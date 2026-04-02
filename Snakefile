configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())

# Helper variables for directories
DIR_RAW = config["directories"]["raw_h5ad"]
DIR_QC = config["directories"]["qc_h5ad"]
DIR_MERGED = config["directories"]["merged"]
DIR_PLOTS = config["directories"]["plots"]
DIR_TABLES = config["directories"]["tables"]
DIR_LOGS = config["directories"]["logs"]
DIR_OPTIMIZATION = config["directories"].get("optimization", os.path.join(DIR_TABLES, "optimization"))

# Python executable
PYTHON = "/home/user/miniforge3/envs/matrix2marker/bin/python"

rule all:
    input:
        f"{DIR_MERGED}/merged.processed.h5ad",
        f"{DIR_TABLES}/qc_single_summary.tsv",
        f"{DIR_TABLES}/all_markers.csv",
        f"{DIR_OPTIMIZATION}/optimal_params.yaml"

# =============================================================================
# Stage 1: Read 10x data
# =============================================================================
rule read_10x:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{DIR_RAW}/{{sample}}.raw.h5ad"
    params:
        group = lambda wildcards: config.get("sample_metadata", {}).get(wildcards.sample, {}).get("group", "NA"),
        time = lambda wildcards: config.get("sample_metadata", {}).get(wildcards.sample, {}).get("time", "NA"),
        tissue = lambda wildcards: config.get("sample_metadata", {}).get(wildcards.sample, {}).get("tissue", "NA"),
        age = lambda wildcards: config.get("sample_metadata", {}).get(wildcards.sample, {}).get("age", "NA")
    log:
        f"{DIR_LOGS}/read_10x/{{sample}}.log"
    shell:
        """
        {PYTHON} stages/01_read_10x_to_h5ad.py \
            --sample-name {wildcards.sample} \
            --input-dir "{input}" \
            --output-file {output} \
            --sample-group {params.group} \
            --sample-time {params.time} \
            --sample-tissue {params.tissue} \
            --sample-age {params.age} > {log} 2>&1
        """

# =============================================================================
# Stage 2: QC and filtering
# =============================================================================
rule qc_filter:
    input:
        f"{DIR_RAW}/{{sample}}.raw.h5ad"
    output:
        h5ad = f"{DIR_QC}/{{sample}}.qc.h5ad",
        stats = f"{DIR_TABLES}/qc_single/{{sample}}_stats.json",
        plot_violin = f"{DIR_PLOTS}/qc_single/{{sample}}_qc_violin.png",
        plot_scatter = f"{DIR_PLOTS}/qc_single/{{sample}}_qc_scatter.png"
    params:
        plot_dir = f"{DIR_PLOTS}/qc_single",
        min_genes = config["filtering"]["min_genes"],
        max_mt = config["filtering"]["max_mt_pct"],
        min_cells = config["filtering"]["min_cells_per_gene"],
        doublet_rate = config["filtering"]["doublet_rate"],
        sim_doublet_ratio = config["filtering"]["sim_doublet_ratio"],
        n_prin_comps = config["filtering"]["n_prin_comps"]
    log:
        f"{DIR_LOGS}/qc_filter/{{sample}}.log"
    shell:
        """
        {PYTHON} stages/02_qc_filter_single.py \
            --input-h5ad {input} \
            --output-h5ad {output.h5ad} \
            --output-plot-dir {params.plot_dir} \
            --output-stats-json {output.stats} \
            --min-genes {params.min_genes} \
            --max-mt-pct {params.max_mt} \
            --min-cells-per-gene {params.min_cells} \
            --doublet-rate {params.doublet_rate} \
            --sim-doublet-ratio {params.sim_doublet_ratio} \
            --n-prin-comps {params.n_prin_comps} > {log} 2>&1
        """

rule aggregate_qc_stats:
    input:
        expand(f"{DIR_TABLES}/qc_single/{{sample}}_stats.json", sample=SAMPLES)
    output:
        f"{DIR_TABLES}/qc_single_summary.tsv"
    run:
        import json
        import pandas as pd
        
        data = []
        for json_file in input:
            with open(json_file) as f:
                data.append(json.load(f))
        
        df = pd.DataFrame(data)
        df.to_csv(output[0], sep="\t", index=False)

# =============================================================================
# Stage 3: Merge and embed (no clustering)
# =============================================================================
rule merge_and_embed:
    input:
        expand(f"{DIR_QC}/{{sample}}.qc.h5ad", sample=SAMPLES)
    output:
        merged_raw = f"{DIR_MERGED}/merged.raw.h5ad",
        embedded = f"{DIR_MERGED}/merged.embedded.h5ad"
    params:
        plot_dir = f"{DIR_PLOTS}/embedding",
        table_dir = DIR_TABLES,
        n_top_genes = config["processing"]["n_top_genes"],
        harmony_key = config["processing"]["harmony_key"],
        target_sum = config["processing"]["target_sum"],
        scale_max_value = config["processing"]["scale_max_value"],
        mt_pattern = config["gene_patterns"]["mt"],
        rp_pattern = config["gene_patterns"]["rp"],
        ncrna_pattern = config["gene_patterns"]["ncrna"],
        linc_pattern = config["gene_patterns"]["linc"]
    log:
        f"{DIR_LOGS}/merge_and_embed.log"
    shell:
        """
        {PYTHON} stages/03_merge_and_embed.py \
            --input-h5ads {input} \
            --output-merged-raw {output.merged_raw} \
            --output-embedded {output.embedded} \
            --output-plot-dir {params.plot_dir} \
            --output-table-dir {params.table_dir} \
            --n-top-genes {params.n_top_genes} \
            --harmony-key {params.harmony_key} \
            --target-sum {params.target_sum} \
            --scale-max-value {params.scale_max_value} \
            --mt-pattern "{params.mt_pattern}" \
            --rp-pattern "{params.rp_pattern}" \
            --ncrna-pattern "{params.ncrna_pattern}" \
            --linc-pattern "{params.linc_pattern}" > {log} 2>&1
        """

# =============================================================================
# Stage 4: Parameter optimization
# =============================================================================
rule optimize_clustering:
    input:
        f"{DIR_MERGED}/merged.embedded.h5ad"
    output:
        params_yaml = f"{DIR_OPTIMIZATION}/optimal_params.yaml"
    params:
        output_dir = DIR_OPTIMIZATION,
        resolutions = config["optimization"]["resolutions"],
        n_neighbors_list = config["optimization"]["n_neighbors_list"],
        skip_flag = "--skip-neighbors-test" if config["optimization"].get("skip_neighbors_test", False) else ""
    log:
        f"{DIR_LOGS}/optimize_clustering.log"
    shell:
        """
        {PYTHON} stages/04_clustering_optimization.py \
            --input-h5ad {input} \
            --output-dir {params.output_dir} \
            --output-params {output.params_yaml} \
            --resolutions {params.resolutions} \
            --n-neighbors-list {params.n_neighbors_list} \
            {params.skip_flag} > {log} 2>&1
        """

# =============================================================================
# Stage 5: Apply clustering and find markers
# =============================================================================
rule apply_clustering:
    input:
        h5ad = f"{DIR_MERGED}/merged.embedded.h5ad",
        params_yaml = f"{DIR_OPTIMIZATION}/optimal_params.yaml"
    output:
        h5ad = f"{DIR_MERGED}/merged.processed.h5ad",
        markers = f"{DIR_TABLES}/all_markers.csv"
    params:
        plot_dir = f"{DIR_PLOTS}/clustering",
        table_dir = DIR_TABLES,
        marker_groupby = config["find_markers"]["groupby"],
        marker_method = config["find_markers"]["method"],
        marker_n_genes = config["find_markers"]["n_genes"],
        marker_n_genes_plot = config["find_markers"]["n_genes_plot"]
    log:
        f"{DIR_LOGS}/apply_clustering.log"
    shell:
        """
        {PYTHON} stages/05_apply_clustering.py \
            --input-h5ad {input.h5ad} \
            --params-yaml {input.params_yaml} \
            --output-h5ad {output.h5ad} \
            --output-markers {output.markers} \
            --output-plot-dir {params.plot_dir} \
            --output-table-dir {params.table_dir} \
            --marker-groupby {params.marker_groupby} \
            --marker-method {params.marker_method} \
            --marker-n-genes {params.marker_n_genes} \
            --marker-n-genes-plot {params.marker_n_genes_plot} > {log} 2>&1
        """