configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())

# Helper variables for directories
DIR_RAW = config["directories"]["raw_h5ad"]
DIR_QC = config["directories"]["qc_h5ad"]
DIR_MERGED = config["directories"]["merged"]
DIR_PLOTS = config["directories"]["plots"]
DIR_TABLES = config["directories"]["tables"]
DIR_LOGS = config["directories"]["logs"]

rule all:
    input:
        f"{DIR_MERGED}/merged.processed.h5ad",
        f"{DIR_TABLES}/qc_single_summary.tsv",
        f"{DIR_TABLES}/all_markers.csv"

rule read_10x:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{DIR_RAW}/{{sample}}.raw.h5ad"
    params:
        group = lambda wildcards: config.get("sample_metadata", {}).get(wildcards.sample, {}).get("group", "NA"),
        time = lambda wildcards: config.get("sample_metadata", {}).get(wildcards.sample, {}).get("time", "NA")
    log:
        f"{DIR_LOGS}/read_10x/{{sample}}.log"
    shell:
        """
        /home/user/miniforge3/envs/matrix2marker/bin/python stages/01_read_10x_to_h5ad.py \
            --sample-name {wildcards.sample} \
            --input-dir "{input}" \
            --output-file {output} \
            --sample-group {params.group} \
            --sample-time {params.time} > {log} 2>&1
        """

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
        /home/user/miniforge3/envs/matrix2marker/bin/python stages/02_qc_filter_single.py \
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

rule merge_and_process:
    input:
        expand(f"{DIR_QC}/{{sample}}.qc.h5ad", sample=SAMPLES)
    output:
        merged_raw = f"{DIR_MERGED}/merged.raw.h5ad",
        merged_processed = f"{DIR_MERGED}/merged.processed.h5ad",
        plot_sample = f"{DIR_PLOTS}/umap/umap_by_sample.png"
    params:
        plot_dir = f"{DIR_PLOTS}/umap",
        table_dir = DIR_TABLES,
        n_top_genes = config["processing"]["n_top_genes"],
        resolution = config["processing"]["resolution"],
        harmony_key = config["processing"]["harmony_key"],
        target_sum = config["processing"]["target_sum"],
        n_neighbors = config["processing"]["n_neighbors"],
        scale_max_value = config["processing"]["scale_max_value"],
        mt_pattern = config["gene_patterns"]["mt"],
        rp_pattern = config["gene_patterns"]["rp"],
        ncrna_pattern = config["gene_patterns"]["ncrna"],
        linc_pattern = config["gene_patterns"]["linc"]
    log:
        f"{DIR_LOGS}/merge_and_process.log"
    shell:
        """
        /home/user/miniforge3/envs/matrix2marker/bin/python stages/03_merge_and_process.py \
            --input-h5ads {input} \
            --output-merged-raw {output.merged_raw} \
            --output-merged-processed {output.merged_processed} \
            --output-plot-dir {params.plot_dir} \
            --output-table-dir {params.table_dir} \
            --n-top-genes {params.n_top_genes} \
            --resolution {params.resolution} \
            --harmony-key {params.harmony_key} \
            --target-sum {params.target_sum} \
            --n-neighbors {params.n_neighbors} \
            --scale-max-value {params.scale_max_value} \
            --mt-pattern "{params.mt_pattern}" \
            --rp-pattern "{params.rp_pattern}" \
            --ncrna-pattern "{params.ncrna_pattern}" \
            --linc-pattern "{params.linc_pattern}" > {log} 2>&1
        """

rule find_markers:
    input:
        f"{DIR_MERGED}/merged.processed.h5ad"
    output:
        table = f"{DIR_TABLES}/all_markers.csv",
        plot_rank = f"{DIR_PLOTS}/markers/rank_genes_{config['find_markers']['groupby']}.png"
    params:
        plot_dir = f"{DIR_PLOTS}/markers",
        groupby = config["find_markers"]["groupby"],
        method = config["find_markers"]["method"],
        n_genes = config["find_markers"]["n_genes"],
        n_genes_plot = config["find_markers"]["n_genes_plot"]
    log:
        f"{DIR_LOGS}/find_markers.log"
    shell:
        """
        /home/user/miniforge3/envs/matrix2marker/bin/python stages/04_find_markers.py \
            --input-h5ad {input} \
            --output-table {output.table} \
            --output-plot-dir {params.plot_dir} \
            --groupby {params.groupby} \
            --method {params.method} \
            --n-genes {params.n_genes} \
            --n-genes-plot {params.n_genes_plot} > {log} 2>&1
        """
