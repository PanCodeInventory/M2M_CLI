# M2M_CLI - Matrix2Markers Pipeline CLI

Config-driven Snakemake workflow for single-cell RNA-seq analysis with automatic parameter optimization.

## Features

- **Automatic parameter optimization** using silhouette score
- **Harmony batch correction** (fixed for harmonypy 0.2.0+)
- **Rich metadata support**: tissue, age, group, time
- **Config-driven**: All parameters in YAML config file
- **Complete pipeline**: QC → Merge → Harmony → Optimize → Cluster → Markers

## Pipeline Overview

```
Stage 1: read_10x        → Load 10x data, add metadata (tissue, age)
Stage 2: qc_filter       → QC filtering, doublet removal
Stage 3: merge_and_embed → Merge samples, Harmony correction
Stage 4: optimize        → Find optimal resolution/n_neighbors
Stage 5: apply_clustering → Cluster with optimal params, find markers
```

## Directory Structure

```
M2M_CLI/
├── m2m.py                    # CLI wrapper
├── Snakefile                 # Snakemake workflow
├── config.user.example.yaml  # Config template
├── stages/
│   ├── 01_read_10x_to_h5ad.py
│   ├── 02_qc_filter_single.py
│   ├── 03_merge_and_embed.py
│   ├── 04_clustering_optimization.py
│   └── 05_apply_clustering.py
└── utils/
    └── io_and_qc_utils.py
```

## Quick Start

### 1. Prepare Config File

```bash
# Copy template
cp config.user.example.yaml config.my_project.yaml

# Edit with your sample paths and metadata
vim config.my_project.yaml
```

### 2. Run Pipeline

```bash
# Set environment
export M2M_PIPELINE_ROOT="/path/to/M2M_CLI"

# Validate config
python m2m.py validate --config config.my_project.yaml

# Run workflow
python m2m.py run \
  --config config.my_project.yaml \
  --cores 8 \
  --snakemake /path/to/snakemake
```

## Configuration

### Sample Metadata

```yaml
samples:
  Liver_Aged: "/path/to/Liver_Aged_10x/"
  Liver_Young: "/path/to/Liver_Young_10x/"

sample_metadata:
  Liver_Aged:
    tissue: "Liver"
    age: "Aged"
  Liver_Young:
    tissue: "Liver"
    age: "Young"
```

### Parameter Optimization

```yaml
optimization:
  resolutions: [0.1, 0.2, 0.3, 0.4, 0.5]
  n_neighbors_list: [5, 10, 15, 20, 30]
  skip_neighbors_test: false
```

## Output Structure

```
results/
├── data/
│   ├── merged/
│   │   ├── merged.raw.h5ad        # Raw merged data
│   │   ├── merged.embedded.h5ad   # After Harmony
│   │   └── merged.processed.h5ad  # Final clustered data
│   ├── single_sample_raw/
│   └── single_sample_qc/
├── plots/
│   ├── clustering/                # UMAP, dotplot, rank_genes
│   └── qc_single/                 # QC violin/scatter plots
├── optimization/
│   ├── optimal_params.yaml        # Best parameters found
│   └── *_optimization.png         # Optimization curves
└── tables/
    ├── all_markers.csv            # Marker genes per cluster
    └── qc_single_summary.tsv      # QC summary
```

## Key Results

| File | Description |
|------|-------------|
| `merged.processed.h5ad` | Final data with leiden, tissue, age annotations |
| `all_markers.csv` | Marker genes for each cluster |
| `optimal_params.yaml` | Optimized resolution and n_neighbors |
| `umap_by_leiden.png` | UMAP colored by cluster |
| `umap_by_sample.png` | UMAP colored by sample |

## Requirements

- Python 3.10+
- snakemake 7.x
- scanpy, anndata, harmonypy, sklearn

```bash
conda create -n matrix2marker python=3.10
conda activate matrix2marker
pip install snakemake scanpy anndata harmonypy scikit-learn pyyaml
```

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.

## License

MIT