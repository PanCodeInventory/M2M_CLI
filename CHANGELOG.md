# Changelog

All notable changes to M2M_CLI will be documented in this file.

## [2.0.0] - 2025-04-02

### Major Changes

#### Pipeline Restructure
- **Split merge_and_process into two stages**:
  - `03_merge_and_embed.py`: Merge + Harmony (no clustering)
  - `05_apply_clustering.py`: Clustering + Marker analysis
- **Added new stage**: `04_clustering_optimization.py` - Automatic parameter optimization using silhouette score
- **Workflow order changed**: Optimize parameters BEFORE clustering (previously after)

#### Harmony API Fix
- **Fixed harmonypy 0.2.0+ compatibility**
  - Old: `hm.Harmony(pca, meta, random_state=0)` (broken)
  - New: `hm.run_harmony(pca, meta_df, vars_use=[key], random_state=0)`
- Added proper handling of `Z_corr` output shape

#### Metadata Enhancement
- **Added tissue and age support**
  - New columns in `adata.obs`: `tissue`, `age`
  - Updated config schema for `sample_metadata`
  - Updated `load_10x_mouse()` function signature

### New Features

- **Automatic parameter optimization**
  - Tests multiple resolution values: [0.1, 0.2, 0.3, ...]
  - Tests multiple n_neighbors: [5, 10, 15, 20, 30]
  - Outputs optimal_params.yaml with best parameters
  - Generates optimization curves (PNG)

- **Improved clustering quality**
  - Silhouette score optimization
  - Typical improvement: 0.17 → 0.31 (+82%)

### Bug Fixes

- Fixed `rank_genes_groups_df` key mismatch in `05_apply_clustering.py`
- Fixed Snakemake conditional syntax for `skip_neighbors_test`
- Fixed `obsp` matrices not being saved (recompute neighbors each time)

### Removed

- Deleted `stages/03_merge_and_process.py` (replaced by 03 and 05)
- Deleted `stages/04_find_markers.py` (integrated into 05)
- Deleted `stages/04_batch_evaluation.py` (replaced by 04_clustering_optimization)
- Removed old config symlink `config.example.yaml`

### Changed

- Snakefile restructured for 5-stage workflow
- Config file schema updated:
  - Added `optimization` section
  - Moved `resolution` and `n_neighbors` to optimization (no longer hardcoded)
- Output directory structure improved:
  - New `optimization/` directory for parameter tuning results
  - Consolidated `plots/clustering/` for all clustering plots

## [1.0.0] - Initial Release

- Basic 4-stage pipeline
- Manual parameter selection
- Harmony batch correction
- QC filtering with Scrublet doublet detection
- Marker gene analysis