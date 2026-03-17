#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
03_merge_and_process.py

功能：
- 合并 QC 后的样本
- 保留所有基因的原始计数 (layers['counts'])
- 保留所有基因的归一化数据 (X)
- 仅基于高变基因 (HVGs) 进行 Scale、PCA、去批次、聚类
- 在 HVG 基础上按基因名模式剔除 MT/RP/ncRNA/LINC/LOC 等基因
- 最终输出包含全基因组表达量 + 基于“净化后 HVG”的聚类结果的 h5ad
"""

import os
import sys
import argparse
import traceback
import json
from datetime import datetime

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse
import importlib
import scanpy.external as sce

# ========================
# Path Setup
# ========================
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
UTILS_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "utils")
if UTILS_DIR not in sys.path:
    sys.path.append(UTILS_DIR)

from io_and_qc_utils import ensure_dir, save_fig, refine_hvgs_by_gene_patterns

def log_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def merge_qc_samples(input_files: list) -> ad.AnnData:
    """Read and merge QC samples."""
    adatas = []
    sample_keys = []

    print(f"[{log_timestamp()}] [INFO] Merging {len(input_files)} samples...")
    for path in input_files:
        sample_name = os.path.basename(path).replace(".qc.h5ad", "").replace(".h5ad", "")
        print(f"  - Loading {sample_name} from {path}")
        adata = ad.read_h5ad(path)
        adatas.append(adata)
        sample_keys.append(sample_name)

    adata_merged = ad.concat(
        adatas,
        join="outer",
        label="sample_from_file",
        keys=sample_keys,
        index_unique="-",
    )
    return adata_merged

def summarize_cell_counts(adata: ad.AnnData, out_dir: str, cluster_key: str = "leiden"):
    """Summarize cell counts."""
    ensure_dir(out_dir)
    if "sample" in adata.obs.columns:
        adata.obs.groupby("sample").size().reset_index(name="n_cells").to_csv(os.path.join(out_dir, "cell_counts_by_sample.tsv"), sep="\t", index=False)
    if "group" in adata.obs.columns:
        adata.obs.groupby("group").size().reset_index(name="n_cells").to_csv(os.path.join(out_dir, "cell_counts_by_group.tsv"), sep="\t", index=False)
    if cluster_key in adata.obs.columns:
        adata.obs.groupby(cluster_key).size().reset_index(name="n_cells").to_csv(os.path.join(out_dir, "cell_counts_by_cluster.tsv"), sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Merge and process scRNA-seq samples.")
    parser.add_argument("--input-h5ads", nargs="+", required=True, help="List of QC h5ad files")
    parser.add_argument("--output-merged-raw", required=True, help="Output path for merged raw h5ad")
    parser.add_argument("--output-merged-processed", required=True, help="Output path for processed h5ad")
    parser.add_argument("--output-plot-dir", required=True, help="Directory for plots")
    parser.add_argument("--output-table-dir", required=True, help="Directory for tables")
    
    # Processing Params
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--harmony-key", type=str, default="sample")
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--scale-max-value", type=float, default=10.0)

    # Gene Patterns
    parser.add_argument("--mt-pattern", type=str, default="^MT-")
    parser.add_argument("--rp-pattern", type=str, default="^RP[SL]")
    parser.add_argument("--ncrna-pattern", type=str, default="^[A-Z][A-Z][0-9].*\\.[0-9]")
    parser.add_argument("--linc-pattern", type=str, default="(^LOC|LINC)[1-9]*")

    args = parser.parse_args()

    ensure_dir(os.path.dirname(args.output_merged_raw))
    ensure_dir(os.path.dirname(args.output_merged_processed))
    ensure_dir(args.output_plot_dir)
    ensure_dir(args.output_table_dir)

    try:
        # 1. Merge
        adata = merge_qc_samples(args.input_h5ads)
        if not adata.obs_names.is_unique:
            adata.obs_names_make_unique()
        if not adata.var_names.is_unique:
            adata.var_names_make_unique()
        print(f"[{log_timestamp()}] [INFO] Merged: cells={adata.n_obs}, genes={adata.n_vars}")

        # 1.1 Save Merged Raw
        adata.write(args.output_merged_raw)
        print(f"[{log_timestamp()}] [INFO] Saved merged.raw.h5ad")

        # 2. Backup counts
        adata.layers["counts"] = adata.X.copy()
        if sparse.issparse(adata.layers["counts"]):
            adata.layers["counts"] = adata.layers["counts"].astype(np.int64)
        else:
            adata.layers["counts"] = adata.layers["counts"].astype(np.int64)

        # 3. Normalize & Log1p
        sc.pp.normalize_total(adata, target_sum=args.target_sum)
        sc.pp.log1p(adata)
        adata.raw = adata
        print(f"[{log_timestamp()}] [INFO] Normalized (target_sum={args.target_sum}) and Log1p transformed.")

        # 4. HVG
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=args.n_top_genes,
            layer="counts",
        )
        
        # 4.1 Refine HVG
        refined = refine_hvgs_by_gene_patterns(
            adata,
            mt_pattern=args.mt_pattern,
            rp_pattern=args.rp_pattern,
            ncrna_pattern=args.ncrna_pattern,
            linc_pattern=args.linc_pattern
        )
        hv_mask = adata.var["highly_variable"].values
        print(f"[{log_timestamp()}] [INFO] Final HVGs: {hv_mask.sum()}")

        # 5. Process on Copy
        adata_hvg = adata[:, hv_mask].copy()
        sc.pp.scale(adata_hvg, max_value=args.scale_max_value)
        sc.pp.pca(adata_hvg, svd_solver="arpack")

        # 6. Harmony (直接使用 harmonypy)
        import harmonypy as hm
        print(f"[{log_timestamp()}] [INFO] Running Harmony on '{args.harmony_key}'...")
        
        # Get batch column and PCA data
        pca_embeddings = adata_hvg.obsm["X_pca"]
        
        # Get metadata for Harmony
        meta_data = adata_hvg.obs[[args.harmony_key]].values
        
        # Run Harmony
        ho = hm.Harmony(pca_embeddings, meta_data, random_state=0)
        
        # Store corrected embeddings
        adata_hvg.obsm["X_pca_harmony"] = ho.Z_corr
        print(f"[{log_timestamp()}] [INFO] Harmony completed. corrected embeddings shape: {ho.Z_corr.shape}")
        
        # Continue with neighbors using corrected embeddings
        sc.pp.neighbors(adata_hvg, n_neighbors=args.n_neighbors, use_rep="X_pca_harmony")
        print(f"[{log_timestamp()}] [INFO] Neighbors done (using X_pca_harmony)")
    except Exception as e:
        print(f"[{log_timestamp()}] [WARN] Harmony failed: {e}")
        # Continue with PCA and neighbors...
        rep_key = "X_pca"
        sc.pp.neighbors(adata_hvg, n_neighbors=args.n_neighbors, use_rep=rep_key)
            sc.tl.umap(adata_hvg, title="UMAP by sample", show=False)
            save_fig(fig, os.path.join(args.output_plot_dir, "umap_by_sample.png"))
        if "group" in adata.obs.columns:
            fig = sc.pl.umap(adata, color="group", title="UMAP by group", show=False)
            save_fig(fig, os.path.join(args.output_plot_dir, "umap_by_group.png"))
        if cluster_key in aata.obs.columns:
            fig = sc.pl.umap(adata, color=cluster_key, title=f"UMAP by {cluster_key}", show=False)
            save_fig(fig, os.path.join(args.output_plot_dir, f"umap_by_{cluster_key}.png")

        # 12. Stats
        summarize_cell_counts(adata, args.output_table_dir, cluster_key)
        print(f"[{log_timestamp()}] [INFO] Done.")

        # 7. Neighbors & UMAP
        sc.pp.neighbors(adata_hvg, n_neighbors=args.n_neighbors, use_rep=rep_key)
        sc.tl.umap(adata_hvg)

        # 8. Clustering
        cluster_key = "leiden"
        try:
            if importlib.util.find_spec("leidenalg") is None:
                raise ImportError("leidenalg not installed")
            sc.tl.leiden(adata_hvg, resolution=args.resolution)
        except ImportError:
            sc.tl.louvain(adata_hvg, resolution=args.resolution)
            cluster_key = "louvain"
        
        # 9. Sync back
        adata.obs[cluster_key] = adata_hvg.obs[cluster_key]
        if "X_pca" in adata_hvg.obsm: adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
        if "X_pca_harmony" in adata_hvg.obsm: adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]
        if "X_umap" in adata_hvg.obsm: adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
        if "neighbors" in adata_hvg.uns: adata.uns["neighbors"] = adata_hvg.uns["neighbors"]
        if "leiden" in adata_hvg.uns: adata.uns["leiden"] = adata_hvg.uns["leiden"]
        if "louvain" in adata_hvg.uns: adata.uns["louvain"] = adata_hvg.uns["louvain"]

        del adata_hvg
        import gc
        gc.collect()

        # 10. Save Processed
        adata.write(args.output_merged_processed)
        print(f"[{log_timestamp()}] [INFO] Saved merged.processed.h5ad")

        # 11. Plots
        sc.settings.figdir = args.output_plot_dir
        if "sample" in adata.obs.columns:
            fig = sc.pl.umap(adata, color="sample", title="UMAP by sample", show=False)
            save_fig(fig, os.path.join(args.output_plot_dir, "umap_by_sample.png"))
        if "group" in adata.obs.columns:
            fig = sc.pl.umap(adata, color="group", title="UMAP by group", show=False)
            save_fig(fig, os.path.join(args.output_plot_dir, "umap_by_group.png"))
        if cluster_key in adata.obs.columns:
            fig = sc.pl.umap(adata, color=cluster_key, title=f"UMAP by {cluster_key}", show=False)
            save_fig(fig, os.path.join(args.output_plot_dir, f"umap_by_{cluster_key}.png"))

        # 12. Stats
        summarize_cell_counts(adata, args.output_table_dir, cluster_key)
        print(f"[{log_timestamp()}] [INFO] Done.")

    except Exception as e:
        print(f"[ERROR] Processing failed: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
