#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
03_merge_and_embed.py

功能：
- 合并 QC 后的样本
- 保留所有基因的原始计数 (layers['counts'])
- 保留所有基因的归一化数据 (X)
- 仅基于高变基因 (HVGs) 进行 Scale、PCA、Harmony 批次校正
- 在 HVG 基础上按基因名模式剔除 MT/RP/ncRNA/LINC/LOC 等基因
- 输出包含嵌入空间的 h5ad，供后续参数优化使用
- 不进行聚类（聚类参数由 04_clustering_optimization.py 优化）
"""

import os
import sys
import argparse
import traceback
from datetime import datetime

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

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


def summarize_cell_counts(adata: ad.AnnData, out_dir: str):
    """Summarize cell counts by sample and group."""
    ensure_dir(out_dir)
    if "sample" in adata.obs.columns:
        adata.obs.groupby("sample").size().reset_index(name="n_cells").to_csv(
            os.path.join(out_dir, "cell_counts_by_sample.tsv"), sep="\t", index=False
        )
    if "group" in adata.obs.columns:
        adata.obs.groupby("group").size().reset_index(name="n_cells").to_csv(
            os.path.join(out_dir, "cell_counts_by_group.tsv"), sep="\t", index=False
        )


def main():
    parser = argparse.ArgumentParser(description="Merge samples and compute embeddings (no clustering).")
    parser.add_argument("--input-h5ads", nargs="+", required=True, help="List of QC h5ad files")
    parser.add_argument("--output-merged-raw", required=True, help="Output path for merged raw h5ad")
    parser.add_argument("--output-embedded", required=True, help="Output path for embedded h5ad")
    parser.add_argument("--output-plot-dir", required=True, help="Directory for plots")
    parser.add_argument("--output-table-dir", required=True, help="Directory for tables")

    parser.add_argument("--n-top-genes", type=int, default=2000, help="Number of HVGs")
    parser.add_argument("--harmony-key", type=str, default="sample", help="Batch key for Harmony")
    parser.add_argument("--target-sum", type=float, default=1e4, help="Normalization target sum")
    parser.add_argument("--scale-max-value", type=float, default=10.0, help="Max value for scaling")

    parser.add_argument("--mt-pattern", type=str, default="^MT-")
    parser.add_argument("--rp-pattern", type=str, default="^RP[SL]")
    parser.add_argument("--ncrna-pattern", type=str, default="^[A-Z][A-Z][0-9].*\\.[0-9]")
    parser.add_argument("--linc-pattern", type=str, default="(^LOC|LINC)[1-9]*")

    args = parser.parse_args()

    ensure_dir(os.path.dirname(args.output_merged_raw))
    ensure_dir(os.path.dirname(args.output_embedded))
    ensure_dir(args.output_plot_dir)
    ensure_dir(args.output_table_dir)

    try:
        # 1. 合并样本
        adata = merge_qc_samples(args.input_h5ads)
        if not adata.obs_names.is_unique:
            adata.obs_names_make_unique()
        if not adata.var_names.is_unique:
            adata.var_names_make_unique()
        print(f"[{log_timestamp()}] [INFO] Merged: cells={adata.n_obs}, genes={adata.n_vars}")

        adata.write(args.output_merged_raw)
        print(f"[{log_timestamp()}] [INFO] Saved merged.raw.h5ad")

        # 2. 归一化
        adata.layers["counts"] = adata.X.copy()
        if sparse.issparse(adata.layers["counts"]):
            adata.layers["counts"] = adata.layers["counts"].astype(np.int64)
        else:
            adata.layers["counts"] = adata.layers["counts"].astype(np.int64)

        sc.pp.normalize_total(adata, target_sum=args.target_sum)
        sc.pp.log1p(adata)
        adata.raw = adata
        print(f"[{log_timestamp()}] [INFO] Normalized (target_sum={args.target_sum}) and Log1p transformed.")

        # 3. HVG 选择
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=args.n_top_genes,
            layer="counts",
        )

        refined = refine_hvgs_by_gene_patterns(
            adata,
            mt_pattern=args.mt_pattern,
            rp_pattern=args.rp_pattern,
            ncrna_pattern=args.ncrna_pattern,
            linc_pattern=args.linc_pattern,
        )
        hv_mask = adata.var["highly_variable"].values
        print(f"[{log_timestamp()}] [INFO] Final HVGs: {hv_mask.sum()}")

        # 4. Scale + PCA
        adata_hvg = adata[:, hv_mask].copy()
        sc.pp.scale(adata_hvg, max_value=args.scale_max_value)
        sc.pp.pca(adata_hvg, svd_solver="arpack")
        print(f"[{log_timestamp()}] [INFO] PCA computed")

        # 5. Harmony 批次校正
        harmony_success = False
        try:
            import harmonypy as hm
            print(f"[{log_timestamp()}] [INFO] Running Harmony on '{args.harmony_key}'...")

            pca_embeddings = adata_hvg.obsm["X_pca"]
            meta_data_df = adata_hvg.obs[[args.harmony_key]]

            ho = hm.run_harmony(
                pca_embeddings,
                meta_data_df,
                vars_use=[args.harmony_key],
                random_state=0
            )

            adata_hvg.obsm["X_pca_harmony"] = ho.Z_corr
            print(f"[{log_timestamp()}] [INFO] Harmony completed. Shape: {ho.Z_corr.shape}")
            harmony_success = True
        except Exception as e:
            print(f"[{log_timestamp()}] [WARN] Harmony failed: {e}. Will use X_pca instead.")

        # 6. 将嵌入空间写回 adata
        if "X_pca" in adata_hvg.obsm:
            adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
        if "X_pca_harmony" in adata_hvg.obsm:
            adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]

        # 记录使用的嵌入
        adata.uns["embedding_used"] = "X_pca_harmony" if harmony_success else "X_pca"

        del adata_hvg
        import gc
        gc.collect()

        # 7. 保存结果
        adata.write(args.output_embedded)
        print(f"[{log_timestamp()}] [INFO] Saved embedded h5ad to {args.output_embedded}")

        # 8. 输出统计
        summarize_cell_counts(adata, args.output_table_dir)
        print(f"[{log_timestamp()}] [INFO] Done.")

    except Exception as e:
        print(f"[ERROR] Processing failed: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()