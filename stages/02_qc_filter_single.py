#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
02_qc_filter_single.py

功能：
- 读取步骤 01 生成的单样本 raw h5ad 文件
- 计算/补充 QC 指标（如果尚未计算）
- 按统一规则进行 QC 过滤（细胞 + 基因）
- 在单样本层面进行 doublet 检测并去除（基于 scrublet）
- 导出 QC+去双胞 后的单样本 h5ad
- 绘制基础 QC 图（violin & scatter）
- 汇总各样本过滤前后细胞数及 doublet 信息，输出表格

使用说明（在项目根目录 /home/harry/NASH/scRNA_CD45/Scanpy 下运行）：
    python 1_DataProcess/Matrix2H5AD_SC/scripts/stages/02_qc_filter_single.py

注意：
- 运行前需确保已安装 scrublet：
    pip install scrublet
"""

import os
import sys
import argparse
import traceback
import json
from datetime import datetime

import scanpy as sc
import pandas as pd
import anndata as ad

# ========================
# Path Setup
# ========================
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
UTILS_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "utils")
if UTILS_DIR not in sys.path:
    sys.path.append(UTILS_DIR)

from io_and_qc_utils import (
    add_mouse_mt_qc,
    basic_qc_filter,
    ensure_dir,
    save_fig,
    run_doublet_detection_with_scrublet,
)

def log_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def qc_and_filter_single_sample(
    h5ad_path: str,
    out_h5ad_path: str,
    out_plot_dir: str,
    min_genes: int,
    max_mt_pct: float,
    min_cells_per_gene: int,
    doublet_rate: float,
    sim_doublet_ratio: float,
    n_prin_comps: int
) -> dict:
    """
    Process a single h5ad file: QC, Filter, Doublet Detection.
    """
    sample_name = os.path.basename(h5ad_path).replace(".raw.h5ad", "").replace(".h5ad", "")

    # Load data
    print(f"[{log_timestamp()}] [INFO] Loading {h5ad_path}...")
    adata = ad.read_h5ad(h5ad_path)
    
    # Make var_names unique to avoid downstream errors
    adata.var_names_make_unique()

    # Ensure QC metrics exist
    need_qc = False
    for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
        if col not in adata.obs.columns:
            need_qc = True
            break
    
    if need_qc or "mt" not in adata.var.columns:
        print(f"[{log_timestamp()}] [INFO] Calculating QC metrics...")
        add_mouse_mt_qc(adata)

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # 1) Basic QC Filter
    print(f"[{log_timestamp()}] [INFO] Filtering: min_genes={min_genes}, max_mt={max_mt_pct}%, min_cells={min_cells_per_gene}")
    adata_qc = basic_qc_filter(
        adata,
        min_genes=min_genes,
        max_mt_pct=max_mt_pct,
        min_cells_per_gene=min_cells_per_gene,
    )

    n_cells_after_qc = adata_qc.n_obs
    n_genes_after_qc = adata_qc.n_vars

    # 2) Doublet Detection
    print(f"[{log_timestamp()}] [INFO] Running Scrublet (expected_rate={doublet_rate}, sim_ratio={sim_doublet_ratio}, n_pcs={n_prin_comps})...")
    try:
        adata_qc_dedup, dbl_stats = run_doublet_detection_with_scrublet(
            adata_qc,
            expected_doublet_rate=doublet_rate,
            sim_doublet_ratio=sim_doublet_ratio,
            n_prin_comps=n_prin_comps,
            doublet_score_key="doublet_score",
            doublet_label_key="doublet_call",
            random_state=0,
        )
        n_doublets = int(dbl_stats.get("n_doublets", 0))
        doublet_rate_observed = float(dbl_stats.get("doublet_rate_observed", 0.0))
    except ImportError:
        print(f"[{log_timestamp()}] [WARN] Scrublet not installed. Skipping doublet detection.")
        adata_qc_dedup = adata_qc
        n_doublets = 0
        doublet_rate_observed = 0.0

    n_cells_after_qc_dedup = adata_qc_dedup.n_obs
    n_genes_after_qc_dedup = adata_qc_dedup.n_vars

    # 3) Save Output
    print(f"[{log_timestamp()}] [INFO] Saving QC h5ad to {out_h5ad_path}...")
    ensure_dir(os.path.dirname(out_h5ad_path))
    adata_qc_dedup.write(out_h5ad_path)

    # 4) Plots
    print(f"[{log_timestamp()}] [INFO] Generating plots in {out_plot_dir}...")
    ensure_dir(out_plot_dir)
    sc.settings.figdir = out_plot_dir

    # Violin
    fig_violin = sc.pl.violin(
        adata, # Plot original data distribution
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    save_fig(fig_violin, os.path.join(out_plot_dir, f"{sample_name}_qc_violin.png"))

    # Scatter
    fig_scatter = sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        show=False,
    )
    save_fig(fig_scatter, os.path.join(out_plot_dir, f"{sample_name}_qc_scatter.png"))

    return {
        "sample": sample_name,
        "n_cells_before": int(n_cells_before),
        "n_cells_after_qc": int(n_cells_after_qc),
        "n_cells_after_qc_dedup": int(n_cells_after_qc_dedup),
        "n_genes_before": int(n_genes_before),
        "n_genes_after_qc": int(n_genes_after_qc),
        "n_genes_after_qc_dedup": int(n_genes_after_qc_dedup),
        "n_doublets": n_doublets,
        "doublet_rate_observed": doublet_rate_observed,
    }

def main():
    parser = argparse.ArgumentParser(description="QC and Filter single h5ad.")
    parser.add_argument("--input-h5ad", required=True, help="Path to input raw h5ad")
    parser.add_argument("--output-h5ad", required=True, help="Path to output qc h5ad")
    parser.add_argument("--output-plot-dir", required=True, help="Directory for QC plots")
    parser.add_argument("--output-stats-json", required=True, help="Path to output stats JSON")
    
    # Params
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-mt-pct", type=float, default=10.0)
    parser.add_argument("--min-cells-per-gene", type=int, default=3)
    parser.add_argument("--doublet-rate", type=float, default=0.06)
    parser.add_argument("--sim-doublet-ratio", type=float, default=2.0)
    parser.add_argument("--n-prin-comps", type=int, default=30)

    args = parser.parse_args()

    try:
        stats = qc_and_filter_single_sample(
            args.input_h5ad,
            args.output_h5ad,
            args.output_plot_dir,
            args.min_genes,
            args.max_mt_pct,
            args.min_cells_per_gene,
            args.doublet_rate,
            args.sim_doublet_ratio,
            args.n_prin_comps
        )

        with open(args.output_stats_json, "w") as f:
            json.dump(stats, f, indent=4)
        
        print(f"[{log_timestamp()}] [INFO] Stats saved to {args.output_stats_json}")

    except Exception as e:
        print(f"[ERROR] Failed to process {args.input_h5ad}: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
