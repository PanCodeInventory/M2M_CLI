#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
04_find_markers.py

功能：
- 读取合并处理后的 h5ad (merged.processed.h5ad)
- 基于指定的列（如 'leiden', 'group'）计算差异表达基因 (Find Markers)
- 默认使用 Wilcoxon rank-sum test (One-vs-Rest)
- 输出：
    1. 包含所有分组差异基因的 Excel/CSV 表格
    2. Rank genes groups plot (Top 基因排序图)
    3. Dotplot / Matrixplot (Top 基因表达气泡图/热图)
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
import matplotlib.pyplot as plt

# ========================
# Path Setup
# ========================
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
UTILS_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "utils")
if UTILS_DIR not in sys.path:
    sys.path.append(UTILS_DIR)

from io_and_qc_utils import ensure_dir, save_fig

def log_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def find_markers(
    adata: ad.AnnData,
    groupby: str,
    method: str = "wilcoxon",
    n_genes: int = 50,
    pts: bool = True,
    use_raw: bool = True
):
    """
    计算差异基因并返回结果 DataFrame
    """
    print(f"[{log_timestamp()}] [INFO] Computing rank_genes_groups (groupby='{groupby}', method='{method}')...")
    
    # 确保 groupby 列存在且为 category 类型
    if groupby not in adata.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in adata.obs")
    
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        print(f"[{log_timestamp()}] [WARN] Converting '{groupby}' to category...")
        adata.obs[groupby] = adata.obs[groupby].astype("category")

    # 计算差异基因
    # 注意：通常建议在 normalized/log1p 后的数据上做 diff exp
    # 如果 adata.raw 存在且 use_raw=True，scanpy 会使用 adata.raw.X
    
    # 强制 n_genes=None 以计算所有基因的分数，满足后续 dotplot 等绘图需求
    print(f"[{log_timestamp()}] [INFO] Running rank_genes_groups with n_genes=None (all genes) for plotting compatibility...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=None,
        use_raw=use_raw,
        pts=pts, # 计算表达比例
    )
    
    # 整理结果为 DataFrame
    # sc.get.rank_genes_groups_df 可以方便提取
    # 我们将所有组的结果合并为一个大表
    groups = adata.obs[groupby].cat.categories
    dfs = []
    
    # 如果传入的 n_genes <= 0，则保留所有基因
    limit_genes = n_genes if n_genes > 0 else None
    
    print(f"[{log_timestamp()}] [INFO] Extracting top {limit_genes if limit_genes else 'ALL'} genes per group for table...")

    for group in groups:
        df = sc.get.rank_genes_groups_df(adata, group=group, key="rank_genes_groups")
        
        # 截取前 n_genes 行
        if limit_genes is not None:
            df = df.head(limit_genes)
            
        df["group"] = group
        df["comparison"] = f"{group} vs Rest"
        dfs.append(df)
    
    result_df = pd.concat(dfs, axis=0)
    
    # 调整列顺序，把 group 放在前面
    cols = ["group", "names", "scores", "pvals", "pvals_adj", "logfoldchanges"]
    if pts:
        cols.extend(["pct_nz_group", "pct_nz_reference"])
    
    # 仅保留存在的列
    cols = [c for c in cols if c in result_df.columns]
    result_df = result_df[cols]
    
    return result_df

def plot_markers(adata: ad.AnnData, groupby: str, plot_dir: str, n_genes_plot: int = 5):
    """
    绘制差异基因相关图形
    """
    ensure_dir(plot_dir)
    sc.settings.figdir = plot_dir
    
    # 1. Rank Genes Groups Plot
    print(f"[{log_timestamp()}] [INFO] Plotting rank_genes_groups...")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    save_fig(plt.gcf(), os.path.join(plot_dir, f"rank_genes_{groupby}.png"))
    
    # 2. Dotplot (Top N distinct genes)
    print(f"[{log_timestamp()}] [INFO] Plotting dotplot...")
    # 获取每一组最显著的基因
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=n_genes_plot,
        values_to_plot="logfoldchanges", 
        cmap="bwr",
        vmin=-4, vmax=4,
        min_logfoldchange=0.5, # 仅展示变化明显的
        colorbar_title="logFC",
        show=False
    )
    save_fig(plt.gcf(), os.path.join(plot_dir, f"dotplot_{groupby}.png"))

    # 3. Matrixplot (Heatmap style)
    print(f"[{log_timestamp()}] [INFO] Plotting matrixplot...")
    sc.pl.rank_genes_groups_matrixplot(
        adata,
        n_genes=n_genes_plot,
        values_to_plot="logfoldchanges",
        cmap="bwr",
        vmin=-4, vmax=4,
        colorbar_title="logFC",
        show=False
    )
    save_fig(plt.gcf(), os.path.join(plot_dir, f"matrixplot_{groupby}.png"))

def main():
    parser = argparse.ArgumentParser(description="Find Markers (Differential Expression).")
    parser.add_argument("--input-h5ad", required=True, help="Path to input processed h5ad")
    parser.add_argument("--output-table", required=True, help="Path to output CSV/Excel table")
    parser.add_argument("--output-plot-dir", required=True, help="Directory for marker plots")
    
    # Params
    parser.add_argument("--groupby", type=str, default="leiden", help="Column to group cells by")
    parser.add_argument("--method", type=str, default="wilcoxon", help="Test method (wilcoxon, t-test, etc.)")
    parser.add_argument("--n-genes", type=int, default=100, help="Number of genes to report in table")
    parser.add_argument("--n-genes-plot", type=int, default=5, help="Number of top genes to plot")

    args = parser.parse_args()

    ensure_dir(os.path.dirname(args.output_table))
    ensure_dir(args.output_plot_dir)

    try:
        # 1. Load Data
        print(f"[{log_timestamp()}] [INFO] Loading {args.input_h5ad}...")
        adata = ad.read_h5ad(args.input_h5ad)
        
        # 2. Find Markers
        result_df = find_markers(
            adata,
            groupby=args.groupby,
            method=args.method,
            n_genes=args.n_genes,
        )
        
        # 3. Save Table
        print(f"[{log_timestamp()}] [INFO] Saving table to {args.output_table}...")
        if args.output_table.endswith(".csv"):
            result_df.to_csv(args.output_table, index=False)
        elif args.output_table.endswith(".xlsx"):
            result_df.to_excel(args.output_table, index=False)
        else:
            result_df.to_csv(args.output_table, sep="\t", index=False)
            
        # 4. Plots
        plot_markers(adata, args.groupby, args.output_plot_dir, n_genes_plot=args.n_genes_plot)
        
        print(f"[{log_timestamp()}] [INFO] Done.")

    except Exception as e:
        print(f"[ERROR] Find Markers failed: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
