#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
05_apply_clustering.py

功能：
- 读取参数优化结果
- 应用最优参数进行聚类
- 计算 UMAP
- 进行 Marker 基因分析
- 生成可视化结果
"""

import os
import sys
import argparse
import traceback
import yaml
import importlib.util
from datetime import datetime

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
UTILS_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "utils")
if UTILS_DIR not in sys.path:
    sys.path.append(UTILS_DIR)

from io_and_qc_utils import ensure_dir, save_fig


def log_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def summarize_cell_counts(adata: ad.AnnData, out_dir: str, cluster_key: str = "leiden"):
    """Summarize cell counts."""
    ensure_dir(out_dir)
    if "sample" in adata.obs.columns:
        adata.obs.groupby("sample").size().reset_index(name="n_cells").to_csv(
            os.path.join(out_dir, "cell_counts_by_sample.tsv"), sep="\t", index=False
        )
    if "group" in adata.obs.columns:
        adata.obs.groupby("group").size().reset_index(name="n_cells").to_csv(
            os.path.join(out_dir, "cell_counts_by_group.tsv"), sep="\t", index=False
        )
    if cluster_key in adata.obs.columns:
        adata.obs.groupby(cluster_key).size().reset_index(name="n_cells").to_csv(
            os.path.join(out_dir, "cell_counts_by_cluster.tsv"), sep="\t", index=False
        )


def main():
    parser = argparse.ArgumentParser(description="Apply optimal clustering parameters and find markers.")
    parser.add_argument("--input-h5ad", required=True, help="Path to embedded h5ad file")
    parser.add_argument("--params-yaml", required=True, help="Path to optimal parameters YAML")
    parser.add_argument("--output-h5ad", required=True, help="Output path for processed h5ad")
    parser.add_argument("--output-markers", required=True, help="Output path for markers CSV")
    parser.add_argument("--output-plot-dir", required=True, help="Directory for plots")
    parser.add_argument("--output-table-dir", required=True, help="Directory for tables")
    
    parser.add_argument("--marker-groupby", default="leiden", help="Groupby key for marker analysis")
    parser.add_argument("--marker-method", default="wilcoxon", help="Method for marker analysis")
    parser.add_argument("--marker-n-genes", type=int, default=20, help="Top genes per cluster")
    parser.add_argument("--marker-n-genes-plot", type=int, default=5, help="Genes to show in plots")
    
    args = parser.parse_args()
    
    # 读取最优参数
    print(f"[{log_timestamp()}] [INFO] Loading optimal parameters from {args.params_yaml}...")
    with open(args.params_yaml, 'r') as f:
        params = yaml.safe_load(f)
    
    n_neighbors = params['n_neighbors']
    resolution = params['resolution']
    embedding_key = params['embedding_key']
    
    print(f"[{log_timestamp()}] [INFO] Optimal parameters:")
    print(f"  - n_neighbors: {n_neighbors}")
    print(f"  - resolution: {resolution}")
    print(f"  - embedding: {embedding_key}")
    
    # 加载数据
    print(f"[{log_timestamp()}] [INFO] Loading {args.input_h5ad}...")
    adata = sc.read_h5ad(args.input_h5ad)
    print(f"[{log_timestamp()}] [INFO] Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 检查 embedding
    if embedding_key not in adata.obsm:
        print(f"[{log_timestamp()}] [WARN] {embedding_key} not found, trying X_pca...")
        if 'X_pca' in adata.obsm:
            embedding_key = 'X_pca'
        else:
            print(f"[{log_timestamp()}] [ERROR] No valid embedding found!")
            sys.exit(1)
    
    ensure_dir(os.path.dirname(args.output_h5ad))
    ensure_dir(args.output_plot_dir)
    ensure_dir(args.output_table_dir)
    
    try:
        # 1. 构建 neighbors graph
        print(f"[{log_timestamp()}] [INFO] Building neighbors graph (n_neighbors={n_neighbors})...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=embedding_key)
        
        # 2. UMAP
        print(f"[{log_timestamp()}] [INFO] Computing UMAP...")
        sc.tl.umap(adata)
        
        # 3. Leiden 聚类
        cluster_key = "leiden"
        print(f"[{log_timestamp()}] [INFO] Running Leiden clustering (resolution={resolution})...")
        try:
            if importlib.util.find_spec("leidenalg") is None:
                raise ImportError("leidenalg not installed")
            sc.tl.leiden(adata, resolution=resolution)
        except ImportError:
            print(f"[{log_timestamp()}] [WARN] leidenalg not found, using louvain")
            sc.tl.louvain(adata, resolution=resolution)
            cluster_key = "louvain"
        
        n_clusters = adata.obs[cluster_key].nunique()
        print(f"[{log_timestamp()}] [INFO] Found {n_clusters} clusters")
        
        # 4. 保存结果
        adata.write(args.output_h5ad)
        print(f"[{log_timestamp()}] [INFO] Saved processed h5ad to {args.output_h5ad}")
        
        # 5. UMAP 可视化
        print(f"[{log_timestamp()}] [INFO] Generating UMAP plots...")
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
        
        # 6. Marker 基因分析
        print(f"[{log_timestamp()}] [INFO] Running marker gene analysis...")
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=args.marker_method,
            key_added='rank_genes'
        )
        
        # 提取结果 - 对每个 cluster 取 top N genes
        all_markers = []
        for cluster in sorted(adata.obs[cluster_key].unique(), key=int):
            cluster_markers = sc.get.rank_genes_groups_df(
                adata, 
                group=cluster, 
                key='rank_genes'
            ).head(args.marker_n_genes)
            cluster_markers['cluster'] = cluster
            all_markers.append(cluster_markers)
        
        markers_df = pd.concat(all_markers, ignore_index=True)
        # 重排列：cluster 在前
        cols = ['cluster'] + [col for col in markers_df.columns if col != 'cluster']
        markers_df = markers_df[cols]
        markers_df.to_csv(args.output_markers, index=False)
        print(f"[{log_timestamp()}] [INFO] Saved {len(markers_df)} markers ({args.marker_n_genes} per cluster) to {args.output_markers}")
        
        # Marker 可视化
        print(f"[{log_timestamp()}] [INFO] Generating marker plots...")
        
        # Rank genes plot
        fig = sc.pl.rank_genes_groups(
            adata,
            n_genes=args.marker_n_genes_plot,
            key='rank_genes',
            show=False
        )
        save_fig(fig, os.path.join(args.output_plot_dir, f"rank_genes_{cluster_key}.png"))
        
        # Dotplot
        try:
            sc.tl.dendrogram(adata, groupby=cluster_key)
            fig = sc.pl.rank_genes_groups_dotplot(
                adata,
                n_genes=args.marker_n_genes_plot,
                key='rank_genes',
                show=False
            )
            save_fig(fig, os.path.join(args.output_plot_dir, f"dotplot_{cluster_key}.png"))
        except Exception as e:
            print(f"[{log_timestamp()}] [WARN] Dotplot failed: {e}")
        
        # 7. 统计
        summarize_cell_counts(adata, args.output_table_dir, cluster_key)
        
        print(f"\n[{log_timestamp()}] === Summary ===")
        print(f"  Total cells: {adata.n_obs}")
        print(f"  Total genes: {adata.n_vars}")
        print(f"  Clusters: {n_clusters}")
        print(f"  Cluster sizes:")
        print(adata.obs[cluster_key].value_counts().sort_index().to_string())
        
        print(f"\n[{log_timestamp()}] [INFO] Done.")
        
    except Exception as e:
        print(f"[ERROR] Processing failed: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()