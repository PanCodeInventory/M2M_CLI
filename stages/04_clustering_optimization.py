#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
04_clustering_optimization.py

功能：
- 使用轮廓系数 (Silhouette Score) 评估聚类参数选择
- 测试不同的 resolution 参数
- 测试不同的 n_neighbors 参数
- 输出最优参数到 YAML 文件
- 生成参数优化报告和可视化
"""

import os
import sys
import argparse
import yaml
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
UTILS_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "utils")
if UTILS_DIR not in sys.path:
    sys.path.append(UTILS_DIR)

from io_and_qc_utils import ensure_dir


def log_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def test_resolutions(adata, emb_key, resolutions, n_neighbors):
    """
    测试不同 resolution 参数的聚类效果。
    """
    results = []
    
    # 构建 neighbors graph
    print(f"[{log_timestamp()}] [INFO] Building neighbors graph (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=emb_key, key_added='temp_neighbors')
    
    for res in resolutions:
        # 聚类
        sc.tl.leiden(adata, resolution=res, neighbors_key='temp_neighbors', key_added=f'temp_r{res}')
        
        # 计算轮廓系数
        labels = adata.obs[f'temp_r{res}'].values
        n_clusters = len(np.unique(labels))
        
        # 采样计算
        n_samples = min(10000, len(labels))
        if n_samples < len(labels):
            idx = np.random.choice(len(labels), n_samples, replace=False)
            emb = adata.obsm[emb_key][idx]
            lab = labels[idx]
        else:
            emb = adata.obsm[emb_key]
            lab = labels
        
        score = silhouette_score(emb, lab)
        
        results.append({
            'resolution': res,
            'n_clusters': n_clusters,
            'silhouette_score': round(score, 4)
        })
        
        print(f"[{log_timestamp()}] [INFO] resolution={res}: n_clusters={n_clusters}, silhouette={score:.4f}")
    
    return pd.DataFrame(results)


def test_n_neighbors(adata, emb_key, n_neighbors_list):
    """
    测试不同 n_neighbors 参数对聚类的影响（使用默认 resolution=0.5）。
    """
    results = []
    
    for n_neigh in n_neighbors_list:
        print(f"[{log_timestamp()}] [INFO] Testing n_neighbors={n_neigh}...")
        
        # 构建 neighbor graph
        sc.pp.neighbors(adata, n_neighbors=n_neigh, use_rep=emb_key, key_added=f'temp_nn{n_neigh}')
        
        # 聚类
        sc.tl.leiden(adata, resolution=0.5, neighbors_key=f'temp_nn{n_neigh}', key_added=f'temp_leiden_nn{n_neigh}')
        
        # 计算轮廓系数
        labels = adata.obs[f'temp_leiden_nn{n_neigh}'].values
        n_clusters = len(np.unique(labels))
        
        n_samples = min(10000, len(labels))
        if n_samples < len(labels):
            idx = np.random.choice(len(labels), n_samples, replace=False)
            emb = adata.obsm[emb_key][idx]
            lab = labels[idx]
        else:
            emb = adata.obsm[emb_key]
            lab = labels
        
        score = silhouette_score(emb, lab)
        
        results.append({
            'n_neighbors': n_neigh,
            'n_clusters': n_clusters,
            'silhouette_score': round(score, 4)
        })
    
    return pd.DataFrame(results)


def plot_resolution_optimization(df, output_path):
    """绘制 resolution 优化图。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 轮廓系数 vs resolution
    ax1 = axes[0]
    ax1.plot(df['resolution'], df['silhouette_score'], 'b-o', linewidth=2, markersize=8)
    ax1.set_xlabel('Resolution', fontsize=12)
    ax1.set_ylabel('Silhouette Score', fontsize=12)
    ax1.set_title('Clustering Quality vs Resolution', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # 标记最优点
    best_idx = df['silhouette_score'].idxmax()
    best_res = df.loc[best_idx, 'resolution']
    best_score = df.loc[best_idx, 'silhouette_score']
    ax1.axvline(best_res, color='red', linestyle='--', alpha=0.5)
    ax1.scatter([best_res], [best_score], color='red', s=150, zorder=5, label=f'Best: res={best_res}')
    ax1.legend()
    
    # n_clusters vs resolution
    ax2 = axes[1]
    ax2.plot(df['resolution'], df['n_clusters'], 'g-s', linewidth=2, markersize=8)
    ax2.set_xlabel('Resolution', fontsize=12)
    ax2.set_ylabel('Number of Clusters', fontsize=12)
    ax2.set_title('Cluster Count vs Resolution', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_n_neighbors_optimization(df, output_path):
    """绘制 n_neighbors 优化图。"""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.plot(df['n_neighbors'], df['silhouette_score'], 'b-o', linewidth=2, markersize=8)
    ax.set_xlabel('n_neighbors', fontsize=12)
    ax.set_ylabel('Silhouette Score', fontsize=12)
    ax.set_title('Clustering Quality vs n_neighbors', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # 标记最优点
    best_idx = df['silhouette_score'].idxmax()
    best_nn = df.loc[best_idx, 'n_neighbors']
    best_score = df.loc[best_idx, 'silhouette_score']
    ax.axvline(best_nn, color='red', linestyle='--', alpha=0.5)
    ax.scatter([best_nn], [best_score], color='red', s=150, zorder=5, label=f'Best: nn={best_nn}')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Optimize clustering parameters using silhouette scores.")
    parser.add_argument("--input-h5ad", required=True, help="Path to embedded h5ad file")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--output-params", required=True, help="Output YAML file for optimal parameters")
    parser.add_argument("--embedding-key", default="X_pca_harmony", help="Embedding to use")
    
    parser.add_argument("--resolutions", type=float, nargs='+',
                        default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5],
                        help="Resolution values to test")
    parser.add_argument("--n-neighbors-list", type=int, nargs='+',
                        default=[5, 10, 15, 20, 30],
                        help="n_neighbors values to test")
    parser.add_argument("--skip-neighbors-test", action='store_true',
                        help="Skip n_neighbors optimization")
    
    args = parser.parse_args()
    
    print(f"[{log_timestamp()}] [INFO] Loading {args.input_h5ad}...")
    adata = sc.read_h5ad(args.input_h5ad)
    print(f"[{log_timestamp()}] [INFO] Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 检查 embedding
    if args.embedding_key not in adata.obsm:
        print(f"[{log_timestamp()}] [WARN] {args.embedding_key} not found, trying X_pca...")
        if 'X_pca' in adata.obsm:
            args.embedding_key = 'X_pca'
        else:
            print(f"[{log_timestamp()}] [ERROR] No valid embedding found!")
            sys.exit(1)
    
    print(f"[{log_timestamp()}] [INFO] Using embedding: {args.embedding_key}")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. 测试不同 n_neighbors（先确定，因为会影响 resolution 测试）
    best_n_neighbors = 10  # 默认值
    if not args.skip_neighbors_test:
        print(f"\n[{log_timestamp()}] === Step 1: Testing n_neighbors ===")
        print(f"[{log_timestamp()}] [INFO] n_neighbors to test: {args.n_neighbors_list}")
        
        nn_df = test_n_neighbors(adata, args.embedding_key, args.n_neighbors_list)
        
        # 保存结果
        nn_output = os.path.join(args.output_dir, "n_neighbors_optimization.tsv")
        nn_df.to_csv(nn_output, sep="\t", index=False)
        print(f"[{log_timestamp()}] [INFO] Saved: {nn_output}")
        
        # 绘图
        plot_n_neighbors_optimization(nn_df, os.path.join(args.output_dir, "n_neighbors_optimization.png"))
        
        # 找最优 n_neighbors
        best_nn_idx = nn_df['silhouette_score'].idxmax()
        best_n_neighbors = int(nn_df.loc[best_nn_idx, 'n_neighbors'])
        
        print(f"\n[{log_timestamp()}] === Best n_neighbors ===")
        print(f"  n_neighbors: {best_n_neighbors}")
        print(f"  Silhouette Score: {nn_df.loc[best_nn_idx, 'silhouette_score']:.4f}")
    
    # 2. 测试不同 resolution
    print(f"\n[{log_timestamp()}] === Step 2: Testing Resolutions ===")
    print(f"[{log_timestamp()}] [INFO] Resolutions to test: {args.resolutions}")
    print(f"[{log_timestamp()}] [INFO] Using n_neighbors: {best_n_neighbors}")
    
    res_df = test_resolutions(adata, args.embedding_key, args.resolutions, best_n_neighbors)
    
    # 保存结果
    res_output = os.path.join(args.output_dir, "resolution_optimization.tsv")
    res_df.to_csv(res_output, sep="\t", index=False)
    print(f"[{log_timestamp()}] [INFO] Saved: {res_output}")
    
    # 绘图
    plot_resolution_optimization(res_df, os.path.join(args.output_dir, "resolution_optimization.png"))
    
    # 找最优 resolution
    best_idx = res_df['silhouette_score'].idxmax()
    best_res = float(res_df.loc[best_idx, 'resolution'])
    best_score = float(res_df.loc[best_idx, 'silhouette_score'])
    best_n_clusters = int(res_df.loc[best_idx, 'n_clusters'])
    
    print(f"\n[{log_timestamp()}] === Best Resolution ===")
    print(f"  Resolution: {best_res}")
    print(f"  Silhouette Score: {best_score:.4f}")
    print(f"  Number of Clusters: {best_n_clusters}")
    
    # 3. 保存最优参数到 YAML
    optimal_params = {
        'n_neighbors': best_n_neighbors,
        'resolution': best_res,
        'embedding_key': args.embedding_key,
        'silhouette_score': round(best_score, 4),
        'n_clusters': best_n_clusters
    }
    
    ensure_dir(os.path.dirname(args.output_params))
    with open(args.output_params, 'w') as f:
        yaml.dump(optimal_params, f, default_flow_style=False)
    print(f"\n[{log_timestamp()}] [INFO] Saved optimal parameters to: {args.output_params}")
    
    # 4. 生成总结报告
    print(f"\n[{log_timestamp()}] === Summary Report ===\n")
    print("=" * 60)
    print("CLUSTERING PARAMETER OPTIMIZATION REPORT")
    print("=" * 60)
    print(f"\nDataset: {adata.n_obs} cells")
    print(f"Embedding: {args.embedding_key}")
    print(f"\n--- Optimal Parameters ---")
    print(f"  n_neighbors: {best_n_neighbors}")
    print(f"  resolution: {best_res}")
    print(f"  Expected clusters: {best_n_clusters}")
    print(f"  Silhouette score: {best_score:.4f}")
    
    # 判断聚类质量
    if best_score > 0.5:
        quality = "Excellent"
    elif best_score > 0.3:
        quality = "Good"
    elif best_score > 0.1:
        quality = "Moderate"
    else:
        quality = "Poor (consider adjusting parameters)"
    
    print(f"\n--- Quality Assessment ---")
    print(f"  Clustering quality: {quality}")
    print("\n" + "=" * 60)
    
    print(f"\n[{log_timestamp()}] [INFO] Done.")


if __name__ == "__main__":
    main()