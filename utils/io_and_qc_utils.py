#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
io_and_qc_utils.py

通用 I/O 与 QC 工具函数：
1. 手动从 10x mtx + features + barcodes 构建 AnnData（不再使用 read_10x_mtx）
2. 计算 QC 指标（n_genes、total_counts、pct_counts_mt 等）
3. 基于简单规则进行细胞与基因过滤
4. 保存图形到指定路径
5. 高变基因过滤（按基因名模式剔除不希望驱动聚类的基因）
6. 单样本层面的 doublet 检测与去除（基于 scrublet）

注意：
- 明确避开 scanpy.read_10x_mtx 内部对 genes[1] 的假设，避免 KeyError: 1。
"""

import os
from typing import Dict, List, Tuple

import scanpy as sc
import anndata as ad
import numpy as np


def load_10x_mouse(
    sample_dir: str, 
    sample_name: str,
    group: str = "NA",
    time: str = "NA"
) -> ad.AnnData:
    """
    从 10x 目录读取小鼠单细胞数据，返回 AnnData 对象。
    
    参数:
    - sample_dir: 10x matrix 所在目录
    - sample_name: 样本名称
    - group: 样本分组 (默认为 "NA")
    - time: 时间点 (默认为 "NA")
    """
    import pandas as pd

    # 1）定位 matrix.mtx(.gz) 或 Alevin quants_mat.mtx
    mtx_path = None
    file_format = "10x"  # "10x" or "alevin"

    # Try 10x
    for fname in ["matrix.mtx.gz", "matrix.mtx"]:
        candidate = os.path.join(sample_dir, fname)
        if os.path.exists(candidate):
            mtx_path = candidate
            file_format = "10x"
            break
    
    # Try Alevin if 10x not found
    if mtx_path is None:
        candidate = os.path.join(sample_dir, "quants_mat.mtx")
        if os.path.exists(candidate):
            mtx_path = candidate
            file_format = "alevin"

    if mtx_path is None:
        raise FileNotFoundError(f"未在 {sample_dir} 中找到 matrix.mtx(.gz) 或 quants_mat.mtx")

    # 2）定位 features / genes 文件
    features_path = None
    if file_format == "10x":
        for fname in ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]:
            candidate = os.path.join(sample_dir, fname)
            if os.path.exists(candidate):
                features_path = candidate
                break
    else: # Alevin
        candidate = os.path.join(sample_dir, "quants_mat_cols.txt")
        if os.path.exists(candidate):
            features_path = candidate

    if features_path is None:
        raise FileNotFoundError(f"未在 {sample_dir} 中找到 features/genes 文件 (或 quants_mat_cols.txt)")

    # 3）定位 barcodes 文件
    barcodes_path = None
    if file_format == "10x":
        for fname in ["barcodes.tsv.gz", "barcodes.tsv"]:
            candidate = os.path.join(sample_dir, fname)
            if os.path.exists(candidate):
                barcodes_path = candidate
                break
    else: # Alevin
        candidate = os.path.join(sample_dir, "quants_mat_rows.txt")
        if os.path.exists(candidate):
            barcodes_path = candidate

    if barcodes_path is None:
        raise FileNotFoundError(f"未在 {sample_dir} 中找到 barcodes 文件 (或 quants_mat_rows.txt)")

    # 4）读取 matrix.mtx 为 AnnData，然后取 X（n_genes x n_cells 或 n_cells x n_genes）
    adata_mtx = sc.read_mtx(mtx_path)
    X = adata_mtx.X

    # 5）读取 features 和 barcodes
    # 10x features 通常无 header (header=None)
    # Alevin cols/rows 通常也无 header (header=None)
    features = pd.read_csv(features_path, header=None, sep="\t")
    # 如果 features 只有一列且看起来像 csv，尝试用逗号分隔 (Alevin txt 可能有时不一样，但通常是单列)
    # 这里保持默认 sep="\t" 也可以兼容单列数据（pandas 会读整行为一列）
    
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")

    # 6）根据列数构造 var_names 和 var 字段
    n_feature_cols = features.shape[1]
    if n_feature_cols == 1:
        # 只有 gene_id 一列
        gene_ids = features.iloc[:, 0].astype(str).values
        var_names = gene_ids
        var_dict = {"gene_ids": gene_ids}
    else:
        # 至少两列：第 1 列 gene_id，第 2 列 gene_symbol
        gene_ids = features.iloc[:, 0].astype(str).values
        gene_symbols = features.iloc[:, 1].astype(str).values
        var_names = gene_symbols
        var_dict = {
            "gene_ids": gene_ids,
            "gene_symbol": gene_symbols,
        }

    # 7）obs_names（细胞条形码）
    obs_names = barcodes.iloc[:, 0].astype(str).values

    # 8）校验矩阵方向并必要时转置为 AnnData 期望的形状：行=细胞(n_cells), 列=基因(n_genes)
    n_genes = len(var_names)
    n_cells = len(obs_names)

    if X.shape == (n_genes, n_cells):
        # 当前为 行=基因, 列=细胞（10x 原生形状），转置为 行=细胞, 列=基因
        X = X.T
    elif X.shape == (n_cells, n_genes):
        # 已是 行=细胞, 列=基因，直接使用
        pass
    else:
        # 其余形状皆为异常，给出详细提示便于排查
        raise ValueError(
            f"matrix.mtx 形状={X.shape} 与 (n_genes={n_genes}, n_cells={n_cells}) 不匹配，"
            "请检查 mtx/features/barcodes 是否对应同一数据集、或文件是否损坏。"
        )

    # 9）构建 AnnData（行=细胞，列=基因）
    adata = ad.AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = var_names

    for k, v in var_dict.items():
        adata.var[k] = v

    # 10）添加样本信息及 group/time 注释
    adata.obs["sample"] = sample_name
    adata.obs["group"] = group
    adata.obs["time"] = time

    # 11）确保基因名唯一（scanpy 要求）
    adata.var_names_make_unique()

    return adata


def add_mouse_mt_qc(adata: ad.AnnData) -> None:
    """
    针对小鼠数据，标记线粒体基因并计算 QC 指标。
    """
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True
    )


def basic_qc_filter(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_mt_pct: float = 10.0,
    min_cells_per_gene: int = 3
) -> ad.AnnData:
    """
    对 AnnData 进行基础 QC 过滤，返回过滤后的新 AnnData 对象。
    """
    cell_mask = adata.obs["n_genes_by_counts"] >= min_genes
    cell_mask &= adata.obs["pct_counts_mt"] < max_mt_pct
    adata_cell_filtered = adata[cell_mask, :].copy()
    sc.pp.filter_genes(adata_cell_filtered, min_cells=min_cells_per_gene)
    return adata_cell_filtered


def ensure_dir(path: str) -> None:
    """
    确保目录存在，如不存在则创建。
    """
    if path == "" or path is None:
        return
    os.makedirs(path, exist_ok=True)


def save_fig(obj, out_path: str) -> None:
    """
    保存各种绘图对象（matplotlib Figure/Axes、seaborn FacetGrid、scanpy 返回对象等）到指定路径，
    自动创建上级目录，并安全关闭图形以节省内存。

    兼容性说明：
    - scanpy.pl.violin 通常返回 seaborn.FacetGrid（有 savefig 方法，但无 clf 方法）
    - scanpy.pl.scatter 通常返回 matplotlib Axes（有 .figure）
    - 也可能返回 list[Axes] 或带 .fig 的对象
    """
    import matplotlib.pyplot as plt

    out_dir = os.path.dirname(out_path)
    ensure_dir(out_dir)

    saved = False

    # Case 1: 对象本身有 savefig（如 seaborn.FacetGrid）
    if hasattr(obj, "savefig"):
        try:
            obj.savefig(out_path, bbox_inches="tight", dpi=150)
            saved = True
        except Exception:
            saved = False

    # Case 2: matplotlib Axes（具有 .figure）
    if not saved and hasattr(obj, "figure"):
        try:
            obj.figure.savefig(out_path, bbox_inches="tight", dpi=150)
            saved = True
        except Exception:
            saved = False

    # Case 3: 具有 .fig（某些对象提供）
    if not saved and hasattr(obj, "fig"):
        try:
            obj.fig.savefig(out_path, bbox_inches="tight", dpi=150)
            saved = True
        except Exception:
            saved = False

    # Case 4: list/tuple of Axes
    if not saved and isinstance(obj, (list, tuple)) and len(obj) > 0 and hasattr(obj[0], "figure"):
        try:
            obj[0].figure.savefig(out_path, bbox_inches="tight", dpi=150)
            saved = True
        except Exception:
            saved = False

    # Fallback：当前激活的全局 Figure
    if not saved:
        plt.gcf().savefig(out_path, bbox_inches="tight", dpi=150)

    # 关闭资源（不同对象关闭方式不同，不能使用 clf()）
    try:
        if hasattr(obj, "fig"):
            plt.close(obj.fig)
        elif hasattr(obj, "figure"):
            plt.close(obj.figure)
        else:
            plt.close("all")
    except Exception:
        plt.close("all")


# =============================================================================
# 高变基因过滤：按基因名模式从 HVG 中剔除不希望驱动聚类的基因
# =============================================================================

def refine_hvgs_by_gene_patterns(
    adata: ad.AnnData,
    mt_pattern: str = r"^MT-",
    rp_pattern: str = r"^RP[SL]",
    ncrna_pattern: str = r"^[A-Z][A-Z][0-9].*\.[0-9]",
    linc_pattern: str = r"(^LOC|LINC)[1-9]*",
    hvg_col: str = "highly_variable",
) -> Dict[str, List[str]]:
    """
    在 adata.var[hvg_col] 已经标记初步 HVG 的基础上，根据基因名模式剔除某些类型的基因，
    以避免这些基因在降维 / 聚类中占主导。

    参数
    ----------
    adata : ad.AnnData
        输入 AnnData 对象，要求 adata.var[hvg_col] 已经存在（bool）。
    mt_pattern : str
        线粒体基因的正则表达式模式（默认匹配类似 "MT-" 开头的基因）。
    rp_pattern : str
        核糖体蛋白基因的正则表达式模式（默认匹配 "RPS"/"RPL" 开头）。
    ncrna_pattern : str
        ncRNA 基因名的正则表达式模式。
    linc_pattern : str
        LINC/LOC 基因名的正则表达式模式。
    hvg_col : str
        存放 HVG 标记的列名（默认 "highly_variable"）。

    返回
    ----------
    removed_genes : Dict[str, List[str]]
        一个字典，记录每一类被剔除的基因列表，键包括：
        - "mt_genes"
        - "rp_genes"
        - "ncRNA_genes"
        - "LINC_genes"
        同时可以用 sum(len(v) for v in removed_genes.values()) 查看总剔除数。

    使用说明
    ----------
    - 通常流程：
        1) 先运行 sc.pp.highly_variable_genes(...) 得到初步 HVG；
        2) 再调用本函数剔除 MT/RP/ncRNA/LINC 等类型；
        3) 之后再根据 adata.var[hvg_col] 进行子集用于 Scale/PCA/聚类。
    """
    # 若没有 hvg_col 列，直接返回空
    if hvg_col not in adata.var.columns:
        return {
            "mt_genes": [],
            "rp_genes": [],
            "ncRNA_genes": [],
            "LINC_genes": [],
        }

    # 只在初步 HVG 中进行进一步剔除
    hvg_mask = adata.var[hvg_col].astype(bool).values
    var_names = adata.var_names

    # 线粒体基因
    mt_genes = list(var_names[hvg_mask & var_names.str.match(mt_pattern)])
    # 核糖体蛋白基因
    rp_genes = list(var_names[hvg_mask & var_names.str.match(rp_pattern)])
    # ncRNA
    ncRNA_genes = list(var_names[hvg_mask & var_names.str.match(ncrna_pattern)])
    # LINC/LOC
    LINC_genes = list(var_names[hvg_mask & var_names.str.match(linc_pattern)])

    remove_genes = set(mt_genes + rp_genes + ncRNA_genes + LINC_genes)

    # 当前 HVG 列表
    hvg_list = var_names[hvg_mask].tolist()
    # 过滤后的 HVG 列表
    hvg_keep = [g for g in hvg_list if g not in remove_genes]

    # 更新 adata.var[hvg_col]
    is_keep = np.isin(var_names, hvg_keep)
    adata.var[hvg_col] = is_keep

    removed = {
        "mt_genes": mt_genes,
        "rp_genes": rp_genes,
        "ncRNA_genes": ncRNA_genes,
        "LINC_genes": LINC_genes,
    }
    return removed


# =============================================================================
# 单样本层面的 doublet 检测与去除（基于 scrublet）
# =============================================================================

def run_doublet_detection_with_scrublet(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.06,
    sim_doublet_ratio: float = 2.0,
    n_prin_comps: int = 30,
    doublet_score_key: str = "doublet_score",
    doublet_label_key: str = "doublet_call",
    random_state: int = 0,
) -> Tuple[ad.AnnData, Dict[str, float]]:
    """
    使用 scrublet 对单个样本的 AnnData 进行 doublet 检测，并返回去除 doublet 后的新 AnnData。

    设计思路：
    ----------
    - scrublet 原始实现基于原始 counts 矩阵，因此本函数默认直接使用 adata.X
      （假设此时 adata.X 仍为计数矩阵；若已归一化，建议在归一化前调用本函数）。
    - 函数会在 adata.obs 中新增：
        * doublet_score_key : float，doublet 打分
        * doublet_label_key : bool，True 表示预测为 doublet
    - 返回的新 AnnData 已经去除了预测 doublet 的细胞，且保留上述 obs 列，方便后续追踪。

    参数
    ----------
    adata : ad.AnnData
        输入 AnnData 对象，最好为单样本 QC 后数据。
    expected_doublet_rate : float
        预期 doublet 比例，默认 0.06（6%）。
    sim_doublet_ratio : float
        模拟 doublet 与观测 singlet 的比例，详见 scrublet 文档。
    n_prin_comps : int
        预留参数：原本用于控制 scrublet 内部 PCA 主成分数量，当前版本的 Scrublet
        构造函数不再接受该参数，因此本函数暂不传递该值，仅保留以便未来扩展。
    doublet_score_key : str
        在 obs 中保存 doublet score 的列名。
    doublet_label_key : str
        在 obs 中保存 doublet 调用结果的列名（bool）。
    random_state : int
        随机种子，保证可复现。

    返回
    ----------
    adata_filtered : ad.AnnData
        去除 doublet 后的 AnnData（深拷贝）。
    stats : Dict[str, float]
        简单统计信息，包括：
        - "n_cells_before"
        - "n_cells_after"
        - "n_doublets"
        - "doublet_rate_observed"
    """
    try:
        import scrublet as scr
    except ImportError as e:
        raise ImportError(
            "需要安装 scrublet 才能运行 run_doublet_detection_with_scrublet，"
            "请使用 `pip install scrublet` 或 `conda install -c bioconda scrublet` 安装。"
        ) from e

    # 确保使用稀疏矩阵或 numpy array 均可
    counts_matrix = adata.X.copy()

    # 初始化 scrublet 对象
    # 注意：当前安装的 scrublet 版本 (0.2.3) 的 Scrublet.__init__ 不再接受 n_prin_comps 参数，
    # 因此这里只传入 counts_matrix / expected_doublet_rate / sim_doublet_ratio / random_state 等通用参数。
    scrub = scr.Scrublet(
        counts_matrix,
        expected_doublet_rate=expected_doublet_rate,
        sim_doublet_ratio=sim_doublet_ratio,
        random_state=random_state,
    )

    # 运行 doublet 检测（使用当前版本的默认参数设置）
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # 将结果写入 obs
    adata.obs[doublet_score_key] = doublet_scores
    adata.obs[doublet_label_key] = predicted_doublets.astype(bool)

    # 统计
    n_cells_before = adata.n_obs
    n_doublets = int(predicted_doublets.sum())
    n_cells_after = int(n_cells_before - n_doublets)
    doublet_rate_observed = n_doublets / float(n_cells_before) if n_cells_before > 0 else 0.0

    stats = {
        "n_cells_before": float(n_cells_before),
        "n_cells_after": float(n_cells_after),
        "n_doublets": float(n_doublets),
        "doublet_rate_observed": float(doublet_rate_observed),
    }

    # 过滤掉 doublet，返回新对象
    keep_mask = ~adata.obs[doublet_label_key].astype(bool).values
    adata_filtered = adata[keep_mask, :].copy()

    return adata_filtered, stats
