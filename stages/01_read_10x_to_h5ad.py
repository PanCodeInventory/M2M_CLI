#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
01_read_10x_to_h5ad.py

功能：
- 从 10x 输出目录（matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz）读取小鼠单细胞数据
- 为每个样本构建 AnnData 对象，添加 sample / group / time 等注释
- 计算基础 QC 指标（线粒体比例等）
- 将结果保存为单样本 raw h5ad，供后续 QC / 合并使用

运行约定（**避免路径堆叠的关键**）：
- 始终在项目根目录运行本脚本：
    /home/harry/NASH/scRNA_CD45/Scanpy
- 命令示例：
    python3 1_DataProcess/Matrix2H5AD_SC/scripts/stages/01_read_10x_to_h5ad.py

注意：
- 本脚本内部所有路径都基于固定的 PROJECT_ROOT，不会再产生 1_DataProcess/1_DataProcess 之类的嵌套。
- 物种：小鼠；线粒体前缀：MT-（在 io_and_qc_utils 中约定）。
"""

import os
import sys
import argparse
import traceback
from datetime import datetime
import anndata as ad

# ========================
# Path Setup
# ========================
# Allow importing from utils/ relative to this script
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
UTILS_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "utils")
if UTILS_DIR not in sys.path:
    sys.path.append(UTILS_DIR)

from io_and_qc_utils import (
    load_10x_mouse,
    add_mouse_mt_qc,
    ensure_dir,
)

def log_timestamp() -> str:
    """Returns current timestamp string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def main():
    parser = argparse.ArgumentParser(description="Read 10x matrix and convert to raw h5ad.")
    parser.add_argument("--sample-name", required=True, help="Name of the sample")
    parser.add_argument("--input-dir", required=True, help="Directory containing 10x matrix files")
    parser.add_argument("--output-file", required=True, help="Output path for .raw.h5ad file")
    parser.add_argument("--sample-group", default="NA", help="Sample Group (e.g. MCD, NCD)")
    parser.add_argument("--sample-time", default="NA", help="Sample Timepoint (e.g. 1w, 0w)")
    
    args = parser.parse_args()
    
    sample_name = args.sample_name
    sample_dir = args.input_dir
    out_path = args.output_file
    group = args.sample_group
    time_point = args.sample_time

    print(f"[{log_timestamp()}] === Processing Sample: {sample_name} ===")
    print(f"[INFO] Input Dir: {sample_dir}")
    print(f"[INFO] Metadata: group={group}, time={time_point}")
    print(f"[INFO] Output File: {out_path}")

    if not os.path.isdir(sample_dir):
        print(f"[ERROR] Directory does not exist: {sample_dir}")
        sys.exit(1)

    try:
        # 1) Read 10x data
        print(f"[{log_timestamp()}] [INFO] Loading 10x data...")
        adata = load_10x_mouse(
            sample_dir=sample_dir, 
            sample_name=sample_name,
            group=group,
            time=time_point
        )
        print(f"[INFO] Loaded: cells={adata.n_obs}, genes={adata.n_vars}")

        # 2) Calculate QC
        print(f"[{log_timestamp()}] [INFO] Calculating QC metrics...")
        add_mouse_mt_qc(adata)

        # 3) Write output
        print(f"[{log_timestamp()}] [INFO] Writing h5ad to {out_path}...")
        ensure_dir(os.path.dirname(out_path))
        adata.write(out_path)
        print(f"[{log_timestamp()}] [INFO] Done.")

    except Exception as e:
        print(f"[ERROR] Failed to process {sample_name}: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
