#!/usr/bin/env python3
"""Matrix2Markers wrapper CLI.

Config-driven commands for Snakemake execution and cellxgene-compatible export.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any


TOOL_DIR = Path(__file__).resolve().parent


def _resolve_pipeline_root() -> Path:
    env_root = os.environ.get("M2M_PIPELINE_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()

    project_root = TOOL_DIR.parent
    candidate = project_root / "Analysis" / "Matrix2Markers"
    if (candidate / "Snakefile").exists() and (candidate / "config.yaml").exists():
        return candidate

    return project_root


PIPELINE_ROOT = _resolve_pipeline_root()
DEFAULT_CONFIG = PIPELINE_ROOT / "config.yaml"
DEFAULT_SNAKEFILE = PIPELINE_ROOT / "Snakefile"


def _load_yaml(config_path: Path) -> dict[str, Any]:
    try:
        yaml = importlib.import_module("yaml")
    except ImportError as exc:
        raise RuntimeError("PyYAML is required. Install with: pip install pyyaml") from exc

    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError("Config must be a YAML mapping at top level.")
    return data


def _validate_config_struct(cfg: dict[str, Any], config_path: Path) -> list[str]:
    errors: list[str] = []
    required_top = ["samples", "directories", "filtering", "processing", "find_markers", "gene_patterns"]
    for key in required_top:
        if key not in cfg:
            errors.append(f"Missing top-level key: {key}")

    samples = cfg.get("samples", {})
    if not isinstance(samples, dict) or not samples:
        errors.append("'samples' must be a non-empty mapping.")
    else:
        for sample, path in samples.items():
            if not isinstance(path, str):
                errors.append(f"samples.{sample} must be a string path")
                continue
            p = Path(path)
            if not p.exists():
                errors.append(f"samples.{sample} path not found: {p}")

    dirs = cfg.get("directories", {})
    required_dirs = ["raw_h5ad", "qc_h5ad", "merged", "plots", "tables", "logs"]
    if not isinstance(dirs, dict):
        errors.append("'directories' must be a mapping.")
    else:
        for key in required_dirs:
            if key not in dirs:
                errors.append(f"Missing directories.{key}")

    if errors:
        errors.insert(0, f"Config validation failed: {config_path}")
    return errors


def _run(cmd: list[str], cwd: Path) -> int:
    print("[CMD]", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=str(cwd), check=False)
    return proc.returncode


def _sanitize_for_cellxgene(adata: Any) -> None:
    """Apply compatibility cleanup for cellxgene readers."""
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if "log1p" in adata.uns and isinstance(adata.uns["log1p"], dict):
        if adata.uns["log1p"].get("base", "__MISSING__") is None:
            del adata.uns["log1p"]["base"]

    if "neighbors" in adata.uns and isinstance(adata.uns["neighbors"], dict):
        conn_key = adata.uns["neighbors"].get("connectivities_key")
        dist_key = adata.uns["neighbors"].get("distances_key")
        if (conn_key and conn_key not in adata.obsp) or (dist_key and dist_key not in adata.obsp):
            del adata.uns["neighbors"]


def cmd_init(args: argparse.Namespace) -> int:
    target = Path(args.out).resolve()
    if target.exists() and not args.force:
        print(f"Target exists: {target}. Use --force to overwrite.")
        return 2

    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(DEFAULT_CONFIG, target)
    print(f"Config template written: {target}")
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Config not found: {config_path}")
        return 2

    try:
        cfg = _load_yaml(config_path)
        errors = _validate_config_struct(cfg, config_path)
    except Exception as exc:  # noqa: BLE001
        if "PyYAML is required" not in str(exc):
            print(f"Validation failed: {exc}")
            return 2

        print(f"Validation fallback: {exc}")
        print("Running Snakemake dry-run as fallback validator...")
        snakemake_exec = getattr(args, "snakemake", "snakemake")
        cmd = [
            snakemake_exec,
            "--snakefile",
            str(DEFAULT_SNAKEFILE),
            "--configfile",
            str(config_path),
            "--directory",
            str(PIPELINE_ROOT),
            "--cores",
            "1",
            "-n",
        ]
        code = _run(cmd, cwd=PIPELINE_ROOT)
        if code == 0:
            print("Config looks runnable (Snakemake dry-run passed).")
            return 0
        return code

    if errors:
        for line in errors:
            print(f"[ERROR] {line}")
        return 2

    print(f"Config OK: {config_path}")
    return 0


def cmd_doctor(args: argparse.Namespace) -> int:
    config_path = Path(args.config).resolve()
    snakemake_path = shutil.which(args.snakemake)
    python_path = shutil.which(args.python)

    ok = True
    if snakemake_path:
        print(f"[OK] snakemake: {snakemake_path}")
    else:
        print(f"[FAIL] snakemake not found: {args.snakemake}")
        ok = False

    if python_path:
        print(f"[OK] python: {python_path}")
    else:
        print(f"[FAIL] python not found: {args.python}")
        ok = False

    if not DEFAULT_SNAKEFILE.exists():
        print(f"[FAIL] Snakefile not found: {DEFAULT_SNAKEFILE}")
        ok = False
    else:
        print(f"[OK] Snakefile: {DEFAULT_SNAKEFILE}")

    if not config_path.exists():
        print(f"[FAIL] config not found: {config_path}")
        ok = False
    else:
        print(f"[OK] config: {config_path}")
        validate_code = cmd_validate(
            argparse.Namespace(config=str(config_path), snakemake=args.snakemake)
        )
        ok = ok and (validate_code == 0)

    return 0 if ok else 2


def cmd_run(args: argparse.Namespace) -> int:
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"Config not found: {config_path}")
        return 2

    cmd = [
        args.snakemake,
        "--snakefile",
        str(DEFAULT_SNAKEFILE),
        "--configfile",
        str(config_path),
        "--directory",
        str(PIPELINE_ROOT),
        "--cores",
        str(args.cores),
    ]
    if args.use_conda:
        cmd.append("--use-conda")
    if args.rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.forceall:
        cmd.append("--forceall")
    if args.dry_run:
        cmd.append("-n")

    return _run(cmd, cwd=PIPELINE_ROOT)


def cmd_export_cellxgene(args: argparse.Namespace) -> int:
    in_path = Path(args.input).resolve()
    out_path = Path(args.output).resolve()

    if not in_path.exists():
        print(f"Input h5ad not found: {in_path}")
        return 2

    try:
        import importlib

        ad = importlib.import_module("anndata")
    except ImportError as exc:
        print("anndata is required. Install in runtime environment first.")
        print(f"Details: {exc}")
        return 2

    adata = ad.read_h5ad(in_path)

    _sanitize_for_cellxgene(adata)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(f"Exported cellxgene-compatible file: {out_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="m2m", description="Matrix2Markers wrapper CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init", help="Create config template")
    p_init.add_argument("--out", default="./config.user.yaml", help="Output config path")
    p_init.add_argument("--force", action="store_true", help="Overwrite existing file")
    p_init.set_defaults(func=cmd_init)

    p_validate = sub.add_parser("validate", help="Validate config fields and sample paths")
    p_validate.add_argument("--config", default=str(DEFAULT_CONFIG), help="Config YAML path")
    p_validate.set_defaults(func=cmd_validate)

    p_doctor = sub.add_parser("doctor", help="Check runtime prerequisites")
    p_doctor.add_argument("--config", default=str(DEFAULT_CONFIG), help="Config YAML path")
    p_doctor.add_argument("--snakemake", default="snakemake", help="Snakemake executable")
    p_doctor.add_argument("--python", default="python3", help="Python executable")
    p_doctor.set_defaults(func=cmd_doctor)

    p_run = sub.add_parser("run", help="Run Matrix2Markers workflow via Snakemake")
    p_run.add_argument("--config", default=str(DEFAULT_CONFIG), help="Config YAML path")
    p_run.add_argument("--cores", type=int, default=4, help="CPU cores")
    p_run.add_argument("--snakemake", default="snakemake", help="Snakemake executable")
    p_run.add_argument("--use-conda", action="store_true", help="Pass --use-conda")
    p_run.add_argument("--rerun-incomplete", action="store_true", help="Pass --rerun-incomplete")
    p_run.add_argument("--forceall", action="store_true", help="Force all jobs")
    p_run.add_argument("--dry-run", action="store_true", help="Dry run only")
    p_run.set_defaults(func=cmd_run)

    p_export = sub.add_parser("export-cellxgene", help="Write a cellxgene-friendly h5ad")
    p_export.add_argument("--input", required=True, help="Input h5ad path")
    p_export.add_argument("--output", required=True, help="Output h5ad path")
    p_export.set_defaults(func=cmd_export_cellxgene)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
