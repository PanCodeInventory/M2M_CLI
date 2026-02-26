---
name: Matrix2Markers Runner
description: This skill should be used when the user asks to "run Matrix2Markers", "validate Matrix2Markers config", "generate h5ad and plots from config", "prepare h5ad for cellxgene", or "troubleshoot Matrix2Markers run failure".
version: 0.1.0
---

# Purpose

Execute Matrix2Markers as a config-driven workflow with stable operational steps.
Use `M2M_CLI/m2m.py` as the single entrypoint for init, validation, execution, and cellxgene export.

# Use Cases

- Run full pipeline after editing config.
- Validate config before submitting long jobs.
- Diagnose missing runtime dependencies or wrong paths.
- Export a cellxgene-friendly h5ad.

# Lightweight Agent Workflow

Use this exact 3-step flow:

1. Ask concise intake questions.
2. Update config based on answers.
3. Return copy-paste run commands.

## Step 1: Intake Questions (ask only these)

1. What is the Matrix2Markers pipeline root path?
2. Which sample names and 10x directories should be used?
3. Keep default QC thresholds or provide overrides (min_genes, max_mt_pct, min_cells_per_gene)?
4. Keep default processing params or provide overrides (n_top_genes, resolution, n_neighbors)?
5. Which grouping key for marker finding (`leiden`/`group`/`sample`)?
6. Do you need a cellxgene export file (`yes`/`no`)?

## Step 2: Config Update Rules

- Update `samples` from provided sample-name to path pairs.
- Update `sample_metadata` for each sample if time/group was provided.
- Apply only explicitly provided parameter overrides.
- Keep untouched keys unchanged.
- Write changes to `M2M_CLI/config.user.yaml`.

## Step 3: Command Output Contract

Always return these commands after config update:

```bash
export M2M_PIPELINE_ROOT="<pipeline_root>"
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py validate --config /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py doctor --config /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml --snakemake /home/user/miniforge3/envs/matrix2markers/bin/snakemake --python /home/user/miniforge3/envs/matrix2markers/bin/python
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py run --config /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml --cores 8 --snakemake /home/user/miniforge3/envs/matrix2markers/bin/snakemake
```

If user selected cellxgene export, append:

```bash
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py export-cellxgene --input <merged.processed.h5ad path> --output <merged.processed.cellxgene.h5ad path>
```

# Commands

```bash
python M2M_CLI/m2m.py init --out M2M_CLI/config.user.yaml
python M2M_CLI/m2m.py validate --config M2M_CLI/config.user.yaml
python M2M_CLI/m2m.py doctor --config M2M_CLI/config.user.yaml
python M2M_CLI/m2m.py run --config M2M_CLI/config.user.yaml --cores 8
python M2M_CLI/m2m.py export-cellxgene --input results/data/merged/merged.processed.h5ad --output results/data/merged/merged.processed.cellxgene.h5ad
```

# Failure Handling

- If `validate` fails, fix config keys and sample paths first.
- If `doctor` fails, resolve missing executables and rerun.
- If run fails, check `results/logs/` and rerun with `--rerun-incomplete`.
- If cellxgene cannot open file, run `export-cellxgene` and test the exported file.

# Output Contract

- Main processed object: `results/data/merged/merged.processed.h5ad`
- Optional cellxgene export: `results/data/merged/*.cellxgene.h5ad`
- Figures: `results/plots/`
- Tables: `results/tables/`
- Logs: `results/logs/`
