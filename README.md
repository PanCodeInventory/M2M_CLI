# Matrix2Markers Tool Wrapper

This directory packages `Matrix2Markers` as a config-driven CLI plus a reusable skill spec.

## Layout

- `m2m.py`: wrapper CLI
- `skill/SKILL.md`: skill definition draft

## Quick Start

Use a Python environment that has `snakemake` and `pyyaml`.

```bash
export M2M_PIPELINE_ROOT="/home/user/Pan Chongshi/Projects/NK_Expansion_sc/Analysis/Matrix2Markers"
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py init --out /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py validate --config /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py doctor --config /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml --snakemake /home/user/miniforge3/envs/matrix2markers/bin/snakemake --python /home/user/miniforge3/envs/matrix2markers/bin/python
python /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/m2m.py run --config /home/user/Pan\ Chongshi/Repository/Software/M2M_CLI/config.user.yaml --cores 8 --snakemake /home/user/miniforge3/envs/matrix2markers/bin/snakemake
```

## Agent Lightweight Mode

For AI-agent usage, keep interaction minimal:

1. Ask intake questions (sample paths, optional parameter overrides, groupby key, need cellxgene export).
2. Update only the relevant keys in `config.user.yaml`.
3. Return direct copy-paste commands for `validate`, `doctor`, `run`, and optional `export-cellxgene`.

## Cellxgene Export

```bash
python M2M_CLI/m2m.py export-cellxgene \
  --input results/data/merged/merged.processed.h5ad \
  --output results/data/merged/merged.processed.cellxgene.h5ad
```

The export command applies compatibility cleanup used during this session:

- remove `uns['log1p']['base']` when it is `None`
- ensure `obs_names` and `var_names` are strings and unique
- remove dangling `uns['neighbors']` if referenced graph keys are missing in `obsp`
