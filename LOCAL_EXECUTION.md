# Phase 2 — Local Execution Guide

How to reproduce the Phase 2 analysis on a local machine.

## Prerequisites

- **Conda environment `humica`** with: scanpy 1.11.5, anndata 0.12.19,
  numpy 2.4.6, scipy 1.17.1, pandas 2.3.3, matplotlib 3.11.0, h5py 3.16.0,
  scikit-learn 1.9.0 (Python 3.11.15).
- **Master checkpoint** `checkpoints/HuMicAtlas_validated.h5ad` present
  (6.3 GB on disk; MD5 `ccc6c9d77cb46ab29f4329e2465d1cf6`). This is the fixed
  Phase 1 output; step 01 verifies its MD5 and will refuse to proceed on a
  mismatch.
- **~8.5 GB free RAM.** Step 01 loads the master once (read-only) to subset
  the 41 candidate genes; steps 02–05 use the 238 MB candidate copy and are
  light.
- **Writable numba cache.** The env mount is read-only, so scripts set
  `NUMBA_CACHE_DIR` to `data/external/.numba_cache` before importing scanpy.
  The orchestrator and `_env.py` do this automatically.

## Recreating the environment (if needed)

```bash
conda create -n humica python=3.11
conda activate humica
pip install "scanpy==1.11.5" "anndata==0.12.19" "numpy==2.4.6" \
            "scipy==1.17.1" "pandas==2.3.3" "matplotlib==3.11.0" \
            "h5py==3.16.0" "scikit-learn==1.9.0"
```

## Running the full pipeline

From the project root:

```bash
conda activate humica
cd "<project root>/LSD insilico"
python scripts/06_run_all.py
```

This runs steps 01→05 in order, logging each step and its elapsed time, and
stops on the first failure. Total runtime is dominated by the one-time master
load in step 01 (a few minutes); steps 02–05 finish in seconds each.

## Partial / resume runs

```bash
python scripts/06_run_all.py --skip-load   # reuse existing candidate h5ad
python scripts/06_run_all.py --from 02      # resume from step 02
python scripts/06_run_all.py --only 05      # regenerate figures only
```

Steps 02–05 depend only on `checkpoints/phase2_candidates.h5ad` (step 01's
output) and the CSV tables, so once step 01 has produced the candidate object
you can iterate on the analysis and figures without reloading the 6.3 GB
master.

## Individual steps

Each script is standalone and can be run directly (same env, same cwd):

```bash
python scripts/01_load_checkpoint.py
python scripts/02_expression_analysis.py
python scripts/03_state_assignment.py
python scripts/04_module_conservation.py
python scripts/05_generate_figures.py
```

## Outputs

Tables land in `tables/`, figures in `figures/` (PNG + PDF). See
`PIPELINE_MANIFEST.md` for the full file-by-file inventory and
`docs/phase2_methods.md` for methodological detail.

## Troubleshooting

- **`numba` cache / permission error on scanpy import** — ensure
  `NUMBA_CACHE_DIR` points at a writable path (the scripts set it; if running a
  bare `python -c "import scanpy"`, export it yourself).
- **MD5 mismatch in step 01** — the master checkpoint differs from the fixed
  Phase 1 output; do not overwrite it. Restore the correct file before
  proceeding.
- **`MemoryError` in step 01** — close other applications; the master load
  needs the full object in RAM once. Steps 02–05 do not.
