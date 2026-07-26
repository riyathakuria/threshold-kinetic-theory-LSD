# HuMicA Phase 1 - Environment Report

_Generated: 2026-07-04T21:40:26_

## System

- **Python**: 3.11.15
- **Platform**: macOS-14.8.7-arm64-arm-64bit
- **Architecture**: arm64
- **CPU cores**: 8
- **Total RAM**: 8.59 GB
- **Available RAM**: 2.43 GB
- **Conda env**: humica

## Package versions

| Package | Version |
|---|---|
| scanpy | 1.11.5 |
| anndata | 0.12.19 |
| h5py | 3.16.0 |
| numba | 0.65.1 |
| scipy | 1.17.1 |
| pandas | 2.3.3 |
| numpy | 2.4.6 |
| scikit-learn | 1.9.0 |
| psutil | 7.2.2 |
| openpyxl | 3.1.5 |
| requests | 2.34.2 |
| leidenalg | 0.12.0 |

## Software issues encountered & resolved

1. **scanpy import failed with numba caching error** (`RuntimeError: cannot cache function '_is_constant_csr_rows': no locator available`). Cause: scanpy JIT functions are decorated with `@njit(cache=True)`, and numba tries to write compiled caches next to the source files inside the conda env, which is mounted read-only. **Fix**: set `NUMBA_CACHE_DIR` to a writable path *before* importing numba/scanpy. All scanpy-importing scripts in this pipeline set `os.environ['NUMBA_CACHE_DIR']` at the top (see `scripts/_env.py`).

2. **RAM constraint**: total RAM is 8.59 GB while the HuMicA v2 `.h5ad` is ~6.8 GB on disk. A full in-memory load would exceed available memory. **Fix**: all structural inspection is done in **backed mode** (`sc.read_h5ad(path, backed='r')`), which reads `.obs`/`.var`/`.uns`/`.obsm` metadata and matrix shapes without materializing the expression matrix in RAM. Gene-presence checks operate on `.var`/`.raw.var` index only.

3. **Base `python` env lacked the single-cell stack** (scanpy/anndata/h5py/numba absent). **Fix**: created a dedicated conda env `humica` (python 3.11) with the full stack pinned.
