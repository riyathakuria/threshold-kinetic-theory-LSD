#!/usr/bin/env python3
"""
Shared environment setup for the HuMicA Phase 1 pipeline.

Import this FIRST in any script that imports numba or scanpy:

    from _env import ROOT, setup_numba_cache
    setup_numba_cache()
    import scanpy as sc

`setup_numba_cache()` points NUMBA_CACHE_DIR at a writable project-local path,
which avoids the read-only-conda-env caching failure that otherwise breaks
`import scanpy`.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def setup_numba_cache():
    cache_dir = os.path.join(ROOT, "data", "external", ".numba_cache")
    os.makedirs(cache_dir, exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = cache_dir
    return cache_dir
