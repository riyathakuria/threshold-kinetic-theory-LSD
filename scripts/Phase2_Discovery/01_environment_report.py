#!/usr/bin/env python3
"""
Phase 1 - Step 3: Software validation & environment report.

import os, sys, platform, datetime
import importlib.metadata as ilm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")

PKGS = ["scanpy", "anndata", "h5py", "numba", "scipy", "pandas",
        "numpy", "scikit-learn", "psutil", "openpyxl", "requests", "leidenalg"]


def main():
    os.makedirs(DOCS, exist_ok=True)
    import psutil
    vm = psutil.virtual_memory()

    versions = {}
    for p in PKGS:
        try:
            versions[p] = ilm.version(p)
        except Exception:
            versions[p] = "NOT FOUND"

    lines = []
    lines.append("# HuMicA Phase 1 - Environment Report\n")
    lines.append(f"_Generated: {datetime.datetime.now().isoformat(timespec='seconds')}_\n")
    lines.append("## System\n")
    lines.append(f"- **Python**: {sys.version.split()[0]}")
    lines.append(f"- **Platform**: {platform.platform()}")
    lines.append(f"- **Architecture**: {platform.machine()}")
    lines.append(f"- **CPU cores**: {psutil.cpu_count()}")
    lines.append(f"- **Total RAM**: {vm.total/1e9:.2f} GB")
    lines.append(f"- **Available RAM**: {vm.available/1e9:.2f} GB")
    lines.append(f"- **Conda env**: humica\n")

    lines.append("## Package versions\n")
    lines.append("| Package | Version |")
    lines.append("|---|---|")
    for p in PKGS:
        lines.append(f"| {p} | {versions[p]} |")
    lines.append("")

    lines.append("## Software issues encountered & resolved\n")
    lines.append(
        "1. **scanpy import failed with numba caching error** "
        "(`RuntimeError: cannot cache function '_is_constant_csr_rows': no locator available`). "
        "Cause: scanpy JIT functions are decorated with `@njit(cache=True)`, and numba tries to "
        "write compiled caches next to the source files inside the conda env, which is mounted "
        "read-only. **Fix**: set `NUMBA_CACHE_DIR` to a writable path *before* importing numba/scanpy. "
        "All scanpy-importing scripts in this pipeline set "
        "`os.environ['NUMBA_CACHE_DIR']` at the top (see `scripts/_env.py`).\n"
    )
    lines.append(
        f"2. **RAM constraint**: total RAM is {vm.total/1e9:.2f} GB while the HuMicA v2 `.h5ad` is "
        "~6.8 GB on disk. A full in-memory load would exceed available memory. **Fix**: all "
        "structural inspection is done in **backed mode** (`sc.read_h5ad(path, backed='r')`), which "
        "reads `.obs`/`.var`/`.uns`/`.obsm` metadata and matrix shapes without materializing the "
        "expression matrix in RAM. Gene-presence checks operate on `.var`/`.raw.var` index only.\n"
    )
    lines.append(
        "3. **Base `python` env lacked the single-cell stack** (scanpy/anndata/h5py/numba absent). "
        "**Fix**: created a dedicated conda env `humica` (python 3.11) with the full stack pinned.\n"
    )

    out = os.path.join(DOCS, "environment_report.md")
    with open(out, "w") as f:
        f.write("\n".join(lines))
    print("Wrote", out)


if __name__ == "__main__":
    main()
