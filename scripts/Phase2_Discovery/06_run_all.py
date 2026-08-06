

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths & environment
# ---------------------------------------------------------------------------
SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent

# Writable numba cache before any scanpy import in child processes.
_NUMBA_CACHE = PROJECT_ROOT / "data" / "external" / ".numba_cache"
_NUMBA_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(_NUMBA_CACHE))

# Ordered pipeline: (step id, script filename, one-line description).
STEPS = [
    ("01", "01_load_checkpoint.py",     "Load master checkpoint & derive candidate object"),
    ("02", "02_expression_analysis.py", "Per-gene expression across states"),
    ("03", "03_state_assignment.py",    "Descriptive state assignment"),
    ("04", "04_module_conservation.py", "STRING-module conservation"),
    ("05", "05_generate_figures.py",    "Publication figures H1-H4"),
]


def _log(msg: str) -> None:
    print(f"[run_all] {msg}", flush=True)


def run_step(step_id: str, script: str, desc: str) -> float:
    """Run one pipeline step as a subprocess. Returns elapsed seconds; raises on failure."""
    path = SCRIPTS_DIR / script
    if not path.exists():
        raise FileNotFoundError(f"missing pipeline script: {path}")
    _log(f"START  step {step_id}  {desc}")
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, str(path)],
        cwd=str(PROJECT_ROOT),
        env=os.environ.copy(),
    )
    dt = time.time() - t0
    if result.returncode != 0:
        _log(f"FAILED step {step_id}  (exit {result.returncode}, {dt:.1f}s)")
        raise SystemExit(result.returncode)
    _log(f"OK     step {step_id}  ({dt:.1f}s)")
    return dt


def main() -> None:
    ap = argparse.ArgumentParser(description="Run the Phase 2 pipeline (steps 01-05).")
    ap.add_argument("--from", dest="start", default=None,
                    help="resume from this step id (e.g. 02)")
    ap.add_argument("--only", dest="only", default=None,
                    help="run only this step id (e.g. 05)")
    ap.add_argument("--skip-load", action="store_true",
                    help="skip step 01 (reuse existing candidate h5ad)")
    args = ap.parse_args()

    steps = list(STEPS)
    if args.only:
        steps = [s for s in steps if s[0] == args.only]
        if not steps:
            raise SystemExit(f"unknown step id: {args.only}")
    else:
        if args.skip_load:
            steps = [s for s in steps if s[0] != "01"]
        if args.start:
            steps = [s for s in steps if s[0] >= args.start]
            if not steps:
                raise SystemExit(f"no steps at/after: {args.start}")

    _log(f"pipeline: {' -> '.join(s[0] for s in steps)}")
    total = 0.0
    for step_id, script, desc in steps:
        total += run_step(step_id, script, desc)
    _log(f"DONE   {len(steps)} step(s) in {total:.1f}s")


if __name__ == "__main__":
    main()
