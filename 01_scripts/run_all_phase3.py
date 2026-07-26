#!/usr/bin/env python3
"""
Phase 3 orchestrator — run the full NPC disease-kinetics pipeline (07→12)
=========================================================================

Runs every Phase-3 step in order, in one logged, standalone command:

  07a  Acquire NPC datasets from GEO (GSE221609 snRNA-seq, GSE152158 bulk)
       + build the human→mouse ortholog map.
  07   Preprocess GSE221609 snRNA-seq & annotate microglia (2,656 cells).
  07b  Preprocess GSE152158 bulk real-time course (31 disease-course samples).
  08   Diffusion pseudotime + QC GATE with a HARD STOP. This step decides,
       from metrics (tables/qc_gate_metrics.csv) not arbitrary cutoffs, which
       axis is the PRIMARY disease axis.
  09   Project the 8 Phase-2 STRING modules onto the disease axis.
  10   Module kinetics — ruptures PELT change-point on the disease-deviation
       signal Δ(week)=NPC−WT, bootstrap CIs, goodness-of-fit.
  11   Stable vs Dynamic remodelling (module-level cross-species test).
  12   Integrated tables + publication figures P3-H1..P3-H7.

QC-GATE / AXIS-SELECTION BRANCH (recorded outcome)
--------------------------------------------------
Step 08's gate is authoritative. In the executed run:
  * the single-cell diffusion-pseudotime axis FAILED the disease-relevance
    gate (DPT-vs-DAM Spearman=0.069; Npc1_KO not shifted, p=1.0), and
  * the bulk real-time axis PASSED (DAM-vs-week Spearman=0.746),
so the PRIMARY DISEASE AXIS = bulk_realtime and steps 10–11 run their onset
analysis on the bulk WEEK axis. The pseudotime axis is retained only as a
FAILING SUPPORTING axis. After step 08 this orchestrator re-reads
tables/qc_gate_metrics.csv and PRINTS which axis the gate selected, so the
branch is visible in the run log rather than silently assumed. The downstream
scripts already read the same gate table, so the pipeline stays correct
whichever way the gate falls.

Usage
-----
    python scripts/run_all_phase3.py                 # full pipeline
    python scripts/run_all_phase3.py --skip-acquire  # reuse cached raw data
    python scripts/run_all_phase3.py --from 09       # resume at a step
    python scripts/run_all_phase3.py --dry-run       # print plan only

Each step is run as a subprocess with the repo root as CWD; stdout/stderr are
streamed and also captured to logs/run_all_phase3.log. A non-zero exit from
any step aborts the pipeline (fail fast). Re-running is safe: every step
regenerates its own outputs.
"""
import argparse
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LOGS_DIR = os.path.join(ROOT, "logs")
TABLES_DIR = os.path.join(ROOT, "tables")
QC_GATE_CSV = os.path.join(TABLES_DIR, "qc_gate_metrics.csv")
RUN_LOG = os.path.join(LOGS_DIR, "run_all_phase3.log")

# Ordered pipeline. (step_id, script, one-line description).
STEPS = [
    ("07a", "07a_acquire_npc.py",     "Acquire NPC datasets from GEO + ortholog map"),
    ("07",  "07_prepare_npc.py",      "Preprocess GSE221609 snRNA-seq & annotate microglia"),
    ("07b", "07b_prepare_bulk.py",    "Preprocess GSE152158 bulk real-time course"),
    ("08",  "08_pseudotime.py",       "Diffusion pseudotime + QC gate (HARD STOP / axis choice)"),
    ("09",  "09_project_modules.py",  "Project Phase-2 STRING modules onto disease axis"),
    ("10",  "10_module_kinetics.py",  "Module kinetics: PELT change-point, CIs, goodness-of-fit"),
    ("11",  "11_stability_remodel.py","Stable vs Dynamic remodelling (module-level)"),
    ("12",  "12_generate_figures.py", "Integrated tables + figures P3-H1..P3-H7"),
]
STEP_IDS = [s[0] for s in STEPS]


def _log(msg, fh=None):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} [run_all] {msg}"
    print(line, flush=True)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def read_axis_decision():
    """Read tables/qc_gate_metrics.csv and derive which axis the gate selected.

    Returns (primary_axis, detail_dict) or (None, {}) if the table is absent.
    The rule mirrors step 08: the single-cell pseudotime axis is usable only if
    it tracks the DAM gradient (|rho|>=0.30) AND Npc1_KO cells shift toward the
    DAM end (p<0.05); otherwise the bulk real-time axis (if it passes B1/B2) is
    primary.
    """
    if not os.path.exists(QC_GATE_CSV):
        return None, {}
    m = {}
    with open(QC_GATE_CSV) as fh:
        for row in csv.DictReader(fh):
            m[row["metric"]] = row
    def _pass(key):
        return key in m and str(m[key]["pass"]).strip().lower() == "true"
    pseudotime_ok = _pass("P4_spearman_dpt_DAM") and _pass("P5_genotype_shift_p")
    bulk_ok = _pass("B1_bulk_timepoints") and _pass("B2_bulk_DAM_vs_week_rho")
    if pseudotime_ok:
        primary = "pseudotime"
    elif bulk_ok:
        primary = "bulk_realtime"
    else:
        primary = "NONE"
    return primary, {"pseudotime_ok": pseudotime_ok, "bulk_ok": bulk_ok, "metrics": m}


def run_step(step_id, script, desc, fh, dry_run=False):
    path = os.path.join(HERE, script)
    if not os.path.exists(path):
        _log(f"MISSING SCRIPT: {path}", fh)
        return 1
    _log(f"=== STEP {step_id}: {desc} ===", fh)
    _log(f"    $ python scripts/{script}", fh)
    if dry_run:
        return 0
    t0 = time.time()
    proc = subprocess.run([sys.executable, path], cwd=ROOT)
    dt = time.time() - t0
    if proc.returncode != 0:
        _log(f"    STEP {step_id} FAILED (exit {proc.returncode}) after {dt:.1f}s", fh)
    else:
        _log(f"    STEP {step_id} OK ({dt:.1f}s)", fh)
    return proc.returncode


def main():
    ap = argparse.ArgumentParser(description="Run the full Phase 3 pipeline (07->12).")
    ap.add_argument("--skip-acquire", action="store_true",
                    help="skip step 07a (reuse cached raw data in data/raw/npc/).")
    ap.add_argument("--from", dest="from_step", default=None, metavar="STEP",
                    help=f"resume at a step id (one of {', '.join(STEP_IDS)}).")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the plan and axis decision; run nothing.")
    args = ap.parse_args()

    os.makedirs(LOGS_DIR, exist_ok=True)
    steps = list(STEPS)
    if args.skip_acquire:
        steps = [s for s in steps if s[0] != "07a"]
    if args.from_step:
        if args.from_step not in STEP_IDS:
            ap.error(f"--from must be one of {STEP_IDS}")
        idx = [s[0] for s in steps].index(args.from_step) if args.from_step in [s[0] for s in steps] else 0
        steps = steps[idx:]

    with open(RUN_LOG, "a") as fh:
        _log("################ Phase 3 pipeline start ################", fh)
        _log(f"repo root: {ROOT}", fh)
        _log(f"steps to run: {', '.join(s[0] for s in steps)}", fh)

        for step_id, script, desc in steps:
            rc = run_step(step_id, script, desc, fh, dry_run=args.dry_run)
            if rc != 0 and not args.dry_run:
                _log(f"ABORTING pipeline at step {step_id}.", fh)
                sys.exit(rc)

            # After the QC gate, surface which axis it selected.
            if step_id == "08" and not args.dry_run:
                primary, detail = read_axis_decision()
                if primary is None:
                    _log("QC gate table not found after step 08 — cannot confirm axis.", fh)
                else:
                    _log(f"QC GATE OUTCOME: primary disease axis = {primary} "
                         f"(pseudotime_ok={detail['pseudotime_ok']}, "
                         f"bulk_ok={detail['bulk_ok']}).", fh)
                    if primary == "pseudotime":
                        _log("    -> single-cell pseudotime PASSED; steps 10-11 use the DPT axis.", fh)
                    elif primary == "bulk_realtime":
                        _log("    -> pseudotime FAILED; steps 10-11 use the BULK WEEK axis "
                             "(pseudotime retained as failing supporting axis).", fh)
                    else:
                        _log("    -> NEITHER axis passed. Onset claims are not defensible; "
                             "downstream kinetics should be treated as exploratory only.", fh)

        # Dry-run: still show the recorded axis decision if the table exists.
        if args.dry_run:
            primary, detail = read_axis_decision()
            if primary is not None:
                _log(f"[dry-run] recorded QC gate outcome: primary axis = {primary}", fh)

        _log("################ Phase 3 pipeline complete ################", fh)
        _log("Deliverables: tables/*.csv, figures/P3-H1..P3-H7.{png,pdf}, "
             "docs/phase3_*.md. STOP before Phase 4 — await approval.", fh)


if __name__ == "__main__":
    main()
