#!/usr/bin/env python3
"""
Phase 4 orchestrator — orthogonal evidence integration & prioritization
=======================================================================

Phase 4 has TWO kinds of step, and this orchestrator is explicit about the
difference because it matters for reproducibility.

(1) EVIDENCE-LAYER ACQUISITION (network-interactive, access date 2026-07-19)
    The five external layers and the internal-evidence assembly were computed
    against live web resources:

      A1/A2/A3  Internal evidence   -> tables/internal_evidence_scores.csv
                                       checkpoints/phase4_internal.pkl
      E1  SMPD1 perturbation (null) -> tables/asm_perturbation_module_scores.csv
                                       checkpoints/phase4_checkpoint1_perturbation.pkl
      E2  ASM functional network    -> tables/asm_network_module_scores.csv
          (STRING+Reactome+KEGG)       checkpoints/phase4_checkpoint2_network.pkl
      E3  Regulatory (ChEA3 TFs)    -> tables/regulatory_module_scores.csv
                                       checkpoints/phase4_checkpoint3_regulatory.pkl
      E4  Brain expression (HPA)    -> tables/brain_expression_module_scores.csv
                                       checkpoints/phase4_checkpoint4_brain.pkl
      E5  Disease (OpenTargets+     -> tables/disease_association_module_scores.csv
          ClinVar)                     checkpoints/phase4_checkpoint5_disease.pkl

    These are NOT re-run here: the remote resources are versioned and rate
    limited, and re-querying at a later date would silently change results
    (e.g. Open Targets scores, ChEA3 libraries). The saved checkpoints and
    per-layer tables ARE the reproducible record of that acquisition, pinned
    to ACCESS_DATE in scripts/_phase4_config.py. Each layer's method, source,
    and any null/gate handling is documented in docs/checkpoint{1..5}_summary.md
    and docs/REPRODUCIBILITY_phase4.md.

(2) INTEGRATION & PRESENTATION (fully deterministic, offline)
    From the consolidated matrices onward, everything is a pure function of the
    saved tables and is re-run here, in order:

      14  Integrated prioritization (APPROVED framework: modules 0.60 internal /
          0.40 external; genes 0.45/0.40/0.15; min-max per layer) + sensitivity
          over 50/50, 60/40, 70/30. Reproduces every ranked table.
      15  Publication figures P4-1..P4-5 (reuses 14's score_modules()).

Usage
-----
    python scripts/run_all_phase4.py            # verify layers, run 14 -> 15
    python scripts/run_all_phase4.py --dry-run  # print plan + layer check only

The orchestrator first VERIFIES that every evidence-layer table and checkpoint
from stage (1) is present (fail fast with a clear message if a checkpoint is
missing), PRINTS the E1 null and E4 gate outcomes so the two non-standard layer
handlings are visible in the run log, then executes stage (2). stdout/stderr are
streamed and captured to logs/run_all_phase4.log.
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LOGS_DIR = os.path.join(ROOT, "logs")
TABLES_DIR = os.path.join(ROOT, "tables")
CKPT_DIR = os.path.join(ROOT, "checkpoints")
RUN_LOG = os.path.join(LOGS_DIR, "run_all_phase4.log")

# Evidence-layer artifacts that stage (2) depends on (verified, not recomputed).
REQUIRED_TABLES = [
    "internal_evidence_scores.csv",
    "asm_perturbation_module_scores.csv",
    "asm_network_module_scores.csv",
    "regulatory_module_scores.csv",
    "brain_expression_module_scores.csv",
    "disease_association_module_scores.csv",
    "consolidated_module_evidence_matrix.csv",
    "consolidated_gene_evidence_matrix.csv",
]
REQUIRED_CKPTS = [
    "phase4_internal.pkl",
    "phase4_checkpoint1_perturbation.pkl",
    "phase4_checkpoint2_network.pkl",
    "phase4_checkpoint3_regulatory.pkl",
    "phase4_checkpoint4_brain.pkl",
    "phase4_checkpoint5_disease.pkl",
]
# Stage (2): deterministic scripts, run in order.
STAGE2 = [
    ("14", "14_integrate_prioritize.py", "Integrated prioritization + sensitivity"),
    ("15", "15_phase4_figures.py", "Publication figures P4-1..P4-5"),
]


def _log(fh, msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n"); fh.flush()


def verify_layers(fh):
    missing = []
    for t in REQUIRED_TABLES:
        if not os.path.exists(os.path.join(TABLES_DIR, t)):
            missing.append(f"tables/{t}")
    for c in REQUIRED_CKPTS:
        if not os.path.exists(os.path.join(CKPT_DIR, c)):
            missing.append(f"checkpoints/{c}")
    if missing:
        _log(fh, "FATAL: evidence-layer artifacts missing (stage 1 not complete):")
        for m in missing:
            _log(fh, f"    - {m}")
        _log(fh, "Re-acquire these layers (see docs/REPRODUCIBILITY_phase4.md) before integrating.")
        sys.exit(2)
    _log(fh, f"Evidence layers verified: {len(REQUIRED_TABLES)} tables + "
             f"{len(REQUIRED_CKPTS)} checkpoints present.")


def report_special_layers(fh):
    """Make the two non-standard layer handlings visible in the run log."""
    import csv
    e1 = os.path.join(TABLES_DIR, "asm_perturbation_module_scores.csv")
    try:
        with open(e1) as f:
            rows = list(csv.DictReader(f))
        ov = sum(int(float(r["genes_in_sig"])) for r in rows) if rows else 0
        _log(fh, f"E1 SMPD1 perturbation: total module-gene overlaps with KO signature = {ov} "
                 f"-> {'NULL soft-fail (zero weight)' if ov == 0 else 'nonzero (check weighting)'}")
    except Exception as e:
        _log(fh, f"E1 report skipped ({e})")
    e4 = os.path.join(TABLES_DIR, "brain_expression_module_scores.csv")
    try:
        with open(e4) as f:
            rows = list(csv.DictReader(f))
        _log(fh, f"E4 brain expression: {len(rows)} modules scored -> used as ELIGIBILITY GATE "
                 f"(non-discriminating positive control, not weighted).")
    except Exception as e:
        _log(fh, f"E4 report skipped ({e})")


def run_step(fh, tag, script, desc, dry):
    path = os.path.join(HERE, script)
    _log(fh, f"STEP {tag}: {desc}  ({script})")
    if dry:
        return
    if not os.path.exists(path):
        _log(fh, f"FATAL: {script} not found."); sys.exit(2)
    t0 = time.time()
    p = subprocess.run([sys.executable, path], cwd=ROOT,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    fh.write(p.stdout); fh.flush()
    sys.stdout.write(p.stdout)
    if p.returncode != 0:
        _log(fh, f"FATAL: step {tag} exited {p.returncode}. Aborting."); sys.exit(p.returncode)
    _log(fh, f"STEP {tag} done in {time.time() - t0:.1f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    os.makedirs(LOGS_DIR, exist_ok=True)
    with open(RUN_LOG, "a") as fh:
        _log(fh, "=" * 66)
        _log(fh, "Phase 4 orchestrator — evidence integration & prioritization")
        _log(fh, "Stage 1 (network-interactive layers) = verified, not recomputed.")
        verify_layers(fh)
        report_special_layers(fh)
        _log(fh, "Stage 2 (deterministic integration + figures):")
        for tag, script, desc in STAGE2:
            run_step(fh, tag, script, desc, args.dry_run)
        _log(fh, "Phase 4 pipeline complete." if not args.dry_run
                 else "[dry-run] plan printed; nothing executed.")


if __name__ == "__main__":
    main()
