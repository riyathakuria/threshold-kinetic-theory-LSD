#!/usr/bin/env python3
"""
Phase 2 · Step 03 — Descriptive gene -> preferred-state assignment
==================================================================

Scale: SCT Pearson residuals (see _phase2_config). No specificity index
(Tau) is computed — the object carries no non-negative detection scale, and
per the Phase 2 decision, gene-level specificity is reported ONLY as simple,
auditable descriptors:

    preferred_state   = argmax over states of the mean SCT residual
    preferred_group   = functional group of that state
    top_mean          = residual at the preferred state
    second_mean       = residual at the 2nd-ranked state
    assignment_margin = top_mean - second_mean   (how clear the preference is)

We deliberately do NOT threshold or bin these margins into
'specific/broad' classes; small margins are reported as-is so the reader
can see that many housekeeping-like genes are broadly represented across
states. Biological interpretation is carried at the module level (step 04).

A per-state composition table summarises how many Stable vs Dynamic genes
prefer each state (feeds the module-centric discussion, not a claim about
state identity).

Inputs
------
  tables/phase2_expression_mean_matrix.csv   (from step 02)
  data/metadata/phase2_candidate_genes.csv

Outputs
-------
  tables/phase2_state_assignment.csv         (per-gene descriptive assignment)
  tables/phase2_state_composition.csv        (per-state Stable/Dynamic counts)
  logs/03_state_assignment.log
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT, setup_numba_cache          # noqa: E402
setup_numba_cache()
from _phase2_config import (                        # noqa: E402
    TABLES_DIR, METADATA_DIR, STATE_ORDER, STATE_GROUP,
    EXPRESSION_SCALE_LABEL, ensure_dirs, get_logger,
)

import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402

LOG = get_logger("03_state_assignment")


def main():
    ensure_dirs()
    mm = pd.read_csv(os.path.join(TABLES_DIR,
                     "phase2_expression_mean_matrix.csv"), index_col=0)
    genes = pd.read_csv(os.path.join(METADATA_DIR, "phase2_candidate_genes.csv"))
    meta = genes.set_index("symbol")

    states = [s for s in STATE_ORDER if s in mm.columns]
    mm = mm[states]

    rows = []
    for sym, r in mm.iterrows():
        ranked = r.sort_values(ascending=False)
        top_state = ranked.index[0]
        top_mean = float(ranked.iloc[0])
        second_state = ranked.index[1]
        second_mean = float(ranked.iloc[1])
        rows.append({
            "symbol": sym,
            "stability_status": meta.loc[sym, "stability_status"],
            "string_cluster_id": int(meta.loc[sym, "string_cluster_id"]),
            "preferred_state": top_state,
            "preferred_group": STATE_GROUP[top_state],
            "top_mean_residual": round(top_mean, 4),
            "second_state": second_state,
            "second_mean_residual": round(second_mean, 4),
            "assignment_margin": round(top_mean - second_mean, 4),
        })
    assign = pd.DataFrame(rows).sort_values(
        ["preferred_group", "preferred_state", "assignment_margin"],
        ascending=[True, True, False]).reset_index(drop=True)
    ap = os.path.join(TABLES_DIR, "phase2_state_assignment.csv")
    assign.to_csv(ap, index=False)
    LOG.info("Wrote per-gene descriptive assignment -> %s", ap)

    # ---- Per-state composition (Stable vs Dynamic) -----------------------
    comp = (assign.groupby(["preferred_state", "stability_status"])
            .size().unstack(fill_value=0))
    for col in ("Stable", "Dynamic"):
        if col not in comp.columns:
            comp[col] = 0
    comp = comp.reindex([s for s in STATE_ORDER if s in comp.index])
    comp["n_genes"] = comp["Stable"] + comp["Dynamic"]
    comp["preferred_group"] = [STATE_GROUP[s] for s in comp.index]
    cp = os.path.join(TABLES_DIR, "phase2_state_composition.csv")
    comp.to_csv(cp)
    LOG.info("Wrote per-state composition -> %s", cp)

    print(json.dumps({
        "expression_scale": EXPRESSION_SCALE_LABEL,
        "n_genes": int(len(assign)),
        "states_used_as_preferred": sorted(assign.preferred_state.unique().tolist()),
        "median_assignment_margin": float(assign.assignment_margin.median()),
        "n_small_margin_lt_0p1": int((assign.assignment_margin < 0.1).sum()),
    }, indent=2))


if __name__ == "__main__":
    main()
