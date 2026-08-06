#!/usr/bin/env python3
Per-gene adult expression across microglial states

import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT, setup_numba_cache          # noqa: E402
setup_numba_cache()
from _phase2_config import (                        # noqa: E402
    CANDIDATE_CHECKPOINT, TABLES_DIR, STATE_ORDER,
    EXPRESSION_MATRIX, EXPRESSION_SCALE_LABEL,
    ensure_dirs, get_logger,
)

import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
import scanpy as sc                                 # noqa: E402
import scipy.sparse as sp                           # noqa: E402

LOG = get_logger("02_expression_analysis")


def main():
    ensure_dirs()
    LOG.info("Loading derived candidate object: %s", CANDIDATE_CHECKPOINT)
    adata = sc.read_h5ad(CANDIDATE_CHECKPOINT)
    LOG.info("  %d cells x %d genes | expression scale = %s (matrix=%s)",
             adata.n_obs, adata.n_vars, EXPRESSION_SCALE_LABEL, EXPRESSION_MATRIX)

    states = [s for s in STATE_ORDER if s in set(adata.obs["microglial_state"])]
    symbols = adata.var["symbol"].to_numpy()

    X = adata.X if EXPRESSION_MATRIX == "X" else adata.layers[EXPRESSION_MATRIX]
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)

    state_vec = adata.obs["microglial_state"].to_numpy()
    state_idx = {s: np.where(state_vec == s)[0] for s in states}
    n_by_state = {s: len(ix) for s, ix in state_idx.items()}
    LOG.info("Cells per state: %s", n_by_state)

    long_rows = []
    for j, sym in enumerate(symbols):
        col = X[:, j]
        for s in states:
            vals = col[state_idx[s]]
            long_rows.append({
                "symbol": sym,
                "ensembl_id": adata.var_names[j],
                "stability_status": adata.var["stability_status"].iloc[j],
                "string_cluster_id": int(adata.var["string_cluster_id"].iloc[j]),
                "microglial_state": s,
                "n_cells": len(vals),
                "mean_residual": float(vals.mean()),
                "median_residual": float(np.median(vals)),
                "sd_residual": float(vals.std(ddof=0)),
            })
    long = pd.DataFrame(long_rows)
    long_path = os.path.join(TABLES_DIR, "phase2_expression_by_state.csv")
    long.to_csv(long_path, index=False)
    LOG.info("Wrote long expression table -> %s (%d rows)", long_path, len(long))

    # ---- Wide mean-residual matrix (gene x state) ------------------------
    mean_mat = long.pivot(index="symbol", columns="microglial_state",
                          values="mean_residual")[states]
    mean_mat.to_csv(os.path.join(TABLES_DIR, "phase2_expression_mean_matrix.csv"))
    LOG.info("Wrote wide mean-residual matrix.")

    # ---- Lightweight descriptive per-gene summary ------------------------
    # Purely descriptive: preferred state = argmax mean residual; across-state
    # range = spread of per-state means. No Tau, no thresholds.
    summ_rows = []
    for sym in mean_mat.index:
        gmean = mean_mat.loc[sym]
        meta = long.loc[long.symbol == sym].iloc[0]
        summ_rows.append({
            "symbol": sym,
            "stability_status": meta["stability_status"],
            "string_cluster_id": int(meta["string_cluster_id"]),
            "preferred_state_max_residual": gmean.idxmax(),
            "max_state_mean_residual": float(gmean.max()),
            "min_state_mean_residual": float(gmean.min()),
            "across_state_residual_range": float(gmean.max() - gmean.min()),
        })
    summ = pd.DataFrame(summ_rows).sort_values(
        ["string_cluster_id", "symbol"]).reset_index(drop=True)
    summ_path = os.path.join(TABLES_DIR, "phase2_gene_expression_summary.csv")
    summ.to_csv(summ_path, index=False)
    LOG.info("Wrote per-gene descriptive summary -> %s", summ_path)

    print(json.dumps({
        "genes": int(len(mean_mat)),
        "states": states,
        "expression_scale": EXPRESSION_SCALE_LABEL,
        "residual_range_global": [float(mean_mat.values.min()),
                                  float(mean_mat.values.max())],
        "median_across_state_range": float(summ["across_state_residual_range"].median()),
    }, indent=2))


if __name__ == "__main__":
    main()
