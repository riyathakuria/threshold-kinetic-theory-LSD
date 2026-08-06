
import os
import sys
import json
import itertools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT, setup_numba_cache          # noqa: E402
setup_numba_cache()
from _phase2_config import (                        # noqa: E402
    CANDIDATE_CHECKPOINT, TABLES_DIR, STATE_ORDER, STRING_CLUSTER_NAMES,
    COHERENCE_METHOD, EXPRESSION_MATRIX, EXPRESSION_SCALE_LABEL,
    RANDOM_SEED, ensure_dirs, get_logger,
)

import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
import scanpy as sc                                 # noqa: E402
import scipy.sparse as sp                           # noqa: E402
from scipy.stats import rankdata                    # noqa: E402

LOG = get_logger("04_module_conservation")


def _dense(mat):
    return np.asarray(mat.todense()) if sp.issparse(mat) else np.asarray(mat)


def mean_pairwise_spearman(mat):
    """mat: cells x genes (>=2 genes). Mean of off-diagonal pairwise
    Spearman correlations (rank-transform then Pearson)."""
    if mat.shape[1] < 2:
        return np.nan, 0
    ranks = np.column_stack([rankdata(mat[:, j]) for j in range(mat.shape[1])])
    c = np.corrcoef(ranks, rowvar=False)
    iu = np.triu_indices_from(c, k=1)
    vals = c[iu]
    vals = vals[np.isfinite(vals)]
    return float(np.mean(vals)), int(len(vals))


def main():
    ensure_dirs()
    sc.settings.verbosity = 0
    LOG.info("Loading candidate object: %s", CANDIDATE_CHECKPOINT)
    adata = sc.read_h5ad(CANDIDATE_CHECKPOINT)
    if EXPRESSION_MATRIX != "X":
        adata.X = adata.layers[EXPRESSION_MATRIX]
    LOG.info("  %d cells x %d genes | scale=%s",
             adata.n_obs, adata.n_vars, EXPRESSION_SCALE_LABEL)

    states = [s for s in STATE_ORDER if s in set(adata.obs["microglial_state"])]
    clusters = sorted(adata.var["string_cluster_id"].unique().tolist())
    LOG.info("Present STRING modules: %s", clusters)

    # ---- 1. ACTIVITY: score_genes per module, mean per state ------------
    activity = pd.DataFrame(index=[f"cluster_{c}" for c in clusters],
                            columns=states, dtype=float)
    for c in clusters:
        members = adata.var_names[adata.var["string_cluster_id"] == c].tolist()
        score_key = f"score_cl{c}"
        if len(members) == 1:
            # score_genes needs a control pool; with 1 gene just use its residual
            gi = adata.var_names.get_loc(members[0])
            col = _dense(adata.X[:, gi]).ravel()
            adata.obs[score_key] = col
        else:
            sc.tl.score_genes(adata, members, score_name=score_key,
                              random_state=RANDOM_SEED, ctrl_size=50)
        for s in states:
            m = adata.obs["microglial_state"] == s
            activity.loc[f"cluster_{c}", s] = float(adata.obs.loc[m, score_key].mean())
    activity.to_csv(os.path.join(TABLES_DIR, "phase2_module_activity_by_state.csv"))
    LOG.info("Wrote module activity-by-state.")

    # ---- 2. COHERENCE: within-module mean pairwise Spearman -------------
    X = _dense(adata.X)
    state_vec = adata.obs["microglial_state"].to_numpy()
    coh_rows = []
    for c in clusters:
        cols = np.where((adata.var["string_cluster_id"] == c).to_numpy())[0]
        members = adata.var_names[cols].tolist()
        coh_all, npairs = mean_pairwise_spearman(X[:, cols])
        # coherence within the module's most-active state
        top_state = activity.loc[f"cluster_{c}"].astype(float).idxmax()
        m = state_vec == top_state
        coh_top, _ = mean_pairwise_spearman(X[np.ix_(m, cols)])
        coh_rows.append({
            "string_cluster_id": c,
            "string_cluster_name": STRING_CLUSTER_NAMES.get(c, f"(cluster {c})"),
            "n_genes": len(members),
            "n_pairs": npairs,
            "coherence_all_cells": None if np.isnan(coh_all) else round(coh_all, 4),
            "most_active_state": top_state,
            "coherence_in_top_state": None if (coh_top is None or np.isnan(coh_top))
                                       else round(coh_top, 4),
            "genes": ",".join(sorted(members)),
        })
    coh = pd.DataFrame(coh_rows)
    coh.to_csv(os.path.join(TABLES_DIR, "phase2_module_coherence.csv"), index=False)
    LOG.info("Wrote module coherence.")

    # ---- 3. COMPOSITION + combined conservation table -------------------
    comp = (adata.var.groupby("string_cluster_id")["stability_status"]
            .value_counts().unstack(fill_value=0))
    for col in ("Stable", "Dynamic"):
        if col not in comp.columns:
            comp[col] = 0
    comp = comp.reindex(clusters)
    comp["n_genes"] = comp["Stable"] + comp["Dynamic"]
    comp["pct_stable"] = (100 * comp["Stable"] / comp["n_genes"]).round(1)

    combined = coh.merge(
        comp[["Stable", "Dynamic", "pct_stable"]],
        left_on="string_cluster_id", right_index=True)
    # peak activity value at most-active state
    combined["peak_activity"] = [
        round(float(activity.loc[f"cluster_{c}", st]), 4)
        for c, st in zip(combined.string_cluster_id, combined.most_active_state)]
    combined = combined[[
        "string_cluster_id", "string_cluster_name", "n_genes",
        "Stable", "Dynamic", "pct_stable",
        "most_active_state", "peak_activity",
        "coherence_all_cells", "coherence_in_top_state", "n_pairs", "genes"]]
    combined.to_csv(os.path.join(TABLES_DIR, "phase2_module_conservation.csv"),
                    index=False)
    LOG.info("Wrote combined module conservation table.")

    print(json.dumps({
        "expression_scale": EXPRESSION_SCALE_LABEL,
        "coherence_method": COHERENCE_METHOD,
        "modules_present": clusters,
        "modules_absent": [2, 3],
        "activity_range": [float(activity.values.min()), float(activity.values.max())],
        "coherence_all_cells": {int(r.string_cluster_id): r.coherence_all_cells
                                for r in combined.itertuples()},
    }, indent=2))


if __name__ == "__main__":
    main()
