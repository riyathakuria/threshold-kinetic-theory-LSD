#!/usr/bin/env python3
"""
Phase 2 · Step 01 — Load master checkpoint & derive candidate object
====================================================================

import os
import sys
import json
import hashlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT, setup_numba_cache          # noqa: E402
setup_numba_cache()
from _phase2_config import (                        # noqa: E402
    MASTER_CHECKPOINT, CANDIDATE_CHECKPOINT, PRESENCE_TABLE, METADATA_DIR,
    STATE_LABELS, STATE_GROUP, STATE_KEY_MARKERS, STATE_CLUSTER_COL,
    ensure_dirs, get_logger,
)

import numpy as np                                  # noqa: E402
import pandas as pd                                 # noqa: E402
import anndata as ad                                # noqa: E402
import scanpy as sc                                 # noqa: E402

LOG = get_logger("01_load_checkpoint")


def md5(path, chunk=8 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ensure_dirs()

    if not os.path.exists(MASTER_CHECKPOINT):
        LOG.error("Master checkpoint not found: %s", MASTER_CHECKPOINT)
        sys.exit(1)

    LOG.info("Hashing master checkpoint (integrity guard, pre-read)...")
    md5_before = md5(MASTER_CHECKPOINT)
    LOG.info("  master MD5 (before) = %s", md5_before)

    # ---- Load candidate table (41 present) -------------------------------
    pres = pd.read_csv(PRESENCE_TABLE)
    present = pres[pres["Present_in_var_X"] == True].copy()   # noqa: E712
    LOG.info("Confirmed present candidates: %d (%s)",
             len(present), present["Stability_Status"].value_counts().to_dict())
    assert len(present) == 41, f"expected 41 present genes, got {len(present)}"

    # ---- Open master READ-ONLY (backed) ----------------------------------
    LOG.info("Opening master in backed='r' mode (read-only)...")
    adata = sc.read_h5ad(MASTER_CHECKPOINT, backed="r")
    LOG.info("  master shape: %d cells x %d genes", adata.n_obs, adata.n_vars)
    assert adata.shape == (90716, 3000), f"unexpected master shape {adata.shape}"
    assert STATE_CLUSTER_COL in adata.obs.columns, \
        f"{STATE_CLUSTER_COL} missing from obs"

    var_ids = adata.var_names.to_numpy()
    ens = present["Ensembl_ID"].to_numpy()
    found_mask = np.isin(ens, var_ids)
    if not found_mask.all():
        missing = present.loc[~found_mask, "Symbol_Used"].tolist()
        LOG.error("Candidate Ensembl IDs absent from master var: %s", missing)
        sys.exit(1)
    LOG.info("All 41 candidate Ensembl IDs located in master var index.")

    # ---- Derive candidate-only object (load only the 41 columns) ---------
    # Positional indexing on the backed object, then bring into memory.
    gene_order = [g for g in var_ids if g in set(ens)]      # keep var order
    LOG.info("Subsetting master to %d candidate genes (derived copy)...",
             len(gene_order))
    sub = adata[:, gene_order].to_memory()
    # Detach from the backing file so nothing can write through.
    sub = sub.copy()
    adata.file.close()
    # Drop .raw: the derived object is a 41-gene efficiency copy holding
    # .X (SCT-normalized) + the SCT layer. The full 3,000-gene .raw matrix
    # lives in the authoritative master and is not duplicated here. (Also
    # avoids anndata's reserved '_index' column conflict when writing raw.var.)
    if sub.raw is not None:
        sub.raw = None
        LOG.info("Dropped .raw from derived object (full matrix stays in master).")
    LOG.info("  derived object: %d cells x %d genes (in memory)",
             sub.n_obs, sub.n_vars)

    # ---- Attach readable gene metadata to DERIVED object -----------------
    meta = present.set_index("Ensembl_ID")
    sub.var["symbol"] = [meta.loc[g, "Symbol_Used"] for g in sub.var_names]
    sub.var["stability_status"] = [meta.loc[g, "Stability_Status"] for g in sub.var_names]
    sub.var["string_cluster_id"] = [int(meta.loc[g, "STRING_Cluster_ID"]) for g in sub.var_names]
    sub.var["string_cluster_name"] = [meta.loc[g, "STRING_Cluster_Name"] for g in sub.var_names]
    sub.var["ensembl_id"] = sub.var_names

    # ---- Attach published state labels to DERIVED object -----------------
    raw_clusters = sub.obs[STATE_CLUSTER_COL].to_numpy()
    cl_int = pd.to_numeric(pd.Series(raw_clusters), errors="coerce").astype("Int64")
    state_name = cl_int.map(STATE_LABELS)
    n_unmapped = int(state_name.isna().sum())
    if n_unmapped:
        LOG.warning("%d cells have unmapped cluster ids -> 'Cluster N' fallback",
                    n_unmapped)
        fb = cl_int.astype("float").astype("Int64").astype(str).radd("Cluster ")
        state_name = state_name.fillna(fb)
    sub.obs["microglial_state"] = pd.Categorical(state_name.astype(str))
    sub.obs["state_group"] = sub.obs["microglial_state"].map(STATE_GROUP).astype("category")
    LOG.info("State composition:\n%s",
             sub.obs["microglial_state"].value_counts().to_string())

    # ---- Provenance in .uns ---------------------------------------------
    sub.uns["phase2_provenance"] = {
        "derived_from": os.path.basename(MASTER_CHECKPOINT),
        "master_md5": md5_before,
        "n_candidate_genes": int(sub.n_vars),
        "state_cluster_col": STATE_CLUSTER_COL,
        "note": ("Derived reduction of the authoritative master checkpoint "
                 "for computational efficiency. Master preserved unchanged."),
    }

    # ---- Write derived object -------------------------------------------
    LOG.info("Writing derived candidate object -> %s", CANDIDATE_CHECKPOINT)
    sub.write_h5ad(CANDIDATE_CHECKPOINT)

    # ---- State label map CSV --------------------------------------------
    slm = pd.DataFrame({
        "V1_cluster": list(STATE_LABELS.keys()),
        "microglial_state": list(STATE_LABELS.values()),
        "state_group": [STATE_GROUP[v] for v in STATE_LABELS.values()],
        "key_markers": [STATE_KEY_MARKERS[v] for v in STATE_LABELS.values()],
    })
    slm_path = os.path.join(METADATA_DIR, "state_label_map.csv")
    slm.to_csv(slm_path, index=False)
    LOG.info("Wrote state label map -> %s", slm_path)

    # ---- Candidate gene table -------------------------------------------
    genes_out = present[["Symbol_Used", "Ensembl_ID", "Stability_Status",
                         "STRING_Cluster_ID", "STRING_Cluster_Name"]].copy()
    genes_out.columns = ["symbol", "ensembl_id", "stability_status",
                         "string_cluster_id", "string_cluster_name"]
    genes_path = os.path.join(METADATA_DIR, "phase2_candidate_genes.csv")
    genes_out.to_csv(genes_path, index=False)
    LOG.info("Wrote candidate gene table -> %s", genes_path)

    # ---- Integrity guard: master unchanged ------------------------------
    LOG.info("Re-hashing master checkpoint (integrity guard, post-read)...")
    md5_after = md5(MASTER_CHECKPOINT)
    if md5_after != md5_before:
        LOG.error("MASTER CHECKPOINT MODIFIED! %s != %s", md5_after, md5_before)
        sys.exit(2)
    LOG.info("  master MD5 (after)  = %s  [UNCHANGED ✓]", md5_after)

    print(json.dumps({
        "master_shape": [90716, 3000],
        "candidate_shape": [int(sub.n_obs), int(sub.n_vars)],
        "master_md5_stable": md5_before == md5_after,
        "n_states": int(sub.obs["microglial_state"].nunique()),
        "candidate_checkpoint": CANDIDATE_CHECKPOINT,
    }, indent=2))


if __name__ == "__main__":
    main()
