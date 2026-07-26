#!/usr/bin/env python3
"""
Phase 1 - Step 5: Locate all 109 candidate genes within the HuMicA v2 atlas.

RULES (from project spec):
  * The 109-gene candidate universe is FIXED. We do NOT reduce it. The goal is to
    record WHERE each gene lives (or that it is absent), never to intersect-and-drop.
  * We never "intersect with .X first". We inspect .X (var), .raw (raw.var) and the
    SCT layer's feature set independently, then classify each candidate.
  * Every candidate is classified into exactly one Status and given a Reason.

Atlas gene namespaces (established in Step 4) - all Ensembl gene IDs:
  * .var index  == 'features'      (3000 HVGs, RNA assay)
  * .var SCT_features               (3000, SCT assay variable features)
  * .raw/var/_index                 (3000)
  All three are the SAME 3000 Ensembl IDs in different orderings, so the atlas has a
  single searchable universe of 3000 highly-variable genes. We still check each slot
  separately and report per-slot presence for full transparency.

Inputs:
  data/raw/HuMicA_v2.0.0.h5ad                 (gene namespaces, via h5py - no matrix load)
  tables/candidate_mapping.csv                (109 genes -> Ensembl IDs, resolved)
  Raw data 2/Gene_Stability_Classification.csv (Stable/Dynamic labels)
  Raw data 2/Top10_PPI_clusters.xlsx           (STRING cluster assignments; 'Gene ' has trailing space)

Outputs:
  tables/candidate_presence.csv        (per-gene: slots present + Stable/Dynamic + STRING cluster)
  tables/candidate_mapping_audit.csv   (Original Symbol, Alias Used, Ensembl ID,
                                        Present in var, Present in raw, Status, Reason)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT

import h5py
import numpy as np
import pandas as pd

RAW = os.path.join(ROOT, "data", "raw", "HuMicA_v2.0.0.h5ad")
TABLES = os.path.join(ROOT, "tables")
RD = os.path.join(ROOT, "Raw data 2")


def decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in arr]


def main():
    # --- atlas gene namespaces (no matrix materialization) ---
    with h5py.File(RAW, "r") as f:
        var_features = decode(f["/var/features"][:])        # == var index
        var_sct = decode(f["/var/SCT_features"][:])
        raw_index = decode(f["/raw/var/_index"][:])
    set_var = set(var_features)
    set_sct = set(var_sct)
    set_raw = set(raw_index)

    # --- candidate mapping (resolved Ensembl IDs) ---
    cm = pd.read_csv(os.path.join(TABLES, "candidate_mapping.csv"))

    # --- Stable/Dynamic classification ---
    cls = pd.read_csv(os.path.join(RD, "Gene_Stability_Classification.csv"))
    # find the symbol + status columns robustly
    sym_col = next(c for c in cls.columns if c.strip().lower() in ("gene", "symbol", "gene_symbol", "gene symbol"))
    stat_col = next(c for c in cls.columns if "stab" in c.lower() and "status" in c.lower())
    cls_map = dict(zip(cls[sym_col].astype(str).str.strip(), cls[stat_col].astype(str).str.strip()))

    # --- STRING clusters (Top10_PPI_clusters.xlsx; 'Gene ' has trailing space) ---
    ppi = pd.read_excel(os.path.join(RD, "Top10_PPI_clusters.xlsx"))
    ppi.columns = [c.strip() for c in ppi.columns]
    gcol = next(c for c in ppi.columns if c.strip().lower() == "gene")
    cid_col = next((c for c in ppi.columns if "cluster" in c.lower() and "id" in c.lower()), None)
    cname_col = next((c for c in ppi.columns if "cluster" in c.lower() and "name" in c.lower()), None)
    ppi[gcol] = ppi[gcol].astype(str).str.strip()
    cluster_id_map, cluster_name_map = {}, {}
    for _, r in ppi.iterrows():
        g = r[gcol]
        cluster_id_map.setdefault(g, r[cid_col] if cid_col else "")
        cluster_name_map.setdefault(g, r[cname_col] if cname_col else "")

    # --- classify each candidate ---
    pres_rows, audit_rows = [], []
    for _, r in cm.iterrows():
        orig = str(r["Original_Symbol"])
        used = str(r["Symbol_Used"])
        alias = "" if pd.isna(r.get("Alias_Used")) else str(r.get("Alias_Used") or "")
        ens = str(r["Ensembl_ID"])

        in_var = ens in set_var
        in_sct = ens in set_sct
        in_raw = ens in set_raw

        # Status logic. In this atlas var==sct==raw, so present-in-any == present-in-all,
        # but we keep the general logic so the script is correct for any object.
        if in_var and in_raw:
            status = "Present in X and raw"
            reason = "Ensembl ID found in .var (X) and .raw.var (both HVG namespaces)."
        elif in_raw and not in_var:
            status = "Present only in raw"
            reason = "Ensembl ID found in .raw.var but not in .var (X)."
        elif in_var and not in_raw:
            status = "Present only in X"
            reason = "Ensembl ID found in .var (X) but not in .raw.var."
        else:
            status = "Not detected"
            reason = ("Ensembl ID not present in the atlas 3000-HVG universe "
                      "(.var/.SCT_features/.raw all identical). Gene was not retained "
                      "as a highly-variable feature in HuMicA v2.")

        stability = cls_map.get(orig, cls_map.get(used, "NOT_FOUND"))
        cid = cluster_id_map.get(orig, cluster_id_map.get(used, ""))
        cname = cluster_name_map.get(orig, cluster_name_map.get(used, ""))

        pres_rows.append({
            "Original_Symbol": orig,
            "Symbol_Used": used,
            "Alias_Used": alias,
            "Ensembl_ID": ens,
            "Present_in_var_X": in_var,
            "Present_in_SCT_features": in_sct,
            "Present_in_raw": in_raw,
            "Stability_Status": stability,
            "STRING_Cluster_ID": cid,
            "STRING_Cluster_Name": cname,
            "Status": status,
        })
        audit_rows.append({
            "Original Symbol": orig,
            "Alias Used": alias,
            "Ensembl ID": ens,
            "Present in var": in_var,
            "Present in raw": in_raw,
            "Status": status,
            "Reason": reason,
        })

    pres = pd.DataFrame(pres_rows)
    audit = pd.DataFrame(audit_rows)
    pres.to_csv(os.path.join(TABLES, "candidate_presence.csv"), index=False)
    audit.to_csv(os.path.join(TABLES, "candidate_mapping_audit.csv"), index=False)

    # --- console summary ---
    n = len(pres)
    n_present = int((pres["Present_in_var_X"] | pres["Present_in_raw"]).sum())
    n_absent = n - n_present
    print(f"candidates total          : {n}")
    print(f"present in atlas (var|raw): {n_present}")
    print(f"not detected              : {n_absent}")
    print("\nby Status:")
    print(pres["Status"].value_counts().to_string())
    print("\nStability x presence:")
    pres["_present"] = pres["Present_in_var_X"] | pres["Present_in_raw"]
    print(pres.groupby("Stability_Status")["_present"].agg(["size", "sum"]).to_string())
    print("\nStability labels found:", sorted(pres["Stability_Status"].unique()))
    miss = pres[pres["Stability_Status"] == "NOT_FOUND"]
    if len(miss):
        print("WARNING - genes with no Stable/Dynamic label:", list(miss["Original_Symbol"]))
    nocluster = pres[pres["STRING_Cluster_ID"].astype(str).isin(["", "nan"])]
    print(f"\ngenes with no STRING cluster: {len(nocluster)}")
    if len(nocluster):
        print(list(nocluster["Original_Symbol"]))


if __name__ == "__main__":
    main()
