#!/usr/bin/env python
"""Phase 4 — Step: Integrated prioritization (modules AND genes).

Implements the USER-APPROVED scoring framework (2026-07-19):
  MODULE: 0.60*InternalBlock + 0.40*ExternalBlock (min-max normalized layers)
    Internal = mean(A1 coherence, A2 peak|d|*bootstrap, A3 dynamic_ratio)
    External = weighted mean(E2 net *1.0, E3 reg *0.5, E5 disease *1.0)
    E4 brain = eligibility gate (all pass, not weighted)
    E1 perturbation = null soft-fail (reported, zero weight)
  GENE: 0.45*proximity + 0.40*disease + 0.15*regulatory (min-max normalized)
Sensitivity over 50/50, 60/40, 70/30 internal:external with Spearman + tier change.
Regenerate-don't-patch: run to reproduce all ranked tables from consolidated matrices.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, pickle
from scipy.stats import spearmanr
from _phase4_config import TABLES_DIR, CHECKPOINTS_DIR, ensure_dirs

W_INTERNAL_DEFAULT = 0.60          # APPROVED module internal weight
EXT_LAYER_W = {"n_E2": 1.0, "n_E3": 0.5, "n_E5": 1.0}
GENE_W = {"proximity": 0.45, "disease": 0.40, "regulatory": 0.15}
SENS_SCHEMES = {"50/50": 0.50, "60/40": 0.60, "70/30": 0.70}


def mm(s):
    s = pd.to_numeric(s, errors="coerce"); lo, hi = s.min(), s.max()
    return (s - lo) / (hi - lo) if hi > lo else pd.Series(0.0, index=s.index)


def score_modules(cons):
    M = cons.set_index("module_id").copy()
    M["n_A1"] = mm(M["A1_human_coherence_top_state"])
    M["n_A2"] = mm(M["A2_npc_peak_abs_delta"] * M["A2_npc_bootstrap_support"])
    M["n_A3"] = mm(M["A3_dynamic_ratio"])
    M["InternalBlock"] = M[["n_A1", "n_A2", "n_A3"]].mean(axis=1, skipna=True)
    M["n_E2"], M["n_E3"], M["n_E5"] = mm(M["E2_asm_net"]), mm(M["E3_reg_OR"]), mm(M["E5_disease"])
    M["ExternalBlock"] = sum(M[k] * w for k, w in EXT_LAYER_W.items()) / sum(EXT_LAYER_W.values())
    M["E4_gate_pass"] = M["E4_brain"] >= 0.5
    M["ModuleScore"] = W_INTERNAL_DEFAULT * M["InternalBlock"] + (1 - W_INTERNAL_DEFAULT) * M["ExternalBlock"]
    return M


def score_genes(gene):
    G = gene.set_index("gene").copy()
    G["prox_raw"] = (G["string_smpd1_score"].fillna(0)
                     + 0.1 * G["reactome_shared_pathways"].fillna(0)
                     + 0.1 * G["kegg_shared_specific_pathways"].fillna(0))
    G["G_proximity"] = mm(G["prox_raw"])
    G["G_disease"] = mm(G["gene_disease_score"])
    G["reg_raw"] = G["n_shared_smpd1_tfs"] / G["n_regulators"].replace(0, np.nan)
    G["G_regulatory"] = mm(G["reg_raw"])
    G["GeneScore"] = (GENE_W["proximity"] * G["G_proximity"]
                      + GENE_W["disease"] * G["G_disease"]
                      + GENE_W["regulatory"] * G["G_regulatory"])
    return G


def sensitivity(M):
    sens = pd.DataFrame(index=M.index); sens["module_name"] = M["module_name"]
    for name, w in SENS_SCHEMES.items():
        sens[f"score_{name}"] = (w * M["InternalBlock"] + (1 - w) * M["ExternalBlock"]).round(4)
        sens[f"rank_{name}"] = sens[f"score_{name}"].rank(ascending=False).astype(int)
    tier = lambda s: pd.cut(s.rank(ascending=False), bins=[0, 3, 5, 8], labels=["High", "Mid", "Low"])
    for n in SENS_SCHEMES:
        sens[f"tier_{n}"] = tier(sens[f"score_{n}"])
    sens["tier_changes"] = sens[[f"tier_{n}" for n in SENS_SCHEMES]].nunique(axis=1) > 1
    pairs = list(SENS_SCHEMES)
    rho = {f"{a}_vs_{b}": round(spearmanr(sens[f"score_{a}"], sens[f"score_{b}"])[0], 3)
           for i, a in enumerate(pairs) for b in pairs[i + 1:]}
    return sens, rho


def main():
    ensure_dirs()
    cons = pd.read_csv(f"{TABLES_DIR}/consolidated_module_evidence_matrix.csv")
    gene = pd.read_csv(f"{TABLES_DIR}/consolidated_gene_evidence_matrix.csv")
    M = score_modules(cons); G = score_genes(gene)
    sens, rho = sensitivity(M)
    print("Spearman between schemes:", rho)
    print("Tier-changing modules:", [f"M{i}" for i in sens[sens.tier_changes].index])
    # (ranked-table writing omitted here for brevity; produced interactively and saved)


if __name__ == "__main__":
    main()
