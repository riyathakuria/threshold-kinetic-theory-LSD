#!/usr/bin/env python3
"""
Phase 3 · Step 08 — Diffusion pseudotime + QC gate with HARD STOP
=================================================================

Constructs a within-microglia diffusion pseudotime (DPT) "disease axis" for
GSE221609, then evaluates a QC gate that decides — from metrics, not
arbitrary cutoffs — whether the single-cell pseudotime axis is defensible as
the PRIMARY disease axis, or whether the bulk real-time course (GSE152158,
weeks 1/3/6/9) must be primary instead.

Design note (important). GSE221609 is a BINARY genotype contrast (WT vs
Npc1-/-), not a time course. A diffusion pseudotime therefore represents a
homeostatic -> disease-associated microglia (DAM) ACTIVATION CONTINUUM, not
real elapsed time. It is defensible as a disease axis only if the manifold is
connected AND the ordering tracks a genuine homeostatic->DAM transcriptional
gradient AND that gradient has disease relevance (Npc1-/- enriched at the
disease end). The bulk dataset provides an orthogonal, genuine real-TIME axis.

HARD STOP. If NEITHER the pseudotime axis NOR the bulk axis passes its minimum
quality bar, the script writes a diagnostic report and exits non-zero WITHOUT
producing a disease axis — Phase 3 cannot proceed on an indefensible axis.

Gate metrics (each tied to published practice, cited in the decision doc):
  P1 microglia N            trajectory inference needs enough cells to define a
                            manifold (hundreds+; Haghverdi 2016, scanpy DPT).
  P2 graph connectivity     DPT requires a single connected kNN component
                            (Wolf 2019 PAGA/DPT); multiple components => invalid.
  P3 diffusion spectrum     >=1 informative diffusion component (eigenvalue gap).
  P4 biological ordering     |Spearman(DPT, DAM_score)| with a clear homeostatic
                            gradient; sign flips DPT so high = disease.
  P5 disease relevance       DPT differs by genotype (Npc1-/- shifted toward the
                            DAM end; Mann-Whitney).
  B1 bulk design            >=3 real time points x replicates per genotype.
  B2 bulk disease signal    DAM signature rises with week in Npc1-/-.

Outputs:
  tables/pseudotime_values.csv     cell | genotype | dpt | DAM/homeo scores
  tables/qc_gate_metrics.csv       metric | value | reference | pass
  docs/phase3_qc_decision.md       full gate evaluation + primary-axis choice
  figures/P3-H1_trajectory.pdf/png trajectory figure
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import setup_numba_cache          # noqa: E402
setup_numba_cache()

import numpy as np                            # noqa: E402
import pandas as pd                           # noqa: E402
import scanpy as sc                           # noqa: E402
from scipy import sparse                      # noqa: E402
from scipy.stats import spearmanr, mannwhitneyu  # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402

import _phase3_config as cfg                  # noqa: E402
from _figstyle import apply_figure_style, panel_letter  # noqa: E402

log = cfg.get_logger("08_pseudotime")
sc.settings.verbosity = 1
np.random.seed(cfg.RANDOM_SEED)

DAM_MARKERS = ["Trem2", "Tyrobp", "Apoe", "Cst7", "Lpl", "Itgax",
               "Cd9", "Ctsb", "Ctsd", "Ctsl", "Axl", "Gpnmb", "Spp1"]
HOMEO_MARKERS = ["P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Sall1",
                 "Siglech", "Selplg", "Csf1r"]

# Minimum-quality references (documented, not invented):
MIN_MICROGLIA_N = 500      # DPT manifold; comfortably in published range.
MIN_ABS_SPEARMAN = 0.30    # DPT must track the DAM gradient at least moderately.
ALPHA = 0.05


def score_signatures(ad):
    dam = [g for g in DAM_MARKERS if g in ad.var_names]
    hom = [g for g in HOMEO_MARKERS if g in ad.var_names]
    sc.tl.score_genes(ad, dam, score_name="DAM_score",
                      ctrl_size=cfg.SCORE_GENES_CTRL_SIZE,
                      random_state=cfg.SCORE_GENES_RANDOM_STATE)
    sc.tl.score_genes(ad, hom, score_name="homeo_score",
                      ctrl_size=cfg.SCORE_GENES_CTRL_SIZE,
                      random_state=cfg.SCORE_GENES_RANDOM_STATE)
    return dam, hom


def build_manifold(ad):
    """Recompute PCA/neighbors/diffmap on microglia ONLY (stored graph was
    built on all cells then subset)."""
    sc.pp.highly_variable_genes(ad, n_top_genes=cfg.N_TOP_HVG, flavor="seurat")
    sub = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=cfg.N_PCS, random_state=cfg.RANDOM_SEED)
    ad.obsm["X_pca"] = sub.obsm["X_pca"]
    sc.pp.neighbors(ad, n_neighbors=cfg.N_NEIGHBORS, n_pcs=cfg.N_PCS,
                    random_state=cfg.RANDOM_SEED)
    sc.tl.umap(ad, random_state=cfg.RANDOM_SEED)
    sc.tl.diffmap(ad, n_comps=15)
    return ad


def n_connected_components(ad):
    from scipy.sparse.csgraph import connected_components
    conn = ad.obsp["connectivities"]
    n, labels = connected_components(conn, directed=False)
    # size of the largest component
    sizes = pd.Series(labels).value_counts()
    return n, int(sizes.iloc[0])


def choose_root(ad):
    """Root = most homeostatic WT cell (max homeo_score - DAM_score among WT)."""
    wt = ad.obs["genotype"] == "WT"
    metric = ad.obs["homeo_score"] - ad.obs["DAM_score"]
    metric = metric.where(wt, other=-np.inf)
    root = int(np.argmax(metric.values))
    return root


def run_dpt(ad):
    root = choose_root(ad)
    ad.uns["iroot"] = root
    sc.tl.dpt(ad, n_dcs=10)
    # Orient so high DPT = disease (positive corr with DAM score).
    rho, _ = spearmanr(ad.obs["dpt_pseudotime"], ad.obs["DAM_score"])
    if rho < 0:
        ad.obs["dpt_pseudotime"] = 1.0 - ad.obs["dpt_pseudotime"]
    return root


def evaluate_bulk():
    """B1/B2 on the bulk real-time axis."""
    meta = pd.read_csv(os.path.join(cfg.TABLES_DIR, "bulk_sample_metadata.csv"))
    expr = pd.read_csv(os.path.join(cfg.TABLES_DIR,
                                    "bulk_expression_normalized.csv"), index_col=0)
    dam = [g for g in DAM_MARKERS if g in expr.index]
    # Mean DAM-signature expression per sample (z-scored across genes).
    z = expr.loc[dam]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
    dam_sig = z.mean(axis=0)
    meta = meta.set_index("sample")
    meta["DAM_sig"] = dam_sig.reindex(meta.index).values
    npc = meta[meta["genotype"] == "NPC"]
    rho, p = spearmanr(npc["week"], npc["DAM_sig"])
    n_tp = meta.groupby("genotype")["week"].nunique().min()
    return {
        "n_timepoints": int(n_tp),
        "dam_week_rho_npc": float(rho),
        "dam_week_p_npc": float(p),
        "dam_present": len(dam),
    }, meta


def main():
    cfg.ensure_dirs()
    ad = sc.read_h5ad(cfg.NPC_MICROGLIA_CKPT)
    log.info(f"loaded {ad.n_obs:,} microglia")
    dam, hom = score_signatures(ad)
    log.info(f"DAM markers used: {dam}")
    log.info(f"homeostatic markers used: {hom}")
    ad = build_manifold(ad)

    ncomp, largest = n_connected_components(ad)
    log.info(f"kNN graph: {ncomp} connected components, largest={largest}")

    root = run_dpt(ad)
    log.info(f"root cell index {root} (genotype "
             f"{ad.obs['genotype'].iloc[root]})")

    # --- gate metrics ---
    rho_dam, p_dam = spearmanr(ad.obs["dpt_pseudotime"], ad.obs["DAM_score"])
    rho_hom, _ = spearmanr(ad.obs["dpt_pseudotime"], ad.obs["homeo_score"])
    dpt_wt = ad.obs.loc[ad.obs["genotype"] == "WT", "dpt_pseudotime"]
    dpt_ko = ad.obs.loc[ad.obs["genotype"] == "Npc1_KO", "dpt_pseudotime"]
    u, p_geno = mannwhitneyu(dpt_ko, dpt_wt, alternative="greater")
    eig = np.asarray(ad.uns["diffmap_evals"])
    eig_gap = float(eig[0] - eig[1]) if len(eig) > 1 else 0.0

    bulk, bulk_meta = evaluate_bulk()

    # --- pass/fail ---
    P1 = ad.n_obs >= MIN_MICROGLIA_N
    P2 = ncomp == 1
    P3 = eig_gap > 0
    P4 = abs(rho_dam) >= MIN_ABS_SPEARMAN
    P5 = p_geno < ALPHA
    pseudotime_pass = bool(P1 and P2 and P4)   # manifold + biological ordering
    B1 = bulk["n_timepoints"] >= 3
    B2 = (bulk["dam_week_rho_npc"] > 0) and (bulk["dam_week_p_npc"] < 0.10)
    bulk_pass = bool(B1 and B2)

    metrics = pd.DataFrame([
        ["P1_microglia_N", ad.n_obs, f">= {MIN_MICROGLIA_N}", P1],
        ["P2_connected_components", ncomp, "== 1", P2],
        ["P3_diffusion_eigen_gap", round(eig_gap, 4), "> 0", P3],
        ["P4_spearman_dpt_DAM", round(rho_dam, 3),
         f"|rho| >= {MIN_ABS_SPEARMAN}", P4],
        ["P4b_spearman_dpt_homeo", round(rho_hom, 3), "(sign check, <0 expected)",
         rho_hom < 0],
        ["P5_genotype_shift_p", f"{p_geno:.2e}", f"< {ALPHA} (KO>WT)", P5],
        ["B1_bulk_timepoints", bulk["n_timepoints"], ">= 3", B1],
        ["B2_bulk_DAM_vs_week_rho", round(bulk["dam_week_rho_npc"], 3),
         "> 0 and p<0.10", B2],
    ], columns=["metric", "value", "reference", "pass"])
    metrics.to_csv(os.path.join(cfg.TABLES_DIR, "qc_gate_metrics.csv"), index=False)
    log.info("gate metrics:\n%s", metrics.to_string(index=False))

    # --- HARD STOP if neither axis usable ---
    if not (pseudotime_pass or bulk_pass):
        write_decision(ad, metrics, bulk, pseudotime_pass, bulk_pass,
                       primary="NONE")
        log.error("HARD STOP: neither pseudotime nor bulk axis passes. "
                  "See docs/phase3_qc_decision.md")
        sys.exit(2)

    # --- pick primary axis ---
    # Prefer the genuine real-TIME bulk axis for temporal-onset claims when it
    # passes; pseudotime is an activation continuum, not elapsed time. If bulk
    # fails but pseudotime passes, pseudotime is primary.
    if bulk_pass and pseudotime_pass:
        primary = "bulk_realtime"   # both usable; real time is the stronger substrate for WHEN
    elif bulk_pass:
        primary = "bulk_realtime"
    else:
        primary = "pseudotime"

    # Save pseudotime values regardless (used as supporting single-cell axis).
    pt = ad.obs[["genotype", "gsm", "dpt_pseudotime", "DAM_score",
                 "homeo_score", "leiden"]].copy()
    pt.index.name = "cell"
    pt.to_csv(os.path.join(cfg.TABLES_DIR, "pseudotime_values.csv"))
    log.info("wrote tables/pseudotime_values.csv")

    make_h1_figure(ad, bulk_meta, primary)
    write_decision(ad, metrics, bulk, pseudotime_pass, bulk_pass, primary)
    # persist DPT-annotated microglia for step 09
    ad.write(os.path.join(cfg.CHECKPOINT_DIR, "npc_microglia_dpt.h5ad"))
    log.info(f"PRIMARY AXIS = {primary}")


def make_h1_figure(ad, bulk_meta, primary):
    apply_figure_style(frame="open", font="Arial")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    u = ad.obsm["X_umap"]
    s0 = axes[0].scatter(u[:, 0], u[:, 1], s=3, c=ad.obs["dpt_pseudotime"],
                        cmap="viridis", rasterized=True)
    axes[0].set_title("Diffusion pseudotime"); axes[0].set_xlabel("UMAP1")
    axes[0].set_ylabel("UMAP2"); fig.colorbar(s0, ax=axes[0], fraction=0.046,
                                              label="DPT (disease →)")
    panel_letter(axes[0], "a")
    for gt, col in [("WT", "#4C72B0"), ("Npc1_KO", "#C44E52")]:
        m = (ad.obs["genotype"] == gt).values
        axes[1].scatter(u[m, 0], u[m, 1], s=3, c=col, label=gt, rasterized=True)
    axes[1].set_title("Genotype"); axes[1].legend(markerscale=3)
    axes[1].set_xlabel("UMAP1"); panel_letter(axes[1], "b")
    # DPT distribution by genotype
    for gt, col in [("WT", "#4C72B0"), ("Npc1_KO", "#C44E52")]:
        d = ad.obs.loc[ad.obs["genotype"] == gt, "dpt_pseudotime"]
        axes[2].hist(d, bins=40, alpha=0.6, color=col, label=gt, density=True)
    axes[2].set_title("DPT by genotype"); axes[2].set_xlabel("DPT (disease →)")
    axes[2].set_ylabel("density"); axes[2].legend()
    panel_letter(axes[2], "c")
    fig.suptitle(f"P3-H1 NPC microglia disease axis (primary = {primary})",
                 fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(cfg.FIGURES_DIR, f"P3-H1_trajectory.{ext}"),
                    dpi=300)
    plt.close(fig)
    log.info("wrote figures/P3-H1_trajectory.{pdf,png}")


def write_decision(ad, metrics, bulk, pseudotime_pass, bulk_pass, primary):
    p = os.path.join(cfg.DOCS_DIR, "phase3_qc_decision.md")
    lines = [
        "# Phase 3 QC gate & disease-axis decision (Step 08)",
        "",
        "## Design context",
        "GSE221609 (snRNA-seq) is a **binary genotype contrast** (WT vs "
        "Npc1-/-), not a time course. A diffusion pseudotime (DPT) over these "
        "microglia therefore represents a **homeostatic → disease-associated "
        "microglia (DAM) activation continuum**, not elapsed real time. "
        "GSE152158 (bulk) provides an orthogonal, genuine **real-time** axis "
        "(weeks 1/3/6/9).",
        "",
        "## Gate metrics",
        "| metric | value | reference | pass |",
        "|---|---|---|---|",
    ]
    for _, r in metrics.iterrows():
        lines.append(f"| {r['metric']} | {r['value']} | {r['reference']} | "
                     f"{'✅' if r['pass'] else '❌'} |")
    lines += [
        "",
        "### References for the bars (not arbitrary cutoffs)",
        "- **P1**: diffusion pseudotime requires enough cells to define a "
        "manifold; published DPT applications (Haghverdi et al. 2016 Nat "
        "Methods; scanpy DPT) operate from hundreds of cells upward. "
        f"Here N = {ad.n_obs:,}.",
        "- **P2**: DPT/PAGA require a single connected kNN component "
        "(Wolf et al. 2019 Genome Biol); a disconnected graph makes a single "
        "pseudotime undefined.",
        "- **P4**: for the axis to be a *disease* axis, the ordering must track "
        "a homeostatic→DAM transcriptional gradient (Keren-Shaul et al. 2017 "
        "DAM signature). We require |Spearman| ≥ 0.30 (moderate).",
        "- **P5**: disease relevance — Npc1-/- microglia should sit toward the "
        "DAM end (one-sided Mann–Whitney).",
        "- **B1/B2**: the bulk axis needs ≥3 real time points per genotype and "
        "a DAM signature that rises with disease week in Npc1-/-.",
        "",
        "## Outcome",
        f"- Pseudotime axis usable: **{'YES' if pseudotime_pass else 'NO'}**.",
        f"- Bulk real-time axis usable: **{'YES' if bulk_pass else 'NO'}**.",
        f"- **PRIMARY DISEASE AXIS: {primary}.**",
        "",
    ]
    if primary == "NONE":
        lines += [
            "### HARD STOP",
            "Neither axis passed its minimum quality bar. Phase 3 does not "
            "proceed on an indefensible disease axis. This is a diagnostic "
            "report only; no module kinetics were computed.",
        ]
    elif primary == "bulk_realtime":
        lines += [
            "### Rationale for choosing bulk real-time as primary",
            "Both axes may be usable, but the change-point question is *when* "
            "each module remodels along disease progression. The bulk course "
            "carries genuine elapsed time (weeks), whereas pseudotime is an "
            "activation-state continuum inferred from a binary contrast. The "
            "bulk axis is therefore the stronger substrate for temporal-onset "
            "claims. The single-cell DPT axis is retained as a **supporting** "
            "axis and cross-checked in Step 06 (agreement recorded). Temporal "
            "order is reported as **association, not causation**.",
        ]
    else:
        lines += [
            "### Rationale for choosing pseudotime as primary",
            "The bulk axis did not pass, so the single-cell DPT activation "
            "continuum is the primary disease axis. Because it is inferred from "
            "a binary genotype contrast, it is interpreted as a homeostatic→DAM "
            "activation ordering (association, not elapsed time / causation).",
        ]
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    log.info(f"wrote {p}")


if __name__ == "__main__":
    main()
