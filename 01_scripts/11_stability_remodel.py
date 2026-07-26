#!/usr/bin/env python3
"""
Phase 3 · Step 11 — Stable vs Dynamic remodelling (MODULE-level test)
=====================================================================

Phase 2 classified each candidate gene as developmentally STABLE or DYNAMIC in
HUMAN microglia (NormFinder-based). Phase 3 asks a MODULE-level cross-species
question:

    Does a module's Phase-2 DYNAMIC character (measured in human microglia)
    predict how strongly that module REMODELS along the mouse NPC disease axis?

This is a conserved-program hypothesis: modules that are transcriptionally
plastic in human homeostatic microglia are predicted to be the ones that move
during mouse NPC disease. It is tested at the MODULE level (n = 8 STRING
modules), NOT gene-by-gene, and framed as human-module property → mouse-module
activity (never human state → mouse state).

Module Phase-2 character:
  dynamic_frac   = fraction of the module's candidate genes called Dynamic
  module_class   = Dynamic-dominant (dynamic_frac > 0.5) vs Stable-dominant
                   (matches Phase-2 Class A/B assignment)

Module NPC remodelling metric (from tables/module_kinetics.csv):
  peak_abs_delta          magnitude of NPC-vs-WT divergence
  bootstrap_frac_remodel  reliability of onset detection
  remodels                confident onset detected (bool)

Tests (small-n, reported honestly):
  * Spearman rho: dynamic_frac vs peak_abs_delta and vs bootstrap_frac_remodel
  * Mann-Whitney U: peak_abs_delta, Dynamic-dominant vs Stable-dominant modules
Interpretation is association, not causation, and n = 8 modules limits power;
this is reported as a consistency check on the conserved-program hypothesis,
not a definitive test.

Outputs:
  tables/stability_vs_pseudotime.csv   (per-module join of Phase-2 char + NPC kinetics)
  tables/stability_remodel_stats.csv   (the two tests)
  figures/P3-H4_stable_dynamic.pdf/png
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                               # noqa: E402
import pandas as pd                              # noqa: E402
from scipy import stats                          # noqa: E402
import matplotlib                                # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                  # noqa: E402

import _phase3_config as cfg                     # noqa: E402
from _figstyle import apply_figure_style         # noqa: E402

log = cfg.get_logger("11_stability_remodel")


def main():
    cfg.ensure_dirs()
    genes = pd.read_csv(os.path.join(cfg.METADATA_DIR,
                                     "phase2_candidate_genes.csv"))
    kin = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_kinetics.csv"))

    # Phase-2 module character
    comp = (genes.groupby("string_cluster_id")["stability_status"]
            .agg(n_genes="size",
                 n_dynamic=lambda s: int((s == "Dynamic").sum()))
            .reset_index())
    comp["dynamic_frac"] = comp["n_dynamic"] / comp["n_genes"]
    comp["module_class"] = np.where(comp["dynamic_frac"] > 0.5,
                                    "Dynamic-dominant", "Stable-dominant")

    m = comp.merge(kin, on="string_cluster_id", how="inner")
    m["module"] = "M" + m["string_cluster_id"].astype(str)
    keep = ["string_cluster_id", "module", "module_name", "n_genes",
            "n_dynamic", "dynamic_frac", "module_class", "onset_week",
            "onset_ci_lo", "onset_ci_hi", "peak_abs_delta",
            "bootstrap_frac_remodel", "remodels", "direction", "timing",
            "pattern", "pelt_r2"]
    m = m[keep].sort_values("dynamic_frac", ascending=False)
    m.to_csv(os.path.join(cfg.TABLES_DIR, "stability_vs_pseudotime.csv"),
             index=False)
    log.info("wrote tables/stability_vs_pseudotime.csv")

    # ---- tests ----
    rows = []
    for metric in ["peak_abs_delta", "bootstrap_frac_remodel"]:
        rho, p = stats.spearmanr(m["dynamic_frac"], m[metric])
        rows.append({"test": "spearman", "x": "dynamic_frac", "y": metric,
                     "n": len(m), "statistic": round(float(rho), 3),
                     "p_value": round(float(p), 4)})
        log.info("Spearman dynamic_frac vs %s: rho=%.3f p=%.4f (n=%d)",
                 metric, rho, p, len(m))

    dyn = m[m["module_class"] == "Dynamic-dominant"]["peak_abs_delta"]
    sta = m[m["module_class"] == "Stable-dominant"]["peak_abs_delta"]
    if len(dyn) and len(sta):
        u, p = stats.mannwhitneyu(dyn, sta, alternative="greater")
        rows.append({"test": "mannwhitney_u_greater",
                     "x": "Dynamic-dominant vs Stable-dominant",
                     "y": "peak_abs_delta",
                     "n": f"{len(dyn)}v{len(sta)}",
                     "statistic": round(float(u), 3),
                     "p_value": round(float(p), 4)})
        log.info("Mann-Whitney peak_abs_delta Dynamic>Stable: U=%.1f p=%.4f "
                 "(n=%dv%d)", u, p, len(dyn), len(sta))

    stat = pd.DataFrame(rows)
    stat.to_csv(os.path.join(cfg.TABLES_DIR, "stability_remodel_stats.csv"),
                index=False)
    log.info("wrote tables/stability_remodel_stats.csv")

    make_h4_figure(m, stat)
    return m, stat


def make_h4_figure(m, stat):
    apply_figure_style(frame="open", font="Arial")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.4))

    # panel a: dynamic_frac vs peak_abs_delta scatter, colored by class
    cmap = {"Dynamic-dominant": "#C44E52", "Stable-dominant": "#4C72B0"}
    for cls, sub in m.groupby("module_class"):
        ax1.scatter(sub["dynamic_frac"], sub["peak_abs_delta"], s=70,
                    color=cmap[cls], label=cls, zorder=3,
                    edgecolor="white", linewidth=0.5)
    for _, r in m.iterrows():
        ax1.annotate(r["module"], (r["dynamic_frac"], r["peak_abs_delta"]),
                     fontsize=6, xytext=(3, 3), textcoords="offset points")
    rho = stat[(stat.test == "spearman") &
               (stat.y == "peak_abs_delta")]["statistic"].iloc[0]
    pv = stat[(stat.test == "spearman") &
              (stat.y == "peak_abs_delta")]["p_value"].iloc[0]
    ax1.set_xlabel("Phase-2 Dynamic fraction (human microglia)")
    ax1.set_ylabel("NPC remodelling magnitude\n(peak |Δ NPC−WT|, mouse)")
    ax1.set_title(f"a  Conserved plasticity → remodelling\nSpearman ρ={rho:.2f}, p={pv:.3f} (n={len(m)})",
                  fontsize=8, loc="left")
    ax1.legend(fontsize=6, loc="upper left")

    # panel b: peak_abs_delta by class (strip + median)
    order = ["Stable-dominant", "Dynamic-dominant"]
    for i, cls in enumerate(order):
        vals = m[m["module_class"] == cls]["peak_abs_delta"].values
        x = np.random.default_rng(0).normal(i, 0.06, size=len(vals))
        ax2.scatter(x, vals, s=60, color=cmap[cls], zorder=3,
                    edgecolor="white", linewidth=0.5)
        ax2.hlines(np.median(vals), i - 0.2, i + 0.2, color=cmap[cls], lw=2)
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(order, fontsize=7)
    ax2.set_ylabel("peak |Δ NPC−WT|")
    prow = stat[stat.test == "mannwhitney_u_greater"]
    pv2 = prow["p_value"].iloc[0] if len(prow) else float("nan")
    ax2.set_title(f"b  Remodelling by module class\nMann–Whitney (Dyn>Stab) p={pv2:.3f}",
                  fontsize=8, loc="left")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(cfg.FIGURES_DIR,
                                 f"P3-H4_stable_dynamic.{ext}"), dpi=300)
    plt.close(fig)
    log.info("wrote figures/P3-H4_stable_dynamic.{pdf,png}")


if __name__ == "__main__":
    main()
