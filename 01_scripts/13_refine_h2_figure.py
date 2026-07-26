#!/usr/bin/env python3
"""
13_refine_h2_figure.py — publication refinement of Figure P3-H2
===============================================================

STRICTLY a presentation/organization refinement of the EXISTING Phase-3
results. This script performs NO new analysis: it reads the already-computed
`tables/module_bulk_timecourse.csv` (the Δ(week)=NPC−WT signal) and
`tables/module_kinetics.csv` (PELT onsets, bootstrap CIs, patterns) and only
re-organizes / re-styles the P3-H2 small-multiples.

What changes vs the original P3-H2:
  * Panels are grouped by BIOLOGICAL remodelling behaviour, classified from the
    observed Δ trajectory (not by STRING module number, and not by the
    kinetics table's single-sign `direction` field which can hide biphasic
    behaviour):
      Group A  Upregulated / sustained activation
      Group B  Downregulated / suppressed
      Group D  Biphasic (early activation then suppression)
      Group C  No statistically supported change-point (remodels == False)
  * Within each group, modules are ordered by earliest PELT onset, then by peak
    week, then by peak |Δ| (documented sort key).
  * Section headings above each group; colourblind-friendly palette
    (up=red, down=blue, biphasic=purple, none=grey); identical panel widths;
    shared y-axis within each group.

The underlying numbers (Δ values, onset weeks, CIs, patterns) are used exactly
as stored; nothing quantitative is recomputed.

Biphasic rule (data-driven, peer-review-defensible)
---------------------------------------------------
A remodelling module is BIPHASIC iff, AFTER its primary peak week (week of
max |Δ|), the trajectory reaches an excursion of OPPOSITE sign whose magnitude
is >= 0.30 x |peak Δ|. Requiring the reversal to occur *after* the peak keeps a
pre-onset baseline blip (e.g. M4's small positive week-1 point preceding a
downward dip) from being mis-called biphasic, while catching a genuine
activation-then-suppression (M1: +0.407 at wk3 -> -0.195 at wk6).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

import _phase3_config as cfg

WEEKS = [1, 3, 6, 9]
MODULE_IDS = cfg.MODULE_IDS if hasattr(cfg, "MODULE_IDS") else [1, 4, 5, 6, 7, 8, 9, 10]

# Colourblind-friendly palette (seaborn 'muted' family; red/blue/purple/grey).
C_UP    = "#C44E52"   # upregulated / sustained activation
C_DOWN  = "#4C72B0"   # downregulated / suppressed
C_BIPH  = "#8172B3"   # biphasic
C_NONE  = "#9A9A9A"   # no supported remodelling
GROUP_COLOR = {"A": C_UP, "B": C_DOWN, "D": C_BIPH, "C": C_NONE}
GROUP_TITLE = {
    "A": "Group A — Upregulated / sustained activation",
    "B": "Group B — Downregulated / suppressed",
    "D": "Group D — Biphasic remodelling (early activation \u2192 suppression)",
    "C": "Group C — No statistically supported change-point",
}
# Fixed, shared y-limits per group (identical scaling within each group).
GROUP_YLIM = {"A": (-0.15, 1.35), "B": (-0.60, 0.35),
              "D": (-0.35, 0.55), "C": (-0.55, 2.80)}
BIPHASIC_REVERSAL_FRAC = 0.30


def _set_style():
    for fam in ("Arial", "Helvetica", "DejaVu Sans"):
        if any(fam in f.name for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = fam
            break
    plt.rcParams.update({
        "font.size": 11, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "figure.dpi": 120,
    })


def delta_series(bulk, cid):
    col = f"M{cid}"
    wt = bulk[bulk.genotype == "WT"].groupby("week")[col].mean().reindex(WEEKS)
    npc = bulk[bulk.genotype == "NPC"].groupby("week")[col].mean().reindex(WEEKS)
    return (npc - wt).values


def classify(cid, d, k):
    """Return group id ('A','B','C','D') from the observed Δ trajectory + PELT."""
    if not bool(k["remodels"]):
        return "C"
    peak_i = int(np.argmax(np.abs(d)))
    peak_val = d[peak_i]
    after = d[peak_i + 1:]
    if peak_val >= 0:
        opp = after[after < 0]
        rev = (-opp.min()) if opp.size else 0.0
    else:
        opp = after[after > 0]
        rev = (opp.max()) if opp.size else 0.0
    if rev >= BIPHASIC_REVERSAL_FRAC * abs(peak_val):
        return "D"
    return "A" if peak_val > 0 else "B"


def build_layout(groups_order, group_members):
    """Compute figure size + per-panel axis rectangles (identical widths).

    Vertical budget per band (top->bottom): section heading, gap for the panel's
    own 2-line title, the axes, then the x-axis label + inter-band gap. Gaps are
    sized so nothing collides.
    """
    ncol = max(len(v) for v in group_members.values())
    nband = len(groups_order)
    left, right = 0.075, 0.985
    col_gap = 0.020
    panel_w = (right - left - (ncol - 1) * col_gap) / ncol

    top_title = 0.912          # first band heading sits below the title block
    bottom = 0.045
    band_h = (top_title - bottom) / nband
    head_drop = 0.022          # heading baseline below band top
    title_gap = 0.100          # room below heading for rule + the panel's 2-line title
    xlabel_gap = 0.050         # room below axes for the x label + spacing
    panel_h = band_h - (title_gap + xlabel_gap)

    layout = {}
    for b, g in enumerate(groups_order):
        band_top = top_title - b * band_h
        head_y = band_top - head_drop
        panel_top = band_top - title_gap
        panel_bot = panel_top - panel_h
        rects = []
        for j, _ in enumerate(group_members[g]):
            x = left + j * (panel_w + col_gap)
            rects.append([x, panel_bot, panel_w, panel_h])
        layout[g] = {"head_y": head_y, "head_x": left, "rects": rects,
                     "rule_x0": left, "rule_x1": right}
    return layout, (ncol, nband)


def make_refined_h2(bulk, kin, outdir):
    _set_style()
    info = {}
    for cid in MODULE_IDS:
        d = delta_series(bulk, cid)
        k = kin[kin.string_cluster_id == cid].iloc[0]
        g = classify(cid, d, k)
        info[cid] = {"d": d, "k": k, "group": g}

    group_members = {"A": [], "B": [], "D": [], "C": []}
    for cid, v in info.items():
        group_members[v["group"]].append(cid)

    def sort_key(cid):
        k = info[cid]["k"]
        onset = k["onset_week"] if not pd.isna(k["onset_week"]) else 99
        pk = k["peak_week"] if not pd.isna(k["peak_week"]) else 99
        return (onset, pk, -float(k["peak_abs_delta"]))

    for g in group_members:
        if g == "C":
            group_members[g].sort(key=lambda c: -float(info[c]["k"]["peak_abs_delta"]))
        else:
            group_members[g].sort(key=sort_key)

    groups_order = [g for g in ["A", "B", "D", "C"] if group_members[g]]

    layout, (ncol, nband) = build_layout(groups_order, group_members)
    fig_w = 2.35 * ncol + 0.9
    fig_h = 2.95 * nband + 1.1
    fig = plt.figure(figsize=(fig_w, fig_h))

    for g in groups_order:
        color = GROUP_COLOR[g]
        lay = layout[g]
        fig.text(lay["head_x"], lay["head_y"], GROUP_TITLE[g],
                 fontsize=13, fontweight="bold", color=color, ha="left", va="center")
        yr = lay["head_y"] - 0.014
        line = plt.Line2D([lay["rule_x0"], lay["rule_x1"]], [yr, yr],
                          color=color, lw=1.4, alpha=0.55,
                          transform=fig.transFigure, zorder=0)
        fig.add_artist(line)

        for cid, rect in zip(group_members[g], lay["rects"]):
            v = info[cid]; d = v["d"]; k = v["k"]
            ax = fig.add_axes(rect)
            ax.axhline(0, color="#cccccc", lw=0.7, zorder=1)
            if g in ("A", "B", "D") and not pd.isna(k["onset_week"]):
                ax.axvspan(k["onset_ci_lo"], k["onset_ci_hi"],
                           color=color, alpha=0.12, zorder=1)
                ax.axvline(k["onset_week"], color=color, ls="--", lw=1.1,
                           alpha=0.8, zorder=2)
            ax.plot(WEEKS, d, "-o", color=color, ms=5, lw=2.0, zorder=3)
            ax.set_ylim(*GROUP_YLIM[g])
            ax.set_xticks(WEEKS)
            ax.set_xlim(0.3, 9.7)
            nm = cfg.STRING_CLUSTER_NAMES[cid]
            nm = nm if len(nm) <= 30 else nm[:28] + "\u2026"
            if g == "C":
                extra = " (large, n.s.)" if cid == 8 else ""
                sub = f"no supported onset{extra}"
            elif g == "D":
                sub = "\u2191 wk3 \u2192 \u2193 wk6 (PELT wk%d)" % int(k["onset_week"])
            else:
                arrow = "\u2191" if g == "A" else "\u2193"
                sub = "%s onset wk%d, %s" % (arrow, int(k["onset_week"]), k["pattern"])
            ax.set_title(f"M{cid}  {nm}\n{sub}", fontsize=10.5, loc="left",
                         linespacing=1.35)
            if rect[0] == layout[g]["rects"][0][0]:
                ax.set_ylabel("\u0394 module score\n(NPC \u2212 WT)", fontsize=10)
            ax.set_xlabel("Disease week", fontsize=10)

    fig.suptitle("P3-H2  Disease-specific microglial module trajectories, grouped by remodelling behaviour",
                 x=0.075, ha="left", y=0.978, fontsize=15, fontweight="bold")
    fig.text(0.075, 0.955,
             "Modules are organized by remodeling behavior rather than STRING module number.",
             ha="left", va="center", fontsize=11, style="italic", color="#333333")
    fig.text(0.075, 0.938,
             "Dashed line = PELT onset change-point;  shaded band = 95% bootstrap CI on onset;  "
             "\u0394 = mean(NPC) \u2212 mean(WT) module score per disease week.",
             ha="left", va="center", fontsize=9.5, color="#666666")

    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, "P3-H2_module_activity.png")
    pdf = os.path.join(outdir, "P3-H2_module_activity.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return {cid: info[cid]["group"] for cid in MODULE_IDS}, group_members


def main():
    bulk = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_bulk_timecourse.csv"))
    kin = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_kinetics.csv"))
    assign, members = make_refined_h2(bulk, kin, cfg.FIGURES_DIR)
    print("group assignment:", assign)
    print("ordered members:", members)


if __name__ == "__main__":
    main()
