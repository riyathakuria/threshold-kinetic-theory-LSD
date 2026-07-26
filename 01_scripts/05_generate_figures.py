#!/usr/bin/env python3
"""
Phase 2 · Step 05 — Publication figures H1-H4
=============================================

All figures are built from the step 02-04 CSV tables (no h5ad reload) on the
native SCT Pearson-residual scale. Arial, 300 dpi PNG + vector PDF.

H1  Gene x state mean-residual heatmap (41 genes, 9 states), rows grouped by
    STRING module with a stability side-strip. The granular expression map.
H2  Module x state activity heatmap (8 modules, 9 states) — module-level
    aggregate showing where each conserved program is most active.
H3  Module conservation properties: (a) within-module coherence (all cells vs
    most-active state), (b) Stable/Dynamic composition per module.
H4  Module-centric synthesis dashboard (the Phase 3 reference figure): one row
    per module aligning peak activity, preferred state, representative gene,
    and Stable/Dynamic composition.

Outputs: figures/H1_gene_state_heatmap.{png,pdf} ... H4_*.{png,pdf}
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT                               # noqa: E402
from _phase2_config import (                         # noqa: E402
    TABLES_DIR, FIGURES_DIR, METADATA_DIR, STATE_ORDER, STATE_GROUP,
    ensure_dirs, get_logger,
)

import numpy as np                                   # noqa: E402
import pandas as pd                                  # noqa: E402
import matplotlib as mpl                             # noqa: E402
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.patches import Rectangle             # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402

from _figstyle import (apply_figure_style, set_frame, panel_letter,  # noqa: E402,F401
                       panel_crops)

LOG = get_logger("05_generate_figures")

# --- shared display constants ---------------------------------------------
MOD_SHORT = {
    1: "PKA / glucagon", 4: "EGFR signaling", 5: "Rho signal transduction",
    6: "Translation elongation", 7: "Apoptosis", 8: "Blood microparticle",
    9: "Unfolded-protein response", 10: "Sphingolipid metabolism",
}
GROUP_COLORS = {"Homeostatic": "#4C72B0", "DAM": "#C44E52",
                "DIM": "#DD8452", "Macrophage": "#8172B3"}
STAB_COLORS = {"Stable": "#0072B2", "Dynamic": "#E69F00"}  # Okabe-Ito, CVD-safe
DIVERGING = "RdBu_r"


def _load():
    mean_mat = pd.read_csv(os.path.join(TABLES_DIR,
                           "phase2_expression_mean_matrix.csv"), index_col=0)[STATE_ORDER]
    act = pd.read_csv(os.path.join(TABLES_DIR,
                      "phase2_module_activity_by_state.csv"), index_col=0)[STATE_ORDER]
    con = pd.read_csv(os.path.join(TABLES_DIR, "phase2_module_conservation.csv"))
    genes = pd.read_csv(os.path.join(METADATA_DIR, "phase2_candidate_genes.csv"))
    return mean_mat, act, con, genes


def _save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGURES_DIR, f"{stem}.{ext}"),
                    dpi=300, bbox_inches="tight")
    LOG.info("Saved %s.{png,pdf}", stem)


def _state_group_spans(states):
    """Return [(group, start_idx, end_idx)] for contiguous state groups."""
    spans, i = [], 0
    while i < len(states):
        g = STATE_GROUP[states[i]]
        j = i
        while j < len(states) and STATE_GROUP[states[j]] == g:
            j += 1
        spans.append((g, i, j))
        i = j
    return spans


# ==========================================================================
# H1 — gene x state mean-residual heatmap
# ==========================================================================
def fig_h1(mean_mat, genes):
    meta = genes.set_index("symbol")
    # order genes by module then symbol
    order = (meta.reset_index()
             .sort_values(["string_cluster_id", "symbol"])["symbol"].tolist())
    M = mean_mat.loc[order]
    vmax = np.abs(M.values).max()

    fig = plt.figure(figsize=(7.4, 9.2))
    # explicit axes: heatmap (left), module-label band handled via text,
    # colorbar in its own far-right axis so nothing overlaps.
    ax = fig.add_axes([0.15, 0.09, 0.50, 0.84])       # heatmap
    cax = fig.add_axes([0.90, 0.30, 0.022, 0.42])     # colorbar, far right
    im = ax.imshow(M.values, aspect="auto", cmap=DIVERGING,
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(STATE_ORDER)))
    ax.set_xticklabels(STATE_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontstyle="italic", fontsize=6)
    ax.set_xlabel("Microglial state")

    # module separators (white rules between modules)
    cl = meta.loc[order, "string_cluster_id"].to_numpy()
    bounds = np.where(np.diff(cl) != 0)[0] + 0.5
    for b in bounds:
        ax.axhline(b, color="white", lw=1.4)
    # module bracket labels in the whitespace right of the heatmap
    seg_start = 0
    for k in range(len(order) + 1):
        if k == len(order) or cl[k] != cl[seg_start]:
            mid = (seg_start + k - 1) / 2
            ax.annotate(MOD_SHORT[cl[seg_start]],
                        xy=(len(STATE_ORDER) - 0.5, mid),
                        xytext=(len(STATE_ORDER) + 0.4, mid),
                        va="center", ha="left", fontsize=6.2, color="0.25",
                        annotation_clip=False)
            seg_start = k
    # state-group separators (vertical)
    for g, s, e in _state_group_spans(STATE_ORDER):
        if s > 0:
            ax.axvline(s - 0.5, color="0.35", lw=0.8)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("mean SCT Pearson residual", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax.set_title("Adult microglial expression of 41 lysosomal candidates,\n"
                 "grouped by STRING module", loc="left", fontsize=8)
    _save(fig, "H1_gene_state_heatmap")
    return fig


# ==========================================================================
# H2 — module x state activity heatmap
# ==========================================================================
def fig_h2(act, con):
    cl_ids = [int(r.split("_")[1]) for r in act.index]
    labels = [f"{MOD_SHORT[c]}" for c in cl_ids]
    A = act.copy()
    A.index = labels
    vmax = np.abs(A.values).max()

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    im = ax.imshow(A.values, aspect="auto", cmap=DIVERGING, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(STATE_ORDER)))
    ax.set_xticklabels(STATE_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Microglial state")
    # print values (72 cells < 200)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            v = A.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                    color="white" if abs(v) > 0.6 * vmax else "0.15")
    for g, s, e in _state_group_spans(STATE_ORDER):
        if s > 0:
            ax.axvline(s - 0.5, color="0.35", lw=0.8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, aspect=25)
    cbar.set_label("module activity\n(score_genes, mean residual)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax.set_title("Module activity across adult microglial states\n"
                 "(score_genes, mean SCT residual)", loc="left", fontsize=8)
    fig.subplots_adjust(left=0.24, right=0.90, top=0.90, bottom=0.22)
    _save(fig, "H2_module_activity_heatmap")
    return fig


# ==========================================================================
# H3 — coherence + composition
# ==========================================================================
def fig_h3(con):
    c = con.sort_values("string_cluster_id").copy()
    c["short"] = c["string_cluster_id"].map(MOD_SHORT)
    y = np.arange(len(c))[::-1]  # top-to-bottom by cluster id

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.6, 4.2),
                                   gridspec_kw=dict(wspace=0.55))
    # -- (a) coherence: all-cells vs top-state --
    all_c = c["coherence_all_cells"].to_numpy(dtype=float)
    top_c = c["coherence_in_top_state"].to_numpy(dtype=float)
    axL.hlines(y, np.minimum(np.nan_to_num(all_c), np.nan_to_num(top_c)),
               np.maximum(np.nan_to_num(all_c), np.nan_to_num(top_c)),
               color="0.8", lw=1.5, zorder=1)
    axL.scatter(all_c, y, s=34, color="0.45", label="all cells", zorder=3)
    axL.scatter(top_c, y, s=34, color="#C44E52",
                label="within most-active state", zorder=3)
    # single-gene module: n.d.
    for yi, (a, t) in zip(y, zip(all_c, top_c)):
        if np.isnan(a) and np.isnan(t):
            axL.text(0.002, yi, "n.d. (1 gene)", va="center", fontsize=6, color="0.5")
    axL.set_yticks(y)
    axL.set_yticklabels(c["short"])
    axL.set_xlabel("mean pairwise Spearman r\n(within-module co-expression)")
    axL.axvline(0, color="0.6", lw=0.8, ls=":")
    axL.set_title("Within-module coherence is modest", loc="left", fontsize=8)
    axL.legend(frameon=False, fontsize=6, loc="lower right")
    set_frame(axL, "open")

    # -- (b) composition: stable/dynamic stacked counts --
    stable = c["Stable"].to_numpy()
    dynamic = c["Dynamic"].to_numpy()
    axR.barh(y, stable, color=STAB_COLORS["Stable"], label="Stable")
    axR.barh(y, dynamic, left=stable, color=STAB_COLORS["Dynamic"], label="Dynamic")
    for yi, s, d in zip(y, stable, dynamic):
        if s > 0:
            axR.text(s / 2, yi, str(int(s)), ha="center", va="center",
                     color="white", fontsize=6)
        if d > 0:
            axR.text(s + d / 2, yi, str(int(d)), ha="center", va="center",
                     color="white", fontsize=6)
    axR.set_yticks(y)
    axR.set_yticklabels(c["short"])
    axR.set_xlabel("number of candidate genes")
    axR.set_title("Developmental composition per module", loc="left", fontsize=8)
    axR.legend(frameon=False, fontsize=6, loc="lower right", ncol=1)
    set_frame(axR, "open")

    panel_letter(axL, "a")
    panel_letter(axR, "b")
    fig.subplots_adjust(left=0.20, right=0.97, top=0.90, bottom=0.18)
    _save(fig, "H3_module_conservation")
    return fig


# ==========================================================================
# H4 — module-centric synthesis dashboard (Phase 3 reference)
# ==========================================================================
def fig_h4(con, genes):
    c = con.sort_values("peak_activity", ascending=True).copy()
    c["short"] = c["string_cluster_id"].map(MOD_SHORT)
    # representative gene = highest across-state range gene in module
    summ = pd.read_csv(os.path.join(TABLES_DIR, "phase2_gene_expression_summary.csv"))
    rep = {}
    for cid in c["string_cluster_id"]:
        sub = summ[summ.string_cluster_id == cid].sort_values(
            "across_state_residual_range", ascending=False)
        rep[cid] = sub.iloc[0]["symbol"] if len(sub) else "-"
    y = np.arange(len(c))

    fig, (ax0, ax1, ax2) = plt.subplots(
        1, 3, figsize=(10.4, 4.6), gridspec_kw=dict(width_ratios=[2.0, 1.3, 1.3], wspace=0.08))

    # -- track 0: peak activity lollipop, colored by preferred-state group --
    groups = [STATE_GROUP[s] for s in c["most_active_state"]]
    colors = [GROUP_COLORS[g] for g in groups]
    ax0.hlines(y, 0, c["peak_activity"], color="0.8", lw=1.5, zorder=1)
    ax0.scatter(c["peak_activity"], y, s=70, color=colors, zorder=3,
                edgecolor="0.3", linewidth=0.5)
    ax0.set_yticks(y)
    ax0.set_yticklabels(c["short"])
    ax0.set_xlabel("peak module activity\n(at most-active state)")
    ax0.axvline(0, color="0.6", lw=0.8, ls=":")
    # annotate preferred state + representative gene at the dot
    for yi, (_, row) in zip(y, c.iterrows()):
        ax0.text(row["peak_activity"] + 0.012, yi,
                 f"{row['most_active_state']}  ·  {rep[row['string_cluster_id']]}",
                 va="center", ha="left", fontsize=6.2, color="0.2")
    ax0.set_xlim(-0.02, 0.72)
    ax0.set_title("Each candidate module peaks in one microglial state",
                  loc="left", fontsize=8)
    set_frame(ax0, "open")

    # -- track 1: pct stable horizontal bar --
    pct_stable = 100 * c["Stable"] / c["n_genes"]
    ax1.barh(y, pct_stable, color=STAB_COLORS["Stable"])
    ax1.barh(y, 100 - pct_stable, left=pct_stable, color=STAB_COLORS["Dynamic"])
    ax1.set_yticks(y)
    ax1.set_yticklabels([])
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("developmental\ncomposition (%)")
    for yi, ps, n in zip(y, pct_stable, c["n_genes"]):
        ax1.text(50, yi, f"n={int(n)}", ha="center", va="center",
                 fontsize=5.8, color="white")
    ax1.set_title("Stable vs Dynamic", loc="left", fontsize=8)
    set_frame(ax1, "open")

    # -- track 2: coherence bar (all cells) --
    coh = c["coherence_all_cells"].to_numpy(dtype=float)
    ax2.barh(y, np.nan_to_num(coh), color="0.55")
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("within-module\ncoherence (Spearman r)")
    for yi, v in zip(y, coh):
        if np.isnan(v):
            ax2.text(0.002, yi, "n.d.", va="center", fontsize=5.8, color="0.5")
    ax2.set_title("Co-expression", loc="left", fontsize=8)
    set_frame(ax2, "open")

    # shared legends
    grp_handles = [Line2D([0], [0], marker="o", ls="", mfc=GROUP_COLORS[g],
                   mec="0.3", ms=7, label=g) for g in
                   ["Homeostatic", "DAM", "DIM", "Macrophage"]]
    stab_handles = [Rectangle((0, 0), 1, 1, fc=STAB_COLORS[k]) for k in
                    ("Stable", "Dynamic")]
    leg1 = ax0.legend(handles=grp_handles, title="preferred-state group",
                      frameon=False, fontsize=6, title_fontsize=6,
                      loc="lower right", ncol=2)
    ax0.add_artist(leg1)
    ax1.legend(stab_handles, ["Stable", "Dynamic"], frameon=False, fontsize=6,
               loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

    fig.suptitle("Module-centric reference framework for Phase 3 (NPC pseudotime)",
                 x=0.02, ha="left", fontsize=9, fontweight="bold")
    fig.subplots_adjust(left=0.19, right=0.985, top=0.86, bottom=0.20)
    _save(fig, "H4_module_synthesis")
    return fig


def main():
    ensure_dirs()
    apply_figure_style(frame="open", font="Arial", sizes=(8, 7, 6))
    mean_mat, act, con, genes = _load()
    fig_h1(mean_mat, genes)
    fig_h2(act, con)
    fig_h3(con)
    fig_h4(con, genes)
    LOG.info("All figures H1-H4 written to %s", FIGURES_DIR)


if __name__ == "__main__":
    main()
