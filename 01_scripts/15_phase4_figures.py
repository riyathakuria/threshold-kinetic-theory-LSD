#!/usr/bin/env python
"""Phase 4 publication figures P4-1..P4-5.
Regenerates all five orthogonal-evidence-integration figures from the
consolidated matrices and ranked outputs. Run from the workspace root.
"""
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from _phase4_config import TABLES_DIR, FIGURES_DIR, ensure_dirs
from _figstyle import apply_figure_style

apply_figure_style(frame="open", font="Arial")
INT_C, EXT_C = "#55A868", "#8172B3"
FULL_NAMES = {1:"PKA/glucagon signalling",4:"EGFR signaling pathway",
              5:"Rho signal transduction",6:"Eukaryotic Translation Elongation",
              7:"Apoptosis",8:"Blood microparticle (A2M)",
              9:"Response to unfolded protein",10:"Sphingolipid metabolism"}
CALL_C = {"Strong":"#C44E52","Moderate":"#DD8452","Ambiguous":"#CCB974","Stable":"#B0B0B0"}


def _save(fig, name):
    fig.savefig(os.path.join(FIGURES_DIR, name + ".png"), dpi=300)
    fig.savefig(os.path.join(FIGURES_DIR, name + ".pdf"))
    plt.close(fig)


def _mm(s):
    s = s.astype(float); lo, hi = s.min(), s.max()
    return (s - lo) / (hi - lo) if hi > lo else s * 0.0


def _score_modules():
    """Reuse the exact approved scoring from 14_integrate_prioritize.py
    (single source of truth for min-max normalized layers)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_integrate14", os.path.join(_HERE, "14_integrate_prioritize.py"))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    cons = pd.read_csv(os.path.join(TABLES_DIR, "consolidated_module_evidence_matrix.csv"))
    return mod.score_modules(cons)


def load():
    M = _score_modules()                       # has n_A1..n_E5, blocks, ModuleScore
    grank = pd.read_csv(os.path.join(TABLES_DIR, "gene_prioritization_ranked.csv"))
    sens = pd.read_csv(os.path.join(TABLES_DIR, "module_weighting_sensitivity.csv")).set_index("module_id")
    fe = pd.read_csv(os.path.join(TABLES_DIR, "final_module_evidence_table.csv"))
    fe["mid"] = fe["module"].str.replace("M", "").astype(int)
    call_map = dict(zip(fe["mid"], fe["evidence_call"]))
    return M, grank, sens, call_map


def fig1_heatmap(M):
    order = M.sort_values("ModuleScore", ascending=False).index.tolist()
    cols = {"A1 coherence":"n_A1","A2 kinetics×conf":"n_A2","A3 dynamics":"n_A3",
            "E2 ASM network":"n_E2","E3 regulatory":"n_E3","E5 disease":"n_E5"}
    L = pd.DataFrame({k: M[v] for k, v in cols.items()})
    L["E4 brain (gate)"] = _mm(M["E4_brain"])
    L = L.loc[order]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    data = L.values.astype(float)
    im = ax.imshow(np.ma.masked_invalid(data), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(L.shape[1])); ax.set_xticklabels(L.columns, rotation=35, ha="right", fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"M{i}  {FULL_NAMES[i]}" for i in order], fontsize=7.5)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, "n.d.", ha="center", va="center", fontsize=6, color="#888888")
            else:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.3,
                        color="white" if v < 0.55 else "black")
    ax.axvline(2.5, color="white", lw=2.2); ax.axvline(5.5, color="white", lw=2.2)
    ax.text(1.0, -0.9, "INTERNAL (Phase 1–3)", ha="center", fontsize=7, color=INT_C, fontweight="bold")
    ax.text(4.0, -0.9, "EXTERNAL (ASM)", ha="center", fontsize=7, color=EXT_C, fontweight="bold")
    ax.text(6.0, -0.9, "gate", ha="center", fontsize=7, color="#666666", fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.02)
    cb.set_label("min–max normalized evidence", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    ax.set_title("Orthogonal evidence integration across 8 conserved lysosomal modules",
                 fontsize=9.5, loc="left", pad=18)
    ax.set_ylabel("Module (ranked by combined score, top→bottom)", fontsize=8)
    fig.subplots_adjust(left=0.30, right=0.98, top=0.86, bottom=0.20)
    _save(fig, "P4-1_evidence_heatmap")


def fig2_module(M):
    mo = M.sort_values("ModuleScore", ascending=False).reset_index()[::-1].reset_index(drop=True)
    y = np.arange(len(mo))
    ic = 0.60 * mo["InternalBlock"]; ec = 0.40 * mo["ExternalBlock"]
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ax.barh(y, ic, color=INT_C, label="Internal block (×0.60)", edgecolor="white", lw=0.5)
    ax.barh(y, ec, left=ic, color=EXT_C, label="External block (×0.40)", edgecolor="white", lw=0.5)
    ax.set_yticks(y); ax.set_yticklabels([f"M{r.module_id}  {FULL_NAMES[r.module_id]}" for r in mo.itertuples()], fontsize=7.5)
    ax.set_xlabel("Combined module priority score", fontsize=8)
    for i, tot in enumerate(mo["ModuleScore"]):
        ax.text(tot + 0.008, i, f"{tot:.2f}", va="center", fontsize=6.8)
    m8 = mo.index[mo.module_id == 8][0]
    ax.text(mo.loc[m8, "ModuleScore"] + 0.07, m8, "* noise-driven\n(Ambiguous, NPC non-remodeler)",
            va="center", fontsize=6, color="#B00000", style="italic")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_title("Module prioritization (60% internal / 40% external)", fontsize=9.5, loc="left")
    ax.set_xlim(0, 0.85)
    fig.subplots_adjust(left=0.33, right=0.97, top=0.90, bottom=0.13)
    _save(fig, "P4-2_module_ranking")


def fig3_gene(grank):
    top = grank.head(15)[::-1].reset_index(drop=True)
    y = np.arange(len(top))
    cp = 0.45 * top["G_proximity"]; cd = 0.40 * top["G_disease"]; cr = 0.15 * top["G_regulatory"]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.barh(y, cp, color="#4C72B0", label="ASM proximity (×0.45)", edgecolor="white", lw=0.4)
    ax.barh(y, cd, left=cp, color="#C44E52", label="Disease assoc. (×0.40)", edgecolor="white", lw=0.4)
    ax.barh(y, cr, left=cp + cd, color="#CCB974", label="Regulatory (×0.15)", edgecolor="white", lw=0.4)
    ax.set_yticks(y); ax.set_yticklabels([f"{r.gene}  (M{r.module_id})" for r in top.itertuples()],
                                          fontsize=7.5, fontstyle="italic")
    ax.set_xlabel("Gene priority score", fontsize=8)
    for i, tot in enumerate(top["GeneScore"]):
        ax.text(tot + 0.006, i, f"{tot:.2f}", va="center", fontsize=6.5)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_title("Gene prioritization — top 15 (0.45 proximity / 0.40 disease / 0.15 regulatory)",
                 fontsize=9, loc="left")
    ax.set_xlim(0, 0.95)
    fig.subplots_adjust(left=0.22, right=0.97, top=0.92, bottom=0.11)
    _save(fig, "P4-3_gene_ranking")


def fig4_sensitivity(sens):
    schemes = ["50/50", "60/40", "70/30"]; x = [0, 1, 2]
    palette = {6:"#DD8452", 10:"#4C72B0"}
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for mid in sens.index:
        ranks = [sens.loc[mid, f"rank_{s}"] for s in schemes]
        chg = bool(sens.loc[mid, "tier_changes"])
        col = palette.get(mid, "#B0B0B0"); lw = 2.2 if chg else 1.0
        ax.plot(x, ranks, "-o", color=col, lw=lw, ms=5 if chg else 3.5,
                zorder=3 if chg else 1, alpha=1.0 if chg else 0.6)
        ax.text(-0.06, ranks[0], f"M{mid}", ha="right", va="center", fontsize=7,
                color=col, fontweight="bold" if chg else "normal")
        ax.text(2.06, ranks[2], f"M{mid} {FULL_NAMES[mid][:18]}", ha="left", va="center",
                fontsize=7, color=col, fontweight="bold" if chg else "normal")
    ax.set_xticks(x); ax.set_xticklabels([f"{s}\ninternal:external" for s in schemes], fontsize=7.5)
    ax.set_yticks(range(1, 9)); ax.invert_yaxis(); ax.set_ylabel("Priority rank", fontsize=8)
    ax.set_xlim(-0.75, 3.3)
    for yb in (3.5, 5.5):
        ax.axhline(yb, color="#DDDDDD", lw=0.8, ls="--", zorder=0)
    ax.set_title("Weighting sensitivity: only M6 and M10 cross tiers (Spearman ρ 0.81–0.95)",
                 fontsize=8.5, loc="left")
    fig.subplots_adjust(left=0.11, right=0.99, top=0.92, bottom=0.13)
    _save(fig, "P4-4_sensitivity_slopegraph")


def fig5_synthesis(M, call_map):
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    xi = M["InternalBlock"]; ye = M["ExternalBlock"]
    sizes = (M["ModuleScore"] / M["ModuleScore"].max() * 520) + 60
    def ccol(x):
        for k, v in CALL_C.items():
            if str(x).startswith(k): return v
        return "#B0B0B0"
    cols = [ccol(call_map[m]) for m in M.index]
    ax.scatter(xi, ye, s=sizes, c=cols, alpha=0.78, edgecolors="white", linewidths=1.0, zorder=3)
    for m in M.index:
        dy = 16 if m == 10 else (15 if m == 8 else -14)
        ax.annotate(f"M{m} {FULL_NAMES[m][:20]}", (xi[m], ye[m]), fontsize=6.8, ha="center",
                    va="center", xytext=(0, dy), textcoords="offset points")
    ax.axvline(xi.median(), color="#E0E0E0", lw=0.8, ls="--", zorder=0)
    ax.axhline(ye.median(), color="#E0E0E0", lw=0.8, ls="--", zorder=0)
    ax.text(0.97, 0.03, "discovery-driven\n(NPC dynamics)", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=6.5, color="#999999", style="italic")
    ax.text(0.03, 0.97, "ASM-mechanistic\n(sphingolipid axis)", transform=ax.transAxes, ha="left",
            va="top", fontsize=6.5, color="#999999", style="italic")
    ax.set_xlabel("Internal evidence block (Phase 1–3 discovery + NPC kinetics)", fontsize=8)
    ax.set_ylabel("External evidence block (ASM network + disease + regulatory)", fontsize=8)
    ax.set_title("Internal discovery and external ASM-centrality are orthogonal axes", fontsize=8.8, loc="left")
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=v, markersize=8, label=k)
           for k, v in CALL_C.items()]
    ax.legend(handles=leg, fontsize=6.8, frameon=False, title="NPC evidence call (Phase 3)",
              title_fontsize=7, loc="upper right")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.05)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.10)
    _save(fig, "P4-5_synthesis_orthogonality")


def main():
    ensure_dirs()
    M, grank, sens, call_map = load()
    fig1_heatmap(M); fig2_module(M); fig3_gene(grank)
    fig4_sensitivity(sens); fig5_synthesis(M, call_map)
    print("P4 figures regenerated: P4-1..P4-5")


if __name__ == "__main__":
    main()
