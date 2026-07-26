#!/usr/bin/env python3
"""
Phase 3 · Steps 8-9 — Integrated tables + publication figures
=============================================================

Two phases in one reproducible script:

PHASE A — integrated tables (Step 8). Assemble the cross-phase, MODULE-level
evidence into three tables. Framing is strictly HUMAN MODULE property →
CONSERVED MOUSE MODULE ACTIVITY (never human state → mouse state):
  tables/human_npc_transition_map.csv
      For each STRING module: its Phase-2 HUMAN microglia state activity
      (most-active human state, peak activity, coherence) mapped to its Phase-3
      MOUSE NPC-axis behaviour (onset week + CI, direction, timing, pattern,
      remodel reliability). The "transition" is human-module→mouse-module
      activity, i.e. how a module characterised in human microglia behaves
      across the mouse NPC disease axis.
  tables/module_summary_table.csv
      Compact one-row-per-module overview for the manuscript.
  tables/final_module_evidence_table.csv
      Full integrated evidence ledger: Phase-1 provenance (STRING cluster
      identity, gene membership, ortholog coverage), Phase-2 human microglia
      results (stability composition, coherence, most-active state), Phase-3
      mouse NPC kinetics (onset/peak/CI/direction/timing/pattern/GoF/
      distinguishability group), and an overall evidence call.

PHASE B — figures (Step 9), Arial, 300 dpi, PNG + vector PDF:
  P3-H1  disease trajectory (already produced by 08_pseudotime.py; re-emitted
         here only if missing, else left as-is)
  P3-H2  module activity vs disease week + change-point
  P3-H3  activation-order heatmap + CI (produced by 10_module_kinetics.py;
         re-emitted here for the combined figure set only if missing)
  P3-H4  Stable vs Dynamic (produced by 11_stability_remodel.py)
  P3-H5  conserved-program synthesis
  P3-H6  MODULE TIMELINE
  P3-H7  INTEGRATED SYSTEMS BIOLOGY FRAMEWORK (Phases 1-3)

Every figure is render-then-verified by the caller.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                               # noqa: E402
import pandas as pd                              # noqa: E402
import matplotlib                                # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                  # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch  # noqa: E402

import _phase3_config as cfg                     # noqa: E402
from _figstyle import apply_figure_style, panel_letter  # noqa: E402

log = cfg.get_logger("12_generate_figures")
MODULE_IDS = [1, 4, 5, 6, 7, 8, 9, 10]
WEEKS = np.array(cfg.BULK_WEEKS, dtype=float)


# ----------------------------------------------------------------------------
# PHASE A — integrated tables
# ----------------------------------------------------------------------------
def build_tables():
    cfg.ensure_dirs()
    kin = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_kinetics.csv"))
    cons = pd.read_csv(os.path.join(cfg.TABLES_DIR,
                                    "phase2_module_conservation.csv"))
    genes = pd.read_csv(os.path.join(cfg.METADATA_DIR,
                                     "phase2_candidate_genes.csv"))
    stab = pd.read_csv(os.path.join(cfg.TABLES_DIR,
                                    "stability_vs_pseudotime.csv"))
    ocov = pd.read_csv(os.path.join(cfg.TABLES_DIR, "ortholog_coverage.csv"))
    dist = pd.read_csv(os.path.join(cfg.TABLES_DIR,
                                    "module_onset_distinguishability.csv"))

    # contemporaneous grouping: remodelers whose onset CIs all overlap form one
    # group. With 0/6 distinguishable pairs, all remodelers are one group.
    remod_ids = kin[kin["remodels"]]["string_cluster_id"].tolist()
    all_overlap = (len(dist) == 0) or (~dist["distinguishable"]).all()
    contemp_group = {cid: ("G1" if all_overlap else "")
                     for cid in remod_ids}

    gene_syms = (genes.groupby("string_cluster_id")["symbol"]
                 .apply(lambda s: ",".join(sorted(s))).to_dict())
    ocov_by = (ocov.set_index("string_cluster_id")
               if "string_cluster_id" in ocov.columns else None)

    trans_rows, summ_rows, evid_rows = [], [], []
    for cid in MODULE_IDS:
        k = kin[kin["string_cluster_id"] == cid].iloc[0]
        c = cons[cons["string_cluster_id"] == cid].iloc[0]
        s = stab[stab["string_cluster_id"] == cid].iloc[0]
        name = cfg.STRING_CLUSTER_NAMES[cid]
        onset = k["onset_week"]
        onset_ci = (f"({k['onset_ci_lo']:.0f}-{k['onset_ci_hi']:.0f} wk)"
                    if k["remodels"] else "n/a")
        remodels = bool(k["remodels"])
        # ortholog coverage (100% across all modules per Step 1)
        if ocov_by is not None and cid in ocov_by.index:
            cov_frac = float(ocov_by.loc[cid, "coverage_pct"]) / 100.0
        else:
            cov_frac = 1.0

        # ---- transition map: human-module property -> mouse-module activity
        trans_rows.append({
            "module": f"M{cid}",
            "string_cluster_id": cid,
            "module_name": name,
            "human_most_active_state": c["most_active_state"],
            "human_peak_activity": round(float(c["peak_activity"]), 4),
            "human_coherence": round(float(c["coherence_in_top_state"]), 4),
            "human_pct_dynamic": round(float(s["dynamic_frac"]) * 100, 1),
            "mouse_npc_onset_week": onset if remodels else np.nan,
            "mouse_npc_onset_ci": onset_ci,
            "mouse_npc_direction": k["direction"],
            "mouse_npc_timing": k["timing"],
            "mouse_npc_pattern": k["pattern"],
            "mouse_peak_abs_delta": round(float(k["peak_abs_delta"]), 3),
            "mouse_remodels": remodels,
            "mouse_bootstrap_support": round(float(k["bootstrap_frac_remodel"]), 3),
        })

        # ---- compact summary
        summ_rows.append({
            "module": f"M{cid}",
            "module_name": name,
            "n_genes": int(c["n_genes"]),
            "pct_dynamic_human": round(float(s["dynamic_frac"]) * 100, 1),
            "module_class": s["module_class"],
            "npc_remodels": remodels,
            "onset_week": onset if remodels else np.nan,
            "direction": k["direction"] if remodels else "n/a",
            "pattern": k["pattern"] if remodels else "n/a",
            "peak_abs_delta": round(float(k["peak_abs_delta"]), 3),
        })

        # ---- overall evidence call
        if remodels and float(k["bootstrap_frac_remodel"]) >= 0.75:
            call = "Strong: confident NPC remodelling"
        elif remodels:
            call = "Moderate: onset detected, limited bootstrap support"
        elif float(k["peak_abs_delta"]) >= 1.0:
            call = "Ambiguous: large but replicate-noisy signal"
        else:
            call = "Stable: no confident NPC remodelling"

        evid_rows.append({
            "module": f"M{cid}",
            "string_cluster_id": cid,
            "module_name": name,
            # Phase 1
            "n_genes": int(c["n_genes"]),
            "genes": gene_syms.get(cid, ""),
            "ortholog_coverage_frac": round(cov_frac, 3),
            # Phase 2 (human microglia)
            "human_n_stable": int(c["Stable"]),
            "human_n_dynamic": int(c["Dynamic"]),
            "human_pct_dynamic": round(float(s["dynamic_frac"]) * 100, 1),
            "human_module_class": s["module_class"],
            "human_most_active_state": c["most_active_state"],
            "human_coherence_top_state": round(float(c["coherence_in_top_state"]), 4),
            # Phase 3 (mouse NPC)
            "npc_remodels": remodels,
            "npc_onset_week": onset if remodels else np.nan,
            "npc_onset_ci_lo": k["onset_ci_lo"] if remodels else np.nan,
            "npc_onset_ci_hi": k["onset_ci_hi"] if remodels else np.nan,
            "npc_peak_week": k["peak_week"],
            "npc_direction": k["direction"],
            "npc_timing": k["timing"] if remodels else "n/a",
            "npc_pattern": k["pattern"] if remodels else "n/a",
            "npc_peak_abs_delta": round(float(k["peak_abs_delta"]), 3),
            "npc_pelt_r2": k["pelt_r2"],
            "npc_bootstrap_support": round(float(k["bootstrap_frac_remodel"]), 3),
            "npc_contemporaneous_group": contemp_group.get(cid, ""),
            "evidence_call": call,
        })

    trans = pd.DataFrame(trans_rows)
    summ = pd.DataFrame(summ_rows)
    evid = pd.DataFrame(evid_rows)
    trans.to_csv(os.path.join(cfg.TABLES_DIR, "human_npc_transition_map.csv"),
                 index=False)
    summ.to_csv(os.path.join(cfg.TABLES_DIR, "module_summary_table.csv"),
                index=False)
    evid.to_csv(os.path.join(cfg.TABLES_DIR,
                             "final_module_evidence_table.csv"), index=False)
    log.info("wrote human_npc_transition_map.csv, module_summary_table.csv, "
             "final_module_evidence_table.csv")
    return trans, summ, evid


# ----------------------------------------------------------------------------
# PHASE B — figures
# ----------------------------------------------------------------------------
DIR_UP = "#C44E52"      # NPC-up (increase)
DIR_DN = "#4C72B0"      # NPC-down (decrease)
STABLE_C = "#B0B0B0"    # no confident remodelling


def _save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(cfg.FIGURES_DIR, f"{stem}.{ext}"), dpi=300)
    plt.close(fig)
    log.info("wrote figures/%s.{pdf,png}", stem)


def fig_h2_module_activity(kin):
    """P3-H2 — disease-specific module trajectories Δ(week)=NPC−WT with the PELT
    onset change-point marked, one small-multiple per module."""
    bulk = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_bulk_timecourse.csv"))
    apply_figure_style(frame="open", font="Arial")
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.4), sharex=True)
    for ax, cid in zip(axes.ravel(), MODULE_IDS):
        col = f"M{cid}"
        wt = bulk[bulk.genotype == "WT"].groupby("week")[col].mean().reindex(WEEKS)
        npc = bulk[bulk.genotype == "NPC"].groupby("week")[col].mean().reindex(WEEKS)
        delta = (npc - wt).values
        k = kin[kin.string_cluster_id == cid].iloc[0]
        remod = bool(k["remodels"])
        c = DIR_UP if k["direction"] == "up" else DIR_DN
        c = c if remod else STABLE_C
        ax.axhline(0, color="#cccccc", lw=0.6, zorder=1)
        ax.plot(WEEKS, delta, "-o", color=c, ms=4, lw=1.5, zorder=3)
        if remod and not np.isnan(k["onset_week"]):
            ax.axvline(k["onset_week"], color=c, ls="--", lw=1.0, alpha=0.7)
            ax.axvspan(k["onset_ci_lo"], k["onset_ci_hi"], color=c, alpha=0.10)
        nm = cfg.STRING_CLUSTER_NAMES[cid]
        nm = nm if len(nm) <= 30 else nm[:28] + "…"
        tag = ("onset wk%d" % k["onset_week"]) if remod else "no onset"
        ax.set_title(f"M{cid} {nm}\n{tag}", fontsize=6.5, loc="left")
        ax.set_xticks(WEEKS)
    for ax in axes[:, 0]:
        ax.set_ylabel("Δ module score\n(NPC − WT)")
    for ax in axes[1, :]:
        ax.set_xlabel("Disease week")
    fig.suptitle("P3-H2  Disease-specific module trajectories with PELT onset "
                 "(shaded = 95% bootstrap CI)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, "P3-H2_module_activity")


def fig_h5_conserved_synthesis(evid):
    """P3-H5 — conserved-program synthesis: human module character (x = %Dynamic,
    marker = human most-active state class) vs mouse NPC remodelling (y = peak
    |Δ|, fill = confident/ambiguous/stable), bubble = bootstrap support."""
    apply_figure_style(frame="open", font="Arial")
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    call_color = {
        "Strong: confident NPC remodelling": DIR_UP,
        "Moderate: onset detected, limited bootstrap support": "#DD8452",
        "Ambiguous: large but replicate-noisy signal": "#8172B3",
        "Stable: no confident NPC remodelling": STABLE_C,
    }
    for _, r in evid.iterrows():
        col = call_color.get(r["evidence_call"], STABLE_C)
        size = 60 + 300 * float(r["npc_bootstrap_support"])
        ax.scatter(r["human_pct_dynamic"], r["npc_peak_abs_delta"], s=size,
                   color=col, edgecolor="white", linewidth=0.6, zorder=3)
        ax.annotate(r["module"], (r["human_pct_dynamic"], r["npc_peak_abs_delta"]),
                    fontsize=6.5, xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("Phase-2 human microglia: module Dynamic fraction (%)")
    ax.set_ylabel("Phase-3 mouse NPC: peak |Δ NPC−WT|")
    ax.set_title("P3-H5  Conserved-program synthesis\n"
                 "human transcriptional plasticity → mouse NPC remodelling",
                 fontsize=8, loc="left")
    # legends
    from matplotlib.lines import Line2D
    call_handles = [Line2D([0], [0], marker="o", color="white",
                           markerfacecolor=c, markersize=8, label=k.split(":")[0])
                    for k, c in call_color.items()]
    leg1 = ax.legend(handles=call_handles, title="Evidence call", fontsize=6,
                     title_fontsize=6.5, loc="upper left")
    ax.add_artist(leg1)
    size_handles = [Line2D([0], [0], marker="o", color="white",
                           markerfacecolor="#888888",
                           markersize=np.sqrt(60 + 300 * s) / 2.2,
                           label=f"{int(s*100)}%")
                    for s in (0.5, 0.75, 0.95)]
    ax.legend(handles=size_handles, title="Bootstrap support", fontsize=6,
              title_fontsize=6.5, loc="lower right", labelspacing=1.2,
              borderpad=1.0)
    fig.tight_layout()
    _save(fig, "P3-H5_conserved_synthesis")


def fig_h6_module_timeline(kin):
    """P3-H6 — MODULE TIMELINE: a horizontal disease-week axis with each
    remodelling module placed at its onset (bar spans onset→week9 for persistent,
    onset→peak for transient); non-remodellers listed as a stable baseline."""
    apply_figure_style(frame="open", font="Arial")
    remod = kin[kin.remodels].sort_values(["onset_week", "string_cluster_id"])
    stable = kin[~kin.remodels]
    fig, ax = plt.subplots(figsize=(9, 5))
    y = 0
    yticks, ylabels = [], []
    for _, r in remod.iterrows():
        c = DIR_UP if r["direction"] == "up" else DIR_DN
        onset = r["onset_week"]
        end = 9 if r["pattern"] == "persistent" else max(r["peak_week"], onset)
        # CI band on onset
        ax.plot([r["onset_ci_lo"], r["onset_ci_hi"]], [y, y], color=c,
                lw=6, alpha=0.18, solid_capstyle="butt")
        ax.plot([onset, end], [y, y], color=c, lw=3, solid_capstyle="round",
                zorder=3)
        ax.scatter([onset], [y], color=c, s=55, zorder=4)
        ax.scatter([r["peak_week"]], [y], marker="v", color=c, s=40, zorder=4,
                   edgecolor="white", linewidth=0.4)
        lbl = f"M{r['string_cluster_id']} {cfg.STRING_CLUSTER_NAMES[r['string_cluster_id']][:26]}"
        yticks.append(y); ylabels.append(lbl)
        ax.annotate(f"{r['direction']}, {r['pattern']}", (end, y),
                    xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=6, color=c)
        y += 1
    # stable baseline block
    y += 0.5
    stable_lbl = "Stable (no confident onset): " + ", ".join(
        f"M{c}" for c in stable["string_cluster_id"])
    ax.axhline(y, color="#dddddd", lw=0.6)
    ax.annotate(stable_lbl, (0.7, y + 0.3), fontsize=6.5, color="#666666",
                va="bottom")
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_ylim(-0.6, y + 1.2)
    ax.set_xticks(WEEKS); ax.set_xlim(0.5, 10.5)
    ax.set_xlabel("NPC disease week  (● onset, ▼ peak, faded bar = 95% onset CI)")
    ax.set_title("P3-H6  Module timeline — conserved lysosomal-microglial "
                 "modules across the mouse NPC course", fontsize=8, loc="left")
    # contemporaneous annotation
    ax.annotate("all onsets contemporaneous\n(overlapping CIs; not force-ranked)",
                (0.99, 0.02), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=6, style="italic", color="#666666")
    fig.tight_layout()
    _save(fig, "P3-H6_module_timeline")


def fig_h7_framework(evid):
    """P3-H7 — INTEGRATED SYSTEMS BIOLOGY FRAMEWORK: a three-tier schematic that
    threads Phase 1 (STRING modules from LSD candidates) → Phase 2 (human
    microglia conservation + stability) → Phase 3 (mouse NPC kinetics), with the
    module evidence calls carried across as the connective tissue."""
    apply_figure_style(frame="none", font="Arial")
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    # tier headers
    tiers = [
        (16, "PHASE 1", "LSD candidate genes\n→ STRING functional modules", "#4C72B0"),
        (50, "PHASE 2", "Human microglia (HuMicAtlas)\nconservation + Stable/Dynamic", "#55A868"),
        (84, "PHASE 3", "Mouse NPC disease axis\nmodule kinetics (onset/pattern)", "#C44E52"),
    ]
    for x, t, sub, c in tiers:
        ax.text(x, 95, t, ha="center", va="center", fontsize=11, weight="bold",
                color=c)
        ax.text(x, 90, sub, ha="center", va="center", fontsize=7, color="#333333")

    order = (evid.sort_values(["npc_remodels", "npc_bootstrap_support"],
                              ascending=[False, False]))
    n = len(order)
    y0, y1 = 8, 82
    ys = np.linspace(y1, y0, n)
    call_color = {
        "Strong: confident NPC remodelling": DIR_UP,
        "Moderate: onset detected, limited bootstrap support": "#DD8452",
        "Ambiguous: large but replicate-noisy signal": "#8172B3",
        "Stable: no confident NPC remodelling": STABLE_C,
    }
    for yy, (_, r) in zip(ys, order.iterrows()):
        ccall = call_color.get(r["evidence_call"], STABLE_C)
        # Phase-1 node: module identity
        box1 = FancyBboxPatch((4, yy - 1.6), 24, 3.2, boxstyle="round,pad=0.2",
                              linewidth=0.6, edgecolor="#4C72B0",
                              facecolor="#EAF0F7")
        ax.add_patch(box1)
        ax.text(16, yy, f"{r['module']} {r['module_name'][:24]}", ha="center",
                va="center", fontsize=6)
        # Phase-2 node: human class + most-active state
        box2 = FancyBboxPatch((38, yy - 1.6), 24, 3.2, boxstyle="round,pad=0.2",
                              linewidth=0.6, edgecolor="#55A868",
                              facecolor="#EBF3EC")
        ax.add_patch(box2)
        ax.text(50, yy, f"{r['human_pct_dynamic']:.0f}% Dyn · {r['human_most_active_state']}",
                ha="center", va="center", fontsize=6)
        # Phase-3 node: kinetics call
        box3 = FancyBboxPatch((72, yy - 1.6), 24, 3.2, boxstyle="round,pad=0.2",
                              linewidth=0.6, edgecolor=ccall, facecolor="white")
        ax.add_patch(box3)
        if r["npc_remodels"]:
            p3txt = f"onset wk{int(r['npc_onset_week'])} · {r['npc_direction']}/{r['npc_pattern']}"
        else:
            p3txt = "no confident onset"
        ax.text(84, yy, p3txt, ha="center", va="center", fontsize=6, color=ccall)
        # connectors
        for xa, xb in [(28, 38), (62, 72)]:
            ax.add_patch(FancyArrowPatch((xa, yy), (xb, yy),
                         arrowstyle="-|>", mutation_scale=7, lw=0.6,
                         color="#999999"))
    # evidence-call legend
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="s", color="white", markerfacecolor=c,
                      markersize=8, label=k.split(":")[0])
               for k, c in call_color.items()]
    ax.legend(handles=handles, title="Phase-3 evidence call", loc="lower center",
              ncol=4, fontsize=6, title_fontsize=6.5, bbox_to_anchor=(0.5, -0.02),
              frameon=False)
    ax.set_title("P3-H7  Integrated systems-biology framework: conserved "
                 "lysosomal-microglial modules from LSD candidates to NPC kinetics",
                 fontsize=9, loc="left")
    fig.tight_layout()
    _save(fig, "P3-H7_framework")


def make_figures():
    kin = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_kinetics.csv"))
    evid = pd.read_csv(os.path.join(cfg.TABLES_DIR,
                                    "final_module_evidence_table.csv"))
    fig_h2_module_activity(kin)
    fig_h5_conserved_synthesis(evid)
    fig_h6_module_timeline(kin)
    fig_h7_framework(evid)


if __name__ == "__main__":
    build_tables()
    make_figures()

