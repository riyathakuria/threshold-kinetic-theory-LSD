#!/usr/bin/env python3
"""
Phase 3 · Step 10 — Module kinetics: change-point, CIs & goodness-of-fit
========================================================================

For each Phase-2 STRING module, determine WHEN it begins to remodel along the
PRIMARY disease axis (GSE152158 bulk real-time course, weeks 1/3/6/9) and
classify its temporal pattern. Change-point detection is ruptures PELT
(PRIMARY); pwlf piecewise-linear is a SENSITIVITY check only.

Disease-specific trajectory. To separate NPC disease dynamics from
developmental/age effects shared with WT, the disease signal for each module
is the per-week difference  Δ(week) = mean_NPC(week) − mean_WT(week)  in the
module score. Onset/peak/pattern are defined on Δ.

Resolution caveat (documented, not hidden). The bulk course samples four weeks
(1,3,6,9). Onset is therefore resolved only to the SAMPLED weeks; we do not
claim finer temporal precision than the design supports. Change-point tools are
applied to the replicate-ordered series so the estimate uses the full replicate
structure, and uncertainty is quantified by bootstrapping replicates.

Reported per module:
  onset_week        first week Δ departs from the week-1 baseline (PELT primary)
  onset_ci          bootstrap 95% CI for onset (percentile)
  peak_week         week of maximum |Δ| in the disease direction
  peak_ci           bootstrap 95% CI for peak
  direction         up / down in NPC vs WT
  timing            early (onset ≤3 wk) / late (onset ≥6 wk)
  pattern           transient (declines after peak) / persistent (holds to wk9)
  pelt_r2           goodness-of-fit of the PELT piecewise-constant model to Δ
  pwlf_break/pwlf_r2  sensitivity (piecewise-linear breakpoint on Δ)
  remodels          whether a change is detected at all
Contemporaneity: modules whose onset CIs overlap are flagged as a
CONTEMPORANEOUS group — temporal order is reported as ASSOCIATION, not causation.

Also records agreement with the SUPPORTING pseudotime axis (DPT-binned), which
did NOT pass the disease-relevance gate and is reported for completeness only.

Outputs:
  tables/module_kinetics.csv
  tables/module_onset_distinguishability.csv
  figures/P3-H3_activation_order.pdf/png   (onset heatmap + CIs)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                            # noqa: E402
import pandas as pd                           # noqa: E402
import ruptures as rpt                        # noqa: E402
import pwlf                                    # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402

import _phase3_config as cfg                  # noqa: E402
from _figstyle import apply_figure_style, panel_letter  # noqa: E402

log = cfg.get_logger("10_module_kinetics")
np.random.seed(cfg.RANDOM_SEED)

WEEKS = np.array(cfg.BULK_WEEKS, dtype=float)   # [1,3,6,9]
MODULE_IDS = [1, 4, 5, 6, 7, 8, 9, 10]


def load_bulk():
    df = pd.read_csv(os.path.join(cfg.TABLES_DIR, "module_bulk_timecourse.csv"))
    return df


def disease_dev(df, cid, samples=None):
    """Per-NPC-sample disease-deviation series d_i = score_i − mean_WT(week_i),
    ordered by week. This is the disease-specific signal (developmental/age
    effects shared with WT are subtracted per week). Returns (d, weeks)."""
    d = df if samples is None else df.loc[samples]
    col = f"M{cid}"
    wtmean = {wk: d[(d["genotype"] == "WT") & (d["week"] == wk)][col].mean()
              for wk in WEEKS}
    npc = d[d["genotype"] == "NPC"].sort_values("week")
    dev = npc[col].values - np.array([wtmean[wk] for wk in npc["week"].values])
    return dev.astype(float), npc["week"].values.astype(float)


def week_delta(df, cid, samples=None):
    """Disease-specific per-week trajectory Δ(week)=mean_NPC-mean_WT (= mean of
    the per-sample disease-deviation d_i within each week)."""
    dev, w = disease_dev(df, cid, samples)
    return np.array([dev[w == wk].mean() if np.any(w == wk) else np.nan
                     for wk in WEEKS], dtype=float)


def _within_week_sigma2(dev, w):
    """Noise model: pooled WITHIN-week replicate variance of the disease-
    deviation series. Replicate scatter at a fixed week is measurement/biological
    noise; between-week shifts are the signal to be detected. Estimating sigma^2
    this way (rather than from lag-1 differences of the ordered series) prevents a
    genuine step from inflating its own penalty and suppressing its detection."""
    res = []
    for wk in np.unique(w):
        x = dev[w == wk]
        if len(x) > 1:
            res.append(x - x.mean())
    if not res:
        return 1e-9
    res = np.concatenate(res)
    return max(float(np.var(res, ddof=1)), 1e-9)


def pelt_onset(dev, w):
    """PELT (l2) on the disease-deviation series ordered by week. BIC-style
    penalty pen = log(n) * sigma^2 (Killick 2012), sigma^2 = within-week variance.
    Returns onset week (first detected mean shift) or np.nan if no breakpoint."""
    order = np.argsort(w, kind="stable")
    sig = dev[order].reshape(-1, 1)
    ww = w[order]
    n = len(sig)
    sigma2 = _within_week_sigma2(dev, w)
    pen = np.log(n) * 1 * sigma2
    algo = rpt.Pelt(model=cfg.RUPTURES_MODEL,
                    min_size=cfg.RUPTURES_MIN_SIZE).fit(sig)
    bkps = [b for b in algo.predict(pen=pen) if b < n]
    if not bkps:
        return np.nan, sig.ravel(), ww, []
    return float(ww[bkps[0]]), sig.ravel(), ww, bkps


def pelt_r2(sig, weeks, bkps):
    """R^2 of the PELT piecewise-constant fit."""
    n = len(sig)
    segs = [0] + [b for b in bkps if b < n] + [n]
    pred = np.empty(n)
    for a, b in zip(segs[:-1], segs[1:]):
        pred[a:b] = sig[a:b].mean()
    ss_res = float(((sig - pred) ** 2).sum())
    ss_tot = float(((sig - sig.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def peak_week(delta):
    """Week of maximum |Δ| in the disease direction (sign = overall direction)."""
    direction = np.sign(delta[np.argmax(np.abs(delta))])
    signed = delta * direction
    return float(WEEKS[int(np.argmax(signed))]), direction


def classify(onset, pk, delta, direction):
    if np.isnan(onset):
        return "none", "no_change"
    timing = "early" if onset <= 3 else "late"
    signed = delta * direction
    # persistent if week-9 level is >=70% of peak level; else transient
    peak_val = signed.max()
    end_val = signed[-1]
    pattern = "persistent" if (peak_val > 0 and end_val >= 0.7 * peak_val) else "transient"
    return timing, pattern


def bootstrap_ci(df, cid, n_boot=cfg.BOOTSTRAP_N):
    """Resample replicates within each (genotype, week) cell; recompute onset &
    peak. Return percentile CIs."""
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    idx_by_cell = {k: g.index.values for k, g in df.groupby(["genotype", "week"])}
    onsets, peaks = [], []
    col = f"M{cid}"
    for _ in range(n_boot):
        samp = []
        for k, idxs in idx_by_cell.items():
            samp.extend(rng.choice(idxs, size=len(idxs), replace=True))
        sub = df.loc[samp].reset_index(drop=True)
        dev, w = disease_dev(sub, cid)
        ow, sig, ww, bkps = pelt_onset(dev, w)
        onsets.append(ow)
        d = week_delta(sub, cid)
        if np.all(np.isnan(d)):
            peaks.append(np.nan)
        else:
            pk, _ = peak_week(np.nan_to_num(d))
            peaks.append(pk)
    onsets = np.array(onsets, float)
    peaks = np.array(peaks, float)
    lo, hi = (100 - cfg.BOOTSTRAP_CI) / 2, 100 - (100 - cfg.BOOTSTRAP_CI) / 2
    def ci(a):
        a = a[~np.isnan(a)]
        if len(a) == 0:
            return (np.nan, np.nan)
        return (float(np.percentile(a, lo)), float(np.percentile(a, hi)))
    frac_remodel = float(np.mean(~np.isnan(onsets)))
    return ci(onsets), ci(peaks), frac_remodel


def main():
    cfg.ensure_dirs()
    df = load_bulk()
    rows = []
    for cid in MODULE_IDS:
        col = f"M{cid}"
        dev, w = disease_dev(df, cid)
        onset, sig, ww, bkps = pelt_onset(dev, w)
        delta = week_delta(df, cid)
        pk, direction = peak_week(delta)
        timing, pattern = classify(onset, pk, delta, direction)
        r2 = pelt_r2(sig, ww, bkps) if bkps else np.nan

        # pwlf sensitivity on Δ(week)
        try:
            m = pwlf.PiecewiseLinFit(WEEKS, delta)
            brk = m.fit(2)
            pwlf_break = float(brk[1])
            pwlf_r2v = float(m.r_squared())
        except Exception:
            pwlf_break, pwlf_r2v = np.nan, np.nan

        (o_lo, o_hi), (p_lo, p_hi), frac = bootstrap_ci(df, cid)
        rows.append({
            "string_cluster_id": cid,
            "module_name": cfg.STRING_CLUSTER_NAMES[cid],
            "onset_week": onset,
            "onset_ci_lo": o_lo, "onset_ci_hi": o_hi,
            "peak_week": pk,
            "peak_ci_lo": p_lo, "peak_ci_hi": p_hi,
            "direction": "up" if direction > 0 else "down",
            "timing": timing,
            "pattern": pattern,
            "delta_week9": round(float(delta[-1]), 3),
            "peak_abs_delta": round(float(np.max(np.abs(delta))), 3),
            "pelt_r2": round(r2, 3) if not np.isnan(r2) else np.nan,
            "pwlf_break_week": round(pwlf_break, 2) if not np.isnan(pwlf_break) else np.nan,
            "pwlf_r2": round(pwlf_r2v, 3) if not np.isnan(pwlf_r2v) else np.nan,
            "bootstrap_frac_remodel": round(frac, 3),
            "remodels": bool(not np.isnan(onset)),
        })
        log.info("M%d %-45s onset=%s CI=(%.1f,%.1f) peak=%s dir=%s %s/%s "
                 "R2pelt=%.2f", cid, cfg.STRING_CLUSTER_NAMES[cid],
                 onset, o_lo, o_hi, pk, "up" if direction > 0 else "down",
                 timing, pattern, r2 if not np.isnan(r2) else -1)

    kin = pd.DataFrame(rows)
    kin.to_csv(os.path.join(cfg.TABLES_DIR, "module_kinetics.csv"), index=False)
    log.info("wrote tables/module_kinetics.csv")

    distinguish(kin)
    make_h3_figure(kin)
    return kin


def distinguish(kin):
    """Pairwise onset-CI overlap → contemporaneous grouping."""
    remod = kin[kin["remodels"]].copy()
    pairs = []
    for i, a in remod.iterrows():
        for j, b in remod.iterrows():
            if a["string_cluster_id"] >= b["string_cluster_id"]:
                continue
            overlap = not (a["onset_ci_hi"] < b["onset_ci_lo"] or
                           b["onset_ci_hi"] < a["onset_ci_lo"])
            pairs.append({
                "module_a": f"M{a['string_cluster_id']}",
                "module_b": f"M{b['string_cluster_id']}",
                "onset_a": a["onset_week"], "onset_b": b["onset_week"],
                "ci_a": f"({a['onset_ci_lo']:.1f},{a['onset_ci_hi']:.1f})",
                "ci_b": f"({b['onset_ci_lo']:.1f},{b['onset_ci_hi']:.1f})",
                "distinguishable": not overlap,
            })
    d = pd.DataFrame(pairs)
    d.to_csv(os.path.join(cfg.TABLES_DIR,
                          "module_onset_distinguishability.csv"), index=False)
    n_dist = int(d["distinguishable"].sum()) if len(d) else 0
    log.info("onset distinguishability: %d/%d module pairs statistically "
             "distinguishable (non-overlapping CIs)", n_dist, len(d))


def make_h3_figure(kin):
    apply_figure_style(frame="open", font="Arial")
    remod = kin[kin["remodels"]].sort_values("onset_week")
    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(remod))
    for i, (_, r) in enumerate(remod.iterrows()):
        ax.plot([r["onset_ci_lo"], r["onset_ci_hi"]], [i, i], color="#888888",
                lw=3, alpha=0.5, solid_capstyle="round")
        c = "#C44E52" if r["direction"] == "up" else "#4C72B0"
        ax.scatter(r["onset_week"], i, s=70, color=c, zorder=3,
                   label=r["direction"])
        ax.scatter(r["peak_week"], i, s=40, marker="D", facecolor="none",
                   edgecolor=c, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"M{r['string_cluster_id']} {r['module_name'][:28]}"
                        for _, r in remod.iterrows()], fontsize=6)
    ax.set_xticks(WEEKS)
    ax.set_xlabel("Disease week (onset = filled circle, peak = open diamond; "
                  "grey bar = 95% bootstrap CI)")
    ax.set_title("P3-H3 Module activation order (bulk real-time axis)")
    ax.set_xlim(0.5, 9.5)
    # dedup legend
    h, l = ax.get_legend_handles_labels()
    seen = dict(zip(l, h))
    ax.legend(seen.values(), [f"NPC {k}" for k in seen], loc="lower right",
              title="direction")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(cfg.FIGURES_DIR,
                                 f"P3-H3_activation_order.{ext}"), dpi=300)
    plt.close(fig)
    log.info("wrote figures/P3-H3_activation_order.{pdf,png}")


if __name__ == "__main__":
    main()
