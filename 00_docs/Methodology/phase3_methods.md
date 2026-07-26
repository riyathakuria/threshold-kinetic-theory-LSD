# Phase 3 methods — NPC disease kinetics of conserved lysosomal microglial modules

Phase 3 asks *when*, along Niemann-Pick type C (NPC) disease progression in
mouse microglia, each Phase-2 STRING functional module remodels. It projects
the eight Phase-2 STRING MCL modules onto an NPC disease axis and applies
change-point detection to locate the onset of remodelling. All interpretation
is at the **module level** (genes illustrate their parent module) and every
temporal statement is **association, not causation**.

All parameters, paths, marker sets, and change-point settings are centralized
in `scripts/_phase3_config.py`; nothing is invented inline. Scripts are
numbered `07`–`12` and orchestrated by `scripts/run_all_phase3.py`.

---

## 1. Datasets (both mouse, *Mus musculus*)

| GEO | Assay | Design | Role |
|---|---|---|---|
| GSE221609 | snRNA-seq (10x) | Npc1−/− vs WT forebrain, **binary genotype** (3 vs 3) | Single-cell **supporting** axis |
| GSE152158 | bulk RNA-seq, FACS-purified microglia | **Real-time** NPC disease course, weeks 1/3/6/9 | **Primary** disease axis |

GSE221609 is an oligodendrocyte-focused study; microglia are a minority subset
that must be extracted. GSE152158 supplies genuine elapsed time.

## 2. Human→mouse ortholog mapping (Step 07a metadata)

All 41 Phase-2 human candidate genes map 1:1 to mouse orthologs
(`tables/ortholog_coverage.csv`): **100% coverage for all 8 modules**.
Orthologs were resolved via Ensembl; two required manual assignment
(`HSPA8`→`Hspa8`, `UGT8`→`Ugt8a`). The map is frozen in
`data/metadata/ortholog_map.csv`. Because coverage is complete, no module is
penalized or flagged for low ortholog detection.

## 3. GSE221609 snRNA-seq preprocessing & microglia annotation (Step 07)

- **Input:** 75,878 nuclei × 32,285 genes (cellranger-filtered).
- **QC:** `min_genes/cell = 200`, `min_cells/gene = 3`, `pct_counts_mt < 5.0%`
  (mito prefix `mt-`). Post-QC: 75,805 nuclei × 27,127 genes. The low mito
  ceiling suits single-nucleus data (little cytoplasmic mito RNA); pre/post
  distributions are shown in `figures/P3_npc_umap_qc`.
- **Normalization / clustering:** CP10k + log1p (counts kept in
  `layers['counts']`); HVG = 2000 (seurat); scale(max=10); PCA 30; neighbors
  k=15; UMAP; Leiden resolution 1.0 (seed 0) → 35 clusters.
- **Genotype:** assigned from barcode suffix (−1/−2/−3 = WT;
  −4/−5/−6 = Npc1_KO).
- **Microglia annotation:** clusters scored against seven cell-type marker sets
  (`sc.tl.score_genes`, ctrl_size=50, seed 0); microglia clusters then
  **confirmed** by canonical pan-microglial marker expression, not accepted on
  de novo score alone. Confirmation (mean log1p CP10k): Hexb 1.97, Cx3cr1 1.29,
  Csf1r 1.23, P2ry12 1.02 in microglia vs ≤0.15 elsewhere.
- **Yield:** **2,656 microglia** (WT 1,473 + Npc1_KO 1,183).
- **Output:** `checkpoints/npc_microglia.h5ad`; marker/count tables.

Full detail: `docs/npc_preprocessing.md`.

## 4. GSE152158 bulk preprocessing (Step 07b)

The supplementary `GSE152158_NormCountData.csv.gz` is **already** normalized, so
TMM/CPM was **not** re-applied (would double-normalize); a log2(x+1) transform
was applied for scoring. The disease-course samples WT/NPC × weeks {1,3,6,9}
were retained (31 samples: WT 4/4/4/4, NPC 4/4/4/3); the 10 CSF1R-inhibitor
arms (`NPC-PLX*`, `NPC-Con*`) were **excluded**. Matrix after dropping all-zero
genes: **18,277 genes × 31 samples**. Output:
`tables/bulk_expression_normalized.csv`, `tables/bulk_sample_metadata.csv`,
`figures/P3_bulk_qc`.

## 5. Diffusion pseudotime + QC gate with hard STOP (Step 08)

A within-microglia diffusion pseudotime (DPT) was built on GSE221609 as a
candidate single-cell disease axis, then subjected to a pre-registered QC gate
(`tables/qc_gate_metrics.csv`) that decides *from metrics, not arbitrary
cutoffs* whether that axis is defensible as PRIMARY.

| Metric | Value | Bar | Pass |
|---|---|---|---|
| Microglia N | 2,656 | ≥ 500 | ✅ |
| kNN connected components | 1 | = 1 | ✅ |
| Diffusion eigen-gap | 0.0445 | > 0 | ✅ |
| Spearman DPT vs DAM signature | **0.069** | \|ρ\| ≥ 0.30 | ❌ |
| Spearman DPT vs homeostatic | 0.176 | <0 expected | ❌ |
| Genotype shift (KO>WT, MWU) | **p = 1.00** | p < 0.05 | ❌ |
| Bulk time points | 4 | ≥ 3 | ✅ |
| Bulk DAM-vs-week Spearman | **0.746** | >0, p<0.10 | ✅ |

**Decision.** The single-cell DPT axis **fails** the disease-relevance gate:
it does not track the homeostatic→DAM gradient (ρ=0.069) and Npc1_KO cells are
not shifted toward the DAM end (p=1.0). The bulk real-time axis **passes**
(DAM signature rises with disease week, ρ=0.746). Because GSE221609 is a binary
genotype contrast, its DPT is an activation-state continuum, not elapsed time —
the wrong substrate for temporal-onset claims regardless. **PRIMARY DISEASE
AXIS = bulk real-time (GSE152158, weeks 1/3/6/9).** The pseudotime axis is
retained as a *failing* SUPPORTING axis and reported transparently, never used
for the primary onset claims. Full rationale: `docs/phase3_qc_decision.md`.

## 6. Project Phase-2 modules onto the disease axis (Step 09)

Each STRING module (mouse orthologs) was scored with
`sc.tl.score_genes(ctrl_size=50, random_state=0)` — **identical** to Phase 2 —
on both axes. STRING module identities are preserved throughout (never
collapsed/renamed). 100% ortholog detection on both axes; no LOW_DETECTION
flags. Outputs: `tables/module_bulk_timecourse.csv` (31 samples × 8 modules +
week/genotype), `tables/module_pseudotime_scores.csv` (2,656 cells,
supporting), `tables/module_mapping_report.csv`.

## 7. Module kinetics — change-point, CIs, goodness-of-fit (Step 10)

**Disease-deviation signal (critical).** To separate NPC dynamics from
developmental/age effects shared with WT, all timing is computed on the
per-week **disease deviation**

  Δ(week) = mean_NPC(week) − mean_WT(week)

of the module score. Onset, peak, direction, and pattern are **all** derived
from this same Δ series (a prior bug where onset used raw NPC while peak used Δ
was fixed).

**Change-point (PRIMARY): ruptures PELT**, `model="l2"` (piecewise-constant
mean shift), `min_size=3`. Penalty = `log(n)·σ²` where **σ² is the within-week
replicate variance** (`_within_week_sigma2`), not the lag-1 difference variance
— using lag-1 differences would let a genuine step inflate its own penalty and
suppress its own detection. Onset week = first detected change-point on Δ; peak
week = argmax |Δ|; direction = sign of Δ at peak.

**Classification.** Timing: early if onset ≤ week 3, late if ≥ week 6. Pattern:
persistent if Δ(week 9) ≥ 70% of peak |Δ|, else transient.

**Confidence intervals: bootstrap** (`n=1000`, percentile 95% CI), resampling
replicates within each genotype×week cell. A module is called a confident
**remodeler** only if PELT finds an NPC-specific onset with bootstrap support.

**Distinguishability.** Pairwise onset-CI overlap (`distinguish()`): overlapping
CIs → onsets are **contemporaneous**, not force-ranked.

**Sensitivity:** pwlf piecewise-linear break week + R² reported alongside (not
primary).

**Results (`tables/module_kinetics.csv`).** Four confident remodelers:

| Module | Onset (95% CI) | Dir | Timing | Pattern | peak\|Δ\| | Bootstrap |
|---|---|---|---|---|---|---|
| **M6** Translation Elongation | wk3 (3–6) | up | early | **persistent** | **1.21** | 0.95 |
| **M7** Apoptosis | wk3 (3–3) | up | early | transient | 0.54 | 0.87 |
| **M4** EGFR | wk3 (3–6) | down | early | transient | 0.46 | 0.66 |
| **M1** PKA/glucagon | wk6 (3–6) | up | late | transient | 0.41 | 0.64 |

Non-remodelers: **M5** Rho, **M8** Blood microparticle (A2M — large peak\|Δ\|=2.57
but replicate-noisy, correctly fails the confidence bar), **M9** UPR, **M10**
Sphingolipid. **Distinguishability:** 0/6 remodeler pairs have non-overlapping
onset CIs → **all remodelers are contemporaneous at week resolution** (single
group G1), not temporally ordered. Outputs:
`tables/module_kinetics.csv`, `tables/module_onset_distinguishability.csv`,
`figures/P3-H3_activation_order`.

## 8. Stable vs Dynamic remodelling (Step 11)

MODULE-level cross-species test: does a module's **human** Phase-2 Dynamic
character predict how strongly it remodels along the **mouse** NPC axis?
`dynamic_frac` = fraction of a module's genes called Dynamic in human;
`module_class` = Dynamic-dominant (>0.5) vs Stable-dominant.

| Test | n | Statistic | p |
|---|---|---|---|
| Spearman dynamic_frac vs peak\|Δ\| | 8 | ρ=0.42 | 0.30 |
| Spearman dynamic_frac vs bootstrap_frac_remodel | 8 | ρ=0.08 | 0.84 |
| Mann–Whitney peak\|Δ\| (Dynamic>Stable) | 4v4 | U=10 | 0.34 |

The trend is **positive and consistent** with the conserved-program hypothesis
(top remodeler M6 is 100% Dynamic) but **non-significant / underpowered at
n=8** — reported as a consistency check, not a definitive test. Outputs:
`tables/stability_vs_pseudotime.csv`, `tables/stability_remodel_stats.csv`,
`figures/P3-H4_stable_dynamic`.

## 9. Integrated tables + figures (Step 12)

**Phase A — three cross-phase tables**, framed strictly HUMAN module property →
CONSERVED MOUSE module activity (never state→state):
`tables/human_npc_transition_map.csv`, `tables/module_summary_table.csv`,
`tables/final_module_evidence_table.csv`. The evidence table is the full ledger
(Phase-1 STRING identity + genes + ortholog coverage; Phase-2 stable/dynamic
composition, class, most-active human state, coherence; Phase-3 onset/CI/peak/
direction/timing/pattern/PELT-R²/bootstrap/contemporaneous group), with an
`evidence_call`: **Strong** (remodels & bootstrap ≥0.75), **Moderate**
(remodels, lower support), **Ambiguous** (peak\|Δ\|≥1.0 but replicate-noisy),
**Stable** (else). Result: M6+M7 = Strong; M4+M1 = Moderate; M8 = Ambiguous;
M5/M9/M10 = Stable.

**Phase B — seven hero figures** (`P3-H1`..`P3-H7`), 300 dpi PNG + vector PDF,
Arial: H1 trajectory/UMAP, H2 module-activity small-multiples of Δ with PELT
onset+CI shading, H3 activation-order heatmap, H4 Stable/Dynamic, H5
conserved-program synthesis bubble, H6 module timeline, H7 integrated
systems-biology framework schematic. Each was rendered then visually verified.

## 10. Software & reproducibility

Conda env `humica`: Python 3.11, scanpy 1.11.5, ruptures, pwlf, pybiomart.
Determinism: `random_state=0`/`seed=0` throughout; bootstrap seed fixed
(`RANDOM_SEED=0`). Numba cache pre-seeded (`scripts/_env.py`) before scanpy
import. Figure style vendored in `scripts/_figstyle.py`. The pipeline is
re-runnable end-to-end via `scripts/run_all_phase3.py` (see
`docs/phase3_qc_decision.md` for the axis-selection branch).

## 11. Standing constraints honored

Master checkpoint untouched; no invented thresholds; module-level primary
interpretation; STRING clusters never collapsed/renamed; parameters centralized;
scripts regenerated not patched; **human module → conserved mouse module**
framing only; temporal order = association not causation; overlapping CIs →
contemporaneous.
