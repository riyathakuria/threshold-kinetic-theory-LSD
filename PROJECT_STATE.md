# PROJECT STATE — LSD *in silico*

_Last updated: end of Phase 3 (awaiting Phase 4 approval)._

## Phase status

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Candidate nomination, STRING modules, NormFinder Stable/Dynamic classes | **Complete** (finalized) |
| 2 | Human microglial validation in HuMicAtlas v2.0.0 (module-centric) | **Complete** |
| 3 | NPC disease kinetics of the conserved module reference (mouse) | **Complete** |
| 4 | (proposed) powered temporal / mechanistic follow-up | **Not started — awaits approval** |

## Phase 3 summary

- **Question:** as NPC disease progresses in mouse microglia, which of the eight
  Phase-2 STRING modules remodel, and *when* along the disease course?
- **Datasets (mouse):** GSE221609 snRNA-seq (binary WT vs Npc1−/−, 2,656
  microglia) — single-cell **supporting** axis; GSE152158 bulk real-time course
  (weeks 1/3/6/9, 31 samples, 18,277 genes) — **primary** axis.
- **Ortholog map:** all 41 human candidates → mouse 1:1, **100% coverage** for
  all 8 modules (`data/metadata/ortholog_map.csv`, `tables/ortholog_coverage.csv`).
- **QC gate (metric-driven axis choice):** single-cell diffusion-pseudotime
  axis **FAILED** the disease-relevance gate (DPT-vs-DAM Spearman 0.069;
  Npc1_KO not shifted, p=1.0); bulk real-time axis **PASSED** (DAM-vs-week
  Spearman 0.746). **PRIMARY AXIS = bulk_realtime.** Pseudotime retained as a
  failing supporting axis. (`tables/qc_gate_metrics.csv`,
  `docs/phase3_qc_decision.md`.)
- **Kinetics:** ruptures PELT (`l2`, min_size=3, penalty log(n)·σ² with σ² =
  within-week replicate variance) on the disease-deviation signal
  Δ(week)=NPC−WT; bootstrap 95% CIs (n=1000).

### Result — four remodelers, all contemporaneous (~wk3)

| Module | STRING identity | Onset (95% CI) | Dir | Pattern | peak\|Δ\| | Evidence |
|---|---|---|---|---|---|---|
| M6 | Eukaryotic Translation Elongation | wk3 (3–6) | ↑ | **persistent** | **1.21** | **Strong** |
| M7 | Apoptosis | wk3 (3–3) | ↑ | transient | 0.54 | **Strong** |
| M4 | EGFR signaling | wk3 (3–6) | ↓ | transient | 0.46 | Moderate |
| M1 | PKA/glucagon | wk6 (3–6) | ↑ | transient | 0.41 | Moderate |
| M8 | Blood microparticle (A2M) | — | (↑ noisy) | — | 2.57 | Ambiguous |
| M5 / M9 / M10 | Rho / UPR / Sphingolipid | — | — | — | ≤0.36 | Stable |

- **All four remodeler onset CIs overlap (0/6 pairs distinguishable) →
  contemporaneous at week resolution, NOT force-ranked into a cascade.**
- **M6 Translation Elongation** (100% Dynamic, Ribo.DAM2-focused in human) is the
  anchor: early, persistent, largest reliable effect, bootstrap 0.95.
- **Human Dynamic → mouse remodelling** is a positive but underpowered trend
  (Spearman dynamic_frac vs peak\|Δ\| ρ=0.42, p=0.30 at n=8;
  `tables/stability_remodel_stats.csv`).
- Temporal statements are **association, not causation**.

## Phase 2 summary

- **Input (fixed):** `checkpoints/HuMicAtlas_validated.h5ad` — 90,716 cells ×
  3,000 HVGs, MD5 `ccc6c9d77cb46ab29f4329e2465d1cf6`. Read-only, MD5-verified,
  never modified. Derived efficiency copy: `checkpoints/phase2_candidates.h5ad`
  (41 genes).
- **Candidates:** 41 v2-present genes = 20 Stable + 21 Dynamic, spanning 8 of
  10 STRING MCL clusters (clusters 2 & 3 have no v2-present members).
- **States:** 9 microglial states from Martins-Ferreira et al. 2025
  (*Nat Commun* 16:739), read from `obs['V1_clusters']`.
- **Scale decision:** object carries SCT Pearson residuals only (no counts /
  log-norm). All quantities computed on the native **mean SCT Pearson
  residual** scale. No percent-expressing / detection proxy / Tau — a
  pre-registered deviation from the approved plan
  (`docs/phase2_decision_memo.md`, `docs/phase2_methods.md`).

### Result — three module classes

- **Class A (state-focused Dynamic):** cl6 Translation → Ribo.DAM2 (100%
  Dynamic); cl10 Sphingolipid → Lipo.DAM (ASAH1, PSAP); cl9 UPR → DIMs (most
  coherent, 0.088; 78% Dynamic). → Phase 3 movement candidates.
- **Class B (broad Stable scaffolds):** cl4 EGFR, cl5 Rho, cl1 PKA/glucagon,
  cl7 Apoptosis. → Phase 3 negative-control backbone.
- **Class C (under-determined):** cl8 (single gene, A2M).
- **Quantitative cautions:** gene-level preferences weak (median margin 0.058;
  26/41 < 0.10); module coherence modest overall (0.017–0.061). Reported
  honestly; used to calibrate Phase 3 expectations.

## Deliverables (all saved as artifacts)

- **Figures:** H1 gene×state heatmap, H2 module×state activity, H3 coherence +
  composition, H4 module synthesis (Phase 3 reference). PNG + PDF, Arial, 300 dpi.
- **Tables (8):** `tables/phase2_*.csv` — expression by/across states, state
  assignment + composition, module activity/coherence/conservation.
- **Scripts (6):** `scripts/01`–`06` + `_phase2_config`, `_figstyle`, `_env`.
- **Docs:** `phase2_methods.md`, `phase2_integrated_interpretation.md`,
  `phase2_decision_memo.md`; `README.md`, `LOCAL_EXECUTION.md`,
  `PIPELINE_MANIFEST.md`, this file.

## Standing constraints (carried into Phase 3)

- Master checkpoint is authoritative and preserved unchanged; work on derived
  copies only.
- No invented thresholds; module-level primary interpretation (genes illustrate
  parent module).
- STRING MCL clusters never collapsed.
- Regenerate (don't patch) scripts; centralize parameters in `_phase2_config.py`.
- **Phase 3 (NPC / pseudotime / disease) not to begin without explicit
  approval.**

## Phase 3 deliverables (all saved as artifacts)

- **Scripts:** `scripts/07a`–`12` + `run_all_phase3.py` (orchestrator, 07→12,
  QC-gate/STOP branch), `_phase3_config.py`, `_figstyle.py`, `_env.py`.
- **Tables:** `pseudotime_values.csv`, `module_pseudotime_scores.csv`,
  `module_bulk_timecourse.csv`, `module_kinetics.csv`,
  `module_onset_distinguishability.csv`, `stability_vs_pseudotime.csv`,
  `stability_remodel_stats.csv`, `human_npc_transition_map.csv`,
  `module_summary_table.csv`, `final_module_evidence_table.csv`,
  `qc_gate_metrics.csv`, `ortholog_coverage.csv`, bulk QC tables.
- **Figures:** `P3-H1`..`P3-H7` (300 dpi PNG + vector PDF, Arial) +
  `P3_npc_umap_qc`, `P3_bulk_qc`.
- **Docs:** `npc_preprocessing.md`, `phase3_qc_decision.md`,
  `phase3_methods.md`, `phase3_interpretation.md`, `phase3_story.md`.

## Next step — STOP before Phase 4; await approval

Phase 3 is complete and self-contained. **Phase 4 is not to begin without
explicit approval.** The recommended Phase 4 (proposed, not started) would
address Phase 3's structural limits: (1) a **single-cell, genuinely
time-resolved** NPC microglial dataset to replace the failed pseudotime axis and
the bulk population average; (2) **denser time sampling** to resolve whether the
week-3 contemporaneous inflection separates into an ordered sequence; (3) a
**properly powered** test of the human-Dynamic → mouse-remodelling link (more
modules / independent cohorts); (4) targeted **mechanistic** follow-up on the
M6 Translation-Elongation anchor (association → causation). See the final
summary and `docs/phase3_story.md`.
