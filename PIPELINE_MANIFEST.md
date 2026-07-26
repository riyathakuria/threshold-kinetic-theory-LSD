# Phase 2 — Pipeline Manifest

Human microglial validation of 41 confirmed lysosomal candidate genes in
HuMicAtlas v2.0.0 (module-centric reference for Phase 3). Every output below is
regenerated deterministically (`RANDOM_SEED = 0`) by `scripts/06_run_all.py`.

## Scripts (execution order)

| Step | Script | Reads | Writes |
|------|--------|-------|--------|
| 01 | `scripts/01_load_checkpoint.py` | `checkpoints/HuMicAtlas_validated.h5ad` (read-only, MD5-verified) | `checkpoints/phase2_candidates.h5ad` |
| 02 | `scripts/02_expression_analysis.py` | candidate h5ad | `tables/phase2_expression_by_state.csv`, `phase2_expression_mean_matrix.csv`, `phase2_gene_expression_summary.csv` |
| 03 | `scripts/03_state_assignment.py` | candidate h5ad | `tables/phase2_state_assignment.csv`, `phase2_state_composition.csv` |
| 04 | `scripts/04_module_conservation.py` | candidate h5ad | `tables/phase2_module_activity_by_state.csv`, `phase2_module_coherence.csv`, `phase2_module_conservation.csv` |
| 05 | `scripts/05_generate_figures.py` | `tables/phase2_*.csv` | `figures/H1-H4` (PNG + PDF) |
| 06 | `scripts/06_run_all.py` | — | orchestrates 01→05 |

Shared modules (not run directly): `_phase2_config.py` (all parameters, paths,
state maps), `_figstyle.py` (vendored publication figure-style helpers),
`_env.py` (numba cache setup).

## Tables (`tables/`)

| File | Contents |
|------|----------|
| `phase2_expression_by_state.csv` | Long gene × state: mean/median/SD SCT residual |
| `phase2_expression_mean_matrix.csv` | Wide gene × state mean-residual matrix |
| `phase2_gene_expression_summary.csv` | Per-gene preferred state + across-state range |
| `phase2_state_assignment.csv` | Per-gene argmax state, group, assignment margin, NormFinder status |
| `phase2_state_composition.csv` | Stable/Dynamic counts per preferred state |
| `phase2_module_activity_by_state.csv` | Module × state score_genes activity |
| `phase2_module_coherence.csv` | Per-module mean off-diagonal Spearman (all cells + top state) |
| `phase2_module_conservation.csv` | Per-module peak activity, most-active state, coherence, Stable/Dynamic composition |

## Figures (`figures/`, 300 dpi PNG + vector PDF, Arial)

| File | Contents |
|------|----------|
| `H1_gene_state_heatmap` | Gene × state residual heatmap, rows grouped by STRING module, stability side-strip |
| `H2_module_activity_heatmap` | Module × state activity (score_genes, mean SCT residual) |
| `H3_module_conservation` | Within-module coherence (all vs top state) + Stable/Dynamic composition |
| `H4_module_synthesis` | Module-centric synthesis dashboard — the Phase 3 reference figure |

## Documentation (`docs/`)

| File | Purpose |
|------|---------|
| `phase2_methods.md` | Full methods: scale decision, parameters, threshold justifications, software |
| `phase2_integrated_interpretation.md` | Module-level interpretation, three module classes, per-module Phase 3 bridge |
| `phase2_decision_memo.md` | Residual-scale decision record (pre-registered deviation) |

## Metadata (`data/metadata/`)

`phase2_candidate_genes.csv` (41 genes, STRING cluster, NormFinder status),
`state_label_map.csv` (V1_clusters → state), plus checkpoint structure exports.

## Provenance guarantees

- Master checkpoint MD5 `ccc6c9d77cb46ab29f4329e2465d1cf6`, opened read-only,
  verified before and after subsetting; never modified.
- All analysis on the native **mean SCT Pearson residual** scale (no counts /
  log-norm slot exists in the object).
- No percent-expressing, detection proxy, or Tau (documented deviation).
- STRING MCL clusters are never collapsed; clusters 2 & 3 absent by
  construction; single-gene module (8/A2M) coherence reported as not defined.

---

# Phase 3 — Pipeline Manifest

NPC disease kinetics of the conserved lysosomal microglial modules (mouse).
Every output is regenerated deterministically (`RANDOM_SEED = 0`) by
`scripts/run_all_phase3.py`, which runs 07a→12 and prints the QC-gate axis
decision after step 08.

## Scripts (execution order)

| Step | Script | Reads | Writes |
|------|--------|-------|--------|
| 07a | `scripts/07a_acquire_npc.py` | GEO (GSE221609, GSE152158) | `data/raw/npc/*`, `data/metadata/ortholog_map.csv`, `tables/ortholog_coverage.csv` |
| 07 | `scripts/07_prepare_npc.py` | GSE221609 10x triplet | `checkpoints/npc_microglia.h5ad`, `tables/npc_celltype_marker_scores.csv`, `npc_microglia_marker_confirm.csv`, `npc_microglia_counts.csv`, `figures/P3_npc_umap_qc` |
| 07b | `scripts/07b_prepare_bulk.py` | GSE152158 norm-count CSV | `tables/bulk_expression_normalized.csv`, `bulk_sample_metadata.csv`, `figures/P3_bulk_qc` |
| 08 | `scripts/08_pseudotime.py` | `npc_microglia.h5ad`, bulk tables | `tables/pseudotime_values.csv`, `qc_gate_metrics.csv`, `figures/P3-H1_trajectory` — **QC GATE / axis choice (HARD STOP)** |
| 09 | `scripts/09_project_modules.py` | bulk + microglia, `ortholog_map.csv` | `tables/module_bulk_timecourse.csv`, `module_pseudotime_scores.csv`, `module_mapping_report.csv` |
| 10 | `scripts/10_module_kinetics.py` | `module_bulk_timecourse.csv` | `tables/module_kinetics.csv`, `module_onset_distinguishability.csv`, `figures/P3-H3_activation_order` |
| 11 | `scripts/11_stability_remodel.py` | kinetics + Phase-2 stability | `tables/stability_vs_pseudotime.csv`, `stability_remodel_stats.csv`, `figures/P3-H4_stable_dynamic` |
| 12 | `scripts/12_generate_figures.py` | all Phase-3 tables | `tables/human_npc_transition_map.csv`, `module_summary_table.csv`, `final_module_evidence_table.csv`; `figures/P3-H2,H5,H6,H7` |
| — | `scripts/run_all_phase3.py` | — | orchestrates 07a→12 with QC-gate/STOP branch, logs to `logs/run_all_phase3.log` |

Shared modules: `_phase3_config.py` (all parameters, paths, marker sets,
disease-axis metadata, change-point settings), `_figstyle.py`, `_env.py`.

## Figures (`figures/`, 300 dpi PNG + vector PDF, Arial)

| File | Contents |
|------|----------|
| `P3-H1_trajectory` | Disease-axis trajectory / microglia UMAP |
| `P3-H2_module_activity` | 2×4 small-multiples of module Δ(week) with PELT onset + bootstrap-CI shading |
| `P3-H3_activation_order` | Module activation-order heatmap (contemporaneous, not force-ranked) |
| `P3-H4_stable_dynamic` | Human Stable/Dynamic character vs mouse remodelling |
| `P3-H5_conserved_synthesis` | %Dynamic vs peak\|Δ\| bubble (size = bootstrap support, colour = evidence call) |
| `P3-H6_module_timeline` | Horizontal week timeline of remodeler onsets → end, ● onset ▼ peak |
| `P3-H7_framework` | Integrated Phase 1→2→3 systems-biology schematic |
| `P3_npc_umap_qc`, `P3_bulk_qc` | Preprocessing QC panels |

## Documentation (`docs/`)

| File | Purpose |
|------|---------|
| `npc_preprocessing.md` | GSE221609 + GSE152158 preprocessing / microglia annotation |
| `phase3_qc_decision.md` | QC gate metrics and metric-driven axis selection |
| `phase3_methods.md` | Full Phase-3 methods (datasets → figures) |
| `phase3_interpretation.md` | Module-level biological interpretation |
| `phase3_story.md` | One-page Discussion narrative (Phases 1–3) |

## Provenance guarantees

- Master Phase-2 checkpoint untouched; Phase 3 works only on NPC-derived data.
- No invented thresholds; the disease axis is chosen by the QC gate
  (`qc_gate_metrics.csv`), not by preference.
- Human→mouse orthologs are 1:1, 100% coverage; STRING clusters never
  collapsed/renamed.
- All timing on the disease-deviation signal Δ(week)=NPC−WT; onsets with
  overlapping CIs are contemporaneous, not force-ranked; association not causation.
