# Phase 4 — reproducibility record

Phase 4 integrates external, ASM(SMPD1)-centered evidence with our internal
Phase 1–3 discoveries. Reproducibility has two tiers because two kinds of step
were run.

## Tier 1 — network-interactive acquisition (pinned, not auto-re-run)

These layers queried live, versioned web resources on **ACCESS_DATE = 2026-07-19**
(`scripts/_phase4_config.py`). Re-querying at a later date would silently change
results (databases update), so the **saved checkpoints + per-layer tables are the
authoritative reproducible record**, not a re-execution. Each layer's method and
source is documented in the matching `docs/checkpoint{N}_summary.md`.

| Layer | Source(s) | Query anchor | Output table | Checkpoint | Notes |
|---|---|---|---|---|---|
| A1/A2/A3 internal | Phase 1–3 tables (offline) | — | `internal_evidence_scores.csv` | `phase4_internal.pkl` | deterministic from prior phases |
| E1 SMPD1 perturbation | CRISPR-KO signature (491-gene union) | SMPD1 | `asm_perturbation_module_scores.csv`, `asm_perturbation_gene_membership.csv` | `phase4_checkpoint1_perturbation.pkl` | **NULL: 0/41 overlap, all p_FDR=1.0 → zero weight** |
| E2 ASM functional network | STRING + Reactome + KEGG (one layer) | SMPD1 → 9606.ENSP00000340409 | `asm_network_module_scores.csv`, `asm_network_gene_scores.csv` | `phase4_checkpoint2_network.pkl` | M10 top (0.863) |
| E3 regulatory | ChEA3 vs SMPD1's 107 upstream regulators | universe = 406 | `regulatory_module_scores.csv`, `regulatory_gene_scores.csv`, `regulatory_module_top_tfs.csv` | `phase4_checkpoint3_regulatory.pkl` | all 8 sig; **down-weighted ×0.5** |
| E4 brain expression | Human Protein Atlas (GTEx rate-limited) | 42 genes | `brain_expression_module_scores.csv`, `brain_expression_gene.csv` | `phase4_checkpoint4_brain.pkl` | 41/41 expressed → **eligibility gate, unweighted** |
| E5 disease association | Open Targets + ClinVar | SMPD1 (assoc 1997; 502 pathogenic) | `disease_association_module_scores.csv`, `disease_association_gene_scores.csv` | `phase4_checkpoint5_disease.pkl` | M10 top (0.601) |

To re-acquire a layer from scratch, re-run its documented queries against the same
resources and expect drift; then rebuild the consolidated matrices below. The
checkpoints let you reproduce every downstream number without any network access.

## Tier 2 — deterministic integration & presentation (fully re-runnable offline)

| Step | Script | Consumes | Produces |
|---|---|---|---|
| Consolidate | (interactive, saved) | all Tier-1 tables | `consolidated_module_evidence_matrix.csv`, `consolidated_gene_evidence_matrix.csv` |
| 14 Integrate | `scripts/14_integrate_prioritize.py` | consolidated matrices | `module_prioritization_ranked.csv`, `gene_prioritization_ranked.csv`, provenance/coverage tables; sensitivity (ρ 0.929 / 0.810 / 0.952; tier-changers M6, M10) |
| 15 Figures | `scripts/15_phase4_figures.py` | 14's `score_modules()` + ranked tables | `P4-1 … P4-5` (PNG + PDF) |

`scripts/15_phase4_figures.py` reuses `14_integrate_prioritize.score_modules()`
as the single source of truth for the min–max normalized layers, so figures and
rankings can never diverge.

## Orchestrator

```
python scripts/run_all_phase4.py            # verify Tier-1 layers, then run 14 → 15
python scripts/run_all_phase4.py --dry-run  # print plan + layer presence check
```

It fail-fasts if any Tier-1 checkpoint/table is missing, prints the **E1 null**
and **E4 gate** handling to the run log, then runs the deterministic tier. Log →
`logs/run_all_phase4.log`.

## Approved scoring framework (frozen)

- Modules: **0.60 InternalBlock + 0.40 ExternalBlock**; Internal = mean(A1, A2
  peak|Δ|×bootstrap, A3 dynamic frac); External = weighted mean(E2 ×1.0, E3 ×0.5,
  E5 ×1.0). E4 = gate; E1 = reported null, zero weight.
- Genes: **0.45 proximity + 0.40 disease + 0.15 regulatory**.
- Min–max normalization per layer. Sensitivity over 50/50, 60/40, 70/30.
- Approved by user before ranking (`docs/phase4_methods.md`,
  `docs/phase4_scoring_framework_PROPOSAL.md`).

## Environment

conda env `humica` (Python 3.11, scanpy 1.11.5, ruptures, pandas, matplotlib).
No R/Rscript and no pyarrow in this env → checkpoints are pickle, not parquet.
