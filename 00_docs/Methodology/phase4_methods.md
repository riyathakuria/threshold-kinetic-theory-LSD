# Phase 4 — Methods: Orthogonal Evidence Integration and Prioritization

## Overview
Phase 4 prioritizes the 8 conserved lysosomal microglial modules (and their 41
member genes) identified in Phases 1–3, using acid-sphingomyelinase
(SMPD1/ASM)-centered orthogonal evidence. **Biological scope note:** Phases 1–3
characterized Niemann–Pick type C (NPC1/NPC2, cholesterol trafficking); the
Phase 4 ASM axis (SMPD1) is the Niemann–Pick A/B gene. The two are integrated as
*orthogonal* lysosomal-disease evidence, not as the same disorder.

Eight evidence layers were assembled: three **internal** (from the project's own
Phase 1–3 results) and five **external** (public databases, accessed 2026-07-19).

## Evidence layers

### Internal (Phase 1–3 discovery)
- **A1 Conservation / coherence.** Within-module human microglial co-expression
  coherence at each module's most-active DAM state (HuMicAtlas v2.0.0), plus
  1:1 mouse ortholog coverage (100% all modules).
- **A2 NPC disease kinetics.** Per-module NPC remodelling in the GSE152158 bulk
  time course (primary axis), summarized as peak\|Δ\| deviation from WT baseline,
  change-point onset week, direction, and a 1000-iteration bootstrap support.
- **A3 Developmental dynamics.** Fraction of member genes classified Dynamic vs
  Stable by NormFinder across the developmental atlas (dynamic_ratio).

### External (Phase 4 ASM-centered)
- **E1 SMPD1 perturbation (soft-fail).** Overlap of module genes with the
  SMPD1 CRISPR-KO consensus signature (LINCS L1000, Enrichr). 0/41 module genes
  overlapped → null layer; recorded as attempted, contributes zero weight.
- **E2 ASM functional network.** Consolidated single score from STRING (SMPD1↔
  gene edge confidence), Reactome (shared pathways with SMPD1), and KEGG (shared
  specific pathways), to avoid double-counting three correlated network sources.
- **E3 Regulatory.** Overlap of each module's enriched upstream TFs (ChEA3) with
  SMPD1's upstream ChIP-X regulators; Fisher odds ratio.
- **E4 Brain expression.** Ordinal detection concordance in brain (Human Protein
  Atlas, single-nucleus + regional). Positive control.
- **E5 Disease association.** Max Open Targets association across an LSD/
  neurodegeneration disease panel + log-scaled ClinVar pathogenic-variant count.

## APPROVED scoring rules (user-approved 2026-07-19, before any ranking)

Layer roles were assigned by **measured discriminating power** (coefficient of
variation of the raw metric across the 8 modules): layers unable to discriminate
are not folded into a diluting average.
- **E4 brain expression** (CV 0.035, range 0.92–1.00) → **eligibility gate**, not
  a weighted term. All 8 modules / 41 genes pass; none penalized.
- **E1 perturbation** (null, soft-fail) → **reported, zero weight**.
- **E3 regulatory** (all 8 significant, weak discriminator) → **down-weighted ×0.5**.

### Framework 1 — MODULE prioritization
Each weighted layer is min–max normalized to [0,1] across the 8 modules. Two
blocks are computed separately, then combined:

- **Internal block** = mean of: A1 coherence; A2 (peak\|Δ\| × bootstrap support);
  A3 dynamic_ratio — each min–max normalized first.
- **External block** = weighted mean of: E2 network (w 1.0); E3 regulatory
  (w 0.5); E5 disease (w 1.0) — each min–max normalized first.

**Combine (APPROVED):** `ModuleScore = 0.60 · InternalBlock + 0.40 · ExternalBlock`.
Both blocks are reported alongside the combined score.

### Framework 2 — GENE prioritization (separate)
Per-gene direct evidence, min–max normalized to [0,1] across the 41 genes:
- `G_proximity` = STRING SMPD1 edge score + Reactome/KEGG shared-pathway flags
- `G_disease` = per-gene OT LSD/neuro score + ClinVar pathogenic count
- `G_regulatory` = shared SMPD1 upstream TFs / gene's regulator count

**Combine (APPROVED):** `GeneScore = 0.45·G_proximity + 0.40·G_disease + 0.15·G_regulatory`.
E1 perturbation membership retained as an inert provenance annotation.
Developmental stability_status and the parent module's NPC evidence_call are
carried as annotation, not scored into the rank.

### Normalization
Min–max per layer across the scored set (approved; z-score/rank alternatives not used).

### Sensitivity analysis (APPROVED requirement)
Module weighting swept over **50/50, 60/40, 70/30** internal:external. Report
Spearman rank correlation between every pair of schemes and flag any module whose
priority tier changes across schemes. Gene weights held at the approved values.

### Provenance and coverage (APPROVED requirement)
Final tables report, for every evidence layer: source + version + access date,
per-layer coverage (fraction of modules/genes with a non-missing value), and
explicit missingness handling. Structurally undefined cells (e.g. single-gene
module M8 coherence) are recorded as such, not imputed.
