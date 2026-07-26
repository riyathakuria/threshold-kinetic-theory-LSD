# LSD *in silico*

**Conserved lysosomal microglial modules: discovery → human validation →
NPC disease kinetics → orthogonal ASM-centered prioritization.** A four-phase
module-centric framework tracing a curated set of lysosomal-storage-disease
candidate genes from STRING functional modules through human-microglia validation,
mouse NPC disease-progression kinetics, and finally a prioritization against
orthogonal acid-sphingomyelinase (SMPD1)-centered evidence.

## What Phase 2 does

Phase 1 nominated a set of lysosomal-relevant candidate genes and organized them
into STRING MCL functional modules and NormFinder developmental classes
(Stable / Dynamic). Phase 2 asks a single validation question against a fixed,
checksum-verified human microglial atlas:

> Across the nine published adult human microglial states, how are these
> candidate genes and their parent functional modules expressed — and which
> modules behave as state-focused programs versus broad scaffolds?

The answer is expressed at the **module level** (genes illustrate their parent
module; they are not interpreted individually), producing a reference against
which Phase 3 can test disease-driven movement.

## Key result

The 41 v2-present candidates fall into **three module classes**:

- **Class A — state-focused Dynamic programs** (Translation Elongation →
  Ribo.DAM2; Sphingolipid metabolism → Lipo.DAM; Unfolded-protein response →
  DIMs). Enriched for Dynamic genes and locally coherent; the movement
  candidates for Phase 3.
- **Class B — broad Stable scaffolds** (EGFR, Rho, PKA/glucagon, Apoptosis
  modules). Distributed across states; the negative-control backbone for
  Phase 3.
- **Class C — under-determined** (single-gene module, A2M).

Gene-level state preferences are mostly weak (median assignment margin 0.058;
26/41 genes < 0.10), and STRING modules show only modest single-cell
co-expression coherence — both reported honestly and used to calibrate Phase 3
expectations.

## Data & scale

- Fixed input: `checkpoints/HuMicAtlas_validated.h5ad` (90,716 cells × 3,000
  HVGs; MD5 `ccc6c9d77cb46ab29f4329e2465d1cf6`). Read-only; never modified.
- The object stores **SCT Pearson residuals only** (no counts, no log-norm
  slot). All Phase 2 quantities are computed on the native **mean SCT Pearson
  residual** scale. Percent-expressing / detection proxies / Tau are **not**
  computed (documented deviation — see `docs/phase2_methods.md`).
- States from Martins-Ferreira et al. 2025, *Nat Commun* 16:739
  (doi:10.1038/s41467-025-56124-1).

## Layout

```
scripts/    01-06 pipeline (06_run_all.py orchestrates); _phase2_config, _figstyle, _env
tables/     phase2_*.csv  (expression, state assignment, module conservation)
figures/    H1-H4  (PNG + vector PDF, Arial, 300 dpi)
docs/       phase2_methods.md, phase2_integrated_interpretation.md, phase2_decision_memo.md
data/metadata/  candidate genes, state label map, checkpoint structure exports
checkpoints/    HuMicAtlas_validated.h5ad (master), phase2_candidates.h5ad (derived)
```

## Reproduce

```bash
conda activate humica
python scripts/06_run_all.py        # Phase 2 (human validation reference)
python scripts/run_all_phase3.py    # Phase 3 (NPC disease kinetics, 07a->12)
python scripts/run_all_phase4.py    # Phase 4 (verify evidence layers -> integrate 14->15)
```

See `LOCAL_EXECUTION.md` for prerequisites and partial-run options,
`PIPELINE_MANIFEST.md` for the file-by-file inventory, and
`docs/phase2_methods.md` / `docs/phase3_methods.md` for full methods.

## Phase 3 — NPC disease kinetics of the conserved modules (mouse)

Projects the eight Phase-2 STRING modules onto an NPC disease axis and uses
change-point detection to locate *when* each remodels.

- **Axis chosen by a QC gate, not preference.** The single-cell diffusion-
  pseudotime axis **failed** the disease-relevance gate (DPT-vs-DAM ρ=0.069;
  Npc1_KO not shifted, p=1.0); the **bulk real-time** course (GSE152158, weeks
  1/3/6/9) **passed** (DAM-vs-week ρ=0.746) and is **primary**.
- **Four modules remodel, all contemporaneous (~week 3):** M6 Translation
  Elongation (↑ persistent, **Strong**, the anchor), M7 Apoptosis (↑ transient,
  **Strong**), M4 EGFR (↓, Moderate), M1 PKA (↑ late, Moderate). M8/A2M is a
  large but replicate-noisy Ambiguous; M5/M9/M10 are Stable. Onset CIs overlap
  (0/6 pairs distinguishable) → **not** force-ranked into a cascade.
- **Human Dynamic → mouse remodelling** is a positive but underpowered trend
  (ρ=0.42, p=0.30 at n=8). Temporal claims are **association, not causation**.

Full result: `docs/phase3_interpretation.md`; one-page Discussion:
`docs/phase3_story.md`; ledger: `tables/final_module_evidence_table.csv`.

## Phase 4 — orthogonal ASM(SMPD1)-centered prioritization

Prioritizes the 8 modules and 41 genes using two evidence blocks — *internal*
(Phase 1–3 coherence, NPC kinetics, dynamics) and *external* ASM-centered
(perturbation, functional network, regulatory, brain-expression gate, disease
association) — under a **user-approved** scoring framework (modules 60/40
internal:external; genes 0.45/0.40/0.15; min–max per layer). The one direct
SMPD1-perturbation layer returned a **clean null** (0/41 overlap → zero weight).

- **Central finding:** internal disease dynamics and external ASM centrality are
  **orthogonal**. **M10 Sphingolipid metabolism** is the most ASM-proximal module
  (network + disease) yet is a Phase-3 **non-remodeler**; the top remodeler **M6
  Translation Elongation** is functionally distant from ASM.
- **Top genes:** the M10 sphingolipid enzymes **ASAH1, PSAP, UGT8** — the most
  defensible single-gene follow-up targets.
- **Caveat carried in every table:** M8/A2M ranks #1 on internal structure alone
  (single-gene outlier), but Phase 3 called it a non-remodeler — "large but
  low-confidence", not a robust target.
- Rankings robust across weightings (Spearman ρ 0.81–0.95); only M6/M10 cross a
  tier. All claims associational.

Full result: `docs/phase4_interpretation.md`; one-page Phases 1–4 Discussion:
`docs/phase4_story.md`; reproducibility: `docs/REPRODUCIBILITY_phase4.md`.

## Status

**All four phases complete:** discovery (P1) → human validation (P2) → disease
kinetics (P3) → orthogonal prioritization (P4). See `PROJECT_STATE.md` for the
full inventory and recommended future work.
