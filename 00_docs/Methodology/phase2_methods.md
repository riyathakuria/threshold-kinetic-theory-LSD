# Phase 2 — Methods

Human microglial validation of confirmed lysosomal candidate genes in
HuMicAtlas v2.0.0 (module-centric reference framework for Phase 3).

## Data source

The analysis uses a fixed, checksum-verified HuMicAtlas v2.0.0 microglial
reference object (`checkpoints/HuMicAtlas_validated.h5ad`; 90,716 cells ×
3,000 highly-variable Ensembl-indexed genes; MD5
`ccc6c9d77cb46ab29f4329e2465d1cf6`). Nine published microglial states are read
from `obs['V1_clusters']` and mapped to the HuMicA state nomenclature of
Martins-Ferreira et al. (2025, *Nat Commun* 16:739, doi:10.1038/s41467-025-56124-1):
0 Homeos1, 1 Inflam.DAM, 2 DIMs, 3 Ribo.DAM1, 4 Homeos2, 5 Ribo.DAM2,
6 Lipo.DAM, 7 MAC, 8 Homeos3. States are grouped as Homeostatic
(Homeos1/2/3), DAM (Inflam.DAM, Ribo.DAM1/2, Lipo.DAM), DIM (DIMs), and
Macrophage (MAC) for display only.

The 41 candidate genes are the Phase 1 confirmed lysosomal candidates present
in v2.0.0, spanning 8 of the 10 STRING MCL functional modules (clusters 2 and 3
have no v2-present members). Each gene carries a fixed NormFinder developmental
classification (20 Stable, 21 Dynamic); this label is used but never
re-derived or altered in Phase 2.

## Expression scale — decision and justification

On inspection, the fixed checkpoint stores the Seurat SCT `scale.data` slot
(Pearson residuals) in `.X`, `layers['SCT']`, and `.raw/X`. `.X` is clipped at
the Seurat `ScaleData` default maximum of 10 and ~2/3 of its values are
negative. The object contains **no raw counts and no log-normalized `data`
slot**; only per-cell count summaries (`nCount_SCT`, `nFeature_SCT`) survive in
`obs`.

Consequently, all Phase 2 quantities are computed directly on the **mean SCT
Pearson residual** scale (`_phase2_config.EXPRESSION_MATRIX = "X"`). We do
**not** compute percent-expressing, detection-fraction proxies, or Tau tissue-
specificity, because each requires a non-negative detection or count scale that
this object does not provide, and fabricating one (e.g. by min-shifting
residuals) would introduce an arbitrary transform outside the fixed input. This
is a deliberate decision to remain strictly within the provided object rather
than import an external count matrix, which would break the reproducibility
guarantee of the checksum-fixed checkpoint. SCT Pearson residuals are variance-
stabilized and are appropriate for the two operations Phase 2 performs:
cross-state relative-activity comparison and rank-based co-expression
estimation.

**Deviation from the pre-approved plan.** The approved plan listed a
percent-expressing readout (H2 dot plot) and a descriptive Tau specificity
index. Both were removed for the reason above and replaced by residual-native,
module-centric readouts. All other plan elements are unchanged.

## Analysis 1 — Per-gene expression across states (`02_expression_analysis.py`)

For each candidate gene and each state, the mean, median, and SD of the SCT
Pearson residual are computed over all cells assigned to that state. Outputs:
a long gene × state table, a wide gene × state mean-residual matrix, and a
per-gene descriptive summary reporting the argmax preferred state and the
across-state residual range (max − min of per-state means). No thresholding or
classification is applied.

## Analysis 2 — Descriptive state assignment (`03_state_assignment.py`)

Each gene is assigned its argmax preferred state (highest mean residual). An
**assignment margin** (top minus second-best state mean residual) is reported
as a transparent, non-index descriptor of preference strength; margins are
reported as continuous values and never binned into specific/broad classes. A
per-state composition table counts Stable vs Dynamic genes preferring each
state. No Tau or specificity index is computed (see scale decision).

## Analysis 3 — Module conservation (`04_module_conservation.py`)

Modules are the original STRING MCL clusters (not collapsed). Three readouts:

1. **Activity** — `scanpy.tl.score_genes` (control-gene-corrected mean residual
   of module members; `ctrl_size=50`, `random_state=0`) computed per cell, then
   averaged per state. Single-gene modules use the member residual directly.
2. **Coherence** — mean of the off-diagonal pairwise Spearman correlations
   among member-gene residual vectors, computed across all cells and,
   separately, within the module's most-active state. Spearman was chosen over
   Pearson because it is rank-based and therefore robust to the `ScaleData`
   clip at 10 and to residual outliers; on variance-stabilized residuals it
   estimates monotonic co-expression without distributional assumptions.
   Single-gene modules (cluster 8) yield no pairs and are reported as not
   defined.
3. **Composition** — Stable vs Dynamic member counts and percent-stable per
   module, from the fixed NormFinder labels.

Clusters 2 and 3 are absent by construction (no v2-present members) and are
reported as such rather than imputed.

## Figures (`05_generate_figures.py`)

Four figures, generated with a vendored subset of the publication figure-style
helpers (`scripts/_figstyle.py`): Arial, open frame, 300 dpi PNG + vector PDF.
Diverging heatmaps (`RdBu_r`) are symmetric about zero (semantic-zero centring)
since residuals are signed. H1 = gene × state residual heatmap grouped by
module; H2 = module × state activity heatmap; H3 = within-module coherence
(all-cells vs top-state) and Stable/Dynamic composition; H4 = module-centric
synthesis dashboard (peak activity, preferred state, representative gene,
composition, coherence) serving as the Phase 3 reference figure. Colour maps:
preferred-state group (Homeostatic/DAM/DIM/Macrophage) and developmental status
(Stable/Dynamic, Okabe–Ito colour-blind-safe).

## Reproducibility

All parameters are centralized in `scripts/_phase2_config.py` (paths, state
maps, STRING cluster names, `RANDOM_SEED = 0`, `COHERENCE_METHOD = "spearman"`,
expression-scale constants). The master checkpoint is opened read-only and its
MD5 is verified before and after subsetting; Phase 2 operates only on a derived
41-gene efficiency copy (`checkpoints/phase2_candidates.h5ad`) and never
modifies the master. Scripts are run in order via `06_run_all.py`.

## Software

Python 3.11.15 · scanpy 1.11.5 · anndata 0.12.19 · numpy 2.4.6 · scipy 1.17.1 ·
pandas 2.3.3 · matplotlib 3.11.0 · h5py 3.16.0 · scikit-learn 1.9.0. Conda
environment `humica`. `NUMBA_CACHE_DIR` is set to a writable path before scanpy
import because the conda environment is mounted read-only.
