# Phase 1 — Data Preparation & Validation (Methods)

**Project:** LSD in-silico — locating 109 candidate genes in the HuMicA v2.0.0 microglia atlas
**Scope of this phase:** data acquisition, atlas structural inspection, candidate gene
resolution, presence classification, and a validated metadata checkpoint. **No downstream
analysis** (subtype scoring, heatmaps, ML, pseudotime) is performed in Phase 1.

---

## 1. Dataset acquisition

- **Source:** Zenodo record **18458280** — *HuMicA v2.0.0* (published 2026-02-02).
- **File:** `HuMicA_v2.0.0.h5ad`, 6,803,451,711 bytes (6.803 GB).
- **Integrity:** MD5 verified against the Zenodo-published checksum
  `0630ebce364221c2d88abcddac50c0a2` — **PASS**. SHA-256 also recorded in
  `data/download_manifest.csv`.
- **Download method:** custom resumable downloader (`scripts/00_acquire_humica_v2.py`).
  The sandbox network proxy does not reliably honour HTTP `Range` boundaries, so a naive
  append-on-resume strategy duplicated bytes and produced an oversized, corrupt file. The
  final logic is **position-authoritative**: it reads the server's `Content-Range` response
  header, `seek()`s the local `.part` file to that exact offset (mode `r+b`), and refuses to
  write past the expected size. Overlapping bytes overwrite in place rather than growing the
  file, making the download correct regardless of proxy behaviour. Final assembled size and
  MD5 both matched.
- **Note on versions:** an earlier record, **14697727** (HuMicA v1, 2025-01-20), is *not*
  used. Only v2.0.0 is authoritative for this project.

## 2. Compute environment

- Execution in an Anthropic cloud sandbox (Linux), **not** the local Mac. Project files live
  on the user's iCloud Drive, mounted read-write.
- Conda env **`humica`** (Python 3.11): scanpy 1.11.5, anndata, h5py, numba, scipy, pandas,
  numpy, scikit-learn, openpyxl, leidenalg.
- **RAM constraint:** total RAM ≈ 8.6 GB < the 6.8 GB matrix, so the AnnData object is
  **never** fully materialised. All inspection uses `backed='r'` mode and direct `h5py`
  reads of index datasets. See `docs/environment_report.md`.
- **numba cache fix:** the conda env is a read-only mount, so numba's default on-disk cache
  raises at `import scanpy`. `scripts/_env.py::setup_numba_cache()` redirects
  `NUMBA_CACHE_DIR` to a writable project path (`data/external/.numba_cache`) before any
  numba/scanpy import.

## 3. Atlas structure (key finding)

`scripts/02_inspect_anndata.py` (backed mode) established:

- **90,716 cells × 3,000 genes.**
- The object is **HVG-reduced**: `.X`, the single `SCT` layer, and `.raw` all carry the
  **same 3,000 highly-variable Ensembl gene IDs** (verified set-equal, only ordering
  differs). **There is no full transcriptome in this object** — `.raw` is *not* a larger
  gene space here.
- Gene identifiers are **Ensembl gene IDs** (`ENSG…`) in all namespaces.
- `.obs` has 21 columns spanning 12 studies, 7 diagnostic groups, tissue metadata, and a
  `V1_clusters` label. Full detail in `docs/humica_structure_report.md`.

## 4. Candidate gene resolution (109 genes)

Multi-stage resolution (details in `tables/candidate_mapping.csv`), never a one-step
HGNC→Ensembl intersect:

1. **Stage 1 — direct HGNC symbol** (MyGene): 108/109 resolved.
2. **Stage 2 — alias/previous-symbol:** 1 case, **GBA → GBA1** (`ENSG00000177628`).
3. **Stage 3 — Ensembl BioMart:** confirmed/canonicalised all 109, resolving 7 genes that
   had multiple Ensembl IDs to their primary-assembly ID.
- Stages 4 (NCBI) and 5 (UniProt) were available as fallbacks but not required.
- `tables/manual_review_required.csv` is **empty** — no gene needed manual arbitration.
- All 109 candidates carry a Stable/Dynamic label and a STRING/PPI cluster assignment.

## 5. Presence classification

`scripts/03_candidate_presence.py` checks each resolved Ensembl ID against the atlas gene
namespaces (`.var`/`features`, `.var`/`SCT_features`, `.raw`/`_index`) **independently** —
never intersecting candidates with `.X` first, and never dropping a gene from the fixed
universe of 109.

**Result: 41 of 109 candidates are present; 68 are not detected.** Because the atlas is
HVG-reduced, "not detected" means the gene was not retained among the 3,000 highly-variable
features — it is a property of the deposited data, not a mapping failure (all 109 resolved
cleanly to Ensembl IDs). Since the three namespaces are identical, "present in `.X`" and
"present in `.raw`" coincide for every gene.

| Class | Present | Total | % present |
|---|---|---|---|
| Stable | 20 | 50 | 40% |
| Dynamic | 21 | 59 | 36% |
| **All** | **41** | **109** | **38%** |

Outputs: `tables/candidate_presence.csv` (per-gene slots + Stable/Dynamic + STRING cluster)
and `tables/candidate_mapping_audit.csv` (Original Symbol, Alias Used, Ensembl ID, Present
in var, Present in raw, Status, Reason). Figure: `figures/candidate_presence_overview.png`.

## 6. Validated checkpoint

`scripts/04_build_validated_checkpoint.py` → `checkpoints/HuMicAtlas_validated.h5ad`
(6.81 GB). Built RAM-safe (disk copy + in-place `h5py` metadata edits; matrix never loaded).
**Genes are NOT subset** — full 90,716 × 3,000 retained. Additions:

- **Corrected annotation:** `/obs/Sex` in the source mixed `Male/Female` with `F/M`;
  14,263 `F`/`M` cells were standardised to `Female`/`Male`. Original values preserved in
  `/obs/Sex_original`.
- **Candidate mapping attached to `.var`:** `lsd_is_candidate` (41 flagged),
  `lsd_candidate_symbol`, `lsd_stability_status`.
- **Provenance in `/uns/lsd_phase1`:** JSON provenance block plus the full candidate mapping
  and presence tables.

## 7. Limitations

- The atlas is HVG-reduced to 3,000 genes; **57 of 109 candidates cannot be studied at the
  expression level in this object** because they were not retained as HVGs. Any downstream
  expression analysis is necessarily restricted to the 41 present genes unless a
  full-transcriptome source is added in a later phase.
- The `SCT` layer is Pearson-residual normalised (Seurat SCTransform); `.X` and `.raw` share
  the same gene set. Which matrix to use for which downstream task is a Phase-2 decision.
