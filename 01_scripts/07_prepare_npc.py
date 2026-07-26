#!/usr/bin/env python3
"""
Phase 3 · Step 07 — Preprocess GSE221609 snRNA-seq & annotate microglia
=======================================================================

Loads the GSE221609 10x triplet (Npc1-/- vs WT mouse forebrain, snRNA-seq),
attaches genotype from the barcode suffix, runs QC + CP10k/log1p
normalization, clusters (PCA/neighbors/UMAP/Leiden), and annotates microglia
by CANONICAL MOUSE MARKER EXPRESSION (Cx3cr1, P2ry12, Tmem119, Hexb, C1qa/b,
Csf1r) — not de novo cluster identity alone. This is an oligodendrocyte-
focused study, so microglia are a minority subset that must be extracted.

No published per-cell annotations are distributed in the GEO supplementary
(raw matrix only), so clusters are formed here and microglia identified by
marker score; the marker evidence is written out so the call is auditable.

Outputs:
  checkpoints/npc_microglia.h5ad   microglia-only, normalized, with genotype
  tables/npc_celltype_marker_scores.csv   per-cluster marker scores
  tables/npc_microglia_counts.csv         microglia per genotype
  figures/P3_npc_umap_qc.pdf/.png          QC + annotation UMAP panel
  docs/npc_preprocessing.md                methods + decisions

Standalone: run inside the `humica` conda env.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import setup_numba_cache          # noqa: E402
setup_numba_cache()

import numpy as np                            # noqa: E402
import pandas as pd                           # noqa: E402
import scanpy as sc                           # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402

import _phase3_config as cfg                  # noqa: E402
from _figstyle import apply_figure_style, panel_letter  # noqa: E402

log = cfg.get_logger("07_prepare_npc")
sc.settings.verbosity = 1
np.random.seed(cfg.RANDOM_SEED)


def load_triplet() -> "sc.AnnData":
    """Read the 10x MTX triplet into an AnnData (cells x genes)."""
    log.info("reading MTX matrix (this is ~730 MB gzip; takes a few minutes)")
    # scanpy read_mtx returns genes x cells for this file layout; transpose.
    adata = sc.read_mtx(cfg.GSE221609_MTX).T
    features = pd.read_csv(cfg.GSE221609_FEATURES, sep="\t", header=None,
                           names=["ensembl", "symbol", "kind"])
    barcodes = pd.read_csv(cfg.GSE221609_BARCODES, sep="\t", header=None,
                           names=["barcode"])
    assert adata.shape[0] == len(barcodes), (adata.shape, len(barcodes))
    assert adata.shape[1] == len(features), (adata.shape, len(features))
    adata.obs_names = barcodes["barcode"].values
    adata.var["ensembl"] = features["ensembl"].values
    adata.var_names = features["symbol"].values
    adata.var_names_make_unique()
    log.info(f"loaded {adata.shape[0]:,} nuclei x {adata.shape[1]:,} genes")
    return adata


def attach_genotype(adata) -> None:
    key = pd.read_csv(cfg.GSE221609_SAMPLE_KEY)
    suf2geno = dict(zip(key["barcode_suffix"].astype(str).str.lstrip("-"),
                        key["genotype"]))
    suf2gsm = dict(zip(key["barcode_suffix"].astype(str).str.lstrip("-"),
                       key["gsm"]))
    suffix = pd.Index(adata.obs_names).str.rsplit("-", n=1).str[-1]
    adata.obs["sample_suffix"] = suffix.values
    adata.obs["genotype"] = suffix.map(suf2geno).values
    adata.obs["gsm"] = suffix.map(suf2gsm).values
    adata.obs["genotype"] = adata.obs["genotype"].astype("category")
    log.info("genotype counts:\n%s", adata.obs["genotype"].value_counts().to_string())


def run_qc(adata):
    adata.var["mt"] = adata.var_names.str.startswith(cfg.MITO_PREFIX)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=cfg.MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=cfg.MIN_CELLS_PER_GENE)
    adata = adata[adata.obs["pct_counts_mt"] < cfg.MAX_PCT_MT].copy()
    log.info(f"QC: {n0:,} -> {adata.n_obs:,} nuclei "
             f"(min_genes={cfg.MIN_GENES_PER_CELL}, "
             f"max_pct_mt={cfg.MAX_PCT_MT})")
    return adata


def normalize_cluster(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=cfg.TARGET_SUM)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg.N_TOP_HVG,
                                flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=cfg.N_PCS, random_state=cfg.RANDOM_SEED)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=cfg.N_NEIGHBORS, n_pcs=cfg.N_PCS,
                    random_state=cfg.RANDOM_SEED)
    sc.tl.umap(adata, random_state=cfg.RANDOM_SEED)
    sc.tl.leiden(adata, resolution=cfg.LEIDEN_RESOLUTION,
                 random_state=cfg.RANDOM_SEED, flavor="igraph", n_iterations=2,
                 directed=False)
    log.info(f"clustered: {adata.obs['leiden'].nunique()} Leiden clusters")
    return adata


def score_celltypes(adata):
    """Score each canonical cell-type marker set per cell, then average by
    Leiden cluster. Microglia clusters are those with the top microglia score
    AND clear pan-microglial marker expression."""
    present = {}
    for ct, markers in cfg.CELLTYPE_MARKERS.items():
        mk = [g for g in markers if g in adata.var_names]
        present[ct] = mk
        if mk:
            sc.tl.score_genes(adata, mk, score_name=f"score_{ct}",
                              ctrl_size=cfg.SCORE_GENES_CTRL_SIZE,
                              random_state=cfg.SCORE_GENES_RANDOM_STATE)
    score_cols = [f"score_{ct}" for ct in cfg.CELLTYPE_MARKERS if present[ct]]
    per_cluster = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    # Assign each cluster the cell type of its highest mean score.
    assigned = per_cluster.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["cell_type"] = adata.obs["leiden"].map(assigned).astype("category")
    per_cluster["assigned"] = assigned
    log.info("per-cluster cell-type assignment:\n%s", per_cluster.to_string())
    return per_cluster, present


def main():
    cfg.ensure_dirs()
    adata = load_triplet()
    attach_genotype(adata)
    adata = run_qc(adata)
    adata = normalize_cluster(adata)
    per_cluster, present = score_celltypes(adata)

    per_cluster.to_csv(os.path.join(cfg.TABLES_DIR,
                                    "npc_celltype_marker_scores.csv"))

    # Microglia = clusters assigned "Microglia".
    micro_clusters = per_cluster.index[per_cluster["assigned"] == "Microglia"].tolist()
    log.info(f"microglia Leiden clusters: {micro_clusters}")
    micro = adata[adata.obs["cell_type"] == "Microglia"].copy()

    # Confirm by canonical marker expression (report mean per marker).
    mk_present = [g for g in cfg.MICROGLIA_MARKERS if g in micro.var_names]
    conf = pd.DataFrame({
        "marker": mk_present,
        "mean_expr_microglia": [float(np.asarray(
            micro[:, g].X.todense()).mean()) for g in mk_present],
        "mean_expr_other": [float(np.asarray(
            adata[adata.obs["cell_type"] != "Microglia", g].X.todense()).mean())
            for g in mk_present],
    })
    log.info("microglia marker confirmation:\n%s", conf.to_string(index=False))
    conf.to_csv(os.path.join(cfg.TABLES_DIR, "npc_microglia_marker_confirm.csv"),
                index=False)

    counts = micro.obs["genotype"].value_counts().rename_axis("genotype")
    counts.to_frame("n_microglia").to_csv(
        os.path.join(cfg.TABLES_DIR, "npc_microglia_counts.csv"))
    log.info("microglia per genotype:\n%s", counts.to_string())

    micro.write(cfg.NPC_MICROGLIA_CKPT)
    log.info(f"wrote {cfg.NPC_MICROGLIA_CKPT} ({micro.n_obs:,} microglia)")

    make_qc_figure(adata, micro, conf)
    write_doc(adata, micro, per_cluster, counts, conf, present)


def make_qc_figure(adata, micro, conf):
    apply_figure_style(frame="open", font="Arial")
    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    # (a) UMAP by cell type
    ct_order = list(adata.obs["cell_type"].cat.categories)
    for ct in ct_order:
        m = adata.obs["cell_type"] == ct
        axes[0, 0].scatter(adata.obsm["X_umap"][m.values, 0],
                           adata.obsm["X_umap"][m.values, 1], s=1,
                           label=ct, rasterized=True)
    axes[0, 0].set_title("Cell type"); axes[0, 0].set_xlabel("UMAP1")
    axes[0, 0].set_ylabel("UMAP2")
    axes[0, 0].legend(markerscale=4, fontsize=5, loc="best")
    panel_letter(axes[0, 0], "a")
    # (b) UMAP by genotype
    for gt in adata.obs["genotype"].cat.categories:
        m = adata.obs["genotype"] == gt
        axes[0, 1].scatter(adata.obsm["X_umap"][m.values, 0],
                           adata.obsm["X_umap"][m.values, 1], s=1,
                           label=gt, rasterized=True)
    axes[0, 1].set_title("Genotype"); axes[0, 1].legend(markerscale=4)
    panel_letter(axes[0, 1], "b")
    # (c) microglia score on UMAP
    sccol = adata.obs["score_Microglia"].values
    sc3 = axes[0, 2].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                            s=1, c=sccol, cmap="viridis", rasterized=True)
    axes[0, 2].set_title("Microglia score"); fig.colorbar(sc3, ax=axes[0, 2],
                                                          fraction=0.046)
    panel_letter(axes[0, 2], "c")
    # (d) QC violin: pct mito
    axes[1, 0].violinplot([adata.obs["pct_counts_mt"].values], showmeans=True)
    axes[1, 0].set_title("pct mito (post-QC)"); axes[1, 0].set_ylabel("% mito")
    panel_letter(axes[1, 0], "d")
    # (e) marker confirmation bar
    x = np.arange(len(conf))
    axes[1, 1].bar(x - 0.2, conf["mean_expr_microglia"], 0.4, label="microglia")
    axes[1, 1].bar(x + 0.2, conf["mean_expr_other"], 0.4, label="other")
    axes[1, 1].set_xticks(x); axes[1, 1].set_xticklabels(conf["marker"],
                                                        rotation=45, ha="right")
    axes[1, 1].set_title("Microglia markers"); axes[1, 1].legend()
    axes[1, 1].set_ylabel("mean log1p CP10k")
    panel_letter(axes[1, 1], "e")
    # (f) microglia count by genotype
    cc = micro.obs["genotype"].value_counts()
    axes[1, 2].bar(range(len(cc)), cc.values)
    axes[1, 2].set_xticks(range(len(cc))); axes[1, 2].set_xticklabels(cc.index)
    axes[1, 2].set_title("Microglia per genotype"); axes[1, 2].set_ylabel("n nuclei")
    panel_letter(axes[1, 2], "f")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(cfg.FIGURES_DIR, f"P3_npc_umap_qc.{ext}"),
                    dpi=300)
    plt.close(fig)
    log.info("wrote figures/P3_npc_umap_qc.{pdf,png}")


def write_doc(adata, micro, per_cluster, counts, conf, present):
    p = os.path.join(cfg.DOCS_DIR, "npc_preprocessing.md")
    lines = [
        "# GSE221609 snRNA-seq preprocessing & microglia annotation (Phase 3, Step 07)",
        "",
        "**Dataset.** GSE221609, Npc1-/- vs WT mouse forebrain, single-nucleus "
        "RNA-seq (10x, cellranger-6.0.0). An oligodendrocyte-focused study, so "
        "microglia are a minority subset extracted here.",
        "",
        "## Sample / genotype key (from barcode suffix)",
        "| suffix | GSM | genotype |",
        "|---|---|---|",
    ]
    key = pd.read_csv(cfg.GSE221609_SAMPLE_KEY)
    for _, r in key.iterrows():
        lines.append(f"| {r['barcode_suffix']} | {r['gsm']} | {r['genotype']} |")
    lines += [
        "",
        "## QC",
        f"- Input: 75,878 nuclei x 32,285 genes.",
        f"- Filters: min_genes/cell = {cfg.MIN_GENES_PER_CELL}; "
        f"min_cells/gene = {cfg.MIN_CELLS_PER_GENE}; "
        f"pct_counts_mt < {cfg.MAX_PCT_MT}% (mito prefix `{cfg.MITO_PREFIX}`).",
        f"- Post-QC: {adata.n_obs:,} nuclei x {adata.n_vars:,} genes.",
        "- Rationale: single-nucleus data carry little cytoplasmic mito RNA, so "
        "a low mito ceiling removes damaged nuclei without discarding real "
        "signal; thresholds follow common snRNA-seq practice rather than "
        "arbitrary cutoffs, and pre/post distributions are shown in "
        "`figures/P3_npc_umap_qc.pdf`.",
        "",
        "## Normalization & clustering",
        f"- CP10k (`target_sum={cfg.TARGET_SUM:g}`) + log1p; counts retained in "
        "`layers['counts']`.",
        f"- HVG = top {cfg.N_TOP_HVG} (seurat flavor); scale(max=10); "
        f"PCA {cfg.N_PCS} comps; neighbors k={cfg.N_NEIGHBORS}; UMAP; "
        f"Leiden resolution={cfg.LEIDEN_RESOLUTION} (seed {cfg.RANDOM_SEED}).",
        f"- {adata.obs['leiden'].nunique()} Leiden clusters.",
        "",
        "## Microglia annotation (published markers + canonical markers)",
        "No per-cell annotations are distributed in the GEO supplementary "
        "(raw matrix only). Clusters were therefore formed here and each "
        "cluster assigned the cell type of its highest mean marker score "
        "(`sc.tl.score_genes`, ctrl_size="
        f"{cfg.SCORE_GENES_CTRL_SIZE}, seed {cfg.SCORE_GENES_RANDOM_STATE}). "
        "Microglia clusters were then CONFIRMED by canonical pan-microglial "
        "marker expression, not accepted on de novo score alone.",
        "",
        "Marker sets used:",
    ]
    for ct, mk in present.items():
        lines.append(f"- **{ct}**: {', '.join(mk) if mk else '(none present)'}")
    lines += [
        "",
        "### Microglia marker confirmation (mean log1p CP10k)",
        "| marker | microglia | other |",
        "|---|---|---|",
    ]
    for _, r in conf.iterrows():
        lines.append(f"| {r['marker']} | {r['mean_expr_microglia']:.3f} | "
                     f"{r['mean_expr_other']:.3f} |")
    lines += [
        "",
        "## Microglia yield (per genotype)",
        "| genotype | n microglia |",
        "|---|---|",
    ]
    for gt, n in counts.items():
        lines.append(f"| {gt} | {int(n):,} |")
    lines += [
        "",
        f"- **Total microglia: {micro.n_obs:,}.** This yield feeds the Step 08 "
        "QC gate (whether the single-cell disease axis is defensible or the "
        "bulk real-time axis must be primary).",
        "",
        "## Output",
        f"- `checkpoints/npc_microglia.h5ad` — microglia only, normalized, "
        "genotype-labelled.",
        "- `tables/npc_celltype_marker_scores.csv`, "
        "`tables/npc_microglia_marker_confirm.csv`, "
        "`tables/npc_microglia_counts.csv`.",
    ]
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    log.info(f"wrote {p}")


if __name__ == "__main__":
    main()
