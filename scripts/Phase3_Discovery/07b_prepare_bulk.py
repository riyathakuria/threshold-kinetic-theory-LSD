#!/usr/bin/env python3
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                            # noqa: E402
import pandas as pd                           # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402

import _phase3_config as cfg                  # noqa: E402
from _figstyle import apply_figure_style, panel_letter  # noqa: E402

log = cfg.get_logger("07b_prepare_bulk")

SAMPLE_RE = re.compile(r"^(WT|NPC)-(\d+)([A-D])$")


def parse_samples(columns):
    """Return {sample: (genotype, week, rep)} for disease-course samples only."""
    keep = {}
    for c in columns:
        m = SAMPLE_RE.match(c)
        if not m:
            continue
        geno, week, rep = m.group(1), int(m.group(2)), m.group(3)
        if any(p in c for p in cfg.BULK_EXCLUDE_PATTERNS):
            continue
        if week not in cfg.BULK_WEEKS:
            continue
        keep[c] = (geno, week, rep)
    return keep


def main():
    cfg.ensure_dirs()
    log.info("reading normalized count matrix")
    df = pd.read_csv(cfg.GSE152158_CSV, index_col=0)
    log.info(f"raw matrix: {df.shape[0]:,} genes x {df.shape[1]} samples")

    keep = parse_samples(df.columns)
    excluded = [c for c in df.columns if c not in keep]
    log.info(f"disease-course samples kept: {len(keep)}; excluded "
             f"(PLX/Con/other): {excluded}")
    df = df[list(keep.keys())]

    # Drop genes with zero signal across all kept samples.
    nz = (df.sum(axis=1) > 0)
    df = df.loc[nz]
    log.info(f"after dropping all-zero genes: {df.shape[0]:,} genes")

    # log2(x+1) — matrix is already normalized counts, so no re-TMM.
    logdf = np.log2(df + 1.0)
    logdf.to_csv(os.path.join(cfg.TABLES_DIR, "bulk_expression_normalized.csv"))
    log.info("wrote tables/bulk_expression_normalized.csv (log2 normalized)")

    meta = pd.DataFrame(
        [{"sample": s, "genotype": g, "week": w, "replicate": r}
         for s, (g, w, r) in keep.items()]
    ).sort_values(["genotype", "week", "replicate"]).reset_index(drop=True)
    meta.to_csv(os.path.join(cfg.TABLES_DIR, "bulk_sample_metadata.csv"),
                index=False)
    log.info("sample counts by genotype x week:\n%s",
             meta.groupby(["genotype", "week"]).size().to_string())

    make_qc_figure(logdf, meta)
    write_doc(df, logdf, meta, excluded)


def make_qc_figure(logdf, meta):
    apply_figure_style(frame="open", font="Arial")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    order = meta["sample"].tolist()
    L = logdf[order]
    # (a) library size (sum of normalized counts) per sample
    libs = L.sum(axis=0)
    colors = ["#4C72B0" if g == "WT" else "#C44E52"
              for g in meta.set_index("sample").loc[order, "genotype"]]
    axes[0].bar(range(len(order)), libs.values, color=colors)
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels(order, rotation=90, fontsize=4)
    axes[0].set_title("Column sum (log2 norm)"); axes[0].set_ylabel("sum")
    panel_letter(axes[0], "a")
    # (b) PCA of samples
    from sklearn.decomposition import PCA
    top = L.loc[L.var(axis=1).sort_values(ascending=False).index[:2000]]
    pcs = PCA(n_components=2, random_state=0).fit_transform(top.T.values)
    for g, col in [("WT", "#4C72B0"), ("NPC", "#C44E52")]:
        m = meta.set_index("sample").loc[order, "genotype"].values == g
        axes[1].scatter(pcs[m, 0], pcs[m, 1], c=col, label=g, s=30)
        for i, s in enumerate(order):
            if m[i]:
                wk = meta.set_index("sample").loc[s, "week"]
                axes[1].annotate(str(wk), (pcs[i, 0], pcs[i, 1]), fontsize=4)
    axes[1].set_title("Sample PCA (top 2000 var)")
    axes[1].set_xlabel("PC1"); axes[1].set_ylabel("PC2"); axes[1].legend()
    panel_letter(axes[1], "b")
    # (c) sample-sample correlation heatmap
    corr = np.corrcoef(top.T.values)
    im = axes[2].imshow(corr, cmap="viridis", aspect="auto")
    axes[2].set_xticks(range(len(order))); axes[2].set_yticks(range(len(order)))
    axes[2].set_xticklabels(order, rotation=90, fontsize=3)
    axes[2].set_yticklabels(order, fontsize=3)
    axes[2].set_title("Sample correlation"); fig.colorbar(im, ax=axes[2],
                                                         fraction=0.046)
    panel_letter(axes[2], "c")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(cfg.FIGURES_DIR, f"P3_bulk_qc.{ext}"), dpi=300)
    plt.close(fig)
    log.info("wrote figures/P3_bulk_qc.{pdf,png}")


def write_doc(df, logdf, meta, excluded):
    p = os.path.join(cfg.DOCS_DIR, "npc_preprocessing.md")
    # Append a bulk section to the same preprocessing doc.
    section = [
        "",
        "---",
        "",
        "# GSE152158 bulk RNA-seq preprocessing (Phase 3, Step 07b)",
        "",
        "**Dataset.** GSE152158, bulk RNA-seq of FACS-purified mouse microglia "
        "across the NPC disease course. The GEO file "
        "`GSE152158_NormCountData.csv.gz` is ALREADY a normalized count matrix, "
        "so no TMM/CPM was re-applied (that would double-normalize); a "
        "log2(x+1) transform was applied for downstream module scoring.",
        "",
        "## Disease-course axis (kept samples)",
        "- Genotypes: WT and NPC (Npc1-/-).",
        f"- Weeks: {cfg.BULK_WEEKS}. Replicates A-D.",
        f"- Kept {meta.shape[0]} samples; EXCLUDED the CSF1R-inhibitor arms "
        f"(pattern {cfg.BULK_EXCLUDE_PATTERNS}): {excluded}.",
        "",
        "### Samples by genotype x week",
        "| genotype | week | n |",
        "|---|---|---|",
    ]
    g = meta.groupby(["genotype", "week"]).size()
    for (geno, week), n in g.items():
        section.append(f"| {geno} | {week} | {n} |")
    section += [
        "",
        f"- Matrix after dropping all-zero genes: {df.shape[0]:,} genes x "
        f"{df.shape[1]} samples.",
        "",
        "## Output",
        "- `tables/bulk_expression_normalized.csv` — log2 normalized, genes x "
        "disease-course samples.",
        "- `tables/bulk_sample_metadata.csv` — sample | genotype | week | replicate.",
        "- `figures/P3_bulk_qc.pdf` — column sums, sample PCA, sample correlation.",
    ]
    with open(p, "a") as fh:
        fh.write("\n".join(section) + "\n")
    log.info(f"appended bulk section to {p}")


if __name__ == "__main__":
    main()
