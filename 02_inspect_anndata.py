#!/usr/bin/env python3
"""
Phase 1 - Step 4: Inspect the HuMicA v2 AnnData object (RAM-safe, backed mode).

Because total RAM (~8.6 GB) is smaller than the on-disk matrix, the object is
opened with backed='r' so the expression matrix is never materialized in RAM.
We read structure from .obs/.var/.raw/.layers/.uns/.obsm and, where matrix stats
are needed, we stream via h5py directly.

Writes:
  docs/humica_structure_report.md
  data/metadata/humica_var_index.csv        (.X var names)
  data/metadata/humica_raw_var_index.csv     (.raw var names, if present)
  data/metadata/humica_obs_head.csv          (first 200 obs rows)
  data/metadata/humica_structure.json        (machine-readable summary)
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import ROOT, setup_numba_cache
setup_numba_cache()

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

RAW = os.path.join(ROOT, "data", "raw", "HuMicA_v2.0.0.h5ad")
DOCS = os.path.join(ROOT, "docs")
META = os.path.join(ROOT, "data", "metadata")


def h5_tree(path, max_depth=2):
    """Top-level structure of the .h5ad without loading matrices."""
    out = {}
    with h5py.File(path, "r") as f:
        def walk(g, prefix, depth):
            for k in g.keys():
                item = g[k]
                key = f"{prefix}/{k}"
                if isinstance(item, h5py.Group):
                    out[key] = {"type": "group", "keys": list(item.keys())[:50]}
                    if depth < max_depth:
                        walk(item, key, depth + 1)
                else:
                    out[key] = {"type": "dataset", "shape": item.shape, "dtype": str(item.dtype)}
        walk(f, "", 0)
    return out


def main():
    os.makedirs(DOCS, exist_ok=True)
    os.makedirs(META, exist_ok=True)

    # --- raw h5 structure (no matrix load) ---
    tree = h5_tree(RAW)

    # --- backed-mode AnnData (metadata only) ---
    adata = sc.read_h5ad(RAW, backed="r")
    n_obs, n_var = adata.shape

    summary = {}
    summary["file"] = RAW
    summary["n_obs_cells"] = int(n_obs)
    summary["n_var_genes_X"] = int(n_var)
    summary["X_dtype"] = str(adata.X.dtype) if adata.X is not None else None

    # .var
    var = adata.var.copy()
    var.to_csv(os.path.join(META, "humica_var_index.csv"))
    summary["var_columns"] = list(var.columns)
    summary["var_index_name"] = var.index.name
    summary["var_index_sample"] = var.index[:10].tolist()

    # .obs
    obs = adata.obs
    summary["obs_columns"] = list(obs.columns)
    summary["n_obs_columns"] = len(obs.columns)
    obs.head(200).to_csv(os.path.join(META, "humica_obs_head.csv"))
    # categorical value counts for likely annotation cols
    obs_summaries = {}
    for c in obs.columns:
        try:
            nun = obs[c].nunique(dropna=True)
        except Exception:
            continue
        if nun <= 60:
            obs_summaries[c] = obs[c].value_counts(dropna=False).head(60).to_dict()
        else:
            obs_summaries[c] = {"_n_unique": int(nun), "_dtype": str(obs[c].dtype)}
    def _coerce(vv):
        try:
            return int(vv)
        except (ValueError, TypeError):
            return str(vv)
    summary["obs_value_counts"] = {k: {str(kk): _coerce(vv) for kk, vv in v.items()}
                                   for k, v in obs_summaries.items()}

    # .raw
    if adata.raw is not None:
        rv = adata.raw.var
        summary["has_raw"] = True
        summary["n_raw_genes"] = int(adata.raw.shape[1])
        summary["raw_var_columns"] = list(rv.columns)
        summary["raw_var_index_name"] = rv.index.name
        summary["raw_var_index_sample"] = rv.index[:10].tolist()
        pd.DataFrame(index=rv.index).to_csv(os.path.join(META, "humica_raw_var_index.csv"))
        rv.to_csv(os.path.join(META, "humica_raw_var_full.csv"))
    else:
        summary["has_raw"] = False
        summary["n_raw_genes"] = 0

    # .layers
    summary["layers"] = list(adata.layers.keys()) if adata.layers is not None else []

    # .obsm (embeddings)
    summary["obsm_keys"] = list(adata.obsm.keys()) if adata.obsm is not None else []
    summary["obsm_shapes"] = {k: list(np.asarray(adata.obsm[k]).shape) for k in summary["obsm_keys"]}

    # .uns
    try:
        summary["uns_keys"] = list(adata.uns.keys())
    except Exception:
        summary["uns_keys"] = []

    # --- heuristic: is .X HVG-only? is .raw the full transcriptome? ---
    summary["X_is_subset_of_raw"] = (
        summary.get("has_raw") and summary["n_raw_genes"] > summary["n_var_genes_X"]
    )

    with open(os.path.join(META, "humica_structure.json"), "w") as f:
        json.dump({"tree": {k: str(v) for k, v in tree.items()}, "summary": summary}, f, indent=1)

    # --- write markdown report ---
    write_report(summary, tree)
    print("Inspection complete.")
    print(f"  cells={summary['n_obs_cells']}  X_genes={summary['n_var_genes_X']}  "
          f"raw_genes={summary['n_raw_genes']}  layers={summary['layers']}")


def write_report(s, tree):
    L = []
    L.append("# HuMicA v2.0.0 - Structure Report\n")
    L.append(f"Source file: `{os.path.basename(s['file'])}`\n")
    L.append("## Dimensions\n")
    L.append(f"- **Cells (obs)**: {s['n_obs_cells']:,}")
    L.append(f"- **Genes in `.X` (var)**: {s['n_var_genes_X']:,}")
    L.append(f"- **Genes in `.raw`**: {s['n_raw_genes']:,}" + ("" if s['has_raw'] else "  (no .raw)"))
    L.append(f"- **`.X` dtype**: {s['X_dtype']}\n")

    L.append("## Interpretation\n")
    if s["has_raw"] and s["X_is_subset_of_raw"]:
        L.append(f"- `.raw` holds **{s['n_raw_genes']:,}** genes vs **{s['n_var_genes_X']:,}** in `.X`; "
                 "`.X` is a **subset (likely HVG / analysis matrix)** and **`.raw` carries the fuller "
                 "transcriptome**. Gene-presence checks MUST consult `.raw.var`, not only `.var`.")
    elif s["has_raw"]:
        L.append(f"- `.raw` present with {s['n_raw_genes']:,} genes (not larger than `.X`); "
                 "`.X` is not a strict HVG subset of `.raw`.")
    else:
        L.append("- No `.raw` present; `.X` var index is the only gene namespace.")
    L.append("")

    L.append("## `.var` (X gene annotation)\n")
    L.append(f"- index name: `{s['var_index_name']}`; columns: {s['var_columns']}")
    L.append(f"- sample gene IDs: {s['var_index_sample']}\n")
    if s["has_raw"]:
        L.append("## `.raw.var` (raw gene annotation)\n")
        L.append(f"- index name: `{s['raw_var_index_name']}`; columns: {s['raw_var_columns']}")
        L.append(f"- sample gene IDs: {s['raw_var_index_sample']}\n")

    L.append("## `.obs` (cell metadata)\n")
    L.append(f"- {s['n_obs_columns']} columns: {s['obs_columns']}\n")
    L.append("### Annotation column value counts (<=60 categories)\n")
    for c, vc in s["obs_value_counts"].items():
        if "_n_unique" in vc:
            L.append(f"- **{c}**: {vc['_n_unique']} unique values (dtype {vc['_dtype']})")
        else:
            top = list(vc.items())[:12]
            L.append(f"- **{c}**: " + ", ".join(f"{k}={v}" for k, v in top)
                     + (" ..." if len(vc) > 12 else ""))
    L.append("")

    L.append("## Embeddings (`.obsm`)\n")
    if s["obsm_keys"]:
        for k in s["obsm_keys"]:
            L.append(f"- `{k}`: shape {s['obsm_shapes'][k]}")
    else:
        L.append("- none")
    L.append("")

    L.append("## Layers\n")
    L.append(f"- {s['layers'] if s['layers'] else 'none'}\n")

    L.append("## `.uns` keys\n")
    L.append(f"- {s['uns_keys'] if s['uns_keys'] else 'none'}\n")

    L.append("## Raw HDF5 top-level structure\n")
    L.append("```")
    for k in sorted(tree.keys()):
        if k.count("/") <= 2:
            L.append(f"{k}: {tree[k]}")
    L.append("```")

    with open(os.path.join(DOCS, "humica_structure_report.md"), "w") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()
