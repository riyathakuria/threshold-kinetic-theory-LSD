import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _env import setup_numba_cache          # noqa: E402
setup_numba_cache()

import numpy as np                            # noqa: E402
import pandas as pd                           # noqa: E402
import scanpy as sc                           # noqa: E402
import anndata as anndata_mod                 # noqa: E402

import _phase3_config as cfg                  # noqa: E402

log = cfg.get_logger("09_project_modules")
sc.settings.verbosity = 1
np.random.seed(cfg.RANDOM_SEED)


def load_modules():
    """Return {cluster_id: {'name':..., 'mouse_genes':[...], 'human_genes':[...],
    'n_human':int, 'n_mapped':int}} from the ortholog map."""
    omap = pd.read_csv(cfg.ORTHOLOG_MAP)
    modules = {}
    for cid, grp in omap.groupby(omap["string_cluster_id"].astype(int)):
        mouse = grp.loc[grp["confidence"] == "1:1", "mouse_symbol"].tolist()
        modules[cid] = {
            "name": cfg.STRING_CLUSTER_NAMES.get(cid, grp["string_cluster_name"].iloc[0]),
            "mouse_genes": mouse,
            "human_genes": grp["human_symbol"].tolist(),
            "n_human": len(grp),
            "n_mapped": int((grp["confidence"] == "1:1").sum()),
        }
    return modules


def score_bulk(modules):
    expr = pd.read_csv(os.path.join(cfg.TABLES_DIR,
                                    "bulk_expression_normalized.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(cfg.TABLES_DIR, "bulk_sample_metadata.csv"))
    # Build AnnData samples x genes so score_genes parity holds.
    ad = anndata_mod.AnnData(X=expr.T.values.astype(np.float32))
    ad.obs_names = expr.columns
    ad.var_names = expr.index
    rep = []
    for cid, m in modules.items():
        genes = [g for g in m["mouse_genes"] if g in ad.var_names]
        present_frac = len(genes) / m["n_human"] if m["n_human"] else 0
        rep.append({"string_cluster_id": cid, "module_name": m["name"],
                    "axis": "bulk", "n_module_genes": m["n_human"],
                    "n_detected": len(genes),
                    "detected_frac": round(present_frac, 3),
                    "flag": "OK" if present_frac >= 0.5 else "LOW_DETECTION"})
        if genes:
            sc.tl.score_genes(ad, genes, score_name=f"M{cid}",
                              ctrl_size=cfg.SCORE_GENES_CTRL_SIZE,
                              random_state=cfg.SCORE_GENES_RANDOM_STATE)
        else:
            ad.obs[f"M{cid}"] = np.nan
    score_cols = [f"M{cid}" for cid in modules]
    out = ad.obs[score_cols].copy()
    out.index.name = "sample"
    out = out.merge(meta.set_index("sample"), left_index=True, right_index=True)
    out.to_csv(os.path.join(cfg.TABLES_DIR, "module_bulk_timecourse.csv"))
    log.info("wrote tables/module_bulk_timecourse.csv (%d samples x %d modules)",
             out.shape[0], len(score_cols))
    return out, rep


def score_pseudotime(modules):
    ckpt = os.path.join(cfg.CHECKPOINT_DIR, "npc_microglia_dpt.h5ad")
    if not os.path.exists(ckpt):
        log.warning("no DPT checkpoint (%s); skipping supporting axis", ckpt)
        return None, []
    ad = sc.read_h5ad(ckpt)
    rep = []
    for cid, m in modules.items():
        genes = [g for g in m["mouse_genes"] if g in ad.var_names]
        present_frac = len(genes) / m["n_human"] if m["n_human"] else 0
        rep.append({"string_cluster_id": cid, "module_name": m["name"],
                    "axis": "pseudotime", "n_module_genes": m["n_human"],
                    "n_detected": len(genes),
                    "detected_frac": round(present_frac, 3),
                    "flag": "OK" if present_frac >= 0.5 else "LOW_DETECTION"})
        if genes:
            sc.tl.score_genes(ad, genes, score_name=f"M{cid}",
                              ctrl_size=cfg.SCORE_GENES_CTRL_SIZE,
                              random_state=cfg.SCORE_GENES_RANDOM_STATE)
        else:
            ad.obs[f"M{cid}"] = np.nan
    score_cols = [f"M{cid}" for cid in modules]
    out = ad.obs[["genotype", "dpt_pseudotime"] + score_cols].copy()
    out.index.name = "cell"
    out.to_csv(os.path.join(cfg.TABLES_DIR, "module_pseudotime_scores.csv"))
    log.info("wrote tables/module_pseudotime_scores.csv (%d cells x %d modules)",
             out.shape[0], len(score_cols))
    return out, rep


def main():
    cfg.ensure_dirs()
    modules = load_modules()
    log.info("modules: %s", {cid: m["name"] for cid, m in modules.items()})

    bulk_out, bulk_rep = score_bulk(modules)
    pt_out, pt_rep = score_pseudotime(modules)

    rep = pd.DataFrame(bulk_rep + pt_rep)
    rep.to_csv(os.path.join(cfg.TABLES_DIR, "module_mapping_report.csv"),
               index=False)
    low = rep[rep["flag"] != "OK"]
    if len(low):
        log.warning("modules flagged LOW_DETECTION:\n%s", low.to_string(index=False))
    else:
        log.info("all modules adequately detected on both axes")

    # Quick per-week module means on the primary axis, for a glance.
    wk = (bulk_out.groupby(["genotype", "week"])[[f"M{c}" for c in modules]]
          .mean())
    log.info("module means by genotype x week (primary bulk axis):\n%s",
             wk.round(3).to_string())


if __name__ == "__main__":
    main()
