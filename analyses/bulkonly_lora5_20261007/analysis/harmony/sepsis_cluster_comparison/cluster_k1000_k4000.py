from __future__ import annotations

import gc
import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[4]
RUNS = [
    ("K1000", ROOT / "bulkonly_lora5_20261006/analysis/harmony/sepsis_k1000/sepsis_k1000_harmony_embedding.h5ad", ROOT / "bulkonly_lora5_20261006/analysis/harmony/sepsis_k1000"),
    ("K4000", ROOT / "bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000/sepsis_k4000_harmony_embedding.h5ad", ROOT / "bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000"),
]
OUT = Path(__file__).parent
SEED = 20261007


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    for name, input_path, output_dir in RUNS:
        print(f"Reading {name}: {input_path}", flush=True)
        a = ad.read_h5ad(input_path)
        if "X_pca_harmony" not in a.obsm:
            raise KeyError(f"{name}: X_pca_harmony is absent")
        if "celltype_l1_pred" not in a.obs:
            raise KeyError(f"{name}: celltype_l1_pred is absent")
        print(f"{name}: {a.n_obs:,} cells, {a.n_vars:,} genes; obs={list(a.obs.columns)}", flush=True)
        a.obs["celltype_l1_pred"] = a.obs["celltype_l1_pred"].astype("category")
        a.uns["paper_comparison"] = {
            "paper": "Ye et al. Nature Immunology 2026, DOI 10.1038/s41590-025-02345-x",
            "input_scope": "bulk-only model-predicted cells, not measured scRNA-seq observations",
            "clustering_space": "first 10 Harmony-corrected PCs",
            "neighbors": 30,
            "leiden_resolutions": [0.5, 1.0],
            "note": "Paper uses fastMNN and top 10 PCs for graph clustering/UMAP; Harmony replaces fastMNN. Paper's exact clustering resolution is not reported in accessible Methods text.",
        }
        graph_key = "neighbors_harmony10"
        sc.pp.neighbors(a, use_rep="X_pca_harmony", n_pcs=10, n_neighbors=30, random_state=SEED, key_added=graph_key)
        for resolution in (0.5, 1.0):
            key = f"leiden_harmony10_r{resolution:g}"
            sc.tl.leiden(a, resolution=resolution, neighbors_key=graph_key, key_added=key, random_state=SEED)
        sc.tl.umap(a, neighbors_key=graph_key, random_state=SEED, min_dist=0.3)
        sc.settings.figdir = str(OUT / "figures")
        Path(sc.settings.figdir).mkdir(exist_ok=True)

        label = a.obs["celltype_l1_pred"].astype(str)
        records = {"K": name, "n_cells": int(a.n_obs), "n_genes": int(a.n_vars), "n_predicted_l1_types": int(label.nunique())}
        for resolution in (0.5, 1.0):
            key = f"leiden_harmony10_r{resolution:g}"
            cluster = a.obs[key].astype(str)
            tab = pd.crosstab(cluster, label, margins=False)
            tab.index.name = "leiden_cluster"
            tab.to_csv(OUT / f"{name.lower()}_r{resolution:g}_cluster_by_predicted_l1.csv")
            sizes = cluster.value_counts().sort_index().rename_axis("leiden_cluster").rename("n_cells").to_frame()
            sizes["dominant_predicted_l1"] = tab.idxmax(axis=1).reindex(sizes.index)
            sizes["dominant_type_fraction"] = (tab.max(axis=1) / tab.sum(axis=1)).reindex(sizes.index)
            sizes.to_csv(OUT / f"{name.lower()}_r{resolution:g}_cluster_summary.csv")
            records[f"clusters_r{resolution:g}"] = int(cluster.nunique())
            records[f"ARI_cluster_vs_predicted_l1_r{resolution:g}"] = float(adjusted_rand_score(label, cluster))
            records[f"NMI_cluster_vs_predicted_l1_r{resolution:g}"] = float(normalized_mutual_info_score(label, cluster))
        label.value_counts().rename_axis("predicted_l1").rename("n_cells").to_csv(OUT / f"{name.lower()}_predicted_l1_counts.csv")
        sc.pl.umap(a, color=["leiden_harmony10_r0.5", "celltype_l1_pred", "draw"], ncols=1, show=False, save=f"_{name.lower()}_leiden_r0.5.png")
        # Save a copy carrying the clustering labels; original Harmony integration files remain unchanged.
        out_h5ad = output_dir / f"sepsis_{name.lower()}_harmony_clustering.h5ad"
        a.write_h5ad(out_h5ad, compression="gzip")
        records["clustered_h5ad"] = str(out_h5ad)
        summary.append(records)
        print(json.dumps(records, ensure_ascii=False), flush=True)
        del a
        gc.collect()
    pd.DataFrame(summary).to_csv(OUT / "k1000_k4000_cluster_metrics.csv", index=False)
    (OUT / "k1000_k4000_cluster_metrics.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
