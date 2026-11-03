from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import harmonypy as hm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
INPUT = ROOT / "ye2026_sepsis"
N_PER_SAMPLE_DRAW = 500
N_HVG = 2000
N_PCS = 50
SEED = 20261007

matplotlib.use("Agg")
sc.settings.verbosity = 2


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pieces: list[ad.AnnData] = []
    input_records: list[dict[str, object]] = []
    sample_counts: dict[str, dict[str, int]] = {}

    for draw in range(5):
        path = INPUT / f"lora5_bulkonly_ye2026_sepsis_draw{draw}.h5ad"
        a = ad.read_h5ad(path)
        if "measured_in_cohort" not in a.var:
            raise KeyError(f"{path.name}: var['measured_in_cohort'] is missing")
        measured = a.var["measured_in_cohort"].astype(bool).to_numpy()
        a = a[:, measured].copy()
        if "sample_id" not in a.obs:
            raise KeyError(f"{path.name}: obs['sample_id'] is missing")
        a.obs["draw"] = f"draw{draw}"
        a.obs["K"] = 1000
        a.obs["cohort"] = "ye2026_sepsis"

        rng = np.random.default_rng(SEED + draw)
        chosen: list[np.ndarray] = []
        sample_counts[f"draw{draw}"] = {}
        for sample_id, positions in a.obs.groupby("sample_id", observed=True, sort=False).indices.items():
            positions = np.asarray(positions, dtype=np.int64)
            n_take = min(N_PER_SAMPLE_DRAW, len(positions))
            chosen.append(rng.choice(positions, size=n_take, replace=False))
            sample_counts[f"draw{draw}"][str(sample_id)] = int(len(positions))
        selected = np.sort(np.concatenate(chosen))
        a = a[selected, :].copy()
        pieces.append(a)
        input_records.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "n_obs_full": int(sum(sample_counts[f"draw{draw}"].values())),
                "n_obs_sampled": int(a.n_obs),
                "n_vars_measured": int(a.n_vars),
                "n_samples": int(a.obs["sample_id"].nunique()),
            }
        )
        print(f"Loaded draw{draw}: {a.n_obs:,} sampled cells × {a.n_vars:,} measured genes", flush=True)

    adata = ad.concat(
        pieces,
        join="inner",
        merge="same",
        index_unique="-",
        pairwise=False,
    )
    adata.obs["draw"] = pd.Categorical(adata.obs["draw"], categories=[f"draw{i}" for i in range(5)], ordered=True)
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    print(f"Combined object: {adata.n_obs:,} cells × {adata.n_vars:,} measured genes", flush=True)
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(N_HVG, adata.n_vars),
        flavor="seurat",
        batch_key="draw",
        subset=False,
    )

    hvg = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    print(f"Selected {hvg.n_vars:,} highly variable genes; scaling and PCA", flush=True)
    sc.pp.scale(hvg, max_value=10)
    sc.pp.pca(hvg, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"].copy()
    del hvg

    print("Running Harmony with draw as the batch covariate", flush=True)
    harmony = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        vars_use=["draw"],
        max_iter_harmony=20,
        random_state=SEED,
        verbose=True,
    )
    adata.obsm["X_pca_harmony"] = np.asarray(harmony.Z_corr).T.astype(np.float32, copy=False)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED, min_dist=0.3)

    adata.uns["harmony_run"] = {
        "cohort": "ye2026_sepsis",
        "K_recipe_size": 1000,
        "draws": [f"draw{i}" for i in range(5)],
        "batch_covariate": "draw",
        "measured_gene_filter": "var['measured_in_cohort'] == True",
        "cells_per_sample_per_draw_for_embedding": N_PER_SAMPLE_DRAW,
        "n_hvg": int(adata.var["highly_variable"].sum()),
        "n_pcs": N_PCS,
        "harmony_version": "0.0.10",
        "scanpy_version": sc.__version__,
        "anndata_version": ad.__version__,
        "random_seed": SEED,
        "counts_preserved_in_layer": "counts",
        "expression_matrix_X": "library-size normalized to 10,000 and log1p transformed",
        "note": "Predicted cells, not measured single-cell observations; embedding subset is stratified by sample_id and draw.",
    }

    out_h5ad = OUT / "sepsis_k1000_harmony_embedding.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"Saved {out_h5ad}", flush=True)

    draw_order = [f"draw{i}" for i in range(5)]
    draw_colors = {
        "draw0": "#1F6F5F",
        "draw1": "#242424",
        "draw2": "#555555",
        "draw3": "#858585",
        "draw4": "#B6B6B0",
    }
    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(7.0, 5.1), constrained_layout=True)
    ax.set_facecolor("white")
    for draw in draw_order:
        mask = adata.obs["draw"].astype(str).to_numpy() == draw
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=1.0,
            alpha=0.42,
            c=draw_colors[draw],
            linewidths=0,
            rasterized=True,
        )
    ax.set_title("Sepsis K=1,000 predicted cells in Harmony space", loc="left", fontsize=11, weight="semibold", color="#171717", pad=10)
    ax.set_xlabel("UMAP 1", fontsize=9, color="#404040")
    ax.set_ylabel("UMAP 2", fontsize=9, color="#404040")
    ax.tick_params(labelsize=8, colors="#5E6462", length=2.5, width=0.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#B6B6B0")
    ax.spines["bottom"].set_color("#B6B6B0")
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=draw_colors[d], markeredgewidth=0, label=d)
        for d in draw_order
    ]
    ax.legend(handles=handles, title="Model draw", frameon=False, fontsize=8, title_fontsize=8, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.savefig(OUT / "umap_by_draw.png", dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(OUT / "umap_by_draw.svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "cohort": "ye2026_sepsis",
        "K_recipe_size": 1000,
        "files": input_records,
        "n_cells_integrated": int(adata.n_obs),
        "n_measured_genes": int(adata.n_vars),
        "n_hvg": int(adata.var["highly_variable"].sum()),
        "n_pcs": N_PCS,
        "cells_per_sample_draw": N_PER_SAMPLE_DRAW,
        "sample_cell_counts_per_draw_min": min(v for draw_counts in sample_counts.values() for v in draw_counts.values()) if sample_counts else None,
        "sample_cell_counts_per_draw_max": max(v for draw_counts in sample_counts.values() for v in draw_counts.values()) if sample_counts else None,
        "celltype_l1_counts": adata.obs["celltype_l1_pred"].astype(str).value_counts().to_dict() if "celltype_l1_pred" in adata.obs else {},
        "draw_cell_counts": adata.obs["draw"].astype(str).value_counts().sort_index().to_dict(),
        "output_h5ad": out_h5ad.name,
        "plots": ["umap_by_draw.png", "umap_by_draw.svg"],
        "interpretation": "Harmony is for embeddings/neighborhoods/clustering; use original sample-level data and paired calibration for inference.",
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
