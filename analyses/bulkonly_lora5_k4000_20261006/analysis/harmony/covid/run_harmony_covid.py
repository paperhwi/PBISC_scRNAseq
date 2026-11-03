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
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chisquare

WORKSPACE = Path(__file__).resolve().parents[4]
ROOTS = {1000: WORKSPACE / "bulkonly_lora5_20261006",
         4000: WORKSPACE / "bulkonly_lora5_k4000_20261006"}
N_PER_SAMPLE_DRAW = 500
N_HVG = 2000
N_PCS = 50
SEED = 20261007
K_NEIGHBORS = 30
TEST_FRACTION = 0.10
REPEATS = 100
NULL_PERMUTATIONS = 20


def kbet_rejection(distances, batch_codes, rng, repeats=REPEATS, alpha=0.05):
    """kBET Pearson chi-square rejection-rate estimate on a fixed kNN graph.

    Repeats random 10% query subsets; null baseline shuffles batch labels while
    preserving global batch sizes. Returns per-neighborhood reject fraction.
    """
    degrees = np.diff(distances.indptr)
    unique_degrees = np.unique(degrees)
    if len(unique_degrees) != 1:
        raise ValueError(f"Expected a fixed neighborhood size; found {unique_degrees}")
    k_actual = int(unique_degrees[0])
    neighbors = distances.indices.reshape(-1, k_actual)
    n_batches = int(batch_codes.max()) + 1
    expected = np.full(n_batches, k_actual / n_batches, dtype=float)
    n_test = max(25, int(round(len(batch_codes) * TEST_FRACTION)))
    observed_rates = []
    null_rates = []

    for _ in range(repeats):
        query = rng.choice(len(batch_codes), size=n_test, replace=False)
        local = batch_codes[neighbors[query]]
        counts = np.stack([(local == b).sum(axis=1) for b in range(n_batches)], axis=1)
        observed_rates.append(float(np.mean(chisquare(counts, f_exp=expected, axis=1).pvalue < alpha)))

    query = rng.choice(len(batch_codes), size=n_test, replace=False)
    for _ in range(NULL_PERMUTATIONS):
        permuted = rng.permutation(batch_codes)
        local = permuted[neighbors[query]]
        counts = np.stack([(local == b).sum(axis=1) for b in range(n_batches)], axis=1)
        null_rates.append(float(np.mean(chisquare(counts, f_exp=expected, axis=1).pvalue < alpha)))

    return {
        "observed_rejection_rate": float(np.mean(observed_rates)),
        "observed_rejection_rate_sd_across_subsamples": float(np.std(observed_rates, ddof=1)),
        "null_expected_rejection_rate": float(np.mean(null_rates)),
        "null_rejection_rate_sd": float(np.std(null_rates, ddof=1)),
        "alpha": alpha,
        "n_neighbors_requested": K_NEIGHBORS,
        "n_neighbors_actual": k_actual,
        "n_test_per_repeat": n_test,
        "test_fraction": TEST_FRACTION,
        "subsample_repeats": repeats,
        "null_permutations": NULL_PERMUTATIONS,
        "method": "Pearson chi-square goodness-of-fit of local draw counts vs global proportions; report fraction p<alpha",
    }


def run_one(k: int) -> None:
    root = ROOTS[k]
    input_dir = root / "imperial_eskd"
    out = root / "analysis" / "harmony" / f"covid_k{k}"
    out.mkdir(parents=True, exist_ok=True)
    pieces = []
    input_records = []
    sample_counts = {}

    for draw in range(5):
        prefix = "lora5k4_" if k == 4000 else "lora5_"
        path = input_dir / f"{prefix}bulkonly_imperial_eskd_draw{draw}.h5ad"
        a = ad.read_h5ad(path)
        measured = a.var["measured_in_cohort"].astype(bool).to_numpy()
        a = a[:, measured].copy()
        a.obs["draw"] = f"draw{draw}"
        a.obs["K"] = k
        a.obs["cohort"] = "imperial_eskd_covid"
        rng = np.random.default_rng(SEED + draw)
        chosen = []
        sample_counts[f"draw{draw}"] = {}
        for sample_id, positions in a.obs.groupby("sample_id", observed=True, sort=False).indices.items():
            positions = np.asarray(positions, dtype=np.int64)
            take = min(N_PER_SAMPLE_DRAW, len(positions))
            chosen.append(rng.choice(positions, size=take, replace=False))
            sample_counts[f"draw{draw}"][str(sample_id)] = int(len(positions))
        a = a[np.sort(np.concatenate(chosen)), :].copy()
        pieces.append(a)
        input_records.append({"file": path.name, "bytes": path.stat().st_size,
                              "n_obs_full": int(sum(sample_counts[f"draw{draw}"].values())),
                              "n_obs_sampled": int(a.n_obs), "n_vars_measured": int(a.n_vars),
                              "n_samples": int(a.obs["sample_id"].nunique())})
        print(f"K={k} draw{draw}: {a.n_obs:,} cells × {a.n_vars:,} measured genes", flush=True)

    adata = ad.concat(pieces, join="inner", merge="same", index_unique="-", pairwise=False)
    adata.obs["draw"] = pd.Categorical(adata.obs["draw"], categories=[f"draw{i}" for i in range(5)], ordered=True)
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(N_HVG, adata.n_vars), flavor="seurat", batch_key="draw", subset=False)
    hvg = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.pp.pca(hvg, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"].copy()
    del hvg

    # Compute pre-Harmony neighborhood graph for a transparent kBET baseline.
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=K_NEIGHBORS, random_state=SEED)
    pre_distances = adata.obsp["distances"].copy().tocsr()

    harmony = hm.run_harmony(adata.obsm["X_pca"], adata.obs, vars_use=["draw"],
                             max_iter_harmony=20, random_state=SEED, verbose=True)
    adata.obsm["X_pca_harmony"] = np.asarray(harmony.Z_corr).T.astype(np.float32, copy=False)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=K_NEIGHBORS, random_state=SEED)
    post_distances = adata.obsp["distances"].copy().tocsr()

    codes = adata.obs["draw"].cat.codes.to_numpy()
    kbet = {
        "pre_harmony_pca": kbet_rejection(pre_distances, codes, np.random.default_rng(SEED + 100)),
        "post_harmony": kbet_rejection(post_distances, codes, np.random.default_rng(SEED + 200)),
    }
    adata.uns["harmony_run"] = {
        "cohort": "imperial_eskd_covid", "K_recipe_size": k,
        "draws": [f"draw{i}" for i in range(5)], "batch_covariate": "draw",
        "measured_gene_filter": "var['measured_in_cohort'] == True",
        "cells_per_sample_per_draw_for_embedding": N_PER_SAMPLE_DRAW,
        "n_hvg": int(adata.var["highly_variable"].sum()), "n_pcs": N_PCS,
        "harmony_version": "0.0.10", "scanpy_version": sc.__version__,
        "anndata_version": ad.__version__, "random_seed": SEED,
        "counts_preserved_in_layer": "counts",
        "expression_matrix_X": "library-size normalized to 10,000 and log1p transformed",
        "kBET": kbet,
        "note": "Predicted cells, not measured single-cell observations; sample × draw stratified subset.",
    }
    out_h5ad = out / f"covid_k{k}_harmony_embedding.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")
    summary = {
        "cohort": "imperial_eskd_covid", "K_recipe_size": k, "files": input_records,
        "n_cells_integrated": int(adata.n_obs), "n_measured_genes": int(adata.n_vars),
        "n_hvg": int(adata.var["highly_variable"].sum()), "n_pcs": N_PCS,
        "cells_per_sample_draw": N_PER_SAMPLE_DRAW,
        "source_cells_per_sample_min": min(v for d in sample_counts.values() for v in d.values()),
        "source_cells_per_sample_max": max(v for d in sample_counts.values() for v in d.values()),
        "draw_cell_counts": adata.obs["draw"].astype(str).value_counts().sort_index().to_dict(),
        "kBET": kbet, "output_h5ad": out_h5ad.name,
    }
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    for setting in (1000, 4000):
        run_one(setting)
