from pathlib import Path

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix


ROOT = Path(__file__).parent
SAMPLES = [
    ("GSM8571091_Abd-S101", "adult_sepsis", "abdominal"),
    ("GSM8571092_Abd-S102", "adult_sepsis", "abdominal"),
    ("GSM8571093_Abd-S103", "adult_sepsis", "abdominal"),
    ("GSM8571082_Abd-PS04", "pediatric_sepsis", "abdominal"),
    ("GSM8571083_Abd-PS05", "pediatric_sepsis", "abdominal"),
]


def load_sample(sample_id, group, site):
    folder = ROOT / sample_id
    matrix = mmread(folder / f"{sample_id}.matrix.mtx.gz").tocsr().T.tocsr()
    features = pd.read_csv(folder / f"{sample_id}.features.tsv.gz", sep="\t", header=None)
    barcodes = pd.read_csv(folder / f"{sample_id}.barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    symbols = features.iloc[:, 1].astype(str) if features.shape[1] > 1 else features.iloc[:, 0].astype(str)
    gene_names = pd.Index(symbols).astype(str)
    gene_names.name = None
    gene_names = ad.utils.make_index_unique(gene_names)
    var = pd.DataFrame(index=gene_names)
    var["gene_id"] = features.iloc[:, 0].astype(str).to_numpy()
    var["mt"] = var.index.str.upper().str.startswith("MT-")
    obs = pd.DataFrame(index=barcodes)
    obs.index.name = None
    obs["sample_id"] = sample_id
    obs["cohort_group"] = group
    obs["infection_site"] = site
    return ad.AnnData(X=matrix, obs=obs, var=var)


adatas = [load_sample(*x) for x in SAMPLES]
adata = ad.concat(adatas, join="inner", index_unique="-")
adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
qc_before = adata.n_obs
adata = adata[adata.obs["n_genes_by_counts"] > 200].copy()
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
adata = adata[:, ~adata.var["mt"]].copy()
qc_after = adata.n_obs

sc.pp.normalize_total(adata, target_sum=10_000)
sc.pp.log1p(adata)
adata.layers["log1p_norm"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2_000, flavor="seurat", batch_key="sample_id")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.regress_out(adata, ["n_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=10, svd_solver="arpack")

ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, vars_use=["sample_id"], max_iter_harmony=20)
adata.obsm["X_pca_harmony"] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=10)
sc.tl.umap(adata, random_state=0)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", random_state=0)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", layer="log1p_norm", use_raw=False)

out = ROOT / "results"
out.mkdir(exist_ok=True)
sc.settings.figdir = str(out / "figures")
Path(sc.settings.figdir).mkdir(exist_ok=True)
adata.write_h5ad(out / "sepsis_primary_5samples_harmony_clusters.h5ad", compression="gzip")
adata.obs.groupby(["sample_id", "cohort_group"], observed=True).size().rename("cells_after_qc").to_csv(out / "cell_counts_by_sample.csv")
pd.DataFrame({"metric": ["cells_before_qc", "cells_after_qc", "n_genes_hvg", "n_pcs", "n_neighbors", "leiden_resolution", "n_clusters"], "value": [qc_before, qc_after, 2000, 10, 15, 0.5, adata.obs["leiden"].nunique()]}).to_csv(out / "integration_parameters.csv", index=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(out / "cluster_markers_wilcoxon.csv", index=False)
sc.pl.umap(adata, color=["sample_id", "cohort_group", "leiden"], ncols=1, show=False, save="_sepsis_primary_5samples.png")
