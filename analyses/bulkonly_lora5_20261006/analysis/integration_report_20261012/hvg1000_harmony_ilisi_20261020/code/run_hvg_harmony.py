"""Joint Harmony integration of real sepsis (Ye 2026 recap) + PBISC K=1,000 + PBISC K=4,000 generated cells,
HVG 1,000 vs HVG 2,000. Same recipe as integration_report code/c_joint.py, only n_top_genes changes:
panel genes shared by both sources, log1p CP10k on panel totals, seurat_v3 HVG on real counts, scale with real mean/sd,
PCA 30 fit on real and generated projected ('pre'), Harmony batch = source ('post'), iLISI (perplexity 30, k = 90, median).
Inputs: patient_pair_ilisi_20261018/_extract (real 800 cells/patient, generated 800 cells/sample, 5 draws pooled);
here 400 real cells per patient (133 patients, as in c_joint.py) and 400 generated cells per bulk-only sample per K.
usage: py -3.12 run_hvg_harmony.py"""
import os, sys, json, time, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, harmonypy
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
IR = os.path.dirname(OUT); EX = IR + "/patient_pair_ilisi_20261018/_extract"
sys.path.insert(0, IR + "/code"); import metrics as M
rng = np.random.default_rng(0); PER = 400


def load(name, by):
    z = np.load(f"{EX}/{name}", allow_pickle=True)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.DataFrame({k: z[k] for k in z.files if k not in ("data", "indices", "indptr", "shape", "genes", "rows")})
    take = np.sort(np.concatenate([rng.permutation(np.where(obs[by] == s)[0])[:PER] for s in pd.unique(obs[by])]))
    return X[take], obs.iloc[take].reset_index(drop=True), z["genes"].astype(str)


Xr, ro, genes = load("sepsis_patients_cells.npz", "sample")
gen = {K: load(f"sepsis_{K}_samples_cells.npz", "sample") for K in ["k1000", "k4000"]}
for K in gen: assert (gen[K][2] == genes).all()
print("cells real / k1000 / k4000:", Xr.shape[0], gen["k1000"][0].shape[0], gen["k4000"][0].shape[0], "genes", len(genes), flush=True)


def lognorm(X):
    X = X.tocsr().astype(np.float32); tot = np.asarray(X.sum(1)).ravel()
    Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y


Lr = lognorm(Xr); Lg = {K: lognorm(gen[K][0]) for K in gen}
src = np.array(["real"] * Xr.shape[0] + ["k1000"] * gen["k1000"][0].shape[0] + ["k4000"] * gen["k4000"][0].shape[0])
rmaj = ro["celltype_major_ye2026"].to_numpy()
res, umaps = [], {}
for nhv in [1000, 2000]:
    t0 = time.time()
    a = ad.AnnData(Xr.copy()); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=nhv)
    hv = np.where(a.var.highly_variable.to_numpy())[0]
    mu = np.asarray(Lr[:, hv].mean(0)).ravel(); sd = np.sqrt(np.asarray(Lr[:, hv].multiply(Lr[:, hv]).mean(0)).ravel() - mu ** 2) + 1e-6
    scale = lambda Y: np.clip((Y[:, hv].toarray() - mu) / sd, -10, 10)
    pca = PCA(30, random_state=0).fit(scale(Lr))
    pre = np.vstack([pca.transform(scale(Lr))] + [pca.transform(scale(Lg[K])) for K in gen]).astype(np.float32)
    ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"src": src}), ["src"], random_state=0, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
    for t, X in [("pre", pre), ("post", post)]:
        res.append(dict(n_hvg=nhv, stage=t, comparison="real + K1k + K4k (3 sources)", iLISI=M.lisi(X, {"s": src}, n_sub=12000)["s"], max=3))
        for K in gen:
            m = np.isin(src, ["real", K])
            res.append(dict(n_hvg=nhv, stage=t, comparison=f"real vs {K.upper().replace('K1000','K1k').replace('K4000','K4k')}",
                            iLISI=M.lisi(X[m], {"s": src[m]}, n_sub=10000)["s"], max=2))
    # transfer major labels to generated cells (kNN k = 30, distance-weighted) for UMAP colouring only
    nr = len(rmaj); nn = NearestNeighbors(n_neighbors=30).fit(post[:nr]); d, idx = nn.kneighbors(post[nr:]); w = 1 / (d + 1e-3)
    cats = pd.Categorical(rmaj); sc_ = np.zeros((len(idx), len(cats.categories)))
    for j in range(idx.shape[1]): np.add.at(sc_, (np.arange(len(idx)), cats.codes[idx[:, j]]), w[:, j])
    major = np.concatenate([rmaj, cats.categories[sc_.argmax(1)]])
    take = np.concatenate([rng.permutation(np.where(src == s)[0])[:30000] for s in ["real", "k1000", "k4000"]])
    aa = ad.AnnData(np.zeros((len(take), 0), np.float32)); aa.obsm["X"] = post[take]
    sc.pp.neighbors(aa, use_rep="X", n_neighbors=30, random_state=0); sc.tl.umap(aa, random_state=0, min_dist=0.3)
    umaps[nhv] = dict(umap=aa.obsm["X_umap"], src=src[take], major=major[take].astype(str))
    np.savez_compressed(f"{OUT}/results/joint_umap_hvg{nhv}.npz", **umaps[nhv])
    np.save(f"{OUT}/results/X_pca_harmony_hvg{nhv}.npy", post)
    pd.Series(a.var_names[hv] if False else genes[hv]).to_csv(f"{OUT}/results/hvg{nhv}_genes.csv", index=False, header=["gene"])
    print(nhv, "done", round(time.time() - t0), "s", flush=True)
R = pd.DataFrame(res); R.to_csv(f"{OUT}/results/ilisi_hvg1000_vs_hvg2000.csv", index=False); print(R.to_string())
pd.DataFrame({"src": src, "sample": np.concatenate([ro["sample"], gen["k1000"][1]["sample"], gen["k4000"][1]["sample"]])}).to_csv(
    f"{OUT}/results/cells_obs.csv.gz", index=False)
json.dump({"n_real": int((src == "real").sum()), "n_k1000": int((src == "k1000").sum()), "n_k4000": int((src == "k4000").sum()),
           "panel_genes": int(len(genes)), "cells_per_unit": PER}, open(f"{OUT}/results/run_info.json", "w"), indent=1)
