"""Integration robustness (2026-10-21). Usage: python integ.py <coh>
Pairs: real A + {real B, floor, prev K4000}, 20,000 + 20,000 cells drawn from the ~60k bundles (seed 3). Same 2,000 real-A HVG.
Methods: Harmony (on real-A-fit PCA 30, = main pipeline), Scanorama (log1p CP10k HVG, 30-d), scVI (raw counts on HVG, batch = source, 30-d latent,
40 epochs, CPU), BBKNN (graph on PCA 30, batch = source). Metrics: iLISI of source (perplexity 30, k = 90; for BBKNN = weighted inverse
Simpson over the BBKNN graph), cLISI of six-type label (1 = pure), label ASW ((silhouette+1)/2, 8,000 cells). Real B pair = ceiling."""
import sys, json, time, numpy as np, pandas as pd, anndata as ad, scanpy as sc, harmonypy, scanorama, bbknn, scvi, torch
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from common import *
torch.set_num_threads(2); coh = sys.argv[1]; rng = np.random.default_rng(3); t0 = time.time()
R = load(coh, "ref"); hv = L.hvg_seurat_v3([R.X.tocsr()]); rows = []
def lisi_on(E, lab, k=90, n=8000):
    return L.lisi_stats(E, lab, n_sub=n, k=k, seed=0)[0]
for s in ("realU", "floor", "prev4k"):
    S = load(coh, s); ia = rng.choice(R.n_obs, 20000, replace=False); ib = rng.choice(S.n_obs, 20000, replace=False)
    A = ad.concat([R[ia], S[ib]], merge="same"); A.obs_names_make_unique(); src = np.r_[np.zeros(20000, int), np.ones(20000, int)].astype(str); A.obs["src"] = src
    six = A.obs.six_type.astype(str).to_numpy(); A = A[:, hv].copy(); A.layers["counts"] = A.X.copy()
    Yl = L.lognorm(A.X.tocsr()); Ar = Yl[:20000].toarray(); mu, sd = Ar.mean(0), Ar.std(0) + 1e-6
    pca = PCA(30, random_state=0).fit(np.clip((Ar - mu) / sd, -10, 10)); P = pca.transform(np.clip((Yl.toarray() - mu) / sd, -10, 10)).astype(np.float32)
    E = {"PCA (none)": P}
    ho = harmonypy.run_harmony(P.astype(np.float64), pd.DataFrame({"s": src}), ["s"], random_state=0, max_iter_harmony=20, verbose=False); Z = np.asarray(ho.Z_corr)
    E["Harmony"] = (Z.T if Z.shape[0] == 30 else Z).astype(np.float32)
    X0 = Yl[:20000].tocsr(); X1 = Yl[20000:].tocsr(); integ, _ = scanorama.integrate([X0, X1], [list(map(str, range(len(hv))))] * 2, dimred=30, verbose=0)
    E["Scanorama"] = np.vstack(integ).astype(np.float32)
    B = ad.AnnData(A.layers["counts"].copy(), obs=pd.DataFrame({"src": src}, index=A.obs_names)); scvi.model.SCVI.setup_anndata(B, batch_key="src")
    m = scvi.model.SCVI(B, n_latent=30, n_layers=2); m.train(max_epochs=40, accelerator="cpu", enable_progress_bar=False, batch_size=512)
    E["scVI"] = m.get_latent_representation().astype(np.float32); print(f"[{time.time()-t0:.0f}s] {coh} {s} scVI done", flush=True)
    for meth, emb in E.items():
        sub = rng.choice(len(src), 8000, replace=False)
        rows.append({"cohort": coh, "pair": s, "method": meth, "iLISI": lisi_on(emb, src), "cLISI": lisi_on(emb, six),
                     "label_ASW": float((silhouette_score(emb[sub], six[sub]) + 1) / 2)})
    G = ad.AnnData(np.zeros((len(src), 0), np.float32), obs=pd.DataFrame({"src": src}, index=A.obs_names)); G.obsm["X_pca"] = P
    bbknn.bbknn(G, batch_key="src", use_rep="X_pca", neighbors_within_batch=15)
    Cn = G.obsp["connectivities"].tocsr(); srcc = src.astype(int); ilisi_b, clisi_b = [], []
    for i in rng.choice(len(src), 8000, replace=False):
        a, b = Cn.indptr[i], Cn.indptr[i + 1]; j = Cn.indices[a:b]; w = Cn.data[a:b] / Cn.data[a:b].sum()
        pb = np.bincount(srcc[j], weights=w, minlength=2); ilisi_b.append(1 / (pb ** 2).sum())
        lab_ = pd.Categorical(six[j]).codes; pc = np.bincount(lab_, weights=w); clisi_b.append(1 / (pc ** 2).sum())
    rows.append({"cohort": coh, "pair": s, "method": "BBKNN (graph)", "iLISI": float(np.median(ilisi_b)), "cLISI": float(np.median(clisi_b)), "label_ASW": np.nan})
    print(pd.DataFrame(rows).tail(6).round(3).to_string(), flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/{coh}_integration_methods.csv", index=False)
