"""iLISI of one source vs the reference real cells after Harmony. Usage: python ilisi_run.py <coh> <source> [V0,V1]
V0 = reference HVG + reference PCA (the report's main pipeline); V1 = joint HVG (batch_key = source) + reference PCA.
Writes results/<coh>_ilisi_<source>.json (one file per run: parallel runs must not share a file) and data/<coh>_joint_<source>_V0.h5ad (PCs pre/post Harmony, V0 only)."""
import sys, os, json, time, numpy as np, pandas as pd, anndata as ad, harmonypy, scipy.sparse as sp
from sklearn.decomposition import PCA
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
coh, name = sys.argv[1], sys.argv[2]; variants = (sys.argv[3] if len(sys.argv) > 3 else "V0,V1").split(",")
O = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); t0 = time.time()
def log(*a): print(f"[{time.time()-t0:5.0f}s]", *a, flush=True)
R = ad.read_h5ad(f"{O}/data/{coh}_ref.h5ad"); S = ad.read_h5ad(f"{O}/data/{coh}_{name}.h5ad")
XR, XS = R.X.tocsr(), S.X.tocsr(); nR, nS = XR.shape[0], XS.shape[0]; YR, YS = L.lognorm(XR), L.lognorm(XS)
jf = f"{O}/results/{coh}_ilisi_{name}.json"; res = json.load(open(jf)) if os.path.exists(jf) else {}
lin = np.r_[R.obs.lineage.to_numpy(), S.obs.lineage.to_numpy()]; src = np.r_[np.zeros(nR, int), np.ones(nS, int)]
for v in variants:
    hv = L.hvg_seurat_v3([XR]) if v == "V0" else L.hvg_seurat_v3([XR, XS])
    A = YR[:, hv].toarray().astype(np.float32); mu = A.mean(0); sd = A.std(0) + 1e-6
    sc_ = lambda Y: np.clip((Y[:, hv].toarray().astype(np.float32) - mu) / sd, -10, 10)
    pca = PCA(30, random_state=0).fit(sc_(YR)); pre = np.vstack([pca.transform(sc_(YR)), pca.transform(sc_(YS))]).astype(np.float32)
    ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"s": src.astype(str)}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
    r = {"n_ref": int(nR), "n_source": int(nS), "hvg_overlap_with_ref_hvg": int(len(np.intersect1d(hv, L.hvg_seurat_v3([XR]))))}
    for tag, Y in (("pre", pre), ("post", post)):
        vals = [L.lisi_stats(Y, src, seed=s) for s in (0, 1, 2)]
        r[f"iLISI_{tag}"] = float(np.mean([x[0] for x in vals])); r[f"iLISI_{tag}_sd_query"] = float(np.std([x[0] for x in vals])); r[f"own_share_{tag}"] = float(np.mean([x[2] for x in vals]))
        for c in ["T", "NK", "B", "Myeloid"]:
            m = lin == c
            if (m & (src == 0)).sum() > 500 and (m & (src == 1)).sum() > 500: r[f"iLISI_{tag}_{c}"] = L.lisi_stats(Y[m], src[m], n_sub=5000)[0]
    res.setdefault(v, {})[name] = r; json.dump(res, open(jf, "w"), indent=1); log(coh, name, v, {k: round(x, 3) for k, x in r.items() if k.startswith(("iLISI_pre", "iLISI_post", "own"))})
    if v == "V0":
        o = pd.concat([R.obs.assign(source="ref"), S.obs.assign(source=name)]); o.index = [f"{s}_{i}" for i, s in enumerate(o.source)]
        a_ = ad.AnnData(np.zeros((len(o), 0), np.float32), obs=o); a_.obsm["X_pca"] = pre; a_.obsm["X_pca_harmony"] = post
        a_.write_h5ad(f"{O}/data/{coh}_joint_{name}_V0.h5ad", compression="gzip")
