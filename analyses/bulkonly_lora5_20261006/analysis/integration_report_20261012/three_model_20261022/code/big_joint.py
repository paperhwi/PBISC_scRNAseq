"""Pairwise Harmony of each source with real half A on the larger bundles (~60,000 cells per source), then
iLISI (robustness re-run on a second subsample), kNN label transfer from real half A, joint Leiden clusters (cluster matching) and UMAP.
Usage: python big_joint.py <coh> <source>   (source = realU | floor | prev4k)
Pipeline = ilisi_run.py V0: 2,000 seurat_v3 HVG on half A, scale with half A mean/sd (clip 10), PCA 30 fit on half A, Harmony (batch = source, max_iter 20, seed 0)."""
import sys, os, json, time, numpy as np, pandas as pd, anndata as ad, harmonypy, scanpy as sc
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
coh, name = sys.argv[1], sys.argv[2]; t0 = time.time()
IN = "/home/claude/work/s3/big"
OUT = "/mnt/user-data/outputs/ilisi_floor_20261021"
def log(*a): print(f"[{time.time()-t0:5.0f}s]", coh, name, *a, flush=True)
R = ad.read_h5ad(f"{IN}/{coh}_ref.h5ad"); S = ad.read_h5ad(f"{IN}/{coh}_{name}.h5ad")
XR, XS = R.X.tocsr(), S.X.tocsr(); nR, nS = XR.shape[0], XS.shape[0]; YR, YS = L.lognorm(XR), L.lognorm(XS)
hv = L.hvg_seurat_v3([XR]); A = YR[:, hv].toarray().astype(np.float32); mu = A.mean(0); sd = A.std(0) + 1e-6; del A
scl = lambda Y: np.clip((Y[:, hv].toarray().astype(np.float32) - mu) / sd, -10, 10)
pca = PCA(30, random_state=0).fit(scl(YR)); pre = np.vstack([pca.transform(scl(YR)), pca.transform(scl(YS))]).astype(np.float32)
src = np.r_[np.zeros(nR, int), np.ones(nS, int)]
ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"s": src.astype(str)}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32); log("harmony done")
res = {"n_ref": int(nR), "n_source": int(nS)}
for tag, Y in (("pre", pre), ("post", post)):
    v = [L.lisi_stats(Y, src, seed=s) for s in (0, 1, 2)]
    res[f"iLISI_{tag}"] = float(np.mean([x[0] for x in v])); res[f"own_share_{tag}"] = float(np.mean([x[2] for x in v]))
log("iLISI", round(res["iLISI_pre"], 3), round(res["iLISI_post"], 3))
# kNN label transfer (k = 30) from half A in the Harmony space
fine_r = R.obs.fine.astype(str).to_numpy(); major_r = R.obs.major.astype(str).to_numpy(); six_r = R.obs.six_type.astype(str).to_numpy()
f2m = pd.crosstab(fine_r, major_r).idxmax(1).to_dict(); f2s = pd.crosstab(fine_r, six_r).idxmax(1).to_dict()
d, ix = NearestNeighbors(n_neighbors=30, n_jobs=2).fit(post[:nR]).kneighbors(post[nR:])
cats, codes = np.unique(fine_r, return_inverse=True); C = np.zeros((nS, len(cats)), np.int16)
for j in range(30): np.add.at(C, (np.arange(nS), codes[ix[:, j]]), 1)
tf = cats[C.argmax(1)]; conf = C.max(1) / 30.0
q = S.obs.copy(); q["fine_tr"] = tf; q["major_tr"] = [f2m[x] for x in tf]; q["six_tr"] = [f2s[x] for x in tf]; q["conf_tr"] = conf
res["conf_fine_median"] = float(np.median(conf)); res["six_agree_with_own_label"] = float((q.six_tr.to_numpy() == q.six_type.astype(str).to_numpy()).mean())
if name == "realU": res["acc_fine_vs_true"] = float((q.fine_tr.to_numpy() == q.fine.astype(str).to_numpy()).mean())
# six-type agreement per own label (generated: frozen scANVI label)
res["six_agree_by_type"] = q.assign(ok=q.six_tr.to_numpy() == q.six_type.astype(str).to_numpy()).groupby("six_type").ok.mean().round(4).to_dict()
log("transfer conf", round(res["conf_fine_median"], 3), "six agree", round(res["six_agree_with_own_label"], 3))
# joint Leiden (cluster matching)
J = ad.AnnData(np.zeros((nR + nS, 0), np.float32)); J.obsm["X_h"] = post
sc.pp.neighbors(J, n_neighbors=15, use_rep="X_h", random_state=0); sc.tl.leiden(J, resolution=1.0, flavor="igraph", n_iterations=2, random_state=0, key_added="leiden")
lc = J.obs.leiden.astype(int).to_numpy(); ct = pd.crosstab(lc, src); ct.columns = ["ref", "src"]
fr = ct / ct.sum(0); w = fr.src / (fr.src + fr.ref)       # size-normalised source share per cluster (0.5 = balanced)
maj = pd.Series(fine_r).groupby(lc[:nR]).agg(lambda x: x.value_counts().index[0]) if nR else None
cl = pd.DataFrame({"n_ref": ct.ref, "n_src": ct.src, "frac_ref": fr.ref, "frac_src": fr.src, "src_share_norm": w})
cl["ref_major_fine"] = [maj.get(i, "") if i in maj.index else "(no real cells)" for i in cl.index]
cl.to_csv(f"{OUT}/results/big/{coh}_{name}_leiden_clusters.csv", index_label="leiden")
res.update({"n_clusters": int(len(cl)), "src_cells_in_mixed_clusters": float(cl.frac_src[(w >= 0.2) & (w <= 0.8)].sum()),
            "ref_cells_in_mixed_clusters": float(cl.frac_ref[(w >= 0.2) & (w <= 0.8)].sum()),
            "src_cells_in_src_only_clusters": float(cl.frac_src[w > 0.9].sum()), "ref_cells_in_ref_only_clusters": float(cl.frac_ref[w < 0.1].sum()),
            "cluster_freq_pearson": float(np.corrcoef(fr.ref, fr.src)[0, 1]), "cluster_freq_JSD": float(0.5 * (np.nansum(fr.ref * np.log2(2 * fr.ref / (fr.ref + fr.src))) + np.nansum(fr.src * np.log2(np.where(fr.src > 0, 2 * fr.src / (fr.ref + fr.src), 1)))))})
q["leiden_joint"] = lc[nR:]
log("leiden", res["n_clusters"], "src in mixed", round(res["src_cells_in_mixed_clusters"], 3), "r", round(res["cluster_freq_pearson"], 3))
# UMAP on a 25k + 25k subsample, post (all sources) and pre (generated only)
rng = np.random.default_rng(0); sub = np.r_[rng.choice(nR, 25000, replace=False), nR + rng.choice(nS, 25000, replace=False)]
um = {}
for tag, Y in (("post", post), ("pre", pre)):
    if tag == "pre" and name == "realU": continue
    U = ad.AnnData(np.zeros((len(sub), 0), np.float32)); U.obsm["X_h"] = Y[sub]
    sc.pp.neighbors(U, n_neighbors=15, use_rep="X_h", random_state=0); sc.tl.umap(U, random_state=0); um[tag] = U.obsm["X_umap"].astype(np.float32)
lab_sub = np.r_[fine_r, tf][sub]; six_sub = np.r_[six_r, S.obs.six_type.astype(str).to_numpy()][sub]
np.savez_compressed(f"{OUT}/results/big/{coh}_{name}_umap.npz", sub=sub, src=src[sub], fine=lab_sub, six=six_sub, leiden=lc[sub], **{f"umap_{k}": v for k, v in um.items()})
log("umap done")
q.to_csv(f"{OUT}/results/big/{coh}_{name}_transfer.csv.gz")
o = pd.concat([R.obs.assign(source="ref", fine_tr=R.obs.fine, conf_tr=1.0), q.assign(source=name)]); o.index = [f"{s}_{i}" for i, s in enumerate(o.source)]
o["leiden_joint"] = lc.astype(str); o = o.astype({c: str for c in o.columns if o[c].dtype == object})
a_ = ad.AnnData(np.zeros((len(o), 0), np.float32), obs=o); a_.obsm["X_pca"] = pre; a_.obsm["X_pca_harmony"] = post
a_.write_h5ad(f"{OUT}/data/big_joint/{coh}_joint_{name}.h5ad", compression="gzip")
json.dump(res, open(f"{OUT}/results/big/{coh}_{name}_summary.json", "w"), indent=1); log("saved", json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in res.items() if not isinstance(v, dict)}))
