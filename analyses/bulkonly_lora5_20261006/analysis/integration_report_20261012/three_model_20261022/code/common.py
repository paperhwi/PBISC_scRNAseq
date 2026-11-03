import os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
IN = "/home/claude/work/s3/big"
BJ = "/mnt/user-data/outputs/ilisi_floor_20261021/data/big_joint"; TR = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"
OUT = "/mnt/user-data/outputs/s3/results"; VAL = "/mnt/user-data/outputs/val/results"
LIN = {"T cell": "T", "NK cell": "NK", "B cell": "B", "Monocyte": "Myeloid", "DC": "Myeloid", "Other": "Other"}
def load(coh, src):
    return ad.read_h5ad(f"{IN}/{coh}_{src}.h5ad")
def thin_to(B, G, seed):
    """B counts thinned per six-type to depths drawn from G cells of the same six-type"""
    rng = np.random.default_rng(seed); X = B.X.tocsr().copy(); tot = np.asarray(X.sum(1)).ravel(); gt = np.asarray(G.X.sum(1)).ravel()
    lab = G.obs.six_type.astype(str).to_numpy(); pool = {c: gt[lab == c] for c in np.unique(lab)}
    tgt = np.array([rng.choice(pool.get(c, gt)) for c in B.obs.six_type.astype(str)], float); keep = np.minimum(1.0, tgt / np.maximum(tot, 1))
    X.data = rng.binomial(X.data.astype(np.int64), np.repeat(keep, np.diff(X.indptr))).astype(np.float32); X.eliminate_zeros()
    return ad.AnnData(X, obs=B.obs.copy(), var=B.var.copy())
def bundle(coh, with_thin=True, sources=("ref", "realU", "floor", "prev4k", "scsemi"), scsemiF=False):
    D = {s: load(coh, s) for s in sources}
    if with_thin:
        if "floor" in D: D["realT1f"] = thin_to(D["realU"], D["floor"], 20261025)
        if "prev4k" in D: D["realT1p"] = thin_to(D["realU"], D["prev4k"], 20261026)
        if "scsemi" in D: D["realT1s"] = thin_to(D["realU"], D["scsemi"], 20261027)
    if scsemiF and "scsemi" in D:
        a = D["scsemi"]; D["scsemiF"] = ad.AnnData(a.layers["floor_rule"].tocsr(), obs=a.obs.copy(), var=a.var.copy())
    for s, a in D.items(): a.obs["lineage"] = [LIN.get(x, "Other") for x in a.obs.six_type.astype(str)]
    return D
def transferred(coh, src):
    return pd.read_csv(f"{TR}/{coh}_{src}_transfer.csv.gz", index_col=0, keep_default_na=False)
