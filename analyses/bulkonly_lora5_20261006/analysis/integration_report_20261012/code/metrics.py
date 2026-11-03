"""Integration / clustering metrics used in the lab-meeting report.
kBET  (Buttner 2019): per-cell chi-square test of local batch composition (k neighbours) vs global; rejection rate.
LISI  (Korsunsky 2019): inverse Simpson index of label composition in a perplexity-30 Gaussian neighbourhood (k=90).
       iLISI on batch (higher = better mixing, max = #batches), cLISI on cell type (lower = better separation, min 1).
"""
import numpy as np, pandas as pd
from scipy.stats import chi2
from sklearn.neighbors import NearestNeighbors


def knn(X, k):
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=2).fit(X)
    d, i = nn.kneighbors(X)
    return d[:, 1:], i[:, 1:]


def kbet(X, batch, k=30, n_query=10000, alpha=0.05, seed=0, null_reps=5):
    rng = np.random.default_rng(seed)
    b = pd.Categorical(batch); codes = b.codes; B = len(b.categories)
    p = np.bincount(codes, minlength=B) / len(codes)
    q = rng.choice(len(X), size=min(n_query, len(X)), replace=False)
    nn = NearestNeighbors(n_neighbors=k, n_jobs=2).fit(X); _, idx = nn.kneighbors(X[q])
    def rate(c):
        obs = np.stack([np.bincount(row, minlength=B) for row in c[idx]])
        exp = k * p; stat = ((obs - exp) ** 2 / exp).sum(1)
        return float((chi2.sf(stat, B - 1) < alpha).mean())
    obs_rate = rate(codes)
    null = np.mean([rate(rng.permutation(codes)) for _ in range(null_reps)])
    return obs_rate, float(null)


def _lisi_from_knn(d, idx, labels, perplexity=30, tol=1e-5, iters=50):
    lab = pd.Categorical(labels).codes; L = lab.max() + 1
    n = len(d); out = np.zeros(n); logU = np.log(perplexity)
    for i in range(n):
        di = d[i] ** 2; beta = 1.0; lo, hi = -np.inf, np.inf
        for _ in range(iters):
            P = np.exp(-di * beta); s = P.sum()
            if s == 0: H = 0.0; P = np.ones_like(di) / len(di)
            else:
                H = np.log(s) + beta * (di * P).sum() / s; P = P / s
            diff = H - logU
            if abs(diff) < tol: break
            if diff > 0: lo = beta; beta = beta * 2 if hi == np.inf else (beta + hi) / 2
            else: hi = beta; beta = beta / 2 if lo == -np.inf else (beta + lo) / 2
        w = np.bincount(lab[idx[i]], weights=P, minlength=L)
        out[i] = 1.0 / (w ** 2).sum()
    return out


def lisi(X, labels_dict, n_sub=15000, k=90, seed=0):
    rng = np.random.default_rng(seed)
    sub = rng.choice(len(X), size=min(n_sub, len(X)), replace=False)
    nn = NearestNeighbors(n_neighbors=k, n_jobs=2).fit(X)
    d, idx = nn.kneighbors(X[sub])
    return {name: float(np.median(_lisi_from_knn(d, idx, np.asarray(lab)))) for name, lab in labels_dict.items()}


def leiden_on(X, k=30, resolutions=(0.5, 1.0), seed=2021):
    import scanpy as sc, anndata as ad, igraph as ig, leidenalg, scipy.sparse as sp
    a = ad.AnnData(np.zeros((len(X), 0), np.float32)); a.obsm["X"] = X
    sc.pp.neighbors(a, use_rep="X", n_neighbors=k, random_state=seed)
    U = sp.triu(a.obsp["connectivities"], k=1).tocoo()
    g = ig.Graph(n=len(X), edges=np.column_stack([U.row, U.col]), directed=False); g.es["weight"] = U.data
    out = {}
    for r in resolutions:
        out[r] = np.array(leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight",
                                                   resolution_parameter=r, seed=seed, n_iterations=2).membership)
    return out
