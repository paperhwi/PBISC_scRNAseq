"""Helpers for the real-vs-generated iLISI (2026-10-21). Methods follow integration_report_20261012/code (c_joint.py, i_hvg_pca.py, metrics.py):
2,000 seurat_v3 HVG on the reference real cells, PCA 30 fitted on the reference (scaled with reference mean/sd, clip 10), pairwise Harmony
(batch = source), iLISI = inverse Simpson index of source labels in a perplexity-30 Gaussian neighbourhood (k = 90), median over cells."""
import numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

def lognorm(X):
    t = np.asarray(X.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(t, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y

def _norm_var(X):
    """scanpy 'seurat_v3' normalised variance of raw counts (one batch)."""
    from skmisc.loess import loess
    X = X.tocsr().astype(np.float64); N, G = X.shape
    mean = np.asarray(X.mean(0)).ravel(); var = (np.asarray(X.multiply(X).mean(0)).ravel() - mean ** 2) * N / (N - 1)
    nc = var > 0; est = np.zeros(G); m = loess(np.log10(mean[nc]), np.log10(var[nc]), span=0.3, degree=2); m.fit()
    est[nc] = m.outputs.fitted_values; reg = np.sqrt(10 ** est)
    clip = reg * np.sqrt(N) + mean; Xc = X.copy(); Xc.data = np.minimum(Xc.data, clip[Xc.indices])
    s2 = np.asarray(Xc.multiply(Xc).sum(0)).ravel(); s1 = np.asarray(Xc.sum(0)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        nv = (1 / ((N - 1) * reg ** 2)) * (N * mean ** 2 + s2 - 2 * s1 * mean)
    nv[~nc] = 0; return np.nan_to_num(nv, nan=0.0, posinf=0.0)

def hvg_seurat_v3(blocks, n_top=2000):
    """blocks = [counts of batch 1, counts of batch 2, ...]; one block = plain seurat_v3, several = batch_key (median rank, n batches)."""
    nv = [_norm_var(b) for b in blocks]; ranks = np.vstack([np.argsort(np.argsort(-v)) for v in nv])
    if len(blocks) == 1: return np.argsort(-nv[0])[:n_top]
    nb = (ranks < n_top).sum(0); med = np.median(ranks, axis=0); order = np.lexsort((med, -nb)); return order[:n_top]

def lisi_vec(d, idx, lab, perplexity=30, tol=1e-5, iters=50):
    """Vectorised twin of metrics._lisi_from_knn (same bisection on beta, same stopping rule)."""
    D = d.astype(np.float64) ** 2; n = len(D); lab = pd.Categorical(lab).codes; L = lab.max() + 1; logU = np.log(perplexity)
    beta = np.ones(n); lo = np.full(n, -np.inf); hi = np.full(n, np.inf); act = np.ones(n, bool); Pf = np.ones_like(D) / D.shape[1]
    for _ in range(iters):
        a = np.where(act)[0]
        if len(a) == 0: break
        Pm = np.exp(-D[a] * beta[a, None]); s = Pm.sum(1); ok = s > 0
        H = np.zeros(len(a)); H[ok] = np.log(s[ok]) + beta[a][ok] * (D[a][ok] * Pm[ok]).sum(1) / s[ok]
        Pm[ok] /= s[ok, None]; Pm[~ok] = 1 / D.shape[1]; Pf[a] = Pm
        diff = H - logU; conv = np.abs(diff) < tol; up = (~conv) & (diff > 0); dn = (~conv) & (diff <= 0)
        ba, la, ha = beta[a], lo[a], hi[a]
        nb = ba.copy(); la2, ha2 = la.copy(), ha.copy()
        la2[up] = ba[up]; nb[up] = np.where(np.isinf(ha[up]), ba[up] * 2, (ba[up] + ha[up]) / 2)
        ha2[dn] = ba[dn]; nb[dn] = np.where(np.isinf(la[dn]), ba[dn] / 2, (ba[dn] + la[dn]) / 2)
        beta[a], lo[a], hi[a] = nb, la2, ha2; act[a[conv]] = False
    w = np.zeros((n, L))
    for l in range(L): w[:, l] = (Pf * (lab[idx] == l)).sum(1)
    return 1.0 / (w ** 2).sum(1)

def lisi_stats(Y, src, n_sub=10000, k=90, seed=0):
    rng = np.random.default_rng(seed); sub = rng.choice(len(Y), size=min(n_sub, len(Y)), replace=False)
    d, ix = NearestNeighbors(n_neighbors=k, n_jobs=2).fit(Y).kneighbors(Y[sub]); v = lisi_vec(d, ix, np.asarray(src))
    own = (np.asarray(src)[ix] == np.asarray(src)[sub][:, None]).mean(1)
    return float(np.median(v)), float(np.mean(v)), float(own.mean())
