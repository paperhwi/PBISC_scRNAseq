"""Real-vs-generated iLISI after pairwise Harmony.

Method (PBISC notes 2026-10-21 / 2026-10-13, ported from three_model_20261022/code/ilisi_lib.py and integration_report_20261012/code):
2,000 seurat_v3 HVG chosen on the reference real cells -> log1p CP10k -> scaling with reference mean/sd (clip 10) -> PCA 30 fitted on the
reference -> Harmony with batch = source -> iLISI = inverse Simpson index of the source labels in a perplexity-30 Gaussian neighbourhood
(k = 90), median over (a subsample of) cells. Maximum for two sources is 2. The ceiling is real held-out half B vs real half A.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from ._utils import align_genes, lognorm, to_csr


# ----------------------------------------------------------------------------- HVG
def _norm_var(X):
    """scanpy 'seurat_v3' normalised variance of raw counts (one batch). Needs scikit-misc."""
    from skmisc.loess import loess
    X = X.tocsr().astype(np.float64)
    N, G = X.shape
    mean = np.asarray(X.mean(0)).ravel()
    var = (np.asarray(X.multiply(X).mean(0)).ravel() - mean ** 2) * N / (N - 1)
    nc = var > 0
    est = np.zeros(G)
    m = loess(np.log10(mean[nc]), np.log10(var[nc]), span=0.3, degree=2)
    m.fit()
    est[nc] = m.outputs.fitted_values
    reg = np.sqrt(10 ** est)
    clip = reg * np.sqrt(N) + mean
    Xc = X.copy()
    Xc.data = np.minimum(Xc.data, clip[Xc.indices])
    s2 = np.asarray(Xc.multiply(Xc).sum(0)).ravel()
    s1 = np.asarray(Xc.sum(0)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        nv = (1 / ((N - 1) * reg ** 2)) * (N * mean ** 2 + s2 - 2 * s1 * mean)
    nv[~nc] = 0
    return np.nan_to_num(nv, nan=0.0, posinf=0.0)


def hvg_seurat_v3(blocks, n_top=2000):
    """Indices of the top highly variable genes. One block = plain seurat_v3; several blocks = batch-aware (n batches, median rank)."""
    blocks = [to_csr(b) for b in blocks]
    nv = [_norm_var(b) for b in blocks]
    if len(blocks) == 1:
        return np.argsort(-nv[0])[:n_top]
    ranks = np.vstack([np.argsort(np.argsort(-v)) for v in nv])
    nb = (ranks < n_top).sum(0)
    med = np.median(ranks, axis=0)
    return np.lexsort((med, -nb))[:n_top]


# ----------------------------------------------------------------------------- embedding
def reference_pca(ref_X, query_X_list, hvg, n_pcs=30, clip=10.0, seed=0):
    """PCA fitted on the reference only. Returns (Z_ref, [Z_query, ...]); all inputs restricted to the same gene columns."""
    def prep(X, mu=None, sd=None):
        D = lognorm(to_csr(X)[:, hvg]).toarray()
        if mu is None:
            mu, sd = D.mean(0), D.std(0)
            sd[sd == 0] = 1.0
        return np.clip((D - mu) / sd, -clip, clip), mu, sd
    A, mu, sd = prep(ref_X)
    pca = PCA(n_components=n_pcs, random_state=seed, svd_solver="randomized").fit(A)
    return pca.transform(A), [pca.transform(prep(Q, mu, sd)[0]) for Q in query_X_list]


def pairwise_harmony(Z, source, seed=0, max_iter=20):
    """Harmony with batch = source. Handles harmonypy returning (cells x PCs) or (PCs x cells)."""
    import harmonypy
    meta = pd.DataFrame({"source": np.asarray(source).astype(str)})
    ho = harmonypy.run_harmony(Z, meta, "source", random_state=seed, max_iter_harmony=max_iter, verbose=False)
    Zc = np.asarray(ho.Z_corr)
    return Zc if Zc.shape == Z.shape else Zc.T


# ----------------------------------------------------------------------------- LISI
def lisi_vec(d, idx, lab, perplexity=30, tol=1e-5, iters=50):
    """Per-cell LISI from kNN distances d and indices idx (same bisection on beta as the Harmony/LISI reference implementation)."""
    D = d.astype(np.float64) ** 2
    n = len(D)
    lab = pd.Categorical(lab).codes
    L = lab.max() + 1
    logU = np.log(perplexity)
    beta = np.ones(n); lo = np.full(n, -np.inf); hi = np.full(n, np.inf)
    act = np.ones(n, bool); Pf = np.ones_like(D) / D.shape[1]
    for _ in range(iters):
        a = np.where(act)[0]
        if len(a) == 0:
            break
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
    for l in range(L):
        w[:, l] = (Pf * (lab[idx] == l)).sum(1)
    return 1.0 / (w ** 2).sum(1)


def lisi_per_cell(Y, src, n_sub=10000, k=90, seed=0):
    """iLISI of a random subsample of cells. Returns (cell_index, lisi_values, own_source_neighbour_share)."""
    src = np.asarray(src)
    rng = np.random.default_rng(seed)
    sub = rng.choice(len(Y), size=min(n_sub, len(Y)), replace=False)
    k = min(k, len(Y))
    d, ix = NearestNeighbors(n_neighbors=k).fit(Y).kneighbors(Y[sub])
    v = lisi_vec(d, ix, src)
    own = (src[ix] == src[sub][:, None]).mean(1)
    return sub, v, own


def lisi_stats(Y, src, n_sub=10000, k=90, seed=0):
    """(median iLISI, mean iLISI, own-source neighbour share)."""
    _, v, own = lisi_per_cell(Y, src, n_sub, k, seed)
    return float(np.median(v)), float(np.mean(v)), float(own.mean())


def patient_bootstrap_median(values, patient, n_boot=1000, seed=0):
    """Median iLISI with a 95% CI from resampling patients (clusters) with replacement. The neighbourhoods stay fixed; only the
    cell-level iLISI values are resampled, so a patient's cells enter or leave together."""
    values, patient = np.asarray(values), np.asarray(patient).astype(str)
    ids = np.unique(patient)
    groups = {p: values[patient == p] for p in ids}
    rng = np.random.default_rng(seed)
    meds = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(ids, size=len(ids), replace=True)
        meds[b] = np.median(np.concatenate([groups[p] for p in pick]))
    return float(np.median(values)), float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


# ----------------------------------------------------------------------------- pipeline
def ilisi_pairwise(ref_X, ref_genes, query_X, query_genes, n_per_source=20000, n_hvg=2000, n_pcs=30, n_draws=3, seed=0,
                   ref_patient=None, query_patient=None, n_boot=0):
    """Reference vs one query source: iLISI before (PCA projection) and after Harmony.

    Both sources are subsampled to the same number of cells (iLISI is only comparable at balanced source sizes). The median iLISI is the
    mean over `n_draws` query draws (seeds seed..seed+n_draws-1) of a 10,000-cell median. With patient ids and n_boot > 0 a patient-bootstrap
    95% CI of the post-Harmony iLISI is added (taken from the first draw).
    """
    ia, ib, names = align_genes(ref_genes, query_genes)
    R, Q = to_csr(ref_X)[:, ia], to_csr(query_X)[:, ib]
    rng = np.random.default_rng(seed)
    n = min(n_per_source, R.shape[0], Q.shape[0])
    ri, qi = np.sort(rng.choice(R.shape[0], n, replace=False)), np.sort(rng.choice(Q.shape[0], n, replace=False))
    R, Q = R[ri], Q[qi]
    hvg = hvg_seurat_v3([R], n_top=min(n_hvg, R.shape[1]))
    Zr, (Zq,) = reference_pca(R, [Q], hvg, n_pcs=n_pcs, seed=seed)
    Z = np.vstack([Zr, Zq])
    src = np.array(["ref"] * n + ["query"] * n)
    Zh = pairwise_harmony(Z, src, seed=seed)
    out = {"n_per_source": n, "n_hvg": len(hvg)}
    for stage, Y in (("pre", Z), ("post", Zh)):
        s = [lisi_stats(Y, src, seed=seed + j) for j in range(n_draws)]
        out[f"ilisi_{stage}"] = float(np.mean([x[0] for x in s]))
        out[f"own_share_{stage}"] = float(np.mean([x[2] for x in s]))
    if n_boot and ref_patient is not None and query_patient is not None:
        pat = np.concatenate([np.asarray(ref_patient)[ri], np.asarray(query_patient)[qi]])
        sub, v, _ = lisi_per_cell(Zh, src, seed=seed)
        med, lo, hi = patient_bootstrap_median(v, pat[sub], n_boot=n_boot, seed=seed)
        out.update({"ilisi_post_boot_median": med, "ilisi_post_ci_lo": lo, "ilisi_post_ci_hi": hi})
    return out


def ilisi_report(ref_X, ref_genes, sources, ceiling, ref_patient=None, **kw):
    """iLISI of several sources against the same real reference, with the fraction of the real-vs-real ceiling reached.

    sources: {name: (X, genes)} or {name: (X, genes, patient)}; ceiling: name of the real held-out source (half B) in `sources`.
    fraction_of_ceiling = (iLISI_source - 1) / (iLISI_ceiling - 1), so 1 means indistinguishable from a second real sample and 0 means
    no mixing at all. Extra keyword arguments go to ilisi_pairwise. Returns a DataFrame with one row per source.
    """
    rows = []
    for name, tup in sources.items():
        X, g = tup[0], tup[1]
        pat = tup[2] if len(tup) > 2 else None
        res = ilisi_pairwise(ref_X, ref_genes, X, g, ref_patient=ref_patient, query_patient=pat, **kw)
        rows.append({"source": name, **res})
    df = pd.DataFrame(rows).set_index("source")
    c = df.loc[ceiling]
    for stage in ("pre", "post"):
        df[f"fraction_of_ceiling_{stage}"] = (df[f"ilisi_{stage}"] - 1) / (c[f"ilisi_{stage}"] - 1)
    return df
