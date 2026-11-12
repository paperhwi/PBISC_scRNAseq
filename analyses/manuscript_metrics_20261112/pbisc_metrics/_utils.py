"""Small helpers shared by the metric modules."""
import numpy as np
import pandas as pd
import scipy.sparse as sp


def to_csr(X, dtype=np.float64):
    """Any array-like / sparse matrix -> CSR with the requested dtype."""
    if sp.issparse(X):
        return X.tocsr().astype(dtype)
    return sp.csr_matrix(np.asarray(X), dtype=dtype)


def depth(X):
    """Per-cell total UMI."""
    return np.asarray(to_csr(X).sum(1)).ravel()


def cp10k(X):
    """Counts -> counts per 10,000 (sparse, per cell)."""
    X = to_csr(X)
    t = np.maximum(depth(X), 1.0)
    return (sp.diags(1e4 / t) @ X).tocsr()


def lognorm(X):
    """log1p(CP10k), sparse."""
    Y = cp10k(X)
    Y.data = np.log1p(Y.data)
    return Y


def dedupe_genes(genes):
    """Boolean keep-mask: first occurrence of each gene symbol is kept (the PBISC panels contain duplicated symbols)."""
    return ~pd.Index(np.asarray(genes).astype(str)).duplicated(keep="first")


def align_genes(genes_a, genes_b):
    """Indices of the shared genes (first occurrence) in a and b, in the order of a. Returns (idx_a, idx_b, names)."""
    ga = pd.Index(np.asarray(genes_a).astype(str))
    gb = pd.Index(np.asarray(genes_b).astype(str))
    ka, kb = ~ga.duplicated(keep="first"), ~gb.duplicated(keep="first")
    pos_b = pd.Series(np.arange(len(gb))[kb], index=gb[kb])
    ia = np.array([i for i in np.arange(len(ga))[ka] if ga[i] in pos_b.index], dtype=int)
    names = ga[ia]
    ib = pos_b.loc[names].to_numpy()
    return ia, ib, np.asarray(names)


def thin_to_depth(X, target, seed=0):
    """Binomial thinning of each cell to a target UMI depth (scalar or per-cell array).

    Cells already below the target are left unchanged (keep probability capped at 1). This is the 'depth-matched real' control
    used to compare real and generated cells at the same sequencing depth.
    """
    X = to_csr(X, np.float64)
    d = depth(X)
    tgt = np.broadcast_to(np.asarray(target, dtype=np.float64), d.shape)
    p = np.minimum(1.0, tgt / np.maximum(d, 1.0))
    rng = np.random.default_rng(seed)
    out = X.copy()
    row_p = np.repeat(p, np.diff(out.indptr))
    out.data = rng.binomial(np.round(out.data).astype(np.int64), row_p).astype(np.float64)
    out.eliminate_zeros()
    return out
