"""Compact on-disk CSR store (uint16 counts, uint16 gene index, int64 indptr) read through memmaps,
plus streaming versions of the Seurat / Scanpy preprocessing used by the two papers."""
import os, json, gc
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, svds


class Store:
    def __init__(self, path):
        self.path = path
        m = json.load(open(f"{path}/meta.json"))
        self.n_cells, self.n_genes = m["n_cells"], m["n_genes"]
        self.data = np.load(f"{path}/data.npy", mmap_mode="r")
        self.indices = np.load(f"{path}/indices.npy", mmap_mode="r")
        self.indptr = np.load(f"{path}/indptr.npy")
        self.var = pd.read_csv(f"{path}/var.csv")
        self.obs = pd.read_parquet(f"{path}/obs.parquet")
        self.ncount = np.load(f"{path}/ncount.npy")

    def rows(self, idx):
        """CSR float32 counts for sorted row indices idx (contiguous runs read efficiently)."""
        idx = np.asarray(idx)
        if len(idx) == 0:
            return sp.csr_matrix((0, self.n_genes), dtype=np.float32)
        if np.all(np.diff(idx) == 1):
            a, b = self.indptr[idx[0]], self.indptr[idx[-1] + 1]
            ip = self.indptr[idx[0]: idx[-1] + 2] - a
            return sp.csr_matrix((self.data[a:b].astype(np.float32), self.indices[a:b].astype(np.int32), ip),
                                 shape=(len(idx), self.n_genes))
        starts = self.indptr[idx]; ends = self.indptr[idx + 1]; lens = ends - starts
        ip = np.r_[0, np.cumsum(lens)]
        sel = np.repeat(starts - ip[:-1], lens) + np.arange(ip[-1])
        return sp.csr_matrix((self.data[sel].astype(np.float32), self.indices[sel].astype(np.int32), ip),
                             shape=(len(idx), self.n_genes))

    def chunks(self, cells=None, size=15000):
        cells = np.arange(self.n_cells) if cells is None else np.sort(np.asarray(cells))
        for i in range(0, len(cells), size):
            c = cells[i:i + size]
            yield c, self.rows(c)


def write_store(path, csr_iter, var, obs, n_genes):
    """csr_iter yields CSR count blocks (rows in obs order)."""
    os.makedirs(path, exist_ok=True)
    dfile = open(f"{path}/data.bin", "wb"); ifile = open(f"{path}/indices.bin", "wb")
    indptr = [0]; ncount = []; nnz = 0
    for X in csr_iter:
        X = X.tocsr(); X.sort_indices()
        assert X.data.max(initial=0) < 65535
        dfile.write(np.rint(X.data).astype(np.uint16).tobytes()); ifile.write(X.indices.astype(np.uint16).tobytes())
        indptr.extend((np.asarray(X.indptr[1:], dtype=np.int64) + nnz).tolist()); nnz += X.nnz
        ncount.append(np.asarray(X.sum(1)).ravel())
    dfile.close(); ifile.close()
    for name in ("data", "indices"):
        raw = np.memmap(f"{path}/{name}.bin", dtype=np.uint16, mode="r")
        out = np.lib.format.open_memmap(f"{path}/{name}.npy", mode="w+", dtype=np.uint16, shape=raw.shape)
        for i in range(0, len(raw), 200_000_000):
            out[i:i + 200_000_000] = raw[i:i + 200_000_000]
        out.flush(); del out, raw; os.remove(f"{path}/{name}.bin")
    np.save(f"{path}/indptr.npy", np.asarray(indptr, dtype=np.int64))
    np.save(f"{path}/ncount.npy", np.concatenate(ncount).astype(np.float32))
    var.to_csv(f"{path}/var.csv", index=False); obs.to_parquet(f"{path}/obs.parquet")
    json.dump({"n_cells": int(len(indptr) - 1), "n_genes": int(n_genes), "nnz": int(nnz)}, open(f"{path}/meta.json", "w"))


def lognorm(X, ncount, target=1e4):
    Y = sp.diags((target / np.maximum(ncount, 1)).astype(np.float32)) @ X
    Y = Y.tocsr(); Y.data = np.log1p(Y.data); return Y


# ---------------------------------------------------------------- HVG selection
def vst_hvg(st, cells=None, n_top=2000, min_cells=3, span=0.3):
    """Seurat v4 FindVariableFeatures(selection.method='vst') on raw counts."""
    from skmisc.loess import loess
    G = st.n_genes; s1 = np.zeros(G); s2 = np.zeros(G); nnz = np.zeros(G); N = 0
    for c, X in st.chunks(cells):
        s1 += np.asarray(X.sum(0)).ravel(); s2 += np.bincount(X.indices, weights=X.data.astype(np.float64) ** 2, minlength=G)
        nnz += np.bincount(X.indices, minlength=G); N += X.shape[0]
    mean = s1 / N; var = (s2 - N * mean ** 2) / (N - 1)
    ok = (nnz >= min_cells) & (var > 0)
    lo = loess(np.log10(mean[ok]), np.log10(var[ok]), span=span, degree=2); lo.fit()
    sd = np.zeros(G); sd[ok] = np.sqrt(10 ** lo.outputs.fitted_values); vmax = np.sqrt(N)
    zz = np.zeros(G)
    for c, X in st.chunks(cells):
        g = X.indices; good = sd[g] > 0
        z = np.zeros(len(g)); z[good] = np.minimum((X.data[good] - mean[g][good]) / sd[g][good], vmax)
        zz += np.bincount(g, weights=z ** 2, minlength=G)
    with np.errstate(divide="ignore", invalid="ignore"):
        zz += (N - nnz) * np.where(sd > 0, (mean / sd) ** 2, 0)
        vstd = np.where(ok, zz / (N - 1), 0)
    hvg = np.sort(np.argsort(-vstd)[:n_top])
    return hvg, pd.DataFrame({"mean": mean, "variance": var, "variance_standardized": vstd, "n_cells": nnz})


def seurat_disp_hvg(st, cells=None, min_mean=0.0125, max_mean=3, min_disp=0.5, n_bins=20, exclude=None):
    """scanpy.pp.highly_variable_genes(flavor='seurat') on log1p(normalize_total 1e4) data."""
    G = st.n_genes; s1 = np.zeros(G); s2 = np.zeros(G); N = 0
    for c, X in st.chunks(cells):
        Y = sp.diags((1e4 / np.maximum(st.ncount[c], 1)).astype(np.float32)) @ X; Y = Y.tocsr()
        s1 += np.asarray(Y.sum(0)).ravel(); s2 += np.bincount(Y.indices, weights=Y.data.astype(np.float64) ** 2, minlength=G); N += Y.shape[0]
    mean = s1 / N; var = (s2 - N * mean ** 2) / (N - 1)
    mean[mean == 0] = 1e-12
    disp = var / mean; disp[disp == 0] = np.nan
    disp = np.log(disp); mean_l = np.log1p(mean)
    df = pd.DataFrame({"means": mean_l, "dispersions": disp})
    df["bin"] = pd.cut(df.means, bins=n_bins)
    g = df.groupby("bin", observed=False)["dispersions"]
    m = g.transform("mean"); s = g.transform("std")
    one = g.transform("count") == 1
    s = s.where(~one, m); m = m.where(~one, 0)   # scanpy single-gene-bin convention
    dn = (df.dispersions - m) / s
    hv = (df.means > min_mean) & (df.means < max_mean) & (dn > min_disp)
    if exclude is not None:
        hv &= ~exclude
    return np.where(hv.fillna(False).to_numpy())[0], df.assign(dispersions_norm=dn, highly_variable=hv)


def lognorm_cols(st, cols, cells=None, target=1e4):
    blocks = []
    for c, X in st.chunks(cells):
        Y = lognorm(X[:, cols], st.ncount[c], target); blocks.append(Y.astype(np.float32))
    return sp.vstack(blocks, format="csr")


# ---------------------------------------------------------------- scaling + PCA
def regress_scale_pca(X, covars=None, n_pcs=50, seed=2021, center_only_regress=False):
    """PCA of the scaled matrix A = (X - C B) / s without forming it (C = [1, covars]).
    covars=None -> plain centering + unit variance (scanpy.pp.scale without clip)."""
    X = X.tocsr(); N, G = X.shape
    C = np.ones((N, 1)) if covars is None else np.column_stack([np.ones(N), covars])
    if C.shape[1] > 1:
        C[:, 1:] = (C[:, 1:] - C[:, 1:].mean(0)) / C[:, 1:].std(0)
    CtX = np.asarray((X.T @ C).T); B = np.linalg.solve(C.T @ C, CtX)
    ss = np.bincount(X.indices, weights=X.data.astype(np.float64) ** 2, minlength=G)
    rss = ss - np.einsum("kg,kg->g", CtX, B)
    inv_s = 1.0 / np.sqrt(np.maximum(rss, 1e-12) / (N - 1))
    def mv(v):
        v = np.asarray(v).reshape(G, -1) * inv_s[:, None]; return X @ v - C @ (B @ v)
    def rmv(u):
        u = np.asarray(u).reshape(N, -1); return (np.asarray(X.T @ u) - B.T @ (C.T @ u)) * inv_s[:, None]
    A = LinearOperator((N, G), matvec=mv, rmatvec=rmv, matmat=mv, rmatmat=rmv, dtype=np.float64)
    U, S, Vt = svds(A, k=n_pcs, v0=np.random.default_rng(seed).standard_normal(min(N, G)))
    o = np.argsort(-S); U, S, Vt = U[:, o], S[o], Vt[o]
    sgn = np.sign(Vt[np.arange(n_pcs), np.abs(Vt).argmax(1)]); U *= sgn; Vt *= sgn[:, None]
    return (U * S).astype(np.float32), Vt.T.astype(np.float32), (S ** 2 / (N - 1)).astype(np.float32)


def run_harmony(pcs, batch, seed=2021, max_iter=20):
    import harmonypy
    meta = pd.DataFrame({"batch": pd.Categorical(np.asarray(batch).astype(str))})
    ho = harmonypy.run_harmony(np.asarray(pcs, dtype=np.float64), meta, ["batch"], max_iter_harmony=max_iter,
                               random_state=seed, verbose=False)
    Z = np.asarray(ho.Z_corr)
    return (Z.T if Z.shape[0] == pcs.shape[1] else Z).astype(np.float32)


# ---------------------------------------------------------------- per-group statistics on all genes
def group_stats(st, cells, labels, target=1e4):
    labels = np.asarray(labels).astype(str); groups = sorted(pd.unique(labels), key=lambda x: (len(x), x))
    gi = {g: i for i, g in enumerate(groups)}; li = np.array([gi[g] for g in labels])
    order = np.argsort(cells); cells = np.asarray(cells)[order]; li = li[order]
    S = np.zeros((len(groups), st.n_genes)); P = np.zeros_like(S); n = np.bincount(li, minlength=len(groups)).astype(float)
    pos = 0
    for c, X in st.chunks(cells):
        idx = li[pos:pos + len(c)]; pos += len(c)
        Y = lognorm(X, st.ncount[c], target)
        M = sp.csr_matrix((np.ones(len(idx)), (idx, np.arange(len(idx)))), shape=(len(groups), len(idx)))
        S += (M @ Y).toarray(); Y.data[:] = 1; P += (M @ Y).toarray()
    return groups, S / n[:, None], P / n[:, None], n


def marker_table(groups, mean, pct, n, gene_names, top=30, min_pct=0.1, lfc=0.25):
    tot = n.sum(); rows = []
    for i, g in enumerate(groups):
        rm = ((mean * n[:, None]).sum(0) - mean[i] * n[i]) / (tot - n[i])
        rp = ((pct * n[:, None]).sum(0) - pct[i] * n[i]) / (tot - n[i])
        l2 = np.log2(np.expm1(mean[i]) + 1) - np.log2(np.expm1(rm) + 1)
        keep = (l2 > lfc) & (np.maximum(pct[i], rp) > min_pct)
        for j in np.argsort(-(np.where(keep, l2, -np.inf)))[:top]:
            if keep[j]:
                rows.append((g, gene_names[j], round(float(l2[j]), 3), round(float(pct[i, j]), 3), round(float(rp[j]), 3)))
    return pd.DataFrame(rows, columns=["cluster", "gene", "avg_log2FC", "pct.1", "pct.2"])
