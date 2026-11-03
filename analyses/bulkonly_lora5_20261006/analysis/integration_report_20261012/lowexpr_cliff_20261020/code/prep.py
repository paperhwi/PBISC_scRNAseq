"""Low-expression cliff: analyses 1 and 2 (2026-10-20). Per-gene statistics for one cohort.

Usage: python prep.py <sepsis|covid> [n_cells=30000]

Sources (all on the same gene panel = generated K1000 ∩ K4000 ∩ real held-out panel, first occurrence of duplicated symbols):
  realU     real held-out cells (cell_transcriptome_validation_20261014), untouched
  realT1    depth-matched binomial thinning to a UMI total drawn from K=1,000 generated cells of the same six-type
  realT2    PBISC-like read-out, identical to integration_report e_checks.py: expected counts = kNN(k=15, PCA30)-smoothed
            profile restricted to the cell's OWN detected genes, scaled to the matched depth, one Poisson draw
  realT2b   same, but support re-drawn as Bernoulli(p), p = share of the 15 neighbours detecting the gene (calibrated support)
  realT2m   same, but support = genes with p >= 0.5 (mode support: what a deterministic, mode-seeking support channel returns)
  k1000/k4000  generated cells (Harmony-set h5ad, counts layer, 5 draws)
Sepsis: adults only on both sides (generated set has one paediatric ICU-control sample, removed).

Outputs (lowexpr_cliff_20261020/): data/<coh>_gene_stats.h5ad (obs = source x group, var = genes; X = mean CP10k,
layers det / mean_umi), results/<coh>_fullfile_gene_totals.csv (all scoped cells, not a sample), results/<coh>_depth.csv
"""
import sys, json, time, numpy as np, pandas as pd, h5py, scipy.sparse as sp, anndata as ad
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
coh = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
A = "/sessions/" if False else None
import os
BASE = os.path.expanduser("~/mnt/PBISC")
AN = f"{BASE}/bulkonly_lora5_20261006/analysis"
OUT = f"{AN}/integration_report_20261012/lowexpr_cliff_20261020"
REAL = {"sepsis": f"{AN}/cell_transcriptome_validation_20261014/inputs/ye2026_sepsis_real_heldout_panel_counts.h5ad",
        "covid": f"{AN}/cell_transcriptome_validation_20261014/inputs/imperial_eskd_real_heldout_panel_counts.h5ad"}[coh]
GEN = {"k1000": f"{AN}/harmony/{coh}_k1000/{coh}_k1000_harmony_{'clustering' if coh == 'sepsis' else 'embedding'}.h5ad",
       "k4000": f"{BASE}/bulkonly_lora5_k4000_20261006/analysis/harmony/{coh}_k4000/{coh}_k4000_harmony_{'clustering' if coh == 'sepsis' else 'embedding'}.h5ad"}
T4 = f"{AN}/integration_report_20261012/file_component_fig_20261018/results/t4_real_label_to_six_type.csv"
rng = np.random.default_rng(20261020); t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)


def col(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object)
        return cats[g["codes"][:]]
    v = g[:]; return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v


def scan(f, key, keep_rows, take_rows, cols, ncol_out, block=20000):
    """One pass over a CSR matrix: returns (sample matrix of take_rows restricted to cols, full-scope per-gene nnz and sum
    over keep_rows). cols maps file column -> output column (-1 = drop)."""
    g = f[key]; ip = g["indptr"][:]; n = len(ip) - 1; ncol = int(f[key].attrs["shape"][1])
    det = np.zeros(ncol_out, np.int64); tot = np.zeros(ncol_out, np.float64); parts = []
    take = np.zeros(n, bool); take[take_rows] = True
    for s in range(0, n, block):
        e = min(n, s + block); kr = keep_rows[s:e]
        if not kr.any(): continue
        a, b = ip[s], ip[e]
        M = sp.csr_matrix((g["data"][a:b].astype(np.float32), g["indices"][a:b], ip[s:e + 1] - a), shape=(e - s, ncol))
        M = M[:, np.where(cols >= 0)[0]]                      # panel columns, in output order (cols sorted below)
        Mk = M[np.where(kr)[0]]
        det += np.bincount(Mk.indices, minlength=ncol_out); tot += np.bincount(Mk.indices, weights=Mk.data, minlength=ncol_out)
        tr = np.where(take[s:e])[0]
        if len(tr): parts.append(M[tr])
    return sp.vstack(parts).tocsr(), det, tot


# ---------------- gene panel
fr = h5py.File(REAL, "r"); fg = {K: h5py.File(p, "r") for K, p in GEN.items()}
sym_r = col(fr, "var/gene_symbol").astype(str); sym_g = {K: col(f, "var/gene_symbol").astype(str) for K, f in fg.items()}
meas_g = fg["k1000"]["var/measured_in_cohort"][:].astype(bool); meas_r = fr["var/measured_in_cohort"][:].astype(bool)
def first(sym): s = pd.Series(np.arange(len(sym)), index=sym); return s[~s.index.duplicated()]
pr, p1, p4 = first(sym_r), first(sym_g["k1000"]), first(sym_g["k4000"])
panel = [g for g in p1.index if g in p4.index and g in pr.index]
log(coh, "panel", len(panel), "| gen", len(sym_g["k1000"]), "real", len(sym_r))
def colmap(pos, nfile):
    """file column -> output position. Output order must equal the order of kept file columns, so sort panel by file column."""
    m = np.full(nfile, -1); m[pos.loc[panel].to_numpy()] = 1; return m
# make output order = order of real file columns; for generated files remap after reading
ord_r = np.sort(pr.loc[panel].to_numpy()); genes = sym_r[ord_r]
cm_r = np.full(len(sym_r), -1); cm_r[ord_r] = np.arange(len(ord_r))

# ---------------- real
grp = col(fr, "obs/group").astype(str); major = col(fr, "obs/major").astype(str); fine = col(fr, "obs/fine").astype(str)
if coh == "sepsis":
    m6 = {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell",
          "NK cells": "NK cell", "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC", "Other cells": "Other"}
    six_r = pd.Series(major).map(m6).to_numpy(); keep_r = pd.Series(grp).str.startswith("adult").to_numpy()
else:
    t4 = pd.read_csv(T4); t4 = t4[t4.cohort == "covid"]; six_r = pd.Series(major).map(dict(zip(t4.real_label, t4.six_type))).to_numpy()
    keep_r = np.ones(len(grp), bool)
log("real unmapped labels:", sorted(set(major[pd.isna(six_r)])))
six_r = np.where(pd.isna(six_r), "Other", six_r).astype(str)
take_r = np.sort(rng.choice(np.where(keep_r)[0], size=min(N, keep_r.sum()), replace=False))
Xr, det_r, tot_r = scan(fr, "X", keep_r, take_r, cm_r, len(genes))
lab_r, fine_r = six_r[take_r], fine[take_r]
log("real read", Xr.shape, Xr.nnz, "scope cells", int(keep_r.sum()))
full = pd.DataFrame({"gene": genes, "measured_real": meas_r[ord_r], "real_n_det": det_r, "real_sum": tot_r}); full_n = {"real": int(keep_r.sum())}

# ---------------- generated
G, lab_g, draw_g = {}, {}, {}
for K, f in fg.items():
    pos = (p1 if K == "k1000" else p4).loc[genes].to_numpy()          # file column of each output gene
    ordk = np.argsort(pos); cm = np.full(len(sym_g[K]), -1); cm[pos[ordk]] = np.arange(len(pos))
    gk = col(f, "obs/group").astype(str); kk = ~np.char.startswith(gk.astype(str), "ped") if coh == "sepsis" else np.ones(len(gk), bool)
    tk = np.sort(rng.choice(np.where(kk)[0], size=min(N, kk.sum()), replace=False))
    M, d, t = scan(f, "layers/counts", kk, tk, cm, len(genes))
    inv = np.empty_like(ordk); inv[ordk] = np.arange(len(ordk))      # scan returns columns in file order -> back to gene order
    G[K] = M[:, inv].tocsr(); full[f"{K}_n_det"] = d[inv]; full[f"{K}_sum"] = t[inv]; full_n[K] = int(kk.sum())
    lab_g[K] = col(f, "obs/celltype_l1_pred").astype(str)[tk]; draw_g[K] = col(f, "obs/draw").astype(str)[tk]
    if K == "k1000": full["measured_gen"] = meas_g[pos]
    log(K, "read", G[K].shape, G[K].nnz, "scope cells", full_n[K])
for k, v in full_n.items(): full[f"{k}_n_cells"] = v
full.to_csv(f"{OUT}/results/{coh}_fullfile_gene_totals.csv", index=False)

# ---------------- transforms of the real sample
tot = np.asarray(Xr.sum(1)).ravel(); gt = np.asarray(G["k1000"].sum(1)).ravel()
pool = {c: gt[lab_g["k1000"] == c] for c in np.unique(lab_g["k1000"])}
target = np.array([rng.choice(pool.get(c, gt)) for c in lab_r], np.float64)
keep = np.minimum(1.0, target / np.maximum(tot, 1))
T1 = Xr.copy(); T1.data = rng.binomial(T1.data.astype(np.int64), np.repeat(keep, np.diff(T1.indptr))).astype(np.float32); T1.eliminate_zeros()
def lognorm(X):
    t = np.asarray(X.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(t, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y
Y = lognorm(Xr); mu = np.asarray(Y.mean(0)).ravel(); var = np.asarray(Y.multiply(Y).mean(0)).ravel() - mu ** 2
detf = np.diff(Xr.tocsc().indptr) / Xr.shape[0]
hv = np.argsort(-np.where(detf >= 0.01, var, -1))[:2000]
Z = np.clip((Y[:, hv].toarray() - mu[hv]) / (np.sqrt(var[hv]) + 1e-6), -10, 10)
P = PCA(30, random_state=0).fit_transform(Z).astype(np.float32); del Z
_, nb = NearestNeighbors(n_neighbors=15).fit(P).kneighbors(P)
n = Xr.shape[0]
W = sp.csr_matrix((np.full(nb.size, 1 / 15, np.float32), nb.ravel(), np.arange(0, nb.size + 1, 15)), shape=(n, n))
prof = (sp.diags(1 / np.maximum(tot, 1)) @ Xr).tocsr(); B = Xr.copy(); B.data[:] = 1
def readout(E, rows):
    rs = np.asarray(E.sum(1)).ravel(); E = (sp.diags(target[rows] / np.maximum(rs, 1e-12)) @ E).tocsr()
    E.data = rng.poisson(E.data).astype(np.float32); E.eliminate_zeros(); return E
T2, T2b, T2m, pmean = [], [], [], []
for i in range(0, n, 2000):
    r = np.arange(i, min(n, i + 2000)); S = (W[r] @ prof).tocsr(); Pd = (W[r] @ B).tocsr()   # same sparsity pattern (union support)
    S.sort_indices(); Pd.sort_indices()
    T2.append(readout(S.multiply(B[r]).tocsr(), r))
    u = rng.random(Pd.nnz).astype(np.float32); mb = Pd.copy(); mb.data = (u < Pd.data).astype(np.float32); mb.eliminate_zeros()
    T2b.append(readout(S.multiply(mb).tocsr(), r))
    mm = Pd.copy(); mm.data = (Pd.data >= 0.5 - 1e-6).astype(np.float32); mm.eliminate_zeros()
    T2m.append(readout(S.multiply(mm).tocsr(), r))
T2, T2b, T2m = [sp.vstack(x).tocsr() for x in (T2, T2b, T2m)]
log("transforms done")

# ---------------- per-gene statistics
SRC = {"realU": (Xr, lab_r), "realT1": (T1, lab_r), "realT2": (T2, lab_r), "realT2b": (T2b, lab_r), "realT2m": (T2m, lab_r),
       "k1000": (G["k1000"], lab_g["k1000"]), "k4000": (G["k4000"], lab_g["k4000"])}
rows, Xm, Ld, Lu, dep = [], [], [], [], []
def add(X, src, level, label):
    t = np.asarray(X.sum(1)).ravel(); C = (sp.diags(1e4 / np.maximum(t, 1)) @ X).tocsr()
    Xm.append(np.asarray(C.mean(0)).ravel()); Lu.append(np.asarray(X.mean(0)).ravel())
    Ld.append(np.diff(X.tocsc().indptr) / X.shape[0]); rows.append(dict(source=src, level=level, label=label, n_cells=X.shape[0]))
for s, (X, lab) in SRC.items():
    add(X, s, "all", "all")
    t = np.asarray(X.sum(1)).ravel(); d = np.diff(X.indptr)
    dep.append(dict(source=s, umi_median=float(np.median(t)), umi_mean=float(t.mean()), detected_median=float(np.median(d)), frac_ones=float((X.data == 1).mean())))
    for c in np.unique(lab):
        m = lab == c
        if m.sum() >= 100: add(X[np.where(m)[0]], s, "six_type", c)
for c, nc in pd.Series(fine_r).value_counts().items():
    if nc >= 100: add(Xr[np.where(fine_r == c)[0]], "realU", "fine", c)
for K in G:
    for d in np.unique(draw_g[K]): add(G[K][np.where(draw_g[K] == d)[0]], K, "draw", d)
obs = pd.DataFrame(rows); obs.index = [f"{a}|{b}|{c}" for a, b, c in zip(obs.source, obs.level, obs.label)]
var = pd.DataFrame({"measured_real": meas_r[ord_r], "measured_gen": full.measured_gen.to_numpy()}, index=pd.Index(genes, name="gene"))
A_ = ad.AnnData(np.vstack(Xm).astype(np.float32), obs=obs, var=var, layers={"det": np.vstack(Ld).astype(np.float32), "mean_umi": np.vstack(Lu).astype(np.float32)})
A_.uns["description"] = ("Per-gene statistics. X = mean CP10k (each cell normalised to 10,000 on its panel total), layers['det'] = "
                         "fraction of cells with >=1 UMI, layers['mean_umi'] = mean raw UMI per cell. obs = source x level x label; "
                         "level 'all', 'six_type' (real: author major label mapped to the six PBISC types; generated: frozen scANVI "
                         "celltype_l1_pred), 'fine' (real author fine labels, n>=100), 'draw' (generated draws). See code/prep.py.")
A_.uns["n_sampled"] = N; A_.uns["seed"] = 20261020; A_.uns["cohort"] = coh
A_.write_h5ad(f"{OUT}/data/{coh}_gene_stats.h5ad", compression="gzip")
pd.DataFrame(dep).to_csv(f"{OUT}/results/{coh}_depth.csv", index=False)
log("done", A_.shape)
