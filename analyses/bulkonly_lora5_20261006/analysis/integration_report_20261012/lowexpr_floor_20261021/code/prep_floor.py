"""Low-expression cliff, floor read-out (2026-10-21). Per-gene statistics of the new inference-style file
(depth256-K4000-floor1e-4, draw0) with the same real sample / panel / depth-matching as lowexpr_cliff_20261020/code/prep.py.
Usage: python prep_floor.py <sepsis|covid> [n_cells=30000]
Sources written to data/<coh>_floor_gene_stats.h5ad:
  realU   real held-out cells (identical sample to the 2026-10-20 run: same seed, first rng draw)
  realT1  real, binomially thinned to the depth of FLOOR cells of the same six-type (depth only)
  floor   new generated cells (X = floor-integerised counts), 30,000 random kept cells
Sepsis: adults only (the file has 1 paediatric ICU-control sample). COVID: the 5 adapter-training samples
(donor_trained_draw0 == True) are excluded.
"""
import sys, os, time, numpy as np, pandas as pd, h5py, scipy.sparse as sp, anndata as ad
coh = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 30000
BASE = os.path.expanduser("~/mnt/PBISC"); AN = f"{BASE}/bulkonly_lora5_20261006/analysis"
IR = f"{AN}/integration_report_20261012"; PREV = f"{IR}/lowexpr_cliff_20261020"; OUT = f"{IR}/lowexpr_floor_20261021"
REAL = f"{AN}/cell_transcriptome_validation_20261014/inputs/" + {"sepsis": "ye2026_sepsis", "covid": "imperial_eskd"}[coh] + "_real_heldout_panel_counts.h5ad"
FL = f"{BASE}/new_inference_depth256_K4000_floor1e-4_20261021/" + {"sepsis": "sepsis/PBISC_ye2026_sepsis", "covid": "covid/PBISC_imperial_eskd"}[coh] + "_draw0_depth256_K4000_floor0p0001.h5ad"
T4 = f"{IR}/file_component_fig_20261018/results/t4_real_label_to_six_type.csv"
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

def col(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object)
        return cats[g["codes"][:]]
    v = g[:]; return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v

def scan(f, key, keep_rows, take_rows, cols, ncol_out, block=20000):
    g = f[key]; ip = g["indptr"][:]; n = len(ip) - 1; ncol = int(f[key].attrs["shape"][1])
    det = np.zeros(ncol_out, np.int64); tot = np.zeros(ncol_out, np.float64); parts = []
    take = np.zeros(n, bool); take[take_rows] = True
    for s in range(0, n, block):
        e = min(n, s + block); kr = keep_rows[s:e]
        if not kr.any(): continue
        a, b = ip[s], ip[e]
        M = sp.csr_matrix((g["data"][a:b].astype(np.float32), g["indices"][a:b], ip[s:e + 1] - a), shape=(e - s, ncol))
        M = M[:, np.where(cols >= 0)[0]]; Mk = M[np.where(kr)[0]]
        det += np.bincount(Mk.indices, minlength=ncol_out); tot += np.bincount(Mk.indices, weights=Mk.data, minlength=ncol_out)
        tr = np.where(take[s:e])[0]
        if len(tr): parts.append(M[tr])
    return sp.vstack(parts).tocsr(), det, tot

def first(sym): s = pd.Series(np.arange(len(sym)), index=sym); return s[~s.index.duplicated()]

# ---- panel = the 2026-10-20 panel (stored var_names), restricted to genes present in the floor file
P0 = ad.read_h5ad(f"{PREV}/data/{coh}_gene_stats.h5ad"); genes0 = np.array(P0.var_names)
fr = h5py.File(REAL, "r"); ff = h5py.File(FL, "r")
sym_r = col(fr, "var/gene_symbol").astype(str); sym_f = col(ff, "var/gene_symbol").astype(str)
pr, pf = first(sym_r), first(sym_f)
genes = np.array([g for g in genes0 if g in pf.index]); miss = [g for g in genes0 if g not in pf.index]
log(coh, "previous panel", len(genes0), "| in floor file", len(genes), "| missing", miss[:10])
ord_r = pr.loc[genes].to_numpy(); assert np.all(np.diff(ord_r) > 0), "real panel order"
cm_r = np.full(len(sym_r), -1); cm_r[ord_r] = np.arange(len(ord_r))

# ---- real (same sample as 2026-10-20: np.random.default_rng(20261020), first draw)
rng = np.random.default_rng(20261020)
grp = col(fr, "obs/group").astype(str); major = col(fr, "obs/major").astype(str)
if coh == "sepsis":
    m6 = {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell",
          "NK cells": "NK cell", "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC", "Other cells": "Other"}
    six_r = pd.Series(major).map(m6).to_numpy(); keep_r = pd.Series(grp).str.startswith("adult").to_numpy()
else:
    t4 = pd.read_csv(T4); t4 = t4[t4.cohort == "covid"]; six_r = pd.Series(major).map(dict(zip(t4.real_label, t4.six_type))).to_numpy()
    keep_r = np.ones(len(grp), bool)
six_r = np.where(pd.isna(six_r), "Other", six_r).astype(str)
take_r = np.sort(rng.choice(np.where(keep_r)[0], size=min(N, keep_r.sum()), replace=False))
Xr, det_r, tot_r = scan(fr, "X", keep_r, take_r, cm_r, len(genes)); lab_r = six_r[take_r]
log("real", Xr.shape, Xr.nnz, "scope cells", int(keep_r.sum()))
# consistency with the previous run
iu = np.where((P0.obs.source == "realU") & (P0.obs.level == "all"))[0][0]
prev_det = pd.Series(P0.layers["det"][iu], index=P0.var_names).loc[genes].to_numpy(); mine = np.diff(Xr.tocsc().indptr) / Xr.shape[0]
log("realU sample reproduces previous run: max |det diff| =", float(np.abs(prev_det - mine).max()))

# ---- floor file
pos = pf.loc[genes].to_numpy(); ordk = np.argsort(pos); cm = np.full(len(sym_f), -1); cm[pos[ordk]] = np.arange(len(pos))
gk = col(ff, "obs/group").astype(str)
if coh == "sepsis": kk = ~np.char.startswith(gk, "ped")
else: kk = ~col(ff, "obs/donor_trained_draw0").astype(bool)
rng2 = np.random.default_rng(20261021)
tk = np.sort(rng2.choice(np.where(kk)[0], size=min(N, kk.sum()), replace=False))
M, dF, tF = scan(ff, "X", kk, tk, cm, len(genes)); inv = np.empty_like(ordk); inv[ordk] = np.arange(len(ordk))
Xf = M[:, inv].tocsr(); dF, tF = dF[inv], tF[inv]; lab_f = col(ff, "obs/celltype_l1_pred").astype(str)[tk]
log("floor", Xf.shape, Xf.nnz, "scope cells", int(kk.sum()))
full = pd.DataFrame({"gene": genes, "measured_real": fr["var/measured_in_cohort"][:].astype(bool)[ord_r],
                     "measured_floor": ff["var/measured_in_cohort"][:].astype(bool)[pos], "real_n_det": det_r, "real_sum": tot_r,
                     "floor_n_det": dF, "floor_sum": tF, "real_n_cells": int(keep_r.sum()), "floor_n_cells": int(kk.sum())})
full.to_csv(f"{OUT}/results/{coh}_floor_fullfile_gene_totals.csv", index=False)

# ---- depth-matched real (T1 against FLOOR depth, per six-type)
tot = np.asarray(Xr.sum(1)).ravel(); gt = np.asarray(Xf.sum(1)).ravel()
pool = {c: gt[lab_f == c] for c in np.unique(lab_f)}
target = np.array([rng2.choice(pool.get(c, gt)) for c in lab_r], np.float64)
keep = np.minimum(1.0, target / np.maximum(tot, 1))
T1 = Xr.copy(); T1.data = rng2.binomial(T1.data.astype(np.int64), np.repeat(keep, np.diff(T1.indptr))).astype(np.float32); T1.eliminate_zeros()

SRC = {"realU": (Xr, lab_r), "realT1": (T1, lab_r), "floor": (Xf, lab_f)}
rows, Xm, Ld, Lu, dep = [], [], [], [], []
def add(X, src, level, label):
    t = np.asarray(X.sum(1)).ravel(); C = (sp.diags(1e4 / np.maximum(t, 1)) @ X).tocsr()
    Xm.append(np.asarray(C.mean(0)).ravel()); Lu.append(np.asarray(X.mean(0)).ravel())
    Ld.append(np.diff(X.tocsc().indptr) / X.shape[0]); rows.append(dict(source=src, level=level, label=label, n_cells=X.shape[0]))
for s, (X, lab) in SRC.items():
    add(X, s, "all", "all"); t = np.asarray(X.sum(1)).ravel(); d = np.diff(X.indptr)
    dep.append(dict(source=s, umi_median=float(np.median(t)), umi_mean=float(t.mean()), detected_median=float(np.median(d)), frac_ones=float((X.data == 1).mean())))
    for c in np.unique(lab):
        m = lab == c
        if m.sum() >= 100: add(X[np.where(m)[0]], s, "six_type", c)
obs = pd.DataFrame(rows); obs.index = [f"{a}|{b}|{c}" for a, b, c in zip(obs.source, obs.level, obs.label)]
var = pd.DataFrame({"measured_real": full.measured_real.to_numpy(), "measured_floor": full.measured_floor.to_numpy()}, index=pd.Index(genes, name="gene"))
A_ = ad.AnnData(np.vstack(Xm).astype(np.float32), obs=obs, var=var, layers={"det": np.vstack(Ld).astype(np.float32), "mean_umi": np.vstack(Lu).astype(np.float32)})
A_.uns["description"] = ("Per-gene statistics of the floor read-out file (draw0, depth256, K4000, floor 1e-4) and matched real. X = mean CP10k, "
                         "layers det = fraction of cells with >=1 UMI, mean_umi = mean raw UMI per cell. realT1 is depth-matched to the FLOOR cells. See code/prep_floor.py.")
A_.uns["n_sampled"] = N; A_.uns["seed_real"] = 20261020; A_.uns["seed_floor"] = 20261021; A_.uns["cohort"] = coh
A_.write_h5ad(f"{OUT}/data/{coh}_floor_gene_stats.h5ad", compression="gzip")
pd.DataFrame(dep).to_csv(f"{OUT}/results/{coh}_floor_depth.csv", index=False)
log("done", A_.shape)
