"""Build balanced input bundles for the real-vs-generated iLISI (2026-10-21). Usage: python ilisi_prep.py <sepsis|covid> [n=20000]
Real held-out cells (cell_transcriptome_validation_20261014) are split BY PATIENT into reference half A and control half B.
Bundles (data/<coh>_<name>.h5ad, counts on the lowexpr panel, n cells each, per-sample-stratified draws):
  ref     real half A (batch partner for every comparison)
  realU   real half B, untouched (upper control)
  realT1  real half B, binomially thinned to the depth of floor cells of the same six-type (depth-matched control)
  floor   new floor read-out file (depth256-K4000-floor1e-4, draw0)
  prev4k  previous K4000 Harmony-set cells (Poisson read-out, 5 draws pooled)
Sepsis: adults only. COVID floor: the 5 adapter-training samples are excluded."""
import sys, os, time, numpy as np, pandas as pd, h5py, scipy.sparse as sp, anndata as ad
coh = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
BASE = os.path.expanduser("~/mnt/PBISC"); AN = f"{BASE}/bulkonly_lora5_20261006/analysis"; IR = f"{AN}/integration_report_20261012"; OUT = f"{IR}/ilisi_floor_20261021"
REAL = f"{AN}/cell_transcriptome_validation_20261014/inputs/" + {"sepsis": "ye2026_sepsis", "covid": "imperial_eskd"}[coh] + "_real_heldout_panel_counts.h5ad"
FL = f"{BASE}/new_inference_depth256_K4000_floor1e-4_20261021/" + {"sepsis": "sepsis/PBISC_ye2026_sepsis", "covid": "covid/PBISC_imperial_eskd"}[coh] + "_draw0_depth256_K4000_floor0p0001.h5ad"
PREV = f"{BASE}/bulkonly_lora5_k4000_20261006/analysis/harmony/{coh}_k4000/{coh}_k4000_harmony_{'clustering' if coh == 'sepsis' else 'embedding'}.h5ad"
T4 = f"{IR}/file_component_fig_20261018/results/t4_real_label_to_six_type.csv"
t0 = time.time(); rng = np.random.default_rng(20261021)
def log(*a): print(f"[{time.time()-t0:5.0f}s]", *a, flush=True)
def col(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object); return cats[g["codes"][:]]
    v = g[:]; return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v
def first(sym): s = pd.Series(np.arange(len(sym)), index=sym); return s[~s.index.duplicated()]
def read_rows(f, key, rows, colpos):
    """rows (sorted file row indices) x columns colpos (file columns in output order) -> csr float32"""
    g = f[key]; ip = g["indptr"][:]; ncol = int(g.attrs["shape"][1]); out = []
    inv = np.full(ncol, -1); inv[colpos] = np.arange(len(colpos)); cm = inv
    for s in range(0, len(rows), 4000):
        r = rows[s:s + 4000]; blocks = []
        for i in r:
            a, b = ip[i], ip[i + 1]; blocks.append((g["indices"][a:b], g["data"][a:b]))
        idx = np.concatenate([x[0] for x in blocks]); dat = np.concatenate([x[1] for x in blocks]).astype(np.float32); lens = [len(x[0]) for x in blocks]
        M = sp.csr_matrix((dat, idx, np.r_[0, np.cumsum(lens)]), shape=(len(r), ncol))[:, colpos]; out.append(M)
    return sp.vstack(out).tocsr()
def strat(samples, ids, n, rng):
    """n cells, as equal as possible per sample"""
    u = pd.unique(samples[ids]); q = int(np.ceil(n / len(u))); pick = []
    for s in u:
        ix = ids[samples[ids] == s]; pick.append(rng.permutation(ix)[:q])
    pick = np.concatenate(pick); return np.sort(rng.permutation(pick)[:n])

P0 = ad.read_h5ad(f"{IR}/lowexpr_floor_20261021/data/{coh}_floor_gene_stats.h5ad"); genes = np.array(P0.var_names)
fr = h5py.File(REAL, "r"); ff = h5py.File(FL, "r"); fp = h5py.File(PREV, "r")
sym_r, sym_f, sym_p = [col(f, "var/gene_symbol").astype(str) for f in (fr, ff, fp)]
pos = {k: first(s).loc[genes].to_numpy() for k, s in (("r", sym_r), ("f", sym_f), ("p", sym_p))}
grp = col(fr, "obs/group").astype(str); samp = col(fr, "obs/sample_id").astype(str); pat = col(fr, "obs/patient_id").astype(str); major = col(fr, "obs/major").astype(str)
if coh == "sepsis":
    m6 = {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell", "NK cells": "NK cell",
          "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC", "Other cells": "Other"}
    six_r = pd.Series(major).map(m6).to_numpy(); keep_r = pd.Series(grp).str.startswith("adult").to_numpy()
else:
    t4 = pd.read_csv(T4); t4 = t4[t4.cohort == "covid"]; six_r = pd.Series(major).map(dict(zip(t4.real_label, t4.six_type))).to_numpy(); keep_r = np.ones(len(grp), bool)
six_r = np.where(pd.isna(six_r), "Other", six_r).astype(str)
LIN = {"T cell": "T", "NK cell": "NK", "B cell": "B", "Monocyte": "Myeloid", "DC": "Myeloid", "Other": "Other"}
pats = np.array(sorted(pd.unique(pat[keep_r]))); rng.shuffle(pats); half = len(pats) // 2; A, B = set(pats[:half]), set(pats[half:])
inA = keep_r & np.isin(pat, list(A)); inB = keep_r & np.isin(pat, list(B))
log(coh, "patients", len(pats), "A", len(A), "B", len(B), "| cells A", int(inA.sum()), "B", int(inB.sum()))
iA = strat(samp, np.where(inA)[0], n, rng); iB = strat(samp, np.where(inB)[0], n, rng)
# floor
gk = col(ff, "obs/group").astype(str); kk = ~np.char.startswith(gk, "ped") if coh == "sepsis" else ~col(ff, "obs/donor_trained_draw0").astype(bool)
sf = col(ff, "obs/sample_id").astype(str); iF = strat(sf, np.where(kk)[0], n, rng)
# previous K4000 (Harmony-set)
gp = col(fp, "obs/group").astype(str); kp = ~np.char.startswith(gp, "ped") if coh == "sepsis" else np.ones(len(gp), bool)
spn = "obs/sample_id" if "sample_id" in fp["obs"] else "obs/sample"; sp_ = col(fp, spn).astype(str); iP = strat(sp_, np.where(kp)[0], n, rng)
log("selected", len(iA), len(iB), len(iF), len(iP))
XA = read_rows(fr, "X", iA, pos["r"]); XB = read_rows(fr, "X", iB, pos["r"]); XF = read_rows(ff, "X", iF, pos["f"]); XP = read_rows(fp, "layers/counts", iP, pos["p"])
log("read", XA.shape, XB.shape, XF.shape, XP.shape)
labF = col(ff, "obs/celltype_l1_pred").astype(str)[iF]; labP = col(fp, "obs/celltype_l1_pred").astype(str)[iP]
# depth-matched control (B thinned to floor depth, per six-type)
tot = np.asarray(XB.sum(1)).ravel(); gt = np.asarray(XF.sum(1)).ravel(); pool = {c: gt[labF == c] for c in np.unique(labF)}
sixB = six_r[iB]; target = np.array([rng.choice(pool.get(c, gt)) for c in sixB], np.float64); keep = np.minimum(1.0, target / np.maximum(tot, 1))
T1 = XB.copy(); T1.data = rng.binomial(T1.data.astype(np.int64), np.repeat(keep, np.diff(T1.indptr))).astype(np.float32); T1.eliminate_zeros()
def save(name, X, sample, patient, six):
    o = pd.DataFrame({"sample": sample, "patient": patient, "six_type": six, "lineage": [LIN.get(s, "Other") for s in six]}, index=[f"{name}_{i}" for i in range(X.shape[0])])
    a = ad.AnnData(X.tocsr(), obs=o, var=pd.DataFrame(index=pd.Index(genes, name="gene"))); a.uns["source"] = name; a.uns["cohort"] = coh
    a.write_h5ad(f"{OUT}/data/{coh}_{name}.h5ad", compression="gzip")
save("ref", XA, samp[iA], pat[iA], six_r[iA]); save("realU", XB, samp[iB], pat[iB], six_r[iB]); save("realT1", T1, samp[iB], pat[iB], six_r[iB])
save("floor", XF, sf[iF], sf[iF], labF); save("prev4k", XP, sp_[iP], sp_[iP], labP)
log("saved; depth (mean UMI) ref/realU/realT1/floor/prev4k:", [round(float(np.asarray(x.sum(1)).mean()), 0) for x in (XA, XB, T1, XF, XP)])
