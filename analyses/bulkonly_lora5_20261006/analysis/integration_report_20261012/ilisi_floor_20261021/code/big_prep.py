"""Larger bundles for UMAP / cluster matching / rare-population checks (2026-10-21).
Usage: python big_prep.py <sepsis|covid> <ref|realU|floor|prev4k> [n=60000]
Same patient split as ilisi_prep.py (rng seed 20261021, first draw = patient shuffle); cells are drawn with seed 20261022
(+ source index), equal per sample. Output data/big/<coh>_<src>.h5ad: raw counts on the lowexpr panel,
obs sample, patient, group, six_type, lineage, fine, major (real: author labels; generated: frozen scANVI l1 label, fine = '')."""
import sys, os, time, numpy as np, pandas as pd, h5py, scipy.sparse as sp, anndata as ad
coh, src = sys.argv[1], sys.argv[2]; n = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
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
def read_rows(f, key, rows, colpos, step=20000):
    """contiguous chunked read; rows sorted"""
    g = f[key]; ip = g["indptr"][:]; ncol = int(g.attrs["shape"][1]); out = []
    for s in range(0, int(rows.max()) + 1, step):
        r = rows[(rows >= s) & (rows < s + step)]
        if not len(r): continue
        e = min(s + step, len(ip) - 1); a, b = ip[s], ip[e]
        M = sp.csr_matrix((g["data"][a:b].astype(np.float32), g["indices"][a:b], ip[s:e + 1] - a), shape=(e - s, ncol))
        out.append(M[r - s][:, colpos])
    return sp.vstack(out).tocsr()
def strat(samples, ids, n, rng):
    u = pd.unique(samples[ids]); q = int(np.ceil(n / len(u))); pick = []
    for s in u:
        ix = ids[samples[ids] == s]; pick.append(rng.permutation(ix)[:q])
    pick = np.concatenate(pick); return np.sort(rng.permutation(pick)[:n])
genes = np.array(ad.read_h5ad(f"{OUT}/data/{coh}_ref.h5ad", backed="r").var_names)
LIN = {"T cell": "T", "NK cell": "NK", "B cell": "B", "Monocyte": "Myeloid", "DC": "Myeloid", "Other": "Other"}
if src in ("ref", "realU"):
    f = h5py.File(REAL, "r"); grp = col(f, "obs/group").astype(str); samp = col(f, "obs/sample_id").astype(str); pat = col(f, "obs/patient_id").astype(str)
    major = col(f, "obs/major").astype(str); fine = col(f, "obs/fine").astype(str)
    if coh == "sepsis":
        m6 = {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell", "NK cells": "NK cell",
              "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC", "Other cells": "Other"}
        six = pd.Series(major).map(m6).to_numpy(); keep = pd.Series(grp).str.startswith("adult").to_numpy()
    else:
        t4 = pd.read_csv(T4); t4 = t4[t4.cohort == "covid"]; six = pd.Series(major).map(dict(zip(t4.real_label, t4.six_type))).to_numpy(); keep = np.ones(len(grp), bool)
    six = np.where(pd.isna(six), "Other", six).astype(str)
    pats = np.array(sorted(pd.unique(pat[keep]))); rng.shuffle(pats); half = len(pats) // 2; A, B = set(pats[:half]), set(pats[half:])
    pool = np.where(keep & np.isin(pat, list(A if src == "ref" else B)))[0]
    rng2 = np.random.default_rng(20261022 + (0 if src == "ref" else 1)); ix = strat(samp, pool, min(n, len(pool)), rng2)
    log(coh, src, "pool", len(pool), "take", len(ix), "patients", len(A if src == "ref" else B))
    X = read_rows(f, "X", ix, first(col(f, "var/gene_symbol").astype(str)).loc[genes].to_numpy())
    obs = pd.DataFrame({"sample": samp[ix], "patient": pat[ix], "group": grp[ix], "six_type": six[ix], "fine": fine[ix], "major": major[ix]})
else:
    f = h5py.File(FL if src == "floor" else PREV, "r"); grp = col(f, "obs/group").astype(str)
    sn = "obs/sample_id" if "sample_id" in f["obs"] else "obs/sample"; samp = col(f, sn).astype(str); lab = col(f, "obs/celltype_l1_pred").astype(str)
    if coh == "sepsis": keep = ~np.char.startswith(grp, "ped")
    else: keep = ~col(f, "obs/donor_trained_draw0").astype(bool) if src == "floor" else np.ones(len(grp), bool)
    rng2 = np.random.default_rng(20261022 + (2 if src == "floor" else 3)); ix = strat(samp, np.where(keep)[0], n, rng2)
    log(coh, src, "pool", int(keep.sum()), "take", len(ix), "samples", len(pd.unique(samp[keep])))
    X = read_rows(f, "X" if src == "floor" else "layers/counts", ix, first(col(f, "var/gene_symbol").astype(str)).loc[genes].to_numpy())
    obs = pd.DataFrame({"sample": samp[ix], "patient": samp[ix], "group": grp[ix], "six_type": lab[ix], "fine": "", "major": ""})
obs["lineage"] = [LIN.get(s, "Other") for s in obs.six_type]; obs.index = [f"{src}_{i}" for i in range(len(obs))]
a = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene"))); a.uns["source"] = src; a.uns["cohort"] = coh
a.write_h5ad(f"{OUT}/data/big/{coh}_{src}.h5ad", compression="gzip"); log("saved", X.shape, "mean UMI", round(float(X.sum() / X.shape[0]), 0), obs.group.value_counts().to_dict())
