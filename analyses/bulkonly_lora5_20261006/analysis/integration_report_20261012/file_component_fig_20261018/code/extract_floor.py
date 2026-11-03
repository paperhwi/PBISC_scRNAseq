"""Floor read-out (depth256-K4000-floor1e-4, draw0) extract in the same schema as extract_file_component.py (2026-10-22).
Usage: python extract_floor.py <sepsis|covid> <out_dir>. The floor file has one draw, so 20,000 random cells are taken from draw0
(= the 4,000 x 5 draws of the K files) and stored under 'draw0'. It has no aifi labels: aifi_l2_sample = celltype_l1_pred."""
import h5py, numpy as np, sys, os, time
ROOT = os.environ.get("PBISC_ROOT", ".")
P = {"sepsis": "new_inference_depth256_K4000_floor1e-4_20261021/sepsis/PBISC_ye2026_sepsis_draw0_depth256_K4000_floor0p0001.h5ad",
     "covid": "new_inference_depth256_K4000_floor1e-4_20261021/covid/PBISC_imperial_eskd_draw0_depth256_K4000_floor0p0001.h5ad"}
GENES = "CD3E IL7R CCR7 CD4 CD8A NKG7 GNLY MS4A1 CD79A JCHAIN CD14 LYZ S100A8 FCGR3A FCER1A CLEC4C PPBP HBB MKI67".split()
N = 20000; coh, out = sys.argv[1], sys.argv[2]; t0 = time.time()
def col(f, grp, c):
    g = f[grp][c]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]); return cats, g["codes"][:]
    v = np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]]); cats, codes = np.unique(v, return_inverse=True); return cats, codes
with h5py.File(os.path.join(ROOT, P[coh]), "r") as f:
    sc, scode = col(f, "var", "gene_symbol"); sym = sc[scode]; meas = f["var"]["measured_in_cohort"][:]
    gidx = np.array([np.where(sym == g)[0][0] for g in GENES])
    ip = f["X"]["indptr"][:]; n = len(ip) - 1; rows = np.sort(np.random.default_rng(20261022).choice(n, N, replace=False))
    D = f["X"]["data"]; I = f["X"]["indices"]
    tot = np.zeros(N, np.int64); ng = np.zeros(N, np.int32); M = np.zeros((N, len(GENES)), np.int32)
    gsum = np.zeros(len(sym)); gdet = np.zeros(len(sym)); gcnt = np.zeros(len(sym))
    # read in contiguous blocks of rows for speed
    for b0 in range(0, N, 500):
        rr = rows[b0:b0 + 500]
        for k, r in zip(range(b0, b0 + len(rr)), rr):
            d = D[ip[r]:ip[r + 1]].astype(np.int64); ix = I[ip[r]:ip[r + 1]]
            t = d.sum(); tot[k] = t; ng[k] = len(d)
            gsum[ix] += np.log1p(d / max(t, 1) * 1e4); gdet[ix] += 1; gcnt[ix] += d
            h = np.searchsorted(ix, gidx); ok = (h < len(ix)) & (ix[np.minimum(h, len(ix) - 1)] == gidx); M[k, ok] = d[h[ok]]
        if b0 % 5000 == 0: print(coh, b0, round(time.time() - t0), "s", flush=True)
    labs = {c: col(f, "obs", c) for c in ("celltype_l1_pred", "group", "sample_id")}
res = {"rows": rows, "total": tot, "ngenes": ng, "marker_counts": M, "gene_logcp10k_sum": gsum, "gene_det": gdet, "gene_count_sum": gcnt, "n_cells": n,
       "n_samples": len(labs["sample_id"][0])}
for c, (cats, codes) in labs.items(): res[f"{c}_sample"] = cats[codes[rows]]
res["aifi_l2_sample"] = res["celltype_l1_pred_sample"]
for c in ("celltype_l1_pred", "group"): res[f"{c}_allcounts"] = np.bincount(labs[c][1], minlength=len(labs[c][0])); res[f"{c}_cats"] = labs[c][0]
flat = {f"draw0__{k}": np.asarray(v) for k, v in res.items()}; flat["gene_symbol"] = sym; flat["measured"] = meas
np.savez_compressed(os.path.join(out, f"{coh}_floor.npz"), **flat); print("saved", coh, round(time.time() - t0), "s")
