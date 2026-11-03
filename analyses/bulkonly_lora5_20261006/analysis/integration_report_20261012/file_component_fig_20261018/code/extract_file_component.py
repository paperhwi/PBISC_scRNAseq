"""Extract a small, figure-ready sample from the PBISC bulk-only generated draw files.
Usage: python extract_file_component.py <dataset> <out_dir>   dataset in sepsis_k1000, sepsis_k4000, covid_k1000, covid_k4000
Per draw: 4,000 random cells -> total UMI, detected genes, raw counts of marker genes, labels; plus
label counts over all cells and the per-gene sum of CP10k / detection over the sampled cells."""
import h5py, numpy as np, sys, json, os, time
ROOT = os.environ.get("PBISC_ROOT", ".")
PATHS = {
 "sepsis_k1000": "bulkonly_lora5_20261006/ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw{j}.h5ad",
 "sepsis_k4000": "bulkonly_lora5_k4000_20261006/ye2026_sepsis/lora5k4_bulkonly_ye2026_sepsis_draw{j}.h5ad",
 "covid_k1000": "bulkonly_lora5_20261006/imperial_eskd/lora5_bulkonly_imperial_eskd_draw{j}.h5ad",
 "covid_k4000": "bulkonly_lora5_k4000_20261006/imperial_eskd/lora5k4_bulkonly_imperial_eskd_draw{j}.h5ad"}
GENES = "CD3E IL7R CCR7 CD4 CD8A NKG7 GNLY MS4A1 CD79A JCHAIN CD14 LYZ S100A8 FCGR3A FCER1A CLEC4C PPBP HBB MKI67".split()
N_PER_DRAW = 4000
ds, out = sys.argv[1], sys.argv[2]
def cat(f, grp, col):
    g = f[grp][col]; c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]); return c, g["codes"][:]
res = {"genes": GENES}; t0 = time.time()
for j in range(5):
    p = os.path.join(ROOT, PATHS[ds].format(j=j))
    with h5py.File(p, "r") as f:
        sc, scode = cat(f, "var", "gene_symbol"); sym = sc[scode]
        meas = f["var"]["measured_in_cohort"][:]
        gidx = np.array([np.where(sym == g)[0][0] for g in GENES])
        ip = f["X"]["indptr"][:]; n = len(ip) - 1
        rows = np.sort(np.random.default_rng(20261018 + j).choice(n, N_PER_DRAW, replace=False))
        D = f["X"]["data"]; I = f["X"]["indices"]
        tot = np.zeros(N_PER_DRAW, np.int64); ng = np.zeros(N_PER_DRAW, np.int32)
        M = np.zeros((N_PER_DRAW, len(GENES)), np.int32)
        gsum = np.zeros(len(sym)); gdet = np.zeros(len(sym)); gcnt = np.zeros(len(sym))
        for k, r in enumerate(rows):
            d = D[ip[r]:ip[r + 1]]; ix = I[ip[r]:ip[r + 1]]
            t = d.sum(); tot[k] = t; ng[k] = len(d)
            gsum[ix] += np.log1p(d / max(t, 1) * 1e4); gdet[ix] += 1; gcnt[ix] += d
            for q, gi in enumerate(gidx):
                h = np.searchsorted(ix, gi)
                if h < len(ix) and ix[h] == gi: M[k, q] = d[h]
        labs = {}
        for col in ["celltype_l1_pred", "aifi_l1", "aifi_l2", "group", "sample_id"]:
            c, code = cat(f, "obs", col); labs[col] = (c, code)
        res[f"draw{j}"] = dict(rows=rows, total=tot, ngenes=ng, marker_counts=M, gene_logcp10k_sum=gsum, gene_det=gdet, gene_count_sum=gcnt,
            **{f"{col}_sample": labs[col][0][labs[col][1][rows]] for col in labs},
            **{f"{col}_allcounts": np.bincount(labs[col][1], minlength=len(labs[col][0])) for col in ["celltype_l1_pred", "aifi_l1", "aifi_l2", "group"]},
            **{f"{col}_cats": labs[col][0] for col in ["celltype_l1_pred", "aifi_l1", "aifi_l2", "group"]},
            n_cells=n, n_samples=len(labs["sample_id"][0]))
        res["gene_symbol"] = sym; res["measured"] = meas
    print(ds, "draw", j, "done", round(time.time() - t0), "s", flush=True)
flat = {}
for k, v in res.items():
    if isinstance(v, dict):
        for kk, vv in v.items(): flat[f"{k}__{kk}"] = np.asarray(vv)
    else: flat[k] = np.asarray(v)
np.savez_compressed(os.path.join(out, f"{ds}.npz"), **flat); print("saved", ds)
