"""Step 1: real T/NK cells from the Ye 2026 recap, on the generated panel, thinned to generated depth.

Outputs (data/):
  real_tnk_panel.h5ad   sampled real T/NK cells; X = thinned counts, layers['raw'] = original panel counts
  real_truth_fractions.csv  per-sample author-label fractions over ALL cells (no subsampling)
  depth.json            library medians used for thinning
"""
import json, os
import numpy as np, pandas as pd, h5py, anndata as ad, scipy.sparse as sp

ROOT = "F:/PBISC/fig6_bulkonly_reproduction_20261020"
REC = "F:/PBISC/bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad"
GEN = "F:/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw{}.h5ad"
os.makedirs(f"{ROOT}/data", exist_ok=True)
rng = np.random.default_rng(20261020)
PER_SAMPLE = 1000

# ---- panel: measured generated genes, deduplicated by first occurrence, present in recap
g0 = ad.read_h5ad(GEN.format(0), backed="r")
gv = g0.var[g0.var.measured_in_cohort.to_numpy()].copy()
gv = gv[~gv.gene_symbol.duplicated()]
rid = ad.read_h5ad(REC, backed="r").var["gene_ids"].astype(str).to_numpy()
rpos = pd.Series(np.arange(len(rid)), index=rid)
rpos = rpos[~rpos.index.duplicated()]
gv = gv[gv.gene_id.isin(rpos.index)]
panel_ids = gv.gene_id.to_numpy(); panel_sym = gv.gene_symbol.to_numpy()
col_real = rpos.loc[panel_ids].to_numpy()
print("panel genes", len(panel_ids))

# ---- real obs and truth fractions over all cells
o = ad.read_h5ad(REC, backed="r").obs
lab = o.celltype_ye2026.astype(str)
lin = np.select([lab.str.startswith("T_CD4"), lab.str.startswith("T_CD8"), lab.str.startswith("NKT_"),
                 lab.str.startswith("NK_"), lab.eq("gdT")], ["CD4T", "CD8T", "NKT", "NK", "other"], "none")
o = o.assign(lab=lab.values, lin=lin)
TARGETS = {"F1": ("CD4T", "T_CD4_c02-GPR183-NR4A2"), "F2": ("CD8T", "T_CD8_c05-GNLY-CCL4hi"),
           "F3": ("NK", "NK_c03-PRF1-CCL4hi"), "F4": ("NKT", "NKT_c05-IFNG-CCL4hi"), "C1": ("CD4T", "T_CD4_c01-LEF1")}
rows = []
for s, d in o.groupby("sample", observed=True):
    for k, (par, st) in TARGETS.items():
        n = int((d.lin == par).sum()); kk = int((d.lab == st).sum())
        rows.append(dict(sample=s, group=d.group.iloc[0], age_group=d.age_group.iloc[0], infection_site=d.infection_site.iloc[0],
                         target=k, parent_n=n, state_n=kk))
pd.DataFrame(rows).to_csv(f"{ROOT}/data/real_truth_fractions.csv", index=False)

# ---- sample real T/NK cells
cand = np.where(o.lin.isin(["CD4T", "CD8T", "NKT", "NK", "other"]).to_numpy())[0]
cs = pd.Series(cand).groupby(o["sample"].to_numpy()[cand])
pick = np.sort(np.concatenate([rng.choice(v.to_numpy(), min(PER_SAMPLE, len(v)), replace=False) for _, v in cs]))
print("real cells sampled", len(pick))

# ---- read selected rows chunk-wise, keep panel columns
colmap = np.full(37474 + 10, -1); colmap[col_real] = np.arange(len(col_real))
data_l, ind_l, ptr = [], [], [0]
with h5py.File(REC, "r") as f:
    X = f["X"]; indptr = X["indptr"][:]
    CH = 40000
    for start in range(0, len(indptr) - 1, CH):
        stop = min(start + CH, len(indptr) - 1)
        sel = pick[(pick >= start) & (pick < stop)]
        if len(sel) == 0: continue
        a, b = indptr[start], indptr[stop]
        dat = X["data"][a:b]; ind = X["indices"][a:b]
        for r in sel:
            s0, s1 = indptr[r] - a, indptr[r + 1] - a
            ii = colmap[ind[s0:s1]]; keep = ii >= 0
            data_l.append(dat[s0:s1][keep]); ind_l.append(ii[keep]); ptr.append(ptr[-1] + int(keep.sum()))
        print("chunk", start, flush=True)
raw = sp.csr_matrix((np.concatenate(data_l).astype(np.float32), np.concatenate(ind_l), np.array(ptr)), shape=(len(pick), len(panel_ids)))
raw.data = np.round(raw.data)

# ---- thinning to generated depth (T-like vs NK-like)
glib = {"T": [], "NK": []}
for j in range(5):
    g = ad.read_h5ad(GEN.format(j), backed="r")
    m = g.var.gene_id.isin(panel_ids).to_numpy() & ~g.var.gene_symbol.duplicated().to_numpy()
    for L, lab1 in [("T", "T cell"), ("NK", "NK cell")]:
        idx = np.where((g.obs.celltype_l1_pred == lab1).to_numpy())[0]
        idx = np.sort(rng.choice(idx, min(5000, len(idx)), replace=False))
        glib[L].append(np.asarray(g.X[idx][:, np.where(m)[0]].sum(1)).ravel())
gmed = {L: float(np.median(np.concatenate(v))) for L, v in glib.items()}
obs = o.iloc[pick].copy()
rlib = np.asarray(raw.sum(1)).ravel()
LL = np.where(obs.lin.eq("NK").to_numpy(), "NK", "T")
rmed = {L: float(np.median(rlib[LL == L])) for L in ["T", "NK"]}
p = np.array([min(1.0, gmed[L] / rmed[L]) for L in LL])
thin = raw.copy(); rep = np.repeat(p, np.diff(thin.indptr))
thin.data = rng.binomial(thin.data.astype(np.int64), rep).astype(np.float32); thin.eliminate_zeros()
json.dump(dict(gen_median=gmed, real_median=rmed, keep_prob={L: min(1, gmed[L] / rmed[L]) for L in rmed}, n_cells=int(len(pick)),
               per_sample=PER_SAMPLE, panel_genes=int(len(panel_ids))), open(f"{ROOT}/data/depth.json", "w"), indent=1)
print(gmed, rmed)

A = ad.AnnData(thin, obs=obs[["sample", "group", "age_group", "infection_site", "lab", "lin"]].astype(str).reset_index(drop=True),
               var=pd.DataFrame({"gene_id": panel_ids}, index=panel_sym))
A.layers["raw"] = raw
A.write_h5ad(f"{ROOT}/data/real_tnk_panel.h5ad")
print(A)
