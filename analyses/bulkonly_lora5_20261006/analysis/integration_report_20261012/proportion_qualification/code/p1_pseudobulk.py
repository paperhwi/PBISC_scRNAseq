"""One pass over the real sepsis reference: summed UMI counts per (sample, fine cell type) over all genes, plus cell counts.
Gives per-sample truth proportions, cell-type signatures (any sample subset) and per-sample pseudobulk."""
import sys, numpy as np, pandas as pd, scipy.sparse as sp, time
sys.path.insert(0, "/home/claude/work"); import store
t0 = time.time(); O = "/home/claude/work/prop"
st = store.Store("/home/claude/work/sepsis/store")
ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                            usecols=["sample", "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme"]) for i in range(3)], ignore_index=True)
assert len(ro) == st.n_cells
ok = ro.in_paper_scheme.astype(str).eq("True").to_numpy()
samp = pd.Categorical(ro["sample"].astype(str)); fine = pd.Categorical(ro.celltype_ye2026.astype(str))
nS, nF = len(samp.categories), len(fine.categories)
key = samp.codes.astype(np.int64) * nF + fine.codes
keys_used = np.unique(key[ok]); kpos = {k: i for i, k in enumerate(keys_used)}
row_of = np.full(nS * nF, -1); row_of[keys_used] = np.arange(len(keys_used))
acc = np.zeros((len(keys_used), st.n_genes), dtype=np.float64)
cells = np.where(ok)[0]
for i, (c, X) in enumerate(st.chunks(cells, size=20000)):
    r = row_of[key[c]]
    A = sp.csr_matrix((np.ones(len(c), np.float32), (r, np.arange(len(c)))), shape=(len(keys_used), len(c)))
    acc += (A @ X).toarray()
    if i % 5 == 0: print(f"[{time.time()-t0:5.0f}s] {c[-1]} / {st.n_cells}", flush=True)
ncell = np.bincount(row_of[key[ok]], minlength=len(keys_used))
idx = pd.DataFrame({"sample": samp.categories[keys_used // nF], "fine": fine.categories[keys_used % nF], "n_cells": ncell})
f2m = ro[ok].groupby("celltype_ye2026").celltype_major_ye2026.agg(lambda s: s.value_counts().index[0])
idx["major"] = idx.fine.map(f2m)
grp = st.obs.assign(sample=ro["sample"].astype(str).values).groupby("sample").group.first()
idx["group"] = idx["sample"].map(grp)
idx.to_csv(f"{O}/sepsis_pb_index.csv", index=False)
np.save(f"{O}/sepsis_pb_counts.npy", acc.astype(np.float32))
st.var.to_csv(f"{O}/sepsis_pb_var.csv", index=False)
print("done", acc.shape, idx.n_cells.sum(), f"{time.time()-t0:.0f}s")
