"""Subsample cells (raw counts + labels) for the real-vs-generated random-forest test. h5py only; run from F:/PBISC.
   python extract_rf_cells.py sepsis_real|covid_real|k1000|k4000 [budget_s]
real files: cell_transcriptome_validation_20261014/inputs/*_real_heldout_panel_counts.h5ad (patients not used for any LoRA adapter)
generated:  Harmony-integrated 305k sets (500 cells / sample / draw), layers['counts']"""
import sys, os, time, json, h5py, numpy as np
T0 = time.time(); src = sys.argv[1]; BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 140
A = "bulkonly_lora5_20261006/analysis"
OUT = f"{A}/integration_report_20261012/rf_realism_20261017/_extract"; os.makedirs(OUT, exist_ok=True)
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
def col(g, c):
    x = g[c]
    if isinstance(x, h5py.Group): return dec(x["categories"][:])[x["codes"][:]] if "categories" in x else dec(x["values"][:])
    return dec(x[:])
CFG = {"sepsis_real": (f"{A}/cell_transcriptome_validation_20261014/inputs/ye2026_sepsis_real_heldout_panel_counts.h5ad", "X", 700, "patient_id",
                       ["major", "fine", "patient_id", "sample_id", "group"]),
       "covid_real": (f"{A}/cell_transcriptome_validation_20261014/inputs/imperial_eskd_real_heldout_panel_counts.h5ad", "X", 220, "sample_id",
                      ["major", "fine", "patient_id", "sample_id", "group"]),
       "k1000": (f"{A}/harmony/sepsis_k1000/sepsis_k1000_harmony_clustering.h5ad", "layers/counts", 210, "sample_id",
                 ["celltype_l1_pred", "aifi_l2", "draw", "sample_id", "group"]),
       "k4000": ("bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000/sepsis_k4000_harmony_clustering.h5ad", "layers/counts", 210, "sample_id",
                 ["celltype_l1_pred", "aifi_l2", "draw", "sample_id", "group"])}
path, xkey, per, by, cols = CFG[src]
h = h5py.File(path, "r"); o = h["obs"]; meta_f = f"{OUT}/{src}_meta.npz"
if not os.path.exists(meta_f):
    key = col(o, by); rng = np.random.default_rng(17)
    rows = np.sort(np.concatenate([rng.permutation(np.where(key == k)[0])[:per] for k in np.unique(key)]))
    sym = col(h["var"], "gene_symbol"); gid = col(h["var"], "gene_id")
    np.savez_compressed(meta_f, rows=rows, gene_symbol=sym, gene_id=gid, **{c: col(o, c)[rows] for c in cols})
m = np.load(meta_f, allow_pickle=True); rows = m["rows"]; X = h[xkey]; ip = X["indptr"][:].astype(np.int64); n = len(ip) - 1; B = 60000
for i, a in enumerate(range(0, n, B)):
    p = f"{OUT}/{src}_block{i:02d}.npz"
    if os.path.exists(p): continue
    if time.time() - T0 > BUDGET: print("budget reached; run again", flush=True); sys.exit(0)
    b = min(a + B, n); r = rows[(rows >= a) & (rows < b)]
    d = X["data"][ip[a]:ip[b]]; ix = X["indices"][ip[a]:ip[b]]
    s, e = ip[r] - ip[a], ip[r + 1] - ip[a]
    sel = np.concatenate([np.arange(x, y) for x, y in zip(s, e)]) if len(r) else np.array([], np.int64)
    np.savez_compressed(p + ".tmp.npz", data=d[sel].astype(np.int32), indices=ix[sel].astype(np.int32), indptr=np.concatenate([[0], np.cumsum(e - s)]), rows=r)
    os.replace(p + ".tmp.npz", p); print(src, i, len(r), f"{time.time() - T0:.0f}s", flush=True)
print(src, "DONE", len(rows), "cells")
