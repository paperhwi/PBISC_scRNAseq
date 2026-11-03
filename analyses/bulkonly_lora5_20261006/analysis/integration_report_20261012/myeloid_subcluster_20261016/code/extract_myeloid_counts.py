"""Extract raw counts of myeloid cells (h5py only, resumable; run from F:/PBISC).
real : GSE279452 recap, in_paper_scheme & celltype_major_ye2026 in {Monocytes, Macrophages, DC}  -> all 37,474 genes
k1000/k4000 : Harmony-integrated 305k generated cells, celltype_l1_pred in {Monocyte, DC}       -> 5,975 panel genes (layers['counts'])
Output: _extract/<src>_block*.npz with CSR pieces (data uint16/int32, indices uint16, indptr, global row ids)."""
import sys, os, time, h5py, numpy as np
OUT = "bulkonly_lora5_20261006/analysis/integration_report_20261012/myeloid_subcluster_20261016/_extract"; T0 = time.time(); BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 140
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
def catcol(o, c):
    g = o[c]
    return dec(g["categories"][:])[g["codes"][:]] if isinstance(g, h5py.Group) else dec(g[:])
def save(path, X, rows, sel_ip, data, idx):
    dt = np.uint16 if data.max(initial=0) < 65535 else np.int32
    np.savez_compressed(path, data=data.astype(dt), indices=idx.astype(np.uint16), indptr=sel_ip, rows=rows)
def take(ip, data_ds, idx_ds, rows, a, b):
    """rows (sorted, global) inside [a, b); read the block once, keep selected rows."""
    off = ip[a]; d = data_ds[ip[a]:ip[b]]; ix = idx_ds[ip[a]:ip[b]]
    starts = ip[rows] - off; ends = ip[rows + 1] - off; lens = ends - starts
    sel = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)]) if len(rows) else np.array([], int)
    return d[sel], ix[sel], np.concatenate([[0], np.cumsum(lens)])
src = sys.argv[1]
if src == "real":
    h = h5py.File("bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad", "r"); o = h["obs"]
    ok = o["in_paper_scheme"][:].astype(bool) & np.isin(catcol(o, "celltype_major_ye2026"), ["Monocytes", "Macrophages", "DC"])
    rows_all = np.where(ok)[0]; X = h["X"]; ip = X["indptr"][:]; B = 50000
    if not os.path.exists(f"{OUT}/real_genes.npy"): np.save(f"{OUT}/real_genes.npy", catcol(h["var"], "gene_symbol")); np.save(f"{OUT}/real_rows.npy", rows_all)
    for i, a in enumerate(range(0, len(ip) - 1, B)):
        p = f"{OUT}/real_block{i:02d}.npz"
        if os.path.exists(p): continue
        if time.time() - T0 > BUDGET: print("budget reached; rerun to continue"); break
        b = min(a + B, len(ip) - 1); r = rows_all[(rows_all >= a) & (rows_all < b)]
        d, ix, sip = take(ip, X["data"], X["indices"], r, a, b); save(p + ".tmp.npz", X, r, sip, d, ix); os.replace(p + ".tmp.npz", p)
        print(i, len(r), "cells", f"{time.time() - T0:.0f}s", flush=True)
    else: print("real DONE", len(rows_all))
else:
    f = {"k1000": "bulkonly_lora5_20261006/analysis/harmony/sepsis_k1000/sepsis_k1000_harmony_clustering.h5ad",
         "k4000": "bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000/sepsis_k4000_harmony_clustering.h5ad"}[src]
    h = h5py.File(f, "r"); o = h["obs"]; r = np.where(np.isin(catcol(o, "celltype_l1_pred"), ["Monocyte", "DC"]))[0]
    C = h["layers/counts"]; ip = C["indptr"][:].astype(np.int64)
    d, ix, sip = take(ip, C["data"], C["indices"], r, 0, len(ip) - 1)
    np.save(f"{OUT}/{src}_genes.npy", catcol(h["var"], "gene_symbol"))
    save(f"{OUT}/{src}_block00.npz", C, r, sip, d, ix); print(src, "DONE", len(r), "cells", f"{time.time() - T0:.0f}s")
