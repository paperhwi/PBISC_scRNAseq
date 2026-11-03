"""Write data/sepsis_<name>_myeloid_subcluster.h5ad with h5py + numpy only (no anndata/scipy needed), following the
AnnData on-disk spec, so it can run where only h5py is installed. Run from F:/PBISC:
    python write_myeloid_h5ad.py real|k1000|k4000
Inputs: _extract/<name>_block*.npz (raw counts, CSR pieces), results/<name>_myeloid_subcluster_* (from myeloid_subcluster.ipynb),
        celltype_umap_20261016/data/sepsis_<name>_umap_obs.npz (global UMAP).
Output: X = log1p(CP10k) on all genes (identical to sc.pp.normalize_total(1e4) + sc.pp.log1p), layers['counts'] = raw counts."""
import sys, os, glob, json, time, h5py, numpy as np
name = sys.argv[1]; T0 = time.time()
REP = os.environ.get("REP", "bulkonly_lora5_20261006/analysis/integration_report_20261012")
W = f"{REP}/myeloid_subcluster_20261016"; E, R, D = f"{W}/_extract", f"{W}/results", f"{W}/data"; os.makedirs(D, exist_ok=True)
STR = h5py.string_dtype("utf-8")
def enc(o, t, v): o.attrs["encoding-type"] = t; o.attrs["encoding-version"] = v
def w_array(g, k, a): d = g.create_dataset(k, data=a); enc(d, "array", "0.2.0")
def w_strings(g, k, a): d = g.create_dataset(k, data=np.asarray(a, dtype=object), dtype=STR); enc(d, "string-array", "0.2.0")
def w_cat(g, k, a):
    a = np.asarray(a).astype(str); cats, codes = np.unique(a, return_inverse=True)
    if all(c.isdigit() for c in cats):            # Leiden labels: numeric order (0, 1, 2, ..., 10) instead of lexicographic
        order = np.argsort(cats.astype(int)); rank = np.empty_like(order); rank[order] = np.arange(len(order)); cats, codes = cats[order], rank[codes]
    c = g.create_group(k); enc(c, "categorical", "0.2.0"); c.attrs["ordered"] = False
    w_array(c, "codes", codes.astype(np.int8 if len(cats) < 127 else np.int16 if len(cats) < 32000 else np.int32)); w_strings(c, "categories", cats)
def w_df(g, k, index, cols):
    d = g.create_group(k); enc(d, "dataframe", "0.2.0"); d.attrs["_index"] = "_index"
    d.attrs.create("column-order", np.array(list(cols), dtype=object), dtype=STR) if cols else d.attrs.create("column-order", np.array([], dtype=np.float64))
    w_strings(d, "_index", index)
    for c, v in cols.items():
        (w_array(d, c, v) if np.asarray(v).dtype.kind in "fiub" else w_cat(d, c, v))
def empty_dict(g, k): x = g.create_group(k); enc(x, "dict", "0.1.0"); return x

blocks = sorted(glob.glob(f"{E}/{name}_block*.npz"))
meta = [np.load(f) for f in blocks]; rows = np.concatenate([m["rows"] for m in meta])
nnz = int(sum(int(m["indptr"][-1]) for m in meta)); n = len(rows)
var_names = np.load(f"{R}/{name}_var_names.npy", allow_pickle=True).astype(str); m_genes = len(var_names)
gid = np.load(f"{E}/{'real_gene_ids' if name == 'real' else 'panel_gene_ids'}.npy", allow_pickle=True).astype(str)
pid = set(np.load(f"{E}/panel_gene_ids.npy", allow_pickle=True).astype(str))
ob = np.load(f"{R}/{name}_myeloid_subcluster_obs.npz", allow_pickle=True)
obs_names = ob["obs_names"].astype(str); assert len(obs_names) == n
expect = np.array([f"{name}_{r}" for r in rows]); assert (obs_names == expect).all(), "row order mismatch between extract and results"
out = f"{D}/sepsis_{name}_myeloid_subcluster.h5ad"; tmp = out + ".tmp"; state_f = out + ".progress.json"
BUDGET = float(os.environ.get("BUDGET", "1e9"))      # seconds; resumable: re-run the same command until it prints "wrote"
state = json.load(open(state_f)) if os.path.exists(state_f) and os.path.exists(tmp) else {"done": 0, "o_nz": 0, "o_r": 0}
with h5py.File(tmp, "a" if state["done"] else "w") as h:
    if state["done"] == 0:
        enc(h, "anndata", "0.1.0")
        X = h.create_group("X"); L = h.create_group("layers"); enc(L, "dict", "0.1.0"); C = L.create_group("counts")
        for g in (X, C):
            enc(g, "csr_matrix", "0.1.0"); g.attrs["shape"] = np.array([n, m_genes], dtype=np.int64)
            ck = (max(1, min(nnz, 4_000_000)),)   # chunked: space is allocated as written (contiguous storage would be zero-filled up front)
            g.create_dataset("data", (nnz,), np.float32, chunks=ck); g.create_dataset("indices", (nnz,), np.int32, chunks=ck); g.create_dataset("indptr", data=np.zeros(n + 1, np.int64))
    X, C = h["X"], h["layers/counts"]; dX, dC, iX, iC = X["data"], C["data"], X["indices"], C["indices"]
    o_nz, o_r = state["o_nz"], state["o_r"]
    for bi, f in enumerate(blocks):
        if bi < state["done"]: continue
        if time.time() - T0 > BUDGET:
            json.dump({"done": bi, "o_nz": o_nz, "o_r": o_r}, open(state_f, "w")); print(f"budget reached at block {bi}/{len(blocks)}; run again"); sys.exit(0)
        z = np.load(f); d = z["data"].astype(np.float32); ix = z["indices"].astype(np.int32); p = z["indptr"].astype(np.int64); k = len(p) - 1
        tot = np.add.reduceat(d, p[:-1]) if len(d) else np.zeros(k, np.float32); tot[np.diff(p) == 0] = 0
        lg = np.log1p(d * np.repeat((1e4 / np.maximum(tot, 1)).astype(np.float32), np.diff(p)))
        dX[o_nz:o_nz + len(d)] = lg; dC[o_nz:o_nz + len(d)] = d; iX[o_nz:o_nz + len(d)] = ix; iC[o_nz:o_nz + len(d)] = ix
        for g in (X, C): g["indptr"][o_r + 1:o_r + k + 1] = p[1:] + o_nz
        o_nz += len(d); o_r += k; h.flush()
        print(os.path.basename(f), k, "cells", f"{time.time() - T0:.0f}s", flush=True)
    assert o_r == n and o_nz == nnz
    cols = {}
    for k in ob.files:
        if k == "obs_names": continue
        cols[k[5:] if k.startswith("num__") else k] = ob[k]
    w_df(h, "obs", obs_names, cols)
    hvg = np.load(f"{R}/{name}_myeloid_subcluster_hvg.npy").astype(bool)
    w_df(h, "var", var_names, {"gene_id": gid, "in_panel": np.isin(gid, list(pid)), "highly_variable": hvg})
    om = empty_dict(h, "obsm")
    for k in ["X_pca", "X_pca_harmony", "X_umap"]: w_array(om, k, np.load(f"{R}/{name}_myeloid_subcluster_{k}.npy"))
    zg = np.load(f"{REP}/celltype_umap_20261016/data/sepsis_{name}_umap_obs.npz", allow_pickle=True); w_array(om, "X_umap_global", zg["umap"][rows])
    for k in ["varm", "obsp", "varp"]: empty_dict(h, k)
    u = empty_dict(h, "uns"); p = json.load(open(f"{R}/{name}_myeloid_subcluster_params.json"))
    pu = empty_dict(u, "subcluster_params")
    for k, v in p.items():
        if isinstance(v, list): w_array(pu, k, np.asarray(v, float))
        elif isinstance(v, str): d = pu.create_dataset(k, data=v, dtype=STR); enc(d, "string", "0.2.0")
        else: d = pu.create_dataset(k, data=v); enc(d, "numeric-scalar", "0.2.0")
os.replace(tmp, out); os.path.exists(state_f) and os.replace(state_f, state_f + ".done"); print("wrote", out, n, "x", m_genes, f"nnz={nnz:,}", f"{time.time() - T0:.0f}s")
