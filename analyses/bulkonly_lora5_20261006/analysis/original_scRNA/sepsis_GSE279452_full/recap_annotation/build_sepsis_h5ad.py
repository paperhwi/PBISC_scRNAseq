"""Assemble the Ye 2026 sepsis (GSE279452) recapitulation .h5ad on the local computer.

Reads the raw 10x matrices downloaded from GEO, keeps the QC-passed singlets of the re-analysis (same order as the
annotation table), drops mitochondrial genes (as the authors did), and streams the raw counts into X.
Resumable: run repeatedly until it prints 'ALL SAMPLES WRITTEN', then run with --finalize.

usage: python build_sepsis_h5ad.py <geo_dir> <annotation_dir> <out.h5ad> [--finalize] [--budget SECONDS]
"""
import sys, os, json, time, gzip, glob
import numpy as np, pandas as pd, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5ad_lite as hl

geo, ann, out = sys.argv[1], sys.argv[2], sys.argv[3]
finalize = "--finalize" in sys.argv
budget = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else 150
t0 = time.time()
state_f = out + ".state.json"
obs = pd.concat([pd.read_csv(p, dtype={"barcode": str}, keep_default_na=False) for p in sorted(glob.glob(f"{ann}/sepsis_obs_part*.csv.gz"))], ignore_index=True)
for c in ["nCount_RNA", "nCount_RNA_raw", "nFeature_RNA", "percent_mito", "doublet_score"]:
    obs[c] = pd.to_numeric(obs[c])
obs["in_paper_scheme"] = obs["in_paper_scheme"].astype(str).eq("True")
samples = list(dict.fromkeys(obs["sample_dir"]))
feat = pd.read_csv(glob.glob(f"{geo}/{samples[0]}/*.features.tsv.gz")[0], sep="\t", header=None,
                   names=["gene_ids", "gene_symbol", "feature_types"])
mt = feat.gene_symbol.str.upper().str.startswith("MT-").to_numpy()
newidx = np.full(len(feat), -1, np.int64); newidx[~mt] = np.arange((~mt).sum())
G = int((~mt).sum())

if not os.path.exists(state_f):
    with h5py.File(out, "w") as f:
        f.attrs["encoding-type"] = "anndata"; f.attrs["encoding-version"] = "0.1.0"
        g = f.create_group("X"); hl._attrs(g, "csr_matrix", "0.1.0")
        g.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(4_000_000,))
        g.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(4_000_000,))
        g.create_dataset("indptr", shape=(1,), maxshape=(None,), dtype=np.int64, chunks=(100_000,), data=[0])
    json.dump({"done": []}, open(state_f, "w"))

state = json.load(open(state_f))
if not finalize:
    for s in samples:
        if s in state["done"]:
            continue
        if time.time() - t0 > budget:
            break
        bc = pd.read_csv(f"{geo}/{s}/{s}.barcodes.tsv.gz", header=None)[0].astype(str).to_numpy()
        df = pd.read_csv(f"{geo}/{s}/{s}.matrix.mtx.gz", sep=" ", header=None, comment="%", dtype=np.int64,
                         compression="gzip").iloc[1:]
        gi = df[0].to_numpy() - 1; ci = df[1].to_numpy() - 1; v = df[2].to_numpy(); del df
        keep_bc = obs.loc[obs.sample_dir == s, "barcode"].to_numpy()
        pos = pd.Series(np.arange(len(bc)), index=bc)
        rows = pos.loc[keep_bc].to_numpy()                       # raw cell index of each kept cell, in obs order
        rank = np.full(len(bc), -1, np.int64); rank[rows] = np.arange(len(rows))
        r = rank[ci]; ng = newidx[gi]; m = (r >= 0) & (ng >= 0)
        r, ng, v = r[m], ng[m], v[m]
        o = np.lexsort((ng, r)); r, ng, v = r[o], ng[o], v[o]
        indptr = np.zeros(len(rows) + 1, np.int64); indptr[1:] = np.cumsum(np.bincount(r, minlength=len(rows)))
        with h5py.File(out, "a") as f:
            X = f["X"]; nnz = X["data"].shape[0]; n_ip = X["indptr"].shape[0]
            X["data"].resize((nnz + len(v),)); X["data"][nnz:] = v.astype(np.float32)
            X["indices"].resize((nnz + len(v),)); X["indices"][nnz:] = ng.astype(np.int32)
            X["indptr"].resize((n_ip + len(rows),)); X["indptr"][n_ip:] = indptr[1:] + nnz
        state["done"].append(s); json.dump(state, open(state_f, "w"))
        print(f"{s}: {len(rows)} cells, {len(v)} nnz, {time.time() - t0:.0f}s", flush=True)
    left = len(samples) - len(state["done"])
    print("ALL SAMPLES WRITTEN" if left == 0 else f"{left} samples left", flush=True)
    sys.exit(0)

# ---- finalize: shape, obs, var, obsm, uns
assert len(state["done"]) == len(samples), "samples missing"
emb = {}
for fpath in sorted(glob.glob(f"{ann}/sepsis_obsm_*.npz")):
    z = np.load(fpath)
    for k in z.files:
        key, part = k.rsplit("__", 1); emb.setdefault(key, []).append((int(part), z[k]))
with h5py.File(out, "a") as f:
    X = f["X"]; assert X["indptr"].shape[0] == len(obs) + 1
    X.attrs["shape"] = np.array([len(obs), G], dtype=np.int64)
    for k in ["obs", "var", "obsm", "uns", "layers", "obsp", "varm", "varp"]:
        if k in f: del f[k]
    cols = [c for c in obs.columns if c not in ("cell_id",)]
    hl.write_dataframe(f, "obs", obs[cols], obs["cell_id"].astype(str).to_numpy())
    var = feat.loc[~mt, ["gene_ids", "feature_types"]].reset_index(drop=True)
    sym = feat.loc[~mt, "gene_symbol"].astype(str)
    dup = sym.duplicated(keep=False)
    names = np.where(dup, sym + "_" + feat.loc[~mt, "gene_ids"].astype(str), sym)
    var["gene_symbol"] = sym.to_numpy()
    hl.write_dataframe(f, "var", var, names)
    om = hl.empty_dict(f, "obsm")
    for key, parts in emb.items():
        arr = np.vstack([p for _, p in sorted(parts, key=lambda t: t[0])]).astype(np.float32)
        assert arr.shape[0] == len(obs), key
        hl.write_array(om, key, arr)
    u = hl.empty_dict(f, "uns")
    meta = json.load(open(f"{ann}/sepsis_uns.json"))
    for k, v in meta.items():
        d = u.create_dataset(k, data=str(v), dtype=hl.STR); hl._attrs(d, "string", "0.2.0")
    for k in ["layers", "obsp", "varm", "varp"]:
        hl.empty_dict(f, k)
os.remove(state_f)
print("FINALIZED", out, len(obs), "cells x", G, "genes", flush=True)
