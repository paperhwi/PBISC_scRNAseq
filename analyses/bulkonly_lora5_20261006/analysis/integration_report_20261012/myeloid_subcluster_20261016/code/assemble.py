"""Assemble extracted CSR blocks into count AnnData files (X = raw UMI / generated counts) with obs."""
import os, glob, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, warnings; warnings.filterwarnings("ignore")
BASE = os.environ.get("PBISC_BASE", "F:/PBISC"); REP = f"{BASE}/bulkonly_lora5_20261006/analysis/integration_report_20261012"
E = f"{REP}/myeloid_subcluster_20261016/_extract"; OUT = f"{REP}/myeloid_subcluster_20261016/data"; os.makedirs(OUT, exist_ok=True)
U = f"{REP}/celltype_umap_20261016/data"
def cat(z, c): return z[c + "__cats"][z[c + "__codes"]].astype(str)
def load_blocks(pat, ncol):
    mats, rows = [], []
    for f in sorted(glob.glob(pat)):
        z = np.load(f); n = len(z["rows"])
        mats.append(sp.csr_matrix((z["data"].astype(np.float32), z["indices"].astype(np.int32), z["indptr"].astype(np.int64)), shape=(n, ncol))); rows.append(z["rows"])
    return sp.vstack(mats, format="csr"), np.concatenate(rows)
# real
g = np.load(f"{E}/real_genes.npy", allow_pickle=True).astype(str)
X, rows = load_blocks(f"{E}/real_block*.npz", len(g)); assert (rows == np.load(f"{E}/real_rows.npy")).all()
zr = np.load(f"{U}/sepsis_real_umap_obs.npz", allow_pickle=True)
obs = pd.DataFrame({c: cat(zr, c)[rows] for c in ["celltype_ye2026", "celltype_major_ye2026", "sample", "group", "age_group", "major_lineage"]})
obs["source"] = "real"; obs.index = [f"real_{r}" for r in rows]
var = pd.DataFrame(index=pd.Index(g).astype(str)); var.index.name = None
a = ad.AnnData(X, obs=obs, var=var); a.var_names_make_unique(); a.obsm["X_umap_global"] = zr["umap"][rows]
print(a, "nnz", X.nnz); a.write_h5ad(f"{OUT}/sepsis_real_myeloid_counts.h5ad", compression="gzip"); del a, X
for K in ["k1000", "k4000"]:
    g = np.load(f"{E}/{K}_genes.npy", allow_pickle=True).astype(str)
    X, rows = load_blocks(f"{E}/{K}_block*.npz", len(g)); z = np.load(f"{U}/sepsis_{K}_umap_obs.npz", allow_pickle=True)
    obs = pd.DataFrame({c: cat(z, c)[rows] for c in ["celltype_l1_pred", "aifi_l1", "aifi_l2", "draw", "sample_id", "group"]})
    obs["aifi_l2_confidence"] = z["aifi_l2_confidence__raw"][rows]; obs["source"] = K; obs.index = [f"{K}_{r}" for r in rows]
    a = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=pd.Index(g).astype(str))); a.var_names_make_unique(); a.obsm["X_umap_global"] = z["umap"][rows]
    print(a, "nnz", X.nnz); a.write_h5ad(f"{OUT}/sepsis_{K}_myeloid_counts.h5ad", compression="gzip")
