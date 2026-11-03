"""Stephenson et al. 2025 (Cell Genomics) compartment re-analysis, following the STAR Methods:
  normalize_total(1e4) + log1p -> highly_variable_genes(flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5)
  -> [T/NK only: regress_out(total_counts, pct_mito)] -> scale(max_value=10) -> PCA -> harmony (myeloid)
  -> neighbors -> leiden -> rank_genes_groups -> manual annotation (iterative re-clustering).
The scaled matrix is never formed (exact matrix-free PCA, clip at 10 omitted).
Usage: python c1_compartment.py <name> <obs-column> <value> [regress] [harmony] [n_pcs_for_neighbors] [resolutions...]
"""
import sys, json, time, numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, "/home/claude/work"); import store

name, col, val = sys.argv[1], sys.argv[2], sys.argv[3]
regress = sys.argv[4] == "1"; harmony = sys.argv[5] == "1"; npc = int(sys.argv[6])
res = [float(r) for r in sys.argv[7:]] or [1.0]
st = store.Store("/home/claude/work/covid/store")
obs = st.obs
cells = np.sort(np.load(val)) if col == "@file" else np.where(obs[col].astype(str).isin(val.split("|")).to_numpy())[0]
t0 = time.time()
mt = st.var.gene.str.upper().str.startswith("MT-").to_numpy()
import os, re
exc = os.environ.get('HVG_EXCLUDE')
exm = st.var.gene.str.contains(exc, regex=True).to_numpy() if exc else None
hvg, hinfo = store.seurat_disp_hvg(st, cells, exclude=exm)
print(name, "cells", len(cells), "hvg", len(hvg), round(time.time() - t0), flush=True)
X = store.lognorm_cols(st, hvg, cells)
cov = None
if regress:
    o = obs.iloc[cells]
    cov = np.column_stack([st.ncount[cells], o["percentMito"].to_numpy()])
pcs, load, ev = store.regress_scale_pca(X, cov, n_pcs=50)
del X
rep = pcs
if harmony:
    rep = store.run_harmony(pcs, obs["orig.ident"].astype(str).to_numpy()[cells])
a = ad.AnnData(obs=obs.iloc[cells].copy())
a.obsm["X_pca"] = pcs; a.obsm["X_rep"] = rep
sc.pp.neighbors(a, use_rep="X_rep", n_pcs=npc, n_neighbors=15, random_state=0)
for r in res:
    sc.tl.leiden(a, resolution=r, key_added=f"leiden_{r:g}", random_state=0, flavor="igraph", n_iterations=2)
sc.tl.umap(a, random_state=0)
a.uns["hvg_genes"] = st.var.gene.to_numpy()[hvg].astype(str)
a.obs["cell_index"] = cells
a.write_h5ad(f"/home/claude/work/covid/{name}.h5ad")
print(name, "done", round(time.time() - t0), {f"leiden_{r:g}": a.obs[f"leiden_{r:g}"].nunique() for r in res}, flush=True)
