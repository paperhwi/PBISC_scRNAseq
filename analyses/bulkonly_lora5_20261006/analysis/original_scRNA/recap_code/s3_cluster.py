"""Steps 3/5 (Ye 2026): HVG (vst, 2,000) -> LogNormalize -> ScaleData(regress nCount_RNA + percent.mito)
-> PCA(50) -> batch correction across samples -> top 10 corrected PCs -> SNN/kNN graph (k=20) -> clustering
-> UMAP. Used for the whole dataset (major lineages) and again per lineage (re-clustering).
fastMNN (paper) is replaced by Harmony on the same PCs, with sample (orig.ident) as batch.
Stages are check-pointed to <name>_*.npz so each heavy step runs in a fresh process.
Usage: python s3_cluster.py <name> <cells.npy | all> <stage: embed|graph|leiden|umap|all> [resolutions...]
"""
import sys, time, gc, os, subprocess, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/sepsis"
name, cfile, stage = sys.argv[1], sys.argv[2], sys.argv[3]; res = [float(x) for x in sys.argv[4:]] or [0.8]
t0 = time.time(); log = lambda *m: print(name, stage, round(time.time() - t0), *m, flush=True)
P = f"{W}/run_{name}"; os.makedirs(P, exist_ok=True)

if stage == "all":
    for s in ["embed", "graph", "leiden", "umap"]:
        r = subprocess.run([sys.executable, __file__, name, cfile, s] + sys.argv[4:])
        if r.returncode: sys.exit(r.returncode)
    sys.exit(0)

if stage == "embed":
    st = store.Store(f"{W}/store")
    cells = np.arange(st.n_cells) if cfile == "all" else np.sort(np.load(cfile))
    hvg, info = store.vst_hvg(st, cells, n_top=2000); log("hvg")
    X = store.lognorm_cols(st, hvg, cells); log("lognorm", X.shape, X.nnz)
    o = st.obs.iloc[cells]
    pcs, loadings, ev = store.regress_scale_pca(X, np.column_stack([st.ncount[cells], o["percent_mito"].to_numpy()]), n_pcs=50)
    del X; gc.collect(); log("pca")
    Z = store.run_harmony(pcs, o["sample"].to_numpy()); log("harmony")
    np.savez(f"{P}/embed.npz", cells=cells, hvg=hvg, pcs=pcs, harmony=Z, ev=ev, loadings=loadings)

elif stage == "graph":
    import scanpy as sc, anndata as ad
    e = np.load(f"{P}/embed.npz")
    a = ad.AnnData(np.zeros((len(e["cells"]), 0), dtype=np.float32)); a.obsm["X_h"] = e["harmony"][:, :10]
    sc.pp.neighbors(a, use_rep="X_h", n_neighbors=20, random_state=2021); log("neighbors")
    sp.save_npz(f"{P}/conn.npz", a.obsp["connectivities"].tocsr()); sp.save_npz(f"{P}/dist.npz", a.obsp["distances"].tocsr())

elif stage == "leiden":
    import igraph as ig, leidenalg
    C = sp.load_npz(f"{P}/conn.npz"); U = sp.triu(C, k=1).tocoo(); del C
    g = ig.Graph(n=U.shape[0], edges=np.column_stack([U.row, U.col]).astype(np.int64), directed=False)
    g.es["weight"] = U.data.astype(float); del U; gc.collect(); log("igraph", g.ecount())
    out = {}
    for r in res:
        part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight",
                                        resolution_parameter=r, seed=2021, n_iterations=2)
        out[f"leiden_{r:g}"] = np.array(part.membership); log("leiden", r, len(set(part.membership)))
    np.savez(f"{P}/leiden.npz", **out)

elif stage == "umap":
    import scanpy as sc, anndata as ad
    e = np.load(f"{P}/embed.npz")
    a = ad.AnnData(np.zeros((len(e["cells"]), 0), dtype=np.float32)); a.obsm["X_h"] = e["harmony"][:, :10]
    a.obsp["connectivities"] = sp.load_npz(f"{P}/conn.npz"); a.obsp["distances"] = sp.load_npz(f"{P}/dist.npz")
    a.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances",
                          "params": {"n_neighbors": 20, "method": "umap", "use_rep": "X_h", "random_state": 2021}}
    init = "spectral"
    if a.n_obs > 600000:      # spectral init on >600k cells exceeds memory here; start from the corrected PCs
        h = e["harmony"][:, :2].astype(np.float64); a.obsm["X_init"] = (h - h.mean(0)) / h.std(0) * 3; init = "X_init"
    sc.tl.umap(a, random_state=2021, min_dist=0.3, init_pos=init); log("umap")
    np.save(f"{P}/umap.npy", a.obsm["X_umap"])
log("done")
