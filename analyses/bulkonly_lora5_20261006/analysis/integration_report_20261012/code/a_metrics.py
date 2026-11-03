"""A: integration + clustering metrics for generated (K1000/K4000 x sepsis/covid) and real reference datasets."""
import sys, json, numpy as np, pandas as pd, anndata as ad, scanpy as sc
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI
sys.path.insert(0, "/home/claude/work/report"); sys.path.insert(0, "/home/claude/work")
import metrics as M, store
U = "/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; R = "/home/claude/work/report"
res = {}
for coh in ["sepsis", "covid"]:
    for K in ["k1000", "k4000"]:
        name = f"{coh}_{K}"; o = pd.read_csv(f"{U}/{name}_obs.csv.gz"); e = np.load(f"{U}/{name}_emb.npz")
        pre, post = e["X_pca"], e["X_pca_harmony"]
        r = {"n_cells": len(o), "n_samples": int(o.sample_id.nunique())}
        for tag, X in [("pre", pre), ("post", post)]:
            kb, kn = M.kbet(X, o.draw.values); li = M.lisi(X, {"iLISI_draw": o.draw.values, "cLISI_l1": o.celltype_l1_pred.values})
            r[f"kBET_{tag}"] = kb; r[f"kBET_null_{tag}"] = kn; r.update({f"{k}_{tag}": v for k, v in li.items()})
        cl = M.leiden_on(post[:, :10])
        for rr, lab in cl.items():
            r[f"clusters_r{rr}"] = int(lab.max() + 1)
            r[f"ARI_l1_r{rr}"] = ARI(o.celltype_l1_pred, lab); r[f"NMI_l1_r{rr}"] = NMI(o.celltype_l1_pred, lab)
            r[f"ARI_l2_r{rr}"] = ARI(o.aifi_l2, lab); r[f"NMI_l2_r{rr}"] = NMI(o.aifi_l2, lab)
            ct = pd.crosstab(lab, o.draw); r[f"draw_dominated_clusters_r{rr}"] = int(((ct.max(1) / ct.sum(1)) > 0.5).sum())
            np.save(f"{R}/{name}_leiden_r{rr}.npy", lab)
        # UMAPs on the extracted subsample (100 cells / sample / draw)
        rows = np.load(f"{U}/{name}_counts.npz")["rows"]
        for tag, X in [("pre", pre), ("post", post)]:
            a = ad.AnnData(np.zeros((len(rows), 0), np.float32)); a.obsm["X"] = X[rows, :30]
            sc.pp.neighbors(a, use_rep="X", n_neighbors=30, random_state=0); sc.tl.umap(a, random_state=0, min_dist=0.3)
            np.save(f"{R}/{name}_umap_{tag}.npy", a.obsm["X_umap"])
        res[name] = r; print(name, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
# ---- real references
def real(name, pre, post, batch, labels):
    r = {"n_cells": len(pre), "n_batches": int(pd.Series(batch).nunique())}
    for tag, X in [("pre", pre), ("post", post)]:
        li = M.lisi(X, {"iLISI_sample": batch, "cLISI_label": labels}); r.update({f"{k}_{tag}": v for k, v in li.items()})
    res[name] = r; print(name, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
st = store.Store("/home/claude/work/sepsis/store")
so = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False) for i in range(3)], ignore_index=True)
e = np.load("/home/claude/work/sepsis/run_major/embed.npz")
real("sepsis_real", e["pcs"][:, :10], e["harmony"][:, :10], so["sample"].values, so.celltype_major_ye2026.values)
cst = store.Store("/home/claude/work/covid/store"); L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
for comp, f, npc in [("myeloid", "myeloid", 30), ("bcell", "bcell", 30), ("tnk", "tnk", 40)]:
    a = ad.read_h5ad(f"/home/claude/work/covid/{f}.h5ad"); idx = a.obs.cell_index.values
    real(f"covid_real_{comp}", a.obsm["X_pca"][:, :npc], a.obsm["X_rep"][:, :npc], a.obs["orig.ident"].astype(str).values,
         L.cell_type_detailed_recap.values[idx])
json.dump(res, open(f"{R}/metrics_A.json", "w"), indent=1)
