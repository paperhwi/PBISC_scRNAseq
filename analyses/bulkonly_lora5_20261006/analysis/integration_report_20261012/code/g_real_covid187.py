"""Real-reference LISI for the COVID compartments restricted to the published 187-sample roster (replaces the
covid_real_* entries of metrics_A.json, which used all 198 local samples)."""
import json, sys, numpy as np, pandas as pd, h5py
sys.path.insert(0, "/home/claude/work/report"); import metrics as M
R = "/home/claude/work/report"
roster = set(pd.read_csv(f"{R}/covid_roster187_local_ids.csv").sample_id.astype(str))
L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
A = json.load(open(f"{R}/metrics_A.json"))
sys.path.insert(0, "/home/claude/work"); import store
SID = store.Store("/home/claude/work/covid/store").obs.sample_id.astype(str).to_numpy()
for comp, npc in [("myeloid", 30), ("bcell", 30), ("tnk", 40)]:
    with h5py.File(f"/home/claude/work/covid/{comp}.h5ad") as h:
        idx = h["obs/cell_index"][:]; pca = h["obsm/X_pca"][:, :npc]; rep = h["obsm/X_rep"][:, :npc]
        g = h["obs/orig.ident"]
        oi = (g["categories"][:][g["codes"][:]] if isinstance(g, h5py.Group) else g[:]).astype(str)
    keep = pd.Series(SID[idx]).isin(roster).to_numpy()          # roster by biological sample; batch = sequencing library (orig.ident)
    lab = L.cell_type_detailed_recap.values[idx][keep]; b = oi[keep]
    r = {"n_cells": int(keep.sum()), "n_batches": int(pd.Series(b).nunique())}
    for tag, X in [("pre", pca[keep]), ("post", rep[keep])]:
        li = M.lisi(X, {"iLISI_sample": b, "cLISI_label": lab}); r.update({f"{k}_{tag}": v for k, v in li.items()})
    A[f"covid_real_{comp}"] = r; print(comp, json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in r.items()}), flush=True)
json.dump(A, open(f"{R}/metrics_A.json", "w"), indent=1)
