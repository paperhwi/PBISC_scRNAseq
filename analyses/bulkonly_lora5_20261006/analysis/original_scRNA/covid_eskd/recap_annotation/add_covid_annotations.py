"""Add the Stephenson-2025 recapitulation annotations to a copy of eskd_covid19.h5ad (in place on the copy).
usage: python add_covid_annotations.py <copy.h5ad> <annotation_dir>"""
import sys, os, json, h5py, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import h5ad_lite as hl
path, ann = sys.argv[1], sys.argv[2]
L = pd.read_csv(f"{ann}/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
E = np.load(f"{ann}/covid_obsm.npz")
with h5py.File(path, "a") as f:
    idx = np.array([x.decode() if isinstance(x, bytes) else x for x in f["obs/_index"][:]])
    assert len(idx) == len(L) and (idx == L.index.astype(str).to_numpy()).all(), "cell order mismatch"
    L["isg_score"] = pd.to_numeric(L["isg_score"].replace("", np.nan))
    L["asc_gmm"] = L["asc_gmm"].astype(str).eq("True"); L["dex_mono_recap"] = L["dex_mono_recap"].astype(str).eq("True")
    hl.add_obs_columns(f, L)
    om = f["obsm"]
    for k in E.files:
        if k in om: del om[k]
        hl.write_array(om, k, E[k])
    u = f["uns"]
    notes = json.load(open(f"{ann}/covid_uns.json"))
    for k, v in notes.items():
        if k in u: del u[k]
        d = u.create_dataset(k, data=str(v), dtype=hl.STR); hl._attrs(d, "string", "0.2.0")
print("annotated", path, L.shape)
