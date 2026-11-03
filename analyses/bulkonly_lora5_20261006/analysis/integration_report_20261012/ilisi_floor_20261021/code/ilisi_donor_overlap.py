"""COVID: does shared donor identity between generated cells and the reference real cells inflate iLISI? Per-cell iLISI (V0, post-Harmony joint space)
of generated cells, grouped by whether the generating donor also has real cells in the reference half A, in the control half B, or neither."""
import os, sys, numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
from sklearn.neighbors import NearestNeighbors
O = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); NF = os.path.expanduser("~/mnt/PBISC/new_inference_depth256_K4000_floor1e-4_20261021"); coh = "covid"
ref = ad.read_h5ad(f"{O}/data/{coh}_ref.h5ad", backed="r").obs; rb = ad.read_h5ad(f"{O}/data/{coh}_realU.h5ad", backed="r").obs
pA, pB = set(ref.patient), set(rb.patient); s = pd.read_csv(f"{NF}/covid/samples.csv").set_index("sample_id").donor_id; rows = []
for name in ["floor", "prev4k"]:
    J = ad.read_h5ad(f"{O}/data/{coh}_joint_{name}_V0.h5ad"); o = J.obs; Y = J.obsm["X_pca_harmony"]; src = (o.source != "ref").to_numpy().astype(int)
    don = s.reindex(o["sample"]).to_numpy() if name == "floor" else pd.Series(o["sample"].to_numpy()).str.split("_").str[0].to_numpy()
    g = np.where(src == 0, "ref", np.where(pd.Series(don).isin(pA).to_numpy(), "donor in ref half A", np.where(pd.Series(don).isin(pB).to_numpy(), "donor in control half B", "donor not in real")))
    q = np.where(src == 1)[0]; d, ix = NearestNeighbors(n_neighbors=90, n_jobs=2).fit(Y).kneighbors(Y[q]); v = L.lisi_vec(d, ix, src)
    for k in ["donor in ref half A", "donor in control half B", "donor not in real"]:
        m = g[q] == k; rows.append(dict(source=name, group=k, n_cells=int(m.sum()), n_donors=int(pd.Series(don[q][m]).nunique()), iLISI_post_median=float(np.median(v[m])) if m.sum() else np.nan,
                                         own_share=float((src[ix[m]] == 1).mean()) if m.sum() else np.nan))
    print(name, "overall median", round(float(np.median(v)), 3), flush=True)
R = pd.DataFrame(rows); R.to_csv(f"{O}/results/covid_donor_overlap_ilisi.csv", index=False); print(R.round(3).to_string(index=False))
