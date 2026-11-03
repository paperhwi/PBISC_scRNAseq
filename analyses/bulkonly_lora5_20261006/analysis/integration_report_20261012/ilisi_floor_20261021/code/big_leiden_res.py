"""Cluster matching at three Leiden resolutions (0.5, 1, 2.5; resolution 1 = big_joint.py) on the saved joint Harmony embeddings.
Usage: python big_leiden_res.py <coh> <source>. Output results/big/<coh>_<source>_leiden_res.csv"""
import sys, numpy as np, pandas as pd, anndata as ad, scanpy as sc
coh, s = sys.argv[1], sys.argv[2]; B = "/mnt/user-data/outputs/ilisi_floor_20261021"
J = ad.read_h5ad(f"{B}/data/big_joint/{coh}_joint_{s}.h5ad"); src = (J.obs.source.to_numpy() != "ref").astype(int)
sc.pp.neighbors(J, n_neighbors=15, use_rep="X_pca_harmony", random_state=0); out = []
for res in (0.5, 1.0, 2.5):
    sc.tl.leiden(J, resolution=res, flavor="igraph", n_iterations=2, random_state=0, key_added="l"); lc = J.obs.l.astype(int).to_numpy()
    ct = pd.crosstab(lc, src); ct.columns = ["ref", "src"]; fr = ct / ct.sum(0); w = fr.src / (fr.src + fr.ref)
    js = 0.5 * (np.nansum(fr.ref * np.log2(np.where(fr.ref > 0, 2 * fr.ref / (fr.ref + fr.src), 1))) + np.nansum(fr.src * np.log2(np.where(fr.src > 0, 2 * fr.src / (fr.ref + fr.src), 1))))
    out.append({"cohort": coh, "source": s, "resolution": res, "n_clusters": len(ct), "src_in_shared": fr.src[(w >= 0.2) & (w <= 0.8)].sum(), "ref_in_shared": fr.ref[(w >= 0.2) & (w <= 0.8)].sum(),
                "n_shared": int(((w >= 0.2) & (w <= 0.8)).sum()), "n_src_dominated": int((w > 0.8).sum()), "n_ref_dominated": int((w < 0.2).sum()),
                "ref_in_ref_dominated": fr.ref[w < 0.2].sum(), "src_in_src_dominated": fr.src[w > 0.8].sum(), "freq_pearson": np.corrcoef(fr.ref, fr.src)[0, 1], "JSD": js})
pd.DataFrame(out).to_csv(f"{B}/results/big/{coh}_{s}_leiden_res.csv", index=False); print(pd.DataFrame(out).round(3).to_string())
