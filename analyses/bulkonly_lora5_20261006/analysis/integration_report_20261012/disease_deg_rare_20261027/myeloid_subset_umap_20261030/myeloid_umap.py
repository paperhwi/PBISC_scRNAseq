import numpy as np, pandas as pd, anndata as ad, umap
B="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big_joint/"
for coh in ["sepsis","covid"]:
    a=ad.read_h5ad(f"{B}{coh}_joint_floor.h5ad")
    o=a.obs.copy()
    six=np.where(o.source=="ref",o.six_type.astype(str),o.six_tr.astype(str))
    m=six=="Monocyte"
    E=a.obsm["X_pca_harmony"][m]
    U=umap.UMAP(n_neighbors=30,min_dist=0.3,random_state=0).fit_transform(E)
    o=o[m].copy(); o["u1"]=U[:,0]; o["u2"]=U[:,1]
    o["sub"]=np.where(o.source=="ref",o.fine.astype(str),o.fine_tr.astype(str))
    o.to_pickle(f"myeloid/{coh}_myeloid_umap.pkl")
    print(coh,o.shape); print(pd.crosstab(o["sub"],o.source))
