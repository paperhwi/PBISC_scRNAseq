import numpy as np, pandas as pd, anndata as ad, umap, pickle, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
B="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big_joint/"
BLUE,RED,GREY,INK="#2B83BA","#D7191C","#7F7F7F","#1F2933"
rng=np.random.default_rng(7)
out={}
for coh in ["sepsis","covid"]:
    a=ad.read_h5ad(f"{B}{coh}_joint_floor.h5ad")
    n=a.n_obs; idx=rng.choice(n,min(40000,n),replace=False); idx.sort()
    E=a.obsm["X_pca_harmony"][idx]
    U=umap.UMAP(n_neighbors=30,min_dist=0.3,random_state=0).fit_transform(E)
    o=a.obs.iloc[idx].copy(); o["u1"]=U[:,0]; o["u2"]=U[:,1]
    o.to_pickle(f"figs/{coh}_umap.pkl")
    print(coh,"done",o.shape)
