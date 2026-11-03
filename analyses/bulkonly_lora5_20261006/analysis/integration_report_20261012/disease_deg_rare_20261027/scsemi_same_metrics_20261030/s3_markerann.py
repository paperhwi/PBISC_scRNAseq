"""Top10-marker de novo annotation of scSemiProfiler cells (same procedure as markerann/s2_annotate.py), then share error and identity r vs real A."""
import sys; sys.path.insert(0,"/home/claude/pbisc/markerann")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, warnings; warnings.filterwarnings("ignore")
from annot_lib import annotate
SS="/mnt/user-data/uploads/integration_report_20261012/three_model_20261022/data/"
for coh in ["sepsis","covid"]:
    a=ad.read_h5ad(f"{SS}{coh}_scsemi.h5ad"); a.obs_names_make_unique()
    tr=pd.read_csv(f"/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/results/big/{coh}_scsemi_transfer.csv.gz",index_col=0,keep_default_na=False)
    a.obs["six_knn"]=tr.six_tr.values; a.obs["fine_knn"]=tr.fine_tr.values
    sc.pp.normalize_total(a,1e4); sc.pp.log1p(a)
    a=annotate(a,coh)
    cols=["sample","patient","group","six_type","aifi_l2","six_knn","fine_knn","cl1","cl2","six_mk_cl","fine_mk_cl","six_mk_cell","fine_mk_cell"]
    a.obs[cols].to_csv(f"scsemi/out/{coh}_scsemi_markerann_labels.csv.gz")
    a.write_h5ad(f"scsemi/out/{coh}_scsemi_annot.h5ad",compression="gzip")
    print(coh,"L1",a.obs.cl1.nunique(),"L2",a.obs.cl2.nunique(),flush=True)
