"""Old S3 realism protocol applied to the matched 60k new-LoRA sample.

Usage: py -3.12 new_lora_realism.py sepsis|covid
"""
from pathlib import Path
import sys
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
sys.path.insert(0, str(BASE / "code"))
import ilisi_lib as L

COH = sys.argv[1]
assert COH in ("sepsis", "covid")
DC_ONLY = "--dc-only" in sys.argv
OUT = Path(__file__).resolve().parent
INPUT = ROOT / "NewLora" / ("ye2026_sepsis_draw1_K4000_ceiling.h5ad" if COH == "sepsis" else "imperial_eskd_draw0_K4000_ceiling.h5ad")
T0 = time.time()


def log(*args):
    print(f"[{time.time()-T0:7.1f}s]", COH, *args, flush=True)


def column(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object)
        return cats[g["codes"][:]]
    arr = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object) if arr.dtype.kind in "SO" else arr


def read_csr_rows(f, rows, cols, step=5000):
    g = f["X"]
    ptr = g["indptr"][:]
    ncol = int(g.attrs["shape"][1])
    blocks = []
    for start in range(0, int(rows.max()) + 1, step):
        selected = rows[(rows >= start) & (rows < start + step)]
        if len(selected) == 0:
            continue
        end = min(start + step, len(ptr)-1)
        a, b = int(ptr[start]), int(ptr[end])
        block = sp.csr_matrix((g["data"][a:b].astype(np.float32), g["indices"][a:b], ptr[start:end+1]-a), shape=(end-start,ncol))
        blocks.append(block[selected-start][:,cols])
    return sp.vstack(blocks, format="csr")


def dc(real, gen, k=5):
    nn = NearestNeighbors(n_neighbors=k+1).fit(real)
    radii = nn.kneighbors(real)[0][:,-1]
    dist, ix = NearestNeighbors(n_neighbors=k).fit(real).kneighbors(gen)
    density = (dist <= radii[ix]).sum()/(k*len(gen))
    nearest_gen = NearestNeighbors(n_neighbors=1).fit(gen).kneighbors(real)[0][:,0]
    coverage = (nearest_gen <= radii).mean()
    return float(density), float(coverage)


ref = ad.read_h5ad(BASE / "data/big" / f"{COH}_ref.h5ad")
labels = pd.read_csv(ROOT / "NewLora/analysis_20261102" / COH / "ceiling_cell_labels.csv.gz", keep_default_na=False)
hv = L.hvg_seurat_v3([ref.X.tocsr()])
genes = np.asarray(ref.var_names)
log("real A and labels loaded; HVGs", len(hv))

with h5py.File(INPUT) as f:
    sample = column(f, "obs/sample_id").astype(str)
    group = column(f, "obs/group").astype(str)
    eligible_flag = column(f, "obs/primary_analysis_eligible").astype(bool)
    keep = (~np.char.startswith(group,"ped")) if COH == "sepsis" else eligible_flag
    eligible = np.where(keep)[0]
    rng = np.random.default_rng(20261024)
    per_sample = int(np.ceil(60000/len(np.unique(sample[eligible]))))
    picked = np.concatenate([rng.choice(eligible[np.asarray(ix)], size=min(per_sample,len(ix)), replace=False)
                             for _,ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
    rows = np.sort(rng.choice(picked,size=min(60000,len(picked)),replace=False))
    assert np.array_equal(sample[rows], labels["sample"].astype(str).to_numpy())
    assert np.array_equal(group[rows], labels["group"].astype(str).to_numpy())
    symbols = column(f,"var/gene_symbol").astype(str)
    positions = pd.Series(np.arange(len(symbols)),index=symbols)
    positions = positions[~positions.index.duplicated()]
    cols = positions.reindex(genes)
    assert cols.notna().all()
    gen_counts = read_csr_rows(f, rows, cols.to_numpy(dtype=int))
log("new LoRA sampled",gen_counts.shape)

YR = L.lognorm(ref.X.tocsr())[:,hv]
YS = L.lognorm(gen_counts)[:,hv]
A = YR.toarray().astype(np.float32)
mu,sd = A.mean(0),A.std(0)+1e-6
A = np.clip((A-mu)/sd,-10,10)
if not DC_ONLY:
    pca = PCA(50,random_state=0).fit(A)
    ZR = pca.transform(A).astype(np.float32)
pca30 = PCA(30,random_state=0).fit(A)
ZR30 = pca30.transform(A).astype(np.float32)
del A,YR,gen_counts
B = YS.toarray().astype(np.float32)
B = np.clip((B-mu)/sd,-10,10)
if not DC_ONLY:
    ZS = pca.transform(B).astype(np.float32)
ZS30 = pca30.transform(B).astype(np.float32)
del B,YS
log("PCA complete", "density/coverage only" if DC_ONLY else "C2ST and density/coverage")

lin_map = {"T cell":"T","NK cell":"NK","B cell":"B","Monocyte":"Myeloid","DC":"Myeloid","Other":"Other"}
lr = np.array([lin_map.get(x,"Other") for x in ref.obs.six_type.astype(str)])
ls = np.array([lin_map.get(x,"Other") for x in labels.six_tr.astype(str)])
sr = ref.obs["sample"].astype(str).to_numpy()
ss = labels["sample"].astype(str).to_numpy()
rng = np.random.default_rng(0)
crows = []
for lin in ("T","NK","B","Myeloid"):
    ir,js = np.where(lr==lin)[0],np.where(ls==lin)[0]
    n = min(4000,len(ir),len(js))
    a,b = rng.choice(ir,n,replace=False),rng.choice(js,n,replace=False)
    if DC_ONLY:
        continue
    features = np.vstack([ZR[a],ZS[b]])
    y = np.r_[np.zeros(n),np.ones(n)]
    patient = np.r_["ref:"+sr[a],"ceiling:"+ss[b]]
    prob = np.zeros(2*n)
    for fold,(tr,te) in enumerate(GroupKFold(5).split(features,y,patient),1):
        model = HistGradientBoostingClassifier(max_iter=200,random_state=0)
        prob[te] = model.fit(features[tr],y[tr]).predict_proba(features[te])[:,1]
    auc = roc_auc_score(y,prob)
    crows.append(dict(cohort=COH,lineage=lin,source="ceiling",n_per_class=n,c2st_auc=auc))
    pd.DataFrame(crows).to_csv(OUT/f"{COH}_new_lora_c2st.csv",index=False)
    log("C2ST",lin,"n",n,"AUROC",round(auc,4))

dr=[]
for lin in ("all","T","NK","B","Myeloid"):
    ir = np.arange(len(lr)) if lin=="all" else np.where(lr==lin)[0]
    js = np.arange(len(ls)) if lin=="all" else np.where(ls==lin)[0]
    n = min(10000,len(ir),len(js))
    a,b = rng.choice(ir,n,replace=False),rng.choice(js,n,replace=False)
    density,coverage = dc(ZR30[a],ZS30[b])
    dr.append(dict(cohort=COH,source="ceiling",space="PCA (no correction)",lineage=lin,n=n,density=density,coverage=coverage))
    pd.DataFrame(dr).to_csv(OUT/f"{COH}_new_lora_density_coverage.csv",index=False)
    log("density/coverage",lin,n,round(density,4),round(coverage,4))
