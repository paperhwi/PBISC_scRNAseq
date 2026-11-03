"""Exploratory disease pseudobulk DEG for the feasible Sepsis rare subtype.

Uses the previously fixed 60k-cell subsets. Minimum 10 subtype cells per
sample/patient; real author fine labels and generated 30-NN transferred labels.
"""
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import anndata as ad
import warnings

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
OLD=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/results/deg_tables"
OUT=ROOT/"NewLora/analysis_20261102/rare_disease_deg"
OUT.mkdir(parents=True,exist_ok=True)
SUB="Mono_c01-CD14-HLA-DRA"
MIN_CELLS=10

def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p)
    v=p[ok];n=len(v)
    if n==0:return q
    o=np.argsort(v);z=np.minimum.accumulate((v[o]*n/np.arange(1,n+1))[::-1])[::-1]
    a=np.empty(n);a[o]=np.minimum(z,1);q[ok]=a
    return q

def h5col(f,path):
    g=f[path]
    if isinstance(g,h5py.Group):
        cats=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in g["categories"][:]],dtype=object)
        return cats[g["codes"][:]]
    a=g[:]
    return np.array([x.decode() if isinstance(x,bytes) else x for x in a],dtype=object) if a.dtype.kind in "SO" else a

def rows_from_csr(f,rows,columns,step=20000):
    g=f["X"];ptr=g["indptr"][:];ncol=int(g.attrs["shape"][1]);blocks=[]
    for start in range(0,int(rows.max())+1,step):
        chosen=rows[(rows>=start)&(rows<start+step)]
        if len(chosen)==0:continue
        end=min(start+step,len(ptr)-1);a,b=int(ptr[start]),int(ptr[end])
        block=sp.csr_matrix((g["data"][a:b].astype(np.float32),g["indices"][a:b],ptr[start:end+1]-a),shape=(end-start,ncol))
        blocks.append(block[chosen-start][:,columns])
    return sp.vstack(blocks,format="csr")

def new_lora(panel):
    path=ROOT/"NewLora/ye2026_sepsis_draw1_K4000_ceiling.h5ad"
    with h5py.File(path) as f:
        sample=h5col(f,"obs/sample_id").astype(str)
        group=h5col(f,"obs/group").astype(str)
        eligible=np.where(~np.char.startswith(group,"ped"))[0]
        rng=np.random.default_rng(20261024)
        per=int(np.ceil(60000/len(np.unique(sample[eligible]))))
        picked=np.concatenate([rng.choice(eligible[np.asarray(ix)],size=min(per,len(ix)),replace=False)
                               for _,ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
        rows=np.sort(rng.choice(picked,size=min(60000,len(picked)),replace=False))
        lab=pd.read_csv(ROOT/"NewLora/analysis_20261102/sepsis/ceiling_cell_labels.csv.gz",index_col=0)
        assert len(rows)==len(lab) and np.array_equal(lab["sample"].astype(str).values,sample[rows])
        mask=lab.fine_tr.astype(str).to_numpy()==SUB
        symbols=h5col(f,"var/gene_symbol").astype(str)
        pos=pd.Series(np.arange(len(symbols)),index=symbols)
        cols=pos[~pos.index.duplicated()].reindex(panel)
        assert cols.notna().all()
        X=rows_from_csr(f,rows[mask],cols.to_numpy(dtype=int))
    return X,lab.loc[mask,["sample","group"]].reset_index(drop=True)

def aggregate(X,obs,genes):
    obs=obs.reset_index(drop=True).copy()
    obs["condition"]=obs.group.astype(str).map({"adult_sepsis":"sepsis","adult_HC":"HC","adult_ICU_ctrl":"ICU"})
    kept=[];arrays=[]
    for sample,ix in obs[obs.condition.notna()].groupby("sample").groups.items():
        if len(ix)<MIN_CELLS:continue
        counts=np.asarray(X[np.asarray(ix,dtype=int)].sum(axis=0)).ravel()
        if counts.sum()==0:continue
        arrays.append(np.log2(counts/counts.sum()*1e6+1))
        kept.append((sample,obs.iloc[int(ix[0])].condition,len(ix)))
    Y=np.vstack(arrays)
    meta=pd.DataFrame(kept,columns=["patient","condition","n_cells"])
    return Y,meta

def deg(Y,meta,genes,case,ctrl="HC"):
    a=Y[meta.condition.eq(case).to_numpy()];b=Y[meta.condition.eq(ctrl).to_numpy()]
    assert len(a)>=3 and len(b)>=3,(case,len(a),len(b))
    lf=a.mean(0)-b.mean(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        p=stats.ttest_ind(a,b,axis=0,equal_var=False).pvalue
    q=bh(p)
    return pd.DataFrame({"gene":genes,"logFC":lf,"p":p,"q":q,
                         "mean_case":a.mean(0),"mean_ctrl":b.mean(0),
                         "n_case":len(a),"n_ctrl":len(b),
                         "deg":(q<.05)&(np.abs(lf)>=1)})

R=ad.read_h5ad(BASE/"data/big/sepsis_ref.h5ad")
U=ad.read_h5ad(BASE/"data/big/sepsis_realU.h5ad")
genes=np.array(R.var_names)
panel=pd.read_csv(OLD/"sepsis_Tcell_real.csv.gz")["gene"].astype(str)
panel=panel[panel.isin(pd.Index(genes))].reset_index(drop=True)
cols=pd.Index(genes).get_indexer(panel)
assert (cols>=0).all()
realX=sp.vstack([R.X,U.X],format="csr")
realobs=pd.concat([R.obs,U.obs],ignore_index=True)
inputs={}
mask=realobs.fine.astype(str).eq(SUB).to_numpy()
inputs["real"]=(realX[mask][:,cols],realobs.loc[mask,["sample","group"]])
for src in ("prev4k","floor"):
    A=ad.read_h5ad(BASE/"data/big"/f"sepsis_{src}.h5ad")
    tr=pd.read_csv(BASE/"results/big"/f"sepsis_{src}_transfer.csv.gz",index_col=0)
    assert A.obs.index.equals(tr.index)
    mask=tr.fine_tr.astype(str).eq(SUB).to_numpy()
    inputs[src]=(A.X.tocsr()[mask][:,cols],tr.loc[mask,["sample","group"]])
inputs["ceiling"]=new_lora(panel)

unit=[];summ=[]
for src,(X,obs) in inputs.items():
    Y,meta=aggregate(X,obs,panel.to_numpy())
    meta.assign(source=src).to_csv(OUT/f"mono_c01_{src}_patient_units.csv",index=False)
    unit.append(meta.assign(source=src))
    for case in ("sepsis","ICU"):
        T=deg(Y,meta,panel.to_numpy(),case)
        T.to_csv(OUT/f"mono_c01_{case}_vs_HC_{src}.csv.gz",index=False,compression="gzip")
        summ.append({"source":src,"contrast":f"{case}_vs_HC","n_case":int(T.n_case.iloc[0]),
                     "n_HC":int(T.n_ctrl.iloc[0]),"n_DEG":int(T.deg.sum())})
summary=pd.DataFrame(summ)
for case in ("sepsis","ICU"):
    rt=pd.read_csv(OUT/f"mono_c01_{case}_vs_HC_real.csv.gz")
    for src in ("prev4k","floor","ceiling"):
        gt=pd.read_csv(OUT/f"mono_c01_{case}_vs_HC_{src}.csv.gz")
        x,y=rt.logFC.to_numpy(),gt.logFC.to_numpy()
        a,b=rt.deg.to_numpy(),gt.deg.to_numpy()
        mask=(summary.source.eq(src)&summary.contrast.eq(f"{case}_vs_HC"))
        summary.loc[mask,"logFC_r"]=float(np.corrcoef(x,y)[0,1])
        summary.loc[mask,"real_n_DEG"]=int(a.sum())
        summary.loc[mask,"n_overlap"]=int((a&b).sum())
        summary.loc[mask,"recall"]=(a&b).sum()/a.sum() if a.sum() else np.nan
        summary.loc[mask,"precision"]=(a&b).sum()/b.sum() if b.sum() else np.nan
        same=(a&b&(np.sign(x)==np.sign(y)))
        summary.loc[mask,"n_overlap_same_direction"]=int(same.sum())
        summary.loc[mask,"recall_same_direction"]=same.sum()/a.sum() if a.sum() else np.nan
summary.to_csv(OUT/"mono_c01_deg_summary.csv",index=False)
print(summary.round(3).to_string(index=False))
