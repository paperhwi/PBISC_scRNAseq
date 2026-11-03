"""Matched 60k original/floor/new-LoRA disease DEG against full-depth real DEG.

Pseudobulk unit: patient (COVID visits averaged), within six major lineages.
Sepsis = adult sepsis vs healthy control; COVID = positive vs negative.
Labels: real author six_type; generated real-A 30-NN-transferred six_tr.
"""
from pathlib import Path
import sys
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
OLD=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/results/deg_tables"
OUT=ROOT/"NewLora/analysis_20261102/disease_deg_matched"
(OUT/"deg_tables").mkdir(parents=True,exist_ok=True)
TYPES=("T cell","B cell","NK cell","Monocyte","DC")
GROUP={"sepsis":{"adult_sepsis":"case","adult_HC":"control"},
       "covid":{"covid_positive":"case","covid_negative":"control"}}
COVID_EXCL={"C126_positive_5","C126_positive_9","C126_positive_11",
            "C141_negative","C58_positive_2"}

def bh(p):
    a=np.asarray(p,float);q=np.full(len(a),np.nan);ok=np.isfinite(a)
    v=a[ok];n=len(v)
    if n==0:return q
    o=np.argsort(v); z=np.minimum.accumulate((v[o]*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n);out[o]=np.minimum(z,1);q[ok]=out
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

def ceiling(coh,genes):
    name={"sepsis":"ye2026_sepsis_draw1_K4000_ceiling.h5ad","covid":"imperial_eskd_draw0_K4000_ceiling.h5ad"}[coh]
    with h5py.File(ROOT/"NewLora"/name) as f:
        sample=h5col(f,"obs/sample_id").astype(str)
        group=h5col(f,"obs/group").astype(str)
        flag=h5col(f,"obs/primary_analysis_eligible").astype(bool)
        keep=(~np.char.startswith(group,"ped")) if coh=="sepsis" else flag
        eligible=np.where(keep)[0];rng=np.random.default_rng(20261024)
        per=int(np.ceil(60000/len(np.unique(sample[eligible]))))
        picked=np.concatenate([rng.choice(eligible[np.asarray(ix)],size=min(per,len(ix)),replace=False)
                               for _,ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
        rows=np.sort(rng.choice(picked,size=min(60000,len(picked)),replace=False))
        symbols=h5col(f,"var/gene_symbol").astype(str)
        pos=pd.Series(np.arange(len(symbols)),index=symbols)
        pos=pos[~pos.index.duplicated()].reindex(genes)
        assert pos.notna().all()
        X=rows_from_csr(f,rows,pos.to_numpy(dtype=int))
    lab=pd.read_csv(ROOT/"NewLora/analysis_20261102"/coh/"ceiling_cell_labels.csv.gz",index_col=0)
    assert X.shape[0]==len(lab) and np.array_equal(lab["sample"].to_numpy(dtype=str),sample[rows])
    return X,lab[["sample","group","six_tr"]].copy()

def aggregate(X,obs,coh):
    obs=obs.copy().reset_index(drop=True)
    obs["condition"]=obs.group.astype(str).map(GROUP[coh])
    if coh=="covid":
        obs.loc[obs["sample"].isin(COVID_EXCL),"condition"]=np.nan
    obs["patient"]=obs["sample"].astype(str) if coh=="sepsis" else obs["sample"].astype(str).str.split("_").str[0]
    out={}
    for celltype in TYPES:
        rows=[]
        sub=obs[(obs.six_tr.astype(str)==celltype)&obs.condition.notna()]
        for sample,idx in sub.groupby("sample").groups.items():
            if len(idx)<20:continue
            block=X[np.asarray(idx,dtype=int)]
            counts=np.asarray(block.sum(axis=0)).ravel()
            total=float(counts.sum())
            if total<=0:continue
            meta=obs.iloc[int(idx[0])]
            rows.append((sample,meta.patient,meta.condition,len(idx),np.log2(counts/total*1e6+1)))
        if not rows:continue
        Y=np.vstack([r[4] for r in rows])
        m=pd.DataFrame({"sample":[r[0] for r in rows],"patient":[r[1] for r in rows],
                        "condition":[r[2] for r in rows],"n_cells":[r[3] for r in rows]})
        # COVID patient may have multiple visits; take the mean of sample logCPM.
        groups=m.groupby(["patient","condition"],sort=True).groups
        Z=np.vstack([Y[np.asarray(ix)].mean(axis=0) for ix in groups.values()])
        keys=list(groups)
        out[celltype]=(Z,pd.DataFrame(keys,columns=["patient","condition"]),m)
    return out

def deg(Z,meta,genes):
    case=Z[meta.condition.to_numpy()=="case"]
    ctrl=Z[meta.condition.to_numpy()=="control"]
    if len(case)<3 or len(ctrl)<3:return None
    logfc=case.mean(0)-ctrl.mean(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        p=stats.ttest_ind(case,ctrl,axis=0,equal_var=False).pvalue
    q=bh(p)
    return pd.DataFrame({"gene":genes,"logFC":logfc,"p":p,"q":q,
                         "mean_case":case.mean(0),"mean_ctrl":ctrl.mean(0),
                         "n_case":len(case),"n_ctrl":len(ctrl),
                         "deg":(q<.05)&(np.abs(logfc)>=1)})

summaries=[];sample_counts=[]
for coh in ("sepsis","covid"):
    print("COHORT",coh,flush=True)
    R=ad.read_h5ad(BASE/"data/big"/f"{coh}_ref.h5ad")
    P=ad.read_h5ad(BASE/"data/big"/f"{coh}_prev4k.h5ad")
    F=ad.read_h5ad(BASE/"data/big"/f"{coh}_floor.h5ad")
    genes=np.array(R.var_names)
    assert np.array_equal(np.array(P.var_names),genes) and np.array_equal(np.array(F.var_names),genes)
    prevtr=pd.read_csv(BASE/"results/big"/f"{coh}_prev4k_transfer.csv.gz",index_col=0)
    assert P.obs.index.equals(prevtr.index)
    prevobs=prevtr[["sample","group","six_tr"]]
    floortr=pd.read_csv(BASE/"results/big"/f"{coh}_floor_transfer.csv.gz",index_col=0)
    assert F.obs.index.equals(floortr.index)
    floorobs=floortr[["sample","group","six_tr"]]
    CX,Cobs=ceiling(coh,genes)
    inputs={"prev4k":(P.X.tocsr(),prevobs),"floor":(F.X.tocsr(),floorobs),"ceiling":(CX,Cobs)}
    panel=pd.read_csv(OLD/f"{coh}_Tcell_real.csv.gz")["gene"].astype(str)
    panel=panel[panel.isin(pd.Index(genes))].reset_index(drop=True)
    pos=pd.Index(genes).get_indexer(panel)
    assert (pos>=0).all()
    aggs={s:aggregate(X[:,pos],obs,coh) for s,(X,obs) in inputs.items()}
    for celltype in TYPES:
        oldpath=OLD/f"{coh}_{celltype.replace(' ','')}_real.csv.gz"
        rt=pd.read_csv(oldpath).set_index("gene").reindex(panel).reset_index()
        assert rt.logFC.notna().all()
        rt.to_csv(OUT/"deg_tables"/f"{coh}_{celltype.replace(' ','')}_real.csv.gz",index=False,compression="gzip")
        tabs={"real":rt}
        for source in ("prev4k","floor","ceiling"):
            if celltype not in aggs[source]:continue
            Z,meta,m=aggs[source][celltype]
            T=deg(Z,meta,panel.to_numpy())
            if T is None:continue
            T.to_csv(OUT/"deg_tables"/f"{coh}_{celltype.replace(' ','')}_{source}.csv.gz",index=False,compression="gzip")
            tabs[source]=T
            sample_counts.append({"cohort":coh,"celltype":celltype,"source":source,
                                  "n_case_patients":int(T.n_case.iloc[0]),"n_control_patients":int(T.n_ctrl.iloc[0]),
                                  "n_case_samples":int((m.condition=="case").sum()),
                                  "n_control_samples":int((m.condition=="control").sum())})
        rt=tabs["real"]
        rank=stats.rankdata(-((rt.mean_case+rt.mean_ctrl)/2).to_numpy(),method="average")
        low=rank>2000
        for source in ("prev4k","floor","ceiling"):
            if source not in tabs:continue
            gt=tabs[source]
            assert np.array_equal(rt.gene.to_numpy(),gt.gene.to_numpy())
            r,g=rt.logFC.to_numpy(),gt.logFC.to_numpy()
            rd,gd=rt.deg.to_numpy(),gt.deg.to_numpy()
            overlap=int((rd&gd).sum())
            def corr(a,b):return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
            summaries.append({"cohort":coh,"celltype":celltype,"source":source,
                              "n_case_real":int(rt.n_case.iloc[0]),"n_control_real":int(rt.n_ctrl.iloc[0]),
                              "n_case_gen":int(gt.n_case.iloc[0]),"n_control_gen":int(gt.n_ctrl.iloc[0]),
                              "n_genes":len(rt),"n_low_expr_genes":int(low.sum()),
                              "logFC_r":corr(r,g),"logFC_r_low_expr":corr(r[low],g[low]),
                              "slope":float(np.polyfit(r,g,1)[0]),
                              "n_DEG_real":int(rd.sum()),"n_DEG_gen":int(gd.sum()),
                              "n_DEG_overlap":overlap,
                              "recall":overlap/rd.sum() if rd.sum() else np.nan,
                              "precision":overlap/gd.sum() if gd.sum() else np.nan,
                              "sign_agree_real_DEG":float((np.sign(r[rd])==np.sign(g[rd])).mean()) if rd.sum() else np.nan,
                              "sign_agree_low_expr":float((np.sign(r[low])==np.sign(g[low])).mean())})
    print("DONE",coh,flush=True)

pd.DataFrame(summaries).to_csv(OUT/"matched_deg_summary.csv",index=False)
pd.DataFrame(sample_counts).to_csv(OUT/"unit_counts.csv",index=False)
print(OUT)
