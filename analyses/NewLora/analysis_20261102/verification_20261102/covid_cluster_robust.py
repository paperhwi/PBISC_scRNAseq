"""COVID DEG sensitivity with patient-clustered standard errors."""
from pathlib import Path
import sys
import gc
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

HERE=Path(__file__).resolve().parent
SCRIPT=HERE.parents[1]/'analyze_disease_deg_matched.py'
ns={'__file__':str(SCRIPT)}
exec(SCRIPT.read_text().split('summaries=[];sample_counts=[]')[0],ns)
OLDROOT=HERE.parents[2]/'bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/code'
sys.path.insert(0,str(OLDROOT))
import deg_lib as D
D.PB=str(HERE.parents[2]/'bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb')+'/'

def robust(Y, patient, condition):
    Y=np.asarray(Y,dtype=float)
    patient=np.asarray(patient,dtype=str)
    condition=np.asarray(condition,dtype=str)
    case=condition=='case'; control=condition=='control'
    n1=case.sum();n0=control.sum();N=len(Y);G=len(np.unique(patient))
    mean1=Y[case].mean(axis=0);mean0=Y[control].mean(axis=0)
    effect=mean1-mean0
    residual=Y-np.where(case[:,None],mean1,mean0)
    influence=residual*np.where(case,1/n1,-1/n0)[:,None]
    cluster=pd.DataFrame({'patient':patient}).groupby('patient').indices
    score=np.vstack([influence[idx].sum(axis=0) for idx in cluster.values()])
    se=np.sqrt((score**2).sum(axis=0)*G/(G-1)*(N-1)/(N-2))
    with np.errstate(divide='ignore',invalid='ignore'):
        p=2*stats.t.sf(np.abs(effect/se),df=G-1)
    p=np.where((se==0)&(effect==0),np.nan,p)
    q=ns['bh'](p)
    return effect,p,q,(q<.05)&(np.abs(effect)>=1)

rows=[]
cohort='covid'
d=D.load('real_covid');des=D.design('real_covid',d)
for ct in ns['TYPES']:
    g=D.cpm_matrix(d,des,ct)
    if g.empty:continue
    Y=g[d['genes']].to_numpy(); meta=pd.DataFrame({'patient':g._u,'condition':g._c.replace({'ctrl':'control'})})
    old=D.deg(g,d['genes'])
    effect,p,q,call=robust(Y,meta.patient,meta.condition)
    assert np.allclose(effect,old.logFC,atol=1e-7)
    both=meta.groupby('patient').condition.nunique()
    rows.append({'source':'real','celltype':ct,'n_case':int((meta.condition=='case').sum()),
                 'n_control':int((meta.condition=='control').sum()),'n_crossover_patients':int((both>1).sum()),
                 'welch_DEG':int(old.deg.sum()),'cluster_robust_DEG':int(call.sum()),
                 'shared_DEG':int((old.deg&call).sum()),
                 'welch_reference_match':int(old.deg.sum())==int(pd.read_csv(ns['OLD']/f'covid_{ct.replace(" ","")}_real.csv.gz').deg.sum())})
del d,des
gc.collect()

real=ad.read_h5ad(ns['BASE']/'data/big/covid_ref.h5ad',backed='r')
genes=np.asarray(real.var_names,dtype=str);real.file.close()
panel=pd.read_csv(ns['OLD']/'covid_Tcell_real.csv.gz').gene.astype(str)
panel=panel[panel.isin(set(genes))].reset_index(drop=True)
pos=pd.Index(genes).get_indexer(panel)
for source in ('prev4k','floor','ceiling'):
    print(source,flush=True)
    if source=='ceiling':
        X,obs=ns['ceiling']('covid',genes)
    else:
        x=ad.read_h5ad(ns['BASE']/'data/big'/f'covid_{source}.h5ad')
        X=x.X.tocsr()
        obs=pd.read_csv(ns['BASE']/'results/big'/f'covid_{source}_transfer.csv.gz',index_col=0)[['sample','group','six_tr']]
        assert x.obs.index.equals(obs.index)
        del x
    aggs=ns['aggregate'](X[:,pos],obs,'covid')
    for ct,(Y,meta,_) in aggs.items():
        old=ns['deg'](Y,meta,panel.to_numpy())
        if old is None:continue
        effect,p,q,call=robust(Y,meta.patient,meta.condition)
        assert np.allclose(effect,old.logFC,atol=1e-4), (source,ct,float(np.max(np.abs(effect-old.logFC))))
        both=meta.groupby('patient').condition.nunique()
        rows.append({'source':source,'celltype':ct,'n_case':int((meta.condition=='case').sum()),
                     'n_control':int((meta.condition=='control').sum()),'n_crossover_patients':int((both>1).sum()),
                     'welch_DEG':int(old.deg.sum()),'cluster_robust_DEG':int(call.sum()),
                     'shared_DEG':int((old.deg&call).sum()),'welch_reference_match':True})
    del X,obs,aggs
    gc.collect()

pd.DataFrame(rows).to_csv(HERE/'covid_cluster_robust_summary.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
