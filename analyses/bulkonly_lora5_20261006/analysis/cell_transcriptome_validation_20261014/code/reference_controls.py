"""Metric calibration on real cells only; these are not PBISC validation results."""
import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1]
PREV=ROOT/'bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014'
sys.path.insert(0,str(PREV/'code'))
from common import read_csr_rows,norm,boot_mean
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy.stats import wasserstein_distance
from scipy.spatial.distance import cdist
rng=np.random.default_rng(20261014);rows_out=[]
def metrics(a,b):
    sd=np.maximum(a.std(0),.1);az=(a-a.mean(0))/sd;bz=(b-a.mean(0))/sd
    marginal=np.mean([wasserstein_distance(az[:,j],bz[:,j]) for j in range(a.shape[1])])
    def correlation(q):
        q=np.asarray(q,dtype=np.float64)
        centered=q-q.mean(0);std=q.std(0);den=np.outer(std,std)
        return np.divide(centered.T@centered/len(q),den,out=np.zeros_like(den),where=den>1e-12)
    ca=correlation(a);cb=correlation(b)
    off=~np.eye(a.shape[1],dtype=bool)
    coexpr=np.abs(ca-cb)[off].mean()
    energy=2*cdist(az,bz).mean()-cdist(az,az).mean()-cdist(bz,bz).mean()
    covariance=np.abs((az-az.mean(0)).T@(az-az.mean(0))/len(az)-(bz-bz.mean(0)).T@(bz-bz.mean(0))/len(bz))[off].mean()
    return dict(marginal_W1=marginal,corr_MAE=coexpr,covariance_MAE=covariance,energy_distance=energy)
for cohort in ['ye2026_sepsis','imperial_eskd']:
    path=OUT/f'inputs/{cohort}_real_heldout_panel_counts.h5ad'
    with h5py.File(path,'r') as f:md=ad.io.read_elem(f['obs'])
    # One alphabetically first visit per patient; up to four largest major types.
    first=md[['patient_id','sample_id']].drop_duplicates().sort_values('sample_id').drop_duplicates('patient_id')
    chosen=[];groups=[]
    for s in first.sample_id:
        counts=md.loc[md.sample_id==s].major.value_counts()
        for typ in counts[counts>=160].head(4).index:
            ix=np.where((md.sample_id==s)&(md.major==typ))[0]
            ix=np.sort(rng.choice(ix,min(300,len(ix)),replace=False));chosen.extend(ix);groups.append((s,typ,ix))
    chosen=np.sort(np.asarray(chosen));x=read_csr_rows(path,chosen);y=norm(x)
    for s,typ,ix in groups:
        yy=y[np.searchsorted(chosen,ix)].toarray();order=rng.permutation(len(yy));n=len(order)//2
        a0,b0=yy[order[:n]],yy[order[n:n*2]]
        cols=np.argsort(a0.var(0))[-50:];a,b=a0[:,cols],b0[:,cols]
        sh=b.copy()
        for j in range(sh.shape[1]):rng.shuffle(sh[:,j])
        assert np.array_equal(np.sort(b,axis=0),np.sort(sh,axis=0))
        collapsed=np.repeat(b.mean(0)[None,:],len(b),axis=0)
        for ctrl,q in [('real_split_half',b),('independent_gene_shuffle',sh),('centroid_collapse',collapsed)]:
            with np.errstate(divide='ignore',invalid='ignore'): vals=metrics(a,q)
            rows_out.append(dict(cohort=cohort,patient_id=md.iloc[ix[0]].patient_id,sample_id=s,type=typ,n_each=n,control=ctrl,**vals))
        with np.errstate(divide='ignore',invalid='ignore'): vals=metrics(b,sh)
        rows_out.append(dict(cohort=cohort,patient_id=md.iloc[ix[0]].patient_id,sample_id=s,type=typ,n_each=n,control='shuffle_vs_its_own_parent',**vals))
    print('Calibrated',cohort,len(groups),'sample-type groups',flush=True)
d=pd.DataFrame(rows_out);d.to_csv(OUT/'results/reference_metric_controls.csv',index=False)
summ=[]
for (cohort,ctrl),z in d.groupby(['cohort','control']):
    for metric in ['marginal_W1','corr_MAE','covariance_MAE','energy_distance']:
        v=z.groupby('patient_id')[metric].mean();est,lo,hi=boot_mean(v,rng)
        summ.append(dict(cohort=cohort,control=ctrl,metric=metric,n_patients=len(v),mean=est,ci_lo=lo,ci_hi=hi))
pd.DataFrame(summ).to_csv(OUT/'results/reference_metric_control_summary.csv',index=False)
delta=[]
for cohort,z in d[d.control.isin(['real_split_half','independent_gene_shuffle','centroid_collapse'])].groupby('cohort'):
    for metric in ['marginal_W1','corr_MAE','covariance_MAE','energy_distance']:
        w=z.pivot(index=['patient_id','sample_id','type'],columns='control',values=metric)
        for ctrl in ['independent_gene_shuffle','centroid_collapse']:
            v=(w[ctrl]-w.real_split_half).groupby('patient_id').mean();est,lo,hi=boot_mean(v,rng)
            delta.append(dict(cohort=cohort,control=ctrl,metric=metric,n_patients=len(v),delta=est,ci_lo=lo,ci_hi=hi))
pd.DataFrame(delta).to_csv(OUT/'results/reference_metric_control_deltas.csv',index=False)
assert (d[d.control=='shuffle_vs_its_own_parent'].marginal_W1==0).all()
json.dump(dict(scope='real-only assay calibration',seed=20261014,genes='50 top variance genes selected separately in real half A; preliminary metric QC, not locked PBISC primary gene sets',normalization='log1p CP10000 on measured panel',patient_weighting='sample-type means within patient, then equal patient weights',n_bootstrap=2000,PBISC_evaluated=False,exact_marginals_preserved=True),open(OUT/'results/reference_control_methods.json','w'),indent=2)
