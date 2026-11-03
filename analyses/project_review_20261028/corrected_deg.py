"""Recompute PBISC floor DEG concordance with patient-cluster resampling.

The 2026-10-27 CSVs are left intact. Run with the local analysis runtime.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path('F:/PBISC')
BASE = ROOT/'bulkonly_lora5_20261006/analysis/integration_report_20261012'
sys.path.insert(0, str(BASE/'disease_deg_rare_20261027/code'))
import deg_lib as old

old.PB = str(BASE/'ilisi_floor_20261021/data/pb') + '/'
OUT = ROOT/'project_review_20261028/corrected_deg'
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20261028)
B = 200
TYPES = ['T cell','B cell','NK cell','Monocyte','DC']

def bh(p):
    p=np.asarray(p,float)
    ok=np.isfinite(p)
    out=np.full(len(p),np.nan)
    v=p[ok]
    if not len(v): return out
    ix=np.argsort(v); q=np.minimum.accumulate((v[ix]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]
    tmp=np.empty(len(v)); tmp[ix]=np.minimum(q,1)
    out[ok]=tmp
    return out

def patient_cluster_deg(y,cond,cluster):
    """OLS case-minus-control with CR1 patient-cluster t inference per gene."""
    n,gene_count=y.shape
    groups,gi=np.unique(cluster,return_inverse=True)
    G=len(groups)
    if G<4 or min((cond==0).sum(),(cond==1).sum())<3:
        return None
    x=np.column_stack([np.ones(n),cond])
    inv=np.linalg.inv(x.T@x)
    beta=inv@x.T@y
    resid=y-x@beta
    z=np.zeros((G,gene_count))
    leverage=inv[1]@x.T
    np.add.at(z,gi,resid*leverage[:,None])
    correction=G/(G-1)*(n-1)/(n-2)
    se=np.sqrt(correction*np.square(z).sum(0))
    with np.errstate(divide='ignore',invalid='ignore'):
        t=np.divide(beta[1],se,out=np.zeros(gene_count),where=se>0)
    p=2*stats.t.sf(abs(t),G-1)
    p[(se==0)&(beta[1]==0)]=1
    p[(se==0)&(beta[1]!=0)]=0
    q=bh(p)
    return beta[1],(q<0.05)&(abs(beta[1])>=1),q

def welch_deg(y,cond):
    if min((cond==0).sum(),(cond==1).sum())<3: return None
    a=y[cond==1]; b=y[cond==0]
    lfc=a.mean(0)-b.mean(0)
    with np.errstate(invalid='ignore'):
        _,p=stats.ttest_ind(a,b,equal_var=False,axis=0)
    q=bh(p)
    return lfc,(q<0.05)&(abs(lfc)>=1),q

def indices_by_patient(g):
    return {u:np.flatnonzero(g._u.astype(str).to_numpy()==u) for u in g._u.astype(str).unique()}

def resample(index,strata=None):
    """Each sampled patient contributes all visits and both conditions; repeats remain repeats."""
    if strata is not None:
        draw=np.concatenate([RNG.choice(u,size=len(u),replace=True) for u in strata])
    else:
        units=np.array(list(index))
        draw=RNG.choice(units,size=len(units),replace=True)
    rows=np.concatenate([index[u] for u in draw])
    cluster=np.concatenate([np.repeat(i,len(index[u])) for i,u in enumerate(draw)])
    return rows,cluster

def metrics(a,b):
    al,ad,_=a; bl,bd,_=b
    good=np.isfinite(al)&np.isfinite(bl)
    r=np.corrcoef(al[good],bl[good])[0,1]
    intersect=(ad&bd).sum()
    rec=intersect/ad.sum() if ad.sum() else np.nan
    pre=intersect/bd.sum() if bd.sum() else np.nan
    sign_all=(np.sign(al[good])==np.sign(bl[good])).mean()
    sign_deg=(np.sign(al[ad])==np.sign(bl[ad])).mean() if ad.sum() else np.nan
    return np.array([r,rec,pre,sign_all,sign_deg],float)

def prepare(coh,roster):
    source={s:old.load(f'{s}_{coh}') for s in ('real','floor')}
    des={s:old.design(f'{s}_{coh}',source[s]) for s in source}
    if coh=='covid' and roster=='published_187':
        x=pd.read_csv(ROOT/'paper_claim_analysis/covid_sample_patient_id_reconciliation.csv')
        keep=set(x.loc[x.in_published_187_sample_table==True,'sample_id'].astype(str))
        des['real']=des['real'][des['real']['sample'].isin(keep)].copy()
        assert len(des['real'])==177  # 187 roster entries include 10 recovery visits
    genes=source['real']['genes']
    assert np.array_equal(genes,source['floor']['genes'])
    measured=(source['real']['S'].sum(0)>0)&(source['floor']['S'].sum(0)>0)
    mats={}
    for ct in TYPES:
        mats[ct]={s:old.cpm_matrix(source[s],des[s],ct) for s in source}
    return genes,measured,mats,des

def run(coh,roster):
    genes,measured,mats,des=prepare(coh,roster)
    report=[]
    for ct in TYPES:
        gR=mats[ct]['real']; gF=mats[ct]['floor']
        cols=list(genes)
        yR=gR[cols].to_numpy(float); yF=gF[cols].to_numpy(float)
        cR=(gR._c=='case').to_numpy(int); cF=(gF._c=='case').to_numpy(int)
        pR=gR._u.astype(str).to_numpy(); pF=gF._u.astype(str).to_numpy()
        test=patient_cluster_deg if coh=='covid' else welch_deg
        def call(y,c,p):
            x=test(y,c,p) if coh=='covid' else test(y,c)
            return None if x is None else (x[0][measured],x[1][measured],x[2][measured])
        fullR=call(yR,cR,pR); fullF=call(yF,cF,pF)
        if fullR is None or fullF is None: continue
        point=metrics(fullR,fullF)
        iR=indices_by_patient(gR); iF=indices_by_patient(gF)
        if coh=='sepsis':
            sr=[np.array([u for u,ix in iR.items() if (gR.iloc[ix]._c==c).all()]) for c in ('case','ctrl')]
            sf=[np.array([u for u,ix in iF.items() if (gF.iloc[ix]._c==c).all()]) for c in ('case','ctrl')]
        else:
            sr=sf=None
        boot=[]
        for _ in range(B):
            rr,pr=resample(iR,sr); ff,pf=resample(iF,sf)
            ar=call(yR[rr],cR[rr],pr); bf=call(yF[ff],cF[ff],pf)
            if ar is not None and bf is not None:
                boot.append(metrics(ar,bf))
        boot=np.asarray(boot)
        row=dict(cohort=coh,roster=roster,celltype=ct,
                 n_real_patients=len(iR),n_floor_patients=len(iF),
                 n_real_case_units=int(cR.sum()),n_real_ctrl_units=int((1-cR).sum()),
                 n_floor_case_units=int(cF.sum()),n_floor_ctrl_units=int((1-cF).sum()),
                 n_genes=int(measured.sum()),n_real_deg=int(fullR[1].sum()),n_floor_deg=int(fullF[1].sum()),
                 n_boot=int(len(boot)),method='patient-cluster OLS CR1' if coh=='covid' else 'Welch')
        for j,k in enumerate(['logfc_r','recall','precision','sign_all','sign_real_deg']):
            row[k]=point[j]
            row[k+'_lo'],row[k+'_hi']=np.nanpercentile(boot[:,j],[2.5,97.5])
        report.append(row)
        print(coh,roster,ct,'r=',round(row['logfc_r'],3),'CI=',tuple(np.round([row['logfc_r_lo'],row['logfc_r_hi']],3)),flush=True)
    return report

def self_check():
    g=pd.DataFrame({'_u':['A','A','B'],'_c':['case','ctrl','case']})
    index=indices_by_patient(g)
    assert len(index['A'])==2
    draw=['A','A','B']
    rows=np.concatenate([index[u] for u in draw])
    assert rows.tolist()==[0,1,0,1,2]
    y=np.array([[1.,2.],[2.,4.],[3.,6.],[4.,8.],[5.,10.],[6.,12.]])
    cond=np.array([0,1,0,1,0,1]); cluster=np.array([0,0,1,1,2,2])
    # With complete pairs, the point estimate is the mean within-patient difference.
    result=patient_cluster_deg(y,cond,cluster)
    assert result is None  # fewer than four patients is rejected
    cluster=np.array([0,0,1,1,2,2,3,3]); y=np.vstack([y,[[7,14],[8,16]]]); cond=np.tile([0,1],4)
    result=patient_cluster_deg(y,cond,cluster)
    assert np.allclose(result[0],[1,2])

if __name__=='__main__':
    self_check()
    rows=[]
    for cohort,roster in [('sepsis','all'),('covid','published_187'),('covid','all_local_198')]:
        rows.extend(run(cohort,roster))
    pd.DataFrame(rows).to_csv(OUT/'corrected_summary.csv',index=False)
    print('saved',OUT/'corrected_summary.csv')
