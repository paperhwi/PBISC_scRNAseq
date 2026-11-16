"""Ver3 (2026-11-15) composition analyses.
1. Patient tables (COVID+ primary scope, Real samples >=200 cells) for 20 features.
2. Association with peak severity: Real (39), PBISC (57), PBISC-only (29, patients absent from Real),
   Extension = Real (39) + PBISC-only (29), each rank-INT within its own source set (no patient counted twice),
   and the ver1/ver2 'averaged' combination (shared patients averaged) for reference.
3. Cohort-extension power: n Real patients drawn without replacement (1,000 draws) tested alone, with the 29 PBISC-only
   patients, and with the 29 PBISC-only patients after shuffling their severity labels.
4. Prediction-powered inference (PPI++, ppi_py): slope of log10(% + 0.05) per severity step; labelled = 28 shared
   patients (Real Y, PBISC prediction), unlabelled = 29 PBISC-only patients.
5. Robustness: published 187-sample roster; age/sex/centre and days-from-symptoms adjustment; first sample only;
   other adapters (final-architecture draw 0, 60k labelled subset; legacy K4000 draws 0-4, 100 cells/sample/draw).
6. Label recoverability (held-out real cells through the PBISC label-transfer route)."""
import os, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import ppi_py
RUN=os.environ.get('RUN','/home/claude/run/'); UP='/mnt/user-data/uploads/PBISC/'; D='data_v3/'; os.makedirs(D,exist_ok=True)
S=pd.read_csv(RUN+'out/sample_level_composition_metadata.csv')
RARE13=['ASDC','B_ASC_IgG','B_ASC_dividing','CD14mono_IFN','CD14mono_dex','DC1','HSPC','ILC','MAIT','NK CD56 bright','T CD4 cycling','T g/d','pDC']
FEATS=['T_pct','B_pct','Mono_pct','NK_pct','Other_pct','DC_within_Other_pct','rare13_total_pct']+[f'rare_{r}_pct' for r in RARE13]
NICE={f:f.replace('rare_','').replace('_pct','').replace('_within_Other','') for f in FEATS}; NICE.update({'T_pct':'T cell','B_pct':'B cell','Mono_pct':'Monocyte','NK_pct':'NK cell','Other_pct':'Other','rare13_total_pct':'Rare (13 total)'})
def rint(x): r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(r))
sev=S[S.condition=='positive'].drop_duplicates('patient_id').set_index('patient_id').severity_ord
ROSTER_DROP_P=['C116','C139']; ROSTER_DROP_S=['C141_0903']

def patient_table(src,roster187=False,first_only=False):
    d=S[(S.source==src)&(S.condition=='positive')&S.in_primary_scope.astype(bool)&S.qc_pass_min200cells.astype(bool)]
    if roster187: d=d[~d.patient_id.isin(ROSTER_DROP_P)&~d.sample_id.isin(ROSTER_DROP_S)]
    if first_only: d=d.sort_values('days_from_first_positive_swab').groupby('patient_id').head(1)
    g=d.groupby('patient_id'); P=g[FEATS].mean()
    P['sev']=g.severity_ord.first(); P['age']=g.age.first(); P['sex']=g.sex.first(); P['centre']=g.centre_or_wave.first()
    P['days']=g.days_from_first_symptoms.mean(); P['fatal']=g.fatal_harmonized.first(); return P
PR,PP=patient_table('Real'),patient_table('PBISC'); PO=PP[~PP.index.isin(PR.index)]; SH=sorted(set(PR.index)&set(PP.index))
print('Real',len(PR),'PBISC',len(PP),'PBISC-only',len(PO),'shared',len(SH))

def ext(a,b,col):  # extension: rank-INT within each source set, no patient twice
    z=pd.concat([pd.Series(rint(a[col].values),index=a.index),pd.Series(rint(b[col].values),index=b.index)]); return z
def avg(a,b,col):
    return pd.concat([pd.Series(rint(a[col].values),index=a.index),pd.Series(rint(b[col].values),index=b.index)]).groupby(level=0).mean()
def sp(x,idx): return stats.spearmanr(x,sev[idx])

rows=[]
for f in FEATS:
    r=dict(feature=f,name=NICE[f])
    for lab,P in [('Real',PR),('PBISC',PP),('PBISC_only',PO)]:
        q=sp(P[f],P.index); r[f'{lab}_rho'],r[f'{lab}_p'],r[f'{lab}_n']=q[0],q[1],len(P)
    z=ext(PR,PO,f); q=sp(z,z.index); r['Extension_rho'],r['Extension_p'],r['Extension_n']=q[0],q[1],len(z)
    z=avg(PR,PP,f); q=sp(z,z.index); r['Averaged_rho'],r['Averaged_p'],r['Averaged_n']=q[0],q[1],len(z)
    sh=[s for s in SH]; r['same_patient_rho']=stats.spearmanr(PR.loc[sh,f],PP.loc[sh,f])[0]
    rows.append(r)
A=pd.DataFrame(rows)
for lab in ['Real','PBISC','PBISC_only','Extension','Averaged']:
    ok=A[f'{lab}_p'].notna(); A[f'{lab}_q']=np.nan; A.loc[ok,f'{lab}_q']=multipletests(A.loc[ok,f'{lab}_p'],method='fdr_bh')[1]
REC=pd.read_csv(UP+'NewLora/analysis_20261102/verification_20261102/holdout_rare_label_precision.csv'); REC=REC[REC.cohort=='covid'].set_index('subset')
A['label_recall']=[REC.recall.get(f.replace('rare_','').replace('_pct',''),np.nan) for f in FEATS]
A['label_precision']=[REC.precision.get(f.replace('rare_','').replace('_pct',''),np.nan) for f in FEATS]
A['recoverable']=np.where(A.label_recall.isna(),True,(A.label_recall>=0.5)&(A.label_precision>=0.5))
A.to_csv(D+'composition_association_v3.csv',index=False)
pd.set_option('display.width',250)
print(A[['name','Real_rho','Real_q','PBISC_rho','PBISC_q','PBISC_only_rho','PBISC_only_p','Extension_rho','Extension_q','Averaged_q','same_patient_rho','recoverable']].round(3).to_string())

# ---- cohort-extension power
rng=np.random.default_rng(20261115); NS=[10,15,20,25,30,35]; REPS=1000
pw=[]
for f in FEATS:
    full=sp(PR[f],PR.index); sg=np.sign(full[0]); real_sig=full[1]<0.05
    for n in NS:
        h=np.zeros(3)
        for _ in range(REPS):
            s=rng.choice(PR.index,n,replace=False); a=PR.loc[s]
            q0=sp(a[f],s)
            z=ext(a,PO,f); q1=sp(z,z.index)
            z=ext(a,PO,f); q2=stats.spearmanr(z.values,np.concatenate([sev[s].values,rng.permutation(sev[PO.index].values)]))
            # q2: PBISC-only values kept, their severity labels permuted (PBISC carries no severity information)
            for k,q in enumerate([q0,q1,q2]):
                h[k]+=(q[1]<0.05)&((np.sign(q[0])==sg) if real_sig else True)
        pw.append(dict(feature=f,name=NICE[f],real_full_rho=full[0],real_full_p=full[1],criterion='p<0.05, sign of full Real' if real_sig else 'p<0.05, any sign (Real null)',
                       n_real=n,real_alone=h[0]/REPS,plus_pbisc_only=h[1]/REPS,plus_shuffled_pbisc=h[2]/REPS))
PW=pd.DataFrame(pw); PW.to_csv(D+'extension_power_v3.csv',index=False)
print(PW[PW.feature.isin(['rare_NK CD56 bright_pct','rare_B_ASC_dividing_pct','DC_within_Other_pct','rare_HSPC_pct','NK_pct','rare_DC1_pct'])].round(2).to_string())

# ---- PPI++
pp=[]
for f in FEATS:
    t=lambda x: np.log10(x+0.05)
    Xl=np.c_[np.ones(len(SH)),sev[SH].values]; Yl=t(PR.loc[SH,f].values); Fl=t(PP.loc[SH,f].values)
    Xu=np.c_[np.ones(len(PO)),sev[PO.index].values]; Fu=t(PO[f].values)
    Xr=np.c_[np.ones(len(PR)),sev[PR.index].values]; Yr=t(PR[f].values)
    for lab,(lo,hi),est in [('classical_shared28',ppi_py.classical_ols_ci(Xl,Yl),None),('classical_real39',ppi_py.classical_ols_ci(Xr,Yr),None),
                            ('PPI++_28+29',ppi_py.ppi_ols_ci(Xl,Yl,Fl,Xu,Fu),ppi_py.ppi_ols_pointestimate(Xl,Yl,Fl,Xu,Fu))]:
        lo,hi=np.asarray(lo)[1],np.asarray(hi)[1]; b=(lo+hi)/2 if est is None else np.asarray(est)[1]; se=(hi-lo)/(2*1.959964)
        pp.append(dict(feature=f,name=NICE[f],estimator=lab,slope=b,ci_lo=lo,ci_hi=hi,se=se,p=2*stats.norm.sf(abs(b/se)),n_labelled=len(SH) if '28' in lab else len(PR),n_unlabelled=len(PO) if 'PPI' in lab else 0))
    # naive: PBISC predictions treated as truth for PBISC-only (biased reference)
    Xn=np.r_[Xl,Xu]; Yn=np.r_[Yl,Fu]; lo,hi=ppi_py.classical_ols_ci(Xn,Yn); lo,hi=lo[1],hi[1]; se=(hi-lo)/3.919928
    pp.append(dict(feature=f,name=NICE[f],estimator='naive_imputed',slope=(lo+hi)/2,ci_lo=lo,ci_hi=hi,se=se,p=2*stats.norm.sf(abs((lo+hi)/2/se)),n_labelled=len(SH),n_unlabelled=len(PO)))
PPI=pd.DataFrame(pp); PPI.to_csv(D+'ppi_v3.csv',index=False)
print(PPI[PPI.feature.isin(['rare_NK CD56 bright_pct','rare_B_ASC_dividing_pct','DC_within_Other_pct','rare_HSPC_pct'])].round(4).to_string())

# ---- robustness
rob=[]
KEY=['rare_NK CD56 bright_pct','rare_B_ASC_dividing_pct','DC_within_Other_pct','rare_DC1_pct','rare_HSPC_pct','NK_pct','T_pct']
for lab,src,kw in [('all deposited samples','Real',{}),('published 187-sample roster','Real',{'roster187':True}),('first sample per patient','Real',{'first_only':True}),
                   ('all samples','PBISC',{}),('first sample per patient','PBISC',{'first_only':True})]:
    P=patient_table(src,**kw)
    for f in KEY:
        q=sp(P[f],P.index); d=P.copy(); d['y']=rint(d[f]); dd=d.dropna(subset=['age'])
        m1=smf.ols('y~sev+age+C(sex)+C(centre)' if d.centre.nunique()>1 else 'y~sev+age+C(sex)',dd).fit()
        d2=d.dropna(subset=['days']); m2=smf.ols('y~sev+days',d2).fit(); m0=smf.ols('y~sev',d2).fit()
        rob.append(dict(source=src,analysis=lab,feature=f,name=NICE[f],n=len(P),rho=q[0],p=q[1],p_adj_age_sex_centre=m1.pvalues['sev'],n_days=len(d2),p_same_subset=m0.pvalues['sev'],p_adj_days=m2.pvalues['sev']))
# other adapters
lm=pd.read_csv(UP+'deliverables/annotated_h5ad_20261112/tables/label_mapping_covid.csv').set_index('fine')
def comp_from_labels(df,sid,lab):
    n=df.groupby(sid).size(); f=pd.crosstab(df[sid],df[lab]); b6=df[lab].map(lm.six_type_broad)
    out=pd.DataFrame(index=n.index)
    for col,bb in [('NK_pct','NK cell'),('T_pct','T cell')]: out[col]=100*(b6==bb).groupby(df[sid]).sum()/n
    out['DC_within_Other_pct']=100*(b6=='DC').groupby(df[sid]).sum()/n
    for r in ['NK CD56 bright','B_ASC_dividing','DC1','HSPC']: out[f'rare_{r}_pct']=100*f.get(r,pd.Series(0,index=f.index)).reindex(n.index).fillna(0)/n
    return out
md=S[S.source=='PBISC'].set_index('sample_id')
def adapter_rows(name,C):
    C=C[C.index.isin(md.index[(md.condition=='positive')&md.in_primary_scope.astype(bool)])]
    pat=md.patient_id.reindex(C.index); P=C.groupby(pat.values).mean()
    for f in KEY:
        if f in P: q=sp(P[f],P.index); rob.append(dict(source='PBISC',analysis=name,feature=f,name=NICE[f],n=len(P),rho=q[0],p=q[1]))
L0=pd.read_csv(UP+'NewLora/analysis_20261102/covid/ceiling_cell_labels.csv.gz'); print(L0.columns.tolist())
labcol=[c for c in L0.columns if 'fine' in c.lower()][0]; adapter_rows('final model, adapter draw 0 (60k labelled subset)',comp_from_labels(L0,'sample',labcol))
o=pd.read_parquet(RUN+'cache/pbisc_obs.parquet')
o60=o[o.in_published_60k.astype(bool)] if 'in_published_60k' in o else None
if o60 is not None: adapter_rows('final model, adapter draw 1 (same 60k subset)',comp_from_labels(o60,'sample_id','fine'))
K0=pd.read_parquet(RUN+'cache/draw0_knn_labels.parquet'); adapter_rows('final model, adapter draw 0, all cells (kNN from draw-1 labels)',comp_from_labels(K0,'sample_id','fine'))
T=pd.read_csv(UP+'bulkonly_lora5_20261006/analysis/integration_report_20261012/results/covid_k4000_transfer.csv.gz')
for dname,g in T.groupby('draw'): adapter_rows(f'legacy K4000 model, {dname} (100 cells/sample)',comp_from_labels(g,'sample_id','label'))
adapter_rows('legacy K4000 model, draws 0-4 pooled (500 cells/sample)',comp_from_labels(T,'sample_id','label'))
ROB=pd.DataFrame(rob); ROB.to_csv(D+'robustness_v3.csv',index=False)
print(ROB.round(4).to_string())
