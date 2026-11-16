"""Figure source data for the ESKD COVID clinical-association figures (2026-11-15).
Reads the 11/15 analysis outputs (RUN/out, RUN/cache) and writes every number a panel draws to data/*.csv."""
import os, numpy as np, pandas as pd, warnings, json
warnings.filterwarnings('ignore')
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_curve, roc_auc_score
RUN=os.environ.get('RUN','/home/claude/run/'); D='data/'
S=pd.read_csv(RUN+'out/sample_level_composition_metadata.csv'); P=pd.read_csv(RUN+'out/patient_level_composition_metadata.csv')
pos=P[P.group=='COVID_positive'].copy()
def intz(x): r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))
def bh(p):
    p=np.asarray(p,float); q=np.full(len(p),np.nan); k=~np.isnan(p)
    if k.sum(): q[k]=multipletests(p[k],method='fdr_bh')[1]
    return q
# ---- copies of analysis tables used directly
for f in ['assoc_patient_level','reliability_split_half_vs_cross_source','rare_subset_severity_DEG_summary','combined_DEG_permutation_control','rare_subclustering_summary','rare_subclustering_joint_cluster_severity','composition_vs_clinical_blood_counts']:
    pd.read_csv(RUN+f'out/{f}.csv').to_csv(D+f+'.csv',index=False)
pos.to_csv(D+'patient_level_covid_positive.csv',index=False)
# ---- 1. ROC: HSPC % vs death
fh=pos.groupby('patient_id').fatal_harmonized.max()
rows=[];aucs=[]
z=pos[['source','patient_id','rare_HSPC_pct']].copy(); z['z']=z.groupby('source').rare_HSPC_pct.transform(intz); comb=z.groupby('patient_id').z.mean()
for lab,x in [('Real',pos[pos.source=='Real'].set_index('patient_id').rare_HSPC_pct),('PBISC',pos[pos.source=='PBISC'].set_index('patient_id').rare_HSPC_pct),('Combined',comb)]:
    y=fh.reindex(x.index).astype(int); fpr,tpr,_=roc_curve(y,x)
    rows.append(pd.DataFrame(dict(analysis=lab,fpr=fpr,tpr=tpr))); aucs.append(dict(analysis=lab,auc=roc_auc_score(y,x),n=len(y),n_fatal=int(y.sum())))
pd.concat(rows).to_csv(D+'roc_hspc_fatal.csv',index=False); pd.DataFrame(aucs).to_csv(D+'roc_hspc_fatal_auc.csv',index=False); print(pd.DataFrame(aucs))
# ---- 2. clinical blood counts (Real metadata) vs composition
lab=S[(S.source=='Real')&(S.condition=='positive')].groupby('patient_id')[['temp_neut','temp_lymph']].mean(); lab['NLR']=lab.temp_neut/lab.temp_lymph
nl=pos.set_index('patient_id')[['source','T_pct','rare_HSPC_pct','NK_pct']].join(lab,how='inner').dropna(subset=['NLR']).reset_index()
nl.to_csv(D+'composition_vs_NLR.csv',index=False); print(nl.groupby('source').size())
# ---- 3. pseudobulk units (same code as deg.py) for gene panels + permutation nulls
ok=S[(S.condition=='positive')&S.in_primary_scope].set_index(['source','sample_id'])
def load(src,tag):
    z=np.load(RUN+f'cache/{tag}_pb.npz',allow_pickle=True); g=pd.read_csv(RUN+f'cache/{tag}_pb_groups.csv'); pb=pd.DataFrame(z['pb'],columns=z['genes'])
    pb=pb.T.groupby(level=0).sum().T; g['patient']=g.sample_id.map(ok.loc[src].patient_id); return pb,g
pbR,gR=load('Real','real'); pbP,gP=load('PBISC','pbisc'); genes=sorted(set(pbR.columns)&set(pbP.columns)); pbR=pbR[genes]; pbP=pbP[genes]
sev=S[S.condition=='positive'].drop_duplicates('patient_id').set_index('patient_id').severity_ord
def units(pb,g,fines):
    m=g.patient.notna() if fines is None else (g.patient.notna()&g.fine.isin(fines)); idx=np.where(m)[0]; pat=g.patient.values[idx]
    cnt=pd.DataFrame(pb.values[idx]).groupby(pat).sum(); cnt.columns=pb.columns; n=g.n.iloc[idx].groupby(pat).sum(); keep=n.index[n>=10]; cnt=cnt.loc[keep]
    return np.log2(cnt.div(cnt.sum(1),axis=0)*1e6+1)
def trend(L,s_):
    x=(s_-s_.mean())/s_.std(ddof=0); Lc=(L-L.mean())/L.std(ddof=0).replace(0,np.nan); r=(Lc.values*x.values[:,None]).mean(0); n=len(x)
    t=r*np.sqrt((n-2)/(1-r**2)); return r,2*stats.t.sf(np.abs(t),n-2)
ASC=['B_ASC_IgG','B_ASC_IgA','B_ASC_IgM','B_ASC_dividing','B_ASC_other']
SUB={'CD14mono_IFN':['CD14mono_IFN'],'HSPC':['HSPC'],'B_ASC_all':ASC,'NK CD56 bright':['NK CD56 bright'],'pDC':['pDC'],'REF_CD14mono':['CD14mono'],'REF_NK':['NK']}
tabs={}
for name in SUB: tabs[name]=pd.read_csv(RUN+f"out/deg/deg_{name.replace(' ','_')}_Combined.csv.gz",index_col=0)
rng=np.random.default_rng(1); prow=[]
Lkeep={}
for name,f in SUB.items():
    LR=units(pbR,gR,f); LP=units(pbP,gP,f); gl=tabs[name].index; Lkeep[name]=(LR,LP)
    zR=(LR[gl]-LR[gl].mean())/LR[gl].std(); zP=(LP[gl]-LP[gl].mean())/LP[gl].std()
    obs=int((tabs[name].trend_q<0.05).sum())
    for i in range(200):
        zs=zP.copy(); zs.index=rng.permutation(zs.index.values); Z=pd.concat([zR,zs]).groupby(level=0).mean(); tr,tp=trend(Z,sev.reindex(Z.index))
        prow.append(dict(subset=name,perm=i,n_DEG=int((bh(tp)<0.05).sum()),observed=obs))
PN=pd.DataFrame(prow); PN.to_csv(D+'combined_DEG_permutation_null.csv',index=False)
print(PN.groupby('subset').agg(obs=('observed','first'),med=('n_DEG','median'),q95=('n_DEG',lambda x:np.quantile(x,.95)),frac=('n_DEG',lambda x:None)))
# ---- 4. CD14mono_IFN gene-level trend: Real vs PBISC vs Combined
g3=[pd.read_csv(RUN+f'out/deg/deg_CD14mono_IFN_{s}.csv.gz',index_col=0)[['trend_r','trend_q']].add_prefix(s+'_') for s in ['Real','PBISC','Combined']]
G=pd.concat(g3,axis=1).dropna(subset=['Combined_trend_r']); G.index.name='gene'; G.to_csv(D+'cd14monoIFN_gene_trend.csv'); print('combined sig',(G.Combined_trend_q<.05).sum())
LR,LP=Lkeep['CD14mono_IFN']; show=['S100A8','S100A12','CD163','FKBP5','HLA-DRA','HLA-DPB1']
rows=[]
for src,L in [('Real',LR),('PBISC',LP)]:
    zz=(L[show]-L[show].mean())/L[show].std()
    for g_ in show:
        for p_,v in zz[g_].items(): rows.append(dict(source=src,patient_id=p_,gene=g_,z=v,severity_ord=sev.get(p_)))
pd.DataFrame(rows).to_csv(D+'cd14monoIFN_gene_patient_z.csv',index=False)
# ---- 5. DEG count summary
R=pd.read_csv(RUN+'out/rare_subset_severity_DEG_summary.csv'); R.to_csv(D+'rare_subset_severity_DEG_summary.csv',index=False)
