"""Ver3 (2026-11-15): reference-based deconvolution of the same bulk RNA-seq that PBISC was conditioned on - baseline for
'does PBISC add anything over deconvolution?'.
Signature: Real scRNA (all deposited samples), mean UMI per cell for 30 groups over the 5,762 shared genes; 30 marker genes
per group (highest log fold change versus the next-highest group, mean >= 0.05 UMI/cell). NNLS on CPM; fractions = weights
normalised to 1. Two variants: uncalibrated, and calibrated with gene-wise bulk/scRNA ratios estimated on the 111 bulk samples
drawn on the same day as a Real scRNA sample (cross-fitted by patient halves for the accuracy check).
Accuracy: estimated vs scRNA composition on those 111 matched samples. Association: patient means over the PBISC input
samples (COVID+ primary scope) vs peak severity."""
import os, numpy as np, pandas as pd, anndata as ad, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from scipy.optimize import nnls
RUN=os.environ.get('RUN','/home/claude/run/'); UP='/mnt/user-data/uploads/PBISC/'; D='data_v3/'
S=pd.read_csv(RUN+'out/sample_level_composition_metadata.csv'); sev=S[S.condition=='positive'].drop_duplicates('patient_id').set_index('patient_id').severity_ord
GROUP={'NK CD56 bright':'NK CD56bright','NK':'NK','NK IFN stim':'NK','B_ASC_dividing':'B_ASC_dividing','B_ASC_IgA':'B_ASC_other','B_ASC_IgG':'B_ASC_other','B_ASC_IgM':'B_ASC_other','B_ASC_other':'B_ASC_other',
 'HSPC':'HSPC','pDC':'pDC','DC1':'DC1','DC2':'cDC2','DC3':'cDC2','DC3_IFN':'cDC2','ASDC':'cDC2','MAIT':'MAIT','T g/d':'T g/d','T CD4 cycling':'T CD4 Helper','ILC':'ILC',
 'CD14mono':'CD14mono','CD14mono_IFN':'CD14mono_IFN','CD14mono_dex':'CD14mono_dex','Int.mono':'CD14mono','CD16mono':'CD16mono','CD16mono_C1':'CD16mono','CD16mono_IFN':'CD16mono',
 'T CD4 Naive':'T CD4 Naive','T CD4 Naive IFN stim':'T CD4 Naive','T CD4 Helper':'T CD4 Helper','T CD4 Helper IFN stim':'T CD4 Helper','T CD4 CTL':'T CD4 CTL','T reg':'T reg','T reg IFN stim':'T reg',
 'T CD8 Naive':'T CD8 Naive','T CD8 Naive IFN stim':'T CD8 Naive','T CD8 CM':'T CD8 CM','T CD8 CTL':'T CD8 CTL','T CD8 CTL IFN stim':'T CD8 CTL','NKT':'NKT',
 'B_naive':'B_naive','B_switched_mem':'B_memory','B_non-switched_mem':'B_memory','B_CD11c':'B_memory','Platelet':'Platelet','RBC':'RBC'}
BROAD={'NK CD56bright':'NK','NK':'NK','B_ASC_dividing':'B','B_ASC_other':'B','B_naive':'B','B_memory':'B','HSPC':'Other','pDC':'DC','DC1':'DC','cDC2':'DC','ILC':'Other','Platelet':'Other','RBC':'Other',
 'CD14mono':'Mono','CD14mono_IFN':'Mono','CD14mono_dex':'Mono','CD16mono':'Mono'}
z=np.load(RUN+'cache/real_pb.npz',allow_pickle=True); g=pd.read_csv(RUN+'cache/real_pb_groups.csv'); genes=pd.Index(z['genes'])
pb=pd.DataFrame(z['pb'],columns=genes); pb=pb.T.groupby(level=0).sum().T; genes=pb.columns
g['grp']=g.fine.map(GROUP); assert g.grp.notna().all()
sig=pb.groupby(g.grp.values).sum().div(g.groupby('grp').n.sum(),axis=0)          # mean UMI per cell
lg=np.log2(sig+0.01); mk=[]
for k in sig.index:
    other=lg.drop(index=k).max(0); fc=(lg.loc[k]-other)[sig.loc[k]>=0.05]; mk+=list(fc.sort_values(ascending=False).index[:30])
mk=sorted(set(mk)); print('groups',len(sig),'marker genes',len(mk))
# bulk
B=pd.read_csv(UP+'bulkonly_lora5_20261006/analysis/bulk_public/zenodo6497251_imperial_eskd/htseq_counts.csv',index_col=0); B=B[~B.index.str.startswith('__')]; B.index=B.index.str.split('.').str[0]; B=B.groupby(level=0).sum()
var=ad.read_h5ad(UP+'deliverables/annotated_h5ad_20261112/pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad',backed='r').var
e2s=pd.Series(var.gene_symbol.values,index=var.index); e2s=e2s[e2s.index.isin(B.index)]; e2s=e2s[~e2s.duplicated()]
cpm=(B.div(B.sum(0),axis=1)*1e6).loc[e2s.index]; cpm.index=e2s.values; cpm=cpm.groupby(level=0).sum()
mk=[m for m in mk if m in cpm.index]; print('markers in bulk',len(mk))
def deconv(M,sigm):
    A=sigm.loc[:,mk].T.values; out={}
    for s in M.columns: w,_=nnls(A,M.loc[mk,s].values); out[s]=w/w.sum() if w.sum()>0 else w
    return pd.DataFrame(out,index=sigm.index).T
def features(F):
    o=pd.DataFrame(index=F.index)
    for f,cols in [('rare_NK CD56 bright_pct',['NK CD56bright']),('rare_B_ASC_dividing_pct',['B_ASC_dividing']),('rare_HSPC_pct',['HSPC']),('rare_DC1_pct',['DC1']),('rare_pDC_pct',['pDC']),('rare_CD14mono_IFN_pct',['CD14mono_IFN']),('rare_MAIT_pct',['MAIT'])]:
        o[f]=100*F[cols].sum(1)
    o['DC_within_Other_pct']=100*F[[c for c in F if BROAD.get(c)=='DC']].sum(1); o['NK_pct']=100*F[['NK','NK CD56bright']].sum(1)
    o['T_pct']=100*F[[c for c in F if c.startswith('T ') or c in ('MAIT','NKT')]].sum(1); o['Mono_pct']=100*F[[c for c in F if BROAD.get(c)=='Mono']].sum(1)
    return o
# matched samples (Real scRNA vs bulk same day)
rm=S[S.source=='Real'].set_index('sample_id')
def bk(r):
    if r.condition=='negative': return f'{r.patient_id}_negative'
    return f'{r.patient_id}_positive_{int(r.days_from_first_positive_swab)}' if pd.notna(r.days_from_first_positive_swab) else None
rm['bulk_id']=[bk(r) if r.condition in('positive','negative') else None for _,r in rm.iterrows()]
mt=rm[rm.bulk_id.isin(cpm.columns)&rm.qc_pass_min200cells.astype(bool)]
truth=rm.loc[mt.index,['rare_NK CD56 bright_pct','rare_B_ASC_dividing_pct','rare_HSPC_pct','rare_DC1_pct','rare_pDC_pct','rare_CD14mono_IFN_pct','rare_MAIT_pct','DC_within_Other_pct','NK_pct','T_pct','Mono_pct']]
# real pseudobulk CPM of matched samples for calibration
idx=g.sample_id.isin(mt.index).values; rpb=pb[idx].groupby(g.sample_id.values[idx]).sum(); rcpm=rpb.div(rpb.sum(1),axis=0)*1e6
pats=np.array(sorted(mt.patient_id.unique())); rng=np.random.default_rng(0); rng.shuffle(pats); half={p:i%2 for i,p in enumerate(pats)}
def calib(samples):
    b=cpm[mt.loc[samples,'bulk_id']].mean(1); r=rcpm.loc[samples].mean(0); return ((b+1)/(r.reindex(b.index)+1)).fillna(1)
acc=[]; EST={}
for variant in ['uncalibrated','calibrated']:
    parts=[]
    for h in [0,1]:
        test=[s for s in mt.index if half[mt.loc[s,'patient_id']]==h]; train=[s for s in mt.index if half[mt.loc[s,'patient_id']]!=h]
        sg=sig.copy() if variant=='uncalibrated' else sig.mul(calib(train).reindex(sig.columns).fillna(1),axis=1)
        e=deconv(cpm[mt.loc[test,'bulk_id']],sg); e.index=test; parts.append(features(e))
    E=pd.concat(parts); EST[variant]=E
    for f in truth.columns:
        acc.append(dict(variant=variant,feature=f,n=len(E),pearson_r=stats.pearsonr(E[f],truth.loc[E.index,f])[0],spearman_r=stats.spearmanr(E[f],truth.loc[E.index,f])[0],mean_est=E[f].mean(),mean_true=truth[f].mean()))
ACC=pd.DataFrame(acc); ACC.to_csv(D+'deconvolution_accuracy_matched_v3.csv',index=False); print(ACC.round(3).to_string())
# association on PBISC input samples
pbm=S[S.source=='PBISC'].set_index('sample_id'); inp=pbm[(pbm.condition=='positive')&pbm.in_primary_scope.astype(bool)].index
res=[]; PATS={}
for variant in ['uncalibrated','calibrated']:
    sg=sig.copy() if variant=='uncalibrated' else sig.mul(calib(list(mt.index)).reindex(sig.columns).fillna(1),axis=1)
    F=features(deconv(cpm[inp],sg)); P=F.groupby(pbm.loc[inp,'patient_id'].values).mean(); PATS[variant]=P
    P.assign(variant=variant).to_csv(D+f'deconvolution_patient_{variant}_v3.csv')
    for f in F.columns:
        q=stats.spearmanr(P[f],sev[P.index]); po=P[~P.index.isin(rm.patient_id.unique())]; q2=stats.spearmanr(po[f],sev[po.index])
        res.append(dict(method=f'deconvolution ({variant})',feature=f,n=len(P),rho=q[0],p=q[1],rho_new_patients=q2[0],p_new_patients=q2[1],n_new=len(po)))
A=pd.read_csv(D+'composition_association_v3.csv').set_index('feature')
for f in truth.columns:
    res.append(dict(method='PBISC',feature=f,n=int(A.loc[f,'PBISC_n']),rho=A.loc[f,'PBISC_rho'],p=A.loc[f,'PBISC_p'],rho_new_patients=A.loc[f,'PBISC_only_rho'],p_new_patients=A.loc[f,'PBISC_only_p'],n_new=int(A.loc[f,'PBISC_only_n'])))
    res.append(dict(method='Real scRNA',feature=f,n=int(A.loc[f,'Real_n']),rho=A.loc[f,'Real_rho'],p=A.loc[f,'Real_p']))
RES=pd.DataFrame(res); RES.to_csv(D+'deconvolution_vs_pbisc_severity_v3.csv',index=False)
pd.set_option('display.width',250); print(RES.round(3).sort_values(['feature','method']).to_string())
