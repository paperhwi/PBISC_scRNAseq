import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
S=pd.read_csv('out/sample_level_composition_metadata.csv'); P=pd.read_csv('out/patient_level_composition_metadata.csv')
FE=[c for c in P.columns if c.endswith('_pct')]
pos=P[P.group=='COVID_positive'].copy(); neg=P[P.group=='COVID_negative_control'].copy()
fh=pos.groupby('patient_id').fatal_harmonized.max()
def intz(x):
    r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))
def combined(df,feats):
    """within-source rank-INT, then average over the sources available for each patient (one row per patient)"""
    z=df[['source','patient_id']+feats].copy()
    for f in feats: z[f]=z.groupby('source')[f].transform(intz)
    return z.groupby('patient_id')[feats].mean()
meta=pos.drop_duplicates('patient_id').set_index('patient_id')[['severity_ord','age','sex']]; meta['fatal']=fh
def bh(p): 
    p=np.asarray(p,float); q=np.full(len(p),np.nan); ok=~np.isnan(p)
    if ok.sum(): q[ok]=multipletests(p[ok],method='fdr_bh')[1]
    return q
def tests(X,m,label):
    rows=[]
    for f in X.columns:
        x=X[f]; d=pd.concat([x.rename('x'),m],axis=1).dropna(subset=['x','severity_ord'])
        n=len(d); varied=d.x.nunique()>2
        rho,p1=stats.spearmanr(d.x,d.severity_ord) if varied else (np.nan,np.nan)
        hi=d.x[d.severity_ord>=2]; lo=d.x[d.severity_ord<2]
        u,p2=stats.mannwhitneyu(hi,lo) if varied else (np.nan,np.nan); auc=u/(len(hi)*len(lo)) if varied else np.nan
        df_=d.dropna(subset=['fatal']); a=df_.x[df_.fatal==1]; b=df_.x[df_.fatal==0]
        u3,p3=stats.mannwhitneyu(a,b) if varied and len(a)>=3 else (np.nan,np.nan)
        # age/sex adjusted trend (rank-INT outcome)
        try:
            d2=d.assign(z=intz(d.x)); fit=smf.ols('z~severity_ord+age+C(sex)',d2).fit(); padj=fit.pvalues['severity_ord']; badj=fit.params['severity_ord']
        except Exception: padj=badj=np.nan
        rows.append(dict(analysis=label,feature=f,n_patients=n,n_severe_critical=len(hi),n_mild_moderate=len(lo),
            median_mild_moderate=lo.median(),median_severe_critical=hi.median(),
            spearman_rho_severity=rho,p_severity_trend=p1,AUC_sevcrit_vs_mildmod=auc,p_sevcrit_vs_mildmod=p2,
            beta_trend_age_sex_adj=badj,p_trend_age_sex_adj=padj,
            n_fatal=len(a),n_survived=len(b),median_fatal=a.median(),median_survived=b.median(),AUC_fatal=(u3/(len(a)*len(b)) if not np.isnan(u3) else np.nan),p_fatal=p3))
    R=pd.DataFrame(rows)
    for c in ['p_severity_trend','p_sevcrit_vs_mildmod','p_trend_age_sex_adj','p_fatal']: R[c.replace('p_','q_')]=bh(R[c])
    return R
out=[]
for src in ['Real','PBISC']:
    d=pos[pos.source==src].set_index('patient_id'); out.append(tests(d[FE],meta,src))
out.append(tests(combined(pos,FE),meta,'Combined_union_patients'))
both=pos.groupby('patient_id').source.nunique().eq(2); ids_both=both.index[both]
d=pos[(pos.source=='PBISC')&~pos.patient_id.isin(ids_both)].set_index('patient_id'); out.append(tests(d[FE],meta,'PBISC_only_patients_not_in_Real(independent)'))
# naive stacking that ignores patient overlap (shown only as the wrong way)
R=pd.concat(out); R.to_csv('out/assoc_patient_level.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
show=R[['analysis','feature','n_patients','spearman_rho_severity','p_severity_trend','q_severity_trend','AUC_sevcrit_vs_mildmod','q_sevcrit_vs_mildmod','p_trend_age_sex_adj','n_fatal','AUC_fatal','p_fatal','q_fatal']]
print(show.round(3).to_string())
# Real vs PBISC agreement within the same patients
rows=[]
w=pos[pos.patient_id.isin(ids_both)].pivot(index='patient_id',columns='source',values=FE)
for f in FE:
    a=w[f]['Real'];b=w[f]['PBISC']
    r,p=stats.spearmanr(a,b) if b.nunique()>2 else (np.nan,np.nan)
    rows.append(dict(feature=f,n_patients=len(a),real_mean=a.mean(),pbisc_mean=b.mean(),spearman_real_vs_pbisc=r,p=p))
A=pd.DataFrame(rows); A.to_csv('out/real_vs_pbisc_same_patient_agreement.csv',index=False); print(A.round(3).to_string())
# effect-size concordance
piv=R.pivot(index='feature',columns='analysis',values='spearman_rho_severity')
print('rho(severity effect Real vs PBISC across features)',stats.spearmanr(piv['Real'],piv['PBISC'],nan_policy='omit'))
print('Real vs PBISC-only',stats.spearmanr(piv['Real'],piv['PBISC_only_patients_not_in_Real(independent)'],nan_policy='omit'))
# case vs control
rows=[]
for src in ['Real','PBISC']:
    a=pos[pos.source==src]; b=neg[neg.source==src]
    for f in FE:
        if a[f].nunique()<3: continue
        u,p=stats.mannwhitneyu(a[f],b[f]); rows.append(dict(source=src,feature=f,n_pos=len(a),n_neg=len(b),median_pos=a[f].median(),median_neg=b[f].median(),AUC=u/len(a)/len(b),p=p))
CC=pd.DataFrame(rows); CC["q"]=CC.groupby("source").p.transform(bh); CC.to_csv('out/assoc_case_vs_control.csv',index=False)
# sample-level: contemporaneous severity (WHO_temp_severity); OLS on standardized arcsine-sqrt proportion, patient-clustered SE
Sq=S[(S.condition=='positive')&S.in_primary_scope&S.qc_pass_min200cells].copy()
Sq['logCRP']=np.log10(Sq.temp_CRP+1)
rows=[]
for lab,d,form,term in [('Real',Sq[Sq.source=='Real'],'y~temp_severity_ord+days_from_first_positive_swab','temp_severity_ord'),('PBISC',Sq[Sq.source=='PBISC'],'y~temp_severity_ord+days_from_first_positive_swab','temp_severity_ord'),
                   ('Combined',Sq,'y~temp_severity_ord+days_from_first_positive_swab+C(source)','temp_severity_ord'),('Real_CRP',Sq[Sq.source=='Real'],'y~logCRP+days_from_first_positive_swab','logCRP')]:
    for f in FE:
        dd=d.dropna(subset=[term,'days_from_first_positive_swab',f]).copy()
        if dd[f].nunique()<5: continue
        dd['y']=dd.groupby('source')[f].transform(lambda x:(np.arcsin(np.sqrt(x/100))-np.arcsin(np.sqrt(x/100)).mean())/np.arcsin(np.sqrt(x/100)).std())
        dd=dd.dropna(subset=['y']).reset_index(drop=True)
        try: m=smf.ols(form,dd).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(dd.patient_id)[0]})
        except Exception as e: print('fail',lab,f,type(e).__name__,str(e)[:80]); continue
        rows.append(dict(analysis=lab,feature=f,n_samples=len(dd),n_patients=dd.patient_id.nunique(),beta_std=m.params[term],se=m.bse[term],z=m.tvalues[term],p=m.pvalues[term]))
L=pd.DataFrame(rows); L['q']=L.groupby('analysis').p.transform(bh); L.to_csv('out/assoc_sample_level_clustered.csv',index=False)
print(L.groupby('analysis')[['n_samples','n_patients']].first())
print(L.pivot(index='feature',columns='analysis',values=['beta_std','q']).round(3).to_string())
# noise ceiling: within-source split-half reliability of patient-level proportions (patients with >=2 samples, 200 random splits)
rng=np.random.default_rng(0); rows=[]
for src in ['Real','PBISC']:
    d=Sq[Sq.source==src]; g=[x for _,x in d.groupby('patient_id') if len(x)>=2]
    for f in FE:
        rs=[]
        for _ in range(200):
            a=[];b=[]
            for x in g:
                idx=rng.permutation(len(x)); h=len(x)//2; a.append(x[f].values[idx[:h]].mean()); b.append(x[f].values[idx[h:]].mean())
            if np.std(a)>0 and np.std(b)>0: rs.append(stats.spearmanr(a,b)[0])
        rows.append(dict(source=src,feature=f,n_patients=len(g),split_half_spearman=np.mean(rs) if rs else np.nan))
N=pd.DataFrame(rows).pivot(index='feature',columns='source',values='split_half_spearman'); N['real_vs_pbisc_same_patient']=A.set_index('feature').spearman_real_vs_pbisc
N.to_csv('out/reliability_split_half_vs_cross_source.csv'); print(N.round(2).to_string())
# is PBISC composition tracking whole-blood counts? clinical lymphocyte / neutrophil counts exist only in Real metadata -> patient level, 28 shared patients
lab=S[(S.source=='Real')&(S.condition=='positive')].groupby('patient_id')[['temp_lymph','temp_neut','temp_WCC','temp_mono','temp_CRP']].mean(); lab['lymph_frac_of_WCC']=lab.temp_lymph/lab.temp_WCC; lab['NLR']=lab.temp_neut/lab.temp_lymph
rows=[]
for src in ['Real','PBISC']:
    d=pos[pos.source==src].set_index('patient_id').join(lab,how='inner')
    for f in ['T_pct','NK_pct','B_pct','Mono_pct','Other_pct','rare_HSPC_pct']:
        for l in ['temp_lymph','lymph_frac_of_WCC','NLR','temp_neut','temp_CRP']:
            x=d[[f,l]].dropna(); r,p=stats.spearmanr(x[f],x[l]); rows.append(dict(source=src,feature=f,lab=l,n=len(x),rho=r,p=p))
W=pd.DataFrame(rows); W.to_csv('out/composition_vs_clinical_blood_counts.csv',index=False)
print(W.pivot(index=['feature','source'],columns='lab',values='rho').round(2).to_string()); print(W.groupby('source').n.max())
