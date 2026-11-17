"""Independent re-implementation (no code shared with run_ver8_pipeline.py) of core Fig. 7 ver8 statistics."""
import os, warnings, numpy as np, pandas as pd, statsmodels.api as sm
warnings.filterwarnings('ignore')
from scipy.stats import rankdata, norm, spearmanr, wilcoxon, mannwhitneyu, pearsonr
U = os.environ.get('PBISC_ROOT', r'F:\PBISC').rstrip('/\\') + os.sep  # folder holding eskd_covid_clinical_assoc_20261115/ and HSPC_.../
S = pd.read_csv(os.path.join(U,'eskd_covid_clinical_assoc_20261115','results','sample_level_composition_metadata.csv'))
F = pd.read_csv(os.path.join(U,'eskd_covid_clinical_assoc_20261115','results','sample_level_fine_label_pct.csv'))
OUT = os.path.join(U,'HSPC_CD56brightNK_COVID_analysis_summary_20261116','15_Figure7_ver8_Main_and_Supp','data') + os.sep
CLL = ['C116','C139']
cells = [c for c in F.columns if c not in ('source','sample_id','RBC','Platelet')]
X = F.set_index(['source','sample_id'])[cells].astype(float).fillna(0)
X = X.div(X.sum(axis=1), axis=0)*100

def table(src, drop=CLL):
    s = S[(S.source==src)&(S.condition=='positive')&(S.in_primary_scope==True)&(S.qc_pass_min200cells==True)&~S.patient_id.isin(drop)]
    v = X.loc[src].loc[s.sample_id].copy(); v['pid'] = s.patient_id.values
    T = v.groupby('pid').mean()
    meta = s.groupby('patient_id').agg(sev=('severity_ord','first'), fatal=('fatal_harmonized','first'), batch=('centre_or_wave','first'))
    return T.join(meta)

def inv_norm(x):  # rank-based INT, Blom-like (r-0.5)/n
    return pd.Series(norm.ppf((rankdata(x)-0.5)/len(x)), index=x.index)

def auc(y, d):  # explicit pairwise AUC
    y=np.asarray(y,float); d=np.asarray(d)==1; a=y[d]; b=y[~d]
    return ((a[:,None]>b[None,:]).sum()+0.5*(a[:,None]==b[None,:]).sum())/(len(a)*len(b))

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); q[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; return np.minimum(q,1)

R = table('Real'); B = table('PBISC'); N = B.loc[~B.index.isin(R.index)]; SH = R.index.intersection(B.index)
res = {'n': (len(R), int(R.fatal.sum()), len(N), int(N.fatal.sum()), len(SH))}
feats = [c for c in cells if (R[c]>0).mean()>=0.30]
rows=[]
for f in feats:
    E = pd.concat([pd.DataFrame({'y':inv_norm(R[f]),'sev':R.sev,'fatal':R.fatal,'batch':R.batch,'src':'R'}),
                   pd.DataFrame({'y':inv_norm(N[f]),'sev':N.sev,'fatal':N.fatal,'batch':N.batch,'src':'N'})])
    r = dict(feature=f)
    for k,D in [('Real',R.assign(y=R[f])),('New',N.assign(y=N[f])),('Ext',E)]:
        r['rho_'+k], r['p_sev_'+k] = spearmanr(D.sev, D.y) if D.y.nunique()>=4 else (np.nan,np.nan)
        r['AUC_'+k] = auc(D.y, D.fatal) if D.y.nunique()>=4 else np.nan
        r['p_death_'+k] = mannwhitneyu(D.y[D.fatal==1], D.y[D.fatal==0]).pvalue if D.y.nunique()>=4 else np.nan
    # centre/wave-adjusted (partial association via OLS residuals; design = batch + source dummies)
    Z = sm.add_constant(pd.get_dummies(E[['batch','src']], drop_first=True).astype(float))
    ey = sm.OLS(E.y.values, Z).fit().resid; es = sm.OLS(E.sev.values.astype(float), Z).fit().resid
    r['p_sev_adj'] = pearsonr(es, ey)[1]
    r['p_death_adj'] = mannwhitneyu(ey[E.fatal.values==1], ey[E.fatal.values==0]).pvalue
    rows.append(r)
T = pd.DataFrame(rows)
for k in ['Real','Ext']:
    T['q_sev_'+k] = bh(T['p_sev_'+k]); T['q_death_'+k] = bh(T['p_death_'+k])
s2s, s2d = T.q_sev_Ext<0.05, T.q_death_Ext<0.10
sg = lambda v: np.sign(v)
s3s = s2s & (sg(T.rho_Real)==sg(T.rho_New)) & (sg(T.rho_New)==sg(T.rho_Ext)) & (np.minimum(T.rho_Real.abs(),T.rho_New.abs())>=0.10)
dA = lambda k: T['AUC_'+k]-0.5
s3d = s2d & (sg(dA('Real'))==sg(dA('New'))) & (sg(dA('New'))==sg(dA('Ext'))) & (np.minimum(dA('Real').abs(),dA('New').abs())>=0.05)
nd = T.feature!='CD14mono_dex'
s4s, s4d = s3s & nd & (T.p_sev_adj<0.05), s3d & nd & (T.p_death_adj<0.05)
s5s, s5d = s4s & (T.q_sev_Real>=0.05), s4d & (T.q_death_Real>=0.10)
res['funnel'] = [len(T), int((s2s|s2d).sum()), int((s3s|s3d).sum()), int((s4s|s4d).sum()), int((s5s|s5d).sum())]
res['S5'] = sorted(T.feature[s5s|s5d])
# compare with pipeline outputs
P = pd.read_csv(OUT+'selection_funnel_43fine.csv').set_index('feature'); Ti = T.set_index('feature')
cmp = {}
for mine, theirs in [('rho_Real','rho_Real'),('rho_New','rho_New'),('rho_Ext','rho_Ext'),('AUC_Real','AUC_Real'),('AUC_New','AUC_New'),('AUC_Ext','AUC_Ext'),
                     ('q_sev_Real','q_sev_Real'),('q_sev_Ext','q_sev_Ext'),('q_death_Real','q_death_Real'),('q_death_Ext','q_death_Ext')]:
    cmp[mine] = float(np.nanmax(np.abs(Ti[mine]-P.loc[Ti.index, theirs])))
adj = P.p_sev_centre_adj.notna()
cmp['p_sev_adj (S3 rows)'] = float(np.nanmax(np.abs(Ti.loc[adj[adj].index,'p_sev_adj']-P.loc[adj,'p_sev_centre_adj'])))
cmp['p_death_adj (S3 rows)'] = float(np.nanmax(np.abs(Ti.loc[adj[adj].index,'p_death_adj']-P.loc[adj,'p_death_centre_adj'])))
res['max_abs_diff_vs_pipeline'] = cmp
# shared-donor Wilcoxon + modality
M = pd.read_csv(OUT+'modality_scale_shared27_and_raw_pct_v8.csv').set_index('feature')
res['wilcoxon'] = {f: (round(wilcoxon(R.loc[SH,f], B.loc[SH,f]).pvalue,9), round(float(M.loc[f,'wilcoxon_paired_p']),9)) for f in ['B_ASC_IgA','B_ASC_dividing','HSPC','NK CD56 bright','CD16mono','T CD8 CTL']}
iga_sign = (B.loc[SH,'B_ASC_IgA']>R.loc[SH,'B_ASC_IgA']).sum(); res['IgA PBISC>Real pairs'] = f'{iga_sign}/27'
# HSPC: Cambridge-only and within WHO2-3
E = pd.concat([pd.DataFrame({'y':inv_norm(R.HSPC),'sev':R.sev,'fatal':R.fatal,'batch':R.batch,'src':'R'}),pd.DataFrame({'y':inv_norm(N.HSPC),'sev':N.sev,'fatal':N.fatal,'batch':N.batch,'src':'N'})])
c = E[~((E.src=='R')&(E.batch=='NCL'))]; w = E[E.sev>=2]
res['HSPC'] = {'all': (round(auc(E.y,E.fatal),4), mannwhitneyu(E.y[E.fatal==1],E.y[E.fatal==0]).pvalue),
               'cambridge': (int(c.fatal.sum()), round(auc(c.y,c.fatal),4), mannwhitneyu(c.y[c.fatal==1],c.y[c.fatal==0]).pvalue),
               'WHO23': (int(w.fatal.sum()), round(auc(w.y,w.fatal),4), mannwhitneyu(w.y[w.fatal==1],w.y[w.fatal==0]).pvalue)}
# severity-stratified AUC (point estimate) + permutation with an independent seed
def strat(D):
    num=den=0
    for _,g in D.groupby('sev'):
        d=g.fatal==1
        if 0<d.sum()<len(g): num+=auc(g.y,g.fatal)*d.sum()*(~d).sum(); den+=d.sum()*(~d).sum()
    return num/den
obs = strat(E); rng=np.random.default_rng(777); hits=0; NP=4000
for _ in range(NP):
    Ep=E.copy(); Ep['fatal']=Ep.groupby('sev').fatal.transform(lambda s: rng.permutation(s.values)); hits += abs(strat(Ep)-.5)>=abs(obs-.5)
res['HSPC']['stratified'] = (round(obs,4), (hits+1)/(NP+1))
# power check (independent seed), n=20
def power(f, outcome, exp_dir, n=20, draws=1000, seed=99):
    rng=np.random.default_rng(seed); Nn=pd.DataFrame({'y':inv_norm(N[f]),'sev':N.sev,'fatal':N.fatal}); h=np.zeros(2); k=0
    while k<draws:
        s=R.loc[rng.choice(R.index,n,replace=False)]
        if outcome=='death' and not (0<s.fatal.sum()<n): continue
        k+=1; sr=pd.DataFrame({'y':inv_norm(s[f]),'sev':s.sev,'fatal':s.fatal})
        for i,D in enumerate([sr, pd.concat([sr,Nn])]):
            if outcome=='sev':
                r,p=spearmanr(D.sev,D.y); h[i]+= (p<.05)&(np.sign(r)==exp_dir)
            else:
                a=auc(D.y,D.fatal); p=mannwhitneyu(D.y[D.fatal==1],D.y[D.fatal==0]).pvalue; h[i]+=(p<.05)&(np.sign(a-.5)==exp_dir)
    return tuple(np.round(h/draws,3))
res['power_n20 (indep seed)'] = {k: power(*k) for k in [('B_ASC_IgG','sev',1),('B_ASC_IgA','sev',1),('T CD8 CM','death',-1),('HSPC','death',1)]}
import pprint; pprint.pprint(res, width=160)
