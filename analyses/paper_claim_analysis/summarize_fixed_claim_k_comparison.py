from pathlib import Path
import pandas as pd, numpy as np, json
from scipy.stats import bootstrap
ROOT=Path('paper_claim_analysis')
OUT=ROOT/'fixed_claim_k_paired_contrasts.csv'
rng=np.random.default_rng(92107)
rows=[]

def summarize(d, cohort, group_col, case, control, keys, metric_map, patient_map=None, draws=range(5)):
    d=d[d[group_col].isin([case,control])].copy()
    if patient_map is not None:
        d['patient_id']=d.sample_id.astype(str).str.split('_').str[0]
        unit='patient_id'
    else: unit='sample_id'
    for claim, metric in metric_map:
        x=d[d['claim_id']==claim]
        if x.empty or metric not in x or x[metric].notna().sum()==0: continue
        # per-draw group effects; average draws as descriptive summary
        effects={}; pairedK={}
        for K in [1000,4000]:
            effects[K]={}
            for draw in draws:
                z=x[(x.K==K)&(x.draw==draw)]
                means=z.groupby([group_col,unit],observed=True)[metric].mean().groupby(level=0).mean()
                if case in means and control in means: effects[K][draw]=means[case]-means[control]
        valid=sorted(set(effects[1000])&set(effects[4000]))
        if not valid: continue
        delta_draw={j: effects[4000][j]-effects[1000][j] for j in valid}
        # Bootstrap samples (Ye) or patients (COVID), separately within each group.
        agg=x.groupby([group_col,unit,'K','draw'],observed=True)[metric].mean().reset_index()
        combo=[(K,j) for K in [1000,4000] for j in valid]
        arrays={}
        units_by_group={}
        for g in [case,control]:
            z=agg[agg[group_col]==g].pivot(index=unit,columns=['K','draw'],values=metric).reindex(columns=pd.MultiIndex.from_tuples(combo))
            z=z.dropna(axis=0,how='any')
            arrays[g]=z.to_numpy(); units_by_group[g]=z.index.to_list()
        reps=[]
        for b in range(1000):
            means={}
            for g in [case,control]:
                a=arrays[g]; sampled=rng.integers(0,len(a),size=len(a)); means[g]=np.mean(a[sampled],axis=0)
            contrast=means[case]-means[control]
            reps.append((np.mean(contrast[:len(valid)]),np.mean(contrast[len(valid):])))
        reps=np.asarray(reps)
        eff1=float(np.mean(list(effects[1000].values()))); eff4=float(np.mean(list(effects[4000].values())))
        diff=eff4-eff1
        lo,hi=np.quantile(reps[:,1]-reps[:,0],[.025,.975]) if len(reps) else (np.nan,np.nan)
        rows.append(dict(cohort=cohort,claim_id=claim,metric=metric,contrast=f'{case} - {control}',n_units_case=len(units_by_group[case]),n_units_control=len(units_by_group[control]),n_draws=len(valid),K1000_effect=eff1,K4000_effect=eff4,K4000_minus_K1000=diff,paired_K_boot_CI_low=lo,paired_K_boot_CI_high=hi,draw_effects_K1000=';'.join(f'{j}:{effects[1000][j]:.4f}' for j in valid),draw_effects_K4000=';'.join(f'{j}:{effects[4000][j]:.4f}' for j in valid),draw_deltas=';'.join(f'{j}:{delta_draw[j]:.4f}' for j in valid),draws_same_delta_sign=sum(np.sign(v)==np.sign(diff) for v in delta_draw.values()),bootstrap_unit=unit))

s=pd.read_csv(ROOT/'fixed_claim_pbisc_sepsis_by_sample_draw.csv')
for c in s.claim_id.unique():
    for metric in ['fraction_marker_positive','marker_score_mean']:
        summarize(s,'Ye sepsis','group','adult_sepsis','adult_HC',None,[(c,metric)])
# COVID patient map from published roster reconciliation
pm=pd.read_csv(ROOT/'covid_sample_patient_id_reconciliation.csv')[['sample_id','patient_id']].drop_duplicates()
c=pd.read_csv(ROOT/'fixed_claim_covid_pbisc_by_sample_draw.csv')
c=c.rename(columns={'fraction_signature_positive':'fraction_marker_positive'})
for cl in c.claim_id.unique():
    for metric in ['fraction_marker_positive','mean_IFN_score','mean_DEX_up_score','mean_DEX_down_score']:
        summarize(c,'Imperial COVID','group','covid_positive','covid_negative',None,[(cl,metric)],patient_map=pm)
pd.DataFrame(rows).to_csv(OUT,index=False)
print('saved',OUT,'rows',len(rows))
print(pd.DataFrame(rows).to_string(index=False))
# Provenance compare all draws for each cohort / K from H5AD uns JSON
import anndata as ad
prov=[]
for cohort, dirs, pats in [('sepsis',[Path('bulkonly_lora5_20261006/ye2026_sepsis'),Path('bulkonly_lora5_k4000_20261006/ye2026_sepsis')],['lora5_bulkonly_ye2026_sepsis_draw','lora5k4_bulkonly_ye2026_sepsis_draw']),('covid',[Path('bulkonly_lora5_20261006/imperial_eskd'),Path('bulkonly_lora5_k4000_20261006/imperial_eskd')],['lora5_bulkonly_imperial_eskd_draw','lora5k4_bulkonly_imperial_eskd_draw'])]:
  for ki,(dr,prefix) in enumerate(zip(dirs,pats)):
    for j in range(5):
      f=dr/f'{prefix}{j}.h5ad'; a=ad.read_h5ad(f,backed='r'); p=json.loads(a.uns['provenance']); a.file.close()
      prov.append({'cohort':cohort,'K':1000 if ki==0 else 4000,'draw':j,'training_donors':';'.join(p.get('training_donors',[])),'checkpoint_sha256':p.get('checkpoint_sha256'),'generation':p.get('generation'),'seed':p.get('seed'),'model':p.get('model')})
pd.DataFrame(prov).to_csv(ROOT/'fixed_claim_lora_provenance.csv',index=False)
print('\nProvenance first rows')
print(pd.DataFrame(prov).to_string(index=False))
