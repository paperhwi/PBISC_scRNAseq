"""Paired composition accuracy, donor-level uncertainty, and held-out nulls.

Primary truth: archived observed scANVI L1 fractions, same frozen labeling as
paired PBISC predictions. No trained donor is included for its adapter.
Secondary comparison: common adult sepsis samples and harmonized recap truth;
deconvolution and PBISC training/reference protocols remain different.
"""
from common import *
import json
from scipy.stats import spearmanr

rng = np.random.default_rng(20261014)
TYPES = ['T cell', 'B cell', 'Monocyte', 'NK cell', 'DC', 'Other']
metrics, errors, types, predictions = [], [], [], []
integrity = []

def score(cohort, scope, level, method, draw, truth, pred, meta, null=None):
    idx = truth.index.intersection(pred.index)
    truth, pred = truth.loc[idx], pred.loc[idx]
    assert np.allclose(truth.sum(axis=1), 1) and np.allclose(pred.sum(axis=1), 1)
    md = meta.loc[idx]
    l1 = np.abs(pred.to_numpy()-truth.to_numpy()).sum(1)
    ai = np.linalg.norm(clr(pred)-clr(truth), axis=1)
    e = pd.DataFrame({'sample_id': idx, 'donor_id': md.donor_id.to_numpy(), 'L1': l1, 'aitchison': ai})
    if null is not None: e['delta_L1_vs_null'] = l1-np.abs(null.loc[idx].to_numpy()-truth.to_numpy()).sum(1)
    donor = e.groupby('donor_id').mean(numeric_only=True)
    row = dict(cohort=cohort, scope=scope, level=level, method=method, draw=draw, n_samples=len(idx), n_donors=len(donor))
    for key in donor:
        row[key], row[key+'_lo'], row[key+'_hi'] = boot_mean(donor[key], rng)
    row['sample_weighted_L1'] = float(l1.mean())
    metrics.append(row)
    for _, r in e.iterrows(): errors.append(dict(cohort=cohort,scope=scope,level=level,method=method,draw=draw,**r.to_dict()))
    # Across-patient correlation separately by cell type; visits averaged first.
    td = truth.assign(donor_id=md.donor_id.to_numpy()).groupby('donor_id').mean()
    pd_ = pred.assign(donor_id=md.donor_id.to_numpy()).groupby('donor_id').mean()
    for c in truth.columns:
        r = spearmanr(td[c], pd_[c]).statistic if td[c].nunique()>1 and pd_[c].nunique()>1 else np.nan
        types.append(dict(cohort=cohort,scope=scope,level=level,method=method,draw=draw,cell_type=c,spearman=r,n_donors=len(td)))
    if level == 'l1' and scope in ('all_paired','adult_paired') and method.startswith('PBISC'):
        for s in idx:
            for c in truth.columns: predictions.append(dict(cohort=cohort,scope=scope,method=method,draw=draw,sample_id=s,donor_id=meta.loc[s,'donor_id'],cell_type=c,pred=pred.loc[s,c],truth=truth.loc[s,c]))

stored = {}
for coh in ['ye2026_sepsis','imperial_eskd']:
    source = {}
    for k in ['k1000','k4000']:
        p = OUT/'inputs'/k/'analysis/features'/coh
        source[k] = (pd.read_parquet(p/'composition.parquet'),pd.read_csv(p/'meta.csv').set_index('sample_id'))
    d, meta = source['k1000']; assert meta.index.is_unique
    assert source['k4000'][1].index.equals(meta.index)
    scopes = {'all_paired': meta.index[meta['set'].eq('paired')]}
    common_test = meta['set'].eq('paired') & meta[[f'paired_test_d{i}' for i in range(5)]].all(axis=1)
    scopes['common_test_paired'] = meta.index[common_test]
    if coh == 'ye2026_sepsis': scopes['adult_paired'] = meta.index[meta['set'].eq('paired') & meta.group.str.startswith('adult')]
    else:
        roster = set(pd.read_csv(REPORT/'results/covid_roster187_local_ids.csv').sample_id.astype(str))
        scopes['published187_paired'] = meta.index[meta['set'].eq('paired') & meta.index.isin(roster)]
        scopes['published187_active'] = meta.index[meta['set'].eq('paired') & meta.index.isin(roster) & meta.group.isin(['covid_positive','covid_negative'])]
    for level in ['l1','l2']:
        cats = sorted(d.loc[d.level.eq(level),'type'].unique())
        def pivot(df, set_, dr):
            return df[df.level.eq(level) & df['set'].eq(set_) & df.draw.eq(dr)].pivot(index='sample_id',columns='type',values='frac').reindex(columns=cats,fill_value=0).fillna(0)
        obs = pivot(d,'observed',-1)
        assert np.allclose(obs, pivot(source['k4000'][0],'observed',-1).loc[obs.index])
        byk = {}
        for k in source:
            byk[k] = {}
            for dr in range(5):
                gen = pivot(source[k][0],'paired_gen',dr)
                assert set(gen.index) <= set(obs.index)
                assert meta.loc[gen.index,f'paired_test_d{dr}'].astype(bool).all()
                assert not meta.loc[gen.index,f'donor_trained_d{dr}'].astype(bool).any()
                train = obs.index[meta.loc[obs.index,f'donor_trained_d{dr}'].astype(bool)]
                train_patients = set(meta.loc[train,'donor_id'])
                assert len(train_patients)==5 and not set(meta.loc[gen.index,'donor_id']) & train_patients
                train_mean = obs.loc[train].assign(patient=meta.loc[train,'donor_id'].values).groupby('patient').mean().mean()
                null = pd.DataFrame(np.tile(train_mean,(len(gen),1)),index=gen.index,columns=cats)
                byk[k][dr] = gen
                for scope, ids in scopes.items():
                    ids = gen.index.intersection(ids)
                    if not len(ids): continue
                    score(coh,scope,level,f'PBISC_{k}',dr,obs.loc[ids],gen.loc[ids],meta,null)
                    if k=='k1000': score(coh,scope,level,'training5_null',dr,obs.loc[ids],null.loc[ids],meta)
                integrity.append(dict(cohort=coh,K=k,level=level,draw=dr,test_samples=len(gen),train_donors=len(train_patients),train_test_patient_overlap=0))
        for scope, ids in scopes.items():
            ids = obs.index.intersection(ids)
            # Cross-fitted optimistic null: every other observed patient's mean.
            om = obs.loc[ids].assign(patient=meta.loc[ids,'donor_id'].values).groupby('patient').mean()
            loo = pd.DataFrame([(om.sum()-om.loc[meta.loc[s,'donor_id']])/(len(om)-1) for s in ids],index=ids,columns=cats)
            score(coh,scope,level,'other_patient_null',-1,obs.loc[ids],loo,meta)
            ens = {}
            for k in byk:
                ens[k] = pd.concat(byk[k].values()).groupby(level=0).mean().reindex(ids)
                assert not ens[k].isna().any().any()
                score(coh,scope,level,f'PBISC_{k}_eligible_ensemble',-1,obs.loc[ids],ens[k],meta,loo)
            ek = np.abs(ens['k4000']-obs.loc[ids]).sum(axis=1)-np.abs(ens['k1000']-obs.loc[ids]).sum(axis=1)
            dd = ek.groupby(meta.loc[ids,'donor_id']).mean()
            mean,lo,hi = boot_mean(dd,rng)
            metrics.append(dict(cohort=coh,scope=scope,level=level,method='paired_K4_minus_K1',draw=-1,n_samples=len(ids),n_donors=len(dd),L1=mean,L1_lo=lo,L1_hi=hi))
            if level=='l1':
                stored[(coh,scope)] = (obs.loc[ids],ens,meta)
                # Does own-donor pairing help after removing clinical-group means?
                cols = cats
                t = obs.loc[ids].copy(); p = ens['k4000'].copy()
                index = pd.MultiIndex.from_arrays([meta.loc[ids,'donor_id'],meta.loc[ids,'group']],names=['donor','group'])
                t.index=index;p.index=index;t=t.groupby(level=[0,1]).mean();p=p.groupby(level=[0,1]).mean()
                tr=t-t.groupby(level='group').transform('mean');pr=p-p.groupby(level='group').transform('mean')
                stat=float(np.mean((tr.to_numpy()-pr.to_numpy())**2))
                # Patient blocks, not visits: each donor has a vector of group-specific residuals.
                # Restrict to patients with exactly one group to keep permutation exchanges valid.
                single=t.index.get_level_values('donor').value_counts();keep=single[single==1].index
                tr=tr[tr.index.get_level_values('donor').isin(keep)];pr=pr.loc[tr.index]
                stat=float(np.mean((tr.to_numpy()-pr.to_numpy())**2));perms=[]
                for _ in range(1000):
                    q=pr.copy()
                    for group in tr.index.get_level_values('group').unique():
                        loc=np.where(tr.index.get_level_values('group')==group)[0];q.iloc[loc]=pr.iloc[rng.permutation(loc)].to_numpy()
                    perms.append(float(np.mean((tr.to_numpy()-q.to_numpy())**2)))
                metrics.append(dict(cohort=coh,scope=scope,level=level,method='K4_within_group_pairing_permutation',draw=-1,n_samples=len(tr),n_donors=len(keep),residual_MSE=stat,permuted_MSE=float(np.mean(perms)),p_pairing=(1+np.sum(np.asarray(perms)<=stat))/1001))
    print('paired complete',coh,flush=True)

# Align adult sepsis evaluation samples and truth, keeping training differences explicit.
truth = pd.read_csv(REPORT/'proportion_qualification/results/sepsis_sc_truth_props.csv.gz')
truth = truth[truth.level.eq('major') & truth.grp.str.startswith('adult')].copy()
truth['type']=truth.label.map(benchmark6)
truth = truth.groupby(['sample','type']).frac.sum().unstack(fill_value=0).reindex(columns=TYPES,fill_value=0)
dec = pd.read_csv(REPORT/'proportion_qualification/results/sepsis_deconv_props.csv.gz')
dec = dec[dec['set'].eq('paired') & dec.level.eq('major')].copy();dec['type']=dec.label.map(benchmark6)
obs,ens,meta=stored[('ye2026_sepsis','adult_paired')]
ids = truth.index.intersection(obs.index).intersection(dec['sample'].unique())
for k in ens:score('ye2026_sepsis','common60_recap_truth','harmonized6',f'PBISC_{k}_eligible_ensemble',-1,truth.loc[ids],ens[k].reindex(columns=TYPES).loc[ids],meta)
for method,dm in dec.groupby('method'):
    p=dm.groupby(['sample','type']).frac.sum().unstack(fill_value=0).reindex(index=ids,columns=TYPES,fill_value=0).fillna(0);p=p.div(p.sum(axis=1),axis=0)
    score('ye2026_sepsis','common60_recap_truth','harmonized6',method,-1,truth.loc[ids],p,meta)
decm=json.load(open(REPORT/'proportion_qualification/results/sepsis_deconv_meta.json'))
null=pd.DataFrame(index=ids,columns=TYPES,dtype=float)
for fold in decm['folds']:
    ev=ids.intersection(fold);ref=truth.index.difference(fold);null.loc[ev]=truth.loc[ref].mean().values
score('ye2026_sepsis','common60_recap_truth','harmonized6','deconv_fold_null',-1,truth.loc[ids],null,meta)
assert null.notna().all().all()
for K in ['k1000','k4000']:
    raw=pd.read_parquet(OUT/'inputs'/K/'analysis/features/ye2026_sepsis/composition.parquet')
    for draw in range(5):
        q=raw[(raw['set']=='paired_gen')&(raw.level=='l1')&(raw.draw==draw)].pivot(index='sample_id',columns='type',values='frac').reindex(columns=TYPES,fill_value=0).fillna(0)
        ii=ids.intersection(q.index)
        score('ye2026_sepsis','common60_recap_truth','harmonized6',f'PBISC_{K}',draw,truth.loc[ii],q.loc[ii],meta,null.loc[ii])
pd.DataFrame(metrics).to_csv(OUT/'results/paired_accuracy_summary.csv',index=False)
pd.DataFrame(errors).to_csv(OUT/'results/paired_sample_errors.csv.gz',index=False)
pd.DataFrame(types).to_csv(OUT/'results/paired_type_correlations.csv',index=False)
pd.DataFrame(predictions).to_csv(OUT/'results/paired_predictions.csv.gz',index=False)
json.dump(integrity,open(OUT/'results/paired_integrity.json','w'),indent=2)
print('All paired outputs saved',flush=True)
