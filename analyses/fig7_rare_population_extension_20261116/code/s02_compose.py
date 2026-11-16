import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
U='/mnt/user-data/uploads/PBISC/'; D=U+'eskd_covid_clinical_assoc_20261115/data/'
RARE13=['ASDC','B_ASC_IgG','B_ASC_dividing','CD14mono_IFN','CD14mono_dex','DC1','HSPC','ILC','MAIT','NK CD56 bright','T CD4 cycling','T g/d','pDC']
B5={'T cell':'T','B cell':'B','Monocyte':'Mono','NK cell':'NK','DC':'Other','Other':'Other'}
lm=pd.read_csv(U+'deliverables/annotated_h5ad_20261112/tables/label_mapping_covid.csv'); FINE=list(lm.fine)
def comp(cells,sid):
    n=cells.groupby(sid,observed=True).size().rename('n_cells')
    b=pd.crosstab(cells[sid],cells.broad6.map(B5)).reindex(columns=['T','B','Mono','NK','Other'],fill_value=0)
    dc=(cells.broad6=='DC').groupby(cells[sid],observed=True).sum().rename('DC')
    f=pd.crosstab(cells[sid],cells.fine).reindex(columns=FINE,fill_value=0)
    out=pd.concat([n,100*b.div(n,axis=0).add_suffix('_pct'),(100*dc/n).rename('DC_within_Other_pct'),
                   100*f[RARE13].div(n,axis=0).add_prefix('rare_').add_suffix('_pct')],axis=1)
    out['rare13_total_pct']=100*f[RARE13].sum(1)/n
    fine=100*f.div(n,axis=0)
    return out,fine
# ---- Real ----
r=pd.read_csv(D+'real_obs_full.csv.gz',low_memory=False); r['fine']=r.cell_type_IFNstim_recap; r['broad6']=r.fine.map(lm.set_index('fine').six_type_broad)
cR,fR=comp(r,'sample_id')
mc=['patient_id','case_control','WHO_severity','WHO_temp_severity','fatal_disease','calc_age','sex','ethnicity','diabetes','ihd','copd','previous_vte','smoking','cause_eskd','centre','days_since_posswab','days_since_symptoms','temp_CRP','temp_WCC','temp_lymph','temp_mono','temp_neut','temp_ddimer','temp_ferritin','temp_trop','sample_date']
mR=r.groupby('sample_id')[mc].first()
assert (r.groupby('sample_id')[mc].nunique(dropna=False).max()<=1).all()
mR=mR.rename(columns={'case_control':'condition','calc_age':'age','days_since_posswab':'days_from_first_positive_swab','days_since_symptoms':'days_from_first_symptoms','centre':'centre_or_wave'})
mR['condition']=mR.condition.str.lower(); mR['fatal_raw']=mR.fatal_disease
mR['fatal']=mR.fatal_disease.map({'yes':1,'indirect':1,'no':0}); mR['fatal_direct_only']=mR.fatal_disease.map({'yes':1,'indirect':0,'no':0})
mR['source']='Real'
# ---- PBISC ----
o=pd.read_parquet('cache/pbisc_obs.parquet'); o['sample_id']=o.sample_id.astype(str)
cP,fP=comp(o,'sample_id')
pc=['donor_id','condition','WHO_severity','WHO_temp_severity','fatal_disease','calc_age','sex','ethnicity','diabetes','ihd','copd','smoking','cause_eskd','wave','time_from_first_positive_swab','time_from_first_symptoms','in_primary_scope']
mP=o.groupby('sample_id',observed=True)[pc].first()
w=pd.concat([pd.read_csv(U+'bulkonly_lora5_20261006/analysis/bulk_public/zenodo6497251_imperial_eskd/w%d_metadata.csv'%k) for k in (1,2)]).drop_duplicates('sample_id').set_index('sample_id')
mP['previous_vte']=w.previous_vte.reindex(mP.index)
mP=mP.rename(columns={'donor_id':'patient_id','calc_age':'age','wave':'centre_or_wave','time_from_first_positive_swab':'days_from_first_positive_swab','time_from_first_symptoms':'days_from_first_symptoms'})
mP['condition']=mP.condition.astype(str).map({'covid_positive':'positive','covid_negative':'negative','recovery':'recovery'})
for c in ['WHO_severity','WHO_temp_severity']: mP[c]=mP[c].astype(str).replace({'':None,'nan':None,'NEGATIVE':None,'RECOVERY':None})
fd=mP.fatal_disease.astype(str); mP['fatal_raw']=fd; mP['fatal']=fd.map({'True':1,'False':0}); mP['fatal_direct_only']=mP['fatal']
mP.loc[mP.condition!='positive',['fatal','fatal_direct_only']]=np.nan if False else mP.loc[mP.condition!='positive',['fatal','fatal_direct_only']]
mP['source']='PBISC'
print('PBISC sev values',mP.WHO_severity.unique(), mP.WHO_temp_severity.unique(), mP.fatal_raw.unique())
S=pd.concat([mR.join(cR),mP.join(cP)]).reset_index().rename(columns={'index':'sample_id'})
S['in_primary_scope']=S.in_primary_scope.fillna(True).astype(bool)
for c in ['days_from_first_positive_swab','days_from_first_symptoms']: S[c]=pd.to_numeric(S[c],errors='coerce')
sevmap={'mild':0,'moderate':1,'severe':2,'critical':3}
S['severity_ord']=S.WHO_severity.map(sevmap); S['temp_severity_ord']=S.WHO_temp_severity.map(sevmap)
S['severity_bin']=S.severity_ord.map(lambda x: np.nan if pd.isna(x) else ('severe/critical' if x>=2 else 'mild/moderate'))
# same-day matched Real<->PBISC sample
S['day_key']=S.patient_id.astype(str)+'|'+S.condition+'|'+S.days_from_first_positive_swab.round().astype('Int64').astype(str)
k=S[S.condition=='positive'].groupby('day_key').source.nunique(); S['same_day_pair']=S.day_key.map(k).eq(2)&S.condition.eq('positive')
pat_src=S.groupby('patient_id').source.nunique(); S['patient_in_both_sources']=S.patient_id.map(pat_src).eq(2)
front=['source','sample_id','patient_id','condition','WHO_severity','severity_ord','severity_bin','WHO_temp_severity','temp_severity_ord','fatal','fatal_raw','fatal_direct_only','days_from_first_positive_swab','days_from_first_symptoms','age','sex','ethnicity','diabetes','ihd','copd','previous_vte','smoking','cause_eskd','centre_or_wave','in_primary_scope','patient_in_both_sources','same_day_pair','n_cells']
pct=[c for c in S.columns if c.endswith('_pct')]
labs=['temp_CRP','temp_WCC','temp_lymph','temp_mono','temp_neut','temp_ddimer','temp_ferritin','temp_trop']
S=S[front+pct+labs+['sample_date','day_key']]
for c in labs: S[c]=pd.to_numeric(S[c],errors='coerce')
S['age']=pd.to_numeric(S.age,errors='coerce')
S.to_csv('out/sample_level_composition_metadata.csv',index=False)
F=pd.concat([fR.assign(source='Real'),fP.assign(source='PBISC')]).reset_index().rename(columns={'index':'sample_id'})
F=F[['source','sample_id']+FINE]; F.to_csv('out/sample_level_fine_label_pct.csv',index=False)
# ---- patient level (COVID-positive samples averaged; controls = negative samples) ----
S['qc_pass_min200cells']=S.n_cells>=200
fh=S[S.condition=='positive'].groupby('patient_id').fatal.max(); S['fatal_harmonized']=np.where(S.condition=='positive',S.patient_id.map(fh),np.nan)
S.to_csv('out/sample_level_composition_metadata.csv',index=False)
Sp=S[S.in_primary_scope&S.qc_pass_min200cells]
def agg(g):
    d=g.iloc[0][['source','patient_id','WHO_severity','severity_ord','severity_bin','fatal','fatal_harmonized','fatal_raw','fatal_direct_only','age','sex','ethnicity','diabetes','ihd','copd','previous_vte','smoking','cause_eskd','centre_or_wave','patient_in_both_sources']].to_dict()
    d.update(n_samples=len(g),total_cells=g.n_cells.sum(),first_sample_day=g.days_from_first_positive_swab.min(),last_sample_day=g.days_from_first_positive_swab.max(),
             max_temp_severity_ord=g.temp_severity_ord.max(),max_CRP=g.temp_CRP.max(),min_lymph=g.temp_lymph.min())
    d.update(g[pct].mean().to_dict()); return pd.Series(d)
pos=Sp[Sp.condition=='positive'].groupby(['source','patient_id'],sort=False).apply(agg).reset_index(drop=True); pos.insert(2,'group','COVID_positive')
neg=Sp[Sp.condition=='negative'].groupby(['source','patient_id'],sort=False).apply(agg).reset_index(drop=True); neg.insert(2,'group','COVID_negative_control')
for c in ['WHO_severity','severity_ord','severity_bin','fatal','fatal_harmonized','fatal_raw','fatal_direct_only','max_temp_severity_ord']: neg[c]=np.nan
P=pd.concat([pos,neg]).sort_values(['group','patient_id','source']); P.to_csv('out/patient_level_composition_metadata.csv',index=False)
print('severity/fatal unique within patient:',S[S.condition=='positive'].groupby(['source','patient_id'])[['WHO_severity','fatal']].nunique().max().to_dict())
print(S.groupby(['source','condition']).agg(samples=('sample_id','size'),patients=('patient_id','nunique')))
print(P.groupby(['group','source']).size())
pp=P[P.group=='COVID_positive']
print(pd.crosstab([pp.source],[pp.WHO_severity,pp.fatal]))
print('union positive patients',pp.patient_id.nunique(),'both',pp.groupby('patient_id').source.nunique().eq(2).sum(),'same-day pairs',S[S.same_day_pair].day_key.nunique())
# metadata consistency for patients in both
w_=pp.pivot(index='patient_id',columns='source',values=['WHO_severity','fatal','age','sex']).dropna()
for c in ['WHO_severity','fatal','age','sex']: print(c,'mismatch',(w_[c]['Real']!=w_[c]['PBISC']).sum(), list(w_.index[w_[c]['Real']!=w_[c]['PBISC']]))
print(S.groupby('source')[['T_pct','B_pct','Mono_pct','NK_pct','Other_pct','rare13_total_pct']].mean().round(2))
print(S.groupby('source')[[c for c in pct if c.startswith('rare_')]].mean().round(3).T)
print(S[S.source=='Real'].n_cells.describe())
