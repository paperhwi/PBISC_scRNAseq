import h5py, numpy as np, pandas as pd, os
p=r'bulkonly_lora5_20261006/analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad'
f=h5py.File(p,'r'); obs=f['obs']; idx=obs['_index'][:]
def col(name):
 o=obs[name]
 if isinstance(o,h5py.Group):
  cats=o['categories'][:]; cats=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in cats]); codes=o['codes'][:].astype(int)
  vals=np.array([cats[x] if x>=0 else '' for x in codes])
 else:
  vals=o[:]
  if vals.dtype.kind in 'SO': vals=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in vals])
 return vals
cols=['sample_id','patient_id','centre','case_control','WHO_severity','sample_date','date_first_symptoms','date_positive_swab','days_since_symptoms','exclude','orig.ident']
df=pd.DataFrame({c:col(c) for c in cols})
print('rows',len(df),'excluded=',pd.Series(df.exclude).value_counts(dropna=False).to_dict())
print('cross-tab centre/case_control\n',df.drop_duplicates('sample_id').groupby(['centre','case_control'],dropna=False).sample_id.nunique().to_string())
sm=df.groupby(['sample_id','patient_id','centre','case_control','WHO_severity','sample_date','date_first_symptoms','date_positive_swab','exclude'],dropna=False).size().reset_index(name='n_cells')
sm.to_csv(r'paper_claim_analysis/covid_sample_patient_manifest_local.csv',index=False)
print('samples',len(sm),'patients',sm.patient_id.nunique())
print('patients per centre and control',sm.groupby(['centre','case_control']).patient_id.nunique().to_string())
print('samples per centre and control',sm.groupby(['centre','case_control']).sample_id.nunique().to_string())
print('first rows\n',sm.head(20).to_string(index=False))
print('exclude sample aggregates\n',sm.groupby('exclude').agg(samples=('sample_id','nunique'),patients=('patient_id','nunique'),cells=('n_cells','sum')).to_string())
print('patient count per case control\n',sm.drop_duplicates('patient_id').groupby(['centre','case_control']).patient_id.nunique().to_string())
f.close()
