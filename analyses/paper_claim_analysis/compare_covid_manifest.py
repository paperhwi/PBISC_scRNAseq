import pandas as pd
from pathlib import Path
base=Path('paper_claim_analysis')
local=pd.read_csv(base/'covid_sample_patient_manifest_local.csv',dtype=str).fillna('')
pub=pd.read_excel(base/'Stephenson_Table_S1_sample_metadata.xlsx',header=2,dtype={'sample_id':str,'patient_id':str})
def key(s):
 s=str(s).strip(); pre,sep,suf=s.partition('_')
 return pre+'_'+suf[:4] if sep and len(suf)>=4 else s
local['sample_key']=local.sample_id.map(key); pub['sample_key']=pub.sample_id.map(key)
print('local dup keys',local[local.sample_key.duplicated(False)][['sample_id','sample_key','patient_id']].to_string(index=False))
print('pub dup keys',pub[pub.sample_key.duplicated(False)][['sample_id','sample_key','patient_id']].to_string(index=False))
L=set(local.sample_key); P=set(pub.sample_key)
print('local not pub',local[~local.sample_key.isin(P)][['sample_id','patient_id','centre','case_control','WHO_severity','sample_date','n_cells']].to_string(index=False))
print('pub not local',pub[~pub.sample_key.isin(L)][['sample_id','patient_id','case_control','year_of_sample']].to_string(index=False))
LP=set(local.patient_id); PP=set(pub.patient_id)
print('local patient not pub',sorted(LP-PP)); print('pub patient not local',sorted(PP-LP))
print('local samples',local.shape[0],'pub samples',len(pub),'local patients',len(LP),'pub patients',len(PP))
