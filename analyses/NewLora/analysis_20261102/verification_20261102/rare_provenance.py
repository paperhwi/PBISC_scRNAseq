"""Recount rare subtype cells from saved transferred labels."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/code'))
from rare_lib import SEPSIS,COVID
OUT=HERE.parent
group_map={'sepsis':{'adult_sepsis':'case','adult_HC':'control','adult_ICU_ctrl':'ICU control'},
           'covid':{'covid_positive':'case','covid_negative':'control'}}
unit=pd.read_csv(OUT/'rare_populations/rare_unit_shares_real_vs_ceiling.csv')
pooled=pd.read_csv(OUT/'figure_realism_rare/rare_parent_share_with_new_lora.csv')
checks=[]
for cohort,defs in [('sepsis',SEPSIS),('covid',COVID)]:
    lab=pd.read_csv(OUT/cohort/'ceiling_cell_labels.csv.gz',usecols=['sample','group','major_tr','fine_tr'])
    lab['condition']=lab.group.map(group_map[cohort])
    lab['patient']=lab['sample'] if cohort=='sepsis' else lab['sample'].str.split('_').str[0]
    for subset,fun in defs.items():
        major=lab.major_tr.astype(str).to_numpy();fine=lab.fine_tr.astype(str).to_numpy();group=lab.group.astype(str).to_numpy()
        parent=fun(major,fine);target=parent&(fine==subset)
        for grp in ('all', 'adult_sepsis' if cohort=='sepsis' else 'covid_positive', 'adult_HC' if cohort=='sepsis' else 'covid_negative'):
            mask=np.ones(len(lab),bool) if grp=='all' else group==grp
            old=pooled[(pooled.cohort==cohort)&(pooled.subset==subset)&(pooled.source=='new LoRA ceiling (transferred)')&(pooled.group==grp)]
            checks.append({'cohort':cohort,'section':'pooled','subset':subset,'group':grp,
                           'passed':len(old)==1 and old.n_target.iloc[0]==int((target&mask).sum()) and old.n_parent.iloc[0]==int((parent&mask).sum())})
        sub=lab[lab.condition.notna()].copy();sub['parent']=parent[lab.condition.notna().to_numpy()];sub['target']=target[lab.condition.notna().to_numpy()]
        calc=sub.groupby(['patient','condition']).agg(n_cells=('target','size'),n_subset=('target','sum'),n_parent=('parent','sum')).reset_index()
        old=unit[(unit.cohort==cohort)&(unit.source=='ceiling')&(unit.subset==subset)].rename(columns={'group':'condition'})
        match=calc.merge(old,on=['patient','condition'],how='outer',suffixes=('_calc','_out'),indicator=True)
        good=(len(match)==len(calc)==len(old) and match['_merge'].eq('both').all() and all((match[f'{c}_calc']==match[f'{c}_out']).all() for c in ('n_cells','n_subset','n_parent')))
        checks.append({'cohort':cohort,'section':'patient','subset':subset,'group':'all','passed':good})
checks=pd.DataFrame(checks)
checks.to_csv(HERE/'rare_provenance_checks.csv',index=False)
print(checks.groupby(['cohort','section']).passed.agg(['sum','count']).to_string())
print('failed',int((~checks.passed).sum()))
