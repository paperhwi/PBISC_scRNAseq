from pathlib import Path
import pandas as pd
import numpy as np
import json

base = Path('F:/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012')
report = {}
pb = base/'ilisi_floor_20261021/data/pb'
for f in pb.glob('*covid_obs.csv'):
    df = pd.read_csv(f)
    report[f.name] = {'rows':len(df),'columns':list(df.columns)}
    if 'sample' in df:
        ids=df['sample'].astype(str)
        donors=ids.str.split('_').str[0]
        train={'C126','C141','C58'}
        report[f.name]['samples_from_named_training_donors'] = sorted(ids[donors.isin(train)].unique().tolist())
    if 'patient' in df and 'case_control' in df:
        d=df[df.case_control.isin(['POSITIVE','NEGATIVE'])]
        both=d.groupby('patient').case_control.nunique()
        report[f.name]['patients_in_both_conditions']=both[both>1].index.astype(str).tolist()

a=pd.read_csv(base/'week1_defensive_20261027/results/a3_power_matched_ceiling.csv')
report['ceiling_definition_comparison'] = a[['cohort','celltype','floor_r_full','floor_r_matchedN','ceil_case_split_r','pct_of_ceiling_r']].assign(full_n_pct=lambda d:100*d.floor_r_full/d.ceil_case_split_r).round(4).to_dict('records')

# Demonstrate the selection semantics used by the saved bootstrap code.
toy=pd.DataFrame({'patient':['A','B','C','D']})
resample=['A','A','B','B']
report['bootstrap_selection_example']={'draw':resample,'isin_selected':toy.loc[toy.patient.isin(resample),'patient'].tolist(),'correct_row_selection':toy.set_index('patient').loc[resample].reset_index().patient.tolist()}

notebooks=[]
for directory in ['disease_deg_rare_20261027','week1_defensive_20261027','tcell_batch_floor_20261027']:
    for p in (base/directory).glob('*.ipynb'):
        nb=json.loads(p.read_text(encoding='utf-8'))
        code=[c for c in nb.get('cells',[]) if c.get('cell_type')=='code']
        notebooks.append({'path':str(p),'code_cells':len(code),'executed_cells':sum(c.get('execution_count') is not None for c in code),'stored_errors':[o.get('ename') for c in code for o in c.get('outputs',[]) if o.get('output_type')=='error']})
report['notebook_saved_execution']=notebooks
out=Path('F:/PBISC/project_review_20261028/table_audit.json')
out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
