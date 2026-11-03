"""Verify acquired counts, patient isolation and prepared intervention inputs."""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1]
PREV=ROOT/'bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014'
sys.path.insert(0,str(PREV/'code'))
from common import read_csr_rows
import numpy as np,pandas as pd,h5py,anndata as ad
checks=[]
for rec in json.load(open(OUT/'results/real_counts_manifest.json')):
    with h5py.File(rec['path'],'r') as f:obs=ad.io.read_elem(f['obs']);var=ad.io.read_elem(f['var'])
    with h5py.File(rec['original_source'],'r') as f:rawvar=ad.io.read_elem(f['var'])
    if rec['cohort']=='ye2026_sepsis':rk=rawvar.gene_ids.astype(str).str.split('.').str[0];key=var.gene_id.astype(str).str.split('.').str[0]
    else:rk=rawvar.index.astype(str);key=var.gene_symbol.astype(str)
    mp=pd.Series(np.arange(len(rawvar)),index=rk);mp=mp[~mp.index.duplicated(keep=False)]
    rows=np.unique(np.linspace(0,len(obs)-1,12).astype(int));src=obs.source_row.iloc[rows].to_numpy()
    xx=read_csr_rows(rec['path'],rows,chunk=256);yy=read_csr_rows(rec['original_source'],src,mp.loc[key].to_numpy(),chunk=256)
    assert (xx!=yy).nnz==0
    roster=pd.read_csv(OUT/f"results/{rec['cohort']}_strict_heldout_roster.csv")
    assert set(obs.sample_id)==set(roster.sample_id)
    assert not obs.patient_id.isna().any() and obs.patient_id.nunique()==rec['patients']
    for j in range(5):assert not roster[f'donor_trained_d{j}'].astype(str).str.lower().isin(['true','1']).any()
    checks.append(dict(cohort=rec['cohort'],counts_equal_original_rows=len(rows),patients=rec['patients'],samples=rec['samples'],adapter_train_patient_overlap=0,pretraining_patient_overlap='not verifiable without base-model training roster'))
s=pd.read_csv(OUT/'results/bulk_swap_assignments.csv');ok=s.status.eq('assignment_ready_inference_not_run')
assert (s.loc[ok,'target_patient']!=s.loc[ok,'source_patient']).all()
assert s.loc[ok].groupby(['cohort','target_sample']).size().eq(3).all()
for cohort,m in s.loc[ok].groupby('cohort'):
    r=pd.read_csv(OUT/f'results/{cohort}_strict_heldout_roster.csv').set_index('sample_id')
    cols=['group'] if cohort=='ye2026_sepsis' else ['group','wave','severity_peak']
    for col in cols:assert np.array_equal(r.loc[m.target_sample,col].fillna('missing').to_numpy(),r.loc[m.source_sample,col].fillna('missing').to_numpy())
p=np.load(OUT/'inputs/rare_panel_pseudobulk_doses.npz',allow_pickle=False);x=p['expected_panel_CPM'];m=pd.read_csv(OUT/'results/rare_panel_pseudobulk_dose_manifest.csv')
assert x.shape[0]==len(m)==105 and np.all(x>=0) and np.allclose(x.sum(1),1e6,rtol=1e-6)
assert m.groupby(['state','sample_id']).parent_fraction_fixed.nunique().eq(1).all()
assert m.groupby(['state','sample_id']).state_fraction_all_cells.min().eq(0).all()
for (_,sample),z in m.groupby(['state','sample_id']):
    z=z.sort_values('state_fraction_all_cells');lo=x[z.row.iloc[0]];hi=x[z.row.iloc[-1]];fmax=z.state_fraction_all_cells.iloc[-1]
    for row in z.itertuples():assert np.allclose(x[row.row],lo+(hi-lo)*(row.state_fraction_all_cells/fmax),atol=.02,rtol=1e-5)
c=pd.read_csv(OUT/'results/reference_metric_controls.csv');assert (c[c.control=='shuffle_vs_its_own_parent'].marginal_W1==0).all()
delta=pd.read_csv(OUT/'results/reference_metric_control_deltas.csv');q=delta[(delta.control=='independent_gene_shuffle')&(delta.metric.isin(['corr_MAE','covariance_MAE','energy_distance']))]
assert (q.ci_lo>0).all()
qa=dict(status='passed_preparation_checks',real_counts=checks,swaps=dict(assignments=int(ok.sum()),patient_disjoint=True,clinical_strata_match=True),rare_inputs=dict(n=105,panel_CPM_sum=1000000,parent_fraction_fixed=True,linearly_increasing_dose_verified=True,PBISC_response_tested=False),metric_calibration=dict(real_sample_type_groups=int(len(c)/4),gene_shuffle_marginals_exactly_preserved=True,joint_metrics_detect_disruption_both_cohorts=True,PBISC_tested=False),all_four_requested_stages_complete=False)
json.dump(qa,open(OUT/'VALIDATION.json','w'),indent=2)
print(json.dumps(qa,indent=2))
