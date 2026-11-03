"""Read-only acquisition audit. Counts are distinguished from feature summaries."""
import os, sys, json, tarfile, zipfile, socket
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parents[1]
PREV=ROOT/'bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014'
sys.path[:0]=[str(PREV/'runtime'),str(ROOT/'.analysis_runtime')]
import numpy as np
import pandas as pd
import h5py, anndata as ad
(OUT/'results').mkdir(parents=True,exist_ok=True)
archives=[]
for branch in ['bulkonly_lora5_20261006','bulkonly_lora5_k4000_20261006']:
    path=ROOT/branch/'analysis/analysis_bundle.tar.gz'
    with tarfile.open(path,'r|gz') as t:
        for m in t:
            if m.isfile(): archives.append(dict(archive=str(path),member=m.name,bytes=m.size))
pd.DataFrame(archives).to_csv(OUT/'results/archive_inventory.csv',index=False)
possible=[r for r in archives if r['member'].lower().endswith(('.npz','.npy','.h5ad','.pt','.pth','.ckpt','.safetensors'))]
pd.DataFrame(possible).to_csv(OUT/'results/archive_matrix_candidates.csv',index=False)
local=[]
for base in [ROOT/'bulkonly_lora5_20261006',ROOT/'bulkonly_lora5_k4000_20261006',ROOT/'bulkonly_lora5_20261007']:
    for p in base.rglob('*'):
        if not p.is_file() or 'runtime' in p.parts: continue
        if p.suffix.lower() in {'.npz','.h5ad','.pt','.pth','.ckpt','.safetensors'}:
            local.append(dict(path=str(p),bytes=p.stat().st_size,kind=p.suffix))
pd.DataFrame(local).to_csv(OUT/'results/local_matrix_inventory.csv',index=False)
for cohort in ['ye2026_sepsis','imperial_eskd']:
    meta=pd.read_csv(PREV/f'inputs/k1000/analysis/features/{cohort}/meta.csv')
    hold=meta.set.eq('paired')
    for j in range(5): hold &= meta[f'paired_test_d{j}'].astype(str).str.lower().isin(['true','1'])
    roster=meta.loc[hold].copy()
    roster.to_csv(OUT/f'results/{cohort}_strict_heldout_roster.csv',index=False)
    print(cohort,'strict heldout',len(roster),'samples',roster.donor_id.nunique(),'patients')
paths={
 'sepsis':ROOT/'bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad',
 'covid':ROOT/'bulkonly_lora5_20261006/analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad',
 'generated_sepsis':ROOT/'bulkonly_lora5_20261006/ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw0.h5ad',
 'generated_covid':ROOT/'bulkonly_lora5_20261006/imperial_eskd/lora5_bulkonly_imperial_eskd_draw0.h5ad'}
schemas={}
for tag,path in paths.items():
    with h5py.File(path,'r') as f:
        obs=ad.io.read_elem(f['obs']);var=ad.io.read_elem(f['var'])
        schemas[tag]=dict(path=str(path),shape=[int(x) for x in f['X'].attrs.get('shape',f['X'].shape if isinstance(f['X'],h5py.Dataset) else [])],obs_columns=list(obs),var_columns=list(var),layers=list(f.get('layers',{})),obs_examples={c:obs[c].astype(str).unique()[:8].tolist() for c in obs if any(s in c.lower() for s in ['sample','patient','donor','recap'])},var_head=var.head(3).astype(str).reset_index().to_dict('records'))
        print(tag,json.dumps(schemas[tag],ensure_ascii=False))
json.dump(schemas,open(OUT/'results/input_schemas.json','w'),indent=2,ensure_ascii=False)
try:
    remote_address=socket.gethostbyname('spark-0110');dns_error=None
except OSError as e:
    remote_address=None;dns_error=str(e)
summary=dict(archive_members=len(archives),matrix_candidates=possible,local_checkpoints=[r for r in local if r['kind'] in {'.pt','.pth','.ckpt','.safetensors'}],remote_host='spark-0110',dns_resolved=remote_address is not None,dns_error=dns_error,ssh_config_present=(Path.home()/'.ssh/config').exists(),paired_generated_full_counts_acquired=False)
json.dump(summary,open(OUT/'results/acquisition_audit.json','w'),indent=2)
print('Archive matrix candidates:',len(possible))
for r in possible: print(r['member'],r['bytes'])
