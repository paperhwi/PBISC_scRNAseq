"""Independent invariants and benchmark reconciliation; fails on leakage."""
from common import *
import json,hashlib,importlib.metadata
checks=[]
s=pd.read_csv(OUT/'results/paired_accuracy_summary.csv')
for coh in ['ye2026_sepsis','imperial_eskd']:
    for draw in range(5):
        raw=pd.read_csv(OUT/'inputs/k4000/scores'/coh/f'draw{draw}/per_donor.csv')
        for K,arm in [('k1000','lora5'),('k4000','lora5_k4000')]:
            a=raw[raw.arm==arm].comp_l1.mean()
            b=s[(s.cohort==coh)&(s.scope=='all_paired')&(s.level=='l1')&(s.method==f'PBISC_{K}')&(s.draw==draw)].iloc[0].sample_weighted_L1
            assert abs(a-b)<1e-10,(coh,K,draw,a,b)
checks.append('All 20 original per-draw paired benchmark means reproduced to 1e-10')
pi=json.load(open(OUT/'results/paired_integrity.json'))
assert all(x['train_test_patient_overlap']==0 and x['train_donors']==5 for x in pi)
checks.append('Every paired adapter excludes all visits of its five training donors')
mm=json.load(open(OUT/'results/covid_mapping_methods.json'))
patient_folds=[x for x in mm['split_rosters'] if x['protocol']=='patient_disjoint']
assert len(patient_folds)==2
for f in patient_folds:assert not set(f['ref_patients'])&set(f['test_patients']) and not f['overlap_patients']
test0,test1=[set(f['test_patients']) for f in patient_folds]
assert not test0&test1 and len(test0|test1)==61
checks.append('COVID 2-fold OOF: all 61 patients tested exactly once; zero shared train/test patients')
assert not any(f['warnings'] for f in mm['classifier_fit_log']),mm['classifier_fit_log']
checks.append('Every format-aware logistic classifier converged')
for f in patient_folds:
    d=pd.read_csv(OUT/f"results/patient_disjoint_fold{f['fold']}_generated_labels.csv.gz",usecols=['patient'])
    assert not set(d.patient)&set(f['ref_patients'])
checks.append('Generated visits never mapped by a classifier trained on that real-reference patient')
m=pd.read_csv(OUT/'results/covid_mapping_patient_summary.csv')
assert m.patient_accuracy.between(0,1).all() and m.composition_L1.between(0,2).all()
assert (m.accuracy_lo<=m.patient_accuracy).all() and (m.patient_accuracy<=m.accuracy_hi).all()
checks.append('Mapping accuracy, composition bounds and bootstrap intervals are valid')
effects=pd.read_csv(OUT/'results/covid_patient_effect_summary.csv')
assert effects.valid_bootstraps.min()>=1900
checks.append('Patient effect bootstrap has at least 1,900 valid draws for every reported comparison')
files=[]
for folder in ['code','results','figs']:
    for p in sorted((OUT/folder).glob('*')):
        if p.is_file():files.append(dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
versions={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','scikit-learn','h5py','anndata','scanpy','pyarrow','harmonypy','matplotlib']}
json.dump({'checks':checks,'versions':versions,'python':sys.version,'files':files},open(OUT/'VALIDATION.json','w'),indent=2)
print('\n'.join(checks))
