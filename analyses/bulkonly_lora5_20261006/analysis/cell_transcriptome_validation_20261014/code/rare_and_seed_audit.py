"""Rare-state eligibility, finite-cell power and fixed-checkpoint seed audit."""
import sys,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1]
PREV=ROOT/'bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014'
sys.path[:0]=[str(PREV/'runtime'),str(ROOT/'.analysis_runtime')]
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy.stats import binom
recover=pd.read_csv(PREV/'results/covid_mapping_recoverability.csv')
recover=recover[(recover.protocol=='patient_disjoint')&(recover.source=='realT2')&(recover.route=='R3')&(recover.level=='fine')]
census=pd.read_csv(OUT/'results/imperial_eskd_real_cell_census.csv');total=census.n.sum()
freq=census.groupby('fine').n.sum().rename('n_heldout_real').reset_index().rename(columns={'fine':'label'})
c= recover.merge(freq,on='label',how='left').fillna({'n_heldout_real':0})
c['heldout_frequency']=c.n_heldout_real/total
c['candidate_for_intervention']=(c.recall>=.5)&(c.precision>=.5)&(c.heldout_frequency<.02)&(c.n_true_cells>=25)
c['eligibility_scope']='preliminary realT2 mapper recoverability; independent markers and held-out precision still required'
c.to_csv(OUT/'results/rare_state_eligibility.csv',index=False)
power=[]
for p in [.0001,.0005,.001,.002,.005,.01,.02]:
    for n in [200,500,1500,4000,10000,30000]:
        power.append(dict(true_frequency=p,n_generated=n,expected_cells=n*p,p_zero=(1-p)**n,p_at_least_10=binom.sf(9,n,p),n_for_95pct_at_least_one=math.ceil(math.log(.05)/math.log1p(-p))))
pd.DataFrame(power).to_csv(OUT/'results/rare_sampling_power.csv',index=False)
provs=[];requests=[]
for K,branch,prefix in [('k1000','bulkonly_lora5_20261006','lora5'),('k4000','bulkonly_lora5_k4000_20261006','lora5k4')]:
    for cohort in ['ye2026_sepsis','imperial_eskd']:
        for j in range(5):
            path=ROOT/branch/cohort/f'{prefix}_bulkonly_{cohort}_draw{j}.h5ad'
            with h5py.File(path,'r') as f:u=ad.io.read_elem(f['uns'])
            p=json.loads(u['provenance'])
            provs.append(dict(K=K,cohort=cohort,draw=j,checkpoint_sha256=p.get('checkpoint_sha256'),lora_run=p.get('lora_run'),generation=p.get('generation'),poisson_seed=p.get('poisson_readout_seed'),training_donors=p.get('training_donors')))
            requests.append(dict(K=K,cohort=cohort,draw=j,checkpoint_dir=p.get('lora_run'),checkpoint_sha256=p.get('checkpoint_sha256'),compact_paired_dir=p.get('plan',f'/home/sj_server_1/PBISC/immun-pbmc-fm/model_8_sparse_data_diffusion/plans/{branch}')+f'/work/{cohort}/paired_gen/draw{j}',required=['full per-cell count or expected-count matrix with integer read-out and seed','row-aligned sample and patient IDs','gene IDs, symbols, measured-gene mask','paired bulk conditions and row mapping','test/train donor manifest','receipt and RUN.json'],insufficient_alone=['composition.parquet','markers.npz','pb_l1.npz','pb_l2.npz','200-cell subsample.npz']))
json.dump(provs,open(OUT/'results/generation_provenance.json','w'),indent=2)
json.dump(dict(remote_host='spark-0110 / server_1',pairs=requests,shared_required=['sample_sdd_standalone.py and its model dependencies','canonical_genes.tsv','exact condition-builder and sampler CLI from plan tools','base-model pretraining donor/sample roster and all model-selection cohorts','frozen annotation model or validated mapper','independent external cohort with bulk+scRNA and disjoint donors'],stages={'1':'real counts acquired, paired generated counts missing','2':'no PBISC coexpression or bulk swap executed','3':'no PBISC removal/spike-in executed','4':'no fixed-checkpoint reseeding or external-cohort inference executed'}),open(OUT/'results/server_acquisition_request.json','w'),indent=2)
summary=[]
for (K,cohort),grp in pd.DataFrame(provs).groupby(['K','cohort']):
    summary.append(dict(K=K,cohort=cohort,models=len(grp),unique_checkpoints=grp.checkpoint_sha256.nunique(),fixed_checkpoint_replicate_seeds=0,independent_seed_experiment_complete=False))
pd.DataFrame(summary).to_csv(OUT/'results/fixed_checkpoint_seed_audit.csv',index=False)
print('Rare preliminary candidates:',c.loc[c.candidate_for_intervention,['label','heldout_frequency','recall','precision']].to_string(index=False))
print(pd.DataFrame(summary).to_string(index=False))
