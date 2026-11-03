"""Prepare bulk-swap assignments and panel pseudobulk interventions, without inference."""
import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1]
PREV=ROOT/'bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014'
sys.path.insert(0,str(PREV/'code'))
from common import lin6
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad
rng=np.random.default_rng(20261014);swaps=[]
for cohort in ['ye2026_sepsis','imperial_eskd']:
    m=pd.read_csv(OUT/f'results/{cohort}_strict_heldout_roster.csv')
    strata=['group'] if cohort=='ye2026_sepsis' else ['group','wave','severity_peak']
    for _,row in m.iterrows():
        z=m[m.donor_id!=row.donor_id].copy()
        for col in strata:z=z[z[col].fillna('missing').astype(str)==str(row[col] if pd.notna(row[col]) else 'missing')]
        for rep in range(3):
            if len(z):
                pid=rng.choice(sorted(z.donor_id.unique()));src=z[z.donor_id==pid].iloc[rng.integers((z.donor_id==pid).sum())]
                swaps.append(dict(cohort=cohort,target_sample=row.sample_id,target_patient=row.donor_id,target_bulk=row.get('bulk_sample',row.sample_id),source_sample=src.sample_id,source_patient=src.donor_id,source_bulk=src.get('bulk_sample',src.sample_id),replicate=rep,eligible_donors=z.donor_id.nunique(),status='assignment_ready_inference_not_run',stratum='|'.join(str(row[c]) for c in strata)))
            else:swaps.append(dict(cohort=cohort,target_sample=row.sample_id,target_patient=row.donor_id,replicate=rep,eligible_donors=0,status='no_other_heldout_patient_in_stratum',stratum='|'.join(str(row[c]) for c in strata)))
s=pd.DataFrame(swaps);s.to_csv(OUT/'results/bulk_swap_assignments.csv',index=False)
ok=s.status.eq('assignment_ready_inference_not_run');assert (s.loc[ok,'target_patient']!=s.loc[ok,'source_patient']).all()
a=ad.read_h5ad(OUT/'inputs/imperial_eskd_real_heldout_panel_counts.h5ad')
md=a.obs.reset_index(drop=True);md['parent']=md.major.map(lin6)
x=a.X.tocsr();depth=np.asarray(x.sum(1)).ravel();profiles=(sp.diags(1/np.maximum(depth,1))@x).tocsr()
mix=[];des=[];states=['MAIT','NK IFN stim','CD14mono_IFN'];doses=[0,.0005,.001,.002,.005,.01,.02]
for state in states:
    eligible=md[md.fine==state].sample_id.value_counts();eligible=eligible[eligible>=20]
    seen=set();selected=[]
    for sample in eligible.index:
        pid=md.loc[md.sample_id==sample,'patient_id'].iloc[0]
        if pid not in seen:selected.append(sample);seen.add(pid)
        if len(selected)==5:break
    for sample in selected:
        sm=md.sample_id.eq(sample).to_numpy()&(depth>0);st=sm&md.fine.eq(state).to_numpy()
        parents=md.loc[st,'parent'].unique();assert len(parents)==1;parent=parents[0]
        bg=sm&md.parent.eq(parent).to_numpy()&~st;other=sm&~md.parent.eq(parent).to_numpy()
        if bg.sum()<100 or not other.any():continue
        parent_fraction=(bg.sum()+st.sum())/sm.sum()
        state_profile=np.asarray(profiles[st].mean(0)).ravel();bg_profile=np.asarray(profiles[bg].mean(0)).ravel();outside=np.asarray(profiles[other].mean(0)).ravel()
        for dose in doses:
            if dose>=parent_fraction:continue
            v=1e6*((1-parent_fraction)*outside+(parent_fraction-dose)*bg_profile+dose*state_profile)
            assert np.isclose(v.sum(),1e6,rtol=1e-6) and (v>=0).all()
            mix.append(v.astype(np.float32));des.append(dict(row=len(mix)-1,state=state,sample_id=sample,patient_id=md.loc[sm,'patient_id'].iloc[0],state_fraction_all_cells=dose,parent=parent,parent_fraction_fixed=parent_fraction,n_real_state=int(st.sum()),n_real_parent_background=int(bg.sum()),status='panel_pseudobulk_input_only_PBISC_not_run'))
np.savez_compressed(OUT/'inputs/rare_panel_pseudobulk_doses.npz',expected_panel_CPM=np.vstack(mix),gene_ids=a.var.gene_id.to_numpy().astype(str),gene_symbols=a.var.gene_symbol.to_numpy().astype(str))
pd.DataFrame(des).to_csv(OUT/'results/rare_panel_pseudobulk_dose_manifest.csv',index=False)
json.dump(dict(scope='preparation only; not a removal/spike-in response result',condition_warning='5745-gene measured panel CPM is not the canonical 20097-gene PBISC condition. Reconstruct full-gene pseudobulks and apply exact server condition-builder before inference.',depth='per-cell panel library normalized; equal RNA contribution per cell',composition='total parent lineage fraction fixed; rare state replaces state-negative cells within parent',state_sources='up to 5 distinct strict-heldout patients per state, >=20 real state cells',seeds=[2026101401,2026101402,2026101403],matched_N_primary=1500,rare_N=30000,swap_generation='use same checkpoint and noise/readout seeds for native and swapped condition; no retraining',independent_seeds='three independent seeds per unchanged checkpoint; keep adapter draw and readout seed effects separate'),open(OUT/'results/experiment_input_methods.json','w'),indent=2)
print('Swap eligible samples',s.loc[ok,'target_sample'].nunique(),'of',s.target_sample.nunique(),flush=True)
print('Intervention inputs',len(des),'rows',len(set(d['patient_id'] for d in des)),'patients',flush=True)
