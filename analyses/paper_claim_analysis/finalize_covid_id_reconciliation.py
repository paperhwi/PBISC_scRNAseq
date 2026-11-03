from pathlib import Path
import h5py, numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'paper_claim_analysis'
LOCAL=OUT/'covid_sample_patient_manifest_local.csv'
PUB=OUT/'Stephenson_Table_S1_sample_metadata.xlsx'
local=pd.read_csv(LOCAL,dtype=str).fillna('')
pub=pd.read_excel(PUB,header=2,dtype={'sample_id':str,'patient_id':str}).fillna('')
def key(s):
    prefix, sep, suffix=str(s).strip().partition('_')
    return f'{prefix}_{suffix[:4]}' if sep and len(suffix)>=4 else str(s).strip()
local['sample_key']=local.sample_id.map(key)
pub['sample_key']=pub.sample_id.map(key)
if local.sample_key.duplicated().any() or pub.sample_key.duplicated().any():
    raise ValueError('Sample ID normalization is not one-to-one; stop rather than silently merge.')
pub_cols=pub[['sample_key','sample_id','patient_id','case_control','year_of_sample','WHO_severity','fatal_disease','Received_Steroid']].rename(columns={
    'sample_id':'published_sample_id','patient_id':'published_patient_id','case_control':'published_case_control','year_of_sample':'published_year','WHO_severity':'published_WHO_severity','fatal_disease':'published_fatal_disease','Received_Steroid':'published_received_steroid'})
m=local.merge(pub_cols,on='sample_key',how='left',validate='one_to_one')
m['in_published_187_sample_table']=m.published_sample_id.fillna('').ne('')
m['patient_id_matches']=m.patient_id.eq(m.published_patient_id)
m['case_control_matches']=m.case_control.eq(m.published_case_control)
m['reconciliation_status']=np.where(m.in_published_187_sample_table,'included in published analysis sample table',np.where(m.patient_id.isin(set(pub.patient_id)),'sample absent from published analysis sample table','patient absent from published analysis sample table'))
m['exclusion_reason']=np.where(m.in_published_187_sample_table,'',np.where(m.patient_id.isin(set(pub.patient_id)),'sample ID absent from final Table S1 roster; sample-level exclusion reason is not specified in Table S1','patient absent from Table S1; two excluded patients are consistent with the article-stated benign CLL removals'))
m.to_csv(OUT/'covid_sample_patient_id_reconciliation.csv',index=False)
# Strong consistency checks for all retained IDs.
ret=m[m.in_published_187_sample_table]
assert len(ret)==187 and ret.patient_id.nunique()==61
assert ret.patient_id_matches.all(), 'patient ID mismatch among published samples'
assert ret.case_control_matches.all(), 'case-control mismatch among published samples'
extra=m[~m.in_published_187_sample_table]
assert len(extra)==11 and extra.patient_id.nunique()==3
assert set(extra.patient_id)=={'C116','C139','C141'}
assert not (m[~m.in_published_187_sample_table].published_sample_id.fillna('').ne('')).any()
# Correct exact cell count from local H5AD's per-sample cell tally.
cells_in=int(ret.n_cells.astype(int).sum()); cells_all=int(local.n_cells.astype(int).sum()); cells_extra=int(extra.n_cells.astype(int).sum())
assert cells_all==588389 and cells_in==578835
assert cells_extra==9554
# Preserve all-object per-sample files while making canonical filenames paper-aligned.
for name in ['original_covid_ifn_frequencies_by_sample_and_parent.csv','original_covid_dex_mono_frequencies_by_sample.csv','original_covid_ifn_marker_scores_by_sample.csv','original_covid_dex_marker_scores_by_sample.csv','original_covid_claim_frequencies_by_sample.csv']:
    p=OUT/name
    if not p.exists(): continue
    d=pd.read_csv(p)
    if 'sample_id' not in d: continue
    bak=OUT/(p.stem+'_all_local_198.csv')
    if not bak.exists(): d.to_csv(bak,index=False)
    d=d[d.sample_id.astype(str).map(key).isin(set(ret.sample_key))].copy()
    d.to_csv(p,index=False)
# Recompute descriptive summaries from publication-cohort subsets.
ifn=pd.read_csv(OUT/'original_covid_ifn_frequencies_by_sample_and_parent.csv')
ifn.groupby(['population','case_control','WHO_severity','centre'],dropna=False).agg(
    n_samples=('sample_id','nunique'),median_ifn_fraction_parent=('ifn_fraction_parent','median'),
    mean_ifn_fraction_parent=('ifn_fraction_parent','mean'),samples_ifn_positive=('ifn_cells',lambda x:int((x>0).sum()))
).reset_index().to_csv(OUT/'original_covid_claim_summary.csv',index=False)
dex=pd.read_csv(OUT/'original_covid_dex_mono_frequencies_by_sample.csv')
dex.groupby(['case_control','WHO_severity','centre'],dropna=False).agg(
    n_samples=('sample_id','nunique'),n_patients=('patient_id','nunique'),samples_with_candidate_dex=('dex_cells',lambda x:int((x>0).sum())),
    median_candidate_dex_fraction=('dex_fraction_monocytes','median'),mean_candidate_dex_fraction=('dex_fraction_monocytes','mean')
).reset_index().to_csv(OUT/'original_covid_dex_mono_summary.csv',index=False)
# Update score summaries using matched sample rows.
for stem, value, title in [('original_covid_ifn_marker_scores_by_sample','IFN_score_difference','IFN'),('original_covid_dex_marker_scores_by_sample','paired_score_difference','DEX')]:
    p=OUT/(stem+'.csv')
    if not p.exists(): continue
    d=pd.read_csv(p)
    if title=='IFN':
        rows=[]
        for parent,g in d.groupby('parent_cell_type'):
            z=g[value].dropna()
            if len(z):
                from scipy.stats import wilcoxon
                rows.append({'parent_cell_type':parent,'n_sample_parent_pairs':len(z),'median_IFN_score_difference':z.median(),'positive_differences':int((z>0).sum()),'paired_wilcoxon_pvalue':wilcoxon(z).pvalue if np.any(z!=0) else np.nan})
        pd.DataFrame(rows).to_csv(OUT/'original_covid_ifn_marker_score_summary.csv',index=False)
    else:
        from scipy.stats import wilcoxon
        if 'DEX_up_score_candidate' in d.columns:
            rows=[]
            for panel, a,b in [('DEX_up','DEX_up_score_candidate','DEX_up_score_background'),('DEX_down','DEX_down_score_candidate','DEX_down_score_background')]:
                z=(d[a]-d[b]).dropna()
                rows.append({'panel':panel,'n_sample_pairs':len(z),'median_candidate_score':d[a].median(),'median_background_score':d[b].median(),'median_paired_difference':z.median(),'positive_differences':int((z>0).sum()),'paired_wilcoxon_pvalue':wilcoxon(z).pvalue if np.any(z!=0) else np.nan})
            pd.DataFrame(rows).to_csv(OUT/'original_covid_dex_marker_score_summary.csv',index=False)
# manifest of published sample IDs is the frozen rerun input.
pub[['sample_id','patient_id','case_control','year_of_sample','WHO_severity','fatal_disease','Received_Steroid']].to_csv(OUT/'covid_published_187_sample_ids.csv',index=False)
audit=pd.DataFrame([
 {'cohort':'local source object after exclude==0','cells':cells_all,'samples':len(local),'patients':local.patient_id.nunique(),'status':'full downloaded object; includes 2 patients not analyzed in article'},
 {'cohort':'published article analysis cohort, matched to Supplementary Table S1','cells':cells_in,'samples':len(ret),'patients':ret.patient_id.nunique(),'status':'187/61 sample-patient roster matches Table S1; local cells remain 1,205 below paper count'},
 {'cohort':'removed from published analysis set','cells':cells_extra,'samples':len(extra),'patients':extra.patient_id.nunique(),'status':'11 sample IDs / 9,554 cells absent: all local samples for C116 and C139, plus C141_0903'},
])
audit.to_csv(OUT/'original_covid_cohort_size_audit.csv',index=False)
summary=f'''# COVID-19/ESKD published cohort ID reconciliation\n\nDate: 2026-10-12\n\n## Result\n\nThe deposited local H5AD metadata contain 198 biological `sample_id` values and 63 `patient_id` values (588,389 cells). Sample IDs were reconciled against the article’s Supplementary Table S1, which contains 187 sample IDs and 61 patient IDs. Where the publisher’s ID appends the two-digit year to the month/day suffix, IDs were normalized to the shared patient + first four date digits; this was one-to-one in both files. Patient IDs and case-control labels matched for every retained sample.\n\nExactly 11 local sample IDs from two patients (`C116`, `C139`) are absent from the published sample table. Their local samples contain {cells_extra:,} cells. Removing the 11 unmatched sample IDs yields 578,835 cells, 187 samples and 61 patients. The sample and patient roster now matches Table S1 exactly, but the local selected cells are 1,205 fewer than the article's 580,040. Sample identity is reconciled; the remaining 1,205-cell difference is a cell-level QC/provenance mismatch and is not assigned to samples without the authors’ cell inclusion manifest. Two absent patients are C116 and C139; the article methods state that two patients with benign chronic lymphocytic leukemia were removed, so this patient-level mapping is consistent with that stated exclusion (the article does not name them in Table S1). The remaining absent sample is C141_0903, while C141 has other retained samples. Its specific exclusion reason is not enumerated in Table S1; the paper states its reported roster is after sample-level QC exclusions.\n\n## Application\n\nCOVID sample-level IFN-state, dexamethasone-candidate frequency, and marker-score tables now use only the 187 IDs in the published Table S1 roster. Previous all-local summaries are retained under filenames ending `_all_local_198.csv`. The cell-level metadata source remains an unfiltered copy of the deposited object; downstream aggregation should use `covid_published_187_sample_ids.csv` or the reconciled manifest.\n\n## Files\n\n- `covid_sample_patient_id_reconciliation.csv`: all 198 local sample IDs, matched article IDs, patient/status checks, cell counts and inclusion flag.\n- `covid_published_187_sample_ids.csv`: frozen article sample roster for reruns.\n- `Stephenson_Table_S1_sample_metadata.xlsx`: publisher’s full sample metadata workbook.\n- `original_covid_cohort_size_audit.csv`: local, matched and removed cohort counts.\n\n## Reproduction caveat\n\nThis fixes sample/patient identity and the cohort denominator. It does not restore authors’ discarded cells, redo doublet/quality filtering, or rerun their original differential abundance and mixed-effects models.\n'''
(OUT/'COVID_COHORT_ID_RECONCILIATION.md').write_text(summary,encoding='utf-8')
print('included',len(ret),ret.patient_id.nunique(),cells_in,'removed',len(extra),extra.patient_id.nunique(),cells_extra)
print(extra[['sample_id','patient_id','centre','case_control','sample_date','n_cells']].to_string(index=False))
print('consistency matches: sample-level patient IDs and case_control all exact')
