"""Acquire real held-out count matrices locally; never substitute bulk-only predictions."""
import sys,os,json,hashlib,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1]
PREV=ROOT/'bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014'
sys.path.insert(0,str(PREV/'code'))
from common import read_csr_rows
import numpy as np,pandas as pd,h5py,anndata as ad
(OUT/'inputs').mkdir(exist_ok=True)
info=[]
for cohort,short,path in [
 ('ye2026_sepsis','sepsis','sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad'),
 ('imperial_eskd','covid','covid_eskd/eskd_covid19_recap.h5ad')]:
    src=ROOT/'bulkonly_lora5_20261006/analysis/original_scRNA'/path
    roster=pd.read_csv(OUT/f'results/{cohort}_strict_heldout_roster.csv')
    generated=ROOT/f'bulkonly_lora5_20261006/{cohort}/lora5_bulkonly_{cohort}_draw0.h5ad'
    with h5py.File(src,'r') as f: obs=ad.io.read_elem(f['obs']);rv=ad.io.read_elem(f['var'])
    with h5py.File(generated,'r') as f: gv=ad.io.read_elem(f['var'])
    if short=='sepsis':
        # Single explicit naming discrepancy: raw GEO/local recap Skin-152,
        # benchmark roster Skin-S152. Keep the original name in provenance.
        sid=obs['sample'].astype(str).replace({'Skin-152':'Skin-S152'});major=obs.celltype_major_ye2026.astype(str);fine=obs.celltype_ye2026.astype(str)
        keep=sid.isin(roster.sample_id)&obs.in_paper_scheme.astype(bool)
        raw_ids=rv.gene_ids.astype(str).str.split('.').str[0];mode='unique ENSG'
    else:
        sid=obs.sample_id.astype(str);major=obs.cell_type_detailed_recap.astype(str);fine=obs.cell_type_IFNstim_recap.astype(str)
        keep=sid.isin(roster.sample_id)
        raw_ids=rv.index.astype(str);mode='unique gene symbol'
    gkey=gv.gene_id.astype(str).str.split('.').str[0] if short=='sepsis' else gv.gene_symbol.astype(str)
    rawmap=pd.Series(np.arange(len(rv)),index=raw_ids)
    rawmap=rawmap[~rawmap.index.duplicated(keep=False)]
    good=gv.measured_in_cohort.astype(bool).to_numpy()&gkey.isin(rawmap.index).to_numpy()&~gkey.duplicated(keep=False).to_numpy()
    genes=gv.loc[good].copy();cols=rawmap.loc[gkey[good]].to_numpy();rows=np.where(keep)[0]
    samplemap=roster.set_index('sample_id')
    md=pd.DataFrame({'sample_id':sid.iloc[rows].values,'major':major.iloc[rows].values,'fine':fine.iloc[rows].values,'source_row':rows},index=obs.index[rows].astype(str))
    md['original_sample_id']=(obs['sample'] if short=='sepsis' else obs.sample_id).iloc[rows].astype(str).values
    md['patient_id']=md.sample_id.map(samplemap.donor_id);md['group']=md.sample_id.map(samplemap.group)
    if short=='covid': assert np.array_equal(md.patient_id.to_numpy(),obs.patient_id.iloc[rows].astype(str).to_numpy())
    assert set(md.sample_id)==set(roster.sample_id)
    for j in range(5): assert not roster[f'donor_trained_d{j}'].astype(str).str.lower().isin(['true','1']).any()
    dest=OUT/f'inputs/{cohort}_real_heldout_panel_counts.h5ad'
    if not dest.exists():
        print('Reading',cohort,len(rows),'cells',len(cols),'genes',flush=True)
        x=read_csr_rows(src,rows,cols)
        assert np.all(x.data>=0) and np.all(x.data==np.round(x.data))
        md['panel_counts']=np.asarray(x.sum(1)).ravel().astype(np.int64)
        a=ad.AnnData(x,obs=md,var=genes)
        a.uns['acquisition']={'original_path':str(src),'gene_matching':mode,'scope':'all five adapters held-out, measured unique panel genes','labels':'local recap, not original benchmark scANVI','full_transcriptome_in_original_file':str(src),'generated_paired_counts_available':False}
        a.write_h5ad(dest,compression='lzf');del a,x
    else:
        with h5py.File(dest,'r') as f: md=ad.io.read_elem(f['obs'])
    print('Saved',cohort,len(md),'cells',md.patient_id.nunique(),'patients',flush=True)
    md.groupby(['patient_id','sample_id','group','major','fine'],observed=True).size().rename('n').reset_index().to_csv(OUT/f'results/{cohort}_real_cell_census.csv',index=False)
    h=hashlib.sha256()
    with open(dest,'rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    info.append(dict(cohort=cohort,path=str(dest),sha256=h.hexdigest(),n_cells=len(md),n_genes=len(genes),samples=md.sample_id.nunique(),patients=md.patient_id.nunique(),gene_matching=mode,original_source=str(src),counts_integer=True))
json.dump(info,open(OUT/'results/real_counts_manifest.json','w'),indent=2)
