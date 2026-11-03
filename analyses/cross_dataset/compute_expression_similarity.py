"""Descriptive centroid similarity: original versus LoRA cells within broad labels.
Uses matched expressed genes and log-normalized counts. Cell centroids are descriptive;
the statistical replicate is not the cell and no cell-level p-values are computed.
"""
from pathlib import Path
import h5py, anndata as ad, numpy as np, pandas as pd
from scipy import sparse, stats

ROOT=Path(r"F:\PBISC")
ORIG_BASE=ROOT/"bulkonly_lora5_20261006"
LORA_BASE=ROOT/"bulkonly_lora5_k4000_20261006"
OUT=ROOT/"paper_claim_analysis"
OUT.mkdir(exist_ok=True)
PY=ROOT/".analysis_runtime"

def decode(x): return x.decode() if isinstance(x,bytes) else str(x)
def read_meta(path):
    with h5py.File(path,'r') as f:
        obs=ad.io.read_elem(f['obs'])
        var=ad.io.read_elem(f['var'])
        if 'gene_symbol' in var: genes=var['gene_symbol'].astype(str).to_numpy()
        else: genes=var.index.astype(str).to_numpy()
        if 'measured_in_cohort' in var: measured=var['measured_in_cohort'].astype(bool).to_numpy()
        else: measured=np.ones(len(genes),bool)
    return obs,genes,measured

def aggregate(path, cohort, dataset, draw=None):
    obs,genes,measured=read_meta(path)
    if cohort=='sepsis':
        if dataset=='original':
            labels=obs['celltype_major_ye2026'].astype(str).to_numpy()
            cond=obs['group'].astype(str).to_numpy()
            age=obs['age_group'].astype(str).to_numpy()
            cond=np.where((age=='adult') & (cond=='sepsis'),'adult_sepsis',cond)
            cond=np.where((age=='adult') & pd.Series(cond).str.contains('healthy',case=False).to_numpy(),'adult_HC',cond)
            cellmaps={'Monocytes':'Monocyte','B cells':'B cell','NK cells':'NK cell','DC':'DC','Plasma cells':'B cell'}
            cell=np.array([cellmaps.get(x,x) for x in labels])
            keep=(age=='adult') & np.isin(cond,['adult_sepsis','adult_HC']) & np.isin(cell,['Monocyte','B cell','NK cell','DC'])
        else:
            cell=obs['celltype_l1_pred'].astype(str).to_numpy(); cond=obs['group'].astype(str).to_numpy()
            keep=np.isin(cell,['Monocyte','B cell','NK cell','DC']) & np.isin(cond,['adult_sepsis','adult_HC'])
    else:
        if dataset=='original':
            labels=obs['cell_type_detailed_recap'].astype(str).to_numpy()
            cond=obs['case_control'].astype(str).str.lower().to_numpy()
            cmap={'CD14mono':'Monocyte','CD16mono':'Monocyte','Int.mono':'Monocyte','NK':'NK cell','NKT':'NK cell','T CD4 Helper':'T cell','T CD4 Naive':'T cell','T CD4 CTL':'T cell','T CD8 CTL':'T cell','T CD8 Naive':'T cell','T CD8 CM':'T cell','B_naive':'B cell','B_switched_mem':'B cell','B_non-switched_mem':'B cell','B_ASC_IgA':'B cell','B_ASC_IgG':'B cell','B_ASC_IgM':'B cell','B_ASC_dividing':'B cell','DC1':'DC','DC2':'DC','DC3':'DC','pDC':'DC'}
            cell=np.array([cmap.get(x,x) for x in labels]); keep=(cond=='positive') & np.isin(cell,['Monocyte','NK cell','T cell','B cell','DC'])
            cond=np.where(cond=='positive','covid_positive',cond)
        else:
            cell=obs['celltype_l1_pred'].astype(str).to_numpy(); cond=obs['group'].astype(str).to_numpy()
            keep=(cond=='covid_positive') & np.isin(cell,['Monocyte','NK cell','T cell','B cell','DC'])
    ix=np.flatnonzero(measured)
    syms=genes[ix]
    # Collapse duplicate symbols by retaining first occurrence to avoid double weighting.
    _,first=np.unique(syms,return_index=True); ix=ix[np.sort(first)]; syms=genes[ix]
    index={g:i for i,g in enumerate(syms)}
    types=sorted(set(cell[keep])); sums={t:np.zeros(len(ix),np.float64) for t in types}; n={t:0 for t in types}
    with h5py.File(path,'r') as f:
        x=f['X']; ip=x['indptr']; ids=x['indices']; dat=x['data']; nrows=len(obs)
        remap=np.full(len(genes),-1,np.int32); remap[ix]=np.arange(len(ix),dtype=np.int32)
        # Process contiguous row blocks; only count matrix is read, never dense logcounts.
        for lo in range(0,nrows,10000):
            hi=min(lo+10000,nrows); rows=np.arange(lo,hi); mask=keep[lo:hi]
            if not mask.any(): continue
            a,b=int(ip[lo]),int(ip[hi]); lens=np.diff(ip[lo:hi+1]); rr=np.repeat(np.arange(hi-lo),lens)
            cc=ids[a:b]; vv=dat[a:b]; mapped=remap[cc]; ok=(mapped>=0)&mask[rr]
            mat=sparse.csr_matrix((vv[ok].astype(np.float64),(rr[ok],mapped[ok])),shape=(hi-lo,len(ix)))
            lib=np.asarray(mat.sum(axis=1)).ravel(); scale=np.divide(1e4,lib,out=np.zeros_like(lib),where=lib>0)
            mat=sparse.diags(scale)@mat
            mat.data=np.log1p(mat.data)
            labs=cell[lo:hi]
            for t in types:
                sel=np.flatnonzero(mask & (labs==t))
                if len(sel):
                    sums[t]+=np.asarray(mat[sel].sum(axis=0)).ravel(); n[t]+=len(sel)
    return {t:sums[t]/max(1,n[t]) for t in types},n,syms

def main():
    rows=[]
    cohorts={'sepsis':('analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad','ye2026_sepsis/lora5k4_bulkonly_ye2026_sepsis_draw{}.h5ad'),
             'covid_eskd':('analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad','imperial_eskd/lora5k4_bulkonly_imperial_eskd_draw{}.h5ad')}
    for cohort,(origrel,genpat) in cohorts.items():
        orig,ocount,genes=aggregate(ORIG_BASE/origrel,cohort,'original')
        for d in range(5):
            gen,gcount,ggenes=aggregate(LORA_BASE/genpat.format(d),cohort,'LoRA',d)
            lookup={g:i for i,g in enumerate(ggenes)}; common=[(i,lookup[g]) for i,g in enumerate(genes) if g in lookup]
            oi,gi=np.array(common).T
            for t in sorted(set(orig)&set(gen)):
                x=orig[t][oi]; y=gen[t][gi]
                pear=float(stats.pearsonr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else np.nan
                spear=float(stats.spearmanr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else np.nan
                rows.append(dict(cohort=cohort,cell_type=t,draw=d,n_genes=len(common),original_cells=ocount[t],generated_cells=gcount[t],pearson_r=pear,spearman_r=spear))
    pd.DataFrame(rows).to_csv(OUT/'expression_centroid_similarity.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__': main()
