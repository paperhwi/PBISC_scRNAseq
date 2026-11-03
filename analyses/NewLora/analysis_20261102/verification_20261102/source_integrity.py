"""Check exact New LoRA input selection and figure point provenance."""
from pathlib import Path
import json
import h5py
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=HERE.parent
rows=[]

def col(f,path):
    g=f[path]
    if isinstance(g,h5py.Group):
        cats=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in g['categories'][:]],dtype=object)
        return cats[g['codes'][:]]
    a=g[:]
    return np.array([x.decode() if isinstance(x,bytes) else x for x in a],dtype=object) if a.dtype.kind in 'SO' else a

for cohort,name in [('sepsis','ye2026_sepsis_draw1_K4000_ceiling.h5ad'),('covid','imperial_eskd_draw0_K4000_ceiling.h5ad')]:
    labels=pd.read_csv(BASE/cohort/'ceiling_cell_labels.csv.gz',index_col=0)
    summary=json.loads((BASE/cohort/'summary.json').read_text())
    with h5py.File(ROOT/'NewLora'/name) as f:
        sample=col(f,'obs/sample_id').astype(str)
        group=col(f,'obs/group').astype(str)
        eligible_flag=col(f,'obs/primary_analysis_eligible').astype(bool)
        keep=(~np.char.startswith(group,'ped')) if cohort=='sepsis' else eligible_flag
        eligible=np.where(keep)[0]
        rng=np.random.default_rng(20261024)
        per=int(np.ceil(60000/len(np.unique(sample[eligible]))))
        picked=np.concatenate([rng.choice(eligible[np.asarray(ix)],size=min(per,len(ix)),replace=False)
                               for _,ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
        selected=np.sort(rng.choice(picked,size=min(60000,len(picked)),replace=False))
        symbols=col(f,'var/gene_symbol').astype(str)
        proto=f['uns/generation_protocol']
        kcond=int(proto['K_cond'][()]);cells=int(proto['cells_per_bulk'][()])
        readout=proto['readout'][()].decode()
        valid=(len(sample)==summary['n_input'] and int(keep.sum())==summary['n_eligible'] and
               len(labels)==len(selected)==summary['n_ceiling'] and
               np.array_equal(sample[selected],labels['sample'].astype(str).to_numpy()) and
               np.array_equal(group[selected],labels['group'].astype(str).to_numpy()) and
               f['X/data'].dtype==np.dtype('int32') and len(symbols)==6000 and
               kcond==1000 and cells==4000 and readout=='support-ceiling.v1')
    rows.append({'section':'input','cohort':cohort,'check':'exact 60k selection and protocol','passed':valid,
                 'detail':f'{len(sample)} input; {int(keep.sum())} eligible; {len(selected)} selected; {len(symbols)} genes; K_cond={kcond}'})
    for typ in ('general','myeloid'):
        p=BASE/'figure_recreation'/f'{typ}_umap_points_{cohort}.csv.gz'
        x=pd.read_csv(p)
        labels_ok=set(x.source.astype(str))=={'real_A','ceiling'}
        nceil=int(x.source.eq('ceiling').sum())
        finite=np.isfinite(x[['u1','u2']].to_numpy()).all()
        cell_ok=x.loc[x.source=='ceiling',['sample','fine','six_type']].astype(str).drop_duplicates().merge(
            labels[['sample','fine_tr','six_tr']].astype(str).drop_duplicates().rename(columns={'fine_tr':'fine','six_tr':'six_type'}),
            on=['sample','fine','six_type'],how='left',indicator=True)['_merge'].eq('both').all()
        rows.append({'section':'figure','cohort':cohort,'check':f'{typ} UMAP point labels','passed':labels_ok and finite and cell_ok,
                     'detail':f'{len(x)} plotted, {nceil} ceiling, finite={finite}'})
pd.DataFrame(rows).to_csv(HERE/'source_integrity_checks.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
