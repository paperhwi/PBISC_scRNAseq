"""Compare archived BH family with one shared, measured-gene BH family."""
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
DIR=HERE.parent/'disease_deg_matched/deg_tables'

def bh(p):
    p=np.asarray(p,float);ok=np.isfinite(p);q=np.full(len(p),np.nan)
    n=ok.sum()
    if not n:return q
    ix=np.argsort(p[ok]);v=p[ok][ix]
    z=np.minimum.accumulate((v*n/np.arange(1,n+1))[::-1])[::-1]
    r=np.empty(n);r[ix]=np.minimum(z,1);q[ok]=r
    return q

rows=[]
for cohort in ('sepsis','covid'):
    for ct in ('Tcell','Bcell','NKcell','Monocyte'):
        rt=pd.read_csv(DIR/f'{cohort}_{ct}_real.csv.gz')
        r_old=rt.deg.to_numpy(bool)
        r_new=(bh(rt.p)<.05)&(rt.logFC.abs().to_numpy()>=1)
        for source in ('prev4k','floor','ceiling'):
            tab=pd.read_csv(DIR/f'{cohort}_{ct}_{source}.csv.gz')
            assert tab.gene.equals(rt.gene)
            g_old=tab.deg.to_numpy(bool)
            g_new=(bh(tab.p)<.05)&(tab.logFC.abs().to_numpy()>=1)
            old_overlap=int((r_old&g_old).sum());new_overlap=int((r_new&g_new).sum())
            rows.append({'cohort':cohort,'celltype':ct,'source':source,
                         'real_DEG_archived_BH':int(r_old.sum()),'real_DEG_common_BH':int(r_new.sum()),
                         'gen_DEG_archived_BH':int(g_old.sum()),'gen_DEG_common_BH':int(g_new.sum()),
                         'overlap_archived_BH':old_overlap,'overlap_common_BH':new_overlap,
                         'recall_archived_BH':old_overlap/max(r_old.sum(),1),
                         'recall_common_BH':new_overlap/max(r_new.sum(),1),
                         'precision_archived_BH':old_overlap/max(g_old.sum(),1),
                         'precision_common_BH':new_overlap/max(g_new.sum(),1)})
out=pd.DataFrame(rows)
out.to_csv(HERE/'bh_family_sensitivity.csv',index=False)
print(out.to_string(index=False))
