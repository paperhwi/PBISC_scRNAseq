"""Recalculate reported post-Harmony iLISI from saved, plotted embeddings."""
from pathlib import Path
import sys,json
import anndata as ad
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
PBISC=HERE.parents[2]
sys.path.insert(0,str(PBISC/'bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/code'))
import ilisi_lib as L

rows=[]
for cohort in ('sepsis','covid'):
    f=HERE.parent/cohort/'joint_harmony_embedding.h5ad'
    a=ad.read_h5ad(f)
    src=(a.obs.source.astype(str).to_numpy()=='ceiling').astype(np.int8)
    target=json.loads((HERE.parent/cohort/'summary.json').read_text())
    assert len(src)==target['n_ref']+target['n_ceiling']
    for stage,col in [('pre','X_pca'),('post','X_pca_harmony')]:
        scores=[L.lisi_stats(a.obsm[col],src,seed=s) for s in (0,1,2)]
        median=float(np.mean([x[0] for x in scores]))
        query_mean=float(np.mean([x[1] for x in scores]))
        own=float(np.mean([x[2] for x in scores]))
        rows.append({'cohort':cohort,'stage':stage,'reported_iLISI':target[f'iLISI_{stage}'],
                     'recomputed_iLISI':median,'absolute_difference':abs(median-target[f'iLISI_{stage}']),
                     'reported_own_neighbor_share':target[f'own_share_{stage}'],
                     'recomputed_own_neighbor_share':own,'query_mean_difference':query_mean-target[f'iLISI_{stage}_query_mean']})
        print(cohort,stage,median,flush=True)
pd.DataFrame(rows).to_csv(HERE/'ilisi_recomputation.csv',index=False)
