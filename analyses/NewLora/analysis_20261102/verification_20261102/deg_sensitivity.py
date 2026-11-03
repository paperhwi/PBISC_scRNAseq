"""Recalculate generated DEG with full-panel CPM and crossover patients removed."""
from pathlib import Path
import gc
import numpy as np
import pandas as pd
import anndata as ad

HERE=Path(__file__).resolve().parent
SCRIPT=HERE.parents[1]/'analyze_disease_deg_matched.py'
ns={'__file__':str(SCRIPT)}
exec(SCRIPT.read_text().split('summaries=[];sample_counts=[]')[0],ns)

BASE=ns['BASE']; OLD=ns['OLD']
rows=[]; crossings=[]
for cohort in ('sepsis','covid'):
    real=ad.read_h5ad(BASE/'data/big'/f'{cohort}_ref.h5ad',backed='r')
    genes=np.asarray(real.var_names,dtype=str)
    real.file.close()
    panel=pd.read_csv(OLD/f'{cohort}_Tcell_real.csv.gz').gene.astype(str)
    panel=panel[panel.isin(set(genes))].reset_index(drop=True)
    pos=pd.Index(genes).get_indexer(panel)
    for source in ('prev4k','floor','ceiling'):
        print(cohort,source,flush=True)
        if source=='ceiling':
            X,obs=ns['ceiling'](cohort,genes)
        else:
            d=ad.read_h5ad(BASE/'data/big'/f'{cohort}_{source}.h5ad')
            X=d.X.tocsr()
            obs=pd.read_csv(BASE/'results/big'/f'{cohort}_{source}_transfer.csv.gz',index_col=0)[['sample','group','six_tr']]
            assert d.obs.index.equals(obs.index)
            del d
        full=ns['aggregate'](X,obs,cohort)
        old=ns['aggregate'](X[:,pos],obs,cohort)
        for celltype in full:
            if celltype not in old:continue
            Z,meta,m=full[celltype]
            oldZ,oldmeta,_=old[celltype]
            assert meta.equals(oldmeta)
            current=ns['deg'](oldZ,meta,panel.to_numpy())
            alternate=ns['deg'](Z[:,pos],meta,panel.to_numpy())
            if current is None or alternate is None:continue
            a=current.deg.to_numpy();b=alternate.deg.to_numpy()
            overlap=meta.groupby('patient').condition.nunique()
            both=set(overlap.index[overlap>1])
            mask=~meta.patient.isin(both)
            exclude=ns['deg'](Z[mask.to_numpy()][:,pos],meta[mask].reset_index(drop=True),panel.to_numpy())
            crossing={'cohort':cohort,'source':source,'celltype':celltype,'n_crossover_patients':len(both),
                      'crossover_ids':','.join(sorted(both)),'n_case_current':int(current.n_case.iloc[0]),
                      'n_control_current':int(current.n_ctrl.iloc[0])}
            if exclude is not None:
                crossing.update(n_case_excluding=int(exclude.n_case.iloc[0]),
                                n_control_excluding=int(exclude.n_ctrl.iloc[0]),
                                n_DEG_fullpanel_excluding=int(exclude.deg.sum()),
                                overlap_with_fullpanel=int((exclude.deg&alternate.deg).sum()),
                                logFC_r_vs_fullpanel=float(np.corrcoef(exclude.logFC,alternate.logFC)[0,1]))
            crossings.append(crossing)
            rows.append({'cohort':cohort,'source':source,'celltype':celltype,
                         'n_DEG_current':int(a.sum()),'n_DEG_fullpanel':int(b.sum()),
                         'n_DEG_shared':int((a&b).sum()),
                         'logFC_r':float(np.corrcoef(current.logFC,alternate.logFC)[0,1]),
                         'median_abs_logFC_delta':float(np.median(np.abs(current.logFC-alternate.logFC))),
                         'max_abs_logFC_delta':float(np.max(np.abs(current.logFC-alternate.logFC)))})
        del full,old,X,obs
        gc.collect()
pd.DataFrame(rows).to_csv(HERE/'full_panel_cpm_sensitivity.csv',index=False)
pd.DataFrame(crossings).to_csv(HERE/'crossover_patient_sensitivity.csv',index=False)
print('COMPLETE',flush=True)
