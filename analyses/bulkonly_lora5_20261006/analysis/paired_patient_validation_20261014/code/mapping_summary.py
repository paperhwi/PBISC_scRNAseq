"""Aggregate truly out-of-fold patient scores and recoverability."""
from common import *
from scipy.stats import spearmanr
import json
rng=np.random.default_rng(1407)
files=list((OUT/'results').glob('*_confusion_counts.csv.gz'))
d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
summary=[];labels=[]
for (protocol,source,route,level),g in d.groupby(['protocol','source','route','level']):
    g['correct_n']=np.where(g['true']==g.pred,g.n,0)
    pt=g.groupby('patient')[['n','correct_n']].sum();v=pt.correct_n/pt.n
    mean,lo,hi=boot_mean(v,rng)
    cats=sorted(set(g['true'])|set(g.pred))
    tr=g.groupby(['patient','sample','true']).n.sum().unstack(fill_value=0).reindex(columns=cats,fill_value=0)
    pr=g.groupby(['patient','sample','pred']).n.sum().unstack(fill_value=0).reindex(columns=cats,fill_value=0)
    tr=tr.div(tr.sum(axis=1),axis=0);pr=pr.div(pr.sum(axis=1),axis=0)
    l1=np.abs(pr-tr).sum(axis=1).groupby(level='patient').mean();lm,ll,lh=boot_mean(l1,rng)
    summary.append(dict(protocol=protocol,source=source,route=route,level=level,patient_accuracy=mean,accuracy_lo=lo,accuracy_hi=hi,
                        composition_L1=lm,L1_lo=ll,L1_hi=lh,n_patients=len(pt),n_samples=g['sample'].nunique(),n_cells=int(g.n.sum())))
    # Patient-balanced confusion counts, not visit-count-weighted precision.
    g['w']=g.n/g.patient.map(pt.n)
    c=g.pivot_table(index='true',columns='pred',values='w',aggfunc='sum',fill_value=0)
    for label in cats:
        tp=float(c.loc[label,label]) if label in c.index and label in c else 0
        ntrue=float(c.loc[label].sum()) if label in c.index else 0
        npred=float(c[label].sum()) if label in c else 0
        recall=tp/ntrue if ntrue else np.nan;precision=tp/npred if npred else np.nan
        labels.append(dict(protocol=protocol,source=source,route=route,level=level,label=label,recall=recall,precision=precision,
                           recoverable=bool(recall>=.5 and precision>=.5),n_true_cells=int(g.loc[g['true']==label,'n'].sum())))
pd.DataFrame(summary).to_csv(OUT/'results/covid_mapping_patient_summary.csv',index=False)
pd.DataFrame(labels).to_csv(OUT/'results/covid_mapping_recoverability.csv',index=False)
print(pd.DataFrame(summary).query("source=='realT2'").to_string(index=False))
