"""Measure 30-NN fine-label transfer against held-out real-B author labels."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021'
sys.path.insert(0,str(ROOT/'bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/code'))
from rare_lib import SEPSIS,COVID

rows=[];overall=[]
for cohort,definitions in [('sepsis',SEPSIS),('covid',COVID)]:
    x=pd.read_csv(BASE/'results/big'/f'{cohort}_realU_transfer.csv.gz')
    overall.append({'cohort':cohort,'n_real_B':len(x),'fine_accuracy':float((x.fine==x.fine_tr).mean()),
                    'six_accuracy':float((x.six_type==x.six_tr).mean()),
                    'median_confidence':float(x.conf_tr.median())})
    for subset in definitions:
        truth=x.fine.astype(str).eq(subset)
        pred=x.fine_tr.astype(str).eq(subset)
        tp=int((truth&pred).sum());n_truth=int(truth.sum());n_pred=int(pred.sum())
        precision=tp/n_pred if n_pred else np.nan
        recall=tp/n_truth if n_truth else np.nan
        rows.append({'cohort':cohort,'subset':subset,'n_author':n_truth,'n_transfer':n_pred,
                     'true_positive':tp,'precision':precision,'recall':recall,
                     'f1':2*precision*recall/(precision+recall) if n_pred and n_truth and precision+recall else np.nan})
pd.DataFrame(overall).to_csv(HERE/'holdout_label_accuracy.csv',index=False)
pd.DataFrame(rows).to_csv(HERE/'holdout_rare_label_precision.csv',index=False)
print(pd.DataFrame(overall).to_string(index=False))
print(pd.DataFrame(rows).sort_values(['cohort','n_author'],ascending=[True,False]).to_string(index=False))
