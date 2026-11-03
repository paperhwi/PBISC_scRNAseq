"""Disease shift of each rare subset within each source: case vs control (per-patient shares), Welch t; then real-vs-floor agreement of the shift."""
import numpy as np, pandas as pd
from scipy import stats
U=pd.read_csv("perpatient/out/rare_unit_shares.csv"); rows=[]
for (coh,sub,src),d in U.groupby(["cohort","subset","source"]):
    a=d[d.group=="case"].pct_all; b=d[d.group=="control"].pct_all
    rows.append(dict(cohort=coh,subset=sub,source=src,mean_case=a.mean(),mean_ctrl=b.mean(),
                     log2FC=np.log2((a.mean()+0.01)/(b.mean()+0.01)),p=stats.ttest_ind(a,b,equal_var=False).pvalue if a.std()+b.std()>0 else np.nan))
T=pd.DataFrame(rows); W=T.pivot_table(index=["cohort","subset"],columns="source",values=["log2FC","p"]); W.columns=[f"{a}_{b}" for a,b in W.columns]
W["same_direction"]=np.sign(W.log2FC_real)==np.sign(W.log2FC_floor)
W["real_sig"]=W.p_real<0.05; W["floor_sig"]=W.p_floor<0.05
W.reset_index().to_csv("perpatient/out/rare_disease_shift.csv",index=False)
pd.set_option("display.width",200); print(W.round(3).to_string())
for coh in ["sepsis","covid"]:
    w=W.loc[coh]; r=stats.spearmanr(w.log2FC_real,w.log2FC_floor)
    print(coh,"shift rho",round(r.statistic,2),"p",round(r.pvalue,3),"| same direction",int(w.same_direction.sum()),"/",len(w),"| real sig",int(w.real_sig.sum()),"floor sig",int(w.floor_sig.sum()),"both sig & same dir",int((w.real_sig&w.floor_sig&w.same_direction).sum()))
