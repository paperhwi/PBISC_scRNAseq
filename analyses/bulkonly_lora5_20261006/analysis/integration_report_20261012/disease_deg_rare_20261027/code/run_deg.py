import numpy as np, pandas as pd, json, os
from scipy import stats
from deg_lib import *
os.makedirs("out",exist_ok=True)
res={}; allde={}
for coh in ["sepsis","covid"]:
    srcs={s:load(f"{s}_{coh}") for s in ["real","floor","prev"]}
    genes=srcs["real"]["genes"]
    # measured genes: non-zero in real somewhere AND non-zero in floor somewhere
    meas=(srcs["real"]["S"].sum(0)>0)&(srcs["floor"]["S"].sum(0)>0)
    gm=genes[meas]
    des={s:design(f"{s}_{coh}",srcs[s]) for s in srcs}
    for ct in SIX:
        tabs={}
        for s in srcs:
            g=cpm_matrix(srcs[s],des[s],ct)
            if len(g)==0: continue
            t=deg(g,list(genes))
            if t is None: continue
            t=t[t.gene.isin(gm)].reset_index(drop=True)
            tabs[s]=t
        if "real" not in tabs: continue
        allde[(coh,ct)]=tabs
        R=tabs["real"]
        # expression rank in real (case+ctrl mean)
        rexp=(R.mean_case+R.mean_ctrl)/2
        rank=pd.Series(stats.rankdata(-rexp),index=R.gene)
        for s in ["floor","prev"]:
            if s not in tabs: continue
            F=tabs[s].set_index("gene").loc[R.gene]
            row=dict(cohort=coh,celltype=ct,source=s,
                     n_case_real=int(R.n_case[0]),n_ctrl_real=int(R.n_ctrl[0]),
                     n_case_gen=int(F.n_case.iloc[0]),n_ctrl_gen=int(F.n_ctrl.iloc[0]),
                     n_genes=len(R))
            row["logFC_r"]=float(np.corrcoef(R.logFC,F.logFC)[0,1])
            row["logFC_rho"]=float(stats.spearmanr(R.logFC,F.logFC).statistic)
            sl=np.polyfit(R.logFC.values,F.logFC.values,1)[0]; row["slope"]=float(sl)
            hi=(rank.values>2000)
            if hi.sum()>50:
                row["logFC_r_rank2000plus"]=float(np.corrcoef(R.logFC[hi],F.logFC[hi])[0,1])
                row["slope_rank2000plus"]=float(np.polyfit(R.logFC[hi],F.logFC[hi],1)[0])
            rd=R.deg.values; fd=F.deg.values
            row["n_DEG_real"]=int(rd.sum()); row["n_DEG_gen"]=int(fd.sum())
            row["recall"]=float((rd&fd).sum()/max(rd.sum(),1)) if rd.sum() else np.nan
            row["precision"]=float((rd&fd).sum()/max(fd.sum(),1)) if fd.sum() else np.nan
            row["jaccard"]=float((rd&fd).sum()/max((rd|fd).sum(),1))
            row["F1"]=float(2*row["recall"]*row["precision"]/(row["recall"]+row["precision"])) if (rd.sum() and fd.sum() and (row["recall"]+row["precision"])>0) else np.nan
            if rd.sum():
                row["sign_agree_realDEG"]=float((np.sign(R.logFC[rd])==np.sign(F.logFC.values[rd])).mean())
            row["sign_agree_all"]=float((np.sign(R.logFC)==np.sign(F.logFC.values)).mean())
            # top-100 by |logFC| among q<0.05
            def top100(T):
                t=T[T.q<0.05]
                return set(t.reindex(t.logFC.abs().sort_values(ascending=False).index).gene.head(100))
            tr=top100(R); tf=top100(F.reset_index())
            row["top100_overlap"]=len(tr&tf)
            res[(coh,ct,s)]=row
D=pd.DataFrame(res.values())
D.to_csv("out/deg_summary.csv",index=False)
pd.set_option("display.width",250)
print(D[["cohort","celltype","source","n_case_real","n_ctrl_real","n_case_gen","n_ctrl_gen","n_DEG_real","n_DEG_gen","logFC_r","logFC_r_rank2000plus","slope","sign_agree_realDEG","recall","precision","top100_overlap"]].round(3).to_string(index=False))
import pickle
pickle.dump(allde,open("out/allde.pkl","wb"))
