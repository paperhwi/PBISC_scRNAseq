import numpy as np, pandas as pd, pickle
from scipy import stats
from deg_lib import *
rng=np.random.default_rng(20261027)
rows=[]
for coh in ["sepsis","covid"]:
    d=load(f"real_{coh}"); des=design(f"real_{coh}",d); genes=d["genes"]
    fl=load(f"floor_{coh}")
    meas=(d["S"].sum(0)>0)&(fl["S"].sum(0)>0); gm=set(genes[meas])
    for ct in SIX:
        g=cpm_matrix(d,des,ct)
        if len(g)==0: continue
        case=g[g._c=="case"]._u.unique(); ctrl=g[g._c=="ctrl"]._u.unique()
        if len(case)<8 or len(ctrl)<6:
            rows.append(dict(cohort=coh,celltype=ct,note=f"too few units case={len(case)} ctrl={len(ctrl)}")); continue
        acc=[]
        for rep in range(20):
            ca=rng.permutation(case); co=rng.permutation(ctrl)
            h=[(ca[:len(ca)//2],co[:len(co)//2]),(ca[len(ca)//2:],co[len(co)//2:])]
            ts=[]
            for cs,cl in h:
                sub=g[g._u.isin(list(cs)+list(cl))]
                t=deg(sub,list(genes))
                if t is None: ts=[]; break
                ts.append(t[t.gene.isin(gm)].reset_index(drop=True))
            if len(ts)!=2: continue
            A,B=ts
            r=np.corrcoef(A.logFC,B.logFC)[0,1]
            rec=(A.deg&B.deg).sum()/max(A.deg.sum(),1) if A.deg.sum() else np.nan
            pre=(A.deg&B.deg).sum()/max(B.deg.sum(),1) if B.deg.sum() else np.nan
            sa=(np.sign(A.logFC[A.deg])==np.sign(B.logFC[A.deg])).mean() if A.deg.sum() else np.nan
            acc.append((r,rec,pre,sa,A.deg.sum(),B.deg.sum()))
        if not acc: continue
        a=np.array(acc,float)
        rows.append(dict(cohort=coh,celltype=ct,n_case=len(case),n_ctrl=len(ctrl),reps=len(acc),
                         logFC_r=np.nanmean(a[:,0]),logFC_r_sd=np.nanstd(a[:,0]),
                         recall=np.nanmean(a[:,1]),precision=np.nanmean(a[:,2]),
                         sign_agree=np.nanmean(a[:,3]),nDEG_halfA=np.nanmean(a[:,4]),nDEG_halfB=np.nanmean(a[:,5])))
C=pd.DataFrame(rows); C.to_csv("out/real_real_ceiling.csv",index=False)
pd.set_option("display.width",250)
print(C.round(3).to_string(index=False))
