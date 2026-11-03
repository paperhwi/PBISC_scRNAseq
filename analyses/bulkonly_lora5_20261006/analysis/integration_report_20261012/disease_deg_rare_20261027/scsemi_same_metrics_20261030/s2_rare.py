"""Rare populations in scSemiProfiler vs floor, same definitions as run_rare.py (f_rare_parent subsets, kNN-transferred labels from real half A
in the pairwise Harmony space; scSemi transfer from three_model big_joint.py). Yield (% of all cells, ratio to real A) and marker score ratio."""
import numpy as np, pandas as pd, anndata as ad, os
from rare_lib import SEPSIS, COVID, PAR_S, PAR_C, lognorm, IN, TR
SS="/mnt/user-data/uploads/integration_report_20261012/three_model_20261022/data/"
rows=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad"); F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad"); S=ad.read_h5ad(f"{SS}{coh}_scsemi.h5ad")
    assert list(S.var_names)==list(R.var_names)
    tr={s:pd.read_csv(f"{TR}{coh}_{s}_transfer.csv.gz",index_col=0,keep_default_na=False) for s in ("realU","floor","scsemi")}
    genes=np.array(R.var_names)
    X={"refA":lognorm(R.X),"realB":lognorm(U.X),"floor":lognorm(F.X),"scsemi":lognorm(S.X)}
    V={"refA":(R.obs.major.astype(str).to_numpy(),R.obs.fine.astype(str).to_numpy()),
       "realB":(tr["realU"].major_tr.to_numpy(),tr["realU"].fine_tr.to_numpy()),
       "floor":(tr["floor"].major_tr.to_numpy(),tr["floor"].fine_tr.to_numpy()),
       "scsemi":(tr["scsemi"].major_tr.to_numpy(),tr["scsemi"].fine_tr.to_numpy())}
    conf={k:tr[k2].conf_tr.to_numpy() for k,k2 in [("realB","realU"),("floor","floor"),("scsemi","scsemi")]}
    DEF,PAR=(SEPSIS,PAR_S) if coh=="sepsis" else (COVID,PAR_C)
    mA,fA=V["refA"]
    for sub,fn in DEF.items():
        parA=fn(mA,fA); tgA=parA&(fA==sub)
        if tgA.sum()<30: continue
        A=X["refA"]
        d=np.asarray(A[tgA].mean(0)).ravel()-np.asarray(A[parA&~tgA].mean(0)).ravel()
        det=np.asarray((A[tgA]>0).mean(0)).ravel(); d2=d.copy(); d2[det<0.1]=-np.inf; top=np.argsort(-d2)[:10]
        r=dict(cohort=coh,subset=sub,parent=PAR[sub],markers=",".join(genes[top]))
        for k in ["refA","realB","floor","scsemi"]:
            m,f=V[k]; t=(fn(m,f)&(f==sub)) if k!="refA" else tgA
            r[f"n_{k}"]=int(t.sum()); r[f"pct_{k}"]=100*t.sum()/len(f)
            if k!="refA": r[f"conf_{k}"]=float(np.median(conf[k][t])) if t.sum() else np.nan
            r[f"mk_{k}"]=float(np.asarray(X[k][t][:,top].mean(1)).ravel().mean()) if t.sum()>=10 else np.nan
        for k in ["realB","floor","scsemi"]:
            r[f"yield_{k}"]=r[f"pct_{k}"]/r["pct_refA"]; r[f"mkratio_{k}"]=r[f"mk_{k}"]/r["mk_refA"]
        rows.append(r)
    print(coh,"done",flush=True)
D=pd.DataFrame(rows); os.makedirs("scsemi/out",exist_ok=True); D.to_csv("scsemi/out/rare_yield_marker.csv",index=False)
pd.set_option("display.width",250)
print(D[["cohort","subset","n_refA","n_floor","n_scsemi","yield_realB","yield_floor","yield_scsemi","mkratio_realB","mkratio_floor","mkratio_scsemi","conf_floor","conf_scsemi"]].round(2).to_string(index=False))
for c in ["realB","floor","scsemi"]:
    print(c,"median yield",D[f"yield_{c}"].median().round(2),"| n subsets <0.5x",(D[f"yield_{c}"]<0.5).sum(),"/",len(D),"| median mk ratio",D[f"mkratio_{c}"].median().round(2))
