import numpy as np, pandas as pd
from rare_lib import *
rng=np.random.default_rng(99)
rows=[]
for coh in ["sepsis","covid"]:
    X,V,SRC,genes,(DEF,PAR)=load_coh(coh)
    mA,fA=V["refA"]; detA=np.asarray((X["refA"]>0).mean(0)).ravel(); gok=detA>=0.01
    for sub,fn in DEF.items():
        idx={}
        parA=fn(mA,fA); idx["refA"]=np.where(parA&(fA==sub))[0]
        for lab in ["realB_tr","floor_tr"]:
            m,f=V[lab]; p=fn(m,f); idx[lab]=np.where(p&(f==sub))[0]
        n=min(len(v) for v in idx.values())
        if n<40: continue
        res=dict(cohort=coh,subset=sub,n_matched=n,n_refA=len(idx["refA"]),n_realB=len(idx["realB_tr"]),n_floor=len(idx["floor_tr"]))
        reps=[]
        for _ in range(10):
            pr={}
            for lab,key in [("refA","refA"),("realB_tr","realB"),("floor_tr","floor")]:
                s=rng.choice(idx[lab],n,replace=False)
                pr[key]=np.asarray(X[SRC.get(lab,"refA")][s].mean(0)).ravel()[gok] if lab!="refA" else np.asarray(X["refA"][s].mean(0)).ravel()[gok]
            # internal real A split for the true noise floor
            s2=rng.permutation(idx["refA"])
            if len(s2)>=2*n:
                a1=np.asarray(X["refA"][s2[:n]].mean(0)).ravel()[gok]; a2=np.asarray(X["refA"][s2[n:2*n]].mean(0)).ravel()[gok]
                reps.append((np.corrcoef(a1,a2)[0,1],(np.abs(a2-a1)>=0.5).sum(),
                             np.corrcoef(pr["refA"],pr["realB"])[0,1],(np.abs(pr["realB"]-pr["refA"])>=0.5).sum(),
                             np.corrcoef(pr["refA"],pr["floor"])[0,1],(np.abs(pr["floor"]-pr["refA"])>=0.5).sum()))
            else:
                reps.append((np.nan,np.nan,np.corrcoef(pr["refA"],pr["realB"])[0,1],(np.abs(pr["realB"]-pr["refA"])>=0.5).sum(),
                             np.corrcoef(pr["refA"],pr["floor"])[0,1],(np.abs(pr["floor"]-pr["refA"])>=0.5).sum()))
        a=np.array(reps,float)
        res.update(r_AA=np.nanmean(a[:,0]),nD_AA=np.nanmean(a[:,1]),r_realB=np.nanmean(a[:,2]),nD_realB=np.nanmean(a[:,3]),
                   r_floor=np.nanmean(a[:,4]),nD_floor=np.nanmean(a[:,5]))
        rows.append(res)
P=pd.DataFrame(rows); P.to_csv("out/rare_power_matched.csv",index=False)
pd.set_option("display.width",250)
print(P.round(3).to_string(index=False))
print("\nmedians: r_AA %.4f r_realB %.4f r_floor %.4f | nD_AA %.0f nD_realB %.0f nD_floor %.0f"%(
 P.r_AA.median(),P.r_realB.median(),P.r_floor.median(),P.nD_AA.median(),P.nD_realB.median(),P.nD_floor.median()))
