"""Specificity law for scSemiProfiler: gene-level (tau vs recovery) and data-driven module level, same definitions as b1b_genelaw.py / b1_modules.py.
Recovery = mean CP10k in the generator / depth-matched real B (thinned to floor depth, identical seed) within the gene's six-type.
Generator six-type = its own label (as for floor)."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.cluster import KMeans
from scipy import stats
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"; SS="/mnt/user-data/uploads/integration_report_20261012/three_model_20261022/data/"
SIX=["T cell","B cell","NK cell","Monocyte","DC"]
gl=[]; ml=[]
def cp10k_mean(X,sel):
    sub=X[sel]; t=np.asarray(sub.sum(1)).ravel(); t[t==0]=1; return np.asarray((sp.diags(1e4/t)@sub).mean(0)).ravel()
for coh in ["sepsis","covid"]:
    rng=np.random.default_rng(20261027)
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad"); F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad"); S=ad.read_h5ad(f"{SS}{coh}_scsemi.h5ad")
    g=np.array(R.var_names); fineA=R.obs.fine.astype(str).to_numpy(); sixA=R.obs.six_type.astype(str).to_numpy()
    XR=sp.csr_matrix(R.X); t=np.asarray(XR.sum(1)).ravel(); t[t==0]=1; CR=sp.diags(1e4/t)@XR
    types=[x for x in pd.unique(fineA) if (fineA==x).sum()>=100]
    M=np.vstack([np.asarray(CR[fineA==x].mean(0)).ravel() for x in types])
    xr=M/np.maximum(M.max(0,keepdims=True),1e-9); tau=(1-xr).sum(0)/(len(types)-1)
    sixF=F.obs.six_type.astype(str).to_numpy(); gt=np.asarray(sp.csr_matrix(F.X).sum(1)).ravel()
    XU=sp.csr_matrix(U.X).astype(np.int64); sixU=U.obs.six_type.astype(str).to_numpy()
    pool={c:gt[sixF==c] for c in np.unique(sixF)}
    tg=np.array([rng.choice(pool.get(c,gt)) for c in sixU],float); kp=np.minimum(1.0,tg/np.maximum(np.asarray(XU.sum(1)).ravel(),1))
    XU.data=rng.binomial(XU.data,np.repeat(kp,np.diff(XU.indptr))); XU.eliminate_zeros(); XU=sp.csr_matrix(XU.astype(np.float32))
    SRC={"real_T1":(XU,sixU),"floor":(sp.csr_matrix(F.X),sixF),"scsemi":(sp.csr_matrix(S.X),S.obs.six_type.astype(str).to_numpy())}
    for ct in SIX:
        mA=np.asarray(CR[sixA==ct].mean(0)).ravel()
        v={s:(cp10k_mean(X,six==ct) if (six==ct).sum()>=100 else None) for s,(X,six) in SRC.items()}
        if any(v[s] is None for s in v): continue
        ok=(v["real_T1"]>=0.2)&(mA>=0.2)
        gl.append(pd.DataFrame(dict(cohort=coh,celltype=ct,gene=g[ok],tau=tau[ok],rank_in_type=stats.rankdata(-mA)[ok],
                                    floor_ratio=v["floor"][ok]/v["real_T1"][ok],scsemi_ratio=v["scsemi"][ok]/v["real_T1"][ok])))
    # modules
    keep=M.max(0)>=0.5; gk=g[keep]; Mk=M[:,keep]; x=Mk/np.maximum(Mk.max(0,keepdims=True),1e-9); tk=(1-x).sum(0)/(len(types)-1)
    Z=(x-x.mean(0))/np.maximum(x.std(0),1e-9); lab=KMeans(n_clusters=30,n_init=10,random_state=0).fit(Z.T).labels_
    six_of_fine=pd.Series(sixA,index=fineA).groupby(level=0).agg(lambda s:s.mode().iloc[0]); gidx={x:i for i,x in enumerate(g)}
    for k in range(30):
        mg=gk[lab==k]
        if len(mg)<5 or len(mg)>400: continue
        ii=[gidx[q] for q in mg]; dom_type=types[int(np.argmax(Mk[:,lab==k].mean(1)))]; dom6=six_of_fine.get(dom_type,"Other")
        if dom6 not in SIX: continue
        vals={}
        for s,(X,six) in SRC.items():
            sel=six==dom6
            if sel.sum()<50: vals[s]=np.nan; continue
            sub=X[sel]; tt=np.asarray(sub.sum(1)).ravel(); tt[tt==0]=1
            vals[s]=float(np.asarray(sp.diags(1e4/tt).dot(sub)[:,ii].sum(1)).ravel().mean())
        if not vals["real_T1"]>0: continue
        ml.append(dict(cohort=coh,module=f"{coh}_M{k:02d}",n_genes=len(mg),dom_six=dom6,tau=float(tk[lab==k].mean()),
                       floor_ratio=vals["floor"]/vals["real_T1"],scsemi_ratio=vals["scsemi"]/vals["real_T1"],top_genes=",".join(mg[:8])))
    print(coh,"done",flush=True)
G=pd.concat(gl,ignore_index=True); Mo=pd.DataFrame(ml)
G.to_csv("scsemi/out/law_gene_level.csv.gz",index=False,compression="gzip"); Mo.to_csv("scsemi/out/law_modules.csv",index=False)
out=[]
for lvl,D in [("gene",G),("module",Mo)]:
    for s in ["floor","scsemi"]:
        d=D.replace([np.inf,-np.inf],np.nan).dropna(subset=[f"{s}_ratio"])
        r=stats.spearmanr(d.tau,d[f"{s}_ratio"]); out.append(dict(level=lvl,generator=s,n=len(d),rho=r.statistic,p=r.pvalue,
            median_ratio=d[f"{s}_ratio"].median(),q1_broad=d[d.tau<=d.tau.quantile(.2)][f"{s}_ratio"].median(),q5_specific=d[d.tau>=d.tau.quantile(.8)][f"{s}_ratio"].median()))
O=pd.DataFrame(out); O.to_csv("scsemi/out/law_summary.csv",index=False); print(O.round(4).to_string(index=False))
