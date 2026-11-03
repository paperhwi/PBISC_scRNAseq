"""B1b - gene-level law: recovery vs cell-type specificity (all measured genes)."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy import stats
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
rng=np.random.default_rng(20261027)
SIX=["T cell","B cell","NK cell","Monocyte","DC"]
out=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad")
    F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad"); P4=ad.read_h5ad(f"{IN}{coh}_prev4k.h5ad")
    g=np.array(R.var_names); fineA=R.obs.fine.astype(str).to_numpy(); sixA=R.obs.six_type.astype(str).to_numpy()
    XR=sp.csr_matrix(R.X); t=np.asarray(XR.sum(1)).ravel(); t[t==0]=1; CR=sp.diags(1e4/t)@XR
    types=[x for x in pd.unique(fineA) if (fineA==x).sum()>=100]
    M=np.vstack([np.asarray(CR[fineA==x].mean(0)).ravel() for x in types])
    xr=M/np.maximum(M.max(0,keepdims=True),1e-9); tau=(1-xr).sum(0)/(len(types)-1)
    sixF=F.obs.six_type.astype(str).to_numpy(); gt=np.asarray(sp.csr_matrix(F.X).sum(1)).ravel()
    XU=sp.csr_matrix(U.X).astype(np.int64); sixU=U.obs.six_type.astype(str).to_numpy()
    pool={c:gt[sixF==c] for c in np.unique(sixF)}
    tg=np.array([rng.choice(pool.get(c,gt)) for c in sixU],float)
    kp=np.minimum(1.0,tg/np.maximum(np.asarray(XU.sum(1)).ravel(),1))
    XU.data=rng.binomial(XU.data,np.repeat(kp,np.diff(XU.indptr))); XU.eliminate_zeros()
    SRC={"real_T1":(sp.csr_matrix(XU.astype(np.float32)),sixU),"floor":(sp.csr_matrix(F.X),sixF),
         "prev":(sp.csr_matrix(P4.X),P4.obs.six_type.astype(str).to_numpy())}
    for ct in SIX:
        mA=np.asarray(CR[sixA==ct].mean(0)).ravel()
        v={}
        for s,(X,six) in SRC.items():
            sel=six==ct
            if sel.sum()<100: v[s]=None; continue
            sub=X[sel]; tt=np.asarray(sub.sum(1)).ravel(); tt[tt==0]=1
            v[s]=np.asarray((sp.diags(1e4/tt)@sub).mean(0)).ravel()
        if v["real_T1"] is None or v["floor"] is None: continue
        ok=(v["real_T1"]>=0.2)&(mA>=0.2)
        out.append(pd.DataFrame(dict(cohort=coh,celltype=ct,gene=g[ok],tau=tau[ok],
            mean_realT1=v["real_T1"][ok],
            rank_in_type=stats.rankdata(-mA)[ok],
            floor_ratio=v["floor"][ok]/v["real_T1"][ok],
            prev_ratio=(v["prev"][ok]/v["real_T1"][ok]) if v["prev"] is not None else np.nan)))
G=pd.concat(out,ignore_index=True)
G["log2_floor_ratio"]=np.log2(np.maximum(G.floor_ratio,1e-3))
G.to_csv("out2/b1b_gene_level_law.csv.gz",index=False,compression="gzip")
pd.set_option("display.width",250)
print("genes x celltype scored:",len(G))
r=stats.spearmanr(G.tau,G.floor_ratio); r2=stats.spearmanr(G.rank_in_type,G.floor_ratio)
print(f"Spearman(tau, floor recovery) = {r.statistic:.3f} (p={r.pvalue:.1e})")
print(f"Spearman(expression rank, floor recovery) = {r2.statistic:.3f} (p={r2.pvalue:.1e})   <- confounder check")
# partial: within expression-rank bands
print("\n=== floor recovery by specificity quintile, within expression-rank bands ===")
G["rank_band"]=pd.cut(G.rank_in_type,[0,500,1500,3000,99999],labels=["1-500","501-1500","1501-3000",">3000"])
G["tau_q"]=pd.qcut(G.tau,5,labels=["Q1 broad","Q2","Q3","Q4","Q5 specific"])
print(G.pivot_table(index="rank_band",columns="tau_q",values="floor_ratio",aggfunc="median",observed=True).round(3).to_string())
print("\n=== same, prev K4000 ===")
print(G.pivot_table(index="rank_band",columns="tau_q",values="prev_ratio",aggfunc="median",observed=True).round(3).to_string())
print("\n=== by cohort / cell type (median floor recovery, Q1 broad vs Q5 specific) ===")
p=G.pivot_table(index=["cohort","celltype"],columns="tau_q",values="floor_ratio",aggfunc="median",observed=True)
p["Q5/Q1"]=p["Q5 specific"]/p["Q1 broad"]; print(p.round(3).to_string())
