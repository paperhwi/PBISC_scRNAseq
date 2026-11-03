"""B1 - data-driven co-expression modules defined on real half A only, scored blind in floor."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.cluster import KMeans
from scipy import stats
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
rng=np.random.default_rng(20261027)
SIX=["T cell","B cell","NK cell","Monocyte","DC"]
rows=[]; modtab=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad"); F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad")
    P4=ad.read_h5ad(f"{IN}{coh}_prev4k.h5ad")
    g=np.array(R.var_names)
    fineA=R.obs.fine.astype(str).to_numpy(); sixA=R.obs.six_type.astype(str).to_numpy()
    XR=sp.csr_matrix(R.X); tR=np.asarray(XR.sum(1)).ravel(); tR[tR==0]=1
    CR=sp.diags(1e4/tR)@XR
    # ---- gene x fine-celltype mean profile on real A (types with >=100 cells)
    types=[t for t in pd.unique(fineA) if (fineA==t).sum()>=100]
    M=np.vstack([np.asarray(CR[fineA==t].mean(0)).ravel() for t in types])     # types x genes
    keep=M.max(0)>=0.5
    gk=g[keep]; Mk=M[:,keep]
    # Tau specificity
    x=Mk/np.maximum(Mk.max(0,keepdims=True),1e-9)
    tau=(1-x).sum(0)/(len(types)-1)
    # ---- modules: cluster genes by their across-celltype profile shape
    Z=(x-x.mean(0))/np.maximum(x.std(0),1e-9)
    km=KMeans(n_clusters=30,n_init=10,random_state=0).fit(Z.T)
    lab=km.labels_
    # ---- score each module in each source, inside its dominant six-type
    six_of_fine=pd.Series(sixA,index=fineA).groupby(level=0).agg(lambda s:s.mode().iloc[0])
    gt=np.asarray(sp.csr_matrix(F.X).sum(1)).ravel(); sixF=F.obs.six_type.astype(str).to_numpy()
    def thin(X,six):
        Xb=sp.csr_matrix(X,copy=True).astype(np.int64); tot=np.asarray(Xb.sum(1)).ravel()
        pool={c:gt[sixF==c] for c in np.unique(sixF)}
        tg=np.array([rng.choice(pool.get(c,gt)) for c in six],float)
        kp=np.minimum(1.0,tg/np.maximum(tot,1)); Xb.data=rng.binomial(Xb.data,np.repeat(kp,np.diff(Xb.indptr))); Xb.eliminate_zeros()
        return sp.csr_matrix(Xb.astype(np.float32))
    XU=sp.csr_matrix(U.X); sixU=U.obs.six_type.astype(str).to_numpy()
    SRC={"real_T1":(thin(XU,sixU),sixU),"floor":(sp.csr_matrix(F.X),sixF),
         "prev":(sp.csr_matrix(P4.X),P4.obs.six_type.astype(str).to_numpy())}
    gidx={x:i for i,x in enumerate(g)}
    for k in range(30):
        mg=gk[lab==k]
        if len(mg)<5 or len(mg)>400: continue
        ii=[gidx[x] for x in mg]
        dom_type=types[int(np.argmax(Mk[:,lab==k].mean(1)))]
        dom6=six_of_fine.get(dom_type,"Other")
        if dom6 not in SIX: continue
        tau_m=float(np.mean(tau[lab==k]))
        # mean rank of the module in real A within dom6
        subA=CR[sixA==dom6]
        mA=np.asarray(subA.mean(0)).ravel()
        rk=stats.rankdata(-mA); rank_m=float(np.median(rk[ii]))
        ref=None; vals={}
        for s,(X,six) in SRC.items():
            sel=six==dom6
            if sel.sum()<50: vals[s]=np.nan; continue
            sub=X[sel]; t=np.asarray(sub.sum(1)).ravel(); t[t==0]=1
            cp=np.asarray(sp.diags(1e4/t).dot(sub)[:,ii].sum(1)).ravel()
            vals[s]=float(cp.mean())
        if not np.isfinite(vals.get("real_T1",np.nan)) or vals["real_T1"]<=0: continue
        rows.append(dict(cohort=coh,module=f"{coh}_M{k:02d}",n_genes=len(mg),dom_type=dom_type,dom_six=dom6,
                         tau=tau_m,median_rank_in_type=rank_m,
                         real_T1=vals["real_T1"],floor=vals["floor"],prev=vals["prev"],
                         floor_ratio=vals["floor"]/vals["real_T1"],prev_ratio=vals["prev"]/vals["real_T1"],
                         top_genes=",".join(mg[np.argsort(-Mk[:,lab==k].mean(0))][:8])))
B=pd.DataFrame(rows); B.to_csv("out2/b1_datadriven_modules.csv",index=False)
pd.set_option("display.width",260); pd.set_option("display.max_colwidth",60)
print(f"modules scored: {len(B)}")
print(B[["cohort","module","n_genes","dom_six","tau","median_rank_in_type","floor_ratio","prev_ratio","top_genes"]].sort_values("tau").round(3).to_string(index=False))
for coh in ["sepsis","covid","ALL"]:
    d=B if coh=="ALL" else B[B.cohort==coh]
    if len(d)<5: continue
    r=stats.spearmanr(d.tau,d.floor_ratio); r2=stats.spearmanr(d.tau,d.prev_ratio)
    print(f"\n{coh}: Spearman(specificity tau, floor recovery) = {r.statistic:.3f} (p={r.pvalue:.1e}, n={len(d)});  prev = {r2.statistic:.3f} (p={r2.pvalue:.1e})")
hi=B[B.tau>=B.tau.quantile(.75)]; lo=B[B.tau<=B.tau.quantile(.25)]
print(f"\nhigh-specificity modules (top quartile tau): median floor recovery {hi.floor_ratio.median():.2f}  (n={len(hi)})")
print(f"low-specificity  modules (bottom quartile tau): median floor recovery {lo.floor_ratio.median():.2f}  (n={len(lo)})")
print("Mann-Whitney p =", f"{stats.mannwhitneyu(hi.floor_ratio,lo.floor_ratio).pvalue:.2e}")
print("\nSpearman(expression rank in type, floor recovery) =",
      f"{stats.spearmanr(B.median_rank_in_type,B.floor_ratio).statistic:.3f}")
