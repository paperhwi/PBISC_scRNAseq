import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
rng=np.random.default_rng(20261027)
SET=[("C1Q / complement","sepsis","Monocyte",["C1QA","C1QB","C1QC"]),
     ("plasma-cell secretory","sepsis","B cell",["JCHAIN","MZB1","IGLL5","DERL3","TXNDC5"]),
     ("pDC / cDC identity","sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C"]),
     ("sepsis myeloid activation","sepsis","Monocyte",["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1"]),
     ("ISG (COVID)","covid","Monocyte",["ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1"])]
rows=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad"); F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad")
    P=ad.read_h5ad(f"{IN}{coh}_prev4k.h5ad")
    g=np.array(R.var_names); gi={x:i for i,x in enumerate(g)}
    XF=sp.csr_matrix(F.X); labF=F.obs.six_type.astype(str).to_numpy()
    Xb=sp.csr_matrix(U.X).copy().astype(np.int64); tot=np.asarray(Xb.sum(1)).ravel()
    gt=np.asarray(XF.sum(1)).ravel(); pool={c:gt[labF==c] for c in np.unique(labF)}
    tgt=np.array([rng.choice(pool.get(c,gt)) for c in U.obs.six_type.astype(str)],float)
    kp=np.minimum(1.0,tgt/np.maximum(tot,1)); Xb.data=rng.binomial(Xb.data,np.repeat(kp,np.diff(Xb.indptr))); Xb.eliminate_zeros()
    SRC={"real A":(sp.csr_matrix(R.X),R.obs.six_type.astype(str).to_numpy()),
         "real (depth-matched)":(sp.csr_matrix(Xb.astype(np.float32)),U.obs.six_type.astype(str).to_numpy()),
         "floor":(XF,labF),
         "prev K4000 (Poisson)":(sp.csr_matrix(P.X),P.obs.six_type.astype(str).to_numpy())}
    for name,c,ct,gs in SET:
        if c!=coh: continue
        ii=[gi[x] for x in gs if x in gi]
        for lab,(X,six) in SRC.items():
            sub=X[six==ct]; t=np.asarray(sub.sum(1)).ravel(); t[t==0]=1
            cp=sp.diags(1e4/t).dot(sub)[:,ii]
            raw=sub[:,ii]
            pos=np.asarray((raw>0).sum(1)).ravel()>0
            sc=np.asarray(cp.sum(1)).ravel()
            umi=np.asarray(raw.sum(1)).ravel()
            rows.append(dict(module=name,celltype=ct,source=lab,n_cells=int((six==ct).sum()),
                pct_pos=100*pos.mean(),
                mean_CP10k_all=float(sc.mean()),
                mean_CP10k_in_pos=float(sc[pos].mean()) if pos.any() else np.nan,
                mean_UMI_in_pos=float(umi[pos].mean()) if pos.any() else np.nan))
D=pd.DataFrame(rows); D.to_csv("out2/c2b_module_load.csv",index=False)
pd.set_option("display.width",230)
print(D.round(2).to_string(index=False))
print("\n=== floor / depth-matched-real ratios ===")
p=D.pivot_table(index=["module"],columns="source",values=["pct_pos","mean_CP10k_in_pos"])
out=pd.DataFrame({"pos% floor/real":p[("pct_pos","floor")]/p[("pct_pos","real (depth-matched)")],
                  "load floor/real":p[("mean_CP10k_in_pos","floor")]/p[("mean_CP10k_in_pos","real (depth-matched)")],
                  "pos% prev/real":p[("pct_pos","prev K4000 (Poisson)")]/p[("pct_pos","real (depth-matched)")],
                  "load prev/real":p[("mean_CP10k_in_pos","prev K4000 (Poisson)")]/p[("mean_CP10k_in_pos","real (depth-matched)")]})
print(out.round(3).to_string())
