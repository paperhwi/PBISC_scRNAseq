import numpy as np, pandas as pd, pickle, anndata as ad, scipy.sparse as sp
from scipy import stats
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
MODG={"C1Q / complement":("sepsis","Monocyte",["C1QA","C1QB","C1QC","RNASE1","CLU"]),
 "plasma-cell secretory":("sepsis","B cell",["JCHAIN","MZB1","IGLL5","DERL3","TNFRSF17","TXNDC5"]),
 "pDC / cDC identity":("sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C"]),
 "sepsis myeloid activation":("sepsis","Monocyte",["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1"]),
 "ISG (COVID)":("covid","Monocyte",["ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1"])}
TR="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/results/big/"
rows=[]; conc=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad")
    U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad")
    rng=np.random.default_rng(20261027)
    Xb=sp.csr_matrix(U.X).copy().astype(np.int64)
    tot=np.asarray(Xb.sum(1)).ravel(); Fx=ad.read_h5ad(f"{IN}{coh}_floor.h5ad")
    gt=np.asarray(sp.csr_matrix(Fx.X).sum(1)).ravel(); labF=Fx.obs.six_type.astype(str).to_numpy()
    pool={c:gt[labF==c] for c in np.unique(labF)}
    tgt=np.array([rng.choice(pool.get(c,gt)) for c in U.obs.six_type.astype(str)],float)
    kp=np.minimum(1.0,tgt/np.maximum(tot,1))
    Xb.data=rng.binomial(Xb.data,np.repeat(kp,np.diff(Xb.indptr))); Xb.eliminate_zeros()
    T1=ad.AnnData(Xb.astype(np.float32),obs=U.obs.copy(),var=U.var.copy())
    genes=np.array(R.var_names); gi={g:i for i,g in enumerate(genes)}
    sixR=R.obs.six_type.astype(str).to_numpy(); sixF=F.obs.six_type.astype(str).to_numpy(); sixT=T1.obs.six_type.astype(str).to_numpy()
    XR=sp.csr_matrix(R.X); XF=sp.csr_matrix(F.X); XT=sp.csr_matrix(T1.X)
    for name,(c,ct,gs) in MODG.items():
        if c!=coh: continue
        mR=sixR==ct; mF=sixF==ct; mT=sixT==ct
        # expression rank within that cell type in real A (CP10k mean)
        sub=XR[mR]; tot=np.asarray(sub.sum(1)).ravel(); tot[tot==0]=1
        cp=np.asarray(sp.diags(1e4/tot).dot(sub).mean(0)).ravel()
        rank=stats.rankdata(-cp)
        for g in gs:
            i=gi.get(g)
            if i is None: continue
            dR=float((XR[mR][:,i]>0).mean()); dT=float((XT[mT][:,i]>0).mean()); dF=float((XF[mF][:,i]>0).mean())
            rows.append(dict(module=name,cohort=coh,celltype=ct,gene=g,
                             rank_in_type=int(rank[i]),rank_pct=float(rank[i]/len(genes)),
                             mean_CP10k_realA=float(cp[i]),
                             det_realA=dR,det_realT1=dT,det_floor=dF,
                             det_ratio_floor_over_T1=float(dF/dT) if dT>0 else np.nan))
    # ---- concentration: where do module UMIs sit inside the parent lineage? ----
    if coh=="sepsis":
        for name,ct,gs in [("plasma-cell secretory","B cell",["JCHAIN","MZB1","IGLL5","DERL3","TXNDC5"]),
                           ("C1Q / complement","Monocyte",["C1QA","C1QB","C1QC"])]:
            ii=[gi[g] for g in gs if g in gi]
            for lab,X,mk in [("real A",XR,sixR==ct),("real T1 (depth-matched)",XT,sixT==ct),("floor",XF,sixF==ct)]:
                sub=X[mk]; tot=np.asarray(sub.sum(1)).ravel(); tot[tot==0]=1
                score=np.asarray(sp.diags(1e4/tot).dot(sub)[:,ii].sum(1)).ravel()
                o=np.sort(score)[::-1]; s=o.sum()
                conc.append(dict(module=name,celltype=ct,source=lab,n_cells=int(mk.sum()),
                                 mean_score=float(score.mean()),
                                 pct_cells_pos=float((score>0).mean()*100),
                                 share_top1pct=float(o[:max(1,len(o)//100)].sum()/s*100) if s>0 else np.nan,
                                 share_top5pct=float(o[:max(1,len(o)//20)].sum()/s*100) if s>0 else np.nan,
                                 gini=float(1-2*np.trapezoid(np.cumsum(np.sort(score))/max(score.sum(),1e-9),dx=1/len(score)))))
C=pd.DataFrame(rows); C.to_csv("out2/c2_module_rank.csv",index=False)
K=pd.DataFrame(conc); K.to_csv("out2/c2_concentration.csv",index=False)
pd.set_option("display.width",250)
print("=== module genes: expression rank inside their own cell type (real A) and detection ===")
print(C[["module","gene","rank_in_type","mean_CP10k_realA","det_realA","det_realT1","det_floor","det_ratio_floor_over_T1"]].round(3).to_string(index=False))
print("\n=== rank summary (panel has 5,961 / 5,927 measured genes) ===")
print(C.groupby("module").agg(median_rank=("rank_in_type","median"),worst_rank=("rank_in_type","max"),
                              median_det_realA=("det_realA","median"),median_det_floor=("det_floor","median"),
                              median_det_ratio=("det_ratio_floor_over_T1","median")).round(3).to_string())
print("\n=== where module UMIs sit inside the lineage ===")
print(K.round(2).to_string(index=False))
