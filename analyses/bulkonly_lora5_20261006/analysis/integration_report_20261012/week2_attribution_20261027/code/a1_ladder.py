"""A1 - control ladder + B2 three models, on the module-load metric (cell level)."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
SC="/mnt/user-data/uploads/integration_report_20261012/three_model_20261022/data/"
rng=np.random.default_rng(20261027)
SET=[("C1Q / complement","sepsis","Monocyte",["C1QA","C1QB","C1QC"]),
     ("plasma-cell secretory","sepsis","B cell",["JCHAIN","MZB1","IGLL5","DERL3","TXNDC5"]),
     ("pDC / cDC identity","sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C"]),
     ("sepsis myeloid activation","sepsis","Monocyte",["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1"]),
     ("ISG (COVID)","covid","Monocyte",["ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1"]),
     ("IFN / IFI27 (COVID)","covid","Monocyte",["IFI27","SIGLEC1"])]

def lognorm(X):
    X=sp.csr_matrix(X,copy=True).astype(np.float64); t=np.asarray(X.sum(1)).ravel(); t[t==0]=1
    X=sp.diags(1e4/t)@X; X.data=np.log1p(X.data); return sp.csr_matrix(X)

def thin(X,six,target_pool,lab_pool):
    Xb=sp.csr_matrix(X,copy=True).astype(np.int64); tot=np.asarray(Xb.sum(1)).ravel()
    pool={c:target_pool[lab_pool==c] for c in np.unique(lab_pool)}
    tgt=np.array([rng.choice(pool.get(c,target_pool)) for c in six],float)
    kp=np.minimum(1.0,tgt/np.maximum(tot,1))
    Xb.data=rng.binomial(Xb.data,np.repeat(kp,np.diff(Xb.indptr))); Xb.eliminate_zeros()
    return sp.csr_matrix(Xb.astype(np.float32)),tgt

def readout_emulate(X,six,tgt,rule):
    """kNN-15 smoothed profile on own support -> target depth -> Poisson or floor rule."""
    L=lognorm(X); P=PCA(50,random_state=0).fit_transform(np.asarray(L.todense()[:, :0]) if False else L.toarray()[:, :0]) # placeholder
    return None

def smooth_profile(X,k=15,npc=50):
    L=lognorm(X)
    P=PCA(npc,random_state=0).fit_transform(L.toarray())
    nn=NearestNeighbors(n_neighbors=k).fit(P); _,idx=nn.kneighbors(P)
    Xc=sp.csr_matrix(X); t=np.asarray(Xc.sum(1)).ravel(); t[t==0]=1
    N=sp.diags(1.0/t)@Xc                      # fractions
    N=np.asarray(N.todense())
    S=N[idx].mean(1)                          # kNN-smoothed profile
    sup=(np.asarray(Xc.todense())>0)          # own support
    S=S*sup                                   # restrict to own support
    s=S.sum(1,keepdims=True); s[s==0]=1
    return S/s

def to_counts(prof,depth,rule):
    mu=prof*depth[:,None]
    if rule=="poisson":
        return sp.csr_matrix(rng.poisson(mu).astype(np.float32))
    x=np.where(mu<1e-4,0.0,np.where(mu<1.0,1.0,np.floor(mu)))
    return sp.csr_matrix(x.astype(np.float32))

rows=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad")
    F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad"); P4=ad.read_h5ad(f"{IN}{coh}_prev4k.h5ad")
    SS=ad.read_h5ad(f"{SC}{coh}_scsemi.h5ad")
    g=np.array(R.var_names); gi={x:i for i,x in enumerate(g)}
    sixR=R.obs.six_type.astype(str).to_numpy(); sixU=U.obs.six_type.astype(str).to_numpy()
    sixF=F.obs.six_type.astype(str).to_numpy(); sixP=P4.obs.six_type.astype(str).to_numpy()
    sc_col=[c for c in ["six_type","celltype_l1_pred"] if c in SS.obs][0]
    sixS=SS.obs[sc_col].astype(str).to_numpy()
    XF=sp.csr_matrix(F.X); gt=np.asarray(XF.sum(1)).ravel()
    XU=sp.csr_matrix(U.X)
    XT1,tgt=thin(XU,sixU,gt,sixF)
    SRC={"1 real A (author)":(sp.csr_matrix(R.X),sixR),
         "2 real B (other patients)":(XU,sixU),
         "3 real B depth-matched":(XT1,sixU)}
    # read-out emulation per cell type of interest (done lineage-wise to keep PCA small)
    for name,c,ct,gs in SET:
        if c!=coh: continue
        mk=sixU==ct
        if mk.sum()<60: continue
        prof=smooth_profile(XU[mk])
        for rule,lab in [("poisson","4 real B + read-out (Poisson)"),("floor","5 real B + read-out (floor rule)")]:
            Xe=to_counts(prof,tgt[mk],rule)
            SRC.setdefault(lab,{}); SRC[lab]=SRC.get(lab,{})
            if not isinstance(SRC[lab],dict): SRC[lab]={}
            SRC[lab][ct]=Xe
    SRC["6 PBISC floor"]=(XF,sixF); SRC["7 PBISC prev K4000"]=(sp.csr_matrix(P4.X),sixP)
    SRC["8 scSemiProfiler"]=(sp.csr_matrix(SS.X),sixS)
    for name,c,ct,gs in SET:
        if c!=coh: continue
        ii=[gi[x] for x in gs if x in gi]
        for lab,v in SRC.items():
            if isinstance(v,dict):
                if ct not in v: continue
                X=v[ct]; sel=np.ones(X.shape[0],bool)
            else:
                X,six=v; sel=six==ct
                if sel.sum()<20: continue
                X=X[sel]
            t=np.asarray(X.sum(1)).ravel(); t[t==0]=1
            cp=np.asarray(sp.diags(1e4/t).dot(X)[:,ii].sum(1)).ravel()
            pos=cp>0
            rows.append(dict(module=name,cohort=coh,celltype=ct,source=lab,n=int(X.shape[0]),
                             pct_pos=100*pos.mean(),load_in_pos=float(cp[pos].mean()) if pos.any() else 0.0,
                             mean_all=float(cp.mean())))
L=pd.DataFrame(rows)
base=L[L.source=="3 real B depth-matched"].set_index(["module"])[["pct_pos","load_in_pos"]]
L["pos_ratio"]=[r.pct_pos/base.loc[r.module,"pct_pos"] for r in L.itertuples()]
L["load_ratio"]=[r.load_in_pos/base.loc[r.module,"load_in_pos"] for r in L.itertuples()]
L.to_csv("out2/a1_ladder.csv",index=False)
pd.set_option("display.width",250)
print("=== A1 control ladder: module load, ratio vs depth-matched real B ===")
print(L.pivot_table(index="source",columns="module",values="load_ratio").round(3).to_string())
print("\n=== % positive cells, same reference ===")
print(L.pivot_table(index="source",columns="module",values="pos_ratio").round(3).to_string())
