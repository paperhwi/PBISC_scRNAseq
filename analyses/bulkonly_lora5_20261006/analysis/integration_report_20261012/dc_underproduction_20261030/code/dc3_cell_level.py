"""DC-3: label-free, cell-level. Are pDC/cDC marker UMIs produced? Are they co-expressed in the same cell (identity)?
Where do marker-positive generated cells sit (labels), and do they co-express other lineages' markers?"""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
TR="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/results/big/"
rng=np.random.default_rng(20261030)
PDC=["LILRA4","CLEC4C","PTCRA","SERPINF1","IL3RA","PLD4"]; PDC_EXCL=["LILRA4","CLEC4C","PTCRA"]
CDC=["CD1C","CLEC10A","CLEC9A","FCER1A"]
BMK=["MS4A1","CD79A","CD79B"]; MMK=["CD14","LYZ","S100A8","S100A9","VCAN"]; TMK=["CD3E","CD3D","IL32"]
rows=[]; lab=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad")
    F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad"); P=ad.read_h5ad(f"{IN}{coh}_prev4k.h5ad")
    gi={g:i for i,g in enumerate(R.var_names)}
    XF=sp.csr_matrix(F.X); sixF=F.obs.six_type.astype(str).to_numpy(); gt=np.asarray(XF.sum(1)).ravel()
    XU=sp.csr_matrix(U.X).astype(np.int64); sixU=U.obs.six_type.astype(str).to_numpy()
    pool={c:gt[sixF==c] for c in np.unique(sixF)}
    tg=np.array([rng.choice(pool.get(c,gt)) for c in sixU],float)
    kp=np.minimum(1.0,tg/np.maximum(np.asarray(XU.sum(1)).ravel(),1))
    XT=XU.copy(); XT.data=rng.binomial(XT.data,np.repeat(kp,np.diff(XT.indptr))); XT.eliminate_zeros()
    trF=pd.read_csv(f"{TR}{coh}_floor_transfer.csv.gz",index_col=0,keep_default_na=False)
    SRC={"real A":(sp.csr_matrix(R.X),R.obs.six_type.astype(str).to_numpy(),R.obs.fine.astype(str).to_numpy()),
         "real B depth-matched":(sp.csr_matrix(XT.astype(np.float32)),sixU,U.obs.fine.astype(str).to_numpy()),
         "PBISC floor":(XF,sixF,trF.fine_tr.astype(str).to_numpy()),
         "PBISC prev K4000":(sp.csr_matrix(P.X),P.obs.six_type.astype(str).to_numpy(),None)}
    def det(X,gs): ii=[gi[g] for g in gs if g in gi]; return (X[:,ii]>0).toarray()
    for s,(X,six,fine) in SRC.items():
        n=X.shape[0]
        Dp=det(X,PDC_EXCL); Dc=det(X,CDC)
        any_p=Dp.any(1); co_p=Dp.sum(1)>=2
        any_c=Dc.any(1); co_c=Dc.sum(1)>=2
        b=det(X,BMK).any(1); mo=det(X,MMK).any(1); t=det(X,TMK).any(1)
        rows.append(dict(cohort=coh,source=s,n_cells=n,
            pct_any_pDCmarker=100*any_p.mean(),pct_pDC_coexpr2=100*co_p.mean(),
            coexpr_given_any_pDC=100*co_p.sum()/max(any_p.sum(),1),
            pct_any_cDCmarker=100*any_c.mean(),pct_cDC_coexpr2=100*co_c.mean(),
            coexpr_given_any_cDC=100*co_c.sum()/max(any_c.sum(),1),
            pDCco_with_B=100*(co_p&b).sum()/max(co_p.sum(),1),
            pDCco_with_Mono=100*(co_p&mo).sum()/max(co_p.sum(),1),
            pDCco_with_T=100*(co_p&t).sum()/max(co_p.sum(),1),
            cDCco_with_B=100*(co_c&b).sum()/max(co_c.sum(),1),
            cDCco_with_Mono=100*(co_c&mo).sum()/max(co_c.sum(),1),
            cDCco_with_T=100*(co_c&t).sum()/max(co_c.sum(),1)))
        # where do single-marker pDC+ cells sit (lineage label)?
        vc=pd.Series(six[any_p]).value_counts(normalize=True)*100
        for k,v in vc.items(): lab.append(dict(cohort=coh,source=s,population="any pDC-exclusive marker",label=k,pct=v))
        vc=pd.Series(six[co_p]).value_counts(normalize=True)*100
        for k,v in vc.items(): lab.append(dict(cohort=coh,source=s,population="≥2 pDC-exclusive markers",label=k,pct=v))
T=pd.DataFrame(rows); T.to_csv("dc/dc3_cell_level.csv",index=False)
L=pd.DataFrame(lab); L.to_csv("dc/dc3_marker_cell_labels.csv",index=False)
pd.set_option("display.width",260)
print("=== % of ALL cells (label-free) ===")
print(T[["cohort","source","pct_any_pDCmarker","pct_pDC_coexpr2","coexpr_given_any_pDC","pct_any_cDCmarker","pct_cDC_coexpr2","coexpr_given_any_cDC"]].round(3).to_string(index=False))
print("\n=== lineage blending: among co-expressing DC-like cells, % that also carry another lineage's markers ===")
print(T[["cohort","source","pDCco_with_B","pDCco_with_Mono","pDCco_with_T","cDCco_with_B","cDCco_with_Mono","cDCco_with_T"]].round(1).to_string(index=False))
print("\n=== label of cells carrying pDC-exclusive marker UMIs ===")
print(L.pivot_table(index=["cohort","population","label"],columns="source",values="pct").round(1).fillna(0).to_string())
