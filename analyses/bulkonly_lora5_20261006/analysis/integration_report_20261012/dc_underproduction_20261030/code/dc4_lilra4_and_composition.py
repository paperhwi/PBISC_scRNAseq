"""DC-4: (a) LILRA4+ cells: are they assembled pDCs (co-express the pDC programme)? (b) per-sample DC % vs bulk DC signal."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, pickle
from scipy import stats
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
PB="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/pb/"
rng=np.random.default_rng(20261030)
PROG=["CLEC4C","IRF7","TCF4","SPIB","IRF8","GZMB","PLD4","SERPINF1","IL3RA","ITM2C","PTCRA"]
rows=[]
for coh in ["sepsis","covid"]:
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad")
    F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad"); P=ad.read_h5ad(f"{IN}{coh}_prev4k.h5ad")
    gi={g:i for i,g in enumerate(R.var_names)}
    XF=sp.csr_matrix(F.X); sixF=F.obs.six_type.astype(str).to_numpy(); gt=np.asarray(XF.sum(1)).ravel()
    XU=sp.csr_matrix(U.X).astype(np.int64); sixU=U.obs.six_type.astype(str).to_numpy()
    pool={c:gt[sixF==c] for c in np.unique(sixF)}; tg=np.array([rng.choice(pool.get(c,gt)) for c in sixU],float)
    kp=np.minimum(1.0,tg/np.maximum(np.asarray(XU.sum(1)).ravel(),1)); XT=XU.copy()
    XT.data=rng.binomial(XT.data,np.repeat(kp,np.diff(XT.indptr))); XT.eliminate_zeros()
    SRC={"real A":(sp.csr_matrix(R.X),R.obs.six_type.astype(str).to_numpy()),
         "real B depth-matched":(sp.csr_matrix(XT.astype(np.float32)),sixU),
         "PBISC floor":(XF,sixF),"PBISC prev K4000":(sp.csr_matrix(P.X),P.obs.six_type.astype(str).to_numpy())}
    for s,(X,six) in SRC.items():
        L=(X[:,gi["LILRA4"]]>0).toarray().ravel()
        prog=(X[:,[gi[g] for g in PROG if g in gi]]>0).toarray()
        nprog=prog.sum(1)
        umi=X[:,gi["LILRA4"]].toarray().ravel()
        rows.append(dict(cohort=coh,source=s,
            pct_LILRA4pos=100*L.mean(),
            LILRA4_UMI_in_pos=float(umi[L].mean()) if L.any() else np.nan,
            pDCprog_genes_in_LILRA4pos=float(nprog[L].mean()) if L.any() else np.nan,
            pDCprog_genes_in_LILRA4neg=float(nprog[~L].mean()),
            pct_LILRA4pos_labelled_DC=100*(six[L]=="DC").mean() if L.any() else np.nan,
            pct_assembled_pDC=100*((L)&(nprog>=4)).mean()))
A=pd.DataFrame(rows); A.to_csv("dc/dc4_lilra4.csv",index=False)
pd.set_option("display.width",240)
print("=== LILRA4 (most pDC-exclusive gene): produced? assembled into a pDC? ===")
print(A.round(3).to_string(index=False))

# ---- (b) composition tracking: per-sample DC % vs bulk DC marker score ----
m=pickle.load(open("out2/sym2ens.pkl","rb"))
def bscore(L,gs):
    ids=[max([i for i in m.get(g,[]) if i in L.index],key=lambda i:L.loc[i].mean()) for g in gs if any(i in L.index for i in m.get(g,[]))]
    return L.loc[ids].mean(0)
out=[]
for coh in ["sepsis"]:
    B=pd.read_pickle(f"out2/{coh}_bulk.pkl"); Lb=np.log2(B/B.sum(0)*1e6+1)
    xs=bscore(Lb,["LILRA4","CLEC4C","PTCRA","CD1C","CLEC10A","CLEC9A"])
    for src in ["real","floor","prev"]:
        z=np.load(f"{PB}{src}_{coh}.npz",allow_pickle=True); gr=pd.Series(z["groups"].astype(str))
        df=pd.DataFrame({"s":gr.str.split("\\|\\|").str[0],"six":gr.str.split("\\|\\|").str[1],"n":z["n"]})
        w=df.pivot_table(index="s",columns="six",values="n",aggfunc="sum",fill_value=0)
        dc=100*w["DC"]/w.sum(1); c=[s for s in dc.index if s in xs.index]
        r=stats.spearmanr(xs[c],dc[c])
        out.append(dict(cohort=coh,source=src,n=len(c),spearman_bulkDCscore_vs_DCpct=r.statistic,p=r.pvalue,
                        DCpct_range_p10_p90=f"{np.percentile(dc[c],10):.2f}–{np.percentile(dc[c],90):.2f}"))
O=pd.DataFrame(out); O.to_csv("dc/dc4_composition_tracking.csv",index=False)
print("\n=== per-sample DC % vs bulk DC marker score (sepsis, same samples) ===")
print(O.round(3).to_string(index=False))
