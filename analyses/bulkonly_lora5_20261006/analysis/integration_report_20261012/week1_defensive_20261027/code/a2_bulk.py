import numpy as np, pandas as pd, pickle
from scipy import stats
from deg_lib import bh
m=pickle.load(open("out2/sym2ens.pkl","rb"))
allde=pickle.load(open("out/allde.pkl","rb"))

def bulk_deg(X,meta,unit=None):
    tot=X.sum(0); L=np.log2(X/tot*1e6+1)
    L=L.T; L["_c"]=meta.loc[L.index,"cond"].values
    if unit is not None:
        L["_u"]=meta.loc[L.index,unit].values; L=L.groupby(["_u","_c"]).mean().reset_index()
    a=L[L._c=="case"].drop(columns=[c for c in L.columns if c.startswith("_")]).to_numpy()
    b=L[L._c=="ctrl"].drop(columns=[c for c in L.columns if c.startswith("_")]).to_numpy()
    gl=[c for c in L.columns if not c.startswith("_")]
    t,p=stats.ttest_ind(a,b,equal_var=False)
    return pd.DataFrame(dict(ensembl=gl,bulk_logFC=a.mean(0)-b.mean(0),bulk_p=p,bulk_q=bh(p),
                             bulk_mean=(a.mean(0)+b.mean(0))/2,n_case=len(a),n_ctrl=len(b)))

S=pd.read_pickle("out2/sepsis_bulk.pkl"); sm=pd.read_csv("out2/sepsis_bulk_meta.csv",index_col=0)
sa=sm[(sm.age=="adult")&(sm.cond.isin(["case","ctrl"]))]
BS=bulk_deg(S[sa.index],sa)
C=pd.read_pickle("out2/covid_bulk.pkl"); cm=pd.read_csv("out2/covid_bulk_meta.csv",index_col=0)
BC=bulk_deg(C,cm,unit="patient")
print("sepsis bulk n:",BS.n_case[0],"/",BS.n_ctrl[0],"  covid bulk patients:",BC.n_case[0],"/",BC.n_ctrl[0])

# symbol-level bulk table (sum over mapped ensembl ids -> take the most expressed)
def to_symbol(B):
    rows=[]
    bi=B.drop_duplicates("ensembl").set_index("ensembl")
    for s,ids in m.items():
        hit=[i for i in dict.fromkeys(ids) if i in bi.index]
        if not hit: continue
        sub=bi.loc[hit]
        r=sub.iloc[int(np.argmax(sub.bulk_mean.to_numpy()))]
        rows.append((s,float(r.bulk_logFC),float(r.bulk_q),float(r.bulk_mean)))
    return pd.DataFrame(rows,columns=["gene","bulk_logFC","bulk_q","bulk_mean"]).set_index("gene")
bs=to_symbol(BS); bc=to_symbol(BC)
bs["bulk_rank"]=stats.rankdata(-bs.bulk_mean.to_numpy())
bc["bulk_rank"]=stats.rankdata(-bc.bulk_mean.to_numpy())
bs.to_csv("out2/sepsis_bulk_deg.csv"); bc.to_csv("out2/covid_bulk_deg.csv")
print("symbols with bulk values:",len(bs),len(bc))

MOD={"C1Q / complement":("sepsis","Monocyte",["C1QA","C1QB","C1QC","RNASE1","CLU"]),
     "plasma-cell secretory":("sepsis","B cell",["JCHAIN","MZB1","IGLL5","DERL3","TNFRSF17","TXNDC5"]),
     "pDC / cDC identity":("sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C"]),
     "sepsis myeloid activation":("sepsis","Monocyte",["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1"]),
     "MHC-II":("sepsis","Monocyte",["HLA-DRA","HLA-DRB1","HLA-DPA1","CD74"]),
     "ISG (COVID)":("covid","Monocyte",["ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1"]),
     "IFI27 / SIGLEC1":("covid","Monocyte",["IFI27","SIGLEC1"])}
rows=[]
for name,(coh,ct,gs) in MOD.items():
    B=bs if coh=="sepsis" else bc
    R=allde[(coh,ct)]["real"].set_index("gene"); F=allde[(coh,ct)]["floor"].set_index("gene")
    for g in gs:
        if g not in B.index or g not in R.index:
            rows.append(dict(module=name,cohort=coh,celltype=ct,gene=g,note="missing")); continue
        rows.append(dict(module=name,cohort=coh,celltype=ct,gene=g,
            bulk_logFC=B.loc[g,"bulk_logFC"],bulk_q=B.loc[g,"bulk_q"],
            bulk_rank=int(B.loc[g,"bulk_rank"]),
            real_logFC=R.loc[g,"logFC"],real_q=R.loc[g,"q"],floor_logFC=F.loc[g,"logFC"],floor_q=F.loc[g,"q"]))
M=pd.DataFrame(rows); M["bulk_sig"]=M.bulk_q<0.05; M["bulk_same_sign"]=np.sign(M.bulk_logFC)==np.sign(M.real_logFC)
M.to_csv("out2/a2_module_bulk.csv",index=False)
pd.set_option("display.width",250)
print("\n=== per-gene: bulk vs real(sc) vs floor ===")
print(M[["module","gene","bulk_logFC","bulk_q","bulk_sig","real_logFC","floor_logFC"]].round(3).to_string(index=False))
print("\n=== per-module summary ===")
g=M.dropna(subset=["bulk_logFC"]).groupby("module")
sumtab=pd.DataFrame({"n_genes":g.size(),
  "bulk_sig_n":g.bulk_sig.sum(),"bulk_median_logFC":g.bulk_logFC.median(),
  "bulk_sign_match":g.bulk_same_sign.mean(),
  "real_median_logFC":g.real_logFC.median(),"floor_median_logFC":g.floor_logFC.median(),
  "floor/real":g.apply(lambda d:(d.floor_logFC/d.real_logFC).median()),
  "bulk/real":g.apply(lambda d:(d.bulk_logFC/d.real_logFC).median())})
print(sumtab.round(3).to_string())

# background: what fraction of real cell-type DEGs are bulk-detectable?
bg=[]
for (coh,ct),t in allde.items():
    B=bs if coh=="sepsis" else bc
    R=t["real"].set_index("gene"); d=R[R.deg]
    common=[x for x in d.index if x in B.index]
    if len(common)<10: continue
    bb=B.loc[common]
    bg.append(dict(cohort=coh,celltype=ct,n_realDEG=len(common),
                   bulk_sig_frac=float((bb.bulk_q<0.05).mean()),
                   bulk_sign_frac=float((np.sign(bb.bulk_logFC)==np.sign(d.loc[common,"logFC"])).mean()),
                   bulk_logFC_r=float(np.corrcoef(bb.bulk_logFC,d.loc[common,"logFC"])[0,1])))
BG=pd.DataFrame(bg); BG.to_csv("out2/a2_bulk_background.csv",index=False)
print("\n=== background: real cell-type DEGs recovered in bulk ===")
print(BG.round(3).to_string(index=False))
