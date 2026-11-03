"""Item 7 - separate the immediate-early / handling programme from the disease contrast."""
import numpy as np, pandas as pd, pickle
from scipy import stats
from deg_lib import *
allde=pickle.load(open("out/allde.pkl","rb"))

# core immediate-early / stress-response genes (transcription factors + canonical primary response)
IEG_CORE=["FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","EGR4","NR4A1","NR4A2","NR4A3",
 "ATF3","IER2","IER3","IER5","DUSP1","DUSP2","DUSP5","KLF2","KLF4","KLF6","KLF10","ZFP36","ZFP36L1","ZFP36L2",
 "CSRNP1","RGS1","RGS2","RGS16","PPP1R15A","BTG1","BTG2","MYC","NFKBIA","NFKBIZ","TNFAIP3","CD69","CD83",
 "RASD1","PHLDA1","TRIB1","PER1","HES1","MAFF","SGK1","SOCS3","SLC7A5","HSPA1A","HSPA1B","DNAJB1","HSPH1",
 "ID1","ID2","ID3","CCNL1","MIDN","SIK1","PLK2","GADD45B","CEBPD","ARC","NPAS4"]
AMBIENT=["HBB","HBA1","HBA2","HBD","PPBP","PF4","NRGN","ALAS2"]   # RBC / platelet carry-over

rows=[]; detail=[]
for (coh,ct),t in allde.items():
    R=t["real"].set_index("gene"); F=t["floor"].set_index("gene").loc[R.index]
    P=t["prev"].set_index("gene").loc[R.index] if "prev" in t else None
    ieg=R.index.isin(IEG_CORE); amb=R.index.isin(AMBIENT)
    drop=ieg|amb
    def m(sel):
        a=R.logFC[sel].values; b=F.logFC[sel].values
        rd=R.deg[sel].values; fd=F.deg[sel].values
        out=dict(n_genes=int(sel.sum()),logFC_r=float(np.corrcoef(a,b)[0,1]),
                 slope=float(np.polyfit(a,b,1)[0]),
                 n_DEG_real=int(rd.sum()),n_DEG_gen=int(fd.sum()),
                 recall=float((rd&fd).sum()/max(rd.sum(),1)) if rd.sum() else np.nan,
                 precision=float((rd&fd).sum()/max(fd.sum(),1)) if fd.sum() else np.nan)
        return out
    full=m(np.ones(len(R),bool)); kept=m(~drop)
    ieg_only=m(ieg) if ieg.sum()>=8 else {k:np.nan for k in full}
    rows.append(dict(cohort=coh,celltype=ct,
        r_all=full["logFC_r"],r_noIEG=kept["logFC_r"],delta_r=kept["logFC_r"]-full["logFC_r"],
        r_IEGonly=ieg_only["logFC_r"],
        recall_all=full["recall"],recall_noIEG=kept["recall"],
        prec_all=full["precision"],prec_noIEG=kept["precision"],
        nDEG_real_all=full["n_DEG_real"],nDEG_real_noIEG=kept["n_DEG_real"],
        pct_realDEG_IEG=100*(R.deg&ieg).sum()/max(R.deg.sum(),1),
        pct_realDEG_ambient=100*(R.deg&amb).sum()/max(R.deg.sum(),1),
        n_IEG_in_panel=int(ieg.sum())))
    # how much of the real signal is IEG, by |logFC| rank
    top=R[R.q<0.05].reindex(R[R.q<0.05].logFC.abs().sort_values(ascending=False).index).head(30)
    detail.append(dict(cohort=coh,celltype=ct,
        pct_top30_IEG=100*top.index.isin(IEG_CORE).mean(),
        top30_IEG=",".join([g for g in top.index if g in IEG_CORE][:10])))
S=pd.DataFrame(rows); DT=pd.DataFrame(detail)
S=S.merge(DT,on=["cohort","celltype"])
S.to_csv("out2/ieg_split.csv",index=False)
pd.set_option("display.width",260)
print("=== disease logFC r before / after removing immediate-early + ambient genes ===")
print(S[["cohort","celltype","n_IEG_in_panel","pct_realDEG_IEG","pct_top30_IEG","r_all","r_noIEG","delta_r","r_IEGonly",
         "recall_all","recall_noIEG","nDEG_real_all","nDEG_real_noIEG"]].round(3).to_string(index=False))
print("\n=== which IEGs dominate the real top-30 ===")
print(DT.to_string(index=False))

# ---- is the IEG programme itself disease-associated in the real data? ----
rows=[]
for coh in ["sepsis","covid"]:
    d=load(f"real_{coh}"); des=design(f"real_{coh}",d); genes=list(d["genes"])
    df=load(f"floor_{coh}"); desf=design(f"floor_{coh}",df)
    for ct in SIX:
        for src,(dd,de) in [("real",(d,des)),("floor",(df,desf))]:
            g=cpm_matrix(dd,de,ct)
            if len(g)==0: continue
            cols=[x for x in genes if x in IEG_CORE]
            sc=g[cols].mean(1)
            a=sc[g._c=="case"]; b=sc[g._c=="ctrl"]
            if len(a)<3 or len(b)<3: continue
            t,p=stats.ttest_ind(a,b,equal_var=False)
            rows.append(dict(cohort=coh,celltype=ct,source=src,n_case=len(a),n_ctrl=len(b),
                             ieg_score_case=a.mean(),ieg_score_ctrl=b.mean(),
                             delta=a.mean()-b.mean(),p=float(p),
                             cohens_d=float((a.mean()-b.mean())/np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2))))
Q=pd.DataFrame(rows); Q.to_csv("out2/ieg_score_contrast.csv",index=False)
print("\n=== IEG programme score, case vs control (is the stress programme itself disease-linked?) ===")
print(Q.round(3).to_string(index=False))
