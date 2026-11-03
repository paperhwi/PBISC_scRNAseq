"""Matched-n DEG: real and floor subsampled to the same case / control unit counts, 50 random draws.
Design M: n_case = min(real, floor), n_ctrl = min(real, floor).  Design B: n_case = n_ctrl = the smaller control arm."""
import numpy as np, pandas as pd, pickle
from deg_lib import *
rng=np.random.default_rng(20261030)
IEG=set(["FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","EGR4","NR4A1","NR4A2","NR4A3","ATF3","IER2","IER3","IER5",
 "DUSP1","DUSP2","DUSP5","KLF2","KLF4","KLF6","KLF10","ZFP36","ZFP36L1","ZFP36L2","CSRNP1","RGS1","RGS2","RGS16","PPP1R15A",
 "BTG1","BTG2","MYC","NFKBIA","NFKBIZ","TNFAIP3","CD69","CD83","RASD1","PHLDA1","TRIB1","PER1","HES1","MAFF","SGK1","SOCS3",
 "SLC7A5","HSPA1A","HSPA1B","DNAJB1","HSPH1","ID1","ID2","ID3","CCNL1","MIDN","SIK1","PLK2","GADD45B","CEBPD","ARC","NPAS4",
 "HBB","HBA1","HBA2","HBD","PPBP","PF4","NRGN","ALAS2"])
CT=["T cell","B cell","NK cell"]; REPS=50
rows=[]; rep_tabs={}
for coh in ["sepsis","covid"]:
    srcs={s:load(f"{s}_{coh}") for s in ["real","floor"]}; genes=list(srcs["real"]["genes"])
    meas=set(srcs["real"]["genes"][(srcs["real"]["S"].sum(0)>0)&(srcs["floor"]["S"].sum(0)>0)])
    G={s:{ct:cpm_matrix(srcs[s],design(f"{s}_{coh}",srcs[s]),ct) for ct in CT} for s in srcs}
    for ct in CT:
        units={s:(G[s][ct][G[s][ct]._c=="case"]._u.unique(),G[s][ct][G[s][ct]._c=="ctrl"]._u.unique()) for s in G}
        nc=min(len(units["real"][0]),len(units["floor"][0])); no=min(len(units["real"][1]),len(units["floor"][1]))
        for dname,(ncase,nctrl) in {"M":(nc,no),"B":(min(nc,no),min(nc,no))}.items():
            store=[]
            for rep in range(REPS):
                res={}
                for s in ["real","floor"]:
                    ca=rng.choice(units[s][0],ncase,replace=False); co=rng.choice(units[s][1],nctrl,replace=False)
                    g=G[s][ct]; sub=g[((g._c=="case")&g._u.isin(ca))|((g._c=="ctrl")&g._u.isin(co))]
                    t=deg(sub,genes); t=t[t.gene.isin(meas)].reset_index(drop=True); res[s]=t
                R=res["real"]; F=res["floor"].set_index("gene").loc[R.gene].reset_index()
                rd=R.deg.values; fd=F.deg.values
                r=dict(cohort=coh,celltype=ct,design=dname,n_case=ncase,n_ctrl=nctrl,rep=rep,
                       real_up=int((rd&(R.logFC>0)).sum()),real_down=int((rd&(R.logFC<0)).sum()),
                       floor_up=int((fd&(F.logFC>0)).sum()),floor_down=int((fd&(F.logFC<0)).sum()),
                       real_IEG=int((rd&R.gene.isin(IEG)).sum()),floor_IEG=int((fd&F.gene.isin(IEG)).sum()),
                       logFC_r=float(np.corrcoef(R.logFC,F.logFC)[0,1]),
                       recall=float((rd&fd).sum()/rd.sum()) if rd.sum() else np.nan,
                       precision=float((rd&fd).sum()/fd.sum()) if fd.sum() else np.nan,
                       real_maxlogq=float(-np.log10(R.q.clip(1e-300)).max()),floor_maxlogq=float(-np.log10(F.q.clip(1e-300)).max()))
                rows.append(r); store.append((r,res))
            # representative draw: total DEG counts closest to the median of both sources
            d=pd.DataFrame([x[0] for x in store]); tr=d.real_up+d.real_down; tf=d.floor_up+d.floor_down
            k=int(((tr-tr.median()).abs()/max(tr.std(),1)+(tf-tf.median()).abs()/max(tf.std(),1)).idxmin())
            rep_tabs[(coh,ct,dname)]=store[k][1]
D=pd.DataFrame(rows); D.to_csv("matched/matched_n_reps.csv",index=False)
pickle.dump(rep_tabs,open("matched/representative_tables.pkl","wb"))
S=D.groupby(["cohort","celltype","design","n_case","n_ctrl"]).agg(
    real_up=("real_up","mean"),real_up_sd=("real_up","std"),real_down=("real_down","mean"),real_down_sd=("real_down","std"),
    floor_up=("floor_up","mean"),floor_up_sd=("floor_up","std"),floor_down=("floor_down","mean"),floor_down_sd=("floor_down","std"),
    real_IEG=("real_IEG","mean"),floor_IEG=("floor_IEG","mean"),
    logFC_r=("logFC_r","mean"),logFC_r_sd=("logFC_r","std"),recall=("recall","mean"),precision=("precision","mean"),
    real_maxlogq=("real_maxlogq","median"),floor_maxlogq=("floor_maxlogq","median")).reset_index()
S.to_csv("matched/matched_n_summary.csv",index=False)
pd.set_option("display.width",260)
print(S.round(2).to_string(index=False))
