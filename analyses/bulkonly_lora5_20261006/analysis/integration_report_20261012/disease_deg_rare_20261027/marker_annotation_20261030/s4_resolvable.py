# Step 4: per-subset marker resolvability on real B (F1 vs author) and key-subset table
import numpy as np, pandas as pd
from sklearn.metrics import precision_recall_fscore_support
P=pd.read_csv("markerann/proportions_long.csv"); E=pd.read_csv("markerann/expression_similarity_long.csv")
R=[]
for coh in ["sepsis","covid"]:
    L=pd.read_csv(f"markerann/{coh}_realB_thin_labels.csv.gz",index_col=0)
    au=L.fine.astype(str).values
    for col in ["fine_mk_cl","fine_mk_cell"]:
        pr=L[col].fillna("").astype(str).values; labs=sorted(set(P[P.cohort==coh].subset))
        p,r,f,_=precision_recall_fscore_support(au,pr,labels=labs,zero_division=0)
        R+= [dict(cohort=coh,label=col,subset=s,precision=a,recall=b,F1=c) for s,a,b,c in zip(labs,p,r,f)]
R=pd.DataFrame(R); R.to_csv("markerann/realB_per_subset_F1.csv",index=False)
F=R[R.label=="fine_mk_cl"].set_index(["cohort","subset"]).F1
print("resolvable subsets (real B F1>=0.5, cluster):",(F>=0.5).groupby(level=0).sum().to_dict(),"of",F.groupby(level=0).size().to_dict())
P["F1_realB"]=[F.get((c,s),np.nan) for c,s in zip(P.cohort,P.subset)]
E["F1_realB"]=[F.get((c,s),np.nan) for c,s in zip(E.cohort,E.subset)]
P["l2r"]=np.log2((P.pct+0.01)/(P.pct_realA+0.01))
res=P[P.F1_realB>=0.5]
print(res.groupby(["cohort","combo"]).l2r.agg(med_abs=lambda x: np.median(np.abs(x)),n_lost=lambda x:(x<-2).sum(),n="size").round(2).to_string())
er=E[E.F1_realB>=0.5]
print(er.groupby(["cohort","combo"])[["r_identity_nonmarker","marker_contrast_ratio"]].median().round(2).to_string())
KEY={"sepsis":["Macrophage","Mono_c06-CD14-MKI67","pDC","cDC","B_c07-MZB1","B_c08-MKI67","T_CD4_c07-MKI67","Megakaryocyte"],
     "covid":["CD14mono_dex","CD16mono_C1","pDC","DC1","DC2","B_ASC_dividing","B_ASC_IgA","MAIT","Platelet"]}
rows=[]
for coh,ks in KEY.items():
    for s in ks:
        q=P[(P.cohort==coh)&(P.subset==s)].set_index("combo")
        e=E[(E.cohort==coh)&(E.subset==s)].set_index("combo")
        d=dict(cohort=coh,subset=s,F1_realB=round(F.get((coh,s),np.nan),2),pct_realA=round(q.pct_realA.iloc[0],2))
        for c,short in [("real B · author","realB_author"),("real B · marker (cluster)","realB_mk"),("floor · kNN transfer","floor_knn"),("floor · marker (cluster)","floor_mk")]:
            d["pct_"+short]=round(q.pct.get(c,np.nan),2)
            d["rID_"+short]=round(e.r_identity_nonmarker.get(c,np.nan),2) if c in e.index else np.nan
            d["mkC_"+short]=round(e.marker_contrast_ratio.get(c,np.nan),2) if c in e.index else np.nan
        rows.append(d)
K=pd.DataFrame(rows); K.to_csv("markerann/key_subsets_table.csv",index=False)
pd.set_option("display.width",250); print(K.to_string())
