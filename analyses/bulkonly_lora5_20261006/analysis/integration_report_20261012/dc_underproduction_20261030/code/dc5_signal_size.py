"""DC-5: how large is the DC signal inside the bulk input, compared with other lineages' markers?"""
import numpy as np, pandas as pd, pickle
m=pickle.load(open("out2/sym2ens.pkl","rb"))
G={"pDC":["LILRA4","CLEC4C","PTCRA"],"cDC":["CD1C","CLEC10A","CLEC9A"],
   "monocyte":["CD14","LYZ","VCAN"],"B cell":["MS4A1","CD79A","CD19"],"T cell":["CD3E","CD3D","CD2"],"NK":["GNLY","KLRF1","NCAM1"]}
rows=[]
for coh in ["sepsis","covid"]:
    B=pd.read_pickle(f"out2/{coh}_bulk.pkl"); C=B/B.sum(0)*1e6
    for k,gs in G.items():
        ids=[max([i for i in m.get(g,[]) if i in C.index],key=lambda i:C.loc[i].mean()) for g in gs if any(i in C.index for i in m.get(g,[]))]
        v=C.loc[ids].sum(0)
        rows.append(dict(cohort=coh,lineage=k,marker_CPM_median=v.median(),
                         share_of_bulk_pct=100*v.median()/1e6,
                         cv_across_samples=v.std()/v.mean()))
S=pd.DataFrame(rows); S.to_csv("dc/dc5_bulk_signal_size.csv",index=False)
pd.set_option("display.width",200)
p=S.pivot_table(index="lineage",columns="cohort",values=["marker_CPM_median","cv_across_samples"]).round(2)
print(p.loc[["T cell","monocyte","NK","B cell","cDC","pDC"]].to_string())
