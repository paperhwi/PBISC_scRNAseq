import numpy as np, pandas as pd, mygene, os, pickle
PB="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/pb/"
BK="/mnt/user-data/uploads/integration_report_20261012/disease_deg_rare_20261027/bulk_cache/"
ZN="/mnt/user-data/uploads/bulk_public/zenodo6497251_imperial_eskd/"
genes=np.load(PB+"real_sepsis.npz",allow_pickle=True)["genes"].astype(str)
mg=mygene.MyGeneInfo()
r=mg.querymany(list(genes),scopes="symbol",fields="ensembl.gene",species="human",verbose=False)
m={}
for x in r:
    if "ensembl" not in x: continue
    e=x["ensembl"]; ids=[e["gene"]] if isinstance(e,dict) else [y["gene"] for y in e]
    m.setdefault(x["query"],[]).extend(ids)
print("mapped symbols:",len(m),"/",len(genes))
pickle.dump(m,open("out2/sym2ens.pkl","wb"))

# ---- sepsis bulk ----
S=pd.read_csv(BK+"sepsis_bulk_counts.csv.gz",index_col=0)
S.index=[i.split(".")[0] for i in S.index]
print("sepsis bulk",S.shape)
ro=pd.read_csv(PB+"real_sepsis_obs.csv"); fo=pd.read_csv(PB+"floor_sepsis_obs.csv")
gmap={}
for _,x in ro.iterrows():
    gmap[x["sample"]]=("case" if x.group=="sepsis" else "ctrl" if x.group=="healthy control" else "icu" if x.group=="ICU control" else None, x.age_group)
for _,x in fo.iterrows():
    g=x.group
    gmap.setdefault(x["sample"],("case" if g=="adult_sepsis" else "ctrl" if g=="adult_HC" else "icu" if "ICU" in g else None,
                                "paediatric" if g.startswith("ped") else "adult"))
cov=[s for s in S.columns if s in gmap]
print("bulk samples with group:",len(cov),"/",S.shape[1])
meta=pd.DataFrame([(s,)+gmap[s] for s in cov],columns=["sample","cond","age"]).set_index("sample")
print(meta.groupby(["age","cond"]).size())
S[cov].to_pickle("out2/sepsis_bulk.pkl"); meta.to_csv("out2/sepsis_bulk_meta.csv")

# ---- covid bulk ----
C=pd.read_csv(ZN+"htseq_counts.csv",index_col=0)
C.index=[i.split(".")[0] for i in C.index]
w=pd.concat([pd.read_csv(ZN+"w1_metadata.csv"),pd.read_csv(ZN+"w2_metadata.csv")],ignore_index=True)
w=w.drop_duplicates("sample_id").set_index("sample_id")
cc=[s for s in C.columns if s in w.index and w.loc[s,"case_control"] in ("POSITIVE","NEGATIVE")]
cm=pd.DataFrame({"cond":["case" if w.loc[s,"case_control"]=="POSITIVE" else "ctrl" for s in cc],
                 "patient":[w.loc[s,"individual_id"] for s in cc]},index=cc)
print("covid bulk",C.shape,"used",len(cc)); print(cm.cond.value_counts())
C[cc].to_pickle("out2/covid_bulk.pkl"); cm.to_csv("out2/covid_bulk_meta.csv")
