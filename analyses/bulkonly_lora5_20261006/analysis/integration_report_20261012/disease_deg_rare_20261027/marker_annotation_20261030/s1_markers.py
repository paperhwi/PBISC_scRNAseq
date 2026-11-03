# Step 1: hierarchical Top10 markers from real A (author labels)
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, pickle, warnings; warnings.filterwarnings("ignore")
B="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
MINC=20
def top10(a,key,groups,ref="rest"):
    sc.tl.rank_genes_groups(a,key,groups=groups,reference=ref,method="wilcoxon",pts=True,use_raw=False)
    out={}
    for g in groups:
        d=sc.get.rank_genes_groups_df(a,g)
        d=d[~d.names.str.match(r"^(RP[SL]|MT-)")]  # drop ribosomal / mito genes
        ok=d[(d.logfoldchanges>0.5)&(d.pvals_adj<0.05)]
        if "pct_nz_group" in d: ok=ok[ok.pct_nz_group>ok.pct_nz_reference]
        genes=ok.names.tolist()[:10]
        if len(genes)<10: genes+= [x for x in d.names if x not in genes][:10-len(genes)]
        out[g]=genes
    return out
M={}
for coh in ["sepsis","covid"]:
    a=ad.read_h5ad(B+f"{coh}_ref.h5ad"); a.obs["six_type"]=a.obs.six_type.astype(str); a.obs["fine"]=a.obs.fine.astype(str)
    sc.pp.normalize_total(a,1e4); sc.pp.log1p(a)
    six=[s for s,n in a.obs.six_type.value_counts().items() if n>=MINC]
    L1=top10(a,"six_type",six)
    L2={}
    for s in six:
        b=a[a.obs.six_type==s].copy()
        fines=[f for f,n in b.obs.fine.value_counts().items() if n>=MINC]
        if len(fines)<2: L2[s]={f:L1[s] for f in fines}; continue
        L2[s]=top10(b,"fine",fines)
    M[coh]=dict(L1=L1,L2=L2)
    print(coh,{s:len(v) for s,v in L2.items()})
    for s in L2:
        for f,g in L2[s].items(): print("  ",s,"|",f,":",",".join(g))
pickle.dump(M,open("markerann/markers.pkl","wb"))
rows=[dict(cohort=c,level="six_type",six_type=s,subset=s,markers=",".join(g)) for c in M for s,g in M[c]["L1"].items()]
rows+=[dict(cohort=c,level="fine",six_type=s,subset=f,markers=",".join(g)) for c in M for s in M[c]["L2"] for f,g in M[c]["L2"][s].items()]
pd.DataFrame(rows).to_csv("markerann/top10_markers_realA.csv",index=False)
