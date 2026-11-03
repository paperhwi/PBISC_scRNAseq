# Step 3: proportions + expression similarity vs real A, kNN vs marker labels
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, pickle, warnings; warnings.filterwarnings("ignore")
from sklearn.metrics import adjusted_rand_score, f1_score
B="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
M=pickle.load(open("markerann/markers.pkl","rb"))
MINC=20
def means(X,lab):
    lab=pd.Series(lab); out={}
    for g,idx in lab.groupby(lab).groups.items():
        if len(idx)>=MINC and g!="": out[g]=np.asarray(X[np.asarray(idx)].mean(0)).ravel()
    return out
P=[];E=[];V=[]
for coh in ["sepsis","covid"]:
    A=ad.read_h5ad(B+f"{coh}_ref.h5ad"); sc.pp.normalize_total(A,1e4); sc.pp.log1p(A)
    fineA=A.obs.fine.astype(str).values; sixA=A.obs.six_type.astype(str).values
    f2s=pd.Series(sixA,index=fineA).groupby(level=0).first().to_dict()
    allmk=set(g for d in M[coh]["L2"].values() for v in d.values() for g in v)|set(g for v in M[coh]["L1"].values() for g in v)
    H=A.copy(); sc.pp.highly_variable_genes(H,n_top_genes=2000,flavor="seurat")
    hv=np.array([g in set(H.var_names[H.var.highly_variable]) and g not in allmk for g in A.var_names]); del H
    genes=np.array(A.var_names); gi={g:i for i,g in enumerate(genes)}
    mA=means(A.X,fineA); muA=np.asarray(A.X.mean(0)).ravel()
    pA=pd.Series(fineA).value_counts(normalize=True)
    fmk={f:v for d in M[coh]["L2"].values() for f,v in d.items()}
    combos=[("floor","fine_knn","floor · kNN transfer"),("floor","fine_mk_cl","floor · marker (cluster)"),("floor","fine_mk_cell","floor · marker (cell)"),
            ("realB_thin","fine","real B · author"),("realB_thin","fine_mk_cl","real B · marker (cluster)"),("realB_thin","fine_mk_cell","real B · marker (cell)"),
            ("realB_full","fine","real B full · author"),("realB_full","fine_mk_cl","real B full · marker (cluster)")]
    cache={}
    for ds,col,lab in combos:
        L=pd.read_csv(f"markerann/{coh}_{ds}_labels.csv.gz",index_col=0)[col].fillna("").astype(str).values
        pD=pd.Series(L).value_counts(normalize=True)
        for s in pA.index:
            P.append(dict(cohort=coh,combo=lab,dataset=ds,label=col,subset=s,six=f2s.get(s),pct_realA=100*pA[s],pct=100*pD.get(s,0.0),n=int((L==s).sum())))
        if ds=="realB_full": continue
        if ds not in cache:
            D=ad.read_h5ad(f"markerann/{coh}_{ds}_annot.h5ad"); cache[ds]=(D.X, np.asarray(D.X.mean(0)).ravel(), D.obs)
        X,muD,obs=cache[ds]; mD=means(X,L)
        for s in mA:
            if s not in mD: continue
            dA=mA[s]-muA; dD=mD[s]-muD
            r_non=np.corrcoef(dA[hv],dD[hv])[0,1]
            mk=[gi[g] for g in fmk.get(s,[]) if g in gi]; mk=[i for i in mk if dA[i]>0.1]
            ratio=np.median(dD[mk]/dA[mk]) if mk else np.nan
            r_mk=np.corrcoef(dA[mk],dD[mk])[0,1] if len(mk)>2 else np.nan
            E.append(dict(cohort=coh,combo=lab,dataset=ds,label=col,subset=s,six=f2s.get(s),n=int((L==s).sum()),pct_realA=100*pA[s],
                          r_identity_nonmarker=r_non,marker_contrast_ratio=ratio,r_marker=r_mk))
        if ds=="realB_thin" and col!="fine":
            au=obs.fine.astype(str).values
            V.append(dict(cohort=coh,label=col,ARI=adjusted_rand_score(au,L),accuracy=(au==L).mean(),
                          macroF1=f1_score(au,L,average="macro",labels=list(mA.keys()),zero_division=0)))
    print(coh,"done",flush=True)
P=pd.DataFrame(P);E=pd.DataFrame(E);V=pd.DataFrame(V)
P.to_csv("markerann/proportions_long.csv",index=False);E.to_csv("markerann/expression_similarity_long.csv",index=False);V.to_csv("markerann/realB_marker_accuracy.csv",index=False)
print(V.round(3).to_string())
