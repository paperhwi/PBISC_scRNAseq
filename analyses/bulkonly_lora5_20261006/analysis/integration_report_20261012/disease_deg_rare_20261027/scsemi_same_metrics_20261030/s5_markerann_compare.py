"""Share error and non-marker identity r vs real A: scSemi (kNN / marker) next to floor and real B, on marker-resolvable subsets (real B F1>=0.5)."""
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, pickle, warnings; warnings.filterwarnings("ignore")
B="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
M=pickle.load(open("markerann/markers.pkl","rb")); F1=pd.read_csv("markerann/realB_per_subset_F1.csv"); F1=F1[F1.label=="fine_mk_cl"].set_index(["cohort","subset"]).F1
P0=pd.read_csv("markerann/proportions_long.csv"); E0=pd.read_csv("markerann/expression_similarity_long.csv")
def means(X,lab):
    lab=pd.Series(lab); return {g:np.asarray(X[np.asarray(i)].mean(0)).ravel() for g,i in lab.groupby(lab).groups.items() if len(i)>=20 and g!=""}
P=[];E=[]
for coh in ["sepsis","covid"]:
    A=ad.read_h5ad(B+f"{coh}_ref.h5ad"); sc.pp.normalize_total(A,1e4); sc.pp.log1p(A)
    fineA=A.obs.fine.astype(str).values; allmk=set(g for d in M[coh]["L2"].values() for v in d.values() for g in v)|set(g for v in M[coh]["L1"].values() for g in v)
    H=A.copy(); sc.pp.highly_variable_genes(H,n_top_genes=2000,flavor="seurat"); hvs=set(H.var_names[H.var.highly_variable]); del H
    hv=np.array([g in hvs and g not in allmk for g in A.var_names])
    mA=means(A.X,fineA); muA=np.asarray(A.X.mean(0)).ravel(); pA=pd.Series(fineA).value_counts(normalize=True)
    D=ad.read_h5ad(f"scsemi/out/{coh}_scsemi_annot.h5ad"); muD=np.asarray(D.X.mean(0)).ravel()
    for col,lab in [("fine_knn","scSemi · kNN transfer"),("fine_mk_cl","scSemi · marker (cluster)")]:
        L=D.obs[col].fillna("").astype(str).values; pD=pd.Series(L).value_counts(normalize=True); mD=means(D.X,L)
        for s in pA.index:
            P.append(dict(cohort=coh,combo=lab,subset=s,pct_realA=100*pA[s],pct=100*pD.get(s,0.0)))
            if s in mA and s in mD:
                E.append(dict(cohort=coh,combo=lab,subset=s,r_identity_nonmarker=np.corrcoef((mA[s]-muA)[hv],(mD[s]-muD)[hv])[0,1]))
P=pd.concat([P0[["cohort","combo","subset","pct_realA","pct"]],pd.DataFrame(P)]); E=pd.concat([E0[["cohort","combo","subset","r_identity_nonmarker"]],pd.DataFrame(E)])
for D_ in (P,E): D_["res"]=[F1.get((c,s),0)>=0.5 for c,s in zip(D_.cohort,D_.subset)]
P["l2r"]=np.log2((P.pct+0.01)/(P.pct_realA+0.01))
S=P[P.res].groupby(["cohort","combo"]).l2r.agg(med_abs=lambda x: np.median(np.abs(x)),n_lost=lambda x:(x<-2).sum(),n="size").join(
  E[E.res].groupby(["cohort","combo"]).r_identity_nonmarker.agg(["median","count"]))
S.to_csv("scsemi/out/markerann_compare.csv"); P.to_csv("scsemi/out/markerann_props_long.csv",index=False); E.to_csv("scsemi/out/markerann_identity_long.csv",index=False)
pd.set_option("display.width",200); print(S.round(2).to_string())
