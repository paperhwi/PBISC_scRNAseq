# Step 2: de novo clustering + Top10-marker annotation (two-level) of floor and real B
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, pickle, sys, warnings; warnings.filterwarnings("ignore")
import scipy.sparse as sp
B="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/"
M=pickle.load(open("markerann/markers.pkl","rb"))
rng=np.random.default_rng(0)
def thin(X,p):
    X=sp.csr_matrix(X).astype(np.int64); X.sort_indices()
    rows=np.repeat(np.arange(X.shape[0]),np.diff(X.indptr))
    X.data=rng.binomial(X.data,p[rows]); X.eliminate_zeros(); return X.astype(np.float32)
def cluster(a,res,npc):
    b=a.copy(); sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor="seurat")
    b=b[:,b.var.highly_variable].copy(); sc.pp.scale(b,max_value=10)
    sc.tl.pca(b,n_comps=min(npc,b.n_obs-1)); sc.pp.neighbors(b,n_neighbors=15)
    sc.tl.leiden(b,resolution=res,key_added="cl",flavor="igraph",n_iterations=2,directed=False)
    return b.obs.cl.astype(str).values
def scores(a,sets):
    S={}
    for k,g in sets.items():
        g=[x for x in g if x in a.var_names]
        sc.tl.score_genes(a,g,score_name="_s",random_state=0,use_raw=False); S[k]=a.obs["_s"].values.copy()
    return pd.DataFrame(S,index=a.obs_names)
def annotate(a,coh):
    L1,L2=M[coh]["L1"],M[coh]["L2"]
    a.obs["cl1"]=cluster(a,1.0,30)
    S1=scores(a,L1)
    cm=S1.groupby(a.obs.cl1.values).mean()
    a.obs["six_mk_cl"]=a.obs.cl1.map(cm.idxmax(1)).values
    a.obs["six_mk_cell"]=S1.idxmax(1).values
    a.obs["fine_mk_cl"]=""; a.obs["fine_mk_cell"]=""; a.obs["cl2"]=""
    for s in L1:
        m=(a.obs.six_mk_cl==s).values
        if m.sum()>=30:
            b=a[m].copy(); b.obs["cl2"]=cluster(b,2.0,20)
            S2=scores(b,L2[s]); cm2=S2.groupby(b.obs.cl2.values).mean()
            a.obs.loc[m,"fine_mk_cl"]=b.obs.cl2.map(cm2.idxmax(1)).values
            a.obs.loc[m,"cl2"]=(s+"_"+b.obs.cl2.astype(str)).values
        mc=(a.obs.six_mk_cell==s).values
        if mc.sum():
            b=a[mc].copy(); S2=scores(b,L2[s]); a.obs.loc[mc,"fine_mk_cell"]=S2.idxmax(1).values
    return a
for coh in ["sepsis","covid"]:
    f=ad.read_h5ad(B+f"big/{coh}_floor.h5ad")
    j=ad.read_h5ad(B+f"big_joint/{coh}_joint_floor.h5ad",backed="r"); jf=j.obs[j.obs.source=="floor"]
    f.obs["six_knn"]=jf.six_tr.astype(str).values; f.obs["fine_knn"]=jf.fine_tr.astype(str).values
    u=ad.read_h5ad(B+f"big/{coh}_realU.h5ad")
    flib=pd.Series(np.asarray(f.X.sum(1)).ravel()).groupby(f.obs.six_knn.values).median()
    ulib=pd.Series(np.asarray(u.X.sum(1)).ravel()).groupby(u.obs.six_type.astype(str).values).median()
    p=(flib/ulib).clip(upper=1).reindex(u.obs.six_type.astype(str).values).fillna(1).values
    print(coh,"thinning p per six:",(flib/ulib).round(2).to_dict())
    ut=u.copy(); ut.X=thin(u.X,p)
    for name,a in [("floor",f),("realB_thin",ut),("realB_full",u)]:
        a.obs_names_make_unique()
        sc.pp.normalize_total(a,1e4); sc.pp.log1p(a)
        a=annotate(a,coh)
        cols=[c for c in ["sample","patient","group","six_type","fine","six_knn","fine_knn","cl1","cl2","six_mk_cl","fine_mk_cl","six_mk_cell","fine_mk_cell"] if c in a.obs]
        a.obs[cols].to_csv(f"markerann/{coh}_{name}_labels.csv.gz")
        if name!="realB_full": a.write_h5ad(f"markerann/{coh}_{name}_annot.h5ad",compression="gzip")
        print(coh,name,"n clusters L1",a.obs.cl1.nunique(),"L2",a.obs.cl2.nunique(),flush=True)
