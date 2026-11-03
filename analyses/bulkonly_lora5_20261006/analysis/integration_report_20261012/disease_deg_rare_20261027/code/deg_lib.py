import numpy as np, pandas as pd
from scipy import stats
PB="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/pb/"
COVID_EXCL={"C126_positive_5","C126_positive_9","C126_positive_11","C141_negative","C58_positive_2"}
SIX=["T cell","B cell","NK cell","Monocyte","DC"]
MINCELL=20

def bh(p):
    p=np.asarray(p,float); ok=~np.isnan(p); q=np.full(p.shape,np.nan)
    pv=p[ok]; n=len(pv); o=np.argsort(pv); r=np.arange(1,n+1)
    qv=np.minimum.accumulate((pv[o]*n/r)[::-1])[::-1]
    out=np.empty(n); out[o]=np.minimum(qv,1); q[ok]=out; return q

def load(name):
    z=np.load(PB+name+".npz",allow_pickle=True)
    S=z["S"].astype(np.float64); n=z["n"]; groups=z["groups"].astype(str); genes=z["genes"].astype(str)
    samp=np.array([g.split("||")[0] for g in groups]); six=np.array([g.split("||")[1] for g in groups])
    obs=pd.read_csv(PB+name+"_obs.csv")
    return dict(S=S,n=n,samp=samp,six=six,genes=genes,obs=obs)

def design(name,d):
    """return sample -> (unit, condition) ; condition in {'case','ctrl'} else None"""
    o=d["obs"].copy()
    if name.endswith("sepsis"):
        if "age_group" in o:  # real
            o=o[o.age_group=="adult"]
            m={"sepsis":"case","healthy control":"ctrl"}
            o["cond"]=o.group.map(m)
        else:
            m={"adult_sepsis":"case","adult_HC":"ctrl"}
            o["cond"]=o.group.map(m)
        o["unit"]=o["sample"]
    else:
        if "case_control" in o:  # real covid
            m={"POSITIVE":"case","NEGATIVE":"ctrl"}
            o["cond"]=o.case_control.map(m); o["unit"]=o["patient"]
        else:
            o=o[~o["sample"].isin(COVID_EXCL)]
            m={"covid_positive":"case","covid_negative":"ctrl"}
            o["cond"]=o.group.map(m)
            o["unit"]=o["sample"].str.split("_").str[0]
    o=o.dropna(subset=["cond"])
    return o[["sample","unit","cond"]]

def cpm_matrix(d,des,ct):
    """samples x genes log2CPM for one six_type, then averaged to unit level"""
    sel=(d["six"]==ct)
    S=d["S"][sel]; n=d["n"][sel]; samp=d["samp"][sel]
    keep=n>=MINCELL
    S,samp=S[keep],samp[keep]
    tot=S.sum(1,keepdims=True); tot[tot==0]=1
    L=np.log2(S/tot*1e6+1)
    df=pd.DataFrame(L,index=samp,columns=d["genes"])
    df=df.loc[df.index.isin(des["sample"])]
    mp=des.set_index("sample")
    df["_u"]=mp.loc[df.index,"unit"].values; df["_c"]=mp.loc[df.index,"cond"].values
    g=df.groupby(["_u","_c"]).mean()
    g=g.reset_index()
    return g

def deg(g,genes,lfc=1.0,q=0.05):
    a=g[g._c=="case"][genes].to_numpy(); b=g[g._c=="ctrl"][genes].to_numpy()
    if len(a)<3 or len(b)<3: return None
    lf=a.mean(0)-b.mean(0)
    t,p=stats.ttest_ind(a,b,equal_var=False)
    qv=bh(p)
    return pd.DataFrame(dict(gene=genes,logFC=lf,p=p,q=qv,
                             mean_case=a.mean(0),mean_ctrl=b.mean(0),
                             n_case=len(a),n_ctrl=len(b),
                             deg=(qv<q)&(np.abs(lf)>=lfc)))
