import numpy as np, pandas as pd, pickle
from scipy import stats
from deg_lib import *
rng=np.random.default_rng(20261027)
CT=["T cell","B cell","NK cell","Monocyte","DC"]
SRC={}
for coh in ["sepsis","covid"]:
    for s in ["real","floor"]:
        d=load(f"{s}_{coh}"); SRC[(coh,s)]=(d,design(f"{s}_{coh}",d))
genes={coh:load(f"real_{coh}")["genes"] for coh in ["sepsis","covid"]}
meas={}
for coh in ["sepsis","covid"]:
    r=load(f"real_{coh}"); f=load(f"floor_{coh}")
    meas[coh]=set(r["genes"][(r["S"].sum(0)>0)&(f["S"].sum(0)>0)])

def mat(coh,src,ct):
    d,des=SRC[(coh,src)]; return cpm_matrix(d,des,ct)
def tab(g,coh,units_case=None,units_ctrl=None):
    if units_case is not None:
        g=g[((g._c=="case")&g._u.isin(units_case))|((g._c=="ctrl")&g._u.isin(units_ctrl))]
    t=deg(g,list(genes[coh]))
    return None if t is None else t[t.gene.isin(meas[coh])].reset_index(drop=True)
def cmp(A,B):
    B=B.set_index("gene").loc[A.gene]
    r=float(np.corrcoef(A.logFC,B.logFC)[0,1])
    rec=float((A.deg.values&B.deg.values).sum()/max(A.deg.sum(),1)) if A.deg.sum() else np.nan
    pre=float((A.deg.values&B.deg.values).sum()/max(B.deg.sum(),1)) if B.deg.sum() else np.nan
    return r,rec,pre

# ================= A3: power-matched ceiling =================
rows=[]
for coh in ["sepsis","covid"]:
    for ct in CT:
        gR=mat(coh,"real",ct); gF=mat(coh,"floor",ct)
        if len(gR)==0 or len(gF)==0: continue
        cR=gR[gR._c=="case"]._u.unique(); oR=gR[gR._c=="ctrl"]._u.unique()
        full=tab(gR,coh); fl=tab(gF,coh)
        if full is None or fl is None: continue
        r_full,rec_full,pre_full=cmp(full,fl)
        acc=[]
        for _ in range(30):
            ca=rng.permutation(cR); h=len(ca)//2
            # V1: split case arm only, controls shared
            A=tab(gR,coh,ca[:h],oR); B=tab(gR,coh,ca[h:],oR)
            if A is None or B is None: continue
            v1=cmp(A,B)
            # V2: floor subsampled to the same unit counts as the halves
            cF=gF[gF._c=="case"]._u.unique(); oF=gF[gF._c=="ctrl"]._u.unique()
            fs=tab(gF,coh,rng.choice(cF,min(h,len(cF)),replace=False),rng.choice(oF,min(len(oR),len(oF)),replace=False))
            v2=cmp(A,fs) if fs is not None else (np.nan,)*3
            acc.append(v1+v2)
        if not acc: continue
        a=np.array(acc,float)
        rows.append(dict(cohort=coh,celltype=ct,n_case=len(cR),n_ctrl=len(oR),
            floor_r_full=r_full,floor_recall_full=rec_full,floor_prec_full=pre_full,
            ceil_case_split_r=np.nanmean(a[:,0]),ceil_case_split_r_sd=np.nanstd(a[:,0]),
            ceil_recall=np.nanmean(a[:,1]),ceil_prec=np.nanmean(a[:,2]),
            floor_r_matchedN=np.nanmean(a[:,3]),floor_recall_matchedN=np.nanmean(a[:,4]),floor_prec_matchedN=np.nanmean(a[:,5])))
A3=pd.DataFrame(rows)
A3["pct_of_ceiling_r"]=100*A3.floor_r_matchedN/A3.ceil_case_split_r
A3.to_csv("out2/a3_power_matched_ceiling.csv",index=False)
pd.set_option("display.width",260)
print("=== A3: control-arm-shared ceiling, floor at matched n ===")
print(A3[["cohort","celltype","n_case","n_ctrl","ceil_case_split_r","ceil_case_split_r_sd","floor_r_matchedN","pct_of_ceiling_r","ceil_recall","floor_recall_matchedN"]].round(3).to_string(index=False))

# ================= D1: patient bootstrap CI =================
rows=[]
for coh in ["sepsis","covid"]:
    for ct in CT:
        gR=mat(coh,"real",ct); gF=mat(coh,"floor",ct)
        if len(gR)==0 or len(gF)==0: continue
        cR=gR[gR._c=="case"]._u.unique(); oR=gR[gR._c=="ctrl"]._u.unique()
        cF=gF[gF._c=="case"]._u.unique(); oF=gF[gF._c=="ctrl"]._u.unique()
        boot=[]
        for _ in range(200):
            A=tab(gR,coh,rng.choice(cR,len(cR),True),rng.choice(oR,len(oR),True))
            B=tab(gF,coh,rng.choice(cF,len(cF),True),rng.choice(oF,len(oF),True))
            if A is None or B is None: continue
            boot.append(cmp(A,B))
        if not boot: continue
        b=np.array(boot,float)
        q=lambda i:(np.nanpercentile(b[:,i],2.5),np.nanpercentile(b[:,i],97.5))
        full=tab(gR,coh); fl=tab(gF,coh); r,rec,pre=cmp(full,fl)
        rows.append(dict(cohort=coh,celltype=ct,logFC_r=r,r_lo=q(0)[0],r_hi=q(0)[1],
                         recall=rec,rec_lo=q(1)[0],rec_hi=q(1)[1],precision=pre,pre_lo=q(2)[0],pre_hi=q(2)[1]))
D1=pd.DataFrame(rows); D1.to_csv("out2/d1_bootstrap_ci.csv",index=False)
print("\n=== D1: patient bootstrap 95% CI (200 reps) ===")
print(D1.round(3).to_string(index=False))

# ================= D2: negative controls =================
rows=[]
for coh in ["sepsis","covid"]:
    for ct in CT:
        gR=mat(coh,"real",ct); gF=mat(coh,"floor",ct)
        if len(gR)==0 or len(gF)==0: continue
        # permutation null: shuffle condition labels in BOTH, 20 reps
        perm=[]
        for _ in range(20):
            a=gR.copy(); a["_c"]=rng.permutation(a._c.values)
            b=gF.copy(); b["_c"]=rng.permutation(b._c.values)
            A=tab(a,coh); B=tab(b,coh)
            if A is None or B is None: continue
            perm.append(cmp(A,B))
        p=np.array(perm,float) if perm else np.full((1,3),np.nan)
        full=tab(gR,coh); fl=tab(gF,coh); r,rec,pre=cmp(full,fl)
        rows.append(dict(cohort=coh,celltype=ct,observed_r=r,perm_r_mean=np.nanmean(p[:,0]),
                         perm_r_p975=np.nanpercentile(p[:,0],97.5),
                         perm_nDEG_real=np.nanmean([np.nan]),observed_recall=rec,perm_recall=np.nanmean(p[:,1])))
D2=pd.DataFrame(rows).drop(columns=["perm_nDEG_real"])
D2.to_csv("out2/d2_permutation_null.csv",index=False)
print("\n=== D2a: label-permutation null ===")
print(D2.round(3).to_string(index=False))

# D2b: sepsis vs ICU control (known hard negative)
import copy
rows=[]
d=load("real_sepsis"); o=d["obs"]; o=o[o.age_group=="adult"]
dr=o.assign(cond=o.group.map({"sepsis":"case","ICU control":"ctrl"})).dropna(subset=["cond"]).assign(unit=lambda x:x["sample"])[["sample","unit","cond"]]
df=load("floor_sepsis"); of=df["obs"]
dfm=of.assign(cond=of.group.map({"adult_sepsis":"case","adult_ICU_ctrl":"ctrl"})).dropna(subset=["cond"]).assign(unit=lambda x:x["sample"])[["sample","unit","cond"]]
for ct in CT:
    A=deg(cpm_matrix(d,dr,ct),list(genes["sepsis"])); B=deg(cpm_matrix(df,dfm,ct),list(genes["sepsis"]))
    if A is None or B is None: continue
    A=A[A.gene.isin(meas["sepsis"])].reset_index(drop=True); B=B[B.gene.isin(meas["sepsis"])].reset_index(drop=True)
    r,rec,pre=cmp(A,B)
    rows.append(dict(celltype=ct,n_realDEG=int(A.deg.sum()),n_floorDEG=int(B.deg.sum()),logFC_r=r,recall=rec,precision=pre))
D2b=pd.DataFrame(rows); D2b.to_csv("out2/d2b_sepsis_vs_icu.csv",index=False)
print("\n=== D2b: sepsis vs ICU control (hard negative) ===")
print(D2b.round(3).to_string(index=False))
