"""scSemiProfiler: disease DEG metrics with the same pseudobulk pipeline as floor (deg_lib).
Designs: (full) generator vs all real samples; (paired) generator vs real restricted to the SAME samples (scSemi regenerates real samples).
Outputs scsemi/out/deg_compare.csv, module_ratios.csv, sd_ratio.csv, allde_scsemi.pkl"""
import numpy as np, pandas as pd, pickle, os
from scipy import stats
import deg_lib as L
from deg_lib import SIX, cpm_matrix, deg, design
os.makedirs("scsemi/out",exist_ok=True)
P3="/mnt/user-data/uploads/integration_report_20261012/three_model_20261022/data/pb/"
def load_any(path,name):
    z=np.load(path+name+".npz",allow_pickle=True)
    S=z["S"].astype(np.float64); groups=z["groups"].astype(str)
    return dict(S=S,n=z["n"],samp=np.array([g.split("||")[0] for g in groups]),six=np.array([g.split("||")[1] for g in groups]),
                genes=z["genes"].astype(str),obs=pd.read_csv(path+name+"_obs.csv"))
def des_scsemi(d,coh,nonA=False):
    o=d["obs"].copy()
    if nonA: o=o[~o.in_A]
    if coh=="sepsis": o["cond"]=o.group.map({"adult_sepsis":"case","adult_HC":"ctrl"}); o["unit"]=o["sample"]
    else:
        o=o[~o["sample"].isin(L.COVID_EXCL)]; o["cond"]=o.group.map({"covid_positive":"case","covid_negative":"ctrl"}); o["unit"]=o["patient"].astype(str)
    return o.dropna(subset=["cond"])[["sample","unit","cond"]]
def compare(R,F):
    F=F.set_index("gene").loc[R.gene]; rd=R.deg.values; fd=F.deg.values
    up_r=(rd&(R.logFC.values>0)).sum(); dn_r=(rd&(R.logFC.values<0)).sum(); up_f=(fd&(F.logFC.values>0)).sum(); dn_f=(fd&(F.logFC.values<0)).sum()
    return dict(logFC_r=float(np.corrcoef(R.logFC,F.logFC)[0,1]),slope=float(np.polyfit(R.logFC.values,F.logFC.values,1)[0]),
                n_DEG_real=int(rd.sum()),n_DEG_gen=int(fd.sum()),up_real=int(up_r),down_real=int(dn_r),up_gen=int(up_f),down_gen=int(dn_f),
                recall=float((rd&fd).sum()/rd.sum()) if rd.sum() else np.nan,precision=float((rd&fd).sum()/fd.sum()) if fd.sum() else np.nan,
                n_case_real=int(R.n_case.iloc[0]),n_ctrl_real=int(R.n_ctrl.iloc[0]),n_case_gen=int(F.n_case.iloc[0]),n_ctrl_gen=int(F.n_ctrl.iloc[0]))
MOD=[("C1Q / complement","sepsis","Monocyte",["C1QA","C1QB","C1QC","RNASE1","CLU"]),("MHC-II","sepsis","Monocyte",["HLA-DRA","HLA-DRB1","HLA-DPA1","CD74"]),
     ("pDC/cDC identity","sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C"]),("sepsis myeloid activation","sepsis","Monocyte",["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1"]),
     ("ISG","covid","Monocyte",["ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1"]),("IFI27 + SIGLEC1","covid","Monocyte",["IFI27","SIGLEC1"])]
rows=[]; allde={}; mods=[]; sds=[]
for coh in ["sepsis","covid"]:
    real=L.load(f"real_{coh}"); floor=L.load(f"floor_{coh}"); sem=load_any(P3,f"scsemi_{coh}")
    genes=list(real["genes"]); assert genes==list(sem["genes"])
    dR=design(f"real_{coh}",real); dF=design(f"floor_{coh}",floor); dS=des_scsemi(sem,coh); dSn=des_scsemi(sem,coh,nonA=True)
    dRp=dR[dR["sample"].isin(dS["sample"])]                       # real, same samples as scSemi
    for ct in SIX:
        T={}
        for nm,(d,ds) in {"real":(real,dR),"real_paired":(real,dRp),"floor":(floor,dF),"scsemi":(sem,dS),"scsemi_nonA":(sem,dSn)}.items():
            g=cpm_matrix(d,ds,ct)
            if len(g)==0: continue
            t=deg(g,genes);
            if t is not None: T[nm]=(t,g)
        if "real" not in T: continue
        meas=(real["S"].sum(0)>0)&(floor["S"].sum(0)>0)&(sem["S"].sum(0)>0); gm=set(np.array(genes)[meas])
        tabs={k:v[0][v[0].gene.isin(gm)].reset_index(drop=True) for k,v in T.items()}; allde[(coh,ct)]=tabs
        for gen,ref in [("floor","real"),("scsemi","real"),("scsemi_nonA","real"),("scsemi","real_paired")]:
            if gen in tabs and ref in tabs:
                rows.append(dict(cohort=coh,celltype=ct,generator=gen,reference=ref,**compare(tabs[ref],tabs[gen])))
        # between-sample SD (expressed genes), generator vs its reference
        for gen,ref in [("floor","real"),("scsemi","real_paired"),("scsemi","real")]:
            if gen not in T or ref not in T or ct=="DC": continue
            for grp in ["case","ctrl"]:
                a=T[ref][1]; b=T[gen][1]; ga=[x for x in genes if x in gm]
                xa=a[a._c==grp][ga]; xb=b[b._c==grp][ga]
                ok=(xa.mean(0)>1)&(xb.mean(0)>1); rr=(xb.std(0)[ok]/xa.std(0)[ok]).replace([np.inf,-np.inf],np.nan).dropna()
                sds.append(dict(cohort=coh,celltype=ct,group=grp,generator=gen,reference=ref,n_units_ref=len(xa),n_units_gen=len(xb),median_SD_ratio=float(rr.median())))
for nm,coh,ct,gs in MOD:
    tabs=allde[(coh,ct)]
    for gen,ref in [("floor","real"),("scsemi","real"),("scsemi","real_paired")]:
        if gen not in tabs or ref not in tabs: continue
        R=tabs[ref].set_index("gene"); F=tabs[gen].set_index("gene"); gg=[g for g in gs if g in R.index and abs(R.logFC[g])>0.2]
        rat=[F.logFC[g]/R.logFC[g] for g in gg]
        mods.append(dict(module=nm,cohort=coh,celltype=ct,generator=gen,reference=ref,n_genes=len(gg),median_ratio=float(np.median(rat)) if rat else np.nan,
                         genes_DEG_real=int(R.loc[gg].deg.sum()),genes_DEG_gen=int(F.loc[gg].deg.sum()),detail=";".join(f"{g}:{R.logFC[g]:+.2f}->{F.logFC[g]:+.2f}" for g in gg)))
D=pd.DataFrame(rows); M=pd.DataFrame(mods); SD=pd.DataFrame(sds)
D.to_csv("scsemi/out/deg_compare.csv",index=False); M.to_csv("scsemi/out/module_ratios.csv",index=False); SD.to_csv("scsemi/out/sd_ratio.csv",index=False)
pickle.dump(allde,open("scsemi/out/allde_scsemi.pkl","wb"))
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",90)
print(D.round(3).to_string(index=False)); print(M.drop(columns="detail").round(3).to_string(index=False))
print(SD.groupby(["cohort","generator","reference"]).median_SD_ratio.agg(["min","median","max"]).round(2))
