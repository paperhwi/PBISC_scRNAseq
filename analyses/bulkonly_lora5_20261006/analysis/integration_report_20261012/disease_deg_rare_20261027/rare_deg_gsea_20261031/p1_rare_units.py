"""Per-patient rare-population shares: real (author labels, real half A + half B) vs PBISC floor (kNN-transferred labels, 60k subset).
Unit = sample for sepsis (sample = patient); patient x group for COVID (visits averaged by pooling cells).
Share = % of all cells of the unit (primary) and % of parent (secondary, parent >= 20 cells).
Stats per subset and group: mean, SD, binomial-corrected SD, Welch t-test real vs floor, Levene test for SD, BH FDR."""
import numpy as np, pandas as pd, anndata as ad, os
from scipy import stats
from rare_lib import SEPSIS, COVID, IN, TR
from deg_lib import bh, COVID_EXCL
os.makedirs("perpatient/out",exist_ok=True)
GRP={"sepsis":{"adult_sepsis":"case","adult_HC":"control","adult_ICU_ctrl":"ICU control"},"covid":{"covid_positive":"case","covid_negative":"control"}}
U=[]
for coh in ["sepsis","covid"]:
    DEF=SEPSIS if coh=="sepsis" else COVID
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad",backed="r").obs; Ub=ad.read_h5ad(f"{IN}{coh}_realU.h5ad",backed="r").obs
    F=pd.read_csv(f"{TR}{coh}_floor_transfer.csv.gz",index_col=0,keep_default_na=False)
    real=pd.concat([R,Ub]); real=pd.DataFrame({"sample":real["sample"].astype(str).values,"patient":real.patient.astype(str).values,"group":real.group.astype(str).values,
                                                "major":real.major.astype(str).values,"fine":real.fine.astype(str).values,"source":"real"})
    flo=pd.DataFrame({"sample":F["sample"].astype(str).values,"group":F.group.astype(str).values,"major":F.major_tr.values,"fine":F.fine_tr.values,"source":"floor"})
    flo=flo[~flo["sample"].isin(COVID_EXCL)]
    flo["patient"]=flo["sample"].str.split("_").str[0] if coh=="covid" else flo["sample"]
    D=pd.concat([real,flo],ignore_index=True); D["grp"]=D.group.map(GRP[coh]); D=D.dropna(subset=["grp"])
    D["unit"]=D.source+"|"+D.patient+"|"+D.grp
    for (unit,src,pat,grp),d in D.groupby(["unit","source","patient","grp"]):
        m=d.major.values; f=d.fine.values; n=len(d)
        for sub,fn in DEF.items():
            par=fn(m,f); k=int((par&(f==sub)).sum()); npar=int(par.sum())
            U.append(dict(cohort=coh,source=src,patient=pat,group=grp,unit=unit,subset=sub,n_cells=n,n_subset=k,n_parent=npar,
                          pct_all=100*k/n,pct_parent=100*k/npar if npar>=20 else np.nan))
U=pd.DataFrame(U); U.to_csv("perpatient/out/rare_unit_shares.csv",index=False)
rows=[]
for (coh,sub,grp),d in U.groupby(["cohort","subset","group"]):
    r=d[d.source=="real"]; f=d[d.source=="floor"]
    if len(r)<2 or len(f)<2: continue
    row=dict(cohort=coh,subset=sub,group=grp,n_real=len(r),n_floor=len(f))
    for nm,x in [("real",r),("floor",f)]:
        p=x.pct_all.values/100; row[f"mean_{nm}"]=100*p.mean(); row[f"sd_{nm}"]=100*p.std(ddof=1)
        binvar=np.mean(p.mean()*(1-p.mean())/x.n_cells.values)                     # expected sampling variance
        row[f"sd_bio_{nm}"]=100*np.sqrt(max(p.var(ddof=1)-binvar,0)); row[f"cells_per_unit_{nm}"]=x.n_cells.median()
    row["fold_floor_vs_real"]=row["mean_floor"]/row["mean_real"] if row["mean_real"]>0 else np.nan
    row["sd_ratio"]=row["sd_floor"]/row["sd_real"] if row["sd_real"]>0 else np.nan
    row["sd_bio_ratio"]=row["sd_bio_floor"]/row["sd_bio_real"] if row["sd_bio_real"]>0 else np.nan
    t=stats.ttest_ind(r.pct_all,f.pct_all,equal_var=False); row["welch_t"]=t.statistic; row["welch_p"]=t.pvalue
    row["mwu_p"]=stats.mannwhitneyu(r.pct_all,f.pct_all).pvalue
    row["levene_p"]=stats.levene(r.pct_all,f.pct_all,center="median").pvalue
    rp=r.pct_parent.dropna(); fp=f.pct_parent.dropna()
    row["n_real_parent"]=len(rp); row["n_floor_parent"]=len(fp)
    if len(rp)>=3 and len(fp)>=3:
        row["mean_parent_real"]=rp.mean(); row["mean_parent_floor"]=fp.mean(); row["welch_p_parent"]=stats.ttest_ind(rp,fp,equal_var=False).pvalue
    rows.append(row)
S=pd.DataFrame(rows)
for c in ["welch_p","levene_p","welch_p_parent","mwu_p"]:
    S[c.replace("_p","_q")]=np.nan
    for g,i in S.groupby(["cohort","group"]).groups.items():
        v=S.loc[i,c]; ok=v.notna(); S.loc[v[ok].index,c.replace("_p","_q")]=bh(v[ok].values)
S.to_csv("perpatient/out/rare_stats.csv",index=False)
pd.set_option("display.width",260)
print(S[["cohort","subset","group","n_real","n_floor","mean_real","mean_floor","fold_floor_vs_real","sd_real","sd_floor","sd_ratio","sd_bio_ratio","welch_p","welch_q","levene_q"]].round(3).to_string(index=False))
