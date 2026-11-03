"""DC-1: how big is the DC deficit, in which group, which DC subtype, under which label route."""
import numpy as np, pandas as pd
from scipy import stats
PB="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/pb/"
TR="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/results/big/"
def grp(c,s,o):
    if c=="sepsis":
        if "age_group" in o:   # real
            m={"sepsis":"case","healthy control":"HC","ICU control":"ICU"}; o=o[o.age_group=="adult"].copy(); o["g"]=o.group.map(m)
        else:
            m={"adult_sepsis":"case","adult_HC":"HC","adult_ICU_ctrl":"ICU"}; o=o.copy(); o["g"]=o.group.map(m)
    else:
        if "case_control" in o: o=o.copy(); o["g"]=o.case_control.map({"POSITIVE":"case","NEGATIVE":"HC"})
        else: o=o.copy(); o["g"]=o.group.map({"covid_positive":"case","covid_negative":"HC"})
    return o.dropna(subset=["g"])[["sample","g"]]
rows=[]
for c in ["sepsis","covid"]:
    for s in ["real","floor","prev"]:
        z=np.load(f"{PB}{s}_{c}.npz",allow_pickle=True); o=pd.read_csv(f"{PB}{s}_{c}_obs.csv")
        gr=pd.Series(z["groups"].astype(str)); n=z["n"]
        df=pd.DataFrame({"sample":gr.str.split("\\|\\|").str[0],"six":gr.str.split("\\|\\|").str[1],"n":n})
        w=df.pivot_table(index="sample",columns="six",values="n",aggfunc="sum",fill_value=0)
        w["tot"]=w.sum(1); w["dc_pct"]=100*w.get("DC",0)/w.tot
        w=w.join(grp(c,s,o).set_index("sample"),how="inner")
        for g,d in w.groupby("g"):
            rows.append(dict(cohort=c,source=s,group=g,n_samples=len(d),
                             DC_pct_median=d.dc_pct.median(),DC_pct_mean=d.dc_pct.mean(),
                             pct_samples_DC0=100*(d.get("DC",0)==0).mean(),
                             DC_cells_per_sample_median=float(np.median(d.get("DC",0)))))
A=pd.DataFrame(rows); A.to_csv("dc/dc1_fraction_by_group.csv",index=False)
pd.set_option("display.width",220)
print("=== DC % of all cells per sample (full files; real = author labels, generated = frozen classifier) ===")
print(A.round(3).to_string(index=False))
# ---- subtype via kNN-transferred fine labels on the 60k sets (same route for real B and floor) ----
rows=[]
for c in ["sepsis","covid"]:
    for s in ["realU","floor"]:
        t=pd.read_csv(f"{TR}{c}_{s}_transfer.csv.gz",index_col=0,keep_default_na=False)
        g=t.group.astype(str)
        if c=="sepsis": t["g"]=g.map({"adult_sepsis":"case","adult_HC":"HC","adult_ICU_ctrl":"ICU"})
        else: t["g"]=g.map({"covid_positive":"case","covid_negative":"HC"})
        t=t.dropna(subset=["g"])
        for gg,d in t.groupby("g"):
            per=d.groupby("sample").apply(lambda x: pd.Series({
                "pDC":100*(x.fine_tr=="pDC").mean() if c=="sepsis" else 100*(x.fine_tr=="pDC").mean(),
                "cDC":100*x.fine_tr.isin(["cDC","DC1","DC2","DC3","DC_prolif"]).mean() if c=="sepsis" else 100*x.fine_tr.str.startswith("DC").mean(),
                "ASDC":100*(x.fine_tr=="ASDC").mean(),
                "DC_six_tr":100*(x.six_tr=="DC").mean(),
                "DC_frozen":100*(x.six_type=="DC").mean()}),include_groups=False)
            rows.append(dict(cohort=c,source=s,group=gg,n_samples=len(per),**{k:per[k].median() for k in per.columns},
                             **{k+"_mean":per[k].mean() for k in ["pDC","cDC"]}))
B=pd.DataFrame(rows); B.to_csv("dc/dc1_subtype_transfer.csv",index=False)
print("\n=== DC subtypes, % of all cells per sample (kNN-transferred labels, 60k balanced sets; median) ===")
print(B.round(3).to_string(index=False))
