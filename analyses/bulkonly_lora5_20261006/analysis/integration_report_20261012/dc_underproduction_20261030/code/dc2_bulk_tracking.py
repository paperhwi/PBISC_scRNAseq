"""DC-2: does the generated pseudobulk follow its own bulk input for DC markers?
Same samples: bulk input (x) vs pseudobulk of all generated cells (y). Reference: real sc pseudobulk vs bulk on matched sepsis samples."""
import numpy as np, pandas as pd, pickle
from scipy import stats
PB="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/pb/"
m=pickle.load(open("out2/sym2ens.pkl","rb"))
SETS={"pDC (LILRA4, CLEC4C, PTCRA)":["LILRA4","CLEC4C","PTCRA"],
      "cDC (CD1C, CLEC10A, CLEC9A)":["CD1C","CLEC10A","CLEC9A"],
      "monocyte (CD14, LYZ, VCAN)":["CD14","LYZ","VCAN"],
      "B cell (MS4A1, CD79A, CD19)":["MS4A1","CD79A","CD19"],
      "T cell (CD3E, CD3D, CD2)":["CD3E","CD3D","CD2"],
      "NK (GNLY, KLRF1, NCAM1)":["GNLY","KLRF1","NCAM1"]}
def bulk_lcpm(B):
    L=np.log2(B/B.sum(0)*1e6+1); return L
def bulk_score(L,genes):
    ids=[];
    for g in genes:
        cand=[i for i in m.get(g,[]) if i in L.index]
        if cand: ids.append(max(cand,key=lambda i:L.loc[i].mean()))
    return L.loc[ids].mean(0)
def pb_all(name):
    z=np.load(f"{PB}{name}.npz",allow_pickle=True); S=z["S"].astype(float); gr=pd.Series(z["groups"].astype(str))
    samp=gr.str.split("\\|\\|").str[0].values; genes=z["genes"].astype(str)
    df=pd.DataFrame(S,columns=genes); df["s"]=samp; P=df.groupby("s").sum()
    L=np.log2(P.div(P.sum(1),axis=0)*1e6+1); return L
def sc_score(L,genes): return L[[g for g in genes if g in L.columns]].mean(1)
rows=[]; pts=[]
for coh in ["sepsis","covid"]:
    B=pd.read_pickle(f"out2/{coh}_bulk.pkl"); Lb=bulk_lcpm(B)
    for src in ["real","floor","prev"]:
        Lg=pb_all(f"{src}_{coh}")
        common=[s for s in Lg.index if s in Lb.columns]
        if len(common)<10: rows.append(dict(cohort=coh,source=src,n=len(common),note="no sample-level match")); continue
        for k,gs in SETS.items():
            x=bulk_score(Lb[common],gs); y=sc_score(Lg.loc[common],gs)
            r=stats.pearsonr(x,y); sl=np.polyfit(x,y,1)[0]
            rows.append(dict(cohort=coh,source=src,markers=k,n=len(common),r=r.statistic,p=r.pvalue,slope=sl,
                             bulk_mean=x.mean(),sc_mean=y.mean(),bulk_sd=x.std(),sc_sd=y.std()))
            pts.append(pd.DataFrame(dict(cohort=coh,source=src,markers=k,sample=common,bulk=x.values,sc=y.values)))
R=pd.DataFrame(rows); R.to_csv("dc/dc2_bulk_tracking.csv",index=False)
pd.concat(pts).to_csv("dc/dc2_bulk_tracking_points.csv.gz",index=False,compression="gzip")
pd.set_option("display.width",230)
print("=== per-sample: bulk marker score (x) vs pseudobulk of all cells (y), same samples ===")
print(R.dropna(subset=["r"]).pivot_table(index=["cohort","markers"],columns="source",values=["r","slope"]).round(2).to_string())
print("\n(note) COVID real sc samples are date-coded and cannot be matched to bulk samples one-to-one")
# COVID real: patient-level reference
B=pd.read_pickle("out2/covid_bulk.pkl"); cm=pd.read_csv("out2/covid_bulk_meta.csv",index_col=0); Lb=bulk_lcpm(B)
Lg=pb_all("real_covid"); ro=pd.read_csv(PB+"real_covid_obs.csv").set_index("sample")
rows=[]
for k,gs in SETS.items():
    xb=bulk_score(Lb,gs).groupby(cm.loc[Lb.columns,"patient"].values).mean()
    ys=sc_score(Lg,gs); ys=ys.groupby(ro.loc[ys.index,"patient"].astype(str).values).mean()
    c=[p for p in xb.index if p in ys.index]
    r=stats.pearsonr(xb[c],ys[c]); rows.append(dict(markers=k,n_patients=len(c),r=r.statistic,slope=np.polyfit(xb[c],ys[c],1)[0]))
print("\n=== COVID real sc vs bulk, patient-level averages (reference) ===")
print(pd.DataFrame(rows).round(3).to_string(index=False))
