import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"font.size":9,"axes.titlesize":10,"legend.fontsize":8})
BLUE,RED,GREY,INK,GREEN="#2B83BA","#D7191C","#7F7F7F","#1F2933","#1A9641"; SEMI="#E08214"
D=pd.read_csv("scsemi/out/deg_compare.csv"); M=pd.read_csv("scsemi/out/module_ratios.csv"); SD=pd.read_csv("scsemi/out/sd_ratio.csv")
A=pd.read_csv("out2/a3_power_matched_ceiling.csv").merge(pd.read_csv("out/real_real_ceiling.csv")[["cohort","celltype","logFC_r"]].rename(columns={"logFC_r":"ceil_both"}),on=["cohort","celltype"]).set_index(["cohort","celltype"])
fig,ax=plt.subplots(2,2,figsize=(12.4,7.6))
# a logFC r
d=D[D.reference=="real"].pivot_table(index=["cohort","celltype"],columns="generator",values="logFC_r")
order=[(c,t) for c in ["sepsis","covid"] for t in ["T cell","B cell","NK cell","Monocyte","DC"]]; d=d.reindex(order); x=np.arange(len(d)); w=.38
a=ax[0,0]; a.bar(x-w/2,d.floor,w,color=BLUE,label="PBISC floor"); a.bar(x+w/2,d.scsemi,w,color=SEMI,label="scSemiProfiler")
for i,k in enumerate(order):
    lo,hi=sorted([A.loc[k,"ceil_case_split_r"],A.loc[k,"ceil_both"]]); a.plot([i-.45,i+.45],[hi,hi],color=INK,lw=1); a.plot([i-.45,i+.45],[lo,lo],color=GREY,lw=1)
a.set_xticks(x); a.set_xticklabels([f"{c[:3]} {t.replace(' cell','')}" for c,t in order],rotation=35,ha="right"); a.set_ylim(0,1.05)
a.set_ylabel("disease logFC r vs real"); a.set_title("a  disease DEG direction (lines = real-real ceilings)",loc="left"); a.legend(frameon=False,loc="upper right",ncol=2)
# b modules
m=M[M.reference=="real"].pivot_table(index="module",columns="generator",values="median_ratio")
mo=["C1Q / complement","pDC/cDC identity","MHC-II","sepsis myeloid activation","ISG","IFI27 + SIGLEC1"]; m=m.reindex(mo); x=np.arange(len(m))
a=ax[0,1]; a.bar(x-w/2,m.floor,w,color=BLUE); a.bar(x+w/2,m.scsemi,w,color=SEMI)
for i,(f,s) in enumerate(zip(m.floor,m.scsemi)):
    a.text(i-w/2,f*1.08,f"{f:.2f}",ha="center",fontsize=7); a.text(i+w/2,s*1.08,f"{s:.2f}",ha="center",fontsize=7)
a.set_yscale("log"); a.set_ylim(0.05,6); a.axhline(1,ls="--",color=INK,lw=.8); a.axhline(.5,ls=":",color=RED,lw=.8)
a.set_xticks(x); a.set_xticklabels(["C1Q\nsep Mono","DC identity\nsep DC","MHC-II\nsep Mono","myeloid act.\nsep Mono","ISG\nCOV Mono","IFI27+SIGLEC1\nCOV Mono"],fontsize=7,rotation=30,ha="right")
a.set_ylabel("generator / real disease log2FC (median gene)"); a.set_title("b  module recovery (1 = perfect)",loc="left")
# c SD ratio & down DEG
a=ax[1,0]
s=SD[((SD.generator=="floor")&(SD.reference=="real"))|((SD.generator=="scsemi")&(SD.reference=="real_paired"))]
for j,(g,col) in enumerate([("floor",BLUE),("scsemi",SEMI)]):
    v=s[s.generator==g]
    for k,coh in enumerate(["sepsis","covid"]):
        vv=v[v.cohort==coh].median_SD_ratio.values; xx=np.full(len(vv),k*3+j)+np.random.default_rng(0).uniform(-.18,.18,len(vv))
        a.scatter(xx,vv,s=22,color=col,edgecolor="white",lw=.4,zorder=3); a.plot([k*3+j-.3,k*3+j+.3],[np.median(vv)]*2,color=INK,lw=1.5)
a.axhline(1,ls="--",color=INK,lw=.8); a.set_yscale("log"); a.set_ylim(0.4,16)
a.set_xticks([0,1,3,4]); a.set_xticklabels(["sepsis\nfloor","sepsis\nscSemi","COVID\nfloor","COVID\nscSemi"])
a.set_ylabel("between-sample SD, generator / real\n(median gene; dots = cell type × group)"); a.set_title("c  sample-to-sample variation",loc="left")
dd=D[D.reference.isin(["real"])&D.generator.isin(["floor","scsemi"])].groupby(["cohort","generator"])[["up_real","down_real","up_gen","down_gen"]].sum()
txt="\n".join(f"{c} {g}: up {int(r.up_real)}→{int(r.up_gen)}, down {int(r.down_real)}→{int(r.down_gen)}" for (c,g),r in dd.iterrows())
a.text(0.02,0.98,"DEG counts, 5 types summed (real → generator)\n"+txt,transform=a.transAxes,va="top",fontsize=7,color=INK)
# d rare transfer-free
a=ax[1,1]; rows=[]
for coh in ["sepsis","covid"]:
    t=pd.read_csv(f"/mnt/user-data/uploads/integration_report_20261012/three_model_20261022/results/{coh}_rare_transfer_free.csv")
    t["ratio"]=t.gen_estimate_pct/t.realA_cv_estimate_pct; t["ratioB"]=t.realB_estimate_pct/t.realA_cv_estimate_pct
    p=t.pivot_table(index="subset",columns="model_depth",values="ratio"); p["realB"]=t.groupby("subset").ratioB.first(); p["cohort"]=coh; rows.append(p.reset_index())
P=pd.concat(rows,ignore_index=True)
for j,(c,col,lab) in enumerate([("realB",GREY,"real B"),("floor",BLUE,"floor"),("scsemi",SEMI,"scSemi")]):
    v=np.clip(P[c].values,0.005,30); xx=j+np.random.default_rng(1).uniform(-.2,.2,len(v))
    a.scatter(xx,v,s=20,color=col,edgecolor="white",lw=.4,zorder=3); a.plot([j-.3,j+.3],[np.median(v)]*2,color=INK,lw=1.5)
    a.text(j,22,f"<0.5×: {(P[c]<0.5).sum()}/{len(v)}",ha="center",fontsize=8)
a.set_yscale("log"); a.set_ylim(0.004,40); a.axhspan(0.5,2,color="#EAF3FA",zorder=0); a.set_xticks([0,1,2]); a.set_xticklabels(["real B","PBISC floor","scSemiProfiler"])
a.set_ylabel("rare subset share / real A (transfer-free, log)"); a.set_title("d  rare populations (27 subsets, both cohorts; band = 0.5–2×)",loc="left")
for k in ["Macrophage","Mono_c06-CD14-MKI67","CD14mono_dex","pDC"]:
    for _,r in P[P.subset==k].iterrows():
        a.annotate(f"{k[:12]} ({r.cohort[:3]})",(2,max(r.scsemi,0.005)),xytext=(2.32,max(r.scsemi,0.005)),fontsize=6.5,color=SEMI,va="center",arrowprops=dict(arrowstyle="-",color=SEMI,lw=.5))
fig.tight_layout(); fig.savefig("scsemi/out/fig_scsemi_vs_floor.png",dpi=220,bbox_inches="tight",facecolor="white"); print("ok")
