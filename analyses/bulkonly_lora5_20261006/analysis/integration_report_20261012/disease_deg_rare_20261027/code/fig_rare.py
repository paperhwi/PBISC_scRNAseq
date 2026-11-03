import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,GREEN,ORANGE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#1A9641","#E08214"
F="figs/"; save=lambda fig,n:(fig.savefig(F+n+".png",dpi=220,bbox_inches="tight",facecolor="white"),plt.close(fig),print("saved",n))
R=pd.read_csv("out/rare_expression_ratios.csv"); P=pd.read_csv("out/rare_power_matched.csv")
RA=pd.read_csv("out/rare_expression.csv")

# ---- FIG 8: marker programme retention ----
d=R.dropna(subset=["mkRatio_floor"]).sort_values("mkRatio_floor")
fig,ax=plt.subplots(figsize=(8.0,4.9))
y=np.arange(len(d))
ax.hlines(y,d.mkRatio_floor,d.mkRatio_realB.fillna(1.0),color="#D5DBE1",lw=2.2,zorder=1)
ax.scatter(d.mkRatio_realB,y,s=46,color=GREY,zorder=3)
ax.scatter(d.mkRatio_floor,y,s=46,color=BLUE,zorder=3)
ax.axvline(1,ls="--",lw=.9,c=INK); ax.axvline(0.5,ls=":",lw=.9,c=RED)
ax.set_yticks(y); ax.set_yticklabels([f"{r.subset}  ({r.cohort})" for r in d.itertuples()],fontsize=8)
ax.set_xlabel("subset marker-score / real A  (top-10 markers, log1p CP10k)")
ax.set_xlim(0,1.5); ax.set_title("Rare-subset identity programme is diluted in floor",loc="left")
fig.legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,label="real B (same label-transfer route)",ms=7),
                    Line2D([],[],marker="o",ls="",color=BLUE,label="PBISC floor",ms=7)],
           frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.02),ncol=2,columnspacing=2.0,handletextpad=.5)
fig.tight_layout(rect=[0,0.055,1,1]); save(fig,"rare_marker")

# ---- FIG 9: power-matched divergence ----
fig,ax=plt.subplots(1,2,figsize=(10.2,3.7))
v=[P.nD_AA.dropna(),P.nD_realB,P.nD_floor]; labs=["real A internal\n(technical floor)","real A vs real B\n(patient + transfer)","real A vs floor\n(model)"]
bp=ax[0].boxplot(v,widths=.55,patch_artist=True,showfliers=False,medianprops=dict(color=INK,lw=1.4))
for p,c in zip(bp["boxes"],["#E4E8EC","#D5DBE1","#A8CDE3"]): p.set_facecolor(c); p.set_edgecolor(GREY)
for i,vv in enumerate(v): ax[0].scatter(np.random.default_rng(1).normal(i+1,0.06,len(vv)),vv,s=14,color=INK,alpha=.5,zorder=3)
for i,vv in enumerate(v): ax[0].text(i+1,np.nanmedian(vv)+6,f"median {np.nanmedian(vv):.0f}",ha="center",fontsize=8,color=INK)
ax[0].set_xticklabels(labs,fontsize=8); ax[0].set_ylabel("genes with |Δ mean expression| ≥ 0.5")
ax[0].set_title("Within-subset divergence, equal n per source",loc="left")
mm=[P.r_AA.dropna(),P.r_realB,P.r_floor]
bp=ax[1].boxplot(mm,widths=.55,patch_artist=True,showfliers=False,medianprops=dict(color=INK,lw=1.4))
for p,c in zip(bp["boxes"],["#E4E8EC","#D5DBE1","#A8CDE3"]): p.set_facecolor(c); p.set_edgecolor(GREY)
for i,vv in enumerate(mm): ax[1].scatter(np.random.default_rng(2).normal(i+1,0.06,len(vv)),vv,s=14,color=INK,alpha=.5,zorder=3)
ax[1].set_xticklabels(labs,fontsize=8); ax[1].set_ylabel("profile Pearson r"); ax[1].set_ylim(0.93,1.0)
ax[1].set_title("Profile correlation hides it",loc="left")
fig.tight_layout(); save(fig,"rare_power")

# ---- FIG 10: under-production ----
d=RA[["cohort","subset","n_refA","n_floor_tr"]].copy()
d["ratio"]=d.n_floor_tr/d.n_refA; d=d.sort_values("ratio")
fig,ax=plt.subplots(figsize=(8.0,4.9))
y=np.arange(len(d)); cols=[RED if r<0.3 else (ORANGE if r<0.6 else BLUE) for r in d.ratio]
ax.barh(y,d.ratio,0.6,color=[c+"55" for c in cols],edgecolor=cols,lw=1.0)
ax.axvline(1,ls="--",lw=.9,c=INK)
for i,r in enumerate(d.itertuples()):
    ax.text(min(r.ratio,2.6)+0.03,i,f"{r.n_floor_tr} / {r.n_refA}",va="center",fontsize=7.5,color=INK)
ax.set_yticks(y); ax.set_yticklabels([f"{r.subset}  ({r.cohort})" for r in d.itertuples()],fontsize=8)
ax.set_xlabel("floor cells / real A cells with the same subset label  (60k-cell balanced sets)")
ax.set_xlim(0,3.0); ax.set_title("Rare-subset yield: five subsets fall below the 30-cell evaluation threshold",loc="left")
fig.tight_layout(); save(fig,"rare_counts")

# ---- FIG 11: driver genes heatmap ----
G=pd.read_csv("out/rare_gene_diff.csv.gz")
SEL={"B_c07-MZB1":["JCHAIN","MZB1","IGLL5","DERL3","CD27","TXNDC5"],
     "B_c08-MKI67":["JCHAIN","MZB1","TNFRSF17","S100A8","IFITM1","JUN"],
     "pDC":["HLA-DRB5","HLA-DQA1","HIST1H1E","HIST1H1D","HIST1H1C","JUNB"],
     "MAIT":["KLRB1","NCR3","DUSP2","HIST1H1E","RPS26","GAPDH"],
     "cDC":["S100A8","S100A9","IL1B","CD14","HLA-DPA1","CD74"],
     "NK CD56 bright":["AREG","DUSP2","GZMB","HIST1H1C","JUN","RPS26"]}
fig,ax=plt.subplots(1,len(SEL),figsize=(12.4,3.4),sharey=False)
for j,(sub,gs) in enumerate(SEL.items()):
    d=G[G.subset==sub].set_index("gene")
    gs=[g for g in gs if g in d.index]
    yy=np.arange(len(gs))
    ax[j].barh(yy+0.19,[d.loc[g,"d_realB"] for g in gs],0.36,color="#D5DBE1",edgecolor=GREY,lw=.6)
    ax[j].barh(yy-0.19,[d.loc[g,"d_floor"] for g in gs],0.36,color="#A8CDE3",edgecolor=BLUE,lw=.6)
    ax[j].axvline(0,lw=.8,c=INK); ax[j].set_yticks(yy); ax[j].set_yticklabels(gs,fontsize=7.5)
    ax[j].set_title(sub,loc="left",fontsize=9); ax[j].set_xlabel("Δ vs real A",fontsize=8)
    ax[j].invert_yaxis()
fig.legend(handles=[plt.Rectangle((0,0),1,1,fc="#D5DBE1",ec=GREY),plt.Rectangle((0,0),1,1,fc="#A8CDE3",ec=BLUE)],
           labels=["real B","PBISC floor"],frameon=False,fontsize=8,loc="lower center",
           bbox_to_anchor=(0.5,-0.015),ncol=2,columnspacing=2.0,handletextpad=.5)
fig.suptitle("Inside the same subset label: which genes move",x=0.005,ha="left",fontsize=10,fontweight="bold")
fig.tight_layout(rect=[0,0.07,1,0.97]); save(fig,"rare_drivers")
