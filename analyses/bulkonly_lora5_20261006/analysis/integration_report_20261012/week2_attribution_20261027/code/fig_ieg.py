import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,ORANGE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#E08214"
F="figs/"; save=lambda fig,n:(fig.savefig(F+n+".png",dpi=220,bbox_inches="tight",facecolor="white"),plt.close(fig),print("saved",n))
S=pd.read_csv("out2/ieg_split.csv"); Q=pd.read_csv("out2/ieg_score_contrast.csv")
CT=["T cell","B cell","NK cell","Monocyte","DC"]
lab=[f"{c[:3]} {t.replace(' cell','')}" for c in ["sepsis","covid"] for t in CT]
S=S.set_index(["cohort","celltype"]).loc[[(c,t) for c in ["sepsis","covid"] for t in CT]].reset_index()
fig,ax=plt.subplots(1,3,figsize=(12.0,3.9))
x=np.arange(len(S))
# a: r before/after
ax[0].vlines(x,S.r_noIEG,S.r_all,color="#D5DBE1",lw=2.4)
ax[0].plot(x,S.r_all,"o",color=GREY,ms=6,label="all genes")
ax[0].plot(x,S.r_noIEG,"o",color=BLUE,ms=6,label="IEG + ambient removed")
for i,d in enumerate(S.delta_r):
    if abs(d)>=0.03: ax[0].text(i,S.r_noIEG[i]-0.055,f"{d:+.02f}",ha="center",fontsize=7,color=RED)
ax[0].set_xticks(x); ax[0].set_xticklabels(lab,rotation=55,ha="right"); ax[0].set_ylim(0,0.85)
ax[0].set_ylabel("disease logFC r vs real"); ax[0].set_title("The correlation is not an IEG artefact",loc="left")
ax[0].legend(frameon=False,loc="lower left")
# b: share of DEGs that are IEG
w=0.38
ax[1].bar(x-w/2,S.pct_realDEG_IEG,w,color="#F7D9AF",edgecolor=ORANGE,lw=.9,label="% of real DEGs")
ax[1].bar(x+w/2,S.pct_top30_IEG,w,color="#F2A7A8",edgecolor=RED,lw=.9,label="% of real top-30 |logFC|")
ax[1].set_xticks(x); ax[1].set_xticklabels(lab,rotation=55,ha="right")
ax[1].set_ylabel("% immediate-early"); ax[1].set_ylim(0,72)
ax[1].set_title("But it dominates the headline gene lists",loc="left"); ax[1].legend(frameon=False,loc="upper right")
# c: IEG score effect size
P=Q.pivot_table(index=["cohort","celltype"],columns="source",values="cohens_d").loc[[(c,t) for c in ["sepsis","covid"] for t in CT]]
ax[2].bar(x-w/2,P["real"],w,color="#D5DBE1",edgecolor=GREY,lw=.9,label="real")
ax[2].bar(x+w/2,P["floor"],w,color="#A8CDE3",edgecolor=BLUE,lw=.9,label="floor")
ax[2].axhline(0,lw=.9,c=INK)
ax[2].set_xticks(x); ax[2].set_xticklabels(lab,rotation=55,ha="right")
ax[2].set_ylabel("IEG score, case vs control (Cohen's d)")
ax[2].set_title("The IEG axis flips sign between cohorts",loc="left"); ax[2].legend(frameon=False,loc="lower right")
fig.suptitle("Immediate-early / handling programme: separated",x=0.005,ha="left",fontsize=10.5,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.95]); save(fig,"w3_ieg")
