import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"font.size":8,"axes.titlesize":8.5})
GREY,BLUE,INK,RED="#7F7F7F","#2B83BA","#1F2933","#D7191C"
U=pd.read_csv("perpatient/out/rare_unit_shares.csv"); S=pd.read_csv("perpatient/out/rare_stats.csv")
def star(q): return "***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else "ns"
for coh in ["sepsis","covid"]:
    groups=["control","ICU control","case"] if coh=="sepsis" else ["control","case"]
    gl={"control":"HC" if coh=="sepsis" else "neg","ICU control":"ICU","case":"sepsis" if coh=="sepsis" else "pos"}
    subs=list(dict.fromkeys(U[U.cohort==coh].subset)); nc=5; nr=int(np.ceil(len(subs)/nc))
    fig,axs=plt.subplots(nr,nc,figsize=(13.2,2.75*nr)); axs=axs.ravel(); rng=np.random.default_rng(0)
    for a,sub in zip(axs,subs):
        d=U[(U.cohort==coh)&(U.subset==sub)]; ymax=d.pct_all.max()
        for gi,g in enumerate(groups):
            for si,(src,col) in enumerate([("real",GREY),("floor",BLUE)]):
                v=d[(d.group==g)&(d.source==src)].pct_all.values; x=gi*2.6+si
                a.scatter(x+rng.uniform(-.28,.28,len(v)),v,s=7,color=col,alpha=.75,lw=0,zorder=2)
                if len(v)>1: a.errorbar(x+.42,v.mean(),yerr=v.std(ddof=1),fmt="o",ms=3.5,color=INK,capsize=2.5,lw=1,zorder=3)
            s=S[(S.cohort==coh)&(S.subset==sub)&(S.group==g)]
            if len(s):
                s=s.iloc[0]; a.text(gi*2.6+.5,ymax*1.13,f"{star(s.welch_q)}\nSD×{s.sd_ratio:.2f}",ha="center",va="bottom",fontsize=6.3,color=RED if s.welch_q<0.05 else INK)
        a.set_ylim(0,ymax*1.45 if ymax>0 else 1); a.set_xticks([gi*2.6+.5 for gi in range(len(groups))]); a.set_xticklabels([gl[g] for g in groups])
        a.set_title(sub,loc="left"); a.set_ylabel("% of all cells",fontsize=7)
    for a in axs[len(subs):]: a.axis("off")
    h=[Line2D([],[],marker="o",ls="",color=GREY,label="real — one dot per patient (author label)"),Line2D([],[],marker="o",ls="",color=BLUE,label="PBISC floor — one dot per generated patient (kNN label)"),
       Line2D([],[],marker="o",ls="-",color=INK,label="mean ± SD"),Line2D([],[],ls="",label="* / ** / *** = Welch t BH-q < 0.05 / 0.01 / 0.001 ; SD× = floor SD / real SD")]
    fig.legend(handles=h,loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(0.5,-0.01),fontsize=8)
    fig.suptitle(f"{'Sepsis (Ye 2026)' if coh=='sepsis' else 'COVID / ESKD (Imperial)'} — rare populations per patient",x=0.005,ha="left",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,0.05,1,0.97]); fig.savefig(f"perpatient/out/rare_dotplot_{coh}.png",dpi=200,bbox_inches="tight",facecolor="white"); plt.close(fig)
print("ok")
