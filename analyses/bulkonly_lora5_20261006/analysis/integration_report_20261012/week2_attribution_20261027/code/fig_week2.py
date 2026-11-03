import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,GREEN,ORANGE,PURPLE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#1A9641","#E08214","#7B3294"
F="figs/"; save=lambda fig,n:(fig.savefig(F+n+".png",dpi=220,bbox_inches="tight",facecolor="white"),plt.close(fig),print("saved",n))

# ---- A1 ladder ----
L=pd.read_csv("out2/a1_ladder.csv")
order=["1 real A (author)","2 real B (other patients)","3 real B depth-matched",
       "4 real B + read-out (Poisson)","5 real B + read-out (floor rule)",
       "6 PBISC floor","7 PBISC prev K4000","8 scSemiProfiler"]
mods=["C1Q / complement","pDC / cDC identity","plasma-cell secretory","sepsis myeloid activation","ISG (COVID)","IFN / IFI27 (COVID)"]
P=L.pivot_table(index="source",columns="module",values="overall_ratio").reindex(order)[mods]
fig,ax=plt.subplots(figsize=(10.0,4.4))
x=np.arange(len(order)); cols=[RED,ORANGE,PURPLE,GREEN,BLUE,"#5FA8D3"]
for j,m in enumerate(mods):
    ax.plot(x,P[m],"-o",color=cols[j],ms=5,lw=1.6,label=m)
ax.axhline(1,ls="--",lw=.9,c=INK); ax.axhline(0.5,ls=":",lw=.9,c=GREY)
ax.axvspan(2.5,4.5,color="#EDF2F6",zorder=0)
ax.text(3.5,2.85,"read-out emulation\n(real cells)",ha="center",fontsize=8,color=GREY)
ax.axvspan(4.5,7.5,color="#FBF2E6",zorder=0); ax.text(6,2.85,"generative models",ha="center",fontsize=8,color=GREY)
ax.set_xticks(x); ax.set_xticklabels([s.split(" ",1)[1] for s in order],rotation=28,ha="right",fontsize=8)
ax.set_ylabel("module recovery  (vs depth-matched real B)"); ax.set_ylim(0,3.15)
ax.set_title("Control ladder: depth → read-out rule → model",loc="left")
fig.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.05),ncol=3,columnspacing=1.6)
fig.tight_layout(rect=[0,0.10,1,1]); save(fig,"w2_ladder")

# ---- C1Q decomposition waterfall ----
c=P["C1Q / complement"]
steps=[("depth-matched real",1.0),("+ read-out rule",c["5 real B + read-out (floor rule)"]),("PBISC floor",c["6 PBISC floor"])]
fig,ax=plt.subplots(figsize=(5.6,3.6))
v=[s[1] for s in steps]
ax.bar(range(3),v,0.55,color=["#D5DBE1","#F7D9AF","#A8CDE3"],edgecolor=[GREY,ORANGE,BLUE],lw=1.1)
for i,y in enumerate(v): ax.text(i,y+0.02,f"{y:.2f}",ha="center",fontsize=9)
ax.annotate("",xy=(1,v[1]+0.08),xytext=(0,v[0]+0.08),arrowprops=dict(arrowstyle="-|>",color=ORANGE,lw=1.2))
ax.text(0.5,v[0]+0.12,f"read-out rule\n−{100*(1-v[1]):.0f}%",ha="center",fontsize=8,color=ORANGE)
ax.annotate("",xy=(2,v[2]+0.08),xytext=(1,v[1]+0.08),arrowprops=dict(arrowstyle="-|>",color=BLUE,lw=1.2))
ax.text(1.5,v[1]+0.12,f"model\n−{100*(1-v[2]/v[1]):.0f}%",ha="center",fontsize=8,color=BLUE)
ax.set_xticks(range(3)); ax.set_xticklabels([s[0] for s in steps],fontsize=8.5)
ax.set_ylabel("C1Q module recovery"); ax.set_ylim(0,1.35)
ax.set_title("C1Q loss: ~1/3 read-out rule, ~4/5 of the rest is the model",loc="left",fontsize=9.5)
fig.tight_layout(); save(fig,"w2_c1q_decomp")

# ---- B2 three models ----
T=P.loc[["6 PBISC floor","7 PBISC prev K4000","8 scSemiProfiler"]]
fig,ax=plt.subplots(figsize=(9.4,3.9))
x=np.arange(len(mods)); w=0.26
for k,(s,col,ec) in enumerate([("6 PBISC floor","#A8CDE3",BLUE),("7 PBISC prev K4000","#F2A7A8",RED),("8 scSemiProfiler","#D6C2E3",PURPLE)]):
    ax.bar(x+(k-1)*w,T.loc[s],w,color=col,edgecolor=ec,lw=1.0,label=s.split(" ",1)[1])
ax.axhline(1,ls="--",lw=.9,c=INK); ax.axhline(0.5,ls=":",lw=.9,c=GREY)
ax.set_xticks(x); ax.set_xticklabels([m.replace(" / ","/\n").replace(" secretory","\nsecretory").replace(" activation","\nactivation").replace(" (COVID)","\n(COVID)") for m in mods],fontsize=7.5)
ax.set_ylabel("module recovery vs depth-matched real"); ax.set_ylim(0,2.1)
ax.set_title("Three models, same metric",loc="left")
fig.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.03),ncol=3,columnspacing=2.2)
fig.tight_layout(rect=[0,0.09,1,1]); save(fig,"w2_three_models")

# ---- B1 specificity law ----
G=pd.read_csv("out2/b1b_gene_level_law.csv.gz"); B=pd.read_csv("out2/b1_datadriven_modules.csv")
fig,ax=plt.subplots(1,2,figsize=(10.6,4.0))
G["rank_band"]=pd.cut(G.rank_in_type,[0,500,1500,3000,99999],labels=["rank 1-500","rank 501-1500","rank 1501-3000",">3000"])
G["tau_q"]=pd.qcut(G.tau,5,labels=["Q1\nbroad","Q2","Q3","Q4","Q5\nspecific"])
pt=G.pivot_table(index="tau_q",columns="rank_band",values="floor_ratio",aggfunc="median",observed=True)
for j,c in enumerate(pt.columns[:3]):
    ax[0].plot(range(5),pt[c],"-o",ms=5,lw=1.6,color=[INK,BLUE,GREY][j],label=c)
ax[0].axhline(1,ls="--",lw=.9,c=INK); ax[0].set_xticks(range(5)); ax[0].set_xticklabels(pt.index,fontsize=8)
ax[0].set_ylabel("median gene recovery (floor / depth-matched real)")
ax[0].set_xlabel("cell-type specificity quintile (Tau, defined on real half A)")
rho=stats.spearmanr(G.tau,G.floor_ratio)
ax[0].set_title(f"Gene level, n={len(G):,}  ρ = {rho.statistic:.2f}",loc="left")
ax[0].legend(frameon=False,loc="lower left")
sc=ax[1].scatter(B.tau,B.floor_ratio,s=26+B.n_genes/4,c=[BLUE if c=="sepsis" else ORANGE for c in B.cohort],alpha=.75,lw=0)
m,b=np.polyfit(B.tau,B.floor_ratio,1); xs=np.linspace(B.tau.min(),B.tau.max(),20)
ax[1].plot(xs,m*xs+b,lw=1.3,c=RED); ax[1].axhline(1,ls="--",lw=.9,c=INK)
r2=stats.spearmanr(B.tau,B.floor_ratio)
ax[1].set_xlabel("module specificity (mean Tau)"); ax[1].set_ylabel("module recovery")
ax[1].set_title(f"Data-driven modules (k=30/cohort), n={len(B)}  ρ = {r2.statistic:.2f}, p = {r2.pvalue:.1e}",loc="left",fontsize=9)
ax[1].legend(handles=[Line2D([],[],marker='o',ls='',color=BLUE,label='sepsis',ms=6),
                      Line2D([],[],marker='o',ls='',color=ORANGE,label='COVID',ms=6)],frameon=False,loc="lower left")
fig.suptitle("Recovery falls as a programme becomes cell-type-exclusive",x=0.005,ha="left",fontsize=10.5,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.95]); save(fig,"w2_specificity_law")
