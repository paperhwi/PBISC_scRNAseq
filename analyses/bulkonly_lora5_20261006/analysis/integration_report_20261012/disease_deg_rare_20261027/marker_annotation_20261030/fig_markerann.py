import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"font.size":9,"axes.titlesize":10,"legend.fontsize":8})
BLUE,RED,GREY,INK,LBLUE,LGREY="#2B83BA","#D7191C","#7F7F7F","#1F2933","#9ECAE1","#C7CDD3"
P=pd.read_csv("markerann/proportions_long.csv"); E=pd.read_csv("markerann/expression_similarity_long.csv")
F=pd.read_csv("markerann/realB_per_subset_F1.csv"); F=F[F.label=="fine_mk_cl"].set_index(["cohort","subset"]).F1
V=pd.read_csv("markerann/realB_marker_accuracy.csv")
for D in (P,E): D["res"]=[F.get((c,s),0)>=0.5 for c,s in zip(D.cohort,D.subset)]
P["l2r"]=np.log2((P.pct+0.01)/(P.pct_realA+0.01))
C=[("real B · author","real B\nauthor",LGREY),("real B · marker (cluster)","real B\nmarker",GREY),
   ("floor · kNN transfer","floor\nkNN",LBLUE),("floor · marker (cluster)","floor\nmarker",BLUE)]
fig,ax=plt.subplots(1,3,figsize=(12.4,3.9),gridspec_kw={"width_ratios":[0.8,1.2,1.2]})
# (a) ruler validity
acc={}
for coh in ["sepsis","covid"]:
    L=pd.read_csv(f"markerann/{coh}_realB_thin_labels.csv.gz",index_col=0)
    acc[coh]=[(L.six_type.astype(str)==L.six_mk_cl).mean(),V[(V.cohort==coh)&(V.label=="fine_mk_cl")].accuracy.iloc[0],
              V[(V.cohort==coh)&(V.label=="fine_mk_cl")].macroF1.iloc[0]]
x=np.arange(3); w=.36
for k,(coh,col) in enumerate([("sepsis",GREY),("covid",INK)]):
    b=ax[0].bar(x+(k-.5)*w,acc[coh],w,color=col,label=coh)
    for xi,v in zip(x+(k-.5)*w,acc[coh]): ax[0].text(xi,v+.015,f"{v:.2f}",ha="center",fontsize=7)
ax[0].set_xticks(x); ax[0].set_xticklabels(["6-type\naccuracy","fine\naccuracy","fine\nmacro-F1"]); ax[0].set_ylim(0,1.05)
ax[0].set_title("a  Top10-marker ruler on real B\n(vs author labels)",loc="left"); ax[0].legend(frameon=False,loc="upper right")
# (b) proportion error, resolvable subsets
for j,(metric,axx,yl,ttl) in enumerate([("l2r",ax[1],"median |log2 (share / real A share)|","b  subset share error"),
                                        ("r_identity_nonmarker",ax[2],"median identity r (non-marker HVGs)","c  expression identity (non-marker genes)")]):
    D=P if metric=="l2r" else E
    for k,coh in enumerate(["sepsis","covid"]):
        for i,(c,lab,col) in enumerate(C):
            d=D[(D.cohort==coh)&(D.combo==c)&D.res][metric].dropna()
            v=np.median(np.abs(d)) if metric=="l2r" else np.median(d)
            xi=k*5+i; axx.bar(xi,v,0.8,color=col,edgecolor="black" if "marker" in c else "none",lw=.6)
            axx.text(xi,v+(0.03 if metric=="l2r" else 0.005),f"{v:.2f}",ha="center",fontsize=7)
    axx.set_xticks([k*5+i for k in range(2) for i in range(4)]); axx.set_xticklabels([l for _ in range(2) for _,l,_ in C],fontsize=7)
    for k,coh in enumerate(["sepsis","covid"]):
        n=int(D[(D.cohort==coh)&(D.combo==C[0][0])&D.res].shape[0])
        axx.text(k*5+1.5,-0.24,f"{coh} (n={n})",transform=axx.get_xaxis_transform(),ha="center",fontsize=8,fontweight="bold")
    axx.set_ylabel(yl); axx.set_title(ttl+"\nmarker-resolvable subsets (real B F1 ≥ 0.5)",loc="left")
    if metric!="l2r": axx.set_ylim(0.8,1.0)
fig.tight_layout(); fig.savefig("markerann/fig1_ruler_share_identity.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig)

# fig 2: key subsets share dot plot
K=pd.read_csv("markerann/key_subsets_table.csv")
fig,ax=plt.subplots(1,2,figsize=(12.4,4.6))
S=[("pct_realA","real A (author)","o",INK,"none"),("pct_realB_author","real B · author","o",LGREY,"none"),("pct_realB_mk","real B · marker","s",GREY,"black"),
   ("pct_floor_knn","floor · kNN","^",LBLUE,"none"),("pct_floor_mk","floor · marker","^",BLUE,"black")]
for j,coh in enumerate(["sepsis","covid"]):
    k=K[K.cohort==coh].reset_index(drop=True); y=np.arange(len(k))[::-1]
    for i,(c,lab,mk,col,ec) in enumerate(S):
        v=k[c].values.astype(float); v=np.where(v<=0,0.004,v)
        ax[j].scatter(v,y+(i-2)*0.13,marker=mk,s=34,c=col,edgecolors=ec,linewidths=.6,label=lab,zorder=3)
    for yy in y: ax[j].axhline(yy,color="#EEF0F2",lw=.6,zorder=0)
    ax[j].set_xscale("log"); ax[j].set_xlim(0.003,12); ax[j].axvline(0.004,color=GREY,ls=":",lw=.7)
    ax[j].text(0.0042,y.max()+0.6,"0 cells",fontsize=6.5,color=GREY)
    ax[j].set_yticks(y); ax[j].set_yticklabels([f"{s}  (F1 {f:.2f})" for s,f in zip(k.subset,k.F1_realB)],fontsize=8)
    ax[j].set_xlabel("% of all cells (log)"); ax[j].set_title(f"{coh}: key subsets (F1 = marker ruler on real B)",loc="left")
ax[1].legend(frameon=False,loc="upper center",bbox_to_anchor=(-0.15,-0.13),ncol=5)
fig.tight_layout(); fig.savefig("markerann/fig2_key_subsets_share.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig)
print("ok")
