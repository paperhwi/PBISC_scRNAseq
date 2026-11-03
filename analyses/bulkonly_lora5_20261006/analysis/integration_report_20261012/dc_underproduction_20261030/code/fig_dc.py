import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"font.size":9,"axes.titlesize":10,"legend.fontsize":8})
BLUE,RED,GREY,INK,LBLUE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#9ECAE1"
G=pd.read_csv("dc/dc1_fraction_by_group.csv"); T=pd.read_csv("dc/dc2_bulk_tracking.csv"); C=pd.read_csv("dc/dc3_cell_level.csv")
fig,ax=plt.subplots(1,3,figsize=(12.6,3.8),gridspec_kw={"width_ratios":[1,1.35,1]})
# a: sepsis DC% by group
g=G[G.cohort=="sepsis"].pivot(index="group",columns="source",values="DC_pct_median").reindex(["HC","ICU","case"])
x=np.arange(3); w=.38
ax[0].bar(x-w/2,g.real,w,color=GREY,label="real (all samples)"); ax[0].bar(x+w/2,g.floor,w,color=BLUE,label="PBISC floor")
for i,k in enumerate(g.index):
    r=g.loc[k,"floor"]/g.loc[k,"real"]; ax[0].text(i+w/2,g.loc[k,"floor"]+0.03,f"{r:.2f}×",ha="center",fontsize=8,color=RED if r<0.7 else INK)
ax[0].set_xticks(x); ax[0].set_xticklabels(["adult HC","adult ICU","adult sepsis"]); ax[0].set_ylabel("DC % of cells (sample median)")
ax[0].set_title("a  sepsis: DC deficit is in the controls",loc="left"); ax[0].legend(frameon=False,loc="upper right")
# b: bulk tracking r
lin=["pDC","cDC","monocyte","B cell","T cell","NK"]
def get(coh,src):
    d=T[(T.cohort==coh)&(T.source==src)].copy(); d["lin"]=d.markers.str.split(" \\(").str[0]
    return d.set_index("lin").r.reindex(lin)
x=np.arange(len(lin)); w=.26
for k,(coh,src,col,lab) in enumerate([("sepsis","real",GREY,"sepsis real"),("sepsis","floor",BLUE,"sepsis floor"),("covid","floor",LBLUE,"COVID floor")]):
    v=get(coh,src); ax[1].bar(x+(k-1)*w,v,w,color=col,label=lab)
    for xi,vi in zip(x+(k-1)*w,v): ax[1].text(xi,(vi if vi>0 else 0)+0.02,f"{vi:.2f}",ha="center",fontsize=6.5,rotation=90)
ax[1].axhline(0,color=INK,lw=.6); ax[1].set_xticks(x); ax[1].set_xticklabels(lin); ax[1].set_ylim(-0.3,1.15)
ax[1].set_ylabel("r (bulk marker score vs\ngenerated pseudobulk, same sample)"); ax[1].set_title("b  does the output follow its own bulk input?",loc="left")
ax[1].legend(frameon=False,loc="upper left",ncol=3,bbox_to_anchor=(0,1.0))
# c: pDC marker co-expression
c=C.set_index(["cohort","source"]).coexpr_given_any_pDC
srcs=["real A","real B depth-matched","PBISC floor","PBISC prev K4000"]; labs=["real A","real B","floor","prev K4000"]; cols=[INK,GREY,BLUE,RED]
x=np.arange(2); w=.2
for k,(s,l,col) in enumerate(zip(srcs,labs,cols)):
    v=[c.loc[("sepsis",s)],c.loc[("covid",s)]]; ax[2].bar(x+(k-1.5)*w,v,w,color=col,label=l)
    for xi,vi in zip(x+(k-1.5)*w,v): ax[2].text(xi,vi+0.6,f"{vi:.0f}",ha="center",fontsize=7)
ax[2].set_xticks(x); ax[2].set_xticklabels(["sepsis","COVID"]); ax[2].set_ylabel("% of pDC-marker+ cells co-expressing ≥2\n(LILRA4 / CLEC4C / PTCRA …)")
ax[2].set_title("c  pDC identity is not assembled",loc="left"); ax[2].legend(frameon=False,loc="upper left",ncol=2)
ax[2].set_ylim(0,max(c.values)*1.35)
fig.tight_layout(); fig.savefig("figs/v4_dc_underproduction.png",dpi=220,bbox_inches="tight",facecolor="white"); print("ok")
