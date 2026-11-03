import numpy as np, pandas as pd, pickle, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8.5,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,ORANGE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#E08214"
import importlib.util,sys
spec=importlib.util.spec_from_file_location("fv","fig_volcano_major.py")
src=open("fig_volcano_major.py").read().split("allde=pickle.load")[0]; ns={}; exec(src,ns); IEG=ns["IEG"]
HL={("sepsis","T cell"):["NR4A1","AREG","SGK1","CDKN1A","HBB","ZNF683","SLC4A10","S100A8"],
    ("sepsis","B cell"):["NR4A1","SOCS3","TNF","LILRB4","JCHAIN","MZB1","S100A8"],
    ("sepsis","NK cell"):["NR4A1","AREG","HBB","LMNA","ZNF683","ZMAT4","CD83","S100A8"],
    ("covid","T cell"):["IFI27","IFI44L","IFIT1","RSAD2","MX1","CISH","S100A8","S100A9"],
    ("covid","B cell"):["IFIT3","IFIT1","IFI44L","OAS3","USP18","S100A8","S100A9","HES1"],
    ("covid","NK cell"):["IFIT1","IFI44L","IFI27","RSAD2","MX1","KRT86","S100A8","S100A9"]}
REP=pickle.load(open("matched/representative_tables.pkl","rb"))
S=pd.read_csv("matched/matched_n_summary.csv").set_index(["cohort","celltype","design"])
def volc(ax,T,title,hl,base):
    q=np.clip(T.q.values,1e-300,1); y=-np.log10(q); x=T.logFC.values
    sig=((T.q<0.05)&(T.logFC.abs()>=1)).values; ieg=T.gene.isin(IEG).values
    ax.scatter(x[~sig],y[~sig],s=2.5,c="#D5DBE1",lw=0,rasterized=True)
    ax.scatter(x[sig&~ieg],y[sig&~ieg],s=4.5,c=base,lw=0,rasterized=True)
    ax.scatter(x[sig&ieg],y[sig&ieg],s=6,c=ORANGE,lw=0,rasterized=True)
    ax.axhline(-np.log10(0.05),ls="--",lw=.7,c=GREY); ax.axvline(1,ls="--",lw=.7,c=GREY); ax.axvline(-1,ls="--",lw=.7,c=GREY)
    txt=[]
    for g in hl:
        r=T[T.gene==g]
        if r.empty: continue
        r=r.iloc[0]; yy=-np.log10(max(r.q,1e-300))
        ax.scatter([r.logFC],[yy],s=26,facecolor="none",edgecolor=RED,lw=1.1,zorder=4)
        txt.append(ax.text(r.logFC,yy,g,fontsize=7,color=RED,zorder=5,bbox=dict(fc="white",ec="none",pad=0.4,alpha=0.75)))
    up=int((sig&(x>0)).sum()); dn=int((sig&(x<0)).sum())
    ax.set_title(f"{title}  (up {up} / down {dn}; IEG {int((sig&ieg).sum())})",loc="left",fontsize=9.3)
    ax.set_xlabel("log2 FC (case vs control)"); ax.set_ylabel("-log10 FDR")
    return txt
for coh in ["sepsis","covid"]:
    fig,ax=plt.subplots(2,3,figsize=(11.6,6.8)); jobs=[]
    for j,ct in enumerate(["T cell","B cell","NK cell"]):
        t=REP[(coh,ct,"M")]; n=S.loc[(coh,ct,"M")]
        jobs.append((ax[0,j],volc(ax[0,j],t["real"],f"{ct} — real",HL[(coh,ct)],GREY)))
        jobs.append((ax[1,j],volc(ax[1,j],t["floor"],f"{ct} — PBISC floor",HL[(coh,ct)],BLUE)))
        lo=min(ax[0,j].get_xlim()[0],ax[1,j].get_xlim()[0]); hi=max(ax[0,j].get_xlim()[1],ax[1,j].get_xlim()[1])
        yh=max(ax[0,j].get_ylim()[1],ax[1,j].get_ylim()[1])
        for a in (ax[0,j],ax[1,j]): a.set_xlim(lo,hi); a.set_ylim(-0.3,yh)
    fig.suptitle({"sepsis":"Sepsis — matched n (real and floor: 68 case / 7 control, one representative of 50 draws)",
                  "covid":"COVID / ESKD — matched n (real and floor: 37–39 case / 18 control patients, representative of 50 draws)"}[coh],
                 x=0.005,ha="left",fontsize=10.5,fontweight="bold")
    fig.legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,ms=5,label="DEG, real"),
                        Line2D([],[],marker="o",ls="",color=BLUE,ms=5,label="DEG, PBISC floor"),
                        Line2D([],[],marker="o",ls="",color=ORANGE,ms=5,label="DEG, immediate-early / ambient gene"),
                        Line2D([],[],marker="o",ls="",color="#D5DBE1",ms=5,label="not DEG"),
                        Line2D([],[],marker="o",ls="",mfc="none",mec=RED,ms=6,label="highlighted gene")],
               frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.01),ncol=5,columnspacing=1.6,handletextpad=.4)
    fig.tight_layout(rect=[0,0.05,1,0.96])
    for a,tx in jobs:
        if tx: adjust_text(tx,ax=a,expand=(1.8,2.2),force_text=(0.8,1.2),force_static=(0.4,0.6),max_move=60,arrowprops=dict(arrowstyle="-",color=RED,lw=.5))
    fig.savefig(f"figs/volcano_lymphoid_{coh}_matchedN.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig); print("saved",coh)

# ---- summary: up / down DEG counts over 50 draws, and the variance ratio ----
D=pd.read_csv("matched/matched_n_reps.csv"); V=pd.read_csv("matched/variance_check.csv")
fig,ax=plt.subplots(1,2,figsize=(11.4,3.9),gridspec_kw=dict(width_ratios=[1.6,1]))
d=D[D.design=="M"]; keys=[(c,t) for c in ["sepsis","covid"] for t in ["T cell","B cell","NK cell"]]
x=np.arange(len(keys)); w=0.2
for k,(col,lab,fc,ec) in enumerate([("real_up","real up","#D5DBE1",GREY),("floor_up","floor up","#A8CDE3",BLUE),
                                     ("real_down","real down","#ECEFF2",GREY),("floor_down","floor down","#DCEBF5",BLUE)]):
    m=[d[(d.cohort==c)&(d.celltype==t)][col].mean() for c,t in keys]; sd=[d[(d.cohort==c)&(d.celltype==t)][col].std() for c,t in keys]
    sign=1 if "up" in col else -1
    ax[0].bar(x+(k%2-0.5)*w*1.1,[sign*v for v in m],w,color=fc,edgecolor=ec,lw=.9,hatch=None if "up" in col else "///",label=lab)
    ax[0].errorbar(x+(k%2-0.5)*w*1.1,[sign*v for v in m],yerr=sd,fmt="none",ecolor=INK,capsize=2.5,lw=.8)
ax[0].axhline(0,lw=.9,c=INK); ax[0].set_xticks(x); ax[0].set_xticklabels([f"{c[:3]} {t.replace(' cell','')}" for c,t in keys])
ax[0].set_ylabel("DEG count (up ↑ / down ↓), mean ± SD of 50 draws"); ax[0].set_title("Matched n: up-DEGs kept or inflated, down-DEGs lost",loc="left")
v=V[V.celltype.isin(["T cell","B cell","NK cell","Monocyte"])]
lab=[f"{r.cohort[:3]} {r.celltype.replace(' cell','')} {r.group}" for r in v.itertuples()]
ax[1].barh(np.arange(len(v)),v.median_SD_ratio,0.6,color=["#A8CDE3" if g=="case" else "#E4E8EC" for g in v.group],edgecolor=BLUE,lw=.8)
ax[1].axvline(1,ls="--",lw=.9,c=INK); ax[1].set_yticks(np.arange(len(v))); ax[1].set_yticklabels(lab,fontsize=7.2); ax[1].invert_yaxis()
ax[1].set_xlabel("between-sample SD, floor / real (median gene)"); ax[1].set_title("Generated samples are too alike",loc="left")
h=[plt.Rectangle((0,0),1,1,fc="#D5DBE1",ec=GREY),plt.Rectangle((0,0),1,1,fc="#A8CDE3",ec=BLUE),
   plt.Rectangle((0,0),1,1,fc="#ECEFF2",ec=GREY,hatch="///"),plt.Rectangle((0,0),1,1,fc="#DCEBF5",ec=BLUE,hatch="///")]
fig.legend(handles=h,labels=["real up","floor up","real down","floor down"],frameon=False,loc="lower center",bbox_to_anchor=(0.33,-0.03),ncol=4)
fig.tight_layout(rect=[0,0.07,1,1]); fig.savefig("figs/matchedN_summary.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig); print("saved summary")
