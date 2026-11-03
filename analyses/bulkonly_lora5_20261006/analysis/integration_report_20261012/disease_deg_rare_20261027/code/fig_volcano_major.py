"""Volcano plots for the major lymphoid populations (T, B, NK), real vs PBISC floor, both cohorts.
Same layout as v3 slide 10: top row real, bottom row floor, shared x-axis per column, highlighted genes circled in red.
Immediate-early / handling genes (the IEG set from ieg_split.py) are coloured orange when called DEG."""
import numpy as np, pandas as pd, pickle, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8.5,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,ORANGE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#E08214"
IEG=set(["FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","EGR4","NR4A1","NR4A2","NR4A3","ATF3","IER2","IER3","IER5",
 "DUSP1","DUSP2","DUSP5","KLF2","KLF4","KLF6","KLF10","ZFP36","ZFP36L1","ZFP36L2","CSRNP1","RGS1","RGS2","RGS16","PPP1R15A",
 "BTG1","BTG2","MYC","NFKBIA","NFKBIZ","TNFAIP3","CD69","CD83","RASD1","PHLDA1","TRIB1","PER1","HES1","MAFF","SGK1","SOCS3",
 "SLC7A5","HSPA1A","HSPA1B","DNAJB1","HSPH1","ID1","ID2","ID3","CCNL1","MIDN","SIK1","PLK2","GADD45B","CEBPD","ARC","NPAS4",
 "HBB","HBA1","HBA2","HBD","PPBP","PF4","NRGN","ALAS2"])
allde=pickle.load(open("out/allde.pkl","rb"))
HL={("sepsis","T cell"):["NR4A1","AREG","SGK1","CDKN1A","HBB","ZNF683","SLC4A10","S100A8"],
    ("sepsis","B cell"):["NR4A1","SOCS3","TNF","LILRB4","JCHAIN","MZB1","S100A8"],
    ("sepsis","NK cell"):["NR4A1","AREG","HBB","LMNA","ZNF683","ZMAT4","CD83","S100A8"],
    ("covid","T cell"):["IFI27","IFI44L","IFIT1","RSAD2","MX1","CISH","S100A8","S100A9"],
    ("covid","B cell"):["IFIT3","IFIT1","IFI44L","OAS3","USP18","S100A8","S100A9","HES1"],
    ("covid","NK cell"):["IFIT1","IFI44L","IFI27","RSAD2","MX1","KRT86","S100A8","S100A9"]}
TITLE={"sepsis":"Sepsis — adult sepsis vs healthy control","covid":"COVID / ESKD — positive vs negative (patient level)"}

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
    n_ieg=int((sig&ieg).sum()); n=int(sig.sum())
    ax.set_title(f"{title}  (DEG n={n}; IEG {n_ieg})",loc="left",fontsize=9.5)
    ax.set_xlabel("log2 FC (case vs control)"); ax.set_ylabel("-log10 FDR")
    return txt

for coh in ["sepsis","covid"]:
    fig,ax=plt.subplots(2,3,figsize=(11.6,6.6))
    alltxt=[]
    for j,ct in enumerate(["T cell","B cell","NK cell"]):
        t=allde[(coh,ct)]
        alltxt.append((ax[0,j],volc(ax[0,j],t["real"],f"{ct} — real",HL[(coh,ct)],GREY)))
        alltxt.append((ax[1,j],volc(ax[1,j],t["floor"],f"{ct} — PBISC floor",HL[(coh,ct)],BLUE)))
        lo=min(ax[0,j].get_xlim()[0],ax[1,j].get_xlim()[0]); hi=max(ax[0,j].get_xlim()[1],ax[1,j].get_xlim()[1])
        ax[0,j].set_xlim(lo,hi); ax[1,j].set_xlim(lo,hi)
    fig.suptitle(TITLE[coh],x=0.005,ha="left",fontsize=11,fontweight="bold")
    fig.legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,ms=5,label="DEG, real"),
                        Line2D([],[],marker="o",ls="",color=BLUE,ms=5,label="DEG, PBISC floor"),
                        Line2D([],[],marker="o",ls="",color=ORANGE,ms=5,label="DEG, immediate-early / ambient gene"),
                        Line2D([],[],marker="o",ls="",color="#D5DBE1",ms=5,label="not DEG"),
                        Line2D([],[],marker="o",ls="",mfc="none",mec=RED,ms=6,label="highlighted gene")],
               frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.01),ncol=5,columnspacing=1.6,handletextpad=.4)
    fig.tight_layout(rect=[0,0.05,1,0.96])
    for a,tx in alltxt:
        if tx: adjust_text(tx,ax=a,expand=(1.8,2.2),force_text=(0.8,1.2),force_static=(0.4,0.6),max_move=60,arrowprops=dict(arrowstyle="-",color=RED,lw=.5))
    fig.savefig(f"figs/volcano_lymphoid_{coh}.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig)
    print("saved",coh)
