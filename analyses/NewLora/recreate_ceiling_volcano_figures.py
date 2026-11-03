"""Recreate floor-era highlighted DEG volcano panels for the new LoRA output."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/results/deg_tables"
NEW=ROOT/"NewLora/analysis_20261102/disease_deg_matched/deg_tables"
OUT=ROOT/"NewLora/analysis_20261102/figure_recreation"
OUT.mkdir(parents=True,exist_ok=True)

plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,
                     "axes.spines.right":False,"font.size":9,"axes.labelsize":9,
                     "axes.titlesize":10,"legend.fontsize":8,"pdf.fonttype":42})
BLUE,RED,GREY,ORANGE="#2B83BA","#D7191C","#7F7F7F","#E08214"
IEG=set("FOS FOSB JUN JUNB JUND EGR1 EGR2 EGR3 EGR4 NR4A1 NR4A2 NR4A3 ATF3 IER2 IER3 IER5 "
        "DUSP1 DUSP2 DUSP5 KLF2 KLF4 KLF6 KLF10 ZFP36 ZFP36L1 ZFP36L2 CSRNP1 RGS1 RGS2 RGS16 "
        "PPP1R15A BTG1 BTG2 MYC NFKBIA NFKBIZ TNFAIP3 CD69 CD83 RASD1 PHLDA1 TRIB1 PER1 "
        "HES1 MAFF SGK1 SOCS3 SLC7A5 HSPA1A HSPA1B DNAJB1 HSPH1 ID1 ID2 ID3 CCNL1 MIDN SIK1 "
        "PLK2 GADD45B CEBPD ARC NPAS4 HBB HBA1 HBA2 HBD PPBP PF4 NRGN ALAS2".split())
LYMPH_HL={
    ("sepsis","T cell"):["NR4A1","AREG","SGK1","CDKN1A","HBB","ZNF683","SLC4A10","S100A8"],
    ("sepsis","B cell"):["NR4A1","SOCS3","TNF","LILRB4","JCHAIN","MZB1","S100A8"],
    ("sepsis","NK cell"):["NR4A1","AREG","HBB","LMNA","ZNF683","ZMAT4","CD83","S100A8"],
    ("covid","T cell"):["IFI27","IFI44L","IFIT1","RSAD2","MX1","CISH","S100A8","S100A9"],
    ("covid","B cell"):["IFIT3","IFIT1","IFI44L","OAS3","USP18","S100A8","S100A9","HES1"],
    ("covid","NK cell"):["IFIT1","IFI44L","IFI27","RSAD2","MX1","KRT86","S100A8","S100A9"]}
MYELOID_PANELS=[("sepsis","Monocyte",["C1QA","C1QB","C1QC","IL1R2","HLA-DRA","CLU"]),
                ("sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","IL1R2"]),
                ("covid","Monocyte",["IFI27","SIGLEC1","IFIT1","MX1"])]

def get_tab(coh,ct,source):
    name=ct.replace(" ","")
    path=(OLD if source=="real" else NEW)/f"{coh}_{name}_{source}.csv.gz"
    return pd.read_csv(path) if path.exists() else None

def volc(ax,T,title,hl,base,ieg=False,min_annot_y=0):
    if T is None:
        ax.axis("off")
        ax.text(.5,.55,"Not estimable\n0 sepsis / 0 HC patients\nwith ≥20 transferred DC cells",
                transform=ax.transAxes,ha="center",va="center",fontsize=10,color="#475569")
        ax.set_title(title,loc="left",fontsize=9.5)
        return []
    x=T.logFC.to_numpy(float)
    y=-np.log10(np.clip(T.q.to_numpy(float),1e-300,1))
    sig=(T.q.to_numpy(float)<.05)&(np.abs(x)>=1)
    is_ieg=T.gene.isin(IEG).to_numpy() if ieg else np.zeros(len(T),bool)
    ax.scatter(x[~sig],y[~sig],s=2.5,c="#D5DBE1",lw=0,rasterized=True)
    ax.scatter(x[sig&~is_ieg],y[sig&~is_ieg],s=4.5,c=base,lw=0,rasterized=True)
    if ieg:ax.scatter(x[sig&is_ieg],y[sig&is_ieg],s=6,c=ORANGE,lw=0,rasterized=True)
    ax.axhline(-np.log10(.05),ls="--",lw=.7,c=GREY)
    ax.axvline(1,ls="--",lw=.7,c=GREY)
    ax.axvline(-1,ls="--",lw=.7,c=GREY)
    texts=[]
    lookup=T.set_index("gene")
    for gene in hl:
        if gene not in lookup.index:continue
        row=lookup.loc[gene]
        xx=float(row.logFC);yy=float(-np.log10(max(row.q,1e-300)))
        if not np.isfinite(xx) or not np.isfinite(yy) or yy<min_annot_y:continue
        ax.scatter([xx],[yy],s=26,facecolor="none",edgecolor=RED,lw=1.1,zorder=4)
        texts.append(ax.text(xx,yy,gene,fontsize=7,color=RED,zorder=5,
                             bbox={"facecolor":"white","edgecolor":"none","pad":.4,"alpha":.75}))
    extra=f"; IEG {int((sig&is_ieg).sum())}" if ieg else ""
    ax.set_title(f"{title}  (DEG n={int(sig.sum())}{extra})",loc="left",fontsize=9.5)
    ax.set_xlabel("log2 FC (case vs control)")
    ax.set_ylabel("−log10 FDR")
    return texts

def adjust(items):
    for ax,txt in items:
        if txt:
            adjust_text(txt,ax=ax,expand=(1.5,1.8),force_text=(.6,.8),
                        force_static=(.3,.4),max_move=50,
                        arrowprops={"arrowstyle":"-","color":RED,"lw":.5})

def save(fig,stem):
    fig.savefig(OUT/f"{stem}.png",dpi=220,bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/f"{stem}.pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)
    print("saved",stem,flush=True)

for coh in ("sepsis","covid"):
    fig,axes=plt.subplots(2,3,figsize=(11.6,6.6))
    items=[]
    for j,ct in enumerate(("T cell","B cell","NK cell")):
        hl=LYMPH_HL[(coh,ct)]
        for i,source in enumerate(("real","ceiling")):
            T=get_tab(coh,ct,source)
            name="real" if source=="real" else "new LoRA"
            txt=volc(axes[i,j],T,f"{ct} — {name}",hl,GREY if i==0 else BLUE,
                     ieg=True,min_annot_y=1.3 if i==1 else 0)
            items.append((axes[i,j],txt))
        lo=min(axes[0,j].get_xlim()[0],axes[1,j].get_xlim()[0])
        hi=max(axes[0,j].get_xlim()[1],axes[1,j].get_xlim()[1])
        axes[0,j].set_xlim(lo,hi);axes[1,j].set_xlim(lo,hi)
    title="Sepsis — adult sepsis vs healthy control" if coh=="sepsis" else "COVID / ESKD — positive vs negative (patient level)"
    fig.suptitle(title,x=.005,ha="left",fontsize=11,fontweight="bold")
    fig.legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,ms=5,label="DEG, real"),
                        Line2D([],[],marker="o",ls="",color=BLUE,ms=5,label="DEG, new LoRA"),
                        Line2D([],[],marker="o",ls="",color=ORANGE,ms=5,label="DEG, immediate-early / ambient gene"),
                        Line2D([],[],marker="o",ls="",color="#D5DBE1",ms=5,label="not DEG"),
                        Line2D([],[],marker="o",ls="",mfc="none",mec=RED,ms=6,label="highlighted gene")],
               frameon=False,loc="lower center",bbox_to_anchor=(.5,-.01),ncol=5,
               columnspacing=1.6,handletextpad=.4)
    fig.tight_layout(rect=[0,.05,1,.96])
    adjust(items)
    save(fig,f"volcano_lymphoid_{coh}_ceiling")

fig,axes=plt.subplots(2,3,figsize=(11.6,6.4))
items=[]
for j,(coh,ct,hl) in enumerate(MYELOID_PANELS):
    for i,source in enumerate(("real","ceiling")):
        T=get_tab(coh,ct,source)
        txt=volc(axes[i,j],T,f"{coh} {ct} — {'real' if i==0 else 'new LoRA'}",hl,
                 GREY if i==0 else BLUE,min_annot_y=1.3 if i==1 else 0)
        items.append((axes[i,j],txt))
    if axes[1,j].axison:
        lo=min(axes[0,j].get_xlim()[0],axes[1,j].get_xlim()[0])
        hi=max(axes[0,j].get_xlim()[1],axes[1,j].get_xlim()[1])
        axes[0,j].set_xlim(lo,hi);axes[1,j].set_xlim(lo,hi)
fig.tight_layout()
adjust(items)
save(fig,"volcano_myeloid_dc_ceiling")
