"""Volcano plots for major lineage and exploratory rare-subtype DEG tables."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/results/deg_tables"
MAJOR=ROOT/"NewLora/analysis_20261102/disease_deg_matched"
RARE=ROOT/"NewLora/analysis_20261102/rare_disease_deg"
SOURCES=("real","prev4k","floor","ceiling")
NAMES={"real":"Real","prev4k":"Original PBISC","floor":"Floor","ceiling":"New LoRA"}
TYPES=("T cell","B cell","NK cell","Monocyte")
PALETTE={"up":"#C7544B","down":"#386A9F","other":"#A6AFB7"}

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,
                     "axes.spines.top":False,"axes.spines.right":False})

def volcano(ax,tab,title=None,annotate=2,xlim=6,ylim=30):
    x=tab.logFC.to_numpy(float)
    y=-np.log10(np.maximum(tab.q.to_numpy(float),1e-300))
    sig=tab.deg.astype(bool).to_numpy()
    for key,mask in (("other",~sig),("down",sig&(x<0)),("up",sig&(x>0))):
        ax.scatter(np.clip(x[mask],-xlim,xlim),np.clip(y[mask],0,ylim),
                   s=2.8 if key=="other" else 5.5,c=PALETTE[key],
                   alpha=.42 if key=="other" else .82,lw=0,rasterized=True)
    ax.axvline(-1,c="#64748B",lw=.6,ls="--")
    ax.axvline(1,c="#64748B",lw=.6,ls="--")
    ax.axhline(-np.log10(.05),c="#64748B",lw=.6,ls="--")
    if annotate and sig.sum():
        good=np.where(sig&np.isfinite(y))[0]
        chosen=good[np.argsort(tab.q.to_numpy()[good])[:annotate]]
        for i in chosen:
            xx=float(np.clip(x[i],-xlim+.15,xlim-.15));yy=float(np.clip(y[i],0,ylim-.2))
            ax.annotate(str(tab.gene.iloc[i]),(xx,yy),xytext=(3,3),textcoords="offset points",
                        fontsize=6.3,color="#263442",clip_on=True)
    ax.set_xlim(-xlim,xlim);ax.set_ylim(0,ylim)
    ax.text(.98,.97,f"DEG {sig.sum()}  |  n={int(tab.n_case.iloc[0])}/{int(tab.n_ctrl.iloc[0])}",
            transform=ax.transAxes,ha="right",va="top",fontsize=7,
            bbox={"facecolor":"white","alpha":.8,"edgecolor":"none","pad":1})
    if title:ax.set_title(title,fontsize=10,fontweight="bold",pad=7)

for coh in ("sepsis","covid"):
    data={}
    for ct in TYPES:
        ctfile=ct.replace(" ","")
        for source in SOURCES:
            path=(OLD if source=="real" else MAJOR/"deg_tables")/f"{coh}_{ctfile}_{source}.csv.gz"
            data[(ct,source)]=pd.read_csv(path)
    vals=np.concatenate([t.logFC.abs().to_numpy() for t in data.values()])
    xlim=max(5,int(np.ceil(np.nanquantile(vals,.999))))
    vals=np.concatenate([-np.log10(np.maximum(t.q.to_numpy(float),1e-300)) for t in data.values()])
    ylim=max(15,int(np.ceil(np.nanquantile(vals,.995)/5)*5))
    ylim=min(ylim,50)
    fig,axes=plt.subplots(4,4,figsize=(14,12),sharex=True,sharey=True)
    for i,ct in enumerate(TYPES):
        for j,source in enumerate(SOURCES):
            ax=axes[i,j]
            volcano(ax,data[(ct,source)],NAMES[source] if i==0 else None,
                    xlim=xlim,ylim=ylim)
            if j==0:ax.set_ylabel(f"{ct}\n−log10(BH q)",fontsize=9)
            if i==3:ax.set_xlabel("Disease log2FC (case − control)")
    disease="Adult sepsis vs healthy control" if coh=="sepsis" else "COVID positive vs negative ESKD"
    fig.suptitle(f"{disease} — disease DEG volcano plots",fontsize=16,fontweight="bold",y=.995)
    fig.text(.5,.008,"Patient pseudobulk; q < 0.05 and |log2FC| ≥ 1. n = case/control patients. Axes clip extreme values for display.",
             ha="center",fontsize=8,color="#475569")
    fig.tight_layout(rect=(0,.025,1,.975),h_pad=1.3,w_pad=1)
    fig.savefig(MAJOR/f"volcano_{coh}_four_sources.png",dpi=190,bbox_inches="tight")
    fig.savefig(MAJOR/f"volcano_{coh}_four_sources.pdf",bbox_inches="tight")
    plt.close(fig)

for case in ("sepsis","ICU"):
    data={src:pd.read_csv(RARE/f"mono_c01_{case}_vs_HC_{src}.csv.gz") for src in SOURCES}
    vals=np.concatenate([t.logFC.abs().to_numpy() for t in data.values()])
    xlim=max(5,int(np.ceil(np.nanquantile(vals,.999))))
    vals=np.concatenate([-np.log10(np.maximum(t.q.to_numpy(float),1e-300)) for t in data.values()])
    ylim=min(40,max(10,int(np.ceil(np.nanquantile(vals,.995)/5)*5)))
    fig,axes=plt.subplots(1,4,figsize=(15,3.7),sharex=True,sharey=True)
    for ax,src in zip(axes,SOURCES):
        volcano(ax,data[src],NAMES[src],annotate=3,xlim=xlim,ylim=ylim)
        ax.set_xlabel(f"{case} − HC log2FC")
    axes[0].set_ylabel("−log10(BH q)")
    fig.suptitle(f"Sepsis Mono_c01-CD14-HLA-DRA: {case} vs HC (exploratory; ≥10 cells/patient)",
                 fontsize=13,fontweight="bold",y=1.03)
    fig.tight_layout()
    fig.savefig(RARE/f"volcano_mono_c01_{case}_vs_HC.png",dpi=200,bbox_inches="tight")
    fig.savefig(RARE/f"volcano_mono_c01_{case}_vs_HC.pdf",bbox_inches="tight")
    plt.close(fig)
