"""Recreate the floor-era general and myeloid UMAP figure layouts for new LoRA."""
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import umap

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"NewLora/analysis_20261102"
OUT=INPUT/"figure_recreation"
OUT.mkdir(parents=True,exist_ok=True)

plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,
                     "axes.spines.right":False,"font.size":9,"axes.labelsize":9,
                     "axes.titlesize":10,"legend.fontsize":8,"pdf.fonttype":42})
BLUE,RED,GREY,GREEN,ORANGE,PURPLE="#2B83BA","#D7191C","#7F7F7F","#1A9641","#E08214","#7B3294"
SIX_COLORS=dict(zip(["T cell","B cell","NK cell","Monocyte","DC","Other"],
                    ["#4C9BD4","#F4A582","#92C5DE","#D6604D","#B2182B","#BABABA"]))
RARE={"sepsis":["B_c07-MZB1","Macrophage","pDC","B_c08-MKI67"],
      "covid":["pDC","DC1","B_ASC_dividing","MAIT"]}
MYELOID_ORDER={"sepsis":["Mono_c01-CD14-HLA-DRA","Mono_c02-CD14-IL1B","Mono_c03-CD14-RETN",
                          "Mono_c04-CD14-TLN1","Mono_c05-CD14-IFITM1","Mono_c06-CD14-MKI67",
                          "Mono_c07-CD16","Mono_c08-CD14CD16","Macrophage"],
               "covid":["CD14mono","CD14mono_IFN","CD14mono_dex","Int.mono",
                        "CD16mono","CD16mono_IFN","CD16mono_C1"]}
MYELOID_KEY={"sepsis":["Macrophage","Mono_c06-CD14-MKI67","Mono_c08-CD14CD16","Mono_c07-CD16"],
             "covid":["CD14mono_dex","CD16mono_C1","Int.mono","CD14mono_IFN"]}
MYELOID_COLORS=["#1F78B4","#A6CEE3","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#33A02C","#B2DF8A","#6A3D9A"]

def arrows(a):
    a.set_xticks([]);a.set_yticks([])
    for s in a.spines.values():s.set_visible(False)
    a.annotate("",xy=(.115,.02),xytext=(.02,.02),xycoords="axes fraction",
               arrowprops={"arrowstyle":"-|>","color":GREY,"lw":.9})
    a.annotate("",xy=(.02,.115),xytext=(.02,.02),xycoords="axes fraction",
               arrowprops={"arrowstyle":"-|>","color":GREY,"lw":.9})
    a.text(.125,.018,"UMAP 1",transform=a.transAxes,fontsize=7,color=GREY,va="center")
    a.text(.015,.125,"UMAP 2",transform=a.transAxes,fontsize=7,color=GREY,
           rotation=90,ha="center",va="bottom")

def save(fig,stem):
    fig.savefig(OUT/f"{stem}.png",dpi=220,bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/f"{stem}.pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)
    print("saved",stem,flush=True)

for coh,label in (("sepsis","Sepsis (Ye 2026)"),("covid","COVID / ESKD (Imperial)")):
    a=ad.read_h5ad(INPUT/coh/"joint_harmony_embedding.h5ad")
    rng=np.random.default_rng(7)
    idx=np.sort(rng.choice(a.n_obs,min(40000,a.n_obs),replace=False))
    print(coh,"general UMAP start",len(idx),flush=True)
    U=umap.UMAP(n_neighbors=30,min_dist=.3,random_state=0,n_jobs=1).fit_transform(a.obsm["X_pca_harmony"][idx])
    o=a.obs.iloc[idx].copy().reset_index(drop=True)
    o["u1"],o["u2"]=U[:,0],U[:,1]
    o[["source","sample","group","fine","six_type","u1","u2"]].to_csv(
        OUT/f"general_umap_points_{coh}.csv.gz",index=False,compression="gzip")
    ref=o.source.astype(str).eq("real_A").to_numpy()
    fine=o.fine.astype(str).to_numpy()
    six=o.six_type.astype(str).to_numpy()
    fig,ax=plt.subplots(1,3,figsize=(11.2,4.1))
    ax[0].scatter(o.u1[~ref],o.u2[~ref],s=.6,c=BLUE,alpha=.25,lw=0,rasterized=True)
    ax[0].scatter(o.u1[ref],o.u2[ref],s=.6,c=GREY,alpha=.25,lw=0,rasterized=True)
    ax[0].set_title("source",loc="left")
    ax[0].legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,label="real (half A)",ms=5),
                          Line2D([],[],marker="o",ls="",color=BLUE,label="new LoRA ceiling",ms=5)],
                 frameon=False,loc="upper center",bbox_to_anchor=(.5,-.02),ncol=2,
                 handletextpad=.3,columnspacing=1.2)
    for k,v in SIX_COLORS.items():
        mask=six==k
        ax[1].scatter(o.u1[mask],o.u2[mask],s=.6,c=v,alpha=.35,lw=0,rasterized=True)
    ax[1].set_title("cell type (author / transferred)",loc="left")
    ax[1].legend(handles=[Line2D([],[],marker="o",ls="",color=v,label=k,ms=5) for k,v in SIX_COLORS.items()],
                 frameon=False,loc="upper center",bbox_to_anchor=(.5,-.02),ncol=3,
                 handletextpad=.3,columnspacing=.9)
    ax[2].scatter(o.u1,o.u2,s=.5,c="#E8EAED",lw=0,rasterized=True)
    pal=[RED,GREEN,ORANGE,PURPLE]
    for i,r in enumerate(RARE[coh]):
        mask=(fine==r)&ref
        ax[2].scatter(o.u1[mask],o.u2[mask],s=13,marker="o",c=pal[i],alpha=.8,
                      lw=0,zorder=2,rasterized=True)
    for i,r in enumerate(RARE[coh]):
        mask=(fine==r)&~ref
        ax[2].scatter(o.u1[mask],o.u2[mask],s=15,marker="^",c=pal[i],
                      edgecolors="black",linewidths=.45,alpha=.9,zorder=3,rasterized=True)
    ax[2].set_title("rare subsets (● real, ▲ new LoRA with black outline)",loc="left")
    h=[Line2D([],[],marker="s",ls="",color=pal[i],label=r,ms=6) for i,r in enumerate(RARE[coh])]
    h += [Line2D([],[],marker="o",ls="",mfc="#9AA5B1",mec="none",ms=5.5,label="real (half A)"),
          Line2D([],[],marker="^",ls="",mfc="#9AA5B1",mec="black",mew=.8,ms=7,label="new LoRA")]
    ax[2].legend(handles=h,frameon=False,loc="upper center",bbox_to_anchor=(.5,-.02),
                 ncol=3,handletextpad=.3,columnspacing=.9)
    for z in ax:arrows(z)
    fig.suptitle(label+" — new LoRA ceiling",x=.005,ha="left",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,0,1,.97])
    save(fig,f"general_umap_{coh}_ceiling")

    mask=a.obs.six_type.astype(str).eq("Monocyte").to_numpy()
    E=a.obsm["X_pca_harmony"][mask]
    print(coh,"myeloid UMAP start",len(E),flush=True)
    V=umap.UMAP(n_neighbors=30,min_dist=.3,random_state=0,n_jobs=1).fit_transform(E)
    m=a.obs.loc[mask].copy().reset_index(drop=True)
    m["u1"],m["u2"]=V[:,0],V[:,1]
    m[["source","sample","group","fine","six_type","u1","u2"]].to_csv(
        OUT/f"myeloid_umap_points_{coh}.csv.gz",index=False,compression="gzip")
    m=m.sample(frac=1,random_state=1)
    ref=m.source.astype(str).eq("real_A").to_numpy()
    sub=m.fine.astype(str).to_numpy()
    fig,ax=plt.subplots(1,3,figsize=(11.6,4.6))
    ax[0].scatter(m.u1[~ref],m.u2[~ref],s=1.2,c=BLUE,alpha=.35,lw=0,rasterized=True)
    ax[0].scatter(m.u1[ref],m.u2[ref],s=1.2,c=GREY,alpha=.35,lw=0,rasterized=True)
    ax[0].set_title("source",loc="left")
    ax[0].legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,
                                  label=f"real (half A, n={ref.sum():,})",ms=5),
                          Line2D([],[],marker="o",ls="",color=BLUE,
                                 label=f"new LoRA (n={(~ref).sum():,})",ms=5)],
                 frameon=False,loc="upper center",bbox_to_anchor=(.5,-.02),ncol=1,handletextpad=.3)
    cols=dict(zip(MYELOID_ORDER[coh],MYELOID_COLORS))
    for k,v in cols.items():
        q=sub==k
        ax[1].scatter(m.u1[q],m.u2[q],s=1.2,c=v,alpha=.5,lw=0,rasterized=True)
    ax[1].set_title("subset (author / transferred)",loc="left")
    ax[1].legend(handles=[Line2D([],[],marker="o",ls="",color=v,label=k,ms=5) for k,v in cols.items()],
                 frameon=False,loc="upper center",bbox_to_anchor=(.5,-.02),ncol=2,
                 handletextpad=.3,columnspacing=.9)
    ax[2].scatter(m.u1,m.u2,s=1,c="#E8EAED",lw=0,rasterized=True)
    pal=[RED,GREEN,ORANGE,PURPLE]
    K=MYELOID_KEY[coh]
    for rk,i in enumerate(sorted(range(len(K)),key=lambda i:-(sub==K[i]).sum())):
        r=K[i];q=(sub==r)&ref
        ax[2].scatter(m.u1[q],m.u2[q],s=11,marker="o",c=pal[i],alpha=.75,
                      lw=0,zorder=2+2*rk,rasterized=True)
        q=(sub==r)&~ref
        ax[2].scatter(m.u1[q],m.u2[q],s=14 if q.sum()>50 else 26,marker="^",c=pal[i],
                      edgecolors="black",linewidths=.4 if q.sum()>50 else .7,
                      alpha=.9,zorder=3+2*rk,rasterized=True)
    ax[2].set_title("key subsets (● real, ▲ new LoRA with black outline)",loc="left")
    h=[Line2D([],[],marker="s",ls="",color=pal[i],ms=6,
              label=f"{r} ({((sub==r)&ref).sum()} / {((sub==r)&~ref).sum()})") for i,r in enumerate(K)]
    h += [Line2D([],[],marker="o",ls="",mfc="#9AA5B1",mec="none",ms=5.5,label="real (half A)"),
          Line2D([],[],marker="^",ls="",mfc="#9AA5B1",mec="black",mew=.8,ms=7,label="new LoRA")]
    ax[2].legend(handles=h,frameon=False,loc="upper center",bbox_to_anchor=(.5,-.02),
                 ncol=2,handletextpad=.3,columnspacing=.9,title="subset (n real / n new LoRA)",title_fontsize=7.5)
    for z in ax:arrows(z)
    suffix="monocyte / macrophage" if coh=="sepsis" else "monocyte"
    fig.suptitle(label+" — "+suffix+" · new LoRA ceiling",x=.005,ha="left",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,0,1,.97])
    save(fig,f"myeloid_umap_{coh}_ceiling")
