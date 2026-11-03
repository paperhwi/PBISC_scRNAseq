import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
  "font.size":9,"axes.titlesize":10,"legend.fontsize":7.5})
BLUE,RED,GREY,GREEN,ORANGE,PURPLE="#2B83BA","#D7191C","#7F7F7F","#1A9641","#E08214","#7B3294"
ORDER={"sepsis":["Mono_c01-CD14-HLA-DRA","Mono_c02-CD14-IL1B","Mono_c03-CD14-RETN","Mono_c04-CD14-TLN1",
                 "Mono_c05-CD14-IFITM1","Mono_c06-CD14-MKI67","Mono_c07-CD16","Mono_c08-CD14CD16","Macrophage"],
       "covid":["CD14mono","CD14mono_IFN","CD14mono_dex","Int.mono","CD16mono","CD16mono_IFN","CD16mono_C1"]}
KEY={"sepsis":["Macrophage","Mono_c06-CD14-MKI67","Mono_c08-CD14CD16","Mono_c07-CD16"],
     "covid":["CD14mono_dex","CD16mono_C1","Int.mono","CD14mono_IFN"]}
PAL=["#1F78B4","#A6CEE3","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#33A02C","#B2DF8A","#6A3D9A"]
LAB={"sepsis":"Sepsis (Ye 2026) — monocyte / macrophage","covid":"COVID / ESKD (Imperial) — monocyte"}
def arrows(a):
    a.set_xticks([]); a.set_yticks([]); [s.set_visible(False) for s in a.spines.values()]
    for xy in [(0.115,0.02),(0.02,0.115)]:
        a.annotate("",xy=xy,xytext=(0.02,0.02),xycoords="axes fraction",arrowprops=dict(arrowstyle="-|>",color=GREY,lw=.9))
    a.text(0.125,0.018,"UMAP 1",transform=a.transAxes,fontsize=7,color=GREY,va="center")
    a.text(0.015,0.125,"UMAP 2",transform=a.transAxes,fontsize=7,color=GREY,rotation=90,ha="center",va="bottom")
for coh in ["sepsis","covid"]:
    o=pd.read_pickle(f"myeloid/{coh}_myeloid_umap.pkl")
    o=o.sample(frac=1,random_state=1)  # shuffle draw order
    ref=o.source.values=="ref"; sub=o["sub"].values
    fig,ax=plt.subplots(1,3,figsize=(11.6,4.6))
    ax[0].scatter(o.u1[~ref],o.u2[~ref],s=1.2,c=BLUE,alpha=.35,lw=0,rasterized=True)
    ax[0].scatter(o.u1[ref],o.u2[ref],s=1.2,c=GREY,alpha=.35,lw=0,rasterized=True)
    ax[0].set_title("source",loc="left")
    ax[0].legend(handles=[Line2D([],[],marker='o',ls='',color=GREY,label=f'real (half A, n={ref.sum():,})',ms=5),
                          Line2D([],[],marker='o',ls='',color=BLUE,label=f'PBISC floor (n={(~ref).sum():,})',ms=5)],
                 frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=1,handletextpad=.3)
    cols=dict(zip(ORDER[coh],PAL))
    for k,v in cols.items():
        s=sub==k; ax[1].scatter(o.u1[s],o.u2[s],s=1.2,c=v,alpha=.5,lw=0,rasterized=True)
    ax[1].set_title("subset (author / transferred)",loc="left")
    ax[1].legend(handles=[Line2D([],[],marker='o',ls='',color=v,label=k,ms=5) for k,v in cols.items()],
                 frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=2,handletextpad=.3,columnspacing=.9)
    ax[2].scatter(o.u1,o.u2,s=1,c="#E8EAED",lw=0,rasterized=True)
    pal=[RED,GREEN,ORANGE,PURPLE]; K=KEY[coh]
    # draw abundant subsets first, rarest last; within each subset real circles then floor triangles on top
    for rk,i in enumerate(sorted(range(len(K)),key=lambda i:-(sub==K[i]).sum())):
        r=K[i]; s=(sub==r)&ref
        ax[2].scatter(o.u1[s],o.u2[s],s=11,marker="o",c=pal[i],alpha=.75,lw=0,zorder=2+2*rk,rasterized=True)
        s=(sub==r)&~ref
        ax[2].scatter(o.u1[s],o.u2[s],s=14 if s.sum()>50 else 26,marker="^",c=pal[i],edgecolors="black",
                      linewidths=0.4 if s.sum()>50 else 0.7,alpha=.9,zorder=3+2*rk,rasterized=True)
    ax[2].set_title("key subsets (● real, ▲ floor with black outline)",loc="left")
    h=[Line2D([],[],marker='s',ls='',color=pal[i],ms=6,
              label=f"{r}  ({((sub==r)&ref).sum()} / {((sub==r)&~ref).sum()})") for i,r in enumerate(K)]
    h+=[Line2D([],[],marker='o',ls='',mfc="#9AA5B1",mec="none",ms=5.5,label="real (half A)"),
        Line2D([],[],marker='^',ls='',mfc="#9AA5B1",mec="black",mew=0.8,ms=7,label="PBISC floor")]
    ax[2].legend(handles=h,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=2,handletextpad=.3,columnspacing=.9,
                 title="subset (n real / n floor)",title_fontsize=7.5)
    for a in ax: arrows(a)
    fig.suptitle(LAB[coh],x=0.005,ha="left",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig(f"myeloid/umap_myeloid_{coh}.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig)

# composition: share within myeloid compartment, real vs floor
fig,ax=plt.subplots(1,2,figsize=(11.6,3.9),gridspec_kw={"width_ratios":[9,7]})
rows=[]
for j,coh in enumerate(["sepsis","covid"]):
    o=pd.read_pickle(f"myeloid/{coh}_myeloid_umap.pkl")
    ct=pd.crosstab(o["sub"],o.source).reindex(ORDER[coh]).fillna(0)
    fr=ct/ct.sum()*100; x=np.arange(len(ct)); w=.38
    ax[j].bar(x-w/2,fr["ref"],w,color=GREY,label="real (half A)")
    ax[j].bar(x+w/2,fr["floor"],w,color=BLUE,label="PBISC floor")
    for i,g in enumerate(ct.index):
        ratio=fr.loc[g,"floor"]/fr.loc[g,"ref"]
        ax[j].text(i,max(fr.loc[g])+0.8,f"{ratio:.2f}×",ha="center",fontsize=7,color=RED if ratio<0.5 else "#3E4C59")
        rows.append(dict(cohort=coh,subset=g,n_real=int(ct.loc[g,"ref"]),n_floor=int(ct.loc[g,"floor"]),
                         pct_real=round(fr.loc[g,"ref"],2),pct_floor=round(fr.loc[g,"floor"],2),floor_over_real=round(ratio,3)))
    ax[j].set_xticks(x); ax[j].set_xticklabels([g.replace("Mono_","") for g in ct.index],rotation=40,ha="right",fontsize=7.5)
    ax[j].set_ylabel("% of monocyte compartment"); ax[j].set_title(f"{coh}: subset share (label = floor/real)",loc="left")
    ax[j].set_ylim(0,fr.values.max()*1.15+2)
ax[1].legend(frameon=False,loc="upper right")
fig.tight_layout(); fig.savefig("myeloid/composition_myeloid.png",dpi=220,bbox_inches="tight",facecolor="white")
pd.DataFrame(rows).to_csv("myeloid/myeloid_composition.csv",index=False); print(pd.DataFrame(rows).to_string())
