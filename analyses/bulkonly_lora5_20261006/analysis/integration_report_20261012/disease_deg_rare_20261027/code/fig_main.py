import numpy as np, pandas as pd, pickle, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
                     "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,GREEN,ORANGE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#1A9641","#E08214"
F="figs/"
def save(fig,name):
    fig.savefig(F+name+".png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig); print("saved",name)

# ---------- FIG 1: joint UMAP, source + six type ----------
for coh,lab in [("sepsis","Sepsis (Ye 2026)"),("covid","COVID / ESKD (Imperial)")]:
    o=pd.read_pickle(F+f"{coh}_umap.pkl")
    fig,ax=plt.subplots(1,3,figsize=(11.2,4.1))
    m=o.source.values=="ref"
    ax[0].scatter(o.u1[~m],o.u2[~m],s=0.6,c=BLUE,alpha=.25,lw=0,rasterized=True)
    ax[0].scatter(o.u1[m],o.u2[m],s=0.6,c=GREY,alpha=.25,lw=0,rasterized=True)
    ax[0].set_title("source",loc="left")
    ax[0].legend(handles=[Line2D([],[],marker='o',ls='',color=GREY,label='real (half A)',ms=5),
                          Line2D([],[],marker='o',ls='',color=BLUE,label='PBISC floor',ms=5)],
                 frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=2,handletextpad=.3,columnspacing=1.2)
    six=pd.Series(np.where(o.source.values=="ref",o.six_type.astype(str).values,o.six_tr.astype(str).values))
    cols=dict(zip(["T cell","B cell","NK cell","Monocyte","DC","Other"],["#4C9BD4","#F4A582","#92C5DE","#D6604D","#B2182B","#BABABA"]))
    for k,v in cols.items():
        s=six.values==k; ax[1].scatter(o.u1[s],o.u2[s],s=0.6,c=v,alpha=.35,lw=0,rasterized=True)
    ax[1].set_title("cell type (author / transferred)",loc="left")
    ax[1].legend(handles=[Line2D([],[],marker='o',ls='',color=v,label=k,ms=5) for k,v in cols.items()],
                 frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=3,handletextpad=.3,columnspacing=.9)
    # rare subsets
    RS={"sepsis":["B_c07-MZB1","Macrophage","pDC","B_c08-MKI67"],"covid":["pDC","DC1","B_ASC_dividing","MAIT"]}[coh]
    fine=pd.Series(np.where(o.source.values=="ref",o.fine.astype(str).values,o.fine_tr.astype(str).values))
    ax[2].scatter(o.u1,o.u2,s=0.5,c="#E8EAED",lw=0,rasterized=True)
    pal=[RED,GREEN,ORANGE,"#7B3294"]
    # real first (filled circles, no outline), then floor on top (triangles, black outline)
    for i,r in enumerate(RS):
        s=(fine.values==r)&(o.source.values=="ref")
        if s.sum(): ax[2].scatter(o.u1[s],o.u2[s],s=13,marker="o",c=pal[i],alpha=.8,lw=0,zorder=2,rasterized=True)
    for i,r in enumerate(RS):
        s=(fine.values==r)&(o.source.values=="floor")
        if s.sum(): ax[2].scatter(o.u1[s],o.u2[s],s=15,marker="^",c=pal[i],edgecolors="black",linewidths=0.45,alpha=0.9,zorder=3,rasterized=True)
    ax[2].set_title("rare subsets (● real, ▲ floor with black outline)",loc="left")
    h=[Line2D([],[],marker='s',ls='',color=pal[i],label=r,ms=6) for i,r in enumerate(RS)]
    h+=[Line2D([],[],marker='o',ls='',mfc="#9AA5B1",mec="none",ms=5.5,label="real (half A)"),
        Line2D([],[],marker='^',ls='',mfc="#9AA5B1",mec="black",mew=0.8,ms=7,label="PBISC floor")]
    ax[2].legend(handles=h,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=3,handletextpad=.3,columnspacing=.9)
    for a in ax:
        a.set_xticks([]); a.set_yticks([]); [sp.set_visible(False) for sp in a.spines.values()]
        a.annotate("",xy=(0.115,0.02),xytext=(0.02,0.02),xycoords="axes fraction",
                   arrowprops=dict(arrowstyle="-|>",color=GREY,lw=.9))
        a.annotate("",xy=(0.02,0.115),xytext=(0.02,0.02),xycoords="axes fraction",
                   arrowprops=dict(arrowstyle="-|>",color=GREY,lw=.9))
        a.text(0.125,0.018,"UMAP 1",transform=a.transAxes,fontsize=7,color=GREY,va="center")
        a.text(0.015,0.125,"UMAP 2",transform=a.transAxes,fontsize=7,color=GREY,rotation=90,ha="center",va="bottom")
    fig.suptitle(lab,x=0.005,ha="left",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.97]); save(fig,f"umap_{coh}")

# ---------- FIG 2: volcano ----------
allde=pickle.load(open("out/allde.pkl","rb"))
def volc(ax,T,title,hl=(),color=BLUE):
    q=T.q.values.copy(); q[q<1e-300]=1e-300; y=-np.log10(q); x=T.logFC.values
    sig=(T.q<0.05)&(np.abs(T.logFC)>=1)
    ax.scatter(x[~sig],y[~sig],s=2.5,c="#D5DBE1",lw=0,rasterized=True)
    ax.scatter(x[sig],y[sig],s=3.5,c=color,lw=0,rasterized=True)
    ax.axhline(-np.log10(0.05),ls="--",lw=.7,c=GREY); ax.axvline(1,ls="--",lw=.7,c=GREY); ax.axvline(-1,ls="--",lw=.7,c=GREY)
    for g in hl:
        if g in set(T.gene):
            r=T[T.gene==g].iloc[0]; yy=-np.log10(max(r.q,1e-300))
            ax.scatter([r.logFC],[yy],s=22,facecolor="none",edgecolor=RED,lw=1.1)
            ax.annotate(g,(r.logFC,yy),fontsize=7,color=RED,xytext=(3,3),textcoords="offset points")
    ax.set_title(f"{title}  (DEG n={int(sig.sum())})",loc="left")
    ax.set_xlabel("log2 FC (case vs control)"); ax.set_ylabel("-log10 FDR")
    return float(np.nanmax(np.abs(x)))
panels=[("sepsis","Monocyte",["C1QA","C1QB","C1QC","IL1R2","HLA-DRA","CLU"]),
        ("sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","IL1R2"]),
        ("covid","Monocyte",["IFI27","SIGLEC1","IFIT1","MX1"])]
fig,ax=plt.subplots(2,3,figsize=(11.6,6.4))
for j,(coh,ct,hl) in enumerate(panels):
    t=allde[(coh,ct)]
    volc(ax[0,j],t["real"],f"{coh} {ct} — real",hl,GREY)
    volc(ax[1,j],t["floor"],f"{coh} {ct} — PBISC floor",hl,BLUE)
    lo=min(ax[0,j].get_xlim()[0],ax[1,j].get_xlim()[0]); hi=max(ax[0,j].get_xlim()[1],ax[1,j].get_xlim()[1])
    ax[0,j].set_xlim(lo,hi); ax[1,j].set_xlim(lo,hi)
fig.tight_layout(); save(fig,"volcano")

# ---------- FIG 3: logFC scatter ----------
fig,ax=plt.subplots(1,3,figsize=(11.2,3.8))
for j,(coh,ct,hl) in enumerate(panels):
    t=allde[(coh,ct)]; R=t["real"].set_index("gene"); Fl=t["floor"].set_index("gene").loc[R.index]
    r=np.corrcoef(R.logFC,Fl.logFC)[0,1]; sl=np.polyfit(R.logFC,Fl.logFC,1)[0]
    ax[j].scatter(R.logFC,Fl.logFC,s=2.5,c="#C9D2DA",lw=0,rasterized=True)
    d=R.deg.values; ax[j].scatter(R.logFC[d],Fl.logFC[d],s=4,c=BLUE,lw=0,rasterized=True)
    lim=np.nanpercentile(np.abs(np.r_[R.logFC,Fl.logFC]),99.8)
    ax[j].plot([-lim,lim],[-lim,lim],ls="--",lw=.8,c=GREY); ax[j].plot([-lim,lim],[-lim*sl,lim*sl],lw=1.1,c=RED)
    for g in hl:
        if g in R.index:
            ax[j].scatter([R.logFC[g]],[Fl.logFC[g]],s=24,facecolor="none",edgecolor=RED,lw=1.1)
            ax[j].annotate(g,(R.logFC[g],Fl.logFC[g]),fontsize=7,color=RED,xytext=(3,3),textcoords="offset points")
    ax[j].set_xlim(-lim,lim); ax[j].set_ylim(-lim,lim)
    ax[j].set_title(f"{coh} {ct}",loc="left"); ax[j].set_xlabel("real log2FC"); ax[j].set_ylabel("floor log2FC")
    ax[j].text(.03,.95,f"r = {r:.2f}\nslope = {sl:.2f}",transform=ax[j].transAxes,va="top",fontsize=8,color=INK)
fig.tight_layout(); save(fig,"logfc_scatter")
