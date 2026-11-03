import numpy as np, pandas as pd, pickle, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,GREEN,ORANGE,PURPLE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#1A9641","#E08214","#7B3294"
F="figs/"
save=lambda fig,n:(fig.savefig(F+n+".png",dpi=220,bbox_inches="tight",facecolor="white"),plt.close(fig),print("saved",n))
D=pd.read_csv("out/deg_summary.csv"); C=pd.read_csv("out/real_real_ceiling.csv")
CT=["T cell","B cell","NK cell","Monocyte","DC"]

# ---- FIG 4: logFC r by cell type with ceiling ----
fig,ax=plt.subplots(1,2,figsize=(10.6,4.2))
for j,coh in enumerate(["sepsis","covid"]):
    d=D[D.cohort==coh].set_index(["celltype","source"]); c=C[C.cohort==coh].set_index("celltype")
    x=np.arange(len(CT)); w=0.36
    pv=[d.loc[(t,"prev"),"logFC_r"] for t in CT]; fv=[d.loc[(t,"floor"),"logFC_r"] for t in CT]
    cv=[c.loc[t,"logFC_r"] for t in CT]; cs=[c.loc[t,"logFC_r_sd"] for t in CT]
    ax[j].bar(x-w/2,pv,w,color="#F2A7A8",edgecolor=RED,lw=.8,label="PBISC prev K4000")
    ax[j].bar(x+w/2,fv,w,color="#A8CDE3",edgecolor=BLUE,lw=.8,label="PBISC floor")
    for i in x:
        ax[j].hlines(cv[i],i-0.46,i+0.46,color=INK,lw=1.6)
        ax[j].add_patch(plt.Rectangle((i-0.46,cv[i]-cs[i]),0.92,2*cs[i],color=INK,alpha=.13,lw=0))
        ax[j].text(i+w/2,fv[i]+0.02,f"{100*fv[i]/cv[i]:.0f}%",ha="center",fontsize=7.5,color=BLUE)
    ax[j].set_xticks(x); ax[j].set_xticklabels(CT,rotation=20,ha="right")
    ax[j].set_ylim(0,1.05); ax[j].set_ylabel("disease logFC Pearson r vs real")
    ax[j].set_title({"sepsis":"Sepsis — sepsis vs healthy control","covid":"COVID/ESKD — positive vs negative"}[coh],loc="left")
fig.legend(handles=[plt.Rectangle((0,0),1,1,fc="#F2A7A8",ec=RED),plt.Rectangle((0,0),1,1,fc="#A8CDE3",ec=BLUE),
                    Line2D([],[],color=INK,lw=1.6)],labels=["PBISC prev K4000","PBISC floor","real-vs-real ceiling (±SD)"],
           frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.015),ncol=3,columnspacing=2.2,handletextpad=.6)
fig.tight_layout(rect=[0,0.07,1,1]); save(fig,"bar_logfc_r")

# ---- FIG 5: recall / precision vs ceiling ----
fig,ax=plt.subplots(1,2,figsize=(10.6,3.9))
for j,m in enumerate(["recall","precision"]):
    rows=[];
    for coh in ["sepsis","covid"]:
        d=D[(D.cohort==coh)&(D.source=="floor")].set_index("celltype"); c=C[C.cohort==coh].set_index("celltype")
        for t in CT: rows.append((f"{coh[:3]} {t.replace(' cell','')}",d.loc[t,m],c.loc[t,m]))
    lab=[r[0] for r in rows]; fv=[r[1] for r in rows]; cv=[r[2] for r in rows]
    x=np.arange(len(rows))
    ax[j].bar(x,fv,0.6,color="#A8CDE3",edgecolor=BLUE,lw=.8)
    ax[j].scatter(x,cv,marker="_",s=320,color=INK,lw=1.8,zorder=3)
    ax[j].set_xticks(x); ax[j].set_xticklabels(lab,rotation=55,ha="right"); ax[j].set_ylim(0,0.85)
    ax[j].set_ylabel(f"{m} of real DEGs" if m=="recall" else "precision of generated DEGs")
    ax[j].set_title({"recall":"DEG recall","precision":"DEG precision"}[m],loc="left")
fig.legend(handles=[plt.Rectangle((0,0),1,1,fc="#A8CDE3",ec=BLUE),Line2D([],[],color=INK,lw=1.8)],
           labels=["PBISC floor","real-vs-real ceiling"],frameon=False,loc="upper center",
           bbox_to_anchor=(0.5,1.015),ncol=2,columnspacing=2.2,handletextpad=.6)
fig.tight_layout(rect=[0,0,1,0.93]); save(fig,"bar_recall_precision")

# ---- FIG 6: module recovery ----
N=pd.read_csv("out/named_genes.csv"); f=N[N.source=="floor"]
MOD=[("C1Q / complement\n(sepsis Mono)","sepsis","Monocyte",["C1QA","C1QB","C1QC","RNASE1","CLU"]),
     ("MHC-II\n(sepsis Mono)","sepsis","Monocyte",["HLA-DRA","HLA-DRB1","HLA-DPA1","CD74"]),
     ("pDC/cDC identity\n(sepsis DC)","sepsis","DC",["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C"]),
     ("sepsis myeloid\nactivation (Mono)","sepsis","Monocyte",["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1"]),
     ("ISG\n(COVID Mono)","covid","Monocyte",["ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1"]),
     ("IFI27 + SIGLEC1\n(COVID Mono)","covid","Monocyte",["IFI27","SIGLEC1"])]
fig,ax=plt.subplots(figsize=(10.4,3.9))
names=[];vals=[];errs=[];cols=[]
for nm,coh,ct,gs in MOD:
    d=f[(f.cohort==coh)&(f.celltype==ct)&(f.gene.isin(gs))].dropna(subset=["ratio"])
    names.append(nm); vals.append(d.ratio.median()); errs.append(d.ratio.std()/np.sqrt(max(len(d),1)))
    cols.append(RED if d.ratio.median()<0.5 else (ORANGE if d.ratio.median()<0.8 else BLUE))
x=np.arange(len(names))
ax.bar(x,vals,0.58,color=[c+"55" for c in cols],edgecolor=cols,lw=1.1)
ax.errorbar(x,vals,yerr=errs,fmt="none",ecolor=INK,capsize=3,lw=.9)
ax.axhline(1,ls="--",lw=.9,c=INK); ax.axhline(0.5,ls=":",lw=.9,c=RED)
for i,v in enumerate(vals): ax.text(i,v+0.05,f"{v:.2f}",ha="center",fontsize=8.5,color=INK)
ax.set_xticks(x); ax.set_xticklabels(names,fontsize=8); ax.set_ylim(0,1.5)
ax.set_ylabel("floor / real  disease log2FC (median ± SEM)")
ax.set_title("Disease-programme recovery by gene module",loc="left")
ax.text(0.995,0.93,"dashed = perfect   dotted = 0.5×",transform=ax.transAxes,ha="right",fontsize=8,color=GREY)
fig.tight_layout(); save(fig,"bar_modules")

# ---- FIG 7: S100A8/9 leakage ----
fig,ax=plt.subplots(figsize=(6.9,4.0))
d=f[(f.cohort=="covid")&(f.gene.isin(["S100A8","S100A9"]))]
order=["Monocyte","DC","T cell","B cell","NK cell"]
x=np.arange(len(order)); w=0.38
for k,(g,off,col) in enumerate([("S100A8",-w/2,BLUE),("S100A9",w/2,"#5FA8D3")]):
    dd=d[d.gene==g].set_index("celltype")
    ax.bar(x+off,[dd.loc[t,"real_logFC"] for t in order],w,color="#D9DEE3",edgecolor=GREY,lw=.8,label="real (bars)" if k==0 else None)
for k,(g,off,col) in enumerate([("S100A8",-w/2,BLUE),("S100A9",w/2,"#5FA8D3")]):
    dd=d[d.gene==g].set_index("celltype")
    ax.plot(x+off,[dd.loc[t,"gen_logFC"] for t in order],"o",color=col,ms=7,label=f"floor {g}")
    for i,t in enumerate(order):
        if dd.loc[t,"gen_DEG"]: ax.text(x[i]+off,dd.loc[t,"gen_logFC"]+0.09,"*",ha="center",fontsize=13,color=RED)
ax.axhline(0,lw=.8,c=INK); ax.set_xticks(x); ax.set_xticklabels(order,rotation=20,ha="right")
ax.set_ylabel("disease log2FC (COVID+ vs COVID−)"); ax.set_ylim(-0.2,1.9)
fig.legend(frameon=False,fontsize=8,loc="lower center",bbox_to_anchor=(0.5,-0.02),ncol=3,
           columnspacing=1.6,handletextpad=.5)
ax.set_title("Alarmin leakage: S100A8/A9 become DEGs in lymphoid lineages\nbars = real, dots = floor,  * = called DEG in floor",loc="left",fontsize=10)
fig.tight_layout(rect=[0,0.09,1,1]); save(fig,"bar_s100_leak")
