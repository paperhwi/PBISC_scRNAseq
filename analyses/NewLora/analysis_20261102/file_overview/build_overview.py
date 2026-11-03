"""Reference-style file overview for New LoRA, with real and floor context.

Panel scopes are stated in the figure: full-file counts for cell totals, 20k
uniform random cells per file for expression/depth, and precomputed matched
transfer-label analysis subsets for composition.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches, lines
from scipy.ndimage import gaussian_filter1d

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
OLD=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/file_component_fig_20261018/_extract"
BASE=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/results/big"
SOURCES=["real","prev","floor","new"]
COL={"real":"#404B5A","prev":"#E96836","floor":"#2B83BA","new":"#008C78"}
LINE={"real":"-","prev":":","floor":"--","new":"-"}
LAB={"real":"Real","prev":"PBISC prev K4k","floor":"PBISC floor","new":"New LoRA ceiling"}
TYPES=["T cell","NK cell","B cell","Monocyte","DC","Other"]
TYPECOL=dict(zip(TYPES,["#1BAE79","#EEA600","#E87CA4","#008B55","#7652B7","#D95A4C"]))
LINEAGE=[("T",["CD3E","IL7R","CCR7","CD4","CD8A"]),("NK",["NKG7","GNLY"]),
         ("B / plasma",["MS4A1","CD79A","JCHAIN"]),("Monocyte",["CD14","LYZ","S100A8","FCGR3A"]),
         ("DC",["FCER1A","CLEC4C"]),("Plt / Ery",["PPBP","HBB"]),("Cycling",["MKI67"])]
MARKERS=[g for _,v in LINEAGE for g in v]
NAME={"sepsis":"Sepsis — Ye 2026 (GSE279452)","covid":"COVID-19 — Imperial ESKD (Stephenson 2025)"}
INK="#17212B";MUTED="#475569";GRID="#D9E1E7"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,
    "axes.edgecolor":"#8795A3","text.color":INK,"axes.labelcolor":MUTED,"xtick.color":MUTED,"ytick.color":MUTED,
    "figure.facecolor":"white","axes.facecolor":"white","savefig.facecolor":"white","pdf.fonttype":42})

def get_data(c):
    r=np.load(OLD/f"real_{c}.npz",allow_pickle=True)
    p=np.load(OLD/f"{c}_k4000.npz",allow_pickle=True)
    f=np.load(OLD/f"{c}_floor.npz",allow_pickle=True)
    n=np.load(HERE/f"new_lora_{c}_20k.npz",allow_pickle=True)
    assert np.array_equal(r["gene_symbol"],p["gene_symbol"])
    assert np.array_equal(r["gene_symbol"],f["gene_symbol"])
    assert np.array_equal(r["gene_symbol"],n["gene_symbol"])
    real_tr=pd.read_csv(BASE/f"{c}_realU_transfer.csv.gz")
    prev_tr=pd.read_csv(BASE/f"{c}_prev4k_transfer.csv.gz")
    floor_tr=pd.read_csv(BASE/f"{c}_floor_transfer.csv.gz")
    new_tr=pd.read_csv(ROOT/"NewLora/analysis_20261102"/c/"ceiling_cell_labels.csv.gz")
    comp=[]
    for source,d in [("real",real_tr),("prev",prev_tr),("floor",floor_tr),("new",new_tr)]:
        counts=d.six_tr.value_counts(normalize=True).reindex(TYPES,fill_value=0)*100
        for ct,v in counts.items():comp.append({"cohort":c,"source":source,"six_type":ct,"percent":float(v),"n_labelled":len(d)})
    comp=pd.DataFrame(comp)
    return r,p,f,n,comp

def title(ax,t,pad=6):ax.set_title(t,loc="left",fontsize=11.5,fontweight="bold",pad=pad,color=INK)

def p_genes(ax,c,r,f,n):
    values=[int(r["n_genes_file"]),len(n["gene_symbol"]),int(n["measured"].sum()),int(r["measured"].sum()),2000]
    labels=["Real file (all)","PBISC panel","Measured in cohort","Found in real file","HVG (Harmony)"]
    yy=np.arange(5)[::-1]
    ax.barh(yy,values,color=[COL["real"],"#AAB4C0","#AAB4C0",COL["real"],"#C9D1D9"],height=.62)
    for y,v in zip(yy,values):ax.text(v*1.08,y,f"{v:,}",va="center",fontsize=8)
    ax.set_xscale("log");ax.set_xlim(1000,150000);ax.set_xticks([1000,10000,100000],["1,000","10,000","100,000"])
    ax.set_yticks(yy,labels,fontsize=8);ax.set_xlabel("genes (log scale)",fontsize=8)
    ax.grid(axis="x",color=GRID,zorder=0);ax.set_axisbelow(True);title(ax,"Genes")

def p_cells(ax,c,r,p,f,n):
    vals=[int(r["n_scope"]),int(p["draw0__n_cells"])*5,int(f["draw0__n_cells"]),int(n["n_cells"])]
    samples=[int(r["n_samples"]),int(p["draw0__n_samples"]),int(f["draw0__n_samples"]),int(n["n_samples"])]
    yy=np.arange(4)[::-1]
    ax.barh(yy,np.array(vals)/1e6,color=[COL[s] for s in SOURCES],height=.6)
    for y,v,s in zip(yy,vals,samples):ax.text(v/1e6+.02,y,f"{v:,} · {s} samples",va="center",fontsize=7.8)
    ax.set_yticks(yy,[LAB[s] for s in SOURCES],fontsize=8)
    ax.set_xlim(0,4.15 if c=="sepsis" else 5.15);ax.set_xlabel("cells (millions)",fontsize=8)
    ax.grid(axis="x",color=GRID,zorder=0);ax.set_axisbelow(True);title(ax,"Cells · full files")

def arrays(r,p,f,n):
    return {"real":(r["total_panel"],r["marker_counts"]),
            "prev":(np.concatenate([p[f"draw{i}__total"] for i in range(5)]),
                    np.concatenate([p[f"draw{i}__marker_counts"] for i in range(5)])),
            "floor":(f["draw0__total"],f["draw0__marker_counts"]),
            "new":(n["total"],n["marker_counts"])}

def p_depth(ax,c,r,p,f,n):
    d=arrays(r,p,f,n);bins=np.linspace(2,4.8,75);mid=(bins[:-1]+bins[1:])/2
    ha,_=np.histogram(np.log10(np.maximum(r["total_all"],1)),bins,density=True)
    ax.plot(mid,gaussian_filter1d(ha,1.5),ls="--",lw=1.3,color="#AAB4C0")
    for s,(tot,_) in d.items():
        h,_=np.histogram(np.log10(np.maximum(tot,1)),bins,density=True)
        ax.plot(mid,gaussian_filter1d(h,1.5),lw=1.9,color=COL[s],ls=LINE[s])
    ax.set_xlim(2,4.8);ax.set_xticks([2,3,4],["100","1,000","10,000"]);ax.set_yticks([])
    ax.spines["left"].set_visible(False);ax.set_xlabel("UMI per cell (log scale)",fontsize=8)
    med=[np.median(d[s][0]) for s in SOURCES]
    text="median panel UMI\n"+"\n".join(f"{LAB[s]}: {v:,.0f}" for s,v in zip(SOURCES,med))+f"\nReal all genes: {np.median(r['total_all']):,.0f} (dash)"
    ax.text(.01,.99,text,transform=ax.transAxes,ha="left",va="top",fontsize=7.2,color=MUTED,
            bbox={"facecolor":"white","edgecolor":"none","alpha":.86,"pad":1.5})
    title(ax,"Depth · 6,000 panel genes")

def gene_table(c,r,p,f,n):
    sym=r["gene_symbol"].astype(str);nr=len(r["rows"]);nf=len(f["draw0__rows"]);nn=len(n["rows"])
    tab=pd.DataFrame({"gene":sym,"in_real_file":r["measured"],"measured_in_cohort":n["measured"],
        "real_mean":np.where(r["measured"],r["gene_count_sum"]/nr,np.nan),
        "prev_mean":sum(p[f"draw{i}__gene_count_sum"] for i in range(5))/sum(len(p[f"draw{i}__rows"]) for i in range(5)),
        "floor_mean":f["draw0__gene_count_sum"]/nf,"new_mean":n["gene_count_sum"]/nn})
    tab.to_csv(HERE/f"{c}_gene_mean_20k.csv",index=False)
    stats={}
    for source in ["prev","floor","new"]:
        m=tab.in_real_file&(tab.real_mean>0)&(tab[f"{source}_mean"]>0)
        stats[source]={"pearson_log10":float(np.corrcoef(np.log10(tab.real_mean[m]),np.log10(tab.loc[m,f"{source}_mean"]))[0,1]),"n_genes":int(m.sum())}
    return tab,stats

def p_genecount(ax,c,tab,stats):
    a=tab[tab.in_real_file].sort_values("real_mean",ascending=False)
    b=tab[~tab.in_real_file].sort_values("new_mean",ascending=False)
    t=pd.concat([a,b],ignore_index=True);x=np.arange(len(t));edge=len(a);eps=2e-5
    for s in ["prev","floor","new"]:
        v=t[f"{s}_mean"].clip(lower=eps)
        ax.scatter(x,v,s=1.1,color=COL[s],alpha=.14,lw=0,rasterized=True)
        med=10**pd.Series(np.log10(v)).rolling(151,center=True,min_periods=20).median()
        ax.plot(x,med,color=COL[s],lw=1.8,ls=LINE[s],label=LAB[s])
    ax.plot(x[:edge],t.real_mean.iloc[:edge].clip(lower=eps),color=COL["real"],lw=1.6,label="Real")
    ax.axvspan(edge-.5,len(t)-.5,color="#F1F4F6",zorder=-1)
    ax.text(.98,.12,f"shaded tail: {len(t)-edge} genes absent in real",transform=ax.transAxes,
            ha="right",va="bottom",fontsize=7.0,color=MUTED,
            bbox={"facecolor":"white","edgecolor":"none","alpha":.86,"pad":1})
    ax.set_xlim(0,6000);ax.set_yscale("log");ax.set_ylim(eps*.85,30)
    ax.set_yticks([1e-4,1e-3,1e-2,1e-1,1,10],["0.0001","0.001","0.01","0.1","1","10"])
    ax.set_xlabel("6,000 panel genes ordered by real mean count (left = highest)",fontsize=8)
    ax.set_ylabel("mean UMI / cell",fontsize=8);ax.grid(axis="y",color=GRID,zorder=0)
    for j,g in enumerate(["LYZ","NKG7","CD3E","CD14","MS4A1","CD8A","MKI67"]):
        inds=np.flatnonzero(t.gene.to_numpy()==g)
        if not len(inds):continue
        i=int(inds[0]);v=float(t.real_mean.iloc[i]);
        if not np.isfinite(v) or v<=0:continue
        ax.annotate(g,(i,v),xytext=(300+j*400,8 if j%2==0 else 3),fontsize=7,fontstyle="italic",
                    arrowprops={"arrowstyle":"-","color":"#9AA5B1","lw":.5})
    ax.legend(frameon=False,loc="lower left",fontsize=7.5,ncol=4)
    ax.text(.92,.82,f"per-gene log-mean r vs Real:\nprev {stats['prev']['pearson_log10']:.3f} ({stats['prev']['n_genes']:,}) · floor {stats['floor']['pearson_log10']:.3f} ({stats['floor']['n_genes']:,})\nNew {stats['new']['pearson_log10']:.3f} ({stats['new']['n_genes']:,})",
            transform=ax.transAxes,ha="right",va="top",fontsize=7.6,color=MUTED)
    title(ax,"Gene counts · per-gene mean and 151-gene running median")

def p_comp(ax,c,comp):
    p=comp.pivot(index="source",columns="six_type",values="percent").loc[SOURCES,TYPES]
    yy=np.arange(4)[::-1];left=np.zeros(4)
    for typ in TYPES:
        v=p[typ].to_numpy();ax.barh(yy,v,left=left,color=TYPECOL[typ],height=.6,edgecolor="white",lw=1)
        for y,l,val in zip(yy,left,v):
            if val>=6:ax.text(l+val/2,y,f"{val:.0f}",ha="center",va="center",fontsize=7.5,color="white" if typ in ["Monocyte","DC","Other"] else INK)
        left+=v
    ax.set_yticks(yy,[LAB[s] for s in SOURCES],fontsize=8);ax.set_xlim(0,100)
    ax.set_xlabel("% cells in labelled analysis subsets",fontsize=8)
    ax.legend([patches.Patch(color=TYPECOL[t]) for t in TYPES],TYPES,ncol=6,fontsize=6.6,loc="lower left",
        bbox_to_anchor=(-.02,1.0),frameon=False,handlelength=.65,handletextpad=.2,columnspacing=.55)
    title(ax,"Six-type composition · 30-NN labels",pad=18)

def marker_table(c,r,p,f,n):
    d=arrays(r,p,f,n);records=[]
    for s,(tot,M) in d.items():
        for j,g in enumerate(MARKERS):
            x=M[:,j];det=x>0;v=np.log1p(x[det]/np.maximum(tot[det],1)*1e4)
            records.append({"cohort":c,"source":s,"gene":g,"detected_pct":100*det.mean(),"n_detected":int(det.sum()),
                "mean_log1p_cp10k_detected":float(v.mean()) if len(v) else np.nan})
    tab=pd.DataFrame(records);tab.to_csv(HERE/f"{c}_marker_profile_20k.csv",index=False);return tab

def p_marker(ax,c,r,p,f,n,profile):
    d=arrays(r,p,f,n);bounds=np.linspace(0,9,100);mid=(bounds[:-1]+bounds[1:])/2
    for j,g in enumerate(MARKERS):
        for s,side in [("real",-1),("prev",1),("floor",1),("new",1)]:
            tot,M=d[s];x=M[:,j];det=x>0
            if det.sum()<5:continue
            v=np.log1p(x[det]/np.maximum(tot[det],1)*1e4)
            hist,_=np.histogram(v,bounds,density=True);smooth=gaussian_filter1d(hist,1.5)
            if smooth.max()<=0:continue
            w=.36*smooth/smooth.max()
            keep=smooth>.025*smooth.max()
            if s=="real":ax.fill_betweenx(mid[keep],j,j-w[keep],color=COL[s],alpha=.25,lw=0)
            ax.plot(j+side*w[keep],mid[keep],color=COL[s],lw=1.1,ls=LINE[s])
            ax.plot([j,j+side*.36],[np.median(v)]*2,color=COL[s],lw=.8)
        pct=[profile[(profile.source==s)&(profile.gene==g)].detected_pct.iloc[0] for s in SOURCES]
        ax.text(j,-.45,"/".join(f"{p:.0f}" for p in pct),ha="center",va="top",fontsize=5.2,color=MUTED)
    x0=0
    for lab,genes in LINEAGE:
        ax.text(x0+(len(genes)-1)/2,9.12,lab,ha="center",va="bottom",fontsize=7.5,fontweight="bold")
        if x0>0:ax.axvline(x0-.5,color=GRID,lw=.8)
        x0+=len(genes)
    ax.set_xlim(-.7,len(MARKERS)-.3);ax.set_ylim(-1.1,9.5)
    ax.set_xticks(range(len(MARKERS)),MARKERS,rotation=35,ha="right",fontsize=7.5)
    ax.tick_params(axis="x",pad=12,length=0)
    ax.set_ylabel("log1p(CP10k)\ndetected cells",fontsize=8)
    ax.grid(axis="y",color=GRID,zorder=0)
    title(ax,"Marker expression · detected cells (left Real; right prev / floor / New; numbers = detection %)",pad=20)

def p_annot(ax,c):
    ax.axis("off");ax.set_xlim(0,1);ax.set_ylim(0,1)
    real_scope=("Ye 2026 adult scheme · 108 patients" if c=="sepsis" else
                "Published ESKD roster · 187 samples")
    blocks=[(.98,"Real data",COL["real"],[
        real_scope,
        "20,000 random cells for profile; Real B labels"]),
        (.73,"Previous PBISC K=4,000",COL["prev"],[
        "five draws, 4,000 cells per bulk sample",
        "20,000 pooled cells for profile; 30-NN labels"]),
        (.48,"PBISC floor",COL["floor"],[
        "one draw, 4,000 cells per bulk sample",
        "20,000 random cells for profile; 30-NN labels"]),
        (.23,"New LoRA ceiling",COL["new"],[
        "K_cond=1,000; support-ceiling.v1 readout",
        "one draw, 4,000 cells per bulk sample",
        "20,000 random cells; 60,000 30-NN labels"])]
    for y,head,color,lines_ in blocks:
        ax.add_patch(patches.Rectangle((0,y-.22),.011,.23,color=color,lw=0))
        ax.text(.04,y,head,fontsize=8.8,fontweight="bold",va="top",color=INK)
        for i,t in enumerate(lines_):ax.text(.04,y-.085-i*.06,t,fontsize=7.4,va="top",color=MUTED)
    title(ax,"Source and label scope")

def figure(c):
    r,p,f,n,comp=get_data(c)
    gene,stats=gene_table(c,r,p,f,n);profile=marker_table(c,r,p,f,n)
    comp.to_csv(HERE/f"{c}_six_type_transfer_composition.csv",index=False)
    fig=plt.figure(figsize=(16,9),dpi=100)
    gs=fig.add_gridspec(3,3,height_ratios=[.86,1.11,1.25],width_ratios=[1.1,1.1,1.0],
        left=.09,right=.98,top=.865,bottom=.075,hspace=.66,wspace=.38)
    fig.text(.09,.953,NAME[c]+": real vs PBISC prev / floor / New LoRA",fontsize=18,fontweight="bold",color=INK)
    fig.text(.09,.92,"File-level QC: same 6,000-gene panel · 20,000 cells per source for expression/depth · PBISC prev: five draws; floor/New: one draw · transferred labels for composition",
        fontsize=9.2,color=MUTED)
    p_genes(fig.add_subplot(gs[0,0]),c,r,f,n)
    p_cells(fig.add_subplot(gs[0,1]),c,r,p,f,n)
    p_depth(fig.add_subplot(gs[0,2]),c,r,p,f,n)
    p_genecount(fig.add_subplot(gs[1,:2]),c,gene,stats)
    p_comp(fig.add_subplot(gs[1,2]),c,comp)
    p_marker(fig.add_subplot(gs[2,:2]),c,r,p,f,n,profile)
    p_annot(fig.add_subplot(gs[2,2]),c)
    for ext in ("png","pdf"):
        fig.savefig(HERE/f"new_lora_overview_{c}.{ext}",dpi=240)
    plt.close(fig)
    summary={"cohort":c,"real_scope_cells":int(r["n_scope"]),"prev_file_cells":int(p["draw0__n_cells"])*5,"floor_file_cells":int(f["draw0__n_cells"]),"new_file_cells":int(n["n_cells"]),
        "real_scope_samples":int(r["n_samples"]),"prev_file_samples":int(p["draw0__n_samples"]),"floor_file_samples":int(f["draw0__n_samples"]),"new_file_samples":int(n["n_samples"]),
        "panel_genes":len(n["gene_symbol"]),"new_measured_genes":int(n["measured"].sum()),"genes_in_real":int(r["measured"].sum()),
        "panel_UMI_median":{s:float(np.median(v[0])) for s,v in arrays(r,p,f,n).items()},"gene_logmean_vs_real":stats,
        "composition_n_labelled":{s:int(comp.loc[comp.source==s,"n_labelled"].iloc[0]) for s in SOURCES}}
    (HERE/f"new_lora_overview_{c}_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(c,summary,flush=True)

if __name__=="__main__":
    for cohort in ("sepsis","covid"):figure(cohort)
