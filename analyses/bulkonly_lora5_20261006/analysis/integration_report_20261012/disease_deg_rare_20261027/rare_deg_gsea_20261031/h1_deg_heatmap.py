"""Disease DEG heatmap per cell population, real vs PBISC floor (same pseudobulk DEG tables as the deck: out/allde.pkl).
Rows = union of top-10 real and top-10 floor DEGs (|log2FC| among FDR<0.05 & |log2FC|>=1) per six-type; columns = six-type x source."""
import numpy as np, pandas as pd, pickle, matplotlib, os
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
os.makedirs("heatgsea/out",exist_ok=True)
A=pickle.load(open("out/allde.pkl","rb")); SIX=["T cell","B cell","NK cell","Monocyte","DC"]; N=10
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7})
met=[]
for coh in ["sepsis","covid"]:
    rows=[]; seen=set()
    for ct in SIX:
        for src in ["real","floor"]:
            t=A[(coh,ct)][src]; d=t[t.deg].reindex(t[t.deg].logFC.abs().sort_values(ascending=False).index).head(N)
            for g in d.gene:
                if g not in seen: rows.append((ct,src,g)); seen.add(g)
    genes=[g for _,_,g in rows]; blk=[ct for ct,_,_ in rows]; origin=[s for _,s,_ in rows]
    L=np.zeros((len(genes),10)); Dg=np.zeros_like(L,bool); cols=[]
    for j,(ct,src) in enumerate([(c,s) for c in SIX for s in ["real","floor"]]):
        t=A[(coh,ct)][src].set_index("gene").reindex(genes); L[:,j]=t.logFC.values; Dg[:,j]=t.deg.fillna(False).values; cols.append(f"{ct.replace(' cell','')}\n{src}")
    pd.DataFrame(L,index=genes,columns=[c.replace("\n"," ") for c in cols]).assign(selected_in=blk,selected_by=origin).to_csv(f"heatgsea/out/deg_heatmap_{coh}.csv")
    fig,ax=plt.subplots(figsize=(6.4,0.118*len(genes)+1.6))
    im=ax.imshow(np.clip(np.nan_to_num(L),-4,4),aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(0,-4,4),interpolation="nearest")
    yy,xx=np.where(Dg); ax.scatter(xx,yy,s=3,c="black",marker="o",lw=0)
    ax.set_xticks(range(10)); ax.set_xticklabels(cols); ax.xaxis.tick_top()
    ax.set_yticks(range(len(genes))); ax.set_yticklabels([f"{g}{'  ◆' if o=='floor' else ''}" for g,o in zip(genes,origin)],fontsize=5.6)
    for j in range(1,5): ax.axvline(2*j-.5,color="white",lw=2)
    b=np.array(blk);
    for ct in SIX[1:]:
        i=np.where(b==ct)[0]
        if len(i): ax.axhline(i[0]-.5,color="black",lw=.8)
    for ct in SIX:
        i=np.where(b==ct)[0]
        if len(i): ax.text(9.62,i.mean(),{"T cell":"T","B cell":"B","NK cell":"NK","Monocyte":"Mono","DC":"DC"}[ct],va="center",ha="left",fontsize=8,fontweight="bold")
    cb=fig.colorbar(im,ax=ax,fraction=0.025,pad=0.015,orientation="horizontal",shrink=.6); cb.set_label("disease log2FC (case vs control), clipped ±4")
    ax.set_title(f"{'Sepsis: adult sepsis vs HC' if coh=='sepsis' else 'COVID: positive vs negative'}\n• = DEG (FDR<0.05 & |log2FC|≥1)   ◆ = gene picked from floor top-10",fontsize=8,loc="left",pad=26)
    fig.savefig(f"heatgsea/out/deg_heatmap_{coh}.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig)
    # agreement metrics
    for ct in SIX:
        R=A[(coh,ct)]["real"].set_index("gene"); F=A[(coh,ct)]["floor"].set_index("gene").reindex(R.index)
        u=R.deg|F.deg.fillna(False); rd=set(R.index[R.deg]); fd=set(F.index[F.deg.fillna(False)])
        tr=list(R[R.deg].logFC.abs().sort_values(ascending=False).index[:50]); tf=list(F[F.deg.fillna(False)].logFC.abs().sort_values(ascending=False).index[:50])
        met.append(dict(cohort=coh,celltype=ct,n_real=len(rd),n_floor=len(fd),overlap=len(rd&fd),jaccard=len(rd&fd)/max(len(rd|fd),1),
                        r_union=np.corrcoef(R.logFC[u],F.logFC[u])[0,1] if u.sum()>3 else np.nan,
                        sign_agree_realDEG=float((np.sign(R.logFC[R.deg])==np.sign(F.logFC[R.deg])).mean()) if len(rd) else np.nan,
                        top50_overlap=len(set(tr)&set(tf)),slope_union=np.polyfit(R.logFC[u],F.logFC[u],1)[0] if u.sum()>3 else np.nan))
    # cross-cell-type similarity of disease response (is floor more 'broadcast'?)
    for src in ["real","floor"]:
        M=pd.DataFrame({ct:A[(coh,ct)][src].set_index("gene").logFC for ct in SIX}).dropna()
        C=M.corr().values; iu=np.triu_indices(5,1)
        met.append(dict(cohort=coh,celltype="cross-type mean r",source=src,cross_type_r=float(C[iu].mean())))
Mt=pd.DataFrame(met); Mt.to_csv("heatgsea/out/deg_heatmap_metrics.csv",index=False)
pd.set_option("display.width",220); print(Mt.round(3).to_string(index=False))
