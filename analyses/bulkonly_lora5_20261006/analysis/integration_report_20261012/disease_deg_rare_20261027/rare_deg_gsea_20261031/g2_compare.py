"""Compare GSEA real vs floor: NES correlation, significant-set overlap, sign agreement, discordant sets; Hallmark NES heatmap + scatter."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7})
G=pd.read_csv("heatgsea/out/gsea_all.csv.gz"); SIX=["T cell","B cell","NK cell","Monocyte","DC"]; Q=0.05
W=G.pivot_table(index=["library","cohort","celltype","term"],columns="source",values=["NES","fdr"]); W.columns=[f"{a}_{b}" for a,b in W.columns]; W=W.reset_index()
rows=[]
for (lib,coh,ct),d in W.groupby(["library","cohort","celltype"]):
    d=d.dropna(subset=["NES_real","NES_floor"]); sr=d.fdr_real<Q; sf=d.fdr_floor<Q
    rows.append(dict(library=lib,cohort=coh,celltype=ct,n_sets=len(d),r_NES=np.corrcoef(d.NES_real,d.NES_floor)[0,1],rho_NES=stats.spearmanr(d.NES_real,d.NES_floor).statistic,
                     sig_real=int(sr.sum()),sig_floor=int(sf.sum()),sig_both=int((sr&sf).sum()),jaccard=(sr&sf).sum()/max((sr|sf).sum(),1),
                     sign_agree_real_sig=float((np.sign(d.NES_real[sr])==np.sign(d.NES_floor[sr])).mean()) if sr.sum() else np.nan,
                     real_sig_lost=int((sr&~sf).sum()),floor_only=int((sf&~sr).sum()),opposite_sig=int((sr&sf&(np.sign(d.NES_real)!=np.sign(d.NES_floor))).sum())))
M=pd.DataFrame(rows); M.to_csv("heatgsea/out/gsea_agreement.csv",index=False)
pd.set_option("display.width",230); print(M.round(2).to_string(index=False))
# discordant sets (Hallmark + Reactome): significant in real, lost or reversed in floor, and floor-only
dis=W[(W.fdr_real<Q)&((W.fdr_floor>=0.25)|(np.sign(W.NES_real)!=np.sign(W.NES_floor)))].assign(kind="real sig, lost/reversed in floor")
fon=W[(W.fdr_floor<Q)&(W.fdr_real>=0.25)].assign(kind="floor only")
pd.concat([dis,fon]).sort_values(["library","cohort","celltype","kind"]).to_csv("heatgsea/out/gsea_discordant.csv",index=False)
# Hallmark heatmap per cohort
H=W[W.library=="hallmark"]
for coh in ["sepsis","covid"]:
    h=H[H.cohort==coh]; keep=h[(h.fdr_real<Q)|(h.fdr_floor<Q)].term.unique()
    N=h[h.term.isin(keep)].pivot_table(index="term",columns="celltype",values=["NES_real","NES_floor"]); F=h[h.term.isin(keep)].pivot_table(index="term",columns="celltype",values=["fdr_real","fdr_floor"])
    cols=[(f"NES_{s}",c) for c in SIX for s in ["real","floor"]]; Mx=N.reindex(columns=pd.MultiIndex.from_tuples(cols)); Fx=F.reindex(columns=pd.MultiIndex.from_tuples([(a.replace("NES","fdr"),c) for a,c in cols]))
    order=Mx[[c for c in cols if c[0]=="NES_real"]].mean(1).sort_values(ascending=False).index; Mx=Mx.loc[order]; Fx=Fx.loc[order]
    fig,ax=plt.subplots(figsize=(6.6,0.17*len(Mx)+1.6))
    im=ax.imshow(Mx.values,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(0,-3.2,3.2),interpolation="nearest")
    yy,xx=np.where(Fx.values<Q); ax.scatter(xx,yy,s=5,c="black",lw=0)
    ax.set_xticks(range(10)); ax.set_xticklabels([f"{c.replace(' cell','')}\n{s[4:]}" for s,c in cols]); ax.xaxis.tick_top()
    ax.set_yticks(range(len(Mx))); ax.set_yticklabels([t.replace("HALLMARK_","").replace("_"," ").lower() for t in Mx.index],fontsize=6.2)
    for j in range(1,5): ax.axvline(2*j-.5,color="white",lw=2)
    cb=fig.colorbar(im,ax=ax,fraction=0.03,pad=0.015,orientation="horizontal",shrink=.6); cb.set_label("GSEA NES (case vs control)")
    ax.set_title(f"{'Sepsis' if coh=='sepsis' else 'COVID'} — Hallmark GSEA, real vs PBISC floor   (• FDR < 0.05)",fontsize=8,loc="left",pad=24)
    fig.savefig(f"heatgsea/out/gsea_hallmark_heatmap_{coh}.png",dpi=220,bbox_inches="tight",facecolor="white"); plt.close(fig)
# NES scatter (Hallmark and Reactome)
fig,axs=plt.subplots(2,2,figsize=(9,8.4)); cols5=dict(zip(SIX,["#4C9BD4","#F4A582","#92C5DE","#D6604D","#B2182B"]))
for i,lib in enumerate(["hallmark","reactome"]):
    for j,coh in enumerate(["sepsis","covid"]):
        a=axs[i,j]; d=W[(W.library==lib)&(W.cohort==coh)]
        for ct in SIX:
            e=d[d.celltype==ct]; a.scatter(e.NES_real,e.NES_floor,s=6 if lib=="reactome" else 14,alpha=.6 if lib=="reactome" else .85,color=cols5[ct],lw=0,label=ct)
        lim=np.nanmax(np.abs(d[["NES_real","NES_floor"]].values))*1.05; a.plot([-lim,lim],[-lim,lim],ls="--",c="grey",lw=.7); a.axhline(0,c="grey",lw=.4); a.axvline(0,c="grey",lw=.4)
        m=M[(M.library==lib)&(M.cohort==coh)]
        a.text(0.02,0.98,"\n".join(f"{r.celltype}: r {r.r_NES:.2f}, sig {r.sig_real}/{r.sig_floor}/{r.sig_both}" for r in m.itertuples()),transform=a.transAxes,va="top",fontsize=6.5)
        a.set_xlabel("NES real"); a.set_ylabel("NES PBISC floor"); a.set_title(f"{lib.capitalize()} · {coh}  (sig = real/floor/both, FDR<0.05)",loc="left",fontsize=8)
axs[0,1].legend(frameon=False,loc="lower right",fontsize=7)
fig.tight_layout(); fig.savefig("heatgsea/out/gsea_nes_scatter.png",dpi=220,bbox_inches="tight",facecolor="white"); print("figs ok")
