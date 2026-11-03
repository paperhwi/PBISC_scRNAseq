"""Heatmap of patient counts for prespecified rare-subtype pseudobulk."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"NewLora/analysis_20261102/rare_disease_deg"
D=pd.read_csv(OUT/"feasibility_by_condition.csv")
SOURCES=("real","prev4k","floor","ceiling")
NAMES={"real":"Real","prev4k":"Original","floor":"Floor","ceiling":"New LoRA"}
fig,axes=plt.subplots(1,2,figsize=(15,7.8),gridspec_kw={"width_ratios":[1,1]})
for ax,coh in zip(axes,("sepsis","covid")):
    sub=D[D.cohort.eq(coh)]
    labels=sorted(sub.subset.unique())
    conds=("case","control")
    cols=[(src,cond) for src in SOURCES for cond in conds]
    p=sub.pivot(index="subset",columns=["source","condition"],values="n_patients_ge10").reindex(index=labels,columns=pd.MultiIndex.from_tuples(cols)).fillna(0).astype(int)
    a=p.to_numpy()
    im=ax.imshow(a,cmap="Blues",vmin=0,vmax=40,aspect="auto")
    ax.set_xticks(range(8),[f"{NAMES[src]}\n{'case' if cond=='case' else ('HC' if coh=='sepsis' else 'neg')}" for src,cond in cols],rotation=45,ha="right")
    ax.set_yticks(range(len(labels)),labels)
    for i in range(len(labels)):
        for j in range(8):
            ax.text(j,i,str(a[i,j]),ha="center",va="center",fontsize=8,color="white" if a[i,j]>=20 else "#1E293B")
    for j in (1.5,3.5,5.5):ax.axvline(j,color="white",lw=2)
    ax.set_title("Sepsis: HC vs sepsis" if coh=="sepsis" else "COVID positive vs negative ESKD",fontsize=12,fontweight="bold")
fig.suptitle("Rare subtype pseudobulk feasibility: patients with ≥10 subtype cells",fontsize=15,fontweight="bold",y=1.01)
fig.text(.5,-.02,"Prespecified rare subtypes; real author labels, generated 30-NN transferred labels. Analysis requires ≥5 patients per arm for exploratory DEG.",
         ha="center",fontsize=9,color="#475569")
fig.tight_layout()
fig.savefig(OUT/"rare_feasibility_heatmap.png",dpi=190,bbox_inches="tight")
fig.savefig(OUT/"rare_feasibility_heatmap.pdf",bbox_inches="tight")
