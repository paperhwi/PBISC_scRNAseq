"""Show how a 496-cell draw quantizes real rare-subtype patient proportions."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"NewLora/analysis_20261102/rare_populations"
SUBSET="T_CD4_c07-MKI67"
SEED=11
TARGET=496

d=pd.read_csv(OUT/"rare_unit_shares_real_vs_ceiling.csv")
sub=d[(d.cohort=="sepsis")&(d.subset==SUBSET)].copy()
real=sub[sub.source=="real"].reset_index(drop=True)
new=sub[sub.source=="ceiling"].sort_values(["group","patient"]).reset_index(drop=True)
assert len(real)==40 and real.n_cells.min()>=TARGET
assert len(new)==121 and new.n_cells.isin([495,496]).all()

rng=np.random.default_rng(SEED)
real["n_thinned"]=TARGET
real["k_thinned"]=[rng.hypergeometric(int(k),int(n-k),TARGET)
                    for n,k in zip(real.n_cells,real.n_subset)]
real["pct_thinned"]=100*real.k_thinned/TARGET
real[["patient","group","n_cells","n_subset","pct_all","n_thinned","k_thinned","pct_thinned"]].to_csv(
    OUT/"real_T_CD4_c07_MKI67_thinned_496_seed11.csv",index=False)

groups=("control","ICU control","case")
names=("HC","ICU","sepsis")
sources=[
    (real,"pct_all","Real: original sample", "#7C8998"),
    (real,"pct_thinned","Same real patients: 496 cells", "#1D6D9C"),
    (new,"pct_all","New LoRA: 495–496 cells", "#D45D48"),
]
fig,axes=plt.subplots(1,3,figsize=(13.6,4.4),sharey=True)
for i,(tab,field,title,color) in enumerate(sources):
    ax=axes[i]
    jitter=np.random.default_rng(99)
    for j,g in enumerate(groups):
        x=tab[tab.group==g]
        xx=j+jitter.uniform(-.24,.24,len(x))
        ax.scatter(xx,x[field],s=20 if i<2 else 14,color=color,alpha=.78,lw=0,zorder=3)
    if i>0:
        for k in range(1,6):
            ax.axhline(100*k/TARGET,color="#CBD5E1",lw=.7,ls="--",zorder=1)
    pos=tab[field].to_numpy()>0
    n_unique=len(np.unique(np.round(tab.loc[pos,field],10)))
    ax.text(.5,.95,f"nonzero {pos.sum()}/{len(tab)}  ·  {n_unique} distinct nonzero %",
            transform=ax.transAxes,ha="center",va="top",fontsize=9,
            bbox={"facecolor":"white","edgecolor":"none","alpha":.9,"pad":2})
    ax.set_xticks(range(3),names)
    ax.set_xlim(-.48,2.48)
    ax.set_ylim(-.035,1.15)
    ax.set_title(title,fontsize=11,fontweight="bold",pad=14)
    ax.grid(axis="y",alpha=.14)
    ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("T_CD4_c07-MKI67 (% of all cells)")
fig.suptitle("The same real patients form 0.2%-spaced bands after 496-cell sampling",
             fontsize=15,fontweight="bold",y=1.02)
fig.text(.5,-.035,"Each dot is one patient; horizontal jitter is for visibility only. Real thinning: hypergeometric draw without replacement, seed 11. Dashed lines: 1–5 cells per 496.",
         ha="center",fontsize=9,color="#475569")
fig.tight_layout()
fig.savefig(OUT/"real_496_downsampling_T_CD4_c07_MKI67.png",dpi=220,bbox_inches="tight")
fig.savefig(OUT/"real_496_downsampling_T_CD4_c07_MKI67.pdf",bbox_inches="tight")
print(f"real original: {sum(real.n_subset>0)}/{len(real)} nonzero, {real.loc[real.n_subset>0,'pct_all'].nunique()} distinct nonzero percentages")
print(f"real thinned: {sum(real.k_thinned>0)}/{len(real)} nonzero, {real.loc[real.k_thinned>0,'pct_thinned'].nunique()} distinct nonzero percentages")
print(f"new LoRA: {sum(new.n_subset>0)}/{len(new)} nonzero, {new.loc[new.n_subset>0,'pct_all'].nunique()} distinct nonzero percentages")
