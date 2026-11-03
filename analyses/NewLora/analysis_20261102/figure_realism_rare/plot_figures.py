"""Figure analogs for old S3 realism and rare-parent-share panels."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
OLD_R = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/three_model_20261022/results"
COLORS = {"realU":"#3E4C59", "realT1s":"#AEB7C2", "scsemi":"#7B3294", "prev4k":"#D7191C", "floor":"#2B83BA", "ceiling":"#009E73"}
NAMES = {"realU":"real B", "realT1s":"real B, model depth", "scsemi":"scSemiProfiler", "prev4k":"PBISC prev K4000", "floor":"PBISC floor", "ceiling":"new LoRA ceiling"}

plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9, "axes.spines.top":False, "axes.spines.right":False,
                     "pdf.fonttype":42, "savefig.facecolor":"white"})


def save(fig,stem):
    fig.savefig(OUT/f"{stem}.png",dpi=220,bbox_inches="tight")
    fig.savefig(OUT/f"{stem}.pdf",bbox_inches="tight")
    plt.close(fig)


c2st = pd.concat([pd.read_csv(OLD_R/f"{c}_c2st.csv") for c in ("sepsis","covid")] +
                 [pd.read_csv(OUT/f"{c}_new_lora_c2st.csv") for c in ("sepsis","covid")],ignore_index=True)
dc = pd.concat([pd.read_csv(OLD_R/f"{c}_density_coverage.csv") for c in ("sepsis","covid")] +
               [pd.read_csv(OUT/f"{c}_new_lora_density_coverage.csv") for c in ("sepsis","covid")],ignore_index=True)
c2st.to_csv(OUT/"c2st_all_sources.csv",index=False)
dc.to_csv(OUT/"density_coverage_all_sources.csv",index=False)

fig,axs=plt.subplots(1,3,figsize=(18.6,4.5),gridspec_kw={"width_ratios":[1.85,1,1]})
left=axs[0]
srcs=("realU","realT1s","scsemi","prev4k","floor","ceiling")
xt,xlab=[],[]
for ci,coh in enumerate(("sepsis","covid")):
    for li,lin in enumerate(("T","NK","B","Myeloid")):
        base=ci*4.45+li
        for j,src in enumerate(srcs):
            v=c2st.query("cohort == @coh and lineage == @lin and source == @src").c2st_auc
            assert len(v)==1,(coh,lin,src)
            left.bar(base+(j-2.5)*0.14,v.iloc[0],width=0.135,color=COLORS[src])
        xt.append(base);xlab.append(f"{lin}\n{coh.capitalize()}")
left.axhline(.5,color="#94A3B8",ls="--",lw=1)
left.set_xticks(xt,xlab,fontsize=8)
left.set_ylim(.3,1.02)
left.set_ylabel("C2ST AUC vs real A")
left.set_title("Classifier two-sample test (0.5 = indistinguishable)",loc="left",fontsize=10.5)

for ax,metric in zip(axs[1:],("coverage","density")):
    for i,coh in enumerate(("sepsis","covid")):
        for j,src in enumerate(("realU","scsemi","prev4k","floor","ceiling")):
            d=dc.query("cohort == @coh and source == @src and lineage == 'all' and space == 'PCA (no correction)'")
            assert len(d)==1,(coh,src)
            v=d[metric].iloc[0]
            x=i+(j-2)*0.16
            ax.bar(x,v,width=.15,color=COLORS[src])
            ax.text(x,v+.01,f"{v:.2f}",ha="center",va="bottom",fontsize=7)
    ax.set_xticks([0,1],["Sepsis","COVID"])
    ax.set_xlim(-.52,1.52)
    ax.set_ylim(0,max(1.0 if metric=="coverage" else .55,ax.get_ylim()[1]))
    ax.set_ylabel(metric)
    ax.set_title(f"{metric.capitalize()} (k = 5, PCA, no correction)",loc="left",fontsize=10.5)
for ax in axs:
    ax.yaxis.grid(True,color="#E6EBEF",lw=.6)
    ax.set_axisbelow(True)
handles=[plt.Rectangle((0,0),1,1,color=COLORS[s]) for s in srcs]
fig.legend(handles,[NAMES[s] for s in srcs],loc="lower center",ncol=6,frameon=False,fontsize=8.5,bbox_to_anchor=(.5,-.035))
fig.tight_layout(rect=(0,.05,1,1),w_pad=2)
save(fig,"realism_new_lora_comparison")

rare=pd.read_csv(OUT/"rare_parent_share_with_new_lora.csv")
rare=rare[rare.group=="all"]
RARE_SRCS=(("real B (author)","#3E4C59"),("real B (transferred)","#AEB7C2"),
           ("prev K4000 (transferred)","#D7191C"),("floor (transferred)","#2B83BA"),
           ("new LoRA ceiling (transferred)","#009E73"))
fig,axs=plt.subplots(1,2,figsize=(15.5,7.1))
for ax,coh,title in zip(axs,("sepsis","covid"),("Sepsis (Ye 2026, adults)","COVID (Imperial ESKD)")):
    d=rare[rare.cohort==coh]
    subs=d[d.source=="real B (author)"].sort_values("pct_of_parent").subset.tolist()
    par=d.drop_duplicates("subset").set_index("subset").parent
    y=np.arange(len(subs));h=.16
    for j,(source,color) in enumerate(RARE_SRCS):
        vals=d[d.source==source].set_index("subset").reindex(subs).pct_of_parent.to_numpy()
        assert np.isfinite(vals).all()
        ax.barh(y+(2-j)*h,np.maximum(vals,.003),height=h,color=color,label=source)
    ax.set_yticks(y,[f"{s}  [{par[s]}]" for s in subs],fontsize=8)
    ax.set_xscale("log");ax.set_xlim(.003,100)
    ax.set_xlabel("% of parent lineage (pooled; bar at 0.003 = none)")
    ax.set_title(title,loc="left",fontsize=11)
    ax.xaxis.grid(True,color="#E6EBEF",lw=.6);ax.set_axisbelow(True)
handles=[plt.Rectangle((0,0),1,1,color=c) for _,c in RARE_SRCS]
fig.legend(handles,[s for s,_ in RARE_SRCS],loc="lower center",ncol=5,frameon=False,fontsize=8.5,bbox_to_anchor=(.5,-.015))
fig.tight_layout(rect=(0,.045,1,1),w_pad=2)
save(fig,"rare_parent_share_new_lora_comparison")

# Numeric audit comparing the transferred labels for real B with generated labels.
wide=rare.pivot(index=["cohort","subset","parent"],columns="source",values="pct_of_parent").reset_index()
for source in ("prev K4000 (transferred)","floor (transferred)","new LoRA ceiling (transferred)"):
    wide[f"{source} / real_B_transferred"]=(wide[source]+.01)/(wide["real B (transferred)" ]+.01)
wide.to_csv(OUT/"rare_parent_share_wide.csv",index=False)
print("Figures and comparison tables saved",OUT,flush=True)
