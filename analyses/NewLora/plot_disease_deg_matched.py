"""Plot disease DEG recovery from three matched generated samples."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent/"analysis_20261102/disease_deg_matched"
D=pd.read_csv(HERE/"matched_deg_summary.csv")
COHORTS=("sepsis","covid")
TYPES=("T cell","B cell","NK cell","Monocyte")
METRICS=(("logFC_r","log2FC correlation"),("recall","Real DEG recall"),("precision","Calls overlapping real DEG"))
COLORS={"prev4k":"#C79B5B","floor":"#7A8792","ceiling":"#277DA8"}
NAMES={"prev4k":"Original PBISC","floor":"Floor","ceiling":"New LoRA"}

fig,axes=plt.subplots(2,3,figsize=(13.8,7.3),sharey="col")
for i,coh in enumerate(COHORTS):
    for j,(metric,label) in enumerate(METRICS):
        ax=axes[i,j]
        for source,shift in (("prev4k",-.25),("floor",0),("ceiling",.25)):
            sub=D[(D.cohort==coh)&(D.source==source)].set_index("celltype").reindex(TYPES)
            ax.bar(np.arange(len(TYPES))+shift,sub[metric],width=.23,
                   label=NAMES[source],color=COLORS[source])
        ax.set_xticks(range(len(TYPES)),TYPES,rotation=20,ha="right")
        ax.set_ylim(0,1.02 if metric!="precision" else .6)
        ax.set_title(label,fontsize=12)
        ax.grid(axis="y",alpha=.22)
        ax.set_axisbelow(True)
        if j==0:ax.set_ylabel("Sepsis" if coh=="sepsis" else "COVID / ESKD",fontweight="bold")
fig.suptitle("Disease DEG recovery against full real patient pseudobulk",fontsize=16,fontweight="bold",y=.98)
handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc="lower center",ncol=3,bbox_to_anchor=(.5,.035),frameon=False)
fig.text(.5,.008,"Generated sources: matched ~60k cells each; Welch t test, BH q < 0.05 and |log2FC| >= 1. DC omitted: insufficient generated patient units.",
         ha="center",fontsize=9,color="#475569")
fig.tight_layout(rect=(0,.075,1,.94),h_pad=2)
fig.savefig(HERE/"disease_deg_comparison.png",dpi=200,bbox_inches="tight")
fig.savefig(HERE/"disease_deg_comparison.pdf",bbox_inches="tight")
