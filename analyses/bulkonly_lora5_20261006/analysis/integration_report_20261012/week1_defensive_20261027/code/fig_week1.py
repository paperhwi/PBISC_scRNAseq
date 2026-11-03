import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,
 "font.size":9,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"xtick.labelsize":8,"ytick.labelsize":8})
BLUE,RED,GREY,INK,GREEN,ORANGE="#2B83BA","#D7191C","#7F7F7F","#1F2933","#1A9641","#E08214"
F="figs/"; save=lambda fig,n:(fig.savefig(F+n+".png",dpi=220,bbox_inches="tight",facecolor="white"),plt.close(fig),print("saved",n))

# --- A2: is the module in the bulk? ---
M=pd.read_csv("out2/a2_module_bulk.csv").dropna(subset=["bulk_logFC"])
g=M.groupby("module")
S=pd.DataFrame({"bulk/real":g.apply(lambda d:(d.bulk_logFC/d.real_logFC).median()),
                "floor/real":g.apply(lambda d:(d.floor_logFC/d.real_logFC).median())})
order=["pDC / cDC identity","C1Q / complement","IFI27 / SIGLEC1","MHC-II","ISG (COVID)","sepsis myeloid activation"]
S=S.loc[[o for o in order if o in S.index]]
fig,ax=plt.subplots(figsize=(8.6,4.0)); x=np.arange(len(S)); w=0.38
ax.bar(x-w/2,S["bulk/real"],w,color="#C9B8D8",edgecolor="#7B3294",lw=1.0,label="bulk input")
ax.bar(x+w/2,S["floor/real"],w,color="#A8CDE3",edgecolor=BLUE,lw=1.0,label="PBISC floor")
ax.axhline(1,ls="--",lw=.9,c=INK); ax.axhline(0.5,ls=":",lw=.9,c=RED)
for i,(b,f) in enumerate(zip(S["bulk/real"],S["floor/real"])):
    ax.text(i-w/2,b+0.03,f"{b:.2f}",ha="center",fontsize=8); ax.text(i+w/2,f+0.03,f"{f:.2f}",ha="center",fontsize=8,color=BLUE)
ax.set_xticks(x); ax.set_xticklabels([s.replace(" / ","/\n") for s in S.index],fontsize=8)
ax.set_ylabel("module disease log2FC / real single-cell"); ax.set_ylim(0,1.6)
ax.set_title("The signal is in the bulk — floor fails to route it to the cell type",loc="left")
fig.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.02),ncol=2,columnspacing=2.4)
fig.tight_layout(rect=[0,0.07,1,1]); save(fig,"w1_bulk_vs_floor")

# --- C2b: positivity vs load ---
L=pd.read_csv("out2/c2b_module_load.csv")
p=L.pivot_table(index="module",columns="source",values=["pct_pos","mean_CP10k_in_pos"])
R=pd.DataFrame({"pos%":p[("pct_pos","floor")]/p[("pct_pos","real (depth-matched)")],
                "load":p[("mean_CP10k_in_pos","floor")]/p[("mean_CP10k_in_pos","real (depth-matched)")],
                "pos%_prev":p[("pct_pos","prev K4000 (Poisson)")]/p[("pct_pos","real (depth-matched)")],
                "load_prev":p[("mean_CP10k_in_pos","prev K4000 (Poisson)")]/p[("mean_CP10k_in_pos","real (depth-matched)")]})
R=R.loc[[o for o in ["C1Q / complement","pDC / cDC identity","plasma-cell secretory","sepsis myeloid activation","ISG (COVID)"] if o in R.index]]
fig,ax=plt.subplots(1,2,figsize=(10.6,4.0))
x=np.arange(len(R)); w=0.38
for j,(ttl,c1,c2) in enumerate([("PBISC floor","pos%","load"),("prev K4000 (Poisson read-out)","pos%_prev","load_prev")]):
    ax[j].bar(x-w/2,R[c1],w,color="#D5DBE1",edgecolor=GREY,lw=1.0)
    ax[j].bar(x+w/2,R[c2],w,color="#A8CDE3" if j==0 else "#F2A7A8",edgecolor=BLUE if j==0 else RED,lw=1.0)
    ax[j].axhline(1,ls="--",lw=.9,c=INK)
    for i,(a,b) in enumerate(zip(R[c1],R[c2])):
        ax[j].text(i-w/2,a+0.03,f"{a:.2f}",ha="center",fontsize=7.5)
        ax[j].text(i+w/2,b+0.03,f"{b:.2f}",ha="center",fontsize=7.5)
    ax[j].set_xticks(x); ax[j].set_xticklabels([s.replace(" / ","/\n").replace(" secretory","\nsecretory").replace(" activation","\nactivation") for s in R.index],fontsize=7.5)
    ax[j].set_ylim(0,2.0); ax[j].set_title(ttl,loc="left"); ax[j].set_ylabel("ratio vs depth-matched real")
fig.legend(handles=[plt.Rectangle((0,0),1,1,fc="#D5DBE1",ec=GREY),plt.Rectangle((0,0),1,1,fc="#A8CDE3",ec=BLUE)],
           labels=["% of cells positive for the module","transcript load inside positive cells"],
           frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.02),ncol=2,columnspacing=2.4)
fig.suptitle("Right number of positive cells, wrong amount inside them",x=0.005,ha="left",fontsize=10.5,fontweight="bold")
fig.tight_layout(rect=[0,0.08,1,0.95]); save(fig,"w1_pos_vs_load")

# --- A3: ceiling bracket ---
A=pd.read_csv("out2/a3_power_matched_ceiling.csv"); C=pd.read_csv("out/real_real_ceiling.csv")
A=A.merge(C[["cohort","celltype","logFC_r"]].rename(columns={"logFC_r":"ceil_bothsplit"}),on=["cohort","celltype"])
fig,ax=plt.subplots(1,2,figsize=(10.6,3.9))
CT=["T cell","B cell","NK cell","Monocyte","DC"]
for j,coh in enumerate(["sepsis","covid"]):
    d=A[A.cohort==coh].set_index("celltype").loc[CT]; x=np.arange(len(CT))
    ax[j].fill_between(x,d.ceil_bothsplit,d.ceil_case_split_r,color=INK,alpha=.12,step="mid")
    ax[j].plot(x,d.ceil_bothsplit,"_",ms=28,color=GREY,mew=2)
    ax[j].plot(x,d.ceil_case_split_r,"_",ms=28,color=INK,mew=2)
    ax[j].bar(x,d.floor_r_full,0.45,color="#A8CDE3",edgecolor=BLUE,lw=1.0)
    for i,t in enumerate(CT):
        lo=100*d.floor_r_full[t]/d.ceil_case_split_r[t]; hi=100*d.floor_r_full[t]/d.ceil_bothsplit[t]
        ax[j].text(i,d.floor_r_full[t]+0.025,f"{lo:.0f}–{min(hi,999):.0f}%",ha="center",fontsize=7.5,color=BLUE)
    ax[j].set_xticks(x); ax[j].set_xticklabels(CT,rotation=20,ha="right"); ax[j].set_ylim(0,1.05)
    ax[j].set_ylabel("disease logFC Pearson r"); ax[j].set_title({"sepsis":"Sepsis","covid":"COVID / ESKD"}[coh],loc="left")
fig.legend(handles=[Line2D([],[],color=GREY,lw=2),Line2D([],[],color=INK,lw=2),plt.Rectangle((0,0),1,1,fc="#A8CDE3",ec=BLUE)],
           labels=["ceiling, both arms split (under-powered → low)","ceiling, controls shared (correlated → high)","PBISC floor"],
           frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.02),ncol=3,columnspacing=1.8)
fig.suptitle("The ceiling is itself uncertain — report floor as a bracket",x=0.005,ha="left",fontsize=10.5,fontweight="bold")
fig.tight_layout(rect=[0,0.09,1,0.95]); save(fig,"w1_ceiling_bracket")

# --- D1 CI forest ---
D=pd.read_csv("out2/d1_bootstrap_ci.csv")
fig,ax=plt.subplots(figsize=(7.4,4.2)); y=np.arange(len(D))[::-1]
ax.hlines(y,D.r_lo,D.r_hi,color=GREY,lw=2.2)
ax.plot(D.logFC_r,y,"o",color=BLUE,ms=7)
ax.axvline(0,lw=.8,c=INK)
ax.set_yticks(y); ax.set_yticklabels([f"{r.cohort[:3]}  {r.celltype}" for r in D.itertuples()],fontsize=8)
ax.set_xlabel("disease logFC r vs real  (point estimate, 95% patient-bootstrap CI)")
ax.set_xlim(-0.05,1.0); ax.set_title("Uncertainty of the headline metric",loc="left")
fig.tight_layout(); save(fig,"w1_bootstrap_ci")
