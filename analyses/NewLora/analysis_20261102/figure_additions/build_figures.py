"""Slide-sized, source-traceable figures for the New LoRA comparison."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)
COL = {"real": "#374151", "prev": "#7B3294", "floor": "#2B83BA", "new": "#008C78", "pale": "#AAB4C0", "grid": "#D9E1E7"}
SOURCES = [("prev4k", "기존 PBISC", COL["prev"]), ("floor", "floor", COL["floor"]),
           ("new_ceiling", "새 LoRA", COL["new"]), ("realU", "real B", COL["real"])]
SC = [("prev4k", "기존 PBISC", COL["prev"]), ("floor", "floor", COL["floor"]),
      ("ceiling", "새 LoRA", COL["new"]), ("realB", "real B", COL["real"])]
plt.rcParams.update({"font.family": "Malgun Gothic", "font.size": 13, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#52606D", "axes.labelcolor": "#1F2937",
                     "text.color": "#17212B", "xtick.color": "#334155", "ytick.color": "#334155",
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})

def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=200, bbox_inches="tight", pad_inches=.12)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=.12)
    plt.close(fig)

def legend(fig, labels, colors, y=.99, ncol=None):
    hs = [Line2D([0], [0], color=c, lw=8, label=l) for l, c in zip(labels, colors)]
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(.5, y), ncol=ncol or len(labels),
               frameon=False, fontsize=12)

# 1. The reference graphic's case-specific parent-lineage denominator, with matched transfer labels.
rare = pd.read_csv(ROOT / "figure_realism_rare/rare_parent_share_with_new_lora.csv")
sets = {
    "sepsis": ("adult_sepsis", "Sepsis · adult cases", ["Mono_c06-CD14-MKI67", "B_c08-MKI67",
       "T_CD4_c07-MKI67", "T_CD8_c06-MKI67", "NK_c04-MKI67", "B_c07-MZB1", "Macrophage", "gdT", "ILC", "HSC"]),
    "covid": ("covid_positive", "COVID / ESKD · positive", ["CD14mono_dex", "B_ASC_IgG", "T CD4 cycling",
       "ASDC", "DC1", "pDC", "HSPC", "ILC", "T g/d", "MAIT", "NK CD56 bright"]),
}
fig, axs = plt.subplots(1, 2, figsize=(15.6, 7.3), gridspec_kw={"left":.19,"right":.99,"top":.88,"bottom":.14,"wspace":.53})
rare_sources = [("real B (transferred)", "real B · transfer", COL["real"]),
                ("floor (transferred)", "floor", COL["floor"]),
                ("new LoRA ceiling (transferred)", "new LoRA", COL["new"])]
for ax, (cohort, (grp, title, names)) in zip(axs, sets.items()):
    yy = np.arange(len(names)); w=.22
    for j, (source, _, color) in enumerate(rare_sources):
        data = rare[(rare.cohort==cohort)&(rare.group==grp)&(rare.source==source)].set_index("subset")
        vals = np.array([data.loc[n,"pct_of_parent"] for n in names], dtype=float)
        yyj = yy+(j-1)*w
        ax.barh(yyj, np.maximum(vals-.003,0), height=w*.92, left=.003, color=color,
                edgecolor="#111827" if j==2 else color, linewidth=.25,
                hatch="///" if j==2 else None, zorder=3)
        for y0, v in zip(yyj, vals):
            if v==0: ax.text(.0042,y0,"0",ha="left",va="center",color=color,fontsize=9,fontweight="bold")
    ax.set_yticks(yy, names, fontsize=11.5)
    ax.invert_yaxis(); ax.set_xscale("log"); ax.set_xlim(.003, 100)
    ax.set_xticks([.01,.1,1,10,100], ["0.01","0.1","1","10","100"])
    ax.grid(axis="x", color=COL["grid"], lw=.8, zorder=0)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=16, pad=15)
    ax.set_xlabel("% of parent lineage (case, pooled)", fontsize=12)
legend(fig,[s[1] for s in rare_sources],[s[2] for s in rare_sources],.98)
fig.text(.5,.035,"Same 30-NN transferred label route across sources; 0 is marked at the left edge. Real author labels can differ.",
         ha="center",fontsize=11,color="#475569")
save(fig,"01_rare_case_parent_new_lora")

# A four-subtype extract from the all-patient table fits the existing numerical slide.
focus=[("covid","pDC","COVID pDC / DC"),("sepsis","pDC","Sepsis pDC / DC"),
       ("sepsis","B_c07-MZB1","Sepsis B_c07-MZB1 / B"),
       ("sepsis","Mono_c06-CD14-MKI67","Sepsis Mono_c06-MKI67 / Mono")]
fig,ax=plt.subplots(figsize=(13.7,5.1),gridspec_kw={"left":.31,"right":.96,"top":.81,"bottom":.19})
for i,(cohort,sub,lab) in enumerate(focus):
    for j,(src,_,col) in enumerate(rare_sources):
        v=rare[(rare.cohort==cohort)&(rare.group=="all")&(rare.subset==sub)&(rare.source==src)].pct_of_parent.iloc[0]
        ax.scatter(max(v,.01),i+(j-1)*.19,s=145,c=col,zorder=3,
                   marker="D" if j==2 else ("s" if j==1 else "o"))
        ax.text(max(v,.01)*1.12,i+(j-1)*.19,f"{v:.3g}%",va="center",fontsize=11,color=col,fontweight="bold")
ax.set_yticks(np.arange(4),[x[2] for x in focus]);ax.invert_yaxis();ax.set_xscale("log");ax.set_xlim(.01,100)
ax.set_xticks([.01,.1,1,10,100],["0.01","0.1","1","10","100"]);ax.grid(axis="x",color=COL["grid"],zorder=0)
ax.set_xlabel("% of parent lineage (all patients pooled)",fontsize=12)
legend(fig,[s[1] for s in rare_sources],[s[2] for s in rare_sources],.98)
fig.text(.5,.055,"Matched transferred labels across sources; author annotation can disagree for rare cells.",ha="center",fontsize=11,color="#475569")
save(fig,"08_rare_focus_slide")

# 2. PCA versus Harmony iLISI. The 1–2 range is kept explicit; bars show excess over the 1.0 no-mixing floor.
h = pd.read_csv(ROOT / "historical_metric_comparison.csv")
fig, axs = plt.subplots(1, 2, figsize=(14.5,5.6), gridspec_kw={"left":.07,"right":.99,"top":.85,"bottom":.19,"wspace":.18})
for ax, cohort, title in zip(axs,["sepsis","covid"],["Sepsis","COVID / ESKD"]):
    d=h[h.cohort==cohort].set_index("source"); xx=np.arange(4); w=.35
    for j, (key,label,color) in enumerate(SOURCES):
        pre=float(d.loc[key,"iLISI_pre"]); post=float(d.loc[key,"iLISI_post"])
        ax.bar(j-w/2,pre-1,bottom=1,width=w,color="#B8C1C8",zorder=3)
        ax.bar(j+w/2,post-1,bottom=1,width=w,color=color,zorder=3)
        ax.text(j+w/2,post+.018,f"{post:.2f}",ha="center",fontsize=11,color=color,fontweight="bold")
    ax.axhline(2,ls="--",color="#536170",lw=1.1)
    ax.set_ylim(1,2.08);ax.set_xlim(-.65,3.65)
    ax.set_xticks(xx,[s[1] for s in SOURCES],fontsize=11)
    ax.set_yticks([1,1.25,1.5,1.75,2])
    ax.grid(axis="y",color=COL["grid"],zorder=0)
    ax.set_title(title,loc="left",fontsize=16,fontweight="bold",pad=16)
    ax.set_ylabel("median iLISI (real A vs source)",fontsize=12 if cohort=="sepsis" else 0)
    if cohort=="covid": ax.set_ylabel("")
fig.legend([Line2D([0],[0],lw=8,color="#B8C1C8"),Line2D([0],[0],lw=8,color=COL["new"])],
           ["PCA projection (before)","after pairwise Harmony (source color)"],loc="upper center",
           bbox_to_anchor=(.5,.99),ncol=2,frameon=False,fontsize=12)
fig.text(.5,.045,"1 = separated · 2 = ideal two-source mixing. Each source has its own real A pairwise integration; no draw-level CI.",
         ha="center",fontsize=11,color="#475569")
save(fig,"02_ilisi_pre_post_new_lora")

# 3. Readable realism figure for the numerical slide.
dc=pd.read_csv(ROOT/"figure_realism_rare/density_coverage_all_sources.csv")
c2=pd.read_csv(ROOT/"figure_realism_rare/c2st_all_sources.csv")
fig,axs=plt.subplots(1,3,figsize=(14.7,5.4),gridspec_kw={"left":.06,"right":.99,"top":.82,"bottom":.18,"wspace":.28})
for ax,metric,title in zip(axs[:2],["coverage","density"],["Coverage · PCA","Density · PCA"]):
    for i,cohort in enumerate(["sepsis","covid"]):
        for j,(src,lab,col) in enumerate([("realU","real B",COL["real"]),("floor","floor",COL["floor"]),("ceiling","new LoRA",COL["new"]) ]):
            v=dc[(dc.cohort==cohort)&(dc.space=="PCA (no correction)")&(dc.lineage=="all")&(dc.source==src)][metric].iloc[0]
            ax.bar(i+(j-1)*.23,v,width=.21,color=col,zorder=3)
            ax.text(i+(j-1)*.23,v+.013,f"{v:.2f}",ha="center",fontsize=10)
    ax.set_ylim(0,1 if metric=="coverage" else .5);ax.set_xticks([0,1],["Sepsis","COVID"])
    ax.set_title(title,loc="left",fontweight="bold",fontsize=15);ax.grid(axis="y",color=COL["grid"],zorder=0)
ax=axs[2]
for i,cohort in enumerate(["sepsis","covid"]):
    for j,(src,lab,col) in enumerate([("realU","real B",COL["real"]),("floor","floor",COL["floor"]),("ceiling","new LoRA",COL["new"]) ]):
        v=c2[(c2.cohort==cohort)&(c2.source==src)&(c2.lineage.isin(["T","NK","B","Myeloid"]))].c2st_auc.mean()
        ax.bar(i+(j-1)*.23,v,width=.21,color=col,zorder=3)
        ax.text(i+(j-1)*.23,v+.01,f"{v:.2f}",ha="center",fontsize=10)
ax.set_ylim(0,1.12);ax.set_xticks([0,1],["Sepsis","COVID"]);ax.axhline(.5,color="#64748B",ls="--",lw=1)
ax.set_title("C2ST AUC · 4 lineages mean",loc="left",fontweight="bold",fontsize=15);ax.grid(axis="y",color=COL["grid"],zorder=0)
legend(fig,["real B","floor","new LoRA"],[COL["real"],COL["floor"],COL["new"]],.98)
fig.text(.5,.035,"Coverage / density: k=5, no correction. C2ST: PCA 50, lower toward 0.5 is less distinguishable; lineage route differs by source.",ha="center",fontsize=10.5,color="#475569")
save(fig,"03_realism_slide")

# 4. Case composition, equal patient weighting, matched transferred broad labels.
comp=pd.read_csv(ROOT/"group_composition_comparison.csv")
fig,axs=plt.subplots(1,2,figsize=(14.3,5.5),gridspec_kw={"left":.16,"right":.99,"top":.82,"bottom":.18,"wspace":.22})
for ax,cohort,title in zip(axs,["sepsis","covid"],["Sepsis cases","COVID positive"]):
    d=comp[(comp.cohort==cohort)&(comp.level=="six_tr")&(comp.condition=="case")].set_index("cluster")
    cats=["T cell","B cell","Monocyte","NK cell","DC","Other"]
    y=np.arange(len(cats)); r=np.array([d.loc[k,"real_B_mean"]*100 for k in cats]); n=np.array([d.loc[k,"ceiling_mean"]*100 for k in cats])
    ax.hlines(y,r,n,color="#AAB4C0",lw=3,zorder=1)
    ax.scatter(r,y,s=95,color=COL["real"],zorder=3);ax.scatter(n,y,s=95,color=COL["new"],zorder=3)
    ax.set_yticks(y,cats);ax.invert_yaxis();ax.set_xlim(0,77);ax.set_xticks([0,20,40,60]);ax.grid(axis="x",color=COL["grid"])
    ax.set_title(title,loc="left",fontweight="bold",fontsize=16)
    ax.set_xlabel("mean patient share (%)",fontsize=12)
legend(fig,["real B","new LoRA"],[COL["real"],COL["new"]],.99)
fig.text(.5,.04,"Each sample contributes equally before group mean. Cases: adult sepsis or COVID positive; matched broad transfer labels.",ha="center",fontsize=11,color="#475569")
save(fig,"04_case_composition_slide")

# 5. Fine-cluster reproduction.
score=pd.read_csv(ROOT/"three_way/three_way_scorecard.csv")
fig,axs=plt.subplots(1,2,figsize=(13.5,5.3),gridspec_kw={"left":.07,"right":.99,"top":.84,"bottom":.17,"wspace":.21})
for ax,cohort,title in zip(axs,["sepsis","covid"],["Sepsis · 35 abundant clusters","COVID · 28 abundant clusters"]):
    d=score[score.cohort==cohort].set_index("source")
    for j,(src,lab,col) in enumerate(SC):
        v=int(d.loc[src,"n_reproduced"]); tot=int(d.loc[src,"n_abundant_real"])
        ax.bar(j,v/ tot*100,color=col,width=.65,zorder=3)
        ax.text(j,v/tot*100+2,f"{v}/{tot}",ha="center",fontsize=12,fontweight="bold")
    ax.set_xticks(range(4),[s[1] for s in SC],fontsize=11);ax.set_ylim(0,105);ax.grid(axis="y",color=COL["grid"],zorder=0)
    ax.set_title(title,loc="left",fontweight="bold",fontsize=16)
    if cohort=="sepsis": ax.set_ylabel("reproduced abundant fine clusters (%)",fontsize=12)
fig.text(.5,.04,"Threshold: real A and source ≥0.5%; source/real A 0.5–2; source ≥100 cells. Labels transferred separately.",ha="center",fontsize=11,color="#475569")
save(fig,"05_fine_cluster_slide")

# 6. Transcriptome similarity in clusters shared across all three generated conditions.
tr=pd.read_csv(ROOT/"three_way/common_cluster_transcriptome_summary.csv")
fig,axs=plt.subplots(1,2,figsize=(14.2,5.4),gridspec_kw={"left":.08,"right":.99,"top":.82,"bottom":.19,"wspace":.25})
for ax,metric,title in zip(axs,["median_pearson_low_expr","median_identity_residual_r"],["Low-expression mean r","Cluster-specific residual r"]):
    for i,cohort in enumerate(["sepsis","covid"]):
        d=tr[tr.cohort==cohort].set_index("source")
        for j,(src,lab,col) in enumerate(SC):
            v=float(d.loc[src,metric]);ax.scatter(i+(j-1.5)*.15,v,s=135,color=col,zorder=3)
            ax.text(i+(j-1.5)*.15,v+.017,f"{v:.3f}",ha="center",fontsize=10,rotation=45)
    ax.set_xticks([0,1],["Sepsis (15)","COVID (13)"])
    ax.set_xlim(-.43,1.43);ax.set_ylim(.25,1.08);ax.grid(axis="y",color=COL["grid"],zorder=0)
    ax.set_title(title,loc="left",fontweight="bold",fontsize=16)
    if metric=="median_pearson_low_expr":ax.set_ylabel("median Pearson r",fontsize=12)
legend(fig,[s[1] for s in SC],[s[2] for s in SC],.99)
fig.text(.5,.045,"Same reproduced clusters in all three generators; real A is reference. Low expression and relative cluster identity are distinct metrics.",ha="center",fontsize=10.5,color="#475569")
save(fig,"06_transcriptome_slide")

# 7. DEG recovery and precision; a scatter exposes over-calling better than recovery alone.
deg=pd.read_csv(ROOT/"disease_deg_matched/matched_deg_summary.csv")
fig,axs=plt.subplots(1,2,figsize=(14.2,5.5),gridspec_kw={"left":.08,"right":.99,"top":.84,"bottom":.19,"wspace":.25})
for ax,cohort,title in zip(axs,["sepsis","covid"],["Sepsis · disease vs HC","COVID positive vs negative"]):
    d=deg[deg.cohort==cohort]
    for src,lab,col in SC[:3]:
        dd=d[d.source==src]
        ax.scatter(dd.recall*100,dd.precision*100,s=100,color=col,alpha=.88,label=lab,zorder=3)
        if src=="ceiling":
            for _,r in dd.iterrows():
                if r.celltype in (["T cell","B cell","NK cell","Monocyte"] if cohort=="sepsis" else ["T cell","Monocyte"]):
                    ax.annotate(r.celltype,(r.recall*100,r.precision*100),xytext=(5,5),textcoords="offset points",fontsize=10,color=col)
    ax.set_xlim(0,84);ax.set_ylim(0,64);ax.grid(color=COL["grid"],zorder=0)
    ax.set_xlabel("real DEG recall (%)",fontsize=12)
    if cohort=="sepsis":ax.set_ylabel("generated DEG precision (%)",fontsize=12)
    ax.set_title(title,loc="left",fontweight="bold",fontsize=16)
legend(fig,[s[1] for s in SC[:3]],[s[2] for s in SC[:3]],.99)
fig.text(.5,.04,"Pseudobulk Welch + BH, q<0.05 and |log2FC|≥1. COVID repeated patients make these q-values provisional.",ha="center",fontsize=11,color="#475569")
save(fig,"07_disease_deg_slide")

spec={"intent":"Add actual figures to numeric-only deck slides and reproduce two user-supplied chart formats for New LoRA.",
      "audience":"PBISC co-investigators, Korean 30-minute talk","figure_size":"slide wide, 13–15.6 in; 200 dpi PNG and vector PDF",
      "data_sources":["historical_metric_comparison.csv","figure_realism_rare/rare_parent_share_with_new_lora.csv",
        "figure_realism_rare/density_coverage_all_sources.csv","figure_realism_rare/c2st_all_sources.csv",
        "group_composition_comparison.csv","three_way/three_way_scorecard.csv",
        "three_way/common_cluster_transcriptome_summary.csv","disease_deg_matched/matched_deg_summary.csv"],
      "comparability":"New LoRA file has K_cond=1000 and cells_per_bulk=4000; it is a single new protocol output, not a K=1000 vs K=4000 ablation.",
      "annotation":"Rare parent comparison uses transferred labels for real B and both generators; author-label discrepancies remain in original figure.",
      "iLISI":"Pairwise real A integration, PCA before versus Harmony after; local two-source ideal 2, separated 1.",
      "limitations":"One draw; no uncertainty intervals. COVID repeated patients compromise independent Welch q-values."}
(OUT/"visual_spec.json").write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding="utf-8")
print("generated",len(list(OUT.glob("*.png"))),"PNGs in",OUT)
