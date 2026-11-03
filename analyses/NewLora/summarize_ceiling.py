"""Summarize group-matched sample composition and historic pairwise iLISI."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "NewLora/analysis_20261102"
OLD = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/results/big"
GMAP = {"adult_sepsis":"case", "adult_HC":"control", "adult_ICU_ctrl":"icu_control",
        "covid_positive":"case", "covid_negative":"control"}
all_metrics=[]
all_group=[]
def markdown_table(df):
    vals = df.astype(str)
    head = "| " + " | ".join(vals.columns) + " |"
    rule = "| " + " | ".join(["---"]*len(vals.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in vals.itertuples(index=False, name=None)]
    return "\n".join([head,rule]+rows)
for coh in ("sepsis","covid"):
    folder = OUT/coh
    result = json.loads((folder/"summary.json").read_text(encoding="utf-8"))
    rows = pd.read_csv(folder/"composition_by_sample.csv")
    # Use the historical held-out real-B transfer for fine labels, matching the
    # exact reference and label route used for the new generated cells.
    transfer = pd.read_csv(OLD/f"{coh}_realU_transfer.csv.gz", index_col=0)
    fine_types = sorted(set(transfer.fine_tr.astype(str)) | set(rows.loc[(rows.source=="ceiling") & (rows.level=="fine_tr"),"cluster"]))
    fine_counts = pd.crosstab(transfer["sample"], transfer["fine_tr"]).reindex(columns=fine_types, fill_value=0)
    fine_rates = fine_counts.div(fine_counts.sum(1), axis=0)
    fine_rates.columns.name = "cluster"
    new_fine = fine_rates.stack().rename("fraction").reset_index()
    new_fine["group"] = new_fine["sample"].map(transfer.drop_duplicates("sample").set_index("sample").group)
    new_fine["source"] = "real_B"
    new_fine["level"] = "fine_tr"
    rows = pd.concat([rows.loc[~((rows.source=="real_B") & (rows.level=="fine_tr"))], new_fine[rows.columns]], ignore_index=True)
    rows.to_csv(folder/"composition_by_sample.csv", index=False)
    rows.groupby(["level","source","cluster"]).fraction.agg(["mean","std","count"]).reset_index().to_csv(folder/"composition_summary.csv",index=False)
    rows["condition"] = rows.group.map(GMAP)
    rows = rows[rows.condition.notna()]
    for level in ("six_tr","fine_tr"):
        part = rows[rows.level==level]
        types = sorted(part.cluster.unique())
        for cond in (["case","control","icu_control"] if coh=="sepsis" else ["case","control"]):
            p = part[part.condition==cond]
            if p.empty: continue
            means = p.groupby(["source","cluster"]).fraction.mean().unstack("source").reindex(types).fillna(0)
            counts = p.groupby("source")["sample"].nunique().to_dict()
            for ty, rec in means.iterrows():
                all_group.append({"cohort":coh,"level":level,"condition":cond,"cluster":ty,
                                  "real_B_mean":rec.get("real_B",0),"ceiling_mean":rec.get("ceiling",0),
                                  "delta_pp":100*(rec.get("ceiling",0)-rec.get("real_B",0)),
                                  "n_real_B_samples":counts.get("real_B",0),"n_ceiling_samples":counts.get("ceiling",0)})
            if level=="six_tr":
                result[f"six_TV_{cond}_vs_real_B"] = float(.5*np.abs(means.get("ceiling",0)-means.get("real_B",0)).sum())
        if level=="six_tr":
            p=part.pivot_table(index=["source","sample"],columns="cluster",values="fraction",fill_value=0)
            # Group-level contrasts use equal sample weights and no cell-level p values.
            g=part[["source","sample","condition"]].drop_duplicates().set_index(["source","sample"]).condition
            for ctrl in (["control","icu_control"] if coh=="sepsis" else ["control"]):
                dr=p.loc[g=="case"].loc["real_B"].mean()-p.loc[g==ctrl].loc["real_B"].mean()
                dg=p.loc[g=="case"].loc["ceiling"].mean()-p.loc[g==ctrl].loc["ceiling"].mean()
                keep=dr.index!="Other"
                result[f"six_shift_r_case_vs_{ctrl}"]=float(np.corrcoef(dr[keep],dg[keep])[0,1])
                pd.DataFrame({"real_B_shift":dr,"ceiling_shift":dg}).to_csv(folder/f"six_shift_case_vs_{ctrl}.csv",index_label="six_type")
    (folder/"summary.json").write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    for source in ("realU","prev4k","floor"):
        old=json.loads((OLD/f"{coh}_{source}_summary.json").read_text())
        all_metrics.append({"cohort":coh,"source":source,"iLISI_pre":old["iLISI_pre"],"iLISI_post":old["iLISI_post"],
                            "cluster_r":old["cluster_freq_pearson"],"cluster_JSD":old["cluster_freq_JSD"],
                            "fraction_in_mixed_clusters":old["src_cells_in_mixed_clusters"]})
    all_metrics.append({"cohort":coh,"source":"new_ceiling","iLISI_pre":result["iLISI_pre"],"iLISI_post":result["iLISI_post"],
                        "cluster_r":result["cluster_fraction_pearson"],"cluster_JSD":result["cluster_fraction_JSD"],
                        "fraction_in_mixed_clusters":result["ceiling_fraction_in_mixed_clusters"]})

group=pd.DataFrame(all_group); group.to_csv(OUT/"group_composition_comparison.csv",index=False)
pd.DataFrame(all_metrics).to_csv(OUT/"historical_metric_comparison.csv",index=False)

fig,axes=plt.subplots(2,2,figsize=(12,8.5))
for i,coh in enumerate(("sepsis","covid")):
    d=group[(group.cohort==coh)&(group.level=="six_tr")]
    for j,cond in enumerate(("case","control")):
        p=d[d.condition==cond].sort_values("cluster")
        x=np.arange(len(p)); w=.38
        axes[i,j].bar(x-w/2,100*p.real_B_mean,width=w,color="#66788E",label="real B")
        axes[i,j].bar(x+w/2,100*p.ceiling_mean,width=w,color="#D45D48",label="new ceiling")
        axes[i,j].set_xticks(x,p.cluster,rotation=35,ha="right")
        axes[i,j].set_ylabel("Mean sample share (%)")
        axes[i,j].set_title(f"{coh} · {cond} · real B {int(p.n_real_B_samples.iloc[0])} / ceiling {int(p.n_ceiling_samples.iloc[0])} samples")
        axes[i,j].legend(frameon=False,fontsize=8)
fig.tight_layout();fig.savefig(OUT/"six_type_composition.png",dpi=220);plt.close(fig)

with (OUT/"REPORT.md").open("w",encoding="utf-8") as f:
    f.write("# New LoRA ceiling vs real: Harmony, iLISI and cell composition\n\n")
    f.write("> 2026-11-02 검증 기록: 계산 재현과 희소 라벨·COVID DEG의 해석 제한은 [검증 보고서](verification_20261102/검증_보고서.md)에 정리했다.\n\n")
    f.write("## Analysis design\n\n")
    f.write("- Real A (previous patient split) is the reference; real B is the held-out observed comparison.\n")
    f.write("- Stratified 60,000 generated cells per cohort by sample, with adult-only sepsis and primary-analysis-eligible COVID cells. The COVID filter removes five support-donor-overlap samples.\n")
    f.write("- Aligned gene symbols to the existing measured panels (sepsis 5,961; COVID 5,745). Real A count-based 2,000 Seurat v3 HVGs, log1p CP10k, real-A scaling and 30-PC PCA; pairwise Harmony on source (max 20 iterations).\n")
    f.write("- iLISI uses perplexity 30 and 90 neighbors; reported value is the mean of three 10,000-query median iLISI estimates. Two-source range is 1–2, higher indicates better mixing.\n")
    f.write("- Fine and six-type labels in new cells were transferred from real A by 30 nearest neighbors in Harmony space. Real B six-type labels are author labels; real B fine labels are the prior 30-neighbor transfers. Six/fine composition summaries give equal weight per sample.\n")
    f.write("- Joint Leiden (resolution 1) was fit to real A plus new cells. Cluster frequencies are cell fractions for those same two groups.\n\n")
    f.write("## Integration results\n\n")
    f.write(markdown_table(pd.DataFrame(all_metrics).round(3))+"\n\n")
    for coh in ("sepsis","covid"):
        r=json.loads((OUT/coh/"summary.json").read_text())
        f.write(f"## {coh}\n\n")
        f.write(f"- New ceiling iLISI: {r['iLISI_pre']:.3f} before → {r['iLISI_post']:.3f} after Harmony.\n")
        f.write(f"- Joint Leiden: {r['n_leiden_clusters']} clusters; fraction correlation {r['cluster_fraction_pearson']:.3f}, JSD {r['cluster_fraction_JSD']:.3f}; {100*r['ceiling_fraction_in_mixed_clusters']:.1f}% of generated cells are in clusters with normalized source share 0.2–0.8.\n")
        for cond in (["case","control","icu_control"] if coh=="sepsis" else ["case","control"]):
            f.write(f"- Six-type total variation vs real B, {cond}: {r[f'six_TV_{cond}_vs_real_B']:.3f}.\n")
        f.write(f"- Six-type case–control shift correlation vs real B: {r['six_shift_r_case_vs_control']:.3f}.\n\n")
        t=group[(group.cohort==coh)&(group.level=="six_tr")&group.condition.isin(["case","control"])].copy()
        for cond in ("case","control"):
            x=t[t.condition==cond][["cluster","real_B_mean","ceiling_mean","delta_pp"]].copy()
            x[["real_B_mean","ceiling_mean"]]*=100
            x.columns=["Cell type","real B %","ceiling %","Δ pp"]
            f.write(f"### {cond} six-type sample means\n\n{markdown_table(x.round(2))}\n\n")
        c=pd.read_csv(OUT/coh/"joint_leiden_clusters.csv")
        top=c.reindex(c.delta_pp.abs().sort_values(ascending=False).index).head(8)[["leiden","dominant_real_fine","fraction_real_A","fraction_ceiling","delta_pp"]]
        top[["fraction_real_A","fraction_ceiling"]]*=100
        top.columns=["Leiden","Dominant real fine label","real A %","ceiling %","Δ pp"]
        f.write("### Largest joint Leiden differences\n\n"+markdown_table(top.round(2))+"\n\n")
    f.write("## Interpretation limits\n\n")
    f.write("- Harmony corrects the embedding and neighbor graph; original count matrices are unchanged. iLISI alone does not establish biological fidelity.\n")
    f.write("- New inputs have no cell-type annotations. Their composition depends on real-A label transfer; the historical floor/previous composition metrics used frozen generated labels, so those composition scores are not directly interchangeable.\n")
    f.write("- Real B contains fewer samples than the bulk-derived cohorts, especially sepsis controls. Composition differences are descriptive and patient/sample is the unit; no cell-level significance claims are made.\n")
    f.write("- New sepsis uses draw 1, COVID draw 0; both are single draws.\n\n")
    f.write("## Files\n\n")
    f.write("- `historical_metric_comparison.csv`, `group_composition_comparison.csv`, `six_type_composition.png` at this level.\n")
    f.write("- Each cohort folder contains `summary.json`, `composition_by_sample.csv`, `composition_summary.csv`, `joint_leiden_clusters.csv`, `ceiling_cell_labels.csv.gz`, `joint_harmony_embedding.h5ad`, `umap_harmony.png`, and disease-shift tables.\n")
print(OUT/"REPORT.md")
