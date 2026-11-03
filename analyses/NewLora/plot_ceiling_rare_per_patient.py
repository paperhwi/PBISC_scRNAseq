"""Recreate the prior rare-population dotplots for the new LoRA ceiling output.

Real points and rare-population/parent definitions come from the original
2026-10-31 analysis. Ceiling cells are the 60k balanced subset used for the
real-A Harmony comparison, with 30-NN-transferred labels.
"""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RARE_BASE = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027"
OLD = RARE_BASE / "rare_deg_gsea_20261031/rare_unit_shares.csv"
sys.path.insert(0, str(RARE_BASE / "code"))
from rare_lib import SEPSIS, COVID
from deg_lib import bh

OUT = ROOT / "NewLora/analysis_20261102/rare_populations"
OUT.mkdir(parents=True, exist_ok=True)
GRP = {"sepsis": {"adult_sepsis": "case", "adult_HC": "control", "adult_ICU_ctrl": "ICU control"},
       "covid": {"covid_positive": "case", "covid_negative": "control"}}
GREY, NEW, INK, RED = "#7F7F7F", "#D45D48", "#1F2933", "#D7191C"

real = pd.read_csv(OLD)
real = real[real.source == "real"].copy()
real["source"] = "real"
rows = []
for cohort, definitions in (("sepsis", SEPSIS), ("covid", COVID)):
    labels = pd.read_csv(ROOT / "NewLora/analysis_20261102" / cohort / "ceiling_cell_labels.csv.gz",
                         usecols=["sample", "group", "major_tr", "fine_tr"])
    labels["condition"] = labels.group.map(GRP[cohort])
    labels = labels[labels.condition.notna()].copy()
    labels["patient"] = labels["sample"] if cohort == "sepsis" else labels["sample"].str.split("_").str[0]
    for (patient, condition), group in labels.groupby(["patient", "condition"], sort=True):
        major = group.major_tr.astype(str).to_numpy()
        fine = group.fine_tr.astype(str).to_numpy()
        n = len(group)
        for subset, parent_fn in definitions.items():
            parent = parent_fn(major, fine)
            k = int((parent & (fine == subset)).sum())
            nparent = int(parent.sum())
            rows.append(dict(cohort=cohort, source="ceiling", patient=patient, group=condition,
                             unit=f"ceiling|{patient}|{condition}", subset=subset, n_cells=n,
                             n_subset=k, n_parent=nparent, pct_all=100*k/n,
                             pct_parent=100*k/nparent if nparent >= 20 else np.nan))

ceiling = pd.DataFrame(rows)
unit = pd.concat([real[ceiling.columns], ceiling], ignore_index=True)
unit.to_csv(OUT / "rare_unit_shares_real_vs_ceiling.csv", index=False)

test_rows = []
for (cohort, subset, condition), data in unit.groupby(["cohort", "subset", "group"], sort=False):
    r = data[data.source == "real"]
    c = data[data.source == "ceiling"]
    if len(r) < 2 or len(c) < 2:
        continue
    mr, mc = float(r.pct_all.mean()), float(c.pct_all.mean())
    sdr, sdc = float(r.pct_all.std(ddof=1)), float(c.pct_all.std(ddof=1))
    with np.errstate(invalid="ignore"):
        p = float(stats.ttest_ind(r.pct_all, c.pct_all, equal_var=False).pvalue)
    test_rows.append(dict(cohort=cohort, subset=subset, group=condition,
                          n_real=len(r), n_ceiling=len(c), mean_real=mr, mean_ceiling=mc,
                          sd_real=sdr, sd_ceiling=sdc, sd_ratio=sdc/sdr if sdr>0 else np.nan,
                          fold_ceiling_vs_real=mc/mr if mr>0 else np.nan, welch_p=p))

tests = pd.DataFrame(test_rows)
tests["welch_q"] = np.nan
for _, ix in tests.groupby(["cohort", "group"]).groups.items():
    values = tests.loc[ix, "welch_p"]
    tests.loc[ix, "welch_q"] = bh(values.to_numpy())
tests.to_csv(OUT / "rare_stats_real_vs_ceiling.csv", index=False)

plt.rcParams.update({"font.family":"DejaVu Sans", "axes.spines.top":False,
                     "axes.spines.right":False, "font.size":8, "axes.titlesize":8.5})
def stars(q):
    if not np.isfinite(q):return "ns"
    return "***" if q<.001 else "**" if q<.01 else "*" if q<.05 else "ns"

for cohort, definitions in (("sepsis", SEPSIS), ("covid", COVID)):
    groups = ["control", "ICU control", "case"] if cohort == "sepsis" else ["control", "case"]
    labels = {"control":"HC" if cohort=="sepsis" else "neg", "ICU control":"ICU",
              "case":"sepsis" if cohort=="sepsis" else "pos"}
    subsets = list(definitions)
    ncols = 5
    nrows = int(np.ceil(len(subsets)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.2, 2.75*nrows))
    axes = axes.ravel()
    rng = np.random.default_rng(0)
    for ax, subset in zip(axes, subsets):
        d = unit[(unit.cohort==cohort) & (unit.subset==subset)]
        ymax = float(d.pct_all.max())
        for gi, condition in enumerate(groups):
            for si, (source, color) in enumerate((("real", GREY), ("ceiling", NEW))):
                values = d[(d.group==condition)&(d.source==source)].pct_all.to_numpy()
                x = gi*2.6 + si
                ax.scatter(x+rng.uniform(-.28,.28,len(values)),values,s=7,color=color,
                           alpha=.75,lw=0,zorder=2)
                if len(values)>1:
                    ax.errorbar(x+.42, values.mean(), yerr=values.std(ddof=1), fmt="o",
                                ms=3.5, color=INK, capsize=2.5, lw=1, zorder=3)
            st = tests[(tests.cohort==cohort)&(tests.subset==subset)&(tests.group==condition)]
            if len(st):
                s=st.iloc[0]
                sd_text=f"SD×{s.sd_ratio:.2f}" if np.isfinite(s.sd_ratio) else "SD×NA"
                ax.text(gi*2.6+.5, ymax*1.13 if ymax>0 else .78,
                        f"{stars(s.welch_q)}\n{sd_text}",ha="center",va="bottom",fontsize=6.3,
                        color=RED if np.isfinite(s.welch_q) and s.welch_q<.05 else INK)
        ax.set_ylim(0,ymax*1.45 if ymax>0 else 1)
        ax.set_xticks([gi*2.6+.5 for gi in range(len(groups))])
        ax.set_xticklabels([labels[g] for g in groups])
        ax.set_title(subset,loc="left")
        ax.set_ylabel("% of all cells",fontsize=7)
    for ax in axes[len(subsets):]:ax.axis("off")
    legend=[Line2D([],[],marker="o",ls="",color=GREY,label="real — one dot per patient (author label)"),
            Line2D([],[],marker="o",ls="",color=NEW,label="new LoRA ceiling — one dot per generated patient (kNN label)"),
            Line2D([],[],marker="o",ls="-",color=INK,label="mean ± SD"),
            Line2D([],[],ls="",label="* / ** / *** = Welch t BH-q < 0.05 / 0.01 / 0.001 ; SD× = ceiling SD / real SD")]
    fig.legend(handles=legend,loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(.5,-.01),fontsize=8)
    fig.suptitle(("Sepsis (Ye 2026)" if cohort=="sepsis" else "COVID / ESKD (Imperial)")+
                 " — rare populations per patient · new LoRA ceiling",x=.005,ha="left",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,.05,1,.97])
    fig.savefig(OUT/f"rare_dotplot_{cohort}_ceiling.png",dpi=240,bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/f"rare_dotplot_{cohort}_ceiling.pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)

(OUT/"README.md").write_text("""# New LoRA ceiling rare-population per-patient plots

The real points are the exact author-labeled observations from the prior
`rare_deg_gsea_20261031/rare_unit_shares.csv`. New LoRA ceiling cells use the
60,000-cell balanced subsample and real-A 30-NN-transferred fine/major labels
from `NewLora/analysis_20261102/{cohort}/ceiling_cell_labels.csv.gz`.

Sepsis unit is sample/patient. COVID unit is patient × condition with visits
pooled, matching the prior plot. Rare subset and parent definitions come from
`disease_deg_rare_20261027/code/rare_lib.py`. The ordinate is subset cells as
a percentage of *all* cells in a unit. Dots are individual units, with
mean ± sample SD. Welch t tests compare real and ceiling within each
condition; BH FDR is applied across the plotted subsets separately for each
cohort and condition. SD× is ceiling sample SD divided by real sample SD;
NA means the real SD is zero.

The ceiling file has no native cell-type annotations, so these proportions
depend on the transferred labels. Small percentages are sensitive to the
~500 generated cells sampled per sepsis patient and varying cells per COVID
patient; zeros need not imply biological absence. The plot gives a descriptive
comparison of one generated draw per cohort and does not treat cells as
independent biological replicates.
""",encoding="utf-8")
print(OUT)
