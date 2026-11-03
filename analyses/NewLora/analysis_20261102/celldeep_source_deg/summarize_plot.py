"""Summarize CellDEEP source DEG and render comparison/volcano figures."""
from pathlib import Path
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

HERE = Path(__file__).resolve().parent
RESULT = Path(os.environ.get("PBISC_CELLDEEP_OUTPUT_ROOT", HERE)).resolve() / "results"
INK = "#1B2733"
MUTED = "#546477"
GRID = "#D9E0E6"
CELL = "#008C78"
PB = "#364A5D"
ONLY = "#E87943"
TYPES = ["T cell", "B cell", "NK cell", "Monocyte"]
LABEL = {"adult_HC": "Healthy", "adult_ICU_ctrl": "ICU control",
         "adult_sepsis": "Sepsis", "covid_negative": "COVID negative",
         "covid_positive": "COVID positive"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#91A0AD", "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "figure.facecolor": "white", "axes.facecolor": "white", "pdf.fonttype": 42})


def read_result(row):
    return pd.read_csv(RESULT / f"{row.stratum}_deg.csv.gz")


def verify_and_recurrence(summary):
    records = []
    confirmed_tables = []
    for row in summary.itertuples():
        tab = read_result(row)
        assert len(tab) > 5000 and tab.gene.is_unique
        assert np.isfinite(tab.avg_log2FC).sum() > 4000
        assert int(tab.deg_CellDEEP.sum()) == int(row.n_DEG_CellDEEP)
        assert int(tab.deg_pseudobulk.sum()) == int(row.n_DEG_pseudobulk)
        overlap = tab.deg_CellDEEP & tab.deg_pseudobulk & (np.sign(tab.avg_log2FC) == np.sign(tab.pb_logFC))
        assert int(overlap.sum()) == int(row.n_overlap_same_direction)
        confirmed_tables.append(tab.loc[overlap, ["cohort", "condition", "celltype", "gene",
            "avg_log2FC", "q_BH", "pb_logFC", "pb_q_BH"]])
        for v in tab.loc[overlap, ["gene", "avg_log2FC", "pb_logFC"]].itertuples(index=False):
            records.append({"cohort": row.cohort, "gene": v.gene,
                            "stratum": row.stratum, "condition": row.condition,
                            "celltype": row.celltype, "direction": "New higher" if v.pb_logFC > 0 else "New lower"})
    rec = pd.DataFrame(records)
    pd.concat(confirmed_tables, ignore_index=True).sort_values(
        ["cohort", "condition", "celltype", "pb_q_BH"]
    ).to_csv(RESULT / "confirmed_degs.csv", index=False)
    if len(rec):
        out = rec.groupby(["cohort", "gene", "direction"]).agg(
            n_strata=("stratum", "nunique"),
            n_celltypes=("celltype", "nunique"),
            conditions=("condition", lambda x: ";".join(sorted(set(x)))),
            celltypes=("celltype", lambda x: ";".join(sorted(set(x))))).reset_index()
        out.sort_values(["cohort", "n_strata", "gene"], ascending=[True, False, True]).to_csv(
            RESULT / "recurrent_confirmed_genes.csv", index=False)
    return rec


def summary_plot(summary):
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [1.15, 1]})
    for ax, cohort, title in zip(axes, ["sepsis", "covid"], ["Sepsis · Ye 2026", "COVID-19 · Imperial ESKD"]):
        sub = summary[summary.cohort == cohort].copy()
        order = {"sepsis": ["adult_HC", "adult_ICU_ctrl", "adult_sepsis"],
                 "covid": ["covid_negative", "covid_positive"]}[cohort]
        sub["gorder"] = sub.condition.map({g: i for i, g in enumerate(order)})
        sub["torder"] = sub.celltype.map({t: i for i, t in enumerate(TYPES)})
        sub = sub.sort_values(["gorder", "torder"])
        yy = np.arange(len(sub))[::-1]
        labels = [f"{LABEL[g]}  |  {t}" for g, t in zip(sub.condition, sub.celltype)]
        ax.hlines(yy, sub.n_DEG_pseudobulk, sub.n_DEG_CellDEEP, color="#B7C4CC", lw=2, zorder=1)
        ax.scatter(sub.n_DEG_CellDEEP, yy, color=CELL, s=48, label="CellDEEP metacells", zorder=3)
        ax.scatter(sub.n_DEG_pseudobulk, yy, color=PB, s=45, marker="s", label="Donor pseudobulk", zorder=3)
        for y, cd, pb, n in zip(yy, sub.n_DEG_CellDEEP, sub.n_DEG_pseudobulk, sub.n_real_samples):
            ax.text(cd + 10, y, str(int(cd)), va="center", fontsize=8, color=CELL)
            ax.text(max(float(pb)-10,0), y-.23, str(int(pb)), va="top", ha="right", fontsize=7.5, color=PB)
        ax.set_yticks(yy, labels, fontsize=9)
        ax.set_xlim(0, 850)
        ax.set_xlabel("Genes with BH q < 0.05 and |log₂FC| ≥ 1", fontsize=9)
        ax.set_title(title, loc="left", fontsize=14, fontweight="bold", color=INK, pad=14)
        ax.grid(axis="x", color=GRID, lw=.7)
        ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, loc="upper right", bbox_to_anchor=(.96, .985), frameon=False, fontsize=8)
    fig.suptitle("Real vs New LoRA: metacell DEG calls and donor-level check", x=.28, y=.99,
                 fontsize=17, fontweight="bold", color=INK)
    fig.text(.04,.025,"Each contrast is within one condition and one broad cell type; positive log₂FC means higher in New LoRA. Sample counts and full gene tables are in summary.csv.",
             fontsize=8.2,color=MUTED)
    fig.subplots_adjust(left=.17,right=.97,bottom=.09,top=.88,wspace=.52)
    for ext in ("png", "pdf"):
        fig.savefig(RESULT / f"celldeep_vs_pseudobulk_summary.{ext}", dpi=240)
    plt.close(fig)


def volcano_plot(summary, cohort, group):
    sub = summary[(summary.cohort == cohort) & (summary.condition == group)]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    for ax, ct in zip(axes.ravel(), TYPES):
        row = sub[sub.celltype == ct]
        if row.empty:
            ax.axis("off")
            continue
        row = row.iloc[0]
        tab = read_result(row)
        x = tab.avg_log2FC.to_numpy(dtype=float)
        y = -np.log10(np.clip(tab.q_BH.to_numpy(dtype=float), 1e-300, 1))
        confirmed = (tab.deg_CellDEEP & tab.deg_pseudobulk &
                     (np.sign(tab.avg_log2FC) == np.sign(tab.pb_logFC))).to_numpy(dtype=bool)
        cd_only = tab.deg_CellDEEP.to_numpy(dtype=bool) & ~confirmed
        rest = ~tab.deg_CellDEEP.to_numpy(dtype=bool)
        valid = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[valid & rest], y[valid & rest], s=4, color="#D1D9DF", alpha=.55, lw=0, rasterized=True)
        ax.scatter(x[valid & cd_only], y[valid & cd_only], s=7, color=ONLY, alpha=.65, lw=0, rasterized=True)
        ax.scatter(x[valid & confirmed], y[valid & confirmed], s=12, color=CELL, alpha=.85, lw=0, rasterized=True)
        ax.axvline(-1, color="#A9B6C0", ls="--", lw=.8)
        ax.axvline(1, color="#A9B6C0", ls="--", lw=.8)
        ax.axhline(-np.log10(.05), color="#A9B6C0", ls="--", lw=.8)
        labels = tab.loc[confirmed & valid].sort_values("pb_q_BH").head(3)
        texts = [ax.text(gene.avg_log2FC, -np.log10(max(gene.q_BH, 1e-300)), gene.gene,
                         fontsize=7.2, color=INK) for gene in labels.itertuples()]
        ax.set_title(f"{ct}  ·  CellDEEP {int(row.n_DEG_CellDEEP)}  /  pseudobulk overlap {int(row.n_overlap_same_direction)}",
                     loc="left", fontsize=11.5, fontweight="bold", color=INK)
        ax.grid(axis="y", color=GRID, lw=.5)
        ax.set_xlim(-6,6)
        ax.set_ylim(0, 17)
        if texts:
            adjust_text(texts, ax=ax, ensure_inside_axes=True,
                        force_text=(.45, .65), force_static=(.15, .25),
                        expand=(1.1, 1.2),
                        arrowprops={"arrowstyle": "-", "color": "#9AA8B2", "lw": .5})
    fig.supxlabel("CellDEEP log₂FC (New LoRA − Real)", fontsize=10, y=.035, color=MUTED)
    fig.supylabel("−log₁₀ BH q", fontsize=10, x=.025, color=MUTED)
    fig.suptitle(f"{'COVID-19' if cohort == 'covid' else 'Sepsis'} · {LABEL[group]} · source-associated genes", x=.10, y=.99,
                 ha="left", fontsize=17, fontweight="bold", color=INK)
    fig.text(.12,.06,"Teal: CellDEEP and donor pseudobulk both significant, same direction   ·   Orange: CellDEEP only   ·   Gray: below CellDEEP threshold",
             fontsize=8.4, color=MUTED)
    fig.subplots_adjust(left=.10,right=.98,bottom=.12,top=.92,wspace=.2,hspace=.27)
    for ext in ("png", "pdf"):
        fig.savefig(RESULT / f"volcano_{cohort}_{group}.{ext}", dpi=240)
    plt.close(fig)


if __name__ == "__main__":
    summary = pd.read_csv(RESULT / "summary.csv")
    assert len(summary) == 19 and summary.error.isna().all() if "error" in summary else len(summary)==19
    verify_and_recurrence(summary)
    summary_plot(summary)
    volcano_plot(summary, "sepsis", "adult_sepsis")
    volcano_plot(summary, "covid", "covid_positive")
