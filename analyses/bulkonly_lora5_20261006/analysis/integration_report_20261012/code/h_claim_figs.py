"""Figures for parent-lineage denominators: claim populations control -> case, and rare subsets within parent."""
import logging, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
R = "/home/claude/work/report"; F = f"{R}/figs"
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 220})
SRC = [("real", "Real", "#555555"), ("k1000", "K = 1,000", "#2b83ba"), ("k4000", "K = 4,000", "#d7191c")]
d = pd.read_csv(f"{R}/rare_parent_denominators.csv")


def claim_fig():
    rows = [("sepsis", "T_CD4_c02-GPR183-NR4A2", "SEP-1  NR4A2+ CD4 T  (in CD4 T)"), ("sepsis", "T_CD8_c05-GNLY-CCL4hi", "SEP-2  GNLY/CCL4 CD8 T  (in CD8 T)"),
            ("sepsis", "NKT_c05-IFNG-CCL4hi", "SEP-3  IFNG/CCL4 NKT  (in NKT)"), ("sepsis", "NK_c03-PRF1-CCL4hi", "SEP-3-NK  PRF1/CCL4 NK  (in NK)"),
            ("sepsis", "Mono_c01-CD14-HLA-DRA", "SEP-4A  HLA-DR CD14 mono  (in mono)"), ("sepsis", "Mono_c03-CD14-RETN", "SEP-4B  RETN CD14 mono  (in mono)"),
            ("sepsis", "Mono_c06-CD14-MKI67", "SEP-5  MKI67 CD14 mono  (in mono), adult"),
            ("covid", "CD14mono_IFN", "COVID  IFN CD14 mono  (in CD14 mono)"), ("covid", "CD14mono_dex", "COVID  dex CD14 mono  (in CD14 mono)")]
    case = {"sepsis": "adult_sepsis", "covid": "covid_positive"}; ctrl = {"sepsis": "adult_HC", "covid": "covid_negative"}
    fig, ax = plt.subplots(figsize=(10, 5.6)); y0 = 0; yt = []
    for coh, sub, lab in rows[::-1]:
        a = d[(d.cohort == coh) & (d.subset == sub) & (d.group == ctrl[coh])].iloc[0]; b = d[(d.cohort == coh) & (d.subset == sub) & (d.group == case[coh])].iloc[0]
        for j, (s, n, c) in enumerate(SRC):
            y = y0 + (1 - j) * 0.24; x0, x1 = max(a[f"{s}_pooled_pct"], 0.03), max(b[f"{s}_pooled_pct"], 0.03)
            low = (s != "real") and (a[f"{s}_n_parent"] < 300 or b[f"{s}_n_parent"] < 300)
            ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.6, alpha=0.45 if low else 1, mutation_scale=10))
            ax.plot([x0], [y], "o", color=c, ms=3.5, alpha=0.45 if low else 1)
            if low: ax.text(max(x0, x1) * 1.15, y, "few parent cells", fontsize=6.5, color=c, va="center")
        yt.append((y0, lab)); y0 += 1
    ax.set_xscale("log"); ax.set_xlim(0.025, 100); ax.set_yticks([p for p, _ in yt], [l for _, l in yt], fontsize=8.5)
    ax.set_xlabel("% of parent lineage (pooled cells; dot = control, arrowhead = case; 0 drawn at 0.03%)")
    ax.legend(handles=[plt.Line2D([], [], color=c, lw=2, label=n) for _, n, c in SRC], loc="lower right", fontsize=8)
    ax.set_title("Paper-claim populations: control → case, within parent lineage (generated labels by kNN transfer)", loc="left", fontsize=10)
    ax.grid(axis="x", color="#eeeeee", lw=0.6); plt.tight_layout(); fig.savefig(f"{F}/claims_parent.png"); plt.close(fig)


def rare_fig():
    sel = {"sepsis": ["Mono_c06-CD14-MKI67", "B_c08-MKI67", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67", "NK_c04-MKI67", "B_c07-MZB1", "Macrophage", "gdT", "ILC", "HSC"],
           "covid": ["CD14mono_dex", "B_ASC_IgG", "T CD4 cycling", "ASDC", "DC1", "pDC", "HSPC", "ILC", "T g/d", "MAIT", "NK CD56 bright"]}
    case = {"sepsis": "adult_sepsis", "covid": "covid_positive"}
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.6))
    for ax, coh in zip(axs, ["sepsis", "covid"]):
        t = d[(d.cohort == coh) & (d.group == case[coh])].set_index("subset").loc[sel[coh]][::-1]; y = np.arange(len(t))
        for j, (s, n, c) in enumerate(SRC):
            v = t[f"{s}_pooled_pct"].clip(lower=0.003)
            ax.barh(y + (1 - j) * 0.26, v, 0.26, color=c, label=n)
            for yi, (vv, nt) in enumerate(zip(t[f"{s}_pooled_pct"], t[f"{s}_n_target"])):
                if s != "real" and nt == 0: ax.text(0.0035, yi + (1 - j) * 0.26, "0", fontsize=6, color=c, va="center")
        lab = [f"{i}  [{p}]" for i, p in zip(t.index, t.parent)]
        ax.set_yticks(y, lab, fontsize=7.5); ax.set_xscale("log"); ax.set_xlim(0.003, 100)
        ax.set_xlabel(f"% of parent lineage, {case[coh].replace('_', ' ')} (pooled)"); ax.set_title({"sepsis": "Sepsis (adult)", "covid": "ESKD COVID (187-sample roster)"}[coh], loc="left")
    axs[1].legend(fontsize=8, loc="lower right"); plt.tight_layout(); fig.savefig(f"{F}/rare_parent.png"); plt.close(fig)


claim_fig(); rare_fig()
