import logging, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 220})
O = "/home/claude/work/report/rescore"
d = pd.read_csv(f"{O}/summary_panel_matched.csv"); d = d[d.variant == "A_panel"].set_index("claim_id")
names = {"SEP-1": "SEP-1 NR4A2+ CD4 T", "SEP-2": "SEP-2 GNLY/CCL4 CD8 T", "SEP-3": "SEP-3 IFNG/CCL4 NKT*", "SEP-3-NK": "SEP-3-NK PRF1/CCL4 NK",
         "SEP-4A": "SEP-4A HLA-DR CD14 mono", "SEP-4B": "SEP-4B RETN CD14 mono", "SEP-5": "SEP-5 MKI67 CD14 mono (adult)"}
order = ["SEP-1", "SEP-2", "SEP-3", "SEP-3-NK", "SEP-4A", "SEP-4B", "SEP-5"][::-1]
SRC = [("real_classifier", "real_score", "Real (same classifier)", "#555555"), ("k1k_classifier", "k1k_score", "K = 1,000", "#2b83ba"), ("k4k_classifier", "k4k_score", "K = 4,000", "#d7191c")]
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for ax, kind, xl in [(axs[0], 0, "Δ classifier-positive fraction in parent (adult sepsis − HC)"), (axs[1], 1, "Δ mean marker score in parent (adult sepsis − HC)")]:
    for i, c in enumerate(order):
        usable = d.loc[c, "cv_auc"] >= 0.65 and not (c == "SEP-3")
        if not usable and kind == 0: ax.axhspan(i - 0.45, i + 0.45, color="#f2f2f2", zorder=0)
        for j, s in enumerate(SRC):
            key = s[kind]; y = i + (1 - j) * 0.22; m, lo, hi = d.loc[c, f"{key}_diff"], d.loc[c, f"{key}_lo"], d.loc[c, f"{key}_hi"]
            ax.errorbar(m, y, xerr=[[m - lo], [hi - m]], fmt="o", color=s[3], ms=4, lw=1.2, capsize=2, label=s[2] if i == 0 else None)
    ax.axvline(0, c="k", lw=0.6); ax.set_xlabel(xl); ax.grid(axis="x", color="#eeeeee", lw=0.6)
axs[0].set_yticks(range(len(order)), [f"{names[c]}  (AUC {d.loc[c, 'cv_auc']:.2f})" for c in order], fontsize=8.5)
axs[0].set_title("Thresholded (panel-matched markers, refitted thresholds)", loc="left", fontsize=10)
axs[1].set_title("Threshold-free score", loc="left", fontsize=10)
axs[1].legend(fontsize=8, loc="upper left")
fig.text(0.01, 0.01, "Grey rows: classifier unusable (AUC < 0.65) or parent not mappable (*generated NKT parent = DN T proxy). Bars: 95% bootstrap CI over samples; generated = mean of 5 draws per sample.", fontsize=7.5, color="#555555")
plt.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(f"{O}/rescore_contrasts.png"); fig.savefig("/home/claude/work/report/figs/rescore_contrasts.png"); plt.close(fig)
