import json, logging, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 220})
O = "/home/claude/work/report/check"
V = ["V0_refHVG_refPCA", "V1_jointHVG_refPCA", "V2_refHVG_jointPCA", "V3_jointHVG_jointPCA"]
VL = ["V0\nref HVG\nref PCA", "V1\njoint HVG\nref PCA", "V2\nref HVG\njoint PCA", "V3\njoint HVG\njoint PCA"]
rows = []
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
for ax, coh in zip(axs, ["sepsis", "covid"]):
    r = json.load(open(f"{O}/{coh}_hvgpca.json")); oth = "covid" if coh == "sepsis" else "sepsis"
    S = [("realU", "Real, held-out samples", "#6a51a3", "-"), ("realT1", "Real, depth-matched", "#9e9ac8", "-"),
         ("realT2", "Real, PBISC-like read-out", "#e69f00", "-"), (oth, "Real, other study", "#5aa469", "-"),
         ("k1000", "Generated K = 1,000", "#2b83ba", "-"), ("k4000", "Generated K = 4,000", "#d7191c", "--")]
    for k, lab, c, ls in S:
        y = [r[v][k]["iLISI_post"] for v in V]; ax.plot(range(4), y, ls, marker="o", color=c, lw=1.8, ms=5, label=lab)
        if k not in ("realT1", "k4000"): ax.text(3.08, y[-1], f"{y[-1]:.2f}", va="center", fontsize=7.5, color=c)
        for v in V: rows.append({"cohort": coh, "variant": v, "source": k, **{m: r[v][k].get(m) for m in ["iLISI_pre", "iLISI_post", "own_share_post", "iLISI_post_Myeloid", "conf_fine_median", "acc_fine", "hvg_overlap_with_ref_hvg"]}})
    ax.set_xticks(range(4), VL, fontsize=8); ax.set_xlim(-0.2, 3.4); ax.set_ylim(0.95, 2.0); ax.axhline(2, c="k", lw=0.6, ls=":")
    ax.set_title({"sepsis": "Sepsis reference", "covid": "ESKD COVID reference (187-sample roster)"}[coh], loc="left")
    ax.grid(axis="y", color="#eeeeee", lw=0.6)
axs[0].set_ylabel("median iLISI vs reference after Harmony\n(1 = separate, 2 = fully mixed)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=6, fontsize=8, bbox_to_anchor=(0.5, 0))
plt.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig("/home/claude/work/report/figs/hvgpca_sensitivity.png"); plt.close(fig)
pd.DataFrame(rows).to_csv(f"{O}/hvgpca_sensitivity_summary.csv", index=False)
