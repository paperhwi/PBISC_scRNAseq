import nbformat as nbf
nb = nbf.v4.new_notebook(); C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("# Summary figure for the RF realism test\nReads `results/rf_auc_whole.csv`, `rf_auc_within_lineage.csv`, `rf_auc_after_harmony.csv` (from `rf_realism.ipynb` and `rf_realism_harmony.ipynb`).")
code(r'''import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, logging
from IPython.display import Image, display
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
BASE = os.environ.get("PBISC_BASE", "F:/PBISC"); W = f"{BASE}/bulkonly_lora5_20261006/analysis/integration_report_20261012/rf_realism_20261017"
R, F = f"{W}/results", f"{W}/figs"
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 14, "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 200,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
U = pd.read_csv(f"{R}/rf_auc_whole.csv"); L = pd.read_csv(f"{R}/rf_auc_within_lineage.csv"); H = pd.read_csv(f"{R}/rf_auc_after_harmony.csv")
order = ["C0 real A vs real A", "C1 real A vs real B", "C2 real A vs real COVID", "C3 real A vs real B (PBISC read-out)",
         "C4 real A vs K = 1,000", "C5 real A vs K = 4,000", "C6 real B (PBISC read-out) vs K = 1,000", "C7 K = 1,000 vs K = 4,000"]
LAB = {"C0": "real A vs real A (same patients)", "C1": "real A vs real B (other patients)", "C2": "real A vs real COVID (other study)",
       "C3": "real A vs real B in PBISC read-out", "C4": "real A vs PBISC K = 1,000", "C5": "real A vs PBISC K = 4,000",
       "C6": "real B in PBISC read-out vs K = 1,000", "C7": "PBISC K = 1,000 vs K = 4,000"}
T = pd.DataFrame({"uncorrected, natural": U[U["mode"] == "natural"].groupby("comparison").auc.mean(),
                  "uncorrected, lineage-matched": U[U["mode"] == "matched"].groupby("comparison").auc.mean(),
                  "after Harmony, lineage-matched": H.groupby("comparison").auc_after_harmony.mean()}).reindex(order)
for l in ["Myeloid", "T", "NK", "B"]: T[f"within {l}"] = L[L["mode"] == f"within {l}"].groupby("comparison").auc.mean().reindex(order)
T.index = [LAB[c[:2]] for c in T.index]; T.round(3).to_csv(f"{R}/rf_auc_summary_table.csv"); display(T.round(3))
nl = L.groupby("mode").n_per_class.min(); print(nl.to_dict())''')
code(r'''fig, axs = plt.subplots(1, 2, figsize=(18, 7.2), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]; y = np.arange(len(T))[::-1]
for dy, k, col in [(0.27, "uncorrected, natural", "#C9CED8"), (0, "uncorrected, lineage-matched", "#9AA3B2"), (-0.27, "after Harmony, lineage-matched", "#2A78D6")]:
    ax.barh(y + dy, T[k] - 0.5, left=0.5, height=0.26, color=col, label=k)
    for yi, v in zip(y + dy, T[k]): ax.text(max(v, 0.5) + 0.008, yi, f"{v:.2f}", va="center", fontsize=11, color="#1F2430")
ax.set_yticks(y, T.index, fontsize=13); ax.set_xlim(0.38, 1.07); ax.axvline(0.5, color="#1F2430", lw=1)
ax.axvline(0.65, color="#C4373A", lw=1.2, ls="--"); ax.text(0.656, len(T) - 0.4, "scGAN vs real (Marouf 2020): 0.65", color="#C4373A", fontsize=12)
ax.set_xlabel("random forest AUC (1,000 trees, 5-fold CV, PCA 50)"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=12)
ax.set_title("Whole samples, 3,000 cells per class", loc="left", fontsize=15)
ax = axs[1]; lin = ["Myeloid", "T", "NK", "B"]; pick = list(T.index[[1, 3, 4, 5, 6]]); cols = ["#4A4F5C", "#9AA3B2", "#2A78D6", "#EB6834", "#1BAF7A"]; w = 0.16
for j, (c, col) in enumerate(zip(pick, cols)):
    ax.bar(np.arange(4) + (j - 2) * w, T.loc[c, [f"within {l}" for l in lin]].to_numpy() - 0.5, bottom=0.5, width=w, color=col, label=c)
ax.set_xticks(range(4), [f"{l}\n{nl[f'within {l}']:,}/class" for l in lin], fontsize=12); ax.set_ylim(0.45, 1.02); ax.axhline(0.5, color="#1F2430", lw=1)
ax.set_ylabel("random forest AUC, uncorrected"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=11.5)
ax.set_title("Within lineage", loc="left", fontsize=15)
plt.tight_layout(); fig.savefig(f"{F}/rf_auc_summary.png", bbox_inches="tight"); fig.savefig(f"{F}/rf_auc_summary.pdf", bbox_inches="tight"); plt.close(fig)
display(Image(f"{F}/rf_auc_summary.png", width=1100))''')
nb["cells"] = C; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/rf_realism_20261017/rf_summary_figure.ipynb")
