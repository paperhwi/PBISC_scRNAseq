import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from deck_common import *
U = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/results"
R = pd.read_csv(f"{U}/ilisi_all.csv"); R = R[R.variant == "V0"]; L = pd.read_csv(f"{U}/ilisi_post_by_lineage_V0.csv")
ORD = [("realU", "Real B\n(other pts)", "#3E4C59"), ("realT1", "Real B\n(depth-m.)", "#AEB7C2"), ("prev4k", "Previous\nK = 4,000", "#D7191C"), ("floor", "Floor\nread-out", "#2B83BA")]
fig, axs = plt.subplots(1, 3, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1, 1, 1.6]})
for ax, coh in zip(axs[:2], ("sepsis", "covid")):
    d = R[R.cohort == coh].set_index("source")
    for i, (s, lab, c) in enumerate(ORD):
        ax.bar(i, d.loc[s, "iLISI_post"], 0.62, color=c); ax.scatter(i, d.loc[s, "iLISI_pre"], s=36, facecolor="white", edgecolor="#1F2933", lw=1.1, zorder=3)
        ax.text(i, d.loc[s, "iLISI_post"] + 0.02, f"{d.loc[s, 'iLISI_post']:.2f}", ha="center", va="bottom", fontsize=10, fontweight="bold", color="#1F2933")
    ax.set_xticks(range(4), [o[1] for o in ORD], fontsize=8.5); ax.set_ylim(0.9, 2.1); ax.axhline(2, color="#9AA5B1", lw=0.8, ls=":")
    ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True); ax.set_title(COH[coh], loc="left", fontsize=11, color="#1F2933")
axs[0].set_ylabel("median iLISI vs real half A (1 = separate, 2 = mixed)"); axs[0].scatter([], [], s=36, facecolor="white", edgecolor="#1F2933", label="before Harmony"); axs[0].legend(frameon=False, fontsize=9, loc="upper right")
ax = axs[2]; LIN = ["T", "NK", "B", "Myeloid"]; k = 0; xt, xl = [], []
for coh, t in (("sepsis", "Sepsis"), ("covid", "COVID")):
    for ln in LIN:
        for j, (s, _, c) in enumerate([o for o in ORD if o[0] != "realT1"]):
            v = L[(L.cohort == coh) & (L.source == s)][ln].iloc[0]; ax.bar(k + j * 0.26, v, 0.24, color=c)
        xt.append(k + 0.26); xl.append(f"{ln}\n{t}"); k += 1
    k += 0.5
ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylim(1, 2); ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
ax.set_title("After Harmony, within lineage", loc="left", fontsize=11, color="#1F2933"); ax.set_ylabel("median iLISI within lineage")
h = [plt.Rectangle((0, 0), 1, 1, color=c, label={"realU": "Real B (other patients)", "prev4k": "Previous K = 4,000", "floor": "Floor read-out"}[s]) for s, l, c in ORD if s != "realT1"]; ax.legend(handles=h, frameon=False, fontsize=8.5, ncol=3, loc="upper center")
fig.tight_layout(); save(fig, "ilisi_bar_deck"); print("ok")
