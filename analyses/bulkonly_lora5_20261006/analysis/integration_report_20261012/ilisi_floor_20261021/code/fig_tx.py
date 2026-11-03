import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from deck_common import *
# A: pseudobulk correlation per six-type (three gene sets) ; B: marker dot plot (sepsis + covid)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
TY = ["T cell", "NK cell", "B cell", "Monocyte", "DC"]; SETS = ["all genes", "2,000 ref HVG", "genes ranked 2,001+"]
for ax, coh in zip(axs, ["sepsis", "covid"]):
    P = pd.read_csv(f"{RB}/{coh}_tx_pb.csv"); x = np.arange(len(SETS)); w = 0.25
    for j, s in enumerate(["realT1f", "floor", "prev4k"]):
        d = P[P.source == s]; m = d.groupby("gene_set").pearson.median().reindex(SETS); lo = d.groupby("gene_set").pearson.min().reindex(SETS); hi = d.groupby("gene_set").pearson.max().reindex(SETS)
        xx = x + (j - 1) * (w + 0.03); ax.bar(xx, m, w, color=SRC[s][1], label=SRC[s][0]); ax.errorbar(xx, m, yerr=[m - lo, hi - m], fmt="none", ecolor="#1F2933", lw=0.9, capsize=2.5)
        for xi, v, h in zip(xx, m, hi): ax.text(xi, h + 0.008, f"{v:.2f}", ha="center", va="bottom", fontsize=8.5, color="#1F2933", fontweight="bold")
    ax.set_xticks(x, SETS); ax.set_ylim(0.3, 1.02); ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
    ax.set_title(COH[coh] + " · pseudobulk r with real half A (median of 5 types; bar = range)", loc="left", fontsize=10, color="#1F2933")
axs[0].set_ylabel("Pearson r, log1p mean CP10k per gene"); axs[0].set_ylim(0.3, 1.06)
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "pseudobulk_r")
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
for ax, coh in zip(axs, ["sepsis", "covid"]):
    D = pd.read_csv(f"{RB}/{coh}_tx_dot.csv"); genes = list(dict.fromkeys(D.gene)); genes = [g for g in genes if g not in ("HBB", "MKI67")]
    rows = [(t, s) for t in TY for s in ["ref", "floor", "prev4k"]]
    for i, (t, s) in enumerate(rows):
        d = D[(D.six_type == t) & (D.source == s)].set_index("gene").reindex(genes)
        ax.scatter(np.arange(len(genes)), np.full(len(genes), -i), s=d.frac.fillna(0) * 90, c=d.mean_log.fillna(0), cmap="Blues", vmin=0, vmax=3,
                   edgecolors=SRC[s][1], linewidths=1.0)
    ax.set_yticks([-i for i in range(len(rows))], [f"{t} · {'real A' if s == 'ref' else ('floor' if s == 'floor' else 'prev')}" for t, s in rows], fontsize=7.5)
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=8); ax.set_title(COH[coh] + " · marker detection (size) and mean log CP10k (colour)", loc="left", fontsize=10, color="#1F2933")
    for k in range(1, len(TY)): ax.axhline(-3 * k + 0.5, color="#E4E7EB", lw=0.8)
    for k in range(4, len(genes), 4): ax.axvline(k - 0.5, color="#E4E7EB", lw=0.8)
for f in (0.25, 0.5, 1.0): axs[1].scatter([], [], s=f * 90, c="#9AA5B1", label=f"{int(f*100)}% detected")
for s in ["ref", "floor", "prev4k"]: axs[1].scatter([], [], s=50, facecolor="white", edgecolors=SRC[s][1], linewidths=1.4, label=SRC[s][0].split(" (")[0])
axs[1].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
fig.tight_layout(); save(fig, "marker_dot")
print("ok")
