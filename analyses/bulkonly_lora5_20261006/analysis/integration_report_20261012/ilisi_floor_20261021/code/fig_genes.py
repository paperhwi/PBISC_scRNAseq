import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from deck_common import *
fig, axs = plt.subplots(2, 2, figsize=(12, 7.2))
for r, coh in enumerate(["sepsis", "covid"]):
    G = pd.read_csv(f"{RB}/{coh}_tx_gene.csv.gz"); ref = G[G.source == "ref"].set_index("gene")
    for s in ["ref", "realT1f", "prev4k", "floor"]:
        d = G[G.source == s].sort_values("real_rank"); y = d.mean_cp10k.rolling(51, center=True, min_periods=10).median()
        axs[r, 0].plot(d.real_rank, np.maximum(y, 1e-5), color=SRC[s][1], lw=1.6 if s in ("floor", "prev4k") else 1.3, ls="--" if s == "realT1f" else "-", label=SRC[s][0])
        if s == "ref": continue
        rr = (d.set_index("gene").mean_cp10k / ref.mean_cp10k.reindex(d.gene).to_numpy()); rr = pd.Series(rr.to_numpy(), index=d.real_rank.to_numpy()).replace(np.inf, np.nan)
        axs[r, 1].plot(rr.index, np.clip(rr.rolling(151, center=True, min_periods=30).median(), 1e-3, None), color=SRC[s][1], lw=1.6, ls="--" if s == "realT1f" else "-", label=SRC[s][0])
    for ax in axs[r]:
        ax.set_yscale("log"); ax.axvspan(2000, G.real_rank.max(), color="#F5F7FA", zorder=0); ax.set_xlim(0, G.real_rank.max())
    axs[r, 0].set_ylim(1e-4, 300); axs[r, 0].set_ylabel("mean CP10k per gene\n(running median, 51 genes)")
    axs[r, 1].set_ylim(1e-3, 5); axs[r, 1].axhline(1, color="#9AA5B1", lw=0.8, ls=":"); axs[r, 1].set_ylabel("source / real half A\n(running median, 151 genes)")
    axs[r, 0].set_title(COH[coh] + " · gene count curve", loc="left", fontsize=11, color="#1F2933"); axs[r, 1].set_title(COH[coh] + " · ratio to real", loc="left", fontsize=11, color="#1F2933")
    axs[r, 1].text(2050, 3.2, "rank > 2,000 (low-expression genes)", fontsize=8.5, color="#7B8794")
axs[1, 0].set_xlabel("genes ordered by real mean (left = highest)"); axs[1, 1].set_xlabel("genes ordered by real mean (left = highest)")
h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "gene_count_curve")
# per-cell detected genes and mean-detection relationship
fig, axs = plt.subplots(2, 3, figsize=(12, 6.6), gridspec_kw={"width_ratios": [1.15, 1, 1]})
for r, coh in enumerate(["sepsis", "covid"]):
    C = pd.read_csv(f"{RB}/{coh}_tx_cell.csv.gz"); G = pd.read_csv(f"{RB}/{coh}_tx_gene.csv.gz")
    order = ["ref", "realT1f", "floor", "realT1p", "prev4k"]; ax = axs[r, 0]
    data = [C[C.source == s].n_genes.to_numpy() for s in order]
    vp = ax.violinplot(data, showextrema=False, widths=0.8)
    for b, s in zip(vp["bodies"], order): b.set_facecolor(SRC[s][1]); b.set_alpha(0.85); b.set_edgecolor("none")
    for i, d in enumerate(data): ax.text(i + 1, np.percentile(d, 99) * 1.02, f"{np.median(d):.0f}", ha="center", fontsize=9, fontweight="bold", color="#1F2933")
    ax.set_xticks(range(1, 6), ["real A", "real B\n(floor depth)", "floor", "real B\n(prev depth)", "prev\nK4000"], fontsize=8.5); ax.set_ylim(0, None)
    ax.set_ylabel("detected genes per cell"); ax.set_title(COH[coh] + " · genes per cell (median)", loc="left", fontsize=10.5, color="#1F2933")
    for s in ["realT1f", "floor", "realT1p", "prev4k"]:
        d = G[G.source == s]; d = d[d.mean_count > 0]; b = pd.cut(np.log10(d.mean_count), np.linspace(-4, 1.5, 34)); m = d.groupby(b, observed=True).agg(x=("mean_count", "median"), y=("det", "median"))
        axs[r, 1].plot(m.x, m.y, color=SRC[s][1], lw=1.6, ls="--" if s.startswith("real") else "-", label=SRC[s][0])
        cv = d.var_count / d.mean_count ** 2; m2 = d.assign(cv=cv).groupby(b, observed=True).agg(x=("mean_count", "median"), y=("cv", "median"))
        axs[r, 2].plot(m2.x, m2.y, color=SRC[s][1], lw=1.6, ls="--" if s.startswith("real") else "-")
    xx = np.logspace(-4, 1.5, 100); axs[r, 1].plot(xx, 1 - np.exp(-xx), color="#1F2933", lw=0.8, ls=":", label="Poisson 1 − e^(−μ)")
    axs[r, 2].plot(xx, 1 / xx, color="#1F2933", lw=0.8, ls=":")
    for ax in axs[r, 1:]: ax.set_xscale("log")
    axs[r, 1].set_ylabel("fraction of cells detecting the gene"); axs[r, 2].set_yscale("log"); axs[r, 2].set_ylabel("CV² of counts")
    axs[r, 1].set_title("mean–detection", loc="left", fontsize=10.5, color="#1F2933"); axs[r, 2].set_title("mean–variance (dotted = Poisson)", loc="left", fontsize=10.5, color="#1F2933")
axs[1, 1].set_xlabel("mean count per cell"); axs[1, 2].set_xlabel("mean count per cell")
h, l = axs[0, 1].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "cell_shape")
print("ok")
