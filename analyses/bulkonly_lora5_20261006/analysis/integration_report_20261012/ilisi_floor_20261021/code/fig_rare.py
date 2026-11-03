import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from deck_common import *
SS = [("real B (author)", "#3E4C59"), ("real B (transferred)", "#AEB7C2"), ("prev K4000 (transferred)", "#D7191C"), ("floor (transferred)", "#2B83BA")]
fig, axs = plt.subplots(1, 2, figsize=(13, 6.4))
for ax, coh in zip(axs, ("sepsis", "covid")):
    D = pd.read_csv(f"{RB}/{coh}_rare_parent_share.csv"); D = D[D.group == "all"]
    subs = D[D.source == "real B (author)"].sort_values("pct_of_parent").subset.tolist(); par = D.drop_duplicates("subset").set_index("subset").parent
    y = np.arange(len(subs)); h = 0.2
    for j, (s, c) in enumerate(SS):
        d = D[D.source == s].set_index("subset").reindex(subs); v = np.maximum(d.pct_of_parent.to_numpy(), 0.003)
        ax.barh(y + (1.5 - j) * h, v, h, color=c, label=s)
    ax.set_xscale("log"); ax.set_xlim(0.003, 100); ax.set_yticks(y, [f"{s}  [{par[s]}]" for s in subs], fontsize=8.5)
    ax.set_xlabel("% of parent lineage (pooled; bar at 0.003 = none)"); ax.xaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
    ax.set_title(COH[coh], loc="left", fontsize=11, color="#1F2933")
h_, l_ = axs[0].get_legend_handles_labels(); fig.legend(h_[::-1], l_[::-1], loc="lower center", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "rare_parent_share"); print("ok")
