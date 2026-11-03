import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from deck_common import *
SIXC = dict(SIX); J = {(c, s): json.load(open(f"{RB}/{c}_{s}_summary.json")) for c in ("sepsis", "covid") for s in ("realU", "floor", "prev4k")}
def sc_(ax, u, c, s=0.15): o = np.random.default_rng(0).permutation(len(u)); ax.scatter(u[o, 0], u[o, 1], s=s, c=np.asarray(c, dtype=object)[o], rasterized=True, linewidths=0)
for coh in ("sepsis", "covid"):
    U = {s: np.load(f"{RB}/{coh}_{s}_umap.npz", allow_pickle=True) for s in ("realU", "floor", "prev4k")}
    fig, axs = plt.subplots(2, 4, figsize=(14, 7.2))
    for k, (s, tag) in enumerate([("prev4k", "pre"), ("floor", "pre"), ("prev4k", "post"), ("floor", "post")]):
        u = U[s][f"umap_{tag}"]; src = U[s]["src"]; ax = axs[0, k]
        sc_(ax, u, np.where(src == 0, "#B8BEC7", SRC[s][1]), s=0.2)
        clean(ax, f"{'Before' if tag == 'pre' else 'After'} Harmony · {SRC[s][0]}\niLISI {J[(coh, s)]['iLISI_' + tag]:.2f}")
    for k, s in enumerate(["realU", "prev4k", "floor"]):
        u = U[s]["umap_post"]; sc_(axs[1, k], u, [SIXC.get(x, "#9AA5B1") for x in U[s]["six"]]); clean(axs[1, k], f"After Harmony · real A + {dict(realU='real B', prev4k='previous K = 4,000', floor='floor read-out')[s]}\ncell type (real: author; generated: own label)")
    u = U["floor"]["umap_post"]; lc = U["floor"]["leiden"]; pal = plt.get_cmap("tab20")(np.linspace(0, 1, 20)); sc_(axs[1, 3], u, [pal[i % 20] for i in lc])
    clean(axs[1, 3], f"After Harmony · real A + floor\njoint Leiden clusters (n = {J[(coh, 'floor')]['n_clusters']})")
    h = [plt.Line2D([], [], marker="o", ls="", ms=6, color="#C9CED6", label="Real half A (reference)")] + [plt.Line2D([], [], marker="o", ls="", ms=6, color=SRC[s][1], label=SRC[s][0]) for s in ("prev4k", "floor")] + \
        [plt.Line2D([], [], marker="o", ls="", ms=6, color=c, label=l) for l, c in SIX]
    fig.legend(handles=h, loc="lower center", ncol=9, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02)); fig.suptitle(COH[coh] + " · 25,000 real A + 25,000 source cells per panel", x=0.01, ha="left", fontsize=11, color="#1F2933")
    fig.tight_layout(rect=(0, 0.04, 1, 0.97)); save(fig, f"umap_{coh}")
# cluster matching: cluster frequency real A vs source
fig, axs = plt.subplots(2, 3, figsize=(12, 7.4))
for r, coh in enumerate(("sepsis", "covid")):
    for k, s in enumerate(("realU", "prev4k", "floor")):
        C = pd.read_csv(f"{RB}/{coh}_{s}_leiden_clusters.csv"); ax = axs[r, k]; e = 3e-4
        w = C.src_share_norm; col = np.where(w > 0.8, SRC[s][1], np.where(w < 0.2, "#7F7F7F", "#1F2933"))
        ax.scatter(C.frac_ref + e, C.frac_src + e, s=18 + 500 * (C.frac_ref + C.frac_src), c=col, alpha=0.75, edgecolors="white", linewidths=0.5)
        ax.plot([e, 0.5], [e, 0.5], color="#9AA5B1", lw=0.8, ls=":"); ax.plot([e, 0.5], [4 * e, 2.0], color="#E4E7EB", lw=0.6); ax.plot([4 * e, 2.0], [e, 0.5], color="#E4E7EB", lw=0.6)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(e * 0.8, 0.5); ax.set_ylim(e * 0.8, 0.5)
        for _, row in C[(C.src_share_norm > 0.8) | (C.src_share_norm < 0.2)].iterrows():
            ax.annotate(str(row.ref_major_fine)[:22], (row.frac_ref + e, row.frac_src + e), fontsize=7, color="#3E4C59", xytext=(3, 2), textcoords="offset points")
        j = J[(coh, s)]
        ax.set_title(f"{COH[coh].split(' (')[0]} · real A vs {dict(realU='real B', prev4k='previous K = 4,000', floor='floor read-out')[s]}\nr = {j['cluster_freq_pearson']:.2f}, JSD = {j['cluster_freq_JSD']:.3f}, "
                     f"{100 * j['src_cells_in_mixed_clusters']:.0f}% in shared clusters", loc="left", fontsize=9.5, color="#1F2933")
        if k == 0: ax.set_ylabel("share of source cells in cluster")
        if r == 1: ax.set_xlabel("share of real half A cells in cluster")
h = [plt.Line2D([], [], marker="o", ls="", ms=7, color="#1F2933", label="shared cluster (size-normalised source share 0.2–0.8)"), plt.Line2D([], [], marker="o", ls="", ms=7, color="#7F7F7F", label="real-dominated (< 0.2)"),
     plt.Line2D([], [], marker="o", ls="", ms=7, color="#2B83BA", label="source-dominated (> 0.8; colour = source)")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02)); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "cluster_matching")
print("ok")
