"""Slide figures for the Korean HVG 1,000 deck, in the integration-report deck palette
(real = grey 7F7F7F, K = 1,000 = 2B83BA, K = 4,000 = D7191C). usage: py -3.12 plot_deck_figs.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES, FIG = OUT + "/results", OUT + "/figs/deck"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#3E4C59", "axes.labelcolor": "#1F2933", "xtick.color": "#3E4C59", "ytick.color": "#3E4C59"})
SRC = {"real": ("Real sepsis (Ye 2026)", "#7F7F7F"), "k1000": ("Generated K = 1,000", "#2B83BA"), "k4000": ("Generated K = 4,000", "#D7191C")}
FOLD = {"B cells": "B / plasma", "Plasma cells": "B / plasma", "DC": "DC / macrophage", "Macrophages": "DC / macrophage"}
TYPES = [("CD4+ T cells", "#2a78d6"), ("CD8+ T cells", "#eb6834"), ("Monocytes", "#1baf7a"), ("B / plasma", "#eda100"),
         ("NK cells", "#e87ba4"), ("NKT cells", "#008300"), ("DC / macrophage", "#6250d6"), ("Other cells", "#e34948")]
U = {n: np.load(f"{RES}/joint_umap_hvg{n}.npz", allow_pickle=True) for n in (1000, 2000)}
R = pd.read_csv(f"{RES}/ilisi_hvg1000_vs_hvg2000.csv")
post = {n: R[(R.n_hvg == n) & (R.stage == "post")].set_index("comparison").iLISI for n in (1000, 2000)}


def clean(ax, title):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10, color="#1F2933")
    for s in ax.spines.values(): s.set_visible(False)


def save(fig, name):
    fig.savefig(f"{FIG}/{name}.png", dpi=250, bbox_inches="tight", facecolor="white"); plt.close(fig)


# A: HVG 1,000 joint UMAP, all sources + each source over grey
u, src = U[1000]["umap"], U[1000]["src"].astype(str); order = np.random.default_rng(0).permutation(len(u))
fig, axes = plt.subplots(1, 4, figsize=(14, 3.9))
axes[0].scatter(u[order, 0], u[order, 1], s=0.25, c=[SRC[s][1] for s in src[order]], rasterized=True, linewidths=0)
clean(axes[0], "All three sources")
for ax, s in zip(axes[1:], ["real", "k1000", "k4000"]):
    m = src == s
    ax.scatter(u[~m, 0], u[~m, 1], s=0.2, c="#E4E7EB", rasterized=True, linewidths=0)
    ax.scatter(u[m, 0], u[m, 1], s=0.25, c=SRC[s][1], rasterized=True, linewidths=0)
    clean(ax, SRC[s][0] + ("" if s == "real" else f"\niLISI vs real = {post[1000]['real vs ' + ('K1k' if s == 'k1000' else 'K4k')]:.2f}"))
h = [plt.Line2D([], [], marker="o", ls="", ms=7, color=c, label=l) for l, c in SRC.values()]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.06))
fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "umap_hvg1000_sources")

# B: HVG 1,000 vs 2,000, all sources and cell type
fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.4))
for r, n in enumerate((1000, 2000)):
    u, src = U[n]["umap"], U[n]["src"].astype(str); maj = pd.Series(U[n]["major"].astype(str)).replace(FOLD).to_numpy()
    order = np.random.default_rng(0).permutation(len(u)); col = dict(TYPES)
    axes[r, 0].scatter(u[order, 0], u[order, 1], s=0.2, c=[SRC[s][1] for s in src[order]], rasterized=True, linewidths=0)
    clean(axes[r, 0], f"HVG {n:,} · source")
    axes[r, 1].scatter(u[order, 0], u[order, 1], s=0.2, c=[col.get(x, "#999999") for x in maj[order]], rasterized=True, linewidths=0)
    clean(axes[r, 1], f"HVG {n:,} · major cell type")
h2 = [plt.Line2D([], [], marker="o", ls="", ms=6, color=c, label=l) for l, c in TYPES]
fig.legend(handles=h + h2, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
fig.tight_layout(); save(fig, "umap_hvg1000_vs_2000")

# C: iLISI bars. Pairwise panel: bar colour = generated set; dark = HVG 1,000, light = HVG 2,000; open circle = before Harmony.
REF = [("Real held-out samples (same study)", 1.74, "#5E4FA2"), ("Real, other study (COVID)", 1.51, "#4DAF4A"), ("Real, PBISC-like read-out", 1.44, "#E69F00")]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={"width_ratios": [1, 2.1]})
light = {"#2B83BA": "#9CC6E2", "#D7191C": "#F09A9B", "#3E4C59": "#AEB7C2"}
for ax, cs, top in [(axs[0], [("real + K1k + K4k (3 sources)", "#3E4C59", "3 sources")], 3),
                    (axs[1], [("real vs K1k", "#2B83BA", "Real vs K = 1,000"), ("real vs K4k", "#D7191C", "Real vs K = 4,000")], 2)]:
    x = np.arange(len(cs)); w = 0.34
    for i, (c, colr, lab) in enumerate(cs):
        for j, n in enumerate((1000, 2000)):
            d = R[R.n_hvg == n].set_index(["comparison", "stage"]).iLISI; xx = i + (j - 0.5) * (w + 0.04)
            ax.bar(xx, d[(c, "post")], w, color=colr if n == 1000 else light[colr])
            ax.scatter(xx, d[(c, "pre")], s=34, facecolor="white", edgecolor="#1F2933", lw=1.1, zorder=3)
            ax.text(xx, d[(c, "post")] + 0.025, f"{d[(c, 'post')]:.2f}", ha="center", va="bottom", fontsize=10, color="#1F2933", fontweight="bold")
            ax.text(xx, 0.93, f"HVG\n{n:,}", ha="center", va="bottom", fontsize=8.5, color="white" if n == 1000 else "#1F2933")
    ax.set_xticks(x, [c[2] for c in cs]); ax.set_xlim(-0.6, len(cs) - 0.4); ax.set_ylim(0.9, top + 0.12)
    ax.axhline(top, color="#9AA5B1", lw=0.8, ls=":"); ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
    ax.set_ylabel(f"median iLISI (1 = separate, {top} = fully mixed)")
for lab, v, c in REF:
    axs[1].axhline(v, color=c, lw=1.2, ls="--"); axs[1].text(1.62, v + 0.012, f"{lab} {v:.2f}", color=c, fontsize=8.5, ha="right", va="bottom")
axs[0].set_title("Three-way mixing (max 3)", loc="left", fontsize=11, color="#1F2933")
axs[1].set_title("Real vs each generated set (max 2); dashed = real-vs-real reference, pooled, HVG 2,000 (results/check)", loc="left", fontsize=11, color="#1F2933")
axs[1].scatter([], [], s=34, facecolor="white", edgecolor="#1F2933", label="before Harmony (PCA projection)")
axs[1].legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout(); save(fig, "ilisi_bar_deck")
print("deck figs done")
