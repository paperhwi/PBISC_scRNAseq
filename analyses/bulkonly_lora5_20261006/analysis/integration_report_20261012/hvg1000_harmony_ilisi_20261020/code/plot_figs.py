"""Figures for run_hvg_harmony.py: joint UMAPs (HVG 1,000 and 2,000) and iLISI bar graph. usage: py -3.12 plot_figs.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES, FIG = OUT + "/results", OUT + "/figs"
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#333331", "xtick.color": "#555553", "ytick.color": "#555553"})
SRC = {"real": ("Real sepsis (Ye 2026)", "#2a78d6"), "k1000": ("PBISC K = 1,000", "#eb6834"), "k4000": ("PBISC K = 4,000", "#1baf7a")}
FOLD = {"B cells": "B / plasma", "Plasma cells": "B / plasma", "DC": "DC / macrophage", "Macrophages": "DC / macrophage"}
TYPES = [("CD4+ T cells", "#2a78d6"), ("CD8+ T cells", "#eb6834"), ("Monocytes", "#1baf7a"), ("B / plasma", "#eda100"),
         ("NK cells", "#e87ba4"), ("NKT cells", "#008300"), ("DC / macrophage", "#6250d6"), ("Other cells", "#e34948")]
U = {n: np.load(f"{RES}/joint_umap_hvg{n}.npz", allow_pickle=True) for n in (1000, 2000)}
R = pd.read_csv(f"{RES}/ilisi_hvg1000_vs_hvg2000.csv")


def clean(ax, title):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9, loc="left", color="#1a1a19")
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xlabel("UMAP1", fontsize=7, loc="left"); ax.set_ylabel("UMAP2", fontsize=7, loc="bottom")


# ---------- UMAP: rows = HVG setting; cols = all sources overlaid, real alone, K1k alone, K4k alone (grey = rest), cell type
fig, axes = plt.subplots(2, 5, figsize=(17, 7.2))
for r, n in enumerate((1000, 2000)):
    u, src = U[n]["umap"], U[n]["src"].astype(str); maj = pd.Series(U[n]["major"].astype(str)).replace(FOLD).to_numpy()
    post = R[(R.n_hvg == n) & (R.stage == "post")].set_index("comparison").iLISI
    order = np.random.default_rng(0).permutation(len(u)); ax = axes[r, 0]
    ax.scatter(u[order, 0], u[order, 1], s=0.3, c=[SRC[s][1] for s in src[order]], rasterized=True, linewidths=0)
    clean(ax, f"HVG {n:,} · all sources\niLISI (3 sources) = {post['real + K1k + K4k (3 sources)']:.2f} / 3")
    for c, s in enumerate(["real", "k1000", "k4000"], start=1):
        ax = axes[r, c]; m = src == s
        ax.scatter(u[~m, 0], u[~m, 1], s=0.2, c="#d9d8d3", rasterized=True, linewidths=0)
        ax.scatter(u[m, 0], u[m, 1], s=0.3, c=SRC[s][1], rasterized=True, linewidths=0)
        sub = "" if s == "real" else f"\niLISI real vs {'K1k' if s == 'k1000' else 'K4k'} = {post['real vs ' + ('K1k' if s == 'k1000' else 'K4k')]:.2f} / 2"
        clean(ax, f"HVG {n:,} · {SRC[s][0]}" + sub)
    ax = axes[r, 4]; col = dict(TYPES)
    ax.scatter(u[order, 0], u[order, 1], s=0.3, c=[col.get(x, "#999999") for x in maj[order]], rasterized=True, linewidths=0)
    clean(ax, f"HVG {n:,} · major cell type\n(real labels; generated = kNN transfer)")
h = [plt.Line2D([], [], marker="o", ls="", ms=6, color=c, label=l) for l, c in SRC.values()]
fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.36, -0.01), ncol=3, frameon=False, fontsize=9)
h2 = [plt.Line2D([], [], marker="o", ls="", ms=6, color=c, label=l) for l, c in TYPES]
fig.legend(handles=h2, loc="lower center", bbox_to_anchor=(0.86, -0.04), ncol=2, frameon=False, fontsize=8)
fig.suptitle("Sepsis: real + PBISC K = 1,000 + K = 4,000, Harmony (batch = source), HVG 1,000 vs HVG 2,000", x=0.01, ha="left", fontsize=11)
fig.tight_layout(rect=(0, 0.06, 1, 0.97)); [fig.savefig(f"{FIG}/umap_hvg1000_vs_hvg2000.{e}", dpi=300, bbox_inches="tight") for e in ("png", "pdf")]
plt.close(fig)

# ---------- iLISI bars: post-Harmony bars, pre-Harmony hollow marker
comps = ["real + K1k + K4k (3 sources)", "real vs K1k", "real vs K4k"]
HV = {1000: ("HVG 1,000", "#1c5cab"), 2000: ("HVG 2,000", "#86b6ef")}
fig, axs = plt.subplots(1, 2, figsize=(9, 3.8), gridspec_kw={"width_ratios": [1, 2]})
for ax, cs, top in [(axs[0], comps[:1], 3), (axs[1], comps[1:], 2)]:
    x = np.arange(len(cs)); w = 0.36
    for j, n in enumerate((1000, 2000)):
        d = R[R.n_hvg == n].set_index(["comparison", "stage"]).iLISI
        post = [d[(c, "post")] for c in cs]; pre = [d[(c, "pre")] for c in cs]; xx = x + (j - 0.5) * (w + 0.03)
        ax.bar(xx, post, w, color=HV[n][1], label=f"{HV[n][0]} · after Harmony", edgecolor="white", linewidth=1)
        ax.scatter(xx, pre, marker="o", s=26, facecolor="white", edgecolor="#333331", linewidth=1, zorder=3,
                   label="before Harmony (PCA projection)" if j == 0 else None)
        for xi, v in zip(xx, post): ax.text(xi, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=8, color="#1a1a19")
    ax.axhline(top, color="#8a8a85", lw=0.8, ls="--"); ax.text(x[-1] + 0.45, top, f"max {top}", fontsize=7, color="#555553", va="bottom", ha="right")
    ax.axhline(1, color="#8a8a85", lw=0.6); ax.set_ylim(0.9, top + 0.25); ax.set_xticks(x, cs); ax.set_xlim(-0.6, len(cs) - 0.4)
    ax.yaxis.grid(True, color="#eeeeec", lw=0.6); ax.set_axisbelow(True); ax.set_ylabel("median iLISI (source)")
axs[0].set_title("Three-way mixing", loc="left", fontsize=10); axs[1].set_title("Real vs each generated set", loc="left", fontsize=10)
axs[1].legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("iLISI: HVG 1,000 vs HVG 2,000 (same cells, same Harmony settings)", x=0.01, ha="left", fontsize=11)
fig.tight_layout(); [fig.savefig(f"{FIG}/ilisi_bar_hvg1000_vs_hvg2000.{e}", dpi=300, bbox_inches="tight") for e in ("png", "pdf")]
print("figs done")
