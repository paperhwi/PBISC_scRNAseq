"""Render consistent, distinguishable draw-color UMAPs for both K settings."""
from pathlib import Path
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("MPLBACKEND", "Agg")

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

matplotlib.use("Agg")
ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = Path(__file__).resolve().parent
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]

for k, folder in [(1000, ROOT / "bulkonly_lora5_20261006/analysis/harmony/sepsis_k1000"),
                  (4000, ANALYSIS / "sepsis_k4000")]:
    path = folder / f"sepsis_k{k}_harmony_embedding.h5ad"
    a = ad.read_h5ad(path, backed="r")
    xy = np.asarray(a.obsm["X_umap"])
    draws = a.obs["draw"].astype(str).to_numpy()
    fig, ax = plt.subplots(figsize=(7.0, 5.1), constrained_layout=True)
    ax.set_facecolor("white")
    for i in range(5):
        draw = f"draw{i}"
        m = draws == draw
        ax.scatter(xy[m, 0], xy[m, 1], s=1.0, alpha=0.48, c=PALETTE[i], linewidths=0, rasterized=True)
    ax.set_title(f"Ye 2026 sepsis K={k:,} · Harmony embedding", loc="left", fontsize=11, weight="semibold", color="#171717", pad=10)
    ax.set_xlabel("UMAP 1", fontsize=9, color="#404040")
    ax.set_ylabel("UMAP 2", fontsize=9, color="#404040")
    ax.tick_params(labelsize=8, colors="#5E6462", length=2.5, width=0.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#B6B6B0")
    ax.spines["bottom"].set_color("#B6B6B0")
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                      markerfacecolor=PALETTE[i], markeredgewidth=0, label=f"draw{i}") for i in range(5)]
    ax.legend(handles=handles, title="Model draw", frameon=False, fontsize=8, title_fontsize=8,
              ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.savefig(folder / "umap_by_draw.png", dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(folder / "umap_by_draw.svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)
    a.file.close()
