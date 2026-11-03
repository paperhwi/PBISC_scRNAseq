import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "axes.edgecolor": "#3E4C59", "axes.labelcolor": "#1F2933", "xtick.color": "#3E4C59", "ytick.color": "#3E4C59"})
RB = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"; FIG = "/mnt/user-data/outputs/ilisi_floor_20261021/figs/deck"
SRC = {"ref": ("Real, half A (reference)", "#7F7F7F"), "realU": ("Real, half B (other patients)", "#3E4C59"),
       "realT1f": ("Real half B, depth-matched to floor", "#AEB7C2"), "realT1p": ("Real half B, depth-matched to prev K4000", "#D9DEE4"),
       "floor": ("Floor read-out (new)", "#2B83BA"), "prev4k": ("Previous K = 4,000", "#D7191C")}
COH = {"sepsis": "Sepsis (Ye 2026, adults)", "covid": "COVID (Imperial ESKD)"}
SIX = [("T cell", "#2a78d6"), ("NK cell", "#e87ba4"), ("B cell", "#eda100"), ("Monocyte", "#1baf7a"), ("DC", "#6250d6"), ("Other", "#9AA5B1")]
def save(fig, name):
    fig.savefig(f"{FIG}/{name}.png", dpi=250, bbox_inches="tight", facecolor="white"); fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight"); plt.close(fig)
def clean(ax, title=""):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10, color="#1F2933")
    for s in ax.spines.values(): s.set_visible(False)
