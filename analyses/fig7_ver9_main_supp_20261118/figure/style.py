"""figures4papers house style (.scientific-figure-making skill): apply_publication_style, PALETTE, finalize_figure."""
import logging
from dataclasses import dataclass
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_1": "#DDF3DE",
    "green_2": "#AADCA9",
    "green_3": "#8BCF8B",
    "red_1": "#F6CFCB",
    "red_2": "#E9A6A1",
    "red_strong": "#B64342",
    "neutral": "#CFCECE",
    "highlight": "#FFD700",
    "yellow_soft": "#FDEBA8",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "gray_1": "#767676",
    "gray_2": "#4D4D4D",
    "gray_3": "#272727",
    "magenta": "#EA84DD",
}

SRC_COL = {"Real": PALETTE["gray_3"], "PBISC": PALETTE["blue_main"], "scSemiProfiler": PALETTE["red_strong"]}
SRC_MK = {"Real": "D", "PBISC": "o", "scSemiProfiler": "^"}
CT_COL = {
    "T cell": PALETTE["blue_secondary"],
    "NK cell": PALETTE["green_3"],
    "B cell": PALETTE["violet"],
    "Monocyte": PALETTE["teal"],
    "DC": PALETTE["red_2"],
    "Other": PALETTE["neutral"],
}


@dataclass(frozen=True)
class FigureStyle:
    font_size: int = 14
    axes_linewidth: float = 1.8
    font_family: tuple = ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans", "sans-serif")


def apply_publication_style(style: FigureStyle = FigureStyle()):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": list(style.font_family),
        "font.size": style.font_size,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": style.axes_linewidth,
        "xtick.major.width": style.axes_linewidth * 0.85,
        "ytick.major.width": style.axes_linewidth * 0.85,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.titlesize": style.font_size,
        "axes.labelsize": style.font_size,
        "legend.fontsize": style.font_size - 2,
        "xtick.labelsize": style.font_size - 1,
        "ytick.labelsize": style.font_size - 1,
        "savefig.facecolor": "white",
    })


def finalize_figure(fig, out_path, formats=("png", "pdf", "svg"), dpi=300, close=True):
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    saved = []
    for f in formats:
        p = out.with_suffix("." + f)
        fig.savefig(p, dpi=dpi, bbox_inches="tight", pad_inches=0.12, facecolor="white")
        saved.append(p)
        if f == "svg":
            s = p.read_text(encoding="utf-8")
            p.write_text(s.replace("Liberation Sans", "Arial"), encoding="utf-8")
    if close:
        plt.close(fig)
    return saved
