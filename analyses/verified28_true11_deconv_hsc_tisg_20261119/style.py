"""Publication styling utilities following .scientific-figure-making/SKILL.md."""
from dataclasses import dataclass
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    "gray_1": "#767676",
    "gray_2": "#4D4D4D",
    "gray_3": "#272727",
    "highlight": "#FFD700",
    "teal": "#42949E",
    "violet": "#9A4D8E",
}

CT_COL = {
    "T cell": "#3775BA",
    "NK cell": "#8BCF8B",
    "B cell": "#9A4D8E",
    "Monocyte": "#42949E",
    "Other": "#CFCECE",
}


@dataclass(frozen=True)
class FigureStyle:
    font_size: int = 13
    axes_linewidth: float = 1.8
    use_tex: bool = False
    font_family: tuple[str, ...] = ("Arial", "Helvetica", "DejaVu Sans", "sans-serif")


def apply_publication_style(style: FigureStyle | None = None) -> None:
    st = style or FigureStyle()
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": list(st.font_family),
        "font.size": st.font_size,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": st.axes_linewidth,
        "xtick.major.width": st.axes_linewidth * 0.85,
        "ytick.major.width": st.axes_linewidth * 0.85,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.usetex": st.use_tex,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def finalize_figure(fig, out_path, formats=("png", "pdf", "svg"), dpi=300, close=True, pad=0.06):
    out_p = pathlib.Path(out_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    saved = []
    for fmt in formats:
        fp = out_p.with_suffix(f".{fmt}")
        fig.savefig(fp, dpi=dpi, bbox_inches="tight", pad_inches=pad, facecolor="white")
        saved.append(fp)
    if close:
        plt.close(fig)
    return saved
