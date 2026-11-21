"""Nature Portfolio 1:1 print-scale house style module (180 mm two-column standard),
strictly implementing .scientific-figure-making/SKILL.md (figures4papers PALETTE, minimalist spines,
frameless legends, crisp black bar edges, vector PDF/SVG + 600-DPI combination-art PNG) with automated
typographic & boundary auditing.
"""
from dataclasses import dataclass
import logging
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# Canonical semantic PALETTE from .scientific-figure-making/references/api.md & design-theory.md
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
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "magenta": "#EA84DD",
    "gray_1": "#767676",
    "gray_2": "#4D4D4D",
    "gray_3": "#272727",
}

SRC_COL = {
    "Real": PALETTE["gray_3"],
    "PBISC": PALETTE["blue_main"],
    "scSemiProfiler": PALETTE["red_strong"],
}
SRC_MK = {"Real": "D", "PBISC": "o", "scSemiProfiler": "^"}

CT_COL = {
    "T cell": PALETTE["blue_secondary"],
    "NK cell": PALETTE["green_3"],
    "B cell": PALETTE["violet"],
    "Monocyte": PALETTE["teal"],
    "DC": PALETTE["red_2"],
    "Other": PALETTE["neutral"],
}

# Nature Portfolio exact two-column width: 180.0 mm = 7.0866 in (max page height 240-247 mm = 9.45-9.72 in)
NATURE_2COL_IN = 180.0 / 25.4  # 7.086614 in


@dataclass(frozen=True)
class FigureStyle:
    font_size: float = 6.5
    axes_linewidth: float = 0.85
    use_tex: bool = False
    font_family: tuple[str, ...] = ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans", "sans-serif")


def apply_publication_style(style: FigureStyle | None = None) -> None:
    st = style or FigureStyle()
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": list(st.font_family),
        "font.size": st.font_size,
        "axes.titlesize": st.font_size + 0.4,
        "axes.labelsize": st.font_size,
        "xtick.labelsize": st.font_size - 0.5,
        "ytick.labelsize": st.font_size - 0.5,
        "legend.fontsize": st.font_size - 0.6,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": st.axes_linewidth,
        "xtick.major.width": st.axes_linewidth * 0.85,
        "ytick.major.width": st.axes_linewidth * 0.85,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "xtick.major.pad": 2.0,
        "ytick.major.pad": 2.0,
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


def panel_label(ax, letter: str, x: float = -0.14, y: float = 1.08, size: float = 8.5) -> None:
    """Nature Portfolio standard: bold lowercase sans-serif letter (8-9 pt at 1x print scale)."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size, fontweight="bold",
            va="bottom", ha="left", color="#000000", clip_on=False)


def fig_panel_label(fig, x_fig: float, y_fig: float, letter: str, size: float = 8.5) -> None:
    fig.text(x_fig, y_fig, letter, fontsize=size, fontweight="bold", va="top", ha="left", color="#000000")


def panel_title(ax, text: str, size: float = 6.9, pad: float = 4.0, color: str | None = None) -> None:
    ax.set_title(text, loc="left", fontsize=size, fontweight="bold", pad=pad,
                 color=color or PALETTE["gray_3"])


def audit_figure_geometry(fig, name: str, min_pt: float = 5.0) -> dict:
    """Inspect all visible in-view Text artists at true 1x print scale for font size and canvas bounds."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w_px, h_px = fig.canvas.get_width_height()
    w_in, h_in = fig.get_size_inches()

    # Collect out-of-view tick label Text objects that Matplotlib instantiates outside xlim/ylim
    off_axis_ticks = set()
    for ax in fig.axes:
        xl0, xl1 = sorted(ax.get_xlim())
        yl0, yl1 = sorted(ax.get_ylim())
        xtol = 1e-6 * max(1.0, xl1 - xl0)
        ytol = 1e-6 * max(1.0, yl1 - yl0)
        for tk in ax.xaxis.get_major_ticks() + ax.xaxis.get_minor_ticks():
            loc = tk.get_loc()
            if loc < xl0 - xtol or loc > xl1 + xtol:
                off_axis_ticks.add(id(tk.label1))
                off_axis_ticks.add(id(tk.label2))
        for tk in ax.yaxis.get_major_ticks() + ax.yaxis.get_minor_ticks():
            loc = tk.get_loc()
            if loc < yl0 - ytol or loc > yl1 + ytol:
                off_axis_ticks.add(id(tk.label1))
                off_axis_ticks.add(id(tk.label2))

    texts = []
    for t in fig.findobj(mpl.text.Text):
        if id(t) in off_axis_ticks or not t.get_visible():
            continue
        s = t.get_text().strip()
        if not s:
            continue
        fs = float(t.get_fontsize())
        bb = t.get_window_extent(renderer=renderer)
        if bb.width <= 0 or bb.height <= 0:
            continue
        texts.append({
            "text": s[:40],
            "pt": round(fs, 2),
            "x0_frac": round(bb.x0 / w_px, 4),
            "x1_frac": round(bb.x1 / w_px, 4),
            "y0_frac": round(bb.y0 / h_px, 4),
            "y1_frac": round(bb.y1 / h_px, 4),
        })
    small = [x for x in texts if x["pt"] < min_pt]
    oob = [x for x in texts if (x["x0_frac"] < -0.002 or x["x1_frac"] > 1.002 or x["y0_frac"] < -0.002 or x["y1_frac"] > 1.002)]
    return {
        "figure": name,
        "width_mm": round(w_in * 25.4, 2),
        "height_mm": round(h_in * 25.4, 2),
        "n_text_spans": len(texts),
        "min_font_pt": min((x["pt"] for x in texts), default=0.0),
        "max_font_pt": max((x["pt"] for x in texts), default=0.0),
        "below_min_pt_count": len(small),
        "below_min_pt_examples": small[:5],
        "out_of_bounds_count": len(oob),
        "out_of_bounds_examples": oob[:5],
    }


def finalize_figure(fig, out_path, formats=("png", "pdf", "svg"), png_dpi: int = 600, vec_dpi: int = 450, close: bool = True):
    """Save exact 1:1 physical print dimensions (180 mm wide) without bbox_inches='tight' distortion."""
    out_p = pathlib.Path(out_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    audit = audit_figure_geometry(fig, out_p.name)
    saved = []
    for fmt in formats:
        fp = out_p.with_suffix(f".{fmt}")
        dpi_use = png_dpi if fmt == "png" else vec_dpi
        fig.savefig(fp, dpi=dpi_use, facecolor="white")
        if fmt == "svg":
            txt = fp.read_text(encoding="utf-8")
            fp.write_text(txt.replace("Liberation Sans", "Arial"), encoding="utf-8")
        saved.append(fp)
    if close:
        plt.close(fig)
    return saved, audit
