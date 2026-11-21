"""Build FigR1: Complete 66-Patient Longitudinal Blood-Draw Schedule Swimmer Plot,
Same-Day vs Peak WHO Severity Discordance, and Real vs PBISC Sampling-Count Asymmetry.

Follows .scientific-figure-making/SKILL.md and nature_style.py (180.0 mm two-column width,
minimalist spines, semantic PALETTE, 600-DPI PNG + vector PDF + editable SVG, oob == 0 audit).
"""
import os
import pathlib
import sys
import tempfile
import warnings

os.environ.setdefault("MPLCONFIGDIR", str(pathlib.Path(tempfile.gettempdir()) / "mpl_pbisc_nature"))
pathlib.Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from nature_style import (
    NATURE_2COL_IN,
    PALETTE,
    FigureStyle,
    apply_publication_style,
    finalize_figure,
    panel_label,
    panel_title,
)

ROOT = pathlib.Path(r"F:\PBISC\deliverables\nature_submission_fig4_6_20261120\draw_sameday_sev_firstdraw_death_20261120")
TAB_DIR = ROOT / "tables"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

WHO_COL = {
    0: "#8BCF8B",  # mild (green_3)
    1: "#3775BA",  # moderate (blue_secondary)
    2: "#E68A3E",  # severe (warm amber-orange, high contrast with blue/red)
    3: "#B64342",  # critical (red_strong)
}
WHO_NAME = {0: "Mild", 1: "Moderate", 2: "Severe", 3: "Critical"}


def build_figR1():
    apply_publication_style(FigureStyle(font_size=6.3, axes_linewidth=0.80))

    t1 = pd.read_csv(TAB_DIR / "T1_all68_patients_draw_schedule.csv")
    t1b = pd.read_csv(TAB_DIR / "T1b_all217_primary_draws_detail.csv")
    t1c = pd.read_csv(TAB_DIR / "T1c_cohort_draw_summary.csv")

    prim_pts = t1[t1["primary_source"] != "Excluded_CLL"].copy()
    group_order = ["Real-only (10)", "Shared-27 (Real used)", "PBISC new (29)"]
    group_labels = {
        "Real-only (10)": "Real-only scRNA (n = 10; 47 draws; 1 death)",
        "Shared-27 (Real used)": "Shared-27 [Real used] (n = 27; 91 draws; 4 deaths)",
        "PBISC new (29)": "PBISC new [bulk-only] (n = 29; 79 draws; 6 deaths)",
    }
    group_bg = {
        "Real-only (10)": "#F5F7FA",
        "Shared-27 (Real used)": "#EEF7F7",
        "PBISC new (29)": "#EDF3FA",
    }
    group_bar = {
        "Real-only (10)": PALETTE["gray_2"],
        "Shared-27 (Real used)": PALETTE["teal"],
        "PBISC new (29)": PALETTE["blue_main"],
    }

    # Build ordered y-positions (top-to-bottom) with header gaps between cohorts
    y_rows = []
    cur_y = 0.0
    block_spans = []
    for g in group_order:
        sub = prim_pts[prim_pts["cohort_group"] == g].copy()
        sub = sub.sort_values(
            ["death_TRUE_PUB11", "peak_WHO_ord", "first_draw_WHO_ord", "first_draw_day", "n_primary_draws"],
            ascending=[False, False, False, True, False],
        ).reset_index(drop=True)
        header_y = cur_y
        cur_y += 1.35
        y_start = cur_y
        for _, r in sub.iterrows():
            y_rows.append((cur_y, g, r))
            cur_y += 1.0
        y_end = cur_y - 1.0
        block_spans.append((g, header_y, y_start, y_end))
        cur_y += 0.9

    total_h_units = cur_y

    # Layout: 180.0 mm x 218.0 mm (7.0866 x 8.58 in)
    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.58))
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        width_ratios=[1.48, 1.0],
        height_ratios=[1.05, 1.05, 1.05],
        left=0.072, right=0.985, top=0.932, bottom=0.068,
        wspace=0.26, hspace=0.40,
    )

    # ------------------------------------------------------------------
    # PANEL a: Full 66-patient longitudinal blood-draw swimmer plot
    # ------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[:, 0])
    panel_label(ax_a, "a", x=-0.13, y=1.032, size=8.5)
    panel_title(ax_a, "Blood-draw schedule across 66 COVID-19+ patients (217 draws)", size=6.5, pad=22.0)

    # Shaded cohort bands + headers
    for g, hdr_y, y0, y1 in block_spans:
        ax_a.axhspan(y0 - 0.55, y1 + 0.55, color=group_bg[g], zorder=0, lw=0)
        ax_a.add_patch(Rectangle((-1.85, y0 - 0.45), 0.65, (y1 - y0) + 0.90,
                                 facecolor=group_bar[g], edgecolor="none", zorder=2))
        ax_a.text(-1.05, hdr_y + 0.25, group_labels[g], fontsize=5.7, fontweight="bold",
                  color=group_bar[g], va="center", ha="left")

    # Day 7 guide
    ax_a.axvline(7.0, color=PALETTE["gray_1"], ls="--", lw=0.7, alpha=0.65, zorder=1)
    ax_a.text(7.35, -0.78, "Day 7 (<=7 d)", fontsize=5.2, color=PALETTE["gray_2"], va="bottom", ha="left")

    # Column headers for right-margin annotation strip
    ax_a.text(32.8, -0.78, "Peak", fontsize=5.3, fontweight="bold", color=PALETTE["gray_3"], ha="center", va="bottom")
    ax_a.text(35.3, -0.78, "n", fontsize=5.3, fontweight="bold", color=PALETTE["gray_3"], ha="center", va="bottom")
    ax_a.text(37.7, -0.78, "Out.", fontsize=5.3, fontweight="bold", color=PALETTE["gray_3"], ha="center", va="bottom")

    # Plot each patient's draws
    draws_by_pt = {pid: grp.sort_values("days") for pid, grp in t1b.groupby("patient_id")}
    yticks, ylabels = [], []
    for y, g, r in y_rows:
        pid = r["patient_id"]
        sub_d = draws_by_pt[pid]
        days = sub_d["days"].to_numpy(float)
        sd_who = sub_d["temp_severity_ord"].to_numpy(int)
        died = int(r["death_TRUE_PUB11"])
        peak_ord = int(r["peak_WHO_ord"])

        # Highlight fatal patient row subtly
        if died == 1:
            ax_a.axhspan(y - 0.46, y + 0.46, color="#FCE8E6", alpha=0.72, zorder=1, lw=0)

        # Horizontal track connecting draws
        if len(days) > 1:
            ax_a.plot([days.min(), days.max()], [y, y], color=PALETTE["gray_1"], lw=0.75, zorder=2)

        # Individual draws colored by same-day WHO; first draw highlighted with thick black ring
        for idx, (d_val, w_val) in enumerate(zip(days, sd_who)):
            is_first = (idx == 0)
            ax_a.scatter(
                d_val, y,
                s=22 if is_first else 14,
                marker="o",
                facecolor=WHO_COL[w_val],
                edgecolor="#000000" if is_first else "#333333",
                linewidth=1.05 if is_first else 0.35,
                zorder=4 if is_first else 3,
            )

        # Right-strip: Peak WHO square, draw count, death marker
        ax_a.scatter(32.8, y, s=16, marker="s", facecolor=WHO_COL[peak_ord],
                     edgecolor="#222222", linewidth=0.4, zorder=4)
        ax_a.text(35.3, y, f"{len(days)}", fontsize=5.2, ha="center", va="center",
                  color=PALETTE["gray_3"], fontweight="bold" if len(days) == 1 else "normal")
        if died == 1:
            ax_a.text(37.7, y, "Died", fontsize=5.1, ha="center", va="center",
                      color=PALETTE["red_strong"], fontweight="bold")
        else:
            ax_a.text(37.7, y, "Surv", fontsize=5.0, ha="center", va="center", color="#777777")

        yticks.append(y)
        ylabels.append(f"{pid}")

    ax_a.set_yticks(yticks)
    ax_a.set_yticklabels(ylabels, fontsize=5.1)
    # Color fatal patient ytick labels red/bold
    for lbl_obj, (_, _, r) in zip(ax_a.get_yticklabels(), y_rows):
        if int(r["death_TRUE_PUB11"]) == 1:
            lbl_obj.set_color(PALETTE["red_strong"])
            lbl_obj.set_fontweight("bold")

    ax_a.set_xlim(-2.1, 39.8)
    ax_a.set_xticks([0, 5, 10, 15, 20, 25, 30])
    ax_a.set_ylim(total_h_units - 0.1, -1.6)
    ax_a.set_xlabel("Days from first positive SARS-CoV-2 swab", fontsize=6.2)
    ax_a.tick_params(axis="y", length=0, pad=3.2)

    # Clean legend sitting cleanly ABOVE plot axes and BELOW panel title
    leg_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=WHO_COL[k], markeredgecolor="#333333",
               markeredgewidth=0.4, markersize=4.3, label=WHO_NAME[k]) for k in range(4)
    ] + [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#FFFFFF", markeredgecolor="#000000",
               markeredgewidth=1.1, markersize=4.8, label="First draw")
    ]
    ax_a.legend(handles=leg_handles, loc="lower left", bbox_to_anchor=(-0.02, 0.998),
                ncol=5, frameon=False, fontsize=5.2, handletextpad=0.25, columnspacing=0.65)

    # ------------------------------------------------------------------
    # PANEL b: Same-Day WHO vs Peak WHO 4x4 discordance matrix (217 draws)
    # ------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    panel_label(ax_b, "b", x=-0.20, y=1.08, size=8.5)
    panel_title(ax_b, "Same-day WHO vs Peak WHO across 217 draws", size=6.4, pad=5.0)

    mat = np.zeros((4, 4), dtype=int)
    for _, row in t1b.iterrows():
        mat[int(row["severity_ord"]), int(row["temp_severity_ord"])] += 1
    n_disc = int((t1b["severity_ord"] != t1b["temp_severity_ord"]).sum())
    pct_disc = 100.0 * n_disc / len(t1b)

    im = ax_b.imshow(mat, cmap="Blues", aspect="auto", vmin=0, vmax=mat.max() * 1.15)
    for i in range(4):
        for j in range(4):
            val = mat[i, j]
            if val == 0 and j > i:
                txt = "—"
                col_t = "#BBBBBB"
            else:
                txt = f"{val}"
                col_t = "white" if val > 55 else "#111111"
            ax_b.text(j, i, txt, ha="center", va="center", fontsize=6.4,
                      fontweight="bold" if (i != j and val > 0) else "normal", color=col_t)
            if i == j:
                ax_b.add_patch(Rectangle((j - 0.48, i - 0.48), 0.96, 0.96,
                                         fill=False, edgecolor="#111111", lw=0.9))

    ax_b.set_xticks(range(4))
    ax_b.set_xticklabels(["Mild", "Mod", "Sev", "Crit"], fontsize=5.8)
    ax_b.set_yticks(range(4))
    ax_b.set_yticklabels(["Mild\n(23 pts)", "Mod\n(14 pts)", "Sev\n(17 pts)", "Crit\n(12 pts)"], fontsize=5.6)
    ax_b.set_xlabel("Same-day WHO severity (at blood draw)", fontsize=6.0)
    ax_b.set_ylabel("Patient Peak WHO severity", fontsize=6.0)
    ax_b.spines["right"].set_visible(True)
    ax_b.spines["top"].set_visible(True)
    ax_b.spines["right"].set_linewidth(0.6)
    ax_b.spines["top"].set_linewidth(0.6)
    ax_b.text(2.48, 0.28, f"{n_disc}/217 draws\n({pct_disc:.1f}%) differ\nfrom Peak WHO",
              ha="center", va="center", fontsize=5.2,
              color=PALETTE["red_strong"], fontweight="bold",
              bbox=dict(boxstyle="round,pad=0.28", facecolor="#FFF5F5", edgecolor=PALETTE["red_strong"], lw=0.6))

    # ------------------------------------------------------------------
    # PANEL c: Draw-count & first-draw-day asymmetry (Real vs PBISC new, Survived vs Died)
    # ------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 1])
    panel_label(ax_c, "c", x=-0.20, y=1.08, size=8.5)
    panel_title(ax_c, "Draw-count asymmetry motivates first-draw death", size=6.4, pad=5.0)

    cats = ["Real\nSurv\n(n=32)", "Real\nDied\n(n=5)", "PBISC\nSurv\n(n=23)", "PBISC\nDied\n(n=6)"]
    subgroups = [
        t1b[(t1b["arm"] == "Real") & (t1b["fat11"] == 0)],
        t1b[(t1b["arm"] == "Real") & (t1b["fat11"] == 1)],
        t1b[(t1b["arm"] == "PBISC new") & (t1b["fat11"] == 0)],
        t1b[(t1b["arm"] == "PBISC new") & (t1b["fat11"] == 1)],
    ]
    bar_cols = [PALETTE["neutral"], PALETTE["red_2"], PALETTE["blue_secondary"], PALETTE["red_strong"]]

    means = []
    tot_draws = []
    rng_j = np.random.default_rng(7)
    for idx, sg in enumerate(subgroups):
        nd_pt = sg.groupby("patient_id").size().to_numpy(float)
        means.append(nd_pt.mean())
        tot_draws.append(int(nd_pt.sum()))
        ax_c.bar(idx, nd_pt.mean(), width=0.58, color=bar_cols[idx], edgecolor="#000000", linewidth=0.75, zorder=2)
        xj = idx + rng_j.uniform(-0.16, 0.16, size=len(nd_pt))
        ax_c.scatter(xj, nd_pt, s=11, color="#111111", alpha=0.65, linewidths=0, zorder=3)
        ax_c.text(idx, nd_pt.max() + 0.35, f"{nd_pt.mean():.1f}/pt\n({int(nd_pt.sum())} dr)",
                  ha="center", va="bottom", fontsize=5.3, fontweight="bold", color=PALETTE["gray_3"])

    ax_c.set_xticks(range(4))
    ax_c.set_xticklabels(cats, fontsize=5.6)
    ax_c.set_ylabel("Draws per patient (mean + dots)", fontsize=6.0)
    ax_c.set_ylim(0, 8.6)
    ax_c.set_yticks([0, 2, 4, 6, 8])
    ax_c.annotate("4/6 died after\nsingle draw!", xy=(3, 1.38), xytext=(2.85, 5.15),
                  fontsize=5.4, fontweight="bold", color=PALETTE["red_strong"], ha="center",
                  arrowprops=dict(arrowstyle="->", color=PALETTE["red_strong"], lw=0.8))

    # ------------------------------------------------------------------
    # PANEL d: First-draw sampling day distribution & framework schematic
    # ------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[2, 1])
    panel_label(ax_d, "d", x=-0.20, y=1.08, size=8.5)
    panel_title(ax_d, "First blood-draw timing (baseline prognosis)", size=6.4, pad=5.0)

    fd = t1b[t1b["is_first_draw"]].copy()
    groups_fd = [
        ("Real Surv\n(n=32)", fd[(fd["arm"] == "Real") & (fd["fat11"] == 0)]["days"].to_numpy(), PALETTE["neutral"]),
        ("Real Died\n(n=5)", fd[(fd["arm"] == "Real") & (fd["fat11"] == 1)]["days"].to_numpy(), PALETTE["red_2"]),
        ("PBISC Surv\n(n=23)", fd[(fd["arm"] == "PBISC new") & (fd["fat11"] == 0)]["days"].to_numpy(), PALETTE["blue_secondary"]),
        ("PBISC Died\n(n=6)", fd[(fd["arm"] == "PBISC new") & (fd["fat11"] == 1)]["days"].to_numpy(), PALETTE["red_strong"]),
    ]

    ax_d.axhline(7.0, color=PALETTE["gray_1"], ls="--", lw=0.75, alpha=0.7, zorder=1)
    ax_d.text(3.48, 7.4, "Day 7\n(53/66 pts;\n10/11 died)", fontsize=5.1, color=PALETTE["gray_2"],
              ha="right", va="bottom")

    for idx, (lbl, vals, col) in enumerate(groups_fd):
        bp = ax_d.boxplot(
            [vals], positions=[idx], widths=0.48, patch_artist=True, showfliers=False,
            medianprops=dict(color="#000000", lw=1.1),
            boxprops=dict(facecolor=col, edgecolor="#000000", lw=0.7),
            whiskerprops=dict(color="#000000", lw=0.7),
            capprops=dict(color="#000000", lw=0.7),
        )
        xj = idx + rng_j.uniform(-0.14, 0.14, size=len(vals))
        ax_d.scatter(xj, vals, s=11, color="#111111", alpha=0.65, linewidths=0, zorder=4)
        ax_d.text(idx, 19.3, f"med {np.median(vals):.1f}d", ha="center", va="bottom",
                  fontsize=5.2, color=PALETTE["gray_3"])

    ax_d.set_xticks(range(4))
    ax_d.set_xticklabels([g[0] for g in groups_fd], fontsize=5.5)
    ax_d.set_ylabel("First blood-draw day (from +swab)", fontsize=6.0)
    ax_d.set_ylim(-1.0, 22.5)
    ax_d.set_yticks([0, 5, 10, 15, 20])
    ax_d.text(0.5, -0.23, "37/37 Real & 16/29 PBISC first draws <=7 d (10/11 deaths <=7 d)",
              transform=ax_d.transAxes, ha="center", va="top", fontsize=5.1, color=PALETTE["gray_2"])

    saved, audit = finalize_figure(fig, FIG_DIR / "FigR1_Patient_Swimmer_Draw_Schedule_66")
    print("FigR1 saved:", saved)
    print("FigR1 audit:", audit)


if __name__ == "__main__":
    build_figR1()
