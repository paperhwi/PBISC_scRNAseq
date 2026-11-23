"""Build FigR2: Side-by-Side Comparison of Previous (Patient-Average) vs NEW Framework
(All-217-Draws Same-Day WHO Severity + Donor First-Draw Final Death) across:
  a. 43 immune subsets: Previous Patient-Mean Peak WHO rho vs NEW Same-Day WHO GEE beta (+wk)
  b. Mundlak within-patient vs between-patient severity slope decomposition
  c. Mortality dissociation: First-Draw (N=66) vs Patient-Mean vs Late (>7 d) draws
  d. Independent cohort replication on First-Draw Mortality (Real n=37 vs PBISC new n=29)
  e. Single-cell transcriptional programs (HSPC modules & T-cell ISG) under NEW framework
  f. Cohort-expansion power rescue (n=20 Real vs +29 PBISC new vs +29 shuffled; 300 sims)

Follows .scientific-figure-making/SKILL.md and nature_style.py (180.0 mm width, 600-DPI PNG +
vector PDF + editable SVG, automated oob == 0 audit).
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
from scipy import stats

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

LBL = {
    "B_ASC_dividing": "B ASC dividing",
    "B_ASC_IgG": "B ASC IgG",
    "B_ASC_IgA": "B ASC IgA",
    "B_ASC_IgM": "B ASC IgM",
    "NK CD56 bright": "NK CD56bright",
    "NK": "NK (CD56dim)",
    "HSPC": "HSPC",
    "T CD8 CM": "T CD8 CM",
    "T CD8 CTL": "T CD8 CTL",
    "CD14mono_dex": "CD14 mono dex",
    "CD14mono_IFN": "CD14 mono IFN",
    "CD16mono": "CD16 mono",
    "CD16mono_C1": "CD16 mono C1",
    "B_CD11c": "B CD11c",
    "T g/d": "T gd",
    "MAIT": "MAIT",
    "NKT": "NKT",
    "T CD4 CTL": "T CD4 CTL",
    "T CD4 cycling": "T CD4 cycling",
    "ASDC": "ASDC",
}


def disp(name: str) -> str:
    return LBL.get(name, name.replace("_", " "))


def build_figR2():
    apply_publication_style(FigureStyle(font_size=6.2, axes_linewidth=0.80))

    t2 = pd.read_csv(TAB_DIR / "T2_comp43_patientmean_vs_samedaySev_firstdrawDeath.csv")
    t2b = pd.read_csv(TAB_DIR / "T2b_comp43_combined_full_comparison.csv")
    t3 = pd.read_csv(TAB_DIR / "T3_new_framework_43_funnel.csv")
    t4 = pd.read_csv(TAB_DIR / "T4_power_rescue_new_framework.csv")
    t5 = pd.read_csv(TAB_DIR / "T5_HSPC_16modules_20genes_new_framework.csv")
    t6 = pd.read_csv(TAB_DIR / "T6_Tcell_ISG_13subsets_new_framework.csv")

    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.85))
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        width_ratios=[1.05, 1.05],
        height_ratios=[1.06, 1.08, 1.04],
        left=0.095, right=0.985, top=0.955, bottom=0.082,
        wspace=0.36, hspace=0.42,
    )

    # ------------------------------------------------------------------
    # PANEL a: 43 subsets Previous Patient-Mean rho vs NEW Same-Day WHO GEE beta
    # ------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    panel_label(ax_a, "a", x=-0.22, y=1.06, size=8.5)
    panel_title(ax_a, "Severity: Patient-mean Peak WHO vs Same-day WHO (43 subsets)", size=6.1, pad=5.0)

    ax_a.axhline(0, color=PALETTE["neutral"], lw=0.7, ls="--", zorder=0)
    ax_a.axvline(0, color=PALETTE["neutral"], lw=0.7, ls="--", zorder=0)

    x_pm = t2b["pm_peakWHO_rho"].to_numpy()
    y_sd = t2b["sd_gee_wk_beta"].to_numpy()
    q_pm = t2b["pm_peakWHO_q43"].to_numpy()
    q_sd = t2b["sd_gee_wk_p_q43"].to_numpy()

    r_align, p_align = stats.pearsonr(x_pm, y_sd)

    for i, row in t2b.iterrows():
        both_sig = (row["pm_peakWHO_q43"] < 0.05) and (row["sd_gee_wk_p_q43"] < 0.05)
        new_only = (row["pm_peakWHO_q43"] >= 0.05) and (row["sd_gee_wk_p_q43"] < 0.05)
        prev_only = (row["pm_peakWHO_q43"] < 0.05) and (row["sd_gee_wk_p_q43"] >= 0.05)
        if both_sig:
            c, s, z = PALETTE["blue_main"], 24, 4
        elif new_only:
            c, s, z = PALETTE["teal"], 22, 3
        elif prev_only:
            c, s, z = PALETTE["red_strong"], 26, 5
        else:
            c, s, z = PALETTE["neutral"], 13, 2
        ax_a.scatter(row["pm_peakWHO_rho"], row["sd_gee_wk_beta"], s=s,
                     facecolor=c, edgecolor="#222222" if z >= 3 else "#888888", lw=0.45, zorder=z)

    # Curated non-overlapping direct labels (fan layout)
    label_offsets_a = {
        "B_ASC_dividing": (-0.08, 0.072),
        "B_ASC_IgG": (-0.25, 0.055),
        "CD14mono_dex": (0.145, -0.015),
        "B_ASC_IgA": (0.145, -0.075),
        "HSPC": (-0.19, 0.035),
        "CD16mono_C1": (-0.17, 0.045),
        "T CD8 CTL": (0.04, 0.115),
        "NK": (-0.01, 0.105),
        "CD16mono": (0.22, -0.145),
        "T CD8 CM": (0.20, -0.025),
        "NK CD56 bright": (0.14, -0.075),
    }
    t2b_idx = t2b.set_index("feature")
    for feat, (dx, dy) in label_offsets_a.items():
        if feat not in t2b_idx.index:
            continue
        xv = float(t2b_idx.loc[feat, "pm_peakWHO_rho"])
        yv = float(t2b_idx.loc[feat, "sd_gee_wk_beta"])
        col_txt = PALETTE["red_strong"] if feat in ("T CD8 CTL", "NK") else PALETTE["gray_3"]
        ax_a.annotate(
            disp(feat), xy=(xv, yv), xytext=(xv + dx, yv + dy),
            fontsize=5.0, fontweight="bold" if feat in ("B_ASC_IgG", "B_ASC_dividing", "HSPC", "NK CD56 bright", "T CD8 CTL", "NK") else "normal",
            color=col_txt, ha="center", va="center",
            arrowprops=dict(arrowstyle="-", color="#666666", lw=0.45, shrinkA=1, shrinkB=2),
        )

    ax_a.set_xlim(-0.62, 0.68)
    ax_a.set_ylim(-0.39, 0.54)
    ax_a.set_xlabel("Previous: Patient-mean rho (vs Peak WHO, n=66)", fontsize=5.8)
    ax_a.set_ylabel("NEW: Same-day WHO GEE beta (+wk, 217 draws)", fontsize=5.8)
    ax_a.text(0.03, 0.96, f"r = {r_align:.2f}\nBoth q<0.05: {int(((q_pm<0.05)&(q_sd<0.05)).sum())}\nNew q<0.05: {int(((q_pm>=0.05)&(q_sd<0.05)).sum())}\nLost (early): {int(((q_pm<0.05)&(q_sd>=0.05)).sum())}",
              transform=ax_a.transAxes, va="top", ha="left", fontsize=5.0, color=PALETTE["gray_3"],
              bbox=dict(boxstyle="round,pad=0.22", facecolor="#F8FAFC", edgecolor="#CBD5E1", lw=0.5))

    leg_a = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["blue_main"], markeredgecolor="#222", ms=4.2, label="Both q<0.05"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["teal"], markeredgecolor="#222", ms=4.2, label="NEW only"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["red_strong"], markeredgecolor="#222", ms=4.2, label="Prev only (early)"),
    ]
    ax_a.legend(handles=leg_a, loc="lower right", frameon=False, fontsize=5.0, handletextpad=0.25)

    # ------------------------------------------------------------------
    # PANEL b: Mundlak Within-Patient vs Between-Patient Severity Slopes
    # ------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    panel_label(ax_b, "b", x=-0.25, y=1.06, size=8.5)
    panel_title(ax_b, "Mundlak decomposition: Within- vs Between-patient WHO slope", size=6.1, pad=5.0)

    mundlak_feats = [
        "B_ASC_IgG", "B_ASC_dividing", "CD14mono_dex", "B_ASC_IgA", "HSPC",
        "CD16mono_C1", "T g/d", "CD16mono", "T CD8 CM", "NK CD56 bright"
    ]
    sub_m = t2b_idx.loc[mundlak_feats].reset_index()
    y_pos = np.arange(len(sub_m))
    h = 0.35

    ax_b.axvline(0, color="#222222", lw=0.7, zorder=1)
    ax_b.barh(y_pos - h / 2, sub_m["mundlak_within_beta"], height=h, color=PALETTE["blue_main"],
              edgecolor="#000000", lw=0.55, label="Within-patient (same pt ΔWHO)", zorder=3)
    ax_b.barh(y_pos + h / 2, sub_m["mundlak_between_beta"], height=h, color=PALETTE["neutral"],
              edgecolor="#000000", lw=0.55, label="Between-patient (pt mean WHO)", zorder=3)

    for idx, row in sub_m.iterrows():
        bw = row["mundlak_within_beta"]
        pw = row["mundlak_within_p"]
        if pw < 0.05:
            star = "***" if pw < 1e-3 else ("**" if pw < 0.01 else "*")
            off = 0.025 if bw >= 0 else -0.025
            ax_b.text(bw + off, idx - h / 2, star, ha="left" if bw >= 0 else "right",
                      va="center", fontsize=5.8, fontweight="bold", color=PALETTE["blue_main"])

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels([disp(f) for f in sub_m["feature"]], fontsize=5.4)
    ax_b.invert_yaxis()
    ax_b.set_xlim(-0.72, 0.58)
    ax_b.set_xlabel("Mundlak Exchangeable GEE slope (+week, 217 draws)", fontsize=5.8)
    ax_b.legend(loc="upper left", frameon=False, fontsize=5.0, handlelength=1.1)

    # ------------------------------------------------------------------
    # PANEL c: Mortality Dissociation: First-Draw (N=66) vs Patient-Mean vs Late (>7 d)
    # ------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    panel_label(ax_c, "c", x=-0.22, y=1.06, size=8.5)
    panel_title(ax_c, "Mortality: First-draw prognosis vs Patient-mean vs Late (>7 d)", size=6.1, pad=5.0)

    dth_feats = [
        "T CD8 CTL", "NK", "B_CD11c", "NK CD56 bright", "T g/d", "MAIT",
        "CD14mono_IFN", "CD16mono_C1", "HSPC", "CD14mono_dex", "B_ASC_dividing", "B_ASC_IgG"
    ]
    sub_d = t2b_idx.loc[dth_feats].reset_index()
    yp_c = np.arange(len(sub_d))

    ax_c.axvline(0.5, color="#333333", lw=0.75, ls="--", zorder=1)
    # Divider separating baseline death predictors vs late/contemporaneous markers
    ax_c.axhline(7.5, color=PALETTE["gray_1"], lw=0.65, ls=":", zorder=1)
    ax_c.text(0.04, 7.22, "Baseline predictors ↑", ha="left", va="bottom", fontsize=5.0,
              color=PALETTE["blue_main"], fontweight="bold")
    ax_c.text(0.04, 7.78, "Late / severity markers ↓", ha="left", va="top", fontsize=5.0,
              color=PALETTE["red_strong"], fontweight="bold")

    for idx, row in sub_d.iterrows():
        a_fd = row["fd_AUC"]
        a_pm = row["pm_death_AUC"]
        a_lt = row["late_gt7_AUC"]
        # Connector segment showing shift from First-Draw to Late (>7d)
        ax_c.plot([a_fd, a_lt], [idx, idx], color="#CBD5E1", lw=1.1, zorder=2)
        ax_c.scatter(a_pm, idx, s=16, marker="D", facecolor=PALETTE["neutral"], edgecolor="#333333", lw=0.45, zorder=3)
        ax_c.scatter(a_lt, idx, s=18, marker="^", facecolor=PALETTE["red_strong"], edgecolor="#222222", lw=0.45, zorder=4)
        ax_c.scatter(a_fd, idx, s=26, marker="o", facecolor=PALETTE["blue_main"], edgecolor="#000000", lw=0.6, zorder=5)

        # Annotate adjWHO-surviving first-draw predictors
        if row["fd_p"] < 0.05:
            tag = "†" if row["fd_adjWHO_p"] < 0.05 else "*"
            xpos = min(a_fd, a_pm, a_lt) - 0.03 if a_fd < 0.5 else max(a_fd, a_pm, a_lt) + 0.03
            ha = "right" if a_fd < 0.5 else "left"
            ax_c.text(xpos, idx, f"P={row['fd_p']:.3f}{tag}", fontsize=5.0, fontweight="bold",
                      color=PALETTE["blue_main"], ha=ha, va="center")

    ax_c.set_yticks(yp_c)
    ax_c.set_yticklabels([disp(f) for f in sub_d["feature"]], fontsize=5.3)
    ax_c.invert_yaxis()
    ax_c.set_xlim(0.01, 0.99)
    ax_c.set_xticks([0.2, 0.35, 0.5, 0.65, 0.8, 0.95])
    ax_c.set_xlabel("Mortality AUC (TRUE_PUB_11; <0.5 depleted, >0.5 elevated)", fontsize=5.7)

    leg_c = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["blue_main"], markeredgecolor="#000", ms=4.5, label="NEW First draw (N=66)"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=PALETTE["neutral"], markeredgecolor="#333", ms=3.8, label="Prev Patient-mean"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=PALETTE["red_strong"], markeredgecolor="#222", ms=4.2, label="Late >7 d draws"),
    ]
    ax_c.legend(handles=leg_c, loc="lower left", frameon=False, fontsize=5.0, handletextpad=0.25)
    ax_c.text(0.02, -0.22, "† survives same-day WHO adjustment (P<0.05); * unadjusted P<0.05",
              transform=ax_c.transAxes, fontsize=5.0, color=PALETTE["gray_2"])

    # ------------------------------------------------------------------
    # PANEL d: Independent Cohort Replication on First-Draw Mortality (Real vs PBISC new)
    # ------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    panel_label(ax_d, "d", x=-0.22, y=1.06, size=8.5)
    panel_title(ax_d, "First-draw mortality replication: Real (n=37) vs PBISC (n=29)", size=6.1, pad=5.0)

    t2_r = t2[t2["set"] == "Real"].set_index("feature")
    t2_n = t2[t2["set"] == "PBISC new"].set_index("feature")
    common_f = sorted(set(t2_r.index) & set(t2_n.index))
    xr = np.array([t2_r.loc[f, "fd_AUC"] for f in common_f])
    yn = np.array([t2_n.loc[f, "fd_AUC"] for f in common_f])

    ax_d.axhline(0.5, color=PALETTE["neutral"], lw=0.7, ls="--", zorder=0)
    ax_d.axvline(0.5, color=PALETTE["neutral"], lw=0.7, ls="--", zorder=0)
    ax_d.plot([0.12, 0.88], [0.12, 0.88], color="#CBD5E1", lw=0.8, ls=":", zorder=0)

    # Shaded concordant quadrants
    ax_d.fill_between([0.10, 0.50], 0.10, 0.50, color="#EEF5FC", alpha=0.55, zorder=0)
    ax_d.fill_between([0.50, 0.90], 0.50, 0.90, color="#EEF5FC", alpha=0.55, zorder=0)

    for f, rv, nv in zip(common_f, xr, yn):
        comb_p = float(t2b_idx.loc[f, "fd_p"])
        if f in ("T CD8 CTL", "NK", "CD14mono_IFN", "B_CD11c", "NK CD56 bright"):
            c, s, z = PALETTE["blue_main"], 28, 5
        elif f == "HSPC":
            c, s, z = PALETTE["red_strong"], 28, 5
        elif comb_p < 0.08:
            c, s, z = PALETTE["teal"], 20, 4
        else:
            c, s, z = PALETTE["neutral"], 12, 2
        ax_d.scatter(rv, nv, s=s, facecolor=c, edgecolor="#222222" if z >= 4 else "#888888", lw=0.5, zorder=z)

    lbl_d = {
        "T CD8 CTL": (0.025, 0.068),
        "NK": (0.16, -0.022),
        "B_CD11c": (0.125, 0.015),
        "NK CD56 bright": (-0.01, 0.095),
        "CD14mono_IFN": (0.01, 0.095),
        "CD16mono_C1": (-0.08, 0.085),
        "HSPC": (0.06, 0.048),
    }
    for feat, (dx, dy) in lbl_d.items():
        rv = float(t2_r.loc[feat, "fd_AUC"])
        nv = float(t2_n.loc[feat, "fd_AUC"])
        col_t = PALETTE["red_strong"] if feat == "HSPC" else PALETTE["gray_3"]
        tag = f"{disp(feat)}\n({rv:.2f}, {nv:.2f})" if feat in ("T CD8 CTL", "NK", "CD14mono_IFN", "HSPC") else disp(feat)
        ax_d.annotate(
            tag, xy=(rv, nv), xytext=(rv + dx, nv + dy),
            fontsize=5.0, fontweight="bold" if feat in ("T CD8 CTL", "NK", "CD14mono_IFN", "HSPC") else "normal",
            color=col_t, ha="center", va="center",
            arrowprops=dict(arrowstyle="-", color="#666666", lw=0.45, shrinkA=1, shrinkB=2),
        )

    ax_d.set_xlim(0.11, 0.89)
    ax_d.set_ylim(0.11, 0.92)
    ax_d.set_xticks([0.2, 0.35, 0.5, 0.65, 0.8])
    ax_d.set_yticks([0.2, 0.35, 0.5, 0.65, 0.8])
    ax_d.set_xlabel("Real first-draw mortality AUC (n=37, 5 deaths)", fontsize=5.8)
    ax_d.set_ylabel("PBISC new first-draw mortality AUC (n=29, 6 deaths)", fontsize=5.8)
    ax_d.text(0.97, 0.04, "Blue quadrants:\nconcordant Real & PBISC\nfirst-draw predictors",
              transform=ax_d.transAxes, va="bottom", ha="right", fontsize=5.0,
              color=PALETTE["blue_main"], fontweight="bold")

    # ------------------------------------------------------------------
    # PANEL e: Single-Cell HSPC Modules & T-cell ISG under NEW Framework
    # ------------------------------------------------------------------
    ax_e = fig.add_subplot(gs[2, 0])
    panel_label(ax_e, "e", x=-0.22, y=1.06, size=8.5)
    panel_title(ax_e, "Single-cell HSPC modules & T-cell ISG (NEW framework)", size=6.1, pad=5.0)

    t5_c = t5[t5["set"] == "Combined"].set_index("feature")
    t6_c = t6[t6["set"] == "Combined"].set_index("subset")

    sc_items = [
        ("HSPC Cycle4", t5_c.loc["Cycle4", "sd_gee_wk_beta"], t5_c.loc["Cycle4", "sd_gee_wk_p"], t5_c.loc["Cycle4", "fd_AUC"], t5_c.loc["Cycle4", "fd_p"]),
        ("HSPC Ery5", t5_c.loc["Ery5", "sd_gee_wk_beta"], t5_c.loc["Ery5", "sd_gee_wk_p"], t5_c.loc["Ery5", "fd_AUC"], t5_c.loc["Ery5", "fd_p"]),
        ("HSPC Alarmin", t5_c.loc["Alarmin_Inflam", "sd_gee_wk_beta"], t5_c.loc["Alarmin_Inflam", "sd_gee_wk_p"], t5_c.loc["Alarmin_Inflam", "fd_AUC"], t5_c.loc["Alarmin_Inflam", "fd_p"]),
        ("HSPC HSC5", t5_c.loc["HSC5", "sd_gee_wk_beta"], t5_c.loc["HSC5", "sd_gee_wk_p"], t5_c.loc["HSC5", "fd_AUC"], t5_c.loc["HSC5", "fd_p"]),
        ("HSPC Myeloid4", t5_c.loc["Myeloid4", "sd_gee_wk_beta"], t5_c.loc["Myeloid4", "sd_gee_wk_p"], t5_c.loc["Myeloid4", "fd_AUC"], t5_c.loc["Myeloid4", "fd_p"]),
        ("ISG: All T cells", t6_c.loc["All T cells", "sd_gee_wk_beta"], t6_c.loc["All T cells", "sd_gee_wk_p"], t6_c.loc["All T cells", "fd_AUC"], t6_c.loc["All T cells", "fd_p"]),
        ("ISG: All CD4 T", t6_c.loc["All CD4 T", "sd_gee_wk_beta"], t6_c.loc["All CD4 T", "sd_gee_wk_p"], t6_c.loc["All CD4 T", "fd_AUC"], t6_c.loc["All CD4 T", "fd_p"]),
        ("ISG: All CD8 T", t6_c.loc["All CD8 T", "sd_gee_wk_beta"], t6_c.loc["All CD8 T", "sd_gee_wk_p"], t6_c.loc["All CD8 T", "fd_AUC"], t6_c.loc["All CD8 T", "fd_p"]),
        ("ISG: T CD8 CTL", t6_c.loc["T CD8 CTL", "sd_gee_wk_beta"], t6_c.loc["T CD8 CTL", "sd_gee_wk_p"], t6_c.loc["T CD8 CTL", "fd_AUC"], t6_c.loc["T CD8 CTL", "fd_p"]),
    ]
    yp_e = np.arange(len(sc_items))
    # Plot standardized effect sizes side-by-side: Same-Day WHO GEE beta vs First-Draw Death (2*(AUC - 0.5))
    sd_vals = [x[1] for x in sc_items]
    fd_eff = [2.0 * (x[3] - 0.5) for x in sc_items]
    h_e = 0.36

    ax_e.axvline(0, color="#222222", lw=0.7, zorder=1)
    ax_e.axhline(4.5, color=PALETTE["gray_1"], lw=0.65, ls=":", zorder=1)
    ax_e.barh(yp_e - h_e / 2, sd_vals, height=h_e, color=PALETTE["teal"], edgecolor="#000000", lw=0.55,
              label="Same-day WHO GEE beta (+wk)", zorder=3)
    ax_e.barh(yp_e + h_e / 2, fd_eff, height=h_e, color=PALETTE["red_strong"], edgecolor="#000000", lw=0.55,
              label="First-draw death 2*(AUC-0.5)", zorder=3)

    for idx, (nm, sb, sp, fa, fp) in enumerate(sc_items):
        if sp < 0.05:
            st = "***" if sp < 1e-3 else ("**" if sp < 0.01 else "*")
            ax_e.text(sb + (0.025 if sb >= 0 else -0.025), idx - h_e / 2, st,
                      ha="left" if sb >= 0 else "right", va="center", fontsize=5.4, fontweight="bold", color=PALETTE["teal"])
        fe = 2.0 * (fa - 0.5)
        if fp < 0.05:
            st = "***" if fp < 1e-3 else ("**" if fp < 0.01 else "*")
            ax_e.text(fe + (0.025 if fe >= 0 else -0.025), idx + h_e / 2, st,
                      ha="left" if fe >= 0 else "right", va="center", fontsize=5.4, fontweight="bold", color=PALETTE["red_strong"])

    ax_e.set_yticks(yp_e)
    ax_e.set_yticklabels([x[0] for x in sc_items], fontsize=5.3)
    ax_e.invert_yaxis()
    ax_e.set_xlim(-0.55, 0.66)
    ax_e.set_xlabel("Effect size (GEE beta or Rank-biserial 2*(AUC-0.5))", fontsize=5.7)
    ax_e.legend(loc="lower left", frameon=False, fontsize=5.0, handlelength=1.0)

    # ------------------------------------------------------------------
    # PANEL f: Cohort-Expansion Power Rescue (n=20 Real vs +29 PBISC new vs +29 Shuffled)
    # ------------------------------------------------------------------
    ax_f = fig.add_subplot(gs[2, 1])
    panel_label(ax_f, "f", x=-0.22, y=1.06, size=8.5)
    panel_title(ax_f, "PBISC power rescue under NEW framework (300 sims)", size=6.1, pad=5.0)

    pw_show = [
        ("B_ASC_IgG", "sameday_sev", "B ASC IgG\n(SD-WHO)", "fdr"),
        ("HSPC", "sameday_sev", "HSPC\n(SD-WHO)", "fdr"),
        ("NK CD56 bright", "sameday_sev", "NK CD56br\n(SD-WHO)", "fdr"),
        ("T CD8 CM", "sameday_sev", "T CD8 CM\n(SD-WHO)", "fdr"),
        ("T CD8 CTL", "firstdraw_dth", "T CD8 CTL\n(1st-Death)", "nom"),
        ("NK", "firstdraw_dth", "NK\n(1st-Death)", "nom"),
        ("CD14mono_IFN", "firstdraw_dth", "CD14 IFN\n(1st-Death)", "nom"),
    ]
    t4_idx = t4.set_index(["feature", "endpoint"])
    x_f = np.arange(len(pw_show))
    w_f = 0.25

    v_real, v_plus, v_shuf = [], [], []
    for feat, ep, _, metric in pw_show:
        row = t4_idx.loc[(feat, ep)]
        v_real.append(100.0 * float(row[f"real_{metric}"]))
        v_plus.append(100.0 * float(row[f"plus_{metric}"]))
        v_shuf.append(100.0 * float(row[f"shuf_{metric}"]))

    ax_f.axvline(3.5, color=PALETTE["gray_1"], lw=0.7, ls=":", zorder=1)
    b1 = ax_f.bar(x_f - w_f, v_real, width=w_f, color=PALETTE["gray_2"], edgecolor="#000000", lw=0.55, label="Real only (n=20)", zorder=3)
    b2 = ax_f.bar(x_f, v_plus, width=w_f, color=PALETTE["blue_main"], edgecolor="#000000", lw=0.55, label="+29 PBISC new", zorder=3)
    b3 = ax_f.bar(x_f + w_f, v_shuf, width=w_f, color=PALETTE["neutral"], edgecolor="#000000", lw=0.55, label="+29 Shuffled null", zorder=3)

    for i, vp in enumerate(v_plus):
        ax_f.text(i, vp + 2.2, f"{vp:.0f}%", ha="center", va="bottom", fontsize=5.0, fontweight="bold", color=PALETTE["blue_main"])

    ax_f.set_xticks(x_f)
    ax_f.set_xticklabels([p[2] for p in pw_show], fontsize=5.0)
    ax_f.set_ylabel("Discovery power (% of 300 sims)", fontsize=5.8)
    ax_f.set_ylim(0, 120)
    ax_f.set_yticks([0, 25, 50, 75, 100])
    ax_f.text(1.5, 114, "Same-day WHO (q43<0.05)", ha="center", va="top", fontsize=5.0, fontweight="bold", color=PALETTE["gray_3"])
    ax_f.text(5.0, 114, "First-draw Death (P<0.05)", ha="center", va="top", fontsize=5.0, fontweight="bold", color=PALETTE["gray_3"])
    ax_f.legend(loc="upper center", bbox_to_anchor=(0.5, -0.21), ncol=3, frameon=False, fontsize=5.0, handlelength=1.0, columnspacing=0.8)

    saved, audit = finalize_figure(fig, FIG_DIR / "FigR2_PatientMean_vs_SameDaySev_FirstDrawDeath")
    # Also sync directly to package figures/SuppFig11_Primary_SameDay_FirstDraw_vs_PatientMean.{png,pdf,svg} & AUDIT_FIG6_7_SFIG6_11.json
    import json, shutil
    out_root = ROOT.parent
    pkg_fig = out_root / "figures"
    pkg_fig.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        shutil.copy2(FIG_DIR / f"FigR2_PatientMean_vs_SameDaySev_FirstDrawDeath.{ext}",
                     pkg_fig / f"SuppFig11_Primary_SameDay_FirstDraw_vs_PatientMean.{ext}")
    audit_pkg = dict(audit)
    audit_pkg["figure"] = "SuppFig11_Primary_SameDay_FirstDraw_vs_PatientMean"
    audit_json = out_root / "AUDIT_FIG6_7_SFIG6_11.json"
    if audit_json.exists():
        recs = json.loads(audit_json.read_text(encoding="utf-8"))
        by_k = {r["figure"]: r for r in recs}
        by_k["SuppFig11_Primary_SameDay_FirstDraw_vs_PatientMean"] = {k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in audit_pkg.items()}
        audit_json.write_text(json.dumps(list(by_k.values()), indent=2), encoding="utf-8")
    print("FigR2 / SuppFig11 saved:", saved)
    print("FigR2 / SuppFig11 audit:", audit_pkg)
    return saved, audit_pkg


if __name__ == "__main__":
    build_figR2()
