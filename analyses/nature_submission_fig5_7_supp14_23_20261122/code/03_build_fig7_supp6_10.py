"""Build Nature Portfolio 180-mm Main Figure 7 and Supplementary Figures 6-10
aligned 1:1 with the finalized methodological framework:
  - Main Figure 7 (a-j): Same-day WHO severity (equal total weight per patient w_ij = 1/n_i across 276 acute samples)
                         and last-eligible-sample mortality (N = 66: 55 survived vs 11 died, Mann-Whitney AUC),
                         each BH-adjusted across 43 cell types (q < 0.05), plus 4 external studies analyzed by
                         the exact same criteria and pooled via DerSimonian-Laird random-effects (k >= 2).
  - Supplementary Fig. 6: Complete Sepsis (255 samples) & ESKD COVID-19 (107 donors / 362 total / 276 acute) metadata tables.
  - Supplementary Fig. 7: 66-patient longitudinal swimmer schedule across all 276 eligible acute draws (marking last draws).
  - Supplementary Fig. 8: 43-cell-type sensitivity (raw % vs INT vs count vs CBC vs CLL) & Official R 4.3.2 deconvolution.
  - Supplementary Fig. 9: Per-study forest plots & tables across the 4 external COVID-19 scRNA-seq cohorts (k >= 2).
  - Supplementary Fig. 10: Single-cell T-cell 7-ISG & HSPC 16-module / 20-gene transcriptional programs.
"""
import json
import os
import pathlib
import sys
import tempfile
import warnings

os.environ.setdefault("MPLCONFIGDIR", str(pathlib.Path(tempfile.gettempdir()) / "mpl_pbisc_fig7"))
pathlib.Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_curve
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

from nature_style import (
    NATURE_2COL_IN,
    PALETTE,
    SRC_COL,
    FigureStyle,
    apply_publication_style,
    fig_panel_label,
    finalize_figure,
    panel_label,
    panel_title,
)

ROOT = pathlib.Path(r"F:\PBISC")
OUT_ROOT = ROOT / "deliverables" / "nature_submission_fig5_7_20261122"
FIG_OUT = OUT_ROOT / "figures"
DATA_OUT = OUT_ROOT / "data"
FIG_OUT.mkdir(parents=True, exist_ok=True)

WHO_COL = {0: "#AADCA9", 1: "#FFD700", 2: "#E9A6A1", 3: "#B64342"}
WHO_LBL = {0: "Mild", 1: "Mod.", 2: "Sev.", 3: "Crit."}
SHORT_CT = {
    "NK CD56 bright": "NK CD56bright",
    "B_ASC_dividing": "B ASC div.",
    "B_ASC_IgG": "B ASC IgG",
    "B_ASC_IgA": "B ASC IgA",
    "B_ASC_IgM": "B ASC IgM",
    "B_ASC_total": "B ASC (total)",
    "CD14mono_dex": "CD14mono dex",
    "CD14mono_IFN": "CD14mono IFN",
    "CD16mono_C1": "CD16mono C1",
    "T CD4 cycling": "CD4 cycling",
    "T CD4 Helper": "CD4 Helper",
    "T CD4 CTL": "CD4 CTL",
    "T CD8 CTL": "CD8 CTL",
    "T CD8 CM": "CD8 CM",
    "T CD8 Naive": "CD8 Naive",
    "T CD4 Naive": "CD4 Naive",
    "B_switched_mem": "B sw. mem",
    "B_non-switched_mem": "B non-sw. mem",
    "B_naive": "B naive",
}


def sh_ct(name: str) -> str:
    return SHORT_CT.get(name, name)


# =============================================================================
# MAIN FIGURE 7 (a-j)
# =============================================================================
def build_figure7():
    df276 = pd.read_csv(DATA_OUT / "fig7_eskd_276_acute_samples_metadata_and_43pct.csv")
    df66 = pd.read_csv(DATA_OUT / "fig7_eskd_66_last_samples.csv")
    eskd = pd.read_csv(DATA_OUT / "fig7_eskd_43celltypes_sameday_sev_lastdraw_death.csv")
    ext_meta = pd.read_csv(DATA_OUT / "fig7_external4_random_effects_meta.csv")

    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.36), facecolor="white")

    # ---- Panel a: Cohort & analytical design (66 haemodialysis patients, 276 acute samples) ----
    ax_a = fig.add_axes([0.022, 0.865, 0.545, 0.120])
    ax_a.set_xlim(0, 100); ax_a.set_ylim(0, 100); ax_a.axis("off")
    ax_a.add_patch(FancyBboxPatch((0.5, 48.0), 99.0, 48.5, boxstyle="round,pad=0.25,rounding_size=2.0",
                                  fc="#F6F9FC", ec=PALETTE["blue_main"], lw=0.85))
    ax_a.text(2.5, 91.5, "Haemodialysis COVID-19: 66 patients · 276 acute draws",
              fontsize=5.9, fontweight="bold", color=PALETTE["blue_main"], va="top")
    ax_a.text(2.5, 78.0,
              "• Real: 138 draws (37 pts) · PBISC: 138 draws (79 new + 59 Shared-27)\n"
              "• Severity: same-day WHO vs %, equal patient weight (1/n per patient)\n"
              "• Mortality: last acute draw/patient (55 surv. vs 11 died, AUC); BH43<0.05",
              fontsize=5.05, color=PALETTE["gray_3"], va="top", linespacing=1.22)

    # Stacked horizontal bar of the 276 acute draws & 66 last draws
    segs_276 = [
        ("Real (35)", 35, PALETTE["gray_3"]),
        ("Shared-27 Real (103)", 103, PALETTE["gray_1"]),
        ("Shared-27 PBISC (59)", 59, PALETTE["blue_secondary"]),
        ("PBISC new 29 pts (79)", 79, PALETTE["blue_main"]),
    ]
    x_left = 2.0
    tot_w = 96.0
    for lbl, cnt, col in segs_276:
        w_seg = tot_w * (cnt / 276.0)
        ax_a.add_patch(Rectangle((x_left, 10.0), w_seg, 21.0, fc=col, ec="white", lw=0.7))
        if w_seg >= 11.0:
            ax_a.text(x_left + w_seg / 2, 20.5, lbl, ha="center", va="center", fontsize=5.05, color="white", fontweight="bold")
        x_left += w_seg
    ax_a.text(2.0, 35.5, "276 acute draws (0 same-day overlap) → 66 unique last draws (24R + 42P)",
              fontsize=5.1, fontweight="bold", color=PALETTE["gray_2"], va="bottom")
    fig_panel_label(fig, 0.008, 0.992, "a")

    # ---- Panel b: BH43 < 0.05 discovery counts (Real vs Combined 276/66) ----
    ax_b = fig.add_axes([0.635, 0.870, 0.340, 0.108])
    r_real = eskd[eskd["stratum"] == "Real_37pts_138draws"]
    r_comb = eskd[eskd["stratum"] == "Combined_66pts_276draws_rawPct"]
    r_int = eskd[eskd["stratum"] == "Combined_66pts_276draws_withinSourceINT"]

    cats = ["Same-day WHO\n(BH43 q < 0.05)", "Last-draw death\n(BH43 q < 0.05)"]
    v_real = [int((r_real["q_sev_bh43"] < 0.05).sum()), int((r_real["q_death_bh43"] < 0.05).sum())]
    v_comb = [int((r_comb["q_sev_bh43"] < 0.05).sum()), int((r_comb["q_death_bh43"] < 0.05).sum())]
    v_int = [int((r_int["q_sev_bh43"] < 0.05).sum()), int((r_int["q_death_bh43"] < 0.05).sum())]

    xb = np.arange(2); wb = 0.25
    b1 = ax_b.bar(xb - wb, v_real, wb, color=SRC_COL["Real"], edgecolor="k", lw=0.6, label="Real (37 pts)")
    b2 = ax_b.bar(xb, v_comb, wb, color=SRC_COL["PBISC"], edgecolor="k", lw=0.6, label="Comb. % (66 pts)")
    b3 = ax_b.bar(xb + wb, v_int, wb, color=PALETTE["blue_secondary"], edgecolor="k", lw=0.6, hatch="//", label="Comb. INT")
    for bars in (b1, b2, b3):
        for rect in bars:
            hv = int(rect.get_height())
            ax_b.text(rect.get_x() + rect.get_width() / 2, hv + 0.5, str(hv), ha="center", va="bottom", fontsize=5.4, fontweight="bold")
    ax_b.set_xticks(xb); ax_b.set_xticklabels(cats, fontsize=5.5)
    ax_b.set_ylim(0, 24.5); ax_b.set_yticks([0, 10, 20]); ax_b.set_ylabel("Significant types (/43)", fontsize=5.6)
    ax_b.legend(loc="upper right", fontsize=5.1, handlelength=0.85, labelspacing=0.14)
    panel_title(ax_b, "FDR discovery rescue (+PBISC)", size=6.2, pad=3)
    fig_panel_label(fig, 0.585, 0.992, "b")

    # ---- Panel c: Joint 43-cell-type landscape (Same-day WHO rho_w vs Last-draw mortality AUC) ----
    ax_c = fig.add_axes([0.082, 0.652, 0.390, 0.148])
    rc = r_comb.copy().reset_index(drop=True)
    sig_s = rc["q_sev_bh43"] < 0.05
    sig_m = rc["q_death_bh43"] < 0.05
    dual = sig_s & sig_m
    only_s = sig_s & (~sig_m)
    ns = (~sig_s) & (~sig_m)

    ax_c.axhline(0.5, color=PALETTE["neutral"], lw=0.65, ls=(0, (2, 2)), zorder=0)
    ax_c.axvline(0.0, color=PALETTE["neutral"], lw=0.65, ls=(0, (2, 2)), zorder=0)
    ax_c.scatter(rc.loc[ns, "rho_w"], rc.loc[ns, "auc_death"], s=15, color=PALETTE["neutral"], edgecolor=PALETTE["gray_1"], lw=0.4, alpha=0.8, label="q ≥ 0.05 (25)", zorder=2)
    ax_c.scatter(rc.loc[only_s, "rho_w"], rc.loc[only_s, "auc_death"], s=24, color=PALETTE["blue_secondary"], edgecolor="k", lw=0.45, label=f"WHO q<0.05 ({int(only_s.sum())})", zorder=3)
    ax_c.scatter(rc.loc[dual, "rho_w"], rc.loc[dual, "auc_death"], s=42, marker="D", color=PALETTE["red_strong"], edgecolor="k", lw=0.65, label=f"Both q<0.05 ({int(dual.sum())})", zorder=5)

    lbl_pos = {
        "HSPC": (0.14, 0.85, "left"),
        "T g/d": (-0.28, 0.11, "center"),
        "CD14mono_dex": (0.50, 0.44, "right"),
        "B_ASC_IgG": (0.50, 0.53, "right"),
        "B_ASC_dividing": (0.15, 0.77, "left"),
        "B_ASC_IgA": (0.03, 0.79, "center"),
        "CD16mono_C1": (-0.06, 0.74, "right"),
        "NK": (-0.42, 0.20, "center"),
        "NK CD56 bright": (-0.42, 0.34, "center"),
        "T CD8 CM": (-0.42, 0.44, "center"),
        "T CD8 CTL": (-0.41, 0.57, "center"),
        "B_CD11c": (-0.02, 0.17, "left"),
    }
    rc_idx = rc.set_index("feature")
    for fn, (tx, ty, ha) in lbl_pos.items():
        if fn in rc_idx.index:
            rw = rc_idx.loc[fn]
            is_d = (rw["q_sev_bh43"] < 0.05) and (rw["q_death_bh43"] < 0.05)
            col_t = PALETTE["red_strong"] if is_d else PALETTE["blue_main"]
            ax_c.annotate(sh_ct(fn), xy=(rw["rho_w"], rw["auc_death"]), xytext=(tx, ty),
                          fontsize=5.2, fontweight="bold" if is_d else "normal", color=col_t, ha=ha, va="center",
                          arrowprops=dict(arrowstyle="-", lw=0.4, color=PALETTE["gray_1"], shrinkA=0, shrinkB=2))
    ax_c.set_xlim(-0.52, 0.52); ax_c.set_ylim(0.06, 0.91)
    ax_c.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4]); ax_c.set_yticks([0.2, 0.5, 0.8])
    ax_c.set_xlabel("Same-day WHO weighted ρ (276 draws)", fontsize=5.6)
    ax_c.set_ylabel("Last-draw death AUC (66 pts)", fontsize=5.6)
    ax_c.legend(loc="lower right", fontsize=5.1, handletextpad=0.15, labelspacing=0.14)
    panel_title(ax_c, "43-cell-type severity vs mortality landscape", size=6.2, pad=3)
    fig_panel_label(fig, 0.008, 0.822, "c")

    # ---- Panel d: Stratified cross-source concordance (Real vs PBISC vs Combined) ----
    ax_d1 = fig.add_axes([0.585, 0.652, 0.178, 0.148])
    ax_d2 = fig.add_axes([0.792, 0.652, 0.180, 0.148])
    foc8 = ["HSPC", "B_ASC_IgG", "B_ASC_dividing", "CD14mono_dex", "T CD8 CTL", "T CD8 CM", "NK CD56 bright", "T g/d"]
    r_pb = eskd[eskd["stratum"] == "PBISC_56pts_138draws"].set_index("feature")
    r_re = r_real.set_index("feature")
    r_co = r_comb.set_index("feature")
    yd = np.arange(len(foc8))

    for off, tag, tab, col, mk in [(-0.22, "Real", r_re, SRC_COL["Real"], "D"),
                                   (0.0, "PBISC", r_pb, PALETTE["blue_secondary"], "s"),
                                   (0.22, "Combined", r_co, SRC_COL["PBISC"], "o")]:
        v = tab.reindex(foc8)["rho_w"].values
        lo = v - tab.reindex(foc8)["ci_lo_rho"].values
        hi = tab.reindex(foc8)["ci_hi_rho"].values - v
        ax_d1.errorbar(v, yd + off, xerr=[np.maximum(0, lo), np.maximum(0, hi)], fmt=mk, ms=3.4, color=col,
                       elinewidth=0.65, capsize=1.0, label=tag)
    ax_d1.axvline(0.0, color="k", lw=0.65, ls=(0, (2, 2)))
    ax_d1.set_yticks(yd); ax_d1.set_yticklabels([sh_ct(x) for x in foc8], fontsize=5.3)
    ax_d1.invert_yaxis(); ax_d1.set_xlim(-0.65, 0.65); ax_d1.set_xticks([-0.4, 0.0, 0.4])
    ax_d1.set_xlabel("Same-day WHO weighted ρ", fontsize=5.5)
    ax_d1.legend(loc="lower left", bbox_to_anchor=(0.0, 1.11), ncol=3, fontsize=5.1, handletextpad=0.15, columnspacing=0.6)
    panel_title(ax_d1, "Severity concordance", size=6.1, pad=3)

    for off, tag, tab, col, mk in [(-0.22, "Real (37)", r_re, SRC_COL["Real"], "D"),
                                   (0.0, "PBISC (56)", r_pb, PALETTE["blue_secondary"], "s"),
                                   (0.22, "Comb. (66)", r_co, SRC_COL["PBISC"], "o")]:
        v = tab.reindex(foc8)["auc_death"].values
        lo = v - tab.reindex(foc8)["ci_lo_auc"].values
        hi = tab.reindex(foc8)["ci_hi_auc"].values - v
        ax_d2.errorbar(v, yd + off, xerr=[np.maximum(0, lo), np.maximum(0, hi)], fmt=mk, ms=3.4, color=col,
                       elinewidth=0.65, capsize=1.0, label=tag)
    ax_d2.axvline(0.5, color="k", lw=0.65, ls=(0, (2, 2)))
    ax_d2.set_yticks(yd); ax_d2.set_yticklabels([])
    ax_d2.invert_yaxis(); ax_d2.set_xlim(0.02, 0.98); ax_d2.set_xticks([0.2, 0.5, 0.8])
    ax_d2.set_xlabel("Last-draw death AUC", fontsize=5.5)
    panel_title(ax_d2, "Mortality concordance", size=6.1, pad=3)
    fig_panel_label(fig, 0.502, 0.822, "d")

    # ---- Panels e & f: Same-day WHO grade distributions across all 276 acute draws ----
    rng = np.random.default_rng(20261122)

    def draw_who_strip(ax, feat, ttl, show_yl=False, show_leg=False):
        rw = r_co.loc[feat]
        w_wts = df276["equal_patient_weight"].to_numpy(float)
        y_raw = df276[feat].to_numpy(float)
        # Use log1p(pct) for skewed rare subsets or raw % cleanly labeled
        use_log = feat in ("HSPC", "B_ASC_IgG", "B_ASC_dividing", "NK CD56 bright", "T g/d")
        yv = np.log10(y_raw + 0.01) if use_log else y_raw
        for g_ord in range(4):
            m = (df276["temp_severity_ord"].to_numpy(int) == g_ord)
            mr = m & (df276["source"].to_numpy(str) == "Real")
            mp = m & (df276["source"].to_numpy(str) == "PBISC")
            ax.scatter(g_ord - 0.12 + rng.uniform(-0.07, 0.07, mr.sum()), yv[mr], s=6.5, color=SRC_COL["Real"], alpha=0.65, lw=0, rasterized=True, label="Real" if (g_ord == 0 and show_leg) else None)
            ax.scatter(g_ord + 0.12 + rng.uniform(-0.07, 0.07, mp.sum()), yv[mp], s=6.5, color=SRC_COL["PBISC"], alpha=0.65, lw=0, rasterized=True, label="PBISC" if (g_ord == 0 and show_leg) else None)
            # Equal-patient-weighted mean across draws in this WHO grade
            wm = float(np.average(yv[m], weights=w_wts[m]))
            ax.plot([g_ord - 0.24, g_ord + 0.24], [wm, wm], color=PALETTE["red_strong"], lw=1.35, zorder=5)
        ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=5.4)
        if use_log:
            ax.set_yticks([-2, -1, 0, 1]); ax.set_yticklabels(["0", "0.1", "1", "10"], fontsize=5.3)
            ax.set_ylim(-2.2, 1.55)
            if show_yl:
                ax.set_ylabel("Cells (% PBMC, log)", fontsize=5.5)
        else:
            ax.set_ylim(-2.0, float(yv.max()) * 1.18)
            if show_yl:
                ax.set_ylabel("Cells (% PBMC)", fontsize=5.5)
        q_str = f"q={rw['q_sev_bh43']:.1e}" if rw["q_sev_bh43"] < 0.001 else f"q={rw['q_sev_bh43']:.3f}"
        ax.text(0.03, 0.96, f"ρ = {rw['rho_w']:+.2f} ({q_str})", transform=ax.transAxes, ha="left", va="top",
                fontsize=5.2, fontweight="bold", color=PALETTE["blue_main"],
                bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.4))
        if show_leg:
            ax.legend(loc="lower right", fontsize=5.1, handletextpad=0.12, markerscale=1.2)
        panel_title(ax, ttl, size=6.0, pad=2.5)

    fig_panel_label(fig, 0.008, 0.602, "e")
    ax_e1 = fig.add_axes([0.074, 0.448, 0.185, 0.128]); draw_who_strip(ax_e1, "HSPC", "HSPC vs WHO", show_yl=True)
    ax_e2 = fig.add_axes([0.288, 0.448, 0.185, 0.128]); draw_who_strip(ax_e2, "B_ASC_IgG", "B ASC IgG vs WHO", show_yl=False, show_leg=True)

    fig_panel_label(fig, 0.502, 0.602, "f")
    ax_f1 = fig.add_axes([0.566, 0.448, 0.188, 0.128]); draw_who_strip(ax_f1, "NK CD56 bright", "NK CD56bright vs WHO", show_yl=True)
    ax_f2 = fig.add_axes([0.785, 0.448, 0.188, 0.128]); draw_who_strip(ax_f2, "T CD8 CTL", "CD8 CTL vs WHO", show_yl=False)

    # ---- Panel g: Last-eligible-sample mortality distributions (55 survived vs 11 died, N=66) ----
    fig_panel_label(fig, 0.008, 0.402, "g")
    for idx_g, (feat, ttl_g) in enumerate([("HSPC", "Last draw · HSPC"), ("T g/d", "Last draw · T g/d")]):
        ax_g = fig.add_axes([0.074 + idx_g * 0.214, 0.252, 0.185, 0.125])
        rw = r_co.loc[feat]
        yv = np.log10(df66[feat].to_numpy(float) + 0.01)
        for d_flag in (0, 1):
            m = df66["fatal_pub11"].to_numpy(int) == d_flag
            mr = m & (df66["source"].to_numpy(str) == "Real")
            mp = m & (df66["source"].to_numpy(str) == "PBISC")
            ax_g.scatter(d_flag - 0.11 + rng.uniform(-0.07, 0.07, mr.sum()), yv[mr], s=11, marker="D", color=SRC_COL["Real"], alpha=0.82, lw=0, label="Real (24)" if (d_flag == 0 and idx_g == 1) else None)
            ax_g.scatter(d_flag + 0.11 + rng.uniform(-0.07, 0.07, mp.sum()), yv[mp], s=13, marker="o", color=SRC_COL["PBISC"], alpha=0.85, lw=0, label="PBISC (42)" if (d_flag == 0 and idx_g == 1) else None)
            med = float(np.median(yv[m]))
            q25, q75 = np.quantile(yv[m], [0.25, 0.75])
            ax_g.plot([d_flag - 0.22, d_flag + 0.22], [med, med], color=PALETTE["red_strong"], lw=1.4, zorder=5)
            ax_g.plot([d_flag, d_flag], [q25, q75], color=PALETTE["red_strong"], lw=1.0, zorder=4)
        ax_g.set_xticks([0, 1]); ax_g.set_xticklabels(["Surv.\n(n=55)", "Died\n(n=11)"], fontsize=5.4)
        ax_g.set_yticks([-2, -1, 0, 1]); ax_g.set_yticklabels(["0", "0.1", "1", "10"], fontsize=5.3)
        ax_g.set_ylim(-2.15, 1.45)
        if idx_g == 0:
            ax_g.set_ylabel("Last-draw % PBMC", fontsize=5.5)
        ax_g.text(0.03, 0.96, f"AUC = {rw['auc_death']:.2f}\nq43 = {rw['q_death_bh43']:.3f}",
                  transform=ax_g.transAxes, ha="left", va="top", fontsize=5.15, fontweight="bold",
                  color=PALETTE["red_strong"], bbox=dict(fc="white", ec="none", alpha=0.84, pad=0.3))
        if idx_g == 1:
            ax_g.legend(loc="upper right", fontsize=5.0, handletextpad=0.12)
        panel_title(ax_g, ttl_g, size=6.0, pad=2.5)

    # ---- Panel h: Last-draw ROC curves (Real n=37 vs Combined n=66) ----
    ax_h = fig.add_axes([0.568, 0.252, 0.405, 0.125])
    y_true66 = df66["fatal_pub11"].to_numpy(int)
    df37_last = df276[df276["source"] == "Real"].sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last()
    y_true37 = df37_last["fatal_pub11"].to_numpy(int)

    fpr_h66, tpr_h66, _ = roc_curve(y_true66, df66["HSPC"].to_numpy(float))
    fpr_h37, tpr_h37, _ = roc_curve(y_true37, df37_last["HSPC"].to_numpy(float))
    fpr_g66, tpr_g66, _ = roc_curve(y_true66, -df66["T g/d"].to_numpy(float))

    ax_h.plot([0, 1], [0, 1], color=PALETTE["neutral"], lw=0.7, ls=(0, (2, 2)))
    ax_h.plot(fpr_h37, tpr_h37, color=SRC_COL["Real"], lw=1.0, ls="--", label="HSPC Real n=37 (AUC 0.73, q=0.64)")
    ax_h.plot(fpr_h66, tpr_h66, color=PALETTE["red_strong"], lw=1.35, label="HSPC Comb. n=66 (AUC 0.80, q=0.039)")
    ax_h.plot(fpr_g66, tpr_g66, color=SRC_COL["PBISC"], lw=1.25, label="Low T g/d n=66 (AUC 0.82, q=0.039)")
    ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1.02); ax_h.set_xticks([0.0, 0.5, 1.0]); ax_h.set_yticks([0.0, 0.5, 1.0])
    ax_h.set_xlabel("1 − Specificity (FPR)", fontsize=5.5); ax_h.set_ylabel("Sensitivity (TPR)", fontsize=5.5)
    ax_h.legend(loc="lower right", fontsize=5.0, handlelength=1.1, labelspacing=0.15)
    panel_title(ax_h, "Last-draw mortality ROC rescue (n=37 → n=66)", size=6.1, pad=2.5)
    fig_panel_label(fig, 0.502, 0.402, "h")

    # ---- Panel i: External 4-study Random-Effects Meta-Analysis — Same-day WHO severity (k >= 2) ----
    ax_i = fig.add_axes([0.148, 0.045, 0.325, 0.145])
    em_exact = ext_meta[ext_meta["is_exact_43_subset"] == True].copy()
    em_sev_top = em_exact[em_exact["q_sev_re"] < 0.05].sort_values("rho_re").reset_index(drop=True)
    yi = np.arange(len(em_sev_top))
    ax_i.axvline(0.0, color="k", lw=0.7, ls=(0, (2, 2)))
    lo_i = em_sev_top["rho_re"] - em_sev_top["ci_lo_rho_re"]
    hi_i = em_sev_top["ci_hi_rho_re"] - em_sev_top["rho_re"]
    ax_i.errorbar(em_sev_top["rho_re"], yi, xerr=[lo_i, hi_i], fmt="o", ms=3.8, color=PALETTE["blue_main"],
                  elinewidth=0.8, capsize=1.3, label="Ext. RE meta (k=2–4)")
    ax_i.scatter(em_sev_top["eskd_rho_w"], yi, marker="D", s=16, facecolor="white", edgecolor=PALETTE["red_strong"],
                 lw=0.85, zorder=5, label="ESKD 276 draws")
    ax_i.set_yticks(yi)
    ax_i.set_yticklabels([f"{sh_ct(r.feature)} (k={int(r.k_studies)})" for r in em_sev_top.itertuples()], fontsize=5.2)
    ax_i.set_xlim(-0.62, 0.65); ax_i.set_xticks([-0.4, 0.0, 0.4])
    ax_i.set_xlabel("Same-day WHO ρ (RE meta, BH q<0.05)", fontsize=5.5)
    ax_i.legend(loc="lower right", fontsize=5.0, handletextpad=0.12, labelspacing=0.14)
    panel_title(ax_i, "External 4 studies · severity meta", size=6.1, pad=2.5)
    fig_panel_label(fig, 0.008, 0.205, "i")

    # ---- Panel j: External 4-study Random-Effects Meta-Analysis — Last-draw mortality AUC (k >= 2) ----
    ax_j = fig.add_axes([0.625, 0.045, 0.348, 0.145])
    foc_m_ext = ["CD16mono", "T g/d", "MAIT", "T CD8 CM", "NK CD56 bright", "T CD8 CTL", "HSPC", "B_ASC_dividing"]
    em_m = em_exact.set_index("feature").reindex(foc_m_ext).reset_index()
    yj = np.arange(len(em_m))
    ax_j.axvline(0.5, color="k", lw=0.7, ls=(0, (2, 2)))
    lo_j = em_m["auc_re"] - em_m["ci_lo_auc_re"]
    hi_j = em_m["ci_hi_auc_re"] - em_m["auc_re"]
    ax_j.errorbar(em_m["auc_re"], yj, xerr=[lo_j, hi_j], fmt="o", ms=3.8, color=PALETTE["blue_main"],
                  elinewidth=0.8, capsize=1.3, label="Ext. RE meta (22 deaths)")
    ax_j.scatter(em_m["eskd_auc_last66"], yj, marker="D", s=16, facecolor="white", edgecolor=PALETTE["red_strong"],
                 lw=0.85, zorder=5, label="ESKD last draw (n=66)")
    for yi_k, r in enumerate(em_m.itertuples()):
        star = "*" if r.q_death_re < 0.05 else ("†" if r.p_death_re < 0.05 else "")
        if star:
            ax_j.text(0.95, yi_k, star, va="center", ha="center", fontsize=5.6, fontweight="bold", color=PALETTE["blue_main"])
    ax_j.set_yticks(yj)
    ax_j.set_yticklabels([f"{sh_ct(r.feature)} (k={int(r.k_studies)})" for r in em_m.itertuples()], fontsize=5.2)
    ax_j.invert_yaxis()
    ax_j.set_xlim(0.02, 0.99); ax_j.set_xticks([0.2, 0.5, 0.8])
    ax_j.set_xlabel("Last-draw death AUC (* q<0.05, † P<0.05)", fontsize=5.5)
    ax_j.legend(loc="lower right", fontsize=5.0, handletextpad=0.12, labelspacing=0.14)
    panel_title(ax_j, "External 4 studies · mortality meta", size=6.1, pad=2.5)
    fig_panel_label(fig, 0.502, 0.205, "j")

    return finalize_figure(fig, FIG_OUT / "Figure7_SameDay_Severity_and_LastDraw_Mortality")


# =============================================================================
# SUPPLEMENTARY FIGURE 6: COHORT & METADATA TABLES (SEPSIS 255 & ESKD 107 / 276)
# =============================================================================
def _draw_pub_table(ax, col_headers, rows, col_widths, title_txt):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    ax.text(0.0, 99.0, title_txt, fontsize=6.4, fontweight="bold", color=PALETTE["gray_3"], va="top")
    n_r = len(rows)
    top_y = 90.0
    row_h = min(8.2, 84.0 / (n_r + 1.2))
    # Header
    ax.add_patch(Rectangle((0, top_y - row_h), 100, row_h, fc=PALETTE["blue_main"], ec="none"))
    x = 0.0
    for ch, cw in zip(col_headers, col_widths):
        ax.text(x + 1.0, top_y - row_h / 2, ch, fontsize=5.35, fontweight="bold", color="white", va="center", ha="left")
        x += cw
    for ri, r in enumerate(rows):
        y0 = top_y - (ri + 2) * row_h
        fc = "#F5F8FB" if ri % 2 == 1 else "#FFFFFF"
        if ri == n_r - 1 and ("Total" in r[0] or "Combined" in r[0]):
            fc = "#E8EFF7"
        ax.add_patch(Rectangle((0, y0), 100, row_h, fc=fc, ec="#DDE3EA", lw=0.45))
        x = 0.0
        is_tot = ("Total" in r[0] or "Combined" in r[0])
        for ci, (val, cw) in enumerate(zip(r, col_widths)):
            ax.text(x + 1.0, y0 + row_h / 2, str(val), fontsize=5.2,
                    fontweight="bold" if (ci == 0 or is_tot) else "normal",
                    color=PALETTE["gray_3"], va="center", ha="left")
            x += cw


def build_suppfig_metadata(snum=19):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.20), facecolor="white")

    # Panel a: Bacterial Sepsis cohort metadata (255 samples, 0 donor overlap)
    ax_a = fig.add_axes([0.03, 0.62, 0.94, 0.34])
    sep_hdr = ["Dataset / Stratum", "Donors", "Samples", "Cells", "Clinical groups (n)", "Age / Sex & QC verification"]
    sep_w = [23, 10, 10, 13, 22, 22]
    sep_rows = [
        ["Real scRNA-seq (Adult)", "108", "108", "633,072", "93 Sepsis / 8 ICU / 7 HC", "Adult primary scope; 82 paired w/ bulk"],
        ["Real scRNA-seq (Paed.)", "25", "25", "147,461", "17 Sepsis / 4 ICU / 4 HC", "Paed. reference audit (Fig. 6c)"],
        ["PBISC bulk-only (Adult)", "121", "121", "484,000", "68 Sepsis / 17 ICU / 36 HC", "K=4,000/sample; 0 overlap w/ Real"],
        ["PBISC bulk-only (Paed.)", "1", "1", "4,000", "1 Paed. ICU (PIC06)", "Excluded from adult primary scope"],
        ["scSemiProfiler baseline", "34", "34", "51,000", "26 Sepsis / 3 ICU / 5 HC", "Adult held-out baseline"],
        ["Total Sepsis (Real+PBISC)", "255", "255", "1,268,533", "0 donor overlap verified", "5-fold LoRA adapted (GSE279448/52)"],
    ]
    _draw_pub_table(ax_a, sep_hdr, sep_rows, sep_w, "Bacterial Sepsis cohort metadata & dataset alignment (Ye et al., Nat. Immunol. 2026)")
    fig_panel_label(fig, 0.008, 0.985, "a")

    # Panel b: Imperial ESKD COVID-19 haemodialysis cohort metadata (107 donors / 362 total / 276 eligible acute draws)
    ax_b = fig.add_axes([0.03, 0.14, 0.94, 0.44])
    cov_hdr = ["ESKD Cohort Stratum", "Patients", "Acute draws", "Last draws", "Deaths (TRUE_PUB_11)", "Same-day WHO & longitudinal scope"]
    cov_w = [24, 10, 12, 12, 16, 26]
    cov_rows = [
        ["Real-only acute COVID+", "10", "35", "10 (Real)", "2 / 10 died", "Wave 2; 100% same-day WHO present"],
        ["Shared-27 Real acute", "27", "103", "14 (Real)", "3 / 27 died", "Different visit days from Shared PBISC"],
        ["Shared-27 PBISC acute", "(same 27)", "59", "13 (PBISC)", "(same 3 / 27)", "Complementary bulk-only acute visits"],
        ["PBISC-only new acute (W1)", "29", "79", "29 (PBISC)", "6 / 29 died", "Wave 1 bulk-only; 0 scRNA-seq"],
        ["Excluded CLL donors", "2", "(7 excl.)", "Excluded", "C116, C139", "B-cell malignancy excluded pre-analysis"],
        ["Negative / recovery / QC<200", "39 neg+", "(79 non-acute)", "Control/rec.", "—", "Used in Fig. 5 case vs control only"],
        ["Combined Fig. 7 Primary", "66", "276", "66 (24R+42P)", "11 died / 55 surv.", "Equal pt weight w_ij=1/n_i; 0 day ties"],
    ]
    _draw_pub_table(ax_b, cov_hdr, cov_rows, cov_w, "Imperial ESKD haemodialysis COVID-19 cohort metadata & 276-draw accounting (Fig. 7)")
    fig_panel_label(fig, 0.008, 0.595, "b")

    fig.text(0.03, 0.04,
             "Notes: Every one of the 276 eligible acute samples (condition=positive, primary 0–30 d scope, ≥200 cells, non-CLL) has complete same-day WHO grade (276/276)\n"
             "and days_from_first_positive_swab (276/276), zero same-day Real/PBISC duplicates, and a unique un-tied last eligible acute sample per patient (N = 66).",
             fontsize=5.3, color=PALETTE["gray_2"], va="bottom")

    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_Cohort_and_Sample_Metadata_Tables")


# =============================================================================
# SUPPLEMENTARY FIGURE 20: 66-PATIENT SWIMMER SCHEDULE ACROSS ALL 276 ACUTE DRAWS
# =============================================================================
def build_suppfig_swimmer(snum=20):
    df276 = pd.read_csv(DATA_OUT / "fig7_eskd_276_acute_samples_metadata_and_43pct.csv")
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.25), facecolor="white")

    # Sort 66 patients: Deceased (11) at top, Survivors (55) below, ordered by max swab day
    pt_sum = df276.groupby("patient_id").agg(
        fatal=("fatal_pub11", "max"),
        max_day=("days_from_first_positive_swab", "max"),
        min_day=("days_from_first_positive_swab", "min"),
        n_draws=("sample_id", "count"),
    ).reset_index()
    pt_sum = pt_sum.sort_values(["fatal", "max_day", "n_draws"], ascending=[False, True, False]).reset_index(drop=True)

    # Split into 2 balanced columns (33 patients left, 33 patients right) for crystal-clear legibility!
    gs = fig.add_gridspec(1, 2, wspace=0.28, left=0.10, right=0.97, top=0.91, bottom=0.07)
    for col_i in range(2):
        ax = fig.add_subplot(gs[0, col_i])
        sub_pts = pt_sum.iloc[col_i * 33 : (col_i + 1) * 33].reset_index(drop=True)
        for yi, pr in sub_pts.iterrows():
            pid = pr["patient_id"]
            g = df276[df276["patient_id"] == pid]
            ax.plot([pr["min_day"], pr["max_day"]], [yi, yi], color="#D0D7DE", lw=0.9, zorder=1)
            for _, dr in g.iterrows():
                mk = "s" if dr["source"] == "Real" else "o"
                fc = WHO_COL[int(dr["temp_severity_ord"])]
                is_last = bool(dr["is_last_sample"])
                ax.scatter(dr["days_from_first_positive_swab"], yi, marker=mk,
                           s=26 if is_last else 16, facecolor=fc,
                           edgecolor="#000000" if is_last else PALETTE["gray_2"],
                           lw=1.15 if is_last else 0.45, zorder=4 if is_last else 3)
            if pr["fatal"] == 1:
                ax.scatter(pr["max_day"] + 1.4, yi, marker="*", s=32, color=PALETTE["red_strong"], zorder=5)
        if col_i == 0:
            ax.axhline(10.5, color=PALETTE["red_strong"], lw=0.75, ls="--")
            ax.text(34.5, 0.3, "Deceased (n=11)", fontsize=5.3, fontweight="bold", color=PALETTE["red_strong"], ha="right", va="center")
        ax.set_yticks(np.arange(len(sub_pts)))
        lbls = [f"{r.patient_id} ({int(r.n_draws)})" for r in sub_pts.itertuples()]
        ax.set_yticklabels(lbls, fontsize=5.2)
        for tick_lbl, r in zip(ax.get_yticklabels(), sub_pts.itertuples()):
            if r.fatal == 1:
                tick_lbl.set_color(PALETTE["red_strong"])
                tick_lbl.set_fontweight("bold")
        ax.invert_yaxis()
        ax.set_xlim(-2, 35.5); ax.set_xticks([0, 7, 14, 21, 28, 35])
        ax.set_xlabel("Days from first positive swab", fontsize=5.7)
        panel_title(ax, f"Patients {col_i*33+1}–{(col_i+1)*33} (draws in parens)", size=6.2, pad=3)
        panel_label(ax, "a" if col_i == 0 else "b", x=-0.22, y=1.02)

    leg_h = [
        Patch(fc=WHO_COL[0], ec="k", lw=0.4, label="WHO Mild"),
        Patch(fc=WHO_COL[1], ec="k", lw=0.4, label="Moderate"),
        Patch(fc=WHO_COL[2], ec="k", lw=0.4, label="Severe"),
        Patch(fc=WHO_COL[3], ec="k", lw=0.4, label="Critical"),
        Line2D([], [], marker="s", ls="", mfc="white", mec=PALETTE["gray_2"], ms=4.2, label="Real (138)"),
        Line2D([], [], marker="o", ls="", mfc="white", mec=PALETTE["gray_2"], ms=4.2, label="PBISC (138)"),
        Line2D([], [], marker="o", ls="", mfc="white", mec="k", mew=1.2, ms=5.0, label="Last draw (66)"),
        Line2D([], [], marker="*", ls="", color=PALETTE["red_strong"], ms=6.0, label="Died (11)"),
    ]
    fig.legend(handles=leg_h, loc="upper center", bbox_to_anchor=(0.52, 0.992), ncol=4, fontsize=5.4,
               handletextpad=0.25, columnspacing=0.8)
    fig.text(0.5, 0.012,
             "All 276 eligible acute draws across 66 non-CLL haemodialysis patients. Bold black-bordered symbol marks each patient's unique last eligible draw used for mortality AUC.",
             ha="center", fontsize=5.3, color=PALETTE["gray_2"])

    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_Patient_Swimmer_276_Acute_Draws")


# =============================================================================
# SUPPLEMENTARY FIGURE 21: 43-CELL-TYPE SENSITIVITY & OFFICIAL R 4.3.2 DECONVOLUTION
# =============================================================================
def build_suppfig_sensitivity(snum=21):
    sens = pd.read_csv(DATA_OUT / "supp8_sensitivity_count_cbc_cll_276framework.csv")
    deconv = pd.read_csv(DATA_OUT / "supp8_officialR_deconv_12methods.csv")
    hspc_dec = pd.read_csv(DATA_OUT / "supp8_officialR_new29_hspc_auc.csv")

    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.95), facecolor="white")

    # Row 1: a (% vs K=4000 cell count severity), b (% vs Real hospital CBC cells/L severity), c (Last-draw AUC: % vs cell count)
    ax_a = fig.add_axes([0.085, 0.725, 0.245, 0.215])
    ax_a.scatter(sens["rho_rawPct_276"], sens["rho_cellCount_276"], s=16, color=PALETTE["blue_main"], alpha=0.8, edgecolor="k", lw=0.35)
    ax_a.plot([-0.5, 0.5], [-0.5, 0.5], color="k", lw=0.65, ls="--")
    r_ac = np.corrcoef(sens["rho_rawPct_276"], sens["rho_cellCount_276"])[0, 1]
    ax_a.text(0.05, 0.94, f"r = {r_ac:.3f}\n(43 types)", transform=ax_a.transAxes, va="top", fontsize=5.4, fontweight="bold")
    ax_a.set_xlabel("Same-day WHO weighted ρ (% PBMC)", fontsize=5.5); ax_a.set_ylabel("Same-day WHO weighted ρ (count)", fontsize=5.5)
    panel_title(ax_a, "% vs cell count (276 draws)", size=6.1); fig_panel_label(fig, 0.010, 0.975, "a")

    ax_b = fig.add_axes([0.415, 0.725, 0.245, 0.215])
    ax_b.scatter(sens["rho_rawPct_276"], sens["rho_CBCabs_134"], s=16, color=PALETTE["teal"], alpha=0.85, edgecolor="k", lw=0.35)
    ax_b.axhline(0, color=PALETTE["neutral"], lw=0.5); ax_b.axvline(0, color=PALETTE["neutral"], lw=0.5)
    r_bc = np.corrcoef(sens["rho_rawPct_276"], sens["rho_CBCabs_134"])[0, 1]
    ax_b.text(0.05, 0.94, f"r = {r_bc:.3f}\n(42 diff-CBC draws)", transform=ax_b.transAxes, va="top", fontsize=5.4, fontweight="bold")
    ax_b.set_xlabel("Same-day WHO weighted ρ (% PBMC)", fontsize=5.5); ax_b.set_ylabel("CBC abs. count weighted ρ (cells/L)", fontsize=5.5)
    panel_title(ax_b, "% vs hospital CBC (cells/L)", size=6.1); fig_panel_label(fig, 0.345, 0.975, "b")

    ax_c = fig.add_axes([0.740, 0.725, 0.235, 0.215])
    ax_c.scatter(sens["auc_rawPct_last66"], sens["auc_cellCount_last66"], s=16, color=PALETTE["red_strong"], alpha=0.82, edgecolor="k", lw=0.35)
    ax_c.plot([0.15, 0.85], [0.15, 0.85], color="k", lw=0.65, ls="--")
    r_cc = np.corrcoef(sens["auc_rawPct_last66"], sens["auc_cellCount_last66"])[0, 1]
    ax_c.text(0.05, 0.94, f"r = {r_cc:.3f}\n(66 last draws)", transform=ax_c.transAxes, va="top", fontsize=5.4, fontweight="bold")
    ax_c.set_xlabel("Last-draw AUC (% PBMC)", fontsize=5.5); ax_c.set_ylabel("Last-draw AUC (cell count)", fontsize=5.5)
    panel_title(ax_c, "Mortality: % vs count", size=6.1); fig_panel_label(fig, 0.675, 0.975, "c")

    # Row 2: d (5-specification consistency matrix across top 14 clinical cell types)
    ax_d = fig.add_axes([0.155, 0.385, 0.805, 0.250])
    top14 = sens.sort_values("q_sev_rawPct_276").head(14).reset_index(drop=True)
    mat_s = np.column_stack([
        top14["rho_rawPct_276"], top14["rho_INT_276"], top14["rho_cellCount_276"],
        top14["rho_CBCabs_134"], top14["rho_withCLL_283"],
        2.0 * (top14["auc_rawPct_last66"] - 0.5), 2.0 * (top14["auc_INT_last66"] - 0.5),
    ])
    im = ax_d.imshow(mat_s.T, cmap="RdBu_r", vmin=-0.65, vmax=0.65, aspect="auto")
    spec_lbls = ["WHO: Raw % (276)", "WHO: Source INT", "WHO: Cell count", "WHO: CBC cells/L", "WHO: +CLL (283)", "Death: 2(AUC−0.5) %", "Death: 2(AUC−0.5) INT"]
    ax_d.set_yticks(np.arange(len(spec_lbls))); ax_d.set_yticklabels(spec_lbls, fontsize=5.3)
    ax_d.set_xticks(np.arange(len(top14))); ax_d.set_xticklabels([sh_ct(x) for x in top14["feature"]], fontsize=5.2, rotation=28, ha="right")
    for r_i in range(mat_s.shape[1]):
        for c_i in range(mat_s.shape[0]):
            val = mat_s[c_i, r_i]
            ax_d.text(c_i, r_i, f"{val:+.2f}", ha="center", va="center", fontsize=5.0,
                      color="white" if abs(val) >= 0.36 else PALETTE["gray_3"])
    panel_title(ax_d, "Effect stability across 7 normalization, count, CBC & CLL specifications", size=6.2)
    fig_panel_label(fig, 0.010, 0.655, "d")

    # Row 3: e & f (Official R 4.3.2 12-method bulk deconvolution benchmark — refreshed Folder 29)
    ax_e = fig.add_axes([0.165, 0.060, 0.335, 0.240])
    d12 = deconv.groupby("method", as_index=False)["passed_concordant_FDR"].sum().sort_values("passed_concordant_FDR", ascending=True)
    ye = np.arange(len(d12))
    cols_e = [SRC_COL["PBISC"] if "PBISC" in str(m) else PALETTE["neutral"] for m in d12["method"]]
    ax_e.barh(ye, d12["passed_concordant_FDR"], color=cols_e, edgecolor="k", lw=0.45, height=0.68)
    for yi, v in zip(ye, d12["passed_concordant_FDR"]):
        ax_e.text(v + 0.18, yi, f"{int(v)}/10", va="center", fontsize=5.1, fontweight="bold" if v == d12["passed_concordant_FDR"].max() else "normal")
    ax_e.set_yticks(ye); ax_e.set_yticklabels(d12["method"], fontsize=5.1)
    ax_e.set_xlim(0, 11.6); ax_e.set_xticks([0, 5, 10])
    ax_e.set_xlabel("Concordant FDR signals (/10)", fontsize=5.5)
    panel_title(ax_e, "Official R 4.3.2 deconv. discovery", size=6.0)
    fig_panel_label(fig, 0.010, 0.315, "e")

    ax_f = fig.add_axes([0.655, 0.060, 0.315, 0.240])
    h12 = hspc_dec.copy()
    h12["hspc_dth_auc_plot"] = h12["hspc_dth_auc"].fillna(0.5)
    h12 = h12.sort_values("hspc_dth_auc_plot", ascending=True)
    yf = np.arange(len(h12))
    cols_f = [SRC_COL["PBISC"] if "PBISC" in str(m) else PALETTE["red_2"] for m in h12["method"]]
    ax_f.barh(yf, h12["hspc_dth_auc_plot"], color=cols_f, edgecolor="k", lw=0.45, height=0.68)
    ax_f.axvline(0.5, color="k", lw=0.65, ls="--")
    ax_f.set_yticks(yf); ax_f.set_yticklabels(h12["method"], fontsize=5.1)
    ax_f.set_xlim(0.0, 1.0); ax_f.set_xticks([0.0, 0.5, 1.0])
    ax_f.set_xlabel("Wave-1 new 29 pts HSPC death AUC", fontsize=5.5)
    panel_title(ax_f, "New-patient HSPC AUC vs R deconv.", size=6.0)
    fig_panel_label(fig, 0.525, 0.315, "f")

    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_43CellType_Sensitivity_and_OfficialR_Deconv")


# =============================================================================
# SUPPLEMENTARY FIGURE 23: EXTERNAL 4-COHORT PER-STUDY FOREST PLOTS & TABLES
# =============================================================================
def build_suppfig_external(snum=23):
    ext_coh = pd.read_csv(DATA_OUT / "fig7_external4_per_cohort_results.csv")
    ext_meta = pd.read_csv(DATA_OUT / "fig7_external4_random_effects_meta.csv")

    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.10), facecolor="white")

    # Panel a: External 4-study sample & endpoint table
    ax_a = fig.add_axes([0.03, 0.74, 0.94, 0.22])
    hdr = ["External COVID-19 Study", "Acute samples", "Donors", "Same-day WHO scale", "Last-draw Surv. / Died", "Exact 1:1 mapped subsets"]
    cw = [22, 13, 10, 23, 16, 16]
    rows_e = []
    scale_map = {
        "Stephenson 2021": "Collection-day status (1–7)",
        "COMBAT 2022": "Sampling source WHO (1–7)",
        "van der Wijst 2021": "Visit NIH ordinal (4–7)",
        "Liu 2021": "Visit clinical severity (4.5–7)",
    }
    for cname in ["Stephenson 2021", "COMBAT 2022", "van der Wijst 2021", "Liu 2021"]:
        g = ext_coh[(ext_coh["cohort"] == cname) & (ext_coh["is_exact_43_subset"] == True)]
        r0 = g.iloc[0]
        rows_e.append([cname, str(int(r0["n_samples_sev"])), str(int(r0["n_donors_sev"])), scale_map[cname],
                       f"{int(r0['n_surv'])} / {int(r0['n_dead'])}", f"{len(g)} / 43"])
    rows_e.append(["Total / Pooled (k ≥ 2 RE)", "300", "231", "Equal donor weight w_ij = 1/n_i", "197 / 22", "23 exact (+B_ASC_total)"])
    _draw_pub_table(ax_a, hdr, rows_e, cw, "External 4 public single-cell COVID-19 studies analyzed by identical Fig. 7 criteria")
    fig_panel_label(fig, 0.008, 0.985, "a")

    # Panels b & c: Complete 23-subset heatmaps (cohort-level + RE meta) for Severity (b) and Last-draw Mortality (c)
    em23 = ext_meta[ext_meta["is_exact_43_subset"] == True].sort_values("rho_re").reset_index(drop=True)
    feats23 = em23["feature"].tolist()
    cohs4 = ["Stephenson 2021", "COMBAT 2022", "van der Wijst 2021", "Liu 2021"]

    mat_sev = np.full((len(feats23), 6), np.nan)
    mat_mrt = np.full((len(feats23), 6), np.nan)
    for fi, f in enumerate(feats23):
        for ci, cn in enumerate(cohs4):
            sub = ext_coh[(ext_coh["cohort"] == cn) & (ext_coh["feature"] == f)]
            if len(sub):
                mat_sev[fi, ci] = sub["rho_w"].iloc[0]
                mat_mrt[fi, ci] = sub["auc_death"].iloc[0]
        mr = em23.iloc[fi]
        mat_sev[fi, 4] = mr["rho_re"]
        mat_sev[fi, 5] = mr["eskd_rho_w"]
        mat_mrt[fi, 4] = mr["auc_re"]
        mat_mrt[fi, 5] = mr["eskd_auc_last66"]

    col_names = ["Steph. 21", "COMBAT 22", "vdWijst 21", "Liu 21", "Ext. RE Meta", "ESKD (276/66)"]
    for idx_p, (mat, vmin, vmax, vmid, ttl, let, xl) in enumerate([
        (mat_sev, -0.55, 0.55, 0.0, "Same-day WHO weighted ρ_w", "b", 0.145),
        (mat_mrt, 0.15, 0.85, 0.5, "Last-draw mortality AUC", "c", 0.625),
    ]):
        ax = fig.add_axes([xl, 0.075, 0.345, 0.615])
        cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#EFEFEF")
        ax.imshow(np.ma.masked_invalid(mat), cmap=cmap, norm=mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=vmid, vmax=vmax), aspect="auto")
        ax.axvline(3.5, color="k", lw=0.85)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.0,
                            fontweight="bold" if j >= 4 else "normal",
                            color="white" if abs(v - vmid) >= 0.24 else PALETTE["gray_3"])
        ax.set_yticks(np.arange(len(feats23)))
        ax.set_yticklabels([sh_ct(x) for x in feats23] if idx_p == 0 else [], fontsize=5.2)
        ax.set_xticks(np.arange(6)); ax.set_xticklabels(col_names, fontsize=5.2, rotation=32, ha="right")
        ax.tick_params(axis="both", length=0, pad=1.5)
        for sp_ in ax.spines.values():
            sp_.set_visible(True); sp_.set_linewidth(0.55); sp_.set_edgecolor(PALETTE["gray_3"])
        panel_title(ax, ttl, size=6.1, pad=3)
        fig_panel_label(fig, 0.008 if idx_p == 0 else 0.515, 0.710, let)

    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_External_4Cohort_SameCriterion_Validation")


# =============================================================================
# SUPPLEMENTARY FIGURE 22: SINGLE-CELL T-ISG & HSPC FUNCTIONAL MODULES
# =============================================================================
def build_suppfig_tisg_hspc(snum=22):
    tisg13 = pd.read_csv(DATA_OUT / "supp10_tcell_isg_13subsets_summary.csv")
    cd8_fl = pd.read_csv(DATA_OUT / "supp10_cd8ctl_isg_sameday_first_last.csv")
    hmod = pd.read_csv(DATA_OUT / "supp10_hspc_16modules_summary.csv")
    hgene = pd.read_csv(DATA_OUT / "supp10_hspc_20genes_summary.csv")

    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.60), facecolor="white")

    # Panel a: T CD8 CTL 7-ISG temporal decoupling (same-day WHO vs First draw vs Last draw mortality)
    ax_a = fig.add_axes([0.10, 0.60, 0.36, 0.31])
    n0_s, n0_p = int(cd8_fl.iloc[0]["n_samples"]), int(cd8_fl.iloc[0]["n_patients"])
    n1_p = int(cd8_fl.iloc[1]["n_patients"])
    n2_p = int(cd8_fl.iloc[2]["n_patients"])
    lbls_a = [f"Same-day WHO\n(weighted ρ, {n0_s} draws)", f"First-draw death\n(AUC−0.5, n={n1_p})", f"Last-draw death\n(AUC−0.5, n={n2_p})"]
    vals_a = [
        cd8_fl.iloc[0]["estimate"],
        cd8_fl.iloc[1]["estimate"] - 0.5,
        cd8_fl.iloc[2]["estimate"] - 0.5,
    ]
    lo_a = [
        vals_a[0] - cd8_fl.iloc[0]["ci_lo"],
        cd8_fl.iloc[1]["estimate"] - cd8_fl.iloc[1]["ci_lo"],
        cd8_fl.iloc[2]["estimate"] - cd8_fl.iloc[2]["ci_lo"],
    ]
    hi_a = [
        cd8_fl.iloc[0]["ci_hi"] - vals_a[0],
        cd8_fl.iloc[1]["ci_hi"] - cd8_fl.iloc[1]["estimate"],
        cd8_fl.iloc[2]["ci_hi"] - cd8_fl.iloc[2]["estimate"],
    ]
    xa = np.arange(3)
    cols_a = [PALETTE["neutral"], PALETTE["blue_main"], PALETTE["red_2"]]
    ax_a.bar(xa, vals_a, 0.52, color=cols_a, edgecolor="k", lw=0.6)
    ax_a.errorbar(xa, vals_a, yerr=[np.maximum(0, lo_a), np.maximum(0, hi_a)], fmt="none", ecolor="k", elinewidth=0.75, capsize=2.0)
    ax_a.axhline(0.0, color="k", lw=0.65)
    for xi, (v, hi_v, pr) in enumerate(zip(vals_a, hi_a, cd8_fl["p"])):
        txt_v = f"ρ={v:+.2f}\nP={pr:.3g}" if xi == 0 else f"AUC={v+0.5:.2f}\nP={pr:.3g}"
        ax_a.text(xi, max(v + max(0, hi_v), 0) + 0.02, txt_v, ha="center", va="bottom", fontsize=5.2, fontweight="bold")
    ax_a.set_xticks(xa); ax_a.set_xticklabels(lbls_a, fontsize=5.3)
    ax_a.set_ylim(-0.22, 0.58); ax_a.set_ylabel("Effect size (centered at 0)", fontsize=5.6)
    panel_title(ax_a, "CD8 CTL 7-ISG: early vs late draws", size=6.1)
    fig_panel_label(fig, 0.010, 0.960, "a")

    # Panel b: T-cell 7-ISG across 13 T-cell subsets (Shared-27 fidelity & clinical decoupling)
    ax_b = fig.add_axes([0.58, 0.60, 0.39, 0.31])
    s_col = "subset" if "subset" in tisg13.columns else tisg13.columns[0]
    r_col = [c for c in tisg13.columns if "sh_rho" in c or "shared" in c.lower()][0]
    yb = np.arange(len(tisg13))
    ax_b.barh(yb, tisg13[r_col], color=PALETTE["blue_secondary"], edgecolor="k", lw=0.45, height=0.68)
    ax_b.axvline(0.0, color="k", lw=0.6)
    ax_b.set_yticks(yb); ax_b.set_yticklabels(tisg13[s_col], fontsize=5.1); ax_b.invert_yaxis()
    ax_b.set_xlim(-0.1, 1.0); ax_b.set_xlabel("Shared-27 Real vs PBISC Spearman ρ", fontsize=5.5)
    panel_title(ax_b, "13-subset T-cell ISG fidelity", size=6.1)
    fig_panel_label(fig, 0.500, 0.960, "b")

    # Panel c: HSPC 16 intrinsic modules (Combined n=62)
    ax_c = fig.add_axes([0.18, 0.08, 0.29, 0.41])
    hm_s = hmod[hmod["cohort"] == "Combined"].sort_values("AUC_dth11", ascending=True).reset_index(drop=True)
    yc = np.arange(len(hm_s))
    ax_c.barh(yc, hm_s["AUC_dth11"], color=PALETTE["teal"], edgecolor="k", lw=0.45, height=0.68)
    ax_c.axvline(0.5, color="k", lw=0.65, ls="--")
    ax_c.set_yticks(yc); ax_c.set_yticklabels(hm_s["module"], fontsize=5.1)
    ax_c.set_xlim(0.2, 0.90); ax_c.set_xlabel("Mortality AUC (62 HSPC+ pts)", fontsize=5.5)
    panel_title(ax_c, "HSPC 16 functional modules", size=6.1)
    fig_panel_label(fig, 0.010, 0.520, "c")

    # Panel d: HSPC 20 canonical marker/program genes (Combined n=62)
    ax_d = fig.add_axes([0.62, 0.08, 0.35, 0.41])
    hg_s = hgene[hgene["cohort"] == "Combined"].sort_values("AUC_dth11", ascending=True).reset_index(drop=True)
    yd = np.arange(len(hg_s))
    ax_d.barh(yd, hg_s["AUC_dth11"], color=PALETTE["violet"], edgecolor="k", lw=0.45, height=0.68)
    ax_d.axvline(0.5, color="k", lw=0.65, ls="--")
    ax_d.set_yticks(yd); ax_d.set_yticklabels(hg_s["gene"], fontsize=5.1, fontstyle="italic")
    ax_d.set_xlim(0.2, 0.90); ax_d.set_xlabel("Mortality AUC (62 HSPC+ pts)", fontsize=5.5)
    panel_title(ax_d, "HSPC 20 canonical genes", size=6.1)
    fig_panel_label(fig, 0.500, 0.520, "d")

    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_SingleCell_TISG_and_HSPC_Programs")


def main():
    apply_publication_style(FigureStyle())
    supp_offset = 13
    cli_args = []
    for a in sys.argv[1:]:
        if a.startswith("--supp-offset="):
            supp_offset = int(a.split("=", 1)[1])
        else:
            cli_args.append(a)

    s_meta = supp_offset + 6
    s_swim = supp_offset + 7
    s_sens = supp_offset + 8
    s_prog = supp_offset + 9
    s_ext = supp_offset + 10

    builders = {
        "fig7": build_figure7,
        "sfig19": lambda: build_suppfig_metadata(s_meta),
        "sfig20": lambda: build_suppfig_swimmer(s_swim),
        "sfig21": lambda: build_suppfig_sensitivity(s_sens),
        "sfig22": lambda: build_suppfig_tisg_hspc(s_prog),
        "sfig23": lambda: build_suppfig_external(s_ext),
    }
    wanted = [a for a in cli_args if a in builders] or list(builders)
    audit_path = OUT_ROOT / f"AUDIT_FIG7_SFIG{s_meta}_{s_ext}.json"
    prev = {a["figure"]: a for a in json.loads(audit_path.read_text(encoding="utf-8"))} if audit_path.exists() else {}
    for key in wanted:
        _, aud = builders[key]()
        prev[aud["figure"]] = aud
        print(f"[OK] {aud['figure']}: {aud['width_mm']}x{aud['height_mm']} mm | min_pt={aud['min_font_pt']} | oob={aud['out_of_bounds_count']}", flush=True)

    order = [
        "Figure7_SameDay_Severity_and_LastDraw_Mortality",
        f"SuppFig{s_meta}_Cohort_and_Sample_Metadata_Tables",
        f"SuppFig{s_swim}_Patient_Swimmer_276_Acute_Draws",
        f"SuppFig{s_sens}_43CellType_Sensitivity_and_OfficialR_Deconv",
        f"SuppFig{s_prog}_SingleCell_TISG_and_HSPC_Programs",
        f"SuppFig{s_ext}_External_4Cohort_SameCriterion_Validation",
    ]
    audit_path.write_text(json.dumps([prev[k] for k in order if k in prev], indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
