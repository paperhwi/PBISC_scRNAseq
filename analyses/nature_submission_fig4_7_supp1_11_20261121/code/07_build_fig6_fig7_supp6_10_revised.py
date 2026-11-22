"""Build revised Nature Portfolio 180-mm figures honoring the user's exact directives:
  - Main Figure 6 (Figure6_Rare_Population_Sepsis_COVID_Reproduction):
      Keeps Sepsis rare-cell data (14 rare shares, 82-donor proliferating-monocyte audit,
      Sepsis rare-cell Harmony UMAP, Sepsis 2,000-HVG specificity & sample ceiling, Sepsis
      Mono_c06-CD14-MKI67 abundance + canonical marker heatmap) TOGETHER WITH COVID-19 rare-cell
      reproduction (13 rare shares, COVID rare-cell Harmony UMAP, COVID 2,000-HVG specificity,
      Shared-27 bridge calibration, Official R 4.3.2 11-method out-of-donor deconvolution
      benchmark, 37-43x HSPC sink inflation, and 3,124 single-cell HSPC 20-gene lineage switch).
  - Main Figure 7 (Figure7_COVID_Clinical_Metadata_CountPct_TISG):
      Focused squarely on the user's 3 priorities:
        (1) Exact patient & dataset metadata matching (107 donors / 362 samples 1:1 audit)
        (2) Population number (raw count & clinical CBC abs) vs percentage (% PBMC & % parent)
            severity (217 draws x same-day WHO GEE +wk) & mortality (66 first draws) check
        (3) Deep T-cell ISG scoring analysis (Shared-27 vs BayesPrism Z, time-bias removal,
            Sepsis T-ISG suppression contrast, first-draw T CD8 CTL abundance vs ISG
            anti-correlation rho=-0.396, triphasic Time x Fatal collapse, and leak-free 50x5 CV).
  - Upgraded SuppFig 6 (exact Sepsis 93/8/7 + COVID 362 1:1 matching metadata tables),
    SuppFig 9 (Population number vs % & Clinical CBC audit + Sepsis 44-subset check + 14-var matrix),
    SuppFig 10 (13-subset T-cell ISG audit + Early/Late + Cross-disease 19-ISG DEG + HSPC/Stephenson/Power).
"""
import importlib.util
import json
import os
import pathlib
import shutil
import sys
import tempfile
import warnings
import zipfile
from xml.sax.saxutils import escape

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
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from scipy import stats
import statsmodels.api as sm

from nature_style import (
    NATURE_2COL_IN,
    PALETTE,
    SRC_COL,
    FigureStyle,
    apply_publication_style,
    fig_panel_label,
    finalize_figure,
    panel_title,
)

# Import verified 1x rare-cell drawing helpers from 02_build_fig4_fig5_supp1_5.py
spec02 = importlib.util.spec_from_file_location("b02", SCRIPT_DIR / "02_build_fig4_fig5_supp1_5.py")
b02 = importlib.util.module_from_spec(spec02)
spec02.loader.exec_module(b02)

ROOT = pathlib.Path(os.environ.get("PBISC_ROOT", r"F:\PBISC"))
D28_DATA = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "28_Verified_True11_Deconv_HSC_TISG_20261119" / "data"
D30_TAB = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "30_DrawLevel_Reanalysis_and_Patient_Draw_Table_20261120" / "tables"
SEP_NK_DATA = ROOT / "sepsis_nk_cd56bright_20261115" / "data"
P_DIR = ROOT / "deliverables" / "manuscript_real_pbisc_scsemi_20261110" / "panels"
OUT_ROOT = ROOT / "deliverables" / "nature_submission_fig4_6_20261120"
NEW_DIR = OUT_ROOT / "draw_sameday_sev_firstdraw_death_20261120"
NEW_TAB = NEW_DIR / "tables"
NEW_FIG = NEW_DIR / "figures"
AUD3 = OUT_ROOT / "audit_focus3_20261121"
FIG_OUT = OUT_ROOT / "figures"
DATA_OUT = OUT_ROOT / "data"
FIG_OUT.mkdir(parents=True, exist_ok=True)

CLL = {"C116", "C139"}
COL = dict(
    Real=PALETTE["gray_3"],
    New=PALETTE["blue_secondary"],
    Ext=PALETTE["blue_main"],
    Shared=PALETTE["teal"],
    Shuf=PALETTE["neutral"],
    All=PALETTE["red_2"],
    CLL=PALETTE["red_strong"],
    HSPC=PALETTE["red_strong"],
    Teal=PALETTE["teal"],
    Violet=PALETTE["violet"],
    Good="#1B7A43",
    Mid="#B86B11",
    Bad="#B64342",
    CardBlue="#DCE6F2",
    CardGray="#EFEFEF",
    CardTeal="#DCF0F2",
    CardRed="#FADBD8",
)
MK = dict(Real="D", New="o", Ext="s")

LBL = {
    "B_ASC_dividing": "B ASC dividing",
    "NK CD56 bright": "NK CD56bright",
    "HSPC": "HSPC",
    "T CD8 CM": "T CD8 CM",
    "T CD8 CTL": "T CD8 CTL",
    "B_ASC_IgG": "B ASC IgG",
    "B_ASC_IgA": "B ASC IgA",
    "B_ASC_IgM": "B ASC IgM",
    "CD14mono": "CD14 mono",
    "CD14mono_IFN": "CD14 mono IFN",
    "CD14mono_dex": "CD14 mono dex",
    "CD16mono": "CD16 mono",
    "CD16mono_C1": "CD16 mono C1",
    "Int.mono": "Int. mono",
    "B_CD11c": "B CD11c",
    "DC1": "DC1",
    "DC2": "DC2",
    "DC3": "DC3",
    "NK": "NK (CD56dim)",
    "NKT": "NKT",
    "T CD4 CTL": "T CD4 CTL",
    "T CD4 cycling": "T CD4 cycling",
    "T g/d": "T γδ",
    "MAIT": "MAIT",
    "ASDC": "ASDC",
}


def fp_s(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"P={p:.1e}".replace("e-0", "e-")
    if p < 0.001:
        return f"P={p:.4f}"
    if p < 0.01:
        return f"P={p:.3f}"
    return f"P={p:.2f}"


def fp_short(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return "<1e-4"
    if p < 0.001:
        return f"{p:.1e}".replace("e-0", "e-")
    if p < 0.01:
        return f"{p:.3f}".lstrip("0")
    return f"{p:.2f}".lstrip("0")


def auc_mw(y, f):
    y = np.asarray(y, float)
    f = np.asarray(f, int)
    m = np.isfinite(y)
    y, f = y[m], f[m]
    a, b = y[f == 1], y[f == 0]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u / (len(a) * len(b))), float(p)


def rbox(ax, x, y, w, h, fc, ec, lw=0.85, r=1.8, pad=0.18, alpha=1.0, ls="-", zorder=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={pad},rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, linestyle=ls,
                                alpha=alpha, zorder=zorder))


def plot_ci_clipped(ax, y, lo, hi, color, xlim, lw=1.0, ls="-", alpha=1.0, zorder=3):
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return
    c_lo, c_hi = max(lo, xlim[0]), min(hi, xlim[1])
    if c_hi >= c_lo:
        ax.plot([c_lo, c_hi], [y, y], color=color, lw=lw, ls=ls, alpha=alpha, zorder=zorder)
    if lo < xlim[0]:
        ax.scatter(xlim[0] + 0.01 * (xlim[1] - xlim[0]), y, marker="<", s=10, c=color, zorder=zorder + 1)
    if hi > xlim[1]:
        ax.scatter(xlim[1] - 0.01 * (xlim[1] - xlim[0]), y, marker=">", s=10, c=color, zorder=zorder + 1)


def load_data():
    d = {}
    d["st_df"] = b02.get_rare_bootstrap_10states()
    d["rcache"] = b02.get_raresim_cache()
    d["T1_PAT68"] = pd.read_csv(NEW_TAB / "T1_all68_patients_draw_schedule.csv")
    d["T1C_SUM"] = pd.read_csv(NEW_TAB / "T1c_cohort_draw_summary.csv")
    d["T1D_NEG39"] = pd.read_csv(NEW_TAB / "T1d_covid_negative_39_donors.csv")
    d["T2_ALL"] = pd.read_csv(NEW_TAB / "T2_comp43_patientmean_vs_samedaySev_firstdrawDeath.csv")
    d["T2B"] = pd.read_csv(NEW_TAB / "T2b_comp43_combined_full_comparison.csv")
    d["T3_FUN"] = pd.read_csv(NEW_TAB / "T3_new_framework_43_funnel.csv")
    d["T4_PW"] = pd.read_csv(NEW_TAB / "T4_power_rescue_new_framework.csv")
    d["T5_HSPC"] = pd.read_csv(NEW_TAB / "T5_HSPC_16modules_20genes_new_framework.csv")
    d["T6_ISG"] = pd.read_csv(NEW_TAB / "T6_Tcell_ISG_13subsets_new_framework.csv")
    d["T7_LABS"] = pd.read_csv(NEW_TAB / "T7_hospital_labs_sameday_pairs.csv")
    d["T8B_WIN"] = pd.read_csv(D30_TAB / "T8b_Tcell_ISG_13subsets_time_windows.csv")
    d["P217_Z"] = pd.read_csv(DATA_OUT / "primary_217draws_rankINT_43subsets.csv")
    d["P66_FD"] = pd.read_csv(DATA_OUT / "primary_66firstdraws_rankINT_43subsets.csv")
    d["H149_Z"] = pd.read_csv(DATA_OUT / "primary_149draws_HSPC_modules.csv")
    d["H61_FD"] = pd.read_csv(DATA_OUT / "primary_61firstdraws_HSPC_modules.csv")
    d["S_meta"] = pd.read_csv(D28_DATA / "eskd_sample_metadata_verified.csv")
    d["pos_217"] = pd.read_csv(D28_DATA / "eskd_longitudinal_217_draws.csv")
    d["MOD"] = pd.read_csv(D28_DATA / "suppS2_modality_scale_summary.csv")
    d["SH27"] = pd.read_csv(D28_DATA / "suppS2_shared27_paired_pct.csv")
    d["PR_REALB"] = pd.read_csv(D28_DATA / "suppS2_realB_transfer_precision_recall.csv").set_index("label")
    d["SA"] = pd.read_csv(D28_DATA / "suppS3_sensitivity_14variants_true11.csv")
    d["NKD"] = pd.read_csv(D28_DATA / "suppS2_nk_decomposition.csv")
    d["PW_UNI"] = pd.read_csv(D28_DATA / "fig7_s3_unified_redetection_true11.csv")
    d["PUB_COH"] = pd.read_csv(D28_DATA / "suppS4_public_cohort_reconciled.csv")
    d["PUB_META"] = pd.read_csv(D28_DATA / "suppS4_public_meta_reconciled.csv")
    d["PUB_HW"] = pd.read_csv(D28_DATA / "suppS4_public_hspc_windows_reconciled.csv")
    d["ESKD_MATCH"] = pd.read_csv(D28_DATA / "suppS4_eskd_matched_severity_ci.csv").set_index("feature")
    d["SG_HSPC"] = pd.read_csv(D28_DATA / "suppS5_single_gene_hspc_new29.csv")
    d["BENCH10_D28"] = pd.read_csv(D28_DATA / "fig8_deconv_10methods_true11.csv")
    d["TISG13"] = pd.read_csv(D28_DATA / "fig8_tcell_isg_13subsets_true11.csv")
    d["TISG_GEE"] = pd.read_csv(D28_DATA / "suppS5_tcell_isg_gee_true11.csv")
    d["CD8_ISG_217"] = pd.read_csv(D28_DATA / "fig8_cd8ctl_isg_longitudinal_draws.csv")
    d["P43_WIN"] = pd.read_csv(D28_DATA / "suppS5_paired43_window_AUC_true11.csv")
    d["CV_500K"] = pd.read_csv(D28_DATA / "fig8_clinical_CV_infold_INT_500k_true11.csv")
    d["OFF_R_SH27"] = pd.read_csv(DATA_OUT / "official_R432_shared27_benchmark.csv")
    d["VER_JSON"] = json.loads((OUT_ROOT / "VERIFICATION_REPORT_20261120.json").read_text(encoding="utf-8"))
    # Audit 3 tables (newly generated!)
    d["A1_SEP"] = pd.read_csv(AUD3 / "A1_sepsis_255_samples_matched_metadata.csv")
    d["A2_COV"] = pd.read_csv(AUD3 / "A2_covid_43subsets_count_vs_pct_vs_parent_vs_CBC.csv")
    d["A2_SEP"] = pd.read_csv(AUD3 / "A2b_sepsis_count_vs_pct_HC_ICU_Sepsis.csv")
    d["A3_CV"] = pd.read_csv(AUD3 / "A3_covid_tisg_cv_benchmark.csv")
    d["A3_SEP_ISG"] = pd.read_csv(AUD3 / "A3b_sepsis_pbisc_tcell_isg.csv")
    d["A3_CROSS_DEG"] = pd.read_csv(AUD3 / "A3c_cross_disease_Tcell_19ISG_deg.csv")
    d["A3_FD66"] = pd.read_csv(AUD3 / "A3d_covid_66firstdraw_CTL_NK_ISG.csv")
    d["A3_JSON"] = json.loads((AUD3 / "AUDIT_FOCUS3_SUMMARY.json").read_text(encoding="utf-8"))
    d["SEP_REAL_NK"] = pd.read_csv(SEP_NK_DATA / "sepsis_real_sample_table.csv")
    d["SEP_PB_NK"] = pd.read_csv(SEP_NK_DATA / "sepsis_pbisc_sample_table.csv")
    return d


def draw_sameday_strip(ax, d, f, rng_local, show_yl=True, title_txt=None, title_col=None, compact=True):
    P217, T2 = d["P217_Z"], d["T2_ALL"]
    zc = f + "__z217"
    for g in range(4):
        for arm_v, k, dx in [("Real", "Real", -0.17), ("PBISC new", "New", 0.17)]:
            vals = P217[(P217["arm"] == arm_v) & (P217["temp_severity_ord"] == g)][zc].dropna().values
            ax.scatter(g + dx + rng_local.uniform(-0.08, 0.08, len(vals)), vals, s=5.0, marker=MK[k],
                       c=COL[k], alpha=0.72, edgecolors="white", linewidths=0.18, zorder=3, rasterized=True)
            if len(vals):
                ax.plot([g + dx - 0.13, g + dx + 0.13], [np.median(vals)] * 2, color="black", lw=0.95, zorder=4)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=5.1)
    ax.set_xlim(-0.55, 3.55); ax.set_ylim(-2.6, 5.25); ax.set_yticks([-2, 0, 2])
    if show_yl:
        ax.set_ylabel("Draw normal score", fontsize=5.3)
    else:
        ax.set_yticklabels([])
    for idx_s, (sname, k, disp_s) in enumerate([("Real", "Real", "R"), ("PBISC new", "New", "N"), ("Combined", "Ext", "C")]):
        r = T2[(T2["feature"] == f) & (T2["set"] == sname)].iloc[0]
        if compact:
            if sname == "Combined":
                txt = f"C:{r.sd_gee_wk_beta:+.2f}, q={fp_short(r.sd_gee_wk_p_q43)}".replace("q=<", "q<")
            else:
                txt = f"{disp_s}:{r.sd_gee_wk_beta:+.2f}, P={fp_short(r.sd_gee_wk_p)}".replace("P=<", "P<")
        else:
            q_s = f", q={fp_short(r.sd_gee_wk_p_q43)}" if sname == "Combined" else ""
            txt = f"{disp_s}: β={r.sd_gee_wk_beta:+.2f}, P={fp_short(r.sd_gee_wk_p)}{q_s}".replace("P=<", "P<").replace("q=<", "q<")
        ax.text(0.02, 0.975 - 0.108 * idx_s, txt, transform=ax.transAxes, fontsize=5.0,
                color=COL[k], va="top", fontweight="bold" if r.sd_gee_wk_p < 0.05 else "normal")
    if title_txt:
        panel_title(ax, title_txt, size=5.8, pad=3, color=title_col or COL["Real"])


def draw_firstdraw_death_strip(ax, d, f, rng_local, show_yl=False, title_txt=None, title_col=None, compact=True):
    P66, T2 = d["P66_FD"], d["T2_ALL"]
    zc = f + "__zfd"
    for g in range(2):
        for arm_v, k, dx in [("Real", "Real", -0.17), ("PBISC new", "New", 0.17)]:
            vals = P66[(P66["arm"] == arm_v) & (P66["fatal_pub11"] == g)][zc].dropna().values
            ax.scatter(g + dx + rng_local.uniform(-0.07, 0.07, len(vals)), vals, s=8.5, marker=MK[k],
                       c=COL[k], edgecolors="white", linewidths=0.22, zorder=3)
            if len(vals):
                ax.plot([g + dx - 0.13, g + dx + 0.13], [np.median(vals)] * 2, color="black", lw=0.95, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Surv.\n(55)", "Died\n(11)"], fontsize=5.1)
    ax.set_xlim(-0.52, 1.52); ax.set_ylim(-2.6, 5.25); ax.set_yticks([-2, 0, 2])
    if show_yl:
        ax.set_ylabel("1st-draw score", fontsize=5.3)
    else:
        ax.set_yticklabels([])
    for idx_s, (sname, k, disp_s) in enumerate([("Real", "Real", "R"), ("PBISC new", "New", "N"), ("Combined", "Ext", "C")]):
        r = T2[(T2["feature"] == f) & (T2["set"] == sname)].iloc[0]
        if compact:
            if sname == "Combined":
                tag = f" aP={fp_short(r.fd_adjWHO_p)}" if r.fd_adjWHO_p < 0.05 else f" q={fp_short(r.fd_q43)}"
                txt = f"C:{r.fd_AUC:.2f} ({fp_short(r.fd_p)},{tag})"
            else:
                txt = f"{disp_s}:{r.fd_AUC:.2f} (P={fp_short(r.fd_p)})"
        else:
            tag = (f", adjP={fp_short(r.fd_adjWHO_p)}" if r.fd_adjWHO_p < 0.05 else f", q={fp_short(r.fd_q43)}") if sname == "Combined" else ""
            txt = f"{disp_s}: {r.fd_AUC:.2f}, P={fp_short(r.fd_p)}{tag}".replace("P=<", "P<")
        ax.text(0.02, 0.975 - 0.108 * idx_s, txt, transform=ax.transAxes, fontsize=5.0,
                color=COL[k], va="top", fontweight="bold" if r.fd_p < 0.05 else "normal")
    if title_txt:
        panel_title(ax, title_txt, size=5.8, pad=3, color=title_col or COL["Real"])


# =============================================================================
# MAIN FIGURE 6 (a - l): DISEASE-SPECIFIC RARE CELL POPULATION REPRODUCTION
#                        ACROSS BACTERIAL SEPSIS & ESKD COVID-19 (WITH UMAPs!)
# =============================================================================
def build_figure6(d):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.42), facecolor="white")
    st_df, rcache = d["st_df"], d["rcache"]

    # ---- Row 1 (y=0.818..0.962): a (Sepsis 14 rare shares), b (COVID 13 rare shares), c (Sepsis 82-donor prolif mono audit) ----
    ax_a = fig.add_axes([0.118, 0.815, 0.172, 0.145])
    b02.draw_rare_parent_1x(ax_a, "sepsis", show_legend=True)
    fig_panel_label(fig, 0.008, 0.990, "a")

    ax_b = fig.add_axes([0.398, 0.815, 0.172, 0.145])
    b02.draw_rare_parent_1x(ax_b, "covid", show_legend=False)
    fig_panel_label(fig, 0.305, 0.990, "b")

    # c: Compact Observed Sepsis proliferating-monocyte (CD14/MKI67/STMN1/MPO/ELANE) heatmap (82 paired donors)
    b02.draw_prolif_mono_1x(fig, [0.592, 0.815, 0.398, 0.155])
    fig_panel_label(fig, 0.585, 0.990, "c")

    # ---- Row 2 (y=0.630..0.768): d (Sepsis rare-cell Harmony UMAP) & e (COVID-19 rare-cell Harmony UMAP) ----
    ax_d = fig.add_axes([0.052, 0.632, 0.390, 0.132])
    b02.draw_umap_rare_1x(ax_d, "sepsis")
    fig_panel_label(fig, 0.008, 0.782, "d")

    ax_e = fig.add_axes([0.548, 0.632, 0.390, 0.132])
    b02.draw_umap_rare_1x(ax_e, "covid")
    fig_panel_label(fig, 0.502, 0.782, "e")

    # ---- Row 3 (y=0.408..0.580): f (Sepsis 2,000-HVG specificity + ceiling) & g (COVID 2,000-HVG specificity + ceiling) ----
    b02.draw_hvg_ceiling_1x(fig, [0.008, 0.408, 0.484, 0.175], rcache["sepsis"], "f")
    b02.draw_hvg_ceiling_1x(fig, [0.504, 0.408, 0.484, 0.175], rcache["covid"], "g")

    # ---- Row 4 (y=0.232..0.388): h (Sepsis Mono_c06-MKI67 abundance + markers) & i (COVID Shared-27 bridge fidelity & scale) ----
    ax_h1 = fig.add_axes([0.062, 0.242, 0.116, 0.130])
    b02.draw_rare_patient_1x(ax_h1, "sepsis", "Mono_c06-CD14-MKI67", st_df, short_title="Sepsis · Mono_MKI67")
    fig_panel_label(fig, 0.008, 0.395, "h")

    b02.draw_pop_marker_heatmap_1x(
        fig, [0.195, 0.228, 0.295, 0.162], rcache["sepsis"],
        "Mono_c06-CD14-MKI67", "Mono_c06-MKI67 markers", ["Monocyte", "Cell cycle"], letter=None
    )

    # i: Paired Shared-27 bridge fidelity & scale calibration (COVID-19 n=27)
    ax_i = fig.add_axes([0.562, 0.258, 0.418, 0.122]); fig_panel_label(fig, 0.502, 0.395, "i")
    FEATS_M = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
    X_SH27_2L = ["NK\n56br", "B ASC\ndiv.", "B ASC\nIgG", "B ASC\nIgA", "T CD8\nCM", "T CD8\nCTL", "HSPC"]
    SH27, MOD, PR_REALB = d["SH27"], d["MOD"], d["PR_REALB"]
    rng_m = np.random.default_rng(19)
    for ii, f in enumerate(FEATS_M):
        vr, vb = SH27[f + "|Real_pct"].values, SH27[f + "|PBISC_pct"].values
        for r_i, b_i in zip(vr, vb):
            ax_i.plot([ii - 0.18, ii + 0.18], [np.log10(r_i + 0.02), np.log10(b_i + 0.02)], color=PALETTE["neutral"], lw=0.4, alpha=0.55)
        ax_i.scatter(ii - 0.18 + rng_m.uniform(-0.04, 0.04, len(vr)), np.log10(vr + 0.02), s=7.0, marker="D", c=COL["Real"], edgecolors="white", lw=0.2, zorder=3)
        ax_i.scatter(ii + 0.18 + rng_m.uniform(-0.04, 0.04, len(vb)), np.log10(vb + 0.02), s=8.0, marker="o", c=COL["New"], edgecolors="white", lw=0.2, zorder=3)
        mrow = MOD[MOD["feature"] == f].iloc[0]
        ax_i.text(ii, 1.84, f"{mrow['fold_PBISC_over_Real']:.1f}×", ha="center", va="top", fontsize=5.1, fontweight="bold", color=COL["Ext"])
    yt_p = [0.0, 0.1, 0.5, 2.0, 10.0, 40.0]
    ax_i.set_yticks([np.log10(v + 0.02) for v in yt_p]); ax_i.set_yticklabels(["0%", "0.1%", "0.5%", "2%", "10%", "40%"], fontsize=5.1)
    ax_i.set_ylim(np.log10(0.008), 1.98); ax_i.set_xticks(range(len(FEATS_M)))
    ax_i.set_xticklabels(X_SH27_2L, fontsize=5.0, fontweight="bold")
    ax_i.set_ylabel("% PBMC (Shared-27)", fontsize=5.3)
    h_pr = PR_REALB.loc["HSPC"]
    ax_i.text(0.99, 0.03, f"Real B 30-NN prec/rec: HSPC {100*h_pr['precision']:.1f}%/{100*h_pr['recall']:.1f}%, B ASC 98.4%/98.6%",
              transform=ax_i.transAxes, ha="right", va="bottom", fontsize=5.0, color=PALETTE["gray_2"],
              bbox=dict(fc="white", ec="none", alpha=0.92, pad=0.2))
    panel_title(ax_i, "COVID-19 · Shared-27 bridge fidelity & scale calibration (n=27)", size=6.0, pad=3)

    # ---- Row 5 (y=0.036..0.195): j (Official R 4.3.2 11-method Deconv), k (37-43x HSPC sink), l (Single-cell HSPC 20-gene switch) ----
    off_r = d["OFF_R_SH27"].set_index("method")
    ax_j = fig.add_axes([0.118, 0.034, 0.285, 0.150]); fig_panel_label(fig, 0.008, 0.208, "j")
    M_ORD = ["PBISC", "BayesPrism-init", "BayesPrism-upd", "BisqueRNA", "CIBERSORTx-Bmode",
             "CIBERSORTx-Smode", "CIBERSORT-R", "Calib-NNLS", "MuSiC", "NNLS", "DWLS"]
    C_ORD = [
        ("med_rho_all43", "All43"), ("med_rho_rare23", "Rare23"),
        ("T_CD8_CTL_rho", "CTL"), ("B_ASC_IgG_rho", "IgG"),
        ("B_ASC_IgA_rho", "IgA"), ("HSPC_rho", "HSPC"),
    ]
    mat = off_r.loc[M_ORD, [c[0] for c in C_ORD]].values.astype(float)
    masked = np.ma.masked_invalid(mat)
    cmap_u = mpl.colormaps["YlGnBu"].copy(); cmap_u.set_bad("#EAEAEA")
    ax_j.imshow(masked, cmap=cmap_u, norm=mpl.colors.Normalize(-0.1, 0.7), aspect="auto")
    for ii in range(mat.shape[0]):
        for jj in range(mat.shape[1]):
            v = mat[ii, jj]
            if np.isnan(v):
                ax_j.text(jj, ii, "0%", ha="center", va="center", fontsize=5.0, color=COL["Bad"], fontweight="bold")
            else:
                ax_j.text(jj, ii, f"{v:+.2f}", ha="center", va="center", fontsize=5.0,
                          fontweight="bold" if (ii == 0 or v >= 0.35) else "normal",
                          color="white" if v > 0.38 else COL["Real"])
    ax_j.set_yticks(range(len(M_ORD)))
    ax_j.set_yticklabels([m.replace("CIBERSORTx-", "CSx-").replace("BayesPrism-", "BP-") for m in M_ORD], fontsize=5.0)
    ax_j.get_yticklabels()[0].set_fontweight("bold"); ax_j.get_yticklabels()[0].set_color(COL["Ext"])
    ax_j.set_xticks(range(len(C_ORD))); ax_j.set_xticklabels([c[1] for c in C_ORD], fontsize=5.0, fontweight="bold")
    ax_j.xaxis.tick_top(); ax_j.tick_params(length=0)
    for sp_ in ax_j.spines.values(): sp_.set_visible(False)
    panel_title(ax_j, "Official R 4.3.2 deconv. (Shared-27 ρ)", size=5.9, pad=12)

    ax_k = fig.add_axes([0.482, 0.042, 0.118, 0.142]); fig_panel_label(fig, 0.418, 0.208, "k")
    real_h_med = d["VER_JSON"]["figure6"]["cohort_counts"]["Real_Real_HSPC_median_pct"]
    h_rows = [("Real", real_h_med), ("PBISC", float(off_r.loc["PBISC", "HSPC_median_pct"])),
              ("BP-init", float(off_r.loc["BayesPrism-init", "HSPC_median_pct"])),
              ("BP-upd", float(off_r.loc["BayesPrism-upd", "HSPC_median_pct"])),
              ("Calib-N", float(off_r.loc["Calib-NNLS", "HSPC_median_pct"])),
              ("DWLS", float(off_r.loc["DWLS", "HSPC_median_pct"])),
              ("Bisque", float(off_r.loc["BisqueRNA", "HSPC_median_pct"])),
              ("MuSiC", float(off_r.loc["MuSiC", "HSPC_median_pct"])),
              ("CSx/NN", 0.01)]
    yb_k = np.arange(len(h_rows))[::-1]
    ax_k.barh(yb_k, [r[1] for r in h_rows], height=0.66,
              color=[COL["Real"], COL["Ext"], COL["Bad"], COL["Bad"], COL["Bad"], COL["Bad"], PALETTE["neutral"], PALETTE["neutral"], PALETTE["neutral"]],
              edgecolor="k", lw=0.45)
    for yi, (mn, hv) in zip(yb_k, h_rows):
        ax_k.text(hv + 0.14, yi, f"{hv:.2f}%" if mn != "CSx/NN" else "≤0.1%", va="center", fontsize=5.0,
                  fontweight="bold" if mn in ("Real", "PBISC") or hv > 1.5 else "normal",
                  color=COL["Bad"] if hv > 1.5 else COL["Real"])
    ax_k.set_yticks(yb_k); ax_k.set_yticklabels([r[0] for r in h_rows], fontsize=5.0)
    ax_k.tick_params(axis="y", length=0); ax_k.spines["left"].set_visible(False)
    ax_k.set_xlim(0, 7.8); ax_k.set_xticks([0, 3, 6]); ax_k.set_xlabel("Med HSPC %", fontsize=5.2)
    panel_title(ax_k, "HSPC sink (37–43×)", size=5.8, pad=9)

    ax_l = fig.add_axes([0.688, 0.056, 0.295, 0.128]); fig_panel_label(fig, 0.615, 0.208, "l")
    T5 = d["T5_HSPC"]; t5_g = T5[T5["family"] == "gene"]
    GENES_20 = ["GATA1", "KLF1", "TFRC", "BLVRB", "APOC1", "GATA2", "ITGA2B", "MKI67", "TOP2A", "STMN1", "TYMS",
                "CD34", "AVP", "HLF", "HOPX", "CRHBP", "MPO", "ELANE", "CSF3R", "CEBPD"]
    cmap_g, nrm_g = mpl.colormaps["RdBu_r"], mpl.colors.Normalize(-0.55, 0.55)
    G20_R = [
        ("Real", "sd_gee_wk_beta", "sd_gee_wk_p", "sd_gee_wk_p_qfam", "SD · Real (74)", 0.0),
        ("PBISC new", "sd_gee_wk_beta", "sd_gee_wk_p", "sd_gee_wk_p_qfam", "SD · PBISC (75)", -1.0),
        ("Combined", "sd_gee_wk_beta", "sd_gee_wk_p", "sd_gee_wk_p_qfam", "SD · Comb (149)", -2.0),
        ("Real", "fd_AUC", "fd_p", "fd_p_qfam", "1st · Real (5d)", -3.3),
        ("PBISC new", "fd_AUC", "fd_p", "fd_p_qfam", "1st · PBISC (6d)", -4.3),
        ("Combined", "fd_AUC", "fd_p", "fd_p_qfam", "1st · Comb (11d)", -5.3),
    ]
    ax_l.axhline(-2.65, color=PALETTE["neutral"], ls="--", lw=0.6)
    for ck, vc, pc, qc, rlb, yv in G20_R:
        sub = t5_g[t5_g["set"] == ck].set_index("feature")
        for gi, gn in enumerate(GENES_20):
            rg = sub.loc[gn]
            eff = float(rg[vc]) if vc == "sd_gee_wk_beta" else 2.0 * (float(rg[vc]) - 0.5)
            pv, qv = float(rg[pc]), float(rg[qc])
            sz = 8 + 25 * min(-np.log10(max(pv, 1e-5)), 4.0) / 4.0
            ax_l.scatter(gi, yv, s=sz, marker="^" if eff >= 0 else "v", c=[cmap_g(nrm_g(eff))],
                         edgecolors="k" if qv < 0.05 else ("#555555" if pv < 0.05 else "none"),
                         linewidths=0.78 if qv < 0.05 else 0.32, zorder=3)
    for dx in [6.5, 10.5, 15.5]:
        ax_l.axvline(dx, color=PALETTE["neutral"], ls="--", lw=0.5)
    ax_l.set_yticks([r[5] for r in G20_R]); ax_l.set_yticklabels([r[4] for r in G20_R], fontsize=5.0)
    ax_l.set_xticks(range(len(GENES_20))); ax_l.set_xticklabels(GENES_20, rotation=50, ha="right", fontsize=5.0, fontstyle="italic", fontweight="bold")
    ax_l.tick_params(length=0)
    for sp_ in ax_l.spines.values(): sp_.set_visible(False)
    ax_l.set_xlim(-0.6, len(GENES_20) - 0.4); ax_l.set_ylim(-5.9, 0.85)
    panel_title(ax_l, "3,124 HSPCs: GATA2↓→GATA1/KLF1/TYMS↑", size=5.8, pad=8, color=COL["HSPC"])

    return finalize_figure(fig, FIG_OUT / "Figure6_Rare_Population_Sepsis_COVID_Reproduction")


# =============================================================================
# MAIN FIGURE 7 (a - k): CLINICAL DISCOVERY FOCUSED ON USER'S 3 CORE PILLARS
#   1) Patient & Dataset Metadata Matching (107 donors / 362 1:1 matched samples)
#   2) Population Number (Count & Clinical CBC) vs Percentage (% PBMC & % Parent)
#      Severity (217 draws x Same-Day WHO GEE +wk) & First-Draw Mortality (66 pts)
#   3) Deep T-Cell ISG Scoring Analysis (13 compartments, BayesPrism Z, Time-bias,
#      Sepsis contrast, CTL-ISG anti-correlation, Time x Fatal collapse, 50x5 CV)
# =============================================================================
def build_figure7(d):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.42), facecolor="white")
    T2 = d["T2_ALL"]
    t2_r = T2[T2["set"] == "Real"].set_index("feature")
    t2_c = d["T2B"].set_index("feature")
    t3_f = d["T3_FUN"].set_index("feature")
    A2 = d["A2_COV"].set_index("feature")

    # ---- Row 1 (y=0.835..0.960): PILLAR 1 — a (1:1 Patient/Dataset Metadata Matching) & b (Event Gains + FDR Rescue) ----
    ax_a = fig.add_axes([0.025, 0.840, 0.485, 0.125]); ax_a.set_xlim(0, 100); ax_a.set_ylim(0, 100); ax_a.axis("off")
    arr_kw = dict(arrowstyle="-|>", lw=0.85, color=PALETTE["gray_2"], mutation_scale=6.5)
    rbox(ax_a, 0.5, 6, 25.0, 86, COL["CardGray"], COL["Real"], lw=0.85)
    ax_a.text(13.0, 80, "1:1 Metadata Audit\n107 ESKD Donors", ha="center", va="center", fontsize=5.5, fontweight="bold", color=COL["Real"])
    ax_a.text(13.0, 50, "362 matched samples:\n172 Real bulk+sc (0d)\n164 bulk→PBISC\n26 scRNA-only", ha="center", va="center", fontsize=5.0, color=PALETTE["gray_2"], linespacing=1.18)
    rbox(ax_a, 1.8, 10, 22.4, 21, COL["CardRed"], COL["CLL"], lw=0.65, r=1.1)
    ax_a.text(13.0, 20.5, "−2 CLL · 163/164 sex OK\n0 shared mismatches", ha="center", va="center", fontsize=5.0, fontweight="bold", color=COL["CLL"])
    ax_a.annotate("", xy=(28.5, 49), xytext=(26.0, 49), arrowprops=arr_kw)

    rbox(ax_a, 29.0, 4, 42.5, 90, "white", PALETTE["gray_1"], lw=0.7, ls="--")
    ax_a.text(50.2, 87, "66 Primary COVID+ (217 dr · 11d)", ha="center", va="center", fontsize=5.2, fontweight="bold", color=COL["Ext"])
    rbox(ax_a, 30.2, 27, 19.2, 51, "#EAEAEA", COL["Real"], lw=0.8)
    ax_a.text(39.8, 68, "Real scRNA", ha="center", va="center", fontsize=5.2, fontweight="bold")
    ax_a.text(39.8, 44, "37 pts (10+Sh27)\n138 dr · 5 deaths\nn_cells 226–9,254", ha="center", va="center", fontsize=5.0, color=PALETTE["gray_2"], linespacing=1.15)
    rbox(ax_a, 50.8, 27, 19.5, 51, COL["CardBlue"], COL["New"], lw=0.85)
    ax_a.text(60.5, 68, "PBISC new", ha="center", va="center", fontsize=5.2, fontweight="bold", color=COL["Ext"])
    ax_a.text(60.5, 44, "+29 Wave-1 pts\n79 dr · 6 deaths\n4,000 cells fixed", ha="center", va="center", fontsize=5.0, color=COL["Real"], linespacing=1.15)
    ax_a.text(50.2, 14.5, "41.9% same-day ≠ Peak WHO\nFatal: 4.0 vs 1.3 draws/pt", ha="center", va="center", fontsize=5.0, fontweight="bold", color=COL["HSPC"], linespacing=1.08)
    ax_a.annotate("", xy=(74.3, 49), xytext=(72.0, 49), arrowprops=arr_kw)

    rbox(ax_a, 74.8, 6, 24.5, 86, COL["CardTeal"], COL["Shared"], lw=0.9)
    ax_a.text(87.0, 77, "Primary Pairing", ha="center", va="center", fontsize=5.4, fontweight="bold", color="#165057")
    ax_a.text(87.0, 40, "1) Severity:\n217 dr × SD-WHO\n(GEE +wk, Mundlak)\n2) Death: 66 1st-dr\n3) Count vs % vs CBC\n4) T-cell ISG 13 sub", ha="center", va="center", fontsize=5.0, color=COL["Real"], linespacing=1.15)
    panel_title(ax_a, "1:1 Patient & dataset metadata matching + longitudinal pairing", size=6.0, pad=2)
    fig_panel_label(fig, 0.008, 0.990, "a")

    # b1: Sample & death event gains + b2: 43-subset FDR rescue
    ax_b1 = fig.add_axes([0.602, 0.848, 0.145, 0.112]); fig_panel_label(fig, 0.522, 0.990, "b")
    cats_b = [("WHO 0", 13, 10), ("WHO 1", 9, 5), ("WHO 2", 10, 7), ("WHO 3", 5, 7), ("PUB 11", 5, 6), ("Dir. 9", 3, 6)]
    yb = np.arange(len(cats_b))[::-1]; yb_adj = np.where(yb < 2, yb - 0.30, yb)
    ax_b1.axhline(1.85, color=PALETTE["neutral"], ls="--", lw=0.6)
    for yi, (lab, nr, nn) in zip(yb_adj, cats_b):
        ax_b1.barh(yi, nr, height=0.60, color=COL["Real"], edgecolor="k", lw=0.5)
        ax_b1.barh(yi, nn, left=nr, height=0.60, color=COL["New"], edgecolor="k", lw=0.5)
        ax_b1.text(nr + nn + 0.5, yi, f"{nr+nn}(+{int(round(100*nn/nr))}%)", va="center", fontsize=5.0, fontweight="bold",
                   color=COL["HSPC"] if "PUB" in lab else COL["Ext"])
    ax_b1.set_yticks(yb_adj); ax_b1.set_yticklabels([c[0] for c in cats_b], fontsize=5.0, fontweight="bold")
    ax_b1.tick_params(axis="y", length=0); ax_b1.spines["left"].set_visible(False)
    ax_b1.set_xlim(0, 39); ax_b1.set_xticks([0, 15, 30]); ax_b1.set_xlabel("COVID+ pts", fontsize=5.2)
    panel_title(ax_b1, "Cohort gains", size=5.9, pad=3)

    ax_b2 = fig.add_axes([0.818, 0.848, 0.158, 0.112])
    cap = 5.5; nl = lambda q: np.minimum(-np.log10(np.maximum(q, 1e-8)), cap); q05 = -np.log10(0.05)
    ax_b2.axhspan(q05, cap, xmin=0, xmax=q05 / cap, color=PALETTE["green_1"], alpha=0.85, zorder=0)
    ax_b2.plot([0, cap], [0, cap], color=PALETTE["neutral"], lw=0.65)
    ax_b2.axhline(q05, color=PALETTE["gray_1"], ls=":", lw=0.65); ax_b2.axvline(q05, color=PALETTE["gray_1"], ls=":", lw=0.65)
    ax_b2.axhline(1.0, color=PALETTE["gray_1"], ls="--", lw=0.5, alpha=0.65)
    for f in t3_f.index:
        x_s, y_s = float(nl(t2_r.loc[f, "sd_gee_wk_p_q43"])), float(nl(t2_c.loc[f, "sd_gee_wk_p_q43"]))
        is_rs = bool(t3_f.loc[f, "new_sdWHO_PBISC_rescued"])
        is_as = (t2_r.loc[f, "sd_gee_wk_p_q43"] < 0.05) and (t2_c.loc[f, "sd_gee_wk_p_q43"] < 0.05) and bool(t3_f.loc[f, "new_sdWHO_same_sign"]) and f != "CD14mono_dex"
        fc_s = COL["HSPC"] if f == "HSPC" else (COL["Ext"] if is_rs else (COL["Real"] if is_as else PALETTE["neutral"]))
        ax_b2.scatter(x_s, y_s, s=17 if is_rs else (12 if is_as else 6), marker="o", c=fc_s,
                      edgecolors="k" if is_rs else ("white" if is_as else "none"), lw=0.45 if is_rs else 0.25, zorder=5 if is_rs else 3)
        x_d, y_d = float(nl(t2_r.loc[f, "fd_q43"])), float(nl(t2_c.loc[f, "fd_q43"]))
        is_rd = bool(t3_f.loc[f, "new_fd_PBISC_rescued"])
        ax_b2.scatter(x_d, y_d, s=18 if is_rd else 5, marker="^", c=COL["Ext"] if is_rd else PALETTE["neutral"],
                      edgecolors="k" if is_rd else "none", lw=0.45 if is_rd else 0, zorder=5 if is_rd else 2)
    for f, mode, txt, xyt, tc in [
        ("HSPC", "sev", "HSPC", (0.12, 4.95), COL["HSPC"]),
        ("NK CD56 bright", "sev", "NK56br", (1.65, 3.05), COL["Ext"]),
        ("T CD8 CTL", "dth", "CTL/NK", (1.95, 0.80), COL["Ext"]),
        ("B_ASC_IgG", "sev", "IgG", (1.65, 5.15), COL["Real"]),
    ]:
        xv = float(nl(t2_r.loc[f, "sd_gee_wk_p_q43" if mode == "sev" else "fd_q43"]))
        yv = float(nl(t2_c.loc[f, "sd_gee_wk_p_q43" if mode == "sev" else "fd_q43"]))
        ax_b2.annotate(txt, xy=(xv, yv), xytext=xyt, fontsize=5.0, fontweight="bold", color=tc, va="center",
                       arrowprops=dict(arrowstyle="-", color=PALETTE["gray_1"], lw=0.4, shrinkA=0, shrinkB=2))
    ax_b2.set_xlim(0, cap); ax_b2.set_ylim(0, 5.95); ax_b2.set_xticks([0, 2, 4]); ax_b2.set_yticks([0, 2, 4])
    ax_b2.set_xlabel("−log10 q43 Real", fontsize=5.1); ax_b2.set_ylabel("Comb", fontsize=5.1, labelpad=1)
    panel_title(ax_b2, "FDR rescue (43)", size=5.9, pad=3)

    # ---- Row 2 (y=0.618..0.768): PILLAR 2 — c (Population Number vs % vs CBC check) + d (Stratified, Mundlak & Power) ----
    # c1: 43-subset % PBMC INT vs Raw Count INT scatter (Same-day WHO r=0.987 & First-draw Death r=0.985)
    ax_c1 = fig.add_axes([0.056, 0.625, 0.132, 0.132]); fig_panel_label(fig, 0.008, 0.790, "c")
    ax_c1.plot([-0.4, 0.5], [-0.4, 0.5], color=PALETTE["neutral"], ls="--", lw=0.65)
    ax_c1.axhline(0, color=PALETTE["neutral"], lw=0.45); ax_c1.axvline(0, color=PALETTE["neutral"], lw=0.45)
    sig_both = (A2["sd_q43_pct_INT"] < 0.05) & (A2["sd_q43_cnt_INT"] < 0.05)
    ax_c1.scatter(A2.loc[~sig_both, "sd_beta_pct_INT"], A2.loc[~sig_both, "sd_beta_cnt_INT"], s=7.5, c=PALETTE["neutral"], lw=0, zorder=2)
    ax_c1.scatter(A2.loc[sig_both, "sd_beta_pct_INT"], A2.loc[sig_both, "sd_beta_cnt_INT"], s=13, c=COL["Ext"], edgecolors="k", lw=0.35, zorder=4)
    ax_c1.scatter([A2.loc["HSPC", "sd_beta_pct_INT"]], [A2.loc["HSPC", "sd_beta_cnt_INT"]], s=18, c=COL["HSPC"], edgecolors="k", lw=0.45, zorder=5)
    r_cp = float(stats.pearsonr(A2["sd_beta_pct_INT"], A2["sd_beta_cnt_INT"])[0])
    r_dp = float(stats.pearsonr(A2["fd_AUC_pct_INT"], A2["fd_AUC_cnt_INT"])[0])
    ax_c1.text(0.03, 0.97, f"SD-WHO r={r_cp:.3f}\n18 vs 18 FDR\n1st-dth r={r_dp:.3f}", transform=ax_c1.transAxes,
               fontsize=5.0, fontweight="bold", color=COL["Ext"], va="top", bbox=dict(fc="white", ec="none", alpha=0.88, pad=0.2))
    ax_c1.set_xlim(-0.38, 0.48); ax_c1.set_ylim(-0.42, 0.48); ax_c1.set_xticks([-0.3, 0.0, 0.3]); ax_c1.set_yticks([-0.3, 0.0, 0.3])
    ax_c1.set_xlabel("% PBMC INT β", fontsize=5.1); ax_c1.set_ylabel("Raw Count INT β", fontsize=5.1, labelpad=1)
    panel_title(ax_c1, "Count vs % (43)", size=5.9, pad=3)

    # c2: Multi-metric check across 6 flagship markers: % PBMC INT, Raw Count INT, % Parent INT, Real CBC abs (cells/uL)
    ax_c2 = fig.add_axes([0.278, 0.625, 0.140, 0.132])
    M6 = [("B_ASC_dividing", "B ASC div", "sev"), ("B_ASC_IgG", "B ASC IgG", "sev"), ("HSPC", "HSPC", "sev"),
          ("NK CD56 bright", "NK56br", "sev"), ("T CD8 CM", "T CD8 CM", "sev"), ("T CD8 CTL", "CTL (1st)*", "dth"), ("NK", "NK (1st)*", "dth")]
    yc2 = np.arange(len(M6))[::-1]
    ax_c2.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=0.65); ax_c2.axhline(1.5, color=PALETTE["neutral"], ls="--", lw=0.55)
    for yi, (fk, lb, kind) in zip(yc2, M6):
        rw = A2.loc[fk]
        if kind == "sev":
            vals = [(rw["sd_beta_pct_INT"], rw["sd_p_pct_INT"], COL["Ext"], "s", +0.24),
                    (rw["sd_beta_cnt_INT"], rw["sd_p_cnt_INT"], COL["New"], "o", +0.08),
                    (rw["sd_beta_par_INT"], rw["sd_p_par_INT"], COL["Teal"], "^", -0.08),
                    (rw["sd_beta_absCBC_Real"], rw["sd_p_absCBC_Real"], COL["HSPC"], "D", -0.24)]
        else:
            # Map first-draw AUC to signed effect 2*(AUC-0.5) so it shares the clean [-0.65, +0.55] axis!
            vals = [(2*(rw["fd_AUC_pct_INT"]-0.5), rw["fd_p_pct_INT"], COL["Ext"], "s", +0.24),
                    (2*(rw["fd_AUC_cnt_INT"]-0.5), rw["fd_p_cnt_INT"], COL["New"], "o", +0.08),
                    (2*(rw["fd_AUC_par_INT"]-0.5), rw["fd_p_par_INT"], COL["Teal"], "^", -0.08),
                    (2*(rw["fd_AUC_absCBC_Real"]-0.5), rw["fd_p_absCBC_Real"], COL["HSPC"], "D", -0.24)]
        for v, pv, cc, mk, dy in vals:
            ax_c2.scatter(v, yi + dy, s=11, marker=mk, c=cc if pv < 0.05 else "white", edgecolors=cc, lw=0.65, zorder=4)
    ax_c2.set_yticks(yc2); ax_c2.set_yticklabels([x[1] for x in M6], fontsize=5.0, fontweight="bold")
    ax_c2.tick_params(axis="y", length=0); ax_c2.spines["left"].set_visible(False)
    ax_c2.set_xlim(-0.68, 0.52); ax_c2.set_xticks([-0.5, 0.0, 0.4]); ax_c2.set_xlabel("β (SD) / 2(AUC−.5)", fontsize=5.1)
    ax_c2.legend(handles=[
        Line2D([], [], marker="s", ls="", ms=3.2, color=COL["Ext"], label="%PBMC"),
        Line2D([], [], marker="o", ls="", ms=3.2, color=COL["New"], label="Count"),
        Line2D([], [], marker="^", ls="", ms=3.2, color=COL["Teal"], label="%Par"),
        Line2D([], [], marker="D", ls="", ms=3.2, color=COL["HSPC"], label="CBC"),
    ], loc="lower left", bbox_to_anchor=(-0.22, 1.00), ncol=4, fontsize=5.0, handletextpad=0.05, columnspacing=0.22)
    panel_title(ax_c2, "Metric check", size=5.9, pad=12)

    # d1: Primary Stratified (SD beta + 1st-dr AUC), d2: Mundlak Within vs Between, d3: 300-sim Power Rescue
    ax_d1 = fig.add_axes([0.516, 0.625, 0.095, 0.132]); fig_panel_label(fig, 0.432, 0.790, "d")
    ax_d2 = fig.add_axes([0.620, 0.625, 0.092, 0.132])
    ax_d3 = fig.add_axes([0.742, 0.625, 0.092, 0.132])
    ax_d4 = fig.add_axes([0.858, 0.625, 0.082, 0.132])
    FEATS_D = ["B_ASC_IgG", "B_ASC_IgA", "HSPC", "NK CD56 bright", "T CD8 CM", "T CD8 CTL", "NK"]
    yd_cur = np.arange(len(FEATS_D))[::-1]
    for axx in (ax_d1, ax_d2):
        axx.axhline(1.5, color=PALETTE["neutral"], ls="--", lw=0.55)
    for yi, f in zip(yd_cur, FEATS_D):
        for sname, k, dy, sz in [("Real", "Real", +0.23, 8), ("PBISC new", "New", 0.0, 8.5), ("Combined", "Ext", -0.23, 10)]:
            rw = T2[(T2["feature"] == f) & (T2["set"] == sname)].iloc[0]
            plot_ci_clipped(ax_d1, yi + dy, rw["sd_gee_wk_lo"], rw["sd_gee_wk_hi"], COL[k], (-0.75, 0.75), lw=0.68)
            ax_d1.scatter(rw["sd_gee_wk_beta"], yi + dy, s=sz, marker=MK[k], c="white" if rw["sd_gee_wk_p"] >= 0.05 else COL[k], edgecolors=COL[k], lw=0.55, zorder=3)
            plot_ci_clipped(ax_d2, yi + dy, rw["fd_AUC_lo"], rw["fd_AUC_hi"], COL[k], (0.04, 0.96), lw=0.68)
            ax_d2.scatter(rw["fd_AUC"], yi + dy, s=sz, marker=MK[k], c="white" if rw["fd_p"] >= 0.05 else COL[k], edgecolors=COL[k], lw=0.55, zorder=3)
        # Mundlak within (blue) vs between (grey) on same rows
        wb, bb, wp = t2_c.loc[f, "mundlak_within_beta"], t2_c.loc[f, "mundlak_between_beta"], t2_c.loc[f, "mundlak_within_p"]
        ax_d3.barh(yi + 0.18, wb, height=0.32, color=COL["Ext"], edgecolor="k", lw=0.4)
        ax_d3.barh(yi - 0.18, bb, height=0.32, color=PALETTE["neutral"], edgecolor="k", lw=0.4)
        if wp < 0.05:
            ax_d3.text(wb + (0.03 if wb >= 0 else -0.03), yi + 0.18, "*" if wp >= 0.01 else "**",
                       ha="left" if wb >= 0 else "right", va="center", fontsize=5.0, fontweight="bold", color=COL["Ext"])
    ax_d1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=0.65); ax_d2.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=0.65)
    ax_d3.axvline(0.0, color="#222222", lw=0.6); ax_d3.axhline(1.5, color=PALETTE["neutral"], ls="--", lw=0.55)
    ax_d1.set_yticks(yd_cur); ax_d1.set_yticklabels([LBL[f].replace("NK CD56bright", "NK56br").replace("NK (CD56dim)", "NK dim") for f in FEATS_D], fontsize=5.0, fontweight="bold")
    ax_d1.get_yticklabels()[2].set_color(COL["HSPC"])
    for axx in (ax_d1, ax_d2, ax_d3):
        axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(-0.55, len(FEATS_D) - 0.45)
    ax_d2.set_yticks(yd_cur); ax_d2.set_yticklabels([]); ax_d3.set_yticks(yd_cur); ax_d3.set_yticklabels([])
    ax_d1.set_xlim(-0.75, 0.75); ax_d1.set_xticks([-0.5, 0, 0.5]); ax_d1.set_xlabel("SD β", fontsize=5.1)
    ax_d2.set_xlim(0.04, 0.96); ax_d2.set_xticks([0.2, 0.5, 0.8]); ax_d2.set_xlabel("1st AUC", fontsize=5.1)
    ax_d3.set_xlim(-0.68, 0.55); ax_d3.set_xticks([-0.5, 0, 0.5]); ax_d3.set_xlabel("Mundlak", fontsize=5.1)
    panel_title(ax_d1, "Stratified", size=5.8, pad=3); panel_title(ax_d2, "1st-dr", size=5.8, pad=3); panel_title(ax_d3, "Within/Betw", size=5.8, pad=3)

    # d4: 300-sim Power Rescue aligned on same 7 rows!
    T4 = d["T4_PW"]
    ep_map = {"B_ASC_IgG": "sameday_sev", "B_ASC_IgA": "sameday_sev", "HSPC": "sameday_sev",
              "NK CD56 bright": "sameday_sev", "T CD8 CM": "sameday_sev", "T CD8 CTL": "firstdraw_dth", "NK": "firstdraw_dth"}
    hb = 0.24
    for yi, f in zip(yd_cur, FEATS_D):
        ep = ep_map[f]
        rr = T4[(T4["feature"] == f) & (T4["endpoint"] == ep)].iloc[0]
        for nom_c, bh_c, fc_b, fc_n, dy, ht in [
            ("real_nom", "real_fdr", PALETTE["gray_2"], PALETTE["neutral"], +hb, ""),
            ("plus_nom", "plus_fdr", COL["Ext"], "#B3CBE6", 0.0, ""),
            ("shuf_nom", "shuf_fdr", PALETTE["gray_1"], "white", -hb, "///"),
        ]:
            ax_d4.barh(yi + dy, rr[nom_c], height=hb, color=fc_n, edgecolor="k", lw=0.4, hatch=ht, zorder=2)
            ax_d4.barh(yi + dy, rr[bh_c], height=hb, color=fc_b, edgecolor="k", lw=0.45, zorder=3)
        ax_d4.text(1.04, yi, f"{int(round(100*rr['plus_fdr']))}/{int(round(100*rr['plus_nom']))}%",
                   va="center", fontsize=5.0, fontweight="bold", color=COL["HSPC"] if f == "HSPC" else COL["Ext"])
    ax_d4.axhline(1.5, color=PALETTE["neutral"], ls="--", lw=0.55); ax_d4.axvline(0.8, color=PALETTE["gray_1"], ls=":", lw=0.6)
    ax_d4.set_yticks(yd_cur); ax_d4.set_yticklabels([]); ax_d4.tick_params(axis="y", length=0); ax_d4.spines["left"].set_visible(False)
    ax_d4.set_xlim(0, 1.02); ax_d4.set_xticks([0, 0.5, 1.0]); ax_d4.set_ylim(-0.55, len(FEATS_D) - 0.45); ax_d4.set_xlabel("Power", fontsize=5.1)
    panel_title(ax_d4, "+New29", size=5.8, pad=3)

    # ---- Row 3 (y=0.398..0.535): e (Same-day WHO & 1st-draw Death Strips) & f (External 4-Cohort Meta + Same-Day Labs) ----
    rng_f7 = np.random.default_rng(7)
    ax_e1 = fig.add_axes([0.052, 0.405, 0.112, 0.128]); fig_panel_label(fig, 0.008, 0.562, "e")
    draw_sameday_strip(ax_e1, d, "B_ASC_IgG", rng_f7, show_yl=True, title_txt="B ASC IgG (SD)", title_col=COL["Ext"], compact=True)
    ax_e2 = fig.add_axes([0.172, 0.405, 0.112, 0.128])
    draw_sameday_strip(ax_e2, d, "HSPC", rng_f7, show_yl=False, title_txt="HSPC (SD)", title_col=COL["HSPC"], compact=True)
    ax_e3 = fig.add_axes([0.292, 0.405, 0.104, 0.128])
    draw_firstdraw_death_strip(ax_e3, d, "T CD8 CTL", rng_f7, show_yl=False, title_txt="T CD8 CTL (1st)", title_col=COL["Ext"], compact=True)
    ax_e4 = fig.add_axes([0.404, 0.405, 0.104, 0.128])
    draw_firstdraw_death_strip(ax_e4, d, "NK", rng_f7, show_yl=False, title_txt="NK dim (1st)", title_col=COL["Ext"], compact=True)

    # f1: External 4-cohort public COVID-19 meta-analysis (Severity N=220 + HSPC >7d Death)
    PUB_META, PUB_HW, ESKD_MATCH = d["PUB_META"], d["PUB_HW"], d["ESKD_MATCH"]
    ax_f1 = fig.add_axes([0.592, 0.405, 0.132, 0.128]); fig_panel_label(fig, 0.520, 0.562, "f")
    EXT_S = [("NK_CD56bright", "NK56br", "NK CD56 bright", "anc"), ("B_ASC_dividing", "B ASC div", "B_ASC_dividing", "anc"),
             ("B_ASC_total", "B ASC tot", "B_ASC_total", "foc"), ("T_CD8_CM", "T CD8 CM", "T CD8 CM", "anc"),
             ("T_CD8_CTL", "T CD8 CTL", "T CD8 CTL", "foc"), ("HSPC", "HSPC", "HSPC", "hspc")]
    ye_i = np.arange(len(EXT_S))[::-1]; xl_s = (-0.75, 0.78); ax_f1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=0.65)
    for yi, (pf, disp, ef, grp) in zip(ye_i, EXT_S):
        em = ESKD_MATCH.loc[ef]
        mr = PUB_META[(PUB_META["feature"] == pf) & (PUB_META["endpoint"] == "severity")].iloc[0]
        pcol = COL["HSPC"] if grp == "hspc" else COL["Teal"]
        ax_f1.plot([em["lo"], em["hi"]], [yi + 0.16, yi + 0.16], color=COL["Ext"], lw=0.7)
        ax_f1.scatter(em["rho"], yi + 0.16, s=9.5, marker="s", c=COL["Ext"], zorder=3, label="ESKD" if yi == ye_i[0] else "")
        plot_ci_clipped(ax_f1, yi - 0.16, float(mr["DL_lo"]), float(mr["DL_hi"]), pcol, xl_s, lw=0.85)
        ax_f1.scatter(mr["estimate"], yi - 0.16, s=11, marker="D", c="white" if mr["DL_p"] >= 0.05 else pcol, edgecolors=pcol, lw=0.65, zorder=4, label="4-coh" if yi == ye_i[0] else "")
    ax_f1.set_yticks(ye_i); ax_f1.set_yticklabels([x[1] for x in EXT_S], fontsize=5.0, fontweight="bold")
    ax_f1.tick_params(axis="y", length=0); ax_f1.spines["left"].set_visible(False)
    ax_f1.set_xlim(*xl_s); ax_f1.set_xticks([-0.5, 0, 0.5]); ax_f1.set_xlabel("Ext. Sev ρ (N=220)", fontsize=5.1)
    ax_f1.legend(loc="lower left", bbox_to_anchor=(-0.08, 1.00), ncol=2, fontsize=5.0, handletextpad=0.08, columnspacing=0.3)
    panel_title(ax_f1, "4-cohort meta", size=5.8, pad=12)

    # f2: Same-day hospital labs GEE (+wk, BH36)
    ax_f2 = fig.add_axes([0.790, 0.405, 0.190, 0.128])
    T7 = d["T7_LABS"]
    IMM6 = ["NK CD56 bright", "T CD8 CM", "T CD8 CTL", "B_ASC_IgG", "B_ASC_IgA", "HSPC"]
    LABS6 = [("temp_CRP", "CRP"), ("temp_lymph", "ALC"), ("NLR", "NLR"), ("temp_ferritin", "Fer"), ("temp_ddimer", "D-d"), ("temp_trop", "Trp")]
    mat_l, p_l, q_l = np.zeros((6, 6)), np.zeros((6, 6)), np.zeros((6, 6))
    for i_f, fk in enumerate(IMM6):
        for j_l, (lk, _) in enumerate(LABS6):
            rw = T7[(T7["feature"] == fk) & (T7["lab"] == lk)].iloc[0]
            mat_l[i_f, j_l], p_l[i_f, j_l], q_l[i_f, j_l] = rw["sameday_gee_wk_beta"], rw["sameday_gee_wk_p"], rw["sameday_gee_wk_q36"]
    ax_f2.imshow(mat_l, cmap="RdBu_r", norm=mpl.colors.Normalize(-0.68, 0.68), aspect="auto")
    for i_f in range(6):
        for j_l in range(6):
            vv, sg, qg = mat_l[i_f, j_l], p_l[i_f, j_l] < 0.05, q_l[i_f, j_l] < 0.05
            star = "**" if qg else ("*" if sg else "")
            ax_f2.text(j_l, i_f, f"{vv:+.2f}{star}", ha="center", va="center", fontsize=5.0,
                       fontweight="bold" if sg else "normal", color="white" if abs(vv) > 0.40 else COL["Real"])
    ax_f2.axhline(4.5, color="k", lw=0.7)
    ax_f2.set_yticks(range(6)); ax_f2.set_yticklabels([LBL[f].replace("NK CD56bright", "NK56br").replace("T CD8 ", "CD8 ") for f in IMM6], fontsize=5.0, fontweight="bold")
    ax_f2.set_xticks(range(6)); ax_f2.set_xticklabels([x[1] for x in LABS6], fontsize=5.0, fontweight="bold")
    ax_f2.xaxis.tick_top(); ax_f2.tick_params(length=0)
    for sp_ in ax_f2.spines.values(): sp_.set_visible(False)
    panel_title(ax_f2, "Same-day labs GEE (**q36)", size=5.8, pad=12)

    # ---- Row 4 (y=0.048..0.312): PILLAR 3 — DEDICATED DEEP T-CELL ISG SCORING ANALYSIS (g, h, i, j) ----
    # Banner divider framing Pillar 3 cleanly
    fig.add_artist(Line2D([0.012, 0.988], [0.352, 0.352], color=COL["Ext"], lw=0.85, ls="-", alpha=0.45))
    fig.text(0.50, 0.356, "T-cell 19-gene ISG scoring: reconstruction, time-bias dissection, CTL-ISG anti-correlation, fatal collapse & prognosis",
             ha="center", va="bottom", fontsize=5.6, fontweight="bold", color=COL["Ext"])

    # g: T-cell ISG multi-compartment audit (left: PBISC vs BayesPrism Z + Time-bias removal; bottom: Sepsis T-ISG suppression contrast)
    ax_g1 = fig.add_axes([0.105, 0.188, 0.140, 0.125]); fig_panel_label(fig, 0.008, 0.342, "g")
    t6_c = d["T6_ISG"][d["T6_ISG"]["set"] == "Combined"].set_index("subset")
    t8b = d["T8B_WIN"]
    TS5 = ["T CD8 CTL", "T CD4 CTL", "All CD8 T", "All CD4 T", "All T cells"]
    yy_5 = np.arange(len(TS5))[::-1]
    ax_g1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=0.65)
    for yi, ts in zip(yy_5, TS5):
        rt = t6_c.loc[ts]
        r_e7 = t8b[(t8b["subset"] == ts) & (t8b["window"] == "≤7 d")].iloc[0]
        ax_g1.plot([rt["sd_gee_wk_beta"], rt["pm_peakWHO_rho"]], [yi, yi], color=PALETTE["neutral"], lw=0.85)
        ax_g1.scatter(rt["pm_peakWHO_rho"], yi, s=11, marker="D", c="white", edgecolors=PALETTE["gray_2"], lw=0.7, zorder=3, label="PM ρ" if yi == yy_5[0] else "")
        ax_g1.scatter(r_e7["gee_sev_beta"], yi, s=11, marker="^", c=COL["Teal"], zorder=4, label="≤7d β*" if yi == yy_5[0] else "")
        ax_g1.scatter(rt["sd_gee_wk_beta"], yi, s=13, marker="o", c=COL["Ext"], zorder=5, label="SD+wk" if yi == yy_5[0] else "")
    ax_g1.set_yticks(yy_5); ax_g1.set_yticklabels([t.replace("All T cells", "All T").replace("All ", "") for t in TS5], fontsize=5.0, fontweight="bold")
    ax_g1.tick_params(axis="y", length=0); ax_g1.spines["left"].set_visible(False)
    ax_g1.set_xlim(-0.12, 0.42); ax_g1.set_xticks([0.0, 0.2, 0.4]); ax_g1.set_xlabel("ISG sev. effect", fontsize=5.1)
    ax_g1.legend(loc="lower left", bbox_to_anchor=(-0.20, 1.00), ncol=3, fontsize=5.0, handletextpad=0.05, columnspacing=0.22)
    panel_title(ax_g1, "Time-bias check", size=5.8, pad=12)

    # g2 (directly below g1): Cross-disease Bacterial Sepsis PBISC T-ISG suppression vs COVID-19 + PBISC vs BayesPrism Z
    ax_g2 = fig.add_axes([0.068, 0.048, 0.178, 0.095])
    sep_isg = d["A3_SEP_ISG"]
    sep_ad = sep_isg[sep_isg["in_primary_scope"].astype(str).isin(["True", "1"])].copy()
    rng_sp = np.random.default_rng(42)
    for xi, (ck, clb, cc) in enumerate([("adult_HC", "HC\n(36)", PALETTE["gray_2"]), ("adult_ICU_ctrl", "ICU\n(17)", COL["Teal"]), ("adult_sepsis", "Sepsis\n(68)", COL["HSPC"])]):
        vv = sep_ad[sep_ad["condition"] == ck]["T_ISG_score"].values
        ax_g2.scatter(xi + rng_sp.uniform(-0.18, 0.18, len(vv)), vv, s=5.5, c=cc, alpha=0.72, lw=0)
        ax_g2.plot([xi - 0.24, xi + 0.24], [np.median(vv)] * 2, color="k", lw=1.0, zorder=4)
    ax_g2.text(0.98, 0.96, "ρ=−0.39, P=1.1e-5\n(Viral-specific ISG)", transform=ax_g2.transAxes, ha="right", va="top",
               fontsize=5.0, fontweight="bold", color=COL["HSPC"], bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15))
    ax_g2.set_xticks([0, 1, 2]); ax_g2.set_xticklabels(["HC(36)", "ICU(17)", "Sep(68)"], fontsize=5.0, fontweight="bold")
    ax_g2.set_ylabel("T-ISG", fontsize=5.1, labelpad=1); ax_g2.set_ylim(0.32, 0.60)
    panel_title(ax_g2, "Sepsis T-ISG suppression", size=5.7, pad=2)

    # h: Baseline Exact First-Draw T CD8 CTL Abundance vs T CD8 CTL ISG Anti-Correlation (n=63 exact draw 1 >=5 CTL, rho=-0.409, P=8.7e-4)
    ax_h = fig.add_axes([0.305, 0.055, 0.175, 0.250]); fig_panel_label(fig, 0.258, 0.342, "h")
    FD63 = d["A3_FD66"]
    for arm_v, sv, mkv in [("Real", "Real", "D"), ("PBISC new", "New", "o")]:
        for fv, cc, szv in [(0, COL[sv], 10), (1, COL["HSPC"], 22)]:
            sb = FD63[(FD63["arm"] == arm_v) & (FD63["fatal_pub11"] == fv)]
            ax_h.scatter(sb["ctl_zfd"], sb["isg_zfd"], s=szv, marker=mkv, c=cc,
                         edgecolors="k" if fv == 1 else "white", lw=0.45, alpha=0.72 if fv == 0 else 0.95, zorder=5 if fv == 1 else 2)
    # Regression line
    xr = np.linspace(-2.2, 2.2, 50)
    sl, ic = np.polyfit(FD63["ctl_zfd"], FD63["isg_zfd"], 1)
    ax_h.plot(xr, sl * xr + ic, color=COL["Ext"], lw=0.95, ls="--", zorder=3)
    ax_h.axhline(0.0, color=PALETTE["neutral"], ls=":", lw=0.6); ax_h.axvline(0.0, color=PALETTE["neutral"], ls=":", lw=0.6)
    ax_h.axvspan(-2.5, 0.0, ymin=0.5, ymax=1.0, color=COL["CardRed"], alpha=0.25, zorder=0)
    tsum = d["A3_JSON"]["tisg_summary"]
    ax_h.text(0.97, 0.97, f"1st-dr ρ={tsum['ctl_vs_isg_63exact_firstdraw_rho']:+.2f}\nP=8.7e-4 (n=63)\nAll 213dr ρ=−0.26",
              transform=ax_h.transAxes, ha="right", va="top", fontsize=5.0, fontweight="bold", color=COL["Ext"],
              bbox=dict(fc="white", ec="none", alpha=0.90, pad=0.2), zorder=7)
    ax_h.text(0.97, 0.03, f"−(CTL+NK)+ISG\nAUC={tsum['firstdraw63_AUC_negCTLNK_posISG']:.3f}\nP=4.3e-4",
              transform=ax_h.transAxes, ha="right", va="bottom", fontsize=5.0, fontweight="bold", color=COL["HSPC"],
              bbox=dict(fc="white", ec="none", alpha=0.90, pad=0.2), zorder=7)
    ax_h.set_xlim(-2.45, 2.45); ax_h.set_ylim(-2.45, 2.55); ax_h.set_xticks([-2, 0, 2]); ax_h.set_yticks([-2, 0, 2])
    ax_h.set_xlabel("1st-draw T CD8 CTL (z)", fontsize=5.2); ax_h.set_ylabel("1st-draw T CD8 CTL ISG (z)", fontsize=5.2, labelpad=1)
    panel_title(ax_h, "CTL loss ↔ ISG coupling (n=63)", size=5.9, pad=3)

    # i: Longitudinal T CD8 CTL ISG fatal collapse (213 draws) + 5-compartment GEE Time x Fatal forest!
    ax_i1 = fig.add_axes([0.535, 0.198, 0.215, 0.118]); fig_panel_label(fig, 0.492, 0.342, "i")
    x_grid = np.linspace(0, 28, 100)
    dm = d["CD8_ISG_217"][["patient_id", "source", "days", "isg_z", "fatal_pub11"]].dropna().rename(columns={"isg_z": "y"})
    ax_i1.axvline(7.0, color=PALETTE["gray_1"], ls="--", lw=0.6)
    for fv, cln, lb_g, al_p in [(0, COL["New"], "Surv (55)", 0.32), (1, COL["HSPC"], "Died (11)", 0.85)]:
        sub = dm[dm["fatal_pub11"] == fv]
        for sk, mkv in [("Real", "D"), ("PBISC", "o")]:
            sb = sub[sub["source"] == sk]
            ax_i1.scatter(sb["days"], sb["y"], s=4.8 if fv == 0 else 7.5, marker=mkv, c=cln, alpha=al_p, lw=0, rasterized=True)
        ols_m = sm.OLS(sub["y"], sm.add_constant(sub["days"])).fit(cov_type="cluster", cov_kwds={"groups": sub["patient_id"]})
        pr = ols_m.get_prediction(sm.add_constant(x_grid)).summary_frame()
        ax_i1.fill_between(x_grid, pr["mean_ci_lower"], pr["mean_ci_upper"], color=cln, alpha=0.16)
        ax_i1.plot(x_grid, pr["mean"], color=cln, lw=1.1, label=lb_g)
    ax_i1.text(0.03, 0.96, "Time×Fatal β=−0.55/wk\nP=3.7e-6 (Real −0.53, New −0.51)", transform=ax_i1.transAxes,
               fontsize=5.0, fontweight="bold", color=COL["HSPC"], va="top", bbox=dict(fc="white", ec="none", alpha=0.86, pad=0.2))
    ax_i1.set_xlim(-0.8, 30); ax_i1.set_ylim(-2.6, 3.6); ax_i1.set_xticks([0, 7, 14, 21, 28]); ax_i1.set_yticks([-2, 0, 2])
    ax_i1.set_xlabel("Days from 1st + swab", fontsize=5.1, labelpad=1.5); ax_i1.set_ylabel("CTL ISG (z)", fontsize=5.1, labelpad=1)
    ax_i1.legend(loc="lower left", fontsize=5.0, handlelength=0.8)
    panel_title(ax_i1, "T CD8 CTL ISG collapse (213 dr)", size=5.8, pad=3)

    # i2 (below i1): 5-compartment GEE Time x Fatal slopes (Real / New / Combined)
    ax_i2 = fig.add_axes([0.585, 0.042, 0.165, 0.088])
    TISG_GEE = d["TISG_GEE"]
    ax_i2.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=0.65)
    for yi, ts in zip(yy_5, TS5):
        sub_tg = TISG_GEE[TISG_GEE["subset"] == ts].set_index("cohort")
        for clb, ks, dy in [("Real", "Real", +0.24), ("PBISC new", "New", 0.0), ("Combined", "HSPC", -0.24)]:
            rg = sub_tg.loc[clb]
            ax_i2.plot([rg["fat_x_wk_lo"], rg["fat_x_wk_hi"]], [yi + dy, yi + dy], color=COL[ks], lw=0.72)
            ax_i2.scatter(rg["fat_x_wk_beta"], yi + dy, s=9.5, marker=MK["Ext" if ks == "HSPC" else ks], c=COL[ks], zorder=4)
    ax_i2.set_yticks(yy_5); ax_i2.set_yticklabels([t.replace("All T cells", "All T").replace("All ", "") for t in TS5], fontsize=5.0, fontweight="bold")
    ax_i2.tick_params(axis="y", length=0); ax_i2.spines["left"].set_visible(False)
    ax_i2.set_xlim(-1.05, 0.08); ax_i2.set_xticks([-0.8, -0.4, 0.0]); ax_i2.set_xlabel("Time×Fatal β/wk", fontsize=5.1)
    panel_title(ax_i2, "5-comp. fatal collapse", size=5.7, pad=2)

    # j: Leak-free In-fold rank-INT 50x5 Canonical Exact First-Draw Mortality CV (n=66 / <=7d n=53 / ISG exact draw-1 n=63)
    ax_j = fig.add_axes([0.840, 0.055, 0.098, 0.250]); fig_panel_label(fig, 0.762, 0.342, "j")
    A3CV = d["A3_CV"].set_index("model")
    CV500 = d["CV_500K"]
    cv_rows_j = [
        ("1st dr (63)\n+CTL ISG", A3CV.loc["Exact 1st-dr (63) + T CD8 CTL ISG", "AUC_base"], A3CV.loc["Exact 1st-dr (63) + T CD8 CTL ISG", "AUC_ext"], COL["Teal"], True),
        ("1st dr (66)\n+T CD8 CTL", A3CV.loc["All 1st-dr (66) + T CD8 CTL", "AUC_base"], A3CV.loc["All 1st-dr (66) + T CD8 CTL", "AUC_ext"], COL["Ext"], True),
        ("1st dr (66)\n+CTL + NK", A3CV.loc["All 1st-dr (66) + CTL + NK", "AUC_base"], A3CV.loc["All 1st-dr (66) + CTL + NK", "AUC_ext"], COL["Ext"], True),
        ("≤7d (53)\n+T CD8 CTL", A3CV.loc["Early <=7d (53) + T CD8 CTL", "AUC_base"], A3CV.loc["Early <=7d (53) + T CD8 CTL", "AUC_ext"], COL["Ext"], True),
        ("WHO 2–3\n+HSPC %", float(CV500.iloc[4]["AUC_base"]), float(CV500.iloc[4]["AUC_ext"]), COL["HSPC"], False),
    ]
    yy_j = np.arange(len(cv_rows_j))[::-1]; ax_j.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=0.65)
    for yi, (lb, ab, ae, cc, is_good) in zip(yy_j, cv_rows_j):
        da = ae - ab
        ax_j.plot([ab, ae], [yi, yi], color=cc, lw=1.0, zorder=2)
        ax_j.scatter(ab, yi, s=13, marker="D", c="white", edgecolors=PALETTE["gray_2"], lw=0.75, zorder=3, label="Base" if yi == yy_j[0] else "")
        ax_j.scatter(ae, yi, s=16, marker="s", c=cc if is_good else "white", edgecolors=cc, lw=0.75, zorder=4, label="+Imm" if yi == yy_j[0] else "")
        ax_j.text(0.95, yi, f"{ae:.2f}\n(Δ{da:+.2f})", va="center", fontsize=5.0,
                  fontweight="bold" if is_good else "normal", color=cc if is_good else PALETTE["gray_2"], linespacing=0.95)
    ax_j.set_yticks(yy_j); ax_j.set_yticklabels([x[0] for x in cv_rows_j], fontsize=5.0, fontweight="bold")
    ax_j.tick_params(axis="y", length=0); ax_j.spines["left"].set_visible(False)
    ax_j.set_xlim(0.45, 0.93); ax_j.set_xticks([0.5, 0.7, 0.9]); ax_j.set_ylim(-0.55, len(cv_rows_j) - 0.35)
    ax_j.set_xlabel("50×5 CV AUC", fontsize=5.1)
    ax_j.legend(loc="lower left", bbox_to_anchor=(-0.15, 1.00), ncol=2, fontsize=5.0, handletextpad=0.08, columnspacing=0.25)
    panel_title(ax_j, "1st-dr CV", size=5.9, pad=12)

    return finalize_figure(fig, FIG_OUT / "Figure7_COVID_Clinical_Metadata_CountPct_TISG")


# =============================================================================
# UPGRADED SUPPLEMENTARY FIGURE 6: EXACT PATIENT & DATASET METADATA MATCHING TABLES
#   (Fixes Sepsis Adult Real to audited 93/8/7 + exact age/sex + 0 overlap + sex_flag
#    + COVID 362-sample 1:1 matching & n_cells + individual clinical audit table!)
# =============================================================================
def draw_pub_table(ax, headers, rows, col_widths, col_aligns, header_bg="#0F4D92", stripe_bg="#F3F7FB",
                   highlight_rows=None, bold_rows=None, row_height=1.0, font_pt=5.0):
    highlight_rows = highlight_rows or {}
    bold_rows = bold_rows or set()
    n_rows = len(rows)
    total_w = sum(col_widths)
    x_edges = [0.0]
    for w in col_widths:
        x_edges.append(x_edges[-1] + w)
    ax.set_xlim(0, total_w)
    ax.set_ylim(-(n_rows + 1.25) * row_height, 0.15 * row_height)
    ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")

    ax.add_patch(Rectangle((0, -1.25 * row_height), total_w, 1.25 * row_height, facecolor=header_bg, edgecolor="none", zorder=1))
    for j, htxt in enumerate(headers):
        cx = x_edges[j] + (0.010 * total_w if col_aligns[j] == "left" else 0.5 * col_widths[j])
        ax.text(cx, -0.625 * row_height, htxt, ha=col_aligns[j], va="center",
                fontsize=font_pt + 0.1, fontweight="bold", color="white", zorder=3)

    for i, rvals in enumerate(rows):
        y_top = -(1.25 + i) * row_height
        y_mid = y_top - 0.5 * row_height
        bg = highlight_rows.get(i, stripe_bg if (i % 2 == 1) else "white")
        ax.add_patch(Rectangle((0, y_top - row_height), total_w, row_height, facecolor=bg, edgecolor="none", zorder=1))
        ax.plot([0, total_w], [y_top - row_height, y_top - row_height], color="#CBD5E1", lw=0.4, zorder=2)
        is_b = i in bold_rows
        for j, val in enumerate(rvals):
            cx = x_edges[j] + (0.010 * total_w if col_aligns[j] == "left" else 0.5 * col_widths[j])
            sval = str(val)
            txt_c = COL["HSPC"] if ("Died" in sval or "CLL" in sval or "Mismatch" in sval) else (COL["Good"] if "ALIVE" in sval else "#111111")
            ax.text(cx, y_mid, sval, ha=col_aligns[j], va="center",
                    fontsize=font_pt, fontweight="bold" if (is_b or j == 0 or "Died" in sval) else "normal",
                    color=txt_c, zorder=3)
    ax.plot([0, total_w], [-(1.25 + n_rows) * row_height, -(1.25 + n_rows) * row_height], color="#0F4D92", lw=0.9, zorder=3)


def build_suppfig6_metadata(d):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.38), facecolor="white")
    t1 = d["T1_PAT68"]; t1c = d["T1C_SUM"].set_index("group"); t1d = d["T1D_NEG39"]

    def age_sex(sub):
        a = sub["age"].dropna()
        m = int((sub["sex"] == "M").sum())
        return f"{a.median():.0f} ({a.quantile(0.25):.0f}–{a.quantile(0.75):.0f})", f"{m}/{len(sub)} ({100*m/max(len(sub),1):.0f}%)"

    # Panel a: Master Patient & Dataset Metadata Matching Table (Sepsis + ESKD COVID-19)
    ax_a = fig.add_axes([0.020, 0.505, 0.962, 0.448]); fig_panel_label(fig, 0.008, 0.985, "a")
    headers_a = [
        "Cohort / Dataset Stratum", "Pts", "Draws", "Dr/pt\n(IQR)", "Cells/dr\nmed (IQR)",
        "Age yr\n(IQR)", "Male\nn (%)", "Sex-QC /\nOverlap", "1st-dr\n(≤7d)",
        "Same-day WHO / Group", "Peak WHO", "Deaths\nPUB11(D9)"
    ]
    cw_a = [22.2, 4.0, 4.8, 7.6, 9.8, 8.2, 8.2, 8.2, 7.8, 12.2, 7.8, 7.2]
    ca_a = ["left", "center", "center", "center", "center", "center", "center", "center", "center", "center", "center", "center"]

    r_only = t1[t1["cohort_group"] == "Real-only (10)"]
    sh27 = t1[t1["cohort_group"] == "Shared-27 (Real used)"]
    real_all = t1[t1["primary_source"] == "Real"]
    pb_new = t1[t1["primary_source"] == "PBISC new"]
    comb_all = t1[t1["primary_source"] != "Excluded_CLL"]
    cll_2 = t1[t1["primary_source"] == "Excluded_CLL"]
    ag_neg, sx_neg = age_sex(t1d); ag_cll, sx_cll = age_sex(cll_2)

    def row_cov(label, key_c, sub_df, d9_cnt, cells_s, qc_s):
        rc = t1c.loc[key_c]
        ag, sx = age_sex(sub_df)
        return [
            label, str(int(rc["n_patients"])), str(int(rc["n_draws"])),
            f"{rc['draws_per_pt_median']:.1f} ({rc['draws_per_pt_iqr']})", cells_s,
            ag, sx, qc_s, f"d{rc['first_draw_day_median']:.0f} ({int(rc['first_draw_le7d_pts'])}/{int(rc['n_patients'])})",
            str(rc["sameday_WHO_draws_m_mo_s_c"]), str(rc["peak_WHO_mild_mod_sev_crit"]),
            f"{int(rc['deaths_pub11'])} ({d9_cnt})"
        ]

    rows_a = [
        # EXACT audited Bacterial Sepsis (Ye et al., GSE279452 / GSE279448)
        ["Sepsis · Real scRNA (adult)", "108", "108", "1.0 (1–1)", "5,664 (4.7–6.8k)", "67 (54–74)", "76/108 (70%)", "GEO GSM OK", "1-dr", "Sep 93 / ICU 8 / HC 7", "—", "—"],
        ["Sepsis · Real scRNA (paed.)", "25", "25", "1.0 (1–1)", "5,460 (4.1–6.5k)", "4 (1–9)", "17/25 (68%)", "GEO GSM OK", "1-dr", "Sep 17 / ICU 4 / HC 4", "—", "—"],
        ["Sepsis · PBISC (adult bulk*)", "121", "121", "1.0 (1–1)", "4,000 (fixed)", "65 (54–74)", "78/121 (64%)", "119/122 · 0 ovl", "1-dr", "Sep 68 / ICU 17 / HC 36", "—", "—"],
        # EXACT audited Imperial ESKD COVID-19 (107 donors / 362 matched samples)
        row_cov("ESKD · Real-only (NCL/Cam)", "Real-only (10 pts)", r_only, int(r_only["death_DIRECT9"].sum()), "2,718 (1.2–4.5k)", "0 mismatches"),
        row_cov("ESKD · Shared-27 (Real used)", "Shared-27 (Real used, 27 pts)", sh27, int(sh27["death_DIRECT9"].sum()), "2,078 (1.1–3.8k)", "59 PBISC val."),
        row_cov("ESKD · Real ALL (Primary)", "Real ALL (37 pts)", real_all, int(real_all["death_DIRECT9"].sum()), "2,204 (1.1–4.1k)", "0 mismatches"),
        row_cov("  ├─ Real · Survived", "Real · Survived (32 pts)", real_all[real_all["death_TRUE_PUB11"] == 0], 0, "2,112 (1.0–3.9k)", "118 dr"),
        row_cov("  └─ Real · Died", "Real · Died (5 pts)", real_all[real_all["death_TRUE_PUB11"] == 1], int(real_all["death_DIRECT9"].sum()), "2,761 (1.7–4.8k)", "4.0 dr/pt"),
        row_cov("ESKD · PBISC new (Wave-1)", "PBISC new ALL (29 pts)", pb_new, int(pb_new["death_DIRECT9"].sum()), "4,000 (fixed)", "79/79 sex OK"),
        row_cov("  ├─ PBISC new · Survived", "PBISC new · Survived (23 pts)", pb_new[pb_new["death_TRUE_PUB11"] == 0], 0, "4,000 (fixed)", "71 dr"),
        row_cov("  └─ PBISC new · Died", "PBISC new · Died (6 pts)", pb_new[pb_new["death_TRUE_PUB11"] == 1], int(pb_new["death_DIRECT9"].sum()), "4,000 (fixed)", "1.3 dr/pt"),
        row_cov("ESKD · Combined Primary", "Combined ALL (66 pts)", comb_all, int(comb_all["death_DIRECT9"].sum()), "2.2k / 4.0k", "163/164 sex OK"),
        ["ESKD · Excluded CLL (2 pts)", "2", "11†", "3–6", "Var. / 4,000", ag_cll, sx_cll, "ROR1/FMOD+", "d3,d7 (2/2)", "Excluded (monoclonal)", "0/0/1/1", "0 (0)"],
        ["ESKD · COVID-neg ESKD ctrl", "39", "39", "1.0 (1–1)", "Var. / 4,000", ag_neg, sx_neg, "38R / 19P‡", "Baseline", "COVID-negative ESKD", "—", "0 (0)"],
    ]
    hl_a = {2: "#EFF6FF", 5: "#E8EFF8", 7: "#FCE8E6", 8: "#E8EFF8", 10: "#FCE8E6", 11: "#DCE8F6", 12: "#FEF3C7"}
    draw_pub_table(ax_a, headers_a, rows_a, cw_a, ca_a, highlight_rows=hl_a, bold_rows={0, 2, 5, 8, 11}, font_pt=5.0)
    panel_title(ax_a, "Exact patient & dataset metadata matching across Bacterial Sepsis (255 samples) and ESKD COVID-19 (107 donors / 362 matched samples)", size=6.0, pad=3)

    # Panel b: Individual Clinical, Dataset-Matching & QC Flag Table (11 fatal + 3 resolved + 2 CLL + 2 sex-QC flags)
    ax_b = fig.add_axes([0.020, 0.056, 0.962, 0.408]); fig_panel_label(fig, 0.008, 0.489, "b")
    headers_b = [
        "Patient / ID", "Dataset / Role", "Site/\nWave", "Age/\nSex", "Dr", "1st\nDay",
        "Blood-Draw Days (+swab)", "Same-Day WHO Sequence", "Peak\nWHO", "Steroid", "Forensic / QC Status"
    ]
    cw_b = [7.6, 13.5, 6.8, 6.2, 3.6, 4.4, 16.5, 22.0, 6.2, 7.6, 12.6]
    ca_b = ["left", "left", "center", "center", "center", "center", "left", "left", "center", "center", "center"]

    target_ids = ["C141", "C21", "C33", "C108", "C145", "C32", "C37", "C48", "C66", "C77", "C84",
                  "C101", "C146", "C147", "C116", "C139"]
    t1_ix = t1.set_index("patient_id")
    rows_b = []; hl_b = {}
    for idx_i, pid in enumerate(target_ids):
        r = t1_ix.loc[pid]
        grp_s = str(r["cohort_group"]).replace(" (Real used)", " (Real)").replace(" (10)", "").replace(" (29)", "").replace(" (2)", "")
        ag_s = f"{int(r['age'])}/{r['sex']}"
        if pid in ("C116", "C139"):
            dr_s = "3R+2P" if pid == "C116" else "6R"; fd_s = "d3" if pid == "C116" else "d7"
            days_s = "d3, d8, d12 (Real)" if pid == "C116" else "d7, d10, d12, d14, d21, d24"
            sd_s = "mod, mild, crit" if pid == "C116" else "mild×4, sev×2"; st_s = "Excl."; out_s = "Excluded (CLL)"; hl_b[idx_i] = "#FEF3C7"
        else:
            dr_s = str(int(r["n_primary_draws"])); fd_s = f"d{int(r['first_draw_day'])}"
            days_s = str(r["draw_days_swab"])
            sd_s = str(r["same_day_WHO_each_draw"]).replace("moderate", "mod").replace("severe", "sev").replace("critical", "crit")
            st_s = "Unrec. (W1)" if "UNREC" in str(r["recorded_steroid"]) else str(r["recorded_steroid"])
            if int(r["death_TRUE_PUB11"]) == 1:
                out_s = "Died (Direct)" if int(r["death_DIRECT9"]) == 1 else "Died (Indirect)"; hl_b[idx_i] = "#FCE8E6"
            else:
                out_s = "ALIVE (2m rec.)" if pid in ("C101", "C146") else "ALIVE (Tab. S1)"; hl_b[idx_i] = "#E6F4EA"
        rows_b.append([pid, grp_s, str(r["centre_or_wave"]), ag_s, dr_s, fd_s, days_s, sd_s, str(r["peak_WHO"]), st_s, out_s])
    # Add molecular sex-QC audit rows explicitly!
    rows_b.append(["C137_pos_2", "COVID Shared-27 (P)", "w2", "68/F", "1P", "d2", "d2 (non-primary PBISC)", "ambiguous bulk sex call", "mild", "N", "Sex-QC: Ambiguous"])
    hl_b[len(rows_b) - 1] = "#EFF6FF"
    rows_b.append(["Res-S34", "Sepsis PBISC bulk", "Chongq.", "71/M", "1P", "—", "Cross-sectional", "GEO Male vs bulk female", "Sepsis", "—", "Sex-QC: Mismatch"])
    hl_b[len(rows_b) - 1] = "#FEF3C7"

    draw_pub_table(ax_b, headers_b, rows_b, cw_b, ca_b, highlight_rows=hl_b, font_pt=5.0)
    panel_title(ax_b, "Individual forensic & molecular QC audit: all 11 fatal COVID-19 cases, 3 resolved survivors, 2 excluded CLL & molecular sex-QC flags", size=5.9, pad=3)

    fig.text(0.020, 0.008,
             "*Sepsis PBISC total = 122 bulk-only samples (121 adult primary + 1 paediatric ICU control PIC06; 0 donor overlap with Real scRNA; 2 HC ambiguous sex calls HC23/HC30 + 1 discordant Res-S34).\n"
             "†C116/C139 excluded for monoclonal CLL. ‡39 COVID-neg donors contribute 38 Real + 19 PBISC neg samples.\n"
             "All 362 matched COVID-19 rows (172 exact 0-day Real bulk+scRNA + 164 bulk→PBISC + 26 scRNA-only) and 255 Sepsis rows are tabulated in PBISC_Metadata_CountVsPct_TISG_Audit_20261121.xlsx.",
             fontsize=5.0, color=PALETTE["gray_2"], va="bottom", ha="left", linespacing=1.18)

    return finalize_figure(fig, FIG_OUT / "SuppFig6_Cohort_Metadata_Tables")


# =============================================================================
# UPGRADED SUPPLEMENTARY FIGURE 9: POPULATION NUMBER (COUNT & CLINICAL CBC) VS
#   PERCENTAGE AUDIT (COVID-43 & SEPSIS-44) + FIRST-DRAW SENS + 14-VAR MATRIX
# =============================================================================
def build_suppfig9_count_pct(d):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.35), facecolor="white")
    A2 = d["A2_COV"].set_index("feature")
    A2S = d["A2_SEP"]

    # a: 43 COVID subsets First-draw mortality AUC (% PBMC INT vs Raw Count INT vs Naive Pooled Raw Count)
    ax_a = fig.add_axes([0.072, 0.775, 0.215, 0.170]); fig_panel_label(fig, 0.012, 0.978, "a")
    ax_a.plot([0.15, 0.85], [0.15, 0.85], color=PALETTE["neutral"], ls="--", lw=0.65)
    ax_a.axhline(0.5, color=PALETTE["neutral"], lw=0.45); ax_a.axvline(0.5, color=PALETTE["neutral"], lw=0.45)
    ax_a.scatter(A2["fd_AUC_pct_INT"], A2["fd_AUC_cnt_INT"], s=11, c=COL["Ext"], alpha=0.85, edgecolors="k", lw=0.25, label="Count INT (r=0.985)")
    ax_a.scatter(A2["fd_AUC_pct_INT"], A2["fd_AUC_rawcnt_pooled"], s=9, marker="^", c=COL["HSPC"], alpha=0.55, lw=0, label="Naive pooled count")
    for fk, xyt in [("T CD8 CTL", (0.28, 0.18)), ("NK", (0.28, 0.26)), ("CD14mono_IFN", (0.52, 0.78))]:
        ax_a.annotate(LBL[fk].replace(" (CD56dim)", ""), xy=(A2.loc[fk, "fd_AUC_pct_INT"], A2.loc[fk, "fd_AUC_cnt_INT"]),
                      xytext=xyt, fontsize=5.0, fontweight="bold", color=COL["Ext"], arrowprops=dict(arrowstyle="-", lw=0.4, color=PALETTE["gray_1"]))
    ax_a.set_xlim(0.16, 0.84); ax_a.set_ylim(0.16, 0.84); ax_a.set_xticks([0.2, 0.5, 0.8]); ax_a.set_yticks([0.2, 0.5, 0.8])
    ax_a.set_xlabel("1st-draw AUC (% PBMC INT)", fontsize=5.3); ax_a.set_ylabel("1st-draw AUC (Count)", fontsize=5.3)
    ax_a.legend(loc="lower right", fontsize=5.0, handletextpad=0.12)
    panel_title(ax_a, "COVID-19: Count vs % (1st-dr death)", size=6.0, pad=3)

    # b: Real sequencing depth n_cells confounder check (226-9,254 cells vs 4,000 fixed PBISC)
    ax_b = fig.add_axes([0.365, 0.775, 0.215, 0.170]); fig_panel_label(fig, 0.305, 0.978, "b")
    meta = d["S_meta"]
    p217 = meta[(meta["condition"] == "positive") & meta["in_primary_scope"].astype(bool) & meta["qc_pass_min200cells"].astype(bool) & (~meta["patient_id"].isin(CLL))]
    r138 = p217[p217["source"] == "Real"]
    rng_n = np.random.default_rng(14)
    for g in range(4):
        vv = r138[r138["temp_severity_ord"] == g]["n_cells"].values / 1000.0
        ax_b.scatter(g + rng_n.uniform(-0.18, 0.18, len(vv)), vv, s=7.5, marker="D", c=COL["Real"], alpha=0.72, lw=0)
        if len(vv):
            ax_b.plot([g - 0.22, g + 0.22], [np.median(vv)] * 2, color="k", lw=1.0, zorder=4)
    ax_b.axhline(4.0, color=COL["New"], ls="--", lw=0.85)
    ax_b.text(3.38, 4.25, "PBISC fixed (4.0k)", ha="right", va="bottom", fontsize=5.0, fontweight="bold", color=COL["New"],
              bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.1))
    ax_b.text(0.03, 0.96, "Real n_cells vs SD-WHO:\nρ=−0.13 (P=0.13)\n1st-dr Death AUC=0.58 (ns)", transform=ax_b.transAxes,
              fontsize=5.0, va="top", color=COL["Real"], bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.15))
    ax_b.set_xticks(range(4)); ax_b.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=5.2)
    ax_b.set_ylabel("Total cells / draw (×1,000)", fontsize=5.3); ax_b.set_ylim(0, 10.2)
    panel_title(ax_b, "Real cell-depth variance (40×)", size=6.0, pad=3)

    # c: Sepsis 44-subset Count vs % concordance (Real r=0.9925, PBISC r=1.000) + NK_c01-GZMK cross-disease depletion
    ax_c = fig.add_axes([0.665, 0.775, 0.295, 0.170]); fig_panel_label(fig, 0.602, 0.978, "c")
    ax_c.plot([-0.6, 0.6], [-0.6, 0.6], color=PALETTE["neutral"], ls="--", lw=0.65)
    ax_c.axhline(0, color=PALETTE["neutral"], lw=0.45); ax_c.axvline(0, color=PALETTE["neutral"], lw=0.45)
    ax_c.scatter(A2S["real_rho_pct"], A2S["real_rho_cnt"], s=11, marker="D", c=COL["Real"], alpha=0.82, edgecolors="white", lw=0.25, label="Real adult (r=0.993)")
    ax_c.scatter(A2S["pbisc_rho_pct"], A2S["pbisc_rho_cnt"], s=9, marker="o", c=COL["New"], alpha=0.55, lw=0, label="PBISC adult (r=1.000)")
    for ct, disp, xyt in [("NK_c01-GZMK", "NK_c01-GZMK", (-0.52, -0.05)), ("B_c07-MZB1", "B_c07-MZB1", (0.22, -0.08)), ("pDC", "pDC", (-0.55, -0.48))]:
        rw = A2S[A2S["celltype"] == ct]
        if len(rw):
            ax_c.annotate(disp, xy=(rw["real_rho_pct"].iloc[0], rw["real_rho_cnt"].iloc[0]), xytext=xyt,
                          fontsize=5.0, fontweight="bold", color=COL["Real"], arrowprops=dict(arrowstyle="-", lw=0.4, color=PALETTE["gray_1"]))
    ax_c.set_xlim(-0.62, 0.78); ax_c.set_ylim(-0.62, 0.78); ax_c.set_xlabel("Sepsis HC→ICU→Sep ρ (% PBMC)", fontsize=5.2); ax_c.set_ylabel("ρ (Raw count)", fontsize=5.2)
    ax_c.legend(loc="lower right", fontsize=5.0, handletextpad=0.12)
    panel_title(ax_c, "Sepsis 44 subsets: Count vs %", size=6.0, pad=3)

    # d: First-draw mortality sensitivity & e: NK-pool decomposition + Sepsis NK_c01-GZMK cross-disease panel
    ax_d = fig.add_axes([0.155, 0.565, 0.255, 0.145]); fig_panel_label(fig, 0.012, 0.735, "d")
    t2_c = d["T2B"].set_index("feature")
    sens_rows = [
        ("T CD8 CTL: 1st dr (66)", 1.0 - t2_c.loc["T CD8 CTL", "fd_AUC"], t2_c.loc["T CD8 CTL", "fd_p"], COL["Ext"], "s"),
        ("T CD8 CTL: adjWHO", 1.0 - t2_c.loc["T CD8 CTL", "fd_adjWHO_AUC"], t2_c.loc["T CD8 CTL", "fd_adjWHO_p"], COL["Ext"], "s"),
        ("T CD8 CTL: ≤7d (53)", 1.0 - t2_c.loc["T CD8 CTL", "fd_le7_AUC"], t2_c.loc["T CD8 CTL", "fd_le7_p"], COL["Ext"], "s"),
        ("T CD8 CTL: WHO 2–3", 1.0 - t2_c.loc["T CD8 CTL", "fd_who23_AUC"], t2_c.loc["T CD8 CTL", "fd_who23_p"], COL["Ext"], "s"),
        ("T CD8 CTL: Real CBC", 1.0 - A2.loc["T CD8 CTL", "fd_AUC_absCBC_Real"], A2.loc["T CD8 CTL", "fd_p_absCBC_Real"], COL["Teal"], "D"),
        ("NK: 1st dr / adjWHO", 1.0 - t2_c.loc["NK", "fd_AUC"], t2_c.loc["NK", "fd_p"], COL["Teal"], "D"),
        ("HSPC: 1st dr (ns)", t2_c.loc["HSPC", "fd_AUC"], t2_c.loc["HSPC", "fd_p"], COL["HSPC"], "o"),
        ("HSPC: Late >7d (50)", t2_c.loc["HSPC", "late_gt7_AUC"], t2_c.loc["HSPC", "late_gt7_p"], COL["HSPC"], "o"),
    ]
    yy_f = np.arange(len(sens_rows))[::-1]
    ax_d.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=0.7); ax_d.axhline(1.5, color=PALETTE["neutral"], ls="--", lw=0.6)
    for yi, (lb, auc_dir, pv, cc, mk) in zip(yy_f, sens_rows):
        ax_d.scatter(auc_dir, yi, s=14, marker=mk, c=cc if pv < 0.05 else "white", edgecolors=cc, lw=0.75, zorder=4)
        ax_d.plot([0.5, auc_dir], [yi, yi], color=cc, lw=0.8, alpha=0.55, zorder=2)
        ax_d.text(0.97, yi, f"{auc_dir:.2f} ({fp_s(pv)})", va="center", fontsize=5.0,
                  fontweight="bold" if pv < 0.05 else "normal", color=cc if pv < 0.05 else PALETTE["gray_2"])
    ax_d.set_yticks(yy_f); ax_d.set_yticklabels([x[0] for x in sens_rows], fontsize=5.0, fontweight="bold")
    ax_d.tick_params(axis="y", length=0); ax_d.spines["left"].set_visible(False)
    ax_d.set_xlim(0.45, 0.95); ax_d.set_xticks([0.5, 0.7, 0.9]); ax_d.set_xlabel("Mortality AUC (aligned >0.5)", fontsize=5.3)
    panel_title(ax_d, "1st-draw & Real CBC sensitivity", size=6.0, pad=3)

    ax_e = fig.add_axes([0.625, 0.565, 0.335, 0.145]); fig_panel_label(fig, 0.552, 0.735, "e")
    NKD = d["NKD"]
    measures = ["CD56bright % PBMC", "CD56bright % of NK", "Total NK % PBMC"]
    for ks, ck in [("Real", "Real"), ("PBISC new", "New")]:
        for ii, mn in enumerate(measures):
            rw = NKD[(NKD["source"] == ks) & (NKD["measure"] == mn)].iloc[0]
            x = ii + (-0.18 if ks == "Real" else 0.18)
            ax_e.bar(x, rw["rho"], width=0.33, color=COL[ck] if rw["p"] < 0.05 else "white", edgecolor=COL[ck], lw=0.8, label=ks if ii == 0 else "")
            ax_e.text(x, rw["rho"] - 0.06 if rw["rho"] < 0 else rw["rho"] + 0.03, "*" if rw["p"] < 0.05 else "ns", ha="center", fontsize=5.2, fontweight="bold")
    ax_e.axhline(0, color="k", lw=0.6); ax_e.set_xticks(range(3)); ax_e.set_xticklabels(["CD56br\n% PBMC", "CD56br\n% of NK", "Total NK\n% PBMC"], fontsize=5.1)
    ax_e.set_ylim(-0.92, 0.28); ax_e.set_ylabel("Severity ρ", fontsize=5.3); ax_e.legend(loc="lower left", fontsize=5.0)
    panel_title(ax_e, "COVID NK pool decomposition", size=6.0, pad=3)

    # f: 14-variant Patient-mean Sensitivity Matrix
    SA = d["SA"]
    FEN = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN"]
    VAR_ORD = SA["variant"].drop_duplicates().tolist()
    ax_s1 = fig.add_axes([0.235, 0.032, 0.340, 0.402]); ax_s2 = fig.add_axes([0.595, 0.032, 0.340, 0.402])
    fig_panel_label(fig, 0.012, 0.515, "f")
    short_f = {f: LBL[f].replace("CD14 mono IFN", "CD14 IFN") for f in FEN}
    for ax, val, pcol, ttl, nrm, show_y in [
        (ax_s1, "rho", "p_sev", "Patient-mean sensitivity: Peak WHO ρ", mpl.colors.Normalize(-0.7, 0.7), True),
        (ax_s2, "AUC", "p_death", "Patient-mean sensitivity: Death AUC", mpl.colors.Normalize(0.15, 0.85), False),
    ]:
        M = SA.pivot_table(index="variant", columns="feature", values=val).reindex(index=VAR_ORD, columns=FEN)
        Pm = SA.pivot_table(index="variant", columns="feature", values=pcol).reindex(index=VAR_ORD, columns=FEN)
        ax.imshow(M.values, cmap="RdBu_r", norm=nrm, aspect="auto")
        for ii in range(M.shape[0]):
            for jj in range(M.shape[1]):
                v = M.values[ii, jj]
                if np.isnan(v):
                    ax.text(jj, ii, "–", ha="center", va="center", fontsize=5.1, color=PALETTE["gray_1"])
                else:
                    sig = Pm.values[ii, jj] < 0.05
                    ax.text(jj, ii, f"{v:.2f}", ha="center", va="center", fontsize=5.1, fontweight="bold" if sig else "normal",
                            color=("white" if abs(nrm(v) - 0.5) > 0.32 else "black") if sig else PALETTE["gray_1"])
        ax.axvline(3.5, color="k", lw=0.8)
        ax.set_xticks(range(len(FEN))); ax.set_xticklabels([short_f[f] for f in FEN], rotation=35, ha="left", fontsize=5.1, fontweight="bold")
        ax.xaxis.tick_top(); ax.set_yticks(range(len(VAR_ORD)))
        sh_v = [v.replace("Published 11d", "Pub11").replace("Direct COVID only", "Direct").replace("Zenodo w1+w2 raw", "Zenodo").replace("G  combined: lineage %, no outliers, adj.", "G combined adj.").replace("combined: lineage %, no outliers, adj.", "combined adj.") for v in VAR_ORD]
        ax.set_yticklabels(sh_v if show_y else [], fontsize=5.1); ax.tick_params(length=0)
        for sp_ in ax.spines.values(): sp_.set_visible(False)
        ax.text(0.5, 1.19, ttl, transform=ax.transAxes, ha="center", va="bottom", fontsize=5.9, fontweight="bold", color=COL["Real"])

    return finalize_figure(fig, FIG_OUT / "SuppFig9_Count_vs_Pct_CBC_and_14VarMatrix")


# =============================================================================
# UPGRADED SUPPLEMENTARY FIGURE 10: COMPLETE 13-SUBSET T-CELL ISG AUDIT +
#   CROSS-DISEASE 19-ISG DEG + HSPC SUPPORTING + STEPHENSON + POWER CURVES
# =============================================================================
def build_suppfig10_tisg_hspc(d):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.25), facecolor="white")
    t6_c = d["T6_ISG"][d["T6_ISG"]["set"] == "Combined"].set_index("subset")
    t8b = d["T8B_WIN"]
    tisg13 = d["TISG13"].set_index("subset")
    bpz = d["VER_JSON"]["figure6"]["official_R432_BayesPrism_Z_shared27"]

    # ---- Row 1: a (Complete 13-subset T-cell ISG audit matrix) & b (PBISC vs BayesPrism Z + Cross-disease 19-ISG DEG) ----
    ax_a = fig.add_axes([0.132, 0.675, 0.435, 0.255]); fig_panel_label(fig, 0.012, 0.975, "a")
    SUB12 = ["All T cells", "All CD8 T", "All CD4 T", "T CD8 CTL", "T CD4 CTL", "T CD4 CM/EM/CTL",
             "T CD8 CM", "T CD8 Naive", "T CD4 Helper", "T CD4 Naive", "T g/d", "MAIT", "T reg"]
    SUB12 = [s for s in SUB12 if s in t6_c.index]
    COLS_T = ["Sh27 ρ", "PM Peak", "≤7d β", "SD+wk", "1st AUC", ">7d AUC"]
    mat_t = np.zeros((len(SUB12), 6)); sig_t = np.zeros((len(SUB12), 6), bool)
    for ii, s in enumerate(SUB12):
        rc = t6_c.loc[s]
        e7 = t8b[(t8b["subset"] == s) & (t8b["window"] == "≤7 d")].iloc[0]
        l7 = t8b[(t8b["subset"] == s) & (t8b["window"] == ">7 d")].iloc[0]
        sh = float(tisg13.loc[s, "sh_rho"]) if s in tisg13.index else np.nan
        mat_t[ii] = [sh, rc["pm_peakWHO_rho"], e7["gee_sev_beta"], rc["sd_gee_wk_beta"], 2*(rc["fd_AUC"]-0.5), 2*(l7["pm_AUC"]-0.5)]
        sig_t[ii] = [sh > 0.5, rc["pm_q13_sev"] < 0.05, e7["gee_sev_p"] < 0.05, rc["sd_gee_wk_p"] < 0.05, rc["fd_p"] < 0.10, l7["pm_p_death"] < 0.10]
    ax_a.imshow(mat_t, cmap="RdBu_r", norm=mpl.colors.Normalize(-0.85, 0.85), aspect="auto")
    for ii in range(len(SUB12)):
        for jj in range(6):
            vv = mat_t[ii, jj]
            disp = f"{0.5 + vv/2:.2f}" if jj >= 4 else f"{vv:+.2f}"
            ax_a.text(jj, ii, disp, ha="center", va="center", fontsize=5.0,
                      fontweight="bold" if sig_t[ii, jj] else "normal",
                      color="white" if abs(vv) > 0.52 else COL["Real"])
    ax_a.set_yticks(range(len(SUB12))); ax_a.set_yticklabels(SUB12, fontsize=5.0, fontweight="bold")
    ax_a.set_xticks(range(6)); ax_a.set_xticklabels(COLS_T, fontsize=5.0, fontweight="bold")
    ax_a.xaxis.tick_top(); ax_a.tick_params(length=0)
    for sp_ in ax_a.spines.values(): sp_.set_visible(False)
    panel_title(ax_a, "All 13 T-cell compartments: ISG reconstruction, severity & mortality", size=6.0, pad=13)

    # b1 (PBISC vs BayesPrism Z) & b2 (Cross-disease COVID vs Sepsis 19-ISG T-cell log2FC)
    ax_b1 = fig.add_axes([0.648, 0.832, 0.325, 0.105]); fig_panel_label(fig, 0.585, 0.975, "b")
    pb_ctl, pb_cm = float(tisg13.loc["T CD8 CTL", "sh_rho"]), float(tisg13.loc["T CD8 CM", "sh_rho"])
    bp_ctl, bp_cm = bpz["T CD8 CTL_ISG"]["official_R_BayesPrism_Z_rho"], bpz["T CD8 CM_ISG"]["official_R_BayesPrism_Z_rho"]
    yb1 = np.arange(2)[::-1]; h1 = 0.34
    ax_b1.barh(yb1 + h1 / 2, [pb_ctl, pb_cm], h1, color=COL["Ext"], edgecolor="k", lw=0.55, label="PBISC scRNA")
    ax_b1.barh(yb1 - h1 / 2, [bp_ctl, bp_cm], h1, color=PALETTE["red_2"], edgecolor="k", lw=0.55, hatch="///", label="BayesPrism Z")
    for yi, vp, vb in zip(yb1, [pb_ctl, pb_cm], [bp_ctl, bp_cm]):
        ax_b1.text(vp + 0.02, yi + h1 / 2, f"{vp:.2f}", va="center", fontsize=5.1, fontweight="bold", color=COL["Ext"])
        ax_b1.text(vb + 0.02, yi - h1 / 2, f"{vb:.2f}", va="center", fontsize=5.1, color=COL["Real"])
    ax_b1.set_yticks(yb1); ax_b1.set_yticklabels(["CD8 CTL\nISG", "CD8 CM\nISG"], fontsize=5.1, fontweight="bold")
    ax_b1.set_xlim(0, 1.02); ax_b1.set_xlabel("Shared-27 ρ vs Real", fontsize=5.2)
    ax_b1.legend(loc="lower right", fontsize=5.0, handlelength=0.85)
    panel_title(ax_b1, "PBISC vs official BayesPrism Z", size=5.9, pad=3)

    ax_b2 = fig.add_axes([0.648, 0.675, 0.325, 0.105])
    cdeg = d["A3_CROSS_DEG"]
    G_ISG_SHOW = ["IFI44L", "IFIT1", "IFIT3", "MX1", "OAS1", "ISG15", "STAT1", "IRF7"]
    x_g = np.arange(len(G_ISG_SHOW))
    for coh, mk, dy in [("covid", "s", 0.0), ("sepsis", "o", 0.0)]:
        for src, dx, cc in [("Real", -0.18, COL["Real"]), ("PBISC", +0.18, COL["Ext"])]:
            sub = cdeg[(cdeg["cohort"] == coh) & (cdeg["source"] == src)].set_index("gene").reindex(G_ISG_SHOW)
            ax_b2.scatter(x_g + dx, sub["logFC"], s=13 if coh == "covid" else 10, marker=mk,
                          c=cc if coh == "covid" else "white", edgecolors=cc, lw=0.75, zorder=4,
                          label=f"{'COVID' if coh=='covid' else 'Sep'} {src}")
    ax_b2.axhline(0, color="k", lw=0.6)
    ax_b2.set_xticks(x_g); ax_b2.set_xticklabels(G_ISG_SHOW, rotation=32, ha="right", fontsize=5.0, fontstyle="italic", fontweight="bold")
    ax_b2.set_ylabel("T log2FC", fontsize=5.1); ax_b2.set_ylim(-1.6, 4.3)
    ax_b2.legend(loc="upper right", ncol=2, fontsize=5.0, handletextpad=0.05, columnspacing=0.25)
    panel_title(ax_b2, "COVID induction vs Sepsis suppression", size=5.8, pad=3)

    # ---- Row 2: c (HSPC Ery5, Alarmin, New-29 vs raw bulk, Paired n=43) ----
    H149, H61, T5 = d["H149_Z"], d["H61_FD"], d["T5_HSPC"]
    t5_c = T5[T5["set"] == "Combined"].set_index("feature")
    ax_c1 = fig.add_axes([0.072, 0.375, 0.132, 0.195]); ax_c2 = fig.add_axes([0.218, 0.375, 0.095, 0.195])
    fig_panel_label(fig, 0.012, 0.605, "c")
    rng_b = np.random.default_rng(28)
    for g in range(4):
        for arm_v, ks, dx in [("Real", "Real", -0.17), ("PBISC new", "New", +0.17)]:
            vals = H149[(H149["arm"] == arm_v) & (H149["temp_severity_ord"] == g)]["Ery5__z149"].values
            ax_c1.scatter(g + dx + rng_b.uniform(-0.07, 0.07, len(vals)), vals, s=6, marker=MK[ks], c=COL[ks], alpha=0.75, edgecolors="white", lw=0.15)
            if len(vals): ax_c1.plot([g + dx - 0.13, g + dx + 0.13], [np.median(vals)] * 2, color="k", lw=0.9)
    ax_c1.set_xticks(range(4)); ax_c1.set_xticklabels(["0", "1", "2", "3"], fontsize=5.1); ax_c1.set_ylim(-2.5, 4.1)
    ax_c1.set_xlabel("Same-day WHO", fontsize=5.2); ax_c1.set_ylabel("Normal score", fontsize=5.3)
    ax_c1.text(0.04, 0.96, f"β=+0.26\nq16=1.7e-3*", transform=ax_c1.transAxes, fontsize=5.0, fontweight="bold", color=COL["HSPC"], va="top")
    panel_title(ax_c1, "HSPC Ery5 (149 dr)", size=5.9, pad=3)

    for gf in [0, 1]:
        for arm_v, ks, dx in [("Real", "Real", -0.17), ("PBISC new", "New", +0.17)]:
            vals = H61[(H61["arm"] == arm_v) & (H61["fatal_pub11"] == gf)]["Alarmin_Inflam__zfd"].values
            ax_c2.scatter(gf + dx + rng_b.uniform(-0.07, 0.07, len(vals)), vals, s=7.5, marker=MK[ks], c=COL[ks], edgecolors="white", lw=0.2)
            if len(vals): ax_c2.plot([gf + dx - 0.13, gf + dx + 0.13], [np.median(vals)] * 2, color="k", lw=0.9)
    ax_c2.set_xticks([0, 1]); ax_c2.set_xticklabels(["Surv", "Died"], fontsize=5.1); ax_c2.set_ylim(-2.5, 4.1); ax_c2.set_yticks([])
    ax_c2.spines["left"].set_visible(False)
    ax_c2.text(0.04, 0.96, f"AUC=0.72\nadjP=.022†", transform=ax_c2.transAxes, fontsize=5.0, fontweight="bold", color=COL["HSPC"], va="top")
    panel_title(ax_c2, "Alarmin (1st)", size=5.9, pad=3)

    # d: New-29 PBISC HSPC vs raw bulk HSPC genes + Paired n=43 HSPC early vs late windows
    ax_d1 = fig.add_axes([0.445, 0.375, 0.185, 0.195]); fig_panel_label(fig, 0.340, 0.605, "d")
    SG = d["SG_HSPC"]
    b8_h = d["BENCH10_D28"][(d["BENCH10_D28"]["method"] == "Bulk-8-markers") & (d["BENCH10_D28"]["feature"] == "HSPC") & (d["BENCH10_D28"]["endpoint"] == "death")].iloc[0]
    sg_ext = pd.concat([SG.iloc[[0]], pd.DataFrame([{"label": "Bulk 6-gene", "is_pbisc": False, "auc": b8_h["new29_eff"], "lo": 0.14, "hi": 0.63, "p": b8_h["new29_p"]}]), SG.iloc[1:]], ignore_index=True)
    yd = np.arange(len(sg_ext))[::-1]; ax_d1.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=0.7)
    for yi, (_, r_sg) in zip(yd, sg_ext.iterrows()):
        isp = bool(r_sg["is_pbisc"]); cc = COL["HSPC"] if isp else PALETTE["gray_2"]
        ax_d1.plot([r_sg["lo"], r_sg["hi"]], [yi, yi], color=cc, lw=1.0 if isp else 0.7)
        ax_d1.scatter(r_sg["auc"], yi, s=14 if isp else 10, marker="s" if isp else "o", c=cc if r_sg["p"] < 0.05 else "white", edgecolors=cc, lw=0.7, zorder=4)
    ax_d1.set_yticks(yd); ax_d1.set_yticklabels([s.replace("PBISC HSPC (%)", "PBISC HSPC").replace(" (CD133)", "") for s in sg_ext["label"]], fontsize=5.0)
    ax_d1.get_yticklabels()[0].set_fontweight("bold"); ax_d1.get_yticklabels()[0].set_color(COL["HSPC"])
    ax_d1.tick_params(axis="y", length=0); ax_d1.spines["left"].set_visible(False)
    ax_d1.set_xlim(0.05, 1.04); ax_d1.set_xlabel("New-29 HSPC AUC", fontsize=5.2)
    panel_title(ax_d1, "PBISC HSPC vs raw bulk", size=5.9, pad=3)

    ax_d2 = fig.add_axes([0.742, 0.375, 0.165, 0.195])
    P43_WIN = d["P43_WIN"].set_index("window")
    WIN_L = [("early", "Early ≤7d\n(n=43)"), ("late", "Late >7d\n(n=43)"), ("delta", "Late − Early\nΔ (n=43)")]
    yy_w = np.arange(3)[::-1]; ax_d2.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=0.7)
    for yi, (wk, _) in zip(yy_w, WIN_L):
        rw = P43_WIN.loc[wk]; cc = COL["HSPC"] if rw["P"] < 0.05 else PALETTE["gray_2"]
        ax_d2.plot([rw["CI_lo"], rw["CI_hi"]], [yi, yi], color=cc, lw=0.95)
        ax_d2.scatter(rw["AUC"], yi, s=16, marker="s" if wk != "early" else "D", c=cc if rw["P"] < 0.05 else "white", edgecolors=cc, lw=0.75, zorder=4)
        ax_d2.text(1.03, yi, f"{rw['AUC']:.2f}\n({fp_s(rw['P'])})", va="center", fontsize=5.0, fontweight="bold" if rw["P"] < 0.05 else "normal", color=cc, linespacing=0.95)
    ax_d2.set_yticks(yy_w); ax_d2.set_yticklabels([x[1] for x in WIN_L], fontsize=5.0, fontweight="bold")
    ax_d2.tick_params(axis="y", length=0); ax_d2.spines["left"].set_visible(False)
    ax_d2.set_xlim(0.10, 1.00); ax_d2.set_ylim(-0.55, 2.55); ax_d2.set_xlabel("Paired AUC (5d)", fontsize=5.2)
    panel_title(ax_d2, "Paired n=43 HSPC shift", size=5.9, pad=3)

    # ---- Row 3: e (Patient-mean sensitivity power curves) & f (Stephenson 2021 B ASC isotypes) ----
    PW_UNI, PUB_COH, ESKD_MATCH = d["PW_UNI"], d["PUB_COH"], d["ESKD_MATCH"]
    ax_e = fig.add_axes([0.088, 0.055, 0.415, 0.235]); fig_panel_label(fig, 0.012, 0.320, "e")
    for fk_p, ep_p, lbl_p, cc_p, has_r in [
        ("B_ASC_IgG", "sev", "B ASC IgG (s)", COL["Ext"], True),
        ("B_ASC_IgA", "sev", "B ASC IgA (s)", COL["New"], True),
        ("HSPC", "sev", "HSPC (s)", PALETTE["teal"], True),
        ("T CD8 CTL", "death", "T CD8 CTL (d)", PALETTE["gray_3"], False),
        ("HSPC", "death", "HSPC (d)", COL["HSPC"], True),
    ]:
        sub = PW_UNI[(PW_UNI["feature"] == fk_p) & (PW_UNI["outcome"] == ep_p)].sort_values("n_real")
        if has_r:
            ax_e.plot(sub["n_real"], sub["real_nom"], color=cc_p, lw=0.72, ls="--", marker="D", ms=2.4, mfc="white")
        ax_e.plot(sub["n_real"], sub["plus_new29_nom"], color=cc_p, lw=1.0, ls="-", marker="s", ms=2.8, label=f"{lbl_p} (+New29)")
    d0 = PW_UNI[(PW_UNI["feature"] == "HSPC") & (PW_UNI["outcome"] == "death")].sort_values("n_real")
    shuf_avg = PW_UNI.groupby("n_real")["plus_shuf_nom"].mean().loc[d0["n_real"]].values
    ax_e.plot(d0["n_real"], shuf_avg, color=PALETTE["gray_1"], lw=0.85, ls=":", label="+New29 shuffled")
    ax_e.axhline(0.8, color=PALETTE["gray_1"], ls=":", lw=0.65)
    ax_e.set_ylim(-0.03, 1.08); ax_e.set_xlim(8.5, 36.5); ax_e.set_xticks([10, 15, 20, 25, 30, 35])
    ax_e.set_xlabel("Sampled Real patients (+29 PBISC new)", fontsize=5.3); ax_e.set_ylabel("Re-detection (P < 0.05)", fontsize=5.3)
    ax_e.legend(loc="center left", bbox_to_anchor=(0.02, 0.56), fontsize=5.0, handlelength=1.0, labelspacing=0.14)
    panel_title(ax_e, "Patient-mean sensitivity power curves", size=6.0, pad=3)

    ax_f = fig.add_axes([0.665, 0.055, 0.285, 0.225]); fig_panel_label(fig, 0.545, 0.320, "f")
    ISO = [("B_ASC_dividing", "Plasmablast div.", "B_ASC_dividing"), ("B_ASC_IgG", "Plasma IgG", "B_ASC_IgG"),
           ("B_ASC_IgA", "Plasma IgA", "B_ASC_IgA"), ("B_ASC_IgM", "Plasma IgM", "B_ASC_IgM"), ("B_ASC_total", "B ASC total", "B_ASC_total")]
    st_b = PUB_COH[PUB_COH["cohort"] == "Stephenson 2021"].set_index("feature")
    yy_b = np.arange(len(ISO))[::-1]; ax_f.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=0.7)
    for yi, (fk, flb, efk) in zip(yy_b, ISO):
        r, em = st_b.loc[fk], ESKD_MATCH.loc[efk]
        ax_f.plot([em["lo"], em["hi"]], [yi + 0.16, yi + 0.16], color=COL["Ext"], lw=0.8)
        ax_f.scatter(em["rho"], yi + 0.16, s=12, marker="s", c=COL["Ext"], zorder=3, label="ESKD (66)" if yi == yy_b[0] else "")
        ax_f.plot([r["rho_lo"], r["rho_hi"]], [yi - 0.16, yi - 0.16], color=COL["Teal"], lw=0.8)
        ax_f.scatter(r["rho"], yi - 0.16, s=13, marker="D", c="white" if r["p_sev"] >= 0.05 else COL["Teal"], edgecolors=COL["Teal"], lw=0.7, zorder=4, label="Stephenson (74)" if yi == yy_b[0] else "")
        ax_f.text(0.76, yi, f"OR={r['OR_glm']:.2f}*", va="center", fontsize=5.0, color=COL["Ext"], fontweight="bold")
    ax_f.set_yticks(yy_b); ax_f.set_yticklabels([x[1] for x in ISO], fontsize=5.1, fontweight="bold")
    ax_f.tick_params(axis="y", length=0); ax_f.spines["left"].set_visible(False)
    ax_f.set_xlim(-0.20, 1.20); ax_f.set_xticks([0.0, 0.4, 0.7]); ax_f.set_xlabel("Severity Spearman ρ", fontsize=5.3)
    ax_f.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=5.0, handletextpad=0.10, columnspacing=0.35)
    panel_title(ax_f, "Stephenson 2021 B ASC isotypes", size=6.0, pad=12)

    return finalize_figure(fig, FIG_OUT / "SuppFig10_Deep_TISG_HSPC_Stephenson_Power")


def write_stdlib_xlsx(sheets: list, out_path: pathlib.Path):
    """Pure standard-library OOXML .xlsx writer (zero openpyxl dependency)."""
    def col_letter(idx: int) -> str:
        res = ""
        n = idx + 1
        while n > 0:
            n, r = divmod(n - 1, 26)
            res = chr(65 + r) + res
        return res

    sheet_xmls = []
    for name, df in sheets:
        rows_xml = []
        # Header row
        h_cells = [f'<c r="{col_letter(j)}1" t="inlineStr"><is><t>{escape(str(c))}</t></is></c>' for j, c in enumerate(df.columns)]
        rows_xml.append(f'<row r="1">{"".join(h_cells)}</row>')
        for i, row in enumerate(df.itertuples(index=False), start=2):
            c_xml = []
            for j, val in enumerate(row):
                if pd.isna(val):
                    continue
                ref = f"{col_letter(j)}{i}"
                if isinstance(val, (int, float, np.integer, np.floating)) and not isinstance(val, bool) and np.isfinite(val):
                    c_xml.append(f'<c r="{ref}"><v>{val}</v></c>')
                else:
                    c_xml.append(f'<c r="{ref}" t="inlineStr"><is><t>{escape(str(val))}</t></is></c>')
            rows_xml.append(f'<row r="{i}">{"".join(c_xml)}</row>')
        sheet_xmls.append((name[:31], (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
            + "".join(rows_xml) + "</sheetData></worksheet>"
        )))

    wb_sheets = "".join([f'<sheet name="{escape(nm)}" sheetId="{i+1}" r:id="rId{i+1}"/>' for i, (nm, _) in enumerate(sheet_xmls)])
    wb_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets>{wb_sheets}</sheets></workbook>'
    )
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join([f'<Relationship Id="rId{i+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i+1}.xml"/>'
                   for i in range(len(sheet_xmls))])
        + "</Relationships>"
    )
    ct_overrides = "".join([f'<Override PartName="/xl/worksheets/sheet{i+1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                            for i in range(len(sheet_xmls))])
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        f'{ct_overrides}</Types>'
    )
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        '</Relationships>'
    )
    with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("xl/workbook.xml", wb_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        for i, (_, sxml) in enumerate(sheet_xmls):
            zf.writestr(f"xl/worksheets/sheet{i+1}.xml", sxml)


def build_excel_audit_bundle(d):
    sheets = [
        ("1a_COVID_107Donors_Match", pd.read_csv(D28_DATA / "table_s3_all_107_donors_overview.csv")),
        ("1b_COVID_362Samples_1to1", pd.read_csv(D28_DATA / "table_s2_sample_1to1_mapping_362.csv")),
        ("1c_Sepsis_255Samples_Match", d["A1_SEP"]),
        ("2a_COVID43_Count_vs_Pct_CBC", d["A2_COV"]),
        ("2b_Sepsis44_Count_vs_Pct", d["A2_SEP"]),
        ("3a_COVID_TISG_13Subsets", d["T6_ISG"]),
        ("3b_COVID_TISG_EarlyVsLate", d["T8B_WIN"]),
        ("3c_COVID_63Exact1stDr_ISG", d["A3_FD66"]),
        ("3d_COVID_TISG_50x5_CV", d["A3_CV"]),
        ("3e_Sepsis_PBISC_TISG", d["A3_SEP_ISG"]),
        ("3f_CrossDisease_19ISG_DEG", d["A3_CROSS_DEG"]),
    ]
    out_xlsx = AUD3 / "PBISC_Metadata_CountVsPct_TISG_Audit_20261121.xlsx"
    write_stdlib_xlsx(sheets, out_xlsx)
    print(f"[OK] Built 11-sheet Excel audit workbook: {out_xlsx.name} ({out_xlsx.stat().st_size:,} bytes)")


def cleanup_old_stems():
    obsolete = [
        "Figure6_COVID_Rare_Population_Reproduction",
        "Figure7_COVID_Clinical_Discovery_Severity_Mortality",
        "SuppFig9_SameDay_Anchors_NK_FirstDrawSens_14VarMatrix",
        "SuppFig10_Power_Curves_Stephenson_HSPC_TISG_Details",
    ]
    for st in obsolete:
        for ext in ("png", "pdf", "svg"):
            p = FIG_OUT / f"{st}.{ext}"
            if p.exists():
                p.unlink()


def main():
    apply_publication_style(FigureStyle())
    d = load_data()
    builders = {
        "fig6": build_figure6,
        "fig7": build_figure7,
        "sfig6": build_suppfig6_metadata,
        "sfig9": build_suppfig9_count_pct,
        "sfig10": build_suppfig10_tisg_hspc,
    }
    wanted = [a for a in sys.argv[1:] if a in builders] or list(builders)
    audit_path = OUT_ROOT / "AUDIT_FIG6_7_SFIG6_11.json"
    prev = {a["figure"]: a for a in json.loads(audit_path.read_text(encoding="utf-8"))} if audit_path.exists() else {}
    for k_drop in ("Figure6_COVID_Rare_Population_Reproduction", "Figure7_COVID_Clinical_Discovery_Severity_Mortality",
                   "SuppFig9_SameDay_Anchors_NK_FirstDrawSens_14VarMatrix", "SuppFig10_Power_Curves_Stephenson_HSPC_TISG_Details"):
        prev.pop(k_drop, None)
    for key in wanted:
        saved, aud = builders[key](d)
        aud = {k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in aud.items()}
        prev[aud["figure"]] = aud
        print(f"[OK] {aud['figure']}: {aud['width_mm']}x{aud['height_mm']} mm | min_pt={aud['min_font_pt']} | below5={aud['below_min_pt_count']} | oob={aud['out_of_bounds_count']}")
    cleanup_old_stems()
    build_excel_audit_bundle(d)
    order = [
        "Figure6_Rare_Population_Sepsis_COVID_Reproduction",
        "Figure7_COVID_Clinical_Metadata_CountPct_TISG",
        "SuppFig6_Cohort_Metadata_Tables",
        "SuppFig8_Audit_CLL_Funnel_Directional_Matrix",
        "SuppFig9_Count_vs_Pct_CBC_and_14VarMatrix",
        "SuppFig10_Deep_TISG_HSPC_Stephenson_Power",
    ]
    audit_path.write_text(json.dumps([prev[k] for k in order if k in prev], indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
