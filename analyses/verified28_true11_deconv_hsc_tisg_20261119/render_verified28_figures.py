"""Stage 2 Standalone Figure Renderer for Folder 28 (Verified True-11 Mortality,
10-Method Deconvolution incl. BayesPrism/CIBERSORTx/MuSiC, HSPC Transcriptome, and T-cell ISG).

100% self-contained `--render-only` script: reads ONLY local CSVs from `--data-dir`
(zero external reads outside `<data-dir>`), renders Main Figures 7 & 8 and Supplementary
Figures 1-5 (`PNG`, `PDF`, `SVG` at 300 DPI) into `--output-dir`, merges the 7-page
vector PDF bundle, and writes `RUN_RECEIPT.json`.
"""
import argparse
import hashlib
import json
import pathlib
import sys
import warnings

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

from style import CT_COL, FigureStyle, PALETTE, apply_publication_style, finalize_figure

DEFAULT_OUT = SCRIPT_DIR
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
    ScoreGood="#1B7A43",
    ScoreMid="#B86B11",
    ScoreBad="#B64342",
    Good="#1B7A43",
    Mid="#B86B11",
    Bad="#B64342",
    CardBlue="#DCE6F2",
    CardGray="#EFEFEF",
    CardTeal="#DCF0F2",
    CardRed="#FADBD8",
    CardGreen=PALETTE["green_1"],
    SoftRed="#FADBD8",
    SoftBlue="#DCE6F2",
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
    "CD14mono_IFN": "CD14 mono IFN",
    "CD14mono_dex": "CD14 mono dex",
    "CD16mono": "CD16 mono",
    "CD16mono|lin": "CD16 mono (% Mono)",
    "CD16mono_C1": "CD16 mono C1",
    "Int.mono": "Int. mono",
    "B_CD11c": "B CD11c",
    "DC1": "DC1",
    "DC2": "DC2",
    "DC3": "DC3",
    "NK": "NK CD56dim",
    "NKT": "NKT",
    "T CD4 CTL": "T CD4 CTL",
    "T g/d": "T γδ",
    "MAIT": "MAIT",
}

ANCHOR_SET = {"NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgM"}
FOCUS_SET = {"B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC"}


def fp(p):
    if p < 1e-4:
        return "P < 1e-4"
    if p < 0.001:
        return f"P = {p:.1e}".replace("e-0", "e-")
    if p < 0.1:
        return f"P = {p:.3f}"
    return f"P = {p:.2f}"


def fp_short(p):
    if p < 1e-4:
        return "<1e-4"
    if p < 0.001:
        return f"{p:.1e}".replace("e-0", "e-")
    if p < 0.1:
        return f"{p:.3f}"
    return f"{p:.2f}"


def p_eq_s(p):
    if p < 1e-4:
        return "P<1e-4"
    if p < 0.001:
        return f"P={p:.1e}".replace("e-0", "e-")
    if p < 0.1:
        return f"P={p:.3f}"
    return f"P={p:.2f}"


def auc_mw(y, f):
    y_a, f_a = np.asarray(y, float), np.asarray(f, int)
    m = np.isfinite(y_a)
    y_a, f_a = y_a[m], f_a[m]
    n1, n0 = int((f_a == 1).sum()), int((f_a == 0).sum())
    if n1 < 1 or n0 < 1 or np.nanstd(y_a) == 0:
        return 0.5, 1.0
    u = stats.mannwhitneyu(y_a[f_a == 1], y_a[f_a == 0])
    return float(u.statistic / (n1 * n0)), float(u.pvalue)


def n80_interp(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if y[0] >= 0.8:
        return x[0], True
    k = np.where(y >= 0.8)[0]
    if not len(k):
        return np.nan, False
    i = k[0]
    return x[i - 1] + (0.8 - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]), False


def rbox(ax, x, y, w, h, fc, ec, lw=1.6, r=2.0, pad=0.20, alpha=1.0, ls="-", zorder=2):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad={pad},rounding_size={r}",
        fc=fc, ec=ec, lw=lw, ls=ls, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(p)
    return p


def plot_ci_clipped(ax, y, lo, hi, color, xlim, lw=1.8, ls="-", alpha=1.0, zorder=3):
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return
    lo_c = max(xlim[0] + 0.008, lo)
    hi_c = min(xlim[1] - 0.008, hi)
    ax.plot([lo_c, hi_c], [y, y], color=color, lw=lw, ls=ls, alpha=alpha, zorder=zorder)
    if lo < xlim[0] + 0.008:
        ax.annotate("", xy=(xlim[0] + 0.004, y), xytext=(xlim[0] + 0.052, y),
                    arrowprops=dict(arrowstyle="-|>", lw=1.0, color=color, alpha=alpha, mutation_scale=7))
    if hi > xlim[1] - 0.008:
        ax.annotate("", xy=(xlim[1] - 0.004, y), xytext=(xlim[1] - 0.052, y),
                    arrowprops=dict(arrowstyle="-|>", lw=1.0, color=color, alpha=alpha, mutation_scale=7))


def plot_dl_and_mkh(ax, y, dl_lo, dl_hi, mkh_lo, mkh_hi, k, color, xlim, lw_dl=2.0):
    plot_ci_clipped(ax, y, dl_lo, dl_hi, color, xlim, lw=lw_dl, ls="-", alpha=1.0, zorder=3)
    if k >= 2 and np.isfinite(mkh_lo) and np.isfinite(mkh_hi):
        plot_ci_clipped(ax, y, mkh_lo, mkh_hi, color, xlim, lw=1.15, ls="--", alpha=0.68, zorder=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=pathlib.Path, default=DEFAULT_OUT / "data")
    ap.add_argument("--output-dir", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    apply_publication_style(FigureStyle(font_size=13, axes_linewidth=1.8))

    # -----------------------------------------------------------------
    # Load 100% of inputs exclusively from local `data_dir`
    # -----------------------------------------------------------------
    S_meta = pd.read_csv(data_dir / "eskd_sample_metadata_verified.csv")
    pos_217 = pd.read_csv(data_dir / "eskd_longitudinal_217_draws.csv")

    pos_samps = S_meta[
        (S_meta["condition"] == "positive")
        & S_meta["in_primary_scope"].astype(bool)
        & S_meta["qc_pass_min200cells"].astype(bool)
        & (~S_meta["patient_id"].isin(CLL))
    ].copy()
    r_samps = pos_samps[pos_samps["source"] == "Real"].copy()
    b_all_samps = pos_samps[pos_samps["source"] == "PBISC"].copy()
    b_new_samps = b_all_samps[~b_all_samps["patient_id"].isin(set(r_samps["patient_id"]))].copy()

    CLL_BULK = pd.read_csv(data_dir / "fig7_cll_molecular_screen.csv")
    SCAN = pd.read_csv(data_dir / "fig7_scan_43fine_true11.csv")
    SCAN_NO = SCAN[SCAN["CLL_excluded"] == True].copy()
    SEL = pd.read_csv(data_dir / "fig7_selection_funnel_43fine_true11.csv")
    SC = pd.read_csv(data_dir / "fig7_credibility_scorecard_true11.csv")
    CI = pd.read_csv(data_dir / "fig7_cohort_effect_sizes_stratified_ci_true11.csv")
    HS = pd.read_csv(data_dir / "fig7_hspc_death_strata_true11.csv")
    HS_COND = pd.read_csv(data_dir / "fig7_exact_source_peakWHO_cond_true11.csv").set_index("feature")
    HS_L1O = pd.read_csv(data_dir / "fig7_hspc_leave_one_death_true11.csv")
    PW_UNI = pd.read_csv(data_dir / "fig7_s3_unified_redetection_true11.csv")
    MORT_AUD = pd.read_csv(data_dir / "metadata_mortality_audit_discrepancies.csv")

    MOD = pd.read_csv(data_dir / "suppS2_modality_scale_summary.csv")
    SH27 = pd.read_csv(data_dir / "suppS2_shared27_paired_pct.csv")
    PT_NS = pd.read_csv(data_dir / "fig7_patient_normal_scores_noCLL.csv", index_col=0)
    PT_RAW = pd.read_csv(data_dir / "fig7_patient_raw_pct_noCLL.csv", index_col=0)
    SA = pd.read_csv(data_dir / "suppS3_sensitivity_14variants_true11.csv")
    NKD = pd.read_csv(data_dir / "suppS2_nk_decomposition.csv")
    PR_REALB = pd.read_csv(data_dir / "suppS2_realB_transfer_precision_recall.csv").set_index("label")

    PUB_COH_CORR = pd.read_csv(data_dir / "suppS4_public_cohort_reconciled.csv")
    PUB_META_CORR = pd.read_csv(data_dir / "suppS4_public_meta_reconciled.csv")
    PUB_HW_CORR = pd.read_csv(data_dir / "suppS4_public_hspc_windows_reconciled.csv")
    ESKD_MATCH_DF = pd.read_csv(data_dir / "suppS4_eskd_matched_severity_ci.csv").set_index("feature")

    BENCH10 = pd.read_csv(data_dir / "fig8_deconv_10methods_true11.csv")
    NW29_CI = pd.read_csv(data_dir / "fig8_new29_hspc_10methods_CI.csv").set_index("method")
    UNB9 = pd.read_csv(data_dir / "suppS5_unbiased_deconv_9methods_summary.csv")
    SG_HSPC = pd.read_csv(data_dir / "suppS5_single_gene_hspc_new29.csv")

    HSPC_MOD16 = pd.read_csv(data_dir / "fig8_hspc_modules_16unique_true11.csv")
    HSPC_PT62 = pd.read_csv(data_dir / "fig8_hspc_patients_62_true11.csv")
    HSPC_G20 = pd.read_csv(data_dir / "suppS5_hspc_20gene_true11.csv")

    TISG13 = pd.read_csv(data_dir / "fig8_tcell_isg_13subsets_true11.csv")
    TISG_GEE = pd.read_csv(data_dir / "suppS5_tcell_isg_gee_true11.csv")
    CD8_ISG_217 = pd.read_csv(data_dir / "fig8_cd8ctl_isg_longitudinal_draws.csv")
    BPZ_COMP = pd.read_csv(data_dir / "suppS5_pbisc_vs_bayesprism_Z_expression.csv")

    P43_WIN = pd.read_csv(data_dir / "suppS5_paired43_window_AUC_true11.csv")
    GEE_SRC = pd.read_csv(data_dir / "suppS5_source_GEE_robust_CI_true11.csv")
    LAB_IND = pd.read_csv(data_dir / "fig8_hospital_labs_true11.csv")
    CV_500K = pd.read_csv(data_dir / "fig8_clinical_CV_infold_INT_500k_true11.csv")

    # =====================================================================
    # 1) MAIN FIGURE 7 (a - q) — VERIFIED TRUE_PUB_11 MORTALITY
    # =====================================================================
    def draw_7a(ax):
        ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
        arr_kw = dict(arrowstyle="-|>", lw=1.9, color=PALETTE["gray_2"], mutation_scale=14)

        rbox(ax, 0.4, 8, 20.6, 80, COL["CardGray"], COL["Real"], lw=1.7)
        ax.text(10.7, 78.5, "Imperial ESKD Cohort\n(107 Donors · 105 Bulk)", ha="center", va="center", fontsize=9.4, fontweight="bold", color=COL["Real"])
        ax.text(1.7, 53.5, "• Hemodialysis (med. 71 yr)\n• 68 COVID-19+ · 39 COVID−", ha="left", va="center", fontsize=8.2, color=PALETTE["gray_2"], linespacing=1.32)
        rbox(ax, 1.8, 13.5, 17.8, 21.0, COL["CardRed"], COL["CLL"], lw=1.3, r=1.5)
        ax.text(10.7, 24.0, "ROR1 / FMOD Screen\n−2 CLL (C116, C139)", ha="center", va="center", fontsize=8.0, fontweight="bold", color=COL["CLL"])
        ax.annotate("", xy=(24.0, 48), xytext=(21.4, 48), arrowprops=arr_kw)

        rbox(ax, 24.4, 5, 49.4, 85, "white", PALETTE["gray_1"], lw=1.4, ls="--")
        ax.text(49.1, 84.0, "66 Unique Non-CLL COVID-19+ Patients (217 longitudinal draws · 11 verified deaths*)", ha="center", va="center", fontsize=9.2, fontweight="bold", color=COL["Ext"])

        rbox(ax, 25.6, 24, 12.6, 50, "#EAEAEA", COL["Real"], lw=1.6)
        ax.text(31.9, 64.5, "Real-only\nscRNA-seq", ha="center", va="center", fontsize=8.8, fontweight="bold", color=COL["Real"])
        ax.text(31.9, 43.0, "10 patients\n\nCITE-seq\ndroplet only\n(Camb. 6, NCL 4)", ha="center", va="center", fontsize=7.9, color=PALETTE["gray_2"], linespacing=1.24)

        rbox(ax, 39.4, 24, 16.2, 50, COL["CardTeal"], COL["Shared"], lw=1.8)
        ax.text(47.5, 64.5, "Shared Bridge\n(Real ∩ PBISC)", ha="center", va="center", fontsize=8.8, fontweight="bold", color="#1F636B")
        ax.text(40.5, 42.5, "• 27 patients\n• Non-contemp.\n  fidelity/scale\n• Counted once", ha="left", va="center", fontsize=7.8, color=COL["Real"], linespacing=1.25)

        rbox(ax, 56.8, 24, 15.8, 50, COL["CardBlue"], COL["New"], lw=1.8)
        ax.text(64.7, 64.5, "PBISC new\n(Bulk → scRNA)", ha="center", va="center", fontsize=8.8, fontweight="bold", color=COL["Ext"])
        ax.text(57.9, 42.5, "• 29 patients\n• 79 Wave-1 draws\n• Pre-RECOVERY era\n• +6 direct deaths", ha="left", va="center", fontsize=7.8, color=COL["Real"], linespacing=1.28)

        ax.plot([25.8, 55.4], [17.5, 17.5], color=COL["Real"], lw=2.2)
        ax.text(40.6, 11.2, "Real: n = 37 (138 draws · 5 verified deaths)", ha="center", va="center", fontsize=8.1, fontweight="bold", color=COL["Real"])
        ax.plot([39.6, 72.4], [20.4, 20.4], color=COL["New"], lw=1.8, ls=":")
        ax.text(64.7, 11.2, "+ PBISC new: +29 (+78%)", ha="center", va="center", fontsize=8.1, fontweight="bold", color=COL["New"])
        ax.annotate("", xy=(77.0, 48), xytext=(74.2, 48), arrowprops=arr_kw)

        rbox(ax, 77.4, 8, 22.0, 80, COL["CardBlue"], COL["Ext"], lw=2.0)
        ax.text(88.4, 77.5, "Harmonized Cohort\nReal + PBISC (n = 66)", ha="center", va="center", fontsize=9.6, fontweight="bold", color=COL["Ext"])
        ax.text(78.6, 43.0, "• Within-source INT\n  standardizes ranks\n• +78% sample size\n• 2.2× verified deaths\n  (5 → 11; *C101/C146/\n   C147 fixed as alive)",
                ha="left", va="center", fontsize=7.9, color=COL["Real"], linespacing=1.30)

    def draw_7b(ax):
        cats_b = [
            ("Mild (WHO 0)", 13, 10), ("Moderate (WHO 1)", 9, 5), ("Severe (WHO 2)", 10, 7), ("Critical (WHO 3)", 5, 7),
            ("Published COVID deaths (11)", 5, 6), ("Direct COVID only (9)", 3, 6), ("Zenodo w1+w2 raw (14)*", 8, 6),
        ]
        yb = np.arange(len(cats_b))[::-1]
        yb_adj = np.where(yb < 3, yb - 0.38, yb)
        ax.plot([0, 32.0], [2.80, 2.80], color=PALETTE["neutral"], ls="--", lw=1.1, zorder=1)
        for yi, (lab, nr, nn) in zip(yb_adj, cats_b):
            ax.barh(yi, nr, height=0.64, color=COL["Real"], edgecolor="black", lw=1.1, label="Real (37)" if yi == yb_adj[0] else "")
            ax.barh(yi, nn, left=nr, height=0.64, color=COL["New"], edgecolor="black", lw=1.1, label="+ PBISC new (29)" if yi == yb_adj[0] else "")
            pct_gain = int(round(100.0 * nn / nr))
            ax.text(nr + nn + 0.6, yi, f"{nr + nn} (+{pct_gain}%)", va="center", ha="left", fontsize=8.3,
                    fontweight="bold", color=PALETTE["red_strong"] if yi == yb_adj[4] else (PALETTE["gray_1"] if yi == yb_adj[6] else COL["Ext"]))
            ax.text(nr / 2.0, yi, f"{nr}", va="center", ha="center", fontsize=8.0, color="white", fontweight="bold")
            ax.text(nr + nn / 2.0, yi, f"{nn}", va="center", ha="center", fontsize=8.0, color="white", fontweight="bold")
        ax.set_yticks(yb_adj); ax.set_yticklabels([c[0] for c in cats_b], fontsize=8.7)
        ax.get_yticklabels()[4].set_fontweight("bold"); ax.get_yticklabels()[4].set_color(PALETTE["red_strong"])
        ax.tick_params(axis="y", length=0, pad=4); ax.spines["left"].set_visible(False)
        ax.set_xlim(0, 34.5); ax.set_ylim(-1.25, 6.55); ax.set_xticks([0, 10, 20, 30])
        ax.set_xlabel("Number of COVID-19+ patients", fontsize=9.8)
        ax.legend(loc="lower right", bbox_to_anchor=(1.02, -0.02), ncol=2, frameon=False, fontsize=7.9, handlelength=1.0, handletextpad=0.3, columnspacing=0.8)

    def draw_7c(ax):
        df_b = CLL_BULK.set_index("patient_id")
        neg = df_b[~df_b["covid_positive_scope"] & (df_b["cohort_group"] != "Excluded CLL (C116, C139)")]
        pos_real = df_b[df_b["covid_positive_scope"] & (df_b["cohort_group"] == "Real / Shared donor")]
        pos_pbisc = df_b[df_b["covid_positive_scope"] & (df_b["cohort_group"] == "PBISC-only donor")]
        cll_pts = df_b[df_b["cohort_group"] == "Excluded CLL (C116, C139)"]
        ax.scatter(neg["ROR1"], neg["FMOD"], s=22, marker="o", c="white", edgecolors=PALETTE["neutral"], linewidths=0.9, alpha=0.8, zorder=1, label=f"COVID− ESKD ({len(neg)})")
        ax.scatter(pos_real["ROR1"], pos_real["FMOD"], s=34, marker="D", c=PALETTE["neutral"], edgecolors=COL["Real"], linewidths=1.0, alpha=0.92, zorder=2, label=f"Real/Shared ({len(pos_real)})")
        ax.scatter(pos_pbisc["ROR1"], pos_pbisc["FMOD"], s=40, marker="^", c=COL["New"], edgecolors="black", linewidths=0.8, alpha=0.92, zorder=3, label=f"PBISC-only ({len(pos_pbisc)})")
        ax.scatter(cll_pts["ROR1"], cll_pts["FMOD"], s=105, marker="D", c=COL["CLL"], edgecolors="black", linewidths=1.2, zorder=5, label="CLL (C116, C139)")
        ax.annotate("C116\n(CLL)", xy=(df_b.loc["C116", "ROR1"], df_b.loc["C116", "FMOD"]), xytext=(3.15, 4.30), fontsize=8.5, fontweight="bold", color=COL["CLL"], arrowprops=dict(arrowstyle="-", color=COL["CLL"], lw=1.2))
        ax.annotate("C139\n(CLL)", xy=(df_b.loc["C139", "ROR1"], df_b.loc["C139", "FMOD"]), xytext=(3.40, 1.95), fontsize=8.5, fontweight="bold", color=COL["CLL"], arrowprops=dict(arrowstyle="-", color=COL["CLL"], lw=1.2))
        ax.axhline(1.0, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.75)
        ax.axvline(2.5, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.75)
        ax.set_xlim(-0.3, 6.2); ax.set_ylim(-0.35, 6.55)
        ax.set_xticks([0, 2, 4, 6]); ax.set_yticks([0, 2, 4, 6])
        ax.set_xlabel("Bulk ROR1 log2(CPM+1)", fontsize=9.8); ax.set_ylabel("Bulk FMOD log2(CPM+1)", fontsize=9.8)
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.00), fontsize=7.5, frameon=False, handletextpad=0.22, labelspacing=0.18, borderpad=0.05)

    def draw_7d(ax):
        steps = [
            ("Tested fine types", len(SEL), "≥30% Real detection"),
            ("FDR-associated\n(Real + PBISC new)", int((SEL["S2_sev"] | SEL["S2_death"]).sum()), "sev q43<0.05 / dth q43<0.10"),
            ("Same direction\nin both sources", int((SEL["S3_sev"] | SEL["S3_death"]).sum()), "− NKT (Real ρ = −0.01)"),
            ("Retained after\nsteroid/dex check", int((SEL["S4_sev"] | SEL["S4_death"]).sum()), "− CD14 mono dex"),
            ("PBISC-enabled\n(not FDR in Real)", int((SEL["S5_sev"] | SEL["S5_death"]).sum()), "− 5 Real-sig. anchors"),
        ]
        y = np.arange(len(steps))[::-1]
        fills = ["#DCE6F2", "#B3CBE6", "#7FA9D6", COL["New"], COL["Ext"]]
        for yi, (lab, cnt, note), fc in zip(y, steps, fills):
            ax.barh(yi, cnt, height=0.66, color=fc, edgecolor=COL["Real"], linewidth=1.3, zorder=2)
            ax.text(cnt + 1.2, yi, f"{cnt}", va="center", ha="left", fontsize=9.8, fontweight="bold", color=COL["Real"])
            ax.text(cnt + (8.2 if cnt >= 10 else 5.8), yi, note, va="center", ha="left", fontsize=7.7,
                    color=PALETTE["red_strong"] if note.startswith("−") else PALETTE["gray_2"])
        ax.set_yticks(y); ax.set_yticklabels([s[0] for s in steps], fontsize=8.8)
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlim(0, 92); ax.set_xticks([0, 20, 40]); ax.set_xlabel("Fine-grained cell types", fontsize=9.8)

    def draw_7e(fig, ax, cax_rect):
        keep_feats = SCAN_NO[(SCAN_NO["q_sev"] < 0.05) | (SCAN_NO["q_death"] < 0.10)]["feature"].unique()
        keep = [f for f in LBL if f in keep_feats and f != "CD16mono|lin"]
        ext = SCAN_NO[SCAN_NO["set"] == "Real + PBISC new"].set_index("feature")
        keep = sorted(keep, key=lambda f: ext.loc[f, "rho"])
        ROWS = [
            ("Real", "rho", "q_sev", 0.05, "Severity · Real (37)", 0.0),
            ("PBISC new", "rho", "q_sev", 0.05, "Severity · PBISC new (29)", -1.0),
            ("Real + PBISC new", "rho", "q_sev", 0.05, "Severity · Combined (66)", -2.0),
            ("Real", "AUC_death", "q_death", 0.10, "Death (5d) · Real", -3.35),
            ("PBISC new", "AUC_death", "q_death", 0.10, "Death (6d) · PBISC new", -4.35),
            ("Real + PBISC new", "AUC_death", "q_death", 0.10, "Death (11d) · Combined", -5.35),
        ]
        cmap = mpl.colormaps["RdBu_r"]; norm = mpl.colors.Normalize(-0.8, 0.8)
        ax.axhline(-2.675, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=1)
        for st, vcol, qcol, q_thr, rlab, y_pos in ROWS:
            for i, f in enumerate(keep):
                r = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == st)]
                if not len(r) or pd.isna(r[vcol].iloc[0]): continue
                r = r.iloc[0]
                val = r[vcol] if vcol == "rho" else 2.0 * (r[vcol] - 0.5)
                is_fdr = r[qcol] < q_thr
                size = 24 + 135 * min(-np.log10(max(r[qcol], 1e-4)), 4.0) / 4.0
                ax.scatter(i, y_pos, s=size, c=[cmap(norm(val))], edgecolors="black" if is_fdr else "none", linewidths=1.4 if is_fdr else 0.0, zorder=3)
        ax.set_yticks([r[5] for r in ROWS]); ax.set_yticklabels([r[4] for r in ROWS], fontsize=8.6)
        ax.set_xticks(range(len(keep))); ax.set_xticklabels([LBL[f] for f in keep], rotation=44, ha="right", fontsize=8.5)
        ax.tick_params(length=0)
        for sp_e in ax.spines.values(): sp_e.set_visible(False)
        ax.set_xlim(-0.6, len(keep) - 0.4); ax.set_ylim(-5.95, 1.05)
        for i, f in enumerate(keep):
            t = ax.get_xticklabels()[i]
            if f in ANCHOR_SET: t.set_fontweight("bold"); t.set_color(COL["Real"])
            elif f in FOCUS_SET: t.set_fontweight("bold"); t.set_color(PALETTE["red_strong"] if f == "HSPC" else COL["Ext"])
            p_het = ext.loc[f, "p_het_sev"]
            flag = "S" if f == "CD14mono_dex" else ("H" if (pd.notna(p_het) and p_het < 0.05) else "")
            if flag:
                ax.text(i, 0.68, flag, ha="center", va="center", fontsize=8.5, color=PALETTE["red_strong"], fontweight="bold")
        cax = fig.add_axes(cax_rect)
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm, cmap), cax=cax, orientation="horizontal")
        cb.set_ticks([-0.8, 0, 0.8]); cb.ax.xaxis.set_ticks_position("top")
        cb.ax.tick_params(labelsize=8.2, length=2, pad=1); cb.outline.set_visible(False)

    def draw_7f(ax):
        ITEMS = [
            ("NK CD56 bright", "sev", "NK CD56bright", "anchor"),
            ("B_ASC_dividing", "sev", "B ASC dividing", "anchor"),
            ("T CD8 CM", "sev", "T CD8 CM (sev.)", "anchor"),
            ("CD16mono", "sev", "CD16 mono", "anchor"),
            ("B_ASC_IgG", "sev", "B ASC IgG (sev.)", "focus"),
            ("B_ASC_IgA", "sev", "B ASC IgA (sev.)", "focus"),
            ("HSPC", "sev", "HSPC (sev.)", "hspc"),
            ("T CD8 CTL", "sev", "T CD8 CTL (sev.)", "focus"),
            ("T CD8 CTL", "death", "T CD8 CTL (death)", "focus"),
            ("HSPC", "death", "HSPC (death)", "hspc"),
        ]
        y_positions, y_cur = [], 0.0
        for idx, (f, ep, disp, grp) in enumerate(ITEMS):
            if idx == 4:
                ax.axhline(y_cur + 0.42, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=0); y_cur -= 0.25
            y_positions.append(y_cur)
            d = PW_UNI[(PW_UNI["feature"] == f) & (PW_UNI["outcome"] == ep)].sort_values("n_real")
            a, _ = n80_interp(d["n_real"], d["real_nom"])
            b, le = n80_interp(d["n_real"], d["plus_new29_nom"])
            ra = a if np.isfinite(a) else 40.0
            ext_col = PALETTE["red_strong"] if grp == "hspc" else (COL["Ext"] if grp == "focus" else COL["Real"])
            ax.plot([b, ra], [y_cur, y_cur], color=PALETTE["neutral"], lw=2.3, zorder=1)
            ax.scatter(ra, y_cur, marker="D", s=40, zorder=3, c=COL["Real"] if np.isfinite(a) else "white", edgecolors=COL["Real"], linewidths=1.5)
            ax.scatter(b, y_cur, marker="s", s=44, zorder=4, c=ext_col if np.isfinite(b) else "white", edgecolors=ext_col)
            ax.text(ra + 1.2, y_cur, f"{a:.0f}" if np.isfinite(a) else "n.r.", va="center", fontsize=8.3, color=COL["Real"])
            ax.text(b - 1.2, y_cur, "≤10" if le else f"{b:.0f}", va="center", ha="right", fontsize=8.3, color=ext_col, fontweight="bold")
            y_cur -= 1.0
        ax.set_yticks(y_positions); ax.set_yticklabels([it[2] for it in ITEMS], fontsize=8.7)
        for idx_t, t in enumerate(ax.get_yticklabels()):
            grp = ITEMS[idx_t][3]
            if grp != "anchor":
                t.set_fontweight("bold"); t.set_color(PALETTE["red_strong"] if grp == "hspc" else COL["Ext"])
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlim(0.5, 46.5); ax.set_xticks([10, 20, 30, 40]); ax.set_xticklabels(["10", "20", "30", "n.r."])
        ax.set_ylim(y_positions[-1] - 0.55, 0.6); ax.set_xlabel("Sampled Real pts for 80% P < 0.05", fontsize=9.5)

    def draw_7g(ax):
        cap = 4.15; nl = lambda q: np.minimum(-np.log10(np.maximum(q, 1e-6)), cap)
        q05 = -np.log10(0.05)
        ax.axhspan(q05, cap, xmin=0, xmax=q05 / cap, color=PALETTE["green_1"], alpha=0.85, zorder=0)
        ax.plot([0, cap], [0, cap], color=PALETTE["neutral"], lw=1.3, zorder=1)
        ax.axhline(q05, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1); ax.axvline(q05, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1)
        ax.axhline(1.0, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.65, zorder=1); ax.axvline(1.0, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.65, zorder=1)
        for out_k, mk, qr, qe in [("sev", "o", "q_sev_Real", "q_sev_Ext"), ("death", "^", "q_death_Real", "q_death_Ext")]:
            for _, r in SEL.iterrows():
                xv, yv = float(nl(r[qr])), float(nl(r[qe]))
                is_s5 = bool(r[f"S5_{out_k}"])
                is_foc = is_s5 and (r["feature"] in FOCUS_SET)
                is_anc = (out_k == "sev") and bool(r["Anchor_sev"])
                fc, ec, sz, zo = (PALETTE["red_strong"] if r["feature"] == "HSPC" else COL["Ext"], "black", 62, 5) if is_foc else ((COL["Real"], "white", 50, 4) if is_anc else ((PALETTE["gray_1"], "none", 32, 3) if is_s5 else (PALETTE["neutral"], "none", 20, 2)))
                ax.scatter(xv, yv, s=sz, marker=mk, c=fc, edgecolors=ec, linewidths=1.0 if ec != "none" else 0.0, zorder=zo)
        callouts = [
            ("B_ASC_IgG", "q_sev_Real", "q_sev_Ext", "B ASC IgG (sev)", (0.10, 3.74), COL["Ext"]),
            ("B_ASC_IgA", "q_sev_Real", "q_sev_Ext", "B ASC IgA (sev)", (0.10, 3.16), COL["Ext"]),
            ("HSPC", "q_sev_Real", "q_sev_Ext", "HSPC (sev)", (0.06, 2.38), PALETTE["red_strong"]),
            ("T CD8 CTL", "q_sev_Real", "q_sev_Ext", "T CD8 CTL (sev)", (0.06, 1.78), COL["Ext"]),
            ("HSPC", "q_death_Real", "q_death_Ext", "HSPC (death)", (1.42, 0.92), PALETTE["red_strong"]),
            ("T CD8 CTL", "q_death_Real", "q_death_Ext", "T CD8 CTL (death)", (1.42, 0.42), COL["Ext"]),
            ("NK CD56 bright", "q_sev_Real", "q_sev_Ext", "NK CD56bright", (2.30, 3.88), COL["Real"]),
            ("T CD8 CM", "q_sev_Real", "q_sev_Ext", "T CD8 CM (sev)", (2.42, 3.05), COL["Real"]),
            ("B_ASC_dividing", "q_sev_Real", "q_sev_Ext", "B ASC dividing", (2.42, 2.48), COL["Real"]),
            ("CD16mono", "q_sev_Real", "q_sev_Ext", "CD16 mono", (2.35, 1.72), COL["Real"]),
        ]
        sel_ix = SEL.set_index("feature")
        for f, qr, qe, txt, xy_t, tcol in callouts:
            rr = sel_ix.loc[f]
            ax.annotate(txt, xy=(float(nl(rr[qr])), float(nl(rr[qe]))), xytext=xy_t, fontsize=7.6, fontweight="bold", color=tcol, va="center", ha="left",
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.84),
                        arrowprops=dict(arrowstyle="-", color=PALETTE["gray_1"], lw=0.9, shrinkA=1, shrinkB=4), zorder=6)
        ax.set_xlim(0, cap); ax.set_ylim(0, cap); ax.set_xticks([0, 1, 2, 3, 4]); ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_xlabel("−log10 q43, Real alone", fontsize=9.8); ax.set_ylabel("−log10 q43, Real + PBISC new", fontsize=9.8)

    def draw_7h(ax):
        cols = ["Real\n(ρ / AUC)", "PBISC new\n(ρ / AUC)", "Diff. P†\n(P_het)", "Shared ρ†\n(n = 27)", "Sensitivity\n(of 14)", "WHO 2–3\ndeath AUC"]
        GOOD, MID, BAD = COL["ScoreGood"], COL["ScoreMid"], COL["ScoreBad"]
        for i, (_, r) in enumerate(SC.iterrows()):
            o = r["outcome"]; fmt = (lambda z: f"{z:+.2f}") if o == "sev" else (lambda z: f"{z:.2f}")
            sh_txt = f"{r['shared_rho']:+.2f}{'*' if r['shared_p'] < 0.05 else ''}"
            rob_txt = f"{int(r['robust_pass'])}/{int(r['robust_total'])}"
            cells = [
                (fmt(r["effect_Real"]), "white", COL["Real"]),
                (fmt(r["effect_New"]), "white", COL["New"]),
                (f"{r['p_het']:.2f}", GOOD if r["p_het"] >= 0.05 else BAD, "white"),
                (sh_txt, GOOD if r["shared_p"] < 0.05 else (MID if r["shared_rho"] >= 0.25 else BAD), "white"),
                (rob_txt, GOOD if r["robust_pass"] >= 12 else (MID if r["robust_pass"] >= 9 else BAD), "white"),
                ((f"{r['AUC_within_WHO23']:.2f}", GOOD if r["p_within_WHO23"] < 0.05 else BAD, "white") if o == "death" else ("–", "white", PALETTE["gray_1"])),
            ]
            if i % 2 == 1:
                ax.add_patch(Rectangle((0.0, i + 0.05), 2.0, 0.90, fc="#F5F6F8", ec="none", zorder=0))
            for j, (txt, fc, tc) in enumerate(cells):
                if fc != "white":
                    ax.add_patch(FancyBboxPatch((j + 0.06, i + 0.10), 0.88, 0.80, boxstyle="round,pad=0.02,rounding_size=0.10", fc=fc, ec="#222222", lw=0.45, zorder=2))
                ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=8.8, color=tc,
                        fontweight="bold" if (fc != "white" or (j in (0, 1) and r["followed"])) else "normal", zorder=3)
        ax.set_xlim(0, len(cols)); ax.set_ylim(len(SC), 0)
        ax.set_yticks(np.arange(len(SC)) + 0.5)
        ax.set_yticklabels([f"{LBL[r.feature]} · {'sev.' if r.outcome == 'sev' else 'death'}" for _, r in SC.iterrows()], fontsize=8.9)
        for idx_t, t in enumerate(ax.get_yticklabels()):
            if SC.iloc[idx_t]["followed"]:
                t.set_fontweight("bold"); t.set_color(PALETTE["red_strong"] if SC.iloc[idx_t]["feature"] == "HSPC" else COL["Ext"])
            else:
                t.set_color(PALETTE["gray_2"])
        ax.set_xticks(np.arange(len(cols)) + 0.5); ax.set_xticklabels(cols, fontsize=8.5, fontweight="bold")
        ax.xaxis.tick_top(); ax.tick_params(length=0)
        for sp_h in ax.spines.values(): sp_h.set_visible(False)
        ax.plot([-0.34, 1.0], [6.0, 6.0], transform=ax.get_yaxis_transform(), color=COL["Real"], lw=1.5, clip_on=False)
        ax.text(-0.33, 3.0, "Prioritized", transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center", fontsize=8.8, fontweight="bold", color=COL["Ext"], clip_on=False)
        ax.text(-0.33, 7.0, "Deprior.", transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center", fontsize=8.3, color=PALETTE["gray_2"], clip_on=False)

    def draw_7i(ax_sev, ax_dth):
        FEATS_E = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
        ci_no = CI[CI["CLL_excluded"] == True]
        y_centers, y_cur = [], 0.0
        for idx, f in enumerate(FEATS_E):
            if idx == 2:
                for axx in (ax_sev, ax_dth): axx.axhline(y_cur + 0.52, color=PALETTE["neutral"], ls="--", lw=1.1, zorder=1)
                y_cur -= 0.35
            y_centers.append(y_cur)
            for sname, k, dy, sz in [("Real", "Real", +0.26, 36), ("PBISC new", "New", 0.0, 38), ("Real + PBISC new", "Ext", -0.26, 42)]:
                row = ci_no[(ci_no["feature"] == f) & (ci_no["set"] == sname)].iloc[0]
                ax_sev.plot([row["rho_lo"], row["rho_hi"]], [y_cur + dy, y_cur + dy], color=COL[k], lw=1.5, alpha=0.8, zorder=2)
                ax_sev.scatter(row["rho"], y_cur + dy, s=sz, marker=MK[k], c="white" if row["p_sev"] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.4, zorder=3)
                ax_dth.plot([row["AUC_lo"], row["AUC_hi"]], [y_cur + dy, y_cur + dy], color=COL[k], lw=1.5, alpha=0.8, zorder=2)
                ax_dth.scatter(row["AUC"], y_cur + dy, s=sz, marker=MK[k], c="white" if row["p_death"] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.4, zorder=3)
            y_cur -= 1.05
        ax_sev.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3); ax_dth.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.3)
        ax_sev.set_yticks(y_centers); ax_sev.set_yticklabels([LBL[f] for f in FEATS_E], fontsize=9.1)
        for idx_t, t in enumerate(ax_sev.get_yticklabels()):
            f = FEATS_E[idx_t]; t.set_fontweight("bold")
            t.set_color(COL["Real"] if idx_t < 2 else (PALETTE["red_strong"] if f == "HSPC" else COL["Ext"]))
        ax_dth.set_yticks(y_centers); ax_dth.set_yticklabels([])
        for axx in (ax_sev, ax_dth):
            axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(y_centers[-1] - 0.55, 0.55)
        ax_sev.set_xlim(-0.92, 0.92); ax_sev.set_xticks([-0.5, 0.0, 0.5]); ax_sev.set_xlabel("Spearman ρ (severity)", fontsize=9.8)
        ax_dth.set_xlim(0.0, 1.05); ax_dth.set_xticks([0.25, 0.50, 0.75, 1.0]); ax_dth.set_xticklabels(["0.25", "0.5", "0.75", "1"]); ax_dth.set_xlabel("AUC (verified death)", fontsize=9.8)

    def draw_strip_ns(ax, f, outcome, rng_local, show_ylabel=True, compact_txt=False):
        if outcome == "sev":
            for g in range(4):
                for k, dx in [("Real", -0.18), ("New", 0.18)]:
                    vals = PT_NS[(PT_NS["src"] == k) & (PT_NS["sev"] == g)][f].dropna().values
                    ax.scatter(g + dx + rng_local.uniform(-0.075, 0.075, len(vals)), vals, s=24, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                    if len(vals): ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.9, zorder=4)
            ax.set_xticks(range(4)); ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=9.2); ax.set_xlim(-0.55, 3.55)
        else:
            for g in range(2):
                for k, dx in [("Real", -0.18), ("New", 0.18)]:
                    vals = PT_NS[(PT_NS["src"] == k) & (PT_NS["fat11"] == g)][f].dropna().values
                    ax.scatter(g + dx + rng_local.uniform(-0.075, 0.075, len(vals)), vals, s=26, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                    if len(vals): ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.9, zorder=4)
            ax.set_xticks([0, 1]); ax.set_xticklabels([f"Survived\n({int((PT_NS['fat11']==0).sum())})", f"Died\n({int((PT_NS['fat11']==1).sum())})"], fontsize=9.1); ax.set_xlim(-0.52, 1.52)
        ax.set_ylim(-2.6, 5.15); ax.set_yticks([-2, 0, 2])
        if show_ylabel: ax.set_ylabel("Normal score (within-source)", fontsize=9.5)
        p_fn = p_eq_s if compact_txt else fp
        fs_t = 7.05 if compact_txt else 7.65
        for idx_s, (sname, k, disp_s) in enumerate([("Real", "Real", "Real"), ("PBISC new", "New", "PBISC new"), ("Real + PBISC new", "Ext", "Comb." if compact_txt else "Combined")]):
            r = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == sname)].iloc[0]
            if outcome == "sev":
                q_s = f" (q43={r.q_sev:.3f})" if (sname == "Real + PBISC new" and r.q_sev >= 0.001) else (f" (q43={r.q_sev:.1e})" if sname == "Real + PBISC new" else "")
                txt, sig = f"{disp_s}: ρ={r.rho:+.2f}, {p_fn(r.p_sev)}{q_s}".replace("e-0", "e-"), r.p_sev < 0.05
            else:
                q_s = f" (q43={r.q_death:.3f})" if sname == "Real + PBISC new" else ""
                txt, sig = f"{disp_s}: AUC={r.AUC_death:.2f}, {p_fn(r.p_death)}{q_s}", r.p_death < 0.05
            ax.text(0.03, 0.975 - 0.092 * idx_s, txt, transform=ax.transAxes, fontsize=fs_t, color=COL[k], va="top", fontweight="bold" if sig else "normal")

    def draw_7p(ax):
        ORDER_H = [
            ("Published COVID deaths (11)", "Published 11 (L1O: .74–.81)"),
            ("Direct COVID only (9)", "Direct COVID only (9)"),
            ("Zenodo w1+w2 raw (14)", "Zenodo w1+w2 raw (14)"),
            ("Within WHO 2-3", "Within WHO 2–3"),
            ("Source_x_PeakWHO", "Source × Peak WHO (exact)"),
            ("Cambridge + PBISC new", "Cambridge + PBISC new"),
            ("Samples > day 7", "Samples > day 7"),
            ("PBISC new alone", "PBISC new alone (6d)"),
        ]
        yy = np.arange(len(ORDER_H))[::-1]
        ax.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1)
        l1o_min, l1o_max = float(HS_L1O["AUC"].min()), float(HS_L1O["AUC"].max())
        cond_r = HS_COND.loc["HSPC"]
        for yi, (an_key, _) in zip(yy, ORDER_H):
            if an_key == "PBISC new alone":
                rn = HS[HS["analysis"] == an_key].iloc[0]
                ax.scatter(rn["AUC"], yi, s=44, marker="o", c=COL["New"], edgecolors=COL["New"], linewidths=1.4, zorder=4)
                ax.text(1.02, yi, f"{int(rn['deaths'])} · {fp_short(rn['p'])}", va="center", ha="left", fontsize=7.6, color=COL["New"], fontweight="bold")
            elif an_key == "Source_x_PeakWHO":
                ax.scatter(float(cond_r["cond_AUC"]), yi, s=46, marker="s", c="white", edgecolors=COL["Ext"], linewidths=1.5, zorder=4)
                ax.text(1.02, yi, f"11 · {float(cond_r['exact_p_2s']):.2f} (ns)", va="center", ha="left", fontsize=7.5, color=PALETTE["gray_2"])
            else:
                rr = HS[(HS["analysis"] == an_key) & (HS["set"] == "Real")].iloc[0]
                re = HS[(HS["analysis"] == an_key) & (HS["set"] == "Real + PBISC new")].iloc[0]
                if an_key == "Published COVID deaths (11)":
                    ax.plot([l1o_min, l1o_max], [yi - 0.16, yi - 0.16], color=PALETTE["red_strong"], lw=2.4, alpha=0.55, zorder=2)
                ax.plot([rr["AUC"], re["AUC"]], [yi + 0.16, yi - 0.16], color=PALETTE["neutral"], lw=1.2, zorder=2)
                ax.scatter(rr["AUC"], yi + 0.16, s=36, marker="D", c="white" if rr["p"] >= 0.05 else COL["Real"], edgecolors=COL["Real"], linewidths=1.3, zorder=3)
                ax.scatter(re["AUC"], yi - 0.16, s=40, marker="s", c="white" if re["p"] >= 0.05 else COL["Ext"], edgecolors=COL["Ext"], linewidths=1.3, zorder=4)
                ax.text(1.02, yi, f"{int(rr['deaths'])}→{int(re['deaths'])} · {fp_short(re['p'])}", va="center", ha="left", fontsize=7.6, color=COL["Ext"], fontweight="bold" if re["p"] < 0.05 else "normal")
        ax.text(1.02, yy[0] + 0.72, "deaths · +PBISC P", va="bottom", ha="left", fontsize=7.5, color=PALETTE["gray_2"], fontweight="bold")
        ax.set_yticks(yy); ax.set_yticklabels([lab for _, lab in ORDER_H], fontsize=8.0)
        ax.get_yticklabels()[4].set_color(PALETTE["red_strong"]); ax.get_yticklabels()[4].set_fontweight("bold")
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlim(0.46, 1.00); ax.set_ylim(-0.6, len(ORDER_H) - 0.25)
        ax.set_xticks([0.5, 0.75, 1.0]); ax.set_xticklabels(["0.5", "0.75", "1.0"])
        ax.set_xlabel("HSPC death AUC", fontsize=9.6)

    def draw_7q(ax):
        PP = [
            ("B_ASC_IgG", "sev", "B ASC IgG · sev."),
            ("B_ASC_IgA", "sev", "B ASC IgA · sev."),
            ("HSPC", "sev", "HSPC · sev."),
            ("T CD8 CTL", "sev", "T CD8 CTL · sev."),
            ("T CD8 CTL", "death", "T CD8 CTL · death"),
            ("HSPC", "death", "HSPC · death"),
        ]
        yy = np.arange(len(PP))[::-1]; hb = 0.25; pw20 = PW_UNI[PW_UNI["n_real"] == 20]
        for nom_col, bh_col, fc_bh, fc_nom, dy, hatch in [
            ("real_nom", "real_bh43", PALETTE["gray_2"], PALETTE["neutral"], +hb, ""),
            ("plus_new29_nom", "plus_new29_bh43", COL["Ext"], "#B3CBE6", 0.0, ""),
            ("plus_shuf_nom", "plus_shuf_bh43", PALETTE["gray_1"], "white", -hb, "///"),
        ]:
            v_nom = [float(pw20[(pw20["feature"] == f) & (pw20["outcome"] == o)][nom_col].iloc[0]) for f, o, _ in PP]
            v_bh = [float(pw20[(pw20["feature"] == f) & (pw20["outcome"] == o)][bh_col].iloc[0]) for f, o, _ in PP]
            ax.barh(yy + dy, v_nom, height=hb, color=fc_nom, edgecolor=COL["Real"], linewidth=0.85, hatch=hatch, zorder=2)
            ax.barh(yy + dy, v_bh, height=hb, color=fc_bh, edgecolor=COL["Real"], linewidth=0.95, zorder=3)
        for yi, (f, o, _) in zip(yy, PP):
            rr = pw20[(pw20["feature"] == f) & (pw20["outcome"] == o)].iloc[0]
            bh_pct, nm_pct = int(round(100.0 * rr["plus_new29_bh43"])), int(round(100.0 * rr["plus_new29_nom"]))
            ax.text(1.03, yi, f"{bh_pct}%/{nm_pct}%", va="center", ha="left", fontsize=7.6, color=PALETTE["red_strong"] if f == "HSPC" else COL["Ext"], fontweight="bold")
        ax.text(1.03, yy[0] + 0.65, "BH43 / nom", va="bottom", ha="left", fontsize=7.4, color=PALETTE["gray_2"], fontweight="bold")
        ax.axvline(0.80, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=4)
        ax.set_yticks(yy); ax.set_yticklabels([lab for _, _, lab in PP], fontsize=8.4)
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlim(0.0, 1.01); ax.set_ylim(-0.55, len(PP) - 0.2)
        ax.set_xticks([0.0, 0.5, 0.8, 1.0]); ax.set_xticklabels(["0", "0.5", "0.8", "1"])
        ax.set_xlabel("Re-detection at 20 Real pts", fontsize=9.4)

    fig_m = plt.figure(figsize=(17.8, 26.0), facecolor="white")
    fig_m.text(0.012, 0.992, "Figure 7. Multi-modal cohort assembly, verified published mortality audit, and PBISC-enabled immune associations in ESKD COVID-19",
               fontsize=14.0, fontweight="bold", va="top", ha="left", color=COL["Real"])
    ax_a = fig_m.add_axes([0.016, 0.838, 0.620, 0.124]); draw_7a(ax_a); ax_a.set_title("Multi-modal ESKD COVID-19 cohort assembly via PBISC", x=0.042, ha="left", fontsize=11.2, fontweight="bold", pad=4); fig_m.text(0.012, 0.974, "a", fontsize=17.5, fontweight="bold", va="top")
    ax_b = fig_m.add_axes([0.775, 0.842, 0.198, 0.116]); draw_7b(ax_b); ax_b.set_title("Severity & verified fatal expansion", loc="left", fontsize=11.0, fontweight="bold", pad=6); fig_m.text(0.652, 0.974, "b", fontsize=17.5, fontweight="bold", va="top")
    ax_c = fig_m.add_axes([0.048, 0.668, 0.156, 0.112]); draw_7c(ax_c); ax_c.set_title("Bulk CLL molecular screen", loc="left", fontsize=11.0, fontweight="bold", pad=6); fig_m.text(0.012, 0.800, "c", fontsize=17.5, fontweight="bold", va="top")
    ax_d = fig_m.add_axes([0.318, 0.668, 0.158, 0.112]); draw_7d(ax_d); ax_d.set_title("Structured 5-step screen (43 types)", loc="left", fontsize=10.9, fontweight="bold", pad=6); fig_m.text(0.224, 0.800, "d", fontsize=17.5, fontweight="bold", va="top")
    ax_e = fig_m.add_axes([0.612, 0.668, 0.372, 0.104]); draw_7e(fig_m, ax_e, cax_rect=[0.612, 0.788, 0.082, 0.007])
    fig_m.text(0.700, 0.791, "ρ or 2·(AUC−0.5) · ring: fixed BH43 FDR · H: P_het<0.05 · S: dex proxy", fontsize=7.8, va="center", color=COL["Real"])
    fig_m.text(0.495, 0.800, "e", fontsize=17.5, fontweight="bold", va="top")

    ax_leg = fig_m.add_axes([0.012, 0.584, 0.976, 0.026]); ax_leg.set_axis_off()
    hd_main = [
        Line2D([0], [0], marker="D", ls="", ms=7.2, mfc=COL["Real"], mec=COL["Real"], label="Real (37 patients, 5 verified deaths)"),
        Line2D([0], [0], marker="o", ls="", ms=7.2, mfc=COL["New"], mec=COL["New"], label="PBISC new (29 patients, 6 deaths)"),
        Line2D([0], [0], marker="s", ls="", ms=7.2, mfc=COL["Ext"], mec=COL["Ext"], label="Real + PBISC new (66 patients, 11 deaths)"),
        Line2D([0], [0], marker="o", ls="", ms=7.2, mfc="white", mec=PALETTE["gray_2"], mew=1.5, label="open: P ≥ 0.05 (f, i, p)"),
        Line2D([0], [0], marker="o", ls="", ms=6.2, mfc=COL["Ext"], mec="black", label="severity (g)"),
        Line2D([0], [0], marker="^", ls="", ms=6.8, mfc=COL["Ext"], mec="black", label="death (g)"),
        Patch(facecolor=COL["Ext"], edgecolor=COL["Real"], lw=0.8, label="inner solid: BH43 / outer: P<0.05 (q)"),
        Patch(facecolor="white", edgecolor=COL["Real"], hatch="///", label="outcome-shuffled New 29 (q)"),
        Patch(facecolor=COL["ScoreGood"], edgecolor="#222222", lw=0.5, label="met (h)"),
        Patch(facecolor=COL["ScoreMid"], edgecolor="#222222", lw=0.5, label="borderline (h)"),
        Patch(facecolor=COL["ScoreBad"], edgecolor="#222222", lw=0.5, label="not met (h; †non-contemp. visits)"),
    ]
    ax_leg.legend(handles=hd_main, loc="center", ncol=6, fontsize=8.4, frameon=False, columnspacing=0.90, handletextpad=0.30, labelspacing=0.32)

    ax_f = fig_m.add_axes([0.114, 0.418, 0.156, 0.122]); draw_7f(ax_f); ax_f.set_title("Conditional re-detection n80", loc="left", fontsize=10.9, fontweight="bold", pad=6); fig_m.text(0.012, 0.562, "f", fontsize=17.5, fontweight="bold", va="top")
    ax_g = fig_m.add_axes([0.330, 0.418, 0.180, 0.122]); draw_7g(ax_g); ax_g.set_title("Associations gained by adding PBISC", loc="left", fontsize=10.9, fontweight="bold", pad=6); fig_m.text(0.288, 0.562, "g", fontsize=17.5, fontweight="bold", va="top")
    ax_h = fig_m.add_axes([0.642, 0.414, 0.342, 0.118]); draw_7h(ax_h); ax_h.set_title("Prioritization of PBISC-enabled candidates", loc="left", fontsize=11.0, fontweight="bold", pad=28); fig_m.text(0.532, 0.562, "h", fontsize=17.5, fontweight="bold", va="top")

    ax_i_sev = fig_m.add_axes([0.095, 0.228, 0.120, 0.122]); ax_i_dth = fig_m.add_axes([0.232, 0.228, 0.120, 0.122])
    draw_7i(ax_i_sev, ax_i_dth); ax_i_sev.set_title("Cohort effect sizes (strat. 95% CI)", loc="left", fontsize=10.8, fontweight="bold", pad=8); fig_m.text(0.012, 0.370, "i", fontsize=17.5, fontweight="bold", va="top")

    rng_main = np.random.default_rng(7)
    for fk, let, show_yl, x_l in [("B_ASC_IgG", "j", True, 0.432), ("B_ASC_IgA", "k", False, 0.628), ("HSPC", "l", False, 0.824)]:
        ax_s = fig_m.add_axes([x_l, 0.228, 0.152, 0.122]); draw_strip_ns(ax_s, fk, "sev", rng_main, show_ylabel=show_yl)
        ax_s.set_title(LBL[fk], loc="left", fontsize=11.0, fontweight="bold", color=COL["HSPC"] if fk == "HSPC" else COL["Real"], pad=5)
        fig_m.text(x_l - (0.044 if show_yl else 0.028), 0.370, let, fontsize=17.5, fontweight="bold", va="top")

    for fk, out_k, let, show_yl, x_l in [("T CD8 CTL", "sev", "m", True, 0.048), ("T CD8 CTL", "death", "n", False, 0.200), ("HSPC", "death", "o", False, 0.352)]:
        ax_s = fig_m.add_axes([x_l, 0.040, 0.122, 0.122]); draw_strip_ns(ax_s, fk, out_k, rng_main, show_ylabel=show_yl, compact_txt=True)
        ax_s.set_title(f"{LBL[fk]} ({'sev.' if out_k=='sev' else 'death'})", loc="left", fontsize=10.6, fontweight="bold", color=COL["HSPC"] if fk == "HSPC" else COL["Real"], pad=5)
        fig_m.text(x_l - (0.036 if show_yl else 0.024), 0.182, let, fontsize=17.5, fontweight="bold", va="top")

    ax_p = fig_m.add_axes([0.618, 0.040, 0.086, 0.120]); draw_7p(ax_p); ax_p.set_title("HSPC death sensitivity", loc="left", fontsize=10.8, fontweight="bold", color=COL["Real"], pad=14); fig_m.text(0.492, 0.182, "p", fontsize=17.5, fontweight="bold", va="top")
    ax_q = fig_m.add_axes([0.868, 0.040, 0.072, 0.120]); draw_7q(ax_q); ax_q.set_title("Conditional re-detection", loc="left", fontsize=10.8, fontweight="bold", color=COL["Real"], pad=14); fig_m.text(0.770, 0.182, "q", fontsize=17.5, fontweight="bold", va="top")

    saved_7 = finalize_figure(fig_m, out_dir / "PBISC_Figure7_Main_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Main Figure 7:", saved_7[0].name)

    # =====================================================================
    # 2) MAIN FIGURE 8 (a - i) — 10-METHOD DECONV + HSPC TRANSCRIPTOME + T-CELL ISG + CV
    # =====================================================================
    fig8 = plt.figure(figsize=(18.2, 24.4), facecolor="white")
    fig8.text(0.012, 0.993, "Figure 8. Expanded 10-method bulk deconvolution benchmark, HSPC and T-cell ISG transcriptomes, longitudinal collapse, and mortality CV",
              fontsize=13.6, fontweight="bold", va="top", ha="left", color=COL["Real"])

    # Panel 8a: 10-method deconvolution benchmark (PBISC, Calib-NNLS, BayesPrism, CIBERSORTx S/B, MuSiC, DWLS, Bulk-8, NNLS)
    ax_a8 = fig8.add_axes([0.132, 0.785, 0.355, 0.148]); ax_a_bar = fig8.add_axes([0.495, 0.785, 0.064, 0.148])
    fig8.text(0.012, 0.968, "a", fontsize=17.5, fontweight="bold", va="top")

    METH10 = [
        ("PBISC", "PBISC (generative)"),
        ("Calib-NNLS", "Calib-NNLS† (paired)"),
        ("BayesPrism-init", "BayesPrism (init φ)"),
        ("BayesPrism-upd", "BayesPrism (upd ψ)"),
        ("Bulk-8-markers", "Bulk 8-markers* (BH8)"),
        ("MuSiC", "MuSiC (W-NNLS)"),
        ("CIBERSORTx-Bmode", "CIBERSORTx B-mode"),
        ("CIBERSORTx-Smode", "CIBERSORTx S-mode"),
        ("DWLS", "DWLS"),
        ("NNLS", "Uncalib. NNLS"),
    ]
    COL10_A = [
        ("NK CD56 bright", "sev", "NK56br"), ("B_ASC_dividing", "sev", "B ASC div"),
        ("T CD8 CM", "sev", "CD8 CM"), ("CD16mono", "sev", "CD16m"),
        ("B_ASC_IgG", "sev", "ASC IgG"), ("B_ASC_IgA", "sev", "ASC IgA"),
        ("HSPC", "sev", "HSPC(s)"), ("T CD8 CTL", "sev", "CTL(s)"),
        ("T CD8 CTL", "death", "CTL(d)"), ("HSPC", "death", "HSPC(d)"),
    ]
    pass_cnts = []
    for im, (mk_nm, _) in enumerate(METH10):
        np_m = 0
        for jc, (fk, ep, _) in enumerate(COL10_A):
            row = BENCH10[(BENCH10["method"] == mk_nm) & (BENCH10["feature"] == fk) & (BENCH10["endpoint"] == ep)].iloc[0]
            ok = bool(row["passed_concordant_FDR"]); np_m += int(ok)
            fc = COL["Good"] if ok else COL["Bad"]
            txt = f"{-np.log10(max(row['comb66_q'], 1e-5)):.1f}" if ok else "ns"
            ax_a8.add_patch(FancyBboxPatch((jc + 0.06, im + 0.09), 0.88, 0.82, boxstyle="round,pad=0.015,rounding_size=0.09",
                                           fc=fc, ec="#222222", lw=0.4, hatch="" if ok else "//", alpha=0.95 if ok else 0.78))
            ax_a8.text(jc + 0.5, im + 0.5, txt, ha="center", va="center", fontsize=7.0, color="white", fontweight="bold")
        pass_cnts.append(np_m)
    ax_a8.axvline(4.0, color="black", lw=1.7); ax_a8.set_xlim(0, len(COL10_A)); ax_a8.set_ylim(len(METH10), 0)
    ax_a8.set_yticks(np.arange(len(METH10)) + 0.5); ax_a8.set_yticklabels([m[1] for m in METH10], fontsize=8.3)
    ax_a8.get_yticklabels()[0].set_fontweight("bold"); ax_a8.get_yticklabels()[0].set_color(COL["Ext"])
    ax_a8.set_xticks(np.arange(len(COL10_A)) + 0.5); ax_a8.set_xticklabels([c[2] for c in COL10_A], rotation=36, ha="left", fontsize=7.8, fontweight="bold")
    for ic, t in enumerate(ax_a8.get_xticklabels()):
        t.set_color(COL["Real"] if ic < 4 else (COL["HSPC"] if ic in (6, 9) else COL["Ext"]))
    ax_a8.xaxis.tick_top(); ax_a8.tick_params(length=0)
    for sp_a in ax_a8.spines.values(): sp_a.set_visible(False)
    ax_a8.set_title("10-method bulk deconvolution screen (−log10 q; TRUE 11 deaths)", loc="left", fontsize=10.2, fontweight="bold", pad=22)

    yy_ab = np.arange(len(METH10)) + 0.5
    ax_a_bar.barh(yy_ab, pass_cnts, height=0.68, color=[COL["Ext"] if i == 0 else PALETTE["neutral"] for i in range(len(METH10))], edgecolor=COL["Real"], lw=1.0)
    for yi, pc in zip(yy_ab, pass_cnts):
        ax_a_bar.text(pc + 0.35, yi, f"{pc}/10", va="center", ha="left", fontsize=7.9, fontweight="bold", color=COL["Ext"] if pc == 10 else COL["Real"])
    ax_a_bar.set_ylim(len(METH10), 0); ax_a_bar.set_xlim(0, 14.8); ax_a_bar.set_yticks([]); ax_a_bar.set_xticks([0, 5, 10])
    ax_a_bar.tick_params(axis="y", length=0); ax_a_bar.spines["left"].set_visible(False); ax_a_bar.set_xlabel("Passed", fontsize=8.6)

    # Panel 8b: Wave-1 New-29 HSPC death AUC (95% CI) + HSPC/B ASC IgG/IgA severity rho across all 10 methods
    ax_b1 = fig8.add_axes([0.692, 0.785, 0.122, 0.148]); ax_b2 = fig8.add_axes([0.850, 0.785, 0.128, 0.148])
    fig8.text(0.576, 0.968, "b", fontsize=17.5, fontweight="bold", va="top")
    ax_b1.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3); ax_b2.axvline(0.00, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for im, (mk_nm, _) in enumerate(METH10):
        yi = len(METH10) - 1 - im
        ci_h = NW29_CI.loc[mk_nm]
        rg = BENCH10[(BENCH10["method"] == mk_nm) & (BENCH10["feature"] == "B_ASC_IgG") & (BENCH10["endpoint"] == "sev")].iloc[0]
        ra = BENCH10[(BENCH10["method"] == mk_nm) & (BENCH10["feature"] == "B_ASC_IgA") & (BENCH10["endpoint"] == "sev")].iloc[0]
        sig_h = ci_h["hspc_dth_p"] < 0.05
        c_h = COL["HSPC"] if im == 0 else (PALETTE["gray_3"] if sig_h else PALETTE["gray_1"])
        ax_b1.plot([ci_h["hspc_dth_lo"], ci_h["hspc_dth_hi"]], [yi, yi], color=c_h, lw=1.7 if sig_h else 1.1, alpha=0.85, zorder=2)
        ax_b1.scatter(ci_h["hspc_dth_auc"], yi, s=44 if im == 0 else 32, marker="s" if im == 0 else "o", c=c_h if sig_h else "white", edgecolors=c_h, linewidths=1.35, zorder=4)
        ax_b1.text(1.03, yi, f"{ci_h['hspc_dth_auc']:.2f}{'*' if sig_h else ''}", va="center", ha="left", fontsize=7.4, fontweight="bold" if sig_h else "normal", color=c_h)
        # Right: HSPC sev (red star), B ASC IgG (square), B ASC IgA (triangle)
        c_b = COL["Ext"] if im == 0 else PALETTE["gray_2"]
        ax_b2.scatter(ci_h["hspc_sev_rho"], yi + 0.20, s=46, marker="*", c=COL["HSPC"] if ci_h["hspc_sev_p"] < 0.05 else "white", edgecolors=COL["HSPC"] if im == 0 else PALETTE["gray_1"], linewidths=1.1, zorder=5)
        ax_b2.scatter(rg["new29_eff"], yi, s=30, marker="s", c=c_b if rg["new29_p"] < 0.05 else "white", edgecolors=c_b, linewidths=1.2, zorder=4)
        ax_b2.scatter(ra["new29_eff"], yi - 0.20, s=32, marker="^", c=COL["New"] if ra["new29_p"] < 0.05 else "white", edgecolors=COL["New"] if im == 0 else c_b, linewidths=1.2, zorder=4)
    ax_b1.set_yticks(np.arange(len(METH10))[::-1]); ax_b1.set_yticklabels([m[1] for m in METH10], fontsize=8.0)
    ax_b1.get_yticklabels()[0].set_fontweight("bold"); ax_b1.get_yticklabels()[0].set_color(COL["Ext"])
    ax_b2.set_yticks(np.arange(len(METH10))[::-1]); ax_b2.set_yticklabels([])
    for axx in (ax_b1, ax_b2):
        axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(-0.55, len(METH10) - 0.35)
    ax_b1.set_xlim(0.02, 1.02); ax_b1.set_xticks([0.25, 0.5, 0.8]); ax_b1.set_xlabel("HSPC death AUC (95% CI)", fontsize=8.8)
    ax_b1.set_title("Wave-1 new (n = 29, 6d): rare-cell recovery across 10 methods", loc="left", fontsize=10.1, fontweight="bold", pad=12)
    ax_b2.set_xlim(-0.52, 0.76); ax_b2.set_xticks([-0.4, 0.0, 0.4]); ax_b2.set_xlabel("Sev. ρ: HSPC (*) · IgG (sq) · IgA (tri)", fontsize=8.3)

    ax_leg8 = fig8.add_axes([0.014, 0.736, 0.972, 0.022]); ax_leg8.set_axis_off()
    ax_leg8.legend(handles=[
        Line2D([0], [0], marker="D", ls="", ms=6.8, mfc=COL["Real"], mec=COL["Real"], label="Real (37; HSPC: 33 [5d])"),
        Line2D([0], [0], marker="o", ls="", ms=6.8, mfc=COL["New"], mec=COL["New"], label="PBISC new (29 Wave-1 [6d])"),
        Line2D([0], [0], marker="s", ls="", ms=6.8, mfc=COL["Ext"], mec=COL["Ext"], label="Combined (66 [11d]; HSPC: 62)"),
        Line2D([0], [0], marker="o", ls="", ms=6.8, mfc="white", mec=PALETTE["gray_2"], mew=1.4, label="open: P ≥ 0.05 / CI incl. 0 (i)"),
        Patch(facecolor=COL["Good"], edgecolor="#222222", lw=0.4, label="concordant & FDR (a; †paired cal.; *BH8)"),
        Patch(facecolor=COL["Bad"], edgecolor="#222222", lw=0.4, hatch="//", label="not passed / discordant (a)"),
    ], loc="center", ncol=6, fontsize=8.2, frameon=False, columnspacing=0.85, handletextpad=0.28)

    # Row 2: Point 5 (8c HSPC 16-module severity & death, 8d HSPC Ery5 patient strips) & Point 6 (8e T-cell ISG 13-compartments!)
    ax_c1 = fig8.add_axes([0.148, 0.530, 0.118, 0.170]); ax_c2 = fig8.add_axes([0.282, 0.530, 0.118, 0.170])
    fig8.text(0.012, 0.725, "c", fontsize=17.5, fontweight="bold", va="top")
    HMOD_SHOW = [
        ("Ery5", "Erythroid (Ery5)", COL["HSPC"]),
        ("Cycle4", "Proliferation (Cycle4)", COL["HSPC"]),
        ("MEP5", "MEP-assoc. (MEP5)", COL["Ext"]),
        ("Alarmin_Inflam", "Alarmin / Inflam.", COL["Ext"]),
        ("ISG_HSPC", "HSPC ISG", COL["Teal"]),
        ("HSC5", "Quiescent HSC (HSC5)", PALETTE["gray_2"]),
        ("Myeloid4", "Granulocyte (Myel4)", PALETTE["gray_2"]),
    ]
    yy_c = np.arange(len(HMOD_SHOW))[::-1]
    ax_c1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.2); ax_c2.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.2)
    for yi, (mk_c, _, mcol) in zip(yy_c, HMOD_SHOW):
        sub_m = HSPC_MOD16[HSPC_MOD16["module"] == mk_c].set_index("cohort")
        for snm, k_s, dy, sz in [("Real", "Real", +0.24, 28), ("PBISC new", "New", 0.0, 30), ("Combined", "Ext", -0.24, 38)]:
            rv = sub_m.loc[snm]
            uc = mcol if (k_s == "Ext" and mk_c in ("Ery5", "Cycle4")) else COL[k_s]
            ax_c1.plot([rv["rho_lo"], rv["rho_hi"]], [yi + dy, yi + dy], color=uc, lw=1.35, alpha=0.8)
            ax_c1.scatter(rv["rho_sev"], yi + dy, s=sz, marker=MK[k_s], c=uc if rv["p_sev"] < 0.05 else "white", edgecolors=uc, linewidths=1.3, zorder=3)
            ax_c2.plot([rv["AUC_lo"], rv["AUC_hi"]], [yi + dy, yi + dy], color=uc, lw=1.35, alpha=0.8)
            ax_c2.scatter(rv["AUC_dth11"], yi + dy, s=sz, marker=MK[k_s], c=uc if rv["p_dth11"] < 0.05 else "white", edgecolors=uc, linewidths=1.3, zorder=3)
            if snm == "Combined":
                q_s = "*" if rv["q16_sev"] < 0.05 else ""
                ax_c2.text(1.02, yi, f"ρ={rv['rho_sev']:+.2f}{q_s} | AUC={rv['AUC_dth11']:.2f}", va="center", ha="left", fontsize=7.0,
                           fontweight="bold" if (rv["q16_sev"] < 0.05 or rv["p_dth11"] < 0.05) else "normal",
                           color=COL["HSPC"] if rv["q16_sev"] < 0.05 else (COL["Ext"] if rv["p_dth11"] < 0.05 else PALETTE["gray_1"]))
    ax_c1.set_yticks(yy_c); ax_c1.set_yticklabels([x[1] for x in HMOD_SHOW], fontsize=8.3)
    for it in (0, 1, 2): ax_c1.get_yticklabels()[it].set_fontweight("bold"); ax_c1.get_yticklabels()[it].set_color(COL["HSPC"] if it < 2 else COL["Ext"])
    ax_c2.set_yticks(yy_c); ax_c2.set_yticklabels([])
    for axx in (ax_c1, ax_c2):
        axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(-0.55, len(HMOD_SHOW) - 0.35)
    ax_c1.set_xlim(-0.58, 0.78); ax_c1.set_xticks([-0.4, 0.0, 0.4]); ax_c1.set_xlabel("Severity ρ (95% CI)", fontsize=8.8)
    ax_c2.set_xlim(0.10, 1.00); ax_c2.set_xticks([0.25, 0.50, 0.75]); ax_c2.set_xlabel("Death AUC (11d)", fontsize=8.8)
    ax_c1.set_title("HSPC transcriptome modules (n = 62; *q16 < 0.05)", loc="left", fontsize=10.2, fontweight="bold", pad=10)

    # Panel 8d: Patient-level HSPC Ery5 across WHO 0-3 (left) & Survived vs Died (right)
    ax_d1 = fig8.add_axes([0.525, 0.530, 0.105, 0.170]); ax_d2 = fig8.add_axes([0.642, 0.530, 0.070, 0.170])
    fig8.text(0.492, 0.725, "d", fontsize=17.5, fontweight="bold", va="top")
    rng_d8 = np.random.default_rng(28)
    for g in range(4):
        for arm_v, k_s, dx in [("Real", "Real", -0.17), ("PBISC new", "New", +0.17)]:
            vals = HSPC_PT62[(HSPC_PT62["arm"] == arm_v) & (HSPC_PT62["sev"] == g)]["Ery5_ns"].values
            ax_d1.scatter(g + dx + rng_d8.uniform(-0.07, 0.07, len(vals)), vals, s=22, marker=MK[k_s], c=COL[k_s], edgecolors="white", lw=0.4, zorder=3)
            if len(vals): ax_d1.plot([g + dx - 0.13, g + dx + 0.13], [np.median(vals)] * 2, color="black", lw=1.7, zorder=4)
    ax_d1.set_xticks(range(4)); ax_d1.set_xticklabels(["0", "1", "2", "3"], fontsize=8.8); ax_d1.set_xlim(-0.55, 3.55)
    ax_d1.set_ylim(-2.5, 3.85); ax_d1.set_yticks([-2, 0, 2]); ax_d1.set_xlabel("Peak WHO", fontsize=8.8); ax_d1.set_ylabel("HSPC Ery5 (Normal score)", fontsize=8.8)
    ery_c = HSPC_MOD16[(HSPC_MOD16["module"] == "Ery5") & (HSPC_MOD16["cohort"] == "Combined")].iloc[0]
    ax_d1.text(0.03, 0.97, f"Comb: ρ={ery_c['rho_sev']:+.2f}\nP={ery_c['p_sev']:.4f} (q16=.044*)", transform=ax_d1.transAxes, fontsize=7.2, fontweight="bold", color=COL["HSPC"], va="top")
    ax_d1.set_title("HSPC Ery5 score (n = 62)", loc="left", fontsize=10.3, fontweight="bold", pad=10)

    for g_f in [0, 1]:
        for arm_v, k_s, dx in [("Real", "Real", -0.17), ("PBISC new", "New", +0.17)]:
            vals = HSPC_PT62[(HSPC_PT62["arm"] == arm_v) & (HSPC_PT62["fat11"] == g_f)]["Ery5_ns"].values
            ax_d2.scatter(g_f + dx + rng_d8.uniform(-0.07, 0.07, len(vals)), vals, s=22, marker=MK[k_s], c=COL[k_s], edgecolors="white", lw=0.4, zorder=3)
            if len(vals): ax_d2.plot([g_f + dx - 0.13, g_f + dx + 0.13], [np.median(vals)] * 2, color="black", lw=1.7, zorder=4)
    ax_d2.set_xticks([0, 1]); ax_d2.set_xticklabels(["Surv.\n(51)", "Died\n(11)"], fontsize=8.5); ax_d2.set_xlim(-0.55, 1.55)
    ax_d2.set_ylim(-2.5, 3.85); ax_d2.set_yticks([-2, 0, 2]); ax_d2.set_yticklabels([])
    ax_d2.text(0.03, 0.97, f"AUC={ery_c['AUC_dth11']:.2f}\nP={ery_c['p_dth11']:.4f}", transform=ax_d2.transAxes, fontsize=7.2, fontweight="bold", color=COL["HSPC"], va="top")

    # Panel 8e: Point 6 — T-cell Type-I IFN / ISG across compartments (Shared-27 fidelity + Severity + Death!)
    ax_e1 = fig8.add_axes([0.805, 0.530, 0.056, 0.170]); ax_e2 = fig8.add_axes([0.869, 0.530, 0.056, 0.170]); ax_e3 = fig8.add_axes([0.933, 0.530, 0.056, 0.170])
    fig8.text(0.728, 0.725, "e", fontsize=17.5, fontweight="bold", va="top")
    TSUB_SHOW = ["T CD8 CTL", "T CD4 CTL", "All CD8 T", "All T cells", "All CD4 T", "T CD8 CM", "T CD4 Naive", "T reg"]
    t_ix = TISG13.set_index("subset")
    yy_e = np.arange(len(TSUB_SHOW))[::-1]
    ax_e1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.0); ax_e2.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.0); ax_e3.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.0)
    for yi, ts in zip(yy_e, TSUB_SHOW):
        rt = t_ix.loc[ts]
        # Shared-27 fidelity
        ax_e1.plot([rt["sh_lo"], rt["sh_hi"]], [yi, yi], color=COL["Shared"], lw=1.5)
        ax_e1.scatter(rt["sh_rho"], yi, s=34, marker="o", c=COL["Shared"] if rt["sh_p"] < 0.05 else "white", edgecolors=COL["Shared"], lw=1.3, zorder=4)
        # Severity rho (q13 < 0.05 = solid blue)
        c_sv = COL["Ext"] if rt["q13_sev"] < 0.05 else PALETTE["gray_2"]
        ax_e2.plot([rt["c_sev_lo"], rt["c_sev_hi"]], [yi, yi], color=c_sv, lw=1.5)
        ax_e2.scatter(rt["c_sev"], yi, s=36, marker="s", c=c_sv if rt["q13_sev"] < 0.05 else "white", edgecolors=c_sv, lw=1.3, zorder=4)
        # Death AUC (p < 0.05 = solid red)
        c_dt = COL["HSPC"] if rt["cp_auc11"] < 0.05 else PALETTE["gray_2"]
        ax_e3.plot([rt["c_auc11_lo"], rt["c_auc11_hi"]], [yi, yi], color=c_dt, lw=1.5)
        ax_e3.scatter(rt["c_auc11"], yi, s=36, marker="s", c=c_dt if rt["cp_auc11"] < 0.05 else "white", edgecolors=c_dt, lw=1.3, zorder=4)
    ax_e1.set_yticks(yy_e); ax_e1.set_yticklabels(TSUB_SHOW, fontsize=8.0)
    for it in range(5): ax_e1.get_yticklabels()[it].set_fontweight("bold"); ax_e1.get_yticklabels()[it].set_color(COL["Ext"])
    for axx in (ax_e2, ax_e3): axx.set_yticks(yy_e); axx.set_yticklabels([])
    for axx in (ax_e1, ax_e2, ax_e3): axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(-0.55, len(TSUB_SHOW) - 0.35)
    ax_e1.set_xlim(-0.15, 0.98); ax_e1.set_xticks([0.0, 0.8]); ax_e1.set_xlabel("Sh27 ρ", fontsize=8.2)
    ax_e2.set_xlim(-0.20, 0.62); ax_e2.set_xticks([0.0, 0.4]); ax_e2.set_xlabel("Sev. ρ*", fontsize=8.2)
    ax_e3.set_xlim(0.20, 0.95); ax_e3.set_xticks([0.5, 0.8]); ax_e3.set_xlabel("Dth AUC", fontsize=8.2)
    ax_e1.set_title("T-cell ISG (19g; *q13=.037)", loc="left", fontsize=10.1, fontweight="bold", pad=10)

    # Row 3: Panel 8f (3 longitudinal trajectories: T CD8 CTL, HSPC, T CD8 CTL ISG!) & Panel 8g (PBISC vs BayesPrism Z!)
    ax_f1 = fig8.add_axes([0.048, 0.290, 0.188, 0.165])
    ax_f2 = fig8.add_axes([0.264, 0.290, 0.188, 0.165])
    ax_f3 = fig8.add_axes([0.480, 0.290, 0.188, 0.165])
    fig8.text(0.012, 0.480, "f", fontsize=17.5, fontweight="bold", va="top")

    x_grid = np.linspace(0, 28, 100)
    gee_h = GEE_SRC[GEE_SRC["feature"] == "HSPC"].set_index("cohort")
    gee_c = GEE_SRC[GEE_SRC["feature"] == "T CD8 CTL"].set_index("cohort")
    gee_i = TISG_GEE[TISG_GEE["subset"] == "T CD8 CTL"].set_index("cohort")
    TRAJ_SPECS = [
        (ax_f1, pos_217, "T CD8 CTL_ns", "fatal_pub11", "T CD8 CTL: early depletion", True,
         f"0–7d AUC = 0.24 (P = 0.011)\nGEE Time×Fatal P = {gee_c.loc['Combined','P']:.3f} (R .017; N .50)"),
        (ax_f2, pos_217, "HSPC_ns", "fatal_pub11", "HSPC: late expansion (>7d)", False,
         f">7d AUC = 0.82 (P = 0.011)\nPaired Δ(43) = 0.77; GEE P = {gee_h.loc['Combined','P']:.3f}"),
        (ax_f3, CD8_ISG_217, "isg_z", "fatal_pub11", "T CD8 CTL ISG: steep fatal collapse", False,
         f"Day-adj Sev P=.0016 · Dth P=.039\nTime×Fatal β=−0.55 (P=3.7e-6; R 9e-5, N .016)"),
    ]
    for ax_f, df_t, ycol, fcol, ttl_f, show_yl, stat_txt in TRAJ_SPECS:
        ax_f.axvspan(7.0, 30.5, color=COL["SoftRed"] if "HSPC" in ttl_f else "#F4F6F9", alpha=0.35, zorder=0)
        ax_f.axvline(7.0, color=COL["HSPC"] if "HSPC" in ttl_f else PALETTE["gray_1"], ls="--", lw=1.2, zorder=1)
        dm = df_t[["patient_id", "source", "days", ycol, fcol]].dropna().rename(columns={ycol: "y", fcol: "fatal"})
        for fat_v, c_ln, lbl_g, al_pt in [(0, COL["New"], "Survived (55 pts)", 0.35), (1, COL["HSPC"], "Died (11 verified)", 0.82)]:
            sub = dm[dm["fatal"] == fat_v]
            for s_k, mk_v in [("Real", "D"), ("PBISC", "o")]:
                sb_s = sub[sub["source"] == s_k]
                ax_f.scatter(sb_s["days"], sb_s["y"], s=16 if fat_v == 0 else 25, marker=mk_v, c=c_ln, alpha=al_pt, edgecolors="none", zorder=2 if fat_v == 0 else 3)
            ols_m = sm.OLS(sub["y"], sm.add_constant(sub["days"])).fit(cov_type="cluster", cov_kwds={"groups": sub["patient_id"]})
            pr = ols_m.get_prediction(sm.add_constant(x_grid)).summary_frame()
            ax_f.fill_between(x_grid, pr["mean_ci_lower"], pr["mean_ci_upper"], color=c_ln, alpha=0.16, zorder=3)
            ax_f.plot(x_grid, pr["mean"], color=c_ln, lw=2.3, label=lbl_g, zorder=4)
        ax_f.text(0.025, 0.97, stat_txt, transform=ax_f.transAxes, fontsize=6.9, fontweight="bold", color=COL["HSPC"] if ("HSPC" in ttl_f or "ISG" in ttl_f) else COL["Ext"], va="top",
                  bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.88))
        ax_f.set_xlim(-0.8, 30.5); ax_f.set_ylim(-2.6, 3.7); ax_f.set_xticks([0, 7, 14, 21, 28]); ax_f.set_yticks([-2, 0, 2])
        ax_f.set_xlabel("Days from first positive swab", fontsize=9.2)
        if show_yl:
            ax_f.set_ylabel("Normal score (longitudinal draws)", fontsize=9.2); ax_f.legend(loc="lower right", frameon=False, fontsize=7.6, handlelength=1.3)
        ax_f.set_title(ttl_f, loc="left", fontsize=10.0, fontweight="bold", pad=6)

    # Panel 8g: PBISC vs BayesPrism Posterior Z cell-type transcriptome recovery!
    ax_g8 = fig8.add_axes([0.750, 0.290, 0.232, 0.165]); fig8.text(0.696, 0.480, "g", fontsize=17.5, fontweight="bold", va="top")
    TGTS_G = ["T CD8 CTL ISG", "T CD8 CM ISG", "HSPC Ery5"]
    METR_G = [("shared_rho", "Sh27 ρ", False), ("new29_sev_rho", "New29 Sev ρ", False), ("new29_dth11_auc", "New29 Dth\n2(AUC−.5)", True)]
    xg = np.arange(len(TGTS_G) * 3)
    pb_vals, bp_vals, x_lbls = [], [], []
    for t_nm in TGTS_G:
        r_pb = BPZ_COMP[(BPZ_COMP["target"] == t_nm) & (BPZ_COMP["method"] == "PBISC")].iloc[0]
        r_bp = BPZ_COMP[(BPZ_COMP["target"] == t_nm) & (BPZ_COMP["method"] == "BayesPrism Z")].iloc[0]
        for mc, ml, is_auc in METR_G:
            vp = 2.0 * (r_pb[mc] - 0.5) if is_auc else r_pb[mc]
            vb = 2.0 * (r_bp[mc] - 0.5) if is_auc else r_bp[mc]
            pb_vals.append(vp); bp_vals.append(vb)
            x_lbls.append(ml.split("\n")[0])
    w_g = 0.36
    ax_g8.axhline(0.0, color="black", lw=1.0)
    ax_g8.bar(xg - w_g / 2, pb_vals, width=w_g, color=COL["Ext"], edgecolor="black", lw=0.8, label="PBISC scRNA")
    ax_g8.bar(xg + w_g / 2, bp_vals, width=w_g, color=PALETTE["red_2"], edgecolor="black", lw=0.8, hatch="///", label="BayesPrism Z")
    for dv in [2.5, 5.5]: ax_g8.axvline(dv, color=PALETTE["neutral"], ls="--", lw=1.1)
    for idx_t, t_nm in enumerate(TGTS_G):
        ax_g8.text(idx_t * 3 + 1.0, 0.90, t_nm, ha="center", va="top", fontsize=7.8, fontweight="bold", color=COL["HSPC"] if "HSPC" in t_nm else COL["Ext"])
    ax_g8.set_xticks(xg); ax_g8.set_xticklabels(x_lbls, rotation=38, ha="right", fontsize=7.2)
    ax_g8.set_ylim(-0.78, 0.98); ax_g8.set_yticks([-0.5, 0.0, 0.5, 0.8])
    ax_g8.set_ylabel("Fidelity ρ / effect", fontsize=8.8)
    ax_g8.legend(loc="lower right", frameon=False, fontsize=7.6, handlelength=1.2)
    ax_g8.set_title("Cell-type transcriptome: PBISC vs BayesPrism Z", loc="left", fontsize=9.8, fontweight="bold", pad=6)

    # Row 4: Panel 8h (Hospital bloods, 5 verified Real deaths) & Panel 8i (500k In-fold CV under TRUE_PUB_11!)
    ax_h1 = fig8.add_axes([0.118, 0.048, 0.132, 0.155]); ax_h2 = fig8.add_axes([0.270, 0.048, 0.132, 0.155])
    fig8.text(0.012, 0.228, "h", fontsize=17.5, fontweight="bold", va="top")
    HOSP_LABS = [("temp_CRP", "CRP (37, 5d)"), ("temp_lymph", "ALC (37, 5d)"), ("temp_ferritin", "Ferritin (33, 5d)"), ("temp_ddimer", "D-dimer (31, 5d)"), ("temp_trop", "Troponin (24, 4d)")]
    yy_h = np.arange(len(HOSP_LABS))[::-1]
    ax_h1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3); ax_h2.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
    hspc_lab = LAB_IND[LAB_IND["feature"] == "HSPC"].set_index("lab")
    for yi, (lk, _) in zip(yy_h, HOSP_LABS):
        r_l = hspc_lab.loc[lk]
        ax_h1.barh(yi, r_l["rho_with_lab"], height=0.55, color="white", edgecolor=PALETTE["gray_2"], lw=1.4)
        ax_h1.text(r_l["rho_with_lab"] + (0.04 if r_l["rho_with_lab"] >= 0 else -0.04), yi, f"{r_l['rho_with_lab']:+.2f} (ns)", va="center", ha="left" if r_l["rho_with_lab"] >= 0 else "right", fontsize=7.5, color=PALETTE["gray_1"])
        sig_l = r_l["resid_p11"] < 0.05
        ax_h2.plot([0.50, r_l["resid_auc11"]], [yi, yi], color=COL["HSPC"] if sig_l else PALETTE["gray_2"], lw=1.8)
        ax_h2.scatter(r_l["resid_auc11"], yi, s=44, marker="D", c=COL["HSPC"] if sig_l else "white", edgecolors=COL["HSPC"], linewidths=1.3, zorder=4)
        ax_h2.text(r_l["resid_auc11"] + 0.025, yi, f"{r_l['resid_auc11']:.2f} (P={r_l['resid_p11']:.2f})", va="center", ha="left", fontsize=7.4, color=PALETTE["gray_2"])
    ax_h1.set_yticks(yy_h); ax_h1.set_yticklabels([x[1] for x in HOSP_LABS], fontsize=8.6, fontweight="bold"); ax_h2.set_yticks(yy_h); ax_h2.set_yticklabels([])
    for axx in (ax_h1, ax_h2): axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(-0.55, len(HOSP_LABS) - 0.35)
    ax_h1.set_xlim(-0.54, 0.64); ax_h1.set_xticks([-0.4, 0.0, 0.4]); ax_h1.set_xlabel("ρ(HSPC, mean lab)", fontsize=9.0)
    ax_h1.set_title("HSPC vs. patient-mean bloods (Real, 5 verified deaths)", loc="left", fontsize=10.1, fontweight="bold", pad=8)
    ax_h2.set_xlim(0.46, 1.04); ax_h2.set_xticks([0.5, 0.75, 1.0]); ax_h2.set_xlabel("Rank-resid. HSPC death AUC", fontsize=9.0)

    ax_i8 = fig8.add_axes([0.652, 0.048, 0.175, 0.155]); fig8.text(0.455, 0.228, "i", fontsize=17.5, fontweight="bold", va="top")
    CV_SPECS = [
        (0, "First available draw (n=66, 11d)\nAge + Sex + Same-day WHO → + CD8 CTL", COL["Ext"]),
        (1, "First draw ≤ 7d subset (n=53, 10d)\nAge + Sex + Same-day WHO → + CD8 CTL", COL["Ext"]),
        (2, "Hospital bloods (Real n=37, 5d)\nAge + CRP + ALC → + HSPC", COL["HSPC"]),
        (3, "Routine labs only (Real n=37, 5d)\nCRP + ALC → + HSPC", COL["HSPC"]),
        (4, "Severe/Critical WHO 2–3 (n=29, 11d)\nAge + Sex + WHO 3 vs 2 → + HSPC", COL["HSPC"]),
    ]
    yy_i = np.arange(len(CV_SPECS))[::-1]
    ax_i8.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1)
    for yi, (r_idx, _, pcol) in zip(yy_i, CV_SPECS):
        rv = CV_500K.iloc[r_idx]
        d_lo, d_hi, d_v = float(rv["delta_lo"]), float(rv["delta_hi"]), float(rv["Delta_AUC"])
        excl_zero = d_lo > 0.0
        ax_i8.plot([rv["base_lo"], rv["base_hi"]], [yi + 0.16, yi + 0.16], color=PALETTE["gray_1"], lw=1.5, zorder=2)
        ax_i8.scatter(rv["AUC_base"], yi + 0.16, s=38, marker="D", c="white", edgecolors=PALETTE["gray_2"], linewidths=1.4, zorder=3, label="Clinical / lab baseline" if yi == yy_i[0] else "")
        ax_i8.plot([rv["ext_lo"], rv["ext_hi"]], [yi - 0.16, yi - 0.16], color=pcol, lw=1.9, zorder=3)
        ax_i8.scatter(rv["AUC_ext"], yi - 0.16, s=46, marker="s", c=pcol if excl_zero else "white", edgecolors=pcol, linewidths=1.5, zorder=4, label="+ Immune biomarker" if yi == yy_i[0] else "")
        ax_i8.text(1.03, yi, f"Δ={d_v:+.2f} [{d_lo:+.2f}, {d_hi:+.2f}]", va="center", ha="left", fontsize=7.4, fontweight="bold" if excl_zero else "normal", color=pcol if excl_zero else PALETTE["gray_2"])
    ax_i8.text(1.03, yy_i[0] + 0.65, "Paired ΔCV-AUC [500k 95% CI]", va="bottom", ha="left", fontsize=7.3, fontweight="bold", color=PALETTE["gray_2"])
    ax_i8.set_yticks(yy_i); ax_i8.set_yticklabels([x[1] for x in CV_SPECS], fontsize=7.95, fontweight="bold")
    for idx_y in (2, 3, 4): ax_i8.get_yticklabels()[idx_y].set_color(PALETTE["gray_2"])
    ax_i8.tick_params(axis="y", length=0); ax_i8.spines["left"].set_visible(False)
    ax_i8.set_xlim(0.32, 1.01); ax_i8.set_xticks([0.4, 0.5, 0.7, 0.9]); ax_i8.set_ylim(-0.58, len(CV_SPECS) - 0.30)
    ax_i8.set_xlabel("Out-of-fold CV AUC (in-fold INT, 50 × 5-fold)", fontsize=9.1)
    ax_i8.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.05), ncol=2, frameon=False, fontsize=7.6, columnspacing=0.8, handlelength=1.1)
    ax_i8.set_title("In-fold preprocessed mortality CV (TRUE 11 deaths)", loc="left", fontsize=10.3, fontweight="bold", pad=18)

    saved_8 = finalize_figure(fig8, out_dir / "PBISC_Figure8_Main_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Main Figure 8:", saved_8[0].name)

    # =====================================================================
    # 3) SUPPLEMENTARY FIGURE 1 (a - d) — COHORT + MORTALITY FORENSIC AUDIT
    # =====================================================================
    fig_s1 = plt.figure(figsize=(16.8, 12.6), facecolor="white")
    fig_s1.text(0.015, 0.986, "Supplementary Figure 1. Longitudinal sampling, Point-1 mortality forensic audit (C101/C146/C147), and analysis architecture",
                fontsize=13.5, fontweight="bold", va="top", ha="left", color=COL["Real"])

    ax_s1a = fig_s1.add_axes([0.140, 0.560, 0.215, 0.345]); fig_s1.text(0.015, 0.938, "a", fontsize=18, fontweight="bold", va="top")
    ax_s1a.axvspan(7.0, 31.5, color=COL["CardRed"], alpha=0.38, zorder=0); ax_s1a.axvline(7.0, color=COL["CLL"], ls="--", lw=1.3, zorder=2)
    ax_s1a.text(19.2, 4.18, "Late window (> day 7)\nHSPC mortality separation", ha="center", va="top", fontsize=8.2, fontweight="bold", color=COL["CLL"])
    ax_s1a.text(3.3, 4.18, "Early\n(≤7d)", ha="center", va="top", fontsize=8.2, color=PALETTE["gray_2"])
    rng_c = np.random.default_rng(12)
    groups_c = [
        ("PBISC new · Died (6)", b_new_samps[b_new_samps["fatal_pub11"] == 1]["days_from_first_positive_swab"], 0.0, "o", COL["CLL"]),
        ("PBISC new · Surv. (23)", b_new_samps[b_new_samps["fatal_pub11"] == 0]["days_from_first_positive_swab"], 1.0, "o", COL["New"]),
        ("Real · Died (5)", r_samps[r_samps["fatal_pub11"] == 1]["days_from_first_positive_swab"], 2.1, "D", COL["CLL"]),
        ("Real · Surv. (32)", r_samps[r_samps["fatal_pub11"] == 0]["days_from_first_positive_swab"], 3.1, "D", COL["Real"]),
    ]
    for lab, s_days, y0, mk, c_dot in groups_c:
        dv = s_days.dropna().values
        ax_s1a.scatter(dv, y0 + rng_c.uniform(-0.18, 0.18, len(dv)), s=25, marker=mk, c=c_dot, edgecolors="white", linewidths=0.45, alpha=0.9, zorder=3)
        q25, med, q75 = np.percentile(dv, [25, 50, 75])
        ax_s1a.plot([q25, q75], [y0 - 0.30, y0 - 0.30], color="black", lw=2.0, zorder=4)
        ax_s1a.scatter([med], [y0 - 0.30], s=30, marker="|", c="black", linewidths=2.2, zorder=5)
    ax_s1a.set_yticks([g[2] for g in groups_c]); ax_s1a.set_yticklabels([g[0] for g in groups_c], fontsize=9.0)
    ax_s1a.tick_params(axis="y", length=0, pad=5); ax_s1a.spines["left"].set_visible(False)
    ax_s1a.set_xlim(-1.0, 31.5); ax_s1a.set_ylim(-0.55, 4.30); ax_s1a.set_xticks([0, 7, 14, 21, 30])
    ax_s1a.set_xlabel("Days from first positive swab", fontsize=10.0)
    ax_s1a.set_title("Longitudinal sampling (217 draws)", loc="left", fontsize=11.0, fontweight="bold", pad=8)

    ax_s1b = fig_s1.add_axes([0.460, 0.560, 0.155, 0.345]); fig_s1.text(0.376, 0.938, "b", fontsize=18, fontweight="bold", va="top")
    r_cnts, n_cnts = r_samps.groupby("patient_id").size(), b_new_samps.groupby("patient_id").size()
    r_pts_meta = r_samps.groupby("patient_id").agg(batch=("centre_or_wave", "first"), ster=("steroid_real_Y", "max"))
    cats_s1b = [
        ("1 draw", int((r_cnts == 1).sum()), int((n_cnts == 1).sum())),
        ("2 draws", int((r_cnts == 2).sum()), int((n_cnts == 2).sum())),
        ("3 draws", int((r_cnts == 3).sum()), int((n_cnts == 3).sum())),
        ("≥4 draws", int((r_cnts >= 4).sum()), int((n_cnts >= 4).sum())),
        ("NCL centre", int((r_pts_meta["batch"] == "NCL").sum()), 0),
        ("Cambridge", int((r_pts_meta["batch"] == "Cambridge").sum()), 0),
        ("Wave-1 bulk", 0, len(n_cnts)),
        ("Steroid (Real)", int(r_pts_meta["ster"].sum()), 0),
    ]
    yb_s1 = np.arange(len(cats_s1b))[::-1]; yb_s1_adj = np.where(yb_s1 < 4, yb_s1 - 0.40, yb_s1)
    ax_s1b.plot([0, 34], [3.78, 3.78], color=PALETTE["neutral"], ls="--", lw=1.1, zorder=1)
    for yi, (lab, nr, nn) in zip(yb_s1_adj, cats_s1b):
        if nr > 0: ax_s1b.barh(yi, nr, height=0.62, color=COL["Real"], edgecolor="black", lw=1.0, label="Real (37)" if yi == yb_s1_adj[0] else "")
        if nn > 0: ax_s1b.barh(yi, nn, left=nr, height=0.62, color=COL["New"], edgecolor="black", lw=1.0, label="PBISC new (29)" if yi == yb_s1_adj[0] else "")
        ax_s1b.text(nr + nn + 0.7, yi, f"{nr + nn}", va="center", ha="left", fontsize=8.7, fontweight="bold", color=COL["Real"])
    ax_s1b.set_yticks(yb_s1_adj); ax_s1b.set_yticklabels([c[0] for c in cats_s1b], fontsize=8.8)
    ax_s1b.tick_params(axis="y", length=0, pad=5); ax_s1b.spines["left"].set_visible(False)
    ax_s1b.set_xlim(0, 36.0); ax_s1b.set_xticks([0, 10, 20, 30]); ax_s1b.set_xlabel("COVID-19+ patients", fontsize=10.0)
    ax_s1b.set_title("Draws & centre strata", loc="left", fontsize=11.0, fontweight="bold", pad=8)
    ax_s1b.legend(loc="lower right", frameon=False, fontsize=8.2, handlelength=1.1, handletextpad=0.35)

    # Panel S1c: Point-1 Mortality Forensic Audit Matrix (ALL 14 Zenodo-flagged patients!)
    ax_s1c = fig_s1.add_axes([0.708, 0.548, 0.275, 0.352]); fig_s1.text(0.638, 0.938, "c", fontsize=18, fontweight="bold", va="top")
    COLS_M = ["Table S1\n(2025)", "Zenodo\nw1/w2", "2-mo\nRecov.", "TRUE\n11d", "Legacy\n13d"]
    for i_m, (_, rm) in enumerate(MORT_AUD.iterrows()):
        s1_t = "yes" if rm["stephenson2025_TableS1_fatal"] == "yes" else ("indir" if rm["stephenson2025_TableS1_fatal"] == "indirect" else ("no" if rm["stephenson2025_TableS1_fatal"] == "no" else "–"))
        s1_c = COL["CLL"] if s1_t in ("yes", "indir") else (COL["ScoreGood"] if s1_t == "no" else "#EAEAEA")
        zen_t = "TRUE"
        rec_t = "YES" if bool(rm["has_2month_recovery_sample"]) else "–"
        rec_c = COL["ScoreGood"] if rec_t == "YES" else "white"
        tr_t = "DIED" if int(rm["TRUE_PUBLISHED_11_PRIMARY"]) == 1 else "ALIVE"
        tr_c = COL["CLL"] if tr_t == "DIED" else COL["ScoreGood"]
        lg_t = "DIED" if int(rm["legacy_13_hybrid"]) == 1 else "ALIVE"
        lg_c = COL["ScoreMid"] if (tr_t == "ALIVE" and lg_t == "DIED") else ("#E5E7EB" if tr_t == lg_t else COL["ScoreGood"])
        row_c = [(s1_t, s1_c, "white" if s1_c != "#EAEAEA" else COL["Real"]), (zen_t, "#FADBD8", COL["CLL"]),
                 (rec_t, rec_c, "white" if rec_t == "YES" else PALETTE["gray_1"]), (tr_t, tr_c, "white"),
                 (lg_t, lg_c, "white" if lg_c == COL["ScoreMid"] else COL["Real"])]
        for j_m, (tx, fc, tc) in enumerate(row_c):
            ax_s1c.add_patch(FancyBboxPatch((j_m + 0.06, i_m + 0.10), 0.88, 0.80, boxstyle="round,pad=0.015,rounding_size=0.08", fc=fc, ec="#333333", lw=0.4))
            ax_s1c.text(j_m + 0.5, i_m + 0.5, tx, ha="center", va="center", fontsize=6.9, color=tc, fontweight="bold")
    ax_s1c.set_xlim(0, len(COLS_M)); ax_s1c.set_ylim(len(MORT_AUD), 0)
    ax_s1c.set_yticks(np.arange(len(MORT_AUD)) + 0.5)
    ax_s1c.set_yticklabels([f"{r.patient_id} ({'Sh' if 'Shared' in r.cohort else ('R' if 'Real' in r.cohort else 'N')})" for _, r in MORT_AUD.iterrows()], fontsize=7.8, fontweight="bold")
    for idx_y, r in enumerate(MORT_AUD.itertuples()):
        if r.patient_id in {"C101", "C146", "C147"}:
            ax_s1c.get_yticklabels()[idx_y].set_color(COL["ScoreGood"])
    ax_s1c.set_xticks(np.arange(len(COLS_M)) + 0.5); ax_s1c.set_xticklabels(COLS_M, fontsize=7.6, fontweight="bold")
    ax_s1c.xaxis.tick_top(); ax_s1c.tick_params(length=0)
    for sp_m in ax_s1c.spines.values(): sp_m.set_visible(False)
    ax_s1c.set_title("Point-1 forensic audit: 3 w2 survivors fixed", loc="left", fontsize=10.3, fontweight="bold", pad=22)

    ax_s1d = fig_s1.add_axes([0.020, 0.040, 0.962, 0.395]); ax_s1d.set_xlim(0, 100); ax_s1d.set_ylim(0, 100); ax_s1d.axis("off")
    fig_s1.text(0.015, 0.455, "d", fontsize=18, fontweight="bold", va="top")
    ax_s1d.text(0.2, 97.2, "End-to-end verified workflow: Published-11 ground truth (Fig. 7), 10-method deconvolution & single-cell transcriptomics (Fig. 8), and external replication (Supp Fig. 4)",
                fontsize=10.7, fontweight="bold", va="top", ha="left", color=COL["Real"])
    arr_kw_s1 = dict(arrowstyle="-|>", lw=2.0, color=PALETTE["gray_2"], mutation_scale=16)
    rbox(ax_s1d, 0.4, 7, 29.0, 81, "#F8FAFC", COL["Ext"], lw=1.7)
    ax_s1d.text(14.9, 81.2, "1. Verified ESKD Screen (n = 66)", ha="center", va="center", fontsize=11.0, fontweight="bold", color=COL["Ext"])
    ax_s1d.text(2.1, 43.0, "• WHO 0–3: 100% exact (0/68 diff.)\n• Mortality: 11 verified COVID deaths\n  (Real 5 + PBISC new 6; C101/C146/\n  C147 fixed as survivors)\n• 5-step screen (fixed BH43):\n  43 → 13 → 12 → 11 → 6 PBISC-enabled",
                ha="left", va="center", fontsize=10.1, color=COL["Real"], linespacing=1.36)
    ax_s1d.annotate("", xy=(31.8, 47.5), xytext=(29.8, 47.5), arrowprops=arr_kw_s1)

    rbox(ax_s1d, 32.2, 7, 33.6, 81, "white", COL["Real"], lw=1.7)
    ax_s1d.text(49.0, 81.2, "2. Prioritized Cellular & Transcriptomic Signals", ha="center", va="center", fontsize=10.8, fontweight="bold", color=COL["Real"])
    rbox(ax_s1d, 33.3, 52.5, 31.4, 21.0, COL["CardGray"], COL["Real"], lw=1.1, r=1.4)
    ax_s1d.text(49.0, 63.0, "Tier 1 · Real-sig Anchors (5) + PBISC-enabled (6 pairs)\nSeverity: B ASC IgG/IgA (↑), HSPC (↑), T CD8 CTL (↓)\nDeath: T CD8 CTL (↓, 0.21; CV Δ=+0.14), HSPC (↑, 0.76)",
                ha="center", va="center", fontsize=8.4, fontweight="bold", color=COL["Real"], linespacing=1.30)
    rbox(ax_s1d, 33.3, 27.5, 31.4, 21.5, COL["CardBlue"], COL["Ext"], lw=1.2, r=1.4)
    ax_s1d.text(49.0, 38.2, "Points 4–6 · Deconvolution & Cell-Type Transcriptome\nPBISC 10/10 vs BayesPrism/CIBERSORTx/MuSiC (2–5/10)\nHSPC Ery5 sev q16=.044* · T-cell ISG sev q13=.037*, collapse P=3.7e-6",
                ha="center", va="center", fontsize=8.2, fontweight="bold", color=COL["Ext"], linespacing=1.30)
    rbox(ax_s1d, 33.3, 10.5, 31.4, 13.5, COL["CardRed"], COL["CLL"], lw=1.0, r=1.2)
    ax_s1d.text(49.0, 17.2, "Honest boundaries: HSPC WHO 2–3 (P=.13), Src×WHO (P=.45), lab CV null",
                ha="center", va="center", fontsize=8.3, color=COL["CLL"], fontweight="bold")
    ax_s1d.annotate("", xy=(68.2, 47.5), xytext=(66.2, 47.5), arrowprops=arr_kw_s1)

    rbox(ax_s1d, 68.6, 7, 30.9, 81, COL["CardTeal"], COL["Shared"], lw=1.8)
    ax_s1d.text(84.05, 80.2, "3. External Adult COVID-19 Check\n(4 Cohorts · N_sev=220 · N_dth=219 [24d])",
                ha="center", va="center", fontsize=10.3, fontweight="bold", color="#165057", linespacing=1.22)
    ax_s1d.text(69.9, 41.0, "• Stephenson 2021 (85 acute; 74 sev / 73 dth [7d])\n• COMBAT 2022 (Wave-1, n=75, 7 deaths)\n• van der Wijst 2021 (UCSF, n=42, 7 deaths)\n• Liu 2021 (NIAID, n=29, 3 deaths)\n→ Anchors, B ASC total & HSPC (DL; mKH wider)\n→ Deprioritized NK CD56dim & CD16 C1 null",
                ha="left", va="center", fontsize=9.8, color=COL["Real"], linespacing=1.38)

    saved_s1 = finalize_figure(fig_s1, out_dir / "PBISC_SuppFig_S1_Cohort_Mortality_Audit_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Supplementary Figure 1:", saved_s1[0].name)

    # =====================================================================
    # 4) SUPPLEMENTARY FIGURE 2 (a - o) — CLL, ANCHORS, NK, MODALITY & RAW % (TRUE_PUB_11!)
    # =====================================================================
    def draw_s2_cll_bars(ax):
        BR = [("T_pct", "T cell"), ("NK_pct", "NK cell"), ("B_pct", "B cell"), ("Mono_pct", "Monocyte"), ("Other_pct", "Other")]
        sc = S_meta[(S_meta["condition"] == "positive") & S_meta["in_primary_scope"].astype(bool) & S_meta["qc_pass_min200cells"].astype(bool)]
        rows = []
        for src in ["Real", "PBISC"]:
            m = sc[(sc["source"] == src) & (~sc["patient_id"].isin(CLL))].groupby("patient_id")[[b for b, _ in BR]].mean().median()
            rows.append((f"{src}: cohort median", m / m.sum() * 100.0, src))
            for p in sorted(CLL):
                q = sc[(sc["source"] == src) & (sc["patient_id"] == p)]
                if len(q):
                    v = q[[b for b, _ in BR]].mean()
                    rows.append((f"{src}: {p} ({len(q)})", v / v.sum() * 100.0, src))
        for i, (lab, v, src) in enumerate(rows):
            left = 0.0
            for b, ct in BR:
                ax.barh(-i, v[b], left=left, color=CT_COL[ct], edgecolor="white", lw=1.1, height=0.68); left += v[b]
            if "C1" in lab:
                ax.text(101.5, -i, f"B {v['B_pct']:.0f}% · Mono {v['Mono_pct']:.0f}% · NK {v['NK_pct']:.0f}%", va="center", fontsize=9.4,
                        color=PALETTE["red_strong"] if src == "Real" else PALETTE["gray_2"], fontweight="bold" if src == "Real" else "normal")
        ax.set_yticks(-np.arange(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=10.0)
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlim(0, 158); ax.set_xticks([0, 50, 100]); ax.set_xlabel("% of PBMCs (COVID+ samples, patient mean)", fontsize=10.3)
        ax.legend(handles=[Rectangle((0, 0), 1, 1, color=CT_COL[c]) for _, c in BR], labels=[c for _, c in BR],
                  ncol=5, loc="lower left", bbox_to_anchor=(-0.02, 1.01), fontsize=9.4, handlelength=1.0, columnspacing=0.8, handletextpad=0.3)

    def draw_s2_cll_shift(ax, val, lo, hi, pcol, ref, xlab, xlim, xt):
        FEATS_SH = ["B_ASC_dividing", "NK CD56 bright", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC"]
        y, ticks = 0.0, []
        for f in FEATS_SH:
            for st, k in [("Real", "Real"), ("Real + PBISC new", "Ext")]:
                a = CI[(CI["feature"] == f) & (CI["set"] == st) & (~CI["CLL_excluded"])].iloc[0]
                b = CI[(CI["feature"] == f) & (CI["set"] == st) & (CI["CLL_excluded"])].iloc[0]
                ax.plot([b[lo], b[hi]], [y, y], color=COL[k], lw=1.8, alpha=0.55)
                ax.annotate("", xy=(b[val], y), xytext=(a[val], y), arrowprops=dict(arrowstyle="-|>", color=PALETTE["gray_1"], lw=1.2, mutation_scale=10, shrinkA=2, shrinkB=4))
                ax.scatter(a[val], y, s=36, marker=MK[k], c="white" if a[pcol] >= 0.05 else COL["All"], edgecolors=COL["All"], linewidths=1.4, zorder=3)
                ax.scatter(b[val], y, s=50, marker=MK[k], c="white" if b[pcol] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.6, zorder=4)
                y -= 1.0
            ticks.append(y + 1.5); y -= 0.55
        ax.axvline(ref, color=PALETTE["gray_1"], ls=":", lw=1.4)
        ax.set_yticks(ticks); ax.set_yticklabels([LBL[f] for f in FEATS_SH], fontsize=9.8)
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_ylim(y + 0.35, 0.55); ax.set_xlim(*xlim); ax.set_xticks(xt); ax.set_xlabel(xlab, fontsize=10.3)

    def draw_s2_nk_decomp(ax):
        measures = ["CD56bright % PBMC", "CD56bright % of NK", "Total NK % PBMC"]
        xlabels = ["CD56bright\n% PBMC", "CD56bright\n% of NK", "Total NK\n% PBMC"]
        for k_src, col_k in [("Real", "Real"), ("PBISC new", "New")]:
            for i, m_name in enumerate(measures):
                row = NKD[(NKD["source"] == k_src) & (NKD["measure"] == m_name)].iloc[0]
                r, p = row["rho"], row["p"]; x = i + (-0.19 if k_src == "Real" else 0.19)
                ax.bar(x, r, width=0.35, color=COL[col_k] if p < 0.05 else "white", edgecolor=COL[col_k], lw=1.6)
                ax.text(x, r - 0.045 if r < 0 else r + 0.03, "*" if p < 0.05 else "ns", ha="center", va="top" if r < 0 else "bottom", fontsize=9.6, fontweight="bold" if p < 0.05 else "normal")
        ax.axhline(0, color="black", lw=1.0); ax.set_xticks(range(3)); ax.set_xticklabels(xlabels, fontsize=9.4)
        ax.set_ylim(-0.92, 0.28); ax.set_yticks([-0.8, -0.4, 0.0]); ax.set_ylabel("Spearman ρ with severity", fontsize=10.0); ax.tick_params(axis="x", length=0)

    def draw_s2_modality_scale(ax):
        FEATS_M = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
        x = np.arange(len(FEATS_M)); rng_m = np.random.default_rng(19)
        for i, f in enumerate(FEATS_M):
            vr, vb = SH27[f + "|Real_pct"].values, SH27[f + "|PBISC_pct"].values
            for r_i, b_i in zip(vr, vb):
                ax.plot([i - 0.18, i + 0.18], [np.log10(r_i + 0.02), np.log10(b_i + 0.02)], color=PALETTE["neutral"], lw=0.6, alpha=0.55, zorder=1)
            ax.scatter(i - 0.18 + rng_m.uniform(-0.04, 0.04, len(vr)), np.log10(vr + 0.02), s=22, marker="D", c=COL["Real"], edgecolors="white", linewidths=0.35, zorder=3)
            ax.scatter(i + 0.18 + rng_m.uniform(-0.04, 0.04, len(vb)), np.log10(vb + 0.02), s=24, marker="o", c=COL["New"], edgecolors="white", linewidths=0.35, zorder=3)
            ax.plot([i - 0.30, i - 0.06], [np.log10(np.median(vr) + 0.02)] * 2, color="black", lw=2.0, zorder=4)
            ax.plot([i + 0.06, i + 0.30], [np.log10(np.median(vb) + 0.02)] * 2, color="black", lw=2.0, zorder=4)
            mrow = MOD[MOD["feature"] == f].iloc[0]; wp, fold = mrow["wilcoxon_paired_p"], mrow["fold_PBISC_over_Real"]
            star_s = "***" if wp < 0.001 else ("**" if wp < 0.01 else ("*" if wp < 0.05 else "ns"))
            ax.text(i, 1.80, f"{fold:.1f}×\n({star_s})", ha="center", va="top", fontsize=8.2, color=PALETTE["red_strong"] if wp < 0.01 else PALETTE["gray_2"], fontweight="bold" if wp < 0.01 else "normal")
        h_pr, g_pr = PR_REALB.loc["HSPC"], PR_REALB.loc["B_ASC_IgG"]
        ax.text(0.99, 0.03, f"Real B 30-NN P/R: HSPC {100*h_pr['precision']:.1f}%/{100*h_pr['recall']:.1f}%, pooled B ASC 98.4%/98.6% (IgG alone {100*g_pr['precision']:.1f}%/{100*g_pr['recall']:.1f}%)",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=PALETTE["gray_2"],
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=PALETTE["neutral"], lw=0.7, alpha=0.92))
        yt_pct = [0.0, 0.1, 0.5, 2.0, 10.0, 40.0]
        ax.set_yticks([np.log10(v + 0.02) for v in yt_pct]); ax.set_yticklabels(["0%", "0.1%", "0.5%", "2%", "10%", "40%"], fontsize=9.4)
        ax.set_ylim(np.log10(0.012), 2.04); ax.set_xticks(x); ax.set_xticklabels([LBL[f] for f in FEATS_M], rotation=22, ha="right", fontsize=9.5)
        ax.set_ylabel("% of PBMCs (27 shared donors, log scale)", fontsize=9.8)

    def draw_strip_raw_pct(ax, f, outcome, rng_local, show_ylabel=True):
        if outcome == "sev":
            for g in range(4):
                for k, dx in [("Real", -0.18), ("New", 0.18)]:
                    vals = PT_RAW[(PT_RAW["src"] == k) & (PT_RAW["sev"] == g)][f].dropna().values
                    ax.scatter(g + dx + rng_local.uniform(-0.075, 0.075, len(vals)), vals, s=22, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                    if len(vals): ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.8, zorder=4)
            ax.set_xticks(range(4)); ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=9.0); ax.set_xlim(-0.55, 3.55)
            rp_r, pp_r = stats.spearmanr(PT_RAW["sev"].values, PT_RAW[f].values)
            t2 = f"Raw % pool: ρ={rp_r:+.2f}, {fp(pp_r)}"
        else:
            for g in range(2):
                for k, dx in [("Real", -0.18), ("New", 0.18)]:
                    vals = PT_RAW[(PT_RAW["src"] == k) & (PT_RAW["fat11"] == g)][f].dropna().values
                    ax.scatter(g + dx + rng_local.uniform(-0.075, 0.075, len(vals)), vals, s=24, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                    if len(vals): ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.8, zorder=4)
            ax.set_xticks([0, 1]); ax.set_xticklabels(["Surv.\n(55)", "Died\n(11)"], fontsize=9.0); ax.set_xlim(-0.52, 1.52)
            ap_r, pp_r = auc_mw(PT_RAW[f].values, PT_RAW["fat11"].values)
            t2 = f"Raw % pool: AUC={ap_r:.2f}, {fp(pp_r)}"
        vmax = float(PT_RAW[f].max()); ax.set_ylim(-0.04 * vmax, vmax * 1.38)
        if show_ylabel: ax.set_ylabel("% of PBMCs (patient mean)", fontsize=9.6)
        r_real = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == "Real")].iloc[0]
        r_new = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == "PBISC new")].iloc[0]
        t1 = f"Real: ρ={r_real.rho:+.2f} | PBISC: {r_new.rho:+.2f}" if outcome == "sev" else f"Real: AUC={r_real.AUC_death:.2f} | PBISC: {r_new.AUC_death:.2f}"
        ax.text(0.03, 0.97, t1, transform=ax.transAxes, fontsize=7.6, color=PALETTE["gray_2"], va="top")
        ax.text(0.03, 0.86, t2, transform=ax.transAxes, fontsize=7.7, color=COL["Ext"], fontweight="bold", va="top")

    fig_s2 = plt.figure(figsize=(15.2, 20.6), facecolor="white")
    fig_s2.text(0.015, 0.990, "Supplementary Figure 2. CLL exclusion audit, concordant anchors, NK decomposition, and between-modality scale offset (TRUE 11 deaths)",
                fontsize=13.3, fontweight="bold", va="top", ha="left", color=COL["Real"])
    ax_s2a = fig_s2.add_axes([0.125, 0.838, 0.330, 0.110]); draw_s2_cll_bars(ax_s2a); fig_s2.text(0.015, 0.962, "a", fontsize=17, fontweight="bold", va="top")
    ax_s2b = fig_s2.add_axes([0.575, 0.832, 0.165, 0.122]); draw_s2_cll_shift(ax_s2b, "rho", "rho_lo", "rho_hi", "p_sev", 0.0, "ρ with severity (CLL incl. → excl.)", (-0.88, 0.88), [-0.5, 0, 0.5]); fig_s2.text(0.482, 0.962, "b", fontsize=17, fontweight="bold", va="top")
    ax_s2c = fig_s2.add_axes([0.815, 0.832, 0.165, 0.122]); draw_s2_cll_shift(ax_s2c, "AUC", "AUC_lo", "AUC_hi", "p_death", 0.5, "AUC for death (CLL incl. → excl.)", (0.05, 1.02), [0.25, 0.5, 0.75, 1.0]); fig_s2.text(0.755, 0.962, "c", fontsize=17, fontweight="bold", va="top")

    ax_s2_leg = fig_s2.add_axes([0.02, 0.780, 0.96, 0.020]); ax_s2_leg.set_axis_off()
    ax_s2_leg.legend(handles=[
        Line2D([0], [0], marker="D", ls="", ms=7.5, mfc=COL["Real"], mec=COL["Real"], label="Real (37, 5 verified deaths)"),
        Line2D([0], [0], marker="o", ls="", ms=7.5, mfc=COL["New"], mec=COL["New"], label="PBISC new (29, 6d) / PBISC shared (i)"),
        Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["Ext"], mec=COL["Ext"], label="Real + PBISC new (66, 11d, strat. CI)"),
        Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["All"], mec=COL["All"], label="all patients, CLL included (b, c)"),
        Line2D([0], [0], marker="o", ls="", ms=7.5, mfc="white", mec=PALETTE["gray_2"], mew=1.5, label="open: P ≥ 0.05"),
    ], loc="center", ncol=5, fontsize=9.3, frameon=False, columnspacing=1.1)

    rng_s2 = np.random.default_rng(17)
    for idx_p, (f, let) in enumerate([("NK CD56 bright", "d"), ("B_ASC_dividing", "e"), ("T CD8 CM", "f"), ("CD16mono", "g")]):
        x_l = 0.062 + idx_p * 0.242
        ax_an = fig_s2.add_axes([x_l, 0.565, 0.182, 0.165]); draw_strip_ns(ax_an, f, "sev", rng_s2, show_ylabel=(idx_p == 0))
        ax_an.set_title(LBL[f], loc="left", fontsize=11.0, fontweight="bold", color=COL["Real"], pad=6); fig_s2.text(x_l - 0.044, 0.748, let, fontsize=17, fontweight="bold", va="top")

    ax_s2h = fig_s2.add_axes([0.070, 0.318, 0.240, 0.168]); draw_s2_nk_decomp(ax_s2h); ax_s2h.set_title("NK compartment decomposition", loc="left", fontsize=11.0, fontweight="bold", pad=6); fig_s2.text(0.015, 0.505, "h", fontsize=17, fontweight="bold", va="top")
    ax_s2i = fig_s2.add_axes([0.405, 0.318, 0.570, 0.168]); draw_s2_modality_scale(ax_s2i); ax_s2i.set_title("Between-modality differences across 27 shared donors (non-contemporaneous visits)", loc="left", fontsize=10.6, fontweight="bold", pad=6); fig_s2.text(0.345, 0.505, "i", fontsize=17, fontweight="bold", va="top")

    for idx_r, (f, out_t, let) in enumerate([("B_ASC_IgG", "sev", "j"), ("B_ASC_IgA", "sev", "k"), ("HSPC", "sev", "l"), ("T CD8 CTL", "sev", "m"), ("T CD8 CTL", "death", "n"), ("HSPC", "death", "o")]):
        x_l = 0.052 + idx_r * 0.158
        ax_rw = fig_s2.add_axes([x_l, 0.050, 0.125, 0.175]); draw_strip_raw_pct(ax_rw, f, out_t, rng_s2, show_ylabel=(idx_r == 0))
        ax_rw.set_title(f"{LBL[f]} ({'s' if out_t=='sev' else 'd'})", loc="left", fontsize=10.4, fontweight="bold", color=COL["HSPC"] if f == "HSPC" else COL["Real"], pad=5); fig_s2.text(x_l - 0.030, 0.242, let, fontsize=17, fontweight="bold", va="top")

    saved_s2 = finalize_figure(fig_s2, out_dir / "PBISC_SuppFig_S2_CLL_Anchors_Modality_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Supplementary Figure 2:", saved_s2[0].name)

    # =====================================================================
    # 5) SUPPLEMENTARY FIGURE 3 (a - c) — RE-DETECTION, CD14 IFN, 14-VARIANT MATRIX (TRUE_PUB_11)
    # =====================================================================
    def draw_s3a(ax):
        for fk_p, ep_p, lbl_p, cc_p, ls_p, has_real in [
            ("B_ASC_IgG", "sev", "B ASC IgG (sev)", COL["Ext"], "-", True),
            ("B_ASC_IgA", "sev", "B ASC IgA (sev)", COL["New"], "-", True),
            ("HSPC", "sev", "HSPC (sev)", PALETTE["teal"], "-", True),
            ("T CD8 CTL", "sev", "T CD8 CTL (sev)", PALETTE["violet"], "-", False),
            ("T CD8 CTL", "death", "T CD8 CTL (death)", PALETTE["gray_3"], "-", False),
            ("HSPC", "death", "HSPC (death)", PALETTE["red_strong"], "-", True),
        ]:
            d = PW_UNI[(PW_UNI["feature"] == fk_p) & (PW_UNI["outcome"] == ep_p)].sort_values("n_real")
            if has_real:
                ax.plot(d["n_real"], d["real_nom"], color=cc_p, lw=1.5, ls="--", marker="D", ms=4.0, mfc="white", label=f"{lbl_p}: Real")
            ax.plot(d["n_real"], d["plus_new29_nom"], color=cc_p, lw=2.1, ls="-", marker="s", ms=4.8, label=f"{lbl_p}: +New29")
        d0 = PW_UNI[(PW_UNI["feature"] == "HSPC") & (PW_UNI["outcome"] == "death")].sort_values("n_real")
        shuf_avg = PW_UNI.groupby("n_real")["plus_shuf_nom"].mean().loc[d0["n_real"]].values
        ax.plot(d0["n_real"], shuf_avg, color=PALETTE["gray_1"], lw=1.8, ls=":", label="+New29 shuffled (mean)")
        pw20 = PW_UNI[PW_UNI["n_real"] == 20]
        for idx_b, (fk_b, out_b, cc_b) in enumerate([("B_ASC_IgG", "sev", COL["Ext"]), ("HSPC", "sev", PALETTE["teal"]), ("T CD8 CTL", "death", PALETTE["gray_3"]), ("HSPC", "death", PALETTE["red_strong"])]):
            bh_v = float(pw20[(pw20["feature"] == fk_b) & (pw20["outcome"] == out_b)]["plus_new29_bh43"].iloc[0])
            ax.scatter(20.0, bh_v, s=72, marker="P", c="white", edgecolors=cc_b, linewidths=1.8, zorder=6, label="Same-run BH43 at n=20 (+New29)" if idx_b == 0 else "")
        ax.axhline(0.8, color=PALETTE["gray_1"], ls=":", lw=1.3)
        ax.set_ylim(-0.03, 1.45); ax.set_yticks([0, 0.4, 0.8, 1.0]); ax.set_xticks([10, 15, 20, 25, 30, 35]); ax.set_xlim(8.5, 36.5)
        ax.set_xlabel("Sampled Real patients (1,000 unified draws; fixed +PBISC new 29)", fontsize=10.0)
        ax.set_ylabel("Conditional re-detection rate (P < 0.05)", fontsize=10.0)
        ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=2, fontsize=7.6, handlelength=1.4, columnspacing=0.8, handletextpad=0.30, labelspacing=0.18)

    def draw_s3c(ax_sev, ax_death):
        FEN = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN"]
        VAR_ORD = SA["variant"].drop_duplicates().tolist()
        for ax, val, pcol, ttl, nrm, show_y in [
            (ax_sev, "rho", "p_sev", "Severity: ρ (Real + PBISC new, n = 66)", mpl.colors.Normalize(-0.7, 0.7), True),
            (ax_death, "AUC", "p_death", "Death: AUC (Real + PBISC new; TRUE 11d baseline)", mpl.colors.Normalize(0.15, 0.85), False),
        ]:
            M = SA.pivot_table(index="variant", columns="feature", values=val).reindex(index=VAR_ORD, columns=FEN)
            Pm = SA.pivot_table(index="variant", columns="feature", values=pcol).reindex(index=VAR_ORD, columns=FEN)
            ax.imshow(M.values, cmap="RdBu_r", norm=nrm, aspect="auto")
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    v = M.values[i, j]
                    if np.isnan(v): ax.text(j, i, "–", ha="center", va="center", fontsize=9.0, color=PALETTE["gray_1"])
                    else:
                        sig = Pm.values[i, j] < 0.05
                        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.8, fontweight="bold" if sig else "normal",
                                color=("white" if abs(nrm(v) - 0.5) > 0.32 else "black") if sig else PALETTE["gray_1"])
            ax.axvline(3.5, color="black", lw=1.6)
            ax.set_xticks(range(len(FEN))); ax.set_xticklabels([LBL[f] for f in FEN], rotation=32, ha="left", fontsize=9.6); ax.xaxis.tick_top()
            ax.set_yticks(range(len(VAR_ORD))); ax.set_yticklabels(VAR_ORD if show_y else [], fontsize=9.2); ax.tick_params(length=0)
            for j, t in enumerate(ax.get_xticklabels()):
                f = FEN[j]; t.set_fontweight("bold")
                t.set_color(COL["Real"] if j < 4 else (PALETTE["red_strong"] if f in {"HSPC", "CD14mono_IFN"} else COL["Ext"]))
            for sp_c in ax.spines.values(): sp_c.set_visible(False)
            ax.set_xticks(np.arange(-0.5, len(FEN)), minor=True); ax.set_yticks(np.arange(-0.5, len(VAR_ORD)), minor=True)
            ax.grid(which="minor", color="white", lw=1.3); ax.tick_params(which="minor", length=0)
            ax.set_title(ttl, loc="left", fontsize=11.0, fontweight="bold", pad=28)

    fig_s3 = plt.figure(figsize=(15.4, 13.6), facecolor="white")
    fig_s3.text(0.015, 0.986, "Supplementary Figure 3. Unified conditional re-detection curves, CD14 mono IFN check, and 14-variant sensitivity matrix (TRUE 11 deaths)",
                fontsize=13.2, fontweight="bold", va="top", ha="left", color=COL["Real"])
    ax_s3a = fig_s3.add_axes([0.075, 0.595, 0.545, 0.320]); draw_s3a(ax_s3a); ax_s3a.set_title("Conditional re-detection curves across Real sample sizes (n = 10–35)", loc="left", fontsize=11.0, fontweight="bold", pad=8); fig_s3.text(0.015, 0.942, "a", fontsize=17.5, fontweight="bold", va="top")
    ax_s3b = fig_s3.add_axes([0.725, 0.595, 0.240, 0.320]); draw_strip_ns(ax_s3b, "CD14mono_IFN", "death", np.random.default_rng(23), show_ylabel=True)
    ax_s3b.text(0.03, 0.67, "Did not pass BH43 (q43 = 0.184)\nunder verified 11 deaths", transform=ax_s3b.transAxes, fontsize=8.4, color=PALETTE["red_strong"], fontweight="bold", va="top")
    ax_s3b.set_title("CD14 mono IFN (death)", loc="left", fontsize=11.2, fontweight="bold", color=COL["Real"], pad=8); fig_s3.text(0.658, 0.942, "b", fontsize=17.5, fontweight="bold", va="top")
    ax_s3c_s = fig_s3.add_axes([0.238, 0.040, 0.360, 0.385]); ax_s3c_d = fig_s3.add_axes([0.616, 0.040, 0.360, 0.385]); draw_s3c(ax_s3c_s, ax_s3c_d); fig_s3.text(0.015, 0.505, "c", fontsize=17.5, fontweight="bold", va="top")

    saved_s3 = finalize_figure(fig_s3, out_dir / "PBISC_SuppFig_S3_Power_IFN_Sensitivity_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Supplementary Figure 3:", saved_s3[0].name)

    # =====================================================================
    # 6) SUPPLEMENTARY FIGURE 4 (a - d) — 4-COHORT EXTERNAL REPLICATION
    # =====================================================================
    fig_s4 = plt.figure(figsize=(16.2, 14.8), facecolor="white")
    fig_s4.text(0.015, 0.988, "Supplementary Figure 4. External single-cell replication across four adult COVID-19 PBMC cohorts (300 acute samples / 231 donors; N_sev = 220, N_death = 219 [24 deaths])",
                fontsize=12.3, fontweight="bold", va="top", ha="left", color=PALETTE["gray_3"])

    ax_4a = fig_s4.add_axes([0.148, 0.535, 0.250, 0.375])
    ORDER_4A = [
        ("NK_CD56bright", "NK CD56bright", "NK CD56 bright", "anchor"),
        ("B_ASC_dividing", "B ASC dividing*", "B_ASC_dividing", "anchor"),
        ("T_CD8_CM", "T CD8 CM / EM*", "T CD8 CM", "anchor"),
        ("CD16_mono_ofMono", "CD16 mono (% Mono)", "CD16mono|lin", "anchor"),
        ("B_ASC_total", "B ASC total", "B_ASC_total", "focus"),
        ("T_CD8_CTL", "T CD8 CTL / Mem*", "T CD8 CTL", "focus"),
        ("HSPC", "HSPC", "HSPC", "hspc"),
        ("NK_CD56dim", "NK CD56dim", "NK", "excl"),
        ("CD16_mono_C1", "CD16 mono C1", "CD16mono_C1", "excl"),
    ]
    xlim_4a = (-0.82, 0.80); y_pos_4, yc_4 = [], 0.0
    for idx, (pub_f, disp, eskd_f, grp) in enumerate(ORDER_4A):
        if idx in (4, 7): ax_4a.axhline(yc_4 + 0.52, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=1); yc_4 -= 0.35
        y_pos_4.append(yc_4); em = ESKD_MATCH_DF.loc[eskd_f]
        mr = PUB_META_CORR[(PUB_META_CORR["feature"] == pub_f) & (PUB_META_CORR["endpoint"] == "severity")].iloc[0]
        sub_c = PUB_COH_CORR[PUB_COH_CORR["feature"] == pub_f].dropna(subset=["rho"])
        ax_4a.plot([em["lo"], em["hi"]], [yc_4 + 0.18, yc_4 + 0.18], color=PALETTE["blue_main"], lw=1.7, zorder=2)
        ax_4a.scatter(em["rho"], yc_4 + 0.18, s=46, marker="s", c=PALETTE["blue_main"], edgecolors=PALETTE["blue_main"], zorder=3)
        p_col = PALETTE["red_strong"] if grp == "hspc" else (PALETTE["gray_1"] if grp == "excl" else PALETTE["teal"])
        ax_4a.scatter(sub_c["rho"], np.repeat(yc_4 - 0.18, len(sub_c)), s=24, marker="o", c=PALETTE["neutral"], edgecolors=PALETTE["gray_1"], linewidths=0.8, zorder=2)
        plot_dl_and_mkh(ax_4a, yc_4 - 0.18, float(mr["DL_lo"]), float(mr["DL_hi"]), float(mr["mHK_lo"]), float(mr["mHK_hi"]), int(mr["k"]), p_col, xlim_4a)
        ax_4a.scatter(mr["estimate"], yc_4 - 0.18, s=52, marker="D", c="white" if mr["DL_p"] >= 0.05 else p_col, edgecolors=p_col, linewidths=1.6, zorder=4)
        k_i, n_i = int(mr["k"]), int(mr["N"]); mkh_s = f", mKH {fp_short(mr['mHK_p'])}" if k_i >= 2 else ""
        ax_4a.text(0.83, yc_4 - 0.02, f"{mr['estimate']:+.2f} (DL {fp_short(mr['DL_p'])}{mkh_s}) · k={k_i}, N={n_i}",
                   va="center", ha="left", fontsize=7.4, color=p_col if mr["DL_p"] < 0.05 else PALETTE["gray_1"], fontweight="bold" if mr["DL_p"] < 0.05 else "normal")
        yc_4 -= 1.05
    ax_4a.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
    ax_4a.set_yticks(y_pos_4); ax_4a.set_yticklabels([x[1] for x in ORDER_4A], fontsize=9.4)
    for idx_t, t in enumerate(ax_4a.get_yticklabels()):
        grp = ORDER_4A[idx_t][3]; t.set_fontweight("bold" if grp != "excl" else "normal")
        t.set_color(PALETTE["gray_3"] if grp == "anchor" else (PALETTE["red_strong"] if grp == "hspc" else (PALETTE["gray_1"] if grp == "excl" else PALETTE["blue_main"])))
    ax_4a.tick_params(axis="y", length=0); ax_4a.spines["left"].set_visible(False)
    ax_4a.set_xlim(*xlim_4a); ax_4a.set_xticks([-0.5, 0.0, 0.5]); ax_4a.set_ylim(y_pos_4[-1] - 0.55, 0.55)
    ax_4a.set_xlabel("Spearman ρ with WHO severity", fontsize=10.0)
    ax_4a.text(0.83, 0.50, "Public Meta ρ (DL; mKH) · k, N", fontsize=7.9, fontweight="bold", color=PALETTE["gray_2"], va="bottom")
    ax_4a.set_title("Matched severity comparison", loc="left", fontsize=11.2, fontweight="bold", pad=8)
    fig_s4.text(0.015, 0.930, "a", fontsize=18, fontweight="bold", va="top")

    ax_4b1 = fig_s4.add_axes([0.676, 0.575, 0.128, 0.335]); ax_4b2 = fig_s4.add_axes([0.832, 0.575, 0.138, 0.335])
    ISO_ORDER = [("B_ASC_dividing", "Plasmablast (dividing)", "B_ASC_dividing"), ("B_ASC_IgG", "Plasma cell IgG", "B_ASC_IgG"), ("B_ASC_IgA", "Plasma cell IgA", "B_ASC_IgA"), ("B_ASC_IgM", "Plasma cell IgM", "B_ASC_IgM"), ("B_ASC_total", "B ASC total", "B_ASC_total")]
    st_b_df = PUB_COH_CORR[PUB_COH_CORR["cohort"] == "Stephenson 2021"].set_index("feature")
    yy_4b = np.arange(len(ISO_ORDER))[::-1]
    ax_4b1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3); ax_4b2.axvline(1.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for yi, (fk, flab, efk) in zip(yy_4b, ISO_ORDER):
        r, em = st_b_df.loc[fk], ESKD_MATCH_DF.loc[efk]
        col_b = PALETTE["blue_main"] if fk in ("B_ASC_IgG", "B_ASC_IgA") else PALETTE["teal"]
        ax_4b1.plot([em["lo"], em["hi"]], [yi + 0.16, yi + 0.16], color=PALETTE["blue_main"], lw=1.6)
        ax_4b1.scatter(em["rho"], yi + 0.16, s=42, marker="s", c=PALETTE["blue_main"], edgecolors=PALETTE["blue_main"], zorder=3)
        ax_4b1.plot([r["rho_lo"], r["rho_hi"]], [yi - 0.16, yi - 0.16], color=PALETTE["teal"], lw=1.8)
        ax_4b1.scatter(r["rho"], yi - 0.16, s=48, marker="D", c="white" if r["p_sev"] >= 0.05 else PALETTE["teal"], edgecolors=PALETTE["teal"], linewidths=1.5, zorder=4)
        ax_4b1.text(0.99, yi - 0.16, f"P={r['p_sev']:.2f}" if r["p_sev"] >= 0.1 else f"P={r['p_sev']:.3f}", va="center", ha="right", fontsize=7.5, color=PALETTE["gray_1"] if r["p_sev"] >= 0.05 else PALETTE["teal"], fontweight="bold" if r["p_sev"] < 0.05 else "normal")
        ax_4b2.plot([r["OR_lo"], r["OR_hi"]], [yi, yi], color=col_b, lw=2.0)
        ax_4b2.scatter(r["OR_glm"], yi, s=50, marker="D", c=col_b, edgecolors="white", linewidths=0.7, zorder=4)
        ax_4b2.text(r["OR_hi"] + 0.08, yi, f"{r['OR_glm']:.2f} (P={r['p_glm']:.3f})", va="center", ha="left", fontsize=7.6, color=col_b, fontweight="bold")
    ax_4b1.set_yticks(yy_4b); ax_4b1.set_yticklabels([x[1] for x in ISO_ORDER], fontsize=9.2)
    for idx_t, t in enumerate(ax_4b1.get_yticklabels()):
        if ISO_ORDER[idx_t][0] in ("B_ASC_IgG", "B_ASC_IgA"): t.set_fontweight("bold"); t.set_color(PALETTE["blue_main"])
    ax_4b2.set_yticks(yy_4b); ax_4b2.set_yticklabels([])
    for axx in (ax_4b1, ax_4b2): axx.tick_params(axis="y", length=0); axx.spines["left"].set_visible(False); axx.set_ylim(-0.55, len(ISO_ORDER) - 0.35)
    ax_4b1.set_xlim(-0.20, 1.04); ax_4b1.set_xticks([0.0, 0.4, 0.8]); ax_4b1.set_xlabel("Spearman ρ (% PBMC)", fontsize=9.6)
    ax_4b1.set_title("Stephenson 2021 B ASC (n = 74): rank ρ vs. quasi-binomial GLM", loc="left", fontsize=10.4, fontweight="bold", pad=8)
    ax_4b2.set_xlim(0.85, 5.10); ax_4b2.set_xticks([1.0, 2.0, 3.0, 4.0]); ax_4b2.set_xlabel("Exploratory logit-GLM OR / step", fontsize=9.4)
    fig_s4.text(0.552, 0.930, "b", fontsize=18, fontweight="bold", va="top")

    ax_4leg = fig_s4.add_axes([0.015, 0.462, 0.97, 0.026]); ax_4leg.set_axis_off()
    ax_4leg.legend(handles=[
        Line2D([0], [0], marker="s", ls="", ms=7.2, mfc=PALETTE["blue_main"], mec=PALETTE["blue_main"], label="ESKD Real + PBISC new (n = 66, 11 verified deaths)"),
        Line2D([0], [0], marker="D", ls="-", lw=2.0, ms=7.2, mfc=PALETTE["teal"], mec=PALETTE["teal"], color=PALETTE["teal"], label="Public Meta DerSimonian–Laird 95% CI (solid)"),
        Line2D([0], [0], ls="--", lw=1.3, color=PALETTE["teal"], label="modified Knapp–Hartung 95% CI (dashed; arrows = clipped)"),
        Line2D([0], [0], marker="o", ls="", ms=6.0, mfc=PALETTE["neutral"], mec=PALETTE["gray_1"], label="Individual cohort (a)"),
        Line2D([0], [0], marker="D", ls="", ms=7.2, mfc="white", mec=PALETTE["gray_2"], mew=1.5, label="open: DL P ≥ 0.05"),
    ], loc="center", ncol=5, fontsize=8.6, frameon=False, columnspacing=0.95)

    ax_4c = fig_s4.add_axes([0.155, 0.098, 0.240, 0.318])
    D_ITEMS_4 = [
        ("HSPC", "HSPC (k=3)", PALETTE["red_strong"]), ("T_CD8_CM", "T CD8 CM / EM* (k=2)", PALETTE["blue_main"]),
        ("T_CD8_CTL", "T CD8 CTL / Mem* (k=4)", PALETTE["blue_main"]), ("B_ASC_total", "B ASC total (k=4)", PALETTE["gray_2"]),
        ("NK_CD56dim", "NK CD56dim (excl., k=4)", PALETTE["gray_1"]), ("CD16_mono_C1", "CD16 mono C1 (excl., k=1)", PALETTE["gray_1"]),
    ]
    xlim_4c = (0.08, 0.96); yy_4c = np.arange(len(D_ITEMS_4))[::-1]; ax_4c.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for yi, (pub_f, disp, col_c) in zip(yy_4c, D_ITEMS_4):
        mr_all = PUB_META_CORR[(PUB_META_CORR["feature"] == pub_f) & (PUB_META_CORR["endpoint"] == "death_all")].iloc[0]
        mr_w6 = PUB_META_CORR[(PUB_META_CORR["feature"] == pub_f) & (PUB_META_CORR["endpoint"] == "death_WHO6")].iloc[0]
        k_c = int(mr_all["k"])
        plot_dl_and_mkh(ax_4c, yi + 0.17, float(mr_all["DL_lo"]), float(mr_all["DL_hi"]), float(mr_all["mHK_lo"]), float(mr_all["mHK_hi"]), k_c, col_c, xlim_4c, lw_dl=1.9)
        ax_4c.scatter(mr_all["estimate"], yi + 0.17, s=48, marker="D", c="white" if mr_all["DL_p"] >= 0.05 else col_c, edgecolors=col_c, linewidths=1.5, zorder=4)
        plot_ci_clipped(ax_4c, yi - 0.17, float(mr_w6["DL_lo"]), float(mr_w6["DL_hi"]), col_c, xlim_4c, lw=1.4, ls=":", alpha=0.82, zorder=3)
        ax_4c.scatter(mr_w6["estimate"], yi - 0.17, s=44, marker="^", c="white" if mr_w6["DL_p"] >= 0.05 else col_c, edgecolors=col_c, linewidths=1.3, zorder=4)
        mkh_str = f", mKH {fp_short(mr_all['mHK_p'])}" if k_c >= 2 else ""
        na, da, nw, dw = int(mr_all["N"]), int(mr_all["deaths"]), int(mr_w6["N"]), int(mr_w6["deaths"])
        ax_4c.text(0.98, yi + 0.17, f"All ({na}/{da}d): {mr_all['estimate']:.2f} (DL {fp_short(mr_all['DL_p'])}{mkh_str})", va="center", fontsize=7.25, color=col_c if mr_all["DL_p"] < 0.05 else PALETTE["gray_2"], fontweight="bold" if mr_all["DL_p"] < 0.05 else "normal")
        ax_4c.text(0.98, yi - 0.17, f"WHO≥6 ({nw}/{dw}d): {mr_w6['estimate']:.2f} (DL {fp_short(mr_w6['DL_p'])})", va="center", fontsize=7.15, color=PALETTE["gray_1"])
    ax_4c.set_yticks(yy_4c); ax_4c.set_yticklabels([x[1] for x in D_ITEMS_4], fontsize=8.9)
    for idx_t in range(3): ax_4c.get_yticklabels()[idx_t].set_fontweight("bold"); ax_4c.get_yticklabels()[idx_t].set_color(D_ITEMS_4[idx_t][2])
    ax_4c.tick_params(axis="y", length=0); ax_4c.spines["left"].set_visible(False)
    ax_4c.set_xlim(*xlim_4c); ax_4c.set_xticks([0.25, 0.50, 0.75]); ax_4c.set_ylim(-0.55, len(D_ITEMS_4) - 0.35)
    ax_4c.set_xlabel("Public Meta Death AUC (diamond: All · tri.: WHO≥6)", fontsize=9.3)
    ax_4c.set_title("Mortality: all-patient vs. within WHO-CPS ≥ 6", loc="left", fontsize=11.0, fontweight="bold", pad=8)
    fig_s4.text(0.015, 0.435, "c", fontsize=18, fontweight="bold", va="top")

    ax_4d = fig_s4.add_axes([0.675, 0.098, 0.175, 0.318])
    COH_4D = [("Stephenson 2021", "Stephenson 2021"), ("van der Wijst 2021", "van der Wijst 2021"), ("COMBAT 2022", "COMBAT 2022"), ("Meta3", "Public Meta (k = 3)")]
    xlim_4d = (0.34, 0.99); yy_4d = np.arange(len(COH_4D))[::-1]; ax_4d.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for yi, (cn, clab) in zip(yy_4d, COH_4D):
        r1 = PUB_HW_CORR[(PUB_HW_CORR["cohort"] == cn) & (PUB_HW_CORR["window"] == "First acute sample") & (PUB_HW_CORR["endpoint"] == "all")].iloc[0]
        r2 = PUB_HW_CORR[(PUB_HW_CORR["cohort"] == cn) & (PUB_HW_CORR["window"] == "Samples > day 7") & (PUB_HW_CORR["endpoint"] == "all")].iloc[0]
        is_meta = cn == "Meta3"
        a1, p1, n1, d1 = float(r1["estimate"] if is_meta else r1["AUC"]), float(r1["DL_p"] if is_meta else r1["p"]), int(r1["n"]), int(r1["deaths"])
        a2, p2, n2, d2 = float(r2["estimate"] if is_meta else r2["AUC"]), float(r2["DL_p"] if is_meta else r2["p"]), int(r2["n"]), int(r2["deaths"])
        if is_meta:
            plot_dl_and_mkh(ax_4d, yi + 0.16, float(r1["DL_lo"]), float(r1["DL_hi"]), float(r1["mHK_lo"]), float(r1["mHK_hi"]), 3, PALETTE["gray_3"], xlim_4d, lw_dl=1.8)
            plot_dl_and_mkh(ax_4d, yi - 0.16, float(r2["DL_lo"]), float(r2["DL_hi"]), float(r2["mHK_lo"]), float(r2["mHK_hi"]), 3, PALETTE["red_strong"], xlim_4d, lw_dl=1.9)
        else:
            ax_4d.plot([a1, a2], [yi + 0.16, yi - 0.16], color=PALETTE["neutral"], lw=1.5, zorder=2)
        ax_4d.scatter(a1, yi + 0.16, s=60 if is_meta else 46, marker="D", c="white" if p1 >= 0.05 else PALETTE["gray_3"], edgecolors=PALETTE["gray_3"], linewidths=1.5, zorder=3)
        ax_4d.scatter(a2, yi - 0.16, s=64 if is_meta else 50, marker="s", c="white" if p2 >= 0.05 else PALETTE["red_strong"], edgecolors=PALETTE["red_strong"], linewidths=1.5, zorder=4)
        t1_s = f"1st ({n1}/{d1}d): {a1:.2f} (DL {p1:.3f}, mKH {r1['mHK_p']:.2f})" if is_meta else f"1st ({n1}/{d1}d): {a1:.2f} (P={p1:.3f})"
        t2_s = f">7d ({n2}/{d2}d): {a2:.2f} (DL {p2:.3f}, mKH {r2['mHK_p']:.2f})" if is_meta else f">7d ({n2}/{d2}d): {a2:.2f} (P={p2:.3f})"
        ax_4d.text(1.01, yi + 0.16, t1_s, va="center", fontsize=7.25, color=PALETTE["gray_3"] if p1 < 0.05 else PALETTE["gray_1"], fontweight="bold" if p1 < 0.05 else "normal")
        ax_4d.text(1.01, yi - 0.16, t2_s, va="center", fontsize=7.25, color=PALETTE["red_strong"] if p2 < 0.05 else PALETTE["gray_1"], fontweight="bold" if p2 < 0.05 else "normal")
    ax_4d.set_yticks(yy_4d); ax_4d.set_yticklabels([c[1] for c in COH_4D], fontsize=9.1)
    ax_4d.get_yticklabels()[-1].set_fontweight("bold"); ax_4d.get_yticklabels()[-1].set_color(PALETTE["red_strong"])
    ax_4d.tick_params(axis="y", length=0); ax_4d.spines["left"].set_visible(False)
    ax_4d.set_xlim(*xlim_4d); ax_4d.set_xticks([0.50, 0.70, 0.90]); ax_4d.set_ylim(-0.55, len(COH_4D) - 0.35)
    ax_4d.set_xlabel("HSPC Death AUC (diamond: 1st · sq.: >7d)", fontsize=9.3)
    ax_4d.set_title("HSPC mortality across cross-sectional windows", loc="left", fontsize=10.9, fontweight="bold", pad=8)
    fig_s4.text(0.552, 0.435, "d", fontsize=18, fontweight="bold", va="top")

    fig_s4.text(0.015, 0.014,
                "*Stephenson 2021 pre-filter COVID+ samples (101) reconciled to published Supp Table 2; acute (<28 d, ≥500 cells) scope comprises 89 samples across 85 donors: "
                "74 WHO-known, 73 outcome-known (7 deaths; 12 outcome-unknown excluded from death), and 63 with both known (AP3 deceased lacks WHO). Liu 2021 lacks HSPC (k=3). Dashed = mKH 95% CIs.",
                fontsize=7.65, color=PALETTE["gray_2"], va="bottom", ha="left")

    saved_s4 = finalize_figure(fig_s4, out_dir / "PBISC_SuppFig_S4_Public_COVID_Replication_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Supplementary Figure 4:", saved_s4[0].name)

    # =====================================================================
    # 7) SUPPLEMENTARY FIGURE 5 (a - h) — UNBIASED 9-METHOD BENCHMARK, HSPC 20-GENE, T-CELL ISG 13-SUBSETS, LABS
    # =====================================================================
    fig_s5 = plt.figure(figsize=(17.6, 20.4), facecolor="white")
    fig_s5.text(0.014, 0.991, "Supplementary Figure 5. Unbiased 9-method deconvolution summary, 20-gene HSPC matrix, 13-subset T-cell ISG dynamics, and hospital labs",
                fontsize=13.0, fontweight="bold", va="top", ha="left", color=COL["Real"])

    # Panel 5a: Unbiased 9-method summary heatmap (All-43 & Rare-23 Shared-27 fidelity + T1 concordance + 10-signal score)
    ax_5a = fig_s5.add_axes([0.128, 0.778, 0.425, 0.145]); fig_s5.text(0.014, 0.955, "a", fontsize=17.5, fontweight="bold", va="top")
    UNB_COLS = [
        ("med_rho_all43", "Sh27 med ρ\n(all 43)"), ("med_rho_rare23", "Sh27 med ρ\n(rare 23)"),
        ("T1_conc_all43", "T1 conc.\n(all 43)"), ("T1_conc_rare23", "T1 conc.\n(rare 23)"),
        ("pass_10_signals", "Clinical\n(/10)"),
    ]
    M_u = UNB9[[c[0] for c in UNB_COLS]].values.astype(float).copy()
    M_u_norm = M_u.copy(); M_u_norm[:, 4] = M_u_norm[:, 4] / 10.0
    ax_5a.imshow(M_u_norm, cmap="YlGnBu", norm=mpl.colors.Normalize(0.0, 0.95), aspect="auto")
    for i_u in range(M_u.shape[0]):
        for j_u in range(M_u.shape[1]):
            v = M_u[i_u, j_u]
            txt = f"{int(v)}/10" if j_u == 4 else f"{v:+.2f}"
            ax_5a.text(j_u, i_u, txt, ha="center", va="center", fontsize=8.4, fontweight="bold" if (i_u == 0 or (j_u == 4 and v >= 6) or (j_u < 4 and v >= 0.30)) else "normal",
                       color="white" if M_u_norm[i_u, j_u] > 0.52 else COL["Real"])
    ax_5a.axvline(3.5, color="black", lw=1.6)
    ax_5a.set_yticks(range(len(UNB9))); ax_5a.set_yticklabels(UNB9["method"].tolist(), fontsize=8.6)
    ax_5a.get_yticklabels()[0].set_fontweight("bold"); ax_5a.get_yticklabels()[0].set_color(COL["Ext"])
    ax_5a.set_xticks(range(len(UNB_COLS))); ax_5a.set_xticklabels([c[1] for c in UNB_COLS], fontsize=8.4, fontweight="bold")
    ax_5a.xaxis.tick_top(); ax_5a.tick_params(length=0)
    for sp_5a in ax_5a.spines.values(): sp_5a.set_visible(False)
    ax_5a.set_title("Unbiased 9-method comparison: generic 43-type ρ vs. rare clinical signals", loc="left", fontsize=10.0, fontweight="bold", pad=22)

    # Panel 5b: Single-gene raw bulk HSPC vs PBISC (Wave-1 new n=29, unclipped [0, 1.02])
    ax_5b = fig_s5.add_axes([0.708, 0.778, 0.195, 0.145]); fig_s5.text(0.585, 0.955, "b", fontsize=17.5, fontweight="bold", va="top")
    yy_5b = np.arange(len(SG_HSPC))[::-1]; ax_5b.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for yi, (_, r_sg) in zip(yy_5b, SG_HSPC.iterrows()):
        is_pb = bool(r_sg["is_pbisc"]); cc = COL["HSPC"] if is_pb else PALETTE["gray_2"]
        ax_5b.plot([r_sg["lo"], r_sg["hi"]], [yi, yi], color=cc, lw=1.9 if is_pb else 1.2)
        ax_5b.scatter(r_sg["auc"], yi, s=48 if is_pb else 36, marker="s" if is_pb else "o", c=cc if r_sg["p"] < 0.05 else "white", edgecolors=cc, linewidths=1.5, zorder=4)
        ax_5b.text(1.03, yi, f"{r_sg['auc']:.2f}* (P={r_sg['p']:.3f})" if r_sg["p"] < 0.05 else f"{r_sg['auc']:.2f} (ns)", va="center", fontsize=7.6, fontweight="bold" if is_pb else "normal", color=cc if r_sg["p"] < 0.05 else PALETTE["gray_1"])
    ax_5b.set_yticks(yy_5b); ax_5b.set_yticklabels(SG_HSPC["label"].tolist(), fontsize=8.6)
    ax_5b.get_yticklabels()[0].set_fontweight("bold"); ax_5b.get_yticklabels()[0].set_color(COL["HSPC"])
    ax_5b.tick_params(axis="y", length=0); ax_5b.spines["left"].set_visible(False)
    ax_5b.set_xlim(0.00, 1.02); ax_5b.set_xticks([0.0, 0.5, 1.0]); ax_5b.set_xlabel("Death AUC (95% CI, n=29, 6d)", fontsize=8.8)
    ax_5b.set_title("Single-gene raw bulk vs. PBISC (unclipped)", loc="left", fontsize=10.3, fontweight="bold", pad=10)

    # Row 2: Panel 5c (20-gene HSPC Severity & TRUE-11 Death matrix!) & Panel 5d (All 13 T-cell compartments Early vs Late + GEE Time x Fatal!)
    ax_5c = fig_s5.add_axes([0.115, 0.552, 0.500, 0.148]); fig_s5.text(0.014, 0.728, "c", fontsize=17.5, fontweight="bold", va="top")
    GENES_20 = ["GATA1", "KLF1", "TFRC", "BLVRB", "APOC1", "GATA2", "ITGA2B", "MKI67", "TOP2A", "STMN1", "TYMS", "CD34", "AVP", "HLF", "HOPX", "CRHBP", "MPO", "ELANE", "CSF3R", "CEBPD"]
    cmap_g, nrm_g = mpl.colormaps["RdBu_r"], mpl.colors.Normalize(-0.7, 0.7)
    G20_ROWS = [
        ("Real", "rho_sev", "p_sev", "q20_sev", "Sev · Real (33)", 0.0),
        ("PBISC new", "rho_sev", "p_sev", "q20_sev", "Sev · New (29)", -1.0),
        ("Combined", "rho_sev", "p_sev", "q20_sev", "Sev · Comb (62)", -2.0),
        ("Real", "AUC_dth11", "p_dth11", "q20_dth11", "Dth · Real (5d)", -3.35),
        ("PBISC new", "AUC_dth11", "p_dth11", "q20_dth11", "Dth · New (6d)", -4.35),
        ("Combined", "AUC_dth11", "p_dth11", "q20_dth11", "Dth · Comb (11d)", -5.35),
    ]
    ax_5c.axhline(-2.675, color=PALETTE["neutral"], ls="--", lw=1.0)
    for c_key, vcol, pcol, qcol, rlbl, yv in G20_ROWS:
        sub_g = HSPC_G20[HSPC_G20["cohort"] == c_key].set_index("gene")
        for g_i, gn in enumerate(GENES_20):
            rg = sub_g.loc[gn]
            eff = float(rg[vcol]) if vcol == "rho_sev" else 2.0 * (float(rg[vcol]) - 0.5)
            pv, qv = float(rg[pcol]), float(rg[qcol])
            sz = 26 + 110 * min(-np.log10(max(pv, 1e-4)), 3.0) / 3.0
            ec = "black" if qv < 0.05 else (PALETTE["gray_1"] if pv < 0.05 else "none")
            ax_5c.scatter(g_i, yv, s=sz, c=[cmap_g(nrm_g(eff))], edgecolors=ec, linewidths=1.7 if qv < 0.05 else (0.9 if pv < 0.05 else 0.0), zorder=3)
    for div_x in [6.5, 10.5, 15.5]: ax_5c.axvline(div_x, color=PALETTE["neutral"], ls="--", lw=1.0)
    ax_5c.text(3.0, 0.72, "Erythroid / MEP", ha="center", fontsize=8.0, fontweight="bold", color=COL["HSPC"])
    ax_5c.text(8.5, 0.72, "Prolif.", ha="center", fontsize=8.0, fontweight="bold", color=COL["HSPC"])
    ax_5c.text(13.0, 0.72, "Quiescent HSC", ha="center", fontsize=8.0, fontweight="bold", color=PALETTE["gray_2"])
    ax_5c.text(17.5, 0.72, "Granulocyte", ha="center", fontsize=8.0, fontweight="bold", color=PALETTE["gray_2"])
    ax_5c.set_yticks([r[5] for r in G20_ROWS]); ax_5c.set_yticklabels([r[4] for r in G20_ROWS], fontsize=8.0, fontweight="bold")
    ax_5c.set_xticks(range(len(GENES_20))); ax_5c.set_xticklabels(GENES_20, rotation=44, ha="right", fontsize=8.0, fontweight="bold")
    for idx_g, t in enumerate(ax_5c.get_xticklabels()): t.set_color(COL["HSPC"] if idx_g < 11 else PALETTE["gray_2"])
    ax_5c.tick_params(length=0)
    for sp_5c in ax_5c.spines.values(): sp_5c.set_visible(False)
    ax_5c.set_xlim(-0.6, len(GENES_20) - 0.4); ax_5c.set_ylim(-5.95, 1.05)
    ax_5c.set_title("20-gene HSPC Severity & Death matrix (bold ring: BH20 q<0.05; thin: P<0.05)", loc="left", fontsize=9.8, fontweight="bold", pad=12)

    # Panel 5d: T-cell ISG Early (0-7d) vs Late (>7d) Severity + Source-specific GEE Time x Fatal collapse!
    ax_5d1 = fig_s5.add_axes([0.708, 0.552, 0.115, 0.148]); ax_5d2 = fig_s5.add_axes([0.858, 0.552, 0.120, 0.148])
    fig_s5.text(0.642, 0.728, "d", fontsize=17.5, fontweight="bold", va="top")
    TS5 = ["T CD8 CTL", "T CD4 CTL", "All CD8 T", "All CD4 T", "All T cells"]
    yy_5d = np.arange(len(TS5))[::-1]
    ax_5d1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.2)
    for yi, ts in zip(yy_5d, TS5):
        rt = t_ix.loc[ts]
        ax_5d1.plot([rt["late_sev"], rt["early_sev"]], [yi, yi], color=PALETTE["neutral"], lw=1.8)
        ax_5d1.scatter(rt["early_sev"], yi, s=42, marker="s", c=COL["Ext"] if rt["early_p_sev"] < 0.05 else "white", edgecolors=COL["Ext"], lw=1.4, zorder=4, label="Early ≤7d" if yi == yy_5d[0] else "")
        ax_5d1.scatter(rt["late_sev"], yi, s=36, marker="o", c="white", edgecolors=PALETTE["gray_2"], lw=1.3, zorder=3, label="Late >7d" if yi == yy_5d[0] else "")
    ax_5d1.set_yticks(yy_5d); ax_5d1.set_yticklabels(TS5, fontsize=8.2, fontweight="bold")
    ax_5d1.tick_params(axis="y", length=0); ax_5d1.spines["left"].set_visible(False)
    ax_5d1.set_xlim(-0.08, 0.44); ax_5d1.set_xticks([0.0, 0.2, 0.4]); ax_5d1.set_xlabel("ISG severity ρ", fontsize=8.5)
    ax_5d1.legend(loc="lower left", bbox_to_anchor=(-0.05, 1.005), ncol=2, frameon=False, fontsize=7.0, handletextpad=0.2, columnspacing=0.6)
    ax_5d1.set_title("T-cell ISG: 0–7d vs >7d & fatal collapse", loc="left", fontsize=9.8, fontweight="bold", pad=20)

    ax_5d2.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.2)
    for yi, ts in zip(yy_5d, TS5):
        sub_tg = TISG_GEE[TISG_GEE["subset"] == ts].set_index("cohort")
        for clab, k_s, dy in [("Real", "Real", +0.24), ("PBISC new", "New", 0.0), ("Combined", "HSPC", -0.24)]:
            rg = sub_tg.loc[clab]
            ax_5d2.plot([rg["fat_x_wk_lo"], rg["fat_x_wk_hi"]], [yi + dy, yi + dy], color=COL[k_s], lw=1.5)
            ax_5d2.scatter(rg["fat_x_wk_beta"], yi + dy, s=32, marker=MK["Ext" if k_s == "HSPC" else k_s], c=COL[k_s] if rg["fat_x_wk_p"] < 0.05 else "white", edgecolors=COL[k_s], lw=1.2, zorder=4)
    ax_5d2.set_yticks(yy_5d); ax_5d2.set_yticklabels([]); ax_5d2.tick_params(axis="y", length=0); ax_5d2.spines["left"].set_visible(False)
    ax_5d2.set_xlim(-1.05, 0.10); ax_5d2.set_xticks([-0.8, -0.4, 0.0]); ax_5d2.set_xlabel("GEE Time×Fatal β/wk", fontsize=8.5)

    # Row 3: Panel 5e (Paired n=43 HSPC windows, 5 deaths) & Panel 5f (Source-split GEE HSPC & CD8 CTL)
    ax_5e = fig_s5.add_axes([0.138, 0.300, 0.255, 0.160]); fig_s5.text(0.014, 0.488, "e", fontsize=17.5, fontweight="bold", va="top")
    WIN_LBLS = [("early", "Early (0–7d, n=43)"), ("late", "Late (>7d, n=43)"), ("delta", "Late − Early Δ (43)")]
    yy_w = np.arange(3)[::-1]; p43_ix = P43_WIN.set_index("window"); ax_5e.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for yi, (wk, _) in zip(yy_w, WIN_LBLS):
        rw = p43_ix.loc[wk]; cc_w = COL["HSPC"] if rw["P"] < 0.05 else PALETTE["gray_2"]
        ax_5e.plot([rw["CI_lo"], rw["CI_hi"]], [yi, yi], color=cc_w, lw=2.0)
        ax_5e.scatter(rw["AUC"], yi, s=48, marker="s" if wk != "early" else "D", c=cc_w if rw["P"] < 0.05 else "white", edgecolors=cc_w, linewidths=1.4, zorder=4)
        ax_5e.text(1.02, yi, f"{rw['AUC']:.2f} ({p_eq_s(rw['P'])})", va="center", fontsize=8.0, fontweight="bold" if rw["P"] < 0.05 else "normal", color=cc_w)
    ax_5e.set_yticks(yy_w); ax_5e.set_yticklabels([x[1] for x in WIN_LBLS], fontsize=8.8, fontweight="bold")
    ax_5e.tick_params(axis="y", length=0); ax_5e.spines["left"].set_visible(False)
    ax_5e.set_xlim(0.10, 1.00); ax_5e.set_xticks([0.2, 0.5, 0.8]); ax_5e.set_ylim(-0.55, 2.55)
    ax_5e.set_xlabel("Same-patient HSPC death AUC (5 verified deaths)", fontsize=8.8)
    ax_5e.set_title("Paired within-patient HSPC windows (n = 43)", loc="left", fontsize=10.5, fontweight="bold", pad=8)

    ax_5f = fig_s5.add_axes([0.640, 0.300, 0.205, 0.160]); fig_s5.text(0.518, 0.488, "f", fontsize=17.5, fontweight="bold", va="top")
    GEE_ROWS = [
        ("HSPC", "Real", "HSPC · Real (37, 5d)", COL["Real"], "D"), ("HSPC", "New", "HSPC · PBISC new (29, 6d)", COL["New"], "o"), ("HSPC", "Combined", "HSPC · Combined (66, 11d)", COL["HSPC"], "s"),
        ("T CD8 CTL", "Real", "CD8 CTL · Real (37, 5d)", COL["Real"], "D"), ("T CD8 CTL", "New", "CD8 CTL · PBISC new (29, 6d)", COL["New"], "o"), ("T CD8 CTL", "Combined", "CD8 CTL · Combined (66, 11d)", COL["Ext"], "s"),
    ]
    yy_5f = np.array([5.3, 4.3, 3.3, 2.0, 1.0, 0.0])
    ax_5f.axhline(2.65, color=PALETTE["neutral"], ls="--", lw=1.0); ax_5f.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
    for yi, (fk_g, coh_g, _, cc_g, mk_g) in zip(yy_5f, GEE_ROWS):
        rg = GEE_SRC[(GEE_SRC["feature"] == fk_g) & (GEE_SRC["cohort"] == coh_g)].iloc[0]; sig_g = rg["P"] < 0.05
        ax_5f.plot([rg["CI_lo"], rg["CI_hi"]], [yi, yi], color=cc_g, lw=1.8 if sig_g else 1.3)
        ax_5f.scatter(rg["beta"], yi, s=46, marker=mk_g, c=cc_g if sig_g else "white", edgecolors=cc_g, linewidths=1.5, zorder=4)
        ax_5f.text(1.02, yi, f"{rg['beta']:+.2f} ({p_eq_s(rg['P'])})" if sig_g else f"{rg['beta']:+.2f} (P={rg['P']:.2f})", va="center", fontsize=7.6, fontweight="bold" if sig_g else "normal", color=cc_g if sig_g else PALETTE["gray_1"])
    ax_5f.set_yticks(yy_5f); ax_5f.set_yticklabels([x[2] for x in GEE_ROWS], fontsize=8.3, fontweight="bold")
    ax_5f.tick_params(axis="y", length=0); ax_5f.spines["left"].set_visible(False)
    ax_5f.set_xlim(-0.78, 1.00); ax_5f.set_xticks([-0.5, 0.0, 0.5, 1.0]); ax_5f.set_ylim(-0.55, 5.95)
    ax_5f.set_xlabel("GEE Time × Fatal β / wk (95% CI)", fontsize=8.8)
    ax_5f.set_title("Source-specific vs. combined GEE slope interaction", loc="left", fontsize=10.5, fontweight="bold", pad=8)

    # Row 4: Panel 5g (Hospital labs) & Panel 5h (Bivariate scatters under TRUE_PUB_11!)
    ax_5g = fig_s5.add_axes([0.115, 0.045, 0.365, 0.155]); fig_s5.text(0.014, 0.232, "g", fontsize=17.5, fontweight="bold", va="top")
    IMM6 = ["NK CD56 bright", "T CD8 CM", "T CD8 CTL", "B_ASC_IgG", "B_ASC_IgA", "HSPC"]
    IMM6_LBL = ["NK CD56bright", "T CD8 CM", "T CD8 CTL", "B ASC IgG", "B ASC IgA", "HSPC"]
    LABS6 = ["temp_CRP", "temp_lymph", "NLR", "temp_ferritin", "temp_ddimer", "temp_trop"]
    LABS6_LBL = ["CRP (37/5d)", "ALC (37/5d)", "NLR (21/1d)", "Ferritin (33/5d)", "D-dimer (31/5d)", "Troponin (24/4d)"]
    mat_lc, p_lc = np.zeros((len(IMM6), len(LABS6))), np.zeros((len(IMM6), len(LABS6)))
    for i_f, fk in enumerate(IMM6):
        for j_l, lk in enumerate(LABS6):
            rw = LAB_IND[(LAB_IND["feature"] == fk) & (LAB_IND["lab"] == lk)].iloc[0]
            mat_lc[i_f, j_l], p_lc[i_f, j_l] = rw["rho_with_lab"], rw["p_with_lab"]
    ax_5g.imshow(mat_lc, cmap="RdBu_r", norm=mpl.colors.Normalize(-0.65, 0.65), aspect="auto")
    for i_f in range(len(IMM6)):
        for j_l in range(len(LABS6)):
            vv, sg = mat_lc[i_f, j_l], p_lc[i_f, j_l] < 0.05
            ax_5g.text(j_l, i_f, f"{vv:+.2f}{'*' if sg else ''}", ha="center", va="center", fontsize=8.0, fontweight="bold" if sg else "normal", color="white" if abs(vv) > 0.42 else COL["Real"])
    ax_5g.axhline(4.5, color="black", lw=1.8); ax_5g.set_yticks(range(len(IMM6))); ax_5g.set_yticklabels(IMM6_LBL, fontsize=8.8, fontweight="bold"); ax_5g.get_yticklabels()[-1].set_color(COL["HSPC"])
    ax_5g.set_xticks(range(len(LABS6))); ax_5g.set_xticklabels(LABS6_LBL, rotation=24, ha="left", fontsize=8.2, fontweight="bold")
    ax_5g.xaxis.tick_top(); ax_5g.tick_params(length=0)
    for sp_5g in ax_5g.spines.values(): sp_5g.set_visible(False)
    ax_5g.set_title("Correlation ρ with patient-mean hospital labs (Real; *P < 0.05)", loc="left", fontsize=10.1, fontweight="bold", pad=24)

    ax_5h1 = fig_s5.add_axes([0.555, 0.045, 0.182, 0.155]); ax_5h2 = fig_s5.add_axes([0.790, 0.045, 0.182, 0.155])
    fig_s5.text(0.505, 0.232, "h", fontsize=17.5, fontweight="bold", va="top")
    pos_r_draws = pos_217[pos_217["source"] == "Real"]
    pt_crp = pos_r_draws.groupby("patient_id")[["temp_CRP"]].mean().join(PT_NS.loc[PT_NS["src"] == "Real", ["HSPC", "fat11", "sev"]])
    pt_crp["logCRP"] = np.log2(pt_crp["temp_CRP"] + 1.0)
    for fv, cc, sz_v, lb in [(0, COL["Real"], 32, "Survived (32)"), (1, COL["HSPC"], 64, "Died (5)")]:
        sb = pt_crp[pt_crp["fat11"] == fv]
        ax_5h1.scatter(sb["logCRP"], sb["HSPC"], s=sz_v, marker="D", c=cc, edgecolors="black" if fv == 1 else "white", linewidths=0.8, zorder=4 if fv == 1 else 2, label=lb)
    ax_5h1.axhline(0.5, color=COL["HSPC"], ls="--", lw=1.1, alpha=0.8); ax_5h1.set_ylim(-2.45, 2.65)
    ax_5h1.set_xlabel("Patient-mean log2(CRP + 1)", fontsize=9.2); ax_5h1.set_ylabel("HSPC (Normal score)", fontsize=9.2)
    ax_5h1.legend(loc="upper left", frameon=False, fontsize=7.8); ax_5h1.set_title("Mean CRP vs. HSPC (Real, n = 37)", loc="left", fontsize=10.4, fontweight="bold", pad=6)

    for src_v, mk_v in [("Real", "D"), ("New", "o")]:
        for fv, cc, sz_v in [(0, COL[src_v], 28), (1, COL["HSPC"], 62)]:
            sb = PT_NS[(PT_NS["src"] == src_v) & (PT_NS["fat11"] == fv)]
            ax_5h2.scatter(sb["T CD8 CTL"], sb["HSPC"], s=sz_v, marker=mk_v, c=cc, edgecolors="black" if fv == 1 else "white", linewidths=0.8, alpha=0.65 if fv == 0 else 0.95, zorder=5 if fv == 1 else 2)
    ax_5h2.plot([-2.2, 2.2], [-2.2, 2.2], color=PALETTE["gray_1"], ls=":", lw=1.2)
    ax_biv_auc, ax_biv_p = auc_mw((PT_NS["HSPC"] - PT_NS["T CD8 CTL"]).values, PT_NS["fat11"].values)
    ax_5h2.text(2.28, -2.28, f"HSPC − CD8 CTL (11d)\nAUC = {ax_biv_auc:.2f} ({p_eq_s(ax_biv_p)})",
                fontsize=7.6, fontweight="bold", color=COL["HSPC"], ha="right", va="bottom", bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.88))
    ax_5h2.set_xlim(-2.4, 2.4); ax_5h2.set_ylim(-2.4, 2.55)
    ax_5h2.set_xlabel("T CD8 CTL (Normal score)", fontsize=9.2); ax_5h2.set_ylabel("HSPC (Normal score)", fontsize=9.2)
    ax_5h2.set_title("Bivariate distribution (n = 66, red = 11 Died)", loc="left", fontsize=10.3, fontweight="bold", pad=6)

    saved_s5 = finalize_figure(fig_s5, out_dir / "PBISC_SuppFig_S5_Deconv_HSC_TISG_Clinical_Verified28_20261119", formats=("png", "pdf", "svg"), dpi=300, close=True)
    print("Rendered Supplementary Figure 5:", saved_s5[0].name)

    # -----------------------------------------------------------------
    # Merge 7-page vector PDF bundle & write RUN_RECEIPT.json
    # -----------------------------------------------------------------
    stems = [
        "PBISC_Figure7_Main_Verified28_20261119",
        "PBISC_Figure8_Main_Verified28_20261119",
        "PBISC_SuppFig_S1_Cohort_Mortality_Audit_Verified28_20261119",
        "PBISC_SuppFig_S2_CLL_Anchors_Modality_Verified28_20261119",
        "PBISC_SuppFig_S3_Power_IFN_Sensitivity_Verified28_20261119",
        "PBISC_SuppFig_S4_Public_COVID_Replication_Verified28_20261119",
        "PBISC_SuppFig_S5_Deconv_HSC_TISG_Clinical_Verified28_20261119",
    ]
    bundle_pdf = out_dir / "PBISC_Verified28_True11_Figures_7_8_S1-S5_20261119.pdf"
    from pypdf import PdfWriter
    writer = PdfWriter()
    for st in stems:
        writer.append(str(out_dir / f"{st}.pdf"))
    with open(bundle_pdf, "wb") as f_pdf:
        writer.write(f_pdf)
    print("Merged 7-page vector PDF:", bundle_pdf.name)

    receipt = {
        "mode": "standalone_render_only",
        "mortality_ground_truth": "TRUE_PUB_11 (Real=5 [C21,C33,C141,C108,C145] + PBISC_new=6 [C32,C37,C48,C66,C77,C84]; fixed w2 survivors C101,C146,C147)",
        "data_dir": str(data_dir),
        "output_dir": str(out_dir),
        "renderer_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "style_sha256": hashlib.sha256((SCRIPT_DIR / "style.py").read_bytes()).hexdigest(),
        "data_csv_count": len(list(data_dir.glob("*.csv"))),
        "outputs": [
            {"file": p.name, "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for st in stems for ext in ("png", "pdf", "svg") for p in [out_dir / f"{st}.{ext}"]
        ] + [{"file": bundle_pdf.name, "bytes": bundle_pdf.stat().st_size, "sha256": hashlib.sha256(bundle_pdf.read_bytes()).hexdigest()}],
    }
    (out_dir / "RUN_RECEIPT.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print("Wrote RUN_RECEIPT.json")


if __name__ == "__main__":
    main()
