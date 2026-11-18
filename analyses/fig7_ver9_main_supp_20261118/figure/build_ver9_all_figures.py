"""Build Reorganized Figure 7 (ver9, a-q) + Updated Supplementary Figures 1-4.

Implements the user's exact reorganization:
  - Main Figure 7 (a-q):
      a, b  : Cohort assembly & PBISC clinical event expansion (from former Fig. 0a, b)
      c, d  : CLL molecular screen & 5-step selection funnel (from former Fig. 7a, b)
      e, f  : 6-row x 22-subset fine screen & n80 sample-size dumbbell justifying T CD8, B ASC, and HSPC
              (from former Supplementary Fig. 2a, b)
      g     : Associations gained by adding PBISC (-log10 q Gain Map; former Fig. 7c)
      h     : Credibility scorecard of PBISC-dependent candidates WITH NON-PASTEL COLORS (former Fig. 7d)
      i     : Cohort effect sizes with 95% CIs (former Fig. 7e)
      j-o   : Patient-level within-source Normal score strip plots, top title = cell population name only
              (j: B ASC IgG, k: B ASC IgA, l: T CD8 CTL [severity]; m: T CD8 CM, n: T CD8 CTL, o: HSPC [death])
      p     : HSPC death robustness across 8 mortality definitions / strata (former Fig. 7l)
      q     : Power at 20 Real patients & n80 sample-size reduction (former Fig. 7m)

  - Supplementary Figure 1 (a-d):
      Renamed & adapted from former Fig. 0 (with former 0a, 0b moved to Main Fig. 7a, b):
      Longitudinal sampling timeline (a), serial draw density & wave/steroid breakdown (b),
      cohort characteristics & confounder table (c), and two-stage external validation architecture (d).

  - Supplementary Figure 2 (a-o):
      Renumbered from former Supplementary Fig. 1:
      CLL PBMC distortion (a), CLL before/after effect shift (b, c), 4 concordant severity anchors (d-g),
      NK compartment decomposition (h), 27-shared-donor modality scale offset proof (i),
      and Raw % of PBMCs patient-level distributions matching Main Fig. 7j-o (j-o).

  - Supplementary Figure 3 (a-c):
      Renumbered & rebuilt from former Supplementary Fig. 2 (after moving former a, b into Main Fig. 7e, f):
      Continuous subsampling power curves n=10..35 (a), CD14 mono IFN Step-4 filter audit (b),
      and 14-variant dataset & confounder sensitivity matrix (c).

  - Supplementary Figure 4 (a-d):
      Renumbered from former Supplementary Fig. 3 (external single-cell replication across 4 adult COVID-19
      PBMC cohorts, N = 219 patients, with cross-references updated to Fig. 7h).
"""
import pathlib
import shutil
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
V8_DIR = pathlib.Path(r"F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\15_Figure7_ver8_Main_and_Supp")
PUB_DIR = pathlib.Path(r"F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\16_Public_COVID_Replication_ver8")
OUT_DIR = pathlib.Path(r"F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\17_Figure7_ver9_Main_and_Supp")
OUT_DIR.mkdir(parents=True, exist_ok=True)
shutil.copy2(V8_DIR / "style.py", OUT_DIR / "style.py")
sys.path.insert(0, str(OUT_DIR))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from scipy import stats

from style import CT_COL, FigureStyle, PALETTE, apply_publication_style, finalize_figure

ROOT = pathlib.Path(r"F:\PBISC")
S_PATH = ROOT / "eskd_covid_clinical_assoc_20261115" / "results" / "sample_level_composition_metadata.csv"
RECONC_PATH = ROOT / "paper_claim_analysis" / "covid_sample_patient_id_reconciliation.csv"
DATA_V8 = V8_DIR / "data"
DATA_PUB = PUB_DIR / "data"

PANEL_DIR = OUT_DIR / "panels"
for sub in ("main", "supp_s1", "supp_s2", "supp_s3", "supp_s4"):
    (PANEL_DIR / sub).mkdir(parents=True, exist_ok=True)

apply_publication_style(FigureStyle(font_size=13, axes_linewidth=1.8))

CLL = {"C116", "C139"}
COL = dict(
    Real=PALETTE["gray_3"],
    New=PALETTE["blue_secondary"],
    Ext=PALETTE["blue_main"],
    Shared=PALETTE["teal"],
    Shuf=PALETTE["neutral"],
    All=PALETTE["red_2"],
    CLL=PALETTE["red_strong"],
    # Non-pastel authoritative scorecard colors for Panel h (former Fig. 7d)
    ScoreGood="#1B7A43",   # Deep forest emerald (white text)
    ScoreMid="#B86B11",    # Rich bronze-amber (white text)
    ScoreBad="#B64342",    # Deep crimson (PALETTE['red_strong'], white text)
    # Clean diagram card accents
    CardBlue="#DCE6F2",
    CardGray="#EFEFEF",
    CardTeal="#DCF0F2",
    CardRed="#FADBD8",
    CardGreen=PALETTE["green_1"],
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

ANCHOR_SET = {"NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono"}
FOCUS_SET = {"B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"}


def fp(p):
    if p < 1e-4:
        return "P < 1e-4"
    if p < 0.001:
        return f"P = {p:.1e}".replace("e-0", "e-")
    if p < 0.1:
        return f"P = {p:.3f}"
    return f"P = {p:.2f}"


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


# Load verified precomputed datasets
S_meta = pd.read_csv(S_PATH)
reconc = pd.read_csv(RECONC_PATH)
ster_map = reconc.set_index("sample_id")["published_received_steroid"].eq("Y").to_dict()
S_meta["steroid"] = S_meta.apply(lambda r: bool(ster_map.get(str(r["sample_id"]), False)) if r["source"] == "Real" else False, axis=1)

pos_samps = S_meta[
    (S_meta["condition"] == "positive")
    & S_meta["in_primary_scope"].astype(bool)
    & S_meta["qc_pass_min200cells"].astype(bool)
    & (~S_meta["patient_id"].isin(CLL))
].copy()
r_samps = pos_samps[pos_samps["source"] == "Real"].copy()
b_all_samps = pos_samps[pos_samps["source"] == "PBISC"].copy()
b_new_samps = b_all_samps[~b_all_samps["patient_id"].isin(set(r_samps["patient_id"]))].copy()

CLL_BULK = pd.read_csv(DATA_V8 / "cll_molecular_screen_bulk_donors.csv")
SCAN = pd.read_csv(DATA_V8 / "cll_scan_43fine_v8.csv")
SCAN_NO = SCAN[SCAN["CLL_excluded"] == True].copy()
SEL = pd.read_csv(DATA_V8 / "selection_funnel_43fine.csv")
SC = pd.read_csv(DATA_V8 / "credibility_scorecard_v8.csv")
CI = pd.read_csv(DATA_V8 / "cohort_effect_sizes_ci_v8.csv")
HS = pd.read_csv(DATA_V8 / "hspc_death_analyses_v8.csv")
PW = pd.read_csv(DATA_V8 / "power_curves_v8.csv")
MOD = pd.read_csv(DATA_V8 / "modality_scale_shared27_and_raw_pct_v8.csv")
SH27 = pd.read_csv(DATA_V8 / "shared27_paired_pct_v8.csv")
PT_NS = pd.read_csv(DATA_V8 / "patient_normal_scores_noCLL_v8.csv", index_col=0)
PT_RAW = pd.read_csv(DATA_V8 / "patient_raw_pct_noCLL_v8.csv", index_col=0)
SA = pd.read_csv(DATA_V8 / "sensitivity_14variants_all_v8.csv")
NKD = pd.read_csv(DATA_V8 / "nk_decomposition_v8.csv")

PUB_COH = pd.read_csv(DATA_PUB / "public_cohort_level_results_v8.csv")
PUB_META = pd.read_csv(DATA_PUB / "public_meta_analysis_v8.csv")
PUB_ST_B = pd.read_csv(DATA_PUB / "stephenson_b_asc_isotype_glm_v8.csv")
PUB_HW = pd.read_csv(DATA_PUB / "public_hspc_timewindow_death_v8.csv")


# =========================================================================
# DRAWING FUNCTIONS FOR MAIN FIGURE 7 (a - q)
# =========================================================================
def draw_panel_a_cohort_assembly(ax):
    """Panel a (from former Fig. 0a): Multi-modal ESKD donor architecture & PBISC addition."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    arr_kw = dict(arrowstyle="-|>", lw=1.9, color=PALETTE["gray_2"], mutation_scale=14)

    # Left source card
    rbox(ax, 0.4, 8, 20.6, 80, COL["CardGray"], COL["Real"], lw=1.7)
    ax.text(10.7, 78.5, "Imperial ESKD Cohort\n(105 Bulk Donors)", ha="center", va="center", fontsize=9.8, fontweight="bold", color=COL["Real"])
    ax.text(
        1.7, 53.5,
        "• Hemodialysis (med. 71 yr)\n• 68 COVID-19+ · 37 COVID−",
        ha="left", va="center", fontsize=8.3, color=PALETTE["gray_2"], linespacing=1.32,
    )
    rbox(ax, 1.8, 13.5, 17.8, 21.0, COL["CardRed"], COL["CLL"], lw=1.3, r=1.5)
    ax.text(10.7, 24.0, "ROR1 / FMOD Gate\n−2 CLL (C116, C139)", ha="center", va="center", fontsize=8.2, fontweight="bold", color=COL["CLL"])

    ax.annotate("", xy=(24.0, 48), xytext=(21.4, 48), arrowprops=arr_kw)

    # Center 66-patient assembly block
    rbox(ax, 24.4, 5, 49.4, 85, "white", PALETTE["gray_1"], lw=1.4, ls="--")
    ax.text(
        49.1, 84.0,
        "66 Unique Non-CLL COVID-19+ Patients (217 longitudinal draws)",
        ha="center", va="center", fontsize=9.6, fontweight="bold", color=COL["Ext"],
    )

    # 1) Real-only (10)
    rbox(ax, 25.6, 24, 12.6, 50, "#EAEAEA", COL["Real"], lw=1.6)
    ax.text(31.9, 64.5, "Real-only\nscRNA-seq", ha="center", va="center", fontsize=8.9, fontweight="bold", color=COL["Real"])
    ax.text(31.9, 43.0, "10 patients\n\nCITE-seq\ndroplet only\n(Cambridge 6, NCL 4)", ha="center", va="center", fontsize=8.0, color=PALETTE["gray_2"], linespacing=1.24)

    # 2) Shared Bridge (27)
    rbox(ax, 39.4, 24, 16.2, 50, COL["CardTeal"], COL["Shared"], lw=1.8)
    ax.text(47.5, 64.5, "Shared Bridge\n(Real ∩ PBISC)", ha="center", va="center", fontsize=8.9, fontweight="bold", color="#1F636B")
    ax.text(40.6, 42.5, "• 27 patients\n• Fidelity check\n• Scale audit\n• Counted once", ha="left", va="center", fontsize=8.0, color=COL["Real"], linespacing=1.28)

    # 3) PBISC new (29)
    rbox(ax, 56.8, 24, 15.8, 50, COL["CardBlue"], COL["New"], lw=1.8)
    ax.text(64.7, 64.5, "PBISC new\n(Bulk → scRNA)", ha="center", va="center", fontsize=8.9, fontweight="bold", color=COL["Ext"])
    ax.text(58.0, 42.5, "• 29 patients\n• 79 Wave-1 draws\n• 0% steroid\n• +6 direct deaths", ha="left", va="center", fontsize=8.0, color=COL["Real"], linespacing=1.28)

    ax.plot([25.8, 55.4], [17.5, 17.5], color=COL["Real"], lw=2.2)
    ax.text(40.6, 11.2, "Real: n = 37 (138 draws · 7 deaths)", ha="center", va="center", fontsize=8.3, fontweight="bold", color=COL["Real"])
    ax.plot([39.6, 72.4], [20.4, 20.4], color=COL["New"], lw=1.8, ls=":")
    ax.text(64.7, 11.2, "+ PBISC new: +29 (+78%)", ha="center", va="center", fontsize=8.3, fontweight="bold", color=COL["New"])

    ax.annotate("", xy=(77.0, 48), xytext=(74.2, 48), arrowprops=arr_kw)

    # Right combined cohort card
    rbox(ax, 77.4, 8, 22.0, 80, COL["CardBlue"], COL["Ext"], lw=2.0)
    ax.text(88.4, 77.5, "Harmonized Cohort\nReal + PBISC (n = 66)", ha="center", va="center", fontsize=9.7, fontweight="bold", color=COL["Ext"])
    ax.text(
        78.7, 43.0,
        "• Within-source INT\n  aligns RNA scales\n• +78% sample size\n• Steroid-naive W1\n  enrichment\n• 2× fatal events",
        ha="left", va="center", fontsize=8.3, color=COL["Real"], linespacing=1.32,
    )


def draw_panel_b_event_expansion(ax):
    """Panel b (from former Fig. 0b): Clinical severity & nested mortality expansion."""
    cats_b = [
        ("Mild (WHO 0)", 13, 10),
        ("Moderate (WHO 1)", 9, 5),
        ("Severe (WHO 2)", 10, 7),
        ("Critical (WHO 3)", 5, 7),
        ("All deaths (incl. follow-up)", 7, 6),
        ("Acute COVID (dir. + indir.)", 5, 6),
        ("Direct COVID only", 3, 6),
    ]
    yb = np.arange(len(cats_b))[::-1]
    yb_adj = np.where(yb < 3, yb - 0.38, yb)
    ax.plot([0, 32.0], [2.80, 2.80], color=PALETTE["neutral"], ls="--", lw=1.1, zorder=1)

    for yi, (lab, nr, nn) in zip(yb_adj, cats_b):
        ax.barh(yi, nr, height=0.64, color=COL["Real"], edgecolor="black", lw=1.1, label="Real (37)" if yi == yb_adj[0] else "")
        ax.barh(yi, nn, left=nr, height=0.64, color=COL["New"], edgecolor="black", lw=1.1, label="+ PBISC new (29)" if yi == yb_adj[0] else "")
        pct_gain = int(round(100.0 * nn / nr))
        ax.text(nr + nn + 0.6, yi, f"{nr + nn} (+{pct_gain}%)", va="center", ha="left", fontsize=8.4,
                fontweight="bold", color=PALETTE["red_strong"] if yi < 3 else COL["Ext"])
        ax.text(nr / 2.0, yi, f"{nr}", va="center", ha="center", fontsize=8.0, color="white", fontweight="bold")
        ax.text(nr + nn / 2.0, yi, f"{nn}", va="center", ha="center", fontsize=8.0, color="white", fontweight="bold")

    ax.set_yticks(yb_adj)
    ax.set_yticklabels([c[0] for c in cats_b], fontsize=8.9)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 33.5)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_xlabel("Number of COVID-19+ patients", fontsize=9.8)
    ax.legend(loc="lower right", bbox_to_anchor=(1.02, -0.04), frameon=False, fontsize=8.1, handlelength=1.1, handletextpad=0.35)


def draw_panel_c_cll_screen(ax):
    """Panel c (former Fig. 7a): Bulk PBMC ROR1 vs FMOD CLL molecular screen."""
    df_b = CLL_BULK.set_index("patient_id")
    neg = df_b[~df_b["covid_positive_scope"] & (df_b["cohort_group"] != "Excluded CLL (C116, C139)")]
    pos_real = df_b[df_b["covid_positive_scope"] & (df_b["cohort_group"] == "Real / Shared donor")]
    pos_pbisc = df_b[df_b["covid_positive_scope"] & (df_b["cohort_group"] == "PBISC-only donor")]
    cll_pts = df_b[df_b["cohort_group"] == "Excluded CLL (C116, C139)"]

    ax.scatter(neg["ROR1"], neg["FMOD"], s=22, marker="o", c="white", edgecolors=PALETTE["neutral"], linewidths=0.9, alpha=0.8, zorder=1,
               label=f"COVID− ESKD ({len(neg)})")
    ax.scatter(pos_real["ROR1"], pos_real["FMOD"], s=34, marker="D", c=PALETTE["neutral"], edgecolors=COL["Real"], linewidths=1.0, alpha=0.92, zorder=2,
               label=f"Real/Shared ({len(pos_real)})")
    ax.scatter(pos_pbisc["ROR1"], pos_pbisc["FMOD"], s=40, marker="^", c=COL["New"], edgecolors="black", linewidths=0.8, alpha=0.92, zorder=3,
               label=f"PBISC-only ({len(pos_pbisc)})")
    ax.scatter(cll_pts["ROR1"], cll_pts["FMOD"], s=105, marker="D", c=COL["CLL"], edgecolors="black", linewidths=1.2, zorder=5,
               label="CLL (C116, C139)")

    ax.annotate("C116\n(CLL)", xy=(df_b.loc["C116", "ROR1"], df_b.loc["C116", "FMOD"]),
                xytext=(3.20, 4.30), fontsize=8.8, fontweight="bold", color=COL["CLL"],
                arrowprops=dict(arrowstyle="-", color=COL["CLL"], lw=1.2))
    ax.annotate("C139\n(CLL)", xy=(df_b.loc["C139", "ROR1"], df_b.loc["C139", "FMOD"]),
                xytext=(3.45, 1.95), fontsize=8.8, fontweight="bold", color=COL["CLL"],
                arrowprops=dict(arrowstyle="-", color=COL["CLL"], lw=1.2))

    ax.axhline(1.0, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.75)
    ax.axvline(2.5, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.75)
    ax.set_xlim(-0.3, 6.2)
    ax.set_ylim(-0.35, 6.55)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_yticks([0, 2, 4, 6])
    ax.set_xlabel("Bulk ROR1 log2(CPM+1)", fontsize=9.8)
    ax.set_ylabel("Bulk FMOD log2(CPM+1)", fontsize=9.8)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.00), fontsize=7.7, frameon=False,
              handletextpad=0.22, labelspacing=0.18, borderpad=0.05)


def draw_panel_d_funnel(ax):
    """Panel d (former Fig. 7b): Pre-specified 5-step selection funnel."""
    steps = [
        ("Tested fine types", len(SEL), "≥30% Real detection"),
        ("Associated\n(Real + PBISC new)", int((SEL["S2_sev"] | SEL["S2_death"]).sum()), "sev q<0.05 / death q<0.10"),
        ("Same direction\nin both sources", int((SEL["S3_sev"] | SEL["S3_death"]).sum()), "− NKT (Real ρ ≈ 0)"),
        ("Not centre- or\ntreatment-driven", int((SEL["S4_sev"] | SEL["S4_death"]).sum()), "− CD14 dex/IFN, MAIT"),
        ("PBISC-dependent\n(not in Real alone)", int((SEL["S5_sev"] | SEL["S5_death"]).sum()), "− 4 Real-sig. anchors"),
    ]
    y = np.arange(len(steps))[::-1]
    fills = ["#DCE6F2", "#B3CBE6", "#7FA9D6", COL["New"], COL["Ext"]]
    for yi, (lab, cnt, note), fc in zip(y, steps, fills):
        ax.barh(yi, cnt, height=0.66, color=fc, edgecolor=COL["Real"], linewidth=1.3, zorder=2)
        ax.text(cnt + 1.2, yi, f"{cnt}", va="center", ha="left", fontsize=9.8, fontweight="bold", color=COL["Real"])
        offset_x = 8.0 if cnt >= 10 else 5.8
        ax.text(cnt + offset_x, yi, note, va="center", ha="left", fontsize=7.9,
                color=PALETTE["red_strong"] if note.startswith("−") else PALETTE["gray_2"])
    ax.set_yticks(y)
    ax.set_yticklabels([s[0] for s in steps], fontsize=9.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 86)
    ax.set_xticks([0, 20, 40])
    ax.set_xlabel("Fine-grained cell types", fontsize=9.8)


def draw_panel_e_screen(fig, ax, cax_rect):
    """Panel e (from former Supp Fig. 2a): 6-row x 22-subset overview screen justifying T CD8, B ASC, HSPC."""
    keep_feats = SCAN_NO[(SCAN_NO["q_sev"] < 0.05) | (SCAN_NO["q_death"] < 0.10)]["feature"].unique()
    keep = [f for f in LBL if f in keep_feats and f != "CD16mono|lin"]
    ext = SCAN_NO[SCAN_NO["set"] == "Real + PBISC new"].set_index("feature")
    keep = sorted(keep, key=lambda f: ext.loc[f, "rho"])

    ROWS = [
        ("Real", "rho", "q_sev", 0.05, "Severity · Real", 0.0),
        ("PBISC new", "rho", "q_sev", 0.05, "Severity · PBISC new", -1.0),
        ("Real + PBISC new", "rho", "q_sev", 0.05, "Severity · Real + PBISC new", -2.0),
        ("Real", "AUC_death", "q_death", 0.10, "Death · Real", -3.35),
        ("PBISC new", "AUC_death", "q_death", 0.10, "Death · PBISC new", -4.35),
        ("Real + PBISC new", "AUC_death", "q_death", 0.10, "Death · Real + PBISC new", -5.35),
    ]
    cmap = mpl.colormaps["RdBu_r"]
    norm = mpl.colors.Normalize(-0.8, 0.8)
    ax.axhline(-2.675, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=1)

    for st, vcol, qcol, q_thr, rlab, y_pos in ROWS:
        for i, f in enumerate(keep):
            r = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == st)]
            if not len(r) or pd.isna(r[vcol].iloc[0]):
                continue
            r = r.iloc[0]
            val = r[vcol] if vcol == "rho" else 2.0 * (r[vcol] - 0.5)
            is_fdr = r[qcol] < q_thr
            size = 24 + 135 * min(-np.log10(max(r[qcol], 1e-4)), 4.0) / 4.0
            ax.scatter(i, y_pos, s=size, c=[cmap(norm(val))], edgecolors="black" if is_fdr else "none",
                       linewidths=1.4 if is_fdr else 0.0, zorder=3)

    ax.set_yticks([r[5] for r in ROWS])
    ax.set_yticklabels([r[4] for r in ROWS], fontsize=8.8)
    ax.set_xticks(range(len(keep)))
    ax.set_xticklabels([LBL[f] for f in keep], rotation=44, ha="right", fontsize=8.5)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.6, len(keep) - 0.4)
    ax.set_ylim(-5.95, 1.05)

    for i, f in enumerate(keep):
        t = ax.get_xticklabels()[i]
        if f in ANCHOR_SET:
            t.set_fontweight("bold")
            t.set_color(COL["Real"])
        elif f in FOCUS_SET:
            t.set_fontweight("bold")
            t.set_color(PALETTE["red_strong"] if f == "HSPC" else COL["Ext"])
        p_het = ext.loc[f, "p_het_sev"]
        flag = "S" if f == "CD14mono_dex" else ("H" if (pd.notna(p_het) and p_het < 0.05) else "")
        if flag:
            ax.text(i, 0.68, flag, ha="center", va="center", fontsize=8.6, color=PALETTE["red_strong"], fontweight="bold")

    cax = fig.add_axes(cax_rect)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm, cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-0.8, 0, 0.8])
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.tick_params(labelsize=8.2, length=2, pad=1)
    cb.outline.set_visible(False)


def draw_panel_f_n80(ax):
    """Panel f (from former Supp Fig. 2b): Sample size for 80% power justifying PBISC rescue."""
    ITEMS = [
        ("NK CD56 bright", "sev", "NK CD56bright", "anchor"),
        ("B_ASC_dividing", "sev", "B ASC dividing", "anchor"),
        ("T CD8 CM", "sev", "T CD8 CM (sev.)", "anchor"),
        ("CD16mono", "sev", "CD16 mono", "anchor"),
        ("B_ASC_IgG", "sev", "B ASC IgG (sev.)", "focus"),
        ("B_ASC_IgA", "sev", "B ASC IgA (sev.)", "focus"),
        ("T CD8 CTL", "sev", "T CD8 CTL (sev.)", "focus"),
        ("T CD8 CM", "death", "T CD8 CM (death)", "focus"),
        ("T CD8 CTL", "death", "T CD8 CTL (death)", "focus"),
        ("HSPC", "death", "HSPC (death)", "hspc"),
    ]
    y_positions = []
    y_cur = 0.0
    for idx, (f, ep, disp, grp) in enumerate(ITEMS):
        if idx == 4:
            ax.axhline(y_cur + 0.42, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=0)
            y_cur -= 0.25
        y_positions.append(y_cur)
        d = PW[(PW["feature"] == f) & (PW["outcome"] == ep)].sort_values("n_real")
        a, _ = n80_interp(d["n_real"], d["real_alone"])
        b, le = n80_interp(d["n_real"], d["plus_pbisc_new"])
        ra = a if np.isfinite(a) else 40.0
        ext_col = PALETTE["red_strong"] if grp == "hspc" else (COL["Ext"] if grp == "focus" else COL["Real"])
        ax.plot([b, ra], [y_cur, y_cur], color=PALETTE["neutral"], lw=2.3, zorder=1)
        ax.scatter(ra, y_cur, marker="D", s=42, zorder=3, c=COL["Real"] if np.isfinite(a) else "white", edgecolors=COL["Real"], linewidths=1.5)
        ax.scatter(b, y_cur, marker="s", s=46, zorder=4, c=ext_col, edgecolors=ext_col)
        ax.text(ra + 1.2, y_cur, f"{a:.0f}" if np.isfinite(a) else "n.r.", va="center", fontsize=8.4, color=COL["Real"])
        ax.text(b - 1.2, y_cur, "≤10" if le else f"{b:.0f}", va="center", ha="right", fontsize=8.4, color=ext_col, fontweight="bold")
        y_cur -= 1.0

    ax.set_yticks(y_positions)
    ax.set_yticklabels([it[2] for it in ITEMS], fontsize=8.8)
    for idx_t, t in enumerate(ax.get_yticklabels()):
        grp = ITEMS[idx_t][3]
        if grp != "anchor":
            t.set_fontweight("bold")
            t.set_color(PALETTE["red_strong"] if grp == "hspc" else COL["Ext"])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0.5, 46.5)
    ax.set_xticks([10, 20, 30, 40])
    ax.set_xticklabels(["10", "20", "30", "n.r."])
    ax.set_ylim(y_positions[-1] - 0.55, 0.6)
    ax.set_xlabel("Real patients for 80% power", fontsize=9.8)


def draw_panel_g_gain_map(ax):
    """Panel g (former Fig. 7c): -log10 q in Real alone vs Real + PBISC new."""
    cap = 4.15
    nl = lambda q: np.minimum(-np.log10(np.maximum(q, 1e-6)), cap)
    q05 = -np.log10(0.05)
    q10 = 1.0

    ax.axhspan(q05, cap, xmin=0, xmax=q05 / cap, color=PALETTE["green_1"], alpha=0.85, zorder=0)
    ax.plot([0, cap], [0, cap], color=PALETTE["neutral"], lw=1.3, zorder=1)
    ax.axhline(q05, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1)
    ax.axvline(q05, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1)
    ax.axhline(q10, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.65, zorder=1)
    ax.axvline(q10, color=PALETTE["gray_1"], ls="--", lw=1.0, alpha=0.65, zorder=1)

    for out_k, mk, qr, qe in [("sev", "o", "q_sev_Real", "q_sev_Ext"), ("death", "^", "q_death_Real", "q_death_Ext")]:
        for _, r in SEL.iterrows():
            xv, yv = float(nl(r[qr])), float(nl(r[qe]))
            is_s5 = bool(r[f"S5_{out_k}"])
            is_foc = is_s5 and (r["feature"] in FOCUS_SET)
            is_anc = (out_k == "sev") and bool(r["Anchor_sev"]) and (r["feature"] in ANCHOR_SET)

            if is_foc:
                fc, ec, sz, zo = COL["Ext"], "black", 62, 5
            elif is_anc:
                fc, ec, sz, zo = COL["Real"], "white", 52, 4
            elif is_s5:
                fc, ec, sz, zo = PALETTE["gray_1"], "none", 32, 3
            else:
                fc, ec, sz, zo = PALETTE["neutral"], "none", 20, 2

            ax.scatter(xv, yv, s=sz, marker=mk, c=fc, edgecolors=ec, linewidths=1.0 if ec != "none" else 0.0, zorder=zo)

    callouts = [
        ("B_ASC_IgG", "sev", "q_sev_Real", "q_sev_Ext", "B ASC IgG (sev)", (0.08, 3.68), COL["Ext"]),
        ("B_ASC_IgA", "sev", "q_sev_Real", "q_sev_Ext", "B ASC IgA (sev)", (0.10, 3.12), COL["Ext"]),
        ("HSPC", "death", "q_death_Real", "q_death_Ext", "HSPC (death)", (0.08, 2.56), PALETTE["red_strong"]),
        ("T CD8 CTL", "sev", "q_sev_Real", "q_sev_Ext", "T CD8 CTL (sev)", (0.08, 1.96), COL["Ext"]),
        ("T CD8 CTL", "death", "q_death_Real", "q_death_Ext", "T CD8 CTL (death)", (0.18, 0.40), COL["Ext"]),
        ("T CD8 CM", "death", "q_death_Real", "q_death_Ext", "T CD8 CM (death)", (1.42, 0.72), COL["Ext"]),
        ("NK CD56 bright", "sev", "q_sev_Real", "q_sev_Ext", "NK CD56bright", (2.32, 3.88), COL["Real"]),
        ("T CD8 CM", "sev", "q_sev_Real", "q_sev_Ext", "T CD8 CM (sev)", (2.48, 3.00), COL["Real"]),
        ("B_ASC_dividing", "sev", "q_sev_Real", "q_sev_Ext", "B ASC dividing", (2.46, 2.45), COL["Real"]),
        ("CD16mono", "sev", "q_sev_Real", "q_sev_Ext", "CD16 mono", (2.35, 1.70), COL["Real"]),
    ]
    sel_ix = SEL.set_index("feature")
    for f, _, qr, qe, txt, xy_t, tcol in callouts:
        rr = sel_ix.loc[f]
        xv, yv = float(nl(rr[qr])), float(nl(rr[qe]))
        ax.annotate(
            txt, xy=(xv, yv), xytext=xy_t, fontsize=7.8, fontweight="bold", color=tcol, va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.82),
            arrowprops=dict(arrowstyle="-", color=PALETTE["gray_1"], lw=0.9, shrinkA=1, shrinkB=4),
            zorder=6,
        )

    ax.set_xlim(0, cap)
    ax.set_ylim(0, cap)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlabel("−log10 q, Real alone", fontsize=9.8)
    ax.set_ylabel("−log10 q, Real + PBISC new", fontsize=9.8)


def draw_panel_h_scorecard_nonpastel(ax):
    """Panel h (former Fig. 7d): Credibility scorecard with NON-PASTEL saturated colors & white text."""
    cols = [
        "Real\n(ρ / AUC)",
        "PBISC new\n(ρ / AUC)",
        "Homog.\nP_het",
        "Shared ρ\n(n = 27)",
        "Robust\n(of 14)",
        "WHO 2–3\ndeath AUC",
    ]
    GOOD, MID, BAD = COL["ScoreGood"], COL["ScoreMid"], COL["ScoreBad"]

    for i, (_, r) in enumerate(SC.iterrows()):
        o = r["outcome"]
        fmt = (lambda z: f"{z:+.2f}") if o == "sev" else (lambda z: f"{z:.2f}")
        sh_txt = f"{r['shared_rho']:+.2f}{'*' if r['shared_p'] < 0.05 else ''}"
        rob_txt = f"{int(r['robust_pass'])}/{int(r['robust_total'])}"
        cells = [
            (fmt(r["effect_Real"]), "white", COL["Real"]),
            (fmt(r["effect_New"]), "white", COL["New"]),
            (f"{r['p_het']:.2f}", GOOD if r["p_het"] >= 0.05 else BAD, "white"),
            (sh_txt, GOOD if r["shared_p"] < 0.05 else (MID if r["shared_rho"] >= 0.25 else BAD), "white"),
            (rob_txt, GOOD if r["robust_pass"] >= 12 else (MID if r["robust_pass"] >= 9 else BAD), "white"),
            (
                (f"{r['AUC_within_WHO23']:.2f}", GOOD if r["p_within_WHO23"] < 0.05 else BAD, "white")
                if o == "death"
                else ("–", "white", PALETTE["gray_1"])
            ),
        ]
        # Subtle alternating neutral row stripe behind cols 0-1 for crisp table structure
        if i % 2 == 1:
            ax.add_patch(Rectangle((0.0, i + 0.05), 2.0, 0.90, fc="#F5F6F8", ec="none", zorder=0))

        for j, (txt, fc, tc) in enumerate(cells):
            if fc != "white":
                ax.add_patch(FancyBboxPatch(
                    (j + 0.06, i + 0.10), 0.88, 0.80,
                    boxstyle="round,pad=0.02,rounding_size=0.10",
                    fc=fc, ec="#222222", lw=0.45, zorder=2,
                ))
            ax.text(
                j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=8.9, color=tc,
                fontweight="bold" if (fc != "white" or (j in (0, 1) and r["followed"])) else "normal",
                zorder=3,
            )

    ax.set_xlim(0, len(cols))
    ax.set_ylim(len(SC), 0)
    ax.set_yticks(np.arange(len(SC)) + 0.5)
    ylabs = [f"{LBL[r.feature]} · {'sev.' if r.outcome == 'sev' else 'death'}" for _, r in SC.iterrows()]
    ax.set_yticklabels(ylabs, fontsize=9.0)
    for idx_t, t in enumerate(ax.get_yticklabels()):
        if SC.iloc[idx_t]["followed"]:
            t.set_fontweight("bold")
            t.set_color(PALETTE["red_strong"] if SC.iloc[idx_t]["feature"] == "HSPC" else COL["Ext"])
        else:
            t.set_color(PALETTE["gray_2"])

    ax.set_xticks(np.arange(len(cols)) + 0.5)
    ax.set_xticklabels(cols, fontsize=8.6, fontweight="bold")
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)

    ax.plot([-0.34, 1.0], [6.0, 6.0], transform=ax.get_yaxis_transform(), color=COL["Real"], lw=1.5, clip_on=False)
    ax.text(-0.33, 3.0, "Followed", transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center",
            fontsize=9.0, fontweight="bold", color=COL["Ext"], clip_on=False)
    ax.text(-0.33, 8.0, "Not followed", transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center",
            fontsize=8.6, color=PALETTE["gray_2"], clip_on=False)


def draw_panel_i_forest(ax_sev, ax_dth):
    """Panel i (former Fig. 7e): Forest plots of Spearman rho & Death AUC."""
    FEATS_E = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
    ci_no = CI[CI["CLL_excluded"] == True]

    y_centers = []
    y_cur = 0.0
    for idx, f in enumerate(FEATS_E):
        if idx == 2:
            for axx in (ax_sev, ax_dth):
                axx.axhline(y_cur + 0.52, color=PALETTE["neutral"], ls="--", lw=1.1, zorder=1)
            y_cur -= 0.35
        y_centers.append(y_cur)
        offsets = [("Real", "Real", +0.26, 36), ("PBISC new", "New", 0.0, 38), ("Real + PBISC new", "Ext", -0.26, 42)]
        for sname, k, dy, sz in offsets:
            row = ci_no[(ci_no["feature"] == f) & (ci_no["set"] == sname)].iloc[0]
            ax_sev.plot([row["rho_lo"], row["rho_hi"]], [y_cur + dy, y_cur + dy], color=COL[k], lw=1.5, alpha=0.8, zorder=2)
            ax_sev.scatter(row["rho"], y_cur + dy, s=sz, marker=MK[k],
                           c="white" if row["p_sev"] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.4, zorder=3)
            ax_dth.plot([row["AUC_lo"], row["AUC_hi"]], [y_cur + dy, y_cur + dy], color=COL[k], lw=1.5, alpha=0.8, zorder=2)
            ax_dth.scatter(row["AUC"], y_cur + dy, s=sz, marker=MK[k],
                           c="white" if row["p_death"] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.4, zorder=3)
        y_cur -= 1.05

    ax_sev.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
    ax_dth.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.3)

    ax_sev.set_yticks(y_centers)
    ax_sev.set_yticklabels([LBL[f] for f in FEATS_E], fontsize=9.2)
    for idx_t, t in enumerate(ax_sev.get_yticklabels()):
        f = FEATS_E[idx_t]
        t.set_fontweight("bold")
        t.set_color(COL["Real"] if idx_t < 2 else (PALETTE["red_strong"] if f == "HSPC" else COL["Ext"]))
    ax_dth.set_yticks(y_centers)
    ax_dth.set_yticklabels([])

    for axx in (ax_sev, ax_dth):
        axx.tick_params(axis="y", length=0)
        axx.spines["left"].set_visible(False)
        axx.set_ylim(y_centers[-1] - 0.55, 0.55)

    ax_sev.set_xlim(-0.92, 0.92)
    ax_sev.set_xticks([-0.5, 0.0, 0.5])
    ax_sev.set_xlabel("Spearman ρ (severity)", fontsize=9.8)

    ax_dth.set_xlim(0.0, 1.05)
    ax_dth.set_xticks([0.25, 0.50, 0.75, 1.0])
    ax_dth.set_xticklabels(["0.25", "0.5", "0.75", "1"])
    ax_dth.set_xlabel("AUC (death)", fontsize=9.8)


def draw_strip_ns(ax, f, outcome, rng_local, show_ylabel=True):
    """Patient-level within-source Normal score strip plots (j-o)."""
    if outcome == "sev":
        for g in range(4):
            for k, dx in [("Real", -0.18), ("New", 0.18)]:
                vals = PT_NS[(PT_NS["src"] == k) & (PT_NS["sev"] == g)][f].dropna().values
                x = g + dx + rng_local.uniform(-0.075, 0.075, len(vals))
                ax.scatter(x, vals, s=24, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                if len(vals):
                    ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.9, zorder=4)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=9.3)
        ax.set_xlim(-0.55, 3.55)
    else:
        for g in range(2):
            for k, dx in [("Real", -0.18), ("New", 0.18)]:
                vals = PT_NS[(PT_NS["src"] == k) & (PT_NS["fatal"] == g)][f].dropna().values
                x = g + dx + rng_local.uniform(-0.075, 0.075, len(vals))
                ax.scatter(x, vals, s=26, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                if len(vals):
                    ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.9, zorder=4)
        n_surv = int((PT_NS["fatal"] == 0).sum())
        n_died = int((PT_NS["fatal"] == 1).sum())
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Survived\n({n_surv})", f"Died\n({n_died})"], fontsize=9.2)
        ax.set_xlim(-0.52, 1.52)

    ax.set_ylim(-2.6, 5.15)
    ax.set_yticks([-2, 0, 2])
    if show_ylabel:
        ax.set_ylabel("Normal score (within-source)", fontsize=9.6)

    for idx_s, (sname, k, disp_s) in enumerate([("Real", "Real", "Real"), ("PBISC new", "New", "PBISC new"), ("Real + PBISC new", "Ext", "Combined")]):
        r = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == sname)].iloc[0]
        if outcome == "sev":
            q_s = f" (q={r.q_sev:.3f})" if (sname == "Real + PBISC new" and r.q_sev >= 0.001) else (f" (q={r.q_sev:.1e})" if sname == "Real + PBISC new" else "")
            txt = f"{disp_s}: ρ={r.rho:+.2f}, {fp(r.p_sev)}{q_s}".replace("e-0", "e-")
            sig = r.p_sev < 0.05
        else:
            q_s = f" (q={r.q_death:.2f})" if sname == "Real + PBISC new" else ""
            txt = f"{disp_s}: AUC={r.AUC_death:.2f}, {fp(r.p_death)}{q_s}"
            sig = r.p_death < 0.05
        ax.text(0.03, 0.975 - 0.092 * idx_s, txt, transform=ax.transAxes, fontsize=8.0,
                color=COL[k], va="top", fontweight="bold" if sig else "normal")


def draw_panel_p_hspc_analyses(ax):
    """Panel p (former Fig. 7l): HSPC mortality AUC across 8 definitions / strata."""
    ORDER_H = [
        ("All deaths", "All deaths (incl. late)"),
        ("Table S1 deaths", "Acute COVID (dir.+indir.)"),
        ("Direct COVID-19 deaths", "Direct COVID only"),
        ("Within WHO 2–3", "Within WHO 2–3"),
        ("Severity-stratified", "Severity-stratified"),
        ("Cambridge only (Real)", "Cambridge only"),
        ("Samples > day 7", "Samples > day 7"),
        ("PBISC new alone", "PBISC new alone"),
    ]
    yy = np.arange(len(ORDER_H))[::-1]
    ax.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=1)

    for yi, (an_key, _) in zip(yy, ORDER_H):
        if an_key == "PBISC new alone":
            rn = HS[HS["analysis"] == an_key].iloc[0]
            ax.scatter(rn["AUC"], yi, s=46, marker="o", c=COL["New"], edgecolors=COL["New"], linewidths=1.4, zorder=4)
            ax.text(1.02, yi, f"{int(rn['deaths'])} · {fp(rn['p']).replace('P = ', '')}", va="center", ha="left",
                    fontsize=7.9, color=COL["New"], fontweight="bold")
        else:
            rr = HS[(HS["analysis"] == an_key) & (HS["set"] == "Real")].iloc[0]
            re = HS[(HS["analysis"] == an_key) & (HS["set"] == "Real + PBISC new")].iloc[0]
            ax.plot([rr["AUC"], re["AUC"]], [yi + 0.16, yi - 0.16], color=PALETTE["neutral"], lw=1.2, zorder=2)
            ax.scatter(rr["AUC"], yi + 0.16, s=40, marker="D",
                       c="white" if rr["p"] >= 0.05 else COL["Real"], edgecolors=COL["Real"], linewidths=1.4, zorder=3)
            ax.scatter(re["AUC"], yi - 0.16, s=44, marker="s",
                       c="white" if re["p"] >= 0.05 else COL["Ext"], edgecolors=COL["Ext"], linewidths=1.4, zorder=4)
            p_str = fp(re["p"]).replace("P = ", "")
            ax.text(1.02, yi, f"{int(rr['deaths'])}→{int(re['deaths'])} · {p_str}", va="center", ha="left",
                    fontsize=7.9, color=COL["Ext"], fontweight="bold" if re["p"] < 0.05 else "normal")

    ax.text(1.02, yy[0] + 0.72, "deaths · +PBISC P", va="bottom", ha="left", fontsize=7.7, color=PALETTE["gray_2"], fontweight="bold")
    ax.set_yticks(yy)
    ax.set_yticklabels([lab for _, lab in ORDER_H], fontsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0.46, 1.00)
    ax.set_ylim(-0.6, len(ORDER_H) - 0.25)
    ax.set_xticks([0.5, 0.75, 1.0])
    ax.set_xticklabels(["0.5", "0.75", "1.0"])
    ax.set_xlabel("HSPC death AUC", fontsize=9.6)


def draw_panel_q_power(ax):
    """Panel q (former Fig. 7m): Power at 20 Real patients & n80 reduction."""
    PP = [
        ("B_ASC_IgG", "sev", "B ASC IgG · sev."),
        ("B_ASC_IgA", "sev", "B ASC IgA · sev."),
        ("T CD8 CTL", "sev", "T CD8 CTL · sev."),
        ("T CD8 CM", "death", "T CD8 CM · death"),
        ("T CD8 CTL", "death", "T CD8 CTL · death"),
        ("HSPC", "death", "HSPC · death"),
    ]
    yy = np.arange(len(PP))[::-1]
    hb = 0.25
    for col_k, fc, dy, hatch in [
        ("real_alone", PALETTE["neutral"], +hb, ""),
        ("plus_pbisc_new", COL["Ext"], 0.0, ""),
        ("plus_shuffled", "white", -hb, "///"),
    ]:
        vals = [PW[(PW["feature"] == f) & (PW["outcome"] == o) & (PW["n_real"] == 20)][col_k].iloc[0] for f, o, _ in PP]
        ax.barh(yy + dy, vals, height=hb, color=fc, edgecolor=COL["Real"], linewidth=1.0, hatch=hatch, zorder=3)

    for yi, (f, o, _) in zip(yy, PP):
        d_sub = PW[(PW["feature"] == f) & (PW["outcome"] == o)].sort_values("n_real")
        n_r, _ = n80_interp(d_sub["n_real"], d_sub["real_alone"])
        n_e, le = n80_interp(d_sub["n_real"], d_sub["plus_pbisc_new"])
        r_txt = f"{n_r:.0f}" if np.isfinite(n_r) else "n.r."
        e_txt = "≤10" if le else f"{n_e:.0f}"
        ax.text(1.03, yi, f"{r_txt}→{e_txt}", va="center", ha="left", fontsize=7.9, color=COL["Ext"], fontweight="bold")

    ax.text(1.03, yy[0] + 0.65, "n80 Real", va="bottom", ha="left", fontsize=7.7, color=PALETTE["gray_2"], fontweight="bold")
    ax.axvline(0.80, color=PALETTE["gray_1"], ls=":", lw=1.3, zorder=2)
    ax.set_yticks(yy)
    ax.set_yticklabels([lab for _, _, lab in PP], fontsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0.0, 1.01)
    ax.set_ylim(-0.55, len(PP) - 0.2)
    ax.set_xticks([0.0, 0.5, 0.8, 1.0])
    ax.set_xticklabels(["0", "0.5", "0.8", "1"])
    ax.set_xlabel("Power at 20 Real pts", fontsize=9.6)


# =========================================================================
# 1) RENDER REORGANIZED MAIN FIGURE 7 (a - q)
# =========================================================================
fig_m = plt.figure(figsize=(17.6, 25.8), facecolor="white")
fig_m.text(
    0.012,
    0.992,
    "Figure 7. Multi-modal cohort extension, pre-specified screening, and validation of PBISC-rescued immune associations in ESKD COVID-19",
    fontsize=14.6,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

# --- Row 1: a (Cohort assembly) & b (Severity & fatal event expansion) ---
ax_a = fig_m.add_axes([0.016, 0.838, 0.620, 0.124])
draw_panel_a_cohort_assembly(ax_a)
ax_a.set_title("Multi-modal ESKD COVID-19 cohort assembly via PBISC", x=0.042, ha="left", fontsize=11.2, fontweight="bold", pad=4)
fig_m.text(0.012, 0.974, "a", fontsize=17.5, fontweight="bold", va="top")

ax_b = fig_m.add_axes([0.775, 0.842, 0.198, 0.116])
draw_panel_b_event_expansion(ax_b)
ax_b.set_title("Severity & fatal event expansion", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_m.text(0.652, 0.974, "b", fontsize=17.5, fontweight="bold", va="top")

# --- Row 2: c (CLL molecular screen), d (Selection steps), e (6-row x 22-subset fine screen) ---
ax_c = fig_m.add_axes([0.048, 0.668, 0.156, 0.112])
draw_panel_c_cll_screen(ax_c)
ax_c.set_title("CLL molecular screen", loc="left", fontsize=11.0, fontweight="bold", pad=6)
fig_m.text(0.012, 0.800, "c", fontsize=17.5, fontweight="bold", va="top")

ax_d = fig_m.add_axes([0.318, 0.668, 0.158, 0.112])
draw_panel_d_funnel(ax_d)
ax_d.set_title("Selection steps (43 cell types)", loc="left", fontsize=11.0, fontweight="bold", pad=6)
fig_m.text(0.224, 0.800, "d", fontsize=17.5, fontweight="bold", va="top")

ax_e = fig_m.add_axes([0.612, 0.668, 0.372, 0.104])
draw_panel_e_screen(fig_m, ax_e, cax_rect=[0.612, 0.788, 0.088, 0.007])
fig_m.text(
    0.708, 0.791,
    "ρ or 2·(AUC−0.5) · size: −log10 q · ring: FDR sig. · H: het. · S: steroid",
    fontsize=8.1, va="center", color=COL["Real"],
)
fig_m.text(0.495, 0.800, "e", fontsize=17.5, fontweight="bold", va="top")

# --- Master Legend Strip between Row 2 and Row 3 ---
ax_leg = fig_m.add_axes([0.012, 0.586, 0.976, 0.024])
ax_leg.set_axis_off()
hd_main = [
    Line2D([0], [0], marker="D", ls="", ms=7.2, mfc=COL["Real"], mec=COL["Real"], label="Real (37 patients, 7 deaths)"),
    Line2D([0], [0], marker="o", ls="", ms=7.2, mfc=COL["New"], mec=COL["New"], label="PBISC new (29 patients, 6 deaths)"),
    Line2D([0], [0], marker="s", ls="", ms=7.2, mfc=COL["Ext"], mec=COL["Ext"], label="Real + PBISC new (66 patients, 13 deaths)"),
    Line2D([0], [0], marker="o", ls="", ms=7.2, mfc="white", mec=PALETTE["gray_2"], mew=1.5, label="open: P ≥ 0.05 (f, i, p)"),
    Line2D([0], [0], marker="o", ls="", ms=6.2, mfc=COL["Ext"], mec="black", label="severity (g)"),
    Line2D([0], [0], marker="^", ls="", ms=6.8, mfc=COL["Ext"], mec="black", label="death (g)"),
    Patch(facecolor="white", edgecolor=COL["Real"], hatch="///", label="PBISC new, outcome shuffled (q)"),
    Patch(facecolor=COL["ScoreGood"], edgecolor="#222222", lw=0.5, label="criterion met (h)"),
    Patch(facecolor=COL["ScoreMid"], edgecolor="#222222", lw=0.5, label="borderline (h)"),
    Patch(facecolor=COL["ScoreBad"], edgecolor="#222222", lw=0.5, label="not met (h)"),
]
ax_leg.legend(handles=hd_main, loc="center", ncol=5, fontsize=9.2, frameon=False, columnspacing=1.2, handletextpad=0.35, labelspacing=0.35)

# --- Row 3: f (Sample size for 80% power), g (-log10 q gain map), h (Non-pastel credibility scorecard) ---
ax_f = fig_m.add_axes([0.114, 0.418, 0.156, 0.122])
draw_panel_f_n80(ax_f)
ax_f.set_title("Sample size for 80% power", loc="left", fontsize=11.0, fontweight="bold", pad=6)
fig_m.text(0.012, 0.562, "f", fontsize=17.5, fontweight="bold", va="top")

ax_g = fig_m.add_axes([0.330, 0.418, 0.180, 0.122])
draw_panel_g_gain_map(ax_g)
ax_g.set_title("Associations gained by adding PBISC", loc="left", fontsize=11.0, fontweight="bold", pad=6)
fig_m.text(0.288, 0.562, "g", fontsize=17.5, fontweight="bold", va="top")

ax_h = fig_m.add_axes([0.642, 0.414, 0.342, 0.118])
draw_panel_h_scorecard_nonpastel(ax_h)
ax_h.set_title("Credibility of PBISC-dependent candidates", loc="left", fontsize=11.0, fontweight="bold", pad=28)
fig_m.text(0.532, 0.562, "h", fontsize=17.5, fontweight="bold", va="top")

# --- Row 4: i (Forest CIs) + j, k, l (Rescued Severity Strip Plots; top title = cell name ONLY) ---
ax_i_sev = fig_m.add_axes([0.095, 0.228, 0.120, 0.122])
ax_i_dth = fig_m.add_axes([0.232, 0.228, 0.120, 0.122])
draw_panel_i_forest(ax_i_sev, ax_i_dth)
ax_i_sev.set_title("Cohort effect sizes (95% CI)", loc="left", fontsize=11.0, fontweight="bold", pad=8)
fig_m.text(0.012, 0.370, "i", fontsize=17.5, fontweight="bold", va="top")

rng_main = np.random.default_rng(7)
for idx_p, (fk, let, show_yl, x_l) in enumerate([
    ("B_ASC_IgG", "j", True, 0.432),
    ("B_ASC_IgA", "k", False, 0.628),
    ("T CD8 CTL", "l", False, 0.824),
]):
    ax_s = fig_m.add_axes([x_l, 0.228, 0.152, 0.122])
    draw_strip_ns(ax_s, fk, "sev", rng_main, show_ylabel=show_yl)
    ax_s.set_title(LBL[fk], loc="left", fontsize=11.0, fontweight="bold", color=COL["Real"], pad=5)
    fig_m.text(x_l - (0.044 if show_yl else 0.028), 0.370, let, fontsize=17.5, fontweight="bold", va="top")

# --- Row 5: m, n, o (Rescued Death Strip Plots; top title = cell name ONLY) + p (HSPC robustness) + q (Power & n80) ---
for idx_d, (fk, let, show_yl, x_l) in enumerate([
    ("T CD8 CM", "m", True, 0.048),
    ("T CD8 CTL", "n", False, 0.205),
    ("HSPC", "o", False, 0.362),
]):
    ax_s = fig_m.add_axes([x_l, 0.040, 0.125, 0.122])
    draw_strip_ns(ax_s, fk, "death", rng_main, show_ylabel=show_yl)
    ax_s.set_title(LBL[fk], loc="left", fontsize=11.0, fontweight="bold", color=COL["Real"], pad=5)
    fig_m.text(x_l - (0.036 if show_yl else 0.024), 0.182, let, fontsize=17.5, fontweight="bold", va="top")

ax_p = fig_m.add_axes([0.614, 0.040, 0.088, 0.122])
draw_panel_p_hspc_analyses(ax_p)
ax_p.set_title("HSPC death robustness", loc="left", fontsize=11.0, fontweight="bold", color=COL["Real"], pad=5)
fig_m.text(0.502, 0.182, "p", fontsize=17.5, fontweight="bold", va="top")

ax_q = fig_m.add_axes([0.872, 0.040, 0.074, 0.122])
draw_panel_q_power(ax_q)
ax_q.set_title("Power & n80 gain", loc="left", fontsize=11.0, fontweight="bold", color=COL["Real"], pad=5)
fig_m.text(0.776, 0.182, "q", fontsize=17.5, fontweight="bold", va="top")

main_out = OUT_DIR / "PBISC_Figure7_Main_ver9_20261118"
saved_main = finalize_figure(fig_m, main_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Main Figure 7 (ver9, a-q):", [str(p) for p in saved_main])


# =========================================================================
# 2) SUPPLEMENTARY FIGURE 1 (Renamed & Adapted from Former Fig. 0)
# =========================================================================
fig_s1 = plt.figure(figsize=(16.6, 12.2), facecolor="white")
fig_s1.text(
    0.015,
    0.985,
    "Supplementary Figure 1. Longitudinal sampling dynamics, confounder profile, and two-stage external validation architecture",
    fontsize=14.0,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

# Panel a: Longitudinal sampling timeline (217 draws)
ax_s1a = fig_s1.add_axes([0.138, 0.560, 0.220, 0.345])
fig_s1.text(0.015, 0.935, "a", fontsize=18, fontweight="bold", va="top")
ax_s1a.axvspan(7.0, 31.5, color=COL["CardRed"], alpha=0.38, zorder=0)
ax_s1a.axvline(7.0, color=COL["CLL"], ls="--", lw=1.3, zorder=2)
ax_s1a.text(19.2, 4.18, "Late window (> day 7)\nsustained HSPC surge", ha="center", va="top", fontsize=8.6, fontweight="bold", color=COL["CLL"])
ax_s1a.text(3.3, 4.18, "Early\n(≤7d)", ha="center", va="top", fontsize=8.4, color=PALETTE["gray_2"])

rng_c = np.random.default_rng(12)
groups_c = [
    ("PBISC new · Died (6 pts)", b_new_samps[b_new_samps["fatal_harmonized"] == 1]["days_from_first_positive_swab"], 0.0, "o", COL["CLL"]),
    ("PBISC new · Survived (23)", b_new_samps[b_new_samps["fatal_harmonized"] == 0]["days_from_first_positive_swab"], 1.0, "o", COL["New"]),
    ("Real · Died (7 pts)", r_samps[r_samps["fatal_harmonized"] == 1]["days_from_first_positive_swab"], 2.1, "D", COL["CLL"]),
    ("Real · Survived (30)", r_samps[r_samps["fatal_harmonized"] == 0]["days_from_first_positive_swab"], 3.1, "D", COL["Real"]),
]
for lab, s_days, y0, mk, c_dot in groups_c:
    dv = s_days.dropna().values
    jit = rng_c.uniform(-0.18, 0.18, len(dv))
    ax_s1a.scatter(dv, y0 + jit, s=26, marker=mk, c=c_dot, edgecolors="white", linewidths=0.45, alpha=0.9, zorder=3)
    med = float(np.median(dv))
    q25, q75 = float(np.percentile(dv, 25)), float(np.percentile(dv, 75))
    ax_s1a.plot([q25, q75], [y0 - 0.30, y0 - 0.30], color="black", lw=2.0, zorder=4)
    ax_s1a.scatter([med], [y0 - 0.30], s=32, marker="|", c="black", linewidths=2.2, zorder=5)

ax_s1a.set_yticks([g[2] for g in groups_c])
ax_s1a.set_yticklabels([g[0] for g in groups_c], fontsize=9.2)
ax_s1a.tick_params(axis="y", length=0, pad=5)
ax_s1a.spines["left"].set_visible(False)
ax_s1a.set_xlim(-1.0, 31.5)
ax_s1a.set_ylim(-0.55, 4.30)
ax_s1a.set_xticks([0, 7, 14, 21, 30])
ax_s1a.set_xlabel("Days from first positive swab", fontsize=10.2)
ax_s1a.set_title("Longitudinal sampling (217 draws)", loc="left", fontsize=11.4, fontweight="bold", pad=8)

# Panel b: Serial draw density & wave/steroid breakdown across the 66 patients
ax_s1b = fig_s1.add_axes([0.468, 0.560, 0.165, 0.345])
fig_s1.text(0.382, 0.935, "b", fontsize=18, fontweight="bold", va="top")

r_cnts = r_samps.groupby("patient_id").size()
n_cnts = b_new_samps.groupby("patient_id").size()
r_pts_meta = r_samps.groupby("patient_id").agg(batch=("centre_or_wave", "first"), ster=("steroid", "max"))

cats_s1b = [
    ("1 draw", int((r_cnts == 1).sum()), int((n_cnts == 1).sum())),
    ("2 draws", int((r_cnts == 2).sum()), int((n_cnts == 2).sum())),
    ("3 draws", int((r_cnts == 3).sum()), int((n_cnts == 3).sum())),
    ("≥4 draws (up to 7)", int((r_cnts >= 4).sum()), int((n_cnts >= 4).sum())),
    ("NCL centre", int((r_pts_meta["batch"] == "NCL").sum()), 0),
    ("Cambridge centre", int((r_pts_meta["batch"] == "Cambridge").sum()), 0),
    ("Wave-1 bulk (w1)", 0, len(n_cnts)),
    ("Steroid-exposed", int(r_pts_meta["ster"].sum()), 0),
]
yb_s1 = np.arange(len(cats_s1b))[::-1]
yb_s1_adj = np.where(yb_s1 < 4, yb_s1 - 0.40, yb_s1)
ax_s1b.plot([0, 34], [3.78, 3.78], color=PALETTE["neutral"], ls="--", lw=1.1, zorder=1)

for yi, (lab, nr, nn) in zip(yb_s1_adj, cats_s1b):
    if nr > 0:
        ax_s1b.barh(yi, nr, height=0.62, color=COL["Real"], edgecolor="black", lw=1.0, label="Real (37)" if yi == yb_s1_adj[0] else "")
    if nn > 0:
        ax_s1b.barh(yi, nn, left=nr, height=0.62, color=COL["New"], edgecolor="black", lw=1.0, label="+ PBISC new (29)" if yi == yb_s1_adj[0] else "")
    ax_s1b.text(nr + nn + 0.7, yi, f"{nr + nn}", va="center", ha="left", fontsize=8.8, fontweight="bold", color=COL["Real"])

ax_s1b.set_yticks(yb_s1_adj)
ax_s1b.set_yticklabels([c[0] for c in cats_s1b], fontsize=9.1)
ax_s1b.tick_params(axis="y", length=0, pad=5)
ax_s1b.spines["left"].set_visible(False)
ax_s1b.set_xlim(0, 36.0)
ax_s1b.set_xticks([0, 10, 20, 30])
ax_s1b.set_xlabel("COVID-19+ patients", fontsize=10.2)
ax_s1b.set_title("Serial draws & centre/steroid strata", loc="left", fontsize=11.4, fontweight="bold", pad=8)
ax_s1b.legend(loc="lower right", frameon=False, fontsize=8.5, handlelength=1.1, handletextpad=0.35)

# Panel c: Cohort characteristics & confounder profile table (bounded inside cells)
ax_s1c = fig_s1.add_axes([0.678, 0.555, 0.304, 0.355])
ax_s1c.set_xlim(0, 10.0)
ax_s1c.set_ylim(0, 7.55)
ax_s1c.axis("off")
fig_s1.text(0.654, 0.935, "c", fontsize=18, fontweight="bold", va="top")
ax_s1c.set_title("Cohort characteristics & confounder profile", loc="left", fontsize=11.4, fontweight="bold", pad=6)

COL1_X, COL1_W = 0.08, 2.26
COL2_X, COL2_W = 2.52, 3.62
COL3_X, COL3_W = 6.32, 3.60

rbox(ax_s1c, 0.08, 6.52, 9.84, 0.80, COL["CardGray"], COL["Real"], lw=1.2, r=0.18, pad=0.06)
for xh, th in [
    (COL1_X + COL1_W / 2, "Dimension"),
    (COL2_X + COL2_W / 2, "Real (n = 37)"),
    (COL3_X + COL3_W / 2, "PBISC new (n = 29)"),
]:
    ax_s1c.text(xh, 6.92, th, ha="center", va="center", fontsize=8.9, fontweight="bold", color=COL["Real"])

rows_d = [
    ("Pandemic wave", "Wave 1 & Wave 2", "100% Wave 1 (Spring 2020)", "white", COL["CardBlue"]),
    ("Centre / batch", "NCL (21) · Cambridge (16)", "Single bulk protocol (29)", "#FBFBFB", COL["CardBlue"]),
    ("Dexamethasone", "7 pts (34 draws) exposed", "0% (pre-steroid era)", COL["CardRed"], COL["CardGreen"]),
    ("Monocyte loss", "15/60 NCL draws <2%", "0% loss (whole-PBMC RNA)", COL["CardRed"], COL["CardGreen"]),
    ("Sampling span", "Days 0–20 (median 9d)", "Days 0–30 (median 12d)", "white", COL["CardBlue"]),
    ("Abundance unit", "Cell-count % of PBMCs", "mRNA-weighted % → INT", "#FBFBFB", COL["CardTeal"]),
]
for idx_r, (dim_t, real_t, new_t, fc_r, fc_n) in enumerate(rows_d):
    yr = 5.48 - idx_r * 1.04
    rbox(ax_s1c, COL1_X, yr, COL1_W, 0.80, "#F7F7F7", PALETTE["neutral"], lw=0.85, r=0.14, pad=0.06)
    rbox(ax_s1c, COL2_X, yr, COL2_W, 0.80, fc_r, PALETTE["neutral"], lw=0.85, r=0.14, pad=0.06)
    rbox(ax_s1c, COL3_X, yr, COL3_W, 0.80, fc_n, PALETTE["neutral"], lw=0.85, r=0.14, pad=0.06)
    ax_s1c.text(COL1_X + COL1_W / 2, yr + 0.40, dim_t, ha="center", va="center", fontsize=8.0, fontweight="bold", color=COL["Real"])
    ax_s1c.text(COL2_X + COL2_W / 2, yr + 0.40, real_t, ha="center", va="center", fontsize=7.9, color=COL["Real"])
    ax_s1c.text(COL3_X + COL3_W / 2, yr + 0.40, new_t, ha="center", va="center", fontsize=7.9, fontweight="bold", color=COL["Real"])

# Panel d: Two-stage discovery & external validation bridge (enlarged bullet points)
ax_s1d = fig_s1.add_axes([0.020, 0.040, 0.962, 0.395])
ax_s1d.set_xlim(0, 100)
ax_s1d.set_ylim(0, 100)
ax_s1d.axis("off")
fig_s1.text(0.015, 0.455, "d", fontsize=18, fontweight="bold", va="top")
ax_s1d.text(
    0.2, 97.2,
    "Two-stage discovery, confounder-aware credibility filtering (Fig. 7), and external single-cell validation (Supplementary Fig. 4, N = 219)",
    fontsize=11.8, fontweight="bold", va="top", ha="left", color=COL["Real"],
)

arr_kw_s1 = dict(arrowstyle="-|>", lw=2.0, color=PALETTE["gray_2"], mutation_scale=16)
rbox(ax_s1d, 0.4, 7, 29.0, 81, "#F8FAFC", COL["Ext"], lw=1.7)
ax_s1d.text(14.9, 81.2, "1. Harmonized ESKD Screen (n = 66)", ha="center", va="center", fontsize=11.4, fontweight="bold", color=COL["Ext"])
ax_s1d.text(
    2.3, 43.0,
    "• Within-source Normal score\n  aligns CITE-seq & PBISC scales\n• 43 fine leukocyte types screened\n• 5-step pre-specified funnel:\n  43 → 16 (FDR) → 15 (same dir.)\n  → 12 (centre/dex) → 8 gained",
    ha="left", va="center", fontsize=11.2, color=COL["Real"], linespacing=1.42,
)
ax_s1d.annotate("", xy=(31.8, 47.5), xytext=(29.8, 47.5), arrowprops=arr_kw_s1)

rbox(ax_s1d, 32.2, 7, 33.6, 81, "white", COL["Real"], lw=1.7)
ax_s1d.text(49.0, 81.2, "2. Credibility-Validated Immune Signatures", ha="center", va="center", fontsize=11.4, fontweight="bold", color=COL["Real"])
rbox(ax_s1d, 33.3, 52.5, 31.4, 21.0, COL["CardGray"], COL["Real"], lw=1.1, r=1.4)
ax_s1d.text(49.0, 63.0, "Tier 1 · Concordant Anchors (Real & PBISC agree)\nNK CD56bright (↓) · B ASC dividing (↑) · T CD8 CM (↓) · CD16 mono (↓)",
            ha="center", va="center", fontsize=9.1, fontweight="bold", color=COL["Real"], linespacing=1.34)
rbox(ax_s1d, 33.3, 27.5, 31.4, 21.5, COL["CardBlue"], COL["Ext"], lw=1.2, r=1.4)
ax_s1d.text(49.0, 38.2, "Tier 2 · PBISC-Rescued Focus Associations (Fig. 7j–o)\nSeverity: B ASC IgG/IgA (↑), T CD8 CTL (↓)\nDeath: HSPC (↑, >7d), T CD8 CM/CTL (↓)",
            ha="center", va="center", fontsize=9.1, fontweight="bold", color=COL["Ext"], linespacing=1.30)
rbox(ax_s1d, 33.3, 10.5, 31.4, 13.5, COL["CardRed"], COL["CLL"], lw=1.0, r=1.2)
ax_s1d.text(49.0, 17.2, "Filtered Out (Fig. 7d, h): CD14 mono dex · NK CD56dim & CD16 mono C1",
            ha="center", va="center", fontsize=8.8, color=COL["CLL"], fontweight="bold")
ax_s1d.annotate("", xy=(68.2, 47.5), xytext=(66.2, 47.5), arrowprops=arr_kw_s1)

rbox(ax_s1d, 68.6, 7, 30.9, 81, COL["CardTeal"], COL["Shared"], lw=1.8)
ax_s1d.text(84.05, 80.2, "3. External Adult COVID-19 Validation\n(4 Cohorts · N = 219 · 23 Deaths · 2.46M Cells)",
            ha="center", va="center", fontsize=11.0, fontweight="bold", color="#165057", linespacing=1.22)
ax_s1d.text(
    70.2, 41.0,
    "• Stephenson 2021 (Nat Med, n=73, 6 deaths)\n• COMBAT 2022 (Cell Wave-1, n=75, 7 deaths)\n• van der Wijst 2021 (UCSF, n=42, 7 deaths)\n• Liu 2021 (Cell NIAID, n=29, 3 deaths)\n→ Anchors, B ASC, T CD8 CM & HSPC (>7d) hold\n→ Not-followed NK CD56dim & CD16 C1 stay null",
    ha="left", va="center", fontsize=10.8, color=COL["Real"], linespacing=1.38,
)

s1_out = OUT_DIR / "PBISC_SuppFig_S1_Cohort_Characteristics_ver9_20261118"
saved_s1 = finalize_figure(fig_s1, s1_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Supplementary Figure 1 (ver9):", [str(p) for p in saved_s1])


# =========================================================================
# 3) SUPPLEMENTARY FIGURE 2 (Renumbered from Former Supplementary Fig. 1)
# =========================================================================
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
            ax.barh(-i, v[b], left=left, color=CT_COL[ct], edgecolor="white", lw=1.1, height=0.68)
            left += v[b]
        if "C1" in lab:
            txt = f"B {v['B_pct']:.0f}% · Mono {v['Mono_pct']:.0f}% · NK {v['NK_pct']:.0f}%"
            ax.text(101.5, -i, txt, va="center", fontsize=9.5,
                    color=PALETTE["red_strong"] if src == "Real" else PALETTE["gray_2"],
                    fontweight="bold" if src == "Real" else "normal")
    ax.set_yticks(-np.arange(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=10.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 158)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("% of PBMCs (COVID+ samples, patient mean)", fontsize=10.5)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=CT_COL[c]) for _, c in BR],
              labels=[c for _, c in BR], ncol=5, loc="lower left", bbox_to_anchor=(-0.02, 1.01),
              fontsize=9.5, handlelength=1.0, columnspacing=0.8, handletextpad=0.3)


def draw_s2_cll_shift(ax, val, lo, hi, pcol, ref, xlab, xlim, xt):
    FEATS_SH = ["B_ASC_dividing", "NK CD56 bright", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "HSPC"]
    y = 0.0
    ticks = []
    for f in FEATS_SH:
        for st, k in [("Real", "Real"), ("Real + PBISC new", "Ext")]:
            a = CI[(CI["feature"] == f) & (CI["set"] == st) & (~CI["CLL_excluded"])].iloc[0]
            b = CI[(CI["feature"] == f) & (CI["set"] == st) & (CI["CLL_excluded"])].iloc[0]
            ax.plot([b[lo], b[hi]], [y, y], color=COL[k], lw=1.8, alpha=0.55)
            ax.annotate("", xy=(b[val], y), xytext=(a[val], y),
                        arrowprops=dict(arrowstyle="-|>", color=PALETTE["gray_1"], lw=1.2, mutation_scale=10, shrinkA=2, shrinkB=4))
            ax.scatter(a[val], y, s=36, marker=MK[k], c="white" if a[pcol] >= 0.05 else COL["All"], edgecolors=COL["All"], linewidths=1.4, zorder=3)
            ax.scatter(b[val], y, s=50, marker=MK[k], c="white" if b[pcol] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.6, zorder=4)
            y -= 1.0
        ticks.append(y + 1.5)
        y -= 0.55
    ax.axvline(ref, color=PALETTE["gray_1"], ls=":", lw=1.4)
    ax.set_yticks(ticks)
    ax.set_yticklabels([LBL[f] for f in FEATS_SH], fontsize=10.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(y + 0.35, 0.55)
    ax.set_xlim(*xlim)
    ax.set_xticks(xt)
    ax.set_xlabel(xlab, fontsize=10.5)


def draw_s2_nk_decomp(ax):
    measures = ["CD56bright % PBMC", "CD56bright % of NK", "Total NK % PBMC"]
    xlabels = ["CD56bright\n% PBMC", "CD56bright\n% of NK", "Total NK\n% PBMC"]
    for k_src, col_k in [("Real", "Real"), ("PBISC new", "New")]:
        for i, m_name in enumerate(measures):
            row = NKD[(NKD["source"] == k_src) & (NKD["measure"] == m_name)].iloc[0]
            r, p = row["rho"], row["p"]
            x = i + (-0.19 if k_src == "Real" else 0.19)
            ax.bar(x, r, width=0.35, color=COL[col_k] if p < 0.05 else "white", edgecolor=COL[col_k], lw=1.6)
            ax.text(x, r - 0.045 if r < 0 else r + 0.03, "*" if p < 0.05 else "ns",
                    ha="center", va="top" if r < 0 else "bottom", fontsize=9.8, fontweight="bold" if p < 0.05 else "normal")
    ax.axhline(0, color="black", lw=1.0)
    ax.set_xticks(range(3))
    ax.set_xticklabels(xlabels, fontsize=9.5)
    ax.set_ylim(-0.92, 0.28)
    ax.set_yticks([-0.8, -0.4, 0.0])
    ax.set_ylabel("Spearman ρ with severity", fontsize=10.2)
    ax.tick_params(axis="x", length=0)


def draw_s2_modality_scale(ax):
    FEATS_M = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
    x = np.arange(len(FEATS_M))
    rng_m = np.random.default_rng(19)
    for i, f in enumerate(FEATS_M):
        vr = SH27[f + "|Real_pct"].values
        vb = SH27[f + "|PBISC_pct"].values
        for r_i, b_i in zip(vr, vb):
            ax.plot([i - 0.18, i + 0.18], [np.log10(r_i + 0.02), np.log10(b_i + 0.02)], color=PALETTE["neutral"], lw=0.6, alpha=0.55, zorder=1)
        ax.scatter(i - 0.18 + rng_m.uniform(-0.04, 0.04, len(vr)), np.log10(vr + 0.02), s=22, marker="D", c=COL["Real"], edgecolors="white", linewidths=0.35, zorder=3)
        ax.scatter(i + 0.18 + rng_m.uniform(-0.04, 0.04, len(vb)), np.log10(vb + 0.02), s=24, marker="o", c=COL["New"], edgecolors="white", linewidths=0.35, zorder=3)
        ax.plot([i - 0.30, i - 0.06], [np.log10(np.median(vr) + 0.02)] * 2, color="black", lw=2.0, zorder=4)
        ax.plot([i + 0.06, i + 0.30], [np.log10(np.median(vb) + 0.02)] * 2, color="black", lw=2.0, zorder=4)

        mrow = MOD[MOD["feature"] == f].iloc[0]
        wp = mrow["wilcoxon_paired_p"]
        fold = mrow["fold_PBISC_over_Real"]
        star_s = "***" if wp < 0.001 else ("**" if wp < 0.01 else ("*" if wp < 0.05 else "ns"))
        ax.text(i, 1.82, f"{fold:.1f}×\n({star_s})", ha="center", va="top", fontsize=8.3,
                color=PALETTE["red_strong"] if wp < 0.01 else PALETTE["gray_2"],
                fontweight="bold" if wp < 0.01 else "normal")

    yt_pct = [0.0, 0.1, 0.5, 2.0, 10.0, 40.0]
    ax.set_yticks([np.log10(v + 0.02) for v in yt_pct])
    ax.set_yticklabels(["0%", "0.1%", "0.5%", "2%", "10%", "40%"], fontsize=9.5)
    ax.set_ylim(np.log10(0.015), 2.02)
    ax.set_xticks(x)
    ax.set_xticklabels([LBL[f] for f in FEATS_M], rotation=22, ha="right", fontsize=9.6)
    ax.set_ylabel("% of PBMCs (27 shared donors, log scale)", fontsize=10.0)


def draw_strip_raw_pct(ax, f, outcome, rng_local, show_ylabel=True):
    if outcome == "sev":
        for g in range(4):
            for k, dx in [("Real", -0.18), ("New", 0.18)]:
                vals = PT_RAW[(PT_RAW["src"] == k) & (PT_RAW["sev"] == g)][f].dropna().values
                x = g + dx + rng_local.uniform(-0.075, 0.075, len(vals))
                ax.scatter(x, vals, s=22, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                if len(vals):
                    ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.8, zorder=4)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=9.2)
        ax.set_xlim(-0.55, 3.55)
    else:
        for g in range(2):
            for k, dx in [("Real", -0.18), ("New", 0.18)]:
                vals = PT_RAW[(PT_RAW["src"] == k) & (PT_RAW["fatal"] == g)][f].dropna().values
                x = g + dx + rng_local.uniform(-0.075, 0.075, len(vals))
                ax.scatter(x, vals, s=24, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.4, zorder=3)
                if len(vals):
                    ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=1.8, zorder=4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Survived", "Died"], fontsize=9.2)
        ax.set_xlim(-0.52, 1.52)

    vmax = float(PT_RAW[f].max())
    ax.set_ylim(-0.04 * vmax, vmax * 1.38)
    if show_ylabel:
        ax.set_ylabel("% of PBMCs (patient mean)", fontsize=9.8)

    mrow = MOD[MOD["feature"] == f].iloc[0]
    r_real = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == "Real")].iloc[0]
    r_new = SCAN_NO[(SCAN_NO["feature"] == f) & (SCAN_NO["set"] == "PBISC new")].iloc[0]
    if outcome == "sev":
        t1 = f"Real: ρ={r_real.rho:+.2f} | PBISC: {r_new.rho:+.2f}"
        t2 = f"Raw % pool: ρ={mrow.combined_rho_raw_pct:+.2f}, {fp(mrow.combined_p_sev_raw_pct)}"
    else:
        t1 = f"Real: AUC={r_real.AUC_death:.2f} | PBISC: {r_new.AUC_death:.2f}"
        t2 = f"Raw % pool: AUC={mrow.combined_AUC_raw_pct:.2f}, {fp(mrow.combined_p_death_raw_pct)}"
    ax.text(0.03, 0.97, t1, transform=ax.transAxes, fontsize=7.9, color=PALETTE["gray_2"], va="top")
    ax.text(0.03, 0.86, t2, transform=ax.transAxes, fontsize=8.0, color=COL["Ext"], fontweight="bold", va="top")


fig_s2 = plt.figure(figsize=(15.2, 20.6), facecolor="white")
fig_s2.text(
    0.015,
    0.990,
    "Supplementary Figure 2. CLL exclusion audit, concordant anchors, NK decomposition, and within-source normal-score justification",
    fontsize=13.8,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

ax_s2a = fig_s2.add_axes([0.125, 0.838, 0.330, 0.110])
draw_s2_cll_bars(ax_s2a)
fig_s2.text(0.015, 0.962, "a", fontsize=17, fontweight="bold", va="top")

ax_s2b = fig_s2.add_axes([0.575, 0.832, 0.165, 0.122])
draw_s2_cll_shift(ax_s2b, "rho", "rho_lo", "rho_hi", "p_sev", 0.0, "ρ with severity (CLL incl. → excl.)", (-0.88, 0.88), [-0.5, 0, 0.5])
fig_s2.text(0.482, 0.962, "b", fontsize=17, fontweight="bold", va="top")

ax_s2c = fig_s2.add_axes([0.815, 0.832, 0.165, 0.122])
draw_s2_cll_shift(ax_s2c, "AUC", "AUC_lo", "AUC_hi", "p_death", 0.5, "AUC for death (CLL incl. → excl.)", (0.05, 1.02), [0.25, 0.5, 0.75, 1.0])
fig_s2.text(0.755, 0.962, "c", fontsize=17, fontweight="bold", va="top")

ax_s2_leg = fig_s2.add_axes([0.02, 0.780, 0.96, 0.020])
ax_s2_leg.set_axis_off()
hd_s2 = [
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc=COL["Real"], mec=COL["Real"], label="Real (37, CLL excluded)"),
    Line2D([0], [0], marker="o", ls="", ms=7.5, mfc=COL["New"], mec=COL["New"], label="PBISC new (29) / PBISC shared (i)"),
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["Ext"], mec=COL["Ext"], label="Real + PBISC new (66)"),
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["All"], mec=COL["All"], label="all patients, CLL included (b, c)"),
    Line2D([0], [0], marker="o", ls="", ms=7.5, mfc="white", mec=PALETTE["gray_2"], mew=1.5, label="open: P ≥ 0.05"),
]
ax_s2_leg.legend(handles=hd_s2, loc="center", ncol=5, fontsize=9.6, frameon=False, columnspacing=1.3)

rng_s2 = np.random.default_rng(17)
for idx_p, (f, let) in enumerate([("NK CD56 bright", "d"), ("B_ASC_dividing", "e"), ("T CD8 CM", "f"), ("CD16mono", "g")]):
    x_l = 0.062 + idx_p * 0.242
    ax_an = fig_s2.add_axes([x_l, 0.565, 0.182, 0.165])
    draw_strip_ns(ax_an, f, "sev", rng_s2, show_ylabel=(idx_p == 0))
    ax_an.set_title(LBL[f], loc="left", fontsize=11.2, fontweight="bold", color=COL["Real"], pad=6)
    fig_s2.text(x_l - 0.044, 0.748, let, fontsize=17, fontweight="bold", va="top")

ax_s2h = fig_s2.add_axes([0.070, 0.318, 0.240, 0.168])
draw_s2_nk_decomp(ax_s2h)
ax_s2h.set_title("NK compartment decomposition", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_s2.text(0.015, 0.505, "h", fontsize=17, fontweight="bold", va="top")

ax_s2i = fig_s2.add_axes([0.405, 0.318, 0.570, 0.168])
draw_s2_modality_scale(ax_s2i)
ax_s2i.set_title("Modality scale shift on 27 shared donors (PBISC / Real median fold & paired Wilcoxon)", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_s2.text(0.345, 0.505, "i", fontsize=17, fontweight="bold", va="top")

for idx_r, (f, out_t, let) in enumerate([
    ("B_ASC_IgG", "sev", "j"),
    ("B_ASC_IgA", "sev", "k"),
    ("T CD8 CTL", "sev", "l"),
    ("T CD8 CM", "death", "m"),
    ("T CD8 CTL", "death", "n"),
    ("HSPC", "death", "o"),
]):
    x_l = 0.052 + idx_r * 0.158
    ax_rw = fig_s2.add_axes([x_l, 0.050, 0.125, 0.175])
    draw_strip_raw_pct(ax_rw, f, out_t, rng_s2, show_ylabel=(idx_r == 0))
    ax_rw.set_title(LBL[f], loc="left", fontsize=10.8, fontweight="bold", color=COL["Real"], pad=5)
    fig_s2.text(x_l - 0.030, 0.242, let, fontsize=17, fontweight="bold", va="top")

s2_out = OUT_DIR / "PBISC_SuppFig_S2_CLL_Anchors_Modality_ver9_20261118"
saved_s2 = finalize_figure(fig_s2, s2_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Supplementary Figure 2 (ver9):", [str(p) for p in saved_s2])


# =========================================================================
# 4) SUPPLEMENTARY FIGURE 3 (Rebuilt from Former Supp Fig. 2 after moving old a, b to Main Fig. 7e, f)
# =========================================================================
def draw_s3_power_curves(ax):
    d_igg = PW[(PW["feature"] == "B_ASC_IgG") & (PW["outcome"] == "sev")].sort_values("n_real")
    d_iga = PW[(PW["feature"] == "B_ASC_IgA") & (PW["outcome"] == "sev")].sort_values("n_real")
    d_ctl = PW[(PW["feature"] == "T CD8 CTL") & (PW["outcome"] == "sev")].sort_values("n_real")
    d_hsp = PW[(PW["feature"] == "HSPC") & (PW["outcome"] == "death")].sort_values("n_real")
    d_cm = PW[(PW["feature"] == "T CD8 CM") & (PW["outcome"] == "death")].sort_values("n_real")

    ax.plot(d_igg["n_real"], d_igg["real_alone"], color=COL["Ext"], lw=1.8, ls="--", marker="D", ms=4.8, label="B ASC IgG (sev): Real")
    ax.plot(d_igg["n_real"], d_igg["plus_pbisc_new"], color=COL["Ext"], lw=2.3, ls="-", marker="s", ms=5.2, label="B ASC IgG (sev): +PBISC")

    ax.plot(d_iga["n_real"], d_iga["real_alone"], color=COL["New"], lw=1.7, ls="--", marker="D", ms=4.5, mfc="white", label="B ASC IgA (sev): Real")
    ax.plot(d_iga["n_real"], d_iga["plus_pbisc_new"], color=COL["New"], lw=2.1, ls="-", marker="s", ms=5.0, label="B ASC IgA (sev): +PBISC")

    ax.plot(d_ctl["n_real"], d_ctl["real_alone"], color=PALETTE["violet"], lw=1.7, ls="--", marker="D", ms=4.5, mfc="white", label="T CD8 CTL (sev): Real")
    ax.plot(d_ctl["n_real"], d_ctl["plus_pbisc_new"], color=PALETTE["violet"], lw=2.1, ls="-", marker="s", ms=5.0, label="T CD8 CTL (sev): +PBISC")

    ax.plot(d_hsp["n_real"], d_hsp["real_alone"], color=PALETTE["red_strong"], lw=1.8, ls="--", marker="D", ms=4.8, mfc="white", label="HSPC (death): Real")
    ax.plot(d_hsp["n_real"], d_hsp["plus_pbisc_new"], color=PALETTE["red_strong"], lw=2.3, ls="-", marker="s", ms=5.2, label="HSPC (death): +PBISC")

    ax.plot(d_cm["n_real"], d_cm["plus_pbisc_new"], color=PALETTE["teal"], lw=2.0, ls="-", marker="s", ms=4.8, label="T CD8 CM (death): +PBISC")

    shuf_avg = (d_igg["plus_shuffled"].values + d_iga["plus_shuffled"].values + d_hsp["plus_shuffled"].values) / 3.0
    ax.plot(d_igg["n_real"], shuf_avg, color=PALETTE["gray_1"], lw=1.8, ls=":", label="+PBISC shuffled (mean)")

    ax.axhline(0.8, color=PALETTE["gray_1"], ls=":", lw=1.3)
    ax.set_ylim(-0.03, 1.38)
    ax.set_yticks([0, 0.4, 0.8, 1.0])
    ax.set_xticks([10, 15, 20, 25, 30, 35])
    ax.set_xlim(8.5, 36.5)
    ax.set_xlabel("Real patients sampled (1,000 random draws per point)", fontsize=10.4)
    ax.set_ylabel("Empirical detection power (P < 0.05)", fontsize=10.4)
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=2, fontsize=8.5, handlelength=1.6, columnspacing=1.0, handletextpad=0.35, labelspacing=0.22)


def draw_s3_sensitivity_matrix(ax_sev, ax_death):
    SA_ext = SA[SA["set"] == "Real + PBISC new"]
    FEN = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN"]
    VN = {
        "0 baseline (no CLL)": "Baseline (CLL excluded)",
        "A1 adjust centre/wave": "Adjust centre / wave",
        "A2 within-lineage %": "% of own lineage",
        "A3 % of non-myeloid cells": "% of non-myeloid cells",
        "A4 drop Real samples with mono <2%": "Drop monocyte-depleted samples",
        "B  also drop C187, C190": "Also drop C187 (NK) & C190 (T) outliers",
        "C1 death = Table S1 (C101, C147 alive)": "Death: acute COVID (dir.+indir., S1)",
        "C2 death = direct only": "Death: direct COVID-19 only",
        "D1 Real steroid-naive samples only": "Real: steroid-naive samples only",
        "D2 adjust patient CD14mono_dex %": "Adjust CD14 mono dex (steroid proxy)",
        "E  samples <= day 20 (both sources)": "Samples ≤ day 20 after first swab",
        "F1 first sample per patient": "First sample per patient",
        "F2 adjust age + sex": "Adjust age and sex",
        "G  combined: lineage %, no C187/C190, adj centre/wave+age+sex": "Combined: lineage %, no outliers, adj.",
    }
    specs = [
        (ax_sev, "rho", "p_sev", "Severity: ρ (Real + PBISC new)", mpl.colors.Normalize(-0.7, 0.7), True),
        (ax_death, "AUC", "p_death", "Death: AUC (Real + PBISC new)", mpl.colors.Normalize(0.15, 0.85), False),
    ]
    for ax, val, pcol, ttl, nrm, show_y in specs:
        M = SA_ext.pivot_table(index="variant", columns="feature", values=val).reindex(index=list(VN), columns=FEN)
        Pm = SA_ext.pivot_table(index="variant", columns="feature", values=pcol).reindex(index=list(VN), columns=FEN)
        ax.imshow(M.values, cmap="RdBu_r", norm=nrm, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M.values[i, j]
                if np.isnan(v):
                    ax.text(j, i, "–", ha="center", va="center", fontsize=9.2, color=PALETTE["gray_1"])
                    continue
                sig = Pm.values[i, j] < 0.05
                ax.text(
                    j, i, f"{v:.2f}", ha="center", va="center", fontsize=9.0,
                    fontweight="bold" if sig else "normal",
                    color=("white" if abs(nrm(v) - 0.5) > 0.32 else "black") if sig else PALETTE["gray_1"],
                )
        ax.axvline(3.5, color="black", lw=1.6, ls="-")
        ax.set_xticks(range(len(FEN)))
        ax.set_xticklabels([LBL[f] for f in FEN], rotation=32, ha="left", fontsize=9.8)
        ax.xaxis.tick_top()
        ax.set_yticks(range(len(VN)))
        ax.set_yticklabels(list(VN.values()) if show_y else [], fontsize=9.8)
        ax.tick_params(length=0)
        for j, t in enumerate(ax.get_xticklabels()):
            f = FEN[j]
            t.set_fontweight("bold")
            t.set_color(COL["Real"] if j < 4 else (PALETTE["red_strong"] if f in {"HSPC", "CD14mono_IFN"} else COL["Ext"]))
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(FEN)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(VN)), minor=True)
        ax.grid(which="minor", color="white", lw=1.3)
        ax.tick_params(which="minor", length=0)
        ax.set_title(ttl, loc="left", fontsize=11.4, fontweight="bold", pad=28)


fig_s3 = plt.figure(figsize=(15.2, 13.4), facecolor="white")
fig_s3.text(
    0.015,
    0.986,
    "Supplementary Figure 3. Continuous subsampling power curves, Step-4 monocyte audit, and 14-variant confounder sensitivity matrix",
    fontsize=13.8,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

# Row 1 of S3: a (Continuous subsampling power curves, spacious!) & b (CD14 mono IFN Step-4 audit)
ax_s3a = fig_s3.add_axes([0.075, 0.595, 0.545, 0.320])
draw_s3_power_curves(ax_s3a)
ax_s3a.set_title("Subsampling power curves across Real sample sizes (n = 10–35)", loc="left", fontsize=11.4, fontweight="bold", pad=8)
fig_s3.text(0.015, 0.942, "a", fontsize=17.5, fontweight="bold", va="top")

ax_s3b = fig_s3.add_axes([0.725, 0.595, 0.240, 0.320])
draw_strip_ns(ax_s3b, "CD14mono_IFN", "death", np.random.default_rng(23), show_ylabel=True)
ax_s3b.text(0.03, 0.68, "1st sample AUC=0.78 (P=0.002)\nCentre-adj P=0.071 (Step 4 filter)",
            transform=ax_s3b.transAxes, fontsize=8.8, color=PALETTE["red_strong"], fontweight="bold", va="top")
ax_s3b.set_title("CD14 mono IFN", loc="left", fontsize=11.4, fontweight="bold", color=COL["Real"], pad=8)
fig_s3.text(0.658, 0.942, "b", fontsize=17.5, fontweight="bold", va="top")

# Row 2 of S3: c (14-variant confounder & denominator sensitivity matrix)
ax_s3c_s = fig_s3.add_axes([0.228, 0.040, 0.365, 0.385])
ax_s3c_d = fig_s3.add_axes([0.612, 0.040, 0.365, 0.385])
draw_s3_sensitivity_matrix(ax_s3c_s, ax_s3c_d)
fig_s3.text(0.015, 0.505, "c", fontsize=17.5, fontweight="bold", va="top")

s3_out = OUT_DIR / "PBISC_SuppFig_S3_Power_IFN_Sensitivity_ver9_20261118"
saved_s3 = finalize_figure(fig_s3, s3_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Supplementary Figure 3 (ver9):", [str(p) for p in saved_s3])


# =========================================================================
# 5) SUPPLEMENTARY FIGURE 4 (Renumbered from Former Supp Fig. 3; updated Fig. 7h refs)
# =========================================================================
pt_raw_eskd = PT_RAW.copy()
b_asc_cols_eskd = [c for c in ["B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM", "B_ASC_other"] if c in pt_raw_eskd.columns]
pt_raw_eskd["B_ASC_total"] = pt_raw_eskd[b_asc_cols_eskd].sum(axis=1)


def rint_s(x: pd.Series) -> pd.Series:
    r = x.rank()
    return pd.Series(stats.norm.ppf((r - 0.5) / x.notna().sum()), index=x.index)


ESKD_MATCH_FEATS = [
    "NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono|lin",
    "B_ASC_total", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM", "T CD8 CTL", "HSPC", "NK", "CD16mono_C1",
]
rng_e = np.random.default_rng(20261118)
eskd_matched = {}
R_e = pt_raw_eskd[pt_raw_eskd["src"] == "Real"]
N_e = pt_raw_eskd[pt_raw_eskd["src"] == "New"]
for fk in ESKD_MATCH_FEATS:
    comb = pd.concat([
        pd.DataFrame({"y": rint_s(R_e[fk]), "sev": R_e["sev"]}),
        pd.DataFrame({"y": rint_s(N_e[fk]), "sev": N_e["sev"]}),
    ]).dropna()
    r_v, p_v = stats.spearmanr(comb["sev"], comb["y"])
    yv, sv = comb["y"].values, comb["sev"].values
    n_c = len(comb)
    bs_e = [stats.spearmanr(sv[ix := rng_e.integers(0, n_c, n_c)], yv[ix])[0] for _ in range(2000)]
    eskd_matched[fk] = dict(rho=float(r_v), lo=float(np.percentile(bs_e, 2.5)), hi=float(np.percentile(bs_e, 97.5)), p=float(p_v))

fig_s4 = plt.figure(figsize=(15.4, 14.4), facecolor="white")
fig_s4.text(
    0.015,
    0.988,
    "Supplementary Figure 4. External single-cell replication of Figure 7 findings across four adult COVID-19 PBMC cohorts (N = 219 patients)",
    fontsize=13.4,
    fontweight="bold",
    va="top",
    ha="left",
    color=PALETTE["gray_3"],
)

ax_4a = fig_s4.add_axes([0.155, 0.535, 0.285, 0.375])
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
y_pos_4 = []
yc_4 = 0.0
for idx, (pub_f, disp, eskd_f, grp) in enumerate(ORDER_4A):
    if idx in (4, 7):
        ax_4a.axhline(yc_4 + 0.52, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=1)
        yc_4 -= 0.35
    y_pos_4.append(yc_4)

    em = eskd_matched[eskd_f]
    mr = PUB_META[PUB_META["feature"] == pub_f].iloc[0]
    sub_c = PUB_COH[PUB_COH["feature"] == pub_f]

    ax_4a.plot([em["lo"], em["hi"]], [yc_4 + 0.18, yc_4 + 0.18], color=PALETTE["blue_main"], lw=1.7, zorder=2)
    ax_4a.scatter(em["rho"], yc_4 + 0.18, s=46, marker="s", c=PALETTE["blue_main"], edgecolors=PALETTE["blue_main"], zorder=3)

    p_col = PALETTE["red_strong"] if grp == "hspc" else (PALETTE["gray_1"] if grp == "excl" else PALETTE["teal"])
    ax_4a.scatter(sub_c["rho"], np.repeat(yc_4 - 0.18, len(sub_c)), s=24, marker="o", c=PALETTE["neutral"], edgecolors=PALETTE["gray_1"], linewidths=0.8, zorder=2)
    ax_4a.plot([mr["rho_lo"], mr["rho_hi"]], [yc_4 - 0.18, yc_4 - 0.18], color=p_col, lw=2.0, zorder=3)
    ax_4a.scatter(mr["rho_meta"], yc_4 - 0.18, s=54, marker="D",
                  c="white" if mr["p_rho_meta"] >= 0.05 else p_col, edgecolors=p_col, linewidths=1.6, zorder=4)

    p_txt = f"{mr['p_rho_meta']:.1e}".replace("e-0", "e-") if mr["p_rho_meta"] < 0.001 else (f"{mr['p_rho_meta']:.3f}" if mr["p_rho_meta"] < 0.1 else f"{mr['p_rho_meta']:.2f}")
    ax_4a.text(0.83, yc_4 - 0.02, f"{mr['rho_meta']:+.2f} (P={p_txt}) · k={int(mr['k_cohorts'])}, N={int(mr['N_patients'])}",
               va="center", ha="left", fontsize=7.8,
               color=p_col if mr["p_rho_meta"] < 0.05 else PALETTE["gray_1"],
               fontweight="bold" if mr["p_rho_meta"] < 0.05 else "normal")
    yc_4 -= 1.05

ax_4a.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
ax_4a.set_yticks(y_pos_4)
ax_4a.set_yticklabels([x[1] for x in ORDER_4A], fontsize=9.5)
for idx_t, t in enumerate(ax_4a.get_yticklabels()):
    grp = ORDER_4A[idx_t][3]
    t.set_fontweight("bold" if grp != "excl" else "normal")
    t.set_color(PALETTE["gray_3"] if grp == "anchor" else (PALETTE["red_strong"] if grp == "hspc" else (PALETTE["gray_1"] if grp == "excl" else PALETTE["blue_main"])))
ax_4a.tick_params(axis="y", length=0)
ax_4a.spines["left"].set_visible(False)
ax_4a.set_xlim(-0.82, 0.80)
ax_4a.set_xticks([-0.5, 0.0, 0.5])
ax_4a.set_ylim(y_pos_4[-1] - 0.55, 0.55)
ax_4a.set_xlabel("Spearman ρ with WHO severity", fontsize=10.2)
ax_4a.text(0.83, 0.50, "Public Meta ρ (P) · cohorts", fontsize=8.1, fontweight="bold", color=PALETTE["gray_2"], va="bottom")
ax_4a.set_title("Matched severity comparison", loc="left", fontsize=11.3, fontweight="bold", pad=8)
fig_s4.text(0.015, 0.930, "a", fontsize=18, fontweight="bold", va="top")

ax_4b1 = fig_s4.add_axes([0.670, 0.575, 0.132, 0.335])
ax_4b2 = fig_s4.add_axes([0.830, 0.575, 0.140, 0.335])
ISO_ORDER = [
    ("B_ASC_dividing", "Plasmablast (dividing)", "B_ASC_dividing"),
    ("B_ASC_IgG", "Plasma cell IgG", "B_ASC_IgG"),
    ("B_ASC_IgA", "Plasma cell IgA", "B_ASC_IgA"),
    ("B_ASC_IgM", "Plasma cell IgM", "B_ASC_IgM"),
    ("B_ASC_total", "B ASC total", "B_ASC_total"),
]
yy_4b = np.arange(len(ISO_ORDER))[::-1]
ax_4b1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
ax_4b2.axvline(1.0, color=PALETTE["gray_1"], ls=":", lw=1.3)

for yi, (fk, flab, efk) in zip(yy_4b, ISO_ORDER):
    r = PUB_ST_B[PUB_ST_B["feature"] == fk].iloc[0]
    em = eskd_matched[efk]
    col_b = PALETTE["blue_main"] if fk in ("B_ASC_IgG", "B_ASC_IgA") else PALETTE["teal"]

    ax_4b1.plot([em["lo"], em["hi"]], [yi + 0.16, yi + 0.16], color=PALETTE["blue_main"], lw=1.6)
    ax_4b1.scatter(em["rho"], yi + 0.16, s=42, marker="s", c=PALETTE["blue_main"], edgecolors=PALETTE["blue_main"], zorder=3)

    ax_4b1.plot([r["rho_lo"], r["rho_hi"]], [yi - 0.16, yi - 0.16], color=PALETTE["teal"], lw=1.8)
    ax_4b1.scatter(r["rho"], yi - 0.16, s=48, marker="D",
                   c="white" if r["p_sev"] >= 0.05 else PALETTE["teal"], edgecolors=PALETTE["teal"], linewidths=1.5, zorder=4)
    ax_4b1.text(0.98, yi - 0.16, f"P={r['p_sev']:.2f}" if r["p_sev"] >= 0.05 else f"P={r['p_sev']:.3f}",
                va="center", ha="right", fontsize=7.6, color=PALETTE["gray_1"] if r["p_sev"] >= 0.05 else PALETTE["teal"],
                fontweight="bold" if r["p_sev"] < 0.05 else "normal")

    ax_4b2.plot([r["OR_lo"], r["OR_hi"]], [yi, yi], color=col_b, lw=2.0)
    ax_4b2.scatter(r["OR_glm"], yi, s=52, marker="D", c=col_b, edgecolors="white", linewidths=0.7, zorder=4)
    ax_4b2.text(r["OR_hi"] + 0.08, yi, f"{r['OR_glm']:.2f} (P={r['p_glm']:.3f})", va="center", ha="left",
                fontsize=7.8, color=col_b, fontweight="bold")

ax_4b1.set_yticks(yy_4b)
ax_4b1.set_yticklabels([x[1] for x in ISO_ORDER], fontsize=9.4)
for idx_t, t in enumerate(ax_4b1.get_yticklabels()):
    if ISO_ORDER[idx_t][0] in ("B_ASC_IgG", "B_ASC_IgA"):
        t.set_fontweight("bold")
        t.set_color(PALETTE["blue_main"])
ax_4b2.set_yticks(yy_4b)
ax_4b2.set_yticklabels([])
for axx in (ax_4b1, ax_4b2):
    axx.tick_params(axis="y", length=0)
    axx.spines["left"].set_visible(False)
    axx.set_ylim(-0.55, len(ISO_ORDER) - 0.35)
ax_4b1.set_xlim(-0.20, 1.02)
ax_4b1.set_xticks([0.0, 0.4, 0.8])
ax_4b1.set_xlabel("Spearman ρ (% PBMC)", fontsize=9.8)
ax_4b1.set_title("Stephenson 2021 B ASC isotypes (n = 73): rank ρ vs. count GLM", loc="left", fontsize=11.1, fontweight="bold", pad=8)
ax_4b2.set_xlim(0.85, 4.95)
ax_4b2.set_xticks([1.0, 2.0, 3.0, 4.0])
ax_4b2.set_xlabel("Adj. count-GLM OR / step", fontsize=9.8)
fig_s4.text(0.545, 0.930, "b", fontsize=18, fontweight="bold", va="top")

ax_4leg = fig_s4.add_axes([0.02, 0.462, 0.96, 0.026])
ax_4leg.set_axis_off()
hd_4 = [
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=PALETTE["blue_main"], mec=PALETTE["blue_main"], label="ESKD Real + PBISC new (n = 66, 13 deaths)"),
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc=PALETTE["teal"], mec=PALETTE["teal"], label="Public COVID-19 Meta-analysis / Stephenson 2021 (N = 219 across 4 cohorts)"),
    Line2D([0], [0], marker="o", ls="", ms=6.2, mfc=PALETTE["neutral"], mec=PALETTE["gray_1"], label="Individual public cohort (a)"),
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc="white", mec=PALETTE["gray_2"], mew=1.6, label="open: P ≥ 0.05"),
]
ax_4leg.legend(handles=hd_4, loc="center", ncol=4, fontsize=9.5, frameon=False, columnspacing=1.3)

ax_4c = fig_s4.add_axes([0.155, 0.095, 0.285, 0.320])
D_ITEMS_4 = [
    ("HSPC", "HSPC (k=3)", PALETTE["red_strong"]),
    ("T_CD8_CM", "T CD8 CM / EM* (k=2)", PALETTE["blue_main"]),
    ("T_CD8_CTL", "T CD8 CTL / Mem* (k=4)", PALETTE["blue_main"]),
    ("B_ASC_total", "B ASC total (k=4)", PALETTE["gray_2"]),
    ("NK_CD56dim", "NK CD56dim (Fig. 7h excl., k=4)", PALETTE["gray_1"]),
    ("CD16_mono_C1", "CD16 mono C1 (Fig. 7h excl., k=1)", PALETTE["gray_1"]),
]
yy_4c = np.arange(len(D_ITEMS_4))[::-1]
ax_4c.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
for yi, (pub_f, disp, col_c) in zip(yy_4c, D_ITEMS_4):
    mr = PUB_META[PUB_META["feature"] == pub_f].iloc[0]
    ax_4c.plot([mr["AUC_all_lo"], mr["AUC_all_hi"]], [yi + 0.17, yi + 0.17], color=col_c, lw=1.9)
    ax_4c.scatter(mr["AUC_death_all_meta"], yi + 0.17, s=50, marker="D",
                  c="white" if mr["p_death_all_meta"] >= 0.05 else col_c, edgecolors=col_c, linewidths=1.5, zorder=4)
    ax_4c.plot([mr["AUC_WHO6_lo"], mr["AUC_WHO6_hi"]], [yi - 0.17, yi - 0.17], color=col_c, lw=1.5, ls="--", alpha=0.8)
    ax_4c.scatter(mr["AUC_death_WHO6_meta"], yi - 0.17, s=46, marker="^",
                  c="white" if mr["p_death_WHO6_meta"] >= 0.05 else col_c, edgecolors=col_c, linewidths=1.4, zorder=4)
    p1_t = f"{mr['p_death_all_meta']:.3f}" if mr["p_death_all_meta"] < 0.1 else f"{mr['p_death_all_meta']:.2f}"
    p2_t = f"{mr['p_death_WHO6_meta']:.2f}"
    ax_4c.text(0.88, yi + 0.17, f"All: {mr['AUC_death_all_meta']:.2f} (P={p1_t})", va="center", fontsize=8.0,
               color=col_c if mr["p_death_all_meta"] < 0.05 else PALETTE["gray_2"],
               fontweight="bold" if mr["p_death_all_meta"] < 0.05 else "normal")
    ax_4c.text(0.88, yi - 0.17, f"WHO≥6: {mr['AUC_death_WHO6_meta']:.2f} (P={p2_t})", va="center", fontsize=7.8, color=PALETTE["gray_1"])

ax_4c.set_yticks(yy_4c)
ax_4c.set_yticklabels([x[1] for x in D_ITEMS_4], fontsize=9.2)
for idx_t, t in enumerate(ax_4c.get_yticklabels()):
    if idx_t < 3:
        t.set_fontweight("bold")
        t.set_color(D_ITEMS_4[idx_t][2])
ax_4c.tick_params(axis="y", length=0)
ax_4c.spines["left"].set_visible(False)
ax_4c.set_xlim(0.12, 0.86)
ax_4c.set_xticks([0.25, 0.50, 0.75])
ax_4c.set_ylim(-0.55, len(D_ITEMS_4) - 0.35)
ax_4c.set_xlabel("Public Meta Death AUC (diamond: All  ·  triangle: WHO-CPS ≥ 6)", fontsize=9.5)
ax_4c.set_title("Mortality: univariate vs. within-severe", loc="left", fontsize=11.3, fontweight="bold", pad=8)
fig_s4.text(0.015, 0.435, "c", fontsize=18, fontweight="bold", va="top")

ax_4d = fig_s4.add_axes([0.665, 0.095, 0.215, 0.320])
COH_4D = ["Stephenson 2021", "van der Wijst 2021", "COMBAT 2022", "Public Meta (3 cohorts)"]
yy_4d = np.arange(len(COH_4D))[::-1]
ax_4d.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)
for yi, cn in zip(yy_4d, COH_4D):
    r1 = PUB_HW[(PUB_HW["cohort"] == cn) & (PUB_HW["window"] == "First acute sample")].iloc[0]
    r2 = PUB_HW[(PUB_HW["cohort"] == cn) & (PUB_HW["window"] == "Samples > day 7")].iloc[0]
    is_meta = "Meta" in cn
    ax_4d.plot([r1["AUC_all"], r2["AUC_all"]], [yi + 0.15, yi - 0.15], color=PALETTE["neutral"], lw=1.5, zorder=2)
    ax_4d.scatter(r1["AUC_all"], yi + 0.15, s=62 if is_meta else 46, marker="D",
                  c="white" if r1["p_all"] >= 0.05 else PALETTE["gray_3"], edgecolors=PALETTE["gray_3"], linewidths=1.5, zorder=3)
    ax_4d.scatter(r2["AUC_all"], yi - 0.15, s=66 if is_meta else 50, marker="s",
                  c="white" if r2["p_all"] >= 0.05 else PALETTE["red_strong"], edgecolors=PALETTE["red_strong"], linewidths=1.5, zorder=4)
    ax_4d.text(0.95, yi + 0.15, f"1st: {r1['AUC_all']:.2f} (P={r1['p_all']:.3f})", va="center", fontsize=8.0,
               color=PALETTE["gray_3"] if r1["p_all"] < 0.05 else PALETTE["gray_1"],
               fontweight="bold" if r1["p_all"] < 0.05 else "normal")
    ax_4d.text(0.95, yi - 0.15, f">7d: {r2['AUC_all']:.2f} (P={r2['p_all']:.3f})", va="center", fontsize=8.0,
               color=PALETTE["red_strong"] if r2["p_all"] < 0.05 else PALETTE["gray_1"],
               fontweight="bold" if r2["p_all"] < 0.05 else "normal")

ax_4d.set_yticks(yy_4d)
ax_4d.set_yticklabels(COH_4D, fontsize=9.5)
ax_4d.get_yticklabels()[-1].set_fontweight("bold")
ax_4d.get_yticklabels()[-1].set_color(PALETTE["red_strong"])
ax_4d.tick_params(axis="y", length=0)
ax_4d.spines["left"].set_visible(False)
ax_4d.set_xlim(0.42, 0.93)
ax_4d.set_xticks([0.50, 0.70, 0.90])
ax_4d.set_ylim(-0.55, len(COH_4D) - 0.35)
ax_4d.set_xlabel("HSPC Death AUC (diamond: 1st  ·  square: >7d)", fontsize=9.5)
ax_4d.set_title("HSPC mortality across windows", loc="left", fontsize=11.3, fontweight="bold", pad=8)
fig_s4.text(0.545, 0.435, "d", fontsize=18, fontweight="bold", va="top")

fig_s4.text(
    0.015,
    0.018,
    "*Cohort annotation proxies: B ASC dividing = Plasmablast (Stephenson), PB_Prolif (van der Wijst); T CD8 CM = CD8.EM (Stephenson), CD8.TCM+TCM.CCL5 (COMBAT); "
    "T CD8 CTL = CD8.TE (Stephenson), CD8.TEMRA+TEFF+TEM (COMBAT), T8_Mem (van der Wijst), CD8 memory T (Liu). All available cohorts per feature are pooled.",
    fontsize=8.3,
    color=PALETTE["gray_2"],
    va="bottom",
    ha="left",
)

s4_out = OUT_DIR / "PBISC_SuppFig_S4_Public_COVID_Replication_ver9_20261118"
saved_s4 = finalize_figure(fig_s4, s4_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Supplementary Figure 4 (ver9):", [str(p) for p in saved_s4])
