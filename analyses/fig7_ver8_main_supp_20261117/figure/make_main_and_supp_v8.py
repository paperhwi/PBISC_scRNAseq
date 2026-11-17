"""Render Harmonized Main Figure 7 (ver8) + Supplementary Figures S1 & S2 + Standalone Subpanels.

Strictly adheres to F:\\PBISC\\AGENTS.md (.scientific-figure-making house style) and transfers
14_Figure7_ver7_CLL_excluded visual identity onto all newly integrated panels:
  - Main Figure 7 (a-m):
      a     : CLL molecular screen (Bulk PBMC ROR1 vs FMOD across bulk donors; user's requested panel)
      b     : Pre-specified 5-step selection funnel (43 -> 16 -> 15 -> 12 -> 8)
      c     : Associations gained by adding PBISC (-log10 q Gain Map, showing Focus hits & Concordant Anchors)
      d     : Credibility scorecard of the 8 PBISC-dependent candidates (Followed vs Not followed)
      e     : Cohort effect sizes with bootstrap 95% CIs (Real, PBISC new, Combined) for Anchors & Focus populations
      f-k   : Patient-level distributions (within-source Normal score); TOP TITLE = ONLY CELL POPULATION NAME
              (f: B ASC IgG, g: B ASC IgA, h: T CD8 CTL [severity]; i: T CD8 CM, j: T CD8 CTL, k: HSPC [death])
      l     : HSPC mortality robustness across 8 definitions / strata (Table S1, Direct COVID, WHO 2-3, >7d, etc.)
      m     : Statistical power at 20 Real patients (Real, +PBISC new, +PBISC shuffled) & n80 sample-size reduction
  - Supplementary Figure S1 (a-o):
      CLL PBMC distortion (a), CLL before/after effect shift (b, c), 4 Concordant Anchor strips (d-g),
      NK compartment decomposition (h), 27-shared-donor modality scale proof for Normal score vs % PBMC (i),
      and Raw % of PBMCs patient-level distributions for all 6 Main focus associations (j-o).
  - Supplementary Figure S2 (a-e):
      Full 6-row Fine-label association screen with H/S flags (a), n80 sample-size dumbbell (b),
      Continuous subsampling power curves n=10..35 (c), CD14 mono IFN Step-4 filter audit (d),
      and 14-variant dataset & confounder sensitivity matrix (e, moved from old Panel n).
"""
import pathlib
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
V8_DIR = pathlib.Path(r"F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\15_Figure7_ver8_Main_and_Supp")
sys.path.insert(0, str(V8_DIR))
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
DATA_DIR = V8_DIR / "data"
PANEL_DIR = V8_DIR / "panels"
for sub in ("main", "supp_s1", "supp_s2"):
    (PANEL_DIR / sub).mkdir(parents=True, exist_ok=True)

apply_publication_style(FigureStyle(font_size=13, axes_linewidth=1.8))

CLL = {"C116", "C139"}
COL = dict(
    Real=PALETTE["gray_3"],
    New=PALETTE["blue_secondary"],
    Ext=PALETTE["blue_main"],
    Shuf=PALETTE["neutral"],
    All=PALETTE["red_2"],
    CLL=PALETTE["red_strong"],
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


# Load all precomputed v8 tables
S_meta = pd.read_csv(S_PATH)
CLL_BULK = pd.read_csv(DATA_DIR / "cll_molecular_screen_bulk_donors.csv")
SCAN = pd.read_csv(DATA_DIR / "cll_scan_43fine_v8.csv")
SCAN_NO = SCAN[SCAN["CLL_excluded"] == True].copy()
SEL = pd.read_csv(DATA_DIR / "selection_funnel_43fine.csv")
SC = pd.read_csv(DATA_DIR / "credibility_scorecard_v8.csv")
CI = pd.read_csv(DATA_DIR / "cohort_effect_sizes_ci_v8.csv")
HS = pd.read_csv(DATA_DIR / "hspc_death_analyses_v8.csv")
PW = pd.read_csv(DATA_DIR / "power_curves_v8.csv")
MOD = pd.read_csv(DATA_DIR / "modality_scale_shared27_and_raw_pct_v8.csv")
SH27 = pd.read_csv(DATA_DIR / "shared27_paired_pct_v8.csv")
PT_NS = pd.read_csv(DATA_DIR / "patient_normal_scores_noCLL_v8.csv", index_col=0)
PT_RAW = pd.read_csv(DATA_DIR / "patient_raw_pct_noCLL_v8.csv", index_col=0)
SA = pd.read_csv(DATA_DIR / "sensitivity_14variants_all_v8.csv")
NKD = pd.read_csv(DATA_DIR / "nk_decomposition_v8.csv")


# =========================================================================
# MAIN FIGURE 7 PANEL FUNCTIONS (a - m)
# =========================================================================
def draw_main_a_cll_screen(ax):
    """Bulk PBMC ROR1 vs FMOD CLL molecular screen transferred to ver7 house style."""
    df_b = CLL_BULK.set_index("patient_id")
    # Plot COVID-negative ESKD bulk controls lightly in background, COVID+ Real/Shared & PBISC-only prominently
    neg = df_b[~df_b["covid_positive_scope"] & (df_b["cohort_group"] != "Excluded CLL (C116, C139)")]
    pos_real = df_b[df_b["covid_positive_scope"] & (df_b["cohort_group"] == "Real / Shared donor")]
    pos_pbisc = df_b[df_b["covid_positive_scope"] & (df_b["cohort_group"] == "PBISC-only donor")]
    cll_pts = df_b[df_b["cohort_group"] == "Excluded CLL (C116, C139)"]

    ax.scatter(neg["ROR1"], neg["FMOD"], s=26, marker="o", c="white", edgecolors=PALETTE["neutral"], linewidths=0.9, alpha=0.8, zorder=1,
               label=f"COVID− ESKD bulk ({len(neg)})")
    ax.scatter(pos_real["ROR1"], pos_real["FMOD"], s=38, marker="D", c=PALETTE["neutral"], edgecolors=COL["Real"], linewidths=1.0, alpha=0.92, zorder=2,
               label=f"Real / Shared COVID+ ({len(pos_real)})")
    ax.scatter(pos_pbisc["ROR1"], pos_pbisc["FMOD"], s=44, marker="^", c=COL["New"], edgecolors="black", linewidths=0.8, alpha=0.92, zorder=3,
               label=f"PBISC-only COVID+ ({len(pos_pbisc)})")
    ax.scatter(cll_pts["ROR1"], cll_pts["FMOD"], s=115, marker="D", c=COL["CLL"], edgecolors="black", linewidths=1.2, zorder=5,
               label="Benign CLL (C116, C139)")

    ax.annotate("C116\n(CLL)", xy=(df_b.loc["C116", "ROR1"], df_b.loc["C116", "FMOD"]),
                xytext=(3.35, 4.35), fontsize=9.5, fontweight="bold", color=COL["CLL"],
                arrowprops=dict(arrowstyle="-", color=COL["CLL"], lw=1.3))
    ax.annotate("C139\n(CLL)", xy=(df_b.loc["C139", "ROR1"], df_b.loc["C139", "FMOD"]),
                xytext=(3.60, 2.10), fontsize=9.5, fontweight="bold", color=COL["CLL"],
                arrowprops=dict(arrowstyle="-", color=COL["CLL"], lw=1.3))

    ax.axhline(1.0, color=PALETTE["gray_1"], ls="--", lw=1.1, alpha=0.75)
    ax.axvline(2.5, color=PALETTE["gray_1"], ls="--", lw=1.1, alpha=0.75)
    ax.set_xlim(-0.3, 6.2)
    ax.set_ylim(-0.35, 6.45)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_yticks([0, 2, 4, 6])
    ax.set_xlabel("Bulk PBMC ROR1 log2(CPM+1)", fontsize=10.8)
    ax.set_ylabel("Bulk PBMC FMOD log2(CPM+1)", fontsize=10.8)
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.99), fontsize=8.3, frameon=False,
              handletextpad=0.25, labelspacing=0.20, borderpad=0.1)


def draw_main_b_funnel(ax):
    """Pre-specified 5-step selection funnel on the 43 fine cell types."""
    steps = [
        ("Tested fine types", len(SEL), "≥30% Real detection"),
        ("Associated\n(Real + PBISC new)", int((SEL["S2_sev"] | SEL["S2_death"]).sum()), "sev q<0.05 or death q<0.10"),
        ("Same direction\nin both sources", int((SEL["S3_sev"] | SEL["S3_death"]).sum()), "− NKT (Real ρ ≈ 0)"),
        ("Not centre- or\ntreatment-driven", int((SEL["S4_sev"] | SEL["S4_death"]).sum()), "− CD14 mono dex/IFN, MAIT"),
        ("PBISC-dependent\n(not in Real alone)", int((SEL["S5_sev"] | SEL["S5_death"]).sum()), "− 4 Real-sig. anchors"),
    ]
    y = np.arange(len(steps))[::-1]
    fills = ["#DCE6F2", "#B3CBE6", "#7FA9D6", COL["New"], COL["Ext"]]
    for yi, (lab, cnt, note), fc in zip(y, steps, fills):
        ax.barh(yi, cnt, height=0.66, color=fc, edgecolor=COL["Real"], linewidth=1.4, zorder=2)
        ax.text(cnt + 1.2, yi, f"{cnt}", va="center", ha="left", fontsize=10.8, fontweight="bold", color=COL["Real"])
        offset_x = 7.4 if cnt >= 10 else 5.4
        ax.text(cnt + offset_x, yi, note, va="center", ha="left", fontsize=8.4,
                color=PALETTE["red_strong"] if note.startswith("−") else PALETTE["gray_2"])
    ax.set_yticks(y)
    ax.set_yticklabels([s[0] for s in steps], fontsize=9.8)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 82)
    ax.set_xticks([0, 20, 40])
    ax.set_xlabel("Fine-grained cell types", fontsize=10.8)


def draw_main_c_gain_map(ax):
    """-log10 q in Real alone vs Real + PBISC new across all 43 fine cell types."""
    cap = 4.15
    nl = lambda q: np.minimum(-np.log10(np.maximum(q, 1e-6)), cap)
    q05 = -np.log10(0.05)
    q10 = 1.0

    # Discovery region (Real q >= 0.05, Combined q < 0.05)
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
                fc, ec, sz, zo = COL["Ext"], "black", 68, 5
            elif is_anc:
                fc, ec, sz, zo = COL["Real"], "white", 56, 4
            elif is_s5:
                fc, ec, sz, zo = PALETTE["gray_1"], "none", 36, 3
            else:
                fc, ec, sz, zo = PALETTE["neutral"], "none", 22, 2

            ax.scatter(xv, yv, s=sz, marker=mk, c=fc, edgecolors=ec, linewidths=1.0 if ec != "none" else 0.0, zorder=zo)

    # Non-overlapping curated callouts for Focus (PBISC-gained) & Concordant Anchors
    callouts = [
        ("B_ASC_IgG", "sev", "q_sev_Real", "q_sev_Ext", "B ASC IgG (sev)", (0.08, 3.68), COL["Ext"]),
        ("B_ASC_IgA", "sev", "q_sev_Real", "q_sev_Ext", "B ASC IgA (sev)", (0.12, 3.12), COL["Ext"]),
        ("HSPC", "death", "q_death_Real", "q_death_Ext", "HSPC (death)", (0.08, 2.56), PALETTE["red_strong"]),
        ("T CD8 CTL", "sev", "q_sev_Real", "q_sev_Ext", "T CD8 CTL (sev)", (0.08, 1.96), COL["Ext"]),
        ("T CD8 CTL", "death", "q_death_Real", "q_death_Ext", "T CD8 CTL (death)", (0.20, 0.42), COL["Ext"]),
        ("T CD8 CM", "death", "q_death_Real", "q_death_Ext", "T CD8 CM (death)", (1.42, 0.72), COL["Ext"]),
        # Concordant anchors in upper-right quadrant
        ("NK CD56 bright", "sev", "q_sev_Real", "q_sev_Ext", "NK CD56bright", (2.38, 3.88), COL["Real"]),
        ("T CD8 CM", "sev", "q_sev_Real", "q_sev_Ext", "T CD8 CM (sev)", (2.52, 3.00), COL["Real"]),
        ("B_ASC_dividing", "sev", "q_sev_Real", "q_sev_Ext", "B ASC dividing", (2.50, 2.45), COL["Real"]),
        ("CD16mono", "sev", "q_sev_Real", "q_sev_Ext", "CD16 mono", (2.35, 1.72), COL["Real"]),
    ]
    sel_ix = SEL.set_index("feature")
    for f, _, qr, qe, txt, xy_t, tcol in callouts:
        rr = sel_ix.loc[f]
        xv, yv = float(nl(rr[qr])), float(nl(rr[qe]))
        ax.annotate(
            txt,
            xy=(xv, yv),
            xytext=xy_t,
            fontsize=8.1,
            fontweight="bold",
            color=tcol,
            va="center",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.82),
            arrowprops=dict(arrowstyle="-", color=PALETTE["gray_1"], lw=0.9, shrinkA=1, shrinkB=4),
            zorder=6,
        )

    ax.set_xlim(0, cap)
    ax.set_ylim(0, cap)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlabel("−log10 q, Real alone", fontsize=10.8)
    ax.set_ylabel("−log10 q, Real + PBISC new", fontsize=10.8)


def draw_main_d_scorecard(ax):
    """Credibility scorecard for the 8 PBISC-dependent candidates (10 rows)."""
    cols = [
        "Real\n(ρ / AUC)",
        "PBISC new\n(ρ / AUC)",
        "Homog.\nP_het",
        "Shared ρ\n(n = 27)",
        "Robust\n(of 14)",
        "WHO 2–3\ndeath AUC",
    ]
    GOOD, MID, BAD = PALETTE["green_2"], PALETTE["yellow_soft"], PALETTE["red_1"]

    for i, (_, r) in enumerate(SC.iterrows()):
        o = r["outcome"]
        fmt = (lambda z: f"{z:+.2f}") if o == "sev" else (lambda z: f"{z:.2f}")
        sh_txt = f"{r['shared_rho']:+.2f}{'*' if r['shared_p'] < 0.05 else ''}"
        rob_txt = f"{int(r['robust_pass'])}/{int(r['robust_total'])}"
        cells = [
            (fmt(r["effect_Real"]), "white", COL["Real"]),
            (fmt(r["effect_New"]), "white", COL["New"]),
            (f"{r['p_het']:.2f}", GOOD if r["p_het"] >= 0.05 else BAD, COL["Real"]),
            (sh_txt, GOOD if r["shared_p"] < 0.05 else (MID if r["shared_rho"] >= 0.25 else BAD), COL["Real"]),
            (rob_txt, GOOD if r["robust_pass"] >= 12 else (MID if r["robust_pass"] >= 9 else BAD), COL["Real"]),
            (
                (f"{r['AUC_within_WHO23']:.2f}", GOOD if r["p_within_WHO23"] < 0.05 else BAD, COL["Real"])
                if o == "death"
                else ("–", "white", PALETTE["gray_1"])
            ),
        ]
        for j, (txt, fc, tc) in enumerate(cells):
            if fc != "white":
                ax.add_patch(FancyBboxPatch((j + 0.06, i + 0.10), 0.88, 0.80, boxstyle="round,pad=0.02,rounding_size=0.10", fc=fc, ec="none"))
            ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=9.2, color=tc,
                    fontweight="bold" if (j in (0, 1) and r["followed"]) or (fc == GOOD) else "normal")

    ax.set_xlim(0, len(cols))
    ax.set_ylim(len(SC), 0)
    ax.set_yticks(np.arange(len(SC)) + 0.5)
    ylabs = [f"{LBL[r.feature]} · {'sev.' if r.outcome == 'sev' else 'death'}" for _, r in SC.iterrows()]
    ax.set_yticklabels(ylabs, fontsize=9.5)
    for idx_t, t in enumerate(ax.get_yticklabels()):
        if SC.iloc[idx_t]["followed"]:
            t.set_fontweight("bold")
            t.set_color(PALETTE["red_strong"] if SC.iloc[idx_t]["feature"] == "HSPC" else COL["Ext"])
        else:
            t.set_color(PALETTE["gray_2"])

    ax.set_xticks(np.arange(len(cols)) + 0.5)
    ax.set_xticklabels(cols, fontsize=9.0, fontweight="bold")
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)

    # Separator between Followed (rows 0..5) and Not followed (rows 6..9)
    ax.plot([-0.34, 1.0], [6.0, 6.0], transform=ax.get_yaxis_transform(), color=COL["Real"], lw=1.5, clip_on=False)
    ax.text(-0.33, 3.0, "Followed", transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center",
            fontsize=9.3, fontweight="bold", color=COL["Ext"], clip_on=False)
    ax.text(-0.33, 8.0, "Not followed", transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center",
            fontsize=9.0, color=PALETTE["gray_2"], clip_on=False)


def draw_main_e_forest(ax_sev, ax_dth):
    """Side-by-side forest plots (Real, PBISC new, Combined) for 2 Anchors + 5 Focus populations."""
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
        offsets = [("Real", "Real", +0.26, 40), ("PBISC new", "New", 0.0, 42), ("Real + PBISC new", "Ext", -0.26, 46)]
        for sname, k, dy, sz in offsets:
            row = ci_no[(ci_no["feature"] == f) & (ci_no["set"] == sname)].iloc[0]
            # Severity
            ax_sev.plot([row["rho_lo"], row["rho_hi"]], [y_cur + dy, y_cur + dy], color=COL[k], lw=1.6, alpha=0.75, zorder=2)
            ax_sev.scatter(row["rho"], y_cur + dy, s=sz, marker=MK[k],
                           c="white" if row["p_sev"] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.5, zorder=3)
            # Death
            ax_dth.plot([row["AUC_lo"], row["AUC_hi"]], [y_cur + dy, y_cur + dy], color=COL[k], lw=1.6, alpha=0.75, zorder=2)
            ax_dth.scatter(row["AUC"], y_cur + dy, s=sz, marker=MK[k],
                           c="white" if row["p_death"] >= 0.05 else COL[k], edgecolors=COL[k], linewidths=1.5, zorder=3)
        y_cur -= 1.05

    ax_sev.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
    ax_dth.axvline(0.5, color=PALETTE["gray_1"], ls=":", lw=1.3)

    ax_sev.set_yticks(y_centers)
    ax_sev.set_yticklabels([LBL[f] for f in FEATS_E], fontsize=9.8)
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
    ax_sev.set_xlabel("Spearman ρ (severity)", fontsize=10.3)

    ax_dth.set_xlim(0.0, 1.05)
    ax_dth.set_xticks([0.25, 0.50, 0.75, 1.0])
    ax_dth.set_xticklabels(["0.25", "0.5", "0.75", "1"])
    ax_dth.set_xlabel("AUC (death)", fontsize=10.3)


def draw_strip_ns(ax, f, outcome, rng_local, show_ylabel=True):
    """Patient-level Normal-score strip plot (severity or death) in ver7 house style."""
    if outcome == "sev":
        for g in range(4):
            for k, dx in [("Real", -0.18), ("New", 0.18)]:
                vals = PT_NS[(PT_NS["src"] == k) & (PT_NS["sev"] == g)][f].dropna().values
                x = g + dx + rng_local.uniform(-0.075, 0.075, len(vals))
                ax.scatter(x, vals, s=26, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.45, zorder=3)
                if len(vals):
                    ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=2.0, zorder=4)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["Mild", "Mod.", "Sev.", "Crit."], fontsize=10.0)
        ax.set_xlim(-0.55, 3.55)
    else:
        for g in range(2):
            for k, dx in [("Real", -0.18), ("New", 0.18)]:
                vals = PT_NS[(PT_NS["src"] == k) & (PT_NS["fatal"] == g)][f].dropna().values
                x = g + dx + rng_local.uniform(-0.075, 0.075, len(vals))
                ax.scatter(x, vals, s=28, marker=MK[k], c=COL[k], edgecolors="white", linewidths=0.45, zorder=3)
                if len(vals):
                    ax.plot([g + dx - 0.14, g + dx + 0.14], [np.median(vals)] * 2, color="black", lw=2.0, zorder=4)
        n_surv = int((PT_NS["fatal"] == 0).sum())
        n_died = int((PT_NS["fatal"] == 1).sum())
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Survived\n({n_surv})", f"Died\n({n_died})"], fontsize=9.8)
        ax.set_xlim(-0.52, 1.52)

    ax.set_ylim(-2.6, 5.05)
    ax.set_yticks([-2, 0, 2])
    if show_ylabel:
        ax.set_ylabel("Normal score (within-source)", fontsize=10.2)

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
        ax.text(0.03, 0.975 - 0.090 * idx_s, txt, transform=ax.transAxes, fontsize=8.6,
                color=COL[k], va="top", fontweight="bold" if sig else "normal")


def draw_main_l_hspc_analyses(ax):
    """HSPC mortality AUC across 8 definitions / strata in ver7 diamond/square/circle style."""
    ORDER_H = [
        ("All deaths", "All deaths"),
        ("Table S1 deaths", "Table S1 deaths"),
        ("Direct COVID-19 deaths", "Direct COVID-19"),
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
            ax.scatter(rn["AUC"], yi, s=52, marker="o", c=COL["New"], edgecolors=COL["New"], linewidths=1.5, zorder=4)
            ax.text(1.02, yi, f"{int(rn['deaths'])} · {fp(rn['p']).replace('P = ', '')}", va="center", ha="left",
                    fontsize=8.2, color=COL["New"], fontweight="bold")
        else:
            rr = HS[(HS["analysis"] == an_key) & (HS["set"] == "Real")].iloc[0]
            re = HS[(HS["analysis"] == an_key) & (HS["set"] == "Real + PBISC new")].iloc[0]
            ax.plot([rr["AUC"], re["AUC"]], [yi + 0.16, yi - 0.16], color=PALETTE["neutral"], lw=1.2, zorder=2)
            ax.scatter(rr["AUC"], yi + 0.16, s=44, marker="D",
                       c="white" if rr["p"] >= 0.05 else COL["Real"], edgecolors=COL["Real"], linewidths=1.5, zorder=3)
            ax.scatter(re["AUC"], yi - 0.16, s=48, marker="s",
                       c="white" if re["p"] >= 0.05 else COL["Ext"], edgecolors=COL["Ext"], linewidths=1.5, zorder=4)
            p_str = fp(re["p"]).replace("P = ", "")
            ax.text(1.02, yi, f"{int(rr['deaths'])}→{int(re['deaths'])} · {p_str}", va="center", ha="left",
                    fontsize=8.2, color=COL["Ext"], fontweight="bold" if re["p"] < 0.05 else "normal")

    ax.text(1.02, yy[0] + 0.72, "deaths · +PBISC P", va="bottom", ha="left", fontsize=8.0, color=PALETTE["gray_2"], fontweight="bold")
    ax.set_yticks(yy)
    ax.set_yticklabels([lab for _, lab in ORDER_H], fontsize=9.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0.46, 1.00)
    ax.set_ylim(-0.6, len(ORDER_H) - 0.25)
    ax.set_xticks([0.5, 0.75, 1.0])
    ax.set_xticklabels(["0.5", "0.75", "1.0"])
    ax.set_xlabel("HSPC death AUC", fontsize=10.3)


def draw_main_m_power(ax):
    """Power at 20 Real patients (Real, +PBISC new, +PBISC shuffled) + n80 sample-size reduction."""
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
        ax.barh(yy + dy, vals, height=hb, color=fc, edgecolor=COL["Real"], linewidth=1.1, hatch=hatch, zorder=3)

    # Annotate n80 Real sample size reduction on the right
    for yi, (f, o, _) in zip(yy, PP):
        d_sub = PW[(PW["feature"] == f) & (PW["outcome"] == o)].sort_values("n_real")
        n_r, _ = n80_interp(d_sub["n_real"], d_sub["real_alone"])
        n_e, le = n80_interp(d_sub["n_real"], d_sub["plus_pbisc_new"])
        r_txt = f"{n_r:.0f}" if np.isfinite(n_r) else "n.r."
        e_txt = "≤10" if le else f"{n_e:.0f}"
        ax.text(1.03, yi, f"{r_txt}→{e_txt}", va="center", ha="left", fontsize=8.3, color=COL["Ext"], fontweight="bold")

    ax.text(1.03, yy[0] + 0.65, "n80 Real", va="bottom", ha="left", fontsize=8.0, color=PALETTE["gray_2"], fontweight="bold")
    ax.axvline(0.80, color=PALETTE["gray_1"], ls=":", lw=1.4, zorder=2)
    ax.set_yticks(yy)
    ax.set_yticklabels([lab for _, _, lab in PP], fontsize=9.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0.0, 1.01)
    ax.set_ylim(-0.55, len(PP) - 0.2)
    ax.set_xticks([0.0, 0.5, 0.8, 1.0])
    ax.set_xticklabels(["0", "0.5", "0.8", "1"])
    ax.set_xlabel("Power at 20 Real patients", fontsize=10.3)


# =========================================================================
# 1) RENDER MAIN FIGURE 7 COMPOSITE (a - m)
# =========================================================================
fig_m = plt.figure(figsize=(15.2, 21.4), facecolor="white")

fig_m.text(
    0.015,
    0.990,
    "Figure 7. Pre-specified screening and validation of PBISC-rescued immune associations in ESKD COVID-19",
    fontsize=15.0,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

# Row 1: a (CLL molecular screen), b (Selection funnel), c (Gain map)
ax_a = fig_m.add_axes([0.055, 0.825, 0.215, 0.122])
draw_main_a_cll_screen(ax_a)
ax_a.set_title("CLL molecular screen", loc="left", fontsize=11.5, fontweight="bold", pad=8)
fig_m.text(0.015, 0.962, "a", fontsize=18, fontweight="bold", va="top")

ax_b = fig_m.add_axes([0.388, 0.825, 0.245, 0.122])
draw_main_b_funnel(ax_b)
ax_b.set_title("Selection steps (43 cell types)", loc="left", fontsize=11.5, fontweight="bold", pad=8)
fig_m.text(0.292, 0.962, "b", fontsize=18, fontweight="bold", va="top")

ax_c = fig_m.add_axes([0.715, 0.825, 0.260, 0.122])
draw_main_c_gain_map(ax_c)
ax_c.set_title("Associations gained by adding PBISC", loc="left", fontsize=11.5, fontweight="bold", pad=8)
fig_m.text(0.665, 0.962, "c", fontsize=18, fontweight="bold", va="top")

# Dedicated Legend Strip below Row 1
ax_leg = fig_m.add_axes([0.015, 0.768, 0.97, 0.024])
ax_leg.set_axis_off()
hd_main = [
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc=COL["Real"], mec=COL["Real"], label="Real (37 patients, 7 deaths)"),
    Line2D([0], [0], marker="o", ls="", ms=7.5, mfc=COL["New"], mec=COL["New"], label="PBISC new (29 patients, 6 deaths)"),
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["Ext"], mec=COL["Ext"], label="Real + PBISC new (66 patients, 13 deaths)"),
    Line2D([0], [0], marker="o", ls="", ms=7.5, mfc="white", mec=PALETTE["gray_2"], mew=1.6, label="open: P ≥ 0.05 (e, l)"),
    Line2D([0], [0], marker="o", ls="", ms=6.5, mfc=COL["Ext"], mec="black", label="severity (c)"),
    Line2D([0], [0], marker="^", ls="", ms=7.0, mfc=COL["Ext"], mec="black", label="death (c)"),
    Patch(facecolor="white", edgecolor=COL["Real"], hatch="///", label="PBISC new, outcome shuffled (m)"),
    Patch(facecolor=PALETTE["green_2"], label="criterion met (d)"),
    Patch(facecolor=PALETTE["yellow_soft"], label="borderline (d)"),
    Patch(facecolor=PALETTE["red_1"], label="not met (d)"),
]
ax_leg.legend(handles=hd_main, loc="center", ncol=5, fontsize=9.6, frameon=False, columnspacing=1.15, handletextpad=0.35, labelspacing=0.35)

# Row 2: d (Credibility scorecard) & e (Forest CIs for Anchors + Focus populations)
ax_d = fig_m.add_axes([0.145, 0.554, 0.425, 0.156])
draw_main_d_scorecard(ax_d)
ax_d.set_title("Credibility of PBISC-dependent candidates", loc="left", fontsize=11.5, fontweight="bold", pad=34)
fig_m.text(0.015, 0.748, "d", fontsize=18, fontweight="bold", va="top")

ax_e_sev = fig_m.add_axes([0.688, 0.558, 0.132, 0.162])
ax_e_dth = fig_m.add_axes([0.845, 0.558, 0.132, 0.162])
draw_main_e_forest(ax_e_sev, ax_e_dth)
ax_e_sev.set_title("Cohort effect sizes (95% CI)", loc="left", fontsize=11.5, fontweight="bold", pad=10)
fig_m.text(0.592, 0.748, "e", fontsize=18, fontweight="bold", va="top")

# Row 3: f, g, h (Severity: B ASC IgG, B ASC IgA, T CD8 CTL), i (Death: T CD8 CM)
# CRITICAL REQUIREMENT: Top line of f-k has ONLY the cell population name!
rng_main = np.random.default_rng(7)
row3_specs = [
    ("B_ASC_IgG", "sev", "f", True),
    ("B_ASC_IgA", "sev", "g", False),
    ("T CD8 CTL", "sev", "h", False),
    ("T CD8 CM", "death", "i", True),
]
for idx_p, (f, out_t, let, show_yl) in enumerate(row3_specs):
    x_l = 0.062 + idx_p * 0.242
    ax_s = fig_m.add_axes([x_l, 0.315, 0.182, 0.162])
    draw_strip_ns(ax_s, f, out_t, rng_main, show_ylabel=show_yl)
    ax_s.set_title(LBL[f], loc="left", fontsize=11.5, fontweight="bold", color=COL["Real"], pad=6)
    fig_m.text(x_l - 0.044, 0.495, let, fontsize=18, fontweight="bold", va="top")

# Row 4: j (T CD8 CTL death), k (HSPC death), l (HSPC death analyses), m (Power at 20 Real & n80)
ax_j = fig_m.add_axes([0.062, 0.058, 0.165, 0.165])
draw_strip_ns(ax_j, "T CD8 CTL", "death", rng_main, show_ylabel=True)
ax_j.set_title("T CD8 CTL", loc="left", fontsize=11.5, fontweight="bold", color=COL["Real"], pad=6)
fig_m.text(0.018, 0.242, "j", fontsize=18, fontweight="bold", va="top")

ax_k = fig_m.add_axes([0.268, 0.058, 0.165, 0.165])
draw_strip_ns(ax_k, "HSPC", "death", rng_main, show_ylabel=False)
ax_k.set_title("HSPC", loc="left", fontsize=11.5, fontweight="bold", color=COL["Real"], pad=6)
fig_m.text(0.232, 0.242, "k", fontsize=18, fontweight="bold", va="top")

ax_l = fig_m.add_axes([0.538, 0.058, 0.118, 0.165])
draw_main_l_hspc_analyses(ax_l)
ax_l.set_title("HSPC death robustness", loc="left", fontsize=11.5, fontweight="bold", color=COL["Real"], pad=6)
fig_m.text(0.448, 0.242, "l", fontsize=18, fontweight="bold", va="top")

ax_m = fig_m.add_axes([0.844, 0.058, 0.100, 0.165])
draw_main_m_power(ax_m)
ax_m.set_title("Power & n80 gain", loc="left", fontsize=11.5, fontweight="bold", pad=6)
fig_m.text(0.754, 0.242, "m", fontsize=18, fontweight="bold", va="top")

main_out = V8_DIR / "PBISC_Figure7_Main_ver8_20261117"
saved_main = finalize_figure(fig_m, main_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Main Figure 7:", [str(p) for p in saved_main])


# =========================================================================
# SUPPLEMENTARY FIGURE S1: CLL Audit, Concordant Anchors, NK Decomposition &
#                          Normal-Score vs Raw % of PBMCs Modality Proof
# =========================================================================
def draw_s1_a_cll_bars(ax):
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


def draw_s1_bc_cll_shift(ax, val, lo, hi, pcol, ref, xlab, xlim, xt):
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


def draw_s1_h_nk_decomp(ax):
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


def draw_s1_i_modality_scale(ax):
    """Paired 27-shared-donor Real CITE-seq vs PBISC bulk % PBMC showing systematic platform scale offsets."""
    FEATS_M = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
    x = np.arange(len(FEATS_M))
    rng_m = np.random.default_rng(19)
    for i, f in enumerate(FEATS_M):
        vr = SH27[f + "|Real_pct"].values
        vb = SH27[f + "|PBISC_pct"].values
        # Paired connecting lines in faint neutral
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
    """Patient-level RAW % of PBMCs strip plot (for direct comparison with Normal score)."""
    # Use log10(% + 0.02) tick transform or linear %? Wait: user asked about % of total cells on y-axis!
    # Let's plot exact % of PBMCs on a clean symlog / square-root or linear % axis!
    # Even better: linear % of PBMCs with clean cohort medians & pooled raw-% stats so values are 100% literal %!
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


fig_s1 = plt.figure(figsize=(15.2, 20.6), facecolor="white")
fig_s1.text(
    0.015,
    0.990,
    "Supplementary Figure 1. CLL exclusion audit, concordant anchors, NK decomposition, and within-source normal-score justification",
    fontsize=13.8,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

# Row 1 of S1: a (CLL composition), b (CLL shift severity), c (CLL shift death)
ax_s1a = fig_s1.add_axes([0.125, 0.838, 0.330, 0.110])
draw_s1_a_cll_bars(ax_s1a)
fig_s1.text(0.015, 0.962, "a", fontsize=17, fontweight="bold", va="top")

ax_s1b = fig_s1.add_axes([0.575, 0.832, 0.165, 0.122])
draw_s1_bc_cll_shift(ax_s1b, "rho", "rho_lo", "rho_hi", "p_sev", 0.0, "ρ with severity (CLL incl. → excl.)", (-0.88, 0.88), [-0.5, 0, 0.5])
fig_s1.text(0.482, 0.962, "b", fontsize=17, fontweight="bold", va="top")

ax_s1c = fig_s1.add_axes([0.815, 0.832, 0.165, 0.122])
draw_s1_bc_cll_shift(ax_s1c, "AUC", "AUC_lo", "AUC_hi", "p_death", 0.5, "AUC for death (CLL incl. → excl.)", (0.05, 1.02), [0.25, 0.5, 0.75, 1.0])
fig_s1.text(0.755, 0.962, "c", fontsize=17, fontweight="bold", va="top")

# Legend strip on S1
ax_s1_leg = fig_s1.add_axes([0.02, 0.780, 0.96, 0.020])
ax_s1_leg.set_axis_off()
hd_s1 = [
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc=COL["Real"], mec=COL["Real"], label="Real (37, CLL excluded)"),
    Line2D([0], [0], marker="o", ls="", ms=7.5, mfc=COL["New"], mec=COL["New"], label="PBISC new (29) / PBISC shared (i)"),
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["Ext"], mec=COL["Ext"], label="Real + PBISC new (66)"),
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=COL["All"], mec=COL["All"], label="all patients, CLL included (b, c)"),
    Line2D([0], [0], marker="o", ls="", ms=7.5, mfc="white", mec=PALETTE["gray_2"], mew=1.5, label="open: P ≥ 0.05"),
]
ax_s1_leg.legend(handles=hd_s1, loc="center", ncol=5, fontsize=9.6, frameon=False, columnspacing=1.3)

# Row 2 of S1: d-g (4 Cross-Cohort Concordant Anchors; top title = cell population name only)
rng_s1 = np.random.default_rng(17)
anchors_s1 = [("NK CD56 bright", "d"), ("B_ASC_dividing", "e"), ("T CD8 CM", "f"), ("CD16mono", "g")]
for idx_p, (f, let) in enumerate(anchors_s1):
    x_l = 0.062 + idx_p * 0.242
    ax_an = fig_s1.add_axes([x_l, 0.565, 0.182, 0.165])
    draw_strip_ns(ax_an, f, "sev", rng_s1, show_ylabel=(idx_p == 0))
    ax_an.set_title(LBL[f], loc="left", fontsize=11.2, fontweight="bold", color=COL["Real"], pad=6)
    fig_s1.text(x_l - 0.044, 0.748, let, fontsize=17, fontweight="bold", va="top")

# Row 3 of S1: h (NK compartment decomposition) & i (27-shared-donor platform scale proof)
ax_s1h = fig_s1.add_axes([0.070, 0.318, 0.240, 0.168])
draw_s1_h_nk_decomp(ax_s1h)
ax_s1h.set_title("NK compartment decomposition", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_s1.text(0.015, 0.505, "h", fontsize=17, fontweight="bold", va="top")

ax_s1i = fig_s1.add_axes([0.405, 0.318, 0.570, 0.168])
draw_s1_i_modality_scale(ax_s1i)
ax_s1i.set_title("Modality scale shift on 27 shared donors (PBISC / Real median fold & paired Wilcoxon)", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_s1.text(0.345, 0.505, "i", fontsize=17, fontweight="bold", va="top")

# Row 4 of S1: j-o (6 Raw % of PBMCs distributions corresponding to Main f-k; top title = cell population name only)
raw_specs = [
    ("B_ASC_IgG", "sev", "j"),
    ("B_ASC_IgA", "sev", "k"),
    ("T CD8 CTL", "sev", "l"),
    ("T CD8 CM", "death", "m"),
    ("T CD8 CTL", "death", "n"),
    ("HSPC", "death", "o"),
]
for idx_r, (f, out_t, let) in enumerate(raw_specs):
    x_l = 0.052 + idx_r * 0.158
    ax_rw = fig_s1.add_axes([x_l, 0.050, 0.125, 0.175])
    draw_strip_raw_pct(ax_rw, f, out_t, rng_s1, show_ylabel=(idx_r == 0))
    ax_rw.set_title(LBL[f], loc="left", fontsize=10.8, fontweight="bold", color=COL["Real"], pad=5)
    fig_s1.text(x_l - 0.030, 0.242, let, fontsize=17, fontweight="bold", va="top")

s1_out = V8_DIR / "PBISC_SuppFig_S1_CLL_Anchors_Modality_20261117"
saved_s1 = finalize_figure(fig_s1, s1_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Supplementary Figure S1:", [str(p) for p in saved_s1])


# =========================================================================
# SUPPLEMENTARY FIGURE S2: Full 6-Row Screen, Continuous Power Curves,
#                          CD14 mono IFN Audit & 14-Variant Robustness Matrix
# =========================================================================
def draw_s2_a_screen(fig, ax, cax_rect):
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
            size = 34 + 195 * min(-np.log10(max(r[qcol], 1e-4)), 4.0) / 4.0
            ax.scatter(i, y_pos, s=size, c=[cmap(norm(val))], edgecolors="black" if is_fdr else "none",
                       linewidths=1.6 if is_fdr else 0.0, zorder=3)

    ax.set_yticks([r[5] for r in ROWS])
    ax.set_yticklabels([r[4] for r in ROWS], fontsize=10.2)
    ax.set_xticks(range(len(keep)))
    ax.set_xticklabels([LBL[f] for f in keep], rotation=42, ha="right", fontsize=10.0)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.6, len(keep) - 0.4)
    ax.set_ylim(-5.95, 0.95)

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
            ax.text(i, 0.65, flag, ha="center", va="center", fontsize=10.0, color=PALETTE["red_strong"], fontweight="bold")

    cax = fig.add_axes(cax_rect)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm, cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-0.8, 0, 0.8])
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.tick_params(labelsize=9.5, length=2, pad=1)
    cb.outline.set_visible(False)


def draw_s2_b_n80(ax):
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
        ax.plot([b, ra], [y_cur, y_cur], color=PALETTE["neutral"], lw=2.5, zorder=1)
        ax.scatter(ra, y_cur, marker="D", s=48, zorder=3, c=COL["Real"] if np.isfinite(a) else "white", edgecolors=COL["Real"], linewidths=1.6)
        ax.scatter(b, y_cur, marker="s", s=52, zorder=4, c=ext_col, edgecolors=ext_col)
        ax.text(ra + 1.2, y_cur, f"{a:.0f}" if np.isfinite(a) else "n.r.", va="center", fontsize=9.2, color=COL["Real"])
        ax.text(b - 1.2, y_cur, "≤10" if le else f"{b:.0f}", va="center", ha="right", fontsize=9.2, color=ext_col, fontweight="bold")
        y_cur -= 1.0

    ax.set_yticks(y_positions)
    ax.set_yticklabels([it[2] for it in ITEMS], fontsize=9.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(1, 45.5)
    ax.set_xticks([10, 20, 30, 40])
    ax.set_xticklabels(["10", "20", "30", "n.r."])
    ax.set_ylim(y_positions[-1] - 0.55, 0.6)
    ax.set_xlabel("Real patients needed for 80% power", fontsize=10.2)


def draw_s2_c_power_curves(ax):
    d_igg = PW[(PW["feature"] == "B_ASC_IgG") & (PW["outcome"] == "sev")].sort_values("n_real")
    d_iga = PW[(PW["feature"] == "B_ASC_IgA") & (PW["outcome"] == "sev")].sort_values("n_real")
    d_hsp = PW[(PW["feature"] == "HSPC") & (PW["outcome"] == "death")].sort_values("n_real")
    d_cm = PW[(PW["feature"] == "T CD8 CM") & (PW["outcome"] == "death")].sort_values("n_real")

    ax.plot(d_igg["n_real"], d_igg["real_alone"], color=COL["Ext"], lw=1.8, ls="--", marker="D", ms=4.5, label="B ASC IgG (sev): Real")
    ax.plot(d_igg["n_real"], d_igg["plus_pbisc_new"], color=COL["Ext"], lw=2.3, ls="-", marker="s", ms=5.0, label="B ASC IgG (sev): +PBISC")

    ax.plot(d_iga["n_real"], d_iga["real_alone"], color=COL["New"], lw=1.7, ls="--", marker="D", ms=4.2, mfc="white", label="B ASC IgA (sev): Real")
    ax.plot(d_iga["n_real"], d_iga["plus_pbisc_new"], color=COL["New"], lw=2.1, ls="-", marker="s", ms=4.8, label="B ASC IgA (sev): +PBISC")

    ax.plot(d_hsp["n_real"], d_hsp["real_alone"], color=PALETTE["red_strong"], lw=1.8, ls="--", marker="D", ms=4.5, mfc="white", label="HSPC (death): Real")
    ax.plot(d_hsp["n_real"], d_hsp["plus_pbisc_new"], color=PALETTE["red_strong"], lw=2.3, ls="-", marker="s", ms=5.0, label="HSPC (death): +PBISC")

    ax.plot(d_cm["n_real"], d_cm["plus_pbisc_new"], color=PALETTE["teal"], lw=2.0, ls="-", marker="s", ms=4.5, label="T CD8 CM (death): +PBISC")

    shuf_avg = (d_igg["plus_shuffled"].values + d_iga["plus_shuffled"].values + d_hsp["plus_shuffled"].values) / 3.0
    ax.plot(d_igg["n_real"], shuf_avg, color=PALETTE["gray_1"], lw=1.6, ls=":", label="+PBISC shuffled (mean)")

    ax.axhline(0.8, color=PALETTE["gray_1"], ls=":", lw=1.2)
    ax.set_ylim(-0.03, 1.44)
    ax.set_yticks([0, 0.4, 0.8, 1.0])
    ax.set_xticks([10, 15, 20, 25, 30, 35])
    ax.set_xlim(8.5, 36.5)
    ax.set_xlabel("Real patients sampled (1,000 draws)", fontsize=10.2)
    ax.set_ylabel("Power (P < 0.05)", fontsize=10.2)
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=2, fontsize=7.8, handlelength=1.4, columnspacing=0.55, handletextpad=0.25, labelspacing=0.18)


def draw_s2_e_sensitivity_matrix(ax_sev, ax_death):
    SA_ext = SA[SA["set"] == "Real + PBISC new"]
    FEN = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN"]
    VN = {
        "0 baseline (no CLL)": "Baseline (CLL excluded)",
        "A1 adjust centre/wave": "Adjust centre / wave",
        "A2 within-lineage %": "% of own lineage",
        "A3 % of non-myeloid cells": "% of non-myeloid cells",
        "A4 drop Real samples with mono <2%": "Drop monocyte-depleted samples",
        "B  also drop C187, C190": "Also drop C187 (NK) & C190 (T) outliers",
        "C1 death = Table S1 (C101, C147 alive)": "Death per Table S1 (acute COVID)",
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
                    ax.text(j, i, "–", ha="center", va="center", fontsize=9.0, color=PALETTE["gray_1"])
                    continue
                sig = Pm.values[i, j] < 0.05
                ax.text(
                    j,
                    i,
                    f"{v:.2f}",
                    ha="center",
                    va="center",
                    fontsize=8.8,
                    fontweight="bold" if sig else "normal",
                    color=("white" if abs(nrm(v) - 0.5) > 0.32 else "black") if sig else PALETTE["gray_1"],
                )
        ax.axvline(3.5, color="black", lw=1.6, ls="-")
        ax.set_xticks(range(len(FEN)))
        ax.set_xticklabels([LBL[f] for f in FEN], rotation=32, ha="left", fontsize=9.6)
        ax.xaxis.tick_top()
        ax.set_yticks(range(len(VN)))
        ax.set_yticklabels(list(VN.values()) if show_y else [], fontsize=9.6)
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
        ax.set_title(ttl, loc="left", fontsize=11.2, fontweight="bold", pad=28)


fig_s2 = plt.figure(figsize=(15.2, 18.2), facecolor="white")
fig_s2.text(
    0.015,
    0.990,
    "Supplementary Figure 2. Full 43-cell-type screen, subsampling power curves, and 14-variant confounder sensitivity audit",
    fontsize=13.8,
    fontweight="bold",
    va="top",
    ha="left",
    color=COL["Real"],
)

# Row 1 of S2: a (Full 6-row Fine-label screen)
ax_s2a = fig_s2.add_axes([0.195, 0.775, 0.775, 0.145])
draw_s2_a_screen(fig_s2, ax_s2a, cax_rect=[0.195, 0.942, 0.125, 0.010])
fig_s2.text(
    0.332,
    0.947,
    "ρ (severity) or 2·(AUC − 0.5) (death)  ·  dot size: −log10 q (BH, 43 fine types)  ·  ring: q < 0.05 (sev) or q < 0.10 (death)  ·  H: het. P < 0.05  ·  S: steroid",
    fontsize=9.2,
    va="center",
    color=COL["Real"],
)
fig_s2.text(0.015, 0.960, "a", fontsize=17, fontweight="bold", va="top")

# Row 2 of S2: b (n80 dumbbell), c (continuous power curves), d (CD14 mono IFN audit)
ax_s2b = fig_s2.add_axes([0.160, 0.435, 0.215, 0.220])
draw_s2_b_n80(ax_s2b)
ax_s2b.set_title("Sample size for 80% power", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_s2.text(0.015, 0.672, "b", fontsize=17, fontweight="bold", va="top")

ax_s2c = fig_s2.add_axes([0.455, 0.435, 0.265, 0.220])
draw_s2_c_power_curves(ax_s2c)
ax_s2c.set_title("Subsampling power curves (n = 10–35)", loc="left", fontsize=11.2, fontweight="bold", pad=6)
fig_s2.text(0.405, 0.672, "c", fontsize=17, fontweight="bold", va="top")

ax_s2d = fig_s2.add_axes([0.795, 0.435, 0.175, 0.220])
draw_strip_ns(ax_s2d, "CD14mono_IFN", "death", np.random.default_rng(23), show_ylabel=True)
ax_s2d.text(0.03, 0.69, "1st sample AUC=0.78 (P=0.002)\nCentre-adj P=0.071 (Step 4 filter)",
            transform=ax_s2d.transAxes, fontsize=8.3, color=PALETTE["red_strong"], fontweight="bold", va="top")
ax_s2d.set_title("CD14 mono IFN", loc="left", fontsize=11.2, fontweight="bold", color=COL["Real"], pad=6)
fig_s2.text(0.742, 0.672, "d", fontsize=17, fontweight="bold", va="top")

# Row 3 of S2: e (14-variant sensitivity matrix)
ax_s2e_s = fig_s2.add_axes([0.225, 0.025, 0.365, 0.275])
ax_s2e_d = fig_s2.add_axes([0.612, 0.025, 0.365, 0.275])
draw_s2_e_sensitivity_matrix(ax_s2e_s, ax_s2e_d)
fig_s2.text(0.015, 0.368, "e", fontsize=17, fontweight="bold", va="top")

s2_out = V8_DIR / "PBISC_SuppFig_S2_Screen_Power_Sensitivity_20261117"
saved_s2 = finalize_figure(fig_s2, s2_out, formats=("png", "pdf", "svg"), dpi=300, close=True)
print("Saved Supplementary Figure S2:", [str(p) for p in saved_s2])


# =========================================================================
# 4) EXPORT STANDALONE MODULAR PANELS (PNG / PDF / SVG)
# =========================================================================
def export_single(sub_dir, name, figsize, draw_fn):
    f_s = plt.figure(figsize=figsize, facecolor="white")
    draw_fn(f_s)
    finalize_figure(f_s, PANEL_DIR / sub_dir / name, formats=("png", "pdf", "svg"), dpi=300, close=True)


export_single("main", "Fig7_a_CLL_molecular_screen", (4.8, 4.4), lambda fg: draw_main_a_cll_screen(fg.add_axes([0.18, 0.16, 0.78, 0.78])))
export_single("main", "Fig7_b_selection_steps", (5.8, 4.0), lambda fg: draw_main_b_funnel(fg.add_axes([0.32, 0.18, 0.64, 0.76])))
export_single("main", "Fig7_c_gain_map", (5.4, 4.6), lambda fg: draw_main_c_gain_map(fg.add_axes([0.16, 0.16, 0.80, 0.78])))
export_single("main", "Fig7_d_credibility_scorecard", (8.6, 4.6), lambda fg: draw_main_d_scorecard(fg.add_axes([0.26, 0.06, 0.72, 0.78])))
export_single("main", "Fig7_e_cohort_effect_sizes", (6.8, 4.8), lambda fg: draw_main_e_forest(fg.add_axes([0.26, 0.16, 0.33, 0.78]), fg.add_axes([0.64, 0.16, 0.33, 0.78])))
for fk, ot, lt in [("B_ASC_IgG", "sev", "f"), ("B_ASC_IgA", "sev", "g"), ("T CD8 CTL", "sev", "h"),
                   ("T CD8 CM", "death", "i"), ("T CD8 CTL", "death", "j"), ("HSPC", "death", "k")]:
    export_single("main", f"Fig7_{lt}_{fk.replace(' ', '_')}_{ot}", (3.8, 3.8),
                  lambda fg, f_in=fk, o_in=ot: (
                      ax := fg.add_axes([0.24, 0.18, 0.72, 0.74]),
                      draw_strip_ns(ax, f_in, o_in, np.random.default_rng(7), True),
                      ax.set_title(LBL[f_in], loc="left", fontsize=11.5, fontweight="bold", pad=6),
                  ))
export_single("main", "Fig7_l_HSPC_death_analyses", (4.8, 4.2), lambda fg: draw_main_l_hspc_analyses(fg.add_axes([0.36, 0.16, 0.38, 0.76])))
export_single("main", "Fig7_m_power_and_n80", (4.8, 4.0), lambda fg: draw_main_m_power(fg.add_axes([0.36, 0.18, 0.42, 0.74])))

export_single("supp_s1", "SuppS1_a_CLL_PBMC_composition", (6.8, 3.4), lambda fg: draw_s1_a_cll_bars(fg.add_axes([0.24, 0.20, 0.72, 0.68])))
export_single("supp_s1", "SuppS1_h_NK_decomposition", (4.2, 3.6), lambda fg: draw_s1_h_nk_decomp(fg.add_axes([0.24, 0.22, 0.72, 0.72])))
export_single("supp_s1", "SuppS1_i_shared27_modality_scale", (8.2, 4.2), lambda fg: draw_s1_i_modality_scale(fg.add_axes([0.14, 0.24, 0.83, 0.70])))
export_single("supp_s2", "SuppS2_a_full_screen", (10.5, 4.6), lambda fg: draw_s2_a_screen(fg, fg.add_axes([0.22, 0.32, 0.75, 0.56]), [0.22, 0.93, 0.18, 0.03]))
export_single("supp_s2", "SuppS2_e_14variant_sensitivity", (13.8, 4.8), lambda fg: draw_s2_e_sensitivity_matrix(fg.add_axes([0.22, 0.05, 0.37, 0.72]), fg.add_axes([0.61, 0.05, 0.37, 0.72])))

print("Saved all standalone panels in", PANEL_DIR)
