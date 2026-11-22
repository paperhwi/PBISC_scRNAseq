"""Build Nature Portfolio 180-mm Main Figure 5, Main Figure 6, and Supplementary Figures 1-5
aligned 1:1 with the finalized methodological framework:
  - Figure 5: Cross-cohort bulk-only comparison (122 Sepsis, 164 COVID-19 vs separate Real cohorts)
              * Patient x condition cell proportions + patient-unit resampling (49,999 bootstrap)
              * Cell-type-specific literature disease gene scores (Ye et al. 2026; Stephenson et al. 2025)
              * Welch's t-test DEGs (BH q < 0.05 and |log2FC| >= 1) across T, B, NK, Monocyte
  - Figure 6: Rare cell state frequency & sample-summed expression (14 Sepsis, 13 COVID-19 states)
              * Case-group parent-lineage shares (total cells denominator for HSC / HSPC)
              * Rare-state Harmony UMAP localization & 82-donor Sepsis monocyte reference audit
              * Per-sample summed state expression + inter-patient (pt_i != pt_j) Pearson r vs Real
              * Per-sample summed canonical marker expression heatmaps
  - Supplementary Figures 1-5: Coverage/depth/mixing/decomposition, 8 rare abundances, 8 summed marker
              heatmaps, and complete 4-lineage x 3-source Welch DEG volcanoes for Sepsis & COVID-19.
"""
import json
import os
import pathlib
import pickle
import sys
import tempfile
import warnings
import zlib

os.environ.setdefault("MPLCONFIGDIR", str(pathlib.Path(tempfile.gettempdir()) / "mpl_pbisc_fig56"))
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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

from nature_style import (
    CT_COL,
    NATURE_2COL_IN,
    PALETTE,
    SRC_COL,
    SRC_MK,
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

P_DIR = ROOT / "deliverables" / "manuscript_real_pbisc_scsemi_20261110" / "panels"
F4P_DIR = ROOT / "deliverables" / "manuscript_figures4papers_20261111"

SRC = ["Real", "PBISC", "scSemiProfiler"]
T6 = ["T cell", "NK cell", "B cell", "Monocyte", "DC", "Other"]
CT4 = ["T cell", "B cell", "NK cell", "Monocyte"]
NAME = {"sepsis": "Sepsis", "covid": "COVID-19"}
L = "abcdefghijklmnopqrstuvwxyz"

CMAP_CORR = LinearSegmentedColormap.from_list(
    "pbisc_corr", ["#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"], N=256
)
CMAP_EXPR = LinearSegmentedColormap.from_list(
    "pbisc_expr", ["#EAEEF3", "#F8FAFC", "#9ECAE1", "#3775BA", "#0F4D92", "#061E42"], N=256
)

MOD_COLOR_MAP = {
    "Plasmablast": PALETTE["violet"],
    "Macrophage": PALETTE["teal"],
    "CD4 T": PALETTE["green_3"],
    "CD8 T": PALETTE["blue_secondary"],
    "Monocyte": PALETTE["magenta"],
    "Cell cycle": PALETTE["gray_2"],
    "IgG ASC": PALETTE["violet"],
    "Dex Mono": PALETTE["teal"],
    "IFN response": PALETTE["blue_secondary"],
    "pDC": PALETTE["green_3"],
    "HSPC": PALETTE["magenta"],
}


def hex_to_rgb(hx: str):
    hx = hx.lstrip("#")
    return [int(hx[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]


def _resolve_1d_positions(targets, min_sep, lo, hi):
    pos = np.array(targets, dtype=float)
    n = len(pos)
    if n <= 1:
        return np.clip(pos, lo, hi)
    for _ in range(80):
        for i in range(n - 1):
            gap = pos[i + 1] - pos[i]
            if gap < min_sep:
                mid = 0.5 * (pos[i] + pos[i + 1])
                pos[i] = mid - 0.5 * min_sep
                pos[i + 1] = mid + 0.5 * min_sep
        pos = np.clip(pos, lo, hi)
    return pos


# =============================================================================
# MAIN FIGURE 5 HELPERS — CROSS-COHORT BULK-ONLY COMPARISON
# =============================================================================
def draw_fig5_schematic(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, text, fc, ec=PALETTE["gray_2"], bold=False, fs=5.15):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=2.0", fc=fc, ec=ec, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", color=PALETTE["gray_3"], linespacing=1.18)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=PALETTE["gray_2"], lw=0.85, mutation_scale=6.5))

    box(0.8, 55.0, 29.5, 39.0, "Bulk-only cohorts\nSepsis: n = 122\nCOVID-19: n = 164", PALETTE["green_1"], bold=True)
    box(35.2, 55.0, 29.5, 39.0, "Adapted PBISC\n4,000 cells / sample\n(488k & 656k cells)", "#DCE6F2", ec=PALETTE["blue_main"], bold=True)
    box(69.6, 55.0, 29.6, 39.0, "Separate Real scRNA\nSepsis: n = 108\nCOVID: n = 187", "#F4F5F7", fs=5.1)
    arrow(30.8, 74.5, 34.7, 74.5)
    arrow(65.2, 74.5, 69.1, 74.5)

    ax.add_patch(FancyBboxPatch((0.8, 5.0), 98.4, 42.5, boxstyle="round,pad=0.25,rounding_size=2.0", fc="#FAFBFD", ec="#D5DCE4", lw=0.7))
    ax.text(3.0, 41.0, "Cross-cohort rules (Fig. 5):", fontsize=5.5, fontweight="bold", color=PALETTE["blue_main"], va="top")
    ax.text(3.0, 29.0,
            "• Patient × condition summary; patient bootstrap (49,999×)\n"
            "• Cell-type literature scores; Welch DEG (q<0.05, |log2FC|≥1)",
            fontsize=5.1, color=PALETTE["gray_3"], va="top", linespacing=1.28)


def draw_gene_counts(ax, coh):
    d = pd.read_csv(P_DIR / f"{coh}_gene_counts_displayed.csv").sort_values("rank")
    m = pd.read_csv(P_DIR / "main_summary_numbers.csv")
    R = {s: m[(m.cohort == coh) & (m.source == s) & (m.metric == "log10_gene_mean_Pearson_r")].value.iloc[0] for s in ["PBISC", "scSemiProfiler"]}
    for src in ["scSemiProfiler", "PBISC", "Real"]:
        y = d[src].values.astype(float)
        ok = y > 0
        if src != "Real":
            ax.scatter(d["rank"][ok], y[ok], s=1.2, color=SRC_COL[src], alpha=0.34, lw=0, rasterized=True)
        med = pd.Series(np.where(ok, y, np.nan)).rolling(151, center=True, min_periods=30).median()
        lab = "Real" if src == "Real" else f"{'scSemi' if src=='scSemiProfiler' else src} (r={R[src]:.2f})"
        ax.plot(d["rank"], med, color=SRC_COL[src], lw=1.2 if src == "Real" else 1.0, label=lab, zorder=5 if src == "Real" else 4)
    ax.set_yscale("log"); ax.set_ylim(1e-4, 300); ax.set_xlim(0, 6000); ax.set_xticks([0, 3000, 6000])
    ax.set_xlabel("Gene rank in Real", fontsize=5.7); ax.set_ylabel("Mean UMI / cell", fontsize=5.7)
    h, l = ax.get_legend_handles_labels()
    order = [l.index("Real")] + [i for i, x in enumerate(l) if x != "Real"]
    ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper right", fontsize=5.3,
              handlelength=1.05, labelspacing=0.18)
    panel_title(ax, f"{NAME[coh]} · gene-level UMI", size=6.4, pad=3)


def _clean_umap(ax, o):
    x0, x1 = float(o.umap1.min()), float(o.umap1.max())
    y0, y1 = float(o.umap2.min()), float(o.umap2.max())
    xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    half = 0.54 * max(x1 - x0, y1 - y0)
    ax.set_xlim(xm - half, xm + half); ax.set_ylim(ym - half, ym + half)
    ax.set_box_aspect(1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.annotate("", xy=(0.20, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color="k", shrinkA=0, shrinkB=0, mutation_scale=5.5))
    ax.annotate("", xy=(0.02, 0.20), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color="k", shrinkA=0, shrinkB=0, mutation_scale=5.5))
    ax.text(0.11, -0.02, "UMAP1", transform=ax.transAxes, fontsize=5.3, ha="center", va="top")
    ax.text(-0.02, 0.11, "UMAP2", transform=ax.transAxes, fontsize=5.3, rotation=90, ha="right", va="center")


def draw_umap_source(ax, coh, ilisi_txt):
    full = pd.read_csv(F4P_DIR / f"umap_{coh}_post.csv.gz", index_col=0)
    idx = [g.sample(min(30000, len(g)), random_state=1).index for _, g in full.groupby("source")]
    o = full.loc[[i for ix in idx for i in ix]].sample(frac=1, random_state=2)
    col = np.where(o.source == "real_A", SRC_COL["Real"], SRC_COL["PBISC"])
    ax.scatter(o.umap1, o.umap2, c=col, s=0.5, alpha=0.45, lw=0, rasterized=True)
    _clean_umap(ax, full)
    ax.text(0.98, 0.98, ilisi_txt, transform=ax.transAxes, ha="right", va="top", fontsize=5.2,
            color=PALETTE["blue_main"], fontweight="bold",
            bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.4))
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=SRC_COL["Real"], label="Real A", ms=3.5),
        Line2D([], [], marker="o", ls="", color=SRC_COL["PBISC"], label="PBISC", ms=3.5),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, fontsize=5.4, handletextpad=0.16, columnspacing=0.6)
    panel_title(ax, f"{NAME[coh]} · source mixing", size=6.3, pad=2.5)


def draw_umap_six(ax, coh, show_legend=False):
    full = pd.read_csv(F4P_DIR / f"umap_{coh}_post.csv.gz", index_col=0)
    idx = [g.sample(min(30000, len(g)), random_state=1).index for _, g in full.groupby("source")]
    o = full.loc[[i for ix in idx for i in ix]].sample(frac=1, random_state=4)
    ax.scatter(o.umap1, o.umap2, c=[CT_COL.get(t, PALETTE["neutral"]) for t in o.six_type.astype(str)],
               s=0.5, alpha=0.60, lw=0, rasterized=True)
    _clean_umap(ax, full)
    if show_legend:
        ax.legend(handles=[Line2D([], [], marker="o", ls="", color=CT_COL[t], label=t.replace(" cell", ""), ms=3.4) for t in T6],
                  loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3, fontsize=5.2,
                  handletextpad=0.10, columnspacing=0.38, labelspacing=0.18)
    panel_title(ax, f"{NAME[coh]} · six lineages", size=6.3, pad=2.5)


def draw_patient_resampled_shift(ax, coh, df_shift6, show_legend=False):
    sub = df_shift6[df_shift6["cohort"] == coh]
    yy = np.arange(len(T6))
    for j, src in enumerate(SRC):
        s = sub[sub["source"] == src].set_index("cell_type").reindex(T6)
        yp = yy + (j - 1) * 0.26
        ax.barh(yp, s["shift_pp"].values, height=0.24, color=SRC_COL[src], edgecolor="k", lw=0.45,
                label=("scSemi" if src == "scSemiProfiler" else src))
        xerr_lo = np.maximum(0.0, s["shift_pp"].values - s["ci_lo"].values)
        xerr_hi = np.maximum(0.0, s["ci_hi"].values - s["shift_pp"].values)
        ax.errorbar(s["shift_pp"].values, yp, xerr=[xerr_lo, xerr_hi], fmt="none", ecolor="k", elinewidth=0.55, capsize=1.2, zorder=4)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(yy); ax.set_yticklabels(["T", "NK", "B", "Mono", "DC", "Other"], fontsize=5.6)
    ax.invert_yaxis(); ax.set_xlim(-22, 22); ax.set_xticks([-15, 0, 15])
    ax.set_xlabel("Case − ctrl (pp, pt boot.)", fontsize=5.6)
    if show_legend:
        ax.legend(loc="lower right", fontsize=5.2, handlelength=0.8, labelspacing=0.16)
    panel_title(ax, f"{NAME[coh]} · resampled Δ%", size=6.1, pad=3)


def draw_literature_scores(ax, coh, df_lit_sum, show_legend=False):
    sub = df_lit_sum[df_lit_sum["cohort"] == coh].copy()
    if coh == "sepsis":
        items = [
            ("Monocyte", "Myeloid S100/IL1R2", "Mono: S100/IL1R2"),
            ("Monocyte", "Myeloid MHC-II", "Mono: MHC-II"),
            ("T cell", "T NR4A/GPR183", "T: NR4A/GPR183"),
            ("B cell", "B NR4A/CD69", "B: NR4A/CD69"),
            ("NK cell", "NK Stress/PRF1", "NK: Stress/PRF1"),
        ]
        xlim = (-2.6, 4.6)
        xticks = [-2, 0, 2, 4]
        ttl = "Sepsis · Ye gene scores"
    else:
        items = [
            ("Monocyte", "Type-I IFN (7-ISG)", "Mono: 7-ISG"),
            ("T cell", "Type-I IFN (7-ISG)", "T: 7-ISG"),
            ("B cell", "Type-I IFN (7-ISG)", "B: 7-ISG"),
            ("NK cell", "Type-I IFN (7-ISG)", "NK: 7-ISG"),
        ]
        xlim = (-0.5, 4.9)
        xticks = [0, 2, 4]
        ttl = "COVID-19 · 7-ISG scores"

    yy = np.arange(len(items))
    for j, src in enumerate(SRC):
        vals, lo_e, hi_e = [], [], []
        for ct, sig, _ in items:
            r = sub[(sub["cell_type"] == ct) & (sub["signature"] == sig) & (sub["source"] == src)].iloc[0]
            vals.append(r["delta_score"])
            lo_e.append(max(0.0, r["delta_score"] - r["ci_lo"]))
            hi_e.append(max(0.0, r["ci_hi"] - r["delta_score"]))
        yp = yy + (j - 1) * 0.26
        ax.barh(yp, vals, height=0.24, color=SRC_COL[src], edgecolor="k", lw=0.45,
                label=("scSemi" if src == "scSemiProfiler" else src))
        ax.errorbar(vals, yp, xerr=[lo_e, hi_e], fmt="none", ecolor="k", elinewidth=0.55, capsize=1.2, zorder=4)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(yy); ax.set_yticklabels([x[2] for x in items], fontsize=5.4)
    ax.invert_yaxis(); ax.set_xlim(*xlim); ax.set_xticks(xticks)
    ax.set_xlabel("Δ score (log2 CPM+1)", fontsize=5.6)
    if show_legend:
        ax.legend(loc="lower right", fontsize=5.2, handlelength=0.8, labelspacing=0.16)
    panel_title(ax, ttl, size=6.2, pad=3)


def draw_welch_calibration(ax, coh):
    F = pd.read_csv(DATA_OUT / "fig5_welch_deg_fidelity.csv")
    C = pd.read_csv(DATA_OUT / "fig5_welch_deg_concordance.csv")
    C["source"] = C["source"].replace({"Ceiling": "PBISC"})
    xx = np.arange(4); w = 0.35
    for j, src in enumerate(["PBISC", "scSemiProfiler"]):
        s = F[(F.cohort == coh) & (F.source == src)].set_index("celltype").reindex(CT4)
        ax.bar(xx + (j - 0.5) * w, s.slope_corrected, w, color=SRC_COL[src], edgecolor="k", lw=0.6,
               label=f"{'scSemi' if src=='scSemiProfiler' else src} slope")
        ax.scatter(xx + (j - 0.5) * w, s.r_disattenuated, marker="D", s=14, color="white", edgecolor="k", lw=0.65,
                   zorder=4, label="Disattenuated r" if j == 0 else None)
    P_ = F[(F.cohort == coh) & (F.source == "PBISC")].set_index("celltype").AUROC
    S_ = F[(F.cohort == coh) & (F.source == "scSemiProfiler")].set_index("celltype").AUROC
    ax.axhline(1.0, color="k", lw=0.65, ls=(0, (3, 2)))
    ax.set_xticks(xx)
    ax.set_xticklabels([f"{c.replace(' cell', '').replace('Monocyte', 'Mono')}\nAUROC {P_[c]:.2f}·{S_[c]:.2f}" for c in CT4], fontsize=5.3)
    ax.set_ylabel("Effect size vs Real", fontsize=5.7); ax.set_ylim(0, 2.08)
    ax.legend(loc="upper right", fontsize=5.2, handlelength=0.85, labelspacing=0.15, ncol=3 if coh == "sepsis" else 1)
    panel_title(ax, f"{NAME[coh]} · Welch DEG calibration (BH q<0.05, |log2FC|≥1)", size=6.2, pad=3)


VOLC_UP = {
    ("sepsis", "T cell"): ["RPS26", "CD69", "CD83", "NR4A1", "SOCS3"],
    ("sepsis", "B cell"): ["RASD1", "RPS26", "CD83", "CD69", "NR4A1"],
    ("sepsis", "NK cell"): ["RPS26", "CD69", "CD83", "AREG", "NR4A1"],
    ("sepsis", "Monocyte"): ["IL1R2", "MCEMP1", "AQP9", "TNFAIP6", "ACSL1", "S100A8"],
    ("covid", "T cell"): ["IFI27", "IFIT1", "IFIT3", "IFI44L", "MX1", "OAS1"],
    ("covid", "B cell"): ["MX1", "OAS1", "OAS3", "IRF7", "USP18", "HERC5"],
    ("covid", "NK cell"): ["OAS1", "IRF7", "MX1", "IFIT1", "IFI44L", "IFI27"],
    ("covid", "Monocyte"): ["IFI27", "IFIT1", "MX1", "CD163", "MCEMP1", "IFI44L"],
}
VOLC_DN = {
    ("sepsis", "T cell"): ["CX3CR1", "EOMES"],
    ("sepsis", "B cell"): ["TLR10", "KLRB1"],
    ("sepsis", "NK cell"): ["ZNF683", "CX3CR1"],
    ("sepsis", "Monocyte"): ["CIITA", "HLA-DMB", "HLA-DRA"],
    ("covid", "T cell"): ["NR4A2", "JUNB"],
    ("covid", "B cell"): ["NR4A2", "CXCR5"],
    ("covid", "NK cell"): ["NR4A2", "ZNF331"],
    ("covid", "Monocyte"): ["VMO1", "CDKN1C"],
}


def draw_volcano_1x(ax, coh, ct, src, show_xl=False, show_yl=False, scsemi_note=None):
    d = pd.read_csv(P_DIR / f"deg_{coh}_{ct.replace(' ', '')}_{src}.csv.gz")
    q = -np.log10(d.q.clip(lower=1e-300)); yclip = 20.0
    deg = d.deg.astype(bool); qq = q.clip(upper=yclip)
    ax.scatter(d.logFC[~deg], qq[~deg], s=1.1, color=PALETTE["neutral"], lw=0, rasterized=True)
    ax.scatter(d.logFC[deg], qq[deg], s=2.1, color=SRC_COL[src], lw=0, rasterized=True)
    ax.axhline(-np.log10(0.05), color=PALETTE["gray_1"], lw=0.55, ls=(0, (2, 2)))
    ax.axvline(-1, color=PALETTE["gray_1"], lw=0.55, ls=(0, (2, 2)))
    ax.axvline(1, color=PALETTE["gray_1"], lw=0.55, ls=(0, (2, 2)))
    ax.set_xlim(-6, 6); ax.set_ylim(0, 26.8); ax.set_xticks([-4, 0, 4]); ax.set_yticks([0, 10, 20])
    if show_xl:
        ax.set_xlabel("log2 FC (case − ctrl)", fontsize=5.6)
    if show_yl:
        ax.set_ylabel("−log10 BH q", fontsize=5.6)
    extra = f"\nscSemi: {scsemi_note}" if scsemi_note is not None else ""
    ax.text(0.97, 0.96, f"DEG: {int(deg.sum())} ({int(d.n_case.iloc[0])}/{int(d.n_ctrl.iloc[0])}){extra}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.2,
            bbox=dict(fc="white", ec="none", alpha=0.84, pad=0.5), zorder=7)

    deg_df = d[deg].copy(); deg_df["qq"] = qq[deg]
    up_g = deg_df[(deg_df.logFC > 0) & deg_df.gene.isin(VOLC_UP[(coh, ct)])].sort_values("qq", ascending=False).head(3)
    dn_g = deg_df[(deg_df.logFC < 0) & deg_df.gene.isin(VOLC_DN[(coh, ct)])].sort_values("qq", ascending=False).head(2)
    for rank_i, (_, r) in enumerate(up_g.iterrows()):
        yt = 22.8 if rank_i == 0 else (17.5 - (rank_i - 1) * 4.2)
        xt = min(1.65, max(1.10, r.logFC * 0.60)) if rank_i == 0 else min(4.5, max(1.8, r.logFC + 0.45))
        ax.annotate(r.gene, xy=(r.logFC, r.qq), xytext=(xt, yt),
                    fontsize=5.2, fontstyle="italic", fontweight="bold", ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.2),
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="#444444"))
    for rank_i, (_, r) in enumerate(dn_g.iterrows()):
        yt = 21.2 - rank_i * 4.5
        ax.annotate(r.gene, xy=(r.logFC, r.qq), xytext=(max(-5.1, min(-1.8, r.logFC - 0.55)), yt),
                    fontsize=5.2, fontstyle="italic", fontweight="bold", ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.78, pad=0.2),
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="#444444"))
    sh_s = "scSemiProfiler" if src == "scSemiProfiler" else src
    panel_title(ax, f"{ct.replace(' cell', '')} · {sh_s}", size=6.0, pad=2.5)


def build_figure5():
    df_shift6 = pd.read_csv(DATA_OUT / "fig5_six_lineage_shift_bootstrap.csv")
    df_lit_sum = pd.read_csv(DATA_OUT / "fig5_literature_gene_scores_summary.csv")

    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.32), facecolor="white")

    # Row 1: a (cross-cohort workflow), b (Sepsis gene counts), c (COVID gene counts)
    ax_a = fig.add_axes([0.025, 0.858, 0.370, 0.125])
    draw_fig5_schematic(ax_a); fig_panel_label(fig, 0.008, 0.992, "a")

    ax_b = fig.add_axes([0.468, 0.862, 0.212, 0.116])
    draw_gene_counts(ax_b, "sepsis"); fig_panel_label(fig, 0.410, 0.992, "b")

    ax_c = fig.add_axes([0.758, 0.862, 0.212, 0.116])
    draw_gene_counts(ax_c, "covid"); fig_panel_label(fig, 0.700, 0.992, "c")

    # Row 2: d-g (Post-Harmony UMAPs + iLISI)
    u_specs = [
        ("d", lambda a: draw_umap_source(a, "sepsis", "iLISI 1.57\n(scSemi 1.00)")),
        ("e", lambda a: draw_umap_six(a, "sepsis", show_legend=True)),
        ("f", lambda a: draw_umap_source(a, "covid", "iLISI 1.60\n(scSemi 1.00)")),
        ("g", lambda a: draw_umap_six(a, "covid", show_legend=True)),
    ]
    for idx, (lbl, fn) in enumerate(u_specs):
        ax = fig.add_axes([0.050 + idx * 0.238, 0.678, 0.194, 0.126])
        fn(ax); panel_label(ax, lbl, x=-0.18, y=1.05)

    # Row 3: h, i (Patient-resampled 6-lineage shift) & j, k (Cell-type literature disease gene scores)
    ax_h = fig.add_axes([0.068, 0.492, 0.172, 0.115])
    draw_patient_resampled_shift(ax_h, "sepsis", df_shift6, show_legend=False); fig_panel_label(fig, 0.008, 0.625, "h")

    ax_i = fig.add_axes([0.298, 0.492, 0.172, 0.115])
    draw_patient_resampled_shift(ax_i, "covid", df_shift6, show_legend=True); fig_panel_label(fig, 0.245, 0.625, "i")

    ax_j = fig.add_axes([0.585, 0.492, 0.158, 0.115])
    draw_literature_scores(ax_j, "sepsis", df_lit_sum, show_legend=False); fig_panel_label(fig, 0.488, 0.625, "j")

    ax_k = fig.add_axes([0.825, 0.492, 0.152, 0.115])
    draw_literature_scores(ax_k, "covid", df_lit_sum, show_legend=True); fig_panel_label(fig, 0.752, 0.625, "k")

    # Row 4: l, m (Welch's t-test DEG calibration & AUROC across T, B, NK, Monocyte)
    ax_l = fig.add_axes([0.080, 0.285, 0.400, 0.132])
    draw_welch_calibration(ax_l, "sepsis"); fig_panel_label(fig, 0.008, 0.436, "l")

    ax_m = fig.add_axes([0.568, 0.285, 0.405, 0.132])
    draw_welch_calibration(ax_m, "covid"); fig_panel_label(fig, 0.502, 0.436, "m")

    # Row 5: n (Sepsis Monocyte Welch DEG Real vs PBISC) & o (COVID T-cell Welch DEG Real vs PBISC)
    sc_sep_mono = int(pd.read_csv(P_DIR / "deg_sepsis_Monocyte_scSemiProfiler.csv.gz")["deg"].sum())
    sc_cov_t = int(pd.read_csv(P_DIR / "deg_covid_Tcell_scSemiProfiler.csv.gz")["deg"].sum())

    fig_panel_label(fig, 0.008, 0.218, "n")
    fig.text(0.036, 0.215, "Sepsis · Monocyte Welch DEGs", fontsize=6.2, fontweight="bold", color=PALETTE["gray_3"], va="top")
    ax_n1 = fig.add_axes([0.072, 0.052, 0.188, 0.132])
    draw_volcano_1x(ax_n1, "sepsis", "Monocyte", "Real", show_xl=True, show_yl=True)
    ax_n2 = fig.add_axes([0.292, 0.052, 0.188, 0.132])
    draw_volcano_1x(ax_n2, "sepsis", "Monocyte", "PBISC", show_xl=True, show_yl=False, scsemi_note=f"{sc_sep_mono} DEGs")

    fig_panel_label(fig, 0.502, 0.218, "o")
    fig.text(0.530, 0.215, "COVID-19 · T-cell Welch DEGs", fontsize=6.2, fontweight="bold", color=PALETTE["gray_3"], va="top")
    ax_o1 = fig.add_axes([0.566, 0.052, 0.188, 0.132])
    draw_volcano_1x(ax_o1, "covid", "T cell", "Real", show_xl=True, show_yl=True)
    ax_o2 = fig.add_axes([0.786, 0.052, 0.188, 0.132])
    draw_volcano_1x(ax_o2, "covid", "T cell", "PBISC", show_xl=True, show_yl=False, scsemi_note=f"{sc_cov_t} DEGs")

    return finalize_figure(fig, FIG_OUT / "Figure5_CrossCohort_BulkOnly_Comparison")


# =============================================================================
# MAIN FIGURE 6 HELPERS — RARE CELL FREQUENCY AND SAMPLE-SUMMED EXPRESSION
# =============================================================================
PROLIF_GROUPS = [
    "Ped. HC (4)", "Ped. ICU (4)", "Ped. abd. (7)", "Ped. pulm. (8)",
    "Adult HC (7)", "Adult ICU (6)", "Adult pulm. (18)", "Adult oth. (28)",
]
PROLIF_DOTS = [
    [0.00, 0.05, 0.18, 0.69], [0.18, 0.25, 0.55, 0.61],
    [0.00, 0.03, 0.13, 0.63, 6.71, 8.93, 12.37], [0.00, 0.00, 0.05, 0.34, 0.66, 1.29, 5.23, 11.13],
    [0.00, 0.00, 0.00, 0.00, 0.03, 0.08, 0.13], [0.00, 0.00, 0.00, 0.03, 0.08, 0.13],
    [0.00]*8 + [0.04, 0.08, 0.15, 0.28, 0.36, 0.44, 0.53, 0.74, 0.92, 2.01],
    [0.00]*15 + [0.04, 0.06, 0.08, 0.12, 0.18, 0.26, 0.31, 0.37, 0.44, 0.53, 0.66, 1.06, 3.33],
]
PROLIF_SQRS = [
    [0.00, 0.00, 0.00, 0.03], [0.00, 0.00, 0.00, 0.03],
    [0.00, 0.00, 0.00, 0.03, 2.93, 3.14, 5.50], [0.00, 0.00, 0.00, 0.00, 0.05, 0.29, 1.29, 2.09],
    [0.00]*7, [0.00]*5 + [0.08],
    [0.00]*12 + [0.05, 0.12, 0.23, 0.42, 0.66, 0.87],
    [0.00]*22 + [0.04, 0.06, 0.08, 0.14, 0.18, 1.19],
]


def draw_rare_parent_shares(ax, coh, show_legend=False):
    r = pd.read_csv(DATA_OUT / "fig6_rare_parent_shares_14sepsis_13covid.csv")
    r = r[r.cohort == coh]
    subs = list(dict.fromkeys(r[r.source == "Real"].subset))
    short_s = {
        "Mono_c06-CD14-MKI67": "Mono_c06-MKI67",
        "Mono_c01-CD14-HLA-DRA": "Mono_c01-HLA-DRA",
        "NKT_c05-IFNG-CCL4hi": "NKT_IFNG-CCL4hi",
        "HSC": "HSC (/all cells)",
        "HSPC": "HSPC (/all cells)",
    }
    y = np.arange(len(subs)); h = 0.26
    for j, src in enumerate(SRC):
        v = r[r.source == src].set_index("subset").percent_parent.reindex(subs).values
        vv = np.where(v > 0, v, np.nan)
        ax.barh(y + (j - 1) * h, vv, h, left=0.004, color=SRC_COL[src], edgecolor="k", lw=0.35,
                label=("scSemi" if src == "scSemiProfiler" else src))
        for yi, vi in zip(y + (j - 1) * h, v):
            if not vi > 0:
                ax.text(0.0045, yi, "0", va="center", fontsize=5.2, color=SRC_COL[src], fontweight="bold")
    ax.set_xscale("log"); ax.set_xlim(0.004, 120); ax.set_yticks(y)
    ax.set_yticklabels([short_s.get(s, s) for s in subs], fontsize=5.3)
    ax.invert_yaxis(); ax.set_xlabel("% of parent (cases)", fontsize=5.7)
    if show_legend:
        ax.legend(loc="lower right", fontsize=5.2, handlelength=0.8, labelspacing=0.16)
    n_st = len(subs)
    panel_title(ax, f"{NAME[coh]} · {n_st} rare states", size=6.3, pad=3)


def draw_sepsis_prolif_compact(ax):
    for xg in [0, 5, 10]:
        ax.axvline(xg, color="#EAEAEA", lw=0.5, zorder=0)
    for i in range(len(PROLIF_GROUPS)):
        vd = np.asarray(PROLIF_DOTS[i], float); vs = np.asarray(PROLIF_SQRS[i], float)
        ax.scatter(vd, np.full_like(vd, i - 0.15), s=7.5, marker="o", facecolor=PALETTE["blue_main"], edgecolor="white", lw=0.25, zorder=3)
        ax.scatter(vs, np.full_like(vs, i + 0.16), s=7.5, marker="s", facecolor="white", edgecolor=PALETTE["gray_1"], lw=0.55, zorder=2)
    ax.set_yticks(np.arange(len(PROLIF_GROUPS))); ax.set_yticklabels(PROLIF_GROUPS, fontsize=5.3)
    ax.set_ylim(len(PROLIF_GROUPS) - 0.4, -0.6); ax.set_xlim(-0.5, 13.2); ax.set_xticks([0, 5, 10])
    ax.set_xlabel("Co-detected mono. (%)", fontsize=5.6)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=PALETTE["blue_main"], mec="white", ms=3.4, label="MKI67+STMN1"),
        Line2D([], [], marker="s", ls="", mfc="white", mec=PALETTE["gray_1"], mew=0.6, ms=3.4, label="+MPO/ELANE"),
    ], loc="lower right", fontsize=5.1, handletextpad=0.15, labelspacing=0.15)
    panel_title(ax, "Sepsis 82d · mono.", size=6.1, pad=3)


RARE_UMAP = {"sepsis": ["B_c07-MZB1", "Macrophage", "pDC", "B_c08-MKI67"], "covid": ["pDC", "DC1", "B_ASC_dividing", "MAIT"]}
RARE_UMAP_COL = [PALETTE["violet"], PALETTE["teal"], PALETTE["highlight"], PALETTE["magenta"]]


def draw_umap_rare_1x(ax, coh):
    o = pd.read_csv(F4P_DIR / f"umap_{coh}_post.csv.gz", index_col=0)
    bg = o.sample(60000, random_state=3)
    ax.scatter(bg.umap1, bg.umap2, c="#E6E6E6", s=0.42, lw=0, rasterized=True)
    h = []
    for st, colr in zip(RARE_UMAP[coh], RARE_UMAP_COL):
        r = o[(o.fine.astype(str) == st) & (o.source == "real_A")]
        g = o[(o.fine.astype(str) == st) & (o.source == "ceiling")]
        ax.scatter(r.umap1, r.umap2, s=3.8, color=colr, lw=0, alpha=0.82, rasterized=True)
        ax.scatter(g.umap1, g.umap2, s=11.0, marker="^", facecolor=colr, edgecolor="k", lw=0.35, zorder=5, rasterized=True)
        h.append(Patch(color=colr, label=f"{st.replace('B_ASC_dividing','ASC_div')} ({len(r)}/{len(g)})"))
    h += [
        Line2D([], [], marker="o", ls="", color=PALETTE["gray_1"], label="Real A", ms=3.0),
        Line2D([], [], marker="^", ls="", mfc=PALETTE["gray_1"], mec="k", mew=0.4, label="PBISC", ms=3.5),
    ]
    _clean_umap(ax, o)
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, fontsize=5.1,
              handlelength=0.75, columnspacing=0.35, labelspacing=0.15)
    panel_title(ax, f"{NAME[coh]} · rare UMAP", size=6.2, pad=2.5)


def draw_rare_patient_1x(ax, coh, sub, st_df, short_title=None):
    u = pd.read_csv(P_DIR / "rare_patient_units.csv.gz")
    u = u[(u.cohort == coh) & (u.subset == sub)]
    st = st_df[(st_df.cohort == coh) & (st_df.subset == sub)]
    conds = [("control", "HC"), ("ICU control", "ICU"), ("case", "Sepsis")] if coh == "sepsis" else [("control", "Neg."), ("case", "Pos.")]
    rng = np.random.default_rng(20261111 + (zlib.crc32(f"{coh}:{sub}".encode()) % 100000))

    if (coh, sub) == ("covid", "CD14mono_dex"):
        fig = ax.figure; bbox = ax.get_position(); ax.axis("off")
        x0, y0, w_tot, h_tot = bbox.x0, bbox.y0, bbox.width, bbox.height
        ax_bot = fig.add_axes([x0, y0, w_tot, 0.62 * h_tot])
        ax_top = fig.add_axes([x0, y0 + 0.69 * h_tot, w_tot, 0.31 * h_tot], sharex=ax_bot)
        y_bot_max, y_top_min, y_top_max = 0.85, 1.5, 13.9
        for i, (cond, _) in enumerate(conds):
            xr = i - 0.27
            for j, src in enumerate(SRC):
                v = u[(u.group == cond) & (u.source == src)].pct_all.values
                xp = i + (j - 1) * 0.27
                if not len(v):
                    continue
                jit = rng.uniform(-0.07, 0.07, len(v))
                mb, mt = v <= y_bot_max, v >= y_top_min
                if mb.any(): ax_bot.scatter(xp + jit[mb], v[mb], s=6, color=SRC_COL[src], alpha=0.78, lw=0)
                if mt.any(): ax_top.scatter(xp + jit[mt], v[mt], s=6, color=SRC_COL[src], alpha=0.78, lw=0)
                m, sd = float(v.mean()), (float(v.std(ddof=1)) if len(v) > 1 else 0.0)
                if m <= y_bot_max:
                    ax_bot.plot([xp - 0.06, xp + 0.06], [m, m], color="k", lw=0.9)
                    ax_bot.plot([xp, xp], [max(0.0, m - sd), min(m + sd, y_bot_max)], color="k", lw=0.75)
                if src != "Real":
                    s = st[(st.condition == cond) & (st.source == src)]
                    if len(s):
                        yb = 8.0 if src == "PBISC" else 11.2
                        ax_top.plot([xr, xr, xp, xp], [yb, yb + 0.5, yb + 0.5, yb], color=SRC_COL[src], lw=0.6)
                        ax_top.text(0.5 * (xr + xp), yb + 0.6, s.stars.iloc[0], ha="center", va="bottom", fontsize=5.4, color=SRC_COL[src])
        ax_bot.set_xlim(-0.55, len(conds) - 0.45); ax_bot.set_ylim(-0.02, y_bot_max); ax_bot.set_yticks([0.0, 0.4, 0.8])
        ax_bot.set_xticks(range(len(conds))); ax_bot.set_xticklabels([c[1] for c in conds], fontsize=5.7)
        ax_bot.set_ylabel("% cells / patient", y=0.80, fontsize=5.7)
        ax_top.set_ylim(y_top_min, y_top_max); ax_top.set_yticks([3, 7]); ax_top.spines["bottom"].set_visible(False)
        ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        panel_title(ax_top, short_title or f"{NAME[coh]} · {sub}", size=6.3, pad=3)
        return

    ymax = max(float(u.pct_all.max()), 0.01)
    for cond, _ in conds:
        for src in SRC:
            v = u[(u.group == cond) & (u.source == src)].pct_all.values
            if len(v) > 1:
                ymax = max(ymax, float(v.mean() + v.std(ddof=1)))
    hh = ymax * 0.035
    for i, (cond, _) in enumerate(conds):
        xr = i - 0.27
        for j, src in enumerate(SRC):
            v = u[(u.group == cond) & (u.source == src)].pct_all.values
            x0 = i + (j - 1) * 0.27
            if not len(v):
                continue
            ax.scatter(x0 + rng.uniform(-0.07, 0.07, len(v)), v, s=5.8, color=SRC_COL[src], alpha=0.75, lw=0, zorder=3, rasterized=True)
            ax.errorbar(x0, v.mean(), yerr=v.std(ddof=1) if len(v) > 1 else 0, fmt="_", color="k", ms=4.8, capsize=1.6, lw=0.75, zorder=4)
            if src != "Real":
                s = st[(st.condition == cond) & (st.source == src)]
                if len(s):
                    yb = ymax * (1.06 if src == "PBISC" else 1.25)
                    ax.plot([xr, xr, x0, x0], [yb, yb + hh, yb + hh, yb], color=SRC_COL[src], lw=0.58, zorder=5)
                    ax.text(0.5 * (xr + x0), yb + hh + ymax * 0.01, s.stars.iloc[0], ha="center", va="bottom",
                            fontsize=5.4, color=SRC_COL[src], fontweight="bold" if s.stars.iloc[0] != "ns" else "normal")
    ax.set_xticks(range(len(conds))); ax.set_xticklabels([c[1] for c in conds], fontsize=5.7)
    ax.set_ylim(0, ymax * 1.48); ax.set_ylabel("% cells / patient", fontsize=5.7)
    panel_title(ax, short_title or f"{NAME[coh]} · {sub}", size=6.2, pad=3)


def draw_summed_interpatient_similarity(fig, bbox, coh, cdata, letter):
    mat_4x5 = cdata["mat_4x5"]
    pair_dists = cdata["pair_dists"]
    focal4 = cdata["focal4"]
    x0, y0, w_box, h_box = bbox
    fig_panel_label(fig, x0, y0 + h_box * 0.99, letter)
    fig.text(x0 + 0.030, y0 + h_box * 0.985, f"{NAME[coh]} · summed sample inter-patient Pearson r",
             fontsize=6.3, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")

    ax_m = fig.add_axes([x0 + 0.235 * w_box, y0 + 0.21 * h_box, 0.335 * w_box, 0.61 * h_box])
    ax_p = fig.add_axes([x0 + 0.610 * w_box, y0 + 0.21 * h_box, 0.280 * w_box, 0.61 * h_box])

    norm = TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=1.0)
    masked = np.ma.masked_invalid(np.clip(mat_4x5, -0.5, 1.0))
    cmap_m = CMAP_CORR.copy(); cmap_m.set_bad("#EAEAEA")
    ax_m.imshow(masked, cmap=cmap_m, norm=norm, aspect="auto", interpolation="nearest")
    ax_m.axvline(3.5, color=PALETTE["gray_3"], lw=0.8)
    for i in range(4):
        for j in range(5):
            v = mat_4x5[i, j]
            if np.isfinite(v):
                tc = "white" if (v >= 0.62 or v <= -0.35) else PALETTE["gray_3"]
                ax_m.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.2,
                          fontweight="bold" if (i == j or j == 4) else "normal", color=tc)
            else:
                ax_m.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor="#EFEFEF", edgecolor="#CCCCCC", hatch="///", lw=0.4))
                ax_m.text(j, i, "0", ha="center", va="center", fontsize=5.2, color=PALETTE["gray_2"], fontweight="bold")
    sh_map = {
        "B_c07-MZB1": "B_c07-MZB1", "Macrophage": "Macrophage",
        "T_CD4_c07-MKI67": "CD4_MKI67", "Mono_c06-CD14-MKI67": "Mono_MKI67",
        "B_ASC_IgG": "B_ASC_IgG", "CD14mono_dex": "CD14mono_dex", "pDC": "pDC", "HSPC": "HSPC",
    }
    sh_lbls = [sh_map[s] for s in focal4]
    ax_m.set_yticks(np.arange(4)); ax_m.set_yticklabels(sh_lbls, fontsize=5.3, fontweight="bold")
    ax_m.set_xticks(np.arange(5)); ax_m.set_xticklabels(sh_lbls + ["scSemi"], fontsize=5.2, rotation=32, ha="right")
    ax_m.tick_params(axis="both", length=0, pad=1.5)
    for sp_ in ax_m.spines.values():
        sp_.set_visible(True); sp_.set_linewidth(0.6); sp_.set_edgecolor(PALETTE["gray_3"])

    ax_p.axvline(0.0, color=PALETTE["neutral"], lw=0.65, ls=(0, (2, 2)))
    for i, sub in enumerate(focal4):
        st = pair_dists[sub]
        for key, yy, col, mk in [("rr_all", i - 0.22, SRC_COL["Real"], "D"), ("rp_all", i + 0.02, SRC_COL["PBISC"], "o"), ("off_sig", i + 0.24, PALETTE["gray_1"], "s")]:
            arr = st[key]; arr = arr[np.isfinite(arr)]
            if len(arr) == 0:
                continue
            q10, q25, med, q75, q90 = np.quantile(arr, [0.10, 0.25, 0.50, 0.75, 0.90])
            ax_p.plot([q10, q90], [yy, yy], color=col, lw=0.65, alpha=0.8)
            ax_p.plot([q25, q75], [yy, yy], color=col, lw=1.45)
            ax_p.scatter([med], [yy], marker=mk, s=10.5, facecolor="white" if key == "off_sig" else col, edgecolor=col, lw=0.6, zorder=4)
        mrp = float(np.nanmedian(st["rp_all"])); mrr = float(np.nanmedian(st["rr_all"]))
        ax_p.text(1.06, i - 0.12, f"{mrp:.2f}", fontsize=5.2, fontweight="bold", color=SRC_COL["PBISC"], va="center")
        ax_p.text(1.06, i + 0.16, f"({mrr:.2f})", fontsize=5.2, color=PALETTE["gray_2"], va="center")
    ax_p.set_ylim(3.55, -0.55); ax_p.set_xlim(-0.42, 1.04); ax_p.set_yticks([])
    ax_p.set_xticks([0.0, 0.5, 1.0]); ax_p.spines["left"].set_visible(False)
    ax_p.set_xlabel("Inter-patient summed r", fontsize=5.4)
    ax_p.legend(handles=[
        Line2D([], [], marker="D", ls="", color=SRC_COL["Real"], ms=3.1, label="Real–Real"),
        Line2D([], [], marker="o", ls="", color=SRC_COL["PBISC"], ms=3.1, label="Real–PBISC"),
        Line2D([], [], marker="s", ls="", mfc="white", mec=PALETTE["gray_1"], ms=3.1, label="Other state"),
    ], loc="upper right", bbox_to_anchor=(1.22, 1.21), ncol=3, fontsize=5.1, handletextpad=0.10, columnspacing=0.28)


def draw_summed_marker_heatmap(fig, bbox, cdata, sub, title_str, letter=None):
    mdf_all = cdata["mdf"]
    Z_all = cdata["Z_mk"]
    C_all = cdata["C_all"]
    mk_g2k = cdata["mk_g2k"]
    mods = cdata["marker_mods"][sub]

    sub_idx = np.flatnonzero((mdf_all["subset"] == sub).to_numpy())
    # Order samples within each source by mean correlation to all samples in the same state
    ordered_sub = []
    for sn in ("Real", "PBISC", "scSemiProfiler"):
        ix = np.flatnonzero((mdf_all["subset"] == sub) & (mdf_all["source"] == sn))
        if len(ix) > 1:
            mr = C_all[np.ix_(ix, sub_idx)].mean(axis=1)
            ix = ix[np.argsort(-mr)]
        ordered_sub.extend(ix.tolist())
    ordered_sub = np.array(ordered_sub, dtype=int)
    mdf = mdf_all.iloc[ordered_sub].reset_index(drop=True)

    col_indices, marker_meta = [], []
    for m_name, m_genes in mods:
        mc = MOD_COLOR_MAP.get(m_name, PALETTE["blue_secondary"])
        for g in m_genes:
            col_indices.append(mk_g2k[g])
            marker_meta.append((m_name, mc, g))
    Z_sub = Z_all[np.ix_(ordered_sub, np.array(col_indices, int))]
    N, G = len(mdf), len(marker_meta)
    nsc = int((mdf["source"] == "scSemiProfiler").sum())

    x0, y0, w_box, h_box = bbox
    if letter is not None:
        fig_panel_label(fig, x0, y0 + h_box * 0.99, letter)
    n_tag = f"n={N}; scSemi: 0" if nsc == 0 else f"n={N}"
    fig.text(x0 + (0.030 if letter else 0.01), y0 + h_box * 0.985, f"{title_str} ({n_tag})",
             fontsize=6.1, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")

    mx = x0 + 0.285 * w_box; my = y0 + 0.08 * h_box; mw = 0.610 * w_box; mh = 0.64 * h_box
    st_h = 0.045 * h_box; gy = 0.015 * h_box; mod_w = 0.022 * w_box; gx = 0.008 * w_box
    ax_hm = fig.add_axes([mx, my, mw, mh])
    ax_src = fig.add_axes([mx, my + mh + gy, mw, st_h])
    ax_mod = fig.add_axes([mx - gx - mod_w, my, mod_w, mh])
    ax_cb = fig.add_axes([mx + mw + 0.02 * w_box, my + 0.12 * mh, 0.022 * w_box, 0.76 * mh])

    im = ax_hm.imshow(np.clip(Z_sub.T, -1.2, 2.2), cmap=CMAP_EXPR, vmin=-1.2, vmax=2.2,
                      aspect="auto", interpolation="nearest", rasterized=True)
    ax_hm.set_xticks([]); ax_hm.set_yticks([])
    for sp_ in ax_hm.spines.values():
        sp_.set_visible(True); sp_.set_linewidth(0.6); sp_.set_edgecolor(PALETTE["gray_3"])

    src_rgb = np.array([hex_to_rgb(SRC_COL[s]) for s in mdf["source"]], float)
    mod_rgb = np.array([hex_to_rgb(mc) for _, mc, _ in marker_meta], float)
    ax_src.imshow(src_rgb[None, :, :], aspect="auto", interpolation="nearest", rasterized=True)
    ax_mod.imshow(mod_rgb[:, None, :], aspect="auto", interpolation="nearest", rasterized=True)
    for a_s in (ax_src, ax_mod):
        a_s.set_xticks([]); a_s.set_yticks([])
        for sp_ in a_s.spines.values():
            sp_.set_visible(True); sp_.set_linewidth(0.45); sp_.set_edgecolor(PALETTE["gray_2"])

    present = []
    for sn in ("Real", "PBISC", "scSemiProfiler"):
        ix = np.flatnonzero((mdf["source"] == sn).to_numpy())
        if len(ix):
            present.append((sn, int(ix[0]), int(ix[-1]), len(ix)))
    for b_i, (_, _, end, _) in enumerate(present):
        if b_i < len(present) - 1:
            ax_hm.axvline(end + 0.5, color="white", lw=0.9, zorder=3)
            ax_hm.axvline(end + 0.5, color=PALETTE["gray_3"], lw=0.55, ls="--", zorder=4)
            ax_src.axvline(end + 0.5, color="white", lw=0.7, zorder=4)

    raw_xf = [(0.5 * (s + e) + 0.5) / N for _, s, e, _ in present]
    adj_xf = _resolve_1d_positions(raw_xf, min_sep=0.37, lo=0.12, hi=0.90)
    y_ts = my + mh + gy + st_h
    for (sn, _, _, cnt), x_lbl in zip(present, adj_xf):
        sh_n = "scSemi" if sn == "scSemiProfiler" else sn
        fig.text(mx + x_lbl * mw, y_ts + 0.015 * h_box, f"{sh_n} ({cnt})",
                 fontsize=5.2, fontweight="bold", color=SRC_COL[sn], ha="center", va="bottom")

    x_ml = mx - gx - mod_w
    for k, (_, _, gene) in enumerate(marker_meta):
        yr = my + (1.0 - (k + 0.5) / G) * mh
        fig.text(x_ml - 0.01 * w_box, yr, gene, fontsize=5.2, fontstyle="italic", fontweight="bold",
                 color=PALETTE["gray_3"], ha="right", va="center")
    cur = 0
    for m_name, mgenes in mods:
        end_r = cur + len(mgenes) - 1
        ym = my + (1.0 - (0.5 * (cur + end_r) + 0.5) / G) * mh
        fig.text(x_ml - 0.14 * w_box, ym, m_name.replace(" response", "\nresp.").replace(" cycle", "\ncycle"),
                 fontsize=5.2, fontweight="bold", color=PALETTE["gray_2"], ha="right", va="center")
        if end_r < G - 1:
            ax_hm.axhline(end_r + 0.5, color="white", lw=0.85, zorder=3)
            ax_hm.axhline(end_r + 0.5, color=PALETTE["gray_3"], lw=0.5, zorder=4)
        cur = end_r + 1

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_ticks([-1, 0, 1, 2]); cb.ax.tick_params(labelsize=5.2, width=0.5, length=1.6, pad=1)
    cb.outline.set_linewidth(0.5); cb.set_label("Z", fontsize=5.3, labelpad=1.5)


def build_figure6(st_df, scache):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.34), facecolor="white")

    # Row 1: a (Sepsis 14 rare shares), b (COVID 13 rare shares), c (Sepsis 82d cycling monocyte audit)
    ax_a = fig.add_axes([0.148, 0.768, 0.228, 0.195])
    draw_rare_parent_shares(ax_a, "sepsis", show_legend=False); fig_panel_label(fig, 0.008, 0.990, "a")

    ax_b = fig.add_axes([0.522, 0.768, 0.212, 0.195])
    draw_rare_parent_shares(ax_b, "covid", show_legend=True); fig_panel_label(fig, 0.395, 0.990, "b")

    ax_c = fig.add_axes([0.840, 0.768, 0.135, 0.195])
    draw_sepsis_prolif_compact(ax_c); fig_panel_label(fig, 0.742, 0.990, "c")

    # Row 2: d, e (Sepsis & COVID rare UMAPs) + f, g (Patient-level frequencies)
    ax_d = fig.add_axes([0.042, 0.558, 0.175, 0.132])
    draw_umap_rare_1x(ax_d, "sepsis"); fig_panel_label(fig, 0.008, 0.710, "d")

    ax_e = fig.add_axes([0.270, 0.558, 0.175, 0.132])
    draw_umap_rare_1x(ax_e, "covid"); fig_panel_label(fig, 0.238, 0.710, "e")

    ax_f = fig.add_axes([0.525, 0.542, 0.195, 0.148])
    draw_rare_patient_1x(ax_f, "sepsis", "Mono_c06-CD14-MKI67", st_df, short_title="Sepsis · Mono_c06-MKI67")
    fig_panel_label(fig, 0.468, 0.710, "f")

    ax_g = fig.add_axes([0.792, 0.542, 0.185, 0.148])
    draw_rare_patient_1x(ax_g, "covid", "CD14mono_IFN", st_df, short_title="COVID-19 · CD14mono_IFN")
    fig_panel_label(fig, 0.736, 0.710, "g")

    # Row 3: h, i (Summed sample inter-patient Pearson r & 2000-HVG specificity)
    draw_summed_interpatient_similarity(fig, [0.012, 0.278, 0.472, 0.182], "sepsis", scache["sepsis"], "h")
    draw_summed_interpatient_similarity(fig, [0.505, 0.278, 0.478, 0.182], "covid", scache["covid"], "i")

    # Row 4: j, k (Per-sample summed canonical marker expression heatmaps)
    draw_summed_marker_heatmap(fig, [0.015, 0.032, 0.468, 0.205], scache["sepsis"], "Mono_c06-CD14-MKI67",
                               "Sepsis · Mono_c06-MKI67 markers", letter="j")
    draw_summed_marker_heatmap(fig, [0.506, 0.032, 0.474, 0.205], scache["covid"], "CD14mono_IFN",
                               "COVID-19 · CD14mono_IFN markers", letter="k")

    return finalize_figure(fig, FIG_OUT / "Figure6_Rare_Cell_Frequency_and_Expression")


# =============================================================================
# SUPPLEMENTARY FIGURES 1-5
# =============================================================================
def build_suppfig_coverage(snum=14):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.25), facecolor="white")
    gs = fig.add_gridspec(5, 2, hspace=0.62, wspace=0.40, left=0.16, right=0.96, top=0.96, bottom=0.05)

    cover_spec = {
        "sepsis": [("Genes in Real file", 37474, SRC_COL["Real"]), ("PBISC panel", 6000, SRC_COL["PBISC"]),
                   ("Cohort bulk genes", 5975, PALETTE["blue_secondary"]), ("Panel in Real", 5967, PALETTE["gray_1"]),
                   ("scSemi genes", 5961, SRC_COL["scSemiProfiler"]), ("Harmony HVGs", 2000, PALETTE["neutral"])],
        "covid": [("Genes in Real file", 33559, SRC_COL["Real"]), ("PBISC panel", 6000, SRC_COL["PBISC"]),
                  ("Cohort bulk genes", 5927, PALETTE["blue_secondary"]), ("Panel in Real", 5791, PALETTE["gray_1"]),
                  ("scSemi genes", 5745, SRC_COL["scSemiProfiler"]), ("Harmony HVGs", 2000, PALETTE["neutral"])],
    }
    s_scope = pd.read_csv(DATA_OUT / "supp1_source_scope_and_depth.csv")
    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[0, j])
        rows = cover_spec[coh]; y = np.arange(len(rows))
        ax.barh(y, [r[1] for r in rows], color=[r[2] for r in rows], edgecolor="k", lw=0.5, height=0.65)
        for yi, r in zip(y, rows):
            ax.text(r[1] * 1.10, yi, f"{r[1]:,}", va="center", fontsize=5.4)
        ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=5.5)
        ax.invert_yaxis(); ax.set_xscale("log"); ax.set_xlim(1000, 2.5e5); ax.set_xlabel("Genes (log10 scale)", fontsize=5.6)
        panel_title(ax, f"{NAME[coh]} · gene coverage", size=6.3); panel_label(ax, L[j], x=-0.32, y=1.06)

    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[1, j])
        s = s_scope[s_scope.cohort == coh].set_index("source").reindex(SRC); y = np.arange(3)
        ax.barh(y, s.n_cells / 1e6, color=[SRC_COL[x] for x in SRC], edgecolor="k", lw=0.5, height=0.58)
        for yi, (_, r) in zip(y, s.iterrows()):
            ax.text(r.n_cells / 1e6 + 0.015, yi, f"{int(r.n_cells):,} ({int(r.n_samples)} samples)", va="center", fontsize=5.3)
        ax.set_yticks(y); ax.set_yticklabels(["Real", "PBISC", "scSemi"], fontsize=5.6)
        ax.invert_yaxis(); ax.set_xlim(0, 1.30); ax.set_xlabel("Cells (millions)", fontsize=5.6)
        panel_title(ax, f"{NAME[coh]} · cells & samples", size=6.3); panel_label(ax, L[2 + j], x=-0.32, y=1.06)

    bins = np.logspace(1.5, 4.7, 90); c_bin = 0.5 * (bins[1:] + bins[:-1])
    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[2, j]); txt = []
        for src in SRC:
            u = np.load(P_DIR / f"{coh}_{src}_all_cell_panel_umi.npy").astype(float)
            h, _ = np.histogram(u[u > 0], bins=bins)
            ax.plot(c_bin, 100 * h / h.sum(), color=SRC_COL[src], lw=1.05, label=("scSemi" if src == "scSemiProfiler" else src))
            txt.append(f"{'scSemi' if src=='scSemiProfiler' else src}: {np.median(u):,.0f}")
        ax.set_xscale("log"); ax.set_xlabel("Panel UMI / cell", fontsize=5.6); ax.set_ylabel("Cells / bin (%)", fontsize=5.6)
        ax.legend(loc="upper left", fontsize=5.2, handlelength=0.9)
        ax.text(0.98, 0.96, "Median UMI\n" + "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=5.2)
        panel_title(ax, f"{NAME[coh]} · UMI depth density", size=6.3); panel_label(ax, L[4 + j], x=-0.32, y=1.06)

    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[3, j])
        for src in SRC:
            u = np.sort(np.load(P_DIR / f"{coh}_{src}_all_cell_panel_umi.npy").astype(float))
            step = max(1, len(u) // 4000)
            idx = np.unique(np.r_[np.arange(0, len(u), step), len(u) - 1])
            ax.plot(u[idx], 100 * (idx + 1) / len(u), color=SRC_COL[src], lw=1.05, label=("scSemi" if src == "scSemiProfiler" else src))
        ax.set_xscale("log"); ax.set_xlim(30, 5e4); ax.set_ylim(0, 101)
        ax.set_xlabel("Panel UMI / cell", fontsize=5.6); ax.set_ylabel("Cum. cells (%)", fontsize=5.6)
        ax.legend(loc="upper left", fontsize=5.2, handlelength=0.9)
        panel_title(ax, f"{NAME[coh]} · UMI depth ECDF", size=6.3); panel_label(ax, L[6 + j], x=-0.32, y=1.06)

    # Row 5 (i, j): Shift decomposition (total, Δ share, Δ expr.) + iLISI
    F_dec = pd.read_csv(DATA_OUT / "fig5_shift_decomposition.csv")
    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[4, j])
        comps = ["total", "composition", "within"]
        alphas = [1.0, 0.62, 0.32]; hatches = ["", "//", ".."]
        pos = []
        for i_s, src in enumerate(["PBISC", "scSemiProfiler"]):
            s = F_dec[(F_dec.cohort == coh) & (F_dec.source == src)].set_index("component")
            for m, (k_c, al, ht) in enumerate(zip(comps, alphas, hatches)):
                x = i_s * 1.18 + m * 0.31
                v = s.loc[k_c, "slope_corrected"]
                ax.bar(x, v, 0.27, color=SRC_COL[src], alpha=al, edgecolor="k", lw=0.6, hatch=ht, zorder=2)
                ax.text(x, v + 0.04, f"{s.loc[k_c, 'r']:.2f}", ha="center", va="bottom", fontsize=5.2,
                        fontweight="bold" if src == "PBISC" else "normal",
                        bbox=dict(fc="white", ec="none", alpha=0.92, pad=0.25), zorder=5)
            pos.append(i_s * 1.18 + 0.31)
        ax.axhline(1.0, color="k", lw=0.7, ls=(0, (3, 2)), zorder=1)
        ax.set_xticks(pos); ax.set_xticklabels(["PBISC", "scSemi"], fontsize=5.6)
        ax.set_ylim(0, 2.02); ax.set_ylabel("Corr. slope vs Real", fontsize=5.6)
        ax.legend(handles=[Patch(fc=PALETTE["gray_2"], alpha=al, ec="k", lw=0.5, hatch=ht, label=lbl)
                           for lbl, al, ht in zip(["Total", "Δ share", "Δ expr."], alphas, hatches)],
                  loc="upper right", fontsize=5.1, handlelength=0.9, ncol=3)
        panel_title(ax, f"{NAME[coh]} · shift decomposition (r above)", size=6.3); panel_label(ax, L[8 + j], x=-0.32, y=1.06)

    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_Coverage_Depth_and_Decomposition")


def build_suppfig_rare_abund(st_df, snum=17):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.85), facecolor="white")
    outer = fig.add_gridspec(1, 2, wspace=0.38, left=0.11, right=0.97, top=0.93, bottom=0.06)
    gs_l = outer[0].subgridspec(4, 1, hspace=0.52)
    gs_r = outer[1].subgridspec(4, 1, hspace=0.52)

    sep_4 = ["B_c07-MZB1", "Macrophage", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67"]
    cov_4 = ["B_ASC_IgG", "CD14mono_dex", "pDC", "HSPC"]
    for i, s in enumerate(sep_4):
        ax = fig.add_subplot(gs_l[i, 0])
        draw_rare_patient_1x(ax, "sepsis", s, st_df)
        panel_label(ax, L[i], x=-0.22, y=1.05)
    for i, s in enumerate(cov_4):
        ax = fig.add_subplot(gs_r[i, 0])
        draw_rare_patient_1x(ax, "covid", s, st_df)
        panel_label(ax, L[4 + i], x=-0.22, y=1.05)

    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=SRC_COL[s], label=s, ms=4.5) for s in SRC],
               loc="upper center", bbox_to_anchor=(0.52, 0.992), ncol=3, fontsize=6.2)
    fig.text(0.5, 0.012,
             "Points: patient × condition units; bars: mean ± s.d. Brackets: centered patient-block bootstrap vs Real (49,999 resamples, BH q; * q<0.05, ** q<0.01, *** q<0.001, ns otherwise).",
             ha="center", fontsize=5.3, color=PALETTE["gray_2"])
    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_Rare_Populations_Across_Patients")


def build_suppfig_rare_heatmaps(scache, snum=18):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.15), facecolor="white")
    specs_l = [
        ("sepsis", "B_c07-MZB1", "Sepsis · Plasmablast (B_c07-MZB1)", "a"),
        ("sepsis", "Macrophage", "Sepsis · Macrophage", "b"),
        ("sepsis", "T_CD4_c07-MKI67", "Sepsis · Cycling CD4 T", "c"),
        ("sepsis", "T_CD8_c06-MKI67", "Sepsis · Cycling CD8 T", "d"),
    ]
    specs_r = [
        ("covid", "B_ASC_IgG", "COVID-19 · IgG ASC (B_ASC_IgG)", "e"),
        ("covid", "CD14mono_dex", "COVID-19 · Dex Monocyte", "f"),
        ("covid", "pDC", "COVID-19 · Plasmacytoid DC (pDC)", "g"),
        ("covid", "HSPC", "COVID-19 · HSPC", "h"),
    ]
    for row_i, (ck, sub, ttl, let) in enumerate(specs_l):
        y0 = 0.755 - row_i * 0.242
        draw_summed_marker_heatmap(fig, [0.028, y0, 0.458, 0.222], scache[ck], sub, ttl, letter=let)
    for row_i, (ck, sub, ttl, let) in enumerate(specs_r):
        y0 = 0.755 - row_i * 0.242
        draw_summed_marker_heatmap(fig, [0.512, y0, 0.462, 0.222], scache[ck], sub, ttl, letter=let)
    return finalize_figure(fig, FIG_OUT / f"SuppFig{snum}_Rare_Population_Summed_Marker_Heatmaps")


def build_suppfig4_5(coh, snum):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.25), facecolor="white")
    bulk_ref = "GSE279448" if coh == "sepsis" else "Zenodo 6497251"

    ax_a = fig.add_axes([0.155, 0.838, 0.315, 0.118])
    d_b = pd.read_csv(DATA_OUT / "fig5_whole_pbmc_vs_bulk.csv"); d_b = d_b[d_b.cohort == coh].reset_index(drop=True)
    yb = np.arange(len(d_b))
    ax_a.barh(yb, d_b.slope_corrected, color=[SRC_COL[s] for s in d_b.source], edgecolor="k", lw=0.6, height=0.58)
    for yi, r in zip(yb, d_b.itertuples()):
        ax_a.text(r.slope_corrected + 0.02, yi, f"{r.slope_corrected:.2f} (r={r.r:.2f}; {r.n_case}/{r.n_ctrl})", va="center", fontsize=5.2)
    ax_a.axvline(d_b.slope_corrected.iloc[0], color=SRC_COL["Real"], lw=0.75, ls=(0, (3, 2)))
    ax_a.set_yticks(yb); ax_a.set_yticklabels(["Real", "PBISC", "scSemi"], fontsize=5.7); ax_a.invert_yaxis()
    ax_a.set_xlim(0, 1.58); ax_a.set_xlabel(f"Corr. slope vs bulk ({bulk_ref})", fontsize=5.7)
    panel_title(ax_a, f"{NAME[coh]} · pseudobulk vs bulk", size=6.4, pad=3)
    fig_panel_label(fig, 0.015, 0.985, "a")

    ax_b = fig.add_axes([0.585, 0.838, 0.385, 0.118])
    F = pd.read_csv(DATA_OUT / "fig5_welch_deg_fidelity.csv"); xx = np.arange(4); w = 0.35
    for j, src in enumerate(["PBISC", "scSemiProfiler"]):
        s = F[(F.cohort == coh) & (F.source == src)].set_index("celltype").reindex(CT4)
        ax_b.bar(xx + (j - 0.5) * w, s.slope_corrected, w, color=SRC_COL[src], edgecolor="k", lw=0.6,
                 label=f"{'scSemi' if src=='scSemiProfiler' else src} slope")
        ax_b.scatter(xx + (j - 0.5) * w, s.r_disattenuated, marker="D", s=14, color="white", edgecolor="k", lw=0.65,
                     zorder=4, label="Disatten. r" if j == 0 else None)
    P_ = F[(F.cohort == coh) & (F.source == "PBISC")].set_index("celltype").AUROC
    S_ = F[(F.cohort == coh) & (F.source == "scSemiProfiler")].set_index("celltype").AUROC
    ax_b.axhline(1.0, color="k", lw=0.65, ls=(0, (3, 2)))
    ax_b.set_xticks(xx); ax_b.set_xticklabels([f"{c.replace(' cell', '').replace('Monocyte', 'Mono')}\n{P_[c]:.2f}·{S_[c]:.2f}" for c in CT4], fontsize=5.3)
    ax_b.set_ylabel("vs Real (corr.)", fontsize=5.7); ax_b.set_ylim(0, 2.05)
    ax_b.legend(loc="upper right", fontsize=5.2, handlelength=0.9, labelspacing=0.15)
    panel_title(ax_b, f"{NAME[coh]} · Welch DEG effects (AUROC below)", size=6.3, pad=3)
    fig_panel_label(fig, 0.515, 0.985, "b")

    gs_v = fig.add_gridspec(4, 3, hspace=0.45, wspace=0.28, left=0.085, right=0.975, top=0.755, bottom=0.048)
    for i, ct in enumerate(CT4):
        for j, src in enumerate(SRC):
            ax_v = fig.add_subplot(gs_v[i, j])
            draw_volcano_1x(ax_v, coh, ct, src, show_xl=(i == 3), show_yl=(j == 0))
            if j == 0:
                panel_label(ax_v, L[2 + i], x=-0.26, y=1.06)

    stem = f"SuppFig{snum}_{'Sepsis' if coh=='sepsis' else 'COVID'}_Welch_DEG_Pseudobulk"
    return finalize_figure(fig, FIG_OUT / stem)


def main():
    apply_publication_style(FigureStyle())
    st_df = pd.read_csv(DATA_OUT / "fig6_rare_patient_bootstrap_10states.csv")
    with open(DATA_OUT / "fig6_rare_summed_cache.pkl", "rb") as fp:
        scache = pickle.load(fp)

    # Default supp_offset=13 aligns after Figures 1-4 (Supplementary Figs. 1-13) in strict Nature citation order:
    #   Fig. 5 -> SuppFig 14 (coverage/depth/decomp), SuppFig 15 (Sepsis Welch), SuppFig 16 (COVID Welch)
    #   Fig. 6 -> SuppFig 17 (rare patient abundance), SuppFig 18 (rare summed marker heatmaps)
    supp_offset = 13
    cli_args = []
    for a in sys.argv[1:]:
        if a.startswith("--supp-offset="):
            supp_offset = int(a.split("=", 1)[1])
        else:
            cli_args.append(a)

    s_cov = supp_offset + 1
    s_sep = supp_offset + 2
    s_covd = supp_offset + 3
    s_rab = supp_offset + 4
    s_rhm = supp_offset + 5

    builders = {
        "fig5": lambda: build_figure5(),
        "fig6": lambda: build_figure6(st_df, scache),
        "sfig14": lambda: build_suppfig_coverage(s_cov),
        "sfig15": lambda: build_suppfig4_5("sepsis", s_sep),
        "sfig16": lambda: build_suppfig4_5("covid", s_covd),
        "sfig17": lambda: build_suppfig_rare_abund(st_df, s_rab),
        "sfig18": lambda: build_suppfig_rare_heatmaps(scache, s_rhm),
    }
    wanted = [a for a in cli_args if a in builders] or list(builders)
    audit_path = OUT_ROOT / f"AUDIT_FIG5_6_SFIG{s_cov}_{s_rhm}.json"
    prev = {a["figure"]: a for a in json.loads(audit_path.read_text(encoding="utf-8"))} if audit_path.exists() else {}
    for key in wanted:
        _, aud = builders[key]()
        prev[aud["figure"]] = aud
        print(f"[OK] {aud['figure']}: {aud['width_mm']}x{aud['height_mm']} mm | min_pt={aud['min_font_pt']} | oob={aud['out_of_bounds_count']}", flush=True)

    order = [
        "Figure5_CrossCohort_BulkOnly_Comparison",
        "Figure6_Rare_Cell_Frequency_and_Expression",
        f"SuppFig{s_cov}_Coverage_Depth_and_Decomposition",
        f"SuppFig{s_sep}_Sepsis_Welch_DEG_Pseudobulk",
        f"SuppFig{s_covd}_COVID_Welch_DEG_Pseudobulk",
        f"SuppFig{s_rab}_Rare_Populations_Across_Patients",
        f"SuppFig{s_rhm}_Rare_Population_Summed_Marker_Heatmaps",
    ]
    audit_path.write_text(json.dumps([prev[k] for k in order if k in prev], indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
