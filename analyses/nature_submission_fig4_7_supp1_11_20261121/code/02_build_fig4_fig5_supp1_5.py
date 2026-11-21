"""Build Nature Portfolio 180-mm Main Figure 4 (a-o), Main Figure 5 (a-n), and Supplementary Figures 1-5.
Follows .scientific-figure-making/SKILL.md (figures4papers house style) and runs automated geometry/font audits.
"""
import json
import os
import pathlib
import pickle
import sys
import tempfile
import warnings
import zlib

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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyBboxPatch, Patch, Rectangle

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
F4P_DIR = ROOT / "deliverables" / "manuscript_figures4papers_20261111"
P_DIR = ROOT / "deliverables" / "manuscript_real_pbisc_scsemi_20261110" / "panels"
OUT_ROOT = ROOT / "deliverables" / "nature_submission_fig4_6_20261120"
FIG_OUT = OUT_ROOT / "figures"
DATA_OUT = OUT_ROOT / "data"
FIG_OUT.mkdir(parents=True, exist_ok=True)
DATA_OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(F4P_DIR))
from make_per_population_marker_heatmaps import (
    COVID_IFN_MODULES,
    SEPSIS_CD8_MODULES,
    compute_extended_state_data,
)
from make_rare_similarity_heatmap import (
    CMAP_CORR,
    CMAP_EXPR,
    MARKER_MODULES,
    SFIG2,
    SHORT_LABEL,
    _resolve_1d_positions,
    compute_cohort_data,
)

SRC = ["Real", "PBISC", "scSemiProfiler"]
T6 = ["T cell", "NK cell", "B cell", "Monocyte", "DC", "Other"]
CT4 = ["T cell", "B cell", "NK cell", "Monocyte"]
NAME = {"sepsis": "Sepsis", "covid": "COVID-19"}
L = "abcdefghijklmnopqrstuvwxyz"


def hex_to_rgb(hx: str):
    hx = hx.lstrip("#")
    return [int(hx[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def get_rare_bootstrap_10states() -> pd.DataFrame:
    """Compute 49,999-resample centered patient-block bootstrap for all 5 Sepsis & 5 COVID rare states."""
    out_csv = DATA_OUT / "suppfig2_bootstrap_10states.csv"
    if out_csv.exists():
        return pd.read_csv(out_csv)
    u = pd.read_csv(P_DIR / "rare_patient_units.csv.gz")
    names = {
        "sepsis": ["B_c07-MZB1", "Macrophage", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67", "Mono_c06-CD14-MKI67"],
        "covid": ["B_ASC_IgG", "CD14mono_dex", "CD14mono_IFN", "pDC", "HSPC"],
    }
    nboot = 49999
    tests = []
    for ci, c in enumerate(["sepsis", "covid"]):
        cu = u[u.cohort == c]
        patients = sorted(cu.patient.unique())
        rng = np.random.default_rng(20261110 + ci)
        weights = rng.multinomial(len(patients), np.full(len(patients), 1.0 / len(patients)), size=nboot).astype(float)
        conds = ["control", "ICU control", "case"] if c == "sepsis" else ["control", "case"]
        for sub in names[c]:
            for g in conds:
                values = cu[(cu.subset == sub) & (cu.group == g)].pivot(index="patient", columns="source", values="pct_all").reindex(patients)
                rv = values.Real.to_numpy(float)
                real = rv[np.isfinite(rv)]
                rden = weights @ np.isfinite(rv).astype(float)
                rmean = (weights @ np.nan_to_num(rv)) / np.maximum(rden, 1.0)
                for src in ["PBISC", "scSemiProfiler"]:
                    gv = values[src].to_numpy(float)
                    gen = gv[np.isfinite(gv)]
                    gden = weights @ np.isfinite(gv).astype(float)
                    gmean = (weights @ np.nan_to_num(gv)) / np.maximum(gden, 1.0)
                    delta = float(gen.mean() - real.mean())
                    valid = (rden > 0) & (gden > 0)
                    boot = (gmean - rmean)[valid]
                    err = boot - delta
                    pv = min(1.0, 2.0 * min((1 + np.sum(err >= delta)) / (len(err) + 1), (1 + np.sum(err <= delta)) / (len(err) + 1)))
                    qlo, qhi = np.quantile(err, [0.025, 0.975])
                    tests.append({
                        "cohort": c, "subset": sub, "condition": g, "source": src,
                        "n_real": len(real), "n_generated": len(gen),
                        "mean_real": float(real.mean()), "mean_generated": float(gen.mean()),
                        "diff_pp": delta, "ci_lo": float(delta - qhi), "ci_hi": float(delta - qlo), "p": pv,
                    })
    df_t = pd.DataFrame(tests)
    df_t["q"] = np.nan
    for cc, ix in df_t.groupby("cohort").groups.items():
        ix = np.asarray(ix)
        pv_arr = df_t.loc[ix, "p"].to_numpy()
        order = np.argsort(pv_arr)
        q = np.minimum.accumulate((pv_arr[order] * len(pv_arr) / np.arange(1, len(pv_arr) + 1))[::-1])[::-1]
        back = np.empty(len(pv_arr))
        back[order] = np.minimum(q, 1.0)
        df_t.loc[ix, "q"] = back
    df_t["stars"] = df_t["q"].map(lambda q: "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else "ns")))
    df_t.to_csv(out_csv, index=False)
    return df_t


def get_raresim_cache() -> dict:
    cache_p = DATA_OUT / "raresim_all10_cache.pkl"
    if cache_p.exists():
        with open(cache_p, "rb") as fp:
            return pickle.load(fp)
    print("Computing rare-population HVG/sample similarity & marker matrices (one-time cache)...", flush=True)
    d = {c: compute_cohort_data(c) for c in ("sepsis", "covid")}
    d["covid_ifn"] = compute_extended_state_data(
        "covid", ["B_ASC_IgG", "CD14mono_dex", "CD14mono_IFN", "pDC", "HSPC"], COVID_IFN_MODULES
    )
    d["sepsis_cd8"] = compute_extended_state_data(
        "sepsis", ["B_c07-MZB1", "Macrophage", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67", "Mono_c06-CD14-MKI67"], SEPSIS_CD8_MODULES
    )
    with open(cache_p, "wb") as fp:
        pickle.dump(d, fp)
    return d


# =============================================================================
# PANEL HELPERS AT 1x NATURE SCALE (7.087 in = 180 mm)
# =============================================================================
def draw_fig4_schematic(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 52); ax.axis("off")
    def box(x, y, w, h, text, fc, ec=PALETTE["gray_2"], bold=False, fs=6.1):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.1", fc=fc, ec=ec, lw=0.85))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", color=PALETTE["gray_3"], linespacing=1.16)
    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=PALETTE["gray_2"], lw=0.95, mutation_scale=8))

    box(0.8, 37.5, 18.5, 11.8, "Paired bulk +\nscRNA-seq", PALETTE["green_1"])
    box(25.2, 37.5, 22.8, 11.8, "5-fold LoRA on\npretrained PBISC", "#DCE6F2")
    box(0.8, 21.5, 18.5, 11.8, "Unpaired bulk-only\npatient RNA-seq", PALETTE["green_1"], bold=True)
    box(25.2, 21.5, 22.8, 11.8, "Cohort-adapted\nPBISC", "#DCE6F2", ec=PALETTE["blue_main"], bold=True)
    box(54.0, 21.5, 19.2, 11.8, "Generated scRNA-seq\n(4,000 cells / sample)", "#DCE6F2", ec=PALETTE["blue_main"], bold=True)
    box(54.0, 37.5, 19.2, 11.8, "scSemiProfiler\n(semi-profiling baseline)", PALETTE["red_1"], ec=PALETTE["red_strong"])
    box(79.0, 21.5, 20.2, 27.8, "Evaluation vs Real\n• Gene-mean profile\n• Harmony iLISI & UMAP\n• 6-lineage composition\n• Δ share vs Δ expr.", "#F4F5F7", fs=5.8)

    arrow(19.8, 43.4, 24.7, 43.4)
    arrow(36.6, 37.0, 36.6, 33.8)
    arrow(19.8, 27.4, 24.7, 27.4)
    arrow(48.5, 27.4, 53.5, 27.4)
    arrow(73.7, 27.4, 78.5, 31.0)
    arrow(73.7, 43.4, 78.5, 40.5)

    coh_info = [
        ("Ye et al. (2026) bacterial sepsis",
         "207 bulk samples (3 excl.); 82 paired donors · 122 unpaired bulk-only (68 sepsis, 54 ctrl)\n"
         "488,000 PBISC cells vs 633,072 Real (108 samples) · 5-fold LoRA L1 = 0.36 (null 0.53)"),
        ("Imperial ESKD COVID-19 (Gisby 2022; Stephenson 2025)",
         "336 bulk samples (15 excl.); 57 paired donors (157 visits) · 164 unpaired bulk-only\n"
         "656,000 PBISC cells vs 578,835 Real (187 samples) · 5-fold LoRA L1 = 0.36 (null 0.52)"),
    ]
    for idx, (hdr, body) in enumerate(coh_info):
        x0 = 0.8 + idx * 50.2
        ax.text(x0, 17.2, hdr, fontsize=6.2, fontweight="bold", color=PALETTE["gray_3"], va="top")
        ax.text(x0, 11.8, body, fontsize=5.5, color=PALETTE["gray_2"], va="top", linespacing=1.28)


def draw_gene_counts(ax, coh):
    d = pd.read_csv(P_DIR / f"{coh}_gene_counts_displayed.csv").sort_values("rank")
    m = pd.read_csv(P_DIR / "main_summary_numbers.csv")
    R = {s: m[(m.cohort == coh) & (m.source == s) & (m.metric == "log10_gene_mean_Pearson_r")].value.iloc[0] for s in ["PBISC", "scSemiProfiler"]}
    for src in ["scSemiProfiler", "PBISC", "Real"]:
        y = d[src].values.astype(float)
        ok = y > 0
        if src != "Real":
            ax.scatter(d["rank"][ok], y[ok], s=1.4, color=SRC_COL[src], alpha=0.35, lw=0, rasterized=True)
        med = pd.Series(np.where(ok, y, np.nan)).rolling(151, center=True, min_periods=30).median()
        lab = "Real" if src == "Real" else f"{src} (r = {R[src]:.2f})"
        ax.plot(d["rank"], med, color=SRC_COL[src], lw=1.25 if src == "Real" else 1.05, label=lab, zorder=5 if src == "Real" else 4)
    ax.set_yscale("log"); ax.set_ylim(1e-4, 300); ax.set_xlim(0, 6000); ax.set_xticks([0, 2000, 4000, 6000])
    ax.set_xlabel("Gene rank in Real"); ax.set_ylabel("Mean UMI per cell")
    h, l = ax.get_legend_handles_labels()
    order = [l.index("Real")] + [i for i, x in enumerate(l) if x != "Real"]
    ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper right", fontsize=5.6,
              handlelength=1.2, labelspacing=0.22, title="log10 Pearson r vs Real", title_fontsize=5.5)
    panel_title(ax, f"{NAME[coh]} · gene-level counts")


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
                arrowprops=dict(arrowstyle="-|>", lw=0.75, color="k", shrinkA=0, shrinkB=0, mutation_scale=6))
    ax.annotate("", xy=(0.02, 0.20), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", lw=0.75, color="k", shrinkA=0, shrinkB=0, mutation_scale=6))
    ax.text(0.11, -0.02, "UMAP1", transform=ax.transAxes, fontsize=5.4, ha="center", va="top")
    ax.text(-0.02, 0.11, "UMAP2", transform=ax.transAxes, fontsize=5.4, rotation=90, ha="right", va="center")


def draw_umap_source(ax, coh):
    full = pd.read_csv(F4P_DIR / f"umap_{coh}_post.csv.gz", index_col=0)
    idx = [g.sample(min(30000, len(g)), random_state=1).index for _, g in full.groupby("source")]
    o = full.loc[[i for ix in idx for i in ix]].sample(frac=1, random_state=2)
    col = np.where(o.source == "real_A", SRC_COL["Real"], SRC_COL["PBISC"])
    ax.scatter(o.umap1, o.umap2, c=col, s=0.55, alpha=0.45, lw=0, rasterized=True)
    _clean_umap(ax, full)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=SRC_COL["Real"], label="Real A", ms=3.8),
        Line2D([], [], marker="o", ls="", color=SRC_COL["PBISC"], label="PBISC", ms=3.8),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize=5.7, handletextpad=0.18, columnspacing=0.7)
    panel_title(ax, f"{NAME[coh]} · Real A vs PBISC", size=6.6, pad=3)


def draw_umap_six(ax, coh):
    full = pd.read_csv(F4P_DIR / f"umap_{coh}_post.csv.gz", index_col=0)
    idx = [g.sample(min(30000, len(g)), random_state=1).index for _, g in full.groupby("source")]
    o = full.loc[[i for ix in idx for i in ix]].sample(frac=1, random_state=4)
    ax.scatter(o.umap1, o.umap2, c=[CT_COL.get(t, PALETTE["neutral"]) for t in o.six_type.astype(str)],
               s=0.55, alpha=0.60, lw=0, rasterized=True)
    _clean_umap(ax, full)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=CT_COL[t], label=t, ms=3.6) for t in T6],
              loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3, fontsize=5.4,
              handletextpad=0.12, columnspacing=0.45, labelspacing=0.20)
    panel_title(ax, f"{NAME[coh]} · six lineages", size=6.6, pad=3)


ILISI_DF = pd.DataFrame([
    ("sepsis", "Real A/B", 1.660, 1.682), ("sepsis", "PBISC", 1.476, 1.569), ("sepsis", "scSemi", 1.000, 1.000),
    ("covid", "Real A/B", 1.598, 1.648), ("covid", "PBISC", 1.429, 1.600), ("covid", "scSemi", 1.000, 1.001),
], columns=["cohort", "source", "pre", "post"])


def draw_ilisi(ax, coh):
    d = ILISI_DF[ILISI_DF.cohort == coh].reset_index(drop=True)
    x = np.arange(3); w = 0.36
    ax.bar(x - w / 2, d.pre, w, color=PALETTE["neutral"], edgecolor="k", lw=0.65, hatch="//", label="PCA (pre)")
    ax.bar(x + w / 2, d.post, w, color=[SRC_COL["Real"], SRC_COL["PBISC"], SRC_COL["scSemiProfiler"]],
           edgecolor="k", lw=0.65, label="Harmony (post)")
    for xi, v in zip(x + w / 2, d.post):
        ax.text(xi, v + 0.025, f"{v:.2f}", ha="center", va="bottom", fontsize=5.6, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(d.source, fontsize=5.8)
    ax.set_ylim(1.0, 2.12); ax.set_yticks([1.0, 1.5, 2.0]); ax.set_ylabel("Median iLISI")
    ax.legend(loc="upper right", fontsize=5.4, handlelength=1.0, labelspacing=0.18)
    ax.text(0.98, 0.32, "1 = separated\n2 = mixed", transform=ax.transAxes, ha="right", va="center",
            fontsize=5.3, color=PALETTE["gray_1"])
    panel_title(ax, f"{NAME[coh]} · manifold mixing", size=6.6, pad=3)


def draw_composition6(ax, coh, show_legend=False):
    a = pd.read_csv(P_DIR / f"{coh}_composition_counts.csv")
    b = pd.read_csv(P_DIR / "scsemi_composition.csv"); b = b[b.cohort == coh]
    tab = pd.concat([a[["source", "cell_type", "percent"]], b[["source", "cell_type", "percent"]]]).pivot_table(
        index="source", columns="cell_type", values="percent"
    ).reindex(index=SRC, columns=T6).fillna(0)
    y = np.arange(3); left = np.zeros(3)
    for t in T6:
        v = tab[t].values
        ax.barh(y, v, left=left, color=CT_COL[t], edgecolor="k", lw=0.55, height=0.60, label=t)
        for yi, (l, w) in enumerate(zip(left, v)):
            if w >= 8.0:
                ax.text(l + w / 2, yi, f"{w:.0f}", ha="center", va="center", fontsize=5.6,
                        fontweight="bold", color="white" if t in ["T cell", "B cell", "Monocyte"] else "k")
        left += v
    ax.set_yticks(y); ax.set_yticklabels(["Real", "PBISC", "scSemi"], fontsize=5.9)
    ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel("Composition (%)")
    if show_legend:
        ax.legend(loc="upper center", bbox_to_anchor=(1.15, -0.22), ncol=6, fontsize=5.4,
                  handlelength=0.85, columnspacing=0.65, handletextpad=0.25)
    panel_title(ax, f"{NAME[coh]} · six-lineage composition", size=6.6, pad=3)


def draw_att_decomp(ax, coh):
    F = pd.read_csv(F4P_DIR / "F2d_decomposition.csv")
    comps = ["total", "composition", "within"]
    alphas = [1.0, 0.62, 0.32]
    hatches = ["", "//", ".."]
    pos = []
    for i, src in enumerate(["PBISC", "scSemiProfiler"]):
        s = F[(F.cohort == coh) & (F.source == src)].set_index("component")
        for m, (k, al, ht) in enumerate(zip(comps, alphas, hatches)):
            x = i * 1.18 + m * 0.31
            v = s.loc[k, "slope_corrected"]
            ax.bar(x, v, 0.27, color=SRC_COL[src], alpha=al, edgecolor="k", lw=0.65, hatch=ht)
            ax.text(x, v + 0.04, f"{s.loc[k, 'r']:.2f}", ha="center", va="bottom", fontsize=5.3,
                    fontweight="bold" if src == "PBISC" else "normal", color=PALETTE["gray_3"])
        pos.append(i * 1.18 + 0.31)
    ax.axhline(1.0, color="k", lw=0.75, ls=(0, (3, 2)))
    ax.set_xticks(pos); ax.set_xticklabels(["PBISC", "scSemi"], fontsize=5.9)
    ax.set_ylim(0, 2.02); ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax.set_ylabel("Corr. slope vs Real")
    ax.legend(handles=[
        Patch(fc=PALETTE["gray_2"], alpha=al, ec="k", lw=0.55, hatch=ht, label=lbl)
        for lbl, al, ht in zip(["Total", "Δ share", "Δ expr."], alphas, hatches)
    ], loc="upper right", fontsize=5.3, handlelength=0.95, labelspacing=0.16)
    ax.text(0.03, 0.96, "Numbers = r", transform=ax.transAxes, fontsize=5.3, color=PALETTE["gray_1"], va="top")
    panel_title(ax, f"{NAME[coh]} · shift decomp.", size=6.5, pad=3)


def draw_att_comp(ax, coh, show_legend=False):
    F = pd.read_csv(F4P_DIR / "F2b_composition_shift.csv")
    yy = np.arange(len(T6))
    for j, src in enumerate(SRC):
        v = F[(F.cohort == coh) & (F.source == src)].set_index("celltype").shift_pp.reindex(T6)
        ax.barh(yy + (j - 1) * 0.26, v, height=0.25, color=SRC_COL[src], edgecolor="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.75)
    ax.set_yticks(yy); ax.set_yticklabels(["T", "NK", "B", "Mono", "DC", "Other"], fontsize=5.7)
    ax.invert_yaxis(); ax.set_xlim(-18, 18); ax.set_xticks([-10, 0, 10])
    ax.set_xlabel("Case − control (pp)")
    if show_legend:
        ax.legend(handles=[Patch(fc=SRC_COL[s], ec="k", lw=0.5, label=("scSemi" if s == "scSemiProfiler" else s)) for s in SRC],
                  loc="upper center", bbox_to_anchor=(1.22, -0.24), ncol=3, fontsize=5.3,
                  handlelength=0.85, columnspacing=0.65)
    panel_title(ax, f"{NAME[coh]} · lineage shift", size=6.5, pad=3)


def build_figure4():
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.32), facecolor="white")

    ax_a = fig.add_axes([0.035, 0.828, 0.945, 0.158])
    draw_fig4_schematic(ax_a); fig_panel_label(fig, 0.012, 0.993, "a")

    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_axes([0.095 + j * 0.485, 0.662, 0.395, 0.128])
        draw_gene_counts(ax, coh); panel_label(ax, L[1 + j], x=-0.18, y=1.06)

    u_specs = [
        ("d", lambda a: draw_umap_source(a, "sepsis")),
        ("e", lambda a: draw_umap_six(a, "sepsis")),
        ("f", lambda a: draw_umap_source(a, "covid")),
        ("g", lambda a: draw_umap_six(a, "covid")),
    ]
    for idx, (lbl, fn) in enumerate(u_specs):
        ax = fig.add_axes([0.052 + idx * 0.238, 0.482, 0.196, 0.132])
        fn(ax); panel_label(ax, lbl, x=-0.17, y=1.06)

    for idx, coh in enumerate(["sepsis", "covid"]):
        yh = 0.348 - idx * 0.114
        ax = fig.add_axes([0.075, yh, 0.215, 0.084])
        draw_ilisi(ax, coh); panel_label(ax, L[7 + idx], x=-0.26, y=1.06)
    for idx, coh in enumerate(["sepsis", "covid"]):
        xj = 0.382 + idx * 0.328
        ax = fig.add_axes([xj, 0.262, 0.228, 0.165])
        draw_composition6(ax, coh, show_legend=(idx == 0)); panel_label(ax, L[9 + idx], x=-0.26, y=1.05)

    r5 = [
        ("l", lambda a: draw_att_decomp(a, "sepsis"), 0.068, 0.175, -0.30),
        ("m", lambda a: draw_att_decomp(a, "covid"), 0.302, 0.175, -0.30),
        ("n", lambda a: draw_att_comp(a, "sepsis", show_legend=True), 0.565, 0.165, -0.32),
        ("o", lambda a: draw_att_comp(a, "covid", show_legend=False), 0.805, 0.165, -0.32),
    ]
    for lbl, fn, xl, wl, xlb in r5:
        ax = fig.add_axes([xl, 0.052, wl, 0.122])
        fn(ax); panel_label(ax, lbl, x=xlb, y=1.06)

    return finalize_figure(fig, FIG_OUT / "Figure4_BulkOnly_Molecular_Fidelity")


# =============================================================================
# MAIN FIGURE 5 (a-n): RARE POPULATIONS & DISEASE DEGs
# =============================================================================
def draw_target_cards_1x(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    def card(y0, h, coh_title, main_title, prog_line, note, cite, is_covid=False):
        ax.add_patch(FancyBboxPatch((1.2, y0), 97.6, h, boxstyle="round,pad=0.25,rounding_size=2.6",
                                    fc="#FFFFFF", ec="#C8D1DA", lw=0.85))
        ax.text(4.8, y0 + h - 5.0, coh_title, fontsize=6.4, fontweight="bold", color=PALETTE["gray_3"], va="top")
        if not is_covid:
            for (cx, cy, w, hh, ang, fc, ec) in [
                (13.0, y0 + h * 0.44, 10.5, 16.0, -12, "#B7C5CC", "#7E929B"),
                (12.8, y0 + h * 0.44, 5.0, 7.8, -12, "#466370", "none"),
                (16.0, y0 + h * 0.54, 8.6, 12.6, 18, "#B7C5CC", "#7E929B"),
                (15.8, y0 + h * 0.54, 3.8, 5.8, 18, "#466370", "none"),
            ]:
                ax.add_patch(Ellipse((cx, cy), w, hh, angle=ang, fc=fc, ec=ec, lw=0.7))
        else:
            ax.add_patch(FancyBboxPatch((9.2, y0 + h * 0.58), 8.6, 4.6, boxstyle="round,pad=0.1,rounding_size=1.8",
                                        fc="#DDD4B0", ec="#9E9574", lw=0.7))
            ax.add_patch(Ellipse((13.5, y0 + h * 0.34), 11.0, 15.0, angle=10, fc="#B7C5CC", ec="#7E929B", lw=0.7))
            ax.add_patch(Ellipse((13.3, y0 + h * 0.34), 5.2, 7.5, angle=10, fc="#466370", ec="none"))
        tx = 24.5
        ax.text(tx, y0 + h - 6.5, main_title, fontsize=6.1, fontweight="bold", color=PALETTE["gray_3"], va="top")
        ax.text(tx, y0 + h - 17.5, prog_line, fontsize=5.6, color=PALETTE["blue_main"], fontweight="bold", va="top")
        ax.text(tx, y0 + h - 27.0, note, fontsize=5.4, color=PALETTE["gray_2"], va="top", linespacing=1.18)
        ax.text(4.8, y0 + 2.8, cite, fontsize=5.3, color=PALETTE["gray_1"], va="bottom")

    card(52.0, 46.5, "Sepsis", "Proliferating monocytes", "CD14 / MKI67 / STMN1",
         "Raw marker co-detection vs\npublished cluster assignment", "Ye et al., Nat. Immunol. 2026", False)
    card(1.5, 46.5, "COVID-19", "Steroid monocyte state", "CD14mono_dex (Wave 2)",
         "Requires visit alignment with\ncorticosteroid treatment", "Stephenson et al., Cell Genom. 2025", True)


PROLIF_GROUPS = [
    "Ped. HC (n=4)", "Ped. ICU (n=4)", "Ped. abdom. (n=7)", "Ped. pulm. (n=8)",
    "Adult HC (n=7)", "Adult ICU (n=6)", "Adult pulm. (n=18)", "Adult other (n=28)",
]
PROLIF_GENES = ["CD14", "MKI67", "STMN1", "MPO", "ELANE"]
PROLIF_HEAT = np.array([
    [65.6, 0.6, 3.7, 1.7, 1.3], [74.5, 1.0, 4.2, 2.3, 2.0],
    [77.2, 6.5, 10.2, 5.7, 5.3], [80.5, 2.9, 5.5, 4.6, 5.2],
    [73.7, 0.3, 4.2, 3.5, 1.1], [86.6, 0.3, 2.7, 2.8, 3.3],
    [86.9, 1.3, 3.5, 4.3, 4.5], [80.9, 1.1, 5.0, 2.6, 4.2],
])
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


def draw_prolif_mono_1x(fig, bbox):
    x0, y0, w_tot, h_tot = bbox
    fig.text(x0, y0 + h_tot * 0.98, "Observed Sepsis monocyte marker detection (82 paired donors)",
             fontsize=6.7, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")
    ax_hm = fig.add_axes([x0 + 0.24 * w_tot, y0 + 0.19 * h_tot, 0.34 * w_tot, 0.62 * h_tot])
    ax_cb = fig.add_axes([x0 + 0.59 * w_tot, y0 + 0.19 * h_tot, 0.018 * w_tot, 0.62 * h_tot])
    ax_dt = fig.add_axes([x0 + 0.70 * w_tot, y0 + 0.19 * h_tot, 0.29 * w_tot, 0.62 * h_tot])

    cmap = LinearSegmentedColormap.from_list("pbisc_greys", ["#F7F7F7", "#D9D9D9", "#767676", "#252525"])
    im = ax_hm.imshow(PROLIF_HEAT, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    for i in range(len(PROLIF_GROUPS)):
        for j in range(len(PROLIF_GENES)):
            v = PROLIF_HEAT[i, j]
            ax_hm.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.3,
                       color="white" if v >= 50 else PALETTE["gray_3"],
                       fontweight="bold" if v >= 50 else "normal")
    ax_hm.set_yticks(np.arange(len(PROLIF_GROUPS))); ax_hm.set_yticklabels(PROLIF_GROUPS, fontsize=5.5)
    ax_hm.set_xticks(np.arange(len(PROLIF_GENES)))
    ax_hm.set_xticklabels(PROLIF_GENES, fontsize=5.4, fontstyle="italic", rotation=28, ha="right")
    ax_hm.tick_params(axis="both", length=0, pad=1.5)
    for sp_ in ax_hm.spines.values():
        sp_.set_visible(True); sp_.set_linewidth(0.55); sp_.set_edgecolor("#CCCCCC")

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_ticks([0, 50, 100]); cb.ax.tick_params(labelsize=5.3, width=0.6, length=1.8, pad=1)
    cb.outline.set_linewidth(0.55); cb.set_label("Donor %", fontsize=5.4, labelpad=1.5)

    for xg in [0, 5, 10]:
        ax_dt.axvline(xg, color="#EAEAEA", lw=0.5, zorder=0)
    for i in range(len(PROLIF_GROUPS)):
        vd = np.asarray(PROLIF_DOTS[i], float); vs = np.asarray(PROLIF_SQRS[i], float)
        ax_dt.scatter(vd, np.full_like(vd, i - 0.15), s=7, marker="o", facecolor=PALETTE["blue_main"], edgecolor="white", lw=0.25, zorder=3)
        ax_dt.scatter(vs, np.full_like(vs, i + 0.16), s=7, marker="s", facecolor="white", edgecolor=PALETTE["gray_1"], lw=0.55, zorder=2)
    ax_dt.set_ylim(len(PROLIF_GROUPS) - 0.4, -0.6); ax_dt.set_xlim(-0.5, 13.2)
    ax_dt.set_xticks([0, 5, 10]); ax_dt.set_yticks([]); ax_dt.spines["left"].set_visible(False)
    ax_dt.set_xlabel("Co-detected mono. (%)", fontsize=5.6)
    ax_dt.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=PALETTE["blue_main"], mec="white", ms=3.5, label="MKI67+STMN1"),
        Line2D([], [], marker="s", ls="", mfc="white", mec=PALETTE["gray_1"], mew=0.6, ms=3.5, label="+MPO/ELANE"),
    ], loc="lower right", fontsize=5.3, handletextpad=0.2, labelspacing=0.15)


def draw_rare_parent_1x(ax, coh, show_legend=False):
    r = pd.read_csv(P_DIR / "rare_parent_shares.csv"); r = r[r.cohort == coh]
    subs = list(dict.fromkeys(r[r.source == "Real"].subset))
    short_s = {
        "Mono_c06-CD14-MKI67": "Mono_c06-MKI67", "Mono_c01-CD14-HLA-DRA": "Mono_c01-HLA-DRA",
        "NKT_c05-IFNG-CCL4hi": "NKT_IFNG-CCL4hi",
    }
    y = np.arange(len(subs)); h = 0.26
    for j, src in enumerate(SRC):
        v = r[r.source == src].set_index("subset").percent_parent.reindex(subs).values
        vv = np.where(v > 0, v, np.nan)
        ax.barh(y + (j - 1) * h, vv, h, left=0.004, color=SRC_COL[src], edgecolor="k", lw=0.35,
                label=("scSemi" if src == "scSemiProfiler" else src))
        for yi, vi in zip(y + (j - 1) * h, v):
            if not vi > 0:
                ax.text(0.0045, yi, "0", va="center", fontsize=5.3, color=SRC_COL[src], fontweight="bold")
    ax.set_xscale("log"); ax.set_xlim(0.004, 120); ax.set_yticks(y)
    ax.set_yticklabels([short_s.get(s, s) for s in subs], fontsize=5.4)
    ax.invert_yaxis(); ax.set_xlabel("% of parent lineage (cases)", fontsize=5.9)
    if show_legend:
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=5.3,
                  handlelength=0.8, columnspacing=0.45, handletextpad=0.22)
    panel_title(ax, f"{NAME[coh]} · rare shares", size=6.5, pad=12 if show_legend else 3)


def draw_deg_fidelity_1x(ax, coh):
    F = pd.read_csv(F4P_DIR / "F1_deg_fidelity.csv"); d = F[F.cohort == coh]; ref = d.drop_duplicates("celltype")
    ax.axhline(1.0, color=PALETTE["gray_1"], lw=0.75, ls=(0, (3, 2)), zorder=0)
    ax.axvline(0.5, color=PALETTE["neutral"], lw=0.65, zorder=0)
    ax.scatter(ref.real_splithalf_AUROC, np.ones(len(ref)), marker="D", s=22, facecolor="white",
               edgecolor=SRC_COL["Real"], lw=0.9, zorder=3, label="Real split-half")
    off_map = {
        ("sepsis", "PBISC", "T cell"): (0.0, 0.16, "center"), ("sepsis", "PBISC", "B cell"): (-0.02, 0.18, "center"),
        ("sepsis", "PBISC", "NK cell"): (-0.05, 0.06, "right"), ("sepsis", "PBISC", "Monocyte"): (0.02, -0.02, "left"),
        ("sepsis", "scSemiProfiler", "T cell"): (0.0, 0.18, "center"), ("sepsis", "scSemiProfiler", "B cell"): (0.04, 0.05, "left"),
        ("sepsis", "scSemiProfiler", "NK cell"): (-0.03, 0.14, "right"), ("sepsis", "scSemiProfiler", "Monocyte"): (0.03, -0.08, "left"),
        ("covid", "PBISC", "B cell"): (-0.04, -0.02, "right"), ("covid", "PBISC", "T cell"): (-0.02, -0.18, "center"),
        ("covid", "PBISC", "NK cell"): (-0.05, -0.03, "right"), ("covid", "PBISC", "Monocyte"): (0.0, -0.22, "center"),
        ("covid", "scSemiProfiler", "NK cell"): (-0.04, 0.06, "right"), ("covid", "scSemiProfiler", "T cell"): (-0.01, 0.22, "center"),
        ("covid", "scSemiProfiler", "Monocyte"): (0.03, 0.0, "left"), ("covid", "scSemiProfiler", "B cell"): (0.03, 0.02, "left"),
    }
    for src in ["PBISC", "scSemiProfiler"]:
        s = d[d.source == src]
        ax.scatter(s.AUROC, s.slope_corrected, marker=SRC_MK[src], s=36, color=SRC_COL[src],
                   edgecolor="white", lw=0.45, zorder=4, label=("scSemi" if src == "scSemiProfiler" else src))
        for _, r in s.iterrows():
            dx, dy, ha = off_map.get((coh, src, r.celltype), (0.02, 0.05, "left"))
            lbl_c = r.celltype.replace(" cell", "").replace("Monocyte", "Mono")
            ax.annotate(lbl_c, (r.AUROC, r.slope_corrected), xytext=(r.AUROC + dx, r.slope_corrected + dy),
                        color=SRC_COL[src], fontsize=5.5, fontweight="bold" if src == "PBISC" else "normal",
                        ha=ha, va="center", arrowprops=dict(arrowstyle="-", color=PALETTE["gray_1"], lw=0.45, shrinkA=0, shrinkB=2.5))
    ax.set_xlim(0.45, 1.03); ax.set_ylim(0, 1.98); ax.set_xticks([0.5, 0.75, 1.0])
    ax.set_xlabel("Real-DEG ranking (AUROC)", fontsize=5.9); ax.set_ylabel("Corr. slope vs Real", fontsize=5.9)
    ax.legend(loc="upper left", fontsize=5.3, handletextpad=0.15, labelspacing=0.15, borderaxespad=0.1)
    panel_title(ax, f"{NAME[coh]} · DEG fidelity", size=6.6, pad=3)


RARE_UMAP = {"sepsis": ["B_c07-MZB1", "Macrophage", "pDC", "B_c08-MKI67"], "covid": ["pDC", "DC1", "B_ASC_dividing", "MAIT"]}
RARE_UMAP_COL = [PALETTE["violet"], PALETTE["teal"], PALETTE["highlight"], PALETTE["magenta"]]


def draw_umap_rare_1x(ax, coh):
    o = pd.read_csv(F4P_DIR / f"umap_{coh}_post.csv.gz", index_col=0)
    bg = o.sample(60000, random_state=3)
    ax.scatter(bg.umap1, bg.umap2, c="#E6E6E6", s=0.45, lw=0, rasterized=True)
    h = []
    for st, colr in zip(RARE_UMAP[coh], RARE_UMAP_COL):
        r = o[(o.fine.astype(str) == st) & (o.source == "real_A")]
        g = o[(o.fine.astype(str) == st) & (o.source == "ceiling")]
        ax.scatter(r.umap1, r.umap2, s=4.0, color=colr, lw=0, alpha=0.82, rasterized=True)
        ax.scatter(g.umap1, g.umap2, s=12.0, marker="^", facecolor=colr, edgecolor="k", lw=0.35, zorder=5, rasterized=True)
        h.append(Patch(color=colr, label=f"{st} ({len(r)}/{len(g)})"))
    h += [
        Line2D([], [], marker="o", ls="", color=PALETTE["gray_1"], label="Real A", ms=3.2),
        Line2D([], [], marker="^", ls="", mfc=PALETTE["gray_1"], mec="k", mew=0.4, label="PBISC", ms=3.8),
    ]
    _clean_umap(ax, o)
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, fontsize=5.3,
              handlelength=0.8, columnspacing=0.45, labelspacing=0.18)
    panel_title(ax, f"{NAME[coh]} · rare localization", size=6.5, pad=3)


def draw_hvg_ceiling_1x(fig, bbox, data, letter):
    coh = data["cohort"]
    mat_4x5 = data["mat_4x5"]
    pair_stats = data["pair_stats"]
    x0, y0, w_box, h_box = bbox
    fig_panel_label(fig, x0, y0 + h_box * 0.99, letter)
    fig.text(x0 + 0.032, y0 + h_box * 0.985, f"{NAME[coh]} · 2,000-HVG specificity & sample ceiling",
             fontsize=6.5, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")

    ax_m = fig.add_axes([x0 + 0.24 * w_box, y0 + 0.20 * h_box, 0.34 * w_box, 0.62 * h_box])
    ax_p = fig.add_axes([x0 + 0.62 * w_box, y0 + 0.20 * h_box, 0.275 * w_box, 0.62 * h_box])

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
                ax_m.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.3,
                          fontweight="bold" if (i == j or j == 4) else "normal", color=tc)
            else:
                ax_m.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor="#EFEFEF", edgecolor="#CCCCCC", hatch="///", lw=0.4))
                ax_m.text(j, i, "0", ha="center", va="center", fontsize=5.3, color=PALETTE["gray_2"], fontweight="bold")
    sh_lbls = [SHORT_LABEL[s].replace("T_CD4_c07-", "CD4_").replace("Mono_c06-", "Mono_") for s in SFIG2[coh]]
    ax_m.set_yticks(np.arange(4)); ax_m.set_yticklabels(sh_lbls, fontsize=5.4, fontweight="bold")
    ax_m.set_xticks(np.arange(5)); ax_m.set_xticklabels(sh_lbls + ["scSemi"], fontsize=5.3, rotation=32, ha="right")
    ax_m.tick_params(axis="both", length=0, pad=1.5)
    for sp_ in ax_m.spines.values():
        sp_.set_visible(True); sp_.set_linewidth(0.6); sp_.set_edgecolor(PALETTE["gray_3"])

    ax_p.axvline(0.0, color=PALETTE["neutral"], lw=0.65, ls=(0, (2, 2)))
    for i, sub in enumerate(SFIG2[coh]):
        st = pair_stats[sub]
        for key, yy, col, mk in [("rr", i - 0.22, SRC_COL["Real"], "D"), ("rp", i + 0.02, SRC_COL["PBISC"], "o"), ("off", i + 0.24, PALETTE["gray_1"], "s")]:
            arr = st[key]; arr = arr[np.isfinite(arr)]
            if len(arr) == 0:
                continue
            q10, q25, med, q75, q90 = np.quantile(arr, [0.10, 0.25, 0.50, 0.75, 0.90])
            ax_p.plot([q10, q90], [yy, yy], color=col, lw=0.65, alpha=0.8)
            ax_p.plot([q25, q75], [yy, yy], color=col, lw=1.45)
            ax_p.scatter([med], [yy], marker=mk, s=11, facecolor="white" if key == "off" else col, edgecolor=col, lw=0.6, zorder=4)
        mrp = float(np.nanmedian(st["rp"])); mrr = float(np.nanmedian(st["rr"]))
        ax_p.text(1.06, i - 0.12, f"{mrp:.2f}", fontsize=5.3, fontweight="bold", color=SRC_COL["PBISC"], va="center")
        ax_p.text(1.06, i + 0.16, f"({mrr:.2f})", fontsize=5.3, color=PALETTE["gray_2"], va="center")
    ax_p.set_ylim(3.55, -0.55); ax_p.set_xlim(-0.45, 1.05); ax_p.set_yticks([])
    ax_p.set_xticks([0.0, 0.5, 1.0]); ax_p.spines["left"].set_visible(False)
    ax_p.set_xlabel("Sample-pair r", fontsize=5.6)
    ax_p.legend(handles=[
        Line2D([], [], marker="D", ls="", color=SRC_COL["Real"], ms=3.2, label="Real–Real"),
        Line2D([], [], marker="o", ls="", color=SRC_COL["PBISC"], ms=3.2, label="Real–PBISC"),
        Line2D([], [], marker="s", ls="", mfc="white", mec=PALETTE["gray_1"], ms=3.2, label="Cross"),
    ], loc="upper right", bbox_to_anchor=(1.22, 1.21), ncol=3, fontsize=5.3, handletextpad=0.10, columnspacing=0.30)


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
                        ax_top.text(0.5 * (xr + xp), yb + 0.6, s.stars.iloc[0], ha="center", va="bottom", fontsize=5.5, color=SRC_COL[src])
        ax_bot.set_xlim(-0.55, len(conds) - 0.45); ax_bot.set_ylim(-0.02, y_bot_max); ax_bot.set_yticks([0.0, 0.4, 0.8])
        ax_bot.set_xticks(range(len(conds))); ax_bot.set_xticklabels([c[1] for c in conds], fontsize=5.8)
        ax_bot.set_ylabel("% cells / patient", y=0.80, fontsize=5.8)
        ax_top.set_ylim(y_top_min, y_top_max); ax_top.set_yticks([3, 7]); ax_top.spines["bottom"].set_visible(False)
        ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        panel_title(ax_top, short_title or f"{NAME[coh]} · {sub}", size=6.5, pad=3)
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
            ax.scatter(x0 + rng.uniform(-0.07, 0.07, len(v)), v, s=6, color=SRC_COL[src], alpha=0.75, lw=0, zorder=3, rasterized=True)
            ax.errorbar(x0, v.mean(), yerr=v.std(ddof=1) if len(v) > 1 else 0, fmt="_", color="k", ms=5, capsize=1.8, lw=0.8, zorder=4)
            if src != "Real":
                s = st[(st.condition == cond) & (st.source == src)]
                if len(s):
                    yb = ymax * (1.06 if src == "PBISC" else 1.25)
                    ax.plot([xr, xr, x0, x0], [yb, yb + hh, yb + hh, yb], color=SRC_COL[src], lw=0.6, zorder=5)
                    ax.text(0.5 * (xr + x0), yb + hh + ymax * 0.01, s.stars.iloc[0], ha="center", va="bottom",
                            fontsize=5.5, color=SRC_COL[src], fontweight="bold" if s.stars.iloc[0] != "ns" else "normal")
    ax.set_xticks(range(len(conds))); ax.set_xticklabels([c[1] for c in conds], fontsize=5.8)
    ax.set_ylim(0, ymax * 1.48); ax.set_ylabel("% cells / patient", fontsize=5.8)
    panel_title(ax, short_title or f"{NAME[coh]} · {sub}", size=6.4, pad=3)


def draw_pop_marker_heatmap_1x(fig, bbox, coh_data, sub, title_str, target_mods, letter=None):
    coh = coh_data["cohort"]
    mdf_all = coh_data["mdf"]
    Z_all = coh_data["Z_mk"]
    marker_meta_all = coh_data["marker_meta"]
    avail_modules = coh_data.get("modules", MARKER_MODULES[coh])

    sub_mask = (mdf_all["subset"] == sub).to_numpy()
    mdf = mdf_all.iloc[np.flatnonzero(sub_mask)].reset_index(drop=True)
    row_indices, marker_meta, selected_modules = [], [], []
    mod_lookup = {m_name: (m_name, m_col, m_genes) for m_name, m_col, m_genes in avail_modules}
    for m_name in target_mods:
        _, m_col, m_genes = mod_lookup[m_name]
        selected_modules.append((m_name, m_col, m_genes))
        for k, (mk_mod, mk_col, mk_gene) in enumerate(marker_meta_all):
            if mk_mod == m_name:
                row_indices.append(k); marker_meta.append((mk_mod, mk_col, mk_gene))
    Z_sub = Z_all[np.ix_(np.flatnonzero(sub_mask), np.array(row_indices, int))]
    N, G = len(mdf), len(marker_meta)
    nsc = int((mdf["source"] == "scSemiProfiler").sum())

    x0, y0, w_box, h_box = bbox
    if letter is not None:
        fig_panel_label(fig, x0, y0 + h_box * 0.99, letter)
    n_tag = f"n={N}; scSemi: 0" if nsc == 0 else f"n={N}"
    fig.text(x0 + (0.032 if letter else 0.01), y0 + h_box * 0.985, f"{title_str} ({n_tag})",
             fontsize=6.2, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")

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
                 fontsize=5.3, fontweight="bold", color=SRC_COL[sn], ha="center", va="bottom")

    x_ml = mx - gx - mod_w
    for k, (_, _, gene) in enumerate(marker_meta):
        yr = my + (1.0 - (k + 0.5) / G) * mh
        fig.text(x_ml - 0.01 * w_box, yr, gene, fontsize=5.3, fontstyle="italic", fontweight="bold",
                 color=PALETTE["gray_3"], ha="right", va="center")
    cur = 0
    for m_name, _, mgenes in selected_modules:
        end_r = cur + len(mgenes) - 1
        ym = my + (1.0 - (0.5 * (cur + end_r) + 0.5) / G) * mh
        fig.text(x_ml - 0.14 * w_box, ym, m_name.replace(" response", "\nresp.").replace(" cycle", "\ncycle"),
                 fontsize=5.3, fontweight="bold", color=PALETTE["gray_2"], ha="right", va="center")
        if end_r < G - 1:
            ax_hm.axhline(end_r + 0.5, color="white", lw=0.85, zorder=3)
            ax_hm.axhline(end_r + 0.5, color=PALETTE["gray_3"], lw=0.5, zorder=4)
        cur = end_r + 1

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_ticks([-1, 0, 1, 2]); cb.ax.tick_params(labelsize=5.3, width=0.5, length=1.6, pad=1)
    cb.outline.set_linewidth(0.5); cb.set_label("Z", fontsize=5.4, labelpad=1.5)


def build_figure5(st_df=None, rcache=None):
    """Main Figure 5 (a-h): Lineage-resolved disease DEG & transcriptional signature fidelity
    across Bacterial Sepsis and ESKD COVID-19 (strictly non-overlapping with rare-cell Figure 6).
    """
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.15), facecolor="white")

    # ---- Row 1: a, b (Whole-PBMC pseudobulk vs cohort bulk RNA-seq) & c, d (AUROC vs slope DEG fidelity) ----
    for idx_c, (coh, bulk_ref, let_b, let_d) in enumerate([
        ("sepsis", "GSE279448", "a", "c"),
        ("covid", "Zenodo 6497251", "b", "d"),
    ]):
        yb_top = 0.855 - idx_c * 0.145
        ax_bk = fig.add_axes([0.115, yb_top, 0.255, 0.100])
        d_b = pd.read_csv(F4P_DIR / "F2a_bulk_reference.csv")
        d_b = d_b[d_b.cohort == coh].reset_index(drop=True)
        yb = np.arange(len(d_b))
        ax_bk.barh(yb, d_b.slope_corrected, color=[SRC_COL[s] for s in d_b.source], edgecolor="k", lw=0.6, height=0.58)
        for yi, r in zip(yb, d_b.itertuples()):
            ax_bk.text(r.slope_corrected + 0.025, yi, f"{r.slope_corrected:.2f} (r={r.r:.2f})", va="center", fontsize=5.2,
                       fontweight="bold" if r.source == "PBISC" else "normal")
        ax_bk.axvline(d_b.slope_corrected.iloc[0], color=SRC_COL["Real"], lw=0.75, ls=(0, (3, 2)))
        ax_bk.set_yticks(yb); ax_bk.set_yticklabels(["Real", "PBISC", "scSemi"], fontsize=5.6); ax_bk.invert_yaxis()
        ax_bk.set_xlim(0, 1.48); ax_bk.set_xlabel(f"Corr. slope vs bulk ({bulk_ref})", fontsize=5.5)
        panel_title(ax_bk, f"{NAME[coh]} · whole-PBMC vs bulk", size=6.2, pad=3)
        fig_panel_label(fig, 0.012, yb_top + 0.125, let_b)

    ax_c = fig.add_axes([0.455, 0.728, 0.228, 0.228]); draw_deg_fidelity_1x(ax_c, "sepsis"); fig_panel_label(fig, 0.392, 0.980, "c")
    ax_d = fig.add_axes([0.748, 0.728, 0.228, 0.228]); draw_deg_fidelity_1x(ax_d, "covid"); fig_panel_label(fig, 0.688, 0.980, "d")

    # ---- Row 2: e, f (Per-lineage reliability-corrected slope, disattenuated r & AUROC) ----
    F = pd.read_csv(F4P_DIR / "F1_deg_fidelity.csv")
    xx = np.arange(4); w = 0.35
    for idx_c, (coh, let_e) in enumerate([("sepsis", "e"), ("covid", "f")]):
        xl = 0.088 + idx_c * 0.485
        ax_e = fig.add_axes([xl, 0.510, 0.395, 0.145])
        for j, src in enumerate(["PBISC", "scSemiProfiler"]):
            s = F[(F.cohort == coh) & (F.source == src)].set_index("celltype").reindex(CT4)
            ax_e.bar(xx + (j - 0.5) * w, s.slope_corrected, w, color=SRC_COL[src], edgecolor="k", lw=0.6,
                     label=f"{'scSemi' if src=='scSemiProfiler' else src} slope")
            ax_e.scatter(xx + (j - 0.5) * w, s.r_disattenuated, marker="D", s=15, color="white", edgecolor="k", lw=0.65,
                         zorder=4, label="Disattenuated r" if j == 0 else None)
        P_ = F[(F.cohort == coh) & (F.source == "PBISC")].set_index("celltype").AUROC
        S_ = F[(F.cohort == coh) & (F.source == "scSemiProfiler")].set_index("celltype").AUROC
        ax_e.axhline(1.0, color="k", lw=0.65, ls=(0, (3, 2)))
        ax_e.set_xticks(xx)
        ax_e.set_xticklabels([f"{c.replace(' cell', '').replace('Monocyte', 'Mono')}\nAUROC {P_[c]:.2f}·{S_[c]:.2f}" for c in CT4], fontsize=5.3)
        ax_e.set_ylabel("Effect size vs Real (corr.)", fontsize=5.6); ax_e.set_ylim(0, 2.08)
        ax_e.legend(loc="upper right", fontsize=5.2, handlelength=0.9, labelspacing=0.15)
        panel_title(ax_e, f"{NAME[coh]} · per-lineage DEG calibration (T, B, NK, Mono)", size=6.3, pad=3)
        fig_panel_label(fig, xl - 0.076, 0.680, let_e)

    # ---- Row 3 (g): Representative Sepsis Monocyte DEG volcanoes (Real / PBISC / scSemi) ----
    fig_panel_label(fig, 0.012, 0.452, "g")
    fig.text(0.042, 0.448, "Sepsis · Monocyte pseudobulk DEGs (IL1R2/ACSL1/TNFAIP6 up; HLA-DMB/CIITA down)",
             fontsize=6.2, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")
    for j, src in enumerate(SRC):
        ax_v = fig.add_axes([0.085 + j * 0.305, 0.278, 0.275, 0.140])
        draw_volcano_1x(ax_v, "sepsis", "Monocyte", src, show_xl=False, show_yl=(j == 0))

    # ---- Row 4 (h): Representative COVID-19 T-cell DEG volcanoes (Real / PBISC / scSemi) ----
    fig_panel_label(fig, 0.012, 0.232, "h")
    fig.text(0.042, 0.228, "COVID-19 · T-cell pseudobulk DEGs (type-I IFN MX1/IFIT1/IFIT3/IFI44L induction)",
             fontsize=6.2, fontweight="bold", color=PALETTE["gray_3"], ha="left", va="top")
    for j, src in enumerate(SRC):
        ax_v = fig.add_axes([0.085 + j * 0.305, 0.052, 0.275, 0.140])
        draw_volcano_1x(ax_v, "covid", "T cell", src, show_xl=True, show_yl=(j == 0))

    return finalize_figure(fig, FIG_OUT / "Figure5_Lineage_Disease_DEG_Fidelity")


# =============================================================================
# SUPPLEMENTARY FIGURES 1-5 (180 mm Nature Extended Data / Supp standard)
# =============================================================================
def build_suppfig1():
    fig = plt.figure(figsize=(NATURE_2COL_IN, 8.40), facecolor="white")
    gs = fig.add_gridspec(4, 2, hspace=0.58, wspace=0.42, left=0.16, right=0.96, top=0.95, bottom=0.06)

    cover_spec = {
        "sepsis": [("Genes in Real file", 37474, SRC_COL["Real"]), ("PBISC panel", 6000, SRC_COL["PBISC"]),
                   ("Cohort bulk genes", 5975, PALETTE["blue_secondary"]), ("Panel in Real", 5967, PALETTE["gray_1"]),
                   ("scSemi genes", 5961, SRC_COL["scSemiProfiler"]), ("Harmony HVGs", 2000, PALETTE["neutral"])],
        "covid": [("Genes in Real file", 33559, SRC_COL["Real"]), ("PBISC panel", 6000, SRC_COL["PBISC"]),
                  ("Cohort bulk genes", 5927, PALETTE["blue_secondary"]), ("Panel in Real", 5791, PALETTE["gray_1"]),
                  ("scSemi genes", 5745, SRC_COL["scSemiProfiler"]), ("Harmony HVGs", 2000, PALETTE["neutral"])],
    }
    s_scope = pd.read_csv(P_DIR / "supplementary_1_source_scope_and_depth.csv")
    k = 0
    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[0, j])
        rows = cover_spec[coh]; y = np.arange(len(rows))
        ax.barh(y, [r[1] for r in rows], color=[r[2] for r in rows], edgecolor="k", lw=0.55, height=0.65)
        for yi, r in zip(y, rows):
            ax.text(r[1] * 1.10, yi, f"{r[1]:,}", va="center", fontsize=5.6)
        ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=5.7)
        ax.invert_yaxis(); ax.set_xscale("log"); ax.set_xlim(1000, 2.2e5); ax.set_xlabel("Genes (log10 scale)")
        panel_title(ax, f"{NAME[coh]} · gene coverage"); panel_label(ax, L[j], x=-0.32, y=1.06)

    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[1, j])
        s = s_scope[s_scope.cohort == coh].set_index("source").reindex(SRC); y = np.arange(3)
        ax.barh(y, s.n_cells / 1e6, color=[SRC_COL[x] for x in SRC], edgecolor="k", lw=0.55, height=0.58)
        for yi, (_, r) in zip(y, s.iterrows()):
            ax.text(r.n_cells / 1e6 + 0.015, yi, f"{int(r.n_cells):,} ({int(r.n_samples)} samples)", va="center", fontsize=5.5)
        ax.set_yticks(y); ax.set_yticklabels(["Real", "PBISC", "scSemi"], fontsize=5.8)
        ax.invert_yaxis(); ax.set_xlim(0, 1.25); ax.set_xlabel("Cells (millions)")
        panel_title(ax, f"{NAME[coh]} · cells & samples"); panel_label(ax, L[2 + j], x=-0.32, y=1.06)

    bins = np.logspace(1.5, 4.7, 90); c_bin = 0.5 * (bins[1:] + bins[:-1])
    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[2, j]); txt = []
        for src in SRC:
            u = np.load(P_DIR / f"{coh}_{src}_all_cell_panel_umi.npy").astype(float)
            h, _ = np.histogram(u[u > 0], bins=bins)
            ax.plot(c_bin, 100 * h / h.sum(), color=SRC_COL[src], lw=1.1, label=("scSemi" if src == "scSemiProfiler" else src))
            txt.append(f"{'scSemi' if src=='scSemiProfiler' else src}: {np.median(u):,.0f}")
        ax.set_xscale("log"); ax.set_xlabel("Panel UMI / cell"); ax.set_ylabel("Cells / log bin (%)")
        ax.legend(loc="upper left", fontsize=5.4, handlelength=1.0)
        ax.text(0.98, 0.96, "Median UMI\n" + "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=5.4)
        panel_title(ax, f"{NAME[coh]} · UMI depth density"); panel_label(ax, L[4 + j], x=-0.32, y=1.06)

    for j, coh in enumerate(["sepsis", "covid"]):
        ax = fig.add_subplot(gs[3, j])
        for src in SRC:
            u = np.sort(np.load(P_DIR / f"{coh}_{src}_all_cell_panel_umi.npy").astype(float))
            step = max(1, len(u) // 4000)
            idx = np.unique(np.r_[np.arange(0, len(u), step), len(u) - 1])
            ax.plot(u[idx], 100 * (idx + 1) / len(u), color=SRC_COL[src], lw=1.1, label=("scSemi" if src == "scSemiProfiler" else src))
        ax.set_xscale("log"); ax.set_xlim(30, 5e4); ax.set_ylim(0, 101)
        ax.set_xlabel("Panel UMI / cell"); ax.set_ylabel("Cumulative cells (%)")
        ax.legend(loc="upper left", fontsize=5.4, handlelength=1.0)
        panel_title(ax, f"{NAME[coh]} · UMI depth ECDF"); panel_label(ax, L[6 + j], x=-0.32, y=1.06)

    return finalize_figure(fig, FIG_OUT / "SuppFig1_Coverage_Cells_Sequencing_Depth")


def build_suppfig2(st_df):
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
             "Points: patients; bars: mean ± s.d. Brackets: two-sided centered patient-block bootstrap vs Real (49,999 resamples, BH q within cohort; * q<0.05, ** q<0.01, *** q<0.001, ns otherwise).",
             ha="center", fontsize=5.4, color=PALETTE["gray_2"])
    return finalize_figure(fig, FIG_OUT / "SuppFig2_Rare_Populations_Across_Patients")


def build_suppfig3(rcache):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.15), facecolor="white")
    specs_l = [
        ("sepsis", "B_c07-MZB1", "Sepsis · Plasmablast (B_c07-MZB1)", ["Plasmablast"], "a"),
        ("sepsis", "Macrophage", "Sepsis · Macrophage", ["Macrophage"], "b"),
        ("sepsis", "T_CD4_c07-MKI67", "Sepsis · Cycling CD4 T", ["CD4 T", "Cell cycle"], "c"),
        ("sepsis_cd8", "T_CD8_c06-MKI67", "Sepsis · Cycling CD8 T", ["CD8 T", "Cell cycle"], "d"),
    ]
    specs_r = [
        ("covid", "B_ASC_IgG", "COVID-19 · IgG ASC (B_ASC_IgG)", ["IgG ASC"], "e"),
        ("covid", "CD14mono_dex", "COVID-19 · Dex Monocyte", ["Dex Mono"], "f"),
        ("covid", "pDC", "COVID-19 · Plasmacytoid DC (pDC)", ["pDC"], "g"),
        ("covid", "HSPC", "COVID-19 · HSPC", ["HSPC"], "h"),
    ]
    for row_i, (ck, sub, ttl, mods, let) in enumerate(specs_l):
        y0 = 0.755 - row_i * 0.242
        draw_pop_marker_heatmap_1x(fig, [0.028, y0, 0.458, 0.222], rcache[ck], sub, ttl, mods, letter=let)
    for row_i, (ck, sub, ttl, mods, let) in enumerate(specs_r):
        y0 = 0.755 - row_i * 0.242
        draw_pop_marker_heatmap_1x(fig, [0.512, y0, 0.462, 0.222], rcache[ck], sub, ttl, mods, letter=let)
    return finalize_figure(fig, FIG_OUT / "SuppFig3_Rare_Population_Marker_Heatmaps")


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


def draw_volcano_1x(ax, coh, ct, src, show_xl=False, show_yl=False):
    d = pd.read_csv(P_DIR / f"deg_{coh}_{ct.replace(' ', '')}_{src}.csv.gz")
    q = -np.log10(d.q.clip(lower=1e-300)); yclip = 20.0
    deg = d.deg.astype(bool); qq = q.clip(upper=yclip)
    ax.scatter(d.logFC[~deg], qq[~deg], s=1.2, color=PALETTE["neutral"], lw=0, rasterized=True)
    ax.scatter(d.logFC[deg], qq[deg], s=2.2, color=SRC_COL[src], lw=0, rasterized=True)
    ax.axhline(-np.log10(0.05), color=PALETTE["gray_1"], lw=0.55, ls=(0, (2, 2)))
    ax.axvline(-1, color=PALETTE["gray_1"], lw=0.55, ls=(0, (2, 2)))
    ax.axvline(1, color=PALETTE["gray_1"], lw=0.55, ls=(0, (2, 2)))
    ax.set_xlim(-6, 6); ax.set_ylim(0, 26.5); ax.set_xticks([-4, 0, 4]); ax.set_yticks([0, 10, 20])
    if show_xl: ax.set_xlabel("log2 FC (case − ctrl)", fontsize=5.8)
    if show_yl: ax.set_ylabel("−log10 BH q", fontsize=5.8)
    ax.text(0.97, 0.96, f"DEG: {int(deg.sum())}\n({int(d.n_case.iloc[0])}/{int(d.n_ctrl.iloc[0])})",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.3,
            bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.6), zorder=7)

    deg_df = d[deg].copy(); deg_df["qq"] = qq[deg]
    up_g = deg_df[(deg_df.logFC > 0) & deg_df.gene.isin(VOLC_UP[(coh, ct)])].sort_values("qq", ascending=False).head(3)
    dn_g = deg_df[(deg_df.logFC < 0) & deg_df.gene.isin(VOLC_DN[(coh, ct)])].sort_values("qq", ascending=False).head(2)
    for rank_i, (_, r) in enumerate(up_g.iterrows()):
        # Keep top rank_i==0 label at x<=1.75 so it never overlaps the top-right DEG count box (x>2.3, y>20)
        yt = 23.0 if rank_i == 0 else (18.2 - (rank_i - 1) * 4.2)
        xt = min(1.75, max(1.15, r.logFC * 0.65)) if rank_i == 0 else min(4.6, max(1.8, r.logFC + 0.5))
        ax.annotate(r.gene, xy=(r.logFC, r.qq), xytext=(xt, yt),
                    fontsize=5.3, fontstyle="italic", fontweight="bold", ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.2),
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="#444444"))
    for rank_i, (_, r) in enumerate(dn_g.iterrows()):
        yt = 21.5 - rank_i * 4.5
        ax.annotate(r.gene, xy=(r.logFC, r.qq), xytext=(max(-5.2, min(-1.8, r.logFC - 0.6)), yt),
                    fontsize=5.3, fontstyle="italic", fontweight="bold", ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.78, pad=0.2),
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="#444444"))
    sh_s = "scSemiProfiler" if src == "scSemiProfiler" else src
    panel_title(ax, f"{ct.replace(' cell', '')} · {sh_s}", size=6.1, pad=2.5)


def build_suppfig4_5(coh, snum):
    fig = plt.figure(figsize=(NATURE_2COL_IN, 9.25), facecolor="white")
    bulk_ref = "GSE279448" if coh == "sepsis" else "Zenodo 6497251"

    # Top row: a (whole-PBMC vs same-sample bulk), b (cell-type corrected slope + disattenuated r + AUROC)
    ax_a = fig.add_axes([0.155, 0.838, 0.315, 0.118])
    d_b = pd.read_csv(F4P_DIR / "F2a_bulk_reference.csv"); d_b = d_b[d_b.cohort == coh].reset_index(drop=True)
    yb = np.arange(len(d_b))
    ax_a.barh(yb, d_b.slope_corrected, color=[SRC_COL[s] for s in d_b.source], edgecolor="k", lw=0.6, height=0.58)
    for yi, r in zip(yb, d_b.itertuples()):
        ax_a.text(r.slope_corrected + 0.02, yi, f"{r.slope_corrected:.2f} (r={r.r:.2f}; {r.n_case}/{r.n_ctrl})", va="center", fontsize=5.3)
    ax_a.axvline(d_b.slope_corrected.iloc[0], color=SRC_COL["Real"], lw=0.75, ls=(0, (3, 2)))
    ax_a.set_yticks(yb); ax_a.set_yticklabels(["Real", "PBISC", "scSemi"], fontsize=5.8); ax_a.invert_yaxis()
    ax_a.set_xlim(0, 1.55); ax_a.set_xlabel(f"Corr. slope vs bulk ({bulk_ref})", fontsize=5.8)
    panel_title(ax_a, f"{NAME[coh]} · pseudobulk vs bulk", size=6.5, pad=3)
    fig_panel_label(fig, 0.015, 0.985, "a")

    ax_b = fig.add_axes([0.585, 0.838, 0.385, 0.118])
    F = pd.read_csv(F4P_DIR / "F1_deg_fidelity.csv"); xx = np.arange(4); w = 0.35
    for j, src in enumerate(["PBISC", "scSemiProfiler"]):
        s = F[(F.cohort == coh) & (F.source == src)].set_index("celltype").reindex(CT4)
        ax_b.bar(xx + (j - 0.5) * w, s.slope_corrected, w, color=SRC_COL[src], edgecolor="k", lw=0.6,
                 label=f"{'scSemi' if src=='scSemiProfiler' else src} slope")
        ax_b.scatter(xx + (j - 0.5) * w, s.r_disattenuated, marker="D", s=14, color="white", edgecolor="k", lw=0.65,
                     zorder=4, label="Disatten. r" if j == 0 else None)
    P_ = F[(F.cohort == coh) & (F.source == "PBISC")].set_index("celltype").AUROC
    S_ = F[(F.cohort == coh) & (F.source == "scSemiProfiler")].set_index("celltype").AUROC
    ax_b.axhline(1.0, color="k", lw=0.65, ls=(0, (3, 2)))
    ax_b.set_xticks(xx); ax_b.set_xticklabels([f"{c.replace(' cell', '').replace('Monocyte', 'Mono')}\n{P_[c]:.2f}·{S_[c]:.2f}" for c in CT4], fontsize=5.4)
    ax_b.set_ylabel("vs Real (corr.)", fontsize=5.8); ax_b.set_ylim(0, 2.05)
    ax_b.legend(loc="upper right", fontsize=5.3, handlelength=0.9, labelspacing=0.15)
    panel_title(ax_b, f"{NAME[coh]} · lineage DEG effects (AUROC below)", size=6.4, pad=3)
    fig_panel_label(fig, 0.515, 0.985, "b")

    # Rows c-f: 4 lineages x 3 sources volcano grid
    gs_v = fig.add_gridspec(4, 3, hspace=0.45, wspace=0.28, left=0.085, right=0.975, top=0.755, bottom=0.048)
    for i, ct in enumerate(CT4):
        for j, src in enumerate(SRC):
            ax_v = fig.add_subplot(gs_v[i, j])
            draw_volcano_1x(ax_v, coh, ct, src, show_xl=(i == 3), show_yl=(j == 0))
            if j == 0:
                panel_label(ax_v, L[2 + i], x=-0.26, y=1.06)

    stem = f"SuppFig{snum}_{'Sepsis' if coh=='sepsis' else 'COVID'}_Disease_DEG_Pseudobulk"
    return finalize_figure(fig, FIG_OUT / stem)


def main():
    apply_publication_style(FigureStyle())
    st_df = get_rare_bootstrap_10states()
    rcache = get_raresim_cache()

    builders = {
        "fig4": lambda: build_figure4(),
        "fig5": lambda: build_figure5(st_df, rcache),
        "sfig1": lambda: build_suppfig1(),
        "sfig2": lambda: build_suppfig2(st_df),
        "sfig3": lambda: build_suppfig3(rcache),
        "sfig4": lambda: build_suppfig4_5("sepsis", 4),
        "sfig5": lambda: build_suppfig4_5("covid", 5),
    }
    wanted = [a for a in sys.argv[1:] if a in builders] or list(builders)
    audit_path = OUT_ROOT / "AUDIT_FIG4_5_SFIG1_5.json"
    prev = {a["figure"]: a for a in json.loads(audit_path.read_text(encoding="utf-8"))} if audit_path.exists() else {}
    # Remove old Figure5 stem if rebuilding fig5
    if "fig5" in wanted:
        prev.pop("Figure5_Rare_Populations_Disease_DEGs", None)
        for ext in ("png", "pdf", "svg"):
            old_f = FIG_OUT / f"Figure5_Rare_Populations_Disease_DEGs.{ext}"
            if old_f.exists():
                old_f.unlink()
    for key in wanted:
        saved, aud = builders[key]()
        prev[aud["figure"]] = aud
        print(f"[OK] {aud['figure']}: {aud['width_mm']}x{aud['height_mm']} mm | min_pt={aud['min_font_pt']} | oob={aud['out_of_bounds_count']}")

    order = [
        "Figure4_BulkOnly_Molecular_Fidelity",
        "Figure5_Lineage_Disease_DEG_Fidelity",
        "SuppFig1_Coverage_Cells_Sequencing_Depth",
        "SuppFig2_Rare_Populations_Across_Patients",
        "SuppFig3_Rare_Population_Marker_Heatmaps",
        "SuppFig4_Sepsis_Disease_DEG_Pseudobulk",
        "SuppFig5_COVID_Disease_DEG_Pseudobulk",
    ]
    audit_path.write_text(json.dumps([prev[k] for k in order if k in prev], indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
