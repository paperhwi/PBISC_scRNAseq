import nbformat as nbf
C=[]
C.append(nbf.v4.new_markdown_cell("""# PBISC generated files vs the original single-cell datasets — one figure per cohort (2026-10-18)
Sepsis: Ye 2026 (GSE279452) recap, adult cells in the paper scheme. COVID: Stephenson 2025 recap, 187-sample published roster.
Both real datasets are restricted to the PBISC 6,000-gene panel so counts are comparable.

Inputs (`_extract/`, made on the PC): `real_<cohort>.npz` from `code/extract_real.py` (20,000 random scoped cells; label counts over all scoped cells)
and `<cohort>_k1000/k4000.npz` from `code/extract_file_component.py` (4,000 random cells per draw × 5 draws).
Outputs: `figs/pbisc_vs_real_<cohort>.png|pdf|svg`, per-panel files `figs/vs_real_<cohort>_*.png|svg`, `results/t4–t6_*`, `data/pbisc_vs_real_marker_cells.h5ad`."""))
C.append(nbf.v4.new_code_cell(r'''import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, anndata as ad, os
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
BASE = os.environ.get("FC_BASE", ".")
EX, FIG, RES, DAT = [os.path.join(BASE, d) for d in ("_extract", "figs", "results", "data")]
for d in (FIG, RES, DAT): os.makedirs(d, exist_ok=True)
COH = ["sepsis", "covid"]
Z = {f"{c}_{k}": np.load(f"{EX}/{c}_{k}.npz", allow_pickle=True) for c in COH for k in ("k1000", "k4000")}
R = {c: np.load(f"{EX}/real_{c}.npz", allow_pickle=True) for c in COH}
GENES = list(R["sepsis"]["genes"])
SRC = ["real", "k1000", "k4000"]
COL = {"real": "#4A4F5C", "k1000": "#2a78d6", "k4000": "#eb6834"}
LAB = {"real": "Real", "k1000": "PBISC K=1,000", "k4000": "PBISC K=4,000"}
CT = ["T cell", "NK cell", "B cell", "Monocyte", "DC", "Other"]
CTCOL = dict(zip(CT, ["#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]))
INK, INK2, MUTED, GRID, SURF = "#1F2430", "#4A5162", "#8A90A0", "#E4E6EB", "#FFFFFF"
NAME = {"sepsis": "Sepsis — Ye 2026 (GSE279452)", "covid": "COVID-19 — Imperial ESKD (Stephenson 2025)"}
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#B8BDC7", "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF, "svg.fonttype": "none", "pdf.fonttype": 42})'''))
C.append(nbf.v4.new_markdown_cell("## 1. Real labels → the same six types as the PBISC gate"))
C.append(nbf.v4.new_code_cell(r'''def six_sepsis(l):   # Ye 2026 Fig. 1c subsets (recap)
    if l.startswith(("T_CD4", "T_CD8", "NKT_")) or l == "gdT": return "T cell"
    if l.startswith("NK_"): return "NK cell"
    if l.startswith("B_"): return "B cell"            # incl. B_c07-MZB1 / B_c08-MKI67 plasma cells (AIFI L1 'B cell' also holds plasma cells)
    if l.startswith("Mono_") or l == "Macrophage": return "Monocyte"
    if l in ("cDC", "pDC"): return "DC"
    return "Other"                                     # Megakaryocyte, HSC, ILC
def six_covid(l):    # Stephenson 2025 cell_type_detailed (recap)
    if l.startswith("T ") or l in ("MAIT", "NKT"): return "T cell"
    if l.startswith("NK"): return "NK cell"
    if l.startswith("B_"): return "B cell"
    if l in ("CD14mono", "CD16mono", "Int.mono"): return "Monocyte"
    if l in ("DC1", "DC2", "DC3", "pDC", "ASDC"): return "DC"
    return "Other"                                     # HSPC, ILC, Platelet, RBC
SIX = {"sepsis": six_sepsis, "covid": six_covid}
maptab = pd.concat([pd.DataFrame({"cohort": c, "real_label": R[c]["fine_cats"], "cells": R[c]["fine_counts"],
                                  "six_type": [SIX[c](x) for x in R[c]["fine_cats"]]}) for c in COH])
maptab.to_csv(f"{RES}/t4_real_label_to_six_type.csv", index=False); maptab.groupby(["cohort", "six_type"]).cells.sum().unstack()'''))
C.append(nbf.v4.new_code_cell(r'''# cells table (sampled cells, raw marker counts) for real and PBISC
def pb_cells(key):
    z = Z[key]; out = []
    for j in range(5):
        k = f"draw{j}__"; df = pd.DataFrame(z[k + "marker_counts"], columns=GENES)
        df["total_panel"] = z[k + "total"]; df["ngenes_panel"] = z[k + "ngenes"]; df["label"] = z[k + "aifi_l2_sample"]
        df["six_type"] = z[k + "celltype_l1_pred_sample"]; df["sample"] = z[k + "sample_id_sample"]; out.append(df)
    return pd.concat(out, ignore_index=True)
rows = []
for c in COH:
    r = R[c]; df = pd.DataFrame(r["marker_counts"], columns=GENES)
    df["total_panel"] = r["total_panel"]; df["total_all"] = r["total_all"]; df["ngenes_panel"] = r["ngenes_panel"]
    df["label"] = r["fine_sample"]; df["six_type"] = [SIX[c](x) for x in r["fine_sample"]]; df["sample"] = r["sample_sample"]
    df["cohort"] = c; df["source"] = "real"; rows.append(df)
    for k in ("k1000", "k4000"):
        d = pb_cells(f"{c}_{k}"); d["cohort"] = c; d["source"] = k; rows.append(d)
cells = pd.concat(rows, ignore_index=True)
cells.groupby(["cohort", "source"]).agg(n=("total_panel", "size"), median_umi_panel=("total_panel", "median"),
                                        median_genes_panel=("ngenes_panel", "median"), median_umi_all_real=("total_all", "median"))'''))
C.append(nbf.v4.new_code_cell(r'''# per-gene mean UMI per cell over all 6,000 panel genes
gene_tabs = {}
for c in COH:
    r = R[c]; sym = r["gene_symbol"]; inreal = r["measured"]
    t = pd.DataFrame({"gene": sym, "in_real_file": inreal,
                      "measured_in_cohort": Z[f"{c}_k1000"]["measured"],
                      "real_mean_umi": np.where(inreal, r["gene_count_sum"] / len(r["rows"]), np.nan),
                      "real_detected_frac": np.where(inreal, r["gene_det"] / len(r["rows"]), np.nan)})
    for k in ("k1000", "k4000"):
        z = Z[f"{c}_{k}"]; n = sum(len(z[f"draw{j}__rows"]) for j in range(5))
        t[f"{k}_mean_umi"] = sum(z[f"draw{j}__gene_count_sum"] for j in range(5)) / n
        t[f"{k}_detected_frac"] = sum(z[f"draw{j}__gene_det"] for j in range(5)) / n
    gene_tabs[c] = t; t.to_csv(f"{RES}/t5_gene_counts_all6000_{c}.csv", index=False)
def logr(t, k):
    m = t.in_real_file & (t.real_mean_umi > 0) & (t[f"{k}_mean_umi"] > 0)
    return np.corrcoef(np.log10(t.real_mean_umi[m]), np.log10(t[f"{k}_mean_umi"][m]))[0, 1], int(m.sum())
STATS = {c: {k: logr(gene_tabs[c], k) for k in ("k1000", "k4000")} for c in COH}
STATS'''))
C.append(nbf.v4.new_code_cell(r'''# marker profile (log1p CP10k over panel genes) and composition
prof = []
for c in COH:
    for s in SRC:
        d = cells[(cells.cohort == c) & (cells.source == s)]
        for g in GENES:
            x = d[g].to_numpy(); det = x > 0; v = np.log1p(x[det] / d.total_panel.to_numpy()[det] * 1e4)
            prof.append(dict(cohort=c, source=s, gene=g, detected_frac=det.mean(), n_detected=int(det.sum()),
                             mean_log1p_cp10k=v.mean() if det.any() else np.nan, sd_log1p_cp10k=v.std(ddof=1) if det.sum() > 1 else np.nan))
prof = pd.DataFrame(prof); prof.to_csv(f"{RES}/t6_marker_profile_vs_real.csv", index=False)
comp = []
for c in COH:
    t = maptab[maptab.cohort == c].groupby("six_type").cells.sum(); tot = t.sum()
    for ct in CT: comp.append(dict(cohort=c, source="real", six_type=ct, pct=100 * t.get(ct, 0) / tot))
    for k in ("k1000", "k4000"):
        z = Z[f"{c}_{k}"]; cats = list(z["draw0__celltype_l1_pred_cats"]); v = sum(z[f"draw{j}__celltype_l1_pred_allcounts"] for j in range(5))
        for ct in CT: comp.append(dict(cohort=c, source=k, six_type=ct, pct=100 * v[cats.index(ct)] / v.sum()))
comp = pd.DataFrame(comp); comp.to_csv(f"{RES}/t7_six_type_composition_vs_real.csv", index=False)
comp.pivot_table(index=["cohort", "source"], columns="six_type", values="pct")[CT].round(1)'''))
C.append(nbf.v4.new_markdown_cell("## 2. Panels"))
C.append(nbf.v4.new_code_cell(r'''LINEAGE = [("T", ["CD3E", "IL7R", "CCR7", "CD4", "CD8A"]), ("NK", ["NKG7", "GNLY"]), ("B / plasma", ["MS4A1", "CD79A", "JCHAIN"]),
           ("Monocyte", ["CD14", "LYZ", "S100A8", "FCGR3A"]), ("DC", ["FCER1A", "CLEC4C"]), ("Plt / Ery", ["PPBP", "HBB"]), ("Cycling", ["MKI67"])]
ORDER = [g for _, gs in LINEAGE for g in gs]
def T(ax, s, pad=6): ax.set_title(s, loc="left", fontsize=12, fontweight="bold", color=INK, pad=pad)

def p_genes(ax, c):
    r = R[c]; z = Z[f"{c}_k1000"]
    items = [("Real file (all)", int(r["n_genes_file"])), ("PBISC panel", len(z["gene_symbol"])),
             ("Measured in cohort", int(z["measured"].sum())), ("Found in real file", int(r["measured"].sum())),
             ("HVG (Harmony)", 2000)]
    y = np.arange(len(items))[::-1]
    ax.barh(y, [v for _, v in items], color=["#4A4F5C", "#9AA3B2", "#9AA3B2", "#4A4F5C", "#C9CED8"], height=0.6, edgecolor=SURF, lw=2)
    for yy, (_, v) in zip(y, items): ax.text(v * 1.08, yy, f"{v:,}", va="center", fontsize=9.5, color=INK)
    ax.set_xscale("log"); ax.set_xlim(1000, 150000); ax.set_xticks([1000, 10000, 100000], ["1,000", "10,000", "100,000"])
    ax.set_yticks(y, [l for l, _ in items], fontsize=9); ax.set_xlabel("genes (log scale)"); ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    T(ax, "Genes")

def p_cells(ax, c):
    r = R[c]; y = np.arange(3)[::-1]; vals = []
    txt = [f"{int(r['n_scope']):,} cells · {int(r['n_samples'])} samples" + (" (adult, paper scheme)" if c == "sepsis" else " (published roster)")]
    for k, per in (("k1000", 1500), ("k4000", 4000)):
        z = Z[f"{c}_{k}"]; n = sum(int(z[f"draw{j}__n_cells"]) for j in range(5)); ns = int(z["draw0__n_samples"])
        vals.append(n); txt.append(f"{n:,} = {ns} bulk samples × {per:,} × 5 draws")
    vals = [int(r["n_scope"])] + vals
    ax.barh(y, np.array(vals) / 1e6, color=[COL[s] for s in SRC], height=0.58, edgecolor=SURF, lw=2)
    for yy, v, t in zip(y, vals, txt): ax.text(v / 1e6 + 0.06, yy, t, va="center", fontsize=8.8, color=INK)
    ax.set_yticks(y, [LAB[s] for s in SRC], fontsize=9.5); ax.set_xlim(0, 6.6 if c == "covid" else 5.0); ax.set_xlabel("cells (millions)")
    ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True); T(ax, "Cells")

def p_depth(ax, c):
    bins = np.linspace(2.0, 4.6, 70); mid = (bins[1:] + bins[:-1]) / 2
    d = cells[(cells.cohort == c) & (cells.source == "real")]
    h, _ = np.histogram(np.log10(d.total_all.clip(lower=1)), bins=bins, density=True); ax.plot(mid, h, color="#9AA3B2", lw=1.5, ls="--")
    for s in SRC:
        d = cells[(cells.cohort == c) & (cells.source == s)]
        h, _ = np.histogram(np.log10(d.total_panel.clip(lower=1)), bins=bins, density=True); ax.plot(mid, h, color=COL[s], lw=2)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.75)
    ax.set_xticks([2, 3, 4], ["100", "1,000", "10,000"]); ax.set_xlabel("UMI per cell (log scale)"); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    med = cells[cells.cohort == c].groupby("source").total_panel.median()
    ax.text(0.0, 0.99, "median, panel genes\n" + "\n".join(f"{LAB[s]}: {med[s]:,.0f}" for s in SRC) +
            f"\nReal, all genes: {cells[(cells.cohort == c) & (cells.source == 'real')].total_all.median():,.0f} (dashed)",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.8, color=INK2)
    T(ax, "Depth (6,000 panel genes)")

def p_genecount(ax, c):
    t = gene_tabs[c].copy(); floor = 2e-5
    a = t[t.in_real_file].sort_values("real_mean_umi", ascending=False); b = t[~t.in_real_file].sort_values("k1000_mean_umi", ascending=False)
    t = pd.concat([a, b], ignore_index=True); x = np.arange(len(t)); nr = len(a)
    for k in ("k1000", "k4000"):
        ax.scatter(x, t[f"{k}_mean_umi"].clip(lower=floor), s=1.2, color=COL[k], alpha=0.25, lw=0, rasterized=True)
        sm = pd.Series(np.log10(t[f"{k}_mean_umi"].clip(lower=floor))).rolling(151, center=True, min_periods=20).median()
        ax.plot(x, 10 ** sm, color=COL[k], lw=1.8)
    ax.plot(x[:nr], a.real_mean_umi.clip(lower=floor), color=COL["real"], lw=2.2)
    ax.axvspan(nr - 0.5, len(t) - 0.5, color="#F1F3F6", lw=0)
    ax.text(nr + (len(t) - nr) / 2, 3, f"not in real\nfile ({len(t) - nr})", ha="center", va="top", fontsize=8, color=MUTED)
    ax.set_yscale("log"); ax.set_ylim(floor * 0.8, 30); ax.set_xlim(-20, len(t) + 20)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10], ["0.0001", "0.001", "0.01", "0.1", "1", "10"])
    ax.set_ylabel("mean UMI per cell"); ax.set_xlabel(f"all {len(t):,} panel genes, ordered by real mean count (left = highest)")
    ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    lab = [g for g in ("LYZ", "NKG7", "CD3E", "CD14", "MS4A1", "CD8A", "MKI67") if g in set(t.gene)]
    for n, g in enumerate(sorted(lab, key=lambda g: int(np.where(t.gene.values == g)[0][0]))):
        i = int(np.where(t.gene.values == g)[0][0]); v = t.real_mean_umi.iloc[i]
        ax.annotate(f"{g} (#{i + 1})", (i, v), xytext=(450 + 330 * n, 6.0 if n % 2 == 0 else 2.2), fontsize=8, fontstyle="italic", color=INK,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    s1, s4 = STATS[c]["k1000"], STATS[c]["k4000"]
    ax.text(0.995, 0.60, f"per-gene log mean, real vs PBISC: r = {s1[0]:.2f} (K=1,000), {s4[0]:.2f} (K=4,000); {s1[1]:,} genes",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK2)
    ax.legend(handles=[mpl.lines.Line2D([], [], color=COL[s], lw=2, label=LAB[s]) for s in SRC], frameon=False, fontsize=8.5,
              loc="lower left", ncol=3)
    T(ax, "Gene counts, all 6,000 panel genes (dots = each gene; PBISC line = running median)")

def p_profile(ax, c):
    yy = np.linspace(0, 8.4, 400)
    for i, g in enumerate(ORDER):
        for s, side in (("real", -1), ("k1000", 1), ("k4000", 1)):
            r = prof[(prof.cohort == c) & (prof.source == s) & (prof.gene == g)].iloc[0]
            if r.n_detected < 5: continue
            dens = np.exp(-0.5 * ((yy - r.mean_log1p_cp10k) / r.sd_log1p_cp10k) ** 2); keep = dens > 0.02; w = 0.40 * dens[keep]
            if s == "real": ax.fill_betweenx(yy[keep], i, i - w, color=COL["real"], alpha=0.25, lw=0)
            else: ax.fill_betweenx(yy[keep], i, i + w, color=COL[s], alpha=0.10, lw=0)
            ax.plot(i + side * w, yy[keep], color=COL[s], lw=1.5)
            ax.plot([i, i + side * 0.40], [r.mean_log1p_cp10k] * 2, color=COL[s], lw=1.1)
        dd = [100 * prof[(prof.cohort == c) & (prof.source == s) & (prof.gene == g)].detected_frac.iloc[0] for s in SRC]
        ax.text(i, -0.4, " ".join(f"{v:.0f}" for v in dd), ha="center", va="top", fontsize=7.2, color=MUTED)
    x0 = 0
    for lab, gs in LINEAGE:
        ax.text(x0 + (len(gs) - 1) / 2, 8.55, lab, ha="center", fontsize=9, color=INK2, fontweight="bold")
        if x0 > 0: ax.axvline(x0 - 0.5, color=GRID, lw=1)
        x0 += len(gs)
    ax.axhline(0, color="#B8BDC7", lw=0.8)
    ax.set_xticks(range(len(ORDER)), ORDER, fontsize=9.5, fontstyle="italic", rotation=30, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", pad=14, length=0); ax.set_xlim(-0.6, len(ORDER) - 0.4); ax.set_ylim(-1.0, 8.9); ax.set_yticks([0, 2, 4, 6, 8])
    ax.set_ylabel("log1p(CP10k), detected cells\ncurve = Gaussian(mean, SD)", fontsize=9); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    T(ax, "Marker gene expression — left half Real, right half PBISC (grey numbers = % cells detected: Real, K1k, K4k)", pad=18)

def p_comp(ax, c):
    p = comp[comp.cohort == c].pivot(index="source", columns="six_type", values="pct").loc[SRC, CT]
    y = np.arange(3)[::-1]; left = np.zeros(3)
    for ct in CT:
        v = p[ct].to_numpy(); ax.barh(y, v, left=left, color=CTCOL[ct], height=0.6, edgecolor=SURF, lw=2)
        for yy, l, vv in zip(y, left, v):
            if vv >= 4: ax.text(l + vv / 2, yy, f"{vv:.0f}", ha="center", va="center", fontsize=8.5, color=SURF if ct in ("Monocyte", "DC", "Other") else INK)
        left += v
    ax.set_yticks(y, [LAB[s] for s in SRC], fontsize=9.5); ax.set_xlim(0, 100); ax.set_xlabel("% of cells (real: all scoped cells; PBISC: all 5 draws)")
    ax.legend(handles=[mpl.patches.Patch(color=CTCOL[ct], label=ct) for ct in CT], ncol=6, frameon=False, fontsize=8.2,
              loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=0.9, columnspacing=0.8, handletextpad=0.4)
    T(ax, "Six-type composition", pad=24)

REALTXT = {"sepsis": ("Real — Ye 2026 author scheme (our recap)",
                      ["raw 10x → QC → Harmony (10 PCs) → Leiden clusters", "lineage re-clustering + marker genes → 44 Fig. 1c subsets",
                       "composition r = 0.98 vs paper · subsets → six types (t4)"]),
           "covid": ("Real — Stephenson 2025 author scheme (recap)",
                     ["CellTypist compartments → per-compartment clustering", "canonical markers → 34 types (79% agree with deposited)",
                      "cell_type_detailed → six types (t4)"])}
PBTXT = ("PBISC — frozen scANVI annotator (f01common)",
         ["same weights for every cohort, K and draw; not refit", "aifi_l2: 29 types (+ confidence) → aifi_l1: 9 types",
          "six-type gate: Platelet, Erythrocyte, ILC, Progenitor → Other"])

def p_annot(ax, c):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    for y0, (head, lines), col in ((0.97, REALTXT[c], COL["real"]), (0.47, PBTXT, COL["k1000"])):
        ax.add_patch(mpl.patches.Rectangle((0, y0 - 0.43), 0.012, 0.43, color=col, lw=0))
        ax.text(0.04, y0 - 0.02, head, fontsize=9.3, fontweight="bold", color=INK, va="top")
        for n, l in enumerate(lines):
            ax.text(0.04, y0 - 0.135 - n * 0.095, ("" if n == 0 else "→ ") + l, fontsize=8.4, color=INK2, va="top")
    T(ax, "Cell-type annotation")'''))
C.append(nbf.v4.new_markdown_cell("## 3. One 16:9 figure per cohort"))
C.append(nbf.v4.new_code_cell(r'''def figure(c):
    fig = plt.figure(figsize=(16, 9), dpi=100)
    gs = fig.add_gridspec(3, 3, height_ratios=[0.85, 1.12, 1.25], width_ratios=[1.12, 1.12, 1.0], hspace=0.66, wspace=0.36,
                          left=0.085, right=0.98, top=0.875, bottom=0.075)
    fig.text(0.085, 0.955, f"{NAME[c]}: real data vs PBISC-generated files", fontsize=18, fontweight="bold", color=INK)
    sub = ("Real = adult cells in the Ye 2026 scheme (recap)" if c == "sepsis" else "Real = 187-sample published roster (recap)") + \
          " · PBISC = bulk-only generation, LoRA-5, 5 draws · all counts restricted to the 6,000-gene PBISC panel · expression from 20,000 random cells per source"
    fig.text(0.085, 0.922, sub, fontsize=10, color=INK2)
    p_genes(fig.add_subplot(gs[0, 0]), c); p_cells(fig.add_subplot(gs[0, 1]), c); p_depth(fig.add_subplot(gs[0, 2]), c)
    p_genecount(fig.add_subplot(gs[1, :2]), c); p_comp(fig.add_subplot(gs[1, 2]), c)
    p_profile(fig.add_subplot(gs[2, :2]), c); p_annot(fig.add_subplot(gs[2, 2]), c)
    for ext in ("png", "pdf", "svg"): fig.savefig(f"{FIG}/pbisc_vs_real_{c}.{ext}", dpi=300)
    plt.show()
figure("sepsis")'''))
C.append(nbf.v4.new_code_cell(r'''figure("covid")'''))
C.append(nbf.v4.new_markdown_cell("## 4. Each panel as its own file"))
C.append(nbf.v4.new_code_cell(r'''def save(fn, size, fn_draw, c):
    f, a = plt.subplots(figsize=size); fn_draw(a, c); f.tight_layout()
    for ext in ("png", "svg"): f.savefig(f"{FIG}/vs_real_{c}_{fn}.{ext}", dpi=300)
    plt.close(f)
for c in COH:
    save("genes", (5.4, 3.0), p_genes, c); save("cells", (6.4, 2.6), p_cells, c); save("depth", (5.4, 3.2), p_depth, c)
    save("genecount_all6000", (11, 3.8), p_genecount, c); save("composition", (6.0, 2.8), p_comp, c)
    save("marker_profile", (11, 4.0), p_profile, c); save("annotation", (5.6, 3.6), p_annot, c)
sorted(f for f in os.listdir(FIG) if f.startswith(("vs_real", "pbisc_vs_real")))'''))
C.append(nbf.v4.new_markdown_cell("## 5. Save the sampled cells (real + PBISC) as h5ad"))
C.append(nbf.v4.new_code_cell(r'''obs = cells[["cohort", "source", "sample", "label", "six_type", "total_panel", "ngenes_panel"]].copy()
for k in ["cohort", "source", "sample", "label", "six_type"]: obs[k] = obs[k].astype(str).astype("category")
obs.index = [f"c{i}" for i in range(len(obs))]
a = ad.AnnData(X=cells[GENES].to_numpy(np.float32), obs=obs, var=pd.DataFrame(index=GENES))
a.uns["note"] = ("Raw UMI counts of 19 markers. real: 20,000 random cells per cohort (sepsis Ye 2026 recap adult in-scheme; COVID Stephenson 2025 roster 187). "
                 "PBISC: 4,000 per draw x 5 draws per K. total_panel = UMI over the 6,000 PBISC panel genes.")
a.write_h5ad(f"{DAT}/pbisc_vs_real_marker_cells.h5ad"); a'''))
nb=nbf.v4.new_notebook(); nb.cells=C
nb.metadata["kernelspec"]={"name":"python3","display_name":"Python 3","language":"python"}
nbf.write(nb,'/home/claude/fcfig2/pbisc_vs_real_figure.ipynb'); print('ok')
