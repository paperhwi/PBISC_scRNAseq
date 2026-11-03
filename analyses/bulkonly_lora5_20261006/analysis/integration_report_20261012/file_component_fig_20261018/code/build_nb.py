import nbformat as nbf
C=[]
C.append(nbf.v4.new_markdown_cell("""# PBISC generated files — component figure for PPT (2026-10-18)
Figure of what is inside the bulk-only generated files (Ye 2026 sepsis and Imperial ESKD COVID, K = 1,000 and K = 4,000):
gene numbers, cell numbers, per-gene expression profile (Gaussian fit), depth, and the cell-type annotation route.

Input: `_extract/<dataset>.npz`, made on the PC by `code/extract_file_component.py` (4,000 random cells per draw × 5 draws = 20,000 cells per dataset; label counts over **all** cells).
Outputs: `figs/` (combined figure + each panel as PNG 300 dpi / PDF / SVG), `results/` (tables), `data/file_component_marker_cells.h5ad`."""))
C.append(nbf.v4.new_code_cell(r'''import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, anndata as ad, os, json
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
BASE = os.environ.get("FC_BASE", ".")
EX, FIG, RES, DAT = [os.path.join(BASE, d) for d in ("_extract", "figs", "results", "data")]
for d in (FIG, RES, DAT): os.makedirs(d, exist_ok=True)
DS = ["sepsis_k1000", "sepsis_k4000", "covid_k1000", "covid_k4000"]
NICE = {"sepsis_k1000": "Sepsis K=1,000", "sepsis_k4000": "Sepsis K=4,000", "covid_k1000": "COVID K=1,000", "covid_k4000": "COVID K=4,000"}
Z = {d: np.load(f"{EX}/{d}.npz", allow_pickle=True) for d in DS}
GENES = list(Z["sepsis_k1000"]["genes"])
# colours: K identity (validated pair) and six cell types (validated, fixed order; segments are direct-labelled)
KCOL = {"k1000": "#2a78d6", "k4000": "#eb6834"}
CT = ["T cell", "NK cell", "B cell", "Monocyte", "DC", "Other"]
CTCOL = dict(zip(CT, ["#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]))
INK, INK2, MUTED, GRID, SURF = "#1F2430", "#4A5162", "#8A90A0", "#E4E6EB", "#FFFFFF"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#B8BDC7", "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF, "svg.fonttype": "none", "pdf.fonttype": 42})'''))
C.append(nbf.v4.new_markdown_cell("## 1. Pool the sampled cells and build the tables"))
C.append(nbf.v4.new_code_cell(r'''def pooled(d):
    z = Z[d]; out = []
    for j in range(5):
        k = f"draw{j}__"
        df = pd.DataFrame(z[k + "marker_counts"], columns=GENES)
        df["total_umi"] = z[k + "total"]; df["n_genes"] = z[k + "ngenes"]; df["draw"] = f"draw{j}"
        for col in ["celltype_l1_pred", "aifi_l1", "aifi_l2", "group", "sample_id"]: df[col] = z[k + col + "_sample"]
        out.append(df)
    df = pd.concat(out, ignore_index=True); df["dataset"] = d; return df
cells = pd.concat([pooled(d) for d in DS], ignore_index=True)

# file-level numbers (from the files themselves)
rows = []
for d in DS:
    z = Z[d]; n = [int(z[f"draw{j}__n_cells"]) for j in range(5)]; ns = int(z["draw0__n_samples"])
    rows.append(dict(dataset=NICE[d], samples=ns, cells_per_sample=n[0] // ns, cells_per_draw_file=n[0], draw_files=5,
                     cells_all_draws=sum(n), panel_genes=len(z["gene_symbol"]), measured_genes=int(z["measured"].sum()),
                     unique_symbols=len(set(z["gene_symbol"])), harmony_subset_cells=500 * ns * 5, harmony_hvg=2000,
                     median_umi=float(cells.loc[cells.dataset == d, "total_umi"].median()),
                     median_genes=float(cells.loc[cells.dataset == d, "n_genes"].median())))
files = pd.DataFrame(rows); files.to_csv(f"{RES}/t1_file_numbers.csv", index=False); files'''))
C.append(nbf.v4.new_code_cell(r'''# composition over ALL cells of the 5 draws (not the sample)
comp = []
for d in DS:
    z = Z[d]; cats = list(z["draw0__celltype_l1_pred_cats"]); tot = sum(z[f"draw{j}__celltype_l1_pred_allcounts"] for j in range(5))
    for c, v in zip(cats, tot): comp.append(dict(dataset=NICE[d], celltype_l1_pred=c, cells=int(v), pct=100 * v / tot.sum()))
comp = pd.DataFrame(comp); comp.to_csv(f"{RES}/t2_celltype_l1_composition_all_cells.csv", index=False)
# check: six-type gate = AIFI L1 with Erythrocyte / ILC / Platelet / Progenitor folded into "Other"
z = Z["sepsis_k1000"]; l1 = dict(zip(z["draw0__aifi_l1_cats"], z["draw0__aifi_l1_allcounts"])); g6 = dict(zip(z["draw0__celltype_l1_pred_cats"], z["draw0__celltype_l1_pred_allcounts"]))
print("Other == Ery+ILC+Plt+Prog:", g6["Other"] == l1["Erythrocyte"] + l1["ILC"] + l1["Platelet"] + l1["Progenitor cell"],
      "| shared types equal:", all(g6[k] == l1[k] for k in ["B cell", "DC", "Monocyte", "NK cell", "T cell"]))
comp.pivot(index="dataset", columns="celltype_l1_pred", values="pct").round(1)'''))
C.append(nbf.v4.new_code_cell(r'''# per-gene expression: log1p(CP10k) of detected cells, Gaussian fit (mean, SD) + detection fraction
prof = []
for d in DS:
    s = cells[cells.dataset == d]
    for g in GENES:
        c = s[g].to_numpy(); det = c > 0; x = np.log1p(c[det] / s["total_umi"].to_numpy()[det] * 1e4)
        prof.append(dict(dataset=d, gene=g, detected_frac=det.mean(), mean_log1p_cp10k=x.mean() if det.any() else np.nan,
                         sd_log1p_cp10k=x.std(ddof=1) if det.sum() > 1 else np.nan, n_detected=int(det.sum()), mean_umi_all=c.mean()))
prof = pd.DataFrame(prof); prof.to_csv(f"{RES}/t3_marker_gene_profile.csv", index=False)
prof.pivot(index="gene", columns="dataset", values="detected_frac").loc[GENES].round(3)'''))
C.append(nbf.v4.new_markdown_cell("## 2. Panels"))
C.append(nbf.v4.new_code_cell(r'''LINEAGE = [("T", ["CD3E", "IL7R", "CCR7", "CD4", "CD8A"]), ("NK", ["NKG7", "GNLY"]), ("B / plasma", ["MS4A1", "CD79A", "JCHAIN"]),
           ("Monocyte", ["CD14", "LYZ", "S100A8", "FCGR3A"]), ("DC", ["FCER1A", "CLEC4C"]), ("Plt / Ery", ["PPBP", "HBB"]), ("Cycling", ["MKI67"])]
ORDER = [g for _, gs in LINEAGE for g in gs]

def panel_genes(ax):
    pan = files[["dataset", "panel_genes", "measured_genes"]].iloc[[0, 2]].copy(); pan["dataset"] = ["Sepsis", "COVID"]
    labels = ["Model panel", "Measured\nin cohort", "Harmony\nHVG"]
    y = np.arange(3)[::-1]
    for i, (_, r) in enumerate(pan.iterrows()):
        vals = [r.panel_genes, r.measured_genes, 2000]; off = 0.18 if i == 0 else -0.18
        ax.barh(y + off, vals, height=0.32, color=["#9AA3B2", "#4A5162"][i], edgecolor=SURF, linewidth=2)
        for yy, v in zip(y + off, vals): ax.text(v + 80, yy, f"{v:,}", va="center", fontsize=9.5, color=INK)
    ax.set_yticks(y, labels, fontsize=9.5); ax.set_xlim(0, 7600); ax.set_xticks([0, 2000, 4000, 6000], ["0", "2,000", "4,000", "6,000"])
    ax.set_xlabel("genes"); ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.legend(handles=[mpl.patches.Patch(color="#9AA3B2", label="Sepsis"), mpl.patches.Patch(color="#4A5162", label="COVID")],
              frameon=False, fontsize=9, loc="lower right")
    ax.set_title("Genes", loc="left", fontsize=12, fontweight="bold", color=INK)

def panel_cells(ax):
    f = files.copy(); y = np.arange(4)[::-1]
    ax.barh(y, f.cells_all_draws / 1e6, height=0.55, color=[KCOL["k1000"], KCOL["k4000"]] * 2, edgecolor=SURF, linewidth=2)
    for yy, (_, r) in zip(y, f.iterrows()):
        ax.text(r.cells_all_draws / 1e6 + 0.05, yy, f"{r.cells_all_draws:,}  =  {r.samples} samples × {r.cells_per_sample:,} cells × 5 draws",
                va="center", fontsize=9, color=INK)
    ax.set_yticks(y, f.dataset, fontsize=9.5); ax.set_xlim(0, 6.2); ax.set_xlabel("cells, all 5 draw files (millions)")
    ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.legend(handles=[mpl.patches.Patch(color=KCOL["k1000"], label="K = 1,000"), mpl.patches.Patch(color=KCOL["k4000"], label="K = 4,000")],
              frameon=False, fontsize=9.5, loc="upper right")
    ax.set_title("Cells", loc="left", fontsize=12, fontweight="bold", color=INK)

def panel_profile(ax, cohort):
    xs = np.arange(len(ORDER)); yy = np.linspace(0, 8.4, 400)
    for i, g in enumerate(ORDER):
        for k, side in (("k1000", -1), ("k4000", 1)):
            r = prof[(prof.dataset == f"{cohort}_{k}") & (prof.gene == g)].iloc[0]
            if r.n_detected < 5: continue
            dens = np.exp(-0.5 * ((yy - r.mean_log1p_cp10k) / r.sd_log1p_cp10k) ** 2); keep = dens > 0.02
            w = 0.40 * dens[keep]; y_ = yy[keep]
            ax.fill_betweenx(y_, i, i + side * w, color=KCOL[k], alpha=0.30, lw=0)
            ax.plot(i + side * w, y_, color=KCOL[k], lw=1.6)
            ax.plot([i, i + side * 0.40], [r.mean_log1p_cp10k] * 2, color=KCOL[k], lw=1.2)
        d1 = prof[(prof.dataset == f"{cohort}_k1000") & (prof.gene == g)].detected_frac.iloc[0]
        d4 = prof[(prof.dataset == f"{cohort}_k4000") & (prof.gene == g)].detected_frac.iloc[0]
        ax.text(i, -0.4, f"{100*d1:.0f} | {100*d4:.0f}", ha="center", va="top", fontsize=8, color=MUTED)
    ax.axhline(0, color="#B8BDC7", lw=0.8)
    x0 = 0
    for lab, gs in LINEAGE:
        ax.text(x0 + (len(gs) - 1) / 2, 8.55, lab, ha="center", fontsize=9, color=INK2, fontweight="bold")
        if x0 > 0: ax.axvline(x0 - 0.5, color=GRID, lw=1)
        x0 += len(gs)
    ax.set_xticks(xs, ORDER, fontsize=9.5, fontstyle="italic", rotation=30, ha="right", rotation_mode="anchor"); ax.tick_params(axis="x", pad=14, length=0)
    ax.set_xlim(-0.6, len(ORDER) - 0.4); ax.set_ylim(-1.0, 8.9); ax.set_yticks([0, 2, 4, 6, 8])
    ax.set_ylabel("log1p(CP10k), detected cells\ncurve = Gaussian(mean, SD)", fontsize=9)
    ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)

def panel_depth(ax):
    bins = np.linspace(2.0, 4.2, 70)
    for i, (cohort, ls) in enumerate((("sepsis", "-"), ("covid", "--"))):
        for k in ("k1000", "k4000"):
            v = np.log10(cells.loc[cells.dataset == f"{cohort}_{k}", "total_umi"].clip(lower=1))
            h, e = np.histogram(v, bins=bins, density=True)
            ax.plot((e[1:] + e[:-1]) / 2, h, color=KCOL[k], lw=1.8, ls=ls)
    ax.set_xticks([2, 3, 4], ["100", "1,000", "10,000"]); ax.set_xlabel("UMI per cell (log scale)"); ax.set_ylabel("density")
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    med = files.set_index("dataset").median_umi
    ax.text(0.02, 0.97, "median UMI\n" + "\n".join(f"{k}: {v:,.0f}" for k, v in med.items()), transform=ax.transAxes, va="top", fontsize=8, color=INK2)
    ax.legend(handles=[mpl.lines.Line2D([], [], color=INK2, ls="-", label="Sepsis"), mpl.lines.Line2D([], [], color=INK2, ls="--", label="COVID")],
              frameon=False, fontsize=8, loc="upper right")
    ax.set_title("Depth  (X = one Poisson draw of expected counts)", loc="left", fontsize=12, fontweight="bold", color=INK)

def box(ax, x, y, w, h, title, body, fc="#F1F3F6", ec="#C9CED8", tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.02", fc=fc, ec=ec, lw=1.2))
    ax.text(x + 0.02, y + h - 0.035, title, fontsize=10, fontweight="bold", color=tc, va="top")
    ax.text(x + 0.02, y + h - 0.105, body, fontsize=8.6, color=INK2, va="top", linespacing=1.35)

def arrow(ax, a, b):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14, color="#6A7080", lw=1.4))

def panel_annotation(ax):
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.02); ax.axis("off")
    B = [("Generated cell", "integer UMI counts on 6,000 panel genes", "#F1F3F6", "#C9CED8"),
         ("Frozen scANVI annotator (f01common)", "same weights for every cohort, K and draw; not refit", "#EEF1F8", "#2a78d6"),
         ("aifi_l2  ·  29 types  (+ aifi_l2_confidence)", "AIFI fine labels, e.g. Memory CD4 T, CD14 monocyte, cDC2", "#FFFFFF", "#C9CED8"),
         ("aifi_l1  ·  9 types", "T, NK, B, Monocyte, DC, Platelet, Erythrocyte, ILC, Progenitor", "#FFFFFF", "#C9CED8"),
         ("celltype_l1_pred  ·  six-type gate", "aifi_l1 with Platelet, Erythrocyte, ILC, Progenitor → Other", "#FFFFFF", "#C9CED8")]
    h = 0.15; gap = (0.985 - 5 * h) / 4
    for n, (t, b, fc, ec) in enumerate(B):
        y = 1 - (n + 1) * h - n * gap
        ax.add_patch(FancyBboxPatch((0.01, y), 0.98, h, boxstyle="round,pad=0.004,rounding_size=0.015", fc=fc, ec=ec, lw=1.3))
        ax.text(0.04, y + h * 0.66, t, fontsize=9.6, fontweight="bold", color=INK, va="center")
        ax.text(0.04, y + h * 0.27, b, fontsize=8.6, color=INK2, va="center")
        if n < 4: ax.add_patch(FancyArrowPatch((0.5, y - 0.004), (0.5, y - gap + 0.004), arrowstyle="-|>", mutation_scale=11, color="#6A7080", lw=1.2))
    ax.set_title("Cell-type annotation (stored in obs)", loc="left", fontsize=12, fontweight="bold", color=INK)

def panel_comp(ax):
    p = comp.pivot(index="dataset", columns="celltype_l1_pred", values="pct").loc[[NICE[d] for d in DS], CT]
    y = np.arange(4)[::-1]; left = np.zeros(4)
    for c in CT:
        v = p[c].to_numpy(); ax.barh(y, v, left=left, color=CTCOL[c], height=0.6, edgecolor=SURF, linewidth=2)
        for yy, l, vv in zip(y, left, v):
            if vv >= 4: ax.text(l + vv / 2, yy, f"{vv:.0f}", ha="center", va="center", fontsize=8.5, color=SURF if c in ("Monocyte", "DC", "Other") else INK)
        left += v
    ax.set_yticks(y, p.index, fontsize=9.5); ax.set_xlim(0, 100); ax.set_xlabel("% of all cells (5 draws), celltype_l1_pred")
    ax.legend(handles=[mpl.patches.Patch(color=CTCOL[c], label=c) for c in CT], ncol=6, frameon=False, fontsize=8.2,
              loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=0.9, columnspacing=0.8, handletextpad=0.4)
    ax.set_title("Resulting composition", loc="left", fontsize=12, fontweight="bold", color=INK, pad=24)'''))
C.append(nbf.v4.new_markdown_cell("## 3. Combined 16:9 figure"))
C.append(nbf.v4.new_code_cell(r'''fig = plt.figure(figsize=(16, 9), dpi=100)
gs = fig.add_gridspec(3, 3, height_ratios=[0.95, 1.3, 1.3], width_ratios=[1, 1.15, 1.05], hspace=0.62, wspace=0.40,
                      left=0.07, right=0.975, top=0.875, bottom=0.07)
fig.text(0.075, 0.955, "PBISC generated single-cell files: what is inside", fontsize=18, fontweight="bold", color=INK)
fig.text(0.075, 0.922, "Bulk-only generation, LoRA-5, 5 draws per cohort · Ye 2026 sepsis (122 bulk samples) and Imperial ESKD COVID (164) · "
         "K = 1,000 vs 4,000 · expression and depth from 20,000 random cells per dataset", fontsize=10.5, color=INK2)
panel_genes(fig.add_subplot(gs[0, 0])); panel_cells(fig.add_subplot(gs[0, 1:]))
gsp = gs[1:, :2].subgridspec(2, 1, hspace=0.62)
a1 = fig.add_subplot(gsp[0]); panel_profile(a1, "sepsis"); a1.set_title("Gene expression profile — Sepsis   (Gaussian fit to detected cells · left half K=1,000, right half K=4,000)", loc="left", fontsize=12, fontweight="bold", color=INK, pad=18)
a2 = fig.add_subplot(gsp[1]); panel_profile(a2, "covid"); a2.set_title("Gene expression profile — COVID   (grey numbers = % of cells with ≥1 UMI, K=1,000 | K=4,000)", loc="left", fontsize=12, fontweight="bold", color=INK, pad=18)

gsr = gs[1:, 2].subgridspec(2, 1, hspace=0.42, height_ratios=[1.35, 1])
panel_annotation(fig.add_subplot(gsr[0])); panel_comp(fig.add_subplot(gsr[1]))
for ext in ("png", "pdf", "svg"): fig.savefig(f"{FIG}/pbisc_file_components.{ext}", dpi=300)
plt.show()'''))
C.append(nbf.v4.new_markdown_cell("## 4. Each panel as its own file (for re-arranging in PowerPoint)"))
C.append(nbf.v4.new_code_cell(r'''def save(fn, size, draw):
    f, a = plt.subplots(figsize=size); draw(a); f.tight_layout()
    for ext in ("png", "svg"): f.savefig(f"{FIG}/{fn}.{ext}", dpi=300)
    plt.close(f)
save("panel_genes", (5.2, 2.8), panel_genes); save("panel_cells", (9, 2.8), panel_cells)
for c in ("sepsis", "covid"):
    def dr(a, c=c):
        panel_profile(a, c); a.set_title(f"Gene expression profile — {c.capitalize()}  (left half K=1,000, right half K=4,000; grey = % cells detected)", loc="left", fontsize=12, fontweight="bold", color=INK, pad=18)
    save(f"panel_profile_{c}", (11, 3.6), dr)
save("panel_depth", (6, 3.4), panel_depth); save("panel_annotation", (5.6, 4.4), panel_annotation); save("panel_composition", (6.4, 3.0), panel_composition := panel_comp)
sorted(os.listdir(FIG))'''))
C.append(nbf.v4.new_markdown_cell("## 5. Save the sampled cells as h5ad"))
C.append(nbf.v4.new_code_cell(r'''obs = cells[["dataset", "draw", "sample_id", "group", "celltype_l1_pred", "aifi_l1", "aifi_l2", "total_umi", "n_genes"]].astype({c: "category" for c in ["dataset", "draw", "sample_id", "group", "celltype_l1_pred", "aifi_l1", "aifi_l2"]})
obs.index = [f"c{i}" for i in range(len(obs))]
a = ad.AnnData(X=cells[GENES].to_numpy(np.float32), obs=obs, var=pd.DataFrame(index=GENES))
a.uns["note"] = "Raw UMI counts of 19 marker genes for 4,000 random cells per draw x 5 draws per dataset (PBISC bulk-only, 2026-10-06). total_umi = all 6,000 panel genes."
a.write_h5ad(f"{DAT}/file_component_marker_cells.h5ad"); a'''))
nb=nbf.v4.new_notebook(); nb.cells=C
nb.metadata["kernelspec"]={"name":"python3","display_name":"Python 3","language":"python"}
nbf.write(nb,'/home/claude/fcfig/file_component_figure.ipynb'); print('ok')
