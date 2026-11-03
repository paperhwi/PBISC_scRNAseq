"""Builds myeloid_subcluster.ipynb (executed afterwards with nbconvert)."""
import nbformat as nbf
nb = nbf.v4.new_notebook(); C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))

md("""# Sepsis myeloid cells: why PBISC has no macrophage label, and Mono / Macrophage / DC sub-clustering
*2026-10-16 · Ye 2026 sepsis (GSE279452) vs PBISC bulk-only LoRA-5, K = 1,000 and K = 4,000*

**Part 1 — is the missing macrophage label caused by the gene panel?** PBISC generates 5,975 genes. We
(a) find macrophage markers *from the real data* (DE, not a curated list) and check how many sit in the panel;
(b) ask whether real macrophages can still be told apart from monocytes/DC **using only panel genes**;
(c) apply that panel-only real-data classifier to generated cells; (d) check what each labelling route can output.

**Part 2 — sub-clustering.** Each dataset's myeloid compartment is re-processed on its own
(HVG → PCA → Harmony → kNN → Leiden → UMAP; batch = sample for real, draw for generated).

| dataset | myeloid cells | definition | genes |
|---|---|---|---|
| real | Ye 2026 major ∈ {Monocytes, Macrophages, DC}, in paper scheme | published labels | 37,474 |
| K = 1,000 / 4,000 | PBISC `celltype_l1_pred` ∈ {Monocyte, DC} in the Harmony-integrated 305k-cell set | model labels | 5,975 panel |

Panel membership is matched on **Ensembl IDs** (113 panel genes carry older symbols, e.g. FAM46C → TENT5C).""")

code(r'''import os, json, warnings, logging, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from scipy.stats import hypergeom
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score, recall_score
from IPython.display import Image, display
warnings.filterwarnings("ignore"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
BASE = os.environ.get("PBISC_BASE", "F:/PBISC"); REP = f"{BASE}/bulkonly_lora5_20261006/analysis/integration_report_20261012"
W = f"{REP}/myeloid_subcluster_20261016"; E = f"{W}/_extract"; D = f"{W}/data"; R = f"{W}/results"; F = f"{W}/figs"
for p in (R, F): os.makedirs(p, exist_ok=True)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "axes.titlesize": 10, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "savefig.dpi": 250, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
sc.settings.verbosity = 0; RNG = np.random.default_rng(0)
def show(p, w=1000): display(Image(p, width=w))''')

md("## Load count matrices and the panel")
code(r'''real = ad.read_h5ad(f"{D}/sepsis_real_myeloid_counts.h5ad")
gen = {K: ad.read_h5ad(f"{D}/sepsis_{K}_myeloid_counts.h5ad") for K in ["k1000", "k4000"]}
real.var["gene_id"] = np.load(f"{E}/real_gene_ids.npy", allow_pickle=True).astype(str)
pid = np.load(f"{E}/panel_gene_ids.npy", allow_pickle=True).astype(str)
for K in gen: gen[K].var["gene_id"] = pid
real.var["in_panel"] = real.var.gene_id.isin(set(pid))
common = [g for g in pid if g in set(real.var.gene_id)]          # 5,967 (the 8 MT genes were removed from the real data)
print(real.shape, {K: a.shape for K, a in gen.items()}, "| panel genes:", len(pid), "| shared with real:", len(common))
print(real.obs.celltype_major_ye2026.value_counts().to_dict())
for K in gen: print(K, gen[K].obs.celltype_l1_pred.value_counts().to_dict())''')

md("## Part 1a — macrophage markers found in the real data")
code(r'''n = real.n_obs; det = real.X.getnnz(axis=0) / n; real.var["det_myeloid"] = det      # real.X = raw counts throughout Part 1
lab = real.obs.celltype_ye2026.astype(str).to_numpy()
grp = np.where(lab == "Macrophage", "Macrophage", np.where(np.isin(lab, ["Mono_c07-CD16", "Mono_c08-CD14CD16"]), "CD16/int mono",
               np.where(np.isin(lab, ["cDC", "pDC"]), "DC", "CD14 mono")))
real.obs["mgroup"] = pd.Categorical(grp, ["Macrophage", "CD14 mono", "CD16/int mono", "DC"])
print(real.obs.mgroup.value_counts().to_dict())
def lognorm_subset(a, rows):   # log CP10k on ALL genes for a row subset (keeps the full matrix as counts)
    s = a[np.sort(rows)].copy(); sc.pp.normalize_total(s, target_sum=1e4); sc.pp.log1p(s); return s
# DE on a subsample: all macrophages + up to 6,000 of each other group; genes detected in >= 1% of myeloid cells
keep_g = det >= 0.01
idx = np.concatenate([RNG.permutation(np.where(grp == g)[0])[: (None if g == "Macrophage" else 6000)] for g in ["Macrophage", "CD14 mono", "CD16/int mono", "DC"]])
de = lognorm_subset(real, idx)[:, keep_g].copy()
de.obs["vs"] = np.where(de.obs.mgroup == "Macrophage", "Macrophage", "rest")
sc.tl.rank_genes_groups(de, "vs", groups=["Macrophage"], reference="rest", method="wilcoxon", pts=True, key_added="all")
sc.tl.rank_genes_groups(de, "mgroup", groups=["Macrophage"], reference="CD16/int mono", method="wilcoxon", pts=True, key_added="c16")
def tab(key):
    t = sc.get.rank_genes_groups_df(de, "Macrophage", key=key).set_index("names")
    return t[["scores", "logfoldchanges", "pvals_adj"]]
T = tab("all").join(tab("c16").add_suffix("_vsCD16"))
pct = lambda m: pd.Series(np.asarray((de.X[m] > 0).mean(0)).ravel(), index=de.var_names)
T["pct_nz_group"] = pct((de.obs.mgroup == "Macrophage").to_numpy()); T["pct_nz_reference"] = pct((de.obs.mgroup != "Macrophage").to_numpy())
T["pct_CD16int_mono"] = pct((de.obs.mgroup == "CD16/int mono").to_numpy())
T["gene_id"] = real.var.loc[T.index, "gene_id"]; T["in_panel"] = real.var.loc[T.index, "in_panel"]
mk = T[(T.logfoldchanges >= 1) & (T.pvals_adj < 0.05) & (T.pct_nz_group >= 0.25)].sort_values("scores", ascending=False)
mk.to_csv(f"{R}/real_macrophage_markers_DE.csv"); T.to_csv(f"{R}/real_macrophage_DE_all_genes.csv.gz"); del de
print(f"{len(mk)} macrophage markers (log2FC >= 1, padj < 0.05, detected in >= 25% of macrophages)")
mk.head(40)[["logfoldchanges", "pct_nz_group", "pct_nz_reference", "logfoldchanges_vsCD16", "in_panel"]].round(2)''')

md("## Part 1b — how many of those markers does the PBISC panel contain?\nBackground = all genes detected in ≥ 1% of real myeloid cells. Over-representation is tested with the hypergeometric distribution.")
code(r'''bg = real.var[keep_g]; Nbg, Kbg = len(bg), int(bg.in_panel.sum())
rows = []
for n in sorted(set([n for n in [10, 25, 50, 100] if n < len(mk)] + [len(mk)])):
    top = mk.head(n); k = int(top.in_panel.sum())
    rows.append({"top_n_markers": n, "in_panel": k, "fraction_in_panel": k / n, "background_fraction": Kbg / Nbg,
                 "enrichment": (k / n) / (Kbg / Nbg), "hypergeom_p_enriched": hypergeom.sf(k - 1, Nbg, Kbg, n)})
cov = pd.DataFrame(rows); cov.to_csv(f"{R}/panel_coverage_of_macrophage_markers.csv", index=False)
print(f"background: {Kbg:,} of {Nbg:,} expressed myeloid genes are in the panel ({Kbg / Nbg:.1%})")
missing_top = mk.head(50)[~mk.head(50).in_panel].index.tolist(); print("top-50 markers NOT in panel:", missing_top)
cov.round(4)''')
code(r'''fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 2.2]})
ax = axs[0]; x = np.arange(len(cov))
ax.bar(x, cov.fraction_in_panel * 100, color="#2a78d6", width=0.6)
ax.axhline(Kbg / Nbg * 100, ls="--", c="#555", lw=1); ax.text(len(cov) - 0.5, Kbg / Nbg * 100 + 1.5, f"background {Kbg / Nbg:.0%}", ha="right", fontsize=8, color="#555")
for xi, r in zip(x, cov.itertuples()): ax.text(xi, r.fraction_in_panel * 100 + 1.5, f"{r.in_panel}/{r.top_n_markers}", ha="center", fontsize=8)
ax.set_xticks(x, [f"top {n}" if n != len(mk) else f"all {n}" for n in cov.top_n_markers]); ax.set_ylim(0, 105)
ax.set_ylabel("% of real macrophage markers in the PBISC panel"); ax.set_title("Panel coverage of data-driven macrophage markers", loc="left")
ax = axs[1]; top = mk.head(40)
ax.bar(np.arange(len(top)), top.logfoldchanges, color=["#2a78d6" if v else "#e34948" for v in top.in_panel], width=0.75)
ax.set_xticks(np.arange(len(top)), top.index, rotation=70, ha="right", fontsize=7.5); ax.set_ylabel("log2FC, macrophage vs other myeloid")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#2a78d6"), plt.Rectangle((0, 0), 1, 1, color="#e34948")], labels=["in panel", "not in panel"], frameon=False, fontsize=8)
ax.set_title("Top 40 real macrophage markers (ranked by Wilcoxon score)", loc="left")
plt.tight_layout(); fig.savefig(f"{F}/macrophage_marker_panel_coverage.png"); plt.close(fig); show(f"{F}/macrophage_marker_panel_coverage.png")''')

md("""## Part 1c — can real macrophages be recognised with panel genes only?
Logistic regression, macrophage vs other myeloid, 5-fold cross-validation **grouped by sample** (no sample in both train and test).
- *all genes*: every gene detected in ≥ 1% of myeloid cells, log CP10k on all genes.
- *panel only*: the 5,967 panel genes, counts re-normalised on panel totals — the format PBISC cells come in.""")
code(r'''cid = np.load(f"{E}/real_rows.npy")
ci = np.concatenate([RNG.permutation(np.where(grp == g)[0])[: (None if g == "Macrophage" else 8000)] for g in ["Macrophage", "CD14 mono", "CD16/int mono", "DC"]]); ci.sort()
y = (grp[ci] == "Macrophage").astype(int); groups = real.obs["sample"].to_numpy()[ci]
gpos = pd.Series(np.arange(real.n_vars), index=real.var.gene_id)
def panel_fmt(Xc):   # counts on the shared panel genes -> log CP10k on panel totals
    tot = np.asarray(Xc.sum(1)).ravel(); Y = sp.diags(1e4 / np.maximum(tot, 1)) @ Xc; Y = sp.csr_matrix(Y, dtype=np.float32); Y.data = np.log1p(Y.data); return Y
Xr_panel_counts = real.X[:, gpos.loc[common].to_numpy()].tocsr()
FEAT = {"all genes": lognorm_subset(real, ci).X[:, np.where(keep_g)[0]], "panel only": panel_fmt(Xr_panel_counts[ci])}
cv = []; oof = {}
for name, Xf in FEAT.items():
    p = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(Xf, y, groups):
        m = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000, solver="liblinear").fit(Xf[tr], y[tr]); p[te] = m.predict_proba(Xf[te])[:, 1]
    oof[name] = p
    cv.append({"features": name, "n_genes": Xf.shape[1], "AUROC": roc_auc_score(y, p), "average_precision": average_precision_score(y, p),
               "precision@0.5": precision_score(y, p >= 0.5), "recall@0.5": recall_score(y, p >= 0.5)})
# same panel-only model, but every real myeloid cell is scored by a model that never saw its sample (natural macrophage prevalence)
Pall = panel_fmt(Xr_panel_counts); del Xr_panel_counts; yall = (grp == "Macrophage").astype(int); gall = real.obs["sample"].to_numpy(); p_all = np.zeros(real.n_obs)
for tr, te in GroupKFold(5).split(np.zeros(real.n_obs), yall, gall):
    trs = np.concatenate([RNG.permutation(tr[grp[tr] == g])[: (None if g == "Macrophage" else 8000)] for g in ["Macrophage", "CD14 mono", "CD16/int mono", "DC"]])
    m = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000, solver="liblinear").fit(Pall[trs], yall[trs]); p_all[te] = m.predict_proba(Pall[te])[:, 1]
real.obs["macro_prob_panelLR"] = p_all; del Pall
cv.append({"features": "panel only, all real myeloid cells (natural prevalence)", "n_genes": len(common), "AUROC": roc_auc_score(yall, p_all),
           "average_precision": average_precision_score(yall, p_all), "precision@0.5": precision_score(yall, p_all >= 0.5), "recall@0.5": recall_score(yall, p_all >= 0.5)})
cv = pd.DataFrame(cv); cv.to_csv(f"{R}/macrophage_classifier_cv_real.csv", index=False); cv.round(3)''')

md("## Part 1d — apply the panel-only real classifier to generated myeloid cells")
code(r'''clf = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000, solver="liblinear").fit(FEAT["panel only"], y)
coef = pd.Series(clf.coef_[0], index=real.var.index[gpos.loc[common].to_numpy()]).sort_values()
coef.to_csv(f"{R}/macrophage_classifier_panel_coefficients.csv")
gpos_g = pd.Series(np.arange(len(pid)), index=pid)
prob = {"real: labelled macrophage": p_all[yall == 1], "real: other myeloid": p_all[yall == 0]}
real_myeloid_rate = (grp == "Macrophage").mean()
rows = [{"dataset": "real", "cells": real.n_obs, "macrophage_fraction": real_myeloid_rate, "basis": "published label"},
        {"dataset": "real", "cells": real.n_obs, "macrophage_fraction": (p_all >= 0.5).mean(), "basis": "panel classifier p >= 0.5 (held-out samples)"},
        {"dataset": "real", "cells": real.n_obs, "macrophage_fraction": (p_all >= 0.9).mean(), "basis": "panel classifier p >= 0.9 (held-out samples)"}]
for K, a in gen.items():
    Pg = clf.predict_proba(panel_fmt(a.X[:, gpos_g.loc[common].to_numpy()].tocsr()))[:, 1]; a.obs["macro_prob_panelLR"] = Pg; prob[K] = Pg
    rows.append({"dataset": K, "cells": a.n_obs, "macrophage_fraction": (Pg >= 0.5).mean(), "basis": "panel classifier p >= 0.5"})
    rows.append({"dataset": K, "cells": a.n_obs, "macrophage_fraction": (Pg >= 0.9).mean(), "basis": "panel classifier p >= 0.9"})
# precision-adjusted expectation for real cells at the same threshold
pr = precision_score(y, oof["panel only"] >= 0.5); rc = recall_score(y, oof["panel only"] >= 0.5)
gp = pd.DataFrame(rows); gp.to_csv(f"{R}/generated_macrophage_prediction.csv", index=False)
print(f"real, held-out samples, p>=0.5: precision {precision_score(yall, p_all >= 0.5):.2f}, recall {recall_score(yall, p_all >= 0.5):.2f}; p>=0.9: precision {precision_score(yall, p_all >= 0.9):.2f}, recall {recall_score(yall, p_all >= 0.9):.2f}")
gp''')
code(r'''# detection of the top in-panel macrophage markers: real groups vs generated cells (generated: all myeloid, and the top 1% by classifier probability)
top_in = mk[mk.in_panel].head(20)
tg = [real.var.index.get_loc(g) for g in top_in.index]; tgid = real.var.gene_id.iloc[tg].to_numpy()
cols = {}
hi_r = p_all >= 0.9; cols[f"real classifier p>=0.9 (n={hi_r.sum()})"] = np.asarray((real.X[np.where(hi_r)[0]][:, tg] > 0).mean(0)).ravel()
for g in ["Macrophage", "CD14 mono", "CD16/int mono", "DC"]:
    cols[f"real {g}"] = np.asarray((real.X[np.where(grp == g)[0]][:, tg] > 0).mean(0)).ravel()
for K, a in gen.items():
    Xg = a.X[:, gpos_g.loc[tgid].to_numpy()]
    cols[f"{K} all myeloid"] = np.asarray((Xg > 0).mean(0)).ravel()
    hi = a.obs.macro_prob_panelLR.to_numpy() >= 0.9
    cols[f"{K} classifier p>=0.9 (n={hi.sum()})"] = np.asarray((Xg[hi] > 0).mean(0)).ravel()
detm = pd.DataFrame(cols, index=top_in.index); detm.to_csv(f"{R}/macrophage_marker_detection_real_vs_generated.csv")
fig, axs = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axs[0]; im = ax.imshow(detm.T.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=1)
ax.set_xticks(range(len(detm)), detm.index, rotation=70, ha="right", fontsize=8); ax.set_yticks(range(detm.shape[1]), detm.columns, fontsize=8)
for i in range(detm.shape[1]):
    for j in range(len(detm)): v = detm.iloc[j, i]; ax.text(j, i, f"{v:.2f}".lstrip("0") if v < 1 else "1", ha="center", va="center", fontsize=5.8, color="white" if v > 0.55 else "#1f1f1f")
plt.colorbar(im, ax=ax, fraction=0.025, label="fraction of cells detecting the gene"); ax.set_title("Top 20 in-panel macrophage markers: detection rate", loc="left")
ax = axs[1]; bins = np.linspace(0, 1, 41)
for k, c in zip(prob, ["#1f1f1f", "#9a9a9a", "#2a78d6", "#e34948"]):
    ax.hist(prob[k], bins=bins, histtype="step", lw=1.6, color=c, density=True, label=k)
ax.set_yscale("log"); ax.set_xlabel("P(macrophage) from the panel-only real-data classifier"); ax.set_ylabel("density (log)"); ax.legend(frameon=False, fontsize=8)
ax.set_title("Macrophage-likeness of generated cells", loc="left")
plt.tight_layout(); fig.savefig(f"{F}/macrophage_markers_generated_vs_real.png"); plt.close(fig); show(f"{F}/macrophage_markers_generated_vs_real.png")''')

md("## Part 1e — what each labelling route can output")
code(r'''aifi = sorted(set(gen["k1000"].obs.aifi_l2) | set(pd.read_csv(f"{REP}/results/sepsis_k1000_transfer.csv.gz", usecols=["aifi_l2"]).aifi_l2))
tr = {K: pd.read_csv(f"{REP}/results/sepsis_{K}_transfer.csv.gz", usecols=["major", "label"], keep_default_na=False) for K in ["k1000", "k4000"]}
route = pd.DataFrame([
    {"route": "PBISC native label (AIFI L1/L2, PBMC reference)", "has a macrophage class": any("acroph" in x for x in aifi), "macrophage cells K1000": 0, "macrophage cells K4000": 0},
    {"route": "kNN transfer from real (joint Harmony, integration report)", "has a macrophage class": True,
     "macrophage cells K1000": int((tr["k1000"].major == "Macrophages").sum()), "macrophage cells K4000": int((tr["k4000"].major == "Macrophages").sum())}])
print("AIFI L2 classes in generated data:", aifi); route''')

md("## Part 2 — sub-clustering of each myeloid dataset")
code(r'''CANON = ["CD14", "LYZ", "S100A8", "S100A12", "VCAN", "FCN1", "IL1B", "CXCL8", "HLA-DRA", "ISG15", "IFITM1", "RETN", "MKI67",
         "FCGR3A", "CDKN1C", "LST1", "C1QA", "C1QB", "C1QC", "CD163", "APOE", "MRC1", "SELENOP", "CD1C", "FCER1A", "CLEC10A", "CLEC9A", "LILRA4", "IL3RA", "PPBP"]
def subcluster(a, batch, n_hvg=2000, seed=0):
    """a.X = raw counts. Returns a with lognorm X, HVG, PCA, Harmony, kNN, Leiden (0.5, 1.0), UMAP."""
    try:    # batch-aware HVG; loess can fail on batches with very few cells -> fall back to pooled HVG (batch is still corrected by Harmony)
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=n_hvg, batch_key=batch, span=0.5); a.uns["hvg_batch_key"] = batch
    except ValueError as e:
        print("batch-aware HVG failed:", str(e)[:80], "-> pooled seurat_v3"); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=n_hvg, span=0.5); a.uns["hvg_batch_key"] = "none"
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    # scaled HVG matrix kept in float32 (sc.pp.scale would densify in float64 -> OOM for 135k real cells in 7 GB)
    Z = a[:, a.var.highly_variable.to_numpy()].X.toarray().astype(np.float32)
    Z -= Z.mean(0); Z /= (Z.std(0) + 1e-8); np.clip(Z, -10, 10, out=Z)
    from sklearn.decomposition import PCA
    a.obsm["X_pca"] = PCA(50, svd_solver="randomized", random_state=seed).fit_transform(Z).astype(np.float32); del Z
    import gc; gc.collect()
    import harmonypy
    ho = harmonypy.run_harmony(a.obsm["X_pca"].astype(np.float64), a.obs[[batch]], [batch], random_state=seed, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); a.obsm["X_pca_harmony"] = (Z.T if Z.shape[0] == 50 else Z).astype(np.float32)
    sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=30, random_state=seed)
    for r in [0.5, 1.0]: sc.tl.leiden(a, resolution=r, key_added=f"leiden_r{r}", flavor="igraph", n_iterations=2, directed=False, random_state=seed)
    sc.tl.umap(a, min_dist=0.3, random_state=seed); return a
def markers(a, key, per=1000, seed=0):
    idx = np.concatenate([np.random.default_rng(seed).permutation(np.where(a.obs[key] == c)[0])[:per] for c in a.obs[key].cat.categories]); s = a[np.sort(idx)].copy()
    s = s[:, s.X.getnnz(axis=0) >= 0.01 * s.n_obs].copy()
    sc.tl.rank_genes_groups(s, key, method="wilcoxon", pts=True); t = sc.get.rank_genes_groups_df(s, None); return t''')
code(r'''# Part 1 never normalised the full real matrix, so it is still raw counts here
del FEAT, oof, Xg, detm
import gc; gc.collect()
rs = real                                   # counts -> normalised in place inside subcluster()
rs = subcluster(rs, "sample"); print(rs)
gs = {K: subcluster(a, "draw") for K, a in gen.items()}; print({K: a.obs["leiden_r0.5"].nunique() for K, a in gs.items()}, rs.obs["leiden_r0.5"].nunique())''')
code(r'''SETS = {"real": (rs, "celltype_ye2026"), "k1000": (gs["k1000"], "aifi_l2"), "k4000": (gs["k4000"], "aifi_l2")}
macro_genes = [g for g in mk[mk.in_panel].head(20).index if g in rs.var_names]
for name, (a, lab) in SETS.items():
    gl = [g for g in macro_genes if g in a.var_names]; sc.tl.score_genes(a, gl, score_name="macro_score", random_state=0)
    t = markers(a, "leiden_r0.5"); t.to_csv(f"{R}/{name}_myeloid_leiden_r0.5_markers.csv.gz", index=False)
    top = t[(t.logfoldchanges > 0.5) & (t.pvals_adj < 0.05)].sort_values("scores", ascending=False).groupby("group").head(8).groupby("group")["names"].apply(lambda s: ", ".join(s))
    ct = pd.crosstab(a.obs["leiden_r0.5"], a.obs[lab]); maj = ct.idxmax(1); pur = (ct.max(1) / ct.sum(1))
    ann = pd.DataFrame({"cells": ct.sum(1), "majority_label": maj, "purity": pur.round(3), "second_label": ct.apply(lambda r: r.sort_values().index[-2] if (r > 0).sum() > 1 else "", axis=1),
                        "macro_score_mean": a.obs.groupby("leiden_r0.5").macro_score.mean().round(3),
                        "C1QA_det": [float((a[a.obs["leiden_r0.5"] == c, "C1QA"].X > 0).mean()) if "C1QA" in a.var_names else np.nan for c in ct.index],
                        "top_markers": top.reindex(ct.index)})
    if name == "real": ann["macrophage_share"] = (ct.get("Macrophage", 0) / ct.sum(1)).round(3)
    else: ann["macro_prob_panelLR_mean"] = a.obs.groupby("leiden_r0.5").macro_prob_panelLR.mean().round(3)
    ann.index.name = "leiden_r0.5"; ann.to_csv(f"{R}/{name}_myeloid_cluster_annotation.csv"); SETS[name] = (a, lab, ann)
    print(f"== {name}: {a.n_obs:,} cells, {len(ann)} clusters (r = 0.5)"); display(ann.drop(columns="top_markers"))''')
code(r'''HALO = [pe.withStroke(linewidth=2.2, foreground="white")]
TAB = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
def scat(ax, U, lab, title, s, cats=None, number=False, legend=True, fs=7):
    lab = np.asarray(lab).astype(str); cats = cats or list(pd.Series(lab).value_counts().index); col = {c: TAB[i % len(TAB)] for i, c in enumerate(cats)}
    o = RNG.permutation(len(U)); ax.scatter(U[o, 0], U[o, 1], c=[col[x] for x in lab[o]], s=s, linewidths=0, rasterized=True)
    for c in cats:
        m = lab == c
        if m.sum() >= 30: ax.text(*np.median(U[m], 0), c if number else "", fontsize=fs, ha="center", va="center", path_effects=HALO)
    if legend:
        n = pd.Series(lab).value_counts()
        ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=4.5, color=col[c], label=f"{c} ({n[c]:,})") for c in cats if c in n], loc="upper left",
                  bbox_to_anchor=(1, 1), frameon=False, fontsize=6.3, handletextpad=0.1, ncol=1 if len(cats) <= 22 else 2)
    ax.set_title(title, loc="left"); ax.set_xticks([]); ax.set_yticks([]); [sp_.set_visible(False) for sp_ in ax.spines.values()]
TITLE = {"real": "Real (GSE279452)", "k1000": "PBISC K = 1,000", "k4000": "PBISC K = 4,000"}
fig, axs = plt.subplots(3, 3, figsize=(20, 17), gridspec_kw={"wspace": 0.55, "hspace": 0.12})
for i, (name, (a, lab, ann)) in enumerate(SETS.items()):
    U = a.obsm["X_umap"]; s = 0.25 if name == "real" else 0.8
    cl = [str(c) for c in a.obs["leiden_r0.5"].cat.categories]
    scat(axs[i, 0], U, a.obs["leiden_r0.5"], f"{TITLE[name]} myeloid · {a.n_obs:,} cells\nLeiden r = 0.5", s, cats=cl, number=True)
    scat(axs[i, 1], U, a.obs[lab], f"{TITLE[name]}: {'Ye 2026 subsets' if name == 'real' else 'PBISC AIFI L2'}", s)
    ax = axs[i, 2]; o = np.argsort(a.obs.macro_score.to_numpy()); scx = ax.scatter(U[o, 0], U[o, 1], c=a.obs.macro_score.to_numpy()[o], s=s, cmap="Reds", vmin=0, vmax=np.quantile(rs.obs.macro_score, 0.995), linewidths=0, rasterized=True)
    plt.colorbar(scx, ax=ax, fraction=0.035, label="macrophage score (top-20 in-panel markers)"); ax.set_xticks([]); ax.set_yticks([]); [sp_.set_visible(False) for sp_ in ax.spines.values()]
    ax.set_title(f"{TITLE[name]}: macrophage signature", loc="left")
fig.savefig(f"{F}/myeloid_subcluster_umaps.png", bbox_inches="tight"); fig.savefig(f"{F}/myeloid_subcluster_umaps.pdf", bbox_inches="tight"); plt.close(fig); show(f"{F}/myeloid_subcluster_umaps.png", 1100)''')
code(r'''# canonical myeloid markers per cluster (genes present in all three datasets)
genes = [g for g in CANON if all(g in a.var_names for a, _, _ in SETS.values())]
fig, axs = plt.subplots(3, 1, figsize=(13, 17))
for ax, (name, (a, lab, ann)) in zip(axs, SETS.items()):
    cats = list(a.obs["leiden_r0.5"].cat.categories); M = np.zeros((len(cats), len(genes))); P = np.zeros_like(M)
    Xg = a[:, genes].X
    for i, c in enumerate(cats):
        m = (a.obs["leiden_r0.5"] == c).to_numpy(); Xi = Xg[m]; M[i] = np.asarray(Xi.mean(0)).ravel(); P[i] = np.asarray((Xi > 0).mean(0)).ravel()
    Mz = (M - M.mean(0)) / (M.std(0) + 1e-9)
    for i in range(len(cats)):
        ax.scatter(np.arange(len(genes)), [i] * len(genes), s=P[i] * 120, c=Mz[i], cmap="RdBu_r", vmin=-2, vmax=2, linewidths=0)
    ax.set_xticks(range(len(genes)), genes, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(cats)), [f"{c}  {ann.loc[c, 'majority_label']} ({ann.loc[c, 'cells']:,})" for c in cats], fontsize=7.5); ax.invert_yaxis()
    ax.set_title(f"{TITLE[name]}: canonical myeloid markers by Leiden r = 0.5 cluster (dot size = % detecting, colour = scaled mean)", loc="left", fontsize=9.5)
    ax.set_xlim(-0.7, len(genes) - 0.3)
plt.tight_layout(); fig.savefig(f"{F}/myeloid_subcluster_dotplots.png"); fig.savefig(f"{F}/myeloid_subcluster_dotplots.pdf"); plt.close(fig); show(f"{F}/myeloid_subcluster_dotplots.png", 1000)''')

md("## Save\n`results/<name>_myeloid_subcluster_*.npy|npz` hold everything computed here; `code/write_myeloid_h5ad.py` turns them plus the extracted counts into `data/sepsis_<name>_myeloid_subcluster.h5ad` (X = log CP10k, `layers['counts']`, obs, `obsm['X_pca' | 'X_pca_harmony' | 'X_umap' | 'X_umap_global']`).")
code(r'''for name, (a, lab, ann) in SETS.items():
    o = a.obs.copy()
    np.savez_compressed(f"{R}/{name}_myeloid_subcluster_obs.npz", obs_names=o.index.to_numpy().astype(str),
                        **{c: o[c].astype(str).to_numpy() for c in o.columns if o[c].dtype.name in ("category", "object")},
                        **{f"num__{c}": o[c].to_numpy(np.float32) for c in o.columns if o[c].dtype.name not in ("category", "object")})
    for k in ["X_pca", "X_pca_harmony", "X_umap"]: np.save(f"{R}/{name}_myeloid_subcluster_{k}.npy", a.obsm[k].astype(np.float32))
    np.save(f"{R}/{name}_myeloid_subcluster_hvg.npy", a.var.highly_variable.to_numpy())
    json.dump({"batch_key": "sample" if name == "real" else "draw", "hvg_batch_key": a.uns["hvg_batch_key"], "n_hvg": 2000, "hvg_flavor": "seurat_v3 (span 0.5)", "n_pcs": 50, "harmony_max_iter": 20,
               "n_neighbors": 30, "leiden": [0.5, 1.0], "umap_min_dist": 0.3, "normalize": "CP10k + log1p", "scanpy": sc.__version__},
              open(f"{R}/{name}_myeloid_subcluster_params.json", "w"), indent=1)
sorted(os.listdir(R))''')
code(r'''# sanity checks
for name, (a, lab, ann) in SETS.items():
    assert a.obsm["X_umap"].shape[0] == a.n_obs and not np.isnan(a.obsm["X_umap"]).any()
    print(name, a.n_obs, "cells |", ann.shape[0], "clusters | largest cluster share", f"{ann.cells.max() / a.n_obs:.1%}")
mac_cl = SETS["real"][2].sort_values("macrophage_share", ascending=False).head(3); print("real clusters richest in Macrophage label:"); display(mac_cl[["cells", "majority_label", "macrophage_share", "macro_score_mean", "C1QA_det"]])''')

nb["cells"] = C; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/myeloid_subcluster_20261016/myeloid_subcluster.ipynb")
