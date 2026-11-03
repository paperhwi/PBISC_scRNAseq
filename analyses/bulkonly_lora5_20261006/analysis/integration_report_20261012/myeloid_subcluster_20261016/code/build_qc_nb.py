import nbformat as nbf
nb = nbf.v4.new_notebook(); C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Myeloid sub-clusters at Leiden r = 1.0, and checks on the saved .h5ad files
Reads the three `data/sepsis_<name>_myeloid_subcluster.h5ad` files written by `code/write_myeloid_h5ad.py`
(from the counts and the results of `myeloid_subcluster.ipynb`).
1. File checks: X equals scanpy log-CP10k of `layers['counts']`, embeddings present.
2. r = 1.0 cluster tables (finer than the r = 0.5 tables in `myeloid_subcluster.ipynb`): does each dataset have a macrophage cluster?""")
code(r'''import os, warnings, logging, scipy.sparse as sp, numpy as np, pandas as pd, scanpy as sc, anndata as ad, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from anndata.io import read_elem
from IPython.display import Image, display
warnings.filterwarnings("ignore"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
BASE = os.environ.get("PBISC_BASE", "F:/PBISC"); W = f"{BASE}/bulkonly_lora5_20261006/analysis/integration_report_20261012/myeloid_subcluster_20261016"
D, R, F = f"{W}/data", f"{W}/results", f"{W}/figs"
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 250, "pdf.fonttype": 42})
def load(name):   # everything except layers['counts'] (keeps memory down for the 135k x 37k real file)
    with h5py.File(f"{D}/sepsis_{name}_myeloid_subcluster.h5ad", "r") as h:
        return ad.AnnData(X=read_elem(h["X"]), obs=read_elem(h["obs"]), var=read_elem(h["var"]), obsm=read_elem(h["obsm"]), uns=read_elem(h["uns"]))
A = {n: load(n) for n in ["real", "k1000", "k4000"]}
for n, a in A.items(): print(n, a.shape, list(a.obsm), "| params:", dict(a.uns["subcluster_params"]))''')
md("## 1. File checks")
code(r'''rng = np.random.default_rng(0); rows = []
for n, a in A.items():
    i = np.sort(rng.choice(a.n_obs, 3000, replace=False))
    with h5py.File(f"{D}/sepsis_{n}_myeloid_subcluster.h5ad", "r") as h:      # read only the 3,000 sampled rows of the counts layer
        g = h["layers/counts"]; ip = g["indptr"][:]; dd = [g["data"][ip[r]:ip[r + 1]] for r in i]; ii = [g["indices"][ip[r]:ip[r + 1]] for r in i]
    cnt = sp.csr_matrix((np.concatenate(dd), np.concatenate(ii), np.concatenate([[0], np.cumsum([len(x) for x in dd])])), shape=(len(i), a.n_vars))
    t = ad.AnnData(cnt.copy()); sc.pp.normalize_total(t, target_sum=1e4); sc.pp.log1p(t)
    rows.append({"dataset": n, "cells": a.n_obs, "genes": a.n_vars, "counts_integer": bool(np.all(np.mod(cnt.data, 1) == 0)),
                 "max|X - lognorm(counts)| (3k cells)": float(abs(a.X[i] - t.X).max()), "UMAP finite": bool(np.isfinite(a.obsm["X_umap"]).all()),
                 "HVG": int(a.var.highly_variable.sum()), "Leiden r0.5 / r1.0": f"{a.obs['leiden_r0.5'].nunique()} / {a.obs['leiden_r1.0'].nunique()}"}); del cnt
pd.DataFrame(rows)''')
md("## 2. r = 1.0 clusters")
code(r'''MAC = ["C1QA", "C1QB", "C1QC", "VSIG4", "CD163", "NR1H3", "MSR1", "MAF", "STAB1", "SELENOP", "APOE", "FCGR3A", "CD14", "CD1C", "LILRA4"]
LAB = {"real": "celltype_ye2026", "k1000": "aifi_l2", "k4000": "aifi_l2"}; TAB = {}
for n, a in A.items():
    k = "leiden_r1.0"; ct = pd.crosstab(a.obs[k], a.obs[LAB[n]]); g = [x for x in MAC if x in a.var_names]
    det = pd.DataFrame({x: [float((a[a.obs[k] == c, x].X > 0).mean()) for c in ct.index] for x in g}, index=ct.index).round(3)
    t = pd.DataFrame({"cells": ct.sum(1), "majority_label": ct.idxmax(1), "purity": (ct.max(1) / ct.sum(1)).round(3),
                      "macro_prob_mean": a.obs.groupby(k).macro_prob_panelLR.mean().round(3),
                      "frac_macro_prob>=0.9": a.obs.assign(h=a.obs.macro_prob_panelLR >= 0.9).groupby(k).h.mean().round(3)})
    if n == "real":
        m = ct.get("Macrophage", 0); t["macrophage_share"] = (m / ct.sum(1)).round(3); t["share_of_all_macrophages"] = (m / m.sum()).round(3)
    t = t.join(det.add_prefix("det_")); t.index.name = k; t.to_csv(f"{R}/{n}_myeloid_cluster_annotation_r1.0.csv"); TAB[n] = t
    print(f"== {n}"); display(t.sort_values("macro_prob_mean", ascending=False).head(6))''')
code(r'''# best macrophage candidate cluster per dataset vs the real macrophage cluster
best = {n: t.sort_values("macro_prob_mean", ascending=False).index[0] for n, t in TAB.items()}
cmp_ = pd.DataFrame({f"{n} cluster {c}": TAB[n].loc[c] for n, c in best.items()}).T
cmp_[["cells", "majority_label", "macro_prob_mean", "frac_macro_prob>=0.9"] + [f"det_{g}" for g in MAC if f"det_{g}" in cmp_.columns]]''')
code(r'''HALO = [pe.withStroke(linewidth=2.2, foreground="white")]; TABC = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
TITLE = {"real": "Real (GSE279452)", "k1000": "PBISC K = 1,000", "k4000": "PBISC K = 4,000"}
fig, axs = plt.subplots(2, 3, figsize=(19, 12))
for j, (n, a) in enumerate(A.items()):
    U = a.obsm["X_umap"]; lab = a.obs["leiden_r1.0"].astype(str).to_numpy(); cats = list(a.obs["leiden_r1.0"].cat.categories); col = {c: TABC[i % 40] for i, c in enumerate(cats)}
    o = np.random.default_rng(0).permutation(a.n_obs); s = 0.2 if n == "real" else 0.7
    ax = axs[0, j]; ax.scatter(U[o, 0], U[o, 1], c=[col[x] for x in lab[o]], s=s, linewidths=0, rasterized=True)
    for c in cats:
        m = lab == c
        if m.sum() > 30: ax.text(*np.median(U[m], 0), c, fontsize=8, ha="center", va="center", path_effects=HALO, weight="bold" if c == best[n] else "normal")
    ax.set_title(f"{TITLE[n]} myeloid — Leiden r = 1.0 ({len(cats)} clusters)\nbold = highest mean P(macrophage): cluster {best[n]}", loc="left", fontsize=9.5)
    ax = axs[1, j]; v = a.obs.macro_prob_panelLR.to_numpy(); oo = np.argsort(v)
    sc_ = ax.scatter(U[oo, 0], U[oo, 1], c=v[oo], s=s, cmap="Reds", vmin=0, vmax=1, linewidths=0, rasterized=True)
    ax.set_title(f"{TITLE[n]}: P(macrophage), panel-only real-data classifier\ncells with P >= 0.9: {(v >= 0.9).mean():.2%}", loc="left", fontsize=9.5)
    if j == 2: plt.colorbar(sc_, ax=ax, fraction=0.03)
    for ax in axs[:, j]: ax.set_xticks([]); ax.set_yticks([]); [x.set_visible(False) for x in ax.spines.values()]
plt.tight_layout(); fig.savefig(f"{F}/myeloid_subcluster_r1.0_macrophage.png"); fig.savefig(f"{F}/myeloid_subcluster_r1.0_macrophage.pdf"); plt.close(fig)
display(Image(f"{F}/myeloid_subcluster_r1.0_macrophage.png", width=1100))''')
nb["cells"] = C; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/myeloid_subcluster_20261016/myeloid_subcluster_r1_checks.ipynb")
