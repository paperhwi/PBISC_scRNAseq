"""Builds macrophage_labeling_figs.ipynb: figures for the macrophage-labelling deck (2026-10-17)."""
import nbformat as nbf
nb = nbf.v4.new_notebook(); C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Sepsis macrophage labelling — figures for the deck (2026-10-17)
Reads stored results only (no new label transfer):
- `proportion_qualification/results/`: `sepsis_realT2_mapping.h5ad`, `sepsis_k{1000,4000}_mapping.h5ad` (per-cell R1–R4 labels), `sepsis_confusion_recall.csv`
- `macrophage_markers/macrophage_marker_counts.h5ad` (40 in-panel macrophage / monocyte genes; generated groups from R3 fine labels)
- `myeloid_subcluster_20261016/data/sepsis_{real,k1000,k4000}_myeloid_subcluster.h5ad` (Leiden r = 1.0 sub-clusters)

| route | how | where it was used |
|---|---|---|
| R1 | Harmony (batch = source) on reference PCA, then kNN (k = 30, 1/d weights) from real reference cells | joint co-embedding, composition, claims |
| R2 | same kNN in reference PCA, no Harmony | robustness check |
| R3 | multinomial logistic regression on 2,000 reference HVGs, trained on reference cells passed through the PBISC-like read-out; no integration | proportion qualification, NK/ILC, monocytes, macrophages |
| R4 | marker gating, six lineages only | robustness check |
Evaluation: real samples split in half (A = reference, 39,600 cells; B = held-out, 40,200 cells, no shared donor); B passed through the PBISC-like read-out = `realT2` (labels known).""")
code(r'''import os, warnings, logging, json, numpy as np, pandas as pd, h5py, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from anndata.io import read_elem
from IPython.display import Image, display
warnings.filterwarnings("ignore"); logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
BASE = os.environ.get("PBISC_BASE", "F:/PBISC"); REP = f"{BASE}/bulkonly_lora5_20261006/analysis/integration_report_20261012"
PQ = f"{REP}/proportion_qualification/results"; MM = f"{REP}/macrophage_markers"; SUB = f"{REP}/myeloid_subcluster_20261016/data"
OUT = f"{REP}/macrophage_labeling_20261017"; F = f"{OUT}/figs"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 200, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False})
INK, MUTED = "#1f2430", "#5b6170"
COL = {"CD14 mono": "#2a78d6", "CD16 / int. mono": "#1baf7a", "Macrophage": "#e34948", "DC": "#9a9a9a", "pDC": "#c9c9c9"}
R1C, R3C, TRUTH = "#2a78d6", "#eb6834", "#4a4f5c"
def show(p, w=1000): display(Image(p, width=w))''')

md("## 1. Myeloid cells and their original labels")
code(r'''def load_obs(n):
    with h5py.File(f"{SUB}/sepsis_{n}_myeloid_subcluster.h5ad", "r") as h: return read_elem(h["obs"]), h["obsm/X_umap"][:]
def coarse_real(x):
    if x == "Macrophage": return "Macrophage"
    if x in ("Mono_c07-CD16", "Mono_c08-CD14CD16"): return "CD16 / int. mono"
    if x == "cDC": return "DC"
    if x == "pDC": return "pDC"
    return "CD14 mono"
def coarse_aifi(x):
    if x == "CD14 monocyte": return "CD14 mono"
    if x in ("CD16 monocyte", "Intermediate monocyte"): return "CD16 / int. mono"
    if x == "pDC": return "pDC"
    return "DC"
OBS = {}
for n in ["real", "k1000", "k4000"]:
    o, U = load_obs(n); o["coarse"] = (o.celltype_ye2026.map(coarse_real) if n == "real" else o.aifi_l2.map(coarse_aifi)).astype(str); OBS[n] = (o, U)
    print(n, o.coarse.value_counts().to_dict())
TITLE = {"real": "Real (Ye 2026 labels)", "k1000": "PBISC K = 1,000 (AIFI L2)", "k4000": "PBISC K = 4,000 (AIFI L2)"}
fig, axs = plt.subplots(1, 3, figsize=(18, 6.4))
for ax, n in zip(axs, OBS):
    o, U = OBS[n]; lab = o.coarse.to_numpy(); rng = np.random.default_rng(0); idx = rng.permutation(len(U))
    idx = np.concatenate([idx[lab[idx] != "Macrophage"], idx[lab[idx] == "Macrophage"]])     # macrophages drawn last
    ax.scatter(U[idx, 0], U[idx, 1], c=[COL[x] for x in lab[idx]], s=0.25 if n == "real" else 0.8, linewidths=0, rasterized=True)
    vc = o.coarse.value_counts()
    ax.set_title(f"{TITLE[n]}\n{len(o):,} myeloid cells · macrophage: {int(vc.get('Macrophage', 0)):,}", loc="left", color=INK)
    ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, label=k) for k, c in COL.items()]
fig.legend(handles=h, loc="lower center", ncol=5, fontsize=13, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(f"{F}/f1_myeloid_umap_labels.png"); plt.close(fig); show(f"{F}/f1_myeloid_umap_labels.png")''')

md("## 2. R1 vs R3 on real held-out cells (`realT2`): where do the 239 true macrophages go?")
code(r'''rt = ad.read_h5ad(f"{PQ}/sepsis_realT2_mapping.h5ad").obs
def dest(x):
    x = str(x)
    if x in ("Macrophage", "Macrophages"): return "Macrophage"
    if x in ("Mono_c07-CD16", "Mono_c08-CD14CD16"): return "CD16 / int. mono"
    if x.startswith("Mono_") or x == "Monocytes": return "CD14 mono"
    if x in ("cDC", "pDC", "DC"): return "DC"
    return "other"
tm = rt[rt.true_fine == "Macrophage"]; cats = ["Macrophage", "CD16 / int. mono", "CD14 mono", "DC", "other"]
D = pd.DataFrame({f"{r} {lv}": tm[f"{r}_{lv}"].map(dest).value_counts().reindex(cats).fillna(0).astype(int) for lv in ["fine", "major"] for r in ["R1", "R3"]}).T
cr = pd.read_csv(f"{PQ}/sepsis_confusion_recall.csv"); cr = cr[cr.label.isin(["Macrophage", "Macrophages"])].set_index(["route", "level"])
D["recall"] = [cr.loc[(i.split()[0], i.split()[1]), "recall"] for i in D.index]; D["precision"] = [cr.loc[(i.split()[0], i.split()[1]), "precision"] for i in D.index]
D["predicted_macrophage"] = [int((rt[f"{i.split()[0]}_{i.split()[1]}"].map(dest) == "Macrophage").sum()) for i in D.index]
D.to_csv(f"{F}/t2_realT2_true_macrophage_destination.csv"); display(D)
fig, ax = plt.subplots(figsize=(15, 4.8)); rows = ["R1 fine", "R3 fine", "R1 major", "R3 major"]; y = np.arange(len(rows))[::-1]
c2 = {**COL, "other": "#e6e6e6"}
for yi, r in zip(y, rows):
    left = 0
    for c in cats:
        v = D.loc[r, c]
        if v: ax.barh(yi, v, left=left, color=c2[c], height=0.62, edgecolor="white", linewidth=2)
        if v >= 14: ax.text(left + v / 2, yi, f"{v}", ha="center", va="center", fontsize=13, color="white" if c not in ("other", "DC") else INK, weight="bold")
        left += v
    ax.text(242, yi, f"recall {D.loc[r, 'recall']:.2f} · precision {D.loc[r, 'precision']:.2f}", va="center", fontsize=13, color=INK)
ax.set_yticks(y, rows, fontsize=14); ax.set_xlim(0, 330); ax.set_xticks([0, 60, 120, 180, 239]); ax.set_xlabel("239 real macrophages in the held-out half, PBISC-like read-out (realT2)")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c2[c]) for c in cats[:4]], labels=["→ Macrophage (correct)", "→ CD16 / int. mono", "→ CD14 mono", "→ DC"], ncol=4, loc="upper center", bbox_to_anchor=(0.4, 1.17), fontsize=12)
plt.tight_layout(); fig.savefig(f"{F}/f2_realT2_macrophage_destination.png"); plt.close(fig); show(f"{F}/f2_realT2_macrophage_destination.png")''')

md("## 3. R1 vs R3 on generated cells: how many cells are called macrophage?")
code(r'''rows = [{"set": "Real held-out", "what": "truth", "n": int((rt.true_fine == "Macrophage").sum()), "cells": len(rt)}]
for r in ["R1", "R3"]: rows.append({"set": "Real held-out", "what": r, "n": int((rt[f"{r}_fine"] == "Macrophage").sum()), "cells": len(rt)})
for K, nm in [("k1000", "PBISC K = 1,000"), ("k4000", "PBISC K = 4,000")]:
    o = ad.read_h5ad(f"{PQ}/sepsis_{K}_mapping.h5ad").obs
    for r in ["R1", "R3"]:
        m = o[f"{r}_fine"] == "Macrophage"
        rows.append({"set": nm, "what": r, "n": int(m.sum()), "n_not_rejected": int((m & ~o[f"{r}_fine_rejected"].astype(bool)).sum()), "cells": len(o)})
    if K == "k1000": o1 = o
    else: o4 = o
G = pd.DataFrame(rows); G["pct"] = 100 * G.n / G.cells; G.to_csv(f"{F}/t3_macrophage_calls_by_route.csv", index=False); display(G)
fig, ax = plt.subplots(figsize=(13, 5.4)); sets = ["Real held-out", "PBISC K = 1,000", "PBISC K = 4,000"]; x0 = 0; ticks = []
for s in sets:
    g = G[G.set == s]
    for _, r in g.iterrows():
        c = {"truth": TRUTH, "R1": R1C, "R3": R3C}[r.what]
        ax.bar(x0, r.pct, color=c, width=0.8); ax.text(x0, r.pct + 0.012, f"{r.n:,}\n({r.pct:.2f}%)" if r.pct >= 0.01 or r.n == 0 else f"{r.n:,}\n(<0.01%)", ha="center", va="bottom", fontsize=12, color=INK)
        ticks.append((x0, r.what)); x0 += 1
    ax.text(x0 - len(g) / 2 - 0.5, -0.16, s, ha="center", va="top", fontsize=14, weight="bold", color=INK, transform=ax.get_xaxis_transform())
    x0 += 0.8
ax.set_xticks([t[0] for t in ticks], [t[1] for t in ticks], fontsize=13); ax.set_ylim(0, 0.85); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
ax.axhline(0, color=INK, lw=0.8); ax.set_ylabel("% of cells labelled Macrophage (fine)")
plt.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(f"{F}/f3_macrophage_calls_by_route.png"); plt.close(fig); show(f"{F}/f3_macrophage_calls_by_route.png")
X = pd.DataFrame({K: o.loc[o.R3_fine == "Macrophage", "R1_fine"].map(dest).value_counts() for K, o in [("k1000", o1), ("k4000", o4)]}).fillna(0).astype(int)
print("generated cells called Macrophage by R3 — what R1 called them:"); display(X)
X.to_csv(f"{F}/t3b_R3_macrophage_R1_label.csv")''')

md("## 4. Macrophage genes by label group (R3 groups for generated cells)\nDetection = fraction of cells with count > 0 (robust to depth). Colour = log2(detection in the group / detection in the same dataset's CD14 monocytes): how specific the gene is to the group.")
code(r'''A = ad.read_h5ad(f"{MM}/macrophage_marker_counts.h5ad"); A = A[A.obs.cohort == "sepsis"]
GENES = ["C1QA", "C1QB", "C1QC", "VSIG4", "CD163", "MSR1", "MAF", "STAB1", "MERTK", "SIGLEC1", "SELENOP", "APOE", "APOC1", "MRC1", "F13A1",
         "CSF1R", "HMOX1", "MAFB", "CD68", "CD14", "LYZ", "S100A8", "VCAN", "FCN1"]
GR = ["Macrophage", "CD16 mono", "CD14 mono"]; SRC = [("real", "Real"), ("k1000", "K = 1,000"), ("k4000", "K = 4,000")]
rows, det, lr = [], [], []
for s, sn in SRC:
    b = A[A.obs.source == s]; Xb = b[:, GENES].X; ref = None
    dd = {g: np.asarray((Xb[(b.obs.group == g).to_numpy()] > 0).mean(0)).ravel() for g in GR}
    for g in GR:
        rows.append(f"{sn} · {g} (n = {int((b.obs.group == g).sum()):,})"); det.append(dd[g]); lr.append(np.log2((dd[g] + 0.01) / (dd["CD14 mono"] + 0.01)))
det, lr = np.array(det), np.clip(np.array(lr), -4, 4)
pd.DataFrame(det, index=rows, columns=GENES).round(3).to_csv(f"{F}/t4_detection_by_label_group.csv")
fig, ax = plt.subplots(figsize=(17, 6.8))
for i in range(len(rows)):
    sc_ = ax.scatter(np.arange(len(GENES)), [i] * len(GENES), s=det[i] * 420 + 2, c=lr[i], cmap="RdBu_r", vmin=-4, vmax=4, edgecolors="#b0b0b0", linewidths=0.4)
for k in [2.5, 5.5]: ax.axhline(k, color="#c8c8c8", lw=1)
ax.axvline(14.5, color="#c8c8c8", lw=1, ls="--"); ax.axvline(18.5, color="#c8c8c8", lw=1, ls="--")
ax.set_xticks(range(len(GENES)), GENES, rotation=60, ha="right", fontsize=13); ax.set_yticks(range(len(rows)), rows, fontsize=13); ax.invert_yaxis()
ax.set_xlim(-0.7, len(GENES) - 0.3)
cb = plt.colorbar(sc_, ax=ax, fraction=0.025, pad=0.01); cb.set_label("log2 detection vs CD14 mono", fontsize=12)
for v in [0.2, 0.5, 1.0]: ax.scatter([], [], s=v * 420 + 2, color="#9a9a9a", label=f"{int(v * 100)}%")
ax.legend(title="detected", loc="upper left", bbox_to_anchor=(1.07, 0.35), fontsize=12, title_fontsize=12, labelspacing=1.4)
ax.text(7, -1.1, "macrophage programme", ha="center", fontsize=13, color=MUTED); ax.text(16.5, -1.1, "shared myeloid", ha="center", fontsize=13, color=MUTED); ax.text(21, -1.1, "monocyte", ha="center", fontsize=13, color=MUTED)
for s in ax.spines.values(): s.set_visible(False)
plt.tight_layout(); fig.savefig(f"{F}/f4_dotplot_label_groups.png", bbox_inches="tight"); plt.close(fig); show(f"{F}/f4_dotplot_label_groups.png")''')

md("## 5. Macrophage genes by sub-cluster (Leiden r = 1.0)\nPer dataset: the 3 monocyte/macrophage clusters with the highest mean P(macrophage) (panel-only real-data classifier) plus the 2 largest CD14 clusters (≥ 300 cells; DC/pDC clusters excluded). Colour = mean log CP10k, scaled 0–1 per gene across ALL rows of the three datasets (so colours compare across datasets); size = detection.")
code(r'''GS = ["C1QA", "C1QB", "C1QC", "VSIG4", "CD163", "NR1H3", "MSR1", "MAF", "STAB1", "CMKLR1", "SELENOP", "MRC1", "F13A1", "FCGR3A", "MS4A7", "CSF1R", "CD14", "LYZ", "S100A8", "VCAN"]
ANN = {n: pd.read_csv(f"{REP}/myeloid_subcluster_20261016/results/{n}_myeloid_cluster_annotation_r1.0.csv").set_index("leiden_r1.0") for n in ["real", "k1000", "k4000"]}
blocks = []
for n in ["real", "k1000", "k4000"]:
    with h5py.File(f"{SUB}/sepsis_{n}_myeloid_subcluster.h5ad", "r") as h:
        var = read_elem(h["var"]); X = read_elem(h["X"]); o = read_elem(h["obs"])
    gi = [var.index.get_loc(g) for g in GS]; Xg = X[:, gi].toarray(); del X
    a = ANN[n]; keep = [c for c in a.index if a.loc[c, "cells"] >= 300 and a.loc[c, "majority_label"] not in ("cDC", "pDC", "cDC2", "cDC1", "ASDC")]
    top = sorted(keep, key=lambda c: -a.loc[c, "macro_prob_mean"])[:3]
    cd14 = [c for c in sorted(keep, key=lambda c: -a.loc[c, "cells"]) if c not in top and "CD16" not in str(a.loc[c, "majority_label"]) and "Intermediate" not in str(a.loc[c, "majority_label"])][:2]
    keep = top + cd14
    cl = o["leiden_r1.0"].astype(str).to_numpy()
    M = np.array([Xg[cl == str(c)].mean(0) for c in keep]); P = np.array([(Xg[cl == str(c)] > 0).mean(0) for c in keep])
    Ms = M
    short = lambda s: s.replace("Mono_", "").replace("monocyte", "mono").replace("Intermediate", "Int.")
    labs = [f"{n.replace('real', 'Real').replace('k1000', 'K1k').replace('k4000', 'K4k')} c{c} · {short(a.loc[c, 'majority_label'])} ({a.loc[c, 'cells']:,})" for c in keep]
    blocks.append((n, labs, Ms, P)); pd.DataFrame(P, index=labs, columns=GS).round(3).to_csv(f"{F}/t5_detection_subclusters_{n}.csv"); pd.DataFrame(M, index=labs, columns=GS).round(3).to_csv(f"{F}/t5_mean_logcp10k_subclusters_{n}.csv")
allM = np.vstack([b[2] for b in blocks]); lo, hi = allM.min(0), allM.max(0)
blocks = [(n, labs, (M - lo) / (hi - lo + 1e-9), P) for n, labs, M, P in blocks]
nrow = sum(len(b[1]) for b in blocks)
fig, ax = plt.subplots(figsize=(17, 0.5 * nrow + 2.6)); i0 = 0
for n, labs, Ms, P in blocks:
    for i in range(len(labs)):
        sc_ = ax.scatter(np.arange(len(GS)), [i0 + i] * len(GS), s=P[i] * 300 + 2, c=Ms[i], cmap="Reds", vmin=0, vmax=1, edgecolors="#b0b0b0", linewidths=0.4)
    i0 += len(labs); ax.axhline(i0 - 0.5, color="#c8c8c8", lw=1)
ax.set_yticks(range(nrow), [l for b in blocks for l in b[1]], fontsize=13); ax.invert_yaxis()
ax.set_xticks(range(len(GS)), GS, rotation=60, ha="right", fontsize=13); ax.set_xlim(-0.7, len(GS) - 0.3)
ax.axvline(12.5, color="#c8c8c8", lw=1, ls="--"); ax.axvline(15.5, color="#c8c8c8", lw=1, ls="--")
cb = plt.colorbar(sc_, ax=ax, fraction=0.02, pad=0.01); cb.set_label("mean log CP10k, scaled per gene across datasets", fontsize=12)
for v in [0.2, 0.5, 1.0]: ax.scatter([], [], s=v * 300 + 2, color="#9a9a9a", label=f"{int(v * 100)}%")
ax.legend(title="detected", loc="upper left", bbox_to_anchor=(1.08, 0.3), fontsize=12, title_fontsize=12, labelspacing=1.3)
for s in ax.spines.values(): s.set_visible(False)
plt.tight_layout(); fig.savefig(f"{F}/f5_dotplot_subclusters.png", bbox_inches="tight"); plt.close(fig); show(f"{F}/f5_dotplot_subclusters.png")
# compact numbers for the slide: best macrophage-like cluster per dataset
for n, labs, Ms, P in blocks: print(n, labs[0], dict(zip(GS[:6], P[0].round(2))))''')
nb["cells"] = C; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
import os; os.makedirs("/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/macrophage_labeling_20261017", exist_ok=True)
nbf.write(nb, "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/macrophage_labeling_20261017/macrophage_labeling_figs.ipynb")
