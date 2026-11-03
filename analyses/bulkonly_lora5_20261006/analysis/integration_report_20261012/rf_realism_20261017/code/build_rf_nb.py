"""Builds rf_realism.ipynb — scGAN-style real-vs-generated random forest AUC with real-real controls (2026-10-17)."""
import nbformat as nbf, os
nb = nbf.v4.new_notebook(); C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Can a classifier tell PBISC cells from real cells? (scGAN protocol, sepsis, 2026-10-17)
**Protocol (Marouf et al. 2020, scGAN):** equal numbers of cells per class → library-size normalisation + log1p → PCA, 50 components →
random forest (1,000 trees) → 5-fold cross-validated ROC AUC. AUC 0.5 = indistinguishable, 1.0 = perfectly separable. Each comparison is repeated
on 3 independent subsamples (mean ± SD). Same gene set for every source (genes shared by all four inputs).

**Inputs** (subsampled on the PC by `code/extract_rf_cells.py`):
- real sepsis = Ye 2026 held-out patients (none used for any LoRA adapter), **adults only**, 700 cells / patient; split by patient into A (reference) and B
- real COVID = Stephenson 2025 ESKD held-out patients (another study, disease and lab), 220 cells / sample
- PBISC K = 1,000 / 4,000 = Harmony-integrated generated sets, 210 cells / bulk sample (122 samples, 5 draws)

| id | comparison | what it measures |
|---|---|---|
| C0 | real A vs real A (random cell halves) | technical floor (should be ≈ 0.5) |
| C1 | real A vs real B (other patients, same study) | patient-to-patient difference |
| C2 | real A vs real COVID (other study) | study / disease / lab difference |
| C3 | real A vs real B in PBISC read-out | what the read-out format alone does |
| C4, C5 | real A vs PBISC K = 1,000 / 4,000 | generated vs real |
| C6 | PBISC B-format vs K = 1,000 | generated vs real when both are in PBISC format |
| C7 | K = 1,000 vs K = 4,000 | generator vs generator |
Each is run twice: **natural** composition (cells drawn at random) and **lineage-matched** (both classes drawn to real A's five-lineage composition, so the
classifier cannot win on cell-type proportions). Then **within lineage** (myeloid, T, NK, B).""")
code(r'''import os, glob, json, time, warnings, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors
from IPython.display import Image, display
warnings.filterwarnings("ignore"); import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
BASE = os.environ.get("PBISC_BASE", "F:/PBISC"); W = f"{BASE}/bulkonly_lora5_20261006/analysis/integration_report_20261012/rf_realism_20261017"
E, R, F = f"{W}/_extract", f"{W}/results", f"{W}/figs"; [os.makedirs(p, exist_ok=True) for p in (R, F)]
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 13, "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 200,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42})
N_TREES, N_PC, FOLDS, REPS = 1000, 50, 5, 3''')
md("## Load and align genes")
code(r'''def lin5(x):
    x = str(x)
    if x in ("CD4+ T cells", "CD8+ T cells", "NKT cells", "T cell", "MAIT", "NKT", "T g/d") or x.startswith("T "): return "T"
    if x in ("NK cells", "NK cell", "NK") or x.startswith("NK "): return "NK"
    if x in ("B cells", "Plasma cells", "B cell") or x.startswith("B_"): return "B"
    if x in ("Monocytes", "Macrophages", "DC", "Monocyte", "CD14mono", "CD16mono", "Int.mono", "pDC", "ASDC") or x.startswith("DC"): return "Myeloid"
    return "Other"
def load(src):
    m = np.load(f"{E}/{src}_meta.npz", allow_pickle=True); mats = []
    for f in sorted(glob.glob(f"{E}/{src}_block*.npz")):
        z = np.load(f); mats.append(sp.csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=(len(z["rows"]), len(m["gene_symbol"]))))
    X = sp.vstack(mats).tocsr(); lab = m["major"] if "major" in m.files else m["celltype_l1_pred"]
    obs = pd.DataFrame({"lineage": [lin5(x) for x in lab], "label": lab.astype(str), "sample": m["sample_id"].astype(str), "group": m["group"].astype(str),
                        "patient": (m["patient_id"] if "patient_id" in m.files else m["sample_id"]).astype(str)})
    if "draw" in m.files: obs["draw"] = m["draw"].astype(str)
    return X, obs, m["gene_symbol"].astype(str)
D = {s: load(s) for s in ["sepsis_real", "covid_real", "k1000", "k4000"]}
def first_idx(sym): s = pd.Series(np.arange(len(sym)), index=sym); return s[~s.index.duplicated()]
common = sorted(set.intersection(*[set(v[2]) for v in D.values()]) - {g for g in D["k1000"][2] if g.startswith("MT-")})
D = {k: (X[:, first_idx(sym).loc[common].to_numpy()].tocsr(), o) for k, (X, o, sym) in D.items()}
print("genes shared by all inputs:", len(common))
for k, (X, o) in D.items(): print(k, X.shape, o.lineage.value_counts().to_dict())''')
code(r'''# real sepsis: adults only (the bulk-only generated cohort is adult except 1 sample); patient split stratified by group
Xs, os_ = D["sepsis_real"]; keep = os_.group.str.startswith("adult").to_numpy(); Xs, os_ = Xs[keep], os_[keep].reset_index(drop=True)
rng = np.random.default_rng(7); A_p, B_p = [], []
for g, d in os_.groupby("group"):
    p = rng.permutation(sorted(d.patient.unique())); A_p += list(p[: len(p) // 2]); B_p += list(p[len(p) // 2:])
ia, ib = np.where(os_.patient.isin(A_p))[0], np.where(os_.patient.isin(B_p))[0]
REAL_A, REAL_B = (Xs[ia], os_.iloc[ia].reset_index(drop=True)), (Xs[ib], os_.iloc[ib].reset_index(drop=True))
print("real A:", len(A_p), "patients", len(ia), "cells | real B:", len(B_p), "patients", len(ib), "cells")
print(pd.DataFrame({"A": REAL_A[1].group.value_counts(), "B": REAL_B[1].group.value_counts()}))''')
md("## PBISC-like read-out of real B (same emulation as the earlier `realT2` checks)\nEach B cell → mean normalised profile of its 15 nearest B neighbours (reference PCA) restricted to the genes the cell itself detected → rescaled to a depth drawn from generated cells of the same lineage → one Poisson draw. Labels are unchanged; only the count format becomes PBISC-like.")
code(r'''def lognorm(X):
    tot = np.asarray(X.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y
Xa = REAL_A[0]; a = ad.AnnData(Xa.copy()); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000); hv = np.where(a.var.highly_variable.to_numpy())[0]
Ya = lognorm(Xa)[:, hv].toarray(); mu, sd = Ya.mean(0), Ya.std(0) + 1e-6; pca_ref = PCA(30, random_state=0).fit((Ya - mu) / sd)
proj = lambda X: pca_ref.transform(np.clip((lognorm(X)[:, hv].toarray() - mu) / sd, -10, 10))
gdepth = {}
for K in ["k1000", "k4000"]:
    X, o = D[K]; tot = np.asarray(X.sum(1)).ravel()
    for c in o.lineage.unique(): gdepth.setdefault(c, []).append(tot[(o.lineage == c).to_numpy()])
gdepth = {c: np.concatenate(v) for c, v in gdepth.items()}
def readout(X, lin, seed=0):
    r = np.random.default_rng(seed); tot = np.asarray(X.sum(1)).ravel()
    target = np.array([r.choice(gdepth[c]) for c in lin], dtype=np.float64)
    P = proj(X); _, nb = NearestNeighbors(n_neighbors=15).fit(P).kneighbors(P)
    prof = (sp.diags(1 / np.maximum(tot, 1)) @ X).tocsr()
    Wm = sp.csr_matrix((np.full(nb.size, 1 / 15, np.float32), nb.ravel(), np.arange(0, nb.size + 1, 15)), shape=(len(nb), len(nb)))
    sup = X.copy(); sup.data[:] = 1
    Em = (Wm @ prof).multiply(sup).tocsr(); Em = (sp.diags(target / np.maximum(np.asarray(Em.sum(1)).ravel(), 1e-12)) @ Em).tocsr()
    T = Em.copy(); T.data = r.poisson(Em.data).astype(np.float32); T.eliminate_zeros(); return T
REAL_BT = (readout(REAL_B[0], REAL_B[1].lineage.to_numpy()), REAL_B[1])
for nm, (X, o) in {"real A": REAL_A, "real B": REAL_B, "real B, PBISC read-out": REAL_BT, "K1000": D["k1000"], "K4000": D["k4000"], "real COVID": D["covid_real"]}.items():
    t = np.asarray(X.sum(1)).ravel(); g = np.asarray((X > 0).sum(1)).ravel(); print(f"{nm:24s} median UMI {np.median(t):7.0f}  median genes {np.median(g):6.0f}")''')
md("## Random forest AUC")
code(r'''LIN = ["T", "NK", "B", "Myeloid", "Other"]
target_comp = REAL_A[1].lineage.value_counts(normalize=True).reindex(LIN).fillna(0)
def draw(o, n, rng, matched, exclude=None, lineage=None):
    pool = np.arange(len(o)) if exclude is None else np.setdiff1d(np.arange(len(o)), exclude)
    if lineage is not None: pool = pool[(o.lineage.to_numpy()[pool] == lineage)]
    if not matched: return rng.choice(pool, min(n, len(pool)), replace=False)
    out = []
    for c in LIN:
        k = int(round(n * target_comp[c])); cand = pool[o.lineage.to_numpy()[pool] == c]
        if k: out.append(rng.choice(cand, min(k, len(cand)), replace=False))
    return np.concatenate(out)
def rf_auc(X1, X2, seed):
    X = sp.vstack([X1, X2]).tocsr(); y = np.r_[np.zeros(X1.shape[0]), np.ones(X2.shape[0])]
    Z = PCA(N_PC, random_state=seed).fit_transform(lognorm(X).toarray())
    aucs = []
    for tr, te in StratifiedKFold(FOLDS, shuffle=True, random_state=seed).split(Z, y):
        m = RandomForestClassifier(N_TREES, n_jobs=2, random_state=seed).fit(Z[tr], y[tr]); aucs.append(roc_auc_score(y[te], m.predict_proba(Z[te])[:, 1]))
    return float(np.mean(aucs))
COMP = {"C0 real A vs real A": ("A", "A"), "C1 real A vs real B": ("A", "B"), "C2 real A vs real COVID": ("A", "covid"),
        "C3 real A vs real B (PBISC read-out)": ("A", "BT"), "C4 real A vs K = 1,000": ("A", "k1000"), "C5 real A vs K = 4,000": ("A", "k4000"),
        "C6 real B (PBISC read-out) vs K = 1,000": ("BT", "k1000"), "C7 K = 1,000 vs K = 4,000": ("k1000", "k4000")}
SRC = {"A": REAL_A, "B": REAL_B, "BT": REAL_BT, "covid": D["covid_real"], "k1000": D["k1000"], "k4000": D["k4000"]}
def run(n, matched, reps, lineage=None):
    rows = []
    for name, (s1, s2) in COMP.items():
        for rep in range(reps):
            rng = np.random.default_rng(100 + rep); (X1, o1), (X2, o2) = SRC[s1], SRC[s2]
            try:
                i1 = draw(o1, n, rng, matched, lineage=lineage); i2 = draw(o2, n, rng, matched, exclude=i1 if s1 == s2 else None, lineage=lineage)
            except ValueError: continue
            if lineage: m = min(len(i1), len(i2)); i1, i2 = i1[:m], i2[:m]
            rows.append({"comparison": name, "mode": "matched" if matched else ("natural" if lineage is None else f"within {lineage}"), "rep": rep,
                         "n_per_class": int(min(len(i1), len(i2))), "auc": rf_auc(X1[i1], X2[i2], seed=rep)})
    return pd.DataFrame(rows)
t0 = time.time(); res = pd.concat([run(3000, False, REPS), run(3000, True, REPS)]); print(f"{time.time() - t0:.0f}s")
res.to_csv(f"{R}/rf_auc_whole.csv", index=False)
summ = res.groupby(["comparison", "mode"]).auc.agg(["mean", "std", "count"]).unstack("mode").round(3); summ''')
code(r'''t0 = time.time(); lin_res = pd.concat([run(1000, False, 2, lineage=L) for L in ["Myeloid", "T", "NK", "B"]]); print(f"{time.time() - t0:.0f}s")
lin_res.to_csv(f"{R}/rf_auc_within_lineage.csv", index=False)
lin_res.groupby(["comparison", "mode"]).auc.mean().unstack("mode").round(3)''')
md("## Figure")
code(r'''order = list(COMP); ylab = [c.split(" ", 1)[1] for c in order]
fig, axs = plt.subplots(1, 2, figsize=(17, 6.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]; y = np.arange(len(order))[::-1]
for dy, mode, col in [(0.18, "natural", "#9AA3B2"), (-0.18, "matched", "#2A78D6")]:
    d = res[res["mode"] == mode].groupby("comparison").auc.agg(["mean", "std"]).reindex(order)
    ax.barh(y + dy, d["mean"] - 0.5, left=0.5, height=0.34, color=col, label={"natural": "natural composition", "matched": "lineage-matched"}[mode], xerr=d["std"], error_kw={"lw": 1})
    for yi, v in zip(y + dy, d["mean"]): ax.text(v + 0.01, yi, f"{v:.2f}", va="center", fontsize=11)
ax.set_yticks(y, ylab, fontsize=12); ax.set_xlim(0.45, 1.05); ax.axvline(0.5, color="#1F2430", lw=1); ax.axvline(0.65, color="#C4373A", lw=1, ls="--")
ax.text(0.655, len(order) - 0.45, "scGAN vs real: 0.65", color="#C4373A", fontsize=11); ax.set_xlabel("random forest AUC (5-fold CV, PCA 50)"); ax.legend(loc="lower right", fontsize=11)
ax.set_title("Whole samples", loc="left")
ax = axs[1]; L = ["Myeloid", "T", "NK", "B"]; pick = [c for c in order if c[:2] in ("C1", "C3", "C4", "C5", "C6")]
cols = {"C1": "#4A4F5C", "C3": "#9AA3B2", "C4": "#2A78D6", "C5": "#EB6834", "C6": "#1BAF7A"}; w = 0.16
for j, c in enumerate(pick):
    v = [lin_res[(lin_res.comparison == c) & (lin_res["mode"] == f"within {l}")].auc.mean() for l in L]
    ax.bar(np.arange(len(L)) + (j - 2) * w, np.array(v) - 0.5, bottom=0.5, width=w, color=cols[c[:2]], label=c.split(" ", 1)[1])
ax.set_xticks(range(len(L)), L); ax.set_ylim(0.45, 1.02); ax.axhline(0.5, color="#1F2430", lw=1); ax.set_ylabel("random forest AUC"); ax.legend(fontsize=10, loc="lower left", ncol=1)
ax.set_title("Within lineage (1,000 cells per class)", loc="left")
plt.tight_layout(); fig.savefig(f"{F}/rf_auc_real_vs_generated.png"); fig.savefig(f"{F}/rf_auc_real_vs_generated.pdf"); plt.close(fig)
display(Image(f"{F}/rf_auc_real_vs_generated.png", width=1100))''')
md("## Save the cells used (one .h5ad, raw counts on the shared genes)")
code(r'''parts = []
for nm, (X, o) in {"real_A": REAL_A, "real_B": REAL_B, "real_B_pbisc_readout": REAL_BT, "real_covid": D["covid_real"], "k1000": D["k1000"], "k4000": D["k4000"]}.items():
    a = ad.AnnData(X.tocsr(), obs=o.assign(source=nm).astype(str)); a.obs_names = [f"{nm}_{i}" for i in range(a.n_obs)]; a.var_names = common; parts.append(a)
A = ad.concat(parts); A.obs = A.obs.astype("category"); A.write_h5ad(f"{W}/rf_realism_cells.h5ad", compression="gzip"); print(A)
json.dump({"n_trees": N_TREES, "n_pcs": N_PC, "folds": FOLDS, "reps": REPS, "n_per_class_whole": 3000, "n_per_class_lineage": 1000, "genes": len(common),
           "real_A_patients": A_p, "real_B_patients": B_p}, open(f"{R}/params.json", "w"), indent=1)''')
nb["cells"] = C; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/rf_realism_20261017"; os.makedirs(out, exist_ok=True)
nbf.write(nb, f"{out}/rf_realism.ipynb")
