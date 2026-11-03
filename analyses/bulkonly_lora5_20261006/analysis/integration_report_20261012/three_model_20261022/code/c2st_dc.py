"""(1) Classifier two-sample test, per lineage: real A vs each source (real B, real B depth-matched, floor, prev K4000).
Features = 2,000 real-A seurat_v3 HVG, log1p CP10k, scaled by real-A mean/sd, PCA 50 fit on real A (no batch correction: we want what is detectable).
HistGradientBoosting, 5-fold CV grouped by sample (patient identity cannot be used to separate), 4,000 cells per class per lineage, AUROC.
Also logistic regression on the 2,000 HVG to list the genes that separate floor from real A (top coefficients).
(2) Density and coverage (Naeem et al. 2020, k = 5): real A = reference manifold, source cells = samples; in the pre-Harmony PCA and in the
pairwise Harmony space (data/big_joint). Real B is the ceiling. Usage: python c2st_dc.py <coh>"""
import sys, json, numpy as np, pandas as pd, anndata as ad
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors
from common import *
coh = sys.argv[1]; D = bundle(coh, scsemiF=True); rng = np.random.default_rng(0)
genes = np.array(D["ref"].var_names); hv = L.hvg_seurat_v3([D["ref"].X.tocsr()])
Y = {s: L.lognorm(a.X.tocsr())[:, hv] for s, a in D.items()}
A = Y["ref"].toarray().astype(np.float32); mu, sd = A.mean(0), A.std(0) + 1e-6; pca = PCA(50, random_state=0).fit(np.clip((A - mu) / sd, -10, 10)); del A
Z = {s: pca.transform(np.clip((y.toarray().astype(np.float32) - mu) / sd, -10, 10)).astype(np.float32) for s, y in Y.items()}
rows, toprows = [], []
for lin in ("T", "NK", "B", "Myeloid"):
    ir = np.where(D["ref"].obs.lineage.to_numpy() == lin)[0]
    for s in ("realU", "realT1f", "realT1p", "realT1s", "floor", "prev4k", "scsemi", "scsemiF"):
        js = np.where(D[s].obs.lineage.to_numpy() == lin)[0]; n = min(4000, len(ir), len(js))
        a, b = rng.choice(ir, n, replace=False), rng.choice(js, n, replace=False)
        X = np.vstack([Z["ref"][a], Z[s][b]]); y = np.r_[np.zeros(n), np.ones(n)]
        g = np.r_["ref:" + D["ref"].obs["sample"].astype(str).to_numpy()[a], f"{s}:" + D[s].obs["sample"].astype(str).to_numpy()[b]]
        p = np.zeros(len(y))
        for tr, te in GroupKFold(5).split(X, y, g):
            p[te] = HistGradientBoostingClassifier(max_iter=200, random_state=0).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        rows.append({"cohort": coh, "lineage": lin, "source": s, "n_per_class": int(n), "c2st_auc": float(roc_auc_score(y, p))})
        if s in ("floor", "prev4k", "scsemi"):
            Xg = np.vstack([Y["ref"][a].toarray(), Y[s][b].toarray()]); lr = LogisticRegression(C=0.05, max_iter=500).fit((Xg - Xg.mean(0)) / (Xg.std(0) + 1e-6), y)
            o = np.argsort(lr.coef_[0]); d = Y[s][b].toarray().mean(0) - Y["ref"][a].toarray().mean(0)
            for k in list(o[:10]) + list(o[-10:]):
                toprows.append({"cohort": coh, "lineage": lin, "source": s, "gene": genes[hv][k], "coef": float(lr.coef_[0][k]), "mean_logdiff_src_minus_real": float(d[k]),
                                "direction": "higher in source" if lr.coef_[0][k] > 0 else "higher in real"})
    print(coh, lin, {r["source"]: round(r["c2st_auc"], 3) for r in rows if r["lineage"] == lin}, flush=True)
C = pd.DataFrame(rows); C.to_csv(f"{OUT}/{coh}_c2st.csv", index=False); pd.DataFrame(toprows).to_csv(f"{OUT}/{coh}_c2st_top_genes.csv", index=False)
# density / coverage
def dc(R, F, k=5):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(R); r = nn.kneighbors(R)[0][:, -1]          # k-th NN radius of each real point (self excluded)
    d, i = NearestNeighbors(n_neighbors=k).fit(R).kneighbors(F)                               # real neighbours of each fake point
    inside = d <= r[i]; density = inside.sum() / (k * len(F))
    # coverage: real points whose nearest fake point lies within their radius
    df, _ = NearestNeighbors(n_neighbors=1).fit(F).kneighbors(R); coverage = (df[:, 0] <= r).mean()
    return float(density), float(coverage)
dr = []
for s in ("realU", "floor", "prev4k", "scsemi"):
    J = ad.read_h5ad(f"{BJ}/{coh}_joint_{s}.h5ad"); src = J.obs.source.to_numpy(); lin = J.obs.lineage.astype(str).to_numpy()
    for space in ("X_pca", "X_pca_harmony"):
        E = J.obsm[space]
        for l in ("all", "T", "NK", "B", "Myeloid"):
            m = np.ones(len(src), bool) if l == "all" else (lin == l)
            ir, js = np.where(m & (src == "ref"))[0], np.where(m & (src != "ref"))[0]; n = min(10000, len(ir), len(js))
            de_, co_ = dc(E[rng.choice(ir, n, replace=False)], E[rng.choice(js, n, replace=False)])
            dr.append({"cohort": coh, "source": s, "space": "PCA (no correction)" if space == "X_pca" else "Harmony", "lineage": l, "n": int(n), "density": de_, "coverage": co_})
DC = pd.DataFrame(dr); DC.to_csv(f"{OUT}/{coh}_density_coverage.csv", index=False)
print(C.pivot_table(index="lineage", columns="source", values="c2st_auc").round(3).to_string())
print(DC.pivot_table(index=["space", "lineage"], columns="source", values=["density", "coverage"]).round(3).to_string())
