"""Step 2: same-instrument classifiers trained on thinned real cells, applied to real (cross-fitted) and generated cells.

Outputs (results/):
  classifier_performance.csv   held-out (sample-grouped) AUC, TPR, FPR per stage/target
  calls_real.csv               per real sample x target: parent_n, state_n (instrument), cross-fitted
  calls_gen.csv                per generated draw x sample x target: parent_n, state_n
"""
import json, os
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = "F:/PBISC/fig6_bulkonly_reproduction_20261020"
GEN = "F:/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw{}.h5ad"
os.makedirs(f"{ROOT}/results", exist_ok=True)
TARGETS = {"F1": ("CD4T", "T_CD4_c02-GPR183-NR4A2"), "F2": ("CD8T", "T_CD8_c05-GNLY-CCL4hi"),
           "F3": ("NK", "NK_c03-PRF1-CCL4hi"), "F4": ("NKT", "NKT_c05-IFNG-CCL4hi"), "C1": ("CD4T", "T_CD4_c01-LEF1")}
LINS = ["CD4T", "CD8T", "NK", "NKT", "other"]

R = ad.read_h5ad(f"{ROOT}/data/real_tnk_panel.h5ad")
panel_ids = R.var.gene_id.to_numpy()
groups = R.obs["sample"].to_numpy()
lin = R.obs.lin.to_numpy(); lab = R.obs.lab.to_numpy()


def lognorm(X):
    X = sp.csr_matrix(X, dtype=np.float32)
    Y = (sp.diags(1e4 / np.maximum(np.asarray(X.sum(1)).ravel(), 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y


class Feat:
    def __init__(self, Xc, n=2000):
        a = ad.AnnData(lognorm(Xc)); sc.pp.highly_variable_genes(a, n_top_genes=n, flavor="seurat")
        self.hv = np.where(a.var.highly_variable.to_numpy())[0]
        Y = a.X[:, self.hv].toarray(); self.mu = Y.mean(0); self.sd = Y.std(0) + 1e-3

    def __call__(self, Xc):
        return np.clip((lognorm(Xc)[:, self.hv].toarray() - self.mu) / self.sd, -10, 10).astype(np.float32)


def stage1_model(): return LogisticRegression(C=0.05, max_iter=300)
def stage2_model(): return HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, class_weight="balanced", random_state=0)


cv = GroupKFold(5)
perf = []
# ---------------- stage 1: lineage
F1 = Feat(R.X)
XR = F1(R.X)
lin_cf = np.empty(len(lin), dtype=object)
for tr, te in cv.split(XR, lin, groups):
    lin_cf[te] = stage1_model().fit(XR[tr], lin[tr]).predict(XR[te])
acc = (lin_cf == lin).mean()
for L in LINS:
    t, p = lin == L, lin_cf == L
    perf.append(dict(stage="lineage", target=L, auc=np.nan, recall=(p & t).sum() / max(t.sum(), 1),
                     precision=(p & t).sum() / max(p.sum(), 1), overall_acc=acc))
m1 = stage1_model().fit(XR, lin)
del XR
print("stage1 acc", round(acc, 3), flush=True)

# ---------------- stage 2: state within parent (train on true parent)
m2, f2, thr = {}, {}, {}
state_cf = {}
for k, (par, st) in TARGETS.items():
    ix = np.where(lin == par)[0]; y = (lab[ix] == st).astype(int)
    f = Feat(R.X[ix]); Xp = f(R.X[ix]); p = np.zeros(len(ix))
    for tr, te in cv.split(Xp, y, groups[ix]):
        p[te] = stage2_model().fit(Xp[tr], y[tr]).predict_proba(Xp[te])[:, 1]
    fpr, tpr, th = roc_curve(y, p); j = np.argmax(tpr - fpr); thr[k] = th[j]
    pred = p >= thr[k]
    perf.append(dict(stage="state", target=k, auc=roc_auc_score(y, p), recall=pred[y == 1].mean(), precision=y[pred].mean(),
                     TPR=pred[y == 1].mean(), FPR=pred[y == 0].mean(), threshold=thr[k], n_train=len(ix), prevalence=y.mean()))
    m2[k] = stage2_model().fit(Xp, y); f2[k] = f
    # cross-fitted calls for ALL real cells with cross-fitted parent == par (same route as generated)
    jx = np.where(lin_cf == par)[0]; Xq = f(R.X[jx]); q = np.zeros(len(jx))
    for tr, te in cv.split(Xp, y, groups[ix]):
        test_samples = set(groups[ix][te]); mm = np.isin(groups[jx], list(test_samples))
        if mm.any():
            q[mm] = stage2_model().fit(Xp[tr], y[tr]).predict_proba(Xq[mm])[:, 1]
    state_cf[k] = (jx, q >= thr[k])
    print(k, "AUC", round(perf[-1]["auc"], 3), flush=True)
pd.DataFrame(perf).to_csv(f"{ROOT}/results/classifier_performance.csv", index=False)

# per-sample real instrument counts
rows = []
for k in TARGETS:
    jx, pos = state_cf[k]
    d = pd.DataFrame({"sample": groups[jx], "pos": pos}).groupby("sample").pos.agg(["count", "sum"])
    for s in np.unique(groups):
        n, kk = (int(d.loc[s, "count"]), int(d.loc[s, "sum"])) if s in d.index else (0, 0)
        rows.append(dict(sample=s, target=k, parent_n=n, state_n=kk))
pd.DataFrame(rows).to_csv(f"{ROOT}/results/calls_real.csv", index=False)

# ---------------- generated cells
rows = []
for j in range(5):
    g = ad.read_h5ad(GEN.format(j))
    gid = pd.Series(np.arange(g.n_vars), index=g.var.gene_id.to_numpy()); gid = gid[~gid.index.duplicated()]
    keep = np.where(g.obs.celltype_l1_pred.isin(["T cell", "NK cell"]).to_numpy())[0]
    Xg = sp.csr_matrix(g.X[keep][:, gid.loc[panel_ids].to_numpy()], dtype=np.float32)
    smp = g.obs.sample_id.astype(str).to_numpy()[keep]
    lg = m1.predict(F1(Xg))
    for k, (par, st) in TARGETS.items():
        jx = np.where(lg == par)[0]
        pos = (m2[k].predict_proba(f2[k](Xg[jx]))[:, 1] >= thr[k]) if len(jx) else np.zeros(0, bool)
        d = pd.DataFrame({"sample": smp[jx], "pos": pos}).groupby("sample").pos.agg(["count", "sum"])
        for s in np.unique(g.obs.sample_id.astype(str)):
            n, kk = (int(d.loc[s, "count"]), int(d.loc[s, "sum"])) if s in d.index else (0, 0)
            rows.append(dict(draw=j, sample=s, target=k, parent_n=n, state_n=kk))
    print("draw", j, "done", flush=True)
    del g, Xg
pd.DataFrame(rows).to_csv(f"{ROOT}/results/calls_gen.csv", index=False)
json.dump({k: float(v) for k, v in thr.items()}, open(f"{ROOT}/results/thresholds.json", "w"), indent=1)
