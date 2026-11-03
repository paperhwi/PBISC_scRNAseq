"""Re-score the seven Ye 2026 sepsis claims with panel-matched markers on both sides and refit thresholds.

Differences from paper_claim_analysis/run_fixed_claims_k_comparison.py (08:55 version), everything else identical:
  * every author marker list (Supplementary Table 6) is restricted to genes MEASURED in the generated 6k panel and present
    in the real matrix, and the SAME restricted list scores real and generated cells;
  * real cells are normalised on the same panel total (log1p CP10k over the panel genes), as generated cells are;
  * generated cells are compared with REAL CELLS SCORED BY THE SAME CLASSIFIER (classifier-positive fraction within parent),
    alongside the real exact-label fraction.
Unchanged: frozen claim definitions (fixed_paper_claim_definitions.csv), real parent = recap label prefix within
in_paper_scheme, target = exact author cluster, best_threshold() = 5-fold sample-held-out CV with sample-balanced Youden
threshold (seed 20261012), generated parent = AIFI L2 mapping YE_PARENTS.
Sensitivity variant B: real cells binomially thinned (panel UMI) to a depth drawn from generated cells of the same lineage
before scoring and threshold fitting."""
import sys, json, time, numpy as np, pandas as pd, h5py, scipy.sparse as sp
from sklearn.metrics import roc_curve, roc_auc_score
sys.path.insert(0, "/home/claude/work"); import store
t0 = time.time(); log = lambda *a: print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)
O = "/home/claude/work/report/rescore"; PCA_ = "/mnt/user-data/uploads/PBISC/paper_claim_analysis"
G1 = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw{d}.h5ad"
G4 = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_k4000_20261006/ye2026_sepsis/lora5k4_bulkonly_ye2026_sepsis_draw{d}.h5ad"
SEED = 20261012
YE_PARENTS = {
    "SEP-1": {"Memory CD4 T cell", "Naive CD4 T cell", "Treg", "Proliferating T cell"},
    "SEP-2": {"Memory CD8 T cell", "Naive CD8 T cell", "CD8aa", "Proliferating T cell"},
    "SEP-3": {"DN T cell"},
    "SEP-3-NK": {"CD56dim NK cell", "CD56bright NK cell", "Proliferating NK cell"},
    "SEP-4A": {"CD14 monocyte", "CD16 monocyte", "Intermediate monocyte"},
    "SEP-4B": {"CD14 monocyte", "CD16 monocyte", "Intermediate monocyte"},
    "SEP-5": {"CD14 monocyte"},
}


def best_threshold(score, y, groups):  # verbatim logic from run_fixed_claims_k_comparison.py
    score = np.asarray(score, dtype=np.float32); y = np.asarray(y, dtype=bool); groups = np.asarray(groups, dtype=str)
    valid = np.isfinite(score); score, y, groups = score[valid], y[valid], groups[valid]
    uniq = np.unique(groups)
    if len(uniq) < 3 or y.sum() < 20 or (~y).sum() < 20:
        return {"threshold": float(np.quantile(score[y], 0.05)) if y.any() else float("inf"), "cv_auc": np.nan,
                "cv_sensitivity": np.nan, "cv_specificity": np.nan, "n_samples": len(uniq), "cv_note": "fallback"}
    rng = np.random.default_rng(SEED)
    gi = {g: np.flatnonzero(groups == g) for g in uniq}
    def sample_for_roc(sg):
        chosen, counts = [], {}
        for g in sg:
            ig = gi[g]; pos = ig[y[ig]]; neg = ig[~y[ig]]
            if len(pos) > 250: pos = rng.choice(pos, 250, replace=False)
            if len(neg) > 1000: neg = rng.choice(neg, 1000, replace=False)
            chosen.extend(pos.tolist()); chosen.extend(neg.tolist()); counts[g] = len(pos) + len(neg)
        return np.asarray(chosen, dtype=np.int64), counts
    perm = rng.permutation(uniq); fold_of = {g: i % min(5, len(uniq)) for i, g in enumerate(perm)}
    pred = np.zeros(len(y), dtype=bool)
    for fold in range(min(5, len(uniq))):
        test = np.isin(groups, [g for g in uniq if fold_of[g] == fold])
        fit_idx, counts = sample_for_roc([g for g in uniq if fold_of[g] != fold])
        if not test.any() or y[fit_idx].sum() == 0 or (~y[fit_idx]).sum() == 0: continue
        w = np.asarray([1.0 / counts[g] for g in groups[fit_idx]])
        fpr, tpr, th = roc_curve(y[fit_idx], score[fit_idx], sample_weight=w); pred[test] = score[test] >= th[np.argmax(tpr - fpr)]
    auc_idx, _ = sample_for_roc(uniq); auc = roc_auc_score(y[auc_idx], score[auc_idx])
    fit_idx, counts = sample_for_roc(uniq); w = np.asarray([1.0 / counts[g] for g in groups[fit_idx]])
    fpr, tpr, th = roc_curve(y[fit_idx], score[fit_idx], sample_weight=w)
    return {"threshold": float(th[np.argmax(tpr - fpr)]), "cv_auc": float(auc), "cv_sensitivity": float(pred[y].mean()),
            "cv_specificity": float((~pred[~y]).mean()), "n_samples": len(uniq), "cv_note": "5-fold sample-held-out"}


def gen_symbols(path):
    with h5py.File(path) as f:
        v = f["var"]; s = v["gene_symbol"]
        sym = (s["categories"][:][s["codes"][:]] if isinstance(s, h5py.Group) else s[:]).astype(str)
        meas = v["measured_in_cohort"][:].astype(bool)
    return sym, meas


def obs_col(f, k):
    x = f["obs"][k]
    if isinstance(x, h5py.Group): return x["categories"][:].astype(str)[x["codes"][:]]
    return x[:].astype(str)


def lineage_gen(l1):
    l1 = np.asarray(l1, str); out = np.full(len(l1), "Other", object)
    for key, lab in [("T cell", "T"), ("NK", "NK"), ("Monocyte", "Myeloid"), ("DC", "Myeloid"), ("B cell", "B"), ("Plasma", "B")]:
        out[np.char.find(l1, key) >= 0] = lab
    return out


def lineage_real(m):
    return pd.Series(m).map({"CD4+ T cells": "T", "CD8+ T cells": "T", "NKT cells": "T", "NK cells": "NK", "Monocytes": "Myeloid",
                             "Macrophages": "Myeloid", "DC": "Myeloid", "B cells": "B", "Plasma cells": "B"}).fillna("Other").to_numpy()


# ---------------- panel and claims
sym1, meas1 = gen_symbols(G1.format(d=0)); sym4, meas4 = gen_symbols(G4.format(d=0))
assert (sym1 == sym4).all() and (meas1 == meas4).all()
st = store.Store("/home/claude/work/sepsis/store"); rsym = st.var.gene_symbol.astype(str).to_numpy()
rpos = pd.Series(np.arange(len(rsym)), index=rsym); rpos = rpos[~rpos.index.duplicated()]
gpos = pd.Series(np.arange(len(sym1)), index=sym1)[meas1]; gpos = gpos[~gpos.index.duplicated()]
panel = [g for g in gpos.index if g in rpos.index]
P = pd.Series(np.arange(len(panel)), index=panel)
defs = pd.read_csv(f"{PCA_}/fixed_paper_claim_definitions.csv", dtype=str).fillna("")
defs = defs[defs.study.eq("Ye et al. 2026")].set_index("id")
claims = list(defs.index)
markers = {c: [g for g in dict.fromkeys(defs.loc[c, "markers"].split(";")) if g in P.index] for c in claims}
cover = pd.DataFrame({"claim": claims, "author_markers": [len(set(defs.loc[c, "markers"].split(";")) - {""}) for c in claims],
                      "panel_markers": [len(markers[c]) for c in claims], "panel_markers_list": [";".join(markers[c]) for c in claims]})
cover.to_csv(f"{O}/panel_marker_coverage.csv", index=False)
log("panel genes", len(panel), dict(zip(claims, cover.panel_markers)))
W = sp.csr_matrix((np.concatenate([np.full(len(markers[c]), 1 / len(markers[c])) for c in claims]),
                   (np.concatenate([P.loc[markers[c]].to_numpy() for c in claims]),
                    np.concatenate([np.full(len(markers[c]), j) for j, c in enumerate(claims)]))), shape=(len(panel), len(claims)))


def score(Xp):
    tot = np.asarray(Xp.sum(1)).ravel()
    Y = (sp.diags((1e4 / np.maximum(tot, 1)).astype(np.float32)) @ Xp).tocsr(); Y.data = np.log1p(Y.data)
    return np.asarray((Y @ W).todense(), dtype=np.float32), tot


# ---------------- generated panel depth by lineage (for variant B)
gdepth = {}
for path in [G1.format(d=0), G4.format(d=0)]:
    with h5py.File(path) as f:
        n = len(f["X/indptr"]) - 1; lin = lineage_gen(obs_col(f, "aifi_l1"))
        Xg = sp.csr_matrix((f["X/data"][:].astype(np.float32), f["X/indices"][:], f["X/indptr"][:]), shape=(n, len(sym1)))[:, gpos.loc[panel].to_numpy()]
        tot = np.asarray(Xg.sum(1)).ravel()
    for c in np.unique(lin): gdepth.setdefault(c, []).append(tot[lin == c])
gdepth = {c: np.concatenate(v) for c, v in gdepth.items()}
log("generated panel depth medians", {c: int(np.median(v)) for c, v in gdepth.items()})

# ---------------- real cells: scores (A: as is, B: depth-matched)
ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                            usecols=["sample", "age_group", "group", "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme"]) for i in range(3)], ignore_index=True)
ok = ro.in_paper_scheme.astype(str).eq("True").to_numpy(); idx = np.where(ok)[0]
rlin = lineage_real(ro.celltype_major_ye2026.to_numpy())
SA = np.zeros((len(ro), len(claims)), np.float32); SB = np.zeros_like(SA); RT = np.zeros(len(ro), np.float32)
rng = np.random.default_rng(SEED); cols = rpos.loc[panel].to_numpy()
for i in range(0, len(idx), 20000):
    c = idx[i:i + 20000]; Xp = st.rows(c)[:, cols].tocsr()
    SA[c], RT[c] = score(Xp)
    tgt = np.array([rng.choice(gdepth.get(l, gdepth["T"])) for l in rlin[c]], np.float64)
    p = np.minimum(1.0, tgt / np.maximum(RT[c], 1)); Xt = Xp.copy()
    Xt.data = rng.binomial(Xt.data.astype(np.int64), np.repeat(p, np.diff(Xt.indptr))).astype(np.float32); Xt.eliminate_zeros()
    SB[c], _ = score(Xt)
    if i % 200000 == 0: log("real scored", i + len(c), "/", len(idx))
np.savez_compressed(f"{O}/real_scores.npz", SA=SA[ok], SB=SB[ok], panel_umi=RT[ok], cells=idx)
log("real median panel UMI", float(np.median(RT[ok])))

# ---------------- thresholds (same procedure, refit per variant)
lab = ro.celltype_ye2026.astype(str).to_numpy(); smp = ro["sample"].astype(str).to_numpy()
cv_rows, thr = [], {}
old = pd.read_csv(f"{PCA_}/fixed_claim_sepsis_thresholds_cv.csv").set_index("claim_id")
for j, c in enumerate(claims):
    par = ok & pd.Series(lab).str.startswith(defs.loc[c, "lineage_prefix"]).to_numpy(); y = lab == defs.loc[c, "target"]
    for v, S in [("A_panel", SA), ("B_panel_depth", SB)]:
        r = best_threshold(S[par, j], y[par], smp[par]); thr[(v, c)] = r["threshold"]
        cv_rows.append({"claim_id": c, "variant": v, **r, "panel_markers": len(markers[c]), "author_markers": int(cover.author_markers[j]),
                        "previous_full_panel_auc": float(old.loc[c, "cv_auc"])})
cv = pd.DataFrame(cv_rows); cv.to_csv(f"{O}/thresholds_cv_panel_matched.csv", index=False)
log("CV\n" + cv[["claim_id", "variant", "cv_auc", "cv_sensitivity", "cv_specificity", "previous_full_panel_auc", "panel_markers"]].round(3).to_string())

# ---------------- real per-sample: exact-label fraction and classifier-positive fraction within parent
grp = (ro.age_group.str[:5].replace({"paedi": "ped"}) + "_" + ro.group.replace({"healthy control": "HC", "ICU control": "ICU_ctrl"})).to_numpy()
real_rows = []
for j, c in enumerate(claims):
    par = ok & pd.Series(lab).str.startswith(defs.loc[c, "lineage_prefix"]).to_numpy(); y = lab == defs.loc[c, "target"]
    d = pd.DataFrame({"sample": smp[par], "group": grp[par], "y": y[par], "A": SA[par, j] >= thr[("A_panel", c)],
                      "B": SB[par, j] >= thr[("B_panel_depth", c)], "s": SA[par, j]})
    a = d.groupby(["sample", "group"]).agg(n_parent=("y", "size"), label_frac=("y", "mean"), posA_frac=("A", "mean"), posB_frac=("B", "mean"), score_mean=("s", "mean")).reset_index()
    a.insert(0, "claim_id", c); real_rows.append(a)
real = pd.concat(real_rows); real.to_csv(f"{O}/real_by_sample.csv", index=False)

# ---------------- generated draws
gen_rows = []
for K, tmpl in [(1000, G1), (4000, G4)]:
    for dr in range(5):
        with h5py.File(tmpl.format(d=dr)) as f:
            n = len(f["X/indptr"]) - 1; l2 = obs_col(f, "aifi_l2"); s_id = obs_col(f, "sample_id"); g = obs_col(f, "group")
            ip = f["X/indptr"][:]; Sg = np.zeros((n, len(claims)), np.float32); gt = np.zeros(n, np.float32)
            for a0 in range(0, n, 50000):
                b0 = min(n, a0 + 50000); lo, hi = ip[a0], ip[b0]
                Xg = sp.csr_matrix((f["X/data"][lo:hi].astype(np.float32), f["X/indices"][lo:hi], ip[a0:b0 + 1] - lo), shape=(b0 - a0, len(sym1)))
                Sg[a0:b0], gt[a0:b0] = score(Xg[:, gpos.loc[panel].to_numpy()].tocsr())
        for j, c in enumerate(claims):
            par = np.isin(l2, list(YE_PARENTS[c]))
            d = pd.DataFrame({"sample": s_id[par], "group": g[par], "A": Sg[par, j] >= thr[("A_panel", c)], "B": Sg[par, j] >= thr[("B_panel_depth", c)], "s": Sg[par, j]})
            a = d.groupby(["sample", "group"]).agg(n_parent=("A", "size"), posA_frac=("A", "mean"), posB_frac=("B", "mean"), score_mean=("s", "mean")).reset_index()
            a.insert(0, "claim_id", c); a.insert(1, "K", K); a.insert(2, "draw", dr); gen_rows.append(a)
        log(f"K={K} draw={dr} scored, {n} cells, median panel UMI {np.median(gt):.0f}")
gen = pd.concat(gen_rows); gen.to_csv(f"{O}/generated_by_sample_draw.csv", index=False)
json.dump({"panel_genes": len(panel), "seed": SEED}, open(f"{O}/run_meta.json", "w"))
log("done")
