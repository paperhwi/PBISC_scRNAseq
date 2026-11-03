"""B. Standard deconvolution of the same sepsis bulk (GSE279448) with the same single-cell reference, plus a no-information null.

Reference: Ye 2026 sepsis scRNA (in_paper_scheme cells), summed counts per (sample, fine type) from p1_pseudobulk.py.
Signature: mean UMI per cell of each type (cell-size aware, so coefficients are cell fractions), on marker genes.
Methods
  NNLS     : non-negative least squares on linear counts
  WNNLS    : NNLS with each gene scaled by 1 / its mean signature value (relative error; dampens high-expression genes)
  SVR      : CIBERSORT-style linear nu-SVR (nu 0.25/0.5/0.75, lowest RMSE), negatives set to 0
  WNNLS+cal: WNNLS after a per-gene bulk -> pseudobulk correction learned from the 5 paired donors that trained each
             LoRA draw (five versions, matching PBISC's LoRA-5 information budget)
Evaluation
  paired   : bulk samples that have a single-cell partner; 2-fold by sample, signature and calibration donors from the
             other fold only; per-sample error vs the sample's own single-cell composition
  bulkonly : the 122 bulk-only samples PBISC generated from; group means vs the real single-cell group composition
Null: pooled single-cell composition of all adult samples, identical for every sample (carries no bulk information)."""
import glob, json, numpy as np, pandas as pd
from scipy.optimize import nnls
from sklearn.svm import NuSVR
O = "/home/claude/work/prop"; rng = np.random.default_rng(20261013)

# ---------------- reference pseudobulk
idx = pd.read_csv(f"{O}/sepsis_pb_index.csv"); PB = np.load(f"{O}/sepsis_pb_counts.npy", mmap_mode="r")
var = pd.read_csv(f"{O}/sepsis_pb_var.csv")
obs = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                             usecols=["sample", "age_group", "group"]) for i in range(3)], ignore_index=True).drop_duplicates("sample")
obs["grp"] = obs.age_group.str[:5].replace({"paedi": "ped"}) + "_" + obs.group.replace({"healthy control": "HC", "ICU control": "ICU_ctrl"})
grp_of = obs.set_index("sample").grp
idx["grp"] = idx["sample"].map(grp_of)


def lin6(x):
    if x == "Plasma cells": return "Plasma"
    if x == "B cells": return "B"
    if x in ("CD4+ T cells", "CD8+ T cells", "NKT cells"): return "T"
    if x == "NK cells": return "NK"
    if x in ("Monocytes", "Macrophages", "DC"): return "Myeloid"
    return "Other"


idx["lin6"] = idx.major.map(lin6)
LEVELS = ["lin6", "major", "fine"]

# ---------------- bulk
files = sorted(glob.glob(f"{O}/bulk/*.htseq.results.gz"))
bulk = {}
for f in files:
    s = f.split("/")[-1].split("_", 1)[1].replace(".htseq.results.gz", "")
    t = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]; bulk[s] = t[~t.index.str.startswith("__")]
bulk = pd.DataFrame(bulk); bulk.index = bulk.index.str.split(".").str[0]
gid = var.gene_id.astype(str).str.split(".").str[0].to_numpy()
common = np.array([g for g in pd.unique(gid) if g in bulk.index])
gpos = pd.Series(np.arange(len(gid)), index=gid); gpos = gpos[~gpos.index.duplicated()]
gi = gpos.loc[common].to_numpy()
Bk = bulk.loc[common].astype(float)
bo = pd.read_csv("/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/samples.csv")
bulkonly = [s for s in bo.sample_id.astype(str) if s in Bk.columns]
sc_samples = set(idx["sample"])
paired = [s for s in Bk.columns if s in sc_samples and grp_of.get(s, "").startswith("adult")]
donors = json.load(open(f"{O}/sepsis_lora_training_donors.json"))
print("genes", len(common), "bulk", Bk.shape[1], "bulk-only", len(bulkonly), "paired adult", len(paired), flush=True)


def type_sums(samples, lv):
    m = idx["sample"].isin(samples).to_numpy()
    sub = idx[m]; rows = np.where(m)[0]
    labs = sub[lv].to_numpy(); cats = sorted(pd.unique(idx[lv]))
    S = np.zeros((len(cats), len(gi))); n = np.zeros(len(cats))
    for j, c in enumerate(cats):
        r = rows[labs == c]
        if len(r): S[j] = np.asarray(PB[r][:, gi]).sum(0); n[j] = sub.n_cells.to_numpy()[labs == c].sum()
    return cats, S, n


def signature(samples, lv, n_mark=50):
    cats, S, n = type_sums(samples, lv)
    keep = n >= 30
    cats = [c for c, k in zip(cats, keep) if k]; S = S[keep]; n = n[keep]
    sig = S / n[:, None]                                     # mean UMI per cell
    cpm = S / S.sum(1, keepdims=True) * 1e6                  # type profile
    lc = np.log2(cpm + 1)
    mk = set()
    for j in range(len(cats)):
        other = np.delete(lc, j, 0).max(0)
        score = lc[j] - other; ok = cpm[j] > 20
        mk |= set(np.where(ok)[0][np.argsort(-score[ok])[:n_mark]])
    mk = np.array(sorted(mk))
    return cats, sig[:, mk].T, mk                            # genes x types


def solve(sig, b, method):
    if method == "NNLS":
        p, _ = nnls(sig, b)
    elif method.startswith("WNNLS"):
        w = 1 / np.maximum(sig.mean(1), 1e-3); p, _ = nnls(sig * w[:, None], b * w)
    else:                                                    # CIBERSORT-style
        X = (sig - sig.mean()) / sig.std(); y = (b - b.mean()) / b.std(); best = None
        for nu in (0.25, 0.5, 0.75):
            m = NuSVR(nu=nu, C=1.0, kernel="linear").fit(X, y); c = np.maximum(m.coef_.ravel(), 0)
            rmse = np.sqrt(np.mean((X @ c - y) ** 2)) if c.sum() > 0 else np.inf
            if best is None or rmse < best[0]: best = (rmse, c)
        p = best[1]
    return p / p.sum() if p.sum() > 0 else np.full(len(p), 1 / len(p))


def calib_factor(cal_samples, mk):
    """per-gene log ratio pseudobulk / bulk over calibration donors (both as CPM over the common genes)."""
    r = []
    for s in cal_samples:
        if s not in Bk.columns or s not in sc_samples: continue
        pbs = np.asarray(PB[np.where(idx["sample"].to_numpy() == s)[0]][:, gi]).sum(0)
        r.append(np.log((pbs / pbs.sum() * 1e6 + 1) / (Bk[s].to_numpy() / Bk[s].sum() * 1e6 + 1)))
    return np.exp(np.mean(r, 0))[mk] if r else np.ones(len(mk))


def run(samples_eval, ref_samples, lv, cal_sets):
    cats, sig, mk = signature(ref_samples, lv); out = []
    facs = {d: calib_factor(c, mk) for d, c in cal_sets.items()}
    for s in samples_eval:
        b = Bk[s].to_numpy()[mk]; b = b / b.sum() * sig.sum(0).mean() * 1  # scale is irrelevant for the fractions
        for m in ["NNLS", "WNNLS", "SVR"]:
            out.append(pd.DataFrame({"sample": s, "method": m, "draw": "-", "label": cats, "frac": solve(sig, b, m)}))
        for d, f in facs.items():
            out.append(pd.DataFrame({"sample": s, "method": "WNNLS+cal", "draw": d, "label": cats, "frac": solve(sig, b * f, "WNNLS")}))
    return pd.concat(out)


adult_sc = [s for s in sc_samples if grp_of.get(s, "").startswith("adult")]
res = []
# ---- paired, 2-fold by sample
perm = rng.permutation(paired); folds = [perm[: len(perm) // 2], perm[len(perm) // 2:]]
for k, ev in enumerate(folds):
    ref = [s for s in sc_samples if s not in set(ev)]
    cal = {"fold_cal": list(rng.choice([s for s in paired if s not in set(ev)], 5, replace=False))}
    for lv in LEVELS:
        r = run(ev, ref, lv, cal); r["set"] = "paired"; r["level"] = lv; res.append(r)
    print("paired fold", k, "done", flush=True)
# ---- bulk-only: full reference; calibration = each LoRA draw's five training donors
for lv in LEVELS:
    r = run(bulkonly, list(sc_samples), lv, donors); r["set"] = "bulkonly"; r["level"] = lv; res.append(r)
    print("bulk-only", lv, "done", flush=True)
D = pd.concat(res, ignore_index=True)
D["grp"] = D["sample"].map(lambda s: grp_of.get(s, None)).fillna(D["sample"].map(bo.set_index("sample_id").group))
D.to_csv(f"{O}/sepsis_deconv_props.csv.gz", index=False)

# ---------------- truth: per-sample single-cell composition (all in-scheme cells)
T = []
for lv in LEVELS:
    t = idx.groupby(["sample", lv]).n_cells.sum().unstack(fill_value=0); t = t.div(t.sum(1), axis=0)
    t = t.stack().rename("frac").reset_index().rename(columns={lv: "label"}); t["level"] = lv; t["grp"] = t["sample"].map(grp_of); T.append(t)
pd.concat(T).to_csv(f"{O}/sepsis_sc_truth_props.csv.gz", index=False)
json.dump({"genes_common": int(len(common)), "n_bulk": int(Bk.shape[1]), "bulkonly": bulkonly, "paired_adult": paired,
           "folds": [list(map(str, f)) for f in folds]}, open(f"{O}/sepsis_deconv_meta.json", "w"), indent=1)
print("done")
