"""Equal-information deconvolution baselines for the paired PBISC ensemble (58 adult sepsis samples, six types).

PBISC's five-model ensemble uses the five LoRA draws' training donors (25 donor sets in total, 5 per draw).
Matching baselines, all with the reference built from scRNA samples outside the 58 test samples:
  WNNLS+cal(draw j) : bulk corrected gene by gene with draw j's 5 training donors (single-model budget)
  WNNLS+cal ensemble: mean of the five draw-specific estimates (ensemble budget, same structure as PBISC)
  WNNLS+cal25       : one correction from the union of the 25 training donors
  BayesPrism+cal25  : BayesPrism (major types, fine states) on the bulk corrected with the 25 donors
Truth: recap scRNA composition mapped to the paired benchmark's six types (benchmark6), as in paired_accuracy.py."""
import json, glob, numpy as np, pandas as pd, subprocess, shutil
from scipy.optimize import nnls
O = "/home/claude/work/prop"; E = f"{O}/eq"
PV = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014"
idx = pd.read_csv(f"{O}/sepsis_pb_index.csv"); PB = np.load(f"{O}/sepsis_pb_counts.npy", mmap_mode="r"); var = pd.read_csv(f"{O}/sepsis_pb_var.csv")
donors = json.load(open(f"{O}/sepsis_lora_training_donors.json")); D25 = sorted(set(sum(donors.values(), [])))
TYPES = ["B cell", "DC", "Monocyte", "NK cell", "Other", "T cell"]
def b6(x):
    return {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell",
            "NK cells": "NK cell", "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC"}.get(x, "Other")
# bulk
bulk = {}
for f in sorted(glob.glob(f"{O}/bulk/*.htseq.results.gz")):
    s = f.split("/")[-1].split("_", 1)[1].replace(".htseq.results.gz", ""); t = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]; bulk[s] = t[~t.index.str.startswith("__")]
bulk = pd.DataFrame(bulk); bulk.index = bulk.index.str.split(".").str[0]
gid = var.gene_id.astype(str).str.split(".").str[0].to_numpy(); sym = var.gene_symbol.astype(str).to_numpy()
gpos = pd.Series(np.arange(len(gid)), index=gid); gpos = gpos[~gpos.index.duplicated()]
common = np.array([g for g in gpos.index if g in bulk.index]); gi = gpos.loc[common].to_numpy(); Bk = bulk.loc[common].astype(float)
# test samples = the 58 used in paired_accuracy.py (adult paired, in deconvolution set, with PBISC predictions)
meta = json.load(open(f"{O}/sepsis_deconv_meta.json"))
truth = pd.read_csv(f"{O}/sepsis_sc_truth_props.csv.gz"); truth = truth[(truth.level == "major") & truth.grp.str.startswith("adult")]
truth["type"] = truth.label.map(b6); truth = truth.groupby(["sample", "type"]).frac.sum().unstack(fill_value=0).reindex(columns=TYPES, fill_value=0)
pred = pd.read_csv(f"{PV}/results/paired_predictions.csv.gz")
print(pred.columns.tolist()[:12], flush=True)
ER = pd.read_csv(f"{PV}/results/paired_sample_errors.csv.gz")
ids = sorted(ER[(ER.scope == "common60_recap_truth") & (ER.method == "PBISC_k4000_eligible_ensemble")].sample_id.unique())
print("test samples", len(ids), flush=True)
ref = [s for s in idx["sample"].unique() if s not in set(ids)]          # reference scRNA: everything outside the 58


def markers_sig(samples, lv="major", n=50):
    m = idx["sample"].isin(samples).to_numpy(); sub = idx[m]; rows = np.where(m)[0]
    cats = [c for c in sorted(sub[lv].unique()) if sub.loc[sub[lv] == c, "n_cells"].sum() >= 30]
    S = np.vstack([np.asarray(PB[rows[sub[lv].to_numpy() == c]][:, gi]).sum(0) for c in cats]); n_ = np.array([sub.loc[sub[lv] == c, "n_cells"].sum() for c in cats])
    lc = np.log2(S / S.sum(1, keepdims=True) * 1e6 + 1); mk = set()
    for j in range(len(cats)):
        sc = lc[j] - np.delete(lc, j, 0).max(0); ok = lc[j] > np.log2(21); mk |= set(np.where(ok)[0][np.argsort(-sc[ok])[:n]])
    mk = np.array(sorted(mk)); return cats, (S / n_[:, None])[:, mk].T, mk


def calib(cal, mk):
    r = []
    for s in cal:
        if s not in Bk.columns or s not in set(idx["sample"]): continue
        pbs = np.asarray(PB[np.where(idx["sample"].to_numpy() == s)[0]][:, gi]).sum(0)
        r.append(np.log((pbs / pbs.sum() * 1e6 + 1) / (Bk[s].to_numpy() / Bk[s].sum() * 1e6 + 1)))
    return np.exp(np.mean(r, 0))[mk], len(r)


def wnnls(sig, b):
    w = 1 / np.maximum(sig.mean(1), 1e-3); p, _ = nnls(sig * w[:, None], b * w); return p / p.sum()


cats, sig, mk = markers_sig(ref)
rows, used = [], {}
for s in ids:
    b = Bk[s].to_numpy()[mk]
    per = []
    for d, dl in donors.items():
        if s in dl: continue                                      # same eligibility rule as the PBISC ensemble
        f, n = calib(dl, mk); p = wnnls(sig, b * f); per.append(p)
        rows.append(dict(sample=s, method=f"WNNLS+cal_{d}", **dict(zip(cats, p))))
    rows.append(dict(sample=s, method="WNNLS+cal ensemble", **dict(zip(cats, np.mean(per, 0)))))
    c25 = [x for x in D25 if x != s]; f, n = calib(c25, mk); used[s] = n
    rows.append(dict(sample=s, method="WNNLS+cal25", **dict(zip(cats, wnnls(sig, b * f)))))
P = pd.DataFrame(rows).fillna(0); P.to_csv(f"{E}/equal_info_major_props.csv", index=False)
print("calibration donors usable per sample (of 25):", sorted(set(used.values())), flush=True)
# BayesPrism + cal25 inputs (one run, all 58 mixtures; the sample itself is never among the 25 calibration donors' rows
# because test samples are not LoRA training donors except where noted — those get their own factor above, here we use 25 minus overlaps)
over = [s for s in ids if s in D25]; print("test samples that are LoRA donors:", over, flush=True)
f25, n25 = calib([x for x in D25 if x not in set(ids)], np.arange(len(gi)))
m = idx["sample"].isin(ref).to_numpy() & (idx.n_cells >= 1).to_numpy(); sub = idx[m].reset_index(drop=True)
# genes: union of fine and major markers on the reference (same rule as b2_bayesprism_prep.py)
_, _, mkM = markers_sig(ref, "major"); _, _, mkF = markers_sig(ref, "fine"); g = np.array(sorted(set(mkM) | set(mkF)))
bad = pd.Series(sym[gi[g]]).str.match(r"^(RPL|RPS|MT-|MRPL|MRPS)").to_numpy() | np.isin(sym[gi[g]], ["XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "MALAT1", "NEAT1"])
g = g[~bad]; names = pd.Series(sym[gi[g]]); keep = ~names.duplicated().to_numpy(); g = g[keep]
R = pd.DataFrame(np.rint(np.asarray(PB[np.where(m)[0]][:, gi[g]])).astype(np.int64), columns=sym[gi[g]]); R.index = [f"r{i}" for i in range(len(R))]
R.to_csv(f"{O}/bp/eq25cal_ref.csv.gz"); sub[["fine", "major", "n_cells"]].assign(row=R.index).to_csv(f"{O}/bp/eq25cal_ref_labels.csv", index=False)
mix = Bk.loc[common[g], ids].T * f25[g][None, :]; mix.columns = sym[gi[g]]; np.rint(mix).astype(int).to_csv(f"{O}/bp/eq25cal_mix.csv.gz")
tot = np.asarray(PB[np.where(m)[0]]).sum(1); pd.DataFrame({"major": sub.major, "fine": sub.fine, "umi": tot, "n": sub.n_cells}).to_csv(f"{O}/bp/eq25cal_size.csv", index=False)
json.dump({"test_samples": ids, "calibration_donors_25": D25, "bayesprism_cal_donors": [x for x in D25 if x not in set(ids)], "n_genes_bp": int(len(g))}, open(f"{E}/equal_info_meta.json", "w"), indent=1)
print("BayesPrism inputs written", len(g), "genes; cal donors", n25, flush=True)
