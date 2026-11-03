"""Score equal-information baselines against recap six-type truth on the 58 samples, alongside the PBISC results
from paired_patient_validation_20261014 (paired_sample_errors.csv.gz). Patient bootstrap 2,000."""
import json, os, numpy as np, pandas as pd
O = "/home/claude/work/prop"; E = f"{O}/eq"; PV = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014"
TYPES = ["B cell", "DC", "Monocyte", "NK cell", "Other", "T cell"]
def b6(x):
    return {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell",
            "NK cells": "NK cell", "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC"}.get(x, "Other")
truth = pd.read_csv(f"{O}/sepsis_sc_truth_props.csv.gz"); truth = truth[(truth.level == "major") & truth.grp.str.startswith("adult")]
truth["type"] = truth.label.map(b6); truth = truth.groupby(["sample", "type"]).frac.sum().unstack(fill_value=0).reindex(columns=TYPES, fill_value=0)
ids = json.load(open(f"{E}/equal_info_meta.json"))["test_samples"]
def clr(p): q = p + 1e-3; q = q / q.sum(1, keepdims=True); l = np.log(q); return l - l.mean(1, keepdims=True)
err = []
P = pd.read_csv(f"{E}/equal_info_major_props.csv")
for meth, d in P.groupby("method"):
    d = d.set_index("sample"); x = d.drop(columns="method")
    x6 = x.T.groupby(x.columns.map(b6)).sum().T.reindex(columns=TYPES, fill_value=0); x6 = x6.div(x6.sum(1), axis=0)
    t = truth.loc[x6.index]
    for s in x6.index: err.append(dict(method=meth, sample=s, L1=float(np.abs(x6.loc[s] - t.loc[s]).sum()), aitchison=float(np.sqrt(((clr(x6.loc[[s]].to_numpy()) - clr(t.loc[[s]].to_numpy())) ** 2).sum()))))
if os.path.exists(f"{O}/bp/eq25cal_theta_type.csv"):
    th = pd.read_csv(f"{O}/bp/eq25cal_theta_type.csv", index_col=0); sz = pd.read_csv(f"{O}/bp/eq25cal_size.csv")
    size = sz.groupby("major").umi.sum() / sz.groupby("major").n.sum(); cf = th.div(size.reindex(th.columns).to_numpy(), axis=1); cf = cf.div(cf.sum(1), axis=0)
    x6 = cf.T.groupby(cf.columns.map(b6)).sum().T.reindex(columns=TYPES, fill_value=0)
    for s in x6.index:
        t = truth.loc[s]; err.append(dict(method="BayesPrism+cal25", sample=s, L1=float(np.abs(x6.loc[s] - t).sum()), aitchison=float(np.sqrt(((clr(x6.loc[[s]].to_numpy()) - clr(t.to_frame().T.to_numpy())) ** 2).sum()))))
ER = pd.read_csv(f"{PV}/results/paired_sample_errors.csv.gz"); ER = ER[ER.scope == "common60_recap_truth"]
for _, r in ER.iterrows():
    meth = r.method if r.draw == -1 else f"{r.method}_draw{int(r.draw)}"
    err.append(dict(method=meth, sample=r.sample_id, L1=r.L1, aitchison=r.aitchison))
err = pd.DataFrame(err); err = err[err["sample"].isin(ids)]; err.to_csv(f"{E}/equal_info_sample_errors.csv", index=False)
rng = np.random.default_rng(20261014); rows = []
for m, d in err.groupby("method"):
    v = d.set_index("sample").L1.reindex(ids).dropna().to_numpy(); a = d.set_index("sample").aitchison.reindex(ids).dropna().to_numpy()
    b = v[rng.integers(0, len(v), (2000, len(v)))].mean(1)
    rows.append(dict(method=m, n=len(v), L1=v.mean(), lo=np.quantile(b, .025), hi=np.quantile(b, .975), aitchison=a.mean()))
S = pd.DataFrame(rows).sort_values("L1"); S.to_csv(f"{E}/equal_info_summary.csv", index=False)
# paired differences vs the PBISC K4 ensemble and vs single draws (same samples)
piv = err.pivot_table(index="sample", columns="method", values="L1")
def diff(a, b):
    x = (piv[a] - piv[b]).dropna().to_numpy(); bb = x[rng.integers(0, len(x), (2000, len(x)))].mean(1); return x.mean(), np.quantile(bb, .025), np.quantile(bb, .975), len(x)
C = []
for a, b in [("PBISC_k4000_eligible_ensemble", "WNNLS+cal ensemble"), ("PBISC_k4000_eligible_ensemble", "WNNLS+cal25"), ("PBISC_k4000_eligible_ensemble", "BayesPrism+cal25"),
             ("PBISC_k1000_eligible_ensemble", "WNNLS+cal ensemble"), ("PBISC_k4000_eligible_ensemble", "deconv_fold_null")]:
    if a in piv and b in piv: m_, lo, hi, n = diff(a, b); C.append(dict(a=a, b=b, mean_diff=m_, lo=lo, hi=hi, n=n))
pd.DataFrame(C).to_csv(f"{E}/equal_info_contrasts.csv", index=False)
# single-model level: mean over draws of PBISC single draws vs mean over draws of WNNLS+cal(draw j)
sd = [c for c in piv.columns if c.startswith("PBISC_k4000_draw")]; wd = [c for c in piv.columns if c.startswith("WNNLS+cal_draw")]
print(S.round(3).to_string()); print(pd.DataFrame(C).round(3).to_string())
print("single-draw PBISC K4 mean of draws:", round(np.nanmean([piv[c].mean() for c in sd]), 3), "range", [round(piv[c].mean(), 3) for c in sd])
print("single WNNLS+cal draw sets:", [round(piv[c].mean(), 3) for c in wd])
