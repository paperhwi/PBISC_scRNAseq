"""Step 4: prespecified sample-level models, permutation null, verdicts, figure data.

Outputs (results/):
  effects.csv     every source x target x contrast x draw: estimate (logit diff, or z diff for B1), SE, p, n
  verdicts.csv    per target: draws positive, draws significant (BH within draw), verdict
  site_effects.csv  per-site vs HC effects
"""
import warnings
import numpy as np, pandas as pd, statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
from statsmodels.stats.multitest import multipletests

ROOT = "F:/PBISC/fig6_bulkonly_reproduction_20261020"
META = "F:/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/samples.csv"
rng = np.random.default_rng(11)
MIN_PARENT = 20
TG = ["F1", "F2", "F3", "F4", "C1"]

meta = pd.read_csv(META).rename(columns={"sample_id": "sample"})
meta = meta[meta.group.str.startswith("adult")].copy()
meta["grp"] = meta.group.map({"adult_HC": "HC", "adult_ICU_ctrl": "ICU", "adult_sepsis": "sepsis"})
meta["site"] = np.where(meta.grp == "sepsis", meta.infection_source, meta.grp)
meta["sex"] = meta.sex_geo.fillna("NA")

real_truth = pd.read_csv(f"{ROOT}/data/real_truth_fractions.csv")
real_truth = real_truth[real_truth.age_group == "adult"]
rmeta = real_truth.drop_duplicates("sample")[["sample", "group", "infection_site"]].copy()
rmeta["grp"] = rmeta.group.map({"healthy control": "HC", "ICU control": "ICU", "sepsis": "sepsis"})
SITE = {"abdominal": "Abdominal", "respiratory": "Respiratory", "urinary": "Urinary", "skin/soft tissue": "Skin", "CNS": "CNS"}
rmeta["site"] = np.where(rmeta.grp == "sepsis", rmeta.infection_site.map(SITE), rmeta.grp)


def logit_frac(k, n): return np.log((k + 0.5) / (n - k + 0.5))


def fit(d, y, covars, ref, contrast):
    d = d[d.grp.isin([ref, contrast])].copy()
    d["is_c"] = (d.grp == contrast).astype(int)
    if d.is_c.sum() < 3 or (1 - d.is_c).sum() < 3:
        return np.nan, np.nan, np.nan, int(d.is_c.sum()), int((1 - d.is_c).sum())
    # drop a covariate that is constant within this subset (e.g. one sex only)
    covars = [c for c in covars if d[c.replace("C(", "").rstrip(")")].nunique() > 1]
    f = f"{y} ~ is_c" + "".join(f" + {c}" for c in covars)
    m = smf.ols(f, data=d).fit(cov_type="HC3")
    return m.params["is_c"], m.bse["is_c"], m.pvalues["is_c"], int((d.is_c == 1).sum()), int((d.is_c == 0).sum())


def perm_p(d, y, covars, ref, contrast, obs, n=1000):
    d = d[d.grp.isin([ref, contrast])].copy(); null = []
    for _ in range(n):
        d2 = d.copy(); d2["grp"] = rng.permutation(d2.grp.to_numpy()); null.append(fit(d2, y, covars, ref, contrast)[0])
    if np.isnan(obs): return np.nan
    null = np.asarray(null); null = null[~np.isnan(null)]
    return (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)


rows, site_rows = [], []
CONTRASTS = [("HC", "sepsis"), ("ICU", "sepsis")]

# ---- real single-cell: author labels and same instrument
cr = pd.read_csv(f"{ROOT}/results/calls_real.csv")
for src, tab in [("real sc, author labels", real_truth[["sample", "target", "parent_n", "state_n"]]), ("real sc, instrument", cr)]:
    d0 = tab.merge(rmeta, on="sample")
    for k in TG:
        d = d0[(d0.target == k) & (d0.parent_n >= MIN_PARENT)].copy(); d["y"] = logit_frac(d.state_n, d.parent_n)
        for ref, c in CONTRASTS:
            e, se, p, n1, n0 = fit(d, "y", [], ref, c)
            rows.append(dict(source=src, target=k, contrast=f"{c} vs {ref}", draw=-1, est=e, se=se, p=p, n_case=n1, n_ref=n0))
        for s in ["Abdominal", "Respiratory", "Urinary", "Skin"]:
            dd = d.copy(); dd.loc[dd.site == s, "grp"] = "S"; dd = dd[(dd.grp == "HC") | (dd.site == s)]
            if (dd.grp == "S").sum() >= 3:
                e, se, p, n1, n0 = fit(dd, "y", [], "HC", "S")
                site_rows.append(dict(source=src, target=k, site=s, draw=-1, est=e, se=se, p=p, n_case=n1, n_ref=n0))

# ---- generated, per draw
cg = pd.read_csv(f"{ROOT}/results/calls_gen.csv").merge(meta, on="sample")
for j in range(5):
    for k in TG:
        d = cg[(cg.draw == j) & (cg.target == k) & (cg.parent_n >= MIN_PARENT)].copy(); d["y"] = logit_frac(d.state_n, d.parent_n)
        for ref, c in CONTRASTS:
            e, se, p, n1, n0 = fit(d, "y", ["age", "C(sex)"], ref, c)
            pp = perm_p(d, "y", ["age", "C(sex)"], ref, c, e, 500) if ref == "HC" else np.nan
            rows.append(dict(source="PBISC generated", target=k, contrast=f"{c} vs {ref}", draw=j, est=e, se=se, p=p, p_perm=pp, n_case=n1, n_ref=n0))
            ds = d[~d.sex_flag]
            e2, se2, p2, *_ = fit(ds, "y", ["age", "C(sex)"], ref, c)
            rows.append(dict(source="PBISC generated (no sex_flag)", target=k, contrast=f"{c} vs {ref}", draw=j, est=e2, se=se2, p=p2, n_case=n1, n_ref=n0))
        for s in ["Abdominal", "Respiratory", "Urinary", "Skin"]:
            dd = d[(d.grp == "HC") | (d.site == s)].copy(); dd["grp"] = np.where(dd.grp == "HC", "HC", "S")
            e, se, p, n1, n0 = fit(dd, "y", ["age", "C(sex)"], "HC", "S")
            site_rows.append(dict(source="PBISC generated", target=k, site=s, draw=j, est=e, se=se, p=p, n_case=n1, n_ref=n0))

# ---- bulk comparators
bc = pd.read_csv(f"{ROOT}/results/bulk_comparators.csv").merge(meta, on="sample")
bc["b2_logit"] = np.log(np.clip(bc.b2_parent_frac, 1e-3, 1 - 1e-3) / (1 - np.clip(bc.b2_parent_frac, 1e-3, 1 - 1e-3)))
for src, y in [("B1 bulk marker score", "b1_score"), ("B2 signature NNLS", "b2_logit")]:
    for k in TG:
        d = bc[(bc.target == k)].dropna(subset=[y])
        for ref, c in CONTRASTS:
            e, se, p, n1, n0 = fit(d, y, ["age", "C(sex)"], ref, c)
            rows.append(dict(source=src, target=k, contrast=f"{c} vs {ref}", draw=-1, est=e, se=se, p=p, n_case=n1, n_ref=n0))
        for s in ["Abdominal", "Respiratory", "Urinary", "Skin"]:
            dd = d[(d.grp == "HC") | (d.site == s)].copy(); dd["grp"] = np.where(dd.grp == "HC", "HC", "S")
            e, se, p, n1, n0 = fit(dd, y, ["age", "C(sex)"], "HC", "S")
            site_rows.append(dict(source=src, target=k, site=s, draw=-1, est=e, se=se, p=p, n_case=n1, n_ref=n0))

E = pd.DataFrame(rows)
# BH within (source, draw, contrast) over F1-F4
E["p_bh"] = np.nan
for key, d in E[E.target.isin(TG[:4]) & E.p.notna()].groupby(["source", "draw", "contrast"]):
    E.loc[d.index, "p_bh"] = multipletests(d.p, method="fdr_bh")[1]
E.to_csv(f"{ROOT}/results/effects.csv", index=False)
pd.DataFrame(site_rows).to_csv(f"{ROOT}/results/site_effects.csv", index=False)

# ---- verdicts (prespecified rule, generated, sepsis vs HC)
V = []
g = E[(E.source == "PBISC generated") & (E.contrast == "sepsis vs HC")]
ref = E[(E.source == "real sc, instrument") & (E.contrast == "sepsis vs HC")].set_index("target").est
for k in TG:
    d = g[g.target == k]
    pos = int((d.est > 0).sum()); sig = int(((d.est > 0) & (d.p_bh < 0.05)).sum()) if k != "C1" else int(((d.est < 0) & (d.p < 0.05)).sum())
    if k == "C1":
        verdict = "direction control"
    else:
        verdict = "reproduced" if (pos >= 4 and sig >= 3) else ("partial" if pos >= 4 else "not reproduced")
    V.append(dict(target=k, draws_positive=pos, draws_significant=sig, median_est=d.est.median(), real_instrument_est=ref.get(k, np.nan),
                  ratio_to_real=d.est.median() / ref.get(k, np.nan), verdict=verdict))
pd.DataFrame(V).to_csv(f"{ROOT}/results/verdicts.csv", index=False)
print(pd.DataFrame(V).to_string(index=False))
print(E[E.contrast == "sepsis vs HC"].pivot_table(index="target", columns=["source", "draw"], values="est").round(2).T.to_string())
