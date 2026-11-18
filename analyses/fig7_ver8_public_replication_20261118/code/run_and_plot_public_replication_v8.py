"""Comprehensive External Public COVID-19 Replication of Harmonized Figure 7 (ver8) Claims (v8.1 Corrected).

Incorporates every verification point from README_Public_Replication_ver8_corrected.md:
  1. Exact sample count N = 219 acute COVID-19 patients (23 deaths) across 4 independent adult cohorts:
     Stephenson 2021 (n=73, 6 deaths; CLL donor CV0198 excluded), COMBAT 2022 (n=75, 7 deaths),
     van der Wijst 2021 (n=42, 7 deaths), and Liu 2021 (n=29, 3 deaths).
  2. Uniform cohort sets: Liu 2021 is consistently included in every meta-analysis for features it resolves
     (NK CD56bright k=4, B ASC total k=4, CD16 mono % Mono k=4, T CD8 CTL k=4, NK CD56dim k=4), with (k, N)
     shown on each row and 1-to-1 matching between individual cohort dots and the meta-analysis diamond.
  3. Exact apples-to-apples ESKD comparators: computes ESKD B ASC total (all 5 B_ASC subsets summed) to pair
     with Public B ASC total, and ESKD CD16mono|lin (% Mono) to pair with Public CD16 mono (% Mono).
  4. Transparent B ASC isotype panel (b): side-by-side display in Stephenson 2021 (n=73) of BOTH unadjusted
     Spearman rho (% PBMC, showing open non-significant markers for IgG P=0.052 & IgA P=0.12 alongside ESKD)
     AND age/sex/site-adjusted quasi-binomial count-GLM ORs (IgG OR=1.95 P=0.004, IgA OR=1.84 P=0.003).
  5. Explicit cohort cell-type proxy footnote and inclusion of both NK CD56dim and CD16 mono C1 among
     the candidates not followed in Fig. 7d.
"""
import pathlib
import re
import shutil
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
V8_STYLE_DIR = pathlib.Path(r"F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\15_Figure7_ver8_Main_and_Supp")
sys.path.insert(0, str(V8_STYLE_DIR))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

from style import FigureStyle, PALETTE, apply_publication_style, finalize_figure

ROOT = pathlib.Path(r"F:\PBISC")
P1 = ROOT / "public_covid_who_hspc_nk_20261115"
P2 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "10_CD56bright_public_replication" / "data"
ESKD_V8_DATA = V8_STYLE_DIR / "data"
OUT_DIR = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "16_Public_COVID_Replication_ver8"
DATA_DIR = OUT_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
shutil.copy2(V8_STYLE_DIR / "style.py", OUT_DIR / "style.py")

apply_publication_style(FigureStyle(font_size=13, axes_linewidth=1.8))


def age_from_stage(s):
    s = str(s)
    m = re.match(r"(\d+)-year-old", s)
    if m:
        return float(m.group(1))
    dec = {"third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75, "ninth": 85, "tenth": 95}
    for k, v in dec.items():
        if s.startswith(k):
            return float(v)
    return np.nan


def who_group(x):
    return pd.cut(x, [0.5, 3.5, 5.75, 9.5, 10.5], labels=[1, 2, 3, 4]).astype(float)


def rint(x: pd.Series) -> pd.Series:
    r = x.rank()
    return pd.Series(stats.norm.ppf((r - 0.5) / x.notna().sum()), index=x.index)


def meta_dl(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]
    if len(est) == 0:
        return np.nan, np.nan, np.nan, np.nan, 0
    if len(est) == 1:
        z = est[0] / se[0]
        return est[0], se[0], 2 * stats.norm.sf(abs(z)), 0.0, 1
    w = 1.0 / (se ** 2)
    fe = (w * est).sum() / w.sum()
    Q = (w * (est - fe) ** 2).sum()
    k = len(est)
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    ws = 1.0 / (se ** 2 + tau2)
    re_est = (ws * est).sum() / ws.sum()
    sre = np.sqrt(1.0 / ws.sum())
    z = re_est / sre
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100.0 if Q > 0 else 0.0
    return re_est, sre, 2.0 * stats.norm.sf(abs(z)), I2, k


def auc_with_se(y, d):
    y = np.asarray(y, float)
    d = np.asarray(d, bool)
    m = np.isfinite(y)
    y, d = y[m], d[m]
    n1, n0 = int(d.sum()), int((~d).sum())
    if n1 < 2 or n0 < 2:
        return np.nan, np.nan, np.nan, n1, n0
    u = stats.mannwhitneyu(y[d], y[~d])
    auc = u.statistic / (n1 * n0)
    q1 = auc / (2.0 - auc) if auc < 1.0 else 1.0
    q2 = 2.0 * (auc ** 2) / (1.0 + auc)
    var = (auc * (1.0 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)
    se = max(np.sqrt(max(var, 1e-6)), 0.035)
    return auc, se, u.pvalue, n1, n0


def count_glm(d, num_col, den_col, sev_col="who_peak"):
    sub = d[(d[den_col] > 0) & d[sev_col].notna()].copy()
    sub["y"] = sub[num_col] / sub[den_col]
    sub["sevg"] = who_group(sub[sev_col]) if sub[sev_col].max() > 7.5 else sub[sev_col]
    covs = [a for a in ("age", "sex", "site") if a in sub and sub[a].nunique() > 1]
    sub = sub.dropna(subset=["sevg"] + [a for a in covs if a != "sex"])
    if len(sub) < 12 or sub[num_col].sum() < 5:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    f = "y ~ sevg" + "".join(f" + C({a})" if a in ("sex", "site") else f" + {a}" for a in covs)
    r = smf.glm(f, sub, family=sm.families.Binomial(), var_weights=sub[den_col]).fit(scale="X2")
    b, se = float(r.params["sevg"]), float(r.bse["sevg"])
    return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), se, float(r.pvalues["sevg"])


# -------------------------------------------------------------------------
# 1. Build Harmonized Sample Tables Across 4 Independent Adult Public Cohorts
# -------------------------------------------------------------------------
# (A) Stephenson et al. 2021 Nat Med (exclude B_malignant CLL donor CV0198)
s_raw = pd.read_parquet(P1 / "obs_stephenson2021.parquet")
s = s_raw[(s_raw["Status"].astype(str) == "Covid") & (s_raw["donor_id"].astype(str) != "CV0198")].copy()
ct = s["author_cell_type"].astype(str)
den_s = ~ct.isin(["RBC", "Platelets", "B_malignant"])
s_masks = {
    "NK_CD56bright": ct == "NK_56hi",
    "NK_CD56dim": ct == "NK_16hi",
    "NK_total": ct.isin(["NK_56hi", "NK_16hi", "NK_prolif"]),
    "B_ASC_dividing": ct == "Plasmablast",
    "B_ASC_IgG": ct == "Plasma_cell_IgG",
    "B_ASC_IgA": ct == "Plasma_cell_IgA",
    "B_ASC_IgM": ct == "Plasma_cell_IgM",
    "B_ASC_switched": ct.isin(["Plasma_cell_IgG", "Plasma_cell_IgA"]),
    "B_ASC_total": ct.isin(["Plasmablast", "Plasma_cell_IgG", "Plasma_cell_IgA", "Plasma_cell_IgM"]),
    "B_total": (ct.str.startswith("B_") & (ct != "B_malignant")) | ct.isin(["Plasmablast", "Plasma_cell_IgG", "Plasma_cell_IgA", "Plasma_cell_IgM"]),
    "T_CD8_CM": ct == "CD8.EM",  # CD8 memory proxy (no separate CD8.CM label in Stephenson)
    "T_CD8_CTL": ct == "CD8.TE",
    "CD16_mono": ct == "CD16_mono",
    "CD16_mono_C1": ct == "C1_CD16_mono",
    "Mono_total": ct.str.contains("mono"),
    "HSPC": ct.str.startswith("HSC"),
}
gs = s.assign(**s_masks, _den=den_s).groupby("sample_id", observed=True)[list(s_masks) + ["_den"]].sum()
ms = s.drop_duplicates("sample_id").set_index("sample_id")
peak_s = {"Asymptomatic": 1, "Mild": 4, "Moderate": 5, "Severe": 6, "Critical ": 7, "Critical": 7, "Death": 10}
gs["cohort"] = "Stephenson 2021"
gs["donor"] = ms["donor_id"].astype(str)
gs["who_peak"] = ms["Worst_Clinical_Status"].astype(str).map(peak_s)
gs["death"] = ms["Worst_Clinical_Status"].astype(str).eq("Death") | ms["Outcome"].astype(str).eq("Death")
gs["age"] = ms["development_stage"].map(age_from_stage)
gs["sex"] = ms["sex"].astype(str)
gs["site"] = ms["Site"].astype(str)
gs["col_day"] = ms["Collection_Day"].astype(str).str.lstrip("D").astype(float)
gs["days_onset"] = pd.to_numeric(ms["Days_from_onset"].astype(str), errors="coerce") + gs["col_day"]
gs["acute"] = gs["col_day"] < 28
gs = gs[gs["_den"] >= 500].copy()

# (B) COMBAT Consortium 2022 Cell
c_raw = pd.read_parquet(P1 / "obs_combat2022.parquet")
c = c_raw[c_raw["Source"].astype(str).isin(["COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"])].copy()
mi = c["minor_subset"].astype(str)
ma = c["major_subset"].astype(str)
den_c = ~ma.isin(["PLT", "RET", "nan"])
c_masks = {
    "NK_CD56bright": mi == "NK.CD56hi.CD16int.XCL1/2",
    "NK_CD56dim": mi.isin(["NK.CD16hi", "NK.CD16int"]),
    "NK_total": ma == "NK",
    "B_ASC_total": ma == "PB",
    "B_total": ma.isin(["B", "PB"]),
    "T_CD8_CM": mi.isin(["CD8.TCM", "CD8.TCM.CCL5"]),
    "T_CD8_CTL": mi.isin(["CD8.TEMRA", "CD8.TEFF", "CD8.TEM"]),
    "CD16_mono": ma == "ncMono",
    "Mono_total": ma.isin(["cMono", "ncMono"]),
    "HSPC": ma == "HSC",
}
gc = c.assign(**c_masks, _den=den_c).groupby("scRNASeq_sample_ID", observed=True)[list(c_masks) + ["_den"]].sum()
mc = c.drop_duplicates("scRNASeq_sample_ID").set_index("scRNASeq_sample_ID")
peak_c = {"Not hospitalized": 2, "Hospitalized, no O2": 4, "Hospitalized, O2": 5, "Non-invasive ventilation": 6, "Intubated, ventilated": 7, "Death": 10}
agebin = {"19-30": 25, "31-40": 35, "41-50": 45, "51-60": 55, "61-70": 65, "71-80": 75, "81-90": 85, ">=91": 92}
gc["cohort"] = "COMBAT 2022"
gc["donor"] = mc["donor_id"].astype(str)
gc["who_peak"] = mc["Outcome"].astype(str).map(peak_c)
asym = (mc["Source"].astype(str) == "COVID_HCW_MILD") & (~mc["Symptomatic"].astype(bool))
gc.loc[asym, "who_peak"] = 1
gc["death"] = mc["Outcome"].astype(str).eq("Death") | mc["Death28"].astype(bool)
gc["age"] = mc["Age"].astype(str).map(agebin)
gc["sex"] = mc["sex"].astype(str)
gc["site"] = mc["Institute"].astype(str)
gc["days_onset"] = pd.to_numeric(mc["TimeSinceOnset"], errors="coerce")
gc["acute"] = gc["days_onset"].isna() | (gc["days_onset"] <= 28)
gc = gc[gc["_den"] >= 500].copy()

# (C) van der Wijst et al. 2021 Sci Transl Med
v_raw = pd.read_parquet(P1 / "obs_vanderwijst2021.parquet")
v = v_raw[v_raw["COVID_status"].astype(str) == "Positive"].copy()
v["sample_id"] = v["donor_id"].astype(str) + "_t" + v["timepoint"].astype(str)
v["nih"] = pd.to_numeric(v["NIH_ordinal"].astype(str), errors="coerce")
pk_v = v.groupby("donor_id", observed=True)["nih"].max().to_dict()
c2 = v["ct2"].astype(str)
c1 = v["ct1"].astype(str)
v_masks = {
    "NK_CD56bright": c2 == "NK_CD56++",
    "NK_CD56dim": c2 == "NK_CD16+",
    "NK_total": c1 == "NK",
    "B_ASC_dividing": c2 == "PB_Prolif",
    "B_ASC_switched": c2 == "PB",
    "B_ASC_total": c1 == "PB",
    "B_total": c1.isin(["B", "PB"]),
    "T_CD8_CTL": c2 == "T8_Mem",  # CD8 memory/effector proxy
    "CD16_mono": c2 == "ncM",
    "Mono_total": c1.isin(["cM", "ncM"]),
    "HSPC": c1 == "Progen",
}
gv = v.assign(**v_masks, _den=True).groupby("sample_id", observed=True)[list(v_masks) + ["_den"]].sum()
mv = v.drop_duplicates("sample_id").set_index("sample_id")
gv["cohort"] = "van der Wijst 2021"
gv["donor"] = mv["donor_id"].astype(str)
gv["death"] = mv["death"].astype(str).eq("yes")
gv["who_peak"] = gv["donor"].map(pk_v)
gv.loc[gv["death"] & gv["who_peak"].notna(), "who_peak"] = 10
gv["age"] = pd.to_numeric(mv["donor_age"].astype(str), errors="coerce")
gv["sex"] = mv["sex"].astype(str)
gv["site"] = "UCSF"
gv["tp_num"] = pd.to_numeric(mv["days_from_start"], errors="coerce").fillna(0)
gv["days_onset"] = pd.to_numeric(mv["onset_to_D0_days-as_float"], errors="coerce") + gv["tp_num"]
gv["acute"] = gv["tp_num"] <= 14
gv = gv[gv["_den"] >= 500].copy()

# (D) Liu et al. 2021 Cell (complete innate+adaptive PBMC samples)
L_raw = pd.read_parquet(P2 / "obs_liu2021.parquet")
L = L_raw[L_raw["timepoint"].astype(str) != "HC"].copy()
L["sample_id"] = L["donor_id"].astype(str) + "_" + L["timepoint"].astype(str)
la = L["author_cell_type"].astype(str)
den_L = ~la.isin(["platelet", "granulocyte"])
L_masks = {
    "NK_CD56bright": la == "NK_CD56hiCD16lo",
    "NK_CD56dim": la == "NK_CD16hi",
    "NK_total": la.str.startswith("NK"),
    "B_ASC_total": la == "plasmablast",
    "T_CD8_CTL": la == "CD8-positive, alpha-beta memory T cell",
    "CD16_mono": la == "non-classical monocyte",
    "Mono_total": la.str.contains("monocyte"),
}
gL = L.assign(**L_masks, _den=den_L).groupby("sample_id", observed=True)[list(L_masks) + ["_den"]].sum()
part_cnt = L.groupby("sample_id", observed=True)["part"].nunique()
gL = gL.loc[part_cnt[part_cnt == 2].index].copy()
mL = L.drop_duplicates("sample_id").set_index("sample_id").loc[gL.index]
sev_L = {"Moderate": 4.5, "Severe": 5.5, "Critical": 7.0}
gL["cohort"] = "Liu 2021"
gL["donor"] = mL["donor_id"].astype(str)
gL["death"] = mL["outcome"].astype(str).eq("deceased")
gL["who_peak"] = mL["severity"].astype(str).map(sev_L)
gL.loc[gL["death"], "who_peak"] = 10.0
gL["age"] = pd.to_numeric(mL["age"].astype(str), errors="coerce")
gL["sex"] = mL["sex"].astype(str)
gL["site"] = "NIAID"
gL["days_onset"] = pd.to_numeric(mL["days_since_onset"].astype(str), errors="coerce")
gL["acute"] = gL["days_onset"].isna() | (gL["days_onset"] <= 28)
gL = gL[gL["_den"] >= 500].copy()

COHORT_DFS = {
    "Stephenson 2021": gs,
    "COMBAT 2022": gc,
    "van der Wijst 2021": gv,
    "Liu 2021": gL,
}

# -------------------------------------------------------------------------
# 2. Compute Patient-Level Acute Means & Cohort-Specific + Uniform Meta-Analytic Stats
# -------------------------------------------------------------------------
FEATS_EVAL = [
    "NK_CD56bright", "NK_CD56dim", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM",
    "B_ASC_switched", "B_ASC_total", "T_CD8_CM", "T_CD8_CTL", "CD16_mono", "CD16_mono_C1", "HSPC"
]

coh_rows = []
for cname, gdf in COHORT_DFS.items():
    sub = gdf[gdf["acute"]].copy()
    present = [f for f in FEATS_EVAL if f in sub.columns]
    for f in present:
        sub[f + "_pct"] = 100.0 * sub[f] / sub["_den"]
    if "CD16_mono" in sub.columns and "Mono_total" in sub.columns:
        sub["CD16_mono_ofMono_pct"] = 100.0 * sub["CD16_mono"] / sub["Mono_total"].replace(0, np.nan)
        present.append("CD16_mono_ofMono")

    pct_cols = [f + "_pct" for f in present]
    raw_cols = [f for f in FEATS_EVAL if f in sub.columns] + ["_den"]
    P = sub.groupby("donor")[pct_cols + raw_cols].mean()
    P["who_peak"] = sub.groupby("donor")["who_peak"].first()
    P["age"] = sub.groupby("donor")["age"].first()
    P["sex"] = sub.groupby("donor")["sex"].first()
    P["site"] = sub.groupby("donor")["site"].first()
    P["death"] = sub.groupby("donor")["death"].max()

    P_sev = P.dropna(subset=["who_peak"])
    P_sc = P_sev[P_sev["who_peak"] >= 6]

    for f in present:
        col = f + "_pct"
        ok = P_sev[col].notna()
        n_s = int(ok.sum())
        rho, p_s = stats.spearmanr(P_sev.loc[ok, "who_peak"], P_sev.loc[ok, col])
        rng = np.random.default_rng(42)
        yv = P_sev.loc[ok, col].values
        sv = P_sev.loc[ok, "who_peak"].values
        bs = [stats.spearmanr(sv[ix := rng.integers(0, n_s, n_s)], yv[ix])[0] for _ in range(1500)]
        rho_lo, rho_hi = float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))

        covs = [a for a in ("age", "sex", "site") if P_sev[a].nunique() > 1]
        e = P_sev.dropna(subset=[col, "who_peak"] + [a for a in covs if a != "sex"]).copy()
        e["ry"] = e[col].rank()
        e["rs"] = e["who_peak"].rank()
        p_adj = float(smf.ols("ry ~ rs" + "".join(f" + C({a})" if a in ("sex", "site") else f" + {a}" for a in covs), e).fit().pvalues["rs"])

        if f in raw_cols:
            or_v, or_lo, or_hi, se_or, p_glm = count_glm(P_sev, f, "_den", "who_peak")
        else:
            or_v = or_lo = or_hi = se_or = p_glm = np.nan

        auc_all, se_all, pd_all, n1a, n0a = auc_with_se(P[col], P["death"])
        auc_sc, se_sc, pd_sc, n1s, n0s = auc_with_se(P_sc[col], P_sc["death"])

        coh_rows.append({
            "cohort": cname,
            "feature": f,
            "n_sev": n_s,
            "rho": rho,
            "rho_lo": rho_lo,
            "rho_hi": rho_hi,
            "p_sev": p_s,
            "p_sev_adj": p_adj,
            "OR_glm": or_v,
            "OR_lo": or_lo,
            "OR_hi": or_hi,
            "p_glm": p_glm,
            "n_death_all": n1a,
            "n_total_all": n1a + n0a,
            "AUC_death_all": auc_all,
            "SE_death_all": se_all,
            "p_death_all": pd_all,
            "n_death_WHO6": n1s,
            "n_total_WHO6": n1s + n0s,
            "AUC_death_WHO6": auc_sc,
            "SE_death_WHO6": se_sc,
            "p_death_WHO6": pd_sc,
        })

df_coh = pd.DataFrame(coh_rows)
df_coh.to_csv(DATA_DIR / "public_cohort_level_results_v8.csv", index=False)

st_b = df_coh[(df_coh["cohort"] == "Stephenson 2021") & (df_coh["feature"].str.startswith("B_ASC"))].copy()
st_b.to_csv(DATA_DIR / "stephenson_b_asc_isotype_glm_v8.csv", index=False)

# UNIFORM COHORT INCLUSION: Every cohort resolving a feature is included in its meta-analysis!
ALL4 = ["Stephenson 2021", "COMBAT 2022", "van der Wijst 2021", "Liu 2021"]
META_SPECS = [
    ("NK_CD56bright", "NK CD56bright", "Concordant anchor", ALL4),
    ("B_ASC_dividing", "B ASC dividing", "Concordant anchor", ["Stephenson 2021", "van der Wijst 2021"]),
    ("T_CD8_CM", "T CD8 CM (EM proxy)", "Anchor (sev) / Focus (death)", ["Stephenson 2021", "COMBAT 2022"]),
    ("CD16_mono_ofMono", "CD16 mono (% Mono)", "Concordant anchor", ALL4),
    ("CD16_mono", "CD16 mono (% PBMC)", "Concordant anchor", ALL4),
    ("B_ASC_IgG", "B ASC IgG (Steph.)", "Focus: severity", ["Stephenson 2021"]),
    ("B_ASC_IgA", "B ASC IgA (Steph.)", "Focus: severity", ["Stephenson 2021"]),
    ("B_ASC_total", "B ASC total", "Focus: severity", ALL4),
    ("T_CD8_CTL", "T CD8 CTL (eff. proxy)", "Focus: sev & death", ALL4),
    ("HSPC", "HSPC", "Focus: sev & death", ["Stephenson 2021", "COMBAT 2022", "van der Wijst 2021"]),
    ("NK_CD56dim", "NK CD56dim", "Not followed in Fig. 7d", ALL4),
    ("CD16_mono_C1", "CD16 mono C1 (Steph.)", "Not followed in Fig. 7d", ["Stephenson 2021"]),
]

meta_rows = []
for fkey, flab, role, clist in META_SPECS:
    sub = df_coh[(df_coh["feature"] == fkey) & (df_coh["cohort"].isin(clist))]
    z_arr = np.arctanh(sub["rho"].values)
    sz_arr = np.sqrt(1.06 / (sub["n_sev"].values - 3))
    zm, szm, pm, I2_s, k_s = meta_dl(z_arr, sz_arr)

    am_all, sam_all, pm_all, _, _ = meta_dl(sub["AUC_death_all"].values - 0.5, sub["SE_death_all"].values)
    am_w6, sam_w6, pm_w6, _, _ = meta_dl(sub["AUC_death_WHO6"].values - 0.5, sub["SE_death_WHO6"].values)

    meta_rows.append({
        "feature": fkey,
        "label": flab,
        "role": role,
        "k_cohorts": k_s,
        "N_patients": int(sub["n_sev"].sum()),
        "N_deaths": int(sub["n_death_all"].sum()),
        "rho_meta": float(np.tanh(zm)),
        "rho_lo": float(np.tanh(zm - 1.96 * szm)),
        "rho_hi": float(np.tanh(zm + 1.96 * szm)),
        "p_rho_meta": pm,
        "I2_rho": I2_s,
        "AUC_death_all_meta": float(am_all + 0.5),
        "AUC_all_lo": float(np.clip(am_all + 0.5 - 1.96 * sam_all, 0, 1)),
        "AUC_all_hi": float(np.clip(am_all + 0.5 + 1.96 * sam_all, 0, 1)),
        "p_death_all_meta": pm_all,
        "AUC_death_WHO6_meta": float(am_w6 + 0.5),
        "AUC_WHO6_lo": float(np.clip(am_w6 + 0.5 - 1.96 * sam_w6, 0, 1)),
        "AUC_WHO6_hi": float(np.clip(am_w6 + 0.5 + 1.96 * sam_w6, 0, 1)),
        "p_death_WHO6_meta": pm_w6,
    })

df_meta = pd.DataFrame(meta_rows)
df_meta.to_csv(DATA_DIR / "public_meta_analysis_v8.csv", index=False)

# HSPC Time-window mortality comparison (First sample vs Late > day 7) + Meta-analysis
hw_rows = []
for win_key, win_label in [("first", "First acute sample"), ("late_gt7d", "Samples > day 7")]:
    ests_all, ses_all, ests_w6, ses_w6 = [], [], [], []
    n_d_tot, n_p_tot = 0, 0
    for cname in ["Stephenson 2021", "COMBAT 2022", "van der Wijst 2021"]:
        gdf = COHORT_DFS[cname].copy()
        gdf["HSPC_pct"] = 100.0 * gdf["HSPC"] / gdf["_den"]
        ac = gdf[gdf["acute"]].sort_values(["donor", "days_onset"])
        if win_key == "first":
            dw = ac.groupby("donor").head(1).set_index("donor")
        else:
            sub7 = ac[ac["days_onset"] > 7]
            dw = sub7.groupby("donor")[["HSPC_pct"]].mean()
            dw["death"] = sub7.groupby("donor")["death"].max()
            dw["who_peak"] = sub7.groupby("donor")["who_peak"].first()
        dw = dw.dropna(subset=["HSPC_pct"])
        a_all, se_all, p_all, n1, n0 = auc_with_se(dw["HSPC_pct"], dw["death"])
        dw6 = dw[dw["who_peak"] >= 6]
        a_w6, se_w6, p_w6, n1s, n0s = auc_with_se(dw6["HSPC_pct"], dw6["death"])
        hw_rows.append(dict(cohort=cname, window=win_label, deaths=n1, n=n1 + n0, AUC_all=a_all, SE_all=se_all, p_all=p_all,
                            deaths_WHO6=n1s, n_WHO6=n1s + n0s, AUC_WHO6=a_w6, SE_WHO6=se_w6, p_WHO6=p_w6))
        ests_all.append(a_all - 0.5)
        ses_all.append(se_all)
        ests_w6.append(a_w6 - 0.5)
        ses_w6.append(se_w6)
        n_d_tot += n1
        n_p_tot += (n1 + n0)
    ma, sa, pa, _, _ = meta_dl(ests_all, ses_all)
    mw, sw, pw, _, _ = meta_dl(ests_w6, ses_w6)
    hw_rows.append(dict(cohort="Public Meta (3 cohorts)", window=win_label, deaths=n_d_tot, n=n_p_tot,
                        AUC_all=ma + 0.5, SE_all=sa, p_all=pa, deaths_WHO6=np.nan, n_WHO6=np.nan, AUC_WHO6=mw + 0.5, SE_WHO6=sw, p_WHO6=pw))

df_hw = pd.DataFrame(hw_rows)
df_hw.to_csv(DATA_DIR / "public_hspc_timewindow_death_v8.csv", index=False)


# -------------------------------------------------------------------------
# 3. Compute Exact Matched ESKD Comparators (including B_ASC_total & CD16mono|lin)
# -------------------------------------------------------------------------
pt_raw_eskd = pd.read_csv(ESKD_V8_DATA / "patient_raw_pct_noCLL_v8.csv", index_col=0)
b_asc_cols_eskd = [c for c in ["B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM", "B_ASC_other"] if c in pt_raw_eskd.columns]
pt_raw_eskd["B_ASC_total"] = pt_raw_eskd[b_asc_cols_eskd].sum(axis=1)

ESKD_MATCH_FEATS = [
    "NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono|lin", "CD16mono",
    "B_ASC_total", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM", "T CD8 CTL", "HSPC", "NK", "CD16mono_C1"
]
rng_e = np.random.default_rng(20261118)
eskd_matched = {}
R_e = pt_raw_eskd[pt_raw_eskd["src"] == "Real"]
N_e = pt_raw_eskd[pt_raw_eskd["src"] == "New"]
for fk in ESKD_MATCH_FEATS:
    comb = pd.concat([
        pd.DataFrame({"y": rint(R_e[fk]), "sev": R_e["sev"], "fatal": R_e["fatal"]}),
        pd.DataFrame({"y": rint(N_e[fk]), "sev": N_e["sev"], "fatal": N_e["fatal"]}),
    ]).dropna()
    r_v, p_v = stats.spearmanr(comb["sev"], comb["y"])
    yv, sv = comb["y"].values, comb["sev"].values
    n_c = len(comb)
    bs_e = [stats.spearmanr(sv[ix := rng_e.integers(0, n_c, n_c)], yv[ix])[0] for _ in range(2000)]
    eskd_matched[fk] = dict(rho=float(r_v), lo=float(np.percentile(bs_e, 2.5)), hi=float(np.percentile(bs_e, 97.5)), p=float(p_v))

print("=== UNIFORM PUBLIC META-ANALYSIS SUMMARY (N = 219) ===")
print(df_meta.round(4).to_string(index=False))


# =========================================================================
# 4. RENDER CORRECTED EXTERNAL REPLICATION FIGURE (Supp Fig S3)
# =========================================================================
fig = plt.figure(figsize=(15.4, 14.4), facecolor="white")
fig.text(
    0.015,
    0.988,
    "Supplementary Figure 3. External single-cell replication of Figure 7 findings across four adult COVID-19 PBMC cohorts (N = 219 patients)",
    fontsize=13.4,
    fontweight="bold",
    va="top",
    ha="left",
    color=PALETTE["gray_3"],
)

# ---------------- Panel a: Matched ESKD vs. Uniform Public Meta-Analysis (Severity) ----------------
ax_a = fig.add_axes([0.155, 0.535, 0.285, 0.375])
ORDER_A = [
    ("NK_CD56bright", "NK CD56bright", "NK CD56 bright", "anchor"),
    ("B_ASC_dividing", "B ASC dividing*", "B_ASC_dividing", "anchor"),
    ("T_CD8_CM", "T CD8 CM / EM*", "T CD8 CM", "anchor"),
    ("CD16_mono_ofMono", "CD16 mono (% Mono)", "CD16mono|lin", "anchor"),
    ("B_ASC_total", "B ASC total", "B_ASC_total", "focus"),
    ("T_CD8_CTL", "T CD8 CTL / Mem*", "T CD8 CTL", "focus"),
    ("HSPC", "HSPC", "HSPC", "hspc"),
    ("NK_CD56dim", "NK CD56dim", "NK", "excl"),
    ("CD16_mono_C1", "CD16 mono C1", "CD16mono_C1", "excl"),
]

y_pos = []
yc = 0.0
for idx, (pub_f, disp, eskd_f, grp) in enumerate(ORDER_A):
    if idx in (4, 7):
        ax_a.axhline(yc + 0.52, color=PALETTE["neutral"], ls="--", lw=1.0, zorder=1)
        yc -= 0.28
    y_pos.append(yc)

    em = eskd_matched[eskd_f]
    ax_a.plot([em["lo"], em["hi"]], [yc + 0.18, yc + 0.18], color=PALETTE["blue_main"], lw=1.7, zorder=2)
    ax_a.scatter(em["rho"], yc + 0.18, s=44, marker="s", c="white" if em["p"] >= 0.05 else PALETTE["blue_main"],
                 edgecolors=PALETTE["blue_main"], linewidths=1.4, zorder=4)

    mr = df_meta[df_meta["feature"] == pub_f].iloc[0]
    sub_c = df_coh[df_coh["feature"] == pub_f]
    for _, cr in sub_c.iterrows():
        ax_a.scatter(cr["rho"], yc - 0.18, s=20, marker="o", c=PALETTE["neutral"], edgecolors=PALETTE["gray_1"], linewidths=0.8, alpha=0.85, zorder=2)

    m_col = PALETTE["red_strong"] if grp == "hspc" else (PALETTE["gray_2"] if grp == "excl" else PALETTE["teal"])
    ax_a.plot([mr["rho_lo"], mr["rho_hi"]], [yc - 0.18, yc - 0.18], color=m_col, lw=2.0, zorder=3)
    ax_a.scatter(mr["rho_meta"], yc - 0.18, s=50, marker="D", c="white" if mr["p_rho_meta"] >= 0.05 else m_col,
                 edgecolors=m_col, linewidths=1.5, zorder=5)

    p_m = mr["p_rho_meta"]
    p_txt = f"P={p_m:.1e}".replace("e-0", "e-") if p_m < 0.001 else (f"P={p_m:.3f}" if p_m < 0.1 else f"P={p_m:.2f}")
    ax_a.text(0.83, yc - 0.10, f"{mr['rho_meta']:+.2f} ({p_txt}) · k={int(mr['k_cohorts'])}, N={int(mr['N_patients'])}",
              va="center", ha="left", fontsize=7.9,
              color=m_col if p_m < 0.05 else PALETTE["gray_1"], fontweight="bold" if p_m < 0.05 else "normal")
    yc -= 1.05

ax_a.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
ax_a.set_yticks(y_pos)
ax_a.set_yticklabels([x[1] for x in ORDER_A], fontsize=9.4)
for idx_t, t in enumerate(ax_a.get_yticklabels()):
    grp = ORDER_A[idx_t][3]
    t.set_fontweight("bold" if grp != "excl" else "normal")
    t.set_color(PALETTE["red_strong"] if grp == "hspc" else (PALETTE["blue_main"] if grp == "focus" else (PALETTE["gray_1"] if grp == "excl" else PALETTE["gray_3"])))
ax_a.tick_params(axis="y", length=0)
ax_a.spines["left"].set_visible(False)
ax_a.set_xlim(-0.82, 0.80)
ax_a.set_ylim(y_pos[-1] - 0.55, 0.58)
ax_a.set_xticks([-0.5, 0.0, 0.5])
ax_a.set_xlabel("Spearman ρ with WHO severity", fontsize=10.2)
ax_a.text(0.83, 0.50, "Public Meta ρ (P) · cohorts", fontsize=8.1, fontweight="bold", color=PALETTE["gray_2"], va="bottom")
ax_a.set_title("Matched severity comparison", loc="left", fontsize=11.3, fontweight="bold", pad=8)
fig.text(0.015, 0.930, "a", fontsize=18, fontweight="bold", va="top")


# ---------------- Panel b: Stephenson 2021 B ASC Isotypes (Spearman rho AND Count-GLM OR) ----------------
ax_b1 = fig.add_axes([0.670, 0.575, 0.132, 0.335])
ax_b2 = fig.add_axes([0.830, 0.575, 0.140, 0.335])

ISO_ORDER = [
    ("B_ASC_dividing", "Plasmablast (dividing)", "B_ASC_dividing"),
    ("B_ASC_IgG", "Plasma cell IgG", "B_ASC_IgG"),
    ("B_ASC_IgA", "Plasma cell IgA", "B_ASC_IgA"),
    ("B_ASC_IgM", "Plasma cell IgM", "B_ASC_IgM"),
    ("B_ASC_total", "B ASC total", "B_ASC_total"),
]
yy_b = np.arange(len(ISO_ORDER))[::-1]
ax_b1.axvline(0.0, color=PALETTE["gray_1"], ls=":", lw=1.3)
ax_b2.axvline(1.0, color=PALETTE["gray_1"], ls=":", lw=1.3)

for yi, (fk, flab, efk) in zip(yy_b, ISO_ORDER):
    r = st_b[st_b["feature"] == fk].iloc[0]
    em = eskd_matched[efk]
    col_b = PALETTE["blue_main"] if fk in ("B_ASC_IgG", "B_ASC_IgA") else PALETTE["teal"]

    # Left sub-axis: ESKD (square) vs Stephenson 2021 (diamond) Spearman rho (% PBMC)
    ax_b1.plot([em["lo"], em["hi"]], [yi + 0.16, yi + 0.16], color=PALETTE["blue_main"], lw=1.6)
    ax_b1.scatter(em["rho"], yi + 0.16, s=42, marker="s", c=PALETTE["blue_main"], edgecolors=PALETTE["blue_main"], zorder=3)

    ax_b1.plot([r["rho_lo"], r["rho_hi"]], [yi - 0.16, yi - 0.16], color=PALETTE["teal"], lw=1.8)
    ax_b1.scatter(r["rho"], yi - 0.16, s=48, marker="D",
                  c="white" if r["p_sev"] >= 0.05 else PALETTE["teal"], edgecolors=PALETTE["teal"], linewidths=1.5, zorder=4)
    ax_b1.text(0.98, yi - 0.16, f"P={r['p_sev']:.2f}" if r["p_sev"] >= 0.05 else f"P={r['p_sev']:.3f}",
               va="center", ha="right", fontsize=7.6, color=PALETTE["gray_1"] if r["p_sev"] >= 0.05 else PALETTE["teal"],
               fontweight="bold" if r["p_sev"] < 0.05 else "normal")

    # Right sub-axis: Stephenson 2021 Age/Sex/Site Count-GLM OR
    ax_b2.plot([r["OR_lo"], r["OR_hi"]], [yi, yi], color=col_b, lw=2.0)
    ax_b2.scatter(r["OR_glm"], yi, s=52, marker="D", c=col_b, edgecolors="white", linewidths=0.7, zorder=4)
    ax_b2.text(r["OR_hi"] + 0.08, yi, f"{r['OR_glm']:.2f} (P={r['p_glm']:.3f})", va="center", ha="left",
               fontsize=7.8, color=col_b, fontweight="bold")

ax_b1.set_yticks(yy_b)
ax_b1.set_yticklabels([x[1] for x in ISO_ORDER], fontsize=9.4)
for idx_t, t in enumerate(ax_b1.get_yticklabels()):
    if ISO_ORDER[idx_t][0] in ("B_ASC_IgG", "B_ASC_IgA"):
        t.set_fontweight("bold")
        t.set_color(PALETTE["blue_main"])
ax_b2.set_yticks(yy_b)
ax_b2.set_yticklabels([])

for axx in (ax_b1, ax_b2):
    axx.tick_params(axis="y", length=0)
    axx.spines["left"].set_visible(False)
    axx.set_ylim(-0.55, len(ISO_ORDER) - 0.35)

ax_b1.set_xlim(-0.20, 1.02)
ax_b1.set_xticks([0.0, 0.4, 0.8])
ax_b1.set_xlabel("Spearman ρ (% PBMC)", fontsize=9.8)
ax_b1.set_title("Stephenson 2021 B ASC isotypes (n = 73): rank ρ vs. count GLM", loc="left", fontsize=11.1, fontweight="bold", pad=8)

ax_b2.set_xlim(0.85, 4.95)
ax_b2.set_xticks([1.0, 2.0, 3.0, 4.0])
ax_b2.set_xlabel("Adj. count-GLM OR / step", fontsize=9.8)
fig.text(0.545, 0.930, "b", fontsize=18, fontweight="bold", va="top")

# Legend bar between Row 1 and Row 2
ax_leg = fig.add_axes([0.02, 0.462, 0.96, 0.026])
ax_leg.set_axis_off()
hd = [
    Line2D([0], [0], marker="s", ls="", ms=7.5, mfc=PALETTE["blue_main"], mec=PALETTE["blue_main"], label="ESKD Real + PBISC new (n = 66, 13 deaths)"),
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc=PALETTE["teal"], mec=PALETTE["teal"], label="Public COVID-19 Meta-analysis / Stephenson 2021 (N = 219 across 4 cohorts)"),
    Line2D([0], [0], marker="o", ls="", ms=6.2, mfc=PALETTE["neutral"], mec=PALETTE["gray_1"], label="Individual public cohort (a)"),
    Line2D([0], [0], marker="D", ls="", ms=7.5, mfc="white", mec=PALETTE["gray_2"], mew=1.6, label="open: P ≥ 0.05"),
]
ax_leg.legend(handles=hd, loc="center", ncol=4, fontsize=9.5, frameon=False, columnspacing=1.3)


# ---------------- Panel c: Death AUC: All COVID+ vs Within Severe/Critical (WHO>=6) ----------------
ax_c = fig.add_axes([0.155, 0.095, 0.285, 0.320])
D_ITEMS = [
    ("HSPC", "HSPC (k=3)", PALETTE["red_strong"]),
    ("T_CD8_CM", "T CD8 CM / EM* (k=2)", PALETTE["blue_main"]),
    ("T_CD8_CTL", "T CD8 CTL / Mem* (k=4)", PALETTE["blue_main"]),
    ("B_ASC_total", "B ASC total (k=4)", PALETTE["gray_2"]),
    ("NK_CD56dim", "NK CD56dim (not fol., k=4)", PALETTE["gray_1"]),
    ("CD16_mono_C1", "CD16 mono C1 (not fol., k=1)", PALETTE["gray_1"]),
]
yy_c = np.arange(len(D_ITEMS))[::-1]
ax_c.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)

for yi, (pub_f, disp, col_c) in zip(yy_c, D_ITEMS):
    mr = df_meta[df_meta["feature"] == pub_f].iloc[0]
    ax_c.plot([mr["AUC_all_lo"], mr["AUC_all_hi"]], [yi + 0.17, yi + 0.17], color=col_c, lw=1.9)
    ax_c.scatter(mr["AUC_death_all_meta"], yi + 0.17, s=50, marker="D",
                 c="white" if mr["p_death_all_meta"] >= 0.05 else col_c, edgecolors=col_c, linewidths=1.5, zorder=4)
    ax_c.plot([mr["AUC_WHO6_lo"], mr["AUC_WHO6_hi"]], [yi - 0.17, yi - 0.17], color=col_c, lw=1.5, ls="--", alpha=0.8)
    ax_c.scatter(mr["AUC_death_WHO6_meta"], yi - 0.17, s=46, marker="^",
                 c="white" if mr["p_death_WHO6_meta"] >= 0.05 else col_c, edgecolors=col_c, linewidths=1.4, zorder=4)

    p1_t = f"{mr['p_death_all_meta']:.3f}" if mr["p_death_all_meta"] < 0.1 else f"{mr['p_death_all_meta']:.2f}"
    p2_t = f"{mr['p_death_WHO6_meta']:.2f}"
    ax_c.text(0.88, yi + 0.17, f"All: {mr['AUC_death_all_meta']:.2f} (P={p1_t})", va="center", fontsize=8.0,
              color=col_c if mr["p_death_all_meta"] < 0.05 else PALETTE["gray_2"],
              fontweight="bold" if mr["p_death_all_meta"] < 0.05 else "normal")
    ax_c.text(0.88, yi - 0.17, f"WHO≥6: {mr['AUC_death_WHO6_meta']:.2f} (P={p2_t})", va="center", fontsize=7.8, color=PALETTE["gray_1"])

ax_c.set_yticks(yy_c)
ax_c.set_yticklabels([x[1] for x in D_ITEMS], fontsize=9.2)
for idx_t, t in enumerate(ax_c.get_yticklabels()):
    if idx_t < 3:
        t.set_fontweight("bold")
        t.set_color(D_ITEMS[idx_t][2])
ax_c.tick_params(axis="y", length=0)
ax_c.spines["left"].set_visible(False)
ax_c.set_xlim(0.12, 0.86)
ax_c.set_xticks([0.25, 0.50, 0.75])
ax_c.set_ylim(-0.55, len(D_ITEMS) - 0.35)
ax_c.set_xlabel("Public Meta Death AUC (diamond: All  ·  triangle: WHO-CPS ≥ 6)", fontsize=9.6)
ax_c.set_title("Mortality: univariate vs. within-severe", loc="left", fontsize=11.3, fontweight="bold", pad=8)
fig.text(0.015, 0.435, "c", fontsize=18, fontweight="bold", va="top")


# ---------------- Panel d: HSPC Mortality Across Time Windows (First acute vs > day 7) ----------------
ax_d = fig.add_axes([0.675, 0.095, 0.215, 0.320])
COH_D = ["Stephenson 2021", "van der Wijst 2021", "COMBAT 2022", "Public Meta (3 cohorts)"]
yy_d = np.arange(len(COH_D))[::-1]
ax_d.axvline(0.50, color=PALETTE["gray_1"], ls=":", lw=1.3)

for yi, cname in zip(yy_d, COH_D):
    r_f = df_hw[(df_hw["cohort"] == cname) & (df_hw["window"] == "First acute sample")].iloc[0]
    r_l = df_hw[(df_hw["cohort"] == cname) & (df_hw["window"] == "Samples > day 7")].iloc[0]
    is_meta = "Meta" in cname
    c_f = PALETTE["gray_3"]
    c_l = PALETTE["red_strong"]

    ax_d.plot([r_f["AUC_all"], r_l["AUC_all"]], [yi + 0.16, yi - 0.16], color=PALETTE["neutral"], lw=1.3)
    ax_d.scatter(r_f["AUC_all"], yi + 0.16, s=62 if is_meta else 46, marker="D",
                 c="white" if r_f["p_all"] >= 0.05 else c_f, edgecolors=c_f, linewidths=1.5, zorder=3)
    ax_d.scatter(r_l["AUC_all"], yi - 0.16, s=66 if is_meta else 50, marker="s",
                 c="white" if r_l["p_all"] >= 0.05 else c_l, edgecolors=c_l, linewidths=1.5, zorder=4)

    ax_d.text(0.95, yi + 0.16, f"1st: {r_f['AUC_all']:.2f} (P={r_f['p_all']:.3f})", va="center", fontsize=8.0,
              color=c_f if r_f["p_all"] < 0.05 else PALETTE["gray_1"], fontweight="bold" if r_f["p_all"] < 0.05 else "normal")
    ax_d.text(0.95, yi - 0.16, f">7d: {r_l['AUC_all']:.2f} (P={r_l['p_all']:.3f})", va="center", fontsize=8.0,
              color=c_l if r_l["p_all"] < 0.05 else PALETTE["gray_1"], fontweight="bold" if r_l["p_all"] < 0.05 else "normal")

ax_d.set_yticks(yy_d)
ax_d.set_yticklabels(COH_D, fontsize=9.5)
ax_d.get_yticklabels()[-1].set_fontweight("bold")
ax_d.get_yticklabels()[-1].set_color(PALETTE["red_strong"])
ax_d.tick_params(axis="y", length=0)
ax_d.spines["left"].set_visible(False)
ax_d.set_xlim(0.42, 0.93)
ax_d.set_xticks([0.50, 0.70, 0.90])
ax_d.set_ylim(-0.55, len(COH_D) - 0.35)
ax_d.set_xlabel("HSPC Death AUC (diamond: 1st  ·  square: >7d)", fontsize=9.5)
ax_d.set_title("HSPC mortality across windows", loc="left", fontsize=11.3, fontweight="bold", pad=8)
fig.text(0.545, 0.435, "d", fontsize=18, fontweight="bold", va="top")

# Explicit Cohort Proxy Footnote at bottom of figure
fig.text(
    0.015,
    0.018,
    "*Cohort annotation proxies: B ASC dividing = Plasmablast (Stephenson), PB_Prolif (van der Wijst); T CD8 CM = CD8.EM (Stephenson), CD8.TCM+TCM.CCL5 (COMBAT); "
    "T CD8 CTL = CD8.TE (Stephenson), CD8.TEMRA+TEFF+TEM (COMBAT), T8_Mem (van der Wijst), CD8 memory T (Liu). All available cohorts per feature are pooled.",
    fontsize=8.3,
    color=PALETTE["gray_2"],
    va="bottom",
    ha="left",
)

out_fig = OUT_DIR / "PBISC_SuppFig_S3_Public_COVID_Replication_20261118"
saved = finalize_figure(fig, out_fig, formats=("png", "pdf", "svg"), dpi=300, close=True)
for ext in ("png", "pdf", "svg"):
    shutil.copy2(OUT_DIR / f"PBISC_SuppFig_S3_Public_COVID_Replication_20261118.{ext}",
                 OUT_DIR / f"PBISC_SuppFig_S3_Public_COVID_Replication_20261117.{ext}")
print("Saved Corrected External Replication Figure:", [str(p) for p in saved])
