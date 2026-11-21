"""Comprehensive forensic audit for the three user priorities:
  1) Exact patient & dataset metadata matching (Bacterial Sepsis + Imperial ESKD COVID-19),
     including sex_flag / bulk_sex_call checks, donor overlap, 1:1 visit mapping, and draw QC.
  2) Population NUMBER (raw cell count, log1p count, rank-INT count, AND clinical CBC absolute count)
     vs PERCENTAGE (% of PBMC, % of parent lineage) severity (same-day WHO GEE +wk) and
     mortality (donor first-draw AUC + adjWHO) correlation check across COVID-19 and Sepsis.
  3) Deep T-cell ISG scoring analysis across all 13 COVID-19 T-cell compartments (Shared-27 vs
     BayesPrism Z, same-day WHO +wk vs unadjusted vs Peak WHO, first-draw death + adjWHO + <=7d + >7d,
     longitudinal Time x Fatal collapse, T CD8 CTL abundance vs ISG independence, and leak-free 50x5 CV)
     PLUS T-cell ISG scoring in Bacterial Sepsis (PBISC 122 samples vs Real).
"""
import json
import os
import pathlib
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

ROOT = pathlib.Path(r"F:\PBISC")
D28 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "28_Verified_True11_Deconv_HSC_TISG_20261119" / "data"
D30 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "30_DrawLevel_Reanalysis_and_Patient_Draw_Table_20261120" / "tables"
OUT_ROOT = ROOT / "deliverables" / "nature_submission_fig4_6_20261120"
AUDIT_OUT = OUT_ROOT / "audit_focus3_20261121"
AUDIT_OUT.mkdir(parents=True, exist_ok=True)

CLL = {"C116", "C139"}
ISG19 = [
    "ISG15", "IFI6", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "XAF1", "STAT1",
    "IRF7", "EPSTI1", "HERC5", "SIGLEC1",
]


def read_h5ad_obs(path: pathlib.Path, cols=None) -> pd.DataFrame:
    out = {}
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        keys = cols if cols is not None else [k for k in obs.keys() if k != "_index"]
        for k in keys:
            if k not in obs:
                continue
            item = obs[k]
            if isinstance(item, h5py.Group) and "codes" in item and "categories" in item:
                cats = [c.decode("utf-8") if isinstance(c, bytes) else c for c in item["categories"][:]]
                codes = item["codes"][:]
                out[k] = pd.Categorical.from_codes(codes, categories=cats)
            else:
                vals = item[:]
                if vals.dtype.kind in ("S", "O"):
                    vals = [v.decode("utf-8") if isinstance(v, bytes) else v for v in vals]
                out[k] = vals
    return pd.DataFrame(out)


def rank_int(vals):
    s = pd.Series(vals, dtype=float)
    m = s.notna()
    res = pd.Series(np.nan, index=s.index)
    if m.sum() >= 3:
        r = stats.rankdata(s[m].values)
        res.loc[m] = stats.norm.ppf((r - 0.5) / len(r))
    return res.values


def auc_mw(y, f):
    y = np.asarray(y, float)
    f = np.asarray(f, int)
    m = np.isfinite(y)
    y, f = y[m], f[m]
    a, b = y[f == 1], y[f == 0]
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u / (len(a) * len(b))), float(p)


def gee_sameday_wk(df, ycol, who_col="temp_severity_ord", days_col="days", id_col="patient_id"):
    sub = df[[id_col, ycol, who_col, days_col]].dropna().copy()
    sub["wk"] = sub[days_col] / 7.0
    sub = sub.sort_values([id_col, days_col])
    if len(sub) < 10 or sub[ycol].std() == 0:
        return np.nan, np.nan, np.nan
    X = sm.add_constant(sub[[who_col, "wk"]])
    try:
        m = sm.GEE(sub[ycol], X, groups=sub[id_col], family=sm.families.Gaussian(),
                   cov_struct=sm.cov_struct.Exchangeable()).fit()
        b = float(m.params[who_col])
        p = float(m.pvalues[who_col])
        try:
            mbr = sm.GEE(sub[ycol], X, groups=sub[id_col], family=sm.families.Gaussian(),
                         cov_struct=sm.cov_struct.Exchangeable()).fit(cov_type="bias_reduced")
            pbr = float(mbr.pvalues[who_col])
        except Exception:
            pbr = np.nan
        return b, p, pbr
    except Exception:
        return np.nan, np.nan, np.nan


def bh_q(pvals):
    p = np.asarray(pvals, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranks = np.arange(1, len(pv) + 1)
    adj = np.minimum.accumulate((pv[order] * len(pv) / ranks)[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def audit_1_metadata_matching():
    print("=" * 80)
    print("AUDIT 1: PATIENT & DATASET METADATA MATCHING (SEPSIS + ESKD COVID-19)")
    print("=" * 80)
    summary = {}

    # --- 1A. Bacterial Sepsis (Ye et al., GSE279448 bulk / GSE279452 scRNA) ---
    geo = pd.read_csv(ROOT / "sepsis_nk_cd56bright_20261115" / "data" / "geo_characteristics.csv")
    sep_real = pd.read_csv(ROOT / "sepsis_nk_cd56bright_20261115" / "data" / "sepsis_real_sample_table.csv")
    sep_pb_obs = read_h5ad_obs(
        ROOT / "deliverables" / "annotated_h5ad_20261112" / "pbisc_sepsis_ye2026_draw1_K4000_ceiling_annotated.h5ad",
        cols=["sample_id", "donor_id", "gsm", "condition", "group", "infection_source", "age", "sex_geo",
              "bulk_sex_call", "sex_flag", "collection_site", "in_primary_scope", "support_donor_overlap"],
    )
    sep_pb = sep_pb_obs.drop_duplicates("sample_id").reset_index(drop=True)
    # Match geo to sep_pb and sep_real
    geo_sc = geo[geo["title"].str.startswith("scRNA ")].copy()
    geo_bk = geo[~geo["title"].str.startswith("scRNA ")].copy()
    sep_real_ids = set(sep_real["sample"])
    sep_pb_ids = set(sep_pb["sample_id"])
    sep_overlap = sep_real_ids & sep_pb_ids

    # Build unified sepsis donor metadata table (Real 133 + PBISC 122 = 255 rows; 108 adult Real + 121 adult PBISC + 1 ped PBISC)
    sep_pb["age_num"] = pd.to_numeric(sep_pb["age"], errors="coerce")
    sep_pb["age_group"] = np.where(sep_pb["condition"].astype(str).str.startswith("ped"), "paediatric", "adult")
    sep_pb["dataset"] = "PBISC_bulk_only"
    sep_real["dataset"] = "Real_scRNA"

    sep_real_ad = sep_real[sep_real["age_group"] == "adult"]
    sep_real_pd = sep_real[sep_real["age_group"] == "paediatric"]
    sep_pb_ad = sep_pb[sep_pb["age_group"] == "adult"]
    sep_pb_pd = sep_pb[sep_pb["age_group"] == "paediatric"]

    sex_flags_sep = sep_pb[sep_pb["sex_flag"].astype(str).isin(["mismatch", "discordant", "True", "1"]) | (
        sep_pb["sex_geo"].notna() & sep_pb["bulk_sex_call"].notna() & (
            sep_pb["sex_geo"].astype(str).str[0].str.upper() != sep_pb["bulk_sex_call"].astype(str).str[0].str.upper()
        )
    )]
    print(f"[Sepsis] GEO rows: {len(geo)} (scRNA={len(geo_sc)}, non-scRNA={len(geo_bk)})")
    print(f"[Sepsis] Real scRNA samples: {len(sep_real)} (adult={len(sep_real_ad)}, ped={len(sep_real_pd)}) | Real cell total med (IQR): "
          f"{sep_real_ad['total'].median():.0f} ({sep_real_ad['total'].quantile(0.25):.0f}-{sep_real_ad['total'].quantile(0.75):.0f})")
    print(f"[Sepsis] PBISC bulk-only samples: {len(sep_pb)} (adult={len(sep_pb_ad)}, ped={len(sep_pb_pd)}: {sep_pb_pd['sample_id'].tolist()})")
    print(f"[Sepsis] Donor overlap Real ∩ PBISC bulk-only: {len(sep_overlap)} (expected 0; support_donor_overlap sum={int(sep_pb['support_donor_overlap'].astype(str).isin(['True','1']).sum())})")
    print(f"[Sepsis] Real adult age med (IQR): {sep_real_ad['age'].median():.0f} ({sep_real_ad['age'].quantile(0.25):.0f}-{sep_real_ad['age'].quantile(0.75):.0f}), "
          f"Male: {(sep_real_ad['gender']=='Male').sum()}/{len(sep_real_ad)} ({100*(sep_real_ad['gender']=='Male').mean():.1f}%)")
    print(f"[Sepsis] PBISC adult age med (IQR): {sep_pb_ad['age_num'].median():.0f} ({sep_pb_ad['age_num'].quantile(0.25):.0f}-{sep_pb_ad['age_num'].quantile(0.75):.0f}), "
          f"Male(sex_geo): {(sep_pb_ad['sex_geo'].astype(str)=='Male').sum()}/{len(sep_pb_ad)} ({100*(sep_pb_ad['sex_geo'].astype(str)=='Male').mean():.1f}%)")
    print(f"[Sepsis] PBISC sex_flag value_counts: {sep_pb['sex_flag'].astype(str).value_counts().to_dict()}")
    if len(sex_flags_sep):
        print("  Sepsis sex_flag rows:", sex_flags_sep[["sample_id", "condition", "age", "sex_geo", "bulk_sex_call", "sex_flag"]].to_dict("records"))
    print(f"[Sepsis] PBISC condition counts: {sep_pb['condition'].astype(str).value_counts().to_dict()}")
    print(f"[Sepsis] Real adult group counts: {sep_real_ad['group'].value_counts().to_dict()}")
    print(f"[Sepsis] Real adult infection_site: {sep_real_ad['infection_site'].fillna('none').value_counts().to_dict()}")
    print(f"[Sepsis] PBISC adult infection_source: {sep_pb_ad['infection_source'].astype(str).value_counts().to_dict()}")

    # Save Sepsis complete matched metadata table
    sep_pb_clean = sep_pb[["sample_id", "donor_id", "gsm", "condition", "group", "age_group", "age_num", "sex_geo", "bulk_sex_call", "sex_flag", "infection_source", "collection_site", "in_primary_scope", "support_donor_overlap"]].copy()
    sep_pb_clean["dataset"] = "PBISC_bulk_only"
    sep_pb_clean["n_cells"] = 4000
    sep_pb_clean.rename(columns={"age_num": "age", "sex_geo": "sex"}, inplace=True)

    sep_real_clean = sep_real[["sample", "gsm", "group", "age_group", "age", "gender", "infection_site", "total", "shared_with_real"]].copy()
    sep_real_clean.rename(columns={"sample": "sample_id", "gender": "sex", "infection_site": "infection_source", "total": "n_cells"}, inplace=True)
    sep_real_clean["donor_id"] = sep_real_clean["sample_id"]
    sep_real_clean["dataset"] = "Real_scRNA"
    sep_real_clean["in_primary_scope"] = sep_real_clean["age_group"] == "adult"
    sep_all = pd.concat([sep_real_clean, sep_pb_clean], ignore_index=True, sort=False)
    sep_all.to_csv(AUDIT_OUT / "A1_sepsis_255_samples_matched_metadata.csv", index=False)

    # --- 1B. Imperial ESKD COVID-19 (107 donors / 362 samples) ---
    s2 = pd.read_csv(D28 / "table_s2_sample_1to1_mapping_362.csv")
    s3 = pd.read_csv(D28 / "table_s3_all_107_donors_overview.csv")
    meta = pd.read_csv(D28 / "eskd_sample_metadata_verified.csv")
    cov_pb_obs = read_h5ad_obs(
        ROOT / "deliverables" / "annotated_h5ad_20261112" / "pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad",
        cols=["sample_id", "donor_id", "condition", "wave", "waves", "WHO_severity", "WHO_temp_severity",
              "fatal_disease", "calc_age", "sex", "bulk_sex_call", "sex_flag", "time_from_first_positive_swab",
              "in_primary_scope", "support_donor_overlap"],
    )
    cov_pb = cov_pb_obs.drop_duplicates("sample_id").reset_index(drop=True)
    sex_flags_cov = cov_pb[cov_pb["sex_flag"].astype(str).isin(["True", "1", "mismatch", "discordant"])]
    print(f"\n[COVID-19] Total donors in Table S3: {len(s3)} | cohort_role counts: {s3['cohort_role'].value_counts().to_dict()}")
    print(f"[COVID-19] 1:1 mapping rows in Table S2: {len(s2)} | mapping_type counts:")
    for k, v in s2["mapping_type"].value_counts().items():
        print(f"    {k}: {v}")
    print(f"[COVID-19] PBISC h5ad samples: {len(cov_pb)} | sex_flag counts: {cov_pb['sex_flag'].astype(str).value_counts().to_dict()}")
    if len(sex_flags_cov):
        print("  COVID PBISC non-match sex_flag:", sex_flags_cov[["sample_id", "donor_id", "condition", "sex", "bulk_sex_call", "sex_flag"]].drop_duplicates("donor_id").to_dict("records"))

    # Check Shared-27 Real vs PBISC metadata concordance (age, sex, peak WHO, fatal)
    r_pt = meta[meta["source"] == "Real"].drop_duplicates("patient_id").set_index("patient_id")
    p_pt = meta[meta["source"] == "PBISC"].drop_duplicates("patient_id").set_index("patient_id")
    shared_ids = sorted(set(r_pt.index) & set(p_pt.index))
    mismatches = []
    for pid in shared_ids:
        for col in ["age", "sex", "WHO_severity", "fatal_pub11"]:
            rv, pv = r_pt.loc[pid, col], p_pt.loc[pid, col]
            if pd.notna(rv) and pd.notna(pv) and str(rv) != str(pv):
                mismatches.append((pid, col, rv, pv))
    print(f"[COVID-19] Shared donors in meta (Real ∩ PBISC): {len(shared_ids)} | demographic/WHO/fatal mismatches: {len(mismatches)}")
    if mismatches:
        print("  Mismatches:", mismatches[:10])

    # Check n_cells distribution in primary 217 draws (Real vs PBISC new)
    p217 = meta[(meta["condition"] == "positive") & meta["in_primary_scope"].astype(bool) &
                meta["qc_pass_min200cells"].astype(bool) & (~meta["patient_id"].isin(CLL))].copy()
    real_pids = set(p217[p217["source"] == "Real"]["patient_id"])
    p217 = p217[(p217["source"] == "Real") | (~p217["patient_id"].isin(real_pids))].copy()
    r217 = p217[p217["source"] == "Real"]
    b217 = p217[p217["source"] == "PBISC"]
    print(f"[COVID-19] Primary 217 draws verified: Total={len(p217)} (Real={len(r217)} across {r217['patient_id'].nunique()} pts, "
          f"PBISC new={len(b217)} across {b217['patient_id'].nunique()} pts)")
    print(f"[COVID-19] Real n_cells/draw med (IQR) [min-max]: {r217['n_cells'].median():.0f} "
          f"({r217['n_cells'].quantile(0.25):.0f}-{r217['n_cells'].quantile(0.75):.0f}) [{r217['n_cells'].min():.0f}-{r217['n_cells'].max():.0f}]")
    print(f"[COVID-19] PBISC n_cells/draw: constant {b217['n_cells'].unique().tolist()}")

    # CRITICAL CHECK: Does Real total n_cells per draw correlate with same-day WHO or death?!
    rho_nc_sev, p_nc_sev = stats.spearmanr(r217["n_cells"], r217["temp_severity_ord"])
    r_fd = r217.sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).first()
    auc_nc_dth, p_nc_dth = auc_mw(r_fd["n_cells"], r_fd["fatal_pub11"])
    print(f"[CRITICAL CONFOUNDER CHECK] Real total n_cells vs Same-day WHO: Spearman rho={rho_nc_sev:+.3f}, P={p_nc_sev:.4g}")
    print(f"[CRITICAL CONFOUNDER CHECK] Real total n_cells (1st draw) vs Death: AUC={auc_nc_dth:.3f}, P={p_nc_dth:.4g}")

    summary["sepsis"] = {
        "real_total": int(len(sep_real)), "real_adult": int(len(sep_real_ad)), "real_ped": int(len(sep_real_pd)),
        "pbisc_total": int(len(sep_pb)), "pbisc_adult": int(len(sep_pb_ad)), "pbisc_ped": int(len(sep_pb_pd)),
        "pbisc_ped_ids": sep_pb_pd["sample_id"].tolist(),
        "real_pbisc_overlap": int(len(sep_overlap)),
        "real_adult_age_med": float(sep_real_ad["age"].median()),
        "real_adult_age_iqr": f"{sep_real_ad['age'].quantile(0.25):.0f}-{sep_real_ad['age'].quantile(0.75):.0f}",
        "real_adult_male_n": int((sep_real_ad["gender"] == "Male").sum()),
        "pbisc_adult_age_med": float(sep_pb_ad["age_num"].median()),
        "pbisc_adult_age_iqr": f"{sep_pb_ad['age_num'].quantile(0.25):.0f}-{sep_pb_ad['age_num'].quantile(0.75):.0f}",
        "pbisc_adult_male_n": int((sep_pb_ad["sex_geo"].astype(str) == "Male").sum()),
        "pbisc_sex_flag_counts": {str(k): int(v) for k, v in sep_pb["sex_flag"].astype(str).value_counts().items()},
    }
    summary["covid"] = {
        "total_donors": int(len(s3)), "total_mapping_rows": int(len(s2)),
        "shared_donors_meta": int(len(shared_ids)), "shared_metadata_mismatches": len(mismatches),
        "pbisc_sex_flag_counts": {str(k): int(v) for k, v in cov_pb["sex_flag"].astype(str).value_counts().items()},
        "real_ncells_med": float(r217["n_cells"].median()),
        "real_ncells_iqr": f"{r217['n_cells'].quantile(0.25):.0f}-{r217['n_cells'].quantile(0.75):.0f}",
        "real_ncells_min": int(r217["n_cells"].min()), "real_ncells_max": int(r217["n_cells"].max()),
        "real_ncells_vs_sdWHO_rho": float(rho_nc_sev), "real_ncells_vs_sdWHO_p": float(p_nc_sev),
        "real_ncells_1st_vs_death_AUC": float(auc_nc_dth), "real_ncells_1st_vs_death_p": float(p_nc_dth),
    }
    return summary, p217


def audit_2_count_vs_percentage(p217_meta):
    print("\n" + "=" * 80)
    print("AUDIT 2: POPULATION NUMBER (COUNT & CLINICAL CBC ABS) VS PERCENTAGE (% PBMC & % PARENT)")
    print("=" * 80)

    # Load primary 217-draw raw percentages + cell counts directly from real_obs_full & cell_labels_covid
    p217_z = pd.read_csv(OUT_ROOT / "data" / "primary_217draws_rankINT_43subsets.csv")
    feat43 = [c.replace("__z217", "") for c in p217_z.columns if c.endswith("__z217")]

    # Merge exact total n_cells & CBC labs (temp_WCC, temp_lymph, temp_mono, temp_neut) onto p217_z
    lab_cols = ["source", "sample_id", "n_cells", "temp_WCC", "temp_lymph", "temp_mono", "temp_neut", "temp_CRP"]
    meta_sub = p217_meta[lab_cols].drop_duplicates(["source", "sample_id"]).copy()
    meta_sub["arm"] = np.where(meta_sub["source"] == "Real", "Real", "PBISC new")
    drop_dup = [c for c in ["n_cells", "temp_WCC", "temp_lymph", "temp_mono", "temp_neut", "temp_CRP"] if c in p217_z.columns]
    df = p217_z.drop(columns=drop_dup).merge(meta_sub.drop(columns=["source"]), on=["arm", "sample_id"], how="left")
    assert df["n_cells"].notna().all(), "Missing n_cells on primary 217 draws!"

    # Exact raw cell count per subset = round(pct / 100 * n_cells)
    # Verify exact integer accuracy: pct was computed as 100 * count / n_cells
    counts_mat = {}
    for f in feat43:
        counts_mat[f] = np.round(df[f].values / 100.0 * df["n_cells"].values).astype(int)
    counts_df = pd.DataFrame(counts_mat, index=df.index)

    # Map each subset to parent lineage for % of parent lineage AND clinical CBC scaling
    def parent_of(f):
        if f.startswith("T ") or f in ("MAIT", "NKT"):
            return "T"
        if f.startswith("B_"):
            return "B"
        if f.startswith("NK"):
            return "NK"
        if "mono" in f.lower() or f.startswith("DC") or f in ("pDC", "ASDC"):
            return "Myeloid"
        return "Other"

    # Compute parent sums per draw directly from raw percentages (across all 43 subsets)
    parent_groups = {p: [f for f in feat43 if parent_of(f) == p] for p in ["T", "B", "NK", "Myeloid"]}
    # Total lymphocyte % = T + B + NK
    lymph_feats = parent_groups["T"] + parent_groups["B"] + parent_groups["NK"]
    df["lymph_pct_sum"] = df[lymph_feats].sum(axis=1)
    for p, flist in parent_groups.items():
        df[f"parent_{p}_pct"] = df[flist].sum(axis=1)

    # First-draw slice (N=66, 11 deaths)
    fd_idx = df.sort_values(["patient_id", "days"]).groupby("patient_id").head(1).index
    df_fd = df.loc[fd_idx].copy()
    counts_fd = counts_df.loc[fd_idx].copy()

    rows = []
    for f in feat43:
        par = parent_of(f)
        # 1. % PBMC (within-source rank-INT — Primary)
        z_pct_217 = df[f + "__z217"].values
        # Raw % PBMC (untransformed)
        raw_pct_217 = df[f].values
        # 2. Raw cell count (population number): raw count & within-source rank-INT(count)
        raw_cnt_217 = counts_df[f].values.astype(float)
        log_cnt_217 = np.log1p(raw_cnt_217)
        z_cnt_217 = df.assign(_v=raw_cnt_217).groupby("arm")["_v"].transform(rank_int).values
        # 3. % of parent lineage (if T/B/NK/Myeloid)
        if par in parent_groups:
            pct_par_217 = 100.0 * df[f].values / np.maximum(df[f"parent_{par}_pct"].values, 1e-6)
            z_par_217 = df.assign(_v=pct_par_217).groupby("arm")["_v"].transform(rank_int).values
        else:
            pct_par_217 = raw_pct_217
            z_par_217 = z_pct_217

        # 4. Clinical absolute blood count in Real (where CBC exists)
        # For T/B/NK: ALC (10^9/L) * (subset% / lymph_pct_sum)
        # For Monocytes: temp_mono * (subset% / parent_Myeloid)
        # For HSPC/Other: temp_WCC * (subset% / 100)
        if par in ("T", "B", "NK"):
            abs_clin_217 = df["temp_lymph"].values * (df[f].values / np.maximum(df["lymph_pct_sum"].values, 1e-6))
        elif "mono" in f.lower():
            abs_clin_217 = df["temp_mono"].values * (df[f].values / np.maximum(df["parent_Myeloid_pct"].values, 1e-6))
        else:
            abs_clin_217 = df["temp_WCC"].values * (df[f].values / 100.0)
        z_abs_real = rank_int(np.where(df["arm"] == "Real", abs_clin_217, np.nan))

        # Run Same-Day WHO GEE (+wk) on Combined (217) and Real (138) and PBISC new (79)
        tmp = df[["patient_id", "arm", "days", "temp_severity_ord", "fatal_pub11"]].copy()
        tmp["z_pct"] = z_pct_217
        tmp["raw_pct"] = raw_pct_217
        tmp["z_cnt"] = z_cnt_217
        tmp["log_cnt"] = log_cnt_217
        tmp["z_par"] = z_par_217
        tmp["z_abs_real"] = z_abs_real

        b_pct_c, p_pct_c, _ = gee_sameday_wk(tmp, "z_pct")
        b_cnt_c, p_cnt_c, _ = gee_sameday_wk(tmp, "z_cnt")
        b_lcn_c, p_lcn_c, _ = gee_sameday_wk(tmp, "log_cnt")
        b_par_c, p_par_c, _ = gee_sameday_wk(tmp, "z_par")
        b_pct_r, p_pct_r, _ = gee_sameday_wk(tmp[tmp["arm"] == "Real"], "z_pct")
        b_cnt_r, p_cnt_r, _ = gee_sameday_wk(tmp[tmp["arm"] == "Real"], "z_cnt")
        b_abs_r, p_abs_r, _ = gee_sameday_wk(tmp[tmp["arm"] == "Real"], "z_abs_real")
        b_pct_n, p_pct_n, _ = gee_sameday_wk(tmp[tmp["arm"] == "PBISC new"], "z_pct")
        b_cnt_n, p_cnt_n, _ = gee_sameday_wk(tmp[tmp["arm"] == "PBISC new"], "z_cnt")

        # Donor first-draw mortality (N=66, 11 deaths)
        tmp_fd = tmp.loc[fd_idx].copy()
        tmp_fd["z_pct_fd"] = df_fd.groupby("arm")[f].transform(rank_int).values
        tmp_fd["raw_cnt_fd"] = counts_fd[f].values.astype(float)
        tmp_fd["z_cnt_fd"] = tmp_fd.groupby("arm")["raw_cnt_fd"].transform(rank_int).values
        tmp_fd["raw_pct_fd"] = df_fd[f].values
        tmp_fd["pct_par_fd"] = pct_par_217[df.index.isin(fd_idx)]
        tmp_fd["z_par_fd"] = tmp_fd.groupby("arm")["pct_par_fd"].transform(rank_int).values
        tmp_fd["abs_clin_fd"] = abs_clin_217[df.index.isin(fd_idx)]

        auc_pct_c, pd_pct_c = auc_mw(tmp_fd["z_pct_fd"], tmp_fd["fatal_pub11"])
        auc_cnt_c, pd_cnt_c = auc_mw(tmp_fd["z_cnt_fd"], tmp_fd["fatal_pub11"])
        auc_rcnt_c, pd_rcnt_c = auc_mw(tmp_fd["raw_cnt_fd"], tmp_fd["fatal_pub11"])
        auc_rpct_c, pd_rpct_c = auc_mw(tmp_fd["raw_pct_fd"], tmp_fd["fatal_pub11"])
        auc_par_c, pd_par_c = auc_mw(tmp_fd["z_par_fd"], tmp_fd["fatal_pub11"])

        r_sub = tmp_fd[tmp_fd["arm"] == "Real"]
        n_sub = tmp_fd[tmp_fd["arm"] == "PBISC new"]
        auc_pct_r, pd_pct_r = auc_mw(r_sub["raw_pct_fd"], r_sub["fatal_pub11"])
        auc_cnt_r, pd_cnt_r = auc_mw(r_sub["raw_cnt_fd"], r_sub["fatal_pub11"])
        auc_abs_r, pd_abs_r = auc_mw(r_sub["abs_clin_fd"], r_sub["fatal_pub11"])
        auc_pct_n, pd_pct_n = auc_mw(n_sub["raw_pct_fd"], n_sub["fatal_pub11"])
        auc_cnt_n, pd_cnt_n = auc_mw(n_sub["raw_cnt_fd"], n_sub["fatal_pub11"])

        # Same-day WHO adjusted first-draw death for z_pct_fd and z_cnt_fd
        def adj_who_auc(sub_f, col):
            d0 = sub_f[[col, "temp_severity_ord", "fatal_pub11"]].dropna()
            res = sm.OLS(d0[col], sm.add_constant(d0[["temp_severity_ord"]])).fit().resid
            return auc_mw(res, d0["fatal_pub11"])

        adj_pct_auc, adj_pct_p = adj_who_auc(tmp_fd, "z_pct_fd")
        adj_cnt_auc, adj_cnt_p = adj_who_auc(tmp_fd, "z_cnt_fd")

        rows.append({
            "feature": f, "parent": par,
            "real_med_cnt": float(np.median(counts_df.loc[df["arm"] == "Real", f])),
            "pbisc_med_cnt": float(np.median(counts_df.loc[df["arm"] == "PBISC new", f])),
            "real_med_pct": float(np.median(df.loc[df["arm"] == "Real", f])),
            "pbisc_med_pct": float(np.median(df.loc[df["arm"] == "PBISC new", f])),
            # Same-day WHO GEE +wk
            "sd_beta_pct_INT": b_pct_c, "sd_p_pct_INT": p_pct_c,
            "sd_beta_cnt_INT": b_cnt_c, "sd_p_cnt_INT": p_cnt_c,
            "sd_beta_logcnt": b_lcn_c, "sd_p_logcnt": p_lcn_c,
            "sd_beta_par_INT": b_par_c, "sd_p_par_INT": p_par_c,
            "sd_beta_pct_Real": b_pct_r, "sd_p_pct_Real": p_pct_r,
            "sd_beta_cnt_Real": b_cnt_r, "sd_p_cnt_Real": p_cnt_r,
            "sd_beta_absCBC_Real": b_abs_r, "sd_p_absCBC_Real": p_abs_r,
            "sd_beta_pct_PBISC": b_pct_n, "sd_p_pct_PBISC": p_pct_n,
            "sd_beta_cnt_PBISC": b_cnt_n, "sd_p_cnt_PBISC": p_cnt_n,
            # First-draw mortality
            "fd_AUC_pct_INT": auc_pct_c, "fd_p_pct_INT": pd_pct_c, "fd_adjWHO_p_pct": adj_pct_p,
            "fd_AUC_cnt_INT": auc_cnt_c, "fd_p_cnt_INT": pd_cnt_c, "fd_adjWHO_p_cnt": adj_cnt_p,
            "fd_AUC_rawcnt_pooled": auc_rcnt_c, "fd_p_rawcnt_pooled": pd_rcnt_c,
            "fd_AUC_rawpct_pooled": auc_rpct_c, "fd_p_rawpct_pooled": pd_rpct_c,
            "fd_AUC_par_INT": auc_par_c, "fd_p_par_INT": pd_par_c,
            "fd_AUC_pct_Real": auc_pct_r, "fd_p_pct_Real": pd_pct_r,
            "fd_AUC_cnt_Real": auc_cnt_r, "fd_p_cnt_Real": pd_cnt_r,
            "fd_AUC_absCBC_Real": auc_abs_r, "fd_p_absCBC_Real": pd_abs_r,
            "fd_n_absCBC_Real": int(r_sub["abs_clin_fd"].notna().sum()),
            "fd_AUC_pct_PBISC": auc_pct_n, "fd_p_pct_PBISC": pd_pct_n,
            "fd_AUC_cnt_PBISC": auc_cnt_n, "fd_p_cnt_PBISC": pd_cnt_n,
        })

    res_df = pd.DataFrame(rows)
    res_df["sd_q43_pct_INT"] = bh_q(res_df["sd_p_pct_INT"])
    res_df["sd_q43_cnt_INT"] = bh_q(res_df["sd_p_cnt_INT"])
    res_df["sd_q43_par_INT"] = bh_q(res_df["sd_p_par_INT"])
    res_df["fd_q43_pct_INT"] = bh_q(res_df["fd_p_pct_INT"])
    res_df["fd_q43_cnt_INT"] = bh_q(res_df["fd_p_cnt_INT"])
    res_df.to_csv(AUDIT_OUT / "A2_covid_43subsets_count_vs_pct_vs_parent_vs_CBC.csv", index=False)

    # Print key comparison across all 43 subsets & focus subsets
    r_sev_cnt_pct, _ = stats.pearsonr(res_df["sd_beta_pct_INT"], res_df["sd_beta_cnt_INT"])
    r_dth_cnt_pct, _ = stats.pearsonr(res_df["fd_AUC_pct_INT"], res_df["fd_AUC_cnt_INT"])
    print(f"[COVID-19 43 subsets] Correlation (% PBMC INT vs Raw Count INT):")
    print(f"  Same-day WHO GEE +wk beta: r = {r_sev_cnt_pct:.4f} | FDR q43<0.05: pct={int((res_df['sd_q43_pct_INT']<0.05).sum())} vs cnt={int((res_df['sd_q43_cnt_INT']<0.05).sum())} vs parent={int((res_df['sd_q43_par_INT']<0.05).sum())}")
    print(f"  First-draw Death AUC:      r = {r_dth_cnt_pct:.4f} | FDR q43<0.10: pct={int((res_df['fd_q43_pct_INT']<0.10).sum())} vs cnt={int((res_df['fd_q43_cnt_INT']<0.10).sum())}")
    print("\n  Key subsets — Same-Day WHO GEE (+wk) β [% PBMC | Raw Count | % Parent | Real CBC abs]:")
    focus = ["B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "HSPC", "NK CD56 bright", "T CD8 CM", "T CD8 CTL", "NK", "CD14mono_IFN", "CD14mono_dex"]
    sub_f = res_df.set_index("feature").loc[focus]
    for f, r in sub_f.iterrows():
        print(f"    {f:16s} | SD β: pct={r.sd_beta_pct_INT:+.3f}(q={r.sd_q43_pct_INT:.1e}) cnt={r.sd_beta_cnt_INT:+.3f}(q={r.sd_q43_cnt_INT:.1e}) "
              f"par={r.sd_beta_par_INT:+.3f}(q={r.sd_q43_par_INT:.1e}) RealCBC={r.sd_beta_absCBC_Real:+.3f}(P={r.sd_p_absCBC_Real:.3g})")
    print("\n  Key subsets — Donor First-Draw Mortality AUC [% PBMC INT | Raw Count INT | Raw Count Pooled | % Parent | Real CBC abs]:")
    for f, r in sub_f.iterrows():
        print(f"    {f:16s} | 1st AUC: pct={r.fd_AUC_pct_INT:.3f}(P={r.fd_p_pct_INT:.4f},aP={r.fd_adjWHO_p_pct:.3f}) "
              f"cntINT={r.fd_AUC_cnt_INT:.3f}(P={r.fd_p_cnt_INT:.4f},aP={r.fd_adjWHO_p_cnt:.3f}) "
              f"rawCnt={r.fd_AUC_rawcnt_pooled:.3f} par={r.fd_AUC_par_INT:.3f} RealCBC={r.fd_AUC_absCBC_Real:.3f}(P={r.fd_p_absCBC_Real:.3g},n={int(r.fd_n_absCBC_Real)})")

    # Also check Sepsis count vs % across HC (0) -> ICU (1) -> Sepsis (2) in Real adult (108) & PBISC adult (121)
    rc = pd.read_csv(ROOT / "sepsis_nk_cd56bright_20261115" / "data" / "real_sample_celltype_counts.csv")
    pc = pd.read_csv(ROOT / "sepsis_nk_cd56bright_20261115" / "data" / "pbisc_sample_celltype_counts.csv")
    rc = rc[rc["age_group"] == "adult"].copy()
    pc = pc[pc["in_primary_scope"].astype(bool)].copy()
    ord_map_r = {"healthy control": 0, "ICU control": 1, "sepsis": 2}
    ord_map_p = {"adult_HC": 0, "adult_ICU_ctrl": 1, "adult_sepsis": 2}
    rc["ord"] = rc["group"].map(ord_map_r)
    pc["ord"] = pc["condition"].map(ord_map_p)
    r_tot = rc.groupby("sample")["n"].sum()
    p_tot = pc.groupby("sample_id")["n"].sum()

    r_piv = rc.pivot_table(index="sample", columns="celltype_ye2026", values="n", aggfunc="sum", fill_value=0)
    p_piv = pc.pivot_table(index="sample_id", columns="celltype_fine", values="n", aggfunc="sum", fill_value=0)
    r_ord = rc.drop_duplicates("sample").set_index("sample")["ord"].reindex(r_piv.index)
    p_ord = pc.drop_duplicates("sample_id").set_index("sample_id")["ord"].reindex(p_piv.index)

    sep_rows = []
    common_ct = sorted(set(r_piv.columns) & set(p_piv.columns))
    for ct in common_ct:
        rcnt = r_piv[ct].values
        rpct = 100.0 * rcnt / r_tot.reindex(r_piv.index).values
        pcnt = p_piv[ct].values
        ppct = 100.0 * pcnt / p_tot.reindex(p_piv.index).values
        rr_c, rp_c = stats.spearmanr(rcnt, r_ord.values)
        rr_p, rp_p = stats.spearmanr(rpct, r_ord.values)
        pr_c, pp_c = stats.spearmanr(pcnt, p_ord.values)
        pr_p, pp_p = stats.spearmanr(ppct, p_ord.values)
        sep_rows.append({
            "celltype": ct,
            "real_rho_cnt": rr_c, "real_p_cnt": rp_c, "real_rho_pct": rr_p, "real_p_pct": rp_p,
            "pbisc_rho_cnt": pr_c, "pbisc_p_cnt": pp_c, "pbisc_rho_pct": pr_p, "pbisc_p_pct": pp_p,
        })
    sep_cnt_df = pd.DataFrame(sep_rows)
    sep_cnt_df.to_csv(AUDIT_OUT / "A2b_sepsis_count_vs_pct_HC_ICU_Sepsis.csv", index=False)
    print(f"\n[Sepsis HC->ICU->Sepsis] Common fine cell types={len(common_ct)} | "
          f"Real count vs % rho corr: r={stats.pearsonr(sep_cnt_df['real_rho_cnt'].fillna(0), sep_cnt_df['real_rho_pct'].fillna(0))[0]:.4f} | "
          f"PBISC count vs % rho corr: r={stats.pearsonr(sep_cnt_df['pbisc_rho_cnt'].fillna(0), sep_cnt_df['pbisc_rho_pct'].fillna(0))[0]:.4f}")
    for ct in ["NK_c01-GZMK", "B_c07-MZB1", "Macrophage", "Mono_c06-CD14-MKI67", "pDC", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67"]:
        if ct in sep_cnt_df["celltype"].values:
            r = sep_cnt_df[sep_cnt_df["celltype"] == ct].iloc[0]
            print(f"    {ct:22s} | Real: pct ρ={r.real_rho_pct:+.3f} cnt ρ={r.real_rho_cnt:+.3f} | PBISC: pct/cnt ρ={r.pbisc_rho_pct:+.3f}")

    return res_df, sep_cnt_df


def audit_3_tcell_isg(p217_meta):
    print("\n" + "=" * 80)
    print("AUDIT 3: DEEP T-CELL ISG SCORING ANALYSIS (ESKD COVID-19 + SEPSIS)")
    print("=" * 80)

    # 3A. COVID-19 T-cell ISG (all 13 subsets, primary + time windows + independence from T CD8 CTL % + leak-free CV)
    t6 = pd.read_csv(OUT_ROOT / "draw_sameday_sev_firstdraw_death_20261120" / "tables" / "T6_Tcell_ISG_13subsets_new_framework.csv")
    t8b = pd.read_csv(D30 / "T8b_Tcell_ISG_13subsets_time_windows.csv")
    tisg13 = pd.read_csv(D28 / "fig8_tcell_isg_13subsets_true11.csv").set_index("subset")
    tisg_gee = pd.read_csv(D28 / "suppS5_tcell_isg_gee_true11.csv")
    cd8_isg_217 = pd.read_csv(D28 / "fig8_cd8ctl_isg_longitudinal_draws.csv")
    p66_fd = pd.read_csv(OUT_ROOT / "data" / "primary_66firstdraws_rankINT_43subsets.csv")
    p217_z = pd.read_csv(OUT_ROOT / "data" / "primary_217draws_rankINT_43subsets.csv")

    print("[COVID-19 T-cell ISG 13 subsets — Combined Primary Summary]:")
    t6_c = t6[t6["set"] == "Combined"].set_index("subset")
    for s in ["All T cells", "All CD8 T", "All CD4 T", "T CD8 CTL", "T CD4 CTL", "T CD8 CM", "T CD8 Naive", "T CD4 Helper", "T CD4 Naive", "T reg", "T g/d", "MAIT"]:
        if s not in t6_c.index:
            continue
        r = t6_c.loc[s]
        sh = tisg13.loc[s, "sh_rho"] if s in tisg13.index else np.nan
        g_row = tisg_gee[(tisg_gee["subset"] == s) & (tisg_gee["cohort"] == "Combined")]
        f_wk = f"{g_row['fat_x_wk_beta'].iloc[0]:+.2f}(P={g_row['fat_x_wk_p'].iloc[0]:.1e})" if len(g_row) else "NA"
        print(f"  {s:14s} | Sh27 ρ={sh:+.2f} | PM Peak ρ={r.pm_peakWHO_rho:+.2f}(q={r.pm_q13_sev:.3f}) | "
              f"SD-unadj β={r.sd_gee_unadj_beta:+.2f}(P={r.sd_gee_unadj_p:.3f}) | SD+wk β={r.sd_gee_wk_beta:+.2f}(P={r.sd_gee_wk_p:.2f}) | "
              f"1st AUC={r.fd_AUC:.3f}(P={r.fd_p:.3f},aP={r.fd_adjWHO_p:.3f},≤7d={r.fd_le7_AUC:.3f}) | Time×Fatal={f_wk}")

    # Check independence between T CD8 CTL abundance (%) and T CD8 CTL ISG at first draw & across 217 draws
    cd8_m = cd8_isg_217.copy()
    cd8_m["arm"] = np.where(cd8_m["source"] == "Real", "Real", "PBISC new")
    m217 = p217_z.merge(cd8_m[["arm", "sample_id", "isg", "isg_z"]], on=["arm", "sample_id"], how="inner")
    r_abund_isg_217, p_abund_isg_217 = stats.spearmanr(m217["T CD8 CTL__z217"], m217["isg_z"])
    m66 = m217.sort_values(["patient_id", "days"]).groupby("patient_id", as_index=False).first()
    m66["ctl_zfd"] = m66.groupby("arm")["T CD8 CTL"].transform(rank_int)
    m66["nk_zfd"] = m66.groupby("arm")["NK"].transform(rank_int)
    m66["isg_zfd"] = m66.groupby("arm")["isg"].transform(rank_int)
    r_abund_isg_fd, p_abund_isg_fd = stats.spearmanr(m66["ctl_zfd"], m66["isg_zfd"])
    print(f"\n[COVID-19 Independence Check] T CD8 CTL abundance vs T CD8 CTL ISG:")
    print(f"  All 213 draws: Spearman rho = {r_abund_isg_217:+.3f}, P = {p_abund_isg_217:.4g}")
    print(f"  First draw (N={len(m66)}): Spearman rho = {r_abund_isg_fd:+.3f}, P = {p_abund_isg_fd:.4g}")
    # Composite first-draw score: -T CD8 CTL + T CD8 CTL ISG, and -(CTL+NK) + ISG
    a_ctl, p_ctl = auc_mw(-m66["ctl_zfd"], m66["fatal_pub11"])
    a_isg, p_isg = auc_mw(m66["isg_zfd"], m66["fatal_pub11"])
    a_ctl_isg, p_ctl_isg = auc_mw(-m66["ctl_zfd"] + m66["isg_zfd"], m66["fatal_pub11"])
    a_tri, p_tri = auc_mw(-m66["ctl_zfd"] - m66["nk_zfd"] + m66["isg_zfd"], m66["fatal_pub11"])
    print(f"  1st-draw Mann-Whitney AUC: -CTL={a_ctl:.3f}(P={p_ctl:.4g}), +ISG={a_isg:.3f}(P={p_isg:.4g}), "
          f"-CTL+ISG={a_ctl_isg:.3f}(P={p_ctl_isg:.4g}), -(CTL+NK)+ISG={a_tri:.3f}(P={p_tri:.4g})")

    # Leak-free in-fold rank-INT 50x5 CV comparing Base vs +CTL vs +ISG vs +CTL+ISG vs +CTL+NK+ISG
    cv_df = m66.copy()
    cv_df["sex_num"] = (cv_df["sex"].astype(str).str.upper().str.startswith("M")).astype(float)
    meta_age = p217_meta.drop_duplicates("patient_id").set_index("patient_id")[["age", "sex"]]
    cv_df["age"] = cv_df["patient_id"].map(meta_age["age"])
    cv_df["sex_num"] = (cv_df["patient_id"].map(meta_age["sex"]).astype(str).str.upper().str.startswith("M")).astype(float)

    def run_infold_cv(sub_df, base_cols, add_raw_cols, n_rep=50):
        y = sub_df["fatal_pub11"].values.astype(int)
        arms = sub_df["arm"].values
        oof_b_list, oof_e_list = [], []
        for rep in range(n_rep):
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=2026 + rep)
            pb = np.zeros(len(sub_df))
            pe = np.zeros(len(sub_df))
            for tr, te in skf.split(sub_df, y):
                Xb_tr = sub_df.iloc[tr][base_cols].values.astype(float)
                Xb_te = sub_df.iloc[te][base_cols].values.astype(float)
                sc_b = StandardScaler().fit(Xb_tr)
                mb = LogisticRegression(C=1.0, max_iter=500).fit(sc_b.transform(Xb_tr), y[tr])
                pb[te] = mb.predict_proba(sc_b.transform(Xb_te))[:, 1]

                # In-fold within-source rank-INT for each added feature
                add_tr_cols, add_te_cols = [], []
                for col in add_raw_cols:
                    v_tr = sub_df.iloc[tr][col].values.astype(float)
                    v_te = sub_df.iloc[te][col].values.astype(float)
                    z_tr = np.zeros_like(v_tr)
                    z_te = np.zeros_like(v_te)
                    for a_val in np.unique(arms):
                        mtr = arms[tr] == a_val
                        mte = arms[te] == a_val
                        sv = np.sort(v_tr[mtr])
                        r_tr = stats.rankdata(v_tr[mtr])
                        z_tr[mtr] = stats.norm.ppf((r_tr - 0.5) / len(sv))
                        # Monotone interpolation onto test fold
                        r_te = np.searchsorted(sv, v_te[mte], side="left") + 0.5
                        r_te = np.clip(r_te, 0.5, len(sv) - 0.5)
                        z_te[mte] = stats.norm.ppf(r_te / len(sv))
                    add_tr_cols.append(z_tr)
                    add_te_cols.append(z_te)
                Xe_tr = np.column_stack([Xb_tr] + add_tr_cols)
                Xe_te = np.column_stack([Xb_te] + add_te_cols)
                sc_e = StandardScaler().fit(Xe_tr)
                me = LogisticRegression(C=1.0, max_iter=500).fit(sc_e.transform(Xe_tr), y[tr])
                pe[te] = me.predict_proba(sc_e.transform(Xe_te))[:, 1]
            oof_b_list.append(auc_mw(pb, y)[0])
            oof_e_list.append(auc_mw(pe, y)[0])
        ab = float(np.mean(oof_b_list))
        ae = float(np.mean(oof_e_list))
        return ab, ae, ae - ab

    cv_rows = []
    for label, sub_cv, add_cols in [
        ("All 1st-dr (66) + T CD8 CTL", cv_df, ["T CD8 CTL"]),
        ("All 1st-dr (66) + T CD8 CTL ISG", cv_df, ["isg"]),
        ("All 1st-dr (66) + CTL + ISG", cv_df, ["T CD8 CTL", "isg"]),
        ("All 1st-dr (66) + CTL + NK", cv_df, ["T CD8 CTL", "NK"]),
        ("All 1st-dr (66) + CTL + NK + ISG", cv_df, ["T CD8 CTL", "NK", "isg"]),
        ("Early <=7d (53) + T CD8 CTL", cv_df[cv_df["days"] <= 7], ["T CD8 CTL"]),
        ("Early <=7d (53) + T CD8 CTL ISG", cv_df[cv_df["days"] <= 7], ["isg"]),
        ("Early <=7d (53) + CTL + ISG", cv_df[cv_df["days"] <= 7], ["T CD8 CTL", "isg"]),
    ]:
        ab, ae, da = run_infold_cv(sub_cv, ["age", "sex_num", "temp_severity_ord"], add_cols)
        cv_rows.append({"model": label, "n": len(sub_cv), "deaths": int(sub_cv["fatal_pub11"].sum()),
                        "AUC_base": ab, "AUC_ext": ae, "Delta_AUC": da})
        print(f"  [50x5 CV] {label:34s}: Base={ab:.3f} -> Ext={ae:.3f} (Delta={da:+.3f})")
    pd.DataFrame(cv_rows).to_csv(AUDIT_OUT / "A3_covid_tisg_cv_benchmark.csv", index=False)

    # 3B. Bacterial Sepsis T-cell ISG scoring in PBISC (122 samples)
    sep_h5 = ROOT / "deliverables" / "annotated_h5ad_20261112" / "pbisc_sepsis_ye2026_draw1_K4000_ceiling_annotated.h5ad"
    with h5py.File(sep_h5, "r") as f:
        def rcol(grp, k):
            it = grp[k]
            if isinstance(it, h5py.Group) and "codes" in it:
                cats = [c.decode("utf-8") if isinstance(c, bytes) else c for c in it["categories"][:]]
                return pd.Categorical.from_codes(it["codes"][:], categories=cats).astype(str)
            v = it[:]
            return np.array([x.decode("utf-8") if isinstance(x, bytes) else x for x in v])

        var_names = list(rcol(f["var"], "gene_symbol"))
        isg_present = [g for g in ISG19 if g in var_names]
        isg_idx = [var_names.index(g) for g in isg_present]
        sid = np.asarray(rcol(f["obs"], "sample_id"), dtype=str)
        cond = np.asarray(rcol(f["obs"], "condition"), dtype=str)
        maj = np.asarray(rcol(f["obs"], "celltype_major"), dtype=str)
        prim = np.asarray(rcol(f["obs"], "in_primary_scope"), dtype=str)
        t_mask = np.isin(maj, ["CD4+ T cells", "CD8+ T cells"])
        X_grp = f["X"]
        X_csr = sparse.csr_matrix((X_grp["data"][:], X_grp["indices"][:], X_grp["indptr"][:]), shape=(len(sid), len(var_names)))
        X_t = X_csr[t_mask]
        lib_t = np.asarray(X_t.sum(axis=1)).ravel()
        X_t_isg = X_t[:, isg_idx].toarray()
        obs_t = pd.DataFrame(X_t_isg, columns=isg_present)
        obs_t["sample_id"] = sid[t_mask]
        obs_t["condition"] = cond[t_mask]
        obs_t["in_primary_scope"] = prim[t_mask]
        obs_t["lib"] = lib_t
        pb_t = obs_t.groupby(["sample_id", "condition", "in_primary_scope"])[["lib"] + isg_present].sum().reset_index()
        pb_t["n_t_cells"] = obs_t.groupby("sample_id").size().reindex(pb_t["sample_id"]).values
        cp10k = np.log1p(1e4 * pb_t[isg_present].values / np.maximum(pb_t["lib"].values[:, None], 1.0))
        pb_t["T_ISG_score"] = cp10k.mean(axis=1)
        pb_ad = pb_t[pb_t["in_primary_scope"].isin(["True", "1"])].copy()
        pb_ad["ord"] = pb_ad["condition"].map({"adult_HC": 0, "adult_ICU_ctrl": 1, "adult_sepsis": 2})
        rho_s, p_s = stats.spearmanr(pb_ad["T_ISG_score"], pb_ad["ord"])
        _, p_sc = stats.mannwhitneyu(pb_ad[pb_ad["ord"] == 2]["T_ISG_score"], pb_ad[pb_ad["ord"] == 0]["T_ISG_score"])
        pb_t.to_csv(AUDIT_OUT / "A3b_sepsis_pbisc_tcell_isg.csv", index=False)
        print(f"\n[Sepsis PBISC T-cell ISG] {len(isg_present)}/19 ISG genes in panel | Adult n={len(pb_ad)}")
        for cnd in ["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]:
            sub = pb_ad[pb_ad["condition"] == cnd]["T_ISG_score"]
            print(f"  {cnd:16s} (n={len(sub):2d}): med={sub.median():.3f}, IQR=({sub.quantile(0.25):.3f}-{sub.quantile(0.75):.3f})")
        print(f"  Ordinal HC->ICU->Sepsis: Spearman rho={rho_s:+.3f}, P={p_s:.4g} | Sepsis vs HC MW P={p_sc:.4g}")

    # Also save cross-disease Sepsis vs COVID Real/PBISC/scSemi 19-ISG T-cell DEG log2FC table
    deg_dir = ROOT / "deliverables" / "annotated_h5ad_20261112" / "support" / "deg_tables"
    isg_deg_rows = []
    for coh in ["covid", "sepsis"]:
        for src in ["Real", "PBISC", "scSemiProfiler"]:
            df_d = pd.read_csv(deg_dir / f"deg_{coh}_Tcell_{src}.csv.gz")
            sub_d = df_d[df_d["gene"].isin(ISG19)].copy()
            sub_d["cohort"] = coh
            sub_d["source"] = src
            isg_deg_rows.append(sub_d[["cohort", "source", "gene", "logFC", "p", "q"]])
    pd.concat(isg_deg_rows, ignore_index=True).to_csv(AUDIT_OUT / "A3c_cross_disease_Tcell_19ISG_deg.csv", index=False)

    # Save first-draw T CD8 CTL abundance vs ISG table for plotting
    m66[["patient_id", "arm", "days", "temp_severity_ord", "fatal_pub11", "T CD8 CTL", "NK", "isg", "ctl_zfd", "nk_zfd", "isg_zfd"]].to_csv(
        AUDIT_OUT / "A3d_covid_66firstdraw_CTL_NK_ISG.csv", index=False
    )

    tisg_meta = {
        "ctl_vs_isg_213draws_rho": float(r_abund_isg_217), "ctl_vs_isg_213draws_p": float(p_abund_isg_217),
        "ctl_vs_isg_66firstdraw_rho": float(r_abund_isg_fd), "ctl_vs_isg_66firstdraw_p": float(p_abund_isg_fd),
        "firstdraw_AUC_negCTL": float(a_ctl), "firstdraw_p_negCTL": float(p_ctl),
        "firstdraw_AUC_posISG": float(a_isg), "firstdraw_p_posISG": float(p_isg),
        "firstdraw_AUC_negCTL_posISG": float(a_ctl_isg), "firstdraw_p_negCTL_posISG": float(p_ctl_isg),
        "firstdraw_AUC_negCTLNK_posISG": float(a_tri), "firstdraw_p_negCTLNK_posISG": float(p_tri),
        "sepsis_pbisc_tisg_rho": float(rho_s), "sepsis_pbisc_tisg_p": float(p_s), "sepsis_pbisc_tisg_sepsis_vs_HC_p": float(p_sc),
    }
    return cv_rows, tisg_meta


def main():
    s1, p217 = audit_1_metadata_matching()
    res_cnt, sep_cnt = audit_2_count_vs_percentage(p217)
    cv_rows, tisg_meta = audit_3_tcell_isg(p217)
    (AUDIT_OUT / "AUDIT_FOCUS3_SUMMARY.json").write_text(json.dumps({
        "metadata_audit": s1,
        "tisg_cv_benchmark": cv_rows,
        "tisg_summary": tisg_meta,
    }, indent=2), encoding="utf-8")

    # Export clean unified Excel workbook for the 3-priority audit
    xlsx_path = AUDIT_OUT / "PBISC_Metadata_CountVsPct_TISG_Audit_20261121.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
        pd.read_csv(D28 / "table_s3_all_107_donors_overview.csv").to_excel(xw, sheet_name="1a_COVID_107Donors_Match", index=False)
        pd.read_csv(D28 / "table_s2_sample_1to1_mapping_362.csv").to_excel(xw, sheet_name="1b_COVID_362Samples_1to1", index=False)
        pd.read_csv(AUDIT_OUT / "A1_sepsis_255_samples_matched_metadata.csv").to_excel(xw, sheet_name="1c_Sepsis_255Samples_Match", index=False)
        res_cnt.to_excel(xw, sheet_name="2a_COVID43_Count_vs_Pct_CBC", index=False)
        sep_cnt.to_excel(xw, sheet_name="2b_Sepsis44_Count_vs_Pct", index=False)
        pd.read_csv(OUT_ROOT / "draw_sameday_sev_firstdraw_death_20261120" / "tables" / "T6_Tcell_ISG_13subsets_new_framework.csv").to_excel(xw, sheet_name="3a_COVID_TISG_13Subsets", index=False)
        pd.read_csv(D30 / "T8b_Tcell_ISG_13subsets_time_windows.csv").to_excel(xw, sheet_name="3b_COVID_TISG_EarlyVsLate", index=False)
        pd.read_csv(AUDIT_OUT / "A3d_covid_66firstdraw_CTL_NK_ISG.csv").to_excel(xw, sheet_name="3c_COVID_66FirstDraw_ISG", index=False)
        pd.DataFrame(cv_rows).to_excel(xw, sheet_name="3d_COVID_TISG_50x5_CV", index=False)
        pd.read_csv(AUDIT_OUT / "A3b_sepsis_pbisc_tcell_isg.csv").to_excel(xw, sheet_name="3e_Sepsis_PBISC_TISG", index=False)
        pd.read_csv(AUDIT_OUT / "A3c_cross_disease_Tcell_19ISG_deg.csv").to_excel(xw, sheet_name="3f_CrossDisease_19ISG_DEG", index=False)
    print(f"\n[DONE] Saved all audit tables + 11-sheet Excel workbook to {xlsx_path}")


if __name__ == "__main__":
    main()
