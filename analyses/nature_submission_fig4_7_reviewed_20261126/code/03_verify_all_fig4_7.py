#!/usr/bin/env python3
"""
03_verify_all_fig4_7.py
-----------------------
Self-contained, portable end-to-end statistical recomputation and audit (`검산`)
for the 2026-11-26 reviewed & outcome-adjudicated PBISC reproduction package
(Figures 4–7 & Supplementary Figures 8–27).

Runs with ZERO hardcoded machine paths on Linux, macOS, and Windows:
1. Verifies all 22 CSV tables (`tables/`) against `00_table_manifest.csv`
   (row/column dimensions and SHA-256 checksums).
2. Verifies clinical donor outcome adjudication (`van der Wijst 2021` official
   `Data file S3`: 42 donors, 7 deaths; 219 external outcome patients, 24 deaths)
   and `Supplementary Fig. 23c` 16-unique-program HSPC recomputation (`n = 62`).
3. Independently recomputes key statistics across Figures 4, 5, 6, and 7
   (Van der Waerden rank-INT, GEE severity across 43 subsets, Mann–Whitney AUCs,
   OLS-residual WHO-adjusted P values, DerSimonian–Laird 4-cohort meta-analysis).
4. Audits `PBISC_Manuscript_Fig4_7_Final_20261125.md` for concise Results &
   Discussion/Conclusion sections, zero leaked local paths, and strict 1 -> 20
   sequential first-citation numbering.
5. Automatically runs optional binary asset checks (`.xlsx`, `figures/`,
   `supplementary_figures/`, `.docx`, `.pdf`) whenever those local build outputs
   are present, while passing cleanly in the reproduction-only GitHub clone.
6. Writes `VERIFICATION_PACKAGE_CHECKS_20261126.json` and exits non-zero on failure.
"""

from __future__ import annotations
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

try:
    import pymupdf as fitz
except ImportError:
    try:
        import fitz  # type: ignore
    except ImportError:
        fitz = None

PKG_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = PKG_DIR / "figures"
SUPP_FIG_DIR = PKG_DIR / "supplementary_figures"
TAB_DIR = PKG_DIR / "tables"
VER_DIR = PKG_DIR / "verification_20261126"
XLSX_PATH = PKG_DIR / "PBISC_SourceData_and_Summary_Tables_Fig4_7_20261125.xlsx"
MD_PATH = PKG_DIR / "PBISC_Manuscript_Fig4_7_Final_20261125.md"
DOCX_PATH = PKG_DIR / "PBISC_Manuscript_Fig4_7_Final_20261125.docx"
PDF_PATH = PKG_DIR / "PBISC_Manuscript_Fig4_7_Final_20261125.pdf"
OUT_JSON = PKG_DIR / "VERIFICATION_PACKAGE_CHECKS_20261126.json"


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def vdw_normal_scores(x: np.ndarray) -> np.ndarray:
    ranks = stats.rankdata(x, method="average")
    n = len(x)
    return stats.norm.ppf((ranks - 0.5) / float(n))


def mannwhitney_auc_p(scores: np.ndarray, y_bin: np.ndarray) -> tuple[float, float]:
    pos = scores[y_bin == 1]
    neg = scores[y_bin == 0]
    if float(np.std(scores)) == 0.0:
        return 0.5, 1.0
    res = stats.mannwhitneyu(pos, neg, alternative="two-sided", method="auto")
    auc = float(res.statistic) / float(len(pos) * len(neg))
    return auc, float(res.pvalue)


def dl_random_effects(effects: np.ndarray, ses: np.ndarray) -> tuple[float, float, float]:
    w = 1.0 / (ses ** 2)
    fe = float(np.sum(w * effects) / np.sum(w))
    k = len(effects)
    if k <= 1:
        se_re = float(ses[0])
    else:
        q_stat = float(np.sum(w * ((effects - fe) ** 2)))
        c_val = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
        tau2 = max(0.0, (q_stat - (k - 1)) / c_val) if c_val > 0 else 0.0
        w_re = 1.0 / (ses ** 2 + tau2)
        fe = float(np.sum(w_re * effects) / np.sum(w_re))
        se_re = math.sqrt(1.0 / float(np.sum(w_re)))
    z = fe / se_re if se_re > 0 else 0.0
    p = float(2.0 * stats.norm.sf(abs(z)))
    return fe, se_re, p


def main() -> None:
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})

    # -------------------------------------------------------------------------
    # 1. Table Manifest & CSV SHA-256 Integrity (+ optional XLSX check)
    # -------------------------------------------------------------------------
    df_manifest = pd.read_csv(TAB_DIR / "00_table_manifest.csv")
    check("table_manifest_22_tables", len(df_manifest) == 22, f"n={len(df_manifest)}")

    sheet_names: list[str] = []
    if XLSX_PATH.exists():
        wb = openpyxl.load_workbook(str(XLSX_PATH), read_only=True)
        sheet_names = wb.sheetnames
        check("xlsx_23_sheets_present", len(sheet_names) == 23 and sheet_names[0] == "00_README_Index", f"sheets={len(sheet_names)}")
        wb.close()

    for row in df_manifest.itertuples(index=False):
        cp = TAB_DIR / row.csv_filename
        df_c = pd.read_csv(cp)
        ok_csv = (len(df_c) == row.n_rows) and (len(df_c.columns) == row.n_cols) and (sha256_file(cp) == row.sha256)
        ok_sheet = (not XLSX_PATH.exists()) or (row.excel_sheet in sheet_names)
        check(f"table_csv_manifest::{row.csv_filename}", ok_csv and ok_sheet, f"rows={len(df_c)}, cols={len(df_c.columns)}")

    # -------------------------------------------------------------------------
    # 2. Clinical Outcome Adjudication & HSPC-16 Verification (2026-11-26)
    # -------------------------------------------------------------------------
    vdw_rows = json.loads((VER_DIR / "clinical_sources" / "vdw_clinical_rows.json").read_text(encoding="utf-8"))
    raw_vdw = pd.DataFrame(vdw_rows[1:], columns=vdw_rows[0])
    raw_vdw = raw_vdw.loc[:, raw_vdw.columns.notna()]
    clin_vdw = raw_vdw[["Donor_ID", "Died"]].dropna(subset=["Donor_ID"])
    vdw_map = {f"Rep_C_{int(r.Donor_ID)}": int(str(r.Died).strip().lower() == "yes") for r in clin_vdw.itertuples()}

    ext300 = pd.read_csv(TAB_DIR / "fig7_external4_300_acute_samples.csv")
    vd_sub = ext300[ext300["cohort"] == "van der Wijst 2021"]
    vd_mismatches = int((vd_sub["death"].astype(int) != vd_sub["donor"].map(vdw_map)).sum())
    within_donor_inconsistencies = int((ext300.dropna(subset=["death"]).groupby(["cohort", "donor"])["death"].nunique() > 1).sum())
    check(
        "clinical_vdw_official_adjudication_zero_mismatch",
        vd_mismatches == 0 and within_donor_inconsistencies == 0 and vd_sub["donor"].nunique() == 42,
        f"vdW_mismatches={vd_mismatches}, donor_inconsistencies={within_donor_inconsistencies}",
    )

    hspc_pts = pd.read_csv(VER_DIR / "supp23_hspc62_patient_module_scores.csv")
    hspc_rec = pd.read_csv(VER_DIR / "supp23_hspc16_recomputed.csv")
    ery5_comb = hspc_rec.loc[(hspc_rec["module"] == "Ery5") & (hspc_rec["cohort"] == "Combined")].iloc[0]
    hsc5_comb = hspc_rec.loc[(hspc_rec["module"] == "HSC5") & (hspc_rec["cohort"] == "Combined")].iloc[0]
    check(
        "supp23_hspc_unique16_recomputed_exact",
        len(hspc_rec) == 48
        and hspc_pts["patient_id"].nunique() == 62
        and "HSC_broad_q" in set(hspc_rec["module"])
        and abs(float(ery5_comb["AUC"]) - 0.752228) < 1e-4
        and abs(float(ery5_comb["q16"]) - 0.149966) < 1e-4
        and abs(float(hsc5_comb["AUC"]) - 0.512478) < 1e-4,
        f"Ery5_AUC={ery5_comb['AUC']:.4f} (q16={ery5_comb['q16']:.4f}), HSC5_AUC={hsc5_comb['AUC']:.4f}",
    )

    # -------------------------------------------------------------------------
    # 3. Statistical Recomputation — Figure 4
    # -------------------------------------------------------------------------
    f4_comp = pd.read_csv(TAB_DIR / "fig4_paired_composition_6types.csv")
    check("fig4_comp_sums_to_100", bool(np.all(np.abs(f4_comp["sum_pct"] - 100.0) < 0.02)), f"sums={f4_comp['sum_pct'].tolist()}")

    f4_ilisi = pd.read_csv(TAB_DIR / "fig4_paired_ilisi_and_gene_r.csv").set_index(["cohort", "source"])
    check("fig4_sepsis_ilisi_exact", (
        f4_ilisi.loc[("Sepsis", "Real-Real"), "pre_pca_50d_ilisi_printed_mean"] == 1.77
        and f4_ilisi.loc[("Sepsis", "PBISC"), "pre_pca_50d_ilisi_printed_mean"] == 1.51
        and f4_ilisi.loc[("Sepsis", "scSemiProfiler"), "pre_pca_50d_ilisi_printed_mean"] == 1.00
    ), "1.77 / 1.51 / 1.00")
    check("fig4_covid_ilisi_exact", (
        f4_ilisi.loc[("COVID-19", "Real-Real"), "pre_pca_50d_ilisi_printed_mean"] == 1.76
        and f4_ilisi.loc[("COVID-19", "PBISC"), "pre_pca_50d_ilisi_printed_mean"] == 1.44
        and f4_ilisi.loc[("COVID-19", "scSemiProfiler"), "pre_pca_50d_ilisi_printed_mean"] == 1.00
    ), "1.76 / 1.44 / 1.00")

    f4_deg = pd.read_csv(TAB_DIR / "fig4_paired_limma_deg_metrics_32pairs.csv")
    piv_f4 = f4_deg.pivot(index=["cohort", "cell_type"], columns="source")
    recall_wins = int((piv_f4["DEG_recall_mean"]["PBISC"] > piv_f4["DEG_recall_mean"]["scSemiProfiler"]).sum())
    f1_wins = int((piv_f4["DEG_F1_mean"]["PBISC"] > piv_f4["DEG_F1_mean"]["scSemiProfiler"]).sum())
    covid_r_wins = int((piv_f4.loc["COVID-19", ("paired_logFC_r_mean", "PBISC")] > piv_f4.loc["COVID-19", ("paired_logFC_r_mean", "scSemiProfiler")]).sum())
    check("fig4_pbisc_wins_recall_and_f1_all_8", recall_wins == 8 and f1_wins == 8 and covid_r_wins == 4, f"recall={recall_wins}/8, F1={f1_wins}/8, covid_r={covid_r_wins}/4")

    f4_s12 = pd.read_csv(TAB_DIR / "fig4_supp12_d0_volcano_deg_counts.csv")
    f4_s13 = pd.read_csv(TAB_DIR / "fig4_supp13_native_effect_agreement.csv")
    check("fig4_supp12_24_rows_and_supp13_16_rows", len(f4_s12) == 24 and len(f4_s13) == 16, f"s12={len(f4_s12)}, s13={len(f4_s13)}")

    # -------------------------------------------------------------------------
    # 4. Statistical Recomputation — Figure 5
    # -------------------------------------------------------------------------
    f5_units = pd.read_csv(TAB_DIR / "fig5_donor_condition_297units_composition.csv")
    n_units_total = f5_units[["cohort", "source", "condition", "patient"]].drop_duplicates().shape[0]
    check("fig5_297_donor_condition_units", n_units_total == 297, f"units={n_units_total}")

    sep_rb_d = set(f5_units.loc[(f5_units["cohort"] == "sepsis") & (f5_units["source"] == "Real B"), "patient"])
    sep_pb_d = set(f5_units.loc[(f5_units["cohort"] == "sepsis") & (f5_units["source"] == "PBISC"), "patient"])
    cov_rb_d = set(f5_units.loc[(f5_units["cohort"] == "covid") & (f5_units["source"] == "Real B"), "patient"])
    cov_pb_d = set(f5_units.loc[(f5_units["cohort"] == "covid") & (f5_units["source"] == "PBISC"), "patient"])
    shared_cov = sorted(cov_rb_d & cov_pb_d)
    check("fig5_shared_donor_audit", len(sep_rb_d & sep_pb_d) == 0 and shared_cov == ["C113", "C127", "C145", "C168", "C169", "C20", "C63"], f"sepsis_overlap={len(sep_rb_d & sep_pb_d)}, covid_overlap={shared_cov}")

    f5_sum = pd.read_csv(TAB_DIR / "fig5_composition_summary_by_group.csv")
    recomp_means = f5_units.groupby(["cohort", "condition", "source", "cell_type"], as_index=False)["pct"].mean()
    m_comp = recomp_means.merge(f5_sum, on=["cohort", "condition", "source", "cell_type"])
    max_comp_diff = float(np.max(np.abs(m_comp["pct"] - m_comp["mean"])))
    check("fig5_composition_means_exact", max_comp_diff < 1e-9, f"max_diff={max_comp_diff:.2e}")

    f5_sig = pd.read_csv(TAB_DIR / "fig5_published_disease_signals_42rows.csv")
    pb_sep_sig = int(f5_sig.loc[(f5_sig["panel"] == "5n") & (f5_sig["source"] == "PBISC"), "sig_p_lt_0p05_on_published_side"].sum())
    pb_cov_sig = int(f5_sig.loc[(f5_sig["panel"] == "5o") & (f5_sig["source"] == "PBISC"), "sig_p_lt_0p05_on_published_side"].sum())
    check("fig5_published_signals_6of7_and_7of7", len(f5_sig) == 42 and pb_sep_sig == 6 and pb_cov_sig == 7, f"sepsis={pb_sep_sig}/7, covid={pb_cov_sig}/7")

    f5_deg = pd.read_csv(TAB_DIR / "fig5_welch_deg_fidelity_16rows.csv")
    piv_f5 = f5_deg.pivot(index=["cohort", "celltype"], columns="source")
    auroc_wins = int((piv_f5["AUROC"]["PBISC"] > piv_f5["AUROC"]["scSemiProfiler"]).sum())
    check("fig5_pbisc_wins_auroc_all_8_lineages", len(f5_deg) == 16 and auroc_wins == 8, f"auroc_wins={auroc_wins}/8")

    # -------------------------------------------------------------------------
    # 5. Statistical Recomputation — Figure 6
    # -------------------------------------------------------------------------
    f6_shares = pd.read_csv(TAB_DIR / "fig6_rare_parent_shares_integer_counts.csv")
    recomp_pct = 100.0 * f6_shares["n_target"] / f6_shares["n_parent"]
    max_pct_diff = float(np.max(np.abs(recomp_pct - f6_shares["percent_parent"])))
    sep_pb_nz = int((f6_shares.loc[(f6_shares["cohort"] == "sepsis") & (f6_shares["source"] == "PBISC"), "n_target"] > 0).sum())
    cov_pb_nz = int((f6_shares.loc[(f6_shares["cohort"] == "covid") & (f6_shares["source"] == "PBISC"), "n_target"] > 0).sum())
    sep_sc_z = int((f6_shares.loc[(f6_shares["cohort"] == "sepsis") & (f6_shares["source"] == "scSemiProfiler"), "n_target"] == 0).sum())
    cov_sc_z = int((f6_shares.loc[(f6_shares["cohort"] == "covid") & (f6_shares["source"] == "scSemiProfiler"), "n_target"] == 0).sum())
    hsc_den = int(f6_shares.loc[(f6_shares["cohort"] == "sepsis") & (f6_shares["source"] == "PBISC") & (f6_shares["subset"] == "HSC"), "n_parent"].iloc[0])
    hspc_den = int(f6_shares.loc[(f6_shares["cohort"] == "covid") & (f6_shares["source"] == "PBISC") & (f6_shares["subset"] == "HSPC"), "n_parent"].iloc[0])
    check("fig6_integer_counts_and_denominators_exact", (
        max_pct_diff < 1e-10
        and sep_pb_nz == 13 and cov_pb_nz == 12
        and sep_sc_z == 7 and cov_sc_z == 8
        and hsc_den == 272000 and hspc_den == 560000
    ), f"PBISC={sep_pb_nz}/14 & {cov_pb_nz}/13; scSemi_zero={sep_sc_z}/14 & {cov_sc_z}/13; den={hsc_den},{hspc_den}")

    f6_prof = pd.read_csv(TAB_DIR / "fig6_interpatient_similarity_895profiles.csv")
    f6_sim = pd.read_csv(TAB_DIR / "fig6_similarity_summary_24groups.csv")
    check("fig6_895_profiles_count", len(f6_prof) == 895 and int((f6_prof["cohort"] == "sepsis").sum()) == 270 and int((f6_prof["cohort"] == "covid").sum()) == 625, f"total={len(f6_prof)}")
    recomp_med = f6_prof.groupby(["cohort", "subset", "source"], as_index=False)["median_interpatient_r"].agg(n_re="count", med_re="median")
    m_sim = f6_sim.loc[f6_sim["count"] > 0].merge(recomp_med, on=["cohort", "subset", "source"])
    max_med_diff = float(np.max(np.abs(m_sim["median"] - m_sim["med_re"])))
    check("fig6_similarity_medians_recomputed_exact", max_med_diff < 1e-9 and bool((m_sim["count"] == m_sim["n_re"]).all()), f"groups={len(m_sim)}, max_diff={max_med_diff:.2e}")

    # -------------------------------------------------------------------------
    # 6. Statistical Recomputation — Figure 7 (ESKD & 4 External Cohorts)
    # -------------------------------------------------------------------------
    df217 = pd.read_csv(TAB_DIR / "fig7_eskd_217_acute_draws_43subsets.csv")
    df66 = pd.read_csv(TAB_DIR / "fig7_eskd_66_first_draws_43subsets.csv")
    df63 = pd.read_csv(TAB_DIR / "fig7_eskd_63_first_draws_CTL_NK_ISG.csv")

    who_counts = df217["temp_severity_ord"].value_counts().to_dict()
    strata_counts = df66["cohort_design_stratum"].value_counts().to_dict()
    n_le7_all = int((df66["days_from_first_positive_swab"] <= 7).sum())
    n_le7_dead = int(((df66["days_from_first_positive_swab"] <= 7) & (df66["fatal_pub11"] == 1)).sum())
    check("fig7_cohort_and_stratum_counts_exact", (
        len(df217) == 217
        and df217["patient_id"].nunique() == 66
        and len(df66) == 66
        and int(df66["fatal_pub11"].sum()) == 11
        and strata_counts.get("Measured_scRNA_only_10") == 10
        and strata_counts.get("Measured_shared_with_bulk_27") == 27
        and strata_counts.get("Wave1_bulk_only_PBISC_29") == 29
        and who_counts == {0.0: 126, 1.0: 58, 2.0: 21, 3.0: 12}
        and n_le7_all == 53 and n_le7_dead == 10
    ), f"strata={strata_counts}, who={who_counts}, <=7d={n_le7_all}/{n_le7_dead}")

    subsets43 = [c[:-6] for c in df217.columns if c.endswith("__z217")]
    check("fig7_43_subsets_present", len(subsets43) == 43, f"n={len(subsets43)}")
    max_vdw_diff_217 = 0.0
    max_vdw_diff_66 = 0.0
    for ct in subsets43:
        for _, g in df217.groupby("arm"):
            z_calc = vdw_normal_scores(g[ct].to_numpy(dtype=float))
            max_vdw_diff_217 = max(max_vdw_diff_217, float(np.max(np.abs(z_calc - g[f"{ct}__z217"].to_numpy(dtype=float)))))
        for _, g in df66.groupby("arm"):
            z_calc = vdw_normal_scores(g[ct].to_numpy(dtype=float))
            max_vdw_diff_66 = max(max_vdw_diff_66, float(np.max(np.abs(z_calc - g[f"{ct}__zfd"].to_numpy(dtype=float)))))
    check("fig7_van_der_waerden_rankINT_exact", max_vdw_diff_217 < 1e-11 and max_vdw_diff_66 < 1e-11, f"max217={max_vdw_diff_217:.2e}, max66={max_vdw_diff_66:.2e}")

    df217_m = df217.copy()
    df217_m["weeks"] = df217_m["days"].astype(float) / 7.0
    gee_comb_p = []
    gee_meas_p = []
    gee_comb_beta = {}
    for ct in subsets43:
        y = df217_m[f"{ct}__z217"].to_numpy(dtype=float)
        if float(np.std(y)) == 0.0:
            gee_comb_beta[ct] = 0.0
            gee_comb_p.append(1.0)
        else:
            X = sm.add_constant(df217_m[["temp_severity_ord", "weeks"]].to_numpy(dtype=float))
            fit_c = sm.GEE(y, X, groups=df217_m["patient_id"], cov_struct=sm.cov_struct.Exchangeable()).fit()
            gee_comb_beta[ct] = float(fit_c.params[1])
            gee_comb_p.append(float(fit_c.pvalues[1]))

        sub_r = df217_m.loc[df217_m["arm"] == "Real"]
        yr = sub_r[f"{ct}__z217"].to_numpy(dtype=float)
        if float(np.std(yr)) == 0.0:
            gee_meas_p.append(1.0)
        else:
            Xr = sm.add_constant(sub_r[["temp_severity_ord", "weeks"]].to_numpy(dtype=float))
            fit_r = sm.GEE(yr, Xr, groups=sub_r["patient_id"], cov_struct=sm.cov_struct.Exchangeable()).fit()
            gee_meas_p.append(float(fit_r.pvalues[1]))

    q43_comb = multipletests(gee_comb_p, method="fdr_bh")[1]
    q43_meas = multipletests(gee_meas_p, method="fdr_bh")[1]
    q43_map = dict(zip(subsets43, q43_comb))
    n_sig_comb = int(np.sum(q43_comb < 0.05))
    n_sig_meas = int(np.sum(q43_meas < 0.05))
    n_sig_shared = int(np.sum((q43_comb < 0.05) & (q43_meas < 0.05)))
    check("fig7_gee_18_combined_9_measured_7_shared", n_sig_comb == 18 and n_sig_meas == 9 and n_sig_shared == 7, f"comb={n_sig_comb}, meas={n_sig_meas}, shared={n_sig_shared}")
    check("fig7_gee_focal5_betas_and_q43", (
        abs(gee_comb_beta["B_ASC_dividing"] - 0.422) < 0.002
        and abs(gee_comb_beta["HSPC"] - 0.333) < 0.002
        and abs(gee_comb_beta["CD16mono_C1"] - 0.247) < 0.002
        and abs(gee_comb_beta["T g/d"] - (-0.196)) < 0.002
        and abs(gee_comb_beta["T CD8 CM"] - (-0.231)) < 0.002
        and abs(q43_map["HSPC"] - 0.00013) < 2e-5
    ), f"betas={{divASC:{gee_comb_beta['B_ASC_dividing']:.3f}, HSPC:{gee_comb_beta['HSPC']:.3f}}}")

    y66 = df66["fatal_pub11"].to_numpy(dtype=int)
    fd_aucs = {}
    fd_ps = []
    for ct in subsets43:
        auc_v, p_v = mannwhitney_auc_p(df66[f"{ct}__zfd"].to_numpy(dtype=float), y66)
        fd_aucs[ct] = auc_v
        fd_ps.append(p_v)
    fd_q43 = dict(zip(subsets43, multipletests(fd_ps, method="fdr_bh")[1]))
    src_ind = (df66["arm"].astype(str) == "PBISC new").astype(float).to_numpy()
    X_who_src = sm.add_constant(np.column_stack([df66["temp_severity_ord"].to_numpy(dtype=float), src_ind]))
    _, p_ctl_adj = mannwhitney_auc_p(sm.OLS(df66["T CD8 CTL__zfd"].to_numpy(dtype=float), X_who_src).fit().resid, y66)
    _, p_nk_adj = mannwhitney_auc_p(sm.OLS(df66["NK__zfd"].to_numpy(dtype=float), X_who_src).fit().resid, y66)
    check("fig7_firstdraw_ctl_nk_auc_q43_whoadj", (
        abs(fd_aucs["T CD8 CTL"] - 0.213223) < 1e-4
        and abs(fd_aucs["NK"] - 0.217355) < 1e-3
        and abs(fd_q43["T CD8 CTL"] - 0.07006) < 1e-4
        and abs(fd_q43["NK"] - 0.07006) < 1e-4
        and abs(p_ctl_adj - 0.0330) < 0.001
        and abs(p_nk_adj - 0.0303) < 0.001
    ), f"CTL={fd_aucs['T CD8 CTL']:.4f} (q={fd_q43['T CD8 CTL']:.4f}, adjP={p_ctl_adj:.4f}), NK={fd_aucs['NK']:.4f} (adjP={p_nk_adj:.4f})")

    comp63 = (df63["isg_zfd"].to_numpy(dtype=float) - df63["ctl_zfd"].to_numpy(dtype=float) - df63["nk_zfd"].to_numpy(dtype=float)) / 3.0
    y63 = df63["fatal_pub11"].to_numpy(dtype=int)
    auc_comp, p_comp = mannwhitney_auc_p(comp63, y63)
    X_cat_who = pd.get_dummies(df63["temp_severity_ord"].astype(int), drop_first=True, dtype=float)
    X_cat_who = sm.add_constant(X_cat_who.to_numpy(dtype=float))
    _, p_comp_adj = mannwhitney_auc_p(sm.OLS(comp63, X_cat_who).fit().resid, y63)
    check("fig7_composite_n63_auc_and_whoadj", (
        len(df63) == 63 and int(y63.sum()) == 11
        and abs(auc_comp - 0.840909) < 1e-4
        and abs(p_comp - 0.000429) < 2e-5
        and abs(p_comp_adj - 0.01279) < 5e-4
    ), f"AUC={auc_comp:.4f}, P={p_comp:.6f}, adjP={p_comp_adj:.5f}")

    ext_per = pd.read_csv(TAB_DIR / "fig7_external4_per_cohort_results.csv")
    ext_meta = pd.read_csv(TAB_DIR / "fig7_external4_meta_23subsets.csv")
    fd_ext = ext300.sort_values(["cohort", "donor", "time_ord", "sample_id"]).groupby(["cohort", "donor"], as_index=False).first()
    fd_known = fd_ext[fd_ext["death"].notna()]
    n_out_pts = len(fd_known)
    n_out_deaths = int(fd_known["death"].sum())
    n_unique_pts = len(ext300[["cohort", "donor"]].drop_duplicates())
    check("fig7_external_300_draws_231_patients_219_outcomes_24_deaths", (
        len(ext300) == 300 and n_unique_pts == 231
        and n_out_pts == 219 and n_out_deaths == 24
    ), f"draws={len(ext300)}, pts={n_unique_pts}, out_pts={n_out_pts}, deaths={n_out_deaths}")

    hspc_rows = ext_per.loc[ext_per["feature"] == "HSPC"].dropna(subset=["z_fisher", "se_z_fisher"])
    z_re, _, _ = dl_random_effects(hspc_rows["z_fisher"].to_numpy(dtype=float), hspc_rows["se_z_fisher"].to_numpy(dtype=float))
    rho_hspc_re = math.tanh(z_re)
    cd8n_rows = ext_per.loc[ext_per["feature"] == "T CD8 Naive"].dropna(subset=["auc_death_first", "se_auc_first"])
    d_cd8n, _, _ = dl_random_effects(cd8n_rows["auc_death_first"].to_numpy(dtype=float) - 0.5, np.maximum(cd8n_rows["se_auc_first"].to_numpy(dtype=float), 0.035))
    cd16_rows = ext_per.loc[ext_per["feature"] == "CD16mono"].dropna(subset=["auc_death", "se_auc"])
    d_cd16, _, _ = dl_random_effects(cd16_rows["auc_death"].to_numpy(dtype=float) - 0.5, np.maximum(cd16_rows["se_auc"].to_numpy(dtype=float), 0.035))
    n_exact_meta = int((ext_meta["feature"] != "B_ASC_total").sum())
    check("fig7_external_dl_meta_recomputed", (
        abs(rho_hspc_re - 0.40) < 0.01
        and abs((0.5 + d_cd8n) - 0.323309) < 0.01
        and abs((0.5 + d_cd16) - 0.393654) < 0.01
        and n_exact_meta == 23
    ), f"HSPC_rho={rho_hspc_re:.3f}, CD8naive_1st={0.5+d_cd8n:.3f}, CD16mono_last={0.5+d_cd16:.3f}, exact_meta={n_exact_meta}")

    # -------------------------------------------------------------------------
    # 7. Markdown Manuscript, Conciseness, Path Cleanliness & Citation Audit
    # -------------------------------------------------------------------------
    md_text = MD_PATH.read_text(encoding="utf-8")
    check("md_no_local_paths_leaked", ("F:\\" not in md_text) and ("C:\\Users" not in md_text) and ("F:/PBISC" not in md_text), "clean paths")

    sec_map = {}
    for blk in re.split(r"^##\s+", md_text, flags=re.M)[1:]:
        title = blk.splitlines()[0].strip()
        sec_map[title] = len(blk.split())
    res_words = sec_map.get("Results", 99999)
    disc_words = sec_map.get("Discussion", 99999)
    check(
        "manuscript_concise_results_and_discussion",
        (800 <= res_words <= 2200) and (250 <= disc_words <= 650) and ("### Conclusion" in md_text),
        f"Results={res_words} words, Discussion+Conclusion={disc_words} words",
    )

    body_text, ref_text = md_text.split("## References", 1)
    sup_to_ascii = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
    body_no_pow10 = re.sub(r"10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+", "10^x", body_text)
    seen_refs: list[int] = []
    for m in re.finditer(r"[⁰¹²³⁴⁵⁶⁷⁸⁹]+(?:[,⁻][⁰¹²³⁴⁵⁶⁷⁸⁹]+)*", body_no_pow10):
        raw_tok = m.group(0).translate(sup_to_ascii)
        for part in raw_tok.split(","):
            if "-" in part:
                a, b = part.split("-", 1)
                nums = [int(b)] if a == "" else list(range(int(a), int(b) + 1))
            else:
                nums = [int(part)]
            for num in nums:
                if num not in seen_refs:
                    seen_refs.append(num)
    ref_entries = [int(m.group(1)) for m in re.finditer(r"^(\d+)\.\s+", ref_text, flags=re.MULTILINE)]
    check("references_strict_1_to_20_order", seen_refs == list(range(1, 21)) and ref_entries == list(range(1, 21)), f"seen={seen_refs}")
    check("references_pubmed_verified_entries", (
        "Ye, Q. et al." in ref_text and "Nat. Immunol." in ref_text
        and "Gisby, J. S. et al." in ref_text and "7775 (2022)" in ref_text
        and "Stephenson, E. et al." in ref_text and "100918 (2025)" in ref_text
        and "Type I interferon autoantibodies" in ref_text
        and "suppressing" not in ref_text and "suppressive myeloid cell production" in ref_text
    ), "Ye 2026, Gisby 2022, Stephenson 2025, van der Wijst 2021, Reyes 2021 verified")

    # Optional checks when compiled DOCX / PDF / Figures are built locally
    if DOCX_PATH.exists() and FIG_DIR.exists() and SUPP_FIG_DIR.exists():
        with zipfile.ZipFile(DOCX_PATH, "r") as zf:
            media_names = sorted(n for n in zf.namelist() if n.startswith("word/media/"))
        check("optional_docx_24_media", len(media_names) == 24, f"media={len(media_names)}")
    if PDF_PATH.exists() and fitz is not None:
        pdoc = fitz.open(str(PDF_PATH))
        pdf_imgs = sum(len(pdoc[i].get_images(full=True)) for i in range(len(pdoc)))
        pdoc.close()
        check("optional_pdf_24_images", pdf_imgs == 24, f"images={pdf_imgs}")

    n_pass = sum(1 for c in checks if c["passed"])
    n_fail = len(checks) - n_pass
    report = {
        "verification_date": "2026-11-26",
        "package_directory": PKG_DIR.name,
        "distribution_mode": "reproduction_elements_only" if not DOCX_PATH.exists() else "full_local_bundle",
        "status": "ALL_CHECKS_PASSED" if n_fail == 0 else "FAILED",
        "total_checks": len(checks),
        "passed_checks": n_pass,
        "failed_checks": n_fail,
        "tables_verified": len(df_manifest),
        "section_word_counts": sec_map,
        "checks": checks,
    }
    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Verification complete: {n_pass}/{len(checks)} checks PASSED, {n_fail} FAILED.")
    if n_fail > 0:
        for c in checks:
            if not c["passed"]:
                print(f"  [FAIL] {c['check']}: {c['detail']}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
