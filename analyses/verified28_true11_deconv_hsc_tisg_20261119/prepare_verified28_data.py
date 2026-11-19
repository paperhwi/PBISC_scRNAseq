"""Stage 1 Data Freezer & Complete Verification Engine (Folder 28: Verified True-11 Mortality,
Expanded BayesPrism/CIBERSORTx/MuSiC Deconvolution, HSPC Transcriptome, and T-cell ISG).

Addresses all 6 user verification points from primary raw files:
  Point 1 (Metadata Forensic Lock):
    - WHO severity (0-3): 100% verified (0 discrepancies across 68 COVID+ patients).
    - Mortality: Corrected from corrupted 13-death hybrid to TRUE PUBLISHED 11 COVID-19 deaths
      (Real n=37: 5 deaths [C21, C33, C141 direct + C108, C145 indirect], verified against BOTH
       Gisby et al. Nat Commun 2022 Fig. 1/p.14 & Stephenson et al. Cell Genomics 2025 Fig. 1B/Table S1;
       PBISC new n=29: 6 direct Wave-1 deaths [C32, C37, C48, C66, C77, C84]).
      Fixes the Zenodo w2_metadata.csv bug that erroneously coded survivors C101, C146, C147 as fatal.
      Also retains Direct-9 (3R+6N) and symmetric Zenodo-raw-14 (8R+6N) as sensitivity variants.
  Point 2 (Clear Patient & Sample Tables):
    - Exports verified 68-COVID patient table, 107-donor overview, 362-sample 1:1 mapping,
      and formatted multi-sheet Excel workbook `Table_S1_Verified_Patient_and_Sample_Metadata_20261119.xlsx`.
  Point 3 (Logical Soundness & 100% Standalone Reproducibility):
    - Recomputes all Figure 7 & 8 statistics under TRUE_PUB_11 (43-subset scan, 5-step funnel,
      3,000-resample source-stratified CIs, exact Source x Peak WHO permutation [8,400 allocations],
      unified 1,000-iteration conditional re-detection, 14-variant sensitivity matrix).
    - Freezes SHA-256 of every upstream input in INPUT_MANIFEST.csv and all local CSVs in data/.
  Point 4 (Expanded Unbiased Bulk Deconvolution Benchmark):
    - Evaluates 10 methods on identical Zenodo 6497251 HTSeq counts & strictly out-of-donor Real reference:
      PBISC, BayesPrism (updated psi & initial phi + posterior Z), CIBERSORTx-style S-mode & B-mode,
      CIBERSORT-SVR, MuSiC, Calibrated NNLS, DWLS, Uncalibrated NNLS, Bulk 8-panel markers (pure BH8).
    - Scores on Shared-27 (all 43 & 23 rare <1%), unbiased T1 (Real-37) & T2 (4 external public cohorts),
      Wave-1 new (n=29) rare-cell clinical signals, Combined (n=66, True 11d) recovery, AND PBISC vs BayesPrism Z.
  Point 5 (HSC / HSPC Transcriptome vs Severity & True-11 Mortality):
    - Disjoint n=62 patients (Real 33 [5d] + PBISC new 29 [6d], 3,124 HSPC cells), exact 16 unique biological
      modules (deduplicating Cycling/Cycle4) + 20 candidate genes + full 5,762-gene scan (KIAA0391).
  Point 6 (T-cell Type-I IFN / ISG vs Severity & True-11 Mortality):
    - All 13 T-cell compartments (All T, CD4, CD8, CD4/CD8 CTL, CM, Naive, Treg, T g/d, MAIT):
      Shared-27 fidelity, COVID+ vs COVID-, patient-level severity (BH13) & True-11 death (BH13),
      early (0-7d) vs late (>7d) windows, sampling-day-adjusted GEE, and source-specific Time x Fatal GEE.
"""
import argparse
import hashlib
import json
import pathlib
import sys
import warnings
from itertools import combinations

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import nnls
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import NuSVR
from statsmodels.stats.multitest import multipletests

DEFAULT_OUT = pathlib.Path(r"F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\28_Verified_True11_Deconv_HSC_TISG_20261119")
ROOT = pathlib.Path(r"F:\PBISC")
ESKD_DIR = ROOT / "eskd_covid_clinical_assoc_20261115"
BULK_DIR = ROOT / "bulkonly_lora5_20261006" / "analysis" / "bulk_public" / "zenodo6497251_imperial_eskd"
SUM_DIR = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116"
D15_DIR = SUM_DIR / "15_Figure7_ver8_Main_and_Supp" / "data"
D21_DIR = SUM_DIR / "21_Fig7_Fig8_Corrected_Publication_20261118" / "data"
D24_DIR = SUM_DIR / "24_OptionA_Audit23_Resolved_20261118"

CLL = {"C116", "C139"}
# Verified Published COVID-19 Mortality Ground Truth (Gisby 2022 Nat Commun + Stephenson 2025 Cell Genomics)
D_PUB11 = {"C21", "C33", "C141", "C108", "C145", "C32", "C37", "C48", "C66", "C77", "C84"}
D_DIR9 = {"C21", "C33", "C141", "C32", "C37", "C48", "C66", "C77", "C84"}
D_ZEN14 = D_PUB11 | {"C101", "C146", "C147"}
D_LEG13 = D_PUB11 | {"C101", "C147"}


def spearman_vec_cols(s_vec: np.ndarray, X: np.ndarray):
    n = len(s_vec)
    rs = stats.rankdata(s_vec) - (n + 1) / 2.0
    rx = stats.rankdata(X, axis=0) - (n + 1) / 2.0
    ss_s = float(np.sum(rs * rs))
    ss_x = np.sum(rx * rx, axis=0)
    denom = np.sqrt(ss_s * ss_x)
    r = np.where(denom > 0, (rx.T @ rs) / np.maximum(denom, 1e-30), 0.0)
    r_c = np.clip(r, -1.0 + 1e-15, 1.0 - 1e-15)
    t_st = r_c * np.sqrt((n - 2.0) / (1.0 - r_c * r_c))
    p = np.where(denom > 0, 2.0 * stats.t.sf(np.abs(t_st), df=n - 2), 1.0)
    return r, p


def mw_vec_cols(f_vec: np.ndarray, X: np.ndarray):
    m1, m0 = f_vec == 1, f_vec == 0
    n1, n0 = int(m1.sum()), int(m0.sum())
    u_res = stats.mannwhitneyu(X[m1], X[m0], axis=0)
    std_x = np.std(X, axis=0)
    aucs = np.where(std_x > 0, u_res.statistic / (n1 * n0), 0.5)
    ps = np.where(std_x > 0, u_res.pvalue, 1.0)
    return aucs, ps


def rint_arr(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    if len(np.unique(v)) < 2:
        return np.zeros_like(v)
    rk = stats.rankdata(v, method="average")
    return stats.norm.ppf((rk - 0.5) / len(v))


def rint(x: pd.Series) -> pd.Series:
    x = pd.Series(x, dtype=float)
    if x.nunique(dropna=True) < 2:
        return pd.Series(0.0, index=x.index)
    r = x.rank()
    return pd.Series(stats.norm.ppf((r - 0.5) / x.notna().sum()), index=x.index)


def auc_mw(y, f):
    y_a, f_a = np.asarray(y, float), np.asarray(f, int)
    m = np.isfinite(y_a)
    y_a, f_a = y_a[m], f_a[m]
    n1, n0 = int((f_a == 1).sum()), int((f_a == 0).sum())
    if n1 < 1 or n0 < 1 or np.nanstd(y_a) == 0:
        return 0.5, 1.0
    u = stats.mannwhitneyu(y_a[f_a == 1], y_a[f_a == 0])
    return float(u.statistic / (n1 * n0)), float(u.pvalue)


def fast_strat_boot_auc(y, f, rng, n_boot=10000):
    a, p = auc_mw(y, f)
    yv, fv = np.asarray(y, float), np.asarray(f, int)
    m = np.isfinite(yv)
    yv, fv = yv[m], fv[m]
    i1, i0 = np.where(fv == 1)[0], np.where(fv == 0)[0]
    n1, n0 = len(i1), len(i0)
    if n1 < 1 or n0 < 1 or np.nanstd(yv) == 0:
        return a, a, a, p
    K = ((yv[i1, None] > yv[None, i0]).astype(np.float64) + 0.5 * (yv[i1, None] == yv[None, i0])) / (n1 * n0)
    w1 = rng.multinomial(n1, np.ones(n1) / n1, size=n_boot).astype(np.float64)
    w0 = rng.multinomial(n0, np.ones(n0) / n0, size=n_boot).astype(np.float64)
    bs = ((w1 @ K) * w0).sum(axis=1)
    return a, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), p


def fast_boot_rho(x, y, rng, n_boot=5000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4 or np.std(x) == 0 or np.std(y) == 0:
        return 0.0, 0.0, 0.0, 1.0
    r0, p0 = stats.spearmanr(x, y)
    n = len(x)
    idx = rng.integers(0, n, size=(n_boot, n))
    xb, yb = x[idx], y[idx]
    rx = stats.rankdata(xb, axis=1) - (n + 1) / 2.0
    ry = stats.rankdata(yb, axis=1) - (n + 1) / 2.0
    num = (rx * ry).sum(axis=1)
    den = np.sqrt((rx * rx).sum(axis=1) * (ry * ry).sum(axis=1))
    bs = np.where(den > 0, num / np.maximum(den, 1e-30), 0.0)
    return float(r0), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), float(p0)


def p_het_fisher_z(r1, n1, r2, n2):
    r1c, r2c = np.clip(r1, -0.999, 0.999), np.clip(r2, -0.999, 0.999)
    z1, z2 = np.arctanh(r1c), np.arctanh(r2c)
    se = np.sqrt(1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0))
    return float(2.0 * stats.norm.sf(abs(z1 - z2) / se))


def p_het_auc_boot(y1, f1, y2, f2, rng, n_boot=5000):
    a1, _ = auc_mw(y1, f1)
    a2, _ = auc_mw(y2, f2)
    def boot_se(yv, fv):
        yv, fv = np.asarray(yv, float), np.asarray(fv, int)
        i1, i0 = np.where(fv == 1)[0], np.where(fv == 0)[0]
        n1, n0 = len(i1), len(i0)
        K = ((yv[i1, None] > yv[None, i0]).astype(np.float64) + 0.5 * (yv[i1, None] == yv[None, i0])) / (n1 * n0)
        w1 = rng.multinomial(n1, np.ones(n1) / n1, size=n_boot).astype(np.float64)
        w0 = rng.multinomial(n0, np.ones(n0) / n0, size=n_boot).astype(np.float64)
        bs = ((w1 @ K) * w0).sum(axis=1)
        return float(np.std(bs, ddof=1))
    se1, se2 = boot_se(y1, f1), boot_se(y2, f2)
    se = np.sqrt(se1 * se1 + se2 * se2)
    return float(2.0 * stats.norm.sf(abs(a1 - a2) / max(se, 1e-12)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    out_dir = args.output_dir.resolve()
    data_out = out_dir / "data"
    data_out.mkdir(parents=True, exist_ok=True)

    read_registry = []
    def reg(p: pathlib.Path, role: str) -> pathlib.Path:
        rp = p.resolve()
        h = hashlib.sha256(rp.read_bytes()).hexdigest()
        read_registry.append({"role": role, "path": str(rp), "bytes": rp.stat().st_size, "sha256": h})
        return rp

    # ------------------------------------------------------------------
    # 1. Load Core Metadata & Lock TRUE PUBLISHED 11-Death Ground Truth
    # ------------------------------------------------------------------
    S_meta = pd.read_csv(reg(ESKD_DIR / "results" / "sample_level_composition_metadata.csv", "eskd_sample_metadata"))
    Fm_raw = pd.read_csv(reg(ESKD_DIR / "results" / "sample_level_fine_label_pct.csv", "eskd_sample_fine_pct"))
    reconc = pd.read_csv(reg(ROOT / "paper_claim_analysis" / "covid_sample_patient_id_reconciliation.csv", "real_clinical_reconciliation"))
    s1_pub = pd.read_excel(reg(ROOT / "paper_claim_analysis" / "Stephenson_Table_S1_sample_metadata.xlsx", "stephenson2025_table_s1"), header=2)
    w1_meta = pd.read_csv(reg(BULK_DIR / "w1_metadata.csv", "gisby2022_w1_metadata"))
    w2_meta = pd.read_csv(reg(BULK_DIR / "w2_metadata.csv", "gisby2022_w2_metadata"))

    ster_map = reconc.set_index("sample_id")["published_received_steroid"].to_dict()
    S_meta["recorded_steroid"] = S_meta.apply(
        lambda r: str(ster_map.get(str(r["sample_id"]), "NA")) if r["source"] == "Real" else "UNRECORDED_W1",
        axis=1,
    )
    S_meta["steroid_real_Y"] = S_meta["recorded_steroid"].eq("Y")
    S_meta["fatal_pub11"] = S_meta["patient_id"].isin(D_PUB11).astype(int)
    S_meta["fatal_dir9"] = S_meta["patient_id"].isin(D_DIR9).astype(int)
    S_meta["fatal_zen14"] = S_meta["patient_id"].isin(D_ZEN14).astype(int)
    S_meta.to_csv(data_out / "eskd_sample_metadata_verified.csv", index=False)

    LAB43_ALL = [c for c in Fm_raw.columns if c not in ("source", "sample_id", "RBC", "Platelet")]
    Fm = Fm_raw.copy()
    Fm[LAB43_ALL] = 100.0 * Fm[LAB43_ALL].div(Fm[LAB43_ALL].sum(axis=1), axis=0)
    S_all = S_meta.merge(Fm[["source", "sample_id"] + LAB43_ALL], on=["source", "sample_id"], how="left")

    pos_all_incl_cll = S_all[(S_all["condition"] == "positive") & S_all["in_primary_scope"].astype(bool) & S_all["qc_pass_min200cells"].astype(bool)].copy()
    pos_samps = pos_all_incl_cll[~pos_all_incl_cll["patient_id"].isin(CLL)].copy()

    real_ids = set(pos_samps[pos_samps["source"] == "Real"]["patient_id"])
    pos_r = pos_samps[pos_samps["source"] == "Real"].copy()
    pos_p_all = pos_samps[pos_samps["source"] == "PBISC"].copy()
    pos_n = pos_p_all[~pos_p_all["patient_id"].isin(real_ids)].copy()
    SHARED27 = sorted(real_ids & set(pos_p_all["patient_id"]))
    NEW29 = sorted(set(pos_n["patient_id"]))
    REAL37 = sorted(real_ids)

    # Exact 43 fine subsets (>=30% detection in Real)
    R_pt_all = pos_r.groupby("patient_id")[LAB43_ALL].mean().loc[REAL37]
    FINE43 = [f for f in LAB43_ALL if R_pt_all[f].gt(0).mean() >= 0.30]
    assert len(FINE43) == 43
    R_pt = R_pt_all[FINE43].copy()
    N_pt = pos_n.groupby("patient_id")[FINE43].mean().loc[NEW29]
    PB_sh27_pt = pos_p_all[pos_p_all["patient_id"].isin(SHARED27)].groupby("patient_id")[FINE43].mean().loc[SHARED27]
    RARE23 = [f for f in FINE43 if (R_pt[f].mean() < 1.0 and R_pt[f].median() < 1.0)]

    # Patient-level raw % and within-source Normal score tables (66 non-CLL patients, TRUE_PUB_11!)
    pt_cl = pos_samps.groupby("patient_id").agg(sev=("severity_ord", "first"), WHO_severity=("WHO_severity", "first"))
    pt_cl["fat11"] = pt_cl.index.isin(D_PUB11).astype(int)
    pt_cl["fat9"] = pt_cl.index.isin(D_DIR9).astype(int)
    pt_cl["fat14"] = pt_cl.index.isin(D_ZEN14).astype(int)
    pt_cl["fat13_legacy"] = pt_cl.index.isin(D_LEG13).astype(int)

    MONO_COLS = [c for c in ["CD14mono", "CD14mono_IFN", "CD14mono_dex", "CD16mono", "CD16mono_C1", "CD16mono_IFN", "Int.mono"] if c in FINE43]
    PT_RAW = pd.concat([R_pt.assign(src="Real"), N_pt.assign(src="New")])
    PT_RAW["CD16mono|lin"] = 100.0 * PT_RAW["CD16mono"] / PT_RAW[MONO_COLS].sum(axis=1).replace(0, np.nan)
    PT_RAW = PT_RAW.join(pt_cl[["sev", "fat11", "fat9", "fat14", "fat13_legacy"]])
    PT_RAW.to_csv(data_out / "fig7_patient_raw_pct_noCLL.csv")

    PT_NS = pd.concat([
        pd.DataFrame({f: rint(R_pt[f]) for f in FINE43}, index=REAL37).assign(src="Real"),
        pd.DataFrame({f: rint(N_pt[f]) for f in FINE43}, index=NEW29).assign(src="New"),
    ])
    PT_NS = PT_NS.join(pt_cl[["sev", "fat11", "fat9", "fat14", "fat13_legacy"]])
    PT_NS.to_csv(data_out / "fig7_patient_normal_scores_noCLL.csv")

    # 217 longitudinal draws table (with verified TRUE_PUB_11!)
    for fk in ["HSPC", "T CD8 CTL", "T CD8 CM", "B_ASC_IgG", "B_ASC_IgA"]:
        pos_r[fk + "_ns"] = rint(pos_r[fk])
        pos_n[fk + "_ns"] = rint(pos_n[fk])
    pos_217 = pd.concat([pos_r, pos_n], ignore_index=True).sort_values(["patient_id", "days_from_first_positive_swab"])
    pos_217["days"] = pd.to_numeric(pos_217["days_from_first_positive_swab"], errors="coerce")
    pos_217[
        ["source", "sample_id", "patient_id", "days", "age", "sex", "severity_ord", "temp_severity_ord",
         "fatal_pub11", "fatal_dir9", "fatal_zen14", "temp_CRP", "temp_lymph",
         "HSPC", "T CD8 CTL", "HSPC_ns", "T CD8 CTL_ns"]
    ].to_csv(data_out / "eskd_longitudinal_217_draws.csv", index=False)

    # Copy verified static reference tables (CLL bulk screen, public cohorts, Real B transfer P/R, NK decomp)
    for src_p, dst_name, role in [
        (D15_DIR / "cll_molecular_screen_bulk_donors.csv", "fig7_cll_molecular_screen.csv", "fig7_cll_screen"),
        (D15_DIR / "modality_scale_shared27_and_raw_pct_v8.csv", "suppS2_modality_scale_summary.csv", "s2_mod_scale"),
        (D15_DIR / "shared27_paired_pct_v8.csv", "suppS2_shared27_paired_pct.csv", "s2_sh27_pct"),
        (D15_DIR / "nk_decomposition_v8.csv", "suppS2_nk_decomposition.csv", "s2_nk_decomp"),
        (D21_DIR / "public_corrected_cohort_endpoint_counts.csv", "suppS4_public_endpoint_counts.csv", "s4_pub_cnts"),
        (D21_DIR / "public_corrected_cohort_results.csv", "suppS4_public_cohort_reconciled.csv", "s4_pub_coh"),
        (D21_DIR / "public_corrected_meta_results.csv", "suppS4_public_meta_reconciled.csv", "s4_pub_meta"),
        (D21_DIR / "public_corrected_hspc_windows.csv", "suppS4_public_hspc_windows_reconciled.csv", "s4_pub_win"),
        (D21_DIR / "S1_annotation_precision_recall.csv", "suppS2_realB_transfer_precision_recall.csv", "s2_realB_pr"),
    ]:
        pd.read_csv(reg(src_p, role)).to_csv(data_out / dst_name, index=False)

    # ------------------------------------------------------------------
    # 2. Point 1 & 2: Verified Mortality Forensic Discrepancy Table & Patient Master Tables
    # ------------------------------------------------------------------
    p68_prev = pd.read_csv(reg(D24_DIR / "Patient_Metadata_and_Dataset_Inventory_20261119.csv", "prev_patient68_inv"))
    s362_prev = pd.read_csv(reg(D24_DIR / "Sample_1to1_Mapping_All_3Datasets_20261119.csv", "prev_sample362_map"))
    u107_prev = pd.read_csv(reg(D24_DIR / "All_107_ESKD_Donors_Dataset_Availability_20261119.csv", "prev_donor107_inv"))

    def true_death_label(pid):
        if pid in D_DIR9:
            return "Died (Direct COVID)"
        if pid in D_PUB11:
            return "Died (Acute indirect, Table S1)"
        if pid in {"C101", "C146", "C147"}:
            return "Survived (Table S1=no; Gisby Fig.1 survivor; w2 bulk metadata bug=TRUE)"
        return "Survived"

    p68_ver = p68_prev.copy()
    p68_ver["verified_outcome_pub11"] = p68_ver["patient_id"].map(true_death_label)
    p68_ver["fatal_pub11_PRIMARY"] = p68_ver["patient_id"].isin(D_PUB11).astype(int)
    p68_ver["fatal_direct9"] = p68_ver["patient_id"].isin(D_DIR9).astype(int)
    p68_ver["fatal_zenodo_w1w2_raw14"] = p68_ver["patient_id"].isin(D_ZEN14).astype(int)
    p68_ver["fatal_legacy13_hybrid_ERRONEOUS"] = p68_ver["patient_id"].isin(D_LEG13).astype(int)
    p68_ver.drop(columns=["outcome_status", "fatal_all", "fatal_acute_S1", "fatal_direct"], inplace=True)
    p68_ver.to_csv(data_out / "table_s1_patients_68_verified.csv", index=False)

    s362_ver = s362_prev.copy()
    s362_ver["outcome_status"] = s362_ver["patient_id"].map(lambda p: true_death_label(p) if p in set(p68_ver["patient_id"]) else "Control (COVID-)")
    s362_ver["fatal_pub11_PRIMARY"] = s362_ver["patient_id"].isin(D_PUB11).astype(int)
    s362_ver.to_csv(data_out / "table_s2_sample_1to1_mapping_362.csv", index=False)
    u107_prev.to_csv(data_out / "table_s3_all_107_donors_overview.csv", index=False)

    # Mortality discrepancy forensic audit table
    disc_rows = []
    s1_pos = s1_pub[s1_pub["case_control"] == "POSITIVE"].drop_duplicates("patient_id").set_index("patient_id")
    w1_pos = w1_meta[w1_meta["case_control"] == "POSITIVE"].drop_duplicates("individual_id").set_index("individual_id")
    w2_pos = w2_meta[w2_meta["case_control"] == "POSITIVE"].drop_duplicates("individual_id").set_index("individual_id")
    w2_rec_ids = set(w2_meta[w2_meta["case_control"] == "RECOVERY"]["individual_id"])
    for pid in sorted(D_ZEN14, key=lambda x: int(x.lstrip("C"))):
        in_r = pid in real_ids
        in_p = pid in set(pos_p_all["patient_id"])
        disc_rows.append(dict(
            patient_id=pid,
            cohort=("Shared (Real ∩ PBISC)" if (in_r and in_p) else ("Real-only" if in_r else "PBISC new")),
            WHO_severity=pt_cl.loc[pid, "WHO_severity"],
            stephenson2025_TableS1_fatal=(str(s1_pos.loc[pid, "fatal_disease"]) if pid in s1_pos.index else "Not in CITE-seq"),
            gisby2022_w1_meta_fatal=(str(w1_pos.loc[pid, "fatal_disease"]) if pid in w1_pos.index else "-"),
            gisby2022_w2_meta_fatal=(str(w2_pos.loc[pid, "fatal_disease"]) if pid in w2_pos.index else "-"),
            has_2month_recovery_sample=(pid in w2_rec_ids),
            TRUE_PUBLISHED_11_PRIMARY=int(pid in D_PUB11),
            direct_COVID_9=int(pid in D_DIR9),
            legacy_13_hybrid=int(pid in D_LEG13),
            zenodo_raw_14=int(pid in D_ZEN14),
            forensic_verdict=(
                "VERIFIED COVID-19 DEATH (Direct)" if pid in D_DIR9 else
                "VERIFIED COVID-19 DEATH (Indirect acute, Table S1 & Gisby Fig.1)" if pid in D_PUB11 else
                "SURVIVOR MISCODED IN w2_metadata.csv (Table S1=no; Gisby Fig.1=13 survivors; donated 2m recovery sample)" if pid in w2_rec_ids else
                "SURVIVOR MISCODED IN w2_metadata.csv (Table S1=no; Gisby Fig.1=13th survivor missing recovery due to logistics)"
            ),
        ))
    df_mort_audit = pd.DataFrame(disc_rows)
    df_mort_audit.to_csv(data_out / "metadata_mortality_audit_discrepancies.csv", index=False)
    print("[Point 1] Verified TRUE_PUB_11 deaths:", len(D_PUB11), "(Real=5, PBISC new=6); fixed C101, C146, C147!")

    # Write formatted multi-sheet Excel Table S1 workbook
    wb = openpyxl.Workbook(); wb.remove(wb.active)
    hdr_fill = PatternFill("solid", fgColor="0F4D92")
    hdr_font = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
    body_font = Font(name="Calibri", size=10, color="1F2937")
    bold_font = Font(name="Calibri", size=10, bold=True, color="0F172A")
    border = Border(*[Side(style="thin", color="D1D5DB")] * 4)
    for sh_nm, df_sh in [
        ("Mortality_Audit_Point1", df_mort_audit),
        ("Patients_68_Verified", p68_ver),
        ("Sample_1to1_Mapping_362", s362_ver),
        ("All_107_Donors", u107_prev),
    ]:
        ws = wb.create_sheet(title=sh_nm)
        cols = list(df_sh.columns)
        for ci, cn in enumerate(cols, 1):
            c = ws.cell(row=1, column=ci, value=cn)
            c.fill, c.font, c.border = hdr_fill, hdr_font, border
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        for ri, (_, r) in enumerate(df_sh.iterrows(), 2):
            for ci, cn in enumerate(cols, 1):
                v = r[cn]
                if pd.isna(v): v = ""
                elif isinstance(v, (np.bool_, bool)): v = "TRUE" if v else "FALSE"
                elif isinstance(v, (np.integer,)): v = int(v)
                elif isinstance(v, (np.floating,)): v = int(v) if float(v).is_integer() else round(float(v), 4)
                c = ws.cell(row=ri, column=ci, value=v)
                c.font = bold_font if ci == 1 else body_font
                c.border = border
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{len(df_sh) + 1}"
        for ci, cn in enumerate(cols, 1):
            ml = max([len(str(cn))] + [len(str(ws.cell(row=r, column=ci).value or "")) for r in range(2, min(len(df_sh) + 2, 100))])
            ws.column_dimensions[get_column_letter(ci)].width = min(max(ml + 3, 11), 46)
    wb.save(out_dir / "Table_S1_Verified_Patient_and_Sample_Metadata_20261119.xlsx")

    # ------------------------------------------------------------------
    # 3. Point 3: Full 43-Subset Scan, 5-Step Funnel, Scorecard & Stratified CIs under TRUE_PUB_11
    # ------------------------------------------------------------------
    rng_main = np.random.default_rng(20261119)
    mr = (PT_NS["src"] == "Real").values
    mn = (PT_NS["src"] == "New").values
    sev_r, sev_n, sev_c = PT_NS.loc[mr, "sev"].values, PT_NS.loc[mn, "sev"].values, PT_NS["sev"].values
    fat_r, fat_n, fat_c = PT_NS.loc[mr, "fat11"].values, PT_NS.loc[mn, "fat11"].values, PT_NS["fat11"].values

    scan_rows = []
    for cll_excl in [True, False]:
        if cll_excl:
            r_df, n_df = R_pt, N_pt
            sr_, sn_ = sev_r, sev_n
            fr_, fn_ = fat_r, fat_n
        else:
            r_all_c = pos_all_incl_cll[pos_all_incl_cll["source"] == "Real"].groupby("patient_id")[FINE43 + ["severity_ord"]].mean()
            r_all_c["fat11"] = r_all_c.index.isin(D_PUB11).astype(int)
            p_all_c = pos_all_incl_cll[pos_all_incl_cll["source"] == "PBISC"].groupby("patient_id")[FINE43 + ["severity_ord"]].mean()
            n_all_c = p_all_c.loc[~p_all_c.index.isin(r_all_c.index)].copy()
            n_all_c["fat11"] = n_all_c.index.isin(D_PUB11).astype(int)
            r_df, n_df = r_all_c[FINE43], n_all_c[FINE43]
            sr_, sn_ = r_all_c["severity_ord"].values, n_all_c["severity_ord"].values
            fr_, fn_ = r_all_c["fat11"].values, n_all_c["fat11"].values
        sc_, fc_ = np.r_[sr_, sn_], np.r_[fr_, fn_]
        c_ns_df = pd.concat([
            pd.DataFrame({f: rint(r_df[f]) for f in FINE43}, index=r_df.index),
            pd.DataFrame({f: rint(n_df[f]) for f in FINE43}, index=n_df.index),
        ])
        for sname, df_s, sv, fv in [("Real", r_df, sr_, fr_), ("PBISC new", n_df, sn_, fn_), ("Real + PBISC new", c_ns_df, sc_, fc_)]:
            sub_r = []
            for fk in FINE43:
                y = df_s[fk].values
                rho_v, ps_v = (stats.spearmanr(sv, y) if np.std(y) > 0 else (0.0, 1.0))
                auc_v, pd_v = auc_mw(y, fv)
                sub_r.append(dict(CLL_excluded=cll_excl, set=sname, feature=fk, rho=float(rho_v), p_sev=float(ps_v), AUC_death=float(auc_v), p_death=float(pd_v)))
            ds = pd.DataFrame(sub_r)
            # Fixed BH43 across all 43 subsets (constant/non-estimable kept at P=1)
            ds["q_sev"] = multipletests(ds["p_sev"], method="fdr_bh")[1]
            ds["q_death"] = multipletests(ds["p_death"], method="fdr_bh")[1]
            scan_rows.append(ds)
    df_scan = pd.concat(scan_rows, ignore_index=True)

    # Add P_het (Fisher-z for severity rho & bootstrap-SE z for death AUC) to Combined noCLL
    sc_no = df_scan[df_scan["CLL_excluded"]].copy()
    r43 = sc_no[sc_no["set"] == "Real"].set_index("feature")
    n43 = sc_no[sc_no["set"] == "PBISC new"].set_index("feature")
    c43 = sc_no[sc_no["set"] == "Real + PBISC new"].set_index("feature")
    phet_s, phet_d = {}, {}
    for fk in FINE43:
        phet_s[fk] = p_het_fisher_z(r43.loc[fk, "rho"], len(REAL37), n43.loc[fk, "rho"], len(NEW29))
        phet_d[fk] = p_het_auc_boot(R_pt[fk].values, fat_r, N_pt[fk].values, fat_n, rng_main, n_boot=3000)
    df_scan["p_het_sev"] = df_scan["feature"].map(phet_s)
    df_scan["p_het_death"] = df_scan["feature"].map(phet_d)
    df_scan.to_csv(data_out / "fig7_scan_43fine_true11.csv", index=False)

    # All-43 threshold transitions (Real -> Combined)
    trans_rows = []
    for fk in FINE43:
        trans_rows.append(dict(
            feature=fk, is_rare23=(fk in RARE23),
            real_rho=r43.loc[fk, "rho"], real_p_sev=r43.loc[fk, "p_sev"], real_q_sev=r43.loc[fk, "q_sev"],
            comb_rho=c43.loc[fk, "rho"], comb_p_sev=c43.loc[fk, "p_sev"], comb_q_sev=c43.loc[fk, "q_sev"],
            sev_q05_status=("retained" if (r43.loc[fk, "q_sev"] < 0.05 and c43.loc[fk, "q_sev"] < 0.05) else
                            "gained" if (r43.loc[fk, "q_sev"] >= 0.05 and c43.loc[fk, "q_sev"] < 0.05) else
                            "lost" if (r43.loc[fk, "q_sev"] < 0.05 and c43.loc[fk, "q_sev"] >= 0.05) else "neither"),
            real_auc11=r43.loc[fk, "AUC_death"], real_p_dth11=r43.loc[fk, "p_death"], real_q_dth11=r43.loc[fk, "q_death"],
            comb_auc11=c43.loc[fk, "AUC_death"], comb_p_dth11=c43.loc[fk, "p_death"], comb_q_dth11=c43.loc[fk, "q_death"],
            dth11_q10_status=("retained" if (r43.loc[fk, "q_death"] < 0.10 and c43.loc[fk, "q_death"] < 0.10) else
                              "gained" if (r43.loc[fk, "q_death"] >= 0.10 and c43.loc[fk, "q_death"] < 0.10) else
                              "lost" if (r43.loc[fk, "q_death"] < 0.10 and c43.loc[fk, "q_death"] >= 0.10) else "neither"),
        ))
    df_trans = pd.DataFrame(trans_rows)
    df_trans.to_csv(data_out / "fig7_all43_threshold_transitions_true11.csv", index=False)

    # 5-step selection funnel under TRUE_PUB_11
    sel_rows = []
    for fk in FINE43:
        s2_s = c43.loc[fk, "q_sev"] < 0.05
        s2_d = c43.loc[fk, "q_death"] < 0.10
        s3_s = bool(s2_s and (np.sign(r43.loc[fk, "rho"]) == np.sign(n43.loc[fk, "rho"])) and (abs(r43.loc[fk, "rho"]) > 0.03))
        s3_d = bool(s2_d and (np.sign(r43.loc[fk, "AUC_death"] - 0.5) == np.sign(n43.loc[fk, "AUC_death"] - 0.5)))
        s4_s = bool(s3_s and (fk not in {"CD14mono_dex", "CD14mono_IFN", "MAIT"}))
        s4_d = bool(s3_d and (fk not in {"CD14mono_dex", "CD14mono_IFN", "MAIT"}))
        s5_s = bool(s4_s and (r43.loc[fk, "q_sev"] >= 0.05))
        s5_d = bool(s4_d and (r43.loc[fk, "q_death"] >= 0.10))
        anc_s = bool(s4_s and (r43.loc[fk, "q_sev"] < 0.05))
        sel_rows.append(dict(
            feature=fk, q_sev_Real=r43.loc[fk, "q_sev"], q_sev_Ext=c43.loc[fk, "q_sev"],
            q_death_Real=r43.loc[fk, "q_death"], q_death_Ext=c43.loc[fk, "q_death"],
            S2_sev=s2_s, S2_death=s2_d, S3_sev=s3_s, S3_death=s3_d,
            S4_sev=s4_s, S4_death=s4_d, S5_sev=s5_s, S5_death=s5_d, Anchor_sev=anc_s,
        ))
    df_sel = pd.DataFrame(sel_rows)
    df_sel.to_csv(data_out / "fig7_selection_funnel_43fine_true11.csv", index=False)

    # 14-variant sensitivity matrix under TRUE_PUB_11 (recomputed from sample-level data!)
    SENS_FEATS = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN", "NK", "CD16mono_C1"]
    LIN_MAP = {}
    for fk in SENS_FEATS:
        if fk.startswith("B_ASC"): LIN_MAP[fk] = [c for c in FINE43 if c.startswith("B_")]
        elif fk.startswith("T "): LIN_MAP[fk] = [c for c in FINE43 if c.startswith("T ") or c in ("MAIT", "NKT")]
        elif fk.startswith("NK"): LIN_MAP[fk] = [c for c in FINE43 if c.startswith("NK")]
        elif "mono" in fk: LIN_MAP[fk] = MONO_COLS
        else: LIN_MAP[fk] = FINE43
    NON_MYEL = [c for c in FINE43 if c not in MONO_COLS and not c.startswith("DC") and c not in ("pDC", "ASDC")]

    def eval_variant(df_samps, fat_col="fatal_pub11", use_lin=False, use_nonmyel=False, adj_cols=None, first_only=False):
        ds = df_samps.copy()
        if first_only:
            ds = ds.sort_values("days_from_first_positive_swab").groupby(["source", "patient_id"], as_index=False).first()
        r_sub = ds[ds["source"] == "Real"].copy()
        n_sub = ds[(ds["source"] == "PBISC") & (~ds["patient_id"].isin(set(r_sub["patient_id"])))].copy()
        res_v = {}
        for fk in SENS_FEATS:
            if use_nonmyel and fk not in NON_MYEL:
                res_v[fk] = (np.nan, np.nan, np.nan, np.nan); continue
            if use_lin and fk != "HSPC":
                denom_cols = LIN_MAP[fk]
                r_val = (100.0 * r_sub[fk] / r_sub[denom_cols].sum(axis=1).replace(0, np.nan)).groupby(r_sub["patient_id"]).mean()
                n_val = (100.0 * n_sub[fk] / n_sub[denom_cols].sum(axis=1).replace(0, np.nan)).groupby(n_sub["patient_id"]).mean()
            elif use_nonmyel:
                r_val = (100.0 * r_sub[fk] / r_sub[NON_MYEL].sum(axis=1).replace(0, np.nan)).groupby(r_sub["patient_id"]).mean()
                n_val = (100.0 * n_sub[fk] / n_sub[NON_MYEL].sum(axis=1).replace(0, np.nan)).groupby(n_sub["patient_id"]).mean()
            else:
                r_val = r_sub.groupby("patient_id")[fk].mean()
                n_val = n_sub.groupby("patient_id")[fk].mean()
            y_comb = pd.concat([rint(r_val), rint(n_val)]).dropna()
            meta_p = pd.concat([
                r_sub.groupby("patient_id").agg(sev=("severity_ord","first"), fat=(fat_col,"first"), age=("age","first"), sex=("sex","first"), cw=("centre_or_wave","first"), dex=("CD14mono_dex","mean")),
                n_sub.groupby("patient_id").agg(sev=("severity_ord","first"), fat=(fat_col,"first"), age=("age","first"), sex=("sex","first"), cw=("centre_or_wave","first"), dex=("CD14mono_dex","mean")),
            ]).loc[y_comb.index]
            y_use = y_comb.values.copy()
            if adj_cols:
                X_a = pd.DataFrame(index=y_comb.index)
                if "cw" in adj_cols:
                    X_a = pd.concat([X_a, pd.get_dummies(meta_p["cw"], drop_first=True, dtype=float)], axis=1)
                if "age" in adj_cols: X_a["age"] = meta_p["age"].astype(float)
                if "sex" in adj_cols: X_a["sex"] = (meta_p["sex"].astype(str).str.upper().str.startswith("M")).astype(float)
                if "dex" in adj_cols: X_a["dex"] = meta_p["dex"].astype(float)
                X_a = sm.add_constant(X_a.fillna(X_a.mean()))
                y_use = sm.OLS(y_use, X_a).fit().resid.values
            rv, ps = stats.spearmanr(meta_p["sev"].values, y_use)
            av, pd_ = auc_mw(y_use, meta_p["fat"].values)
            res_v[fk] = (float(rv), float(ps), float(av), float(pd_))
        return res_v

    pos_mono_ok = pos_samps[~((pos_samps["source"] == "Real") & (pos_samps["Mono_pct"] < 2.0))]
    pos_no_out = pos_samps[~pos_samps["patient_id"].isin({"C187", "C190"})]
    pos_ster_no = pos_samps[~((pos_samps["source"] == "Real") & pos_samps["steroid_real_Y"])]
    pos_d20 = pos_samps[pd.to_numeric(pos_samps["days_from_first_positive_swab"], errors="coerce") <= 20]

    VAR_SPECS = [
        ("0 baseline (Published 11d)", pos_samps, "fatal_pub11", False, False, None, False),
        ("A1 adjust centre / wave", pos_samps, "fatal_pub11", False, False, ["cw"], False),
        ("A2 % of own lineage", pos_samps, "fatal_pub11", True, False, None, False),
        ("A3 % of non-myeloid cells", pos_samps, "fatal_pub11", False, True, None, False),
        ("A4 drop Real mono <2%", pos_mono_ok, "fatal_pub11", False, False, None, False),
        ("B  drop C187 & C190 outliers", pos_no_out, "fatal_pub11", False, False, None, False),
        ("C1 death = Direct COVID only (9d)", pos_samps, "fatal_dir9", False, False, None, False),
        ("C2 death = Zenodo w1+w2 raw (14d)", pos_samps, "fatal_zen14", False, False, None, False),
        ("D1 Real recorded steroid-No draws", pos_ster_no, "fatal_pub11", False, False, None, False),
        ("D2 adjust CD14 mono dex (proxy)", pos_samps, "fatal_pub11", False, False, ["dex"], False),
        ("E  samples <= day 20", pos_d20, "fatal_pub11", False, False, None, False),
        ("F1 first sample per patient", pos_samps, "fatal_pub11", False, False, None, True),
        ("F2 adjust age + sex", pos_samps, "fatal_pub11", False, False, ["age", "sex"], False),
        ("G  combined: lineage %, no outliers, adj.", pos_no_out, "fatal_pub11", True, False, ["cw", "age", "sex"], False),
    ]
    sa_rows = []
    for vname, df_v, fcol, ulin, unm, adjc, fonly in VAR_SPECS:
        rv_dict = eval_variant(df_v, fcol, ulin, unm, adjc, fonly)
        for fk, (rv, ps, av, pd_) in rv_dict.items():
            sa_rows.append(dict(variant=vname, feature=fk, rho=rv, p_sev=ps, AUC=av, p_death=pd_))
    df_sa = pd.DataFrame(sa_rows)
    df_sa.to_csv(data_out / "suppS3_sensitivity_14variants_true11.csv", index=False)

    # Credibility Scorecard (Fig. 7h) for all 8 Step-5 PBISC-enabled subset-endpoint pairs!
    SC_SPECS = [
        ("B_ASC_IgG", "sev", True),
        ("B_ASC_IgA", "sev", True),
        ("HSPC", "sev", True),
        ("HSPC", "death", True),
        ("T CD8 CTL", "sev", True),
        ("T CD8 CTL", "death", True),
        ("NK", "sev", False),
        ("CD16mono_C1", "death", False),
    ]
    pt_w23 = PT_NS[PT_NS["sev"] >= 2]
    sc_rows = []
    for fk, out_k, fol in SC_SPECS:
        sh_r, sh_p = stats.spearmanr(R_pt.loc[SHARED27, fk], PB_sh27_pt.loc[SHARED27, fk])
        sa_f = df_sa[df_sa["feature"] == fk]
        if out_k == "sev":
            r_eff, n_eff, ph = r43.loc[fk, "rho"], n43.loc[fk, "rho"], phet_s[fk]
            rob_p = int((sa_f["p_sev"] < 0.05).sum())
            rob_t = int(sa_f["p_sev"].notna().sum())
            w23_a, w23_p = np.nan, np.nan
        else:
            r_eff, n_eff, ph = r43.loc[fk, "AUC_death"], n43.loc[fk, "AUC_death"], phet_d[fk]
            rob_p = int((sa_f["p_death"] < 0.05).sum())
            rob_t = int(sa_f["p_death"].notna().sum())
            w23_a, w23_p = auc_mw(pt_w23[fk].values, pt_w23["fat11"].values)
        sc_rows.append(dict(
            feature=fk, outcome=out_k, followed=fol,
            effect_Real=r_eff, effect_New=n_eff, p_het=ph,
            shared_rho=float(sh_r), shared_p=float(sh_p),
            robust_pass=rob_p, robust_total=rob_t,
            AUC_within_WHO23=w23_a, p_within_WHO23=w23_p,
        ))
    df_sc = pd.DataFrame(sc_rows)
    df_sc.to_csv(data_out / "fig7_credibility_scorecard_true11.csv", index=False)

    # 3,000-resample source-stratified per-draw INT-refit CIs for Fig. 7i & S2b, c (TRUE_PUB_11!)
    FEATS_CI = ["NK CD56 bright", "B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC"]
    rng_ci = np.random.default_rng(20261119)
    ci7_rows = []
    for cll_excl in [False, True]:
        sub_all = pos_samps if cll_excl else pos_all_incl_cll
        r_p = sub_all[sub_all["source"] == "Real"].groupby("patient_id")[FEATS_CI + ["severity_ord"]].mean()
        r_p["fat11"] = r_p.index.isin(D_PUB11).astype(int)
        b_p = sub_all[sub_all["source"] == "PBISC"].groupby("patient_id")[FEATS_CI + ["severity_ord"]].mean()
        n_p = b_p.loc[~b_p.index.isin(r_p.index)].copy()
        n_p["fat11"] = n_p.index.isin(D_PUB11).astype(int)
        nr, nn = len(r_p), len(n_p)
        Xr, sr, fr = r_p[FEATS_CI].values, r_p["severity_ord"].values, r_p["fat11"].values.astype(int)
        Xn, sn, fn = n_p[FEATS_CI].values, n_p["severity_ord"].values, n_p["fat11"].values.astype(int)
        Zr0 = stats.norm.ppf((stats.rankdata(Xr, axis=0) - 0.5) / nr)
        Zn0 = stats.norm.ppf((stats.rankdata(Xn, axis=0) - 0.5) / nn)
        Xc0, sc0, fc0 = np.vstack([Zr0, Zn0]), np.r_[sr, sn], np.r_[fr, fn]
        rr0, psr0 = spearman_vec_cols(sr, Xr); ar0, pdr0 = mw_vec_cols(fr, Xr)
        rn0, psn0 = spearman_vec_cols(sn, Xn); an0, pdn0 = mw_vec_cols(fn, Xn)
        rc0, psc0 = spearman_vec_cols(sc0, Xc0); ac0, pdc0 = mw_vec_cols(fc0, Xc0)

        idx_r = rng_ci.integers(0, nr, size=(3000, nr))
        idx_n = rng_ci.integers(0, nn, size=(3000, nn))
        bs_rr, bs_ar = [[] for _ in FEATS_CI], [[] for _ in FEATS_CI]
        bs_rn, bs_an = [[] for _ in FEATS_CI], [[] for _ in FEATS_CI]
        bs_rc, bs_ac = [[] for _ in FEATS_CI], [[] for _ in FEATS_CI]
        for b in range(3000):
            ir, in_ = idx_r[b], idx_n[b]
            Xrb, srb, frb = Xr[ir], sr[ir], fr[ir]
            Xnb, snb, fnb = Xn[in_], sn[in_], fn[in_]
            if len(np.unique(srb)) > 1:
                v, _ = spearman_vec_cols(srb, Xrb)
                for j in range(len(FEATS_CI)):
                    if np.std(Xrb[:, j]) > 0: bs_rr[j].append(v[j])
            if 0 < frb.sum() < nr:
                v, _ = mw_vec_cols(frb, Xrb)
                for j in range(len(FEATS_CI)):
                    if np.std(Xrb[:, j]) > 0: bs_ar[j].append(v[j])
            if len(np.unique(snb)) > 1:
                v, _ = spearman_vec_cols(snb, Xnb)
                for j in range(len(FEATS_CI)):
                    if np.std(Xnb[:, j]) > 0: bs_rn[j].append(v[j])
            if 0 < fnb.sum() < nn:
                v, _ = mw_vec_cols(fnb, Xnb)
                for j in range(len(FEATS_CI)):
                    if np.std(Xnb[:, j]) > 0: bs_an[j].append(v[j])
            Zrb = stats.norm.ppf((stats.rankdata(Xrb, axis=0) - 0.5) / nr); Zrb[:, np.std(Xrb, axis=0) == 0] = 0.0
            Znb = stats.norm.ppf((stats.rankdata(Xnb, axis=0) - 0.5) / nn); Znb[:, np.std(Xnb, axis=0) == 0] = 0.0
            Xcb, scb, fcb = np.vstack([Zrb, Znb]), np.r_[srb, snb], np.r_[frb, fnb]
            if len(np.unique(scb)) > 1:
                v, _ = spearman_vec_cols(scb, Xcb)
                for j in range(len(FEATS_CI)):
                    if np.std(Xcb[:, j]) > 0: bs_rc[j].append(v[j])
            if 0 < fcb.sum() < (nr + nn):
                v, _ = mw_vec_cols(fcb, Xcb)
                for j in range(len(FEATS_CI)):
                    if np.std(Xcb[:, j]) > 0: bs_ac[j].append(v[j])
        for j, fk in enumerate(FEATS_CI):
            for snm, re_, bsr, pse, ae_, bsa, pde in [
                ("Real", rr0[j], bs_rr[j], psr0[j], ar0[j], bs_ar[j], pdr0[j]),
                ("PBISC new", rn0[j], bs_rn[j], psn0[j], an0[j], bs_an[j], pdn0[j]),
                ("Real + PBISC new", rc0[j], bs_rc[j], psc0[j], ac0[j], bs_ac[j], pdc0[j]),
            ]:
                ci7_rows.append(dict(
                    CLL_excluded=cll_excl, set=snm, feature=fk,
                    rho=float(re_), rho_lo=float(np.percentile(bsr, 2.5)), rho_hi=float(np.percentile(bsr, 97.5)), p_sev=float(pse),
                    AUC=float(ae_), AUC_lo=float(np.percentile(bsa, 2.5)), AUC_hi=float(np.percentile(bsa, 97.5)), p_death=float(pde),
                ))
    df_ci7 = pd.DataFrame(ci7_rows)
    df_ci7.to_csv(data_out / "fig7_cohort_effect_sizes_stratified_ci_true11.csv", index=False)

    # Fig. 7p: HSPC mortality sensitivity across 9 definitions/strata under TRUE_PUB_11 (+ L1O + Exact Source x Peak WHO!)
    strata_rows = []
    def add_hspc_stratum(an_key, r_y, r_f, n_y, n_f, is_new_only=False):
        if is_new_only:
            an, pn = auc_mw(n_y, n_f)
            strata_rows.append(dict(analysis=an_key, set="PBISC new", n=len(n_y), deaths=int(np.sum(n_f)), AUC=an, p=pn))
        else:
            ar, pr = auc_mw(r_y, r_f)
            yc = np.r_[rint_arr(r_y), rint_arr(n_y)]
            fc = np.r_[r_f, n_f]
            ac, pc = auc_mw(yc, fc)
            strata_rows.append(dict(analysis=an_key, set="Real", n=len(r_y), deaths=int(np.sum(r_f)), AUC=ar, p=pr))
            strata_rows.append(dict(analysis=an_key, set="Real + PBISC new", n=len(yc), deaths=int(np.sum(fc)), AUC=ac, p=pc))

    add_hspc_stratum("Published COVID deaths (11)", R_pt["HSPC"].values, fat_r, N_pt["HSPC"].values, fat_n)
    add_hspc_stratum("Direct COVID only (9)", R_pt["HSPC"].values, PT_NS.loc[mr, "fat9"].values, N_pt["HSPC"].values, PT_NS.loc[mn, "fat9"].values)
    add_hspc_stratum("Zenodo w1+w2 raw (14)", R_pt["HSPC"].values, PT_NS.loc[mr, "fat14"].values, N_pt["HSPC"].values, PT_NS.loc[mn, "fat14"].values)
    rw23 = PT_RAW[(PT_RAW["src"] == "Real") & (PT_RAW["sev"] >= 2)]
    nw23 = PT_RAW[(PT_RAW["src"] == "New") & (PT_RAW["sev"] >= 2)]
    add_hspc_stratum("Within WHO 2-3", rw23["HSPC"].values, rw23["fat11"].values, nw23["HSPC"].values, nw23["fat11"].values)
    # Cambridge only (Real) + PBISC new
    cam_ids = sorted(set(pos_r[pos_r["centre_or_wave"] == "Cambridge"]["patient_id"]))
    add_hspc_stratum("Cambridge + PBISC new", R_pt.loc[cam_ids, "HSPC"].values, PT_RAW.loc[cam_ids, "fat11"].values, N_pt["HSPC"].values, fat_n)
    # Samples > day 7
    r_l7 = pos_r[pd.to_numeric(pos_r["days_from_first_positive_swab"], errors="coerce") > 7].groupby("patient_id")["HSPC"].mean()
    n_l7 = pos_n[pd.to_numeric(pos_n["days_from_first_positive_swab"], errors="coerce") > 7].groupby("patient_id")["HSPC"].mean()
    add_hspc_stratum("Samples > day 7", r_l7.values, pt_cl.loc[r_l7.index, "fat11"].values, n_l7.values, pt_cl.loc[n_l7.index, "fat11"].values)
    # Samples <= day 20
    r_d20 = pos_r[pd.to_numeric(pos_r["days_from_first_positive_swab"], errors="coerce") <= 20].groupby("patient_id")["HSPC"].mean()
    n_d20 = pos_n[pd.to_numeric(pos_n["days_from_first_positive_swab"], errors="coerce") <= 20].groupby("patient_id")["HSPC"].mean()
    add_hspc_stratum("Samples <= day 20", r_d20.values, pt_cl.loc[r_d20.index, "fat11"].values, n_d20.values, pt_cl.loc[n_d20.index, "fat11"].values)
    add_hspc_stratum("PBISC new alone", None, None, N_pt["HSPC"].values, fat_n, is_new_only=True)
    pd.DataFrame(strata_rows).to_csv(data_out / "fig7_hspc_death_strata_true11.csv", index=False)

    # Leave-one-death-out range (11 deaths)
    l1o_rows = []
    for pid_d in sorted(D_PUB11):
        sub_l = PT_NS.drop(index=pid_d)
        al, pl = auc_mw(sub_l["HSPC"].values, sub_l["fat11"].values)
        l1o_rows.append(dict(dropped_death=pid_d, AUC=al, p=pl))
    df_l1o = pd.DataFrame(l1o_rows)
    df_l1o.to_csv(data_out / "fig7_hspc_leave_one_death_true11.csv", index=False)

    # Exact Source x Peak WHO conditional permutation (8,400 allocations, 33 pairs)
    strata_inf = [(src, sev, g.index.tolist(), int(g["fat11"].sum()), len(g)) for (src, sev), g in PT_NS.groupby(["src", "sev"]) if 0 < int(g["fat11"].sum()) < len(g)]
    cond_rows = []
    for fk in ["HSPC", "T CD8 CTL"]:
        obs_sum, tot_p = 0.0, 0
        s_dists = []
        for _, _, pids, d, n in strata_inf:
            vals, fo = PT_NS.loc[pids, fk].values, PT_NS.loc[pids, "fat11"].values
            tot_p += d * (n - d)
            rks = stats.rankdata(vals)
            obs_sum += float(rks[fo == 1].sum() - d * (d + 1) / 2.0)
            s_dists.append(np.array([float(rks[list(cb)].sum() - d * (d + 1) / 2.0) for cb in combinations(range(n), d)]))
        dist = s_dists[0]
        for sd in s_dists[1:]:
            dist = (dist[:, None] + sd[None, :]).ravel()
        c_auc = obs_sum / tot_p
        p_2s = float(np.mean(np.abs(dist / tot_p - 0.5) >= abs(c_auc - 0.5) - 1e-12))
        cond_rows.append(dict(feature=fk, cond_AUC=c_auc, pairs=tot_p, allocations=len(dist), exact_p_2s=p_2s))
    pd.DataFrame(cond_rows).to_csv(data_out / "fig7_exact_source_peakWHO_cond_true11.csv", index=False)

    # Unified 1,000-iteration conditional re-detection simulation (n_Real = 10..35, TRUE_PUB_11!)
    ITEMS_RD = [
        ("NK CD56 bright", "sev", -1), ("B_ASC_dividing", "sev", +1), ("T CD8 CM", "sev", -1), ("CD16mono", "sev", -1),
        ("B_ASC_IgG", "sev", +1), ("B_ASC_IgA", "sev", +1), ("HSPC", "sev", +1), ("T CD8 CTL", "sev", -1),
        ("T CD8 CTL", "death", -1), ("HSPC", "death", +1),
    ]
    z_new_43 = np.column_stack([rint_arr(N_pt[fk].values) for fk in FINE43])
    R_mat_43 = R_pt[FINE43].values
    rng_rd = np.random.default_rng(2026111917)
    rd_rows = []
    for n_sub in [10, 15, 20, 25, 30, 35]:
        cnt = {(fk, ep, cond, mode): 0 for fk, ep, _ in ITEMS_RD for cond in ("real", "new29", "shuf") for mode in ("nom", "bh43")}
        valid = 0
        while valid < 1000:
            ix = rng_rd.choice(len(REAL37), size=n_sub, replace=False)
            sr_, fr_ = sev_r[ix], fat_r[ix]
            if fr_.sum() < 1 or (fr_ == 0).sum() < 1 or len(np.unique(sr_)) < 2:
                continue
            valid += 1
            pm = rng_rd.permutation(len(NEW29))
            zr = stats.norm.ppf((stats.rankdata(R_mat_43[ix], axis=0) - 0.5) / n_sub)
            zr[:, np.std(R_mat_43[ix], axis=0) == 0] = 0.0
            zc = np.vstack([zr, z_new_43])
            for clbl, zm, sv, fv in [("real", zr, sr_, fr_), ("new29", zc, np.r_[sr_, sev_n], np.r_[fr_, fat_n]), ("shuf", zc, np.r_[sr_, sev_n[pm]], np.r_[fr_, fat_n[pm]])]:
                ra, psa = spearman_vec_cols(sv, zm)
                aa, pda = mw_vec_cols(fv, zm)
                qsl = multipletests(psa, method="fdr_bh")[1]
                qdl = multipletests(pda, method="fdr_bh")[1]
                for fk, ep, sgn in ITEMS_RD:
                    j = FINE43.index(fk)
                    if ep == "sev":
                        if psa[j] < 0.05 and np.sign(ra[j]) == sgn: cnt[(fk, ep, clbl, "nom")] += 1
                        if qsl[j] < 0.05 and np.sign(ra[j]) == sgn: cnt[(fk, ep, clbl, "bh43")] += 1
                    else:
                        if pda[j] < 0.05 and np.sign(aa[j] - 0.5) == sgn: cnt[(fk, ep, clbl, "nom")] += 1
                        if qdl[j] < 0.10 and np.sign(aa[j] - 0.5) == sgn: cnt[(fk, ep, clbl, "bh43")] += 1
        for fk, ep, _ in ITEMS_RD:
            rd_rows.append(dict(
                feature=fk, outcome=ep, n_real=n_sub, repetitions=1000,
                real_nom=cnt[(fk, ep, "real", "nom")] / 1000.0, real_bh43=cnt[(fk, ep, "real", "bh43")] / 1000.0,
                plus_new29_nom=cnt[(fk, ep, "new29", "nom")] / 1000.0, plus_new29_bh43=cnt[(fk, ep, "new29", "bh43")] / 1000.0,
                plus_shuf_nom=cnt[(fk, ep, "shuf", "nom")] / 1000.0, plus_shuf_bh43=cnt[(fk, ep, "shuf", "bh43")] / 1000.0,
            ))
    df_rd = pd.DataFrame(rd_rows)
    df_rd.to_csv(data_out / "fig7_s3_unified_redetection_true11.csv", index=False)
    print("[Point 3] Unified n=20 conditional re-detection (TRUE_PUB_11):\n", df_rd[df_rd["n_real"] == 20][["feature","outcome","real_bh43","plus_new29_bh43","plus_new29_nom","plus_shuf_bh43"]].to_string(index=False))

    # Supp Fig. 4a ESKD matched severity bootstrap CIs
    b_asc_cols = [c for c in ["B_ASC_dividing", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM", "B_ASC_other"] if c in PT_RAW.columns]
    PT_RAW["B_ASC_total"] = PT_RAW[b_asc_cols].sum(axis=1)
    R_em, N_em = PT_RAW[PT_RAW["src"] == "Real"], PT_RAW[PT_RAW["src"] == "New"]
    em_rows = []
    for fk in ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono|lin", "B_ASC_total", "B_ASC_IgG", "B_ASC_IgA", "B_ASC_IgM", "T CD8 CTL", "HSPC", "NK", "CD16mono_C1"]:
        vr, vn = R_em[fk].values, N_em[fk].values
        yc = np.r_[rint_arr(vr), rint_arr(vn)]
        r0, p0 = stats.spearmanr(sev_c, yc)
        bs = []
        for _ in range(3000):
            ir = rng_ci.integers(0, len(R_em), len(R_em))
            in_ = rng_ci.integers(0, len(N_em), len(N_em))
            bs.append(stats.spearmanr(np.r_[sev_r[ir], sev_n[in_]], np.r_[rint_arr(vr[ir]), rint_arr(vn[in_])])[0])
        em_rows.append(dict(feature=fk, rho=float(r0), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), p=float(p0)))
    pd.DataFrame(em_rows).to_csv(data_out / "suppS4_eskd_matched_severity_ci.csv", index=False)

    # ------------------------------------------------------------------
    # 4. Point 4: Comprehensive 10-Method Unbiased Bulk Deconvolution Benchmark
    #    (PBISC, BayesPrism-upd, BayesPrism-init, CIBERSORTx-Smode, CIBERSORTx-Bmode,
    #     CIBERSORT-SVR, MuSiC, Calib-NNLS, DWLS, NNLS, Bulk 8-panel markers)
    # ------------------------------------------------------------------
    print("\n[Point 4] Running 10-method bulk deconvolution benchmark (incl. BayesPrism, CIBERSORTx S/B-mode, MuSiC)...")
    B_raw = pd.read_csv(reg(BULK_DIR / "htseq_counts.csv", "bulk_htseq_counts"), index_col=0)
    B_raw = B_raw[~B_raw.index.str.startswith("__")]
    B_raw.index = B_raw.index.str.split(".").str[0]
    B_raw = B_raw.groupby(level=0).sum()

    pb_h5 = ad.read_h5ad(reg(ROOT / "deliverables" / "annotated_h5ad_20261112" / "pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad", "pbisc_annotated_h5ad"))
    ens2sym = pd.Series(pb_h5.var["gene_symbol"].values, index=pb_h5.var.index)
    for ek, sv in {
        "ENSG00000211896": "IGHG1", "ENSG00000211893": "IGHG2", "ENSG00000211897": "IGHG3", "ENSG00000211892": "IGHG4",
        "ENSG00000211895": "IGHA1", "ENSG00000211890": "IGHA2", "ENSG00000211899": "IGHM",
        "ENSG00000174059": "CD34", "ENSG00000007062": "PROM1", "ENSG00000149294": "NCAM1",
    }.items():
        if ek in B_raw.index and ek not in ens2sym.index:
            ens2sym.loc[ek] = sv
    ens_keep = ens2sym[ens2sym.index.isin(B_raw.index)]
    B_cnt_sym = B_raw.loc[ens_keep.index].copy()
    B_cnt_sym.index = ens_keep.values
    B_cnt_sym = B_cnt_sym.groupby(level=0).sum()
    B_cpm_sym = B_cnt_sym.div(B_cnt_sym.sum(axis=0), axis=1) * 1e6

    z_pb = np.load(reg(ESKD_DIR / "data" / "real_pseudobulk_sample_x_fine_panel.npz", "real_pseudobulk_npz"), allow_pickle=True)
    r_genes = pd.Index(z_pb["genes"].astype(str))
    r_mat = z_pb["pb"]
    pb_df = pd.DataFrame(r_mat, columns=r_genes).T.groupby(level=0).sum().T
    grps = pd.Series(z_pb["groups"].astype(str)).str.split("|", expand=True)
    grps.columns = ["sample_id", "fine"]
    grps["ncells"] = z_pb["ncells"].astype(int)
    rm_all = S_meta[(S_meta["source"] == "Real") & (~S_meta["patient_id"].isin(CLL)) & S_meta["qc_pass_min200cells"].astype(bool)].set_index("sample_id")
    grps["patient_id"] = grps["sample_id"].map(rm_all["patient_id"])
    gv = grps.dropna(subset=["patient_id"]).copy()
    gv = gv[gv["fine"].isin(FINE43)]

    common_g = sorted(set(pb_df.columns) & set(B_cnt_sym.index))
    BP_G = [g for g in common_g if not (g.startswith("RPS") or g.startswith("RPL") or g.startswith("MT-") or g in ("MALAT1", "NEAT1"))]
    pb_df_c = pb_df.loc[gv.index, common_g]

    def build_ref_stats(allowed_pts, top_k=25):
        m = gv["patient_id"].isin(allowed_pts).values
        sub_pb, sub_g = pb_df_c.iloc[m], gv.iloc[m]
        tot_u = sub_pb.groupby(sub_g["fine"].values).sum()
        tot_nc = sub_g.groupby("fine")["ncells"].sum()
        mean_u = tot_u.div(tot_nc, axis=0).reindex(FINE43).fillna(0.0)
        phi_raw = tot_u.reindex(FINE43).fillna(0.0)[BP_G].values.astype(np.float64)
        phi_p = np.maximum(phi_raw / np.maximum(phi_raw.sum(axis=1, keepdims=True), 1.0), 1e-8)
        phi_p /= phi_p.sum(axis=1, keepdims=True)
        pt_u = sub_pb.groupby([sub_g["patient_id"].values, sub_g["fine"].values]).sum()
        pt_nc = sub_g.groupby(["patient_id", "fine"])["ncells"].sum()
        var_ct = pt_u.div(pt_nc, axis=0).groupby(level=1).var().reindex(FINE43).fillna(0.0)
        csz = tot_u.sum(axis=1).div(tot_nc).reindex(FINE43).fillna(1.0)
        lg = np.log2(mean_u + 0.01)
        mks = []
        for ct in FINE43:
            fc = (lg.loc[ct] - lg.drop(index=ct).max(axis=0))[mean_u.loc[ct] >= 0.03]
            mks.extend(list(fc.sort_values(ascending=False).index[:top_k]))
        return mean_u, phi_p, var_ct, csz, sorted(set(mks))

    def run_bayesprism(M_cnt, phi_p, csz):
        X = M_cnt.loc[BP_G].values.T.astype(np.float64)
        N, C = X.shape[0], phi_p.shape[0]
        th = np.full((N, C), 1.0 / C)
        for _ in range(200):
            rt = X / np.maximum(th @ phi_p, 1e-30)
            th_n = th * (rt @ phi_p.T)
            th_n /= np.maximum(th_n.sum(axis=1, keepdims=True), 1e-30)
            if np.max(np.abs(th_n - th)) < 1e-6: th = th_n; break
            th = th_n
        Z_cg = (th.T @ (X / np.maximum(th @ phi_p, 1e-30))) * phi_p
        psi = np.maximum((Z_cg + 10000.0 * phi_p) / np.maximum((Z_cg + 10000.0 * phi_p).sum(axis=1, keepdims=True), 1e-30), 1e-8)
        psi /= psi.sum(axis=1, keepdims=True)
        th2 = th.copy()
        for _ in range(120):
            rt = X / np.maximum(th2 @ psi, 1e-30)
            th_n = th2 * (rt @ psi.T)
            th_n /= np.maximum(th_n.sum(axis=1, keepdims=True), 1e-30)
            if np.max(np.abs(th_n - th2)) < 1e-6: th2 = th_n; break
            th2 = th_n
        cs = np.maximum(csz.values.astype(np.float64), 1.0)
        p1 = (th / cs[None, :]); p1 = 100.0 * p1 / np.maximum(p1.sum(axis=1, keepdims=True), 1e-30)
        p2 = (th2 / cs[None, :]); p2 = 100.0 * p2 / np.maximum(p2.sum(axis=1, keepdims=True), 1e-30)
        rt_f = X / np.maximum(th2 @ psi, 1e-30)
        bp_z = {}
        for ct in ["HSPC", "T CD8 CM", "T CD8 CTL"]:
            ci = FINE43.index(ct)
            Z_ng = (th2[:, ci:ci+1] * rt_f) * psi[ci:ci+1, :]
            bp_z[ct] = pd.DataFrame(np.log1p(Z_ng / np.maximum(Z_ng.sum(axis=1, keepdims=True), 1.0) * 1e4), index=M_cnt.columns, columns=BP_G)
        return pd.DataFrame(p1, index=M_cnt.columns, columns=FINE43), pd.DataFrame(p2, index=M_cnt.columns, columns=FINE43), bp_z

    def fit_opt_nusvr(A_std, b_std):
        best_w, best_rmse = None, np.inf
        for nu_v in (0.25, 0.50, 0.75):
            svr = NuSVR(nu=nu_v, kernel="linear", C=1.0, max_iter=3000)
            svr.fit(A_std, b_std)
            w = np.maximum(svr.coef_.ravel(), 0.0)
            if w.sum() > 0:
                rmse = float(np.sqrt(np.mean((b_std - A_std @ w) ** 2)))
                if rmse < best_rmse:
                    best_rmse, best_w = rmse, w
        return best_w if best_w is not None else np.zeros(A_std.shape[1])

    def run_cibersortx(M_cpm, sig_m, mks):
        A = sig_m.loc[FINE43, mks].T.values.astype(np.float64)
        B = M_cpm.loc[mks].values.astype(np.float64)
        N = B.shape[1]
        A0 = (A - A.mean(axis=0, keepdims=True)) / np.maximum(A.std(axis=0, keepdims=True), 1e-6)
        W0 = np.zeros((N, len(FINE43)))
        for i in range(N):
            b = B[:, i]
            if np.std(b) == 0: continue
            w = fit_opt_nusvr(A0, (b - b.mean()) / np.std(b))
            W0[i] = w / w.sum() if w.sum() > 0 else 0.0
        M_star = A @ W0.T
        lgB, lgM = np.log2(B + 1.0), np.log2(np.maximum(M_star, 0.0) + 1.0)
        B_bm = np.maximum(2.0 ** np.clip((lgB - lgB.mean(axis=1, keepdims=True)) * (np.maximum(lgM.std(axis=1, keepdims=True), 0.1) / np.maximum(lgB.std(axis=1, keepdims=True), 0.1)) + lgM.mean(axis=1, keepdims=True), 0.0, 22.0) - 1.0, 0.0)
        Wb = np.zeros((N, len(FINE43)))
        for i in range(N):
            b = B_bm[:, i]
            if np.std(b) == 0: continue
            w = fit_opt_nusvr(A0, (b - b.mean()) / np.std(b))
            Wb[i] = 100.0 * w / w.sum() if w.sum() > 0 else 0.0
        ratio_g = np.clip((B.mean(axis=1, keepdims=True) + 1.0) / (M_star.mean(axis=1, keepdims=True) + 1.0), 0.1, 10.0)
        As = A * ratio_g
        As0 = (As - As.mean(axis=0, keepdims=True)) / np.maximum(As.std(axis=0, keepdims=True), 1e-6)
        Ws = np.zeros((N, len(FINE43)))
        for i in range(N):
            b = B[:, i]
            if np.std(b) == 0: continue
            w = fit_opt_nusvr(As0, (b - b.mean()) / np.std(b))
            Ws[i] = 100.0 * w / w.sum() if w.sum() > 0 else 0.0
        return pd.DataFrame(100.0 * W0, index=M_cpm.columns, columns=FINE43), pd.DataFrame(Ws, index=M_cpm.columns, columns=FINE43), pd.DataFrame(Wb, index=M_cpm.columns, columns=FINE43)

    def run_music_dwls_nnls(M_cpm, sig_m, var_ct, mks, cal_ratios=None):
        A = sig_m.loc[FINE43, mks].T.values.astype(np.float64)
        V = var_ct.loc[FINE43, mks].T.values.astype(np.float64)
        A_cal = A * cal_ratios[mks].values[:, None] if cal_ratios is not None else A
        o_mu, o_dw, o_nn, o_cn = {}, {}, {}, {}
        for s in M_cpm.columns:
            b = M_cpm.loc[mks, s].values.astype(np.float64)
            w0, _ = nnls(A, b)
            o_nn[s] = (100.0 * w0 / w0.sum()) if w0.sum() > 0 else np.zeros(len(FINE43))
            wc, _ = nnls(A_cal, b)
            o_cn[s] = (100.0 * wc / wc.sum()) if wc.sum() > 0 else np.zeros(len(FINE43))
            if w0.sum() == 0:
                o_mu[s] = o_dw[s] = np.zeros(len(FINE43)); continue
            wm, wd = w0.copy(), w0.copy()
            for _ in range(5):
                vg = (V @ (wm ** 2)) + (b - A @ wm) ** 2 + 1e-2
                sqw = np.minimum(1.0 / np.sqrt(np.maximum(vg, 1e-4)), np.percentile(1.0 / np.sqrt(np.maximum(vg, 1e-4)), 95))
                wn, _ = nnls(A * sqw[:, None], b * sqw)
                if wn.sum() > 0: wm = wn
            for _ in range(4):
                wts = 1.0 / (np.maximum(A @ wd, 1e-3) ** 2)
                sqw = np.sqrt(np.minimum(wts, np.percentile(wts, 95)))
                wn, _ = nnls(A * sqw[:, None], b * sqw)
                if wn.sum() > 0: wd = wn
            o_mu[s] = 100.0 * wm / wm.sum()
            o_dw[s] = 100.0 * wd / wd.sum()
        return pd.DataFrame(o_mu, index=FINE43).T, pd.DataFrame(o_dw, index=FINE43).T, pd.DataFrame(o_nn, index=FINE43).T, pd.DataFrame(o_cn, index=FINE43).T

    # Same-visit paired Real/Bulk ratio for Calibrated NNLS (strictly fit on training donors only!)
    rm_match = rm_all.copy()
    rm_match["bulk_id"] = [
        f"{r['patient_id']}_negative" if r["condition"] == "negative" else
        (f"{r['patient_id']}_positive_{int(r['days_from_first_positive_swab'])}" if (r["condition"] == "positive" and pd.notna(r["days_from_first_positive_swab"])) else None)
        for _, r in rm_match.iterrows()
    ]
    mt_pairs = rm_match[rm_match["bulk_id"].isin(B_cpm_sym.columns)].copy()
    idx_mt = gv["sample_id"].isin(mt_pairs.index).values
    rpb_mt = pb_df_c.iloc[idx_mt].groupby(gv["sample_id"].values[idx_mt]).sum()
    rcpm_mt = rpb_mt.div(rpb_mt.sum(axis=1), axis=0) * 1e6

    def fit_cal_ratios(train_pts):
        sub_mt = mt_pairs[mt_pairs["patient_id"].isin(train_pts)]
        s_real = [s for s in sub_mt.index if s in rcpm_mt.index]
        b_ids = sub_mt.loc[s_real, "bulk_id"].tolist()
        rat = (B_cpm_sym.loc[common_g, b_ids].mean(axis=1) + 1.0) / (rcpm_mt.loc[s_real, common_g].mean(axis=0) + 1.0)
        return np.clip(rat, 0.1, 10.0)

    BULK_MARKERS_8 = {
        "NK CD56 bright": ["NCAM1", "XCL1", "XCL2", "KLRC1", "GZMK", "SELL"],
        "B_ASC_dividing": ["MKI67", "TOP2A", "JCHAIN", "MZB1", "PRDM1"],
        "T CD8 CM": ["CD8A", "CCR7", "SELL", "TCF7", "IL7R", "GPR183"],
        "CD16mono": ["FCGR3A", "LST1", "MS4A7", "TCF7L2", "AIF1"],
        "B_ASC_IgG": ["IGHG1", "IGHG2", "IGHG3", "IGHG4", "MZB1", "JCHAIN"],
        "B_ASC_IgA": ["IGHA1", "IGHA2", "JCHAIN", "MZB1"],
        "T CD8 CTL": ["CD8A", "CD8B", "GZMB", "PRF1", "GNLY", "FGFBP2"],
        "HSPC": ["CD34", "PROM1", "AVP", "SPINK2", "CRHBP", "HLF"],
    }
    def run_bulk_markers(M_cpm):
        lg = np.log2(M_cpm + 1.0)
        out = pd.DataFrame(0.0, index=M_cpm.columns, columns=FINE43)
        for fk, gl in BULK_MARKERS_8.items():
            gu = [g for g in gl if g in lg.index]
            out[fk] = lg.loc[gu].mean(axis=0)
        return out

    sh_samps = pos_p_all[pos_p_all["patient_id"].isin(SHARED27)].set_index("sample_id")
    nw_samps = pos_n.set_index("sample_id")
    all_real_noncll = sorted(gv["patient_id"].unique().tolist())
    fold_a, fold_b = set(SHARED27[::2]), set(SHARED27[1::2])
    non_sh = set(all_real_noncll) - set(SHARED27)

    ALL_METHODS = ["PBISC", "BayesPrism-upd", "BayesPrism-init", "CIBERSORTx-Smode", "CIBERSORTx-Bmode", "MuSiC", "Calib-NNLS", "DWLS", "Bulk-8-markers", "NNLS"]
    sh_buf = {m: [] for m in ALL_METHODS if m != "PBISC"}
    bp_z_sh_buf = {ct: [] for ct in ["HSPC", "T CD8 CM", "T CD8 CTL"]}

    for te_f, tr_f in [(fold_a, fold_b | non_sh), (fold_b, fold_a | non_sh)]:
        sids = sh_samps[sh_samps["patient_id"].isin(te_f)].index.tolist()
        sig_m, phi_p, var_ct, csz, mks = build_ref_stats(tr_f)
        cal_r = fit_cal_ratios(tr_f)
        bp1, bp2, bpz = run_bayesprism(B_cnt_sym[sids], phi_p, csz)
        _, css, csb = run_cibersortx(B_cpm_sym[sids], sig_m, mks)
        mu_d, dw_d, nn_d, cn_d = run_music_dwls_nnls(B_cpm_sym[sids], sig_m, var_ct, mks, cal_ratios=cal_r)
        bm_d = run_bulk_markers(B_cpm_sym[sids])
        for k, df_k in [("BayesPrism-upd", bp2), ("BayesPrism-init", bp1), ("CIBERSORTx-Smode", css), ("CIBERSORTx-Bmode", csb),
                        ("MuSiC", mu_d), ("Calib-NNLS", cn_d), ("DWLS", dw_d), ("Bulk-8-markers", bm_d), ("NNLS", nn_d)]:
            sh_buf[k].append(df_k)
        for ct in bp_z_sh_buf:
            bp_z_sh_buf[ct].append(bpz[ct])

    sh_pt = {"PBISC": PB_sh27_pt.loc[SHARED27]}
    for k, lst in sh_buf.items():
        sh_pt[k] = pd.concat(lst).join(sh_samps["patient_id"]).groupby("patient_id")[FINE43].mean().loc[SHARED27]

    sids_nw = nw_samps.index.tolist()
    sig_f, phi_f, var_f, csz_f, mks_f = build_ref_stats(all_real_noncll)
    cal_f = fit_cal_ratios(all_real_noncll)
    bp1_nw, bp2_nw, bp_z_nw = run_bayesprism(B_cnt_sym[sids_nw], phi_f, csz_f)
    _, css_nw, csb_nw = run_cibersortx(B_cpm_sym[sids_nw], sig_f, mks_f)
    mu_nw, dw_nw, nn_nw, cn_nw = run_music_dwls_nnls(B_cpm_sym[sids_nw], sig_f, var_f, mks_f, cal_ratios=cal_f)
    bm_nw = run_bulk_markers(B_cpm_sym[sids_nw])

    nw_pt = {"PBISC": N_pt.loc[NEW29]}
    for k, df_k in [("BayesPrism-upd", bp2_nw), ("BayesPrism-init", bp1_nw), ("CIBERSORTx-Smode", css_nw), ("CIBERSORTx-Bmode", csb_nw),
                    ("MuSiC", mu_nw), ("Calib-NNLS", cn_nw), ("DWLS", dw_nw), ("Bulk-8-markers", bm_nw), ("NNLS", nn_nw)]:
        nw_pt[k] = df_k.join(nw_samps["patient_id"]).groupby("patient_id")[FINE43].mean().loc[NEW29]

    # Build comprehensive benchmark table (Fig. 8a, 8b, S5a, S5b)
    BENCH_10 = [
        ("NK CD56 bright", "sev", -1), ("B_ASC_dividing", "sev", +1), ("T CD8 CM", "sev", -1), ("CD16mono", "sev", -1),
        ("B_ASC_IgG", "sev", +1), ("B_ASC_IgA", "sev", +1), ("HSPC", "sev", +1), ("T CD8 CTL", "sev", -1),
        ("T CD8 CTL", "death", -1), ("HSPC", "death", +1),
    ]
    bench_rows = []
    for m in ALL_METHODS:
        # Compute full 43-type (or pure BH8 for Bulk-8-markers!) Combined q-values under TRUE_PUB_11
        feat_fam = list(BULK_MARKERS_8.keys()) if m == "Bulk-8-markers" else FINE43
        ps_s_m, ps_d_m, eff_s_m, eff_d_m = [], [], {}, {}
        for f in feat_fam:
            yc = np.r_[rint_arr(R_pt[f].values), rint_arr(nw_pt[m][f].values)]
            rc, psc = stats.spearmanr(sev_c, yc)
            ac, pdc = auc_mw(yc, fat_c)
            ps_s_m.append(psc); ps_d_m.append(pdc)
            eff_s_m[f] = float(rc); eff_d_m[f] = float(ac)
        qs_s_m = dict(zip(feat_fam, multipletests(ps_s_m, method="fdr_bh")[1]))
        qs_d_m = dict(zip(feat_fam, multipletests(ps_d_m, method="fdr_bh")[1]))

        for f, ep, sgn in BENCH_10:
            sh_r = float(stats.spearmanr(R_pt.loc[SHARED27, f], sh_pt[m].loc[SHARED27, f])[0]) if sh_pt[m].loc[SHARED27, f].nunique() > 1 else np.nan
            if ep == "sev":
                ne, np_ = (stats.spearmanr(sev_n, nw_pt[m][f]) if nw_pt[m][f].nunique() > 1 else (0.0, 1.0))
                ce, cq = eff_s_m[f], qs_s_m[f]
                ok = bool((cq < 0.05) and (np.sign(ce) == sgn) and (np.sign(ne) == sgn))
            else:
                ne, np_ = auc_mw(nw_pt[m][f].values, fat_n)
                ce, cq = eff_d_m[f], qs_d_m[f]
                ok = bool((cq < 0.10) and (np.sign(ce - 0.5) == sgn) and (np.sign(ne - 0.5) == sgn))
            bench_rows.append(dict(
                method=m, feature=f, endpoint=ep, expected_sign=sgn,
                shared27_rho=sh_r, new29_eff=float(ne), new29_p=float(np_),
                comb66_eff=float(ce), comb66_q=float(cq), passed_concordant_FDR=ok,
                fdr_family=("BH8" if m == "Bulk-8-markers" else "BH43"),
            ))
    df_bench = pd.DataFrame(bench_rows)
    df_bench.to_csv(data_out / "fig8_deconv_10methods_true11.csv", index=False)

    # Wave-1 New-29 HSPC death AUC with 10,000-resample stratified CI + paired Delta-AUC vs PBISC!
    nw29_ci_rows = []
    pb_h_v = nw_pt["PBISC"]["HSPC"].values
    i1_n, i0_n = np.where(fat_n == 1)[0], np.where(fat_n == 0)[0]
    w1_n = rng_main.multinomial(len(i1_n), np.ones(len(i1_n)) / len(i1_n), size=20000).astype(np.float64)
    w0_n = rng_main.multinomial(len(i0_n), np.ones(len(i0_n)) / len(i0_n), size=20000).astype(np.float64)
    K_pb = ((pb_h_v[i1_n, None] > pb_h_v[None, i0_n]).astype(np.float64) + 0.5 * (pb_h_v[i1_n, None] == pb_h_v[None, i0_n])) / (len(i1_n) * len(i0_n))
    bs_pb = ((w1_n @ K_pb) * w0_n).sum(axis=1)

    for m in ALL_METHODS:
        hv = nw_pt[m]["HSPC"].values
        a_h, lo_h, hi_h, p_h = fast_strat_boot_auc(hv, fat_n, rng_main, n_boot=10000)
        rs_h, ps_h = (stats.spearmanr(sev_n, hv) if np.std(hv) > 0 else (0.0, 1.0))
        K_m = ((hv[i1_n, None] > hv[None, i0_n]).astype(np.float64) + 0.5 * (hv[i1_n, None] == hv[None, i0_n])) / (len(i1_n) * len(i0_n))
        bs_m = ((w1_n @ K_m) * w0_n).sum(axis=1)
        dbs = bs_pb - bs_m
        nw29_ci_rows.append(dict(
            method=m, hspc_sev_rho=float(rs_h), hspc_sev_p=float(ps_h),
            hspc_dth_auc=a_h, hspc_dth_lo=lo_h, hspc_dth_hi=hi_h, hspc_dth_p=p_h,
            delta_pbisc_minus_m=float(auc_mw(pb_h_v, fat_n)[0] - a_h),
            delta_lo=float(np.percentile(dbs, 2.5)), delta_hi=float(np.percentile(dbs, 97.5)),
            hspc_median_pct=float(np.median(hv)),
        ))
    df_nw29_ci = pd.DataFrame(nw29_ci_rows)
    df_nw29_ci.to_csv(data_out / "fig8_new29_hspc_10methods_CI.csv", index=False)

    # Unbiased 43-type & Rare-23 Summary across the 9 full-43 methods (Shared-27 fidelity + T1 concordance vs Real-37)
    vec_r43, vec_r23 = [], []
    for f in FINE43:
        rr = stats.spearmanr(sev_r, R_pt[f])[0]
        ar = auc_mw(R_pt[f], fat_r)[0]
        vec_r43.extend([rr, 2.0 * (ar - 0.5)])
        if f in RARE23: vec_r23.extend([rr, 2.0 * (ar - 0.5)])
    unbiased_rows = []
    for m in [x for x in ALL_METHODS if x != "Bulk-8-markers"]:
        rhos43 = [stats.spearmanr(R_pt.loc[SHARED27, f], sh_pt[m].loc[SHARED27, f])[0] if sh_pt[m].loc[SHARED27, f].nunique() > 1 else np.nan for f in FINE43]
        s_r43 = pd.Series(rhos43, index=FINE43)
        vn43, vn23 = [], []
        for f in FINE43:
            rn = stats.spearmanr(sev_n, nw_pt[m][f])[0] if nw_pt[m][f].nunique() > 1 else 0.0
            an = auc_mw(nw_pt[m][f], fat_n)[0]
            vn43.extend([rn, 2.0 * (an - 0.5)])
            if f in RARE23: vn23.extend([rn, 2.0 * (an - 0.5)])
        n_pass10 = int(df_bench[df_bench["method"] == m]["passed_concordant_FDR"].sum())
        unbiased_rows.append(dict(
            method=m, med_rho_all43=float(s_r43.median()), mean_rho_all43=float(s_r43.mean()),
            med_rho_rare23=float(s_r43.loc[RARE23].median()), mean_rho_rare23=float(s_r43.loc[RARE23].mean()),
            rare23_nonzero=int(s_r43.loc[RARE23].notna().sum()),
            T1_conc_all43=float(stats.spearmanr(vec_r43, vn43)[0]),
            T1_conc_rare23=float(stats.spearmanr(vec_r23, vn23)[0]),
            pass_10_signals=n_pass10,
        ))
    df_unb = pd.DataFrame(unbiased_rows)
    df_unb.to_csv(data_out / "suppS5_unbiased_deconv_9methods_summary.csv", index=False)
    print("[Point 4] Unbiased 9-method summary:\n", df_unb.round(3).to_string(index=False))

    # Single-gene raw bulk check in Wave-1 New-29 (for S5b reference)
    def bulk_gene_new29(ens_id):
        b_cpm_ens = B_raw.div(B_raw.sum(axis=0), axis=1) * 1e6
        return np.log2(b_cpm_ens.loc[ens_id, nw_samps.index] + 1.0).groupby(nw_samps["patient_id"].values).mean().loc[NEW29]
    sg_rows = []
    for glab, gvec, is_pb in [
        ("PBISC HSPC (%)", N_pt["HSPC"], True),
        ("Bulk PROM1 (CD133)", bulk_gene_new29("ENSG00000007062"), False),
        ("Bulk AVP", bulk_gene_new29("ENSG00000101200"), False),
        ("Bulk CD34", bulk_gene_new29("ENSG00000174059"), False),
        ("Bulk SPINK2", bulk_gene_new29("ENSG00000128040"), False),
    ]:
        av, lo, hi, pv = fast_strat_boot_auc(gvec.values, fat_n, rng_main, n_boot=10000)
        sg_rows.append(dict(label=glab, is_pbisc=is_pb, auc=av, lo=lo, hi=hi, p=pv))
    pd.DataFrame(sg_rows).to_csv(data_out / "suppS5_single_gene_hspc_new29.csv", index=False)

    # ------------------------------------------------------------------
    # 5. Point 5: HSC / HSPC Transcriptome (Disjoint n=62: Real 33 [5d] + PBISC new 29 [6d])
    #    16 Unique Biological Modules + 20 Candidate Genes + Full 5,762-Gene Scan + BayesPrism Z comparison!
    # ------------------------------------------------------------------
    print("\n[Point 5] Computing HSPC transcriptome (n=62 disjoint) under Severity & TRUE_PUB_11 Death...")
    r_h5 = ad.read_h5ad(reg(ESKD_DIR / "data" / "real_rare_subset_cells_counts.h5ad", "real_rare_cells_h5ad"))
    p_h5 = ad.read_h5ad(reg(ESKD_DIR / "data" / "pbisc_rare_subset_cells_counts.h5ad", "pbisc_rare_cells_h5ad"))
    p_sym = pd.Series(p_h5.var["gene_symbol"].values, index=np.arange(p_h5.n_vars)); p_sym = p_sym[~p_sym.duplicated()]
    r_sym = pd.Series(r_h5.var_names.values, index=np.arange(r_h5.n_vars)); r_sym = r_sym[~r_sym.duplicated()]
    shared_g5762 = sorted(set(p_sym.values) & set(r_sym.values))
    r_gmap = {s: i for i, s in r_sym.items()}
    p_gmap = {s: i for i, s in p_sym.items()}

    def hspc_pt_expr(adata, idx_c, src_nm, pid_col, min_cells=3):
        m = (adata.obs["fine"] == "HSPC").values
        obs = adata.obs.loc[m]
        vs = set(pos_samps[pos_samps["source"] == src_nm]["sample_id"])
        keep = obs["sample_id"].astype(str).isin(vs).values
        obs = obs.loc[keep]
        X = adata.X[m][keep][:, idx_c]
        if sp.issparse(X): X = X.toarray()
        X = np.log1p(np.divide(X, np.maximum(X.sum(axis=1, keepdims=True), 1.0)) * 1e4)
        pids = obs[pid_col].astype(str).values
        out, nc = {}, {}
        for pid in sorted(set(pids)):
            pm = pids == pid
            if pm.sum() >= min_cells:
                out[pid] = X[pm].mean(axis=0); nc[pid] = int(pm.sum())
        return pd.DataFrame.from_dict(out, orient="index", columns=shared_g5762), nc

    hr_pt, hr_nc = hspc_pt_expr(r_h5, [r_gmap[g] for g in shared_g5762], "Real", "patient_id", min_cells=3)
    hp_pt_all, hp_nc = hspc_pt_expr(p_h5, [p_gmap[g] for g in shared_g5762], "PBISC", "donor_id", min_cells=3)
    hp_pt = hp_pt_all.loc[[p for p in NEW29 if p in hp_pt_all.index]]
    sh24_hspc = sorted(set(hr_pt.index) & set(hp_pt_all.index) & set(SHARED27))

    # Exact 16 unique biological HSPC modules (matching Folder 25 unique-16 set, deduplicating Cycling/Cycle4!)
    MOD16_DEF = {
        "Ery5": ["GATA1", "KLF1", "TFRC", "APOC1", "BLVRB"],
        "Cycle4": ["MKI67", "TOP2A", "STMN1", "TYMS"],
        "MEP5": ["GATA1", "GATA2", "KLF1", "TFRC", "ITGA2B"],
        "MEP_EryMk": ["GATA1", "GATA2", "KLF1", "TFRC", "ITGA2B", "PF4", "PPBP", "GP9"],
        "MkP_Platelet": ["ITGA2B", "PF4", "PPBP", "GP9", "MPL", "ITGB3", "TUBB1", "PLEK"],
        "HSC5": ["CD34", "PROM1", "AVP", "HLF", "CRHBP"],
        "Stem_HSC": ["CD34", "PROM1", "AVP", "HLF", "CRHBP", "SPINK2", "HOPX", "MEIS1"],
        "Myeloid4": ["MPO", "ELANE", "CSF3R", "CEBPD"],
        "GMP_Myeloid": ["MPO", "ELANE", "CTSG", "LYZ", "S100A8", "S100A9", "S100A12", "CSF3R", "CEBPD", "MS4A3", "RETN", "MNDA"],
        "Granule_GMP": ["MPO", "ELANE", "CTSG", "AZU1", "PRTN3", "MS4A3"],
        "Alarmin_Inflam": ["S100A8", "S100A9", "S100A12", "LYZ", "RETN", "MNDA"],
        "ISG_HSPC": ["ISG15", "IFI6", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "STAT1", "IRF7", "XAF1"],
        "MEP_minus_Stem": None,
        "GMP_minus_Stem": None,
        "Alarmin_minus_Stem": None,
        "HSC_broad_q": ["CD34", "PROM1", "AVP", "HLF", "CRHBP", "SPINK2", "HOPX", "MEIS1", "MLLT3"],
    }
    def calc_hspc_mod(df_e, mname):
        if mname == "MEP_minus_Stem":
            return calc_hspc_mod(df_e, "MEP_EryMk") - calc_hspc_mod(df_e, "Stem_HSC")
        if mname == "GMP_minus_Stem":
            return calc_hspc_mod(df_e, "GMP_Myeloid") - calc_hspc_mod(df_e, "Stem_HSC")
        if mname == "Alarmin_minus_Stem":
            return calc_hspc_mod(df_e, "Alarmin_Inflam") - calc_hspc_mod(df_e, "Stem_HSC")
        gu = [g for g in MOD16_DEF[mname] if g in df_e.columns]
        return df_e[gu].mean(axis=1)

    hmod_rows = []
    for mn_k in MOD16_DEF:
        vr = calc_hspc_mod(hr_pt, mn_k)
        vn = calc_hspc_mod(hp_pt, mn_k)
        v_sh_p = calc_hspc_mod(hp_pt_all.loc[sh24_hspc], mn_k)
        vc = pd.concat([rint(vr), rint(vn)])
        sh_r, sh_p = stats.spearmanr(vr.loc[sh24_hspc], v_sh_p)
        for coh_nm, vv in [("Real", vr), ("PBISC new", vn), ("Combined", vc)]:
            sv = pt_cl.loc[vv.index, "sev"].values
            fv = pt_cl.loc[vv.index, "fat11"].values
            r_s, r_lo, r_hi, p_s = fast_boot_rho(sv, vv.values, rng_main, n_boot=5000)
            a_d, a_lo, a_hi, p_d = fast_strat_boot_auc(vv.values, fv, rng_main, n_boot=10000)
            hmod_rows.append(dict(
                module=mn_k, cohort=coh_nm, n=len(vv), deaths=int(fv.sum()),
                rho_sev=r_s, rho_lo=r_lo, rho_hi=r_hi, p_sev=p_s,
                AUC_dth11=a_d, AUC_lo=a_lo, AUC_hi=a_hi, p_dth11=p_d,
                shared24_rho=float(sh_r), shared24_p=float(sh_p),
            ))
    df_hmod = pd.DataFrame(hmod_rows)
    for coh_nm in ["Real", "PBISC new", "Combined"]:
        m_c = df_hmod["cohort"] == coh_nm
        df_hmod.loc[m_c, "q16_sev"] = multipletests(df_hmod.loc[m_c, "p_sev"], method="fdr_bh")[1]
        df_hmod.loc[m_c, "q16_dth11"] = multipletests(df_hmod.loc[m_c, "p_dth11"], method="fdr_bh")[1]
    df_hmod.to_csv(data_out / "fig8_hspc_modules_16unique_true11.csv", index=False)

    # Save 62-patient HSPC Ery5 + KIAA0391 patient-level table (for Fig. 8d!)
    ery5_r, ery5_n = calc_hspc_mod(hr_pt, "Ery5"), calc_hspc_mod(hp_pt, "Ery5")
    df_hspc62 = pd.DataFrame({
        "patient_id": list(hr_pt.index) + list(hp_pt.index),
        "arm": ["Real"] * len(hr_pt) + ["PBISC new"] * len(hp_pt),
        "n_hspc_cells": [hr_nc[p] for p in hr_pt.index] + [hp_nc[p] for p in hp_pt.index],
        "Ery5_raw": np.r_[ery5_r.values, ery5_n.values],
        "Ery5_ns": np.r_[rint(ery5_r).values, rint(ery5_n).values],
        "KIAA0391_ns": np.r_[rint(hr_pt["KIAA0391"]).values, rint(hp_pt["KIAA0391"]).values],
    })
    df_hspc62["sev"] = df_hspc62["patient_id"].map(pt_cl["sev"])
    df_hspc62["fat11"] = df_hspc62["patient_id"].map(pt_cl["fat11"])
    df_hspc62.to_csv(data_out / "fig8_hspc_patients_62_true11.csv", index=False)

    # 20 candidate genes in HSPC (Severity & TRUE_PUB_11 Death)
    GENES_20 = ["GATA1", "KLF1", "TFRC", "BLVRB", "APOC1", "GATA2", "ITGA2B", "MKI67", "TOP2A", "STMN1", "TYMS", "CD34", "AVP", "HLF", "HOPX", "CRHBP", "MPO", "ELANE", "CSF3R", "CEBPD"]
    g20_rows = []
    for g in GENES_20:
        vr, vn = hr_pt[g], hp_pt[g]
        vc = pd.concat([rint(vr), rint(vn)])
        for coh_nm, vv in [("Real", vr), ("PBISC new", vn), ("Combined", vc)]:
            rs, ps = stats.spearmanr(pt_cl.loc[vv.index, "sev"], vv)
            ad_, pd_ = auc_mw(vv.values, pt_cl.loc[vv.index, "fat11"].values)
            g20_rows.append(dict(gene=g, cohort=coh_nm, rho_sev=float(rs), p_sev=float(ps), AUC_dth11=float(ad_), p_dth11=float(pd_)))
    df_g20 = pd.DataFrame(g20_rows)
    for coh_nm in ["Real", "PBISC new", "Combined"]:
        mc = df_g20["cohort"] == coh_nm
        df_g20.loc[mc, "q20_sev"] = multipletests(df_g20.loc[mc, "p_sev"], method="fdr_bh")[1]
        df_g20.loc[mc, "q20_dth11"] = multipletests(df_g20.loc[mc, "p_dth11"], method="fdr_bh")[1]
    df_g20.to_csv(data_out / "suppS5_hspc_20gene_true11.csv", index=False)

    # ------------------------------------------------------------------
    # 6. Point 6: T-cell Type-I IFN / ISG Across All 13 T-cell Compartments
    #    + Day-Adjusted GEE + Source-Specific Time x Fatal Collapse + BayesPrism Z Comparison!
    # ------------------------------------------------------------------
    print("\n[Point 6] Computing T-cell ISG across all 13 T-cell compartments + day-adjusted GEE + BayesPrism Z comparison...")
    ISG_USE = [g for g in ["ISG15", "IFI6", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "XAF1", "STAT1", "IRF7", "EPSTI1", "HERC5", "SIGLEC1"] if g in r_genes]
    r_isg_idx = [r_genes.get_loc(g) for g in ISG_USE]
    pb_sym_all = pd.Index(pb_h5.var["gene_symbol"].astype(str))
    pb_isg_idx = [int(np.where(pb_sym_all == g)[0][0]) for g in ISG_USE]
    X_pb = pb_h5.X
    tot_umi_pb = np.asarray(X_pb.sum(axis=1)).ravel().astype(np.float64)
    X_pb_isg = X_pb[:, pb_isg_idx]
    if sp.issparse(X_pb_isg): X_pb_isg = X_pb_isg.toarray()
    pb_tot_s = pd.DataFrame({
        "sample_id": pb_h5.obs["sample_id"].astype(str).values,
        "fine": pb_h5.obs["celltype_fine"].astype(str).values,
        "tot_umi": tot_umi_pb, "ncells": 1,
    })
    for j, g in enumerate(ISG_USE): pb_tot_s[g] = X_pb_isg[:, j]

    all_r_fines = sorted(grps["fine"].unique().tolist())
    t_fines = [f for f in all_r_fines if f.startswith("T ") or f in ("MAIT", "NKT")]
    GROUPS_T = [
        ("All T cells", set(t_fines) - {"NKT"}),
        ("All CD4 T", {f for f in t_fines if "CD4" in f or "T reg" in f}),
        ("All CD8 T", {f for f in t_fines if "CD8" in f}),
        ("T CD4 Naive", {"T CD4 Naive"}),
        ("T CD4 Helper", {"T CD4 Helper"}),
        ("T CD4 CM/EM/CTL", {f for f in t_fines if f in ("T CD4 CM", "T CD4 EM", "T CD4 CTL", "T CD4 Helper")}),
        ("T CD4 CTL", {"T CD4 CTL"}),
        ("T reg", {"T reg"}),
        ("T CD8 Naive", {"T CD8 Naive"}),
        ("T CD8 CM", {"T CD8 CM"}),
        ("T CD8 CTL", {"T CD8 CTL"}),
        ("T g/d", {"T g/d"}),
        ("MAIT", {"MAIT"}),
    ]

    r_pos_m = pos_r.set_index("sample_id")
    p_pos_m = pos_p_all.set_index("sample_id")
    neg_samps = S_meta[(S_meta["condition"] == "negative") & S_meta["in_primary_scope"].astype(bool) & S_meta["qc_pass_min200cells"].astype(bool) & (~S_meta["patient_id"].isin(CLL))]
    r_neg_m = neg_samps[neg_samps["source"] == "Real"].set_index("sample_id")
    p_neg_m = neg_samps[neg_samps["source"] == "PBISC"].set_index("sample_id")

    tisg_rows, tisg_gee_rows = [], []
    cd8_ctl_isg_217 = None
    for gname, fset in GROUPS_T:
        m_r = grps["fine"].isin(fset).values
        df_rs = pd.DataFrame(r_mat[m_r][:, r_isg_idx], columns=ISG_USE)
        df_rs["sample_id"] = grps.loc[m_r, "sample_id"].values
        df_rs["tot_umi"] = r_mat[m_r].sum(axis=1); df_rs["ncells"] = grps.loc[m_r, "ncells"].values
        df_rs = df_rs.groupby("sample_id").sum(); df_rs = df_rs[df_rs["ncells"] >= 5].copy()
        df_rs["isg"] = np.log1p(df_rs[ISG_USE].div(np.maximum(df_rs["tot_umi"], 1.0), axis=0) * 1e4).mean(axis=1)

        df_ps = pb_tot_s[pb_tot_s["fine"].isin(fset)].groupby("sample_id")[["tot_umi", "ncells"] + ISG_USE].sum()
        df_ps = df_ps[df_ps["ncells"] >= 5].copy()
        df_ps["isg"] = np.log1p(df_ps[ISG_USE].div(np.maximum(df_ps["tot_umi"], 1.0), axis=0) * 1e4).mean(axis=1)

        r_ps = df_rs.reindex(r_pos_m.index).dropna().join(r_pos_m[["patient_id", "days_from_first_positive_swab"]])
        p_ps = df_ps.reindex(p_pos_m.index).dropna().join(p_pos_m[["patient_id", "days_from_first_positive_swab"]])
        r_pt_i = r_ps.groupby("patient_id")["isg"].mean()
        p_pt_i_all = p_ps.groupby("patient_id")["isg"].mean()
        n_pt_i = p_pt_i_all.reindex([p for p in NEW29 if p in p_pt_i_all.index])
        c_ns_i = pd.concat([rint(r_pt_i), rint(n_pt_i)])

        sh_c = sorted(set(r_pt_i.index) & set(p_pt_i_all.index) & set(SHARED27))
        shr, shlo, shhi, shp = fast_boot_rho(r_pt_i.loc[sh_c].values, p_pt_i_all.loc[sh_c].values, rng_main, n_boot=3000)

        rs_r, prs_r = stats.spearmanr(pt_cl.loc[r_pt_i.index, "sev"], r_pt_i)
        rs_n, prs_n = stats.spearmanr(pt_cl.loc[n_pt_i.index, "sev"], n_pt_i)
        rs_c, rslo_c, rshi_c, prs_c = fast_boot_rho(pt_cl.loc[c_ns_i.index, "sev"].values, c_ns_i.values, rng_main, n_boot=5000)
        ad_r, pad_r = auc_mw(r_pt_i.values, pt_cl.loc[r_pt_i.index, "fat11"].values)
        ad_n, pad_n = auc_mw(n_pt_i.values, pt_cl.loc[n_pt_i.index, "fat11"].values)
        ad_c, adlo_c, adhi_c, pad_c = fast_strat_boot_auc(c_ns_i.values, pt_cl.loc[c_ns_i.index, "fat11"].values, rng_main, n_boot=10000)

        # Early <=7d vs Late >7d
        r_e = r_ps[pd.to_numeric(r_ps["days_from_first_positive_swab"], errors="coerce") <= 7].groupby("patient_id")["isg"].mean()
        n_e = p_ps[p_ps["patient_id"].isin(NEW29) & (pd.to_numeric(p_ps["days_from_first_positive_swab"], errors="coerce") <= 7)].groupby("patient_id")["isg"].mean()
        c_e = pd.concat([rint(r_e), rint(n_e)])
        es, eps = stats.spearmanr(pt_cl.loc[c_e.index, "sev"], c_e)
        r_l = r_ps[pd.to_numeric(r_ps["days_from_first_positive_swab"], errors="coerce") > 7].groupby("patient_id")["isg"].mean()
        n_l = p_ps[p_ps["patient_id"].isin(NEW29) & (pd.to_numeric(p_ps["days_from_first_positive_swab"], errors="coerce") > 7)].groupby("patient_id")["isg"].mean()
        c_l = pd.concat([rint(r_l), rint(n_l)])
        ls, lps = stats.spearmanr(pt_cl.loc[c_l.index, "sev"], c_l)

        # Day-residualized patient check
        pt_d_m = pos_217.groupby("patient_id")["days"].mean().loc[c_ns_i.index].values
        slp_d, ic_d = np.polyfit(stats.rankdata(pt_d_m), c_ns_i.values, 1)
        res_i = c_ns_i.values - (slp_d * stats.rankdata(pt_d_m) + ic_d)
        r_dadj, p_dadj = stats.spearmanr(pt_cl.loc[c_ns_i.index, "sev"].values, res_i)
        a_dadj, pa_dadj = auc_mw(res_i, pt_cl.loc[c_ns_i.index, "fat11"].values)

        tisg_rows.append(dict(
            subset=gname, n_comb=len(c_ns_i), n_sh=len(sh_c),
            sh_rho=shr, sh_lo=shlo, sh_hi=shhi, sh_p=shp,
            r_sev=float(rs_r), n_sev=float(rs_n), c_sev=rs_c, c_sev_lo=rslo_c, c_sev_hi=rshi_c, cp_sev=prs_c,
            r_auc11=float(ad_r), n_auc11=float(ad_n), c_auc11=ad_c, c_auc11_lo=adlo_c, c_auc11_hi=adhi_c, cp_auc11=pad_c,
            early_sev=float(es), early_p_sev=float(eps), late_sev=float(ls), late_p_sev=float(lps),
            day_adj_sev=float(r_dadj), day_adj_p_sev=float(p_dadj), day_adj_auc11=float(a_dadj), day_adj_p_auc11=float(pa_dadj),
        ))

        # Source-specific GEE for key T compartments
        if gname in ("T CD8 CTL", "T CD4 CTL", "All CD8 T", "All CD4 T", "All T cells"):
            d_r = pos_217[pos_217["source"] == "Real"].merge(df_rs[["isg"]], left_on="sample_id", right_index=True)
            d_p = pos_217[pos_217["source"] == "PBISC"].merge(df_ps[["isg"]], left_on="sample_id", right_index=True)
            d_r["isg_z"] = rint(d_r["isg"]).values
            d_p["isg_z"] = rint(d_p["isg"]).values
            if gname == "T CD8 CTL":
                cd8_ctl_isg_217 = pd.concat([d_r, d_p], ignore_index=True)[["source", "sample_id", "patient_id", "days", "severity_ord", "fatal_pub11", "isg", "isg_z"]]
            for clab, sub in [("Real", d_r), ("PBISC new", d_p), ("Combined", pd.concat([d_r, d_p], ignore_index=True))]:
                sub = sub.dropna(subset=["days", "isg_z", "severity_ord", "fatal_pub11"]).copy()
                sub["wk"] = sub["days"] / 7.0; sub["fat_wk"] = sub["fatal_pub11"] * sub["wk"]
                ms = sm.GEE(sub["isg_z"], sm.add_constant(sub[["severity_ord", "wk"]]), groups=sub["patient_id"], cov_struct=sm.cov_struct.Exchangeable()).fit()
                mf = sm.GEE(sub["isg_z"], sm.add_constant(sub[["fatal_pub11", "wk"]]), groups=sub["patient_id"], cov_struct=sm.cov_struct.Exchangeable()).fit()
                mfi = sm.GEE(sub["isg_z"], sm.add_constant(sub[["fatal_pub11", "wk", "fat_wk"]]), groups=sub["patient_id"], cov_struct=sm.cov_struct.Exchangeable()).fit()
                ci_fi = mfi.conf_int().loc["fat_wk"]
                tisg_gee_rows.append(dict(
                    subset=gname, cohort=clab, draws=len(sub),
                    day_adj_sev_beta=float(ms.params["severity_ord"]), day_adj_sev_p=float(ms.pvalues["severity_ord"]),
                    day_adj_fat11_beta=float(mf.params["fatal_pub11"]), day_adj_fat11_p=float(mf.pvalues["fatal_pub11"]),
                    fat_x_wk_beta=float(mfi.params["fat_wk"]), fat_x_wk_lo=float(ci_fi[0]), fat_x_wk_hi=float(ci_fi[1]), fat_x_wk_p=float(mfi.pvalues["fat_wk"]),
                ))

    df_tisg = pd.DataFrame(tisg_rows)
    df_tisg["q13_sev"] = multipletests(df_tisg["cp_sev"], method="fdr_bh")[1]
    df_tisg["q13_auc11"] = multipletests(df_tisg["cp_auc11"], method="fdr_bh")[1]
    df_tisg.to_csv(data_out / "fig8_tcell_isg_13subsets_true11.csv", index=False)
    pd.DataFrame(tisg_gee_rows).to_csv(data_out / "suppS5_tcell_isg_gee_true11.csv", index=False)
    cd8_ctl_isg_217.to_csv(data_out / "fig8_cd8ctl_isg_longitudinal_draws.csv", index=False)

    # Head-to-head PBISC vs BayesPrism Posterior Z cell-type expression table (Shared-27 fidelity + New-29 clinical!)
    ERY5 = ["GATA1", "KLF1", "TFRC", "APOC1", "BLVRB"]
    bpz_comp_rows = []
    for ct, gl, mnm in [("HSPC", ERY5, "Ery5"), ("T CD8 CTL", ISG_USE, "ISG"), ("T CD8 CM", ISG_USE, "ISG")]:
        gl_bp = [g for g in gl if g in BP_G]
        z_sh_bp = pd.concat(bp_z_sh_buf[ct])[gl_bp].mean(axis=1).to_frame("z").join(sh_samps["patient_id"]).groupby("patient_id")["z"].mean()
        z_nw_bp = bp_z_nw[ct][gl_bp].mean(axis=1).to_frame("z").join(nw_samps["patient_id"]).groupby("patient_id")["z"].mean().loc[NEW29]
        # Real ground truth for Shared-27
        if ct == "HSPC":
            r_gt = calc_hspc_mod(hr_pt, "Ery5")
            pb_sh_v = calc_hspc_mod(hp_pt_all.loc[sh24_hspc], "Ery5")
            pb_nw_v = calc_hspc_mod(hp_pt, "Ery5")
        else:
            m_r = grps["fine"].eq(ct).values
            df_rs = pd.DataFrame(r_mat[m_r][:, r_isg_idx], columns=ISG_USE)
            df_rs["sample_id"] = grps.loc[m_r, "sample_id"].values; df_rs["tot_umi"] = r_mat[m_r].sum(axis=1); df_rs["ncells"] = grps.loc[m_r, "ncells"].values
            df_rs = df_rs.groupby("sample_id").sum(); df_rs = df_rs[df_rs["ncells"] >= 5]
            df_rs["sc"] = np.log1p(df_rs[ISG_USE].div(np.maximum(df_rs["tot_umi"], 1.0), axis=0) * 1e4).mean(axis=1)
            r_gt = df_rs.reindex(r_pos_m.index).dropna().join(r_pos_m["patient_id"]).groupby("patient_id")["sc"].mean()
            df_ps = pb_tot_s[pb_tot_s["fine"].eq(ct)].groupby("sample_id")[["tot_umi", "ncells"] + ISG_USE].sum()
            df_ps = df_ps[df_ps["ncells"] >= 5]
            df_ps["sc"] = np.log1p(df_ps[ISG_USE].div(np.maximum(df_ps["tot_umi"], 1.0), axis=0) * 1e4).mean(axis=1)
            p_gt = df_ps.reindex(p_pos_m.index).dropna().join(p_pos_m["patient_id"]).groupby("patient_id")["sc"].mean()
            pb_sh_v = p_gt.reindex([p for p in SHARED27 if p in p_gt.index])
            pb_nw_v = p_gt.reindex([p for p in NEW29 if p in p_gt.index])
        sh_idx_pb = sorted(set(r_gt.index) & set(pb_sh_v.index))
        sh_idx_bp = sorted(set(r_gt.index) & set(z_sh_bp.index))
        for m_lbl, sh_r_v, nw_v in [
            ("PBISC", stats.spearmanr(r_gt.loc[sh_idx_pb], pb_sh_v.loc[sh_idx_pb])[0], pb_nw_v),
            ("BayesPrism Z", stats.spearmanr(r_gt.loc[sh_idx_bp], z_sh_bp.loc[sh_idx_bp])[0], z_nw_bp),
        ]:
            rs_n, ps_n = stats.spearmanr(pt_cl.loc[nw_v.index, "sev"], nw_v)
            ad_n, pd_n = auc_mw(nw_v.values, pt_cl.loc[nw_v.index, "fat11"].values)
            bpz_comp_rows.append(dict(
                target=f"{ct} {mnm}", method=m_lbl, shared_rho=float(sh_r_v),
                new29_sev_rho=float(rs_n), new29_sev_p=float(ps_n),
                new29_dth11_auc=float(ad_n), new29_dth11_p=float(pd_n),
            ))
    df_bpz_comp = pd.DataFrame(bpz_comp_rows)
    df_bpz_comp.to_csv(data_out / "suppS5_pbisc_vs_bayesprism_Z_expression.csv", index=False)

    # ------------------------------------------------------------------
    # 7. Longitudinal Paired n=43, Source-Split GEE, Hospital Labs, & 500k Clinical CV (TRUE_PUB_11!)
    # ------------------------------------------------------------------
    # Paired n=43 Early vs Late HSPC windows
    e_h = pos_217[pos_217["days"] <= 7].groupby(["patient_id", "source"]).agg(v=("HSPC", "mean"), f=("fatal_pub11", "first")).reset_index()
    l_h = pos_217[pos_217["days"] > 7].groupby(["patient_id", "source"]).agg(v=("HSPC", "mean"), f=("fatal_pub11", "first")).reset_index()
    both43 = sorted(set(e_h["patient_id"]) & set(l_h["patient_id"]))
    eb = e_h[e_h["patient_id"].isin(both43)].set_index("patient_id").loc[both43]
    lb = l_h[l_h["patient_id"].isin(both43)].set_index("patient_id").loc[both43]
    eb["z"] = eb.groupby("source")["v"].transform(rint)
    lb["z"] = lb.groupby("source")["v"].transform(rint)
    p43_rows = []
    for wk_k, vec_w in [("early", eb["z"].values), ("late", lb["z"].values), ("delta", (lb["z"] - eb["z"]).values)]:
        a_, lo_, hi_, p_ = fast_strat_boot_auc(vec_w, eb["f"].values, rng_main, n_boot=10000)
        p43_rows.append(dict(window=wk_k, n=len(both43), deaths=int(eb["f"].sum()), AUC=a_, CI_lo=lo_, CI_hi=hi_, P=p_))
    pd.DataFrame(p43_rows).to_csv(data_out / "suppS5_paired43_window_AUC_true11.csv", index=False)

    # Source-specific GEE Time x Fatal for HSPC & T CD8 CTL (TRUE_PUB_11)
    gee_rows = []
    for fk in ["HSPC", "T CD8 CTL"]:
        for clab, sub in [("Real", pos_217[pos_217["source"] == "Real"]), ("New", pos_217[pos_217["source"] == "PBISC"]), ("Combined", pos_217)]:
            sub = sub.dropna(subset=["days", fk + "_ns", "fatal_pub11"]).copy()
            sub["wk"] = sub["days"] / 7.0; sub["intx"] = sub["wk"] * sub["fatal_pub11"]
            m = sm.GEE(sub[fk + "_ns"], sm.add_constant(sub[["wk", "fatal_pub11", "intx"]]), groups=sub["patient_id"], cov_struct=sm.cov_struct.Exchangeable()).fit()
            ci_i = m.conf_int().loc["intx"]
            gee_rows.append(dict(feature=fk, cohort=clab, draws=len(sub), beta=float(m.params["intx"]), CI_lo=float(ci_i[0]), CI_hi=float(ci_i[1]), P=float(m.pvalues["intx"])))
    pd.DataFrame(gee_rows).to_csv(data_out / "suppS5_source_GEE_robust_CI_true11.csv", index=False)

    # Patient-mean hospital labs vs immune subsets & rank-residualized HSPC death AUC (Real n=37, 5 deaths)
    LABS_6 = ["temp_CRP", "temp_lymph", "NLR", "temp_ferritin", "temp_ddimer", "temp_trop"]
    pos_r_lab = pos_r.copy()
    pos_r_lab["NLR"] = pos_r_lab["temp_neut"] / pos_r_lab["temp_lymph"].replace(0, np.nan)
    lab_rows = []
    for fk in ["NK CD56 bright", "T CD8 CM", "T CD8 CTL", "B_ASC_IgG", "B_ASC_IgA", "HSPC"]:
        for lk in LABS_6:
            pm = pos_r_lab.groupby("patient_id")[lk].mean().dropna()
            pids = sorted(set(pm.index) & set(REAL37))
            x = pm.loc[pids].values
            y = PT_NS.loc[pids, fk].values
            f11 = PT_NS.loc[pids, "fat11"].values
            r_l, p_l = stats.spearmanr(x, y)
            slp, ic = np.polyfit(stats.rankdata(x), y, 1)
            a_res, p_res = auc_mw(y - (slp * stats.rankdata(x) + ic), f11)
            lab_rows.append(dict(feature=fk, lab=lk, n=len(pids), deaths=int(f11.sum()), rho_with_lab=float(r_l), p_with_lab=float(p_l), resid_auc11=float(a_res), resid_p11=float(p_res)))
    pd.DataFrame(lab_rows).to_csv(data_out / "fig8_hospital_labs_true11.csv", index=False)

    # 500,000-resample converged in-fold preprocessed 50x5-fold CV under TRUE_PUB_11!
    pos_217["sex_m"] = (pos_217["sex"].astype(str).str.upper().str.startswith("M")).astype(float)
    fd_all = pos_217.sort_values("days").groupby("patient_id", as_index=False).first()
    pt_r_cv = pos_217[pos_217["source"] == "Real"].groupby("patient_id").agg(
        age=("age", "first"), sex_m=("sex_m", "first"), sev=("severity_ord", "first"), fat11=("fatal_pub11", "first"),
        CRP=("temp_CRP", "mean"), ALC=("temp_lymph", "mean"), HSPC=("HSPC", "mean"),
    ).assign(source="Real")
    pt_r_cv["logCRP"] = np.log2(pt_r_cv["CRP"] + 1.0)
    pt_n_cv = pos_217[pos_217["source"] == "PBISC"].groupby("patient_id").agg(
        age=("age", "first"), sex_m=("sex_m", "first"), sev=("severity_ord", "first"), fat11=("fatal_pub11", "first"),
        HSPC=("HSPC", "mean"),
    ).assign(source="PBISC")
    pt_66_cv = pd.concat([pt_r_cv[["age", "sex_m", "sev", "fat11", "HSPC", "source"]], pt_n_cv[["age", "sex_m", "sev", "fat11", "HSPC", "source"]]])
    pt_w23_cv = pt_66_cv[pt_66_cv["sev"] >= 2].copy(); pt_w23_cv["who3"] = (pt_w23_cv["sev"] == 3).astype(float)

    def fit_rint_1d(tr_v, te_v):
        nt = len(tr_v)
        zt = stats.norm.ppf((pd.Series(tr_v).rank().values - 0.5) / nt)
        ss = np.sort(tr_v)
        ql, qr = np.searchsorted(ss, te_v, side="left"), np.searchsorted(ss, te_v, side="right")
        pu = np.clip((0.5 * (ql + qr) + 0.5) / (nt + 1.0), 1.0 / (2 * (nt + 1.0)), 1.0 - 1.0 / (2 * (nt + 1.0)))
        return zt, stats.norm.ppf(pu)

    def eval_cv_500k(df_in, y_col, clin_cols, bio_col):
        d = df_in[[y_col, "source"] + clin_cols + [bio_col]].dropna().copy()
        y = d[y_col].values.astype(int)
        rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=50, random_state=42)
        pb, pe, ct = np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d))
        for tr, te in rskf.split(d, y):
            tr_d, te_d = d.iloc[tr], d.iloc[te]
            Xb_tr, Xb_te = tr_d[clin_cols].values.astype(float), te_d[clin_cols].values.astype(float)
            sc_b = StandardScaler().fit(Xb_tr)
            pb[te] += LogisticRegression(C=1.0).fit(sc_b.transform(Xb_tr), y[tr]).predict_proba(sc_b.transform(Xb_te))[:, 1]
            zt, ze = np.zeros(len(tr_d)), np.zeros(len(te_d))
            for sv in tr_d["source"].unique():
                mt, me = (tr_d["source"] == sv).values, (te_d["source"] == sv).values
                zv_t, zv_e = fit_rint_1d(tr_d.loc[mt, bio_col].values, te_d.loc[me, bio_col].values)
                zt[mt] = zv_t
                if me.any(): ze[me] = zv_e
            Xe_tr, Xe_te = np.hstack([Xb_tr, zt[:, None]]), np.hstack([Xb_te, ze[:, None]])
            sc_e = StandardScaler().fit(Xe_tr)
            pe[te] += LogisticRegression(C=1.0).fit(sc_e.transform(Xe_tr), y[tr]).predict_proba(sc_e.transform(Xe_te))[:, 1]
            ct[te] += 1
        pb /= ct; pe /= ct
        ab, ae = float(roc_auc_score(y, pb)), float(roc_auc_score(y, pe))
        i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
        n1, n0 = len(i1), len(i0)
        Kb = ((pb[i1, None] > pb[None, i0]).astype(np.float64) + 0.5 * (pb[i1, None] == pb[None, i0])) / (n1 * n0)
        Ke = ((pe[i1, None] > pe[None, i0]).astype(np.float64) + 0.5 * (pe[i1, None] == pe[None, i0])) / (n1 * n0)
        rg = np.random.default_rng(20261119)
        w1 = rg.multinomial(n1, np.ones(n1) / n1, size=500000).astype(np.float64)
        w0 = rg.multinomial(n0, np.ones(n0) / n0, size=500000).astype(np.float64)
        bs_b = ((w1 @ Kb) * w0).sum(axis=1)
        bs_e = ((w1 @ Ke) * w0).sum(axis=1)
        dbs = bs_e - bs_b
        return dict(
            n=len(d), deaths=n1, AUC_base=ab, base_lo=float(np.percentile(bs_b, 2.5)), base_hi=float(np.percentile(bs_b, 97.5)),
            AUC_ext=ae, ext_lo=float(np.percentile(bs_e, 2.5)), ext_hi=float(np.percentile(bs_e, 97.5)),
            Delta_AUC=ae - ab, delta_lo=float(np.percentile(dbs, 2.5)), delta_hi=float(np.percentile(dbs, 97.5)),
            Pr_delta_le_0=float(np.mean(dbs <= 0.0)), bootstrap_resamples=500000,
        )

    cv_rows = []
    for m_nm, df_m, y_c, cl_c, b_c in [
        ("First available draw (n=66, 11d): Age+Sex+SameDayWHO -> + T CD8 CTL", fd_all, "fatal_pub11", ["age", "sex_m", "temp_severity_ord"], "T CD8 CTL"),
        ("First draw <=7d subset (n=53, 10d): Age+Sex+SameDayWHO -> + T CD8 CTL", fd_all[fd_all["days"] <= 7], "fatal_pub11", ["age", "sex_m", "temp_severity_ord"], "T CD8 CTL"),
        ("Hospital bloods (Real n=37, 5d): Age+logCRP+ALC -> + HSPC", pt_r_cv, "fat11", ["age", "logCRP", "ALC"], "HSPC"),
        ("Routine labs only (Real n=37, 5d): logCRP+ALC -> + HSPC", pt_r_cv, "fat11", ["logCRP", "ALC"], "HSPC"),
        ("Within WHO 2-3 (n=29, 11d): Age+Sex+WHO3 -> + HSPC", pt_w23_cv, "fat11", ["age", "sex_m", "who3"], "HSPC"),
    ]:
        r_cv = eval_cv_500k(df_m, y_c, cl_c, b_c)
        r_cv["model"] = m_nm
        cv_rows.append(r_cv)
    df_cv = pd.DataFrame(cv_rows)[["model", "n", "deaths", "AUC_base", "base_lo", "base_hi", "AUC_ext", "ext_lo", "ext_hi", "Delta_AUC", "delta_lo", "delta_hi", "Pr_delta_le_0", "bootstrap_resamples"]]
    df_cv.to_csv(data_out / "fig8_clinical_CV_infold_INT_500k_true11.csv", index=False)
    print("\n[CV 500k TRUE_PUB_11]:\n", df_cv[["model", "AUC_base", "AUC_ext", "Delta_AUC", "delta_lo", "delta_hi", "Pr_delta_le_0"]].round(4).to_string(index=False))

    # Freeze INPUT_MANIFEST.csv & OUTPUT_DATA_MANIFEST.csv
    df_in_man = pd.DataFrame(read_registry).drop_duplicates(subset=["path"])
    df_in_man.to_csv(out_dir / "INPUT_MANIFEST.csv", index=False)
    out_csvs = sorted(data_out.glob("*.csv"))
    pd.DataFrame([{"file": f"data/{p.name}", "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in out_csvs]).to_csv(out_dir / "OUTPUT_DATA_MANIFEST.csv", index=False)
    print(f"\nDONE Stage 1: Froze {len(df_in_man)} upstream inputs in INPUT_MANIFEST.csv and wrote {len(out_csvs)} standalone CSVs to {data_out}.")


if __name__ == "__main__":
    main()
