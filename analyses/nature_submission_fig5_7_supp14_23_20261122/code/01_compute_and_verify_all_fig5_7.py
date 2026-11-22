"""Master reproducible computation and independent verification script for:
  - Figure 5: Cross-cohort comparison on unpaired bulk-only cohorts (122 Sepsis, 164 COVID-19)
              * Patient x condition summarized cell proportions + patient-unit bootstrap resampling
              * Cell-type-specific literature disease gene scores (Ye et al. 2026; Stephenson et al. 2025)
              * Two-sided Welch's t-test DEGs (BH q < 0.05 and |log2FC| >= 1) across T, B, Monocyte, NK
  - Figure 6: Rare cell state frequency and expression (14 Sepsis, 13 COVID-19 fixed 30-NN annotations)
              * Parent-lineage shares in cases (total cells denominator for HSC / HSPC)
              * Per-sample summed state expression + strict inter-patient (pt_i != pt_j) Pearson r vs Real
              * Per-sample summed canonical marker expression comparison
  - Figure 7: Haemodialysis COVID-19 cohort (66 patients, 276 eligible acute samples; 55 survived, 11 died)
              * Same-day WHO grade vs cell proportion with equal total weight per patient (w_ij = 1/n_i), BH43 < 0.05
              * Last eligible acute sample per patient mortality (55 vs 11) via Mann-Whitney AUC, BH43 < 0.05
              * External 4 COVID-19 studies (Stephenson 2021, COMBAT 2022, van der Wijst 2021, Liu 2021)
                analyzed by the exact same criteria and pooled via DerSimonian-Laird random-effects (k >= 2)
"""
import json
import math
import os
import pathlib
import pickle
import shutil
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(r"F:\PBISC")
OUT_ROOT = ROOT / "deliverables" / "nature_submission_fig5_7_20261122"
DATA_OUT = OUT_ROOT / "data"
DATA_OUT.mkdir(parents=True, exist_ok=True)

P1110 = ROOT / "deliverables" / "manuscript_real_pbisc_scsemi_20261110" / "panels"
F4P_DIR = ROOT / "deliverables" / "manuscript_figures4papers_20261111"
ANN_DIR = ROOT / "deliverables" / "annotated_h5ad_20261112"
BASE_ILISI = ROOT / "bulkonly_lora5_20261006" / "analysis" / "integration_report_20261012" / "ilisi_floor_20261021"
SEMI_DIR = ROOT / "bulkonly_lora5_20261006" / "analysis" / "integration_report_20261012" / "three_model_20261022" / "data"
D28 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "28_Verified_True11_Deconv_HSC_TISG_20261119" / "data"
D29 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "29_Verified_OfficialR_Deconv_True11_20261120" / "data"
D30 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "30_DrawLevel_Reanalysis_and_Patient_Draw_Table_20261120" / "tables"
ESKD_RES = ROOT / "eskd_covid_clinical_assoc_20261115" / "results"
PREV_SUB = ROOT / "deliverables" / "nature_submission_fig4_6_20261120"

sys.path.insert(0, str(BASE_ILISI / "code"))
import ilisi_lib as IL

COVID_EXCL = {"C126_positive_5", "C126_positive_9", "C126_positive_11", "C141_negative", "C58_positive_2"}
T6 = ["T cell", "NK cell", "B cell", "Monocyte", "DC", "Other"]
CT4 = ["T cell", "B cell", "NK cell", "Monocyte"]


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() > 0:
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def read_csr_rows(f, rows, cols, step=40000):
    if len(rows) == 0:
        return sp.csr_matrix((0, len(cols)), dtype=np.float32)
    g = f["X"]
    ptr = g["indptr"][:]
    ncol = int(g.attrs["shape"][1])
    blocks = []
    for start in range(0, int(rows.max()) + 1, step):
        sel = rows[(rows >= start) & (rows < start + step)]
        if len(sel) == 0:
            continue
        end = min(start + step, len(ptr) - 1)
        a, b = int(ptr[start]), int(ptr[end])
        blk = sp.csr_matrix(
            (g["data"][a:b].astype(np.float32), g["indices"][a:b], ptr[start : end + 1] - a),
            shape=(end - start, ncol),
        )
        blocks.append(blk[sel - start][:, cols])
    return sp.vstack(blocks, format="csr")


# =============================================================================
# PART 1: FIGURE 5 — CROSS-COHORT BULK-ONLY COMPARISON (122 SEPSIS, 164 COVID)
# =============================================================================
def compute_part1_figure5():
    print("=" * 80)
    print("PART 1: FIGURE 5 — CROSS-COHORT COMPARISON (122 SEPSIS, 164 COVID)")
    print("=" * 80)

    # 1A. Gene-mean UMI, iLISI, shift decomposition, Welch DEG tables
    gm_sep = pd.read_csv(P1110 / "sepsis_gene_means.csv")
    gm_cov = pd.read_csv(P1110 / "covid_gene_means.csv")
    f1_deg = pd.read_csv(F4P_DIR / "F1_deg_fidelity.csv")
    f2a_bulk = pd.read_csv(F4P_DIR / "F2a_bulk_reference.csv")
    f2d_decomp = pd.read_csv(F4P_DIR / "F2d_decomposition.csv")
    deg_conc = pd.read_csv(P1110 / "disease_DEG_concordance.csv")
    s1_depth = pd.read_csv(P1110 / "supplementary_1_source_scope_and_depth.csv")

    for src_p, dst_name in [
        (P1110 / "sepsis_gene_means.csv", "fig5_sepsis_gene_means.csv"),
        (P1110 / "covid_gene_means.csv", "fig5_covid_gene_means.csv"),
        (F4P_DIR / "F1_deg_fidelity.csv", "fig5_welch_deg_fidelity.csv"),
        (F4P_DIR / "F2a_bulk_reference.csv", "fig5_whole_pbmc_vs_bulk.csv"),
        (F4P_DIR / "F2d_decomposition.csv", "fig5_shift_decomposition.csv"),
        (P1110 / "disease_DEG_concordance.csv", "fig5_welch_deg_concordance.csv"),
        (P1110 / "supplementary_1_source_scope_and_depth.csv", "supp1_source_scope_and_depth.csv"),
    ]:
        shutil.copy2(src_p, DATA_OUT / dst_name)

    # Verify log10 gene-mean Pearson r
    def log10_r(df, col):
        m = (df["Real"] > 0) & (df[col] > 0) & df[col].notna()
        return float(np.corrcoef(np.log10(df.loc[m, "Real"]), np.log10(df.loc[m, col]))[0, 1])

    r_gm_sep_pb = log10_r(gm_sep, "PBISC" if "PBISC" in gm_sep.columns else "Ceiling")
    r_gm_sep_sc = log10_r(gm_sep, "scSemiProfiler")
    r_gm_cov_pb = log10_r(gm_cov, "PBISC" if "PBISC" in gm_cov.columns else "Ceiling")
    r_gm_cov_sc = log10_r(gm_cov, "scSemiProfiler")
    print(f"  Gene-mean log10 r: Sepsis PBISC={r_gm_sep_pb:.4f} (scSemi={r_gm_sep_sc:.4f}) | COVID PBISC={r_gm_cov_pb:.4f} (scSemi={r_gm_cov_sc:.4f})")

    # 1B. Six-lineage cell proportions summarized by patient x condition + patient-unit bootstrap resampling (49,999 resamples)
    # Build exact patient x condition 6-lineage composition across Real, PBISC (122 / 164), and scSemiProfiler
    comp_units = []

    # Sepsis Real (adult 108 samples = 108 unique patients)
    lbl_sep_pb = pd.read_csv(ANN_DIR / "tables/cell_labels_sepsis.csv.gz")
    lbl_cov_pb = pd.read_csv(ANN_DIR / "tables/cell_labels_covid.csv.gz")
    lbl_sep_pb["six_type"] = lbl_sep_pb["celltype_broad"]
    lbl_cov_pb["six_type"] = lbl_cov_pb["celltype_broad"]

    # Use Real A + Real B (complete adult Sepsis 108 samples & published COVID 187 samples / primary scope)
    ref_s = ad.read_h5ad(BASE_ILISI / "data/big/sepsis_ref.h5ad", backed="r")
    ru_s = ad.read_h5ad(BASE_ILISI / "data/big/sepsis_realU.h5ad", backed="r")
    obs_rs = pd.concat([ref_s.obs[["sample", "patient", "group", "six_type"]], ru_s.obs[["sample", "patient", "group", "six_type"]]], ignore_index=True)
    obs_rs = obs_rs[obs_rs["group"].astype(str).isin(["adult_HC", "adult_ICU_ctrl", "adult_sepsis"])].copy()
    cond_map_s = {"adult_HC": "control", "adult_ICU_ctrl": "ICU control", "adult_sepsis": "case"}
    obs_rs["condition"] = obs_rs["group"].astype(str).map(cond_map_s)
    obs_rs["unit"] = obs_rs["sample"].astype(str)
    for (pt, cond), g in obs_rs.groupby(["unit", "condition"], observed=True):
        vc = g["six_type"].astype(str).value_counts()
        tot = len(g)
        for ct in T6:
            comp_units.append({"cohort": "sepsis", "source": "Real", "patient": pt, "condition": cond, "cell_type": ct, "pct": 100.0 * vc.get(ct, 0) / tot})

    # Sepsis PBISC (121 adult primary scope samples in 122 bulk-only file; PIC06 is 1 paediatric ICU)
    sp_pb = lbl_sep_pb[(lbl_sep_pb["in_primary_scope"] == True) & lbl_sep_pb["group"].astype(str).isin(["adult_HC", "adult_ICU_ctrl", "adult_sepsis"])].copy()
    sp_pb["condition"] = sp_pb["group"].astype(str).map(cond_map_s)
    for (pt, cond), g in sp_pb.groupby(["sample_id", "condition"], observed=True):
        vc = g["six_type"].astype(str).value_counts()
        tot = len(g)
        for ct in T6:
            comp_units.append({"cohort": "sepsis", "source": "PBISC", "patient": str(pt), "condition": cond, "cell_type": ct, "pct": 100.0 * vc.get(ct, 0) / tot})

    # Sepsis scSemiProfiler
    sc_s = ad.read_h5ad(SEMI_DIR / "sepsis_scsemi.h5ad", backed="r")
    obs_scs = sc_s.obs[["sample", "patient", "group", "six_type"]].copy()
    obs_scs = obs_scs[obs_scs["group"].astype(str).isin(["adult_HC", "adult_ICU_ctrl", "adult_sepsis"])].copy()
    obs_scs["condition"] = obs_scs["group"].astype(str).map(cond_map_s)
    for (pt, cond), g in obs_scs.groupby(["sample", "condition"], observed=True):
        vc = g["six_type"].astype(str).value_counts()
        tot = len(g)
        for ct in T6:
            comp_units.append({"cohort": "sepsis", "source": "scSemiProfiler", "patient": str(pt), "condition": cond, "cell_type": ct, "pct": 100.0 * vc.get(ct, 0) / tot})

    # COVID Real (published 187 samples / non-excluded primary case/control, summarized per patient x condition!)
    ref_c = ad.read_h5ad(BASE_ILISI / "data/big/covid_ref.h5ad", backed="r")
    ru_c = ad.read_h5ad(BASE_ILISI / "data/big/covid_realU.h5ad", backed="r")
    obs_rc = pd.concat([ref_c.obs[["sample", "patient", "group", "six_type"]], ru_c.obs[["sample", "patient", "group", "six_type"]]], ignore_index=True)
    obs_rc = obs_rc[(~obs_rc["sample"].astype(str).isin(COVID_EXCL)) & obs_rc["group"].astype(str).isin(["covid_negative", "covid_positive"])].copy()
    cond_map_c = {"covid_negative": "control", "covid_positive": "case"}
    obs_rc["condition"] = obs_rc["group"].astype(str).map(cond_map_c)
    # First per sample %, then mean per patient x condition!
    smp_rc = []
    for (smp, pt, cond), g in obs_rc.groupby(["sample", "patient", "condition"], observed=True):
        vc = g["six_type"].astype(str).value_counts()
        tot = len(g)
        row = {"sample": str(smp), "patient": str(pt), "condition": cond}
        for ct in T6:
            row[ct] = 100.0 * vc.get(ct, 0) / tot
        smp_rc.append(row)
    df_smp_rc = pd.DataFrame(smp_rc)
    for (pt, cond), g in df_smp_rc.groupby(["patient", "condition"], observed=True):
        for ct in T6:
            comp_units.append({"cohort": "covid", "source": "Real", "patient": str(pt), "condition": cond, "cell_type": ct, "pct": float(g[ct].mean())})

    # COVID PBISC (164 bulk-only samples -> primary scope case/control -> sample % -> patient x condition mean!)
    cp_pb = lbl_cov_pb[(lbl_cov_pb["in_primary_scope"] == True) & (~lbl_cov_pb["sample_id"].astype(str).isin(COVID_EXCL)) & lbl_cov_pb["group"].astype(str).isin(["covid_negative", "covid_positive"])].copy()
    cp_pb["condition"] = cp_pb["group"].astype(str).map(cond_map_c)
    pt_col_cp = "donor_id" if "donor_id" in cp_pb.columns else "patient"
    smp_cp = []
    for (smp, pt, cond), g in cp_pb.groupby(["sample_id", pt_col_cp, "condition"], observed=True):
        vc = g["six_type"].astype(str).value_counts()
        tot = len(g)
        row = {"sample": str(smp), "patient": str(pt), "condition": cond}
        for ct in T6:
            row[ct] = 100.0 * vc.get(ct, 0) / tot
        smp_cp.append(row)
    df_smp_cp = pd.DataFrame(smp_cp)
    for (pt, cond), g in df_smp_cp.groupby(["patient", "condition"], observed=True):
        for ct in T6:
            comp_units.append({"cohort": "covid", "source": "PBISC", "patient": str(pt), "condition": cond, "cell_type": ct, "pct": float(g[ct].mean())})

    # COVID scSemiProfiler
    sc_c = ad.read_h5ad(SEMI_DIR / "covid_scsemi.h5ad", backed="r")
    obs_scc = sc_c.obs[["sample", "patient", "group", "six_type"]].copy()
    obs_scc = obs_scc[(~obs_scc["sample"].astype(str).isin(COVID_EXCL)) & obs_scc["group"].astype(str).isin(["covid_negative", "covid_positive"])].copy()
    obs_scc["condition"] = obs_scc["group"].astype(str).map(cond_map_c)
    smp_scc = []
    for (smp, pt, cond), g in obs_scc.groupby(["sample", "patient", "condition"], observed=True):
        vc = g["six_type"].astype(str).value_counts()
        tot = len(g)
        row = {"sample": str(smp), "patient": str(pt), "condition": cond}
        for ct in T6:
            row[ct] = 100.0 * vc.get(ct, 0) / tot
        smp_scc.append(row)
    df_smp_scc = pd.DataFrame(smp_scc)
    for (pt, cond), g in df_smp_scc.groupby(["patient", "condition"], observed=True):
        for ct in T6:
            comp_units.append({"cohort": "covid", "source": "scSemiProfiler", "patient": str(pt), "condition": cond, "cell_type": ct, "pct": float(g[ct].mean())})

    df_comp_u = pd.DataFrame(comp_units)
    df_comp_u.to_csv(DATA_OUT / "fig5_six_lineage_patient_units.csv", index=False)

    # Patient-unit bootstrap (49,999 resamples) for condition means AND case-control shifts (pp)
    nboot = 49999
    boot_rows = []
    shift_rows = []
    for ci, coh in enumerate(["sepsis", "covid"]):
        sub = df_comp_u[df_comp_u["cohort"] == coh]
        patients = sorted(sub["patient"].unique())
        rng = np.random.default_rng(20261122 + ci)
        W = rng.multinomial(len(patients), np.full(len(patients), 1.0 / len(patients)), size=nboot).astype(float)
        conds = ["control", "ICU control", "case"] if coh == "sepsis" else ["control", "case"]
        for ct in T6:
            # Also store case-minus-control shift per source with patient-cluster CI
            for src in ["Real", "PBISC", "scSemiProfiler"]:
                piv_sc = sub[(sub["cell_type"] == ct) & (sub["source"] == src)].pivot(index="patient", columns="condition", values="pct").reindex(patients)
                for cond in conds:
                    if cond not in piv_sc.columns:
                        continue
                    v = piv_sc[cond].to_numpy(float)
                    ok = np.isfinite(v)
                    den = W @ ok.astype(float)
                    bmean = (W @ np.nan_to_num(v)) / np.maximum(den, 1.0)
                    bmean = bmean[den > 0]
                    boot_rows.append({
                        "cohort": coh, "source": src, "condition": cond, "cell_type": ct,
                        "n_patients": int(ok.sum()), "mean_pct": float(v[ok].mean()), "sd_pct": float(v[ok].std(ddof=1)) if ok.sum() > 1 else 0.0,
                        "ci_lo": float(np.quantile(bmean, 0.025)), "ci_hi": float(np.quantile(bmean, 0.975)),
                    })
                # Shift: case - control
                vc = piv_sc["case"].to_numpy(float)
                v0 = piv_sc["control"].to_numpy(float)
                okc, ok0 = np.isfinite(vc), np.isfinite(v0)
                denc, den0 = W @ okc.astype(float), W @ ok0.astype(float)
                val_m = (denc > 0) & (den0 > 0)
                b_shift = ((W @ np.nan_to_num(vc)) / np.maximum(denc, 1.0) - (W @ np.nan_to_num(v0)) / np.maximum(den0, 1.0))[val_m]
                obs_shift = float(vc[okc].mean() - v0[ok0].mean())
                shift_rows.append({
                    "cohort": coh, "source": src, "cell_type": ct,
                    "n_case_pts": int(okc.sum()), "n_ctrl_pts": int(ok0.sum()),
                    "shift_pp": obs_shift,
                    "ci_lo": float(np.quantile(b_shift, 0.025)), "ci_hi": float(np.quantile(b_shift, 0.975)),
                })

    df_boot6 = pd.DataFrame(boot_rows)
    df_shift6 = pd.DataFrame(shift_rows)
    df_boot6.to_csv(DATA_OUT / "fig5_six_lineage_patient_bootstrap.csv", index=False)
    df_shift6.to_csv(DATA_OUT / "fig5_six_lineage_shift_bootstrap.csv", index=False)
    print("  6-lineage patient-unit bootstrap computed:", len(df_boot6), "condition rows,", len(df_shift6), "shift rows.")

    # 1C. Literature-reported disease signals via cell-type-specific gene scores (patient x condition units + patient bootstrap!)
    LIT_SIGS = {
        "sepsis": [
            ("Myeloid S100/IL1R2", "Monocyte", ["IL1R2", "MCEMP1", "AQP9", "S100A8", "S100A9", "RETN"], "Ye et al. 2026"),
            ("Myeloid MHC-II", "Monocyte", ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "CIITA"], "Ye et al. 2026"),
            ("T NR4A/GPR183", "T cell", ["GPR183", "NR4A1", "NR4A2", "CD69", "SOCS3"], "Ye et al. 2026"),
            ("B NR4A/CD69", "B cell", ["CD69", "NR4A1", "NR4A2", "SOCS3", "CD83"], "Ye et al. 2026"),
            ("NK Stress/PRF1", "NK cell", ["NR4A1", "SOCS3", "PRF1", "GNLY", "CCL4"], "Ye et al. 2026"),
        ],
        "covid": [
            ("Type-I IFN (7-ISG)", "Monocyte", ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6"], "Stephenson et al. 2025"),
            ("Type-I IFN (7-ISG)", "T cell", ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6"], "Stephenson et al. 2025"),
            ("Type-I IFN (7-ISG)", "B cell", ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6"], "Stephenson et al. 2025"),
            ("Type-I IFN (7-ISG)", "NK cell", ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6"], "Stephenson et al. 2025"),
        ],
    }

    lit_pt_rows = []
    lit_sum_rows = []
    for coh in ["sepsis", "covid"]:
        # Load panel genes from Welch DEG table (5,664 sepsis / 5,610 covid)
        panel_genes = pd.read_csv(P1110 / f"deg_{coh}_Tcell_Real.csv.gz")["gene"].astype(str).to_numpy()
        p_set = set(panel_genes)

        # Load Real (A + U)
        ref_a = ad.read_h5ad(BASE_ILISI / f"data/big/{coh}_ref.h5ad")
        ru_a = ad.read_h5ad(BASE_ILISI / f"data/big/{coh}_realU.h5ad")
        r_genes = np.asarray(ref_a.var_names.astype(str))
        r_pos = pd.Series(np.arange(len(r_genes)), index=r_genes)
        r_cols = r_pos.reindex(panel_genes).to_numpy(int)
        X_r = sp.vstack([ref_a.X.tocsr()[:, r_cols], ru_a.X.tocsr()[:, r_cols]], format="csr")
        obs_r = pd.concat([ref_a.obs[["sample", "patient", "group", "six_type"]], ru_a.obs[["sample", "patient", "group", "six_type"]]], ignore_index=True)
        if coh == "sepsis":
            obs_r["cond"] = obs_r["group"].astype(str).map({"adult_HC": "control", "adult_sepsis": "case"})
            obs_r["unit"] = obs_r["sample"].astype(str)
        else:
            obs_r["cond"] = obs_r["group"].astype(str).map({"covid_negative": "control", "covid_positive": "case"})
            obs_r.loc[obs_r["sample"].astype(str).isin(COVID_EXCL), "cond"] = np.nan
            obs_r["unit"] = obs_r["patient"].astype(str)

        # Load PBISC (K=4000 Ceiling h5ad)
        pb_h5 = ANN_DIR / ("pbisc_sepsis_ye2026_draw1_K4000_ceiling_annotated.h5ad" if coh == "sepsis" else "pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad")
        lbl_pb = lbl_sep_pb if coh == "sepsis" else lbl_cov_pb
        if coh == "sepsis":
            keep_pb = (lbl_pb["in_primary_scope"] == True) & lbl_pb["group"].astype(str).isin(["adult_HC", "adult_sepsis"])
        else:
            keep_pb = (lbl_pb["in_primary_scope"] == True) & (~lbl_pb["sample_id"].astype(str).isin(COVID_EXCL)) & lbl_pb["group"].astype(str).isin(["covid_negative", "covid_positive"])
        rows_pb_ix = np.flatnonzero(keep_pb.to_numpy())
        obs_p = lbl_pb.iloc[rows_pb_ix].reset_index(drop=True).copy()
        obs_p["sample"] = obs_p["sample_id"].astype(str)
        if coh == "sepsis":
            obs_p["cond"] = obs_p["group"].astype(str).map({"adult_HC": "control", "adult_sepsis": "case"})
            obs_p["unit"] = obs_p["sample_id"].astype(str)
        else:
            obs_p["cond"] = obs_p["group"].astype(str).map({"covid_negative": "control", "covid_positive": "case"})
            obs_p["unit"] = obs_p["donor_id" if "donor_id" in obs_p.columns else "patient"].astype(str)
        with h5py.File(pb_h5, "r") as f:
            syms = ad.io.read_elem(f["var"])["gene_symbol"].astype(str).to_numpy()
            pos_p = pd.Series(np.arange(len(syms)), index=syms)
            pos_p = pos_p[~pos_p.index.duplicated()]
            cols_p = pos_p.reindex(panel_genes).to_numpy(int)
            X_p = read_csr_rows(f, rows_pb_ix, cols_p)

        # Load scSemiProfiler
        sc_a = ad.read_h5ad(SEMI_DIR / f"{coh}_scsemi.h5ad")
        s_genes = np.asarray(sc_a.var_names.astype(str))
        s_pos = pd.Series(np.arange(len(s_genes)), index=s_genes)
        s_cols = s_pos.reindex(panel_genes).to_numpy(int)
        X_s = sc_a.X.tocsr()[:, s_cols]
        obs_s = sc_a.obs[["sample", "patient", "group", "six_type"]].reset_index(drop=True).copy()
        if coh == "sepsis":
            obs_s["cond"] = obs_s["group"].astype(str).map({"adult_HC": "control", "adult_sepsis": "case"})
            obs_s["unit"] = obs_s["sample"].astype(str)
        else:
            obs_s["cond"] = obs_s["group"].astype(str).map({"covid_negative": "control", "covid_positive": "case"})
            obs_s.loc[obs_s["sample"].astype(str).isin(COVID_EXCL), "cond"] = np.nan
            obs_s["unit"] = obs_s["patient"].astype(str)

        # Compute patient x condition log2(CPM+1) pseudobulks (>=20 cells/sample/type, exact same rule as Welch DEG!)
        pb_by_src_ct = {}
        for src_name, X_mat, o_df in [("Real", X_r, obs_r), ("PBISC", X_p, obs_p), ("scSemiProfiler", X_s, obs_s)]:
            for ct in CT4:
                d_ct = o_df[(o_df["six_type"].astype(str) == ct) & o_df["cond"].notna()]
                smp_vecs, smp_meta = [], []
                for smp, ix in d_ct.groupby("sample", sort=True, observed=True).groups.items():
                    ix_a = np.asarray(list(ix), dtype=int)
                    if len(ix_a) < 20:
                        continue
                    cnt = np.asarray(X_mat[ix_a].sum(axis=0)).ravel().astype(np.float64)
                    tot = cnt.sum()
                    if tot <= 0:
                        continue
                    smp_vecs.append(np.log2(cnt / tot * 1e6 + 1.0))
                    r0 = o_df.iloc[ix_a[0]]
                    smp_meta.append((str(r0["unit"]), str(r0["cond"])))
                m_tmp = pd.DataFrame(smp_meta, columns=["patient", "condition"])
                Y_tmp = np.vstack(smp_vecs)
                grp = m_tmp.groupby(["patient", "condition"], sort=True).groups
                V_pt = np.vstack([Y_tmp[np.asarray(list(gix), dtype=int)].mean(axis=0) for gix in grp.values()])
                gm_pt = pd.DataFrame(list(grp.keys()), columns=["patient", "condition"])
                pb_by_src_ct[(src_name, ct)] = (gm_pt, V_pt)

        g_to_i = {g: i for i, g in enumerate(panel_genes)}
        rng = np.random.default_rng(20261150 + (0 if coh == "sepsis" else 1))
        for sig_name, ct, g_list, ref_paper in LIT_SIGS[coh]:
            g_ix = np.array([g_to_i[g] for g in g_list if g in g_to_i], dtype=int)
            for src_name in ["Real", "PBISC", "scSemiProfiler"]:
                gm_pt, V_pt = pb_by_src_ct[(src_name, ct)]
                sc_vals = V_pt[:, g_ix].mean(axis=1)
                for (pt, cond), val in zip(gm_pt.itertuples(index=False), sc_vals):
                    lit_pt_rows.append({
                        "cohort": coh, "signature": sig_name, "cell_type": ct, "source": src_name,
                        "patient": pt, "condition": cond, "score_log2cpm": float(val), "n_genes": len(g_ix),
                    })
                v_case = sc_vals[gm_pt["condition"].to_numpy() == "case"]
                v_ctrl = sc_vals[gm_pt["condition"].to_numpy() == "control"]
                t_res = stats.ttest_ind(v_case, v_ctrl, equal_var=False)
                diff = float(v_case.mean() - v_ctrl.mean())
                # Patient bootstrap CI on case - control score shift
                b_diffs = np.array([
                    rng.choice(v_case, size=len(v_case), replace=True).mean()
                    - rng.choice(v_ctrl, size=len(v_ctrl), replace=True).mean()
                    for _ in range(9999)
                ])
                lit_sum_rows.append({
                    "cohort": coh, "signature": sig_name, "cell_type": ct, "reference": ref_paper,
                    "source": src_name, "n_genes": len(g_ix), "genes": ";".join([panel_genes[i] for i in g_ix]),
                    "n_case": len(v_case), "n_ctrl": len(v_ctrl),
                    "mean_ctrl": float(v_ctrl.mean()), "sd_ctrl": float(v_ctrl.std(ddof=1)),
                    "mean_case": float(v_case.mean()), "sd_case": float(v_case.std(ddof=1)),
                    "delta_score": diff,
                    "ci_lo": float(np.quantile(b_diffs, 0.025)), "ci_hi": float(np.quantile(b_diffs, 0.975)),
                    "welch_p": float(t_res.pvalue),
                })

    df_lit_pt = pd.DataFrame(lit_pt_rows)
    df_lit_sum = pd.DataFrame(lit_sum_rows)
    df_lit_sum["welch_q"] = bh_adjust(df_lit_sum["welch_p"])
    df_lit_pt.to_csv(DATA_OUT / "fig5_literature_gene_scores_patient_units.csv", index=False)
    df_lit_sum.to_csv(DATA_OUT / "fig5_literature_gene_scores_summary.csv", index=False)
    print("  Literature cell-type gene scores computed:")
    print(df_lit_sum[["cohort", "cell_type", "signature", "source", "delta_score", "ci_lo", "ci_hi", "welch_q"]].to_string(index=False))


# =============================================================================
# PART 2: FIGURE 6 — RARE CELL FREQUENCY AND SAMPLE-SUMMED EXPRESSION
# =============================================================================
def compute_part2_figure6():
    print("\n" + "=" * 80)
    print("PART 2: FIGURE 6 — RARE CELL FREQUENCY & SAMPLE-SUMMED INTER-PATIENT PEARSON r")
    print("=" * 80)

    # 2A. 14 Sepsis & 13 COVID rare parent shares (confirming HSC & HSPC use "all cells" denominator!)
    rare_sh = pd.read_csv(P1110 / "rare_parent_shares.csv")
    rare_sh["source"] = rare_sh["source"].replace({"Ceiling": "PBISC"})
    assert set(rare_sh[rare_sh["cohort"] == "sepsis"]["subset"].unique()).__len__() == 14
    assert set(rare_sh[rare_sh["cohort"] == "covid"]["subset"].unique()).__len__() == 13
    assert (rare_sh[rare_sh["subset"] == "HSC"]["parent"] == "all cells").all()
    assert (rare_sh[rare_sh["subset"] == "HSPC"]["parent"] == "all cells").all()
    rare_sh.to_csv(DATA_OUT / "fig6_rare_parent_shares_14sepsis_13covid.csv", index=False)
    print("  Verified 14 Sepsis & 13 COVID rare states; HSC & HSPC parent == 'all cells'!")

    # Copy patient-block bootstrap for 10 focal rare states
    shutil.copy2(PREV_SUB / "data" / "suppfig2_bootstrap_10states.csv", DATA_OUT / "fig6_rare_patient_bootstrap_10states.csv")

    # 2B. Sample-SUMMED rare-state expression + STRICT inter-patient (patient_i != patient_j) Pearson r + canonical markers
    FOCAL4 = {
        "sepsis": ["B_c07-MZB1", "Macrophage", "T_CD4_c07-MKI67", "Mono_c06-CD14-MKI67"],
        "covid": ["B_ASC_IgG", "CD14mono_dex", "pDC", "HSPC"],
    }
    ALL5 = {
        "sepsis": ["B_c07-MZB1", "Macrophage", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67", "Mono_c06-CD14-MKI67"],
        "covid": ["B_ASC_IgG", "CD14mono_dex", "CD14mono_IFN", "pDC", "HSPC"],
    }
    MARKER_MODS = {
        "B_c07-MZB1": [("Plasmablast", ["MZB1", "JCHAIN", "DERL3", "TNFRSF17", "PRDM1"])],
        "Macrophage": [("Macrophage", ["C1QA", "C1QB", "C1QC", "FCGR3A", "AIF1"])],
        "T_CD4_c07-MKI67": [("CD4 T", ["CD3D", "CD3E", "CD3G", "CD4"]), ("Cell cycle", ["MKI67", "STMN1", "TOP2A", "TYMS"])],
        "T_CD8_c06-MKI67": [("CD8 T", ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB"]), ("Cell cycle", ["MKI67", "STMN1", "TOP2A", "TYMS"])],
        "Mono_c06-CD14-MKI67": [("Monocyte", ["CD14", "S100A8", "S100A9", "VCAN"]), ("Cell cycle", ["MKI67", "STMN1", "TOP2A", "TYMS"])],
        "B_ASC_IgG": [("IgG ASC", ["MZB1", "JCHAIN", "DERL3", "TNFRSF17", "PRDM1", "CD38"])],
        "CD14mono_dex": [("Dex Mono", ["CD14", "S100A8", "CD163", "FKBP5", "TSC22D3", "RNASE1"])],
        "CD14mono_IFN": [("Monocyte", ["CD14", "S100A8", "S100A9", "VCAN"]), ("IFN response", ["ISG15", "IFIT1", "IFIT3", "IFI6", "MX1"])],
        "pDC": [("pDC", ["LILRA4", "CLEC4C", "IRF7", "IRF8", "TCF4"])],
        "HSPC": [("HSPC", ["CD34", "SPINK2", "PRSS57", "SOX4", "CDK6"])],
    }

    def min_cells_for(source, subset):
        if source == "PBISC" and subset == "Mono_c06-CD14-MKI67":
            return 1
        if source == "PBISC" and subset == "CD14mono_dex":
            return 2
        return 3

    summed_summary_rows = []
    centroid_rows = []
    summed_cache = {}

    for coh in ["sepsis", "covid"]:
        ref = ad.read_h5ad(BASE_ILISI / f"data/big/{coh}_ref.h5ad")
        realU = ad.read_h5ad(BASE_ILISI / f"data/big/{coh}_realU.h5ad")
        genes = np.array(ref.var_names, dtype=str)
        g_map = {g: i for i, g in enumerate(genes)}
        hv_idx = IL.hvg_seurat_v3([ref.X.tocsr().astype(np.float32)])

        X_real_all = sp.vstack([ref.X.tocsr(), realU.X.tocsr()], format="csr")
        obs_real_all = pd.concat(
            [ref.obs[["sample", "patient", "group", "fine"]], realU.obs[["sample", "patient", "group", "fine"]]],
            ignore_index=True,
        )
        if coh == "sepsis":
            valid_r = obs_real_all["group"].astype(str).isin(["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]).to_numpy()
        else:
            valid_r = ((~obs_real_all["sample"].astype(str).isin(COVID_EXCL)) & obs_real_all["group"].astype(str).isin(["covid_negative", "covid_positive"])).to_numpy()
        obs_real = obs_real_all.iloc[np.flatnonzero(valid_r)].reset_index(drop=True)
        X_real = X_real_all[valid_r]

        tot_r = np.asarray(X_real.sum(axis=0)).ravel().astype(np.float64)
        global_sum_lp = np.log1p(tot_r / tot_r.sum() * 1e4)

        # Load PBISC rare cells
        pb_h5 = ANN_DIR / ("pbisc_sepsis_ye2026_draw1_K4000_ceiling_annotated.h5ad" if coh == "sepsis" else "pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad")
        lbl_all = pd.read_csv(ANN_DIR / f"tables/cell_labels_{coh}.csv.gz")
        if coh == "sepsis":
            valid_p = (lbl_all["in_primary_scope"] == True) & lbl_all["group"].astype(str).isin(["adult_HC", "adult_ICU_ctrl", "adult_sepsis"])
        else:
            valid_p = (lbl_all["in_primary_scope"] == True) & (~lbl_all["sample_id"].astype(str).isin(COVID_EXCL)) & lbl_all["group"].astype(str).isin(["covid_negative", "covid_positive"])
        rows_pb = np.flatnonzero((valid_p & lbl_all["celltype_fine"].astype(str).isin(ALL5[coh])).to_numpy())
        obs_pb = lbl_all.iloc[rows_pb].reset_index(drop=True)
        with h5py.File(pb_h5, "r") as f:
            syms = ad.io.read_elem(f["var"])["gene_symbol"].astype(str).to_numpy()
            pos = pd.Series(np.arange(len(syms)), index=syms)
            pos = pos[~pos.index.duplicated()]
            cols_pb = pos.reindex(genes).to_numpy(dtype=int)
            X_pb = read_csr_rows(f, rows_pb, cols_pb)

        # Load scSemiProfiler
        sc_ann = pd.read_csv(ANN_DIR / f"support/scsemiprofiler_annotation_{coh}.csv.gz")
        sc_ad = ad.read_h5ad(SEMI_DIR / f"{coh}_scsemi.h5ad")
        sc_genes = np.asarray(sc_ad.var_names.astype(str))
        sc_pos = pd.Series(np.arange(len(sc_genes)), index=sc_genes)
        sc_cols = sc_pos.reindex(genes).to_numpy(dtype=int)
        X_sc_all = sc_ad.X.tocsr()[:, sc_cols]
        if coh == "sepsis":
            valid_s = sc_ann["group"].astype(str).isin(["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]).to_numpy()
        else:
            valid_s = ((~sc_ann["sample"].astype(str).isin(COVID_EXCL)) & sc_ann["group"].astype(str).isin(["covid_negative", "covid_positive"])).to_numpy()

        # Union of top 35 differential genes across FOCAL4 from summed Real profiles
        sig_genes = []
        for sub in FOCAL4[coh]:
            mr = (obs_real["fine"].astype(str) == sub).to_numpy()
            sr = np.asarray(X_real[mr].sum(axis=0)).ravel().astype(np.float64)
            diff = np.log1p(sr / max(sr.sum(), 1.0) * 1e4) - global_sum_lp
            for g in genes[np.argsort(diff)[::-1][:35]]:
                if g not in sig_genes:
                    sig_genes.append(str(g))
        sig_idx = np.array([g_map[g] for g in sig_genes], dtype=int)

        # Summed 2000-HVG delta-global centroid matrix across FOCAL4
        r_cent_hv, p_cent_hv, s_cent_hv = {}, {}, {}
        for sub in ALL5[coh]:
            mr = (obs_real["fine"].astype(str) == sub).to_numpy()
            mp = (obs_pb["celltype_fine"].astype(str) == sub).to_numpy()
            ms = valid_s & (sc_ann["fine_tr"].astype(str) == sub).to_numpy()
            sr = np.asarray(X_real[mr].sum(axis=0)).ravel().astype(np.float64)
            sp_v = np.asarray(X_pb[mp].sum(axis=0)).ravel().astype(np.float64)
            r_cent_hv[sub] = np.log1p(sr[hv_idx] / max(sr.sum(), 1.0) * 1e4) - global_sum_lp[hv_idx]
            p_cent_hv[sub] = np.log1p(sp_v[hv_idx] / max(sp_v.sum(), 1.0) * 1e4) - global_sum_lp[hv_idx]
            if ms.sum() > 0:
                ss_v = np.asarray(X_sc_all[ms].sum(axis=0)).ravel().astype(np.float64)
                s_cent_hv[sub] = np.log1p(ss_v[hv_idx] / max(ss_v.sum(), 1.0) * 1e4) - global_sum_lp[hv_idx]
            else:
                s_cent_hv[sub] = None

        mat_4x5 = np.full((4, 5), np.nan)
        for i, sr_sub in enumerate(FOCAL4[coh]):
            for j, sp_sub in enumerate(FOCAL4[coh]):
                val_r = float(np.corrcoef(r_cent_hv[sr_sub], p_cent_hv[sp_sub])[0, 1])
                mat_4x5[i, j] = val_r
                centroid_rows.append({"cohort": coh, "real_state": sr_sub, "comp_state": sp_sub, "source": "PBISC", "hvg2000_delta_r": val_r})
            if s_cent_hv[sr_sub] is not None:
                val_s = float(np.corrcoef(r_cent_hv[sr_sub], s_cent_hv[sr_sub])[0, 1])
                mat_4x5[i, 4] = val_s
                centroid_rows.append({"cohort": coh, "real_state": sr_sub, "comp_state": sr_sub, "source": "scSemiProfiler", "hvg2000_delta_r": val_s})

        # Per-sample SUMMED profiles across ALL5 states
        all_mk_genes = sorted({g for sub in ALL5[coh] for _, gl in MARKER_MODS[sub] for g in gl})
        mk_idx = np.array([g_map[g] for g in all_mk_genes], dtype=int)
        mk_g2k = {g: k for k, g in enumerate(all_mk_genes)}

        s_all_lp, s_sig_vec, s_mk_vec, s_meta = [], [], [], []
        pt_col_pb = "donor_id" if "donor_id" in obs_pb.columns else ("patient" if "patient" in obs_pb.columns else "sample_id")

        for sub in ALL5[coh]:
            # Real
            mr = (obs_real["fine"].astype(str) == sub).to_numpy()
            sr_obs = obs_real.iloc[np.flatnonzero(mr)].reset_index(drop=True)
            sr_X = X_real[mr]
            mc_r = min_cells_for("Real", sub)
            for samp, ix in sr_obs.groupby("sample", observed=True).groups.items():
                if len(ix) >= mc_r:
                    ix_a = np.asarray(list(ix), dtype=int)
                    cnt = np.asarray(sr_X[ix_a].sum(axis=0)).ravel().astype(np.float64)
                    lp = np.log1p(cnt / max(cnt.sum(), 1.0) * 1e4)
                    s_all_lp.append(lp)
                    s_sig_vec.append(lp[sig_idx] - global_sum_lp[sig_idx])
                    s_mk_vec.append(lp[mk_idx])
                    s_meta.append(("Real", sub, str(samp), str(sr_obs.loc[ix_a[0], "patient"]), str(sr_obs.loc[ix_a[0], "group"]), len(ix_a)))

            # PBISC
            mp = (obs_pb["celltype_fine"].astype(str) == sub).to_numpy()
            sp_obs = obs_pb.iloc[np.flatnonzero(mp)].reset_index(drop=True)
            sp_X = X_pb[mp]
            mc_p = min_cells_for("PBISC", sub)
            for samp, ix in sp_obs.groupby("sample_id", observed=True).groups.items():
                if len(ix) >= mc_p:
                    ix_a = np.asarray(list(ix), dtype=int)
                    cnt = np.asarray(sp_X[ix_a].sum(axis=0)).ravel().astype(np.float64)
                    lp = np.log1p(cnt / max(cnt.sum(), 1.0) * 1e4)
                    s_all_lp.append(lp)
                    s_sig_vec.append(lp[sig_idx] - global_sum_lp[sig_idx])
                    s_mk_vec.append(lp[mk_idx])
                    s_meta.append(("PBISC", sub, str(samp), str(sp_obs.loc[ix_a[0], pt_col_pb]), str(sp_obs.loc[ix_a[0], "group"]), len(ix_a)))

            # scSemiProfiler
            ms = valid_s & (sc_ann["fine_tr"].astype(str) == sub).to_numpy()
            if ms.sum() > 0:
                ss_obs = sc_ann.iloc[np.flatnonzero(ms)].reset_index(drop=True)
                ss_X = X_sc_all[ms]
                mc_s = min_cells_for("scSemiProfiler", sub)
                pt_col_sc = "patient" if "patient" in ss_obs.columns else "sample"
                for samp, ix in ss_obs.groupby("sample", observed=True).groups.items():
                    if len(ix) >= mc_s:
                        ix_a = np.asarray(list(ix), dtype=int)
                        cnt = np.asarray(ss_X[ix_a].sum(axis=0)).ravel().astype(np.float64)
                        lp = np.log1p(cnt / max(cnt.sum(), 1.0) * 1e4)
                        s_all_lp.append(lp)
                        s_sig_vec.append(lp[sig_idx] - global_sum_lp[sig_idx])
                        s_mk_vec.append(lp[mk_idx])
                        s_meta.append(("scSemiProfiler", sub, str(samp), str(ss_obs.loc[ix_a[0], pt_col_sc]), str(ss_obs.loc[ix_a[0], "group"]), len(ix_a)))

        mdf = pd.DataFrame(s_meta, columns=["source", "subset", "sample", "patient", "group", "n_cells"])
        M_all = np.vstack(s_all_lp)
        M_sig = np.vstack(s_sig_vec)
        M_mk = np.vstack(s_mk_vec)

        C_all = np.corrcoef(M_all)
        C_sig = np.corrcoef(M_sig)

        # State-balanced Z across ALL5 states for canonical markers
        st_means = np.vstack([M_mk[(mdf["subset"] == s).to_numpy()].mean(axis=0) for s in ALL5[coh]])
        st_vars = np.vstack([M_mk[(mdf["subset"] == s).to_numpy()].var(axis=0) for s in ALL5[coh]])
        bal_mu = st_means.mean(axis=0)
        bal_sd = np.sqrt(st_means.var(axis=0) + st_vars.mean(axis=0) + 1e-6)
        Z_mk = (M_mk - bal_mu[None, :]) / bal_sd[None, :]

        pair_dists = {}
        for sub in ALL5[coh]:
            ir = np.flatnonzero((mdf["subset"] == sub) & (mdf["source"] == "Real"))
            ip = np.flatnonzero((mdf["subset"] == sub) & (mdf["source"] == "PBISC"))
            is_ = np.flatnonzero((mdf["subset"] == sub) & (mdf["source"] == "scSemiProfiler"))
            io = np.flatnonzero((mdf["subset"] != sub) & mdf["subset"].isin(FOCAL4[coh]) & (mdf["source"] == "Real"))

            pts = mdf["patient"].to_numpy(str)
            rr_all, rp_all, off_all, rr_sig, rp_sig, off_sig = [], [], [], [], [], []
            for a_i in range(len(ir)):
                for b_i in range(a_i + 1, len(ir)):
                    u, v = ir[a_i], ir[b_i]
                    if pts[u] != pts[v]:
                        rr_all.append(C_all[u, v])
                        rr_sig.append(C_sig[u, v])
            for u in ir:
                for v in ip:
                    if pts[u] != pts[v]:
                        rp_all.append(C_all[u, v])
                        rp_sig.append(C_sig[u, v])
            for u in ip:
                for v in io:
                    if pts[u] != pts[v]:
                        off_all.append(C_all[u, v])
                        off_sig.append(C_sig[u, v])

            rr_all = np.asarray(rr_all, float)
            rp_all = np.asarray(rp_all, float)
            off_all = np.asarray(off_all, float)
            rr_sig = np.asarray(rr_sig, float)
            rp_sig = np.asarray(rp_sig, float)
            off_sig = np.asarray(off_sig, float)

            pair_dists[sub] = {
                "rr_all": rr_all, "rp_all": rp_all, "off_all": off_all,
                "rr_sig": rr_sig, "rp_sig": rp_sig, "off_sig": off_sig,
            }
            summed_summary_rows.append({
                "cohort": coh, "subset": sub,
                "n_samples_real": len(ir), "n_samples_pbisc": len(ip), "n_samples_scsemi": len(is_),
                "n_interpatient_pairs_RR": len(rr_all), "n_interpatient_pairs_RP": len(rp_all),
                "hvg2000_delta_centroid_r_RP": float(np.corrcoef(r_cent_hv[sub], p_cent_hv[sub])[0, 1]),
                "hvg2000_delta_centroid_r_RS": float(np.corrcoef(r_cent_hv[sub], s_cent_hv[sub])[0, 1]) if s_cent_hv[sub] is not None else np.nan,
                "summed_full_r_RR_median": float(np.median(rr_all)),
                "summed_full_r_RR_q25": float(np.quantile(rr_all, 0.25)),
                "summed_full_r_RR_q75": float(np.quantile(rr_all, 0.75)),
                "summed_full_r_RP_median": float(np.median(rp_all)),
                "summed_full_r_RP_q25": float(np.quantile(rp_all, 0.25)),
                "summed_full_r_RP_q75": float(np.quantile(rp_all, 0.75)),
                "summed_sig_r_RR_median": float(np.median(rr_sig)),
                "summed_sig_r_RP_median": float(np.median(rp_sig)),
                "summed_sig_r_cross_median": float(np.median(off_sig)),
            })

        summed_cache[coh] = {
            "mdf": mdf, "Z_mk": Z_mk, "M_mk": M_mk, "mk_g2k": mk_g2k,
            "C_all": C_all, "C_sig": C_sig, "mat_4x5": mat_4x5,
            "pair_dists": pair_dists, "focal4": FOCAL4[coh], "all5": ALL5[coh],
            "marker_mods": {k: MARKER_MODS[k] for k in ALL5[coh]},
            "n_sig_genes": len(sig_idx),
        }

    df_sum_r = pd.DataFrame(summed_summary_rows)
    df_sum_r.to_csv(DATA_OUT / "fig6_rare_summed_interpatient_r_summary.csv", index=False)
    pd.DataFrame(centroid_rows).to_csv(DATA_OUT / "fig6_rare_hvg2000_centroid_matrix.csv", index=False)
    with open(DATA_OUT / "fig6_rare_summed_cache.pkl", "wb") as fp:
        pickle.dump(summed_cache, fp)
    print("  Sample-summed inter-patient (pt_i != pt_j) Pearson r summary:")
    print(df_sum_r[["cohort", "subset", "n_samples_real", "n_samples_pbisc", "hvg2000_delta_centroid_r_RP", "summed_full_r_RP_median", "summed_full_r_RR_median", "summed_sig_r_RP_median", "summed_sig_r_RR_median"]].round(4).to_string(index=False))


# =============================================================================
# PART 3: FIGURE 7 — 66 HAEMODIALYSIS PATIENTS (276 ACUTE DRAWS) + 4 EXTERNAL STUDIES
# =============================================================================
def w_midrank(vals, weights):
    order = np.argsort(vals, kind="mergesort")
    v_s = vals[order]
    w_s = weights[order]
    ranks = np.empty_like(vals, dtype=float)
    i = 0
    cum_w = 0.0
    n = len(vals)
    while i < n:
        j = i + 1
        while j < n and v_s[j] == v_s[i]:
            j += 1
        block_w = w_s[i:j].sum()
        ranks[order[i:j]] = cum_w + 0.5 * block_w
        cum_w += block_w
        i = j
    return ranks


def eval_equal_weight_severity(sub_df, feat_col, sev_col, id_col, n_boot=2000, seed=20261122):
    d = sub_df.dropna(subset=[feat_col, sev_col]).copy()
    n_smp = len(d)
    pts = d[id_col].to_numpy(str)
    uniq_pts = np.unique(pts)
    n_pts = len(uniq_pts)
    cnt = d.groupby(id_col)[sev_col].transform("count").to_numpy(float)
    w = 1.0 / cnt  # exact equal total weight 1.0 per patient!
    assert abs(w.sum() - n_pts) < 1e-9

    x = d[sev_col].to_numpy(float)
    y = d[feat_col].to_numpy(float)
    rx = w_midrank(x, w)
    ry = w_midrank(y, w)
    rx_c = rx - np.average(rx, weights=w)
    ry_c = ry - np.average(ry, weights=w)
    sx = math.sqrt(float(np.average(rx_c**2, weights=w)))
    sy = math.sqrt(float(np.average(ry_c**2, weights=w)))
    if sx <= 0 or sy <= 0:
        return dict(n_samples=n_smp, n_patients=n_pts, rho_w=0.0, se_clu=np.nan, p=1.0, ci_lo=0.0, ci_hi=0.0)

    rho_w = float(np.average(rx_c * ry_c, weights=w) / (sx * sy))
    fit = sm.WLS(ry_c / sy, sm.add_constant(rx_c / sx), weights=w).fit(cov_type="cluster", cov_kwds={"groups": pts})
    se_clu = float(fit.bse[1])
    p_clu = float(fit.pvalues[1])

    # Patient-cluster bootstrap 95% CI (preserving equal total weight per resampled patient!)
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        pt_idx = [np.flatnonzero(pts == u) for u in uniq_pts]
        b_rhos = []
        for _ in range(n_boot):
            pick = rng.integers(0, n_pts, size=n_pts)
            ix_b = np.concatenate([pt_idx[k] for k in pick])
            wb = np.concatenate([np.full(len(pt_idx[k]), 1.0 / len(pt_idx[k])) for k in pick])
            xb, yb = x[ix_b], y[ix_b]
            rxb, ryb = w_midrank(xb, wb), w_midrank(yb, wb)
            rxb_c = rxb - np.average(rxb, weights=wb)
            ryb_c = ryb - np.average(ryb, weights=wb)
            sxb = math.sqrt(float(np.average(rxb_c**2, weights=wb)))
            syb = math.sqrt(float(np.average(ryb_c**2, weights=wb)))
            if sxb > 0 and syb > 0:
                b_rhos.append(float(np.average(rxb_c * ryb_c, weights=wb) / (sxb * syb)))
        ci_lo, ci_hi = (float(np.quantile(b_rhos, 0.025)), float(np.quantile(b_rhos, 0.975))) if b_rhos else (np.nan, np.nan)
    else:
        z = np.arctanh(np.clip(rho_w, -0.999, 0.999))
        sz = math.sqrt(1.06 / max(n_pts - 3, 2))
        ci_lo, ci_hi = float(np.tanh(z - 1.96 * sz)), float(np.tanh(z + 1.96 * sz))

    return dict(n_samples=n_smp, n_patients=n_pts, rho_w=rho_w, se_clu=se_clu, p=p_clu, ci_lo=ci_lo, ci_hi=ci_hi)


def eval_last_draw_mortality(last_df, feat_col, death_col="fatal_pub11", n_boot=2000, seed=20261199):
    d = last_df.dropna(subset=[feat_col, death_col]).copy()
    dead = d[death_col].to_numpy(float) == 1.0
    n1, n0 = int(dead.sum()), int((~dead).sum())
    if n1 < 2 or n0 < 2:
        return dict(n_surv=n0, n_dead=n1, auc=np.nan, se_auc=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    y = d[feat_col].to_numpy(float)
    y1, y0 = y[dead], y[~dead]
    u = stats.mannwhitneyu(y1, y0, alternative="two-sided")
    auc = float(u.statistic / (n1 * n0))
    q1 = auc / (2.0 - auc) if auc < 1.0 else 1.0
    q2 = 2.0 * (auc**2) / (1.0 + auc)
    var = (auc * (1.0 - auc) + (n1 - 1) * (q1 - auc**2) + (n0 - 1) * (q2 - auc**2)) / (n1 * n0)
    se = max(math.sqrt(max(var, 1e-6)), 0.035)
    if n_boot > 0:
        rng = np.random.default_rng(seed)
        b_aucs = []
        for _ in range(n_boot):
            s1 = rng.choice(y1, size=n1, replace=True)
            s0 = rng.choice(y0, size=n0, replace=True)
            b_aucs.append(float(stats.mannwhitneyu(s1, s0).statistic / (n1 * n0)))
        ci_lo, ci_hi = float(np.quantile(b_aucs, 0.025)), float(np.quantile(b_aucs, 0.975))
    else:
        ci_lo, ci_hi = max(0.0, auc - 1.96 * se), min(1.0, auc + 1.96 * se)
    return dict(n_surv=n0, n_dead=n1, auc=auc, se_auc=se, p=float(u.pvalue), ci_lo=ci_lo, ci_hi=ci_hi)


def dl_meta(est, se):
    est = np.asarray(est, float)
    se = np.asarray(se, float)
    k = len(est)
    w = 1.0 / (se**2)
    fe = (w * est).sum() / w.sum()
    Q = (w * (est - fe)**2).sum()
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w**2).sum() / w.sum())) if k > 1 else 0.0
    ws = 1.0 / (se**2 + tau2)
    re = float((ws * est).sum() / ws.sum())
    sre = float(math.sqrt(1.0 / ws.sum()))
    p = float(2.0 * stats.norm.sf(abs(re / sre)))
    I2 = float(max(0.0, (Q - (k - 1)) / Q) * 100.0) if (k > 1 and Q > 0) else 0.0
    return re, sre, p, I2, float(tau2), k


def compute_part3_figure7():
    print("\n" + "=" * 80)
    print("PART 3: FIGURE 7 — 66 HAEMODIALYSIS PATIENTS (276 ACUTE DRAWS) + 4 EXTERNAL STUDIES")
    print("=" * 80)

    meta_v = pd.read_csv(D28 / "eskd_sample_metadata_verified.csv")
    fine_pct = pd.read_csv(ESKD_RES / "sample_level_fine_label_pct.csv")
    scan29 = pd.read_csv(D29 / "fig7_scan_43fine_true11.csv")
    ct43 = list(dict.fromkeys(scan29["feature"].tolist()))
    assert len(ct43) == 43

    # Exact 276 eligible acute samples across 66 non-CLL haemodialysis patients
    sub276 = meta_v[
        (meta_v["condition"] == "positive")
        & (meta_v["in_primary_scope"] == True)
        & (meta_v["qc_pass_min200cells"] == True)
        & (~meta_v["patient_id"].isin(["C116", "C139"]))
    ].copy()
    assert len(sub276) == 276
    assert sub276["patient_id"].nunique() == 66
    assert (sub276["source"] == "Real").sum() == 138
    assert (sub276["source"] == "PBISC").sum() == 138
    assert sub276["temp_severity_ord"].isna().sum() == 0
    assert sub276["days_from_first_positive_swab"].isna().sum() == 0

    df276 = sub276.merge(fine_pct[["source", "sample_id"] + ct43], on=["source", "sample_id"], how="left")
    assert df276[ct43].isna().sum().sum() == 0

    # Classify PBISC samples into PBISC_new29 (79 draws) vs Shared27_PBISC (59 draws)
    real_pts = set(df276[df276["source"] == "Real"]["patient_id"].unique())
    pb_pts = set(df276[df276["source"] == "PBISC"]["patient_id"].unique())
    shared27_pts = real_pts & pb_pts
    new29_pts = pb_pts - real_pts
    real10_pts = real_pts - pb_pts
    assert len(real_pts) == 37 and len(pb_pts) == 56 and len(shared27_pts) == 27 and len(new29_pts) == 29 and len(real10_pts) == 10

    df276["sample_stratum"] = np.where(
        df276["source"] == "Real",
        np.where(df276["patient_id"].isin(shared27_pts), "Shared-27 Real", "Real-only (10 pts)"),
        np.where(df276["patient_id"].isin(shared27_pts), "Shared-27 PBISC", "PBISC-only (29 new pts)"),
    )
    df276["equal_patient_weight"] = 1.0 / df276.groupby("patient_id")["sample_id"].transform("count")
    assert abs(df276["equal_patient_weight"].sum() - 66.0) < 1e-9

    # Mark unique last eligible acute sample per patient (0 ties verified!)
    df276 = df276.sort_values(["patient_id", "days_from_first_positive_swab"]).reset_index(drop=True)
    max_day_by_pt = df276.groupby("patient_id")["days_from_first_positive_swab"].transform("max")
    df276["is_last_sample"] = df276["days_from_first_positive_swab"] == max_day_by_pt
    assert int(df276["is_last_sample"].sum()) == 66
    df66_last = df276[df276["is_last_sample"]].copy().reset_index(drop=True)
    assert int((df66_last["fatal_pub11"] == 1).sum()) == 11
    assert int((df66_last["fatal_pub11"] == 0).sum()) == 55

    df276.to_csv(DATA_OUT / "fig7_eskd_276_acute_samples_metadata_and_43pct.csv", index=False)
    df66_last.to_csv(DATA_OUT / "fig7_eskd_66_last_samples.csv", index=False)

    print(f"  Verified 276 acute samples: Real={int((df276['source']=='Real').sum())} (37 pts), PBISC={int((df276['source']=='PBISC').sum())} (29 new pts / 79 draws + 27 shared pts / 59 draws)")
    print(f"  Verified 66 unique last samples: 55 survived, 11 died | Last-sample source: {df66_last['source'].value_counts().to_dict()}")

    # Within-source rank-INT version for sensitivity check
    def r_int(s):
        r = stats.rankdata(s, method="average")
        return stats.norm.ppf((r - 0.5) / len(s))

    df276_int = df276.copy()
    for c in ct43:
        df276_int[c] = df276_int.groupby("source")[c].transform(r_int)
    df66_last_int = df276_int[df276_int["is_last_sample"]].copy().reset_index(drop=True)

    # Evaluate across strata: Real (138 draws / 37 pts), PBISC all (138 draws / 56 pts), PBISC new29 (79 draws / 29 pts), Combined (276 draws / 66 pts), Combined_INT
    strata_specs = [
        ("Real_37pts_138draws", df276[df276["source"] == "Real"], df276[df276["source"] == "Real"].sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last(), 1500),
        ("PBISC_56pts_138draws", df276[df276["source"] == "PBISC"], df276[df276["source"] == "PBISC"].sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last(), 1500),
        ("PBISC_new29pts_79draws", df276[df276["patient_id"].isin(new29_pts)], df276[df276["patient_id"].isin(new29_pts)].sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last(), 1500),
        ("Combined_66pts_276draws_rawPct", df276, df66_last, 2500),
        ("Combined_66pts_276draws_withinSourceINT", df276_int, df66_last_int, 1500),
    ]

    eskd_rows = []
    for st_name, d_sev, d_last, nb in strata_specs:
        tmp = []
        for fi, c in enumerate(ct43):
            sv = eval_equal_weight_severity(d_sev, c, "temp_severity_ord", "patient_id", n_boot=nb, seed=20261200 + fi)
            mt = eval_last_draw_mortality(d_last, c, "fatal_pub11", n_boot=nb, seed=20261500 + fi)
            tmp.append({
                "stratum": st_name, "feature": c,
                "n_samples_sev": sv["n_samples"], "n_patients_sev": sv["n_patients"],
                "rho_w": sv["rho_w"], "se_rho_clu": sv["se_clu"], "ci_lo_rho": sv["ci_lo"], "ci_hi_rho": sv["ci_hi"], "p_sev": sv["p"],
                "n_surv": mt["n_surv"], "n_dead": mt["n_dead"],
                "auc_death": mt["auc"], "se_auc": mt["se_auc"], "ci_lo_auc": mt["ci_lo"], "ci_hi_auc": mt["ci_hi"], "p_death": mt["p"],
            })
        tdf = pd.DataFrame(tmp)
        tdf["q_sev_bh43"] = bh_adjust(tdf["p_sev"])
        tdf["q_death_bh43"] = bh_adjust(tdf["p_death"])
        eskd_rows.append(tdf)
        print(f"  [{st_name}] BH43<0.05 counts -> Severity: {int((tdf['q_sev_bh43']<0.05).sum())}/43 | Last-draw Mortality: {int((tdf['q_death_bh43']<0.05).sum())}/43")

    df_eskd_res = pd.concat(eskd_rows, ignore_index=True)
    df_eskd_res.to_csv(DATA_OUT / "fig7_eskd_43celltypes_sameday_sev_lastdraw_death.csv", index=False)

    # Sensitivity checks on the exact 276-sample / 66-last-draw framework:
    # (1) Raw cell counts (pct * n_cells / 100), (2) Real hospital CBC absolute counts (cells/L = pct/100 * (temp_lymph + temp_mono)),
    # (3) CLL-included (68 patients, 283 acute draws)
    sub283_cll = meta_v[
        (meta_v["condition"] == "positive")
        & (meta_v["in_primary_scope"] == True)
        & (meta_v["qc_pass_min200cells"] == True)
    ].merge(fine_pct[["source", "sample_id"] + ct43], on=["source", "sample_id"], how="left")
    last68_cll = sub283_cll.sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last()

    df276_cnt = df276.copy()
    for c in ct43:
        df276_cnt[c] = df276_cnt[c] * df276_cnt["n_cells"] / 100.0
    last66_cnt = df276_cnt[df276_cnt["is_last_sample"]].copy()

    df_cbc = df276[(df276["source"] == "Real") & df276["temp_lymph"].notna() & df276["temp_mono"].notna()].copy()
    pbmc_abs = df_cbc["temp_lymph"].to_numpy(float) + df_cbc["temp_mono"].to_numpy(float)
    for c in ct43:
        df_cbc[c] = (df_cbc[c].to_numpy(float) / 100.0) * pbmc_abs
    last_cbc = df_cbc.sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last()

    sens_rows = []
    for c in ct43:
        sv_pct = df_eskd_res[(df_eskd_res["stratum"] == "Combined_66pts_276draws_rawPct") & (df_eskd_res["feature"] == c)].iloc[0]
        sv_int = df_eskd_res[(df_eskd_res["stratum"] == "Combined_66pts_276draws_withinSourceINT") & (df_eskd_res["feature"] == c)].iloc[0]
        sv_cnt = eval_equal_weight_severity(df276_cnt, c, "temp_severity_ord", "patient_id", n_boot=0)
        mt_cnt = eval_last_draw_mortality(last66_cnt, c, "fatal_pub11", n_boot=0)
        sv_cbc = eval_equal_weight_severity(df_cbc, c, "temp_severity_ord", "patient_id", n_boot=0)
        mt_cbc = eval_last_draw_mortality(last_cbc, c, "fatal_pub11", n_boot=0)
        sv_cll = eval_equal_weight_severity(sub283_cll, c, "temp_severity_ord", "patient_id", n_boot=0)
        mt_cll = eval_last_draw_mortality(last68_cll, c, "fatal_pub11", n_boot=0)
        sens_rows.append({
            "feature": c,
            "rho_rawPct_276": sv_pct["rho_w"], "q_sev_rawPct_276": sv_pct["q_sev_bh43"],
            "rho_INT_276": sv_int["rho_w"], "q_sev_INT_276": sv_int["q_sev_bh43"],
            "rho_cellCount_276": sv_cnt["rho_w"], "p_sev_cellCount_276": sv_cnt["p"],
            "rho_CBCabs_134": sv_cbc["rho_w"], "p_sev_CBCabs_134": sv_cbc["p"],
            "rho_withCLL_283": sv_cll["rho_w"], "p_sev_withCLL_283": sv_cll["p"],
            "auc_rawPct_last66": sv_pct["auc_death"], "q_death_rawPct_last66": sv_pct["q_death_bh43"],
            "auc_INT_last66": sv_int["auc_death"], "q_death_INT_last66": sv_int["q_death_bh43"],
            "auc_cellCount_last66": mt_cnt["auc"], "p_death_cellCount_last66": mt_cnt["p"],
            "auc_CBCabs_last37": mt_cbc["auc"], "p_death_CBCabs_last37": mt_cbc["p"],
            "auc_withCLL_last68": mt_cll["auc"], "p_death_withCLL_last68": mt_cll["p"],
        })
    df_sens = pd.DataFrame(sens_rows)
    df_sens["q_sev_cellCount_276"] = bh_adjust(df_sens["p_sev_cellCount_276"])
    df_sens["q_death_cellCount_last66"] = bh_adjust(df_sens["p_death_cellCount_last66"])
    df_sens["q_sev_withCLL_283"] = bh_adjust(df_sens["p_sev_withCLL_283"])
    df_sens["q_death_withCLL_last68"] = bh_adjust(df_sens["p_death_withCLL_last68"])
    df_sens.to_csv(DATA_OUT / "supp8_sensitivity_count_cbc_cll_276framework.csv", index=False)
    print(f"  Sensitivity correlation (276 rawPct vs cellCount): severity r={np.corrcoef(df_sens['rho_rawPct_276'], df_sens['rho_cellCount_276'])[0,1]:.4f}, mortality AUC r={np.corrcoef(df_sens['auc_rawPct_last66'], df_sens['auc_cellCount_last66'])[0,1]:.4f}")
    print(f"  Sensitivity correlation (276 rawPct vs Real CBC cells/L n={len(df_cbc)}): severity r={np.corrcoef(df_sens['rho_rawPct_276'], df_sens['rho_CBCabs_134'])[0,1]:.4f}")

    # Copy Official R 4.3.2 deconvolution benchmark CSVs from Folder 29 (refreshed 11/20)
    shutil.copy2(D29 / "fig8_deconv_12methods_officialR_true11.csv", DATA_OUT / "supp8_officialR_deconv_12methods.csv")
    shutil.copy2(D29 / "fig8_new29_hspc_12methods_officialR_CI.csv", DATA_OUT / "supp8_officialR_new29_hspc_auc.csv")
    shutil.copy2(PREV_SUB / "audit_focus3_20261121" / "A1_sepsis_255_samples_matched_metadata.csv", DATA_OUT / "supp6_sepsis_255_samples_metadata.csv")
    shutil.copy2(D29 / "table_s3_all_107_donors_overview.csv", DATA_OUT / "supp6_eskd_107_donors_overview.csv")

    # 3B. External 4 Studies (300 acute samples across 231 donors) — EXACT SAME CRITERION + Random-Effects Meta-Analysis (k >= 2)
    ext_src = pathlib.Path(r"C:\Users\scott\.gemini\antigravity\brain\add4945a-39a0-4f8d-acd4-3287c21c6121\scratch\external4_acute_samples_300.csv")
    ext_dst = DATA_OUT / "fig7_external4_300_acute_samples.csv"
    shutil.copy2(ext_src, ext_dst)
    ext = pd.read_csv(ext_dst)

    # Exact 1:1 mapped subsets from the 43 ESKD cell types (+ B_ASC_total documented as composite pan-ASC sensitivity)
    ext_feats_43 = [c[4:] for c in ext.columns if c.startswith("pct_") and c[4:] in set(ct43)]
    all_ext_feats = ext_feats_43 + ["B_ASC_total"]

    ext_coh_rows = []
    for coh in ["Stephenson 2021", "COMBAT 2022", "van der Wijst 2021", "Liu 2021"]:
        sub = ext[ext["cohort"] == coh].copy()
        sub_last = sub.sort_values(["donor", "time_ord"]).groupby("donor", as_index=False).last()
        for f in all_ext_feats:
            col = f"pct_{f}"
            if not sub[col].notna().any():
                continue
            sv = eval_equal_weight_severity(sub, col, "who_sameday", "donor", n_boot=1000, seed=20262000)
            mt = eval_last_draw_mortality(sub_last, col, "death", n_boot=1000, seed=20263000)
            z_f = float(np.arctanh(np.clip(sv["rho_w"], -0.999, 0.999)))
            se_zf = math.sqrt(1.06 / max(sv["n_patients"] - 3, 2))
            ext_coh_rows.append({
                "cohort": coh, "feature": f, "is_exact_43_subset": f in set(ct43),
                "n_samples_sev": sv["n_samples"], "n_donors_sev": sv["n_patients"],
                "rho_w": sv["rho_w"], "se_rho_clu": sv["se_clu"], "ci_lo_rho": sv["ci_lo"], "ci_hi_rho": sv["ci_hi"], "p_sev": sv["p"],
                "z_fisher": z_f, "se_z_fisher": se_zf,
                "n_surv": mt["n_surv"], "n_dead": mt["n_dead"],
                "auc_death": mt["auc"], "se_auc": mt["se_auc"], "ci_lo_auc": mt["ci_lo"], "ci_hi_auc": mt["ci_hi"], "p_death": mt["p"],
            })
    df_ext_coh = pd.DataFrame(ext_coh_rows)
    df_ext_coh.to_csv(DATA_OUT / "fig7_external4_per_cohort_results.csv", index=False)

    # Pool via DerSimonian-Laird random-effects whenever k >= 2 studies have matching annotations!
    eskd_comb_map = df_eskd_res[df_eskd_res["stratum"] == "Combined_66pts_276draws_rawPct"].set_index("feature")
    meta_rows = []
    for f, g in df_ext_coh.groupby("feature", sort=False):
        if len(g) < 2:
            continue
        zr, szr, pz, I2_s, tau2_s, k_s = dl_meta(g["z_fisher"], g["se_z_fisher"])
        ar, sar, pa, I2_m, tau2_m, k_m = dl_meta(g["auc_death"] - 0.5, g["se_auc"])
        auc_re = ar + 0.5
        e_rho = float(eskd_comb_map.loc[f, "rho_w"]) if f in eskd_comb_map.index else np.nan
        e_q_s = float(eskd_comb_map.loc[f, "q_sev_bh43"]) if f in eskd_comb_map.index else np.nan
        e_auc = float(eskd_comb_map.loc[f, "auc_death"]) if f in eskd_comb_map.index else np.nan
        e_q_m = float(eskd_comb_map.loc[f, "q_death_bh43"]) if f in eskd_comb_map.index else np.nan
        meta_rows.append({
            "feature": f,
            "is_exact_43_subset": bool(g["is_exact_43_subset"].iloc[0]),
            "k_studies": k_s,
            "cohorts": "; ".join(g["cohort"].tolist()),
            "n_samples_ext": int(g["n_samples_sev"].sum()),
            "n_donors_ext": int(g["n_donors_sev"].sum()),
            "rho_re": float(np.tanh(zr)),
            "ci_lo_rho_re": float(np.tanh(zr - 1.96 * szr)),
            "ci_hi_rho_re": float(np.tanh(zr + 1.96 * szr)),
            "p_sev_re": pz,
            "I2_sev": I2_s,
            "n_surv_ext": int(g["n_surv"].sum()),
            "n_dead_ext": int(g["n_dead"].sum()),
            "auc_re": auc_re,
            "ci_lo_auc_re": max(0.0, auc_re - 1.96 * sar),
            "ci_hi_auc_re": min(1.0, auc_re + 1.96 * sar),
            "p_death_re": pa,
            "I2_death": I2_m,
            "eskd_rho_w": e_rho,
            "eskd_q_sev_bh43": e_q_s,
            "eskd_auc_last66": e_auc,
            "eskd_q_death_bh43": e_q_m,
            "same_sign_sev": bool(np.sign(np.tanh(zr)) == np.sign(e_rho)) if np.isfinite(e_rho) else np.nan,
            "same_sign_death": bool(np.sign(ar) == np.sign(e_auc - 0.5)) if np.isfinite(e_auc) else np.nan,
        })
    df_ext_meta = pd.DataFrame(meta_rows)
    # Apply BH adjustment across the 23 exact 1:1 mapped subsets (and separately report B_ASC_total composite)
    m_exact = df_ext_meta["is_exact_43_subset"].to_numpy(bool)
    df_ext_meta["q_sev_re"] = np.nan
    df_ext_meta["q_death_re"] = np.nan
    df_ext_meta.loc[m_exact, "q_sev_re"] = bh_adjust(df_ext_meta.loc[m_exact, "p_sev_re"])
    df_ext_meta.loc[m_exact, "q_death_re"] = bh_adjust(df_ext_meta.loc[m_exact, "p_death_re"])
    # For B_ASC_total composite row, compute BH if included in 24-test family
    df_ext_meta["q_sev_re_all24"] = bh_adjust(df_ext_meta["p_sev_re"])
    df_ext_meta["q_death_re_all24"] = bh_adjust(df_ext_meta["p_death_re"])
    df_ext_meta.to_csv(DATA_OUT / "fig7_external4_random_effects_meta.csv", index=False)

    n_exact_k2 = int(m_exact.sum())
    print(f"\n  External 4-study RE Meta-Analysis (k>=2): {n_exact_k2} exact 1:1 mapped subsets (+ B_ASC_total composite)")
    print(f"    Severity q_sev_re < 0.05 (among {n_exact_k2} exact 1:1): {int((df_ext_meta.loc[m_exact, 'q_sev_re'] < 0.05).sum())}")
    print(f"    Mortality q_death_re < 0.05 (among {n_exact_k2} exact 1:1): {int((df_ext_meta.loc[m_exact, 'q_death_re'] < 0.05).sum())} | nominal P < 0.05: {int((df_ext_meta.loc[m_exact, 'p_death_re'] < 0.05).sum())}")

    # 3C. Single-cell T-cell ISG & HSPC functional modules under the unified framework (Supplementary Fig 10)
    shutil.copy2(D29 / "fig8_tcell_isg_13subsets_true11.csv", DATA_OUT / "supp10_tcell_isg_13subsets_summary.csv")
    shutil.copy2(D29 / "fig8_hspc_modules_16unique_true11.csv", DATA_OUT / "supp10_hspc_16modules_summary.csv")
    shutil.copy2(D29 / "suppS5_hspc_20gene_true11.csv", DATA_OUT / "supp10_hspc_20genes_summary.csv")
    shutil.copy2(D29 / "fig8_cd8ctl_isg_longitudinal_draws.csv", DATA_OUT / "supp10_cd8ctl_isg_draws.csv")

    # Evaluate CD8 CTL 7-ISG on equal-patient-weight same-day WHO, first-draw mortality, AND last-draw mortality!
    cd8_draws = pd.read_csv(D29 / "fig8_cd8ctl_isg_longitudinal_draws.csv")
    cd8_draws = cd8_draws.merge(
        df276[["source", "sample_id", "patient_id", "days_from_first_positive_swab", "temp_severity_ord", "fatal_pub11", "T CD8 CTL"]],
        on=["source", "sample_id", "patient_id"],
        how="inner",
        suffixes=("_old", ""),
    )
    cd8_first = cd8_draws.sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).first()
    cd8_last = cd8_draws.sort_values(["patient_id", "days_from_first_positive_swab"]).groupby("patient_id", as_index=False).last()
    sv_isg = eval_equal_weight_severity(cd8_draws, "isg_z", "temp_severity_ord", "patient_id", n_boot=1500)
    mt_isg_first = eval_last_draw_mortality(cd8_first, "isg_z", "fatal_pub11", n_boot=1500)
    mt_isg_last = eval_last_draw_mortality(cd8_last, "isg_z", "fatal_pub11", n_boot=1500)
    pd.DataFrame([
        {"analysis": "Same-day WHO (equal patient weight)", "metric": "rho_w", "estimate": sv_isg["rho_w"], "ci_lo": sv_isg["ci_lo"], "ci_hi": sv_isg["ci_hi"], "p": sv_isg["p"], "n_samples": sv_isg["n_samples"], "n_patients": sv_isg["n_patients"]},
        {"analysis": "First-draw mortality (early IFN)", "metric": "AUC", "estimate": mt_isg_first["auc"], "ci_lo": mt_isg_first["ci_lo"], "ci_hi": mt_isg_first["ci_hi"], "p": mt_isg_first["p"], "n_samples": mt_isg_first["n_surv"] + mt_isg_first["n_dead"], "n_patients": mt_isg_first["n_surv"] + mt_isg_first["n_dead"]},
        {"analysis": "Last-draw mortality (late IFN collapse)", "metric": "AUC", "estimate": mt_isg_last["auc"], "ci_lo": mt_isg_last["ci_lo"], "ci_hi": mt_isg_last["ci_hi"], "p": mt_isg_last["p"], "n_samples": mt_isg_last["n_surv"] + mt_isg_last["n_dead"], "n_patients": mt_isg_last["n_surv"] + mt_isg_last["n_dead"]},
    ]).to_csv(DATA_OUT / "supp10_cd8ctl_isg_sameday_first_last.csv", index=False)
    print(f"  T CD8 CTL ISG: Same-day WHO rho_w={sv_isg['rho_w']:+.3f} (P={sv_isg['p']:.3f}) | First-draw AUC={mt_isg_first['auc']:.3f} (P={mt_isg_first['p']:.4f}) | Last-draw AUC={mt_isg_last['auc']:.3f} (P={mt_isg_last['p']:.3f})")


if __name__ == "__main__":
    compute_part1_figure5()
    compute_part2_figure6()
    compute_part3_figure7()
    print("\nALL PART 1-3 COMPUTATIONS & ASSERTIONS COMPLETED SUCCESSFULLY.")
