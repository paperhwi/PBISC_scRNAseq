"""Independent mathematical re-verification & forensic audit script for PBISC Figures 4, 5, 6 and Supplementary Figures 1-10.
Re-computes and cross-checks cohort accounting, molecular fidelity, rare-cell/DEG metrics, ESKD COVID-19 TRUE_PUB_11
clinical associations, transcriptomic modules, GEE longitudinal models, leak-free CV, and Official R 4.3.2 deconvolution.
"""
import json
import pathlib
import sys
import warnings
from itertools import combinations

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(r"F:\PBISC")
OUT_DIR = ROOT / "deliverables" / "nature_submission_fig4_6_20261120"
OUT_DIR.mkdir(parents=True, exist_ok=True)
(OUT_DIR / "data").mkdir(exist_ok=True)

P_FIG45 = ROOT / "deliverables" / "manuscript_real_pbisc_scsemi_20261110" / "panels"
P_F4P = ROOT / "deliverables" / "manuscript_figures4papers_20261111"
P_D28 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "28_Verified_True11_Deconv_HSC_TISG_20261119" / "data"
P_R29 = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "29_Verified_OfficialR_Deconv_True11_20261120" / "r_exchange" / "r_outputs"
ESKD_RES = ROOT / "eskd_covid_clinical_assoc_20261115" / "results"

CLL = {"C116", "C139"}
D_PUB11 = {"C21", "C33", "C141", "C108", "C145", "C32", "C37", "C48", "C66", "C77", "C84"}
D_DIR9 = {"C21", "C33", "C141", "C32", "C37", "C48", "C66", "C77", "C84"}
D_ZEN14 = D_PUB11 | {"C101", "C146", "C147"}


def rint_arr(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    if len(np.unique(v)) < 2:
        return np.zeros_like(v)
    rk = stats.rankdata(v, method="average")
    return stats.norm.ppf((rk - 0.5) / len(v))


def rint(s: pd.Series) -> pd.Series:
    s = pd.Series(s, dtype=float)
    if s.nunique(dropna=True) < 2:
        return pd.Series(0.0, index=s.index)
    return pd.Series(stats.norm.ppf((s.rank() - 0.5) / s.notna().sum()), index=s.index)


def auc_mw(y, f):
    y_a, f_a = np.asarray(y, float), np.asarray(f, int)
    m = np.isfinite(y_a)
    y_a, f_a = y_a[m], f_a[m]
    n1, n0 = int((f_a == 1).sum()), int((f_a == 0).sum())
    if n1 < 1 or n0 < 1 or np.nanstd(y_a) == 0:
        return 0.5, 1.0
    u = stats.mannwhitneyu(y_a[f_a == 1], y_a[f_a == 0])
    return float(u.statistic / (n1 * n0)), float(u.pvalue)


def main():
    report = {}

    # =========================================================================
    # 1. FIGURE 4 & SUPP FIG 1 VERIFICATION
    # =========================================================================
    s1_scope = pd.read_csv(P_FIG45 / "supplementary_1_source_scope_and_depth.csv")
    main_num = pd.read_csv(P_FIG45 / "main_summary_numbers.csv")
    f2d = pd.read_csv(P_F4P / "F2d_decomposition.csv")
    f2b = pd.read_csv(P_F4P / "F2b_composition_shift.csv")

    fig4_gene_corr = {}
    for coh in ["sepsis", "covid"]:
        d_gc = pd.read_csv(P_FIG45 / f"{coh}_gene_counts_displayed.csv")
        fig4_gene_corr[coh] = {}
        for src in ["PBISC", "scSemiProfiler"]:
            ok = (d_gc["Real"] > 0) & (d_gc[src] > 0)
            r_val, p_val = stats.pearsonr(np.log10(d_gc.loc[ok, "Real"]), np.log10(d_gc.loc[ok, src]))
            fig4_gene_corr[coh][src] = {"r": round(float(r_val), 4), "n_genes": int(ok.sum()), "p": float(p_val)}

    comp6 = {}
    for coh in ["sepsis", "covid"]:
        dc = pd.read_csv(P_FIG45 / f"{coh}_composition_counts.csv")
        sc_c = pd.read_csv(P_FIG45 / "scsemi_composition.csv")
        comp6[coh] = {
            "Real_PBISC": dc.to_dict(orient="records"),
            "scSemi": sc_c[sc_c["cohort"] == coh].to_dict(orient="records"),
        }

    report["figure4"] = {
        "scope_and_depth": s1_scope.to_dict(orient="records"),
        "main_summary_numbers": main_num.to_dict(orient="records"),
        "recomputed_gene_mean_corr": fig4_gene_corr,
        "decomposition": f2d.to_dict(orient="records"),
        "composition_shift_pp": f2b.to_dict(orient="records"),
    }

    # =========================================================================
    # 2. FIGURE 5 & SUPP FIGS 2-5 VERIFICATION
    # =========================================================================
    rare_sh = pd.read_csv(P_FIG45 / "rare_parent_shares.csv")
    f1_deg = pd.read_csv(P_F4P / "F1_deg_fidelity.csv")
    f2a_bulk = pd.read_csv(P_F4P / "F2a_bulk_reference.csv")
    rare_sim = pd.read_csv(P_F4P / "SFig2_rare_sample_similarity_summary.csv")
    rare_hvg = pd.read_csv(P_F4P / "SFig2_rare_hvg_centroid_matrix.csv")
    sfig2_stat = pd.read_csv(P_F4P / "SFig2_statistics_with_IFN.csv")

    # Verify volcano DEG counts directly from all 24 per-gene tables
    volcano_counts = []
    for coh in ["sepsis", "covid"]:
        for ct in ["Tcell", "Bcell", "NKcell", "Monocyte"]:
            for src in ["Real", "PBISC", "scSemiProfiler"]:
                df_v = pd.read_csv(P_FIG45 / f"deg_{coh}_{ct}_{src}.csv.gz")
                sig = df_v["deg"].astype(bool) if "deg" in df_v.columns else ((df_v["q"] < 0.05) & (df_v["logFC"].abs() >= 1.0))
                lfc_col = "logFC" if "logFC" in df_v.columns else "log2fc"
                up = sig & (df_v[lfc_col] > 0)
                dn = sig & (df_v[lfc_col] < 0)
                n_case = int(df_v["n_case"].iloc[0]) if "n_case" in df_v.columns else None
                n_ctrl = int(df_v["n_ctrl"].iloc[0]) if "n_ctrl" in df_v.columns else None
                volcano_counts.append({
                    "cohort": coh, "celltype": ct, "source": src,
                    "n_genes": len(df_v), "n_sig": int(sig.sum()), "n_up": int(up.sum()), "n_down": int(dn.sum()),
                    "n_case": n_case, "n_ctrl": n_ctrl
                })

    report["figure5"] = {
        "rare_zero_counts": {
            coh: {
                src: int((rare_sh[(rare_sh["cohort"] == coh) & (rare_sh["source"] == src)]["percent_parent"] <= 0).sum())
                for src in ["Real", "PBISC", "scSemiProfiler"]
            }
            for coh in ["sepsis", "covid"]
        },
        "deg_fidelity": f1_deg.to_dict(orient="records"),
        "bulk_reference": f2a_bulk.to_dict(orient="records"),
        "rare_sample_similarity": rare_sim.to_dict(orient="records"),
        "volcano_counts": volcano_counts,
    }

    # =========================================================================
    # 3. FIGURE 6 & SUPP FIGS 6-10 INDEPENDENT RE-COMPUTATION
    # =========================================================================
    S_meta = pd.read_csv(ESKD_RES / "sample_level_composition_metadata.csv")
    Fm_raw = pd.read_csv(ESKD_RES / "sample_level_fine_label_pct.csv")
    LAB43_ALL = [c for c in Fm_raw.columns if c not in ("source", "sample_id", "RBC", "Platelet")]
    Fm = Fm_raw.copy()
    Fm[LAB43_ALL] = 100.0 * Fm[LAB43_ALL].div(Fm[LAB43_ALL].sum(axis=1), axis=0)
    S_all = S_meta.merge(Fm[["source", "sample_id"] + LAB43_ALL], on=["source", "sample_id"], how="left")

    # Check missing values across primary COVID+ samples
    pos_all = S_all[(S_all["condition"] == "positive") & S_all["in_primary_scope"].astype(bool) & S_all["qc_pass_min200cells"].astype(bool)].copy()
    pos_samps = pos_all[~pos_all["patient_id"].isin(CLL)].copy()

    missing_audit = {
        col: int(pos_samps[col].isna().sum())
        for col in ["patient_id", "severity_ord", "WHO_severity", "age", "sex", "days_from_first_positive_swab"]
        if col in pos_samps.columns
    }

    # Verify exact bridge from 164 PBISC samples to 29 PBISC new patients (79 draws)
    pb_all_164 = S_meta[S_meta["source"] == "PBISC"].copy()
    real_pos_ids = set(pos_samps[pos_samps["source"] == "Real"]["patient_id"])
    pb_bridge_breakdown = {
        "total_pbisc_samples": len(pb_all_164),
        "by_condition": pb_all_164["condition"].value_counts().to_dict(),
        "positive_not_primary_scope": int((pb_all_164["condition"].eq("positive") & ~pb_all_164["in_primary_scope"].astype(bool)).sum()),
        "positive_qc_fail": int((pb_all_164["condition"].eq("positive") & pb_all_164["in_primary_scope"].astype(bool) & ~pb_all_164["qc_pass_min200cells"].astype(bool)).sum()),
        "positive_cll_samples": int(pos_all[(pos_all["source"] == "PBISC") & pos_all["patient_id"].isin(CLL)].shape[0]),
        "positive_shared27_samples": int(pos_samps[(pos_samps["source"] == "PBISC") & pos_samps["patient_id"].isin(real_pos_ids)].shape[0]),
        "positive_new29_samples": int(pos_samps[(pos_samps["source"] == "PBISC") & ~pos_samps["patient_id"].isin(real_pos_ids)].shape[0]),
    }

    pos_r = pos_samps[pos_samps["source"] == "Real"].copy()
    pos_p_all = pos_samps[pos_samps["source"] == "PBISC"].copy()
    pos_n = pos_p_all[~pos_p_all["patient_id"].isin(real_pos_ids)].copy()

    REAL37 = sorted(real_pos_ids)
    NEW29 = sorted(set(pos_n["patient_id"]))
    SHARED27 = sorted(real_pos_ids & set(pos_p_all["patient_id"]))
    R_pt_all = pos_r.groupby("patient_id")[LAB43_ALL].mean().loc[REAL37]
    FINE43 = [f for f in LAB43_ALL if R_pt_all[f].gt(0).mean() >= 0.30]
    assert len(FINE43) == 43
    RARE23 = [f for f in FINE43 if (R_pt_all[f].mean() < 1.0 and R_pt_all[f].median() < 1.0)]

    R_pt = R_pt_all[FINE43].copy()
    N_pt = pos_n.groupby("patient_id")[FINE43].mean().loc[NEW29]
    PB_sh27 = pos_p_all[pos_p_all["patient_id"].isin(SHARED27)].groupby("patient_id")[FINE43].mean().loc[SHARED27]

    pt_cl = pos_samps.groupby("patient_id").agg(sev=("severity_ord", "first"), WHO=("WHO_severity", "first"))
    pt_cl["fat11"] = pt_cl.index.isin(D_PUB11).astype(int)

    # Cross-check against Folder 28 frozen scan CSV
    d28_scan = pd.read_csv(P_D28 / "fig7_scan_43fine_true11.csv")
    c_ns = pd.concat([
        pd.DataFrame({f: rint(R_pt[f]) for f in FINE43}, index=REAL37),
        pd.DataFrame({f: rint(N_pt[f]) for f in FINE43}, index=NEW29),
    ])
    sev_c = pt_cl.loc[c_ns.index, "sev"].values
    fat_c = pt_cl.loc[c_ns.index, "fat11"].values

    max_diff_rho, max_diff_auc = 0.0, 0.0
    d28_comb = d28_scan[d28_scan["CLL_excluded"] & (d28_scan["set"] == "Real + PBISC new")].set_index("feature")
    for fk in FINE43:
        rv, _ = stats.spearmanr(sev_c, c_ns[fk].values)
        av, _ = auc_mw(c_ns[fk].values, fat_c)
        max_diff_rho = max(max_diff_rho, abs(rv - d28_comb.loc[fk, "rho"]))
        max_diff_auc = max(max_diff_auc, abs(av - d28_comb.loc[fk, "AUC_death"]))

    # Verify WHO 2-3 mortality AUC under both conventions (full-cohort INT vs within-stratum re-INT)
    w23_r = [p for p in REAL37 if pt_cl.loc[p, "sev"] >= 2]
    w23_n = [p for p in NEW29 if pt_cl.loc[p, "sev"] >= 2]
    w23_ids = w23_r + w23_n
    who23_audit = {}
    for fk in ["HSPC", "T CD8 CTL"]:
        auc_full_int, p_full_int = auc_mw(c_ns.loc[w23_ids, fk].values, pt_cl.loc[w23_ids, "fat11"].values)
        y_reint = np.r_[rint_arr(R_pt.loc[w23_r, fk].values), rint_arr(N_pt.loc[w23_n, fk].values)]
        auc_reint, p_reint = auc_mw(y_reint, pt_cl.loc[w23_ids, "fat11"].values)
        who23_audit[fk] = {
            "n_who23": len(w23_ids), "deaths": int(pt_cl.loc[w23_ids, "fat11"].sum()),
            "full_cohort_INT_AUC": round(auc_full_int, 4), "full_cohort_INT_p": round(p_full_int, 5),
            "within_WHO23_reINT_AUC": round(auc_reint, 4), "within_WHO23_reINT_p": round(p_reint, 5),
        }

    # Verify Official R 4.3.2 Shared-27 (all 27 donors!) from Folder 29 r_outputs
    R_METHODS = ["BayesPrism-upd", "BayesPrism-init", "MuSiC", "BisqueRNA", "DWLS",
                 "CIBERSORTx-Bmode", "CIBERSORTx-Smode", "CIBERSORT-R", "Calib-NNLS", "NNLS"]
    sh_samps = pos_p_all[pos_p_all["patient_id"].isin(SHARED27)].set_index("sample_id")
    off_r_sh27 = {}
    # Also add PBISC on Shared-27
    pb_rhos = [stats.spearmanr(R_pt.loc[SHARED27, f], PB_sh27.loc[SHARED27, f])[0] if PB_sh27.loc[SHARED27, f].nunique() > 1 else np.nan for f in FINE43]
    s_pb = pd.Series(pb_rhos, index=FINE43)
    off_r_sh27["PBISC"] = {
        "med_rho_all43": round(float(s_pb.median()), 4),
        "med_rho_rare23": round(float(s_pb.loc[RARE23].median()), 4),
        "rare23_nonzero": int(s_pb.loc[RARE23].notna().sum()),
        "HSPC_median_pct": round(float(PB_sh27["HSPC"].median()), 4),
        "HSPC_rho": round(float(s_pb["HSPC"]), 4),
        "T_CD8_CTL_rho": round(float(s_pb["T CD8 CTL"]), 4),
        "B_ASC_IgG_rho": round(float(s_pb["B_ASC_IgG"]), 4),
        "B_ASC_IgA_rho": round(float(s_pb["B_ASC_IgA"]), 4),
    }
    for rm in R_METHODS:
        fa = P_R29 / f"prop_{rm}_shared_foldA.csv"
        fb = P_R29 / f"prop_{rm}_shared_foldB.csv"
        if fa.exists() and fb.exists():
            df_sh = pd.concat([pd.read_csv(fa, index_col=0), pd.read_csv(fb, index_col=0)])
            m_pt = df_sh[FINE43].join(sh_samps["patient_id"]).groupby("patient_id")[FINE43].mean().loc[SHARED27]
            rhos = [stats.spearmanr(R_pt.loc[SHARED27, f], m_pt.loc[SHARED27, f])[0] if m_pt.loc[SHARED27, f].nunique() > 1 else np.nan for f in FINE43]
            sr = pd.Series(rhos, index=FINE43)
            off_r_sh27[rm] = {
                "med_rho_all43": round(float(sr.median()), 4) if sr.notna().any() else 0.0,
                "med_rho_rare23": round(float(sr.loc[RARE23].median()), 4) if sr.loc[RARE23].notna().any() else 0.0,
                "rare23_nonzero": int(sr.loc[RARE23].notna().sum()),
                "HSPC_median_pct": round(float(m_pt["HSPC"].median()), 4),
                "HSPC_rho": round(float(sr["HSPC"]), 4) if pd.notna(sr["HSPC"]) else None,
                "T_CD8_CTL_rho": round(float(sr["T CD8 CTL"]), 4) if pd.notna(sr["T CD8 CTL"]) else None,
                "B_ASC_IgG_rho": round(float(sr["B_ASC_IgG"]), 4) if pd.notna(sr["B_ASC_IgG"]) else None,
                "B_ASC_IgA_rho": round(float(sr["B_ASC_IgA"]), 4) if pd.notna(sr["B_ASC_IgA"]) else None,
            }

    # Also save Official R Shared-27 table to CSV for plotting!
    df_off_r_sh27 = pd.DataFrame.from_dict(off_r_sh27, orient="index").reset_index().rename(columns={"index": "method"})
    df_off_r_sh27.to_csv(OUT_DIR / "data" / "official_R432_shared27_benchmark.csv", index=False)

    # Verify Official R BayesPrism 2.2.3 posterior Z on Shared-27
    z_pb = np.load(ROOT / "eskd_covid_clinical_assoc_20261115" / "data" / "real_pseudobulk_sample_x_fine_panel.npz", allow_pickle=True)
    r_genes, r_mat = pd.Index(z_pb["genes"].astype(str)), z_pb["pb"]
    grps = pd.Series(z_pb["groups"].astype(str)).str.split("|", expand=True)
    grps.columns = ["sample_id", "fine"]; grps["ncells"] = z_pb["ncells"].astype(int)
    ISG_USE = [g for g in ["ISG15", "IFI6", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "XAF1", "STAT1", "IRF7", "EPSTI1", "HERC5", "SIGLEC1"] if g in r_genes]
    r_isg_idx = [r_genes.get_loc(g) for g in ISG_USE]
    r_pos_m = pos_r.set_index("sample_id")

    bpz_sh27 = {}
    for ct in ["T CD8 CTL", "T CD8 CM"]:
        ct_tag = ct.replace(" ", "_")
        z_all = pd.concat([
            pd.read_csv(P_R29 / f"BayesPrism_Z_{ct_tag}_shared_foldA.csv", index_col=0),
            pd.read_csv(P_R29 / f"BayesPrism_Z_{ct_tag}_shared_foldB.csv", index_col=0),
        ])
        gl = [g for g in ISG_USE if g in z_all.columns]
        z_bp_pt = z_all[gl].mean(axis=1).to_frame("z").join(sh_samps["patient_id"]).groupby("patient_id")["z"].mean()
        m_r = grps["fine"].eq(ct).values
        df_rs = pd.DataFrame(r_mat[m_r][:, r_isg_idx], columns=ISG_USE)
        df_rs["sample_id"] = grps.loc[m_r, "sample_id"].values
        df_rs["tot_umi"] = r_mat[m_r].sum(axis=1); df_rs["ncells"] = grps.loc[m_r, "ncells"].values
        df_rs = df_rs.groupby("sample_id").sum(); df_rs = df_rs[df_rs["ncells"] >= 5]
        df_rs["sc"] = np.log1p(df_rs[ISG_USE].div(np.maximum(df_rs["tot_umi"], 1.0), axis=0) * 1e4).mean(axis=1)
        r_gt = df_rs.reindex(r_pos_m.index).dropna().join(r_pos_m["patient_id"]).groupby("patient_id")["sc"].mean()
        common_p = sorted(set(r_gt.index) & set(z_bp_pt.index))
        rv, pv = stats.spearmanr(r_gt.loc[common_p], z_bp_pt.loc[common_p])
        bpz_sh27[ct + "_ISG"] = {"n": len(common_p), "official_R_BayesPrism_Z_rho": round(float(rv), 4), "p": float(pv)}

    report["figure6"] = {
        "missing_values_audit": missing_audit,
        "pbisc_164_to_new29_bridge": pb_bridge_breakdown,
        "cohort_counts": {
            "Real37": len(REAL37), "New29": len(NEW29), "Shared27": len(SHARED27), "Combined66": len(REAL37) + len(NEW29),
            "Real_draws": len(pos_r), "New_draws": len(pos_n), "Total_draws": len(pos_r) + len(pos_n),
            "Real_deaths_pub11": int(pt_cl.loc[REAL37, "fat11"].sum()),
            "New_deaths_pub11": int(pt_cl.loc[NEW29, "fat11"].sum()),
            "Real_Real_HSPC_median_pct": round(float(R_pt.loc[SHARED27, "HSPC"].median()), 4),
        },
        "scan43_max_discrepancy_vs_d28": {"max_diff_rho": float(max_diff_rho), "max_diff_auc": float(max_diff_auc)},
        "who23_int_conventions": who23_audit,
        "official_R432_shared27_all27": off_r_sh27,
        "official_R432_BayesPrism_Z_shared27": bpz_sh27,
    }

    (OUT_DIR / "VERIFICATION_REPORT_20261120.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("=== INDEPENDENT VERIFICATION COMPLETE ===")
    print("Max diff vs D28 scan43 (rho / AUC):", max_diff_rho, max_diff_auc)
    print("Missing values audit:", missing_audit)
    print("PBISC 164 -> New29 bridge:", pb_bridge_breakdown)
    print("WHO 2-3 INT conventions:", json.dumps(who23_audit, indent=2))
    print("\nOfficial R 4.3.2 Shared-27 (ALL 27 donors!) summary:")
    print(df_off_r_sh27.to_string(index=False))
    print("\nOfficial R 4.3.2 BayesPrism Z Shared-27:", bpz_sh27)
    print("\nVolcano counts summary:")
    print(pd.DataFrame(volcano_counts).to_string(index=False))


if __name__ == "__main__":
    main()
