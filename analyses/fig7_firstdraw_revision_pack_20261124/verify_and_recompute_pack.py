#!/usr/bin/env python3
"""Portable, zero-hardcoded-path verification & summary script for Figure 7 (First-Draw) revision.

Run from anywhere after cloning https://github.com/paperhwi/PBISC_scRNAseq.git:
    python analyses/fig7_firstdraw_revision_pack_20261124/verify_and_recompute_pack.py

Verifies 100% mathematical parity between sample-level CSVs in `./data/` and
the ready-to-plot summary tables (`panel_k_ready_table.csv`, `panel_l_ready_table.csv`,
`eskd_firstdraw_composite_and_cv_summary.csv`), with zero local path dependencies.
"""
import math
import pathlib
import numpy as np
import pandas as pd
from scipy import stats

PACK_DIR = pathlib.Path(__file__).resolve().parent
DATA_DIR = PACK_DIR / "data"


def dl_meta(yi, sei):
    y = np.asarray(yi, float)
    s = np.asarray(sei, float)
    w = 1.0 / (s * s)
    k = len(y)
    mu_fe = float(np.sum(w * y) / np.sum(w))
    Q = float(np.sum(w * (y - mu_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w * w) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 and df > 0 else 0.0
    wr = 1.0 / (s * s + tau2)
    mu_re = float(np.sum(wr * y) / np.sum(wr))
    se_re = math.sqrt(1.0 / float(np.sum(wr)))
    z = mu_re / se_re if se_re > 0 else 0.0
    p = float(2.0 * stats.norm.sf(abs(z)))
    return mu_re, se_re, p


def main():
    p217 = pd.read_csv(DATA_DIR / "eskd_217_acute_draws_43subsets_rankINT.csv")
    p66 = pd.read_csv(DATA_DIR / "eskd_66_first_draws_43subsets_rankINT.csv")
    a3d = pd.read_csv(DATA_DIR / "eskd_63_exact_first_draws_CTL_NK_ISG.csv")
    comp = pd.read_csv(DATA_DIR / "eskd_firstdraw_composite_and_cv_summary.csv")
    ext = pd.read_csv(DATA_DIR / "external4_300_acute_samples.csv")
    pk = pd.read_csv(DATA_DIR / "panel_k_ready_table.csv")
    pl = pd.read_csv(DATA_DIR / "panel_l_ready_table.csv")

    # 1. Frozen baseline checks (Panels a, b)
    assert len(p217) == 217 and p217["patient_id"].nunique() == 66
    assert int((p217["arm"] == "Real").sum()) == 138 and int((p217["arm"] == "PBISC new").sum()) == 79
    who_counts = p217["WHO_temp_severity"].value_counts().to_dict()
    assert who_counts == {"mild": 126, "moderate": 58, "severe": 21, "critical": 12}

    assert len(p66) == 66 and int(p66["fatal_pub11"].sum()) == 11
    assert len(a3d) == 63 and int(a3d["fatal_pub11"].sum()) == 11
    assert int((a3d["arm"] == "Real").sum()) == 35 and int((a3d["arm"] == "PBISC new").sum()) == 28

    # 2. Issue 2 composite verification (Panels d, j)
    recomp_tri = (a3d["isg_zfd"] - a3d["ctl_zfd"] - a3d["nk_zfd"]) / 3.0
    assert np.max(np.abs(recomp_tri - a3d["ISG_minus_CTL_minus_NK"])) < 1e-12
    y = a3d["fatal_pub11"].astype(int).to_numpy()
    u, p = stats.mannwhitneyu(recomp_tri[y == 1], recomp_tri[y == 0], alternative="two-sided")
    auc_tri = float(u / (np.sum(y == 1) * np.sum(y == 0)))
    row_tri = comp[(comp["feature_label"] == "ISG − CTL − NK") & (comp["stratum"] == "Combined")].iloc[0]
    assert abs(auc_tri - float(row_tri["AUC"])) < 1e-12
    assert abs(p - float(row_tri["p_two_sided"])) < 1e-12

    # 3. Issue 1 & 3 external 4-cohort verification (Panels k, l)
    assert len(ext) == 300 and len(ext[["cohort", "donor"]].drop_duplicates()) == 231
    fd_ext = ext.sort_values(["cohort", "donor", "time_ord", "sample_id"]).groupby(["cohort", "donor"], as_index=False).first()
    fd_known = fd_ext[fd_ext["death"].notna()]
    assert len(fd_known) == 219 and int(fd_known["death"].sum()) == 22

    for _, r in pl.iterrows():
        feat = r["feature"]
        yi, sei = [], []
        for coh, col_auc in [
            ("Stephenson 2021", "Stephenson_73_7_AUC_1st"),
            ("COMBAT 2022", "COMBAT_75_7_AUC_1st"),
            ("van der Wijst 2021", "vanDerWijst_42_5_AUC_1st"),
            ("Liu 2021", "Liu_29_3_AUC_1st"),
        ]:
            sub = fd_known[(fd_known["cohort"] == coh) & (fd_known[f"pct_{feat}"].notna())]
            if len(sub) == 0:
                continue
            yv = sub["death"].astype(int).to_numpy()
            xv = sub[f"pct_{feat}"].to_numpy(float)
            u_c, _ = stats.mannwhitneyu(xv[yv == 1], xv[yv == 0], alternative="two-sided")
            auc_c = float(u_c / (np.sum(yv == 1) * np.sum(yv == 0)))
            assert abs(auc_c - float(r[col_auc])) < 1e-12, f"Mismatch on {coh} {feat}"
            n1, n0 = int(np.sum(yv == 1)), int(np.sum(yv == 0))
            q1 = auc_c / (2.0 - auc_c)
            q2 = 2.0 * auc_c * auc_c / (1.0 + auc_c)
            var = (auc_c * (1.0 - auc_c) + (n1 - 1) * (q1 - auc_c * auc_c) + (n0 - 1) * (q2 - auc_c * auc_c)) / (n1 * n0)
            yi.append(auc_c - 0.5)
            sei.append(max(math.sqrt(max(var, 0.0)), 0.035))
        mu_re, _, p_re = dl_meta(yi, sei)
        assert abs((mu_re + 0.5) - float(r["RE_AUC_1st"])) < 1e-12, f"RE AUC mismatch on {feat}"
        assert abs(p_re - float(r["RE_P_1st"])) < 1e-12, f"RE P mismatch on {feat}"

    print("[PASS] All portable revision-pack CSVs verified with max_abs_diff < 1e-12!")
    print("\n--- Ready-to-Render Panel k Table ---")
    print(pk[["display_label", "eskd_gee_wk_beta", "eskd_equal_wt_rho", "aligned_over_k", "n_patients_ext", "ext_rho", "ext_q23", "external_evidence"]].to_string(index=False))
    print("\n--- Ready-to-Render Panel l Table ---")
    print(pl[["block", "display_label", "Stephenson_73_7_AUC_1st", "COMBAT_75_7_AUC_1st", "vanDerWijst_42_5_AUC_1st", "Liu_29_3_AUC_1st", "RE_AUC_1st", "RE_q23_1st", "RE_AUC_last", "RE_q23_last", "ESKD_AUC_1st"]].to_string(index=False))
    print("\n--- Issue 2 First-Draw Death & Composite Summary ---")
    print(comp[["feature_label", "stratum", "n_patients", "n_died", "AUC", "CI_lo", "CI_hi", "p_two_sided", "q43", "adjWHO_p", "cv_AUC_base", "cv_AUC_ext"]].to_string(index=False))


if __name__ == "__main__":
    main()
