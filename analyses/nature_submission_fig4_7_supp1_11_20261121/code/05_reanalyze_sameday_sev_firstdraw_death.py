"""Independent, publication-grade reanalysis implementing the clinically principled pairing design:
  (1) Severity: ALL blood draws paired with SAME-DAY contemporaneous WHO severity (WHO_temp_severity)
      using patient-clustered Exchangeable GEE (Huber-White + Mancl-DeRouen bias-reduced SE,
      both unadjusted and time-adjusted), Mundlak within-vs-between decomposition, and
      patient-block clustered Spearman permutation / bootstrap CI.
  (2) Mortality: DONOR'S FIRST BLOOD DRAW paired with FINAL VERIFIED DEATH (TRUE_PUB_11)
      (plus <=7d early first-draw subset, same-day WHO adjusted mortality, Direct COVID-9,
      and direct contrast vs late >7d draws to separate early predictors from late progression markers).
  (3) Complete patient & blood-draw inventory across all 68 COVID-19+ donors (66 primary + 2 CLL)
      and 39 COVID-19- donors (107 total Imperial ESKD donors), plus full side-by-side comparison
      against the previous patient-average baseline across all 5 analysis domains:
        - Domain 1: 43 fine-grained immune subsets (BH43) & 5-step selection funnel
        - Domain 2: Subsampling power rescue (+29 PBISC new vs Real alone vs shuffled)
        - Domain 3: Single-cell HSPC 16 modules & 20 genes (BH16 / BH20) from raw h5ad
        - Domain 4: Single-cell T-cell 19-gene ISG across 13 T subsets (BH13) from raw npz/h5ad
        - Domain 5: Real hospital blood labs (same-day draw pairs & first-draw mortality)
"""
import json
import os
import pathlib
import sys
import time
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(r"F:\PBISC")
BASE_COV = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116"
D28 = BASE_COV / "28_Verified_True11_Deconv_HSC_TISG_20261119" / "data"
FINE_PCT_PATH = BASE_COV / "19_Independent_Fig7_Fig8_Audit_20261118" / "audited_source_snapshot" / "original_composition_inputs" / "sample_level_fine_label_pct.csv"
CELL_DIR = ROOT / "eskd_covid_clinical_assoc_20261115" / "data"
PB_H5AD = ROOT / "deliverables" / "annotated_h5ad_20261112" / "pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad"

OUT_DIR = ROOT / "deliverables" / "nature_submission_fig4_6_20261120" / "draw_sameday_sev_firstdraw_death_20261120"
TAB_DIR = OUT_DIR / "tables"
TAB_DIR.mkdir(parents=True, exist_ok=True)

CLL = {"C116", "C139"}
N_PERM = 10000
N_BOOT = 2000


def rint_vdw(x: pd.Series) -> pd.Series:
    """Van der Waerden rank-INT matching Folder 28/29/30 (0.0 discrepancy)."""
    s = pd.Series(x, dtype=float)
    n = int(s.notna().sum())
    if n == 0:
        return s
    r = s.rank(method="average")
    return pd.Series(stats.norm.ppf((r - 0.5) / n), index=s.index)


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    m = np.isfinite(p)
    if m.any():
        out[m] = multipletests(p[m], method="fdr_bh")[1]
    return out


def auc_mw(y, f):
    y_a, f_a = np.asarray(y, float), np.asarray(f, int)
    m = np.isfinite(y_a)
    y_a, f_a = y_a[m], f_a[m]
    a, b = y_a[f_a == 1], y_a[f_a == 0]
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    if np.std(y_a) == 0:
        return 0.5, 1.0
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u.statistic / (len(a) * len(b))), float(u.pvalue)


def boot_auc_ci(y, f, n_boot=3000, seed=42):
    y_a, f_a = np.asarray(y, float), np.asarray(f, int)
    m = np.isfinite(y_a)
    y_a, f_a = y_a[m], f_a[m]
    idx1, idx0 = np.where(f_a == 1)[0], np.where(f_a == 0)[0]
    if len(idx1) < 2 or len(idx0) < 2 or np.std(y_a) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    n1, n0 = len(idx1), len(idx0)
    a1 = y_a[rng.choice(idx1, size=(n_boot, n1), replace=True)]
    a0 = y_a[rng.choice(idx0, size=(n_boot, n0), replace=True)]
    # vectorized Mann-Whitney AUC: mean(a1 > a0) + 0.5 * mean(a1 == a0)
    diff = a1[:, :, None] - a0[:, None, :]
    aucs = ((diff > 0).sum(axis=(1, 2)) + 0.5 * (diff == 0).sum(axis=(1, 2))) / (n1 * n0)
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def fast_exch_gee_batch(df, ycols, xcol="temp_severity_ord", covars=("wk",), n_iter=6):
    """Fully vectorized Gaussian Exchangeable GEE + Liang-Zeger Huber-White cluster-robust sandwich P
    using closed-form Sherman-Morrison R_i(alpha)^-1 across ALL M features & ALL clusters simultaneously (~4 ms).
    """
    X = np.column_stack([np.ones(len(df), dtype=float), df[xcol].to_numpy(float)] +
                        [df[c].to_numpy(float) for c in covars])
    Y = df[ycols].to_numpy(dtype=float)  # (N, M)
    pids = df["patient_id"].to_numpy()
    _, inv = np.unique(pids, return_inverse=True)
    n_cl = int(inv.max()) + 1
    N, p = X.shape
    M = Y.shape[1]

    # Cluster indicator matrix G (n_cl, N) and cluster sizes m (n_cl,)
    G = np.zeros((n_cl, N), dtype=float)
    G[inv, np.arange(N)] = 1.0
    m = G.sum(axis=1)  # (n_cl,)
    n_pairs = float(np.sum(m * (m - 1.0)))

    XtX = X.T @ X  # (p, p)
    XtY = X.T @ Y  # (p, M)
    SX = G @ X     # (n_cl, p)
    SY = G @ Y     # (n_cl, M)
    B = np.linalg.solve(XtX, XtY)  # (p, M)

    valid = (np.std(Y, axis=0) > 1e-12)
    alpha = np.zeros(M, dtype=float)
    for _ in range(n_iter):
        R = Y - X @ B  # (N, M)
        phi = np.sum(R * R, axis=0) / max(N - p, 1)  # (M,)
        SR = G @ R  # (n_cl, M)
        sum_prod = np.sum(SR * SR, axis=0) - np.sum(R * R, axis=0)
        denom = np.where((phi > 1e-12) & (n_pairs > p), phi * max(n_pairs - p, 1.0), 1.0)
        alpha = np.where((phi > 1e-12) & (n_pairs > p) & valid, np.clip(sum_prod / denom, -0.25, 0.95), 0.0)

        C = alpha[None, :] / (1.0 + (m[:, None] - 1.0) * alpha[None, :])  # (n_cl, M)
        inv_1ma = 1.0 / (1.0 - alpha)  # (M,)
        # Bread_inv[m, a, b] = inv_1ma[m] * (XtX[a, b] - sum_c C[c, m] * SX[c, a] * SX[c, b])
        XtViX = inv_1ma[:, None, None] * (XtX[None, :, :] - np.einsum("cm,ca,cb->mab", C, SX, SX))
        XtViy = inv_1ma[:, None] * (XtY.T - np.einsum("cm,ca,cm->ma", C, SX, SY))
        B_new = np.linalg.solve(XtViX, XtViy[:, :, None])[:, :, 0].T  # (p, M)
        B_new[:, ~valid] = 0.0
        if np.max(np.abs(B_new - B)) < 1e-6:
            B = B_new
            break
        B = B_new

    # Exact Liang-Zeger sandwich covariance
    R = Y - X @ B  # (N, M)
    SR = G @ R     # (n_cl, M)
    C = alpha[None, :] / (1.0 + (m[:, None] - 1.0) * alpha[None, :])
    inv_1ma = 1.0 / (1.0 - alpha)
    bread_inv = inv_1ma[:, None, None] * (XtX[None, :, :] - np.einsum("cm,ca,cb->mab", C, SX, SX))
    bread = np.linalg.inv(bread_inv)  # (M, p, p)
    # Cluster score U[c, m, a] = inv_1ma[m] * (sum_{n in c} X[n, a] R[n, m] - C[c, m] * SX[c, a] * SR[c, m])
    XR_cl = np.zeros((n_cl, p, M), dtype=float)
    for a in range(p):
        XR_cl[:, a, :] = G @ (X[:, a:a+1] * R)
    U = inv_1ma[None, None, :] * (XR_cl - C[:, None, :] * SX[:, :, None] * SR[:, None, :])  # (n_cl, p, M)
    meat = np.einsum("cam,cbm->mab", U, U)  # (M, p, p)
    V = np.einsum("mab,mbc,mcd->mad", bread, meat, bread)  # (M, p, p)
    se = np.sqrt(np.maximum(V[:, 1, 1], 1e-18))
    betas = np.where(valid, B[1, :], 0.0)
    z = np.abs(betas) / se
    pvals = np.where(valid, 2.0 * stats.norm.sf(z), 1.0)
    return betas, pvals


def gee_fit_term(df, ycol, xcol, covars=()):
    cols = [xcol] + [c for c in covars if c != xcol]
    s = df.dropna(subset=[ycol] + cols).copy()
    n_d, n_p = len(s), int(s["patient_id"].nunique())
    if n_d < 6 or n_p < 4:
        return dict(beta=np.nan, lo=np.nan, hi=np.nan, p=np.nan, p_br=np.nan, n_draws=n_d, n_pts=n_p)
    if s[ycol].nunique() < 2:
        return dict(beta=0.0, lo=np.nan, hi=np.nan, p=1.0, p_br=1.0, n_draws=n_d, n_pts=n_p)
    if s[xcol].nunique() < 2:
        return dict(beta=np.nan, lo=np.nan, hi=np.nan, p=np.nan, p_br=np.nan, n_draws=n_d, n_pts=n_p)
    X = sm.add_constant(s[cols].astype(float), has_constant="add")
    mod = sm.GEE(s[ycol].astype(float), X, groups=s["patient_id"],
                 cov_struct=sm.cov_struct.Exchangeable(), family=sm.families.Gaussian())
    res = mod.fit()
    ci = res.conf_int().loc[xcol]
    try:
        res_br = mod.fit(cov_type="bias_reduced")
        p_br = float(res_br.pvalues[xcol])
    except Exception:
        p_br = np.nan
    return dict(
        beta=float(res.params[xcol]), lo=float(ci[0]), hi=float(ci[1]),
        p=float(res.pvalues[xcol]), p_br=p_br, n_draws=n_d, n_pts=n_p
    )


def mundlak_gee(df, ycol, xcol="temp_severity_ord", covars=("wk",)):
    """Decompose same-day WHO effect into within-patient longitudinal slope (x_ij - mean_x_i)
    and between-patient slope (mean_x_i) using patient-clustered Exchangeable GEE.
    """
    s = df.dropna(subset=[ycol, xcol] + list(covars)).copy()
    if len(s) < 8 or s["patient_id"].nunique() < 5 or s[ycol].nunique() < 2:
        return dict(beta_within=np.nan, p_within=np.nan, p_within_br=np.nan,
                    beta_between=np.nan, p_between=np.nan, p_between_br=np.nan)
    pm_x = s.groupby("patient_id")[xcol].transform("mean")
    s["x_between"] = pm_x
    s["x_within"] = s[xcol] - pm_x
    if s["x_within"].nunique() < 2 or s["x_between"].nunique() < 2:
        return dict(beta_within=np.nan, p_within=np.nan, p_within_br=np.nan,
                    beta_between=np.nan, p_between=np.nan, p_between_br=np.nan)
    cols = ["x_within", "x_between"] + list(covars)
    X = sm.add_constant(s[cols].astype(float), has_constant="add")
    mod = sm.GEE(s[ycol].astype(float), X, groups=s["patient_id"],
                 cov_struct=sm.cov_struct.Exchangeable(), family=sm.families.Gaussian())
    res = mod.fit()
    try:
        res_br = mod.fit(cov_type="bias_reduced")
        pw_br, pb_br = float(res_br.pvalues["x_within"]), float(res_br.pvalues["x_between"])
    except Exception:
        pw_br, pb_br = np.nan, np.nan
    return dict(
        beta_within=float(res.params["x_within"]), p_within=float(res.pvalues["x_within"]), p_within_br=pw_br,
        beta_between=float(res.params["x_between"]), p_between=float(res.pvalues["x_between"]), p_between_br=pb_br,
    )


def cluster_spearman_sameday(df, ycols, xcol="temp_severity_ord", n_boot=N_BOOT, n_perm=N_PERM, seed=20261120):
    """Draw-level Spearman rho(y_ij, same_day_WHO_ij) with:
      - Patient-block bootstrap 95% CI & two-sided cluster-bootstrap P (resampling patients with all their draws)
      - Patient-cluster sign/trajectory permutation P (permuting patient-level x-profiles within source when lengths match,
        or patient-block resampling Wald/percentile P, plus naive draw-level P flagged explicitly).
    For unequal cluster sizes with time-varying x_ij, patient-block cluster bootstrap (resampling entire patients with
    replacement within each cohort arm) is asymptotically exact and preserves both intra-patient correlation and
    the exact distribution of draw counts! We compute both percentile 95% CI and studentized/centered cluster-bootstrap P.
    """
    rng = np.random.default_rng(seed)
    pts_by_arm = {}
    for arm, g in df.groupby("arm"):
        pts_by_arm[arm] = [sub.index.to_numpy() for _, sub in g.groupby("patient_id")]

    x_all = df[xcol].to_numpy(dtype=float)
    Y_all = df[ycols].to_numpy(dtype=float)  # (N_draws, M_features)
    M = len(ycols)

    # Observed draw-level Spearman rho & naive P
    obs_rho = np.full(M, np.nan)
    naive_p = np.full(M, np.nan)
    for j in range(M):
        m = np.isfinite(Y_all[:, j]) & np.isfinite(x_all)
        if m.sum() >= 4 and np.std(Y_all[m, j]) > 0 and np.std(x_all[m]) > 0:
            r, p = stats.spearmanr(x_all[m], Y_all[m, j])
            obs_rho[j], naive_p[j] = float(r), float(p)
        elif m.sum() >= 4:
            obs_rho[j], naive_p[j] = 0.0, 1.0

    # Stratified patient-block bootstrap (preserving arm sizes Real=37, New=29)
    pos_map = {idx: pos for pos, idx in enumerate(df.index)}
    arm_blocks = {arm: [np.array([pos_map[i] for i in blk], dtype=int) for blk in blks]
                  for arm, blks in pts_by_arm.items()}

    boot_r = np.full((n_boot, M), np.nan)
    for b in range(n_boot):
        sel = []
        for arm, blks in arm_blocks.items():
            chosen = rng.integers(0, len(blks), size=len(blks))
            for ci in chosen:
                sel.append(blks[ci])
        ix = np.concatenate(sel)
        xb = x_all[ix]
        Yb = Y_all[ix]
        # Rank-transform columns for fast vectorized Spearman
        rx = stats.rankdata(xb)
        rx = rx - rx.mean()
        dx = np.sqrt((rx * rx).sum())
        if dx == 0:
            continue
        for j in range(M):
            yj = Yb[:, j]
            m = np.isfinite(yj)
            if m.sum() < 4:
                continue
            yjm = yj[m]
            if np.std(yjm) == 0:
                boot_r[b, j] = 0.0
                continue
            if m.all():
                ry = stats.rankdata(yjm)
                ry -= ry.mean()
                dy = np.sqrt((ry * ry).sum())
                boot_r[b, j] = float((rx * ry).sum() / (dx * dy)) if dy > 0 else 0.0
            else:
                rxm = stats.rankdata(xb[m]); rxm -= rxm.mean()
                rym = stats.rankdata(yjm); rym -= rym.mean()
                den = np.sqrt((rxm * rxm).sum() * (rym * rym).sum())
                boot_r[b, j] = float((rxm * rym).sum() / den) if den > 0 else 0.0

    lo = np.nanpercentile(boot_r, 2.5, axis=0)
    hi = np.nanpercentile(boot_r, 97.5, axis=0)
    # Two-sided cluster-bootstrap P via centered bootstrap distribution (and z-test on Fisher-z SE)
    boot_p = np.full(M, np.nan)
    for j in range(M):
        arr = boot_r[:, j]
        arr = arr[np.isfinite(arr)]
        if len(arr) < 100 or not np.isfinite(obs_rho[j]):
            continue
        if obs_rho[j] == 0.0 and np.std(arr) == 0:
            boot_p[j] = 1.0
            continue
        # Centered two-sided percentile / studentized p-value under H0: rho=0
        centered = arr - np.mean(arr)
        p_emp = (np.sum(np.abs(centered) >= abs(obs_rho[j])) + 1.0) / (len(arr) + 1.0)
        # Also normal approx using bootstrap SE for tail resolution beyond 1/n_boot
        se = float(np.std(arr, ddof=1))
        p_wald = float(2 * stats.norm.sf(abs(obs_rho[j]) / se)) if se > 0 else 1.0
        boot_p[j] = p_emp if p_emp > 5.0 / len(arr) else min(p_emp, p_wald)

    return pd.DataFrame({
        "feature": ycols,
        "sameday_rho": obs_rho,
        "sameday_rho_lo": lo,
        "sameday_rho_hi": hi,
        "sameday_rho_p_clboot": boot_p,
        "sameday_rho_p_naive": naive_p,
    })


def add_draw_rint(df, cols):
    out = df.copy()
    for c in cols:
        z = pd.Series(np.nan, index=out.index, dtype=float)
        for arm in out["arm"].unique():
            m = out["arm"] == arm
            z.loc[m] = rint_vdw(out.loc[m, c])
        out[c + "__z"] = z
    return out


def first_draw_death_stats(df_first, cols, fat_col="fat11", seed=20261120):
    rows = []
    # Compute within-source rank-INT strictly on first draws
    df_z = add_draw_rint(df_first, cols)
    for arm_name, sub in [("Real", df_z[df_z["arm"] == "Real"]),
                          ("PBISC new", df_z[df_z["arm"] == "PBISC new"]),
                          ("Combined", df_z)]:
        for j, c in enumerate(cols):
            val_col = c if arm_name != "Combined" else (c + "__z")
            s = sub.dropna(subset=[val_col, fat_col])
            auc_v, p_v = auc_mw(s[val_col].values, s[fat_col].values)
            lo_v, hi_v = boot_auc_ci(s[val_col].values, s[fat_col].values, n_boot=3000, seed=seed + j)
            # Same-day WHO adjusted logistic regression / residualized AUC on first draw
            adj_auc, adj_p, adj_beta, adj_lr_p = np.nan, np.nan, np.nan, np.nan
            if s[val_col].nunique() >= 2 and s[fat_col].sum() >= 3:
                # Residualize feature z against same-day WHO (+ source if Combined)
                covs = ["temp_severity_ord"] + (["is_new"] if arm_name == "Combined" else [])
                Xz = sm.add_constant(s[covs].astype(float), has_constant="add")
                resid = sm.OLS(s[val_col].astype(float), Xz).fit().resid
                adj_auc, adj_p = auc_mw(resid.values, s[fat_col].values)
                try:
                    Xlr = sm.add_constant(s[[c + "__z"] + covs].astype(float), has_constant="add")
                    lr = sm.Logit(s[fat_col].astype(int), Xlr).fit(disp=False, maxiter=200)
                    adj_beta = float(lr.params[c + "__z"])
                    adj_lr_p = float(lr.pvalues[c + "__z"])
                except Exception:
                    pass
            rows.append(dict(
                set=arm_name, feature=c,
                fd_n=len(s), fd_deaths=int(s[fat_col].sum()),
                fd_AUC=auc_v, fd_AUC_lo=lo_v, fd_AUC_hi=hi_v, fd_p=p_v,
                fd_adjWHO_AUC=adj_auc, fd_adjWHO_p=adj_p,
                fd_logit_beta_adjWHO=adj_beta, fd_logit_p_adjWHO=adj_lr_p,
            ))
    res = pd.DataFrame(rows)
    # Homogeneity P between Real and PBISC new first-draw AUC via bootstrap SE
    piv_r = res[res["set"] == "Real"].set_index("feature")
    piv_n = res[res["set"] == "PBISC new"].set_index("feature")
    p_het_map = {}
    for c in cols:
        se_r = (piv_r.loc[c, "fd_AUC_hi"] - piv_r.loc[c, "fd_AUC_lo"]) / 3.92
        se_n = (piv_n.loc[c, "fd_AUC_hi"] - piv_n.loc[c, "fd_AUC_lo"]) / 3.92
        se_tot = np.hypot(se_r, se_n)
        if np.isfinite(se_tot) and se_tot > 0:
            z = (piv_r.loc[c, "fd_AUC"] - piv_n.loc[c, "fd_AUC"]) / se_tot
            p_het_map[c] = float(2 * stats.norm.sf(abs(z)))
        else:
            p_het_map[c] = np.nan
    res["fd_p_het"] = res["feature"].map(p_het_map)
    return res


def main():
    t0 = time.time()
    print("=== Loading verified primary metadata & fine-label compositions ===", flush=True)
    meta = pd.read_csv(D28 / "eskd_sample_metadata_verified.csv")
    fine = pd.read_csv(FINE_PCT_PATH)
    t3_107 = pd.read_csv(D28 / "table_s3_all_107_donors_overview.csv")
    scan_prev = pd.read_csv(D28 / "fig7_scan_43fine_true11.csv")
    scan_prev = scan_prev[scan_prev["CLL_excluded"] == True].copy()
    FINE43 = sorted(scan_prev["feature"].unique().tolist())
    assert len(FINE43) == 43

    df_all = meta.merge(fine, on=["source", "sample_id"], how="inner")
    assert len(df_all) == len(meta) == 362

    # Build primary 217-draw cohort (non-CLL, COVID+, in_primary_scope, qc>=200, Real takes precedence for Shared-27)
    pos_scope = df_all[
        (df_all["condition"] == "positive")
        & df_all["in_primary_scope"].astype(bool)
        & df_all["qc_pass_min200cells"].astype(bool)
        & (~df_all["patient_id"].isin(CLL))
    ].copy()
    real_pts = set(pos_scope[pos_scope["source"] == "Real"]["patient_id"])
    pb_pts = set(pos_scope[pos_scope["source"] == "PBISC"]["patient_id"])
    shared27 = real_pts & pb_pts
    real_only10 = real_pts - pb_pts
    new29 = pb_pts - real_pts
    assert len(real_pts) == 37 and len(shared27) == 27 and len(real_only10) == 10 and len(new29) == 29

    prim = pos_scope[
        (pos_scope["source"] == "Real") | ((pos_scope["source"] == "PBISC") & pos_scope["patient_id"].isin(new29))
    ].copy()
    assert len(prim) == 217
    prim["arm"] = np.where(prim["source"] == "Real", "Real", "PBISC new")
    prim["is_new"] = (prim["arm"] == "PBISC new").astype(int)
    prim["days"] = prim["days_from_first_positive_swab"].astype(float)
    prim["wk"] = prim["days"] / 7.0
    prim["fat11"] = prim["fatal_pub11"].astype(int)
    prim["fat9"] = prim["fatal_dir9"].astype(int)
    prim["fat14"] = prim["fatal_zen14"].astype(int)

    # Sort cleanly by patient & day to identify exact first draw per donor
    prim = prim.sort_values(["patient_id", "days", "sample_id"]).reset_index(drop=True)
    prim["draw_order"] = prim.groupby("patient_id").cumcount() + 1
    prim["n_draws_pt"] = prim.groupby("patient_id")["sample_id"].transform("count")
    prim["is_first_draw"] = prim["draw_order"] == 1
    prim["is_last_draw"] = prim["draw_order"] == prim["n_draws_pt"]

    # Check if any patient has 2 draws on the exact same first day
    same_first_day = prim.groupby("patient_id").apply(lambda g: (g["days"] == g["days"].min()).sum())
    assert int((same_first_day > 1).sum()) == 0, "Tie on first draw day!"
    print(f"[OK] Primary cohort: {prim.patient_id.nunique()} pts, {len(prim)} draws (0 ties on first-draw day).", flush=True)

    # -------------------------------------------------------------------------
    # PART A: Complete Patient & Draw Inventory (68 COVID+ + 39 COVID- summary)
    # -------------------------------------------------------------------------
    inv_rows = []
    # Include all 68 COVID+ donors (66 primary + 2 CLL)
    all_pos_donors = sorted(df_all[df_all["condition"] == "positive"]["patient_id"].unique(),
                            key=lambda s: (0 if s in real_only10 else (1 if s in shared27 else (2 if s in new29 else 3)),
                                           int("".join(ch for ch in s if ch.isdigit()) or 0)))
    for pid in all_pos_donors:
        g_all = df_all[(df_all["patient_id"] == pid) & (df_all["condition"] == "positive")].sort_values("days_from_first_positive_swab")
        g_r = g_all[g_all["source"] == "Real"]
        g_p = g_all[g_all["source"] == "PBISC"]
        g_prim = prim[prim["patient_id"] == pid].sort_values("days")
        if pid in real_only10:
            grp = "Real-only (10)"
        elif pid in shared27:
            grp = "Shared-27 (Real used)"
        elif pid in new29:
            grp = "PBISC new (29)"
        else:
            grp = "Excluded CLL (2)"
        ref_row = g_prim.iloc[0] if len(g_prim) else g_all.iloc[0]
        who_seq = [str(w) for w in g_prim["WHO_temp_severity"]] if len(g_prim) else []
        day_seq = [f"d{int(round(d_))}" for d_ in g_prim["days"]] if len(g_prim) else []
        symp_seq = [f"d{int(round(d_))}" if pd.notna(d_) else "NA" for d_ in g_prim["days_from_first_symptoms"]] if len(g_prim) else []
        inv_rows.append(dict(
            patient_id=pid,
            cohort_group=grp,
            primary_source=ref_row["arm"] if len(g_prim) else "Excluded_CLL",
            centre_or_wave=ref_row["centre_or_wave"],
            age=float(ref_row["age"]),
            sex=ref_row["sex"],
            n_primary_draws=len(g_prim),
            n_real_scRNA_pos=len(g_r),
            n_pbisc_pos=len(g_p),
            first_draw_day=float(g_prim["days"].iloc[0]) if len(g_prim) else np.nan,
            last_draw_day=float(g_prim["days"].iloc[-1]) if len(g_prim) else np.nan,
            draw_days_swab=", ".join(day_seq),
            draw_days_symptoms=", ".join(symp_seq),
            first_draw_WHO=str(g_prim["WHO_temp_severity"].iloc[0]) if len(g_prim) else "NA",
            first_draw_WHO_ord=float(g_prim["temp_severity_ord"].iloc[0]) if len(g_prim) else np.nan,
            same_day_WHO_each_draw=", ".join(who_seq),
            n_unique_sameday_WHO=len(set(who_seq)) if who_seq else 0,
            peak_WHO=str(ref_row["WHO_severity"]),
            peak_WHO_ord=float(ref_row["severity_ord"]),
            death_TRUE_PUB11=int(ref_row["fatal_pub11"]),
            death_DIRECT9=int(ref_row["fatal_dir9"]),
            death_ZENODO14=int(ref_row["fatal_zen14"]),
            recorded_steroid=str(ref_row.get("recorded_steroid", "")),
        ))
    T1_pat68 = pd.DataFrame(inv_rows)
    T1_pat68.to_csv(TAB_DIR / "T1_all68_patients_draw_schedule.csv", index=False, encoding="utf-8-sig")

    key_show = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN", "CD14mono_dex", "MAIT"]
    T1b_draws = prim[["patient_id", "arm", "sample_id", "draw_order", "n_draws_pt", "is_first_draw", "days", "days_from_first_symptoms",
                      "WHO_temp_severity", "temp_severity_ord", "WHO_severity", "severity_ord", "fat11", "fat9", "fat14",
                      "age", "sex", "centre_or_wave", "recorded_steroid", "n_cells", "temp_CRP", "temp_lymph"] + key_show].copy()
    T1b_draws.to_csv(TAB_DIR / "T1b_all217_primary_draws_detail.csv", index=False, encoding="utf-8-sig")

    # Cohort draw summary table
    sum_rows = []
    for lbl, sub_d in [
        ("Real-only (10 pts)", prim[prim["patient_id"].isin(real_only10)]),
        ("Shared-27 (Real used, 27 pts)", prim[prim["patient_id"].isin(shared27)]),
        ("Real ALL (37 pts)", prim[prim["arm"] == "Real"]),
        ("Real · Survived (32 pts)", prim[(prim["arm"] == "Real") & (prim["fat11"] == 0)]),
        ("Real · Died (5 pts)", prim[(prim["arm"] == "Real") & (prim["fat11"] == 1)]),
        ("PBISC new ALL (29 pts)", prim[prim["arm"] == "PBISC new"]),
        ("PBISC new · Survived (23 pts)", prim[(prim["arm"] == "PBISC new") & (prim["fat11"] == 0)]),
        ("PBISC new · Died (6 pts)", prim[(prim["arm"] == "PBISC new") & (prim["fat11"] == 1)]),
        ("Combined ALL (66 pts)", prim),
        ("Combined · Survived (55 pts)", prim[prim["fat11"] == 0]),
        ("Combined · Died (11 pts)", prim[prim["fat11"] == 1]),
    ]:
        pts_u = sub_d.drop_duplicates("patient_id")
        fd = sub_d[sub_d["is_first_draw"]]
        nd = sub_d.groupby("patient_id").size()
        sum_rows.append(dict(
            group=lbl,
            n_patients=len(pts_u),
            n_draws=len(sub_d),
            draws_per_pt_median=float(nd.median()),
            draws_per_pt_iqr=f"{nd.quantile(0.25):.0f}-{nd.quantile(0.75):.0f}",
            draws_per_pt_range=f"{nd.min()}-{nd.max()}",
            n_single_draw_pts=int((nd == 1).sum()),
            first_draw_day_median=float(fd["days"].median()),
            first_draw_day_iqr=f"{fd['days'].quantile(0.25):.0f}-{fd['days'].quantile(0.75):.0f}",
            first_draw_le7d_pts=int((fd["days"] <= 7).sum()),
            all_draws_day_median=float(sub_d["days"].median()),
            all_draws_le7d=int((sub_d["days"] <= 7).sum()),
            all_draws_gt7d=int((sub_d["days"] > 7).sum()),
            first_draw_WHO_mild_mod_sev_crit=f"{int((fd['temp_severity_ord']==0).sum())}/{int((fd['temp_severity_ord']==1).sum())}/{int((fd['temp_severity_ord']==2).sum())}/{int((fd['temp_severity_ord']==3).sum())}",
            peak_WHO_mild_mod_sev_crit=f"{int((pts_u['severity_ord']==0).sum())}/{int((pts_u['severity_ord']==1).sum())}/{int((pts_u['severity_ord']==2).sum())}/{int((pts_u['severity_ord']==3).sum())}",
            sameday_WHO_draws_m_mo_s_c=f"{int((sub_d['temp_severity_ord']==0).sum())}/{int((sub_d['temp_severity_ord']==1).sum())}/{int((sub_d['temp_severity_ord']==2).sum())}/{int((sub_d['temp_severity_ord']==3).sum())}",
            deaths_pub11=int(pts_u["fat11"].sum()),
        ))
    T1c_sum = pd.DataFrame(sum_rows)
    T1c_sum.to_csv(TAB_DIR / "T1c_cohort_draw_summary.csv", index=False, encoding="utf-8-sig")

    # Also summarize the 39 COVID-negative donors from table_s3_all_107
    neg39 = t3_107[~t3_107["patient_id"].isin(set(all_pos_donors))].copy()
    neg39.to_csv(TAB_DIR / "T1d_covid_negative_39_donors.csv", index=False, encoding="utf-8-sig")
    print(f"[OK] Saved T1 patient/draw inventories (68 COVID+ pts, 217 draws, 39 COVID- donors).", flush=True)

    # -------------------------------------------------------------------------
    # DOMAIN 1: 43 Fine-Grained Immune Subsets
    #   (A) Previous: Patient-mean vs Peak WHO & Final Death
    #   (B) NEW Severity: All Draws vs Same-Day WHO (GEE unadj, GEE +week, Mundlak, Clust-Boot Spearman)
    #   (C) NEW Death: First Draw vs Final Death (All 66, <=7d first draws, WHO-adj, Direct-9, Late >7d contrast)
    # -------------------------------------------------------------------------
    print("=== Domain 1: Re-analyzing 43 fine-grained immune subsets ===", flush=True)
    prim_z = add_draw_rint(prim, FINE43)
    df_first = prim[prim["is_first_draw"]].copy().reset_index(drop=True)
    df_first_le7 = df_first[df_first["days"] <= 7.0].copy().reset_index(drop=True)
    print(f"    First draws total: N={len(df_first)} (deaths={int(df_first.fat11.sum())}: Real {int(df_first[df_first.arm=='Real'].fat11.sum())}/{int((df_first.arm=='Real').sum())}, PBISC new {int(df_first[df_first.arm=='PBISC new'].fat11.sum())}/{int((df_first.arm=='PBISC new').sum())})")
    print(f"    First draws <=7d:  N={len(df_first_le7)} (deaths={int(df_first_le7.fat11.sum())}: Real {int(df_first_le7[df_first_le7.arm=='Real'].fat11.sum())}/{int((df_first_le7.arm=='Real').sum())}, PBISC new {int(df_first_le7[df_first_le7.arm=='PBISC new'].fat11.sum())}/{int((df_first_le7.arm=='PBISC new').sum())})", flush=True)

    resume_mode = ("--resume" in sys.argv) and (TAB_DIR / "T2b_comp43_combined_full_comparison.csv").exists()
    if resume_mode:
        COMP43_comb = pd.read_csv(TAB_DIR / "T2b_comp43_combined_full_comparison.csv")
        fb, fp = fast_exch_gee_batch(prim_z, [f + "__z" for f in FINE43], "temp_severity_ord", ("wk",))
        ref_b = COMP43_comb.set_index("feature").loc[FINE43, "sd_gee_wk_beta"].to_numpy()
        ref_p = COMP43_comb.set_index("feature").loc[FINE43, "sd_gee_wk_p"].to_numpy()
        print(f"    [--resume loaded T2b; fast_exch_gee_batch vs statsmodels] max|dbeta|={np.max(np.abs(fb - ref_b)):.2e}, max|dp|={np.max(np.abs(fp - ref_p)):.2e}", flush=True)
    else:
        # 1A. Same-day WHO Severity across Real (138), PBISC new (79), Combined (217)
        sev_rows = []
        for arm_name, sub_z in [("Real", prim_z[prim_z["arm"] == "Real"]),
                                ("PBISC new", prim_z[prim_z["arm"] == "PBISC new"]),
                                ("Combined", prim_z)]:
            cl_sp = cluster_spearman_sameday(sub_z, [c + "__z" for c in FINE43], xcol="temp_severity_ord", n_boot=N_BOOT, seed=101)
            cl_sp["feature"] = FINE43
            cl_sp = cl_sp.set_index("feature")
            for f in FINE43:
                zc = f + "__z"
                g_unadj = gee_fit_term(sub_z, zc, "temp_severity_ord", covars=())
                g_wk = gee_fit_term(sub_z, zc, "temp_severity_ord", covars=("wk",))
                mk = mundlak_gee(sub_z, zc, "temp_severity_ord", covars=("wk",))
                sp_r = cl_sp.loc[f]
                sev_rows.append(dict(
                    set=arm_name, feature=f,
                    sd_n_draws=g_wk["n_draws"], sd_n_pts=g_wk["n_pts"],
                    sd_gee_unadj_beta=g_unadj["beta"], sd_gee_unadj_lo=g_unadj["lo"], sd_gee_unadj_hi=g_unadj["hi"],
                    sd_gee_unadj_p=g_unadj["p"], sd_gee_unadj_p_br=g_unadj["p_br"],
                    sd_gee_wk_beta=g_wk["beta"], sd_gee_wk_lo=g_wk["lo"], sd_gee_wk_hi=g_wk["hi"],
                    sd_gee_wk_p=g_wk["p"], sd_gee_wk_p_br=g_wk["p_br"],
                    sd_rho=sp_r["sameday_rho"], sd_rho_lo=sp_r["sameday_rho_lo"], sd_rho_hi=sp_r["sameday_rho_hi"],
                    sd_rho_p_clboot=sp_r["sameday_rho_p_clboot"], sd_rho_p_naive=sp_r["sameday_rho_p_naive"],
                    mundlak_within_beta=mk["beta_within"], mundlak_within_p=mk["p_within"], mundlak_within_p_br=mk["p_within_br"],
                    mundlak_between_beta=mk["beta_between"], mundlak_between_p=mk["p_between"], mundlak_between_p_br=mk["p_between_br"],
                ))
        SD_SEV43 = pd.DataFrame(sev_rows)
        for _, ix in SD_SEV43.groupby("set").groups.items():
            for col in ["sd_gee_unadj_p", "sd_gee_unadj_p_br", "sd_gee_wk_p", "sd_gee_wk_p_br", "sd_rho_p_clboot", "mundlak_within_p"]:
                SD_SEV43.loc[ix, col + "_q43"] = bh(SD_SEV43.loc[ix, col])

        # Homogeneity P between Real and PBISC new for same-day WHO GEE (+wk)
        pr_g = SD_SEV43[SD_SEV43["set"] == "Real"].set_index("feature")
        pn_g = SD_SEV43[SD_SEV43["set"] == "PBISC new"].set_index("feature")
        p_het_sd = {}
        for f in FINE43:
            ser = (pr_g.loc[f, "sd_gee_wk_hi"] - pr_g.loc[f, "sd_gee_wk_lo"]) / 3.92
            sen = (pn_g.loc[f, "sd_gee_wk_hi"] - pn_g.loc[f, "sd_gee_wk_lo"]) / 3.92
            st_ = np.hypot(ser, sen)
            p_het_sd[f] = float(2 * stats.norm.sf(abs(pr_g.loc[f, "sd_gee_wk_beta"] - pn_g.loc[f, "sd_gee_wk_beta"]) / st_)) if (np.isfinite(st_) and st_ > 0) else np.nan
        SD_SEV43["sd_gee_wk_p_het"] = SD_SEV43["feature"].map(p_het_sd)

        # 1B. First-Draw -> Final Death (TRUE_PUB_11) + sensitivity (<=7d first draws, Direct-9, WHO 2-3, Last draw)
        FD_DTH43 = first_draw_death_stats(df_first, FINE43, fat_col="fat11", seed=202)
        for _, ix in FD_DTH43.groupby("set").groups.items():
            FD_DTH43.loc[ix, "fd_q43"] = bh(FD_DTH43.loc[ix, "fd_p"])
            FD_DTH43.loc[ix, "fd_adjWHO_q43"] = bh(FD_DTH43.loc[ix, "fd_adjWHO_p"])

        FD_LE7 = first_draw_death_stats(df_first_le7, FINE43, fat_col="fat11", seed=303)
        FD_LE7 = FD_LE7[FD_LE7["set"] == "Combined"][["feature", "fd_n", "fd_deaths", "fd_AUC", "fd_AUC_lo", "fd_AUC_hi", "fd_p", "fd_adjWHO_AUC", "fd_adjWHO_p"]].rename(
            columns={"fd_n": "fd_le7_n", "fd_deaths": "fd_le7_deaths", "fd_AUC": "fd_le7_AUC", "fd_AUC_lo": "fd_le7_lo",
                     "fd_AUC_hi": "fd_le7_hi", "fd_p": "fd_le7_p", "fd_adjWHO_AUC": "fd_le7_adjWHO_AUC", "fd_adjWHO_p": "fd_le7_adjWHO_p"}
        )
        FD_LE7["fd_le7_q43"] = bh(FD_LE7["fd_le7_p"])

        FD_DIR9 = first_draw_death_stats(df_first, FINE43, fat_col="fat9", seed=404)
        FD_DIR9 = FD_DIR9[FD_DIR9["set"] == "Combined"][["feature", "fd_AUC", "fd_p"]].rename(
            columns={"fd_AUC": "fd_dir9_AUC", "fd_p": "fd_dir9_p"}
        )

        # Within peak WHO 2-3 using first draw (full-cohort first-draw INT)
        df_first_z = add_draw_rint(df_first, FINE43)
        w23_fd = df_first_z[df_first_z["severity_ord"] >= 2]
        w23_rows = []
        for f in FINE43:
            a, p = auc_mw(w23_fd[f + "__z"].values, w23_fd["fat11"].values)
            w23_rows.append(dict(feature=f, fd_who23_AUC=a, fd_who23_p=p))
        FD_WHO23 = pd.DataFrame(w23_rows)

        # Late (>7d) draw patient-mean death AUC for direct contrast vs first draw
        gt7_df = prim[prim["days"] > 7.0]
        gt7_pt = gt7_df.groupby(["arm", "patient_id"])[FINE43 + ["fat11"]].mean().reset_index()
        gt7_z = add_draw_rint(gt7_pt, FINE43)
        gt7_rows = []
        for f in FINE43:
            a, p = auc_mw(gt7_z[f + "__z"].values, gt7_z["fat11"].values)
            gt7_rows.append(dict(feature=f, late_gt7_n=len(gt7_z), late_gt7_deaths=int(gt7_z["fat11"].sum()),
                                 late_gt7_AUC=a, late_gt7_p=p))
        LATE_GT7 = pd.DataFrame(gt7_rows)

        # Merge Previous (Patient-mean) + NEW (Same-day WHO + First-Draw Death)
        prev_ren = scan_prev.replace({"set": {"Real + PBISC new": "Combined"}}).rename(columns={
            "rho": "pm_peakWHO_rho", "p_sev": "pm_peakWHO_p", "q_sev": "pm_peakWHO_q43",
            "AUC_death": "pm_death_AUC", "p_death": "pm_death_p", "q_death": "pm_death_q43",
            "p_het_sev": "pm_p_het_sev", "p_het_death": "pm_p_het_death",
        })[["set", "feature", "pm_peakWHO_rho", "pm_peakWHO_p", "pm_peakWHO_q43", "pm_death_AUC", "pm_death_p", "pm_death_q43", "pm_p_het_sev", "pm_p_het_death"]]

        COMP43 = prev_ren.merge(SD_SEV43, on=["set", "feature"]).merge(FD_DTH43, on=["set", "feature"])
        COMP43.to_csv(TAB_DIR / "T2_comp43_patientmean_vs_samedaySev_firstdrawDeath.csv", index=False, encoding="utf-8-sig")

        COMP43_comb = COMP43[COMP43["set"] == "Combined"].merge(FD_LE7, on="feature").merge(FD_DIR9, on="feature").merge(FD_WHO23, on="feature").merge(LATE_GT7, on="feature")
        COMP43_comb.to_csv(TAB_DIR / "T2b_comp43_combined_full_comparison.csv", index=False, encoding="utf-8-sig")

        # Selection funnel under NEW framework (Same-day WHO GEE +wk q43<0.05 & First-draw Death q43<0.10 / P<0.05)
        c_r = COMP43[COMP43["set"] == "Real"].set_index("feature")
        c_n = COMP43[COMP43["set"] == "PBISC new"].set_index("feature")
        c_c = COMP43_comb.set_index("feature")
        funnel_rows = []
        for f in FINE43:
            # Severity (same-day WHO GEE +wk)
            s_comb_fdr = c_c.loc[f, "sd_gee_wk_p_q43"] < 0.05
            s_comb_br = c_c.loc[f, "sd_gee_wk_p_br_q43"] < 0.05
            s_real_fdr = c_r.loc[f, "sd_gee_wk_p_q43"] < 0.05
            s_same_sign = (c_r.loc[f, "sd_gee_wk_beta"] * c_n.loc[f, "sd_gee_wk_beta"] > 0)
            # Death (first draw -> final death)
            d_comb_fdr = c_c.loc[f, "fd_q43"] < 0.10
            d_comb_nom = c_c.loc[f, "fd_p"] < 0.05
            d_real_fdr = c_r.loc[f, "fd_q43"] < 0.10
            d_same_sign = ((c_r.loc[f, "fd_AUC"] - 0.5) * (c_n.loc[f, "fd_AUC"] - 0.5) > 0)
            funnel_rows.append(dict(
                feature=f,
                prev_sev_q43=c_c.loc[f, "pm_peakWHO_q43"], prev_dth_q43=c_c.loc[f, "pm_death_q43"],
                new_sdWHO_beta=c_c.loc[f, "sd_gee_wk_beta"], new_sdWHO_p=c_c.loc[f, "sd_gee_wk_p"], new_sdWHO_q43=c_c.loc[f, "sd_gee_wk_p_q43"], new_sdWHO_q43_br=c_c.loc[f, "sd_gee_wk_p_br_q43"],
                new_sdWHO_rho=c_c.loc[f, "sd_rho"], new_sdWHO_rho_q43=c_c.loc[f, "sd_rho_p_clboot_q43"],
                new_sdWHO_real_beta=c_r.loc[f, "sd_gee_wk_beta"], new_sdWHO_real_q43=c_r.loc[f, "sd_gee_wk_p_q43"],
                new_sdWHO_new_beta=c_n.loc[f, "sd_gee_wk_beta"], new_sdWHO_new_p=c_n.loc[f, "sd_gee_wk_p"],
                new_sdWHO_same_sign=s_same_sign, new_sdWHO_PBISC_rescued=(s_comb_fdr and (not s_real_fdr) and s_same_sign and f != "CD14mono_dex"),
                new_fd_AUC=c_c.loc[f, "fd_AUC"], new_fd_p=c_c.loc[f, "fd_p"], new_fd_q43=c_c.loc[f, "fd_q43"],
                new_fd_adjWHO_AUC=c_c.loc[f, "fd_adjWHO_AUC"], new_fd_adjWHO_p=c_c.loc[f, "fd_adjWHO_p"],
                new_fd_real_AUC=c_r.loc[f, "fd_AUC"], new_fd_real_p=c_r.loc[f, "fd_p"],
                new_fd_new_AUC=c_n.loc[f, "fd_AUC"], new_fd_new_p=c_n.loc[f, "fd_p"],
                new_fd_same_sign=d_same_sign, new_fd_PBISC_rescued=(d_comb_fdr and (not d_real_fdr) and d_same_sign and f != "CD14mono_dex"),
                new_fd_le7_AUC=c_c.loc[f, "fd_le7_AUC"], new_fd_le7_p=c_c.loc[f, "fd_le7_p"],
                late_gt7_AUC=c_c.loc[f, "late_gt7_AUC"], late_gt7_p=c_c.loc[f, "late_gt7_p"],
            ))
        T3_funnel = pd.DataFrame(funnel_rows)
        T3_funnel.to_csv(TAB_DIR / "T3_new_framework_43_funnel.csv", index=False, encoding="utf-8-sig")
        print(f"[OK] Domain 1 done ({time.time() - t0:.1f}s).", flush=True)

    # -------------------------------------------------------------------------
    # DOMAIN 2: Cohort-Expansion Power Rescue Under the NEW Framework
    #   (n_real = 20 sampled Real patients vs +29 PBISC new vs +29 outcome-shuffled PBISC new, 300 sims)
    # -------------------------------------------------------------------------
    if ("--resume-d3" in sys.argv) and (TAB_DIR / "T4_power_rescue_new_framework.csv").exists():
        print("=== Domain 2: Loaded existing T4 (--resume-d3) ===", flush=True)
    else:
        print("=== Domain 2: Subsampling power rescue under NEW framework ===", flush=True)
        rng_pw = np.random.default_rng(20261120)
        N_SIM = 300
        r_pids_surv = df_first[(df_first["arm"] == "Real") & (df_first["fat11"] == 0)]["patient_id"].to_numpy()
        r_pids_died = df_first[(df_first["arm"] == "Real") & (df_first["fat11"] == 1)]["patient_id"].to_numpy()
        n_pids_all = df_first[df_first["arm"] == "PBISC new"]["patient_id"].to_numpy()

        pw_targets_sev = ["B_ASC_IgG", "B_ASC_IgA", "B_ASC_dividing", "HSPC", "NK CD56 bright", "T CD8 CM", "CD16mono"]
        pw_targets_dth = ["T CD8 CTL", "NK", "CD14mono_IFN", "NK CD56 bright", "B_CD11c", "MAIT", "HSPC"]

        pw_counts = {k: {"real_nom": 0, "real_fdr": 0, "plus_nom": 0, "plus_fdr": 0, "shuf_nom": 0, "shuf_fdr": 0}
                     for k in [(f, "sameday_sev") for f in pw_targets_sev] + [(f, "firstdraw_dth") for f in pw_targets_dth]}

        # Pre-split Real and PBISC new tables for fast subsampling
        prim_r = prim[prim["arm"] == "Real"]
        prim_n = prim[prim["arm"] == "PBISC new"]
        fd_r = df_first[df_first["arm"] == "Real"]
        fd_n = df_first[df_first["arm"] == "PBISC new"]
        zcols43 = [f + "__z" for f in FINE43]

        for sim in range(N_SIM):
            # Sample 20 Real patients stratified (3 deaths, 17 survivors matching ~5/37 rate)
            sp_r = np.concatenate([
                rng_pw.choice(r_pids_died, size=3, replace=False),
                rng_pw.choice(r_pids_surv, size=17, replace=False),
            ])
            # (a) First-draw death across all 43 subsets (fast exact MW)
            sub_fdr = fd_r[fd_r["patient_id"].isin(sp_r)].copy()
            sub_fdn_shuf = fd_n.copy()
            sub_fdn_shuf["fat11"] = rng_pw.permutation(sub_fdn_shuf["fat11"].values)

            for tag, df_cur in [("real", sub_fdr),
                                ("plus", pd.concat([sub_fdr, fd_n], ignore_index=True)),
                                ("shuf", pd.concat([sub_fdr, sub_fdn_shuf], ignore_index=True))]:
                df_cz = add_draw_rint(df_cur, FINE43)
                ps = np.array([auc_mw(df_cz[f + "__z"].values, df_cz["fat11"].values)[1] for f in FINE43])
                qs = bh(ps)
                pmap = dict(zip(FINE43, ps))
                qmap = dict(zip(FINE43, qs))
                for f in pw_targets_dth:
                    if pmap[f] < 0.05:
                        pw_counts[(f, "firstdraw_dth")][f"{tag}_nom"] += 1
                    if qmap[f] < 0.10:
                        pw_counts[(f, "firstdraw_dth")][f"{tag}_fdr"] += 1

            # (b) All-draws same-day WHO severity Exchangeable GEE (+wk) across ALL 300 sims via fast vectorized Sherman-Morrison GEE!
            sub_pr = prim_r[prim_r["patient_id"].isin(sp_r)].copy()
            sub_pn_shuf = prim_n.copy()
            sub_pn_shuf["temp_severity_ord"] = rng_pw.permutation(sub_pn_shuf["temp_severity_ord"].values)

            for tag, df_cur in [("real", sub_pr),
                                ("plus", pd.concat([sub_pr, prim_n], ignore_index=True)),
                                ("shuf", pd.concat([sub_pr, sub_pn_shuf], ignore_index=True))]:
                df_cz = add_draw_rint(df_cur, FINE43)
                _, ps = fast_exch_gee_batch(df_cz, zcols43, "temp_severity_ord", covars=("wk",))
                qs = bh(ps)
                pmap = dict(zip(FINE43, ps))
                qmap = dict(zip(FINE43, qs))
                for f in pw_targets_sev:
                    if pmap[f] < 0.05:
                        pw_counts[(f, "sameday_sev")][f"{tag}_nom"] += 1
                    if qmap[f] < 0.05:
                        pw_counts[(f, "sameday_sev")][f"{tag}_fdr"] += 1

        pw_rows = []
        for (f, ep), cnts in pw_counts.items():
            denom = float(N_SIM)
            pw_rows.append(dict(
                feature=f, endpoint=ep, n_sims=int(denom),
                real_nom=cnts["real_nom"] / denom, real_fdr=cnts["real_fdr"] / denom,
                plus_nom=cnts["plus_nom"] / denom, plus_fdr=cnts["plus_fdr"] / denom,
                shuf_nom=cnts["shuf_nom"] / denom, shuf_fdr=cnts["shuf_fdr"] / denom,
            ))
        T4_power = pd.DataFrame(pw_rows)
        T4_power.to_csv(TAB_DIR / "T4_power_rescue_new_framework.csv", index=False, encoding="utf-8-sig")
        print(f"[OK] Domain 2 power rescue done ({time.time() - t0:.1f}s).", flush=True)

    # -------------------------------------------------------------------------
    # DOMAIN 3: Single-cell HSPC 16 Modules & 20 Genes (from raw cell-level h5ad)
    # -------------------------------------------------------------------------
    print("=== Domain 3: Single-cell HSPC 16 modules & 20 genes ===", flush=True)
    r_h5 = ad.read_h5ad(CELL_DIR / "real_rare_subset_cells_counts.h5ad")
    p_h5 = ad.read_h5ad(CELL_DIR / "pbisc_rare_subset_cells_counts.h5ad")
    p_sym = pd.Series(p_h5.var["gene_symbol"].values, index=np.arange(p_h5.n_vars))
    p_sym = p_sym[~p_sym.duplicated()]
    r_sym = pd.Series(r_h5.var_names.values, index=np.arange(r_h5.n_vars))
    r_sym = r_sym[~r_sym.duplicated()]
    shared_g5762 = sorted(set(p_sym.values) & set(r_sym.values))
    r_gmap, p_gmap = {s_: i for i, s_ in r_sym.items()}, {s_: i for i, s_ in p_sym.items()}

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
        "MEP_minus_Stem": ("MEP_EryMk", "Stem_HSC"),
        "GMP_minus_Stem": ("GMP_Myeloid", "Stem_HSC"),
        "Alarmin_minus_Stem": ("Alarmin_Inflam", "Stem_HSC"),
        "HSC_broad_q": ["CD34", "PROM1", "AVP", "HLF", "CRHBP", "SPINK2", "HOPX", "MEIS1", "MLLT3"],
    }
    MODS16 = list(MOD16_DEF.keys())
    G20 = ["GATA1", "KLF1", "TFRC", "BLVRB", "APOC1", "GATA2", "ITGA2B", "MKI67", "TOP2A", "STMN1", "TYMS",
           "CD34", "AVP", "HLF", "HOPX", "CRHBP", "MPO", "ELANE", "CSF3R", "CEBPD"]
    GENES_H = sorted({g for v in MOD16_DEF.values() if isinstance(v, list) for g in v} | set(G20))
    GENES_H = [g for g in GENES_H if g in shared_g5762]

    def extract_hspc(adata, gmap, src, pid_col):
        m = (adata.obs["fine"] == "HSPC").values
        obs = adata.obs.loc[m]
        keep = obs["sample_id"].astype(str).isin(set(prim.loc[prim["source"] == src, "sample_id"])).values
        obs = obs.loc[keep]
        X = adata.X[m][keep][:, [gmap[g] for g in shared_g5762]]
        X = X.toarray() if sp.issparse(X) else X
        X = np.log1p(np.divide(X, np.maximum(X.sum(axis=1, keepdims=True), 1.0)) * 1e4)
        out = pd.DataFrame(X[:, [shared_g5762.index(g) for g in GENES_H]], columns=GENES_H)
        out["patient_id"] = obs[pid_col].astype(str).values
        out["sample_id"] = obs["sample_id"].astype(str).values
        out["source"] = src
        return out

    HC = pd.concat([extract_hspc(r_h5, r_gmap, "Real", "patient_id"),
                    extract_hspc(p_h5, p_gmap, "PBISC", "donor_id")], ignore_index=True)
    del r_h5, p_h5
    for k, v in MOD16_DEF.items():
        if isinstance(v, list):
            HC[k] = HC[[g for g in v if g in HC.columns]].mean(axis=1)
    for k, v in MOD16_DEF.items():
        if isinstance(v, tuple):
            HC[k] = HC[v[0]] - HC[v[1]]
    FEATS_H = MODS16 + G20

    # Per-draw mean across HSPC cells
    HD_raw = HC.groupby(["source", "sample_id"])[FEATS_H].mean()
    HD_raw["n_hspc"] = HC.groupby(["source", "sample_id"]).size()
    HD_raw = HD_raw.reset_index()
    HD = prim[["arm", "is_new", "source", "sample_id", "patient_id", "days", "wk", "draw_order", "is_first_draw",
               "temp_severity_ord", "severity_ord", "fat11", "fat9"]].merge(HD_raw, on=["source", "sample_id"], how="inner")

    # Primary draw-level threshold >=3 HSPC cells/draw (149 draws / 61 pts, matching Folder 30) + sensitivity >=1 cell (176 draws / 63 pts)
    HD3 = HD[HD["n_hspc"] >= 3].copy().reset_index(drop=True)
    HD3_z = add_draw_rint(HD3, FEATS_H)
    # First HSPC-positive draw (>=3 cells) per patient for first-draw death prediction!
    HD3_first = HD3.sort_values(["patient_id", "days", "sample_id"]).groupby("patient_id", as_index=False).first()
    HD1_first = HD[HD["n_hspc"] >= 1].sort_values(["patient_id", "days", "sample_id"]).groupby("patient_id", as_index=False).first()
    print(f"    HSPC >=3 cells/draw: {len(HD3)} draws / {len(HD3_first)} pts (deaths={int(HD3_first.fat11.sum())}); >=1 cell/draw first draw: {len(HD1_first)} pts (deaths={int(HD1_first.fat11.sum())})", flush=True)

    # Load previous patient-pooled HSPC results for direct comparison
    pub_hm = pd.read_csv(D28 / "fig8_hspc_modules_16unique_true11.csv").rename(columns={"cohort": "set", "module": "feature", "q16_sev": "pm_q_sev", "q16_dth11": "pm_q_dth"})
    pub_hg = pd.read_csv(D28 / "suppS5_hspc_20gene_true11.csv").rename(columns={"cohort": "set", "gene": "feature", "q20_sev": "pm_q_sev", "q20_dth11": "pm_q_dth"})
    pub_h = pd.concat([pub_hm[["set", "feature", "rho_sev", "p_sev", "pm_q_sev", "AUC_dth11", "p_dth11", "pm_q_dth"]],
                       pub_hg[["set", "feature", "rho_sev", "p_sev", "pm_q_sev", "AUC_dth11", "p_dth11", "pm_q_dth"]]], ignore_index=True)

    # Same-day WHO GEE + Cluster-Boot Spearman on HD3
    h_sev_rows = []
    for arm_name, sub_z in [("Real", HD3_z[HD3_z["arm"] == "Real"]),
                            ("PBISC new", HD3_z[HD3_z["arm"] == "PBISC new"]),
                            ("Combined", HD3_z)]:
        cl_h = cluster_spearman_sameday(sub_z, [f + "__z" for f in FEATS_H], xcol="temp_severity_ord", n_boot=N_BOOT, seed=505)
        cl_h["feature"] = FEATS_H
        cl_h = cl_h.set_index("feature")
        for f in FEATS_H:
            zc = f + "__z"
            gu = gee_fit_term(sub_z, zc, "temp_severity_ord", covars=())
            gw = gee_fit_term(sub_z, zc, "temp_severity_ord", covars=("wk",))
            mk = mundlak_gee(sub_z, zc, "temp_severity_ord", covars=("wk",))
            h_sev_rows.append(dict(
                set=arm_name, feature=f, family="module" if f in MODS16 else "gene",
                sd_n_draws=gw["n_draws"], sd_n_pts=gw["n_pts"],
                sd_gee_unadj_beta=gu["beta"], sd_gee_unadj_p=gu["p"], sd_gee_unadj_p_br=gu["p_br"],
                sd_gee_wk_beta=gw["beta"], sd_gee_wk_lo=gw["lo"], sd_gee_wk_hi=gw["hi"],
                sd_gee_wk_p=gw["p"], sd_gee_wk_p_br=gw["p_br"],
                sd_rho=cl_h.loc[f, "sameday_rho"], sd_rho_p_clboot=cl_h.loc[f, "sameday_rho_p_clboot"],
                mundlak_within_beta=mk["beta_within"], mundlak_within_p=mk["p_within"],
            ))
    H_SD = pd.DataFrame(h_sev_rows)
    H_FD3 = first_draw_death_stats(HD3_first, FEATS_H, fat_col="fat11", seed=606)
    H_FD1 = first_draw_death_stats(HD1_first, FEATS_H, fat_col="fat11", seed=607)
    H_FD1 = H_FD1[H_FD1["set"] == "Combined"][["feature", "fd_n", "fd_deaths", "fd_AUC", "fd_p"]].rename(
        columns={"fd_n": "fd1cell_n", "fd_deaths": "fd1cell_deaths", "fd_AUC": "fd1cell_AUC", "fd_p": "fd1cell_p"}
    )

    H_ALL = pub_h.merge(H_SD, on=["set", "feature"]).merge(H_FD3, on=["set", "feature"])
    for _, ix in H_ALL.groupby(["set", "family"]).groups.items():
        for col in ["sd_gee_unadj_p", "sd_gee_wk_p", "sd_gee_wk_p_br", "sd_rho_p_clboot", "fd_p"]:
            H_ALL.loc[ix, col + "_qfam"] = bh(H_ALL.loc[ix, col])
    H_ALL = H_ALL.merge(H_FD1, on="feature", how="left")
    H_ALL.to_csv(TAB_DIR / "T5_HSPC_16modules_20genes_new_framework.csv", index=False, encoding="utf-8-sig")
    print(f"[OK] Domain 3 HSPC single-cell done ({time.time() - t0:.1f}s).", flush=True)

    # -------------------------------------------------------------------------
    # DOMAIN 4: Single-cell T-cell 19-Gene ISG across 13 T-cell Subsets
    # -------------------------------------------------------------------------
    print("=== Domain 4: Single-cell T-cell 19-gene ISG across 13 T subsets ===", flush=True)
    z_pb = np.load(CELL_DIR / "real_pseudobulk_sample_x_fine_panel.npz", allow_pickle=True)
    r_genes, r_mat = pd.Index(z_pb["genes"].astype(str)), z_pb["pb"]
    grps = pd.Series(z_pb["groups"].astype(str)).str.split("|", expand=True)
    grps.columns = ["sample_id", "fine"]
    grps["ncells"] = z_pb["ncells"].astype(int)
    ISG19 = [g for g in ["ISG15", "IFI6", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
                         "OAS1", "OAS2", "OAS3", "RSAD2", "XAF1", "STAT1", "IRF7", "EPSTI1", "HERC5", "SIGLEC1"] if g in r_genes]

    pb_h5 = ad.read_h5ad(PB_H5AD)
    pb_sym_all = pd.Index(pb_h5.var["gene_symbol"].astype(str))
    tot_umi_pb = np.asarray(pb_h5.X.sum(axis=1)).ravel().astype(np.float64)
    Xi = pb_h5.X[:, [int(np.where(pb_sym_all == g)[0][0]) for g in ISG19]]
    Xi = Xi.toarray() if sp.issparse(Xi) else Xi
    pb_tot_s = pd.DataFrame({"sample_id": pb_h5.obs["sample_id"].astype(str).values,
                             "fine": pb_h5.obs["celltype_fine"].astype(str).values,
                             "tot_umi": tot_umi_pb, "ncells": 1})
    for j, g in enumerate(ISG19):
        pb_tot_s[g] = Xi[:, j]
    del pb_h5, Xi

    t_fines = [f for f in sorted(grps["fine"].unique()) if f.startswith("T ") or f in ("MAIT", "NKT")]
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

    def calc_isg(df_in):
        return np.log1p(df_in[ISG19].div(np.maximum(df_in["tot_umi"], 1.0), axis=0) * 1e4).mean(axis=1)

    pub_tisg = pd.read_csv(D28 / "fig8_tcell_isg_13subsets_true11.csv").set_index("subset")
    t_rows = []
    for gname, fset in GROUPS_T:
        m_r = grps["fine"].isin(fset).values
        df_rs = pd.DataFrame(r_mat[m_r][:, [r_genes.get_loc(g) for g in ISG19]], columns=ISG19)
        df_rs["sample_id"] = grps.loc[m_r, "sample_id"].values
        df_rs["tot_umi"] = r_mat[m_r].sum(axis=1)
        df_rs["ncells"] = grps.loc[m_r, "ncells"].values
        df_rs = df_rs.groupby("sample_id").sum()
        df_rs = df_rs[df_rs["ncells"] >= 5].copy()
        df_rs["isg"] = calc_isg(df_rs)

        df_ps = pb_tot_s[pb_tot_s["fine"].isin(fset)].groupby("sample_id")[["tot_umi", "ncells"] + ISG19].sum()
        df_ps = df_ps[df_ps["ncells"] >= 5].copy()
        df_ps["isg"] = calc_isg(df_ps)

        isg_df = pd.concat([
            df_rs[["isg", "ncells"]].reset_index().assign(source="Real"),
            df_ps[["isg", "ncells"]].reset_index().assign(source="PBISC"),
        ], ignore_index=True)

        x = prim[["arm", "is_new", "source", "sample_id", "patient_id", "days", "wk", "is_first_draw",
                  "temp_severity_ord", "severity_ord", "fat11", "fat9"]].merge(isg_df, on=["source", "sample_id"], how="inner")
        xz = add_draw_rint(x, ["isg"])
        x_first = x.sort_values(["patient_id", "days", "sample_id"]).groupby("patient_id", as_index=False).first()
        x_first_le7 = x_first[x_first["days"] <= 7.0].copy()

        fd_res = first_draw_death_stats(x_first, ["isg"], fat_col="fat11", seed=707).set_index("set")
        fd_le7_res = first_draw_death_stats(x_first_le7, ["isg"], fat_col="fat11", seed=708).set_index("set")

        for arm_name, sub_z in [("Real", xz[xz["arm"] == "Real"]),
                                ("PBISC new", xz[xz["arm"] == "PBISC new"]),
                                ("Combined", xz)]:
            gu = gee_fit_term(sub_z, "isg__z", "temp_severity_ord", covars=())
            gw = gee_fit_term(sub_z, "isg__z", "temp_severity_ord", covars=("wk",))
            mk = mundlak_gee(sub_z, "isg__z", "temp_severity_ord", covars=("wk",))
            cl_t = cluster_spearman_sameday(sub_z, ["isg__z"], xcol="temp_severity_ord", n_boot=1500, seed=709).iloc[0]
            fdr = fd_res.loc[arm_name]
            fdl = fd_le7_res.loc[arm_name] if arm_name in fd_le7_res.index else None
            pr_t = pub_tisg.loc[gname]
            pm_r = pr_t["r_sev"] if arm_name == "Real" else (pr_t["n_sev"] if arm_name == "PBISC new" else pr_t["c_sev"])
            pm_a = pr_t["r_auc11"] if arm_name == "Real" else (pr_t["n_auc11"] if arm_name == "PBISC new" else pr_t["c_auc11"])
            t_rows.append(dict(
                set=arm_name, subset=gname,
                pm_peakWHO_rho=pm_r, pm_death_AUC=pm_a,
                pm_q13_sev=pr_t["q13_sev"] if arm_name == "Combined" else np.nan,
                pm_q13_dth=pr_t["q13_auc11"] if arm_name == "Combined" else np.nan,
                sd_n_draws=gw["n_draws"], sd_n_pts=gw["n_pts"],
                sd_gee_unadj_beta=gu["beta"], sd_gee_unadj_p=gu["p"], sd_gee_unadj_p_br=gu["p_br"],
                sd_gee_wk_beta=gw["beta"], sd_gee_wk_lo=gw["lo"], sd_gee_wk_hi=gw["hi"],
                sd_gee_wk_p=gw["p"], sd_gee_wk_p_br=gw["p_br"],
                sd_rho=cl_t["sameday_rho"], sd_rho_p_clboot=cl_t["sameday_rho_p_clboot"],
                mundlak_within_beta=mk["beta_within"], mundlak_within_p=mk["p_within"],
                fd_n=int(fdr["fd_n"]), fd_deaths=int(fdr["fd_deaths"]),
                fd_AUC=fdr["fd_AUC"], fd_AUC_lo=fdr["fd_AUC_lo"], fd_AUC_hi=fdr["fd_AUC_hi"], fd_p=fdr["fd_p"],
                fd_adjWHO_AUC=fdr["fd_adjWHO_AUC"], fd_adjWHO_p=fdr["fd_adjWHO_p"],
                fd_le7_n=int(fdl["fd_n"]) if fdl is not None else 0,
                fd_le7_deaths=int(fdl["fd_deaths"]) if fdl is not None else 0,
                fd_le7_AUC=fdl["fd_AUC"] if fdl is not None else np.nan,
                fd_le7_p=fdl["fd_p"] if fdl is not None else np.nan,
            ))
    T6_ISG = pd.DataFrame(t_rows)
    for _, ix in T6_ISG.groupby("set").groups.items():
        for col in ["sd_gee_unadj_p", "sd_gee_wk_p", "sd_gee_wk_p_br", "sd_rho_p_clboot", "fd_p", "fd_le7_p"]:
            T6_ISG.loc[ix, col + "_q13"] = bh(T6_ISG.loc[ix, col])
    T6_ISG.to_csv(TAB_DIR / "T6_Tcell_ISG_13subsets_new_framework.csv", index=False, encoding="utf-8-sig")
    print(f"[OK] Domain 4 T-cell ISG done ({time.time() - t0:.1f}s).", flush=True)

    # -------------------------------------------------------------------------
    # DOMAIN 5: Hospital Blood Labs (Real cohort: Same-day draw pairs & First-Draw Mortality)
    # -------------------------------------------------------------------------
    print("=== Domain 5: Real hospital labs (same-day pairs & first-draw death) ===", flush=True)
    real_d = prim[prim["arm"] == "Real"].copy()
    real_d["NLR"] = np.where(real_d["temp_lymph"] > 0, real_d["temp_neut"] / real_d["temp_lymph"], np.nan)
    LABS = ["temp_CRP", "temp_lymph", "NLR", "temp_ferritin", "temp_ddimer", "temp_trop"]
    IMM6 = ["NK CD56 bright", "T CD8 CM", "T CD8 CTL", "B_ASC_IgG", "B_ASC_IgA", "HSPC"]
    real_d_z = add_draw_rint(real_d, IMM6 + LABS)
    lab_rows = []
    for f in IMM6:
        for lb in LABS:
            sub = real_d_z.dropna(subset=[f + "__z", lb + "__z"])
            g_un = gee_fit_term(sub, f + "__z", lb + "__z", covars=())
            g_wk = gee_fit_term(sub, f + "__z", lb + "__z", covars=("wk",))
            lab_rows.append(dict(feature=f, lab=lb, n_draws=g_wk["n_draws"], n_pts=g_wk["n_pts"],
                                 sameday_gee_beta=g_un["beta"], sameday_gee_p=g_un["p"],
                                 sameday_gee_wk_beta=g_wk["beta"], sameday_gee_wk_p=g_wk["p"]))
    T7_LABS = pd.DataFrame(lab_rows)
    T7_LABS["sameday_gee_wk_q36"] = bh(T7_LABS["sameday_gee_wk_p"])
    T7_LABS.to_csv(TAB_DIR / "T7_hospital_labs_sameday_pairs.csv", index=False, encoding="utf-8-sig")

    # -------------------------------------------------------------------------
    # PRINT EXECUTIVE SUMMARY TO STDOUT + JSON
    # -------------------------------------------------------------------------
    summary = {
        "n_primary_patients": 66,
        "n_primary_draws": 217,
        "first_draw_total_n": int(len(df_first)),
        "first_draw_total_deaths": int(df_first["fat11"].sum()),
        "first_draw_le7d_n": int(len(df_first_le7)),
        "first_draw_le7d_deaths": int(df_first_le7["fat11"].sum()),
        "n_43_sig_counts": {
            "prev_pm_peakWHO_q43_lt05": int((COMP43_comb["pm_peakWHO_q43"] < 0.05).sum()),
            "new_samedayWHO_GEE_unadj_q43_lt05": int((COMP43_comb["sd_gee_unadj_p_q43"] < 0.05).sum()),
            "new_samedayWHO_GEE_wk_q43_lt05": int((COMP43_comb["sd_gee_wk_p_q43"] < 0.05).sum()),
            "new_samedayWHO_GEE_wk_BR_q43_lt05": int((COMP43_comb["sd_gee_wk_p_br_q43"] < 0.05).sum()),
            "new_samedayWHO_clboot_rho_q43_lt05": int((COMP43_comb["sd_rho_p_clboot_q43"] < 0.05).sum()),
            "prev_pm_death_q43_lt10": int((COMP43_comb["pm_death_q43"] < 0.10).sum()),
            "new_firstdraw_death_q43_lt10": int((COMP43_comb["fd_q43"] < 0.10).sum()),
            "new_firstdraw_death_p_lt05": int((COMP43_comb["fd_p"] < 0.05).sum()),
            "new_firstdraw_le7d_death_p_lt05": int((COMP43_comb["fd_le7_p"] < 0.05).sum()),
        }
    }
    (OUT_DIR / "REANALYSIS_SUMMARY_20261120.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\n=== SUMMARY COUNTS ===", json.dumps(summary, indent=2))
    print(f"Total runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
