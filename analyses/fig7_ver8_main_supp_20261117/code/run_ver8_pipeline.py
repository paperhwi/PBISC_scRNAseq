"""Complete reproducible statistical pipeline for Harmonized Figure 7 (ver8) + Supplementary Figures S1 & S2.

Integrates:
  1. Bulk PBMC ROR1 / FMOD molecular screen across all 88 bulk donors (57 Real/Shared, 29 PBISC-only, 2 Benign CLL: C116, C139).
  2. Clean 43-fine-label screen & pre-specified 5-step selection funnel (43 -> 16 -> 15 -> 12 -> 8).
  3. Credibility scorecard for all 8 PBISC-dependent candidates (cross-source homogeneity P, 27-shared-donor rho, 14-variant pass rate, WHO 2-3 death AUC).
  4. Bootstrap 95% CIs (with and without CLL) across Concordant Anchors & 5 Focus Populations.
  5. HSPC mortality robustness across 8 definitions / strata (All deaths, Table S1, Direct COVID-19, Within WHO 2-3, Severity-stratified, Cambridge only, >day 7, PBISC new alone).
  6. Subsampling power curves (n = 10..35 Real patients) + n=20 power + n80 sample-size reduction + shuffled-outcome controls.
  7. Modality baseline scale audit on the 27 shared donors (Real CITE-seq % PBMC vs PBISC bulk % PBMC) + raw % PBMC vs within-source Normal score comparison.
  8. 14-variant confounder, denominator, and dataset sensitivity matrix + NK compartment decomposition.
"""
import pathlib
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(r"F:\PBISC")
RES_DIR = ROOT / "eskd_covid_clinical_assoc_20261115" / "results"
RECONC_PATH = ROOT / "paper_claim_analysis" / "covid_sample_patient_id_reconciliation.csv"
CLL_BULK_PATH = (
    ROOT
    / "HSPC_CD56brightNK_COVID_analysis_summary_20261116"
    / "12_CLL_excluded_severity_death_discovery"
    / "data"
    / "cll_molecular_verification_bulk.csv"
)
V8_DIR = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116" / "15_Figure7_ver8_Main_and_Supp"
DATA_DIR = V8_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

CLL = {"C116", "C139"}
EXTRA_OUTLIERS = {"C187", "C190"}

# -------------------------------------------------------------------------
# 0. Copy / Verify Bulk PBMC ROR1 / FMOD CLL Molecular Screen (105 total bulk; 68 COVID+)
# -------------------------------------------------------------------------
cll_bulk_all = pd.read_csv(CLL_BULK_PATH)
valid_groups = {"Real / Shared donor", "PBISC-only donor", "Excluded CLL (C116, C139)"}
cll_bulk = cll_bulk_all[cll_bulk_all["cohort_group"].isin(valid_groups)].copy()
S_tmp = pd.read_csv(RES_DIR / "sample_level_composition_metadata.csv")
covid_pos_ids = set(
    S_tmp[(S_tmp["condition"] == "positive") & S_tmp["in_primary_scope"].astype(bool)]["patient_id"]
)
cll_bulk["covid_positive_scope"] = cll_bulk["patient_id"].isin(covid_pos_ids)
cll_bulk.to_csv(DATA_DIR / "cll_molecular_screen_bulk_donors.csv", index=False)

# -------------------------------------------------------------------------
# 1. Load Sample Metadata & Fine-Label Compositions (% of Leukocytes / PBMCs)
# -------------------------------------------------------------------------
S = pd.read_csv(RES_DIR / "sample_level_composition_metadata.csv")
Fm_raw = pd.read_csv(RES_DIR / "sample_level_fine_label_pct.csv")
reconc = pd.read_csv(RECONC_PATH)

ster_map = reconc.set_index("sample_id")["published_received_steroid"].eq("Y").to_dict()
S["steroid"] = S.apply(
    lambda r: bool(ster_map.get(str(r["sample_id"]), False)) if r["source"] == "Real" else np.nan,
    axis=1,
)
S["batch"] = S["centre_or_wave"]

LAB = [c for c in Fm_raw.columns if c not in ("source", "sample_id", "RBC", "Platelet")]
Fm = Fm_raw.copy()
Fm[LAB] = Fm[LAB].apply(pd.to_numeric, errors="coerce").fillna(0.0)
Fm[LAB] = 100.0 * Fm[LAB].div(Fm[LAB].sum(axis=1), axis=0)

LIN = {
    "B": [c for c in LAB if c.startswith("B_")],
    "NK": ["NK", "NK CD56 bright", "NK IFN stim"],
    "Mono": ["CD14mono", "CD14mono_IFN", "CD14mono_dex", "CD16mono", "CD16mono_C1", "CD16mono_IFN", "Int.mono"],
    "T": [c for c in LAB if c.startswith("T ") or c in ("MAIT", "NKT")],
}
MYELOID_ALL = LIN["Mono"] + ["DC1", "DC2", "DC3", "DC3_IFN", "ASDC", "pDC"]

AUDIT_FEATS = {
    "NK CD56 bright": "NK",
    "B_ASC_dividing": "B",
    "T CD8 CM": "T",
    "CD16mono": "Mono",
    "B_ASC_IgM": "B",
    "B_ASC_IgG": "B",
    "B_ASC_IgA": "B",
    "T CD8 CTL": "T",
    "HSPC": None,
    "CD14mono_IFN": "Mono",
    "NK": "NK",
    "T g/d": "T",
    "CD16mono_C1": "Mono",
}
for f, l in AUDIT_FEATS.items():
    Fm[f + "|lin"] = (100.0 * Fm[f] / Fm[LIN[l]].sum(axis=1).replace(0, np.nan)) if l else np.nan
    Fm[f + "|nomyeloid"] = (100.0 * Fm[f] / (100.0 - Fm[MYELOID_ALL].sum(axis=1)).replace(0, np.nan)) if l != "Mono" else np.nan

BROAD_COLS = ["T_pct", "B_pct", "Mono_pct", "NK_pct", "Other_pct"]


def rint(x: pd.Series) -> pd.Series:
    r = x.rank()
    return pd.Series(stats.norm.ppf((r - 0.5) / x.notna().sum()), index=x.index)


def auc_mw(y, f):
    y_a = np.asarray(y, dtype=float)
    f_a = np.asarray(f, dtype=float)
    m = np.isfinite(y_a) & np.isfinite(f_a)
    y_a, f_a = y_a[m], f_a[m]
    n1 = int((f_a == 1).sum())
    n0 = int((f_a == 0).sum())
    if n1 < 1 or n0 < 1:
        return np.nan, np.nan
    u = stats.mannwhitneyu(y_a[f_a == 1], y_a[f_a == 0])
    return u.statistic / (n1 * n0), u.pvalue


def resid(y, X):
    X_mat = np.column_stack([np.ones(len(y))] + [X]) if (X is not None and np.size(X)) else np.ones((len(y), 1))
    b = np.linalg.lstsq(X_mat, y, rcond=None)[0]
    return y - X_mat @ b


def ptab(src, excl=CLL, first=False, maxday=None, steroid_naive=False, mono_ok=False):
    d = S[
        (S["source"] == src)
        & (S["condition"] == "positive")
        & S["in_primary_scope"].astype(bool)
        & S["qc_pass_min200cells"].astype(bool)
    ].copy()
    if excl:
        d = d[~d["patient_id"].isin(excl)]
    if maxday is not None:
        d = d[d["days_from_first_positive_swab"] <= maxday]
    if steroid_naive and src == "Real":
        d = d[d["steroid"] == False]
    if mono_ok and src == "Real":
        d = d[d["Mono_pct"] >= 2.0]
    if first:
        d = d.sort_values("days_from_first_positive_swab").groupby("patient_id").head(1)
    d = d.merge(Fm[Fm["source"] == src].drop(columns="source"), on="sample_id", how="left")
    fm_num_cols = [c for c in Fm.columns if c not in ("source", "sample_id")]
    g = d.groupby("patient_id")
    T = g[BROAD_COLS + fm_num_cols].mean()
    T["sev"] = g["severity_ord"].first()
    T["fatal"] = g["fatal_harmonized"].first()
    T["fatal_direct"] = g["fatal_direct_only"].first()
    T["age"] = g["age"].first()
    T["batch"] = g["batch"].first()
    T["male"] = g["sex"].first().eq("M").astype(float)
    T["steroid_any"] = g["steroid"].max().astype(float) if src == "Real" else 0.0
    T["src"] = src
    return T


def evaluate_variant(R_v, N_v, fcol, cov=(), death_col="fatal"):
    out_v = {}
    for lab, parts in [("Real", [R_v]), ("Real + PBISC new", [R_v, N_v])]:
        D = pd.concat([
            pd.DataFrame({"y": rint(p[fcol]), "sev": p["sev"], "fatal": p[death_col], "src": p["src"], **{c: p[c] for c in cov}})
            for p in parts
        ]).dropna(subset=["y", "sev"] + list(cov))
        cov_cols = list(cov) + (["src"] if len(parts) > 1 else [])
        Xc = pd.get_dummies(D[cov_cols], drop_first=True).astype(float).values if (cov or len(parts) > 1) else None
        if Xc is not None and Xc.shape[1] == 0:
            Xc = None
        if len(cov) == 0:
            # Unadjusted baseline uses rank-INT pooled Spearman / Mann-Whitney directly
            r, p_s = stats.spearmanr(D["sev"], D["y"])
            auc, p_d = auc_mw(D["y"], D["fatal"])
        else:
            ry = resid(D["y"].values, Xc)
            rs = resid(D["sev"].values.astype(float), Xc)
            r, p_s = stats.pearsonr(rs, ry)
            auc, p_d = auc_mw(ry, D["fatal"].values)
        out_v[lab] = dict(n=len(D), deaths=int(D["fatal"].dropna().sum()), rho=r, p_sev=p_s, AUC=auc, p_death=p_d)
    return out_v


# -------------------------------------------------------------------------
# 2. Clean 43-Fine-Label Scan & Pre-Specified Selection Funnel
# -------------------------------------------------------------------------
scan_rows = []
for drop_cll in [False, True]:
    R_cur = ptab("Real", excl=CLL if drop_cll else set())
    B_cur = ptab("PBISC", excl=CLL if drop_cll else set())
    N_cur = B_cur.loc[~B_cur.index.isin(R_cur.index)]

    fine_keep = [f for f in LAB if R_cur[f].gt(0).mean() >= 0.30]
    rng_h = np.random.default_rng(12345)

    for f in fine_keep:
        if R_cur[f].nunique() >= 4 and N_cur[f].nunique() >= 4:
            rr, _ = stats.spearmanr(R_cur["sev"], R_cur[f])
            rn, _ = stats.spearmanr(N_cur["sev"], N_cur[f])
            z = (np.arctanh(rr) - np.arctanh(rn)) / np.sqrt(1.0 / (len(R_cur) - 3) + 1.0 / (len(N_cur) - 3))
            p_het_s = 2.0 * stats.norm.sf(abs(z))
            ur, _ = auc_mw(R_cur[f], R_cur["fatal"])
            un, _ = auc_mw(N_cur[f], N_cur["fatal"])
            r1, r0 = R_cur[f][R_cur["fatal"] == 1].values, R_cur[f][R_cur["fatal"] == 0].values
            n1, n0 = N_cur[f][N_cur["fatal"] == 1].values, N_cur[f][N_cur["fatal"] == 0].values
            diffs_b = [
                auc_mw(rng_h.choice(r1, len(r1), True), np.ones(len(r1)))[0]  # dummy unused
                for _ in range(0)
            ]
            bs_d = []
            for _ in range(2000):
                ar = stats.mannwhitneyu(rng_h.choice(r1, len(r1), True), rng_h.choice(r0, len(r0), True)).statistic / (len(r1) * len(r0))
                an = stats.mannwhitneyu(rng_h.choice(n1, len(n1), True), rng_h.choice(n0, len(n0), True)).statistic / (len(n1) * len(n0))
                bs_d.append(ar - an)
            p_het_d = 2.0 * stats.norm.sf(abs(ur - un) / max(np.std(bs_d, ddof=1), 1e-6))
        else:
            p_het_s, p_het_d = np.nan, np.nan

        sets = {
            "Real": pd.DataFrame({"y": R_cur[f], "sev": R_cur["sev"], "fatal": R_cur["fatal"]}),
            "PBISC new": pd.DataFrame({"y": N_cur[f], "sev": N_cur["sev"], "fatal": N_cur["fatal"]}),
            "Real + PBISC new": pd.concat([
                pd.DataFrame({"y": rint(R_cur[f]), "sev": R_cur["sev"], "fatal": R_cur["fatal"]}),
                pd.DataFrame({"y": rint(N_cur[f]), "sev": N_cur["sev"], "fatal": N_cur["fatal"]}),
            ]),
        }
        for sname, df_s in sets.items():
            df_s = df_s.dropna(subset=["y", "sev"])
            if df_s["y"].nunique() < 4:
                continue
            r, p_s = stats.spearmanr(df_s["sev"], df_s["y"])
            auc_d, p_d = auc_mw(df_s["y"], df_s["fatal"])
            scan_rows.append({
                "CLL_excluded": drop_cll,
                "set": sname,
                "feature": f,
                "n": len(df_s),
                "deaths": int(df_s["fatal"].dropna().sum()),
                "rho": r,
                "p_sev": p_s,
                "AUC_death": auc_d,
                "p_death": p_d,
                "p_het_sev": p_het_s,
                "p_het_death": p_het_d,
            })

df_scan = pd.DataFrame(scan_rows)
for pcol, qcol in [("p_sev", "q_sev"), ("p_death", "q_death")]:
    df_scan[qcol] = np.nan
    for _, idx in df_scan.groupby(["CLL_excluded", "set"]).groups.items():
        df_scan.loc[idx, qcol] = multipletests(df_scan.loc[idx, pcol].values, method="fdr_bh")[1]
df_scan.to_csv(DATA_DIR / "cll_scan_43fine_v8.csv", index=False)

# Build Selection Funnel Table on CLL-excluded 43 fine labels
R = ptab("Real", excl=CLL)
B = ptab("PBISC", excl=CLL)
NEW = B.loc[~B.index.isin(R.index)]
SHARED = R.index.intersection(B.index)

sc_no = df_scan[df_scan["CLL_excluded"] == True]
W = sc_no.pivot_table(index="feature", columns="set", values=["rho", "AUC_death", "q_sev", "q_death", "p_sev", "p_death", "p_het_sev", "p_het_death"])
W.columns = [f"{a}|{b}" for a, b in W.columns]
W = W.reset_index()

sel_rows = []
for _, r in W.iterrows():
    f = r["feature"]
    E = "Real + PBISC new"
    s2_sev = bool(r[f"q_sev|{E}"] < 0.05)
    s2_dth = bool(r[f"q_death|{E}"] < 0.10)
    rho_r, rho_n, rho_e = r["rho|Real"], r.get("rho|PBISC new", np.nan), r[f"rho|{E}"]
    auc_r, auc_n, auc_e = r["AUC_death|Real"], r.get("AUC_death|PBISC new", np.nan), r[f"AUC_death|{E}"]

    s3_sev = bool(
        s2_sev
        and pd.notna(rho_n)
        and (np.sign(rho_r) == np.sign(rho_n) == np.sign(rho_e))
        and (min(abs(rho_r), abs(rho_n)) >= 0.10)
    )
    s3_dth = bool(
        s2_dth
        and pd.notna(auc_n)
        and (np.sign(auc_r - 0.5) == np.sign(auc_n - 0.5) == np.sign(auc_e - 0.5))
        and (min(abs(auc_r - 0.5), abs(auc_n - 0.5)) >= 0.05)
    )
    if s3_sev or s3_dth:
        adj_res = evaluate_variant(R, NEW, f, cov=("batch",))["Real + PBISC new"]
        ps_c, pd_c = adj_res["p_sev"], adj_res["p_death"]
    else:
        ps_c, pd_c = np.nan, np.nan

    not_dex = f != "CD14mono_dex"
    s4_sev = bool(s3_sev and not_dex and (ps_c < 0.05))
    s4_dth = bool(s3_dth and not_dex and (pd_c < 0.05))
    s5_sev = bool(s4_sev and (r["q_sev|Real"] >= 0.05))
    s5_dth = bool(s4_dth and (r["q_death|Real"] >= 0.10))
    anchor_sev = bool(s4_sev and (r["q_sev|Real"] < 0.05))

    sel_rows.append({
        "feature": f,
        "rho_Real": rho_r,
        "rho_New": rho_n,
        "rho_Ext": rho_e,
        "p_sev_Real": r["p_sev|Real"],
        "p_sev_New": r.get("p_sev|PBISC new", np.nan),
        "p_sev_Ext": r[f"p_sev|{E}"],
        "q_sev_Real": r["q_sev|Real"],
        "q_sev_Ext": r[f"q_sev|{E}"],
        "AUC_Real": auc_r,
        "AUC_New": auc_n,
        "AUC_Ext": auc_e,
        "p_death_Real": r["p_death|Real"],
        "p_death_New": r.get("p_death|PBISC new", np.nan),
        "p_death_Ext": r[f"p_death|{E}"],
        "q_death_Real": r["q_death|Real"],
        "q_death_Ext": r[f"q_death|{E}"],
        "p_het_sev": r[f"p_het_sev|{E}"],
        "p_het_death": r[f"p_het_death|{E}"],
        "p_sev_centre_adj": ps_c,
        "p_death_centre_adj": pd_c,
        "S2_sev": s2_sev,
        "S2_death": s2_dth,
        "S3_sev": s3_sev,
        "S3_death": s3_dth,
        "S4_sev": s4_sev,
        "S4_death": s4_dth,
        "S5_sev": s5_sev,
        "S5_death": s5_dth,
        "Anchor_sev": anchor_sev,
    })

df_sel = pd.DataFrame(sel_rows)
df_sel.to_csv(DATA_DIR / "selection_funnel_43fine.csv", index=False)

# -------------------------------------------------------------------------
# 3. 14-Variant Sensitivity Matrix across Anchors + All S5 Candidates + CD14mono_IFN
# -------------------------------------------------------------------------
VARIANTS = [
    ("0 baseline (no CLL)", {}, {}, "pbmc"),
    ("A1 adjust centre/wave", {}, {"cov": ("batch",)}, "pbmc"),
    ("A2 within-lineage %", {}, {}, "lin"),
    ("A3 % of non-myeloid cells", {}, {}, "nomyeloid"),
    ("A4 drop Real samples with mono <2%", {"mono_ok": True}, {}, "pbmc"),
    ("B  also drop C187, C190", {"excl": CLL | EXTRA_OUTLIERS}, {}, "pbmc"),
    ("C1 death = Table S1 (C101, C147 alive)", {}, {"death_col": "fatal_s1"}, "pbmc"),
    ("C2 death = direct only", {}, {"death_col": "fatal_direct"}, "pbmc"),
    ("D1 Real steroid-naive samples only", {"steroid_naive": True}, {}, "pbmc"),
    ("D2 adjust patient CD14mono_dex %", {}, {"cov": ("CD14mono_dex",)}, "pbmc"),
    ("E  samples <= day 20 (both sources)", {"maxday": 20}, {}, "pbmc"),
    ("F1 first sample per patient", {"first": True}, {}, "pbmc"),
    ("F2 adjust age + sex", {}, {"cov": ("age", "male")}, "pbmc"),
    ("G  combined: lineage %, no C187/C190, adj centre/wave+age+sex", {"excl": CLL | EXTRA_OUTLIERS}, {"cov": ("batch", "age", "male")}, "lin"),
]

base_real_ids = R.index
sa_rows = []
for vname, pk, ek, den in VARIANTS:
    Rv = ptab("Real", **pk)
    Bv = ptab("PBISC", **pk)
    Nv = Bv.loc[~Bv.index.isin(base_real_ids)]
    for T_part in (Rv, Nv):
        T_part["fatal_s1"] = T_part["fatal"].copy()
        T_part.loc[T_part.index.isin(["C101", "C147"]), "fatal_s1"] = 0.0
    for f, l in AUDIT_FEATS.items():
        col = f if den == "pbmc" else f"{f}|{den}"
        if den != "pbmc" and Rv[col].isna().all():
            if den == "lin" and l is None:
                col = f + "|nomyeloid"
            else:
                continue
        for lab, res_d in evaluate_variant(Rv, Nv, col, **ek).items():
            sa_rows.append(dict(variant=vname, feature=f, set=lab, **res_d))

df_sa = pd.DataFrame(sa_rows)
df_sa.to_csv(DATA_DIR / "sensitivity_14variants_all_v8.csv", index=False)

# -------------------------------------------------------------------------
# 4. Credibility Scorecard for All 8 S5 Candidates (10 Population x Outcome rows)
# -------------------------------------------------------------------------
SCORECARD_SPECS = [
    ("B_ASC_IgG", "sev", True),
    ("B_ASC_IgA", "sev", True),
    ("T CD8 CTL", "sev", True),
    ("T CD8 CM", "death", True),
    ("T CD8 CTL", "death", True),
    ("HSPC", "death", True),
    ("NK", "sev", False),
    ("NK", "death", False),
    ("T g/d", "death", False),
    ("CD16mono_C1", "death", False),
]
rng_sc = np.random.default_rng(11)
sc_rows = []
for f, out_type, followed in SCORECARD_SPECS:
    rs_sh, ps_sh = stats.spearmanr(R.loc[SHARED, f], B.loc[SHARED, f])
    E_df = pd.concat([
        pd.DataFrame({"y": rint(R[f]), "sev": R["sev"], "fatal": R["fatal"]}),
        pd.DataFrame({"y": rint(NEW[f]), "sev": NEW["sev"], "fatal": NEW["fatal"]}),
    ])
    sub23 = E_df[E_df["sev"] >= 2]
    a23, p23 = auc_mw(sub23["y"], sub23["fatal"])

    if out_type == "sev":
        eff_r, _ = stats.spearmanr(R["sev"], R[f])
        eff_n, _ = stats.spearmanr(NEW["sev"], NEW[f])
        z_h = (np.arctanh(eff_r) - np.arctanh(eff_n)) / np.sqrt(1.0 / (len(R) - 3) + 1.0 / (len(NEW) - 3))
        p_het = 2.0 * stats.norm.sf(abs(z_h))
        vcol, pcol, ref = "rho", "p_sev", 0.0
    else:
        eff_r, _ = auc_mw(R[f], R["fatal"])
        eff_n, _ = auc_mw(NEW[f], NEW["fatal"])
        # Bootstrap SE of AUC in each cohort
        se_list = []
        for T_sub in (R, NEW):
            vals_b = []
            for _ in range(3000):
                sb = T_sub.sample(len(T_sub), replace=True, random_state=rng_sc.integers(1_000_000_000))
                if 0 < sb["fatal"].sum() < len(sb):
                    vals_b.append(auc_mw(sb[f], sb["fatal"])[0])
            se_list.append(np.std(vals_b))
        z_h = (eff_r - eff_n) / np.sqrt(se_list[0] ** 2 + se_list[1] ** 2)
        p_het = 2.0 * stats.norm.sf(abs(z_h))
        vcol, pcol, ref = "AUC", "p_death", 0.5

    sa_sub = df_sa[(df_sa["feature"] == f) & (df_sa["set"] == "Real + PBISC new")]
    s0 = np.sign(sa_sub[sa_sub["variant"].str.startswith("0")][vcol].iloc[0] - ref)
    n_tot_var = int(sa_sub[vcol].notna().sum())
    n_rob = int(((sa_sub[pcol] < 0.05) & (np.sign(sa_sub[vcol] - ref) == s0)).sum())

    sc_rows.append({
        "feature": f,
        "outcome": out_type,
        "followed": followed,
        "effect_Real": eff_r,
        "effect_New": eff_n,
        "p_het": p_het,
        "shared_rho": rs_sh,
        "shared_p": ps_sh,
        "n_shared": len(SHARED),
        "robust_pass": n_rob,
        "robust_total": n_tot_var,
        "AUC_within_WHO23": a23 if out_type == "death" else np.nan,
        "p_within_WHO23": p23 if out_type == "death" else np.nan,
    })

df_sc = pd.DataFrame(sc_rows)
df_sc.to_csv(DATA_DIR / "credibility_scorecard_v8.csv", index=False)

# -------------------------------------------------------------------------
# 5. Bootstrap 95% CIs (With vs Without CLL) across Anchors & Focus Populations
# -------------------------------------------------------------------------
CI_FEATS = ["NK CD56 bright", "B_ASC_dividing", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CM", "T CD8 CTL", "HSPC", "CD14mono_IFN"]
rng_ci = np.random.default_rng(20261117)
ci_rows = []
for drop_cll in [False, True]:
    Rc = ptab("Real", excl=CLL if drop_cll else set())
    Bc = ptab("PBISC", excl=CLL if drop_cll else set())
    Nc = Bc.loc[~Bc.index.isin(Rc.index)]
    for f in CI_FEATS:
        Ec = pd.concat([
            pd.DataFrame({"y": rint(Rc[f]), "sev": Rc["sev"], "fatal": Rc["fatal"]}),
            pd.DataFrame({"y": rint(Nc[f]), "sev": Nc["sev"], "fatal": Nc["fatal"]}),
        ]).dropna()
        for sname, d in [
            ("Real", pd.DataFrame({"y": Rc[f], "sev": Rc["sev"], "fatal": Rc["fatal"]}).dropna()),
            ("PBISC new", pd.DataFrame({"y": Nc[f], "sev": Nc["sev"], "fatal": Nc["fatal"]}).dropna()),
            ("Real + PBISC new", Ec),
        ]:
            r, ps = stats.spearmanr(d["sev"], d["y"])
            auc, pd_ = auc_mw(d["y"], d["fatal"])
            bs_r, bs_a = [], []
            y_arr, s_arr, f_arr = d["y"].values, d["sev"].values, d["fatal"].values
            n_d = len(d)
            for _ in range(2000):
                ix = rng_ci.integers(0, n_d, n_d)
                bs_r.append(stats.spearmanr(s_arr[ix], y_arr[ix])[0])
                fb = f_arr[ix]
                if 0 < fb.sum() < n_d:
                    bs_a.append(auc_mw(y_arr[ix], fb)[0])
            ci_rows.append({
                "CLL_excluded": drop_cll,
                "feature": f,
                "set": sname,
                "n": n_d,
                "deaths": int(f_arr.sum()),
                "rho": r,
                "rho_lo": np.nanpercentile(bs_r, 2.5),
                "rho_hi": np.nanpercentile(bs_r, 97.5),
                "p_sev": ps,
                "AUC": auc,
                "AUC_lo": np.percentile(bs_a, 2.5),
                "AUC_hi": np.percentile(bs_a, 97.5),
                "p_death": pd_,
            })

df_ci = pd.DataFrame(ci_rows)
df_ci.to_csv(DATA_DIR / "cohort_effect_sizes_ci_v8.csv", index=False)

# -------------------------------------------------------------------------
# 6. HSPC Death Robustness Across 8 Mortality Definitions / Strata
# -------------------------------------------------------------------------
rng_hs = np.random.default_rng(2026)


def strat_auc(D):
    num = den = 0.0
    for _, g in D.groupby("sev"):
        if 0 < g["fatal"].sum() < len(g):
            num += stats.mannwhitneyu(g["y"][g["fatal"] == 1], g["y"][g["fatal"] == 0]).statistic
            den += g["fatal"].sum() * (len(g) - g["fatal"].sum())
    return num / den if den else np.nan


def strat_perm(D, nperm=4000):
    obs = strat_auc(D)
    idx = D.groupby("sev").indices
    c = 0
    fat0 = D["fatal"].values
    for _ in range(nperm):
        fat = fat0.copy()
        for _, ii in idx.items():
            fat[ii] = rng_hs.permutation(fat[ii])
        c += abs(strat_auc(D.assign(fatal=fat)) - 0.5) >= abs(obs - 0.5)
    return obs, (c + 1) / (nperm + 1)


def late_hspc(src):
    d = S[
        (S["source"] == src)
        & (S["condition"] == "positive")
        & S["in_primary_scope"].astype(bool)
        & S["qc_pass_min200cells"].astype(bool)
        & (~S["patient_id"].isin(CLL))
        & (S["days_from_first_positive_swab"] > 7)
    ].merge(Fm[Fm["source"] == src].drop(columns="source"), on="sample_id")
    g = d.groupby("patient_id")
    T = pd.DataFrame({
        "HSPC": g["HSPC"].mean(),
        "sev": g["severity_ord"].first(),
        "fatal": g["fatal_harmonized"].first(),
        "batch": g["batch"].first(),
    })
    return T.loc[~T.index.isin(R.index)] if src == "PBISC" else T


R_h = R.copy()
N_h = NEW.copy()
for Tp in (R_h, N_h):
    Tp["fatal_s1"] = Tp["fatal"].copy()
    Tp.loc[Tp.index.isin(["C101", "C147"]), "fatal_s1"] = 0.0

hs_rows = []
for lab, parts in [("Real", [("Real", R_h)]), ("Real + PBISC new", [("Real", R_h), ("New", N_h)])]:
    D = pd.concat([pd.DataFrame({"y": rint(p["HSPC"]), "sev": p["sev"], "fatal": p["fatal"], "batch": p["batch"], "src": nm}) for nm, p in parts])
    Ds = pd.concat([pd.DataFrame({"y": rint(p["HSPC"]), "sev": p["sev"], "fatal": p["fatal_s1"], "batch": p["batch"], "src": nm}) for nm, p in parts])
    Dd = pd.concat([pd.DataFrame({"y": rint(p["HSPC"]), "sev": p["sev"], "fatal": p["fatal_direct"], "batch": p["batch"], "src": nm}) for nm, p in parts])
    sub23 = D[D["sev"] >= 2]
    so, sp = strat_perm(D)
    nc = D[~((D["src"] == "Real") & (D["batch"] == "NCL"))]
    Lt = pd.concat([pd.DataFrame({"y": rint(late_hspc("Real" if nm == "Real" else "PBISC")["HSPC"]),
                                 "fatal": late_hspc("Real" if nm == "Real" else "PBISC")["fatal"]}) for nm, _ in parts])
    for an_name, (a_val, p_val), nd in [
        ("All deaths", auc_mw(D["y"], D["fatal"]), D["fatal"].sum()),
        ("Table S1 deaths", auc_mw(Ds["y"], Ds["fatal"]), Ds["fatal"].sum()),
        ("Direct COVID-19 deaths", auc_mw(Dd["y"], Dd["fatal"]), Dd["fatal"].sum()),
        ("Within WHO 2–3", auc_mw(sub23["y"], sub23["fatal"]), sub23["fatal"].sum()),
        ("Severity-stratified", (so, sp), D["fatal"].sum()),
        ("Cambridge only (Real)", auc_mw(nc["y"], nc["fatal"]), nc["fatal"].sum()),
        ("Samples > day 7", auc_mw(Lt["y"], Lt["fatal"]), Lt["fatal"].sum()),
    ]:
        hs_rows.append(dict(analysis=an_name, set=lab, deaths=int(nd), AUC=a_val, p=p_val))

a_new, p_new = auc_mw(rint(N_h["HSPC"]), N_h["fatal"])
hs_rows.append(dict(analysis="PBISC new alone", set="PBISC new", deaths=int(N_h["fatal"].sum()), AUC=a_new, p=p_new))
df_hs = pd.DataFrame(hs_rows)
df_hs.to_csv(DATA_DIR / "hspc_death_analyses_v8.csv", index=False)

# -------------------------------------------------------------------------
# 7. Subsampling Power Curves (n = 10..35) + Shuffled Null Controls
# -------------------------------------------------------------------------
POWER_SPECS = [
    ("NK CD56 bright", "sev", -1),
    ("B_ASC_dividing", "sev", +1),
    ("T CD8 CM", "sev", -1),
    ("CD16mono", "sev", -1),
    ("CD16mono|lin", "sev", -1),
    ("B_ASC_IgG", "sev", +1),
    ("B_ASC_IgA", "sev", +1),
    ("T CD8 CTL", "sev", -1),
    ("T CD8 CM", "death", -1),
    ("T CD8 CTL", "death", -1),
    ("HSPC", "death", +1),
    ("CD14mono_IFN", "death", +1),
]
rng_pw = np.random.default_rng(2026)
pw_rows = []
for f, out_type, exp_dir in POWER_SPECS:
    Nn = pd.DataFrame({"y": rint(NEW[f]), "sev": NEW["sev"], "fatal": NEW["fatal"]}).dropna()
    R_valid = R[[f, "sev", "fatal"]].dropna()
    for n_sub in [10, 15, 20, 25, 30, 35]:
        hits = np.zeros(3)
        valid_draws = 0
        attempts = 0
        while valid_draws < 1000 and attempts < 5000:
            attempts += 1
            s = R_valid.loc[rng_pw.choice(R_valid.index, n_sub, replace=False)]
            if out_type == "death" and (s["fatal"].sum() < 1 or (len(s) - s["fatal"].sum()) < 1):
                continue
            valid_draws += 1
            sr = pd.DataFrame({"y": rint(s[f]), "sev": s["sev"], "fatal": s["fatal"]})
            sh = Nn.copy()
            sh["sev"] = rng_pw.permutation(sh["sev"].values)
            sh["fatal"] = rng_pw.permutation(sh["fatal"].values)
            for k, D_sub in enumerate([sr, pd.concat([sr, Nn]), pd.concat([sr, sh])]):
                if out_type == "sev":
                    rr, pp = stats.spearmanr(D_sub["sev"], D_sub["y"])
                    if pp < 0.05 and np.sign(rr) == exp_dir:
                        hits[k] += 1
                else:
                    aa, pp = auc_mw(D_sub["y"], D_sub["fatal"])
                    if np.isfinite(aa) and pp < 0.05 and np.sign(aa - 0.5) == exp_dir:
                        hits[k] += 1
        pw_rows.append({
            "feature": f,
            "outcome": out_type,
            "n_real": n_sub,
            "real_alone": hits[0] / valid_draws,
            "plus_pbisc_new": hits[1] / valid_draws,
            "plus_shuffled": hits[2] / valid_draws,
        })
df_pw = pd.DataFrame(pw_rows)
df_pw.to_csv(DATA_DIR / "power_curves_v8.csv", index=False)

# -------------------------------------------------------------------------
# 8. Modality Baseline Scale Audit (27 Shared Donors) & Patient Normal/Raw Tables
# -------------------------------------------------------------------------
MOD_FEATS = ["NK CD56 bright", "B_ASC_dividing", "T CD8 CM", "CD16mono", "B_ASC_IgG", "B_ASC_IgA", "T CD8 CTL", "HSPC", "CD14mono_IFN"]
mod_rows = []
for f in MOD_FEATS:
    r_sh = R.loc[SHARED, f]
    b_sh = B.loc[SHARED, f]
    w_p = stats.wilcoxon(r_sh, b_sh).pvalue
    rs_sh, ps_sh = stats.spearmanr(r_sh, b_sh)
    raw_c = pd.concat([R[[f, "sev", "fatal"]], NEW[[f, "sev", "fatal"]]])
    ns_c = pd.concat([
        pd.DataFrame({f: rint(R[f]), "sev": R["sev"], "fatal": R["fatal"]}),
        pd.DataFrame({f: rint(NEW[f]), "sev": NEW["sev"], "fatal": NEW["fatal"]}),
    ])
    rho_raw, p_sev_raw = stats.spearmanr(raw_c["sev"], raw_c[f])
    rho_ns, p_sev_ns = stats.spearmanr(ns_c["sev"], ns_c[f])
    auc_raw, p_d_raw = auc_mw(raw_c[f], raw_c["fatal"])
    auc_ns, p_d_ns = auc_mw(ns_c[f], ns_c["fatal"])
    mod_rows.append({
        "feature": f,
        "shared_median_Real_pct": r_sh.median(),
        "shared_IQR_Real_pct": r_sh.quantile(0.75) - r_sh.quantile(0.25),
        "shared_median_PBISC_pct": b_sh.median(),
        "shared_IQR_PBISC_pct": b_sh.quantile(0.75) - b_sh.quantile(0.25),
        "fold_PBISC_over_Real": b_sh.median() / max(r_sh.median(), 1e-4),
        "wilcoxon_paired_p": w_p,
        "shared_spearman_rho": rs_sh,
        "shared_spearman_p": ps_sh,
        "combined_rho_raw_pct": rho_raw,
        "combined_p_sev_raw_pct": p_sev_raw,
        "combined_rho_normal_score": rho_ns,
        "combined_p_sev_normal_score": p_sev_ns,
        "combined_AUC_raw_pct": auc_raw,
        "combined_p_death_raw_pct": p_d_raw,
        "combined_AUC_normal_score": auc_ns,
        "combined_p_death_normal_score": p_d_ns,
    })
df_mod = pd.DataFrame(mod_rows)
df_mod.to_csv(DATA_DIR / "modality_scale_shared27_and_raw_pct_v8.csv", index=False)

# Save shared-donor paired values for plotting in Supp Fig S1
sh_df = pd.DataFrame({"patient_id": SHARED})
for f in MOD_FEATS:
    sh_df[f + "|Real_pct"] = R.loc[SHARED, f].values
    sh_df[f + "|PBISC_pct"] = B.loc[SHARED, f].values
sh_df.to_csv(DATA_DIR / "shared27_paired_pct_v8.csv", index=False)

# Save patient-level Normal scores AND Raw % PBMC for all 66 patients
PT_ns = pd.concat([
    pd.DataFrame({**{f: rint(R[f]) for f in LAB}, "CD16mono|lin": rint(R["CD16mono|lin"]), "sev": R["sev"], "fatal": R["fatal"], "src": "Real"}),
    pd.DataFrame({**{f: rint(NEW[f]) for f in LAB}, "CD16mono|lin": rint(NEW["CD16mono|lin"]), "sev": NEW["sev"], "fatal": NEW["fatal"], "src": "New"}),
])
PT_ns.to_csv(DATA_DIR / "patient_normal_scores_noCLL_v8.csv")

PT_raw = pd.concat([
    pd.DataFrame({**{f: R[f] for f in LAB}, "CD16mono|lin": R["CD16mono|lin"], "sev": R["sev"], "fatal": R["fatal"], "src": "Real"}),
    pd.DataFrame({**{f: NEW[f] for f in LAB}, "CD16mono|lin": NEW["CD16mono|lin"], "sev": NEW["sev"], "fatal": NEW["fatal"], "src": "New"}),
])
PT_raw.to_csv(DATA_DIR / "patient_raw_pct_noCLL_v8.csv")

# -------------------------------------------------------------------------
# 9. NK Compartment Decomposition Table
# -------------------------------------------------------------------------
NKL = ["NK", "NK CD56 bright", "NK IFN stim"]
nkd_rows = []
for k, T_nk in [("Real", R), ("PBISC new", NEW)]:
    for m_name, vec in [
        ("CD56bright % PBMC", T_nk["NK CD56 bright"]),
        ("CD56bright % of NK", 100.0 * T_nk["NK CD56 bright"] / T_nk[NKL].sum(axis=1)),
        ("Total NK % PBMC", T_nk[NKL].sum(axis=1)),
    ]:
        r_nk, p_nk = stats.spearmanr(T_nk["sev"], vec)
        nkd_rows.append(dict(source=k, measure=m_name, rho=r_nk, p=p_nk))
pd.DataFrame(nkd_rows).to_csv(DATA_DIR / "nk_decomposition_v8.csv", index=False)

print("=== FUNNEL SUMMARY ===")
print("S1 Tested:", len(df_sel))
print("S2 Assoc (Real+PBISC):", int((df_sel.S2_sev | df_sel.S2_death).sum()), sorted(df_sel[df_sel.S2_sev | df_sel.S2_death].feature.tolist()))
print("S3 Same dir:", int((df_sel.S3_sev | df_sel.S3_death).sum()), sorted(df_sel[df_sel.S3_sev | df_sel.S3_death].feature.tolist()))
print("S4 Not centre/dex:", int((df_sel.S4_sev | df_sel.S4_death).sum()), sorted(df_sel[df_sel.S4_sev | df_sel.S4_death].feature.tolist()))
print("S5 PBISC-dependent:", int((df_sel.S5_sev | df_sel.S5_death).sum()), sorted(df_sel[df_sel.S5_sev | df_sel.S5_death].feature.tolist()))
print("\n=== CREDIBILITY SCORECARD ===")
print(df_sc.round(3).to_string(index=False))
print("\n=== MODALITY SCALE AUDIT (27 SHARED DONORS) ===")
print(df_mod[["feature", "shared_median_Real_pct", "shared_median_PBISC_pct", "fold_PBISC_over_Real", "wilcoxon_paired_p", "combined_rho_raw_pct", "combined_rho_normal_score", "combined_AUC_raw_pct", "combined_AUC_normal_score"]].round(4).to_string(index=False))
