#!/usr/bin/env python3
"""
01_reproduce_adjudicated_tables.py
----------------------------------
Portable, zero-hardcoded-path reproduction script for the 2026-11-26 clinical
outcome adjudication and statistical revisions (Figures 4–7 & Supplementary
Figures 8–27).

Using only repository-bundled plain-text CSV/JSON inputs, this script:
1. Reproduces all 178 baseline external first/last-draw AUCs and P values from
   `analyses/nature_submission_fig4_7_final_20261125/tables/` (if present).
2. Adjudicates all 42 `van der Wijst (2021)` donors against the authors' official
   `Data file S3` (`ClinicalCharacteristics`: `Donor_ID`, `Died` -> 7 deaths),
   resolving 10 first-draw and 12 last-draw discordances (8 donors with visit-
   inconsistent labels).
3. Recomputes the 4-cohort external mortality analysis (219 outcome patients,
   195 survivors, 24 deaths), 10,000-iteration within-outcome bootstrap CIs
   (seed 20261126), and DerSimonian–Laird random-effects meta-analysis (`BH q23`
   across 23 pooled subsets), verifying exact agreement with `tables/`.
4. Recomputes the 16-unique-program HSPC mortality screen (`Supplementary Fig. 23c`:
   33 Real + 29 PBISC new = 62 patients, 11 deaths) from patient module scores.
5. Recomputes the Figure 5f `scSemiProfiler` gene-mean Pearson correlation
   (`r = 0.9365` on full 51,000-cell means vs `r = 0.9363` on display table).
6. Optionally exports the 23-worksheet Source Data Excel workbook when invoked
   with `--export-xlsx`.
"""

from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PKG_DIR = Path(__file__).resolve().parent.parent
TAB_DIR = PKG_DIR / "tables"
VER_DIR = PKG_DIR / "verification_20261126"
BASE_TAB_DIR = PKG_DIR.parent / "nature_submission_fig4_7_final_20261125" / "tables"


def mannwhitney_auc_p_se(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=int)
    pos = x[y == 1]
    neg = x[y == 0]
    if float(np.std(x)) == 0.0:
        return 0.5, 1.0, 0.035
    res = stats.mannwhitneyu(pos, neg, alternative="two-sided", method="auto")
    auc = float(res.statistic / (len(pos) * len(neg)))
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (auc * (1.0 - auc) + (len(pos) - 1) * (q1 - auc * auc) + (len(neg) - 1) * (q2 - auc * auc)) / (
        len(pos) * len(neg)
    )
    se = max(math.sqrt(max(var, 0.0)), 0.035)
    return auc, float(res.pvalue), se


def dl_random_effects(effects: np.ndarray, ses: np.ndarray) -> tuple[float, float, float]:
    effects = np.asarray(effects, dtype=float)
    ses = np.asarray(ses, dtype=float)
    w = 1.0 / (ses ** 2)
    k = len(effects)
    fe = float((w * effects).sum() / w.sum())
    q_stat = float((w * ((effects - fe) ** 2)).sum())
    c_val = float(w.sum() - (w * w).sum() / w.sum())
    tau2 = max(0.0, (q_stat - (k - 1)) / c_val) if (k > 1 and c_val > 0) else 0.0
    w_re = 1.0 / (ses ** 2 + tau2)
    mu = float((w_re * effects).sum() / w_re.sum())
    se = math.sqrt(1.0 / float(w_re.sum()))
    p = float(2.0 * stats.norm.sf(abs(mu / se)))
    return mu, se, p


def export_xlsx_from_csvs(xlsx_out: Path) -> None:
    manifest = pd.read_csv(TAB_DIR / "00_table_manifest.csv")
    readme_rows = json.loads((VER_DIR / "workbook_readme_original.json").read_text(encoding="utf-8"))
    wb = openpyxl.Workbook()
    ws0 = wb.active
    ws0.title = "00_README_Index"
    for r in readme_rows:
        ws0.append([None if v == "None" else v for v in r])
    for row in manifest.itertuples(index=False):
        ws = wb.create_sheet(title=row.excel_sheet)
        df = pd.read_csv(TAB_DIR / row.csv_filename)
        ws.append(list(df.columns))
        for rec in df.astype(object).where(pd.notna(df), None).values.tolist():
            ws.append(rec)
    wb.save(str(xlsx_out))
    print(f"[XLSX] Exported 23-sheet workbook -> {xlsx_out.name} ({xlsx_out.stat().st_size:,} B)")


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce 2026-11-26 adjudicated tables and statistical revisions.")
    parser.add_argument("--export-xlsx", action="store_true", help="Optionally rebuild the 23-sheet Source Data .xlsx workbook.")
    parser.add_argument("--write-tables", action="store_true", help="Overwrite tables/*.csv with recomputed outputs (default: verify in-memory).")
    args = parser.parse_args()

    # 1. Load official van der Wijst clinical rows
    vdw_rows = json.loads((VER_DIR / "clinical_sources" / "vdw_clinical_rows.json").read_text(encoding="utf-8"))
    raw_vdw = pd.DataFrame(vdw_rows[1:], columns=vdw_rows[0])
    raw_vdw = raw_vdw.loc[:, raw_vdw.columns.notna()]
    clin_vdw = raw_vdw[["Donor_ID", "Died"]].dropna(subset=["Donor_ID"])
    assert not clin_vdw["Donor_ID"].duplicated().any()
    vdw_map = {f"Rep_C_{int(r.Donor_ID)}": int(str(r.Died).strip().lower() == "yes") for r in clin_vdw.itertuples()}
    assert len(vdw_map) >= 42

    # 2. If baseline 2026-11-25 tables are available in sibling directory, verify 178 historical effects & discordances
    if BASE_TAB_DIR.exists():
        old_ext = pd.read_csv(BASE_TAB_DIR / "fig7_external4_300_acute_samples.csv")
        old_per = pd.read_csv(BASE_TAB_DIR / "fig7_external4_per_cohort_results.csv")
        old_meta = pd.read_csv(BASE_TAB_DIR / "fig7_external4_meta_23subsets.csv")
        old_panel = pd.read_csv(BASE_TAB_DIR / "fig7_panel_l_summary.csv")

        v_old = old_ext.loc[old_ext["cohort"] == "van der Wijst 2021"].sort_values(["donor", "time_ord", "sample_id"])
        rec = v_old.drop_duplicates("donor", keep="first")[["donor", "death"]].rename(columns={"death": "old_first_death"})
        rec = rec.merge(
            v_old.drop_duplicates("donor", keep="last")[["donor", "death"]].rename(columns={"death": "old_last_death"}),
            on="donor",
        )
        rec["official_death"] = rec["donor"].map(vdw_map)
        rec["n_old_distinct_death_labels"] = rec["donor"].map(v_old.groupby("donor")["death"].nunique())
        assert len(rec) == 42 and int(rec["official_death"].sum()) == 7
        assert int((rec["old_first_death"] != rec["official_death"]).sum()) == 10
        assert int((rec["old_last_death"] != rec["official_death"]).sum()) == 12
        assert int((rec["n_old_distinct_death_labels"] > 1).sum()) == 8

        n_repro_178 = 0
        for _, r in old_per.iterrows():
            for kind, keep in (("first", "first"), ("last", "last")):
                ss = (
                    old_ext[old_ext["cohort"] == r["cohort"]]
                    .dropna(subset=[f"pct_{r['feature']}", "who_sameday"])
                    .sort_values(["donor", "time_ord", "sample_id"])
                    .drop_duplicates("donor", keep=keep)
                    .dropna(subset=["death"])
                )
                a, p, _ = mannwhitney_auc_p_se(ss[f"pct_{r['feature']}"].to_numpy(), ss["death"].to_numpy())
                ac = "auc_death_first" if kind == "first" else "auc_death"
                pc = "p_death_first" if kind == "first" else "p_death"
                assert abs(a - r[ac]) < 1e-12 and abs(p - r[pc]) < 1e-12
                n_repro_178 += 1
        print(f"[1/4] Baseline audit verified: {n_repro_178} historical effects reproduced (max diff < 1e-12); "
              f"10 first-draw & 12 last-draw vdW discordances (8 visit-inconsistent donors) confirmed.")
        e_work = old_ext.copy()
        per_work = old_per.copy()
        meta_work = old_meta.copy()
        panel_work = old_panel.rename(columns={"vanDerWijst_42_5_AUC_1st": "vanDerWijst_42_7_AUC_1st"}).copy()
    else:
        e_work = pd.read_csv(TAB_DIR / "fig7_external4_300_acute_samples.csv")
        per_work = pd.read_csv(TAB_DIR / "fig7_external4_per_cohort_results.csv")
        meta_work = pd.read_csv(TAB_DIR / "fig7_external4_meta_23subsets.csv")
        panel_work = pd.read_csv(TAB_DIR / "fig7_panel_l_summary.csv")
        print("[1/4] Baseline sibling folder not present; reproducing directly from bundled adjudicated inputs.")

    # 3. Apply official van der Wijst outcomes and recompute per-cohort + DL meta-analysis
    vd_mask = e_work["cohort"] == "van der Wijst 2021"
    assert e_work.loc[vd_mask, "donor"].isin(vdw_map).all()
    e_work.loc[vd_mask, "death"] = e_work.loc[vd_mask, "donor"].map(vdw_map).astype(float)

    rng = np.random.default_rng(20261126)
    for i, r in per_work[per_work["cohort"] == "van der Wijst 2021"].iterrows():
        sub = e_work[e_work["cohort"] == r["cohort"]].dropna(subset=[f"pct_{r['feature']}", "who_sameday"])
        for kind, keep in (("first", "first"), ("last", "last")):
            ss = sub.sort_values(["donor", "time_ord", "sample_id"]).drop_duplicates("donor", keep=keep).dropna(subset=["death"])
            x = ss[f"pct_{r['feature']}"].to_numpy(dtype=float)
            y = ss["death"].to_numpy(dtype=int)
            a, p, se = mannwhitney_auc_p_se(x, y)
            suf = "_first" if kind == "first" else ""
            per_work.loc[i, [f"n_surv{suf}", f"n_dead{suf}", f"auc_death{suf}", f"se_auc{suf}", f"p_death{suf}"]] = [
                int(np.sum(y == 0)),
                int(np.sum(y == 1)),
                a,
                se,
                p,
            ]
            xx, yy = x[y == 1], x[y == 0]
            kern = (xx[:, None] > yy).astype(float) + 0.5 * (xx[:, None] == yy)
            w1 = rng.multinomial(len(xx), np.ones(len(xx)) / len(xx), size=10_000)
            w0 = rng.multinomial(len(yy), np.ones(len(yy)) / len(yy), size=10_000)
            boot = ((w1 @ kern) * w0).sum(axis=1) / float(len(xx) * len(yy))
            lo, hi = np.quantile(boot, [0.025, 0.975])
            per_work.loc[i, [f"ci_lo_auc{suf}", f"ci_hi_auc{suf}"]] = [float(lo), float(hi)]

    for i, r in meta_work.iterrows():
        pp = per_work[per_work["feature"] == r["feature"]]
        if len(pp) < 2:
            continue
        meta_work.loc[i, ["n_surv_ext", "n_dead_ext"]] = [int(pp["n_surv"].sum()), int(pp["n_dead"].sum())]
        for kind in ("first", "last"):
            suf = "_first" if kind == "first" else ""
            mu, se, p = dl_random_effects(pp[f"auc_death{suf}"].to_numpy() - 0.5, pp[f"se_auc{suf}"].to_numpy())
            pref = f"ext_{kind}draw_"
            meta_work.loc[i, [f"{pref}AUC_re", f"{pref}ci_lo", f"{pref}ci_hi", f"{pref}p"]] = [
                mu + 0.5,
                mu + 0.5 - 1.96 * se,
                mu + 0.5 + 1.96 * se,
                p,
            ]

    exact_mask = meta_work["feature"] != "B_ASC_total"
    for end in ("first", "last"):
        meta_work.loc[exact_mask, f"ext_{end}draw_q23"] = multipletests(
            meta_work.loc[exact_mask, f"ext_{end}draw_p"].to_numpy(dtype=float), method="fdr_bh"
        )[1]

    col_map = {
        "RE_AUC_1st": "ext_firstdraw_AUC_re",
        "RE_CI_lo_1st": "ext_firstdraw_ci_lo",
        "RE_CI_hi_1st": "ext_firstdraw_ci_hi",
        "RE_P_1st": "ext_firstdraw_p",
        "RE_q23_1st": "ext_firstdraw_q23",
        "RE_AUC_last": "ext_lastdraw_AUC_re",
        "RE_CI_lo_last": "ext_lastdraw_ci_lo",
        "RE_CI_hi_last": "ext_lastdraw_ci_hi",
        "RE_P_last": "ext_lastdraw_p",
        "RE_q23_last": "ext_lastdraw_q23",
    }
    for i, r in panel_work.iterrows():
        mm = meta_work[meta_work["feature"] == r["feature"]].iloc[0]
        vv = per_work[(per_work["cohort"] == "van der Wijst 2021") & (per_work["feature"] == r["feature"])].iloc[0]
        for k_dst, k_src in col_map.items():
            panel_work.loc[i, k_dst] = mm[k_src]
        panel_work.loc[i, "vanDerWijst_42_7_AUC_1st"] = vv["auc_death_first"]

    # Verify exact numerical match against bundled adjudicated CSVs
    for fn, df_calc in (
        ("fig7_external4_300_acute_samples.csv", e_work),
        ("fig7_external4_per_cohort_results.csv", per_work),
        ("fig7_external4_meta_23subsets.csv", meta_work),
        ("fig7_panel_l_summary.csv", panel_work),
    ):
        df_ref = pd.read_csv(TAB_DIR / fn)
        pd.testing.assert_frame_equal(df_calc, df_ref, check_dtype=False, rtol=0.0, atol=1e-11)
        if args.write_tables:
            df_calc.to_csv(TAB_DIR / fn, index=False)

    print("[2/4] External clinical adjudication & 4-cohort meta-analysis (219 patients, 24 deaths) reproduced (max diff < 1e-11).")

    # 4. Reproduce HSPC unique-16 module screen from bundled 62-patient module scores
    hspc_pts = pd.read_csv(VER_DIR / "supp23_hspc62_patient_module_scores.csv")
    hspc_ref = pd.read_csv(VER_DIR / "supp23_hspc16_recomputed.csv")
    hspc_defs = json.loads((VER_DIR / "supp23_hspc16_definitions.json").read_text(encoding="utf-8"))

    def vdw_norm(s: pd.Series) -> pd.Series:
        return pd.Series(stats.norm.ppf((s.rank(method="average") - 0.5) / len(s)), index=s.index)

    hspc_rows = []
    for mod in hspc_defs:
        sub_m = hspc_pts[hspc_pts["module"] == mod]
        s_real = sub_m[sub_m["arm"] == "Real"].set_index("patient_id")["score"]
        s_new = sub_m[sub_m["arm"] == "PBISC new"].set_index("patient_id")["score"]
        y_map = sub_m.set_index("patient_id")["fatal"].to_dict()
        s_comb = pd.concat([vdw_norm(s_real), vdw_norm(s_new)])
        for arm, s in (("Real", s_real), ("PBISC new", s_new), ("Combined", s_comb)):
            yb = np.array([y_map[pid] for pid in s.index], dtype=int)
            auc_v, p_v, _ = mannwhitney_auc_p_se(s.to_numpy(dtype=float), yb)
            hspc_rows.append({"module": mod, "cohort": arm, "n": len(s), "deaths": int(yb.sum()), "AUC": auc_v, "p": p_v})
    df_hspc = pd.DataFrame(hspc_rows)
    df_hspc["q16"] = np.nan
    for _, ix in df_hspc.groupby("cohort").groups.items():
        df_hspc.loc[ix, "q16"] = multipletests(df_hspc.loc[ix, "p"].to_numpy(dtype=float), method="fdr_bh")[1]
    pd.testing.assert_frame_equal(df_hspc, hspc_ref, check_dtype=False, rtol=0.0, atol=1e-7)
    print("[3/4] Supplementary Fig. 23c 16-unique-program HSPC screen (48 estimates across 62 patients) reproduced (max diff < 1e-7).")

    # 5. Reproduce Figure 5f scSemiProfiler correlation check (log1p gene means)
    disp = pd.read_csv(VER_DIR / "sepsis_gene_display_source.csv")
    full_m = pd.read_csv(VER_DIR / "sepsis_scsemi_full_means.csv")
    m_sc = disp.merge(full_m, on="gene")
    r_full = float(np.corrcoef(np.log1p(m_sc["Real"]), np.log1p(m_sc["full_scsemi_mean"]))[0, 1])
    r_disp = float(np.corrcoef(np.log1p(disp["Real"]), np.log1p(disp["scSemiProfiler"]))[0, 1])
    assert abs(r_full - 0.9365457) < 1e-5 and abs(r_disp - 0.9363051) < 1e-5
    print(f"[4/4] Figure 5f scSemiProfiler log1p gene-mean correlation reproduced: full={r_full:.4f} (0.937), display={r_disp:.4f} (0.936).")

    if args.export_xlsx:
        export_xlsx_from_csvs(PKG_DIR / "PBISC_SourceData_and_Summary_Tables_Fig4_7_20261125.xlsx")


if __name__ == "__main__":
    main()
