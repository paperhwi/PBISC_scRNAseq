"""Deterministically export the four primary draw-level and first-draw normal-score matrices
used by 06_audit_metadata_count_pct_tisg.py and 07_build_fig6_fig7_supp6_10_revised.py:
  - data/primary_217draws_rankINT_43subsets.csv
  - data/primary_66firstdraws_rankINT_43subsets.csv
  - data/primary_149draws_HSPC_modules.csv
  - data/primary_61firstdraws_HSPC_modules.csv
"""
import pathlib
import sys
import warnings

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

ROOT = pathlib.Path(r"F:\PBISC")
BASE_COV = ROOT / "HSPC_CD56brightNK_COVID_analysis_summary_20261116"
D28 = BASE_COV / "28_Verified_True11_Deconv_HSC_TISG_20261119" / "data"
FINE_PCT_PATH = BASE_COV / "19_Independent_Fig7_Fig8_Audit_20261118" / "audited_source_snapshot" / "original_composition_inputs" / "sample_level_fine_label_pct.csv"
CELL_DIR = ROOT / "eskd_covid_clinical_assoc_20261115" / "data"
DATA_OUT = ROOT / "deliverables" / "nature_submission_fig4_6_20261120" / "data"
DATA_OUT.mkdir(parents=True, exist_ok=True)
CLL = {"C116", "C139"}


def rint_vdw(x: pd.Series) -> pd.Series:
    s = pd.Series(x, dtype=float)
    n = int(s.notna().sum())
    if n == 0:
        return s
    r = s.rank(method="average")
    return pd.Series(stats.norm.ppf((r - 0.5) / n), index=s.index)


def add_within_arm_rint(df: pd.DataFrame, cols: list, suffix: str) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        z = pd.Series(np.nan, index=out.index, dtype=float)
        for arm in out["arm"].unique():
            m = out["arm"] == arm
            z.loc[m] = rint_vdw(out.loc[m, c])
        out[f"{c}__{suffix}"] = z
    return out


def main():
    meta = pd.read_csv(D28 / "eskd_sample_metadata_verified.csv")
    fine = pd.read_csv(FINE_PCT_PATH)
    scan_prev = pd.read_csv(D28 / "fig7_scan_43fine_true11.csv")
    FINE43 = sorted(scan_prev[scan_prev["CLL_excluded"] == True]["feature"].unique().tolist())
    assert len(FINE43) == 43

    df_all = meta.merge(fine[["source", "sample_id"] + FINE43], on=["source", "sample_id"], how="inner")
    pos_scope = df_all[
        (df_all["condition"] == "positive")
        & df_all["in_primary_scope"].astype(bool)
        & df_all["qc_pass_min200cells"].astype(bool)
        & (~df_all["patient_id"].isin(CLL))
    ].copy()
    real_pts = set(pos_scope[pos_scope["source"] == "Real"]["patient_id"])
    pb_pts = set(pos_scope[pos_scope["source"] == "PBISC"]["patient_id"])
    new29 = pb_pts - real_pts

    prim = pos_scope[
        (pos_scope["source"] == "Real") | ((pos_scope["source"] == "PBISC") & pos_scope["patient_id"].isin(new29))
    ].copy()
    prim["arm"] = np.where(prim["source"] == "Real", "Real", "PBISC new")
    prim["days"] = prim["days_from_first_positive_swab"].astype(float)
    prim = prim.sort_values(["patient_id", "days", "sample_id"]).reset_index(drop=True)
    prim["is_first_draw"] = prim.groupby("patient_id").cumcount() == 0

    p217 = add_within_arm_rint(prim, FINE43, "z217")
    p217.to_csv(DATA_OUT / "primary_217draws_rankINT_43subsets.csv", index=False)

    p66 = p217[p217["is_first_draw"]].copy().reset_index(drop=True)
    p66 = add_within_arm_rint(p66, FINE43, "zfd")
    p66.to_csv(DATA_OUT / "primary_66firstdraws_rankINT_43subsets.csv", index=False)

    r_h5 = ad.read_h5ad(CELL_DIR / "real_rare_subset_cells_counts.h5ad")
    p_h5 = ad.read_h5ad(CELL_DIR / "pbisc_rare_subset_cells_counts.h5ad")
    p_sym = pd.Series(p_h5.var["gene_symbol"].values, index=np.arange(p_h5.n_vars))
    p_sym = p_sym[~p_sym.duplicated()]
    r_sym = pd.Series(r_h5.var_names.values, index=np.arange(r_h5.n_vars))
    r_sym = r_sym[~r_sym.duplicated()]
    shared_g = sorted(set(p_sym.values) & set(r_sym.values))
    r_gmap, p_gmap = {s_: i for i, s_ in r_sym.items()}, {s_: i for i, s_ in p_sym.items()}

    MOD4 = {
        "Ery5": ["GATA1", "KLF1", "TFRC", "APOC1", "BLVRB"],
        "Cycle4": ["MKI67", "TOP2A", "STMN1", "TYMS"],
        "Alarmin_Inflam": ["S100A8", "S100A9", "S100A12", "LYZ", "RETN", "MNDA"],
        "HSC5": ["CD34", "PROM1", "AVP", "HLF", "CRHBP"],
    }
    genes_need = sorted({g for gl in MOD4.values() for g in gl if g in shared_g})

    def extract_hspc(adata, gmap, src, pid_col):
        m = (adata.obs["fine"] == "HSPC").values
        obs = adata.obs.loc[m]
        keep = obs["sample_id"].astype(str).isin(set(prim.loc[prim["source"] == src, "sample_id"])).values
        obs = obs.loc[keep]
        X = adata.X[m][keep][:, [gmap[g] for g in shared_g]]
        X = (X.toarray() if sp.issparse(X) else X).astype(np.float64)
        X = np.log1p(np.divide(X, np.maximum(X.sum(axis=1, keepdims=True), 1.0)) * 1e4)
        out = pd.DataFrame(X[:, [shared_g.index(g) for g in genes_need]], columns=genes_need)
        out["patient_id"] = obs[pid_col].astype(str).values
        out["sample_id"] = obs["sample_id"].astype(str).values
        out["source"] = src
        return out

    hc = pd.concat([extract_hspc(r_h5, r_gmap, "Real", "patient_id"),
                    extract_hspc(p_h5, p_gmap, "PBISC", "donor_id")], ignore_index=True)
    for k, gl in MOD4.items():
        hc[k] = hc[[g for g in gl if g in hc.columns]].mean(axis=1)

    hd_raw = hc.groupby(["source", "sample_id"])[list(MOD4.keys())].mean()
    hd_raw["n_hspc"] = hc.groupby(["source", "sample_id"]).size()
    hd_raw = hd_raw.reset_index()

    hd = prim[["source", "arm", "sample_id", "patient_id", "days", "temp_severity_ord", "severity_ord", "fatal_pub11"]].merge(
        hd_raw, on=["source", "sample_id"], how="inner"
    )
    hd3 = hd[hd["n_hspc"] >= 3].sort_values(["patient_id", "days", "sample_id"]).reset_index(drop=True)
    h149 = add_within_arm_rint(hd3, list(MOD4.keys()), "z149")
    h149.to_csv(DATA_OUT / "primary_149draws_HSPC_modules.csv", index=False)

    h61 = h149.sort_values(["patient_id", "days", "sample_id"]).groupby("patient_id", as_index=False).first()
    h61 = add_within_arm_rint(h61, list(MOD4.keys()), "zfd")
    h61.to_csv(DATA_OUT / "primary_61firstdraws_HSPC_modules.csv", index=False)
    print(f"[OK] Exported primary_217draws ({len(p217)}), primary_66firstdraws ({len(p66)}), primary_149draws ({len(h149)}), primary_61firstdraws ({len(h61)}).")


if __name__ == "__main__":
    main()
