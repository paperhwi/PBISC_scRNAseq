"""Rare / disease-associated subsets with parent-lineage denominators (review item 2).

For each subset: parent lineage defined explicitly below; per group (case, control) and source (real reference, K1000, K4000):
  pooled within-parent fraction = target cells / parent cells (numerator and denominator reported),
  median per-sample within-parent fraction over samples with >= 20 parent cells,
  samples with >= 1 target cell / samples with >= 20 parent cells.
Generated labels are kNN-transferred (c_joint.py); generated sample = bulk sample, pooled over the five draws
(100 cells per sample per draw in the joint subset). Real = full recapitulated reference (sepsis: in_paper_scheme;
COVID: published 187-sample roster)."""
import json, numpy as np, pandas as pd
R = "/home/claude/work/report"; MINP = 20

SEPSIS = {  # subset: (parent description, test on major label)
    "Mono_c06-CD14-MKI67": ("Monocytes", lambda m, l: m == "Monocytes"),
    "B_c08-MKI67": ("B + plasma cells", lambda m, l: np.isin(m, ["B cells", "Plasma cells"])),
    "T_CD4_c07-MKI67": ("CD4+ T cells", lambda m, l: m == "CD4+ T cells"),
    "T_CD8_c06-MKI67": ("CD8+ T cells", lambda m, l: m == "CD8+ T cells"),
    "NK_c04-MKI67": ("NK cells", lambda m, l: m == "NK cells"),
    "T_CD4_c02-GPR183-NR4A2": ("CD4+ T cells", lambda m, l: m == "CD4+ T cells"),
    "T_CD8_c05-GNLY-CCL4hi": ("CD8+ T cells", lambda m, l: m == "CD8+ T cells"),
    "NKT_c05-IFNG-CCL4hi": ("NKT cells", lambda m, l: m == "NKT cells"),
    "NK_c03-PRF1-CCL4hi": ("NK cells", lambda m, l: m == "NK cells"),
    "Mono_c01-CD14-HLA-DRA": ("Monocytes", lambda m, l: m == "Monocytes"),
    "Mono_c03-CD14-RETN": ("Monocytes", lambda m, l: m == "Monocytes"),
    "B_c07-MZB1": ("B + plasma cells", lambda m, l: np.isin(m, ["B cells", "Plasma cells"])),
    "pDC": ("DC", lambda m, l: m == "DC"),
    "cDC": ("DC", lambda m, l: m == "DC"),
    "Macrophage": ("Monocytes + macrophages", lambda m, l: np.isin(m, ["Monocytes", "Macrophages"])),
    "gdT": ("all T (CD4, CD8, NKT, gdT)", lambda m, l: np.isin(m, ["CD4+ T cells", "CD8+ T cells", "NKT cells"]) | (l == "gdT")),
    "ILC": ("NK + ILC", lambda m, l: (m == "NK cells") | (l == "ILC")),
    "HSC": ("all PBMC", lambda m, l: np.ones(len(m), bool)),
}
SEP_CLAIM = {"T_CD4_c02-GPR183-NR4A2": "SEP-1", "T_CD8_c05-GNLY-CCL4hi": "SEP-2", "NKT_c05-IFNG-CCL4hi": "SEP-3", "NK_c03-PRF1-CCL4hi": "SEP-3-NK",
             "Mono_c01-CD14-HLA-DRA": "SEP-4A", "Mono_c03-CD14-RETN": "SEP-4B", "Mono_c06-CD14-MKI67": "SEP-5"}
st = lambda m, *p: np.char.startswith(np.asarray(m, dtype=str), p) if len(p) == 1 else np.any([np.char.startswith(np.asarray(m, dtype=str), q) for q in p], axis=0)
COVID = {
    "CD14mono_dex": ("CD14 monocytes", lambda m, l: st(m, "CD14mono")),
    "CD14mono_IFN": ("CD14 monocytes", lambda m, l: st(m, "CD14mono")),
    "B_ASC_IgG": ("plasmablasts (B_ASC)", lambda m, l: st(m, "B_ASC")),
    "B_ASC_dividing": ("plasmablasts (B_ASC)", lambda m, l: st(m, "B_ASC")),
    "T CD4 cycling": ("CD4 T (incl. Treg)", lambda m, l: st(m, "T CD4", "T reg")),
    "ASDC": ("DC (DC1-3, ASDC, pDC)", lambda m, l: st(m, "DC", "ASDC", "pDC")),
    "DC1": ("DC (DC1-3, ASDC, pDC)", lambda m, l: st(m, "DC", "ASDC", "pDC")),
    "pDC": ("DC (DC1-3, ASDC, pDC)", lambda m, l: st(m, "DC", "ASDC", "pDC")),
    "HSPC": ("all PBMC", lambda m, l: np.ones(len(m), bool)),
    "ILC": ("NK + ILC", lambda m, l: st(m, "NK", "ILC")),
    "T g/d": ("all T (incl. MAIT, NKT)", lambda m, l: st(m, "T ", "MAIT", "NKT")),
    "MAIT": ("all T (incl. MAIT, NKT)", lambda m, l: st(m, "T ", "MAIT", "NKT")),
    "NK CD56 bright": ("NK", lambda m, l: st(m, "NK")),
}


def load(coh):
    if coh == "sepsis":
        ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                                    usecols=["sample", "age_group", "group", "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme"]) for i in range(3)], ignore_index=True)
        ro = ro[ro.in_paper_scheme.astype(str).eq("True")]
        real = pd.DataFrame({"sample": ro["sample"].astype(str), "label": ro.celltype_ye2026, "major": ro.celltype_major_ye2026,
                             "grp": ro.age_group.str[:5].replace({"paedi": "ped"}) + "_" + ro.group.replace({"healthy control": "HC", "ICU control": "ICU_ctrl"})})
        return real, SEPSIS, "adult_sepsis", "adult_HC"
    real = pd.read_csv(f"{R}/covid_real_labels.csv.gz", keep_default_na=False)
    return real, COVID, "covid_positive", "covid_negative"


def summarise(df, subset, parent_fn, grp, sample_col):
    d = df[df.grp == grp]
    par = parent_fn(d.major.to_numpy(), d.label.to_numpy()); tgt = par & (d.label.to_numpy() == subset)
    per = pd.DataFrame({"s": d[sample_col].to_numpy(), "p": par, "t": tgt}).groupby("s")[["p", "t"]].sum()
    per = per[per.p >= MINP]
    return {"n_target": int(tgt.sum()), "n_parent": int(par.sum()), "pooled_pct": 100 * tgt.sum() / max(par.sum(), 1),
            "median_sample_pct": float(100 * (per.t / per.p).median()) if len(per) else np.nan,
            "samples_with_target": int((per.t > 0).sum()), "samples_evaluable": int(len(per))}


rows = []
for coh in ["sepsis", "covid"]:
    real, defs, case, ctrl = load(coh)
    gens = {K: pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False).rename(columns={"group": "grp"}) for K in ["k1000", "k4000"]}
    for subset, (pdesc, fn) in defs.items():
        for grp in [case, ctrl]:
            row = {"cohort": coh, "subset": subset, "claim": SEP_CLAIM.get(subset, ""), "parent": pdesc, "group": grp}
            for src, df, sc in [("real", real, "sample")] + [(K, g, "sample_id") for K, g in gens.items()]:
                row.update({f"{src}_{k}": v for k, v in summarise(df, subset, fn, grp, sc).items()})
            rows.append(row)
out = pd.DataFrame(rows); out.to_csv(f"{R}/rare_parent_denominators.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
show = ["cohort", "subset", "group", "real_pooled_pct", "real_median_sample_pct", "real_samples_with_target", "real_samples_evaluable",
        "k1000_pooled_pct", "k1000_n_target", "k1000_n_parent", "k1000_samples_with_target", "k1000_samples_evaluable", "k4000_pooled_pct", "k4000_samples_with_target"]
print(out[show].round(3).to_string())
