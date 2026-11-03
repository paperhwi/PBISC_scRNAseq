"""Sample-level validation of paper-defined populations in the original cohorts.

This freezes claim labels and marker panels before any K=1,000/4,000 PBISC
comparison; Ye marker panels are populated from Supplementary Table 6 when imported. Frequencies are aggregated per biological sample; cell counts
are never treated as independent replicates.
"""
from pathlib import Path
import gzip
import h5py
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "bulkonly_lora5_20261006" / "analysis" / "original_scRNA"
OUT = ROOT / "paper_claim_analysis"
OUT.mkdir(exist_ok=True)

# Frozen Ye article claim labels and fallback core markers. If the full Supplementary
# Table 6 import exists, write_fixed_claims() expands each Ye panel from that table.
SEPSIS_CLAIMS = [
    dict(id="SEP-1", population="T_CD4_c02-GPR183-NR4A2", parent="CD4 T",
         markers="GPR183;NR4A2", age="adult", contrast="infection_site vs healthy control",
         target="T_CD4_c02-GPR183-NR4A2", lineage_prefix="T_CD4_"),
    dict(id="SEP-2", population="T_CD8_c05-GNLY-CCL4hi", parent="CD8 T",
         markers="GNLY;CCL3;CCL4;TNF", age="adult", contrast="abdominal/respiratory vs healthy control",
         target="T_CD8_c05-GNLY-CCL4hi", lineage_prefix="T_CD8_"),
    dict(id="SEP-3", population="NKT_c05-IFNG-CCL4hi", parent="NKT",
         markers="IFNG;CCL3;CCL4;GNLY", age="adult", contrast="infection_site vs healthy control",
         target="NKT_c05-IFNG-CCL4hi", lineage_prefix="NKT_"),
    dict(id="SEP-3-NK", population="NK_c03-PRF1-CCL4hi", parent="NK",
         markers="PRF1;GNLY;CCL3;CCL4;TNF", age="adult", contrast="abdominal/respiratory vs healthy control",
         target="NK_c03-PRF1-CCL4hi", lineage_prefix="NK_"),
    dict(id="SEP-4A", population="Mono_c01-CD14-HLA-DRA", parent="Monocyte",
         markers="HLA-DRA;HLA-DRB1", age="adult", contrast="infection_site vs healthy control",
         target="Mono_c01-CD14-HLA-DRA", lineage_prefix="Mono_"),
    dict(id="SEP-4B", population="Mono_c03-CD14-RETN", parent="Monocyte",
         markers="RETN;S100A8;S100A9;IL1R2", age="adult", contrast="infection_site vs healthy control",
         target="Mono_c03-CD14-RETN", lineage_prefix="Mono_"),
    dict(id="SEP-5", population="Mono_c06-CD14-MKI67", parent="CD14 monocyte",
         markers="MKI67;TOP2A;TYMS", age="paediatric", contrast="respiratory sepsis vs paediatric healthy control",
         target="Mono_c06-CD14-MKI67", lineage_prefix="Mono_"),
]


def read_sepsis_obs():
    files = sorted((BASE / "sepsis_GSE279452_full" / "recap_annotation").glob("sepsis_obs_part*.csv.gz"))
    use = ["cell_id", "sample", "group", "age_group", "infection_site", "major_lineage",
           "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme"]
    return pd.concat([pd.read_csv(f, usecols=use, low_memory=False) for f in files], ignore_index=True)


def write_fixed_claims():
    rows = []
    for d in SEPSIS_CLAIMS:
        rows.append({"study":"Ye et al. 2026", **d,
                     "definition_source":"article-defined subset name + primary reported marker genes; validate exact gene roster against Supplementary Table 6"})
    rows.extend([
        {"study":"Stephenson et al. 2025", "id":"ESKD-IFN", "population":"IFN-stimulated states within each parent cell type",
         "parent":"each lineage", "markers":"ISG15;IFIT1;IFIT3;MX1;OAS1;OAS2;IFI6", "age":"all",
         "contrast":"sample-level COVID status / time from symptom onset; original trajectory claim requires repeated samples",
         "target":"cell_type_IFNstim label", "lineage_prefix":"within original cell_type_detailed",
         "definition_source":"fixed core ISG module; full paper Supplementary Tables 7-8 / authors’ complete pathway signatures to be added before K comparison"},
        {"study":"Stephenson et al. 2025", "id":"ESKD-DEX", "population":"dexamethasone-associated monocytes (dex. mono)",
         "parent":"monocytes", "markers":"CD163;ADAMTS2;FKBP5;TSC22D3;SLC1A3;MFGE8;FPR1;MTSS1;MRC1;SESN1",
         "age":"all", "contrast":"within-patient pre/post glucocorticoid; assess cohort/wave and severity",
         "target":"monocyte compartment cluster; dex_mono_recap only as secondary label",
         "lineage_prefix":"CD14/CD16/intermediate monocytes",
         "definition_source":"Figure 5 reported marker direction; local recap cluster-23 rule is not independent and must not be used as primary validation"},
    ])
    definitions = pd.DataFrame(rows)
    definitions["core_markers"] = definitions["markers"]
    definitions["author_cluster_name"] = ""
    definitions["author_marker_count"] = ""
    marker_file = OUT / "ye_supplementary_table6_all_markers.csv"
    if marker_file.exists():
        marker_rows = pd.read_csv(marker_file)
        claim_cluster = {
            "SEP-1":"T_CD4_c02-GPR183-NR4A2", "SEP-2":"T_CD8_c05-GNLY-CCL4",
            "SEP-3":"NKT_c05-IFNG-CCL4", "SEP-3-NK":"NK_c03-PRF1-CCL4",
            "SEP-4A":"Mono-c01-CD14-HLA-DRA", "SEP-4B":"Mono-c03-CD14-RETN",
            "SEP-5":"Mono-c06-CD14-MKI67",
        }
        for i, row in definitions[definitions["study"].eq("Ye et al. 2026")].iterrows():
            cluster = claim_cluster[row["id"]]
            genes = marker_rows.loc[marker_rows["cluster"].eq(cluster), "gene"].astype(str)
            definitions.at[i, "author_cluster_name"] = cluster
            definitions.at[i, "author_marker_count"] = str(genes.nunique())
            definitions.at[i, "markers"] = ";".join(genes)
            definitions.at[i, "definition_source"] = "full author marker roster from Supplementary Table 6; exact publisher cluster label mapped; previous core panel retained in core_markers"
    definitions.to_csv(OUT / "fixed_paper_claim_definitions.csv", index=False)


def sample_frequency(data, claim, keep_mask=None):
    d = data if keep_mask is None else data.loc[keep_mask]
    parent = d["celltype_ye2026"].astype(str).str.startswith(claim["lineage_prefix"])
    # For the reported rare cell-state tests, denominator is the annotated parent lineage.
    p = d.loc[parent].copy()
    p["target"] = (p["celltype_ye2026"].astype(str) == claim["target"])
    out = p.groupby(["sample", "group", "age_group", "infection_site"], dropna=False).agg(
        parent_cells=("target", "size"), target_cells=("target", "sum"),
    ).reset_index()
    out["fraction_parent"] = out["target_cells"] / out["parent_cells"].replace(0, np.nan)
    out["claim_id"] = claim["id"]
    out["population"] = claim["population"]
    return out


def sepsis_analysis():
    obs = read_sepsis_obs()
    # Match original authors' reported 44-subset scheme where possible.
    obs = obs[obs["in_paper_scheme"].astype(bool)].copy()
    all_rows = []
    for c in SEPSIS_CLAIMS:
        mask = obs["age_group"].astype(str).eq(c["age"])
        all_rows.append(sample_frequency(obs, c, mask))
    sample = pd.concat(all_rows, ignore_index=True)
    sample.to_csv(OUT / "original_sepsis_claim_frequencies_by_sample.csv", index=False)
    summaries = []
    for claim in SEPSIS_CLAIMS:
        d = sample[sample.claim_id == claim["id"]].copy()
        for (age, group, site), g in d.groupby(["age_group", "group", "infection_site"], dropna=False):
            summaries.append({"claim_id":claim["id"], "population":claim["population"], "age_group":age,
                             "group":group, "infection_site":site, "n_samples":g["sample"].nunique(),
                             "median_fraction_parent":g["fraction_parent"].median(),
                             "mean_fraction_parent":g["fraction_parent"].mean(),
                             "samples_with_cells":int((g["target_cells"] > 0).sum()),
                             "median_target_cells":g["target_cells"].median()})
        # Disease/site contrast vs age-matched healthy controls; Mann-Whitney at sample level.
        healthy = d[d.group.eq("healthy control")]["fraction_parent"].dropna()
        for site in sorted(d.loc[d.group.eq("sepsis"), "infection_site"].dropna().unique()):
            x = d.loc[d.group.eq("sepsis") & d.infection_site.eq(site), "fraction_parent"].dropna()
            if len(x) and len(healthy):
                summaries.append({"claim_id":claim["id"], "population":claim["population"],
                                  "age_group":claim["age"], "group":"sepsis_vs_HC_test", "infection_site":site,
                                  "n_samples":len(x), "median_fraction_parent":x.median(),
                                  "mean_fraction_parent":x.mean(), "samples_with_cells":int((x>0).sum()),
                                  "median_target_cells":np.nan, "HC_n":len(healthy),
                                  "HC_median_fraction_parent":healthy.median(),
                                  "MWU_pvalue":mannwhitneyu(x, healthy, alternative="two-sided").pvalue})
    summary = pd.DataFrame(summaries)
    tests = summary["MWU_pvalue"].notna()
    summary["BH_qvalue_across_all_sepsis_claim_site_tests"] = np.nan
    if tests.any():
        summary.loc[tests, "BH_qvalue_across_all_sepsis_claim_site_tests"] = multipletests(
            summary.loc[tests, "MWU_pvalue"].astype(float), method="fdr_bh")[1]
    summary.to_csv(OUT / "original_sepsis_claim_summary.csv", index=False)
    return obs, sample


def decode_cat(group, n):
    codes = group["codes"][:]
    cats = group["categories"][:]
    vals = np.full(n, "", dtype=object)
    ok = codes >= 0
    decoded = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in cats], dtype=object)
    vals[ok] = decoded[codes[ok]]
    return vals


def covid_analysis():
    h5 = BASE / "covid_eskd" / "eskd_covid19_recap.h5ad"
    ann = pd.read_csv(BASE / "covid_eskd" / "recap_annotation" / "covid_obs_recap.csv.gz")
    with h5py.File(h5, "r") as f:
        n = len(f["obs/_index"])
        md = {}
        for k in ["patient_id", "sample_id", "case_control", "WHO_severity", "cell_type_detailed",
                  "cell_type_IFNstim", "cell_type_compartment", "orig.ident", "centre", "days_since_symptoms",
                  "exclude", "sample_date", "days_since_posswab", "fatal_disease"]:
            o = f["obs"][k]
            if isinstance(o, h5py.Group):
                md[k] = decode_cat(o, n)
            else:
                md[k] = o[:]
        md["cell_id"] = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/_index"][:]], dtype=object)
    orig = pd.DataFrame(md)
    rec = ann.rename(columns={"_index":"cell_id"})
    d = orig.merge(rec[["cell_id", "cell_type_detailed_recap", "cell_type_IFNstim_recap", "dex_mono_recap", "recap_subcluster"]],
                   on="cell_id", how="left", validate="one_to_one")
    d = d[pd.to_numeric(d["exclude"], errors="coerce").fillna(1).eq(0)].copy()
    local_cells, local_samples, local_patients = len(d), d.sample_id.nunique(), d.patient_id.nunique()
    # Restrict COVID summaries to the exact publication sample roster (Table S1).
    roster = pd.read_csv(OUT / "covid_published_187_sample_ids.csv", dtype={"sample_id":str})
    def sample_key(s):
        prefix, sep, suffix = str(s).strip().partition("_")
        return f"{prefix}_{suffix[:4]}" if sep and len(suffix) >= 4 else str(s).strip()
    keep_keys = set(roster["sample_id"].map(sample_key))
    d = d[d["sample_id"].astype(str).map(sample_key).isin(keep_keys)].copy()
    d.to_csv(OUT / "original_covid_claim_cells_metadata.csv.gz", index=False, compression="gzip")
    # Raw labels retain the deposited source annotations. These results do not use generated cells.
    # IFN-stimulated frequency among each original detailed population, grouped per sample.
    # The deposited cell_type_IFNstim field stores IFN-high/low states for B cells,
    # explicit IFN labels for myeloid and T/NK cells, and the unsplit label otherwise.
    d["ifn"] = d["cell_type_IFNstim"].astype(str).str.contains(r"IFN|_High$", case=False, na=False, regex=True)
    d["dex_label"] = d["cell_type_detailed_recap"].eq("CD14mono_dex") | d["dex_mono_recap"].fillna(False).astype(bool)
    # IFN states are tested within their deposited detailed parent type, not across
    # the entire compartment. Use sample_id as the biological unit; orig.ident/pool
    # is only a library batch and may contain multiple libraries per clinical sample.
    out = []
    keys = ["sample_id", "patient_id", "case_control", "WHO_severity", "centre", "cell_type_detailed"]
    for key, gc in d.groupby(keys, dropna=False, observed=True):
        sid, pid, status, sev, centre, dtype = key
        n_ifn = int(gc["ifn"].sum())
        out.append({"sample_id":sid,"patient_id":pid,"case_control":status,"WHO_severity":sev,
                    "centre":centre,"population":dtype,"n_parent_cells":len(gc),"ifn_cells":n_ifn,
                    "ifn_fraction_parent":n_ifn/len(gc) if len(gc) else np.nan})
    ifn_sf = pd.DataFrame(out)
    ifn_sf.to_csv(OUT / "original_covid_ifn_frequencies_by_sample_and_parent.csv", index=False)
    # Dexamethasone-associated monocytes are a monocyte subpopulation. The local
    # cluster-23 annotation is secondary evidence; treatment metadata are absent.
    d["is_monocyte"] = d["cell_type_detailed"].astype(str).str.contains("mono", case=False, na=False)
    mono = d[d.is_monocyte].copy()
    dex_sf = mono.groupby(["sample_id", "patient_id", "case_control", "WHO_severity", "centre"], dropna=False).agg(
        mono_cells=("dex_label", "size"), dex_cells=("dex_label", "sum")
    ).reset_index()
    dex_sf["dex_fraction_monocytes"] = dex_sf["dex_cells"] / dex_sf["mono_cells"].replace(0, np.nan)
    dex_sf.to_csv(OUT / "original_covid_dex_mono_frequencies_by_sample.csv", index=False)
    dex_summary = dex_sf.groupby(["case_control", "WHO_severity", "centre"], dropna=False).agg(
        n_samples=("sample_id", "nunique"), n_patients=("patient_id", "nunique"),
        samples_with_candidate_dex=("dex_cells", lambda x: int((x > 0).sum())),
        median_candidate_dex_fraction=("dex_fraction_monocytes", "median"),
        mean_candidate_dex_fraction=("dex_fraction_monocytes", "mean"),
    ).reset_index()
    dex_summary.to_csv(OUT / "original_covid_dex_mono_summary.csv", index=False)
    sf = ifn_sf
    sf.to_csv(OUT / "original_covid_claim_frequencies_by_sample.csv", index=False)
    # per sample, by source deposited categories, frequency summary for IFN states / monocyte recapped dex state
    summary = sf.groupby(["population","case_control","WHO_severity","centre"], dropna=False).agg(
        n_samples=("sample_id","nunique"), median_ifn_fraction_parent=("ifn_fraction_parent","median"),
        mean_ifn_fraction_parent=("ifn_fraction_parent","mean"), samples_ifn_positive=("ifn_cells",lambda x:int((x>0).sum()))
    ).reset_index()
    summary.to_csv(OUT / "original_covid_claim_summary.csv", index=False)
    # Cohort-size audit documents why this local portal object is not yet an exact
    # published-cohort match (the paper reports 187 samples / 61 individuals).
    pd.DataFrame([
        {"cohort":"local source object after exclude==0", "cells":local_cells, "samples":local_samples,
         "patients":local_patients, "status":"full local roster before Table S1 filtering"},
        {"cohort":"published Table S1 roster after sample ID match", "cells":len(d), "samples":d.sample_id.nunique(),
         "patients":d.patient_id.nunique(), "status":"187/61 IDs; cell count can differ from paper due cell-level QC/provenance"},
        {"cohort":"paper reported analysis cohort", "cells":580040, "samples":187, "patients":61,
         "status":"article-reported counts"},
    ]).to_csv(OUT / "original_covid_cohort_size_audit.csv", index=False)
    return d, ifn_sf


if __name__ == "__main__":
    write_fixed_claims()
    obs, sepsis = sepsis_analysis()
    covid, covid_samples = covid_analysis()
    print(f"Sepsis cells={len(obs):,}; samples={obs['sample'].nunique()}; claims={len(SEPSIS_CLAIMS)}")
    print(f"COVID cells={len(covid):,}; samples={covid['sample_id'].nunique()}; patients={covid['patient_id'].nunique()}")
    print("Wrote frozen definition and original sample-level frequency tables under paper_claim_analysis/")
