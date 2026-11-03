"""Sample-level abundance comparison for paper-specific cell claims.

Uses obs only from backed AnnData objects; writes per-sample proportions for
each original cohort and each LoRA draw without treating cells as replicates.
"""
from pathlib import Path
import anndata as ad
import h5py
import pandas as pd

ROOT = Path(r"F:\PBISC")
OUT = ROOT / "paper_claim_analysis"
OUT.mkdir(exist_ok=True)
BASE = ROOT / "bulkonly_lora5_20261006"

def read_obs(path, cohort, dataset, draw=None):
    # Read only annotations. This avoids materializing the 9 GB float64 COVID
    # logcounts layer just to compute cell fractions.
    with h5py.File(path, "r") as h5:
        obs = ad.io.read_elem(h5["obs"])
    obs["cohort"] = cohort
    obs["dataset"] = dataset
    obs["draw"] = draw if draw is not None else "original"
    return obs.reset_index(drop=True)

def sepsis():
    orig_path = BASE / "analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad"
    orig = read_obs(orig_path, "sepsis", "original")
    gen = []
    for d in range(5):
        p = BASE / f"ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw{d}.h5ad"
        gen.append(read_obs(p, "sepsis", "LoRA", d))
    df = pd.concat([orig, *gen], ignore_index=True, sort=False)
    # Paper-defined subset abundance is available only for original labels;
    # map synthetic broad AIFI labels to the closest parent where possible.
    df["sample"] = df["sample"].astype(object).where(df["sample"].notna(), df["sample_id"].astype(object))
    df["donor"] = df["donor_id"].astype(object).where(df["donor_id"].notna(), df["sample"].astype(object))
    df["condition"] = df["group"].astype(str)
    if "age_group" in df:
        adult = df["age_group"].astype(str).eq("adult")
        df.loc[adult & df["group"].astype(str).eq("sepsis"), "condition"] = "adult_sepsis"
        df.loc[adult & df["group"].astype(str).str.contains("healthy", case=False), "condition"] = "adult_HC"
        df.loc[adult & df["group"].astype(str).str.contains("ICU", case=False), "condition"] = "adult_ICU_ctrl"
    labels = {
        "NR4A2_CD4_parent_or_cluster": ("celltype_ye2026", "T_CD4_c02-GPR183-NR4A2", "aifi_l2", "Memory CD4 T cell"),
        "CCL4hi_CD8_parent_or_cluster": ("celltype_ye2026", "T_CD8_c05-GNLY-CCL4hi", "aifi_l2", "Memory CD8 T cell"),
        "RETN_mono_parent_or_cluster": ("celltype_ye2026", "Mono_c03-CD14-RETN", "aifi_l2", "CD14 monocyte"),
        "HLAII_mono_parent_or_cluster": ("celltype_ye2026", "Mono_c01-CD14-HLA-DRA", "aifi_l2", "CD14 monocyte"),
        "MAIT": ("celltype_ye2026", "T_CD8_c07-TRAV1-2", "aifi_l2", "MAIT"),
    }
    out=[]
    for key,(orig_col,orig_val,gen_col,gen_val) in labels.items():
        for (cohort,dataset,draw,sample,donor,condition),g in df.groupby(["cohort","dataset","draw","sample","donor","condition"],dropna=False,observed=True):
            if pd.isna(sample): continue
            col,val=(orig_col,orig_val) if dataset=="original" else (gen_col,gen_val)
            mask=g[col].astype(str).eq(val)
            out.append(dict(cohort=cohort,dataset=dataset,draw=draw,sample=sample,donor=donor,condition=condition,claim=key,match_level="paper_subset" if dataset=="original" else "closest_parent",n_cells=len(g),fraction=float(mask.mean())))
    pd.DataFrame(out).to_csv(OUT/"sepsis_claim_abundance_by_sample.csv",index=False)

def covid():
    orig_path = BASE / "analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad"
    orig = read_obs(orig_path, "covid_eskd", "original")
    gen=[]
    for d in range(5):
        p=BASE/f"imperial_eskd/lora5_bulkonly_imperial_eskd_draw{d}.h5ad"
        gen.append(read_obs(p,"covid_eskd","LoRA",d))
    df=pd.concat([orig,*gen],ignore_index=True,sort=False)
    df["sample"] = df["sample_id"]
    df["donor"] = df["patient_id"].astype(object).where(df["patient_id"].notna(), df["donor_id"].astype(object))
    df["condition"] = df["case_control"].astype(object).where(df["case_control"].notna(), df["group"].astype(object)).astype(str).str.lower()
    df["condition"] = df["condition"].replace({"positive":"covid_positive","negative":"covid_negative","recovery":"recovery"})
    labels={
      "CD14_monocyte":("cell_type_detailed_recap","CD14mono","aifi_l2","CD14 monocyte"),
      "CD16_monocyte":("cell_type_detailed_recap","CD16mono","aifi_l2","CD16 monocyte"),
      "intermediate_monocyte":("cell_type_detailed_recap","Int.mono","aifi_l2","Intermediate monocyte"),
      "plasma_or_ASC":("cell_type_detailed_recap","B_ASC_IgA","aifi_l2","Plasma cell"),
      "MAIT":("cell_type_detailed_recap","MAIT","aifi_l2","MAIT"),
    }
    out=[]
    for key,(orig_col,orig_val,gen_col,gen_val) in labels.items():
        for (cohort,dataset,draw,sample,donor,condition),g in df.groupby(["cohort","dataset","draw","sample","donor","condition"],dropna=False,observed=True):
            if pd.isna(sample): continue
            col,val=(orig_col,orig_val) if dataset=="original" else (gen_col,gen_val)
            out.append(dict(cohort=cohort,dataset=dataset,draw=draw,sample=sample,donor=donor,condition=condition,claim=key,match_level="paper_subset" if dataset=="original" else "closest_parent",n_cells=len(g),fraction=float(g[col].astype(str).eq(val).mean())))
    pd.DataFrame(out).to_csv(OUT/"covid_claim_abundance_by_sample.csv",index=False)

if __name__=="__main__":
    sepsis(); covid()
    print("Wrote sample-level abundance tables to",OUT)
