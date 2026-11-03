"""Audit patient pseudobulk feasibility for prespecified rare fine subtypes."""
from pathlib import Path
import sys
import anndata as ad
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
RARE=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027"
sys.path.insert(0,str(RARE/"code"))
from rare_lib import SEPSIS,COVID
from deg_lib import COVID_EXCL
OUT=ROOT/"NewLora/analysis_20261102/rare_disease_deg"
OUT.mkdir(parents=True,exist_ok=True)
GRP={"sepsis":{"adult_sepsis":"case","adult_HC":"control","adult_ICU_ctrl":"ICU"},
     "covid":{"covid_positive":"case","covid_negative":"control"}}
sources={}
for coh,defs in (("sepsis",SEPSIS),("covid",COVID)):
    R=ad.read_h5ad(BASE/"data/big"/f"{coh}_ref.h5ad",backed="r")
    U=ad.read_h5ad(BASE/"data/big"/f"{coh}_realU.h5ad",backed="r")
    real=pd.concat([R.obs,U.obs],ignore_index=True)
    real=real[["sample","group","major","fine"]].copy()
    sources[(coh,"real")]=real
    R.file.close();U.file.close()
    for src in ("prev4k","floor"):
        tab=pd.read_csv(BASE/"results/big"/f"{coh}_{src}_transfer.csv.gz",usecols=["sample","group","major_tr","fine_tr"],keep_default_na=False)
        sources[(coh,src)]=tab.rename(columns={"major_tr":"major","fine_tr":"fine"})
    tab=pd.read_csv(ROOT/"NewLora/analysis_20261102"/coh/"ceiling_cell_labels.csv.gz",usecols=["sample","group","major_tr","fine_tr"],keep_default_na=False)
    sources[(coh,"ceiling")]=tab.rename(columns={"major_tr":"major","fine_tr":"fine"})

out=[];unit_rows=[]
for (coh,src),obs in sources.items():
    obs=obs.copy()
    if coh=="covid" and src!="real":obs=obs[~obs["sample"].isin(COVID_EXCL)]
    obs["condition"]=obs.group.map(GRP[coh])
    obs=obs[obs.condition.notna()].copy()
    obs["patient"]=obs["sample"] if coh=="sepsis" else obs["sample"].str.split("_").str[0]
    defs=SEPSIS if coh=="sepsis" else COVID
    m=obs.major.astype(str).to_numpy();f=obs.fine.astype(str).to_numpy()
    for sub,parfn in defs.items():
        mask=parfn(m,f)&(f==sub)
        count=obs.loc[mask].groupby(["patient","condition"],observed=True).size()
        all_units=obs.groupby(["patient","condition"],observed=True).size()
        for (patient,condition),n in all_units.items():
            unit_rows.append({"cohort":coh,"source":src,"subset":sub,"patient":patient,
                              "condition":condition,"n_subset":int(count.get((patient,condition),0)),
                              "n_all":int(n)})
        for cond in ("case","control"):
            k=count.xs(cond,level="condition") if cond in count.index.get_level_values("condition") else pd.Series(dtype=int)
            out.append({"cohort":coh,"source":src,"subset":sub,"condition":cond,
                        "n_cells_total":int(k.sum()),"n_patients_ge1":int((k>=1).sum()),
                        "n_patients_ge5":int((k>=5).sum()),"n_patients_ge10":int((k>=10).sum()),
                        "n_patients_ge20":int((k>=20).sum()),
                        "median_cells_detected":float(k.median()) if len(k) else 0})
pd.DataFrame(out).to_csv(OUT/"feasibility_by_condition.csv",index=False)
pd.DataFrame(unit_rows).to_csv(OUT/"rare_patient_counts.csv.gz",index=False,compression="gzip")
print(pd.DataFrame(out).pivot(index=["cohort","subset"],columns=["source","condition"],values="n_patients_ge20").fillna(0).astype(int).to_string())
