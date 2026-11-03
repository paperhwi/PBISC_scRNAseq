"""Which bulk samples were PBISC floor generated from, and how do they relate to the bulk cohort, the real scRNA samples and the adapter donors?"""
import pandas as pd, numpy as np
U="/mnt/user-data/uploads/new_inference_depth256_K4000_floor1e-4_20261021/"
PB="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/pb/"
ADAPT_COVID={"C126_positive_5","C126_positive_9","C126_positive_11","C141_negative","C58_positive_2"}
out=[]
# ---------------- sepsis ----------------
f=pd.read_csv(U+"sepsis/samples.csv"); b=pd.read_pickle("out2/sepsis_bulk.pkl").columns.tolist()
bm=pd.read_csv("out2/sepsis_bulk_meta.csv",index_col=0); r=pd.read_csv(PB+"real_sepsis_obs.csv")
fs=set(f.sample_id); bs=set(b); rs=set(r["sample"])
print("SEPSIS")
print(f"  floor samples {len(fs)} | in bulk GSE279448 {len(fs&bs)} / bulk total {len(bs)}")
print(f"  bulk samples that also have scRNA: {len(bs&rs)};  bulk-only: {len(bs-rs)}")
print(f"  floor ∩ real scRNA samples (sample = patient here): {len(fs&rs)}")
print(f"  bulk-only samples NOT generated: {sorted(bs-rs-fs)}")
print("  floor groups:",f.group.value_counts().to_dict())
f["cohort"]="sepsis"; f["in_bulk"]=f.sample_id.isin(bs); f["has_real_scRNA_same_sample"]=f.sample_id.isin(rs)
f["same_patient_has_scRNA"]=f.donor_id.isin(rs); f["adapter_training_donor"]=False
out.append(f)
# ---------------- covid ----------------
f=pd.read_csv(U+"covid/samples.csv"); cb=pd.read_pickle("out2/covid_bulk.pkl").columns.tolist()
from_meta=pd.read_csv("out2/covid_bulk_meta.csv",index_col=0)
import os
w=pd.concat([pd.read_csv("/mnt/user-data/uploads/bulk_public/zenodo6497251_imperial_eskd/w1_metadata.csv"),
             pd.read_csv("/mnt/user-data/uploads/bulk_public/zenodo6497251_imperial_eskd/w2_metadata.csv")]).drop_duplicates("sample_id").set_index("sample_id")
r=pd.read_csv(PB+"real_covid_obs.csv")
fs=set(f.sample_id); bs=set(cb); rp=set(r.patient.astype(str)); fp=set(f.donor_id.astype(str))
print("\nCOVID / ESKD")
print(f"  floor samples {len(fs)} ({len(fp)} patients) | in bulk Zenodo 6497251 {len(fs&bs)} / bulk total {len(bs)}")
print(f"  real scRNA: {r['sample'].nunique()} samples, {len(rp)} patients")
print(f"  floor patients that ALSO have real scRNA (other visits): {len(fp&rp)} / {len(fp)}")
print(f"  floor patients with no scRNA at all: {len(fp-rp)}")
print("  floor groups:",f.group.value_counts().to_dict())
print("  adapter-training samples inside floor:",sorted(fs&ADAPT_COVID))
# which bulk samples were NOT generated?
nb=sorted(bs-fs); print(f"  bulk samples not generated: {len(nb)}  (e.g. {nb[:6]})")
cc=w.reindex(nb).case_control.value_counts(dropna=False).to_dict(); print("   their case_control:",cc)
f["cohort"]="covid"; f["in_bulk"]=f.sample_id.isin(bs); f["has_real_scRNA_same_sample"]=False
f["same_patient_has_scRNA"]=f.donor_id.astype(str).isin(rp); f["adapter_training_donor"]=f.sample_id.isin(ADAPT_COVID)
f=f.join(w[["case_control","WHO_severity","time_from_first_positive_swab"]],on="sample_id")
out.append(f)
T=pd.concat(out,ignore_index=True)
T.to_csv("prov/floor_source_samples.csv",index=False)
print("\nsaved prov/floor_source_samples.csv", T.shape)
