"""scSemiProfiler pseudobulk (2026-10-22): sample x six-type sums of the expected counts over ALL 1,500 rows of every test sample
(rates, no Poisson draw: a pseudobulk sum of Poisson draws has the rate sum as its expectation), on the pseudobulk gene axis of data/pb.
Also six-type and aifi_l2 composition per sample. Output /home/claude/work/s3/data/pb/scsemi_<coh>.npz (+ _obs.csv with in_A flag, _comp_six.csv).
Usage: python scsemi_pb.py <coh>"""
import sys, zipfile, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
coh = sys.argv[1]
S = f"/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/scsemi_eval_20261022/stage_{coh}"
PB = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb"
OUT = "/home/claude/work/s3/data/pb"; import os; os.makedirs(OUT, exist_ok=True)
genes = np.load(f"{PB}/real_{coh}.npz", allow_pickle=True)["genes"]
sm = np.load(f"{S}/scsemi.npz", allow_pickle=False); sid = sm["sample_ids"].astype(str); bi = sm["bulk_index"]; sym = sm["gene_symbols"].astype(str)
avail = sm["output_availability_mask"].astype(bool)
s_ = pd.Series(np.arange(len(sym)), index=sym); s_ = s_[~s_.index.duplicated()]; gcol = np.array([s_.get(g, -1) for g in genes]); ok = gcol >= 0
lab = pd.read_parquet(f"{S}/labels.parquet").sort_values("npz_row"); six = lab.celltype_l1_pred.astype(str).to_numpy()
key = pd.Categorical(np.char.add(np.char.add(sid[bi], "||"), six)); K = len(key.categories); Ssum = np.zeros((K, len(genes))); n = np.bincount(key.codes, minlength=K)
z = zipfile.ZipFile(f"{S}/scsemi.npz"); f = z.open("generated_counts.npy"); ver = np.lib.format.read_magic(f)
shape, fort, dt = (np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0)(f); rb = shape[1] * dt.itemsize
for s0 in range(0, shape[0], 5000):
    e = min(s0 + 5000, shape[0]); buf = np.frombuffer(f.read((e - s0) * rb), dtype=dt).reshape(e - s0, shape[1])
    M = np.zeros((e - s0, len(genes))); M[:, ok] = np.where(avail[bi[s0:e]][:, gcol[ok]], buf[:, gcol[ok]], 0)
    Ssum += sp.csr_matrix((np.ones(e - s0), (key.codes[s0:e], np.arange(e - s0))), shape=(K, e - s0)) @ M
np.savez_compressed(f"{OUT}/scsemi_{coh}.npz", S=Ssum.astype(np.float32), n=n, groups=np.array(key.categories, dtype=object), genes=genes)
R = ad.read_h5ad(f"/home/claude/work/s3/big/{coh}_ref.h5ad", backed="r"); A_pat = set(R.obs.patient.astype(str))
if coh == "sepsis":
    ro = pd.read_csv(f"{PB}/real_sepsis_obs.csv").set_index("sample").reindex(np.unique(sid))
    gm = {"sepsis": "adult_sepsis", "healthy control": "adult_HC", "ICU control": "adult_ICU_ctrl"}
    grp = np.where(ro.age_group.eq("adult"), ro.group.map(gm).fillna("other"), "paediatric"); grp = np.where(ro.group.isna(), "other", grp); pat = ro.index.to_numpy()
else:
    ro = pd.read_csv(f"{PB}/real_covid_obs.csv").set_index("sample").reindex(np.unique(sid))
    grp = ro.case_control.map({"POSITIVE": "covid_positive", "NEGATIVE": "covid_negative"}).fillna("recovery").to_numpy(); grp = np.where(ro.case_control.isna(), "other", grp)
    pat = ro.patient.astype(str).to_numpy()
obs = pd.DataFrame({"sample": ro.index, "group": grp, "patient": pat, "in_A": np.isin(pat, list(A_pat))}); obs.to_csv(f"{OUT}/scsemi_{coh}_obs.csv", index=False)
comp = pd.DataFrame({"sample": sid[bi], "six": six}).groupby(["sample", "six"]).size().rename("n").reset_index(); comp.to_csv(f"{OUT}/scsemi_{coh}_comp_six.csv", index=False)
fine = pd.DataFrame({"sample": sid[bi], "aifi_l2": lab.aifi_l2.astype(str).to_numpy()}).groupby(["sample", "aifi_l2"]).size().rename("n").reset_index(); fine.to_csv(f"{OUT}/scsemi_{coh}_comp_aifi_l2.csv", index=False)
print(coh, Ssum.shape, "genes found", int(ok.sum()), "/", len(genes), obs.groupby(["group", "in_A"]).size().to_dict())
