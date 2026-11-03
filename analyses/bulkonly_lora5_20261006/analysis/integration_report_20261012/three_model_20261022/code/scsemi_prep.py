"""scSemiProfiler bundle on the PBISC panel (2026-10-22). Usage: python scsemi_prep.py <sepsis|covid> [n=60000]
Source: scSemiProfiler strict_5pairs draw0 export (expected counts = Poisson rates, 1,500 rows per test sample, frozen scoring labels).
Fairness: the real held-out panel IS the scSemi test set, so scSemi cells are restricted to samples whose patient is NOT in real half A
(the reference) and, for sepsis, to adults. scSemi training (support) donors are not test samples by construction.
Counts: one seeded Poisson draw from the rates (seed 20261030); a floor-rule copy (x<1e-4 -> 0, x<1 -> 1, else floor(x)) is kept as layer 'floor_rule'
for a read-out sensitivity check. Genes joined to the panel by symbol; panel genes scSemi does not model (availability mask 0) stay 0.
Cells: equal per sample (stratified), seed 20261031. Labels: celltype_l1_pred (six-type) and aifi_l2 of the frozen scoring annotator.
Output: /home/claude/work/s3/data/<coh>_scsemi.h5ad (obs sample, patient, group, six_type, lineage, fine='', major='', aifi_l2, label_conf)."""
import sys, zipfile, time, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
coh = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 60000; t0 = time.time()
S = f"/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/scsemi_eval_20261022/stage_{coh}"
BIG = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/big"
PB = BIG.replace("/big", "/pb"); OUT = "/home/claude/work/s3/data"
LIN = {"T cell": "T", "NK cell": "NK", "B cell": "B", "Monocyte": "Myeloid", "DC": "Myeloid", "Other": "Other"}
def log(*a): print(f"[{time.time()-t0:5.0f}s]", coh, *a, flush=True)
z = zipfile.ZipFile(f"{S}/scsemi.npz")
small = np.load(f"{S}/scsemi.npz", allow_pickle=False)
sid, bi, sym, avail = small["sample_ids"].astype(str), small["bulk_index"], small["gene_symbols"].astype(str), small["output_availability_mask"].astype(bool)
row_sample = sid[bi]
R = ad.read_h5ad(f"{BIG}/{coh}_ref.h5ad", backed="r"); U = ad.read_h5ad(f"{BIG}/{coh}_realU.h5ad", backed="r"); panel = np.array(R.var_names)
A_pat = set(R.obs.patient.astype(str))
if coh == "sepsis":
    ro = pd.read_csv(f"{PB}/real_sepsis_obs.csv").set_index("sample")
    gmap = {"sepsis": "adult_sepsis", "healthy control": "adult_HC", "ICU control": "adult_ICU_ctrl"}
    meta = pd.DataFrame({"sample": sid}); meta["patient"] = meta["sample"]
    r2 = ro.reindex(sid); meta["adult"] = r2["age_group"].eq("adult").to_numpy(); meta["group"] = r2["group"].map(gmap).fillna("other").to_numpy()
    print("test samples missing from real obs (dropped):", [s for s in sid if s not in ro.index])
else:
    ro = pd.read_csv(f"{PB}/real_covid_obs.csv").set_index("sample")
    known = pd.concat([R.obs[["sample", "group"]], U.obs[["sample", "group"]]]).drop_duplicates().set_index("sample").group.astype(str)
    meta = pd.DataFrame({"sample": sid}); r2 = ro.reindex(sid); meta["patient"] = r2["patient"].astype(str).to_numpy(); meta["adult"] = r2["patient"].notna().to_numpy()
    cc = r2["case_control"].astype(str).to_numpy(); print("test samples missing from real obs (dropped):", [s for s in sid if s not in ro.index])
    # group names as in the generated files; learn case_control -> group from the samples already in the bundles
    lm = pd.crosstab(ro.loc[known.index, "case_control"].astype(str), known.to_numpy()).idxmax(1).to_dict()
    meta["group"] = [known.get(s, lm.get(c, "other")) for s, c in zip(sid, cc)]
meta["in_A"] = meta.patient.isin(A_pat); meta["in_B"] = meta["sample"].isin(set(U.obs["sample"].astype(str)))
keep_s = meta[(~meta.in_A) & meta.adult]["sample"].to_numpy()
log("test samples", len(sid), "kept (adult, patient not in A)", len(keep_s), "of which real B samples", int(meta.set_index("sample").loc[keep_s, "in_B"].sum()),
    meta.set_index("sample").loc[keep_s, "group"].value_counts().to_dict())
pool = np.where(np.isin(row_sample, keep_s))[0]
rng = np.random.default_rng(20261031); q = int(np.ceil(n / len(keep_s))); pick = []
for s in keep_s:
    ix = pool[row_sample[pool] == s]; pick.append(rng.permutation(ix)[:q])
pick = np.sort(np.concatenate(pick)); pick = np.sort(rng.permutation(pick)[:n]) if len(pick) > n else pick
# gene join by symbol (first occurrence), only modeled columns used
sp_ = pd.Series(np.arange(len(sym)), index=sym); sp_ = sp_[~sp_.index.duplicated()]
gcol = np.array([sp_.get(g, -1) for g in panel]); ok = gcol >= 0
log("panel genes", len(panel), "found in scSemi axis", int(ok.sum()))
# stream generated_counts.npy rows
f = z.open("generated_counts.npy"); ver = np.lib.format.read_magic(f)
shape, fort, dt = (np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0)(f); assert not fort and len(shape) == 2
ncol = shape[1]; rb = ncol * dt.itemsize; chunk = 5000; out_p, out_f = [], []; prng = np.random.default_rng(20261030); want = set(pick.tolist())
for s0 in range(0, shape[0], chunk):
    e = min(s0 + chunk, shape[0]); buf = np.frombuffer(f.read((e - s0) * rb), dtype=dt).reshape(e - s0, ncol)
    r = pick[(pick >= s0) & (pick < e)] - s0
    if len(r):
        M = np.zeros((len(r), len(panel)), np.float32); M[:, ok] = buf[r][:, gcol[ok]]
        # unmodeled (mask 0) columns -> 0
        am = avail[bi[r + s0]][:, gcol[ok]]; M[:, ok] = np.where(am, M[:, ok], 0)
        out_p.append(sp.csr_matrix(prng.poisson(np.clip(M, 0, None)).astype(np.float32)))
        Fm = np.where(M < 1e-4, 0, np.where(M < 1, 1, np.floor(M))).astype(np.float32); out_f.append(sp.csr_matrix(Fm))
    if s0 % 50000 == 0: log("rows", e, "/", shape[0])
Xp = sp.vstack(out_p).tocsr(); Xf = sp.vstack(out_f).tocsr()
lab = pd.read_parquet(f"{S}/labels.parquet").set_index("npz_row").loc[pick]
m = meta.set_index("sample").loc[row_sample[pick]]
obs = pd.DataFrame({"sample": row_sample[pick], "patient": m.patient.to_numpy(), "group": m.group.to_numpy(), "six_type": lab.celltype_l1_pred.astype(str).to_numpy(),
                    "fine": "", "major": "", "aifi_l2": lab.aifi_l2.astype(str).to_numpy(), "label_conf": lab.label_confidence.to_numpy(), "npz_row": pick,
                    "in_realB": m.in_B.to_numpy()})
obs["lineage"] = [LIN.get(s, "Other") for s in obs.six_type]; obs.index = [f"scsemi_{i}" for i in range(len(obs))]
a = ad.AnnData(Xp, obs=obs, var=pd.DataFrame(index=pd.Index(panel, name="gene"))); a.layers["floor_rule"] = Xf
a.var["modeled"] = ok & np.array([avail[:, c].any() if c >= 0 else False for c in gcol])
a.uns["source"] = "scsemi"; a.uns["cohort"] = coh
a.write_h5ad(f"{OUT}/{coh}_scsemi.h5ad", compression="gzip")
t = np.asarray(Xp.sum(1)).ravel()
log("saved", Xp.shape, "mean UMI", round(float(t.mean())), "median", float(np.median(t)), "six", obs.six_type.value_counts().to_dict(), "groups", obs.group.value_counts().to_dict())
