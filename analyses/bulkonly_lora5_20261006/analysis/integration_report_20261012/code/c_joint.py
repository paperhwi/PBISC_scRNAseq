"""C: co-embedding of real reference cells with generated K=1,000 and K=4,000 cells; label transfer; composition.
Shared feature space = generated 6k panel (measured genes) intersected with the real gene set; both sources normalised
on panel totals (log1p CP10k). PCA is fitted on the real reference (scaled with reference mean/sd) and the generated
cells are projected into it ('pre'); Harmony with batch = source then aligns them ('post'). Labels transferred with a
distance-weighted kNN (k = 30) from real cells in the post space."""
import sys, json, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, "/home/claude/work/report"); sys.path.insert(0, "/home/claude/work")
import metrics as M, store
coh = sys.argv[1]
U = "/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; R = "/home/claude/work/report"
rng = np.random.default_rng(0)
# ---------------- real reference
if coh == "sepsis":
    st = store.Store("/home/claude/work/sepsis/store"); sym = st.var.gene_symbol.astype(str).to_numpy()
    ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False) for i in range(3)], ignore_index=True)
    ro["label"] = ro.celltype_ye2026; ro["major"] = ro.celltype_major_ye2026
    ro["grp"] = ro.age_group.str[:5].replace({"paedi": "ped"}) + "_" + ro.group.replace({"healthy control": "HC", "ICU control": "ICU_ctrl"})
    ok = ro.in_paper_scheme.astype(str).eq("True").to_numpy(); per = 400; bkey = ro["sample"]
else:
    st = store.Store("/home/claude/work/covid/store"); sym = st.var.gene.astype(str).to_numpy()
    L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
    ro = st.obs[["case_control", "WHO_severity_simple", "sample_id", "patient_id"]].reset_index(drop=True).copy()
    ro["label"] = L.cell_type_IFNstim_recap.values; ro["major"] = L.cell_type_detailed_recap.values
    ro["grp"] = ro.case_control.astype(str).str.lower().map({"positive": "covid_positive", "negative": "covid_negative", "recovery": "recovery"})
    roster = set(pd.read_csv(f"{R}/covid_roster187_local_ids.csv").sample_id.astype(str))   # published Table S1 roster (187 samples / 61 patients)
    ok = ro["sample_id"].astype(str).isin(roster).to_numpy(); per = 300; bkey = ro["sample_id"].astype(str)
ref_idx = np.sort(np.concatenate([rng.permutation(np.where(ok & (bkey == s).to_numpy())[0])[:per] for s in pd.unique(bkey[ok])]))
# ---------------- generated
gen = {}
for K in ["k1000", "k4000"]:
    z = np.load(f"{U}/{coh}_{K}_counts.npz", allow_pickle=True); o = pd.read_csv(f"{U}/{coh}_{K}_obs.csv.gz").iloc[z["rows"]].reset_index(drop=True)
    X = sp.csr_matrix((z["data"].astype(np.float32), z["indices"].astype(np.int32), z["indptr"]), shape=(len(z["rows"]), len(z["var"])))
    gen[K] = (X, o, z["var"].astype(str))
gvar = gen["k1000"][2]; assert (gvar == gen["k4000"][2]).all()
sym_pos = pd.Series(np.arange(len(sym)), index=sym); sym_pos = sym_pos[~sym_pos.index.duplicated()]
panel = [g for g in gvar if g in sym_pos.index]; gi_gen = np.array([np.where(gvar == g)[0][0] for g in panel]); gi_ref = sym_pos.loc[panel].to_numpy()
print(coh, "panel genes shared:", len(panel), "of", len(gvar), flush=True)
def lognorm(X):
    X = X.tocsr().astype(np.float32); tot = np.asarray(X.sum(1)).ravel()
    Y = sp.diags(1e4 / np.maximum(tot, 1)) @ X; Y = Y.tocsr(); Y.data = np.log1p(Y.data); return Y, tot
Xr_counts = st.rows(ref_idx)[:, gi_ref]
Xr, _ = lognorm(Xr_counts)
a = ad.AnnData(Xr_counts.tocsr()); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000)
hv = np.where(a.var.highly_variable.to_numpy())[0]
mu = np.asarray(Xr[:, hv].mean(0)).ravel(); sd = np.sqrt(np.asarray(Xr[:, hv].multiply(Xr[:, hv]).mean(0)).ravel() - mu ** 2) + 1e-6
def scale(Y): return np.clip((Y[:, hv].toarray() - mu) / sd, -10, 10)
pca = PCA(30, random_state=0).fit(scale(Xr))
parts = [("real", pca.transform(scale(Xr)), ro.iloc[ref_idx].reset_index(drop=True))]
for K in ["k1000", "k4000"]:
    Yg, _ = lognorm(gen[K][0][:, gi_gen]); parts.append((K, pca.transform(scale(Yg)), gen[K][1]))
src = np.concatenate([[s] * len(p) for s, p, _ in parts]); pre = np.vstack([p for _, p, _ in parts]).astype(np.float32)
post = M.__dict__.get("harm") or None
import harmonypy
ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"src": src}), ["src"], random_state=0, max_iter_harmony=20, verbose=False)
Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
mix = {f"iLISI_source_{t}": M.lisi(X, {"s": src}, n_sub=12000)["s"] for t, X in [("pre", pre), ("post", post)]}
# per-K pairwise mixing (real vs one K)
for K in ["k1000", "k4000"]:
    m = np.isin(src, ["real", K])
    for t, X in [("pre", pre), ("post", post)]:
        mix[f"iLISI_real_vs_{K}_{t}"] = M.lisi(X[m], {"s": src[m]}, n_sub=10000)["s"]
print(coh, mix, flush=True)
# label transfer
rlab = parts[0][2]["label"].to_numpy(); rmaj = parts[0][2]["major"].to_numpy(); nr = len(rlab)
nn = NearestNeighbors(n_neighbors=30, n_jobs=2).fit(post[:nr]); out = {"mix": mix}
transfer = {}
for K in ["k1000", "k4000"]:
    m = src == K; d, idx = nn.kneighbors(post[m]); w = 1 / (d + 1e-3)
    labs = pd.Categorical(rlab); codes = labs.codes[idx]
    score = np.zeros((m.sum(), len(labs.categories)))
    for j in range(idx.shape[1]): np.add.at(score, (np.arange(m.sum()), codes[:, j]), w[:, j])
    best = score.argmax(1); conf = score.max(1) / score.sum(1)
    t = pd.DataFrame({"label": labs.categories[best], "conf": conf})
    t["major"] = pd.Series(rmaj, index=rlab).groupby(level=0).first().reindex(t.label).values
    transfer[K] = pd.concat([gen[K][1].reset_index(drop=True), t], axis=1)
    transfer[K].to_csv(f"{R}/{coh}_{K}_transfer.csv.gz", index=False)
# real composition from the FULL reference (not the subsample)
ref_full = ro[ok].copy()
ref_full["sample"] = bkey[ok].values
ref_full.to_csv(f"{R}/{coh}_real_labels.csv.gz", index=False, columns=["grp", "label", "major", "sample"])
# UMAP of the joint space (post) on a balanced subsample for display
take = np.concatenate([rng.permutation(np.where(src == s)[0])[:40000] for s in ["real", "k1000", "k4000"]])
aa = ad.AnnData(np.zeros((len(take), 0), np.float32)); aa.obsm["X"] = post[take]
sc.pp.neighbors(aa, use_rep="X", n_neighbors=30, random_state=0); sc.tl.umap(aa, random_state=0, min_dist=0.3)
lab_all = np.concatenate([rlab, transfer["k1000"].label.values, transfer["k4000"].label.values])
maj_all = np.concatenate([rmaj, transfer["k1000"].major.values, transfer["k4000"].major.values])
np.savez_compressed(f"{R}/{coh}_joint_umap.npz", umap=aa.obsm["X_umap"], src=src[take], label=lab_all[take].astype(str), major=maj_all[take].astype(str))
json.dump({"mix": mix, "panel_genes": len(panel), "n_real_sub": int(nr), "n_real_samples": int(pd.unique(bkey[ok]).size)}, open(f"{R}/{coh}_joint.json", "w"), indent=1)
if coh == "sepsis":   # the COVID pre-Harmony UMAP exceeds the 8 GB container; metrics above are unaffected
    pre_take = ad.AnnData(np.zeros((len(take), 0), np.float32)); pre_take.obsm["X"] = pre[take]
    sc.pp.neighbors(pre_take, use_rep="X", n_neighbors=30, random_state=0); sc.tl.umap(pre_take, random_state=0, min_dist=0.3)
    np.save(f"{R}/{coh}_joint_umap_pre.npy", pre_take.obsm["X_umap"])
print(coh, "done", flush=True)
