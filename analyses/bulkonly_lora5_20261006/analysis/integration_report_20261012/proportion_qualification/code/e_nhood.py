"""E. Label-free neighbourhood analysis (Milo-style, implemented directly; no labels are used for testing).

Joint space as in c_joint.py: panel genes, 2,000 reference HVGs, PCA 30 fitted on real cells, generated cells projected,
Harmony with source as batch (real + one K at a time). kNN graph k = 30 on all cells.
Neighbourhoods: 2 % of cells sampled, each moved to the cell nearest the mean position of its 30 neighbours (Milo
"refined" sampling, duplicates dropped); neighbourhood = that cell and its 299 nearest neighbours (size ~ units, as
Milo advises for many samples).
Units: real samples; generated samples (five draws pooled, 500 cells each).
E1 abundance  : on a balanced subset (generated cells subsampled per clinical group to the real count; index cells drawn
                equally from both sources), per-unit fraction of the unit's cells in each neighbourhood, generated vs
                real samples (log2 ratio of means, Mann-Whitney, BH FDR).
E2 disease    : case vs control log2 ratio computed separately from real units and from generated units, per neighbourhood;
                compared across neighbourhoods that hold >= 30 cells of each source."""
import sys, json, time, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, harmonypy
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, "/home/claude/work"); import store
coh = sys.argv[1]
U = "/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; R = "/home/claude/work/report"; O = "/home/claude/work/prop"
rng = np.random.default_rng(11); t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)


if coh == "sepsis":
    st = store.Store("/home/claude/work/sepsis/store"); sym = st.var.gene_symbol.astype(str).to_numpy()
    ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                                usecols=["sample", "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme", "age_group", "group"]) for i in range(3)], ignore_index=True)
    ro["grp"] = ro.age_group.str[:5].replace({"paedi": "ped"}) + "_" + ro.group.replace({"healthy control": "HC", "ICU control": "ICU_ctrl"})
    ro["label"] = ro.celltype_ye2026; ro["major"] = ro.celltype_major_ye2026
    ok = ro.in_paper_scheme.astype(str).eq("True").to_numpy() & ro.grp.str.startswith("adult").to_numpy()
    CONTRASTS = [("adult_sepsis", "adult_HC"), ("adult_sepsis", "adult_ICU_ctrl")]; GROUPS = ["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]
else:
    st = store.Store("/home/claude/work/covid/store"); sym = st.var.gene.astype(str).to_numpy()
    L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
    ro = pd.DataFrame({"sample": st.obs.sample_id.astype(str).values, "label": L.cell_type_IFNstim_recap.values, "major": L.cell_type_detailed_recap.values})
    ro["grp"] = st.obs.case_control.astype(str).str.lower().map({"positive": "covid_positive", "negative": "covid_negative", "recovery": "recovery"}).values
    roster = set(pd.read_csv(f"{R}/covid_roster187_local_ids.csv").sample_id.astype(str))
    ok = ro["sample"].isin(roster).to_numpy() & ro.grp.isin(["covid_positive", "covid_negative"]).to_numpy()
    CONTRASTS = [("covid_positive", "covid_negative")]; GROUPS = ["covid_negative", "covid_positive"]
per = 300; NH = 300; out = []
for s in pd.unique(ro["sample"][ok]):
    ix = np.where(ok & (ro["sample"].to_numpy() == s))[0]
    if len(ix) <= per: out.append(ix); continue
    j = rng.integers(0, len(ix) - per); out.append(ix[j:j + per])
ir = np.sort(np.concatenate(out)); lr = ro.iloc[ir].reset_index(drop=True)
log("real cells", len(ir), "samples", lr["sample"].nunique())

zk = {}
KS = sys.argv[2:] or ["k1000", "k4000"]
for K in KS:
    z = np.load(f"{U}/{coh}_{K}_counts.npz", allow_pickle=True)
    t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False)
    zk[K] = (sp.csr_matrix((z["data"].astype(np.float32), z["indices"].astype(np.int32), z["indptr"]), shape=(len(z["rows"]), len(z["var"]))), t, z["var"].astype(str))
gvar = zk[KS[0]][2]
p1 = pd.Series(np.arange(len(sym)), index=sym); p1 = p1[~p1.index.duplicated()]
pg = pd.Series(np.arange(len(gvar)), index=gvar); pg = pg[~pg.index.duplicated()]
panel = list(dict.fromkeys(g for g in gvar if g in p1.index))
Xr = sp.vstack([st.rows(ir[i:i + 8000])[:, p1.loc[panel].to_numpy()] for i in range(0, len(ir), 8000)]).tocsr(); log("read real")


def lognorm(X):
    tot = np.asarray(X.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y


a = ad.AnnData(Xr.copy()); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000); hv = np.where(a.var.highly_variable.to_numpy())[0]
Yr = lognorm(Xr); mu = np.asarray(Yr[:, hv].mean(0)).ravel(); sd = np.sqrt(np.asarray(Yr[:, hv].multiply(Yr[:, hv]).mean(0)).ravel() - mu ** 2) + 1e-6
def scale(Y): return np.clip((Y[:, hv].toarray() - mu) / sd, -10, 10).astype(np.float32)
pca = PCA(30, random_state=0).fit(scale(Yr)); Pr = pca.transform(scale(Yr)).astype(np.float32)
import os
res = json.load(open(f"{O}/{coh}_nhood_summary.json")) if os.path.exists(f"{O}/{coh}_nhood_summary.json") else {}
for K in KS:
    X, t, _ = zk.pop(K); Pg = pca.transform(scale(lognorm(X[:, pg.loc[panel].to_numpy()]))).astype(np.float32)
    pre = np.vstack([Pr, Pg]); src = np.r_[np.full(len(Pr), "real"), np.full(len(Pg), "gen")]
    ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"s": src}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
    meta = pd.DataFrame({"source": src, "sample": np.r_[lr["sample"].astype(str), t.sample_id.astype(str)],
                         "draw": np.r_[np.full(len(Pr), "real"), t.draw.astype(str)], "grp": np.r_[lr.grp.astype(str), t.group.astype(str)],
                         "label": np.r_[lr.label.astype(str), t.label.astype(str)], "major": np.r_[lr.major.astype(str), t.major.astype(str)]})
    meta = meta[meta.grp.isin(GROUPS)].copy(); keep = meta.index.to_numpy(); post_k = post[keep]; meta = meta.reset_index(drop=True)
    nn = NearestNeighbors(n_neighbors=31, n_jobs=2).fit(post_k); _, knn = nn.kneighbors(post_k)
    # refined index sampling
    init = rng.choice(len(post_k), size=int(0.02 * len(post_k)), replace=False)
    cen = post_k[knn[init]].mean(1); _, ref_ix = nn.kneighbors(cen, n_neighbors=1); index = np.unique(ref_ix.ravel())
    _, nh = nn.kneighbors(post_k[index], n_neighbors=NH)                # n_nhood x NH cell indices
    unit = (meta.source + "|" + meta["sample"]).to_numpy(); ucat = pd.Categorical(unit); ucodes = ucat.codes
    cnt = np.zeros((len(index), len(ucat.categories)), np.float32)
    for j in range(nh.shape[1]): np.add.at(cnt, (np.arange(len(index)), ucodes[nh[:, j]]), 1)
    utot = np.bincount(ucodes, minlength=len(ucat.categories)).astype(np.float32)
    frac = cnt / utot[None, :]
    umeta = pd.DataFrame({"unit": ucat.categories}); umeta["source"] = umeta.unit.str.split("|").str[0]; umeta["sample"] = umeta.unit.str.split("|").str[1]
    umeta["grp"] = umeta["sample"].map(meta.drop_duplicates("sample").set_index("sample").grp)
    # nhood annotation: majority real label / major, majority transferred label of generated cells, source composition
    nhd = pd.DataFrame({"index_cell": index})
    srcm = meta.source.to_numpy()[nh]; nhd["n_real"] = (srcm == "real").sum(1); nhd["n_gen"] = (srcm == "gen").sum(1)
    for col in ["major", "label"]:
        labm = meta[col].to_numpy()[nh]
        nhd[f"real_{col}"] = [pd.Series(l[s == "real"]).mode().iloc[0] if (s == "real").any() else "(no real cells)" for l, s in zip(labm, srcm)]
        nhd[f"gen_{col}"] = [pd.Series(l[s == "gen"]).mode().iloc[0] if (s == "gen").any() else "(no generated cells)" for l, s in zip(labm, srcm)]
    eps = 1e-4
    # E1 (balanced): generated cells subsampled per group to the real cell count of that group; neighbourhood index
    # cells drawn equally from both sources, so coverage does not favour the larger source. Unit = sample.
    keepb = []
    for g in GROUPS:
        ir_ = np.where((meta.source == "real") & (meta.grp == g))[0]; ig_ = np.where((meta.source == "gen") & (meta.grp == g))[0]
        keepb += [ir_, rng.choice(ig_, size=min(len(ir_), len(ig_)), replace=False)]
    keepb = np.sort(np.concatenate(keepb)); mb = meta.iloc[keepb].reset_index(drop=True); pb_ = post_k[keepb]
    nnb = NearestNeighbors(n_neighbors=31, n_jobs=2).fit(pb_); _, knnb = nnb.kneighbors(pb_)
    nidx = int(0.02 * len(mb) / 2)
    init = np.r_[rng.choice(np.where(mb.source == "real")[0], nidx, replace=False), rng.choice(np.where(mb.source == "gen")[0], nidx, replace=False)]
    _, rix = nnb.kneighbors(pb_[knnb[init]].mean(1), n_neighbors=1); indexb = np.unique(rix.ravel())
    _, nhb = nnb.kneighbors(pb_[indexb], n_neighbors=NH)
    ub = pd.Categorical((mb.source + "|" + mb["sample"]).to_numpy()); cb_ = np.zeros((len(indexb), len(ub.categories)), np.float32)
    for j in range(nhb.shape[1]): np.add.at(cb_, (np.arange(len(indexb)), ub.codes[nhb[:, j]]), 1)
    fb = cb_ / np.bincount(ub.codes, minlength=len(ub.categories)).astype(np.float32)[None, :]
    isg = np.array([c.startswith("gen|") for c in ub.categories])
    e1 = pd.DataFrame({"index_cell": indexb})
    sm = mb.source.to_numpy()[nhb]; e1["n_real"] = (sm == "real").sum(1); e1["n_gen"] = (sm == "gen").sum(1)
    for col in ["major", "label"]:
        lm = mb[col].to_numpy()[nhb]
        e1[f"real_{col}"] = [pd.Series(l[s_ == "real"]).mode().iloc[0] if (s_ == "real").any() else "(no real cells)" for l, s_ in zip(lm, sm)]
    e1["lfc_gen_vs_real"] = np.log2((fb[:, isg].mean(1) + eps) / (fb[:, ~isg].mean(1) + eps))
    p = np.array([mannwhitneyu(fb[i, isg], fb[i, ~isg]).pvalue for i in range(len(indexb))])
    e1["fdr_gen_vs_real"] = multipletests(p, method="fdr_bh")[1]
    e1.to_csv(f"{O}/{coh}_{K}_nhoods_balanced.csv.gz", index=False)
    # E2: disease effect within each source
    for case, ctrl in CONTRASTS:
        for s in ["real", "gen"]:
            uc = np.where((umeta.source == s) & (umeta.grp == case))[0]; uk = np.where((umeta.source == s) & (umeta.grp == ctrl))[0]
            fc, fk = frac[:, uc], frac[:, uk]
            nhd[f"lfc_{case}_vs_{ctrl}_{s}"] = np.log2((fc.mean(1) + eps) / (fk.mean(1) + eps))
            p = np.array([mannwhitneyu(fc[i], fk[i]).pvalue if (fc[i].sum() + fk[i].sum()) > 0 else 1.0 for i in range(len(index))])
            nhd[f"fdr_{case}_vs_{ctrl}_{s}"] = multipletests(p, method="fdr_bh")[1]
    nhd.to_csv(f"{O}/{coh}_{K}_nhoods.csv.gz", index=False)
    # summaries
    r = {"n_cells": int(len(meta)), "n_nhoods": int(len(index)), "units_real": int((umeta.source == "real").sum()), "units_gen": int((umeta.source == "gen").sum())}
    both = (nhd.n_real >= 30) & (nhd.n_gen >= 30); r["share_nhoods_both_sources"] = float(both.mean())
    r["share_nhoods_real_only"] = float((nhd.n_gen == 0).mean()); r["share_nhoods_gen_only"] = float((nhd.n_real == 0).mean())
    sig = e1.fdr_gen_vs_real < 0.1
    r["E1_n_nhoods"] = int(len(e1)); r["E1_enriched_gen"] = float((sig & (e1.lfc_gen_vs_real > 1)).mean())
    r["E1_depleted_gen"] = float((sig & (e1.lfc_gen_vs_real < -1)).mean()); r["E1_median_gen_share"] = float((e1.n_gen / NH).median())
    r["E1_gen_only"] = float((e1.n_real == 0).mean()); r["E1_real_only"] = float((e1.n_gen == 0).mean())
    for case, ctrl in CONTRASTS:
        a_, b_ = nhd[f"lfc_{case}_vs_{ctrl}_real"], nhd[f"lfc_{case}_vs_{ctrl}_gen"]
        m = both.to_numpy(); r[f"{case}_vs_{ctrl}_spearman"] = float(spearmanr(a_[m], b_[m])[0])
        sigr = m & (nhd[f"fdr_{case}_vs_{ctrl}_real"] < 0.1).to_numpy()
        r[f"{case}_vs_{ctrl}_n_sig_real"] = int(sigr.sum())
        r[f"{case}_vs_{ctrl}_sign_agree_sig_real"] = float((np.sign(a_[sigr]) == np.sign(b_[sigr])).mean()) if sigr.sum() else None
        x, y = a_[m].to_numpy(), b_[m].to_numpy(); r[f"{case}_vs_{ctrl}_slope_gen_on_real"] = float(np.polyfit(x, y, 1)[0])
        r[f"{case}_vs_{ctrl}_slope_sig_real"] = float(np.polyfit(a_[sigr], b_[sigr], 1)[0]) if sigr.sum() > 5 else None
    res[K] = r; log(K, json.dumps(r)); json.dump(res, open(f"{O}/{coh}_nhood_summary.json", "w"), indent=1)
    import os
    if K == "k1000" and not os.path.exists(f"{O}/{coh}_k1000_nhood_embedding.h5ad"):
        aa = ad.AnnData(np.zeros((len(meta), 0), np.float32), obs=meta.astype(str)); aa.obsm["X_harmony"] = post_k
        sub = rng.choice(len(meta), size=min(60000, len(meta)), replace=False); sub.sort()
        b = aa[sub].copy(); sc.pp.neighbors(b, use_rep="X_harmony", n_neighbors=30, random_state=0); sc.tl.umap(b, random_state=0, min_dist=0.3)
        b.uns["nhood_index_cells_in_full"] = index
        ui = np.full(len(meta), -1); ui[sub] = np.arange(len(sub)); b.uns["nhood_index_cells_sub"] = ui[index]
        for k_ in list(b.obsp.keys()): del b.obsp[k_]
        b.write_h5ad(f"{O}/{coh}_k1000_nhood_embedding.h5ad", compression="gzip"); log("umap done")
json.dump(res, open(f"{O}/{coh}_nhood_summary.json", "w"), indent=1)
log("done")
