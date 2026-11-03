"""Calibration checks for the real-vs-generated domain gap (same joint pipeline as c_joint.py).

Reference (ref)  = real cells from half A of the cohort's samples.
Held-out real    = real cells from the other half B of samples, in four forms:
   realU  untouched
   realT1 depth-matched: binomial thinning to a UMI total drawn from the generated cells of the same lineage
   realT2 output-emulated: expected counts = kNN-smoothed (k = 15) profile of the cell, restricted to the cell's own
          detected genes (its 'support'), scaled to the matched depth, then one Poisson draw (PBISC's read-out format)
Generated        = K = 1,000 and K = 4,000 cells (100 per sample per draw), as in c_joint.py
Other study      = real cells of the other cohort (sepsis <-> COVID): check 1, cross-study real-vs-real baseline

Feature space: generated panel genes present in BOTH real datasets; log CP10k on panel totals; 2,000 HVG (seurat_v3, ref);
PCA 30 fitted on ref. For each source X: pairwise Harmony (batch = source) on ref + X, then
  iLISI (ref vs X), global and within coarse lineage (T, NK, B, myeloid), pre and post Harmony;
  kNN label transfer (k = 30) from ref -> accuracy where X has true labels; median confidence."""
import sys, json, time, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, harmonypy
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, "/home/claude/work/report"); sys.path.insert(0, "/home/claude/work")
import metrics as M, store
coh = sys.argv[1]; other = "covid" if coh == "sepsis" else "sepsis"
U = "/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; R = "/home/claude/work/report"; O = f"{R}/check"
rng = np.random.default_rng(7); t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)


def coarse(x):
    x = str(x)
    if x in ("B cells", "Plasma cells") or x.startswith("B_"): return "B"
    if x in ("CD4+ T cells", "CD8+ T cells", "NKT cells", "MAIT", "NKT") or x.startswith("T "): return "T"
    if x == "NK cells" or x == "NK" or x.startswith("NK "): return "NK"
    if x in ("Monocytes", "Macrophages", "DC", "Int.mono", "ASDC", "pDC") or x.startswith(("CD14mono", "CD16mono", "DC")): return "Myeloid"
    return "Other"


def load_real(c):
    if c == "sepsis":
        st = store.Store("/home/claude/work/sepsis/store"); sym = st.var.gene_symbol.astype(str).to_numpy()
        ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                                    usecols=["sample", "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme"]) for i in range(3)], ignore_index=True)
        o = pd.DataFrame({"sample": ro["sample"].astype(str), "fine": ro.celltype_ye2026, "major": ro.celltype_major_ye2026})
        ok = ro.in_paper_scheme.astype(str).eq("True").to_numpy()
    else:
        st = store.Store("/home/claude/work/covid/store"); sym = st.var.gene.astype(str).to_numpy()
        L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
        o = pd.DataFrame({"sample": st.obs.sample_id.astype(str).values, "fine": L.cell_type_IFNstim_recap.values, "major": L.cell_type_detailed_recap.values})
        roster = set(pd.read_csv(f"{R}/covid_roster187_local_ids.csv").sample_id.astype(str))   # published 187-sample roster
        ok = ~np.isin(o.major, ["RBC", "Platelet"]) & o["sample"].isin(roster).to_numpy()
    o["coarse"] = o.major.map(coarse)
    return st, sym, o, ok


def sample_blocks(o, ok, samples, per):
    """per cells from each sample as one contiguous run of that sample's eligible rows (fast memmap reads)."""
    out = []
    for s in samples:
        ix = np.where(ok & (o["sample"].to_numpy() == s))[0]
        if len(ix) <= per: out.append(ix); continue
        j = rng.integers(0, len(ix) - per); out.append(ix[j:j + per])
    return np.sort(np.concatenate(out))


# ---------------- data
st, sym, ro, ok = load_real(coh); st2, sym2, ro2, ok2 = load_real(other)
samples = np.array(sorted(pd.unique(ro["sample"][ok]))); rng.shuffle(samples); half = len(samples) // 2
A, B = samples[:half], samples[half:]
per = 600 if coh == "sepsis" else 400
ia = sample_blocks(ro, ok, A, per); ib = sample_blocks(ro, ok, B, per)
s2 = np.array(sorted(pd.unique(ro2["sample"][ok2]))); i2 = sample_blocks(ro2, ok2, s2, 40000 // len(s2) + 1)
log(coh, "ref", len(ia), "held-out", len(ib), other, len(i2))
gen = {}
for K in ["k1000", "k4000"]:
    z = np.load(f"{U}/{coh}_{K}_counts.npz", allow_pickle=True)
    t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False)
    gen[K] = (sp.csr_matrix((z["data"].astype(np.float32), z["indices"].astype(np.int32), z["indptr"]), shape=(len(z["rows"]), len(z["var"]))),
              t.major.map(coarse).to_numpy(), z["var"].astype(str))
gvar = gen["k1000"][2]
p1 = pd.Series(np.arange(len(sym)), index=sym); p1 = p1[~p1.index.duplicated()]
p2 = pd.Series(np.arange(len(sym2)), index=sym2); p2 = p2[~p2.index.duplicated()]
pg = pd.Series(np.arange(len(gvar)), index=gvar); pg = pg[~pg.index.duplicated()]   # 4-6 gene names occur twice in the panel
panel = list(dict.fromkeys(g for g in gvar if g in p1.index and g in p2.index))
log("panel genes in generated panel and both real sets:", len(panel), "of", len(gvar))
def read_panel(S_, idx, cols):
    return sp.vstack([S_.rows(idx[i:i + 8000])[:, cols] for i in range(0, len(idx), 8000)]).tocsr()
Xa = read_panel(st, ia, p1.loc[panel].to_numpy()); log("read ref")
Xb = read_panel(st, ib, p1.loc[panel].to_numpy()); log("read held-out")
X2 = read_panel(st2, i2, p2.loc[panel].to_numpy()); log("read other study")
G = {K: gen[K][0][:, pg.loc[panel].to_numpy()].tocsr() for K in gen}
la, lb, l2 = ro.iloc[ia].reset_index(drop=True), ro.iloc[ib].reset_index(drop=True), ro2.iloc[i2].reset_index(drop=True)


def lognorm(X):
    tot = np.asarray(X.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y


a = ad.AnnData(Xa.copy()); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000); hv = np.where(a.var.highly_variable.to_numpy())[0]
Ya = lognorm(Xa); mu = np.asarray(Ya[:, hv].mean(0)).ravel(); sd = np.sqrt(np.asarray(Ya[:, hv].multiply(Ya[:, hv]).mean(0)).ravel() - mu ** 2) + 1e-6
def scale(Y): return np.clip((Y[:, hv].toarray() - mu) / sd, -10, 10)
pca = PCA(30, random_state=0).fit(scale(Ya)); Pa = pca.transform(scale(Ya)).astype(np.float32)
proj = lambda X: pca.transform(scale(lognorm(X))).astype(np.float32)

# ---------------- transforms of the held-out real cells
gdepth = {}
for K in G:
    tot = np.asarray(G[K].sum(1)).ravel()
    for c in np.unique(gen[K][1]): gdepth.setdefault(c, []).append(tot[gen[K][1] == c])
gdepth = {c: np.concatenate(v) for c, v in gdepth.items()}
alld = np.concatenate(list(gdepth.values()))
cb = lb.coarse.to_numpy(); totb = np.asarray(Xb.sum(1)).ravel()
target = np.array([rng.choice(gdepth.get(c, alld)) for c in cb]).astype(np.float64)
# T1: binomial thinning (cells already shallower than the target are kept as they are)
keep = np.minimum(1.0, target / np.maximum(totb, 1))
T1 = Xb.copy(); T1.data = rng.binomial(T1.data.astype(np.int64), np.repeat(keep, np.diff(T1.indptr))).astype(np.float32); T1.eliminate_zeros()
# T2: expected counts from kNN-smoothed profile on the cell's own support, scaled to target depth, one Poisson draw
Pb = proj(Xb); nn = NearestNeighbors(n_neighbors=15, n_jobs=2).fit(Pb); _, nb = nn.kneighbors(Pb)
prof = (sp.diags(1 / np.maximum(totb, 1)) @ Xb).tocsr()
W = sp.csr_matrix((np.full(nb.size, 1 / 15, np.float32), nb.ravel(), np.arange(0, nb.size + 1, 15)), shape=(len(nb), len(nb)))
sup = Xb.copy(); sup.data[:] = 1                                # the cell's own detected genes
E = sp.vstack([(W[i:i + 2000] @ prof).multiply(sup[i:i + 2000]).tocsr() for i in range(0, W.shape[0], 2000)]).tocsr()  # smoothed profile on own support
rs = np.asarray(E.sum(1)).ravel()
E = (sp.diags(target / np.maximum(rs, 1e-12)) @ E).tocsr()
T2 = E.copy(); T2.data = rng.poisson(E.data).astype(np.float32); T2.eliminate_zeros()
log("transforms done")


def lisi_full(Y, lab, n_sub=10000, k=90):
    sub = rng.choice(len(Y), size=min(n_sub, len(Y)), replace=False)
    d, ix = NearestNeighbors(n_neighbors=k, n_jobs=2).fit(Y).kneighbors(Y[sub])
    v = M._lisi_from_knn(d, ix, np.asarray(lab))
    own = (np.asarray(lab)[ix] == np.asarray(lab)[sub][:, None]).mean(1)       # share of the 90 neighbours from the cell's own source
    return {"median": float(np.median(v)), "mean": float(np.mean(v)), "own_share": float(own.mean())}


def stats(X):
    tot = np.asarray(X.sum(1)).ravel(); det = np.diff(X.tocsr().indptr)
    return {"umi_median": float(np.median(tot)), "detected_median": float(np.median(det)), "frac_ones": float((X.data == 1).mean())}


sources = {"realU": (Xb, cb, lb), "realT1": (T1, cb, lb), "realT2": (T2, cb, lb), other: (X2, l2.coarse.to_numpy(), l2),
           "k1000": (G["k1000"], gen["k1000"][1], None), "k4000": (G["k4000"], gen["k4000"][1], None)}
res = {"panel_genes": len(panel), "n_ref": int(len(ia)), "ref_samples": len(A), "heldout_samples": len(B),
       "stats": {"ref": stats(Xa), **{k: stats(v[0]) for k, v in sources.items()}}}
log(json.dumps(res["stats"]))
ca = la.coarse.to_numpy(); nA = len(Pa)
umap_store = {}
for name, (X, cX, lX) in sources.items():
    P = proj(X); pre = np.vstack([Pa, P]); src = np.r_[np.zeros(nA, int), np.ones(len(P), int)]
    ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"s": src.astype(str)}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
    lin = np.r_[ca, cX]; r = {}
    for tag, Y in [("pre", pre), ("post", post)]:
        f = lisi_full(Y, src); r[f"iLISI_{tag}"] = f["median"]; r[f"iLISI_mean_{tag}"] = f["mean"]; r[f"own_share_{tag}"] = f["own_share"]
        for c in ["T", "NK", "B", "Myeloid"]:
            m = lin == c
            if (m & (src == 0)).sum() > 500 and (m & (src == 1)).sum() > 500:
                r[f"iLISI_{tag}_{c}"] = M.lisi(Y[m], {"s": src[m]}, n_sub=5000)["s"]
    # label transfer from ref in the post space
    nnr = NearestNeighbors(n_neighbors=30, n_jobs=2).fit(post[:nA]); d, idx = nnr.kneighbors(post[nA:]); w = 1 / (d + 1e-3)
    for col in ["coarse", "major", "fine"]:
        cats = pd.Categorical(la[col]); codes = cats.codes[idx]; score = np.zeros((len(P), len(cats.categories)))
        for j in range(idx.shape[1]): np.add.at(score, (np.arange(len(P)), codes[:, j]), w[:, j])
        pred = np.asarray(cats.categories)[score.argmax(1)]; conf = score.max(1) / score.sum(1)
        r[f"conf_{col}_median"] = float(np.median(conf))
        if lX is not None and (col == "coarse" or name != other):
            r[f"acc_{col}"] = float((pred == lX[col].to_numpy()).mean())
    res[name] = r; log(name, json.dumps({k: round(v, 3) for k, v in r.items()}))
    if name in ("realU", "realT2", "k1000", other):
        take = np.r_[rng.permutation(nA)[:15000], nA + rng.permutation(len(P))[:15000]]
        umap_store[name] = (post[take], src[take], lin[take])
# c_joint-style variant: one Harmony over ref + K1000 + K4000, then ref-vs-K mixing on the pair
P1, P4 = proj(G["k1000"]), proj(G["k4000"]); pre = np.vstack([Pa, P1, P4]); src3 = np.r_[np.zeros(nA, int), np.ones(len(P1), int), np.full(len(P4), 2)]
ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"s": src3.astype(str)}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
for K, code in [("k1000", 1), ("k4000", 2)]:
    m = np.isin(src3, [0, code]); f = lisi_full(post[m], src3[m]); g = lisi_full(pre[m], src3[m])
    res[f"{K}_3source_harmony"] = {"iLISI_post": f["median"], "iLISI_mean_post": f["mean"], "own_share_post": f["own_share"], "iLISI_pre": g["median"], "iLISI_mean_pre": g["mean"]}
log("3-source", json.dumps({k: v for k, v in res.items() if "3source" in k}))
json.dump(res, open(f"{O}/{coh}_checks.json", "w"), indent=1)
# UMAPs for display (post-Harmony, ref vs source)
for name, (Y, s, l) in umap_store.items():
    aa = ad.AnnData(np.zeros((len(Y), 0), np.float32)); aa.obsm["X"] = Y
    sc.pp.neighbors(aa, use_rep="X", n_neighbors=30, random_state=0); sc.tl.umap(aa, random_state=0, min_dist=0.3)
    np.savez_compressed(f"{O}/{coh}_umap_{name}.npz", umap=aa.obsm["X_umap"], src=s, lin=l)
log("done")
