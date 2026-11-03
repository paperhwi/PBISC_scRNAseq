"""C + D. Cell proportions by four mapping routes, with a reject option, and a confusion-matrix correction.

Same split and read-out emulation as report/check/e_checks.py (seed 7): reference = real cells of half A of the samples,
held-out = real cells of half B, shown untouched (realU) and through a PBISC-like read-out (realT2, known labels).
Generated = K = 1,000 and K = 4,000 cells (100 per sample per draw).

Routes
  R1 harmony_knn : pairwise Harmony (batch = source) on reference PCA, kNN (k = 30, 1/d weights) from reference cells
  R2 refpca_knn  : same kNN directly in the reference PCA space, no Harmony
  R3 t2_logreg   : multinomial logistic regression on the 2,000 reference HVGs, trained on reference cells passed through
                   the PBISC-like read-out (no integration; it learns the read-out format)
  R4 markers     : marker gating at the six-lineage level only (no reference cells used except for marker scaling)
Reject option: distance (R1, R2) or confidence (R3, R4) beyond the value that rejects 5 % of PBISC-like real cells.
Levels: lin6 (T, NK, B, Plasma, Myeloid, Other), major, fine.
D: confusion matrix P(pred | true) from realT2 per route / level; NNLS-corrected group proportions; per-label recall."""
import sys, json, time, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, harmonypy
from scipy.optimize import nnls
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
sys.path.insert(0, "/home/claude/work/report"); sys.path.insert(0, "/home/claude/work")
import store
coh = sys.argv[1]
U = "/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; R = "/home/claude/work/report"; O = "/home/claude/work/prop"
rng = np.random.default_rng(7); t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)


def lin6(x):
    x = str(x)
    if x in ("Plasma cells",) or x.startswith("B_ASC"): return "Plasma"
    if x == "B cells" or x.startswith("B_"): return "B"
    if x in ("CD4+ T cells", "CD8+ T cells", "NKT cells", "MAIT", "NKT", "T g/d") or x.startswith("T "): return "T"
    if x in ("NK cells", "NK") or x.startswith("NK "): return "NK"
    if x in ("Monocytes", "Macrophages", "DC", "Int.mono", "ASDC", "pDC") or x.startswith(("CD14mono", "CD16mono", "DC")): return "Myeloid"
    return "Other"


def load_real(c):
    if c == "sepsis":
        st = store.Store("/home/claude/work/sepsis/store"); sym = st.var.gene_symbol.astype(str).to_numpy()
        ro = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False,
                                    usecols=["sample", "celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme"]) for i in range(3)], ignore_index=True)
        o = pd.DataFrame({"sample": ro["sample"].astype(str), "fine": ro.celltype_ye2026, "major": ro.celltype_major_ye2026})
        o["group"] = st.obs.group.astype(str).values
        ok = ro.in_paper_scheme.astype(str).eq("True").to_numpy()
    else:
        st = store.Store("/home/claude/work/covid/store"); sym = st.var.gene.astype(str).to_numpy()
        L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False)
        o = pd.DataFrame({"sample": st.obs.sample_id.astype(str).values, "fine": L.cell_type_IFNstim_recap.values, "major": L.cell_type_detailed_recap.values})
        o["group"] = st.obs.case_control.astype(str).str.lower().map({"positive": "covid_positive", "negative": "covid_negative", "recovery": "recovery"}).values
        roster = set(pd.read_csv(f"{R}/covid_roster187_local_ids.csv").sample_id.astype(str))
        ok = ~np.isin(o.major, ["RBC", "Platelet"]) & o["sample"].isin(roster).to_numpy()
    o["lin6"] = o.major.map(lin6)
    return st, sym, o, ok


def sample_blocks(o, ok, samples, per):
    out = []
    for s in samples:
        ix = np.where(ok & (o["sample"].to_numpy() == s))[0]
        if len(ix) <= per: out.append(ix); continue
        j = rng.integers(0, len(ix) - per); out.append(ix[j:j + per])
    return np.sort(np.concatenate(out))


st, sym, ro, ok = load_real(coh)
samples = np.array(sorted(pd.unique(ro["sample"][ok]))); rng.shuffle(samples); half = len(samples) // 2
A, B = samples[:half], samples[half:]
per = 600 if coh == "sepsis" else 400
ia = sample_blocks(ro, ok, A, per); ib = sample_blocks(ro, ok, B, per)
log(coh, "ref cells", len(ia), "held-out cells", len(ib))
gen = {}
for K in ["k1000", "k4000"]:
    z = np.load(f"{U}/{coh}_{K}_counts.npz", allow_pickle=True)
    t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False)
    gen[K] = (sp.csr_matrix((z["data"].astype(np.float32), z["indices"].astype(np.int32), z["indptr"]), shape=(len(z["rows"]), len(z["var"]))),
              t, z["var"].astype(str))
gvar = gen["k1000"][2]
p1 = pd.Series(np.arange(len(sym)), index=sym); p1 = p1[~p1.index.duplicated()]
pg = pd.Series(np.arange(len(gvar)), index=gvar); pg = pg[~pg.index.duplicated()]
panel = list(dict.fromkeys(g for g in gvar if g in p1.index))
log("panel genes:", len(panel))
def read_panel(idx):
    return sp.vstack([st.rows(idx[i:i + 8000])[:, p1.loc[panel].to_numpy()] for i in range(0, len(idx), 8000)]).tocsr()
Xa = read_panel(ia); Xb = read_panel(ib); log("read real")
G = {K: gen[K][0][:, pg.loc[panel].to_numpy()].tocsr() for K in gen}
la, lb = ro.iloc[ia].reset_index(drop=True), ro.iloc[ib].reset_index(drop=True)


def lognorm(X):
    tot = np.asarray(X.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y


a = ad.AnnData(Xa.copy()); sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000); hv = np.where(a.var.highly_variable.to_numpy())[0]
Ya = lognorm(Xa); mu = np.asarray(Ya[:, hv].mean(0)).ravel(); sd = np.sqrt(np.asarray(Ya[:, hv].multiply(Ya[:, hv]).mean(0)).ravel() - mu ** 2) + 1e-6
def scale(Y): return np.clip((Y[:, hv].toarray() - mu) / sd, -10, 10).astype(np.float32)
pca = PCA(30, random_state=0).fit(scale(Ya)); Pa = pca.transform(scale(Ya)).astype(np.float32)
proj = lambda X: pca.transform(scale(lognorm(X))).astype(np.float32)

# ---------------- PBISC-like read-out (as e_checks realT2), applied to held-out B and to reference A (for R3 training)
gdepth = {}
for K in G:
    tot = np.asarray(G[K].sum(1)).ravel(); lin = gen[K][1].major.map(lin6).to_numpy()
    for c in np.unique(lin): gdepth.setdefault(c, []).append(tot[lin == c])
gdepth = {c: np.concatenate(v) for c, v in gdepth.items()}; alld = np.concatenate(list(gdepth.values()))


def t2(X, lin):
    tot = np.asarray(X.sum(1)).ravel()
    target = np.array([rng.choice(gdepth.get(c, alld)) for c in lin]).astype(np.float64)
    P = proj(X); _, nb = NearestNeighbors(n_neighbors=15, n_jobs=2).fit(P).kneighbors(P)
    prof = (sp.diags(1 / np.maximum(tot, 1)) @ X).tocsr()
    W = sp.csr_matrix((np.full(nb.size, 1 / 15, np.float32), nb.ravel(), np.arange(0, nb.size + 1, 15)), shape=(len(nb), len(nb)))
    sup = X.copy(); sup.data[:] = 1
    E = sp.vstack([(W[i:i + 2000] @ prof).multiply(sup[i:i + 2000]).tocsr() for i in range(0, W.shape[0], 2000)]).tocsr()
    E = (sp.diags(target / np.maximum(np.asarray(E.sum(1)).ravel(), 1e-12)) @ E).tocsr()
    T = E.copy(); T.data = rng.poisson(E.data).astype(np.float32); T.eliminate_zeros(); return T


T2b = t2(Xb, lb.lin6.to_numpy()); T2a = t2(Xa, la.lin6.to_numpy()); log("read-out emulation done")

LEVELS = ["lin6", "major", "fine"]
cats = {lv: pd.Categorical(la[lv]).categories for lv in LEVELS}


def knn_vote(ref, Q, labels):
    nn = NearestNeighbors(n_neighbors=30, n_jobs=2).fit(ref); d, idx = nn.kneighbors(Q); w = 1 / (d + 1e-3)
    out = {}
    for lv in LEVELS:
        c = pd.Categorical(labels[lv], categories=cats[lv]).codes[idx]; score = np.zeros((len(Q), len(cats[lv])))
        for j in range(idx.shape[1]): np.add.at(score, (np.arange(len(Q)), c[:, j]), w[:, j])
        out[lv] = np.asarray(cats[lv])[score.argmax(1)]; out[lv + "_conf"] = score.max(1) / score.sum(1)
    out["dist"] = d.mean(1); return out


# R3: logistic regression on read-out-emulated reference cells (one model per level)
Sa2 = scale(lognorm(T2a))
LR = {}
for lv in LEVELS:
    LR[lv] = LogisticRegression(C=0.05, max_iter=300, n_jobs=2).fit(Sa2, la[lv].to_numpy()); log("LR fitted", lv)
del Sa2

# R4: marker gating, six lineages; markers scaled by their mean in the matching reference lineage
MK = {"T": ["CD3E", "CD3D", "CD3G"], "NK": ["KLRF1", "NCAM1", "SH2D1B"], "B": ["MS4A1", "CD79A", "CD79B", "BANK1"],
      "Plasma": ["JCHAIN", "MZB1", "TNFRSF17"], "Myeloid": ["LYZ", "CST3", "S100A9", "FCN1", "FCER1A", "LILRA4"]}
pidx = {g: i for i, g in enumerate(panel)}
mref = {k: np.array([np.asarray(Ya[np.where(la.lin6 == k)[0]][:, pidx[g]].mean()).item() for g in v if g in pidx]) for k, v in MK.items()}


def markers(X):
    Y = lognorm(X); S = np.column_stack([np.asarray(Y[:, [pidx[g] for g in v if g in pidx]].toarray() / mref[k]).mean(1) for k, v in MK.items()])
    k = np.array(list(MK)); best = S.argmax(1); mx = S.max(1)
    return np.where(mx < 0.15, "Other", k[best]), mx


sources = {"realU": (Xb, lb), "realT2": (T2b, lb), "k1000": (G["k1000"], gen["k1000"][1]), "k4000": (G["k4000"], gen["k4000"][1])}
pred = {}
for name, (X, meta) in sources.items():
    P = proj(X); r = {}
    r["R2"] = knn_vote(Pa, P, la)
    pre = np.vstack([Pa, P]); s = np.r_[np.zeros(len(Pa), int), np.ones(len(P), int)].astype(str)
    ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"s": s}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); post = (Z.T if Z.shape[0] == pre.shape[1] else Z).astype(np.float32)
    r["R1"] = knn_vote(post[:len(Pa)], post[len(Pa):], la)
    Sx = scale(lognorm(X)); r["R3"] = {}
    for lv in LEVELS:
        pr = LR[lv].predict_proba(Sx); r["R3"][lv] = LR[lv].classes_[pr.argmax(1)]; r["R3"][lv + "_conf"] = pr.max(1)
    del Sx
    m, mx = markers(X); r["R4"] = {"lin6": m, "lin6_conf": mx}
    pred[name] = (r, meta, P, post[len(Pa):]); log("mapped", name)

# reject thresholds from realT2 (5 % of PBISC-like real cells rejected)
thr = {"R1": np.quantile(pred["realT2"][0]["R1"]["dist"], 0.95), "R2": np.quantile(pred["realT2"][0]["R2"]["dist"], 0.95),
       "R3": {lv: np.quantile(pred["realT2"][0]["R3"][lv + "_conf"], 0.05) for lv in LEVELS},
       "R4": np.quantile(pred["realT2"][0]["R4"]["lin6_conf"], 0.05)}


def rejected(route, r, lv):
    if route in ("R1", "R2"): return r["dist"] > thr[route]
    if route == "R3": return r[lv + "_conf"] < thr["R3"][lv]
    return r["lin6_conf"] < thr["R4"]


# ---------------- long table: per sample (and draw) counts of predicted labels, with reject flag
rows = []
for name, (r, meta, _, _) in pred.items():
    smp = meta["sample"].astype(str).to_numpy() if "sample" in meta else meta.sample_id.astype(str).to_numpy()
    drw = meta["draw"].astype(str).to_numpy() if "draw" in meta else np.full(len(meta), "real")
    grp = meta["group"].astype(str).to_numpy()
    for route in ["R1", "R2", "R3", "R4"]:
        for lv in (["lin6"] if route == "R4" else LEVELS):
            rej = rejected(route, r[route], lv)
            df = pd.DataFrame({"sample": smp, "draw": drw, "group": grp, "label": r[route][lv], "rejected": rej})
            g = df.groupby(["sample", "draw", "group", "label", "rejected"]).size().rename("n").reset_index()
            g["source"] = name; g["route"] = route; g["level"] = lv; rows.append(g)
    if name in ("realU", "realT2"):      # true labels of the held-out cells, same layout
        for lv in LEVELS:
            df = pd.DataFrame({"sample": smp, "draw": drw, "group": grp, "label": meta[lv].to_numpy(), "rejected": False})
            g = df.groupby(["sample", "draw", "group", "label", "rejected"]).size().rename("n").reset_index()
            g["source"] = name; g["route"] = "truth"; g["level"] = lv; rows.append(g)
long = pd.concat(rows, ignore_index=True); long.to_csv(f"{O}/{coh}_mapping_counts.csv.gz", index=False)

# ---------------- D: confusion from realT2 and corrected generated proportions
truthB = lb
conf_rows, corr_rows, acc = [], [], {}
for route in ["R1", "R2", "R3", "R4"]:
    for lv in (["lin6"] if route == "R4" else LEVELS):
        p = pred["realT2"][0][route][lv]; t = truthB[lv].to_numpy()
        C = pd.crosstab(pd.Categorical(t, categories=cats[lv]), pd.Categorical(p, categories=sorted(set(cats[lv]) | set(p))), dropna=False)
        Cn = C.div(C.sum(1).replace(0, np.nan), axis=0).fillna(0)
        acc[f"{route}_{lv}"] = {"acc_realT2": float((p == t).mean()), "acc_realU": float((pred["realU"][0][route][lv] == truthB[lv].to_numpy()).mean())}
        for lab in cats[lv]:
            conf_rows.append({"route": route, "level": lv, "label": lab, "n_true": int(C.loc[lab].sum()),
                              "recall": float(Cn.loc[lab, lab]) if lab in Cn.columns else 0.0,
                              "precision": float(C.loc[lab, lab] / max(C[lab].sum(), 1)) if lab in C.columns else np.nan})
        # correction: q_obs = Cn^T p_true  ->  NNLS for p_true (per source, per group)
        Mx = Cn.to_numpy().T                     # pred x true
        for name in ["realT2", "k1000", "k4000"]:
            r, meta = pred[name][0], pred[name][1]
            for grp_ in pd.unique(meta["group"].astype(str)):
                m = meta["group"].astype(str).to_numpy() == grp_
                q = pd.Series(r[route][lv][m]).value_counts(normalize=True).reindex(Cn.columns, fill_value=0).to_numpy()
                pt, _ = nnls(np.vstack([Mx, np.ones(Mx.shape[1]) * 10]), np.r_[q, 10]); pt = pt / pt.sum()
                obs = pd.Series(r[route][lv][m]).value_counts(normalize=True)
                for j, lab in enumerate(cats[lv]):
                    corr_rows.append({"route": route, "level": lv, "source": name, "group": grp_, "label": lab,
                                      "observed": float(obs.get(lab, 0.0)), "corrected": float(pt[j])})
pd.DataFrame(conf_rows).to_csv(f"{O}/{coh}_confusion_recall.csv", index=False)
pd.DataFrame(corr_rows).to_csv(f"{O}/{coh}_corrected_props.csv", index=False)
json.dump({"thresholds": {k: (v if not isinstance(v, dict) else v) for k, v in thr.items()}, "accuracy": acc,
           "n_ref": int(len(ia)), "n_heldout": int(len(ib)), "ref_samples": list(map(str, A)), "heldout_samples": list(map(str, B)),
           "panel_genes": len(panel)}, open(f"{O}/{coh}_mapping_meta.json", "w"), indent=1, default=float)
# per-cell predictions as h5ad (embeddings + labels; no expression) for inspection
for name in ["realT2", "k1000", "k4000"]:
    r, meta, P, post = pred[name]
    obs = pd.DataFrame({"source": name}, index=np.arange(len(P)).astype(str))
    for c in ["sample", "sample_id", "draw", "group", "lin6", "major", "fine", "label"]:
        if c in meta: obs[("true_" if name == "realT2" and c in LEVELS else "orig_") + c] = meta[c].astype(str).values
    for route in ["R1", "R2", "R3", "R4"]:
        for lv in (["lin6"] if route == "R4" else LEVELS):
            obs[f"{route}_{lv}"] = r[route][lv]; obs[f"{route}_{lv}_rejected"] = rejected(route, r[route], lv)
    aa = ad.AnnData(np.zeros((len(P), 0), np.float32), obs=obs); aa.obsm["X_refpca"] = P; aa.obsm["X_harmony_pair"] = post
    aa.write_h5ad(f"{O}/{coh}_{name}_mapping.h5ad", compression="gzip")
log("done", json.dumps(acc))
