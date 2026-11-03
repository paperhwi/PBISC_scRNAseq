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

# ---- save inputs for the monocyte / claim-classifier analyses (identical cells and read-out draws as c_mapping.py)
import os
O2 = "/home/claude/work/mono"; os.makedirs(O2, exist_ok=True)
obs_r = pd.concat([la.assign(half="A"), lb.assign(half="B")], ignore_index=True).astype(str)
Xr = sp.vstack([Xa, Xb]).tocsr(); T2 = sp.vstack([T2a, T2b]).tocsr()
r = ad.AnnData(Xr, obs=obs_r, var=pd.DataFrame(index=panel)); r.layers["t2"] = T2
r.obs_names = [f"real{i}" for i in range(r.n_obs)]; r.write_h5ad(f"{O2}/{coh}_real_AB.h5ad", compression="gzip")
gl = []
for K in ["k1000", "k4000"]:
    m = ad.read_h5ad(f"{O}/{coh}_{K}_mapping.h5ad", backed="r").obs
    t = gen[K][1]
    ob = pd.DataFrame({"K": K, "sample": m.orig_sample_id.astype(str).values, "draw": m.orig_draw.astype(str).values,
                       "group": m.orig_group.astype(str).values, "R1_fine": m.R1_fine.astype(str).values, "R3_fine": m.R3_fine.astype(str).values,
                       "R3_major": m.R3_major.astype(str).values, "R3_lin6": m.R3_lin6.astype(str).values})
    g = ad.AnnData(G[K], obs=ob, var=pd.DataFrame(index=panel)); g.obs_names = [f"{K}_{i}" for i in range(g.n_obs)]; gl.append(g)
ad.concat(gl).write_h5ad(f"{O2}/{coh}_gen.h5ad", compression="gzip")
log("saved", Xr.shape, T2.shape)
