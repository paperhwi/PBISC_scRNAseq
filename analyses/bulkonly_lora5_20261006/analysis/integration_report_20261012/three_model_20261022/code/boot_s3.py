"""(4) Patient-bootstrap CIs for donor identity and per-sample mixing (2026-10-22). Usage: python boot_s3.py
(a) COVID donor fingerprint (floor vs prev; definitions as val/code/pb_analysis.py part A, top-100 donor-stable genes): per-donor success
    (own donor = best match) and margin (own r - median r to other donors); donors resampled with replacement (B = 2,000), paired floor - prev.
(b) scSemiProfiler sample self-identification: scSemi regenerates real held-out samples from their own bulk, so its pseudobulk can be matched to
    the same real sample (sample level, all six types pooled, same donor-stable gene rule). A different (easier) task from (a): reported separately.
(c) Per-sample iLISI: per-cell LISI of source labels (perplexity 30, k = 90) for 12,000 generated cells in each pairwise Harmony space, median per
    generated sample; samples resampled (B = 2,000) -> CI of the median over samples."""
import numpy as np, pandas as pd, anndata as ad, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
from sklearn.neighbors import NearestNeighbors
P = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb"; P2 = "/home/claude/work/s3/data/pb"
BJ = "/mnt/user-data/outputs/ilisi_floor_20261021/data/big_joint"; OUT = "/mnt/user-data/outputs/s3/results"; NB = 2000; rng = np.random.default_rng(11)
def load(n):
    d = P2 if n.startswith("scsemi") else P; z = np.load(f"{d}/{n}.npz", allow_pickle=True); g = pd.Series(z["groups"]).str.split("\\|\\|", regex=True)
    return pd.DataFrame({"sample": g.str[0], "six": g.str[1], "n": z["n"]}), z["S"].astype(np.float64), np.array(z["genes"])
def cpm_log(M): return np.log2(M / np.maximum(M.sum(1, keepdims=True), 1) * 1e6 + 1)
def by_sample(meta, S):
    u = pd.Categorical(meta["sample"]); A = np.zeros((len(u.categories), S.shape[1])); np.add.at(A, u.codes, S[meta.index.to_numpy()]); return pd.Index(u.categories), A
out = []
# ---- (a)
mr, Sr, genes = load("real_covid"); ro = pd.read_csv(f"{P}/real_covid_obs.csv"); sid, Ar = by_sample(mr, Sr); Lr = cpm_log(Ar)
don_r = ro.set_index("sample").loc[sid, "patient"].to_numpy(); expr = (2 ** Lr - 1).mean(0) >= 10
cnt = pd.Series(don_r).value_counts(); rep = np.isin(don_r, cnt[cnt >= 2].index); X = Lr[rep][:, expr]; d = don_r[rep]; gm = X.mean(0)
dm = pd.DataFrame(X).groupby(d).transform("mean").to_numpy(); icc = ((dm - gm) ** 2).sum(0) / np.maximum(((X - gm) ** 2).sum(0), 1e-12)
cols = np.where(expr)[0][np.argsort(-icc)[:100]]
Dr = pd.DataFrame(Ar).groupby(don_r).sum(); real_d, real_A = Dr.index.to_numpy(), Dr.to_numpy()
z = lambda M: (M[:, cols] - M[:, cols].mean(0)) / (M[:, cols].std(0) + 1e-9)
per = {}
for src in ("floor", "prev"):
    m, S, _ = load(f"{src}_covid"); s, A = by_sample(m, S); dn = pd.Series(s).str.split("_").str[0].to_numpy(); Dg = pd.DataFrame(A).groupby(dn).sum()
    common = np.intersect1d(Dg.index, real_d); G = cpm_log(Dg.loc[common].to_numpy()); R = cpm_log(real_A[pd.Index(real_d).get_indexer(common)])
    C = np.corrcoef(z(G), z(R))[: len(common), len(common):]; rank = (C > np.diag(C)[:, None]).sum(1)
    marg = np.diag(C) - np.array([np.median(np.delete(C[i], i)) for i in range(len(common))]); per[src] = pd.DataFrame({"hit": rank == 0, "margin": marg}, index=common)
dd = per["floor"].join(per["prev"], lsuffix="_floor", rsuffix="_prev"); n = len(dd); bs = {k: [] for k in ("hit_floor", "hit_prev", "margin_floor", "margin_prev", "dhit", "dmargin")}
for _ in range(NB):
    i = rng.integers(0, n, n); b = dd.iloc[i]
    for k in ("hit_floor", "hit_prev", "margin_floor", "margin_prev"): bs[k].append(b[k].mean())
    bs["dhit"].append(b.hit_floor.mean() - b.hit_prev.mean()); bs["dmargin"].append(b.margin_floor.mean() - b.margin_prev.mean())
for k, v in bs.items():
    v = np.array(v); full = {"hit_floor": dd.hit_floor.mean(), "hit_prev": dd.hit_prev.mean(), "margin_floor": dd.margin_floor.mean(), "margin_prev": dd.margin_prev.mean(),
                             "dhit": dd.hit_floor.mean() - dd.hit_prev.mean(), "dmargin": dd.margin_floor.mean() - dd.margin_prev.mean()}[k]
    out.append({"analysis": "COVID donor fingerprint (n donors = %d, chance top-1 = %.3f)" % (n, 1 / n), "metric": k, "estimate": full, "lo": np.quantile(v, .025), "hi": np.quantile(v, .975), "p_le0": float((v <= 0).mean())})
# ---- (b) scSemi sample self-identification
for coh in ("covid", "sepsis"):
    mr, Sr, genes = load(f"real_{coh}"); sid, Ar = by_sample(mr, Sr); m, S, _ = load(f"scsemi_{coh}"); s, A = by_sample(m, S)
    common = np.intersect1d(s, sid); G = cpm_log(A[s.get_indexer(common)]); R = cpm_log(Ar[sid.get_indexer(common)]); Lr2 = cpm_log(Ar)
    if coh == "covid": cc = cols
    else:
        # sepsis has one sample per donor: no ICC; use the 100 most variable expressed genes across real samples instead
        ex = (2 ** Lr2 - 1).mean(0) >= 10; cc = np.where(ex)[0][np.argsort(-Lr2[:, ex].var(0))[:100]]
    zz = lambda M: (M[:, cc] - M[:, cc].mean(0)) / (M[:, cc].std(0) + 1e-9)
    C = np.corrcoef(zz(G), zz(R))[: len(common), len(common):]; rank = (C > np.diag(C)[:, None]).sum(1); hit = rank == 0
    v = np.array([hit[rng.integers(0, len(hit), len(hit))].mean() for _ in range(NB)])
    out.append({"analysis": f"scSemi sample self-identification, {coh} (n = {len(common)}, chance = {1/len(common):.3f})", "metric": "top1", "estimate": hit.mean(), "lo": np.quantile(v, .025), "hi": np.quantile(v, .975),
                "own_r_median": float(np.median(np.diag(C))), "other_r_median": float(np.median(C[~np.eye(len(common), dtype=bool)]))})
# ---- (c) per-sample iLISI
for coh in ("sepsis", "covid"):
    for src in ("realU", "floor", "prev4k", "scsemi"):
        J = ad.read_h5ad(f"{BJ}/{coh}_joint_{src}.h5ad"); E = J.obsm["X_pca_harmony"]; so = J.obs.source.astype(str).to_numpy(); smp = J.obs["sample"].astype(str).to_numpy()
        g = np.where(so != "ref")[0]; sub = rng.choice(g, min(12000, len(g)), replace=False)
        dist, ix = NearestNeighbors(n_neighbors=90, n_jobs=2).fit(E).kneighbors(E[sub]); v = L.lisi_vec(dist, ix, so)
        ps = pd.DataFrame({"s": smp[sub], "v": v}).groupby("s").v.median(); ps = ps[pd.DataFrame({"s": smp[sub]}).groupby("s").size().loc[ps.index] >= 20]
        bm = np.array([np.median(ps.to_numpy()[rng.integers(0, len(ps), len(ps))]) for _ in range(NB)])
        ps.rename("iLISI_median").to_frame().assign(cohort=coh, source=src).to_csv(f"{OUT}/{coh}_{src}_per_sample_ilisi.csv")
        out.append({"analysis": f"per-sample iLISI {coh}", "metric": src, "estimate": float(np.median(ps)), "lo": np.quantile(bm, .025), "hi": np.quantile(bm, .975),
                    "n_samples": int(len(ps)), "share_samples_iLISI_gt_1.5": float((ps > 1.5).mean()), "cell_level_median": float(np.median(v))})
        print(out[-1], flush=True)
O = pd.DataFrame(out); O.to_csv(f"{OUT}/bootstrap_identity_mixing.csv", index=False); pd.set_option("display.width", 250); print(O.round(3).to_string())
