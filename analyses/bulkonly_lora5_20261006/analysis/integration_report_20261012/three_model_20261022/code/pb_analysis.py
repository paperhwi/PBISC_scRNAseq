"""Pseudobulk-level validation (2026-10-21): (A) COVID donor fingerprint, (B) M6 disease effect per six-type.
Inputs data/pb/<source>_<coh>.npz (sample x six-type sums, all cells of each file) from pb_build.py.
(A) 32 COVID donors have generated cells (bulk-only samples) and real scRNA (other time points mostly). Donor-stable genes are chosen on real
    samples only (between-donor share of variance over donors with >= 2 real samples, ICC); donors' pseudobulks (all cells) are z-scored per
    gene within source; identification = correlation of each generated donor with each real donor; top-1 accuracy and rank percentile.
(B) per six-type, per-sample pseudobulk (>= 20 cells), log2(CPM + 1) on the panel; Welch t per gene for case vs control in each source;
    real vs generated logFC agreement (Pearson r, OLS slope gen ~ real, sign agreement and top-100 overlap among real FDR < 0.05 genes),
    split by real expression rank (<= 2,000 vs > 2,000 by mean CP10k in real)."""
import numpy as np, pandas as pd, json
from scipy import stats
P = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb"; OUT = "/mnt/user-data/outputs/val/results"
def load(n):
    z = np.load(f"{P}/{n}.npz", allow_pickle=True); g = pd.Series(z["groups"]).str.split("\\|\\|", regex=True)
    return pd.DataFrame({"sample": g.str[0], "six": g.str[1], "n": z["n"]}), z["S"].astype(np.float64), np.array(z["genes"])
def cpm_log(M): return np.log2(M / np.maximum(M.sum(1, keepdims=True), 1) * 1e6 + 1)
def by_sample(meta, S, types=None):
    m = meta if types is None else meta[meta.six.isin(types)]
    u = pd.Categorical(m["sample"]); A = np.zeros((len(u.categories), S.shape[1]))
    np.add.at(A, u.codes, S[m.index.to_numpy()]); n = np.bincount(u.codes, weights=m.n.to_numpy())
    return pd.Index(u.categories), A, n
res = {}
# ---------------- (A) COVID donor fingerprint
mr, Sr, genes = load("real_covid"); ro = pd.read_csv(f"{P}/real_covid_obs.csv")
sid, Ar, nr = by_sample(mr, Sr); Lr = cpm_log(Ar); don_r = ro.set_index("sample").loc[sid, "patient"].to_numpy()
expr = (2 ** Lr - 1).mean(0) >= 10
cnt = pd.Series(don_r).value_counts(); rep = np.isin(don_r, cnt[cnt >= 2].index)
X = Lr[rep][:, expr]; d = don_r[rep]; gm = X.mean(0)
dm = pd.DataFrame(X).groupby(d).transform("mean").to_numpy()
icc = ((dm - gm) ** 2).sum(0) / np.maximum(((X - gm) ** 2).sum(0), 1e-12)
gexp = genes[expr]; order = np.argsort(-icc)
fp = {"top100_icc_genes": gexp[order[:100]].tolist(), "icc_top100_median": float(np.median(icc[order[:100]])), "icc_all_median": float(np.median(icc))}
def donor_pb(src):
    m, S, _ = load(f"{src}_covid"); ob = pd.read_csv(f"{P}/{src}_covid_obs.csv")
    s, A, n = by_sample(m, S); dn = pd.Series(s).str.split("_").str[0].to_numpy()
    D = pd.DataFrame(A).groupby(dn).sum(); return D.index.to_numpy(), D.to_numpy()
Dr = pd.DataFrame(Ar).groupby(don_r).sum(); real_d, real_A = Dr.index.to_numpy(), Dr.to_numpy()
rows = []
for src in ("floor", "prev"):
    gd, gA = donor_pb(src); common = np.intersect1d(gd, real_d)
    G = cpm_log(gA[pd.Index(gd).get_indexer(common)]); R = cpm_log(real_A[pd.Index(real_d).get_indexer(common)])
    for gs_name, gsel in (("top-100 donor-stable genes", order[:100]), ("top-300 donor-stable genes", order[:300]), ("all expressed genes", np.arange(expr.sum()))):
        cols = np.where(expr)[0][gsel]
        z = lambda M: (M[:, cols] - M[:, cols].mean(0)) / (M[:, cols].std(0) + 1e-9)
        C = np.corrcoef(z(G), z(R))[: len(common), len(common):]
        rank = (C > np.diag(C)[:, None]).sum(1)          # 0 = own donor is best match
        rows.append({"source": src, "gene_set": gs_name, "n_donors": len(common), "top1_acc": float((rank == 0).mean()), "top3_acc": float((rank < 3).mean()),
                     "median_rank_percentile": float(np.median(1 - rank / (len(common) - 1))), "own_r_median": float(np.median(np.diag(C))),
                     "other_r_median": float(np.median(C[~np.eye(len(common), dtype=bool)]))})
    # sex genes: Y-gene detection per donor, real vs generated
    yg = [g for g in ("RPS4Y1", "DDX3Y", "KDM5D", "UTY") if g in genes]; yi = [list(genes).index(g) for g in yg]
    ysum = lambda A: np.log2(A[:, yi].sum(1) / A.sum(1) * 1e6 + 1)
    yr = ysum(real_A[pd.Index(real_d).get_indexer(common)]); yg_ = ysum(gA[pd.Index(gd).get_indexer(common)])
    male_r = yr > 4; rows[-1]["sex_agreement"] = float(((yg_ > np.median(yg_[male_r]) / 2) == male_r).mean()); rows[-1]["n_male_real"] = int(male_r.sum())
    rows[-1]["y_cpm_log_real_median_male_female"] = [float(np.median(yr[male_r])), float(np.median(yr[~male_r]))]
    rows[-1]["y_cpm_log_gen_median_male_female"] = [float(np.median(yg_[male_r])), float(np.median(yg_[~male_r]))]
FP = pd.DataFrame(rows); FP.to_csv(f"{OUT}/covid_donor_fingerprint.csv", index=False); res["fingerprint"] = fp
print(FP.drop(columns=[c for c in FP.columns if c.startswith("y_")]).round(3).to_string()); print("donor-stable genes (top 20):", fp["top100_icc_genes"][:20])
# ---------------- (B) M6
SIX = ["T cell", "NK cell", "B cell", "Monocyte", "DC"]
def groups(coh, src, samples):
    if src == "real" and coh == "sepsis":
        o = pd.read_csv(f"{P}/real_sepsis_obs.csv").set_index("sample"); g = o.loc[samples]
        return np.where(g.age_group.eq("adult") & g.group.eq("sepsis"), "case", np.where(g.age_group.eq("adult") & g.group.eq("healthy control"), "ctrl",
                        np.where(g.age_group.eq("adult") & g.group.eq("ICU control"), "icu", "other")))
    if src == "real":
        o = pd.read_csv(f"{P}/real_covid_obs.csv").set_index("sample"); g = o.loc[samples].case_control
        return np.where(g.eq("POSITIVE"), "case", np.where(g.eq("NEGATIVE"), "ctrl", "other"))
    o = pd.read_csv(f"{P}/{src}_{coh}_obs.csv").set_index("sample"); g = o.loc[samples].group
    if coh == "sepsis": return np.where(g.eq("adult_sepsis"), "case", np.where(g.eq("adult_HC"), "ctrl", np.where(g.eq("adult_ICU_ctrl"), "icu", "other")))
    excl = {"C126_positive_5", "C126_positive_9", "C126_positive_11", "C141_negative", "C58_positive_2"}   # adapter-training samples
    return np.where(np.isin(samples, list(excl)), "other", np.where(g.eq("covid_positive"), "case", np.where(g.eq("covid_negative"), "ctrl", "other")))
def de(coh, src, t, ctrl="ctrl"):
    m, S, genes = load(f"{src}_{coh}"); mm = m[(m.six == t) & (m.n >= 20)]
    smp = mm["sample"].to_numpy(); L = cpm_log(S[mm.index.to_numpy()]); g = groups(coh, src, smp)
    a, b = L[g == "case"], L[g == ctrl]
    t_, p = stats.ttest_ind(a, b, equal_var=False); return pd.DataFrame({"gene": genes, "logFC": a.mean(0) - b.mean(0), "t": t_, "p": p, "mean": L.mean(0)}), (len(a), len(b))
def bh(p): p = np.nan_to_num(p, nan=1.0); o = np.argsort(p); q = p[o] * len(p) / (np.arange(len(p)) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
rows, keep = [], {}
for coh in ("sepsis", "covid"):
    for contrast, ctrl in (("case vs control", "ctrl"),) + ((("sepsis vs ICU control", "icu"),) if coh == "sepsis" else ()):
        # real rank by overall mean CP10k (all cells, real)
        m0, S0, genes = load(f"real_{coh}"); tot = S0.sum(0); rank = pd.Series(tot).rank(ascending=False).to_numpy()
        for t in SIX:
            R, nR = de(coh, "real", t, ctrl); R["q"] = bh(R.p.to_numpy()); R["rank"] = rank
            for src in ("floor", "prev"):
                Gd, nG = de(coh, src, t, ctrl); keep[(coh, contrast, t, src)] = pd.concat([R.add_prefix("real_"), Gd.add_prefix("gen_")], axis=1)
                for band, bm in (("all", np.ones(len(R), bool)), ("rank <= 2,000", rank <= 2000), ("rank > 2,000", rank > 2000)):
                    ok = bm & (R["mean"].to_numpy() > 1) & np.isfinite(R.logFC) & np.isfinite(Gd.logFC)
                    sig = ok & (R.q.to_numpy() < 0.05)
                    x, y = R.logFC.to_numpy()[ok], Gd.logFC.to_numpy()[ok]
                    top_r = R[sig].assign(a=lambda d: d.t.abs()).nlargest(100, "a").gene; top_g = Gd[ok].assign(a=lambda d: d.t.abs()).nlargest(100, "a").gene
                    xs, ys = R.logFC.to_numpy()[sig], Gd.logFC.to_numpy()[sig]
                    rows.append({"cohort": coh, "contrast": contrast, "six_type": t, "source": src, "genes": band, "n_real_case_ctrl": f"{nR[0]}/{nR[1]}", "n_gen_case_ctrl": f"{nG[0]}/{nG[1]}",
                                 "n_genes": int(ok.sum()), "n_real_DEG": int(sig.sum()), "logFC_r": float(np.corrcoef(x, y)[0, 1]) if ok.sum() > 10 else np.nan,
                                 "slope_on_DEG": float(np.polyfit(xs, ys, 1)[0]) if sig.sum() > 10 else np.nan, "sign_agree_DEG": float((np.sign(xs) == np.sign(ys)).mean()) if sig.sum() else np.nan,
                                 "top100_overlap": float(len(set(top_r) & set(top_g)) / max(len(top_r), 1)) if len(top_r) else np.nan})
M6 = pd.DataFrame(rows); M6.to_csv(f"{OUT}/m6_disease_effect.csv", index=False)
pd.concat({"|".join(k): v for k, v in keep.items()}).rename_axis(["key", "i"]).reset_index().to_csv(f"{OUT}/m6_gene_tables.csv.gz", index=False); json.dump(res, open(f"{OUT}/pb_meta.json", "w"), indent=1)
pd.set_option("display.width", 250)
print(M6[M6.genes == "all"].pivot_table(index=["cohort", "contrast", "six_type"], columns="source", values=["logFC_r", "slope_on_DEG", "sign_agree_DEG", "top100_overlap"]).round(2).to_string())
print(M6.groupby(["cohort", "contrast", "genes", "source"])[["logFC_r", "slope_on_DEG", "sign_agree_DEG", "top100_overlap", "n_real_DEG"]].median().round(2).to_string())
