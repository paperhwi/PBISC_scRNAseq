"""M6 disease effect, three models + bulk-broadcast baseline + patient bootstrap (2026-10-22). Usage: python m6_s3.py
Same definitions as val/code/pb_analysis.py part (B): per six-type, per-sample pseudobulk (>= 20 cells), log2(CPM + 1) on the panel, Welch t case vs
control; agreement with real (all real samples) = Pearson r of logFC (genes with real mean > 1), OLS slope and sign agreement on real FDR < 0.05 genes,
top-100 overlap; bands by real expression rank.
Sources: floor (depth256-K4000-floor), prev (Poisson K4000), scsemi (scSemiProfiler strict_5pairs draw0, adult test samples whose patient is not in
real half A; rates summed over all 1,500 rows), scsemi_all (all adult test samples: scSemi regenerates the very samples the truth is computed on),
broadcast (bulk RNA-seq logFC of the bulk-only samples the floor model was given, the same value for every cell type: sepsis GSE279448 121 samples;
COVID Imperial ESKD, Gisby 2022 Zenodo 6497251 htseq counts, 164 samples, same case/control rule and exclusions as floor).
Partial r of each model given broadcast. Bootstrap (B = 300): patients resampled within case and control, real and generated drawn independently,
the same real draw shared by all sources (paired differences); DEG set and gene masks fixed from the full data."""
import numpy as np, pandas as pd, glob, gzip, json, zlib
from scipy import stats
P = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb"; P2 = "/home/claude/work/s3/data/pb"
OUT = "/mnt/user-data/outputs/s3/results"; SIX = ["T cell", "NK cell", "B cell", "Monocyte", "DC"]; NB = 300
def load(n):
    d = P2 if n.startswith("scsemi") else P; z = np.load(f"{d}/{n}.npz", allow_pickle=True); g = pd.Series(z["groups"]).str.split("\\|\\|", regex=True)
    return pd.DataFrame({"sample": g.str[0], "six": g.str[1], "n": z["n"]}), z["S"].astype(np.float64), np.array(z["genes"])
def cpm_log(M): return np.log2(M / np.maximum(M.sum(1, keepdims=True), 1) * 1e6 + 1)
EXCL = {"C126_positive_5", "C126_positive_9", "C126_positive_11", "C141_negative", "C58_positive_2"}
def groups(coh, src, samples):
    """returns (group label case/ctrl/icu/other, patient id)"""
    if src == "real" and coh == "sepsis":
        o = pd.read_csv(f"{P}/real_sepsis_obs.csv").set_index("sample").loc[samples]
        g = np.where(o.age_group.eq("adult") & o.group.eq("sepsis"), "case", np.where(o.age_group.eq("adult") & o.group.eq("healthy control"), "ctrl",
                     np.where(o.age_group.eq("adult") & o.group.eq("ICU control"), "icu", "other"))); return g, np.asarray(samples)
    if src == "real":
        o = pd.read_csv(f"{P}/real_covid_obs.csv").set_index("sample").loc[samples]
        return np.where(o.case_control.eq("POSITIVE"), "case", np.where(o.case_control.eq("NEGATIVE"), "ctrl", "other")), o.patient.astype(str).to_numpy()
    if src.startswith("scsemi"):
        o = pd.read_csv(f"{P2}/scsemi_{coh}_obs.csv").set_index("sample").loc[samples]; g = o.group.astype(str)
        g = np.where(g.str.endswith("sepsis") | g.eq("covid_positive"), "case", np.where(g.eq("adult_HC") | g.eq("covid_negative"), "ctrl", np.where(g.eq("adult_ICU_ctrl"), "icu", "other")))
        if src == "scsemi": g = np.where(o.in_A.to_numpy(), "other", g)
        return g, o.patient.astype(str).to_numpy()
    o = pd.read_csv(f"{P}/{src}_{coh}_obs.csv").set_index("sample").loc[samples]; g = o.group
    if coh == "sepsis":
        return np.where(g.eq("adult_sepsis"), "case", np.where(g.eq("adult_HC"), "ctrl", np.where(g.eq("adult_ICU_ctrl"), "icu", "other"))), np.asarray(samples)
    return (np.where(np.isin(samples, list(EXCL)), "other", np.where(g.eq("covid_positive"), "case", np.where(g.eq("covid_negative"), "ctrl", "other"))),
            pd.Series(samples).str.split("_").str[0].to_numpy())
CACHE = {}
def per_sample(coh, src, t):
    k = (coh, src, t)
    if k not in CACHE:
        m, S, genes = load(f"{'scsemi' if src.startswith('scsemi') else src}_{coh}"); mm = m[(m.six == t) & (m.n >= 20)]
        smp = mm["sample"].to_numpy(); g, pat = groups(coh, src, smp); CACHE[k] = (cpm_log(S[mm.index.to_numpy()]), g, pat, genes)
    return CACHE[k]
def bulk_sepsis(genes):
    o = pd.read_csv(f"{P}/floor_sepsis_obs.csv").set_index("sample")
    sm = np.load("/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/scsemi_eval_20261022/stage_sepsis/scsemi.npz")
    e2s = dict(zip(sm["gene_ids"].astype(str), sm["gene_symbols"].astype(str))); cols = pd.Series(np.arange(len(genes)), index=genes)
    rows, names = [], []
    for f in sorted(glob.glob("/home/claude/work/bulk/*.htseq.results.gz")):
        s = f.split("/")[-1].split("_", 1)[1].replace(".htseq.results.gz", "")
        if s not in o.index: continue
        d = pd.read_csv(f, sep="\t", header=None, names=["id", "c"]); d = d[d.id.str.startswith("ENSG")]; d["sym"] = d.id.map(e2s); d = d.dropna().groupby("sym").c.sum()
        v = np.zeros(len(genes)); ix = cols.reindex(d.index).dropna().astype(int); v[ix.to_numpy()] = d.loc[ix.index].to_numpy(); rows.append(v); names.append(s)
    B = cpm_log(np.array(rows)); g = o.loc[names].group
    return B, np.where(g.eq("adult_sepsis"), "case", np.where(g.eq("adult_HC"), "ctrl", np.where(g.eq("adult_ICU_ctrl"), "icu", "other"))), np.array(names)
def bulk_covid(genes):
    """Imperial ESKD whole-blood/PBMC bulk RNA-seq (Gisby 2022, Zenodo 6497251, htseq counts) of the 164 bulk-only samples the floor model was given"""
    Z = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/bulk_public/zenodo6497251_imperial_eskd/htseq_counts.csv"
    o = pd.read_csv(f"{P}/floor_covid_obs.csv")["sample"].to_numpy(); C = pd.read_csv(Z, index_col=0, usecols=lambda c: c == "gencode_id" or c in set(o))
    sm = np.load("/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/scsemi_eval_20261022/stage_covid/scsemi.npz")
    e2s = dict(zip(sm["gene_ids"].astype(str), sm["gene_symbols"].astype(str))); C.index = C.index.str.split(".").str[0]
    C = C[C.index.isin(e2s)]; C = C.groupby(C.index.map(e2s)).sum().reindex(genes).fillna(0)
    names = np.array(C.columns); g, pat = groups("covid", "floor", names); return cpm_log(C.to_numpy().T.astype(float)), g, pat
def bh(p): p = np.nan_to_num(p, nan=1.0); o = np.argsort(p); q = p[o] * len(p) / (np.arange(len(p)) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
def welch(L, g, ctrl):
    a, b = L[g == "case"], L[g == ctrl]; t, p = stats.ttest_ind(a, b, equal_var=False); return a.mean(0) - b.mean(0), t, p, (len(a), len(b))
def boot_idx(g, pat, ctrl, rng):
    """cluster bootstrap by patient within case and within control; returns index arrays"""
    out = {}
    for lab in ("case", ctrl):
        ix = np.where(g == lab)[0]; up = np.unique(pat[ix]); pick = rng.choice(up, len(up), replace=True)
        out[lab] = np.concatenate([ix[pat[ix] == u] for u in pick])
    return out
def partial_r(y, x, z):
    rxy, rxz, ryz = np.corrcoef(x, y)[0, 1], np.corrcoef(x, z)[0, 1], np.corrcoef(y, z)[0, 1]
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
rows, prow, brow, keep = [], [], [], {}
for coh in ("sepsis", "covid"):
    m0, S0, genes = load(f"real_{coh}"); rank = pd.Series(S0.sum(0)).rank(ascending=False).to_numpy()
    SRC = ["floor", "prev", "scsemi", "scsemi_all"]
    if coh == "sepsis": Bk, gB, nB = bulk_sepsis(genes)
    else: Bk, gB, nB = bulk_covid(genes)
    SRC.append("broadcast")
    for contrast, ctrl in (("case vs control", "ctrl"),) + ((("sepsis vs ICU control", "icu"),) if coh == "sepsis" else ()):
        for t in SIX:
            Lr, gr, pr, _ = per_sample(coh, "real", t); xr, tr_, pv, nR = welch(Lr, gr, ctrl); q = bh(pv); mean_r = Lr.mean(0)
            gen = {}
            for src in SRC:
                if src == "broadcast": Lg, gg, pg = Bk, gB, nB
                else: Lg, gg, pg, _ = per_sample(coh, src, t)
                gen[src] = (Lg, gg, pg); xg, tg, _, nG = welch(Lg, gg, ctrl)
                keep[(coh, contrast, t, src)] = pd.DataFrame({"gene": genes, "real_logFC": xr, "real_q": q, "real_mean": mean_r, "rank": rank, "gen_logFC": xg, "gen_t": tg})
                for band, bm in (("all", np.ones(len(genes), bool)), ("rank <= 2,000", rank <= 2000), ("rank > 2,000", rank > 2000)):
                    ok = bm & (mean_r > 1) & np.isfinite(xr) & np.isfinite(xg); sig = ok & (q < 0.05)
                    top_r = pd.Series(np.abs(tr_), index=genes)[sig].nlargest(100).index; top_g = pd.Series(np.abs(tg), index=genes)[ok].nlargest(100).index
                    rows.append({"cohort": coh, "contrast": contrast, "six_type": t, "source": src, "genes": band, "n_real_case_ctrl": f"{nR[0]}/{nR[1]}", "n_gen_case_ctrl": f"{nG[0]}/{nG[1]}",
                                 "n_genes": int(ok.sum()), "n_real_DEG": int(sig.sum()), "logFC_r": float(np.corrcoef(xr[ok], xg[ok])[0, 1]) if ok.sum() > 10 else np.nan,
                                 "slope_on_DEG": float(np.polyfit(xr[sig], xg[sig], 1)[0]) if sig.sum() > 10 else np.nan,
                                 "sign_agree_DEG": float((np.sign(xr[sig]) == np.sign(xg[sig])).mean()) if sig.sum() else np.nan,
                                 "top100_overlap": float(len(set(top_r) & set(top_g)) / max(len(top_r), 1)) if len(top_r) else np.nan})
            if True:
                xb = keep[(coh, contrast, t, "broadcast")].gen_logFC.to_numpy()
                for band, bm in (("all", np.ones(len(genes), bool)), ("rank > 2,000", rank > 2000)):
                    for src in ("floor", "prev", "scsemi", "scsemi_all"):
                        xg = keep[(coh, contrast, t, src)].gen_logFC.to_numpy(); ok = bm & (mean_r > 1) & np.isfinite(xr) & np.isfinite(xg) & np.isfinite(xb)
                        X = np.column_stack([np.ones(ok.sum()), xb[ok]]); r2b = 1 - np.sum((xr[ok] - X @ np.linalg.lstsq(X, xr[ok], rcond=None)[0]) ** 2) / np.sum((xr[ok] - xr[ok].mean()) ** 2)
                        X2 = np.column_stack([X, xg[ok]]); r2f = 1 - np.sum((xr[ok] - X2 @ np.linalg.lstsq(X2, xr[ok], rcond=None)[0]) ** 2) / np.sum((xr[ok] - xr[ok].mean()) ** 2)
                        prow.append({"cohort": coh, "contrast": contrast, "six_type": t, "genes": band, "source": src, "r_broadcast_real": float(np.corrcoef(xb[ok], xr[ok])[0, 1]),
                                     "r_model_real": float(np.corrcoef(xg[ok], xr[ok])[0, 1]), "r_model_broadcast": float(np.corrcoef(xg[ok], xb[ok])[0, 1]),
                                     "partial_r_model_given_broadcast": float(partial_r(xr[ok], xg[ok], xb[ok])), "R2_broadcast": float(r2b), "R2_broadcast_plus_model": float(r2f)})
            # bootstrap (paired over sources), primary contrast only
            if contrast != "case vs control": continue
            rng = np.random.default_rng(zlib.crc32(f"{coh}|{t}".encode()))
            masks = {b: (bm & (mean_r > 1) & np.isfinite(xr)) for b, bm in (("all", np.ones(len(genes), bool)), ("rank > 2,000", rank > 2000))}
            sig = {b: m & (q < 0.05) for b, m in masks.items()}
            for it in range(NB):
                ir = boot_idx(gr, pr, ctrl, rng); xrb = Lr[ir["case"]].mean(0) - Lr[ir[ctrl]].mean(0)
                for src, (Lg, gg, pg) in gen.items():
                    ig = boot_idx(gg, pg, ctrl, rng); xgb = Lg[ig["case"]].mean(0) - Lg[ig[ctrl]].mean(0)
                    for b in masks:
                        ok = masks[b] & np.isfinite(xgb); s_ = sig[b] & np.isfinite(xgb)
                        brow.append((coh, t, src, b, it, np.corrcoef(xrb[ok], xgb[ok])[0, 1],
                                     np.polyfit(xrb[s_], xgb[s_], 1)[0] if s_.sum() > 10 else np.nan, (np.sign(xrb[s_]) == np.sign(xgb[s_])).mean() if s_.sum() else np.nan))
            print(coh, t, "boot done", flush=True)
M6 = pd.DataFrame(rows); M6.to_csv(f"{OUT}/m6_three_models.csv", index=False); PR = pd.DataFrame(prow); PR.to_csv(f"{OUT}/m6_broadcast_partial.csv", index=False)
BT = pd.DataFrame(brow, columns=["cohort", "six_type", "source", "genes", "it", "logFC_r", "slope_on_DEG", "sign_agree_DEG"]); BT.to_csv(f"{OUT}/m6_bootstrap_draws.csv.gz", index=False)
pd.concat({"|".join(k): v for k, v in keep.items()}).rename_axis(["key", "i"]).reset_index().to_csv(f"{OUT}/m6_three_models_gene_tables.csv.gz", index=False)
# CI of the median over six-types, and paired differences
med = BT.groupby(["cohort", "genes", "source", "it"])[["logFC_r", "slope_on_DEG", "sign_agree_DEG"]].median().reset_index()
ci = []
for (coh, band), d in med.groupby(["cohort", "genes"]):
    W = {m: d.pivot(index="it", columns="source", values=m) for m in ("logFC_r", "slope_on_DEG", "sign_agree_DEG")}
    for m, w in W.items():
        for s in w.columns:
            ci.append({"cohort": coh, "genes": band, "metric": m, "comparison": s, "lo": w[s].quantile(.025), "med": w[s].median(), "hi": w[s].quantile(.975)})
        for a, b in (("floor", "prev"), ("floor", "scsemi"), ("floor", "scsemi_all"), ("prev", "scsemi"), ("floor", "broadcast"), ("prev", "broadcast")):
            dd = w[a] - w[b]; ci.append({"cohort": coh, "genes": band, "metric": m, "comparison": f"{a} - {b}", "lo": dd.quantile(.025), "med": dd.median(), "hi": dd.quantile(.975),
                                         "p_le0": float((dd <= 0).mean())})
CI = pd.DataFrame(ci); CI.to_csv(f"{OUT}/m6_bootstrap_ci.csv", index=False)
pd.set_option("display.width", 250)
print(M6.groupby(["cohort", "contrast", "genes", "source"])[["logFC_r", "slope_on_DEG", "sign_agree_DEG", "top100_overlap"]].median().round(2).to_string())
print(PR.groupby(["genes", "source"])[["r_broadcast_real", "r_model_real", "partial_r_model_given_broadcast", "R2_broadcast", "R2_broadcast_plus_model"]].median().round(3).to_string())
print(CI.round(3).to_string())
