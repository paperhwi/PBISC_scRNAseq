"""(3) Group-level composition and disease shifts, three models (2026-10-22). Usage: python comp_s3.py
Level six-type: per-sample proportions from all cells of each file (real: author labels mapped to six types; generated: frozen scoring label).
Level fine: real = author fine labels over all real cells (data/pb/real_<coh>_comp_fine.csv); generated and real B (fair comparator, same route)
= kNN-transferred fine labels from real half A (k = 30, Harmony space of each pairwise integration; bundles of ~60k cells, equal per sample).
Fine states with real mean share >= 0.5 % of all cells are scored.
Per source: group-mean composition error (L1 / 2 = total variation, case and control separately), disease shift per type
(mean share case - mean share control) agreement with real: Pearson r, sign agreement, slope; log2 fold of mean share (pseudocount 1e-3).
Bootstrap B = 500: patients resampled within groups, real and generated independently, same real draw for every source (paired differences)."""
import numpy as np, pandas as pd, anndata as ad
P = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb"; P2 = "/home/claude/work/s3/data/pb"
TR = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"; OUT = "/mnt/user-data/outputs/s3/results"; NB = 500
SIXM = {"sepsis": {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell", "NK cells": "NK cell",
                   "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC", "Other cells": "Other"}}
EXCL = {"C126_positive_5", "C126_positive_9", "C126_positive_11", "C141_negative", "C58_positive_2"}
def grp_real(coh):
    if coh == "sepsis":
        o = pd.read_csv(f"{P}/real_sepsis_obs.csv"); g = np.where(o.age_group.eq("adult"), o.group.map({"sepsis": "case", "healthy control": "ctrl", "ICU control": "icu"}), "other")
        return pd.DataFrame({"sample": o["sample"], "g": g, "patient": o["sample"]}).set_index("sample")
    o = pd.read_csv(f"{P}/real_covid_obs.csv"); g = o.case_control.map({"POSITIVE": "case", "NEGATIVE": "ctrl"}).fillna("other")
    return pd.DataFrame({"sample": o["sample"], "g": g, "patient": o.patient.astype(str)}).set_index("sample")
def grp_gen(coh, src):
    if src == "scsemi":
        o = pd.read_csv(f"{P2}/scsemi_{coh}_obs.csv"); g = o.group.map({"adult_sepsis": "case", "adult_HC": "ctrl", "adult_ICU_ctrl": "icu", "covid_positive": "case", "covid_negative": "ctrl"}).fillna("other")
        g = np.where(o.in_A, "other", g); return pd.DataFrame({"sample": o["sample"], "g": g, "patient": o.patient.astype(str)}).set_index("sample")
    o = pd.read_csv(f"{P}/{src}_{coh}_obs.csv"); g = o.group.map({"adult_sepsis": "case", "adult_HC": "ctrl", "adult_ICU_ctrl": "icu", "covid_positive": "case", "covid_negative": "ctrl"}).fillna("other")
    if coh == "covid": g = np.where(o["sample"].isin(EXCL), "other", g)
    pat = o["sample"] if coh == "sepsis" else o["sample"].str.split("_").str[0]
    return pd.DataFrame({"sample": o["sample"], "g": g, "patient": pat.to_numpy()}).set_index("sample")
def share_table(long, key):
    W = long.pivot_table(index="sample", columns=key, values="n", aggfunc="sum", fill_value=0); return W.div(W.sum(1), axis=0)
def tables(coh):
    T = {}
    # six-type
    r = pd.read_csv(f"{P}/real_{coh}_comp_six.csv"); T[("six", "real")] = share_table(r, "six")
    for s in ("floor", "prev"): T[("six", s)] = share_table(pd.read_csv(f"{P}/{s}_{coh}_comp_six.csv"), "six")
    T[("six", "scsemi")] = share_table(pd.read_csv(f"{P2}/scsemi_{coh}_comp_six.csv"), "six")
    # fine
    T[("fine", "real")] = share_table(pd.read_csv(f"{P}/real_{coh}_comp_fine.csv"), "fine")
    for s, f in (("realB_tr", "realU"), ("floor", "floor"), ("prev", "prev4k"), ("scsemi", "scsemi")):
        q = pd.read_csv(f"{TR}/{coh}_{f}_transfer.csv.gz", index_col=0, keep_default_na=False, usecols=["Unnamed: 0", "sample", "fine_tr"] if False else None)
        T[("fine", s)] = share_table(q.groupby(["sample", "fine_tr"]).size().rename("n").reset_index().rename(columns={"fine_tr": "fine"}), "fine")
    return T
def metrics(Wr, gr, Wg, gg, types, ctrl):
    out = {}
    for lab in ("case", ctrl):
        mr = Wr.loc[gr == lab].reindex(columns=types, fill_value=0).mean(); mg = Wg.loc[gg == lab].reindex(columns=types, fill_value=0).mean()
        out[f"TV_{lab if lab == 'case' else 'ctrl'}"] = float(0.5 * np.abs(mr - mg).sum())
    dr = Wr.loc[gr == "case"].reindex(columns=types, fill_value=0).mean() - Wr.loc[gr == ctrl].reindex(columns=types, fill_value=0).mean()
    dg = Wg.loc[gg == "case"].reindex(columns=types, fill_value=0).mean() - Wg.loc[gg == ctrl].reindex(columns=types, fill_value=0).mean()
    lr_ = np.log2((Wr.loc[gr == "case"].reindex(columns=types, fill_value=0).mean() + 1e-3) / (Wr.loc[gr == ctrl].reindex(columns=types, fill_value=0).mean() + 1e-3))
    lg_ = np.log2((Wg.loc[gg == "case"].reindex(columns=types, fill_value=0).mean() + 1e-3) / (Wg.loc[gg == ctrl].reindex(columns=types, fill_value=0).mean() + 1e-3))
    out.update({"shift_r": float(np.corrcoef(dr, dg)[0, 1]), "shift_sign_agree": float((np.sign(dr) == np.sign(dg)).mean()), "shift_slope": float(np.polyfit(dr, dg, 1)[0]),
                "log2ratio_r": float(np.corrcoef(lr_, lg_)[0, 1]), "log2ratio_sign_agree": float((np.sign(lr_) == np.sign(lg_)).mean())})
    return out, dr, dg
def bsamp(G, ctrl, rng):
    keep = []
    for lab in ("case", ctrl):
        d = G[G.g == lab]; up = d.patient.unique(); pick = rng.choice(up, len(up), replace=True)
        keep += [d.index[d.patient == u] for u in pick]
    idx = np.concatenate(keep); return idx
rows, trow, brow = [], [], []
for coh in ("sepsis", "covid"):
    T = tables(coh); Gr = grp_real(coh)
    for contrast, ctrl in (("case vs control", "ctrl"),) + ((("sepsis vs ICU control", "icu"),) if coh == "sepsis" else ()):
        for level in ("six", "fine"):
            Wr = T[(level, "real")]; Wr = Wr.loc[Wr.index.intersection(Gr.index)]; gr = Gr.loc[Wr.index].g.to_numpy()
            types = [c for c in Wr.columns if c not in ("Other", "", "nan")] if level == "six" else list(Wr.columns[(Wr.loc[np.isin(gr, ["case", ctrl])].mean() >= 0.005).to_numpy()])
            srcs = ["floor", "prev", "scsemi"] + (["realB_tr"] if level == "fine" else [])
            for s in srcs:
                Wg = T[(level, s)]
                if s == "realB_tr": Gg = Gr.loc[Wg.index.intersection(Gr.index)]
                else: Gg = grp_gen(coh, s); Gg = Gg.loc[Gg.index.intersection(Wg.index)]
                Wg = Wg.loc[Gg.index]; gg = Gg.g.to_numpy()
                m, dr, dg = metrics(Wr, gr, Wg, gg, types, ctrl)
                rows.append({"cohort": coh, "contrast": contrast, "level": level, "source": s, "n_types": len(types), "n_real_case_ctrl": f"{(gr=='case').sum()}/{(gr==ctrl).sum()}",
                             "n_gen_case_ctrl": f"{(gg=='case').sum()}/{(gg==ctrl).sum()}", **m})
                for ty in types: trow.append({"cohort": coh, "contrast": contrast, "level": level, "source": s, "type": ty, "real_shift": dr[ty], "gen_shift": dg[ty],
                                              "real_case": Wr.loc[gr == "case"].reindex(columns=[ty], fill_value=0)[ty].mean(), "gen_case": Wg.loc[gg == "case"].reindex(columns=[ty], fill_value=0)[ty].mean(),
                                              "real_ctrl": Wr.loc[gr == ctrl].reindex(columns=[ty], fill_value=0)[ty].mean(), "gen_ctrl": Wg.loc[gg == ctrl].reindex(columns=[ty], fill_value=0)[ty].mean()})
            if contrast != "case vs control": continue
            rng = np.random.default_rng(7 + (coh == "covid") + 2 * (level == "fine")); GrS = Gr.loc[Wr.index]
            gens = {}
            for s in srcs:
                Wg = T[(level, s)]
                Gg = Gr.loc[Wg.index.intersection(Gr.index)] if s == "realB_tr" else grp_gen(coh, s)
                Gg = Gg.loc[Gg.index.intersection(Wg.index)]; gens[s] = (Wg.loc[Gg.index], Gg)
            for it in range(NB):
                ir = bsamp(GrS, ctrl, rng); Wrb = Wr.loc[ir]; grb = GrS.loc[ir].g.to_numpy()
                for s, (Wg, Gg) in gens.items():
                    ig = bsamp(Gg, ctrl, rng); m, _, _ = metrics(Wrb, grb, Wg.loc[ig], Gg.loc[ig].g.to_numpy(), types, ctrl)
                    brow.append({"cohort": coh, "level": level, "source": s, "it": it, **m})
        print(coh, contrast, "done", flush=True)
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/composition_shift.csv", index=False); TT = pd.DataFrame(trow); TT.to_csv(f"{OUT}/composition_shift_types.csv", index=False)
B = pd.DataFrame(brow); B.to_csv(f"{OUT}/composition_bootstrap_draws.csv.gz", index=False)
ci = []
for (coh, level), d in B.groupby(["cohort", "level"]):
    for m in ("TV_case", "TV_ctrl", "shift_r", "shift_sign_agree", "log2ratio_r"):
        w = d.pivot(index="it", columns="source", values=m)
        for s in w.columns: ci.append({"cohort": coh, "level": level, "metric": m, "comparison": s, "lo": w[s].quantile(.025), "med": w[s].median(), "hi": w[s].quantile(.975)})
        for a, b in (("floor", "prev"), ("floor", "scsemi"), ("prev", "scsemi")):
            dd = w[a] - w[b]; ci.append({"cohort": coh, "level": level, "metric": m, "comparison": f"{a} - {b}", "lo": dd.quantile(.025), "med": dd.median(), "hi": dd.quantile(.975),
                                         "p_le0": float((dd <= 0).mean())})
CI = pd.DataFrame(ci); CI.to_csv(f"{OUT}/composition_bootstrap_ci.csv", index=False)
pd.set_option("display.width", 250); print(R.round(3).to_string()); print(CI.round(3).to_string())
