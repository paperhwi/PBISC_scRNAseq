"""F. Per-sample compositional statistics: real single-cell samples vs PBISC-generated samples vs deconvolution of the same bulk.

Unit = sample. Proportions per sample -> CLR (pseudocount 0.001, renormalised). Disease effect per cell type =
mean CLR(case) - mean CLR(control); 95 % CI by resampling samples within group (1,000 bootstrap replicates).
Generated: draws pooled per sample (500 cells per sample) and also each draw on its own (100 cells per sample).
Attenuation: slope of the generated (or deconvolved) effect on the real effect across cell types (types with mean real
share >= 0.5 %), with a bootstrap CI that resamples samples in both sources. Also: Aitchison distance between group
compositions, and the spread across samples (SD of CLR within group) relative to real."""
import json, numpy as np, pandas as pd
O = "/home/claude/work/prop"; R = "/home/claude/work/report"; rng = np.random.default_rng(5); NB = 1000
PC = 1e-3


def lin6(x):
    x = str(x)
    if x == "Plasma cells" or x.startswith("B_ASC"): return "Plasma"
    if x == "B cells" or x.startswith("B_"): return "B"
    if x in ("CD4+ T cells", "CD8+ T cells", "NKT cells", "MAIT", "NKT", "T g/d") or x.startswith("T "): return "T"
    if x in ("NK cells", "NK") or x.startswith("NK "): return "NK"
    if x in ("Monocytes", "Macrophages", "DC", "Int.mono", "ASDC", "pDC") or x.startswith(("CD14mono", "CD16mono", "DC")): return "Myeloid"
    return "Other"


def clr(P):
    Q = P + PC; Q = Q / Q.sum(1, keepdims=True); L = np.log(Q); return L - L.mean(1, keepdims=True)


def props_from_cells(df, unit_cols, col, cats):
    t = df.groupby(unit_cols + [col]).size().unstack(fill_value=0).reindex(columns=cats, fill_value=0)
    return t.div(t.sum(1), axis=0)


def effect(Pc, Pk, idx_c=None, idx_k=None):
    Lc, Lk = clr(Pc), clr(Pk)
    if idx_c is not None: Lc, Lk = Lc[idx_c], Lk[idx_k]
    return Lc.mean(0) - Lk.mean(0)


def boot_effect(Pc, Pk):
    out = np.array([effect(Pc, Pk, rng.integers(0, len(Pc), len(Pc)), rng.integers(0, len(Pk), len(Pk))) for _ in range(NB)])
    return np.percentile(out, 2.5, 0), np.percentile(out, 97.5, 0), out


def aitchison(p, q):
    a, b = clr(p[None])[0], clr(q[None])[0]; return float(np.sqrt(((a - b) ** 2).sum()))


summary, rows = {}, []
for coh in ["sepsis", "covid"]:
    if coh == "sepsis":
        tr = pd.read_csv(f"{O}/sepsis_sc_truth_props.csv.gz"); tr = tr[tr.grp.astype(str).str.startswith("adult")]
        real = {lv: tr[tr.level == lv].pivot_table(index=["sample", "grp"], columns="label", values="frac", fill_value=0) for lv in ["lin6", "major", "fine"]}
        GR = ["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]; CON = [("adult_sepsis", "adult_HC"), ("adult_sepsis", "adult_ICU_ctrl")]
        dec = pd.read_csv(f"{O}/sepsis_deconv_props.csv.gz"); dec = dec[dec.set == "bulkonly"]
    else:
        rl = pd.read_csv(f"{R}/covid_real_labels.csv.gz").rename(columns={"label": "fine"})
        rl = rl[rl.grp.isin(["covid_negative", "covid_positive"])]; rl["lin6"] = rl.major.map(lin6)
        real = {lv: props_from_cells(rl, ["sample", "grp"], lv, sorted(rl[lv].unique())) for lv in ["lin6", "major", "fine"]}
        GR = ["covid_negative", "covid_positive"]; CON = [("covid_positive", "covid_negative")]; dec = None
    gen = {}
    for K in ["k1000", "k4000"]:
        t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False).rename(columns={"label": "fine", "group": "grp", "sample_id": "sample"})
        t = t[t.grp.isin(GR)]; t["lin6"] = t.major.map(lin6); gen[K] = t
    for lv in ["lin6", "major", "fine"]:
        cats = list(real[lv].columns)
        src = {"real": real[lv]}
        for K, t in gen.items():
            src[K] = props_from_cells(t, ["sample", "grp"], lv, cats)
            for d, td in t.groupby("draw"): src[f"{K}_{d}"] = props_from_cells(td, ["sample", "grp"], lv, cats)
        mc = pd.read_csv(f"{O}/{coh}_mapping_counts.csv.gz")
        for K in ["k1000", "k4000"]:
            for rt in ["R3"]:
                d = mc[(mc.source == K) & (mc.route == rt) & (mc.level == lv) & mc.group.isin(GR)]
                p = d.groupby(["sample", "group", "label"]).n.sum().unstack(fill_value=0).reindex(columns=cats, fill_value=0)
                p.index = p.index.set_names(["sample", "grp"]); src[f"{K}_{rt}"] = p.div(p.sum(1), axis=0)
        if dec is not None:
            for m, dm in dec[(dec.level == lv)].groupby("method"):
                p = dm.groupby(["sample", "grp", "label"]).frac.mean().unstack(fill_value=0).reindex(columns=cats, fill_value=0)
                src[f"deconv_{m}"] = p.div(p.sum(1), axis=0)
            # null: pooled adult composition for every bulk-only sample
            pooled = real[lv].mean(0).to_numpy()
            ref_idx = src["deconv_WNNLS"].index
            src["null"] = pd.DataFrame(np.tile(pooled, (len(ref_idx), 1)), index=ref_idx, columns=cats)
        R_ = src["real"]; main = R_.groupby(level="grp").mean().reindex(GR)
        abund = [c for c in cats if main[c].max() >= 0.005]
        for case, ctrl in CON:
            Pc, Pk = R_.xs(case, level="grp").to_numpy(), R_.xs(ctrl, level="grp").to_numpy()
            e_real = effect(Pc, Pk); lo, hi, B_real = boot_effect(Pc, Pk)
            for name, S in src.items():
                if name == "real" or case not in S.index.get_level_values("grp") or ctrl not in S.index.get_level_values("grp"): continue
                Sc, Sk = S.xs(case, level="grp").to_numpy(), S.xs(ctrl, level="grp").to_numpy()
                e = effect(Sc, Sk); l2, h2, B_s = boot_effect(Sc, Sk)
                ai = [cats.index(c) for c in abund]
                slope = np.polyfit(e_real[ai], e[ai], 1)[0] if name != "null" else 0.0
                bs = [np.polyfit(B_real[b, ai], B_s[b, ai], 1)[0] for b in range(200)] if name != "null" else [0.0]
                within = np.mean((e[ai] >= lo[ai]) & (e[ai] <= hi[ai]))
                sign = np.mean(np.sign(e[ai]) == np.sign(e_real[ai]))
                gm = S.groupby(level="grp").mean().reindex(GR)
                key = f"{coh}|{lv}|{case}_vs_{ctrl}|{name}"
                summary[key] = {"slope": float(slope), "slope_ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                                "spearman": float(pd.Series(e[ai]).corr(pd.Series(e_real[ai]), method="spearman")) if name != "null" else None,
                                "sign_agree": float(sign) if name != "null" else None, "share_in_real_CI": float(within), "n_types": len(ai),
                                "aitchison_case": aitchison(gm.loc[case].to_numpy(), main.loc[case].to_numpy()),
                                "aitchison_ctrl": aitchison(gm.loc[ctrl].to_numpy(), main.loc[ctrl].to_numpy()),
                                "L1_case": float(np.abs(gm.loc[case] - main.loc[case]).sum()), "L1_ctrl": float(np.abs(gm.loc[ctrl] - main.loc[ctrl]).sum()),
                                "spread_ratio_case": float(np.median(clr(Sc)[:, ai].std(0) / np.maximum(clr(Pc)[:, ai].std(0), 1e-6))),
                                "n_case": int(len(Sc)), "n_ctrl": int(len(Sk))}
                for j, c in enumerate(cats):
                    rows.append({"cohort": coh, "level": lv, "contrast": f"{case}_vs_{ctrl}", "source": name, "label": c,
                                 "real_share_case": float(main.loc[case, c]), "real_share_ctrl": float(main.loc[ctrl, c]),
                                 "share_case": float(gm.loc[case, c]), "share_ctrl": float(gm.loc[ctrl, c]),
                                 "effect_real": float(e_real[j]), "real_lo": float(lo[j]), "real_hi": float(hi[j]),
                                 "effect": float(e[j]), "lo": float(l2[j]), "hi": float(h2[j]), "abundant": c in abund})
        print(coh, lv, "done", flush=True)
pd.DataFrame(rows).to_csv(f"{O}/compstats_effects.csv.gz", index=False)
json.dump(summary, open(f"{O}/compstats_summary.json", "w"), indent=1)
print("done")
