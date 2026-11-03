"""Figures and summary tables for the proportion qualification (B-F)."""
import json, logging, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 220})
O = "/home/claude/work/prop"; F = f"{O}/figs"; R = "/home/claude/work/report"
import os; os.makedirs(F, exist_ok=True)
COL = {"real": "#1F2933", "k1000": "#2b83ba", "k4000": "#d7191c", "R1": "#2b83ba", "R2": "#9e9ac8", "R3": "#e69f00", "R4": "#5aa469",
       "BayesPrism": "#5aa469", "BayesPrism+cal": "#2e7d32", "WNNLS+cal": "#8c6d31", "null": "#9aa5b1"}
RN = {"R1": "Harmony + kNN (current)", "R2": "kNN, no Harmony", "R3": "format-aware classifier", "R4": "marker gating"}


def clr(p):
    q = p + 1e-3; q = q / q.sum(); l = np.log(q); return l - l.mean()


# ------------------------------------------------------------------ B1: paired per-sample error
D = pd.read_csv(f"{O}/sepsis_deconv_props.csv.gz"); T = pd.read_csv(f"{O}/sepsis_sc_truth_props.csv.gz")
rows = []
for lv in ["lin6", "major", "fine"]:
    tt = T[T.level == lv].pivot_table(index="sample", columns="label", values="frac", fill_value=0)
    P = D[(D.set == "paired") & (D.level == lv)]
    for m, d in P.groupby("method"):
        dd = d.groupby(["sample", "label"]).frac.mean().unstack(fill_value=0)
        cats = sorted(set(tt.columns) | set(dd.columns)); t2 = tt.reindex(columns=cats, fill_value=0); dd = dd.reindex(columns=cats, fill_value=0)
        s = dd.index.intersection(t2.index); dd, t2 = dd.loc[s], t2.loc[s]
        ab = [c for c in cats if t2[c].mean() >= 0.01]
        rows.append({"level": lv, "method": m, "L1": float((dd - t2).abs().sum(1).median()),
                     "aitchison": float(np.median([np.sqrt(((clr(a) - clr(b)) ** 2).sum()) for a, b in zip(dd.values, t2.values)])),
                     "r_across_samples": float(np.nanmedian([np.corrcoef(dd[c], t2[c])[0, 1] for c in ab])), "n": len(s)})
    adult = T[(T.level == lv) & T.grp.astype(str).str.startswith("adult")].pivot_table(index="sample", columns="label", values="frac", fill_value=0)
    pooled = adult.mean(0); s = tt.index.intersection(P["sample"].unique()); t2 = tt.loc[s].reindex(columns=pooled.index, fill_value=0)
    rows.append({"level": lv, "method": "null", "L1": float((t2 - pooled).abs().sum(1).median()),
                 "aitchison": float(np.median([np.sqrt(((clr(pooled.values) - clr(b)) ** 2).sum()) for b in t2.values])), "r_across_samples": np.nan, "n": len(s)})
B1 = pd.DataFrame(rows); B1.to_csv(f"{O}/B1_paired_error.csv", index=False)

# ------------------------------------------------------------------ F + B2 summary table
S = json.load(open(f"{O}/compstats_summary.json"))
rows = []
for k, v in S.items():
    coh, lv, con, src = k.split("|")
    if "_draw" in src: continue
    rows.append({"cohort": coh, "level": lv, "contrast": con, "source": src, **{x: v[x] for x in ["slope", "spearman", "sign_agree", "L1_case", "L1_ctrl", "aitchison_case", "aitchison_ctrl", "spread_ratio_case", "n_types"]},
                 "slope_lo": v["slope_ci"][0], "slope_hi": v["slope_ci"][1]})
FS = pd.DataFrame(rows); FS.to_csv(f"{O}/F_summary.csv", index=False)

# figure B: two panels
fig, axs = plt.subplots(1, 3, figsize=(13, 4.1))
order = ["null", "NNLS", "SVR", "WNNLS", "BayesPrism", "WNNLS+cal", "BayesPrism+cal"]
for ax, lv in zip(axs[:2], ["lin6", "major"]):
    b = B1[B1.level == lv].set_index("method").reindex([o for o in order if o in set(B1.method)])
    ax.barh(range(len(b)), b.L1, color=[COL.get(m, "#7B8794") for m in b.index]); ax.set_yticks(range(len(b)), b.index)
    ax.axvline(b.loc["null", "L1"], color="#9aa5b1", ls=":", lw=1)
    ax.set_xlabel("median per-sample L1 vs own scRNA (0 = exact)"); ax.set_title(f"Paired samples, {'six lineages' if lv == 'lin6' else 'Ye major types'}", loc="left")
    for i, v in enumerate(b.L1): ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=7.5)
ax = axs[2]
fs = FS[(FS.cohort == "sepsis") & (FS.level == "major") & (FS.contrast == "adult_sepsis_vs_adult_HC")].set_index("source")
srcs = [s for s in ["null", "deconv_WNNLS+cal", "deconv_BayesPrism", "k1000", "k4000", "k1000_R3", "k4000_R3"] if s in fs.index]
lab = {"null": "null", "deconv_WNNLS+cal": "WNNLS+cal", "deconv_BayesPrism": "BayesPrism", "k1000": "PBISC K1k", "k4000": "PBISC K4k", "k1000_R3": "PBISC K1k (R3)", "k4000_R3": "PBISC K4k (R3)"}
x = np.arange(len(srcs)); w = 0.38
ax.bar(x - w / 2, fs.loc[srcs, "L1_ctrl"], w, color="#9aa5b1", label="HC group, L1")
ax.bar(x + w / 2, fs.loc[srcs, "L1_case"], w, color="#3E4C59", label="sepsis group, L1")
ax2 = ax.twinx(); ax2.errorbar(x, fs.loc[srcs, "slope"], yerr=[fs.loc[srcs, "slope"] - fs.loc[srcs, "slope_lo"], fs.loc[srcs, "slope_hi"] - fs.loc[srcs, "slope"]],
                               fmt="o", color="#e69f00", capsize=3, label="disease-effect slope"); ax2.axhline(1, color="#e69f00", ls=":", lw=0.8)
ax2.set_ylabel("slope of effect vs real (1 = exact)", color="#b07800"); ax2.set_ylim(-0.2, 1.8); ax2.spines["right"].set_visible(True)
ax.set_xticks(x, [lab[s] for s in srcs], rotation=30, ha="right", fontsize=8); ax.set_ylabel("group composition L1 vs real")
ax.set_title("Bulk-only samples, Ye major types", loc="left")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper center", bbox_to_anchor=(0.5, -0.36), ncol=3, fontsize=7.5)
plt.tight_layout(); fig.savefig(f"{F}/prop_B_baselines.png"); plt.close(fig)

# ------------------------------------------------------------------ C: routes, major types, case group
def persample_mean(d):
    t = d.groupby(["sample", "label"]).n.sum().unstack(fill_value=0); return t.div(t.sum(1), axis=0).mean(0)


rl = pd.read_csv(f"{R}/covid_real_labels.csv.gz")
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [10, 14]})
Ctab = []
for ax, (coh, case) in zip(axs, [("sepsis", "adult_sepsis"), ("covid", "covid_positive")]):
    L = pd.read_csv(f"{O}/{coh}_mapping_counts.csv.gz")
    if coh == "sepsis": tr = T[(T.level == "major") & (T.grp == case)].groupby("label").frac.mean()
    else:
        r = rl[rl.grp == case]; t = r.groupby(["sample", "major"]).size().unstack(fill_value=0); tr = t.div(t.sum(1), axis=0).mean(0)
    tr = tr[tr >= 0.005].sort_values(ascending=False)
    if coh == "covid": tr = tr.drop([c for c in ["Platelet", "RBC"] if c in tr.index])
    for j, rt in enumerate(["R1", "R2", "R3"]):
        g = persample_mean(L[(L.source == "k1000") & (L.route == rt) & (L.level == "major") & (L.group == case)]).reindex(tr.index, fill_value=0)
        y = np.log2((g + 1e-3) / tr)
        ax.scatter(np.arange(len(tr)) + (j - 1) * 0.22, y, color=COL[rt], s=22, label=RN[rt], zorder=3)
        # positive control: realT2 through the same route vs its truth
        for c in tr.index: Ctab.append({"cohort": coh, "route": rt, "label": c, "real": tr[c], "gen": g[c], "log2_ratio": float(np.log2((g[c] + 1e-3) / tr[c]))})
    lt = L[(L.source == "realT2") & (L.level == "major")]
    tru = persample_mean(lt[lt.route == "truth"]); pr = persample_mean(lt[lt.route == "R1"])
    yc = np.log2((pr.reindex(tr.index, fill_value=0) + 1e-3) / (tru.reindex(tr.index, fill_value=0) + 1e-3))
    ax.scatter(np.arange(len(tr)), yc, marker="_", s=160, color="#7B8794", label="real cells, PBISC-like read-out, route 1 (control)", zorder=2)
    ax.axhline(0, color="k", lw=0.6); ax.axhspan(-0.585, 0.585, color="#F5F7FA", zorder=0)
    ax.set_xticks(range(len(tr)), tr.index, rotation=45, ha="right", fontsize=7.5); ax.set_ylim(-5.5, 2.5)
    ax.set_title({"sepsis": "Sepsis, adult sepsis group (K = 1,000)", "covid": "COVID, COVID+ group (K = 1,000)"}[coh], loc="left")
axs[0].set_ylabel("log2 generated / real share\n(shaded = within 1.5×)")
h, l = axs[1].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, bbox_to_anchor=(0.5, 0))
plt.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(f"{F}/prop_C_routes.png"); plt.close(fig)
pd.DataFrame(Ctab).to_csv(f"{O}/C_routes_major.csv", index=False)

# ------------------------------------------------------------------ D: recoverability of fine labels and ratio on recoverable ones
fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
Dtab = []
for ax, (coh, case) in zip(axs, [("sepsis", "adult_sepsis"), ("covid", "covid_positive")]):
    rc = pd.read_csv(f"{O}/{coh}_confusion_recall.csv"); rc = rc[rc.level == "fine"]
    L = pd.read_csv(f"{O}/{coh}_mapping_counts.csv.gz")
    if coh == "sepsis": tr = T[(T.level == "fine") & (T.grp == case)].groupby("label").frac.mean()
    else:
        r = rl[rl.grp == case]; t = r.groupby(["sample", "label"]).size().unstack(fill_value=0); tr = t.div(t.sum(1), axis=0).mean(0)
    for rt, mk in [("R1", "o"), ("R3", "s")]:
        q = rc[rc.route == rt].set_index("label")
        g = persample_mean(L[(L.source == "k1000") & (L.route == rt) & (L.level == "fine") & (L.group == case)])
        for c in q.index:
            if c not in tr.index or tr[c] < 0.003: continue
            ok = (q.loc[c, "recall"] >= 0.5) and (q.loc[c, "precision"] >= 0.5)
            Dtab.append({"cohort": coh, "route": rt, "label": c, "real_share": tr[c], "gen_share": g.get(c, 0.0), "recall": q.loc[c, "recall"],
                         "precision": q.loc[c, "precision"], "recoverable": ok, "log2_ratio": float(np.log2((g.get(c, 0) + 1e-3) / tr[c]))})
    dt = pd.DataFrame([d for d in Dtab if d["cohort"] == coh])
    for rt in ["R1", "R3"]:
        d = dt[dt.route == rt]
        ax.scatter(d.recall, d.log2_ratio, s=np.clip(d.real_share * 600, 8, 120), c=np.where(d.recoverable, COL[rt], "white"),
                   edgecolors=COL[rt], linewidths=1, label=f"{RN[rt]}", alpha=0.9)
    ax.axvline(0.5, color="#9aa5b1", ls=":", lw=1); ax.axhline(0, color="k", lw=0.6); ax.axhspan(-0.585, 0.585, color="#F5F7FA", zorder=0)
    ax.set_xlabel("recall on real cells with PBISC-like read-out"); ax.set_ylim(-6, 3)
    ax.set_title({"sepsis": "Sepsis fine labels (adult sepsis, K = 1,000)", "covid": "COVID fine labels (COVID+, K = 1,000)"}[coh], loc="left")
    for _, rr in dt[(dt.route == "R3") & dt.recoverable & (dt.log2_ratio.abs() > 1)].iterrows():
        ax.annotate(rr.label, (rr.recall, rr.log2_ratio), fontsize=6.5, xytext=(3, 2), textcoords="offset points", color="#8a5a00")
axs[0].set_ylabel("log2 generated / real share\n(filled = recall and precision ≥ 0.5)"); axs[1].legend(loc="lower right", fontsize=8)
plt.tight_layout(); fig.savefig(f"{F}/prop_D_recover.png"); plt.close(fig)
pd.DataFrame(Dtab).to_csv(f"{O}/D_recoverability.csv", index=False)

# ------------------------------------------------------------------ E2 + F: disease effects
E = pd.read_csv(f"{O}/compstats_effects.csv.gz")
fig, axs = plt.subplots(1, 3, figsize=(13, 4.0))
for ax, (coh, con, lv) in zip(axs[:2], [("sepsis", "adult_sepsis_vs_adult_HC", "fine"), ("covid", "covid_positive_vs_covid_negative", "fine")]):
    for K in ["k1000", "k4000"]:
        e = E[(E.cohort == coh) & (E.contrast == con) & (E.level == lv) & (E.source == K) & E.abundant]
        ax.errorbar(e.effect_real, e.effect, xerr=[e.effect_real - e.real_lo, e.real_hi - e.effect_real], fmt="o", ms=3.5, color=COL[K], ecolor=COL[K], elinewidth=0.4, alpha=0.8,
                    label={"k1000": "K = 1,000", "k4000": "K = 4,000"}[K])
    v = FS[(FS.cohort == coh) & (FS.level == lv) & (FS.contrast == con) & (FS.source == "k1000")].iloc[0]
    lim = max(abs(e.effect_real).max(), abs(e.effect).max()) * 1.15; ax.plot([-lim, lim], [-lim, lim], color="#9aa5b1", lw=0.8, ls=":")
    xs = np.linspace(-lim, lim, 10); ax.plot(xs, v.slope * xs, color=COL["k1000"], lw=1)
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("real effect (Δ mean CLR, 95% CI)"); ax.set_ylabel("generated effect")
    ax.set_title(f"{'Sepsis vs HC' if coh == 'sepsis' else 'COVID+ vs COVID−'}, fine labels: slope {v.slope:.2f}", loc="left"); ax.legend(fontsize=7.5, loc="upper left")
ax = axs[2]
for coh, con, c in [("sepsis", "adult_sepsis_vs_adult_HC", "#2b83ba"), ("covid", "covid_positive_vs_covid_negative", "#d7191c")]:
    nh = pd.read_csv(f"{O}/{coh}_k1000_nhoods.csv.gz"); m = (nh.n_real >= 30) & (nh.n_gen >= 30)
    ax.scatter(nh.loc[m, f"lfc_{con}_real"], nh.loc[m, f"lfc_{con}_gen"], s=4, alpha=0.4, color=c, label={"sepsis": "sepsis vs HC", "covid": "COVID+ vs COVID−"}[coh])
ax.plot([-4, 4], [-4, 4], color="#9aa5b1", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_xlabel("real log2 case / control (neighbourhood)"); ax.set_ylabel("generated log2 case / control")
ax.set_title("Label-free neighbourhoods (K = 1,000)", loc="left"); ax.legend(fontsize=7.5, markerscale=3, loc="upper left")
plt.tight_layout(); fig.savefig(f"{F}/prop_EF_effects.png"); plt.close(fig)
print("figs done")
