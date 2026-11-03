"""Step 5: Fig 6 draft panels + post-hoc sensitivity (MIN_PARENT) for the generated arm.

Outputs:
  figures/fig6_panels.png / .pdf
  results/sensitivity_min_parent.csv   (post-hoc, labelled as such in REPORT.md)
"""
import os, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")

ROOT = "F:/PBISC/fig6_bulkonly_reproduction_20261020"
os.makedirs(f"{ROOT}/figures", exist_ok=True)
TG = ["F1", "F2", "F3", "F4", "C1"]
NAME = {"F1": "F1 CD4 GPR183-NR4A2", "F2": "F2 CD8 GNLY-CCL4hi", "F3": "F3 NK PRF1-CCL4hi", "F4": "F4 NKT IFNG-CCL4hi", "C1": "C1 CD4 naive LEF1"}
COL = {"real sc, author labels": "#5B6B85", "real sc, instrument": "#1B2A4A", "PBISC generated": "#0F7C8C",
       "B1 bulk marker score": "#E0A526", "B2 signature NNLS": "#D9643F"}
E = pd.read_csv(f"{ROOT}/results/effects.csv")
S = pd.read_csv(f"{ROOT}/results/site_effects.csv")
P = pd.read_csv(f"{ROOT}/results/classifier_performance.csv")
V = pd.read_csv(f"{ROOT}/results/verdicts.csv")
e = E[E.contrast == "sepsis vs HC"]

fig = plt.figure(figsize=(15, 9.5))
gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1], hspace=0.38, wspace=0.42)

# ---- A: sepsis vs HC effect, logit-fraction sources
ax = fig.add_subplot(gs[:, 0])
srcs = ["real sc, author labels", "real sc, instrument", "PBISC generated", "B2 signature NNLS"]
off = dict(zip(srcs, [0.27, 0.09, -0.09, -0.27]))
for i, k in enumerate(TG):
    y0 = len(TG) - 1 - i
    for s in srcs:
        d = e[(e.source == s) & (e.target == k) & e.est.notna()]
        if s == "PBISC generated":
            ax.scatter(d.est, np.full(len(d), y0 + off[s]), color=COL[s], s=22, alpha=.75, zorder=3)
            if len(d): ax.plot([d.est.median()] * 2, [y0 + off[s] - .07, y0 + off[s] + .07], color=COL[s], lw=3)
        elif len(d):
            r = d.iloc[0]
            ax.errorbar(r.est, y0 + off[s], xerr=1.96 * r.se, fmt="o", color=COL[s], ms=5, capsize=2, zorder=3)
ax.axvline(0, color="grey", lw=.8, ls="--")
ax.set_yticks(range(len(TG))); ax.set_yticklabels([NAME[k] for k in TG[::-1]])
ax.set_xlabel("sepsis vs HC, logit(state / parent) difference")
ax.set_title("A  Direction and size: real sc vs bulk-only generated", loc="left", fontsize=11)
for s in srcs:
    ax.scatter([], [], color=COL[s], label=s + (" (5 draws; bar = median)" if s == "PBISC generated" else " (95% CI)"))
ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2, frameon=False)
for i, k in enumerate(TG):
    v = V[V.target == k].iloc[0]
    ax.text(1.02, (len(TG) - 1 - i + .5) / len(TG), f"{v.verdict}\n{v.draws_positive}/5 pos, {v.draws_significant} sig",
            transform=ax.transAxes, fontsize=7.5, ha="left", va="center", color="#0F7C8C")

# ---- B: B1 marker score (z units)
ax = fig.add_subplot(gs[0, 1])
d = e[e.source == "B1 bulk marker score"].set_index("target").reindex(TG)
ax.barh(range(len(TG))[::-1], d.est, xerr=1.96 * d.se, color=COL["B1 bulk marker score"], capsize=2)
ax.axvline(0, color="grey", lw=.8, ls="--"); ax.set_yticks(range(len(TG))[::-1]); ax.set_yticklabels(TG)
ax.set_xlabel("sepsis vs HC, mean marker z difference")
ax.set_title("B  Bulk marker score (comparator B1)", loc="left", fontsize=11)

# ---- C: site pattern, real instrument vs generated median
ax = fig.add_subplot(gs[0, 2])
sites = ["Abdominal", "Respiratory", "Urinary", "Skin"]
for j, k in enumerate(TG[:4]):
    g = S[(S.source == "PBISC generated") & (S.target == k)].groupby("site").est.median().reindex(sites)
    r = S[(S.source == "real sc, instrument") & (S.target == k)].set_index("site").est.reindex(sites)
    ax.plot(sites, g.values, "-o", color=plt.cm.viridis(j / 4), label=f"{k} generated")
    ax.plot(sites, r.values, "--x", color=plt.cm.viridis(j / 4), alpha=.6)
ax.axhline(0, color="grey", lw=.8)
ax.set_ylabel("site vs HC effect"); ax.tick_params(axis="x", labelsize=8)
ax.set_title("C  Infection-site pattern (solid gen, dashed real)", loc="left", fontsize=11)
ax.legend(fontsize=7, frameon=False)

# ---- D: instrument held-out AUC
ax = fig.add_subplot(gs[1, 1])
p = P[P.stage == "state"].set_index("target").reindex(TG)
ax.bar(TG, p.auc, color="#5B6B85"); ax.set_ylim(.5, 1)
for i, v in enumerate(p.auc): ax.text(i, v + .01, f"{v:.2f}", ha="center", fontsize=8)
ax.set_ylabel("held-out AUC (sample-grouped)")
ax.set_title("D  Same-instrument classifier, thinned real", loc="left", fontsize=11)

# ---- E: precision of the effect (SE) real (n HC=7) vs generated (n HC=36)
ax = fig.add_subplot(gs[1, 2])
r = e[e.source == "real sc, instrument"].set_index("target").se.reindex(TG)
g = e[e.source == "PBISC generated"].groupby("target").se.median().reindex(TG)
x = np.arange(len(TG))
ax.bar(x - .2, r, .4, color=COL["real sc, instrument"], label="real sc (HC n=7)")
ax.bar(x + .2, g, .4, color=COL["PBISC generated"], label="generated (HC n=36), median of draws")
ax.set_xticks(x); ax.set_xticklabels(TG); ax.set_ylabel("SE of sepsis vs HC effect")
ax.set_title("E  Precision gained from bulk-only cohort", loc="left", fontsize=11)
ax.legend(fontsize=7.5, frameon=False)

fig.savefig(f"{ROOT}/figures/fig6_panels.png", dpi=170, bbox_inches="tight")
fig.savefig(f"{ROOT}/figures/fig6_panels.pdf", bbox_inches="tight")

# ---- post-hoc sensitivity: MIN_PARENT threshold, generated arm, sepsis vs HC
meta = pd.read_csv("F:/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/samples.csv").rename(columns={"sample_id": "sample"})
meta = meta[meta.group.isin(["adult_HC", "adult_sepsis"])]
meta["is_c"] = (meta.group == "adult_sepsis").astype(int); meta["sex"] = meta.sex_geo
cg = pd.read_csv(f"{ROOT}/results/calls_gen.csv").merge(meta, on="sample")
rows = []
for mp in [5, 10, 20]:
    for j in range(5):
        for k in TG:
            d = cg[(cg.draw == j) & (cg.target == k) & (cg.parent_n >= mp)].copy()
            n1, n0 = int(d.is_c.sum()), int((1 - d.is_c).sum())
            if n1 < 3 or n0 < 3:
                rows.append(dict(min_parent=mp, draw=j, target=k, est=np.nan, p=np.nan, n_case=n1, n_ref=n0)); continue
            d["y"] = np.log((d.state_n + .5) / (d.parent_n - d.state_n + .5))
            m = smf.ols("y ~ is_c + age + C(sex)", data=d).fit(cov_type="HC3")
            rows.append(dict(min_parent=mp, draw=j, target=k, est=m.params.is_c, p=m.pvalues.is_c, n_case=n1, n_ref=n0))
sens = pd.DataFrame(rows)
sens.to_csv(f"{ROOT}/results/sensitivity_min_parent.csv", index=False)
print(sens[sens.target == "F4"].round(3).to_string(index=False))
print(sens.groupby(["min_parent", "target"]).apply(lambda d: f"{int((d.est > 0).sum())}/5 pos, {int(((d.est > 0) & (d.p < .05)).sum())} p<.05").unstack())
