"""Figure for analysis A: paired per-patient composition accuracy (from paired_patient_validation_20261014) with the
ensemble-size curve and the equal-information deconvolution baselines on the 58 adult sepsis samples."""
import logging, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 220})
PV = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014/results"; E = "/home/claude/work/prop/eq"
S = pd.read_csv(f"{PV}/paired_accuracy_summary.csv"); S = S[(S.scope == "common_test_paired") & (S.level == "l1")]
C = pd.read_csv(f"{E}/ensemble_size_curve.csv"); Q = pd.read_csv(f"{E}/equal_info_summary.csv").set_index("method")
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.0), gridspec_kw={"width_ratios": [1.1, 1, 1.25]})
ax = axs[0]; y = 0; ticks = []
for coh, lab in [("ye2026_sepsis", "Sepsis, 60 patients"), ("imperial_eskd", "COVID, 35 patients")]:
    s = S[S.cohort == coh].set_index(["method", "draw"])
    items = [("Other-patient mean", s.loc[("other_patient_null", -1)], "#9aa5b1", "s")]
    for K, c in [("k1000", "#2b83ba"), ("k4000", "#d7191c")]:
        d = s.loc[f"PBISC_{K}"]; items.append((f"K = {K[1]},{K[2:]} single draws", None, c, "|", d.L1.min(), d.L1.max(), d.L1.mean()))
        items.append((f"K = {K[1]},{K[2:]} 5-draw average", s.loc[(f"PBISC_{K}_eligible_ensemble", -1)], c, "D"))
    for it in items:
        if it[1] is not None:
            r = it[1]; ax.errorbar(r.L1, y, xerr=[[r.L1 - r.L1_lo], [r.L1_hi - r.L1]], fmt=it[3], color=it[2], ms=6, capsize=3)
        else:
            ax.plot([it[4], it[5]], [y, y], color=it[2], lw=2.5, alpha=0.45); ax.plot(it[6], y, "o", mfc="white", color=it[2], ms=6)
        ticks.append((y, it[0])); y -= 1
    ax.text(0.205, y + len(items) + 0.35, lab, fontsize=9, fontweight="bold"); y -= 1.2
ax.set_yticks([t[0] for t in ticks], [t[1] for t in ticks], fontsize=8); ax.set_xlim(0.2, 0.65); ax.set_xlabel("per-patient composition L1 (six types, lower = better)")
ax.set_title("A. Held-out paired patients", loc="left")
ax = axs[1]
for coh, ls, lab in [("ye2026_sepsis", "-", "sepsis"), ("imperial_eskd", "--", "COVID")]:
    for K, c in [("k1000", "#2b83ba"), ("k4000", "#d7191c")]:
        g = C[(C.cohort == coh) & (C.K == K)].groupby("k").L1.agg(["mean", "min", "max"])
        ax.plot(g.index, g["mean"], ls, marker="o", color=c, ms=4, label=f"{lab}, K = {K[1]},{K[2:]}")
    nul = S[(S.cohort == coh) & (S.method == "other_patient_null")].L1.iloc[0]
    ax.axhline(nul, color="#9aa5b1", ls=ls, lw=1); ax.text(3.0 if lab == "sepsis" else 4.4, nul + 0.003, f"{lab}: other-patient mean", fontsize=7, va="bottom", ha="center", color="#7B8794")
ax.set_xticks(range(1, 6)); ax.set_xlabel("number of adapters averaged"); ax.set_ylabel("per-patient L1"); ax.set_xlim(0.8, 5.2)
ax.set_title("B. Gain from averaging adapters", loc="left"); ax.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 0.92), handlelength=2.6)
ax = axs[2]
rows = [("Average of other patients", "deconv_fold_null", "#9aa5b1"), ("WNNLS + one draw's 5 donors (mean of 5)", None, "#c9a96b"),
        ("PBISC one draw, K = 4,000 (mean of 5)", None, "#f4a09c"),
        ("WNNLS + 5-draw average", "WNNLS+cal ensemble", "#8c6d31"), ("WNNLS + all 25 donors", "WNNLS+cal25", "#8c6d31"),
        ("BayesPrism + all 25 donors", "BayesPrism+cal25", "#2e7d32"), ("PBISC 5-draw average, K = 1,000", "PBISC_k1000_eligible_ensemble", "#2b83ba"),
        ("PBISC 5-draw average, K = 4,000", "PBISC_k4000_eligible_ensemble", "#d7191c")]
singles = {"WNNLS + one draw's 5 donors (mean of 5)": [m for m in Q.index if m.startswith("WNNLS+cal_draw")],
           "PBISC one draw, K = 4,000 (mean of 5)": [m for m in Q.index if m.startswith("PBISC_k4000_draw")]}
for i, (lab, key, c) in enumerate(rows):
    if key is None:
        v = Q.loc[singles[lab], "L1"]; ax.barh(i, v.mean(), color=c); ax.plot([v.min(), v.max()], [i, i], color="k", lw=1); ax.text(v.max() + 0.01, i, f"{v.mean():.2f}", va="center", fontsize=7.5)
    elif key in Q.index:
        r = Q.loc[key]; ax.barh(i, r.L1, color=c); ax.errorbar(r.L1, i, xerr=[[r.L1 - r.lo], [r.hi - r.L1]], fmt="none", ecolor="k", capsize=2, lw=0.8)
        ax.text(r.hi + 0.01, i, f"{r.L1:.2f}", va="center", fontsize=7.5)
ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=7.5); ax.invert_yaxis(); ax.set_xlim(0, 0.62)
ax.set_xlabel("per-sample L1 vs recap scRNA (six types)"); ax.set_title("C. 58 adult sepsis samples, same donors", loc="left")
plt.tight_layout(); fig.savefig("/home/claude/work/prop/figs/prop_A_paired.png"); plt.close(fig); print("ok")
