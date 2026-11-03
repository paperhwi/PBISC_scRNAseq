"""Slide-sized figures for the Result 6 deck (deck palette: real = grey/slate, prev K4000 = D7191C, floor = 2B83BA)."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "axes.edgecolor": "#3E4C59", "axes.labelcolor": "#1F2933", "xtick.color": "#3E4C59", "ytick.color": "#3E4C59"})
R = "/mnt/user-data/outputs/val/results"; F = "/mnt/user-data/outputs/val/figs/deck"
C = {"real": "#3E4C59", "realT": "#AEB7C2", "shuffle": "#E4E7EB", "floor": "#2B83BA", "prev": "#D7191C", "chance": "#9AA5B1"}
LAB = {"real": "Real B (other patients)", "realT": "Real B, depth-matched", "shuffle": "Real B, genes shuffled", "prev": "Previous K = 4,000", "floor": "Floor read-out"}
def save(fig, n): fig.savefig(f"{F}/{n}.png", dpi=250, bbox_inches="tight", facecolor="white"); plt.close(fig)
def grid(ax): ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
def legend(fig, keys, ncol=None, y=-0.05):
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=C[k], label=LAB[k]) for k in keys], loc="lower center", ncol=ncol or len(keys), frameon=False, bbox_to_anchor=(0.5, y))
rd = lambda n: pd.read_csv(f"{R}/{n}")
# A donor fingerprint
FP = rd("covid_donor_fingerprint.csv"); fig, axs = plt.subplots(1, 2, figsize=(8.2, 4.1))
for ax, met, yl in ((axs[0], "top1_acc", "top-1 donor identification"), (axs[1], "own_r_median", "correlation with own donor (median)")):
    for i, gs in enumerate(["top-100 donor-stable genes", "top-300 donor-stable genes", "all expressed genes"]):
        for j, s in enumerate(("prev", "floor")):
            v = FP[(FP.gene_set == gs) & (FP.source == s)][met].item(); ax.bar(i + j * 0.38, v, 0.36, color=C[s]); ax.text(i + j * 0.38, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_xticks([0.19, 1.19, 2.19], ["stable\ntop 100", "stable\ntop 300", "all\nexpressed"]); ax.set_ylabel(yl); grid(ax)
axs[0].axhline(1 / 32, color=C["chance"], ls="--", lw=1); axs[0].text(-0.25, 1 / 32 + 0.012, "chance 1/32", color=C["chance"], fontsize=9, ha="left"); axs[0].set_ylim(0, 0.5)
axs[1].axhline(np.median(FP.other_r_median), color=C["chance"], ls="--", lw=1); axs[1].text(-0.25, np.median(FP.other_r_median) - 0.035, "other donors (median)", color=C["chance"], fontsize=9, ha="left"); axs[1].set_ylim(-0.08, 0.75)
legend(fig, ["prev", "floor"]); fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "d_donor")
# B M6
M = rd("m6_disease_effect.csv"); M = M[M.contrast == "case vs control"]
fig, axs = plt.subplots(1, 3, figsize=(12, 4.0))
for ax, met, yl in zip(axs, ["logFC_r", "slope_on_DEG", "sign_agree_DEG"], ["logFC r (real vs generated)", "slope on real DEGs", "sign agreement on real DEGs"]):
    k = 0; xt, xl = [], []
    for coh in ("sepsis", "covid"):
        for band in ("rank <= 2,000", "rank > 2,000"):
            for j, s in enumerate(("prev", "floor")):
                d = M[(M.cohort == coh) & (M.genes == band) & (M.source == s)][met]; ax.bar(k + j * 0.38, d.median(), 0.36, color=C[s]); ax.scatter(np.full(len(d), k + j * 0.38), d, s=9, color="#1F2933", zorder=3)
            xt.append(k + 0.19); xl.append(f"{'Sepsis' if coh == 'sepsis' else 'COVID'}\n{band.replace('rank ', 'rank ')}"); k += 1
    ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel(yl); grid(ax); ax.set_ylim(0, 1.05 if met != "slope_on_DEG" else 0.9)
legend(fig, ["prev", "floor"]); fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "d_m6")
# C c2st
fig, ax = plt.subplots(figsize=(8.4, 4.1)); k = 0; xt, xl = [], []
for coh in ("sepsis", "covid"):
    T = rd(f"{coh}_c2st.csv")
    for lin in ("T", "NK", "B", "Myeloid"):
        for j, (s, kk) in enumerate((("realU", "real"), ("realT1f", "realT"), ("prev4k", "prev"), ("floor", "floor"))):
            ax.bar(k + j * 0.2, T[(T.lineage == lin) & (T.source == s)].c2st_auc.item(), 0.19, color=C[kk])
        xt.append(k + 0.3); xl.append(f"{lin}\n{'Sepsis' if coh == 'sepsis' else 'COVID'}"); k += 1
    k += 0.4
ax.axhline(0.5, color=C["chance"], ls="--", lw=1); ax.set_ylim(0.3, 1.02); ax.set_xticks(xt, xl, fontsize=9); ax.set_ylabel("classifier AUC vs real A\n(0.5 = indistinguishable)"); grid(ax)
legend(fig, ["real", "realT", "prev", "floor"], ncol=4, y=-0.08); fig.tight_layout(rect=(0, 0.06, 1, 1)); save(fig, "d_c2st")
# C coverage / density
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0))
for ax, met in zip(axs, ("coverage", "density")):
    k = 0; xt, xl = [], []
    for coh in ("sepsis", "covid"):
        T = rd(f"{coh}_density_coverage.csv"); T = T[T.lineage == "all"]
        for sp_ in ("PCA (no correction)", "Harmony"):
            for j, (s, kk) in enumerate((("realU", "real"), ("prev4k", "prev"), ("floor", "floor"))):
                v = T[(T.space == sp_) & (T.source == s)][met].item(); ax.bar(k + j * 0.27, v, 0.25, color=C[kk]); ax.text(k + j * 0.27, v + 0.01, f"{v:.2f}", ha="center", fontsize=8)
            xt.append(k + 0.27); xl.append(f"{'Sepsis' if coh == 'sepsis' else 'COVID'}\n{'before Harmony' if sp_.startswith('PCA') else 'after Harmony'}"); k += 1
    ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel(met + " (k = 5)"); grid(ax); ax.set_title(met.capitalize(), loc="left", fontsize=11)
legend(fig, ["real", "prev", "floor"]); fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "d_dc")
# D co-expression
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0))
for ax, met, yl in ((axs[0], "corr_of_corr", "r of gene–gene correlations vs real A"), (axs[1], "mean_abs_corr_src", "mean |r| between genes")):
    k = 0; xt, xl = [], []
    for coh in ("sepsis", "covid"):
        T = rd(f"{coh}_coexpression.csv")
        for band in ("top 300 genes", "ranks 1,001-1,300"):
            keys = (("realU", "real"), ("shuffle", "shuffle"), ("prev4k", "prev"), ("floor", "floor")) if met == "corr_of_corr" else (("ref", None), ("realU", "real"), ("prev4k", "prev"), ("floor", "floor"))
            for j, (s, kk) in enumerate(keys):
                if s == "ref": v = T[(T.genes == band) & (T.source == "realU")].mean_abs_corr_real; col = "#7F7F7F"
                else: v = T[(T.genes == band) & (T.source == s)][met]; col = C[kk]
                ax.bar(k + j * 0.2, v.median(), 0.19, color=col); ax.scatter(np.full(len(v), k + j * 0.2), v, s=7, color="#1F2933", zorder=3)
            xt.append(k + 0.3); xl.append(f"{'Sepsis' if coh == 'sepsis' else 'COVID'}\n{'top 300' if band.startswith('top') else 'ranks 1,001–1,300'}"); k += 1
    ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel(yl); grid(ax)
axs[0].set_title("Structure (dots = 4 lineages)", loc="left", fontsize=11); axs[1].set_title("Strength (leftmost grey = real A)", loc="left", fontsize=11)
legend(fig, ["real", "shuffle", "prev", "floor"]); fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "d_coexpr")
# E LR
fig, axs = plt.subplots(1, 3, figsize=(12, 4.0))
for ax, met, yl in zip(axs[:2], ("detected_jaccard", "recall_top200_realA"), ("detected interactions, Jaccard vs real A", "recall of real A top-200 interactions")):
    for k, coh in enumerate(("sepsis", "covid")):
        T = rd(f"{coh}_lr_recovery.csv").set_index("source")
        for j, (s, kk) in enumerate((("realU", "real"), ("realT1f", "realT"), ("prev4k", "prev"), ("floor", "floor"))):
            ax.bar(k + j * 0.2, T.loc[s, met], 0.19, color=C[kk]); ax.text(k + j * 0.2, T.loc[s, met] + 0.01, f"{T.loc[s, met]:.2f}", ha="center", fontsize=8)
    ax.set_xticks([0.3, 1.3], ["Sepsis", "COVID"]); ax.set_ylabel(yl); grid(ax); ax.set_ylim(0, 1.08)
ax = axs[2]
for k, coh in enumerate(("sepsis", "covid")):
    T = rd(f"{coh}_lr_disease.csv").set_index("source")
    for j, (s, kk) in enumerate((("prev4k", "prev"), ("floor", "floor"))):
        v = T.loc[s, "disease_spearman_all_detected"]; ax.bar(k + j * 0.38, v, 0.36, color=C[kk]); ax.text(k + j * 0.38, v + 0.01, f"{v:.2f}", ha="center", fontsize=8)
ax.set_xticks([0.19, 1.19], ["Sepsis", "COVID"]); ax.set_ylabel("disease change of scores, Spearman"); grid(ax); ax.set_ylim(0, 1.0)
legend(fig, ["real", "realT", "prev", "floor"]); fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "d_lr")
# F claims
CL = pd.concat([rd(f"{c}_claims.csv") for c in ("sepsis", "covid")]); CL = CL[(CL.auc_realA >= 0.6) & (CL.claim != "ESKD-IFN B cell")]
cl = CL.claim.unique(); x = np.arange(len(cl)); fig, ax = plt.subplots(figsize=(8.4, 4.2))
for j, (s, kk) in enumerate((("realU", "real"), ("prev4k", "prev"), ("floor", "floor"))):
    ax.bar(x + (j - 1) * 0.27, [CL[(CL.claim == k_) & (CL.source == s)].delta_pos_frac.item() for k_ in cl], 0.25, color=C[kk])
ax.axhline(0, color="#3E4C59", lw=0.8); ax.axvline(4.5, color="#C9CED6", lw=0.8, ls=":"); ax.set_ylim(-0.32, 0.34)
ax.set_xticks(x, [k_.replace("ESKD-IFN ", "IFN · ") for k_ in cl], rotation=30, ha="right", fontsize=9)
ax.text(2, 0.31, "Sepsis (Ye 2026)", ha="center", fontsize=9.5, color="#3E4C59"); ax.text(6.5, 0.31, "COVID (Stephenson 2025)", ha="center", fontsize=9.5, color="#3E4C59")
ax.set_ylabel("Δ classifier-positive fraction\n(case − control, sample means)"); grid(ax)
legend(fig, ["real", "prev", "floor"], y=-0.08); fig.tight_layout(rect=(0, 0.06, 1, 1)); save(fig, "d_claims")
# G rare transfer-free
RR = pd.concat([rd(f"{c}_rare_transfer_free.csv") for c in ("sepsis", "covid")]); RR = RR[RR.cv_auc_realA >= 0.8]
RR["ratio"] = RR.gen_estimate_pct / RR.realB_estimate_pct; W = RR.pivot_table(index=["cohort", "subset"], columns="model_depth", values="ratio").reset_index().sort_values(["cohort", "floor"])
fig, ax = plt.subplots(figsize=(8.0, 5.0)); y = np.arange(len(W))
for i, r in enumerate(W.itertuples()): ax.plot([np.clip(r.prev4k, 0.01, 20), np.clip(r.floor, 0.01, 20)], [i, i], color="#C9CED6", lw=1.2, zorder=1)
ax.scatter(np.clip(W.prev4k, 0.01, 20), y, color=C["prev"], s=34, label=LAB["prev"], zorder=3); ax.scatter(np.clip(W.floor, 0.01, 20), y, color=C["floor"], s=34, label=LAB["floor"], zorder=3)
ax.axvspan(0.5, 2, color="#EAF3FA", zorder=0); ax.axvline(1, color="#9AA5B1", lw=0.8, ls=":"); ax.set_xscale("log"); ax.set_xlim(0.01, 20)
ax.set_yticks(y, [f"{r.subset}  ({'S' if r.cohort == 'sepsis' else 'C'})" for r in W.itertuples()], fontsize=8.5); ax.set_xlabel("generated / real B share, same classifier (shaded = 0.5–2×)")
ax.legend(frameon=False, fontsize=9, loc="lower right"); fig.tight_layout(); save(fig, "d_rare")
# H integration
I = pd.concat([rd(f"{c}_integration_methods.csv") for c in ("sepsis", "covid")]); meths = ["PCA (none)", "Harmony", "Scanorama", "scVI", "BBKNN (graph)"]
fig, ax = plt.subplots(figsize=(9.0, 4.1)); k = 0; xt, xl = [], []
for coh in ("sepsis", "covid"):
    for m_ in meths:
        for j, (s, kk) in enumerate((("realU", "real"), ("prev4k", "prev"), ("floor", "floor"))):
            ax.bar(k + j * 0.27, I[(I.cohort == coh) & (I.method == m_) & (I.pair == s)].iLISI.item(), 0.25, color=C[kk])
        xt.append(k + 0.27); xl.append(f"{m_.replace(' (graph)', '').replace('PCA (none)', 'none')}\n{'Sepsis' if coh == 'sepsis' else 'COVID'}"); k += 1
    k += 0.4
ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylim(1, 2); ax.set_ylabel("median iLISI of source (2 = mixed)"); grid(ax)
legend(fig, ["real", "prev", "floor"], y=-0.08); fig.tight_layout(rect=(0, 0.06, 1, 1)); save(fig, "d_integ")
print("deck figs ok")
