"""Figures for the three-model comparison (2026-10-22). Usage: python figs_s3.py"""
import json, numpy as np, pandas as pd, anndata as ad, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "axes.edgecolor": "#3E4C59", "axes.labelcolor": "#1F2933", "xtick.color": "#3E4C59", "ytick.color": "#3E4C59"})
R = "/mnt/user-data/outputs/s3/results"; F = "/mnt/user-data/outputs/s3/figs"; TR = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"
C = {"real": "#3E4C59", "realT": "#AEB7C2", "floor": "#2B83BA", "prev": "#D7191C", "scsemi": "#7B3294", "broadcast": "#E08214", "chance": "#9AA5B1"}
NM = {"scsemi": "scSemiProfiler", "prev": "PBISC prev K4000", "floor": "PBISC floor", "broadcast": "bulk broadcast", "real": "real B", "realT": "real B, model depth"}
def save(fig, n): fig.savefig(f"{F}/{n}.png", dpi=220, bbox_inches="tight", facecolor="white"); fig.savefig(f"{F}/{n}.pdf", bbox_inches="tight"); plt.close(fig)
def grid(ax): ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
def legend(fig, keys, y=-0.05): fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=C[k], label=NM[k]) for k in keys], loc="lower center", ncol=len(keys), frameon=False, bbox_to_anchor=(0.5, y))
# 1. sparsity: genes per cell and UMI per cell, all cells of each bundle
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
for j, coh in enumerate(("sepsis", "covid")):
    ax = axs[j]
    for k, (s, c) in enumerate((("realU", "real"), ("prev4k", "prev"), ("floor", "floor"), ("scsemi", "scsemi"))):
        a = ad.read_h5ad(f"/home/claude/work/s3/big/{coh}_{s}.h5ad"); ng = np.diff(a.X.tocsr().indptr)
        p = ax.violinplot([ng], positions=[k], widths=0.8, showmedians=True, showextrema=False)
        for b in p["bodies"]: b.set_facecolor(C[c]); b.set_alpha(0.75)
        p["cmedians"].set_color("#1F2933"); ax.text(k, np.median(ng) * 1.08, f"{np.median(ng):.0f}", ha="center", fontsize=8.5)
    ax.set_yscale("log"); ax.set_xticks(range(4), [NM[x] for x in ("real", "prev", "floor", "scsemi")], fontsize=8.5); ax.set_ylabel("genes detected per cell (log)"); grid(ax)
    ax.set_title(f"{coh.capitalize()}: genes per cell (panel genes)", loc="left", fontsize=10.5)
fig.tight_layout(); save(fig, "s3_genes_per_cell")
# 2. M6 bars with bootstrap CI
M = pd.read_csv(f"{R}/m6_three_models.csv"); M = M[M.contrast == "case vs control"]; CI = pd.read_csv(f"{R}/m6_bootstrap_ci.csv")
fig, axs = plt.subplots(1, 3, figsize=(15, 3.9))
for ax, met, yl in zip(axs, ["logFC_r", "slope_on_DEG", "sign_agree_DEG"], ["logFC r (real vs generated)", "slope on real DEGs", "sign agreement on real DEGs"]):
    k = 0; xt, xl = [], []
    for coh in ("sepsis", "covid"):
        for band in ("all", "rank > 2,000"):
            srcs = ("broadcast", "scsemi", "prev", "floor")
            for j, s in enumerate(srcs):
                d = M[(M.cohort == coh) & (M.genes == band) & (M.source == s)][met]; x = k + j * 0.2
                ax.bar(x, d.median(), 0.19, color=C[s]); ax.scatter(np.full(len(d), x), d, s=6, color="#1F2933", zorder=3)
                ci = CI[(CI.cohort == coh) & (CI.genes == band) & (CI.metric == met) & (CI.comparison == s)]
                if len(ci): ax.errorbar(x, ci.med.item(), yerr=[[ci.med.item() - ci.lo.item()], [ci.hi.item() - ci.med.item()]], color="#1F2933", lw=1, capsize=2, zorder=4, fmt="none")
            xt.append(k + 0.3); xl.append(f"{coh.capitalize()}\n{band}"); k += 1.1
    ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel(yl); grid(ax); ax.axhline(0, color="#3E4C59", lw=0.6)
axs[0].set_title("Disease effect per cell type (bar = median of 5 types, whisker = patient-bootstrap 95 % CI)", loc="left", fontsize=10.5)
legend(fig, ("broadcast", "scsemi", "prev", "floor")); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "s3_m6_bootstrap")
# 3. broadcast partial r (both cohorts)
PR = pd.read_csv(f"{R}/m6_broadcast_partial.csv"); PR = PR[(PR.contrast == "case vs control")]
fig, axs = plt.subplots(2, 2, figsize=(13, 7.2), sharey=True)
types = ["T cell", "NK cell", "B cell", "Monocyte", "DC"]
for i, coh in enumerate(("sepsis", "covid")):
    for j, band in enumerate(("all", "rank > 2,000")):
        ax = axs[i, j]; d = PR[(PR.genes == band) & (PR.cohort == coh)]
        for k, t in enumerate(types):
            b = d[(d.six_type == t)].set_index("source"); ax.bar(k - 0.3, b.r_broadcast_real.iloc[0], 0.19, color=C["broadcast"])
            for jj, s in enumerate(("scsemi", "prev", "floor")): ax.bar(k - 0.1 + jj * 0.2, b.loc[s, "partial_r_model_given_broadcast"], 0.19, color=C[s])
        ax.set_xticks(range(5), types); ax.axhline(0, color="#3E4C59", lw=0.6); grid(ax); ax.set_ylabel("r with real logFC")
        ax.set_title(f"{coh.capitalize()}, {band} genes: broadcast r | model partial r given broadcast", loc="left", fontsize=10)
legend(fig, ("broadcast", "scsemi", "prev", "floor"), y=-0.03); fig.tight_layout(rect=(0, 0.03, 1, 1)); save(fig, "s3_broadcast_partial")
# 4. composition shifts: real vs generated
T = pd.read_csv(f"{R}/composition_shift_types.csv"); T = T[T.contrast == "case vs control"]; CP = pd.read_csv(f"{R}/composition_shift.csv"); CP = CP[CP.contrast == "case vs control"]
fig, axs = plt.subplots(2, 3, figsize=(13, 8))
for i, (coh, lvl) in enumerate((("sepsis", "fine"), ("covid", "fine"))):
    for j, s in enumerate(("scsemi", "prev", "floor")):
        ax = axs[i, j]; d = T[(T.cohort == coh) & (T.level == lvl) & (T.source == s)]; b = T[(T.cohort == coh) & (T.level == lvl) & (T.source == "realB_tr")]
        ax.scatter(100 * b.real_shift, 100 * b.gen_shift, s=14, color=C["realT"], label="real B (transferred)")
        ax.scatter(100 * d.real_shift, 100 * d.gen_shift, s=18, color=C[s])
        lim = 100 * np.nanmax(np.abs(np.r_[d.real_shift, d.gen_shift, b.gen_shift])) * 1.1; ax.plot([-lim, lim], [-lim, lim], color="#9AA5B1", lw=0.8, ls="--"); ax.axhline(0, color="#CBD2D9", lw=0.6); ax.axvline(0, color="#CBD2D9", lw=0.6)
        for _, r in d.reindex(d.real_shift.abs().sort_values(ascending=False).index).head(4).iterrows(): ax.annotate(r["type"][:18], (100 * r.real_shift, 100 * r.gen_shift), fontsize=7, color="#3E4C59")
        m = CP[(CP.cohort == coh) & (CP.level == lvl) & (CP.source == s)].iloc[0]
        ax.set_title(f"{coh.capitalize()} {lvl}: {NM[s]}\nshift r = {m.shift_r:.2f}, sign {m.shift_sign_agree:.2f}", loc="left", fontsize=9.5)
        ax.set_xlabel("real shift (case − control, % of cells)"); ax.set_ylabel("generated shift")
fig.tight_layout(); save(fig, "s3_composition_shift")
# 5. realism: C2ST per lineage, coverage/density
fig, axs = plt.subplots(1, 3, figsize=(16, 3.9), gridspec_kw={"width_ratios": [1.6, 1, 1]})
ax = axs[0]; k = 0; xt, xl = [], []
for coh in ("sepsis", "covid"):
    t = pd.read_csv(f"{R}/{coh}_c2st.csv")
    for lin in ("T", "NK", "B", "Myeloid"):
        for j, (s, c) in enumerate((("realU", "real"), ("realT1s", "realT"), ("scsemi", "scsemi"), ("prev4k", "prev"), ("floor", "floor"))):
            ax.bar(k + j * 0.16, t[(t.lineage == lin) & (t.source == s)].c2st_auc.item(), 0.15, color=C[c])
        xt.append(k + 0.32); xl.append(f"{lin}\n{coh.capitalize()}"); k += 1
    k += 0.4
ax.axhline(0.5, color=C["chance"], ls="--", lw=1); ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylim(0.3, 1.02); ax.set_ylabel("C2ST AUC vs real A"); grid(ax)
ax.set_title("Classifier two-sample test (0.5 = indistinguishable)", loc="left", fontsize=10.5)
DC = pd.concat([pd.read_csv(f"{R}/{c}_density_coverage.csv") for c in ("sepsis", "covid")]); DC = DC[(DC.lineage == "all") & (DC.space == "PCA (no correction)")]
for ax, met in zip(axs[1:], ("coverage", "density")):
    for i, coh in enumerate(("sepsis", "covid")):
        for j, (s, c) in enumerate((("realU", "real"), ("scsemi", "scsemi"), ("prev4k", "prev"), ("floor", "floor"))):
            v = DC[(DC.cohort == coh) & (DC.source == s)][met].item(); ax.bar(i + j * 0.2, v, 0.19, color=C[c]); ax.text(i + j * 0.2, v + 0.01, f"{v:.2f}", ha="center", fontsize=7.5)
    ax.set_xticks([0.3, 1.3], ["Sepsis", "COVID"]); ax.set_ylabel(met); grid(ax); ax.set_title(f"{met.capitalize()} (k = 5, PCA, no correction)", loc="left", fontsize=10.5)
legend(fig, ("real", "realT", "scsemi", "prev", "floor")); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "s3_realism")
# 6. myeloid diagnosis
MY = pd.concat([pd.read_csv(f"{R}/{c}_myeloid_identity.csv") for c in ("sepsis", "covid")]); G = pd.read_csv(f"{R}/sepsis_myeloid_model_specific_genes.csv")
fig, axs = plt.subplots(1, 3, figsize=(16, 4.1), gridspec_kw={"width_ratios": [1, 1, 1.3]})
srcs = (("realU", "real"), ("realT1f", "realT"), ("scsemi", "scsemi"), ("prev4k", "prev"), ("floor", "floor"))
for ax, met, yl in zip(axs[:2], ("T_marker_any", "mhc_shift_case_minus_ctrl_sdA"), ("monocytes with any T-marker UMI", "MHC-II score shift, case − control (real-A SD)")):
    for i, coh in enumerate(("sepsis", "covid")):
        for j, (s, c) in enumerate(srcs):
            v = MY[(MY.cohort == coh) & (MY.source == s)][met].item(); ax.bar(i + j * 0.16, v, 0.15, color=C[c]); ax.text(i + j * 0.16, v + (0.004 if v >= 0 else -0.08), f"{v:.2f}", ha="center", fontsize=7)
    ax.set_xticks([0.32, 1.32], ["Sepsis", "COVID"]); ax.set_ylabel(yl); grid(ax); ax.axhline(0, color="#3E4C59", lw=0.6)
axs[0].set_title("T/NK co-detection in monocytes", loc="left", fontsize=10.5); axs[1].set_title("Monocyte HLA class II in disease", loc="left", fontsize=10.5)
ax = axs[2]; d = G[G.model == "floor"].head(20).iloc[::-1]
ax.barh(range(len(d)), d.residual, color=[C["floor"] if r > 0 else "#ABD9E9" for r in d.residual]); ax.barh(range(len(d)), d.null_residual, color="#9AA5B1", height=0.35)
ax.set_yticks(range(len(d)), [f"{g} ({c})" if c != "other" else g for g, c in zip(d.gene, d.category)], fontsize=7.5); ax.axvline(0, color="#3E4C59", lw=0.6)
ax.set_xlabel("mean log1p CP10k: floor − real B (same depth); grey = real A − real B"); ax.set_title("Sepsis monocytes: largest floor residuals", loc="left", fontsize=10.5)
legend(fig, ("real", "realT", "scsemi", "prev", "floor")); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "s3_myeloid")
# 7. integration: iLISI per method + per-sample iLISI
I = pd.concat([pd.read_csv(f"{R}/{c}_integration_methods.csv") for c in ("sepsis", "covid")]).drop_duplicates(["cohort", "pair", "method"], keep="last")
BI = pd.read_csv(f"{R}/bootstrap_identity_mixing.csv"); PS = BI[BI.analysis.str.startswith("per-sample iLISI")]
fig, axs = plt.subplots(1, 2, figsize=(14, 3.9), gridspec_kw={"width_ratios": [1.6, 1]})
ax = axs[0]; k = 0; xt, xl = [], []
for coh in ("sepsis", "covid"):
    for meth in ("Harmony", "Scanorama", "scVI", "BBKNN (graph)"):
        for j, (s, c) in enumerate((("realU", "real"), ("scsemi", "scsemi"), ("prev4k", "prev"), ("floor", "floor"))):
            v = I[(I.cohort == coh) & (I.method == meth) & (I.pair == s)]
            if len(v): ax.bar(k + j * 0.2, v.iLISI.item(), 0.19, color=C[c])
        xt.append(k + 0.3); xl.append(f"{meth.replace(' (graph)', '')}\n{coh.capitalize()}"); k += 1
    k += 0.4
ax.set_xticks(xt, xl, fontsize=8); ax.set_ylim(1, 2); ax.set_ylabel("iLISI (2 = perfect mixing)"); grid(ax); ax.set_title("Integration robustness: four methods", loc="left", fontsize=10.5)
ax = axs[1]
for i, coh in enumerate(("sepsis", "covid")):
    for j, (s, c) in enumerate((("realU", "real"), ("scsemi", "scsemi"), ("prev4k", "prev"), ("floor", "floor"))):
        r = PS[(PS.analysis.str.endswith(coh)) & (PS.metric == s)].iloc[0]; ax.bar(i + j * 0.2, r.estimate, 0.19, color=C[c])
        ax.errorbar(i + j * 0.2, r.estimate, yerr=[[r.estimate - r.lo], [r.hi - r.estimate]], color="#1F2933", lw=1, capsize=2, fmt="none")
ax.set_xticks([0.3, 1.3], ["Sepsis", "COVID"]); ax.set_ylim(1, 2); ax.set_ylabel("per-sample iLISI (median over samples)"); grid(ax)
ax.set_title("Per-sample mixing, generated cells (95 % CI)", loc="left", fontsize=10.5)
legend(fig, ("real", "scsemi", "prev", "floor")); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "s3_integration")
# 8. UMAP scSemi vs floor (post-Harmony)
fig, axs = plt.subplots(2, 2, figsize=(10, 9.5))
for i, coh in enumerate(("sepsis", "covid")):
    for j, s in enumerate(("scsemi", "floor")):
        z = np.load(f"{TR}/{coh}_{s}_umap.npz", allow_pickle=True); U = z["umap_post"]; src = z["src"]; ax = axs[i, j]; o = np.random.default_rng(0).permutation(len(src))
        ax.scatter(U[o, 0], U[o, 1], s=0.3, c=np.where(src[o] == 0, "#CBD2D9", C[s]), rasterized=True); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{coh.capitalize()}: real A (grey) + {NM[s]}", loc="left", fontsize=10)
        for sp_ in ax.spines.values(): sp_.set_visible(False)
fig.tight_layout(); save(fig, "s3_umap")
print("figures saved")
