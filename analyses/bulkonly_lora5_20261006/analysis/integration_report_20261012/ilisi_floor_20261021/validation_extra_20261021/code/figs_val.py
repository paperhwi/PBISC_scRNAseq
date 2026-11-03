import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
                     "axes.edgecolor": "#3E4C59", "axes.labelcolor": "#1F2933", "xtick.color": "#3E4C59", "ytick.color": "#3E4C59"})
R = "/mnt/user-data/outputs/val/results"; F = "/mnt/user-data/outputs/val/figs"
C = {"real": "#3E4C59", "realT": "#AEB7C2", "floor": "#2B83BA", "prev": "#D7191C", "chance": "#9AA5B1"}
def save(fig, n): fig.savefig(f"{F}/{n}.png", dpi=220, bbox_inches="tight", facecolor="white"); fig.savefig(f"{F}/{n}.pdf", bbox_inches="tight"); plt.close(fig)
def grid(ax): ax.yaxis.grid(True, color="#EEF0F2", lw=0.6); ax.set_axisbelow(True)
# 1. M6 + donor fingerprint
M = pd.read_csv(f"{R}/m6_disease_effect.csv"); M = M[M.contrast == "case vs control"]; FP = pd.read_csv(f"{R}/covid_donor_fingerprint.csv")
fig, axs = plt.subplots(1, 4, figsize=(15, 3.9), gridspec_kw={"width_ratios": [1, 1, 1, 1.1]})
for ax, met, yl in zip(axs[:3], ["logFC_r", "slope_on_DEG", "sign_agree_DEG"], ["logFC r (real vs generated)", "slope on real DEGs", "sign agreement on real DEGs"]):
    k = 0; xt, xl = [], []
    for coh in ("sepsis", "covid"):
        for band in ("rank <= 2,000", "rank > 2,000"):
            for j, s in enumerate(("prev", "floor")):
                d = M[(M.cohort == coh) & (M.genes == band) & (M.source == s)][met]
                ax.bar(k + j * 0.38, d.median(), 0.36, color=C[s]); ax.scatter(np.full(len(d), k + j * 0.38), d, s=9, color="#1F2933", zorder=3)
            xt.append(k + 0.19); xl.append(f"{coh.capitalize()}\n{band.replace('rank ', '')}"); k += 1
    ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel(yl); grid(ax); ax.set_ylim(0, 1.05 if met != "slope_on_DEG" else 0.9)
axs[0].set_title("Disease effect per cell type (bar = median of 5 types)", loc="left", fontsize=10.5)
ax = axs[3]; d = FP[FP.gene_set != "top-300 donor-stable genes"]
for i, gs in enumerate(["top-100 donor-stable genes", "all expressed genes"]):
    for j, s in enumerate(("prev", "floor")):
        v = d[(d.gene_set == gs) & (d.source == s)]; ax.bar(i + j * 0.38, v.top1_acc.item(), 0.36, color=C[s]); ax.text(i + j * 0.38, v.top1_acc.item() + 0.01, f"{v.top1_acc.item():.2f}", ha="center", fontsize=9)
ax.axhline(1 / 32, color=C["chance"], ls="--", lw=1); ax.text(-0.15, 1 / 32 + 0.012, "chance 1/32", fontsize=8.5, color=C["chance"], ha="left")
ax.set_xticks([0.19, 1.19], ["donor-stable\ngenes (top 100)", "all expressed\ngenes"], fontsize=8.5); ax.set_ylabel("top-1 donor identification"); grid(ax); ax.set_ylim(0, 0.5)
ax.set_title("COVID: 32 donors, generated vs real", loc="left", fontsize=10.5)
h = [plt.Rectangle((0, 0), 1, 1, color=C[s], label=l) for s, l in (("prev", "Previous K = 4,000"), ("floor", "Floor read-out"))]; fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.06))
fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "val_m6_donor")
# 2. C2ST + density/coverage
fig, axs = plt.subplots(1, 3, figsize=(15, 3.9), gridspec_kw={"width_ratios": [1.5, 1, 1]})
ax = axs[0]; k = 0; xt, xl = [], []
for coh in ("sepsis", "covid"):
    T = pd.read_csv(f"{R}/{coh}_c2st.csv")
    for lin in ("T", "NK", "B", "Myeloid"):
        for j, (s, c) in enumerate((("realU", C["real"]), ("realT1f", C["realT"]), ("prev4k", C["prev"]), ("floor", C["floor"]))):
            ax.bar(k + j * 0.2, T[(T.lineage == lin) & (T.source == s)].c2st_auc.item(), 0.19, color=c)
        xt.append(k + 0.3); xl.append(f"{lin}\n{coh.capitalize()}"); k += 1
    k += 0.4
ax.axhline(0.5, color=C["chance"], ls="--", lw=1); ax.set_ylim(0.3, 1.02); ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel("classifier AUC vs real A (0.5 = indistinguishable)"); grid(ax)
ax.set_title("Classifier two-sample test (sample-grouped CV)", loc="left", fontsize=10.5)
for ax, met in zip(axs[1:], ("coverage", "density")):
    k = 0; xt, xl = [], []
    for coh in ("sepsis", "covid"):
        T = pd.read_csv(f"{R}/{coh}_density_coverage.csv"); T = T[T.lineage == "all"]
        for sp_ in ("PCA (no correction)", "Harmony"):
            for j, (s, c) in enumerate((("realU", C["real"]), ("prev4k", C["prev"]), ("floor", C["floor"]))):
                ax.bar(k + j * 0.27, T[(T.space == sp_) & (T.source == s)][met].item(), 0.25, color=c)
            xt.append(k + 0.27); xl.append(f"{coh.capitalize()}\n{'PCA' if sp_.startswith('PCA') else 'Harmony'}"); k += 1
    ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel(met); grid(ax); ax.set_title(f"{met.capitalize()} (k = 5, real A manifold)", loc="left", fontsize=10.5)
h = [plt.Rectangle((0, 0), 1, 1, color=c, label=l) for c, l in ((C["real"], "Real B (other patients)"), (C["realT"], "Real B depth-matched"), (C["prev"], "Previous K = 4,000"), (C["floor"], "Floor read-out"))]
fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.06)); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "val_c2st_dc")
# 3. co-expression + LR
fig, axs = plt.subplots(1, 3, figsize=(15, 3.9), gridspec_kw={"width_ratios": [1.4, 1, 1]})
ax = axs[0]; k = 0; xt, xl = [], []
for coh in ("sepsis", "covid"):
    T = pd.read_csv(f"{R}/{coh}_coexpression.csv")
    for band in ("top 300 genes", "ranks 1,001-1,300"):
        for j, (s, c) in enumerate((("realU", C["real"]), ("shuffle", "#E4E7EB"), ("prev4k", C["prev"]), ("floor", C["floor"]))):
            v = T[(T.genes == band) & (T.source == s)].corr_of_corr; ax.bar(k + j * 0.2, v.median(), 0.19, color=c); ax.scatter(np.full(len(v), k + j * 0.2), v, s=8, color="#1F2933", zorder=3)
        xt.append(k + 0.3); xl.append(f"{coh.capitalize()}\n{band}"); k += 1
ax.set_xticks(xt, xl, fontsize=8.5); ax.set_ylabel("r of gene–gene correlations vs real A"); grid(ax); ax.set_title("Co-expression within lineage (dots = 4 lineages)", loc="left", fontsize=10.5)
for ax, met, yl in zip(axs[1:], ("detected_jaccard", "recall_top200_realA"), ("detected-interaction Jaccard vs real A", "recall of real A top-200 interactions")):
    k = 0
    for coh in ("sepsis", "covid"):
        T = pd.read_csv(f"{R}/{coh}_lr_recovery.csv").set_index("source")
        for j, (s, c) in enumerate((("realU", C["real"]), ("realT1f", C["realT"]), ("prev4k", C["prev"]), ("floor", C["floor"]))):
            ax.bar(k + j * 0.2, T.loc[s, met], 0.19, color=c); ax.text(k + j * 0.2, T.loc[s, met] + 0.01, f"{T.loc[s, met]:.2f}", ha="center", fontsize=8)
        k += 1
    ax.set_xticks([0.3, 1.3], ["Sepsis", "COVID"]); ax.set_ylabel(yl); grid(ax); ax.set_ylim(0, 1.08); ax.set_title("Ligand–receptor (LIANA consensus)", loc="left", fontsize=10.5)
h = [plt.Rectangle((0, 0), 1, 1, color=c, label=l) for c, l in ((C["real"], "Real B"), ("#E4E7EB", "Real B, genes shuffled"), (C["realT"], "Real B depth-matched"), (C["prev"], "Previous K = 4,000"), (C["floor"], "Floor read-out"))]
fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.06)); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "val_coexpr_lr")
# 4. claims + transfer-free rare
fig, axs = plt.subplots(1, 2, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1, 1.3]})
ax = axs[0]; CL = pd.concat([pd.read_csv(f"{R}/{c}_claims.csv") for c in ("sepsis", "covid")]); CL = CL[(CL.auc_realA >= 0.6) & (CL.claim != "ESKD-IFN B cell")]
cl = CL.claim.unique(); x = np.arange(len(cl))
for j, (s, c) in enumerate((("realU", C["real"]), ("prev4k", C["prev"]), ("floor", C["floor"]))):
    ax.bar(x + (j - 1) * 0.27, [CL[(CL.claim == k_) & (CL.source == s)].delta_pos_frac.item() for k_ in cl], 0.25, color=c)
ax.axhline(0, color="#3E4C59", lw=0.8); ax.set_ylim(-0.32, 0.38); ax.set_xticks(x, [k_.replace("ESKD-IFN ", "IFN ") for k_ in cl], rotation=30, ha="right", fontsize=8.5)
ax.set_ylabel("Δ classifier-positive fraction\n(case − control, sample means)"); grid(ax); ax.set_title("Paper claims, same classifier on every source (real-A AUC ≥ 0.6)", loc="left", fontsize=10.5)
ax = axs[1]; RR = pd.concat([pd.read_csv(f"{R}/{c}_rare_transfer_free.csv") for c in ("sepsis", "covid")]); RR = RR[RR.cv_auc_realA >= 0.8]
RR["ratio"] = RR.gen_estimate_pct / RR.realB_estimate_pct; W = RR.pivot_table(index=["cohort", "subset"], columns="model_depth", values="ratio").reset_index()
W = W.sort_values(["cohort", "floor"]); y = np.arange(len(W))
ax.scatter(np.clip(W.prev4k, 0.01, 20), y, color=C["prev"], s=30, label="Previous K = 4,000", zorder=3); ax.scatter(np.clip(W.floor, 0.01, 20), y, color=C["floor"], s=30, label="Floor read-out", zorder=3)
for i, r in enumerate(W.itertuples()): ax.plot([np.clip(r.prev4k, 0.01, 20), np.clip(r.floor, 0.01, 20)], [i, i], color="#C9CED6", lw=1, zorder=1)
ax.axvspan(0.5, 2, color="#EAF3FA", zorder=0); ax.axvline(1, color="#9AA5B1", lw=0.8, ls=":"); ax.set_xscale("log"); ax.set_xlim(0.01, 20)
ax.set_yticks(y, [f"{r.subset} ({r.cohort[:3]})" for r in W.itertuples()], fontsize=8); ax.set_xlabel("generated / real B share, same classifier (shaded = 0.5–2×)")
ax.set_title("Rare states without label transfer (classifier AUC ≥ 0.8)", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5, loc="lower right")
h = [plt.Rectangle((0, 0), 1, 1, color=C[s], label=l) for s, l in (("real", "Real B"), ("prev", "Previous K = 4,000"), ("floor", "Floor read-out"))]; axs[0].legend(handles=h, frameon=False, fontsize=8.5, loc="upper right", ncol=3)
fig.tight_layout(); save(fig, "val_claims_rare"); print("figs ok")
# 5. integration methods (run after integ.py)
if all(os.path.exists(f"{R}/{c}_integration_methods.csv") for c in ("sepsis", "covid")):
    I = pd.concat([pd.read_csv(f"{R}/{c}_integration_methods.csv") for c in ("sepsis", "covid")])
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.0), gridspec_kw={"width_ratios": [1.6, 1]})
    meths = ["PCA (none)", "Harmony", "Scanorama", "scVI", "BBKNN (graph)"]
    for ax, met, yl in ((axs[0], "iLISI", "median iLISI of source (2 = mixed)"), (axs[1], "label_ASW", "six-type label ASW (higher = types kept apart)")):
        k = 0; xt, xl = [], []
        for coh in ("sepsis", "covid"):
            for m_ in (meths if met == "iLISI" else meths[:4]):
                for j, (s, c) in enumerate((("realU", C["real"]), ("prev4k", C["prev"]), ("floor", C["floor"]))):
                    v = I[(I.cohort == coh) & (I.method == m_) & (I.pair == s)][met].item(); ax.bar(k + j * 0.27, v, 0.25, color=c)
                xt.append(k + 0.27); xl.append(f"{m_.replace(' (graph)', '').replace(' (none)', '')}\n{coh.capitalize()}"); k += 1
            k += 0.4
        ax.set_xticks(xt, xl, fontsize=8); ax.set_ylabel(yl); grid(ax)
        ax.set_ylim((1, 2) if met == "iLISI" else (0.4, 0.75))
    axs[0].set_title("Real A + source, 20k + 20k cells, same HVG (BBKNN = graph inverse Simpson)", loc="left", fontsize=10.5); axs[1].set_title("Biology kept", loc="left", fontsize=10.5)
    h = [plt.Rectangle((0, 0), 1, 1, color=c, label=l) for c, l in ((C["real"], "Real B"), (C["prev"], "Previous K = 4,000"), (C["floor"], "Floor read-out"))]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.06)); fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "val_integration"); print("integration fig ok")
