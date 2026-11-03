"""One-table summary of the extra validation (sepsis / COVID): real control, previous K = 4,000, floor."""
import os, numpy as np, pandas as pd
R = "/mnt/user-data/outputs/val/results"; rd = lambda n: pd.read_csv(f"{R}/{n}"); rows = []
def add(area, metric, better, real, prev, floor): rows.append({"area": area, "metric": metric, "better": better, "real_control": real, "prev_K4000": prev, "floor": floor})
f2 = lambda a, b: f"{a:.2f} / {b:.2f}"
FP = rd("covid_donor_fingerprint.csv").set_index(["source", "gene_set"])
add("A patient", "COVID donor top-1 identification, donor-stable genes (chance 0.03)", "higher", "–", f"{FP.loc[('prev','top-100 donor-stable genes'),'top1_acc']:.2f}", f"{FP.loc[('floor','top-100 donor-stable genes'),'top1_acc']:.2f}")
add("A patient", "COVID own-donor r, all expressed genes (other donors ≈ 0)", "higher", "–", f"{FP.loc[('prev','all expressed genes'),'own_r_median']:.2f}", f"{FP.loc[('floor','all expressed genes'),'own_r_median']:.2f}")
M = rd("m6_disease_effect.csv"); M = M[M.contrast == "case vs control"]; g = M.groupby(["genes", "source", "cohort"])
for band in ("all", "rank > 2,000"):
    for met, lab in (("logFC_r", "logFC r"), ("slope_on_DEG", "slope on real DEGs"), ("sign_agree_DEG", "sign agreement on real DEGs")):
        v = g[met].median(); add("B disease effect", f"{lab}, {band} genes (median of 5 types)", "higher", "1", f2(v[(band, 'prev', 'sepsis')], v[(band, 'prev', 'covid')]), f2(v[(band, 'floor', 'sepsis')], v[(band, 'floor', 'covid')]))
P = rd("m6_precision_top100.csv"); P = P[P.contrast == "case vs control"].groupby(["source", "cohort"]).precision_top100_gen.median()
add("B disease effect", "precision of generated top-100 DEGs", "higher", "–", f2(P[('prev','sepsis')], P[('prev','covid')]), f2(P[('floor','sepsis')], P[('floor','covid')]))
C2 = pd.concat([rd(f"{c}_c2st.csv") for c in ("sepsis", "covid")]).groupby(["source", "cohort"]).c2st_auc.median()
add("C realism", "classifier AUC vs real A (0.5 = indistinguishable; median of 4 lineages)", "lower", f2(C2[('realU','sepsis')], C2[('realU','covid')]), f2(C2[('prev4k','sepsis')], C2[('prev4k','covid')]), f2(C2[('floor','sepsis')], C2[('floor','covid')]))
DC = pd.concat([rd(f"{c}_density_coverage.csv") for c in ("sepsis", "covid")]); DC = DC[DC.lineage == "all"].set_index(["space", "source", "cohort"])
for sp_ in ("PCA (no correction)", "Harmony"):
    for met in ("coverage", "density"):
        q = lambda s: f2(DC.loc[(sp_, s, 'sepsis'), met], DC.loc[(sp_, s, 'covid'), met])
        add("C realism", f"{met}, {sp_}", "closer to real", q("realU"), q("prev4k"), q("floor"))
CO = pd.concat([rd(f"{c}_coexpression.csv") for c in ("sepsis", "covid")]).groupby(["genes", "source", "cohort"]).corr_of_corr.median()
for band in ("top 300 genes", "ranks 1,001-1,300"):
    q = lambda s: f2(CO[(band, s, 'sepsis')], CO[(band, s, 'covid')])
    add("D co-expression", f"r of gene–gene correlations vs real A, {band}", "higher", q("realU"), q("prev4k"), q("floor"))
LR = pd.concat([rd(f"{c}_lr_recovery.csv") for c in ("sepsis", "covid")]).set_index(["source", "cohort"]); LD = pd.concat([rd(f"{c}_lr_disease.csv") for c in ("sepsis", "covid")]).set_index(["source", "cohort"])
for met, lab in (("detected_jaccard", "detected-interaction Jaccard"), ("recall_top200_realA", "recall of real A top-200 interactions"), ("spearman_score_on_realA_detected", "Spearman of interaction scores")):
    q = lambda s: f2(LR.loc[(s, 'sepsis'), met], LR.loc[(s, 'covid'), met]); add("E ligand–receptor", lab, "higher", q("realU"), q("prev4k"), q("floor"))
add("E ligand–receptor", "disease change of interaction scores, Spearman", "higher", "1", f2(LD.loc[('prev4k','sepsis'),'disease_spearman_all_detected'], LD.loc[('prev4k','covid'),'disease_spearman_all_detected']),
    f2(LD.loc[('floor','sepsis'),'disease_spearman_all_detected'], LD.loc[('floor','covid'),'disease_spearman_all_detected']))
CL = pd.concat([rd(f"{c}_claims.csv") for c in ("sepsis", "covid")]); CL = CL[(CL.auc_realA >= 0.6) & (CL.claim != "ESKD-IFN B cell")]
W = CL.pivot_table(index=["cohort", "claim"], columns="source", values="delta_pos_frac")
for coh in ("sepsis", "covid"):
    w = W.loc[coh]; s_ = lambda s: f"{int((np.sign(w[s]) == np.sign(w['realU'])).sum())}/{len(w)} sign, size {np.median(w[s] / w['realU']):.2f}×"
    add("F paper claims", f"{coh}: same-classifier contrast vs real B", "closer to real", "–", s_("prev4k"), s_("floor"))
RR = pd.concat([rd(f"{c}_rare_transfer_free.csv") for c in ("sepsis", "covid")]); RR = RR[RR.cv_auc_realA >= 0.8]; RR["ok"] = (RR.gen_estimate_pct / RR.realB_estimate_pct).between(0.5, 2)
k = RR.groupby(["model_depth", "cohort"]).ok.agg(["sum", "size"])
add("G rare (no transfer)", "subsets within 0.5–2× of real B, same classifier (AUC ≥ 0.8)", "higher", "–",
    f"{k.loc[('prev4k','sepsis'),'sum']}/{k.loc[('prev4k','sepsis'),'size']} / {k.loc[('prev4k','covid'),'sum']}/{k.loc[('prev4k','covid'),'size']}", f"{k.loc[('floor','sepsis'),'sum']}/{k.loc[('floor','sepsis'),'size']} / {k.loc[('floor','covid'),'sum']}/{k.loc[('floor','covid'),'size']}")
if os.path.exists(f"{R}/covid_integration_methods.csv"):
    I = pd.concat([rd(f"{c}_integration_methods.csv") for c in ("sepsis", "covid")]).set_index(["method", "pair", "cohort"])
    for meth in ("Harmony", "Scanorama", "scVI", "BBKNN (graph)"):
        q = lambda s: f2(I.loc[(meth, s, 'sepsis'), 'iLISI'], I.loc[(meth, s, 'covid'), 'iLISI']); add("H integration", f"iLISI, {meth}", "higher", q("realU"), q("prev4k"), q("floor"))
S = pd.DataFrame(rows); S.to_csv(f"{R}/scorecard.csv", index=False); pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90); print(S.to_string())
