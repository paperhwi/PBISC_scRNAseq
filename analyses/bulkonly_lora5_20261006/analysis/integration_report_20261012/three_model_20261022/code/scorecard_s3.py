"""Three-model scorecard (2026-10-22): scSemiProfiler (strict 5 pairs, draw0) / previous PBISC K = 4,000 (Poisson read-out) / PBISC floor read-out,
with the real control (real half B, other patients). Each cell = sepsis / COVID. Usage: python scorecard_s3.py"""
import os, json, numpy as np, pandas as pd
R = "/mnt/user-data/outputs/s3/results"; TR = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"; rd = lambda n: pd.read_csv(f"{R}/{n}"); rows = []
def add(area, metric, better, real, scsemi, prev, floor, note=""):
    rows.append({"area": area, "metric": metric, "better": better, "real_control": real, "scSemiProfiler": scsemi, "prev_K4000": prev, "floor": floor, "note": note})
f2 = lambda a, b: f"{a:.2f} / {b:.2f}"
both = lambda fn: f2(fn("sepsis"), fn("covid"))
# A integration / matching (big_joint summaries)
J = {(c, s): json.load(open(f"{TR}/{c}_{s}_summary.json")) for c in ("sepsis", "covid") for s in ("realU", "floor", "prev4k", "scsemi")}
for key, lab in (("iLISI_post", "iLISI after pairwise Harmony (max 2)"), ("src_cells_in_mixed_clusters", "share of cells in mixed joint clusters (20–80 %)"),
                 ("cluster_freq_pearson", "joint-cluster frequency r with real A"), ("six_agree_with_own_label", "kNN six-type = own label")):
    g = lambda s: both(lambda c: J[(c, s)][key]); add("A integration", lab, "higher", g("realU"), g("scsemi"), g("prev4k"), g("floor"))
BI = rd("bootstrap_identity_mixing.csv"); PS = BI[BI.analysis.str.startswith("per-sample iLISI")].assign(c=lambda d: d.analysis.str.split().str[-1]).set_index(["c", "metric"])
g = lambda s: both(lambda c: PS.loc[(c, s), "estimate"]); ci = lambda s: "; ".join(f"{c} [{PS.loc[(c, s), 'lo']:.2f}, {PS.loc[(c, s), 'hi']:.2f}]" for c in ("sepsis", "covid"))
add("A integration", "per-sample iLISI of generated cells, median over samples [95 % CI]", "higher", g("realU"), g("scsemi"), g("prev4k"), g("floor"), "floor " + ci("floor"))
I = pd.concat([rd(f"{c}_integration_methods.csv") for c in ("sepsis", "covid")]).drop_duplicates(["cohort", "pair", "method"], keep="last").set_index(["method", "pair", "cohort"])
for meth in ("Harmony", "Scanorama", "scVI", "BBKNN (graph)"):
    q = lambda s: f2(I.loc[(meth, s, 'sepsis'), 'iLISI'], I.loc[(meth, s, 'covid'), 'iLISI']) if (meth, s, 'sepsis') in I.index else "–"
    add("A integration", f"iLISI, {meth} (20k + 20k)", "higher", q("realU"), q("scsemi"), q("prev4k"), q("floor"))
# B realism
C2 = pd.concat([rd(f"{c}_c2st.csv") for c in ("sepsis", "covid")]).groupby(["source", "cohort"]).c2st_auc.median()
q = lambda s: f2(C2[(s, 'sepsis')], C2[(s, 'covid')])
add("B realism", "classifier AUC vs real A (0.5 = indistinguishable; median of 4 lineages)", "lower", q("realU"), q("scsemi"), q("prev4k"), q("floor"),
    f"real B at model depth: scSemi {q('realT1s')}, prev {q('realT1p')}, floor {q('realT1f')}; scSemi floor-rule copy {q('scsemiF')}")
DC = pd.concat([rd(f"{c}_density_coverage.csv") for c in ("sepsis", "covid")]); DC = DC[DC.lineage == "all"].set_index(["space", "source", "cohort"])
for sp_ in ("PCA (no correction)", "Harmony"):
    for met in ("coverage", "density"):
        q = lambda s: f2(DC.loc[(sp_, s, 'sepsis'), met], DC.loc[(sp_, s, 'covid'), met]); add("B realism", f"{met}, {sp_}", "closer to real", q("realU"), q("scsemi"), q("prev4k"), q("floor"))
MY = pd.concat([rd(f"{c}_myeloid_identity.csv") for c in ("sepsis", "covid")]).set_index(["source", "cohort"])
q = lambda s: f2(MY.loc[(s, 'sepsis'), 'median_genes'], MY.loc[(s, 'covid'), 'median_genes'])
add("B realism", "genes detected per monocyte (median)", "closer to real", q("realU"), q("scsemi"), q("prev4k"), q("floor"), f"real B at model depth: scSemi {q('realT1s')}, floor {q('realT1f')}")
CO = pd.concat([rd(f"{c}_coexpression.csv") for c in ("sepsis", "covid")]).groupby(["genes", "source", "cohort"]).corr_of_corr.median()
for band in ("top 300 genes", "ranks 1,001-1,300"):
    q = lambda s: f2(CO[(band, s, 'sepsis')], CO[(band, s, 'covid')])
    add("C co-expression", f"r of gene–gene correlations vs real A, {band}", "higher", q("realU"), q("scsemi"), q("prev4k"), q("floor"), f"gene-shuffled real B {q('shuffle')}")
# D disease effect
M = rd("m6_three_models.csv"); M = M[M.contrast == "case vs control"]; g = M.groupby(["genes", "source", "cohort"])
for band in ("all", "rank > 2,000"):
    for met, lab in (("logFC_r", "logFC r"), ("slope_on_DEG", "slope on real DEGs"), ("sign_agree_DEG", "sign agreement on real DEGs")):
        v = g[met].median(); q = lambda s: f2(v[(band, s, 'sepsis')], v[(band, s, 'covid')])
        add("D disease effect", f"{lab}, {band} genes (median of 5 types)", "higher", "1", q("scsemi"), q("prev"), q("floor"), f"bulk broadcast {q('broadcast')}")
GT = pd.read_csv(f"{R}/m6_three_models_gene_tables.csv.gz"); prec = {}
for k, d in GT.groupby("key"):
    coh, con, t, s = k.split("|")
    if con != "case vs control": continue
    d = d[(d.real_mean > 1) & np.isfinite(d.gen_t)]; top = d.assign(a=d.gen_t.abs()).nlargest(100, "a")
    prec.setdefault((coh, s), []).append(((top.real_q < 0.05) & (np.sign(top.real_logFC) == np.sign(top.gen_logFC))).mean())
pq = lambda s: f2(np.median(prec[('sepsis', s)]), np.median(prec[('covid', s)]))
add("D disease effect", "precision of generated top-100 DEGs (real FDR < 0.05, same sign)", "higher", "–", pq("scsemi"), pq("prev"), pq("floor"), f"bulk broadcast {pq('broadcast')}")
PR = rd("m6_broadcast_partial.csv"); PR = PR[PR.contrast == "case vs control"].groupby(["genes", "source", "cohort"])
PP = PR.partial_r_model_given_broadcast.median(); RB = PR.r_broadcast_real.median()
for band in ("all", "rank > 2,000"):
    q = lambda s: f2(PP[(band, s, 'sepsis')], PP[(band, s, 'covid')])
    add("D disease effect", f"partial r of model logFC given bulk broadcast, {band} genes", "higher", "–", q("scsemi"), q("prev"), q("floor"),
        f"broadcast alone r {f2(RB[(band, 'floor', 'sepsis')], RB[(band, 'floor', 'covid')])}; single-cell information beyond the bulk input")
CI = rd("m6_bootstrap_ci.csv")
def cis(m, band, comp):
    return "; ".join(f"{c} {r.med:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}]" for c in ("sepsis", "covid") for _, r in CI[(CI.cohort == c) & (CI.genes == band) & (CI.metric == m) & (CI.comparison == comp)].iterrows())
add("D disease effect", "bootstrap Δ logFC r, rank > 2,000 (median [95 % CI], sepsis; covid)", "> 0", "–", "floor − scSemi: " + cis("logFC_r", "rank > 2,000", "floor - scsemi"), "", "floor − prev: " + cis("logFC_r", "rank > 2,000", "floor - prev"))
add("D disease effect", "bootstrap Δ logFC r vs bulk broadcast, all genes (median [95 % CI], sepsis; covid)", "> 0", "–", "", "prev − broadcast: " + cis("logFC_r", "all", "prev - broadcast"), "floor − broadcast: " + cis("logFC_r", "all", "floor - broadcast"))
add("D disease effect", "bootstrap Δ sign agreement vs bulk broadcast, all genes", "> 0", "–", "", "prev − broadcast: " + cis("sign_agree_DEG", "all", "prev - broadcast"), "floor − broadcast: " + cis("sign_agree_DEG", "all", "floor - broadcast"))
# E composition
CP = rd("composition_shift.csv"); CP = CP[CP.contrast == "case vs control"].set_index(["level", "source", "cohort"])
for lvl, lab in (("six", "six-type"), ("fine", "fine (transferred)")):
    for met, ml, bt in (("TV_case", "group-mean composition error, case (TV)", "lower"), ("shift_r", "disease shift r vs real", "higher"), ("log2ratio_r", "disease log2 ratio r vs real", "higher")):
        q = lambda s: f2(CP.loc[(lvl, s, 'sepsis'), met], CP.loc[(lvl, s, 'covid'), met])
        add("E composition", f"{lab}: {ml}", bt, q("realB_tr") if lvl == "fine" else "–", q("scsemi"), q("prev"), q("floor"))
# F ligand-receptor
LR = pd.concat([rd(f"{c}_lr_recovery.csv") for c in ("sepsis", "covid")]).set_index(["source", "cohort"]); LD = pd.concat([rd(f"{c}_lr_disease.csv") for c in ("sepsis", "covid")]).set_index(["source", "cohort"])
for met, lab in (("detected_jaccard", "detected-interaction Jaccard"), ("recall_top200_realA", "recall of real A top-200 interactions"), ("spearman_score_on_realA_detected", "Spearman of interaction scores")):
    q = lambda s: f2(LR.loc[(s, 'sepsis'), met], LR.loc[(s, 'covid'), met]); add("F ligand–receptor", lab, "higher", q("realU"), q("scsemi"), q("prev4k"), q("floor"), f"real B at scSemi depth {q('realT1s')}")
q = lambda s: f2(LD.loc[(s, 'sepsis'), 'disease_spearman_all_detected'], LD.loc[(s, 'covid'), 'disease_spearman_all_detected'])
add("F ligand–receptor", "disease change of interaction scores, Spearman", "higher", "1", q("scsemi"), q("prev4k"), q("floor"))
# G claims / rare
CL = pd.concat([rd(f"{c}_claims.csv") for c in ("sepsis", "covid")]); CL = CL[(CL.auc_realA >= 0.6) & (CL.claim != "ESKD-IFN B cell")]
W = CL.pivot_table(index=["cohort", "claim"], columns="source", values="delta_pos_frac")
for coh in ("sepsis", "covid"):
    w = W.loc[coh]; s_ = lambda s: f"{int((np.sign(w[s]) == np.sign(w['realU'])).sum())}/{len(w)} sign, size {np.nanmedian(w[s] / w['realU']):.2f}×"
    add("G paper claims", f"{coh}: same-classifier contrast vs real B", "closer to real", "–", s_("scsemi"), s_("prev4k"), s_("floor"),
        "scSemi sepsis: marker-positive fraction ~0 in both groups (sparse cells); SEP-3 parent has no scSemi control sample" if coh == "sepsis" else "")
RR = pd.concat([rd(f"{c}_rare_transfer_free.csv") for c in ("sepsis", "covid")]); RR = RR[RR.cv_auc_realA >= 0.8]; RR["ok"] = (RR.gen_estimate_pct / RR.realB_estimate_pct).between(0.5, 2)
k = RR.groupby(["model_depth", "cohort"]).ok.agg(["sum", "size"]); kk = lambda s: f"{k.loc[(s,'sepsis'),'sum']}/{k.loc[(s,'sepsis'),'size']} / {k.loc[(s,'covid'),'sum']}/{k.loc[(s,'covid'),'size']}"
add("G rare (no transfer)", "subsets within 0.5–2× of real B, same classifier (AUC ≥ 0.8)", "higher", "–", kk("scsemi"), kk("prev4k"), kk("floor"))
# H patient identity
FP = BI[BI.analysis.str.startswith("COVID donor")].set_index("metric"); SI = BI[BI.analysis.str.startswith("scSemi sample")]
add("H patient identity", "COVID donor top-1 (32 donors, mostly other time points; chance 0.03) [95 % CI]", "higher", "–", "n/a (see note)",
    f"{FP.loc['hit_prev','estimate']:.2f} [{FP.loc['hit_prev','lo']:.2f}, {FP.loc['hit_prev','hi']:.2f}]", f"{FP.loc['hit_floor','estimate']:.2f} [{FP.loc['hit_floor','lo']:.2f}, {FP.loc['hit_floor','hi']:.2f}]",
    "scSemi same-sample self-ID (easier task): " + "; ".join(f"{r.analysis.split(',')[1].split('(')[0].strip()} {r.estimate:.2f} [{r.lo:.2f}, {r.hi:.2f}]" for _, r in SI.iterrows()))
add("H patient identity", "floor − prev donor top-1 [95 % CI]", "> 0", "–", "", "", f"{FP.loc['dhit','estimate']:+.2f} [{FP.loc['dhit','lo']:+.2f}, {FP.loc['dhit','hi']:+.2f}]")
# I myeloid
q = lambda s, m: f2(MY.loc[(s, 'sepsis'), m], MY.loc[(s, 'covid'), m])
add("I myeloid", "monocytes with any T/NK marker UMI", "closer to real", q("realT1f", "T_or_NK_any"), q("scsemi", "T_or_NK_any"), q("prev4k", "T_or_NK_any"), q("floor", "T_or_NK_any"), "real = real B at floor depth")
add("I myeloid", "monocytes with any T marker UMI (CD3E/CD3D/CD2/IL7R/TRAC)", "closer to real", q("realT1f", "T_marker_any"), q("scsemi", "T_marker_any"), q("prev4k", "T_marker_any"), q("floor", "T_marker_any"))
add("I myeloid", "MHC-II score shift, case − control (real-A SD)", "closer to real", q("realU", "mhc_shift_case_minus_ctrl_sdA"), q("scsemi", "mhc_shift_case_minus_ctrl_sdA"), q("prev4k", "mhc_shift_case_minus_ctrl_sdA"), q("floor", "mhc_shift_case_minus_ctrl_sdA"))
S = pd.DataFrame(rows); S.to_csv(f"{R}/scorecard_three_models.csv", index=False); pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 80); print(S.drop(columns="note").to_string())
