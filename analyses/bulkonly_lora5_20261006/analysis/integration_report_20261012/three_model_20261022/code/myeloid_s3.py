"""(5) Sepsis myeloid gene-level diagnosis (2026-10-22). Usage: python myeloid_s3.py [sepsis|covid]
Cells: six-type 'Monocyte' (real: author labels; generated: frozen scoring label). Depth control: real B thinned per six-type to each model's depth
(realT1f / realT1p / realT1s); null = real half A thinned to the same depth vs real B thinned (between-patient + thinning noise, no model).
(a) Per gene mean log1p CP10k and detection fraction; model residual = model - depth-matched real B; 'model-specific' genes = |residual| above the
    99th percentile of |null| over genes; categories of those genes (T/NK, ribosomal, mitochondrial, MHC-II, ISG, S100/alarmin, other).
(b) T/NK co-detection in monocytes: share of cells with >= 1 UMI of a T marker (CD3E, CD3D, CD2, IL7R, TRAC) or NK marker (NKG7, GNLY, GZMB, PRF1, KLRD1).
(c) Identity programme: classical-monocyte score (CD14, LYZ, S100A8, S100A9, FCN1, VCAN, CSF1R, CST3, MNDA) and MHC-II score (HLA-DRA, HLA-DRB1, HLA-DPA1,
    HLA-DPB1, CD74) per cell; share below real-B 5th percentile; sepsis vs HC shift of the MHC-II score (sample means, in real-A SD) - known biology:
    monocyte HLA-DR falls in sepsis.
(d) Monocyte M6 genes (m6_three_models_gene_tables): real DEGs whose sign each model misses, by category."""
import sys, numpy as np, pandas as pd
from common import *
coh = sys.argv[1] if len(sys.argv) > 1 else "sepsis"; D = bundle(coh)
genes = np.array(D["ref"].var_names); gi = pd.Series(np.arange(len(genes)), index=genes); gi = gi[~gi.index.duplicated()]
TA = {"floor": thin_to(D["ref"], D["floor"], 7), "prev4k": thin_to(D["ref"], D["prev4k"], 8), "scsemi": thin_to(D["ref"], D["scsemi"], 9)}
MB = {"floor": "realT1f", "prev4k": "realT1p", "scsemi": "realT1s"}
def mono(a): return np.where(a.obs.six_type.astype(str).to_numpy() == "Monocyte")[0]
def stats(a):
    i = mono(a); X = a.X.tocsr()[i]; Y = L.lognorm(X); return np.asarray(Y.mean(0)).ravel(), np.asarray((X > 0).mean(0)).ravel(), X, Y, a.obs.iloc[i]
ST = {s: stats(a) for s, a in D.items()}; STA = {k: stats(v) for k, v in TA.items()}
def cat(g):
    if g.startswith(("RPL", "RPS")): return "ribosomal"
    if g.startswith("MT-"): return "mitochondrial"
    if g.startswith("HLA-D") or g == "CD74": return "MHC-II"
    if g in {"ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI44L", "MX1", "OAS1", "OAS2", "OAS3", "RSAD2", "IFITM3", "LY6E", "XAF1"}: return "ISG"
    if g.startswith("S100") or g in {"LCN2", "CAMP", "LTF", "MMP8", "MMP9", "ARG1", "RETN", "CD177"}: return "alarmin/granule"
    if g in TNK or g.startswith(("TRA", "TRB", "KLR", "GZM", "CD3", "CD8")): return "T/NK"
    if g.startswith(("IGH", "IGK", "IGL")) or g in {"JCHAIN", "MZB1"}: return "Ig/plasma"
    if g in {"HBB", "HBA1", "HBA2", "PPBP", "PF4", "GP9", "ITGA2B"}: return "platelet/erythroid"
    return "other"
TM = ["CD3E", "CD3D", "CD2", "IL7R", "TRAC"]; NM = ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1"]; TNK = set(TM + NM + ["CCL5", "CD247", "LCK", "CST7", "GZMA", "GZMH", "CTSW"])
rows, grows = [], []
for m, b in MB.items():
    mm, dm, *_ = ST[m]; mb, db, *_ = ST[b]; ma, da, *_ = STA[m]
    res = mm - mb; null = ma - mb; thr = np.quantile(np.abs(null), 0.99)
    sel = np.abs(res) > thr; order = np.argsort(-np.abs(res))
    cc = pd.Series([cat(g) for g in genes[sel]]).value_counts()
    rows.append({"cohort": coh, "model": m, "null_abs_q99": float(thr), "n_model_specific": int(sel.sum()), "n_up": int((res > thr).sum()), "n_down": int((res < -thr).sum()),
                 "null_n_over_thr": int((np.abs(null) > thr).sum()), "r_mean_model_vs_realB": float(np.corrcoef(mm, mb)[0, 1]), "r_mean_null": float(np.corrcoef(ma, mb)[0, 1]),
                 **{f"cat_{k}": int(v) for k, v in cc.items()}})
    for k in order[:60]:
        grows.append({"cohort": coh, "model": m, "gene": genes[k], "category": cat(genes[k]), "model_mean_log": mm[k], "realB_dm_mean_log": mb[k], "residual": res[k],
                      "null_residual": null[k], "model_det": dm[k], "realB_dm_det": db[k]})
R1 = pd.DataFrame(rows).fillna(0); R1.to_csv(f"{OUT}/{coh}_myeloid_model_specific_summary.csv", index=False); pd.DataFrame(grows).to_csv(f"{OUT}/{coh}_myeloid_model_specific_genes.csv", index=False)
# (b)(c)
mono_sc = [g for g in ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN", "CSF1R", "CST3", "MNDA"] if g in gi.index]
mhc = [g for g in ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74"] if g in gi.index]
CASE, CTRL = ("adult_sepsis", "adult_HC") if coh == "sepsis" else ("covid_positive", "covid_negative")
_, _, XA, YA, _ = ST["ref"]; sdA = {n: np.asarray(YA[:, gi[g].to_numpy()].mean(1)).ravel().std() for n, g in (("mono", mono_sc), ("mhc", mhc))}
_, _, XB, YB, _ = ST["realU"]; p5 = {n: np.quantile(np.asarray(YB[:, gi[g].to_numpy()].mean(1)).ravel(), 0.05) for n, g in (("mono", mono_sc), ("mhc", mhc))}
crow = []
for s in ("ref", "realU", "realT1f", "floor", "realT1p", "prev4k", "realT1s", "scsemi"):
    _, _, X, Y, ob = ST[s]
    t_any = np.asarray((X[:, gi[[g for g in TM if g in gi.index]].to_numpy()] > 0).sum(1)).ravel() > 0
    n_any = np.asarray((X[:, gi[[g for g in NM if g in gi.index]].to_numpy()] > 0).sum(1)).ravel() > 0
    sc_ = {n: np.asarray(Y[:, gi[g].to_numpy()].mean(1)).ravel() for n, g in (("mono", mono_sc), ("mhc", mhc))}
    d = pd.DataFrame({"s": ob["sample"].astype(str).to_numpy(), "g": ob.group.astype(str).to_numpy(), "mhc": sc_["mhc"] / sdA["mhc"], "mono": sc_["mono"] / sdA["mono"]})
    sm = d.groupby(["s", "g"])[["mhc", "mono"]].mean().reset_index()
    tot = np.asarray(X.sum(1)).ravel(); ngen = np.asarray((X > 0).sum(1)).ravel()
    crow.append({"cohort": coh, "source": s, "n_mono": int(X.shape[0]), "median_UMI": float(np.median(tot)), "median_genes": float(np.median(ngen)),
                 "T_marker_any": float(t_any.mean()), "NK_marker_any": float(n_any.mean()), "T_or_NK_any": float((t_any | n_any).mean()),
                 "mono_score_median_sdA": float(np.median(sc_["mono"] / sdA["mono"])), "mono_below_realB_p5": float((sc_["mono"] < p5["mono"]).mean()),
                 "mhc_score_median_sdA": float(np.median(sc_["mhc"] / sdA["mhc"])), "mhc_below_realB_p5": float((sc_["mhc"] < p5["mhc"]).mean()),
                 "mhc_shift_case_minus_ctrl_sdA": float(sm[sm.g == CASE].mhc.mean() - sm[sm.g == CTRL].mhc.mean()),
                 "mono_shift_case_minus_ctrl_sdA": float(sm[sm.g == CASE].mono.mean() - sm[sm.g == CTRL].mono.mean())})
R2 = pd.DataFrame(crow); R2.to_csv(f"{OUT}/{coh}_myeloid_identity.csv", index=False)
# (d) M6 monocyte missed-sign genes
G = pd.read_csv(f"{OUT}/m6_three_models_gene_tables.csv.gz"); mrow = []
for src in ("floor", "prev", "scsemi"):
    d = G[G.key == f"{coh}|case vs control|Monocyte|{src}"]; sig = d[(d.real_q < 0.05) & (d.real_mean > 1)].copy()
    sig["category"] = [cat(g) for g in sig.gene]; sig["miss"] = np.sign(sig.real_logFC) != np.sign(sig.gen_logFC)
    t = sig.groupby("category").agg(n_DEG=("miss", "size"), miss_rate=("miss", "mean")).reset_index(); t["source"] = src; t["cohort"] = coh; mrow.append(t)
    sig.assign(source=src).sort_values("real_logFC", key=np.abs, ascending=False).head(0)
R3 = pd.concat(mrow); R3.to_csv(f"{OUT}/{coh}_myeloid_m6_missed_by_category.csv", index=False)
pd.set_option("display.width", 250); print(R1.round(3).to_string()); print(R2.round(3).to_string()); print(R3.pivot_table(index="category", columns="source", values=["n_DEG", "miss_rate"]).round(2).to_string())
print(pd.DataFrame(grows).groupby("model").head(15)[["model", "gene", "category", "residual", "null_residual", "model_det", "realB_dm_det"]].round(3).to_string())
