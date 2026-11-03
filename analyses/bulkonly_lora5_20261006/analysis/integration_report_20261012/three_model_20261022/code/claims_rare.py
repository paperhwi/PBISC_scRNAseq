"""(1) Paper-claim re-test on the floor read-out (2026-10-21) and (2) transfer-free rare-state estimates. Usage: python claims_rare.py <coh>
(1) Sepsis Ye claims SEP-1/2/3/3-NK/4A/4B: panel-matched author markers (results/rescore/panel_marker_coverage.csv), score = mean log1p CP10k;
    parent = author fine-label prefix (real) or kNN-transferred fine label prefix (generated); threshold = Youden on real A within parent
    (target = exact author cluster); adult sepsis vs adult HC contrast of sample-mean positive fraction and of mean score (in real-A SD units).
    COVID IFN (7 ISG core, frozen definition): per six-type parent, COVID+ vs COVID- contrast; truth label = author 'IFN' state.
(2) Rare subsets (big_rare.py definitions): HistGradientBoosting within parent on [PCA 50 of real-A HVG space + top-50 subset markers],
    trained on real A thinned to the depth of the source (floor or prev, per six-type), 5-fold sample-grouped CV AUC on real A; share =
    mean predicted probability; checked on real B thinned the same way (author labels known); then applied to generated parents."""
import sys, json, numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, roc_curve
from common import *
coh = sys.argv[1]; D = bundle(coh); genes = np.array(D["ref"].var_names); gpos = pd.Series(np.arange(len(genes)), index=genes); gpos = gpos[~gpos.index.duplicated()]
TRN = {s: transferred(coh, s) for s in ("realU", "floor", "prev4k", "scsemi")}
lab = {"ref": (D["ref"].obs.fine.astype(str).to_numpy(), D["ref"].obs.major.astype(str).to_numpy()),
       "realU": (D["realU"].obs.fine.astype(str).to_numpy(), D["realU"].obs.major.astype(str).to_numpy())}
lab["realT1f"] = lab["realU"]; lab["realT1p"] = lab["realU"]; lab["realT1s"] = lab["realU"]
for s in ("floor", "prev4k", "scsemi"): lab[s] = (TRN[s].fine_tr.to_numpy(), TRN[s].major_tr.to_numpy())
lab["realU_tr"] = (TRN["realU"].fine_tr.to_numpy(), TRN["realU"].major_tr.to_numpy())
grp = {s: a.obs.group.astype(str).to_numpy() for s, a in D.items()}; smp = {s: a.obs["sample"].astype(str).to_numpy() for s, a in D.items()}
Y = {s: L.lognorm(a.X.tocsr()) for s, a in D.items()}
CASE, CTRL = ("adult_sepsis", "adult_HC") if coh == "sepsis" else ("covid_positive", "covid_negative")
def youden(s, y):
    f, t, th = roc_curve(y, s); return th[np.argmax(t - f)]
def contrast(vals, g, sm):
    d = pd.DataFrame({"v": vals, "g": g, "s": sm}); m = d.groupby(["s", "g"]).v.mean().reset_index()
    a, b = m[m.g == CASE].v, m[m.g == CTRL].v; return float(a.mean() - b.mean()), int(len(a)), int(len(b))
crow = []
if coh == "sepsis":
    cov = pd.read_csv("/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/results/rescore/panel_marker_coverage.csv")
    defs = pd.read_csv("/mnt/user-data/uploads/PBISC/paper_claim_analysis/fixed_paper_claim_definitions.csv", dtype=str).fillna("").set_index("id")
    claims = [(c, defs.loc[c, "lineage_prefix"], defs.loc[c, "target"], [g for g in cov.set_index("claim").loc[c, "panel_markers_list"].split(";") if g in gpos.index]) for c in ("SEP-1", "SEP-2", "SEP-3", "SEP-3-NK", "SEP-4A", "SEP-4B")]
else:
    isg = [g for g in "ISG15;IFIT1;IFIT3;MX1;OAS1;OAS2;IFI6".split(";") if g in gpos.index]
    claims = [(f"ESKD-IFN {t}", t, "IFN", isg) for t in ("T cell", "NK cell", "B cell", "Monocyte", "DC")]
for cid, pref, tgt, mk in claims:
    w = gpos[mk].to_numpy(); sc = {s: np.asarray(Y[s][:, w].mean(1)).ravel() for s in Y}
    def parent(s):
        if coh == "sepsis": f = lab[s][0]; return np.char.startswith(f.astype(str), pref)
        return D[s].obs.six_type.astype(str).to_numpy() == pref
    def target(s):
        f = lab[s][0]; return (f == tgt) if coh == "sepsis" else np.char.find(f.astype(str), "IFN") >= 0
    pa = parent("ref"); th = youden(sc["ref"][pa], target("ref")[pa]); auc = roc_auc_score(target("ref")[pa], sc["ref"][pa]); sdA = sc["ref"][pa].std()
    for s in ("ref", "realU", "realT1f", "realT1s", "floor", "prev4k", "scsemi"):
        p = parent(s); dpos = contrast((sc[s][p] >= th).astype(float), grp[s][p], smp[s][p]); dsc = contrast(sc[s][p] / sdA, grp[s][p], smp[s][p])
        r = {"cohort": coh, "claim": cid, "n_markers": len(mk), "auc_realA": float(auc), "source": s, "n_parent": int(p.sum()), "delta_pos_frac": dpos[0], "delta_score_sdA": dsc[0],
             "n_case_samples": dpos[1], "n_ctrl_samples": dpos[2]}
        if s in ("ref", "realU"): r["delta_label_frac"] = contrast(target(s)[p].astype(float), grp[s][p], smp[s][p])[0]
        crow.append(r)
CL = pd.DataFrame(crow); CL.to_csv(f"{OUT}/{coh}_claims.csv", index=False)
pd.set_option("display.width", 250)
print(CL.pivot_table(index="claim", columns="source", values=["delta_pos_frac", "delta_score_sdA"]).round(3).to_string())
print(CL[CL.source.isin(["ref", "realU"])][["claim", "source", "auc_realA", "delta_label_frac", "n_case_samples", "n_ctrl_samples"]].round(3).to_string())
# ---------------- (2) transfer-free rare
sys.path.insert(0, "/home/claude/work/big/code")
exec(open("/home/claude/work/s3/code/big_rare.py").read().split("DEF, PAR =")[0].split("coh = sys.argv[1]")[1])   # SEPSIS/COVID dicts and PAR maps
DEF, PAR = (SEPSIS, PAR_S) if coh == "sepsis" else (COVID, PAR_C)
OUT = "/mnt/user-data/outputs/s3/results"   # the exec above redefines OUT
hv = L.hvg_seurat_v3([D["ref"].X.tocsr()]); A0 = Y["ref"][:, hv].toarray().astype(np.float32); mu, sd = A0.mean(0), A0.std(0) + 1e-6
pca = PCA(50, random_state=0).fit(np.clip((A0 - mu) / sd, -10, 10)); del A0
emb = lambda M: pca.transform(np.clip((M[:, hv].toarray().astype(np.float32) - mu) / sd, -10, 10))
TA = {"floor": thin_to(D["ref"], D["floor"], 7), "prev4k": thin_to(D["ref"], D["prev4k"], 8), "scsemi": thin_to(D["ref"], D["scsemi"], 9)}
YT = {k: L.lognorm(v.X.tocsr()) for k, v in TA.items()}; ET = {k: emb(v) for k, v in YT.items()}
E = {s: emb(Y[s]) for s in ("realT1f", "realT1p", "realT1s", "floor", "prev4k", "scsemi")}
rrows = []
fA, mA = lab["ref"]
for sub, fn in DEF.items():
    parA = fn(mA, fA); tg = parA & (fA == sub)
    if tg.sum() < 15: continue
    YA = Y["ref"]; d = np.asarray(YA[tg].mean(0)).ravel() - np.asarray(YA[parA & ~tg].mean(0)).ravel(); top = np.argsort(-d)[:50]
    for depth, (tsrc, bsrc) in {"floor": ("floor", "realT1f"), "prev4k": ("prev4k", "realT1p"), "scsemi": ("scsemi", "realT1s")}.items():
        Xtr = np.hstack([ET[depth][parA], YT[depth][parA][:, top].toarray()]); ytr = (fA[parA] == sub).astype(int); gtr = smp["ref"][parA]
        p = np.zeros(len(ytr))
        for a, b in GroupKFold(5).split(Xtr, ytr, gtr):
            if ytr[a].sum() < 5: continue
            p[b] = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.1, random_state=0).fit(Xtr[a], ytr[a]).predict_proba(Xtr[b])[:, 1]
        auc = roc_auc_score(ytr, p) if 0 < ytr.sum() < len(ytr) else np.nan
        clf = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.1, random_state=0).fit(Xtr, ytr)
        out = {"cohort": coh, "subset": sub, "parent": PAR[sub], "model_depth": depth, "cv_auc_realA": float(auc), "realA_truth_pct": float(100 * ytr.mean()),
               "realA_cv_estimate_pct": float(100 * p.mean())}
        fB, mB = lab[bsrc]; pB = fn(mB, fB)
        out["realB_truth_pct"] = float(100 * (fB[pB] == sub).mean()); out["realB_estimate_pct"] = float(100 * clf.predict_proba(np.hstack([E[bsrc][pB], Y[bsrc][pB][:, top].toarray()]))[:, 1].mean())
        fG, mG = lab[tsrc]; pG = fn(mG, fG)
        out["gen_estimate_pct"] = float(100 * clf.predict_proba(np.hstack([E[tsrc][pG], Y[tsrc][pG][:, top].toarray()]))[:, 1].mean()); out["gen_n_parent"] = int(pG.sum())
        rrows.append(out)
    print(sub, [(r["model_depth"], round(r["cv_auc_realA"], 2), round(r["realB_truth_pct"], 2), round(r["realB_estimate_pct"], 2), round(r["gen_estimate_pct"], 2)) for r in rrows[-3:]], flush=True)
RR = pd.DataFrame(rrows); RR.to_csv(f"{OUT}/{coh}_rare_transfer_free.csv", index=False)
