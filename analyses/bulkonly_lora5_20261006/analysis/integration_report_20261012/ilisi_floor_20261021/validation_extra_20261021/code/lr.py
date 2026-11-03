"""Ligand-receptor interaction recovery (2026-10-21). Usage: python lr.py <coh>
Resource: LIANA 'consensus' (ligand, receptor; complexes '_' = all subunits, value = min over subunits). Panel genes only.
Per source and six-type (T, NK, B, Monocyte, DC): mean CP10k and detection fraction per gene. Interaction (sender, receiver, L, R) is
'detected' if L is detected in >= 10% of sender cells and R in >= 10% of receiver cells; score = sqrt(mean L x mean R) (expression product).
Compared with real A: number detected, recall of real A's top-200 interactions (detected and in the source's top 500), Spearman of scores
over interactions detected in real A; disease change Delta score (case - control) vs the same in real (A + B), Spearman and sign agreement
over real's 100 largest changes."""
import sys, numpy as np, pandas as pd, liana as li
from scipy.stats import spearmanr
from common import *
coh = sys.argv[1]; D = bundle(coh)
D["realAB"] = ad.concat([D["ref"], D["realU"]], merge="same"); D["realAB"].obs["lineage"] = [LIN.get(x, "Other") for x in D["realAB"].obs.six_type.astype(str)]
genes = np.array(D["ref"].var_names); gpos = pd.Series(np.arange(len(genes)), index=genes); gpos = gpos[~gpos.index.duplicated()]
res = li.resource.select_resource("consensus")
def parts(x): return str(x).split("_")
res = res[[all(p in gpos.index for p in parts(l)) and all(p in gpos.index for p in parts(r)) for l, r in zip(res.ligand, res.receptor)]].drop_duplicates().reset_index(drop=True)
TY = ["T cell", "NK cell", "B cell", "Monocyte", "DC"]
CASE, CTRL = ("adult_sepsis", "adult_HC") if coh == "sepsis" else ("covid_positive", "covid_negative")
def stats_by_type(a, mask=None):
    X = a.X.tocsr(); tot = np.asarray(X.sum(1)).ravel(); C = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); six = a.obs.six_type.astype(str).to_numpy()
    m = np.ones(a.n_obs, bool) if mask is None else mask; out = {}
    for t in TY:
        i = np.where(m & (six == t))[0]
        if len(i) < 30: out[t] = (np.zeros(len(genes)), np.zeros(len(genes))); continue
        out[t] = (np.asarray(C[i].mean(0)).ravel(), np.asarray((X[i] > 0).mean(0)).ravel())
    return out
def lr_table(st):
    rows = []
    for s in TY:
        for r in TY:
            mL, dL = st[s]; mR, dR = st[r]
            for l, rc in zip(res.ligand, res.receptor):
                li_, ri_ = [gpos[p] for p in parts(l)], [gpos[p] for p in parts(rc)]
                ml, dl = mL[li_].min(), dL[li_].min(); mr, dr = mR[ri_].min(), dR[ri_].min()
                rows.append((s, r, l, rc, np.sqrt(ml * mr), dl >= 0.1 and dr >= 0.1))
    T = pd.DataFrame(rows, columns=["sender", "receiver", "ligand", "receptor", "score", "detected"]); T.index = T.sender + "|" + T.receiver + "|" + T.ligand + "|" + T.receptor
    return T
tabs = {s: lr_table(stats_by_type(D[s])) for s in ("ref", "realU", "realT1f", "floor", "prev4k")}
grp = lambda a, g: a.obs.group.astype(str).to_numpy() == g
dis = {}
for s in ("realAB", "floor", "prev4k"):
    a = D[s]; dis[s] = lr_table(stats_by_type(a, grp(a, CASE))).score - lr_table(stats_by_type(a, grp(a, CTRL))).score
R = tabs["ref"]; top200 = R[R.detected].score.nlargest(200).index; rows = []
for s in ("realU", "realT1f", "floor", "prev4k"):
    S = tabs[s]; det = R.detected
    top500 = S[S.detected].score.nlargest(500).index
    row = {"cohort": coh, "source": s, "n_interactions_tested": int(len(R)), "n_detected_realA": int(R.detected.sum()), "n_detected_source": int(S.detected.sum()),
           "detected_jaccard": float((R.detected & S.detected).sum() / max((R.detected | S.detected).sum(), 1)),
           "recall_top200_realA": float(np.isin(top200, top500).mean()), "spearman_score_on_realA_detected": float(spearmanr(R.score[det], S.score[det]).statistic)}
    if s in dis:
        pass
    rows.append(row)
dr = dis["realAB"]; big = dr.abs().nlargest(100).index; drows = []
for s in ("floor", "prev4k"):
    d = dis[s]; ok = R.detected | tabs[s].detected
    drows.append({"cohort": coh, "source": s, "disease_spearman_all_detected": float(spearmanr(dr[ok], d[ok]).statistic),
                  "disease_sign_agree_top100_real": float((np.sign(dr[big]) == np.sign(d[big])).mean()),
                  "disease_spearman_top100_real": float(spearmanr(dr[big], d[big]).statistic)})
O1 = pd.DataFrame(rows); O2 = pd.DataFrame(drows); O1.to_csv(f"{OUT}/{coh}_lr_recovery.csv", index=False); O2.to_csv(f"{OUT}/{coh}_lr_disease.csv", index=False)
pd.concat({s: t for s, t in tabs.items()}).to_csv(f"{OUT}/{coh}_lr_scores.csv.gz")
pd.set_option("display.width", 250); print(O1.round(3).to_string()); print(O2.round(3).to_string())
print("top real disease changes:", dr[big][:10].round(3).to_dict())
