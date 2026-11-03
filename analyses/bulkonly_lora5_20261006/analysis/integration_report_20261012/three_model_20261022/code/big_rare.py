"""Rare-population matching (2026-10-21). Usage: python big_rare.py <coh>
Same subset and parent-lineage definitions as integration_report_20261012/code/f_rare_parent.py (slide 23 of the 14 Oct report).
Share = target cells / parent cells (pooled), per group. Real half A uses its author labels; real half B is shown with author labels and
with kNN-transferred labels (same route as generated cells, the fair comparator); generated cells use kNN-transferred labels (k = 30,
Harmony space of each pairwise integration, big_joint.py). Marker check: top-10 marker genes of each subset (half A author labels, subset vs
rest of parent); AUROC of the marker score separating labelled subset cells from the rest of the parent, in every source."""
import sys, os, numpy as np, pandas as pd, anndata as ad
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
coh = sys.argv[1]
IN = "/home/claude/work/s3/big"; OUT = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"
st = lambda m, *p: np.any([np.char.startswith(np.asarray(m, dtype=str), q) for q in p], axis=0)
SEPSIS = {"Mono_c06-CD14-MKI67": lambda m, l: m == "Monocytes", "B_c08-MKI67": lambda m, l: np.isin(m, ["B cells", "Plasma cells"]),
          "T_CD4_c07-MKI67": lambda m, l: m == "CD4+ T cells", "T_CD8_c06-MKI67": lambda m, l: m == "CD8+ T cells", "NK_c04-MKI67": lambda m, l: m == "NK cells",
          "B_c07-MZB1": lambda m, l: np.isin(m, ["B cells", "Plasma cells"]), "pDC": lambda m, l: m == "DC", "cDC": lambda m, l: m == "DC",
          "Macrophage": lambda m, l: np.isin(m, ["Monocytes", "Macrophages"]), "gdT": lambda m, l: np.isin(m, ["CD4+ T cells", "CD8+ T cells", "NKT cells"]) | (l == "gdT"),
          "ILC": lambda m, l: (m == "NK cells") | (l == "ILC"), "HSC": lambda m, l: np.ones(len(m), bool),
          "Mono_c01-CD14-HLA-DRA": lambda m, l: m == "Monocytes", "NKT_c05-IFNG-CCL4hi": lambda m, l: m == "NKT cells"}
PAR_S = {"Mono_c06-CD14-MKI67": "Monocytes", "B_c08-MKI67": "B + plasma", "T_CD4_c07-MKI67": "CD4 T", "T_CD8_c06-MKI67": "CD8 T", "NK_c04-MKI67": "NK",
         "B_c07-MZB1": "B + plasma", "pDC": "DC", "cDC": "DC", "Macrophage": "Mono + macrophage", "gdT": "all T", "ILC": "NK + ILC", "HSC": "all cells",
         "Mono_c01-CD14-HLA-DRA": "Monocytes", "NKT_c05-IFNG-CCL4hi": "NKT"}
COVID = {"CD14mono_dex": lambda m, l: st(m, "CD14mono"), "CD14mono_IFN": lambda m, l: st(m, "CD14mono"), "B_ASC_IgG": lambda m, l: st(m, "B_ASC"),
         "B_ASC_dividing": lambda m, l: st(m, "B_ASC"), "T CD4 cycling": lambda m, l: st(m, "T CD4", "T reg"), "ASDC": lambda m, l: st(m, "DC", "ASDC", "pDC"),
         "DC1": lambda m, l: st(m, "DC", "ASDC", "pDC"), "pDC": lambda m, l: st(m, "DC", "ASDC", "pDC"), "HSPC": lambda m, l: np.ones(len(m), bool),
         "ILC": lambda m, l: st(m, "NK", "ILC"), "T g/d": lambda m, l: st(m, "T ", "MAIT", "NKT"), "MAIT": lambda m, l: st(m, "T ", "MAIT", "NKT"),
         "NK CD56 bright": lambda m, l: st(m, "NK")}
PAR_C = {"CD14mono_dex": "CD14 mono", "CD14mono_IFN": "CD14 mono", "B_ASC_IgG": "plasmablasts", "B_ASC_dividing": "plasmablasts", "T CD4 cycling": "CD4 T",
         "ASDC": "DC", "DC1": "DC", "pDC": "DC", "HSPC": "all cells", "ILC": "NK + ILC", "T g/d": "all T", "MAIT": "all T", "NK CD56 bright": "NK"}
DEF, PAR = (SEPSIS, PAR_S) if coh == "sepsis" else (COVID, PAR_C)
CASE, CTRL = ("adult_sepsis", "adult_HC") if coh == "sepsis" else ("covid_positive", "covid_negative")
R = ad.read_h5ad(f"{IN}/{coh}_ref.h5ad"); genes = np.array(R.var_names)
tr = {s: pd.read_csv(f"/mnt/user-data/outputs/ilisi_floor_20261021/results/big/{coh}_{s}_transfer.csv.gz", index_col=0, keep_default_na=False) for s in ("realU", "floor", "prev4k")}
X = {"ref": L.lognorm(R.X.tocsr())}
for s in tr: X[s] = L.lognorm(ad.read_h5ad(f"{IN}/{coh}_{s}.h5ad").X.tocsr())
V = {"real A (author)": ("ref", R.obs.major.astype(str).to_numpy(), R.obs.fine.astype(str).to_numpy(), R.obs.group.astype(str).to_numpy()),
     "real B (author)": ("realU", tr["realU"].major.to_numpy(), tr["realU"].fine.to_numpy(), tr["realU"].group.to_numpy()),
     "real B (transferred)": ("realU", tr["realU"].major_tr.to_numpy(), tr["realU"].fine_tr.to_numpy(), tr["realU"].group.to_numpy()),
     "floor (transferred)": ("floor", tr["floor"].major_tr.to_numpy(), tr["floor"].fine_tr.to_numpy(), tr["floor"].group.to_numpy()),
     "prev K4000 (transferred)": ("prev4k", tr["prev4k"].major_tr.to_numpy(), tr["prev4k"].fine_tr.to_numpy(), tr["prev4k"].group.to_numpy())}
# depth-matched real B (thinned to floor depth per six-type) for the marker-called shares
import scipy.sparse as sp
rng = np.random.default_rng(20261024); XU = ad.read_h5ad(f"{IN}/{coh}_realU.h5ad"); XF = ad.read_h5ad(f"{IN}/{coh}_floor.h5ad")
Xb = XU.X.tocsr().copy(); tot = np.asarray(Xb.sum(1)).ravel(); gt = np.asarray(XF.X.sum(1)).ravel(); labF = XF.obs.six_type.astype(str).to_numpy()
pool = {c: gt[labF == c] for c in np.unique(labF)}; tgt = np.array([rng.choice(pool.get(c, gt)) for c in XU.obs.six_type.astype(str)], float)
kp = np.minimum(1.0, tgt / np.maximum(tot, 1)); Xb.data = rng.binomial(Xb.data.astype(np.int64), np.repeat(kp, np.diff(Xb.indptr))).astype(np.float32); Xb.eliminate_zeros()
X["realT1f"] = L.lognorm(Xb); del XU, XF
# marker-called (label-free inside the parent): parent from labels (author for real A / real B, transferred for generated); subset = marker score above the
# real-A threshold that maximises F1 against the author label
MC = {"real A (marker)": ("ref", V["real A (author)"][1], V["real A (author)"][2]), "real B (marker)": ("realU", V["real B (author)"][1], V["real B (author)"][2]),
      "real B depth-matched (marker)": ("realT1f", V["real B (author)"][1], V["real B (author)"][2]),
      "floor (marker)": ("floor", V["floor (transferred)"][1], V["floor (transferred)"][2]), "prev K4000 (marker)": ("prev4k", V["prev K4000 (transferred)"][1], V["prev K4000 (transferred)"][2])}
crows = []
rows, mrows = [], []
ma, fa = V["real A (author)"][1], V["real A (author)"][2]
for sub, fn in DEF.items():
    par = fn(ma, fa); tg = par & (fa == sub)
    if tg.sum() < 5: continue
    A = X["ref"]; d = np.asarray(A[tg].mean(0)).ravel() - np.asarray(A[par & ~tg].mean(0)).ravel(); detg = np.asarray((A[tg] > 0).mean(0)).ravel()
    d[detg < 0.1] = -np.inf; top = np.argsort(-d)[:10]
    for lab, (sx, m, l, g) in V.items():
        p = fn(m, l); t = p & (l == sub)
        for grp in ("all", CASE, CTRL):
            gm = np.ones(len(g), bool) if grp == "all" else (g == grp)
            rows.append({"cohort": coh, "subset": sub, "parent": PAR[sub], "source": lab, "group": grp, "n_target": int((t & gm).sum()), "n_parent": int((p & gm).sum()),
                         "pct_of_parent": 100 * (t & gm).sum() / max((p & gm).sum(), 1)})
        sc_ = np.asarray(X[sx][:, top].mean(1)).ravel()
        auc = roc_auc_score(t[p], sc_[p]) if 0 < t[p].sum() < p.sum() else np.nan
        mrows.append({"cohort": coh, "subset": sub, "source": lab, "n_target": int(t.sum()), "marker_auroc": auc, "markers": ",".join(genes[top])})
    sA = np.asarray(X["ref"][:, top].mean(1)).ravel(); yA = tg[par]; xs = sA[par]; best = (0, None)
    for th in np.quantile(xs, np.linspace(0.5, 0.9995, 400)):
        pr = xs > th; tp = (pr & yA).sum(); f1 = 2 * tp / max(pr.sum() + yA.sum(), 1)
        if f1 > best[0]: best = (f1, th)
    for lab, (sx, m, l) in MC.items():
        p = fn(m, l); scv = np.asarray(X[sx][:, top].mean(1)).ravel(); call = p & (scv > best[1])
        crows.append({"cohort": coh, "subset": sub, "parent": PAR[sub], "source": lab, "threshold": best[1], "F1_realA": best[0], "n_called": int(call.sum()), "n_parent": int(p.sum()),
                      "pct_of_parent": 100 * call.sum() / max(p.sum(), 1)})
Cm = pd.DataFrame(crows); Cm.to_csv(f"{OUT}/{coh}_rare_marker_called.csv", index=False)
out = pd.DataFrame(rows); out.to_csv(f"{OUT}/{coh}_rare_parent_share.csv", index=False); M = pd.DataFrame(mrows); M.to_csv(f"{OUT}/{coh}_rare_marker_auroc.csv", index=False)
pd.set_option("display.width", 250)
w = out[out.group == "all"].pivot_table(index="subset", columns="source", values="pct_of_parent").round(2); n = out[out.group == "all"].pivot_table(index="subset", columns="source", values="n_target")
print(w.to_string()); print(n.to_string()); print(Cm.pivot_table(index="subset", columns="source", values="pct_of_parent").round(2).to_string()); print(Cm.drop_duplicates("subset").set_index("subset").F1_realA.round(2).to_string())
print(M.pivot_table(index="subset", columns="source", values="marker_auroc").round(2).to_string())
