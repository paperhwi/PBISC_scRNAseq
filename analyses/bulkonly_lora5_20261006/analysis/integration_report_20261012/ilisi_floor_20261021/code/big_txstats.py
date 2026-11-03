"""Transcriptome characteristics of the larger bundles vs real (2026-10-21). Usage: python big_txstats.py <coh>
Sources: ref (real half A, full depth), realU (real half B), realT1f / realT1p (half B binomially thinned, per six-type, to the
depth of floor / prev K4000 cells), floor, prev4k. Outputs (results/big/<coh>_tx_*.csv):
  gene: per-gene mean CP10k, detection fraction, mean count, variance; rank by real half A mean
  cell: per-cell UMI and detected genes
  pb:   per-six-type pseudobulk log1p mean CP10k correlation with real half A (all genes / 2,000 ref HVG / genes ranked 2,001+)
  dot:  marker panel fraction detected and mean log CP10k per six-type"""
import sys, os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ilisi_lib as L
coh = sys.argv[1]
IN = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/big"; OUT = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"
D = {s: ad.read_h5ad(f"{IN}/{coh}_{s}.h5ad") for s in ("ref", "realU", "floor", "prev4k")}
genes = np.array(D["ref"].var_names); rng = np.random.default_rng(20261023)
def thin(B, G):
    X = B.X.tocsr().copy(); tot = np.asarray(X.sum(1)).ravel(); gt = np.asarray(G.X.sum(1)).ravel(); lab = G.obs.six_type.astype(str).to_numpy()
    pool = {c: gt[lab == c] for c in np.unique(lab)}; six = B.obs.six_type.astype(str).to_numpy()
    tgt = np.array([rng.choice(pool.get(c, gt)) for c in six], float); keep = np.minimum(1.0, tgt / np.maximum(tot, 1))
    X.data = rng.binomial(X.data.astype(np.int64), np.repeat(keep, np.diff(X.indptr))).astype(np.float32); X.eliminate_zeros()
    return ad.AnnData(X, obs=B.obs.copy(), var=B.var.copy())
D["realT1f"] = thin(D["realU"], D["floor"]); D["realT1p"] = thin(D["realU"], D["prev4k"])
ORDER = ["ref", "realU", "realT1f", "realT1p", "floor", "prev4k"]
G, C, P, DOT = [], [], [], []
hv = L.hvg_seurat_v3([D["ref"].X.tocsr()]); refmean = None
MK = {"T cell": ["CD3E", "CD3D", "IL7R", "CD8A"], "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"], "B cell": ["MS4A1", "CD79A", "CD79B", "TCL1A"],
      "Monocyte": ["CD14", "LYZ", "S100A8", "FCGR3A"], "DC": ["FCER1A", "CD1C", "CLEC10A", "LILRA4"], "Other": ["PPBP", "PF4", "HBB", "MKI67"]}
mk = [g for v in MK.values() for g in v if g in set(genes)]
pb = {}
for s in ORDER:
    X = D[s].X.tocsr().astype(np.float32); Y = L.lognorm(X); cp = (sp.diags(1e4 / np.maximum(np.asarray(X.sum(1)).ravel(), 1)) @ X).tocsr()
    m = np.asarray(cp.mean(0)).ravel(); det = np.asarray((X > 0).mean(0)).ravel(); mc = np.asarray(X.mean(0)).ravel(); var = np.asarray(X.multiply(X).mean(0)).ravel() - mc ** 2
    G.append(pd.DataFrame({"gene": genes, "source": s, "mean_cp10k": m, "det": det, "mean_count": mc, "var_count": var}))
    C.append(pd.DataFrame({"source": s, "umi": np.asarray(X.sum(1)).ravel(), "n_genes": np.diff(X.indptr), "six_type": D[s].obs.six_type.astype(str).to_numpy()}))
    six = D[s].obs.six_type.astype(str).to_numpy()
    for t in ["T cell", "NK cell", "B cell", "Monocyte", "DC"]:
        mm = six == t
        if mm.sum() < 50: continue
        pb[(s, t)] = np.log1p(np.asarray(cp[mm].mean(0)).ravel())
        gi = [list(genes).index(g) for g in mk]; sub = Y[mm][:, gi]
        DOT.append(pd.DataFrame({"source": s, "six_type": t, "gene": mk, "frac": np.asarray((sub > 0).mean(0)).ravel(), "mean_log": np.asarray(sub.mean(0)).ravel(), "n_cells": int(mm.sum())}))
G = pd.concat(G); rk = G[G.source == "ref"].set_index("gene").mean_cp10k.rank(ascending=False, method="first"); G["real_rank"] = G.gene.map(rk).astype(int)
gset = {"all genes": np.ones(len(genes), bool), "2,000 ref HVG": np.isin(np.arange(len(genes)), hv), "genes ranked 2,001+": (pd.Series(genes).map(rk) > 2000).to_numpy()}
for (s, t), v in pb.items():
    if s == "ref": continue
    for k, msk in gset.items():
        r0 = pb[("ref", t)]; P.append({"source": s, "six_type": t, "gene_set": k, "pearson": float(np.corrcoef(r0[msk], v[msk])[0, 1]), "n_genes": int(msk.sum())})
G.to_csv(f"{OUT}/{coh}_tx_gene.csv.gz", index=False); pd.concat(C).to_csv(f"{OUT}/{coh}_tx_cell.csv.gz", index=False)
pd.DataFrame(P).to_csv(f"{OUT}/{coh}_tx_pb.csv", index=False); pd.concat(DOT).to_csv(f"{OUT}/{coh}_tx_dot.csv", index=False)
Cc = pd.concat(C); print(coh, Cc.groupby("source")[["umi", "n_genes"]].median().loc[ORDER].round(0).to_string())
print(pd.DataFrame(P).pivot_table(index=["gene_set", "source"], columns="six_type", values="pearson").round(3).to_string())
