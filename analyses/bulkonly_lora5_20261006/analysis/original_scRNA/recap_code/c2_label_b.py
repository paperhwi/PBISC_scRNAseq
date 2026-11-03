"""B/ASC annotation following Stephenson 2025 STAR Methods (Annotation - B cell):
ASC vs non-ASC from a two-component Gaussian mixture on a plasma-cell signature score; non-ASC clusters annotated
with TCL1A/IGHD/CD27/ITGAX/FCRL5; isotype from scBCR constant gene (fallback: dominant constant-region mRNA)
so that naive / non-switched memory / IgM ASC carry IgM/IgD and switched memory / IgA / IgG ASC carry IgG/IgA;
ASC with IgD/IgE -> B_ASC_other. Cycling ASC (MKI67/TOP2A) -> B_ASC_dividing.
Integration used Harmony instead of scVI (no GPU / torch in this environment)."""
import sys, h5py, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.mixture import GaussianMixture
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/covid"
a = ad.read_h5ad(f"{W}/bcell.h5ad"); st = store.Store(f"{W}/store"); o = a.obs.copy()
cells = o.cell_index.to_numpy(); g = st.var.gene.to_numpy()
def expr(genes):
    gi = [int(np.where(g == x)[0][0]) for x in genes]
    M = []
    for c, X in st.chunks(cells):
        M.append(store.lognorm(X, st.ncount[c])[:, gi].toarray())
    return pd.DataFrame(np.vstack(M), columns=genes, index=o.index)
E = expr(["JCHAIN", "MZB1", "XBP1", "PRDM1", "TNFRSF17", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4",
          "IGHM", "IGHD", "IGHE", "MKI67", "TOP2A", "MS4A1", "CD19", "TCL1A", "CD27", "ITGAX", "FCRL5"])
asc_score = E[["JCHAIN", "MZB1", "XBP1", "PRDM1", "TNFRSF17"]].mean(1) - E[["MS4A1", "CD19"]].mean(1)
gm = GaussianMixture(2, covariance_type="full", max_iter=1000, random_state=0).fit(asc_score.to_numpy()[:, None])
comp = gm.predict(asc_score.to_numpy()[:, None]); asc = comp == np.argmax(gm.means_.ravel())
# BCR constant gene from the deposited object (scirpy/dandelion fields)
with h5py.File(f"{W}/eskd_covid19.h5ad", "r") as f:
    gg = f["obs/IR_VDJ_1_c_gene_bcr"]; cats = np.array([x.decode() for x in gg["categories"][:]]); codes = gg["codes"][:][cells]
bcr = np.where(codes >= 0, cats[np.maximum(codes, 0)], "")
def iso_from(cg):
    cg = cg.upper()
    if cg.startswith("IGHA"): return "IgA"
    if cg.startswith("IGHG"): return "IgG"
    if cg.startswith("IGHM"): return "IgM"
    if cg.startswith("IGHD"): return "IgD"
    if cg.startswith("IGHE"): return "IgE"
    return ""
iso = np.array([iso_from(x) for x in bcr])
mrna = pd.DataFrame({"IgA": E[["IGHA1", "IGHA2"]].max(1), "IgG": E[["IGHG1", "IGHG2", "IGHG3", "IGHG4"]].max(1),
                     "IgM": E["IGHM"], "IgD": E["IGHD"], "IgE": E["IGHE"]})
iso_m = mrna.idxmax(1).to_numpy(); iso_m[mrna.max(1).to_numpy() == 0] = ""
iso_final = np.where(iso != "", iso, iso_m)
cyc = (E[["MKI67", "TOP2A"]].max(1) > 0.5).to_numpy()
cl = o.leiden_1.astype(str).to_numpy()
naive_cl = {"3", "4", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "19", "20", "21"}
lab = np.empty(len(o), dtype=object)
for i in range(len(o)):
    if asc[i]:
        if cyc[i] or cl[i] == "0": lab[i] = "B_ASC_dividing"
        else: lab[i] = {"IgA": "B_ASC_IgA", "IgG": "B_ASC_IgG", "IgM": "B_ASC_IgM"}.get(iso_final[i], "B_ASC_other")
    else:
        if cl[i] == "6" and (E.ITGAX.iat[i] > 0 or E.FCRL5.iat[i] > 0): lab[i] = "B_CD11c"
        elif cl[i] in naive_cl and iso_final[i] in ("IgM", "IgD", ""): lab[i] = "B_naive"
        elif iso_final[i] in ("IgG", "IgA", "IgE"): lab[i] = "B_switched_mem"
        else: lab[i] = "B_non-switched_mem"
out = pd.DataFrame({"cell_type_detailed_recap": lab, "recap_subcluster": "bcell_" + o.leiden_1.astype(str).to_numpy(),
                    "asc_gmm": asc, "isotype_recap": iso_final, "bcr_c_gene": bcr}, index=cells)
out.to_csv(f"{W}/labels_b.csv")
auth = o.cell_type_detailed.astype(str).to_numpy()
print("agreement:", round(float((auth == lab).mean()), 3))
print({c: (int((auth == c).sum()), round(float((lab[auth == c] == c).mean()), 3)) for c in sorted(set(auth))})
print(pd.Series(lab).value_counts().to_dict())
