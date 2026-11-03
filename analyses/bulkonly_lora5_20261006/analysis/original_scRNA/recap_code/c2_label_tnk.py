"""T/NK annotation (Stephenson 2025 scheme) from the re-clustered compartment (regress total counts + %mito,
scale, PCA, 40 PCs, Leiden) plus canonical-marker rules applied per cell, as the authors describe
("manually annotated using canonical marker genes through an iterative process of reclustering").
IFN-stimulated states: Leiden (res 3) clusters with mean ISG score > 1.0 (break in the cluster-mean distribution).
Not reproducible from RNA alone: CD45RA/CD45RO protein (authors used CITE-seq) for T CD8 EMRA and T CD8 EM."""
import sys, numpy as np, pandas as pd, anndata as ad
from sklearn.mixture import GaussianMixture
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/covid"
a = ad.read_h5ad(f"{W}/tnk.h5ad"); st = store.Store(f"{W}/store"); o = a.obs
cells = o.cell_index.to_numpy(); g = st.var.gene.to_numpy()
genes = ["CTLA4", "IL2RA", "TYMS", "CD3E", "CD3D", "CD4", "CD8A", "CD8B", "CCR7", "SELL", "LEF1", "TCF7", "FOXP3", "GZMH", "GZMB", "GNLY", "PRF1",
         "NKG7", "KLRB1", "SLC4A10", "TRAV1-2", "TRDC", "NCAM1", "FCGR3A", "XCL1", "KIT", "IL7R", "MKI67", "TOP2A",
         "ISG15", "IFI6", "IFI44L", "MX1", "IFIT3", "RSAD2", "IFI27", "KLRF1"]
gi = [int(np.where(g == x)[0][0]) for x in genes]
E = np.vstack([store.lognorm(X, st.ncount[c])[:, gi].toarray() for c, X in st.chunks(cells)])
E = pd.DataFrame(E, columns=genes)
cl = o.leiden_2.astype(str).to_numpy()
NK = {"17", "18", "19", "20", "21", "30", "31", "36"}; NKb = {"22"}
CTL = {"9", "12", "13", "14", "15", "23", "24", "25", "28", "29", "32", "33", "34", "35", "38", "40", "26"}
NAIVE4 = {"5", "8", "10", "39"}; HELPER4 = {"3", "4", "6", "7", "11"}
isg = E[["ISG15", "IFI6", "IFI44L", "MX1", "IFIT3", "RSAD2"]].mean(1).to_numpy()
cyc = ((E.MKI67 > 0) & ((E.TOP2A > 0) | (E.TYMS > 0))).to_numpy()
cd8 = (E.CD8A + E.CD8B).to_numpy(); cd4 = E.CD4.to_numpy()
lin = np.empty(len(o), dtype=object)
for i, c in enumerate(cl):
    if c == "16": lin[i] = "MAIT"
    elif c == "37": lin[i] = "T g/d"
    elif c == "27": lin[i] = "NKT"
    elif c in NKb: lin[i] = "NK CD56 bright"
    elif c in NK: lin[i] = "NK"
    elif c == "0": lin[i] = "T reg" if (E.FOXP3.iat[i] > 0 or E.CTLA4.iat[i] > 0 or E.IL2RA.iat[i] > 0) else "T CD4 Helper"
    elif c == "1": lin[i] = "T CD8 Naive" if cd8[i] >= cd4[i] else "T CD4 Naive"
    elif c == "2": lin[i] = "T CD8 CM" if cd8[i] >= cd4[i] else "T CD4 Helper"
    elif c in CTL: lin[i] = "T CD8 CTL" if cd8[i] > 0 else "T CD4 CTL"
    elif c in NAIVE4: lin[i] = "T CD8 Naive" if cd8[i] > cd4[i] + 0.5 else "T CD4 Naive"
    elif c in HELPER4: lin[i] = "T reg" if E.FOXP3.iat[i] > 0 else ("T CD8 Naive" if cd8[i] > cd4[i] + 0.5 else "T CD4 Helper")
    else: lin[i] = "T CD4 Helper"
lin = pd.Series(lin)
# rare populations by marker rules
ilc = (E.CD3E == 0) & (E.CD3D == 0) & ((E.KIT > 0) | ((E.IL7R > 0) & (E.KLRB1 > 0) & (E.NKG7 < 1))) & (E.GNLY < 1.5)
lin[ilc.to_numpy()] = "ILC"
cyc4 = cyc & lin.str.startswith("T CD4").to_numpy() | cyc & (lin == "T reg").to_numpy()
lin[cyc4] = "T CD4 cycling"
detailed = lin.copy()
# IFN-stimulated states per lineage (authors' cell_type_IFNstim)
ifn_lineages = {"NK": "NK IFN stim", "T CD4 Helper": "T CD4 Helper IFN stim", "T CD4 Naive": "T CD4 Naive IFN stim",
                "T CD8 CTL": "T CD8 CTL IFN stim", "T CD8 Naive": "T CD8 Naive IFN stim", "T reg": "T reg IFN stim"}
ifnstim = lin.copy()
# IFN-stimulated states defined at cluster level (finest Leiden, res 3): clusters whose mean ISG score lies above the
# natural break (1.0) of the cluster-mean distribution; applied to lineages that have an IFN-stimulated state.
cl3 = o.leiden_3.astype(str).to_numpy()
cm = pd.Series(isg).groupby(cl3).mean(); ifn_cl = set(cm[cm > 1.0].index)
hi = np.isin(cl3, list(ifn_cl))
for k, v in ifn_lineages.items():
    m = (lin == k).to_numpy() & hi; ifnstim[m] = v
out = pd.DataFrame({"cell_type_detailed_recap": detailed.to_numpy(), "cell_type_IFNstim_recap": ifnstim.to_numpy(),
                    "recap_subcluster": "tnk_" + o.leiden_2.astype(str).to_numpy(), "isg_score": isg}, index=cells)
out.to_csv(f"{W}/labels_tnk.csv")
for acol, rcol in [("cell_type_detailed", "cell_type_detailed_recap"), ("cell_type_IFNstim", "cell_type_IFNstim_recap")]:
    au = o[acol].astype(str).to_numpy(); rc = out[rcol].to_numpy()
    print(rcol, "agreement:", round(float((au == rc).mean()), 3))
    print({c: (int((au == c).sum()), round(float((rc[au == c] == c).mean()), 3)) for c in sorted(set(au))})
print(out.cell_type_detailed_recap.value_counts().to_dict())
