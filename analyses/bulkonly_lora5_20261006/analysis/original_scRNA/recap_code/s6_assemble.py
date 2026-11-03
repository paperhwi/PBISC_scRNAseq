"""Assemble per-cell annotation for the Ye 2026 sepsis recapitulation and compare composition with the paper (Fig. 1d,e)."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/sepsis"; OUTD = "/home/claude/work/final/sepsis_ann"
import os; os.makedirs(OUTD, exist_ok=True)
st = store.Store(f"{W}/store"); g = st.var.gene_symbol.astype(str).to_numpy()
obs = st.obs.copy(); N = len(obs)
major = pd.read_csv(f"{W}/run_major/major_lineage.csv")
assert (major.cell.to_numpy() == np.arange(N)).all()
obs["major_lineage"] = major["major"].to_numpy(); obs["major_leiden_1.5"] = major["leiden_1.5"].astype(str).to_numpy()
L = pd.concat([pd.read_csv(f"{W}/{f}") for f in ["labels_b_my.csv", "labels_cd4.csv", "labels_cyto.csv", "labels_other.csv"]])
assert L.cell.is_unique and len(L) == N, (L.cell.duplicated().sum(), len(L), N)
L = L.set_index("cell").sort_index()
obs["celltype_ye2026"] = L.celltype_ye2026.to_numpy(); obs["recap_subcluster"] = L.subcluster.to_numpy()
obs["qc_flag"] = L.qc_flag.fillna("").to_numpy()
# ILC: CD3- KIT+ (or IL7R+ KLRB1+ cytotoxic-low) cells within T/NK
tnk = np.where(obs.major_lineage.eq("T/NK").to_numpy())[0]
genes = ["CD3D", "CD3E", "KIT", "IL7R", "KLRB1", "GNLY", "NKG7"]
gi = [int(np.where(g == x)[0][0]) for x in genes]
E = pd.DataFrame(np.vstack([store.lognorm(X, st.ncount[c])[:, gi].toarray() for c, X in st.chunks(tnk)]), columns=genes)
ilc = ((E.CD3D == 0) & (E.CD3E == 0) & ((E.KIT > 0) | ((E.IL7R > 0) & (E.KLRB1 > 0) & (E.NKG7 < 1))) & (E.GNLY < 1.5)).to_numpy()
obs.iloc[tnk[ilc], obs.columns.get_loc("celltype_ye2026")] = "ILC"
# paper's coarse categories (Fig. 1c-e)
def coarse(x):
    if x.startswith("T_CD4"): return "CD4+ T cells"
    if x.startswith("T_CD8"): return "CD8+ T cells"
    if x.startswith("NKT"): return "NKT cells"
    if x.startswith("NK_"): return "NK cells"
    if x in ("B_c07-MZB1", "B_c08-MKI67"): return "Plasma cells"
    if x.startswith("B_"): return "B cells"
    if x.startswith("Mono_"): return "Monocytes"
    if x == "Macrophage": return "Macrophages"
    if x in ("cDC", "pDC"): return "DC"
    if x in ("gdT", "HSC", "Megakaryocyte", "ILC"): return "Other cells"
    return "Not in paper scheme"
obs["celltype_major_ye2026"] = obs.celltype_ye2026.map(coarse)
obs["in_paper_scheme"] = obs.celltype_major_ye2026.ne("Not in paper scheme") & ~obs.qc_flag.str.contains("doublet|complex|low_quality")
# composition vs paper (fraction of cells in paper-scheme categories, by cohort)
paper = {"adult": {"CD4+ T cells": 33.07, "CD8+ T cells": 19.86, "B cells": 13.49, "Plasma cells": 0.45, "NK cells": 5.47,
                   "NKT cells": 6.84, "Other cells": 2.77, "Macrophages": 0.77, "Monocytes": 17.09, "DC": 0.20},
         "paediatric": {"CD4+ T cells": 28.55, "CD8+ T cells": 19.90, "B cells": 26.53, "Plasma cells": 1.59, "NK cells": 3.34,
                        "NKT cells": 1.54, "Other cells": 3.49, "Macrophages": 1.10, "Monocytes": 13.70, "DC": 0.24}}
rows = []
for coh in ["adult", "paediatric"]:
    sub = obs[obs.in_paper_scheme & obs.age_group.eq(coh)]
    fr = (sub.celltype_major_ye2026.value_counts(normalize=True) * 100).round(2)
    for k, v in paper[coh].items():
        rows.append((coh, k, v, float(fr.get(k, 0.0))))
    print(coh, "cells in paper scheme:", len(sub))
comp = pd.DataFrame(rows, columns=["cohort", "category", "paper_pct", "recap_pct"]); print(comp.to_string(index=False))
comp.to_csv(f"{OUTD}/composition_vs_paper.csv", index=False)
print("Pearson r (adult, paediatric):", [round(float(np.corrcoef(d.paper_pct, d.recap_pct)[0, 1]), 3) for _, d in comp.groupby("cohort")])
obs["sample_dir"] = obs["gsm"] + "_" + obs["sample"]
obs["cell_id"] = obs.index.astype(str)
keep = ["cell_id", "sample_dir", "barcode", "sample", "gsm", "group", "age_group", "infection_site", "nCount_RNA",
        "nCount_RNA_raw", "nFeature_RNA", "percent_mito", "doublet_score", "major_lineage", "major_leiden_1.5",
        "celltype_major_ye2026", "celltype_ye2026", "recap_subcluster", "qc_flag", "in_paper_scheme"]
obs[keep].to_csv(f"{OUTD}/sepsis_obs.csv.gz", index=False)
print(obs.celltype_ye2026.value_counts().to_string())
