"""Marker-based annotation of the re-clustered myeloid compartment (Stephenson 2025 label scheme).
Cluster -> label decisions are made from the marker tables (c1 runs); author labels are only used afterwards
to report concordance."""
import numpy as np, pandas as pd, anndata as ad
W = "/home/claude/work/covid"
mye = ad.read_h5ad(f"{W}/myeloid.h5ad"); mono = ad.read_h5ad(f"{W}/mono.h5ad"); dc = ad.read_h5ad(f"{W}/dc.h5ad")

# myeloid round 1 (leiden_2): non-monocyte, non-DC populations
mye_map = {"21": "RBC", "24": "HSPC", "27": "Platelet", "29": "Platelet", "26": "CD14mono", "30": "CD14mono"}
# monocyte re-clustering (leiden_2)
mono_map = {**{k: "CD14mono" for k in ["0", "1", "2", "3", "4", "5", "6", "10", "11", "12", "17", "20", "25"]},
            **{k: "CD14mono_IFN" for k in ["7", "8", "9"]},
            **{k: "CD16mono" for k in ["13", "16", "19", "21", "22"]},
            **{k: "CD16mono_IFN" for k in ["14", "24"]},
            "15": "CD16mono_C1", "18": "Int.mono",
            "23": "CD14mono_dex"}           # dexamethasone-associated monocytes (paper Fig. 5)
# DC re-clustering (leiden_1)
dc_map = {"0": "pDC", "1": "pDC", "17": "ASDC", "12": "DC1", "13": "DC3_IFN", "11": "DC3_IFN",
          **{k: "DC3" for k in ["2", "3"]},                      # CD1C+ CD14+ S100A8+ (DC3)
          **{k: "DC2" for k in ["4", "5", "10", "6", "7", "8", "9", "14", "15", "16", "18", "19"]}}  # CD1C+ FCER1A+ CD14-low

lab = pd.Series(index=mye.obs.cell_index.to_numpy(), dtype=object)
r1 = mye.obs.leiden_2.astype(str).map(mye_map)
lab.loc[mye.obs.cell_index.to_numpy()[r1.notna().to_numpy()]] = r1.dropna().to_numpy()
lab.loc[mono.obs.cell_index.to_numpy()] = mono.obs.leiden_2.astype(str).map(mono_map).to_numpy()
lab.loc[dc.obs.cell_index.to_numpy()] = dc.obs.leiden_1.astype(str).map(dc_map).to_numpy()
assert lab.notna().all(), lab.isna().sum()
sub = pd.Series(index=mye.obs.cell_index.to_numpy(), dtype=object)
sub.loc[mono.obs.cell_index.to_numpy()] = ("mono_" + mono.obs.leiden_2.astype(str)).to_numpy()
sub.loc[dc.obs.cell_index.to_numpy()] = ("dc_" + dc.obs.leiden_1.astype(str)).to_numpy()
m = sub.isna(); sub[m] = ("myeloid_" + mye.obs.set_index("cell_index").loc[sub.index[m], "leiden_2"].astype(str)).to_numpy()
out = pd.DataFrame({"cell_type_IFNstim_recap": lab, "recap_subcluster": sub})
out["cell_type_detailed_recap"] = out.cell_type_IFNstim_recap.replace(
    {"CD14mono_IFN": "CD14mono", "CD16mono_IFN": "CD16mono", "CD16mono_C1": "CD16mono", "DC3_IFN": "DC3", "CD14mono_dex": "CD14mono"})
out.to_csv(f"{W}/labels_myeloid.csv")
auth = mye.obs.set_index("cell_index").loc[out.index]
for a_col, r_col in [("cell_type_IFNstim", "cell_type_IFNstim_recap"), ("cell_type_detailed", "cell_type_detailed_recap")]:
    agree = (auth[a_col].astype(str).to_numpy() == out[r_col].astype(str).to_numpy())
    print(r_col, "overall agreement with deposited labels:", round(agree.mean(), 3))
    t = pd.crosstab(auth[a_col].astype(str), out[r_col].astype(str))
    rec = (np.diag(t.reindex(index=t.columns, fill_value=0).loc[t.columns, t.columns]) if False else None)
    per = {}
    for c in sorted(set(auth[a_col].astype(str))):
        mm = auth[a_col].astype(str).to_numpy() == c
        per[c] = (int(mm.sum()), round(float((out[r_col].to_numpy()[mm] == c).mean()), 3))
    print(per)
print(out.cell_type_IFNstim_recap.value_counts().to_dict())
