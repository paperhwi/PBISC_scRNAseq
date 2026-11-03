"""Cluster -> Ye 2026 subset maps used for the CD4, CD8/NK/NKT and 'other' re-clusterings (Leiden 1.5 / 2.5 / 1.5).
Decisions from the marker tables printed by s4_markers.py (see validation notebook, section 3)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/sepsis"
# ---- CD4 (run_cd4, leiden_1.5)
P = f"{W}/run_cd4"; cells = np.load(f"{P}/embed.npz")["cells"]; cl = np.load(f"{P}/leiden.npz")["leiden_1.5"].astype(str)
m = {**{k: "T_CD4_c01-LEF1" for k in ["3","4","8","13","16","18","21","14","22"]},          # CCR7/SELL/LEF1, S100A4 low
     **{k: "T_CD4_c02-GPR183-NR4A2" for k in ["0","2","10","12","15","19"]},               # NR4A2 > ANXA2 (NR4A2/AREG/ARID5A)
     **{k: "T_CD4_c03-GPR183-ANXA2" for k in ["1","5","7"]},                               # ANXA2/S100A4 > NR4A2
     "11": "T_CD4_c04-GZMK", "17": "T_CD4_c05-GNLY", "6": "T_CD4_c06-FOXP3", "20": "T_CD4_c07-MKI67",
     "9": "T_CD8_c01-LEF1", "23": "HSC"}                                                    # 9: CD8B/LINC02446; 23: PRSS57/SPINK2
pd.DataFrame({"cell": cells, "lineage": "T/NK", "celltype_ye2026": [m[c] for c in cl], "subcluster": ["CD4_" + c for c in cl],
              "qc_flag": ["B/myeloid_doublet_suspect" if c == "14" else "" for c in cl]}).to_csv(f"{W}/labels_cd4.csv", index=False)
# ---- CD8 / NK / NKT / gdT / MAIT (run_cyto, leiden_2.5)
st = store.Store(f"{W}/store"); g = st.var.gene_symbol.astype(str).to_numpy()
P = f"{W}/run_cyto"; cells = np.load(f"{P}/embed.npz")["cells"]; c25 = np.load(f"{P}/leiden.npz")["leiden_2.5"].astype(str)
m = {"16": "NK_c01-GZMK", "2": "NK_c02-PRF1", "11": "NK_c02-PRF1", "3": "NK_c02-PRF1", "22": "NK_c02-PRF1", "30": "NK_c02-PRF1",
     "4": "NK_c03-PRF1-CCL4hi", "26": "NKT_c01-CX3CR1", "23": "NKT_c02-CXCR4", "21": "NKT_c03-HLA-DRA", "15": "NKT_c04-IFNG",
     "17": "NKT_c04-IFNG", "19": "NKT_c05-IFNG-CCL4hi",
     **{k: "T_CD8_c01-LEF1" for k in ["0", "28", "31", "33"]}, **{k: "T_CD8_c02-GPR183" for k in ["7", "10", "27", "32"]},
     **{k: "T_CD8_c03-GZMK" for k in ["6", "18", "14"]}, "1": "T_CD8_c04-GNLY", "5": "T_CD8_c04-GNLY",
     **{k: "T_CD8_c05-GNLY-CCL4hi" for k in ["8", "12", "29"]}, "9": "T_CD8_c07-TRAV1-2", "20": "gdT",
     "13": "Unresolved (low-quality T/NK)"}
lab = np.array([m.get(c, "") for c in c25], dtype=object); cyc = np.isin(c25, ["24", "25"])   # cycling: split by CD3
gi = [int(np.where(g == x)[0][0]) for x in ["CD3D", "CD3E"]]
E = np.vstack([store.lognorm(X, st.ncount[c])[:, gi].toarray() for c, X in st.chunks(cells[cyc])])
lab[np.where(cyc)[0]] = np.where(E.sum(1) > 0, "T_CD8_c06-MKI67", "NK_c04-MKI67")
pd.DataFrame({"cell": cells, "lineage": "T/NK", "celltype_ye2026": lab, "subcluster": ["Cyto_" + c for c in c25],
              "qc_flag": np.where(c25 == "13", "low_quality", "")}).to_csv(f"{W}/labels_cyto.csv", index=False)
# ---- other (run_other, leiden_1.5): megakaryocyte / RBC / HSC / neutrophil / platelet complexes
P = f"{W}/run_other"; cells = np.load(f"{P}/embed.npz")["cells"]; cl = np.load(f"{P}/leiden.npz")["leiden_1.5"].astype(str)
mk = {"0","1","2","3","7","8","12","15","17","19","21","23","25"}; rbc = {"4","11","14","16","26"}
lab = ["HSC" if c == "22" else "Megakaryocyte" if c in mk else "Erythrocyte" if c in rbc else "Neutrophil" if c in {"10","20"}
       else "B_c01-TCL1A" if c == "18" else "Platelet-T/NK complex" for c in cl]
flag = ["" if c == "22" or c in mk else "not_in_paper_scheme(RBC filter)" if c in rbc else "not_in_paper_scheme" if c in {"10","20"}
        else "platelet_complex" for c in cl]
pd.DataFrame({"cell": cells, "lineage": "Other", "celltype_ye2026": lab, "subcluster": ["Oth_" + c for c in cl], "qc_flag": flag}).to_csv(f"{W}/labels_other.csv", index=False)
