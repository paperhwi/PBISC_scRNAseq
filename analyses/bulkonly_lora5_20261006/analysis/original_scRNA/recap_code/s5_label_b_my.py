"""Ye 2026 labels for the B/plasma and myeloid lineages (re-clustered, Leiden 1.5).
Cluster -> subset decisions follow the paper's subset definitions (Fig. 1c names; Results text for marker logic):
 B: naive TCL1A+ (c01), TCL1A+CD83+ (c02), TCL1A+TNF+ (c03), TLN1+ FCRL5+/ITGAX+ (c04), AIM2+ memory (c05),
    AIM2+ CD27-high memory (c06), MZB1+ plasma (c07), MKI67+ plasmablast (c08).
 CD14 monocytes: ordered rule on cluster-mean programme scores built from the paper's defining genes:
    c02 IL1B/CXCL8/CCL3/CCL4 > 1; c05 IFITM1/ISG15 > 0.9; c01 HLA-DRA/DPA1/DPB1 > 1.25 with S100 < 3.1;
    c03 S100A8/S100A12/RETN >= 3.8 with HLA < 0.8; remaining CD14 clusters (VCAN/ITGA4/TLN1) -> c04. MKI67+ CD14 -> c06; CD16 -> c07; CD14+CD16+ -> c08;
    C1QA+ -> Macrophage; CD1C/FCER1A/CLEC10A -> cDC; LILRA4 -> pDC."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/sepsis"
st = store.Store(f"{W}/store"); g = st.var.gene_symbol.astype(str).to_numpy()
rows = []
# ---------------- B / plasma
P = f"{W}/run_b"; cells = np.load(f"{P}/embed.npz")["cells"]; cl = np.load(f"{P}/leiden.npz")["leiden_1.5"].astype(str)
bmap = {**{k: "B_c01-TCL1A" for k in ["1", "3", "5", "8", "15", "20", "21", "10"]},
        "2": "B_c02-TCL1A-CD83", "6": "B_c02-TCL1A-CD83", "0": "B_c03-TCL1A-TNF",
        "9": "B_c04-TLN1", "4": "B_c05-AIM2", "19": "B_c05-AIM2", "7": "B_c06-AIM2-CD27hi",
        **{k: "B_c07-MZB1" for k in ["11", "17", "16", "12", "14"]}, "13": "B_c08-MKI67", "18": "pDC"}
flag = {"10": "platelet_complex", "16": "platelet_complex", "12": "T_doublet_suspect", "14": "myeloid_doublet_suspect"}
rows.append(pd.DataFrame({"cell": cells, "lineage": "B/Plasma", "celltype_ye2026": [bmap[c] for c in cl],
                          "subcluster": ["B_" + c for c in cl], "qc_flag": [flag.get(c, "") for c in cl]}))
# ---------------- myeloid
P = f"{W}/run_my"; cells = np.load(f"{P}/embed.npz")["cells"]; cl = np.load(f"{P}/leiden.npz")["leiden_1.5"].astype(str)
z = np.load(f"{P}/stats_leiden_1.5.npz"); groups = list(z["groups"]); mean = z["mean"]
def gm(genes): return mean[:, [int(np.where(g == x)[0][0]) for x in genes]].mean(1)
sig = {"Mono_c01-CD14-HLA-DRA": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1"], "Mono_c02-CD14-IL1B": ["IL1B", "CXCL8", "CCL3", "CCL4"],
       "Mono_c03-CD14-RETN": ["S100A8", "S100A12", "RETN"], "Mono_c04-CD14-TLN1": ["VCAN", "ITGA4", "TLN1"],
       "Mono_c05-CD14-IFITM1": ["IFITM1", "ISG15"]}
fixed = {"18": "Mono_c06-CD14-MKI67", "6": "Mono_c07-CD16", "13": "Mono_c07-CD16", "8": "Mono_c08-CD14CD16",
         "14": "Macrophage", "17": "cDC", "21": "pDC"}
cd14 = [c for c in groups if c not in fixed]
T = pd.DataFrame({"HLA": gm(["HLA-DRA", "HLA-DPA1", "HLA-DPB1"]), "IL1B": gm(["IL1B", "CXCL8", "CCL3", "CCL4"]),
                  "S100": gm(["S100A8", "S100A12", "RETN"]), "IFN": gm(["IFITM1", "ISG15"])}, index=groups).loc[cd14]
def rule(r):      # ordered: dominant inflammatory / IFN programme first, then MHC-II-high, then S100-high/HLA-low
    if r.IL1B > 1.0: return "Mono_c02-CD14-IL1B"
    if r.IFN > 0.9: return "Mono_c05-CD14-IFITM1"
    if r.HLA > 1.25 and r.S100 < 3.1: return "Mono_c01-CD14-HLA-DRA"
    if r.S100 >= 3.8 and r.HLA < 0.8: return "Mono_c03-CD14-RETN"
    return "Mono_c04-CD14-TLN1"
auto = {c: rule(T.loc[c]) for c in cd14}
mymap = {**auto, **fixed}
flag = {"10": "platelet_complex", "16": "platelet_complex", "11": "T_doublet_suspect", "12": "T_doublet_suspect",
        "23": "doublet_suspect"}
print("CD14 cluster assignment:", {k: v.split("-")[0] + "-" + v.split("-")[-1] for k, v in auto.items()})
rows.append(pd.DataFrame({"cell": cells, "lineage": "Myeloid", "celltype_ye2026": [mymap[c] for c in cl],
                          "subcluster": ["My_" + c for c in cl], "qc_flag": [flag.get(c, "") for c in cl]}))
out = pd.concat(rows, ignore_index=True); out.to_csv(f"{W}/labels_b_my.csv", index=False)
print(out.celltype_ye2026.value_counts().to_dict())
