"""Extract stored UMAP coordinates + label columns from the large source .h5ad files into small .npz files
(h5py only; reads obs/obsm, never X). Run from F:/PBISC. Used by sepsis_celltype_umap.ipynb."""
import h5py, numpy as np, os
OUT = "bulkonly_lora5_20261006/analysis/integration_report_20261012/celltype_umap_20261016/data"; os.makedirs(OUT, exist_ok=True)
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])
def grab(path, obsm, cols, extra_obsm=()):
    h = h5py.File(path, "r"); o = h["obs"]; out = {"umap": h[f"obsm/{obsm}"][:].astype(np.float32)}
    for k, name in extra_obsm: out[name] = h[f"obsm/{k}"][:].astype(np.float32)
    for c in cols:
        g = o[c]
        if isinstance(g, h5py.Group) and "codes" in g: out[c + "__codes"] = g["codes"][:]; out[c + "__cats"] = dec(g["categories"][:])
        elif isinstance(g, h5py.Group): out[c + "__raw"] = g["values"][:]; out[c + "__mask"] = g["mask"][:]
        else: v = g[:]; out[c + "__raw"] = dec(v) if v.dtype.kind in "OS" else v
    return out
np.savez_compressed(f"{OUT}/sepsis_real_umap_obs.npz", **grab(
    "bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad", "X_umap",
    ["celltype_ye2026", "celltype_major_ye2026", "in_paper_scheme", "sample", "group", "age_group", "major_lineage"], [("X_pca_harmony", "pca_harmony")]))
for K, p in [("k1000", "bulkonly_lora5_20261006/analysis/harmony/sepsis_k1000/sepsis_k1000_harmony_clustering.h5ad"),
             ("k4000", "bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000/sepsis_k4000_harmony_clustering.h5ad")]:
    np.savez_compressed(f"{OUT}/sepsis_{K}_umap_obs.npz", **grab(p, "X_umap",
        ["celltype_l1_pred", "aifi_l1", "aifi_l2", "aifi_l2_confidence", "draw", "sample_id", "group", "cell_id", "leiden_harmony10_r0.5"]))
