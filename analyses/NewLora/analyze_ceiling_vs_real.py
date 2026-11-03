"""Pairwise real-A/NewLora Harmony, iLISI, label transfer and composition.

Run: python analyze_ceiling_vs_real.py sepsis|covid
The input H5ADs are left untouched. Cell sampling uses seed 20261024;
PCA, Harmony, Leiden and UMAP use seed 0, and iLISI uses seeds 0–2.
"""
import json
import sys
import time
from pathlib import Path

import anndata as ad
import harmonypy
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
sys.path.insert(0, str(BASE / "code"))
import ilisi_lib as L

COH = sys.argv[1]
assert COH in ("sepsis", "covid")
INPUT = ROOT / "NewLora" / ({"sepsis": "ye2026_sepsis_draw1_K4000_ceiling.h5ad", "covid": "imperial_eskd_draw0_K4000_ceiling.h5ad"}[COH])
OUT = ROOT / "NewLora/analysis_20261102" / COH
OUT.mkdir(parents=True, exist_ok=True)
START = time.time()

def log(*msg):
    print(f"[{time.time()-START:7.1f}s]", COH, *msg, flush=True)

def column(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        categories = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object)
        return categories[g["codes"][:]]
    a = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object) if a.dtype.kind in "SO" else a

def read_csr_rows(f, rows, cols, step=20000):
    g = f["X"]
    ptr = g["indptr"][:]
    ncol = int(g.attrs["shape"][1])
    blocks = []
    for start in range(0, int(rows.max()) + 1, step):
        selected = rows[(rows >= start) & (rows < start + step)]
        if len(selected) == 0:
            continue
        end = min(start + step, len(ptr) - 1)
        a, b = int(ptr[start]), int(ptr[end])
        block = sp.csr_matrix((g["data"][a:b].astype(np.float32), g["indices"][a:b], ptr[start:end + 1] - a), shape=(end-start, ncol))
        blocks.append(block[selected-start][:, cols])
    return sp.vstack(blocks, format="csr")

ref = ad.read_h5ad(BASE / "data/big" / f"{COH}_ref.h5ad")
real_b = ad.read_h5ad(BASE / "data/big" / f"{COH}_realU.h5ad", backed="r")
with h5py.File(INPUT) as f:
    sample = column(f, "obs/sample_id").astype(str)
    group = column(f, "obs/group").astype(str)
    eligible_flag = column(f, "obs/primary_analysis_eligible").astype(bool)
    keep = (~np.char.startswith(group, "ped")) if COH == "sepsis" else eligible_flag
    eligible = np.where(keep)[0]
    rng = np.random.default_rng(20261024)
    per_sample = int(np.ceil(60000 / len(np.unique(sample[eligible]))))
    picked = np.concatenate([rng.choice(eligible[np.asarray(ix)], size=min(per_sample, len(ix)), replace=False)
                             for _, ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
    rows = np.sort(rng.choice(picked, size=min(60000, len(picked)), replace=False))
    symbols = column(f, "var/gene_symbol").astype(str)
    positions = pd.Series(np.arange(len(symbols)), index=symbols)
    positions = positions[~positions.index.duplicated()]
    cols = positions.reindex(ref.var_names)
    assert cols.notna().all(), f"Missing reference genes: {cols.isna().sum()}"
    X = read_csr_rows(f, rows, cols.to_numpy(dtype=int))
    obs = pd.DataFrame({"sample": sample[rows], "group": group[rows]}, index=[f"ceiling_{i}" for i in range(len(rows))])
    n_input = len(sample)
log("sampled", len(rows), "of", n_input, "cells; genes", X.shape[1], "mean UMI", round(float(X.sum()/X.shape[0]), 1))

XR = ref.X.tocsr().astype(np.float32)
XS = X.tocsr()
YR, YS = L.lognorm(XR), L.lognorm(XS)
hv = L.hvg_seurat_v3([XR])
ref_scaled = YR[:, hv].toarray().astype(np.float32)
mu, sd = ref_scaled.mean(0), ref_scaled.std(0) + 1e-6
ref_scaled = np.clip((ref_scaled-mu)/sd, -10, 10)
pca = PCA(30, random_state=0).fit(ref_scaled)
del ref_scaled
pre_r = pca.transform(np.clip((YR[:, hv].toarray().astype(np.float32)-mu)/sd, -10, 10)).astype(np.float32)
pre_s = pca.transform(np.clip((YS[:, hv].toarray().astype(np.float32)-mu)/sd, -10, 10)).astype(np.float32)
pre = np.vstack((pre_r, pre_s)); del pre_r, pre_s, YR, YS, XR, XS, X
nR, nS = ref.n_obs, len(rows)
src = np.r_[np.zeros(nR, dtype=np.int8), np.ones(nS, dtype=np.int8)]
log("PCA complete; Harmony start")
harmony = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({"source": src.astype(str)}), ["source"], random_state=0, max_iter_harmony=20, verbose=False)
z = np.asarray(harmony.Z_corr)
post = (z.T if z.shape[0] == pre.shape[1] else z).astype(np.float32)
del harmony, z
log("Harmony complete")

result = {"cohort": COH, "input": str(INPUT), "n_input": n_input, "n_ref": nR, "n_ceiling": nS,
          "n_eligible": int(keep.sum()),
          "n_genes": ref.n_vars, "hvg": 2000, "pca": 30, "harmony_iterations": 20,
          "seed": 20261024, "label_method": "30-NN transfer from real A in pairwise Harmony space"}
for name, space in (("pre", pre), ("post", post)):
    scores = [L.lisi_stats(space, src, seed=s) for s in (0,1,2)]
    result[f"iLISI_{name}"] = float(np.mean([s[0] for s in scores]))
    result[f"iLISI_{name}_query_mean"] = float(np.mean([s[1] for s in scores]))
    result[f"own_share_{name}"] = float(np.mean([s[2] for s in scores]))
    log(name, "iLISI", round(result[f"iLISI_{name}"], 3))

fine = ref.obs.fine.astype(str).to_numpy()
six = ref.obs.six_type.astype(str).to_numpy()
major = ref.obs.major.astype(str).to_numpy()
f2six = pd.crosstab(fine, six).idxmax(axis=1).to_dict()
f2major = pd.crosstab(fine, major).idxmax(axis=1).to_dict()
_, nearest = NearestNeighbors(n_neighbors=30, n_jobs=4).fit(post[:nR]).kneighbors(post[nR:])
cats, codes = np.unique(fine, return_inverse=True)
vote = np.zeros((nS, len(cats)), dtype=np.uint8)
for j in range(30):
    np.add.at(vote, (np.arange(nS), codes[nearest[:, j]]), 1)
tf = cats[vote.argmax(1)]
obs["fine_tr"] = tf
obs["major_tr"] = [f2major[t] for t in tf]
obs["six_tr"] = [f2six[t] for t in tf]
obs["confidence"] = vote.max(1)/30
result["median_transfer_confidence"] = float(obs.confidence.median())
log("label transfer complete")

joint = ad.AnnData(np.empty((nR+nS, 0), dtype=np.float32))
joint.obsm["X_pca"] = pre
joint.obsm["X_pca_harmony"] = post
sc.pp.neighbors(joint, use_rep="X_pca_harmony", n_neighbors=15, random_state=0)
sc.tl.leiden(joint, resolution=1.0, flavor="igraph", n_iterations=2, random_state=0)
lc = joint.obs.leiden.astype(int).to_numpy()
ct = pd.crosstab(lc, src).reindex(columns=[0,1], fill_value=0)
ct.columns = ["n_real_A", "n_ceiling"]
fractions = ct.div(ct.sum(0), axis=1)
cl = ct.copy()
cl["fraction_real_A"] = fractions.n_real_A
cl["fraction_ceiling"] = fractions.n_ceiling
cl["delta_pp"] = 100*(cl.fraction_ceiling-cl.fraction_real_A)
cl["source_share_normalized"] = cl.fraction_ceiling/(cl.fraction_ceiling+cl.fraction_real_A)
dominant = pd.DataFrame({"cluster": lc[:nR], "fine": fine}).groupby("cluster").fine.agg(lambda x: x.value_counts().index[0])
cl["dominant_real_fine"] = cl.index.map(dominant).fillna("no real cells")
cl.to_csv(OUT / "joint_leiden_clusters.csv", index_label="leiden")
obs["leiden_joint"] = lc[nR:]
result["n_leiden_clusters"] = len(cl)
result["cluster_fraction_pearson"] = float(np.corrcoef(cl.fraction_real_A, cl.fraction_ceiling)[0,1])
p, q = cl.fraction_real_A.to_numpy(), cl.fraction_ceiling.to_numpy()
m = (p+q)/2
with np.errstate(divide="ignore", invalid="ignore"):
    result["cluster_fraction_JSD"] = float((np.sum(np.where(p>0,p*np.log2(p/m),0))+np.sum(np.where(q>0,q*np.log2(q/m),0)))/2)
result["ceiling_fraction_in_mixed_clusters"] = float(cl.loc[cl.source_share_normalized.between(.2,.8), "fraction_ceiling"].sum())
log("Leiden complete", len(cl), "clusters")

# Composition uses equal weight per sample for cohort summaries; raw cell fractions are also retained.
real_obs = real_b.obs.copy()
real_obs["six_tr"] = real_obs.six_type.astype(str)
real_transfer = pd.read_csv(BASE / "results/big" / f"{COH}_realU_transfer.csv.gz", index_col=0)
matched_fine = real_transfer.reindex(real_obs.index)["fine_tr"]
assert matched_fine.notna().all()
real_obs["fine_tr"] = matched_fine.astype(str)
real_obs["sample"] = real_obs["sample"].astype(str)
real_obs["group"] = real_obs["group"].astype(str)
obs.to_csv(OUT / "ceiling_cell_labels.csv.gz", compression="gzip")
all_types = sorted(set(obs.six_tr) | set(real_obs.six_tr))
all_fine = sorted(set(obs.fine_tr) | set(real_obs.fine_tr))
def comp(frame, field, categories, source):
    sample_group = frame[["sample","group"]].drop_duplicates()
    counts = pd.crosstab(frame["sample"], frame[field]).reindex(columns=categories, fill_value=0)
    rates = counts.div(counts.sum(1), axis=0)
    long = rates.stack().rename("fraction").reset_index().rename(columns={field:"cluster"})
    long["group"] = long["sample"].map(sample_group.set_index("sample").group)
    long["source"] = source
    long["level"] = field
    return long
by_sample = pd.concat([comp(obs, "six_tr", all_types, "ceiling"), comp(real_obs, "six_tr", all_types, "real_B"),
                       comp(obs, "fine_tr", all_fine, "ceiling"), comp(real_obs, "fine_tr", all_fine, "real_B")], ignore_index=True)
by_sample.to_csv(OUT / "composition_by_sample.csv", index=False)
summary = by_sample.groupby(["level","source","cluster"]).fraction.agg(["mean","std","count"]).reset_index()
summary.to_csv(OUT / "composition_summary.csv", index=False)
six_pivot = summary[summary.level=="six_tr"].pivot(index="cluster", columns="source", values="mean").fillna(0)
result["six_type_TV_vs_real_B"] = float(.5*np.abs(six_pivot.ceiling-six_pivot.real_B).sum())
for source, frame in (("ceiling", obs), ("real_B", real_obs)):
    six_counts = frame.six_tr.value_counts(normalize=True).reindex(all_types, fill_value=0)
    result[f"six_cell_fraction_{source}"] = {k: float(v) for k,v in six_counts.items()}
log("composition complete")

# Export embedding without expression counts and plot a fixed subsample.
joint.obs = pd.concat([ref.obs.reset_index(drop=True)[["sample","group","fine","six_type"]].assign(source="real_A"),
                       obs.reset_index(drop=True)[["sample","group","fine_tr","six_tr"]].rename(columns={"fine_tr":"fine","six_tr":"six_type"}).assign(source="ceiling")], ignore_index=True)
joint.write_h5ad(OUT / "joint_harmony_embedding.h5ad", compression="gzip")
prng = np.random.default_rng(20261024)
subset = np.r_[prng.choice(nR, min(15000,nR), False), nR+prng.choice(nS, min(15000,nS), False)]
plot_ad = joint[subset].copy()
sc.pp.neighbors(plot_ad, use_rep="X_pca_harmony", n_neighbors=15, random_state=0)
sc.tl.umap(plot_ad, random_state=0)
u = plot_ad.obsm["X_umap"]
fig, axes = plt.subplots(1,2,figsize=(11,4.5))
for source, color in (("real_A", "#8895A7"), ("ceiling", "#D45D48")):
    mask = plot_ad.obs.source.to_numpy()==source
    axes[0].scatter(u[mask,0],u[mask,1],s=.6,alpha=.55,c=color,label=source,rasterized=True)
axes[0].legend(frameon=False, markerscale=6)
axes[0].set_title(f"{COH}: real A + ceiling, Harmony\niLISI {result['iLISI_pre']:.2f} → {result['iLISI_post']:.2f}")
palette = dict(zip(all_types, plt.cm.tab10(np.linspace(0,1,len(all_types)))))
for label, color in palette.items():
    mask = plot_ad.obs.six_type.to_numpy()==label
    axes[1].scatter(u[mask,0],u[mask,1],s=.6,alpha=.55,c=[color],label=label,rasterized=True)
axes[1].legend(frameon=False,markerscale=6,fontsize=7,bbox_to_anchor=(1.01,1),loc="upper left")
axes[1].set_title("Real author labels / ceiling transferred labels")
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(OUT / "umap_harmony.png", dpi=220, bbox_inches="tight"); plt.close(fig)
real_b.file.close()
(OUT / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
log("saved", OUT)
