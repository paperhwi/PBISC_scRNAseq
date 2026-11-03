"""Deterministic 20k-cell file-level extract for the New LoRA overview figures.

The real/floor inputs are the 2026-10-22 extracts used by the two reference
figures. New LoRA sampling is uniform over all cells in its single output draw.
"""
from pathlib import Path
import sys
import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
OLD = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/file_component_fig_20261018/_extract"
FILES = {"sepsis": "ye2026_sepsis_draw1_K4000_ceiling.h5ad",
         "covid": "imperial_eskd_draw0_K4000_ceiling.h5ad"}
MARKERS = "CD3E IL7R CCR7 CD4 CD8A NKG7 GNLY MS4A1 CD79A JCHAIN CD14 LYZ S100A8 FCGR3A FCER1A CLEC4C PPBP HBB MKI67".split()

def col(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        cat = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
        return cat[g["codes"][:]]
    a = g[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else x for x in a]) if a.dtype.kind in "SO" else a

def extract(cohort):
    real = np.load(OLD / f"real_{cohort}.npz", allow_pickle=True)
    floor = np.load(OLD / f"{cohort}_floor.npz", allow_pickle=True)
    with h5py.File(ROOT / "NewLora" / FILES[cohort], "r") as f:
        sym = col(f, "var/gene_symbol").astype(str)
        assert np.array_equal(sym, real["gene_symbol"])
        assert np.array_equal(sym, floor["gene_symbol"])
        measured = f["uns/assay_availability/mask"][:].any(0)
        assert np.array_equal(measured, floor["measured"])
        ip = f["X/indptr"][:]
        n = len(ip)-1
        sample = col(f, "obs/sample_id").astype(str)
        group = col(f, "obs/group").astype(str)
        rng = np.random.default_rng(20261103)
        rows = np.sort(rng.choice(n, size=20000, replace=False))
        gidx = np.array([np.flatnonzero(sym==g)[0] for g in MARKERS])
        totals=[]; detected=[]; marker=[]
        gene_sum=np.zeros(len(sym), dtype=np.float64)
        gene_det=np.zeros(len(sym), dtype=np.int64)
        x=f["X"]
        for start in range(0,n,5000):
            end=min(start+5000,n)
            pick=rows[(rows>=start)&(rows<end)]
            if not len(pick): continue
            a,b=int(ip[start]),int(ip[end])
            block=sp.csr_matrix((x["data"][a:b],x["indices"][a:b],ip[start:end+1]-a),shape=(end-start,len(sym)))
            sub=block[pick-start].tocsr()
            totals.append(np.asarray(sub.sum(axis=1)).ravel())
            detected.append(np.asarray((sub>0).sum(axis=1)).ravel())
            marker.append(sub[:,gidx].toarray())
            gene_sum+=np.asarray(sub.sum(axis=0)).ravel()
            gene_det+=np.asarray((sub>0).sum(axis=0)).ravel().astype(np.int64)
            if start%100000==0: print(cohort,start,"/",n,flush=True)
    total=np.concatenate(totals); ngenes=np.concatenate(detected); markers=np.concatenate(marker)
    assert len(total)==20000 and len(ngenes)==20000 and markers.shape==(20000,len(MARKERS))
    assert np.isclose(total.sum(),gene_sum.sum())
    np.savez_compressed(OUT/f"new_lora_{cohort}_20k.npz",rows=rows,total=total,ngenes=ngenes,
        marker_counts=markers,gene_count_sum=gene_sum,gene_det=gene_det,gene_symbol=sym,
        measured=measured,n_cells=n,n_samples=len(np.unique(sample)),sample_sample=sample[rows],group_sample=group[rows],
        marker_genes=np.asarray(MARKERS))
    print(cohort,"saved",n,"cells",len(np.unique(sample)),"samples; median UMI",np.median(total),flush=True)

if __name__=="__main__":
    cohorts=sys.argv[1:] or ["sepsis","covid"]
    for c in cohorts:
        assert c in FILES
        extract(c)
