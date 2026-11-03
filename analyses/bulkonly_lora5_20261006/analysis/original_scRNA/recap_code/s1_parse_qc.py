"""Step 1 (Ye 2026 sepsis, GSE279452): per-sample parse, doublet calling and QC.

Follows the authors' pipeline (github.com/rnacenter/Human-Sepsis: Doubletfinder.R, Seurat4Norm.R) and Methods:
  * doublets per sample on the raw Cell Ranger matrix, expected rate 10 %, homotypic-adjusted
    (DoubletFinder nExp_poi.adj = round(0.10*n*(1-homotypic))); Scrublet scores replace DoubletFinder pANN
  * cells kept: 200 < nFeature < 10,000, percent.mito < 20 %
  * mitochondrial genes removed from the expression table
Outputs per sample: <out>/cells/<sample>.npz (CSR counts, QC-passed singlets, all non-MT genes) and
<out>/cells/<sample>.obs.csv; per-gene sums for the streaming HVG step.
"""
import gzip, os, sys, glob, json, time
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc, anndata as ad

SRC = sys.argv[1]            # staged sepsis_GSE279452_full folder
OUT = sys.argv[2]
os.makedirs(f"{OUT}/cells", exist_ok=True)
sc.settings.verbosity = 0

feat_file = sorted(glob.glob(f"{SRC}/GSM8571082_*/*.features.tsv.gz"))[0]
feat = pd.read_csv(feat_file, sep="\t", header=None, names=["gene_id", "gene_symbol", "type"])
mt = feat.gene_symbol.str.upper().str.startswith("MT-").to_numpy()
np.save(f"{OUT}/mt_mask.npy", mt)
feat.to_csv(f"{OUT}/features.csv", index=False)
G = len(feat)

def read_mtx(path):
    with gzip.open(path, "rt") as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        ng, nc, nnz = map(int, line.split())
    df = pd.read_csv(path, sep=" ", header=None, comment="%", skiprows=0, dtype=np.int32,
                     engine="c", compression="gzip")
    df = df.iloc[1:]  # drop the size line
    g = df[0].to_numpy() - 1; c = df[1].to_numpy() - 1; v = df[2].to_numpy()
    del df
    if np.all(np.diff(c) >= 0):
        ip = np.zeros(nc + 1, np.int64); ip[1:] = np.cumsum(np.bincount(c, minlength=nc))
        X = sp.csr_matrix((v.astype(np.float32), g.astype(np.int32), ip), shape=(nc, ng)); X.sort_indices(); return X
    return sp.csr_matrix((v.astype(np.float32), (c, g)), shape=(nc, ng))

WORKER, NW = int(sys.argv[3]), int(sys.argv[4])
def todo():
    ss = sorted(d for d in os.listdir(SRC) if d.startswith("GSM857") and os.path.isdir(f"{SRC}/{d}"))
    return [d for i, d in enumerate(ss) if i % NW == WORKER]
while True:
  done_all = True
  for d in todo():
      if os.path.exists(f"{OUT}/cells/{d}.npz"):
          continue
      mpath = f"{SRC}/{d}/{d}.matrix.mtx.gz"; bpath = f"{SRC}/{d}/{d}.barcodes.tsv.gz"
      if not (os.path.exists(mpath) and os.path.exists(bpath)):
          done_all = False; continue
      try:
          if os.path.getsize(mpath) < 1000: done_all = False; continue
      except OSError: continue
      t0 = time.time()
      try:
          X = read_mtx(mpath); bc = pd.read_csv(bpath, header=None)[0].to_numpy(); assert X.shape == (len(bc), G)
      except Exception as e:
          print('retry later', d, repr(e), flush=True); done_all = False; continue
      bc = pd.read_csv(bpath, header=None)[0].to_numpy()
      assert X.shape == (len(bc), G), (d, X.shape)
      n = X.shape[0]
      # --- doublets (raw matrix, as Doubletfinder.R) ---
      a = ad.AnnData(X.copy())
      sc.pp.filter_genes(a, min_cells=3)
      sc.pp.scrublet(a, expected_doublet_rate=0.10, n_prin_comps=10, random_state=2021)
      score = a.obs["doublet_score"].to_numpy()
      b = ad.AnnData(a.X.copy()); sc.pp.normalize_total(b, 1e4); sc.pp.log1p(b)
      sc.pp.highly_variable_genes(b, flavor="seurat", n_top_genes=2000)
      b = b[:, b.var.highly_variable].copy(); sc.pp.scale(b, max_value=10)
      sc.pp.pca(b, n_comps=10, random_state=2021); sc.pp.neighbors(b, n_pcs=10, n_neighbors=20, random_state=2021)
      sc.tl.leiden(b, resolution=0.8, random_state=2021, flavor="igraph", n_iterations=2)
      p = b.obs.leiden.value_counts(normalize=True).to_numpy(); homotypic = float((p ** 2).sum())
      n_exp = int(round(0.10 * n * (1 - homotypic)))
      order = np.argsort(-score); is_doublet = np.zeros(n, bool); is_doublet[order[:n_exp]] = True
      del a, b
      # --- QC (Seurat4Norm defaults stated in Methods) ---
      n_counts_all = np.asarray(X.sum(1)).ravel()
      mito = np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(n_counts_all, 1)
      Xn = X[:, ~mt]                                   # mitochondrial genes removed
      nfeat = np.diff(X.indptr)                        # nFeature_RNA on the full matrix (before MT removal)
      keep = (~is_doublet) & (nfeat > 200) & (nfeat < 10000) & (mito < 0.20)
      Xk = Xn[keep].tocsr(); Xk.sort_indices()
      obs = pd.DataFrame({
          "barcode": bc, "sample": d.split("_", 1)[1], "gsm": d.split("_", 1)[0],
          "nCount_RNA_raw": n_counts_all, "nFeature_RNA": nfeat, "percent_mito": mito * 100,
          "doublet_score": score, "is_doublet": is_doublet})
      obs["qc_pass"] = keep
      obs.to_csv(f"{OUT}/cells/{d}.obs_all.csv", index=False)
      ok = obs[keep].copy(); ok["nCount_RNA"] = np.asarray(Xk.sum(1)).ravel()
      ok.to_csv(f"{OUT}/cells/{d}.obs.csv", index=False)
      sp.save_npz(f"{OUT}/cells/{d}.tmp.npz", Xk.astype(np.float32), compressed=True); os.replace(f"{OUT}/cells/{d}.tmp.npz", f"{OUT}/cells/{d}.npz")
      print(json.dumps({"sample": d, "n_raw": n, "n_doublet": int(is_doublet.sum()), "homotypic": round(homotypic, 3),
                        "n_keep": int(keep.sum()), "sec": round(time.time() - t0, 1)}), flush=True)

  n_done = len(glob.glob(f"{OUT}/cells/*.npz"))
  if n_done >= 133: break
  time.sleep(30)
