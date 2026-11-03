"""Step 3: comparators computed from the same 122 real bulk profiles.

B1  bulk marker score: mean z of log2(CPM+1) over each state's top Supplementary Table 6 markers.
B2  signature deconvolution (CIBERSORTx-style approximation): 100 cells per fine subset -> mean CP10k signature;
    bulk CPM corrected gene-wise by the median bulk/pseudobulk ratio over paired adult samples; NNLS.
Outputs: results/bulk_comparators.csv (sample x target: b1_score, b2_parent_frac)
"""
import gzip, tarfile, io
import numpy as np, pandas as pd, h5py, anndata as ad
from scipy.optimize import nnls

ROOT = "F:/PBISC/fig6_bulkonly_reproduction_20261020"
TAR = "F:/PBISC/bulkonly_lora5_20261006/analysis/bulk_public/GSE279448/GSE279448_RAW.tar"
REC = "F:/PBISC/bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad"
T6 = "F:/PBISC/paper_claim_analysis/41590_2025_2345_MOESM9_ESM.xlsx"
META = "F:/PBISC/bulkonly_lora5_20261006/ye2026_sepsis/samples.csv"
rng = np.random.default_rng(7)
TARGETS = {"F1": ("CD4T", "T_CD4_c02-GPR183-NR4A2"), "F2": ("CD8T", "T_CD8_c05-GNLY-CCL4hi"),
           "F3": ("NK", "NK_c03-PRF1-CCL4hi"), "F4": ("NKT", "NKT_c05-IFNG-CCL4hi"), "C1": ("CD4T", "T_CD4_c01-LEF1")}
PREFIX = {"CD4T": "T_CD4", "CD8T": "T_CD8", "NK": "NK_c", "NKT": "NKT_"}

# ---- bulk counts (all 207 samples)
cols = {}
with tarfile.open(TAR) as t:
    for m in t.getmembers():
        name = m.name.split("_", 1)[1].replace(".htseq.results.gz", "")
        df = pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())), sep="\t", header=None, index_col=0)[1]
        cols[name] = df[~df.index.str.startswith("__")]
B = pd.DataFrame(cols)
cpm = B / B.sum(0) * 1e6

# ---- symbols
rv = ad.read_h5ad(REC, backed="r").var
id2sym = pd.Series(rv.gene_symbol.astype(str).to_numpy(), index=rv.gene_ids.astype(str).to_numpy())
id2sym = id2sym[~id2sym.index.duplicated()]
cpm = cpm.loc[cpm.index.isin(id2sym.index)]
cpm.index = id2sym.loc[cpm.index].to_numpy()
cpm = cpm.groupby(level=0).sum()

meta = pd.read_csv(META)
bo = meta.sample_id.astype(str).tolist()
assert set(bo) <= set(cpm.columns), set(bo) - set(cpm.columns)

# ---- Table 6 markers
t6 = pd.read_excel(T6, header=1)
t6.columns = ["gene", "p_val", "avg_log2FC", "pct1", "pct2", "p_adj", "cluster"]
T6NAMES = set(t6.cluster.astype(str))
def t6name(cl):  # Table 6 drops the trailing "hi" used in the recap labels (e.g. NK_c03-PRF1-CCL4hi -> NK_c03-PRF1-CCL4)
    return cl if cl in T6NAMES else (cl[:-2] if cl.endswith("hi") and cl[:-2] in T6NAMES else cl)
def top(cl, n=30):
    d = t6[(t6.cluster == t6name(cl)) & (t6.pct1 >= 0.25) & (t6.avg_log2FC > 0)].sort_values("avg_log2FC", ascending=False)
    return [g for g in d.gene if g in cpm.index][:n]

lc = np.log2(cpm[bo] + 1)
z = lc.sub(lc.mean(1), axis=0).div(lc.std(1) + 1e-6, axis=0)
out = []
for k, (par, st) in TARGETS.items():
    gs = top(st)
    sc_ = z.loc[gs].mean(0)
    for s in bo: out.append(dict(sample=s, target=k, b1_score=float(sc_[s]), b1_n_genes=len(gs)))
B1 = pd.DataFrame(out)

# ---- signature + pseudobulk from the recap, one chunked pass
o = ad.read_h5ad(REC, backed="r").obs
lab = o.celltype_ye2026.astype(str).to_numpy(); smp = o["sample"].astype(str).to_numpy()
subsets = [c for c in pd.unique(lab) if c not in ("nan", "Unresolved (low-quality T/NK)", "Platelet-T/NK complex", "Erythrocyte")]
sig_genes = sorted(set(g for c in subsets for g in top(c, 25)))
rsym = rv.gene_symbol.astype(str).to_numpy()
gcol = {g: i for i, g in enumerate(rsym) if g in set(sig_genes)}
sig_genes = [g for g in sig_genes if g in gcol]
colmap = np.full(len(rsym), -1); colmap[[gcol[g] for g in sig_genes]] = np.arange(len(sig_genes))
pick = set(np.concatenate([rng.choice(np.where(lab == c)[0], min(100, (lab == c).sum()), replace=False) for c in subsets]).tolist())
pb_samples = [s for s in np.unique(smp) if s in cpm.columns]
pb_idx = {s: i for i, s in enumerate(pb_samples)}
PB = np.zeros((len(pb_samples), len(sig_genes))); LIB = np.zeros(len(pb_samples))
SIG = {c: np.zeros(len(sig_genes)) for c in subsets}; SIGN = {c: 0 for c in subsets}
with h5py.File(REC, "r") as f:
    X = f["X"]; indptr = X["indptr"][:]
    CH = 40000
    for start in range(0, len(indptr) - 1, CH):
        stop = min(start + CH, len(indptr) - 1); a, b = indptr[start], indptr[stop]
        dat = X["data"][a:b]; ind = X["indices"][a:b]
        for r in range(start, stop):
            s0, s1 = indptr[r] - a, indptr[r + 1] - a
            v = dat[s0:s1]; ii = colmap[ind[s0:s1]]; keep = ii >= 0
            lib = v.sum()
            si = pb_idx.get(smp[r])
            if si is not None:
                np.add.at(PB[si], ii[keep], v[keep]); LIB[si] += lib
            if r in pick:
                vec = np.zeros(len(sig_genes)); vec[ii[keep]] = v[keep] / max(lib, 1) * 1e4
                SIG[lab[r]] += vec; SIGN[lab[r]] += 1
S = pd.DataFrame({c: SIG[c] / max(SIGN[c], 1) for c in subsets}, index=sig_genes)
pbcpm = pd.DataFrame(PB / LIB[:, None] * 1e6, index=pb_samples, columns=sig_genes).T
adult_paired = [s for s in pb_samples if not ("-PS" in s or s.startswith("PHC") or s.startswith("PIC"))]
ratio = (cpm.loc[sig_genes, adult_paired] + 1) / (pbcpm[adult_paired] + 1)
corr = np.exp(np.log(ratio).median(1))
print("signature genes", len(sig_genes), "paired adult samples for correction", len(adult_paired))

out = []
for s in bo:
    y = (cpm.loc[sig_genes, s] / corr).to_numpy()
    w, _ = nnls(S.to_numpy() / 1e4 * 1e6, y)
    fr = pd.Series(w / max(w.sum(), 1e-12), index=subsets)
    for k, (par, st) in TARGETS.items():
        members = [c for c in subsets if c.startswith(PREFIX[par])]
        den = fr[members].sum()
        out.append(dict(sample=s, target=k, b2_parent_frac=float(fr[st] / den) if den > 0 else np.nan, b2_parent_total=float(den)))
B2 = pd.DataFrame(out)
B1.merge(B2, on=["sample", "target"]).to_csv(f"{ROOT}/results/bulk_comparators.csv", index=False)
S.to_csv(f"{ROOT}/results/b2_signature.csv")
print("done")
