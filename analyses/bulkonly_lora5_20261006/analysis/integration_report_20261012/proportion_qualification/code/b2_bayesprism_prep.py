"""Inputs for BayesPrism (Chu et al. 2022): reference = per (sample, fine type) pseudobulk rows of the real sepsis
scRNA (summing rows per state equals summing cells, which is what BayesPrism's reference step does), cell.type = major,
cell.state = fine. Genes: union of the top 50 markers per fine type and per major type, chosen on the run's own
reference samples; ribosomal, mitochondrial and sex-chromosome marker genes removed. Runs: paired fold 0, paired fold 1
(same folds as b_deconv.py), bulk-only (full reference)."""
import json, numpy as np, pandas as pd, glob
O = "/home/claude/work/prop"
idx = pd.read_csv(f"{O}/sepsis_pb_index.csv"); PB = np.load(f"{O}/sepsis_pb_counts.npy", mmap_mode="r"); var = pd.read_csv(f"{O}/sepsis_pb_var.csv")
meta = json.load(open(f"{O}/sepsis_deconv_meta.json"))
files = sorted(glob.glob(f"{O}/bulk/*.htseq.results.gz")); bulk = {}
for f in files:
    s = f.split("/")[-1].split("_", 1)[1].replace(".htseq.results.gz", ""); t = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    bulk[s] = t[~t.index.str.startswith("__")]
bulk = pd.DataFrame(bulk); bulk.index = bulk.index.str.split(".").str[0]
gid = var.gene_id.astype(str).str.split(".").str[0].to_numpy(); sym = var.gene_symbol.astype(str).to_numpy()
bad = pd.Series(sym).str.match(r"^(RPL|RPS|MT-|MRPL|MRPS)").to_numpy() | np.isin(sym, ["XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "MALAT1", "NEAT1"])
inbulk = np.isin(gid, bulk.index) & ~pd.Series(gid).duplicated().to_numpy() & ~bad


def markers(samples, n=50):
    m = idx["sample"].isin(samples).to_numpy(); sub = idx[m]; rows = np.where(m)[0]; mk = set()
    for lv in ["fine", "major"]:
        cats = [c for c in sorted(sub[lv].unique()) if sub.loc[sub[lv] == c, "n_cells"].sum() >= 30]
        S = np.vstack([np.asarray(PB[rows[sub[lv].to_numpy() == c]]).sum(0) for c in cats])
        lc = np.log2(S / S.sum(1, keepdims=True) * 1e6 + 1)
        for j in range(len(cats)):
            sc_ = lc[j] - np.delete(lc, j, 0).max(0); ok = (lc[j] > np.log2(21)) & inbulk
            mk |= set(np.where(ok)[0][np.argsort(-sc_[ok])[:n]])
    return np.array(sorted(mk))


sc_samples = sorted(idx["sample"].unique())
runs = {"fold0": (meta["folds"][0], [s for s in sc_samples if s not in set(meta["folds"][0])]),
        "fold1": (meta["folds"][1], [s for s in sc_samples if s not in set(meta["folds"][1])]),
        "bulkonly": (meta["bulkonly"], sc_samples)}
for name, (ev, ref) in runs.items():
    g = markers(ref); m = idx["sample"].isin(ref).to_numpy() & (idx.n_cells >= 1).to_numpy()
    sub = idx[m].reset_index(drop=True)
    R = pd.DataFrame(np.rint(np.asarray(PB[np.where(m)[0]][:, g])).astype(np.int64), columns=sym[g])
    R.index = [f"r{i}" for i in range(len(R))]
    R.to_csv(f"{O}/bp/{name}_ref.csv.gz")
    sub[["fine", "major", "n_cells"]].assign(row=R.index).to_csv(f"{O}/bp/{name}_ref_labels.csv", index=False)
    Bm = bulk.loc[gid[g], ev].T; Bm.columns = sym[g]; Bm.to_csv(f"{O}/bp/{name}_mix.csv.gz")
    # mean UMI per cell of each type / state in this reference (for the cell-size correction of theta)
    tot = np.asarray(PB[np.where(m)[0]]).sum(1)
    pd.DataFrame({"major": sub.major, "fine": sub.fine, "umi": tot, "n": sub.n_cells}).to_csv(f"{O}/bp/{name}_size.csv", index=False)
    print(name, "genes", len(g), "ref rows", len(R), "mixtures", len(ev), flush=True)
