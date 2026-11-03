"""BayesPrism + 5-donor calibration on the paired folds: the bulk mixture is multiplied, gene by gene, by the
pseudobulk/bulk ratio of the same 5 calibration donors that WNNLS+cal used in b_deconv.py (donor draw reproduced from
its seed and recorded here), then rounded to counts. Reference and genes as in b2_bayesprism_prep.py."""
import json, glob, numpy as np, pandas as pd
O = "/home/claude/work/prop"
idx = pd.read_csv(f"{O}/sepsis_pb_index.csv"); PB = np.load(f"{O}/sepsis_pb_counts.npy", mmap_mode="r"); var = pd.read_csv(f"{O}/sepsis_pb_var.csv")
meta = json.load(open(f"{O}/sepsis_deconv_meta.json"))
files = sorted(glob.glob(f"{O}/bulk/*.htseq.results.gz")); bulk = {}
for f in files:
    s = f.split("/")[-1].split("_", 1)[1].replace(".htseq.results.gz", ""); t = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    bulk[s] = t[~t.index.str.startswith("__")]
bulk = pd.DataFrame(bulk); bulk.index = bulk.index.str.split(".").str[0]
gid = var.gene_id.astype(str).str.split(".").str[0].to_numpy(); common = np.array([g for g in pd.unique(gid) if g in bulk.index])
gpos = pd.Series(np.arange(len(gid)), index=gid); gpos = gpos[~gpos.index.duplicated()]; gi = gpos.loc[common].to_numpy(); Bk = bulk.loc[common].astype(float)
# reproduce b_deconv.py's random stream: permutation of paired samples, then one 5-donor draw per fold
rng = np.random.default_rng(20261013); paired = meta["paired_adult"]; perm = rng.permutation(paired)
folds = [list(perm[: len(perm) // 2]), list(perm[len(perm) // 2:])]; assert [list(map(str, f)) for f in folds] == meta["folds"]
cal = {}
for k, ev in enumerate(folds):
    cal[f"fold{k}"] = list(map(str, rng.choice([s for s in paired if s not in set(ev)], 5, replace=False)))
json.dump(cal, open(f"{O}/bp/cal_donors.json", "w"), indent=1)
sym_by_gid = dict(zip(gid, var.gene_symbol.astype(str)))
for name, donors in cal.items():
    r = []
    for s in donors:
        pbs = np.asarray(PB[np.where(idx["sample"].to_numpy() == s)[0]][:, gi]).sum(0)
        r.append(np.log((pbs / pbs.sum() * 1e6 + 1) / (Bk[s].to_numpy() / Bk[s].sum() * 1e6 + 1)))
    fac = pd.Series(np.exp(np.mean(r, 0)), index=[sym_by_gid[g] for g in common]); fac = fac[~fac.index.duplicated()]
    mix = pd.read_csv(f"{O}/bp/{name}_mix.csv.gz", index_col=0)
    mix = np.rint(mix * fac.reindex(mix.columns).fillna(1).to_numpy()[None, :]).astype(int)
    mix.to_csv(f"{O}/bp/{name}cal_mix.csv.gz")
    for suf in ["_ref.csv.gz", "_ref_labels.csv", "_size.csv"]:
        import shutil; shutil.copy(f"{O}/bp/{name}{suf}", f"{O}/bp/{name}cal{suf}")
    print(name, donors, flush=True)
