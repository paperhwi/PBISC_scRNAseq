"""Step 2: pack the per-sample QC-passed singlet matrices into one memmap store (uint16 counts / gene index)."""
import sys, glob, os, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "/home/claude/work"); import store
OUT = "/home/claude/work/sepsis/out"
samples = sorted(os.path.basename(f)[:-4] for f in glob.glob(f"{OUT}/cells/GSM*.npz") if ".tmp." not in f)
assert len(samples) == 133, len(samples)
meta = pd.read_csv("/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/samples.csv")
obs = []
for d in samples:
    o = pd.read_csv(f"{OUT}/cells/{d}.obs.csv"); obs.append(o)
obs = pd.concat(obs, ignore_index=True)
obs = obs.merge(meta[["gsm", "group", "age_group", "infection_site"]], on="gsm", how="left")
obs.index = (obs["sample"] + "_" + obs["barcode"]).to_numpy()
f = pd.read_csv(f"{OUT}/features.csv"); mt = np.load(f"{OUT}/mt_mask.npy")
var = f[~mt].reset_index(drop=True)
def it():
    for d in samples:
        yield sp.load_npz(f"{OUT}/cells/{d}.npz")
store.write_store("/home/claude/work/sepsis/store", it(), var, obs, int((~mt).sum()))
print("cells", len(obs))
