import sys, h5py, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "/home/claude/work"); import store
f = h5py.File("/home/claude/work/covid/eskd_covid19.h5ad", "r")
X = f["X"]; N, G = X.attrs["shape"]; ip = X["indptr"][:]
def dec(a): return np.array([x.decode() if isinstance(x, bytes) else x for x in a])
var = pd.DataFrame({"gene": dec(f["var/_index"][:])})
obs = {}
for c in f["obs"].keys():
    g = f["obs"][c]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = dec(g["categories"][:]); codes = g["codes"][:]
        obs[c] = pd.Categorical.from_codes(codes, cats)
    elif isinstance(g, h5py.Group) and "values" in g:
        v = g["values"][:].astype(float); v[g["mask"][:]] = np.nan; obs[c] = v
    else:
        v = g[:]; obs[c] = dec(v) if v.dtype.kind in "OS" else v
obs = pd.DataFrame(obs); idx = obs.pop("_index"); obs.index = idx
nonint = 0
def it():
    global nonint
    step = 50000
    for i in range(0, N, step):
        a, b = ip[i], ip[min(i + step, N)]
        d = X["data"][a:b]; ind = X["indices"][a:b]; p = ip[i:min(i + step, N) + 1] - a
        nonint += int((np.abs(d - np.rint(d)) > 1e-3).sum())
        yield sp.csr_matrix((d, ind, p), shape=(len(p) - 1, G))
store.write_store("/home/claude/work/covid/store", it(), var, obs, G)
print("non-integer values:", nonint, "N", N, "G", G)
