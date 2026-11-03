"""Extract 800 random generated cells per bulk-only sepsis sample (all 5 draws pooled) for K = 1,000 and K = 4,000,
on the same panel-gene order as extract_sepsis_patients.py. h5py + numpy only. usage: python extract_generated_samples.py k1000|k4000"""
import os, sys, json, time, h5py, numpy as np
ROOT = os.path.expanduser('~/mnt/PBISC'); A = ROOT + '/bulkonly_lora5_20261006/analysis'
OUT = A + '/integration_report_20261012/patient_pair_ilisi_20261018/_extract'
K = sys.argv[1]; N = 800
path = {'k1000': A + '/harmony/sepsis_k1000/sepsis_k1000_harmony_embedding.h5ad',
        'k4000': ROOT + '/bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000/sepsis_k4000_harmony_embedding.h5ad'}[K]
def col(g):
    if isinstance(g, h5py.Group) and 'codes' in g:
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g['categories'][:]] + ['NA'], dtype=object); return cats[g['codes'][:]].astype(str)
    v = g[:] if not isinstance(g, h5py.Group) else g['values'][:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v]).astype(str)
t0 = time.time()
ref = np.load(OUT + '/sepsis_patients_cells.npz', allow_pickle=True)['genes'].astype(str)
f = h5py.File(path, 'r'); X = f['layers']['counts']; gvar = col(f['var']['gene_symbol'])
first = {}
for j, s in enumerate(gvar): first.setdefault(s, j)
colmap = np.full(len(gvar), -1, np.int64)
for k, g in enumerate(ref): colmap[first[g]] = k
sample = col(f['obs']['sample_id']); draw = col(f['obs']['draw']); group = col(f['obs']['group'])
rng = np.random.default_rng(20261018 + (1 if K == 'k4000' else 0))
pick = np.sort(np.concatenate([rng.choice(np.nonzero(sample == s)[0], N, replace=False) for s in sorted(set(sample))]))
sel = np.zeros(len(sample), bool); sel[pick] = True
indptr = X['indptr'][:].astype(np.int64); D, I, P = [], [], [0]
for r0 in range(0, len(sample), 20000):
    r1 = min(len(sample), r0 + 20000); a, b = indptr[r0], indptr[r1]
    d = X['data'][a:b]; ii = X['indices'][a:b]
    for r in np.nonzero(sel[r0:r1])[0] + r0:
        s0, e0 = indptr[r] - a, indptr[r + 1] - a; k = colmap[ii[s0:e0]]; m = k >= 0
        D.append(d[s0:e0][m].astype(np.float32)); I.append(k[m].astype(np.int32)); P.append(P[-1] + int(m.sum()))
np.savez_compressed(f'{OUT}/sepsis_{K}_samples_cells.npz', data=np.concatenate(D), indices=np.concatenate(I), indptr=np.array(P, np.int64),
                    shape=np.array([len(pick), len(ref)]), genes=ref, rows=pick, sample=sample[pick], draw=draw[pick], group=group[pick])
print(K, len(pick), time.time() - t0)
