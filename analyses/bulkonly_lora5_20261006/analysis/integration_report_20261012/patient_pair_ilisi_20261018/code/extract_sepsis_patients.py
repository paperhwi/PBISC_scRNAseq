"""Extract N random in-paper-scheme cells per sepsis sample (Ye 2026 recap) on the PBISC sepsis panel genes.
Runs on the PC with h5py + numpy. Output: _extract/sepsis_patients_cells.npz (CSR counts, obs, genes)."""
import os, time, json, h5py, numpy as np
ROOT = os.path.expanduser('~/mnt/PBISC'); A = ROOT + '/bulkonly_lora5_20261006/analysis'
OUT = A + '/integration_report_20261012/patient_pair_ilisi_20261018/_extract'; os.makedirs(OUT, exist_ok=True)
N = int(os.environ.get('NPER', '800'))
def col(g):
    if isinstance(g, h5py.Group) and 'codes' in g:
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g['categories'][:]] + ['NA'], dtype=object); return cats[g['codes'][:]].astype(str)
    v = g[:] if not isinstance(g, h5py.Group) else g['values'][:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v]).astype(str)
t0 = time.time()
f = h5py.File(A + '/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad', 'r'); o = f['obs']; X = f['X']
var = f['var']; sym = col(var['gene_symbol']) if 'gene_symbol' in var else col(var[var.attrs.get('_index', '_index')])
g = h5py.File(A + '/harmony/sepsis_k1000/sepsis_k1000_harmony_embedding.h5ad', 'r'); gvar = col(g['var']['gene_symbol'])
first = {}
for j, s in enumerate(sym): first.setdefault(s, j)
panel = list(dict.fromkeys(x for x in gvar if x in first))
colmap = np.full(len(sym), -1, np.int64)
for k, x in enumerate(panel): colmap[first[x]] = k
sample = col(o['sample']); ok = col(o['in_paper_scheme']) == 'True'
meta = {k: col(o[k]) for k in ['group', 'age_group', 'celltype_major_ye2026', 'celltype_ye2026', 'infection_site']}
indptr = X['indptr'][:].astype(np.int64)
rng = np.random.default_rng(20261018)
D, I, P, rows = [], [], [0], []
for s in sorted(set(sample)):
    ix = np.nonzero(sample == s)[0]; lo, hi = ix.min(), ix.max() + 1
    assert hi - lo == len(ix), f'{s} not contiguous'
    el = ix[ok[ix]]; pick = np.sort(rng.choice(el, min(N, len(el)), replace=False))
    a, b = indptr[lo], indptr[hi]
    d = X['data'][a:b]; ii = X['indices'][a:b]
    for r in pick:
        s0, e0 = indptr[r] - a, indptr[r + 1] - a
        k = colmap[ii[s0:e0]]; m = k >= 0
        D.append(d[s0:e0][m].astype(np.float32)); I.append(k[m].astype(np.int32)); P.append(P[-1] + int(m.sum()))
    rows.append(pick)
rows = np.concatenate(rows)
np.savez_compressed(f'{OUT}/sepsis_patients_cells.npz', data=np.concatenate(D), indices=np.concatenate(I), indptr=np.array(P, np.int64),
                    shape=np.array([len(rows), len(panel)]), genes=np.array(panel), rows=rows, sample=sample[rows],
                    **{k: v[rows] for k, v in meta.items()})
json.dump({'n_cells': int(len(rows)), 'n_samples': int(len(set(sample))), 'n_genes': len(panel), 'N_per_sample': N, 'sec': time.time() - t0},
          open(f'{OUT}/sepsis_patients_done.json', 'w'))
print('done', len(rows), len(panel), time.time() - t0)
