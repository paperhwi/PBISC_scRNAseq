"""Local, read-only input helpers for the 2026-10-14 validation."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'bulkonly_lora5_20261006/analysis/integration_report_20261012'
sys.path[:0] = [str(OUT / 'runtime'), str(ROOT / '.analysis_runtime')]
for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(key, '4')
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad

(OUT / 'results').mkdir(exist_ok=True)
(OUT / 'figs').mkdir(exist_ok=True)

def lin6(x):
    x = str(x)
    if x == 'Plasma cells' or x.startswith('B_ASC'): return 'Plasma'
    if x == 'B cells' or x.startswith('B_'): return 'B'
    if x in ('CD4+ T cells', 'CD8+ T cells', 'NKT cells', 'MAIT', 'NKT', 'T g/d') or x.startswith('T '): return 'T'
    if x in ('NK cells', 'NK') or x.startswith('NK '): return 'NK'
    if x in ('Monocytes', 'Macrophages', 'DC', 'Int.mono', 'ASDC', 'pDC') or x.startswith(('CD14mono', 'CD16mono', 'DC')): return 'Myeloid'
    return 'Other'

def benchmark6(x):
    x = str(x)
    if x in ('B cells', 'Plasma cells'): return 'B cell'
    if x in ('CD4+ T cells', 'CD8+ T cells', 'NKT cells'): return 'T cell'
    if x == 'NK cells': return 'NK cell'
    if x in ('Monocytes', 'Macrophages'): return 'Monocyte'
    if x == 'DC': return 'DC'
    return 'Other'

def read_csr_rows(path, rows, cols=None, node='X', chunk=20000):
    """Stream CSR blocks; do not load full expression or the logcounts layer."""
    rows = np.asarray(rows, dtype=np.int64)
    assert np.all(np.diff(rows) > 0)
    blocks = []
    with h5py.File(path, 'r') as f:
        x = f[node]; ip = x['indptr'][:]; shape = tuple(x.attrs['shape'])
        for start in range(0, shape[0], chunk):
            sel = rows[(rows >= start) & (rows < start + chunk)]
            if not len(sel): continue
            end = min(start + chunk, shape[0]); a, b = int(ip[start]), int(ip[end])
            block = sp.csr_matrix((x['data'][a:b].astype(np.float32), x['indices'][a:b], ip[start:end+1]-a), shape=(end-start, shape[1]))
            block = block[sel-start]
            if cols is not None: block = block[:, cols]
            blocks.append(block)
    return sp.vstack(blocks, format='csr')

def norm(x):
    y = (sp.diags(1e4 / np.maximum(np.asarray(x.sum(1)).ravel(), 1)) @ x).tocsr()
    y.data = np.log1p(y.data)
    return y

def boot_mean(values, rng, n=2000):
    v = np.asarray(values, dtype=float)
    b = v[rng.integers(0, len(v), (n, len(v)))].mean(1)
    return float(v.mean()), *np.quantile(b, [.025, .975]).tolist()

def clr(p):
    q = np.asarray(p) + .001
    l = np.log(q / q.sum(axis=-1, keepdims=True))
    return l - l.mean(axis=-1, keepdims=True)

def sample_props(df, units, label, cats):
    t = df.groupby(units + [label], observed=True).size().unstack(fill_value=0).reindex(columns=cats, fill_value=0)
    return t.div(t.sum(1), axis=0)
