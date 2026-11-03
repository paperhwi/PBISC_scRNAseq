"""Extract cells and genes for the COVID CD4-cycling analysis (runs on the PC with h5py + numpy only).
usage: python extract_cd4cyc.py <source> ; source in gen_k1000, gen_k4000, real_heldout, recap
Writes _extract/<source>_obs.npz (labels, lib, ngenes), _extract/<source>_G.npz (dense counts of the gene set G, all cells)
and _extract/<source>_rows.npz (full CSR rows of selected cells)."""
import sys, os, json, time, h5py, numpy as np
ROOT = os.path.expanduser('~/mnt/PBISC')
A = ROOT + '/bulkonly_lora5_20261006/analysis'
OUT = A + '/integration_report_20261012/cd4_cycling_20261017/_extract'
os.makedirs(OUT, exist_ok=True)
SRC = {
 'gen_k1000': (A + '/harmony/covid_k1000/covid_k1000_harmony_embedding.h5ad', 'layers/counts', 'gene_symbol',
               ['sample_id','draw','group','aifi_l1','aifi_l2','aifi_l2_confidence','celltype_l1_pred','donor_id']),
 'gen_k4000': (ROOT + '/bulkonly_lora5_k4000_20261006/analysis/harmony/covid_k4000/covid_k4000_harmony_embedding.h5ad', 'layers/counts', 'gene_symbol',
               ['sample_id','draw','group','aifi_l1','aifi_l2','aifi_l2_confidence','celltype_l1_pred','donor_id']),
 'real_heldout': (A + '/cell_transcriptome_validation_20261014/inputs/imperial_eskd_real_heldout_panel_counts.h5ad', 'X', 'gene_symbol',
               ['fine','major','group','sample_id','patient_id']),
 'recap': (A + '/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad', 'X', '_index',
               ['cell_type_detailed','cell_type_detailed_recap','cell_type_IFNstim_recap','case_control','sample_id','patient_id','recap_compartment']),
}
G = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gene_set.json')))['genes']

def col(g):
    if isinstance(g, h5py.Group) and 'values' in g:
        v = g['values'][:]
        if 'mask' in g:
            v = v.astype(object); v[g['mask'][:]] = None
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v]).astype(str) if v.dtype.kind in 'SO' else v
    if isinstance(g, h5py.Group):
        c = g['codes'][:]
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g['categories'][:]] + ['NA'], dtype=object)
        return cats[c].astype(str)
    v = g[:]
    if v.dtype.kind in 'SO':
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v]).astype(str)
    return v

def main(src):
    path, xkey, symkey, obskeys = SRC[src]
    t0 = time.time()
    f = h5py.File(path, 'r')
    X = f[xkey]
    n, p = [int(x) for x in X.attrs['shape']]
    sym = col(f['var'][symkey])
    # first occurrence of each symbol
    first = {}
    for j, s in enumerate(sym):
        first.setdefault(s, j)
    gcols = [first.get(g, -1) for g in G]
    colmap = np.full(p, -1, np.int64)
    for k, j in enumerate(gcols):
        if j >= 0:
            colmap[j] = k
    obs = {k: col(f['obs'][k]) for k in obskeys if k in f['obs']}
    obs['obs_names'] = col(f['obs'][f['obs'].attrs.get('_index', '_index')])
    indptr = X['indptr'][:].astype(np.int64)
    lib = np.zeros(n, np.float64); ngen = np.zeros(n, np.int32)
    Gm = np.zeros((n, len(G)), np.float32)
    # selected rows for full CSR
    rng = np.random.default_rng(20261017)
    sel = select_rows(src, obs, n, rng)
    selset = np.zeros(n, bool); selset[sel] = True
    rows_data, rows_ind, rows_ptr = [], [], [0]
    B = int(os.environ.get("BLOCK", "20000"))
    for r0 in range(0, n, B):
        r1 = min(n, r0 + B)
        a, b = indptr[r0], indptr[r1]
        d = X['data'][a:b].astype(np.float32); ix = X['indices'][a:b].astype(np.int64)
        cnt = np.diff(indptr[r0:r1 + 1])
        rid = np.repeat(np.arange(r0, r1), cnt)
        lib[r0:r1] = np.bincount(rid - r0, weights=d, minlength=r1 - r0)
        ngen[r0:r1] = cnt
        k = colmap[ix]; m = k >= 0
        Gm[rid[m], k[m]] = d[m]
        for r in np.nonzero(selset[r0:r1])[0] + r0:
            s, e = indptr[r] - a, indptr[r + 1] - a
            rows_data.append(d[s:e].copy()); rows_ind.append(ix[s:e].astype(np.int32)); rows_ptr.append(rows_ptr[-1] + (e - s))
    np.savez_compressed(f'{OUT}/{src}_obs.npz', lib=lib, ngenes=ngen, **{k: v for k, v in obs.items()})
    np.savez_compressed(f'{OUT}/{src}_G.npz', G=Gm, genes=np.array(G), present=np.array([j >= 0 for j in gcols]))
    np.savez_compressed(f'{OUT}/{src}_rows.npz', data=np.concatenate(rows_data), indices=np.concatenate(rows_ind),
                        indptr=np.array(rows_ptr, np.int64), rows=sel, var=sym, shape=np.array([len(sel), p]))
    json.dump({'src': src, 'n': n, 'p': p, 'n_sel': int(len(sel)), 'G_present': int(sum(j >= 0 for j in gcols)),
               'sec': time.time() - t0}, open(f'{OUT}/{src}_done.json', 'w'))
    print(src, n, p, len(sel), time.time() - t0)

def select_rows(src, obs, n, rng):
    if src.startswith('gen'):
        l2 = obs['aifi_l2']; l1 = obs['celltype_l1_pred']
        prolif = np.nonzero(np.isin(l2, ['Proliferating T cell', 'Proliferating NK cell']))[0]
        t = np.setdiff1d(np.nonzero(l1 == 'T cell')[0], prolif); nk = np.setdiff1d(np.nonzero(l1 == 'NK cell')[0], prolif)
        bb = np.nonzero(l2 == 'Plasma cell')[0]
        sel = np.concatenate([prolif, rng.choice(t, min(12000, len(t)), replace=False), rng.choice(nk, min(4000, len(nk)), replace=False), bb])
    elif src == 'real_heldout':
        fine = obs['fine']; maj = obs['major']
        cyc = np.nonzero(np.isin(maj, ['T CD4 cycling', 'B_ASC_dividing']))[0]
        cd4 = np.nonzero(np.isin(maj, ['T CD4 Helper', 'T CD4 Naive', 'T CD4 CTL', 'T reg']))[0]
        cd8 = np.nonzero(np.isin(maj, ['T CD8 CTL', 'T CD8 CM', 'T CD8 Naive']))[0]
        nk = np.nonzero(np.isin(maj, ['NK', 'NK CD56 bright']))[0]
        sel = np.concatenate([cyc, rng.choice(cd4, 12000, replace=False), rng.choice(cd8, 4000, replace=False), rng.choice(nk, 4000, replace=False)])
    else:
        d = obs['cell_type_detailed']; r = obs['cell_type_detailed_recap']
        cyc = np.nonzero((d == 'T CD4 cycling') | (r == 'T CD4 cycling'))[0]
        cd4set = ['T CD4 Helper', 'T CD4 Naive', 'T CD4 CTL', 'T reg']
        cd4 = np.nonzero(np.isin(r, cd4set) & (d != 'T CD4 cycling'))[0]
        oth = np.nonzero(np.isin(r, ['T CD8 CTL', 'T CD8 CM', 'T CD8 Naive', 'NK', 'B_ASC_dividing']))[0]
        sel = np.concatenate([cyc, rng.choice(cd4, 6000, replace=False), rng.choice(oth, 3000, replace=False)])
    return np.unique(sel)

if __name__ == '__main__':
    main(sys.argv[1])
