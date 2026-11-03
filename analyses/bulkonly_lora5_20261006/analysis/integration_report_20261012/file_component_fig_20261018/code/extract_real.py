"""Sample the ORIGINAL (real) single-cell datasets used for comparison, on the PBISC 6,000-gene panel.
Usage: python extract_real.py <sepsis|covid> <out_dir>
Scope: sepsis = Ye 2026 recap, adult, in_paper_scheme (as in the integration report);
       covid  = Stephenson 2025 recap, 187-sample published roster (results/covid_roster187_local_ids.csv).
20,000 random scoped cells -> raw counts of the 19 markers, per-gene count sum / detection over all 6,000 panel genes,
total UMI (all genes and panel genes); label counts over ALL scoped cells."""
import h5py, numpy as np, sys, os, time, csv
ds, out = sys.argv[1], sys.argv[2]; t0 = time.time()
A = "bulkonly_lora5_20261006/analysis/"
GENES = "CD3E IL7R CCR7 CD4 CD8A NKG7 GNLY MS4A1 CD79A JCHAIN CD14 LYZ S100A8 FCGR3A FCER1A CLEC4C PPBP HBB MKI67".split()
N = 20000
def dec(a): return np.array([x.decode() if isinstance(x, bytes) else x for x in a])
def cat(g): return dec(g["categories"][:]), g["codes"][:]
PB = {"sepsis": "bulkonly_lora5_20261006/ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw0.h5ad",
      "covid": "bulkonly_lora5_20261006/imperial_eskd/lora5_bulkonly_imperial_eskd_draw0.h5ad"}[ds]
with h5py.File(PB, "r") as p:
    c, k = cat(p["var"]["gene_symbol"]); psym = c[k]
    pid = dec(p["var"]["gene_id"][:]); pmeas = p["var"]["measured_in_cohort"][:]
if ds == "sepsis":
    f = h5py.File(A + "original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad", "r"); o = f["obs"]
    c, k = cat(f["var"]["gene_ids"]); rkey = c[k]; pkey = pid
    agec, agek = cat(o["age_group"]); scope = (agec[agek] == "adult") & o["in_paper_scheme"][:]
    finec, finek = cat(o["celltype_ye2026"]); majc, majk = cat(o["celltype_major_ye2026"])
    sc_, sk_ = cat(o["sample"]); grc, grk = cat(o["group"])
    fine = finec[finek]; major = majc[majk]; sample = sc_[sk_]; group = grc[grk]
else:
    f = h5py.File(A + "original_scRNA/covid_eskd/eskd_covid19_recap.h5ad", "r"); o = f["obs"]
    rkey = dec(f["var"]["_index"][:]); pkey = psym.copy()
    # COVID var has symbols only. Bridge PBISC Ensembl ids through the sepsis file's 10x Ensembl->symbol table when the
    # PBISC symbol is not in the COVID file (renamed genes, AL/AC clone names).
    with h5py.File(A + "original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad", "r") as sf:
        ec, ek = cat(sf["var"]["gene_ids"]); yc, yk = cat(sf["var"]["gene_symbol"]); e2s = dict(zip(ec[ek], yc[yk]))
    rset = set(rkey); nb = 0
    for j in range(len(pkey)):
        if pkey[j] not in rset and pid[j] in e2s and e2s[pid[j]] in rset: pkey[j] = e2s[pid[j]]; nb += 1
    print("bridged via Ensembl:", nb)
    roster = set(r["sample_id"] for r in csv.DictReader(open(A + "integration_report_20261012/results/covid_roster187_local_ids.csv")))
    sc_, sk_ = cat(o["sample_id"]); sample = sc_[sk_]; scope = np.isin(sample, list(roster))
    finec, finek = cat(o["cell_type_detailed_recap"]); fine = finec[finek]
    authc, authk = cat(o["cell_type_detailed"]); major = authc[authk]          # deposited author label
    grc, grk = cat(o["case_control"]); group = grc[grk]
first = {}
for i, g in enumerate(rkey):
    first.setdefault(g, i)
pmap = np.array([first.get(g, -1) for g in pkey]); pmap[~pmeas] = -1         # real column per panel gene, -1 = not measured
print("panel genes found in real:", (pmap >= 0).sum(), "| scoped cells:", scope.sum(), flush=True)
col2panel = -np.ones(len(rkey), int)
for j, r in enumerate(pmap):
    if r >= 0 and col2panel[r] < 0: col2panel[r] = j
mk = np.array([GENES.index(g) if g in GENES else -1 for g in psym])
mcol = {q: pmap[np.where(psym == g)[0][0]] for q, g in enumerate(GENES)}
print("markers in real:", sum(v >= 0 for v in mcol.values()))
ip = f["X"]["indptr"]; D = f["X"]["data"]; I = f["X"]["indices"]
cand = np.where(scope)[0]; rows = np.sort(np.random.default_rng(20261018).choice(cand, N, replace=False))
tot_all = np.zeros(N); tot_pan = np.zeros(N); ng_pan = np.zeros(N, int); M = np.zeros((N, len(GENES)), np.int32)
gcnt = np.zeros(len(psym)); gdet = np.zeros(len(psym))
for n, r in enumerate(rows):
    a, b = ip[r], ip[r + 1]; d = D[a:b]; ix = I[a:b]
    tot_all[n] = d.sum(); pj = col2panel[ix]; m = pj >= 0
    tot_pan[n] = d[m].sum(); ng_pan[n] = m.sum()
    np.add.at(gcnt, pj[m], d[m]); np.add.at(gdet, pj[m], 1)
    for q in range(len(GENES)):
        h = np.searchsorted(ix, mcol[q]) if mcol[q] >= 0 else len(ix)
        if h < len(ix) and ix[h] == mcol[q]: M[n, q] = d[h]
    if n % 5000 == 0: print(n, round(time.time() - t0), "s", flush=True)
# duplicated panel symbols share one real column: copy counts to the duplicate positions
for j, r in enumerate(pmap):
    if r >= 0 and col2panel[r] != j: gcnt[j] = gcnt[col2panel[r]]; gdet[j] = gdet[col2panel[r]]
assert np.allclose(M.astype(float), np.round(M)), "non-integer counts"
fu, fc = np.unique(fine[scope], return_counts=True); mu, mc = np.unique(major[scope], return_counts=True)
np.savez_compressed(os.path.join(out, f"real_{ds}.npz"), genes=GENES, gene_symbol=psym, measured=pmap >= 0, rows=rows,
    total_all=tot_all, total_panel=tot_pan, ngenes_panel=ng_pan, marker_counts=M, gene_count_sum=gcnt, gene_det=gdet,
    fine_sample=fine[rows], major_sample=major[rows], sample_sample=sample[rows], group_sample=group[rows],
    fine_cats=fu, fine_counts=fc, major_cats=mu, major_counts=mc, n_scope=int(scope.sum()),
    n_samples=len(np.unique(sample[scope])), n_genes_file=len(rkey), n_cells_file=len(scope))
print("saved", ds, round(time.time() - t0), "s")
