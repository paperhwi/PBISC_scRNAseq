"""Per-sample x six-type pseudobulk (sum of raw counts) and cell counts from full files, resumable (2026-10-21).
Usage: python pb_build.py <real_sepsis|real_covid|floor_sepsis|floor_covid|prev_sepsis|prev_covid> [seconds=150]
Genes = floor-file panel (gene_symbol, first occurrence). Real six-type: sepsis celltype_major_ye2026 (in_paper_scheme cells), COVID
cell_type_detailed via t4_real_label_to_six_type.csv. Generated six-type = frozen celltype_l1_pred. Output data/pb/<name>.npz
(S = groups x genes float64, groups = sample||six_type, n = cells per group, genes) + <name>_obs.csv (per-sample metadata)."""
import sys, os, time, numpy as np, pandas as pd, h5py
name = sys.argv[1]; budget = float(sys.argv[2]) if len(sys.argv) > 2 else 150; t0 = time.time()
BASE = os.path.expanduser("~/mnt/PBISC"); AN = f"{BASE}/bulkonly_lora5_20261006/analysis"; IR = f"{AN}/integration_report_20261012"; OUT = f"{IR}/ilisi_floor_20261021/data/pb"
FL = {c: f"{BASE}/new_inference_depth256_K4000_floor1e-4_20261021/" + {"sepsis": "sepsis/PBISC_ye2026_sepsis", "covid": "covid/PBISC_imperial_eskd"}[c] + "_draw0_depth256_K4000_floor0p0001.h5ad" for c in ("sepsis", "covid")}
PATH = {"real_sepsis": f"{AN}/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad", "real_covid": f"{AN}/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad",
        "floor_sepsis": FL["sepsis"], "floor_covid": FL["covid"],
        "prev_sepsis": f"{BASE}/bulkonly_lora5_k4000_20261006/analysis/harmony/sepsis_k4000/sepsis_k4000_harmony_clustering.h5ad",
        "prev_covid": f"{BASE}/bulkonly_lora5_k4000_20261006/analysis/harmony/covid_k4000/covid_k4000_harmony_embedding.h5ad"}[name]
coh = name.split("_")[1]
def col(f, p):
    g = f[p]
    if isinstance(g, h5py.Group):
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object); return c[g["codes"][:]]
    v = g[:]; return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v
ff = h5py.File(FL[coh], "r"); sym = col(ff, "var/gene_symbol").astype(str); genes = pd.unique(sym)
f = h5py.File(PATH, "r")
st_path = f"{OUT}/{name}_state.npz"
if not os.path.exists(st_path):
    if name.startswith("real"):
        vs = col(f, "var/_index").astype(str)
    else:
        vs = col(f, "var/gene_symbol").astype(str)
    pos = pd.Series(np.arange(len(genes)), index=genes); gmap = np.array([pos.get(g, -1) for g in vs], dtype=np.int64)
    if name == "real_sepsis":
        m6 = {"B cells": "B cell", "Plasma cells": "B cell", "CD4+ T cells": "T cell", "CD8+ T cells": "T cell", "NKT cells": "T cell", "NK cells": "NK cell",
              "Monocytes": "Monocyte", "Macrophages": "Monocyte", "DC": "DC", "Other cells": "Other"}
        samp = col(f, "obs/sample").astype(str); six = pd.Series(col(f, "obs/celltype_major_ye2026").astype(str)).map(m6).fillna("Other").to_numpy()
        keep = pd.Series(col(f, "obs/in_paper_scheme")).astype(str).isin(["True", "true", "1"]).to_numpy()
        meta = pd.DataFrame({"sample": samp, "group": col(f, "obs/group").astype(str), "age_group": col(f, "obs/age_group").astype(str), "fine": col(f, "obs/celltype_ye2026").astype(str), "major": col(f, "obs/celltype_major_ye2026").astype(str)})
    elif name == "real_covid":
        t4 = pd.read_csv(f"{IR}/file_component_fig_20261018/results/t4_real_label_to_six_type.csv"); t4 = t4[t4.cohort == "covid"]
        lab = col(f, "obs/cell_type_detailed").astype(str); mp = dict(zip(t4.real_label, t4.six_type)); mp.update({"T CD8 EM": "T cell", "T CD8 EMRA": "T cell"})   # absent from t4 (not in the held-out panel subset)
        six = pd.Series(lab).map(mp).fillna("Other").to_numpy()
        samp = col(f, "obs/sample_id").astype(str); keep = np.ones(len(samp), bool)
        meta = pd.DataFrame({"sample": samp, "patient": col(f, "obs/patient_id").astype(str), "case_control": col(f, "obs/case_control").astype(str),
                             "days_since_posswab": col(f, "obs/days_since_posswab"), "fine": col(f, "obs/cell_type_IFNstim_recap").astype(str), "major": lab})
        print("unmapped covid labels:", sorted(set(lab) - set(mp))[:20])
    else:
        sk = "obs/sample_id" if "sample_id" in f["obs"] else "obs/sample"; samp = col(f, sk).astype(str); six = col(f, "obs/celltype_l1_pred").astype(str)
        keep = np.ones(len(samp), bool); meta = pd.DataFrame({"sample": samp, "group": col(f, "obs/group").astype(str)})
    grp = pd.Series(samp + "||" + six); gk = np.where(keep, pd.Categorical(grp).codes, -1).astype(np.int64); gnames = np.array(pd.Categorical(grp).categories)
    # composition tables (cheap, obs only)
    mk = meta[keep].copy(); mk["six"] = six[keep]
    mk.groupby(["sample", "six"]).size().rename("n").to_csv(f"{OUT}/{name}_comp_six.csv")
    if "fine" in mk: mk.groupby(["sample", "major", "fine"]).size().rename("n").to_csv(f"{OUT}/{name}_comp_fine.csv")
    mk.drop(columns=[c for c in ("fine", "major", "six") if c in mk]).drop_duplicates("sample").to_csv(f"{OUT}/{name}_obs.csv", index=False)
    np.savez(st_path, nxt=0, acc=np.zeros(len(gnames) * len(genes)), gk=gk, gmap=gmap, gnames=gnames)
S = np.load(st_path, allow_pickle=True); nxt = int(S["nxt"]); acc = S["acc"]; gk = S["gk"]; gmap = S["gmap"]; gnames = S["gnames"]; G = len(genes)
key = "X" if (name.startswith("real") or name.startswith("floor")) else "layers/counts"
X = f[key]; ip = X["indptr"]; N = len(ip) - 1; step = 40000 if name == "real_sepsis" else 30000
while nxt < N and time.time() - t0 < budget:
    e = min(nxt + step, N); p = ip[nxt:e + 1]; a, b = int(p[0]), int(p[-1])
    ind = X["indices"][a:b]; dat = X["data"][a:b].astype(np.float64)
    rows = np.repeat(np.arange(nxt, e), np.diff(p)); g = gk[rows]; c = gmap[ind]; m = (g >= 0) & (c >= 0)
    acc += np.bincount(g[m] * G + c[m], weights=dat[m], minlength=len(acc)); nxt = e
np.savez(st_path, nxt=nxt, acc=acc, gk=gk, gmap=gmap, gnames=gnames)
print(f"{name}: rows {nxt:,} / {N:,} ({100*nxt/N:.0f}%) in {time.time()-t0:.0f}s")
if nxt >= N:
    ng = np.bincount(gk[gk >= 0], minlength=len(gnames))
    np.savez_compressed(f"{OUT}/{name}.npz", S=acc.reshape(len(gnames), G).astype(np.float32), n=ng, groups=gnames, genes=genes); print("DONE", name, acc.reshape(len(gnames), G).shape)
