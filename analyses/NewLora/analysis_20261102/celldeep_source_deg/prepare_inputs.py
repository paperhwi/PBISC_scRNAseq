"""Create donor-balanced raw-count inputs for CellDEEP Real vs New LoRA.

The real source is the scoped full original cohort, using author cell types.
The New LoRA source is the exact 60k eligible label-transfer subset used in
the prior New LoRA analysis. There are no shared donor/sample IDs across them.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "bulkonly_lora5_20261006/analysis"
BASE = ANALYSIS / "integration_report_20261012/ilisi_floor_20261021"
NEW = ROOT / "NewLora/analysis_20261102"
TYPES = ("T cell", "B cell", "NK cell", "Monocyte")
GROUPS = {"sepsis": ("adult_HC", "adult_ICU_ctrl", "adult_sepsis"),
          "covid": ("covid_negative", "covid_positive")}
COVID_EXCL = {"C126_positive_5", "C126_positive_9", "C126_positive_11",
              "C141_negative", "C58_positive_2"}
CELLS_PER_SAMPLE = 30
MAX_SAMPLES_PER_SOURCE = 20
MIN_SAMPLES_PER_SOURCE = 5


def h5values(group):
    if isinstance(group, h5py.Group):
        return np.asarray(group["categories"].asstr()[:])[group["codes"][:]]
    if group.dtype.kind in "SO":
        return group.asstr()[:]
    return group[:]


def original_path(coh):
    if coh == "sepsis":
        return ANALYSIS / "original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad"
    return ANALYSIS / "original_scRNA/covid_eskd/eskd_covid19_recap.h5ad"


def generated_path(coh):
    if coh == "sepsis":
        return ROOT / "NewLora/ye2026_sepsis_draw1_K4000_ceiling.h5ad"
    return ROOT / "NewLora/imperial_eskd_draw0_K4000_ceiling.h5ad"


def real_metadata(coh, path, fine_to_six):
    with h5py.File(path) as f:
        if coh == "sepsis":
            sample = h5values(f["obs/sample"]).astype(str)
            raw_group = h5values(f["obs/group"]).astype(str)
            group = pd.Series(raw_group).map({"sepsis": "adult_sepsis",
                "healthy control": "adult_HC", "ICU control": "adult_ICU_ctrl"}).to_numpy()
            fine = h5values(f["obs/celltype_ye2026"]).astype(str)
            scope = (h5values(f["obs/age_group"]) == "adult") & f["obs/in_paper_scheme"][:].astype(bool)
        else:
            sample = h5values(f["obs/sample_id"]).astype(str)
            raw_group = h5values(f["obs/case_control"]).astype(str)
            group = pd.Series(raw_group).map({"POSITIVE": "covid_positive",
                "NEGATIVE": "covid_negative", "RECOVERY": "recovery"}).to_numpy()
            fine = h5values(f["obs/cell_type_detailed_recap"]).astype(str)
            roster = set(pd.read_csv(ANALYSIS / "integration_report_20261012/results/covid_roster187_local_ids.csv")["sample_id"])
            scope = np.isin(sample, list(roster)) & ~np.isin(sample, list(COVID_EXCL))
        six = pd.Series(fine).map(fine_to_six).fillna("Other").to_numpy()
    rows = np.flatnonzero(scope)
    return pd.DataFrame({"row": rows, "sample": sample[rows],
                         "patient": sample[rows] if coh == "sepsis" else np.array([x.split("_")[0] for x in sample[rows]]),
                         "group": group[rows], "celltype": six[rows]})


def new_metadata(coh, path):
    with h5py.File(path) as f:
        sample = h5values(f["obs/sample_id"]).astype(str)
        group = h5values(f["obs/group"]).astype(str)
        flag = h5values(f["obs/primary_analysis_eligible"]).astype(bool)
        keep = ~np.char.startswith(group, "ped") if coh == "sepsis" else flag
        eligible = np.flatnonzero(keep)
        rng = np.random.default_rng(20261024)
        per = int(np.ceil(60000 / len(np.unique(sample[eligible]))))
        picked = np.concatenate([rng.choice(eligible[np.asarray(ix)],
            size=min(per, len(ix)), replace=False)
            for _, ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
        rows = np.sort(rng.choice(picked, size=min(60000, len(picked)), replace=False))
    labels = pd.read_csv(NEW / coh / "ceiling_cell_labels.csv.gz", index_col=0)
    assert len(rows) == len(labels)
    assert np.array_equal(sample[rows], labels["sample"].to_numpy(dtype=str))
    assert np.array_equal(group[rows], labels["group"].to_numpy(dtype=str))
    obs = pd.DataFrame({"row": rows, "sample": sample[rows],
        "patient": sample[rows] if coh == "sepsis" else np.array([x.split("_")[0] for x in sample[rows]]),
        "group": group[rows], "celltype": labels["six_tr"].to_numpy(dtype=str)})
    if coh == "covid":
        obs = obs[~obs["sample"].isin(COVID_EXCL)].reset_index(drop=True)
    return obs


def gene_map(coh, real_path, new_path):
    target = ad.read_h5ad(BASE / "data/big" / f"{coh}_ref.h5ad", backed="r").var_names.to_numpy(dtype=str)
    with h5py.File(real_path) as r, h5py.File(new_path) as n:
        nsym = h5values(n["var/gene_symbol"]).astype(str)
        nids = h5values(n["var/gene_id"]).astype(str)
        rkey = h5values(r["var/gene_ids"]).astype(str) if coh == "sepsis" else h5values(r["var/_index"]).astype(str)
    npos = dict(zip(nsym[::-1], np.arange(len(nsym))[::-1]))
    rpos = dict(zip(rkey[::-1], np.arange(len(rkey))[::-1]))
    symbols, rc, nc = [], [], []
    for gene in target:
        j = npos.get(gene)
        if j is None:
            continue
        k = rpos.get(nids[j] if coh == "sepsis" else gene)
        if k is None:
            continue
        symbols.append(gene)
        rc.append(k)
        nc.append(j)
    assert len(symbols) > 5000 and len(symbols) == len(set(symbols))
    return np.asarray(symbols), np.asarray(rc, dtype=int), np.asarray(nc, dtype=int)


def choose_visit_per_patient(df):
    counts = df.groupby(["patient", "sample"], sort=True).size().reset_index(name="n")
    counts = counts.sort_values(["patient", "n", "sample"], ascending=[True, False, True])
    best = counts.drop_duplicates("patient")["sample"]
    return df[df["sample"].isin(best)].copy()


def select_units(coh, real, new):
    specs, unit_rows = [], []
    for group in GROUPS[coh]:
        for ct in TYPES:
            key = f"{coh}__{group}__{ct.replace(' ', '')}"
            rng = np.random.default_rng(int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "little"))
            chosen = {}
            for source, df in (("real", real), ("new", new)):
                sub = df[(df.group == group) & (df.celltype == ct)].copy()
                if coh == "covid":
                    sub = choose_visit_per_patient(sub)
                counts = sub.groupby("sample").size()
                valid = sorted(counts.index[counts >= CELLS_PER_SAMPLE])
                chosen[source] = (sub, valid)
            n = min(MAX_SAMPLES_PER_SOURCE, len(chosen["real"][1]), len(chosen["new"][1]))
            if n < MIN_SAMPLES_PER_SOURCE:
                specs.append({"stratum": key, "cohort": coh, "group": group, "celltype": ct,
                              "status": "skipped", "n_real_eligible": len(chosen["real"][1]),
                              "n_new_eligible": len(chosen["new"][1]), "n_samples_each": n})
                continue
            specs.append({"stratum": key, "cohort": coh, "group": group, "celltype": ct,
                          "status": "ready", "n_real_eligible": len(chosen["real"][1]),
                          "n_new_eligible": len(chosen["new"][1]), "n_samples_each": n})
            for source in ("real", "new"):
                sub, valid = chosen[source]
                samples = rng.choice(valid, n, replace=False)
                for sample_id in samples:
                    rows = sub.loc[sub["sample"] == sample_id, "row"].to_numpy(dtype=int)
                    selected = rng.choice(rows, CELLS_PER_SAMPLE, replace=False)
                    patient = str(sub.loc[sub["sample"] == sample_id, "patient"].iloc[0])
                    for row in selected:
                        unit_rows.append({"stratum": key, "source": source, "row": int(row),
                                          "sample_id": str(sample_id), "patient": patient,
                                          "group": group, "celltype": ct})
    return pd.DataFrame(specs), pd.DataFrame(unit_rows)


def extract_csr_rows(path, rows, cols, progress_label):
    """Read only selected rows from on-disk CSR, preserving chosen gene order."""
    unique_rows = np.unique(rows)
    with h5py.File(path) as f:
        x = f["X"]
        ptr = x["indptr"][:]
        mapping = np.full(int(x.attrs["shape"][1]), -1, dtype=np.int32)
        mapping[cols] = np.arange(len(cols), dtype=np.int32)
        out_ptr = np.zeros(len(unique_rows) + 1, dtype=np.int64)
        idx_parts, data_parts = [], []
        t0 = time.time()
        for i, row in enumerate(unique_rows):
            a, b = int(ptr[row]), int(ptr[row + 1])
            ix = x["indices"][a:b]
            vals = x["data"][a:b]
            mapped = mapping[ix]
            keep = mapped >= 0
            idx_parts.append(mapped[keep])
            data_parts.append(vals[keep])
            out_ptr[i + 1] = out_ptr[i] + int(keep.sum())
            if (i + 1) % 5000 == 0:
                print(progress_label, i + 1, len(unique_rows), round(time.time() - t0, 1), flush=True)
    indices = np.concatenate(idx_parts).astype(np.int32, copy=False)
    values = np.concatenate(data_parts)
    assert np.all(np.isfinite(values)) and np.allclose(values, np.rint(values))
    mat = sp.csr_matrix((np.rint(values).astype(np.int32), indices, out_ptr),
                        shape=(len(unique_rows), len(cols)))
    mat.sum_duplicates(); mat.sort_indices()
    return unique_rows, mat


def run(coh):
    real_path, new_path = original_path(coh), generated_path(coh)
    ref = ad.read_h5ad(BASE / "data/big" / f"{coh}_ref.h5ad", backed="r")
    fine_to_six = ref.obs.groupby("fine").six_type.first().to_dict()
    real = real_metadata(coh, real_path, fine_to_six)
    new = new_metadata(coh, new_path)
    specs, selected = select_units(coh, real, new)
    out = OUT / coh
    out.mkdir(parents=True, exist_ok=True)
    specs.to_csv(out / "strata.csv", index=False)
    print(coh, specs.to_string(index=False), flush=True)
    genes, real_cols, new_cols = gene_map(coh, real_path, new_path)
    pd.Series(genes).to_csv(out / "genes.txt", index=False, header=False)
    selected.to_csv(out / "selected_cells.csv", index=False)
    matrices = {}
    for source, path, cols in (("real", real_path, real_cols), ("new", new_path, new_cols)):
        rows, mat = extract_csr_rows(path, selected.loc[selected.source == source, "row"].to_numpy(dtype=int),
                                     cols, f"{coh}/{source}")
        matrices[source] = (rows, mat)
    for key, sub in selected.groupby("stratum", sort=True):
        pieces = []
        meta = []
        for source in ("real", "new"):
            part = sub[sub.source == source].copy()
            rows, mat = matrices[source]
            pos = np.searchsorted(rows, part.row.to_numpy(dtype=int))
            assert np.array_equal(rows[pos], part.row.to_numpy(dtype=int))
            pieces.append(mat[pos])
            meta.append(part)
        X = sp.vstack(pieces, format="csr")
        md = pd.concat(meta, ignore_index=True)
        md.insert(0, "cell_id", [f"cell_{i:05d}" for i in range(len(md))])
        target = out / key
        target.mkdir(exist_ok=True)
        mmwrite(str(target / "counts.mtx"), X.T)
        md.to_csv(target / "metadata.csv", index=False)
        assert X.shape[0] == len(md) and X.shape[1] == len(genes)
    (out / "input_summary.json").write_text(json.dumps({"cohort": coh, "n_genes": len(genes),
        "n_strata_ready": int((specs.status == "ready").sum()), "n_selected_cells": len(selected),
        "cells_per_sample": CELLS_PER_SAMPLE, "max_samples_per_source": MAX_SAMPLES_PER_SOURCE}, indent=2))


if __name__ == "__main__":
    for cohort in ("sepsis", "covid"):
        run(cohort)
