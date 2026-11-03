"""Run locked paper-claim signature checks on original and PBISC K datasets.

Scores cells using uncorrected UMI counts normalized on the measured shared gene
panel. Original cluster/state labels define the reference positives; generated
cells are scored with the same frozen marker lists and fixed thresholds. The
generated parent mapping is explicit and intentionally reported as approximate
where AIFI labels do not resolve the paper lineage exactly.
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from sklearn.metrics import roc_curve, roc_auc_score, balanced_accuracy_score
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
BASE = ROOT / "bulkonly_lora5_20261006"
K4 = ROOT / "bulkonly_lora5_k4000_20261006"
CHUNK = 8192
SEED = 20261012

YE_PARENTS = {
    "SEP-1": {"Memory CD4 T cell", "Naive CD4 T cell", "Treg", "Proliferating T cell"},
    "SEP-2": {"Memory CD8 T cell", "Naive CD8 T cell", "CD8aa", "Proliferating T cell"},
    "SEP-3": {"DN T cell"},  # NKT has no dedicated AIFI L2 class; DN T is a coarse candidate only.
    "SEP-3-NK": {"CD56dim NK cell", "CD56bright NK cell", "Proliferating NK cell"},
    "SEP-4A": {"CD14 monocyte", "CD16 monocyte", "Intermediate monocyte"},
    "SEP-4B": {"CD14 monocyte", "CD16 monocyte", "Intermediate monocyte"},
    "SEP-5": {"CD14 monocyte"},
}

COVID_PARENT_MAP = {
    "CD14mono": {"CD14 monocyte"},
    "CD16mono": {"CD16 monocyte"},
    "NK": {"CD56dim NK cell", "CD56bright NK cell"},
    "T CD8 CTL": {"Memory CD8 T cell"},
    "T CD4 Helper": {"Memory CD4 T cell"},
    "T CD4 Naive": {"Naive CD4 T cell"},
}


def decode(a):
    return np.asarray([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def categorical_obs(h5, key, n):
    x = h5["obs"][key]
    if isinstance(x, h5py.Group) and "codes" in x:
        codes = x["codes"][:]
        cats = decode(x["categories"][:])
        out = np.full(n, "", dtype=object)
        valid = codes >= 0
        out[valid] = cats[codes[valid]]
        return out
    if isinstance(x, h5py.Dataset):
        raw = x[:]
        if raw.dtype.kind in "biuf":
            return raw
        return decode(raw)
    raise ValueError(f"Unsupported obs field {key}")


def var_symbols(h5):
    var = h5["var"]
    for key in ["gene_symbol", "gene", "gene_name"]:
        if key in var:
            x = var[key]
            if isinstance(x, h5py.Group) and "codes" in x:
                codes = x["codes"][:]
                cats = decode(x["categories"][:])
                out = np.full(len(codes), "", dtype=object)
                ok = codes >= 0
                out[ok] = cats[codes[ok]]
                return out
            return decode(x[:])
    idx = var["_index"][:] if "_index" in var else h5["var/_index"][:]
    return decode(idx)


def read_panel(h5, panels, counts_mode=True, denominator_symbols=None):
    """Return per-cell module scores for panels, chunked from a sparse CSR matrix."""
    xgroup = h5["X"]
    if not isinstance(xgroup, h5py.Group):
        raise ValueError("Expected CSR sparse X")
    n = len(xgroup["indptr"]) - 1
    symbols = var_symbols(h5)
    # First occurrence is the frozen mapping when a symbol is duplicated.
    gene_to_col = {}
    for j, g in enumerate(symbols):
        gene_to_col.setdefault(str(g), j)
    nvars = len(symbols)
    denominator_cols = None
    if denominator_symbols is not None:
        denominator_cols = np.asarray([gene_to_col[g] for g in denominator_symbols if g in gene_to_col], dtype=np.int32)
    panel_indices = {}
    panel_coverage = {}
    for name, genes in panels.items():
        found = [g for g in genes if g in gene_to_col]
        panel_indices[name] = np.asarray([gene_to_col[g] for g in found], dtype=np.int32)
        panel_coverage[name] = {"n_requested": len(genes), "n_present": len(found), "missing": [g for g in genes if g not in gene_to_col]}
    out = {k: np.zeros(n, dtype=np.float32) for k in panels}
    wr, wc, wv = [], [], []
    names = list(panels)
    for j, name in enumerate(names):
        cols = panel_indices[name]
        if len(cols):
            wr.extend(cols.tolist())
            wc.extend([j] * len(cols))
            wv.extend([1.0 / len(cols)] * len(cols))
    weight = coo_matrix((np.asarray(wv, dtype=np.float32), (wr, wc)), shape=(nvars, len(names))).tocsr()
    indptr = xgroup["indptr"]
    data_ds, idx_ds = xgroup["data"], xgroup["indices"]
    for start in range(0, n, CHUNK):
        stop = min(start + CHUNK, n)
        ptr = indptr[start:stop + 1].astype(np.int64)
        lo, hi = int(ptr[0]), int(ptr[-1])
        ix = idx_ds[lo:hi].astype(np.int32)
        vals = data_ds[lo:hi].astype(np.float32)
        local_ptr = (ptr - lo).astype(np.int64)
        block = csr_matrix((vals, ix, local_ptr), shape=(stop - start, nvars))
        if counts_mode:
            totals = np.asarray((block if denominator_cols is None else block[:, denominator_cols]).sum(axis=1)).ravel().astype(np.float32)
            row_ids = np.repeat(np.arange(stop - start), np.diff(local_ptr))
            block.data = np.log1p(block.data * (10000.0 / np.maximum(totals[row_ids], 1.0)))
        block_scores = (block @ weight).toarray()
        for j, name in enumerate(names):
            out[name][start:stop] = block_scores[:, j]
    return out, panel_coverage


def best_threshold(score, y, groups):
    """Sample-balanced threshold and held-out sample CV performance."""
    score = np.asarray(score, dtype=np.float32)
    y = np.asarray(y, dtype=bool)
    groups = np.asarray(groups, dtype=str)
    valid = np.isfinite(score)
    score, y, groups = score[valid], y[valid], groups[valid]
    uniq = np.unique(groups)
    if len(uniq) < 3 or y.sum() < 20 or (~y).sum() < 20:
        return {"threshold": float(np.quantile(score[y], 0.05)) if y.any() else float("inf"),
                "cv_auc": np.nan, "cv_sensitivity": np.nan, "cv_specificity": np.nan,
                "n_samples": len(uniq), "cv_note": "insufficient positive/negative groups; fallback target 5th percentile"}
    rng = np.random.default_rng(SEED)
    group_indices = {g: np.flatnonzero(groups == g) for g in uniq}
    def sample_for_roc(sample_groups):
        chosen, counts = [], {}
        for grp in sample_groups:
            in_grp = group_indices[grp]
            pos = in_grp[y[in_grp]]
            neg = in_grp[~y[in_grp]]
            if len(pos) > 250:
                pos = rng.choice(pos, 250, replace=False)
            if len(neg) > 1000:
                neg = rng.choice(neg, 1000, replace=False)
            chosen.extend(pos.tolist())
            chosen.extend(neg.tolist())
            counts[grp] = len(pos) + len(neg)
        return np.asarray(chosen, dtype=np.int64), counts
    perm = rng.permutation(uniq)
    fold_of = {g: i % min(5, len(uniq)) for i, g in enumerate(perm)}
    pred = np.zeros(len(y), dtype=bool)
    for fold in range(min(5, len(uniq))):
        test_groups = [g for g in uniq if fold_of[g] == fold]
        train_groups = [g for g in uniq if fold_of[g] != fold]
        test = np.isin(groups, test_groups)
        fit_idx, counts = sample_for_roc(train_groups)
        if not test.any() or y[fit_idx].sum() == 0 or (~y[fit_idx]).sum() == 0:
            continue
        w = np.asarray([1.0 / counts[g] for g in groups[fit_idx]], dtype=np.float64)
        fpr, tpr, th = roc_curve(y[fit_idx], score[fit_idx], sample_weight=w)
        thbest = th[np.argmax(tpr - fpr)]
        pred[test] = score[test] >= thbest
    auc_idx, _ = sample_for_roc(uniq)
    auc = roc_auc_score(y[auc_idx], score[auc_idx])
    sens = float(pred[y].mean()) if y.any() else np.nan
    spec = float((~pred[~y]).mean()) if (~y).any() else np.nan
    fit_idx, counts = sample_for_roc(uniq)
    w = np.asarray([1.0 / counts[g] for g in groups[fit_idx]], dtype=np.float64)
    fpr, tpr, th = roc_curve(y[fit_idx], score[fit_idx], sample_weight=w)
    threshold = float(th[np.argmax(tpr - fpr)])
    return {"threshold": threshold, "cv_auc": float(auc), "cv_sensitivity": sens,
            "cv_specificity": spec, "n_samples": len(uniq), "cv_note": "5-fold sample-held-out; AUC pooled, threshold chosen within training folds"}


def grouped_sums(group_codes, n_groups, mask, values=None):
    codes = group_codes[mask]
    counts = np.bincount(codes, minlength=n_groups)
    sums = None if values is None else np.bincount(codes, weights=np.asarray(values)[mask], minlength=n_groups)
    return counts, sums


def run_ye():
    print("Loading Ye reference and scoring all frozen panels", flush=True)
    defs = pd.read_csv(OUT / "fixed_paper_claim_definitions.csv", dtype=str).fillna("")
    defs = defs[defs.study.eq("Ye et al. 2026")].copy()
    panels = {r.id: list(dict.fromkeys(g for g in r.markers.split(";") if g)) for r in defs.itertuples()}
    source_path = BASE / "analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad"
    src = __import__("anndata").read_h5ad(source_path, backed="r")
    obs = src.obs
    target_defs = {r.id: r.target for r in defs.itertuples()}
    prefix_defs = {r.id: r.lineage_prefix for r in defs.itertuples()}
    inpaper = np.asarray(obs["in_paper_scheme"].astype(bool))
    src_labels = np.asarray(obs["celltype_ye2026"].astype(str), dtype=str)
    src_sample = np.asarray(obs["sample"].astype(str), dtype=str)
    src_age = np.asarray(obs["age_group"].astype(str), dtype=str)
    src_group = np.asarray(obs["group"].astype(str), dtype=str)
    den_path = BASE / "ye2026_sepsis/lora5_bulkonly_ye2026_sepsis_draw0.h5ad"
    with h5py.File(den_path, "r") as gf:
        den_symbols = var_symbols(gf)
        if "measured_in_cohort" in gf["var"]:
            measured = gf["var/measured_in_cohort"][:].astype(bool)
            den_symbols = den_symbols[measured]
    scores, coverage = read_panel(src.file._file, panels, counts_mode=True, denominator_symbols=den_symbols)
    print("Ye reference panel scores complete", flush=True)
    src_rows, cv_rows = [], []
    src_sample_codes, src_sample_names = pd.factorize(src_sample, sort=False)
    src_group_by_sample = {}
    src_age_by_sample = {}
    for code, sid in enumerate(src_sample_names):
        row = np.flatnonzero(src_sample_codes == code)[0]
        src_group_by_sample[sid] = src_group[row]
        src_age_by_sample[sid] = src_age[row]
    for claim in panels:
        parent = inpaper & pd.Series(src_labels).str.startswith(prefix_defs[claim]).to_numpy()
        target = src_labels == target_defs[claim]
        cal = best_threshold(scores[claim][parent], target[parent], src_sample[parent])
        cal.update({"claim_id": claim, "panel_coverage": coverage[claim]["n_present"], "panel_size": coverage[claim]["n_requested"]})
        cv_rows.append(cal)
        npg = len(src_sample_names)
        n_parent, _ = grouped_sums(src_sample_codes, npg, parent)
        n_target, target_sum = grouped_sums(src_sample_codes, npg, parent & target, scores[claim])
        n_bg, bg_sum = grouped_sums(src_sample_codes, npg, parent & ~target, scores[claim])
        n_det, _ = grouped_sums(src_sample_codes, npg, parent & (scores[claim] >= cal["threshold"]))
        _, all_sum = grouped_sums(src_sample_codes, npg, parent, scores[claim])
        for j, sid in enumerate(src_sample_names):
            if not n_parent[j]: continue
            src_rows.append({"study":"Ye et al. 2026", "claim_id":claim, "sample_id":sid,
                             "group":src_group_by_sample[sid], "age_group":src_age_by_sample[sid],
                             "n_parent":int(n_parent[j]), "n_target":int(n_target[j]),
                             "fraction_target":float(n_target[j]/n_parent[j]),
                             "target_score_mean":float(target_sum[j]/n_target[j]) if n_target[j] else np.nan,
                             "background_score_mean":float(bg_sum[j]/n_bg[j]) if n_bg[j] else np.nan,
                             "target_detection_score_fraction":float(n_det[j]/n_parent[j]),
                             "marker_score_parent_mean":float(all_sum[j]/n_parent[j])})
    del scores
    src.file.close()
    source_df = pd.DataFrame(src_rows)
    source_df.to_csv(OUT / "fixed_claim_original_sepsis_parent_frequencies.csv", index=False)
    cv_df = pd.DataFrame(cv_rows)
    cv_df.to_csv(OUT / "fixed_claim_sepsis_thresholds_cv.csv", index=False)
    threshold = dict(zip(cv_df.claim_id, cv_df.threshold))
    gen_rows, gen_group = [], []
    cov_all = []
    for k, root, prefix in [(1000, BASE, "lora5_bulkonly_"), (4000, K4, "lora5k4_bulkonly_")]:
        for draw in range(5):
            p = root / "ye2026_sepsis" / f"{prefix}ye2026_sepsis_draw{draw}.h5ad"
            print(f"Scoring generated Ye K={k} draw={draw}", flush=True)
            with h5py.File(p, "r") as f:
                n = len(f["X/indptr"]) - 1
                obs_cols = {c: categorical_obs(f,c,n) for c in ["sample_id","aifi_l1","aifi_l2","group"]}
                sample_codes, sample_names = pd.factorize(obs_cols["sample_id"], sort=False)
                n_groups = len(sample_names)
                first_by_sample = np.full(n_groups, n, dtype=np.int64)
                np.minimum.at(first_by_sample, sample_codes, np.arange(n, dtype=np.int64))
                scores, coverage = read_panel(f, panels, counts_mode=True)
                cov_all.append({"cohort":"sepsis","K":k,"draw":draw, **{cid: v["n_present"] for cid,v in coverage.items()}})
                for claim in panels:
                    parents = np.isin(obs_cols["aifi_l2"], list(YE_PARENTS[claim]))
                    # SEP-3 has only a coarse DN-T proxy; retain it as an explicitly uncertain estimate.
                    th = threshold[claim]
                    n_parent, _ = grouped_sums(sample_codes, n_groups, parents)
                    n_positive, _ = grouped_sums(sample_codes, n_groups, parents & (scores[claim] >= th))
                    _, score_sum = grouped_sums(sample_codes, n_groups, parents, scores[claim])
                    for j, sid in enumerate(sample_names):
                        if not n_parent[j]: continue
                        first = first_by_sample[j]
                        gen_rows.append({"study":"Ye et al. 2026", "claim_id":claim,"K":k,"draw":draw,"sample_id":sid,
                                         "group":str(obs_cols["group"][first]),
                                         "parent_mapping":";".join(sorted(YE_PARENTS[claim])),"n_parent":int(n_parent[j]),"n_marker_positive":int(n_positive[j]),
                                         "fraction_marker_positive":float(n_positive[j]/n_parent[j]),"marker_score_mean":float(score_sum[j]/n_parent[j]),
                                         "marker_score_median":np.nan})
                del scores
    gen_df = pd.DataFrame(gen_rows)
    gen_df.to_csv(OUT / "fixed_claim_pbisc_sepsis_by_sample_draw.csv", index=False)
    pd.DataFrame(cov_all).to_csv(OUT / "fixed_claim_gene_coverage_by_K_draw.csv", index=False)
    return source_df, cv_df, gen_df


def pubkey(value):
    p, sep, suffix = str(value).strip().partition("_")
    return f"{p}_{suffix[:4]}" if sep and len(suffix) >= 4 else str(value).strip()


def run_covid():
    print("Loading COVID reference on the publisher's 187-sample roster", flush=True)
    source_path = BASE / "analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad"
    gen_ref = BASE / "imperial_eskd/lora5_bulkonly_imperial_eskd_draw0.h5ad"
    with h5py.File(gen_ref, "r") as gf:
        den_symbols = var_symbols(gf)
        if "measured_in_cohort" in gf["var"]:
            den_symbols = den_symbols[gf["var/measured_in_cohort"][:].astype(bool)]
    roster = pd.read_csv(OUT / "covid_published_187_sample_ids.csv", dtype=str)
    roster_keys = set(roster.sample_id.map(pubkey))
    s1 = pd.read_excel(OUT / "Stephenson_Table_S1_sample_metadata.xlsx", sheet_name=0, header=2, dtype=str)
    s1["sample_key"] = s1.sample_id.map(pubkey)
    steroid_by_key = s1.set_index("sample_key")["Received_Steroid"].to_dict()
    IFN = ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6"]
    DEX_UP = ["CD163", "ADAMTS2", "FKBP5", "TSC22D3", "SLC1A3", "MFGE8", "FPR1", "MTSS1", "MRC1", "SESN1"]
    DEX_DOWN = ["JUN", "CXCL8", "IL1B", "HLA-DRA"]
    panels = {"IFN": IFN, "DEX_UP": DEX_UP, "DEX_DOWN": DEX_DOWN}
    src = __import__("anndata").read_h5ad(source_path, backed="r")
    obs = src.obs
    n = src.n_obs
    src_sample = np.asarray(obs["sample_id"].astype(str), dtype=str)
    src_patient = np.asarray(obs["patient_id"].astype(str), dtype=str)
    src_detailed = np.asarray(obs["cell_type_detailed_recap"].astype(str), dtype=str)
    src_ifn_label = np.asarray(obs["cell_type_IFNstim_recap"].astype(str), dtype=str)
    src_case = np.asarray(obs["case_control"].astype(str), dtype=str)
    src_severity = np.asarray(obs["WHO_severity_simple"].astype(str), dtype=str)
    sample_key = np.asarray([pubkey(x) for x in src_sample], dtype=str)
    included = np.isin(sample_key, np.asarray(sorted(roster_keys), dtype=str))
    ifn_positive = pd.Series(src_ifn_label).str.contains(r"IFN|_High$", case=False, regex=True, na=False).to_numpy()
    dex_flag = np.asarray(obs["dex_mono_recap"].fillna(False), dtype=bool)
    mono_mask = pd.Series(src_detailed).str.contains("mono", case=False, regex=False).to_numpy() & included
    scores, coverage = read_panel(src.file._file, panels, counts_mode=True, denominator_symbols=den_symbols)
    signed_dex = scores["DEX_UP"] - scores["DEX_DOWN"]
    src_rows, cv_rows = [], []
    sample_codes, sample_names = pd.factorize(src_sample, sort=False)
    n_samples = len(sample_names)
    first = np.full(n_samples, n, dtype=np.int64)
    np.minimum.at(first, sample_codes, np.arange(n, dtype=np.int64))
    for parent in sorted(COVID_PARENT_MAP):
        parent_mask = included & (src_detailed == parent)
        if not parent_mask.any():
            continue
        cal = best_threshold(scores["IFN"][parent_mask], ifn_positive[parent_mask], src_sample[parent_mask])
        cal.update({"claim_id":"ESKD-IFN", "parent_cell_type":parent, "panel_coverage":coverage["IFN"]["n_present"], "panel_size":len(IFN)})
        cv_rows.append(cal)
        n_parent, _ = grouped_sums(sample_codes, n_samples, parent_mask)
        n_ifn, _ = grouped_sums(sample_codes, n_samples, parent_mask & ifn_positive)
        n_sig, _ = grouped_sums(sample_codes, n_samples, parent_mask & (scores["IFN"] >= cal["threshold"]))
        _, sc_sum = grouped_sums(sample_codes, n_samples, parent_mask, scores["IFN"])
        _, positive_sum = grouped_sums(sample_codes, n_samples, parent_mask & ifn_positive, scores["IFN"])
        _, negative_sum = grouped_sums(sample_codes, n_samples, parent_mask & ~ifn_positive, scores["IFN"])
        for j, sid in enumerate(sample_names):
            if not n_parent[j] or pubkey(sid) not in roster_keys: continue
            i = first[j]
            src_rows.append({"claim_id":"ESKD-IFN", "parent_cell_type":parent, "sample_id":sid, "sample_key":pubkey(sid),
                             "patient_id":src_patient[i], "case_control":src_case[i], "severity":src_severity[i],
                             "received_steroid":steroid_by_key.get(pubkey(sid), ""), "n_parent":int(n_parent[j]), "n_IFN_label":int(n_ifn[j]),
                             "fraction_IFN_label":float(n_ifn[j]/n_parent[j]), "fraction_signature_positive":float(n_sig[j]/n_parent[j]),
                             "mean_IFN_score":float(sc_sum[j]/n_parent[j]),
                             "IFN_score_label_mean":float(positive_sum[j]/n_ifn[j]) if n_ifn[j] else np.nan,
                             "IFN_score_nonlabel_mean":float(negative_sum[j]/(n_parent[j]-n_ifn[j])) if n_parent[j]>n_ifn[j] else np.nan})
    dex_cal = best_threshold(signed_dex[mono_mask], dex_flag[mono_mask], src_sample[mono_mask])
    dex_cal.update({"claim_id":"ESKD-DEX", "parent_cell_type":"all monocytes", "panel_coverage":min(coverage["DEX_UP"]["n_present"],coverage["DEX_DOWN"]["n_present"]),
                    "panel_size":len(DEX_UP)+len(DEX_DOWN)})
    cv_rows.append(dex_cal)
    n_mono, _ = grouped_sums(sample_codes, n_samples, mono_mask)
    n_dex, _ = grouped_sums(sample_codes, n_samples, mono_mask & dex_flag)
    n_dexsig, _ = grouped_sums(sample_codes, n_samples, mono_mask & (signed_dex >= dex_cal["threshold"]))
    _, up_sum = grouped_sums(sample_codes, n_samples, mono_mask, scores["DEX_UP"])
    _, down_sum = grouped_sums(sample_codes, n_samples, mono_mask, scores["DEX_DOWN"])
    for j, sid in enumerate(sample_names):
        if not n_mono[j] or pubkey(sid) not in roster_keys: continue
        i = first[j]
        src_rows.append({"claim_id":"ESKD-DEX", "parent_cell_type":"all monocytes", "sample_id":sid,"sample_key":pubkey(sid),
                         "patient_id":src_patient[i],"case_control":src_case[i],"severity":src_severity[i],
                         "received_steroid":steroid_by_key.get(pubkey(sid), ""),"n_parent":int(n_mono[j]),"n_IFN_label":int(n_dex[j]),
                         "fraction_IFN_label":float(n_dex[j]/n_mono[j]),"fraction_signature_positive":float(n_dexsig[j]/n_mono[j]),
                         "mean_IFN_score":float(up_sum[j]/n_mono[j]),"DEX_down_mean":float(down_sum[j]/n_mono[j])})
    src_df = pd.DataFrame(src_rows)
    cv_df = pd.DataFrame(cv_rows)
    src_df.to_csv(OUT / "fixed_claim_covid_reference_by_sample.csv", index=False)
    cv_df.to_csv(OUT / "fixed_claim_covid_thresholds_cv.csv", index=False)
    thresholds = {(r.claim_id, getattr(r,"parent_cell_type","")):r.threshold for r in cv_df.itertuples()}
    del scores
    src.file.close()

    gen_rows = []
    mono_aifi = {"CD14 monocyte", "CD16 monocyte", "Intermediate monocyte"}
    for k, root, prefix in [(1000, BASE, "lora5_bulkonly_"), (4000, K4, "lora5k4_bulkonly_")]:
        for draw in range(5):
            path = root / "imperial_eskd" / f"{prefix}imperial_eskd_draw{draw}.h5ad"
            print(f"Scoring generated COVID K={k} draw={draw}", flush=True)
            with h5py.File(path, "r") as f:
                nn = len(f["X/indptr"]) - 1
                o = {c:categorical_obs(f,c,nn) for c in ["sample_id","aifi_l2","group"]}
                scs, cov = read_panel(f, panels, counts_mode=True)
                sample_codes, sample_names = pd.factorize(o["sample_id"], sort=False)
                ng = len(sample_names)
                first = np.full(ng, nn, dtype=np.int64)
                np.minimum.at(first, sample_codes, np.arange(nn, dtype=np.int64))
                for parent, gen_parent_set in COVID_PARENT_MAP.items():
                    if parent not in set(src_detailed[included]): continue
                    parent_score = scs["IFN"]
                    threshold = thresholds[("ESKD-IFN", parent)]
                    pm = np.isin(o["aifi_l2"], list(gen_parent_set))
                    npar, _ = grouped_sums(sample_codes, ng, pm)
                    nsig, _ = grouped_sums(sample_codes, ng, pm & (parent_score >= threshold))
                    _, sums = grouped_sums(sample_codes, ng, pm, parent_score)
                    for j,sid in enumerate(sample_names):
                        if not npar[j]: continue
                        q = first[j]
                        gen_rows.append({"claim_id":"ESKD-IFN","parent_cell_type":parent,"K":k,"draw":draw,"sample_id":sid,
                                         "group":o["group"][q],"parent_mapping":";".join(sorted(gen_parent_set)),"n_parent":int(npar[j]),
                                         "n_signature_positive":int(nsig[j]),"fraction_signature_positive":float(nsig[j]/npar[j]),
                                         "mean_IFN_score":float(sums[j]/npar[j])})
                pm = np.isin(o["aifi_l2"], list(mono_aifi))
                signed = scs["DEX_UP"] - scs["DEX_DOWN"]
                th = thresholds[("ESKD-DEX","all monocytes")]
                npar, _ = grouped_sums(sample_codes,ng,pm)
                nsig, _ = grouped_sums(sample_codes,ng,pm & (signed >= th))
                _, su = grouped_sums(sample_codes,ng,pm,scs["DEX_UP"])
                _, sd = grouped_sums(sample_codes,ng,pm,scs["DEX_DOWN"])
                for j,sid in enumerate(sample_names):
                    if not npar[j]: continue
                    q=first[j]
                    gen_rows.append({"claim_id":"ESKD-DEX","parent_cell_type":"all monocytes","K":k,"draw":draw,"sample_id":sid,
                                     "group":o["group"][q],"parent_mapping":";".join(sorted(mono_aifi)),"n_parent":int(npar[j]),
                                     "n_signature_positive":int(nsig[j]),"fraction_signature_positive":float(nsig[j]/npar[j]),
                                     "mean_DEX_up_score":float(su[j]/npar[j]),"mean_DEX_down_score":float(sd[j]/npar[j])})
                del scs
    gen_df = pd.DataFrame(gen_rows)
    gen_df.to_csv(OUT / "fixed_claim_covid_pbisc_by_sample_draw.csv", index=False)
    return src_df,cv_df,gen_df


def main():
    s, cv, g = run_ye()
    print("Ye original rows", len(s), "CV", cv[["claim_id","threshold","cv_auc","cv_sensitivity","cv_specificity","panel_coverage","panel_size"]].to_dict("records"))
    print("Ye generated rows",len(g),"groups",g.groupby(["claim_id","K"]).fraction_marker_positive.agg(["median","mean","std"]).reset_index().to_dict("records"))
    cs, ccv, cg = run_covid()
    print("COVID reference",cs.shape,"thresholds",ccv[["claim_id","parent_cell_type","threshold","cv_auc","cv_sensitivity","cv_specificity"]].to_dict("records"))
    print("COVID generated",cg.shape,"summary",cg.groupby(["claim_id","parent_cell_type","K"]).fraction_signature_positive.agg(["median","mean","std"]).reset_index().to_dict("records"))


if __name__ == "__main__":
    main()
