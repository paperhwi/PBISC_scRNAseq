"""Stream original sparse matrices and score frozen paper marker panels."""
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "bulkonly_lora5_20261006" / "analysis" / "original_scRNA"
OUT = ROOT / "paper_claim_analysis"
CHUNK = 8192

_ye_defs = pd.read_csv(OUT / "fixed_paper_claim_definitions.csv", dtype=str).fillna("")
SEP = {
    r["id"]: (r["target"], r["lineage_prefix"], [g for g in r["markers"].split(";") if g])
    for _, r in _ye_defs[_ye_defs["study"].eq("Ye et al. 2026")].iterrows()
}
IFN_GENES = ["ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6"]
DEX_UP = ["CD163", "ADAMTS2", "FKBP5", "TSC22D3", "SLC1A3", "MFGE8", "FPR1", "MTSS1", "MRC1", "SESN1"]
DEX_DOWN = ["JUN", "CXCL8", "IL1B", "HLA-DRA"]


def decode(values):
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values], dtype=object)


def publication_sample_key(value):
    prefix, sep, suffix = str(value).strip().partition("_")
    return f"{prefix}_{suffix[:4]}" if sep and len(suffix) >= 4 else str(value).strip()


def chunk_scores(f, matrix_key, genes, n, counts_mode):
    """Yield (row start, row stop, selected-gene scores) in H5AD order."""
    matrix = f[matrix_key]
    gene_names = decode(f["var/_index"][:])
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    present = [g for g in genes if g in gene_to_idx]
    indices = np.asarray([gene_to_idx[g] for g in present], dtype=np.int64)
    indptr = matrix["indptr"]
    data_ds, idx_ds = matrix["data"], matrix["indices"]
    for start in range(0, n, CHUNK):
        stop = min(n, start + CHUNK)
        ptr = indptr[start:stop + 1].astype(np.int64)
        lo, hi = int(ptr[0]), int(ptr[-1])
        block_data = data_ds[lo:hi]
        block_indices = idx_ds[lo:hi]
        block_ptr = (ptr - lo).astype(np.int64)
        # Only selected marker columns are materialized as dense values.
        vals = np.zeros((stop - start, len(present)), dtype=np.float64)
        if len(present):
            block = csr_matrix((block_data, block_indices, block_ptr), shape=(stop-start, len(gene_names)))
            vals = block[:, indices].toarray()
        if counts_mode:
            totals = np.zeros(stop-start, dtype=np.float64)
            for row in range(stop-start):
                a, b = block_ptr[row], block_ptr[row+1]
                if b > a:
                    totals[row] = np.asarray(block_data[a:b]).sum()
            denom = np.maximum(totals, 1.0)
            vals = np.log1p(vals * (10000.0 / denom[:, None]))
        yield start, stop, present, vals


def score_groups(group_arrays, score, panel, kind, study):
    """Aggregate cell scores per biological sample for two within-parent groups."""
    rows = []
    samples = group_arrays["sample"]
    for sample in pd.unique(samples):
        sm = samples == sample
        pm = sm & group_arrays["parent"]
        tmask = pm & group_arrays["target"]
        bmask = pm & ~group_arrays["target"]
        if not tmask.any() or not bmask.any():
            continue
        rows.append({"study":study,"claim_id":kind,"sample":sample,
                     "parent_cells":int(pm.sum()),"target_cells":int(tmask.sum()),
                     "target_marker_score":float(score[tmask].mean()),
                     "parent_background_marker_score":float(score[bmask].mean()),
                     "paired_score_difference":float(score[tmask].mean()-score[bmask].mean()),
                     "marker_panel":";".join(panel)})
    return rows


def sepsis_scores():
    obsdir = BASE / "sepsis_GSE279452_full" / "recap_annotation"
    all_obs = pd.concat([pd.read_csv(p, usecols=["sample","age_group","celltype_ye2026","in_paper_scheme"], low_memory=False)
                         for p in sorted(obsdir.glob("sepsis_obs_part*.csv.gz"))], ignore_index=True)
    claims = {k:(v[0],v[1],v[2]) for k,v in SEP.items()}
    genes = sorted({g for _,_,gs in claims.values() for g in gs})
    h5 = BASE / "sepsis_GSE279452_full" / "sepsis_GSE279452_recap.h5ad"
    full_scores = {k:np.zeros(len(all_obs),dtype=np.float32) for k in claims}
    with h5py.File(h5,"r") as f:
        if len(f["obs/_index"]) != len(all_obs):
            raise RuntimeError("Sepsis obs/H5AD row counts cannot be aligned")
        # The recap has in-paper cells interspersed with excluded rows; process all
        # rows, but retain scores only for rows in the paper scheme.
        sample = all_obs["sample"].astype(str).to_numpy()
        labels = all_obs["celltype_ye2026"].astype(str).to_numpy()
        inpaper = all_obs["in_paper_scheme"].astype(bool).to_numpy()
        for start,stop,present,vals in chunk_scores(f,"X",genes,len(sample),counts_mode=True):
            if not present:
                continue
            idx = {g:i for i,g in enumerate(present)}
            for claim,(target,prefix,panel) in claims.items():
                score = vals[:,[idx[g] for g in panel if g in idx]].mean(axis=1)
                # Store cell scores in a full-length vector for sample-level aggregation.
                full_scores[claim][start:stop] = score.astype(np.float32)
    # Aggregate only original target-lineage rows after the streaming pass.
    out=[]
    all_obs = pd.concat([pd.read_csv(p,usecols=["sample","age_group","celltype_ye2026","in_paper_scheme"],low_memory=False)
                         for p in sorted(obsdir.glob("sepsis_obs_part*.csv.gz"))],ignore_index=True)
    labels=all_obs.celltype_ye2026.astype(str).to_numpy(dtype=str)
    ip=all_obs.in_paper_scheme.astype(bool).to_numpy()
    for claim,(target,prefix,panel) in claims.items():
        arrays={"sample":all_obs["sample"].astype(str).to_numpy(),
                "parent":ip & np.asarray([x.startswith(prefix) for x in labels]),
                "target":labels == target}
        out.extend(score_groups(arrays,full_scores[claim],panel,claim,"Ye et al. 2026"))
    result=pd.DataFrame(out)
    result.to_csv(OUT/"original_sepsis_marker_scores_by_sample.csv",index=False)
    summaries=[]
    for claim,g in result.groupby("claim_id"):
        dif=g.paired_score_difference.dropna()
        p=wilcoxon(dif).pvalue if len(dif)>0 and np.any(dif!=0) else np.nan
        summaries.append({"claim_id":claim,"n_samples_with_target_and_background":len(g),
                          "median_target_score":g.target_marker_score.median(),
                          "median_parent_background_score":g.parent_background_marker_score.median(),
                          "median_paired_difference":dif.median(),"positive_sample_differences":int((dif>0).sum()),
                          "paired_wilcoxon_pvalue":p,"marker_panel":g.marker_panel.iloc[0]})
    pd.DataFrame(summaries).to_csv(OUT/"original_sepsis_marker_score_summary.csv",index=False)


def covid_scores():
    h5=BASE/"covid_eskd"/"eskd_covid19_recap.h5ad"
    rec=pd.read_csv(BASE/"covid_eskd"/"recap_annotation"/"covid_obs_recap.csv.gz")
    roster=pd.read_csv(OUT/"covid_published_187_sample_ids.csv",dtype={"sample_id":str})
    keep_sample_keys=set(roster["sample_id"].map(publication_sample_key))
    with h5py.File(h5,"r") as f:
        n=len(f["obs/_index"])
        def cat(k):
            x=f["obs"][k]
            codes=x["codes"][:]; cats=decode(x["categories"][:]); out=np.full(n,"",dtype=object); ok=codes>=0; out[ok]=cats[codes[ok]]; return out
        detailed=cat("cell_type_detailed")
        sample=cat("sample_id")
        patient=cat("patient_id")
        disease=cat("case_control")
        severity=cat("WHO_severity")
        rec_ids=rec["_index"].astype(str).to_numpy()
        ids=decode(f["obs/_index"][:])
        if not np.array_equal(rec_ids,ids):
            raise RuntimeError("COVID recap annotations do not align to source H5AD rows")
        ifn=rec["cell_type_IFNstim_recap"].astype(str).str.contains(r"IFN|_High$",case=False,regex=True,na=False).to_numpy()
        dex=rec["dex_mono_recap"].fillna(False).astype(bool).to_numpy()
        mono=np.char.find(np.char.lower(detailed.astype(str)),"mono")>=0
        genes=IFN_GENES+DEX_UP+DEX_DOWN
        scores={name:np.zeros(n,dtype=np.float32) for name in ["IFN","DEX_UP","DEX_DOWN"]}
        for start,stop,present,vals in chunk_scores(f,"layers/logcounts",genes,n,counts_mode=False):
            ix={g:i for i,g in enumerate(present)}
            for name,panel in [("IFN",IFN_GENES),("DEX_UP",DEX_UP),("DEX_DOWN",DEX_DOWN)]:
                scores[name][start:stop]=vals[:,[ix[g] for g in panel if g in ix]].mean(axis=1)
        out=[]
        for parent in pd.unique(detailed):
            pm=detailed==parent
            if not pm.any(): continue
            for sid in pd.unique(sample[pm]):
                if publication_sample_key(sid) not in keep_sample_keys:
                    continue
                sm=pm & (sample==sid)
                t=sm & ifn; b=sm & ~ifn
                if t.any() and b.any():
                    out.append({"sample_id":sid,"parent_cell_type":parent,"IFN_cells":int(t.sum()),"nonIFN_cells":int(b.sum()),
                                "IFN_score_in_IFN_cells":float(scores["IFN"][t].mean()),"IFN_score_in_nonIFN_cells":float(scores["IFN"][b].mean()),
                                "IFN_score_difference":float(scores["IFN"][t].mean()-scores["IFN"][b].mean()),
                                "case_control":disease[np.flatnonzero(sm)[0]],"WHO_severity":severity[np.flatnonzero(sm)[0]]})
        ifn_scores = pd.DataFrame(out)
        ifn_scores.to_csv(OUT/"original_covid_ifn_marker_scores_by_sample.csv",index=False)
        isum=[]
        for parent,g in ifn_scores.groupby("parent_cell_type"):
            z=g.IFN_score_difference.dropna()
            if len(z):
                isum.append({"parent_cell_type":parent,"n_sample_parent_pairs":len(z),"median_IFN_score_difference":z.median(),
                             "positive_differences":int((z>0).sum()),"paired_wilcoxon_pvalue":wilcoxon(z).pvalue if np.any(z!=0) else np.nan})
        pd.DataFrame(isum).to_csv(OUT/"original_covid_ifn_marker_score_summary.csv",index=False)
        mono_out=[]
        for sid in pd.unique(sample[mono]):
            if publication_sample_key(sid) not in keep_sample_keys:
                continue
            sm=mono & (sample==sid); t=sm & dex; b=sm & ~dex
            if t.any() and b.any():
                i=np.flatnonzero(sm)[0]
                mono_out.append({"sample_id":sid,"patient_id":patient[i],"dex_candidate_cells":int(t.sum()),"other_monocytes":int(b.sum()),
                                 "DEX_up_score_candidate":float(scores["DEX_UP"][t].mean()),"DEX_up_score_background":float(scores["DEX_UP"][b].mean()),
                                 "DEX_down_score_candidate":float(scores["DEX_DOWN"][t].mean()),"DEX_down_score_background":float(scores["DEX_DOWN"][b].mean()),
                                 "case_control":disease[i],"WHO_severity":severity[i]})
        dex_scores=pd.DataFrame(mono_out)
        dex_scores.to_csv(OUT/"original_covid_dex_marker_scores_by_sample.csv",index=False)
        if len(dex_scores):
            up=dex_scores.DEX_up_score_candidate-dex_scores.DEX_up_score_background
            down=dex_scores.DEX_down_score_candidate-dex_scores.DEX_down_score_background
            pd.DataFrame([
                {"panel":"DEX_up", "n_sample_pairs":len(up),"median_candidate_score":dex_scores.DEX_up_score_candidate.median(),
                 "median_background_score":dex_scores.DEX_up_score_background.median(),"median_paired_difference":up.median(),
                 "positive_differences":int((up>0).sum()),"paired_wilcoxon_pvalue":wilcoxon(up).pvalue if np.any(up!=0) else np.nan},
                {"panel":"DEX_down", "n_sample_pairs":len(down),"median_candidate_score":dex_scores.DEX_down_score_candidate.median(),
                 "median_background_score":dex_scores.DEX_down_score_background.median(),"median_paired_difference":down.median(),
                 "positive_differences":int((down>0).sum()),"paired_wilcoxon_pvalue":wilcoxon(down).pvalue if np.any(down!=0) else np.nan},
            ]).to_csv(OUT/"original_covid_dex_marker_score_summary.csv",index=False)


if __name__=="__main__":
    print("Scoring original Ye/sepsis panels from sparse raw counts; this streams the full CSR matrix.",flush=True)
    sepsis_scores()
    print("Scoring original COVID/ESKD panels from the deposited logcounts layer.",flush=True)
    covid_scores()
    print("Original marker score tables written.",flush=True)
