"""Matched real-A fine-cluster composition and transcriptome comparison.

Uses the saved 60k pairwise integrations/30-NN transfer labels for previous
K4000, floor and new ceiling. Real B is the held-out positive control.
"""
import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
OUT = ROOT / "NewLora/analysis_20261102/three_way"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE / "code"))
import ilisi_lib as L

SEED = 20261024
SOURCES = ("realB", "prev4k", "floor", "ceiling")
clock = time.time()
def log(*msg): print(f"[{time.time()-clock:6.1f}s]", *msg, flush=True)

def h5col(f, path):
    g = f[path]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object)
        return cats[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v

def rows_from_csr(f, rows, columns, step=20000):
    g = f["X"]; ip = g["indptr"][:]; ncol = int(g.attrs["shape"][1]); blocks=[]
    for start in range(0, int(rows.max())+1, step):
        chosen = rows[(rows>=start)&(rows<start+step)]
        if len(chosen)==0: continue
        end = min(start+step,len(ip)-1); a,b=int(ip[start]),int(ip[end])
        block=sp.csr_matrix((g["data"][a:b].astype(np.float32),g["indices"][a:b],ip[start:end+1]-a),shape=(end-start,ncol))
        blocks.append(block[chosen-start][:,columns])
    return sp.vstack(blocks,format="csr")

def ceiling_counts(coh, genes):
    path=ROOT/"NewLora"/({"sepsis":"ye2026_sepsis_draw1_K4000_ceiling.h5ad","covid":"imperial_eskd_draw0_K4000_ceiling.h5ad"}[coh])
    with h5py.File(path) as f:
        sample=h5col(f,"obs/sample_id").astype(str)
        group=h5col(f,"obs/group").astype(str)
        flag=h5col(f,"obs/primary_analysis_eligible").astype(bool)
        keep=(~np.char.startswith(group,"ped")) if coh=="sepsis" else flag
        eligible=np.where(keep)[0]
        rng=np.random.default_rng(SEED)
        per_sample=int(np.ceil(60000/len(np.unique(sample[eligible]))))
        picked=np.concatenate([rng.choice(eligible[np.asarray(ix)],size=min(per_sample,len(ix)),replace=False)
                               for _,ix in pd.Series(eligible).groupby(sample[eligible]).groups.items()])
        rows=np.sort(rng.choice(picked,size=min(60000,len(picked)),replace=False))
        symbols=h5col(f,"var/gene_symbol").astype(str)
        positions=pd.Series(np.arange(len(symbols)),index=symbols)
        positions=positions[~positions.index.duplicated()]
        cols=positions.reindex(genes)
        assert cols.notna().all()
        X=rows_from_csr(f,rows,cols.to_numpy(dtype=int))
    labels=pd.read_csv(ROOT/"NewLora/analysis_20261102"/coh/"ceiling_cell_labels.csv.gz",index_col=0)
    assert len(labels)==len(rows)
    assert np.array_equal(labels["sample"].to_numpy(dtype=str),sample[rows])
    return X,labels["fine_tr"].astype(str).to_numpy(),labels["sample"].astype(str).to_numpy()

def cp10k(X):
    total=np.asarray(X.sum(axis=1)).ravel()
    return (sp.diags(1e4/np.maximum(total,1))@X).tocsr()

def corr(x,y):
    if len(x)<3 or np.std(x)==0 or np.std(y)==0:return float("nan")
    return float(np.corrcoef(x,y)[0,1])

fine_rows=[]; summary=[]; expr_rows=[]; depth_rows=[]
for coh in ("sepsis","covid"):
    log(coh,"load real A")
    R=ad.read_h5ad(BASE/"data/big"/f"{coh}_ref.h5ad")
    genes=np.array(R.var_names)
    rlab=R.obs.fine.astype(str).to_numpy()
    RX=R.X.tocsr().astype(np.float32)
    hv=L.hvg_seurat_v3([RX]); hv_mask=np.zeros(len(genes),bool);hv_mask[hv]=True
    Rcp=cp10k(RX)
    rmean=np.asarray(Rcp.mean(axis=0)).ravel()
    low_mask=pd.Series(rmean).rank(ascending=False,method="first").to_numpy()>2000
    rglobal=np.log1p(rmean)
    cats=np.unique(rlab)
    nr={t:int(np.sum(rlab==t)) for t in cats}
    rprofile={t:np.log1p(np.asarray(Rcp[rlab==t].mean(axis=0)).ravel()) for t in cats}
    rfrac={t:nr[t]/len(rlab) for t in cats}
    del Rcp,RX
    for source in SOURCES:
        log(coh,source,"load counts")
        if source=="ceiling":
            X,labels,samples=ceiling_counts(coh,genes)
        else:
            src_name="realU" if source=="realB" else source
            D=ad.read_h5ad(BASE/"data/big"/f"{coh}_{src_name}.h5ad")
            assert np.array_equal(np.array(D.var_names),genes)
            X=D.X.tocsr().astype(np.float32)
            transfer=pd.read_csv(BASE/"results/big"/f"{coh}_{src_name}_transfer.csv.gz",index_col=0)
            assert D.obs.index.equals(transfer.index)
            labels=transfer.fine_tr.astype(str).to_numpy()
            samples=D.obs["sample"].astype(str).to_numpy()
            del D,transfer
        umi=np.asarray(X.sum(axis=1)).ravel()
        detected=np.diff(X.indptr)
        depth_rows.append({"cohort":coh,"source":source,"n_cells":X.shape[0],
                           "median_UMI":float(np.median(umi)),"median_detected_genes":float(np.median(detected))})
        C=cp10k(X); del X
        global_profile=np.log1p(np.asarray(C.mean(axis=0)).ravel())
        ns={t:int(np.sum(labels==t)) for t in cats}
        source_frac={t:ns[t]/len(labels) for t in cats}
        retained=[]
        for t in cats:
            fr=rfrac[t]; fs=source_frac[t]; ratio=fs/fr if fr>0 else np.nan
            # Abundant real-A cluster with a matching, reasonably sized cluster.
            reproduced=(fr>=.005) and (fs>=.005) and (.5<=ratio<=2.0) and (ns[t]>=100)
            fine_rows.append({"cohort":coh,"source":source,"fine_cluster":t,
                              "n_real_A":nr[t],"n_source":ns[t],"fraction_real_A":fr,
                              "fraction_source":fs,"source_to_real_ratio":ratio,"reproduced":reproduced})
            if ns[t]<100 or nr[t]<100:continue
            profile=np.log1p(np.asarray(C[labels==t].mean(axis=0)).ravel())
            real=rprofile[t]
            residual=profile-global_profile; rres=rprofile[t]-rglobal
            top_real=set(np.argsort(rres)[-50:]);top_src=set(np.argsort(residual)[-50:])
            row={"cohort":coh,"source":source,"fine_cluster":t,"reproduced":reproduced,
                 "n_real_A":nr[t],"n_source":ns[t],"fraction_real_A":fr,"fraction_source":fs,
                 "pearson_all":corr(real,profile),"pearson_HVG":corr(real[hv_mask],profile[hv_mask]),
                 "pearson_low_expr":corr(real[low_mask],profile[low_mask]),
                 "identity_residual_r":corr(rres[hv_mask],residual[hv_mask]),
                 "marker_top50_jaccard":len(top_real&top_src)/len(top_real|top_src),
                 "rmse_log1p_cp10k":float(np.sqrt(np.mean((real-profile)**2)))}
            expr_rows.append(row)
            if reproduced:retained.append(row)
        summary.append({"cohort":coh,"source":source,"n_fine_total":len(cats),
                        "n_abundant_real":int(sum(v>=.005 for v in rfrac.values())),
                        "n_reproduced":len(retained),
                        "real_A_mass_reproduced":float(sum(x["fraction_real_A"] for x in retained)),
                        "source_mass_reproduced":float(sum(x["fraction_source"] for x in retained)),
                        "median_pearson_all":float(np.median([x["pearson_all"] for x in retained])) if retained else np.nan,
                        "median_pearson_HVG":float(np.median([x["pearson_HVG"] for x in retained])) if retained else np.nan,
                        "median_pearson_low_expr":float(np.median([x["pearson_low_expr"] for x in retained])) if retained else np.nan,
                        "median_identity_residual_r":float(np.median([x["identity_residual_r"] for x in retained])) if retained else np.nan,
                        "median_marker_top50_jaccard":float(np.median([x["marker_top50_jaccard"] for x in retained])) if retained else np.nan,
                        "median_rmse_log1p_cp10k":float(np.median([x["rmse_log1p_cp10k"] for x in retained])) if retained else np.nan})
        del C
        log(coh,source,"done",len(retained),"clusters")
    del R

F=pd.DataFrame(fine_rows);E=pd.DataFrame(expr_rows);S=pd.DataFrame(summary)
F.to_csv(OUT/"fine_cluster_reproduction.csv",index=False)
E.to_csv(OUT/"per_cluster_transcriptome_similarity.csv",index=False)
S.to_csv(OUT/"three_way_summary.csv",index=False)
pd.DataFrame(depth_rows).to_csv(OUT/"cell_depth_summary.csv",index=False)
comps=[]; commons=[]; pairs=[]
metrics=["pearson_all","pearson_HVG","pearson_low_expr","identity_residual_r","marker_top50_jaccard","rmse_log1p_cp10k"]
for coh in ("sepsis","covid"):
    for source in SOURCES:
        x=F[(F.cohort==coh)&(F.source==source)]
        p=x.fraction_real_A.to_numpy();q=x.fraction_source.to_numpy();m=(p+q)/2
        with np.errstate(divide="ignore",invalid="ignore"):
            jsd=float((np.where(p>0,p*np.log2(p/m),0).sum()+np.where(q>0,q*np.log2(q/m),0).sum())/2)
        comps.append({"cohort":coh,"source":source,"fine_TV":float(.5*np.abs(p-q).sum()),
                      "fine_JSD":jsd,"fine_fraction_r":corr(p,q)})
    ff=F[(F.cohort==coh)&F.reproduced]
    common=set.intersection(*[set(ff[ff.source==source].fine_cluster) for source in ("prev4k","floor","ceiling")])
    for source in SOURCES:
        ee=E[(E.cohort==coh)&(E.source==source)&E.fine_cluster.isin(common)]
        row={"cohort":coh,"source":source,"n_common_three":len(common)}
        for metric in metrics:row[f"median_{metric}"]=float(ee[metric].median())
        commons.append(row)
    for metric in metrics:
        w=E[(E.cohort==coh)&E.fine_cluster.isin(common)].pivot(index="fine_cluster",columns="source",values=metric)
        for t,v in w.iterrows():
            pairs.append({"cohort":coh,"fine_cluster":t,"metric":metric,
                          "ceiling_minus_floor":v.ceiling-v.floor,"floor_minus_prev4k":v.floor-v.prev4k,
                          "realB":v.realB,"prev4k":v.prev4k,"floor":v.floor,"ceiling":v.ceiling})
pd.DataFrame(comps).to_csv(OUT/"fine_composition_summary.csv",index=False)
pd.DataFrame(commons).to_csv(OUT/"common_cluster_transcriptome_summary.csv",index=False)
pd.DataFrame(pairs).to_csv(OUT/"common_cluster_paired_metrics.csv",index=False)
log("saved",OUT)
