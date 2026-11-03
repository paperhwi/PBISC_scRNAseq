"""Read-only arithmetic and provenance audit of the New LoRA result tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
checks = []

def record(section, check, passed, detail):
    checks.append(dict(section=section, check=check, passed=bool(passed), detail=str(detail)))

def close(a, b, atol=1e-9):
    return np.allclose(np.asarray(a, float), np.asarray(b, float), atol=atol, rtol=1e-7, equal_nan=True)

def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    n = ok.sum()
    if not n:
        return q
    ix = np.argsort(p[ok]); sorted_p = p[ok][ix]
    sorted_q = np.minimum.accumulate((sorted_p*n/np.arange(1,n+1))[::-1])[::-1]
    z = np.empty(n); z[ix] = np.minimum(sorted_q,1)
    q[ok] = z
    return q

for cohort in ('sepsis','covid'):
    s = json.loads((BASE/cohort/'summary.json').read_text())
    lab = pd.read_csv(BASE/cohort/'ceiling_cell_labels.csv.gz', index_col=0)
    record(cohort, '60k label rows, unique IDs', len(lab)==60000 and lab.index.is_unique, len(lab))
    six = lab.six_tr.value_counts(normalize=True)
    record(cohort, 'six-type cell fractions', all(close(s['six_cell_fraction_ceiling'][k],v) for k,v in six.items()), six.to_dict())
    cl = pd.read_csv(BASE/cohort/'joint_leiden_clusters.csv')
    nreal = cl.n_real_A.sum(); ngen = cl.n_ceiling.sum()
    cl_ok = (nreal==s['n_ref'] and ngen==s['n_ceiling'] and
             close(cl.fraction_real_A, cl.n_real_A/nreal) and
             close(cl.fraction_ceiling, cl.n_ceiling/ngen) and
             close(cl.delta_pp, 100*(cl.fraction_ceiling-cl.fraction_real_A)))
    record(cohort, 'joint cluster counts and fractions', cl_ok, f'{nreal}, {ngen}, {len(cl)}')
    r = np.corrcoef(cl.fraction_real_A,cl.fraction_ceiling)[0,1]
    p = cl.fraction_real_A.to_numpy(); q = cl.fraction_ceiling.to_numpy(); m=(p+q)/2
    with np.errstate(divide='ignore',invalid='ignore'):
        jsd=(np.nansum(np.where(p>0,p*np.log2(p/m),0))+np.nansum(np.where(q>0,q*np.log2(q/m),0)))/2
    mixed=cl.loc[cl.source_share_normalized.between(.2,.8),'fraction_ceiling'].sum()
    record(cohort, 'joint cluster summary metrics', close([r,jsd,mixed,len(cl)],[s['cluster_fraction_pearson'],s['cluster_fraction_JSD'],s['ceiling_fraction_in_mixed_clusters'],s['n_leiden_clusters']]), f'r={r:.6f}; JSD={jsd:.6f}; mixed={mixed:.6f}')
    comp = pd.read_csv(BASE/cohort/'composition_by_sample.csv')
    summ = pd.read_csv(BASE/cohort/'composition_summary.csv')
    calc = comp.groupby(['level','source','cluster']).fraction.agg(['mean','std','count']).reset_index()
    c=summ.merge(calc,on=['level','source','cluster'],suffixes=('_out','_calc'),validate='one_to_one')
    good=len(c)==len(summ)==len(calc) and all(close(c[f'{col}_out'],c[f'{col}_calc']) for col in ('mean','std','count'))
    record(cohort,'sample-weighted composition summary',good,f'{len(comp)} rows, {len(summ)} summaries')

score=pd.read_csv(BASE/'three_way/three_way_scorecard.csv')
fine=pd.read_csv(BASE/'three_way/fine_cluster_reproduction.csv')
prof=pd.read_csv(BASE/'three_way/per_cluster_transcriptome_similarity.csv')
common=pd.read_csv(BASE/'three_way/common_cluster_transcriptome_summary.csv')
pair=pd.read_csv(BASE/'three_way/common_cluster_paired_metrics.csv')
for (cohort,source), g in fine.groupby(['cohort','source']):
    sc=score[(score.cohort==cohort)&(score.source==source)].iloc[0]
    abundant=g.fraction_real_A>=.005
    record('three_way',f'{cohort}/{source} fine reproduction',
           int((g.reproduced&abundant).sum())==sc.n_reproduced and int(abundant.sum())==sc.n_abundant_real,
           f"reproduced={int((g.reproduced&abundant).sum())}; abundant={int(abundant.sum())}")
    pg=prof[(prof.cohort==cohort)&(prof.source==source)]
    record('three_way',f'{cohort}/{source} expression profile count',len(pg)>0 and pg.fine_cluster.is_unique,len(pg))
    criterion=(g.fraction_real_A>=.005)&(g.fraction_source>=.005)&g.source_to_real_ratio.between(.5,2)&(g.n_source>=100)
    record('three_way',f'{cohort}/{source} reproduction criterion',np.array_equal(g.reproduced.to_numpy(bool),criterion.to_numpy(bool)),len(g))
    names=set.intersection(*[set(fine[(fine.cohort==cohort)&(fine.source==src)&fine.reproduced].fine_cluster) for src in ('prev4k','floor','ceiling')])
    row=common[(common.cohort==cohort)&(common.source==source)].iloc[0]
    p=pg[pg.fine_cluster.isin(names)]
    metrics=['pearson_all','pearson_HVG','pearson_low_expr','identity_residual_r','marker_top50_jaccard','rmse_log1p_cp10k']
    okay=(len(p)==row.n_common_three and all(close(p[col].median(),row['median_'+col]) and close(p[col].median(),sc['median_'+col]) for col in metrics))
    record('three_way',f'{cohort}/{source} common profile medians',okay,len(names))
record('three_way','paired cluster differences',
       close(pair.ceiling_minus_floor,pair.ceiling-pair.floor) and close(pair.floor_minus_prev4k,pair.floor-pair.prev4k),len(pair))

deg_dir=BASE/'disease_deg_matched/deg_tables'
old_dir=BASE.parents[1]/'bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/results/deg_tables'
ds=pd.read_csv(BASE/'disease_deg_matched/matched_deg_summary.csv')
for p in sorted(deg_dir.glob('*.csv.gz')):
    tab=pd.read_csv(p)
    stem=p.name.removesuffix('.csv.gz')
    n=tab.gene.nunique()==len(tab)
    if stem.endswith('_real'):
        original=pd.read_csv(old_dir/p.name).set_index('gene').reindex(tab.gene).reset_index()
        # Archived real q was calculated before the measured-gene filter.
        qok=close(tab.q,original.q,atol=1e-8)
    else:
        qok=close(tab.q,bh(tab.p),atol=1e-8)
    flag=np.asarray((tab.q<.05)&(tab.logFC.abs()>=1))
    degok=np.array_equal(tab.deg.to_numpy(bool),flag)
    fcok=close(tab.logFC,tab.mean_case-tab.mean_ctrl,atol=2e-5)
    record('disease_deg',stem,n and qok and degok and fcok,
           f'n={len(tab)}, q={qok}, flag={degok}, logFC={fcok}, DEG={flag.sum()}')
for row in ds.itertuples(index=False):
    a=pd.read_csv(old_dir/f'{row.cohort}_{row.celltype.replace(" ","")}_real.csv.gz').set_index('gene')
    b=pd.read_csv(deg_dir/f'{row.cohort}_{row.celltype.replace(" ","")}_{row.source}.csv.gz')
    a=a.reindex(b.gene).reset_index()
    overlap=int((a.deg&b.deg).sum())
    ok=(a.gene.equals(b.gene) and len(a)==row.n_genes and
        int(a.deg.sum())==row.n_DEG_real and int(b.deg.sum())==row.n_DEG_gen and
        overlap==row.n_DEG_overlap and close(row.recall,overlap/a.deg.sum()) and
        close(row.precision,overlap/b.deg.sum()))
    record('disease_deg',f'{row.cohort}/{row.celltype}/{row.source} recall',ok,f'{overlap}/{a.deg.sum()}, {overlap}/{b.deg.sum()}')

rare=pd.read_csv(BASE/'rare_populations/rare_unit_shares_real_vs_ceiling.csv')
record('rare_populations','patient percentages',
       close(rare.pct_all,100*rare.n_subset/rare.n_cells) and
       close(rare.pct_parent,np.where(rare.n_parent>=20,100*rare.n_subset/rare.n_parent,np.nan)),
       f'{len(rare)} rows')
stats=pd.read_csv(BASE/'rare_populations/rare_stats_real_vs_ceiling.csv')
for row in stats.itertuples(index=False):
    rr=rare[(rare.cohort==row.cohort)&(rare.subset==row.subset)&(rare.group==row.group)&(rare.source=='real')]
    cc=rare[(rare.cohort==row.cohort)&(rare.subset==row.subset)&(rare.group==row.group)&(rare.source=='ceiling')]
    ok=(len(rr)==row.n_real and len(cc)==row.n_ceiling and
        close(rr.pct_all.mean(),row.mean_real) and close(cc.pct_all.mean(),row.mean_ceiling) and
        close(rr.pct_all.std(),row.sd_real) and close(cc.pct_all.std(),row.sd_ceiling) and
        close(scipy_stats.ttest_ind(rr.pct_all,cc.pct_all,equal_var=False).pvalue,row.welch_p,atol=1e-8))
    record('rare_populations',f'{row.cohort}/{row.subset}/{row.group}',ok,f'n={len(rr)}/{len(cc)}')
for (cohort,group), x in stats.groupby(['cohort','group']):
    record('rare_populations',f'{cohort}/{group} Welch BH family',close(x.welch_q,bh(x.welch_p),atol=1e-8),len(x))
th=pd.read_csv(BASE/'rare_populations/real_T_CD4_c07_MKI67_thinned_496_seed11.csv')
record('rare_populations','496-cell thinning lattice',
       len(th)==40 and th.n_thinned.eq(496).all() and
       (th.k_thinned<=th.n_subset).all() and close(th.pct_thinned,100*th.k_thinned/496),
       f"positive distinct {th.loc[th.pct_all>0,'pct_all'].nunique()} -> {th.loc[th.pct_thinned>0,'pct_thinned'].nunique()}")

rp=pd.read_csv(BASE/'figure_realism_rare/rare_parent_share_with_new_lora.csv')
record('rare_parent','pooled parent percentages',close(rp.pct_of_parent,100*rp.n_target/rp.n_parent),f'{len(rp)} rows')
for filename,col in [('c2st_all_sources.csv','c2st_auc'),('density_coverage_all_sources.csv','density')]:
    x=pd.read_csv(BASE/'figure_realism_rare'/filename)
    good=(x[col].between(0,1).all() and (x.coverage.between(0,1).all() if 'coverage' in x else True))
    record('realism',filename,good,f'{len(x)} rows')

feas=pd.read_csv(BASE/'rare_disease_deg/feasibility_by_condition.csv')
cnt=pd.read_csv(BASE/'rare_disease_deg/rare_patient_counts.csv.gz')
for row in feas.itertuples(index=False):
    x=cnt[(cnt.cohort==row.cohort)&(cnt.source==row.source)&(cnt.subset==row.subset)&(cnt.condition==row.condition)]
    good=len(x)>0 and int(x.n_subset.sum())==row.n_cells_total and all(int((x.n_subset>=k).sum())==getattr(row,f'n_patients_ge{k}') for k in (1,5,10,20))
    record('rare_deg',f'{row.cohort}/{row.source}/{row.subset}/{row.condition} feasibility',good,len(x))
for p in (BASE/'rare_disease_deg').glob('mono_c01_*_vs_HC_*.csv.gz'):
    x=pd.read_csv(p)
    record('rare_deg',p.name,close(x.q,bh(x.p),atol=1e-8) and
           np.array_equal(x.deg.to_numpy(bool),np.asarray((x.q<.05)&(x.logFC.abs()>=1))),
           f'{len(x)} genes, {x.deg.sum()} DEG')

df=pd.DataFrame(checks)
df.to_csv(OUT/'arithmetic_checks.csv',index=False)
print(df.groupby('section').passed.agg(['count','sum']).to_string())
print('FAILED:',len(df)-df.passed.sum())
if (~df.passed).any():
    print(df.loc[~df.passed].to_string(index=False))
