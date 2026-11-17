"""Independent recompute (no code reused from run_and_plot_public_replication_v8.py) of the core public-cohort
statistics: per-patient acute % (mean over acute samples), Spearman vs peak WHO, Mann-Whitney death AUC,
and DerSimonian-Laird random-effects meta via statsmodels.combine_effects. Compares to data/*.csv."""
import os, warnings, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.meta_analysis import combine_effects
warnings.filterwarnings('ignore')
ROOT = os.environ.get('PBISC_ROOT', r'F:\PBISC')
P1 = os.path.join(ROOT, 'public_covid_who_hspc_nk_20261115')
LIU = os.path.join(ROOT, 'HSPC_CD56brightNK_COVID_analysis_summary_20261116', '10_CD56bright_public_replication', 'data', 'obs_liu2021.parquet')
OUT = os.path.join(ROOT, 'HSPC_CD56brightNK_COVID_analysis_summary_20261116', '16_Public_COVID_Replication_ver8', 'data')

def per_sample(df, sid, feats, den):
    X = pd.DataFrame({k: v.values for k, v in feats.items()}, index=df.index); X['_den'] = den.values; X['sid'] = df[sid].astype(str).values
    return X.groupby('sid').sum()

# Stephenson
s = pd.read_parquet(os.path.join(P1, 'obs_stephenson2021.parquet')); s = s[(s.Status.astype(str) == 'Covid') & (s.donor_id.astype(str) != 'CV0198')]
a = s.author_cell_type.astype(str)
S = per_sample(s, 'sample_id', {'NK_CD56bright': a.eq('NK_56hi'), 'HSPC': a.str.startswith('HSC'), 'T_CD8_CTL': a.eq('CD8.TE'), 'T_CD8_CM': a.eq('CD8.EM'),
               'B_ASC_total': a.isin(['Plasmablast', 'Plasma_cell_IgG', 'Plasma_cell_IgA', 'Plasma_cell_IgM']), 'CD16': a.eq('CD16_mono'), 'Mono': a.str.contains('mono')},
               ~a.isin(['RBC', 'Platelets', 'B_malignant']))
m = s.drop_duplicates('sample_id').set_index(s.drop_duplicates('sample_id').sample_id.astype(str))
S['donor'] = m.donor_id.astype(str); wcs = m.Worst_Clinical_Status.astype(str).str.strip()
S['who'] = wcs.map({'Asymptomatic': 1, 'Mild': 4, 'Moderate': 5, 'Severe': 6, 'Critical': 7, 'Death': 10})
S['death'] = wcs.eq('Death') | m.Outcome.astype(str).eq('Death'); S['acute'] = m.Collection_Day.astype(str).str.lstrip('D').astype(float) < 28
# COMBAT
c = pd.read_parquet(os.path.join(P1, 'obs_combat2022.parquet')); c = c[c.Source.astype(str).isin(['COVID_HCW_MILD', 'COVID_MILD', 'COVID_SEV', 'COVID_CRIT'])]
mi, ma = c.minor_subset.astype(str), c.major_subset.astype(str)
C = per_sample(c, 'scRNASeq_sample_ID', {'NK_CD56bright': mi.eq('NK.CD56hi.CD16int.XCL1/2'), 'HSPC': ma.eq('HSC'), 'T_CD8_CTL': mi.isin(['CD8.TEMRA', 'CD8.TEFF', 'CD8.TEM']),
               'T_CD8_CM': mi.isin(['CD8.TCM', 'CD8.TCM.CCL5']), 'B_ASC_total': ma.eq('PB'), 'CD16': ma.eq('ncMono'), 'Mono': ma.isin(['cMono', 'ncMono'])}, ~ma.isin(['PLT', 'RET', 'nan']))
m = c.drop_duplicates('scRNASeq_sample_ID'); m = m.set_index(m.scRNASeq_sample_ID.astype(str))
C['donor'] = m.donor_id.astype(str)
C['who'] = m.Outcome.astype(str).map({'Not hospitalized': 2, 'Hospitalized, no O2': 4, 'Hospitalized, O2': 5, 'Non-invasive ventilation': 6, 'Intubated, ventilated': 7, 'Death': 10})
C.loc[(m.Source.astype(str) == 'COVID_HCW_MILD') & ~m.Symptomatic.astype(bool), 'who'] = 1
C['death'] = m.Outcome.astype(str).eq('Death') | m.Death28.astype(bool); t = pd.to_numeric(m.TimeSinceOnset, errors='coerce'); C['acute'] = t.isna() | (t <= 28)
# van der Wijst
v = pd.read_parquet(os.path.join(P1, 'obs_vanderwijst2021.parquet')); v = v[v.COVID_status.astype(str) == 'Positive'].copy()
v['sid'] = v.donor_id.astype(str) + '_t' + v.timepoint.astype(str); c1, c2 = v.ct1.astype(str), v.ct2.astype(str)
V = per_sample(v, 'sid', {'NK_CD56bright': c2.eq('NK_CD56++'), 'HSPC': c1.eq('Progen'), 'T_CD8_CTL': c2.eq('T8_Mem'), 'B_ASC_total': c1.eq('PB'), 'CD16': c2.eq('ncM'),
               'Mono': c1.isin(['cM', 'ncM'])}, pd.Series(True, index=v.index))
m = v.drop_duplicates('sid').set_index('sid'); V['donor'] = m.donor_id.astype(str)
pk = pd.to_numeric(v.NIH_ordinal.astype(str), errors='coerce').groupby(v.donor_id.astype(str)).max()
V['death'] = m.death.astype(str).eq('yes'); V['who'] = V.donor.map(pk); V.loc[V.death & V.who.notna(), 'who'] = 10
V['acute'] = pd.to_numeric(m.days_from_start, errors='coerce').fillna(0) <= 14
# Liu
L = pd.read_parquet(LIU); L = L[L.timepoint.astype(str) != 'HC'].copy(); L['sid'] = L.donor_id.astype(str) + '_' + L.timepoint.astype(str); la = L.author_cell_type.astype(str)
Lq = per_sample(L, 'sid', {'NK_CD56bright': la.eq('NK_CD56hiCD16lo'), 'T_CD8_CTL': la.eq('CD8-positive, alpha-beta memory T cell'), 'B_ASC_total': la.eq('plasmablast'),
                'CD16': la.eq('non-classical monocyte'), 'Mono': la.str.contains('monocyte')}, ~la.isin(['platelet', 'granulocyte']))
both = L.groupby('sid').part.nunique(); Lq = Lq.loc[both[both == 2].index]
m = L.drop_duplicates('sid').set_index('sid').loc[Lq.index]
Lq['donor'] = m.donor_id.astype(str); Lq['death'] = m.outcome.astype(str).eq('deceased')
Lq['who'] = m.severity.astype(str).map({'Moderate': 4.5, 'Severe': 5.5, 'Critical': 7.0}); Lq.loc[Lq.death, 'who'] = 10.0
d = pd.to_numeric(m.days_since_onset.astype(str), errors='coerce'); Lq['acute'] = d.isna() | (d <= 28)

COH = {'Stephenson 2021': S, 'COMBAT 2022': C, 'van der Wijst 2021': V, 'Liu 2021': Lq}
rows = []
for name, T in COH.items():
    T = T[(T._den >= 500) & T.acute].copy()
    for f in ['NK_CD56bright', 'HSPC', 'T_CD8_CTL', 'T_CD8_CM', 'B_ASC_total']:
        if f in T: T[f + '_pct'] = 100 * T[f] / T._den
    T['CD16_mono_ofMono_pct'] = 100 * T.CD16 / T.Mono.replace(0, np.nan)
    pcols = [x for x in T.columns if x.endswith('_pct')]
    P = T.groupby('donor')[pcols].mean(); P['who'] = T.groupby('donor').who.first(); P['death'] = T.groupby('donor').death.max()
    for col in pcols:
        q = P.dropna(subset=['who', col]); r, p = spearmanr(q.who, q[col])
        z = P.dropna(subset=[col]); dd = z.death.astype(bool)
        u = mannwhitneyu(z[col][dd], z[col][~dd])
        rows.append(dict(cohort=name, feature=col[:-4], n_sev=len(q), rho=r, p_sev=p, n_all=len(z), deaths=int(dd.sum()), AUC=u.statistic / (dd.sum() * (~dd).sum()), p_death=u.pvalue))
R = pd.DataFrame(rows)
ref = pd.read_csv(os.path.join(OUT, 'public_cohort_level_results_v8.csv'))
J = R.merge(ref, on=['cohort', 'feature'], suffixes=('', '_ref'))
print('cohort-level rows compared:', len(J))
for a_, b_ in [('n_sev', 'n_sev'), ('rho', 'rho'), ('p_sev', 'p_sev'), ('AUC', 'AUC_death_all'), ('p_death', 'p_death_all'), ('deaths', 'n_death_all'), ('n_all', 'n_total_all')]:
    print(f'  max|{a_} - ref| = {np.nanmax(np.abs(J[a_] - J[b_ + ("" if b_ != a_ else "_ref")])):.2e}')
print('patients with WHO grade per cohort:', R.groupby('cohort').n_sev.max().to_dict(), '-> total', int(R.groupby('cohort').n_sev.max().sum()))
print('deaths per cohort:', R.groupby('cohort').deaths.max().to_dict(), '-> total', int(R.groupby('cohort').deaths.max().sum()))
# DL meta via statsmodels on Fisher z (SE = sqrt(1.06/(n-3)))
meta_ref = pd.read_csv(os.path.join(OUT, 'public_meta_analysis_v8.csv')).set_index('feature')
specs = {'NK_CD56bright': ['Stephenson 2021', 'COMBAT 2022', 'van der Wijst 2021'], 'HSPC': ['Stephenson 2021', 'COMBAT 2022', 'van der Wijst 2021'],
         'T_CD8_CTL': list(COH), 'CD16_mono_ofMono': list(COH), 'B_ASC_total': ['Stephenson 2021', 'COMBAT 2022', 'van der Wijst 2021'], 'T_CD8_CM': ['Stephenson 2021', 'COMBAT 2022']}
for f, cl in specs.items():
    q = R[(R.feature == f) & R.cohort.isin(cl)]
    res = combine_effects(np.arctanh(q.rho.values), 1.06 / (q.n_sev.values - 3), method_re='dl')
    sf = res.summary_frame()
    # statsmodels does not truncate the DL tau^2 at 0; when Q <= k-1 the DL estimate equals the fixed-effect estimate
    row = 'random effect' if res.q > len(q) - 1 else 'fixed effect'
    est = sf.loc[row, 'eff']; ci = sf.loc[row, ['ci_low', 'ci_upp']].values
    from scipy.stats import norm
    pz = 2 * norm.sf(abs(est / sf.loc[row, 'sd_eff']))
    mr = meta_ref.loc[f]
    print(f'  meta {f:17s} k={len(q)} rho={np.tanh(est):+.4f} [{np.tanh(ci[0]):+.3f},{np.tanh(ci[1]):+.3f}] P={pz:.3g} I2={100*max(res.i2,0):.1f}% | ref rho={mr.rho_meta:+.4f} P={mr.p_rho_meta:.3g} I2={mr.I2_rho:.1f}%')
