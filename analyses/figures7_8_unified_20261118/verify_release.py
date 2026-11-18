"""Independent numeric checks against bundled public plotting/result tables.

This verifies frozen results and figure inputs, not the original generator or
all raw expression values. Original raw-data audits are recorded separately.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'

def bh(values):
    p = np.asarray(values, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1.0)
    return result

def close(a, b, label, atol=2e-12):
    np.testing.assert_allclose(np.asarray(a, dtype=float), np.asarray(b, dtype=float),
                               rtol=2e-11, atol=atol, err_msg=label)

def scalar_tests(frame, column):
    x = frame[column].to_numpy(float)
    if np.ptp(x) == 0:
        return [0.0, 1.0, .5, 1.0]
    rho, psev = stats.spearmanr(x, frame.severity)
    d, s = x[frame.death.eq(1)], x[frame.death.eq(0)]
    u, pdeath = stats.mannwhitneyu(d, s, alternative='two-sided', method='auto')
    return [rho, psev, u/(len(d)*len(s)), pdeath]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--figures', action='store_true', help='Also require all seven rendered figures')
    args = parser.parse_args()
    checks = {}
    original = json.loads((ROOT/'provenance/CANONICAL_INPUTS.json').read_text(encoding='utf-8'))
    required = set(original['bundled_canonical_tables'])
    assert len(required) == 31, 'This release requires exactly31 canonical source tables'
    for rel in required:
        assert (DATA/rel).is_file(), f'Missing canonical source table: {rel}'
    matched = 0
    for item in original['inputs']:
        path = DATA/item['canonical_relative_path']
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], str(path)
            matched += 1
    assert matched == len(required), (matched, len(required))
    checks['canonical_table_hashes'] = matched

    raw = pd.read_csv(DATA/'plotting/composition_patients.csv')
    scores = pd.read_csv(DATA/'plotting/composition_scores.csv')
    meta = ['patient_id','arm','severity','death','draws']
    features = [c for c in raw if c not in meta]
    assert len(raw) == 66 and len(features) == 43 and raw.death.sum() == 13
    assert raw.groupby('arm').size().to_dict() == {'New':29,'Real':37}
    assert raw.groupby('arm').death.sum().to_dict() == {'New':6,'Real':7}
    assert raw.draws.sum() == 217 and raw.patient_id.is_unique
    assert raw.patient_id.str.fullmatch(r'[RN]\d{3}').all()
    assert not {'age','sex'}.intersection(raw.columns)
    close(raw[features].sum(axis=1), 100.0, 'leukocyte closure')
    pd.testing.assert_frame_equal(raw[meta], scores[meta])
    for arm, d in raw.groupby('arm'):
        for feature in features:
            expected = stats.norm.ppf((stats.rankdata(d[feature])-.5)/len(d))
            close(expected, scores.loc[d.index,feature], f'composition source INT {arm}/{feature}')
    comp = pd.read_csv(DATA/'composition/all43_associations.csv')
    for arm in ['Real','New','Combined']:
        frame = scores if arm == 'Combined' else scores[scores.arm.eq(arm)]
        t = comp[comp.arm.eq(arm)].set_index('feature').loc[features]
        observed = np.array([scalar_tests(frame, f) for f in features])
        close(observed, t[['rho','p_severity','AUC','p_death']], f'43 scalar associations {arm}')
        for endpoint in ['severity','death']:
            close(bh(t['p_'+endpoint]),t['q43_'+endpoint],f'BH43 {arm}/{endpoint}')
        joint = bh(np.r_[t.p_severity,t.p_death])
        close(joint[:43],t.q86_severity,f'BH86 severity {arm}')
        close(joint[43:],t.q86_death,f'BH86 death {arm}')
    checks['composition'] = '129 scalar feature records; source INT, BH43, jointBH86, counts and closure pass'

    full = pd.read_csv(DATA/'transcriptome/primary_all_gene_associations.csv.gz')
    assert len(full) == 155574 and full.gene.nunique() == 5762 and full.celltype.nunique() == 9
    comparisons = 0
    for cohort, t in full.groupby('cohort', sort=False):
        assert len(t) == 51858
        for celltype, f in t.groupby('celltype', sort=False):
            assert len(f) == 5762
            for endpoint in ['severity','death']:
                close(bh(f['p_'+endpoint]),f['q_gene_'+endpoint],f'RNA BH5762 {cohort}/{celltype}/{endpoint}')
                comparisons += 1
            joint = bh(np.r_[f.p_severity,f.p_death])
            close(joint[:len(f)],f.q_joint_gene_severity,'RNA joint gene severity')
            close(joint[len(f):],f.q_joint_gene_death,'RNA joint gene death')
            comparisons += 2
        for endpoint in ['severity','death']:
            close(bh(t['p_'+endpoint]),t['q_global_'+endpoint],f'RNA global {cohort}/{endpoint}')
            comparisons += 1
        joint = bh(np.r_[t.p_severity,t.p_death])
        close(joint[:len(t)],t.q_global_joint_severity,'RNA global joint severity')
        close(joint[len(t):],t.q_global_joint_death,'RNA global joint death')
        comparisons += 2
    checks['RNA_BH'] = {'records':len(full),'column_comparisons':comparisons,'status':'PASS'}

    patients = pd.read_csv(DATA/'plotting/RNA_candidate_patients.csv')
    assert patients.patient_id.str.fullmatch(r'[RN]\d{3}').all()
    for (celltype,gene), group in patients.groupby(['celltype','gene']):
        for arm,d in group.groupby('arm'):
            expected = stats.norm.ppf((stats.rankdata(d.mean_logCP10K)-.5)/len(d))
            close(expected,d.score,'RNA source INT')
        for cohort in ['Real','New','Combined']:
            d = group if cohort=='Combined' else group[group.arm.eq(cohort)]
            t = full[full.celltype.eq(celltype)&full.gene.eq(gene)&full.cohort.eq(cohort)].iloc[0]
            assert len(d)==t.n and d.death.sum()==t.deaths
            close(scalar_tests(d,'score'),t[['rho','p_severity','AUC','p_death']],'RNA candidate scalar tests')
    bridge = pd.read_csv(DATA/'plotting/KIAA0391_shared_bridge.csv')
    assert len(bridge)==24
    rho,p = stats.spearmanr(bridge.Real_expression,bridge.PBISC_expression)
    cand = pd.read_csv(DATA/'transcriptome/fixed_primary_candidates_with_rarity_and_composition.csv')
    k = cand[cand.celltype.eq('HSPC')&cand.gene.eq('KIAA0391')].iloc[0]
    close([rho,p],[k.bridge_rho,k.bridge_P],'shared bridge')
    checks['RNA_candidate_plotting_vectors'] = 'nine scalar cohort/gene checks, source INT and shared24 bridge pass'

    modules = pd.read_csv(DATA/'modules/HSPC_all_module_associations.csv')
    unique = modules[modules.module.eq(modules.representative)]
    for _, f in unique.groupby(['configuration','cohort']):
        assert len(f)==16
        close(bh(f.p_severity),f.q_unique16_severity,'module BH16 severity')
        close(bh(f.p_death),f.q_unique16_death,'module BH16 death')
        joint=bh(np.r_[f.p_severity,f.p_death])
        close(joint[:16],f.q_joint32_severity,'module joint32 severity')
        close(joint[16:],f.q_joint32_death,'module joint32 death')
    checks['modules'] = '16 unique hypotheses per configuration/cohort; endpoint and joint correction pass'

    perm = pd.read_csv(DATA/'transcriptome_verification/KIAA0391_source_permutation_10M.csv')
    assert perm.B.eq(9999999).all()
    close((perm.exceedances+1)/(perm.B+1),perm.P_plus1,'permutation plus-one P')
    bf = pd.read_csv(DATA/'transcriptome_verification/KIAA0391_source_permutation_10M_Bonferroni_sensitivity.csv')
    close(np.minimum(bf.m*bf.P_plus1,1),bf.Bonferroni_plus1,'full-family Bonferroni')
    for cohort in ['New','Combined']:
        t=bf[bf.cohort.eq(cohort)&bf.m.eq(11524)].iloc[0]
        assert t.Bonferroni_tail_CP95_high < .05
    checks['10M_precision'] = 'plus-one formula, full-family factors and two-endpoint HSPC MCupper threshold pass'

    if args.figures:
        names=['Figure_7','Figure_8']+[f'Supplementary_Figure_S{i}' for i in range(1,6)]
        for name in names:
            for ext in ['pdf','svg','png','artists.json']:
                assert (ROOT/'figures'/f'{name}.{ext}').stat().st_size > 100
            artists=json.loads((ROOT/'figures'/f'{name}.artists.json').read_text(encoding='utf-8'))
            assert abs(artists['width_mm']-178)<.01
            assert min(t['fontsize_pt'] for t in artists['texts']) >= 7.99
            tree=ET.parse(ROOT/'figures'/f'{name}.svg')
            assert tree.getroot().attrib.get('viewBox')
        checks['rendered_artifact_structure']='seven figures; PDF/SVG/PNG and artists records;178mmwidth;minimum8pt'

    checks['scope']='Independent checks of frozen tables/plotting vectors, not a new raw-data or generator validation'
    checks['status']='PASS'
    (ROOT/'provenance/RELEASE_VALIDATION.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':
    main()
