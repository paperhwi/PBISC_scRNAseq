"""Export frozen canonical25 results for this publication release.

Rendering needs only the bundled data; this optional exporter needs the original
canonical directory. It does not refit associations or rerun the generator.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--canonical', type=Path, required=True)
    args = parser.parse_args()
    src = args.canonical.resolve()
    data = ROOT / 'data'
    data.mkdir(parents=True, exist_ok=True)
    inputs = []

    def register(path):
        rel = str(path.relative_to(src)).replace('\\', '/')
        if not any(row['canonical_relative_path'] == rel for row in inputs):
            inputs.append({'canonical_relative_path': rel, 'bytes': path.stat().st_size,
                           'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})

    selected = {
        '': ['PRIMARY_CANDIDATES.csv'],
        'composition': ['all43_associations.csv', 'all43_fixed_score_bootstrap_CI.csv',
            'all43_gain_loss_and_support.csv', 'all43_source_heterogeneity.csv',
            'all43_source_stratified_permutation.csv', 'candidate_evidence_all_rescues.csv',
            'cohort_summary.csv', 'rarity_all43.csv', 'rare_threshold_sensitivity.csv'],
        'transcriptome': ['all14_primary_coverage_and_composition_rarity.csv',
            'coverage_all14_labels.csv', 'fixed_primary_candidates_with_rarity_and_composition.csv',
            'fixed_primary_candidate_CIs_50000.csv', 'fixed_primary_candidate_sensitivity_comparison.csv',
            'screen_counts_all_configurations.csv', 'screen_counts_overall.csv',
            'primary_all_gene_associations.csv.gz', 'all_gene_shared_patient_fidelity.csv.gz',
            'shared_unique_gene_universe.csv', 'raw_expression_QC.csv', 'duplicated_gene_symbols.csv',
            'historical_HSPC20_both_endpoints.csv'],
        'modules': ['HSPC_all_module_associations.csv', 'HSPC_all_modules_rescue_flags.csv'],
        'transcriptome_verification': ['HSPC_historical20_all_tests.csv',
            'KIAA0391_source_permutation_10M.csv',
            'KIAA0391_source_permutation_10M_Bonferroni_sensitivity.csv',
            'selected_candidates_source_permutation_1M.csv'],
        'verification': ['Bulk_KIAA0391_context.csv', 'Bulk_KIAA0391_geneID_mapping.csv'],
    }
    for folder, files in selected.items():
        for name in files:
            p = src / folder / name
            assert p.exists(), p
            register(p)
            target = data / folder / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, target)

    clinical_path = src / 'transcriptome/clinical_patient_table.csv'
    clinical = pd.read_csv(clinical_path).set_index('patient_id')
    register(clinical_path)
    # Stable release-local pseudonyms; no mapping or age/sex is published.
    identifiers = {}
    for arm, prefix in [('Real', 'R'), ('New', 'N')]:
        for i, pid in enumerate(sorted(clinical.index[clinical.arm.eq(arm)]), 1):
            identifiers[pid] = f'{prefix}{i:03}'
    plot = data / 'plotting'
    plot.mkdir(exist_ok=True)
    for source_name, dest_name in [('patient_raw_percent.csv', 'composition_patients.csv'),
                                   ('patient_fixed_INT_scores.csv', 'composition_scores.csv')]:
        p = src / 'composition' / source_name
        register(p)
        frame = pd.read_csv(p).drop(columns=['age', 'sex'])
        frame['patient_id'] = frame.patient_id.map(identifiers)
        assert frame.patient_id.notna().all()
        frame.to_csv(plot / dest_name, index=False)

    rows = []
    for celltype, genes in [('HSPC', ['KIAA0391']), ('T g/d', ['ADK', 'MXRA7'])]:
        slug = celltype.replace(' ', '_').replace('/', '_')
        for source, arm in [('Real', 'Real'), ('PBISC', 'New')]:
            p = src / 'transcriptome/patient_matrices' / f'{source}__{slug}.npz'
            register(p)
            z = np.load(p, allow_pickle=True)
            for i, pid in enumerate(z['patient_ids'].astype(str)):
                if pid not in clinical.index or clinical.loc[pid, 'arm'] != arm or z['ncells'][i] < 3:
                    continue
                for gene in genes:
                    j = list(z['genes'].astype(str)).index(gene)
                    rows.append(dict(patient_id=identifiers[pid], arm=arm, celltype=celltype, gene=gene,
                        ncells=int(z['ncells'][i]), nvisits=int(z['nvisits'][i]),
                        severity=clinical.loc[pid, 'severity'], death=int(clinical.loc[pid, 'death']),
                        mean_logCP10K=z['cell_mean_logCP10K'][i, j],
                        equal_visit_logCP10K=z['equal_visit_mean_logCP10K'][i, j]))
    patients = pd.DataFrame(rows)
    patients['score'] = np.nan
    for _, frame in patients.groupby(['arm', 'celltype', 'gene']):
        patients.loc[frame.index, 'score'] = stats.norm.ppf((stats.rankdata(frame.mean_logCP10K) - .5) / len(frame))
    patients.to_csv(plot / 'RNA_candidate_patients.csv', index=False)

    r = np.load(src / 'transcriptome/patient_matrices/Real__HSPC.npz', allow_pickle=True)
    p = np.load(src / 'transcriptome/patient_matrices/PBISC__HSPC.npz', allow_pickle=True)
    ri = {pid: i for i, pid in enumerate(r['patient_ids'].astype(str)) if r['ncells'][i] >= 3}
    pi = {pid: i for i, pid in enumerate(p['patient_ids'].astype(str)) if p['ncells'][i] >= 3}
    shared = sorted(set(ri) & set(pi))
    gr = list(r['genes']).index('KIAA0391')
    gp = list(p['genes']).index('KIAA0391')
    bridge = pd.DataFrame([dict(patient_id=identifiers[pid], Real_expression=r['cell_mean_logCP10K'][ri[pid], gr],
        PBISC_expression=p['cell_mean_logCP10K'][pi[pid], gp],
        ncellsReal=int(r['ncells'][ri[pid]]), ncellsPBISC=int(p['ncells'][pi[pid]])) for pid in shared])
    assert len(bridge) == 24
    bridge.to_csv(plot / 'KIAA0391_shared_bridge.csv', index=False)

    bp = src / 'verification/Bulk_KIAA0391_same29patient_values.csv'
    register(bp)
    b = pd.read_csv(bp)
    b['patient_id'] = b.patient_id.map(identifiers)
    assert len(b) == 29 and b.patient_id.notna().all()
    b.to_csv(plot / 'bulk_KIAA0391_patients.csv', index=False)

    docs = ROOT / 'docs'
    docs.mkdir(exist_ok=True)
    for name in ['Positive_Narrative_KO.md', 'Results_Draft_EN.md']:
        original = src / 'narrative_draft' / name
        register(original)
        text = original.read_text(encoding='utf-8')
        # Preserve the approved narrative snapshot while using release-local links.
        text = text.replace('../PRIMARY_CANDIDATES.csv', '../data/PRIMARY_CANDIDATES.csv')
        text = text.replace('../transcriptome/', '../data/transcriptome/')
        text = text.replace('../RESULTS_KO.md', 'Canonical_Results_KO.md')
        text = text.replace('../transcriptome_verification/KIAA0391_10M_INTERPRETATION_KO.md', 'KIAA0391_10M_Interpretation_KO.md')
        text = text.replace('../verification/Bulk_KIAA0391_context.md', 'Bulk_KIAA0391_Context.md')
        text = text.replace('NARRATIVE_QA.json', '../provenance/NARRATIVE_QA.json')
        (docs / name).write_text(text, encoding='utf-8')
    for a, b in [('RESULTS_KO.md', 'Canonical_Results_KO.md'),
                 ('transcriptome_verification/KIAA0391_10M_INTERPRETATION_KO.md', 'KIAA0391_10M_Interpretation_KO.md'),
                 ('verification/Bulk_KIAA0391_context.md', 'Bulk_KIAA0391_Context.md')]:
        register(src / a)
        shutil.copyfile(src / a, docs / b)
    provenance = ROOT / 'provenance'
    provenance.mkdir(exist_ok=True)
    for a, b in [('ANALYSIS_CONTRACT.json', 'ANALYSIS_CONTRACT.json'),
                 ('VALIDATION.json', 'CANONICAL_VALIDATION.json'),
                 ('SOFTWARE_VERSIONS.json', 'CANONICAL_SOFTWARE_VERSIONS.json'),
                 ('narrative_draft/NARRATIVE_QA.json', 'NARRATIVE_QA.json'),
                 ('transcriptome_verification/KIAA0391_10M_validation.json', 'KIAA0391_10M_validation.json'),
                 ('transcriptome_verification/KIAA0391_10M_RUN_MANIFEST.json', 'KIAA0391_10M_RUN_MANIFEST.json')]:
        register(src / a)
        shutil.copyfile(src / a, provenance / b)
    (provenance / 'CANONICAL_INPUTS.json').write_text(json.dumps({
        'canonical_release': src.name, 'inputs': inputs,
        'bundled_canonical_tables': [str(Path(folder)/name).replace('\\','/') for folder,files in selected.items() for name in files],
        'patient_data': 'Release-local pseudonyms; original patient IDs, age, sex, and mapping not published',
        'scope': 'Frozen results and plotting vectors; original counts/generator not bundled',
    }, indent=2), encoding='utf-8')
    print(json.dumps({'tables': sum(len(v) for v in selected.values()), 'patient_rows': len(patients),
                      'bridge_patients': len(bridge), 'registered_inputs': len(inputs)}))

if __name__ == '__main__':
    main()
