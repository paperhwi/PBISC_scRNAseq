"""Write a readable comparison and provenance for corrected_deg.py results."""
from pathlib import Path
import hashlib, json, sys
import numpy as np
import pandas as pd
import scipy

root=Path('F:/PBISC')
review=root/'project_review_20261028'
out=review/'corrected_deg'
a=pd.read_csv(out/'corrected_summary.csv')
old=pd.read_csv(root/'bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/results/deg_summary.csv')
old=old[old.source=='floor']
joined=a.merge(old[['cohort','celltype','logFC_r','recall','precision','n_DEG_real','n_DEG_gen']],on=['cohort','celltype'],suffixes=('','_old'))
sepsis_check=joined[joined.cohort=='sepsis']
assert len(sepsis_check)==5
assert np.allclose(sepsis_check.logfc_r,sepsis_check.logFC_r)
assert np.allclose(sepsis_check.recall,sepsis_check.recall_old)
assert (sepsis_check.n_real_deg==sepsis_check.n_DEG_real).all()
assert (sepsis_check.n_floor_deg==sepsis_check.n_DEG_gen).all()
outside=((a.logfc_r<a.logfc_r_lo)|(a.logfc_r>a.logfc_r_hi)).sum()

def f(x): return '—' if pd.isna(x) else f'{x:.3f}'
lines=[
    '# PBISC DEG 긴급 재계산 결과 — 2026-10-28',
    '',
    '10/27의 bootstrap 구현 오류와 COVID 환자 의존성/roster 차이를 재계산했다. 기존 CSV와 발표자료는 보존했으며, 이 결과는 별도 수정 분석이다.',
    '',
    '## 설계',
    '',
    '- 입력: 10/27과 같은 sample×six-type pseudobulk raw-count 합과 log2(CPM+1), 같은 공통 측정 유전자. FDR<0.05, |log2FC|≥1.',
    '- Sepsis: 성인 sepsis 대 HC, Welch DEG. 환자(샘플)를 두 군 안에서 각각 복원추출해 200회 반복한다.',
    '- COVID: 환자×상태 평균을 분석 단위로 두고, 각 환자의 양성·음성 행을 함께 묶은 OLS cluster-robust(CR1) DEG를 계산한다. bootstrap은 환자를 복원추출해 모든 상태를 함께 반복한다.',
    '- 주분석은 출판 명단 187샘플 중 실제 양성/음성 대비가 가능한 177샘플이다. Recovery 10샘플은 이 대비에서 제외된다. 198 local 중 양성/음성 187샘플을 민감도로 계산한다. 생성 COVID는 기존과 동일하게 학습 관련 5개 샘플을 제외했다.',
    f'- 표의 범위는 200 patient-cluster bootstrap 점추정치 분포의 2.5/97.5 분위수다. logFC r에서 {outside}/{len(a)}개 행은 원래 점추정치가 이 범위 밖이다. 편향과 작은 Sepsis HC군(7명)을 고려하면 이를 확정적인 95% 신뢰구간으로 인용하지 않는다.',
    '- 실제와 생성 데이터는 같은 샘플의 일대일 쌍이 아니다. COVID 주분석에서 실제 177개 샘플과 생성 157개 샘플은 sample ID 교집합이 0이고 환자 ID 교집합은 26명이다. 두 자료의 환자를 각각 재표집했으며, 결과는 환자별 paired 재현도가 아닌 집단 질병 대비의 일치도다.',
    '',
    '## 새 결과',
    '',
    '| 코호트 | 실제 명단 | 세포형 | 실제 환자 | 실제 양성/대조 unit | 생성 환자 | logFC r [bootstrap 분위수] | DEG 회수율 [bootstrap 분위수] | 실제/생성 DEG |',
    '|---|---|---|---:|---:|---:|---:|---:|---:|',
]
for _,r in a.iterrows():
    lines.append(f"| {r.cohort} | {r.roster} | {r.celltype} | {int(r.n_real_patients)} | {int(r.n_real_case_units)}/{int(r.n_real_ctrl_units)} | {int(r.n_floor_patients)} | {f(r.logfc_r)} [{f(r.logfc_r_lo)}, {f(r.logfc_r_hi)}] | {f(r.recall)} [{f(r.recall_lo)}, {f(r.recall_hi)}] | {int(r.n_real_deg)}/{int(r.n_floor_deg)} |")
lines += [
    '',
    'Sepsis 점추정치와 DEG 개수는 10/27 원 결과와 정확히 일치한다. bootstrap 범위는 환자 중복 추출이 유지되어 바뀌었다. COVID는 발표 명단과 환자 내 상관을 반영하면서 비교 대상과 DEG 호출이 달라졌다.',
    '',
    '## 10/27 결과와의 차이',
    '',
    '| 코호트 | 명단 | 세포형 | 이전 r | 새 r | 이전 recall | 새 recall | 이전 실제/생성 DEG | 새 실제/생성 DEG |',
    '|---|---|---|---:|---:|---:|---:|---:|---:|',
]
for _,r in joined.iterrows():
    if r.cohort=='covid' and r.roster=='all_local_198': continue
    lines.append(f"| {r.cohort} | {r.roster} | {r.celltype} | {f(r.logFC_r)} | {f(r.logfc_r)} | {f(r.recall_old)} | {f(r.recall)} | {int(r.n_DEG_real)}/{int(r.n_DEG_gen)} | {int(r.n_real_deg)}/{int(r.n_floor_deg)} |")
lines += [
    '',
    '이전 COVID 결과는 198 local 자료의 양성/음성 부분집합에서 patient×condition 평균을 Welch 검정한 값이다. 새 주분석과는 검정법·명단이 모두 다르므로 차이를 한 원인으로 해석하지 않는다. `corrected_summary.csv`에 198 local 민감도도 있다.',
    '',
    '## 보고 규칙과 남은 제한',
    '',
    '- `sign_all`은 공통 측정 유전자 전체 중 실제와 생성의 logFC 부호가 같은 비율이다. `sign_real_deg`는 **실제에서 DEG로 판정된 유전자**를 분모로 쓴다. 둘을 같은 sign agreement로 표시하지 않는다.',
    '- `recall` 분모는 실제 DEG 수, `precision` 분모는 생성 DEG 수다. 분모 0은 결측으로 기록한다. 수정 결과의 sign·top-100 지표는 [`METRIC_DEFINITIONS.md`](METRIC_DEFINITIONS.md)와 [`corrected_gene_metrics.csv`](corrected_gene_metrics.csv)에 분모를 고정했다. 이전 `top100_overlap`은 q<0.05만 사용한 교집합 개수라 새 정밀도와 다르다.',
    '- 기존 상한선 표의 Sepsis T 64%와 슬라이드의 72%는 서로 다른 floor 분자다. 64%는 matched-n floor, 72%는 full-n floor이다. 두 방식의 수치가 섞인 ‘달성률 범위’를 원고의 확정 결론으로 쓰지 않는다.',
    '- COVID 환자 내 상관에 대한 CR1은 같은 환자의 반복 관측을 고려한다. 시점·중증도·wave의 교란과 동일 환자 전후 비교라는 인과 주장을 해결하지는 않는다.',
    '- 반복 200회 분위수는 특히 작은 DEG 집합에서 불안정할 수 있다. 제출 전 표본 단위·명단·DEG 규칙을 확정하고 bootstrap 편향/구간 방법 및 반복 수의 민감도를 추가 점검한다. 현재 범위는 결과 검토용이다.',
    '',
    '## 확인 및 재실행',
    '',
    '- [`corrected_deg.py`](../corrected_deg.py): 환자 중복 보존, 두 상태 동시 표집, CR1 DEG, Welch 및 bootstrap.',
    '- Sepsis 전체 데이터 점추정치·DEG 개수는 기존 CSV와 정확히 일치했고, 보고서 생성 시 이를 자동 검증한다. 작은 장난감 자료에서 환자 중복 보존과 CR1 계수를 확인했다. 10/27 원본 CSV를 변경하지 않았다.',
    '- 명령: `PYTHONPATH=F:/PBISC/.analysis_runtime python F:/PBISC/project_review_20261028/corrected_deg.py` (Windows PowerShell에서는 `$env:PYTHONPATH=...`).',
    '',
    '절차 참고: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM. *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. 2026. [arXiv:2609.00065](https://arxiv.org/abs/2609.00065). [SciPy bootstrap](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html), [statsmodels cluster covariance](https://www.statsmodels.org/stable/generated/statsmodels.stats.sandwich_covariance.cov_cluster.html).',
]
(out/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

inputs=[review/'corrected_deg.py',review/'corrected_gene_metrics.py',review/'summarize_corrected.py',root/'paper_claim_analysis/covid_sample_patient_id_reconciliation.csv']
pb=root/'bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021/data/pb'
inputs += [pb/f'{source}_{coh}.npz' for source in ('real','floor') for coh in ('sepsis','covid')]
inputs += [pb/f'{source}_{coh}_obs.csv' for source in ('real','floor') for coh in ('sepsis','covid')]
outputs=[out/'corrected_summary.csv',out/'corrected_gene_metrics.csv',out/'RESULTS.md',out/'METRIC_DEFINITIONS.md']
manifest={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'bootstrap_reps':200,'seed':20261028,'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'outputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs}}
(out/'run_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('wrote',out/'RESULTS.md')
