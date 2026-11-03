"""Reconcile A3 full-n and matched-n ratios without treating either as a ceiling."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path('F:/PBISC')
SOURCE=ROOT/'bulkonly_lora5_20261006/analysis/integration_report_20261012/week1_defensive_20261027/results/a3_power_matched_ceiling.csv'
OUT=ROOT/'project_review_20261028'
a=pd.read_csv(SOURCE)
a['full_n_to_split_pct']=100*a.floor_r_full/a.ceil_case_split_r
a['matched_n_to_split_pct']=100*a.floor_r_matchedN/a.ceil_case_split_r
assert np.allclose(a.matched_n_to_split_pct,a.pct_of_ceiling_r)
cols=['cohort','celltype','n_case','n_ctrl','floor_r_full','floor_r_matchedN','ceil_case_split_r','full_n_to_split_pct','matched_n_to_split_pct']
a[cols].to_csv(OUT/'a3_ratio_reconciliation.csv',index=False)

lines=[
    '# A3 비교비율 정리 — 2026-10-28',
    '',
    '10/27 A3 원본 CSV에서 두 floor 분자와 같은 real–real 분모를 분리해 표시한다. 원본 값은 변경하지 않았다.',
    '',
    '| 코호트 | 세포형 | real case/ctrl | real–real 분할 r | floor full-n r | floor matched-n r | full-n/분할 | matched-n/분할 |',
    '|---|---|---:|---:|---:|---:|---:|---:|',
]
for _,r in a.iterrows():
    lines.append(f'| {r.cohort} | {r.celltype} | {int(r.n_case)}/{int(r.n_ctrl)} | {r.ceil_case_split_r:.3f} | {r.floor_r_full:.3f} | {r.floor_r_matchedN:.3f} | {r.full_n_to_split_pct:.1f}% | {r.matched_n_to_split_pct:.1f}% |')
lines += [
    '',
    '원본 `pct_of_ceiling_r`은 **matched-n floor / real–real case 분할**이다. 예컨대 Sepsis T-cell의 64.0%와 슬라이드의 71.9%는 동일 성능의 불확실성 범위가 아니라 서로 다른 floor 표본 수에서 나온 별도 비율이다.',
    '',
    'real–real 분할은 case만 둘로 나누고 control은 공유한다. 따라서 분할 양쪽이 독립적이지 않으며 수학적 성능 상한이나 신뢰구간으로 해석할 수 없다. Sepsis HC는 7명뿐이다. COVID도 일부 환자가 양성·음성 양쪽에 있어 환자 의존성을 다시 고려해야 한다.',
    '',
    '원고에서는 이 표를 탐색적 참고 대조로만 사용하고, 비교 목적에 맞는 한 가지 표본 수 방식을 사전 지정한 뒤 patient-aware 재계산을 수행한다. 현재 수정 DEG 결과와 숫자를 교체하거나 합산하지 않는다.',
    '',
    f'입력: [`{SOURCE.name}`]({SOURCE.as_posix()}). 산출: [`a3_ratio_reconciliation.csv`](a3_ratio_reconciliation.csv).',
]
(OUT/'a3_ratio_reconciliation.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('wrote',OUT/'a3_ratio_reconciliation.csv',OUT/'a3_ratio_reconciliation.md')
