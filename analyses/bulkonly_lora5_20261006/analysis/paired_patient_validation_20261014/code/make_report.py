"""Generate an auditable Korean report directly from final saved tables."""
from common import *
import json

s=pd.read_csv(OUT/'results/paired_accuracy_summary.csv')
ps=pd.read_csv(OUT/'results/covid_patient_effect_summary.csv')
m=pd.read_csv(OUT/'results/covid_mapping_patient_summary.csv')
cf=pd.read_csv(OUT/'results/covid_crossfit_effect_summary.csv')
contrasts=pd.read_csv(OUT/'results/paired_method_contrasts.csv')
methods=json.load(open(OUT/'results/covid_mapping_methods.json'))
def interval(r,col='L1',lo='L1_lo',hi='L1_hi'):return f'{r[col]:.3f} [{r[lo]:.3f}, {r[hi]:.3f}]'
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(str(c) for c in row)+' |' for row in rows])
primary=s[(s.scope=='common_test_paired')&(s.level=='l1')&(s.draw==-1)]
ordered=['other_patient_null','PBISC_k1000_eligible_ensemble','PBISC_k4000_eligible_ensemble']
display={'other_patient_null':'다른 평가 환자의 평균 비율','PBISC_k1000_eligible_ensemble':'PBISC K1 ensemble','PBISC_k4000_eligible_ensemble':'PBISC K4 ensemble'}
rows=[];display_rows=[]
for coh,title in [('ye2026_sepsis','Sepsis'),('imperial_eskd','ESKD COVID')]:
    for method in ordered:
        r=primary[(primary.cohort==coh)&(primary.method==method)].iloc[0]
        rows.append([title,display[method],int(r.n_donors),int(r.n_samples),interval(r)])
        display_rows.append(dict(cohort=title,method=display[method],n_patients=int(r.n_donors),n_samples=int(r.n_samples),
                                 L1_lower_is_better=f'{r.L1:.3f}',CI95_lo=f'{r.L1_lo:.3f}',CI95_hi=f'{r.L1_hi:.3f}'))
pd.DataFrame(display_rows).to_csv(OUT/'results/key_paired_results.csv',index=False,encoding='utf-8-sig')
text=['# PBISC paired 정확도 및 환자 단위 검증',
      '검증일: 2026년 10월 14일. 기존 분석 파일을 보존하고 저장된 paired 예측과 실제 단일세포 자료로 재평가했다.',
      'Sepsis에서는 학습 환자를 모두 제외한 다섯 모델 ensemble이 다른 환자의 평균 비율보다 정확했다. COVID에서는 K4가 K1보다 소폭 개선됐지만 평균 비율 기준값 대비 우위는 확정되지 않았다. 환자 단위 재표집 후에도 기존 매핑의 COVID 질병 효과는 실제 효과보다 작게 추정된다. 아래 환자 분리 검증은 라벨 전달 성능을 평가하며, 생성 세포의 실제 정답을 입증하는 검사는 아니다.',
      '## Paired 정확도',
      '주 평가의 정답은 기존 paired benchmark의 frozen scANVI 관측 L1 비율이다. 각 adapter의 학습 donor와 그 donor의 모든 방문을 제외했다. 엄격한 common-test 집합에서는 평가 환자가 다섯 adapter 중 어느 것의 학습에도 사용되지 않았다. paired 생성량은 샘플당 1,500세포이며 K1과 K4 모두 같은 관측 정답과 여섯 유형을 사용한다.',
      table(['코호트','방법','환자 수','샘플 수','평균 L1 및 95% 환자 bootstrap CI'],rows),
      'L1은 비율 차이의 절댓값 합으로 범위는 0–2이고 낮을수록 좋다. 방문별 오차를 환자 내 평균한 후 환자를 동일 비중으로 평가했다. CI는 환자 2,000회 bootstrap이며 주어진 checkpoint와 생성 결과에 조건부다. 다른 환자 평균 기준값은 자기 환자를 제외한 나머지 평가 환자의 관측 비율 평균이다. 다른 평가 환자의 scRNA 정보를 이용하는 reference-rich 기준값이므로 실제 추론 환경과 동일한 정보 예산의 모델은 아니다.',
      '다섯 모델 ensemble은 서로 다른 다섯 학습 donor 집합을 활용한다. 단일 LoRA-5 모델과 같은 정보량으로 취급하지 않는다. 단일 draw별로 동일한 5명에서 만든 평균 비율 기준값도 별도 평가했고 결과는 paired_accuracy_summary.csv에 저장했다. 단일 draw의 기준값 대비 우위는 일관되지 않는다.',
      '### K 차이와 개인 추적',]
krows=[]
for coh,title in [('ye2026_sepsis','Sepsis'),('imperial_eskd','ESKD COVID')]:
    r=primary[(primary.cohort==coh)&(primary.method=='paired_K4_minus_K1')].iloc[0]
    q=primary[(primary.cohort==coh)&(primary.method=='K4_within_group_pairing_permutation')].iloc[0]
    krows.append([title,interval(r),int(q.n_donors),f'{q.p_pairing:.4f}'])
text += [table(['코호트','K4−K1 L1 및 95% CI','개인 추적 permutation 환자 수','permutation p'],krows),
         '개인 추적 검사는 관측 및 생성 비율에서 임상 그룹 평균을 제거하고 같은 그룹 내 환자 예측을 1,000회 섞는 검증이다. 방문을 먼저 환자·상태별로 평균했고 하나의 임상 상태만 존재하는 환자만 교환했다. p 값은 해당 조건에서 자기 환자와의 연결이 무작위 연결보다 유리한지를 나타낸다. 정확한 절대 비율이나 미세한 세포 상태 재현을 보증하지 않는다. K 차이는 seed와 checkpoint에 조건부이며 sampler seed가 맞춰진 인과적 K 실험이 아니다.',
         '전체 paired 집합의 결과도 제공한다: Sepsis 82환자·82샘플, COVID 57환자·157샘플. common-test는 각각 60환자·60샘플과 35환자·97샘플이다. Sepsis 평가에는 소아도 포함된다. 성인 59샘플, COVID 회복기를 제외한 집합 및 출판된 로컬 표본 ID 집합은 별도 sensitivity scope로 제공한다.',
         '## 같은 성인 Sepsis 샘플에서 deconvolution 비교',
         '기존 deconvolution의 성인 60샘플 중 paired PBISC·recap truth와 모두 일치하는 58샘플을 평가했다. recap의 major 유형을 여섯 유형으로 합쳐 모든 방법을 같은 정답·샘플·분모로 평가했다. PBISC는 저장된 frozen scANVI 라벨을 이 어휘에 대응시킨다. 라벨 정의의 차이는 남는다. deconvolution의 기존 2-fold reference·5명 calibration과 PBISC의 adapter 학습 donor 집합은 동일하지 않으므로 이 비교는 평가 조건을 맞춘 보조 분석이며 학습 조건까지 통제한 우위 실험은 아니다.']
dec=s[(s.scope=='common60_recap_truth')&(s.draw==-1)]
text += [table(['방법','샘플 수','L1 및 95% 환자 bootstrap CI'],[[r.method,int(r.n_samples),interval(r)] for _,r in dec.iterrows()]),
         'ensemble과 단일 draw는 구분해야 한다. 단일 K4 adapter의 이 공통 집합 L1은 0.360–0.515로 draw에 따라 달랐다. paired_method_contrasts.csv에는 같은 환자의 오차 차이와 조건부 CI를 저장했다. Wilcoxon p 값은 탐색적이며 다중비교 보정된 확증 검정으로 사용하지 않는다.',
         '## 환자가 겹치지 않는 COVID 라벨 검증',
         '187샘플·61환자의 실제 counts에서 샘플당 최대 400세포를 무작위 선택했다. 같은 고정 세포 집합으로 기존 sample split과 patient split을 재평가했다. 기존 split은 reference/test에 환자 31명이 겹친다. 새 2-fold split은 각각 학습·평가 환자의 교집합이 0이며 61환자 모두 정확히 한 번 out-of-fold 평가했다. HVG 2,000개, scaling, PCA 30개 성분과 분류기 학습은 reference에서만 수행했다.',
         'R1은 Harmony+kNN, R2는 reference PCA+kNN, R3는 PBISC read-out을 흉내 낸 reference로 학습한 logistic regression이다. C=0.05를 고정했고 충분한 반복으로 수렴을 확인했다. R1은 각 protocol의 fold 0에서 PBISC형 실제 세포 대조군을 평가했다. R2/R3는 patient 2-fold 전체를 평가했으며 원래 실제 세포와 PBISC형 실제 세포를 모두 검사했다. 아래 정확도는 정답을 알고 있는 PBISC형 실제 세포의 환자별 평균이다. sample-overlap은 49 test 환자, patient R1은 한 fold의 환자, patient R2/R3는 61 OOF 환자이므로 교환 가능한 동일 테스트 집합의 인과 효과로 해석하지 않는다.']
mr=m[(m.source=='realT2')&m.level.isin(['major','fine'])]
text += [table(['분할','경로','라벨','평가 환자','정확도 및 95% CI','비율 L1'],
               [[r.protocol,r.route,r.level,int(r.n_patients),interval(r,'patient_accuracy','accuracy_lo','accuracy_hi'),f'{r.composition_L1:.3f}'] for _,r in mr.iterrows()]),
         'labels의 정답은 기존 recap annotation으로, 저자 원래 annotation이나 독립 실험의 생물학적 정답과 동일하다고 가정하지 않는다. read-out emulation은 형식에 대한 대조군이며 생성 모델의 정확한 복제는 아니다. precision/recall과 recoverability는 환자 비중을 맞춰 covid_mapping_recoverability.csv에 저장했다. 분류기 reject 기준을 새로 튜닝하지 않았으며 본 검증은 강제 라벨 할당의 결과다.',
         '## 환자 단위 COVID 질병 효과',
         '샘플별 조성에서 환자·임상 상태별 평균을 만든 후 CLR(case)−CLR(control)을 계산했다. CLR pseudocount는 0.001이다. 같은 환자의 음성·양성 방문 및 실제·생성 source는 bootstrap에서 동일한 재표집 multiplicity를 받아 공분산을 보존한다. 실제 환자는 양성 37명, 음성 37명이고 일부가 양쪽에 포함된다. 생성 환자는 양성 58명, 음성 19명이다. 주 비교는 기존 Harmony 전달 라벨을 유지해 통계 단위 변경의 영향을 평가했다.']
p=ps[(ps.contrast=='infection')&(ps.wave=='all')&ps.level.isin(['major','fine'])&ps.source.isin(['k1000','k4000'])]
text += [table(['라벨','K','효과 기울기 및 95% 환자 CI','방향 일치율'],[[r.level,r.source,interval(r,'slope','slope_lo','slope_hi'),f'{100*r.sign_agree:.1f}%'] for _,r in p.iterrows()]),
         '기울기 1은 실제와 같은 효과 크기다. 방향 일치율은 실제 평균 비중이 한쪽 그룹에서 0.5% 이상인 유형 전체에 대한 값이며 기존 보고서의 유의한 neighbourhood만 고른 97–100%와 다른 지표다. adapter 학습 환자·sex-flag 표본 제외 sensitivity도 저장했다.',
         '새 reference classifier를 적용한 효과는 아래 표에 별도로 제시한다. 생성 방문은 자신의 환자가 reference 학습에 포함되지 않은 fold에서만 매핑했다. 실제 reference에 없는 환자는 두 reference fold의 조성 추정치를 동일 비중으로 평균했다. 참조 환자 누출을 제거해도 generated state의 정답을 알 수 있는 것은 아니다.']
c=cf[cf.level.isin(['major','fine'])&cf.source.isin(['k1000','k4000'])]
text += [table(['경로','라벨','K','환자 효과 기울기 및 95% CI','방향 일치율'],[[r.route,r.level,r.source,interval(r,'slope','slope_lo','slope_hi'),f'{100*r.sign_agree:.1f}%'] for _,r in c.iterrows()]),
         'wave별 sensitivity와 중증도 비교는 covid_patient_effect_summary.csv에 함께 저장했다. 실제 wave 2에는 양성/회복기만 있고 비교할 음성 표본이 없어 wave 2 감염 효과를 독립적으로 추정하지 않았다. 중증도는 active positive 방문에서 환자의 peak WHO 상태로 정의하고 회복기를 제외했다. 전체 infection 효과에 post-infection severity를 공변량으로 넣지는 않았다. site·wave·질병 효과의 분리는 별도 연구 설계가 필요하다.',
         '## 이번 결과로 가능한 주장과 남는 한계',
         '- Sepsis의 대계통 비율에는 평균적인 질병 그룹 조성 이상으로 개별 샘플을 추적하는 정보가 있다. 특히 다섯 모델 ensemble에서 평균 비율 기준값보다 낮은 paired 오차가 확인된다.',
         '- COVID는 개인 연결 정보와 K4의 소폭 개선이 있지만, 평균 비율 기준값보다 절대 비율 정확도가 확실히 좋다는 결론은 유지하기 어렵다.',
         '- 라벨 전달의 환자 누출은 실제 RNA를 다시 학습·평가하는 방식으로 제거했다. 새 환자 성능은 위 표의 경로·라벨·대조군 조건에서만 해석한다.',
         '- 고정 마커 상태 검증, pediatric sepsis의 동일 임상 비교, DEX 노출/마커, 독립 sampler seed 및 새로운 LoRA donor set의 training uncertainty는 이번 검증으로 해결되지 않는다.',
         '## 결과 파일과 재실행',
         'results에는 paired 정확도·개별 오차·유형별 상관·방법 간 차이, 환자 분할 명부·매핑 confusion counts·recoverability, 환자 단위 및 새 cross-fit 질병 효과를 저장했다. figs에는 PNG와 editable SVG를 저장했다. code는 이 폴더를 기준으로 경로를 해석하고 기존 입력을 읽기 전용으로 사용한다. staged_manifest.json은 압축 번들에서 복원한 입력의 SHA256, VALIDATION.json은 독립 검증 결과와 실행 라이브러리 버전을 기록한다.',
         '재실행 순서: paired_accuracy.py → paired_contrasts.py → patient_statistics.py → patient_mapping.py → mapping_summary.py → crossfit_effects.py → figures.py → make_report.py → validate.py. 번들 Python 실행파일을 사용한다. patient_mapping.py는 고정 counts cache와 완료된 fold를 재사용한다.',
         '![Paired 정확도]('+str(OUT/'figs/paired_accuracy.png').replace('\\','/')+')',
         '![환자 단위 질병 효과]('+str(OUT/'figs/covid_patient_effects.png').replace('\\','/')+')',
         '![환자 분리 매핑 검증]('+str(OUT/'figs/covid_patient_mapping.png').replace('\\','/')+')']
(OUT/'REPORT.md').write_text('\n\n'.join(text)+'\n',encoding='utf-8')
print('Report saved:',OUT/'REPORT.md')
