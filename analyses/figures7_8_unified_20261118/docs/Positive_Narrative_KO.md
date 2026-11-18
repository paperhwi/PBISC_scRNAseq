# PBISC의 임상 연관성 탐색 범위를 강조하는 서술안

작성 기준: 2026-11-18, canonical 25번 통일 재분석. 이 문서는 원고 문장과 Figure 7·8의 패널 구성안이다. 기존 연관성 계산과 그림 파일은 바꾸지 않았다. 주요 후보의 효과크기·P·q는 기존 CSV에서 직접 불러왔다.

## 1. 권장 중심 주장과 제목

**PBISC를 통해 실제 bulk-only 환자를 세포 유형별 분석에 포함하자, Real-only에서 검정 기준에 미달했던 희귀 HSPC의 중증도 연관성과 HSPC 내 KIAA0391의 전사체 연관성이 검출되었다. 동시에 HSPC의 사망 연관성과 ASC의 중증도 연관성은 다중검정 기준에서 추가로 뒷받침되었다.**

이 논문의 기여는 scRNA-seq가 있는 환자에 제한됐던 임상 연관성 탐색을 실제 추가 환자에게 확장하는 데 둔다. 논리의 순서는 **실제 환자 편입 → 희귀 세포 조성의 임상 연관성 → 세포 유형별 유전자 연관성 → 직접 검증할 후보의 도출**이다.

권장 통합 제목:

> **PBISC-enabled cohort expansion reveals clinical associations of rare HSPCs and candidate transcriptomic correlates**
>
> **PBISC 기반 코호트 확장으로 드러나는 희귀 HSPC와 전사체 후보의 임상 연관성**

Figure 7 제목:

> **PBISC-enabled cohort expansion reveals severity-associated rare HSPCs and strengthens evidence for outcome-associated cell populations.**

Figure 8 제목:

> **Cell-type-resolved profiling identifies KIAA0391 within HSPCs as a candidate correlate of COVID-19 severity.**

## 2. ‘새로 보였다’를 세 가지 성과로 명확히 나누기

1. **Real-only에서 명목 기준에도 미달했던 연관성의 검출:** HSPC 비율–중증도, HSPC KIAA0391–중증도. 이 분석에서 strict Real-missed는 Real P≥.05이면서 Combined의 해당 primary family q<.05라는 조작적 정의다. 생물학적 현상의 부재를 뜻하지 않는다. 특히 KIAA0391의 Real P=0.05593는 경계값이다.
2. **기존 명목 연관의 다중검정 근거 강화:** HSPC 비율–사망, IgG/IgA ASC 비율–중증도, Tγδ ADK–중증도, Tγδ MXRA7–사망. 이들은 Real에서도 P<.05였고 추가 후 해당 family의 q<.05가 됐다. FDR rescue라고 명시한다. Strict Real-missed도 넓은 FDR rescue 정의에 포함되므로 두 범주를 중복 합산하지 않는다.
3. **세포 유형별 가설 생성의 확장:** bulk-only 환자에서 HSPC에 할당된 발현값과 임상 중증도의 연관 후보를 도출했다. 세포 유형별 분석이 가능해졌다는 성과와, 해당 개인의 실제 HSPC 발현을 정확히 복원했는지는 별도 평가 대상이다.

‘희귀’는 Real 환자 비율의 평균·중앙값이 모두1% 미만이라는 이번 분석의 정의다. HSPC는 Real 평균 0.2666%, 중앙값 0.1343%로 이 조건을 충족한다. HSPC 자체는 Real 35/37명에서 관측되었으므로, 최초 세포 정체성 발견이라는 문장은 쓰지 않는다. NK CD56bright와 CD14mono_IFN은 이1% 기준의 희귀 세포에 포함하지 않는다.

## 3. 성과의 정확한 수치와 보정 범위

중증도 효과는 Spearman ρ, 사망 효과는 AUC다. AUC는 사망자에서 값이 높을수록1에 가까우며 위험 방향에 맞춰 임의 반전하지 않았다. 아래 q는 각 행의 family에 해당하고 서로 같은 보정 범위가 아니다. 모든 q는 endpoint별 primary 결과다.

| 결과 | Real 효과 / P | 추가군 효과 | Combined 효과 / q | q의 family | 성과 분류 |
|---|---|---|---|---|---|
| HSPC 비율–중증도 | 0.095 / 0.5742 | 0.659 | 0.354 / 0.0168 | 43 populations per endpoint | Real P≥.05 → Combined q<.05 |
| IgG ASC 비율–중증도 | 0.358 / 0.02946 | 0.422 | 0.391 / 0.007077 | 43 populations per endpoint | Real 명목 연관의 FDR rescue |
| IgA ASC 비율–중증도 | 0.379 / 0.02074 | 0.388 | 0.384 / 0.007818 | 43 populations per endpoint | Real 명목 연관의 FDR rescue |
| HSPC 비율–사망 | 0.810 / 0.01237 | 0.804 | 0.808 / 0.01356 | 43 populations per endpoint | Real 명목 연관의 FDR rescue |
| HSPC KIAA0391–중증도 | 0.336 / 0.05593 | 0.803 | 0.571 / 0.007249 | 5762 shared genes within HSPC | Real P≥.05 → Combined q<.05 |
| Tγδ ADK–중증도 | -0.505 / 0.002343 | -0.672 | -0.567 / 0.01959 | 5762 shared genes within T g/d | Real 명목 연관의 FDR rescue |
| Tγδ MXRA7–사망 | 0.021 / 0.0001182 | 0.119 | 0.041 / 0.03266 | 5762 shared genes within T g/d | Real 명목 연관의 FDR rescue |
| HSPC Ery5–중증도 | 0.473 / 0.005488 | 0.280 | 0.374 / 0.04445 | 16 unique historical HSPC modules | Real 명목 연관의 FDR rescue |

RNA의 endpoint별 global BH는9개 eligible cell types×5,762genes=51,858검정이다. KIAA0391과 ADK의 Combined global q는 각각 0.04408, 0.04408다. 두 endpoint까지 합친 global-joint q는 각각 0.08816, 0.08816로 .05를 넘는다. MXRA7는 해당 Tγδ gene family q=0.03266를 통과하지만 global q=0.266는 통과하지 않는다.

## 4. 본문에 사용할 한국어 Results 초안

### PBISC를 통한 실제 환자 코호트 확장이 희귀 세포의 임상 연관성을 드러낸다

PBISC를 이용하여 기존 단일세포 측정 환자37명에 bulk-only 신규 환자29명을 추가하고, 총66명에서 세포 조성과 임상 결과의 연관성을 평가하였다. 반복 채혈과 생성 세포 수가 아닌 환자를 독립 분석 단위로 사용했으며, 겹치는 환자는 한 번만 포함하였다. 모든 arm에 동일한43개 세포군과 endpoint별 다중검정 보정을 적용하여, 신규 환자의 편입 후 어떤 연관성이 검정 기준을 충족하는지 평가하였다.

이 접근은 희귀 HSPC의 중증도 연관성을 드러냈다. HSPC 비율과 환자별 최고 WHO 중증도의 연관성은 Real-only에서 검출되지 않았으나(ρ=0.095, P=0.5742), 추가군에서 뚜렷하게 관측되었으며(ρ=0.659), 통합 코호트에서 다중검정 기준을 충족하였다(ρ=0.354, q43=0.0168). Source 내부에서 임상 라벨을 순열한 검정에서도 이 연관성이 유지되었다(q43=0.01816). 따라서 기존 Real-only 분석에 추가 환자 정보를 연결함으로써, 희귀 HSPC의 중증도 관련 변화를 임상 연관 후보로 포착할 수 있었다.

코호트 확장은 HSPC의 사망 연관성에 대해서도 일관된 근거를 제공하였다. HSPC 비율의 사망 관련 AUC는 Real, 추가군, 통합군에서 각각 0.810, 0.804, 0.808로 유사하였고, 통합 분석에서 q43=0.01356를 보였다. 이는 Real에서도 관측된 명목 연관이 추가 환자의 정보와 함께 FDR 기준을 충족한 결과이다. IgG ASC와 IgA ASC 역시 두 source에서 같은 방향의 중증도 연관성을 보였으며, 통합 분석에서 각각 ρ=0.391, q43=0.007077 및 ρ=0.384, q43=0.007818로 다중검정 기준을 충족하였다.

### 세포 유형별 발현 분석은 HSPC 내 중증도 관련 전사체 후보를 제시한다

조성 분석에 이어, 공통14개 세포 라벨 중 환자 수·사망 수 조건을 충족한9개 라벨에서 공통5,762개 유전자를 조사하였다. 이 중7개 라벨은 위 희귀도 정의를 충족하였다. HSPC의 발현 분석에는 Real33명과 추가29명, 총62명이 포함되었다. HSPC 내 KIAA0391은 Real-only에서 명목 유의수준에 미달했지만(ρ=0.336, P=0.05593), 추가군에서는 강한 양의 중증도 연관성을 보였으며(ρ=0.803), 통합 분석에서 해당 HSPC 전체 유전자 family의 보정 기준을 충족하였다(ρ=0.571, q5762=0.007249). Endpoint별9개 세포 유형 전체의 근사 검정 기반 global BH에서도 q=0.04408였다.

KIAA0391의 통합 연관성은 source를 보존한 고정10M 규모 순열 검산에서도 지지되었다(P=2.0e-06). 이 순열 P에 HSPC 전체5,762개 유전자와 두 임상 endpoint를 모두 포함한11,524개 검정의 Bonferroni 배수를 적용해도 보정 P=0.023048로 기준을 충족하였다. 이 결과는 후보 한 개만을 분모로 삼지 않고 전체 HSPC family를 유지한 보수적 민감도 분석에서도 남는 전사체 연관 후보임을 보여준다.

추가적으로, Tγδ 세포의 ADK는 중증도와 음의 연관성을 보였고(통합 ρ=-0.567, 해당 gene-family q=0.01959; endpoint-global q=0.04408), 같은 세포군의 MXRA7는 사망자에서 더 낮은 발현 방향의 연관성을 보였다(AUC=0.041, 해당 gene-family q=0.03266). 이들은 Real에서 관측된 명목 연관의 FDR rescue이며, 특히 MXRA7는 추가군 사망자가3명이고 전역 보정 기준을 충족하지 않아 보조 후보로 제시하였다.

이 결과들은 PBISC를 통한 코호트 확장이 희귀 세포 조성과 세포 유형별 전사체의 임상 연관성 탐색을 넓힐 수 있음을 뒷받침한다. HSPC 중증도 효과의 source별 크기 차이와 KIAA0391 shared-patient 상관(ρ=-0.181, n=24)을 고려하면, 현재의 성과는 후향적 연관 후보의 도출이며 독립 표본의 직접 세포 측정으로 검증할 가치가 있다.

## 5. Discussion에서 성과를 가장 강하게 연결하는 방식

**첫째, PBISC의 실질적 가치는 임상 연관성 분석에 참여할 수 있는 환자 범위를 넓히는 데 있다.** 현재 결과는 희귀 HSPC의 조성 연관성과 HSPC 내 유전자 연관성을 한 분석 흐름에서 제시한다. 조성과 전사체는 서로 다른 관측량이며, 두 결과를 동일한 분화 경로나 기전의 인과적 연쇄로 단정하지 않는다. 그럼에도 어느 희귀 세포군과 어느 분자 후보에 후속 측정을 집중할지 제시한다는 점에서 연구 설계상 가치가 있다.

**둘째, 새 후보의 식별과 기존 신호의 증거 강화가 함께 관측된다.** HSPC–중증도와 KIAA0391–중증도는 Real-only의 명목 기준에 미달했던 후보이고, HSPC–사망은 두 source에서 약0.81의 일관된 효과와 통합 FDR 통과를 보인다. ‘모든 결과가 새로 발견됐다’는 하나의 주장보다, 서로 다른 두 성과를 병렬로 설명하는 편이 데이터의 장점을 더 잘 드러낸다.

**셋째, 세포 유형별 문맥을 이용할 후속 연구의 근거를 제공한다.** 동일 신규29명의 단일 유전자 bulk KIAA0391 점수에서는 중증도 연관성을 검출하지 못했다(ρ=0.167, P=0.3867). 이 비교는 세포 유형별 분석을 추가로 탐구할 이유가 된다. 다만 bulk와 추론 HSPC 상관의 차이를 직접 검정하지 않았고, 독립된 실제 HSPC 발현 측정도 없으므로 세포 특이성이나 bulk 대비 우월성의 확증으로 쓰지 않는다.

한계 문단은 다음처럼 묶어 쓸 수 있다:

> 본 분석은 실제 신규 환자가 추가된 후향적 코호트 비교이며, 단일세포 측정 유무가 같은 환자·시점에서 무작위로 결정된 대조 실험은 아니다. HSPC의 중증도 연관 효과는 Real과 추가 source에서 크기가 달랐으며(서술적 이질성 P=0.007643,43개 보정 q=0.06515), source 내 순열 검정은 모든 임상 교란을 제거하지 않는다. 또한 KIAA0391의 shared-patient 상관(ρ=-0.181, P=0.3961, n=24)은 높은 개인별 복원 정확성을 뒷받침하지 않으며, 방문 시점의 차이 등을 고려한 직접 검증이 필요하다. Real 비유의와 Combined 유의의 차이 자체는 효과 차이 또는 PBISC 고유의 성능 이득에 대한 검정이 아니다. 사망 연관성 역시 중증도 독립 예후 예측이나 임상적 유용성의 확증으로 해석하지 않는다.

## 6. Figure 7 권장 구성: 희귀 세포의 임상 연관성

본문 메시지는 HSPC를 중심으로 하고, ASC는 성과를 넓히는 보조 축으로 둔다. 다섯 패널 구성안을 권장한다.

| 패널 | 보여줄 내용 | 가장 중요한 읽기 포인트 |
|---|---|---|
| 7A | Real37 + 실제 신규29 =66명; 반복 채혈·shared 중복 제외;43개 동일 family | PBISC가 실제 환자의 세포 유형별 분석 참여를 가능하게 함 |
| 7B | HSPC 비율과 최고 WHO 중증도: 환자 점, source별 구분 | 희귀도와 임상 분포를 직접 제시; 점의 수는 환자 수 |
| 7C | HSPC–중증도 Real/New/Combined ρ와50k bootstrap95% CI | strict Real-missed, source별 효과 차이, Combined q43와 source순열 q43 |
| 7D | HSPC–사망 세 arm AUC와95% CI | 약.81의 효과 일관성과 통합 FDR 통과; 예후 모델 성능으로 명명하지 않음 |
| 7E | IgG/IgA ASC–중증도 세 arm 효과와 CI | 같은 방향의 명목 신호가 FDR 기준을 충족하는 두 사례 |

7B에서 raw 비율을 보여주면 통계는 source 내 INT 기반이라는 점을 축·범례에서 구분한다. Raw pooled 산점도의 기울기를 보정 분석의 회귀선처럼 그리지 않는다. 모든 패널에서 Real/New/Combined 순서와 색을 고정한다. CI는 fixed-score 환자 bootstrap이며 생성 모델·세포 라벨 불확실성을 포함하지 않는다고 legend에 쓴다.

전체43개에서는 severity q<.05 결과가10→12개로 변했고,6개 증가·4개 소실이 함께 있었다. Death는1→2개였다. 이 전수 맥락을 Results/legend 한 문장에 남기고, 전체 effect·P·q·gain/loss 표와 rare23개 결과는 보충자료에 보존한다. 일부 예만으로 전체 탐지 정확도가 향상됐다고 쓰지 않는다.

권장 Figure 7 legend 핵심 문장:

> Associations were tested at the patient level using a fixed43-population family within each endpoint. Strict Real-missed denotes Real P≥.05 and Combined q<.05; FDR rescue denotes Real q≥.05 and Combined q<.05. These labels describe threshold transitions rather than a formal test of augmentation benefit. Mortality AUCs summarize univariable associations and do not establish severity-independent prediction.

## 7. Figure 8 권장 구성: HSPC 전사체 후보와 검증 범위

기존 erythroid/20-marker 중심의 광범위한 기전 주장보다, 전수 공유 유전자 screen에서 남은 KIAA0391을 중심으로 구성한다. KIAA0391을 erythroid marker나 기전 인자로 명명하지 않는다.

| 패널 | 보여줄 내용 | 가장 중요한 읽기 포인트 |
|---|---|---|
| 8A | 공통14 labels→eligible9(그중rare7),5,762genes; HSPC33+29=62명 | 분석 가능한 범위와 실제 검정 분모 |
| 8B | KIAA0391–중증도 Real/New/Combined ρ와50k CI | Real P=.0559의 기준 미달 후보가 통합 gene-family 보정 통과 |
| 8C | 추가29명의 환자별 HSPC KIAA0391와 최고 WHO 중증도 | 추가 코호트에서 관측된 강한 연관; 독립 검증군이라는 라벨은 사용하지 않음 |
| 8D | 근사 BH와 고정10M 순열 기반 full-HSPC Bonferroni를 분리한 검증 요약 | q5762, endpoint-global q, joint-global q, Bonf11524 각각의 범위를 명시 |
| 8E | Shared24명의 관측/추론 HSPC KIAA0391 비교 | 집단의 임상 연관성 도출과 개인별 발현 fidelity의 차이를 투명하게 제시 |

8D는 유리한 보정 하나만 제시하지 않는다. Main의 요약에 local q=0.007249, analytic endpoint-global q=0.04408, analytic global-joint q=0.08816,10M within-HSPC two-endpoint Bonferroni P=0.023048를 별도 열로 둔다.10M은 B=9,999,999회 재표본화와 plus-one 식으로 계산되었다. 모든 유전자에 대한 permutation-global BH를 계산한 것으로 표기하지 않는다.

동일 방문 가중, 최소10cells, 중복 symbol 합산의 결과는 고정 후보의 민감도 표로 연결한다. 최소10cells 설정에서 local q는 통과하지만 endpoint-global q는 0.05327으로 .05를 넘는 점도 남긴다. ‘모든 설정에서 전역 보정을 통과’라는 문장은 쓰지 않는다.

## 8. 연관 보충 Figure와 Table 배치

| 보충자료 | 내용 | 역할 |
|---|---|---|
| Supplementary Figure A | 전체43 조성, gain/loss, rarity cutoffs, source별 효과 | 주요 양성 결과를 전수 맥락에 위치시키기 |
| Supplementary Figure B | RNA14-label coverage,9-label full screen, 전체 null·opposite 결과 | 가용성에 따른 분석 범위와 선택 과정을 공개 |
| Supplementary Figure C | Tγδ ADK–severity, MXRA7–death; source별 환자·사망 수와 효과 | 희귀 세포의 보조 전사체 성과; MXRA7 local-family 및3 New deaths 한계 표시 |
| Supplementary Figure D | 같은 신규29명의 bulk KIAA0391, 추가된 세포 문맥, shared24 timing | bulk 단일 점수와 다른 분석 문맥; 세포 특이성 검증과는 구별 |
| Supplementary Figure E | Ery5와 과거20-marker 패널, family별 q와 민감도 | 가설 기반 탐색적 지원 근거 |
| Supplementary Tables | 조성43,9×5,762 genes, 모든module, 보정family·CI·순열·coverage | 선택적 보고를 막는 완전한 결과 보존 |

Ery5–중증도는 고유16-module family q=0.04445이나 두 endpoint joint32 q=0.06206이고 방문 동일 가중에서는 q=0.1538다. Ery5–사망의 primary q16=0.06206은 .05를 넘는다. 과거20-marker의 사망6개 BH20 양성은 전체5,762개 유전자 발견으로 바꾸어 쓰지 않는다. Ery5 구성 유전자를 모두 TF로 부르거나 ambient-free lineage identity로 단정하지 않는다.

이 보충 Figure 이름은 구성상의 가칭이다. 기존 원고의 S-number를 조사하지 않고 실제 번호를 배정하지 않았다.

## 9. 다음 분석: 중심 주장을 더 강하게 만드는 우선순위

1. **환자 단위 masking 실험:** bulk와 scRNA가 모두 있는 환자를 대상으로, reference군과 실제 세포 측정을 숨길 검증군을 분리한다. 검증군 환자가 학습·라벨 튜닝에 섞이지 않도록 generator와 전처리를 fold 안에서 고정한다. Reference-only, reference+PBISC, reference+실측 세포의 oracle을 비교하고 bulk/deconvolution 비교군도 둔다. 여러 mask에서 효과 추정·후보 회수·null 위양성률을 비교하면 실제 환자 수 증가와 PBISC 고유 기여를 분리할 수 있다. Null 검증은 후보 선택을 포함한 전체 평가 절차를 반복한다.
2. **같은 시점의 실제 HSPC 검증:** 신규 환자의 같은 채혈에서 HSPC 비율을 직접 계수하고, 분리 HSPC의 KIAA0391을 측정한다. 현재의 관측/추론 shared-patient 비교는 같은 시점의 독립 진실값을 대신하지 않는다. 이 검증이 성공하면 집단 연관성에서 실제 세포 상태 해석으로 전진할 수 있다.
3. **사망과 중증도·채혈 시점 구별:** 사전에 정의한 초기 채혈/landmark와 소수의 임상 공변량을 사용한다. 현재 사망13명으로 고차원 조정 모델을 과도하게 적합하지 않는다. 최고 WHO 중증도를 baseline 예후 정보로 취급하지 않는다.
4. **모델 및 라벨 불확실성 전파:** 현재 fixed-score CI와 별개로 generator seed/reference resampling, 반복 채혈의 환자 cluster, fine ASC 라벨 변화를 평가한다. ASC 세부 isotype은 직접 라벨 검증 전까지 해당 annotation의 연관성으로 제시한다.
5. **독립 코호트의 후보 고정 검증:** HSPC 조성–severity/death, HSPC KIAA0391–severity를 소수의 고정 가설로 검증한다. Tγδ 사망 후보는 더 많은 실제 사망 표본에서 재평가한다. 현재 발견에 사용한 신규29명을 독립 재현군으로 다시 세지 않는다.

## 10. 발표·초록용 짧은 문장

> **PBISC를 통한 코호트 확장은 희귀 HSPC의 중증도 연관성과 HSPC 내 KIAA0391 전사체 후보를 드러냈고, HSPC 사망 및 ASC 중증도 연관성의 다중검정 근거를 강화했다. 이는 bulk-only 환자를 세포 유형별 임상 가설 생성에 포함하는 전략을 뒷받침한다.**

‘통계적으로 유의한 임상 결과 연관성’이라고 표현할 수 있다. ‘임상적으로 유용한 biomarker’, ‘중증도 독립 사망 예측’, ‘새 세포 정체성의 발견’, ‘정확한 개인별 전사체 복원’은 현재 결과보다 강하다.

## 근거 파일

- [Canonical 통합 결과](Canonical_Results_KO.md)
- [8개 후보의 계산표](../data/PRIMARY_CANDIDATES.csv)
- [RNA 고정 후보의 보정 범위·bridge·희귀도](../data/transcriptome/fixed_primary_candidates_with_rarity_and_composition.csv)
- [KIAA0391 고정10M 순열 해석](KIAA0391_10M_Interpretation_KO.md)
- [동일 신규29명의 bulk 문맥](Bulk_KIAA0391_Context.md)
- [원고용 영문 Results](Results_Draft_EN.md)
- [서술용 원자료와 hash](../provenance/NARRATIVE_QA.json)

이 서술안은 새로운 외부 생물학적 기능 주장이나 문헌상 기전 결론을 추가하지 않는다.
