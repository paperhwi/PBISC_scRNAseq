# 통일 재분석: PBISC 추가 후 실제로 검출되는 연관성

이 보고서의 숫자·후보 표는 `build_unified_summary.py`가 공통 계산표에서 생성한다. 기존15–24번의 그림·이메일 문장을 복사한 결과표가 아니다. 원본은 수정하지 않았다.

## 답부터

**실제로 남는 후보가 있다.** 가장 분명한 엄격한 Real-missed 희귀 조성 결과는 **HSPC–중증도**다. 전사체에서는 **HSPC KIAA0391–중증도**가 공통 가용5,762개 유전자 screen과 고정밀 조건부 순열 검산에서도 후보로 남는다. 다만 이 두 결과의 통계적 의미와 개인별 cell-state 복원 정확성은 다르다.

IgG/IgA ASC–중증도, HSPC–사망은 Real에서도 이미 명목 신호가 있었다. 정확한 표현은 **PBISC 신규 환자 추가 후 FDR 기준을 넘은 연관성**이다. 모든 새로운 결과를 “Real에 전혀 없던 생물학”으로 부르면 안 된다.

## 왜 폴더마다 오류와 다른 숫자가 반복됐는가

**모든 폴더가 같은 계산을 하고 있지는 않았다.** 원시 patient abundance 자체는 기존15번과25번에서 부동소수점 수준으로 일치한다. 문제는 다음 네 가지가 섞인 것이다.

1. **실제 오류:** HSPC가 부족한 기존 Real 환자2명을 신규 환자 arm에 넣었던31/29명 혼동, source/분모/본문의 잘못된 수치 인용, 예전 고정 SE, count와 실제 분석 cell 수 혼용.
2. **합법적인 다른 estimand:** 환자 평균 후 INT 대 sample INT 후 환자 평균, cell 가중 평균 대 방문 동일 가중 평균, pooled bootstrap 대 source-stratified bootstrap, 서로 다른 time window와 유효 환자 수. 이런 숫자가 다르다고 모두 계산 오류는 아니다.
3. **다른 검정 family와 수치 정밀도:** New BH38 대 fixed BH43, 표적20 genes 대5,762 genes 대모든 celltype×gene, 중복 module17 대고유16, CV row order/seed, 부족한 bootstrap 횟수. 같은 q/P라는 이름으로 비교하면 잘못된 결론이 된다.
4. **결과와 문서의 분리:** 계산표를 고쳐도 제목·본문·범례·이메일의 수치를 수작업으로 복사하면 옛 숫자가 남는다. 이전 감사 역시 일부 수치 재현을 확인한 것을 전체 해석 검증 완료처럼 넓게 보고한 부분이 있었다.

이전 제 감사23번 A08의 Cochran Q 명칭 지적은 철회한다. 실제 두 source에서는 사용된 분산 하에 Q=Z²로 같은 P를 주므로 단순 이름 오류가 아니었다. C116 B-cell 수치와 Calibrated NNLS의 부호 등 이전 감사에서 놓친 본문 오류도 별도 독립 확인했다. 정정 기록 (canonical source: `Errata_to_Prior_Audits_KO.md`), 전체 방법 계보 (canonical source: `methods/Method_Drift_Audit_KO.md`)를 참고하면 된다. 24번의 실제 개선도 계보에 기록했으며, 모든24번 분기가 잘못됐다고 판정하지 않았다.

## 이번 계산은 무엇을 고정했는가

- 환자: Real37 + 겹치지 않는 신규29 =66명/13사망/217채혈. Shared는 Real을 한 번 사용한다. Generated cell과 반복 채혈은 독립 n가 아니다.
- 조성: 동일43개 fine population, RBC/platelet 제외 후 leukocyte100% 재정규화, 환자별 방문 평균, source 내부 INT. **모든 arm에서 constant도 P=1로 남겨 BH43**를 적용한다. Primary는 endpoint별 q<.05, 기존 mortality q<.10과 jointBH86은 별도 비교다.
- 희귀도: Real 환자 평균과 중앙값이 모두1% 미만인 23개를 primary rare로 정의했다. .1/.5/2%도 별도 라벨로 공개하며, 희귀도에 맞춰 FDR 분모를 사후 축소하지 않았다. 이 정의는 보편적 세포 분류법이 아니라 이번 분석의 명시적 operational definition이다.
- RNA: raw rare-subset H5AD에 공통으로 있는14 labels 중 두 arm의 환자·사망 수 조건을 충족한 **9개**를 완전 검사했다. 그중 primary rare는 **7개**다. 공통 **5,762개 symbols**, cell별 CP10K/log1p 후 환자 평균, 환자별 최소3cells다. 전체 인간 유전자나 전체43 celltypes의 RNA를 모두 측정한 분석은 아니다. 미충족5개 label의 coverage도 보존했다.
- 검정: 같은 `analysis_core.py`의 Spearman, feature별 MWU auto, BH, CI를 사용한다. Constant/nonestimable을 조용히 버리지 않는다. RNA P가 matrix 전체의 tie 여부 때문에 잘못 같은 방식으로 계산되지 않도록 feature별 exact/asymptotic rule을 고정했다.
- CI: 같은 helper로50,000회 fixed-score 환자 bootstrap. rho는 source, AUC는 source×outcome 층화다. Generator·라벨·후보선정 불확실성을 모두 포함한 구간은 아니다.
- 후보: **엄격한 Real-missed=Real P≥.05, Combined q<.05**, **FDR rescue=Real q≥.05, Combined q<.05**를 분리한다. 전자는 실제 현상의 부재를 입증한다는 뜻이 아니다.

**Real 비유의→Combined 유의라는 문턱 변화 자체는 두 효과의 차이나 formal augmentation gain을 검정한 결과가 아니다.** 추가 환자·source 구성이 달라지는 현재 비교와, 같은 환자·시점에서 PBISC가 제공한 고유 이득을 검증하는 대조 실험은 구분한다.

계약은 기존 결과를 알고 수행하는 소급 재분석 규칙이며 prospective preregistration으로 부르지 않는다. [분석 계약](../provenance/ANALYSIS_CONTRACT.json), 공통 엔진 (canonical source: `analysis_core.py`)을 보존했다.

## 1. 희귀 세포 조성: 무엇이 새로 보이는가

아래 효과는 중증도 Spearmanρ, 사망 AUC다. AUC는 임의로 위험 방향으로 뒤집지 않았다.

| Population | Endpoint | Real 효과 / P | New 효과 / P | Combined 효과 / q43 | Source 순열 q43 | 분류 |
|---|---|---|---|---|---|---|
| HSPC | 중증도 | 0.09544 / 0.5742 | 0.6588 / 0.000102 | 0.3543 / 0.0168 | 0.01816 | 엄격한 Real-missed |
| B_ASC_IgG | 중증도 | 0.3583 / 0.02946 | 0.4221 / 0.02254 | 0.3915 / 0.007077 | 0.008907 | FDR rescue |
| B_ASC_IgA | 중증도 | 0.3789 / 0.02074 | 0.3884 / 0.03735 | 0.3841 / 0.007818 | 0.01156 | FDR rescue |
| HSPC | 사망 | 0.8095 / 0.01237 | 0.8043 / 0.02248 | 0.8084 / 0.01356 | 0.01397 | FDR rescue |

HSPC는 Real의 35/37명에서 이미 관측됐고 Real 평균 비율은 0.2666%다. IgG/IgA도 대부분 Real에서 이미 관측됐다. 따라서 새 cell identity 자체를 발견한 결과가 아니다.

엄격한 HSPC–중증도 결과는 source별 효과가 다르다. Real에서 약한 rho가 신규군에서 강하므로 단순히 “같은 효과가 n 부족 때문에 숨었다”라고 단정하지 않는다. Source/wave/case-mix 차이를 포함한 확장 코호트의 연관성이다. Source-restricted permutation은 source를 보존하지만 WHO 독립성을 검증하는 검정은 아니다.

HSPC–사망의 전체 연관성은 유지되지만, 이전 원자료 검산의 source×peakWHO exact P≈.077은 중증도 독립성의 확증을 주지 않는다. 전체43개에서는 일부 조성 신호가 소실되기도 하므로 검출 수 증가를 정확도나 모든 생물학의 개선으로 부르지 않는다. 전체 null/opposite/gain/loss는 조성 전수 결과 (canonical source: `composition/COMPOSITION_RESULTS_KO.md`)에 있다.

## 2. 세포 내 RNA: 새 후보가 있는가

### HSPC KIAA0391–중증도

| Cohort | N | ρ | 근사 P | HSPC 5762-gene q | 9-celltype global q | 10M 순열 P | 5762 Bonferroni P |
|---|---|---|---|---|---|---|---|
| Real | 33 | 0.336 | 0.05593 | 0.7312 | 0.7336 | 0.05659 | 1 |
| New | 29 | 0.8029 | 1.602e-07 | 0.0009228 | 0.008306 | 1.4e-06 | 0.008067 |
| Combined | 62 | 0.571 | 1.258e-06 | 0.007249 | 0.04408 | 2e-06 | 0.01152 |

Real의 rho는 중간 크기이고 P는.05 근처다. 그러므로 “Real에는 발현·현상이 없었다”가 아니라 “Real-only에서는 통계 기준에 미달했다”가 맞다. 신규와 Combined의 고정10M 순열 결과에서도 작은 tail probability가 남는다. 이 검산은 발견된 한 gene만으로 family를 축소한 보정이 아니다. **HSPC 내 전체5,762개**, 두 endpoint를 함께 보면 **11,524개**의 Bonferroni 배수를 적용해도 New/Combined는 통과하며 Monte Carlo 상한도 확인했다.

그러나9개 celltype 전체 family의 Bonferroni 통과나 permutation-global BH를 검증한 것은 아니다. 전체9×5,762 분석의 근사 global q와 고정밀 순열의 보정 범위를 혼동하지 않는다. 10M 순열의 정확한 범위 (canonical source: `transcriptome_verification/KIAA0391_10M_INTERPRETATION_KO.md`)에 tail count·seed·MC CI가 있다.

같은 알고리즘으로 재현됐다는 사실만으로 작은 표본의 근사 P까지 정확하다고 보증할 수는 없다. 이 gene의 New 근사 P=1.602e-07와10M 순열 P=1.4e-06는 차이가 있고, 후자의 Monte Carlo 구간으로도 그 차이를 확인했다. 그래서 근사 P를 그대로 확정값으로 받아들이지 않고 조건부 순열과 full-family Bonferroni를 별도로 보존했다. Combined의 근사 P는 이번 순열 MC 구간 안에 있어 오류라고 단정하지 않는다.

| 설정 | Combined N | ρ | HSPC gene q | Global q |
|---|---|---|---|---|
| primary | 62 | 0.571 | 0.007249 | 0.04408 |
| equal_visit | 62 | 0.5585 | 0.01375 | 0.03095 |
| minimum10cells | 53 | 0.5798 | 0.03103 | 0.05327 |
| duplicate_sum | 62 | 0.571 | 0.007249 | 0.04408 |

가중 방식과 최소 cells를 바꿔도 HSPC 내부 gene-family 신호는 남지만 global q는 일부 설정에서.05를 넘는다. 설정마다 표본과 estimand가 달라지므로 가장 유리한 설정을 새 primary로 고르지 않았다.

가장 중요한 제한은 공유 환자의 발현 일치도다. KIAA0391의 Real–PBISC bridge rho는 -0.1815로 음수이며, 같은 날 측정된 truth를 확보한 검증도 아니다. **임상 outcome association을 얻은 것과 환자별 세포 내 발현을 정확히 맞힌 것은 별도 기준**이다. 유전자 기능·rare-cell-specific 기전·예후검사로 확대하지 않는다.

같은 New29/79draw의 원 bulk KIAA0391와 비교한 결과는 bulk 문맥 검증 (canonical source: `verification/Bulk_KIAA0391_context.md`)에 있다. 단순 같은 gene의 bulk 값을 옮긴 설명을 지지하지 않는 결과였지만, 두 상관의 차이가 유의하다는 검정이나 세포 특이성의 증명으로 바꾸지 않았다. PBISC는 다유전자 bulk와 reference를 사용한다.

### 다른 primary rare RNA 결과

| Celltype | Gene | Endpoint | Real P | Combined gene q | Global q | 엄격한 Real-missed |
|---|---|---|---|---|---|---|
| HSPC | KIAA0391 | severity | 0.05593 | 0.007249 | 0.04408 | True |
| T g/d | ADK | severity | 0.002343 | 0.01959 | 0.04408 | False |
| T g/d | MXRA7 | death | 0.0001182 | 0.03266 | 0.266 | False |

Tγδ의 ADK·MXRA7는 Real에서도 명목 P가 이미 작았던 FDR rescue다. MXRA7의 whole-celltype global q는.05를 통과하지 않는다. Primary rare가 아닌 CD14mono_IFN의 유전자 결과는 rare-cell 발견으로 세지 않았다. NK CD56bright도 이번1% 기준에는 해당하지 않지만 전체 분석에서 제외하지 않았다.

모든 eligible celltype을 합친 Combined 결과는 celltype별 gene-family 기준 중증도 20개/사망 1개, global family에서는 4개/0개다. Global family에 두 endpoint를 함께 넣은 경우는 0개/0개다. 서로 다른 family의 q를 하나의 유의성처럼 합쳐 쓰지 않는다.

## 3. 기존 HSPC module·20-gene 주장은 어떻게 되는가

공유 universe에서 PCNA/MCM5가 빠져 기존 Cycling과 Cycle4의 실제 gene weights가 완전히 같았다. 동일 P라는 이유만이 아니라 동일 hypothesis이므로 alias로 묶어 **고유16개 module**을 같은 환자에서 전부 재계산했다. 이름17개를 그대로 세는 legacy 값도 함께 보존했다.

Ery5–중증도는 Real P=0.005488/q16=0.08781에서 Combined P=0.002778/q16=0.04445로 FDR rescue가 된다. 엄격한 Real-missed는 아니며, equal-visit 가중에서는 Combined q16=0.1538다. 두 endpoint joint32 q=0.06206도 함께 제시한다.

Ery5–사망 primary q16=0.06206는.05를 넘는다. 조건부 순열이나 최소 cells 민감도에서 더 작은 값이 나와도 그것을 골라 primary 확증이라고 하지 않는다. 기존 표적20개에서는 사망 관련 6개(GATA1, KLF1, TFRC, BLVRB, ITGA2B, TYMS)가 BH20을 통과하지만, 이들은 모든5,762개 gene 보정과는 다르다. [module 전체](../data/modules/HSPC_all_module_associations.csv), [표적20개 전체](../data/transcriptome_verification/HSPC_historical20_all_tests.csv)를 모두 공개했다.

## 무엇을 주장할 수 있는가

1. **지지됨:** 실제 신규 bulk 환자 추가로 일부 희귀 population의 임상 연관성을 더 명확히 검출했다. HSPC–중증도는 엄격한 nominal Real-missed에 해당한다.
2. **후보 수준에서 지지됨:** 추론된 HSPC 내부의 KIAA0391–중증도 연관성은 같은 common-gene universe와 조건부 순열에서도 남는다. 환자별 복원 정확성과 외부 생물학 검증은 별도로 필요하다.
3. **주장하지 않음:** 새로운 cell identity 확정, 모든 희귀 세포 프로그램의 회복, severity-independent death biomarker, 검증된 임상 예후검사, 인과 기전, generator 전체 학습 독립성.

추가로 가장 직접적인 검증은 같은 방문의 bulk/scRNA에서 test donor의 모든 관측 scRNA를 adaptation/reference에서 가린 masking 실험이다. Association recovery와 개인별 expression fidelity를 분리하여 평가해야 한다. 이번 재분석은 새 환자 수집이나 그 실험을 완료한 것이 아니다.

## 검산과 재현

동일 core에 대해 실제 composition/합성322 feature를 scalar SciPy와 비교했고, CI의16개 열×513회도 scalar rank/쌍별 AUC와 대조했다. RNA는9×5,762×3=155,574 feature record의 rho/P/AUC/P와 독립 BH를 전수 비교했다. 최대 차이는 부동소수점 수준이었다. HSPC module204행/408 endpoint 검정도 독립 확인했다.

`PRIMARY_CANDIDATES.csv`, `RESULTS_SUMMARY.json`와 이 문서는 같은 입력 표에서 자동 생성된다. 방법 계보 (canonical source: `methods/Method_Drift_Audit_KO.md`), RNA 전수 검산 (canonical source: `methods/Actual_RNA_Independent_Verification.json`), 조성 검산 (canonical source: `composition/validation.json`), module 검산 (canonical source: `transcriptome/INDEPENDENT_HSPC_MODULE_VERIFICATION.json`)을 함께 보존했다. 입력·코드 hash와 실행 명령은 README 및 run manifest를 따른다.

환자 단위 추론의 문헌 배경은 [Zimmerman et al.2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7854630/)과 [Squair et al.2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8479118/)에서 확인했다. 이 문헌이 본 PBISC 결과를 외부 검증하는 것은 아니다. 분석 절차 지침 출처도 기록했다: Kassis T,Agarwal V,He Y,Patel D,Brueckner AM(2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065).
