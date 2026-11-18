# KIAA0391 중증도 연관성: 고정1천만 규모 순열 검산

2026-11-18. 대상은 이미 선택된 HSPC KIAA0391 한 유전자다. 이 작업의 목적은 작은 permutation tail의 정밀도를 높이는 것이며, 새 유의 결과를 찾거나 유의해질 때까지 반복하는 분석이 아니다. 실행 전 `KIAA0391_10M_LOCKED_PLAN.json`에 각 arm의 **B=9,999,999**, seed와 방법을 고정했고, 정해진 횟수에서 종료했다.

## 관측 효과와 순열 결과

원 NPZ의 patient cell-mean logCP10K에서 canonical cohort 및 최소3세포 기준으로 Real33명, 실제 신규29명, 통합62명을 추출했다. Source INT 후 관측 global rank를 고정하고, globally ranked WHO severity를 source 안에서만 섞었다. 양측 통계량은 절대 centered-rank cross-product이며 동률 이상을 모두 포함했다.

| arm | n | rho | 원 Spearman 근사P | 극단 순열 K/B | raw K/B | +1 순열P | tail의95% MC구간 |
|---|---:|---:|---:|---:|---:|---:|---|
|Real|33|.335969|.0559324|565,941/9,999,999|.056594106|.0565942|[.05645097,.05673750]|
|New|29|.802863|1.601592×10⁻⁷|13/9,999,999|1.300000×10⁻⁶|**1.4×10⁻⁶**|[6.92195×10⁻⁷,2.22304×10⁻⁶]|
|Combined|62|.571011|1.258105×10⁻⁶|19/9,999,999|1.900000×10⁻⁶|**2.0×10⁻⁶**|[1.14392×10⁻⁶,2.96708×10⁻⁶]|

구간은 K~Binomial(B,p)에 대한 beta quantile의 양측 Clopper–Pearson95% 구간이다. 이는 고정된 expression/환자/WHO의 **순열 tail probability를 추정하는 MC 불확실성**이며 임상 효과의 CI가 아니다. +1P는 `(K+1)/(B+1)`이다. 원 rho/근사P는 기존1M 표와 공통 core의 값을 재현했다.

**New에서는 근사 P와 순열 P의 차이가 MC 해상도만의 문제는 아니다.** +1 순열P는 근사P의8.74배이며 MC95% 하한도 근사P의4.32배다. 작은 표본과 강한 WHO 동률이 있는 이 자료에서 Spearman의 t 근사 tail과 조건부 permutation tail을 같은 숫자로 취급하면 안 된다. 다만 기존1M의4×10⁻⁶은 극단 순열3개에 기반한 값이고, 그 넓은 구간 안에 이번1.4×10⁻⁶이 들어가므로 이전 순열 구현의 오류를 뜻하지는 않는다.

Combined의 순열P는 근사P의1.59배지만, 근사P가 이번 MC95% 구간 안에 있어 두 tail의 차이를 이 정밀도로 확정할 수 없다. Real은 두 방법 모두 .05보다 크다. “Real P≥.05”라는 분류는 유지되지만 .056 수준의 경계 P가 생물학적 효과 부재를 뜻하지 않는다.

## 전체 family를 유지한 Bonferroni 민감도

아래는 선택한 유전자1개로 family를 줄이지 않고, 명시된 전체 family 배수를 +1 순열P와 MC구간에 적용한 **보수적 FWER 민감도**다. 모든 값은1에서 cap했다. 완전한 gene-family permutation BH q를 재계산한 표가 아니다.

| 전체 family m | New 보정P[확대한 MC95%구간] | Combined 보정P[확대한 MC95%구간] |
|---|---|---|
|5,762: HSPC×한 endpoint|.008067[.003988,.012809]|.011524[.006591,.017096]|
|11,524: HSPC×두 endpoint|.016134[.007977,.025618]|.023048[.013183,.034193]|
|51,858:9 celltypes×한 endpoint|.072601[.035896,.115282]|.103716[.059322,.153867]|
|103,716:9 celltypes×두 endpoint|.145202[.071792,.230565]|.207432[.118643,.307734]|

따라서 **HSPC 내 전체5,762유전자 및 두 endpoint11,524검정**까지는 New/Combined 모두 보수적 Bonferroni .05 기준을 통과하며 MC상한도 .05 아래다.9개 세포군을 모두 포함한 family에서는 Bonferroni 기준을 통과하지 않는다. New의51,858검정 MC구간은 .05를 가로지르므로 해당 FWER 경계의 정밀도는 여전히 제한적이다.

**Bonferroni 미통과는 원 global BH q의 유의성 소실을 증명하지 않는다.** BH는 전체 유전자 P의 순서와 분포를 사용한다. 이 작업은 KIAA0391만의 tail을 계산했으므로 원 global analytic BH q를 그대로 permutation-validated q라고 부를 수 없고, 반대로 이 보수적 민감도만으로 global BH가 실패했다고 단정할 수도 없다. 선택한 유전자끼리만 BH 보정해서 이를 대신해서도 안 된다.

## 수치·재현 검증과 해석 범위

2배 centered average ranks를 정수로 표현했다. 모든 가능한 dot-product의 절댓값 상한은 Real15,776/New9,240/Combined97,920으로2²⁴보다 작으므로 float32 정수 곱/합이 정확히 표현된다. 첫50,000개 모든 순열 및 **관측값 이상인 모든 순열**을 int64로 다시 계산해 최대 오차0을 확인했다. Real 동률768개를 포함했으며 New/Combined의 절대 관측값 동률은0개였다. Random-key 동률 재추출은0회다. 총 소요 약105초다.

이 검산은 원 source-normalized score에 조건부인 임상 association의 강건성을 평가한다. PBISC 생성/annotation 불확실성, source의 절대 expression calibration, cohort/time/treatment confounding 또는 formal augmentation gain을 검증하는 실험은 아니다. New가 pooled 결과의 일부이므로 독립 외부 검증으로 서술하지 않는다.

재현 파일:

- `verify_KIAA0391_source_permutation_10M.py`:실제 실행 코드.
- `KIAA0391_10M_LOCKED_PLAN.json`, `KIAA0391_10M_RUN_MANIFEST.json`:사전 고정 횟수·seed·plan/input/code hash·실행 시각.
- `KIAA0391_source_permutation_10M.csv`:세 arm의 효과, analyticP, K, raw/+1P, MC구간, 정수 대조.
- `KIAA0391_source_permutation_10M_Bonferroni_sensitivity.csv`:네 전체 family의 보정P와 MC상하한.
- `KIAA0391_source_permutation_10M_patient_ranks.csv`:원 환자 membership/원 logCP10K/sourceINT/관측 순위.
- `KIAA0391_10M_validation.json`, `KIAA0391_10M_OUTPUT_MANIFEST.json`:검산 및 최종 산출물 SHA.
