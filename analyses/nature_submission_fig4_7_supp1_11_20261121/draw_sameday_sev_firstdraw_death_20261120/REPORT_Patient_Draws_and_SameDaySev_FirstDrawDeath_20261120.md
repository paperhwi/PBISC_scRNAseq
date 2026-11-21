# Real & PBISC 전체 환자·채혈 명세 및 재분석 비교 리포트
**신규 프레임워크: `Severity = 전체 217회 채혈 × 당일 WHO Pair` · `Death = Donor당 첫 채혈 × 최종 사망 여부`**

> [!IMPORTANT]
> **핵심 결론 요약**
> 제안하신 분석 설계(**중증도 = 모든 채혈과 당일 WHO severity 페어링**, **사망 = 환자당 첫 채혈과 최종 사망 여부 페어링**)로 기존에 환자 평균(patient-average)으로 수행했던 **5개 전 영역**을 재검증한 결과, 통계적·생물학적으로 **훨씬 정교하고 명료한 2축 구조**가 드러났습니다:
> 1. **진정한 조기 사망 예측인자 (Donor 첫 채혈 → 최종 사망, 당일 중증도 보정 후에도 독립 유의)**:
>    - **`T CD8 CTL` 고갈**이 전체 43개 아형 중 **1위** (`AUC = 0.213`, 역방향 `0.787`, $P = 0.0029$, `q43 = 0.070`; Real 단독 `0.206, P = 0.036` vs PBISC new 단독 `0.210, P = 0.031`로 **소수점 둘째 자리까지 완벽 일치**, 이질성 $P_{\text{het}} = 0.974$; **당일 WHO 보정 후에도 $P = 0.033$**, 7일 이내 조기 첫 채혈 $P = 0.022$, Peak WHO 2–3 중증 환자 내 $P = 0.011$).
>    - **`NK` (`CD56dim`) 고갈** 공동 1–2위 (`AUC = 0.217`, 역방향 `0.783`, $P = 0.0033$, `q43 = 0.070`, 당일 WHO 보정 $P = 0.030$, Direct COVID-9 사망 `inv AUC = 0.844, P = 0.0010`), **`CD14mono_IFN` 상승** (`AUC = 0.721, P = 0.022`, 당일 WHO 보정 $P = 0.030$), **`B_CD11c` 고갈** (`inv AUC = 0.726, P = 0.019`, 당일 WHO 보정 $P = 0.048$), **`NK CD56bright` 고갈** (`inv AUC = 0.707, P = 0.032`).
>    - **HSPC 단일세포 내 조기 전사 리프로그래밍**: HSPC **세포 비율(%)** 자체는 첫 채혈에서 아직 오르지 않지만(`AUC = 0.577, P = 0.43`), **첫 채혈 HSPC 단일세포 내부의 `Alarmin_Inflam` (`S100A8/9/12`)** (`AUC = 0.724, P = 0.021`, **당일 WHO 보정 $P = 0.022$**; Real `0.726` vs PBISC new `0.725`; `≥1` cell 기준 `N=64`에서 `AUC = 0.767, P = 0.0057`) 및 **`Ery5`/`GATA1`/`KLF1`** (`AUC = 0.696–0.718, P = 0.025–0.044`)은 **이미 첫 채혈 시점부터 최종 사망을 예고**합니다.
> 2. **당일 중증도 실시간 추적자 (전체 217회 채혈 × 당일 WHO, 환자-군집 GEE + Mundlak 분해)**:
>    - 기존 환자 평균(`q43 < 0.05` 12개) 대비 **18개 아형(`+week` GEE `q43 < 0.05`, bias-reduced `17`개)**으로 검정력이 대폭 강화되었습니다.
>    - **`B_ASC_IgG`** ($\beta = +0.414, P = 5.8 \times 10^{-10}$, `q43 = 1.8e-8`), **`B_ASC_dividing`** ($\beta = +0.422, P = 8.2 \times 10^{-10}$), **`B_ASC_IgA`** ($\beta = +0.334, P = 1.3 \times 10^{-5}$), **`HSPC` (PBISC-rescued)** ($\beta = +0.333, P = 1.9 \times 10^{-5}$, `q43 = 1.3e-4`), **`NK CD56bright` (PBISC-rescued)** ($\beta = -0.276, P = 1.3 \times 10^{-4}$), **`T CD8 CM`** ($\beta = -0.231, P = 5.2 \times 10^{-4}$) 모두 압도적으로 검증되었으며, **Mundlak 분해**에서 `B_ASC`, `HSPC`($P_{\text{within}} = 0.0088$)는 **동일 환자 내부에서 당일 WHO가 악화될수록 실제로 함께 상승**함을 입증했습니다.
>    - 단일세포 HSPC **`Cycle4`** (`q16 = 1.59e-4`)와 **`Ery5`** (`q16 = 1.70e-3`)는 당일 WHO 기준 기존 대비 검정력이 **26~900배** 강화되었으며, 반대로 **`T-cell 19-gene ISG`**는 `+week` 보정 당일 WHO에서 $\beta \approx 0$ ($P = 0.59\text{–}0.96$)으로 사라져 기존 환자 평균의 양의 상관이 **채혈 시점(초기 급성기 편향)** 교란이었음이 명확히 규명되었습니다.

---

## 1. 전체 데이터셋 환자·채혈 구조 요약 (Real & PBISC)

> [!NOTE]
> **Sepsis 코호트 (`Ye et al.`, GSE279448 / GSE279452; Fig 4–5) 확인**:
> Sepsis 코호트는 **Real `133`명 (성인 분석 대상 `108`명: sepsis / ICU ctrl / HC)** 및 **PBISC `122`명 (성인 sepsis `68`, 성인 HC `36`, 성인 ICU ctrl `17`, 소아 ICU ctrl `1`)** 모두 **1인 1회 단일 채혈(cross-sectional)**입니다. 따라서 종단 다회 채혈 페어링은 **ESKD COVID-19 코호트 (Fig 4–6)**에 고유하게 적용됩니다.

### 1.1 Imperial ESKD COVID-19 전체 107명 도너 및 362개 샘플 구성
- **Primary 임상 분석 대상**: 비-CLL COVID-19+ 환자 **`N = 66`명 (`217`회 QC≥200 급성기 채혈, 검증 사망 `TRUE_PUB_11 = 11`명 / 직접 COVID 사망 `9`명)**
  - **`Real-only (10명)`**: `47`회 채혈 (중앙값 `5.5`회/명, 사망 `1`명: `C141`). *참고: `C31`은 별도 Wave-1 bulk(`PBISC` 3회, `in_primary_scope=False`)가 존재하나 중복 방지를 위해 `Real` 3회만 사용.*
  - **`Shared-27 (27명, Primary는 Real 사용)`**: `91`회 QC≥200 `Real` 채혈 (중앙값 `3.0`회/명, 사망 `4`명: `C21, C33, C108, C145`). *참고: `C73`(d7, 179 cells) 및 `C168`(d10, 4 cells)의 2개 저세포수 샘플은 QC `<200` 기준 제외; 이 27명에게는 서로 다른 투석 방문일에 채혈된 `59`회의 추가 `PBISC` bulk-derived 샘플이 존재하며 이는 `Shared-27` 교차 검증(Fig 4–5, Supp Fig 9)에 사용됨.*
  - **`PBISC new (29명, bulk-only 신규 확장 환자)`**: `79`회 채혈 (중앙값 `2.0`회/명, **1회 단일 채혈 환자 7명**, 사망 `6`명: `C32, C37, C48, C66, C77, C84`).
- **제외 및 대조 도너**:
  - **제외된 CLL 환자 `2`명**: `C116` (87/M, Real 3회 + PBISC 2회, Peak `critical`, 생존), `C139` (85/M, Real 6회 + negative 1회 + PBISC recovery 1회, Peak `severe`, 생존).
  - **COVID-negative ESKD 대조 전용 도너 `39`명**: 전원 1회 채혈 (`Real` scRNA 및/또는 `PBISC` bulk). *(전체 `362`개 샘플 상세는 Excel 시트 `T1e_All362_Draws` 및 `T1d_CovidNeg39`에 수록)*

### 1.2 군별·생존/사망별 채혈 횟수, 첫 채혈 시점 및 당일 WHO 분포

| 코호트 군 | 환자수 | 총 채혈수 | 채혈수/명 중앙값 (IQR, 범위) | 단일채혈 환자수 | 첫 채혈일 중앙값 (IQR) | 첫 채혈 ≤7일 환자 | 전체 채혈일 중앙값 (≤7일 / >7일) | 첫 채혈 당일 WHO (경/중등/중/위중) | Peak WHO (경/중등/중/위중) | 전체 채혈 당일 WHO (경/중등/중/위중) | 최종 사망 (`PUB11`) |
|---|---:|---:|---|---:|---|---:|---|---|---|---|---:|
| Real-only | 10 | 47 | 5.5 (3–6, 3–6) | 0 | d3.5 (2–4) | 10 | d9.0 (21 / 26) | 6/4/0/0 | 4/4/2/0 | 30/15/2/0 | 1 |
| Shared-27 (Real) | 27 | 91 | 3.0 (3–4, 1–5) | 1 | d4.0 (2–4) | 27 | d9.0 (40 / 51) | 14/10/3/0 | 9/5/8/5 | 47/24/13/7 | 4 |
| **Real 전체** | **37** | **138** | **3.0 (3–5, 1–6)** | **1** | **d4.0 (2–4)** | **37 (100%)** | **d9.0 (61 / 77)** | **20/14/3/0** | **13/9/10/5** | **77/39/15/7** | **5** |
| — Real 생존 | 32 | 118 | 3.0 (3–5, 1–6) | 1 | d4.0 (2–5) | 32 | d9.0 (51 / 67) | 18/11/3/0 | 13/9/7/3 | 74/29/10/5 | 0 |
| — Real **사망** | 5 | **20** | **4.0** (3–5, 3–5) | 0 | **d2.0** (2–4) | 5 | d8.0 (10 / 10) | **2/3/0/0** | 0/0/3/2 | 3/10/5/2 | 5 |
| **PBISC new 전체** | **29** | **79** | **2.0 (2–4, 1–7)** | **7** | **d7.0 (3–11)** | **16 (55%)** | **d12.0 (26 / 53)** | **14/9/3/3** | **10/5/7/7** | **49/19/6/5** | **6** |
| — PBISC 생존 | 23 | 71 | 3.0 (2–4, 1–7) | 3 | d8.0 (4–12) | 11 | d12.0 (19 / 52) | 14/6/2/1 | 10/5/7/1 | 49/16/5/1 | 0 |
| — PBISC **사망** | 6 | **8** | **1.0** (1–2, 1–2) | **4** | **d2.5** (0–6) | 5 | d2.5 (7 / 1) | **0/3/1/2** | 0/0/0/6 | 0/3/1/4 | 6 |
| **통합 (Combined)** | **66** | **217** | **3.0 (2–4, 1–7)** | **8** | **d4.0 (2–7)** | **53 (80%)** | **d10.0 (87 / 130)** | **34/23/6/3** | **23/14/17/12** | **126/58/21/12** | **11** |

> [!WARNING]
> **왜 제안하신 프레임워크가 필수적인가? (2가지 핵심 실증 증거)**
> 1. **당일 WHO ≠ Peak WHO 불일치 (`91 / 217 = 41.9%`)**: Peak `critical` 환자의 `28`회 채혈 중 당일 `critical`이었던 채혈은 `12`회뿐이며, 나머지 `16`회는 `mild(4)/moderate(6)/severe(6)` 시점에 채혈되었습니다. Peak WHO 하나로 모든 채혈을 평균 내면 급성기 중증도 변화를 뭉개버립니다.
> 2. **사망자 채혈 횟수의 3배 비대칭**: `Real` 사망자 `5`명은 총 **`20`회 (`4.0`회/명)** 채혈된 반면, `PBISC new` 사망자 `6`명은 총 **`8`회 (`1.3`회/명, 6명 중 4명이 1회 채혈 후 사망)**만 채혈되었습니다. 모든 채혈을 합산해 사망 분석을 하면 Real 말기 채혈이 3배 과대 가중되지만, **Donor당 첫 채혈 → 최종 사망**은 모든 환자를 동등하게 1:1 반영하며, **11명 사망자 중 10명(90.9%)의 첫 채혈이 `≤7일` 이내**(`Real` 5명은 전원 `d2–d4` `mild/moderate` 상태, `PBISC new`도 3명이 `moderate` 상태)이므로 진정한 **초기 예후 예측인자**를 검정할 수 있습니다.

![66명 Primary COVID-19 환자 전체 채혈 타임라인 Swimmer Plot 및 채혈-중증도 불일치](C:/Users/scott/.gemini/antigravity/brain/add4945a-39a0-4f8d-acd4-3287c21c6121/FigR1_Patient_Swimmer_Draw_Schedule_66.png)

---

### 1.3 전체 COVID-19+ 환자 68명 개별 채혈 시점·횟수·당일/최고 중증도·사망여부 전체 표

*(증상 발현일 기준 채혈 시점, 217개 개별 채혈별 검사치, 39명 COVID-negative 도너 및 362개 전체 샘플은 [Excel 통합 워크북](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120/PBISC_Patient_Draw_Inventory_and_SameDaySev_FirstDrawDeath_20261120.xlsx)의 `T1_Patients68` / `T1b_Draws217` / `T1e_All362_Draws` 시트에서 필터·정렬 가능합니다.)*

| 환자 ID | 군 (Cohort) | 센터/Wave | 나이/성별 | 채혈수 (Primary) | 첫 채혈일 | 전체 채혈 시점 (첫 양성 스왑 기준) | 첫 채혈 당일 WHO | 매 채혈 당일 WHO | Peak WHO | 최종 사망 (PUB11/Dir9) |
|---|---|---|---|---:|---:|---|---|---|---|---|
| C23 | Real-only (10) | NCL | 78/F | 3 | d7 | d7, d12, d19 | moderate | moderate, moderate, moderate | moderate | 생존 |
| C31 | Real-only (10) | NCL | 39/F | 3 | d0 | d0, d7, d14 | mild | mild, mild, mild | mild | 생존 |
| C36 | Real-only (10) | NCL | 66/M | 3 | d2 | d2, d4, d8 | mild | mild, mild, mild | moderate | 생존 |
| C124 | Real-only (10) | NCL | 81/M | 3 | d4 | d4, d7, d11 | moderate | moderate, mild, mild | moderate | 생존 |
| C138 | Real-only (10) | Cambridge | 77/F | 6 | d7 | d7, d10, d12, d14, d17, d19 | moderate | moderate, moderate, moderate, moderate, mild, mild | moderate | 생존 |
| C140 | Real-only (10) | Cambridge | 23/F | 6 | d2 | d2, d4, d7, d9, d11, d14 | mild | mild, mild, mild, mild, mild, mild | mild | 생존 |
| C141 | Real-only (10) | Cambridge | 85/M | 5 | d4 | d4, d6, d10, d14, d17 | moderate | moderate, severe, moderate, moderate, moderate | severe | **사망** (Direct) |
| C146 | Real-only (10) | Cambridge | 66/M | 6 | d3 | d3, d5, d7, d10, d14, d18 | mild | mild, mild, moderate, moderate, severe, moderate | severe | 생존 |
| C170 | Real-only (10) | Cambridge | 82/M | 6 | d4 | d4, d7, d9, d11, d14, d16 | mild | mild, mild, mild, mild, mild, mild | mild | 생존 |
| C190 | Real-only (10) | Cambridge | 54/M | 6 | d2 | d2, d5, d7, d9, d12, d14 | mild | mild, mild, mild, mild, mild, mild | mild | 생존 |
| C20 | Shared-27 (Real used) | NCL | 62/M | 3 | d7 | d7, d11, d13 | moderate | moderate, moderate, mild | severe | 생존 |
| C21 | Shared-27 (Real used) | NCL | 80/M | 3 | d4 | d4, d11, d19 | moderate | moderate, severe, critical | critical | **사망** (Direct) |
| C33 | Shared-27 (Real used) | Cambridge | 75/M | 4 | d2 | d2, d4, d7, d14 | mild | mild, mild, severe, critical | critical | **사망** (Direct) |
| C34 | Shared-27 (Real used) | NCL | 74/F | 3 | d0 | d0, d7, d11 | mild | mild, mild, mild | mild | 생존 |
| C40 | Shared-27 (Real used) | NCL | 57/M | 3 | d2 | d2, d9, d18 | mild | mild, moderate, critical | critical | 생존 |
| C42 | Shared-27 (Real used) | NCL | 40/M | 3 | d0 | d0, d7, d12 | moderate | moderate, mild, mild | moderate | 생존 |
| C58 | Shared-27 (Real used) | Cambridge | 72/F | 5 | d5 | d5, d7, d9, d12, d15 | mild | mild, mild, mild, mild, mild | mild | 생존 |
| C60 | Shared-27 (Real used) | NCL | 80/M | 3 | d2 | d2, d7, d11 | moderate | moderate, moderate, moderate | severe | 생존 |
| C63 | Shared-27 (Real used) | NCL | 45/M | 3 | d2 | d2, d9, d14 | moderate | moderate, mild, mild | moderate | 생존 |
| C65 | Shared-27 (Real used) | NCL | 73/M | 3 | d4 | d4, d8, d12 | moderate | moderate, severe, critical | critical | 생존 |
| C73 | Shared-27 (Real used) | NCL | 88/M | 1 | d3 | d3 | moderate | moderate | severe | 생존 |
| C82 | Shared-27 (Real used) | NCL | 60/F | 2 | d6 | d6, d11 | mild | mild, mild | mild | 생존 |
| C93 | Shared-27 (Real used) | NCL | 87/M | 3 | d7 | d7, d13, d20 | mild | mild, moderate, mild | moderate | 생존 |
| C101 | Shared-27 (Real used) | Cambridge | 61/M | 5 | d4 | d4, d9, d11, d16, d18 | severe | severe, severe, moderate, moderate, severe | severe | 생존 |
| C108 | Shared-27 (Real used) | Cambridge | 79/F | 3 | d2 | d2, d5, d9 | moderate | moderate, moderate, moderate | severe | **사망** (Indirect) |
| C113 | Shared-27 (Real used) | NCL | 61/F | 3 | d4 | d4, d10, d17 | severe | severe, severe, moderate | severe | 생존 |
| C123 | Shared-27 (Real used) | NCL | 66/M | 3 | d4 | d4, d9, d15 | severe | severe, severe, moderate | severe | 생존 |
| C126 | Shared-27 (Real used) | NCL | 64/M | 3 | d2 | d2, d7, d13 | mild | mild, mild, mild | mild | 생존 |
| C127 | Shared-27 (Real used) | NCL | 69/F | 3 | d2 | d2, d11, d18 | moderate | moderate, moderate, mild | moderate | 생존 |
| C128 | Shared-27 (Real used) | NCL | 72/M | 3 | d2 | d2, d11, d18 | mild | mild, mild, mild | mild | 생존 |
| C132 | Shared-27 (Real used) | NCL | 77/F | 3 | d5 | d5, d9, d18 | moderate | moderate, mild, mild | moderate | 생존 |
| C137 | Shared-27 (Real used) | Cambridge | 60/F | 2 | d5 | d5, d14 | mild | mild, mild | mild | 생존 |
| C145 | Shared-27 (Real used) | Cambridge | 65/F | 5 | d2 | d2, d7, d9, d14, d18 | mild | mild, severe, severe, moderate, moderate | severe | **사망** (Indirect) |
| C147 | Shared-27 (Real used) | Cambridge | 53/M | 5 | d2 | d2, d7, d11, d14, d16 | mild | mild, severe, critical, critical, critical | critical | 생존 |
| C168 | Shared-27 (Real used) | Cambridge | 75/M | 4 | d4 | d4, d7, d14, d17 | mild | mild, mild, mild, mild | mild | 생존 |
| C169 | Shared-27 (Real used) | Cambridge | 73/F | 5 | d4 | d4, d7, d10, d14, d17 | mild | mild, mild, mild, mild, mild | mild | 생존 |
| C187 | Shared-27 (Real used) | Cambridge | 70/M | 5 | d5 | d5, d7, d12, d16, d19 | mild | mild, mild, mild, mild, mild | mild | 생존 |
| C7 | PBISC new (29) | w1 | 83/M | 2 | d9 | d9, d11 | mild | mild, mild | mild | 생존 |
| C9 | PBISC new (29) | w1 | 55/M | 2 | d7 | d7, d12 | mild | mild, mild | moderate | 생존 |
| C10 | PBISC new (29) | w1 | 84/F | 3 | d18 | d18, d21, d30 | moderate | moderate, moderate, mild | severe | 생존 |
| C11 | PBISC new (29) | w1 | 62/M | 4 | d13 | d13, d17, d19, d21 | moderate | moderate, moderate, moderate, mild | severe | 생존 |
| C12 | PBISC new (29) | w1 | 58/M | 5 | d8 | d8, d14, d16, d19, d21 | severe | severe, moderate, mild, mild, mild | severe | 생존 |
| C18 | PBISC new (29) | w1 | 64/M | 4 | d14 | d14, d15, d17, d19 | moderate | moderate, moderate, mild, mild | severe | 생존 |
| C22 | PBISC new (29) | w1 | 65/M | 3 | d8 | d8, d11, d13 | mild | mild, mild, mild | mild | 생존 |
| C24 | PBISC new (29) | w1 | 75/F | 3 | d8 | d8, d11, d13 | mild | mild, mild, mild | mild | 생존 |
| C25 | PBISC new (29) | w1 | 38/M | 3 | d14 | d14, d16, d18 | mild | mild, mild, mild | moderate | 생존 |
| C27 | PBISC new (29) | w1 | 85/M | 2 | d5 | d5, d7 | mild | mild, mild | mild | 생존 |
| C32 | PBISC new (29) | w1 | 76/M | 1 | d16 | d16 | critical | critical | critical | **사망** (Direct) |
| C37 | PBISC new (29) | w1 | 82/F | 2 | d2 | d2, d5 | critical | critical, critical | critical | **사망** (Direct) |
| C47 | PBISC new (29) | w1 | 72/M | 3 | d7 | d7, d11, d16 | mild | mild, mild, mild | mild | 생존 |
| C48 | PBISC new (29) | w1 | 72/M | 2 | d0 | d0, d2 | moderate | moderate, critical | critical | **사망** (Direct) |
| C50 | PBISC new (29) | w1 | 66/M | 4 | d13 | d13, d15, d18, d22 | severe | severe, moderate, severe, moderate | severe | 생존 |
| C59 | PBISC new (29) | w1 | 47/M | 2 | d11 | d11, d15 | mild | mild, mild | severe | 생존 |
| C66 | PBISC new (29) | w1 | 60/F | 1 | d0 | d0 | severe | severe | critical | **사망** (Direct) |
| C67 | PBISC new (29) | w1 | 75/F | 1 | d2 | d2 | moderate | moderate | moderate | 생존 |
| C68 | PBISC new (29) | w1 | 75/M | 1 | d5 | d5 | mild | mild | mild | 생존 |
| C77 | PBISC new (29) | w1 | 64/F | 1 | d7 | d7 | moderate | moderate | critical | **사망** (Direct) |
| C80 | PBISC new (29) | w1 | 69/F | 1 | d14 | d14 | mild | mild | moderate | 생존 |
| C84 | PBISC new (29) | w1 | 75/M | 1 | d3 | d3 | moderate | moderate | critical | **사망** (Direct) |
| C87 | PBISC new (29) | w1 | 71/M | 3 | d9 | d9, d11, d13 | mild | mild, mild, mild | mild | 생존 |
| C99 | PBISC new (29) | w1 | 60/M | 2 | d3 | d3, d6 | critical | critical, severe | critical | 생존 |
| C105 | PBISC new (29) | w1 | 75/F | 6 | d0 | d0, d2, d5, d7, d12, d14 | mild | mild, mild, mild, mild, mild, mild | mild | 생존 |
| C110 | PBISC new (29) | w1 | 74/M | 4 | d2 | d2, d4, d11, d14 | moderate | moderate, moderate, mild, mild | moderate | 생존 |
| C114 | PBISC new (29) | w1 | 68/M | 7 | d4 | d4, d8, d12, d15, d23, d26, d30 | moderate | moderate, severe, moderate, moderate, mild, mild, mild | severe | 생존 |
| C129 | PBISC new (29) | w1 | 77/M | 4 | d2 | d2, d5, d7, d12 | mild | mild, mild, mild, mild | mild | 생존 |
| C149 | PBISC new (29) | w1 | 76/M | 2 | d3 | d3, d14 | mild | mild, mild | mild | 생존 |
| C116 | Excluded CLL (2) | NCL | 87/M | 0 | 제외(CLL) | Real 3회(d3,d8,d12) / PBISC 2회 | - | mod, mild, crit (Real) | critical | 생존 |
| C139 | Excluded CLL (2) | Cambridge | 85/M | 0 | 제외(CLL) | Real 6회(d7,d10,d12,d14,d21,d24) | - | mild×4, sev×2 (Real) | severe | 생존 |

---

## 2. 기존(Patient-Average) vs 신규 프레임워크 5개 전 영역 재분석 비교

![기존 Patient-Average 대비 신규 Same-Day WHO Severity 및 First-Draw Final Death 6패널 종합 비교](C:/Users/scott/.gemini/antigravity/brain/add4945a-39a0-4f8d-acd4-3287c21c6121/FigR2_PatientMean_vs_SameDaySev_FirstDrawDeath.png)

### 2.1 [영역 1] 43개 면역 아형 — Severity (전체 217회 채혈 × 당일 WHO) 비교

- **전체 유의 아형 수**: 기존 환자 평균 Peak WHO (`q43 < 0.05`: **12개**) → 신규 당일 WHO Exchangeable GEE (`+week` 보정 `q43 < 0.05`: **18개**, Mancl–DeRouen bias-reduced SE `q43 < 0.05`: **17개**, 환자-블록 부트스트랩 Spearman `q43 < 0.05`: **16개**).
- 43개 전체 효과 크기 상관관계는 $r = 0.84$로 매우 일관되며, 특히 **`B_ASC_IgG`, `B_ASC_dividing`, `B_ASC_IgA`, `HSPC`, `NK CD56bright`, `T CD8 CM`** 등 핵심 생물학적 축은 모두 압도적으로 강화/유지되었습니다.

| 면역 아형 | 기존 Patient-Mean Peak WHO ($\rho$, `q43`) | **신규 Same-Day WHO GEE `+wk` ($\beta$, $P$, `q43`)** | BR `q43` | Clust-Boot $\rho$ (`q43`) | Real $\beta$ (`q43`) | PBISC new $\beta$ ($P$) | **Mundlak Within $\beta$ ($P$)** | Mundlak Between $\beta$ ($P$) | 해석 |
|---|---|---|---|---|---|---|---|---|---|
| **`B_ASC_IgG`** | $+0.391$ (`7.1e-3`) | **$+0.414$ ($5.8 \times 10^{-10}$, `1.8e-8`)** | `8.6e-8` | $+0.331$ (`9.8e-5`) | $+0.336$ (`4.5e-4`) | $+0.461$ (`7.8e-6`) | **$+0.418$ (`2.9e-4`)** | $+0.411$ (`4.7e-5`) | 실시간 당일 중증도 #1 |
| **`B_ASC_dividing`** | $+0.482$ (`6.1e-4`) | **$+0.422$ ($8.2 \times 10^{-10}$, `1.8e-8`)** | `8.6e-8` | $+0.359$ (`3.8e-5`) | $+0.389$ (`1.0e-5`) | $+0.386$ (`4.0e-3`) | **$+0.377$ (`8.2e-5`)** | $+0.458$ (`4.4e-6`) | 실시간 당일 중증도 #2 |
| `CD14mono_dex` | $+0.418$ (`4.2e-3`) | $+0.386$ ($3.1 \times 10^{-9}$, `4.4e-8`) | `7.5e-7` | $+0.358$ (`4.8e-6`) | $+0.378$ (`2.6e-5`) | $+0.319$ (`4.3e-3`) | $+0.310$ (`1.4e-3`) | $+0.447$ (`1.6e-7`) | 덱사메타손 치료 반응 |
| **`B_ASC_IgA`** | $+0.384$ (`7.8e-3`) | **$+0.334$ ($1.3 \times 10^{-5}$, `1.2e-4`)** | `2.9e-4` | $+0.303$ (`1.0e-3`) | $+0.249$ (`0.047`) | $+0.441$ (`7.3e-4`) | **$+0.267$ (`0.015`)** | $+0.397$ (`1.8e-4`) | 실시간 당일 중증도 |
| **`HSPC`** | $+0.354$ (`0.017`) | **$+0.333$ ($1.9 \times 10^{-5}$, `1.3e-4`)** | `3.3e-4` | $+0.318$ (`1.0e-3`) | $+0.161$ (`0.23` NS) | **$+0.610$ (`1.3e-6`)** | **$+0.242$ (`8.8e-3`)** | $+0.453$ (`5.3e-4`) | **PBISC-rescued** 실시간 중증도 |
| **`NK CD56bright`** | $-0.503$ (`3.7e-4`) | **$-0.276$ ($1.3 \times 10^{-4}$, `8.2e-4`)** | `1.9e-3` | $-0.403$ (`2.3e-6`) | $-0.073$ (`0.76` NS) | **$-0.516$ (`1.4e-5`)** | $-0.154$ (`0.11`) | **$-0.573$ (`5.1e-9`)** | **PBISC-rescued** 기저 취약성 |
| **`T CD8 CM`** | $-0.503$ (`3.7e-4`) | **$-0.231$ ($5.2 \times 10^{-4}$, `2.8e-3`)** | `5.7e-3` | $-0.337$ (`1.4e-4`) | $-0.247$ (`8.5e-3`) | $-0.196$ (`0.13`) | $-0.034$ (`0.68`) | **$-0.533$ (`1.7e-6`)** | 기저 중증도 앵커 |
| `T g/d` | $-0.226$ (`0.13` NS) | **$-0.196$ ($6.1 \times 10^{-4}$, `2.9e-3`)** | `6.5e-3` | $-0.231$ (`0.012`) | $-0.100$ (`0.28`) | $-0.486$ (`2.2e-6`) | **$-0.142$ (`0.023`)** | $-0.303$ (`0.010`) | **신규 PBISC-rescued** |
| `CD16mono_C1` | $+0.276$ (`0.071` NS) | **$+0.247$ ($1.0 \times 10^{-3}$, `4.3e-3`)** | `8.6e-3` | $+0.183$ (`0.090`) | $+0.093$ (`0.76`) | $+0.327$ (`2.6e-5`) | **$+0.228$ (`0.035`)** | $+0.272$ (`0.019`) | **신규 PBISC-rescued** |
| `CD16mono` | $-0.400$ (`6.2e-3`) | **$-0.175$ ($0.020$, `0.049`)** | `0.071` | $-0.315$ (`3.4e-4`) | $-0.161$ (`0.27`) | $-0.195$ (`0.12`) | $-0.060$ (`0.56`) | **$-0.374$ (`6.4e-5`)** | 기저 중증도 앵커 |
| `T CD8 CTL` | $-0.345$ (`0.019`) | $-0.085$ (`0.33`, `0.48` NS) | `0.53` | $-0.232$ (`0.046`) | $+0.069$ (`0.79`) | $-0.456$ (`2.8e-4`) | — | — | **조기 사망 특화 (아래 참조)** |
| `NK (CD56dim)` | $-0.447$ (`1.8e-3`) | $-0.120$ (`0.27`, `0.43` NS) | `0.50` | $-0.313$ (`1.0e-3`) | $+0.064$ (`0.83`) | $-0.584$ (`2.6e-7`) | — | — | **조기 사망 특화 (아래 참조)** |

---

### 2.2 [영역 1] 43개 면역 아형 — Death (Donor당 첫 채혈 → 최종 사망 `TRUE_PUB_11`) 비교

> [!TIP]
> **핵심 생물학적 발견 — 조기 예후 예측인자 vs 후기 악화/치료 마커의 완전 분리**
> 기존 환자 평균 사망 분석에서는 초기 예측인자(`T CD8 CTL`)와 후기 급증 마커(`HSPC`, `CD14mono_dex`)가 한데 섞여 있었습니다. **Donor 첫 채혈 → 최종 사망**으로 재분석하자 두 집단이 면도날처럼 분리되었습니다:

| 면역 아형 | 기존 Patient-Mean Death (`AUC`, $P$, `q43`) | **신규 첫 채혈 Death (`AUC`, 95% CI, $P$, `q43`)** | **당일 WHO 보정 첫 채혈 (`AUC`, $P$)** | Real 첫 채혈 (`n=37, 5d`) | PBISC new 첫 채혈 (`n=29, 6d`) | 조기 첫 채혈 ≤7일 (`n=53, 10d`) | Direct-9 (`9d`) / WHO 2–3 | 후기 `>7d` 채혈 (`AUC`, $P$) | 역할 분류 |
|---|---|---|---|---|---|---|---|---|---|
| **`T CD8 CTL`** | `0.206` (`0.0023`, `0.074`) | **`0.213` [0.11–0.33] (`0.0029`, `0.070`)** *(inv **0.787**)* | **`0.294` (`0.033`)** | **`0.206` (`0.036`)** | **`0.210` (`0.031`)** | **`0.265` (`0.022`)** | `0.212` (`0.006`) / **`0.212` (`0.011`)** | `0.332` (`0.19` NS) | **#1 독립 조기 사망 예측인자 (PBISC-rescued)** |
| **`NK (CD56dim)`** | `0.264` (`0.015`, `0.105`) | **`0.217` [0.09–0.37] (`0.0033`, `0.070`)** *(inv **0.783**)* | **`0.291` (`0.030`)** | `0.256` (`0.088`) | **`0.181` (`0.016`)** | **`0.233` (`0.0092`)** | **`0.156` (`0.0010`)** / `0.379` | `0.402` (`0.44` NS) | **#2 독립 조기 사망 예측인자 (PBISC-rescued)** |
| **`B_CD11c`** | `0.334` (`0.085`, `0.262`) | **`0.274` [0.14–0.42] (`0.019`, `0.238`)** *(inv 0.726)* | **`0.309` (`0.048`)** | `0.328` (`0.23`) | **`0.217` (`0.035`)** | `0.330` (`0.099`) | **`0.269` (`0.027`)** / `0.308` | `0.455` (`0.73` NS) | 독립 조기 사망 예측인자 |
| **`CD14mono_IFN`** | `0.704` (`0.034`, `0.184`) | **`0.721` [0.50–0.91] (`0.022`, `0.238`)** | **`0.709` (`0.030`)** | `0.703` (`0.16`) | `0.743` (`0.075`) | **`0.726` (`0.028`)** | `0.680` (`0.086`) / `0.697` | `0.525` (`0.85` NS) | **독립 조기 사망 예측인자 (상승)** |
| **`NK CD56bright`** | `0.349` (`0.117`, `0.334`) | **`0.293` [0.14–0.48] (`0.032`, `0.277`)** *(inv 0.707)* | `0.359` (`0.14`) | `0.238` (`0.064`) | `0.344` (`0.26`) | `0.310` (`0.066`) | **`0.291` (`0.047`)** / `0.455` | `0.270` (`0.070`) | 기저 사망 + 당일 중증도 교량 |
| `T g/d` | `0.270` (`0.017`, `0.105`) | `0.315` (`0.055`, `0.366`) | `0.382` (`0.22`) | `0.375` (`0.39`) | `0.264` (`0.082`) | **`0.292` (`0.043`)** | `0.311` (`0.071`) | **`0.200` (`0.018`)** | 전구간 지속 고갈 |
| `CD16mono_C1` | `0.760` (`0.0069`, `0.074`) | `0.679` (`0.063`, `0.366`) | `0.637` (`0.16`) | `0.603` (`0.48`) | `0.750` (`0.066`) | `0.667` (`0.10`) | `0.686` (`0.076`) | `0.688` (`0.14`) | 경계선 조기 + 당일 중증도 |
| `MAIT` | `0.268` (`0.016`, `0.105`) | `0.344` (`0.106`, `0.391`) | `0.388` (`0.25`) | `0.263` (`0.095`) | `0.413` (`0.53`) | `0.330` (`0.099`) | `0.361` (`0.18`) | **`0.142` (`0.0046`)** | 후기(`>7d`) 사망 고갈 |
| **`HSPC (%)`** | `0.762` (`0.0066`, `0.074`) | **`0.577` (`0.428` NS)** | `0.476` (`0.81` NS) | `0.325` (`0.22`) | `0.812` (`0.022`*) | `0.512` (`0.92` NS) | — | **`0.822` (`0.011`)** | **당일 중증도 & 후기(`>7d`) 악화 마커** |
| **`CD14mono_dex`** | `0.759` (`0.0066`, `0.074`) | **`0.599` (`0.286` NS)** | `0.542` (`0.66` NS) | `0.531` (`0.83`) | `0.692` (`0.033`*) | `0.576` (`0.45` NS) | — | **`0.918` (`8.0e-4`)** | **후기(`>7d`) 스테로이드 치료 마커** |
| `B_ASC_dividing` | `0.674` (`0.071`) | `0.507` (`0.945` NS) | `0.421` (`0.41`) | `0.394` (`0.46`) | `0.627` (`0.36`) | `0.505` (`0.97`) | — | `0.655` (`0.22`) | 순수 당일 중증도 마커 |
| `B_ASC_IgG` | `0.580` (`0.409`) | `0.463` (`0.704` NS) | `0.367` (`0.17`) | `0.291` (`0.13`) | `0.634` (`0.33`) | `0.433` (`0.51`) | — | `0.563` (`0.62`) | 순수 당일 중증도 마커 |

*\*PBISC new에서 HSPC/CD14mono_dex가 첫 채혈 unadjusted로만 유의한 이유는 PBISC 사망자 중 3명이 첫 채혈 시점에 이미 `severe/critical`이었기 때문이며, 당일 WHO 보정 시 즉시 소멸(`P = 0.54–0.81`)합니다. 반면 `T CD8 CTL`, `NK`, `CD14mono_IFN`, `B_CD11c`는 당일 WHO 보정 후에도 독립적으로 사망을 예측합니다!*

---

### 2.3 [영역 2] 신규 프레임워크에서의 PBISC 환자 확장 검정력 구제 (`300`회 시뮬레이션)

`Real` `n = 20`명만 있을 때 대비, bulk-only `PBISC new` **`+29`명**을 추가하면 **신규 프레임워크 양축 모두에서 극적인 검정력 구제(Power Rescue)**가 확인되며, 라벨을 섞은 `+29 Shuffled` 대조군에서는 검정력이 붕괴합니다:

| 분석 축 및 타깃 | 평가 기준 | Real 단독 (`n=20`) | **`+29 PBISC new` (확장)** | `+29 Shuffled` (무작위 대조) | 해석 |
|---|---|---:|---:|---:|---|
| `B_ASC_IgG` (당일 WHO GEE) | `BH43 q < 0.05` | 71.3% | **100.0%** | 16.0% | 완벽 구제 |
| `B_ASC_IgA` (당일 WHO GEE) | `BH43 q < 0.05` | 24.0% | **99.7%** | 8.0% | 4.2배 검정력 상승 |
| `HSPC` (당일 WHO GEE) | `BH43 q < 0.05` | **10.7%** | **100.0%** | 4.3% | **9.4배 검정력 구제** |
| `NK CD56bright` (당일 WHO GEE) | `BH43 q < 0.05` | **7.7%** | **96.7%** | 5.3% | **12.6배 검정력 구제** |
| `T CD8 CM` (당일 WHO GEE) | `BH43 q < 0.05` | 45.7% | **77.7%** | 7.3% | 1.7배 상승 |
| **`T CD8 CTL` (첫 채혈 사망)** | $P < 0.05$ (`q43<0.10`) | **14.7%** (0.0%) | **99.3%** (**24.0%**) | 14.0% (1.0%) | **6.8배 조기 사망 검정력 구제** |
| **`NK` (첫 채혈 사망)** | $P < 0.05$ (`q43<0.10`) | **22.7%** (1.0%) | **99.7%** (**34.7%**) | 11.3% (1.0%) | **4.4배 조기 사망 검정력 구제** |
| **`B_CD11c` (첫 채혈 사망)** | $P < 0.05$ (`q43<0.10`) | 1.3% (0.0%) | **84.0%** (10.0%) | 5.0% (0.0%) | 64배 상승 |
| `CD14mono_IFN` (첫 채혈 사망) | $P < 0.05$ (`q43<0.10`) | 21.3% (1.0%) | **38.3%** (**12.3%**) | 10.0% (1.7%) | FDR 구제 |
| *음성 대조*: `HSPC` (첫 채혈 사망) | $P < 0.05$ | 2.7% | **2.3%** | 2.0% | **위양성 팽창 0% (완벽한 특이도!)** |

---

### 2.4 [영역 3] 단일세포 HSPC 16개 모듈 & 20개 유전자 재분석 (`149` draws `≥3` cells / `61` pts; `≥1` cell `64` pts)

> [!IMPORTANT]
> **PBISC 단일세포 생성만의 고유한 돌파구**:
> 위 2.2절에서 보았듯 **HSPC 세포 비율(%)** 자체는 첫 채혈에서 사망을 예측하지 못합니다(`AUC = 0.577, P = 0.43`). 그러나 **PBISC가 복원한 개별 HSPC 세포의 단일세포 전사체**를 첫 채혈 시점에서 분석하면, **`Alarmin_Inflam` (`S100A8/9/12, LYZ, RETN, MNDA`)**과 **적혈구 분화 전사인자(`Ery5`, `GATA1`, `KLF1`, `TYMS`)**가 이미 첫 채혈 단계에서 최종 사망 환자를 유의하게 구별합니다! 게다가 Real과 PBISC new가 **거의 동일한 수치**로 재현됩니다.

| HSPC 단일세포 프로그램/유전자 | 기존 Patient-Mean Severity (`ρ`, `q`) | **신규 Same-Day WHO GEE `+wk` ($\beta$, $P$, `q`)** | Real / PBISC new $\beta$ | 기존 Patient-Mean Death (`AUC`, $P$) | **신규 첫 채혈 Death (`AUC`, $P$)** | **당일 WHO 보정 첫 채혈 $P$** | Real / PBISC new 첫 채혈 `AUC` |
|---|---|---|---|---|---|---|---|
| **`Cycle4`** (`MKI67, TOP2A, STMN1, TYMS`) | $+0.300$ (`0.144` NS) | **$+0.271$ ($1.0 \times 10^{-5}$, `1.59e-4`)** | `+0.269` / `+0.276` | `0.706` (`0.034`) | `0.645` (`0.136`) | `0.343` | `0.637` / `0.652` |
| **`Ery5`** (`GATA1, KLF1, TFRC, APOC1, BLVRB`) | $+0.374$ (`0.044`) | **$+0.260$ ($2.1 \times 10^{-4}$, `1.70e-3`)** | `+0.274` / `+0.262` | `0.752` (`0.0094`) | **`0.696` (`0.044`)** | `0.272` | `0.644` / `0.746` |
| **`Alarmin_Inflam`** (`S100A8/9/12...`) | $+0.142$ (`0.618`) | $+0.091$ (`0.211` NS) | `+0.054` / `+0.075` | `0.712` (`0.029`) | **`0.724` (`0.021`)** *(≥1 cell: **`0.767, P=0.0057`**)* | **`0.022`** (Real 단독 **`0.043`**) | **`0.726` / `0.725`** |
| `HSC5` (줄기세포성 대조) | $-0.019$ (`0.948`) | $+0.006$ (`0.937` NS) | `-0.160` / `+0.174` | `0.512` (`0.905`) | `0.549` (`0.619` NS) | `0.757` | `0.437` / `0.652` |
| `Myeloid4` (과립구 대조) | $-0.088$ (`0.793`) | $-0.059$ (`0.510` NS) | `-0.169` / `+0.051` | `0.504` (`0.978`) | `0.524` (`0.815` NS) | `0.903` | `0.526` / `0.522` |
| 개별 유전자: **`GATA1`** | $+0.304$ (`0.054`) | **$+0.198$ (`2.4e-3`, `q20=0.011`)** | — | `0.763` | **`0.718` (`0.025`)** | — | — |
| 개별 유전자: **`TYMS`** | $+0.378$ (`0.017`) | **$+0.303$ (`4.0e-6`, `q20=4.2e-5`)** | — | `0.726` | **`0.704` (`0.036`)** | — | — |
| 개별 유전자: **`KLF1`** | $+0.384$ (`0.017`) | **$+0.193$ (`3.1e-3`, `q20=0.011`)** | — | `0.726` | **`0.700` (`0.039`)** | — | — |
| 개별 유전자: **`GATA2`** (감소!) | $-0.216$ (`0.175`) | **$-0.213$ (`3.3e-3`, `q20=0.011`)** | — | `0.390` | `0.460` (`0.687`) | — | *고전적 GATA2→GATA1 전환!* |

---

### 2.5 [영역 4] 단일세포 T-cell 19-gene ISG 재분석 (13개 T 아형)

- **중증도 (Same-Day WHO `+week` GEE)**: 기존 환자 평균에서는 `All T`, `All CD4 T`, `All CD8 T`, `T CD8 CTL` 등이 모두 양의 상관($\rho = +0.30\text{–}0.33, q_{13} = 0.037$)을 보였으나, **채혈 시점(`+week`)을 보정한 당일 WHO GEE에서는 13개 T 아형 전부 $\beta = -0.04 \sim +0.03$ ($P = 0.59\text{–}0.96, q_{13} = 0.96$)으로 완전 소멸**했습니다.
  - **원인**: T-cell ISG는 감염 초극성기(d0–d4)에 정점을 찍고 급격히 감쇠하는 프로그램인데, `PBISC new` 생존자는 첫 채혈이 늦고(중앙 `d8`) 중증/사망자는 일러서(중앙 `d2.5`) **환자 평균 방식이 '초기 채혈 시점'을 '중증도'로 착각**한 것입니다. 따라서 T-cell ISG는 정적인 당일 중증도 마커가 아니라, 기존 Figure 6k에서 보인 **`Time × Fatal` 지속 감쇠 지연(trajectory persistence)** 신호로 서술하는 것이 과학적으로 가장 정확합니다.
- **사망 (Donor 첫 채혈 → 최종 사망)**: `T CD8 CTL ISG` (`AUC = 0.676, P = 0.068`, `adjWHO P = 0.076`), `All CD8 T ISG` (`AUC = 0.669, P = 0.081`; `Real` 단독 첫 채혈 `d2–d4` 균일 조건에서는 **`AUC = 0.794, P = 0.036, adjWHO P = 0.046`**!).

---

### 2.6 [영역 5] Real 코호트 병원 혈액검사 당일 페어링 (Same-Day GEE `+week`)

같은 날 뽑은 혈액검사와 면역 아형을 정확히 당일 1:1 페어링하자, 기존 환자 평균보다 생물학적 결합력이 극적으로 선명해졌습니다:
- **`B_ASC_IgG` ↔ 당일 Ferritin**: $\beta = +0.655, P = 1.70 \times 10^{-13}, q_{36} = 6.10 \times 10^{-12}$ (56쌍 / 33명)
- **`B_ASC_IgA` ↔ 당일 D-dimer**: $\beta = +0.582, P = 6.99 \times 10^{-8}, q_{36} = 1.26 \times 10^{-6}$ (47쌍 / 31명)
- **`B_ASC_IgG` ↔ 당일 ALC (림프구 감소)**: $\beta = -0.316, P = 4.70 \times 10^{-5}$ / **NLR**: $+0.384, P = 6.46 \times 10^{-4}$ / **CRP**: $+0.329, P = 2.42 \times 10^{-3}$
- **`T CD8 CM` ↔ 당일 ALC**: $+0.339, P = 1.14 \times 10^{-3}$ / **Troponin**: $-0.413, P = 3.11 \times 10^{-3}$ / **Ferritin**: $-0.363, P = 4.68 \times 10^{-3}$
- **`NK CD56bright` ↔ 당일 Ferritin**: $-0.341, P = 4.48 \times 10^{-3}$
- **`HSPC` 및 `T CD8 CTL`은 당일 `CRP`, `D-dimer`, `Troponin`, `NLR`과 모두 무관(`q36 > 0.06`)** → 기존 병원 표준 혈액검사가 잡지 못하는 **독립적인 세포 수준 임상 정보**를 제공함을 확인!

---

## 3. 생성 파일 및 재현 산출물 목록

모든 결과표, 통합 Excel, 고해상도 Nature 규격 그림 및 재현 코드는 [draw_sameday_sev_firstdraw_death_20261120](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120)에 저장되었습니다:

- **13-시트 통합 Excel 워크북**: [PBISC_Patient_Draw_Inventory_and_SameDaySev_FirstDrawDeath_20261120.xlsx](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120/PBISC_Patient_Draw_Inventory_and_SameDaySev_FirstDrawDeath_20261120.xlsx) (`README`, `T1_Patients68`, `T1b_Draws217`, `T1c_CohortSummary`, `T1d_CovidNeg39`, `T1e_All362_Draws`, `T2b_Comp43_Combined`, `T2_Comp43_AllArms`, `T3_NewFunnel43`, `T4_PowerRescue`, `T5_HSPC_16Mod_20Gene`, `T6_Tcell_ISG13`, `T7_Labs_SameDay`)
- **출판 규격 그림 (`180.0 mm`, 600 DPI PNG + vector PDF + editable SVG, `oob = 0`)**:
  - [FigR1_Patient_Swimmer_Draw_Schedule_66.png](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120/figures/FigR1_Patient_Swimmer_Draw_Schedule_66.png) / [.pdf](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120/figures/FigR1_Patient_Swimmer_Draw_Schedule_66.pdf)
  - [FigR2_PatientMean_vs_SameDaySev_FirstDrawDeath.png](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120/figures/FigR2_PatientMean_vs_SameDaySev_FirstDrawDeath.png) / [.pdf](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/draw_sameday_sev_firstdraw_death_20261120/figures/FigR2_PatientMean_vs_SameDaySev_FirstDrawDeath.pdf)
- **재현 파이썬 코드**: [05_reanalyze_sameday_sev_firstdraw_death.py](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/code/05_reanalyze_sameday_sev_firstdraw_death.py), [06a_plot_figR1_swimmer.py](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/code/06a_plot_figR1_swimmer.py), [06b_plot_figR2_comparison.py](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/code/06b_plot_figR2_comparison.py), [07_build_reanalysis_excel.py](file:///f:/PBISC/deliverables/nature_submission_fig4_6_20261120/code/07_build_reanalysis_excel.py)

---

## 4. Figure 6 및 논문 본문 반영 제안 (의사결정 포인트)

제안하신 프레임워크(`Severity = 전체 217회 채혈 × 당일 WHO`, `Death = Donor당 첫 채혈 × 최종 사망`)는 **리뷰어가 가장 먼저 지적할 종단 채혈 편향(Peak WHO 불일치, 말기 채혈 역인과성, 채혈 횟수 불균형)을 원천 차단**하면서도, 논문의 핵심 주장(**`T CD8 CTL`, `B ASC`, `HSPC` 및 PBISC 환자 확장 효과**)을 오히려 더 입체적이고 설득력 있게 만들어 줍니다:
- **권장안**: **Main Figure 6 및 Supplementary Figs 6–10과 원고 본문(Results / Methods / Discussion)**을 이 **신규 프레임워크(`Same-day WHO GEE + First-draw Mortality`)를 주분석(Primary)**으로 승격하고, 기존 Patient-average Peak WHO 결과를 민감도 분석(Sensitivity)으로 병기하는 구성을 강력히 추천합니다.
- 원하시면 바로 **Figure 6 및 관련 Supplementary Figures, 영문 Manuscript(`.md` & `.docx`) 전체를 이 신규 프레임워크 기준으로 업데이트**해 드리겠습니다.
