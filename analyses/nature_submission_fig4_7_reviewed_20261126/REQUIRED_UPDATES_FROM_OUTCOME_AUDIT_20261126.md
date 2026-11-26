# 외부 코호트 임상 사망 라벨 정정에 따른 패키지 수정 반영 명세서 (`2026-11-26`)

> **근거 문서**: `Donor-identified outcome reconciliation audit note` (*Donor-identified outcome reconciliation and supplemental keys*)  
> **대상 패키지**: `analyses/nature_submission_fig4_7_final_20261125`  
> **핵심 원인 요약**: 외부 4개 검증 코호트 중 **`van der Wijst 2021` (`n = 42`)**의 사망(`death`) 컬럼이 기존 `D121` 테이블에서 `5`명 사망으로 잘못 코딩되어 있었으며, 원 논문 `Data file S3` 및 CELLxGENE 원본(`death == 'yes'`) 기준 실제 사망자는 **`7`명**(`Rep_C_1029`, `Rep_C_1037`, `Rep_C_1047`, `Rep_C_1050`, `Rep_C_1130`, `Rep_C_1154`, `Rep_C_1161`; 기존과 `10`명 불일치 — 위양성 `4`명 해제 + 위음성 `6`명 복원, `Rep_C_1130` `1`명 유지)임.

---

## 0. 수정이 필요 없는 범위 (영향 없음 확인)

아래 항목들은 외부 코호트 사망 라벨을 전혀 사용하지 않으므로 **수치·그래픽·해석이 100% 그대로 유지**됩니다.
- **Figures 4, 5, 6 및 Supplementary Figs. 8–19 전체** (모든 그래픽·표·본문·캡션·Methods 불변)
- **Figure 7 내부 ESKD 코호트 전체 (`Panels a–j`, `Supplementary Figs. 20a,d–o, 23–26`)**: `66`명 환자(`37` Measured + `29` PBISC), `217`건 급성 채혈, 사망 `11`명, 당일 `WHO` 중증도 `GEE` (`18/43`, `BH q43 < 0.05`), 첫 채혈 사망 `AUC` (`CD8 CTL 0.213`, `NK 0.217`, 각 `q43 = 0.070`), 복합 점수 `ISG – CTL – NK` (`AUC = 0.841, WHO-adj. P = 0.013`), `50×5` 내부 교차검증(`0.694 → 0.846 / 0.858`) 전부 불변
- **Figure 7k 하단 외부 당일 중증도 메타분석 표 (`fig7_panel_k_summary.csv`)**: 사망 라벨과 무관한 채혈 당일 서열 척도 `Spearman ρ`이므로 `Dividing ASC*` (`+0.26`), `HSPC` (`+0.40`), `MAIT` (`−0.29`), `CD56-bright NK` (`−0.28`), `γδ T` (`−0.27`), `T CD8 CM*` (`−0.29`) 및 `Ext. q23` 전수 불변

---

## 1. 데이터 테이블 (`tables/` 및 통합 `.xlsx`) 수정 항목 (`4`개 CSV + 워크북 재생성)

1. **`tables/fig7_external4_300_acute_samples.csv`**
   - `cohort == "van der Wijst 2021"`인 `75`개 행(`42`명)에 대해 `death` 값을 실제 사망 `7`명(`Rep_C_1029`, `Rep_C_1037`, `Rep_C_1047`, `Rep_C_1050`, `Rep_C_1130`, `Rep_C_1154`, `Rep_C_1161`)만 `1.0`, 나머지 `35`명은 `0.0`으로 정정 (위양성이었던 `Rep_C_1053`, `Rep_C_1109`, `Rep_C_1113`, `Rep_C_1131`은 `1.0 → 0.0`).
   - 결과 기록 보유 환자 수: `219`명 유지 (`197` 생존 / `22` 사망 → **`195` 생존 / `24` 사망**: `Stephenson 73/7`, `COMBAT 75/7`, `van der Wijst 42/7`, `Liu 29/3`).

2. **`tables/fig7_external4_per_cohort_results.csv`**
   - `cohort == "van der Wijst 2021"` 행들의 `n_surv` (`37 → 35`), `n_dead` (`5 → 7`), 첫 채혈 사망 `auc_death_first` / `se_auc_first` / `p_death_first` 및 마지막 채혈 사망 `auc_death` / `se_auc` / `p_death` / `ci_lo_auc` / `ci_hi_auc` 재계산 반영 (`Stephenson`, `COMBAT`, `Liu` 및 모든 `rho_w` 열은 불변).

3. **`tables/fig7_external4_meta_23subsets.csv`**
   - `n_surv_ext` / `n_dead_ext` 및 `23`개 세포군(`+ B_ASC_total`)의 첫 채혈·마지막 채혈 `DerSimonian–Laird` 메타분석 추정치(`AUC`, `95% CI`, `P`, `I²`, `BH q23`) 갱신.

4. **`tables/fig7_panel_l_summary.csv`**
   - 열 이름 `vanDerWijst_42_5_AUC_1st` → **`vanDerWijst_42_7_AUC_1st`** (또는 기존 열명 유지 시 메타데이터 명시) 및 아래 정정 수치 반영 (필요 시 마지막 채혈 유의 아형인 `DC2` 행 추가 여부 선택):

| 세포군 (`display_label`) | `vdW` 1st AUC (`기존 → 정정`) | Pooled First AUC (`기존 → 정정`) | Ext. `q23` First (`기존 → 정정`) | Pooled Last AUC (`기존 → 정정`) | Ext. `q23` Last (`기존 → 정정`) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`CD8 CTL*`** | `0.465 → 0.290` | `0.477 → 0.411` | `0.926 → 0.343` | `0.423 → 0.423` | `0.482 → 0.449` |
| **`NK`** | `0.459 → 0.392` | `0.496 → 0.476` | `0.965 → 0.834` | `0.552 → 0.504` | `0.926 → 0.966` |
| **`γδ T`** | `0.249 → 0.302` | `0.217 → 0.231` | `0.051 → 0.088` | `0.255 → 0.220` | `0.314 → 0.039†` |
| **`MAIT`** | `0.311 → 0.259` | `0.221 → 0.211` | `0.0064† → 0.0024†` | `0.270 → 0.213` | `0.314 → 0.017†` |
| **`CD8 naïve T`** | `0.211 → 0.469` | `0.266 → 0.323` | `3.4e-5† → 0.0076†` | `0.369 → 0.331` | `0.427 → 0.017†` |
| **`CD16 mono`** | `0.659 → 0.365` | `0.482 → 0.412` | `0.938 → 0.343` | `0.320 → 0.394` | **`0.020† → 0.217` (유의성 상실)** |
| **`HSPC`** *(참고행)* | `0.497 → 0.776` | `0.569 → 0.660` | `0.623 → 0.088` | `0.599 → 0.640` | `0.482 → 0.209` |
| **`DC2`** *(Last 유의군)* | `(0.202)` | `(0.286)` | `(0.100)` | **`0.276`** | **`0.024†`** |

5. **`PBISC_SourceData_and_Summary_Tables_Fig4_7_20261125.xlsx` 및 `tables/00_table_manifest.csv`**
   - 위 `4`개 CSV 갱신 후 `code/01_export_figures_and_tables.py`로 매니페스트 해시 및 `23`-시트 통합 엑셀 파일 재생성.

---

## 2. 피규어 그래픽 수정 항목 (`figures/` 및 `supplementary_figures/`)

1. **본문 `Figure 7` (`figures/Figure7_ESKD_COVID19_Severity_and_FirstDraw_Mortality.{pdf,png,svg}`)**
   - **`Panel k` 상단 워크플로 텍스트/박스**: `HSPC` 외부 첫 채혈 사망 `AUC` 표기를 `0.57` → **`0.66`** (`q23 = 0.088`; 코호트별 `vdW` `0.50 → 0.78`)으로 수정 (`severity ρ = +0.40`은 동일).
   - **`Panel l` 부제(subtitle)**: `219 outcome-documented patients · 197 survived / 22 died` → **`219 outcome-documented patients · 195 survived / 24 died`** (`vdW 42/5` → `vdW 42/7`).
   - **`Panel l` 표/플롯 수치 및 유의성 기호(`†`)**:
     - 위 표의 정정된 `vdW` 1st-draw `AUC`, `Pooled AUC (First / Last)`, `Ext. q23 (First / Last)` 반영.
     - 마지막 채혈(`Last`) 유의 표시(`†`, `q23 < 0.05`): `CD16 mono`에서 `†` 제거(`0.39, q23 = 0.217`), **`γδ T` (`0.22, q23 = 0.039†`), `MAIT` (`0.21, q23 = 0.017†`), `CD8 naïve T` (`0.33, q23 = 0.017†`)**(및 포함 시 `DC2` `0.28, q23 = 0.024†`)에 `†` 부여.

2. **보충 피규어 (`supplementary_figures/` 내 `Supplementary Figs. 20, 21, 22, 27`)**
   - **`Supplementary Fig. 20b,c`**: 범례(key) 및 캡션에서 (i) `CLL` 제외 후 추정치의 소스별 색상 `95% bootstrap CI` 가로선과 (ii) 제외 전·후 추정치를 잇는 회색 연결선(`shift`) 두 요소를 각각 명확히 구분 표기.
   - **`Supplementary Fig. 21b,c`**: 축/패널 라벨에서 `prospective power` 표현을 제거하고 **`observed resampling detection rate`**(관측된 코호트 내 비복원 리샘플링 검출률)로 표기 통일.
   - **`Supplementary Fig. 22`**: 범례에 내부 `ESKD` 코호트 사망 `11`명, 외부 peak-grade 분석 코호트(`n = 219`) 사망 `23`명(`6 + 7 + 7 + 3`)을 명시하고, 내부 ESKD `Table S1` 출처가 외부 `Stephenson 2021`(*Nat. Med.*)과 별개인 `Stephenson et al. 2025`(*Cell Genom.*)임을 명시.
   - **`Supplementary Fig. 27`**: 외부 마지막 채혈(`last-draw`) 사망 후보 검증 열에서 `van der Wijst` `7`명 사망 정정에 따른 후보별 메타분석 추정치 갱신 확인.

---

## 3. 원고 (`PBISC_Manuscript_Fig4_7_Final_20261125.{md,docx,pdf}`) 및 `README.md` 수정 항목

1. **`Results` — Figure 7 마지막 문단**
   - 외부 결과 기록 환자 분모: `219 outcome patients, 22 deaths` → **`219 outcome patients (195 survivors, 24 deaths)`**.
   - 첫 채혈 사망 비교 수치 갱신:
     - `CD8 CTL*`: `pooled First AUC = 0.48` → **`0.41` (`q23 = 0.343`)**; `NK`: `0.50` → **`0.48` (`q23 = 0.834`)** (비투석 코호트에서 비유의하다는 핵심 결론 동일).
     - 비투석 초기 유의 감소군: **`MAIT` (`First AUC = 0.21, q23 = 0.0024†`)** 및 **`CD8 naïve T` (`First AUC = 0.32, q23 = 0.0076†`)** (및 명목 감소 `γδ T` `0.23, P = 0.013, q23 = 0.088`; `HSPC` 첫 채혈 `AUC = 0.66, P = 0.017, q23 = 0.088`로 다중보정 비유의).
   - **마지막 채혈(`Last draw`) 서술 정정 (핵심 내용 변경점)**:
     - 기존 *"양측 코호트가 말기 `CD16 mono` 유의 감소(`Last AUC = 0.32, q23 = 0.020†`)로 수렴했다"*는 문장을 수정하여, 외부 코호트 마지막 채혈에서는 **`CD8 naïve T` (`0.33, q23 = 0.017†`), `MAIT` (`0.21, q23 = 0.017†`), `DC2` (`0.28, q23 = 0.024†`), `γδ T` (`0.22, q23 = 0.039†`)** 감소가 `q23 < 0.05`를 통과하며 지속되고, 총 **`CD16 mono`** 감소는 같은 방향의 경향성(**`Last AUC = 0.394 [0.280–0.507], P = 0.066, q23 = 0.217`**; `ESKD Last AUC = 0.269`)에 머물렀음을 정확히 기술.

2. **`Discussion` — 두 번째 문단 마지막 문장**
   - *"prior to shared late `CD16` monocyte depletion"* 부분을 정정하여, 비투석 코호트에서는 초기부터 말기까지 `CD8 naïve T` 및 선천성 T 세포(`MAIT`, `γδ T`)·`DC2` 고갈이 우세하게 지속되는 반면(`CD16` 단핵구 말기 감소는 다중보정 전 경향 수준, `AUC = 0.394, q23 = 0.217`), 요독증 기저 미감작 T 세포 수축이 선행된 `ESKD` 환자에서는 초기 세포독성-인터페론 취약성(`CD8 CTL`·`NK` 감소 + `ISG` 상승)이 두드러진다는 대조 구도로 정밀화.

3. **`Main Figures and Figure Legends` — `Fig. 7k, l` 범례**
   - `k`: 상단 예시 `HSPC` 첫 채혈 사망 `AUC` `0.57` → **`0.66`**.
   - `l`: `219 outcome-documented patients (197 survivors, 22 deaths: Steph. 73/7, COMBAT 75/7, vdW 42/5, Liu 29/3)` → **`(195 survivors, 24 deaths: Steph. 73/7, COMBAT 75/7, vdW 42/7, Liu 29/3)`**.

4. **`Methods` — `External four-cohort COVID-19 validation` 섹션**
   - `van der Wijst et al., 2021` 사망자 수를 `5 deaths` → **`7 deaths`** (원 논문 `Data file S3` 및 CELLxGENE `death` 필드 일치 명시), 총계 `197 survivors, 22 deaths` → **`195 survivors, 24 deaths`**로 수정.
   - `Stephenson et al., 2021`에서 런던(`AP3`) 사망을 포함한 결과 기록 `73`명(`7`명 사망) 기준과 `Supplementary Fig. 22`의 peak-grade `73`명(`6`명 사망, 외부 합계 `23`명 사망) 기준 차이를 1문장으로 명시.

5. **`Supplementary Figures and Legends` (`Supplementary Figs. 20, 21, 22, 27`) 및 패키지 `README.md`**
   - `Supp. Fig. 20`: 제외 후 추정치 `95% CI` 가로선과 제외 전·후 연결 회색선 구분 명시.
   - `Supp. Fig. 21`: `power curves` 대신 `observed resampling detection-rate curves`로 용어 통일.
   - `Supp. Fig. 22`: 외부 peak-grade 코호트 `219`명·`23`명 사망(`6 + 7 + 7 + 3`) 및 ESKD `Table S1` = `Stephenson et al. 2025` 출처 명시.
   - 패키지 `README.md`: `Figure 7` 설명 줄의 외부 사망자 수 `22사망` → **`24사망`** 수정.

---

## 4. 자동 재현·검산 스크립트 (`code/`) 수정 항목

1. **`code/01_export_figures_and_tables.py`**
   - 외부 샘플 테이블 내보내기 시 `van der Wijst 2021`의 원 논문 검증 `7`명 사망 ID 매핑을 적용하고 `per_cohort`, `meta_23subsets`, `panel_l_summary` CSV 및 `.xlsx`를 일관되게 재생성하도록 반영.
2. **`code/03_verify_all_fig4_7.py`**
   - `fig7_external_300_draws_231_patients_219_outcomes` 검증 조건의 `n_out_deaths == 22` → **`n_out_deaths == 24`** (`vdW == 7`)로 갱신.
   - `fig7_external_dl_meta_recomputed` 검증 기준값을 정정 후 수치(`CD8 naive T` 1st `0.323`, `MAIT` 1st `0.211`, `CD16 mono` last `0.394`)로 갱신 후 `02_build_manuscript_docx_pdf.py` 및 `03_verify_all_fig4_7.py` 재실행.
