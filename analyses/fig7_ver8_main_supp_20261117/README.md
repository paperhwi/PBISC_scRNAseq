# Figure 7 ver8 — 본 그림과 보충 그림 1·2 (2026-11-17)

ESKD COVID-19 코호트에서 Real scRNA-seq 37명(CLL 환자 C116·C139 제외)과 bulk RNA-seq만 있는 새 환자 29명(`PBISC new`)을 합쳐 43개 세포형의 중증도·사망 연관을 선별하고 검증한 코드입니다. 최종 그림은 `PBISC_Figure7_Main_ver8_20261117`, `PBISC_SuppFig_S1_CLL_Anchors_Modality_20261117`, `PBISC_SuppFig_S2_Screen_Power_Sensitivity_20261117`이고, 결과·범례 원고는 `docs/Results_and_Legends_ver8.md`입니다.

## 파일

| 파일 | 내용 |
| --- | --- |
| `code/run_ver8_pipeline.py` | 통계 전체: CLL 분자 스크린, 43개 세포형 스캔과 5단계 선별 퍼널(43 → 16 → 15 → 12 → 8), 14가지 민감도 변형, 신뢰도 스코어카드, bootstrap 95% CI, HSPC 사망 강건성, 검정력 곡선, 27명 공유 환자의 모달리티 스케일 비교, NK 분해. `data/*.csv` 13개를 씁니다 |
| `figure/make_main_and_supp_v8.py`, `figure/style.py` | `data/*.csv`로 본 그림(a–m), 보충 그림 1(a–o), 보충 그림 2(a–e)와 패널별 PNG/PDF/SVG를 그립니다 |
| `code/verify_independent_recompute.py` | 검산용 독립 재구현. 파이프라인 코드를 쓰지 않고 환자 표, 43개 스캔, BH, 퍼널(centre 보정 포함), Wilcoxon, HSPC 분석, 검정력을 다시 계산해 파이프라인 출력과 비교합니다 |
| `PBISC_Figure7_ver8_verification_20261117.ipynb` | 검산 노트북. 원고 문장·범례의 수치를 `data/*.csv`와 대조하고(FAIL은 원고 오류를 뜻함) 마지막에 독립 재구현을 실행합니다. 저장소 규칙에 따라 출력은 지웠습니다 |
| `docs/Results_and_Legends_ver8.md` | 결과·범례 원고. 아래 검산에서 찾은 오류 6건을 고친 판입니다 |

## 실행

1. 입력(저장소에 없음): `F:\PBISC\eskd_covid_clinical_assoc_20261115\results\sample_level_composition_metadata.csv`, `sample_level_fine_label_pct.csv`, `F:\PBISC\paper_claim_analysis\covid_sample_patient_id_reconciliation.csv`, `F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\12_CLL_excluded_severity_death_discovery\data\cll_molecular_verification_bulk.csv`.
2. 두 스크립트의 `ROOT`와 `V8_DIR`(`F:\PBISC\...`)을 자기 환경에 맞게 바꿉니다. `sys.path`의 `F:\PBISC\.analysis_runtime`은 원래 PC의 패키지 경로이며, pandas·scipy·statsmodels·matplotlib만 있으면 됩니다.
3. `python run_ver8_pipeline.py` (약 9분) → `python make_main_and_supp_v8.py`.
4. 검산: 노트북의 `V8_DATA`(ver8 `data/` 폴더), `PBISC_ROOT`(위 입력 폴더들의 상위 폴더), `V8_INDEP`(독립 재구현 스크립트 경로) 환경 변수를 지정하고 실행합니다.

## 분석 정의

- **환자 값:** COVID+, primary scope, 세포 200개 이상인 샘플만 씁니다. RBC·Platelet을 빼고 100%로 다시 맞춘 뒤 환자별 평균을 냅니다.
- **합친 코호트 (`Real + PBISC new`, 66명, 사망 13명):** source 안에서 rank-based inverse-normal score를 매긴 뒤 합칩니다. Real과 PBISC에 모두 있는 27명은 Real 값만 씁니다.
- **검정:** 중증도는 Spearman ρ(peak WHO 0–3), 사망은 Mann–Whitney AUC입니다. BH 보정은 43개 세포형 안에서, CLL 포함/제외 × 집단별로 따로 합니다.
- **퍼널:** (2) 합친 코호트에서 severity q < 0.05 또는 death q < 0.10 → (3) 세 집단의 방향이 같고 |ρ| ≥ 0.10 또는 |AUC − 0.5| ≥ 0.05 → (4) centre/wave 보정 후 P < 0.05이며 CD14 mono dex가 아님 → (5) 같은 endpoint가 Real 단독으로는 FDR 유의하지 않음.

## 검산 (2026-11-17)

- **재실행:** 원본 입력으로 `run_ver8_pipeline.py`를 다시 돌린 결과, 출력 CSV 13개가 저장본과 최대 1.8×10⁻¹⁵ 차이로 같았습니다.
- **독립 재구현:** 퍼널 43/16/15/12/8, S5 후보 8개, 43개 세포형의 ρ·AUC·q, centre 보정 P, 27명 Wilcoxon P, HSPC 분석(Cambridge only, WHO 2–3, 중증도 층화 permutation)이 모두 일치했습니다(최대 차이 10⁻¹⁶). 다른 seed로 돌린 n = 20 검정력도 몬테카를로 오차 안에서 같았습니다.
- **원고에서 고친 오류 6건** (`docs/`에 반영):
  1. paired Wilcoxon P를 고쳤습니다. B ASC IgA 1.5×10⁻⁸ → **7.5×10⁻⁸**(27쌍 중 25쌍이 PBISC > Real), B ASC dividing 1.1×10⁻⁴ → **9.6×10⁻⁵**, HSPC 2.7×10⁻⁴ → **3.0×10⁻⁴**.
  2. B ASC IgA raw % 범위를 고쳤습니다. PBISC new는 **1.3–21.9%**(mild 1.3–8.1%), critical Real은 **0.4–4.5%**입니다.
  3. HSPC의 "WHO 2–3 안에서 43개 세포형 중 유일하게 사망 예측"은 틀렸습니다. 43개 중 6개가 P < 0.05입니다(T CD4 Helper IFN stim, CD14 mono dex·IFN, T CD4 Helper, ASDC 포함). 유일한 것은 **step 4를 통과한 12개 안에서**입니다.
  4. 보충 그림 2a 범례: Int. mono(Real 단독 severity q = 0.0498)는 그림 코드의 `LBL`에 없어 그려지지 않습니다. 범례에 이를 적었습니다.
- **해석 주의:** HSPC의 공유 27명 Real–PBISC 일치도는 ρ +0.20(P 0.31)입니다. 추적하지 않기로 한 후보 3개의 기준(ρ 0.07–0.19, P ≥ 0.34)과 거의 같은 수준이고, Fig. 7d에서도 빨간 칩입니다. B ASC IgG/IgA도 노란 칩(ρ 0.26–0.28)입니다. "다섯 집단이 신뢰도 기준을 충족"이라는 표현은 이 점을 감안해 다듬어야 합니다.
