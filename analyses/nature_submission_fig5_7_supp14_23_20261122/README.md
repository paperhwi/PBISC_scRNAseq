# PBISC Nature Submission Package — Figures 5–7, Methods (Figures 4–7), and Supplementary Figures 14–23 (`2026-11-22`)

## 1. Methodological Alignment Summary (확정 분석 기준 1:1 반영 요약)

- **Figure 4 (Methods 반영 완료 · 본문 그림 및 보충 그림 `Supplementary Figs. 8–13` 번호 매핑 완비)**:
  - 패혈증(`Ye et al., 2026`) 및 혈액투석 COVID-19(`Gisby et al., 2022; Stephenson et al., 2025`) 각 코호트에서 paired donor `5`명으로 모델을 적응시키고, 평가용 held-out donor의 bulk RNA-seq에서 세포 `4,000`개를 생성해 동일 donor의 실측(`Real`) 세포와 비교.
  - 3가지 참조 환자 조합에서 `T`, `B`, `Monocyte`, `NK` 세포의 질환별 발현 차이를 `limma-voom`으로 추정 (`BH q < 0.05`).
  - 전체 논문(`Figures 1–7`) *Nature Portfolio* 인용 순서 기준 `Figures 1–4`의 보충 그림 `13`종(`Fig. 1: S1–S2`, `Fig. 2: S3–S4`, `Fig. 3: S5–S7`, `Fig. 4: S8–S13`) 뒤를 이어 본 패키지의 `Figures 5–7` 보충 그림 `10`종을 **`Supplementary Figures 14–23`**으로 일괄 정렬 (`SUPPLEMENTARY_NUMBERING_MAP.md` 참조; `--supp-offset` 인자로 언제든 재조정 가능).
- **Figure 5 · 다른 환자 집단 간 비교 (`Figure5_CrossCohort_BulkOnly_Comparison` + `Supplementary Figs. 14–16`)**:
  - 패혈증 `122`개(`488,000` cells) 및 COVID-19 `164`개(`656,000` cells)의 unpaired bulk-only 검체에서 생성한 세포를 별도 환자 집단의 실측 세포(`108`개 `633,072` cells / `187`개 `578,835` cells) 및 `scSemiProfiler`(`34`개 / `102`개)와 비교.
  - 6대 면역 세포 비율은 `환자 × 상태`(`patient × condition`) 단위로 요약 후 환자 단위로 `49,999`회 부트스트랩 재표집.
  - 문헌 보고 질환 신호는 세포형별 유전자 점수(`Ye et al.` 5개 모듈, `Stephenson et al.` 7-ISG 모듈)로 비교.
  - 세포형별 DEG는 양측 **Welch 검정**을 사용하고 **`BH q < 0.05`이면서 `|log₂ fold change| ≥ 1`**인 유전자를 선택 (`AUROC`: 패혈증 `0.78–0.96`, COVID-19 `0.85–0.99`).
- **Figure 6 · 희귀 세포의 빈도와 발현 (`Figure6_Rare_Cell_Frequency_and_Expression` + `Supplementary Figs. 17–18`)**:
  - 고정된 `30-NN` 세포 주석으로 패혈증 **`14`종**, COVID-19 **`13`종**의 희귀 세포 상태를 구분.
  - 환자군(`cases`)에서 각 상태가 상위 세포 계통(`parent lineage`)에서 차지한 비율을 산출하되, **`HSC` 및 `HSPC`는 전체 세포 수(`all cells`)를 분모**로 사용 (`PBISC` `13/14` 및 `12/13` 재현 vs `scSemiProfiler` `7/14` 및 `8/13` 0-cell 소실).
  - 검체별로 같은 상태의 세포 발현을 합산(`log₁ₚ(CP10k)`)하여 **다른 환자(`patient_i ≠ patient_j`)**의 같은 상태 실측 세포와 **Pearson 상관(`r`)**을 구하고(`중앙값 r = 0.73–0.92` 및 `0.81–0.89`, `Real–Real` 환자 간 상한 근접/일치) 대표 마커 발현 히트맵을 비교.
- **Figure 7 · 채혈 당일 중증도와 마지막 채혈 기준 사망 비교 (`Figure7_SameDay_Severity_and_LastDraw_Mortality` + `Supplementary Figs. 19–23`)**:
  - 비-CLL 혈액투석 급성 COVID-19 환자 **`66`명**의 적격 급성기 검체 **`276`건**(`138 Real` + `138 PBISC`, 당일 중복 `0`건) 사용.
  - **중증도**: 채혈 당일 WHO 등급과 세포 비율(`% PBMC`)의 연관성을 환자별 총가중치를 같게 두어(`w_ij = 1 / n_i`, 환자당 합 `1.0`, 총합 `66.0`) 평가 → `BH43 q < 0.05` 유의 세포형 **`11/43` (`Real`) → `18/43` (`Combined raw %`, `rank-INT` 시 `20/43`)**로 확대.
  - **사망**: 환자별 **마지막 적격 검체 1건**(`66`건 = `24 Real` + `42 PBISC`, 최대 채혈일 동률 `0`건)으로 **생존 `55`명 vs 사망 `11`명**의 세포 비율을 비교(`Mann–Whitney AUC`) → `BH43 q < 0.05` 유의 예측인자 **`0/43` (`Real`) → `2/43` (`Combined raw %`: `T g/d AUC = 0.183, q = 0.0393` 및 `HSPC AUC = 0.800, q = 0.0393`; `rank-INT` 시 `CD16mono_C1` 포함 `3/43`)** 발굴.
  - **외부 4개 연구 메타분석**: `Stephenson 2021`, `COMBAT 2022`, `van der Wijst 2021`, `Liu 2021` (총 `300`개 급성기 검체 / `231`명, 사망 추적 코호트 생존 `197`명 vs 사망 `22`명)에서 **동일 기준**으로 분석 후 대응 연구 `k ≥ 2`인 세포형을 **DerSimonian–Laird 무작위효과 모형**으로 합성하여 재현 확인.

---

## 2. Directory Contents

- `PBISC_Nature_Manuscript_Fig5_7_and_SuppFig14_23_20261122.md` & `.docx`: Nature Portfolio 형식 원고 (`Results Figs. 5–7`, `Main Legends Figs. 5–7`, `Discussion Figs. 5–7`, `Methods Figs. 4–7`, `Supplementary Legends 14–23`).
- `SUPPLEMENTARY_NUMBERING_MAP.md`: 전체 논문 `Main Figures 1–7` 대응 `Supplementary Figures 1–23` 인용 순서 매핑표.
- `PBISC_SourceData_and_Verification_Fig5_7_20261122.xlsx`: 전 패널 핵심 통계/검산 통합 엑셀 워크북 (20개 시트).
- `VERIFICATION_REPORT_20261122.json`: 독립 재계산, 보충 그림 인용 단조 증가 순서(`14..23`), 그림 규격(`180.0 mm`, `min_pt ≥ 5.0`, `oob = 0`, `600 DPI`), `.docx` OOXML 무결성 전수 검산 리포트 (`81 / 81 PASS`, `FAIL = 0`).
- `figures/`: Main Figures `5–7` 및 Supplementary Figures `14–23` (`13`개 그림 × `PNG [600 DPI] / PDF [450 DPI] / SVG [editable Arial]` = 총 `39`개 파일).
- `data/`: 재현성 100% 원천 및 요약 CSV 데이터.
- `code/`: 데이터 산출(`01_...`), 본문/보충 그림 생성(`02_...`, `03_...`, `nature_style.py`), 문서 빌드 및 독립 검산(`04_...`) 스크립트 일체.
