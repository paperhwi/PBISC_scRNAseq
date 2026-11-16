# Figure 7 — 희귀 세포 집단과 bulk-only 환자 확장 (2026-11-16)

ESKD COVID-19 코호트에서 Real scRNA-seq 39명에 bulk RNA-seq만 있는 새 환자 29명(PBISC로 단일세포 생성)을 더해 B ASC dividing, NK CD56bright, HSPC의 중증도·사망 연관을 분석한 코드입니다. 최종 Figure 7 ver5(`PBISC_Figure7_ASC_NK_HSPC_ver5_20261116.pdf`)의 1쪽 본 그림(a–m)과 2쪽 보충 그림(a–e)을 만드는 단계만 모았습니다.

## 그림과 코드의 대응

**본 그림 (1쪽, `figure/make_panels_v5.py`)**

| 패널 | 내용 | 수치 출처 |
| --- | --- | --- |
| a | 연구 설계 | — |
| b | 세 집단 × Real / PBISC new / 합친 코호트의 Spearman ρ (peak WHO severity), bootstrap 95% CI | `s08` → `gain_severity_rho.csv` |
| c | 사망 AUC, bootstrap 95% CI | `s08` → `gain_death_auc.csv` |
| d | 80% 검정력에 필요한 Real 환자 수 | `s07` 검정력 곡선 → `s08` → `gain_power_curves.csv` |
| e–g / h–j / k–m | B ASC dividing / NK CD56bright / HSPC: severity, 사망, 검정력 곡선 | `s02` 환자 표 + `s08` |

**보충 그림 (2쪽, `figure/export_fig7_v3.py`가 만드는 ver3 패널)**

| 패널 | `panels_fig7_v3/` 파일 | 함수 | 수치 출처 |
| --- | --- | --- | --- |
| a | `SFig6_c_composition_by_severity` | `panels.comp_stack` | `s04` → `patient_level_covid_positive.csv` |
| b | `Fig7_b_severity_dotgrid` | `panels_v3.dotgrid_v3` | `s07` → `composition_association_v3.csv` |
| c | `SFig13_a_extension_power_all_features` | `panels_v3.ext_power_all` | `s07` → `extension_power_v3.csv` |
| d | `SFig7_b_severity_effect_concordance` | `panels.effect_concordance` | `s03` → `assoc_patient_level.csv` |
| e | `SFig6_e_HSPC_death_ROC` | `panels.roc_hspc` | `s04` → `roc_hspc_fatal*.csv` |

## 실행 순서

| 단계 | 파일 (원래 이름) | 출력 |
| --- | --- | --- |
| 1 | `code/s01_prep_cells.py` (`work/prep.py`) | `cache/`: Real·PBISC 세포 표와 pseudobulk |
| 2 | `code/s02_compose.py` (`work/compose.py`) | `out/sample_level_composition_metadata.csv`, `out/patient_level_composition_metadata.csv` |
| 3 | `code/s03_assoc.py` (`work/assoc.py`) | `out/assoc_patient_level.csv` 등 |
| 4 | `code/s04_compute_fig_data.py` (`fig/compute_fig_data.py`) | `data/`: 보충 그림 a, e 입력 |
| 5 | `code/s05_v3_draw0_labels.py` | `cache/`: 두 번째 adapter(draw 0) 라벨 (`s07`의 robustness용) |
| 6 | `code/s06_v3_deconvolution.py` | `data_v3/deconvolution_*`: 같은 bulk의 NNLS deconvolution 비교 |
| 7 | `code/s07_v3_composition.py` | `data_v3/composition_association_v3.csv`, `extension_power_v3.csv` (1,000회 subsampling, severity shuffle 대조군), PPI, robustness |
| 8 | `code/s08_gain_analysis.py` (`gain/gain_analysis.py`) | `data/gain_*.csv`: 본 그림 b–d 수치 |
| 9 | `code/s09_shared_patient_check.py` | `data/shared_patient_check.csv`: 두 플랫폼에 모두 있는 28명에서 Real과 PBISC 비교 |
| 그림 | `figure/make_panels_v5.py` → `figure/build_pptx_v5.js` | `panels_out/` 패널별 PNG/PDF/SVG → 178 × 258 mm PPTX |
| 그림 | `figure/export_fig7_v3.py` (+ `make_figures.py`, `panels*.py`) | `panels_fig7_v3/`: 보충 그림 패널 |

`PBISC_Figure7_ASC_NK_HSPC_ver5.ipynb`는 8단계, 본 그림 패널, 9단계를 차례로 실행합니다. 저장소 규칙에 따라 출력은 지웠습니다.

## 경로와 입력

- `s01`, `s02`, `s05`–`s07`의 `/mnt/user-data/uploads/PBISC/`는 분석 PC의 `F:\PBISC\`입니다. 다음 입력이 필요하고 저장소에는 없습니다.
  - `eskd_covid_clinical_assoc_20261115/data/`: Real obs·pseudobulk·희귀 세포 count
  - `deliverables/annotated_h5ad_20261112/`: PBISC 주석 h5ad와 `label_mapping_covid.csv`
  - `bulkonly_lora5_20261006/analysis/bulk_public/zenodo6497251_imperial_eskd/`: bulk count와 메타데이터
  - `NewLora/`: draw 0 생성본, ceiling 라벨, holdout 라벨 정밀도
- 작업 폴더는 환경 변수로 지정합니다.
  - `RUN`: `s01`–`s03`의 작업 폴더. `cache/`와 `out/` 하위 폴더가 생깁니다.
  - `FIG`: `s04`, `s06`, `s07`, `export_fig7_v3.py`를 실행하는 폴더. `data/`와 `data_v3/`가 생깁니다.
  - `GAIN`: `s08`을 실행하는 폴더. `data/`가 생깁니다.

  기본값은 원래 분석 환경 경로(`/home/claude/...`)이니 자기 환경에 맞게 바꾸세요. `s01`–`s03`은 `RUN` 안에서, 그 밖의 스크립트는 각자 출력 폴더 안에서 실행합니다.
- `s04`의 뒷부분(CD14 mono IFN DEG 표)은 이 그림에 쓰지 않습니다. 이 부분은 별도 DEG 출력(`out/deg/`)이 있어야 하므로, 그 출력이 없으면 보충 그림 입력을 쓴 뒤에 오류로 멈춥니다.
- `build_pptx_v5.js`는 `pptxgenjs`가 필요합니다. 테마 색 적용(`apply_theme.js`)은 `PPTX_SKILL`이 있을 때만 하며, 없어도 패널 배치는 같습니다.

## 분석 정의

- **환자 값:** COVID+ primary scope 샘플 중 세포가 200개 이상인 샘플을 쓰고, 환자별 % of PBMCs 평균을 냈습니다.
- **합친 코호트:** source 안에서 rank-based inverse-normal transform 했습니다. 두 플랫폼에 모두 있는 28명은 Real 값만 씁니다(같은 환자를 두 번 세지 않음).
- **검정:** severity는 Spearman, 사망은 Mann–Whitney AUC입니다. CI는 bootstrap 2,000회, BH 보정은 세 집단에 걸쳐 했습니다.
- **검정력:** Real n명(10–35)을 1,000회 비복원 추출해 P < 0.05 비율을 셉니다. 방향은 B ASC dividing과 NK CD56bright는 전체 Real과 같은 방향, HSPC는 양방향입니다. shuffled 대조군은 매 반복마다 PBISC new 29명의 severity 라벨만 섞습니다.

## 해석 주의

HSPC는 두 플랫폼에 모두 있는 28명에서 Real ρ −0.04, PBISC ρ 0.49로 엇갈립니다(`s09`). 환자별 값 일치도도 낮고(ρ 0.14), 같은 bulk의 deconvolution도 비슷한 추세를 보입니다(`s06`). 그래서 HSPC 결과는 "기존 보고의 방향을 재현"으로만 서술하세요. B ASC dividing과 NK CD56bright는 28명에서도 Real과 PBISC의 방향이 같습니다.
