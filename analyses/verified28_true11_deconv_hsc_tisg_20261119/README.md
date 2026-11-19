> **GitHub 배포 안내 (2026-11-19):** 이 배포에는 생성된 PNG·PDF·SVG 및 합본 그림 PDF를 포함하지 않습니다. 아래 그림 목록·로컬 링크·`RUN_RECEIPT.json`은 원본 로컬 산출물의 기록입니다.
>
> **근사 구현의 한계 — 원문보다 우선 적용:** 아래의 BayesPrism-init/upd·BayesPrism Z, CIBERSORTx S/B-mode, MuSiC 및 DWLS 명칭은 원본 코드·표의 식별자입니다. 공식 패키지를 실행한 결과가 아니며 공식 구현과의 동등성도 검증하지 않았습니다. BayesPrism 표기는 결정론적 EM/기대 할당량, CIBERSORTx 표기는 자체 보정+NuSVR, MuSiC/DWLS 표기는 자체 반복 weighted NNLS를 가리킵니다. 이 비교는 해당 근사 기준선에 한정되며 공식 방법의 실패·PBISC의 우월성 또는 공식 사후분포 복원을 입증하지 않습니다. 원문 중 “수식 그대로/충실 구현/결정적 우위”라는 표현에도 이 제한이 적용됩니다. [근사 구현의 상세 한계](APPROXIMATION_LIMITATIONS.md)를 먼저 읽어 주십시오.
>
> 원본 분석 코드·통계값·방법명과 아래 문서는 보존했으며, 이 배포가 원문 전체의 완전 검증 또는 공식 방법 비교의 확증을 의미하지는 않습니다.

---

# 28_Verified_True11_Deconv_HSC_TISG_20261119

Fully verified, 100% standalone-reproducible publication package resolving all 6 verification checkpoints for ESKD COVID-19 Figures 7–8 and Supplementary Figures 1–5.

## Key Deliverables
1. **Verified Patient & 1:1 Sample Metadata Workbook (Points 1 & 2)**:
   - `Table_S1_Verified_Patient_and_Sample_Metadata_20261119.xlsx` (4 sheets: `Mortality_Audit_Point1`, `Patients_68_Verified`, `Sample_1to1_Mapping_362`, `All_107_Donors`)
   - Fixes the published Zenodo `w2_metadata.csv` mortality bug (`C101, C146, C147` verified as survivors against BOTH `Gisby et al., Nat Commun 2022` and `Stephenson et al., Cell Genomics 2025`), locking primary mortality at **`TRUE_PUB_11` (`Real 5 + PBISC new 6 = 11`)** while retaining `Direct 9` and `Zenodo raw 14` sensitivity checks.
2. **Publication Figures (300 DPI `PNG`, `PDF`, `SVG` + 7-page merged PDF)**:
   - `PBISC_Figure7_Main_Verified28_20261119.{png,pdf,svg}`
   - `PBISC_Figure8_Main_Verified28_20261119.{png,pdf,svg}`
   - `PBISC_SuppFig_S1_Cohort_Mortality_Audit_Verified28_20261119.{png,pdf,svg}`
   - `PBISC_SuppFig_S2_CLL_Anchors_Modality_Verified28_20261119.{png,pdf,svg}`
   - `PBISC_SuppFig_S3_Power_IFN_Sensitivity_Verified28_20261119.{png,pdf,svg}`
   - `PBISC_SuppFig_S4_Public_COVID_Replication_Verified28_20261119.{png,pdf,svg}`
   - `PBISC_SuppFig_S5_Deconv_HSC_TISG_Clinical_Verified28_20261119.{png,pdf,svg}`
   - Merged bundle: `PBISC_Verified28_True11_Figures_7_8_S1-S5_20261119.pdf`
3. **Korean Academic Prose Results + Detailed Per-Analysis Provenance**:
   - `Results_Prose_KO_and_Analysis_Provenance_20261119.md`
4. **Full Reproducibility & Integrity Chain**:
   - Stage 1 data freezer: `prepare_verified28_data.py` -> `INPUT_MANIFEST.csv` (23 upstream inputs with SHA-256) + `data/` (43 frozen CSVs) + `OUTPUT_DATA_MANIFEST.csv`
   - Stage 2 standalone renderer (`0 external reads` AST-verified): `render_verified28_figures.py` + `style.py` -> `RUN_RECEIPT.json`

## Standalone Reproduction Commands
```powershell
# Stage 1: Rebuild all frozen CSVs and Table S1 Excel from upstream raw data
py -3.12 "F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\28_Verified_True11_Deconv_HSC_TISG_20261119\prepare_verified28_data.py"

# Stage 2: Render all 7 figures + merged PDF strictly from local data/*.csv (zero external reads)
py -3.12 "F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\28_Verified_True11_Deconv_HSC_TISG_20261119\render_verified28_figures.py"
```
