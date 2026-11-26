# PBISC Figures 4–7 검산·임상 재판정 재현 패키지 (`2026-11-26`)

본 폴더는 [`analyses/nature_submission_fig4_7_final_20261125`](../nature_submission_fig4_7_final_20261125/README.md)에 대한 **2026-11-26 외부 임상 사망 라벨 재판정 및 통계·방법론 정밀 검산 사항**을 `.docx`·`.pdf`·렌더링 이미지 등 대용량 문서 산출물 없이 **독립 재현에 필요한 데이터 테이블·코드·검산 문서만** 경량 구성한 패키지입니다.

- 상세 검산 내역 및 데이터 한계 범위: [`AUDIT_REPORT_KO_20261126.md`](AUDIT_REPORT_KO_20261126.md)
- 외부 사망 판정 수정 명세: [`REQUIRED_UPDATES_FROM_OUTCOME_AUDIT_20261126.md`](REQUIRED_UPDATES_FROM_OUTCOME_AUDIT_20261126.md)

## 수록 구성 (재현 필수 요소 전용)

1. **`tables/` (`22`개 소스/요약 CSV + `00_table_manifest.csv`)**
   - `van der Wijst et al. (2021)` 저자 공식 임상 원본(`Data file S3` / `ClinicalCharacteristics`: `42`명 중 `7`사망)에 맞춰 환자(`Donor_ID`) 단위로 일관되게 정정된 외부 코호트 데이터(`fig7_external4_300_acute_samples.csv`, 생존 기록 `219`명: 생존 `195` / 사망 `24`)와 재산출된 코호트별·무작위효과 메타분석 표(`fig7_external4_per_cohort_results.csv`, `fig7_external4_meta_23subsets.csv`, `fig7_panel_l_summary.csv`) 및 Figures 4–7 전체 `22`종 CSV 테이블을 수록합니다.
2. **`verification_20261126/` (임상 대조 원본 행·HSPC 모듈 점수·검산 결과표)**
   - `clinical_sources/`: 공식 공개 임상 출처 URL·`SHA-256` 매니페스트(`source_manifest.json`) 및 추출된 환자 임상 행(`vdw_clinical_rows.json`, `steph2021_patient_rows.json`)
   - `supp23_hspc62_patient_module_scores.csv`, `supp23_hspc16_definitions.json`, `supp23_hspc16_recomputed.csv`: `Supplementary Fig. 23c`의 중복 제거된 `16`개 고유 HSPC 프로그램(`HSC_broad_q` 포함, `62`명 환자, `11`사망) 재산출 입력 및 결과표
   - `fig5_split_half_recomputed.csv`, `sepsis_scsemi_full_means.csv`, `vdw_donor_outcome_reconciliation.csv` 등 수치 대조 근거 CSV/JSON
3. **`code/` (무의존성 독립 재현 및 검산 스크립트)**
   - `code/01_reproduce_adjudicated_tables.py`: 임상 사망 재판정, `10,000`회 부트스트랩(`seed 20261126`), `DerSimonian–Laird` 메타분석(`BH q23`), `HSPC` 16개 프로그램 검정, `Fig. 5f` 상관계수를 직접 재계산하여 `tables/`와 오차 `< 1e-11` 일치를 검증합니다(`--export-xlsx` 옵션으로 `23`개 시트 `.xlsx` 생성 지원).
   - `code/02_download_clinical_sources.py`: 공개 URL에서 임상 원본 아카이브를 직접 내려받아 `SHA-256` 해시를 검증하는 선택적 스크립트입니다.
   - `code/03_verify_all_fig4_7.py`: `22`개 CSV 해시·차원, Figures 4–7 통계량(`Van der Waerden` 정규화, `GEE`, `Mann–Whitney`, `OLS` 잔차 보정, `DL` 메타분석), 마크다운 원고(`PBISC_Manuscript_Fig4_7_Final_20261125.md`) 및 참고문헌(`1 -> 20`) 순서를 전수 검증합니다.
4. **마크다운 원문 및 검산 기록**
   - `PBISC_Manuscript_Fig4_7_Final_20261125.md`, `AUDIT_REPORT_KO_20261126.md`, `VERIFICATION_REVIEWED_20261126.json`, `VERIFICATION_PACKAGE_CHECKS_20261126.json`
   - *(조판용 `.docx`, 인쇄용 `.pdf`, 렌더링 이미지 바이너리, 대용량 외부 `.zip` 파일은 제외했습니다.)*

## 빠른 실행

```bash
pip install -r ../../requirements.txt
python analyses/nature_submission_fig4_7_reviewed_20261126/code/01_reproduce_adjudicated_tables.py
python analyses/nature_submission_fig4_7_reviewed_20261126/code/03_verify_all_fig4_7.py
```
