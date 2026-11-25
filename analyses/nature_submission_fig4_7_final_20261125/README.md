# PBISC Unified Nature Portfolio Submission Package — Figures 4–7 (`2026-11-25`)

본 폴더는 **`PBISC-Diffusion_20261123.pdf` 마스터 Figure 4–6 (`Supplementary Figs. 8–19`)**, **최종 확정 `Figure 7` (`Supplementary Figs. 20–27`)**, 그리고 **`audit_20261125`의 모든 통계·방법론·참고문헌 정밀 검산 사항**을 하나로 통합한 최종 제출 패키지입니다.

## 폴더 구성 (4대 산출물 통합)

1. **`figures/` — 최종 본문 Figure 4–7 파일 (`11`개 파일)**
   - `Figure4_Paired_Reconstruction_and_Disease_Signal_Recovery.{pdf,png,svg}`: 동일 환자(`paired`) `5`-donor `LoRA` 재구성, 비통합 `50-D PRE-PCA` `iLISI`, `TMM`/`limma-voom` 질병 신호 복원 (`PBISC-Diffusion_20261123.pdf` Page 8 원본 벡터 보존)
   - `Figure5_BulkOnly_Composition_and_Disease_Programs.{pdf,png,svg}`: 패혈증(`122`건)·ESKD COVID-19(`164`건) `bulk-only` 생성, `pre-PCA iLISI`, 환자×상태 조성 부트스트랩, 문헌 질병 프로그램(`14`종) 및 `Welch` pseudobulk DEG 충실도 (Page 11 원본 벡터 보존)
   - `Figure6_Assigned_Rare_Populations_and_Molecular_Profiles.{pdf,png,svg}`: 패혈증 `14`종·COVID-19 `13`종 희귀 면역 상태 출현율(정수 count 기반), `UMAP` 위치화(8개 희귀 상태 컬러 범례 보존), 타 환자 간 `2,000-HVG` 기저선 차감 유사도(`270` + `625` = `895` 프로파일) 및 마커 히트맵 (Page 14 원본 벡터 보존)
   - `Figure7_ESKD_COVID19_Severity_and_FirstDraw_Mortality.{pdf,png}`: 비-CLL 혈액투석 COVID-19 `66`명(`10 scRNA-only + 27 shared + 29 PBISC bulk-only`)·`217` 급성기 채혈의 당일 `WHO` 중증도(`18/43`, `BH q43 < 0.05`), 첫 채혈 사망 예측(`CD8 CTL`, `NK`, `ISG – CTL – NK` 및 `50×5` 내부 교차검증), 외부 `4`개 코호트(`300`검체/`231`명; 생존 기록 `219`명/`22`사망) 메타분석 (`.png`와 단일 페이지 `.pdf` `100%` 동기화 완료)

2. **`tables/` 및 `PBISC_SourceData_and_Summary_Tables_Fig4_7_20261125.xlsx` — 최종 표 파일 (`22`개 CSV + 통합 `23`-시트 Excel 워크북)**
   - `tables/00_table_manifest.csv`: 전체 `22`개 CSV의 행·열 수, 바이트 크기, `SHA256` 해시 및 Excel 시트 매핑 목록
   - **Figure 4 (`5`종)**: `fig4_paired_composition_6types.csv`, `fig4_paired_ilisi_and_gene_r.csv`, `fig4_paired_limma_deg_metrics_32pairs.csv`, `fig4_supp12_d0_volcano_deg_counts.csv`, `fig4_supp13_native_effect_agreement.csv`
   - **Figure 5 (`4`종)**: `fig5_donor_condition_297units_composition.csv` (`7`명 COVID 공통 donor 명시 플래그 포함), `fig5_composition_summary_by_group.csv`, `fig5_published_disease_signals_42rows.csv`, `fig5_welch_deg_fidelity_16rows.csv`
   - **Figure 6 (`3`종)**: `fig6_rare_parent_shares_integer_counts.csv` (`272,000` 패혈증 / `560,000` COVID case 정수 분모·분자 완비), `fig6_interpatient_similarity_895profiles.csv`, `fig6_similarity_summary_24groups.csv`
   - **Figure 7 (`10`종)**: `fig7_eskd_217_acute_draws_43subsets.csv`, `fig7_eskd_66_first_draws_43subsets.csv` (`10 scRNA-only / 27 shared / 29 PBISC` 구분열 포함), `fig7_eskd_63_first_draws_CTL_NK_ISG.csv`, `fig7_eskd_43subsets_severity_and_mortality_screen.csv`, `fig7_firstdraw_composite_and_50x5_cv.csv`, `fig7_external4_300_acute_samples.csv`, `fig7_external4_per_cohort_results.csv`, `fig7_external4_meta_23subsets.csv`, `fig7_panel_k_summary.csv`, `fig7_panel_l_summary.csv`

3. **통합 원고 파일 (실제 Figure 4–7 본문 삽입 완료)**
   - `PBISC_Manuscript_Fig4_7_Final_20261125.docx`: *Nature Portfolio* 서식(실제 `Figure 4–7` 고해상도 이미지 4장 직접 삽입, 고유 OOXML 위첨자 참고문헌 `1–20`)
   - `PBISC_Manuscript_Fig4_7_Final_20261125.pdf`: 인쇄·열람용 A4 PDF (실제 `Figure 4–7` 4장 삽입)
   - `PBISC_Manuscript_Fig4_7_Final_20261125.md`: 원문 Markdown 소스 (**Results**, **Main Figures & Legends**, **Discussion**, **Methods**, **Supplementary Figs. 8–27 Overview**, **PubMed 검증 References `1–20`**)

4. **`code/` 및 `VERIFICATION_FINAL_FIG4_7_20261125.json` — 독립 재현 및 전수 검산 코드**
   - `code/01_export_figures_and_tables.py`: 마스터 PDF/원시 표 추출 및 `23`-시트 `.xlsx` 생성 (단독 클론 환경 호환)
   - `code/02_build_manuscript_docx_pdf.py`: `.md` + `figures/`로부터 `.docx` 및 `.pdf` 자동 빌드 및 OOXML 검증
   - `code/03_verify_all_fig4_7.py`: 외부 절대경로 의존성 없이 동봉된 `tables/`, `figures/`, 원고만으로 `Van der Waerden` 정규화, `statsmodels GEE`, `Mann–Whitney AUC`, `OLS 잔차 WHO-adj. P`, `DerSimonian–Laird` 메타분석, 정수 희귀세포 분율, 참고문헌 `1→20` 순서 및 이미지 해시를 독립 재계산·검증하여 `VERIFICATION_FINAL_FIG4_7_20261125.json` 생성

## 빠른 재현 실행법

```bash
python code/01_export_figures_and_tables.py
python code/02_build_manuscript_docx_pdf.py
python code/03_verify_all_fig4_7.py
```
