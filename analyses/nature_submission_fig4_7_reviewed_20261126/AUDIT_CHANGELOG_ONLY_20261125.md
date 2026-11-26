> 2026-11-26 정정: 이 문서는 11월 25일 변경 이력입니다. 최신 검산 범위와 미완료 항목은 AUDIT_REPORT_KO_20261126.md가 우선합니다. 아래 정수 count 오기는 원본 표에 맞춰 수정했습니다.

# 검산을 통해 변경·정정된 사항 요약서 (`2026-11-25`)

> **문서 목적**: 기존 초안(`2026-11-23` Figure 4 원고 및 `2026-11-24` Figures 5–7 원고·표·피규어)과 중간 감사본(`audit_20261125` / `corrected_20261125`)에 대해 원시 데이터·원본 코드·마스터 PDF 벡터·PubMed 서지 대조 검산을 수행한 결과, **실제로 수정·반영된 항목만** 발췌해 정리한 변경 내역서입니다.  
> **반영 대상 패키지**: `analyses/nature_submission_fig4_7_final_20261125/`

---

## 1. 한눈에 보는 변경 요약표 (총 20개 변경 항목)

| 구분 | 항목 | 변경 전 (기존 문서/파일) | 변경 후 (최종 반영본) | 검산 근거 |
| :--- | :--- | :--- | :--- | :--- |
| **Fig. 4 / Supp 13** | 1. `Supp Fig. 13` 수치 누락 보완 | COVID-19 `r = 0.80–0.87`만 기재 | **패혈증(`PBISC r = 0.58–0.72, CCC = 0.54–0.63` vs `scSemi r = 0.62–0.66, CCC = 0.58–0.62`) + COVID-19(`r = 0.80–0.87, CCC = 0.79–0.85` vs `0.41–0.68, 0.37–0.64`) 8개 전 패널** 완비 | `PBISC-Diffusion_20261123.pdf` Pages 42–43 텍스트·벡터 대조 |
| **Fig. 4** | 2. `iLISI` 기준선 및 참조셋 해석 명시 | `Real–Real` 및 3-reference 평균만 약술 | `Real–Real`(`1.77`/`1.76`)은 동일 방문 내 **경험적 반분 참조치**(독립 donor 상한 아님), 3-reference 변동은 **지지 표본 선택 민감도**(독립 생물학적 반복 아님)임을 명시 | 마스터 PDF Pages 9–10 레전드 정합화 |
| **Fig. 5a,j–m** | 3. 비대응(`unpaired`) 비교 범위 정밀화 | 서로 완전히 다른 donor 간 비교로 서술 | **검체 수준 비대응(`specimen-level unpaired`)**으로 정정: 패혈증은 공통 donor `0`명이나, COVID-19 `Real B`와 `PBISC` 조성 집단 간 **종단 환자 `7`명(`C113, C127, C145, C168, C169, C20, C63`)** 공유 명시 및 표 플래그 추가 | `fig5_donor_condition_297units_composition.csv` 교집합 검산 |
| **Fig. 5b–e,l** | 4. `Real A` / `Real B` 모집단 범위 명시 | 전체 실측 코호트를 양분한 것처럼 서술 | 전체 코호트(`633,072` / `578,835` 세포)에서 추출한 **고정 donor-disjoint 평가 부분표본(`58,195 / 59,457` 및 `57,613 / 56,171` 세포)**임을 명시 | 원시 `.h5ad` 및 평가 부분표본 세포 수 대조 |
| **Fig. 5f** | 5. 패혈증 `scSemiProfiler` 유전자 상관계수 | `r = 0.937` 단독 표기 (`corrected`는 `0.936`으로만 변경) | 마스터 그림 인쇄값 **`r = 0.937`**(전체 소스 패널 `0.9365`) 유지 + 보관 표시용 요약표 기준 **`0.936`**(`0.9363`) 괄호 병기 | `sepsis_scsemi.h5ad` vs `sepsis_gene_counts_displayed.csv` 대조 |
| **Fig. 5m** | 6. 조성 부트스트랩 비유의(`ns`) 해석 | `"without significant bias"` (무편향 입증 어조) | **유의한 조성 왜곡이 검출되지 않음(`no statistically significant distortion detected`)**으로 정정 (비유의 ≠ 동등성 증명) | 통계적 가설검정 해석 기준 정합화 |
| **Fig. 5p–t** | 7. Welch pseudobulk DEG 최소 세포 수 및 집계 | `≥ 5 cells/sample` 기재 | 원 코드(`lib2.py`) 기준 **`≥ 20 cells per visit and cell type`**, 방문 합산 → `log2(CPM + 1)`(`Real` 전체 유전자 합 vs 생성 공유 패널 합) → COVID 환자×조건 평균(`COVID_EXCL` 5개 방문 제외)으로 정정 | `final_real_pbisc_scsemi_20261111/code/lib2.py` 대조 |
| **Fig. 5p–t** | 8. Welch DEG 참조 집단·`BH` 범위·반분 횟수 | `Real B` 대조 및 동일 패널 내 `BH`, `50`회 반분으로만 기재 | **보관된 전체 `Real` pseudobulk(`93/7`, `86/7`; `39/38`, `37/34`, `37/36`, `39/37`)** 대조 명시; `Real` genome-wide `BH` vs 생성 공유 패널 `BH` 공개(공유 패널 재보정 시 이진 DEG 수 차이 `0`건 확인); COVID 교차 donor 미보정 한계 명시; `real_splithalf_AUROC`(`50`회) vs `split_cov` 신뢰도(`200`회) 구분 | `f1_deg_metrics.py` 및 `deg_*.csv.gz` 전수 재계산 |
| **Fig. 6a,b** | 9. `in_primary_scope` 및 case 분모 공개 | 생성된 `656,000` COVID 세포 전량 사용처럼 기재 | 동결 `in_primary_scope`가 **`C31` 3회 방문(`12,000` 세포) 제외** → 조건 필터 전 `644,000`개, pooled case 분모는 **패혈증 `272,000` / COVID-19 `560,000`개**임을 명시 | `ceiling_cell_labels.csv.gz` 및 `in_primary_scope` 검산 |
| **Fig. 6a,b,e,f** | 10. 희귀 상태 백분율 정수 count 정밀화 및 문구 수정 | 일부 수치가 PDF 막대 좌표 역산 반올림값(`13.22%`, `34.17%`, `90.54%`, `24.46%`, `25.18%` 등); `"negative by construction"` | **원시 정수 분자/분모 기반 정확한 백분율(`13.24%`, `34.29%`, `90.47%`, `24.52%`, `25.26%` 등)**로 본문·표 전수 교체; 점선 설명을 **`"negative in these data after baseline subtraction"`**으로 정정 | `fig6_rare_parent_shares_integer_counts.csv` 원시 집계 대조 |
| **Fig. 7 파일** | 11. `Figure 7` `.png`–`.pdf` 파일 버전 동기화 | `.png`는 최종본이나 `.pdf`는 구버전(`CTL–NK–ISG` 등) 잔존 | 최종 확정 `Figure 7` 원본으로부터 **`.png`와 단일 페이지 `.pdf`를 100% 동일 내용·비율로 동기화** | 시각 대조 및 종횡비·스트림 검증(`ratio_diff < 10⁻⁴`) |
| **Fig. 7a,b** | 12. `Measured` `37`명 내부 구성 명시 | `37`명 전체가 `Bulk + scRNA-seq`인 것처럼 읽힘 | **`10`명 `scRNA-seq-only`(`47` draws, `1` death) + `27`명 `shared with bulk`(`91` measured draws, `4` deaths)**로 본문·레전드·Methods 및 `66`명 표(`cohort_design_stratum` 열)에 분리 명시 | `276` 급성기 검체 및 `66`명 첫 채혈 메타데이터 대조 |
| **Fig. 7** | 13. 순위 기반 정규화 공식 정정 | `Blom` 공식 `Φ⁻¹((rank − 3/8)/(N + 1/4))` 기재 | 실제 저장 점수와 `< 10⁻¹²`로 일치하는 **`Van der Waerden` 공식 `Φ⁻¹((rank − 0.5)/N)` (동점 평균 순위)**로 정정 (`Blom`과 최대 `~0.088` 차이) | `eskd_217_*` 및 `eskd_66_*` `__z217`/`__zfd` 열 재계산 |
| **Fig. 7c** | 14. `Measured` 단독 당일 WHO `GEE` 유의 개수 | `5`개(주요 앵커 개수)로 기재 | 전체 43개 세포형 중 `Measured` 단독 `GEE +wk` `BH q43 < 0.05`는 **`9`개**(`Combined` `18`개와 `7`개 공유)로 정정 | `statsmodels.GEE` 43개 전 세포형 독립 재적합 |
| **Fig. 7d** | 15. 첫 채혈 CI 부트스트랩 횟수 및 `WHO-adj. P` 계산법 | `2,000`회 부트스트랩; `"exact stratified permutation / logistic regression"` | **`3,000`회 층화 부트스트랩**으로 정정; `WHO-adj. P`는 **OLS 회귀 잔차의 양측 Mann–Whitney 검정**(단일 분율 `.033/.030`은 수치형 WHO + Combined source 항; 복합점수 `.013`은 범주형 WHO 단독, source 항 없음)임을 명시하고 명목 탐색적 보정임을 고지 | 원 코드 및 OLS 잔차 독립 재계산 일치 확인 |
| **Fig. 7d, Discussion** | 16. `50 × 5` 교차검증 성격 및 사망 분석 한계 명시 | `"prespecified q43 < 0.10"`, `"leak-free 50 × 5 CV"` | **`exploratory q43 < 0.10` 스크린** 및 전처리만 폴드 내부로 제한한 **`fixed-feature 50 × 5` 내부 교차검증**(변수 선택은 CV 이전 수행)으로 정정; Discussion에 `11` 사망·`q43 = 0.070`·비투석 코호트 미재현·독립 검증 필요성 단락 신설 | `recompute_cv.py` 및 통계적 검증 범위 명시 |
| **Fig. 7 / Supp 23** | 17. `HSPC Ery5` 분석 시점 분리 및 `q16` 병기 | 첫 채혈 문맥 인근에서 `AUC = 0.752, P = 0.0094`만 인용 | 첫 채혈 결과가 아닌 **별도의 전 방문 환자 평균(`all-visit patient-mean`, `62`명/`11`사망) 분석**임을 분리 명시하고 **`q16 = 0.150`** 추가 | `supp10_hspc_16modules_summary.csv` 대조 |
| **References** | 18. 참고문헌(`1–20`) PubMed 원본 전수 교정 | Ref `1`(`Ye, P.`), Ref `2`(*eLife* 오기재), Ref `3`(`Zenodo` 혼동), Ref `8`·`13` 제목 부정확 | **PubMed/DOI 원 출판 기록 대조 전수 정정**: `Ye, Q. et al. Nat. Immunol. 27, 150–165 (2026)`, `Gisby, J. S. et al. Nat. Commun. 13, 7775 (2022)`(`Zenodo 6497251` 원 출처), `Stephenson, E. et al. Cell Genom. 5, 100918 (2025)`, `van der Wijst`·`Reyes` 원제 복원, `TMM`/`voom`/`limma`(Ref `5–7`) 추가 및 `1 → 20` 순번 정렬 | PubMed PMID `41339491`, `36522333`, `40532694`, `34103408` 검증 |
| **시각·조판 회귀 복원** | 19. `corrected_20261125`의 피규어·문장 회귀 오류 제거 | `Fig. 6` 재작성 시 UMAP 컬러 스와치 8개 소실; `Fig. 7` 재작성 시 `Panel l` Forest Plot 소실; 치환 중복 어구 발생 | **마스터 원본 `Fig. 6` 보존**(벡터 오차 `< 0.042 pt`, 스와치 8개 유지), **사용자 확정 `Fig. 7` 디자인 보존**(`Panel l` 첫 채혈 vs 마지막 채혈 Forest Plot 유지), 중복 어구 교열 제거 | PDF 벡터 좌표 정밀 계측 및 원고 전문 교열 |
| **PDF·GitHub 배포** | 20. PDF 특수기호 깨짐 및 GitHub 푸시 차단 해결 | 기본 PDF 폰트에서 `≤`, `α` 등이 `■`로 출력; `.gitignore` 누락 및 SVG 내장 JPEG Base64 문자열 오탐(`GH013`) | PDF 빌더에 **TrueType 유니코드 폰트 임베딩** 적용(전 페이지 기호 정상 출력); `.gitignore` 화이트리스트 등록 및 SVG 내장 래스터 **무손실 PNG 재인코딩 + 48자 줄바꿈**으로 44개 전 파일 푸시 성공 | 페이지 렌더링 시각 QA 및 `git ls-files` 확인 |

---

## 2. 세부 수치 변경 대조 (원고 본문 기준)

### A. Figure 6 희귀 면역 상태 백분율 (PDF 좌표 역산값 → 원시 정수 분자/분모 정확한 백분율)
> 결론(`PBISC`: 패혈증 `13/14`, COVID-19 `12/13` 출현; `scSemiProfiler`: 패혈증 `7/14`, COVID-19 `8/13` 완전 `0` 세포)은 변함없으며, 소수점 둘째 자리 반올림 오차만 원시 정수 집계(`fig6_rare_parent_shares_integer_counts.csv`) 기준으로 정밀화했습니다.

- **패혈증 (`Sepsis`, case 분모 `PBISC = 272,000` cells)**:
  - `B · MZB1`: `PBISC` `13.22%` → **`13.24%`** (`Real B 3.15%` 유지)
  - `Mono · MKI67`: `Real B` `1.27%` (`1.274%`, `92 / 7,219`), `PBISC` `0.019%` (`9 / 47,650`)
  - `pDC`: `PBISC` `9.50%` → **`9.53%`** (`115 / 1,207`), `Real B` `34.17%` → **`34.29%`** (`48 / 140`)
  - `cDC`: `PBISC` `90.54%` → **`90.47%`** (`1,092 / 1,207`), `Real B` `65.75%` → **`65.71%`** (`92 / 140`)
  - `Macrophage`: `Real B` `2.53%` → **`2.54%`** (`188 / 7,407`), `PBISC` `1.18%` 유지
  - `HSC`: `PBISC` `0.19%` → **`0.192%`** (`522 / 272,000`), `Real B` **`0.048%`** (`20 / 41,457`)
  - `Mono HLA-DRA`: `PBISC` `10.83%` → **`10.82%`**, `Real B` `8.30%` → **`8.31%`**
  - `NKT CCL4hi`: `PBISC` `21.88%` → **`21.91%`**, `Real B` `29.69%` → **`29.66%`**
  - `ILC`: `Real B` `0.76%` → **`0.77%`** (`PBISC 0%`)
- **COVID-19 (case 분모 `PBISC = 560,000` cells, `C31` 3회 방문 제외 후)**:
  - `CD14 Mono · dex`: `Real B` `6.64%` → **`6.66%`** (`PBISC 0.12%` 유지)
  - `CD14 Mono · IFN`: `PBISC` `24.46%` → **`24.52%`** (`15,007 / 61,202`), `Real B` `15.21%` → **`15.24%`**
  - `ASC IgG`: `PBISC` `5.71%` → **`5.72%`**, `Real B` `19.57%` → **`19.62%`**
  - `Cycling ASC`: `PBISC` `25.18%` → **`25.26%`** (`15,151 / 59,971`), `Real B` `34.84%` → **`34.86%`**
  - `DC1`: `Real B` `3.59%` → **`3.58%`** (`PBISC 2.34%` 유지)
  - `pDC`: `PBISC` `19.38%` → **`19.37%`**, `Real B` `25.42%` → **`25.44%`**
  - `CD56bright NK`: `Real B` `12.41%` → **`12.42%`** (`PBISC 5.18%` 유지)

### B. Figure 7 및 보충 분석 수치·표기 변경점
- `Measured` 단독 당일 WHO 중증도 `GEE +wk` 유의 세포형(`BH q43 < 0.05`): `5`개 → **`9`개** (`Combined` `18`개는 동일, 두 목록의 교집합 `7`개).
- `HSPC Ery5` 세포 내 발현 사망 연관(`Supplementary Fig. 23b,c`): `AUC = 0.752, P = 0.0094` → 전 방문 환자 평균(`62`명, `11`사망) 분석임을 명시하고 **`q16 = 0.150`** 추가.
- 첫 채혈 `AUC` 95% 신뢰구간 부트스트랩 반복 수: `2,000`회 → **`3,000`회**.
- 그 외 모든 핵심 통계치(`217` 급성기 채혈/`66`명, WHO `126/58/21/12`, `CD8 CTL AUC = 0.213 [q43 = 0.070, WHO-adj. P = 0.033]`, `NK AUC = 0.217 [q43 = 0.070, WHO-adj. P = 0.030]`, `ISG – CTL – NK AUC = 0.841 [P = 4.3 × 10⁻⁴, WHO-adj. P = 0.013]`, `50 × 5 CV` `0.694 → 0.846 / 0.858`, 외부 4개 코호트 `300` 채혈/`231`명·생존 기록 `219`명/`22`사망 메타분석 전 수치)는 독립 재계산에서 `10⁻¹²` 이내로 **원래 값과 완전히 일치하여 그대로 유지**되었습니다.

---

## 3. 참고문헌 변경 전·후 상세 대조

1. **Ref 1 (패혈증 원 코호트)**
   - 변경 전: `Ye, P. et al. Integrated single-cell and bulk blood transcriptomics across adult and paediatric bacterial sepsis. GEO accessions GSE279448 and GSE279452 (2026).`
   - 변경 후: `Ye, Q. et al. Single-cell multi-omic landscape reveals anatomical-specific immune features in adult and pediatric sepsis. Nat. Immunol. 27, 150–165 (2026). https://doi.org/10.1038/s41590-025-02345-x`
2. **Ref 2 (ESKD COVID-19 bulk RNA-seq / Zenodo `6497251` 원 출처)**
   - 변경 전: `Gisby, J. et al. Multi-omics longitudinal profiling of acute COVID-19 in patients with end-stage kidney disease receiving haemodialysis. eLife 11, e74879 (2022).`
   - 변경 후: `Gisby, J. S. et al. Multi-omics identify falling LRRC15 as a COVID-19 severity marker and persistent pro-thrombotic signals in convalescence. Nat. Commun. 13, 7775 (2022). https://doi.org/10.1038/s41467-022-35454-4`
3. **Ref 3 (ESKD COVID-19 종단 scRNA-seq 원 코호트)**
   - 변경 전: `Stephenson, E. et al. Single-cell and bulk peripheral blood transcriptomic atlas of SARS-CoV-2 infection in patients with end-stage kidney disease. Zenodo https://doi.org/10.5281/zenodo.6497251 (2025).`
   - 변경 후: `Stephenson, E. et al. Temporal multi-omics analysis of COVID-19 in end-stage kidney disease. Cell Genom. 5, 100918 (2025). https://doi.org/10.1016/j.xgen.2025.100918`
4. **Ref 5–7 (`Figure 4` 통합에 따른 `TMM` / `limma-voom` 방법론 원전 추가)**
   - 추가: `Robinson, M. D. & Oshlack, A. Genome Biol. 11, R25 (2010)` (`Ref 5`), `Law, C. W. et al. Genome Biol. 15, R29 (2014)` (`Ref 6`), `Ritchie, M. E. et al. Nucleic Acids Res. 43, e47 (2015)` (`Ref 7`).
5. **Ref 11 (`van der Wijst et al., 2021`, 구 Ref 8)**
   - 변경 전: `Longitudinal single-cell epitope and RNA-sequencing reveals the immunological landscape of COVID-19.`
   - 변경 후: `Type I interferon autoantibodies are associated with systemic immune alterations in patients with COVID-19. Sci. Transl. Med. 13, eabh2624 (2021). https://doi.org/10.1126/scitranslmed.abh2624`
6. **Ref 16 (`Reyes et al., 2021`, 구 Ref 13)**
   - 변경 전: `...induces suppressive myelopoiesis from human hematopoietic progenitors in vitro.`
   - 변경 후: `...induces suppressive myeloid cell production from hematopoietic progenitors in vitro. Sci. Transl. Med. 13, eabe9599 (2021). https://doi.org/10.1126/scitranslmed.abe9599`
