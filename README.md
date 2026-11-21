# PBISC scRNA-seq analysis code

PBISC가 생성한 단일세포 데이터의 검증·비교 분석 코드를 모은 저장소입니다. 로컬에 남아 있던 분석 스크립트와 노트북을 분석 시기와 주제별로 `analyses/` 아래에 정리했습니다. 주된 대상 코호트는 Ye 2026 sepsis와 Imperial ESKD/COVID입니다.

최신 **Nature Portfolio 투고 패키지(Main Figures 4–7 & Supplementary Figures 1–11, 2026-11-21)**는 [nature_submission_fig4_7_supp1_11_20261121](analyses/nature_submission_fig4_7_supp1_11_20261121/README.md)에 있습니다. 비중복 4-Main-Figure 구조(**Fig. 4** bulk-only 분자·Harmony 충실도, **Fig. 5** 계통별 질병 DEG, **Fig. 6** Sepsis+COVID 희소 집단 재현·희소 세포 UMAP·공식 R 4.3.2 11개 디컨볼루션·3,124개 단일세포 HSPC 20-유전자 스위치, **Fig. 7** 임상 발견: 1:1 환자/데이터셋 메타데이터 매칭, 세포 수/임상 CBC vs 비율 검산, T-cell 19-유전자 ISG 심층 분석)와 보충 그림 1–11(PNG/PDF/SVG 45개 파일, `0 OOB`), 원고 초안(`.md`/`.docx`), 11-시트 감사 Excel 워크북 및 재현 코드를 모두 포함합니다.

2026-11-19 `TRUE_PUB_11` 분석의 코드·동결 CSV·메타데이터 워크북·원문 결과는 [verified28_true11_deconv_hsc_tisg_20261119](analyses/verified28_true11_deconv_hsc_tisg_20261119/README.md)에 추가했습니다. 이 배포에는 생성된 그림 파일을 포함하지 않습니다. 디컨볼루션 비교는 공식 패키지와 동등성이 검증되지 않은 자체 근사 구현이므로 [근사 구현의 한계](analyses/verified28_true11_deconv_hsc_tisg_20261119/APPROXIMATION_LIMITATIONS.md)를 먼저 읽어야 합니다.

현재 통일 재분석에 기반한 **Figure 7·8 및 보충 Figure S1–S5**는 [figures7_8_unified_20261118](analyses/figures7_8_unified_20261118/README.md)에 있습니다. 희귀 HSPC–중증도와 HSPC KIAA0391–중증도 후보, ASC·HSPC의 FDR rescue를 구분하며, PDF/SVG/PNG·원 계산표·환자 가명 plotting 값·재현 코드·독립 검산을 포함합니다. 이전 Figure 7 버전은 당시 방법과 해석을 확인할 수 있는 보관본으로 유지합니다.

## 분석 폴더

| 폴더 | 내용 |
| --- | --- |
| [`nature_submission_fig4_7_supp1_11_20261121`](analyses/nature_submission_fig4_7_supp1_11_20261121/README.md) | **최신 Nature 투고 패키지 (Main Figs 4–7 + SuppFigs 1–11)**: 비중복 계통 DEG(Fig 5), Sepsis+COVID 희소 집단·UMAP·공식 R 4.3.2 11개 디컨볼루션·HSPC 20-유전자(Fig 6), ESKD COVID-19 임상 3대 축—1:1 메타데이터 매칭(255 sepsis / 362 COVID), 세포 수/CBC vs 비율 검산(`r=0.987/0.985`), T-cell 19-ISG 심층 분석(Fig 7), 원고(`.md`/`.docx`), 11-시트 Excel 감사표 및 45개 그림 파일 |
| [`verified28_true11_deconv_hsc_tisg_20261119`](analyses/verified28_true11_deconv_hsc_tisg_20261119/README.md) | TRUE_PUB_11 분석의 원본 코드·CSV·메타데이터·결과 문서 보존본. 자체 디컨볼루션 근사의 한계를 명시했으며 생성 그림은 제외 |
| [`figures7_8_unified_20261118`](analyses/figures7_8_unified_20261118/README.md) | 현재 Figure 7·8 + 보충 S1–S5: 고정 BH43, 공유 5,762개 유전자 전수 screen, HSPC KIAA0391 고정 10M 순열, 효과·CI·보정 범위와 fidelity, 전체 결과표 및 portable 그림 재현 |
| [`bulkonly_lora5_20261006`](analyses/bulkonly_lora5_20261006/README.md) | LoRA-5 K=1,000 생성 데이터, 원본 scRNA-seq 재처리, 통합·조성·희귀 세포·질병 DEG·환자별 검증 및 후속 분석 |
| [`bulkonly_lora5_k4000_20261006`](analyses/bulkonly_lora5_k4000_20261006/README.md) | K=4,000 Harmony 분석과 K=1,000 비교 |
| [`bulkonly_lora5_20261007`](analyses/bulkonly_lora5_20261007/README.md) | sepsis K=1,000/K=4,000 클러스터 비교 |
| [`new_inference_depth256_K4000_floor1e-4_20261021`](analyses/new_inference_depth256_K4000_floor1e-4_20261021/README.md) | floor readout 데이터 다운로드·출처 확인 코드 |
| [`NewLora`](analyses/NewLora/README.md) | 새 LoRA ceiling, real/이전 생성본/floor 비교, 희귀 집단·DEG 및 재검증 |
| [`fig6_bulkonly_reproduction_20261020`](analyses/fig6_bulkonly_reproduction_20261020/README.md) | bulk-only Figure 6 재현 단계별 분석 |
| [`paper_claim_analysis`](analyses/paper_claim_analysis/README.md) | 원 논문의 세포 집단 주장과 PBISC 생성본 비교 |
| [`project_review_20261028`](analyses/project_review_20261028/README.md) | 정정된 DEG·유전자 지표와 결과표 감사 |
| [`cross_dataset`](analyses/cross_dataset/README.md) | 코호트 간 조성·발현 유사도와 H5AD 구조 확인 |
| [`manuscript_metrics_20261112`](analyses/manuscript_metrics_20261112/README.md) | gene-mean 곡선, iLISI, composition/within-type 분해, 희귀세포 통계의 재사용 모듈·테스트·데모 노트북 (문서 정의 기반 재구현) |
| [`fig7_rare_population_extension_20261116`](analyses/fig7_rare_population_extension_20261116/README.md) | Figure 7 ver5: ESKD COVID-19에서 bulk-only 환자 29명을 PBISC로 더했을 때 B ASC dividing·NK CD56bright·HSPC의 중증도·사망 연관과 검정력, 보충 그림 a–e |
| [`fig7_ver8_main_supp_20261117`](analyses/fig7_ver8_main_supp_20261117/README.md) | Figure 7 ver8과 보충 그림 1·2: CLL 제외 후 43개 세포형 5단계 선별(43 → 16 → 15 → 12 → 8), 신뢰도 스코어카드, HSPC 사망 강건성, 27명 모달리티 스케일 비교, 검산 노트북과 독립 재구현 |
| [`fig7_ver8_public_replication_20261118`](analyses/fig7_ver8_public_replication_20261118/README.md) | Figure 7 ver8 주장의 공개 COVID-19 코호트 재현 v8.1(Stephenson, COMBAT, van der Wijst, Liu; 219명, 코호트 구성 통일): 코호트별·메타분석 ρ와 사망 AUC, Stephenson B ASC isotype GLM, HSPC 사망 시간창, 보충 그림 3, 검산 노트북과 독립 재구현 |
| [`fig7_ver9_main_supp_20261118`](analyses/fig7_ver9_main_supp_20261118/README.md) | Figure 7 ver9 재구성(본 그림 a–q, 보충 그림 1–4): 코호트 구성·사망 확장 패널을 본 그림으로 옮기고 22개 세포형 스크린과 n80을 본문으로 이동, 공개 코호트 재현은 보충 그림 4. 그림 스크립트, 검산 노트북 |

## 실행 전 준비

1. 원본 분석을 실행하려면 원본·PBISC 생성 `.h5ad`와 관련 메타데이터를 별도로 준비합니다. 원본 데이터는 저장소에 포함하지 않았습니다. 현재 Figure 7·8 패키지의 그림 재현과 표 검산은 동봉된 자료만 사용하며, 해당 폴더 README의 실행법을 따릅니다.
2. [requirements.txt](requirements.txt)의 공통 패키지를 설치하고, 사용하는 분석의 추가 패키지를 해당 폴더 README에서 확인합니다.
3. 스크립트에 남아 있는 원래 환경의 절대 경로와 입력 파일 위치를 자신의 환경에 맞게 수정합니다. 일부 스크립트는 원래 작업 디렉터리의 중간 산출물 또는 다른 코드 폴더를 참조합니다. 폴더별 README와 코드의 입력 경로를 먼저 확인하세요.
4. 단계형 분석은 번호가 붙은 스크립트를 순서대로 실행합니다. 노트북은 원래 저장된 출력과 실행 번호를 제거해 코드와 설명 셀만 보관했습니다(`manuscript_metrics_20261112`의 합성 데이터 데모는 예외로 결과를 포함).

이 저장소는 당시 분석 코드의 **보관본**입니다. 모든 분석을 하나의 새 환경에서 즉시 재실행할 수 있도록 경로와 의존성을 통합한 파이프라인은 아닙니다. 결과 수치와 해석은 원래의 입력 파일, 필터, 모델·라벨 버전 및 실행 환경을 확인한 뒤 사용해야 합니다.

## 포함 범위

- Python, R, shell 스크립트와 Jupyter 노트북을 분석 시기와 주제별로 보존하며 후속 분석 패키지를 추가했습니다.
- 노트북 출력·첨부 이미지, H5AD/NPZ/Parquet, 발표자료, 복사된 Python 런타임은 제외했습니다. New LoRA CellDEEP 결과와 현재 Figure 7·8 통일 재분석 패키지에는 검토한 결과표·그림을 포함합니다.

