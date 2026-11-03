# PBISC scRNA-seq analysis code

PBISC가 생성한 단일세포 데이터의 검증·비교 분석 코드를 모은 저장소입니다. 로컬에 남아 있던 분석 스크립트와 노트북을 분석 시기와 주제별로 `analyses/` 아래에 정리했습니다. 주된 대상 코호트는 Ye 2026 sepsis와 Imperial ESKD/COVID입니다.

## 분석 폴더

| 폴더 | 내용 |
| --- | --- |
| [`bulkonly_lora5_20261006`](analyses/bulkonly_lora5_20261006/README.md) | LoRA-5 K=1,000 생성 데이터, 원본 scRNA-seq 재처리, 통합·조성·희귀 세포·질병 DEG·환자별 검증 및 후속 분석 |
| [`bulkonly_lora5_k4000_20261006`](analyses/bulkonly_lora5_k4000_20261006/README.md) | K=4,000 Harmony 분석과 K=1,000 비교 |
| [`bulkonly_lora5_20261007`](analyses/bulkonly_lora5_20261007/README.md) | sepsis K=1,000/K=4,000 클러스터 비교 |
| [`new_inference_depth256_K4000_floor1e-4_20261021`](analyses/new_inference_depth256_K4000_floor1e-4_20261021/README.md) | floor readout 데이터 다운로드·출처 확인 코드 |
| [`NewLora`](analyses/NewLora/README.md) | 새 LoRA ceiling, real/이전 생성본/floor 비교, 희귀 집단·DEG 및 재검증 |
| [`fig6_bulkonly_reproduction_20261020`](analyses/fig6_bulkonly_reproduction_20261020/README.md) | bulk-only Figure 6 재현 단계별 분석 |
| [`paper_claim_analysis`](analyses/paper_claim_analysis/README.md) | 원 논문의 세포 집단 주장과 PBISC 생성본 비교 |
| [`project_review_20261028`](analyses/project_review_20261028/README.md) | 정정된 DEG·유전자 지표와 결과표 감사 |
| [`cross_dataset`](analyses/cross_dataset/README.md) | 코호트 간 조성·발현 유사도와 H5AD 구조 확인 |

## 실행 전 준비

1. 분석하려는 원본·PBISC 생성 `.h5ad`와 관련 메타데이터를 별도로 준비합니다. 데이터 파일과 생성된 표·그림은 이 저장소에 포함하지 않았습니다.
2. [requirements.txt](requirements.txt)의 공통 패키지를 설치하고, 사용하는 분석의 추가 패키지를 해당 폴더 README에서 확인합니다.
3. 스크립트에 남아 있는 원래 환경의 절대 경로와 입력 파일 위치를 자신의 환경에 맞게 수정합니다. 일부 스크립트는 원래 작업 디렉터리의 중간 산출물 또는 다른 코드 폴더를 참조합니다. 폴더별 README와 코드의 입력 경로를 먼저 확인하세요.
4. 단계형 분석은 번호가 붙은 스크립트를 순서대로 실행합니다. 노트북은 원래 저장된 출력과 실행 번호를 제거해 코드와 설명 셀만 보관했습니다.

이 저장소는 당시 분석 코드의 **보관본**입니다. 모든 분석을 하나의 새 환경에서 즉시 재실행할 수 있도록 경로와 의존성을 통합한 파이프라인은 아닙니다. 결과 수치와 해석은 원래의 입력 파일, 필터, 모델·라벨 버전 및 실행 환경을 확인한 뒤 사용해야 합니다.

## 포함 범위

- Python, R, shell 스크립트와 Jupyter 노트북 356개를 원래 분석 폴더 구조에 가깝게 보존했습니다.
- 노트북 출력·첨부 이미지, H5AD/NPZ/Parquet, 생성 CSV, 그림, 발표자료, 복사된 Python 런타임은 제외했습니다.

