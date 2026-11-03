# LoRA-5 K=1,000 및 후속 검증

원본 위치: `F:/PBISC/bulkonly_lora5_20261006/analysis/`. Ye 2026 sepsis와 Imperial ESKD/COVID의 bulk-only PBISC 생성 데이터와 관측 scRNA-seq를 비교한 코드입니다. 이 폴더의 하위 구조는 원본을 유지했습니다.

| 하위 폴더 | 분석 |
| --- | --- |
| `analysis/original_scRNA/` | 원본 scRNA-seq 다운로드, QC, 재주석, H5AD 구성 |
| `analysis/harmony/` | sepsis 생성 세포 Harmony 통합 |
| `analysis/integration_report_20261012/code/` | 통합 지표, 조성, 깊이, 희귀 집단, 그림 생성의 초기 분석 |
| `analysis/integration_report_20261012/*/` | 세포형·myeloid·T 세포 라벨, iLISI, HVG/Harmony, floor readout, 저발현 유전자, 세 모델 비교, 질병 DEG·희귀 집단, 비율 및 배치 효과 등 주제별 후속 분석 |
| `analysis/paired_patient_validation_20261014/code/` | 동일 환자 기준 생성·관측 비교 |
| `analysis/cell_transcriptome_validation_20261014/code/` | 세포 수준 전사체 유사도와 대조군 검증 |
| `analysis/new_inference_depth256_K4000_floor_20261021/` | 새 floor 생성본을 가져오는 보조 코드 |

`analysis/integration_report_20261012/` 안의 각 주제 폴더에는 해당 분석의 코드와 노트북이 있습니다. `code/`의 로컬 모듈을 서로 import하므로 같은 폴더 안의 파일을 함께 사용하세요. 원본 스크립트는 `F:/PBISC`와 과거 sandbox의 `/home/claude`, `/mnt/user-data` 경로를 참조할 수 있습니다. 데이터 위치를 수정한 뒤 개별 주제를 실행해야 합니다. 노트북 출력은 제거했습니다.

주요 외부 의존성: `scanpy`, `anndata`, `harmonypy`, `scikit-learn`, `statsmodels`. 선택 분석에는 `scvi-tools`, `bbknn`, `scanorama`, `gseapy`, `python-igraph`, `leidenalg`, `liana`, `BayesPrism`(R)이 필요합니다.
