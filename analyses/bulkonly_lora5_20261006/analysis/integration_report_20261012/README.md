# 통합 분석 및 후속 검증 목차

이 폴더는 K=1,000 LoRA-5 생성 세포를 관측 세포와 비교한 초기 통합 분석과 이후 확장 분석을 담습니다. 각 주제의 스크립트와 노트북은 원본 폴더 구성을 보존했습니다. 노트북의 저장 출력은 제거했습니다.

| 폴더 | 분석 주제 |
| --- | --- |
| `code/` | 기본 통합 지표, 조성, 발현 깊이, 희귀 집단, 그림 |
| `celltype_umap_20261016/` | sepsis 세포형 UMAP |
| `myeloid_subcluster_20261016/` | 골수계 세포 세분화 |
| `macrophage_labeling_20261017/`, `macrophage_markers/` | 대식세포 라벨과 마커 |
| `cd4_cycling_20261017/`, `nk_ilc_markers/` | CD4 cycling 및 NK/ILC 마커 |
| `rf_realism_20261017/` | Random Forest 기반 진위 판별 |
| `file_component_fig_20261018/` | 파일 성분·floor·관측 데이터 비교 |
| `patient_pair_ilisi_20261018/`, `ilisi_explainer_20261018/` | 환자 쌍 및 iLISI 해석 |
| `hvg1000_harmony_ilisi_20261020/`, `hvg_harmony_diag_20261022/` | HVG 선택과 Harmony/iLISI 진단 |
| `lowexpr_cliff_20261020/`, `lowexpr_floor_20261021/` | 저발현 유전자 절벽과 floor 영향 |
| `ilisi_floor_20261021/` | floor 생성본의 통합·발현·희귀 집단 검증 |
| `three_model_20261022/`, `scsemi_eval_20261022/` | PBISC·floor·scSemiProfiler 등 모델 비교 |
| `proportion_qualification/` | 세포 조성, deconvolution, pseudobulk 대조 |
| `week1_defensive_20261027/`, `week2_attribution_20261027/` | 질병 신호의 방어적 검증과 기여 분해 |
| `disease_deg_rare_20261027/` | 질병 DEG와 희귀 집단, matched-N, 마커 주석 |
| `tcell_batch_floor_20261027/` | T 세포의 배치·floor 효과 |
| `dc_underproduction_20261030/` | DC 세포의 과소 생성 검증 |
| `monocyte_gene_level/`, `labeling_deck/`, `claim_classifiers/` | 보조 유전자·라벨·주장 분석 |

일부 폴더의 `build_nb*.py` 및 `build_deck*.py`는 분석 결과를 노트북·발표자료로 조립한 코드입니다. 실제 수치 분석은 같은 폴더의 지표 계산 스크립트와 노트북을 확인하세요. 과거 sandbox 절대 경로를 사용하는 파일은 데이터 위치를 먼저 수정해야 합니다.
