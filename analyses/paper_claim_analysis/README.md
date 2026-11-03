# 원 논문 주장 검증

원본 위치: `F:/PBISC/paper_claim_analysis/`. Ye sepsis와 Imperial ESKD/COVID 원 논문의 세포 집단·마커 주장을 관측 데이터와 PBISC 생성본에서 확인하는 코드입니다.

`import_ye_table6.py`, `compute_ye_marker_coverage.py`로 원 논문의 마커 정의를 준비하고, `fixed_original_claim_validation.py`로 관측 코호트에서 주장과 표본 구성을 확인합니다. `run_fixed_claims_k_comparison.py`와 `summarize_fixed_claim_k_comparison.py`는 고정된 정의로 K=1,000/K=4,000 생성본을 비교합니다. COVID 표본 ID 확인에는 `reconcile_covid_ids.py`와 `compare_covid_manifest.py`를 사용합니다.

논문 보조표에서 유래한 마커 CSV, 환자별 결과표, 원본·생성 H5AD는 포함하지 않았습니다. `fixed_paper_claim_definitions.csv` 등 선행 산출물을 만들고 스크립트의 입력 경로를 맞춘 뒤 실행하세요.
