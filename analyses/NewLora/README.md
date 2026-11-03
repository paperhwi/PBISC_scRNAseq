# New LoRA ceiling 분석

원본 위치: `F:/PBISC/NewLora/`. 2026-11-02 새 LoRA ceiling 생성본을 관측 세포, 이전 K=4,000 생성본, floor readout과 비교한 코드입니다.

- 루트의 `analyze_ceiling_vs_real.py`, `compare_three_transcriptomes.py`, `analyze_disease_deg_matched.py`, `analyze_rare_disease_deg.py`, `rare_deg_feasibility.py`가 주요 수치 분석입니다.
- 루트의 `plot_*.py`, `recreate_*.py`, `summarize_*.py`는 표와 그림을 만듭니다.
- `analysis_20261102/file_overview/`는 새 파일 구조와 발현 개요를 추출합니다.
- `analysis_20261102/verification_20261102/`는 원본 무결성, label transfer, iLISI, DEG, 희귀 세포 출처 및 민감도 분석을 재검증합니다.
- `analysis_20261102/figure_realism_rare/`, `figure_additions/`, `celldeep_source_deg/`는 후속 그림과 비교 입력을 만듭니다.

입력 H5AD, 메타데이터, 생성 결과표는 포함하지 않았습니다. 이 코드는 원본 `F:/PBISC` 폴더와 앞선 분석의 `ilisi_lib.py` 등을 참조하므로 경로와 필요한 앞 단계 산출물을 확인하세요.
