# Real vs New LoRA: CellDEEP source DEG

이 폴더는 **같은 코호트·질병군·큰 세포형 안에서 New LoRA와 Real의 발현 차이**를 분석한 결과입니다. 질병군 간 DEG 분석이 아닙니다. 양쪽에 동일 환자 ID가 없으므로 짝지은 검정이나 환자별 재현성 검정으로 해석할 수 없습니다.

## 입력과 비교 설계

- Sepsis: Ye 2026 원자료의 성인·논문 분류 대상 Real 세포와 `NewLora/ye2026_sepsis_draw1_K4000_ceiling.h5ad`.
- COVID: Imperial ESKD 원자료의 발표된 187-sample roster와 `NewLora/imperial_eskd_draw0_K4000_ceiling.h5ad`. 기존 분석에서 제외한 5개 표본은 양쪽에서 제외했습니다.
- New LoRA는 기존 분석의 각 코호트 60,000개 eligible 세포와 이전에 계산한 Real A 기반 30-NN 전이 `six_tr` 주석을 사용했습니다. Real은 원 논문 세포형을 6개 큰 세포형으로 매핑했습니다. 주석 방식이 서로 다르므로 source 차이에는 주석 오차가 포함될 수 있습니다.
- T, B, NK, Monocyte 네 세포형을 질병군별로 따로 분석했습니다. 각 질병군·세포형에서 양쪽 환자 수를 맞추고, 환자당 세포 30개를 무작위 추출했습니다. 환자 수는 최대 20명/쪽이며 5명/쪽 미만은 제외했습니다. COVID에서 반복 방문 환자는 질병군별 eligible 세포가 가장 많은 방문 하나를 선택했습니다.
- 총 19개 비교: Sepsis 12개, COVID 7개. COVID 음성 B cell은 New LoRA에서 기준을 만족하는 환자가 2명뿐이어서 제외했습니다. `sepsis/strata.csv`, `covid/strata.csv`에 대상과 제외 사유가 있습니다.
- Real과 New LoRA의 공통 유전자: Sepsis 5,961개, COVID 5,745개. 원시 정수 UMI count를 사용했습니다.

## 분석법

`run_celldeep.R`에서 R **CellDEEP 1.0.1**의 `FindMarker.CellDEEP`을 사용했습니다. 환자·source·세포형 안에서 무작위 세포 10개씩 합산하여 환자당 metacell 3개를 만들고 Wilcoxon 검정을 수행했습니다. 모든 반환 유전자의 p값에 Benjamini–Hochberg 보정을 따로 적용했습니다. CellDEEP의 Seurat `p_val_adj`는 Bonferroni이므로 최종 판정에는 `q_BH`를 사용합니다.

같은 입력 세포를 환자별로 합산하여 **edgeR `filterByExpr` + TMM + limma-voom**으로 독립 환자 단위 검정을 추가했습니다. 두 방법 모두 `BH q < 0.05` 및 `|log2FC| ≥ 1`을 DEG 기준으로 했습니다. log2FC의 양수는 **New LoRA > Real**입니다. `confirmed_degs.csv`는 두 방법에서 모두 기준을 통과하고 방향이 일치한 호출입니다. 동일 유전자가 다른 질병군·세포형에서 여러 번 등장할 수 있습니다.

## 결과 요약

| 비교 | CellDEEP 호출 | 환자 pseudobulk 호출 | 양쪽에서 같은 방향 |
|---|---:|---:|---:|
| Sepsis · sepsis · T | 359 | 43 | 38 |
| Sepsis · sepsis · B | 439 | 66 | 61 |
| Sepsis · sepsis · NK | 421 | 54 | 48 |
| Sepsis · sepsis · Monocyte | 736 | 107 | 86 |
| COVID 양성 · T | 302 | 29 | 27 |
| COVID 양성 · B | 596 | 84 | 76 |
| COVID 양성 · NK | 412 | 49 | 44 |
| COVID 양성 · Monocyte | 431 | 63 | 55 |

19개 비교 전체에서 CellDEEP 5,964건, 환자 pseudobulk 1,115건, 양쪽 일치 968건입니다. 이는 **고유 유전자 수가 아닌 유전자×질병군×세포형 호출 수**입니다. pseudobulk 호출의 86.8%가 CellDEEP에서도 같은 방향으로 검출되었지만, CellDEEP 호출 중 환자 단위로 확인된 비율은 16.2%였습니다. CellDEEP의 metacell은 동일 환자에서 여러 개 생성되므로 독립 환자 수를 늘려주지 않습니다. 따라서 CellDEEP 단독 유의성은 환자 수준의 확정적 증거로 쓰지 않는 것이 타당합니다.

예시로 sepsis Monocyte에서 양쪽 방법이 확인한 유전자는 New LoRA에서 높은 `GNG2`, `MICU2`와 낮은 `PLBD1`, `S100A8`, `S100A9`, `S100A12`입니다. COVID 양성 B cell에서는 `TANK`, `SSH2`, `UVRAG`, `BCL2`가 New LoRA에서 높았습니다. 이 차이는 환자 집단, 포착·생성 과정, 주석 방식의 차이를 분리하지 못하므로 생성 모델의 인과적 효과로 단정할 수 없습니다.

## 파일

이 저장소에는 아래 요약표·유전자별 결과·그림을 포함합니다. 원본 H5AD, 환자별 선택 목록, 추출된 count 행렬과 세포별 metadata는 포함하지 않습니다.

- `results/summary.csv`: 19개 비교의 환자·세포·metacell 수와 DEG 수.
- `results/*_deg.csv.gz`: 유전자별 전체 CellDEEP 및 pseudobulk 통계. `q_BH`, `avg_log2FC`, `pb_q_BH`, `pb_logFC`, 두 DEG 판정 플래그를 포함합니다.
- `results/confirmed_degs.csv`: 양쪽 방법에서 유의하고 방향이 일치한 968건.
- `results/recurrent_confirmed_genes.csv`: 질병군·세포형에 걸친 반복 확인 유전자.
- `results/celldeep_vs_pseudobulk_summary.png` 및 `.pdf`: 비교별 DEG 수.
- `results/volcano_sepsis_adult_sepsis.png` 및 `.pdf`: Sepsis 환자 네 큰 세포형의 source-associated volcano plot.
- `results/volcano_covid_covid_positive.png` 및 `.pdf`: COVID 양성 환자의 source-associated volcano plot.

## 재현

원본 H5AD와 기존 중간 산출물을 이 분석 당시의 `F:/PBISC` 구조로 준비한 뒤, 저장소 루트에서 다음 순서로 실행합니다. 다른 위치라면 `PBISC_DATA_ROOT`에 해당 폴더를 지정합니다. `prepare_inputs.py`는 큰 원자료를 다시 읽어 입력을 생성하므로 수 분 이상 걸릴 수 있습니다. 결과와 임시 입력은 기본적으로 이 README가 있는 폴더에 기록되며, `PBISC_CELLDEEP_OUTPUT_ROOT`로 출력 위치를 바꿀 수 있습니다.

```powershell
$env:PBISC_DATA_ROOT = 'F:\PBISC'
$env:PYTHONPATH = 'F:\PBISC\.analysis_runtime'
py -3.12 analyses/NewLora/analysis_20261102/celldeep_source_deg/prepare_inputs.py
& 'C:\Program Files\R\R-4.3.2\bin\Rscript.exe' analyses/NewLora/analysis_20261102/celldeep_source_deg/run_celldeep.R
py -3.12 analyses/NewLora/analysis_20261102/celldeep_source_deg/summarize_plot.py
```

R에는 `CellDEEP`, `Seurat`, `Matrix`, `edgeR`, `limma`, `data.table`이 필요합니다. Python에는 `anndata`, `h5py`, `numpy`, `pandas`, `scipy`, `matplotlib`, `adjustText`가 필요합니다. 무작위 추출과 metacell 분할은 코드에 고정된 시드를 사용합니다.

## 해석상 제한

이 분석은 source-associated 차이를 찾는 기술 비교입니다. Real과 New LoRA 환자가 겹치지 않고, Real은 저자 주석, New LoRA는 전이 주석을 사용합니다. 질병군 안에서도 나이·성별·중증도·처치 등의 환자 특성이 맞춰진 것은 아닙니다. 환자당 30개 세포와 최대 20명의 단일 무작위 추출이므로 효과 크기와 DEG 수는 선택에 따라 달라질 수 있습니다. 희귀 세포형은 이 분석의 대상이 아닙니다.

CellDEEP 문서: <https://cran.r-universe.dev/CellDEEP/doc/manual.html>
CellDEEP 개발 저장소: <https://github.com/sii-scRNA-Seq/CellDEEP>
