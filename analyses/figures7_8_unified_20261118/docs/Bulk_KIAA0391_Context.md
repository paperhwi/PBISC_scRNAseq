# HSPC KIAA0391 후보의 원래 bulk RNA 문맥 확인

2026-11-18. 이미 선택된 **한 후보**에 대한 해석 범위 검사다. 후보를 추가하거나 새로운 FDR 발견 family를 만들지 않았다. Primary RNA 스크린의 genuine New 29명과 동일한 bulk 79 draw를 고정했다. 모든 draw에 PBISC HSPC가 존재했고 총 2,255개 세포였다. 모든 상관은 환자 29명을 독립 단위로 계산했다.

## ID와 정규화

PBISC rare H5AD의 gene_symbol `KIAA0391`에 매핑되는 var는 정확히 1개, 0-based 위치 4,177의 **ENSG00000100890**이다. Raw bulk CSV의 대응 행은 **ENSG00000100890.16** 1개다. Ensembl version suffix를 제거한 뒤 같은 Ensembl ID 행을 합산했다. 이 raw CSV에는 `__`로 시작하는 HTSeq summary 행이 없었다. Target ID를 다른 symbol 또는 다른 Ensembl로 바꾸지 않았다. `Bulk_KIAA0391_geneID_mapping.csv`와 manifest에 원래 ID를 보존했다.

Bulk draw마다 `log2(1 + 10^6 × KIAA0391 count / library sum)`을 구했다. Primary bulk 요약은 79 draw를 각 환자 안에서 동일 가중치로 평균한 값이다. 분모는 (1) raw bulk의 모든 60,605 unique Ensembl 행 또는 (2) primary shared 5,762 symbol의 first-occurrence Ensembl 중 bulk에서 관측되는 **5,756개 ID**의 합이다. Shared 6개 Ensembl ID가 raw CSV에 없어 (2)는 available intersection이다. 이 부재를 생물학적 측정 0으로 해석하지 않았다. 추가 duplicate-symbol-sum 분모 민감도도 계산했다. 그 경우 요청된 raw Ensembl 5,772개 중 5,762개가 bulk에 존재한다. 모든 ID별 포함 여부와 누락 목록은 `Bulk_shared_gene_mapping_availability.csv`와 manifest에 있다.

같은 79 draw의 HSPC KIAA0391 점수는 원래 rare H5AD에서 별도로 다시 계산했다. Shared first-occurrence 5,762개 유전자에 대한 세포별 CP10K→log1p 후, 환자의 모든 HSPC 세포를 동일 가중치로 평균한 것이 primary HSPC 점수다. Equal-visit HSPC 점수도 함께 보존했다. Bulk와 HSPC의 단위·정규화·가중치는 같다고 가정하지 않았다. 이를 확인하기 위해 bulk draw를 그 draw의 HSPC cell count로 가중한 민감도도 추가했다.

## 환자 단위 결과

| Bulk 분모 및 환자 내 draw 가중치 | Bulk–peak WHO rho | P | Bulk–HSPC primary rho | P |
|---|---:|---:|---:|---:|
| 모든 60,605 Ensembl, 동일 draw 가중치 | 0.16696 | 0.38669 | 0.12210 | 0.52807 |
| Shared first-occurrence 중 관측 5,756 ID, 동일 draw 가중치 | 0.13625 | 0.48096 | 0.10915 | 0.57300 |
| 모든 60,605 Ensembl, HSPC cell count 가중치 | 0.17770 | 0.35640 | 0.17114 | 0.37474 |
| Shared first-occurrence 중 관측 5,756 ID, HSPC cell count 가중치 | 0.11834 | 0.54093 | 0.11563 | 0.55032 |

Duplicate-symbol-sum 분모 민감도는 이 후보의 환자 순위·상관·P가 shared first-occurrence 결과와 같았다. Bulk와 equal-visit HSPC의 상관도 모두 CSV에 보존했으며 명목상 유의하지 않았다. Primary HSPC 점수와 peak WHO의 상관은 rho=0.80286, P=1.60×10^-7이고 equal-visit HSPC 점수에서는 rho=0.79045, P=3.37×10^-7이다. 이 두 값은 기존 primary RNA 스크린과 일치한다.

## 이 검사가 허용하는 해석

이 환자 집단에서는 **동일 KIAA0391 유전자의 bulk 발현 점수**만으로 강한 HSPC 점수–중증도 연관을 단순 설명한다는 증거가 없다. 그러나 bulk와 HSPC 상관의 비유의성을 세포 특이성 또는 bulk 신호와의 차이 유의성으로 바꾸어 해석할 수 없다. N=29의 비유의 결과는 효과가 정확히 0이라는 동등성 증거도 아니다.

PBISC는 다유전자 bulk와 참조 정보를 이용해 생성되며 HSPC 상태는 직접 측정되지 않았다. 이 단일 유전자 bulk 검사가 PBISC HSPC 점수의 생물학적 진실성·발현의 세포 귀속·임상적 효용·외부 재현을 검증하지 않는다. 특히 기존 shared 환자에서 이 후보의 all-visit Real–PBISC HSPC 상관이 음수라는 fidelity 한계는 그대로 남는다. 반대로 어떤 정규화에서 bulk P가 유의하게 나왔더라도 그것만으로 inferred HSPC 연관의 수학적 계산 오류가 되는 것은 아니다. Bulk와 inferred HSPC는 서로 다른 estimand다.

이 검사는 이미 선택된 후보의 후향적 문맥 검사이므로 위 P를 새 발견·다중검정 보정 완료·독립 validation으로 제시하지 않는다. RNA 전수 스크린의 BH5,762/BH51,858/BH103,716 및 source-stratified permutation과 분리해 해석해야 한다.

## 재현 자료

- `run_bulk_kiaa0391_context.py`: core 코호트 사용, 원자료 읽기 전용.
- `Bulk_KIAA0391_context.csv`: 분모 3개×draw 가중치 2개의 모든 환자 단위 상관.
- `Bulk_KIAA0391_same79draw_values.csv`: sample ID, 환자 ID, 임상 시점, raw count, library sum, 정규화 bulk 점수, HSPC 점수/세포 수.
- `Bulk_KIAA0391_same29patient_values.csv`: 같은 순서의 29명에 대한 모든 환자 벡터.
- `Bulk_KIAA0391_HSPC_reference_associations.csv`: primary/equal-visit HSPC–중증도 값.
- `Bulk_KIAA0391_manifest.json`: raw bulk/H5AD/고정 환자 행렬/core/코드의 SHA256·크기·수정 시각과 정확한 mapping/누락 ID.

원래 bulk CSV, H5AD 및 15–24 폴더를 수정하지 않았다.
