# Figure 4–7 및 보충그림 8–27 검산 보고서

검산일: 2026-11-26. 원본 폴더 `nature_submission_fig4_7_final_20261125`는 보존하고, 수정본을 `nature_submission_fig4_7_reviewed_20261126`에 작성했습니다. 이 보고서와 `VERIFICATION_REVIEWED_20261126.json`이 11월 25일의 전체 검증 완료 표기보다 우선합니다.

## 판단

이전 방법론 수정은 대부분 원고에 반영되어 있었습니다. 그러나 외부 사망 라벨 정정 명세서는 실제 데이터·그림·원고에 반영되지 않았습니다. 단순 재계산 일치와 임상 원자료의 정확성은 다릅니다. 원본 패키지의 87개 자체 점검은 모두 통과했지만, 잘못된 사망 라벨을 그대로 사용하고 있었습니다.

가장 중요한 수정은 외부 van der Wijst 코호트입니다. 저자 임상 원본에서 42명을 환자 ID로 대조해 사망 7명을 확인했고, 모든 방문에 같은 최종 사망 라벨을 적용했습니다. 외부 CD16 단핵구의 마지막 채혈 감소는 정정 후 다중보정을 통과하지 않습니다. 따라서 기존의 “양측 코호트에서 말기 CD16 단핵구 감소로 수렴한다”는 주장을 삭제했습니다.

## 1. 사망 라벨과 외부 메타분석

[저자 원본 Data files S1–S9](https://pmc-oa-opendata.s3.amazonaws.com/PMC8601717.2/scitranslmed.abh2624_data_files_s1_to_s9.zip)의 `Data_file_S3.xlsx`, `ClinicalCharacteristics`, `Donor_ID`/`Died`를 사용했습니다. 사망자는 `Rep_C_1029, Rep_C_1037, Rep_C_1047, Rep_C_1050, Rep_C_1130, Rep_C_1154, Rep_C_1161`입니다.

기존 표는 같은 환자의 방문 사이에서도 최종 사망 라벨이 달랐습니다(8명). 기존 첫 채혈 라벨과 원본은 10명, 마지막 채혈 라벨과 원본은 12명이 불일치했습니다. 따라서 전달된 두 설명서의 “10명”과 “12명”은 선택한 채혈 시점에 따라 각각 성립합니다. 환자를 임의로 제외하지 않고 원본의 최종 결과를 모든 방문에 고정했습니다.

[Stephenson 2021 저자 Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-021-01329-2/MediaObjects/41591_2021_1329_MOESM3_ESM.xlsx)도 확보해 89개 검체 ID·환자 ID·사망 결과가 기존 표와 맞는지 대조했습니다. AP3의 `Deceased`를 포함한 73명/7사망을 유지했습니다. COMBAT와 Liu의 라벨은 이번 수정에서 바꾸지 않았으며, 두 코호트 임상 원자료를 새로 수집해 전수 판정한 것은 아닙니다.

| 항목 | 수정 전 | 수정 후 | 의미 |
| --- | --- | --- | --- |
| van der Wijst | 42명, 5사망; 방문 간 불일치 | 42명, 7사망; 환자별 고정 | 저자 원본으로 정정 |
| 외부 결과 기록 집단 | 219명, 197생존/22사망 | 219명, 195생존/24사망 | 검체·환자·방문 선택은 유지 |
| CD8 CTL proxy 첫 채혈 | AUC 0.477, q23 0.926 | 0.411, 0.343 | 여전히 외부 재현되지 않음 |
| NK 첫 채혈 | 0.496, 0.965 | 0.476, 0.834 | 여전히 외부 재현되지 않음 |
| HSPC 첫 채혈 | 0.569, 0.623 | 0.660, 0.088 | BH23 유의하지 않음 |
| MAIT 첫 채혈 | 0.221, 0.00642 | 0.211, 0.00241 | 감소 유지 |
| CD8 naïve T 첫 채혈 | 0.266, 0.0000345 | 0.323, 0.00760 | 감소 유지 |
| CD16 mono 마지막 채혈 | 0.320, 0.0195 | 0.394, 0.217 | 유의성 상실; 해석 수정 |
| 마지막 채혈 q23 < 0.05 | 기존 그림과 상이 | MAIT, CD8 naïve T, DC2, γδ T | 표·그림·본문 동기화 |

정정 전 89개 코호트×세포군의 첫/마지막 AUC와 P(총 178개 결과)를 먼저 재현했습니다. 정정 후 AUC는 별도의 모든 사망–생존 쌍 비교로 다시 확인했고, DL random effects와 정렬 기반 BH23도 별도로 대조했습니다. 공유 발현·세포 분율·중증도·시간 순서는 바꾸지 않았습니다. 코호트별 정정된 vdW CI는 같은 결과층별 환자 bootstrap 방식으로 10,000회(seed 20261126) 다시 산출했습니다. 메타분석 CI는 기존 Wald 방식입니다.

반영 파일: 외부 검체·코호트별 결과·메타분석·panel l CSV 4개, 매니페스트, Excel 4개 대응 시트와 인덱스, Figure 7k 예시와 7l forest/table, 원고 Results/Discussion/Methods/범례. 내부 ESKD의 66명/11사망·GEE·복합점수·내부 CV는 외부 라벨을 사용하지 않아 유지됩니다.

## 2. 추가로 고친 부분

| 위치 | 확인한 문제 | 수정과 근거 |
| --- | --- | --- |
| Fig. 4 본문/Methods | `within-donor case − control`은 같은 사람이 case/control이라는 뜻으로 읽힘 | Real/생성 source가 같은 held-out donor를 공유하는 paired 평가임을 명시. 질병군 대조 자체를 환자 내 변화로 기술하지 않음 |
| Fig. 5a | 그림에 `other donors, unpaired` 잔존 | `unmatched specimens`로 국소 벡터 수정. COVID Real B/PBISC에 공통 환자 7명 존재 |
| Fig. 5 DEG/분할 신뢰도 | PBISC 전체 488k/656k를 사용한 것으로 읽힐 수 있음 | 실제 60,000-cell 평가 부분표본, seed 20261024 및 frozen six_tr를 Methods에 명시 |
| Fig. 5 Real BH 차이 | 최대 Δq 약 0.0025 | 8개 표 재계산의 최대 0.00344685로 정정. DEG 이진 판정 변화는 0개 |
| scSemiProfiler 설명 | `representative-donor copying` | 실제 생성 모형을 donor-copying으로 단정하지 않도록 동일 support의 comparator로 기술 |
| Fig. 6 해석 | `genuine` 상태 인증 및 예측 일반화처럼 읽히는 문구 | 전이 annotation과 타 환자 분자 유사도라는 실제 검증 범위로 정밀화 |
| Fig. 7a | 37명 모두 `Bulk + scRNA-seq`인 도식 | 37 scRNA-seq 환자 = 10 only + 27 shared로 표시 |
| Fig. 7 내보내기 | 원본 768×1024가 기본 렌더에서 576×768로 축소 | 원본 픽셀을 PDF에 직접 삽입하고 국소 표기를 수정. 기존 forest plot 구성을 보존하고 정정된 외부 결과만 갱신 |
| Supp. 20 | CI와 pre/post shift 선 구분이 부족 | source 색 CI 가로선과 회색 이동 연결선의 key/범례 추가 |
| Supp. 21 | `power`라는 명칭 | 관측 코호트 내 비복원 resampling 검출률로 수정. prospective power로 해석하지 않음 |
| Supp. 22 | 내부 ESKD key가 13사망 | 11사망으로 정정. 별도의 peak-grade 외부 219명/23사망과 주분석 219명/24사망의 집합 차이, 내부 ESKD 논문 2025와 외부 Stephenson 2021을 구분 |
| Supp. 23c | 초기 HSPC 표에 누락값·중복 cycling·구버전 HSC5 값이 혼재 | 원시 HSPC 행렬에서 unique16 프로그램을 재계산해 16행을 갱신하고 모든 q16를 표시. 나머지 24 Hallmark 값은 기존 분석으로 유지 |
| Supp. 26 | 범례의 마지막 채혈 표적이 `CD16 mono` | 실제 패널인 `CD16 mono C1`로 정정. 276-draw sensitivity와 주 217-draw 분석을 분리 |
| Supp. 27 | 다른 외부 mapping/대상 수/BH family가 주분석과 섞여 설명됨 | 상단 severity는 historical direct mapping/BH18로 명시. 하단 mortality는 이번 정정된 주분석 mapping/BH23 결과로 갱신: HSPC q23 .209, γδ T .039, C1은 1개 코호트라 pooled q 없음 |
| 변경 이력 | 정수 count 분자/분모 오기 | 9/48,255→9/47,650; 521/272,000→522/272,000; 20/41,720→20/41,457. 실제 원시 표는 이미 맞았음 |

scSemiProfiler는 [원 논문](https://www.nature.com/articles/s41467-024-50150-1)의 deep generative inference 방법입니다. 위 문구 수정은 비교 성능값을 바꾸지 않습니다.

Supp. 23은 33 measured + 29 bulk-only, 총 62명/11사망, 공통 5,762 유전자, 환자당 HSPC ≥3개 기준입니다. 모든 방문의 cell-level log1p(CP10K) 환자 평균 → source별 normal score → Mann–Whitney/BH16을 재계산했습니다. 48개 집단별 AUC의 최대 차이는 0, P/q 차이는 약 3.1×10⁻⁸ 이내였습니다. 예를 들어 HSC5 Combined AUC/P는 초기 표의 약 .41/.333 대신 최종 unique16 분석의 .512/.905입니다. 16개 HSPC CI는 보관된 10,000회 bootstrap 값이며, 동일 RNG 경로까지 새로 재현한 것은 아닙니다. Ery5 AUC=.752, nominal P=.0094, q16=.150의 해석은 유지됩니다.

## 3. 이전 감사 결과 중 이번에 정정한 판단

제가 앞서 Fig. 5f의 scSemiProfiler r을 .936으로 일괄 수정한 것은 입력 집합을 충분히 구분하지 못한 판단이었습니다. 같은 5,961개 유전자에서 보관된 Real 평균과 전체 51,000개 scSemiProfiler 세포의 평균을 비교하면 **r=.936545748 → .937**입니다. 양쪽을 보관 display 표에서 읽으면 **.936305082 → .936**입니다. 따라서 최신 그림의 .937을 유지하고, Methods에 두 계산 입력을 정확하게 설명했습니다. Real 전체 파일의 평균까지 새로 계산한 값이라고 표현하지 않았습니다.

## 4. 실제 검산 범위

- 원 패키지 자체 점검: 87/87 통과. 이는 원본 라벨의 임상적 진실성을 보장하지 않았습니다.
- 처음 수행한 추가 점검: 75개 통과(24개 PDF/PNG 비율, Excel/CSV, 전일 검산 결과의 동일성, 새로 제공된 pseudobulk split-half 재현 포함).
- 수정 후 독립 점검: 56개, **Excel 데이터 92,537셀**을 대응 CSV와 전수 대조. 22개 CSV 해시 확인. 외부 AUC 178개를 별도 pairwise 방식으로 재확인, DL/BH23 독립 대조.
- Fig. 5 split-half: 새로 제공된 실제 pseudobulk 6개에서 16개 cohort×cell-type×method 행의 reliability, disattenuated r, corrected slope, AUROC/recall을 재현(최대 수치 차이 약 2.2×10⁻¹⁶).
- Fig. 5 Real genome-wide vs panel BH: 8개 Real 표 재검산, DEG 판정 차이 0개.
- Fig. 5 조성: 297개 donor×condition, 1,782개 유형별 값은 전일 원시 집계와 동일. 공유 환자 multiplicity를 보존한 별도 49,999회 bootstrap 민감도 분석에서 T/B의 유의/비유의 분류는 같았습니다. 이것이 최종 그림의 모든 정확한 q·별 개수 재현을 뜻하지는 않습니다.
- Fig. 6: 원시 정수 분자/분모 및 895개 타 환자 유사도 결과가 전일 독립 계산과 동일. Figure 4·6 원본 PDF는 byte 단위로 보존했습니다.
- Fig. 7 내부: 기존 217-draw/66-first-draw/63-expression 데이터와 통계 표가 전일 검산본과 동일. 원자료 행렬 기반 내부 GEE·rank INT·사망 AUC·CV에 관한 전일 결과를 그대로 유지했습니다.
- 본문 24개 이미지, 개별 보충그림 20종 및 40-page 합본은 최신 파일로 동기화했습니다. Word/PDF 렌더와 파일별 구조 점검 결과는 JSON에 기록합니다.

## 5. 남는 한계 — 제출 전 추가 원자료가 필요한 부분

1. **Figure 4 원래 limma-voom fit과 per-gene 결과**: 로컬 handoff 설명서도 원래 fit이 외부 실행 환경에 있고 원본 per-gene 테이블/실행 스크립트가 이관되지 않았다고 명시합니다. 현재 16행 요약표(파일명은 32pairs)는 PDF 막대·점 좌표를 읽어 만든 근사치입니다. 점에서 계산한 평균과 저장된 막대 평균의 최대 차이는 .001033이며, 소수 넷째 자리까지 정확한 원시 통계로 취급할 수 없습니다. 이를 점 평균으로 덮어써 정확한 통계인 것처럼 만들지 않았습니다. iLISI dot 값도 원자료가 아니라 그래픽 추출값입니다.
2. **Figure 5n/o 14개 programme**: 42행 표가 추가됐지만 PDF에 표시된 효과·유의성의 전사입니다. 실제 최종 gene-set 목록과 환자별 점수·P·CI가 아니므로 6/7·7/7 재현을 독립적인 원자료 검산 완료로 인정할 수 없습니다.
3. **Figure 5m 및 Supp. 15의 최종 bootstrap**: 마지막 제작 코드/RNG 경로가 확인되지 않았습니다. 이번 별도 patient-block 민감도 분석은 T/B 유의/비유의 분류를 지지하지만, 예를 들어 COVID B case scSemiProfiler는 q24≈.0012(**)로 기존 ***와 차이가 있습니다. 원래 최종 q를 확정할 입력/코드 없이 기존 별을 임의의 새 분석으로 교체하지 않았습니다. 비유의는 동등성 증명이 아닙니다.
4. **Supp. 27 상단과 276-draw 내부 sensitivity**: 주 217/300 분석과 다른 역사적 표본 선정·mapping을 유지합니다. 상단 severity의 원 study-level 입력을 확보하지 못해 완전 재현하지 않았습니다. 하단 외부 mortality는 새 주분석으로 교체했고 BH family 변경을 명시했습니다. 두 블록을 하나의 동일 분석으로 인용하면 안 됩니다.
5. **Figure 7 원본 품질**: 기존 최종 디자인의 실제 원본이 768×1024 래스터입니다. 183 mm 폭에서 약 107 ppi이며 전체 벡터/고해상도 원본으로 복원된 것은 아닙니다. 정정한 표·forest plot·일부 글자는 벡터이지만 나머지 래스터 영역의 세부정보는 원본 한계가 남습니다.
6. 다른 모든 보충 패널의 원래 raw-expression·bootstrap·모형 적합을 이번에 새로 전수 실행한 것은 아닙니다. 점검 통과를 “모든 패널의 완전 독립 재현”으로 확대하지 않았습니다.

## 증거와 재현 파일

`verification_20261126/`에 임상 원본 출처·해시, 환자별 사망 라벨 대조, 수정 후 원자료 재계산, 검사 결과 및 코드가 있습니다. 이전 변경 이력과 기준 exporter는 `history_20261125/`에 보존했습니다. 수정 패키지의 baseline exporter는 임상 정정을 덮어쓰지 못하도록 중지시켰습니다. 최신 manuscript builder와 package checker는 `code/02...`, `code/03...`입니다.

최종 조판 확인: Word 원고 37페이지와 Microsoft Word 내용 전용 PDF 37페이지를 렌더링하여 전 페이지를 시각 점검했습니다. 대체 ReportLab PDF도 37페이지를 점검했습니다. 보충 합본은 40페이지이며 20개 전체 캡션을 별도로 렌더링했습니다. 검사상 대체 문자·검은 사각형·페이지 밖 텍스트 블록은 0개였습니다. 배포 PDF는 Word의 조판과 동일한 판본입니다.
