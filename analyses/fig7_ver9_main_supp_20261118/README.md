# Figure 7 ver9 — 본 그림 a–q와 보충 그림 1–4 재구성 (2026-11-18)

ver8 분석을 원고용으로 다시 배치한 그림판입니다. 새 통계 분석은 없고, 검산을 마친 ver8 출력(`../fig7_ver8_main_supp_20261117/`)과 공개 코호트 v8.1 출력(`../fig7_ver8_public_replication_20261118/`)을 읽어 그림 5개를 그립니다. 코호트 구성 수치(채혈 수, 센터, 스테로이드)만 원본 메타데이터에서 직접 셉니다. 결과·범례 원고는 `docs/Results_and_Legends_ver9.md`입니다.

## 그림 구성

| 그림 | 패널 |
| --- | --- |
| 본 그림 7 (a–q) | a–b 코호트 구성과 중증도·사망 확장(이전 Fig. 0a–b), c CLL 분자 스크린, d 5단계 퍼널, e 22개 세포형 스크린(이전 보충 2a), f n80(이전 보충 2b), g FDR 이득 지도, h 신뢰도 스코어카드, i 효과 크기, j–o 환자별 normal score, p HSPC 사망 강건성, q n = 20 검정력 |
| 보충 그림 1 (a–d) | 이전 Fig. 0: 채혈 타임라인(217회), 채혈 횟수·센터·스테로이드, 코호트 특성, 2단계 검증 구조 |
| 보충 그림 2 (a–o) | 이전 보충 1: CLL 영향, 앵커 4개, NK 분해, 27명 모달리티 스케일, raw % 패널 |
| 보충 그림 3 (a–c) | 이전 보충 2: 검정력 곡선, CD14 mono IFN step 4 점검, 14가지 민감도 |
| 보충 그림 4 (a–d) | 이전 보충 3: 공개 코호트 4개 재현(N = 219) |

## 파일

| 파일 | 내용 |
| --- | --- |
| `figure/build_ver9_all_figures.py`, `figure/style.py` | 그림 5개와 패널별 PNG/PDF/SVG를 그립니다. 원래 판에서 그림 라벨 2개만 고쳤습니다(아래 검산 참고) |
| `PBISC_Figure7_ver9_verification_20261118.ipynb` | 검산 노트북. ver9에서 새로 나온 수치를 원본 메타데이터로 다시 세고, Results의 공개 코호트 수치를 대조하고, 공개 코호트 독립 재구현을 실행합니다(FAIL은 원래 판의 오류). 저장소 규칙에 따라 출력은 지웠습니다 |
| `docs/Results_and_Legends_ver9.md` | 결과·범례 원고. n80 문장을 고치고 4절에 검산 메모를 덧붙인 판입니다 |

## 실행

1. 입력(저장소에 없음): `F:\PBISC\eskd_covid_clinical_assoc_20261115\results\sample_level_composition_metadata.csv`, `F:\PBISC\paper_claim_analysis\covid_sample_patient_id_reconciliation.csv`, ver8 `data/` CSV 13개(`15_Figure7_ver8_Main_and_Supp\data`), 공개 코호트 `data/` CSV 4개(`16_Public_COVID_Replication_ver8\data`).
2. 스크립트의 `V8_DIR`, `PUB_DIR`, `OUT_DIR`, `ROOT`를 자기 환경에 맞게 바꾸고 `python build_ver9_all_figures.py`를 실행합니다.

## 검산 (2026-11-18)

- **입력:** ver8 CSV는 2026-11-17 검산 이후 바뀌지 않았습니다. 공개 코호트 CSV는 v8.1로 다시 만들어졌는데, 재실행(최대 차이 3.4×10⁻¹⁴)과 독립 재구현(코호트별 25행, 메타분석 7개)으로 확인했습니다.
- **그림:** 원래 스크립트를 원본 입력으로 다시 돌린 결과, 그림 5개의 SVG 안 글자(수치·라벨 1,194개)가 원래 그림과 모두 같았습니다.
- **새 수치 확인:** Real 37명·138회 채혈, PBISC new 29명·79회, 합계 217회, 공유 27명·Real 단독 10명, WHO 3 5 → 12명, 사망 7 → 13 / 5 → 11 / 3 → 9, 스테로이드 노출 Real 7명·34회, 채혈 횟수 분포, 채혈 기간, 앵커 n80 19–24, 공개 코호트 메타 수치가 모두 맞았습니다.
- **고친 곳 4건 (중대한 오류 없음, 결론 영향 없음):**
  1. Results의 n80 문장: Real 단독 n80은 "26 – n.r."이 아니라 **31 – n.r.**이고, PBISC 추가 후 T CD8 CM 사망은 **≤ 10이 아니라 26**입니다(그림 7f·q와 일치하도록 문장 수정).
  2. 같은 내용의 1절 Fig. 7f 설명도 고쳤습니다.
  3. 그림 7a "Real-only … (Wave 2)": Real 단독 10명은 Cambridge 6명·NCL 4명이라 웨이브로 묶이지 않습니다 → "(Cambridge 6, NCL 4)".
  4. 보충 그림 1c "17/88 NCL draws <2%"는 음성·회복기 샘플을 포함한 전체 Real 기준입니다. 분석 코호트 기준인 **15/60**으로 고쳤습니다.
- **참고:** PBISC new의 스테로이드 노출은 기록이 없고, "0%"는 1차 유행 시기에서 추정한 값입니다. 스테로이드 표지인 CD14 mono dex가 PBISC new에서 사실상 0(중앙값 0.000%, Real 0.037%, Real 최대 5.2%)이라 데이터와 맞습니다. 그림 7a의 "med. 71 yr"는 코드에 고정된 값이며, 분석 환자 66명의 중앙값은 71.5세입니다.
