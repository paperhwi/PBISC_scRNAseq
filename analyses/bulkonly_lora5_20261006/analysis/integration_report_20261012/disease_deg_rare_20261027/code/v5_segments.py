# -*- coding: utf-8 -*-
# New / replaced segments for deck v4 (2026-10-30). Imported as text by make_v4.py.
NEW_TAG2 = 's.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/31 추가", RED)\n'
NEW_TAG = 's.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)\n'

TITLE = r'''# ========================= 1 title
s = new(); box(s, 0, 0, 0.17, 7.5, BLUE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 1.92, 11.9, 0.3, "PBISC floor read-out · 질병 DEG 및 희귀 집단 재검증 · v5", size=14, color=BLUE, bold=True, spc=150)
title_box(s, "생성 세포는 질병 신호를 어디까지 가져오고,\n어디서 체계적으로 놓치는가", y=2.5, size=40)
text(s, 0.89, 4.15, 11.9, 0.8, "질병 DEG 재계산 · 매칭 n 검정 · 희귀 집단 내부 발현 · DC 과소생성 원인 · 단핵구/대식세포 subset · "
     "marker 재annotation · scSemiProfiler 동일 지표 비교 (패혈증 Ye 2026 / COVID Imperial ESKD)", size=18, color=BODY, line=1.15)
pill(s, 0.89, 5.25, 1.75, "PBISC floor", BLUE); pill(s, 2.8, 5.25, 2.1, "PBISC prev K4000", RED)
pill(s, 5.06, 5.25, 2.25, "Real (held-out)", GREY); pill(s, 7.47, 5.25, 2.6, "real-real 상한선", INK); pill(s, 10.23, 5.25, 2.2, "bulk 입력", PURPLE)
text(s, 0.89, 6.82, 11.9, 0.26, "Jong Hwi Kim · Seoul National University · 2026년 10월 31일 · 약 50분 (부록 별도)", size=12, color=MUTED)
notes(s, """[0:40] 표지
"v5는 v4에 비교 모델 scSemiProfiler를 같은 지표로 평가한 결과를 더한 판입니다. v3 이후 들어온 질문 일곱 가지에 답합니다. 질병 DEG를 case·control 수를 맞춰 다시 본 것, 림프구 볼케이노, DC가 왜 적게 생성되는가, floor가 어느 bulk 샘플에서 왔는가, 단핵구·대식세포 subset 지형, 그리고 real marker로 floor를 직접 annotation하면 달라지는가입니다."
"10/30에 추가된 슬라이드는 주황색, 10/31에 추가된 scSemiProfiler 슬라이드는 빨간색 표시가 오른쪽 위에 있습니다. 나머지는 v3와 같고, 요약·종합·제안·한계·논의만 갱신했습니다."
"모든 수치는 로컬 데이터로 새로 계산했고 새로 생성한 세포는 없습니다." """)
'''

SUMMARY = r'''# ========================= 2 summary
s = new(); header(s, "요약", "신호는 bulk에 있다 — 모델이 세포형으로 보내지 못한다", NEXT())
three_cards(s, 1.75, 3.0, [
 ("1", "방향은 상한의 2/3", f"패혈증 T {brk('sepsis','T cell')} · 단핵구 {brk('sepsis','Monocyte')}, COVID 단핵구 {brk('covid','Monocyte')}.\n"
      "case·control 수를 맞춰도 림프구 logFC r 0.59~0.66 (패혈증) 유지."),
 ("2", "원인을 특정했다", f"C1Q 모듈은 bulk 입력에 온전히 존재(bulk/real {A2[A2.module=='C1Q / complement'].eval('bulk_logFC/real_logFC').median():.2f})."
      f"\n대조 사다리: read-out 규칙 −31% → 나머지 −81%가 모델. DC는 자기 bulk 입력의 pDC 신호를 따라가지 않는다(r −0.18)."),
 ("3", "법칙으로 일반화", f"세포형 전용성이 높을수록 회복률이 낮다 (유전자 ρ −0.25, 모듈 ρ = {stats.spearmanr(BM.tau,BM.floor_ratio).statistic:.2f}). "
      "real marker로 재annotation해도 회복되지 않는다.")],
 accents=[BLUE, RED, ORANGE])
text(s, 0.89, 5.0, 11.55, 1.75, [
  [("한 문장  ", {"bold": True, "color": INK}),
   ("PBISC는 전신 프로그램은 옮기고 세포형 전용 프로그램은 옮기지 못한다. 저발현·깊이·read-out·IEG·라벨링 방식 모두 주원인이 아니다 — "
    "bulk 조건에서 세포형으로의 귀속이 병목이다.", {})],
  [("10/30 추가  ", {"bold": True, "color": ORANGE}),
   ("① 매칭 n: up-DEG는 유지, down-DEG는 소실, 생성 샘플 간 SD는 real의 0.61~0.73× ② DC 결손은 대조군에서 발생(HC 0.64×, ICU 0.33×) "
    "③ 패혈증 real 단핵구의 23%는 환자 고유 영역 — floor에 없음 ④ marker 재annotation은 비율 오차를 키움 (패혈증 0.35 → 1.53)", {})],
  [("10/31 추가  ", {"bold": True, "color": RED}),
   ("scSemiProfiler를 같은 지표로 보면 logFC r 0.05~0.33, 샘플 간 SD 2.5~3.1×, 희귀 subset 18/27 결손, 전용성 기울기 없음 — 전역적으로 실패. floor의 결함은 선택적이다.", {})]],
  size=12, color=BODY, line=1.1, after=5)
notes(s, """[1:50] 요약
"세 카드는 v3와 같은 결론이고, 오늘 추가 분석이 각 카드를 보강합니다."
"첫째, 숫자. case와 control 샘플 수를 real과 floor에서 똑같이 맞춰도 패혈증 림프구 logFC r이 0.59에서 0.66으로 유지됩니다. 즉 v3의 상관값은 floor 쪽 샘플이 많아서 나온 게 아닙니다."
"둘째, 원인. C1Q에 이어 DC에서도 같은 구조를 확인했습니다. 같은 샘플의 bulk 입력에서 pDC 마커 점수와 생성 결과의 pDC 점수 상관이 real에서는 0.57인데 floor는 마이너스 0.18입니다. 모델이 자기 입력을 따라가지 않습니다."
"셋째, 법칙. 라벨링 방식을 바꿔도 회복되지 않는다는 걸 확인했습니다. real의 Top10 marker로 floor를 직접 annotation하면 오히려 비율 오차가 커집니다."
"아래 주황 줄이 10/30에 추가된 네 가지, 빨간 줄이 10/31에 추가된 비교 모델 결과입니다. 각각 뒤 슬라이드에서 설명드리겠습니다." """)
'''

PROV = r'''# ========================= 4b provenance (NEW)
s = new(); header(s, "방법 · 10/30", "floor는 어느 bulk 샘플에서 왔는가", NEXT(), "prov/floor_source_samples.csv · 샘플당 4,000세포 생성")
''' + NEW_TAG + r'''rows = [("그룹 (패혈증)", "real scRNA 샘플", "PBISC floor 샘플", "floor 출처"),
        ("adult sepsis", "93", "68", "GSE279448 bulk-only"),
        ("adult healthy control", "7", "36", "GSE279448 bulk-only"),
        ("adult ICU control", "8", "17", "GSE279448 bulk-only"),
        ("paediatric ICU control", "4", "1 (PIC06)", "GSE279448 bulk-only"),
        ("합계 (floor 대상 그룹)", "112", "122", "real scRNA와 샘플 중복 0")]
table(s, rows, 0.9, 1.72, [2.9, 2.0, 2.0, 2.8], h_row=0.42, bold_cols=(0,), font=12, head_font=12, hl_rows=(2,))
side_cards(s, 10.95, 1.72, 1.65, [("122", "패혈증 floor\n= bulk-only\n전부"), ("0", "real scRNA와\n겹치는 샘플")], h=1.35, accent=[BLUE, INK], body=11, head=20)
text(s, 0.89, 4.65, 11.55, 2.1, [
  [("COVID/ESKD  ", {"bold": True, "color": INK}),
   ("floor 164샘플 = 환자 76명 (positive 143 · negative 19 · recovery 2). 이 중 ", {}),
   ("32명은 real scRNA가 다른 방문으로 존재", {"bold": True, "color": RED}),
   ("하고 44명은 scRNA가 전혀 없다. 같은 방문(샘플)이 겹치는 경우는 0. adapter 학습 donor 5샘플은 DEG에서 제외.", {})],
  [("왜 중요한가  ", {"bold": True, "color": INK}),
   ("패혈증 floor는 real과 환자가 완전히 다르다 → real-floor 차이에는 '생성 오차'와 '환자 구성 차이'가 함께 들어간다. "
    "real B(다른 환자, 같은 라벨 전이)를 대조로 둔 이유가 이것이다. 또한 floor는 HC를 36샘플 갖지만 real은 7샘플 → "
    "샘플 수가 다르면 DEG 수·유의성이 달라지므로 다음 슬라이드들에서 n을 맞춰 다시 본다.", {})]],
  size=13, color=BODY, line=1.15, after=6)
notes(s, """[1:40] floor 출처 (추가)
"floor가 어느 bulk 샘플에서 추론됐는지 메타데이터로 추적했습니다."
"패혈증은 깔끔합니다. GSE279448 bulk 206샘플 중 scRNA가 함께 있는 84개를 빼고 남은 bulk-only 122개 전부를 생성했습니다. 그래서 real scRNA와 겹치는 샘플이 하나도 없습니다. 그룹별로 성인 패혈증 68, 성인 건강 대조 36, 성인 ICU 대조 17, 소아 ICU 대조 1입니다."
"같은 그룹의 real scRNA는 93, 7, 8, 4입니다. 표에서 파랗게 칠한 건강 대조 행을 봐주십시오. real은 7샘플, floor는 36샘플입니다. 다섯 배 차이입니다."
"COVID는 환자 76명의 164샘플입니다. 이 중 32명은 다른 방문에서 real scRNA가 있습니다. 같은 방문이 겹치지는 않지만, 같은 사람이 양쪽에 있다는 점은 리뷰어에게 명시해야 합니다."
"이 슬라이드의 결론은 두 가지입니다. 패혈증에서 real과 floor의 차이에는 환자 구성 차이도 섞이므로 real B 대조가 필요하고, 샘플 수가 크게 다르므로 n을 맞춘 비교가 필요합니다." """)
'''

LYMPH = r'''# ========================= 9b lymphoid volcano (NEW)
s = new(); header(s, "결과 2 · 10/30", "림프구 볼케이노 — 방향은 맞고, 유의성은 부풀려진다", NEXT(),
       "패혈증 성인 sepsis vs HC · 위 real, 아래 floor · 주황 = IEG/ambient")
''' + NEW_TAG + r'''pic(s, "volcano_lymphoid_sepsis", 0.75, 1.62, w=8.75)
side_cards(s, 9.75, 1.7, 2.85, [
    ("DEG 수는 비슷", "T 104 → 75\nB 139 → 136\nNK 118 → 105\nIEG 비중 27~34%로 양쪽 동일"),
    ("유의성 과장", "floor −log10 FDR 25~30\nvs real 8~20.\nfloor HC 36샘플 + 샘플 간\nSD가 작아서 (→ 매칭 n)"),
    ("개별 불일치", "B: JCHAIN real 비DEG → floor +2.4\nB: TNF real +2.8 → floor 0.7\nS100A8: 림프구로 번짐")], h=1.6, accent=[GREY, RED, ORANGE], body=11.5, head=14)
notes(s, """[2:10] 림프구 볼케이노 (추가)
"v3 10번 슬라이드의 단핵구·DC 볼케이노를 T, B, NK로 확장했습니다. 위가 real, 아래가 floor, 같은 x축입니다."
"첫째, 방향은 맞습니다. NR4A1, AREG, CDKN1A, SOCS3 같은 상향 유전자와 ZNF683, SLC4A10 같은 하향 유전자가 같은 쪽에 있습니다. DEG 개수도 T 104 대 75, B 139 대 136, NK 118 대 105로 비슷합니다."
"둘째, 주황 점을 봐주십시오. immediate-early와 ambient 유전자입니다. 양쪽 모두 DEG의 27에서 34퍼센트를 차지합니다. v3 결과 6에서 말씀드린 IEG 문제가 림프구 전반에서 그대로 보입니다."
"셋째, y축입니다. floor는 마이너스 로그 FDR이 25에서 30까지 올라가는데 real은 8에서 20입니다. 같은 효과 크기에서 훨씬 유의하게 나옵니다. 이유는 두 가지입니다. floor는 건강 대조가 36샘플이고 real은 7샘플입니다. 그리고 생성 샘플끼리 너무 비슷합니다. 다음 슬라이드에서 n을 맞춰 확인합니다."
"넷째, 개별 불일치입니다. B세포의 JCHAIN은 real에서는 변화가 없는데 floor에서는 플러스 2.4로 DEG입니다. 형질세포 프로그램이 B세포로 섞여 들어온 것으로 보입니다. 반대로 TNF는 real에서 2.8인데 floor에서는 0.7입니다."
"COVID 림프구 볼케이노는 부록에 있습니다. ISG는 같은 방향으로 보존되지만 IFI27이 눌리고, B세포에서 약한 상향 DEG가 106개로 부풀려집니다." """)
'''

MATCHED = r'''# ========================= 9c matched n (NEW)
s = new(); header(s, "결과 2b · 10/30", "n을 맞추면 — up은 유지, down은 소실, 샘플은 너무 비슷", NEXT(),
       "같은 case/control 수로 subsample × 50회 · 오른쪽 = 샘플 간 SD 비")
''' + NEW_TAG + r'''pic(s, "matchedN_summary", 1.67, 1.62, w=10.0)
rows = [("매칭 n (50회 평균)", "패혈증 T", "패혈증 B", "패혈증 NK", "COVID T", "COVID B", "COVID NK"),
        ("logFC r", "0.66", "0.59", "0.60", "0.57", "0.44", "0.39"),
        ("down-DEG real → floor", "49 → 22", "46 → 41", "51 → 37", "2.9 → 0.04", "5.4 → 1.1", "15.6 → 0.6"),
        ("7/7 균형 설계 정밀도", "0.15", "0.04", "0.08", "0.27", "0.03", "0.22")]
table(s, rows, 0.9, 5.3, [2.6, 1.5, 1.5, 1.5, 1.5, 1.5, 1.45], h_row=0.32, bold_cols=(0,), font=11, head_font=11, warn_rows=(2, 3))
notes(s, """[2:30] 매칭 n (추가)
"질문은 '샘플 수를 똑같이 맞추면 어떻게 되나'였습니다. real과 floor를 같은 case·control 수로 50번 무작위 추출해 DEG를 다시 계산했습니다. 패혈증은 68 대 7, COVID는 환자 37~39 대 18입니다."
"왼쪽 위로 향한 막대가 상향 DEG입니다. 패혈증 T 60 대 60, B 92 대 109, NK 73 대 80으로 floor가 유지하거나 조금 더 만듭니다. COVID B는 31 대 95로 세 배 부풀려집니다."
"아래로 향한 빗금 막대가 하향 DEG입니다. 패혈증 T가 49에서 22로 절반이 되고, COVID는 사실상 0입니다. 하향 신호가 체계적으로 사라집니다."
"오른쪽이 원인입니다. 생성 샘플들 사이의 표준편차가 real의 0.61에서 0.73배입니다. 생성 샘플들이 서로 너무 비슷하다는 뜻입니다. 패혈증 대조군만 예외인데 real 대조군이 7명뿐이라 추정이 불안정합니다."
"표 마지막 줄이 경고입니다. case와 control을 7 대 7로 맞추면 floor는 real보다 DEG를 2에서 5배 더 부르고, 정밀도는 0.03에서 0.27로 떨어집니다. 분산이 작아서 가짜 DEG가 쉽게 나온다는 뜻입니다."
"좋은 소식은 logFC r입니다. n을 맞춰도 패혈증 0.59에서 0.66으로 전체 n의 0.66에서 0.69와 거의 같습니다. 그러니까 v3의 상관 결과는 n 차이 때문이 아닙니다."
"원고 함의는 두 가지입니다. DEG 개수나 FDR을 비교할 때는 반드시 매칭 n으로 보고하고, 하향 DEG는 따로 보고해야 합니다." """)
'''

DC = r'''# ========================= 15b DC underproduction (NEW)
s = new(); header(s, "결과 6b · 10/30", "DC 과소생성 — 모델이 bulk의 pDC 신호를 따르지 않는다", NEXT(),
       "b: 같은 샘플의 bulk vs 생성 pseudobulk 마커 점수 (COVID real은 매칭 없음)")
''' + NEW_TAG + r'''pic(s, "v4_dc_underproduction", 0.85, 1.7, w=11.65)
text(s, 0.89, 5.3, 11.55, 1.5, [
  [("a  ", {"bold": True, "color": INK}), ("결손은 대조군에서: HC 0.64×, ICU 0.33×, sepsis 1.11× → 질병에 따른 DC 감소(real HC 1.11% → sepsis 0.29%)가 압축된다. COVID는 비율 자체는 맞다(0.9~1.1×).", {})],
  [("b  ", {"bold": True, "color": INK}), ("T·B·NK·단핵구는 입력을 따라가지만(r 0.67~0.87) pDC는 −0.18 / 0.06, cDC 0.26 / 0.24 (real 0.57 / 0.69). ", {}),
   ("bulk에서 pDC 마커는 전체 CPM의 0.001%, cDC 0.004%", {"bold": True, "color": RED}), (" — 단핵구(0.49%)의 1/100~1/400.", {})],
  [("c  ", {"bold": True, "color": INK}), ("pDC 마커를 가진 세포 중 2개 이상을 함께 켜는 비율: COVID real 31~36% vs floor 11%. 마커는 흩어져 있고 pDC로 조립되지 않는다.", {})]],
  size=12, color=BODY, line=1.1, after=4)
notes(s, """[2:40] DC 과소생성 원인 (추가)
"v3에서 패혈증 DC가 회수율 1.6퍼센트로 사실상 실패했다고 말씀드렸습니다. 왜 그런지 다섯 가지 분석을 했고, 세 가지로 요약했습니다."
"a. 어디서 부족한가. 그룹별로 나누면 패혈증 환자 샘플에서는 floor가 오히려 1.1배로 맞습니다. 부족한 건 대조군입니다. 건강 대조 0.64배, ICU 대조 0.33배입니다. real에서는 건강 대조 1.1퍼센트에서 패혈증 0.3퍼센트로 DC가 크게 줄어드는데, floor는 대조군을 적게 만들어서 이 감소폭이 눌립니다. 질병 DEG에서 DC가 실패한 직접적인 이유입니다."
"b. 모델이 입력을 따라가는가. 같은 샘플의 bulk에서 계통별 마커 점수를 재고, 생성된 세포의 pseudobulk 점수와 상관을 봤습니다. T, B, NK, 단핵구는 0.67에서 0.87로 잘 따라갑니다. 그런데 pDC는 마이너스 0.18입니다. real scRNA와 bulk 사이는 0.57입니다. 모델이 자기 입력의 pDC 신호를 무시하고 있습니다."
"이유는 신호 크기로 보입니다. bulk 전체 CPM에서 pDC 마커가 차지하는 비중이 0.001퍼센트입니다. 단핵구 마커의 400분의 1입니다. 조건 벡터에서 이 신호가 묻힌다고 해석합니다."
"c. 세포 단위로 보면, pDC 마커를 하나라도 가진 세포 중 두 개 이상을 같이 켠 비율이 COVID real은 31에서 36퍼센트인데 floor는 11퍼센트입니다. 마커가 여러 세포에 흩어져 있고 한 세포에 모여 pDC를 이루지 않습니다. v3에서 말한 '정체성 프로그램 희석'의 DC 버전입니다."
"정리하면, DC 과소생성은 라벨 문제가 아니라 bulk 내 신호가 너무 작아 모델이 그 세포형으로 귀속하지 못하는 문제입니다. 전용성 법칙의 극단 사례입니다." """)
'''

MYE_S = r'''# ========================= 20c myeloid sepsis (NEW)
s = new(); header(s, "결과 12 · 10/30", "단핵구 subset — 흔한 상태는 겹치고, 환자 고유 상태는 없다", NEXT(),
       "패혈증 · 단핵구 구획 joint Harmony UMAP · ● real ▲ floor")
''' + NEW_TAG + r'''pic(s, "umap_myeloid_sepsis", 1.0, 1.6, w=11.3)
text(s, 0.89, 6.1, 11.55, 0.7, [[("읽는 법  ", {"bold": True, "color": INK}),
     ("Mono_c06-MKI67 221 → 3, Macrophage 400 → 60 (0.25×). real 세포의 23%는 주변 30개 이웃의 90% 이상이 real인 '환자 고유 영역'에 있고 "
      "(Abd-S123·IC16·Skin-S140·Abd-S51), 대식세포 400개 중 225개가 이 영역에 속한다. HC 세포는 0.", {})]], size=12, color=BODY, line=1.12)
notes(s, """[2:00] 단핵구·대식세포 subset 패혈증 (추가)
"v3 슬라이드 5, 6의 희귀 집단 UMAP과 같은 형식으로 단핵구 구획만 따로 그렸습니다. joint Harmony 공간을 그대로 쓰고 단핵구 세포만 남겨 UMAP을 다시 계산했습니다."
"왼쪽 소스 패널을 보면 대부분 겹치는데, 아래쪽 꼬리와 왼쪽 위에 회색만 있는 영역이 있습니다."
"오른쪽 패널에서 증식 단핵구는 real 221개 대 floor 3개, 대식세포는 400 대 60입니다. CD16과 CD14CD16 단핵구는 floor 삼각형이 real 원과 같은 자리를 덮습니다."
"회색만 있는 영역을 정량했습니다. real 세포 주변 30개 이웃의 90퍼센트 이상이 real인 경우를 세면 real의 23퍼센트입니다. 이 세포들은 네 명의 환자에 몰려 있고 전부 패혈증 또는 ICU 대조입니다. 건강 대조는 0입니다. 대식세포 400개 중 225개가 이 영역에 있습니다."
"해석은 신중해야 합니다. 패혈증 floor는 real과 환자가 전혀 겹치지 않으니, 이 영역 일부는 환자 구성 차이입니다. 그래도 '환자 고유의 극단적 골수계 상태를 bulk에서 생성하지 못한다'는 건 분명합니다." """)
'''

MYE_C = r'''# ========================= 20d myeloid covid (NEW)
s = new(); header(s, "결과 12b · 10/30", "COVID 단핵구 — 잘 섞이지만 dex·C1 상태가 빠진다", NEXT(),
       "COVID/ESKD · 같은 방식 · 오른쪽 범례 괄호 = n real / n floor")
''' + NEW_TAG + r'''pic(s, "umap_myeloid_covid", 1.0, 1.6, w=11.3)
text(s, 0.89, 6.1, 11.55, 0.7, [[("비교  ", {"bold": True, "color": INK}),
     ("환자 고유 영역은 1.2%로 거의 없다(패혈증 23%). 비율로는 CD14mono 0.98×, IFN 0.94×로 맞지만 CD14mono_dex 57 → 2 (0.03×), "
      "CD16mono_C1 0.38×, 반대로 CD16mono는 1.47× 과생산 — 흔한 상태로 쏠린다. (subset 비율 막대는 부록)", {})]], size=12, color=BODY, line=1.12)
notes(s, """[1:40] COVID 단핵구 (추가)
"COVID에서는 지형이 훨씬 잘 섞입니다. real만 있는 영역이 1.2퍼센트입니다. 일부 환자가 floor 쪽에도 다른 방문으로 들어가 있는 것도 이유일 수 있습니다."
"비율로 보면 흔한 subset은 맞습니다. CD14 단핵구 0.98배, IFN 단핵구 0.94배입니다."
"빠지는 건 두 가지입니다. 덱사메타손 처리 단핵구는 57개 대 2개, C1Q를 켜는 CD16 단핵구 C1은 0.38배입니다. C1은 v3에서 본 C1Q 모듈 결손과 같은 이야기입니다."
"반대로 일반 CD16 단핵구는 1.47배로 과생산됩니다. 희귀 상태가 빠진 자리를 흔한 상태가 채웁니다. 패혈증에서도 RETN과 CD16이 1.4배였습니다." """)
'''

MK1 = r'''# ========================= 20e marker re-annotation 1 (NEW)
s = new(); header(s, "결과 13 · 10/30", "marker 재annotation — floor는 real에 가까워지지 않는다", NEXT(),
       "real A Top10 marker → de novo Leiden + marker 점수 · b,c = 구분 가능 subset")
''' + NEW_TAG + r'''pic(s, "v4_marker_ruler", 0.85, 1.66, w=11.65)
text(s, 0.89, 5.3, 11.55, 1.5, [
  [("a  자 자체가 약하다  ", {"bold": True, "color": INK}),
   ("같은 marker를 저자 라벨이 있는 real B에 적용하면 6-type 0.84/0.92지만 fine subset 정확도 0.31/0.35. 저자 subcluster의 약 1/3(9/46, 14/45)만 marker로 재현된다.", {})],
  [("b  비율  ", {"bold": True, "color": INK}),
   ("구분 가능한 subset에서도 floor marker 라벨은 오차가 커진다 (패혈증 0.35 → 1.53, COVID 0.69 → 1.10). real B는 같은 자로 0.42 / 0.79 — 자 탓이 아니다.", {})],
  [("c  발현  ", {"bold": True, "color": INK}),
   ("marker에 쓰지 않은 HVG의 identity r: floor kNN 0.91/0.95 → marker 0.88/0.92. 개선 없음.", {})]],
  size=12, color=BODY, line=1.1, after=4)
notes(s, """[2:30] marker 기반 재annotation (추가)
"질문은 이랬습니다. 희귀 집단 재현을 kNN 라벨 전이로 봤는데, real이 주장하는 subcluster의 Top10 marker로 floor를 직접 군집하고 annotation하면 비율이나 발현 유사도가 비슷해지지 않느냐."
"방법은 real A의 저자 라벨에서 6개 세포형 marker와 세포형 안의 subset marker를 각각 Top10으로 뽑고, floor와 real B를 각자 처음부터 군집한 뒤 marker 점수가 가장 높은 subset을 붙였습니다."
"a가 가장 중요한 대조입니다. 같은 marker 자를 저자 라벨이 있는 real B에 대보면 6개 세포형은 0.84에서 0.92로 맞지만, fine subset은 0.31에서 0.35만 맞습니다. 진짜 real에서도 저자 subcluster의 3분의 1만 Top10 marker로 다시 찾아집니다. CD8, NKT, NK subset이 GZMH, NKG7, CCL5를 공유하기 때문입니다."
"그래서 공정 비교를 위해 marker로 구분 가능한 subset, 즉 real B에서 F1이 0.5 이상인 subset만 썼습니다."
"b, 비율 오차입니다. real B는 marker로 annotation해도 0.42, 0.79로 저자 라벨과 거의 같습니다. floor는 kNN 0.35에서 marker 1.53으로 네 배 나빠집니다."
"c, marker에 쓰지 않은 유전자에서의 발현 identity도 0.91에서 0.88로 오히려 떨어집니다." """)
'''

MK2 = r'''# ========================= 20f marker re-annotation 2 (NEW)
s = new(); header(s, "결과 13b · 10/30", "marker 라벨이 희귀 subset을 '늘리는' 방식 — 이웃 세포를 끌어온다", NEXT(),
       "F1 = 같은 marker 자를 real B에 적용한 저자 라벨 회복력")
''' + NEW_TAG + r'''pic(s, "v4_marker_keysubsets", 0.95, 1.6, w=11.4)
text(s, 0.89, 5.82, 11.55, 1.0, [
  [("부풀려진 것  ", {"bold": True, "color": RED}),
   ("패혈증 Macrophage 3.76% (real A 0.69%) = 대부분 CD16/CD14 단핵구 · COVID CD16mono_C1 4.44% = 2,667개 중 2,150개가 CD16mono · "
    "패혈증 pDC = HSC 유사 세포, identity r 0.31.  ", {}),
   ("없는 것은 그대로  ", {"bold": True, "color": INK}),
   ("MKI67 계열·CD14mono_dex·B_ASC_dividing은 marker로도 0개.", {})]],
  size=12, color=BODY, line=1.12)
notes(s, """[1:50] 희귀 subset 사례 (추가)
"점 하나하나가 subset 비율입니다. 검은 원이 real A, 진한 파랑 삼각형이 floor marker 라벨입니다."
"패혈증 대식세포를 보십시오. marker로 붙이면 3.76퍼센트로 real A의 0.69퍼센트보다 다섯 배 많아집니다. 그런데 그 세포들의 kNN 라벨을 보면 대부분 CD16과 CD14 단핵구입니다. C1Q가 아니라 FCGR3A, FTL, FTH1 같은 공유 marker로 점수가 올라간 겁니다. COVID CD16mono_C1도 같습니다. 2,667개 중 2,150개가 일반 CD16 단핵구입니다."
"패혈증 pDC는 marker로 167개가 잡히는데 대부분 HSC 유사 세포이고, pDC identity 상관이 0.31입니다."
"반대로 증식 계열, dex 단핵구, 분열 형질모세포는 marker로도 0개입니다. 만들어지지 않은 상태는 어떤 라벨링으로도 나타나지 않습니다."
"결론은 이렇습니다. 라벨링은 경계만 옮길 뿐 생성 세포에 없는 identity 신호를 만들지 못합니다. kNN 전이가 floor에 오히려 유리하고 보수적인 방법입니다. 리뷰어가 marker 기반 재annotation을 요구하면 real B 검증과 함께 이 결과를 제시하면 됩니다." """)
'''

SYN = r'''# ========================= 21 synthesis
s = new(); header(s, "종합", "원인이 어디에 있는지 말할 수 있게 되었다", NEXT())
rows = [("단계", "확인한 것", "근거", "판정"),
        ("bulk 입력", "모듈이 입력에 존재하는가", "C1Q bulk/real 0.86, 5/5 FDR<0.001", "정보 있음"),
        ("깊이", "시퀀싱 깊이 차이인가", "깊이 맞춘 real 1.00", "원인 아님"),
        ("read-out (Poisson)", "출력 규칙 일반의 문제인가", "0.84", "거의 영향 없음"),
        ("read-out (floor 규칙)", "floor 규칙 고유의 손실", "0.69", "−31%"),
        ("생성 모델", "나머지", "C1Q 0.12 · pDC 입력 추종 r −0.18", "−81% (지배적)"),
        ("표현 방식", "support인가 magnitude인가", "양성 세포 1.07× · 양성 내 양 0.11×", "magnitude"),
        ("샘플 수 (10/30)", "n 차이가 결과를 만드는가", "매칭 n에서 r 0.59~0.66 유지 · SD 0.61~0.73×", "r 유지, down 소실"),
        ("라벨링 (10/30)", "kNN 전이 탓인가", "marker 재annotation 시 비율 오차 0.35→1.53", "원인 아님")]
table(s, rows, 0.9, 1.68, [2.4, 3.3, 4.2, 2.0], h_row=0.39, bold_cols=(0,), font=11, head_font=11, warn_rows=(5, 6))
text(s, 0.89, 5.45, 11.55, 1.35, [
  [("함의  ", {"bold": True, "color": INK, "size": 13.5})],
  [("남은 결함은 ", {}), ("세포형 전용 모듈의 magnitude", {"bold": True, "color": RED}),
   ("와 ", {}), ("샘플 간 변동의 압축", {"bold": True, "color": RED}),
   ("이다. 게이트와 수정안을 그쪽(A1/C3 count·NB head, 샘플 간 분산 보존)으로 옮겨야 하며, 라벨링·깊이·n 조정으로는 해결되지 않는다.", {})]],
  size=12.5, color=BODY, line=1.15, after=5)
notes(s, """[2:30] 종합
"v3의 원인 분리 표에 오늘 두 줄을 더했습니다."
"위 여섯 줄은 같습니다. bulk에 정보가 있고, 깊이와 Poisson read-out은 원인이 아니며, floor 규칙이 31퍼센트, 나머지 81퍼센트가 모델입니다. 모델 줄에 DC 결과를 붙였습니다. 모델이 자기 입력의 pDC 신호를 따라가지 않습니다."
"아래 두 줄이 새로 들어간 것입니다. 샘플 수를 맞춰도 상관은 유지되므로 n 차이가 결과를 만든 건 아닙니다. 다만 하향 DEG가 사라지고 생성 샘플 간 분산이 압축되어 있다는 새 결함이 드러났습니다."
"라벨링도 원인이 아닙니다. marker로 다시 붙이면 오히려 나빠집니다."
"그래서 함의에 한 가지를 더합니다. magnitude뿐 아니라 샘플 간 분산 보존도 게이트와 수정안의 대상이 되어야 합니다." """)
'''

ACT = r'''# ========================= 22 actions
s = new(); header(s, "제안", "원고와 다음 실험에 반영할 것", NEXT())
rows = [("#", "항목", "근거", "상태"),
        ("1", "Result 6B sign agreement · top-100 정밀도 분모 명시 또는 재계산 값 교체", "분모 없이 0.90은 방어 불가 (0.5 바닥 지표)", "원고 수정"),
        ("2", "달성률을 구간으로 보고 (상한선 정의 2종 병기)", "패혈증 T 64% ↔ 91%", "원고 수정"),
        ("3", "게이트 M7 '모듈 magnitude ≥ 0.5×' 추가", "C1Q 0.11× — 기존 M1~M6 통과해도 안 잡힘", "지표 산출됨"),
        ("4", "수정안 우선순위 B1(support) → A1/C3(count·NB head)", "양성 비율 1.07×, 양성 내 양 0.11×", "방향 전환"),
        ("5", "모듈 집합을 데이터 주도 정의로 고정", "유전자 ρ −0.25, 모듈 ρ −0.46", "완료"),
        ("6", "패혈증 DC는 '생성 실패'로 명시 + 원인(대조군 결손·입력 미추종) 서술", "recall 0.016 · HC 0.64× · pDC 추종 r −0.18", "원고 수정"),
        ("7", "COVID 수치는 CI와 함께 보고", "COVID T r 0.58 [0.10, 0.71]", "원고 수정"),
        ("8", "패혈증 림프구 IEG 제외값 병기, 유전자 목록 헤드라인 제외", "IEG가 DEG의 26~34%", "완료"),
        ("9", "DEG 개수·FDR 비교는 매칭 n으로, 하향 DEG는 별도 보고", "down 49→22 (T), COVID ≈0 · 7/7 정밀도 0.03~0.27", "원고 수정 (10/30)"),
        ("10", "게이트 M8 '샘플 간 SD 비 ≥ 0.8×' 검토", "생성 SD 0.61~0.73× → 가짜 DEG·과도한 유의성", "제안 (10/30)"),
        ("11", "희귀 subset 판정은 kNN 라벨 유지, marker 재annotation은 보충 근거로", "real B fine 정확도 0.31~0.35 · floor 오차 악화", "결정 필요 (10/30)")]
table(s, rows, 0.9, 1.62, [0.45, 5.35, 4.55, 1.95], h_row=0.405, bold_cols=(0,), font=9.5, head_font=10, hl_rows=(9, 10, 11))
notes(s, """[2:20] 제안
"1번부터 8번은 v3와 같습니다. 6번에 DC 원인 서술을 더했습니다."
"새로 들어간 건 9, 10, 11번입니다."
"9번, DEG 개수나 FDR을 real과 비교할 때는 매칭 n으로 보고하고 하향 DEG를 따로 보고합니다. 지금 원고처럼 전체 n으로 비교하면 floor가 건강 대조 36샘플을 가진 덕에 유의성이 부풀려집니다."
"10번은 새 게이트 제안입니다. 생성 샘플 간 표준편차가 real의 0.8배 이상이어야 한다는 조건입니다. 지금은 0.6에서 0.7배라서 통과하지 못합니다. 이게 가짜 DEG의 직접 원인이라 게이트로 둘 가치가 있습니다."
"11번은 오늘 정할 것입니다. 희귀 subset 판정은 kNN 라벨을 유지하고, marker 기반 재annotation은 reviewer 방어용 보충 자료로 두자는 제안입니다." """)
'''

LIM = r'''# ========================= 23 limitations
s = new(); header(s, "한계", "이 분석이 말할 수 없는 것", NEXT())
side_cards(s, 0.89, 1.66, 5.65, [
    ("paired 조건은 확인 못 함", "paired 생성 카운트가 server_1/spark-0110에만 있음. scSemiProfiler가 대리 증거지만 모델이 다르다."),
    ("draw0 단일", "10/26 보고의 adapter 변동 성분이 오늘 수치에 얼마나 섞였는지 측정되지 않았다."),
    ("패혈증 성인 HC 7샘플", "상한선 간격과 매칭 n 설계의 한계. 7/7 설계는 양쪽 모두 불안정 — bulk(43명) 결과는 이 제약을 받지 않는다."),
    ("환자 구성이 다르다 (10/30)", "패혈증 floor는 real과 환자 중복 0, COVID는 76명 중 32명이 다른 방문으로 중복. '환자 고유 영역'은 일부 구성 차이다.")],
    h=1.17, gap=0.1, accent=[RED, RED, ORANGE, ORANGE], body=11, head=13.5)
side_cards(s, 6.78, 1.66, 5.65, [
    ("read-out 에뮬레이터는 근사", "kNN-15 평활 + 깊이 + 규칙. 실제 μ̂·p̂가 아니므로 사다리 4·5번은 상한 추정."),
    ("음성 대조 (패혈증 vs ICU)", "floor recall 0이 4/5 세포형. 단 real에서도 DEG 7~36개로 적다."),
    ("marker 자의 해상도 (10/30)", "Top10 marker는 real B에서도 fine 정확도 0.31~0.35. 재annotation 결과는 구분 가능한 subset(9개·14개)에 한정해 해석."),
    ("인과의 마지막 한 칸", "조건부 설계인지 LoRA 적응인지는 서버 체크포인트의 μ̂·p̂ 저장 없이 구분 불가.")],
    h=1.17, gap=0.1, accent=[GREY, GREY, GREY, INK], body=11, head=13.5)
notes(s, """[1:50] 한계
"여덟 가지로 늘었습니다. 새로 들어간 건 두 가지입니다."
"왼쪽 넷째, 환자 구성 차이입니다. 패혈증 floor는 real과 환자가 하나도 겹치지 않습니다. 그래서 단핵구 UMAP에서 본 '환자 고유 영역'에는 생성 실패뿐 아니라 환자 구성 차이도 섞여 있습니다. COVID는 반대로 32명이 다른 방문으로 양쪽에 있습니다."
"오른쪽 셋째, marker 자의 해상도입니다. Top10 marker로는 real에서도 fine subset의 3분의 1만 재현되므로, marker 재annotation 결과는 구분 가능한 subset에서만 해석해야 합니다."
"그리고 셋째 카드에 매칭 n의 한계를 붙였습니다. 패혈증 대조군이 7명이라 7 대 7 설계는 real 쪽도 불안정합니다."
"나머지는 v3와 같습니다." """)
'''

DISC = r'''# ========================= 24 discussion
s = new(); header(s, "논의", "오늘 결정하고 싶은 것", NEXT())
three_cards(s, 1.9, 4.3, [
 ("Q1", "수정안 우선순위", "support(B1)에서 count·NB head(A1/C3)로 순서를 바꿀 것인가?\n\n"
   "추가 근거: 생성 샘플 간 SD 0.61~0.73× — 분산 보존도 수정 목표에 넣을지"),
 ("Q2", "게이트 확정", "M7(모듈 magnitude ≥ 0.5×) + M8(샘플 간 SD 비 ≥ 0.8×)을 데이터 주도 모듈 집합으로 고정할 것인가?\n\n"
   "지표 산출 완료, 임계값만 결정"),
 ("Q3", "보고 방식과 서버 요청", "DEG 비교는 매칭 n, 희귀 subset은 kNN 라벨 + marker 보충으로 확정할 것인가?\n\n"
   "서버 요청: ① 다른 draw ② μ̂·p̂ 저장 ③ paired 생성 카운트")],
 accents=[RED, BLUE, ORANGE])
text(s, 0.89, 6.35, 11.55, 0.42, [[("그 외  ", {"bold": True, "color": INK}),
     ("sign agreement 분모 · COVID CI · IEG 제외값 · DC 실패 서술 · COVID 환자 중복 명시 — 모두 원고 수정만 남음", {})]], size=13, color=BODY)
notes(s, """[2:00] 논의
"세 가지입니다."
"Q1, 수정안 순서. v3와 같은 질문인데 근거가 하나 늘었습니다. 생성 샘플 간 분산이 real의 0.6에서 0.7배입니다. 카운트 분포 헤드로 바꿀 때 샘플 간 분산 보존도 목표에 넣을지 같이 정하면 좋겠습니다."
"Q2, 게이트. M7 모듈 magnitude에 M8 샘플 간 표준편차 비를 더하자는 제안입니다. 둘 다 지표는 산출됐고 임계값만 정하면 됩니다."
"Q3, 보고 방식. DEG 비교는 매칭 n으로, 희귀 subset 판정은 kNN 라벨을 기본으로 하고 marker 재annotation은 보충으로 두는 것으로 확정할지입니다. 서버 요청 세 가지는 v3와 같습니다."
"질문 받겠습니다.""")
'''

APPX = r'''# ========================= A appendix divider
s = new()
box(s, 0, 0, 13.333, 7.5, CARD, shape=MSO_SHAPE.RECTANGLE); box(s, 0, 0, 0.17, 7.5, GREY, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 2.3, 11.9, 0.3, "APPENDIX", size=14, color=GREY, bold=True, spc=150)
title_box(s, "부록 — 질의응답용 보충 그림", y=2.9, size=36)
text(s, 0.89, 4.0, 11.9, 1.0, "A1 COVID 림프구 볼케이노 · A2 패혈증 매칭 n 볼케이노 · A3 COVID 매칭 n 볼케이노 · A4 단핵구 subset 비율", size=18, color=BODY, line=1.2)
text(s, 0.89, 6.88, 11.9, 0.26, f"{FOOT} · {NEXT()}", size=12, color=MUTED)
notes(s, """[부록] 질문이 나올 때만 사용합니다.""")

# ========================= A1 covid lymphoid volcano
s = new(); header(s, "부록 A1", "COVID 림프구 볼케이노 — ISG 보존, IFI27 압축, B 과다 호출", NEXT(),
       "COVID/ESKD positive vs negative (환자 단위) · 위 real, 아래 floor")
pic(s, "volcano_lymphoid_covid", 0.75, 1.62, w=8.75)
side_cards(s, 9.75, 1.7, 2.85, [
    ("DEG 수", "T 35 → 33\nB 49 → 106\nNK 59 → 32"),
    ("보존", "MX1 · IFI44L · IFIT1 · RSAD2\n방향·순위 보존"),
    ("손실·누출", "IFI27 압축 · 하향 DEG(KRT86 등) 소실\nS100A8/A9가 T·NK에서 DEG")], h=1.6, accent=[GREY, BLUE, RED], body=11.5, head=14)
notes(s, """[부록 A1] COVID 림프구
"ISG는 방향과 순위가 보존됩니다. IFI27은 T세포에서 2.0이 1.15로 눌립니다. B세포는 logFC 1에서 1.3 사이의 약한 상향 DEG가 106개로 real 49개의 두 배입니다. 하향 DEG는 거의 없습니다. S100A8/A9는 v3 결과 4와 같이 T와 NK에서 DEG로 호출됩니다." """)

# ========================= A2 sepsis matched volcano
s = new(); header(s, "부록 A2", "패혈증 매칭 n 볼케이노 (68 case / 7 control, 50회 중 대표 1회)", NEXT(), "위 real, 아래 floor · y축 공유")
pic(s, "volcano_lymphoid_sepsis_matchedN", 0.75, 1.55, w=8.6)
side_cards(s, 9.6, 1.7, 3.0, [
    ("유의성", "n을 맞추면 floor −log10 FDR이\nreal 범위(≤ 18)로 내려온다"),
    ("up / down", "T up 61 vs 56 · down 45 vs 23\nNK down 54 vs 33"),
    ("남는 불일치", "B: TNF·LILRB4 약화\nJCHAIN floor에서만 DEG")], h=1.6, accent=[GREY, RED, ORANGE], body=11.5, head=14)
notes(s, """[부록 A2] 패혈증 매칭 n
"n을 맞추면 floor의 과도한 유의성이 real 범위로 내려옵니다. 상향 DEG 수는 같고 하향 DEG가 절반 수준입니다." """)

# ========================= A3 covid matched volcano
s = new(); header(s, "부록 A3", "COVID 매칭 n 볼케이노 (환자 37~39 / 18, 대표 1회)", NEXT(), "위 real, 아래 floor")
pic(s, "volcano_lymphoid_covid_matchedN", 0.75, 1.55, w=8.6)
side_cards(s, 9.6, 1.7, 3.0, [
    ("logFC r", "T 0.57 · B 0.44 · NK 0.39"),
    ("down-DEG", "real 2.9~15.6 → floor 0.04~1.1\n사실상 소실"),
    ("B 과다 호출", "up 31 → 95, 정밀도 0.17")], h=1.6, accent=[GREY, RED, ORANGE], body=11.5, head=14)
notes(s, """[부록 A3] COVID 매칭 n
"COVID에서 하향 DEG는 사실상 0이고, B세포 상향 DEG가 세 배로 부풀려집니다." """)

# ========================= A4 myeloid composition
s = new(); header(s, "부록 A4", "단핵구 구획 내 subset 비율 — real A vs floor (kNN 라벨)", NEXT(), "라벨 = floor/real 비율, 빨강 < 0.5×")
pic(s, "composition_myeloid", 0.9, 1.7, w=11.5)
text(s, 0.89, 5.75, 11.55, 1.0, [[("패턴  ", {"bold": True, "color": INK}),
     ("희귀 상태 결손(MKI67 0.02×, dex 0.03×, Macrophage 0.25×, CD16mono_C1 0.38×)과 흔한 상태 과생산(RETN 1.36×, CD16 1.45×, CD16mono 1.47×)이 짝을 이룬다.", {})]],
     size=12.5, color=BODY, line=1.15)
notes(s, """[부록 A4] 단핵구 subset 비율
"희귀 상태가 빠진 만큼 흔한 상태가 과생산됩니다. 전체 비율의 합은 1이므로 당연한 면도 있지만, 어느 쪽이 빠지는지가 일관됩니다." """)
'''

# ======================= v5 additions (2026-10-31) =======================
SEMI1 = r'''# ========================= 20g scSemi scorecard (NEW 10/31)
s = new(); header(s, "비교 모델 · 10/31", "scSemiProfiler와 같은 지표로 비교 — floor가 앞선다", NEXT(),
       "같은 파이프라인 · scSemi = strict_5pairs draw0, 정답 real 샘플을 재생성")
''' + NEW_TAG2 + r'''rows = [("지표 (패혈증 / COVID)", "PBISC floor", "scSemiProfiler", "기준 · 비고"),
        ("질병 logFC r — 패혈증 T·B·NK·Mono", "0.68~0.72", "0.20~0.29", "real-real 상한 0.59~0.97"),
        ("질병 logFC r — COVID", "0.47~0.66", "0.16~0.33", "상한 0.65~0.83"),
        ("상한 대비 달성률 (보수적 하한)", "60~83% (패혈증 DC 31%)", "21~47% (패혈증 DC 6%)", "같은 샘플 paired 설계도 동일"),
        ("DEG recall / precision (중앙값)", "0.44/0.50 · 0.34/0.58", "0.06/0.20 · 0.35/0.13", "COVID scSemi: DEG 과다 호출"),
        ("하향 DEG 보존 (5개 세포형 합)", "226/921 · 5/94", "179/921 · 34/94", "둘 다 하향 손실"),
        ("샘플 간 SD / real (중앙값)", "0.79 · 0.64 (너무 비슷)", "3.09 · 2.49 (너무 다름)", "방향이 반대"),
        ("모듈: C1Q · DC 정체성", "0.10 · 0.09", "0.39 · 0.26", "scSemi = 실제 대표 세포 복사"),
        ("모듈: 골수계 활성화 · ISG", "1.03 · 0.80", "0.24 · 2.23", "scSemi 결손 또는 과대"),
        ("희귀 subset < 0.5× (transfer-free, 27개)", "11/27", "18/27", "real B 6/27"),
        ("marker 재annotation 비율 오차", "1.53 · 1.10", "3.42 · 2.71", "구분 가능 subset만"),
        ("전용성 법칙 ρ (유전자 · 모듈)", "−0.24 · −0.46", "+0.11 · +0.10 (n.s.)", "scSemi는 기울기 없음"),
        ("세포당 검출 유전자 (중앙값)", "468 · 414", "98 · 131", "real B 412 · 340")]
table(s, rows, 0.9, 1.62, [3.6, 2.55, 2.55, 3.15], h_row=0.345, bold_cols=(0,), font=10.5, head_font=11, hl_rows=(3, 6, 11))
notes(s, """[2:10] scSemiProfiler 동일 지표 비교 (추가)
"질문은 '위와 같은 수치들이 scSemiProfiler에서는 어떻게 나오느냐'였습니다. floor에 쓴 것과 똑같은 pseudobulk, DEG, 희귀 집단, marker 파이프라인을 그대로 돌렸습니다."
"먼저 조건 차이를 말씀드려야 합니다. scSemi는 정답이 측정된 real 샘플을 다시 생성합니다. floor는 scRNA가 없는 bulk-only 샘플에서 생성합니다. scSemi 쪽이 훨씬 쉬운 조건입니다. 그래서 real을 같은 샘플로 제한한 paired 설계도 같이 돌렸는데 결론은 같았습니다."
"위 세 줄이 질병 DEG입니다. 패혈증 logFC r이 floor 0.68에서 0.72, scSemi 0.20에서 0.29입니다. 상한 대비로는 floor 60에서 83퍼센트, scSemi 21에서 47퍼센트입니다."
"파랗게 칠한 샘플 간 SD 줄이 흥미롭습니다. floor는 real의 0.6에서 0.8배로 샘플끼리 너무 비슷하고, scSemi는 2.5에서 3배로 너무 다릅니다. 방향이 반대입니다."
"모듈에서 scSemi가 floor보다 나은 건 C1Q와 DC 정체성 두 개뿐입니다. scSemi가 실제 대표 세포를 복사하기 때문입니다. 대신 골수계 활성화는 0.24로 무너지고 ISG는 2.2배로 과대입니다."
"희귀 subset은 27개 중 18개가 0.5배 미만입니다. floor는 11개, real B는 6개입니다."
"마지막으로 전용성 법칙입니다. floor는 유전자 마이너스 0.24, 모듈 마이너스 0.46인데 scSemi는 기울기가 없습니다. 다음 두 장에서 그림으로 보겠습니다." """)
'''

SEMI2 = r'''# ========================= 20h scSemi DEG + modules (NEW 10/31)
s = new(); header(s, "비교 모델 · 10/31", "질병 DEG — floor는 상한의 2/3, scSemi는 1/4", NEXT(),
       "a 선 = real-real 상한 두 정의 · b 모듈 회복 (log 축, 1 = 완벽, 점선 0.5)")
''' + NEW_TAG2 + r'''pic(s, "v5_scsemi_ab", 0.85, 1.62, w=11.65)
text(s, 0.89, 5.45, 11.55, 1.35, [
  [("a  ", {"bold": True, "color": INK}), ("10개 코호트×세포형 모두 floor > scSemi. scSemi는 패혈증 NK DEG 90 → 1, COVID T·NK는 34 → 210, 55 → 304로 과다 호출 (precision 0.09~0.11).", {})],
  [("b  ", {"bold": True, "color": INK}), ("두 모델의 실패 모듈이 다르다: floor는 세포형 전용 모듈(C1Q·DC)에서, scSemi는 전신 모듈(골수계 활성화 0.24)과 과대(MHC-II 3.07×, ISG 2.23×)에서 틀린다.", {})]],
  size=12, color=BODY, line=1.1, after=4)
notes(s, """[1:50] scSemi DEG와 모듈 (추가)
"a는 질병 logFC 상관입니다. 검은 선과 회색 선이 real-real 상한의 두 정의입니다. 파랑 floor는 대부분 회색 선 근처에 있고, 주황 scSemi는 0.2에서 0.3에 머뭅니다. 패혈증 DC는 0.04로 사실상 0입니다."
"DEG 개수도 이상합니다. 패혈증 NK에서는 real 90개 중 1개만 부르고, COVID T와 NK에서는 real의 6배를 부릅니다. 정밀도가 0.1 수준이니 대부분 가짜입니다."
"b가 개념적으로 중요합니다. floor는 C1Q와 DC 정체성 같은 세포형 전용 모듈에서 무너지고, 골수계 활성화나 ISG 같은 전신 모듈은 잘 맞습니다. scSemi는 거꾸로 전용 모듈은 상대적으로 낫고 전신 모듈에서 무너지거나 과대입니다."
"즉 두 모델은 서로 다른 곳에서 틀립니다. floor에서 본 패턴이 평가 파이프라인이 만든 인공물이 아니라는 근거가 됩니다." """)
'''

SEMI3 = r'''# ========================= 20i scSemi variance + rare (NEW 10/31)
s = new(); header(s, "비교 모델 · 10/31", "샘플 간 변동과 희귀 집단 — 실패 방향이 반대다", NEXT(),
       "c scSemi는 같은 real 샘플과 비교 · d 라벨 전이 없는 추정")
''' + NEW_TAG2 + r'''pic(s, "v5_scsemi_cd", 0.85, 1.62, w=11.65)
text(s, 0.89, 5.3, 11.55, 1.5, [
  [("c  ", {"bold": True, "color": INK}), ("floor 0.6~0.8× (생성 샘플이 서로 너무 비슷 → 유의성 과장) vs scSemi 2.0~4.6× (너무 다름 → 잡음). scSemi pseudobulk는 샘플링 전 rate의 합이라 Poisson 잡음 없이도 이만큼 크다.", {})],
  [("d  ", {"bold": True, "color": INK}), ("scSemi는 Macrophage 0.06×, CD14mono_dex 0, Mono_c06-MKI67 0.15×. 생성된 희귀 세포의 marker 점수는 오히려 1.1 (실제 세포 복사) — ", {}),
   ("'적게 만들되 만든 것은 진짜'", {"bold": True, "color": RED}), (". floor는 '대체로 만들되 정체성이 희석'.", {})]],
  size=12, color=BODY, line=1.1, after=4)
notes(s, """[1:50] 샘플 간 변동과 희귀 집단 (추가)
"c는 샘플 간 표준편차 비입니다. floor는 1 아래, scSemi는 1 위 2에서 4.6배입니다. floor는 샘플끼리 너무 비슷해서 유의성이 과장되고, scSemi는 너무 달라서 잡음이 큽니다. scSemi pseudobulk는 Poisson 추출 전의 rate를 합한 것이라, 샘플링 잡음 없이도 이 정도라는 점이 중요합니다."
"d는 희귀 subset 27개를 라벨 전이 없이 추정한 비율입니다. kNN 라벨을 쓰지 않은 이유는 scSemi가 real과 통합되지 않기 때문입니다. Harmony 이후 iLISI가 1.00입니다."
"real B는 대부분 0.5에서 2배 띠 안에 있고, floor는 11개, scSemi는 18개가 띠 아래입니다. 대식세포 0.06배, dex 단핵구는 0입니다."
"그런데 scSemi가 만든 희귀 세포의 marker 점수는 1.1로 오히려 real보다 높습니다. 실제 대표 세포를 복사하기 때문입니다. 한 문장으로 하면 scSemi는 '적게 만들되 만든 것은 진짜'이고, floor는 '대체로 만들되 정체성이 희석된다'입니다." """)
'''

SEMI4 = r'''# ========================= 20j scSemi interpretation (NEW 10/31)
s = new(); header(s, "비교 모델 · 10/31", "두 모델은 다르게 실패한다 — floor의 결함은 선택적이다", NEXT())
''' + NEW_TAG2 + r'''for i, (hh, col, items) in enumerate([
    ("PBISC floor — 선택적 실패", BLUE, ["전신 프로그램은 보존 (골수계 활성화 1.03×, ISG 0.80×)",
                                         "세포형 전용 모듈의 magnitude 손실 (C1Q 0.10×, DC 0.09×)",
                                         "전용성이 높을수록 회복률 낮음 (ρ −0.24 / −0.46)",
                                         "샘플 간 변동 압축 (0.6~0.8×) → 하향 DEG 소실·유의성 과장",
                                         "희귀 상태는 대체로 생성되나 정체성 희석 (marker 0.86)"]),
    ("scSemiProfiler — 전역적 실패", ORANGE, ["세포당 ~100 유전자, 유전자×세포형 51%가 0",
                                               "질병 방향 r 0.05~0.33, DEG 과다·과소 호출",
                                               "전용성 기울기 없음 (ρ +0.11 / +0.10) — 전부 무너짐",
                                               "샘플 간 변동 과대 (2.0~4.6×)",
                                               "희귀 상태 결손 18/27, 생성된 것은 실제 세포 복사 (marker 1.1)"])]):
    x = 0.89 + i * 5.9; card(s, x, 1.72, 5.65, 2.85, col)
    text(s, x + 0.25, 1.92, 5.2, 0.4, hh, size=17, color=INK, bold=True)
    text(s, x + 0.25, 2.45, 5.2, 2.0, [[("•  " + t, {})] for t in items], size=12.5, color=BODY, line=1.12, after=6)
text(s, 0.89, 4.85, 11.55, 1.9, [
  [("원고에서의 활용  ", {"bold": True, "color": INK}),
   ("(1) 비교표: 정답 샘플을 재생성하는 더 쉬운 조건의 scSemi보다 bulk-only floor가 DEG·희귀 지표 전부에서 앞선다. "
    "(2) 전용성 법칙과 분산 압축은 평가 파이프라인의 인공물이 아니다 — 같은 파이프라인에서 scSemi는 다른 양상을 보인다. "
    "(3) 단, C1Q·DC 모듈은 scSemi가 더 낫다는 점을 숨기지 않는다 (대표 세포 복사의 이점).", {})]],
  size=13.5, color=BODY, line=1.2)
notes(s, """[1:40] 해석 (추가)
"정리하면 두 모델은 다르게 실패합니다."
"floor는 선택적입니다. 전신 프로그램은 옮기고, 세포형 전용 모듈의 크기를 잃고, 그 정도가 전용성에 비례하며, 샘플 간 변동이 압축됩니다."
"scSemi는 전역적입니다. 세포당 유전자 100개 정도로 희소해서 거의 모든 것이 무너지고, 그래서 전용성 기울기도 없습니다. 샘플 간 변동은 오히려 과대입니다."
"원고에서는 세 가지로 쓸 수 있습니다. 첫째, 비교표입니다. scSemi는 정답 샘플을 재생성하는 쉬운 조건인데도 floor가 모든 DEG와 희귀 지표에서 앞섭니다. 둘째, 우리가 찾은 전용성 법칙과 분산 압축은 평가 파이프라인이 만든 인공물이 아닙니다. 같은 파이프라인에서 scSemi는 다른 양상을 보입니다. 셋째, C1Q와 DC 모듈은 scSemi가 낫다는 점을 숨기지 말고 대표 세포 복사의 이점으로 설명하는 게 좋겠습니다." """)
'''

SYN = SYN.replace('''        ("라벨링 (10/30)", "kNN 전이 탓인가", "marker 재annotation 시 비율 오차 0.35→1.53", "원인 아님")]
table(s, rows, 0.9, 1.68, [2.4, 3.3, 4.2, 2.0], h_row=0.39, bold_cols=(0,), font=11, head_font=11, warn_rows=(5, 6))
text(s, 0.89, 5.45, 11.55, 1.35, [''', '''        ("라벨링 (10/30)", "kNN 전이 탓인가", "marker 재annotation 시 비율 오차 0.35→1.53", "원인 아님"),
        ("비교 모델 (10/31)", "평가 파이프라인의 인공물인가", "scSemi: r 0.05~0.33, SD 2.5~3.1×, 기울기 없음", "아님 — floor 고유 양상")]
table(s, rows, 0.9, 1.66, [2.4, 3.3, 4.2, 2.0], h_row=0.37, bold_cols=(0,), font=10.5, head_font=11, warn_rows=(5, 6))
text(s, 0.89, 5.62, 11.55, 1.2, [''')
SYN = SYN.replace('"라벨링도 원인이 아닙니다. marker로 다시 붙이면 오히려 나빠집니다."',
                  '"라벨링도 원인이 아닙니다. marker로 다시 붙이면 오히려 나빠집니다."\n"마지막 줄은 10/31 추가입니다. 같은 파이프라인에서 scSemiProfiler는 전혀 다른 양상을 보이므로, 우리가 본 패턴은 평가 방법이 만든 인공물이 아닙니다."')

ACT = ACT.replace('''        ("11", "희귀 subset 판정은 kNN 라벨 유지, marker 재annotation은 보충 근거로", "real B fine 정확도 0.31~0.35 · floor 오차 악화", "결정 필요 (10/30)")]
table(s, rows, 0.9, 1.62, [0.45, 5.35, 4.55, 1.95], h_row=0.405, bold_cols=(0,), font=9.5, head_font=10, hl_rows=(9, 10, 11))''',
'''        ("11", "희귀 subset 판정은 kNN 라벨 유지, marker 재annotation은 보충 근거로", "real B fine 정확도 0.31~0.35 · floor 오차 악화", "결정 필요 (10/30)"),
        ("12", "scSemiProfiler 동일 지표 비교표를 Results/Supp.에 추가 (paired 설계 포함)", "DEG r 0.68~0.72 vs 0.20~0.29 · 희귀 결손 11 vs 18/27", "원고 추가 (10/31)")]
table(s, rows, 0.9, 1.6, [0.45, 5.35, 4.55, 1.95], h_row=0.38, bold_cols=(0,), font=9.5, head_font=10, hl_rows=(9, 10, 11, 12))''')
ACT = ACT.replace('"11번은 오늘 정할 것입니다.','"12번은 10/31 추가입니다. scSemiProfiler를 같은 지표로 평가한 비교표를 본문이나 보충 자료에 넣자는 제안입니다. 같은 샘플 paired 설계까지 포함하면 리뷰어가 조건 차이를 지적해도 방어됩니다."\n"11번은 오늘 정할 것입니다.')

LIM = LIM.replace('("음성 대조 (패혈증 vs ICU)", "floor recall 0이 4/5 세포형. 단 real에서도 DEG 7~36개로 적다."),',
                  '("음성 대조 · scSemi (10/31)", "음성 대조는 real DEG도 7~36개로 적다. scSemi는 real과 통합 실패(iLISI 1.00) → 희귀 결론은 transfer-free 추정, draw0·대조군 6샘플."),')
LIM = LIM.replace('"그리고 셋째 카드에 매칭 n의 한계를 붙였습니다.','"오른쪽 둘째 카드에 scSemi 비교의 한계를 붙였습니다. scSemi는 real과 통합이 안 되므로 kNN 라벨을 믿을 수 없어 희귀 집단은 라벨 전이 없는 추정으로 판단했고, draw0 하나에 패혈증 대조군 6샘플입니다."\n"그리고 셋째 카드에 매칭 n의 한계를 붙였습니다.')

DISC = DISC.replace('("sign agreement 분모 · COVID CI · IEG 제외값 · DC 실패 서술 · COVID 환자 중복 명시 — 모두 원고 수정만 남음", {})]], size=13, color=BODY)',
                    '("sign agreement 분모 · COVID CI · IEG 제외값 · DC 실패 서술 · COVID 환자 중복 · scSemi 비교표 — 원고 수정만 남음", {})]], size=12.5, color=BODY)')
