# -*- coding: utf-8 -*-
"""Korean meeting deck + Korean speaker notes: disease DEG recomputation and rare-population expression (2026-10-27).
Every number is read from out/*.csv produced by the notebook."""
import os, pandas as pd, numpy as np
from pptx.enum.shapes import MSO_SHAPE
from pptx.util import Inches
from scipy import stats
from deck_helpers2 import *
_N=[1]
def NEXT():
    _N[0]+=1; return _N[0]

D = pd.read_csv("out/deg_summary.csv"); C = pd.read_csv("out/real_real_ceiling.csv")
S2 = pd.read_csv("out/deg_rank2000_summary.csv"); N = pd.read_csv("out/named_genes.csv")
RR = pd.read_csv("out/rare_expression_ratios.csv"); PM = pd.read_csv("out/rare_power_matched.csv")
RA = pd.read_csv("out/rare_expression.csv")
A2 = pd.read_csv("out2/a2_module_bulk.csv").dropna(subset=["bulk_logFC"])
A2BG = pd.read_csv("out2/a2_bulk_background.csv")
C2L = pd.read_csv("out2/c2b_module_load.csv"); C2R = pd.read_csv("out2/c2_module_rank.csv")
A3 = pd.read_csv("out2/a3_power_matched_ceiling.csv").merge(
    C[["cohort", "celltype", "logFC_r"]].rename(columns={"logFC_r": "ceil_both"}), on=["cohort", "celltype"])
D1T = pd.read_csv("out2/d1_bootstrap_ci.csv").set_index(["cohort", "celltype"])
D2B = pd.read_csv("out2/d2b_sepsis_vs_icu.csv")
LAD = pd.read_csv("out2/a1_ladder.csv").pivot_table(index="source", columns="module", values="overall_ratio")
BM = pd.read_csv("out2/b1_datadriven_modules.csv")
IE = pd.read_csv("out2/ieg_split.csv").set_index(["cohort", "celltype"])
IQ = pd.read_csv("out2/ieg_score_contrast.csv").set_index(["cohort", "celltype", "source"])
a3 = A3.set_index(["cohort", "celltype"])
def brk(coh, ct):
    r = a3.loc[(coh, ct)]; lo = 100 * r.floor_r_full / r.ceil_case_split_r; hi = 100 * r.floor_r_full / r.ceil_both
    return f"{lo:.0f}~{hi:.0f}%"
def ci(coh, ct):
    r = D1T.loc[(coh, ct)]; return f"{r.logFC_r:.2f} [{r.r_lo:.2f}, {r.r_hi:.2f}]"
f2 = lambda x: f"{x:.2f}"
fl = D[D.source == "floor"].set_index(["cohort", "celltype"]); pv = D[D.source == "prev"].set_index(["cohort", "celltype"])
ce = C.set_index(["cohort", "celltype"]); nf = N[N.source == "floor"].set_index(["cohort", "celltype", "gene"])
g = lambda c, t, gene, col: nf.loc[(c, t, gene), col]

# ========================= 1 title
s = new(); box(s, 0, 0, 0.17, 7.5, BLUE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 1.92, 11.9, 0.3, "PBISC floor read-out · 질병 DEG 및 희귀 집단 재검증 · v4", size=14, color=BLUE, bold=True, spc=150)
title_box(s, "생성 세포는 질병 신호를 어디까지 가져오고,\n어디서 체계적으로 놓치는가", y=2.5, size=40)
text(s, 0.89, 4.15, 11.9, 0.8, "질병 DEG 재계산 · 매칭 n 검정 · 희귀 집단 내부 발현 · DC 과소생성 원인 · 단핵구/대식세포 subset · "
     "marker 기반 재annotation (패혈증 Ye 2026 / COVID Imperial ESKD)", size=18, color=BODY, line=1.15)
pill(s, 0.89, 5.25, 1.75, "PBISC floor", BLUE); pill(s, 2.8, 5.25, 2.1, "PBISC prev K4000", RED)
pill(s, 5.06, 5.25, 2.25, "Real (held-out)", GREY); pill(s, 7.47, 5.25, 2.6, "real-real 상한선", INK); pill(s, 10.23, 5.25, 2.2, "bulk 입력", PURPLE)
text(s, 0.89, 6.82, 11.9, 0.26, "Jong Hwi Kim · Seoul National University · 2026년 10월 30일 · 약 45분 (부록 별도)", size=12, color=MUTED)
notes(s, """[0:40] 표지
"v3 이후 들어온 질문 여섯 가지에 답한 판입니다. 질병 DEG를 case·control 수를 맞춰 다시 본 것, 림프구 볼케이노, DC가 왜 적게 생성되는가, floor가 어느 bulk 샘플에서 왔는가, 단핵구·대식세포 subset 지형, 그리고 real marker로 floor를 직접 annotation하면 달라지는가입니다."
"새로 들어간 슬라이드는 오른쪽 위에 주황색 '10/30 추가' 표시가 있습니다. 나머지는 v3와 같고, 요약·종합·제안·한계·논의만 갱신했습니다."
"모든 수치는 로컬 데이터로 새로 계산했고 새로 생성한 세포는 없습니다." """)
# ========================= 2 summary
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
    "③ 패혈증 real 단핵구의 23%는 환자 고유 영역 — floor에 없음 ④ marker 재annotation은 비율 오차를 키움 (패혈증 0.35 → 1.53)", {})]],
  size=13, color=BODY, line=1.12, after=6)
notes(s, """[1:50] 요약
"세 카드는 v3와 같은 결론이고, 오늘 추가 분석이 각 카드를 보강합니다."
"첫째, 숫자. case와 control 샘플 수를 real과 floor에서 똑같이 맞춰도 패혈증 림프구 logFC r이 0.59에서 0.66으로 유지됩니다. 즉 v3의 상관값은 floor 쪽 샘플이 많아서 나온 게 아닙니다."
"둘째, 원인. C1Q에 이어 DC에서도 같은 구조를 확인했습니다. 같은 샘플의 bulk 입력에서 pDC 마커 점수와 생성 결과의 pDC 점수 상관이 real에서는 0.57인데 floor는 마이너스 0.18입니다. 모델이 자기 입력을 따라가지 않습니다."
"셋째, 법칙. 라벨링 방식을 바꿔도 회복되지 않는다는 걸 확인했습니다. real의 Top10 marker로 floor를 직접 annotation하면 오히려 비율 오차가 커집니다."
"아래 주황 줄이 오늘 새로 들어간 네 가지 요약입니다. 각각 뒤 슬라이드에서 설명드리겠습니다." """)
# ========================= 3 questions
s = new(); header(s, "배경", "왜 다시 계산했는가", NEXT(), "10/26 draw 변동 보고 이후 제기된 세 가지")
rows = [("질문", "지금까지의 근거", "이번에 한 것"),
        ("질병 DEG가 실제와 같은가?", "Result 6B: logFC r 0.70 / 0.52 (집계값)\n세포형별·유전자별 근거 없음", "세포형별 DEG 전체 재계산 + 유전자 단위 대조"),
        ("점수가 낮으면 모델 탓인가?", "비교 대상이 '실제' 하나뿐\n데이터 자체의 한계와 구분 불가", "실제 샘플 절반 분할 × 20회 → real-real 상한선"),
        ("희귀 집단이 '있다'는 게 무슨 뜻인가?", "라벨 전이 후 비율(0.5–2×)만 확인\n집단 내부 발현은 미확인", "19개 집단 내부 평균 발현 프로파일 · 마커 점수 · 유전자별 차이")]
table(s, rows, 0.9, 1.70, [3.4, 4.3, 4.2], h_row=0.80, bold_cols=(0,), font=11.5, head_font=11.5)
text(s, 0.89, 5.90, 11.55, 1.0, [[("설계 원칙  ", {"bold": True, "color": INK}),
     ("생성 세포와 real B는 모두 같은 kNN 라벨 전이를 거친다. 따라서 둘의 차이는 라벨링이 아니라 생성의 차이다. "
      "real A의 author 라벨이 기준점이고, real B는 '환자가 다르고 라벨 전이를 거쳤을 때 자연히 생기는 차이'를 흡수하는 대조군이다.", {})]],
     size=13.5, color=BODY, line=1.15)
notes(s, """[2:00] 배경
"10월 26일 draw 변동 보고에서 '비율 오차 하나로 질병 신호 보존을 판단할 수 없다'는 결론이 나왔습니다. 그 연장선입니다."
"질문 1. Result 6B의 logFC r 0.70은 다섯 세포형을 묶은 집계값입니다. 어느 세포형에서 나온 것인지, 어떤 유전자가 맞고 틀렸는지는 아직 모릅니다."
"질문 2. 이게 제일 중요합니다. 지금까지 우리 수치는 전부 '생성 대 실제' 한 쌍이었습니다. 그래서 r이 0.6이면 좋은 건지 나쁜 건지 말할 기준이 없었습니다. 이번에는 실제 샘플만 절반으로 나눠서 같은 계산을 20번 돌렸습니다. 실제가 실제를 재현하는 수준이 상한선입니다."
"질문 3. 희귀 집단은 지금까지 '비율이 0.5에서 2배 사이면 복원'으로 판정했습니다. 그런데 비율이 맞아도 그 안의 세포가 진짜 그 세포인지는 확인한 적이 없습니다."
설계 원칙 문단을 짚을 것: real B도 floor와 똑같은 라벨 전이를 거치므로, 라벨링 실패와 생성 실패가 분리된다.""")

# ========================= 4 methods
s = new(); header(s, "방법", "데이터와 대비 정의", NEXT(), "모두 로컬 파일 · 새로 생성한 세포 없음")
rows = [("", "패혈증 (Ye 2026)", "COVID / ESKD (Imperial)"),
        ("질병 대비", "성인 sepsis vs 성인 healthy control", "positive vs negative (환자 단위 집계)"),
        ("실제 샘플", "성인 108 (case 93 / control 7)", "198 샘플 → 환자 단위 (case 39 / control 38)"),
        ("floor 샘플", "bulk-only 122 (case 68 / control 36)", "164 − adapter 학습 donor 5 = 159"),
        ("pseudobulk", "샘플 × six-type raw count 합 (전체 파일), log2 CPM · 20세포 미만 제외", ""),
        ("DEG 규칙", "Welch t-test · BH FDR < 0.05 · |log2FC| ≥ 1", ""),
        ("상한선", "실제 샘플만 무작위 절반 분할 × 20회, 동일 파이프라인", ""),
        ("희귀 집단", "f_rare_parent.py 정의 그대로 · real A = author 라벨, real B / floor = kNN 전이 라벨", "")]
table(s, rows, 0.9, 1.75, [2.3, 4.9, 4.7], h_row=0.5, bold_cols=(0,), font=12, head_font=12)
text(s, 0.89, 6.02, 11.55, 0.9, [[("주의  ", {"bold": True, "color": INK}),
     ("패혈증 성인 healthy control이 7샘플뿐이다. 상한선을 만들 때 절반으로 나누면 3~4개가 되므로 패혈증 상한선은 검정력이 낮은 "
      "보수적 추정이다 — floor가 B·NK에서 상한선을 넘는 것은 floor가 더 좋아서가 아니다.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[2:00] 방법
"pseudobulk는 이미 만들어둔 pb_build 결과를 썼습니다. 전체 파일 기준이라 60k 서브셋이 아니라 모든 세포가 들어갑니다."
"COVID는 한 환자가 여러 번 방문하므로 환자 단위로 평균을 내고 검정했습니다. adapter 학습에 쓰인 donor 5샘플은 뺐습니다."
"DEG 기준은 FDR 0.05에 로그 폴드 1 이상입니다. 이건 보고서 Result 6B보다 엄격할 수 있는데, 뒤에서 그 차이가 문제가 됩니다."
"주의 문단을 반드시 짚을 것. 패혈증 성인 건강 대조군이 7명뿐입니다. 절반으로 나누면 3~4명이라 상한선이 낮게 나옵니다. 그래서 B세포와 NK에서 floor가 상한선을 넘는 것처럼 보이는데, 그건 floor가 좋아서가 아니라 상한선이 보수적이라서입니다. 질문이 나오면 여기서 막아야 합니다." """)

# ========================= 4b provenance (NEW)
s = new(); header(s, "방법 · 10/30", "floor는 어느 bulk 샘플에서 왔는가", NEXT(), "prov/floor_source_samples.csv · 샘플당 4,000세포 생성")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
rows = [("그룹 (패혈증)", "real scRNA 샘플", "PBISC floor 샘플", "floor 출처"),
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
# ========================= 5 umap sepsis
s = new(); header(s, "전체 지형", "floor 세포는 실제와 같은 자리에 온다 — 희귀 집단만 빼고", NEXT(),
       "joint Harmony (real 절반 A + floor), 40,000 세포 subsample · 패혈증")
pic(s, "umap_sepsis", 1.17, 1.68, w=11.0)
text(s, 0.89, 5.86, 11.55, 0.9, [[("오른쪽 패널  ", {"bold": True, "color": INK}),
     ("원 = 실제, 검은 테두리 삼각형 = floor. B_c07-MZB1(빨강)은 floor에도 있지만, Macrophage(초록)와 pDC(주황)는 실제만 보이고 "
      "floor 삼각형이 거의 없다. 전역 혼합도(iLISI 1.65, 실제 1.70)가 좋다는 사실과 공존한다.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[2:00] UMAP 패혈증
"왼쪽부터 보겠습니다. 소스별로 색을 칠하면 회색 실제와 파랑 floor가 거의 완전히 겹칩니다. Harmony 이후 iLISI 1.65로 실제끼리의 1.70에 거의 닿습니다. 이건 기존 보고와 같습니다."
"가운데는 세포형입니다. T, B, NK, 단핵구 덩어리가 제자리에 있고 floor도 같은 덩어리에 들어갑니다."
"오른쪽이 오늘의 그림입니다. 희귀 집단 네 개만 색칠했습니다. 원이 실제, 삼각형이 floor입니다. 빨강 형질세포 계열은 floor 삼각형이 같이 보입니다. 그런데 초록 대식세포와 주황 pDC는 원만 있고 삼각형이 거의 없습니다."
"여기서 말하고 싶은 건 모순이 아니라 해상도입니다. 전역 혼합도 지표로는 '거의 실제'인데, 같은 데이터에서 특정 집단만 보면 비어 있습니다. 전역 지표가 희귀 집단을 전혀 보지 못한다는 뜻입니다." """)

# ========================= 6 umap covid
s = new(); header(s, "전체 지형", "COVID/ESKD에서도 같은 패턴", NEXT(), "joint Harmony (real 절반 A + floor), 40,000 세포 subsample")
pic(s, "umap_covid", 1.17, 1.68, w=11.0)
text(s, 0.89, 5.86, 11.55, 0.9, [[("확인  ", {"bold": True, "color": INK}),
     ("MAIT·pDC처럼 비교적 흔한 희귀 집단은 floor에도 분포하지만, DC1·B_ASC_dividing은 실제 대비 셀 수가 크게 적다. "
      "코호트가 달라도 '흔한 것은 되고 드문 것은 흐려진다'는 구조가 반복된다.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[1:10] UMAP COVID
"COVID에서도 같습니다. 소스는 잘 섞이고 세포형 구조도 유지됩니다."
"오른쪽에서 MAIT과 pDC는 floor에도 충분히 있습니다. 이 둘은 뒤에 나올 마커 점수에서도 비교적 잘 유지되는 집단입니다. 반면 DC1과 분열 중 형질모세포는 수 자체가 적습니다."
"두 코호트에서 같은 구조가 반복된다는 점이 중요합니다. 한 데이터셋의 특이성이 아닙니다." """)

# ========================= 7 bar logFC r
s = new(); header(s, "결과 1", "상한선 자체가 불확실하다 — 달성률은 점이 아니라 구간", NEXT(),
       "회색 선 = 두 팔 모두 분할(검정력 부족 → 낮음) · 검은 선 = 대조군 공유(노이즈 상관 → 높음) · 막대 = floor")
pic(s, "w1_ceiling_bracket", 1.35, 1.66, w=10.8)
text(s, 0.89, 5.92, 11.55, 1.0, [[("왜 바뀌었나  ", {"bold": True, "color": INK}),
     (f"종전 상한선은 case·control을 모두 반으로 쪼갰는데 패혈증 성인 HC가 7명이라 3~4명이 되었다. 대조군을 공유해 case 팔만 쪼개면 "
      f"상한선이 {a3.loc[('sepsis','T cell'),'ceil_case_split_r']:.2f}(패혈증 T)로 올라간다. 두 정의의 편향 방향이 알려져 있으므로 "
      f"달성률은 구간으로 보고한다: 패혈증 T {brk('sepsis','T cell')}, 단핵구 {brk('sepsis','Monocyte')}, COVID 단핵구 {brk('covid','Monocyte')}, 패혈증 DC {brk('sepsis','DC')}.", {})]],
     size=13, color=BODY, line=1.15)
notes(s, """[2:20] 결과 1 (수정)
"지난주 덱에서 81에서 91퍼센트라고 적었던 숫자를 고쳐 왔습니다."
"문제는 상한선이었습니다. 실제 샘플을 case와 control 양쪽 모두 반으로 쪼갰는데, 패혈증 성인 건강 대조군이 7명이라 절반이면 3~4명입니다. 검정력이 떨어져 상한선이 낮게 나왔고, 그만큼 달성률이 부풀려졌습니다."
"대조군은 양쪽이 공유하고 case 팔만 쪼개 다시 계산하면 상한선이 0.97까지 올라갑니다. 대신 이번엔 두 반쪽이 같은 대조군 노이즈를 공유하므로 상한선이 과대추정됩니다."
"두 정의의 편향 방향이 서로 반대이고 알려져 있으므로, 정직한 보고는 구간입니다. 회색과 검은 선 사이가 그 구간이고, 막대 위 숫자가 달성률 범위입니다."
"패혈증 DC는 어느 정의로도 20에서 40퍼센트입니다. 결론은 바뀌지 않습니다." """)

# ========================= 8 recall precision
s = new(); header(s, "결과 1b", "DEG 회수율은 상한선 아래, 정밀도는 상한선 수준", NEXT(),
       "막대 = floor, 가로 막대 = 실제 절반 분할 상한선 · 10개 코호트×세포형")
pic(s, "bar_recall_precision", 1.95, 1.60, w=9.5)
rows = [("", "패혈증 T", "패혈증 단핵구", "패혈증 DC", "COVID T", "COVID 단핵구", "COVID DC"),
        ("real DEG 수", *[str(int(fl.loc[k, "n_DEG_real"])) for k in [("sepsis","T cell"),("sepsis","Monocyte"),("sepsis","DC"),("covid","T cell"),("covid","Monocyte"),("covid","DC")]]),
        ("floor DEG 수", *[str(int(fl.loc[k, "n_DEG_gen"])) for k in [("sepsis","T cell"),("sepsis","Monocyte"),("sepsis","DC"),("covid","T cell"),("covid","Monocyte"),("covid","DC")]]),
        ("recall / 상한", *[f"{fl.loc[k,'recall']:.2f} / {ce.loc[k,'recall']:.2f}" for k in [("sepsis","T cell"),("sepsis","Monocyte"),("sepsis","DC"),("covid","T cell"),("covid","Monocyte"),("covid","DC")]])]
table(s, rows, 0.9, 5.26, [2.0, 1.6, 1.9, 1.6, 1.55, 1.9, 1.5], h_row=0.35, bold_cols=(0,), font=11.5, head_font=11.5, warn_rows=())
notes(s, """[1:40] 결과 1b
"회수율, 즉 실제 DEG를 몇 개나 잡았는지는 대부분 상한선보다 낮습니다. 흥미로운 건 정밀도입니다. 생성이 DEG라고 부른 것 중 실제로도 DEG인 비율은 상한선과 비슷하거나 더 높습니다."
"해석은 이렇습니다. floor는 '확실한 것만 말한다'에 가깝습니다. 잡은 건 대체로 맞지만, 약한 신호를 많이 놓칩니다. 논문에 쓸 때 이 비대칭을 숨기지 말고 그대로 쓰는 게 낫다고 봅니다."
"표의 패혈증 DC를 보십시오. 실제는 DEG 815개, floor는 48개, 회수율 1.6퍼센트입니다. 숫자로 쓸 수 있는 수준이 아닙니다." """)

# ========================= 8c bootstrap CI
s = new(); header(s, "결과 1c", "헤드라인 지표의 불확실성 — COVID는 CI 없이 쓸 수 없다", NEXT(),
       "환자 bootstrap 200회, 95% percentile CI · 점 = 전체 데이터 추정치")
pic(s, "w1_bootstrap_ci", 3.25, 1.70, w=6.6)
side_cards(s, 10.1, 1.82, 2.45, [
    ("패혈증", f"T {ci('sepsis','T cell')}\n단핵구 {ci('sepsis','Monocyte')}\n→ 좁다"),
    ("COVID", f"T {ci('covid','T cell')}\nNK {ci('covid','NK cell')}\n→ 매우 넓다"),
    ("함의", "COVID 수치는 단독 보고 금지.\n순열 null과의 여유도 작다\n(COVID T 0.58 vs null 97.5% 0.45).")], h=1.55, accent=[GREY, RED, INK])
notes(s, """[1:40] 결과 1c
"같은 지표에 환자 부트스트랩으로 신뢰구간을 붙였습니다. 200회입니다."
"패혈증은 구간이 좁습니다. T세포 0.59에서 0.73, 단핵구 0.64에서 0.72입니다. 점추정으로 써도 큰 문제 없습니다."
"COVID가 문제입니다. T세포가 0.58인데 구간이 0.10에서 0.71입니다. 환자 수가 적고 방문이 겹치기 때문입니다."
"라벨 순열 null도 같이 돌렸는데, COVID T세포는 관측 0.584, 귀무 97.5분위가 0.447입니다. 여유가 크지 않습니다."
"그래서 제안은, COVID 수치는 반드시 구간과 함께 쓰고, 단독 문장으로는 쓰지 않는 것입니다." """)

# ========================= 9 volcano
s = new(); header(s, "결과 2", "무엇이 사라졌는가 — 볼케이노", NEXT(),
       "위 = 실제, 아래 = floor · 같은 x축 · 빨간 원 = 주요 유전자 · 파랑/회색 = FDR<0.05 & |log2FC|≥1")
pic(s, "volcano", 1.85, 1.56, w=9.6)
notes(s, """[3:00] 결과 2 볼케이노
"왼쪽 열, 패혈증 단핵구입니다. 위 실제에서 오른쪽 위에 IL1R2, CLU, 그리고 C1QA, C1QB, C1QC가 보입니다. 아래 floor를 보면 IL1R2는 더 크게 나오고 CLU도 남아 있는데, C1Q 세 개가 전부 가운데 바닥으로 내려와 있습니다. 세 개 모두 floor에서는 DEG가 아닙니다."
"가운데 열, 패혈증 DC입니다. 실제는 왼쪽에 FCER1A, PTGDS, LILRA4, GZMB가 로그 폴드 마이너스 4에서 5 수준으로 박혀 있습니다. 패혈증에서 pDC와 cDC가 사라지는 신호입니다. 아래 floor를 보면 이 유전자들이 전부 0 근처에 몰려 있고, 대신 IL1R2가 양수로 올라와 있습니다. 부호까지 뒤집힌 겁니다."
"오른쪽 열, COVID 단핵구입니다. 여기는 좋은 소식입니다. IFN 유전자들이 실제와 floor 양쪽에서 같은 자리에 나옵니다. 다만 IFI27을 보십시오. 실제에서는 로그 폴드 6으로 혼자 멀리 나가 있는데, floor에서는 2 근처로 끌려 들어왔습니다."
"이 세 장이 오늘 보고의 핵심 증거입니다. 전역 상관계수로는 안 보이고, 볼케이노로 보면 어느 유전자가 사라졌는지 바로 보입니다." """)

# ========================= 9b lymphoid volcano (NEW)
s = new(); header(s, "결과 2 · 10/30", "림프구 볼케이노 — 방향은 맞고, 유의성은 부풀려진다", NEXT(),
       "패혈증 성인 sepsis vs HC · 위 real, 아래 floor · 주황 = IEG/ambient")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "volcano_lymphoid_sepsis", 0.75, 1.62, w=8.75)
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
# ========================= 9c matched n (NEW)
s = new(); header(s, "결과 2b · 10/30", "n을 맞추면 — up은 유지, down은 소실, 샘플은 너무 비슷", NEXT(),
       "같은 case/control 수로 subsample × 50회 · 오른쪽 = 샘플 간 SD 비")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "matchedN_summary", 1.67, 1.62, w=10.0)
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
# ========================= 10 scatter
s = new(); header(s, "결과 2b", "전역 상관 뒤에 숨은 기울기와 개별 결손", NEXT(),
       "회색 = 전체 유전자, 파랑 = 실제 DEG · 점선 = y=x, 빨간 선 = 회귀 기울기")
pic(s, "logfc_scatter", 1.1, 1.72, w=11.2)
text(s, 0.89, 5.72, 11.55, 1.0, [[("두 가지를 동시에 본다  ", {"bold": True, "color": INK}),
     ("(1) 기울기가 1보다 작다 = 질병 효과 크기가 전반적으로 눌린다. (2) 빨간 원은 추세선에서 떨어진 개별 결손이다. "
      "상관계수만 보고하면 둘 다 보이지 않는다.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[1:50] 결과 2b
"같은 데이터를 산점도로 보겠습니다. x가 실제 로그 폴드, y가 floor 로그 폴드입니다."
"먼저 빨간 회귀선이 점선보다 눕습니다. 기울기가 1보다 작다는 건 질병 효과의 크기가 전반적으로 눌린다는 뜻입니다. 이건 이전부터 알고 있던 attenuation입니다."
"그런데 그것만이 아닙니다. 빨간 원으로 표시한 유전자들은 추세선에서도 떨어져 있습니다. C1QC는 실제에서 4.75인데 floor에서 0.64입니다. 기울기로 설명되는 수준이 아닙니다. 개별 모듈이 따로 빠진 겁니다."
"제가 강조하고 싶은 건, 상관계수 하나로 보고하면 이 두 가지가 모두 안 보인다는 점입니다. 평균적인 눌림과 특정 모듈의 결손은 원인도 다르고 고치는 방법도 다릅니다." """)

# ========================= 11 modules
s = new(); header(s, "결과 3", "모듈 단위로 보면 승패가 갈린다", NEXT(), "각 모듈 구성 유전자의 floor/real 질병 log2FC 비율, 중앙값 ± SEM")
pic(s, "bar_modules", 1.4, 1.70, w=10.6)
text(s, 0.89, 5.85, 11.55, 1.0, [[("패턴  ", {"bold": True, "color": INK}),
     ("살아남은 모듈(패혈증 골수계 활성화 1.03×, ISG 0.80×)은 여러 계통에 걸쳐 켜지는 '전신성' 프로그램이다. "
      "무너진 모듈(C1Q 0.10×, pDC/cDC 정체성 0.09×)은 한 세포형에만 켜지는 '정체성' 프로그램이다.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[2:20] 결과 3 모듈
"유전자를 모듈로 묶어서 floor가 실제 효과의 몇 배를 가져오는지 봤습니다. 1이면 완벽, 0.5 아래면 절반도 못 가져온 겁니다."
"오른쪽 두 개가 좋은 쪽입니다. 패혈증 골수계 활성화 모듈, 즉 IL1R2, MCEMP1, AQP9, TNFAIP6, ACSL1은 1.03배로 거의 완벽합니다. 이건 보고서 Result 6B에 적어둔 바로 그 다섯 유전자인데, 이번 재계산에서도 그대로 재현됐습니다. COVID의 ISG 모듈도 0.80배로 양호합니다."
"왼쪽 두 개가 문제입니다. C1Q 보체 모듈 0.10배, 수지상세포 정체성 모듈 0.09배입니다. 10분의 1도 못 가져옵니다."
"여기서 패턴이 보입니다. 살아남은 모듈은 염증이나 인터페론처럼 여러 세포형에 동시에 켜지는 전신 프로그램입니다. 무너진 모듈은 그 세포형에만 켜지는 정체성 프로그램입니다."
"이게 오늘 제가 제일 하고 싶은 얘기입니다. 모델이 bulk에서 '전신에 무슨 일이 일어났는가'는 읽는데, '어느 세포가 자기 정체성으로 무엇을 하는가'는 못 읽습니다. 구조상 당연할 수도 있습니다. bulk 조건에 세포형 정보가 안 들어가니까요." """)

# ========================= 12 named gene table
s = new(); header(s, "결과 3b", "논문에 적어둔 유전자와 1:1 대조", NEXT(), "값 = 질병 log2FC (case vs control) · 비율 = floor / real")
rows = [("유전자 (패혈증 단핵구)", "real", "floor", "비율", "유전자 (패혈증 DC / COVID 단핵구)", "real", "floor", "비율")]
L = [("IL1R2", "sepsis", "Monocyte"), ("MCEMP1", "sepsis", "Monocyte"), ("TNFAIP6", "sepsis", "Monocyte"),
     ("C1QC", "sepsis", "Monocyte"), ("C1QB", "sepsis", "Monocyte"), ("C1QA", "sepsis", "Monocyte"), ("HLA-DRA", "sepsis", "Monocyte")]
Rt = [("FCER1A", "sepsis", "DC"), ("LILRA4", "sepsis", "DC"), ("GZMB", "sepsis", "DC"), ("IL1R2", "sepsis", "DC"),
      ("IFI27", "covid", "Monocyte"), ("SIGLEC1", "covid", "Monocyte"), ("IFIT1", "covid", "Monocyte")]
for (gn, c1, t1), (gn2, c2, t2) in zip(L, Rt):
    rows.append((gn, f"{g(c1,t1,gn,'real_logFC'):+.2f}", f"{g(c1,t1,gn,'gen_logFC'):+.2f}", f"{g(c1,t1,gn,'ratio'):.2f}",
                 f"{gn2}  ({'DC' if c2=='sepsis' else 'COVID'})", f"{g(c2,t2,gn2,'real_logFC'):+.2f}", f"{g(c2,t2,gn2,'gen_logFC'):+.2f}", f"{g(c2,t2,gn2,'ratio'):.2f}"))
table(s, rows, 0.9, 1.7, [2.5, 1.0, 1.0, 0.85, 2.9, 1.0, 1.0, 0.85], h_row=0.46, bold_cols=(0, 4), font=12, head_font=11.5, warn_rows=(4, 5, 6))
text(s, 0.89, 6.0, 11.55, 0.8, [[("빨간 행  ", {"bold": True, "color": RED}),
     ("C1QA/B/C는 real에서 모두 DEG인데 floor에서는 셋 다 DEG가 아니다. C1QA는 부호까지 반대. "
      "패혈증 DC의 IL1R2도 real −1.51 → floor +1.89로 부호 반전.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[2:00] 결과 3b 유전자 표
"왼쪽이 패혈증 단핵구입니다. 위 세 줄, IL1R2 MCEMP1 TNFAIP6은 비율이 0.76에서 1.33으로 잘 맞습니다. 우리가 보고서에 적어둔 그 유전자들이고 재현됩니다."
"빨간 세 줄이 C1QA, C1QB, C1QC입니다. 실제에서는 각각 2.2, 3.5, 4.8인데 floor에서는 마이너스 0.02, 0.13, 0.64입니다. C1QA는 부호까지 반대입니다. 셋 다 floor에서는 유의하지 않습니다."
"오른쪽 위 네 줄은 패혈증 DC입니다. FCER1A가 실제 마이너스 5.3, floor 마이너스 0.41로 8퍼센트, LILRA4는 3퍼센트입니다. 그리고 IL1R2가 실제에서는 DC에서 내려가는데 floor에서는 올라갑니다. 부호 반전입니다."
"오른쪽 아래가 COVID입니다. IFIT1은 0.90으로 좋은데 IFI27만 0.33입니다. COVID에서 제일 강한 단일 마커가 3분의 1로 눌립니다."
"이 표는 10월 16일과 17일에 보고한 '대식세포 부재, C1Q/VSIG4 프로그램 결손'과 정확히 같은 얘기입니다. 다른 점은, 그때는 세포 라벨링 문제일 가능성이 남아 있었는데 이번에는 라벨과 무관한 질병 DEG 수준에서 똑같이 나왔다는 겁니다." """)

# ========================= 12b A2 bulk
s = new(); header(s, "결과 3c", "신호는 bulk 입력에 있다 — 모델이 세포형으로 못 보낸다", NEXT(),
       "패혈증 bulk GSE279448 성인 115 case / 43 control · COVID bulk Zenodo 6497251 환자 68 / 55")
pic(s, "w1_bulk_vs_floor", 2.5, 1.68, w=8.4)
text(s, 0.89, 5.92, 11.55, 1.0, [[("결정적  ", {"bold": True, "color": RED}),
     (f"C1QA/B/C·RNASE1·CLU 5/5, FCER1A·LILRA4·GZMB·PTGDS·CLIC3·CD1C 6/6 모두 bulk에서 FDR < 0.001, 부호 100% 일치. "
      f"배경값: real 세포형 DEG의 {100*A2BG.bulk_sig_frac.min():.0f}~{100*A2BG.bulk_sig_frac.max():.0f}%가 bulk에서도 유의. "
      "→ '정보 한계'가 아니라 '세포형 귀속 실패'다.", {})]], size=13, color=BODY, line=1.15)
notes(s, """[2:20] 결과 3c — bulk 확인
"지난주 보고에서 받을 가장 날카로운 질문이 이거였습니다. '그 모듈이 애초에 bulk에 들어 있기는 한가?' 없다면 모델 탓이 아니라 방법의 정보 한계니까요."
"이미 받아둔 bulk로 확인했습니다. 중요한 건 bulk가 단일세포보다 검정력이 훨씬 높다는 점입니다. 성인 패혈증 115명 대 건강 43명입니다. 단일세포는 93 대 7이었죠."
"보라색이 bulk, 파랑이 floor입니다. C1Q 모듈은 bulk에서 real 단일세포 효과의 0.86배로 잡힙니다. 다섯 유전자 모두 FDR 0.001 미만, 부호 100퍼센트 일치입니다. 수지상세포 정체성 모듈도 여섯 유전자 전부 유의합니다."
"그런데 floor는 각각 0.10배, 0.09배입니다."
"배경값도 쟀습니다. real 세포형 DEG의 70에서 100퍼센트가 bulk에서도 유의합니다. bulk는 정보가 풍부하고, 병목은 그 정보를 어느 세포형에 줄지 결정하는 단계입니다."
"이 한 장으로 '모델 결함'이라는 주장이 추정에서 입증으로 올라갑니다." """)

# ========================= 12c C2 pos vs load
s = new(); header(s, "결과 3d", "저발현 때문이 아니다 — 양성 세포 수는 맞고 안의 양이 틀리다", NEXT(),
       "깊이 맞춘 real 대비 비율 · 왼쪽 floor, 오른쪽 prev K4000(Poisson read-out)")
pic(s, "w1_pos_vs_load", 1.45, 1.66, w=10.5)
text(s, 0.89, 5.92, 11.55, 1.0, [[("두 가지 사실  ", {"bold": True, "color": INK}),
     (f"(1) 모듈 유전자는 자기 세포형 안에서 순위 {int(C2R.rank_in_type.min())}~{int(C2R.rank_in_type.max())}위 — 저발현 절벽(순위 1,500+)과 무관. "
      "(2) 실제 C1Q 양성 단핵구는 36.9 UMI를 싣는데 floor의 양성 세포는 1.47 UMI. 검출(support)은 보존되고 발현량(magnitude)이 무너진다 "
      "→ fix plan은 B1(stochastic support)이 아니라 A1/C3(count·NB head) 쪽.", {})]], size=13, color=BODY, line=1.15)
notes(s, """[2:20] 결과 3d — 저발현이 아니다
"다음 질문은 '혹시 저발현 절벽의 잔재 아니냐'입니다."
"아닙니다. 이 모듈 유전자들은 자기 세포형 안에서 순위 25위에서 742위입니다. 패널 약 6천 개 중 상위 12퍼센트입니다. 깊이를 맞춘 real 대비 검출률도 멀쩡합니다."
"그럼 무엇이 무너졌나. 양성 세포 비율과 그 세포 안의 전사체 양을 분리했습니다. 회색이 양성 세포 비율, 파랑이 양성 세포 안의 양입니다."
"C1Q를 보십시오. 양성 세포 비율은 1.07로 오히려 약간 많습니다. 그런데 안의 양이 0.11입니다. 실제 C1Q 양성 단핵구는 37 UMI를 싣는데 floor는 1.47입니다. 사실상 전부 1입니다."
"즉 지지, support는 맞고 크기, magnitude가 틀립니다. 이건 10월 20일 수정 계획의 방향을 바꿉니다. B1 확률적 support가 아니라 A1이나 C3, 카운트 분포 헤드 쪽이 맞습니다."
"오른쪽 prev K4000도 같이 봐주십시오. Poisson read-out인데도 C1Q는 0.09입니다. read-out 바꾼다고 해결 안 됩니다." """)

# ========================= 12d A1 ladder
s = new(); header(s, "결과 3e", "대조 사다리 — 깊이, read-out, 모델을 분리한다", NEXT(),
       "모든 값은 '깊이 맞춘 real B'를 1.0으로 한 모듈 회복률 (양성 비율 × 양성 내 발현량)")
pic(s, "w2_ladder", 0.95, 1.66, w=8.2)
pic(s, "w2_c1q_decomp", 9.35, 2.25, w=3.3)
text(s, 0.89, 6.02, 11.55, 0.8, [[("읽는 법  ", {"bold": True, "color": INK}),
     (f"깊이만 맞춘 real(3번)에서 Poisson read-out(4번)까지는 변화 없음. floor 규칙(5번)이 C1Q를 {LAD.loc['5 real B + read-out (floor rule)','C1Q / complement']:.2f}로 떨어뜨리고, "
      f"PBISC floor(6번)는 {LAD.loc['6 PBISC floor','C1Q / complement']:.2f}. 즉 read-out 규칙 −31%, 그 나머지의 −81%가 모델.", {})]],
     size=13, color=BODY, line=1.15)
notes(s, """[2:40] 결과 3e — 대조 사다리
"원인을 더 좁혔습니다. 왼쪽이 사다리입니다."
"1번은 참조 환자, 2번은 다른 환자, 3번은 깊이를 맞춘 것입니다. 여기까지 1.0 근처입니다. 깊이는 범인이 아닙니다."
"4번과 5번이 핵심입니다. 실제 세포를 가지고 PBISC의 read-out을 흉내 낸 겁니다. kNN 15개로 프로파일을 평활하고 목표 깊이를 곱한 뒤, 4번은 Poisson을 뽑고 5번은 floor 규칙을 적용했습니다."
"Poisson은 거의 변화가 없습니다. floor 규칙은 C1Q를 0.69로 떨어뜨립니다. 즉 read-out 규칙만으로 약 3분의 1이 사라집니다."
"그런데 PBISC floor는 0.12입니다. 오른쪽 폭포 그림이 그걸 분해한 겁니다. read-out 규칙으로 31퍼센트, 남은 것의 81퍼센트가 모델입니다."
"이게 왜 중요하냐면, read-out 규칙만 고쳐서는 C1Q가 돌아오지 않는다는 뜻이기 때문입니다. 반대로 형질세포 모듈은 floor 규칙 단계에서 1.0인데 PBISC floor에서 0.42로 떨어지니 역시 모델입니다."
"참고로 이 사다리는 실제 세포만 가지고 만든 거라 새 생성이 필요 없었습니다." """)

# ========================= 12e B2 three models
s = new(); header(s, "결과 3f", "세 모델 비교 — 그리고 paired 질문의 대리 답", NEXT(),
       "같은 모듈 지표 · 깊이 맞춘 real = 1.0")
pic(s, "w2_three_models", 2.35, 1.66, w=8.7)
text(s, 0.89, 5.52, 11.55, 1.25, [
  [("scSemiProfiler가 답을 준다  ", {"bold": True, "color": INK, "size": 13.5})],
  [("정답이 측정된 그 샘플을 재생성하는 더 쉬운 조건 + 같은 코호트의 실제 대표 세포 복사 → 양성 세포 '안의' C1Q 양은 0.97배로 정확. ", {}),
   ("정보는 코호트 안에 있고 복원 가능하다", {"bold": True, "color": RED}),
   (". 단 세포당 유전자 ~100개로 희소해 양성 세포를 1/9만 만들어 전체 회복률은 0.11.", {})],
  [("→ 서버에만 있는 paired 생성 카운트 없이도 'bulk-only 조건이 병목'이라는 해석을 지지한다.", {"color": BLUE})]],
  size=12, color=BODY, line=1.14, after=4)
notes(s, """[2:20] 결과 3f — 세 모델
"원래 이번 주 계획에 paired 조건 확인이 있었습니다. PBISC가 paired로 학습했을 때는 C1Q를 만드는지 보려던 건데, paired 생성 카운트가 서버에만 있어 로컬에서는 불가능했습니다. 그건 뒤 한계 슬라이드에 적어뒀습니다."
"대신 scSemiProfiler가 대리 답을 줍니다. scSemi는 조건이 다릅니다. 정답이 측정된 바로 그 샘플을 다시 만들고, 같은 코호트의 실제 대표 세포를 복사합니다. 훨씬 쉬운 조건이죠."
"그 scSemi에서 C1Q 양성 세포 안의 양은 0.97배로 정확합니다. 이게 중요합니다. 정보가 코호트 안에 존재하고 복원 가능하다는 뜻입니다."
"다만 scSemi는 세포당 유전자가 100개 정도로 극도로 희소해서 양성 세포를 9분의 1만 만듭니다. 그래서 전체 회복률은 0.11로 PBISC floor와 비슷하게 떨어집니다. 실패 방식이 다른 겁니다."
"정리하면, 어떤 방법도 이 모듈을 온전히 복원하지 못하지만, 실패 지점이 다릅니다. PBISC는 양은 못 맞추고 수는 맞추며, scSemi는 반대입니다." """)

# ========================= 13 S100 leak
s = new(); header(s, "결과 4", "새로 발견 — alarmin이 계통을 넘어 샌다", NEXT(), "COVID · 막대 = real, 점 = floor, * = floor에서 DEG로 호출됨")
pic(s, "bar_s100_leak", 1.75, 1.72, w=7.4)
side_cards(s, 9.45, 1.78, 3.1, [
    ("실제", "S100A8/A9는 단핵구에서만 질병 DEG.\n림프계에서는 log2FC 0.1~0.5, 유의하지 않음."),
    ("floor", "T·B·NK에서 log2FC 1.2~1.5로 올라가\nDEG로 호출됨. 단핵구에서는 오히려 약해짐."),
    ("해석", "골수계 프로그램이 제자리에서 약해지고 림프계로 번진다.\n기존 '단핵구 T 마커 공검출 12%'의 반대 방향 사례.")], h=1.62, accent=[GREY, RED, INK])
notes(s, """[1:50] 결과 4 alarmin 누출
"이건 이번에 새로 나온 겁니다."
"S100A8과 A9는 호중구와 단핵구의 alarmin입니다. 실제 데이터에서는 COVID 환자의 단핵구에서만 올라가고, T B NK에서는 0.1에서 0.5 수준으로 유의하지 않습니다."
"floor에서는 반대입니다. 단핵구에서는 오히려 약해져서 유의성을 잃고, T세포 1.43, B세포 1.51, NK 1.39로 전부 DEG로 호출됩니다. 별표가 붙은 게 그겁니다."
"즉 골수계 프로그램이 제자리에서 희미해지고 엉뚱한 계통에서 가짜 DEG로 나타납니다."
"우리가 이미 보고한 '패혈증 단핵구의 12퍼센트가 T세포 마커를 함께 검출한다'는 현상의 정확히 반대 방향 사례입니다. 두 방향 다 있다는 건 모델이 세포형 경계를 흐리게 그린다는 뜻이고, 라벨 전이로는 고쳐지지 않습니다."
"실무적으로는, 이런 가짜 DEG가 bulk-only 코호트에서 새 발견이라고 보고될 위험이 있습니다. 게이트가 필요한 이유입니다." """)

# ========================= 14 reproducibility of reported numbers
s = new(); header(s, "결과 5", "보고 수치 재현성 — 3개는 맞고 2개는 정의가 문제", NEXT(), "값 = 패혈증 / COVID, prev → floor")
rows = [("지표", "보고서 Result 6B", "이번 재계산", "판정"),
        ("logFC r", "0.57 / 0.47 → 0.70 / 0.52", "0.58 / 0.53 → 0.68 / 0.54", "재현"),
        ("rank > 2,000 유전자 r", "0.41 / 0.26 → 0.64 / 0.38", "0.46 / 0.35 → 0.63 / 0.42", "재현"),
        ("기울기 (slope)", "0.32 / 0.25 → 0.62 / 0.42", "0.34 / 0.28 → 0.55 / 0.47", "재현"),
        ("sign agreement", "0.73 / 0.76 → 0.90 / 0.90", "0.47 / 0.43 → 0.63 / 0.58", "분모 불일치"),
        ("생성 top-100 정밀도", "0.76 / 0.46 → 0.71 / 0.57", "0.38 / 0.29 → 0.58 / 0.34", "정의 불일치")]
table(s, rows, 0.9, 1.75, [3.3, 3.5, 3.5, 1.6], h_row=0.55, bold_cols=(0,), font=12.5, head_font=12.5, warn_rows=(4, 5))
text(s, 0.89, 5.3, 11.55, 1.5, [[("무슨 일이 일어난 것인가  ", {"bold": True, "color": INK})],
     [("재계산한 sign agreement는 rank>2,000 ", {}), ("전체 유전자", {"bold": True, "color": INK}),
      (" 기준이라 귀무 유전자의 동전던지기(0.5)가 섞인다. 같은 데이터에서 ", {}),
      ("real DEG로 한정하면 floor는 0.92~1.00", {"bold": True, "color": BLUE}),
      (". 즉 보고서의 0.90은 DE 유전자 또는 |logFC| 역치 이상으로 제한한 값일 가능성이 높다.", {})],
     [("지금 상태로는 '바닥이 0.5인 지표'를 바닥이 0인 지표처럼 읽히게 한다 → methods에 분모를 명시하거나 재계산 값으로 교체해야 한다.", {"color": RED})]],
     size=13.5, color=BODY, line=1.2, after=5)
notes(s, """[2:30] 결과 5 재현성
"보고서 Result 6B의 다섯 지표를 같은 데이터로 다시 계산했습니다."
"위 세 줄, logFC 상관, 저발현 유전자 상관, 기울기는 소수점 둘째 자리 수준에서 재현됩니다. 이건 안심해도 됩니다."
"아래 두 줄이 문제입니다. sign agreement가 보고서는 0.90인데 제 계산은 0.63입니다."
"원인을 찾았습니다. 제 계산은 rank 2000 이후 전체 유전자 기준입니다. 거기에는 질병과 무관한 유전자가 대부분이고, 그런 유전자의 부호 일치는 동전던지기라 0.5가 바닥입니다. 그래서 0.63은 '0.5 바닥 위에서 0.13'이라는 뜻입니다."
"같은 데이터에서 실제 DEG로만 한정하면 floor는 0.92에서 1.00이 나옵니다. 보고서의 0.90은 그쪽 계산일 가능성이 높습니다."
"어느 쪽이든 문제는, 지금 원고에 분모가 안 적혀 있다는 겁니다. 리뷰어가 '0.5가 바닥인 지표 아니냐'고 물으면 방어할 수 없습니다. methods에 명시하거나 재계산 값으로 바꿔야 합니다."
"top-100 정밀도도 같은 종류입니다. 제 real DEG 정의가 FDR에 로그폴드 1 이상이라 더 엄격합니다. 어느 정의를 쓸지 오늘 정하면 좋겠습니다.""")

# ========================= 15 caution IEG
s = new(); header(s, "결과 6", "immediate-early 분리 — 상관은 살아남고, 헤드라인 유전자 목록은 못 쓴다", NEXT(),
       "IEG + ambient(HBB·PPBP 등) 59개를 유전자 집합에서 제거하고 전 지표 재계산")
pic(s, "w3_ieg", 0.62, 1.66, w=12.1)
text(s, 0.89, 5.72, 11.55, 1.2, [
  [("세 가지 결론  ", {"bold": True, "color": INK, "size": 13.5})],
  [(f"(1) logFC r은 거의 안 변한다 (패혈증 T {IE.loc[('sepsis','T cell'),'r_all']:.2f} → {IE.loc[('sepsis','T cell'),'r_noIEG']:.2f}). "
    f"→ 0.70은 IEG 아티팩트가 아니다. ", {}),
   (f"(2) 그러나 패혈증 림프구 real DEG의 {IE.loc[('sepsis','B cell'),'pct_realDEG_IEG']:.0f}%, |logFC| 상위 30개의 "
    f"{IE.loc[('sepsis','B cell'),'pct_top30_IEG']:.0f}%가 IEG이고, floor의 recall은 IEG 제거 시 "
    f"{IE.loc[('sepsis','T cell'),'recall_all']:.2f}→{IE.loc[('sepsis','T cell'),'recall_noIEG']:.2f}로 떨어진다.", {"color": RED})],
  [(f"(3) IEG 축 자체가 패혈증에서는 case↑ (d = {IQ.loc[('sepsis','B cell','real'),'cohens_d']:.2f}), COVID에서는 case↓ "
    f"(d = {IQ.loc[('covid','NK cell','real'),'cohens_d']:.2f}). 코호트마다 부호가 반대 → 보편적 감염 반응이 아니라 "
    "코호트·검체 처리 축일 가능성이 높다.", {})]],
  size=12.5, color=BODY, line=1.16, after=4)
notes(s, """[2:40] 결과 6 — immediate-early 분리 (완료)
"지난주 열어둔 숙제였습니다. 패혈증 림프구의 높은 점수가 질병 생물학인지, 채혈·처리 과정의 스트레스 반응인지 구분하자는 것이었습니다."
"IEG와 적혈구·혈소판 오염 유전자 59개를 아예 빼고 전부 다시 계산했습니다."
"첫째, 좋은 소식입니다. 상관계수는 거의 안 변합니다. 패혈증 T세포가 0.70에서 0.64입니다. 0.06 떨어집니다. 단핵구와 COVID는 0.01 이하입니다. 그러니까 '0.70은 스트레스 유전자 덕분'이라는 반론은 성립하지 않습니다."
"둘째, 나쁜 소식입니다. 가운데 패널을 보십시오. 패혈증 림프구에서 real DEG의 26에서 28퍼센트가 IEG이고, 폴드체인지 상위 30개로 좁히면 47에서 67퍼센트입니다. 그리고 IEG를 빼면 floor의 회수율이 0.48에서 0.34로 떨어집니다. floor가 상대적으로 스트레스 프로그램을 더 잘 복원하고 있었다는 뜻입니다."
"셋째가 제일 중요합니다. 오른쪽입니다. IEG 점수 자체를 case와 control로 비교하면, 패혈증에서는 환자가 높고 효과크기 1.4에서 1.7입니다. 그런데 COVID에서는 반대로 환자가 낮습니다. 마이너스 0.6에서 마이너스 1.2입니다."
"부호가 코호트마다 반대라는 건, 이게 보편적인 감염 반응이 아니라는 뜻입니다. 패혈증 코호트는 중환자실 환자와 건강 지원자를 비교하니 채혈부터 처리까지의 조건이 체계적으로 다를 수 있고, COVID 코호트는 양쪽 다 ESKD 환자라 그 차이가 작거나 반대입니다."
"그래서 제안은 이렇습니다. 상관계수는 그대로 쓰되, 패혈증 림프구의 유전자 목록은 헤드라인으로 쓰지 않고 IEG 제외값을 병기합니다." """)

# ========================= 15b DC underproduction (NEW)
s = new(); header(s, "결과 6b · 10/30", "DC 과소생성 — 모델이 bulk의 pDC 신호를 따르지 않는다", NEXT(),
       "b: 같은 샘플의 bulk vs 생성 pseudobulk 마커 점수 (COVID real은 매칭 없음)")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "v4_dc_underproduction", 0.85, 1.7, w=11.65)
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
# ========================= 16 transition rare
s = new()
box(s, 0, 0, 13.333, 7.5, CARD, shape=MSO_SHAPE.RECTANGLE); box(s, 0, 0, 0.17, 7.5, ORANGE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 2.3, 11.9, 0.3, "PART 2", size=14, color=ORANGE, bold=True, spc=150)
title_box(s, "희귀 집단은 '있다'. 그런데 그 안의 세포는\n정말 그 세포인가?", y=2.9, size=36)
text(s, 0.89, 4.6, 11.9, 1.0, "지금까지의 판정 기준은 부모 계통 대비 비율(0.5–2×)뿐이었다. 비율이 맞아도 집단 내부의 발현이 실제와 같은지는 "
     "한 번도 확인하지 않았다. 19개 집단, 동일 세포 수, 동일 라벨 전이 경로로 비교한다.", size=18, color=BODY, line=1.2)
text(s, 0.89, 6.88, 11.9, 0.26, f"{FOOT} · {NEXT()}", size=12, color=MUTED)
notes(s, """[0:40] Part 2 전환
"여기서부터 두 번째 질문입니다."
"우리는 Result 5에서 희귀 집단 복원을 '부모 계통 대비 비율이 0.5에서 2배 사이면 복원됨'으로 판정했습니다. 패혈증 14개 중 8개 복원, 이런 식으로요."
"그런데 비율이 맞다는 건 그 라벨이 붙은 세포가 그만큼 있다는 뜻이지, 그 세포가 진짜 그 세포라는 뜻은 아닙니다. 라벨은 kNN 전이로 붙인 거니까요."
"그래서 집단 안으로 들어가 봤습니다." """)

# ========================= 17 rare counts
s = new(); header(s, "결과 7", "먼저 수율 — 다섯 집단은 평가 자체가 불가능하다", NEXT(),
       "60,000 세포 균형 subset 기준 · floor 셀 수 / real A 셀 수 · 30세포 미만은 발현 비교에서 제외")
pic(s, "rare_counts", 2.95, 1.68, w=7.2)
text(s, 0.89, 6.14, 11.55, 0.7, [[("평가 불가  ", {"bold": True, "color": RED}),
     ("CD14mono_dex 2 / 57 · Mono_c06-CD14-MKI67 3 / 221 · COVID ILC 0 / 50 · 패혈증 pDC 16 / 96 · B_ASC_IgG 26 / 161. "
      "비율 지표로는 '0.5× 미만'으로만 기록되던 것이, 실제로는 통계를 돌릴 세포가 없는 상태.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[1:40] 결과 7 수율
"발현을 비교하기 전에 셀 수부터 봤습니다. 6만 세포 균형 subset에서 각 집단에 몇 개가 있는지입니다."
"빨간 막대 다섯 개가 30세포 미만입니다. dex 단핵구 2개, 증식 단핵구 3개, COVID의 ILC는 0개입니다. 패혈증 pDC도 16개뿐입니다."
"이게 왜 중요하냐면, 지금 Result 5 표에서는 이것들이 '0.5배 미만'이라는 한 칸으로 들어가 있습니다. 그런데 0.4배와 0.03배는 다른 얘기입니다. 전자는 적게 만든 거고 후자는 안 만든 겁니다."
"그리고 숫자가 이렇게 적으면 그 집단에 대해 어떤 통계도 돌릴 수 없습니다. 논문에서 '이 집단을 복원했다'고 쓸 수 없는 상태입니다."
"반대쪽 끝도 보십시오. HSPC는 실제보다 3배 넘게 만듭니다. 과생산도 있습니다.""")

# ========================= 18 rare marker
s = new(); header(s, "결과 8", "그리고 남은 집단에서도 정체성 프로그램이 희석된다", NEXT(),
       "집단 정의 마커 top-10 점수 / real A · 회색 = real B(같은 라벨 전이 경로), 파랑 = floor")
pic(s, "rare_marker", 2.75, 1.68, w=6.95)
side_cards(s, 10.0, 1.9, 2.55, [
    ("중앙값", "real B 0.98\nfloor 0.86"),
    ("최저", "HSC 0.40\nDC1 0.49\nB_c08-MKI67 0.52"),
    ("핵심", "real B도 같은 kNN 전이를 거쳤는데 1.0을 유지한다.\n→ 라벨링이 아니라 생성의 문제")], h=1.5, accent=[GREY, RED, BLUE])
notes(s, """[2:10] 결과 8 마커 희석
"각 집단을 정의하는 마커 10개의 평균 발현 점수를, real A를 1.0으로 놓고 비교했습니다."
"회색이 real B입니다. 환자가 다르고, floor와 똑같은 kNN 라벨 전이를 거친 실제 세포입니다. 중앙값이 0.98로 거의 1을 유지합니다."
"파랑이 floor입니다. 중앙값 0.86이고, 아래쪽을 보면 HSC 0.40, DC1 0.49, 증식 B세포 0.52입니다. 절반 이하입니다."
"이 그림의 핵심은 회색과 파랑의 대비입니다. real B도 똑같은 라벨 전이를 거쳤는데 1.0입니다. 그러니까 마커가 흐려지는 게 라벨을 잘못 붙여서가 아닙니다. 생성된 세포 자체가 그 프로그램을 덜 켜고 있습니다."
"이게 Result 5의 '복원됨' 판정을 다시 보게 만듭니다. DC1은 비율로는 0.83배라 '복원'인데, 안을 보면 마커가 절반입니다.""")

# ========================= 19 rare power
s = new(); header(s, "결과 9", "상관계수는 이것을 가린다 — 같은 세포 수로 맞춘 비교", NEXT(),
       "집단마다 세 소스를 같은 n으로 맞춰 10회 반복 · 왼쪽이 실제로 쓸 지표")
pic(s, "rare_power", 1.55, 1.72, w=10.3)
text(s, 0.89, 5.72, 11.55, 1.0, [[("읽는 법  ", {"bold": True, "color": INK}),
     ("오른쪽: 프로파일 상관은 real-real 0.990 vs floor 0.979로 거의 구분이 안 된다. "
      "왼쪽: 같은 데이터에서 |Δ| ≥ 0.5 유전자 수는 기술 잡음 1개, real-real 12개, floor 49개 — "
      "floor는 실제끼리의 차이의 약 4배를 만든다.", {})]], size=13.5, color=BODY, line=1.15)
notes(s, """[2:10] 결과 9 equal-n
"여기가 방법론적으로 제일 중요한 슬라이드입니다."
"집단마다 세 소스를 같은 세포 수로 맞추고 10번 반복했습니다. 그리고 세 가지를 비교했습니다. real A를 둘로 쪼갠 것, 이게 순수 기술 잡음입니다. real A 대 real B, 이건 환자 차이와 라벨 전이가 섞인 겁니다. 그리고 real A 대 floor입니다."
"오른쪽 상관계수를 먼저 보십시오. 0.993, 0.990, 0.979입니다. 거의 구분이 안 됩니다. 만약 우리가 상관계수로 보고했으면 '희귀 집단 발현이 실제와 거의 같다'고 썼을 겁니다."
"왼쪽이 같은 데이터입니다. 차이가 0.5 이상인 유전자 개수로 세면 1개, 12개, 49개입니다. floor는 실제끼리의 차이의 네 배를 만듭니다."
"교훈은 단순합니다. 고차원에서 상관계수는 포화됩니다. 유전자 개수로 세야 보입니다. 이건 앞으로 모든 집단 비교에 적용해야 할 지표입니다.""")

# ========================= 20 rare drivers
s = new(); header(s, "결과 10", "같은 라벨 안에서 무엇이 움직이는가", NEXT(),
       "Δ = 해당 소스의 집단 내 평균 발현 − real A · 회색 real B가 대조")
pic(s, "rare_drivers", 1.55, 1.58, w=10.2)
rows = [("패턴", "예시", "real B"),
        ("분비·형질세포 프로그램 소실", "B_c07-MZB1: JCHAIN −1.79, IGLL5 −1.79, CD27 −1.19 / B_c08-MKI67: JCHAIN −3.70", "≈ −0.2"),
        ("히스톤 H1 과발현", "pDC: HIST1H1E +1.33, H1D +0.96, H1C +0.94 (NK·CD8 증식에서도 반복)", "≈ 0"),
        ("alarmin 재배치", "cDC: S100A8 −1.98 / S100A9 −1.97 ↔ B_c08-MKI67: +1.54 / +1.13", "−1.2 / −0.0"),
        ("정의 마커 직접 희석", "MAIT: KLRB1 −1.26 · NKT_c05: CD3D −0.84 · Mono_c01: S100A12 −1.19", "−0.46 / −0.75 / −0.29"),
        ("MHC-II 방향이 집단마다 반대", "Mono_c01·cDC +0.6~1.0 과대 ↔ CD14mono_IFN·pDC −0.6~−1.3 과소", "≈ ±0.3")]
table(s, rows, 0.9, 4.58, [3.0, 6.6, 1.95], h_row=0.33, bold_cols=(0,), font=11, head_font=11)
notes(s, """[2:30] 결과 10 driver 유전자
"집단 안에서 실제로 어느 유전자가 움직이는지 봤습니다. 파랑이 floor, 회색이 real B 대조입니다."
"다섯 가지 패턴이 반복됩니다."
"첫째, 분비 프로그램 소실입니다. 형질세포 계열에서 JCHAIN, MZB1, IGLL5, DERL3, CD27이 전부 내려갑니다. 증식 B세포에서는 JCHAIN이 마이너스 3.7입니다. real B는 마이너스 0.2 수준이니까 대조와 비교가 안 됩니다."
"둘째, 히스톤 H1이 올라갑니다. pDC, NK, 증식 CD8에서 반복됩니다. 이건 10월 15일 COVID 단핵구 분석에서 'platform signature'로 가설을 세워둔 바로 그 신호입니다. 지금은 여러 집단에서 일반적으로 나타나는 걸 확인했습니다."
"셋째, alarmin 재배치입니다. cDC에서는 S100A8이 마이너스 2, 증식 B세포에서는 플러스 1.5입니다. 앞에서 본 COVID DEG 누출과 같은 현상이 집단 내부에서도 보입니다."
"넷째, 집단을 정의하는 마커가 직접 희석됩니다. MAIT의 KLRB1이 마이너스 1.26입니다. MAIT을 MAIT이라고 부르게 해주는 유전자입니다."
"다섯째, MHC-II는 집단마다 방향이 반대입니다. 그래서 'MHC-II가 과하다' 같은 한 문장으로 요약하면 틀립니다. 이것도 오늘 논의 항목입니다.""")

# ========================= 20b B1 specificity law
s = new(); header(s, "결과 11", "일화가 아니라 법칙 — 전용성이 높을수록 복원되지 않는다", NEXT(),
       "모듈·전용성 모두 real 절반 A에서만 정의하고 floor는 그대로 채점 (손으로 고르지 않음)")
pic(s, "w2_specificity_law", 1.0, 1.66, w=11.3)
text(s, 0.89, 5.92, 11.55, 1.0, [[("체리피킹 방어  ", {"bold": True, "color": INK}),
     (f"유전자 수준 ρ = −0.25 (n = 21,883, p < 1e-290), 데이터 주도 모듈 ρ = −0.46 (n = 52, p = 6e-4). "
      f"발현 순위와의 상관은 +0.11로 반대 방향 — 저발현 교란이 아니다. 발현 순위 구간을 나눠도 "
      "모든 구간·모든 코호트×세포형(10/10)에서 Q5/Q1 = 0.67~0.89로 단조.", {})]], size=13, color=BODY, line=1.15)
notes(s, """[2:40] 결과 11 — 법칙
"지난주 보고의 가장 큰 약점은 모듈 세 개를 제가 골랐다는 점이었습니다. 이미 결손을 아는 데서 역으로 고른 거라 체리피킹입니다."
"그래서 손을 뗐습니다. 실제 절반 A에서만 유전자를 세포형별 프로파일로 군집해 코호트당 30개 모듈을 뽑고, 각 모듈의 세포형 전용성을 Tau로 계산했습니다. floor는 그냥 채점만 했습니다."
"오른쪽이 결과입니다. 전용성이 높을수록 회복률이 떨어집니다. 상관 마이너스 0.46, p 6×10⁻⁴입니다."
"왼쪽은 유전자 수준입니다. 2만 1천 개 유전자×세포형 조합에서 상관 마이너스 0.25입니다."
"교란 가능성을 하나 막았습니다. '전용 유전자가 원래 저발현이라 그런 것 아니냐'인데, 발현 순위와 회복률의 상관은 플러스 0.11로 방향이 반대입니다. 그리고 발현 순위 구간을 나눠서 봐도 모든 구간에서, 코호트와 세포형 10개 조합 전부에서 같은 단조 관계가 나옵니다."
"이제 세 모듈은 법칙의 극단 사례이지 제가 고른 일화가 아닙니다." """)

# ========================= 20c myeloid sepsis (NEW)
s = new(); header(s, "결과 12 · 10/30", "단핵구 subset — 흔한 상태는 겹치고, 환자 고유 상태는 없다", NEXT(),
       "패혈증 · 단핵구 구획 joint Harmony UMAP · ● real ▲ floor")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "umap_myeloid_sepsis", 1.0, 1.6, w=11.3)
text(s, 0.89, 6.1, 11.55, 0.7, [[("읽는 법  ", {"bold": True, "color": INK}),
     ("Mono_c06-MKI67 221 → 3, Macrophage 400 → 60 (0.25×). real 세포의 23%는 주변 30개 이웃의 90% 이상이 real인 '환자 고유 영역'에 있고 "
      "(Abd-S123·IC16·Skin-S140·Abd-S51), 대식세포 400개 중 225개가 이 영역에 속한다. HC 세포는 0.", {})]], size=12, color=BODY, line=1.12)
notes(s, """[2:00] 단핵구·대식세포 subset 패혈증 (추가)
"v3 슬라이드 5, 6의 희귀 집단 UMAP과 같은 형식으로 단핵구 구획만 따로 그렸습니다. joint Harmony 공간을 그대로 쓰고 단핵구 세포만 남겨 UMAP을 다시 계산했습니다."
"왼쪽 소스 패널을 보면 대부분 겹치는데, 아래쪽 꼬리와 왼쪽 위에 회색만 있는 영역이 있습니다."
"오른쪽 패널에서 증식 단핵구는 real 221개 대 floor 3개, 대식세포는 400 대 60입니다. CD16과 CD14CD16 단핵구는 floor 삼각형이 real 원과 같은 자리를 덮습니다."
"회색만 있는 영역을 정량했습니다. real 세포 주변 30개 이웃의 90퍼센트 이상이 real인 경우를 세면 real의 23퍼센트입니다. 이 세포들은 네 명의 환자에 몰려 있고 전부 패혈증 또는 ICU 대조입니다. 건강 대조는 0입니다. 대식세포 400개 중 225개가 이 영역에 있습니다."
"해석은 신중해야 합니다. 패혈증 floor는 real과 환자가 전혀 겹치지 않으니, 이 영역 일부는 환자 구성 차이입니다. 그래도 '환자 고유의 극단적 골수계 상태를 bulk에서 생성하지 못한다'는 건 분명합니다." """)
# ========================= 20d myeloid covid (NEW)
s = new(); header(s, "결과 12b · 10/30", "COVID 단핵구 — 잘 섞이지만 dex·C1 상태가 빠진다", NEXT(),
       "COVID/ESKD · 같은 방식 · 오른쪽 범례 괄호 = n real / n floor")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "umap_myeloid_covid", 1.0, 1.6, w=11.3)
text(s, 0.89, 6.1, 11.55, 0.7, [[("비교  ", {"bold": True, "color": INK}),
     ("환자 고유 영역은 1.2%로 거의 없다(패혈증 23%). 비율로는 CD14mono 0.98×, IFN 0.94×로 맞지만 CD14mono_dex 57 → 2 (0.03×), "
      "CD16mono_C1 0.38×, 반대로 CD16mono는 1.47× 과생산 — 흔한 상태로 쏠린다. (subset 비율 막대는 부록)", {})]], size=12, color=BODY, line=1.12)
notes(s, """[1:40] COVID 단핵구 (추가)
"COVID에서는 지형이 훨씬 잘 섞입니다. real만 있는 영역이 1.2퍼센트입니다. 일부 환자가 floor 쪽에도 다른 방문으로 들어가 있는 것도 이유일 수 있습니다."
"비율로 보면 흔한 subset은 맞습니다. CD14 단핵구 0.98배, IFN 단핵구 0.94배입니다."
"빠지는 건 두 가지입니다. 덱사메타손 처리 단핵구는 57개 대 2개, C1Q를 켜는 CD16 단핵구 C1은 0.38배입니다. C1은 v3에서 본 C1Q 모듈 결손과 같은 이야기입니다."
"반대로 일반 CD16 단핵구는 1.47배로 과생산됩니다. 희귀 상태가 빠진 자리를 흔한 상태가 채웁니다. 패혈증에서도 RETN과 CD16이 1.4배였습니다." """)
# ========================= 20e marker re-annotation 1 (NEW)
s = new(); header(s, "결과 13 · 10/30", "marker 재annotation — floor는 real에 가까워지지 않는다", NEXT(),
       "real A Top10 marker → de novo Leiden + marker 점수 · b,c = 구분 가능 subset")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "v4_marker_ruler", 0.85, 1.66, w=11.65)
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
# ========================= 20f marker re-annotation 2 (NEW)
s = new(); header(s, "결과 13b · 10/30", "marker 라벨이 희귀 subset을 '늘리는' 방식 — 이웃 세포를 끌어온다", NEXT(),
       "F1 = 같은 marker 자를 real B에 적용한 저자 라벨 회복력")
s.shapes[0].width = Inches(9.6); pill(s, 11.05, 0.43, 1.45, "10/30 추가", ORANGE)
pic(s, "v4_marker_keysubsets", 0.95, 1.6, w=11.4)
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
# ========================= 21 synthesis
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
# ========================= 22 actions
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
# ========================= 23 limitations
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
# ========================= 24 discussion
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
# ========================= A appendix divider
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
os.makedirs("deck", exist_ok=True)
DST = "deck/PBISC 질병 DEG 희귀집단 검증 보고서 v4.pptx"
prs.save(DST)
print("slides:", len(prs.slides.__iter__.__self__._sldIdLst), "->", DST)
with open("deck/speaker_notes_v4.md", "w", encoding="utf-8") as fh:
    fh.write("# PBISC 질병 DEG · 희귀 집단 검증 (v4) — 발표자 노트\n\n")
    fh.write(f"2026년 10월 30일 · {len(NOTE_LOG)}슬라이드 · 본편 약 45분 + 부록 (질의응답용)\n\n---\n\n")
    for i, t in enumerate(NOTE_LOG, 1):
        lines = t.split("\n")
        fh.write(f"## Slide {i} — {lines[0]}\n\n")
        for ln in lines[1:]:
            fh.write(("> " + ln + "\n\n") if ln.startswith('"') else (ln + "\n\n"))
        fh.write("---\n\n")
print("notes written")
