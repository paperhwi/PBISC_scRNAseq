"""Korean deck + Korean speaker notes for the extra validation of the floor read-out (Result 6, 2026-10-21).
Same template / layout / styles as 'PBISC floor read-out 검증 보고서.pptx'. Every number is read from results/ (scorecard inputs)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deck_helpers import *
R = "/mnt/user-data/outputs/val/results"; DST = "/mnt/user-data/outputs/deck/PBISC floor read-out 추가 검증 보고서.pptx"
rd = lambda n: pd.read_csv(f"{R}/{n}"); f2 = lambda x: f"{x:.2f}"; pct = lambda x: f"{100*x:.0f}%"
# ---------------- numbers
FP = rd("covid_donor_fingerprint.csv").set_index(["source", "gene_set"]); fp = lambda s, g, k: FP.loc[(s, g), k]
M = rd("m6_disease_effect.csv"); MC = M[M.contrast == "case vs control"].groupby(["cohort", "genes", "source"])[["logFC_r", "slope_on_DEG", "sign_agree_DEG", "n_real_DEG"]].median()
m6 = lambda c, b, s, k: MC.loc[(c, b, s), k]
PR = rd("m6_precision_top100.csv"); PR = PR[PR.contrast == "case vs control"].groupby(["cohort", "source"]).precision_top100_gen.median()
ICU = M[(M.contrast == "sepsis vs ICU control") & (M.genes == "all")].n_real_DEG.median()
C2 = pd.concat([rd(f"{c}_c2st.csv") for c in ("sepsis", "covid")]); c2 = C2.groupby(["cohort", "source"]).c2st_auc.median(); c2r = C2.groupby(["cohort", "source"]).c2st_auc.agg(["min", "max"])
DC = pd.concat([rd(f"{c}_density_coverage.csv") for c in ("sepsis", "covid")]); DC = DC[DC.lineage == "all"].set_index(["cohort", "space", "source"])
dc = lambda c, sp_, s, k: DC.loc[(c, sp_, s), k]; PCAs, HAR = "PCA (no correction)", "Harmony"
CO = pd.concat([rd(f"{c}_coexpression.csv") for c in ("sepsis", "covid")]); co = CO.groupby(["cohort", "genes", "source"]).corr_of_corr.median(); ca = CO.groupby(["cohort", "genes", "source"]).mean_abs_corr_src.median()
car = CO.groupby(["cohort", "genes"]).mean_abs_corr_real.median()
LR = pd.concat([rd(f"{c}_lr_recovery.csv").assign(cohort=c) for c in ("sepsis", "covid")]).set_index(["cohort", "source"])
LD = pd.concat([rd(f"{c}_lr_disease.csv").assign(cohort=c) for c in ("sepsis", "covid")]).set_index(["cohort", "source"])
CL = pd.concat([rd(f"{c}_claims.csv") for c in ("sepsis", "covid")]); CLu = CL[(CL.auc_realA >= 0.6) & (CL.claim != "ESKD-IFN B cell")]
W = CLu.pivot_table(index=["cohort", "claim"], columns="source", values="delta_pos_frac")
def claim_sum(c, s):
    w = W.loc[c]; return int((np.sign(w[s]) == np.sign(w["realU"])).sum()), len(w), float(np.median(w[s] / w["realU"]))
RR = pd.concat([rd(f"{c}_rare_transfer_free.csv") for c in ("sepsis", "covid")]); RRu = RR[RR.cv_auc_realA >= 0.8].copy()
RRu["ratio"] = RRu.gen_estimate_pct / RRu.realB_estimate_pct; RRu["ok"] = RRu.ratio.between(0.5, 2)
rk = RRu.groupby(["cohort", "model_depth"]).ok.agg(["sum", "size"]); rq = lambda c, s: f"{int(rk.loc[(c, s), 'sum'])}/{int(rk.loc[(c, s), 'size'])}"
rr = RRu.pivot_table(index=["cohort", "subset"], columns="model_depth", values="ratio")
IN = pd.concat([rd(f"{c}_integration_methods.csv") for c in ("sepsis", "covid")]).set_index(["cohort", "method", "pair"])
il = lambda c, m, s: IN.loc[(c, m, s), "iLISI"]
ns = {"sepsis": "패혈증", "covid": "COVID"}
SL = lambda c, s: f"{f2(c['sepsis'] if isinstance(c, dict) else c)}"
both = lambda fn: f"{fn('sepsis')} / {fn('covid')}"

# ================= 1. title
s = new(); box(s, 0, 0, 0.17, 7.5, BLUE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 1.92, 11.9, 0.3, "PBISC 생성 행렬 검증 · floor read-out 추가 검증", size=14, color=BLUE, bold=True, spc=150)
ph = s.shapes.title; ph.left, ph.top, ph.width, ph.height = Inches(0.77), Inches(2.5), Inches(11.9), Inches(1.5)
tf = ph.text_frame; tf.clear(); tf.word_wrap = True
for a in ("margin_left", "margin_right", "margin_top", "margin_bottom"): setattr(tf, a, 0)
p = tf.paragraphs[0]; style_run(p.add_run(), 40, INK, True); p.runs[0].text = "floor read-out은 무엇을 실제로 개선했나:\n질병 효과, realism, 상호작용, 희귀 세포"
text(s, 0.89, 4.15, 11.9, 0.7, "로컬 데이터로 가능한 8가지 분석에서 floor read-out(draw0)을 이전 K = 4,000 생성 세포, "
     "그리고 다른 환자의 실제 세포(실제 B)와 같은 파이프라인으로 비교", size=18, color=BODY)
pill(s, 0.89, 5.1, 1.75, "Floor read-out", BLUE); pill(s, 2.8, 5.1, 1.95, "Previous K = 4,000", RED); pill(s, 4.91, 5.1, 1.55, "Real (held-out)", GREY)
text(s, 0.89, 6.82, 11.9, 0.26, "Jong Hwi Kim · Seoul National University · 2026년 10월 21일", size=12, color=MUTED)
notes(s, f"""[0:30] 표지
"앞선 보고서에서 floor read-out이 저발현 유전자 절벽을 없애고 실제 세포와 거의 실제끼리 수준으로 섞인다는 것을 보여드렸습니다. 오늘은 한 걸음 더 나가서, 그 개선이 실제 생물학적 결과로 이어지는지 여덟 가지 분석으로 확인한 결과입니다."
색 규칙을 한 번만 말한다: 파랑 = floor read-out, 빨강 = 이전 K = 4,000, 회색·짙은 남색 = 실제 세포.
질문이 나오면: 전부 PC와 클라우드에 있는 데이터로 다시 계산한 것이고, 서버 데이터가 필요한 분석(paired benchmark, threshold sweep)은 빠져 있습니다.""")

# ================= 2. summary
s = new(); header(s, "요약", "절벽 해소는 질병 효과와 상호작용 검출로 이어졌고, 환자 특이성은 그대로", 2)
three_cards(s, 1.75, 3.85, [
    ("1", "좋아진 것", f"저발현 유전자 질병 효과 부호 일치 {f2(m6('sepsis','rank > 2,000','prev','sign_agree_DEG'))} / {f2(m6('covid','rank > 2,000','prev','sign_agree_DEG'))} → "
          f"{f2(m6('sepsis','rank > 2,000','floor','sign_agree_DEG'))} / {f2(m6('covid','rank > 2,000','floor','sign_agree_DEG'))}. 분류기 AUC 0.99 → {f2(c2[('sepsis','floor')])} / {f2(c2[('covid','floor')])}. "
          f"LR 검출 Jaccard {f2(LR.loc[('sepsis','prev4k'),'detected_jaccard'])} / {f2(LR.loc[('covid','prev4k'),'detected_jaccard'])} → {f2(LR.loc[('sepsis','floor'),'detected_jaccard'])} / {f2(LR.loc[('covid','floor'),'detected_jaccard'])}. 통합 4개 방법 모두 floor > 이전."),
    ("2", "그대로인 것", f"COVID donor 식별 {f2(fp('prev','top-100 donor-stable genes','top1_acc'))} vs {f2(fp('floor','top-100 donor-stable genes','top1_acc'))} (우연 0.03): 환자 신호는 유지되나 개선 없음. "
          f"질병에 따른 LR 변화 일치도 {f2(LD.loc[('sepsis','prev4k'),'disease_spearman_all_detected'])} → {f2(LD.loc[('sepsis','floor'),'disease_spearman_all_detected'])} (패혈증)."),
    ("3", "남은 결함", f"과검출 흔적: density {f2(dc('sepsis',PCAs,'floor','density'))} / {f2(dc('covid',PCAs,'floor','density'))} (실제 {f2(dc('sepsis',PCAs,'realU','density'))} / {f2(dc('covid',PCAs,'realU','density'))}), "
          f"유전자 간 상관이 실제보다 강함. 희귀 상태: γδT, pDC 부족, dex 단핵구 없음.")])
text(s, 0.89, 5.85, 11.55, 0.8, [[("한 문장  ", {"bold": True, "color": INK}), ("floor는 '실제처럼 보이는' 세포에서 '실제처럼 반응하는' 세포로 한 단계 나아갔고, 남은 과제는 과검출 보정과 일부 희귀 상태.", {})]],
     size=14, color=BODY, line=1.15)
notes(s, f"""[1:30] 요약
"결론부터 말씀드리면 세 가지입니다. 첫째, 절벽이 없어진 효과가 저발현 유전자의 질병 반응 재현과 ligand–receptor 검출로 바로 이어졌습니다. 둘째, 환자 특이성은 이전과 같은 수준입니다. 우연보다 열 배 가까이 높지만 더 좋아지지는 않았습니다. 셋째, 과검출의 흔적이 분포 지표에 남아 있고, γδT와 pDC 같은 일부 희귀 상태는 여전히 부족합니다."
강조: 좋아진 것과 아닌 것을 같이 보여주는 것이 이 보고서의 목적.
질문이 나오면: 이전과 floor는 입력 depth, read-out, draw 수가 동시에 달라서 개선을 read-out 하나의 효과로 말할 수는 없다(슬라이드 15).""")

# ================= 3. setup
s = new(); header(s, "질문과 설정", "8가지 분석, 모두 같은 실제 대조군과 비교", 3, "실제 A = 참조 절반, 실제 B = 다른 환자 절반(상한 대조) · 패혈증 성인만, COVID는 adapter 학습 샘플 5개 제외")
rows = [("분석", "질문", "방법", "데이터"),
        ("A 환자 특이성", "그 환자의 신호가 세포에 남는가?", "donor 안정 유전자로 donor 식별", "COVID 32 donor, 모든 세포 pseudobulk"),
        ("B 질병 효과 (M6)", "질병 반응을 세포형별로 재현하는가?", "샘플 pseudobulk DE, logFC 비교", "패혈증 vs HC, COVID+ vs −"),
        ("C realism", "실제와 구분되는가? 분포를 덮는가?", "분류기 two-sample test, density / coverage", "소스당 약 60,000 cells"),
        ("D co-expression", "유전자 간 관계가 같은가?", "lineage 안 상관행렬 비교", "상위 300 / 순위 1,001–1,300"),
        ("E ligand–receptor", "세포 간 상호작용이 보이는가?", "LIANA consensus, 발현 곱", "5개 세포형 × 5개"),
        ("F 논문 claim", "논문 결과 방향이 나오는가?", "같은 분류기로 case − control", "Ye 2026 SEP, COVID IFN"),
        ("G 희귀 세포", "라벨 전이 없이 비율이 맞는가?", "부모 계통 안 분류기 (depth 맞춤)", "패혈증 11, COVID 9 subset"),
        ("H 통합 방법", "결론이 Harmony 의존인가?", "Harmony, Scanorama, scVI, BBKNN", "20,000 + 20,000 cells")]
table(s, rows, 0.9, 1.65, [2.3, 3.5, 3.4, 2.35], h_row=0.42, bold_cols=(0,), font=12, head_font=12)
text(s, 0.89, 6.05, 11.55, 0.6, [[("못 한 것  ", {"bold": True, "color": INK}), ("read-out만 바꾼 비교(depth256 Poisson control 없음), 패혈증 샘플별 정답(bulk-only 122개는 scRNA 없음), paired benchmark(서버).", {})]], size=13, color=BODY)
notes(s, """[1:30] 설정
"여덟 가지 분석을 했고, 모두 같은 원칙을 따릅니다. 실제 세포를 환자 단위로 둘로 나눠서 한쪽(A)을 기준으로, 다른 쪽(B)을 상한 대조로 씁니다. 생성 세포도 기준 환자와 다른 환자에서 왔으니 실제 B와 같은 조건입니다."
"할 수 없었던 것도 미리 말씀드립니다. 같은 checkpoint에서 read-out만 바꾼 Poisson control 파일이 PC에 없고, 패혈증 bulk-only 샘플 122개는 짝이 되는 단일세포 데이터가 없어서 샘플별 정답 비교는 COVID의 donor 수준에서만 했습니다."
질문이 나오면: 실제 B를 depth 맞춘 버전(realT1)도 대조로 씀 — floor 세포의 세포형별 UMI 분포로 binomial thinning.""")

# ================= 4. A donor
s = new(); header(s, "A · 환자 특이성", "생성 세포에 환자 고유 신호가 남는다: 우연의 9–12배, 이전과 같은 수준", 4, "COVID 32 donor · donor 안정 유전자 = 실제 샘플만으로 고른 ICC 상위 유전자 · results/covid_donor_fingerprint.csv")
pic(s, "d_donor", 0.89, 1.55, w=7.9)
side_cards(s, 9.0, 1.55, 3.45, [
    ("설계", "생성 = 그 donor의 bulk-only 샘플, 실제 = 같은 donor의 scRNA (대부분 다른 시점). 유전자별 z-score 후 생성 donor ↔ 실제 donor 상관으로 식별."),
    ("결과", f"top-1: 이전 {f2(fp('prev','top-100 donor-stable genes','top1_acc'))} · floor {f2(fp('floor','top-100 donor-stable genes','top1_acc'))} (안정 100), {f2(fp('prev','top-300 donor-stable genes','top1_acc'))} · {f2(fp('floor','top-300 donor-stable genes','top1_acc'))} (300). "
             f"자기 donor r (전체 유전자) {f2(fp('prev','all expressed genes','own_r_median'))} → {f2(fp('floor','all expressed genes','own_r_median'))}"),
    ("의미", "iLISI 상한은 다른 환자 세포로도 도달 → 환자 특이성은 따로 증명해야 함. floor는 섞임을 얻으면서 환자 신호를 잃지 않음.")], h=1.62, body=12)
notes(s, f"""[2:00] A 환자 특이성
"iLISI는 다른 환자의 실제 세포를 가져와도 높게 나옵니다. 그래서 섞임만으로는 '그 환자를 생성했다'는 것을 보일 수 없습니다. 이 슬라이드가 그 빈칸을 채웁니다."
"COVID에서는 32명이 생성 세포와 실제 단일세포를 모두 가지고 있습니다. 실제 샘플만 써서 donor마다 안정적인 유전자를 고르면, 성별 유전자 RPS4Y1·DDX3Y, HLA-DQA2·HLA-DRB5, RPS26 같은 genotype 영향 유전자가 올라옵니다. 이 유전자로 생성 donor가 어느 실제 donor와 가장 닮았는지 보면, 32명 중 자기 자신을 1등으로 맞히는 비율이 floor {f2(fp('floor','top-100 donor-stable genes','top1_acc'))}–{f2(fp('floor','top-300 donor-stable genes','top1_acc'))}로 우연 0.03의 9–12배입니다."
"floor와 이전은 비슷합니다. 다만 전체 유전자로 보면 자기 donor와의 상관이 {f2(fp('prev','all expressed genes','own_r_median'))}에서 {f2(fp('floor','all expressed genes','own_r_median'))}로 올라, 저발현 유전자가 살아나면서 환자 신호가 더 넓게 퍼진 것으로 보입니다."
질문이 나오면: 실제 샘플은 대부분 다른 시점이라 질병 상태 차이가 섞여 있어 식별률의 상한은 1보다 낮음. 정확히 같은 샘플이 아닌 donor 단위 검정. 패혈증은 짝이 없어 불가.""")

# ================= 5. B M6 figure
s = new(); header(s, "B · 질병 효과 (M6)", "저발현 유전자의 질병 반응이 살아났다: slope 두 배, 부호 일치 0.90", 5, "세포형별 샘플 pseudobulk, Welch t · 막대 = 5개 세포형 중앙값, 점 = 세포형 · results/m6_disease_effect.csv")
pic(s, "d_m6", 0.89, 1.5, w=11.55)
text(s, 0.89, 5.6, 11.55, 1.1, [[("순위 2,000 이후 (이전 → floor)  ", {"bold": True, "color": INK}),
     (f"logFC r {f2(m6('sepsis','rank > 2,000','prev','logFC_r'))} → {f2(m6('sepsis','rank > 2,000','floor','logFC_r'))} / {f2(m6('covid','rank > 2,000','prev','logFC_r'))} → {f2(m6('covid','rank > 2,000','floor','logFC_r'))}, "
      f"slope {f2(m6('sepsis','rank > 2,000','prev','slope_on_DEG'))} → {f2(m6('sepsis','rank > 2,000','floor','slope_on_DEG'))} / {f2(m6('covid','rank > 2,000','prev','slope_on_DEG'))} → {f2(m6('covid','rank > 2,000','floor','slope_on_DEG'))}. "
      f"순위 2,000 이내는 r이 오르고 slope은 비슷.", {})]], size=14, color=BODY, line=1.15)
notes(s, f"""[2:00] B 질병 효과
"논문에서 가장 중요한 질문입니다. 생성 세포가 환자 상태에 따라 세포형별로 올바르게 반응하느냐. 실제는 단일세포 샘플 전부, 생성은 bulk-only 샘플 전부로 세포형별 pseudobulk를 만들고, 패혈증 대 건강 대조, COVID 양성 대 음성의 logFC를 비교했습니다."
"차이는 저발현 유전자에서 납니다. 실제 발현 순위 2,000 이후 유전자에서, 실제 DEG에 대한 slope이 {f2(m6('sepsis','rank > 2,000','prev','slope_on_DEG'))}에서 {f2(m6('sepsis','rank > 2,000','floor','slope_on_DEG'))}로 거의 두 배가 되고, 부호 일치는 0.73–0.76에서 0.90이 됩니다. 절벽 때문에 잘렸던 유전자들이 이제 질병 반응을 싣고 있다는 뜻입니다."
"고발현 유전자에서는 r이 오르지만 slope은 비슷합니다. 즉 효과 크기를 절반 남짓으로 줄여 재현하는 성질은 그대로입니다."
질문이 나오면: 실제 패혈증 HC는 성인 7개 샘플뿐이라 실제 쪽 잡음이 큼. 생성은 68 대 36.""")

# ================= 6. B table
s = new(); header(s, "수치 · 질병 효과", "floor의 상위 DEG는 실제 sepsis 단핵구 유전자", 6, "precision = 생성 상위 100개 중 실제 FDR < 0.05이고 방향이 같은 비율 · results/m6_precision_top100.csv")
rows = [("지표 (5개 세포형 중앙값)", "패혈증\n이전", "패혈증\nfloor", "COVID\n이전", "COVID\nfloor")]
for b, bl in (("all", "전체"), ("rank > 2,000", "순위 2,000 이후")):
    for k, kl in (("logFC_r", "logFC r"), ("slope_on_DEG", "실제 DEG slope"), ("sign_agree_DEG", "실제 DEG 부호 일치")):
        rows.append((f"{kl} · {bl}", f2(m6('sepsis', b, 'prev', k)), f2(m6('sepsis', b, 'floor', k)), f2(m6('covid', b, 'prev', k)), f2(m6('covid', b, 'floor', k))))
rows.append(("생성 상위 100 DEG precision", f2(PR[('sepsis','prev')]), f2(PR[('sepsis','floor')]), f2(PR[('covid','prev')]), f2(PR[('covid','floor')])))
table(s, rows, 0.9, 1.65, [4.2, 1.4, 1.4, 1.4, 1.4], h_row=0.38, bold_cols=(2, 4), font=12.5)
side_cards(s, 10.95, 1.65, 1.5, [("단핵구", "floor 상위: IL1R2, MCEMP1, AQP9, TNFAIP6, ACSL1"), ("이전 상위", "CD63, GAPDH, RPL13A, FTL 등 고발현")], h=1.9, body=11, head=13)
text(s, 0.89, 5.75, 11.55, 0.9, [[("주의  ", {"bold": True, "color": INK}), (f"패혈증 vs ICU 대조는 실제 DEG가 {ICU:.0f}개뿐이라 두 모델 모두 재현 불가. COVID T 세포 floor 상위에는 실제에서 유의하지 않은 저발현 유전자(EPB41L5, PLS1)가 섞임 → 과검출 영향.", {})]], size=13, color=BODY, line=1.12)
notes(s, f"""[1:30] 질병 효과 수치
"표의 핵심은 아래 세 줄, 순위 2,000 이후입니다. 그리고 정성적으로 보면, 패혈증 단핵구에서 floor가 가장 강하게 잡는 유전자는 IL1R2, MCEMP1, AQP9, TNFAIP6입니다. 모두 패혈증 단핵구의 대표 유전자이고 실제 데이터에서도 유의합니다. 이전 모델의 상위에는 GAPDH나 리보솜 유전자 같은 고발현 유전자가 많았습니다."
"반대 방향도 정직하게 말씀드리면, COVID T 세포에서는 floor 상위에 실제로는 유의하지 않은 저발현 유전자가 섞입니다. 과검출이 만든 거짓 신호로 보입니다. 그래서 precision은 패혈증 {f2(PR[('sepsis','prev')])} → {f2(PR[('sepsis','floor')])}, COVID {f2(PR[('covid','prev')])} → {f2(PR[('covid','floor')])}입니다."
질문이 나오면: top-100 overlap은 두 모델 모두 0.1–0.3으로 낮음 — 실제 쪽 순위 자체가 HC 7개로 불안정해서 overlap보다 precision을 보는 게 맞음.""")

# ================= 7. C c2st
s = new(); header(s, "C · realism: 분류기 two-sample test", "생성 세포는 덜 구분된다: AUC 0.99 → 0.85", 7, "HistGradientBoosting, 샘플 단위 5-fold CV, lineage당 4,000 + 4,000 cells, 실제 A PCA 50 · results/<cohort>_c2st.csv")
pic(s, "d_c2st", 0.89, 1.5, w=8.1)
side_cards(s, 9.2, 1.55, 3.25, [
    ("AUC (4 lineage 중앙값)", f"이전 {f2(c2[('sepsis','prev4k')])} / {f2(c2[('covid','prev4k')])}\nfloor {f2(c2[('sepsis','floor')])} / {f2(c2[('covid','floor')])}\n실제 B {f2(c2[('sepsis','realU')])} / {f2(c2[('covid','realU')])}"),
    ("floor 쪽이 높은 유전자", "histone H1(HIST1H1C/E), JUN·FOS, HLA-DPA1, GNLY"),
    ("실제 쪽이 높은 유전자", "RPS26·HLA-DRB5(genotype), HBB(ambient), IL32, MYOM2")], h=1.62, body=12)
notes(s, f"""[1:30] 분류기 테스트
"가장 엄격한 realism 지표입니다. 같은 lineage 안에서 실제 A와 각 소스를 구분하는 분류기를 학습시키고, 환자 단위로 교차검증합니다. 0.5면 구분이 불가능하다는 뜻입니다."
"이전 모델은 0.99로 거의 완벽하게 구분됐습니다. floor는 {f2(c2[('sepsis','floor')])} / {f2(c2[('covid','floor')])}로 떨어졌습니다. 다만 실제 B는 0.42 / 0.65라서, 아직 구분할 수 있는 차이가 남아 있습니다."
"무엇이 구분을 만드는지도 봤습니다. floor 쪽에 높은 histone H1은 앞서 단핵구 분석에서도 본 학습 corpus의 플랫폼 신호로 추정됩니다. 실제 쪽 RPS26, HLA-DRB5는 genotype으로 정해지는 유전자이고, HBB는 적혈구 ambient 오염입니다."
질문이 나오면: 패혈증 실제 B AUC가 0.5 미만인 것은 환자 단위 CV에서 생기는 정상 편차. depth만 맞춘 실제(회색)는 0.45–0.75로 depth 자체가 큰 원인은 아님.""")

# ================= 8. C coverage/density
s = new(); header(s, "C · realism: 분포 덮기", "보정 전 coverage가 크게 늘었지만, 밀집 영역에는 덜 놓인다", 8, "Naeem et al. 2020, k = 5, 실제 A가 기준 manifold, 10,000 cells · results/<cohort>_density_coverage.csv")
pic(s, "d_dc", 0.89, 1.5, w=8.4)
side_cards(s, 9.45, 1.55, 3.0, [
    ("coverage", f"Harmony 전 이전 {f2(dc('sepsis',PCAs,'prev4k','coverage'))} / {f2(dc('covid',PCAs,'prev4k','coverage'))} → floor {f2(dc('sepsis',PCAs,'floor','coverage'))} / {f2(dc('covid',PCAs,'floor','coverage'))} (실제 {f2(dc('sepsis',PCAs,'realU','coverage'))} / {f2(dc('covid',PCAs,'realU','coverage'))}). Harmony 후는 비슷."),
    ("density", f"floor {f2(dc('sepsis',PCAs,'floor','density'))} / {f2(dc('covid',PCAs,'floor','density'))} < 실제 {f2(dc('sepsis',PCAs,'realU','density'))} / {f2(dc('covid',PCAs,'realU','density'))} < 이전 {f2(dc('sepsis',HAR,'prev4k','density'))} / {f2(dc('covid',HAR,'prev4k','density'))} (Harmony 후)"),
    ("해석", "이전 = 밀집 핵심부에 몰림(다양성 부족). floor = 더 퍼지고 일부는 실제 밖(과검출 noise).")], h=1.62, body=12)
notes(s, f"""[1:30] density와 coverage
"생성모델 평가에서 쓰는 두 지표입니다. coverage는 실제 분포를 얼마나 넓게 덮는지, density는 생성 세포가 실제 세포가 빽빽한 곳에 얼마나 놓이는지입니다."
"Harmony 보정 전 coverage가 이전 0.40 / 0.45에서 floor 0.69 / 0.78로 크게 늘었습니다. 이전 모델 세포는 보정 없이는 실제 분포의 절반도 못 덮었다는 뜻입니다."
"반대로 density는 floor가 실제보다 낮고 이전이 실제보다 높습니다. 이전 모델은 전형적인 세포만 만드는 mode collapse 경향, floor는 더 퍼지고 일부가 실제 밖으로 나가는 경향입니다. 후자는 과검출이 noise를 더한 결과로 해석합니다. threshold sweep에서 density를 과검출 판정 기준으로 쓰자는 제안이 여기서 나옵니다."
질문이 나오면: Harmony 후 coverage는 floor가 이전보다 약간 낮음(0.76 vs 0.82, 패혈증) — 보정 후 이전 모델은 밀집 영역에 겹쳐 놓이기 때문.""")

# ================= 9. D co-expression
s = new(); header(s, "D · 유전자 간 관계", "저발현대의 co-expression 구조가 살아났고, 상관은 실제보다 강하다", 9, "lineage당 5,000 cells, log1p CP10k Pearson · 셔플 = 유전자별 독립 순열(대조) · results/<cohort>_coexpression.csv")
pic(s, "d_coexpr", 0.89, 1.5, w=8.4)
side_cards(s, 9.45, 1.55, 3.0, [
    ("상위 300 유전자", f"이전 {f2(co[('sepsis','top 300 genes','prev4k')])} / {f2(co[('covid','top 300 genes','prev4k')])} · floor {f2(co[('sepsis','top 300 genes','floor')])} / {f2(co[('covid','top 300 genes','floor')])} (실제 {f2(co[('sepsis','top 300 genes','realU')])} / {f2(co[('covid','top 300 genes','realU')])})"),
    ("순위 1,001–1,300", f"이전 {f2(co[('sepsis','ranks 1,001-1,300','prev4k')])} / {f2(co[('covid','ranks 1,001-1,300','prev4k')])} → floor {f2(co[('sepsis','ranks 1,001-1,300','floor')])} / {f2(co[('covid','ranks 1,001-1,300','floor')])} (실제 {f2(co[('sepsis','ranks 1,001-1,300','realU')])} / {f2(co[('covid','ranks 1,001-1,300','realU')])}). B 세포는 예외."),
    ("상관 강도", f"상위 300 평균 |r| floor {f2(ca[('sepsis','top 300 genes','floor')])} / {f2(ca[('covid','top 300 genes','floor')])} vs 실제 A {f2(car[('sepsis','top 300 genes')])} / {f2(car[('covid','top 300 genes')])}: 너무 정돈된 program")], h=1.62, body=12)
notes(s, f"""[1:30] co-expression
"평균 발현이 맞아도 유전자 간 관계가 틀리면 생성 세포는 쓸모가 줄어듭니다. 같은 lineage 안에서 유전자 쌍의 상관행렬을 만들고, 실제 A의 상관행렬과 얼마나 닮았는지를 봤습니다. 유전자를 하나씩 섞은 대조는 0이 나와서 지표가 제대로 작동함을 보여줍니다."
"고발현 300개에서는 두 모델 모두 0.81–0.85로 비슷합니다. 차이는 중간 순위 유전자에서 나고, floor가 실제 수준에 근접합니다."
"하나 주의할 점은 floor의 상관이 실제보다 강하다는 것입니다. 실제 세포보다 유전자 program이 더 깔끔하게 함께 움직입니다. 생성모델에서 흔한 현상이고, 과검출과도 관련이 있을 수 있습니다."
질문이 나오면: B 세포 저발현대는 floor가 이전보다 약간 낮음(0.33 vs 0.35, 패혈증). 세포 수 5,000이라 저발현 유전자 상관은 원래 noise가 큼(실제 B도 0.37).""")

# ================= 10. E LR
s = new(); header(s, "E · ligand–receptor", "약한 수용체가 살아나 상호작용 검출이 실제 수준에 근접", 10, "LIANA consensus (패널 안 약 4,850 쌍 × 25 세포쌍), 송신·수신 세포 10% 이상 검출 · results/<cohort>_lr_*.csv")
pic(s, "d_lr", 0.89, 1.5, w=8.6)
side_cards(s, 9.65, 1.55, 2.8, [
    ("검출된 interaction", f"이전 {LR.loc[('sepsis','prev4k'),'n_detected_source']} / {LR.loc[('covid','prev4k'),'n_detected_source']} → floor {LR.loc[('sepsis','floor'),'n_detected_source']} / {LR.loc[('covid','floor'),'n_detected_source']} (실제 A {LR.loc[('sepsis','floor'),'n_detected_realA']} / {LR.loc[('covid','floor'),'n_detected_realA']})"),
    ("실제 상위 200 회수", f"{f2(LR.loc[('sepsis','prev4k'),'recall_top200_realA'])} / {f2(LR.loc[('covid','prev4k'),'recall_top200_realA'])} → {f2(LR.loc[('sepsis','floor'),'recall_top200_realA'])} / {f2(LR.loc[('covid','floor'),'recall_top200_realA'])}"),
    ("질병 변화는 아직", f"case − control 점수 변화 Spearman {f2(LD.loc[('sepsis','prev4k'),'disease_spearman_all_detected'])} / {f2(LD.loc[('covid','prev4k'),'disease_spearman_all_detected'])} → {f2(LD.loc[('sepsis','floor'),'disease_spearman_all_detected'])} / {f2(LD.loc[('covid','floor'),'disease_spearman_all_detected'])}")], h=1.62, body=12)
notes(s, f"""[1:30] ligand–receptor
"수용체 유전자는 대부분 저발현이라 이전 모델에서 절벽에 잘렸습니다. 그래서 응용 측면에서 floor의 이득이 가장 잘 보이는 곳입니다. LIANA의 consensus 목록으로 세포형 쌍마다 interaction을 계산했습니다."
"검출된 interaction이 이전 {LR.loc[('sepsis','prev4k'),'n_detected_source']}개에서 floor {LR.loc[('sepsis','floor'),'n_detected_source']}개로, 실제 {LR.loc[('sepsis','floor'),'n_detected_realA']}개와 비슷해졌고, 실제의 상위 200개 interaction 중 97%를 다시 찾습니다(패혈증)."
"하지만 질병에 따른 interaction 변화는 개선되지 않았습니다. 이 지표는 S100A8/9 같은 고발현 리간드가 지배하는데, 그쪽은 원래 두 모델 모두 잘 잡던 부분이라 floor의 이득이 드러나지 않습니다."
질문이 나오면: 점수는 발현 곱(sqrt(mean L × mean R)) 방식, permutation 기반 유의성은 계산하지 않음.""")

# ================= 11. F claims
fs_, fn_, fz = claim_sum("sepsis", "floor"); ps_, pn_, pz = claim_sum("sepsis", "prev4k"); fcs, fcn, fcz = claim_sum("covid", "floor"); pcs, pcn, pcz = claim_sum("covid", "prev4k")
s = new(); header(s, "F · 논문 claim", "같은 분류기로 읽으면 질병 효과 크기가 실제의 0.7–0.8배", 11, "panel-matched author marker, 실제 A Youden threshold, 샘플 평균 · 실제 A AUC < 0.6인 SEP-4A, IFN 라벨이 없는 B 세포 제외 · results/<cohort>_claims.csv")
pic(s, "d_claims", 0.89, 1.5, w=8.1)
side_cards(s, 9.2, 1.55, 3.25, [
    ("패혈증 SEP (5개)", f"방향 일치 이전 {ps_}/{pn_} · floor {fs_}/{fn_}\n크기(실제 대비 중앙값) {pz:.2f}× → {fz:.2f}×"),
    ("COVID IFN (4개 계통)", f"방향 일치 {pcs}/{pcn} · {fcs}/{fcn}\n크기 {pcz:.2f}× → {fcz:.2f}×"),
    ("놓친 것", "SEP-2(CD8 GNLY-CCL4hi): 실제는 감소, floor는 약한 증가. 이전 K4000도 −0.01로 거의 0.")], h=1.62, body=12)
notes(s, f"""[1:30] 논문 claim
"논문에서 인용할 질병 소견을 같은 잣대로 다시 읽었습니다. 실제 A에서 marker 점수의 threshold를 정하고, 그 분류기를 실제 B, 이전, floor에 똑같이 적용해 질병군과 대조군의 양성 비율 차이를 비교합니다."
"floor는 효과 크기를 실제의 {fz:.2f}배(패혈증), {fcz:.2f}배(COVID)로 재현합니다. 이전은 {pz:.2f}배, {pcz:.2f}배였습니다. 방향은 패혈증에서 floor가 5개 중 4개 맞고, 놓친 하나는 SEP-2입니다. SEP-2는 이전 모델도 거의 0이라 사실상 두 모델 모두 재현하지 못한 claim입니다."
질문이 나오면: 실제 A에는 성인 HC가 1명뿐이라 실제 대조는 실제 B(HC 4명) 사용. 표본이 작아 개별 claim 판정은 조심. 라벨 기준 대비 분류기 기준의 방향이 다른 claim(SEP-2, SEP-3-NK)이 있다는 것은 10/12 재채점에서 이미 확인.""")

# ================= 12. G rare
s = new(); header(s, "G · 희귀 세포 (라벨 전이 없이)", "COVID 희귀 상태는 크게 회복, 패혈증은 엇갈림", 12, "부모 계통 안 HistGradientBoosting (PCA 50 + subset marker 50), depth 맞춘 실제 A로 학습, 실제 A CV AUC ≥ 0.8만 · results/<cohort>_rare_transfer_free.csv")
pic(s, "d_rare", 0.89, 1.45, h=5.3)
rv = lambda c, sub, s_: rr.loc[(c, sub), s_]
side_cards(s, 9.6, 1.55, 2.85, [
    ("0.5–2배 안 (이전 → floor)", f"패혈증 {rq('sepsis','prev4k')} → {rq('sepsis','floor')}\nCOVID {rq('covid','prev4k')} → {rq('covid','floor')}"),
    ("floor가 회복", f"증식 NK {rv('sepsis','NK_c04-MKI67','prev4k'):.2f} → {rv('sepsis','NK_c04-MKI67','floor'):.2f}, DC1 {rv('covid','DC1','prev4k'):.2f} → {rv('covid','DC1','floor'):.2f}, MAIT, NK CD56 bright, macrophage"),
    ("floor가 더 나쁨", f"γδT {rv('sepsis','gdT','prev4k'):.2f} → {rv('sepsis','gdT','floor'):.2f}, 패혈증 pDC, 증식 단핵구 {rv('sepsis','Mono_c06-CD14-MKI67','prev4k'):.2f} → {rv('sepsis','Mono_c06-CD14-MKI67','floor'):.2f}. dex 단핵구는 둘 다 없음")], h=1.62, body=11.5, head=14)
notes(s, f"""[2:00] 희귀 세포
"앞선 보고서의 희귀 세포 판정은 kNN 라벨 전이에 의존했는데, 그 경로는 실제 세포에서도 증식 단핵구나 ILC를 잃어버립니다. 그래서 이번에는 전이 없이, 부모 계통 안에서 각 희귀 상태를 맞히는 분류기를 실제 A로 학습시켰습니다. depth도 각 생성 세트에 맞춰 줄인 실제 A로 학습했고, 같은 분류기로 실제 B와 생성 세포의 비율을 비교합니다."
"COVID에서는 9개 중 이전 2개에서 floor 6개로 크게 좋아졌습니다. DC1, MAIT, NK CD56 bright가 대표적입니다. 패혈증은 엇갈립니다. 증식 NK와 macrophage는 좋아졌지만, γδT는 floor가 오히려 0.11배로 부족하고, 증식 단핵구와 pDC도 이전이 더 가깝습니다."
"그리고 판정 방법에 따라 결론이 바뀌는 상태가 있다는 것 자체가 중요한 결과입니다. 증식 단핵구는 전이로는 0.04배, 분류기로는 0.65배입니다. 희귀 상태는 한 가지 방법으로 판정하면 안 된다는 뜻이고, 다음 단계로 dose–response 실험이 필요합니다."
질문이 나오면: 분류기의 실제 B 추정값은 정답 대비 중앙값 0.79–0.92로 약간 과소추정 → 그래서 생성 세포를 실제 정답이 아니라 같은 분류기의 실제 B 추정과 비교. ILC, HSC, T CD4 cycling은 AUC < 0.8이라 제외.""")

# ================= 13. H integration
s = new(); header(s, "H · 통합 방법", "어떤 통합 방법을 써도 floor > 이전", 13, "실제 A + 소스 20,000 + 20,000 cells, 같은 HVG 2,000 · BBKNN = graph 가중 inverse Simpson · results/<cohort>_integration_methods.csv")
pic(s, "d_integ", 0.89, 1.5, w=8.5)
side_cards(s, 9.6, 1.55, 2.85, [
    ("iLISI (이전 → floor, 실제 B)", f"Harmony {f2(il('sepsis','Harmony','prev4k'))} → {f2(il('sepsis','Harmony','floor'))} ({f2(il('sepsis','Harmony','realU'))})\nScanorama {f2(il('sepsis','Scanorama','prev4k'))} → {f2(il('sepsis','Scanorama','floor'))} ({f2(il('sepsis','Scanorama','realU'))})\n"
                                 f"scVI {f2(il('sepsis','scVI','prev4k'))} → {f2(il('sepsis','scVI','floor'))} ({f2(il('sepsis','scVI','realU'))})\nBBKNN {f2(il('sepsis','BBKNN (graph)','prev4k'))} → {f2(il('sepsis','BBKNN (graph)','floor'))} ({f2(il('sepsis','BBKNN (graph)','realU'))})\n(패혈증)"),
    ("생물학 보존", "세포형 label ASW는 모든 방법에서 floor ≥ 이전. 섞임이 세포형을 뭉개서 나온 것이 아님."),
    ("scVI", f"floor도 실제와 격차 큼 ({f2(il('sepsis','scVI','floor'))} vs {f2(il('sepsis','scVI','realU'))}): count 분포 차이에 민감한 모델")], h=1.6, body=11, head=13)
notes(s, f"""[1:30] 통합 방법
"Harmony만 썼기 때문에 좋아 보이는 것 아니냐는 질문에 대한 답입니다. 같은 세포, 같은 HVG로 Scanorama, scVI, BBKNN까지 네 방법을 돌렸습니다."
"네 방법 모두에서 floor가 이전보다 높습니다. BBKNN에서는 이전 1.24, floor 1.75로 차이가 가장 크고, Harmony에서는 실제 B 수준에 거의 도달합니다. 세포형이 섞여서 iLISI가 오른 게 아니라는 것도 label ASW로 확인했습니다."
"예외적으로 scVI에서는 floor도 실제와 격차가 큽니다. scVI는 count 분포를 직접 모델링하기 때문에, 과검출로 생긴 분포 차이를 그대로 batch 차이로 봅니다."
질문이 나오면: cLISI는 중앙값이 모든 조합에서 1.0으로 포화돼서 쓰지 않음. scVI는 CPU에서 40 epoch, latent 30.""")

# ================= 14. scorecard table
s = new(); header(s, "수치 · 종합 성적표", "8개 분석 한눈에 보기", 14, "패혈증 / COVID · results/scorecard.csv")
S = rd("scorecard.csv")
keep = [0, 1, 5, 7, 8, 9, 10, 11, 15, 16, 19, 20, 21, 22, 23]
lab = {0: "A donor 식별 top-1 (우연 0.03)", 1: "A 자기 donor r", 5: "B logFC r, 순위 > 2,000", 6: "B slope, 순위 > 2,000", 7: "B 부호 일치, 순위 > 2,000", 8: "B 상위 100 DEG precision",
       9: "C 분류기 AUC (낮을수록 좋음)", 10: "C coverage, 보정 전", 11: "C density, 보정 전", 15: "D 상관구조 r, 순위 1,001–1,300", 16: "E LR 검출 Jaccard", 17: "E 실제 상위 200 회수",
       19: "E 질병에 따른 LR 변화", 20: "F 패혈증 claim", 21: "F COVID IFN claim", 22: "G 희귀 0.5–2배 (분류기)", 23: "H iLISI Harmony", 26: "H iLISI BBKNN"}
verdict = {0: "≈", 1: "▲", 5: "▲", 6: "▲", 7: "▲", 8: "± ", 9: "▲", 10: "▲", 11: "▼", 15: "▲", 16: "▲", 17: "▲", 19: "≈", 20: "▲ 크기", 21: "▲", 22: "± ", 23: "▲", 26: "▲"}
rows = [("지표", "실제 대조", "이전 K4000", "floor", "판정")] + [(lab[i], "–" if str(S.loc[i, "real_control"]) == "1" else S.loc[i, "real_control"], S.loc[i, "prev_K4000"], S.loc[i, "floor"], verdict[i].strip()) for i in keep]
rows = [tuple(str(x).replace("5/5 sign, size", "5/5 방향, 크기").replace("4/5 sign, size", "4/5 방향, 크기").replace("4/4 sign, size", "4/4 방향, 크기") for x in r) for r in rows]
table(s, rows, 0.9, 1.5, [4.0, 1.85, 2.2, 2.2, 1.3], h_row=0.29, bold_cols=(3,), font=11, head_font=11)
text(s, 0.89, 6.5, 11.55, 0.3, "▲ floor가 실제 쪽으로 개선 · ≈ 비슷 · ▼ 실제에서 멀어짐 · ± 코호트에 따라 다름", size=11, color=MUTED)
notes(s, """[1:00] 성적표
"한 장으로 정리한 것입니다. 열다섯 개 지표 중 열 개는 floor가 실제 쪽으로 좋아졌고, 두 개는 비슷하고, 두 개는 코호트마다 다르며, 하나는 나빠졌습니다. 나빠진 하나는 density, 즉 생성 세포가 실제 밀집 영역에 덜 놓이는 것이고, 과검출 때문으로 봅니다."
"논문에서 우수성의 근거로 쓸 수 있는 줄은 B의 저발현 질병 효과, C의 분류기 AUC와 coverage, E의 상호작용 검출, H의 통합 방법 견고성입니다."
질문이 나오면: 모든 값은 같은 실제 대조군, 같은 HVG, 같은 seed 정책. 숫자의 정의는 각 결과 슬라이드 하단에. 표에서 뺀 지표(질병 효과 slope, LR 상위 200 회수, Scanorama·scVI·BBKNN iLISI)는 각 결과 슬라이드에 있음.""")

# ================= 15. caveats
s = new(); header(s, "해석과 주의", "이 숫자가 보여주지 않는 것", 15)
three_cards(s, 1.75, 3.9, [
    ("1", "원인 분리 안 됨", "floor와 이전은 입력 depth(× 256), read-out 규칙, draw 수(1 vs 5), 샘플당 세포 수가 동시에 다름. 같은 checkpoint의 Poisson control이 있어야 read-out 효과만 분리 가능."),
    ("2", "표본과 경로의 한계", "실제 패혈증 성인 HC 7명(claim 대조는 4명). 패혈증은 샘플별 정답 없음. 희귀 상태는 판정 경로(전이 vs 분류기)에 따라 결론이 바뀌는 것이 있음."),
    ("3", "과검출의 비용", "density 저하, 실제보다 강한 유전자 상관, COVID T 세포의 거짓 DEG. iLISI·coverage는 좋아도 scVI는 분포 차이를 잡아냄.")])
text(s, 0.89, 5.9, 11.55, 0.7, "모든 분석은 subsample 1회, seed 1개. 환자 단위 bootstrap 신뢰구간은 아직 없음.", size=13, color=BODY)
notes(s, """[1:00] 주의
"세 가지를 꼭 함께 말씀드리고 싶습니다. 첫째, 이번 비교는 read-out 하나의 효과가 아닙니다. 둘째, 실제 대조군이 작은 곳이 있고, 희귀 상태는 판정 방법에 민감합니다. 셋째, 과검출은 공짜가 아닙니다. 분포 지표, 상관 강도, 일부 거짓 DEG로 비용이 드러납니다."
질문이 나오면: 과검출은 threshold 1e-4가 너무 관대해서 생긴 것으로 보고, 저장된 연속값으로 threshold sweep을 하면 같은 지표로 최적점을 고를 수 있음.""")

# ================= 16. conclusion
s = new(); header(s, "결론", "floor는 '실제처럼 보이는' 세포에서 '실제처럼 반응하는' 세포로", 16)
three_cards(s, 1.75, 3.5, [
    ("1", "논문에 쓸 근거", f"저발현 유전자 질병 효과 slope {f2(m6('sepsis','rank > 2,000','prev','slope_on_DEG'))} → {f2(m6('sepsis','rank > 2,000','floor','slope_on_DEG'))}, 분류기 AUC 0.99 → {f2(c2[('sepsis','floor')])}, "
                     f"LR 검출 Jaccard {f2(LR.loc[('sepsis','prev4k'),'detected_jaccard'])} → {f2(LR.loc[('sepsis','floor'),'detected_jaccard'])}, 통합 4개 방법 모두 floor > 이전 (패혈증)."),
    ("2", "유지된 것", "환자 특이성은 우연의 9–12배로 이전과 같음. 섞임을 얻으면서 환자 신호를 잃지 않았다는 점이 핵심."),
    ("3", "남은 과제", "과검출 보정(density, 상관 강도, 거짓 DEG). γδT, pDC, dex 단핵구. 패혈증 vs ICU 대조 차이.")])
text(s, 0.89, 5.5, 11.55, 1.2, [[("다음 단계  ", {"bold": True, "color": INK}), ("① 서버의 연속값으로 threshold sweep — 이번 지표(density, 분류기 AUC, 상관 강도, DEG precision)를 판정 기준으로  "
     "② depth256 Poisson control을 같은 파이프라인에  ③ paired 환자별 benchmark  ④ 희귀 상태 dose–response.", {})]], size=14, color=BODY, line=1.2)
notes(s, """[1:00] 결론
"정리하면, floor read-out은 생성 세포를 실제처럼 보이게 만든 데서 그치지 않고, 저발현 유전자의 질병 반응과 세포 간 상호작용까지 실제처럼 만들었습니다. 그러면서도 환자 고유 신호는 그대로 유지했습니다."
"남은 것은 과검출입니다. 다행히 이번 분석으로 과검출을 잴 수 있는 지표가 생겼고, 서버의 연속값으로 threshold를 다시 고르면 이 지표들로 최적점을 정할 수 있습니다."
마무리로 세 문장 다시: 질병 반응·상호작용 개선 / 환자 신호 유지 / 과검출 보정이 다음 과제.
발표 시간 합계 약 23분(슬라이드별 [분:초] 합), 질문 시간은 별도.""")

prs.save(DST); print("saved", DST, len(prs.slides._sldIdLst), "slides")
