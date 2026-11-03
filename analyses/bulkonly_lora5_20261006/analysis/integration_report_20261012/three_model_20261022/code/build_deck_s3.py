"""Korean deck + Korean speaker notes: PBISC floor read-out vs the earlier data sets (previous PBISC K4000, scSemiProfiler), strengths and weaknesses
(2026-10-22). Same template / layout / styles as the earlier floor decks. Every number is read from s3/results (and the big_joint summaries)."""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deck_helpers import *
R = "/mnt/user-data/outputs/s3/results"; TR = "/mnt/user-data/outputs/ilisi_floor_20261021/results/big"
DST = "/mnt/user-data/outputs/deck/PBISC floor vs 기존 데이터 세트 비교 보고서.pptx"
rd = lambda n: pd.read_csv(f"{R}/{n}"); f2 = lambda x: f"{x:.2f}"; pct = lambda x: f"{100*x:.0f}%"
CO_ = ("sepsis", "covid"); sv = lambda fn: f"{fn('sepsis')} / {fn('covid')}"
# ---------------- numbers
J = {(c, s): json.load(open(f"{TR}/{c}_{s}_summary.json")) for c in CO_ for s in ("realU", "floor", "prev4k", "scsemi")}
jj = lambda s, k: sv(lambda c: f2(J[(c, s)][k]))
I = pd.concat([rd(f"{c}_integration_methods.csv") for c in CO_]).drop_duplicates(["cohort", "pair", "method"], keep="last").set_index(["cohort", "method", "pair"])
il = lambda m, s: sv(lambda c: f2(I.loc[(c, m, s), "iLISI"]))
BI = rd("bootstrap_identity_mixing.csv"); PS = BI[BI.analysis.str.startswith("per-sample iLISI")].assign(c=lambda d: d.analysis.str.split().str[-1]).set_index(["c", "metric"])
ps = lambda s: sv(lambda c: f2(PS.loc[(c, s), "estimate"]))
FPb = BI[BI.analysis.str.startswith("COVID donor")].set_index("metric"); SI = BI[BI.analysis.str.startswith("scSemi sample")]
C2 = pd.concat([rd(f"{c}_c2st.csv") for c in CO_]); c2 = C2.groupby(["cohort", "source"]).c2st_auc.median(); cc2 = lambda s: sv(lambda c: f2(c2[(c, s)]))
DC = pd.concat([rd(f"{c}_density_coverage.csv") for c in CO_]); DC = DC[DC.lineage == "all"].set_index(["cohort", "space", "source"])
dc = lambda sp_, s, k: sv(lambda c: f2(DC.loc[(c, sp_, s), k])); PCAs, HAR = "PCA (no correction)", "Harmony"
CO = pd.concat([rd(f"{c}_coexpression.csv") for c in CO_]).groupby(["cohort", "genes", "source"]).corr_of_corr.median(); co = lambda b, s: sv(lambda c: f2(CO[(c, b, s)]))
LR = pd.concat([rd(f"{c}_lr_recovery.csv") for c in CO_]).set_index(["cohort", "source"]); LD = pd.concat([rd(f"{c}_lr_disease.csv") for c in CO_]).set_index(["cohort", "source"])
lr = lambda s, k: sv(lambda c: f2(LR.loc[(c, s), k])); ld = lambda s: sv(lambda c: f2(LD.loc[(c, s), "disease_spearman_all_detected"]))
M = rd("m6_three_models.csv"); M = M[M.contrast == "case vs control"].groupby(["cohort", "genes", "source"])[["logFC_r", "slope_on_DEG", "sign_agree_DEG"]].median()
m6 = lambda b, s, k: sv(lambda c: f2(M.loc[(c, b, s), k])); m6v = lambda c, b, s, k: M.loc[(c, b, s), k]
CI = rd("m6_bootstrap_ci.csv").set_index(["cohort", "genes", "metric", "comparison"])
def ci(c, b, m, comp):
    r = CI.loc[(c, b, m, comp)]; return f"{r.med:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}]"
PR = rd("m6_broadcast_partial.csv"); PR = PR[PR.contrast == "case vs control"].groupby(["cohort", "genes", "source"])[["r_broadcast_real", "partial_r_model_given_broadcast"]].median()
pr = lambda b, s: sv(lambda c: f2(PR.loc[(c, b, s), "partial_r_model_given_broadcast"])); rb = lambda b: sv(lambda c: f2(PR.loc[(c, b, "floor"), "r_broadcast_real"]))
CP = rd("composition_shift.csv"); CP = CP[CP.contrast == "case vs control"].set_index(["cohort", "level", "source"])
cp = lambda l, s, k: sv(lambda c: f2(CP.loc[(c, l, s), k]))
CC = rd("composition_bootstrap_ci.csv").set_index(["cohort", "level", "metric", "comparison"])
CL = pd.concat([rd(f"{c}_claims.csv") for c in CO_]); CL = CL[(CL.auc_realA >= 0.6) & (CL.claim != "ESKD-IFN B cell")]
W = CL.pivot_table(index=["cohort", "claim"], columns="source", values="delta_pos_frac")
def cl(c, s):
    w = W.loc[c]; m = np.nanmedian(w[s] / w['realU']); return f"{int((np.sign(w[s]) == np.sign(w['realU'])).sum())}/{len(w)}, " + ("~0×" if abs(m) < 0.005 else f"{m:.2f}×")
RR = pd.concat([rd(f"{c}_rare_transfer_free.csv") for c in CO_]); RR = RR[RR.cv_auc_realA >= 0.8].copy(); RR["ratio"] = RR.gen_estimate_pct / RR.realB_estimate_pct
RR["ok"] = RR.ratio.between(0.5, 2); rk = RR.groupby(["cohort", "model_depth"]).ok.agg(["sum", "size"])
rq = lambda s: sv(lambda c: f"{int(rk.loc[(c, s), 'sum'])}/{int(rk.loc[(c, s), 'size'])}"); rr = RR.pivot_table(index=["cohort", "subset"], columns="model_depth", values="ratio")
MY = pd.concat([rd(f"{c}_myeloid_identity.csv") for c in CO_]).set_index(["cohort", "source"]); my = lambda s, k, f=f2: sv(lambda c: f(MY.loc[(c, s), k]))
MG = rd("sepsis_myeloid_model_specific_genes.csv"); mg = MG[MG.model == "floor"].set_index("gene")
G = {}
import anndata as ad
for c in CO_:
    for s in ("realU", "prev4k", "floor", "scsemi"):
        a = ad.read_h5ad(f"/home/claude/work/s3/big/{c}_{s}.h5ad"); G[(c, s)] = float(np.median(np.diff(a.X.tocsr().indptr)))
gpc = lambda s: sv(lambda c: f"{G[(c, s)]:.0f}")
ns = {"sepsis": "패혈증", "covid": "COVID"}

# ================= 1. title
s = new(); box(s, 0, 0, 0.17, 7.5, BLUE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 1.92, 11.9, 0.3, "PBISC 생성 행렬 검증 · 데이터 세트 비교", size=14, color=BLUE, bold=True, spc=150)
ph = s.shapes.title; ph.left, ph.top, ph.width, ph.height = Inches(0.77), Inches(2.5), Inches(11.9), Inches(1.5)
tf = ph.text_frame; tf.clear(); tf.word_wrap = True
for a in ("margin_left", "margin_right", "margin_top", "margin_bottom"): setattr(tf, a, 0)
p = tf.paragraphs[0]; style_run(p.add_run(), 40, INK, True); p.runs[0].text = "새 floor PBISC 데이터 세트는 기존보다\n무엇이 낫고, 무엇이 아직 부족한가"
text(s, 0.89, 4.15, 11.9, 0.7, "floor read-out(draw0)을 이전 PBISC K = 4,000, scSemiProfiler, bulk 그대로 복사한 기준선과 같은 실제 대조군·같은 파이프라인으로 비교 "
     "(패혈증 Ye 2026, COVID Imperial ESKD)", size=18, color=BODY)
pill(s, 0.89, 5.1, 1.75, "PBISC floor", BLUE); pill(s, 2.8, 5.1, 2.1, "PBISC prev K4000", RED); pill(s, 5.06, 5.1, 1.95, "scSemiProfiler", "7B3294")
pill(s, 7.17, 5.1, 1.75, "bulk broadcast", "E08214"); pill(s, 9.08, 5.1, 1.55, "Real (held-out)", GREY)
text(s, 0.89, 6.82, 11.9, 0.26, "Jong Hwi Kim · Seoul National University · 2026년 10월 22일", size=12, color=MUTED)
notes(s, """[0:40] 표지
"오늘은 새로 만든 floor read-out PBISC 데이터 세트를 기존 데이터 세트들과 한 자리에서 비교합니다. 비교 대상은 세 가지입니다. 이전 PBISC K = 4,000, 경쟁 방법인 scSemiProfiler, 그리고 bulk의 질병 변화를 모든 세포형에 그대로 복사한 가장 단순한 기준선입니다."
색 규칙: 파랑 floor, 빨강 이전 K4000, 보라 scSemiProfiler, 주황 bulk broadcast, 회색·남색 실제 세포.
"모든 수치는 같은 실제 대조군, 즉 참조 환자 절반과 겹치지 않는 다른 환자 절반에 대해 같은 파이프라인으로 다시 계산한 것입니다."
""")

# ================= 2. summary
s = new(); header(s, "요약", "floor는 대부분 앞서지만, 과검출과 단핵구 결함이 남았다", 2)
three_cards(s, 1.75, 3.95, [
    ("1", "가장 크게 나아진 것", f"저발현 유전자 질병 효과 부호 일치 {m6('rank > 2,000','prev','sign_agree_DEG')} → {m6('rank > 2,000','floor','sign_agree_DEG')}. "
          f"iLISI {jj('prev4k','iLISI_post')} → {jj('floor','iLISI_post')} (실제 {jj('realU','iLISI_post')}). LR 검출 Jaccard {lr('prev4k','detected_jaccard')} → {lr('floor','detected_jaccard')}."),
    ("2", "기존 대비 확실한 우위", f"scSemiProfiler는 세포당 유전자 {gpc('scsemi')}개로 희소해 C2ST {cc2('scsemi')}, iLISI 1.00. "
          f"bulk 복사 기준선 대비 floor의 추가 정보: partial r {pr('all','floor')}."),
    ("3", "남은 약점", f"과검출(세포당 유전자 {gpc('floor')} vs 실제 {gpc('realU')}), density {dc(PCAs,'floor','density')} (실제 {dc(PCAs,'realU','density')}). "
          f"패혈증 단핵구 MHC-II 과다, T 마커 혼입 {pct(MY.loc[('sepsis','floor'),'T_marker_any'])} (실제 {pct(MY.loc[('sepsis','realT1f'),'T_marker_any'])}). 환자 식별은 개선 없음.")])
text(s, 0.89, 5.95, 11.55, 0.8, [[("한 문장  ", {"bold": True, "color": INK}), ("floor는 '실제처럼 보이고 실제처럼 반응하는' 세포에 가장 가까운 데이터 세트이며, 다음 과제는 과검출 보정과 단핵구 정체성.", {})]],
     size=14, color=BODY, line=1.15)
notes(s, f"""[1:30] 요약
"결론은 세 줄입니다. 첫째, floor는 저발현 유전자의 질병 반응과 세포 간 상호작용, 통합 혼합도에서 이전 PBISC를 크게 앞섭니다. 둘째, scSemiProfiler와 bulk 복사 기준선보다 거의 모든 축에서 낫고, 특히 bulk를 통제하고도 floor만의 세포형별 정보가 남습니다. 셋째, 약점도 분명합니다. 유전자를 실제보다 많이 검출하고, 패혈증 단핵구에서 HLA class II가 과하게 나오며 T 마커가 섞입니다. 환자 식별 능력은 이전과 같습니다."
질문이 나오면: floor와 이전 K4000은 input depth, read-out, draw 수가 동시에 다르므로 개선을 read-out 하나의 효과로 말하지 않는다.""")

# ================= 3. data sets
s = new(); header(s, "비교 대상", "네 가지 생성 데이터 세트와 실제 대조군", 3, "실제 A = 참조 절반, 실제 B = 겹치지 않는 다른 환자 절반(상한 대조) · 패혈증 성인만")
rows = [("데이터 세트", "입력 / read-out", "샘플 (패혈증 / COVID)", "분석 세포", "세포당 유전자 (중앙값)"),
        ("PBISC floor (신규)", "bulk ×256, K4000, draw0\nx<1e-4→0, <1→1, 그 외 floor", "121 / 159 (bulk-only)", "60,000 / 60,000", gpc("floor")),
        ("PBISC prev K4000 (기존)", "K4000, Poisson read-out\n5 draw 합", "121 / 164 (bulk-only)", "60,000 / 60,000", gpc("prev4k")),
        ("scSemiProfiler", "strict 5 pairs, draw0\n기댓값 → Poisson 1회", "34 / 102 (실제 A 환자 제외)", "51,000 / 60,000", gpc("scsemi")),
        ("bulk broadcast", "bulk logFC를 모든 세포형에 복사", "121 / 164 (floor와 같은 bulk)", "– (pseudobulk만)", "–"),
        ("실제 B (대조)", "held-out scRNA, 다른 환자", "20 / 57", "59,457 / 56,171", gpc("realU"))]
table(s, rows, 0.9, 1.6, [2.6, 3.3, 2.6, 1.6, 1.75], h_row=0.66, bold_cols=(0,), font=12, head_font=12, hl_rows=(1,))
text(s, 0.89, 5.85, 11.55, 0.9, [[("주의  ", {"bold": True, "color": INK}),
     ("scSemiProfiler는 정답이 측정된 바로 그 샘플을 다시 만든다(더 쉬운 조건). floor·prev는 scRNA가 없는 bulk-only 샘플에서 생성. "
      "scSemi 원 출력은 세포당 비영 유전자 약 93개(native zero fraction 0.984)로 원래부터 희소.", {})]], size=13, color=BODY, line=1.1)
notes(s, """[1:30] 비교 대상
"floor와 이전 K4000은 같은 bulk-only 샘플에서 만든 세포입니다. 이 샘플들은 단일세포 정답이 없어서, 실제 대조군은 다른 환자들의 실제 단일세포입니다."
"scSemiProfiler는 구조가 다릅니다. 실제 단일세포가 있는 샘플을 다시 만들기 때문에 오히려 쉬운 조건입니다. 그래서 참조 절반 환자와 겹치는 샘플은 뺐습니다."
"bulk broadcast는 모델이 아니라 기준선입니다. bulk에서 본 질병 변화를 모든 세포형에 똑같이 붙인 것으로, 모델이 bulk 이상의 정보를 주는지 확인하는 하한입니다. 이번에 COVID bulk 원본도 Zenodo에서 받아 두 코호트 모두 계산했습니다."
마지막 열이 중요: scSemi는 세포당 유전자가 100개 안팎으로, 실제의 4분의 1 수준.""")

# ================= 4. read-out / gene count
s = new(); header(s, "read-out", "절벽은 사라졌지만, 이제는 실제보다 많이 검출한다", 4, "왼쪽: 세포당 검출 유전자 (panel 유전자, 로그 축) · 오른쪽: 저발현 게이트 (Result 1)")
pic(s, "s3_genes_per_cell", 0.89, 1.65, w=7.4)
side_cards(s, 8.6, 1.65, 3.85, [
    ("저발현 절벽 해소", "순위 4,000+ 유전자 CP10k 비율\n이전 0.00 / 0.00 → floor 0.91 / 0.63\n생성 안 된 유전자 1,053 / 1,100 → 308 / 299"),
    ("대신 과검출", "검출률 1.2–1.7× (깊이 맞춘 실제 대비)\n세포당 유전자 floor 468 vs 깊이 맞춘 실제 343 (패혈증)"),
    ("scSemi는 반대 극단", f"세포당 유전자 {gpc('scsemi')} (실제 {gpc('realU')})\n원 출력부터 희소: 최소 양수 rate 2.04")], h=1.55)
notes(s, """[1:40] read-out
"floor read-out의 출발점입니다. 이전 K4000은 발현 순위 2,000 이후 유전자가 거의 생성되지 않는 절벽이 있었고, floor는 이걸 없앴습니다."
"대가는 과검출입니다. 깊이를 맞춘 실제 세포보다 유전자를 1.2–1.7배 많이 검출합니다. 0.0001 이상의 아주 작은 값까지 1로 올리는 규칙 때문입니다."
"scSemiProfiler는 반대 극단입니다. 원 출력에서 세포당 0이 아닌 유전자가 약 90개라, 그래프 오른쪽처럼 100개 안팎에 몰려 있습니다."
"이 한 장이 뒤에 나오는 대부분의 차이를 설명합니다." """)

# ================= 5. integration
s = new(); header(s, "통합", "floor는 네 통합 방법 모두에서 실제 수준으로 섞인다", 5, "iLISI = 소스 혼합도(최대 2) · 오른쪽 = 생성 세포의 샘플별 iLISI 중앙값, 환자 bootstrap 95% CI")
pic(s, "s3_integration", 0.89, 1.6, w=11.55)
rows = [("iLISI", "Harmony", "Scanorama", "scVI", "BBKNN", "샘플별 (생성 세포)"),
        ("PBISC floor", il("Harmony", "floor"), il("Scanorama", "floor"), il("scVI", "floor"), il("BBKNN (graph)", "floor"), ps("floor")),
        ("PBISC prev K4000", il("Harmony", "prev4k"), il("Scanorama", "prev4k"), il("scVI", "prev4k"), il("BBKNN (graph)", "prev4k"), ps("prev4k")),
        ("scSemiProfiler", il("Harmony", "scsemi"), il("Scanorama", "scsemi"), il("scVI", "scsemi"), il("BBKNN (graph)", "scsemi"), ps("scsemi")),
        ("실제 B", il("Harmony", "realU"), il("Scanorama", "realU"), il("scVI", "realU"), il("BBKNN (graph)", "realU"), ps("realU"))]
table(s, rows, 0.9, 4.95, [2.4, 1.8, 1.8, 1.8, 1.8, 1.95], h_row=0.33, bold_cols=(0,), font=11.5, head_font=11.5, hl_rows=(1,))
notes(s, f"""[1:40] 통합
"통합 결과가 Harmony 하나에 기대는지 확인하려고 Scanorama, scVI, BBKNN까지 네 가지로 돌렸습니다. floor는 네 방법 모두에서 이전 K4000보다 높고, BBKNN과 Harmony에서는 실제끼리 수준에 가깝습니다."
"scSemiProfiler는 어느 방법으로도 1.0, 즉 전혀 섞이지 않습니다. 희소한 세포가 실제 세포와 다른 공간에 따로 모이기 때문입니다."
"오른쪽 막대는 샘플마다 따로 본 혼합도입니다. floor는 {ps('floor')}로 신뢰구간이 매우 좁고, 거의 모든 샘플이 1.5를 넘습니다. 특정 샘플 몇 개가 평균을 끌어올린 게 아니라는 뜻입니다."
주의: iLISI는 '실제처럼 보이는가'이지 '그 환자를 닮았는가'가 아님(슬라이드 13).""")

# ================= 6. UMAP + cluster
s = new(); header(s, "클러스터 매칭", "floor는 실제 세포 지도 위에 겹치고, scSemi는 따로 떨어진다", 6, "pairwise Harmony 후 UMAP (실제 A 회색) · Leiden 해상도 1")
pic(s, "s3_umap", 0.89, 1.55, h=5.2)
side_cards(s, 6.75, 1.6, 5.7, [
    ("클러스터 빈도 r (실제 A 대비)", f"floor {jj('floor','cluster_freq_pearson')} · prev {jj('prev4k','cluster_freq_pearson')}\nscSemi {jj('scsemi','cluster_freq_pearson')} · 실제 B {jj('realU','cluster_freq_pearson')}"),
    ("섞인 클러스터(20–80%)에 든 생성 세포", f"floor {jj('floor','src_cells_in_mixed_clusters')} · prev {jj('prev4k','src_cells_in_mixed_clusters')}\nscSemi {jj('scsemi','src_cells_in_mixed_clusters')} (패혈증 절반 이상이 생성 세포만의 클러스터)"),
    ("라벨 전이 일치 (6개 세포형)", f"floor {jj('floor','six_agree_with_own_label')} · prev {jj('prev4k','six_agree_with_own_label')} · scSemi {jj('scsemi','six_agree_with_own_label')}\nscSemi 패혈증 DC 일치 {J[('sepsis','scsemi')]['six_agree_by_type']['DC']:.2f}")], h=1.5)
notes(s, """[1:30] 클러스터 매칭
"UMAP에서 회색이 실제 참조 세포입니다. 오른쪽 열의 floor는 회색 위에 거의 그대로 겹치고, 왼쪽 열의 scSemi는 따로 떨어진 섬을 만듭니다."
"숫자로 보면 클러스터별 세포 비율의 상관이 floor는 실제 B와 비슷한 수준이고, scSemi는 패혈증에서 0에 가깝습니다. 패혈증 scSemi 세포의 절반 이상이 생성 세포만으로 된 클러스터에 들어갑니다."
"라벨 전이도 floor가 실제 B 수준입니다. scSemi는 패혈증 DC가 거의 전부 다른 세포형으로 옮겨집니다." """)

# ================= 7. realism
s = new(); header(s, "realism", "floor는 실제와 가장 덜 구분되지만, 아직 분류기는 구분해 낸다", 7, "C2ST = 실제 A vs 소스 분류기 AUC (0.5 = 구분 불가, 샘플 단위 CV) · density / coverage k = 5, 보정 전 PCA")
pic(s, "s3_realism", 0.89, 1.6, w=11.55)
rows = [("지표", "scSemi", "prev K4000", "floor", "실제 B"),
        ("C2ST AUC (중앙값)", cc2("scsemi"), cc2("prev4k"), cc2("floor"), cc2("realU")),
        ("coverage (높을수록 분포를 덮음)", dc(PCAs, "scsemi", "coverage"), dc(PCAs, "prev4k", "coverage"), dc(PCAs, "floor", "coverage"), dc(PCAs, "realU", "coverage")),
        ("density (실제와 가까울수록 좋음)", dc(PCAs, "scsemi", "density"), dc(PCAs, "prev4k", "density"), dc(PCAs, "floor", "density"), dc(PCAs, "realU", "density"))]
table(s, rows, 0.9, 5.0, [3.6, 1.95, 1.95, 1.95, 1.95], h_row=0.36, bold_cols=(0,), font=12, head_font=12)
notes(s, f"""[1:40] realism
"분류기로 실제와 생성 세포를 가려보면, 이전 K4000과 scSemi는 거의 완벽하게 구분됩니다. floor는 {cc2('floor')}로 내려왔지만 실제 B의 {cc2('realU')}와는 아직 차이가 있습니다."
"coverage는 실제 분포를 얼마나 덮는지입니다. floor가 {dc(PCAs,'floor','coverage')}로 가장 높고, scSemi는 거의 덮지 못합니다."
"density는 반대 방향의 경고입니다. floor는 실제보다 낮습니다. 실제 세포가 몰려 있는 곳보다 바깥쪽에 넓게 퍼져 있다는 뜻이고, 과검출과 같은 방향입니다."
참고: 실제 B를 각 모델 깊이로 낮춘 대조도 0.5–0.7이라, 깊이 차이만으로는 0.85를 설명하지 못한다.""")

# ================= 8. co-expression & LR
s = new(); header(s, "유전자 관계와 세포 간 상호작용", "floor는 저발현 유전자의 상관 구조와 LR 상호작용을 실제 수준으로 회복", 8, "상관행렬 유사도 = 실제 A와 상관계수 행렬의 r · LR = LIANA consensus, 두 세포형 모두 10% 이상 검출")
rows = [("지표 (패혈증 / COVID)", "scSemi", "prev K4000", "floor", "실제 B"),
        ("유전자 상관 유사도, 상위 300 유전자", co("top 300 genes", "scsemi"), co("top 300 genes", "prev4k"), co("top 300 genes", "floor"), co("top 300 genes", "realU")),
        ("유전자 상관 유사도, 순위 1,001–1,300", co("ranks 1,001-1,300", "scsemi"), co("ranks 1,001-1,300", "prev4k"), co("ranks 1,001-1,300", "floor"), co("ranks 1,001-1,300", "realU")),
        ("LR 검출 Jaccard", lr("scsemi", "detected_jaccard"), lr("prev4k", "detected_jaccard"), lr("floor", "detected_jaccard"), lr("realU", "detected_jaccard")),
        ("실제 상위 200 LR 회수율", lr("scsemi", "recall_top200_realA"), lr("prev4k", "recall_top200_realA"), lr("floor", "recall_top200_realA"), lr("realU", "recall_top200_realA")),
        ("LR 점수 Spearman", lr("scsemi", "spearman_score_on_realA_detected"), lr("prev4k", "spearman_score_on_realA_detected"), lr("floor", "spearman_score_on_realA_detected"), lr("realU", "spearman_score_on_realA_detected")),
        ("질병에 따른 LR 변화 Spearman", ld("scsemi"), ld("prev4k"), ld("floor"), "1")]
table(s, rows, 0.9, 1.65, [4.15, 1.85, 1.85, 1.85, 1.85], h_row=0.5, bold_cols=(0,), font=13, head_font=12.5, hl_rows=(2, 3))
text(s, 0.89, 5.75, 11.55, 0.9, [[("읽는 법  ", {"bold": True, "color": INK}),
     ("저발현 유전자 상관과 LR 검출은 floor가 실제 B에 근접. 다만 질병에 따른 LR 변화는 이전 K4000과 같은 수준으로, 개선되지 않았다.", {})]], size=14, color=BODY, line=1.1)
notes(s, """[1:30] co-expression, LR
"유전자끼리 함께 움직이는 구조를 보면, 많이 발현되는 유전자는 세 모델이 비슷하지만 순위 1,000 이후 유전자에서 floor가 실제 B 수준으로 올라옵니다. 절벽이 없어진 결과입니다."
"ligand–receptor 상호작용도 같은 이야기입니다. 검출된 상호작용 집합의 겹침이 0.5 안팎에서 0.75–0.8로 올라왔습니다. scSemi는 희소성 때문에 상호작용을 거의 검출하지 못합니다."
"정직하게 말씀드릴 부분은 마지막 줄입니다. 질병에서 상호작용이 어떻게 변하는지는 floor가 이전보다 낫지 않습니다." """)

# ================= 9. M6
s = new(); header(s, "질병 효과 (M6)", "저발현 유전자의 질병 반응은 floor만 실제 방향을 따른다", 9, "세포형별 샘플 pseudobulk DE (case vs control) · 막대 = 5개 세포형 중앙값, 수염 = 환자 bootstrap 95% CI")
pic(s, "s3_m6_bootstrap", 0.89, 1.6, w=11.55)
rows = [("순위 2,000 이후 유전자", "scSemi", "prev K4000", "floor", "bulk broadcast"),
        ("logFC r", m6("rank > 2,000", "scsemi", "logFC_r"), m6("rank > 2,000", "prev", "logFC_r"), m6("rank > 2,000", "floor", "logFC_r"), m6("rank > 2,000", "broadcast", "logFC_r")),
        ("방향 일치 (실제 DEG)", m6("rank > 2,000", "scsemi", "sign_agree_DEG"), m6("rank > 2,000", "prev", "sign_agree_DEG"), m6("rank > 2,000", "floor", "sign_agree_DEG"), m6("rank > 2,000", "broadcast", "sign_agree_DEG"))]
table(s, rows, 0.9, 5.0, [3.0, 2.0, 2.0, 2.0, 2.15], h_row=0.36, bold_cols=(0,), font=12, head_font=12)
text(s, 0.89, 6.25, 11.55, 0.5, f"floor − prev logFC r: 패혈증 {ci('sepsis','rank > 2,000','logFC_r','floor - prev')}, COVID {ci('covid','rank > 2,000','logFC_r','floor - prev')} · "
     f"floor − scSemi: {ci('sepsis','rank > 2,000','logFC_r','floor - scsemi')}, {ci('covid','rank > 2,000','logFC_r','floor - scsemi')}", size=12, color=BODY)
notes(s, f"""[2:00] 질병 효과
"논문에서 가장 중요한 질문입니다. 생성 세포가 질병에 따른 세포형별 유전자 변화를 재현하는가."
"많이 발현되는 유전자에서는 floor와 이전 모델이 비슷합니다. 차이는 순위 2,000 이후 저발현 유전자에서 납니다. 실제 DEG와 방향이 같은 비율이 floor는 0.90으로, 이전 K4000보다 높습니다."
"환자를 다시 뽑는 bootstrap으로 차이를 검정했습니다. 패혈증은 floor − prev가 {ci('sepsis','rank > 2,000','logFC_r','floor - prev')}로 확실하고, COVID는 {ci('covid','rank > 2,000','logFC_r','floor - prev')}로 경계선입니다."
"scSemi는 저발현 유전자의 logFC가 정확히 0인 경우가 많아 방향 일치가 매우 낮습니다. 같은 샘플을 다시 만드는 유리한 조건인데도 그렇습니다." """)

# ================= 10. broadcast
s = new(); header(s, "bulk 이상의 정보", "bulk를 통제해도 floor의 세포형별 신호는 남는다", 10,
                  "bulk broadcast = 같은 bulk-only 샘플의 bulk logFC를 모든 세포형에 동일하게 사용 · 패혈증 GSE279448, COVID Gisby 2022 (Zenodo 6497251)")
pic(s, "s3_broadcast_partial", 0.89, 1.55, h=5.15)
side_cards(s, 9.95, 1.6, 2.5, [
    ("bulk만으로", f"logFC r {rb('all')}\n방향 일치 {m6('all','broadcast','sign_agree_DEG')}"),
    ("floor − bulk", f"r: 패혈증 {ci('sepsis','all','logFC_r','floor - broadcast')}\nCOVID {ci('covid','all','logFC_r','floor - broadcast')}\n방향: {ci('sepsis','all','sign_agree_DEG','floor - broadcast')}\n{ci('covid','all','sign_agree_DEG','floor - broadcast')}"),
    ("bulk 통제 partial r", f"floor {pr('all','floor')}\nprev {pr('all','prev')}\nscSemi {pr('all','scsemi')}")], h=1.62, body=11, head=13)
notes(s, f"""[2:00] bulk broadcast
"리뷰어가 가장 먼저 물을 질문입니다. 생성 세포가 그냥 bulk의 변화를 모든 세포형에 퍼뜨린 것 아니냐."
"그래서 bulk logFC를 모든 세포형에 똑같이 붙인 기준선을 만들고, 그걸 통제한 뒤에도 모델이 실제 세포형별 logFC를 설명하는지 partial r로 봤습니다. 이번에 COVID bulk 원본(Gisby 2022, Zenodo)을 받아 COVID까지 계산했습니다."
"결과는 두 코호트에서 다르게 나타납니다. 패혈증에서는 bulk만으로도 방향은 0.83까지 맞지만 크기가 맞지 않고, floor가 r을 {ci('sepsis','all','logFC_r','floor - broadcast')}만큼 올립니다. COVID에서는 반대로 bulk의 방향이 0.62에 그치고, floor가 방향 일치를 {ci('covid','all','sign_agree_DEG','floor - broadcast')} 올립니다."
"bulk를 통제한 partial r은 floor가 {pr('all','floor')}로 세 모델 중 가장 높습니다. 즉 floor는 bulk에 없는 세포형별 정보를 만든다고 말할 수 있습니다."
약점: 패혈증 DC의 partial r은 0.16으로 낮다(그림 왼쪽 위).""")

# ================= 11. composition
s = new(); header(s, "세포 구성", "6개 세포형 구성과 질병에 따른 구성 변화는 floor가 가장 정확", 11, "TV = 질병군 평균 구성의 total variation 오차 · 산점도 = 세부 세포 상태(전이 라벨)의 case − control 변화, 회색 = 실제 B")
pic(s, "s3_composition_shift", 0.89, 1.55, h=5.2)
rows = [("(패혈증 / COVID)", "scSemi", "prev", "floor"),
        ("6형 구성 오차 TV", cp("six", "scsemi", "TV_case"), cp("six", "prev", "TV_case"), cp("six", "floor", "TV_case")),
        ("6형 질병 변화 r", cp("six", "scsemi", "shift_r"), cp("six", "prev", "shift_r"), cp("six", "floor", "shift_r")),
        ("세부 상태 변화 r", cp("fine", "scsemi", "shift_r"), cp("fine", "prev", "shift_r"), cp("fine", "floor", "shift_r"))]
table(s, rows, 9.45, 1.65, [1.25, 0.92, 0.92, 0.92], h_row=0.52, bold_cols=(0,), font=10.5, head_font=10.5)
text(s, 9.45, 4.3, 3.0, 2.2, [f"세부 상태 실제 B (같은 전이 경로) {cp('fine','realB_tr','shift_r')}",
                              f"COVID 세부 floor − prev {CC.loc[('covid','fine','shift_r','floor - prev'),'med']:+.2f} [{CC.loc[('covid','fine','shift_r','floor - prev'),'lo']:+.2f}, {CC.loc[('covid','fine','shift_r','floor - prev'),'hi']:+.2f}]: 유의차 없음",
                              "패혈증 vs ICU 대조 구성 변화: 세 모델 모두 재현 못 함"], size=11.5, color=BODY, line=1.1, after=6)
notes(s, """[1:40] 세포 구성
"bulk에서 세포 구성을 맞추는 것은 deconvolution도 하는 일이라, 여기서는 질병군 평균 구성과 질병에 따른 변화가 맞는지를 봤습니다."
"6개 세포형 수준에서는 floor가 구성 오차와 질병 변화 모두 가장 정확합니다. scSemi는 패혈증에서 T 세포가 70%를 차지하는 등 구성이 크게 틀어집니다."
"세부 상태로 내려가면 패혈증은 floor가 좋고, COVID는 floor가 이전 모델보다 약간 낮지만 신뢰구간이 겹쳐 차이가 없다고 봐야 합니다."
"패혈증과 ICU 대조군의 구성 차이는 실제 대조 샘플이 8개뿐이라 어떤 모델도 재현하지 못합니다." """)

# ================= 12. claims & rare
s = new(); header(s, "논문 결과와 희귀 세포", "floor는 논문 결과의 크기를 가장 잘 재현하고, 희귀 상태는 혼재", 12, "같은 분류기로 case − control 비교 (실제 B 대비) · 희귀 = 부모 계통 안 분류기, 라벨 전이 없음, AUC ≥ 0.8")
rows = [("", "scSemi", "prev K4000", "floor"),
        ("패혈증 SEP claim (방향, 크기)", cl("sepsis", "scsemi"), cl("sepsis", "prev4k"), cl("sepsis", "floor")),
        ("COVID IFN claim (방향, 크기)", cl("covid", "scsemi"), cl("covid", "prev4k"), cl("covid", "floor")),
        ("희귀 상태 0.5–2× (패혈증 / COVID)", rq("scsemi"), rq("prev4k"), rq("floor"))]
table(s, rows, 0.9, 1.6, [3.6, 2.6, 2.6, 2.6], h_row=0.46, bold_cols=(0,), font=13, head_font=12.5)
rs = lambda c, sub: " · ".join(f"{k} {rr.loc[(c, sub), k]:.2f}×" for k, lab in (("prev4k", ""), ("floor", ""), ("scsemi", "")) if not np.isnan(rr.loc[(c, sub), k]))
side_cards(s, 0.9, 3.85, 5.6, [
    ("floor가 나아진 희귀 상태", f"NK_c04-MKI67 {rr.loc[('sepsis','NK_c04-MKI67'),'prev4k']:.2f}× → {rr.loc[('sepsis','NK_c04-MKI67'),'floor']:.2f}× · 대식세포 {rr.loc[('sepsis','Macrophage'),'prev4k']:.2f}× → {rr.loc[('sepsis','Macrophage'),'floor']:.2f}×\n"
                              f"DC1 {rr.loc[('covid','DC1'),'prev4k']:.2f}× → {rr.loc[('covid','DC1'),'floor']:.2f}× · MAIT {rr.loc[('covid','MAIT'),'prev4k']:.2f}× → {rr.loc[('covid','MAIT'),'floor']:.2f}×")], h=1.3, accent=BLUE, body=12)
side_cards(s, 6.75, 3.85, 5.7, [
    ("floor가 나빠지거나 여전히 없는 상태", f"γδT {rr.loc[('sepsis','gdT'),'prev4k']:.2f}× → {rr.loc[('sepsis','gdT'),'floor']:.2f}× · Mono_c06-MKI67 {rr.loc[('sepsis','Mono_c06-CD14-MKI67'),'prev4k']:.2f}× → {rr.loc[('sepsis','Mono_c06-CD14-MKI67'),'floor']:.2f}×\n"
                                           f"pDC (패혈증) {rr.loc[('sepsis','pDC'),'floor']:.2f}× · CD14mono_dex {rr.loc[('covid','CD14mono_dex'),'floor']:.2f}× (모든 모델 부재)")], h=1.3, accent=RED, body=12)
text(s, 0.89, 5.45, 11.55, 1.2, [[("scSemi  ", {"bold": True, "color": INK}), ("패혈증 marker 양성 비율이 두 군 모두 ~0 (희소성) → claim 1/5. 반면 COVID IFN claim은 크기까지 잘 재현 (ISG는 고발현 유전자).", {})],
                                 [("해석  ", {"bold": True, "color": INK}), ("floor는 크기를 0.30× → 0.80× (패혈증)로 끌어올렸지만 SEP-2 방향 하나를 놓침. 희귀 상태는 '증식형 NK·DC1·MAIT는 회복, γδT·증식형 단핵구는 후퇴'.", {})]],
     size=13, color=BODY, line=1.12, after=6)
notes(s, """[1:50] claims, 희귀 세포
"논문에 나온 질병 결과를 같은 분류기로 다시 재현해 봤습니다. floor는 크기를 실제의 0.8배까지 재현해 이전의 0.3배보다 훨씬 가깝지만, 다섯 개 중 SEP-2 하나의 방향을 놓칩니다."
"scSemi는 COVID 인터페론 반응은 잘 재현합니다. 인터페론 유전자가 많이 발현되는 유전자라 희소성의 영향을 덜 받기 때문입니다. 반면 패혈증 claim은 marker가 거의 검출되지 않아 재현이 안 됩니다."
"희귀 상태는 혼재입니다. 증식형 NK, DC1, MAIT, 대식세포는 floor에서 나아졌고, γδT와 증식형 단핵구는 오히려 이전보다 나빠졌습니다. dex 단핵구는 어떤 모델에도 없습니다." """)

# ================= 13. patient identity
s = new(); header(s, "환자 특이성", "환자 고유 신호는 우연의 9배지만, floor에서 더 좋아지지는 않았다", 13, "COVID 32 donor: 생성(bulk-only 시점) vs 실제 scRNA(대부분 다른 시점), donor 안정 유전자 100개 · donor bootstrap 2,000회")
three_cards(s, 1.75, 3.3, [
    (f"{FPb.loc['hit_floor','estimate']:.2f}", "floor donor 식별 top-1", f"95% CI [{FPb.loc['hit_floor','lo']:.2f}, {FPb.loc['hit_floor','hi']:.2f}]\n우연 0.03"),
    (f"{FPb.loc['hit_prev','estimate']:.2f}", "prev K4000 donor 식별", f"95% CI [{FPb.loc['hit_prev','lo']:.2f}, {FPb.loc['hit_prev','hi']:.2f}]\nfloor − prev {FPb.loc['dhit','estimate']:+.2f} [{FPb.loc['dhit','lo']:+.2f}, {FPb.loc['dhit','hi']:+.2f}]"),
    (" / ".join(f"{r.estimate:.2f}" for _, r in SI.iterrows()), "scSemi 같은 샘플 재식별", "COVID / 패혈증, 우연 ≤ 0.013\n같은 샘플을 다시 만든 것이라 더 쉬운 과제: 직접 비교 불가")], numsize=30)
text(s, 0.89, 5.35, 11.55, 1.3, [[("의미  ", {"bold": True, "color": INK}), ("iLISI·C2ST가 좋아진 것은 '실제처럼 보이는가'의 개선이지 '그 환자를 닮았는가'의 개선이 아니다. "
     "환자 특이성은 floor read-out이 아니라 condition 경로(bulk 인코딩, LoRA)에서 다뤄야 할 과제.", {})]], size=14, color=BODY, line=1.15)
notes(s, f"""[1:30] 환자 특이성
"COVID에는 생성 세포와 실제 단일세포가 모두 있는 donor가 32명 있습니다. 대부분 다른 시점이지만, donor마다 안정적인 유전자로 '이 생성 세포가 누구 것인지'를 맞혀 봤습니다."
"floor는 32명 중 {FPb.loc['hit_floor','estimate']:.2f}, 이전은 {FPb.loc['hit_prev','estimate']:.2f}로 우연(0.03)보다 아홉 배 높습니다. 하지만 floor − prev 차이의 신뢰구간이 0을 크게 포함합니다. 환자 신호는 유지되지만 좋아지지는 않았습니다."
"scSemi의 0.13–0.15는 같은 샘플을 다시 찾는 과제라 더 쉬운 문제이고, 그래서 직접 비교하지 않았습니다."
메시지: read-out 개선은 환자 특이성을 바꾸지 않는다. 이건 모델 condition 쪽 과제.""")

# ================= 14. myeloid
s = new(); header(s, "단핵구 진단", "패혈증 단핵구: MHC-II 과다와 T 마커 혼입이 핵심 결함", 14, "깊이 맞춘 실제 B 대비 유전자별 잔차 · 회색 막대 = 실제 A − 실제 B (모델 없는 잡음)")
pic(s, "s3_myeloid", 0.89, 1.6, w=11.55)
side_cards(s, 0.9, 5.0, 3.7, [("MHC-II 과다", f"HLA-DRA +{mg.loc['HLA-DRA','residual']:.2f}, CD74 +{mg.loc['CD74','residual']:.2f} (log1p CP10k)\n질병 감소 {MY.loc[('sepsis','floor'),'mhc_shift_case_minus_ctrl_sdA']:.2f} SD vs 실제 {MY.loc[('sepsis','realU'),'mhc_shift_case_minus_ctrl_sdA']:.2f} SD")], h=1.55, accent=RED, body=11.5, head=13)
side_cards(s, 4.82, 5.0, 3.7, [("T 마커 혼입", f"T 마커 UMI 보유 단핵구\n패혈증 floor {pct(MY.loc[('sepsis','floor'),'T_marker_any'])} vs 실제 {pct(MY.loc[('sepsis','realT1f'),'T_marker_any'])}\nCOVID {pct(MY.loc[('covid','floor'),'T_marker_any'])} vs {pct(MY.loc[('covid','realT1f'),'T_marker_any'])} (차이 작음)")], h=1.55, accent=RED, body=11.5, head=13)
side_cards(s, 8.74, 5.0, 3.7, [("그 외", f"ISG(IFI6, ISG15), AP-1(FOS, JUN) 상승\nscSemi 단핵구: 유전자 {MY.loc[('sepsis','scsemi'),'median_genes']:.0f}개, 정체성 점수 하위 5% 미만 {pct(MY.loc[('sepsis','scsemi'),'mono_below_realB_p5'])}")], h=1.55, accent=GREY, body=11.5, head=13)
notes(s, f"""[2:00] 단핵구 진단
"이전 보고서에서 패혈증 골수계가 floor의 약한 곳이었습니다. 깊이를 맞춘 실제 세포와 유전자 하나씩 비교해 원인을 찾았습니다."
"첫째, HLA class II 유전자가 실제보다 높게 나옵니다. HLA-DRA가 약 0.9 log 단위 높습니다. 패혈증의 대표적 특징이 단핵구 HLA-DR 감소인데, floor에서는 이 감소가 실제의 70% 정도로 약하게 나옵니다."
"둘째, 단핵구의 {pct(MY.loc[('sepsis','floor'),'T_marker_any'])}가 T 세포 마커를 갖고 있습니다. 실제는 2% 미만입니다. 이전 모델도 같은 문제가 있어 floor가 새로 만든 문제는 아니고, COVID에서는 나타나지 않습니다."
"셋째, 인터페론 유전자와 FOS·JUN 같은 즉시 반응 유전자도 높습니다."
"이 세 가지가 다음 threshold sweep에서 같이 봐야 할 구체적인 판독 지표입니다." """)

# ================= 15. scorecard
s = new(); header(s, "한 장 요약", "3개 데이터 세트 스코어카드 (패혈증 / COVID)", 15, "강조 행 = floor가 세 데이터 세트 중 최고 · results/scorecard_three_models.csv")
rows = [("지표", "scSemi", "prev K4000", "floor", "실제 / 기준"),
        ("iLISI, Harmony", jj("scsemi", "iLISI_post"), jj("prev4k", "iLISI_post"), jj("floor", "iLISI_post"), jj("realU", "iLISI_post")),
        ("C2ST AUC (낮을수록 좋음)", cc2("scsemi"), cc2("prev4k"), cc2("floor"), cc2("realU")),
        ("coverage, PCA", dc(PCAs, "scsemi", "coverage"), dc(PCAs, "prev4k", "coverage"), dc(PCAs, "floor", "coverage"), dc(PCAs, "realU", "coverage")),
        ("유전자 상관, 순위 1,001–1,300", co("ranks 1,001-1,300", "scsemi"), co("ranks 1,001-1,300", "prev4k"), co("ranks 1,001-1,300", "floor"), co("ranks 1,001-1,300", "realU")),
        ("M6 logFC r, 전체", m6("all", "scsemi", "logFC_r"), m6("all", "prev", "logFC_r"), m6("all", "floor", "logFC_r"), "bulk " + m6("all", "broadcast", "logFC_r")),
        ("M6 방향 일치, 순위 > 2,000", m6("rank > 2,000", "scsemi", "sign_agree_DEG"), m6("rank > 2,000", "prev", "sign_agree_DEG"), m6("rank > 2,000", "floor", "sign_agree_DEG"), "bulk " + m6("rank > 2,000", "broadcast", "sign_agree_DEG")),
        ("bulk 통제 partial r", pr("all", "scsemi"), pr("all", "prev"), pr("all", "floor"), "–"),
        ("6형 질병 구성 변화 r", cp("six", "scsemi", "shift_r"), cp("six", "prev", "shift_r"), cp("six", "floor", "shift_r"), "–"),
        ("LR 검출 Jaccard", lr("scsemi", "detected_jaccard"), lr("prev4k", "detected_jaccard"), lr("floor", "detected_jaccard"), lr("realU", "detected_jaccard")),
        ("질병 LR 변화 Spearman", ld("scsemi"), ld("prev4k"), ld("floor"), "1"),
        ("희귀 상태 0.5–2×", rq("scsemi"), rq("prev4k"), rq("floor"), "–"),
        ("donor 식별 (COVID)", "–", f2(FPb.loc['hit_prev', 'estimate']), f2(FPb.loc['hit_floor', 'estimate']), "우연 0.03"),
        ("단핵구 T 마커 (패혈증 / COVID)", my("scsemi", "T_marker_any"), my("prev4k", "T_marker_any"), my("floor", "T_marker_any"), my("realT1f", "T_marker_any"))]
table(s, rows, 0.9, 1.5, [3.45, 2.0, 2.0, 2.0, 2.1], h_row=0.355, bold_cols=(0,), font=11.5, head_font=11.5, hl_rows=(1, 2, 3, 4, 5, 6, 7, 8, 9))
notes(s, """[1:30] 스코어카드
"지금까지를 한 표로 모았습니다. 위의 아홉 줄은 floor가 세 데이터 세트 중 최고인 지표입니다."
"floor가 최고가 아닌 곳은 아래 네 줄입니다. 질병에 따른 LR 변화는 이전 모델이 약간 높고, 희귀 상태는 패혈증에서 이전 모델이 더 많이 맞힙니다. donor 식별은 비슷하고, 단핵구 T 마커 혼입은 패혈증에서 이전 모델과 floor 모두 실제보다 높습니다."
"scSemiProfiler는 COVID 인터페론 claim을 빼면 모든 줄에서 가장 낮습니다." """)

# ================= 16. strengths
s = new(); header(s, "모델의 강점", "floor PBISC가 기존 데이터 세트보다 확실히 나은 다섯 가지", 16)
items = [("1  실제 세포와 섞인다", f"iLISI {jj('floor','iLISI_post')}: 실제끼리 {jj('realU','iLISI_post')}에 근접, 4개 통합 방법 모두 일관. 샘플별로도 {ps('floor')}로 고르게 높음."),
         ("2  저발현 유전자가 살아 있다", f"절벽 해소 → 순위 1,000 이후 유전자 상관 {co('ranks 1,001-1,300','floor')} (실제 {co('ranks 1,001-1,300','realU')}), LR 검출 Jaccard {lr('floor','detected_jaccard')}."),
         ("3  저발현 유전자의 질병 반응이 맞는다", f"방향 일치 {m6('rank > 2,000','floor','sign_agree_DEG')} (prev {m6('rank > 2,000','prev','sign_agree_DEG')}, scSemi {m6('rank > 2,000','scsemi','sign_agree_DEG')}). 패혈증 차이 bootstrap 유의."),
         ("4  bulk 이상의 세포형별 정보", f"bulk 통제 partial r {pr('all','floor')} (prev {pr('all','prev')}, scSemi {pr('all','scsemi')}). COVID에선 bulk 방향 {m6('all','broadcast','sign_agree_DEG').split(' / ')[1]} → floor {m6('all','floor','sign_agree_DEG').split(' / ')[1]}."),
         ("5  세포 구성과 논문 결과의 크기", f"6형 구성 오차 {cp('six','floor','TV_case')} (prev {cp('six','prev','TV_case')}), 질병 변화 r {cp('six','floor','shift_r')}. 패혈증 claim 크기 0.30× → 0.80×.")]
for i, (h, b) in enumerate(items):
    y = 1.6 + i * 1.0; card(s, 0.89, y, 11.55, 0.88)
    text(s, 1.1, y + 0.14, 4.1, 0.6, h, size=15, color=INK, bold=True); text(s, 5.25, y + 0.14, 7.0, 0.7, b, size=12.5, color=BODY, line=1.08)
notes(s, """[1:40] 강점
"강점 다섯 가지를 정리합니다. 첫째, 실제 세포와 섞입니다. 어떤 통합 방법을 써도, 어떤 샘플을 봐도 그렇습니다. 둘째, 이전 모델에서 사라졌던 저발현 유전자가 살아나서 유전자 간 관계와 세포 간 상호작용이 회복됐습니다. 셋째, 그 저발현 유전자들이 질병에 실제와 같은 방향으로 반응합니다. 넷째, 이건 bulk를 복사한 결과가 아닙니다. bulk를 통제해도 세포형별 신호가 남습니다. 다섯째, 세포 구성과 논문 결과의 크기도 가장 정확합니다."
"논문에서는 3번과 4번이 핵심 주장, 1번과 2번이 그 전제 조건이 됩니다." """)

# ================= 17. weaknesses
s = new(); header(s, "약점과 한계", "아직 남은 다섯 가지, 그리고 이번 분석의 한계", 17)
items = [("1  과검출", f"세포당 유전자 {gpc('floor')} vs 실제 {gpc('realU')}; density {dc(PCAs,'floor','density')} < 실제 {dc(PCAs,'realU','density')}; C2ST {cc2('floor')}. 1e-4 → 1 규칙이 원인."),
         ("2  패혈증 단핵구", f"MHC-II 과다(HLA-DRA +{mg.loc['HLA-DRA','residual']:.2f}), HLA-DR 감소 약화, T 마커 혼입 {pct(MY.loc[('sepsis','floor'),'T_marker_any'])}. 패혈증 DC partial r도 낮음."),
         ("3  환자 특이성", f"donor 식별 {f2(FPb.loc['hit_floor','estimate'])} vs prev {f2(FPb.loc['hit_prev','estimate'])}, 차이 {FPb.loc['dhit','estimate']:+.2f} [{FPb.loc['dhit','lo']:+.2f}, {FPb.loc['dhit','hi']:+.2f}]: read-out으로는 개선 불가."),
         ("4  일부 희귀 상태와 질병 상호작용", f"γδT {rr.loc[('sepsis','gdT'),'floor']:.2f}×, 증식형 단핵구 {rr.loc[('sepsis','Mono_c06-CD14-MKI67'),'floor']:.2f}×, dex 단핵구 부재; 질병 LR 변화 {ld('floor')} (prev {ld('prev4k')})."),
         ("5  분석 한계", "floor와 prev는 depth·read-out·draw 수가 동시에 다름(read-out 단독 효과 미분리); draw0 1개; 패혈증 성인 HC 7개; 패혈증 vs ICU 비교 불가.")]
for i, (h, b) in enumerate(items):
    y = 1.6 + i * 1.0; card(s, 0.89, y, 11.55, 0.88, accent=RED if i < 4 else GREY)
    text(s, 1.1, y + 0.14, 4.1, 0.6, h, size=15, color=INK, bold=True); text(s, 5.25, y + 0.14, 7.0, 0.7, b, size=12.5, color=BODY, line=1.08)
notes(s, """[1:40] 약점과 한계
"약점도 다섯 가지로 정리했습니다. 첫째, 과검출입니다. 아주 작은 값까지 1로 올리는 규칙 때문에 유전자를 실제보다 많이 검출하고, 분류기가 아직 실제와 구분합니다. 둘째, 패혈증 단핵구에서 HLA class II가 과다하고 T 마커가 섞입니다. 셋째, 환자 특이성은 read-out으로는 바뀌지 않았습니다. 넷째, γδT와 증식형 단핵구 같은 일부 희귀 상태와 질병에 따른 상호작용 변화는 아직 약합니다."
"그리고 이번 분석 자체의 한계가 있습니다. floor와 이전 모델은 여러 조건이 동시에 달라서 read-out 하나의 효과로 말할 수 없고, draw도 하나만 썼습니다." """)

# ================= 18. next
s = new(); header(s, "다음 단계", "약점을 판독 지표로 삼아 threshold를 고르고, 논문 주장으로 묶는다", 18)
three_cards(s, 1.75, 3.6, [
    ("1", "threshold sweep (서버)", "1e-4 … 0.5 연속값에서 재판독\n판독 지표: 세포당 유전자, density, C2ST, 단핵구 MHC-II·T 마커, M6 방향 일치"),
    ("2", "read-out 단독 효과", "depth256 Poisson control을 같은 파이프라인으로\n→ floor 개선 중 read-out 몫 분리\n+ draw 1–4 재현성"),
    ("3", "논문 구성", "주장: 저발현 질병 반응 + bulk 이상의 정보\n전제: 혼합·realism\n한계 명시: 환자 특이성, 단핵구")])
text(s, 0.89, 5.6, 11.55, 1.0, [[("산출물  ", {"bold": True, "color": INK}), ("F:\\PBISC\\bulkonly_lora5_20261006\\analysis\\integration_report_20261012\\three_model_20261022\\ (three_model_comparison.ipynb, results/, figs/) · "
     "COVID bulk: bulk_public\\zenodo6497251_imperial_eskd\\", {})]], size=12, color=BODY, line=1.1)
notes(s, """[1:20] 다음 단계
"다음 단계는 세 가지입니다. 첫째, 서버에 저장된 연속값으로 threshold를 바꿔가며 다시 판독합니다. 오늘 찾은 약점들, 즉 세포당 유전자 수, density, 단핵구 MHC-II와 T 마커 혼입을 판독 지표로 씁니다. 둘째, depth256 Poisson control을 같은 파이프라인에 넣어 read-out만의 효과를 분리합니다. 셋째, 논문은 저발현 유전자의 질병 반응과 bulk 이상의 정보를 핵심 주장으로, 혼합과 realism을 전제로, 환자 특이성과 단핵구는 한계로 명시하는 구조를 제안합니다."
"감사합니다." """)
os.makedirs(os.path.dirname(DST), exist_ok=True); prs.save(DST); print("saved", DST, len(prs.slides))
