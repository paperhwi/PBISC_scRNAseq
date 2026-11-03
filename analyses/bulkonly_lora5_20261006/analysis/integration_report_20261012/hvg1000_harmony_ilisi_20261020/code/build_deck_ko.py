"""Korean deck for the HVG 1,000 Harmony / iLISI analysis, in the style of 'PBISC K1000 vs K4000 Integration Report_latest.pptx'
(its master and 'Title and body' layout are reused; all of its slides are dropped). usage: py -3.12 build_deck_ko.py"""
import os, copy, pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG = OUT + "/figs/deck"
TEMPLATE = "F:/PBISC/PBISC K1000 vs K4000 Integration Report_latest.pptx"
DST = "F:/PBISC/PBISC HVG1000 Harmony iLISI 보고서.pptx"
BLUE, RED, GREY, INK, BODY, MUTED, CARD = "2B83BA", "D7191C", "7F7F7F", "1F2933", "3E4C59", "7B8794", "F5F7FA"
KO = "Malgun Gothic"
R = pd.read_csv(OUT + "/results/ilisi_hvg1000_vs_hvg2000.csv")
v = lambda n, st, c: R[(R.n_hvg == n) & (R.stage == st) & (R.comparison == c)].iLISI.item()

prs = Presentation(TEMPLATE)
sld = prs.slides._sldIdLst
for s in list(sld):   # drop every template slide; unreachable parts are not saved
    prs.part.drop_rel(s.get(qn("r:id"))); sld.remove(s)
LAYOUT = [l for l in prs.slide_layouts if l.name == "Title and body"][0]
FOOT = "PBISC HVG 1,000 통합 보고서"


def style_run(r, size, color, bold=False, spc=None):
    r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = RGBColor.from_string(color); r.font.name = "Arial"
    rPr = r._r.get_or_add_rPr()
    for tag in ("a:ea", "a:cs"):
        e = rPr.find(qn(tag))
        if e is None: e = etree.SubElement(rPr, qn(tag))
        e.set("typeface", KO if tag == "a:ea" else "Arial")
    if spc: rPr.set("spc", str(spc))


def text(slide, x, y, w, h, paras, size=14, color=BODY, bold=False, align=PP_ALIGN.LEFT, spc=None, anchor=MSO_ANCHOR.TOP, line=None, after=0):
    """paras: list of str or list of [(text, {overrides})] runs."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0; tf.vertical_anchor = anchor
    for i, p in enumerate(paras if isinstance(paras, list) else [paras]):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); para.alignment = align
        if line: para.line_spacing = line
        if after: para.space_after = Pt(after)
        for t, o in (p if isinstance(p, list) else [(p, {})]):
            style_run(para.add_run(), o.get("size", size), o.get("color", color), o.get("bold", bold), spc); para.runs[-1].text = t
    return tb


def box(slide, x, y, w, h, fill, radius=0.035, shape=MSO_SHAPE.ROUNDED_RECTANGLE):
    s = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid(); s.fill.fore_color.rgb = RGBColor.from_string(fill); s.line.fill.background(); s.shadow.inherit = False
    if shape == MSO_SHAPE.ROUNDED_RECTANGLE: s.adjustments[0] = radius
    return s


def card(slide, x, y, w, h, accent=BLUE):
    box(slide, x, y, w, h, CARD)
    t = box(slide, x, y, w, 0.06, accent, radius=0.5)   # top edge as in the template cards
    return t


def header(slide, kicker, title, n, note="", sub=None):
    text(slide, 0.89, 0.50, 11.9, 0.26, kicker, size=12, color=BLUE, bold=True, spc=100)
    ph = slide.shapes.title; ph.left, ph.top, ph.width, ph.height = Inches(0.89), Inches(0.86), Inches(11.9), Inches(0.5)
    tf = ph.text_frame; tf.clear(); tf.word_wrap = True
    for a in ("margin_left", "margin_right", "margin_top", "margin_bottom"): setattr(tf, a, 0)
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT; style_run(p.add_run(), 28, INK, True); p.runs[0].text = title
    if sub: text(slide, 0.89, 1.48, 11.9, 0.5, sub, size=14, color=BODY)
    text(slide, 0.89, 6.88, 11.9, 0.26, f"{FOOT} · {n}" + (f" · {note}" if note else ""), size=12, color=MUTED)


def pill(slide, x, y, w, label, fill):
    s = box(slide, x, y, w, 0.34, fill, radius=0.5); tf = s.text_frame
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER; style_run(p.add_run(), 12, "FFFFFF", True); p.runs[0].text = label


def new():
    return prs.slides.add_slide(LAYOUT)


def notes(s, t):
    tf = s.notes_slide.notes_text_frame
    if tf is not None: tf.text = t


# ---------------- 1. title
s = new(); box(s, 0, 0, 0.17, 7.5, BLUE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 1.92, 11.9, 0.3, "PBISC 생성 행렬 검증 · HVG 민감도", size=14, color=BLUE, bold=True, spc=150)
ph = s.shapes.title; ph.left, ph.top, ph.width, ph.height = Inches(0.77), Inches(2.5), Inches(11.9), Inches(1.5)
tf = ph.text_frame; tf.clear(); tf.word_wrap = True
for a in ("margin_left", "margin_right", "margin_top", "margin_bottom"): setattr(tf, a, 0)
p = tf.paragraphs[0]; style_run(p.add_run(), 40, INK, True); p.runs[0].text = "HVG 1,000으로 다시 한 Harmony 통합:\n실제 패혈증 세포와 PBISC 생성 세포의 섞임"
text(s, 0.89, 4.15, 11.9, 0.7, "Ye 2026 패혈증 실제 단일세포 + PBISC K = 1,000 + K = 4,000 생성 세포를 HVG 1,000개로 통합하고, "
     "같은 세포·같은 설정의 HVG 2,000 결과와 iLISI를 비교", size=18, color=BODY)
pill(s, 0.89, 5.1, 1.01, "K = 1,000", BLUE); pill(s, 2.06, 5.1, 1.01, "K = 4,000", RED); pill(s, 3.23, 5.1, 1.28, "Real (public)", GREY)
text(s, 0.89, 6.82, 11.9, 0.26, "Jong Hwi Kim · Seoul National University · 2026년 10월 20일", size=12, color=MUTED)
notes(s, "이 보고서는 통합 보고서 p.18 공동 임베딩을 HVG 개수만 바꿔 다시 계산한 결과입니다.")

# ---------------- 2. question / setup (three cards)
s = new(); header(s, "질문과 설정", "HVG를 1,000개로 줄이면 실제–생성 세포가 더 잘 섞이는가?", 2)
cards = [("1", "데이터", "실제 133명 × 400 cells = 53,200\nK = 1,000 / K = 4,000: 벌크 단독 122 샘플 × 400 cells (5 draw 합산) = 각 48,800\n공통 패널 유전자 5,853개"),
         ("2", "파이프라인", "log1p CP10k → 실제 세포에서 seurat_v3 HVG 선택 (1,000 또는 2,000) → 실제 평균·SD로 표준화 → "
                          "실제로 PCA 30 적합, 생성 세포는 투영 → Harmony (batch = source)"),
         ("3", "지표", "iLISI (perplexity 30, k = 90)의 중앙값\n· 3개 소스 전체 (최대 3)\n· 실제 vs K별 쌍 (최대 2)\nHarmony 전(PCA 투영)과 후를 모두 계산")]
for i, (num, h, b) in enumerate(cards):
    x = 0.89 + i * 3.93; card(s, x, 1.95, 3.7, 3.9)
    text(s, x + 0.28, 2.25, 3.2, 0.5, num, size=32, color=BLUE, bold=True)
    text(s, x + 0.28, 2.85, 3.2, 0.45, h, size=18, color=INK, bold=True)
    text(s, x + 0.28, 3.4, 3.2, 2.3, b.split("\n"), size=14, color=BODY, line=1.15)
text(s, 0.89, 6.1, 11.9, 0.5, "두 설정은 같은 세포, 같은 seed, 같은 Harmony 설정(max_iter 20)을 사용. 바뀐 것은 HVG 개수뿐.", size=14, color=BODY)
notes(s, "c_joint.py와 같은 방법. 실제 세포는 patient_pair_ilisi_20261018 추출본(환자당 800)에서 400개, 생성 세포는 샘플당 800에서 400개를 뽑음.")

# ---------------- 3. UMAP HVG 1,000
s = new(); header(s, "결과 · HVG 1,000 공동 UMAP", "HVG 1,000: 같은 계통에 놓이지만 완전히 섞이지는 않는다", 3,
                  "소스당 30,000 cells 표시", sub=None)
pic = s.shapes.add_picture(FIG + "/umap_hvg1000_sources.png", Inches(0.89), Inches(1.75), width=Inches(11.55))
text(s, 0.89, 5.35, 11.55, 1.2, [[("읽는 법  ", {"bold": True, "color": INK}), ("회색 배경 = 다른 소스. Harmony 후 T/NK, 단핵구, B, 혈소판 섬은 세 소스가 공유하지만, "
     "섬 안에서 생성 세포가 몰리는 하위 영역이 따로 있음(쌍 iLISI 1.24–1.25 / 2).", {})]], size=14, color=BODY, line=1.15)

# ---------------- 4. UMAP HVG 1,000 vs 2,000 + cell type
s = new(); header(s, "결과 · HVG 1,000 vs 2,000", "UMAP 모양은 달라져도 계통 구성과 섞임 정도는 거의 같다", 4, "생성 세포 세포형 = 실제 라벨 kNN 전이 (k = 30)")
s.shapes.add_picture(FIG + "/umap_hvg1000_vs_2000.png", Inches(0.89), Inches(1.5), height=Inches(5.25))
x0 = 8.3
for i, (h, b) in enumerate([("같은 계통 구조", "두 설정 모두 CD4 T, CD8 T/NK/NKT, 단핵구, B/형질세포, 기타(혈소판 등)가 같은 섬으로 나뉨."),
                            ("다른 좌표", "UMAP 배치는 HVG에 따라 바뀌지만 이것은 시각화 차이이며, iLISI는 Harmony 임베딩(30 PC)에서 계산."),
                            ("소스 섞임", "섬 안에서 생성 세포가 몰리는 하위 영역이 두 설정 모두에서 남음 (iLISI 차이 ≤ 0.03) → HVG 수는 격차를 설명하지 못함.")]):
    y = 1.6 + i * 1.72; card(s, x0, y, 4.15, 1.55)
    text(s, x0 + 0.25, y + 0.2, 3.7, 0.35, h, size=16, color=INK, bold=True)
    text(s, x0 + 0.25, y + 0.6, 3.7, 0.9, b, size=13, color=BODY, line=1.1)

# ---------------- 5. iLISI bar
s = new(); header(s, "결과 · iLISI", "HVG 1,000 vs 2,000 차이 ≤ 0.03, 실제–실제 기준과는 큰 격차", 5,
                  "막대 = Harmony 후, 흰 원 = Harmony 전 · 점선 = results/check 실제–실제 기준(HVG 2,000)")
s.shapes.add_picture(FIG + "/ilisi_bar_deck.png", Inches(0.89), Inches(1.55), width=Inches(11.55))
d1, d2 = v(1000, "post", "real vs K1k") - v(2000, "post", "real vs K1k"), v(1000, "post", "real vs K4k") - v(2000, "post", "real vs K4k")
text(s, 0.89, 5.75, 11.55, 0.9, [[("Harmony 후 HVG 1,000 − 2,000: ", {"bold": True, "color": INK}),
     (f"실제 vs K1k {d1:+.2f}, 실제 vs K4k {d2:+.2f}, 3개 소스 {v(1000, 'post', 'real + K1k + K4k (3 sources)') - v(2000, 'post', 'real + K1k + K4k (3 sources)'):+.2f}. ".replace("-", "−") +
      f"Harmony 전에는 HVG 1,000이 +0.03 높음 (1.10 vs 1.06). 같은 연구 실제 세포 1.74, 다른 연구 1.51, PBISC 형식 실제 1.44보다 모두 낮음.", {})]],
     size=14, color=BODY, line=1.15)

# ---------------- 6. numbers table
s = new(); header(s, "수치", "iLISI 전체 값", 6, "results/ilisi_hvg1000_vs_hvg2000.csv")
rows = [("비교", "최대", "HVG 1,000\nHarmony 전", "HVG 1,000\nHarmony 후", "HVG 2,000\nHarmony 전", "HVG 2,000\nHarmony 후", "후 차이\n(1,000 − 2,000)")]
for c, lab, mx in [("real + K1k + K4k (3 sources)", "실제 + K1k + K4k (3개 소스)", 3), ("real vs K1k", "실제 vs K = 1,000", 2), ("real vs K4k", "실제 vs K = 4,000", 2)]:
    rows.append((lab, str(mx), f"{v(1000, 'pre', c):.2f}", f"{v(1000, 'post', c):.2f}", f"{v(2000, 'pre', c):.2f}", f"{v(2000, 'post', c):.2f}",
                 f"{v(1000, 'post', c) - v(2000, 'post', c):+.2f}".replace("-", "−")))
tbl = s.shapes.add_table(len(rows), 7, Inches(0.9), Inches(1.75), Inches(11.54), Inches(2.4)).table
widths = [3.3, 0.9, 1.45, 1.45, 1.45, 1.45, 1.54]
for j, w in enumerate(widths): tbl.columns[j].width = Inches(w)
for i, row in enumerate(rows):
    for j, val in enumerate(row):
        cell = tbl.cell(i, j); cell.fill.solid(); cell.fill.fore_color.rgb = RGBColor.from_string("E4ECF3" if i == 0 else ("FFFFFF" if i % 2 else CARD))
        cell.margin_left = cell.margin_right = Inches(0.1); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        tf = cell.text_frame; tf.clear(); lines = val.split("\n")
        for k, ln in enumerate(lines):
            p = tf.paragraphs[0] if k == 0 else tf.add_paragraph(); p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
            style_run(p.add_run(), 13 if i == 0 else 14, INK if (i == 0 or j in (3, 5)) else BODY, i == 0 or j in (3, 5)); p.runs[0].text = ln
text(s, 0.89, 4.5, 11.55, 1.9, [
    [("3개 소스 값이 높은 이유  ", {"bold": True, "color": INK}), ("K1k와 K4k가 서로 잘 섞이기 때문(같은 모델, 같은 벌크). 실제–생성 섞임은 쌍 값으로 판단.", {})],
    [("p.18과의 차이  ", {"bold": True, "color": INK}), ("통합 보고서 p.18 (HVG 2,000)은 1.22 / 1.22. 이번 재계산 1.27 / 1.26은 생성 세포 부분집합(샘플당 400)이 달라서 생긴 차이로, 같은 범위.", {})]],
    size=14, color=BODY, line=1.15, after=8)

# ---------------- 7. conclusion
s = new(); header(s, "결론", "HVG 수는 병목이 아니다: 실제–생성 격차는 HVG 선택과 무관", 7)
for i, (num, h, b) in enumerate([("1", "HVG 1,000 ≈ HVG 2,000", "Harmony 후 쌍 iLISI 1.24–1.25 vs 1.26–1.27. 차이 −0.01 ~ −0.03으로 무시할 수준."),
                                 ("2", "K = 1,000 ≈ K = 4,000", "두 HVG 설정 모두 K 간 차이 ≤ 0.01. K 증가는 실제 세포와의 섞임을 바꾸지 않음."),
                                 ("3", "격차는 그대로", "같은 연구 실제 1.74, 다른 연구 1.51과의 차이 0.24–0.50 유지. 기존 HVG/PCA 민감도(V0–V3) 결론과 일치.")]):
    x = 0.89 + i * 3.93; card(s, x, 1.75, 3.7, 3.4)
    text(s, x + 0.28, 2.05, 3.2, 0.5, num, size=32, color=BLUE, bold=True)
    text(s, x + 0.28, 2.65, 3.2, 0.7, h, size=18, color=INK, bold=True)
    text(s, x + 0.28, 3.3, 3.2, 1.8, b, size=14, color=BODY, line=1.15)
text(s, 0.89, 5.45, 11.55, 1.2, [[("다음 단계  ", {"bold": True, "color": INK}),
     ("남은 격차는 유전자 수준 차이(단핵구의 T/NK 유전자 동시 검출, 스트레스 유전자 감소 등)로 추적. HVG는 기존 2,000개를 유지해도 결론이 같음.", {})]],
     size=14, color=BODY, line=1.15)

prs.save(DST); print("saved", DST, len(prs.slides._sldIdLst), "slides")
