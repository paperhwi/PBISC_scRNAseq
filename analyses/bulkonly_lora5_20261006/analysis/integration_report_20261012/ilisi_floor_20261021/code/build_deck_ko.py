"""Korean deck for the floor read-out validation (2026-10-21), in the style of 'PBISC HVG1000 Harmony iLISI 보고서.pptx' / 'PBISC 저발현 유전자 절벽 보고서.pptx'
(master and 'Title and body' layout of 'PBISC K1000 vs K4000 Integration Report_latest.pptx'; all template slides dropped)."""
import os, json, numpy as np, pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree

B = "/mnt/user-data/outputs/ilisi_floor_20261021"; FIG = B + "/figs/deck"; RB = B + "/results/big"
UP = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/integration_report_20261012"
TEMPLATE = "/mnt/user-data/uploads/PBISC/PBISC K1000 vs K4000 Integration Report_latest.pptx"
DST = "/mnt/user-data/outputs/deck/PBISC floor read-out 검증 보고서.pptx"
BLUE, RED, GREY, INK, BODY, MUTED, CARD = "2B83BA", "D7191C", "7F7F7F", "1F2933", "3E4C59", "7B8794", "F5F7FA"
KO = "Malgun Gothic"; FOOT = "PBISC floor read-out 보고서"
f2 = lambda x: f"{x:.2f}"; mn = lambda s: s.replace("-", "−")

# ---------------- numbers
IL = pd.read_csv(f"{UP}/ilisi_floor_20261021/results/ilisi_all.csv"); IL0 = IL[IL.variant == "V0"].set_index(["cohort", "source"])
LIN = pd.read_csv(f"{UP}/ilisi_floor_20261021/results/ilisi_post_by_lineage_V0.csv").set_index(["cohort", "source"])
DO = pd.read_csv(f"{UP}/ilisi_floor_20261021/results/covid_donor_overlap_ilisi.csv")
BJ = pd.read_csv(f"{RB}/big_joint_summary.csv").set_index(["cohort", "source"])
LR = pd.concat([pd.read_csv(f"{RB}/{c}_{s}_leiden_res.csv") for c in ("sepsis", "covid") for s in ("realU", "floor", "prev4k")]).set_index(["cohort", "source", "resolution"])
MG = pd.read_csv(f"{UP}/lowexpr_floor_20261021/results/m_gate_values.csv").set_index(["cohort", "source"])
RS = pd.read_csv(f"{RB}/rare_summary.csv")
PB = pd.concat([pd.read_csv(f"{RB}/{c}_tx_pb.csv").assign(cohort=c) for c in ("sepsis", "covid")]).groupby(["cohort", "gene_set", "source"]).pearson.median()
CELL = pd.concat([pd.read_csv(f"{RB}/{c}_tx_cell.csv.gz").assign(cohort=c) for c in ("sepsis", "covid")]).groupby(["cohort", "source"]).n_genes.median()
GN = {}
for c in ("sepsis", "covid"):
    G = pd.read_csv(f"{RB}/{c}_tx_gene.csv.gz"); ref = G[G.source == "ref"].set_index("gene").mean_cp10k
    for s in ("realT1f", "floor", "prev4k"):
        d = G[G.source == s].set_index("gene"); r = d.mean_cp10k / ref.reindex(d.index)
        GN[(c, s)] = r.groupby(pd.cut(d.real_rank, [0, 500, 1500, 2000, 3000, 4000, 10 ** 5]), observed=True).median().to_numpy()
        GN[(c, s, "zero")] = int((d.mean_cp10k == 0).sum())
    GN[(c, "ref", "zero")] = int((ref == 0).sum())
def cat(v):
    return "none" if v == 0 or np.isnan(v) else ("low" if v < 0.5 else ("over" if v > 2 else "ok"))
for k in ("floor", "prev"): RS[f"cat_{k}"] = RS[f"{k}_vs_realB_author"].map(cat)
RC = {(c, k): RS[RS.cohort == c][f"cat_{k}"].value_counts().to_dict() for c in ("sepsis", "covid") for k in ("floor", "prev")}
NS = {c: int((RS.cohort == c).sum()) for c in ("sepsis", "covid")}
r_ = lambda c, s: RS[(RS.cohort == c) & (RS.subset == s)].iloc[0]

prs = Presentation(TEMPLATE); sld = prs.slides._sldIdLst
for s in list(sld): prs.part.drop_rel(s.get(qn("r:id"))); sld.remove(s)
LAYOUT = [l for l in prs.slide_layouts if l.name == "Title and body"][0]


def style_run(r, size, color, bold=False, spc=None):
    r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = RGBColor.from_string(color); r.font.name = "Arial"
    rPr = r._r.get_or_add_rPr()
    for tag in ("a:ea", "a:cs"):
        e = rPr.find(qn(tag))
        if e is None: e = etree.SubElement(rPr, qn(tag))
        e.set("typeface", KO if tag == "a:ea" else "Arial")
    if spc: rPr.set("spc", str(spc))


def text(slide, x, y, w, h, paras, size=14, color=BODY, bold=False, align=PP_ALIGN.LEFT, spc=None, anchor=MSO_ANCHOR.TOP, line=None, after=0):
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
    box(slide, x, y, w, h, CARD); box(slide, x, y, w, 0.06, accent, radius=0.5)


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


def new(): return prs.slides.add_slide(LAYOUT)


NOTES_SP = ('<p:sp xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
            '<p:nvSpPr><p:cNvPr id="3" name="Notes Placeholder 2"/><p:cNvSpPr><a:spLocks noGrp="1"/></p:cNvSpPr><p:nvPr><p:ph type="body" idx="1"/></p:nvPr></p:nvSpPr>'
            '<p:spPr><a:xfrm><a:off x="685800" y="4400550"/><a:ext cx="5486400" cy="3600450"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>'
            '<p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:endParaRPr lang="ko-KR"/></a:p></p:txBody></p:sp>')


def notes(s, t):
    """Speaker notes. The template's notes master has no body placeholder, so one is added when missing (otherwise notes would be dropped)."""
    ns = s.notes_slide
    if ns.notes_placeholder is None:
        ns.shapes._spTree.append(etree.fromstring(NOTES_SP))
    tf = ns.notes_text_frame; tf.clear()
    for i, line in enumerate(t.strip().split("\n")):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); r = para.add_run(); r.text = line
        r.font.size = Pt(12); rPr = r._r.get_or_add_rPr()
        e = rPr.find(qn("a:ea"))
        if e is None: e = etree.SubElement(rPr, qn("a:ea"))
        e.set("typeface", KO)


def side_cards(s, x0, y0, w, items, h=1.55, gap=0.17, accent=BLUE, body=12.5, head=15):
    for i, (hh, b) in enumerate(items):
        y = y0 + i * (h + gap); card(s, x0, y, w, h, accent if not isinstance(accent, list) else accent[i])
        text(s, x0 + 0.2, y + 0.18, w - 0.4, 0.35, hh, size=head, color=INK, bold=True)
        text(s, x0 + 0.2, y + 0.55, w - 0.4, h - 0.6, b, size=body, color=BODY, line=1.08)


def three_cards(s, y, h, items, numsize=32):
    for i, (num, hh, b) in enumerate(items):
        x = 0.89 + i * 3.93; card(s, x, y, 3.7, h)
        text(s, x + 0.28, y + 0.3, 3.2, 0.5, num, size=numsize, color=BLUE, bold=True)
        text(s, x + 0.28, y + 0.9, 3.2, 0.45, hh, size=18, color=INK, bold=True)
        text(s, x + 0.28, y + 1.45, 3.2, h - 1.55, b if isinstance(b, list) else b.split("\n"), size=13.5, color=BODY, line=1.12)


def table(s, rows, x, y, widths, h_row=0.36, bold_cols=(), font=13, head_font=12, hl_rows=()):
    tbl = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(sum(widths)), Inches(h_row * len(rows))).table
    for j, w in enumerate(widths): tbl.columns[j].width = Inches(w)
    for i, row in enumerate(rows):
        tbl.rows[i].height = Inches(h_row if i else h_row * 1.35)
        for j, val in enumerate(row):
            cell = tbl.cell(i, j); cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor.from_string("E4ECF3" if i == 0 else ("EAF3FA" if i in hl_rows else ("FFFFFF" if i % 2 else CARD)))
            cell.margin_left = cell.margin_right = Inches(0.08); cell.margin_top = cell.margin_bottom = Inches(0.02); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            tf = cell.text_frame; tf.clear()
            for k, ln in enumerate(str(val).split("\n")):
                p = tf.paragraphs[0] if k == 0 else tf.add_paragraph(); p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
                b = i == 0 or j in bold_cols or i in hl_rows
                style_run(p.add_run(), head_font if i == 0 else font, INK if b else BODY, b); p.runs[0].text = ln
    return tbl


def pic(s, name, x, y, w=None, h=None):
    return s.shapes.add_picture(f"{FIG}/{name}.png", Inches(x), Inches(y), width=Inches(w) if w else None, height=Inches(h) if h else None)


iv = lambda c, s, k="iLISI_post": IL0.loc[(c, s), k]
# ================= 1. title
s = new(); box(s, 0, 0, 0.17, 7.5, BLUE, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.89, 1.92, 11.9, 0.3, "PBISC 생성 행렬 검증 · floor read-out", size=14, color=BLUE, bold=True, spc=150)
ph = s.shapes.title; ph.left, ph.top, ph.width, ph.height = Inches(0.77), Inches(2.5), Inches(11.9), Inches(1.5)
tf = ph.text_frame; tf.clear(); tf.word_wrap = True
for a in ("margin_left", "margin_right", "margin_top", "margin_bottom"): setattr(tf, a, 0)
p = tf.paragraphs[0]; style_run(p.add_run(), 40, INK, True); p.runs[0].text = "새 floor read-out 생성 세포 검증:\n유전자 곡선, 전사체, Harmony 통합, 희귀 세포"
text(s, 0.89, 4.15, 11.9, 0.7, "depth256-K4000-floor1e-4 (draw0) 생성 세포를 실제 held-out 세포(환자 단위 절반 분할) 및 이전 K = 4,000 생성 세포와 "
     "유전자 count 곡선, 전사체 특성, Harmony 후 iLISI · cluster matching, 희귀 세포 비율로 비교", size=18, color=BODY)
pill(s, 0.89, 5.1, 1.75, "Floor read-out", BLUE); pill(s, 2.8, 5.1, 1.95, "Previous K = 4,000", RED); pill(s, 4.91, 5.1, 1.55, "Real (held-out)", GREY)
text(s, 0.89, 6.82, 11.9, 0.26, "Jong Hwi Kim · Seoul National University · 2026년 10월 21일", size=12, color=MUTED)
notes(s, "새 inference 파일(floor read-out)을 저발현 절벽 보고서(10/20)와 HVG 1,000 보고서의 방법으로 평가. 패혈증은 성인만, COVID는 adapter 학습 donor 5 샘플 제외.")

# ================= 2. summary
s = new(); header(s, "요약", "절벽은 사라졌고 실제 세포와의 섞임은 실제–실제 수준에 근접", 2)
fs, fc = RC[("sepsis", "floor")].get("ok", 0), RC[("covid", "floor")].get("ok", 0); ps, pc = RC[("sepsis", "prev")].get("ok", 0), RC[("covid", "prev")].get("ok", 0)
three_cards(s, 1.75, 3.75, [
    ("1", "유전자 · 전사체", f"순위 2,000 이후 유전자의 CP10k 비율 0.0–0.1 → 0.6–1.2 (이전 → floor). 저발현 유전자 pseudobulk r {f2(PB[('sepsis','genes ranked 2,001+','prev4k')])} / {f2(PB[('covid','genes ranked 2,001+','prev4k')])} → "
          f"{f2(PB[('sepsis','genes ranked 2,001+','floor')])} / {f2(PB[('covid','genes ranked 2,001+','floor')])}. 대신 세포당 검출 유전자가 실제보다 많음 (과검출)."),
    ("2", "Harmony · iLISI · cluster", f"iLISI (최대 2) floor {f2(iv('sepsis','floor'))} / {f2(iv('covid','floor'))}, 이전 {f2(iv('sepsis','prev4k'))} / {f2(iv('covid','prev4k'))}, "
          f"실제–실제 {f2(iv('sepsis','realU'))} / {f2(iv('covid','realU'))}. Harmony 전에도 floor가 높음. cluster 빈도 r도 실제–실제 범위."),
    ("3", "희귀 세포", f"실제 절반 B 대비 0.5–2배 안에 드는 희귀 subset: 패혈증 {ps} → {fs} / {NS['sepsis']}, COVID {pc} → {fc} / {NS['covid']} (이전 → floor). "
          f"증식 NK · CD4 T, DC1, NK CD56 bright 회복, macrophage는 0.02 → 0.31배로 개선. 증식 · dex 단핵구, γδT, 패혈증 pDC, ILC는 여전히 부족.")])
text(s, 0.89, 5.75, 11.55, 0.9, [[("다음 단계  ", {"bold": True, "color": INK}), ("과검출 보정(threshold 조정, 양방향 gate)과 함께, 남은 희귀 상태(증식 단핵구, dex 단핵구, γδT, ILC)가 "
     "모델 결손인지 mapping 한계인지 marker 기반으로 분리.", {})]], size=14, color=BODY, line=1.15)
notes(s, "모든 수치는 같은 파이프라인. 실제 대조군은 환자 단위로 나눈 다른 절반(real B)이라 생성 세포와 같은 '다른 환자' 상황.")

# ================= 3. setup
s = new(); header(s, "질문과 설정", "새 read-out은 실제 세포와 얼마나 같아졌는가?", 3)
three_cards(s, 1.75, 4.0, [
    ("1", "데이터", f"Floor: draw0, bulk × 256, K = 4,000, 샘플당 4,000 cells, x < 1e-4 → 0 · < 1 → 1 · 그 외 floor(x)\n이전 K = 4,000: Poisson read-out, 5 draw 합산\n"
               f"실제: held-out 세포를 환자 단위로 A(참조) / B(대조)로 분할 (패혈증 성인 20 / 20명, COVID 17 / 18명)"),
    ("2", "파이프라인", "log1p CP10k → 실제 A에서 seurat_v3 HVG 2,000 → A 평균·SD로 표준화 → A로 PCA 30 → 소스별로 A와 쌍 Harmony (batch = source). "
               "라벨: A의 author 라벨을 kNN (k = 30)으로 전이"),
    ("3", "지표", "유전자: 순위별 CP10k 비율, 검출·분산\n전사체: 세포형 pseudobulk r, marker\n통합: iLISI (perplexity 30, k = 90), joint Leiden cluster 빈도\n희귀: 부모 계통 안 비율")])
text(s, 0.89, 6.0, 11.55, 0.6, "세포 수: iLISI 본 분석은 소스당 20,000 cells, UMAP · cluster · 희귀 세포는 소스당 약 60,000 cells (두 번째 독립 subsample, iLISI 재확인 포함).",
     size=14, color=BODY)
notes(s, "code: ilisi_floor_20261021/code (ilisi_prep.py, ilisi_run.py, big_prep.py, big_joint.py, big_txstats.py, big_rare.py, big_leiden_res.py).")

# ================= 4. gene count curve
s = new(); header(s, "결과 1 · 유전자 count 곡선", "순위 2,000 이후의 절벽이 사라지고 실제 곡선을 따라간다", 4,
                  "소스당 약 60,000 cells · results/big/<cohort>_tx_gene.csv.gz")
pic(s, "gene_count_curve", 0.89, 1.5, w=8.35)
a, b = GN[("sepsis", "floor")], GN[("sepsis", "prev4k")]; c_, d_ = GN[("covid", "floor")], GN[("covid", "prev4k")]
side_cards(s, 9.45, 1.55, 3.0, [
    ("순위 2,001 이후", f"순위 2–3k · 3–4k · 4k+ (÷ 실제 A)\n패혈증 이전 {f2(b[3])} · {f2(b[4])} · {f2(b[5])}\n패혈증 floor {f2(a[3])} · {f2(a[4])} · {f2(a[5])}\nCOVID 이전 {f2(d_[3])} · {f2(d_[4])} · {f2(d_[5])}\nCOVID floor {f2(c_[3])} · {f2(c_[4])} · {f2(c_[5])}"),
    ("한 번도 안 나온 유전자", f"이전 {GN[('sepsis','prev4k','zero')]:,} / {GN[('covid','prev4k','zero')]:,}개 → floor {GN[('sepsis','floor','zero')]} / {GN[('covid','floor','zero')]}개 (실제 A {GN[('sepsis','ref','zero')]} / {GN[('covid','ref','zero')]}, depth 맞춘 B {GN[('sepsis','realT1f','zero')]} / {GN[('covid','realT1f','zero')]})"),
    ("남은 결함", f"순위 501–1,500 과다 {f2(a[1])} / {f2(c_[1])}배 그대로 (이전 {f2(b[1])} / {f2(d_[1])}). COVID 순위 4,000 이후는 {f2(c_[5])}로 약간 부족.")], h=1.62, body=11.5)
notes(s, "곡선 = 실제 A 평균 순으로 정렬한 유전자별 평균 CP10k의 running median (51 / 151 genes). 회색 배경 = 순위 > 2,000.")

# ================= 5. cell shape
s = new(); header(s, "결과 1 · 세포 수준 모양", "세포당 검출 유전자는 실제보다 많다: 절벽 대신 과검출", 5, "점선 = depth를 맞춘 실제 B · results/big/<cohort>_tx_cell.csv.gz")
pic(s, "cell_shape", 0.89, 1.5, w=8.45)
side_cards(s, 9.5, 1.55, 2.95, [
    ("검출 유전자 / 세포", f"패혈증: 실제 B(depth 맞춤) {CELL[('sepsis','realT1f')]:.0f} vs floor {CELL[('sepsis','floor')]:.0f} (이전 {CELL[('sepsis','prev4k')]:.0f} vs {CELL[('sepsis','realT1p')]:.0f})\n"
                       f"COVID: {CELL[('covid','realT1f')]:.0f} vs {CELL[('covid','floor')]:.0f} (이전 {CELL[('covid','prev4k')]:.0f} vs {CELL[('covid','realT1p')]:.0f})"),
    ("mean–detection", "평균 0.1–1 UMI 구간에서 floor가 실제보다 더 자주 검출 → threshold 1e-4가 너무 관대 (M1 1.5–1.7배와 일치)."),
    ("mean–variance", "이전 K = 4,000은 평균 0.01–1 구간에서 CV²가 실제보다 큼 (검출당 UMI 과다). floor는 실제 곡선과 겹침.")], h=1.62, body=12)

# ================= 6. M gates table
s = new(); header(s, "수치 · 저발현 gate (M1–M4)", "하한 gate는 모두 통과, 상한은 열려 있다", 6, "lowexpr_floor_20261021/results/m_gate_values.csv · 소스당 30,000 cells, depth 맞춘 실제 대비")
g = lambda c, s_, k: MG.loc[(c, s_), k]
rows = [("Gate (기준)", "패혈증\n이전 K4000", "패혈증\nfloor", "COVID\n이전 K4000", "COVID\nfloor", "판정")]
rows += [("M1 검출 비율, 실제 2–5% (≥ 0.7)", f2(g('sepsis','prev K4000','M1_ratio_2to5pct')), f2(g('sepsis','floor','M1_ratio_2to5pct')), f2(g('covid','prev K4000','M1_ratio_2to5pct')), f2(g('covid','floor','M1_ratio_2to5pct')), "통과 (과검출)"),
         ("M2 검출 비율, 실제 0.5–2% (≥ 0.5)", f2(g('sepsis','prev K4000','M2_ratio_0p5to2pct')), f2(g('sepsis','floor','M2_ratio_0p5to2pct')), f2(g('covid','prev K4000','M2_ratio_0p5to2pct')), f2(g('covid','floor','M2_ratio_0p5to2pct')), "통과 (과검출)"),
         (f"M3 검출당 UMI (1.0–1.4; 실제 {f2(g('sepsis','floor','M3_umi_per_det_real'))} / {f2(g('covid','floor','M3_umi_per_det_real'))})", f2(g('sepsis','prev K4000','M3_umi_per_det_src')), f2(g('sepsis','floor','M3_umi_per_det_src')), f2(g('covid','prev K4000','M3_umi_per_det_src')), f2(g('covid','floor','M3_umi_per_det_src')), "통과"),
         ("M4 CP10k 비율, 순위 500–1,500 (0.85–1.15)", f2(g('sepsis','prev K4000','M4_ratio_rank500_1500')), f2(g('sepsis','floor','M4_ratio_rank500_1500')), f2(g('covid','prev K4000','M4_ratio_rank500_1500')), f2(g('covid','floor','M4_ratio_rank500_1500')), "실패"),
         ("M4 CP10k 비율, 순위 2,000 이후 (≥ 0.5)", f2(g('sepsis','prev K4000','M4_ratio_gt2000_all')), f2(g('sepsis','floor','M4_ratio_gt2000_all')), f2(g('covid','prev K4000','M4_ratio_gt2000_all')), f2(g('covid','floor','M4_ratio_gt2000_all')), "통과")]
table(s, rows, 0.9, 1.7, [4.35, 1.4, 1.4, 1.4, 1.4, 1.6], h_row=0.5, bold_cols=(2, 4))
text(s, 0.89, 5.2, 11.55, 1.4, [[("읽는 법  ", {"bold": True, "color": INK}), ("M1–M3는 하한만 정한 gate라 1.5배 과검출도 통과. 실제 대비 0.7–1.3배 같은 양방향 gate와 '세포당 유전자 수가 depth 맞춘 실제의 ±10%'를 추가해야 함.", {})],
     [("M4 실패  ", {"bold": True, "color": INK}), ("순위 500–1,500 과다는 검출 과다(1.24–1.29배)에서 오며 검출당 UMI 문제가 아님 (저발현 절벽 보고서 결함 2).", {})]], size=14, color=BODY, line=1.15, after=8)

# ================= 7. pseudobulk
s = new(); header(s, "결과 2 · 전사체 특성", "저발현 유전자까지 세포형 평균 발현이 실제에 근접", 7, "세포형 = 실제 author → six-type, 생성 = 고정 scANVI 라벨 · results/big/<cohort>_tx_pb.csv")
pic(s, "pseudobulk_r", 0.89, 1.5, w=11.55)
q = lambda c, gs, s_: f2(PB[(c, gs, s_)])
text(s, 0.89, 5.95, 11.55, 0.8, [[("전체 유전자  ", {"bold": True, "color": INK}), (f"r {q('sepsis','all genes','prev4k')} → {q('sepsis','all genes','floor')} (패혈증), {q('covid','all genes','prev4k')} → {q('covid','all genes','floor')} (COVID); 실제–실제 {q('sepsis','all genes','realT1f')} / {q('covid','all genes','realT1f')}.  ", {}),
     ("순위 2,001 이후  ", {"bold": True, "color": INK}), (f"{q('sepsis','genes ranked 2,001+','prev4k')} → {q('sepsis','genes ranked 2,001+','floor')}, {q('covid','genes ranked 2,001+','prev4k')} → {q('covid','genes ranked 2,001+','floor')} (실제 {q('sepsis','genes ranked 2,001+','realT1f')} / {q('covid','genes ranked 2,001+','realT1f')}).", {})]],
     size=14, color=BODY, line=1.15)

# ================= 8. markers
s = new(); header(s, "결과 2 · marker 유전자", "강한 marker는 세 소스가 같고, 약한 marker가 floor에서 살아난다", 8, "점 크기 = 검출 세포 비율, 색 = 평균 log CP10k (같은 색 척도) · results/big/<cohort>_tx_dot.csv")
pic(s, "marker_dot", 0.89, 1.5, w=11.55)
text(s, 0.89, 5.75, 11.55, 0.9, [[("읽는 법  ", {"bold": True, "color": INK}), ("세 소스 모두 marker가 맞는 계통에서만 켜지고, 강한 marker(CD3E, CD79A, LYZ, NKG7)의 검출 비율과 평균은 셋이 비슷. 차이는 약한 marker: "
     "NK의 NCAM1 검출 실제 A 23% / 20%, 이전 7% / 7%, floor 29% / 34%; DC의 CD1C 44% / 49%, 25% / 34%, 54% / 60% (패혈증 / COVID). floor는 실제 이상으로 검출 (과검출과 일치).", {})]], size=14, color=BODY, line=1.15)

# ================= 9/10. UMAPs
for k, (coh, nm) in enumerate((("sepsis", "패혈증"), ("covid", "COVID")), start=9):
    s = new(); header(s, f"결과 3 · Harmony 통합 · {nm}", f"{nm}: floor 세포는 Harmony 전부터 실제 세포 사이에 놓인다", k, "실제 A 25,000 + 소스 25,000 cells, 무작위 순서로 그림 · iLISI = 소스당 약 60,000 cells 재계산")
    pic(s, f"umap_{coh}", 0.89, 1.5, h=5.0)
    bj = lambda s_, kk: BJ.loc[(coh, s_), kk]
    side_cards(s, 10.5, 1.55, 1.95, [("Harmony 전", f"이전 {f2(bj('prev4k','iLISI_pre'))}\nfloor {f2(bj('floor','iLISI_pre'))}\n실제 B {f2(bj('realU','iLISI_pre'))}"),
                                       ("Harmony 후", f"이전 {f2(bj('prev4k','iLISI_post'))}\nfloor {f2(bj('floor','iLISI_post'))}\n실제 B {f2(bj('realU','iLISI_post'))}"),
                                       ("6-type 일치", f"전이 라벨 = 자기 라벨\n이전 {100*bj('prev4k','six_agree_with_own_label'):.0f}%\nfloor {100*bj('floor','six_agree_with_own_label'):.0f}%\n실제 B {100*bj('realU','six_agree_with_own_label'):.0f}%")],
               h=1.6, body=11.5, head=13.5)
    notes(s, "윗줄: 회색 = 실제 A, 색 = 소스. 섞인 곳은 점무늬, 분리된 곳은 단색. 아랫줄: 세포형 색과 floor joint Leiden cluster.")

# ================= 11. iLISI
s = new(); header(s, "결과 4 · iLISI", "floor의 섞임은 실제–실제 기준에 근접, 이전 K = 4,000과 큰 차이", 11, "막대 = Harmony 후, 흰 원 = Harmony 전 · 소스당 20,000 cells, query 3회 평균 (SD ≤ 0.004)")
pic(s, "ilisi_bar_deck", 0.89, 1.5, w=11.55)
sh = lambda c, s_: (iv(c, s_) - 1) / (iv(c, "realU") - 1)
for i, (h_, b_) in enumerate([(f"{100*sh('sepsis','floor'):.0f}% / {100*sh('covid','floor'):.0f}%", "실제–실제 상한 대비 floor ((iLISI − 1) 기준, 패혈증 / COVID)"),
                              (f"{100*sh('sepsis','prev4k'):.0f}% / {100*sh('covid','prev4k'):.0f}%", "같은 기준의 이전 K = 4,000"),
                              (f"{f2(LIN.loc[('sepsis','floor'),'Myeloid'])} vs {f2(LIN.loc[('sepsis','realU'),'Myeloid'])}", "약한 지점: 패혈증 골수계 (floor vs 실제 B)")]):
    x = 0.89 + i * 3.93; card(s, x, 5.55, 3.7, 1.15)
    text(s, x + 0.25, 5.7, 3.3, 0.45, h_, size=20, color=INK, bold=True); text(s, x + 0.25, 6.18, 3.3, 0.5, b_, size=12, color=BODY)

# ================= 12. iLISI table
s = new(); header(s, "수치 · iLISI", "두 번의 독립 subsample에서 순서가 같다", 12, "results/ilisi_all.csv (20,000) · results/big/big_joint_summary.csv (약 60,000)")
rows = [("소스", "패혈증\nHarmony 전", "패혈증\nHarmony 후", "패혈증 후\n(60k 재계산)", "COVID\nHarmony 전", "COVID\nHarmony 후", "COVID 후\n(60k 재계산)")]
for s_, lab in (("realU", "실제 B (다른 환자, 상한 대조)"), ("realT1", "실제 B, depth 맞춤"), ("floor", "Floor read-out"), ("prev4k", "이전 K = 4,000")):
    rr = [lab]
    for c in ("sepsis", "covid"):
        rr += [f2(iv(c, s_, "iLISI_pre")), f2(iv(c, s_)), f2(BJ.loc[(c, s_), "iLISI_post"]) if (c, s_) in BJ.index else "–"]
    rows.append(tuple(rr))
table(s, rows, 0.9, 1.7, [3.45, 1.35, 1.35, 1.38, 1.35, 1.35, 1.32], h_row=0.48, bold_cols=(2, 5), hl_rows=(3,))
dd = DO.set_index(["source", "group"]) if "group" in DO.columns else None
text(s, 0.89, 4.6, 11.55, 2.0, [[("60k 재계산  ", {"bold": True, "color": INK}), ("세포가 3배 많으면 k = 90 이웃이 더 국소적이라 값이 0.00–0.09 낮지만 floor ≈ 실제 B ≫ 이전 순서는 같음.", {})],
     [("HVG 공동 선택(V1)  ", {"bold": True, "color": INK}), (f"floor {f2(IL.set_index(['variant','cohort','source']).loc[('V1','sepsis','floor'),'iLISI_post'])} / {f2(IL.set_index(['variant','cohort','source']).loc[('V1','covid','floor'),'iLISI_post'])}, "
                                                          f"이전 {f2(IL.set_index(['variant','cohort','source']).loc[('V1','sepsis','prev4k'),'iLISI_post'])} / {f2(IL.set_index(['variant','cohort','source']).loc[('V1','covid','prev4k'),'iLISI_post'])}: 결론 같음.", {})],
     [("donor 겹침 (COVID)  ", {"bold": True, "color": INK}), ("생성 donor 9 / 8명이 실제 A / B에도 있지만, 그 donor의 생성 세포가 더 잘 섞이지 않음 (세포별 iLISI 1.65 / 1.66 vs 겹치지 않는 donor 1.63). 패혈증은 겹침 없음.", {})],
     [("이전 보고서와 비교 불가  ", {"bold": True, "color": INK}), ("통합 보고서 p.15의 K = 4,000 값 1.22 / 1.17은 다른 실제 절반, 5 draw 합산, 다른 세포 수.", {})]], size=13.5, color=BODY, line=1.12, after=6)

# ================= 13. cluster matching figure
s = new(); header(s, "결과 5 · cluster matching", "joint Leiden cluster의 빈도가 실제–실제 범위에서 맞는다", 13, "joint Leiden (resolution 1, k = 15, Harmony 30 PC) · 점 크기 = cluster 크기 · results/big/<cohort>_<source>_leiden_clusters.csv")
pic(s, "cluster_matching", 0.89, 1.5, h=5.2)
lr = lambda c, s_, r, k: LR.loc[(c, s_, r), k]
side_cards(s, 9.35, 1.55, 3.1, [
    ("cluster 빈도 r", f"패혈증: 이전 {f2(lr('sepsis','prev4k',1.0,'freq_pearson'))} · floor {f2(lr('sepsis','floor',1.0,'freq_pearson'))} · 실제 B {f2(lr('sepsis','realU',1.0,'freq_pearson'))}\n"
                     f"COVID: {f2(lr('covid','prev4k',1.0,'freq_pearson'))} · {f2(lr('covid','floor',1.0,'freq_pearson'))} · {f2(lr('covid','realU',1.0,'freq_pearson'))}"),
    ("생성 세포가 없는 실제 cluster", f"해상도 2.5에서 실제 A 세포 중 실제 쪽 cluster(소스 < 20%)에 남는 비율: 이전 {100*lr('sepsis','prev4k',2.5,'ref_in_ref_dominated'):.1f}% / {100*lr('covid','prev4k',2.5,'ref_in_ref_dominated'):.1f}% → "
                         f"floor {100*lr('sepsis','floor',2.5,'ref_in_ref_dominated'):.1f}% / {100*lr('covid','floor',2.5,'ref_in_ref_dominated'):.1f}% (실제 B {100*lr('sepsis','realU',2.5,'ref_in_ref_dominated'):.1f}% / {100*lr('covid','realU',2.5,'ref_in_ref_dominated'):.1f}%)"),
    ("해상도 2.5", f"cluster {lr('sepsis','floor',2.5,'n_clusters')} / {lr('covid','floor',2.5,'n_clusters')}개에서도 floor 세포의 {100*lr('sepsis','floor',2.5,'src_in_shared'):.0f}% / {100*lr('covid','floor',2.5,'src_in_shared'):.0f}%가 공유 cluster "
                    f"(이전 {100*lr('sepsis','prev4k',2.5,'src_in_shared'):.0f}% / {100*lr('covid','prev4k',2.5,'src_in_shared'):.0f}%)")], h=1.62, body=12)

# ================= 14. cluster + label table
s = new(); header(s, "수치 · cluster matching과 라벨 전이", "floor는 cluster 빈도와 라벨 일치에서 실제에 더 가깝다", 14, "공유 cluster = 크기 보정한 소스 비율 0.2–0.8 · JSD = cluster 빈도 분포 차이 (0 = 동일)")
rows = [("지표", "패혈증\n실제 B", "패혈증\n이전", "패혈증\nfloor", "COVID\n실제 B", "COVID\n이전", "COVID\nfloor")]
def rowf(lab, fn): return tuple([lab] + [fn(c, s_) for c in ("sepsis", "covid") for s_ in ("realU", "prev4k", "floor")])
rows += [rowf("cluster 빈도 Pearson r (res 1)", lambda c, s_: f2(lr(c, s_, 1.0, "freq_pearson"))),
         rowf("cluster 빈도 JSD (res 1)", lambda c, s_: f"{lr(c, s_, 1.0, 'JSD'):.3f}"),
         rowf("cluster 빈도 Pearson r (res 2.5)", lambda c, s_: f2(lr(c, s_, 2.5, "freq_pearson"))),
         rowf("공유 cluster 안 소스 세포 (res 2.5)", lambda c, s_: f"{100*lr(c, s_, 2.5, 'src_in_shared'):.0f}%"),
         rowf("실제만 있는 cluster의 실제 A 세포 (res 2.5)", lambda c, s_: f"{100*lr(c, s_, 2.5, 'ref_in_ref_dominated'):.1f}%"),
         rowf("6-type 일치 (전이 = 자기 라벨)", lambda c, s_: f"{100*BJ.loc[(c, s_), 'six_agree_with_own_label']:.1f}%"),
         rowf("fine 라벨 전이 신뢰도 (중앙값)", lambda c, s_: f2(BJ.loc[(c, s_), "conf_fine_median"]))]
table(s, rows, 0.9, 1.7, [4.15, 1.23, 1.23, 1.23, 1.23, 1.23, 1.23], h_row=0.44, bold_cols=(3, 6))
text(s, 0.89, 5.6, 11.55, 1.1, [[("읽는 법  ", {"bold": True, "color": INK}), (f"COVID 실제 B의 r이 낮은 것은 환자 구성 차이(A와 B는 다른 환자) 때문. 실제 B의 fine 라벨 전이 정확도 {100*BJ.loc[('sepsis','realU'),'acc_fine_vs_true']:.0f}% / {100*BJ.loc[('covid','realU'),'acc_fine_vs_true']:.0f}%가 "
     "이 전이 경로의 상한.", {})]], size=14, color=BODY, line=1.15)

# ================= 15. rare figure
s = new(); header(s, "결과 6 · 희귀 세포", "희귀 subset의 절반 이상이 실제 비율 2배 안으로 돌아왔다", 15, "부모 계통 안 비율, pooled · 정의 = integration_report code/f_rare_parent.py (p.23) · results/big/rare_summary.csv")
pic(s, "rare_parent_share", 0.89, 1.5, w=8.7)
side_cards(s, 9.8, 1.55, 2.65, [
    ("실제 B 대비 0.5–2배", f"패혈증 {ps} → {fs} / {NS['sepsis']}\nCOVID {pc} → {fc} / {NS['covid']}\n(이전 → floor)"),
    ("0.5배 미만", f"패혈증 {RC[('sepsis','prev')].get('low',0)+RC[('sepsis','prev')].get('none',0)} → {RC[('sepsis','floor')].get('low',0)+RC[('sepsis','floor')].get('none',0)}\nCOVID {RC[('covid','prev')].get('low',0)+RC[('covid','prev')].get('none',0)} → {RC[('covid','floor')].get('low',0)+RC[('covid','floor')].get('none',0)}"),
    ("mapping 한계", "실제 B도 전이하면 증식 단핵구 1.08 → 0.18%, ILC 0.58 → 0%로 줄어듦 (연회색 막대)")], h=1.6, body=12, head=14)
notes(s, "생성 세포 라벨 = 실제 A author 라벨의 kNN 전이. 진한 회색 = 실제 B author 라벨(정답), 연회색 = 같은 경로로 전이한 실제 B(경로 상한).")

# ================= 16. rare verdict table
s = new(); header(s, "결과 6 · 희귀 세포 판정", "회복된 상태와 여전히 부족한 상태", 16, "비율 = 생성 / 실제 B (author 라벨), 부모 계통 안 % · n = floor에서 전이된 세포 수 · AUROC = 전이된 세포의 marker 점수")
def lines(c, cats):
    d = RS[(RS.cohort == c) & (RS.cat_floor.isin(cats))].sort_values("floor_vs_realB_author", ascending=False)
    return [f"{x.subset}: {f2(x.prev_vs_realB_author)} → {f2(x.floor_vs_realB_author)} (n {int(x.n_floor)})" for x in d.itertuples()]
cols = [("패혈증 · 회복 (0.5–2배)", lines("sepsis", ["ok"]), BLUE), ("패혈증 · 부족 / 과다", lines("sepsis", ["low", "none", "over"]), RED),
        ("COVID · 회복 (0.5–2배)", lines("covid", ["ok"]), BLUE), ("COVID · 부족 / 과다", lines("covid", ["low", "none", "over"]), RED)]
for i, (h_, ls, acc) in enumerate(cols):
    x = 0.89 + i * 2.93; card(s, x, 1.7, 2.78, 3.85, acc)
    text(s, x + 0.18, 1.88, 2.5, 0.35, h_, size=13.5, color=INK, bold=True)
    text(s, x + 0.18, 2.3, 2.5, 3.2, ls, size=10.5, color=BODY, line=1.05, after=3)
au = RS[RS.n_floor >= 10].auroc_floor
text(s, 0.89, 5.75, 11.55, 1.0, [[("숫자 = 이전 → floor (실제 B 대비 배수).  ", {"bold": True, "color": INK}),
     (f"floor에서 전이된 희귀 세포는 해당 marker를 실제처럼 가짐 (n ≥ 10 subset의 marker AUROC 중앙값 {f2(au.median())}, 범위 {f2(au.min())}–{f2(au.max())}). "
      "HSC/HSPC는 두 생성 세트 모두 2–4배 과다. ILC는 실제 B도 전이로는 0이라 판정 불가.", {})]], size=12.5, color=BODY, line=1.12)

# ================= 17. caveats
s = new(); header(s, "해석과 주의", "이 숫자가 보여주지 않는 것", 17)
three_cards(s, 1.75, 3.6, [
    ("1", "원인이 하나가 아님", "floor와 이전 K = 4,000은 입력 depth(× 256), read-out 규칙, draw 수(1 vs 5), 샘플당 세포 수가 동시에 다름 → 개선을 read-out 하나로 돌릴 수 없음."),
    ("2", "섞임 ≠ 충실도", "iLISI와 cluster 빈도는 전역 지표. 세포당 검출 유전자 약 1.35배 과검출이 있어도 높게 나옴. DEG, 질병 효과 크기는 별도 확인 필요 (M6 미시험)."),
    ("3", "희귀 세포 판정의 한계", "라벨은 kNN 전이라 증식 단핵구, γδT, ILC처럼 실제에서도 전이로 줄어드는 상태는 모델 결손과 mapping 한계가 섞임. subsample 1회, 일부 subset은 n < 20.")])
text(s, 0.89, 5.6, 11.55, 0.4, "Harmony seed 1개(0), 환자 단위 resampling 없음. 실제 대조군은 다른 환자이므로 생성 세포와 같은 조건.", size=13, color=BODY)

# ================= 18. conclusion
s = new(); header(s, "결론", "절벽과 통합 격차는 해결, 과검출과 일부 희귀 상태는 남음", 18)
three_cards(s, 1.75, 3.5, [
    ("1", "유전자 · 전사체", f"순위 2,000 이후 절벽 해소, 저발현 유전자 pseudobulk r 0.93–0.94 (이전 0.46–0.53). 대가는 세포당 검출 유전자 과다."),
    ("2", "통합", f"Harmony 후 iLISI {f2(iv('sepsis','floor'))} / {f2(iv('covid','floor'))}, 실제–실제 {f2(iv('sepsis','realU'))} / {f2(iv('covid','realU'))} (이전 {f2(iv('sepsis','prev4k'))} / {f2(iv('covid','prev4k'))}). cluster 빈도와 라벨 일치도 실제 범위."),
    ("3", "희귀 세포", f"0.5–2배 안 subset 패혈증 {ps} → {fs}, COVID {pc} → {fc}. 증식 NK · CD4 T, DC1, NK CD56 bright 회복; macrophage 개선(0.31배); 증식 · dex 단핵구, γδT, 패혈증 pDC, ILC 미회복.")])
text(s, 0.89, 5.5, 11.55, 1.2, [[("다음 단계  ", {"bold": True, "color": INK}), ("① threshold sweep (1e-4 … 0.5)으로 과검출 보정 후 같은 패널 재실행  ② 양방향 gate (검출 0.7–1.3배, 세포당 유전자 ±10%)  "
     "③ 미회복 희귀 상태는 marker 기반(전이 없이) 판정으로 모델 결손 여부 분리  ④ 패혈증 골수계 iLISI 저하 원인 분석.", {})]], size=14, color=BODY, line=1.2)

prs.save(DST); print("saved", DST, len(prs.slides._sldIdLst), "slides")
