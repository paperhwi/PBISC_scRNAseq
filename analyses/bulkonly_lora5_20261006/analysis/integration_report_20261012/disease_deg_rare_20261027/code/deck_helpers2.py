"""Same template/layout/styles as ilisi_floor_20261021/code/build_deck_ko.py and three_model_20261022/code/build_deck_s3.py.
The original template file is not reachable from this session, so an equivalent 13.333 x 7.5 deck is built from a blank
presentation with the identical palette, fonts, header/footer geometry and helper API."""
import numpy as np
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree
BLUE, RED, GREY, INK, BODY, MUTED, CARD = "2B83BA", "D7191C", "7F7F7F", "1F2933", "3E4C59", "7B8794", "F5F7FA"
GREEN, ORANGE, PURPLE = "1A9641", "E08214", "7B3294"
KO = "Malgun Gothic"; FOOT = "PBISC 질병 DEG · 희귀 집단 발현 검증"; FIG = "figs"
prs = Presentation(); prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
LAYOUT = prs.slide_layouts[6]   # blank


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


def title_box(slide, txt, y=0.86, size=28):
    tb = slide.shapes.add_textbox(Inches(0.89), Inches(y), Inches(11.9), Inches(0.55))
    tf = tb.text_frame; tf.word_wrap = True
    for a in ("margin_left", "margin_right", "margin_top", "margin_bottom"): setattr(tf, a, 0)
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT; style_run(p.add_run(), size, INK, True); p.runs[0].text = txt
    return tb


def header(slide, kicker, title, n, note="", sub=None):
    text(slide, 0.89, 0.50, 11.9, 0.26, kicker, size=12, color=BLUE, bold=True, spc=100)
    title_box(slide, title)
    if sub: text(slide, 0.89, 1.48, 11.9, 0.5, sub, size=14, color=BODY)
    text(slide, 0.89, 6.88, 11.9, 0.26, f"{FOOT} · {n}" + (f" · {note}" if note else ""), size=12, color=MUTED)


def pill(slide, x, y, w, label, fill):
    s = box(slide, x, y, w, 0.34, fill, radius=0.5); tf = s.text_frame
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER; style_run(p.add_run(), 12, "FFFFFF", True); p.runs[0].text = label


def new(): return prs.slides.add_slide(LAYOUT)

NOTE_LOG = []

def notes(s, t):
    ns = s.notes_slide; tf = ns.notes_text_frame; tf.clear()
    for i, line in enumerate(t.strip().split("\n")):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); r = para.add_run(); r.text = line
        r.font.size = Pt(12); rPr = r._r.get_or_add_rPr()
        e = rPr.find(qn("a:ea"))
        if e is None: e = etree.SubElement(rPr, qn("a:ea"))
        e.set("typeface", KO)
    NOTE_LOG.append(t.strip())


def side_cards(s, x0, y0, w, items, h=1.55, gap=0.17, accent=BLUE, body=12.5, head=15):
    for i, (hh, b) in enumerate(items):
        y = y0 + i * (h + gap); card(s, x0, y, w, h, accent if not isinstance(accent, list) else accent[i])
        text(s, x0 + 0.2, y + 0.18, w - 0.4, 0.35, hh, size=head, color=INK, bold=True)
        text(s, x0 + 0.2, y + 0.55, w - 0.4, h - 0.6, b, size=body, color=BODY, line=1.08)


def three_cards(s, y, h, items, numsize=32, accents=None):
    for i, (num, hh, b) in enumerate(items):
        x = 0.89 + i * 3.93; card(s, x, y, 3.7, h, (accents[i] if accents else BLUE))
        text(s, x + 0.28, y + 0.3, 3.2, 0.5, num, size=numsize, color=(accents[i] if accents else BLUE), bold=True)
        text(s, x + 0.28, y + 0.9, 3.2, 0.45, hh, size=18, color=INK, bold=True)
        text(s, x + 0.28, y + 1.45, 3.2, h - 1.55, b if isinstance(b, list) else b.split("\n"), size=13.5, color=BODY, line=1.12)


def table(s, rows, x, y, widths, h_row=0.36, bold_cols=(), font=13, head_font=12, hl_rows=(), warn_rows=()):
    tbl = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(sum(widths)), Inches(h_row * len(rows))).table
    tbl.first_row = True; tbl.horz_banding = False
    for j, w in enumerate(widths): tbl.columns[j].width = Inches(w)
    for i, row in enumerate(rows):
        tbl.rows[i].height = Inches(h_row if i else h_row * 1.35)
        for j, val in enumerate(row):
            cell = tbl.cell(i, j); cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor.from_string(
                "E4ECF3" if i == 0 else ("FBE3E4" if i in warn_rows else ("EAF3FA" if i in hl_rows else ("FFFFFF" if i % 2 else CARD))))
            cell.margin_left = cell.margin_right = Inches(0.08); cell.margin_top = cell.margin_bottom = Inches(0.02); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            tf = cell.text_frame; tf.clear()
            for k, ln in enumerate(str(val).split("\n")):
                p = tf.paragraphs[0] if k == 0 else tf.add_paragraph(); p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
                b = i == 0 or j in bold_cols or i in hl_rows or i in warn_rows
                style_run(p.add_run(), head_font if i == 0 else font, INK if b else BODY, b); p.runs[0].text = ln
    return tbl


def pic(s, name, x, y, w=None, h=None):
    return s.shapes.add_picture(f"{FIG}/{name}.png", Inches(x), Inches(y), width=Inches(w) if w else None, height=Inches(h) if h else None)
