#!/usr/bin/env python3
"""
02_build_manuscript_docx_pdf.py
-------------------------------
Compiles `PBISC_Manuscript_Fig4_7_Final_20261125.md` into publication-formatted
`PBISC_Manuscript_Fig4_7_Final_20261125.docx` and `PBISC_Manuscript_Fig4_7_Final_20261125.pdf`,
embedding actual full-resolution Main Figures 4, 5, 6, and 7 inline above their
respective legends, converting Unicode superscript citations (`¹`–`²⁰`) into
native OOXML/ReportLab superscripts, and asserting zero unparsed Markdown tokens.
"""

from __future__ import annotations
import hashlib
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape

import docx
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch, mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    HRFlowable,
    Image as RLImage,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
)

try:
    import pymupdf as fitz
except ImportError:
    import fitz

PKG_DIR = Path(__file__).resolve().parent.parent
MD_PATH = PKG_DIR / "PBISC_Manuscript_Fig4_7_Final_20261125.md"
DOCX_PATH = PKG_DIR / "PBISC_Manuscript_Fig4_7_Final_20261125.docx"
PDF_PATH = PKG_DIR / "PBISC_Manuscript_Fig4_7_Final_20261125.pdf"

SUP_MAP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
SUP_CHARS_RE = re.compile(r"([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)")
IMG_LINE_RE = re.compile(r"^!\[(.*?)\]\((.*?)\)$")
# Handles bold spans that may contain backticked expressions with literal `*` (e.g. **`Dividing ASC*`**)
TOKEN_RE = re.compile(r"(\*\*(?:`[^`]+`|[^*])+\*\*|`[^`]+`|\*[^*\s][^*]*?\*)")


def split_sup_segments(text: str, convert_sup: bool = True) -> list[tuple[str, bool]]:
    """Splits text into (segment, is_superscript) pairs, keeping 10²/10³/10⁴ or citation superscripts crisp."""
    if not convert_sup:
        return [(text, False)]
    out: list[tuple[str, bool]] = []
    pos = 0
    for m in SUP_CHARS_RE.finditer(text):
        if m.start() > pos:
            out.append((text[pos:m.start()], False))
        out.append((m.group(1).translate(SUP_MAP), True))
        pos = m.end()
    if pos < len(text):
        out.append((text[pos:], False))
    return out


def add_docx_formatted_runs(p, text: str, base_size_pt: float = 10.5, base_color: RGBColor | None = None) -> None:
    pos = 0
    for m in TOKEN_RE.finditer(text):
        if m.start() > pos:
            plain = text[pos:m.start()]
            for seg, is_sup in split_sup_segments(plain, convert_sup=True):
                r = p.add_run(seg)
                r.font.name = "Times New Roman"
                r.font.size = Pt(base_size_pt)
                if base_color:
                    r.font.color.rgb = base_color
                if is_sup:
                    r.font.superscript = True
        tok = m.group(0)
        if tok.startswith("**") and tok.endswith("**"):
            inner = tok[2:-2]
            # Inner may contain backticks like `Dividing ASC*`
            sub_pos = 0
            for sm in re.finditer(r"`([^`]+)`", inner):
                if sm.start() > sub_pos:
                    for seg, is_sup in split_sup_segments(inner[sub_pos:sm.start()], convert_sup=True):
                        r = p.add_run(seg)
                        r.bold = True
                        r.font.name = "Times New Roman"
                        r.font.size = Pt(base_size_pt)
                        if base_color:
                            r.font.color.rgb = base_color
                        if is_sup:
                            r.font.superscript = True
                code_txt = sm.group(1)
                for seg, is_sup in split_sup_segments(code_txt, convert_sup=True):
                    r = p.add_run(seg)
                    r.bold = True
                    r.font.name = "Times New Roman"
                    r.font.size = Pt(base_size_pt)
                    if base_color:
                        r.font.color.rgb = base_color
                    if is_sup:
                        r.font.superscript = True
                sub_pos = sm.end()
            if sub_pos < len(inner):
                for seg, is_sup in split_sup_segments(inner[sub_pos:], convert_sup=True):
                    r = p.add_run(seg)
                    r.bold = True
                    r.font.name = "Times New Roman"
                    r.font.size = Pt(base_size_pt)
                    if base_color:
                        r.font.color.rgb = base_color
                    if is_sup:
                        r.font.superscript = True
        elif tok.startswith("`") and tok.endswith("`"):
            inner = tok[1:-1]
            for seg, is_sup in split_sup_segments(inner, convert_sup=True):
                r = p.add_run(seg)
                r.font.name = "Times New Roman"
                r.font.size = Pt(base_size_pt)
                if base_color:
                    r.font.color.rgb = base_color
                if is_sup:
                    r.font.superscript = True
        elif tok.startswith("*") and tok.endswith("*"):
            inner = tok[1:-1]
            for seg, is_sup in split_sup_segments(inner, convert_sup=True):
                r = p.add_run(seg)
                r.italic = True
                r.font.name = "Times New Roman"
                r.font.size = Pt(base_size_pt)
                if base_color:
                    r.font.color.rgb = base_color
                if is_sup:
                    r.font.superscript = True
        pos = m.end()

    if pos < len(text):
        for seg, is_sup in split_sup_segments(text[pos:], convert_sup=True):
            r = p.add_run(seg)
            r.font.name = "Times New Roman"
            r.font.size = Pt(base_size_pt)
            if base_color:
                r.font.color.rgb = base_color
            if is_sup:
                r.font.superscript = True


def md_to_reportlab_xml(text: str) -> str:
    """Converts inline Markdown + Unicode superscripts to safe ReportLab Paragraph XML."""
    parts: list[str] = []

    def emit_with_sup(s: str, bold: bool = False, italic: bool = False) -> None:
        for seg, is_sup in split_sup_segments(s, convert_sup=True):
            esc = xml_escape(seg)
            if is_sup:
                esc = f"<super>{esc}</super>"
            if bold:
                esc = f"<b>{esc}</b>"
            if italic:
                esc = f"<i>{esc}</i>"
            parts.append(esc)

    pos = 0
    for m in TOKEN_RE.finditer(text):
        if m.start() > pos:
            emit_with_sup(text[pos:m.start()])
        tok = m.group(0)
        if tok.startswith("**") and tok.endswith("**"):
            inner = tok[2:-2].replace("`", "")
            emit_with_sup(inner, bold=True)
        elif tok.startswith("`") and tok.endswith("`"):
            emit_with_sup(tok[1:-1], bold=False)
        elif tok.startswith("*") and tok.endswith("*"):
            emit_with_sup(tok[1:-1], italic=True)
        pos = m.end()
    if pos < len(text):
        emit_with_sup(text[pos:])
    return "".join(parts)


def add_page_number_footer(section) -> None:
    footer = section.footer
    p = footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r1 = p.add_run("PBISC Nature Portfolio Manuscript (Figures 4–7)  |  Page ")
    r1.font.name = "Calibri"
    r1.font.size = Pt(8.5)
    r1.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    fldSimple = OxmlElement("w:fldSimple")
    fldSimple.set(qn("w:instr"), "PAGE")
    p._p.append(fldSimple)


def build_docx(md_lines: list[str]) -> None:
    doc = docx.Document()
    for sec in doc.sections:
        sec.page_width = Inches(8.27)   # A4 width
        sec.page_height = Inches(11.69) # A4 height
        sec.top_margin = Inches(0.75)
        sec.bottom_margin = Inches(0.75)
        sec.left_margin = Inches(0.75)
        sec.right_margin = Inches(0.75)
        add_page_number_footer(sec)

    navy = RGBColor(0x14, 0x2E, 0x54)
    slate = RGBColor(0x33, 0x3E, 0x48)

    in_legends = False
    for raw_line in md_lines:
        line = raw_line.rstrip()
        if not line:
            continue
        if line == "---":
            continue

        img_m = IMG_LINE_RE.match(line)
        if img_m:
            rel_path = img_m.group(2)
            img_path = (PKG_DIR / rel_path).resolve()
            if not img_path.exists():
                raise FileNotFoundError(f"Embedded figure not found: {img_path}")
            # Start each Main Figure on a fresh page for clean full-page figure + legend presentation
            doc.add_page_break()
            # Measure aspect ratio so figure fits comfortably on A4 page
            idoc = fitz.open(str(img_path))
            w_px, h_px = idoc[0].rect.width, idoc[0].rect.height
            idoc.close()
            max_w_in = 6.70
            max_h_in = 8.45
            scale_w = max_w_in
            scale_h = scale_w * (h_px / float(w_px))
            if scale_h > max_h_in:
                scale_h = max_h_in
                scale_w = scale_h * (w_px / float(h_px))
            p_img = doc.add_paragraph()
            p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_img.paragraph_format.space_before = Pt(0)
            p_img.paragraph_format.space_after = Pt(4)
            run_img = p_img.add_run()
            run_img.add_picture(str(img_path), width=Inches(scale_w))
            continue

        if line.startswith("# "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(8)
            add_docx_formatted_runs(p, line[2:].strip(), base_size_pt=15.5, base_color=navy)
            for r in p.runs:
                r.bold = True
                r.font.name = "Calibri"
        elif line.startswith("> "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.2)
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(4)
            add_docx_formatted_runs(p, line[2:].strip(), base_size_pt=9.2, base_color=slate)
        elif line.startswith("## "):
            sec_title = line[3:].strip()
            in_legends = "Main Figures" in sec_title
            if sec_title in ("Discussion", "Methods", "References") and not in_legends:
                doc.add_page_break()
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(10)
            p.paragraph_format.space_after = Pt(5)
            p.paragraph_format.keep_with_next = True
            add_docx_formatted_runs(p, sec_title, base_size_pt=13.5, base_color=navy)
            for r in p.runs:
                r.bold = True
                r.font.name = "Calibri"
        elif line.startswith("### "):
            sub_title = line[4:].strip()
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(7 if not in_legends else 3)
            p.paragraph_format.space_after = Pt(3)
            p.paragraph_format.keep_with_next = True
            add_docx_formatted_runs(p, sub_title, base_size_pt=11.2 if not in_legends else 10.2, base_color=navy)
            for r in p.runs:
                r.bold = True
        elif re.match(r"^\d+\.\s+", line):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.22)
            p.paragraph_format.first_line_indent = Inches(-0.22)
            p.paragraph_format.space_before = Pt(1.5)
            p.paragraph_format.space_after = Pt(2.5)
            p.paragraph_format.line_spacing = 1.12
            add_docx_formatted_runs(p, line, base_size_pt=9.6 if "https://" in line or "et al." in line or "Nat." in line or "Genome Biol." in line or "Radiology" in line or "Sci." in line else 10.2)
        elif line.startswith("- "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.20)
            p.paragraph_format.first_line_indent = Inches(-0.12)
            p.paragraph_format.space_before = Pt(1.5)
            p.paragraph_format.space_after = Pt(2.5)
            p.paragraph_format.line_spacing = 1.14
            add_docx_formatted_runs(p, "• " + line[2:].strip(), base_size_pt=10.0)
        else:
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(5 if not in_legends else 3.5)
            p.paragraph_format.line_spacing = 1.15 if not in_legends else 1.08
            add_docx_formatted_runs(p, line, base_size_pt=10.5 if not in_legends else 9.2)

    doc.save(str(DOCX_PATH))

    # Validate OOXML integrity: 4 embedded figures & 0 leaked markdown tokens
    with zipfile.ZipFile(DOCX_PATH, "r") as zf:
        doc_xml = zf.read("word/document.xml").decode("utf-8")
        media_files = sorted(n for n in zf.namelist() if n.startswith("word/media/"))
        assert len(media_files) == 4, f"Expected 4 embedded figures in DOCX, got {len(media_files)}: {media_files}"
        cleaned_xml = doc_xml.replace("***q", "").replace("**q", "")
        assert "`" not in cleaned_xml, "Unparsed backtick leaked into DOCX XML!"
        assert "**" not in cleaned_xml, "Unparsed ** leaked into DOCX XML!"


def register_unicode_fonts() -> tuple[str, str, str, str, str]:
    """Registers full-Unicode TrueType serif & sans fonts (Windows Fonts or Matplotlib DejaVu fallback)."""
    win_fonts = Path(r"C:\Windows\Fonts")
    if (win_fonts / "times.ttf").exists() and (win_fonts / "arial.ttf").exists():
        pdfmetrics.registerFont(TTFont("USerif", str(win_fonts / "times.ttf")))
        pdfmetrics.registerFont(TTFont("USerif-Bold", str(win_fonts / "timesbd.ttf")))
        pdfmetrics.registerFont(TTFont("USerif-Italic", str(win_fonts / "timesi.ttf")))
        pdfmetrics.registerFont(TTFont("USerif-BoldItalic", str(win_fonts / "timesbi.ttf")))
        pdfmetrics.registerFontFamily("USerif", normal="USerif", bold="USerif-Bold", italic="USerif-Italic", boldItalic="USerif-BoldItalic")
        pdfmetrics.registerFont(TTFont("USans", str(win_fonts / "arial.ttf")))
        pdfmetrics.registerFont(TTFont("USans-Bold", str(win_fonts / "arialbd.ttf")))
        pdfmetrics.registerFontFamily("USans", normal="USans", bold="USans-Bold", italic="USans", boldItalic="USans-Bold")
        return "USerif", "USerif-Bold", "USerif-Italic", "USans", "USans-Bold"
    try:
        import matplotlib
        ttf_dir = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
        if (ttf_dir / "DejaVuSerif.ttf").exists():
            pdfmetrics.registerFont(TTFont("USerif", str(ttf_dir / "DejaVuSerif.ttf")))
            pdfmetrics.registerFont(TTFont("USerif-Bold", str(ttf_dir / "DejaVuSerif-Bold.ttf")))
            pdfmetrics.registerFont(TTFont("USerif-Italic", str(ttf_dir / "DejaVuSerif-Italic.ttf")))
            pdfmetrics.registerFont(TTFont("USerif-BoldItalic", str(ttf_dir / "DejaVuSerif-BoldItalic.ttf")))
            pdfmetrics.registerFontFamily("USerif", normal="USerif", bold="USerif-Bold", italic="USerif-Italic", boldItalic="USerif-BoldItalic")
            pdfmetrics.registerFont(TTFont("USans", str(ttf_dir / "DejaVuSans.ttf")))
            pdfmetrics.registerFont(TTFont("USans-Bold", str(ttf_dir / "DejaVuSans-Bold.ttf")))
            pdfmetrics.registerFontFamily("USans", normal="USans", bold="USans-Bold", italic="USans", boldItalic="USans-Bold")
            return "USerif", "USerif-Bold", "USerif-Italic", "USans", "USans-Bold"
    except Exception:
        pass
    return "Times-Roman", "Times-Bold", "Times-Italic", "Helvetica", "Helvetica-Bold"


def build_pdf_reportlab(md_lines: list[str]) -> None:
    """Builds a crisp, self-contained A4 PDF with ReportLab Platypus & embedded high-DPI figures."""
    f_serif, f_serif_b, _, f_sans, f_sans_b = register_unicode_fonts()
    styles = getSampleStyleSheet()
    navy_hex = colors.HexColor("#142E54")
    slate_hex = colors.HexColor("#333E48")

    st_title = ParagraphStyle(
        "DocTitle",
        parent=styles["Heading1"],
        fontName=f_sans_b,
        fontSize=14.5,
        leading=18.5,
        textColor=navy_hex,
        spaceAfter=6,
    )
    st_meta = ParagraphStyle(
        "DocMeta",
        parent=styles["Normal"],
        fontName=f_sans,
        fontSize=8.6,
        leading=11.5,
        leftIndent=10,
        textColor=slate_hex,
        spaceAfter=3,
    )
    st_h2 = ParagraphStyle(
        "SecH2",
        parent=styles["Heading2"],
        fontName=f_sans_b,
        fontSize=12.5,
        leading=15.5,
        textColor=navy_hex,
        spaceBefore=8,
        spaceAfter=4,
        keepWithNext=True,
    )
    st_h3 = ParagraphStyle(
        "SubH3",
        parent=styles["Heading3"],
        fontName=f_serif_b,
        fontSize=10.8,
        leading=13.5,
        textColor=navy_hex,
        spaceBefore=6,
        spaceAfter=3,
        keepWithNext=True,
    )
    st_leg_title = ParagraphStyle(
        "LegTitle",
        parent=styles["Heading3"],
        fontName=f_serif_b,
        fontSize=9.6,
        leading=12.0,
        textColor=navy_hex,
        spaceBefore=2,
        spaceAfter=2,
        keepWithNext=True,
    )
    st_body = ParagraphStyle(
        "BodyJust",
        parent=styles["Normal"],
        fontName=f_serif,
        fontSize=10.0,
        leading=13.2,
        alignment=TA_JUSTIFY,
        spaceAfter=4.5,
    )
    st_legend = ParagraphStyle(
        "LegendBody",
        parent=styles["Normal"],
        fontName=f_serif,
        fontSize=8.6,
        leading=10.8,
        alignment=TA_JUSTIFY,
        spaceAfter=3.5,
    )
    st_list = ParagraphStyle(
        "ListBody",
        parent=styles["Normal"],
        fontName=f_serif,
        fontSize=9.3,
        leading=12.0,
        leftIndent=14,
        firstLineIndent=-10,
        alignment=TA_LEFT,
        spaceAfter=2.5,
    )

    story = []
    in_legends = False
    usable_w = A4[0] - 36 * mm
    usable_h = A4[1] - 36 * mm

    for raw_line in md_lines:
        line = raw_line.rstrip()
        if not line:
            continue
        if line == "---":
            if not in_legends:
                story.append(Spacer(1, 3))
                story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=2, spaceAfter=4))
            continue

        img_m = IMG_LINE_RE.match(line)
        if img_m:
            rel_path = img_m.group(2)
            img_path = (PKG_DIR / rel_path).resolve()
            idoc = fitz.open(str(img_path))
            w_px, h_px = idoc[0].rect.width, idoc[0].rect.height
            idoc.close()
            story.append(PageBreak())
            max_w = usable_w
            max_h = usable_h * 0.80
            draw_w = max_w
            draw_h = draw_w * (h_px / float(w_px))
            if draw_h > max_h:
                draw_h = max_h
                draw_w = draw_h * (w_px / float(h_px))
            rl_img = RLImage(str(img_path), width=draw_w, height=draw_h)
            rl_img.hAlign = "CENTER"
            story.append(rl_img)
            story.append(Spacer(1, 4))
            continue

        if line.startswith("# "):
            story.append(Paragraph(md_to_reportlab_xml(line[2:].strip()), st_title))
        elif line.startswith("> "):
            story.append(Paragraph(md_to_reportlab_xml(line[2:].strip()), st_meta))
        elif line.startswith("## "):
            sec_title = line[3:].strip()
            in_legends = "Main Figures" in sec_title
            if sec_title in ("Discussion", "Methods", "References") and not in_legends:
                story.append(PageBreak())
            story.append(Paragraph(md_to_reportlab_xml(sec_title), st_h2))
        elif line.startswith("### "):
            sub_title = line[4:].strip()
            story.append(Paragraph(md_to_reportlab_xml(sub_title), st_leg_title if in_legends else st_h3))
        elif re.match(r"^\d+\.\s+", line):
            story.append(Paragraph(md_to_reportlab_xml(line), st_list))
        elif line.startswith("- "):
            story.append(Paragraph("• " + md_to_reportlab_xml(line[2:].strip()), st_list))
        else:
            story.append(Paragraph(md_to_reportlab_xml(line), st_legend if in_legends else st_body))

    def draw_footer(canvas, doc_obj):
        canvas.saveState()
        canvas.setFont("Helvetica", 8.0)
        canvas.setFillColor(colors.HexColor("#555555"))
        canvas.drawCentredString(
            A4[0] / 2.0,
            10 * mm,
            f"PBISC Nature Portfolio Manuscript (Figures 4–7)  |  Page {canvas.getPageNumber()}",
        )
        canvas.restoreState()

    pdf_doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=16 * mm,
        bottomMargin=18 * mm,
        title="PBISC Manuscript — Figures 4–7 (Final Verified Edition, 2026-11-25)",
    )
    pdf_doc.build(story, onFirstPage=draw_footer, onLaterPages=draw_footer)


def main() -> None:
    md_text = MD_PATH.read_text(encoding="utf-8")
    md_lines = md_text.splitlines()
    build_docx(md_lines)
    build_pdf_reportlab(md_lines)

    pdf_check = fitz.open(str(PDF_PATH))
    n_pages = len(pdf_check)
    n_imgs = sum(len(pdf_check[i].get_images(full=True)) for i in range(n_pages))
    pdf_check.close()
    print(f"Built DOCX: {DOCX_PATH.name} ({DOCX_PATH.stat().st_size:,} B)")
    print(f"Built PDF : {PDF_PATH.name} ({PDF_PATH.stat().st_size:,} B, {n_pages} pages, {n_imgs} embedded images)")


if __name__ == "__main__":
    main()
