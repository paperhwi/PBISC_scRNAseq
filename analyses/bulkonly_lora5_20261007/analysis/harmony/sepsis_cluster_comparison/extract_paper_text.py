import sys
from pathlib import Path

from pypdf import PdfReader

sys.stdout.reconfigure(encoding="utf-8")
pdf = Path("F:/PBISC/2. Reference article/Single-cell multi-omic landscape reveals anatomical-specific immune features in adult and pediatric sepsis.pdf")
out = Path(__file__).parent / "sepsis_paper_text.txt"
reader = PdfReader(pdf)
out.write_text("\n\n".join(f"[PAGE {i + 1}]\n{page.extract_text() or ''}" for i, page in enumerate(reader.pages)), encoding="utf-8")
print(f"pages={len(reader.pages)} text={out}")
