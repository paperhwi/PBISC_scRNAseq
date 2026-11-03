from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import json, re
from pypdf import PdfReader

root = Path('F:/PBISC')
out = root / 'project_review_20261028' / 'extracted'
out.mkdir(parents=True, exist_ok=True)
sources = list(root.glob('*.pptx')) + list(root.glob('*.docx')) + list((root/'deliverables').glob('*.pptx')) + list((root/'1. PBISC Manuscript').glob('*'))
manifest = []
for p in sources:
    if p.suffix.lower() not in ('.pdf', '.pptx', '.docx'):
        continue
    parts = []
    if p.suffix.lower() == '.pdf':
        for i, page in enumerate(PdfReader(p).pages, 1):
            parts.append(f'\n## Page {i}\n' + (page.extract_text() or ''))
    else:
        with ZipFile(p) as z:
            if p.suffix.lower() == '.docx':
                names = ['word/document.xml']
            else:
                names = sorted([n for n in z.namelist() if re.match(r'ppt/(slides/slide|notesSlides/notesSlide)\d+\.xml$', n)], key=lambda n: ('notes' in n, int(re.search(r'(\d+)\.xml$',n).group(1))))
            for n in names:
                tree = ET.fromstring(z.read(n))
                paragraphs = []
                for element in tree.iter():
                    if element.tag.rsplit('}',1)[-1] == 'p':
                        s = ''.join(t.text or '' for t in element.iter() if t.tag.rsplit('}',1)[-1] == 't')
                        if s: paragraphs.append(s)
                parts.append('\n## ' + n + '\n' + '\n'.join(paragraphs))
    target = out / (p.stem + '.txt')
    target.write_text('\n'.join(parts), encoding='utf-8')
    manifest.append({'source':str(p), 'extracted':str(target), 'sections':len(parts), 'characters':sum(map(len,parts))})
(out.parent/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(f'Extracted {len(manifest)} source documents; manifest saved.')
