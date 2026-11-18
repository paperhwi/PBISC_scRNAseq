"""Rebuild all seven figures from the bundled frozen tables."""
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
FIGURE_NAMES=['Figure_7','Figure_8']+[f'Supplementary_Figure_S{i}' for i in range(1,6)]

def main():
    for script in ['verify_release.py','build_figure7.py','build_figure8.py',
                   'build_supplementary.py','write_documentation.py']:
        subprocess.run([sys.executable,str(ROOT/script)],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/'verify_release.py'),'--figures'],cwd=ROOT,check=True)
    from pypdf import PdfWriter
    merged=PdfWriter()
    for name in FIGURE_NAMES:merged.append(str(ROOT/'figures'/f'{name}.pdf'))
    with (ROOT/'figures/Figures_7_8_and_S1_S5.pdf').open('wb') as f:merged.write(f)
    merged.close()
    versions={'python':sys.version.split()[0]}
    for name in ['numpy','pandas','scipy','matplotlib','pillow','pypdf']:
        versions[name]=importlib.metadata.version(name)
    (ROOT/'provenance/RENDER_SOFTWARE_VERSIONS.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
    write_manifest()
    print('Seven figures and combined PDF rebuilt; numerical checks passed.')

def write_manifest():
    supplementary=[ROOT/'build_supplementary.py']
    supplementary+=list((ROOT/'figures').glob('Supplementary_Figure_*'))
    supplementary+=list((ROOT/'source_data').glob('S*.csv'))
    supplementary+=list((ROOT/'supplementary_specs').glob('*.json'))
    supplementary+=list((ROOT/'supplementary_specs').glob('*.md'))
    with (ROOT/'supplementary_specs/supplementary_file_manifest.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);writer.writeheader()
        for p in sorted(set(supplementary)):
            writer.writerow(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,
                                 sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    rows=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or any(part in ['__pycache__','qa','qa_preview'] for part in p.relative_to(ROOT).parts):continue
        if p.name=='OUTPUT_MANIFEST.csv':continue
        rows.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,
                         sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    with (ROOT/'provenance/OUTPUT_MANIFEST.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256'])
        writer.writeheader();writer.writerows(rows)

if __name__=='__main__':main()
