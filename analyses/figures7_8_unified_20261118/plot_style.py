"""Shared final-size publication style; no statistical calculations."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
FIGURES = ROOT / 'figures'
COHORT_ORDER = ['Real', 'New', 'Combined']
COLORS = {'Real': '#5E6462', 'New': '#0F4D92', 'Combined': '#1F6F5F'}
MARKERS = {'Real': 'o', 'New': '^', 'Combined': 's'}
INK = '#171717'
MUTED = '#5E6462'
MM = 1 / 25.4

def setup_style():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8.5, 'axes.titlesize': 9, 'axes.labelsize': 8.5,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .7, 'xtick.major.width': .6, 'ytick.major.width': .6,
        'text.color': INK, 'axes.labelcolor': INK,
        'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'savefig.facecolor': 'white', 'figure.facecolor': 'white',
    })

def label_panel(ax, letter, title='', x=-.18, y=1.10):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
    if title:
        ax.set_title(title, loc='left', pad=9, fontweight='semibold')

def save_figure(fig, stem, *, metadata=None):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            texts.append({'text': t.get_text(), 'fontsize_pt': float(t.get_fontsize())})
    if any(t['fontsize_pt'] < 7.99 for t in texts):
        raise ValueError('Final-size text below 8 pt')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(FIGURES / f'{stem}.{ext}', dpi=300)
    record = {'width_mm': fig.get_figwidth()*25.4, 'height_mm': fig.get_figheight()*25.4,
              'texts': texts, 'metadata': metadata or {}}
    (FIGURES / f'{stem}.artists.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    plt.close(fig)
    return record

setup_style()
