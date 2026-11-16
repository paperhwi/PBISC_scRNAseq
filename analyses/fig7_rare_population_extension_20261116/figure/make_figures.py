"""Build the ESKD COVID clinical-association figure set (figures4papers style).
Each panel is drawn once by a function in panels.py and exported (1) as its own fixed-size file in panels_out/
and (2) placed on a full figure page in pages/. Canvases are 2x print size (page 14 in = 178 mm column width),
font 15 -> ~7.5 pt when placed; PNG 300 dpi, PDF and SVG vector with editable text."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from style import apply_publication_style
import panels as p
apply_publication_style()
OUT=Path('panels_out'); PAGES=Path('pages'); OUT.mkdir(exist_ok=True); PAGES.mkdir(exist_ok=True)

def simple(fn, l, r, b, t):
    def draw(fig, place, w, h): fn(place(l, b, w-l-r, h-b-t))
    return draw
def dotgrid(fig, place, w, h):
    ax=place(2.6, 2.45, w-2.6-0.25, h-2.45-1.75); cax=place(2.6, 1.3, (w-2.85)*0.62, 0.22); sax=place(2.6, 0.0, (w-2.85), 0.55)
    p.assoc_dotgrid(ax, cax, sax)
def dimdots(fig, place, w, h):
    ax=place(1.35, 1.05, w-1.35-1.65, h-1.05-0.35); cax=place(w-1.25, 1.05, 0.22, h-1.05-1.15); p.dim_marker_dots(ax, cax)
def markers(fig, place, w, h):
    ax=place(3.3, 1.75, w-3.3-1.9, h-1.75-1.05); cax=place(w-1.55, 1.95, 0.25, h-1.75-1.25); p.marker_dots(ax, cax)

XY=(1.35,0.25,0.95,0.35); YL=(2.6,0.25,0.95,0.35)
FIGS={
 'ClinFig1':dict(title='Cell-type composition vs COVID-19 severity: Real, PBISC and combined',size=(14,18.7),panels=[
   ('a','schematic',(0,0,14,3.1),simple(p.schematic,0.05,0.05,0.05,0.05)),
   ('b','assoc_dotgrid',(0,3.35,6.3,10.05),dotgrid),
   ('c','strip_NKCD56bright',(6.5,3.35,7.5,3.35),simple(lambda a:p.strip_severity(a,'rare_NK CD56 bright_pct',show_legend=True),1.35,0.2,0.65,0.3)),
   ('d','strip_DC',(6.5,6.7,7.5,3.35),simple(lambda a:p.strip_severity(a,'DC_within_Other_pct'),1.35,0.2,0.65,0.3)),
   ('e','strip_HSPC',(6.5,10.05,7.5,3.35),simple(lambda a:p.strip_severity(a,'rare_HSPC_pct'),1.35,0.2,0.65,0.3)),
   ('f','power_dumbbell',(0,13.65,8.4,5.05),simple(p.power_dumbbell,2.6,0.25,0.95,0.3)),
   ('g','roc_HSPC_death',(8.6,13.65,5.4,5.05),simple(p.roc_hspc,1.3,0.25,0.95,0.6))]),
 'ClinFig2':dict(title='Real and PBISC carry different composition signals',size=(14,14.2),panels=[
   ('a','reliability',(0,0,7.0,9.2),simple(p.reliability_dots,2.6,0.25,0.95,0.3)),
   ('b','severity_effect_concordance',(7.2,0,6.8,9.2),simple(p.effect_concordance,1.45,0.3,2.3,1.3)),
   ('c','HSPC_vs_NLR',(0,9.5,7.0,4.7),simple(lambda a:p.nlr_scatter(a,'rare_HSPC_pct','HSPC (% of cells)'),1.35,0.3,0.95,0.3)),
   ('d','T_vs_NLR',(7.2,9.5,6.8,4.7),simple(lambda a:p.nlr_scatter(a,'T_pct','T cell (% of cells)',(30,112)),1.35,0.3,0.95,0.3))]),
 'ClinFig3':dict(title='Severity DEG inside rare populations: Real, PBISC and combined',size=(14,19.4),panels=[
   ('a','DEG_counts',(0,0,5.7,8.2),simple(p.deg_count_dots,2.95,0.15,0.85,1.15)),
   ('b','CD14monoIFN_gene_trend',(5.9,0,8.1,8.2),simple(p.cd14ifn_gene_scatter,1.45,0.3,0.95,0.3)),
   ('c','permutation_control',(0,8.45,7.0,6.2),simple(p.perm_null,2.65,0.2,2.15,0.3)),
   ('d','broadcast',(7.2,8.45,6.8,6.2),simple(p.broadcast_dots,2.95,0.2,2.15,0.3)),
   ('e','S100A8',(0,14.9,4.67,4.5),simple(lambda a:p.gene_strip(a,'S100A8'),1.2,0.15,0.75,0.55)),
   ('f','CD163',(4.67,14.9,4.67,4.5),simple(lambda a:p.gene_strip(a,'CD163'),1.2,0.15,0.75,0.55)),
   ('g','HLA-DRA',(9.33,14.9,4.67,4.5),simple(lambda a:p.gene_strip(a,'HLA-DRA'),1.2,0.15,0.75,0.55))]),
 'ClinFig4':dict(title='Rare-population subclustering: reproducibility and the NK CD56bright cluster 17',size=(14,14.0),panels=[
   ('a','ARI_split_HVG',(0,0,6.0,8.8),simple(p.ari_dots,2.6,0.25,0.95,0.3)),
   ('b','NK_umap_source',(6.2,0,3.9,4.3),simple(lambda a:p.nk_umap(a,'source'),0.45,0.1,0.45,0.15)),
   ('c','NK_umap_cluster17',(10.1,0,3.9,4.3),simple(lambda a:p.nk_umap(a,'c17'),0.45,0.1,0.45,0.15)),
   ('d','cluster17_vs_severity',(6.2,4.4,7.8,4.4),simple(p.c17_strip,1.45,0.2,0.7,0.2)),
   ('e','cluster17_markers',(0,9.1,14,4.9),markers)]),
 'ClinSFig2':dict(title='Supplementary: NK CD56bright cluster 17 is not explained by CD56dim contamination (RNA; ADT not deposited)',size=(14,11.2),panels=[
   ('a','requested_markers',(0,0,5.6,5.0),dimdots),
   ('b','lineage_scores_Real',(5.8,0,4.1,5.0),simple(lambda a:p.score_scatter(a,'Real',legend=True),1.15,0.15,0.9,0.45)),
   ('c','lineage_scores_PBISC',(9.9,0,4.1,5.0),simple(lambda a:p.score_scatter(a,'PBISC'),1.15,0.15,0.9,0.45)),
   ('d','lineage_axis_ecdf',(0,5.3,4.5,5.9),simple(p.lineage_ecdf,1.25,0.15,0.95,1.1)),
   ('e','classifier_pdim_ecdf',(4.5,5.3,4.5,5.9),simple(p.pdim_ecdf,1.25,0.15,0.95,1.1)),
   ('f','severity_after_removing_dimlike',(9.0,5.3,5.0,5.9),simple(p.c17_robust,2.3,0.1,1.75,1.1))]),
 'ClinSFig1':dict(title='Supplementary: composition by peak severity and rare-population scale',size=(14,7.2),panels=[
   ('a','composition_by_severity',(0,0,7.8,7.2),simple(p.comp_stack,1.35,2.1,1.4,0.55)),
   ('b','rare_population_scale',(8.0,0,6.0,7.2),simple(p.rare_scale_bars,2.6,0.3,0.95,0.3))]),
}

def save(fig, stem):
    for ext in ('png','pdf','svg'):
        f=Path(f'{stem}.{ext}'); fig.savefig(f, dpi=300)
        if ext=='svg': f.write_text(f.read_text().replace('Liberation Sans','Arial'))

def build(key):
    spec=FIGS[key]; W,H=spec['size']; manifest=[]
    page=plt.figure(figsize=(W,H))
    for letter,name,(x,y,w,h),draw in spec['panels']:
        # own canvas
        f=plt.figure(figsize=(w,h)); draw(f, lambda l,b,ww,hh: f.add_axes([l/w,b/h,ww/w,hh/h]), w, h)
        stem=OUT/f'{key}_{letter}_{name}'; save(f, stem); plt.close(f)
        # on the page
        X0,Y0=x,H-y-h
        draw(page, lambda l,b,ww,hh,X0=X0,Y0=Y0: page.add_axes([(X0+l)/W,(Y0+b)/H,ww/W,hh/H]), w, h)
        page.text((x+0.05)/W, 1-(y+0.02)/H, letter, fontsize=24, fontweight='bold', va='top', ha='left')
        manifest.append(dict(figure=key,panel=letter,file=stem.name,box_in_2x=[x,y,w,h],print_mm=[round(w/2*25.4,1),round(h/2*25.4,1)]))
    save(page, PAGES/key); plt.close(page)
    return manifest

if __name__=='__main__':
    keys=sys.argv[1:] or list(FIGS)
    allm=[]
    for k in keys: allm+=build(k); print('built',k)
    old=json.loads(Path('panels_out/manifest.json').read_text()) if Path('panels_out/manifest.json').exists() else []
    old=[m for m in old if m['figure'] not in keys]
    Path('panels_out/manifest.json').write_text(json.dumps(old+allm,indent=1))
