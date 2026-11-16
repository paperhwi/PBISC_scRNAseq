"""Ver3 (2026-11-15): Figure 7 + Supplementary Figures 6-13 (ESKD COVID clinical association), figures4papers style.
Every panel is drawn once (panels.py) and exported as its own file (panels_fig7/, PNG 300 dpi at 2x print size + PDF + SVG),
placed on an assembled page (pages_fig7/) and listed in panels_fig7/manifest.json (slide inches, 178 x 258 mm page) for build_pptx_fig7.js."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from style import apply_publication_style
import panels as p
import panels_v2 as q
import panels_v3 as v
from make_figures import simple, save, markers, dimdots
apply_publication_style()
OUT=Path('panels_fig7_v3'); PAGES=Path('pages_fig7_v3'); OUT.mkdir(exist_ok=True); PAGES.mkdir(exist_ok=True)
PAGE_W, PAGE_H, MX, TOP = 7.0087, 10.158, 0.15, 0.42          # slide (1x, inches)
S=(PAGE_W-2*MX)*2/14.0                                          # design canvas (2x) is 14 in wide -> placed scale
def dotgrid7(fig, place, w, h):
    ax=place(2.6, 2.0, w-2.6-0.2, h-2.0-1.6); cax=place(2.6, 1.15, (w-2.8)*0.8, 0.2); sax=place(2.0, 0.0, w-2.05, 0.5); p.assoc_dotgrid(ax, cax, sax)
STRIP=(1.4,0.2,0.65,0.3)
FIGS={
 'Fig7':dict(title='Figure 7. What bulk-inferred single cells add to patient-level clinical association in ESKD COVID-19',panels=[
   ('a','design',(0,0,14,2.4),simple(v.schematic_v3,0.05,0.05,0.05,0.05)),
   ('b','severity_dotgrid',(0,2.55,5.6,8.6),v.dotgrid_v3),
   ('c','NKCD56bright_vs_severity',(5.8,2.55,8.2,3.1),simple(lambda a:p.strip_severity(a,'rare_NK CD56 bright_pct',show_legend=True),*STRIP)),
   ('d','source_comparison',(5.8,5.8,8.2,5.35),simple(v.source_compare,2.0,0.15,1.55,0.15)),
   ('e','cohort_extension_power',(0,11.35,7.6,4.2),v.ext_power),
   ('f','prediction_powered_inference',(7.8,11.35,6.2,4.2),simple(v.ppi_panel,1.75,0.15,2.3,0.15)),
   ('g','GZMB_controls',(0,15.75,4.6,4.2),simple(lambda a:v.gzmb_design_strip(a,'GZMB'),1.25,0.1,0.85,0.1)),
   ('h','PRF1_controls',(4.6,15.75,4.6,4.2),simple(lambda a:v.gzmb_design_strip(a,'PRF1'),1.25,0.1,0.85,0.1)),
   ('i','signal_from_cells_to_bulk',(9.2,15.75,4.8,4.2),simple(v.signal_ladder,2.55,0.1,0.85,0.3))]),
 'SFig6':dict(title='Supplementary Figure 6. Cell-type composition by peak severity and HSPC fraction as a classifier of death',panels=[
   ('a','DC_vs_severity',(0,0,7.0,3.6),simple(lambda a:p.strip_severity(a,'DC_within_Other_pct',show_legend=True),*STRIP)),
   ('b','HSPC_vs_severity',(7.0,0,7.0,3.6),simple(lambda a:p.strip_severity(a,'rare_HSPC_pct'),*STRIP)),
   ('c','composition_by_severity',(0,3.85,8.0,6.8),simple(p.comp_stack,1.35,2.1,1.4,0.55)),
   ('d','rare_population_scale',(8.2,3.85,5.8,6.8),simple(p.rare_scale_bars,2.6,0.3,0.95,0.3)),
   ('e','HSPC_death_ROC',(0,10.9,5.4,5.0),simple(p.roc_hspc,1.3,0.25,0.95,0.6))]),
 'SFig7':dict(title='Supplementary Figure 7. Patient-level agreement between Real and PBISC compositions',panels=[
   ('a','reliability_vs_cross_source',(0,0,7.0,9.2),simple(p.reliability_dots,2.6,0.25,0.95,0.3)),
   ('b','severity_effect_concordance',(7.2,0,6.8,9.2),simple(p.effect_concordance,1.45,0.3,2.3,1.3)),
   ('c','HSPC_vs_NLR',(0,9.5,7.0,4.7),simple(lambda a:p.nlr_scatter(a,'rare_HSPC_pct','HSPC (% of cells)'),1.35,0.3,0.95,0.3)),
   ('d','T_vs_NLR',(7.2,9.5,6.8,4.7),simple(lambda a:p.nlr_scatter(a,'T_pct','T cell (% of cells)',(30,112)),1.35,0.3,0.95,0.3))]),
 'SFig8':dict(title='Supplementary Figure 8. Severity-associated expression inside rare populations',panels=[
   ('a','DEG_counts',(0,0,6.0,8.2),simple(p.deg_count_dots,2.95,0.15,0.85,1.15)),
   ('b','permutation_control',(6.2,0,7.8,8.2),simple(p.perm_null,2.65,0.2,2.15,0.3)),
   ('c','broadcast',(0,8.45,7.0,6.2),simple(p.broadcast_dots,2.95,0.2,2.15,0.3)),
   ('d','CD14monoIFN_gene_trend',(7.2,8.45,6.8,6.2),simple(p.cd14ifn_gene_scatter,1.45,0.3,0.95,0.3)),
   ('e','S100A8',(0,14.9,4.67,4.5),simple(lambda a:p.gene_strip(a,'S100A8'),1.2,0.15,0.75,0.55)),
   ('f','CD163',(4.67,14.9,4.67,4.5),simple(lambda a:p.gene_strip(a,'CD163'),1.2,0.15,0.75,0.55)),
   ('g','HLA-DRA',(9.33,14.9,4.67,4.5),simple(lambda a:p.gene_strip(a,'HLA-DRA'),1.2,0.15,0.75,0.55))]),
 'SFig9':dict(title='Supplementary Figure 9. NK-lineage subclustering of Real and PBISC cells and the cohort-extension test',panels=[
   ('a','NKlineage_umap_clusters',(0,0,3.4,3.9),simple(lambda a:q.nkv2_umap(a,'cluster'),0.45,0.1,0.45,0.1)),
   ('b','NKlineage_umap_source',(3.5,0,3.4,3.9),simple(lambda a:q.nkv2_umap(a,'source'),0.45,0.1,0.45,0.1)),
   ('c','NKlineage_umap_label',(7.0,0,3.4,3.9),simple(lambda a:q.nkv2_umap(a,'fine'),0.45,0.1,0.45,0.1)),
   ('d','NKlineage_cluster_composition',(10.5,0,3.5,3.9),simple(q.nkv2_composition,0.85,0.75,1.4,0.45)),
   ('e','NKlineage_cluster_markers',(0,4.1,14,4.5),q.nkv2_markers),
   ('f','NKlineage_cluster_rho_extension',(0,8.8,5.6,4.0),simple(v.nk_rho_ext,3.1,0.15,1.6,0.15)),
   ('g','NKlineage_extension_power',(5.7,8.8,4.5,4.0),lambda fig,place,w,h: v.ext_power(fig,place,w,h,feats=[('4','cl. 4 ISG-high'),('6','cl. 6 proliferating')],src='nk',compact=True)),
   ('h','NKlineage_ARI',(10.3,8.8,3.7,4.0),simple(v.nk_ari_compact,1.45,0.2,0.85,0.3)),
   ('i','cluster4_by_group',(0,13.0,4.67,4.2),simple(lambda a:q.nkv2_frac_strip(a,'4',legend=False),1.15,0.1,0.55,0.45)),
   ('j','cluster6_by_group',(4.67,13.0,4.67,4.2),simple(lambda a:q.nkv2_frac_strip(a,'6'),1.15,0.1,0.55,0.45)),
   ('k','cluster7_by_group',(9.33,13.0,4.67,4.2),simple(lambda a:q.nkv2_frac_strip(a,'7',legend=True),1.15,0.1,0.55,0.45)),
   ('l','NKlineage_resolution',(0,17.4,7.0,2.5),simple(q.nkv2_resolution,1.3,0.3,0.9,0.1))]),
 'SFig10':dict(title='Supplementary Figure 10. NK CD56bright subclustering and the cytotoxic cluster 17',panels=[
   ('a','NKbright_umap_cluster17',(0,0,4.6,4.3),simple(lambda a:p.nk_umap(a,'c17'),0.45,0.1,0.45,0.15)),
   ('b','NKbright_umap_source',(4.7,0,4.6,4.3),simple(lambda a:p.nk_umap(a,'source'),0.45,0.1,0.45,0.15)),
   ('c','NKbright_umap_severity',(9.4,0,4.6,4.3),simple(p.nk_umap_severity,0.45,0.1,0.45,0.15)),
   ('d','cluster17_vs_severity',(0,4.5,6.0,4.4),simple(p.c17_strip,1.4,0.15,0.65,0.2)),
   ('e','cluster_rho_by_source',(6.2,4.5,7.8,4.4),simple(p.cluster_rho_scatter,1.35,0.2,0.9,0.2)),
   ('f','cluster17_markers',(0,9.1,14,4.9),markers)]),
 'SFig11':dict(title='Supplementary Figure 11. NK CD56bright cluster 17 is not explained by CD56dim contamination',panels=[
   ('a','requested_markers',(0,0,5.6,5.0),dimdots),
   ('b','lineage_scores_Real',(5.8,0,4.1,5.0),simple(lambda a:p.score_scatter(a,'Real',legend=True),1.15,0.15,0.9,0.45)),
   ('c','lineage_scores_PBISC',(9.9,0,4.1,5.0),simple(lambda a:p.score_scatter(a,'PBISC'),1.15,0.15,0.9,0.45)),
   ('d','lineage_axis_ecdf',(0,5.3,4.5,5.9),simple(p.lineage_ecdf,1.25,0.15,0.95,1.1)),
   ('e','classifier_pdim_ecdf',(4.5,5.3,4.5,5.9),simple(p.pdim_ecdf,1.25,0.15,0.95,1.1)),
   ('f','severity_after_removing_dimlike',(9.0,5.3,5.0,5.9),simple(p.c17_robust,2.3,0.1,1.75,1.1))]),
 'SFig12':dict(title='Supplementary Figure 12. NK GZMB and PRF1: control designs, NK specificity, within-patient change and the bulk input',panels=[
   ('a','GZMB_design_AUC',(0,0,7.0,5.2),simple(lambda a:v.design_auc(a,'GZMB'),2.6,0.1,1.75,0.35)),
   ('b','PRF1_design_AUC',(7.0,0,7.0,5.2),simple(lambda a:v.design_auc(a,'PRF1'),2.6,0.1,1.75,0.35)),
   ('c','GZMB_pre_vs_infection',(0,5.4,3.5,4.4),simple(lambda a:v.paired_lines(a,'GZMB'),1.2,0.1,1.0,0.1)),
   ('d','PRF1_pre_vs_infection',(3.5,5.4,3.5,4.4),simple(lambda a:v.paired_lines(a,'PRF1'),1.2,0.1,1.0,0.1)),
   ('e','GZMB_bulk_vs_scRNA_same_day',(7.0,5.4,3.5,4.4),simple(lambda a:v.bulk_vs_real(a,'GZMB'),1.05,0.1,0.85,0.1)),
   ('f','PRF1_bulk_vs_scRNA_same_day',(10.5,5.4,3.5,4.4),simple(lambda a:v.bulk_vs_real(a,'PRF1'),1.05,0.1,0.85,0.1)),
   ('g','NKbright_GZMB_trend_COVIDpos',(0,10.0,7.0,3.6),simple(lambda a:p.nk_gene_strip(a,'GZMB'),1.45,0.2,0.65,0.2)),
   ('h','NKbright_PRF1_trend_COVIDpos',(7.0,10.0,7.0,3.6),simple(lambda a:p.nk_gene_strip(a,'PRF1'),1.45,0.2,0.65,0.2)),
   ('i','same_patient_GZMB',(0,13.8,3.6,4.6),simple(lambda a:q.nkbright_same_patient(a,'GZMB'),1.1,0.15,0.75,0.15)),
   ('j','same_patient_PRF1',(3.65,13.8,3.6,4.6),simple(lambda a:q.nkbright_same_patient(a,'PRF1'),1.1,0.15,0.75,0.15)),
   ('k','AUC_vs_controls_by_severity',(7.4,13.8,6.6,4.6),simple(v.auc_by_severity_indep,1.75,0.15,1.35,0.15))]),
 'SFig13':dict(title='Supplementary Figure 13. Robustness of the composition associations and of what PBISC adds',panels=[
   ('a','extension_power_all_features',(0,0,6.4,9.4),simple(v.ext_power_all,2.0,0.4,1.65,0.15)),
   ('b','robustness_rho',(6.6,0,7.4,4.6),simple(v.robustness_rho,1.7,2.55,0.85,0.15)),
   ('c','deconvolution_accuracy',(6.6,4.8,7.4,4.6),simple(v.deconv_accuracy,1.85,0.2,1.15,0.15)),
   ('d','covariate_adjustment',(0,9.6,8.0,6.4),simple(v.covariate_p,2.75,0.6,0.15,1.0))]),
}
def build(key):
    spec=FIGS[key]; H=max(y+h for _,_,(x,y,w,h),_ in spec['panels']); W=14.0
    page=plt.figure(figsize=(W,H)); items=[]
    for letter,name,(x,y,w,h),draw in spec['panels']:
        f=plt.figure(figsize=(w,h)); draw(f, lambda l,b,ww,hh: f.add_axes([l/w,b/h,ww/w,hh/h]), w, h)
        stem=OUT/f'{key}_{letter}_{name}'; save(f, stem); plt.close(f)
        X0,Y0=x,H-y-h; draw(page, lambda l,b,ww,hh,X0=X0,Y0=Y0: page.add_axes([(X0+l)/W,(Y0+b)/H,ww/W,hh/H]), w, h)
        page.text((x+0.05)/W, 1-(y+0.02)/H, letter, fontsize=24, fontweight='bold', va='top', ha='left')
        items.append(dict(file=stem.name+'.png',label=letter,x=round(MX+x*S/2,4),y=round(TOP+y*S/2,4),w=round(w*S/2,4),h=round(h*S/2,4),print_mm=[round(w*S/2*25.4,1),round(h*S/2*25.4,1)]))
    save(page, PAGES/key); plt.close(page)
    return dict(key=key,title=spec['title'],items=items,bottom=round(TOP+H*S/2,4))
if __name__=='__main__':
    keys=sys.argv[1:] or list(FIGS); man=OUT/'manifest.json'
    old={m['key']:m for m in json.loads(man.read_text())} if man.exists() else {}
    for k in keys: old[k]=build(k); print('built',k,'bottom',old[k]['bottom'])
    man.write_text(json.dumps([old[k] for k in FIGS if k in old],indent=1))
