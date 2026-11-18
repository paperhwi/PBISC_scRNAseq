"""Figure 8: fixed-candidate evidence from the unified patient-level analysis.

Only staged CSVs are read. This script does not refit a generator, recompute
multiple-testing families, or select additional candidates. All uncertainty
and correction values come from their canonical staged result tables.
"""
from pathlib import Path
import sys,json,hashlib,datetime
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
from plot_style import DATA,FIGURES,MM,COLORS,MARKERS,setup_style,label_panel,save_figure

SOURCE=ROOT/'source_data';SOURCE.mkdir(exist_ok=True)
QA=ROOT/'qa';QA.mkdir(exist_ok=True)
ORDER=['Real','New','Combined'];W,H=178,206
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
paths={
 'coverage':DATA/'transcriptome'/'all14_primary_coverage_and_composition_rarity.csv',
 'genes':DATA/'transcriptome'/'shared_unique_gene_universe.csv',
 'candidates':DATA/'transcriptome'/'fixed_primary_candidates_with_rarity_and_composition.csv',
 'CI':DATA/'transcriptome'/'fixed_primary_candidate_CIs_50000.csv',
 'permutation':DATA/'transcriptome_verification'/'KIAA0391_source_permutation_10M.csv',
 'Bonferroni':DATA/'transcriptome_verification'/'KIAA0391_source_permutation_10M_Bonferroni_sensitivity.csv',
 'patients':DATA/'plotting'/'RNA_candidate_patients.csv',
 'bridge':DATA/'plotting'/'KIAA0391_shared_bridge.csv'}
tab={k:pd.read_csv(p) for k,p in paths.items()}
coverage=tab['coverage'];cand=tab['candidates'];cand=cand[(cand.celltype=='HSPC')&(cand.gene=='KIAA0391')&(cand.endpoint=='severity')]
assert len(cand)==1;cand=cand.iloc[0]
ci=tab['CI'];ci=ci[(ci.celltype=='HSPC')&(ci.gene=='KIAA0391')&(ci.endpoint=='severity')].set_index('cohort').loc[ORDER]
perm=tab['permutation'];perm=perm[(perm.celltype=='HSPC')&(perm.gene=='KIAA0391')&(perm.endpoint=='severity')].set_index('cohort').loc[ORDER]
bonf=tab['Bonferroni'];bonf=bonf[(bonf.cohort=='Combined')&(bonf.family=='one_celltype_joint_endpoints')];assert len(bonf)==1;bonf=bonf.iloc[0]
patients=tab['patients'];patients=patients[(patients.celltype=='HSPC')&(patients.gene=='KIAA0391')].copy()
bridge=tab['bridge'].copy();nlabels=len(coverage);neligible=int(coverage.eligible.sum());nrare=int((coverage.eligible&coverage.rare_primary).sum());ngenes=len(tab['genes'])
ns=patients.groupby('arm').size().to_dict();ns['Combined']=len(patients)
assert (nlabels,neligible,nrare,ngenes)==(14,9,7,5762)
assert ns=={'New':29,'Real':33,'Combined':62} and len(bridge)==24
assert patients.patient_id.nunique()==62 and (patients.ncells>=3).all()
assert (bridge[['ncellsReal','ncellsPBISC']]>=3).all().all()
assert ci['B'].eq(50000).all() and perm['B'].eq(9999999).all() and int(bonf.m)==11524
for arm in ORDER:
 assert int(ci.loc[arm,'n'])==ns[arm]
 assert abs(ci.loc[arm,'rho']-perm.loc[arm,'rho'])<1e-12
 assert abs(ci.loc[arm,'rho']-cand[arm+'_effect'])<1e-12
 mask=np.ones(len(patients),bool) if arm=='Combined' else patients.arm.eq(arm)
 r=stats.spearmanr(patients.loc[mask,'score'],patients.loc[mask,'severity'])
 assert abs(r.statistic-ci.loc[arm,'rho'])<1e-12
 assert abs(r.pvalue-ci.loc[arm,'p_severity'])<1e-12
br=stats.spearmanr(bridge.Real_expression,bridge.PBISC_expression)
assert abs(br.statistic-cand.bridge_rho)<1e-12 and abs(br.pvalue-cand.bridge_P)<1e-12

# Evidence lock and complete panel-level source arrays are written before plotting.
pd.DataFrame([{'shared_labels':nlabels,'eligible_labels':neligible,'eligible_rare_labels':nrare,'shared_unique_genes':ngenes,'HSPC_Real':ns['Real'],'HSPC_New':ns['New'],'HSPC_Combined':ns['Combined'],'minimum_cells_per_patient':3,'eligible_minimum_patients_per_arm':10,'eligible_minimum_deaths_per_arm':3,'eligible_minimum_survivors_per_arm':3}]).to_csv(SOURCE/'Figure_8_A_search_space.csv',index=False)
b=ci.reset_index().merge(perm.reset_index()[['cohort','P_plus1','B','seed']],on='cohort',suffixes=('_bootstrap','_permutation'),validate='one_to_one')
b.to_csv(SOURCE/'Figure_8_B_KIAA0391_effects_and_CI.csv',index=False)
new=patients[patients.arm.eq('New')].sort_values('patient_id').copy()
# Deterministic x-only beeswarm reveals seven exactly-zero mild-patient scores.
# Distances use the final panel geometry; y values and association inputs stay fixed.
new['severity_display']=new.severity
xp=(W*.290)*72/25.4/4.;yp=(H*.205)*72/25.4/(.145+.006)
for severity,gg in new.groupby('severity'):
 placed=[]
 for idx,row in gg.sort_values(['mean_logCP10K','patient_id']).iterrows():
  for shift in [0,.135,-.135,.270,-.270,.405,-.405]:
   if all(((shift-ss)*xp)**2+((row.mean_logCP10K-yy)*yp)**2>=4.8**2 for ss,yy in placed):break
  else:raise ValueError('Beeswarm cannot place a new patient without overlap')
  new.loc[idx,'severity_display']=severity+shift;placed.append((shift,row.mean_logCP10K))
new.to_csv(SOURCE/'Figure_8_C_New29_patient_points.csv',index=False)
pd.DataFrame([{'cohort':'New','n':len(new),'rho':ci.loc['New','rho'],'P_analytic':ci.loc['New','p_severity'],'P_source_permutation_10M':perm.loc['New','P_plus1'],'B_source_permutation':int(perm.loc['New','B']),'peak_severity_mapping':'0=Mild;1=Moderate;2=Severe;3=Critical','display_jitter':'deterministic x-only beeswarm; no jitter in statistical tests'}]).to_csv(SOURCE/'Figure_8_C_New29_statistics.csv',index=False)
drows=[
 {'method':'Analytic BH q','scope':'HSPC, severity','m':ngenes,'value':cand.Combined_q_gene},
 {'method':'Analytic BH q','scope':'9 types, severity','m':neligible*ngenes,'value':cand.Combined_q_global},
 {'method':'Analytic BH q','scope':'9 types, 2 endpoints','m':neligible*ngenes*2,'value':cand.Combined_q_global_joint},
 {'method':'10M source permutation Bonferroni P','scope':'HSPC, 2 endpoints','m':int(bonf.m),'value':bonf.Bonferroni_plus1},
 {'method':'10M source permutation Bonferroni tail CP95 upper','scope':'HSPC, 2 endpoints','m':int(bonf.m),'value':bonf.Bonferroni_tail_CP95_high}]
pd.DataFrame(drows).to_csv(SOURCE/'Figure_8_D_correction_scopes.csv',index=False)
bridge.to_csv(SOURCE/'Figure_8_E_shared_patient_points.csv',index=False)
pd.DataFrame([{'celltype':'HSPC','gene':'KIAA0391','n_shared':len(bridge),'rho':cand.bridge_rho,'P':cand.bridge_P,'timing':'all-visit patient means; not same-day matching','primary_threshold':'at least3 HSPC cells per patient in each source','scope':'high reconstruction fidelity is not established'}]).to_csv(SOURCE/'Figure_8_E_shared_patient_statistics.csv',index=False)
panels=[
 {'panel':'A','form':'quantitative diagram','source':'Figure_8_A_search_space.csv','message':'Complete shared-label/shared-gene universe and actual HSPC coverage'},
 {'panel':'B','form':'dot and interval','source':'Figure_8_B_KIAA0391_effects_and_CI.csv','message':'Real/New/Combined KIAA0391 severity correlations; fixed-score50k percentile95% CI'},
 {'panel':'C','form':'discrete-x patient scatter','source':'Figure_8_C_New29_patient_points.csv','message':'All29 genuine new patients; raw meanlogCP10K; deterministic x-only beeswarm is visual only'},
 {'panel':'D','form':'grouped exact-value table','source':'Figure_8_D_correction_scopes.csv','message':'Analytic BH and10M full-HSPC Bonferroni are distinct correction scopes'},
 {'panel':'E','form':'paired-source patient scatter','source':'Figure_8_E_shared_patient_points.csv','message':'All24 shared patients; negative observed rho does not support strong fidelity'}]
pd.DataFrame(panels).to_csv(SOURCE/'Figure_8_panel_manifest.csv',index=False)
spec={
 'claim':'PBISC nominates KIAA0391 as a severity-associated transcriptomic candidate in inferred HSPCs; shared-patient fidelity remains unresolved.',
 'audience':'expert','visual_form':'small_multiples',
 'metric':{'name':'patient-level KIAA0391 severity association','unit':'Spearman rho; patient meanlog1p(CP10K); explicitly scoped correctedP/q','direction':'neutral','aggregation':'cell-weighted meanlog1p(CP10K) per patient over shared5762 genes; source-specific rankINT for pooled correlations','denominator':'Real33, genuineNew29, Combined62; shared24 bridge','precision_decimals':3,'scale':'linear'},
 'data_source':{'path_or_identifier':'; '.join(str(p.relative_to(ROOT)) for p in paths.values()),'immutable_fields':['patient membership','rho','P','q','CI limits','B','gene universe','correction family'],'checksum':sha(paths['candidates']),'notes':'All values are from staged canonical CSVs. Full9type permutation BH was not recomputed. The10M MonteCarlo tail upper bound is not a clinical effect CI.'},
 'focal_series':'Combined','required_baselines':['Real','New'],'series_order':ORDER,
 'category_order':[p['panel']+': '+p['message'] for p in panels],'ordering_rationale':'baseline_to_focal',
 'uncertainty':{'available':True,'type':'bootstrap_interval','level':.95,'sample_or_seed_count':50000,'note':'B50,000 fixed-score patient bootstrap; rho source-stratified. Source rankINT maps, generator, annotation and selection fixed.10M permutation B9,999,999 with+1P; Clopper-Pearson MC interval applies only to tail probability.'},
 'target':{'medium':'paper','width':W,'height':H,'unit':'mm','minimum_text_size_pt':8,'background':'white'},
 'annotations':[{'purpose':'threshold','text':'P/q threshold0.05; correction families are printed explicitly','target':'D'},{'purpose':'counterexample','text':'Shared-patient fidelity remains unresolved','target':'E'}],
 'selected_reference_ids':['OAI-VIS-005','OAI-VIS-008','OAI-VIS-013','OAI-VIS-014'],
 'output_formats':['pdf','svg','png'],'caption_requirements':['sample','estimator','uncertainty','exclusions','metric_direction'],
 'exceptions':[{'rule':'Prefer1-4 panels for a main figure','reason':'Five panels were explicitly approved to show discovery scope, effect uncertainty, patient distribution, correction family and negative fidelity evidence together.'}]}
(ROOT/'figure8_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')

setup_style();fig=plt.figure(figsize=(W*MM,H*MM))
fig.text(.04,.98,'KIAA0391 is a severity-associated candidate in inferred HSPCs',fontsize=10.5,fontweight='semibold',ha='left',va='top')
axs={
 'A':fig.add_axes([.075,.790,.875,.125]),
 'B':fig.add_axes([.180,.490,.330,.205]),
 'C':fig.add_axes([.660,.490,.290,.205]),
 'D':fig.add_axes([.075,.125,.445,.240]),
 'E':fig.add_axes([.660,.125,.290,.240])}
def head(ax,letter,title,column='left'):
 p=ax.get_position();anchor=.04 if column=='left' else .565;title_x=.075 if column=='left' else .615
 label_panel(ax,letter,x=(anchor-p.x0)/p.width,y=1.18)
 ax.text((title_x-p.x0)/p.width,1.18,title,transform=ax.transAxes,fontsize=9,fontweight='semibold',va='bottom',ha='left')
def fmtp(p):return f'{p:.4f}' if p>=.001 else f'{p:.2e}'

# A: no proportional areas; all counts remain literal and source-backed.
ax=axs['A'];ax.axis('off');head(ax,'A','Complete shared-gene search')
centers=[.08,.34,.60,.86]
nums=[str(nlabels),str(neligible),f'{ngenes:,}','HSPC']
descs=['shared cell labels',f'eligible labels\n{nrare} meet rarity criterion','shared genes',f'Real {ns["Real"]} + New {ns["New"]}\nCombined {ns["Combined"]} patients']
for x,n,t in zip(centers,nums,descs):
 ax.text(x,.68,n,ha='center',va='center',fontsize=14,fontweight='semibold',color=COLORS['Combined'] if n=='HSPC' else '#171717')
 ax.text(x,.30,t,ha='center',va='center',fontsize=8.5,linespacing=1.35)
for x in [.21,.73]:ax.annotate('',xy=(x+.035,.64),xytext=(x-.035,.64),arrowprops={'arrowstyle':'->','lw':.8,'color':'#5E6462'})
ax.text(.47,.65,'×',ha='center',va='center',fontsize=12,color='#5E6462')
ax.text(.5,-.20,'Primary RNA coverage: at least 3 cells per patient; all eligible labels and shared genes retained',ha='center',va='top',fontsize=8,color='#5E6462')

# B: actual effect estimates and actual bootstrap intervals, never pseudo-CIs.
ax=axs['B'];head(ax,'B','Severity association')
yy=[2,1,0]
for y,arm in zip(yy,ORDER):
 r=ci.loc[arm];ax.errorbar(r.rho,y,xerr=[[r.rho-r.rho_lo],[r.rho_hi-r.rho]],fmt=MARKERS[arm],markersize=5.2,color=COLORS[arm],ecolor=COLORS[arm],elinewidth=.9,capsize=2.2,markeredgewidth=.8)
 ax.text(r.rho,y+.22,f'{r.rho:.3f}',ha='center',va='bottom',fontsize=8,color=COLORS[arm])
ax.axvline(0,color='#8D8D88',lw=.7,ls='--',zorder=0)
ax.set_yticks(yy,[f'Real (n={ns["Real"]})\nAnalytic P = {fmtp(ci.loc["Real","p_severity"])}',f'New (n={ns["New"]})',f'Combined (n={ns["Combined"]})'])
ax.tick_params(axis='y',length=0,pad=8);ax.spines['left'].set_visible(False)
ax.set(xlim=(-.10,1.01),ylim=(-.48,2.60),xlabel='Spearman rho (95% bootstrap CI)')
ax.set_xticks([0,.5,1]);ax.tick_params(axis='x',pad=3)

# C: every genuine new patient shown at the untransformed expression-score scale.
ax=axs['C'];head(ax,'C',f'New patients (n={len(new)})',column='right')
ax.scatter(new.severity_display,new.mean_logCP10K,s=18,marker=MARKERS['New'],c=COLORS['New'],alpha=.95,linewidths=.4,edgecolors='white')
ax.set(xlim=(-.50,3.50),ylim=(-.006,.145),xlabel='Peak severity',ylabel='Mean log1p(CP10K)')
ax.set_xticks([0,1,2,3],['Mild','Moderate','Severe','Critical']);ax.set_yticks([0,.05,.10]);ax.text(.5,1.015,f'rho = {ci.loc["New","rho"]:.3f}\n10M source-permutation P = {perm.loc["New","P_plus1"]:.1e}',ha='center',va='bottom',transform=ax.transAxes,fontsize=8,linespacing=1.25)

# D: distinct analytic and permutation sections; no false common-statistic axis.
ax=axs['D'];head(ax,'D','Correction scope');ax.axis('off')
ax.text(0,1.01,'Analytic BH q · Combined',fontsize=8.5,fontweight='semibold',ha='left',va='top')
ax.plot([0,1],[.91,.91],color='#8D8D88',lw=.6,transform=ax.transAxes)
for y,row in zip([.84,.69,.54],drows[:3]):
 ax.text(0,y,f'{row["scope"]} ({row["m"]:,})',fontsize=8,ha='left',va='center')
 ax.text(1,y,f'{row["value"]:.4g}',fontsize=8.5,ha='right',va='center')
ax.text(0,.37,'10M permutation · HSPC, 2 endpoints',fontsize=8.5,fontweight='semibold',ha='left',va='center')
ax.plot([0,1],[.29,.29],color='#8D8D88',lw=.6,transform=ax.transAxes)
for y,label,row in zip([.20,.05],['Bonferroni P (11,524)','95% Monte Carlo upper bound'],drows[3:]):
 ax.text(0,y,label,fontsize=8,ha='left',va='center');ax.text(1,y,f'{row["value"]:.4g}',fontsize=8.5,ha='right',va='center',color=COLORS['Combined'])
ax.set(xlim=(0,1),ylim=(0,1.05))
ax.text(0,-.15,'Threshold 0.05; global permutation BH not computed',fontsize=8,ha='left',va='top',transform=ax.transAxes,color='#5E6462')

# E: the complete paired-source patient scatter is deliberately retained.
ax=axs['E'];head(ax,'E','Shared-patient fidelity',column='right')
ax.scatter(bridge.Real_expression,bridge.PBISC_expression,s=24,facecolors='none',edgecolors='#5E6462',linewidths=.8)
ax.set(xlim=(-.012,.30),ylim=(-.012,.30),xlabel='Real mean log1p(CP10K)',ylabel='PBISC mean log1p(CP10K)')
ax.set_xticks([0,.15,.30]);ax.set_yticks([0,.15,.30]);ax.set_aspect('equal',adjustable='box')
ax.text(.5,1.015,f'n = {len(bridge)}; rho = {cand.bridge_rho:.3f}; P = {cand.bridge_P:.3f}',fontsize=8,ha='center',va='bottom',transform=ax.transAxes)
ax.text(.5,-.33,'Strong fidelity not established',fontsize=8,ha='center',va='top',transform=ax.transAxes,color='#5E6462')
fig.text(.075,.026,'Retrospective candidate association; generated expression is not external validation',fontsize=8,color='#5E6462',ha='left',va='bottom')

fig.canvas.draw();renderer=fig.canvas.get_renderer();canvas=fig.bbox
texts=[];outside=[]
for t in fig.findobj(plt.matplotlib.text.Text):
 if t.get_visible() and t.get_text():
  bb=t.get_window_extent(renderer);item={'text':t.get_text(),'fontsize_pt':float(t.get_fontsize()),'bounds_px':list(bb.bounds)};texts.append(item)
  if bb.x0<canvas.x0-.5 or bb.y0<canvas.y0-.5 or bb.x1>canvas.x1+.5 or bb.y1>canvas.y1+.5:outside.append(item)
assert not outside,{'clipped_text':outside}
assert min(t['fontsize_pt'] for t in texts)>=8
metadata={'panels':panels,'source_hashes':{k:sha(p) for k,p in paths.items()},'code_sha256':sha(__file__),'font_minimum_pt':min(t['fontsize_pt'] for t in texts),'canvas_clipped_text':outside,'rendered_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scientific_checks':'All displayed correlations/n/cohort values and patient points asserted against staged canonical tables; full9type permutation BH not recomputed.'}
# Keep a final-size screen preview and grayscale QA companion without resizing exports.
FIGURES.mkdir(exist_ok=True)
fig.savefig(QA/'Figure_8.qa_final_size.png',dpi=96)
from PIL import Image
Image.open(QA/'Figure_8.qa_final_size.png').convert('L').save(QA/'Figure_8.qa_grayscale.png')
(QA/'Figure_8.qa.json').write_text(json.dumps({'status':'programmatic_pass_pending_visual_review','metadata':metadata,'text_bounds':texts},ensure_ascii=False,indent=2),encoding='utf-8')
record=save_figure(fig,'Figure_8',metadata=metadata)
print(json.dumps({k:record[k] for k in ['width_mm','height_mm']},indent=2))
print('Figure8 exported; all scientific/coverage and text-clipping assertions passed.')
