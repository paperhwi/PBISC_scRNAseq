"""Build five publication supplements from this package's frozen CSVs only.

No association, adjustment, confidence interval, or gene-score model is fitted.
Only table selection, counts of existing flags, and display transformations occur.
"""
from pathlib import Path
import json,hashlib,datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from plot_style import ROOT,DATA,FIGURES,MM,COLORS,MARKERS,COHORT_ORDER,setup_style,save_figure

setup_style()
SRC=ROOT/'source_data';SRC.mkdir(exist_ok=True)
SPECS=ROOT/'supplementary_specs';SPECS.mkdir(exist_ok=True)
QA=ROOT/'qa'/'supplementary';QA.mkdir(parents=True,exist_ok=True)
INPUTS={};RECORDS=[];CAPTIONS=[]
INK='#171717';MUTED='#5E6462';RULE='#D8D8D4';TINT='#E7EFEC'

def read(rel):
    p=DATA/rel
    INPUTS[str(rel)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p)

def source(df,name):
    df.to_csv(SRC/(name+'.csv'),index=False)

def lab(s):
    return (str(s).replace('B_ASC_','ASC ').replace('CD14mono','CD14 mono')
            .replace('CD16mono','CD16 mono').replace('Int.mono','Int. mono')
            .replace('T CD4 ','CD4 ').replace('T CD8 ','CD8 ')
            .replace('T g/d','T γδ').replace('_',' ').replace('IFN stim','IFN'))

def qfmt(v):
    if pd.isna(v):return '—'
    if v<.001:return f'{v:.1e}'
    if v<.1:return f'{v:.3f}'
    return f'{v:.2f}'

def pfmt(v):
    if pd.isna(v):return '—'
    return f'{v:.1e}' if v<.001 else f'{v:.3f}'

def fignew(height,title):
    fig=plt.figure(figsize=(178*MM,height*MM))
    fig.text(8/178,(height-7)/height,title,fontsize=10,fontweight='semibold',va='top')
    return fig

def axes_mm(fig,x,y,w,h):
    W,H=fig.get_size_inches()*25.4
    return fig.add_axes([x/W,y/H,w/W,h/H])

def text_mm(fig,x,y,s,**kw):
    W,H=fig.get_size_inches()*25.4
    kw.setdefault('fontsize',8)
    return fig.text(x/W,y/H,s,**kw)

def title(ax,letter,s):
    ax.set_title(f'{letter}  {s}',loc='left',fontsize=9,fontweight='semibold',pad=9)

def cohort_legend(fig,y,height):
    handles=[Line2D([],[],color=COLORS[a],marker=MARKERS[a],lw=.8,ms=4,label=a) for a in COHORT_ORDER]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.58,y/height),ncol=3,handlelength=1.4,columnspacing=1.8,frameon=False)

def record(fig,stem,spec,caption):
    """Retain exact final geometry and diagnostic previews; do not mask failures."""
    (SPECS/(stem+'.json')).write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
    fig.canvas.draw();renderer=fig.canvas.get_renderer();fw,fh=fig.canvas.get_width_height()
    boxes=[];clipped=[];not_drawn=set()
    # Matplotlib retains Text objects for axis-off and out-of-view ticks even
    # though Axis.draw never renders them. Exclude only those defined cases.
    for ax in fig.axes:
        for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            low,high=sorted(limits)
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not ax.axison or not axis.get_visible() or tick.get_loc()<low-1e-10 or tick.get_loc()>high+1e-10:
                    not_drawn.update([tick.label1,tick.label2])
            if not ax.axison or not axis.get_visible():not_drawn.update([axis.label,axis.offsetText])
    for t in fig.findobj(plt.matplotlib.text.Text):
        if t in not_drawn or not(t.get_visible() and t.get_text()):continue
        bb=t.get_window_extent(renderer)
        if bb.width==0 or bb.height==0:continue
        entry=dict(text=t.get_text(),bounds=[float(bb.x0),float(bb.y0),float(bb.x1),float(bb.y1)],font=float(t.get_fontsize()))
        boxes.append(entry)
        if bb.x0<-.7 or bb.y0<-.7 or bb.x1>fw+.7 or bb.y1>fh+.7:clipped.append(entry)
    overlaps=[]
    for i,a in enumerate(boxes):
        for b in boxes[i+1:]:
            aa=a['bounds'];bb=b['bounds'];w=min(aa[2],bb[2])-max(aa[0],bb[0]);h=min(aa[3],bb[3])-max(aa[1],bb[1])
            if w>1.2 and h>1.2:
                overlaps.append(dict(a=a['text'],b=b['text'],overlap_pixels=[w,h]))
    # Overlap list is a candidate diagnostic; data marks are inspected visually.
    report=dict(stem=stem,width_mm=178,height_mm=float(fig.get_figheight()*25.4),min_font_pt=min(z['font'] for z in boxes),clipped_text=clipped,text_overlap_candidates=overlaps)
    (QA/(stem+'_geometry.json')).write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    r=save_figure(fig,stem,metadata=spec)
    from PIL import Image
    im=Image.open(FIGURES/(stem+'.png')).convert('RGB')
    # 96 dpi is the final physical-width view used for screen inspection.
    preview=im.resize((673,round(im.height*673/im.width)),Image.Resampling.LANCZOS)
    preview.save(QA/(stem+'_96dpi.png'));preview.convert('L').save(QA/(stem+'_96dpi_gray.png'))
    RECORDS.append(report);CAPTIONS.append((stem,caption))

def supplement1():
    A=read('composition/all43_associations.csv')
    C=read('composition/all43_fixed_score_bootstrap_CI.csv')
    G=read('composition/all43_gain_loss_and_support.csv')
    R=read('composition/rarity_all43.csv')
    # Frozen fine-label order groups major lineages; never sort by significance.
    features=A[A.arm.eq('Real')].feature.to_list();assert len(features)==43
    rare=R.set_index('feature').rare_primary
    source(A,'S1_all43_associations');source(C,'S1_all43_confidence_intervals');source(G,'S1_all43_gain_loss_and_rarity')
    H=230;fig=fignew(H,'Complete composition screen retains gains, losses and null results')
    cohort_legend(fig,215,H)
    axs=[axes_mm(fig,48,48,55,159),axes_mm(fig,115,48,55,159)]
    for ax,endpoint,metric,lo,hi,null,letter in [(axs[0],'severity','rho','rho_CI_low','rho_CI_high',0,'a'),(axs[1],'death','AUC','AUC_CI_low','AUC_CI_high',.5,'b')]:
        title(ax,letter,'Peak WHO severity' if endpoint=='severity' else 'Recorded mortality')
        for j,f in enumerate(features):
            y=42-j
            if j%2==0:ax.axhspan(y-.49,y+.49,color='#F7F8F7',zorder=0)
            for arm,off in zip(COHORT_ORDER,[.22,0,-.22]):
                d=A[A.arm.eq(arm)&A.feature.eq(f)].iloc[0];z=C[C.arm.eq(arm)&C.feature.eq(f)].iloc[0]
                q=d['q43_'+endpoint];yy=y+off
                valid=d['severity_status' if endpoint=='severity' else 'death_status']=='estimable'
                if valid:
                    x=d[metric]
                    if np.isfinite(z[lo]) and np.isfinite(z[hi]):ax.plot([z[lo],z[hi]],[yy,yy],color=COLORS[arm],lw=.55,alpha=.78,zorder=2)
                    ax.plot(x,yy,marker=MARKERS[arm],ms=2.7,mec=COLORS[arm],mew=.65,mfc=COLORS[arm] if q<.05 else 'white',ls='',zorder=3)
                else:
                    # Undefined rho is not displayed as an estimated biological zero.
                    ax.plot(null,yy,marker='|',ms=3,color=COLORS[arm],mew=.8,zorder=3)
        ax.axvline(null,color='#9DA49F',lw=.7,ls='--',zorder=1)
        ax.set_ylim(-.7,42.7);ax.set_yticks(range(43));ax.tick_params(axis='y',length=0)
        if endpoint=='severity':
            ax.set_xlim(-1,1);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_xlabel('Spearman ρ')
            ax.set_yticklabels([lab(f)+(' †' if rare[f] else '') for f in features[::-1]],fontsize=8)
        else:
            ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.set_xlabel('AUC · higher in deaths →');ax.set_yticklabels([])
        ax.spines['left'].set_visible(False)
    # Exact transition counts accompany all-source effect and CI plots.
    ctx=axes_mm(fig,8,10,162,11);ctx.set_axis_off()
    ctx.text(0,1.10,'c  Threshold transitions at endpoint-wise BH43 q < .05',fontsize=9,fontweight='semibold',transform=ctx.transAxes)
    for y,ep,label in [(.58,'severity','Severity'),(.05,'death','Mortality')]:
        z=G[G.endpoint.eq(ep)]
        vals=[int((z.Real_q43<.05).sum()),int((z.Combined_q43<.05).sum()),int(z.retained.sum()),int(z.FDR_rescue.sum()),int(z.lost.sum())]
        ctx.text(0,y,label,fontsize=8.5,transform=ctx.transAxes)
        ctx.text(.22,y,f'Real {vals[0]} → Combined {vals[1]}     retained {vals[2]}     gained {vals[3]}     lost {vals[4]}',fontsize=8.5,transform=ctx.transAxes)
    text_mm(fig,8,36,'Filled marks: q43 < .05; open: q43 ≥ .05; vertical tick: not estimable.',va='center')
    text_mm(fig,8,31,'† Rare: Real mean and median <1% of leukocytes. All 43 remain in each correction family.',va='center')
    text_mm(fig,8,4,'Intervals: 50,000 fixed-score patient bootstraps; source-stratified ρ, source × outcome AUC.',va='bottom')
    spec=dict(claim='Full composition results include both gains and losses after cohort expansion.',visual_form='two aligned all43 dot-and-interval panels plus exact transition counts',metric=['Spearman rho','AUC'],unit='patient-level source rankINT associations',direction='AUC high means higher in deaths; never flipped',cohort_order=COHORT_ORDER,dimensions_mm=[178,H],uncertainty='50000 fixed-score patient percentile bootstrap; rho source-stratified/AUC source-by-outcome',multiplicity='BH43 per endpoint per cohort; q<.05',rarity='Real patient mean and median below1% of leukocytes, applied after full43 correction',nonestimable='vertical tick at reference, no fabricated CI',source_data=['S1_all43_associations.csv','S1_all43_confidence_intervals.csv','S1_all43_gain_loss_and_rarity.csv'],reference_grammar=['Astra uncertainty comparator','shared-axis forest plots','style_contract dot-and-interval'])
    record(fig,'Supplementary_Figure_S1',spec,'All 43 leukocyte-label associations. Patient scores are equal-visit means followed by within-source rankINT; the full fixed43 family is retained in Real37/New29/Combined66. Marks identify each source and filling indicates endpoint-specific BH43 q<.05. Whiskers are fixed-score 95% percentile patient-bootstrap intervals (50,000 replicates; source strata for rho and source×outcome strata for AUC). A vertical tick denotes an undefined constant/nonestimable association. Rarity (dagger) requires both Real patient mean and median abundance below1% after RBC/platelet exclusion and leukocyte closure. Transition counts are threshold changes, not formal augmentation-effect tests. All null/opposite and lost signals remain shown.')

def supplement2():
    V=read('transcriptome/all14_primary_coverage_and_composition_rarity.csv')
    N=read('transcriptome/screen_counts_all_configurations.csv');N=N[N.configuration.eq('primary')].copy()
    # Summarize existing frozen q columns for every cohort; no test is rerun.
    full=read('transcriptome/primary_all_gene_associations.csv.gz')
    count_rows=[]
    for (arm,celltype),z in full.groupby(['cohort','celltype'],sort=False):
        for ep in ['severity','death']:
            count_rows.append(dict(cohort=arm,celltype=celltype,endpoint=ep,genes=len(z),q_local05=int((z['q_gene_'+ep]<.05).sum()),q_global05=int((z['q_global_'+ep]<.05).sum()),q_global_joint05=int((z['q_global_joint_'+ep]<.05).sum())))
    allcounts=pd.DataFrame(count_rows);source(allcounts,'S2_all9_all_cohort_screen_counts')
    joint_n=int(allcounts.loc[allcounts.cohort.eq('Combined'),'q_global_joint05'].sum())
    assert joint_n==0
    assert len(V)==14 and int(V.eligible.sum())==9
    source(V,'S2_all14_coverage_rarity');source(N,'S2_all9_primary_screen_counts')
    H=192;fig=fignew(H,'RNA coverage and complete shared-gene testing define the discovery scope')
    ax=axes_mm(fig,8,100,162,73);ax.set_axis_off();title(ax,'a','Coverage of all 14 common labels · patient N / deaths')
    xs=[0,.40,.58,.77,.94];heads=['Cell label','Real','New','Eligible','Rare†']
    for x,s in zip(xs,heads):ax.text(x,1,s,fontsize=8.5,fontweight='semibold',ha='left' if x==0 else 'center',transform=ax.transAxes)
    ax.plot([0,1],[.965,.965],transform=ax.transAxes,color=RULE,lw=.7)
    for i,z in V.reset_index(drop=True).iterrows():
        y=.91-i*.067
        vals=[lab(z.celltype),f'{int(z.Real_n)} / {int(z.Real_deaths)}',f'{int(z.New_n)} / {int(z.New_deaths)}','Yes' if z.eligible else 'No','Yes' if z.rare_primary else 'No']
        for x,s in zip(xs,vals):ax.text(x,y,s,fontsize=8,ha='left' if x==0 else 'center',transform=ax.transAxes,color=INK if z.eligible else MUTED)
    labs=V[V.eligible].celltype.to_list();a=axes_mm(fig,8,28,162,42);a.set_axis_off()
    text_mm(fig,8,80,'b  Combined discoveries among all 9 eligible labels',fontsize=9,fontweight='semibold',va='bottom')
    xx=[0,.39,.57,.77,.94]
    for x,s in zip(xx,['Cell label','Severity\nlocal','Severity\nglobal','Death\nlocal','Death\nglobal']):a.text(x,1.01,s,fontsize=8,fontweight='semibold',ha='left' if x==0 else 'center',va='bottom',transform=a.transAxes)
    a.plot([0,1],[.99,.99],transform=a.transAxes,color=RULE,lw=.7)
    for i,l in enumerate(labs):
        y=.91-i*.099;sv=N[N.celltype.eq(l)&N.endpoint.eq('severity')].iloc[0];de=N[N.celltype.eq(l)&N.endpoint.eq('death')].iloc[0]
        vals=[lab(l)+(' †' if bool(V.set_index('celltype').loc[l,'rare_primary']) else ''),str(int(sv.Combined_q_gene05)),str(int(sv.Combined_q_global05)),str(int(de.Combined_q_gene05)),str(int(de.Combined_q_global05))]
        for x,s in zip(xx,vals):a.text(x,y,s,fontsize=8.5 if x else 8,ha='left' if x==0 else 'center',transform=a.transAxes,color=COLORS['Combined'] if x and s!='0' else INK,fontweight='semibold' if x and s!='0' else 'normal')
    text_mm(fig,8,91,'Eligible: each arm ≥10 patients, ≥3 deaths and ≥3 survivors; each patient has ≥3 cells.',va='center')
    text_mm(fig,8,86,'Coverage failure is not a biological null. † Real mean and median abundance <1%.',va='center')
    text_mm(fig,8,17,'Local: BH5,762 within label/endpoint. Global: BH51,858 across 9 labels/endpoint.',va='center')
    text_mm(fig,8,12,f'Global joint: BH103,716 across both endpoints; Combined discoveries = {joint_n}.',va='center')
    text_mm(fig,8,5,'The 9 eligible labels include 7 rare labels; all common genes and eligible labels remain in scope.',va='bottom')
    spec=dict(claim='Coverage and multiplicity scope constrain available RNA discovery.',visual_form='two aligned exact count tables',metric='patients/deaths and count of q<.05 gene tests',unit='patient and gene-label hypotheses',cohort_order=COHORT_ORDER,dimensions_mm=[178,H],multiplicity=['BH5762 each celltype/endpoint','BH51858 all9labels/endpoint','BH103716 joint endpoints sensitivity'],nonestimable='ineligible labels retained in coverage with explicit No; never biological null',source_data=['S2_all14_coverage_rarity.csv','S2_all9_primary_screen_counts.csv','S2_all9_all_cohort_screen_counts.csv'],reference_grammar=['Astra sparse count table','booktabs without vertical rules'])
    record(fig,'Supplementary_Figure_S2',spec,'Coverage of every common label and complete primary RNA screen. Patients require at least three cells of the label across all eligible visits. Eligibility additionally requires ≥10 patients, ≥3 deaths and ≥3 survivors in each source; nine labels qualify and seven satisfy the separately defined rarity criterion. Counts in panel b are Combined gene tests passing q<.05, including every qualifying label and all null counts. Local BH uses5,762 common symbols within each label/endpoint; global BH uses9×5,762=51,858 per endpoint. No Combined test passes the global joint103,716 family. Scores are patient cell-weighted means of log1p(CP10K) over the shared genes and source-wise rankINT. Coverage failures are not biological negative findings, and rare-subset filenames do not make every label rare.')

def supplement3():
    F=read('transcriptome/fixed_primary_candidates_with_rarity_and_composition.csv')
    C=read('transcriptome/fixed_primary_candidate_CIs_50000.csv')
    targets=[('ADK','severity'),('MXRA7','death')]
    focus=F[F.celltype.eq('T g/d')&F.gene.isin([z[0] for z in targets])].copy()
    cis=C[C.celltype.eq('T g/d')&((C.gene.eq('ADK')&C.endpoint.eq('severity'))|(C.gene.eq('MXRA7')&C.endpoint.eq('death')))].copy()
    source(focus,'S3_Tgamma_delta_FDR_rescue_effects');source(cis,'S3_Tgamma_delta_confidence_intervals')
    H=155;fig=fignew(H,'Tγδ RNA signals support FDR rescue of existing nominal associations')
    for k,(gene,ep) in enumerate(targets):
        z=F[F.celltype.eq('T g/d')&F.gene.eq(gene)&F.endpoint.eq(ep)].iloc[0]
        ax=axes_mm(fig,28+89*k,88,53,40);title(ax,chr(97+k),gene+(' · severity' if ep=='severity' else ' · mortality'))
        for i,arm in enumerate(COHORT_ORDER):
            r=cis[cis.gene.eq(gene)&cis.endpoint.eq(ep)&cis.cohort.eq(arm)].iloc[0];y=2-i
            ax.plot([r.selected_endpoint_lo,r.selected_endpoint_hi],[y,y],color=COLORS[arm],lw=1)
            ax.plot(r.selected_endpoint_effect,y,marker=MARKERS[arm],ms=4.5,color=COLORS[arm],ls='')
        ax.set_yticks([2,1,0],[f'{a}\nN={int(z[a+"_n"])}' for a in COHORT_ORDER]);ax.set_ylim(-.5,2.5)
        ax.set_xlim((-1,.1) if ep=='severity' else (0,1));ax.set_xticks([-1,-.5,0] if ep=='severity' else [0,.5,1]);ax.axvline(0 if ep=='severity' else .5,color=RULE,lw=.8,ls='--')
        ax.set_xlabel('Spearman ρ' if ep=='severity' else 'AUC · higher in deaths →');ax.tick_params(axis='y',length=0)
        t=axes_mm(fig,8+89*k,31,73,33);t.set_axis_off()
        for x,s in zip([0,.29,.57,.88],['Cohort','P','q local','q global']):t.text(x,1,s,fontsize=8,fontweight='semibold',ha='left' if x==0 else 'center',transform=t.transAxes)
        t.plot([0,1],[.91,.91],color=RULE,lw=.6,transform=t.transAxes)
        for i,arm in enumerate(COHORT_ORDER):
            y=.70-i*.27
            for x,s in zip([0,.29,.57,.88],[arm,pfmt(z[arm+'_P']),qfmt(z[arm+'_q_gene']),qfmt(z[arm+'_q_global'])]):t.text(x,y,s,fontsize=8,ha='left' if x==0 else 'center',transform=t.transAxes)
    text_mm(fig,8,76,'Local: BH5,762 per label/endpoint. Global: BH51,858 across all 9 eligible labels.',va='center')
    text_mm(fig,8,20,'Real P < .05 for both targets: these are FDR rescues, not strict Real-missed results.',va='center')
    text_mm(fig,8,14,'MXRA7: lower values in deaths; New has only 3 deaths (N=24); global q ≥ .05.',va='center')
    text_mm(fig,8,5,'95% CI: 50,000 fixed-score patient bootstraps; ρ by source, AUC by source × outcome.',va='bottom')
    spec=dict(claim='Tgamma-delta ADK and MXRA7 are FDR rescues, with sparse New mortality support.',visual_form='two cohort forests plus exact P/local/global q tables',metric=['Spearman rho','AUC'],unit='patient-level source rankINT gene scores',dimensions_mm=[178,H],cohort_order=COHORT_ORDER,uncertainty='50000 fixed-score patient bootstrap95% percentile CIs',multiplicity='local BH5762; endpoint-global BH51858',direction='ADK negative severity rho; MXRA7 AUC below.5 means lower values in deaths, no flip',source_data=['S3_Tgamma_delta_FDR_rescue_effects.csv','S3_Tgamma_delta_confidence_intervals.csv'],reference_grammar=['Astra forest comparator','sparse exact value table'])
    record(fig,'Supplementary_Figure_S3',spec,'Source-specific Tγδ ADK severity and MXRA7 mortality associations with50,000 fixed-score patient-bootstrap95% CIs. Both genes were nominally associated in Real data before augmentation, making them FDR rescues. Local q corrects5,762 common genes within label and endpoint; global q corrects51,858 label–gene pairs within endpoint. ADK additionally passes endpoint-global correction, whereas MXRA7 does not. MXRA7 AUCs are retained in their original direction: lower expression scores occur in deaths. The New Tγδ cohort includes24 patients and only3 deaths. Candidate selection and generation uncertainty are not included in the CIs; support from discovery New patients is not independent replication.')

def supplement4():
    B=read('plotting/bulk_KIAA0391_patients.csv');P=read('plotting/KIAA0391_shared_bridge.csv')
    R=read('verification/Bulk_KIAA0391_context.csv');F=read('transcriptome/fixed_primary_candidates_with_rarity_and_composition.csv')
    ref=R[R.normalization.eq('all_raw_Ensembl')&R.aggregation.eq('equal_draw_mean')].iloc[0]
    kr=F[F.celltype.eq('HSPC')&F.gene.eq('KIAA0391')&F.endpoint.eq('severity')].iloc[0]
    assert len(B)==29 and len(P)==24
    source(B,'S4_same29_bulk_HSPC_patient_scores');source(P,'S4_shared24_HSPC_patient_scores');source(R,'S4_bulk_context_all_normalizations')
    H=173;fig=fignew(H,'Bulk context and shared-patient agreement address distinct questions')
    ax=axes_mm(fig,24,101,60,46);title(ax,'a','Bulk gene score versus severity')
    jitter=np.linspace(-.075,.075,len(B));order=np.argsort(B.patient_id.astype(str).to_numpy());jj=np.empty(len(B));jj[order]=jitter
    ax.scatter(B.peak_severity+jj,B.all_raw_Ensembl__equal_draw_mean,s=16,marker=MARKERS['New'],facecolors='white',edgecolors=COLORS['New'],linewidths=.7)
    ax.set_xticks([0,1,2,3]);ax.set_xlim(-.35,3.35);ax.set_xlabel('Peak severity (recoded 0–3)');ax.set_ylabel('Bulk log2(CPM + 1)\npatient mean')
    ax.text(.02,.98,f'N=29; ρ={ref.rho_bulk_vs_peak_severity:.3f}\nP={pfmt(ref.P_bulk_vs_peak_severity)}',transform=ax.transAxes,va='top',fontsize=8)
    ax=axes_mm(fig,112,101,58,46);title(ax,'b','Bulk versus inferred HSPC score')
    ax.scatter(B.all_raw_Ensembl__equal_draw_mean,B.PBISC_HSPC_primary_cellmean,s=16,marker=MARKERS['New'],facecolors='white',edgecolors=COLORS['New'],linewidths=.7)
    ax.set_xlabel('Bulk log2(CPM + 1)\npatient mean');ax.set_ylabel('Inferred HSPC\nmean ln(CP10K + 1)')
    ax.text(.02,.98,f'N=29; ρ={ref.rho_bulk_vs_HSPC_primary:.3f}\nP={pfmt(ref.P_bulk_vs_HSPC_primary)}',transform=ax.transAxes,va='top',fontsize=8)
    ax=axes_mm(fig,24,25,60,46);title(ax,'c','Shared-patient HSPC expression')
    ax.scatter(P.Real_expression,P.PBISC_expression,s=17,color=COLORS['Combined'],marker='o',edgecolors='white',linewidths=.4)
    ax.set_xlabel('Observed Real HSPC\nmean ln(CP10K + 1)');ax.set_ylabel('Inferred PBISC HSPC\nmean ln(CP10K + 1)')
    ax.text(.98,.98,f'N=24; ρ={kr.bridge_rho:.3f}\nP={pfmt(kr.bridge_P)}',transform=ax.transAxes,va='top',ha='right',fontsize=8)
    note=axes_mm(fig,104,20,66,54);note.set_axis_off();title(note,'d','Interpretation boundaries')
    lines=['Same New29 / 79 draws for bulk.',
           'Bulk: all60,605 Ensembl rows;',
           'equal draw weights per patient.',
           '',
           'Shared scores summarize all visits.',
           'These are not matched-day pairs;',
           'individual timing is unavailable here.',
           '',
           'Non-significant bulk correlation is',
           'not proof of absence or specificity.',
           'No formal correlation-difference test.']
    note.text(0,.94,'\n'.join(lines),transform=note.transAxes,va='top',fontsize=8,linespacing=1.32)
    text_mm(fig,8,5,'All plotted points are patients. No fitted line, same-day pairing or unreported uncertainty is imputed.',va='bottom')
    spec=dict(claim='A weak single-gene bulk correlation and low shared-patient fidelity do not establish cell specificity.',visual_form='three patient scatterplots and explicit scope note',metric='raw normalized patient gene scores and peak WHO',unit='bulk log2(CPM+1); HSPC mean natural-log(CP10K+1)',dimensions_mm=[178,H],uncertainty='No resampling CI provided for these scatter correlations; raw patient values shown',cohort='sameNew29/79draws; shared24 HSPC patient means',timing='No timing fields staged; no same-day claims or invented lag distribution',source_data=['S4_same29_bulk_HSPC_patient_scores.csv','S4_shared24_HSPC_patient_scores.csv','S4_bulk_context_all_normalizations.csv'],reference_grammar=['Astra sparse scatter','explicit raw point and unit labels'])
    record(fig,'Supplementary_Figure_S4',spec,'Bulk KIAA0391 context in exactly the same29 New patients and79 eligible draws, and all-visit shared-patient HSPC scores in24 patients. Severity is recoded as0=mild,1=moderate,2=severe,3=critical, not the WHO0–8 clinical progression scale. Bulk values are draw-level log2(CPM+1), using all60,605 unique Ensembl rows in the denominator, then averaged with equal draw weights within patient. HSPC values are cell-weighted patient means of log1p(CP10K) normalized over the shared5,762-symbol universe. Panel a uses deterministic horizontal jitter solely to show tied ordinal severity; no statistic is recomputed from jittered coordinates. Panel b compares different bulk and inferred-cell estimands. Panel c is not a same-day validation: exact lag data are not available in the plotted package. Correlations/P are existing descriptive canonical statistics; no CI was available for these scatter summaries. Non-significance does not establish an absent bulk effect, cell specificity, or superiority; negative shared concordance does not by itself distinguish biological timing from reconstruction error.')

def q_table(ax,rows,rowlabels,columns,header_groups=None,row_height=None):
    """Sparse numeric table with no vertical rules; precise q values remain primary."""
    ax.set_axis_off();n=len(rows);nc=len(columns);left=.32;centers=np.linspace(left+.05,.975,nc)
    ax.set_xlim(0,1);ax.set_ylim(-.65,n+1.1)
    for j,(x,col) in enumerate(zip(centers,columns)):ax.text(x,n+.23,col,fontsize=8,ha='center',va='bottom',fontweight='semibold')
    ax.plot([0,1],[n+.06,n+.06],color=RULE,lw=.7)
    for i,(name,values) in enumerate(zip(rowlabels,rows)):
        y=n-1-i
        ax.text(0,y,name,fontsize=8,ha='left',va='center')
        for x,v in zip(centers,values):
            if np.isfinite(v) and v<.05:ax.add_patch(Rectangle((x-.044,y-.44),.088,.88,facecolor=TINT,edgecolor='none',zorder=0))
            ax.text(x,y,qfmt(v),fontsize=8,ha='center',va='center',fontweight='semibold' if np.isfinite(v) and v<.05 else 'normal')
    ax.plot([0,1],[-.55,-.55],color=RULE,lw=.7)

def supplement5():
    M=read('modules/HSPC_all_module_associations.csv');H=read('transcriptome/historical_HSPC20_both_endpoints.csv')
    p=M[M.configuration.eq('primary')&M.module.eq(M.representative)].copy()
    names=p[p.cohort.eq('Real')].module.to_list();assert len(names)==16
    genes=['GATA1','KLF1','TFRC','BLVRB','APOC1','GATA2','ITGA2B','MKI67','TOP2A','STMN1','TYMS','CD34','AVP','HLF','HOPX','CRHBP','MPO','ELANE','CSF3R','CEBPD']
    hh=H[H.cohort.eq('Combined')].set_index('gene').loc[genes].reset_index()
    source(M,'S5_all_HSPC_module_associations');source(H,'S5_historical20_all_cohorts')
    Ht=223;fig=fignew(Ht,'Module and historical-marker results depend on their declared families')
    a=axes_mm(fig,8,134,162,65);title(a,'a','All 16 unique HSPC modules · endpoint-specific q16')
    a.set_title('a  All 16 unique HSPC modules · endpoint-specific q16',loc='left',fontsize=9,fontweight='semibold',pad=18)
    vals=[]
    for name in names:
        z=p[p.module.eq(name)].set_index('cohort')
        vals.append([z.loc[arm,'q_unique16_'+ep] for ep in ['severity','death'] for arm in COHORT_ORDER])
    q_table(a,vals,[lab(n).replace('Alarmin Inflam','Alarmin/inflam.') for n in names],['Real','New','Comb.','Real','New','Comb.'])
    text_mm(fig,87.5,200,'Severity',ha='center',va='center',fontweight='semibold')
    text_mm(fig,146.3,200,'Mortality',ha='center',va='center',fontweight='semibold')
    b=axes_mm(fig,8,30,162,81);title(b,'b','Historical 20-marker panel · Combined q20 versus all-gene q5,762')
    b.set_title('b  Historical 20-marker panel · Combined q20 versus all-gene q5,762',loc='left',fontsize=9,fontweight='semibold',pad=18)
    vv=hh[['q_target20_severity','q_gene_severity','q_target20_death','q_gene_death']].to_numpy()
    q_table(b,vv,genes,['q20','q5,762','q20','q5,762'])
    text_mm(fig,84.3,112,'Severity',ha='center',va='center',fontweight='semibold')
    text_mm(fig,149.6,112,'Mortality',ha='center',va='center',fontweight='semibold')
    ery=p[p.module.eq('Ery5')&p.cohort.eq('Combined')].iloc[0]
    text_mm(fig,8,126,f'Ery5 severity q16={ery.q_unique16_severity:.3f}; joint32 q={ery.q_joint32_severity:.3f}. Joint results are supplied in full CSV.',va='center')
    text_mm(fig,8,20,'Tint/bold: q < .05 within the labeled family. No historical marker passes HSPC-wide q5,762.',va='center')
    text_mm(fig,8,14,'Cycling and Cycle4 share identical available gene weights; the duplicate is counted only once.',va='center')
    text_mm(fig,8,5,'Modules and markers were selected retrospectively. Their q values are not interchangeable with global-gene q.',va='bottom')
    spec=dict(claim='Historical module/marker FDR is narrower than available-gene discovery.',visual_form='two exact q tables with sparse threshold tint',metric='BH adjusted P values',unit='q',dimensions_mm=[178,Ht],multiplicity=['16 unique modules percohort/endpoint','historical20 markers percohort/endpoint','all5762 sharedHSPCgenes percohort/endpoint'],alias='Cycling and Cycle4 have exactly identical available gene weights, one representative retained',cohort_order=COHORT_ORDER,uncertainty='No new effect uncertainty implied by q; complete effects andexistingCIs retained in package tables',source_data=['S5_all_HSPC_module_associations.csv','S5_historical20_all_cohorts.csv'],reference_grammar=['Astra sparse table','booktabs','family-specific exact values'])
    record(fig,'Supplementary_Figure_S5',spec,'All16 unique HSPC module hypotheses and the complete historical20-marker panel. Panel a shows primary endpoint-wise BH16 q values for Real33/New29/Combined62, retaining null results and the complete module set. Cycling and Cycle4 have identical available gene weights after absent shared genes are removed and are therefore treated as one exact alias; the full17-name and joint32 outputs remain in source data. Ery5 severity passes endpoint-wise q16 but not joint32. Panel b contrasts Combined historical BH20 with HSPC-wide BH5,762, separately for severity and mortality. The20 markers were selected retrospectively; none passes the displayed all-gene correction. Tint/bold denotes only q<.05 within the explicitly labeled family and does not claim external validation or selection-wide error control.')

def main():
    supplement1();supplement2();supplement3();supplement4();supplement5()
    (SPECS/'captions.md').write_text('\n\n'.join('## '+s+'\n\n'+c for s,c in CAPTIONS)+'\n',encoding='utf-8')
    (SPECS/'source_inputs.json').write_text(json.dumps(INPUTS,indent=2),encoding='utf-8')
    summary=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),figures=RECORDS,style='Shared package plot_style.py; all ordinary text≥8pt',inputs_sha256=INPUTS,scope='Only existing canonical package tables selected and displayed; no association or interval recalculation',visual_review='Pending visual inspection of latest color and grayscale renders')
    (SPECS/'SUPPLEMENTARY_QA.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps([dict(stem=r['stem'],clipped=len(r['clipped_text']),text_overlap_candidates=len(r['text_overlap_candidates'])) for r in RECORDS],indent=2))

if __name__=='__main__':main()
