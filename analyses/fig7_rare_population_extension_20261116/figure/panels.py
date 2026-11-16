"""Panel functions for the ESKD COVID clinical-association figures (figures4papers style, see style.py).
Each function draws one panel on the axis it is given; numbers come from data/*.csv (compute_fig_data.py, compute_nk_umap.py)."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from style import PALETTE, SRC_COL, SRC_MK, CT_COL
D='data/'
SRC_COL=dict(SRC_COL); SRC_COL['Combined']=PALETTE['teal']; SRC_COL['PBISC-only']=PALETTE['blue_secondary']
SRC_MK=dict(SRC_MK); SRC_MK['Combined']='s'; SRC_MK['PBISC-only']='o'
SEV=['mild','moderate','severe','critical']
DIV=LinearSegmentedColormap.from_list('div',[PALETTE['blue_main'],'#9DBBDD','#FFFFFF',PALETTE['red_2'],PALETTE['red_strong']])
SEQ=LinearSegmentedColormap.from_list('seq',['#F4F4F4',PALETTE['blue_secondary'],PALETTE['blue_main']])
ANALYSES=[('Real','Real'),('PBISC','PBISC'),('PBISC_only_patients_not_in_Real(independent)','PBISC-only'),('Combined_union_patients','Combined')]
FEAT=[('T_pct','T cell'),('B_pct','B cell'),('Mono_pct','Monocyte'),('NK_pct','NK cell'),('Other_pct','Other'),('DC_within_Other_pct','DC'),
      ('rare13_total_pct','Rare (13 total)'),('rare_NK CD56 bright_pct','NK CD56bright'),('rare_B_ASC_dividing_pct','B ASC dividing'),('rare_B_ASC_IgG_pct','B ASC IgG'),
      ('rare_HSPC_pct','HSPC'),('rare_DC1_pct','DC1'),('rare_ASDC_pct','ASDC'),('rare_pDC_pct','pDC'),('rare_CD14mono_IFN_pct','CD14 mono IFN'),
      ('rare_CD14mono_dex_pct','CD14 mono dex'),('rare_T g/d_pct','T γδ'),('rare_MAIT_pct','MAIT'),('rare_T CD4 cycling_pct','T CD4 cycling'),('rare_ILC_pct','ILC')]
FNAME=dict(FEAT)
def _assoc(): return pd.read_csv(D+'assoc_patient_level.csv')
def _pos(): return pd.read_csv(D+'patient_level_covid_positive.csv')
def fmt_q(q):
    if pd.isna(q): return 'n.a.'
    return f'{q:.0e}'.replace('e-0','e-') if q<0.001 else f'{q:.2g}'
def src_legend(ax, srcs, loc='upper left', ncol=1, bbox=None, ms=9, **kw):
    h=[Line2D([0],[0],marker=SRC_MK[s],color='none',markerfacecolor=SRC_COL[s],markeredgecolor='black',markeredgewidth=0.8,markersize=ms,label=s) for s in srcs]
    return ax.legend(handles=h,loc=loc,ncol=ncol,bbox_to_anchor=bbox,handletextpad=0.3,columnspacing=1.0,borderaxespad=0.2,**kw)

# ---------------- Figure 1 ----------------
def schematic(ax):
    ax.set_xlim(0,100); ax.set_ylim(0,20); ax.axis('off')
    def box(x,y,w,h,txt,fc,ec,fs=14,bold=False):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.3,rounding_size=1.2',fc=fc,ec=ec,lw=2))
        ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs,fontweight='bold' if bold else 'normal',linespacing=1.35)
    box(4,11,21,8,'Real scRNA-seq\n39 COVID+ patients\n149 samples, 588k cells','#EDEDED',SRC_COL['Real'])
    box(4,1,21,8,'PBISC from bulk\n57 COVID+ patients\n4,000 cells/sample','#E3ECF6',SRC_COL['PBISC'])
    box(30,6,17,8,'28 patients in both\n0 same-day pairs','white',PALETTE['gray_1'],fs=13)
    for y in (15,5): ax.annotate('',xy=(30,10),xytext=(25.8,y),arrowprops=dict(arrowstyle='-|>',lw=2,color=PALETTE['gray_2']))
    box(53,11,21,8,'Per source\nReal n=39 | PBISC n=57\nPBISC-only n=29','white',PALETTE['gray_2'],fs=13)
    box(53,1,21,8,'Combined (n=68)\nrank-INT within source,\nmean per patient','#E2F1F2',SRC_COL['Combined'],fs=13)
    ax.annotate('',xy=(52.6,15),xytext=(47.4,11),arrowprops=dict(arrowstyle='-|>',lw=2,color=PALETTE['gray_2']))
    ax.annotate('',xy=(52.6,5),xytext=(47.4,9),arrowprops=dict(arrowstyle='-|>',lw=2,color=PALETTE['gray_2']))
    box(80,11,19,8,'Cell-type fractions\nvs peak WHO severity\nand death','white',PALETTE['gray_2'],fs=13)
    box(80,1,19,8,'Rare-subset pseudobulk\nseverity DEG\nand subclusters','white',PALETTE['gray_2'],fs=13)
    for y in (15,5): ax.annotate('',xy=(79.6,y),xytext=(74.4,y),arrowprops=dict(arrowstyle='-|>',lw=2,color=PALETTE['gray_2']))

def assoc_dotgrid(ax, cax=None, sax=None):
    A=_assoc(); feats=[f for f,_ in FEAT]
    for j,(key,lab) in enumerate(ANALYSES):
        d=A[A.analysis==key].set_index('feature').reindex(feats)
        q=d.q_severity_trend.values; r=d.spearman_rho_severity.values
        s=np.clip(-np.log10(np.where(np.isnan(q),1,q)),0,5)
        y=np.arange(len(feats))
        sc=ax.scatter(np.full(len(feats),j),y,s=40+s*95,c=np.nan_to_num(r),cmap=DIV,norm=TwoSlopeNorm(0,-0.7,0.7),edgecolors=np.where(q<0.05,'black','#BBBBBB'),linewidths=np.where(q<0.05,2.0,0.8),zorder=3)
    ax.set_xticks(range(len(ANALYSES))); ax.set_xticklabels([f'{lab} (n={int(A[A.analysis==k].n_patients.iloc[0])})' for k,lab in ANALYSES],fontsize=13,rotation=35,ha='left',rotation_mode='anchor')
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([FNAME[f] for f in feats]); ax.set_ylim(len(feats)-0.5,-0.6); ax.set_xlim(-0.6,len(ANALYSES)-0.4)
    ax.axhline(6.5,color=PALETTE['neutral'],lw=1.2,ls='--',zorder=1)
    ax.spines['bottom'].set_visible(False); ax.tick_params(axis='x',length=0)
    ax.xaxis.tick_top()
    if cax is not None:
        cb=plt.colorbar(sc,cax=cax,orientation='horizontal'); cb.set_label('Spearman ρ with peak severity',fontsize=13); cb.ax.tick_params(labelsize=12); cb.outline.set_linewidth(1)
    if sax is not None:
        sax.axis('off')
        sax.text(-0.6,0,'BH q',fontsize=13,ha='right',va='center')
        for x,(v,l) in zip([0,1.75,3.6],[(1,'0.1'),(2,'0.01'),(4,'1e-4')]):
            sax.scatter(x,0,s=40+v*95,color='white',edgecolors='black',linewidths=1.0)
            sax.text(x+0.38,0,l,ha='left',va='center',fontsize=12)
        sax.scatter(5.6,0,s=40+2*95,color='white',edgecolors='black',linewidths=2.0); sax.text(5.95,0,'q < 0.05',ha='left',va='center',fontsize=12)
        sax.set_xlim(-0.9,7.6); sax.set_ylim(-1,1)

def strip_severity(ax, feat, ylabel=None, show_legend=False, title=None):
    P=_pos(); A=_assoc(); rng=np.random.default_rng(3)
    for k,src in enumerate(['Real','PBISC']):
        d=P[P.source==src]; off=-0.2 if src=='Real' else 0.2
        for i,s in enumerate(SEV):
            v=d.loc[d.WHO_severity==s,feat].values
            ax.scatter(i+off+rng.uniform(-0.07,0.07,len(v)),v,s=46,marker=SRC_MK[src],color=SRC_COL[src],alpha=0.85,edgecolors='white',linewidths=0.6,zorder=3)
            if len(v): ax.plot([i+off-0.13,i+off+0.13],[np.median(v)]*2,color='black',lw=2.6,zorder=4)
    ax.set_xticks(range(4)); ax.set_xticklabels(['Mild','Mod.','Severe','Crit.']); ax.set_xlim(-0.55,3.55)
    ax.set_ylabel(ylabel or f'{FNAME[feat]} (% of cells)')
    ymax=np.nanmax(P[feat].values); ax.set_ylim(bottom=0,top=ymax*1.38)
    for k,src in enumerate(['Real','PBISC']):
        r=A[(A.analysis==src)&(A.feature==feat)].iloc[0]
        ax.text(0.02+0.5*k,0.97,f'{src}: ρ={r.spearman_rho_severity:+.2f}\nq={fmt_q(r.q_severity_trend)}',transform=ax.transAxes,va='top',ha='left',fontsize=12.5,color=SRC_COL[src])
    if title: ax.set_title(title,loc='left',fontsize=14,pad=6)
    if show_legend: src_legend(ax,['Real','PBISC'],loc='upper right',ms=8,fontsize=12)

def power_dumbbell(ax):
    A=_assoc(); feats=['rare_NK CD56 bright_pct','rare_B_ASC_dividing_pct','rare_HSPC_pct','NK_pct','T_pct','rare_T g/d_pct','rare_MAIT_pct','DC_within_Other_pct','rare_DC1_pct','rare_ASDC_pct']
    keys=[('Real','Real'),('PBISC','PBISC'),('Combined_union_patients','Combined')]
    for i,f in enumerate(feats):
        xs=[]
        for k,lab in keys:
            q=A[(A.analysis==k)&(A.feature==f)].q_severity_trend.iloc[0]; x=-np.log10(q); xs.append(x)
        ax.plot([min(xs),max(xs)],[i,i],color=PALETTE['neutral'],lw=3,zorder=1)
        for (k,lab),x in zip(keys,xs):
            ax.scatter(x,i,s=95,marker=SRC_MK[lab],color=SRC_COL[lab],edgecolors='black',linewidths=0.8,zorder=3)
    ax.axvline(-np.log10(0.05),color=PALETTE['red_strong'],ls='--',lw=1.6,zorder=0); ax.set_xlim(-0.25,7.0)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([FNAME[f] for f in feats]); ax.set_ylim(len(feats)-0.4,-0.6)
    ax.set_xlabel('Severity association, −log$_{10}$ BH q')
    ax.axhspan(-0.5,2.5,color=PALETTE['green_1'],alpha=0.6,zorder=0,lw=0); ax.axhspan(6.5,9.5,color=PALETTE['red_1'],alpha=0.5,zorder=0,lw=0)
    ax.text(0.99,0.04,'Real-only\n(diluted)',transform=ax.transAxes,ha='right',va='bottom',fontsize=12,color=PALETTE['red_strong'])
    ax.text(0.99,0.97,'Same sign\nin both',transform=ax.transAxes,ha='right',va='top',fontsize=12,color='#2E7D32')
    h=src_legend(ax,['Real','PBISC','Combined'],loc='center right',ms=9,fontsize=12); hh=h.legend_handles+[Line2D([0],[0],color=PALETTE['red_strong'],ls='--',lw=1.6,label='q = 0.05')]; ax.legend(handles=hh,loc='center right',fontsize=12,handletextpad=0.4)

def roc_hspc(ax):
    R=pd.read_csv(D+'roc_hspc_fatal.csv'); A=pd.read_csv(D+'roc_hspc_fatal_auc.csv').set_index('analysis')
    ax.plot([0,1],[0,1],color=PALETTE['neutral'],ls='--',lw=1.6)
    for lab in ['Real','PBISC','Combined']:
        d=R[R.analysis==lab]; a=A.loc[lab]
        ax.step(d.fpr,d.tpr,where='post',color=SRC_COL[lab],lw=2.8,label=f'{lab} AUC {a.auc:.2f} ({int(a.n_fatal)}/{int(a.n)})')
    ax.set_xlim(-0.02,1.02); ax.set_ylim(-0.02,1.04); ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1]); ax.legend(loc='lower right',fontsize=12,handlelength=1.2); ax.set_title('HSPC % → death',loc='left',fontsize=14,pad=6)

# ---------------- Figure 2 ----------------
def reliability_dots(ax):
    N=pd.read_csv(D+'reliability_split_half_vs_cross_source.csv').set_index('feature')
    feats=[f for f,_ in FEAT if f in N.index and f!='rare_ILC_pct']
    sty=[('Real','Real split-half',SRC_COL['Real'],'D'),('PBISC','PBISC split-half',SRC_COL['PBISC'],'o'),('real_vs_pbisc_same_patient','Real vs PBISC,\nsame patient',PALETTE['red_strong'],'s')]
    for i,f in enumerate(feats):
        v=[N.loc[f,c] for c,_,_,_ in sty]; ax.plot([np.nanmin(v),np.nanmax(v)],[i,i],color=PALETTE['neutral'],lw=2.5,zorder=1)
        for (c,lab,col,mk),x in zip(sty,v): ax.scatter(x,i,s=80,marker=mk,color=col,edgecolors='black',linewidths=0.7,zorder=3,label=lab if i==0 else None)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([FNAME[f] for f in feats]); ax.set_ylim(len(feats)-0.5,-0.6); ax.set_xlim(-0.7,1.02)
    ax.set_xlabel('Patient-level Spearman ρ'); ax.legend(loc='upper left',fontsize=12,bbox_to_anchor=(0.0,1.0),borderaxespad=0.1)
    med=N.real_vs_pbisc_same_patient.median(); ax.text(0.01,0.01,f'median cross-source\nρ = {med:.2f}',transform=ax.transAxes,ha='left',va='bottom',fontsize=12,color=PALETTE['red_strong'])

def effect_concordance(ax):
    A=_assoc(); feats=[f for f,_ in FEAT if f!='rare_ILC_pct']
    r=A[A.analysis=='Real'].set_index('feature').reindex(feats); p=A[A.analysis=='PBISC'].set_index('feature').reindex(feats)
    cat=np.select([(r.q_severity_trend<.05)&(p.q_severity_trend<.05),r.q_severity_trend<.05,p.q_severity_trend<.05],['both','Real only','PBISC only'],'neither')
    cols={'both':PALETTE['teal'],'Real only':SRC_COL['Real'],'PBISC only':SRC_COL['PBISC'],'neither':PALETTE['neutral']}
    ax.axhline(0,color=PALETTE['gray_1'],lw=1); ax.axvline(0,color=PALETTE['gray_1'],lw=1); ax.plot([-0.8,0.8],[-0.8,0.8],color=PALETTE['neutral'],ls='--',lw=1.4)
    for c in ['neither','Real only','PBISC only','both']:
        m=cat==c; ax.scatter(r.spearman_rho_severity[m],p.spearman_rho_severity[m],s=110,color=cols[c],edgecolors='black',linewidths=0.8,label=f'q<0.05: {c}',zorder=3)
    lab={'rare_NK CD56 bright_pct':(8,-14),'DC_within_Other_pct':(8,6),'rare_HSPC_pct':(-10,8),'rare_B_ASC_dividing_pct':(8,-12),'T_pct':(8,-2),'NK_pct':(12,-17),'rare_DC1_pct':(8,6),'rare_T g/d_pct':(-12,-16)}
    for f,(dx,dy) in lab.items():
        ax.annotate(FNAME[f],(r.loc[f,'spearman_rho_severity'],p.loc[f,'spearman_rho_severity']),xytext=(dx,dy),textcoords='offset points',fontsize=12,ha='left' if dx>0 else 'right')
    rho=pd.Series(r.spearman_rho_severity.values).corr(pd.Series(p.spearman_rho_severity.values),method='spearman')
    ax.set_xlim(-0.75,0.75); ax.set_ylim(-0.75,0.75); ax.set_xlabel('Real: ρ with severity'); ax.set_ylabel('PBISC: ρ with severity')
    ax.text(0.03,0.97,f'across features ρ = {rho:.2f}',transform=ax.transAxes,va='top',fontsize=12.5)
    ax.legend(loc='lower right',fontsize=11.5,handletextpad=0.2,borderaxespad=0.1)

def nlr_scatter(ax, feat, ylab, ylim=None):
    n=pd.read_csv(D+'composition_vs_NLR.csv'); from scipy import stats
    for k,src in enumerate(['Real','PBISC']):
        d=n[n.source==src]; ax.scatter(d.NLR,d[feat],s=70,marker=SRC_MK[src],color=SRC_COL[src],edgecolors='white',linewidths=0.6,alpha=0.9,zorder=3)
        r,p=stats.spearmanr(d.NLR,d[feat])
        ax.text(0.03,0.97-0.12*k,f'{src} (n={len(d)}): ρ={r:+.2f}, p={p:.1g}',transform=ax.transAxes,va='top',fontsize=12.5,color=SRC_COL[src])
    ax.set_xscale('log'); ax.set_xticks([1,2,5,10]); ax.set_xticklabels(['1','2','5','10']); ax.minorticks_off(); ax.set_xlabel('Clinical neutrophil / lymphocyte ratio'); ax.set_ylabel(ylab)
    ax.set_ylim(*(ylim or (0, n[feat].max()*1.35)))

# ---------------- Figure 3 ----------------
DEG_SUB=[('CD14mono_IFN','CD14 mono IFN'),('pDC','pDC'),('NK CD56 bright','NK CD56bright'),('HSPC','HSPC'),('T g/d','T γδ'),('MAIT','MAIT'),('B_ASC_dividing','B ASC dividing'),('B_ASC_all','B ASC (all)'),
         ('REF_CD14mono','CD14 mono (ref.)'),('REF_NK','NK (ref.)'),('REF_T CD4 Helper','T CD4 helper (ref.)'),('REF_B_naive','B naive (ref.)')]
def deg_count_dots(ax):
    R=pd.read_csv(D+'rare_subset_severity_DEG_summary.csv').set_index('subset')
    cols=[('Real','Real'),('PBISC','PBISC'),('Combined','Combined')]
    for i,(s,lab) in enumerate(DEG_SUB):
        for j,(c,cl) in enumerate(cols):
            v=R.loc[s,f'{c}_trendDEG_q05'] if s in R.index else np.nan
            if pd.isna(v): ax.text(j,i,'–',ha='center',va='center',fontsize=13,color=PALETTE['gray_1']); continue
            ax.scatter(j,i,s=30+np.sqrt(v)*38,marker=SRC_MK[cl],color=SRC_COL[cl] if v>0 else 'white',edgecolors='black' if v>0 else PALETTE['neutral'],linewidths=0.8,zorder=3)
            ax.text(j+0.33,i,f'{int(v)}',va='center',ha='left',fontsize=12.5)
    ax.axhline(7.5,color=PALETTE['neutral'],ls='--',lw=1.2)
    ax.set_xticks(range(3)); ax.set_xticklabels([c for c,_ in cols],rotation=35,ha='left',rotation_mode='anchor'); ax.xaxis.tick_top(); ax.tick_params(axis='x',length=0); ax.spines['bottom'].set_visible(False)
    ax.set_yticks(range(len(DEG_SUB))); ax.set_yticklabels([l for _,l in DEG_SUB]); ax.set_ylim(len(DEG_SUB)-0.5,-0.6); ax.set_xlim(-0.5,2.95)
    ax.set_xlabel('Severity-trend DEGs (BH q < 0.05)',fontsize=13)

def cd14ifn_gene_scatter(ax, compact=False):
    G=pd.read_csv(D+'cd14monoIFN_gene_trend.csv').dropna(subset=['Real_trend_r','PBISC_trend_r'])
    sig=G.Combined_trend_q<0.05
    ax.scatter(G.Real_trend_r[~sig],G.PBISC_trend_r[~sig],s=10,color=PALETTE['neutral'],alpha=0.6,lw=0,zorder=1,rasterized=True)
    up=sig&(G.Combined_trend_r>0); dn=sig&(G.Combined_trend_r<0)
    ax.scatter(G.Real_trend_r[up],G.PBISC_trend_r[up],s=55,color=PALETTE['red_strong'],edgecolors='black',linewidths=0.6,zorder=3,label=f'Combined q<0.05, up ({up.sum()})')
    ax.scatter(G.Real_trend_r[dn],G.PBISC_trend_r[dn],s=55,color=PALETTE['blue_main'],edgecolors='black',linewidths=0.6,zorder=3,label=f'Combined q<0.05, down ({dn.sum()})')
    from adjustText import adjust_text
    lim=max(G.Real_trend_r.abs().max(),G.PBISC_trend_r.abs().max())*1.08; ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    txt=[]; xs=[]; ys=[]
    for g_ in (['S100A8','CD163','FKBP5','HLA-DRA'] if compact else ['S100A8','S100A12','CD163','FKBP5','MCEMP1','ADAMTS2','HLA-DRA','HLA-DPB1','HLA-DMB']):
        r=G[G.gene==g_]
        if len(r) and bool(sig[r.index[0]]):
            txt.append(ax.text(r.Real_trend_r.iloc[0],r.PBISC_trend_r.iloc[0],g_,fontsize=12,fontstyle='italic',zorder=5))
    adjust_text(txt,ax=ax,x=G.Real_trend_r[sig].values,y=G.PBISC_trend_r[sig].values,expand=(1.3,1.6),arrowprops=dict(arrowstyle='-',color='black',lw=0.8))
    ax.axhline(0,color=PALETTE['gray_1'],lw=1); ax.axvline(0,color=PALETTE['gray_1'],lw=1)
    k=G[['Real_trend_r','PBISC_trend_r']].corr().iloc[0,1]
    ax.set_xlabel('Real: r with peak severity'); ax.set_ylabel('PBISC: r with peak severity')
    ax.text(0.03,0.97,(f'r = {k:.2f}, {len(G):,} genes' if compact else f'all genes r = {k:.2f}\nReal 0 / PBISC 0 / Combined {int(sig.sum())} DEGs'),transform=ax.transAxes,va='top',fontsize=12 if compact else 12.5)
    if compact:
        for h_ in ax.collections[1:3]: h_.set_label(h_.get_label().replace('Combined q<0.05, ',''))
    ax.legend(loc='lower right',fontsize=11 if compact else 11.5,handletextpad=0.2,borderaxespad=0.1,title='Combined q < 0.05' if compact else None,title_fontsize=11)

PERM_SUB=[('CD14mono_IFN','CD14 mono IFN'),('HSPC','HSPC'),('B_ASC_all','B ASC (all)'),('NK CD56 bright','NK CD56bright'),('pDC','pDC'),('REF_CD14mono','CD14 mono (ref.)'),('REF_NK','NK (ref.)')]
def perm_null(ax):
    N=pd.read_csv(D+'combined_DEG_permutation_null.csv'); rng=np.random.default_rng(5)
    for i,(s,lab) in enumerate(PERM_SUB):
        d=N[N.subset==s]; v=d.n_DEG.values; obs=d.observed.iloc[0]
        ax.scatter(v+1+rng.uniform(-0.08,0.08,len(v))*(v+1),i+rng.uniform(-0.22,0.22,len(v)),s=12,color=PALETTE['neutral'],alpha=0.8,lw=0,zorder=2)
        ax.scatter(obs+1,i,s=170,marker='D',color=SRC_COL['Combined'],edgecolors='black',linewidths=1.0,zorder=4)
        fr=(v>=obs).mean(); ax.text(max(obs+1,2)*1.35,i,(f'{obs}  (p<0.005)' if fr==0 else f'{obs}  (p={fr:.2g})') if obs>0 else '0',va='center',fontsize=12)
    ax.set_xscale('log'); ax.set_xlim(0.7,4000); ax.set_xticks([1,2,11,101,1001]); ax.set_xticklabels(['0','1','10','100','1000'])
    ax.set_yticks(range(len(PERM_SUB))); ax.set_yticklabels([l for _,l in PERM_SUB]); ax.set_ylim(len(PERM_SUB)-0.5,-0.6)
    ax.set_xlabel('Combined severity DEGs (q < 0.05)')
    h=[Line2D([0],[0],marker='D',color='none',markerfacecolor=SRC_COL['Combined'],markeredgecolor='black',markersize=10,label='observed'),
       Line2D([0],[0],marker='o',color='none',markerfacecolor=PALETTE['neutral'],markersize=7,label='PBISC patients shuffled (200×)')]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.4,-0.17),fontsize=12,borderaxespad=0)

def gene_strip(ax, gene):
    Z=pd.read_csv(D+'cd14monoIFN_gene_patient_z.csv'); Z=Z[Z.gene==gene]; G=pd.read_csv(D+'cd14monoIFN_gene_trend.csv').set_index('gene').loc[gene]
    rng=np.random.default_rng(7)
    for src in ['Real','PBISC']:
        d=Z[Z.source==src]; off=-0.2 if src=='Real' else 0.2
        for i in range(4):
            v=d.loc[d.severity_ord==i,'z'].values
            ax.scatter(i+off+rng.uniform(-0.07,0.07,len(v)),v,s=38,marker=SRC_MK[src],color=SRC_COL[src],alpha=0.85,edgecolors='white',linewidths=0.5,zorder=3)
            if len(v): ax.plot([i+off-0.13,i+off+0.13],[np.median(v)]*2,color='black',lw=2.4,zorder=4)
    ax.axhline(0,color=PALETTE['neutral'],lw=1)
    ax.set_xticks(range(4)); ax.set_xticklabels(['Mild','Mod.','Sev.','Crit.']); ax.set_xlim(-0.55,3.55)
    ax.set_title(gene,loc='left',fontsize=14,fontstyle='italic',pad=4); ax.set_ylabel('z (log CPM)')
    ax.text(0.98,0.97,f"Comb. q={fmt_q(G.Combined_trend_q)}",transform=ax.transAxes,ha='right',va='top',fontsize=12,color=SRC_COL['Combined'])
    yl=np.nanmax(np.abs(Z.z))*1.25; ax.set_ylim(-yl,yl*1.15)

def broadcast_dots(ax):
    R=pd.read_csv(D+'rare_subset_severity_DEG_summary.csv').set_index('subset'); subs=DEG_SUB
    sty=[('logFC_r_real_vs_pbisc','Real vs PBISC (same subset)',PALETTE['red_strong'],'s'),('PBISC_logFC_r_with_own_whole_sample','PBISC subset vs PBISC whole sample',SRC_COL['PBISC'],'o'),('Real_logFC_r_with_own_whole_sample','Real subset vs Real whole sample',SRC_COL['Real'],'D')]
    for i,(s,lab) in enumerate(subs):
        if s not in R.index: continue
        v=[R.loc[s,c] for c,_,_,_ in sty]; vv=[x for x in v if not pd.isna(x)]
        if vv: ax.plot([min(vv),max(vv)],[i,i],color=PALETTE['neutral'],lw=2.5,zorder=1)
        for (c,l,col,mk),x in zip(sty,v):
            if not pd.isna(x): ax.scatter(x,i,s=80,marker=mk,color=col,edgecolors='black',linewidths=0.7,zorder=3,label=l if i==0 else None)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.axhline(7.5,color=PALETTE['neutral'],ls='--',lw=1.2)
    ax.set_yticks(range(len(subs))); ax.set_yticklabels([l for _,l in subs]); ax.set_ylim(len(subs)-0.5,-0.6); ax.set_xlim(-0.25,1.0)
    ax.set_xlabel('Pearson r of severity log$_2$FC'); ax.legend(loc='upper center',bbox_to_anchor=(0.4,-0.17),fontsize=12,ncol=1,borderaxespad=0)

# ---------------- Figure 4 ----------------
SUBCL=[('CD14mono_IFN','CD14 mono IFN'),('NK CD56 bright','NK CD56bright'),('B_ASC_all','B ASC (all)'),('pDC','pDC'),('T g/d','T γδ'),('MAIT','MAIT'),('HSPC','HSPC')]
def ari_dots(ax):
    R=pd.read_csv(D+'rare_subclustering_summary.csv'); off={'Real':-0.22,'PBISC':0,'Joint':0.22}; col={'Real':SRC_COL['Real'],'PBISC':SRC_COL['PBISC'],'Joint':SRC_COL['Combined']}; mk={'Real':'D','PBISC':'o','Joint':'s'}
    for i,(s,lab) in enumerate(SUBCL):
        for src in ['Real','PBISC','Joint']:
            v=R[(R.subset==s)&(R.source==src)].ARI_split_HVG.iloc[0]
            ax.scatter(v,i+off[src],s=95,marker=mk[src],color=col[src],edgecolors='black',linewidths=0.8,zorder=3,label=src if i==0 else None)
    ax.axvspan(0.5,1.0,color=PALETTE['green_1'],alpha=0.7,lw=0,zorder=0); ax.text(0.52,len(SUBCL)-0.75,'ARI ≥ 0.5',fontsize=12,color='#2E7D32',va='bottom')
    ax.set_xlim(0,0.75); ax.set_yticks(range(len(SUBCL))); ax.set_yticklabels([l for _,l in SUBCL]); ax.set_ylim(len(SUBCL)-0.5,-0.6)
    ax.set_xlabel('ARI between two random halves of HVGs'); ax.legend(loc='upper right',fontsize=12,title='clustered on',title_fontsize=12)
    r=R[(R.subset=='B_ASC_all')&(R.source=='Real')].ARI_split_HVG.iloc[0]

def nk_umap(ax, mode):
    U=pd.read_csv(D+'nk_cd56bright_joint_umap.csv.gz',index_col=0); rng=np.random.default_rng(0); U=U.iloc[rng.permutation(len(U))]
    if mode=='source':
        for src in ['Real','PBISC']:
            d=U[U.source==src]; ax.scatter(d.UMAP1,d.UMAP2,s=5,color=SRC_COL[src],alpha=0.55,lw=0,rasterized=True,label=f'{src} ({len(d):,})')
        ax.legend(loc='upper left',fontsize=12,markerscale=3.5,handletextpad=0.2,borderaxespad=0.1)
    else:
        m=U.leiden.astype(str)=='17'
        ax.scatter(U.UMAP1[~m],U.UMAP2[~m],s=5,color=PALETTE['neutral'],alpha=0.5,lw=0,rasterized=True,label='other clusters')
        ax.scatter(U.UMAP1[m],U.UMAP2[m],s=8,color=PALETTE['red_strong'],alpha=0.9,lw=0,rasterized=True,label=f'cluster 17 ({m.sum():,})')
        ax.legend(loc='upper left',fontsize=12,markerscale=3,handletextpad=0.2,borderaxespad=0.1)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    for sp in ('left','bottom'): ax.spines[sp].set_linewidth(1.5)

def c17_strip(ax):
    F=pd.read_csv(D+'nk_cluster17_fraction_patient.csv'); J=pd.read_csv(D+'rare_subclustering_joint_cluster_severity.csv')
    j=J[(J.subset=='NK CD56 bright')&(J.cluster.astype(str)=='17')].iloc[0]; rng=np.random.default_rng(11); from scipy import stats
    for k,src in enumerate(['Real','PBISC']):
        d=F[F.source==src]; off=-0.2 if src=='Real' else 0.2
        for i in range(4):
            v=d.loc[d.severity_ord==i,'frac17_pct'].values
            ax.scatter(i+off+rng.uniform(-0.07,0.07,len(v)),v,s=44,marker=SRC_MK[src],color=SRC_COL[src],alpha=0.85,edgecolors='white',linewidths=0.6,zorder=3)
            if len(v): ax.plot([i+off-0.13,i+off+0.13],[np.median(v)]*2,color='black',lw=2.6,zorder=4)
        r,p=stats.spearmanr(d.frac17_pct,d.severity_ord); ax.text(0.02,0.97-0.11*k,f'{src} (n={len(d)}): ρ={r:+.2f}, p={p:.2g}',transform=ax.transAxes,va='top',fontsize=12.5,color=SRC_COL[src])
    ax.text(0.02,0.75,f'Combined (n={int(j.n)}): ρ={j.rho:+.2f}, q={j.q:.3f}',transform=ax.transAxes,va='top',fontsize=12.5,color=SRC_COL['Combined'])
    ax.set_xticks(range(4)); ax.set_xticklabels(['Mild','Mod.','Severe','Crit.']); ax.set_xlim(-0.55,3.55)
    ax.set_ylabel('Cluster 17 (% of NK CD56bright)'); ax.set_ylim(0,F.frac17_pct.max()*1.55)

MKG=[('cluster 17 top',['NKG7','GZMA','CCL4','CST7']),('CD56dim-type',['GZMB','PRF1','SPON2','FCGR3A','CX3CR1','FGFBP2']),('CD56bright-type',['GZMK','XCL1','KLRC1','NCAM1','SELL','IL7R','TCF7'])]
def marker_dots(ax, cax=None):
    M=pd.read_csv(D+'nk_cluster17_markers.csv'); genes=[g for _,gs in MKG for g in gs]
    rows=[('Real','cluster 17'),('Real','other clusters'),('PBISC','cluster 17'),('PBISC','other clusters')]
    piv=M.pivot_table(index=['source','group'],columns='gene',values='mean_lognorm')[genes]; fr=M.pivot_table(index=['source','group'],columns='gene',values='frac_expr')[genes]
    zc=(piv-piv.mean())/piv.std(ddof=0)
    for i,r in enumerate(rows):
        sc=ax.scatter(np.arange(len(genes)),np.full(len(genes),i),s=20+fr.loc[r].values*260,c=zc.loc[r].values,cmap=DIV,vmin=-1.6,vmax=1.6,edgecolors='black',linewidths=0.6,zorder=3)
    x0=0
    for name,gs in MKG:
        ax.plot([x0-0.35,x0+len(gs)-0.65],[-0.95,-0.95],color='black',lw=1.6,clip_on=False); ax.text(x0+(len(gs)-1)/2,-1.1,name,ha='center',va='bottom',fontsize=12.5)
        x0+=len(gs)
        if x0<len(genes): ax.axvline(x0-0.5,color=PALETTE['neutral'],lw=1.2,ls='--')
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes,rotation=55,ha='right',fontstyle='italic',fontsize=13)
    ax.set_yticks(range(4)); ax.set_yticklabels([f'{s}: {g.replace(" clusters","")}' for s,g in rows]); ax.set_ylim(3.5,-0.6); ax.set_xlim(-0.6,len(genes)-0.4)
    ax.axhline(1.5,color=PALETTE['neutral'],lw=1.2)
    if cax is not None:
        cb=plt.colorbar(sc,cax=cax); cb.set_label('mean expr. (z across rows)',fontsize=12); cb.ax.tick_params(labelsize=11); cb.outline.set_linewidth(1)
        cax.text(0.5,-0.18,'size = % cells\nexpressing',transform=cax.transAxes,ha='center',va='top',fontsize=11.5)

# ---------------- Supplementary (bars) ----------------
def comp_stack(ax):
    P=_pos(); cts=[('T_pct','T cell'),('B_pct','B cell'),('NK_pct','NK cell'),('Mono_pct','Monocyte'),('Other_pct','Other')]
    xs=[]; labs=[]; x=0
    for src in ['Real','PBISC']:
        for s in SEV:
            d=P[(P.source==src)&(P.WHO_severity==s)]; bottom=0
            for c,lab in cts:
                v=d[c].mean(); ax.bar(x,v,bottom=bottom,width=0.78,color=CT_COL[lab],edgecolor='black',linewidth=1.2,label=lab if x==0 else None); bottom+=v
            ax.text(x,101,f'n={len(d)}',ha='center',va='bottom',fontsize=11.5)
            xs.append(x); labs.append({'mild':'Mild','moderate':'Mod.','severe':'Sev.','critical':'Crit.'}[s]); x+=1
        x+=0.6
    ax.set_xticks(xs); ax.set_xticklabels(labs,rotation=0,fontsize=12.5); ax.set_ylim(0,108); ax.set_ylabel('% of cells (mean of patients)')
    ax.text(1.5,-15,'Real',ha='center',fontsize=14,fontweight='bold',color=SRC_COL['Real']); ax.text(6.1,-15,'PBISC',ha='center',fontsize=14,fontweight='bold',color=SRC_COL['PBISC'])
    ax.legend(loc='upper left',bbox_to_anchor=(1.0,1.0),fontsize=12.5)

def rare_scale_bars(ax):
    P=_pos(); feats=[f for f,_ in FEAT[7:] if f!='rare_ILC_pct']+['rare_ILC_pct']
    y=np.arange(len(feats)); w=0.38
    for k,src in enumerate(['Real','PBISC']):
        v=[P[P.source==src][f].mean() for f in feats]
        ax.barh(y+(-w/2 if k==0 else w/2),v,height=w,color=SRC_COL[src],edgecolor='black',linewidth=1.0,label=src,hatch='' if k else '//')
    ax.set_yticks(y); ax.set_yticklabels([FNAME[f] for f in feats]); ax.set_ylim(len(feats)-0.4,-0.6); ax.set_xscale('symlog',linthresh=0.05); ax.set_xlim(0,5)
    ax.set_xticks([0,0.1,0.5,1,5]); ax.set_xticklabels(['0','0.1','0.5','1','5'])
    ax.set_xlabel('% of cells (mean of COVID+ patients)'); ax.legend(loc='lower right',fontsize=12.5)

# ---------------- Supplementary C2: is cluster 17 CD56dim contamination? ----------------
POPC={'CD56bright (other)':PALETTE['gray_1'],'cluster 17':PALETTE['red_strong'],'CD56dim (NK)':PALETTE['violet']}
POPL={'CD56bright (other)':'CD56bright, other','cluster 17':'Cluster 17','CD56dim (NK)':'CD56dim NK'}
def dim_marker_dots(ax, cax=None):
    M=pd.read_csv(D+'nk_c17_vs_bright_dim_markers.csv'); genes=['FCGR3A','GZMB','PRF1','KLRC1','SELL','NCAM1']
    cols=[(s,p) for s in ['Real','PBISC'] for p in ['bright','c17','dim']]
    lab={'bright':'bright','c17':'cl. 17','dim':'dim'}
    for i,g in enumerate(genes):
        r=M[M.gene==g].set_index('source'); v=np.array([r.loc[s,f'mean_{p}'] for s,p in cols]); f=np.array([r.loc[s,f'pct_{p}'] for s,p in cols])
        zz=(v-v.mean())/(v.std()+1e-9)
        sc=ax.scatter(np.arange(6),np.full(6,i),s=20+f*2.6,c=zz,cmap=DIV,vmin=-1.6,vmax=1.6,edgecolors='black',linewidths=0.6,zorder=3)
    ax.set_xticks(range(6)); ax.set_xticklabels([lab[p] for _,p in cols],fontsize=12.5,rotation=40,ha='right',rotation_mode='anchor'); ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes,fontstyle='italic')
    ax.set_ylim(len(genes)-0.5,-1.4); ax.set_xlim(-0.6,5.6); ax.axvline(2.5,color=PALETTE['neutral'],lw=1.2)
    for x,s in [(1,'Real'),(4,'PBISC')]: ax.text(x,-1.05,s,ha='center',va='center',fontsize=14,fontweight='bold',color=SRC_COL[s])
    ax.text(-0.95,4.0,'CD56 RNA',ha='right',va='center',fontsize=11,color=PALETTE['gray_1'],transform=ax.transData,alpha=0)
    if cax is not None:
        cb=plt.colorbar(sc,cax=cax); cb.set_label('mean expr. (z across columns)',fontsize=12); cb.ax.tick_params(labelsize=11); cb.outline.set_linewidth(1)
        cax.text(0.5,1.04,'size = %\nexpressing',transform=cax.transAxes,ha='center',va='bottom',fontsize=11.5)

def score_scatter(ax, src, legend=False):
    S=pd.read_csv(D+'nk_c17_cell_scores.csv.gz'); S=S[S.source==src]; rng=np.random.default_rng(0)
    for p,sz,al in [('CD56bright (other)',5,0.25),('CD56dim (NK)',5,0.25),('cluster 17',16,0.9)]:
        d=S[S['pop']==p]
        if len(d)>3000: d=d.iloc[rng.choice(len(d),3000,replace=False)]
        ax.scatter(d.bright_score,d.dim_score,s=sz,color=POPC[p],alpha=al,lw=0,rasterized=True,label=POPL[p],zorder=3 if p=='cluster 17' else 2)
    for p in POPC:
        d=S[S['pop']==p]; ax.scatter(d.bright_score.median(),d.dim_score.median(),s=170,marker='X',color=POPC[p],edgecolors='black',linewidths=1.2,zorder=5)
    ax.set_xlabel('CD56bright lineage score'); ax.set_ylabel('CD56dim lineage score'); ax.set_title(src,loc='left',fontsize=14,color=SRC_COL[src],fontweight='bold',pad=4)
    if legend: ax.legend(loc='upper right',fontsize=11.5,markerscale=2.2,handletextpad=0.2,borderaxespad=0.1)
    ax.set_xlim(-1.2,2.6); ax.set_ylim(-1.0,2.2)

def _ecdf(ax,x,**kw):
    x=np.sort(x); ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',**kw)
def lineage_ecdf(ax):
    S=pd.read_csv(D+'nk_c17_cell_scores.csv.gz'); X=pd.read_csv(D+'nk_c17_mixture_test.csv').set_index('source')
    for src,ls in [('Real','-'),('PBISC','--')]:
        for p in POPC: _ecdf(ax,S[(S.source==src)&(S['pop']==p)].lineage_axis.values,color=POPC[p],lw=2.6,ls=ls)
    ax.set_xlabel('Lineage axis (dim − bright score)'); ax.set_ylabel('Cumulative fraction of cells'); ax.set_xlim(-3,3); ax.set_ylim(0,1.02)
    h=[Line2D([0],[0],color=POPC[p],lw=2.6,label=POPL[p]) for p in POPC]+[Line2D([0],[0],color='black',lw=2,ls='-',label='Real'),Line2D([0],[0],color='black',lw=2,ls='--',label='PBISC')]
    ax.legend(handles=h,loc='lower right',fontsize=11.5,handlelength=1.6,borderaxespad=0.1)
    t=f"cl. 17 position (0=bright, 1=dim):\nReal {X.loc['Real','c17_position_median']:.2f}, PBISC {X.loc['PBISC','c17_position_median']:.2f}\nbright+dim mixture: KS D {X.loc['Real','KS_c17_vs_mixture_D']:.2f} / {X.loc['PBISC','KS_c17_vs_mixture_D']:.2f}"
    ax.text(0.0,1.03,t,transform=ax.transAxes,va='bottom',fontsize=11.5)

def pdim_ecdf(ax):
    C=pd.read_csv(D+'nk_c17_classifier_pdim.csv.gz')
    for src,ls in [('Real','-'),('PBISC','--')]:
        for p in POPC: _ecdf(ax,C[(C.source==src)&(C['pop']==p)].p_dim.values,color=POPC[p],lw=2.6,ls=ls)
    ax.axvline(0.5,color=PALETTE['gray_1'],ls=':',lw=1.6)
    ax.set_xlabel('Classifier P(CD56dim)'); ax.set_ylabel('Cumulative fraction of cells'); ax.set_xlim(-0.02,1.02); ax.set_ylim(0,1.02)
    k={s:100*(C[(C.source==s)&(C['pop']=='cluster 17')].p_dim>0.5).mean() for s in ['Real','PBISC']}
    b={s:100*(C[(C.source==s)&(C['pop']=='CD56bright (other)')].p_dim>0.5).mean() for s in ['Real','PBISC']}
    ax.text(0.0,1.03,f"cells with P > 0.5:\ncluster 17: Real {k['Real']:.0f}%, PBISC {k['PBISC']:.0f}%\nbright other: Real {b['Real']:.0f}%, PBISC {b['PBISC']:.0f}%",transform=ax.transAxes,ha='left',va='bottom',fontsize=11.5)

def c17_robust(ax):
    R=pd.read_csv(D+'nk_c17_robustness.csv'); feats=['cluster 17 (all)','cluster 17, bright-like (P(dim)≤0.5)','cluster 17, dim-like (P(dim)>0.5)']
    lab=['Cluster 17, all','bright-like\n(P(dim) ≤ 0.5)','dim-like\n(P(dim) > 0.5)']; off={'Real':-0.22,'PBISC':0,'Combined':0.22}
    for i,f in enumerate(feats):
        for a in ['Real','PBISC','Combined']:
            r=R[(R.feature==f)&(R.analysis==a)].iloc[0]
            ax.scatter(r.rho,i+off[a],s=110,marker=SRC_MK[a],color=SRC_COL[a],edgecolors='black' if r.p<0.05 else PALETTE['neutral'],linewidths=1.6 if r.p<0.05 else 0.8,zorder=3,label=a if i==0 else None)
            ax.text(0.68,i+off[a],f'p={r.p:.2g}',va='center',fontsize=11.5,color=SRC_COL[a])
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.set_xlim(-0.1,0.95); ax.set_xticks([0,0.25,0.5])
    ax.set_yticks(range(3)); ax.set_yticklabels(lab); ax.set_ylim(2.5,-0.5); ax.set_xlabel('ρ with peak severity')
    ax.legend(loc='upper center',fontsize=11.5,borderaxespad=0,bbox_to_anchor=(0.35,-0.2),ncol=3,columnspacing=0.8,handletextpad=0.2)

# ---------------- Figure 7 additions ----------------
SEVC={0:PALETTE['green_2'],1:PALETTE['red_1'],2:PALETTE['red_2'],3:PALETTE['red_strong']}
def nk_umap_severity(ax):
    U=pd.read_csv(D+'nk_cd56bright_joint_umap.csv.gz',index_col=0); rng=np.random.default_rng(0); U=U.iloc[rng.permutation(len(U))]
    d=U[U.severity_ord.isna()]; ax.scatter(d.UMAP1,d.UMAP2,s=4,color=PALETTE['neutral'],alpha=0.4,lw=0,rasterized=True,label='COVID− / recovery')
    for k,l in zip(range(4),['Mild','Moderate','Severe','Critical']):
        d=U[U.severity_ord==k]; ax.scatter(d.UMAP1,d.UMAP2,s=5,color=SEVC[k],alpha=0.7,lw=0,rasterized=True,label=l)
    m=U.leiden.astype(str)=='17'; from scipy.spatial import ConvexHull
    pts=U.loc[m,['UMAP1','UMAP2']].values; c=np.median(pts,0); keep=np.linalg.norm(pts-c,axis=1)<np.quantile(np.linalg.norm(pts-c,axis=1),0.9); pts=pts[keep]
    h=ConvexHull(pts); poly=np.r_[pts[h.vertices],pts[h.vertices][:1]]; ax.plot(poly[:,0],poly[:,1],color='black',lw=1.6,ls='--'); ax.text(poly[:,0].max(),poly[:,1].max(),' cl. 17',fontsize=12,va='bottom')
    ax.legend(loc='lower left',fontsize=11,markerscale=3,handletextpad=0.1,borderaxespad=0.05,ncol=1,labelspacing=0.2)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    for sp in ('left','bottom'): ax.spines[sp].set_linewidth(1.5)

def nk_gene_strip(ax, gene):
    P=pd.read_csv(D+'nkbright_gzmb_prf1_patient.csv'); P=P[P.gene==gene]; T=pd.read_csv(D+'nkbright_gzmb_prf1_stats.csv'); T=T[T.gene==gene].set_index('analysis')
    rng=np.random.default_rng(13)
    for src in ['Real','PBISC']:
        d=P[P.source==src]; off=-0.2 if src=='Real' else 0.2
        for i in range(4):
            v=d.loc[d.severity_ord==i,'z'].values
            ax.scatter(i+off+rng.uniform(-0.07,0.07,len(v)),v,s=40,marker=SRC_MK[src],color=SRC_COL[src],alpha=0.85,edgecolors='white',linewidths=0.5,zorder=3)
            if len(v): ax.plot([i+off-0.13,i+off+0.13],[np.median(v)]*2,color='black',lw=2.4,zorder=4)
    ax.axhline(0,color=PALETTE['neutral'],lw=1)
    ax.set_xticks(range(4)); ax.set_xticklabels(['Mild','Mod.','Severe','Crit.']); ax.set_xlim(-0.55,3.55)
    ax.set_ylabel(f'{gene} in NK CD56bright\n(z, log CPM)'); yl=np.nanmax(np.abs(P.z))*1.1; ax.set_ylim(-yl,yl*1.55)
    for k,a in enumerate(['Real','PBISC','Combined']):
        r=T.loc[a]; ax.text(0.02+0.33*k,0.97,f'{a}\nr={r.trend_r:+.2f}, p={r.trend_p:.2g}',transform=ax.transAxes,va='top',fontsize=12,color=SRC_COL[a])

def cluster_rho_scatter(ax):
    W=pd.read_csv(D+'nk_joint_cluster_rho_by_source.csv'); W['cluster']=W.cluster.astype(str)
    sz=20+np.sqrt(W.n_real+W.n_pbisc)*6
    sig=W.q_Real<0.05
    ax.axhline(0,color=PALETTE['gray_1'],lw=1); ax.axvline(0,color=PALETTE['gray_1'],lw=1); ax.plot([-0.8,0.8],[-0.8,0.8],color=PALETTE['neutral'],ls='--',lw=1.3)
    ax.scatter(W.rho_Real,W.rho_PBISC,s=sz,color=np.where(sig,PALETTE['gray_3'],PALETTE['neutral']),edgecolors='white',linewidths=0.6,alpha=0.9,zorder=3)
    for c,col in [('17',PALETTE['red_strong']),('21',PALETTE['violet'])]:
        r=W[W.cluster==c].iloc[0]; ax.scatter(r.rho_Real,r.rho_PBISC,s=sz[W.cluster==c].iloc[0],color=col,edgecolors='black',linewidths=1.0,zorder=4)
        ax.annotate(f'cl. {c}\nReal q={r.q_Real:.2g}, PBISC q={r.q_PBISC:.2g}',(r.rho_Real,r.rho_PBISC),xytext=(-10,14) if c=='17' else (10,-2),ha='right' if c=='17' else 'left',textcoords='offset points',fontsize=11.5,color=col)
    ax.set_xlim(-0.85,0.85); ax.set_ylim(-0.6,0.6); ax.set_xlabel('Real: ρ with peak severity'); ax.set_ylabel('PBISC: ρ with peak severity')
    h=[Line2D([0],[0],marker='o',color='none',markerfacecolor=PALETTE['gray_3'],markersize=9,label='Real BH q < 0.05'),Line2D([0],[0],marker='o',color='none',markerfacecolor=PALETTE['neutral'],markersize=9,label='Real q ≥ 0.05')]
    ax.legend(handles=h,loc='lower right',fontsize=11.5,borderaxespad=0.1); ax.text(0.02,0.97,'52 joint NK CD56bright clusters\n(size = cells)',transform=ax.transAxes,va='top',fontsize=11.5)
