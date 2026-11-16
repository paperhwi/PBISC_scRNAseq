"""Ver3 panels (2026-11-15): cohort-extension and PPI tests, comparison with deconvolution and a second adapter,
GZMB/PRF1 control designs and the bulk input, robustness. Same figures4papers style as panels.py / panels_v2.py."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import TwoSlopeNorm
from panels import PALETTE, SRC_COL, SRC_MK, DIV, SEVC, fmt_q, FNAME
D3='data_v3/'
COL=dict(Real=SRC_COL['Real'],PBISC=SRC_COL['PBISC'],PBISC_only=PALETTE['blue_secondary'],Extension=PALETTE['teal'],
         Deconv=PALETTE['red_strong'],Draw0='#9DBBDD',Bulk=PALETTE['violet'],Shuffled=PALETTE['neutral'])
MK=dict(Real='D',PBISC='o',PBISC_only='o',Extension='s',Deconv='^',Draw0='o',Bulk='v',Shuffled='s')
KEY=[('rare_NK CD56 bright_pct','NK CD56bright'),('rare_B_ASC_dividing_pct','B ASC dividing'),('DC_within_Other_pct','DC'),('rare_DC1_pct','DC1'),
     ('rare_HSPC_pct','HSPC'),('NK_pct','NK cell'),('T_pct','T cell')]
def _p(p): return 'n.a.' if pd.isna(p) else (f'{p:.0e}'.replace('e-0','e-') if p<0.001 else f'{p:.2g}')
def _h(c,m,lab,fc=None,ms=9): return Line2D([0],[0],marker=m,color='none',markerfacecolor=fc or c,markeredgecolor='black' if fc is None else c,markeredgewidth=0.8 if fc is None else 1.6,markersize=ms,label=lab)

def schematic_v3(ax):
    ax.set_xlim(0,100); ax.set_ylim(0,20); ax.axis('off')
    def box(x,y,w,h,txt,fc,ec,fs=12):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.3,rounding_size=1.2',fc=fc,ec=ec,lw=2)); ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs,linespacing=1.25)
    arr=dict(arrowstyle='-|>',lw=2,color=PALETTE['gray_2'])
    box(2.5,11,21.5,8,'Real scRNA-seq: 39 COVID+\n24 independent controls\n14 pre-infection baselines','#EDEDED',SRC_COL['Real'],fs=11.5)
    box(2.5,1,21.5,8,'PBISC from bulk: 57 COVID+\n(28 also in Real, 29 new)\n19 controls; same bulk','#E3ECF6',SRC_COL['PBISC'],fs=11.5)
    box(29,11,20,8,'Per source\nReal | PBISC |\nPBISC-only (29 new)','white',PALETTE['gray_2'],fs=11.5)
    box(29,1,20,8,'Extension: Real 39 +\n29 new PBISC patients;\nPPI on 28 shared','#E2F1F2',COL['Extension'],fs=11.5)
    for y in (15,5): ax.annotate('',xy=(28.6,y),xytext=(24.4,y),arrowprops=arr)
    box(54,11,21,8,'Baselines: deconvolution\nof the same bulk;\nsecond adapter (draw 0)','#FBEDEC',COL['Deconv'],fs=11.5)
    box(54,1,21,8,'Controls: contemporaneous\n(NCL 2020) and\nwithin-patient (Cambridge)','white',PALETTE['gray_2'],fs=11.5)
    for y in (15,5): ax.annotate('',xy=(53.6,y),xytext=(49.4,y),arrowprops=arr)
    box(80,11,19,8,'Cell-type fractions\nvs peak severity','white',PALETTE['gray_2'],fs=11.5)
    box(80,1,19,8,'NK GZMB / PRF1\nfrom cells to bulk','white',PALETTE['gray_2'],fs=11.5)
    for y in (15,5): ax.annotate('',xy=(79.6,y),xytext=(75.4,y),arrowprops=arr)

def dotgrid_v3(fig, place, w, h):
    A=pd.read_csv(D3+'composition_association_v3.csv').set_index('feature'); feats=[f for f in FNAME]
    ax=place(2.7,2.0,w-2.7-0.2,h-2.0-1.6); cax=place(2.6,1.15,(w-2.8)*0.8,0.2); sax=place(2.0,0.0,w-2.05,0.5)
    cols=[('Real','Real'),('PBISC','PBISC'),('PBISC_only','PBISC-only'),('Extension','Extension')]
    for j,(k,lab) in enumerate(cols):
        d=A.reindex(feats); q=d[f'{k}_q'].values; r=d[f'{k}_rho'].values; s=np.clip(-np.log10(np.where(np.isnan(q),1,q)),0,5)
        sc=ax.scatter(np.full(len(feats),j),np.arange(len(feats)),s=40+s*95,c=np.nan_to_num(r),cmap=DIV,norm=TwoSlopeNorm(0,-0.8,0.8),
                      edgecolors=np.where(q<0.05,'black','#BBBBBB'),linewidths=np.where(q<0.05,2.0,0.8),zorder=3)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([f'{lab} (n={int(A[f"{k}_n"].max())})' for k,lab in cols],fontsize=13,rotation=35,ha='left',rotation_mode='anchor')
    rec=A.reindex(feats).recoverable.values
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([FNAME[f]+('' if ok else ' †') for f,ok in zip(feats,rec)])
    for t,ok in zip(ax.get_yticklabels(),rec): t.set_color('black' if ok else PALETTE['gray_1'])
    ax.set_ylim(len(feats)-0.5,-0.6); ax.set_xlim(-0.6,len(cols)-0.4); ax.axhline(6.5,color=PALETTE['neutral'],lw=1.2,ls='--',zorder=1)
    ax.spines['bottom'].set_visible(False); ax.tick_params(axis='x',length=0); ax.xaxis.tick_top()
    cb=plt.colorbar(sc,cax=cax,orientation='horizontal'); cb.set_label('Spearman ρ with peak severity',fontsize=13); cb.ax.tick_params(labelsize=12); cb.outline.set_linewidth(1)
    sax.axis('off'); sax.text(-0.6,0,'BH q',fontsize=13,ha='right',va='center')
    for x,(v,l) in zip([0,1.75,3.6],[(1,'0.1'),(2,'0.01'),(4,'1e-4')]):
        sax.scatter(x,0,s=40+v*95,color='white',edgecolors='black',linewidths=1.0); sax.text(x+0.38,0,l,ha='left',va='center',fontsize=12)
    sax.scatter(5.6,0,s=40+2*95,color='white',edgecolors='black',linewidths=2.0); sax.text(5.95,0,'q < 0.05',ha='left',va='center',fontsize=12); sax.set_xlim(-0.9,7.6); sax.set_ylim(-1,1)

def source_compare(ax):
    """Severity rho of key features: Real vs PBISC (adapter draw 1, draw 0) vs deconvolution of the same bulk."""
    A=pd.read_csv(D3+'composition_association_v3.csv').set_index('feature'); R=pd.read_csv(D3+'robustness_v3.csv'); DC=pd.read_csv(D3+'deconvolution_vs_pbisc_severity_v3.csv')
    meth=[('Real','Real scRNA (39)',lambda f:(A.loc[f,'Real_rho'],A.loc[f,'Real_p'])),
          ('PBISC','PBISC, adapter 1 (57)',lambda f:(A.loc[f,'PBISC_rho'],A.loc[f,'PBISC_p'])),
          ('Draw0','PBISC, adapter 0 (57)',lambda f:tuple(R[(R.analysis.str.contains('draw 0, all cells'))&(R.feature==f)][['rho','p']].iloc[0])),
          ('Deconv','Deconvolution, same bulk (57)',lambda f:tuple(DC[(DC.method=='deconvolution (calibrated)')&(DC.feature==f)][['rho','p']].iloc[0]))]
    off=np.linspace(-0.27,0.27,len(meth))
    for i,(f,lab) in enumerate(KEY):
        vals=[m[2](f) for m in meth]; xs=[v[0] for v in vals]; ax.plot([min(xs),max(xs)],[i,i],color=PALETTE['neutral'],lw=2,zorder=1)
        for k,((key,_,_),(r,p)) in enumerate(zip(meth,vals)):
            sig=p<0.05; ax.scatter(r,i+off[k],s=80,marker=MK[key],color=COL[key] if sig else 'white',edgecolors='black' if sig else COL[key],linewidths=0.8 if sig else 1.6,zorder=3)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.set_xlim(-0.75,0.75)
    ax.set_yticks(range(len(KEY))); ax.set_yticklabels([l for _,l in KEY]); ax.set_ylim(len(KEY)-0.5,-0.6)
    ax.set_xlabel('Spearman ρ with peak severity (COVID+ patients)')
    h=[_h(COL[k],MK[k],lab) for k,lab,_ in meth]+[_h(PALETTE['gray_2'],'o','open: P ≥ 0.05',fc='white')]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.42,-0.2),fontsize=11,borderaxespad=0,handletextpad=0.2,ncol=3,columnspacing=0.8)

POWF=[('rare_NK CD56 bright_pct','NK CD56bright'),('rare_B_ASC_dividing_pct','B ASC dividing'),('DC_within_Other_pct','DC'),('rare_HSPC_pct','HSPC (Real null)')]
def ext_power(fig, place, w, h, feats=POWF, src='composition', compact=False):
    P=pd.read_csv(D3+('extension_power_v3.csv' if src=='composition' else 'nk_clusters_extension_power_v3.csv'))
    n=len(feats); gap=0.35; aw=(w-1.35-0.15-gap*(n-1))/n
    for k,(f,lab) in enumerate(feats):
        yb=2.45 if compact else 1.9; ax=place(1.35+k*(aw+gap),yb,aw,h-yb-0.55)
        d=P[P.feature==f] if src=='composition' else P[P.cluster.astype(str)==f]
        for col,key,ls in [('real_alone','Real','-'),('plus_pbisc_only','Extension','-'),('plus_shuffled_pbisc','Shuffled','--')]:
            ax.plot(d.n_real,d[col],color=COL[key],lw=2.3,ls=ls,marker=MK[key],ms=6,mec='black',mew=0.5)
        ax.set_ylim(-0.03,1.05); ax.set_xlim(7,38); ax.set_xticks([10,20,30]); ax.axhline(0.8,color=PALETTE['gray_1'],ls=':',lw=1.2)
        ax.set_title(lab,fontsize=13,loc='left',pad=4)
        if k: ax.set_yticklabels([])
        else: ax.set_ylabel('Detection rate\n(P < 0.05)')
        ax.tick_params(labelsize=12)
    fig.text(0,0,'',fontsize=1)
    lax=place(0.2 if compact else 1.35,0.0 if compact else 0.45,w-0.3 if compact else w-1.5,1.05 if compact else 0.55); lax.axis('off')
    h_=[Line2D([0],[0],color=COL['Real'],marker='D',lw=2.3,mec='black',label='Real alone'),Line2D([0],[0],color=COL['Extension'],marker='s',lw=2.3,mec='black',label='+ 29 new PBISC patients'),
        Line2D([0],[0],color=COL['Shuffled'],marker='s',lw=2.3,ls='--',mec='black',label='+ same, severity shuffled')]
    lax.legend(handles=h_,loc='center',ncol=1 if compact else 3,fontsize=11.5 if compact else 12,frameon=False,handlelength=2.2,columnspacing=1.0,labelspacing=0.25)
    xl=place(1.35,1.6 if compact else 1.1,w-1.5,0.4); xl.axis('off'); xl.text(0.5,0.5,'Real COVID+ patients sampled',ha='center',va='center',fontsize=13.5,transform=xl.transAxes)

def ppi_panel(ax, feats=POWF):
    P=pd.read_csv(D3+'ppi_v3.csv'); est=[('classical_real39','Real, all 39','Real'),('classical_shared28','Real, 28 shared','Real'),('PPI++_28+29','PPI++ (28 + 29 new)','Extension'),('naive_imputed','PBISC imputed as truth','Deconv')]
    y=0; yt=[]; yl=[]
    for f,lab in feats:
        for k,(e,el,c) in enumerate(est):
            r=P[(P.feature==f)&(P.estimator==e)].iloc[0]; sig=r.p<0.05; m=['D','D','s','x'][k]
            ax.plot([r.ci_lo,r.ci_hi],[y,y],color=COL[c],lw=2.6 if k!=1 else 1.6,alpha=1 if k!=1 else 0.6)
            ax.scatter(r.slope,y,s=55,marker=m if m!='x' else 'X',color=COL[c] if sig else 'white',edgecolors='black' if sig else COL[c],linewidths=0.8 if sig else 1.5,zorder=3)
            y+=1
        yt.append(y-2.5); yl.append(lab.replace(' (Real null)','')); y+=0.8
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_ylim(y-0.5,-0.8)
    ax.set_xlabel('Slope per severity step\nlog10(% of cells), 95% CI')
    h=[Line2D([0],[0],color=COL[c],marker=['D','D','s','X'][k],lw=2.4,mec='black',label=el) for k,(e,el,c) in enumerate(est)]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.35,-0.42),fontsize=11,borderaxespad=0,handlelength=1.6,ncol=2,columnspacing=0.8)

GRP=['ctrl','0','1','2','3']; GRPL=['Indep.\nctrl','Mild','Mod.','Sev.','Crit.']
def gzmb_design_strip(ax, gene, level='NK CD56bright', legend=False):
    U=pd.read_csv(D3+'gzmb_prf1_units_v3.csv'); R=pd.read_csv(D3+'gzmb_prf1_designs_v3.csv'); U=U[U.level==level]
    rng=np.random.default_rng(5)
    for src,off in [('Real',-0.2),('PBISC',0.2)]:
        d=U[U.source==src]
        ctrl=d[(d.control_type=='independent control')&(d.centre==('NCL' if src=='Real' else 'w1'))]
        pos=d[(d.covid==1)&(d.centre==('NCL' if src=='Real' else 'w1'))] if False else d[d.covid==1]
        groups=[ctrl[gene].values]+[pos[pos.severity_ord==k][gene].values for k in range(4)]
        for i,v in enumerate(groups):
            ax.scatter(i+off+rng.uniform(-0.07,0.07,len(v)),v,s=34,marker=SRC_MK[src],color=SRC_COL[src],alpha=0.85,edgecolors='white',linewidths=0.5,zorder=3)
            if len(v): ax.plot([i+off-0.13,i+off+0.13],[np.median(v)]*2,color='black',lw=2.4,zorder=4)
    # within-patient pairs (Real, Cambridge): pre-infection -> infection
    w=U[(U.source=='Real')&(U.patient_id.isin(U[(U.source=='Real')&(U.control_type=='pre-infection baseline')].patient_id))].pivot_table(index='patient_id',columns='covid',values=gene).dropna()
    for _,r in w.iterrows():
        ax.plot([5.5,6.4],[r[0],r[1]],color=PALETTE['gray_1'],lw=1.1,alpha=0.8,zorder=2)
    ax.scatter(np.full(len(w),5.5),w[0],s=30,marker='D',color='white',edgecolors=SRC_COL['Real'],linewidths=1.3,zorder=3)
    ax.scatter(np.full(len(w),6.4),w[1],s=30,marker='D',color=SRC_COL['Real'],edgecolors='white',linewidths=0.5,zorder=3)
    ax.axvspan(-0.5,0.5,color='#F2F2F2',zorder=0,lw=0); ax.axvline(4.95,color=PALETTE['neutral'],lw=1.2)
    ax.set_xticks(list(range(5))+[5.5,6.4]); ax.set_xticklabels(GRPL+['Pre-\ninf.','Infec-\nted'],fontsize=12); ax.set_xlim(-0.55,6.85)
    ax.set_ylabel(f'{gene}, NK CD56bright\nlog$_2$(CPM+1)')
    lo,hi=np.nanmin(U[gene]),np.nanmax(U[gene]); rg=hi-lo; ax.set_ylim(lo-0.05*rg,hi+0.85*rg)
    gR=R[(R.level==level)&(R.gene==gene)]
    a=gR[(gR.source=='Real')&gR.design.str.startswith('A')].iloc[0]; b=gR[(gR.source=='Real')&gR.design.str.startswith('B')].iloc[0]; c=gR[(gR.source=='PBISC')&gR.design.str.startswith('A')].iloc[0]
    t=gR[(gR.source=='Real')&gR.design.str.startswith('ordered')].iloc[0]
    ax.text(0.01,0.985,f'Real, NCL 2020: AUC {a.AUC:.2f}, P={_p(a.p)}\nReal, within patient (n={int(b.n_ctrl)}): P={_p(b.p)}\nReal, trend ρ={t.AUC:+.2f}, P={_p(t.p)}',transform=ax.transAxes,va='top',fontsize=11,color=SRC_COL['Real'])
    ax.text(0.01,0.755,f'PBISC, wave 1: AUC {c.AUC:.2f}, P={_p(c.p)}',transform=ax.transAxes,va='top',fontsize=11,color=SRC_COL['PBISC'])
    if legend:
        h=[_h(SRC_COL[s],SRC_MK[s],s) for s in ['Real','PBISC']]; ax.legend(handles=h,loc='upper right',fontsize=11,borderaxespad=0.1,handletextpad=0.2)

LADDER=[('Real','NK CD56bright','A','Real · NK CD56bright'),('Real','All NK','A','Real · all NK'),('Real','Whole sample','A','Real · whole sample'),
        ('SS','Real scRNA Whole sample','','Real · whole sample, same day*'),('SS','bulk, same samples','','Bulk · same samples*'),
        ('Bulk','Whole sample (bulk)','A','Bulk · PBISC input'),('PBISC','Whole sample','A','PBISC · whole sample'),('PBISC','All NK','A','PBISC · all NK'),('PBISC','NK CD56bright','A','PBISC · NK CD56bright')]
def signal_ladder(ax):
    R=pd.read_csv(D3+'gzmb_prf1_designs_v3.csv'); SS=pd.read_csv(D3+'same_sample_real_vs_bulk_v3.csv'); BB=pd.read_csv(D3+'bulk_input_gzmb_prf1_stats_v3.csv')
    cmap={'Real':COL['Real'],'SS':COL['Real'],'Bulk':COL['Bulk'],'PBISC':COL['PBISC']}
    for i,(src,lev,des,lab) in enumerate(LADDER):
        for gene,m,off in [('GZMB','o',-0.15),('PRF1','s',0.15)]:
            if src in ('Real','PBISC'): r=R[(R.source==src)&(R.level==lev)&(R.gene==gene)&R.design.str.startswith('A')].iloc[0]
            elif src=='Bulk': r=BB[(BB.gene==gene)&BB.design.str.startswith('A')].iloc[0]
            else: r=SS[(SS.level==lev)&(SS.gene==gene)].iloc[0]
            c=cmap[src] if not (src=='SS' and lev.startswith('bulk')) else COL['Bulk']; sig=r.p<0.05
            ax.scatter(r.AUC,i+off,s=75,marker=m,color=c if sig else 'white',edgecolors='black' if sig else c,linewidths=0.8 if sig else 1.6,zorder=3)
    ax.axvline(0.5,color=PALETTE['gray_1'],lw=1.2,ls='--'); ax.axhspan(2.5,4.5,color='#F4F4F4',zorder=0,lw=0)
    ax.set_yticks(range(len(LADDER))); ax.set_yticklabels([l for *_,l in LADDER],fontsize=12); ax.set_ylim(len(LADDER)-0.5,-0.6); ax.set_xlim(0.3,1.0); ax.set_xticks([0.4,0.6,0.8])
    ax.set_xlabel('AUC, COVID+ vs control')
    h=[_h(PALETTE['gray_2'],'o','GZMB'),_h(PALETTE['gray_2'],'s','PRF1'),_h(PALETTE['gray_2'],'o','open: P ≥ 0.05',fc='white')]
    ax.legend(handles=h,loc='upper left',fontsize=11,borderaxespad=0.1,handletextpad=0.1)

# ---------------- supplementary
def nk_rho_ext(ax):
    T=pd.read_csv(D3+'nk_clusters_extension_v3.csv'); T['cluster']=T.cluster.astype(str)
    N=pd.read_csv('data/nk_v2_cluster_table.csv',index_col=0); N.index=N.index.astype(str)
    off={'Real':-0.24,'PBISC_only':0,'Extension':0.24}
    for i,c in enumerate(T.cluster):
        r=T.set_index('cluster').loc[c]; xs=[r[f'{k}_rho'] for k in off]; ax.plot([min(xs),max(xs)],[i,i],color=PALETTE['neutral'],lw=2,zorder=1)
        for k in off:
            sig=r[f'{k}_q']<0.1; ax.scatter(r[f'{k}_rho'],i+off[k],s=80,marker=MK[k],color=COL[k] if sig else 'white',edgecolors='black' if sig else COL[k],linewidths=0.8 if sig else 1.6,zorder=3)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.set_xlim(-0.65,0.65)
    ax.set_yticks(range(len(T))); ax.set_yticklabels([f'{c}  {N.loc[c,"name"]}' for c in T.cluster],fontsize=12); ax.set_ylim(len(T)-0.5,-0.6)
    ax.set_xlabel('ρ with peak severity (% of NK cells)')
    h=[_h(COL['Real'],'D','Real (39)'),_h(COL['PBISC_only'],'o','PBISC-only (29)'),_h(COL['Extension'],'s','Extension (68)'),_h(PALETTE['gray_2'],'o','open: BH q ≥ 0.1',fc='white')]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.3,-0.26),ncol=2,fontsize=11,borderaxespad=0,columnspacing=0.8,handletextpad=0.2)

def design_auc(ax, gene):
    R=pd.read_csv(D3+'gzmb_prf1_designs_v3.csv'); R=R[R.gene==gene]
    levs=['NK CD56bright','CD56dim NK','All NK','Whole sample']; rows=[(s,l) for s in ['Real','PBISC'] for l in levs]
    des=[('A','Contemporaneous controls','o'),('all','All independent controls','s'),('B','Within patient','^')]
    for i,(s,l) in enumerate(rows):
        for k,(d,lab,m) in enumerate(des):
            x=R[(R.source==s)&(R.level==l)&R.design.str.startswith(d)]
            if not len(x): continue
            x=x.iloc[0]
            if d=='B' and s=='PBISC': continue
            sig=x.p<0.05; ax.scatter(x.AUC,i+(k-1)*0.22,s=70,marker=m,color=SRC_COL[s] if sig else 'white',edgecolors='black' if sig else SRC_COL[s],linewidths=0.8 if sig else 1.5,zorder=3)
    ax.axvline(0.5,color=PALETTE['gray_1'],ls='--',lw=1.2); ax.axhline(3.5,color=PALETTE['neutral'],lw=1.2)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{s} · {l}' for s,l in rows],fontsize=12); ax.set_ylim(len(rows)-0.5,-0.6); ax.set_xlim(0.25,1.02)
    for t,(s,_) in zip(ax.get_yticklabels(),rows): t.set_color(SRC_COL[s])
    ax.set_xlabel('AUC (within patient: fraction of pairs up)'); ax.set_title(gene,loc='left',fontsize=14,pad=4)
    h=[_h(PALETTE['gray_2'],m,lab) for _,lab,m in des]+[_h(PALETTE['gray_2'],'o','open: P ≥ 0.05',fc='white')]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.4,-0.2),ncol=2,fontsize=11,borderaxespad=0,columnspacing=0.8,handletextpad=0.2)

def paired_lines(ax, gene, level='NK CD56bright'):
    U=pd.read_csv(D3+'gzmb_prf1_units_v3.csv'); U=U[(U.level==level)&(U.source=='Real')]
    pre=U[U.control_type=='pre-infection baseline'].patient_id
    w=U[U.patient_id.isin(pre)].pivot_table(index='patient_id',columns='covid',values=gene).dropna()
    sv=U[U.covid==1].drop_duplicates('patient_id').set_index('patient_id').severity_ord
    for p_,r in w.iterrows():
        c=SEVC.get(int(sv.get(p_,0)),PALETTE['gray_1']); ax.plot([0,1],[r[0],r[1]],color=c,lw=1.8,marker='D',ms=6,mec='black',mew=0.5)
    from scipy import stats
    t=stats.wilcoxon(w[1],w[0]); ax.text(0.03,0.97,f'n={len(w)}; median Δ {np.median(w[1]-w[0]):+.2f}\nWilcoxon P={_p(t.pvalue)}',transform=ax.transAxes,va='top',fontsize=11.5)
    ax.set_xticks([0,1]); ax.set_xticklabels(['Pre-infection\n(Apr–May 2020)','Infection\n(Jan–Feb 2021)'],fontsize=11.5); ax.set_xlim(-0.35,1.35)
    ax.set_ylabel(f'{gene}, Real NK CD56bright\nlog$_2$(CPM+1)'); lo,hi=w.min().min(),w.max().max(); ax.set_ylim(lo-0.2,hi+0.35*(hi-lo))

def bulk_vs_real(ax, gene):
    L=pd.read_csv(D3+'same_sample_real_bulk_Whole_sample.csv',index_col=0); SS=pd.read_csv(D3+'same_sample_real_vs_bulk_v3.csv')
    for cv,c,lab in [(0,PALETTE['neutral'],'COVID−'),(1,PALETTE['red_strong'],'COVID+')]:
        d=L[L.covid==cv]; ax.scatter(d[gene],d[gene+'_bulk'],s=26,color=c,edgecolors='black',linewidths=0.4,alpha=0.85,label=lab,zorder=3)
    a=SS[(SS.gene==gene)&(SS.level=='Real scRNA Whole sample')].iloc[0]; b=SS[(SS.gene==gene)&(SS.level=='bulk, same samples')].iloc[0]
    ax.text(0.03,0.98,f'{int(a.n_pos)} COVID+ / {int(a.n_neg)} COVID− samples\nρ(scRNA, bulk) = {a.rho_real_vs_bulk:.2f}\nAUC scRNA {a.AUC:.2f} (P={_p(a.p)})\nAUC bulk {b.AUC:.2f} (P={_p(b.p)})',transform=ax.transAxes,va='top',fontsize=10.5)
    lo,hi=L[gene+'_bulk'].min(),L[gene+'_bulk'].max(); ax.set_ylim(lo-0.1,hi+0.75*(hi-lo))
    ax.set_xlabel(f'{gene}, scRNA pseudobulk'); ax.set_ylabel(f'{gene}, bulk (same day)')
    if gene=='PRF1': ax.legend(loc='lower right',fontsize=11,borderaxespad=0.1,handletextpad=0.1)

def ext_power_all(ax):
    P=pd.read_csv(D3+'extension_power_v3.csv'); P=P[P.n_real==20]; A=pd.read_csv(D3+'composition_association_v3.csv').set_index('feature')
    feats=[f for f in FNAME if f in set(P.feature)]
    for i,f in enumerate(feats):
        r=P[P.feature==f].iloc[0]
        for col,k in [('real_alone','Real'),('plus_pbisc_only','Extension'),('plus_shuffled_pbisc','Shuffled')]:
            ax.scatter(r[col],i,s=70,marker=MK[k],color=COL[k],edgecolors='black',linewidths=0.6,zorder=3)
        ax.plot([r.real_alone,r.plus_pbisc_only],[i,i],color=PALETTE['neutral'],lw=1.8,zorder=1)
        ax.text(1.08,i,'Real' if r.real_full_p<0.05 else 'null',fontsize=10.5,va='center',color=PALETTE['gray_2'] if r.real_full_p<0.05 else PALETTE['gray_1'])
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([FNAME[f] for f in feats],fontsize=12); ax.set_ylim(len(feats)-0.5,-0.6); ax.set_xlim(-0.03,1.2)
    ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xlabel('Detection rate with 20 Real patients\n(Real-significant: same sign; Real null: any sign)')
    h=[Line2D([0],[0],marker=MK[k],color='none',markerfacecolor=COL[k],markeredgecolor='black',markersize=9,label=l) for k,l in [('Real','Real alone'),('Extension','+ 29 new PBISC'),('Shuffled','+ shuffled')]]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.45,-0.11),ncol=3,fontsize=11,borderaxespad=0,columnspacing=0.8,handletextpad=0.1)

def robustness_rho(ax):
    R=pd.read_csv(D3+'robustness_v3.csv')
    an=[('Real','all deposited samples','Real, all samples','D',COL['Real'],'full'),('Real','published 187-sample roster','Real, 187 roster','D','#7A7A7A','full'),('Real','first sample per patient','Real, first sample','D','#BDBDBD','full'),
        ('PBISC','all samples','PBISC adapter 1','o',COL['PBISC'],'full'),('PBISC','final model, adapter draw 0, all cells (kNN from draw-1 labels)','PBISC adapter 0','o',COL['Draw0'],'full'),
        ('PBISC','legacy K4000 model, draws 0-4 pooled (500 cells/sample)','Legacy model, 5 adapters','v',PALETTE['gray_1'],'full')]
    off=np.linspace(-0.3,0.3,len(an))
    for i,(f,lab) in enumerate(KEY):
        for k,(s,a,l,m,c,fill) in enumerate(an):
            r=R[(R.source==s)&(R.analysis==a)&(R.feature==f)]
            if not len(r) or pd.isna(r.rho.iloc[0]): continue
            r=r.iloc[0]; sig=r.p<0.05
            ax.scatter(r.rho,i+off[k],s=60,marker=m,color=c if (sig and fill=='full') else 'white',edgecolors='black' if (sig and fill=='full') else c,linewidths=0.7 if (sig and fill=='full') else 1.5,zorder=3)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.set_xlim(-0.75,0.75)
    ax.set_yticks(range(len(KEY))); ax.set_yticklabels([l for _,l in KEY],fontsize=12); ax.set_ylim(len(KEY)-0.5,-0.6); ax.set_xlabel('Spearman ρ with peak severity')
    h=[_h(c,m,l) for s,a,l,m,c,fill in an]+[_h(PALETTE['gray_2'],'o','open: P ≥ 0.05',fc='white')]
    ax.legend(handles=h,loc='center left',bbox_to_anchor=(1.01,0.5),fontsize=11,borderaxespad=0,handletextpad=0.2)

def covariate_p(ax):
    R=pd.read_csv(D3+'robustness_v3.csv'); R=R[R.analysis.isin(['all deposited samples','all samples'])]
    cols=[('p','Unadjusted'),('p_adj_age_sex_centre','+ age, sex,\ncentre'),('p_adj_days','+ days from\nsymptoms')]
    rows=[(s,f,l) for s in ['Real','PBISC'] for f,l in KEY]
    for i,(s,f,l) in enumerate(rows):
        r=R[(R.source==s)&(R.feature==f)].iloc[0]
        for j,(c,_) in enumerate(cols):
            v=-np.log10(r[c]); ax.scatter(j,i,s=30+min(v,5)*70,color=SRC_COL[s] if r[c]<0.05 else 'white',edgecolors=SRC_COL[s],linewidths=1.2,zorder=3)
            ax.text(j+0.18,i,_p(r[c]),fontsize=10,va='center')
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([l for _,l in cols],fontsize=12); ax.xaxis.tick_top(); ax.tick_params(axis='x',length=0)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{s} · {l}' for s,_,l in rows],fontsize=11.5); ax.set_ylim(len(rows)-0.5,-0.6); ax.set_xlim(-0.4,len(cols)-0.2)
    for t,(s,_,_) in zip(ax.get_yticklabels(),rows): t.set_color(SRC_COL[s])
    ax.axhline(len(KEY)-0.5,color=PALETTE['neutral'],lw=1.2); ax.spines['bottom'].set_visible(False)

def deconv_accuracy(ax):
    A=pd.read_csv(D3+'deconvolution_accuracy_matched_v3.csv'); names={'rare_NK CD56 bright_pct':'NK CD56bright','rare_B_ASC_dividing_pct':'B ASC dividing','DC_within_Other_pct':'DC','rare_DC1_pct':'DC1','rare_HSPC_pct':'HSPC','NK_pct':'NK cell','T_pct':'T cell','Mono_pct':'Monocyte','rare_pDC_pct':'pDC','rare_CD14mono_IFN_pct':'CD14 mono IFN','rare_MAIT_pct':'MAIT'}
    feats=list(names)
    for i,f in enumerate(feats):
        for k,(v,m,c) in enumerate([('uncalibrated','^','white'),('calibrated','^',COL['Deconv'])]):
            r=A[(A.variant==v)&(A.feature==f)].iloc[0]; ax.scatter(r.spearman_r,i+(k-0.5)*0.3,s=70,marker=m,color=c,edgecolors=COL['Deconv'] if c=='white' else 'black',linewidths=1.4 if c=='white' else 0.7,zorder=3)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2); ax.set_xlim(-0.3,1.0)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([names[f] for f in feats],fontsize=12); ax.set_ylim(len(feats)-0.5,-0.6)
    ax.set_xlabel('Spearman r, deconvolution vs scRNA\n(111 same-day samples)')
    h=[_h(COL['Deconv'],'^','uncalibrated',fc='white'),_h(COL['Deconv'],'^','calibrated (cross-fitted)')]
    ax.legend(handles=h,loc='upper right',fontsize=11,borderaxespad=0.1,handletextpad=0.2)

def auc_by_severity_indep(ax):
    """AUC of each peak-severity group versus independent controls (NCL 2020 for Real, wave-1 for PBISC), NK CD56bright."""
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    U=pd.read_csv(D3+'gzmb_prf1_units_v3.csv'); U=U[U.level=='NK CD56bright']; rows=[]
    for src in ['Real','PBISC']:
        d=U[U.source==src]; ctrl=d[d.control_type=='independent control']
        for gene in ['GZMB','PRF1']:
            res=[]
            for k in range(4):
                x=d[(d.covid==1)&(d.severity_ord==k)][gene]; u=stats.mannwhitneyu(x,ctrl[gene]); res.append((k,u.statistic/len(x)/len(ctrl),u.pvalue,len(x)))
            q=multipletests([r[2] for r in res],method='fdr_bh')[1]
            for (k,a,p,n),qq in zip(res,q): rows.append(dict(source=src,gene=gene,sev=k,AUC=a,p=p,q=qq,n=n,n_ctrl=len(ctrl)))
    R=pd.DataFrame(rows); R.to_csv(D3+'nkbright_auc_by_severity_indep_v3.csv',index=False)
    lab=['Mild','Mod.','Sev.','Crit.']; order=[(s,k) for s in ['Real','PBISC'] for k in range(4)]
    for i,(s,k) in enumerate(order):
        for gene,m,off in [('GZMB','o',-0.13),('PRF1','s',0.13)]:
            r=R[(R.source==s)&(R.sev==k)&(R.gene==gene)].iloc[0]; sig=r.q<0.05
            ax.scatter(r.AUC,i+off,s=80,marker=m,color=SRC_COL[s] if sig else 'white',edgecolors='black' if sig else SRC_COL[s],linewidths=0.8 if sig else 1.5,zorder=3)
    ax.axvline(0.5,color=PALETTE['gray_1'],ls='--',lw=1.2); ax.axhline(3.5,color=PALETTE['neutral'],lw=1.2)
    ax.set_yticks(range(8)); ax.set_yticklabels([f'{s}: {lab[k]}' for s,k in order]); ax.set_ylim(7.6,-0.6); ax.set_xlim(0.2,1.02)
    for t,(s,_) in zip(ax.get_yticklabels(),order): t.set_color(SRC_COL[s])
    ax.set_xlabel('AUC vs independent controls (NK CD56bright)')
    h=[_h(PALETTE['gray_2'],'o','GZMB'),_h(PALETTE['gray_2'],'s','PRF1'),_h(PALETTE['gray_2'],'o','open: BH q ≥ 0.05',fc='white')]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.4,-0.2),ncol=3,fontsize=11,borderaxespad=0,columnspacing=0.8,handletextpad=0.2)

def nk_ari_compact(ax):
    A=pd.read_csv('data/nk_v2_clustering_summary.csv'); col={'Real':SRC_COL['Real'],'PBISC':SRC_COL['PBISC'],'Joint':PALETTE['teal']}; mk={'Real':'D','PBISC':'o','Joint':'s'}
    for i,r in A.iterrows():
        ax.scatter(r.ARI_split_HVG,i,s=110,marker=mk[r.clustering],color=col[r.clustering],edgecolors='black',linewidths=0.8,zorder=3)
        ax.text(r.ARI_split_HVG+0.04,i,f'{r.ARI_split_HVG:.2f} ({int(r.n_clusters)} cl.)',va='center',fontsize=11.5)
    ax.axvspan(0.5,1.0,color=PALETTE['green_1'],alpha=0.7,lw=0,zorder=0); ax.text(0.53,-0.5,'≥ 0.5',fontsize=11,color='#2E7D32')
    ax.set_xlim(0,0.8); ax.set_xticks([0,0.4,0.8]); ax.set_yticks(range(len(A))); ax.set_yticklabels(list(A.clustering)); ax.set_ylim(len(A)-0.5,-0.8)
    ax.set_xlabel('ARI, HVG halves')
