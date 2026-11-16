"""Ver2 panels (2026-11-15): NK-lineage-only subclustering (Fig. 7e-g, Supp. Fig. 9) and NK CD56bright GZMB/PRF1 in
COVID-negative controls (Fig. 7h,i, Supp. Fig. 12). Same figures4papers style and source colours as panels.py."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from panels import PALETTE, SRC_COL, SRC_MK, DIV, SEVC, fmt_q
D='data/'
GRP=['Neg','0','1','2','3']; GRPL=['COVID−','Mild','Mod.','Sev.','Crit.']
NKC={'0':'#9DBBDD','1':PALETTE['blue_secondary'],'2':PALETTE['green_2'],'3':PALETTE['neutral'],'4':PALETTE['red_strong'],
     '5':PALETTE['green_3'],'6':'#E39B17','7':PALETTE['violet']}
def _T(): T=pd.read_csv(D+'nk_v2_cluster_table.csv',index_col=0); T.index=T.index.astype(str); return T
def _p(p): return f'{p:.0e}'.replace('e-0','e-') if p<0.001 else f'{p:.2g}'
def _umap_axes(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    for sp in ('left','bottom'): ax.spines[sp].set_linewidth(1.5)

def schematic_v2(ax):
    ax.set_xlim(0,100); ax.set_ylim(0,20); ax.axis('off')
    def box(x,y,w,h,txt,fc,ec,fs=13):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.3,rounding_size=1.2',fc=fc,ec=ec,lw=2))
        ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs,linespacing=1.3)
    arr=dict(arrowstyle='-|>',lw=2,color=PALETTE['gray_2'])
    box(3,11,22,8,'Real scRNA-seq\n39 COVID+ | 38 COVID−\n588k cells','#EDEDED',SRC_COL['Real'])
    box(3,1,22,8,'PBISC from bulk\n57 COVID+ | 19 COVID−\n4,000 cells/sample','#E3ECF6',SRC_COL['PBISC'])
    box(30,6,16,8,'28 COVID+ patients\nin both sources\n0 same-day pairs','white',PALETTE['gray_1'],fs=12)
    for y in (15,5): ax.annotate('',xy=(30,10),xytext=(25.8,y),arrowprops=arr)
    box(51,11,21,8,'Per source\nReal | PBISC','white',PALETTE['gray_2'])
    box(51,1,21,8,'Combined (n=68)\nrank-INT within source,\nmean per patient','#E2F1F2',SRC_COL['Combined'],fs=12)
    ax.annotate('',xy=(50.6,15),xytext=(46.4,11),arrowprops=arr); ax.annotate('',xy=(50.6,5),xytext=(46.4,9),arrowprops=arr)
    box(77,11,22,8,'Cell-type fractions\nvs peak WHO severity','white',PALETTE['gray_2'],fs=12)
    box(77,1,22,8,'NK-lineage subclusters;\nNK CD56bright GZMB/PRF1\nvs COVID− controls','white',PALETTE['gray_2'],fs=12)
    for y in (15,5): ax.annotate('',xy=(76.6,y),xytext=(72.4,y),arrowprops=arr)

def nkv2_umap(ax, mode='cluster'):
    U=pd.read_csv(D+'nk_v2_umap.csv.gz',index_col=0); U['leiden']=U.leiden.astype(str); U=U.iloc[np.random.default_rng(0).permutation(len(U))]
    if mode=='cluster':
        ax.scatter(U.umap1,U.umap2,s=1.6,c=U.leiden.map(NKC),alpha=0.6,lw=0,rasterized=True)
        for c in NKC:
            d=U[U.leiden==c]; x,y=d.umap1.median(),d.umap2.median()
            ax.text(x,y,c,fontsize=14,fontweight='bold',ha='center',va='center',color='black',bbox=dict(boxstyle='circle,pad=0.18',fc='white',ec=NKC[c],lw=1.8,alpha=0.9))
        ax.text(0.01,0.99,f'{len(U):,} NK cells\n(Real + PBISC)',transform=ax.transAxes,va='top',fontsize=11.5)
    elif mode=='source':
        for src in ['Real','PBISC']:
            d=U[U.source==src]; ax.scatter(d.umap1,d.umap2,s=1.6,color=SRC_COL[src],alpha=0.45,lw=0,rasterized=True,label=f'{src} ({len(d):,})')
        ax.legend(loc='upper left',fontsize=11.5,markerscale=6,handletextpad=0.2,borderaxespad=0.1)
    elif mode=='fine':
        col={'NK':PALETTE['neutral'],'NK IFN stim':PALETTE['red_strong'],'NK CD56 bright':PALETTE['green_3']}; lab={'NK':'NK (CD56dim)','NK IFN stim':'NK IFN stim','NK CD56 bright':'NK CD56bright'}
        for f in ['NK','NK IFN stim','NK CD56 bright']:
            d=U[U.fine==f]; ax.scatter(d.umap1,d.umap2,s=1.6,color=col[f],alpha=0.55,lw=0,rasterized=True,label=lab[f])
        ax.legend(loc='upper left',fontsize=11.5,markerscale=6,handletextpad=0.2,borderaxespad=0.1,title='label (source annotation)',title_fontsize=11)
    _umap_axes(ax)

def nkv2_rho_dots(ax):
    T=_T(); cl=list(T.index); off={'Real':-0.24,'PBISC':0,'Combined':0.24}
    for i,c in enumerate(cl):
        xs=[T.loc[c,f'{a}_rho'] for a in off]; ax.plot([min(xs),max(xs)],[i,i],color=PALETTE['neutral'],lw=2.2,zorder=1)
        for a in off:
            r,q=T.loc[c,f'{a}_rho'],T.loc[c,f'{a}_q']; sig=q<0.1
            ax.scatter(r,i+off[a],s=85,marker=SRC_MK[a],color=SRC_COL[a] if sig else 'white',edgecolors=SRC_COL[a] if not sig else 'black',linewidths=1.6 if not sig else 0.9,zorder=3)
    ax.axvline(0,color=PALETTE['gray_1'],lw=1.2)
    ax.set_yticks(range(len(cl))); ax.set_yticklabels([f'{c}  {T.loc[c,"name"]}' for c in cl],fontsize=12)
    for t,c in zip(ax.get_yticklabels(),cl): t.set_color('black' if c in ('4','6','7') else PALETTE['gray_1'])
    ax.set_ylim(len(cl)-0.5,-0.6); ax.set_xlim(-0.6,0.6); ax.set_xlabel('ρ with peak severity\n(% of patient NK cells)')
    h=[Line2D([0],[0],marker=SRC_MK[a],color='none',markerfacecolor=SRC_COL[a],markeredgecolor='black',markersize=9,label=a) for a in off]
    h.append(Line2D([0],[0],marker='o',color='none',markerfacecolor='white',markeredgecolor=PALETTE['gray_2'],markeredgewidth=1.6,markersize=9,label='open: BH q ≥ 0.1'))
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.3,-0.42),ncol=2,fontsize=11.5,columnspacing=0.8,handletextpad=0.2,borderaxespad=0)

PWC={'4':PALETTE['red_strong'],'6':'#E39B17','7':PALETTE['violet']}
def nkv2_power(ax, clusters=('4','6','7'), legend=True):
    P=pd.read_csv(D+'nk_v2_power_subsample.csv'); P['cluster']=P.cluster.astype(str); T=_T()
    for c in clusters:
        d=P[P.cluster==c]
        ax.plot(d.n_real,d.power_real,color=PWC.get(c,PALETTE['gray_2']),lw=2.4,marker='D',ms=7,mec='black',mew=0.6)
        ax.plot(d.n_real,d.power_real_plus_pbisc,color=PWC.get(c,PALETTE['gray_2']),lw=2.4,ls='--',marker='s',ms=7,mec='black',mew=0.6)
    ax.axhline(0.8,color=PALETTE['gray_1'],ls=':',lw=1.4)
    ax.set_ylim(-0.03,1.05); ax.set_xlim(8,37); ax.set_xticks([10,20,30]); ax.set_xlabel('Real COVID+ patients sampled')
    ax.set_ylabel('Detection rate\n(p < 0.05, same sign)')
    if legend:
        h=[Line2D([0],[0],color=PWC[c],lw=3,label=f'cl. {c}') for c in clusters]
        l1=ax.legend(handles=h,loc='upper left',bbox_to_anchor=(-0.3,-0.27),ncol=3,fontsize=11.5,columnspacing=0.7,handletextpad=0.3,handlelength=1.3,borderaxespad=0); ax.add_artist(l1)
        h2=[Line2D([0],[0],color='black',marker='D',ms=7,lw=2,label='Real alone'),Line2D([0],[0],color='black',marker='s',ms=7,lw=2,ls='--',label='Real + all PBISC')]
        ax.legend(handles=h2,loc='upper left',bbox_to_anchor=(-0.3,-0.42),ncol=2,fontsize=11.5,columnspacing=0.7,handletextpad=0.3,handlelength=2.2,borderaxespad=0)

def _grp_strip(ax, P, col, ylabel, seed=13, jit=0.07, ms=40):
    rng=np.random.default_rng(seed)
    for src in ['Real','PBISC']:
        d=P[P.source==src]; off=-0.2 if src=='Real' else 0.2
        for i,g in enumerate(GRP):
            v=d.loc[d.group==g,col].values
            ax.scatter(i+off+rng.uniform(-jit,jit,len(v)),v,s=ms,marker=SRC_MK[src],color=SRC_COL[src],alpha=0.85,edgecolors='white',linewidths=0.5,zorder=3)
            if len(v): ax.plot([i+off-0.13,i+off+0.13],[np.median(v)]*2,color='black',lw=2.4,zorder=4)
    ax.axvspan(-0.5,0.5,color='#F2F2F2',zorder=0,lw=0)
    ax.set_xticks(range(5)); ax.set_xticklabels(GRPL); ax.set_xlim(-0.55,4.55); ax.set_ylabel(ylabel)

def nkbright_neg_strip(ax, gene, legend=False):
    P=pd.read_csv(D+'nkbright_gzmb_prf1_neg_vs_covid_patient.csv'); P['group']=P.group.astype(str)
    S=pd.read_csv(D+'nkbright_gzmb_prf1_neg_vs_covid_stats.csv'); S=S[S.gene==gene].set_index('comparison')
    _grp_strip(ax,P,gene,f'{gene} in NK CD56bright\nlog$_2$(CPM+1)')
    lo,hi=np.nanmin(P[gene]),np.nanmax(P[gene]); rng_=hi-lo; ax.set_ylim(lo-0.05*rng_,hi+0.8*rng_)
    for k,src in enumerate(['Real','PBISC']):
        a=S.loc[f'{src}: COVID+ vs Neg']; t=S.loc[f'{src}: trend Neg<0<1<2<3']
        ax.text(0.01+0.42*k,0.98,f'{src}\n+ vs −: AUC {a.stat:.2f}, p={_p(a.p)}\ntrend ρ={t.stat:+.2f}, p={_p(t.p)}',transform=ax.transAxes,va='top',fontsize=11.5,color=SRC_COL[src])
    i=S.loc['COVID effect x source (OLS)']; n=S.loc['PBISC vs Real: Neg']
    ax.text(0.01,0.765,f'COVID− PBISC vs Real: p={_p(n.p)};  COVID × source interaction: p={_p(i.p)}',transform=ax.transAxes,va='top',ha='left',fontsize=11.5,color=PALETTE['gray_2'])
    if legend:
        h=[Line2D([0],[0],marker=SRC_MK[s],color='none',markerfacecolor=SRC_COL[s],markeredgecolor='black',markersize=8,label=s) for s in ['Real','PBISC']]
        ax.legend(handles=h,loc='lower right',fontsize=11.5,borderaxespad=0.1,handletextpad=0.2)

def nkv2_markers(fig, place, w, h):
    M=pd.read_csv(D+'nk_v2_cluster_markers.csv'); M['cluster']=M.cluster.astype(str); T=_T(); cl=list(T.index)
    # pooled across sources (cell-weighted)
    n=T[['Real','PBISC']]
    M['w']=[n.loc[c,s] for c,s in zip(M.cluster,M.source)]
    g=M.groupby(['cluster','gene']).apply(lambda d: pd.Series(dict(m=np.average(d.mean_lognorm,weights=d.w),f=np.average(d.frac_expr,weights=d.w)))).reset_index()
    genes=list(dict.fromkeys(M.gene)); piv=g.pivot(index='cluster',columns='gene',values='m').loc[cl,genes]; fr=g.pivot(index='cluster',columns='gene',values='f').loc[cl,genes]
    z=(piv-piv.mean())/piv.std(ddof=0)
    ax=place(3.6,1.35,w-3.6-1.5,h-1.35-0.25); cax=place(w-1.15,1.6,0.22,h-1.6-1.3)
    for i,c in enumerate(cl):
        sc=ax.scatter(np.arange(len(genes)),np.full(len(genes),i),s=15+fr.loc[c].values*230,c=z.loc[c].values,cmap=DIV,vmin=-2,vmax=2,edgecolors='black',linewidths=0.5,zorder=3)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes,rotation=55,ha='right',fontstyle='italic',fontsize=12.5)
    ax.set_yticks(range(len(cl))); ax.set_yticklabels([f'{c}  {T.loc[c,"name"]}' for c in cl],fontsize=12.5); ax.set_ylim(len(cl)-0.5,-0.6); ax.set_xlim(-0.6,len(genes)-0.4)
    cb=plt.colorbar(sc,cax=cax); cb.set_label('z across clusters',fontsize=12); cb.ax.tick_params(labelsize=11); cb.outline.set_linewidth(1)
    cax.text(0.5,1.05,'size = %\nexpressing',transform=cax.transAxes,ha='center',va='bottom',fontsize=11.5)

def nkv2_composition(ax):
    T=_T(); cl=list(T.index); y=np.arange(len(cl))
    parts=[('frac_NK','NK (CD56dim)',PALETTE['neutral']),('frac_NK IFN stim','NK IFN stim',PALETTE['red_strong']),('frac_NK CD56 bright','NK CD56bright',PALETTE['green_3'])]
    left=np.zeros(len(cl))
    for c,l,col in parts:
        v=T[c].values*100; ax.barh(y,v,left=left,height=0.7,color=col,edgecolor='black',linewidth=0.8,label=l); left+=v
    for i,c in enumerate(cl):
        ax.text(102,i,f'{T.loc[c,"real_share"]*100:.0f}%',va='center',fontsize=11.5)
    ax.text(102,-0.95,'Real\nshare',fontsize=11.5,va='bottom')
    ax.set_yticks(y); ax.set_yticklabels(cl); ax.set_ylim(len(cl)-0.5,-0.6); ax.set_xlim(0,100); ax.set_xlabel('% of cluster (source label)'); ax.set_ylabel('Joint NK cluster')
    ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.24),ncol=2,fontsize=11.5,borderaxespad=0,columnspacing=0.8)

def nkv2_ari(ax):
    A=pd.read_csv(D+'nk_v2_clustering_summary.csv'); R=pd.read_csv(D+'rare_subclustering_summary.csv')
    col={'Real':SRC_COL['Real'],'PBISC':SRC_COL['PBISC'],'Joint':SRC_COL['Combined']}; mk={'Real':'D','PBISC':'o','Joint':'s'}
    for i,r in A.iterrows():
        ax.scatter(r.ARI_split_HVG,i,s=120,marker=mk[r.clustering],color=col[r.clustering],edgecolors='black',linewidths=0.8,zorder=3)
        ax.text(r.ARI_split_HVG+0.03,i,f'{r.ARI_split_HVG:.2f} ({int(r.n_clusters)} cl., {int(r.n_cells):,} cells)',va='center',fontsize=11.5)
    ax.axvspan(0.5,1.0,color=PALETTE['green_1'],alpha=0.7,lw=0,zorder=0); ax.text(0.52,-0.45,'ARI ≥ 0.5',fontsize=11.5,color='#2E7D32')
    ax.set_xlim(0,0.75); ax.set_yticks(range(len(A))); ax.set_yticklabels([f'{c} clustering' for c in A.clustering]); ax.set_ylim(len(A)-0.5,-0.7)
    ax.set_xlabel('ARI, two random halves of HVGs')

def nkv2_resolution(ax):
    R=pd.read_csv(D+'nk_v2_resolution_sensitivity.csv')
    for a,c in [('Real','Real_q_lt_0p1'),('PBISC','PBISC_q_lt_0p1'),('Combined','Combined_q_lt_0p1')]:
        ax.plot(range(len(R)),R[c],color=SRC_COL[a],lw=2.2,marker=SRC_MK[a],ms=9,mec='black',mew=0.7,label=a)
    ax.set_xticks(range(len(R))); ax.set_xticklabels([f'res {r.resolution} ({int(r.n_clusters)} cl.)\nρ Real~PBISC {r.rho_concordance_Real_vs_PBISC:.2f}' for _,r in R.iterrows()],fontsize=12); ax.set_ylim(-0.3,R[['Real_q_lt_0p1','PBISC_q_lt_0p1','Combined_q_lt_0p1']].values.max()+1.2)
    ax.set_xlim(-0.4,len(R)-0.6); ax.set_ylabel('Clusters with BH q < 0.1'); ax.legend(loc='upper left',fontsize=11.5,borderaxespad=0.1,handletextpad=0.3)

def nkv2_frac_strip(ax, c, legend=False):
    F=pd.read_csv(D+'nk_v2_cluster_frac_patient.csv'); F['group']=F.group.astype(str); T=_T(); R=pd.read_csv(D+'nk_v2_cluster_neg_vs_covid.csv'); R['cluster']=R.cluster.astype(str)
    _grp_strip(ax,F,c,'% of patient NK cells',seed=int(c)+3,ms=32)
    ax.set_title(f'cl. {c}  {T.loc[c,"name"]}',loc='left',fontsize=13.5,pad=4,color=PWC.get(c,'black'))
    top=np.nanmax(F[c]); ax.set_ylim(0,top*1.85)
    for k,src in enumerate(['Real','PBISC']):
        r=R[(R.cluster==c)&(R.source==src)].iloc[0]; t=T.loc[c]
        ax.text(0.01,0.98-0.15*k,f'{src}: + vs − AUC {r.AUC_pos_vs_neg:.2f}, q={fmt_q(r.q)}\nseverity ρ={t[f"{src}_rho"]:+.2f}, q={fmt_q(t[f"{src}_q"])}',transform=ax.transAxes,va='top',fontsize=11,color=SRC_COL[src])
    if legend:
        h=[Line2D([0],[0],marker=SRC_MK[s],color='none',markerfacecolor=SRC_COL[s],markeredgecolor='black',markersize=8,label=s) for s in ['Real','PBISC']]
        ax.legend(handles=h,loc='center right',fontsize=11.5,borderaxespad=0.1,handletextpad=0.2)

def nkbright_auc_dots(ax):
    S=pd.read_csv(D+'nkbright_gzmb_prf1_neg_vs_covid_stats.csv'); off={'GZMB':-0.13,'PRF1':0.13}
    rows=[(src,lv) for src in ['Real','PBISC'] for lv in ['0','1','2','3']]
    for gene,mk_ in [('GZMB','o'),('PRF1','s')]:
        for i,(src,lv) in enumerate(rows):
            r=S[(S.gene==gene)&(S.comparison==f'{src}: sev {lv} vs Neg')].iloc[0]; sig=r.q_within<0.05
            ax.scatter(r.stat,i+off[gene],s=95,marker=mk_,color=SRC_COL[src] if sig else 'white',edgecolors='black' if sig else SRC_COL[src],linewidths=0.9 if sig else 1.6,zorder=3)
    ax.axvline(0.5,color=PALETTE['gray_1'],lw=1.2,ls='--'); ax.axhline(3.5,color=PALETTE['neutral'],lw=1.2)
    ax.set_yticks(range(8)); ax.set_yticklabels([f'{src}: {GRPL[int(lv)+1]}' for src,lv in rows]); ax.set_ylim(7.6,-0.6); ax.set_xlim(0.25,1.02)
    for t,(src,_) in zip(ax.get_yticklabels(),rows): t.set_color(SRC_COL[src])
    ax.set_xlabel('AUC vs COVID− controls (NK CD56bright)')
    h=[Line2D([0],[0],marker='o',color='none',markerfacecolor=PALETTE['gray_2'],markeredgecolor='black',markersize=9,label='GZMB'),Line2D([0],[0],marker='s',color='none',markerfacecolor=PALETTE['gray_2'],markeredgecolor='black',markersize=9,label='PRF1'),
       Line2D([0],[0],marker='o',color='none',markerfacecolor='white',markeredgecolor=PALETTE['gray_2'],markeredgewidth=1.6,markersize=9,label='open: BH q ≥ 0.05')]
    ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.4,-0.2),ncol=3,fontsize=11.5,borderaxespad=0,columnspacing=0.8,handletextpad=0.2)

def nkbright_same_patient(ax, gene):
    P=pd.read_csv(D+'nkbright_gzmb_prf1_neg_vs_covid_patient.csv'); S=pd.read_csv(D+'nkbright_gzmb_prf1_neg_vs_covid_stats.csv')
    w=P[P.covid==1].pivot_table(index='patient_id',columns='source',values=gene).dropna()
    sv=P[P.covid==1].drop_duplicates('patient_id').set_index('patient_id').severity_ord.reindex(w.index)
    for k in range(4):
        m=sv==k; ax.scatter(w.Real[m],w.PBISC[m],s=70,color=SEVC[k],edgecolors='black',linewidths=0.7,zorder=3,label=GRPL[k+1])
    lo=min(w.min().min(),w.min().min())-0.3; hi=max(w.max().max(),w.max().max())+0.3; ax.plot([lo,hi],[lo,hi],color=PALETTE['neutral'],ls='--',lw=1.4)
    ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_xlabel(f'Real {gene}'); ax.set_ylabel(f'PBISC {gene}')
    r=S[(S.gene==gene)&(S.comparison=='PBISC vs Real same patient: COVID+')].iloc[0]
    ax.text(0.03,0.97,f'COVID+, both sources (n={r.n})\nρ={r.stat:.2f}\nPBISC − Real {r["diff"]:+.2f}\np={_p(r.p)}',transform=ax.transAxes,va='top',fontsize=11.5)
    if gene=='PRF1': ax.legend(loc='lower right',fontsize=11,borderaxespad=0.1,handletextpad=0.1,title='peak',title_fontsize=11)
