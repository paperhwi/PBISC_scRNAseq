"""Figure 7 (ver5): B ASC dividing, NK CD56bright and HSPC after adding 29 bulk-only patients through PBISC.
Style: figures4papers (github.com/ChenLiu-1996/figures4papers, scientific-figure-making skill: PALETTE, spines off,
frameless legends, dedicated legend panel, trend lines with uncertainty bands, dpi 300) and page layout of
PBISC_Figure7_SuppFig6-13_panels_ver4.pptx (178 x 258 mm page, full-width design banner, dense varied panels).
Each panel = own canvas at K x print size (font 15 -> ~6 pt placed, as in ver4); positions in panels_out/manifest.json."""
import os, sys, json; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from scipy import stats
from style import apply_publication_style, PALETTE, finalize_figure
apply_publication_style()
G=os.environ.get('GAIN','/home/claude/gain/')+'data/'; O='panels_out/'; K=2.5  # G = s08 outputs
os.makedirs('data',exist_ok=True); os.makedirs(O,exist_ok=True)
COL=dict(Real=PALETTE['gray_3'],New=PALETTE['blue_secondary'],Ext=PALETTE['blue_main'],Shuf=PALETTE['neutral'])
MK=dict(Real='D',New='o',Ext='s',Shuf='s')
NAME=dict(Real='Real',New='PBISC new',Ext='Real + PBISC new')
SETS={'Real (39)':'Real','PBISC new (29)':'New','Real + PBISC new (68)':'Ext'}
POP=[('rare_B_ASC_dividing_pct','B ASC dividing','B_ASC_dividing'),('rare_NK CD56 bright_pct','NK CD56bright','NK CD56 bright'),('rare_HSPC_pct','HSPC','HSPC')]
def fp(p): return 'P < 1e-4' if p<1e-4 else (f'P = {p:.0e}'.replace('e-0','e-') if p<0.001 else f'P = {p:.3f}' if p<0.1 else f'P = {p:.2f}')
SV=pd.read_csv(G+'gain_severity_rho.csv'); AU=pd.read_csv(G+'gain_death_auc.csv'); PW=pd.read_csv(G+'gain_power_curves.csv')
S=pd.read_csv(os.environ.get('RUN','/home/claude/run/')+'out/sample_level_composition_metadata.csv'); FEATS=[p[0] for p in POP]
def patient_table(src):  # identical to gain_analysis.py
    d=S[(S.source==src)&(S.condition=='positive')&S.in_primary_scope.astype(bool)&S.qc_pass_min200cells.astype(bool)]
    g=d.groupby('patient_id'); T=g[FEATS].mean(); T['sev']=g.severity_ord.first(); T['fatal']=g.fatal_harmonized.first(); return T
R=patient_table('Real'); B=patient_table('PBISC'); NEW=B.loc[~B.index.isin(R.index)]
def rint(x): r=x.rank(); return pd.Series(stats.norm.ppf((r-0.5)/x.notna().sum()),index=x.index)
PT=pd.concat([pd.DataFrame({**{f:rint(R[f]) for f in FEATS},'sev':R.sev,'fatal':R.fatal,'src':'Real'}),
              pd.DataFrame({**{f:rint(NEW[f]) for f in FEATS},'sev':NEW.sev,'fatal':NEW.fatal,'src':'New'})])
PT.to_csv('data/patient_normal_scores.csv'); rng=np.random.default_rng(7); man=[]
def canvas(w,h): return plt.figure(figsize=(K*w,K*h))
def save(fig,name,letter,x,y,w,h): finalize_figure(fig,O+name); man.append(dict(file=name+'.png',label=letter,x=x,y=y,w=w,h=h))
def strip(ax,x0,vals,key):
    x=x0+rng.uniform(-0.08,0.08,len(vals)); ax.scatter(x,vals,s=38,marker=MK[key],c=COL[key],edgecolors='white',linewidths=0.6,zorder=3)
    if len(vals): ax.plot([x0-0.15,x0+0.15],[np.median(vals)]*2,color='black',lw=2.4,zorder=4)
def stat_text(ax,rows,kind):
    for i,(k,r) in enumerate(rows):
        v=f"ρ = {r.rho:.2f}" if kind=='rho' else f"AUC = {r.auc:.2f}"
        ax.text(0.02,0.985-0.085*i,f"{NAME[k]}: {v}, {fp(r.p)}",transform=ax.transAxes,fontsize=12.5,color=COL[k],va='top',fontweight='bold' if r.p<0.05 else 'normal')
# ================= a: design banner (ver4 style: two rows of rounded boxes) =================
AW,AH=6.71,1.10
fig=canvas(AW,AH); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis('off')
def box(x,y,w,h,t,fc,ec):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.5,rounding_size=1.6',fc=fc,ec=ec,lw=2.2)); ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=14,linespacing=1.3)
arr=dict(arrowstyle='-|>',lw=2.2,color=PALETTE['gray_2'],mutation_scale=22)
bw,bh=21.5,38; xs=[1,26,51,76]; yt,yb=57,6
box(xs[0],yt,bw,bh,'Real scRNA-seq\n39 COVID+ patients, 7 deaths','#EDEDED',COL['Real'])
box(xs[1],yt,bw,bh,'Real + PBISC new: 68 patients\n13 deaths; rank-INT within source\n28 shared patients: Real data only','#DCE6F2',COL['Ext'])
box(xs[2],yt,bw,bh,'Patient-level association\npeak WHO severity (Spearman ρ)\ndeath (AUC)','white',PALETTE['gray_2'])
box(xs[3],yt,bw,bh,'Power by subsampling\n1,000 draws of 10–35 Real patients\n± all 29 PBISC new patients','white',PALETTE['gray_2'])
box(xs[0],yb,bw,bh,'Bulk RNA-seq only\n29 new COVID+ patients, 6 deaths','#E3ECF6',COL['New'])
box(xs[1],yb,bw,bh,'PBISC\nbulk → generated single cells','white',COL['New'])
box(xs[2],yb,bw,bh,'B ASC dividing · NK CD56bright · HSPC\n(% of PBMCs, patient mean of\nsamples ≥ 200 cells)','white',PALETTE['gray_2'])
box(xs[3],yb,bw,bh,'Control: severity labels of the\nPBISC new patients shuffled','white',PALETTE['red_strong'])
for i in range(3): ax.annotate('',xy=(xs[i+1]-0.8,yt+bh/2),xytext=(xs[i]+bw+0.8,yt+bh/2),arrowprops=arr)
ax.annotate('',xy=(xs[1]-0.8,yb+bh/2),xytext=(xs[0]+bw+0.8,yb+bh/2),arrowprops=arr)
ax.annotate('',xy=(xs[1]+bw/2,yt-0.8),xytext=(xs[1]+bw/2,yb+bh+0.8),arrowprops=arr)
ax.annotate('',xy=(xs[2]+bw/2,yt-0.8),xytext=(xs[2]+bw/2,yb+bh+0.8),arrowprops=arr)
ax.annotate('',xy=(xs[3]+bw/2,yt-0.8),xytext=(xs[3]+bw/2,yb+bh+0.8),arrowprops=arr)
save(fig,'Fig7_a_design','a',0.15,0.42,AW,AH)
# ================= b, c: forests =================
def forest(val,lo,hi,T,xlab,ref,name,letter,x,w,h):
    fig=canvas(w,h); ax=fig.add_axes([0.40,0.17,0.55,0.79]); y=0; ticks=[]
    for f,lab,_ in POP:
        for s_,k in SETS.items():
            r=T[(T.feature==lab)&(T.set==s_)].iloc[0]
            ax.plot([r[lo],r[hi]],[y,y],color=COL[k],lw=2.6); ax.scatter(r[val],y,s=80,marker=MK[k],zorder=3,c=COL[k] if r.p<0.05 else 'white',edgecolors=COL[k],linewidths=2)
            y-=1
        ticks.append(y+2); y-=0.8
    ax.axvline(ref,color=PALETTE['gray_1'],ls=':',lw=1.6); ax.set_yticks(ticks); ax.set_yticklabels([p[1] for p in POP]); ax.set_xlabel(xlab)
    ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False); ax.set_ylim(y+0.3,0.8)
    ax.set_xlim(*((-0.02,1.04) if ref==0.5 else (-0.85,0.95))); ax.set_xticks([0,0.5,1] if ref==0.5 else [-0.5,0,0.5])
    save(fig,name,letter,x,1.64,w,h)
forest('rho','lo','hi',SV,'Spearman ρ with peak WHO severity',0,'Fig7_b_severity_forest','b',0.10,2.25,2.0)
forest('auc','lo','hi',AU,'AUC for death',0.5,'Fig7_c_death_forest','c',2.35,2.25,2.0)
# ================= d: Real patients needed for 80% power (dumbbell, ver4 g style) =================
w,h=2.30,2.0; fig=canvas(w,h); ax=fig.add_axes([0.36,0.17,0.58,0.79])
def n80(x,y):
    x=np.asarray(x); y=np.asarray(y)
    if y[0]>=0.8: return x[0],True
    k=np.where(y>=0.8)[0]
    if not len(k): return np.nan,False
    i=k[0]; return x[i-1]+(0.8-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]),False
for i,(f,lab,pn) in enumerate(POP):
    d=PW[PW.name==pn]; yv=-i
    a,_=n80(d.n_real,d.real_alone); b,le=n80(d.n_real,d.plus_pbisc_only)
    ra=a if np.isfinite(a) else 40; ax.plot([b,ra],[yv,yv],color=PALETTE['neutral'],lw=3,zorder=1)
    ax.scatter(ra,yv,marker='D',s=90,zorder=3,c=COL['Real'] if np.isfinite(a) else 'white',edgecolors=COL['Real'],linewidths=2)
    ax.scatter(b,yv,marker='s',s=90,zorder=3,c=COL['Ext'],edgecolors=COL['Ext'])
    ax.text(ra,yv+0.22,(f'{a:.0f}' if np.isfinite(a) else ('> 35' if d.real_alone.max()>0.5 else 'not reached')),ha='center',fontsize=12.5,color=COL['Real'])
    ax.text(b,yv+0.22,'≤ 10' if le else f'{b:.0f}',ha='center',fontsize=12.5,color=COL['Ext'],fontweight='bold')
ax.set_yticks([0,-1,-2]); ax.set_yticklabels([p[1] for p in POP]); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
ax.set_xlim(5,44); ax.set_xticks([10,20,30,40]); ax.set_xticklabels(['10','20','30','>35']); ax.set_ylim(-2.6,0.6)
ax.set_xlabel('Real patients for 80% power')
save(fig,'Fig7_d_n_for_80pct_power','d',4.60,1.64,w,h)
# ================= dedicated legend panel =================
fig=canvas(6.7,0.2); ax=fig.add_axes([0,0,1,1]); ax.set_axis_off()
hd=[Line2D([0],[0],marker=MK[k],ls='',ms=10,mfc=COL[k],mec=COL[k],label=l) for k,l in [('Real','Real (39 patients)'),('New','PBISC new (29)'),('Ext','Real + PBISC new (68)')]]
hd+=[Line2D([0],[0],marker='o',ls='',ms=10,mfc='white',mec=PALETTE['gray_2'],mew=2,label='open: P ≥ 0.05 (b, c) / not reached (d)'),
     Line2D([0],[0],color=COL['Shuf'],ls='--',marker='s',ms=8,mec='black',mew=0.5,lw=2.4,label='+ PBISC new, severity shuffled (g, j, m)')]
ax.legend(handles=hd,loc='center',ncol=5,fontsize=14,frameon=False,columnspacing=1.5,handletextpad=0.4)
finalize_figure(fig,O+'Fig7_legend_panel'); man.append(dict(file='Fig7_legend_panel.png',label='',x=0.15,y=3.70,w=6.7,h=0.2))
# ================= per-population rows: severity | death | power =================
Y0=[4.10,6.06,8.02]; RH=1.76; letters=iter('efghijklm')
for (f,lab,pname),y0 in zip(POP,Y0):
    tag=lab.replace(' ','_')
    w=2.60; fig=canvas(w,RH); ax=fig.add_axes([0.21,0.19,0.77,0.79])
    for g in range(4):
        for k,dx in [('Real',-0.19),('New',0.19)]: strip(ax,g+dx,PT[(PT.src==k)&(PT.sev==g)][f].values,k)
    ax.set_xticks(range(4)); ax.set_xticklabels(['Mild','Moderate','Severe','Critical']); ax.set_xlim(-0.55,3.55); ax.set_ylim(-2.6,4.3); ax.set_yticks([-2,-1,0,1,2])
    ax.set_ylabel(f'{lab}\n(normal score)'); stat_text(ax,[(SETS[s_],SV[(SV.feature==lab)&(SV.set==s_)].iloc[0]) for s_ in SETS],'rho')
    save(fig,f'Fig7_{tag}_severity',next(letters),0.10,y0,w,RH)
    w=1.85; fig=canvas(w,RH); ax=fig.add_axes([0.25,0.19,0.72,0.79])
    for g in range(2):
        for k,dx in [('Real',-0.18),('New',0.18)]: strip(ax,g+dx,PT[(PT.src==k)&(PT.fatal==g)][f].values,k)
    n0=int((PT.fatal==0).sum()); n1=int((PT.fatal==1).sum())
    ax.set_xticks([0,1]); ax.set_xticklabels([f'Survived\n({n0})',f'Died\n({n1})']); ax.set_xlim(-0.5,1.5); ax.set_ylim(-2.6,4.3); ax.set_yticks([-2,-1,0,1,2])
    ax.set_ylabel('Normal score'); stat_text(ax,[(SETS[s_],AU[(AU.feature==lab)&(AU.set==s_)].iloc[0]) for s_ in SETS],'auc')
    save(fig,f'Fig7_{tag}_death',next(letters),2.72,y0,w,RH)
    w=2.30; fig=canvas(w,RH); ax=fig.add_axes([0.22,0.19,0.75,0.79]); d=PW[PW.name==pname]
    for col,k,ls in [('real_alone','Real','-'),('plus_pbisc_only','Ext','-'),('plus_shuffled_pbisc','Shuf','--')]:
        p=d[col].values; se=1.96*np.sqrt(p*(1-p)/1000)
        ax.fill_between(d.n_real,np.clip(p-se,0,1),np.clip(p+se,0,1),color=COL[k],alpha=0.18,lw=0)
        ax.plot(d.n_real,p,color=COL[k],lw=2.6,ls=ls,marker=MK[k],ms=7,mec='black',mew=0.5)
    ax.axhline(0.8,color=PALETTE['gray_1'],ls=':',lw=1.5); ax.set_ylim(-0.03,1.42); ax.set_yticks([0,0.5,1]); ax.set_xticks([10,20,30]); ax.set_xlim(8,37)
    ax.set_xlabel('Real patients sampled'); ax.set_ylabel('Power (P < 0.05)')
    r20=d[d.n_real==20].iloc[0]
    for i,(t,c,b) in enumerate([(f"n = 20, Real alone: {r20.real_alone:.0%}",COL['Real'],False),(f"+ PBISC new: {r20.plus_pbisc_only:.0%}",COL['Ext'],True),(f"+ shuffled severity: {r20.plus_shuffled_pbisc:.0%}",PALETTE['gray_1'],False)]):
        ax.text(0.02,0.985-0.085*i,t,transform=ax.transAxes,fontsize=12.5,color=c,va='top',fontweight='bold' if b else 'normal')
    save(fig,f'Fig7_{tag}_power',next(letters),4.62,y0,w,RH)
json.dump(man,open(O+'manifest.json','w'),indent=1); print(len(man),'panels')
