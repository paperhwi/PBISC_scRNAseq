"""Publication-sized interval plots; values read from saved result tables."""
from common import *
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageOps

plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8.5,'ytick.labelsize':9,
                     'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','savefig.facecolor':'white'})
teal='#1F6F5F';ink='#171717';gray='#7B8280'
spec={'width_mm':178,'forms':['horizontal dots and 95% patient bootstrap intervals'],
      'paired_source':'results/paired_accuracy_summary.csv','paired_scope':'common_test_paired, L1, five eligible models in every ensemble',
      'unit':'L1 distance of fractions, 0-2; lower is better','uncertainty':'2,000 bootstrap replicates over unique patients; fixed adapters',
      'style_anchors':['Astra governing tokens','OAI-VIS-004','OAI-VIS-013','OAI-VIS-008'],
      'output_formats':['png','svg'],'palette':{'focal':teal,'baseline':gray,'ink':ink}}
(OUT/'visual_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')

def save(fig,name):
    fig.savefig(OUT/f'figs/{name}.png',dpi=300)
    fig.savefig(OUT/f'figs/{name}.svg')
    ImageOps.grayscale(Image.open(OUT/f'figs/{name}.png')).save(OUT/f'figs/{name}_grayscale.png')
    plt.close(fig)

s=pd.read_csv(OUT/'results/paired_accuracy_summary.csv')
fig,axes=plt.subplots(1,2,figsize=(178/25.4,94/25.4),sharex=True,sharey=True)
methods=['other_patient_null','PBISC_k1000_eligible_ensemble','PBISC_k4000_eligible_ensemble'];names=['Other-patient mean','PBISC K = 1,000','PBISC K = 4,000']
for ax,coh,title in zip(axes,['ye2026_sepsis','imperial_eskd'],['Sepsis','ESKD COVID']):
    ds=s[(s.scope=='common_test_paired')&(s.level=='l1')&(s.draw==-1)&s.cohort.eq(coh)].set_index('method')
    for i,m in enumerate(methods):
        r=ds.loc[m];color=teal if i==2 else ink if i==1 else gray
        ax.errorbar(r.L1,2-i,xerr=[[r.L1-r.L1_lo],[r.L1_hi-r.L1]],fmt=['s','o','D'][i],color=color,ms=4.5,capsize=2,lw=1)
    ax.set_yticks([2,1,0],names);ax.set_ylim(-.5,2.5);ax.set_xlim(.20,.65);ax.set_xlabel('Composition L1 (lower is better)')
    ax.set_title(title+f'\n{int(ds.iloc[0].n_donors)} held-out patients',loc='left',pad=12)
fig.subplots_adjust(left=.25,right=.98,bottom=.30,top=.79,wspace=.30)
fig.text(.25,.95,'Paired composition accuracy',fontsize=12,weight='semibold')
fig.text(.25,.055,'Five-model ensembles; no test patient trained any adapter.\nWhiskers: 95% patient bootstrap intervals; model training uncertainty excluded.',fontsize=8)
save(fig,'paired_accuracy')

ps=pd.read_csv(OUT/'results/covid_patient_effect_summary.csv')
old=pd.read_csv(REPORT/'proportion_qualification/results/F_summary.csv')
fig,ax=plt.subplots(figsize=(178/25.4,90/25.4));positions=[];labels=[]
for i,(level,K) in enumerate([('major','k1000'),('major','k4000'),('fine','k1000'),('fine','k4000')]):
    y=3-i;labels.append(('Major' if level=='major' else 'Fine')+f' K = {int(K[1:]):,}');positions.append(y)
    a=old[(old.cohort=='covid')&(old.level==level)&(old.source==K)].iloc[0]
    b=ps[(ps.level==level)&(ps.source==K)&(ps.contrast=='infection')&(ps.wave=='all')].iloc[0]
    ax.errorbar(a.slope,y+.10,xerr=[[a.slope-a.slope_lo],[a.slope_hi-a.slope]],fmt='s',ms=4,color=gray,capsize=2,lw=.9,label='Sample bootstrap' if i==0 else None)
    ax.errorbar(b.slope,y-.10,xerr=[[b.slope-b.slope_lo],[b.slope_hi-b.slope]],fmt='D',ms=4,color=teal,capsize=2,lw=.9,label='Patient bootstrap' if i==0 else None)
ax.axvline(1,color=gray,linestyle='--',linewidth=.8);ax.set_yticks(positions,labels);ax.set_ylim(-.5,3.6);ax.set_xlim(0,1.1)
ax.set_xlabel('Generated-on-real CLR effect slope (1 = equal magnitude)');ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(.58,1.18),fontsize=8)
fig.subplots_adjust(left=.23,right=.97,bottom=.24,top=.75)
fig.text(.23,.95,'COVID effects remain attenuated after patient weighting',fontsize=11,weight='semibold')
fig.text(.23,.07,'Positive vs negative; existing transferred labels. Patient means within each state.\nPatient bootstrap preserves shared patients across states and sources; 2,000 replicates.',fontsize=8)
save(fig,'covid_patient_effects')

if (OUT/'results/covid_mapping_patient_summary.csv').exists():
    m=pd.read_csv(OUT/'results/covid_mapping_patient_summary.csv');fig,ax=plt.subplots(figsize=(178/25.4,100/25.4))
    pos=[];labs=[];i=0
    for route in ['R1','R2','R3']:
        for lv in ['major','fine']:
            y=5-i;pos.append(y);labs.append(f'{route} '+('Major' if lv=='major' else 'Fine'))
            for proto,off,color,mark in [('sample_overlap',.1,gray,'s'),('patient_disjoint',-.1,teal,'D')]:
                z=m[(m.source=='realT2')&(m.route==route)&(m.level==lv)&(m.protocol==proto)]
                if z.empty:continue
                r=z.iloc[0];ax.errorbar(r.patient_accuracy,y+off,xerr=[[r.patient_accuracy-r.accuracy_lo],[r.accuracy_hi-r.patient_accuracy]],fmt=mark,color=color,ms=4,capsize=2,lw=.9,label=proto.replace('_',' ') if i==0 else None)
            i+=1
    ax.set_yticks(pos,labs);ax.set_ylim(-.5,5.7);ax.set_xlim(.4,.95);ax.set_xlabel('Patient-averaged label accuracy (higher is better)');ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(.58,1.17),fontsize=8)
    fig.subplots_adjust(left=.20,right=.97,bottom=.27,top=.78)
    fig.text(.20,.95,'COVID mapping tested on unseen patients',fontsize=12,weight='semibold')
    fig.text(.20,.06,'PBISC-like real cells with known recap labels; 95% patient bootstrap intervals.\nR1 Harmony + kNN: fold 0. R2 PCA + kNN and R3 format-aware classifier: 2-fold OOF.\nSample-overlap protocol: 31 shared patients; patient protocol: zero shared patients.',fontsize=8)
    save(fig,'covid_patient_mapping')
print('Figures saved')
