"""Render Figure7 only from package-relative, immutable canonical result tables.

No statistical model is fit here. Pseudonymized patient points and all displayed
statistics are copied to source_data/Figure_7_*.csv for direct inspection.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from plot_style import ROOT,DATA,FIGURES,MM,COLORS,MARKERS,COHORT_ORDER,setup_style,label_panel,save_figure
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

SOURCE=ROOT/'source_data'

def read(rel):return pd.read_csv(DATA/rel)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def fmtq(x):return f'{x:.4f}'

def row(data,arm,feature):
    q=data[data.arm.eq(arm)&data.feature.eq(feature)]
    assert len(q)==1,(arm,feature)
    return q.iloc[0]

def forest(ax,intervals,feature,metric='rho',ys=(2,1,0),labels=True):
    lo,hi=('rho_CI_low','rho_CI_high') if metric=='rho' else ('AUC_CI_low','AUC_CI_high')
    ax.axvline(0 if metric=='rho' else .5,color='#B8BEBA',lw=.7,ls=(0,(3,3)),zorder=0)
    records=[]
    for arm,y in zip(COHORT_ORDER,ys):
        t=row(intervals,arm,feature)
        assert np.isfinite([t[metric],t[lo],t[hi]]).all()
        ax.errorbar(t[metric],y,xerr=[[t[metric]-t[lo]],[t[hi]-t[metric]]],
                    fmt=MARKERS[arm],ms=5 if arm=='Combined' else 4.7,color=COLORS[arm],
                    mfc=COLORS[arm],mec='white',mew=.45,elinewidth=1.05 if arm=='Combined' else .85,
                    capsize=2,capthick=.8,zorder=3)
        records.append(dict(arm=arm,feature=feature,metric=metric,y=y,effect=t[metric],CI_low=t[lo],CI_high=t[hi],
                            B=int(t.B),seed=int(t.seed),fixed_score=bool(t.fixed_score),
                            bootstrap_strata=t.rho_strata if metric=='rho' else t.AUC_strata))
    if labels:
        ax.set_yticks(ys,COHORT_ORDER)
        ax.tick_params(axis='y',length=0,pad=5)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#8B918E')
    if metric=='rho':
        ax.set_xlim(-.5,1);ax.set_xticks([-.5,0,.5,1]);ax.set_xlabel('Spearman ρ (95% CI)',labelpad=4)
    else:
        ax.set_xlim(-.02,1.02);ax.set_xticks([0,.5,1]);ax.set_xlabel('AUC for death (95% CI)',labelpad=4)
    return records

def text_bounds(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height()
    outside=[]
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text():continue
        bb=t.get_window_extent(renderer)
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>w+.5 or bb.y1>h+.5:
            outside.append(dict(text=t.get_text(),bounds=list(bb.extents)))
    return outside

def main():
    setup_style();SOURCE.mkdir(exist_ok=True)
    required=['plotting/composition_patients.csv','plotting/composition_scores.csv',
              'composition/all43_associations.csv','composition/all43_fixed_score_bootstrap_CI.csv',
              'composition/all43_source_stratified_permutation.csv','composition/cohort_summary.csv']
    for rel in required:assert (DATA/rel).exists(),rel
    patients=read(required[0]);scores=read(required[1]);assoc=read(required[2]);ci=read(required[3]);perm=read(required[4]);coh=read(required[5])
    assert patients.patient_id.nunique()==66 and len(patients)==66
    assert patients.groupby('arm').size().to_dict()=={'New':29,'Real':37}
    assert int(patients.death.sum())==13 and int(patients.draws.sum())==217
    assert patients.patient_id.tolist()==scores.patient_id.tolist()
    assert ci.B.eq(50000).all()
    for feature in ['HSPC','B_ASC_IgG','B_ASC_IgA']:
        for arm in COHORT_ORDER:
            aa=row(assoc,arm,feature);bb=row(ci,arm,feature)
            assert abs(aa.rho-bb.rho)<1e-12 and abs(aa.AUC-bb.AUC)<1e-12
    fig=plt.figure(figsize=(178*MM,210*MM))
    # Fixed geometry is defined at manuscript size, rather than resized later.
    a=fig.add_axes([.06,.835,.90,.115])
    b=fig.add_axes([.115,.50,.35,.245])
    c=fig.add_axes([.65,.50,.315,.245])
    d=fig.add_axes([.115,.11,.35,.245])
    e=fig.add_axes([.65,.11,.315,.245])
    label_panel(a,'A','Actual bulk-only patients extend the cohort',x=-.025,y=1.08)
    a.set_axis_off()
    positions=[.14,.48,.83]
    for arm,x in zip(COHORT_ORDER,positions):
        t=coh[coh.arm.eq(arm)].iloc[0]
        a.text(x,.84,{'Real':'Real','New':'New bulk-only','Combined':'Combined'}[arm],ha='center',color=COLORS[arm],fontsize=9,fontweight='semibold')
        a.text(x,.39,f'{int(t.patients)}',ha='center',color=COLORS[arm],fontsize=20,fontweight='semibold')
        a.text(x,.17,'patients',ha='center',fontsize=8.5)
        a.text(x,-.07,f'{int(t.deaths)} deaths',ha='center',fontsize=8.5)
    for x,symbol in [(.30,'+'),(.65,'=')]:a.text(x,.40,symbol,ha='center',fontsize=16,color='#5E6462')
    coh.to_csv(SOURCE/'Figure_7_A_cohort.csv',index=False)
    # Patient-level, untransformed composition. No fitted pooled trend line.
    label_panel(b,'B','HSPC abundance across severity',x=-.075,y=1.05)
    rng=np.random.default_rng(71118);raw=[]
    for arm in ['Real','New']:
        t=patients[patients.arm.eq(arm)].copy()
        offset=-.14 if arm=='Real' else .14
        xx=t.severity.to_numpy()+offset+rng.uniform(-.075,.075,len(t))
        b.scatter(xx,t.HSPC,s=17 if arm=='Real' else 22,marker=MARKERS[arm],c=COLORS[arm],alpha=.85,
                  edgecolors='white',linewidths=.35,zorder=3,label=arm)
        for (_,p),x in zip(t.iterrows(),xx):raw.append(dict(patient_id=p.patient_id,arm=arm,severity=int(p.severity),death=int(p.death),HSPC_pct=p.HSPC,display_x=x))
    b.set_xlim(-.5,3.5);b.set_ylim(-.08,3.75)
    b.set_xticks([0,1,2,3],['Mild','Moderate','Severe','Critical'])
    b.set_yticks([0,1,2,3]);b.set_ylabel('HSPCs (% of leukocytes)',labelpad=4)
    b.set_xlabel('Peak WHO severity',labelpad=4)
    b.legend(loc='upper left',ncol=2,handletextpad=.35,columnspacing=1,frameon=False,borderaxespad=0)
    b.text(0,-.30,'Raw patient means; tests use source INT',transform=b.transAxes,fontsize=8,ha='left',va='top')
    pd.DataFrame(raw).to_csv(SOURCE/'Figure_7_B_HSPC_raw_patients.csv',index=False)
    # Source effects and uncertainty stay visible, including the weak Real effect.
    label_panel(c,'C','HSPC severity signal emerges',x=-.085,y=1.05)
    cr=forest(c,ci,'HSPC');c.set_ylim(-.6,2.65)
    rr=row(assoc,'Real','HSPC');cc=row(assoc,'Combined','HSPC');pp=row(perm,'Combined','HSPC')
    c.text(0,-.30,f'Real P={rr.p_severity:.3f}; combined q43={cc.q43_severity:.4f}\nSource-permutation q43={pp.q43_perm_severity:.4f}',
           transform=c.transAxes,fontsize=8,ha='left',va='top',linespacing=1.4)
    for z in cr:
        q=row(assoc,z['arm'],'HSPC');z.update(n=int(q.n),deaths=int(q.deaths),P=q.p_severity,q43=q.q43_severity)
    pd.DataFrame(cr).to_csv(SOURCE/'Figure_7_C_HSPC_severity.csv',index=False)
    # AUC direction is never flipped; Combined is the union, not another cohort.
    label_panel(d,'D','Consistent HSPC death association',x=-.075,y=1.05)
    dr=forest(d,ci,'HSPC',metric='AUC');d.set_ylim(-.6,2.8)
    d.text(.02,.98,f'Combined q43={cc.q43_death:.4f}',transform=d.transAxes,fontsize=8,va='top',
           bbox=dict(facecolor='white',edgecolor='none',pad=1))
    for z in dr:
        q=row(assoc,z['arm'],'HSPC');z.update(n=int(q.n),deaths=int(q.deaths),P=q.p_death,q43=q.q43_death)
    pd.DataFrame(dr).to_csv(SOURCE/'Figure_7_D_HSPC_death.csv',index=False)
    # Aligned ASC small groups use the same rho axis as C.
    label_panel(e,'E','ASC severity evidence strengthens',x=-.085,y=1.05)
    er=[]
    for feature,ys,title,hy in [('B_ASC_IgG',(5.8,4.8,3.8),'IgG ASC',6.6),('B_ASC_IgA',(1.8,.8,-.2),'IgA ASC',2.6)]:
        er.extend(forest(e,ci,feature,ys=ys,labels=False))
        q=row(assoc,'Combined',feature).q43_severity
        e.text(0,hy,f'{title} · combined q43={q:.4f}',fontsize=8.5,va='center',
               transform=e.get_yaxis_transform(),bbox=dict(facecolor='white',edgecolor='none',pad=.5))
    e.set_yticks([5.8,4.8,3.8,1.8,.8,-.2],COHORT_ORDER*2);e.tick_params(axis='y',length=0,pad=5)
    e.set_ylim(-.9,7.15)
    for z in er:
        q=row(assoc,z['arm'],z['feature']);z.update(n=int(q.n),deaths=int(q.deaths),P=q.p_severity,q43=q.q43_severity)
    pd.DataFrame(er).to_csv(SOURCE/'Figure_7_E_ASC_severity.csv',index=False)
    # Full precision source tables; rows are selected only by the approved figure.
    selected=['HSPC','B_ASC_IgG','B_ASC_IgA']
    perm[perm.feature.isin(selected)].to_csv(SOURCE/'Figure_7_source_permutation.csv',index=False)
    scores[['patient_id','arm','severity','death']+selected].to_csv(SOURCE/'Figure_7_patient_INT_scores.csv',index=False)
    fig.text(.115,.026,'95% fixed-score bootstrap CIs (50,000 draws); q43: BH across 43 populations per endpoint',fontsize=8,ha='left')
    outside=text_bounds(fig)
    if outside:raise ValueError('Text outside canvas: '+json.dumps(outside))
    spec=json.loads((ROOT/'figure7_spec.json').read_text(encoding='utf-8'))
    metadata=dict(spec='figure7_spec.json',input_sha256={rel:digest(DATA/rel) for rel in required},
                  no_statistical_refitting=True,point_n=66,CI_B=50000,
                  raw_scale='linear leukocyte percentage; jitter only affects x display',text_outside_canvas=outside,
                  notes=['Combined patients overlap the two source arms; these are not independent replicates.',
                         'Strict/FDR threshold transitions are not formal augmentation benefit tests.',
                         'AUC is an univariable association, not severity-independent prognosis.'])
    record=save_figure(fig,'Figure_7',metadata=metadata)
    # Fixed-size QA derivatives are kept separate from publication exports.
    from PIL import Image,ImageOps
    qa=ROOT/'qa';qa.mkdir(exist_ok=True)
    im=Image.open(FIGURES/'Figure_7.png').convert('RGB')
    preview=im.resize((round(178/25.4*120),round(210/25.4*120)),Image.Resampling.LANCZOS)
    preview.save(qa/'Figure_7_finalsize_120dpi.png')
    ImageOps.grayscale(preview).save(qa/'Figure_7_grayscale_120dpi.png')
    arr=np.asarray(preview,dtype=float)/255
    # Machado deuteranopia simulation, sufficient here with marker/position redundancy.
    matrix=np.array([[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.01182,.04294,.968881]])
    sim=np.clip(arr@matrix.T,0,1)
    Image.fromarray(np.uint8(np.round(sim*255))).save(qa/'Figure_7_deuteranopia_120dpi.png')
    print(json.dumps({'figure':'Figure_7','size_mm':[record['width_mm'],record['height_mm']],
                      'minimum_font_pt':min(t['fontsize_pt'] for t in record['texts']),'source_data_files':len(list(SOURCE.glob('Figure_7_*.csv')))},indent=2))

if __name__=='__main__':main()
