"""Add consistently defined sign and top-k metrics to corrected DEG analysis."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import corrected_deg as d

out=Path('F:/PBISC/project_review_20261028/corrected_deg')
summary=pd.read_csv(out/'corrected_summary.csv')
rows=[]
for cohort,roster in [('sepsis','all'),('covid','published_187'),('covid','all_local_198')]:
    genes,measured,mats,_=d.prepare(cohort,roster)
    names=genes[measured]
    for ct in d.TYPES:
        gR=mats[ct]['real']; gF=mats[ct]['floor']
        yr=gR[list(genes)].to_numpy(float); yf=gF[list(genes)].to_numpy(float)
        cr=(gR._c=='case').to_numpy(int); cf=(gF._c=='case').to_numpy(int)
        test=d.patient_cluster_deg if cohort=='covid' else d.welch_deg
        R=test(yr,cr,gR._u.astype(str).to_numpy()) if cohort=='covid' else test(yr,cr)
        F=test(yf,cf,gF._u.astype(str).to_numpy()) if cohort=='covid' else test(yf,cf)
        lr,rd,qr=(x[measured] for x in R)
        lf,fd,qf=(x[measured] for x in F)
        rank=stats.rankdata(-(yr[cr==1].mean(0)[measured]+yr[cr==0].mean(0)[measured])/2)
        low=rank>2000
        real_set=set(names[rd])
        def top(lfc,deg):
            ix=np.flatnonzero(deg)
            ix=ix[np.argsort(-abs(lfc[ix]),kind='stable')[:100]]
            return list(names[ix])
        tr=set(top(lr,rd)); tf=top(lf,fd)
        rows.append(dict(cohort=cohort,roster=roster,celltype=ct,
            n_genes=len(names),n_rank2000plus=int(low.sum()),
            sign_all=float((np.sign(lr)==np.sign(lf)).mean()),
            sign_real_deg=float((np.sign(lr[rd])==np.sign(lf[rd])).mean()) if rd.sum() else np.nan,
            sign_rank2000plus=float((np.sign(lr[low])==np.sign(lf[low])).mean()) if low.sum() else np.nan,
            real_top100_n=len(tr),gen_top100_n=len(tf),
            top100_overlap_n=len(tr&set(tf)),
            gen_top100_precision=float(sum(g in real_set for g in tf)/len(tf)) if tf else np.nan))

a=pd.DataFrame(rows)
joined=a.merge(summary,on=['cohort','roster','celltype'])
assert len(joined)==15
assert np.allclose(joined.sign_all_x,joined.sign_all_y,equal_nan=True)
assert np.allclose(joined.sign_real_deg_x,joined.sign_real_deg_y,equal_nan=True)
assert (joined.n_genes_x==joined.n_genes_y).all()
a.to_csv(out/'corrected_gene_metrics.csv',index=False)

lines=[
    '# DEG·부호·top-100 지표 정의 — 2026-10-28',
    '',
    '수정 DEG 주분석과 동일한 유전자·환자 명단·검정으로 계산했다. `corrected_gene_metrics.csv`의 15개 행은 `corrected_summary.csv`와 1:1로 일치하며 sign 수치도 자동 교차 검증했다.',
    '',
    '- **공통 유전자:** 실제와 생성 전체 pseudobulk에서 모두 nonzero인 유전자. DEG FDR은 이 필터 전에 canonical 전체 유전자에서 BH 보정한다.',
    '- **DEG:** q<0.05 및 |log2FC|≥1을 동시에 만족. Sepsis는 Welch, COVID는 환자 cluster CR1. 실제와 생성에 같은 규칙을 적용한다.',
    '- **sign_all:** 공통 유전자 전체의 logFC 부호 일치 수 / 공통 유전자 수.',
    '- **sign_real_deg:** 실제 DEG에 한정한 부호 일치 수 / 실제 DEG 수. 실제 DEG가 없으면 결측.',
    '- **sign_rank2000plus:** 실제 case·control 평균 발현으로 내림차순 순위를 매긴 뒤 2000위보다 낮은 공통 유전자의 부호 일치 수 / 해당 유전자 수. 이 값은 `sign_all`이나 `sign_real_deg`와 다른 지표다.',
    '- **real/gen top100:** 각각의 DEG 집합을 |log2FC| 내림차순으로 정렬해 최대 100개 선택. `top100_overlap_n`은 두 목록의 교집합 **개수**다.',
    '- **gen_top100_precision:** 생성 top100 유전자 중 실제 DEG인 수 / 생성 top100의 실제 길이. 생성 DEG가 100개 미만이면 그 개수를 분모로, 0개면 결측으로 둔다.',
    '',
    '이전 10/27 `run_deg.py`의 top100 목록은 q<0.05만으로 고르고 |log2FC|≥1을 강제하지 않았다. 따라서 새 top100 값과 이전 표·슬라이드 수치를 직접 교체해 변화량으로 해석하지 않는다. 원고에는 위 새 정의와 해당 CSV의 수치를 함께 사용한다.',
    '',
    '| 코호트 | 명단 | 세포형 | sign all | sign real DEG | sign rank>2000 | gen top N | gen top precision | top overlap N |',
    '|---|---|---|---:|---:|---:|---:|---:|---:|',
]
for _,r in a.iterrows():
    fmt=lambda x:'—' if pd.isna(x) else f'{x:.3f}'
    lines.append(f'| {r.cohort} | {r.roster} | {r.celltype} | {fmt(r.sign_all)} | {fmt(r.sign_real_deg)} | {fmt(r.sign_rank2000plus)} | {r.gen_top100_n} | {fmt(r.gen_top100_precision)} | {r.top100_overlap_n} |')
(out/'METRIC_DEFINITIONS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('wrote',out/'corrected_gene_metrics.csv',out/'METRIC_DEFINITIONS.md')
