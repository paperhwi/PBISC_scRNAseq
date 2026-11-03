"""Summarize three PBISC outputs against real A and held-out real B."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Malgun Gothic"
plt.rcParams["axes.unicode_minus"]=False

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
DIR=ROOT/"NewLora/analysis_20261102/three_way"
HIST=ROOT/"NewLora/analysis_20261102/historical_metric_comparison.csv"
GMAP={"adult_sepsis":"case","adult_HC":"control","adult_ICU_ctrl":"icu_control",
      "covid_positive":"case","covid_negative":"control"}
SRC={"realB":"실측 real B","prev4k":"기존 PBISC K4000","floor":"floor 수정","ceiling":"새 LoRA ceiling"}
ORDER=["prev4k","floor","ceiling","realB"]

def transfer(coh,source):
    if source=="ceiling":
        path=ROOT/"NewLora/analysis_20261102"/coh/"ceiling_cell_labels.csv.gz"
    else:
        old="realU" if source=="realB" else source
        path=BASE/"results/big"/f"{coh}_{old}_transfer.csv.gz"
    d=pd.read_csv(path,usecols=["sample","group","fine_tr"])
    d["condition"]=d.group.map(GMAP)
    return d[d.condition.notna()]

records=[]; detail=[]
for coh in ("sepsis","covid"):
    frames={source:transfer(coh,source) for source in ORDER}
    cats=sorted(set().union(*[set(x.fine_tr.astype(str)) for x in frames.values()]))
    groupmeans={}
    for source,frame in frames.items():
        counts=pd.crosstab(frame["sample"],frame["fine_tr"]).reindex(columns=cats,fill_value=0)
        rates=counts.div(counts.sum(1),axis=0)
        group=frame.drop_duplicates("sample").set_index("sample").condition.reindex(rates.index)
        groupmeans[source]={}
        for cond in (["case","control","icu_control"] if coh=="sepsis" else ["case","control"]):
            subset=rates[group==cond]
            if subset.empty:continue
            groupmeans[source][cond]=subset.mean()
            for cluster,value in subset.mean().items():
                detail.append({"cohort":coh,"condition":cond,"source":source,"fine_cluster":cluster,
                               "mean_sample_fraction":value,"n_samples":len(subset)})
    for source in ORDER:
        for cond in (["case","control","icu_control"] if coh=="sepsis" else ["case","control"]):
            if cond not in groupmeans[source] or cond not in groupmeans["realB"]:continue
            a=groupmeans["realB"][cond].to_numpy();b=groupmeans[source][cond].to_numpy();m=(a+b)/2
            with np.errstate(divide="ignore",invalid="ignore"):
                jsd=(np.where(a>0,a*np.log2(a/m),0).sum()+np.where(b>0,b*np.log2(b/m),0).sum())/2
            records.append({"cohort":coh,"source":source,"condition":cond,"n_realB_samples":frames["realB"].loc[frames["realB"].condition==cond,"sample"].nunique(),
                            "n_source_samples":frames[source].loc[frames[source].condition==cond,"sample"].nunique(),
                            "fine_TV_vs_realB":.5*np.abs(a-b).sum(),"fine_JSD_vs_realB":jsd,
                            "fine_fraction_r_vs_realB":np.corrcoef(a,b)[0,1]})

group=pd.DataFrame(records);group.to_csv(DIR/"group_matched_fine_similarity.csv",index=False)
pd.DataFrame(detail).to_csv(DIR/"group_fine_cluster_proportions.csv",index=False)

H=pd.read_csv(HIST).replace({"realU":"realB","new_ceiling":"ceiling"})
R=pd.read_csv(DIR/"three_way_summary.csv")
R=R[["cohort","source","n_fine_total","n_abundant_real","n_reproduced","real_A_mass_reproduced","source_mass_reproduced"]]
C=pd.read_csv(DIR/"common_cluster_transcriptome_summary.csv")
F=pd.read_csv(DIR/"fine_composition_summary.csv")
D=pd.read_csv(DIR/"cell_depth_summary.csv")
P=pd.read_csv(DIR/"common_cluster_paired_metrics.csv")
score=H.merge(R,on=["cohort","source"]).merge(C,on=["cohort","source"]).merge(F,on=["cohort","source"]).merge(D,on=["cohort","source"])
for cond in ("case","control"):
    t=group[group.condition==cond][["cohort","source","fine_TV_vs_realB","fine_fraction_r_vs_realB"]]
    t=t.rename(columns={"fine_TV_vs_realB":f"fine_TV_{cond}_vs_realB","fine_fraction_r_vs_realB":f"fine_r_{cond}_vs_realB"})
    score=score.merge(t,on=["cohort","source"])
score["source_order"]=score.source.map({s:i for i,s in enumerate(ORDER)})
score=score.sort_values(["cohort","source_order"]).drop(columns="source_order")
score.to_csv(DIR/"three_way_scorecard.csv",index=False)

fig,axs=plt.subplots(2,3,figsize=(13,7.3))
colors=["#8F9299","#4D82A8","#D45D48","#3F765A"]
for i,coh in enumerate(("sepsis","covid")):
    t=score[score.cohort==coh].set_index("source").loc[ORDER]
    labels=[SRC[s] for s in ORDER];x=np.arange(4)
    for j,(col,title) in enumerate((("iLISI_post","Harmony 후 iLISI (높을수록 혼합)"),
                                    ("n_reproduced","재현 fine cluster 수"),
                                    ("median_pearson_low_expr","공통 cluster 내 저발현 유전자 r"))):
        axs[i,j].bar(x,t[col],color=colors,width=.68)
        axs[i,j].set_xticks(x,labels,rotation=25,ha="right")
        axs[i,j].set_title(f"{coh.upper()} · {title}",fontsize=10,loc="left")
        axs[i,j].grid(axis="y",alpha=.15)
        axs[i,j].set_axisbelow(True)
        if col=="iLISI_post":axs[i,j].set_ylim(1,2)
        if col=="median_pearson_low_expr":axs[i,j].set_ylim(0,1)
fig.tight_layout();fig.savefig(DIR/"three_way_comparison.png",dpi=220,bbox_inches="tight");plt.close(fig)

def md(df):
    s=df.astype(str)
    return "\n".join(["| "+" | ".join(s.columns)+" |","| "+" | ".join(["---"]*len(s.columns))+" |"]+
                     ["| "+" | ".join(row)+" |" for row in s.itertuples(index=False,name=None)])

with (DIR/"비교_보고서.md").open("w",encoding="utf-8") as f:
    f.write("# 기존 PBISC · floor 수정 · 새 LoRA ceiling: real 대비 통합 및 cluster 전사체 비교\n\n")
    f.write("## 방법\n\n")
    f.write("기존 분석의 real A/B 환자 분할과 측정 유전자 패널을 그대로 사용했다. 각 출력에서 약 60,000개 세포를 샘플별 균형 추출했다. Sepsis는 성인, COVID는 새 데이터의 학습 지원 donor 중복 5개 샘플을 제외했다. 이전 출력의 저장된 동일 규모 표본을 사용했다.\n\n")
    f.write("**생성 조건 주의:** 새 파일명의 `K4000`은 샘플당 생성 세포 수 4,000을 가리킨다. 새 H5AD의 `generation_protocol.K_cond=1000`, `training_factor=1`, `inference_factor=1`, `readout=support-ceiling.v1`이다. Floor 생성본의 provenance는 `K_condition=4000`, bulk 입력 256배, 임계값 1e-4를 사용한 deterministic floor readout을 기록한다. Sepsis draw도 floor 0, 새 데이터 1로 다르다. 따라서 아래 차이를 LoRA 교체만의 효과로 분리할 수 없다.\n\n")
    f.write("혼합은 real A와 각 출력의 **별도 pairwise Harmony**에서 얻은 iLISI다. Fine cluster 라벨은 각 통합의 real A에서 30-NN으로 전달했다. 재현 cluster는 real A에서 ≥0.5%, 출력에서 ≥0.5%, 출력/real A 비율 0.5–2, 출력 세포 ≥100이라는 분석용 기술 기준을 모두 만족하는 군이다.\n\n")
    f.write("전사체는 각 cluster의 평균 CP10k를 log1p 변환하여 real A의 같은 fine cluster와 비교했다. 저발현 유전자는 real A의 평균 CP10k 순위 2,001위 이하 유전자다. `identity residual r`은 cluster 평균에서 각 source의 전체 평균 발현을 뺀 벡터를 HVG에서 비교하므로 세포형 특이 발현 패턴을 본다. 전사체 비교의 요약값은 세 생성 출력이 **모두 재현한 같은 cluster**(Sepsis 15개, COVID 13개)의 중앙값이다. Real B는 실측 양성 대조군이다.\n\n")
    f.write("## 전체 비교\n\n")
    cols=["cohort","source","iLISI_post","n_reproduced","real_A_mass_reproduced","fine_TV_case_vs_realB","fine_TV_control_vs_realB",
          "median_pearson_all","median_pearson_low_expr","median_identity_residual_r","median_rmse_log1p_cp10k"]
    t=score[cols].copy();t.source=t.source.map(SRC)
    f.write(md(t.round(3))+"\n\n")
    f.write("`real_A_mass_reproduced`는 재현 기준을 통과한 cluster들이 real A 세포 중 차지하는 비율이다. `fine_TV_case/control_vs_realB`는 각 환자/대조군의 **샘플별 fine cluster 비율 평균**을 실측 real B와 비교한 total variation이다(낮을수록 유사).\n\n")
    for coh in ("sepsis","covid"):
        t=score[score.cohort==coh].set_index("source")
        f.write(f"## {coh.upper()}\n\n")
        f.write(f"- Harmony 후 iLISI: 기존 {t.loc['prev4k','iLISI_post']:.3f}, floor {t.loc['floor','iLISI_post']:.3f}, ceiling {t.loc['ceiling','iLISI_post']:.3f}, real B {t.loc['realB','iLISI_post']:.3f}.\n")
        f.write(f"- 재현 fine cluster: 기존 {int(t.loc['prev4k','n_reproduced'])}, floor {int(t.loc['floor','n_reproduced'])}, ceiling {int(t.loc['ceiling','n_reproduced'])}; abundant real A cluster {int(t.loc['ceiling','n_abundant_real'])}개 중.\n")
        f.write(f"- 공통 cluster 내 전체 유전자 r 중앙값: 기존 {t.loc['prev4k','median_pearson_all']:.3f}, floor {t.loc['floor','median_pearson_all']:.3f}, ceiling {t.loc['ceiling','median_pearson_all']:.3f}, real B {t.loc['realB','median_pearson_all']:.3f}.\n")
        f.write(f"- 저발현 유전자 r 중앙값: 기존 {t.loc['prev4k','median_pearson_low_expr']:.3f}, floor {t.loc['floor','median_pearson_low_expr']:.3f}, ceiling {t.loc['ceiling','median_pearson_low_expr']:.3f}, real B {t.loc['realB','median_pearson_low_expr']:.3f}.\n")
        f.write(f"- Cluster 특이 발현 residual r 중앙값: 기존 {t.loc['prev4k','median_identity_residual_r']:.3f}, floor {t.loc['floor','median_identity_residual_r']:.3f}, ceiling {t.loc['ceiling','median_identity_residual_r']:.3f}, real B {t.loc['realB','median_identity_residual_r']:.3f}.\n")
        f.write(f"- 환자군 fine 구성 TV vs real B: 기존 {t.loc['prev4k','fine_TV_case_vs_realB']:.3f}, floor {t.loc['floor','fine_TV_case_vs_realB']:.3f}, ceiling {t.loc['ceiling','fine_TV_case_vs_realB']:.3f}.\n\n")
        n=int(t.loc['ceiling','n_common_three'])
        low=P[(P.cohort==coh)&(P.metric=="pearson_low_expr")]
        identity=P[(P.cohort==coh)&(P.metric=="identity_residual_r")]
        f.write(f"- 같은 {n}개 cluster를 짝지으면 저발현 유전자 r에서 floor가 ceiling보다 높은 cluster는 {int((low.ceiling_minus_floor<0).sum())}/{n}개다. 반대로 cluster 특이 residual r에서 ceiling이 floor보다 높은 cluster는 {int((identity.ceiling_minus_floor>0).sum())}/{n}개다.\n")
        f.write(f"- 세포당 검출 유전자 중앙값: 기존 {t.loc['prev4k','median_detected_genes']:.0f}, floor {t.loc['floor','median_detected_genes']:.0f}, ceiling {t.loc['ceiling','median_detected_genes']:.0f}, real B {t.loc['realB','median_detected_genes']:.0f}.\n\n")
    f.write("## 해석\n\n")
    f.write("Floor 수정은 두 cohort에서 재현 cluster 내 **절대 발현 프로필**(전체/저발현 유전자 상관 및 log-CP10k RMSE)이 가장 real에 가깝다. 새 ceiling은 기존 K4000보다 발현 유사도가 높지만 floor에는 못 미친다. 반면 ceiling의 cluster 특이 발현 residual r은 floor보다 높다. 따라서 ceiling이 세포형별 상대 패턴을 더 선명하게 유지하더라도 발현량과 저발현 유전자 분포에는 더 큰 오차가 남을 수 있다.\n\n")
    f.write("COVID 환자군의 real B 대비 fine 구성 TV는 기존 0.179, floor 0.225, ceiling 0.259다. 즉 COVID에서 기존 출력은 구성 비율이 더 가깝지만, 혼합과 cluster 내부 발현 유사도는 크게 낮다. 이 세 축은 별도로 판단해야 한다.\n\n")
    f.write("COVID에서 ceiling의 real A 기준 fine 비율 TV는 좋아 보일 수 있다. 그러나 real B와 **질환군별**로 비교하면 T cell 과다·Monocyte 부족이 드러난다. 기준 real A와 real B 사이에도 환자 구성 차이가 있으므로 두 기준을 함께 제시했다.\n\n")
    f.write("## 제한\n\n")
    f.write("- 새 파일에는 세포형 라벨이 없어, 모든 생성 출력의 fine 비교에 발현 기반 30-NN 전달 라벨을 사용했다. 같은 발현 정보를 라벨 선정과 전사체 비교에 사용하므로 유사도가 낙관적일 수 있다.\n")
    f.write("- 세 출력의 생성 설정과 draw가 다르며, 신뢰구간이나 반복 draw 검정 없이 수치 순위를 통계적 우열로 해석하지 않는다.\n")
    f.write("- Sepsis real B의 건강 대조군은 4개, ICU 대조군은 2개 샘플이다. 샘플/환자가 추론 단위이며 생성 세포는 독립 생물학적 반복이 아니다.\n")
    f.write("- iLISI는 혼합 정도, 구성 TV는 cluster 수량, 전사체 r/RMSE는 cluster 내부 발현을 각각 측정한다. 어느 한 값도 전체 생성 품질을 대표하지 않는다.\n\n")
    f.write("## 파일\n\n")
    f.write("- `three_way_scorecard.csv`: 핵심 지표 종합\n- `fine_cluster_reproduction.csv`: cluster별 비율·재현 판정\n- `per_cluster_transcriptome_similarity.csv`: cluster별 발현 상관·RMSE·marker Jaccard\n- `common_cluster_transcriptome_summary.csv`, `common_cluster_paired_metrics.csv`: 동일 재현 cluster의 중앙값과 cluster별 짝비교\n- `group_matched_fine_similarity.csv`, `group_fine_cluster_proportions.csv`: real B 대비 조건별 fine 구성\n- `cell_depth_summary.csv`, `three_way_comparison.png`: readout 깊이와 그래프\n- 재현 스크립트: `../../compare_three_transcriptomes.py`, `../../summarize_three_way.py`\n")
print(DIR/"비교_보고서.md")
