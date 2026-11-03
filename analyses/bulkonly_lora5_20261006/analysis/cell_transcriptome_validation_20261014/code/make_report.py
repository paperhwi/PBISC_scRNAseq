"""Report actual execution separately from prepared, unexecuted experiments."""
import json,hashlib,csv
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]
real=json.load(open(OUT/'results/real_counts_manifest.json'))
qa=json.load(open(OUT/'VALIDATION.json'))
d=list(csv.DictReader(open(OUT/'results/reference_metric_control_deltas.csv')))
def delta(cohort,metric):
    r=next(x for x in d if x['cohort']==cohort and x['control']=='independent_gene_shuffle' and x['metric']==metric)
    return f"{float(r['delta']):.4f} [{float(r['ci_lo']):.4f}, {float(r['ci_hi']):.4f}]"
report=f'''# PBISC 세포 전사체 재현 검증 — 실행 결과와 남은 입력

2026-10-14. 요청한 순서: ① held-out paired counts → ② 유형 내 공동 발현·분포와 bulk 교환 → ③ 희귀 상태 제거·추가 → ④ 독립 seed·외부 코호트.

**현재 완료 범위는 실제 held-out counts 확보와 검증용 입력·대조군 점검이다. paired 생성 counts와 모델 실행환경이 로컬에 없어 PBISC에 대한 ②–④의 새 생성 실험은 실행하지 못했다. 아래 준비 결과를 PBISC 가설의 지지 결과로 해석하면 안 된다.**

| 단계 | 실제 수행 | 아직 수행하지 못한 작업 |
|---|---|---|
| ① | 다섯 LoRA 모두에서 제외된 환자 명단 고정, 실제 counts 686,543세포 확보, 원본 일치 점검 | 동일 환자 paired 생성 counts 확보, base-model 사전학습 환자 중복 확인 |
| ② | 실제 세포 352개 sample×major-type 그룹에서 지표 대조군 실행, bulk 교환 배정 471건 생성 | PBISC 유형 내 공동 발현·분포 비교, 교환 bulk로 새 생성 및 비교 |
| ③ | 판별 가능 후보 선별, MAIT/NK IFN/CD14 IFN 제거·추가 입력 105개와 유한 세포 수 검출력 계산 | canonical bulk condition 재구성, PBISC 재생성, 부재 시 위양성·dose response·de novo cluster 확인 |
| ④ | 20개 생성 파일의 체크포인트·seed 출처 감사 | 고정 체크포인트 독립 seed 반복, 독립 외부 코호트 생성·검증 |

## ① 실제 counts 확보 결과

| 코호트 | 환자 | 샘플 | 세포 | 일대일 대응하는 측정 패널 유전자 |
|---|---:|---:|---:|---:|
| Sepsis | {real[0]['patients']} | {real[0]['samples']} | {real[0]['n_cells']:,} | {real[0]['n_genes']:,} |
| COVID/ESKD | {real[1]['patients']} | {real[1]['samples']} | {real[1]['n_cells']:,} | {real[1]['n_genes']:,} |

`inputs/ye2026_sepsis_real_heldout_panel_counts.h5ad`, `inputs/imperial_eskd_real_heldout_panel_counts.h5ad`의 X는 정수값 raw count이며, source_row, original_sample_id, sample_id, patient_id, major/fine recap labels와 패널 유전자를 저장했다. 정규화 X나 logcounts를 counts로 사용하지 않았다. 원래 전체 전사체 counts는 기존 original_scRNA h5ad에 보존돼 있다. 이 패널 기반 검증만으로 패널 밖 전사체까지 재현했다고 주장할 수 없다.

Sepsis는 ENSG, COVID는 유전자 symbol로 일대일 대응했다. 측정되지 않거나 중복 대응되는 유전자는 제외했다. 원본 `Skin-152`와 benchmark `Skin-S152`는 명시적 별칭으로 연결했고 원래 이름을 보존했다. 서버의 paired row manifest로 이 대응을 마지막으로 대조해야 한다. Sepsis에는 소아가 포함되므로 adult-only 분석을 따로 보고해야 한다.

다섯 LoRA 학습 환자와 검증 환자의 중복은 0이다. **이 제외는 LoRA 단계에 관한 것이며 base model 사전학습과의 독립성은 아직 확인하지 못했다.** 원본에서 코호트별 12개 행을 다시 읽어 모든 보존 유전자의 count가 정확히 일치함을 확인했다. SHA256은 `results/real_counts_manifest.json`에 기록했다.

로컬 두 분석 압축파일의 5,157개 member를 조사했다. 행렬 후보 3개는 PCA 좌표이며 raw counts가 아니다. 로컬 생성 h5ad/NPZ는 bulk-only 샘플용이다. paired feature summaries, marker summaries, pseudobulk sums 또는 PCA 좌표를 full paired 세포 counts로 대체하지 않았다. 새 체크포인트 파일도 로컬에 없다.

## ② 실제 세포에서 공동 발현·분포 측정법 점검

PBISC 결과가 아닌 **실제 세포만 사용한 지표 점검**이다. 환자별 첫 샘플의 큰 major type 최대 4개를 선택했다. 그룹당 160세포 이상, 최대 300세포를 반으로 나누었다. real half A에서 분산 상위 50유전자를 선택했고 패널 기준 log1p CP10,000을 사용했다. Sepsis 60명·233그룹, COVID 33명·119그룹이 조건을 만족했다. 모든 환자를 분석 대상으로 확보했지만 이 점검에는 충분한 세포가 없는 COVID 2명이 제외됐다.

대조군은 real split-half, half B의 유전자별 독립 shuffle, half B의 centroid collapse다. shuffle은 유전자별 분포를 정확히 보존하며 공동 발현만 깨뜨린다. 이 shuffle 대조군은 **bulk 교환 대조군과 다른 실험**이다.

| real split-half 대비 gene shuffle의 추가 오차 | Sepsis: 환자 평균 [95% CI] | COVID: 환자 평균 [95% CI] |
|---|---:|---:|
| 유전자별 분포 W1 | 0 [0, 0] | 0 [0, 0] |
| 공동 발현 correlation MAE | {delta('ye2026_sepsis','corr_MAE')} | {delta('imperial_eskd','corr_MAE')} |
| 표준화 covariance MAE | {delta('ye2026_sepsis','covariance_MAE')} | {delta('imperial_eskd','covariance_MAE')} |
| 다변량 energy distance | {delta('ye2026_sepsis','energy_distance')} | {delta('imperial_eskd','energy_distance')} |

sample×type를 환자 내 평균한 뒤 환자를 같은 비중으로 취급했다. CI는 환자 bootstrap 2,000회다. 단변량 분포가 동일해도 공동 발현·다변량 지표가 악화를 포착했다. Centroid collapse도 분포·energy 지표에서 크게 악화했다. 모든 경우를 한 지표만으로 판단하면 안 된다. 이 half-A 특성 선택은 지표 점검용이며 PBISC의 확증적 검증에서 고정할 gene set/HVG 선택을 대신하지 않는다.

157 held-out 샘플마다 3회, 총 471건 bulk 교환 배정을 저장했다. 모든 교환은 다른 환자의 bulk를 사용한다. Sepsis는 group, COVID는 group×wave×severity_peak를 맞추었다. COVID 방문 시점까지 맞춘 민감도 분석은 별도 필요하다. 이 배정에 따른 재생성은 아직 없다.

생성 counts가 확보되면 ①을 마친 뒤 같은 환자·같은 유형의 real/PBISC를 동일 N·유전자·depth 조건에서 비교한다. detection, mean–variance, tail, marginal distance와 covariance/coexpression, multivariate distance/coverage를 같이 평가하고 Harmony로 source gap을 지우지 않는다. real split-half, readout-matched real, reference resampling, type-mean+noise 대조군을 함께 둔다. Bulk 교환에는 같은 체크포인트·noise/readout seed를 사용하여 native bulk의 추가 정보를 환자 단위로 검정한다. 최종 평가 margin과 marker/gene 목록은 새 생성 결과를 보기 전에 고정해야 한다.

## ③ 희귀 상태 제거·추가 입력과 검출력

기존 patient-disjoint R3 realT2 mapping의 recall/precision≥0.5, n_true_cells≥25, 현재 held-out 실제 빈도<2%를 **준비 단계의 후보 기준**으로 사용했다. 15개 후보 중 MAIT, NK IFN stim, CD14mono_IFN을 먼저 준비했다. 이것은 생성 세포의 정답률이나 독립 marker identity 검증 결과가 아니다. IFN 상태는 반드시 해당 lineage의 독립 marker와 함께 판정해야 한다. DC3, DC3 IFN, gdT, dex-monocyte 등 식별이 약한 상태를 1차 성공 판정에 쓰지 않는다.

각 상태를 가진 환자 5명씩, 합계 11명의 15개 sample×state에서 0%, 0.05%, 0.1%, 0.2%, 0.5%, 1%, 2% dose를 만들었다. 총 105개 입력이다. 희귀 상태를 같은 parent lineage의 state-negative 세포로 교체해 parent 비율을 고정했다. panel CPM 총합은 1,000,000이며 per-cell RNA 기여를 같게 맞추었다. 환자별 20개 이상의 실제 state cell을 가진 샘플을 우선 선택했으므로 희귀 상태가 거의 없는 환자에게 결과를 일반화하면 안 된다.

**현재 NPZ는 5,745유전자 패널의 pseudobulk 입력이다. 모델의 20,097 canonical-gene condition이 아니므로 그대로 sampler에 넣으면 안 된다.** 서버의 canonical gene table과 원본 전체 유전자를 사용해 동일 dose의 full-gene pseudobulk를 재구성하고 원래 condition builder로 변환해야 한다. 생성 세포의 제거·추가 반응은 아직 검정하지 않았다.

진짜 빈도가 0.1%일 때 생성 세포가 하나도 없을 확률은 N=200에서 81.9%, N=1,500에서 22.3%, N=4,000에서 1.83%다. 200-cell subsample만으로 희귀 상태 부재를 입증할 수 없다. 성공 기준을 10세포 이상 회수로 두면 N=30,000에서 기대 state cell은 30개다. 이 계산은 완벽하게 빈도를 재현하는 독립 표본 생성이라는 가정 아래의 binomial 표본 수 계산이며, PBISC의 실제 검출력 측정이 아니다.

서버 실행 시에는 0 dose의 위양성, dose 증가에 따른 회수 빈도·marker identity·within-lineage de novo cluster 안정성을 환자 단위로 평가한다. 라벨 transfer로 희귀 cluster를 강제로 만들지 않는다. 분류기 불확실성, seed, 총 N, depth와 adapter 차이를 분리해서 보고한다.

## ④ 독립 seed와 외부 코호트

각 코호트×K의 5개 draw는 체크포인트 SHA가 모두 다르다. 현재 파일에서는 한 체크포인트를 여러 독립 sampler seed로 반복한 실험이 없다. K1의 공통 sampler seed와 K4의 draw별 seed 차이는 기존 README에 문서화돼 있으나 다른 adapter 효과와 분리할 수 없다. 따라서 기존 5-draw 일치를 독립 seed 재현성으로 보고하지 않는다.

②–③ 다음에 고정 체크포인트별 sampler seed `2026101401`, `2026101402`, `2026101403`을 사용하고 Poisson readout 반복을 따로 둔다. K 비교의 1차 N은 1,500으로 같게, 희귀 상태 dose 실험은 N=30,000으로 같게 준비했다. Recipe K와 생성 세포 수 N은 다른 변수다.

외부 코호트에서 새 생성·검증은 실행하지 않았다. 기존 sepsis/COVID 데이터 둘을 다시 쓰는 것은 외부 검증 완료가 아니다. base-model 사전학습, LoRA, marker/mapper 선택과 최종 평가에 사용되지 않은 paired bulk+scRNA 코호트 및 원본 환자 manifest가 필요하다. 임의의 공개 cohort를 독립 데이터라고 간주하지 않는다.

## 실행을 이어가기 위해 필요한 입력

`spark-0110`은 이 PC에서 DNS로 해석되지 않았고 해당 서버 SSH 설정·키는 확인되지 않았다. 로컬 `.ssh`에는 GitHub known_hosts만 있다. 서버 문서에는 checkpoints와 paired counts가 server_1/spark-0110에 있다고 기록돼 있다. 접속 주소·계정·인증 경로 또는 counts/체크포인트 다운로드 위치가 필요하다. 암호나 비밀 키 내용을 보고서에 기록하지 않는다.

`results/server_acquisition_request.json`에는 20개 코호트×K×draw의 필요한 counts, 환자·샘플 행 매핑, 유전자·측정 mask, checkpoint SHA, condition, train/test manifest, receipt와 서버 디렉터리를 정리했다. K1 checkpoint 디렉터리의 정확한 경로는 로컬 provenance에 없어 확인이 필요하다. 공통으로 sampler/code dependencies, condition-builder, canonical gene table, base-model 학습 명단, annotation 모델 또는 검증된 mapper, 외부 코호트가 필요하다.

이번 작업으로 가설의 PBISC 지지 수준을 새로 높였다고 결론내릴 수는 없다. 확보한 것은 재현 가능한 평가 대상·counts와 대조군·실험 입력이다. 생성 데이터와 실행환경이 확보되면 ①의 매핑·독립성 점검을 마치고 ②→③→④ 순서로 실제 결과를 추가한다.

## 저장 파일과 검증

- `results/real_counts_manifest.json`: counts provenance·SHA256.
- `results/archive_inventory.csv`, `local_matrix_inventory.csv`, `acquisition_audit.json`: 입력 감사.
- `results/reference_metric_controls.csv`, `reference_metric_control_summary.csv`, `reference_metric_control_deltas.csv`: 실제 세포 지표 점검 결과.
- `results/bulk_swap_assignments.csv`: 환자가 겹치지 않는 471건 교환 배정.
- `inputs/rare_panel_pseudobulk_doses.npz`, `results/rare_panel_pseudobulk_dose_manifest.csv`: 105개 dose 입력과 출처.
- `results/rare_state_eligibility.csv`, `rare_sampling_power.csv`: 희귀 상태 후보와 표본 수 계산.
- `results/generation_provenance.json`, `fixed_checkpoint_seed_audit.csv`, `server_acquisition_request.json`: 체크포인트·seed·필요 서버 입력.
- `code/`: 모든 수행 코드. `VALIDATION.json`: 원본 count 일치, LoRA 환자 중복 0, 교환 환자·임상 strata 일치, dose 선형성·CPM 합·parent 비율, 공동 발현 파괴 대조군 점검 통과. **PBISC response validation은 미실행이다.**
'''
(OUT/'REPORT.md').write_text(report,encoding='utf-8')
hashes=[]
for folder in ['code','results','inputs']:
    for p in sorted((OUT/folder).rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        h=hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
        hashes.append(dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=h.hexdigest()))
json.dump(hashes,open(OUT/'FILE_MANIFEST.json','w'),indent=2)
print('Report saved; all_four_requested_stages_complete=False')
