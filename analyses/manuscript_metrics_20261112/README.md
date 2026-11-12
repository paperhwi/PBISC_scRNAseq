# Manuscript metrics (2026-11-12)

gene-mean 곡선, iLISI, composition / within-type 분해, 희귀세포 통계를 계산하는 재사용 가능한 Python 모듈입니다. 다른 분석 폴더와 달리 `F:/` 경로를 쓰지 않고, 모든 함수가 count 행렬과 `obs` 표를 인자로 받습니다.

> **출처 주의.** 11/10–11/12 최종 figure/분해 코드는 이 저장소에 없었고, PC 원본에 접근할 수 없어서 프로젝트 문서에 적힌 정의를 바탕으로 **새로 구현**했습니다. 방법은 기존 `three_model_20261022/code/ilisi_lib.py`, `f_rare_parent.py`, `analyze_ceiling_vs_real.py`와 같은 계열이지만 수치가 원본과 비트 단위로 같다고 보장하지 않습니다. 원본 노트북 수치와 대조한 뒤 논문에 사용하세요.

## 구성

| 파일 | 내용 |
| --- | --- |
| `pbisc_metrics/gene_mean.py` | 유전자별 평균 UMI·CP10k·검출률, real 순위/평균 구간별 generated/real 비, log1p 상관, depth-matched real 대조(binomial thinning), 산점도·순위 곡선 |
| `pbisc_metrics/ilisi.py` | reference-fitted 2,000 seurat_v3 HVG → PCA 30 → pairwise Harmony(batch = source) → iLISI(perplexity 30, k = 90, median). real B 상한 대비 `fraction_of_ceiling`, 환자 단위 bootstrap CI |
| `pbisc_metrics/decomposition.py` | 질병 효과를 `total = composition + within`으로 정확 분해(midpoint/Kitagawa, 잔차 없음). 환자 split-half(200회, 결정적)로 잡음 보정한 slope와 disattenuated r, 타입별 효과 |
| `pbisc_metrics/rare_stats.py` | parent lineage 분모의 fine label share, real A 대비 ratio와 판정(0.5× 미만 / 0.5–2× / 2× 초과 / 평가 불가), 환자 bootstrap CI, 소수 세포 플래그, 샘플 검출·이항 검정력 |
| `pbisc_metrics/synthetic.py` | 알려진 정답이 있는 합성 코호트(테스트·데모용) |
| `tests/` | pytest 16개 (항등 입력, 상한/하한, 분해 합 일치, 정답 약화 복원, 결정성 등) |
| `demo_manuscript_metrics.ipynb` | 합성 데이터 데모. **실행 결과(표·그림)를 포함**합니다 (이 저장소의 다른 노트북은 출력을 제거했지만 이것은 결과 확인용으로 남김) |

## 정의 요약

- **gene-mean**: real 평균 CP10k 순위로 구간을 나눠(1–500, 501–1,500, 1,501–2,000, 2,001–3,000, 3,001–4,000, 4,000+) `ratio = Σ생성 CP10k / Σreal CP10k`, `det_ratio`, 한 번도 생성되지 않은 유전자 수를 계산합니다. `real_dm_*` 열은 real을 생성본 깊이로 thinning한 대조입니다.
- **iLISI**: 두 source를 같은 세포 수로 맞춰 계산합니다(최댓값 2). `fraction_of_ceiling = (iLISI − 1)/(iLISI_real B − 1)`.
- **분해**: 샘플 pseudobulk `R_sg = Σ_t p_st μ_stg`. 그룹 수준 `p_t`, `μ_tg`로 `composition = Σ Δp_t · mean(μ_t)`, `within = Σ mean(p_t) · Δμ_t`. 환자를 그룹 안에서 반으로 나눠 real 효과를 두 번(A, B) 측정하고 `slope = cov(g, r_B)/cov(r_A, r_B)`, `r_dis = corr(g, r_B)/√corr(r_A, r_B)`. 이 조성 성분이 가장 불확실하므로 `split_cv`를 함께 봅니다.
- **희귀세포**: share는 전체가 아니라 parent lineage 안의 비율입니다. `ratio = share_gen / share_real A`. real B(`obs_ctrl`)에서 같은 상태가 매핑으로 사라지면 `not_assessable`.

## 사용 예

```python
import sys; sys.path.insert(0, "analyses/manuscript_metrics_20261112")
from pbisc_metrics import gene_mean as gm, ilisi as il, decomposition as dc, rare_stats as rs

t   = gm.compare_gene_means(real.X, real.var_names, gen.X, gen.var_names)       # counts 행렬
bands = gm.band_table(t, by="rank")
rep = il.ilisi_report(real_A.X, real_A.var_names, {"real_B": (...), "PBISC": (...)}, ceiling="real_B", n_per_source=20000)
rt  = rs.rare_table(obs_realA, obs_gen, fine="fine", parent="parent", unit="patient", sample="sample", obs_ctrl=obs_realB)
```
`dc.aggregate(...)` 입력에는 샘플 메타(`group`, `patient`)가 필요합니다. 자세한 흐름은 노트북을 보세요.

## 요구사항

공통 `requirements.txt`에 `scikit-misc`(seurat_v3 HVG), `harmonypy`, `pytest`가 필요합니다. 테스트: `cd analyses/manuscript_metrics_20261112 && python -m pytest -q tests`

## 한계

- 합성 데모의 수치는 논문 결과가 아닙니다.
- 분해와 split-half 보정은 생성본 자체의 잡음은 보정하지 않습니다(real 측 잡음만 보정).
- iLISI의 HVG는 reference real에서 고르고 Harmony는 한 seed입니다. 큰 코호트에서는 `n_per_source=20000`, `n_draws=3`을 권장합니다.
- 희귀 판정은 라벨 transfer 경로에 의존합니다(`obs_ctrl`로 real B 상한을 함께 보세요).
