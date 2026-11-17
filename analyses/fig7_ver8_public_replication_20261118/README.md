# Figure 7 ver8 — 공개 COVID-19 코호트 재현, 보충 그림 3 (2026-11-18)

Figure 7 ver8의 주장(앵커 4개, focus population 4개, 제외 대조 2개)을 ESKD와 겹치는 환자가 없는 공개 성인 COVID-19 PBMC 단일세포 코호트 4개에서 다시 확인한 코드입니다. 대상은 Stephenson 2021, COMBAT 2022, van der Wijst 2021, Liu 2021이고, WHO 등급이 있는 급성기 환자 219명(사망 23명)입니다. 결과 그림은 `PBISC_SuppFig_S3_Public_COVID_Replication_20261118`, 원고용 요약은 `docs/Public_Replication_ver8.md`입니다.

## 파일

| 파일 | 내용 |
| --- | --- |
| `code/run_and_plot_public_replication_v8.py` | 코호트별 표본 표 작성 → 환자별 급성기 평균 % → Spearman ρ(peak WHO), 사망 Mann–Whitney AUC(Hanley–McNeil SE, 하한 0.035), 나이·성별·기관 보정 순위 OLS, Stephenson count GLM → DerSimonian–Laird 메타분석 → HSPC 사망 시간창(첫 샘플 vs > 7일) → 보충 그림 3 a–d. `data/*.csv` 4개를 씁니다. 원래 판(`..._20261117`)에서 그림 제목·범례의 N만 222 → 219로 고쳤습니다 |
| `code/verify_public_independent.py` | 검산용 독립 재구현. obs parquet에서 환자 표, ρ, AUC, 환자·사망 수, 메타분석(`statsmodels.combine_effects`)을 다시 계산해 파이프라인 출력과 비교합니다 |
| `PBISC_public_replication_ver8_verification_20261118.ipynb` | 검산 노트북. 원본 README의 수치를 CSV와 대조하고(FAIL은 원본 README 오류), Liu 포함 민감도와 독립 재구현을 실행합니다. 저장소 규칙에 따라 출력은 지웠습니다 |
| `docs/Public_Replication_ver8.md` | 코호트별·주장별 재현 요약. 아래 오류를 고치고 3절에 검산 메모를 덧붙인 판입니다 |

## 실행

1. 입력(저장소에 없음): `F:\PBISC\public_covid_who_hspc_nk_20261115\obs_{stephenson2021,combat2022,vanderwijst2021}.parquet`, `F:\PBISC\HSPC_CD56brightNK_COVID_analysis_summary_20261116\10_CD56bright_public_replication\data\obs_liu2021.parquet`, 그리고 ESKD 비교값용 `15_Figure7_ver8_Main_and_Supp\data\` CSV와 `style.py`(`../fig7_ver8_main_supp_20261117/figure/style.py`와 같음).
2. 스크립트의 `ROOT`, `V8_STYLE_DIR`, `OUT_DIR`을 자기 환경에 맞게 바꿉니다. pandas, pyarrow, scipy, statsmodels, matplotlib이 필요합니다.
3. 검산 노트북은 `PUB_DATA`(16번 폴더 `data/`), `PBISC_ROOT`(위 입력 폴더들의 상위 폴더), `PUB_INDEP`(독립 재구현 스크립트 경로)를 환경 변수로 받습니다.

## 검산 (2026-11-18)

- **재실행:** 원본 obs로 다시 실행한 결과 CSV 4개가 저장본과 최대 3.4×10⁻¹⁴ 차이로 같았습니다.
- **독립 재구현:** 코호트별 ρ·AUC·P·환자 수·사망 수 21행과 메타분석 6개가 모두 일치했습니다.
- **원본 README 오류 8건** (`docs/`에 반영):
  - 전체 N 222 → **219**, Liu n 32 → **29**. 그림 제목과 범례도 고쳤습니다.
  - NK CD56bright 메타 P 1.0×10⁻⁴ → **8.0×10⁻⁵**.
  - Stephenson IgG의 "`% B` ρ +0.26, P 0.024"는 출력에 없는 값입니다. B ASC total 값으로 보입니다. 실제 IgG Spearman은 ρ +0.23, P 0.052, IgA는 +0.19, P 0.12입니다.
  - Stephenson T CD8 CM 사망 AUC "0.21–0.26" → **0.24, P 0.031**.
  - van der Wijst T CD8 CTL ρ −0.36, P 0.022 → **−0.34, P 0.030**. 사망 AUC 0.26, P 0.049 → **0.32, P 0.14 (유의하지 않음)**.
  - Stephenson T CD8 CTL ρ +0.05 → **+0.03**.
- **해석 주의:**
  - **세포형 대체 정의:** Stephenson에는 CD8 central-memory 라벨이 없어 `T CD8 CM` 자리에 `CD8.EM`을 썼습니다. T CD8 CTL은 코호트마다 `CD8.TE`, `TEMRA+TEFF+TEM`, `T8_Mem`, `CD8 memory`로 다릅니다.
  - **Liu 2021 제외:** Liu는 NK CD56bright·B ASC total 메타에서만 빠졌습니다. 넣어도 결론은 같습니다(ρ −0.25, P 0.0014; +0.28, P 0.0028).
  - **그림 a의 ESKD 비교값:** "B ASC (mature/total)" 줄의 ESKD 값은 B ASC IgG, "CD16 mono (% Mono)" 줄의 ESKD 값은 % PBMC입니다.
