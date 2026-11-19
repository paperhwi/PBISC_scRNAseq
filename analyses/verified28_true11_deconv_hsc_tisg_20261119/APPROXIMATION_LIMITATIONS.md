# 디컨볼루션 근사 구현의 한계

이 자료의 `BayesPrism-init/upd`, `BayesPrism Z`, `CIBERSORTx-Smode/Bmode`, `MuSiC`, `DWLS` 명칭은 원본 로컬 코드와 결과표의 식별자를 보존한 것이다. **공식 패키지 실행 결과가 아니며, 공식 구현과의 알고리즘·수치적 동등성도 검증하지 않았다.** 결과표의 수치와 기존 분석 코드는 변경하지 않았다.

| 보존된 명칭 | 실제 로컬 계산 | 해석 한계 |
| --- | --- | --- |
| `BayesPrism-init/upd` | 결정론적 multiplicative EM과 자체 참조 갱신 | 공식 BayesPrism의 Gibbs sampling 및 사후분포 주변화를 재현한 것으로 볼 수 없음 |
| `BayesPrism Z` | 위 근사 해에서 계산한 기대 발현 할당량 | 공식 BayesPrism의 사후 발현 복원 결과가 아님 |
| `CIBERSORTx-Smode/Bmode` | 자체 평균·분산/비율 보정과 선형 NuSVR | 공식 CIBERSORTx의 S/B-mode와 동일한 보정·최적화를 수행했다는 검증이 없음 |
| `MuSiC` | 자체 분산 가중 반복 NNLS | 공식 MuSiC 패키지와 동등성이 검증되지 않은 근사 기준선 |
| `DWLS` | 자체 가중 반복 NNLS | 공식 DWLS 패키지의 실행 또는 동일 구현으로 검증된 결과가 아님 |

따라서 모든 비교 결과는 **이 코드에 구현된 자체 근사 기준선**에 한정된다. 이 결과로 공식 BayesPrism·CIBERSORTx·MuSiC·DWLS의 실패, 공식 방법 대비 PBISC의 우월성, 또는 공식 베이지안 사후분포 복원을 입증할 수 없다. Python으로 구현했다는 사실만으로 원논문과 동등한 구현이 되는 것은 아니다.

원본 README·한국어 결과·코드의 설명에 남아 있는 “원논문 수식 그대로”, “충실 구현”, “결정적 우위” 등의 표현에는 **이 제한을 우선 적용**해야 한다. CSV의 방법명과 통계값은 원본과의 대응 및 해시 보존을 위해 유지했다. `PBISC 10/10` 역시 동일 발견 데이터에서 선택된 표적에 대한 탐색적 회복 결과이며, 독립 평가 자료를 사용한 방법 우월성 검정으로 읽어서는 안 된다. 공식 방법을 비교하려면 공식 구현 또는 동등성이 검증된 구현으로 별도 재평가해야 한다.

원 방법의 설명: [BayesPrism 원논문](https://www.nature.com/articles/s43018-022-00356-3), [BayesPrism 저자 저장소](https://github.com/Danko-Lab/BayesPrism), [CIBERSORTx 저자 프로토콜](https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/), [MuSiC 저자 문서](https://xuranw.github.io/MuSiC/articles/MuSiC.html).

이 GitHub 배포에는 생성된 PNG·PDF·SVG 및 합본 그림 PDF를 포함하지 않는다. 그림을 만드는 Python 코드와 기존 `RUN_RECEIPT.json`은 보존한다. 해당 영수증의 그림 경로·해시는 원래 로컬 실행의 기록이며, 그림 파일이 저장소에 포함되어 있다는 뜻은 아니다. Stage 1은 저장소 밖의 원자료와 원래 경로/의존성을 요구한다. 아래 원문 문서의 로컬 산출물 링크도 원본 기록으로 유지한다.
