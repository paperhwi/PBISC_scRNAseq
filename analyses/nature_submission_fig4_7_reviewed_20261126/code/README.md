# PBISC Figures 4–7 임상 재판정·검산 재현 코드 (`2026-11-26`)

모든 스크립트는 저장소 상대경로(`Path(__file__).resolve().parent.parent`)만 사용하며 Windows, macOS, Linux 어디서나 로컬 절대경로 의존성 없이 실행됩니다.

## 실행 순서

```bash
pip install -r ../../../requirements.txt

# 1. van der Wijst 공식 임상 사망 판정(42명/7사망), 외부 4개 코호트 메타분석(219명/24사망),
#    Supp. Fig. 23c 16개 고유 HSPC 프로그램(62명), Fig. 5f 상관계수 독립 재현 검증
python code/01_reproduce_adjudicated_tables.py

# (선택) 23개 시트 Source Data Excel 통합문서(.xlsx) 로컬 생성이 필요한 경우
python code/01_reproduce_adjudicated_tables.py --export-xlsx

# (선택) 저자 공개 임상 원본(PMC8601717.2 / Nature Medicine MOESM3) 직접 재다운로드 및 SHA-256 대조
python code/02_download_clinical_sources.py

# 2. 22개 CSV 테이블 매니페스트·SHA-256, Figures 4–7 핵심 통계량, 원고/참고문헌(1->20) 전수 검산
python code/03_verify_all_fig4_7.py
```
