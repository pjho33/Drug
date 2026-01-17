# 2026-01-01_Phase3_CHARMM_Preparation

**시작일:** 2026-01-01

## 📋 프로젝트 개요

CHARMM36 시스템 준비

## 📁 폴더 구조

```
2026-01-01_Phase3_CHARMM_Preparation/
├── scripts/      # 이 프로젝트 전용 스크립트
├── data/         # 입력 데이터 (작은 파일만, 대용량은 제외)
├── results/      # 출력 결과 (작은 파일만)
├── docs/         # 프로젝트 문서
└── README.md     # 이 파일
```

## 🔧 스크립트 목록

| 파일명 | 설명 |
|--------|------|
| `extract_tripod_coords.py` | TODO: 설명 추가 |
| `fix_pdb_format.py` | TODO: 설명 추가 |
| `merge_tripod_to_membrane.py` | TODO: 설명 추가 |
| `omm_barostat.py` | TODO: 설명 추가 |
| `omm_barostat_from_organized.py` | TODO: 설명 추가 |
| `omm_readinputs.py` | TODO: 설명 추가 |
| `omm_readinputs_from_organized.py` | TODO: 설명 추가 |
| `omm_readparams.py` | TODO: 설명 추가 |
| `omm_readparams_from_organized.py` | TODO: 설명 추가 |
| `omm_restraints.py` | TODO: 설명 추가 |
| `omm_restraints_from_organized.py` | TODO: 설명 추가 |
| `omm_rewrap.py` | TODO: 설명 추가 |
| `omm_rewrap_from_organized.py` | TODO: 설명 추가 |
| `omm_vfswitch.py` | TODO: 설명 추가 |
| `omm_vfswitch_from_organized.py` | TODO: 설명 추가 |
| `openmm_run.py` | TODO: 설명 추가 |
| `openmm_run_from_organized.py` | TODO: 설명 추가 |
| `prepare_simulation.py` | TODO: 설명 추가 |
| `run_control_gpu1.py` | TODO: 설명 추가 |
| `run_control_gpu1_from_organized.py` | TODO: 설명 추가 |
| `run_experimental_gpu0.py` | TODO: 설명 추가 |
| `run_glyco_gpu0.py` | TODO: 설명 추가 |
| `run_md.py` | TODO: 설명 추가 |
| `run_md_from_organized.py` | TODO: 설명 추가 |
| `run_phase3_md.py` | TODO: 설명 추가 |
| `update_psf_for_tripod.py` | TODO: 설명 추가 |
| `check_progress.sh` | TODO: 설명 추가 |
| `check_progress_from_organized.sh` | TODO: 설명 추가 |
| `run_both_gpus.sh` | TODO: 설명 추가 |
| `run_both_simulations.sh` | TODO: 설명 추가 |
| `run_phase3_amber.sh` | TODO: 설명 추가 |
| `run_phase3_hybrid.sh` | TODO: 설명 추가 |

## 📊 데이터

- **입력:** `data/` 폴더
- **출력:** `results/` 폴더
- **대용량 파일:** Git에서 제외 (`.gitignore` 참조)

## 🚀 실행 방법

```bash
cd scripts/
python run_*.py
```

## 📝 노트

- 추가 정보 및 메모

---

**최종 수정:** {datetime.now().strftime("%Y-%m-%d")}
