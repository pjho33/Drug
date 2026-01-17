#!/usr/bin/env python3
"""
프로젝트 중심 폴더 구조 재구성 스크립트

기존 구조:
  Drug/
  ├── scripts/
  └── organized_project/

새 구조:
  Drug/
  ├── 2025-12-01_SMILES_Ligand_Design/
  │   ├── scripts/
  │   ├── data/
  │   └── README.md
  ├── 2025-12-05_SDF_Structures/
  └── ...
"""

import os
import shutil
from pathlib import Path
from datetime import datetime
import json

# 기존 폴더 경로
OLD_ROOT = Path("/home/pjho3/projects/Drug_local_20260117")
OLD_SCRIPTS = OLD_ROOT / "scripts"
OLD_ORGANIZED = OLD_ROOT / "organized_project"

# 새 폴더 경로 (Drug 폴더에 생성)
NEW_ROOT = Path("/home/pjho3/projects/Drug")

# 프로젝트 매핑 (기존 폴더 → 새 폴더)
PROJECT_MAPPING = {
    "01_SMILES_Ligand_Design_2025-12": "2025-12-01_SMILES_Ligand_Design",
    "02_SDF_Structures_2025-12": "2025-12-05_SDF_Structures",
    "07_DiffDock_2025-12": "2025-12-15_Docking_DiffDock",
    "09_Phase2_Results_2025-12_2026-01": "2025-12-20_Phase2_MD_Simulation",
    "10_Phase3_CHARMM_Preparation_2026-01": "2026-01-01_Phase3_CHARMM_Preparation",
    "11_Control_MD_2026-01-08": "2026-01-08_Control_MD",
    "12_Final_Complex_2026-01": "2026-01-10_Final_Complex",
    "13_Validation_2026-01": "2026-01-12_Validation",
    "00_Scripts_Organized": "Shared_Utils",
}

# 스크립트 파일명 패턴 → 프로젝트 매핑
SCRIPT_TO_PROJECT = {
    # MMPBSA 관련
    "mmpbsa": "2026-01-11_MMPBSA_Analysis",
    "ante_mmpbsa": "2026-01-11_MMPBSA_Analysis",
    "calculate_mmpbsa": "2026-01-11_MMPBSA_Analysis",
    
    # MD 시뮬레이션
    "run_glut1": "2025-12-20_Phase2_MD_Simulation",
    "run_control": "2026-01-08_Control_MD",
    "run_simple_md": "2025-12-20_Phase2_MD_Simulation",
    "analyze_control": "2026-01-08_Control_MD",
    
    # 복합체 생성
    "create_final_complex": "2026-01-10_Final_Complex",
    "create_proper_complex": "2026-01-10_Final_Complex",
    "merge_sdg": "2026-01-10_Final_Complex",
    "combine_pdb": "2026-01-10_Final_Complex",
    
    # Topology 관련
    "create_topology": "2026-01-10_Final_Complex",
    "create_psf": "2026-01-10_Final_Complex",
    "check_topology": "2026-01-10_Final_Complex",
    "fix_topology": "2026-01-10_Final_Complex",
    
    # 구조 검증
    "check_ligand": "2026-01-10_Final_Complex",
    "check_complex": "2026-01-10_Final_Complex",
    "check_solvent": "2026-01-10_Final_Complex",
    "diagnose_structure": "2026-01-11_MMPBSA_Analysis",
    "minimize_structure": "2026-01-11_MMPBSA_Analysis",
    
    # PDB 처리
    "fix_pdb": "2026-01-10_Final_Complex",
    
    # SDF/SMILES
    "create_sdg_sdf": "2025-12-05_SDF_Structures",
    "create_sdf": "2025-12-05_SDF_Structures",
    
    # 유틸리티
    "run_bg": "Shared_Utils",
    "check_logs": "Shared_Utils",
    "monitor": "Shared_Utils",
    "BACKGROUND": "Shared_Utils",
}


def create_project_structure(project_path: Path):
    """프로젝트 폴더 구조 생성"""
    subdirs = ["scripts", "data", "results", "docs"]
    for subdir in subdirs:
        (project_path / subdir).mkdir(parents=True, exist_ok=True)
    print(f"  ✅ Created: {project_path.name}")


def classify_script(script_name: str) -> str:
    """스크립트 파일명으로 프로젝트 분류"""
    script_lower = script_name.lower()
    
    for pattern, project in SCRIPT_TO_PROJECT.items():
        if pattern in script_lower:
            return project
    
    # 기본값: Shared_Utils
    return "Shared_Utils"


def generate_readme(project_path: Path, project_name: str):
    """프로젝트 README.md 생성"""
    date_str = project_path.name.split('_')[0] if '_' in project_path.name else "N/A"
    
    # 프로젝트 설명
    descriptions = {
        "SMILES_Ligand_Design": "리간드 SMILES 설계 및 검증",
        "SDF_Structures": "3D 구조 생성 및 최적화",
        "Receptor_Preparation": "GLUT1 수용체 준비 및 글리코실화",
        "Docking_DiffDock": "DiffDock을 이용한 도킹 시뮬레이션",
        "Phase2_MD_Simulation": "Phase2 MD 시뮬레이션 (OpenMM)",
        "Phase3_CHARMM_Preparation": "CHARMM36 시스템 준비",
        "Control_MD": "Control 시뮬레이션",
        "Final_Complex": "최종 복합체 생성 및 검증",
        "MMPBSA_Analysis": "MM/PBSA 결합 자유 에너지 계산",
        "Validation": "결과 검증 및 분석",
        "Shared_Utils": "공통 유틸리티 및 헬퍼 함수",
    }
    
    desc = "프로젝트 설명"
    for key, value in descriptions.items():
        if key in project_name:
            desc = value
            break
    
    readme_content = f"""# {project_name}

**시작일:** {date_str}

## 📋 프로젝트 개요

{desc}

## 📁 폴더 구조

```
{project_path.name}/
├── scripts/      # 이 프로젝트 전용 스크립트
├── data/         # 입력 데이터 (작은 파일만, 대용량은 제외)
├── results/      # 출력 결과 (작은 파일만)
├── docs/         # 프로젝트 문서
└── README.md     # 이 파일
```

## 🔧 스크립트 목록

"""
    
    # 스크립트 목록 추가
    scripts_dir = project_path / "scripts"
    if scripts_dir.exists():
        script_files = sorted(scripts_dir.glob("*.py")) + sorted(scripts_dir.glob("*.sh"))
        if script_files:
            readme_content += "| 파일명 | 설명 |\n|--------|------|\n"
            for script in script_files:
                readme_content += f"| `{script.name}` | TODO: 설명 추가 |\n"
        else:
            readme_content += "(스크립트 없음)\n"
    
    readme_content += """
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
"""
    
    readme_path = project_path / "README.md"
    readme_path.write_text(readme_content)
    print(f"  ✅ Generated: {readme_path.relative_to(NEW_ROOT)}")


def main():
    print("=" * 80)
    print("프로젝트 중심 폴더 구조 재구성")
    print("=" * 80)
    print()
    print(f"기존 경로: {OLD_ROOT}")
    print(f"새 경로: {NEW_ROOT}")
    print()
    
    # 1. 새 프로젝트 폴더 생성
    print("Step 1: 새 프로젝트 폴더 생성")
    print("-" * 80)
    
    for old_name, new_name in PROJECT_MAPPING.items():
        new_project_path = NEW_ROOT / new_name
        create_project_structure(new_project_path)
    
    # MMPBSA_Analysis 폴더 추가 생성
    mmpbsa_path = NEW_ROOT / "2026-01-11_MMPBSA_Analysis"
    create_project_structure(mmpbsa_path)
    
    print()
    
    # 2. 스크립트 분류 및 복사
    print("Step 2: 스크립트 분류 및 복사")
    print("-" * 80)
    
    script_classification = {}
    
    if OLD_SCRIPTS.exists():
        for script_file in OLD_SCRIPTS.glob("*"):
            if script_file.is_file() and (script_file.suffix in ['.py', '.sh', '.md'] or script_file.name.startswith('BACKGROUND')):
                # Zone.Identifier 파일 제외
                if 'Zone.Identifier' in script_file.name:
                    continue
                
                project_name = classify_script(script_file.name)
                
                if project_name not in script_classification:
                    script_classification[project_name] = []
                script_classification[project_name].append(script_file)
    
    # 스크립트 복사
    for project_name, scripts in script_classification.items():
        dest_dir = NEW_ROOT / project_name / "scripts"
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n{project_name}:")
        for script in scripts:
            dest_file = dest_dir / script.name
            shutil.copy2(script, dest_file)
            print(f"  ✅ {script.name}")
    
    print()
    
    # 3. organized_project 내 스크립트 복사
    print("Step 3: organized_project 내 스크립트 복사")
    print("-" * 80)
    
    for old_name, new_name in PROJECT_MAPPING.items():
        old_path = OLD_ORGANIZED / old_name
        if not old_path.exists():
            continue
        
        new_path = NEW_ROOT / new_name
        scripts_dest = new_path / "scripts"
        
        # 스크립트 파일 찾기
        script_files = list(old_path.rglob("*.py")) + list(old_path.rglob("*.sh")) + list(old_path.rglob("*.md"))
        
        if script_files:
            print(f"\n{new_name}:")
            for script in script_files:
                # .git 폴더 제외
                if '.git' in str(script):
                    continue
                
                dest_file = scripts_dest / script.name
                # 중복 방지
                if dest_file.exists():
                    dest_file = scripts_dest / f"{script.stem}_from_organized{script.suffix}"
                
                shutil.copy2(script, dest_file)
                print(f"  ✅ {script.name}")
    
    print()
    
    # 4. README 생성
    print("Step 4: README.md 생성")
    print("-" * 80)
    
    for project_dir in sorted(NEW_ROOT.glob("*")):
        if project_dir.is_dir() and not project_dir.name.startswith('.'):
            generate_readme(project_dir, project_dir.name)
    
    print()
    
    # 5. 분류 결과 저장
    print("Step 5: 분류 결과 저장")
    print("-" * 80)
    
    classification_report = {
        "timestamp": datetime.now().isoformat(),
        "projects": {},
    }
    
    for project_dir in sorted(NEW_ROOT.glob("*")):
        if project_dir.is_dir() and not project_dir.name.startswith('.'):
            scripts_dir = project_dir / "scripts"
            if scripts_dir.exists():
                scripts = [f.name for f in scripts_dir.glob("*") if f.is_file()]
                classification_report["projects"][project_dir.name] = {
                    "script_count": len(scripts),
                    "scripts": scripts,
                }
    
    report_path = NEW_ROOT / "RESTRUCTURE_REPORT.json"
    with open(report_path, 'w') as f:
        json.dump(classification_report, f, indent=2)
    
    print(f"✅ 분류 보고서: {report_path}")
    print()
    
    # 6. 요약
    print("=" * 80)
    print("재구성 완료!")
    print("=" * 80)
    print()
    print("생성된 프로젝트:")
    for project_name in sorted(classification_report["projects"].keys()):
        info = classification_report["projects"][project_name]
        print(f"  - {project_name}: {info['script_count']}개 스크립트")
    
    print()
    print("다음 단계:")
    print("  1. 새 폴더 구조 확인: ls -la /home/pjho3/projects/Drug/")
    print("  2. Git 상태 확인: cd /home/pjho3/projects/Drug && git status")
    print("  3. Git 추가: git add .")
    print("  4. 커밋: git commit -m 'Restructure: Project-centric folder organization'")
    print()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ 오류 발생: {e}")
        import traceback
        traceback.print_exc()
