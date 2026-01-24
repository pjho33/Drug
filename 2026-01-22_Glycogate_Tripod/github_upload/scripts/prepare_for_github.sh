#!/usr/bin/env bash
# GitHub 업로드를 위한 프로젝트 정리 스크립트
# 대용량 데이터 파일 제외, 스크립트와 결과 보고서 위주

set -euo pipefail

PROJECT_DIR="/home/pjho3tr/projects/Drug/2026-01-22_Glycogate_Tripod"
ARCHIVE_DIR="$PROJECT_DIR/github_upload"

echo "=========================================="
echo "GitHub 업로드용 프로젝트 정리"
echo "=========================================="

# 기존 아카이브 디렉토리 삭제
if [ -d "$ARCHIVE_DIR" ]; then
    rm -rf "$ARCHIVE_DIR"
fi

# 새 디렉토리 구조 생성
mkdir -p "$ARCHIVE_DIR"/{scripts,results,docs,data}

echo "1. 스크립트 복사..."
# 모든 Python 및 Shell 스크립트 복사
cp -r "$PROJECT_DIR/scripts/"*.py "$ARCHIVE_DIR/scripts/" 2>/dev/null || true
cp -r "$PROJECT_DIR/scripts/"*.sh "$ARCHIVE_DIR/scripts/" 2>/dev/null || true

echo "2. 결과 파일 복사 (CSV, PNG, DOCX만)..."
# 1ns 결과
mkdir -p "$ARCHIVE_DIR/results/md_tripod_1ns"
cp "$PROJECT_DIR/results/md_tripod_1ns/"*.csv "$ARCHIVE_DIR/results/md_tripod_1ns/" 2>/dev/null || true
cp "$PROJECT_DIR/results/md_tripod_1ns/"*.png "$ARCHIVE_DIR/results/md_tripod_1ns/" 2>/dev/null || true

# 10ns 결과
mkdir -p "$ARCHIVE_DIR/results/md_tripod_10ns"
cp "$PROJECT_DIR/results/md_tripod_10ns/"*.csv "$ARCHIVE_DIR/results/md_tripod_10ns/" 2>/dev/null || true
cp "$PROJECT_DIR/results/md_tripod_10ns/"*.png "$ARCHIVE_DIR/results/md_tripod_10ns/" 2>/dev/null || true

# 100ns 결과
mkdir -p "$ARCHIVE_DIR/results/md_tripod_100ns"
cp "$PROJECT_DIR/results/md_tripod_100ns/"*.csv "$ARCHIVE_DIR/results/md_tripod_100ns/" 2>/dev/null || true
cp "$PROJECT_DIR/results/md_tripod_100ns/"*.png "$ARCHIVE_DIR/results/md_tripod_100ns/" 2>/dev/null || true

# 비교 플롯
cp "$PROJECT_DIR/results/"*.png "$ARCHIVE_DIR/results/" 2>/dev/null || true

# 워드 보고서
cp "$PROJECT_DIR/results/"*.docx "$ARCHIVE_DIR/results/" 2>/dev/null || true

echo "3. 입력 파일 복사 (작은 파일만)..."
# OpenMM 입력 파일들
mkdir -p "$ARCHIVE_DIR/data/openmm"
cp "$PROJECT_DIR/data/solution_builder/openmm/"*.inp "$ARCHIVE_DIR/data/openmm/" 2>/dev/null || true
cp "$PROJECT_DIR/data/solution_builder/openmm/openmm_run.py" "$ARCHIVE_DIR/data/openmm/" 2>/dev/null || true

# PSF, CRD는 크기 확인 후 복사 (50MB 이하만)
for file in "$PROJECT_DIR/data/solution_builder/openmm/"*.{psf,crd}; do
    if [ -f "$file" ]; then
        size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0)
        if [ "$size" -lt 52428800 ]; then  # 50MB
            cp "$file" "$ARCHIVE_DIR/data/openmm/" 2>/dev/null || true
        fi
    fi
done

echo "4. README 생성..."
cat > "$ARCHIVE_DIR/README.md" << 'EOF'
# Tripod MD Simulation Project

## Overview
Molecular Dynamics simulation study of Tripod structure using OpenMM.

## Simulations Performed
- **1ns**: Initial short simulation (250,000 steps)
- **10ns**: Medium-length simulation (2,500,000 steps)
- **100ns**: Long-term equilibrium simulation (25,000,000 steps)

## Key Findings
- Tripod structure shows strong preference for compact conformations
- Maximum arm length decreases from 74.3 Å (1ns) to 30.4 Å (100ns)
- Only 5.4% of conformations exceed 50 Å in 100ns simulation
- Radius of gyration: 20.5 ± 7.3 Å (100ns)

## Directory Structure
```
.
├── scripts/           # Analysis and execution scripts
├── results/           # Analysis results (CSV, PNG, DOCX)
│   ├── md_tripod_1ns/
│   ├── md_tripod_10ns/
│   └── md_tripod_100ns/
├── data/              # Input files
│   └── openmm/
└── docs/              # Documentation
```

## Files Included
- **Scripts**: Python analysis scripts, shell execution scripts
- **Results**: CSV data, PNG plots, Word report
- **Data**: OpenMM input files (.inp), topology files (if < 50MB)

## Files Excluded (Large Data)
- Trajectory files (.dcd) - ~2GB for 100ns
- Restart files (.rst) - ~26MB each
- Checkpoint files (.chk)
- Large topology files (> 50MB)

## Requirements
- Python 3.10+
- OpenMM 8.0+
- MDAnalysis
- pandas, matplotlib
- python-docx (for report generation)

## Usage

### Running Simulations
```bash
# 1ns simulation
./scripts/run_1ns_only.sh

# 10ns simulation
./scripts/run_10ns_only.sh

# 100ns simulation
./scripts/run_100ns_only.sh
```

### Analysis
```bash
# Analyze results
conda activate drug-md
python scripts/analyze_tripod_1ns.py
python scripts/analyze_tripod_10ns.py
python scripts/analyze_tripod_100ns.py

# Generate comparison plots
python scripts/compare_all_tripod.py

# Generate Word report
python scripts/generate_report.py
```

## Results Summary

| Metric | 1ns | 10ns | 100ns |
|--------|-----|------|-------|
| Max Arm (Å) | 74.3 ± 8.2 | 34.2 ± 12.6 | 30.4 ± 11.1 |
| Rg (Å) | 44.1 ± 0.9 | 37.7 ± 6.4 | 20.5 ± 7.3 |
| P(≥40Å) | 100% | 35% | 18% |
| P(≥50Å) | 100% | 13% | 5.4% |

## Citation
If you use this data, please cite:
- OpenMM: http://openmm.org
- MDAnalysis: https://www.mdanalysis.org

## Contact
Project Date: January 2026

## License
This project is for research purposes.
EOF

echo "5. .gitignore 생성..."
cat > "$ARCHIVE_DIR/.gitignore" << 'EOF'
# Large data files
*.dcd
*.rst
*.chk
*.crd
*.psf
*.pdb

# Except small files
!data/openmm/*.inp
!data/openmm/openmm_run.py

# Log files
*.log
nohup.out

# Python
__pycache__/
*.pyc
*.pyo
*.egg-info/

# OS
.DS_Store
Thumbs.db

# IDE
.vscode/
.idea/
*.swp
EOF

echo "6. 파일 크기 확인..."
du -sh "$ARCHIVE_DIR"
echo ""
echo "디렉토리별 크기:"
du -sh "$ARCHIVE_DIR"/*

echo ""
echo "=========================================="
echo "✅ GitHub 업로드 준비 완료!"
echo "=========================================="
echo "위치: $ARCHIVE_DIR"
echo ""
echo "다음 단계:"
echo "1. cd $ARCHIVE_DIR"
echo "2. git init"
echo "3. git add ."
echo "4. git commit -m 'Initial commit: Tripod MD simulation results'"
echo "5. git remote add origin <your-repo-url>"
echo "6. git push -u origin main"
echo "=========================================="
