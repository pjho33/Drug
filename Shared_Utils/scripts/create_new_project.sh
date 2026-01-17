#!/bin/bash
# 새 프로젝트 생성 헬퍼 스크립트
# 사용법: ./create_new_project.sh "Project_Name" "프로젝트 설명"

set -e

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 인자 확인
if [ $# -lt 1 ]; then
    echo -e "${RED}사용법: $0 \"Project_Name\" [\"프로젝트 설명\"]${NC}"
    echo ""
    echo "예시:"
    echo "  $0 \"Metadynamics_Analysis\" \"Metadynamics 시뮬레이션 분석\""
    echo "  $0 \"FEP_Calculation\""
    exit 1
fi

PROJECT_NAME="$1"
PROJECT_DESC="${2:-프로젝트 설명}"
PROJECT_DATE=$(date +%Y-%m-%d)
FULL_NAME="${PROJECT_DATE}_${PROJECT_NAME}"

DRUG_ROOT="$HOME/projects/Drug"

echo "=========================================="
echo "새 프로젝트 생성"
echo "=========================================="
echo "프로젝트명: $PROJECT_NAME"
echo "전체 폴더명: $FULL_NAME"
echo "설명: $PROJECT_DESC"
echo "날짜: $PROJECT_DATE"
echo ""

# Drug 폴더로 이동
cd "$DRUG_ROOT"

# 최신 상태로 업데이트
echo -e "${YELLOW}Step 1: Git Pull (최신 상태로 업데이트)${NC}"
git pull
echo ""

# 폴더 생성
echo -e "${YELLOW}Step 2: 폴더 구조 생성${NC}"
mkdir -p "$FULL_NAME"/{scripts,data,results,docs}
echo "  ✅ $FULL_NAME/scripts/"
echo "  ✅ $FULL_NAME/data/"
echo "  ✅ $FULL_NAME/results/"
echo "  ✅ $FULL_NAME/docs/"
echo ""

# README 생성
echo -e "${YELLOW}Step 3: README.md 생성${NC}"
cat > "$FULL_NAME/README.md" << EOF
# ${PROJECT_NAME}

**시작일:** ${PROJECT_DATE}

## 📋 프로젝트 개요

${PROJECT_DESC}

## 📁 폴더 구조

\`\`\`
${FULL_NAME}/
├── scripts/      # 이 프로젝트 전용 스크립트
├── data/         # 입력 데이터 (작은 파일만, 대용량은 제외)
├── results/      # 출력 결과 (작은 파일만)
├── docs/         # 프로젝트 문서
└── README.md     # 이 파일
\`\`\`

## 🔧 스크립트 목록

(작업 중)

## 📊 데이터

- **입력:** \`data/\` 폴더
- **출력:** \`results/\` 폴더
- **대용량 파일:** Git에서 제외 (\`.gitignore\` 참조)

## 🚀 실행 방법

\`\`\`bash
cd scripts/
python run_*.py
\`\`\`

## 📝 노트

- 추가 정보 및 메모

---

**최종 수정:** ${PROJECT_DATE}
EOF

echo "  ✅ README.md 생성 완료"
echo ""

# 예시 스크립트 생성
echo -e "${YELLOW}Step 4: 예시 스크립트 생성${NC}"
cat > "$FULL_NAME/scripts/run_example.py" << 'EOF'
#!/usr/bin/env python3
"""
예시 스크립트 - 실제 작업에 맞게 수정하세요
"""

import os
from pathlib import Path

def main():
    print("=" * 80)
    print("프로젝트 실행")
    print("=" * 80)
    print()
    
    # 작업 코드 작성
    print("작업 시작...")
    
    print()
    print("✅ 완료!")

if __name__ == "__main__":
    main()
EOF

chmod +x "$FULL_NAME/scripts/run_example.py"
echo "  ✅ scripts/run_example.py 생성 완료"
echo ""

# Git 추가
echo -e "${YELLOW}Step 5: Git 커밋 및 Push${NC}"
git add "$FULL_NAME/"
git commit -m "New Project: Initialize ${FULL_NAME}

- Created project structure
- Added README.md
- Added example script"

git push
echo ""

# 완료
echo "=========================================="
echo -e "${GREEN}✅ 프로젝트 생성 완료!${NC}"
echo "=========================================="
echo ""
echo "프로젝트 경로: $DRUG_ROOT/$FULL_NAME"
echo ""
echo "다음 단계:"
echo "  1. cd $FULL_NAME/scripts"
echo "  2. 스크립트 작성"
echo "  3. git add . && git commit -m 'Update' && git push"
echo ""
echo "다른 컴퓨터에서 동기화:"
echo "  cd ~/projects/Drug && git pull"
echo ""
