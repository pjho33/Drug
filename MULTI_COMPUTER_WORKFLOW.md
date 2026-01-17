# 다중 컴퓨터 Git 워크플로우 가이드

## 🖥️ 시나리오

**3대의 컴퓨터에서 작업:**
- 컴퓨터 A: 현재 컴퓨터 (메인)
- 컴퓨터 B: 두 번째 컴퓨터
- 컴퓨터 C: 세 번째 컴퓨터

**작업 방식:**
- 같은 프로젝트를 여러 컴퓨터에서 작업
- 서로 다른 프로젝트를 동시에 진행

---

## 🚀 초기 설정 (각 컴퓨터마다 한 번만)

### 컴퓨터 B, C에서 처음 시작할 때

```bash
# 1. GitHub에서 클론
cd ~/projects
git clone https://github.com/pjho33/Drug.git
cd Drug

# 2. Git 사용자 정보 설정 (처음 한 번만)
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"

# 3. 현재 구조 확인
ls -la
```

---

## 📋 일상적인 워크플로우

### 시나리오 1: 기존 프로젝트 작업 (여러 컴퓨터)

**예시: 컴퓨터 A에서 MMPBSA 작업 → 컴퓨터 B에서 계속**

#### 컴퓨터 A에서 (작업 완료 후)

```bash
cd ~/projects/Drug

# 1. 변경사항 확인
git status

# 2. 변경된 파일 추가
git add 2026-01-11_MMPBSA_Analysis/scripts/new_script.py

# 또는 해당 프로젝트 전체 추가
git add 2026-01-11_MMPBSA_Analysis/

# 3. 커밋 (의미있는 메시지)
git commit -m "MMPBSA: Add new analysis script for binding energy"

# 4. GitHub에 Push
git push
```

#### 컴퓨터 B에서 (작업 시작 전)

```bash
cd ~/projects/Drug

# 1. 최신 변경사항 가져오기
git pull

# 2. 작업 시작
cd 2026-01-11_MMPBSA_Analysis/scripts
# 작업...

# 3. 완료 후 다시 커밋 & Push
git add .
git commit -m "MMPBSA: Update analysis parameters"
git push
```

---

### 시나리오 2: 새 프로젝트 시작

**예시: 컴퓨터 B에서 새 프로젝트 시작**

```bash
cd ~/projects/Drug

# 1. 최신 상태로 업데이트
git pull

# 2. 새 프로젝트 폴더 생성 (날짜 포함)
mkdir 2026-01-17_New_Project_Name
cd 2026-01-17_New_Project_Name

# 3. 표준 폴더 구조 생성
mkdir scripts data results docs

# 4. README.md 생성
cat > README.md << 'EOF'
# New_Project_Name

**시작일:** 2026-01-17

## 📋 프로젝트 개요

프로젝트 설명...

## 📁 폴더 구조

```
2026-01-17_New_Project_Name/
├── scripts/      # 스크립트
├── data/         # 입력 데이터
├── results/      # 출력 결과
├── docs/         # 문서
└── README.md
```

## 🔧 스크립트 목록

(작업 중)

---

**최종 수정:** 2026-01-17
EOF

# 5. 첫 스크립트 작성
cd scripts
# 스크립트 작성...

# 6. Git에 추가
cd ~/projects/Drug
git add 2026-01-17_New_Project_Name/
git commit -m "New Project: Initialize 2026-01-17_New_Project_Name"
git push
```

---

### 시나리오 3: 서로 다른 프로젝트 동시 작업

**컴퓨터 A: MMPBSA 작업**
**컴퓨터 B: Validation 작업**

#### 컴퓨터 A

```bash
cd ~/projects/Drug

# 작업 전 항상 pull
git pull

# MMPBSA 작업
cd 2026-01-11_MMPBSA_Analysis/scripts
# 작업...

# 커밋 & Push
cd ~/projects/Drug
git add 2026-01-11_MMPBSA_Analysis/
git commit -m "MMPBSA: Add new feature"
git push
```

#### 컴퓨터 B (동시에)

```bash
cd ~/projects/Drug

# 작업 전 항상 pull
git pull

# Validation 작업
cd 2026-01-12_Validation/scripts
# 작업...

# 커밋 & Push
cd ~/projects/Drug
git add 2026-01-12_Validation/
git commit -m "Validation: Add test cases"
git push
```

**충돌 없음!** 서로 다른 프로젝트 폴더이므로 자동으로 병합됩니다.

---

## ⚠️ 충돌 해결 (같은 파일 수정 시)

**상황:** 컴퓨터 A와 B에서 같은 파일을 동시에 수정

```bash
# 컴퓨터 B에서 push 시도
git push

# 오류 발생: Updates were rejected
# 해결:

# 1. 최신 변경사항 가져오기
git pull

# 2. 충돌 발생 시 메시지 확인
# CONFLICT (content): Merge conflict in 2026-01-11_MMPBSA_Analysis/scripts/script.py

# 3. 충돌 파일 수동 수정
nano 2026-01-11_MMPBSA_Analysis/scripts/script.py

# 충돌 마커 찾기:
# <<<<<<< HEAD
# 내 변경사항
# =======
# 다른 컴퓨터의 변경사항
# >>>>>>> origin/main

# 4. 원하는 버전으로 수정 후 저장

# 5. 충돌 해결 완료 표시
git add 2026-01-11_MMPBSA_Analysis/scripts/script.py
git commit -m "Merge: Resolve conflict in script.py"
git push
```

---

## 🎯 베스트 프랙티스

### 1. 작업 시작 전 항상 Pull

```bash
cd ~/projects/Drug
git pull
```

### 2. 자주 커밋 & Push

```bash
# 작은 단위로 자주 커밋
git add specific_file.py
git commit -m "Add feature X"
git push

# 하루 작업 끝날 때 반드시 Push
```

### 3. 의미있는 커밋 메시지

```bash
# ✅ 좋은 예
git commit -m "MMPBSA: Fix GB calculation error in run_mmpbsa_gb.py"
git commit -m "Validation: Add RMSD analysis script"
git commit -m "Final_Complex: Update topology generation for SDG ligand"

# ❌ 나쁜 예
git commit -m "update"
git commit -m "fix"
git commit -m "test"
```

**형식:** `[프로젝트명]: [간단한 설명]`

### 4. 프로젝트별로 작업 분리

```bash
# 한 번에 하나의 프로젝트만 수정
git add 2026-01-11_MMPBSA_Analysis/
git commit -m "MMPBSA: ..."

# 여러 프로젝트 동시 수정 시 분리 커밋
git add 2026-01-11_MMPBSA_Analysis/
git commit -m "MMPBSA: ..."

git add 2026-01-12_Validation/
git commit -m "Validation: ..."

git push
```

### 5. 대용량 파일 주의

```bash
# .gitignore가 자동으로 제외하지만, 확인
git status

# 대용량 파일이 추가되려고 하면
git rm --cached large_file.dcd
```

---

## 📊 상태 확인 명령어

```bash
# 현재 상태
git status

# 최근 커밋 히스토리
git log --oneline -10

# 원격 저장소와 차이
git fetch
git status

# 특정 파일 변경 이력
git log --follow -- 2026-01-11_MMPBSA_Analysis/scripts/script.py

# 누가 언제 수정했는지
git blame 2026-01-11_MMPBSA_Analysis/scripts/script.py
```

---

## 🔧 유용한 Git 명령어

### 변경사항 취소

```bash
# 작업 디렉토리 변경 취소 (커밋 전)
git restore script.py

# 스테이징 취소 (add 취소)
git restore --staged script.py

# 마지막 커밋 수정 (push 전)
git commit --amend -m "New message"
```

### 브랜치 사용 (고급)

```bash
# 새 기능 개발 시 브랜치 생성
git checkout -b feature/new-analysis

# 작업 후 커밋
git add .
git commit -m "Add new analysis"

# 메인으로 병합
git checkout main
git merge feature/new-analysis
git push

# 브랜치 삭제
git branch -d feature/new-analysis
```

---

## 📝 새 프로젝트 추가 템플릿

```bash
#!/bin/bash
# 새 프로젝트 생성 스크립트

PROJECT_DATE="2026-01-17"
PROJECT_NAME="New_Project_Name"
FULL_NAME="${PROJECT_DATE}_${PROJECT_NAME}"

cd ~/projects/Drug
git pull

# 폴더 생성
mkdir -p "$FULL_NAME"/{scripts,data,results,docs}

# README 생성
cat > "$FULL_NAME/README.md" << EOF
# ${PROJECT_NAME}

**시작일:** ${PROJECT_DATE}

## 📋 프로젝트 개요

프로젝트 설명...

## 📁 폴더 구조

\`\`\`
${FULL_NAME}/
├── scripts/      # 스크립트
├── data/         # 입력 데이터
├── results/      # 출력 결과
├── docs/         # 문서
└── README.md
\`\`\`

## 🔧 스크립트 목록

(작업 중)

---

**최종 수정:** ${PROJECT_DATE}
EOF

# Git 추가
git add "$FULL_NAME/"
git commit -m "New Project: Initialize ${FULL_NAME}"
git push

echo "✅ 프로젝트 생성 완료: $FULL_NAME"
```

---

## 🎉 요약

### 일상 워크플로우 (3단계)

1. **작업 시작 전**
   ```bash
   git pull
   ```

2. **작업 중**
   ```bash
   # 스크립트 작성, 수정...
   ```

3. **작업 완료 후**
   ```bash
   git add .
   git commit -m "Project: Description"
   git push
   ```

### 핵심 원칙

- ✅ 작업 시작 전 항상 `git pull`
- ✅ 자주 커밋, 자주 Push
- ✅ 의미있는 커밋 메시지
- ✅ 프로젝트별 폴더 분리
- ✅ 대용량 파일 제외 (.gitignore)

---

**작성일:** 2026-01-17
**버전:** 1.0
