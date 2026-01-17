#!/usr/bin/env python3
"""
Git 동기화 대상 파일 확인 스크립트
"""

import os
from pathlib import Path
import subprocess

ROOT = Path("/home/pjho3/projects/Drug_local_20260117")

print("=" * 80)
print("Git 동기화 대상 파일 분석")
print("=" * 80)
print()

# 1. 스크립트 파일 통계
print("1. 스크립트 파일 통계")
print("-" * 80)

script_extensions = ['.py', '.sh', '.md']
config_extensions = ['.inp', '.in', '.yaml', '.yml', '.json', '.toml', '.cfg']

script_files = []
config_files = []

for ext in script_extensions:
    files = list(ROOT.rglob(f"*{ext}"))
    # .git 폴더 제외
    files = [f for f in files if '.git' not in str(f)]
    script_files.extend(files)

for ext in config_extensions:
    files = list(ROOT.rglob(f"*{ext}"))
    files = [f for f in files if '.git' not in str(f)]
    config_files.extend(files)

print(f"스크립트 파일: {len(script_files)}개")
print(f"  - Python (.py): {len([f for f in script_files if f.suffix == '.py'])}개")
print(f"  - Shell (.sh): {len([f for f in script_files if f.suffix == '.sh'])}개")
print(f"  - Markdown (.md): {len([f for f in script_files if f.suffix == '.md'])}개")
print()
print(f"설정 파일: {len(config_files)}개")
print()

# 2. 총 크기 계산
print("2. 파일 크기 분석")
print("-" * 80)

total_size = 0
for f in script_files + config_files:
    if f.exists():
        total_size += f.stat().st_size

print(f"스크립트 + 설정 파일 총 크기: {total_size / 1024 / 1024:.2f} MB")
print()

# 3. 제외될 대용량 파일
print("3. 제외될 대용량 파일 (100MB 이상)")
print("-" * 80)

large_files = []
for f in ROOT.rglob("*"):
    if f.is_file() and '.git' not in str(f):
        size = f.stat().st_size
        if size > 100 * 1024 * 1024:  # 100MB
            large_files.append((f, size))

large_files.sort(key=lambda x: x[1], reverse=True)

print(f"총 {len(large_files)}개 파일")
for f, size in large_files[:10]:
    rel_path = f.relative_to(ROOT)
    print(f"  {size / 1024 / 1024 / 1024:.2f} GB - {rel_path}")

print()

# 4. Git 상태 확인 (git이 초기화되어 있다면)
print("4. Git 상태")
print("-" * 80)

os.chdir(ROOT)

try:
    # git status
    result = subprocess.run(
        ["git", "status", "--short"],
        capture_output=True,
        text=True,
        check=True
    )
    
    if result.stdout:
        lines = result.stdout.strip().split('\n')
        print(f"변경된 파일: {len(lines)}개")
        print()
        print("최근 변경 파일 (최대 20개):")
        for line in lines[:20]:
            print(f"  {line}")
    else:
        print("변경사항 없음 (모두 커밋됨)")
    
except subprocess.CalledProcessError:
    print("Git이 초기화되지 않았거나 오류 발생")

print()

# 5. 권장사항
print("=" * 80)
print("권장 Git 명령어")
print("=" * 80)
print()
print("1. 추가될 파일 미리보기:")
print("   git add -n .")
print()
print("2. 스크립트만 추가:")
print("   git add '*.py' '*.sh' '*.md'")
print()
print("3. 모든 추적 대상 추가:")
print("   git add .")
print()
print("4. 커밋:")
print("   git commit -m 'Initial commit: Add scripts and configs'")
print()
print("5. 원격 저장소 연결 및 Push:")
print("   git remote add origin <repository-url>")
print("   git push -u origin main")
print()
print("=" * 80)
