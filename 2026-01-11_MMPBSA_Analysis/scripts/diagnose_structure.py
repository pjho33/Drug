#!/usr/bin/env python3
"""
구조 진단: PBC 이미징 및 클래시 확인
"""

import subprocess
from pathlib import Path

AMBER_DIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber")
RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")
RST7 = AMBER_DIR / "step5_input.rst7"
OUTDIR = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/diagnosis")

OUTDIR.mkdir(exist_ok=True)

print("=" * 80)
print("구조 진단: PBC 이미징 및 클래시 확인")
print("=" * 80)
print()

# Step 1: 리간드-단백질 최소 거리 확인 (클래시 체크)
print("Step 1: 리간드-단백질 최소 거리 확인")
print("-" * 80)

cpptraj_input = f"""parm {RADII_PARM}
trajin {RST7}

# 리간드-단백질 최소 거리 (원자 간)
nativecontacts :306 :1-305,307-806 mindist out {OUTDIR}/min_distance.dat

# 리간드-단백질 중심 거리 (질량중심 간)
distance center_dist :306 :1-305,307-806 out {OUTDIR}/center_distance.dat

run
"""

with open(OUTDIR / "check_distance.in", "w") as f:
    f.write(cpptraj_input)

print("cpptraj 실행 중...")
result = subprocess.run(
    ["cpptraj", "-i", str(OUTDIR / "check_distance.in")],
    capture_output=True,
    text=True
)

print(result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr)
    exit(1)

print()

# 결과 읽기
print("=" * 80)
print("진단 결과")
print("=" * 80)
print()

# 최소 거리
min_dist_file = OUTDIR / "min_distance.dat"
if min_dist_file.exists():
    with open(min_dist_file, 'r') as f:
        lines = f.readlines()
        if len(lines) > 1:
            data = lines[1].split()
            min_dist = float(data[1])
            print(f"리간드-단백질 최소 거리: {min_dist:.2f} Å")
            print()
            
            if min_dist < 1.0:
                print("❌ 심각한 클래시 감지! (< 1.0 Å)")
                print("   → Minimization 필요")
            elif min_dist < 2.0:
                print("⚠️  클래시 가능성 있음 (< 2.0 Å)")
                print("   → 구조 확인 권장")
            else:
                print("✅ 클래시 없음 (≥ 2.0 Å)")

# 중심 거리
center_dist_file = OUTDIR / "center_distance.dat"
if center_dist_file.exists():
    with open(center_dist_file, 'r') as f:
        lines = f.readlines()
        if len(lines) > 1:
            data = lines[1].split()
            center_dist = float(data[1])
            print()
            print(f"리간드-단백질 중심 거리: {center_dist:.2f} Å")
            print()
            
            if center_dist > 50.0:
                print("❌ PBC 이미징 문제 가능성 높음! (> 50 Å)")
                print("   → autoimage/center 필요")
            elif center_dist > 30.0:
                print("⚠️  거리가 멀어 보임 (> 30 Å)")
                print("   → 구조 시각화 확인 권장")
            else:
                print("✅ 정상 범위 (≤ 30 Å)")

print()
print("=" * 80)

# Step 2: PDB 생성 (시각화용)
print()
print("Step 2: 시각화용 PDB 생성")
print("-" * 80)

cpptraj_pdb = f"""parm {RADII_PARM}
trajin {RST7}
strip :POT,CLA,WAT
trajout {OUTDIR}/complex_visual.pdb pdb
run
"""

with open(OUTDIR / "make_pdb.in", "w") as f:
    f.write(cpptraj_pdb)

result = subprocess.run(
    ["cpptraj", "-i", str(OUTDIR / "make_pdb.in")],
    capture_output=True,
    text=True
)

if result.returncode == 0:
    print(f"✅ PDB 생성: {OUTDIR}/complex_visual.pdb")
    print()
    print("VMD/PyMOL로 확인:")
    print(f"  vmd {OUTDIR}/complex_visual.pdb")
    print()
    print("확인 사항:")
    print("  1. 리간드가 단백질 결합부위에 정상적으로 위치하는가?")
    print("  2. 리간드와 단백질이 박스 반대편으로 분리되어 있지 않은가?")
    print("  3. 원자들이 겹쳐 파고들어 있지 않은가?")
else:
    print("❌ PDB 생성 실패")
    print(result.stderr)

print()
print("=" * 80)
print(f"진단 완료. 결과 파일: {OUTDIR}")
print("=" * 80)
