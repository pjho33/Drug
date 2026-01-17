#!/usr/bin/env python3
"""
4단계 Restraint Minimization (클래시 제거)
Stage 1: SDG만 움직임 (restraint 50 kcal/mol·Å²)
Stage 2: restraint 10
Stage 3: restraint 1
Stage 4: 완전 자유 (restraint 0)
"""

import subprocess
from pathlib import Path
import os
import sys

RADII_PARM = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/step5_input_radii.parm7")
RST7      = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/GLUT1SDGComplex260110/amber/step5_input.rst7")
OUTDIR    = Path("/home/pjho3/projects/Drug/organized_project/12_Final_Complex_2026-01/final_complex/mmpbsa_amber/minimized_4stage")

LIG_MASK = ":SDG"          # 리간드 residue name
REST_MASK = f"!{LIG_MASK}" # SDG를 제외한 전부 restraint

STAGES = [
    ("stage1_w50", 50.0, 2000, 1000),
    ("stage2_w10", 10.0, 2000, 1000),
    ("stage3_w1",   1.0, 2000, 1000),
    ("stage4_w0",   0.0, 3000, 1500),
]

def run_cmd(cmd, cwd):
    print(f"실행: {' '.join(str(c) for c in cmd)}")
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if p.returncode != 0:
        print("❌ FAILED:", " ".join(str(c) for c in cmd))
        print("STDOUT:\n", p.stdout)
        print("STDERR:\n", p.stderr)
        sys.exit(p.returncode)
    return p

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    os.chdir(OUTDIR)

    if not RADII_PARM.exists():
        print("❌ parm not found:", RADII_PARM)
        sys.exit(1)
    if not RST7.exists():
        print("❌ rst7 not found:", RST7)
        sys.exit(1)

    print("=" * 80)
    print("4단계 Restraint Minimization")
    print("=" * 80)
    print()
    print(f"Topology: {RADII_PARM}")
    print(f"Initial coords: {RST7}")
    print(f"Output dir: {OUTDIR}")
    print(f"Ligand mask: {LIG_MASK}")
    print(f"Restraint mask: {REST_MASK}")
    print()

    current_rst = RST7

    for name, wt, maxcyc, ncyc in STAGES:
        stage_dir = OUTDIR / name
        stage_dir.mkdir(exist_ok=True)

        min_in = stage_dir / "min.in"
        
        if wt > 0:
            title = f"{name}: restrain {REST_MASK} wt={wt}"
            content = f"""{title}
 &cntrl
  imin=1,
  maxcyc={maxcyc},
  ncyc={ncyc},
  ntb=0,
  igb=8,
  cut=999.0,
  ntpr=100,
  ntr=1,
  restraint_wt={wt},
  restraintmask='{REST_MASK}',
 /
"""
        else:
            title = f"{name}: no restraints"
            content = f"""{title}
 &cntrl
  imin=1,
  maxcyc={maxcyc},
  ncyc={ncyc},
  ntb=0,
  igb=8,
  cut=999.0,
  ntpr=100,
  ntr=0,
 /
"""
        min_in.write_text(content)

        out = stage_dir / "min.out"
        rst_out = stage_dir / "min.rst7"

        cmd = [
            "sander", "-O",
            "-i", str(min_in),
            "-o", str(out),
            "-p", str(RADII_PARM),
            "-c", str(current_rst),
            "-r", str(rst_out),
        ]
        
        # restraint 사용 시 reference 좌표 필요
        if wt > 0:
            cmd.extend(["-ref", str(current_rst)])
        
        print("\n" + "="*80)
        print(f"RUN: {name}")
        print("="*80)
        
        result = run_cmd(cmd, cwd=stage_dir)
        
        # 에너지 확인
        if out.exists():
            with open(out, 'r') as f:
                lines = f.readlines()
            
            # 최종 에너지 찾기
            for i, line in enumerate(lines):
                if "FINAL RESULTS" in line:
                    print("\n최종 에너지:")
                    for j in range(i, min(i+10, len(lines))):
                        if lines[j].strip():
                            print(lines[j].strip())
                    break
        
        print(f"✅ {name} 완료")
        current_rst = rst_out

    # 최종 결과 링크
    final_link = OUTDIR / "minimized_final.rst7"
    if final_link.exists():
        final_link.unlink()
    final_link.symlink_to(current_rst.name)

    print("\n" + "="*80)
    print("✅ 4단계 Minimization 완료")
    print("="*80)
    print(f"최종 minimized restart: {current_rst}")
    print(f"Symlink: {final_link}")
    print()
    
    # 경로 저장
    with open('/tmp/minimized_final.txt', 'w') as f:
        f.write(str(final_link))

if __name__ == "__main__":
    main()
