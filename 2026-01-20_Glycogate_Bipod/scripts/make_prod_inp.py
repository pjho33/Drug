#!/usr/bin/env python3
import re
from pathlib import Path

BASE_INP = Path("step5_production.inp")  # CHARMM-GUI 기본 production inp
DT_PS = 0.002  # 2 fs = 0.002 ps (보통 CHARMM-GUI 기본값)
SAVE_EVERY_PS_10NS = 2.0    # 10ns 테스트는 2 ps 간격 저장(프레임 수 충분)
SAVE_EVERY_PS_100NS = 10.0  # 100ns는 10 ps 간격 저장(파일 관리)

def ns_to_nsteps(ns: float, dt_ps: float) -> int:
    total_ps = ns * 1000.0
    return int(round(total_ps / dt_ps))

def ps_to_steps(ps: float, dt_ps: float) -> int:
    return int(round(ps / dt_ps))

def patch_inp(text: str, nsteps: int, dcd_freq_steps: int, rst_freq_steps: int, out_prefix: str):
    """
    CHARMM-GUI inp 포맷은 케이스가 다양해서, 보통 존재하는 키들을 찾아 바꿉니다.
    - nsteps (또는 nstep)
    - dcdfreq / dcdFreq / DCDfreq 류
    - restartfreq / rstfreq 류
    - outputName / prefix 류가 있으면 바꿈(없으면 무시)
    """
    def sub_key(pattern, repl):
        nonlocal text
        text_new, n = re.subn(pattern, repl, text, flags=re.IGNORECASE | re.MULTILINE)
        if n > 0:
            text = text_new
        return n

    # 1) nsteps
    sub_key(r'(^\s*nsteps\s*=\s*)\d+(\s*$)', rf'\g<1>{nsteps}\2')
    sub_key(r'(^\s*nstep\s*=\s*)\d+(\s*$)',  rf'\g<1>{nsteps}\2')

    # 2) DCD frequency
    sub_key(r'(^\s*dcdfreq\s*=\s*)\d+(\s*$)', rf'\g<1>{dcd_freq_steps}\2')
    sub_key(r'(^\s*dcdFreq\s*=\s*)\d+(\s*$)', rf'\g<1>{dcd_freq_steps}\2')

    # 3) Restart / checkpoint frequency
    sub_key(r'(^\s*restartfreq\s*=\s*)\d+(\s*$)', rf'\g<1>{rst_freq_steps}\2')
    sub_key(r'(^\s*rstfreq\s*=\s*)\d+(\s*$)',     rf'\g<1>{rst_freq_steps}\2')
    sub_key(r'(^\s*checkpointfreq\s*=\s*)\d+(\s*$)', rf'\g<1>{rst_freq_steps}\2')

    # 4) output prefix/name (있으면)
    sub_key(r'(^\s*outputname\s*=\s*).+?(\s*$)', rf'\g<1>"{out_prefix}"\2')
    sub_key(r'(^\s*prefix\s*=\s*).+?(\s*$)',     rf'\g<1>"{out_prefix}"\2')

    return text

def main():
    if not BASE_INP.exists():
        raise SystemExit(f"ERROR: {BASE_INP} not found in current directory.")

    base = BASE_INP.read_text()

    # ---- 10 ns ----
    nsteps_10 = ns_to_nsteps(10.0, DT_PS)
    dcd_10 = ps_to_steps(SAVE_EVERY_PS_10NS, DT_PS)
    rst_10 = ps_to_steps(10.0, DT_PS)  # 10 ps마다 checkpoint
    out10 = patch_inp(base, nsteps_10, dcd_10, rst_10, out_prefix="prod_10ns")
    Path("step5_production_10ns.inp").write_text(out10)

    # ---- 100 ns ----
    nsteps_100 = ns_to_nsteps(100.0, DT_PS)
    dcd_100 = ps_to_steps(SAVE_EVERY_PS_100NS, DT_PS)
    rst_100 = ps_to_steps(50.0, DT_PS)  # 50 ps마다 checkpoint
    out100 = patch_inp(base, nsteps_100, dcd_100, rst_100, out_prefix="prod_100ns")
    Path("step5_production_100ns.inp").write_text(out100)

    print("Generated:")
    print(" - step5_production_10ns.inp")
    print(" - step5_production_100ns.inp")
    print()
    print("Settings:")
    print(f" dt = {DT_PS} ps (2 fs)")
    print(f" 10ns: nsteps={nsteps_10}, save every {SAVE_EVERY_PS_10NS} ps (steps={dcd_10})")
    print(f"100ns: nsteps={nsteps_100}, save every {SAVE_EVERY_PS_100NS} ps (steps={dcd_100})")

if __name__ == "__main__":
    main()
