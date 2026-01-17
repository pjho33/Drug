# scripts/05_analyze_movement.py
import sys
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os

def analyze_movement(trajectory_file, output_plot):
    print(f"📉 [Step 5] Analyzing Ligand Movement (RMSD): {trajectory_file}")
    
    # 1. 궤적 로딩
    print("   📂 Loading trajectory (this might take a moment)...")
    try:
        t = md.load(trajectory_file)
    except OSError:
        print(f"   ❌ Error: File not found - {trajectory_file}")
        return

    # 2. 리간드 찾기
    topology = t.topology
    # 물(HOH), 이온(NA, CL 등)이 아닌 마지막 잔기를 리간드로 가정
    # (일반적으로 PDB 파일 맨 끝에 위치함)
    ligand_candidates = [r for r in topology.residues if r.name not in ('HOH', 'WAT', 'NA', 'CL', 'K', 'MG')]
    
    if not ligand_candidates:
        print("   ❌ Error: No ligand found in topology.")
        return

    ligand = ligand_candidates[-1]
    print(f"   🎯 Target Ligand identified: {ligand.name} (Residue Index: {ligand.index})")
    
    # ✅ [핵심 수정] mdtraj에서 atom index를 가져오는 올바른 방법
    ligand_indices = [atom.index for atom in ligand.atoms]
    
    # 3. RMSD 계산
    # 단백질 Backbone을 기준으로 정렬(Superpose)
    prot_indices = topology.select("protein and backbone")
    t.superpose(t, frame=0, atom_indices=prot_indices)
    
    # 리간드 RMSD 계산 (nm -> Angstrom 변환)
    rmsd = md.rmsd(t, t, frame=0, atom_indices=ligand_indices) * 10.0
    
    # 4. 통계 출력
    avg_rmsd = np.mean(rmsd)
    std_rmsd = np.std(rmsd)
    
    print("-" * 40)
    print(f"   📊 Ligand RMSD Statistics:")
    print(f"      - Average Movement: {avg_rmsd:.2f} Å")
    print(f"      - Fluctuation (Std): {std_rmsd:.2f} Å")
    print("-" * 40)
    
    if avg_rmsd < 2.0:
        print("   🔒 Interpretation: [Stable] 꽉 끼어있거나 안정적입니다.")
    elif avg_rmsd > 5.0:
        print("   💃 Interpretation: [Unstable] 포켓 밖으로 나가거나, 심하게 움직입니다.")
    else:
        print("   ⚖️ Interpretation: [Moderate] 포켓 안에서 자리를 잡으려 움직입니다.")

    # 5. 그래프 그리기
    plt.figure(figsize=(10, 6))
    plt.plot(rmsd, label=f'Ligand ({ligand.name})', color='royalblue', linewidth=1.5)
    plt.axhline(y=avg_rmsd, color='red', linestyle='--', label=f'Avg: {avg_rmsd:.2f} Å')
    
    plt.title(f"Ligand Stability Analysis (RMSD)\nSource: {os.path.basename(trajectory_file)}")
    plt.xlabel("Simulation Time (Frames)")
    plt.ylabel("RMSD (Å)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(output_plot)
    print(f"   📈 Plot saved to: {output_plot}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 05_analyze_movement.py <trajectory.pdb> <output_plot.png>")
        sys.exit(1)
    
    analyze_movement(sys.argv[1], sys.argv[2])