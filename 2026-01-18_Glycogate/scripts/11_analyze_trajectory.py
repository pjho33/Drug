#!/usr/bin/env python3
"""
1-Arm PEG24-Glc MD Trajectory 분석

분석 항목:
1. End-to-end 거리 (TRIS 중심 - Glucose 말단)
2. Radius of gyration
3. 시간에 따른 변화
4. 분포 분석
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_trajectory(results_dir, replica=1):
    """Trajectory 분석"""
    
    print("=" * 80)
    print("1-Arm PEG24-Glc Trajectory 분석")
    print("=" * 80)
    print()
    
    results_path = Path(results_dir)
    
    # MDTraj 사용
    try:
        import mdtraj as md
        print("✅ MDTraj 사용")
        use_mdtraj = True
    except ImportError:
        print("⚠️  MDTraj 없음, MDAnalysis 시도")
        use_mdtraj = False
    
    if not use_mdtraj:
        try:
            import MDAnalysis as mda
            print("✅ MDAnalysis 사용")
        except ImportError:
            print("❌ MDTraj 또는 MDAnalysis 설치 필요")
            print("   conda install -c conda-forge mdtraj")
            return
    
    print()
    
    # Topology 파일 찾기
    topology_file = None
    for topo in [
        "/home/pjho3/projects/Drug/2026-01-18_Glycogate/data/solution builder/openmm/step3_input.pdb",
        results_path / "step4_equilibration_rep1.rst"
    ]:
        if Path(topo).exists():
            topology_file = str(topo)
            break
    
    if not topology_file:
        print("❌ Topology 파일을 찾을 수 없습니다.")
        return
    
    print(f"📁 Topology: {topology_file}")
    
    # Trajectory 파일 찾기
    traj_files = sorted(results_path.glob(f"step5_*_rep{replica}.dcd"))
    
    if not traj_files:
        print("❌ Trajectory 파일이 없습니다.")
        return
    
    print(f"📁 Trajectory 파일: {len(traj_files)}개")
    print()
    
    # 분석 시작
    if use_mdtraj:
        analyze_with_mdtraj(topology_file, traj_files, results_path, replica)
    else:
        analyze_with_mdanalysis(topology_file, traj_files, results_path, replica)


def analyze_with_mdtraj(topology_file, traj_files, results_path, replica):
    """MDTraj을 사용한 분석"""
    
    import mdtraj as md
    
    print("Step 1: Trajectory 로드")
    print("-" * 80)
    
    # 첫 번째 파일로 topology 로드
    print(f"로딩 중: {traj_files[0]}")
    traj = md.load(str(traj_files[0]), top=topology_file)
    
    # 나머지 파일 병합
    for traj_file in traj_files[1:]:
        print(f"로딩 중: {traj_file}")
        t = md.load(str(traj_file), top=topology_file)
        traj = traj.join(t)
    
    print(f"✅ 총 프레임: {traj.n_frames}")
    print(f"   시간: {traj.n_frames * 0.1:.1f} ns (100 ps 간격)")
    print()
    
    # Ligand 선택
    print("Step 2: Ligand 원자 선택")
    print("-" * 80)
    
    # Ligand residue 찾기
    ligand = traj.topology.select('resname LIG or resname UNK or resname MOL')
    
    if len(ligand) == 0:
        print("⚠️  Ligand를 찾을 수 없습니다. 모든 원자 사용")
        ligand = traj.topology.select('all')
    
    print(f"✅ Ligand 원자: {len(ligand)}개")
    
    # Ligand만 추출
    traj_lig = traj.atom_slice(ligand)
    print()
    
    # 분석 1: Radius of gyration
    print("Step 3: Radius of Gyration 분석")
    print("-" * 80)
    
    rg = md.compute_rg(traj_lig)
    
    print(f"평균 Rg: {np.mean(rg) * 10:.2f} ± {np.std(rg) * 10:.2f} Å")
    print(f"최소 Rg: {np.min(rg) * 10:.2f} Å")
    print(f"최대 Rg: {np.max(rg) * 10:.2f} Å")
    
    # 저장
    np.save(results_path / f"rg_rep{replica}.npy", rg)
    print(f"✅ 저장: rg_rep{replica}.npy")
    print()
    
    # 분석 2: End-to-end 거리
    print("Step 4: End-to-End 거리 분석")
    print("-" * 80)
    
    # TRIS 중심 원자 찾기 (첫 번째 탄소)
    try:
        tris_atoms = traj_lig.topology.select('name C1 or name C2 or name C3')
        if len(tris_atoms) == 0:
            tris_atoms = [0]  # 첫 번째 원자
        tris_center = tris_atoms[0]
        
        # Glucose 말단 원자 찾기 (마지막 산소)
        glc_atoms = traj_lig.topology.select('name O6 or name O5 or name O4')
        if len(glc_atoms) == 0:
            glc_atoms = [traj_lig.n_atoms - 1]  # 마지막 원자
        glc_end = glc_atoms[-1]
        
        print(f"TRIS 중심: atom {tris_center}")
        print(f"Glucose 말단: atom {glc_end}")
        
        # 거리 계산
        distances = md.compute_distances(traj_lig, [[tris_center, glc_end]])
        distances = distances.flatten()
        
        print(f"평균 거리: {np.mean(distances) * 10:.2f} ± {np.std(distances) * 10:.2f} Å")
        print(f"최소 거리: {np.min(distances) * 10:.2f} Å")
        print(f"최대 거리: {np.max(distances) * 10:.2f} Å")
        
        # 저장
        np.save(results_path / f"end_to_end_rep{replica}.npy", distances)
        print(f"✅ 저장: end_to_end_rep{replica}.npy")
        
    except Exception as e:
        print(f"⚠️  End-to-end 거리 계산 실패: {e}")
        distances = None
    
    print()
    
    # 시각화
    print("Step 5: 결과 시각화")
    print("-" * 80)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Rg 시계열
    axes[0, 0].plot(np.arange(len(rg)) * 0.1, rg * 10, linewidth=0.5)
    axes[0, 0].set_xlabel('Time (ns)')
    axes[0, 0].set_ylabel('Radius of Gyration (Å)')
    axes[0, 0].set_title('Rg vs Time')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Rg 분포
    axes[0, 1].hist(rg * 10, bins=50, density=True, alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('Radius of Gyration (Å)')
    axes[0, 1].set_ylabel('Probability Density')
    axes[0, 1].set_title('Rg Distribution')
    axes[0, 1].axvline(np.mean(rg) * 10, color='red', linestyle='--', label='Mean')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    if distances is not None:
        # End-to-end 시계열
        axes[1, 0].plot(np.arange(len(distances)) * 0.1, distances * 10, linewidth=0.5)
        axes[1, 0].set_xlabel('Time (ns)')
        axes[1, 0].set_ylabel('End-to-End Distance (Å)')
        axes[1, 0].set_title('End-to-End Distance vs Time')
        axes[1, 0].grid(True, alpha=0.3)
        
        # End-to-end 분포
        axes[1, 1].hist(distances * 10, bins=50, density=True, alpha=0.7, edgecolor='black')
        axes[1, 1].set_xlabel('End-to-End Distance (Å)')
        axes[1, 1].set_ylabel('Probability Density')
        axes[1, 1].set_title('End-to-End Distance Distribution')
        axes[1, 1].axvline(np.mean(distances) * 10, color='red', linestyle='--', label='Mean')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(results_path / f"analysis_rep{replica}.png", dpi=300)
    print(f"✅ 저장: analysis_rep{replica}.png")
    print()
    
    # 요약 저장
    print("Step 6: 요약 저장")
    print("-" * 80)
    
    summary = {
        'n_frames': traj.n_frames,
        'time_ns': traj.n_frames * 0.1,
        'rg_mean': float(np.mean(rg) * 10),
        'rg_std': float(np.std(rg) * 10),
        'rg_min': float(np.min(rg) * 10),
        'rg_max': float(np.max(rg) * 10),
    }
    
    if distances is not None:
        summary.update({
            'end_to_end_mean': float(np.mean(distances) * 10),
            'end_to_end_std': float(np.std(distances) * 10),
            'end_to_end_min': float(np.min(distances) * 10),
            'end_to_end_max': float(np.max(distances) * 10),
        })
    
    import json
    with open(results_path / f"summary_rep{replica}.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"✅ 저장: summary_rep{replica}.json")
    print()
    
    print("=" * 80)
    print("✅ 분석 완료!")
    print("=" * 80)


def analyze_with_mdanalysis(topology_file, traj_files, results_path, replica):
    """MDAnalysis를 사용한 분석"""
    
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, distances
    
    print("MDAnalysis 분석은 아직 구현되지 않았습니다.")
    print("MDTraj을 설치해주세요: conda install -c conda-forge mdtraj")


def main():
    """메인 함수"""
    
    results_dir = "/home/pjho3/projects/Drug/2026-01-18_Glycogate/results/md_1arm_openmm"
    replica = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    
    analyze_trajectory(results_dir, replica)


if __name__ == "__main__":
    main()
