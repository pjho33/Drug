# scripts/04_calc_energy.py
import sys
import os
from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
from rdkit import Chem
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule

def calculate_binding_energy(trajectory_file, ligand_file, output_log):
    print(f" [Step 4] Calculating Binding Energy (MM-GBSA): {trajectory_file}")

    # ========================================================
    # 1. 리간드 파라미터 생성 (GAFF)
    # ========================================================
    print(f"   💊 Loading Ligand Reference: {ligand_file}")
    suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
    rdkit_mol = next(suppl)
    
    if rdkit_mol.GetNumAtoms() == sum([a.GetAtomicNum() != 1 for a in rdkit_mol.GetAtoms()]):
         rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
    
    # Stereochemistry 제거 (UndefinedStereochemistryError 우회)
    for atom in rdkit_mol.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in rdkit_mol.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
         
    off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
    
    # Step 3에서 생성된 캐시 파일 재사용
    output_dir = os.path.dirname(trajectory_file)
    cache_file = os.path.join(output_dir, 'gaff_cache.json')
    if os.path.exists(cache_file):
        print(f"   ♻️ Reusing GAFF cache from: {cache_file}")
        gaff = GAFFTemplateGenerator(molecules=[off_mol], cache=cache_file)
    else:
        print(f"   ⚠️ No cache found, generating new GAFF parameters...")
        gaff = GAFFTemplateGenerator(molecules=[off_mol])

    # ========================================================
    # 2. 시스템 구축 (Implicit Solvent)
    # ========================================================
    print("   Loading trajectory...")
    pdb = PDBFile(trajectory_file)
    
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # ✅ [핵심 수정] 물(Water) + 이온(Ions) 모두 제거!
    print("   🧹 Stripping solvent (Water) and Ions (NA, CL)...")
    modeller.deleteWater()
    
    # 이온 제거 로직 추가
    ions = [r for r in modeller.topology.residues() if r.name in ['NA', 'CL', 'K', 'MG']]
    if ions:
        print(f"      - Removing {len(ions)} ions...")
        modeller.delete(ions)
    
    # ForceField 등록
    forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    forcefield.registerTemplateGenerator(gaff.generator)
    
    # 2-1. Complex System
    system_complex = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    # 2-2. Receptor & Ligand Splitting (Chain 기반)
    chains = list(modeller.topology.chains())
    print(f"   Found {len(chains)} chains in trajectory")
    for i, ch in enumerate(chains):
        n_atoms = len(list(ch.atoms()))
        print(f"      Chain {i}: {ch.id} with {n_atoms} atoms")
    
    # Receptor Only (모든 단백질 chain, 리간드 chain 제외)
    modeller_rec = Modeller(pdb.topology, pdb.positions)
    modeller_rec.deleteWater()
    ions_rec = [r for r in modeller_rec.topology.residues() if r.name in ['NA', 'CL', 'K', 'MG']]
    if ions_rec: modeller_rec.delete(ions_rec)
    
    chains_rec = list(modeller_rec.topology.chains())
    # 리간드 chain (마지막 chain) 삭제 - Trimer에서는 Chain D가 리간드
    lig_chain = [chains_rec[-1]]
    if lig_chain: modeller_rec.delete(lig_chain)
    system_rec = forcefield.createSystem(modeller_rec.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    print(f"   Receptor atoms: {system_rec.getNumParticles()}")
    
    # Ligand Only (마지막 chain = 리간드)
    modeller_lig = Modeller(pdb.topology, pdb.positions)
    modeller_lig.deleteWater()
    ions_lig = [r for r in modeller_lig.topology.residues() if r.name in ['NA', 'CL', 'K', 'MG']]
    if ions_lig: modeller_lig.delete(ions_lig)
    
    chains_lig = list(modeller_lig.topology.chains())
    # 단백질 chain들 삭제 (마지막 chain 제외)
    rec_chains = chains_lig[:-1]
    if rec_chains: modeller_lig.delete(rec_chains)
    system_lig = forcefield.createSystem(modeller_lig.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    print(f"   Ligand atoms: {system_lig.getNumParticles()}")

    # ========================================================
    # 3. 계산 준비 (Context)
    # ========================================================
    platform = Platform.getPlatformByName('CUDA')
    prop = {'DeviceIndex': '0', 'Precision': 'mixed'}
    
    int_c = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    int_r = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    int_l = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    
    ctx_complex = Context(system_complex, int_c, platform, prop)
    ctx_rec = Context(system_rec, int_r, platform, prop)
    ctx_lig = Context(system_lig, int_l, platform, prop)

    # ========================================================
    # 4. 프레임별 에너지 계산
    # ========================================================
    num_frames = pdb.getNumFrames()
    stride = max(1, num_frames // 100)
    print(f"   🏃 Processing {num_frames} frames (stride={stride})...")
    
    energies = []
    
    # PDB 파일에서 물/이온이 아닌 원자의 인덱스 찾기
    atoms = list(pdb.topology.atoms())
    valid_indices = [a.index for a in atoms if a.residue.name not in ('HOH', 'WAT', 'NA', 'CL', 'K', 'MG')]
    
    n_rec = system_rec.getNumParticles()
    n_lig = system_lig.getNumParticles()
    
    for i in range(0, num_frames, stride):
        all_pos = pdb.getPositions(frame=i)
        
        # 유효한 좌표만 추출
        complex_pos = [all_pos[j] for j in valid_indices]
        
        if len(complex_pos) != (n_rec + n_lig):
             continue

        # 1. Complex
        ctx_complex.setPositions(complex_pos)
        e_complex = ctx_complex.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        
        # 2. Receptor
        ctx_rec.setPositions(complex_pos[:n_rec])
        e_rec = ctx_rec.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        
        # 3. Ligand
        ctx_lig.setPositions(complex_pos[n_rec:])
        e_lig = ctx_lig.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        
        ddg = e_complex - (e_rec + e_lig)
        energies.append(ddg)
        
        if i % (stride*10) == 0:
            print(f"      Frame {i:4d}: dG = {ddg:.2f} kcal/mol")

    # 결과 저장
    avg_energy = np.mean(energies)
    std_energy = np.std(energies)
    
    print("-" * 50)
    print(f"   📊 Final Binding Energy: {avg_energy:.2f} +/- {std_energy:.2f} kcal/mol")
    print("-" * 50)
    
    with open(output_log, 'w') as f:
        f.write("Frame,BindingEnergy(kcal/mol)\n")
        for idx, e in enumerate(energies):
            f.write(f"{idx*stride},{e:.4f}\n")
        f.write(f"AVERAGE,{avg_energy:.4f}\n")
        f.write(f"STD_DEV,{std_energy:.4f}\n")
    
    print(f"   💾 Score saved to: {output_log}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 04_calc_energy.py <trajectory.pdb> <ligand.sdf> <output_score.csv>")
        sys.exit(1)
        
    calculate_binding_energy(sys.argv[1], sys.argv[2], sys.argv[3])