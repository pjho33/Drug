import sys
import os
import time
import numpy as np
from rdkit import Chem
from openmm import *
from openmm.app import *
from openmm.unit import *
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule
from pdbfixer import PDBFixer

# ==========================================
# 1. 설정 (Configuration)
# ==========================================
protein_file = '1a0q.pdb'
ligand_file = 'rank1.sdf'
output_traj = 'validation_run.pdb' 

print("🚀 [Validation] 검증 파이프라인 시작 (Cleaning Heterogens)...")

# ==========================================
# 2. 단백질 청소 및 수리 (PDBFixer)
# ==========================================
print(f"📦 Fixing Protein Structure: {protein_file}")
fixer = PDBFixer(filename=protein_file)

# ✅ [핵심 수정] 잡동사니(HEP 등) 제거!
# keepWater=False: 결정화 물분자도 제거하고, 나중에 깨끗한 물로 다시 채웁니다.
print("   🧹 Removing artifacts (HEP, waters, ions)...")
fixer.removeHeterogens(keepWater=False) 

# 1. 빠진 잔기 찾기
fixer.findMissingResidues()
# 2. 빠진 원자(Side chain) 복구
fixer.findMissingAtoms()
print("   👉 Adding missing heavy atoms...")
fixer.addMissingAtoms()
# 3. 수소 추가 (pH 7.0)
print("   👉 Adding Hydrogens (pH 7.0)...")
fixer.addMissingHydrogens(7.0)

# 수리된 단백질 저장
print("   ✅ Protein Cleaned & Fixed! Saving to 'fixed_protein.pdb'")
PDBFile.writeFile(fixer.topology, fixer.positions, open('fixed_protein.pdb', 'w'))

# ==========================================
# 3. 리간드 처리 및 합치기
# ==========================================
print(f"📦 Processing Ligand: {ligand_file}")

forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# 리간드 로딩
suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
rdkit_mol = next(suppl)

if rdkit_mol.GetNumAtoms() == sum([a.GetAtomicNum() != 1 for a in rdkit_mol.GetAtoms()]):
     print("   ⚠️ 리간드에 수소가 부족해 보입니다. 추가합니다.")
     rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)

print("   👉 Generating Ligand Parameters (OpenFF)...")
off_mol = Molecule.from_rdkit(rdkit_mol)
gaff = GAFFTemplateGenerator(molecules=[off_mol])
forcefield.registerTemplateGenerator(gaff.generator)

# 병합
print("   👉 Merging Protein and Ligand...")
Chem.MolToPDBFile(rdkit_mol, 'temp_lig.pdb')

lines_prot = [l for l in open('fixed_protein.pdb').readlines() if not l.startswith('END') and not l.startswith('CONECT')]
lines_lig = [l for l in open('temp_lig.pdb').readlines() if l.startswith('HETATM') or l.startswith('CONECT')]

with open('complex_final.pdb', 'w') as f_out:
    f_out.writelines(lines_prot)
    f_out.writelines(lines_lig)
    f_out.write("END\n")

# ==========================================
# 4. 시스템 구축 및 시뮬레이션
# ==========================================
pdb_complex = PDBFile('complex_final.pdb')
modeller = Modeller(pdb_complex.topology, pdb_complex.positions)

print("   💧 Adding Water Box (Solvation)...")
# 이제 HEP이 없으므로 에러가 나지 않습니다!
modeller.addSolvent(forcefield, padding=1.0*nanometer, ionicStrength=0.15*molar)

print("   🔌 Creating OpenMM System...")
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, 
                                 nonbondedCutoff=1.0*nanometer, constraints=HBonds)

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')
prop = {'DeviceIndex': '0', 'Precision': 'mixed'}

simulation = Simulation(modeller.topology, system, integrator, platform, prop)
simulation.context.setPositions(modeller.positions)

print(f"   ✅ System Created on GPU: {platform.getPropertyValue(simulation.context, 'DeviceName')}")
print(f"   🧱 Total Atoms: {modeller.topology.getNumAtoms()}")

print("📉 Minimizing Energy...")
initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"   Initial Energy: {initial_energy.value_in_unit(kilojoules_per_mole):.2f} kJ/mol")

simulation.minimizeEnergy()
final_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"   Final Energy:   {final_energy.value_in_unit(kilojoules_per_mole):.2f} kJ/mol")

print("🏃 Running Stability Test (50,000 steps)...")
simulation.reporters.append(PDBReporter(output_traj, 1000))
simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True, speed=True))

start_time = time.time()
simulation.step(50000)
end_time = time.time()

print("-" * 50)
print(f"🎉 Validation Finished!")
print(f"⏱️ Time taken: {end_time - start_time:.2f} seconds")
print(f"📂 Result saved to: {output_traj}")