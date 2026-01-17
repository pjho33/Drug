import sys
import time
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import openmm as mm
from openmm import app
from openmm import unit

print("🧪 [1/3] RDKit: 분자 생성 및 3D 구조 최적화 시작...")

# 1. RDKit: 아스피린(Aspirin) 생성
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# 2. 3D 좌표 생성
AllChem.EmbedMolecule(mol, AllChem.ETKDG()) 
AllChem.MMFFOptimizeMolecule(mol)

print(f"   ✅ Molecule Created: Aspirin (Atoms: {mol.GetNumAtoms()})")
print("   ✅ 3D Coordinates Generated successfully.")

print("-" * 50)
print("⚙️ [2/3] OpenMM: 시뮬레이션 엔진 가동 (with RTX 3090)...")

# 3. OpenMM System Test
pdb_path = os.path.join(os.path.dirname(app.__file__), 'data', 'test.pdb')
pdb = app.PDBFile(pdb_path)

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# ✅ 수정된 부분: app.PDB -> app.NoCutoff
# (테스트용 단백질은 물 상자(Periodic Box)가 없으므로 Cutoff를 쓰지 않습니다)
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds)

integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)

# 4. 시뮬레이션 객체 생성 및 플랫폼(GPU) 설정
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    prop = {'DeviceIndex': '0', 'Precision': 'mixed'} 
    simulation = app.Simulation(pdb.topology, system, integrator, platform, prop)
    
    device_name = platform.getPropertyValue(simulation.context, 'DeviceName')
    print(f"   ✅ Platform: {platform.getName()}")
    print(f"   ✅ GPU Device: {device_name}")
    
except Exception as e:
    print(f"   ⚠️ GPU Error: {e}")
    print("   CPU로 대체합니다 (느림).")
    simulation = app.Simulation(pdb.topology, system, integrator)

print("-" * 50)
print("🚀 [3/3] Performance Test: 1000 Steps Run...")

simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

start_time = time.time()
simulation.step(1000)
end_time = time.time()

print(f"   ✅ Simulation Completed!")
print(f"   ⏱️ Time taken: {end_time - start_time:.4f} seconds")
print(f"   🎉 Congratulations! Your Pipeline is READY.")