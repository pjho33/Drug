# scripts/03_run_simulation.py
import sys
import os
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from openmm import *
from openmm.app import *
from openmm.unit import *
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule

def run_validation(protein_file, ligand_file, output_traj):
    print(f"🔬 [Step 3] Running MD Validation (Explicit Solvent / High-End): {protein_file} + {ligand_file}")

    # ======================================================
    # 1. 리간드 준비
    # ======================================================
    suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
    rdkit_mol = next(suppl)

    if rdkit_mol.GetNumAtoms() == sum([a.GetAtomicNum() != 1 for a in rdkit_mol.GetAtoms()]):
         print("   ⚠️ No hydrogens found. Adding hydrogens...")
         rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
    
    try:
        conf = rdkit_mol.GetConformer()
        if not conf.Is3D(): raise ValueError
    except:
        print("   ⚠️ 2D Detected -> Embedding 3D...")
        AllChem.EmbedMolecule(rdkit_mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(rdkit_mol)

    # ======================================================
    # 2. 시스템 구축 (Explicit Solvent: 진짜 물!)
    # ======================================================
    off_mol = Molecule.from_rdkit(rdkit_mol)
    gaff = GAFFTemplateGenerator(molecules=[off_mol])

    # ✅ 원래대로 복귀: amber14/tip3p.xml (TIP3P 물 모델 사용)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    forcefield.registerTemplateGenerator(gaff.generator)
    
    pdb_prot = PDBFile(protein_file)
    
    Chem.MolToPDBFile(rdkit_mol, 'temp_lig.pdb')
    lines_prot = [l for l in open(protein_file).readlines() if not l.startswith('END') and not l.startswith('CONECT')]
    lines_lig = [l for l in open('temp_lig.pdb').readlines() if l.startswith('HETATM') or l.startswith('CONECT')]
    
    with open('complex_merged.pdb', 'w') as f:
        f.writelines(lines_prot)
        f.writelines(lines_lig)
        f.write("END\n")
        
    pdb_complex = PDBFile('complex_merged.pdb')
    modeller = Modeller(pdb_complex.topology, pdb_complex.positions)
    
    print("   💧 Adding Water Box (Explicit Solvent)...")
    # ✅ [중요] 700GB 램이 있으니 여유 있게 잡으세요.
    # 일반 약물: 1.0 ~ 1.2 nm
    # Tripod: 2.0 ~ 3.0 nm (여기를 고치면 됩니다!)
    modeller.addSolvent(forcefield, padding=2.0*nanometer, ionicStrength=0.15*molar)
    
    print(f"   🧱 Total Atoms: {modeller.topology.getNumAtoms()} (Ready for High-RAM)")

    # PME (Particle Mesh Ewald) 사용 - 정밀 계산
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, 
                                     nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    
    # ======================================================
    # 3. 시뮬레이션 실행
    # ======================================================
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    platform = Platform.getPlatformByName('CUDA')
    # 700GB 머신에 GPU도 좋다면 'mixed' 사용. 만약 GPU가 옛날 거라면 'single' 고려.
    prop = {'DeviceIndex': '0', 'Precision': 'mixed'}
    
    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)
    
    print("   📉 Minimizing Energy...")
    simulation.minimizeEnergy()
    
    print("   🏃 Running MD Simulation (50,000 steps)...")
    
    # 저장 간격: 메모리가 많아도 파일 관리 편의성을 위해 5000 유지 추천
    simulation.reporters.append(PDBReporter(output_traj, 5000))
    simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True))
    
    simulation.step(50000)
    print(f"   ✅ Validation Done. Trajectory saved to {output_traj}")
    
    if os.path.exists('temp_lig.pdb'): os.remove('temp_lig.pdb')
    if os.path.exists('complex_merged.pdb'): os.remove('complex_merged.pdb')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 03_run_simulation.py <clean_protein> <ligand_file> <output_traj>")
        sys.exit(1)
    
    run_validation(sys.argv[1], sys.argv[2], sys.argv[3])