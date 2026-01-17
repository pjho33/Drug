# scripts/03_run_simulation2.py (Final, OpenMM Template Generator Version)
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from openmm import *
from openmm.app import *
from openmm.unit import *
from openff.toolkit.topology import Molecule as OFFMolecule
from openmmforcefields.generators import GAFFTemplateGenerator # 템플릿 생성기 사용

# ==============================================================================
# 📌 Force Field 로드는 OpenMM 표준 파일명만 사용
# ==============================================================================

def prepare_ligand(ligand_file):
    """Load ligand and return both RDKit mol (for PDB) and OpenFF Molecule (for GAFF)."""
    suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
    rdkit_mol = next(suppl)
    
    if rdkit_mol is None:
        raise ValueError("Ligand could not be loaded by RDKit.")

    # UndefinedStereochemistryError 우회 - 모든 stereochemistry 제거
    for atom in rdkit_mol.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in rdkit_mol.GetBonds():
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            
    # 3D 좌표가 없으면 생성 및 간단한 MM 최적화
    if rdkit_mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(rdkit_mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(rdkit_mol)

    # OpenFF Molecule 생성 (GAFFTemplateGenerator용)
    # allow_undefined_stereo=True로 stereochemistry 에러 우회
    off_mol = OFFMolecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)

    return rdkit_mol, off_mol

def run_validation(protein_file, ligand_file, output_traj):
    print(f"🔬 [Step 3] Running MD Validation (OpenMM GAFF Template): {protein_file} + {ligand_file}")

    temp_lig_pdb = 'temp_lig.pdb'
    cleanup_files = ['complex_merged.pdb', temp_lig_pdb]

    try:
        # 1. Ligand Preparation 및 GAFF 템플릿 생성
        rdkit_mol, off_mol = prepare_ligand(ligand_file)
        
        print("   ⚙️ Generating GAFF Template XML for Tripod...")
        # OpenMM ForceFields의 Generator를 사용하여 Tripod의 템플릿을 생성
        # GAFFTemplateGenerator는 OpenFF Molecule 객체를 필요로 함
        # 캐시 파일 경로 설정 (Step 4에서 재사용)
        output_dir = os.path.dirname(output_traj)
        cache_file = os.path.join(output_dir, 'gaff_cache.json')
        gaff_generator = GAFFTemplateGenerator(molecules=[off_mol], cache=cache_file)
        print(f"   💾 GAFF cache will be saved to: {cache_file}")
        
        # 2. PDB 파일 병합 (Receptor + Ligand)
        Chem.MolToPDBFile(rdkit_mol, temp_lig_pdb)
        
        lines_prot = [l for l in open(protein_file).readlines() if not l.startswith('END') and not l.startswith('CONECT')]
        lines_lig = [l for l in open(temp_lig_pdb).readlines()
                     if l.startswith('HETATM') or l.startswith('ATOM') or l.startswith('CONECT')]
        
        with open('complex_merged.pdb', 'w') as f:
            f.writelines(lines_prot)
            f.writelines(lines_lig)
            f.write("END\n")
            
        pdb_complex = PDBFile('complex_merged.pdb')
        modeller = Modeller(pdb_complex.topology, pdb_complex.positions)
        
        # 3. Force Field 설정 및 시스템 생성
        # amber14-all.xml만 로드하고, GAFFTemplateGenerator가 리간드 파라미터를 자동 생성
        forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
        forcefield.registerTemplateGenerator(gaff_generator.generator) # Tripod의 GAFF 파라미터 등록

        print("   🌐 Using Implicit Solvent (No water box for speed/stability)...")
        
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, 
                                         constraints=HBonds)
        
        print(f"   🧱 Total Atoms: {modeller.topology.getNumAtoms()} (Implicit Solvent Mode)")

        # 4. 시뮬레이션 실행 (이하 동일)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        
        try:
            platform = Platform.getPlatformByName('CUDA')
            prop = {'DeviceIndex': '0', 'Precision': 'mixed'}
        except Exception:
            platform = Platform.getPlatformByName('CPU')
            prop = {}
            print("   ⚠️ CUDA not available. Falling back to CPU.")

        simulation = Simulation(modeller.topology, system, integrator, platform, prop)
        simulation.context.setPositions(modeller.positions)
        
        print("   📉 Minimizing Energy...")
        simulation.minimizeEnergy(maxIterations=1000)
        
        print("   🏃 Running MD Simulation (50,000 steps)...")
        
        simulation.reporters.append(PDBReporter(output_traj, 5000))
        simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True))
        
        simulation.step(50000)
        print(f"   ✅ Validation Done. Trajectory saved to {output_traj}")
    
    except Exception as e:
        print(f"\nFATAL ERROR during MD Simulation: {e}")
        sys.exit(1)
        
    finally:
        # 5. Cleanup
        for f in cleanup_files:
            if os.path.exists(f): 
                try:
                    os.remove(f)
                except OSError:
                    pass


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 03_run_simulation2.py <clean_protein> <ligand_file> <output_traj>")
        sys.exit(1)
    
    run_validation(sys.argv[1], sys.argv[2], sys.argv[3])
