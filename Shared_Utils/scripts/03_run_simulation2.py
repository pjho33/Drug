# scripts/03_run_simulation2.py (Final, OpenMM Template Generator Version)
import sys
import os
import time
print("🔬 [Step 3] Script starting...", flush=True)
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
        t0 = time.time()
        def _mark(stage: str):
            dt = time.time() - t0
            print(f"   ⏱️  [{dt:8.1f}s] {stage}", flush=True)

        # 1. Ligand Preparation 및 GAFF 템플릿 생성
        _mark("Preparing ligand (RDKit/OpenFF)...")
        rdkit_mol, off_mol = prepare_ligand(ligand_file)

        # For large ligands, AM1BCC can take hours. Allow overriding via env var.
        charge_method = os.environ.get('LIGAND_CHARGE_METHOD')
        if charge_method is None:
            charge_method = 'gasteiger' if off_mol.n_atoms > 150 else 'am1bcc'
        _mark(f"Assigning partial charges (method={charge_method})...")
        off_mol.assign_partial_charges(partial_charge_method=charge_method)
        
        _mark("Creating GAFFTemplateGenerator...")
        print("   ⚙️ Generating GAFF Template XML for Tripod...", flush=True)
        # OpenMM ForceFields의 Generator를 사용하여 Tripod의 템플릿을 생성
        # GAFFTemplateGenerator는 OpenFF Molecule 객체를 필요로 함
        # 캐시 파일 경로 설정 (Step 4에서 재사용)
        output_dir = os.path.dirname(output_traj)
        cache_file = os.path.join(output_dir, 'gaff_cache.json')
        gaff_generator = GAFFTemplateGenerator(molecules=[off_mol], cache=cache_file)
        print(f"   💾 GAFF cache will be saved to: {cache_file}", flush=True)
        _mark("GAFFTemplateGenerator created.")
        
        # 2. PDB 파일 병합 (Receptor + Ligand)
        _mark("Writing temporary ligand PDB and merging complex...")
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
        _mark("Loading ForceField and registering GAFF template...")
        forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
        forcefield.registerTemplateGenerator(gaff_generator.generator) # Tripod의 GAFF 파라미터 등록

        print("   🌐 Using Implicit Solvent (No water box for speed/stability)...")

        _mark("Creating OpenMM System (forcefield.createSystem)...")
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, 
                                         constraints=HBonds)
        
        print(f"   🧱 Total Atoms: {modeller.topology.getNumAtoms()} (Implicit Solvent Mode)")

        # 4. 시뮬레이션 실행 (이하 동일)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        
        try:
            platform_name = os.environ.get('OPENMM_PLATFORM', 'CUDA')
            platform = Platform.getPlatformByName(platform_name)
            prop = {'DeviceIndex': '0', 'Precision': 'mixed'}
            _mark(f"Selected OpenMM platform: {platform.getName()}")
        except Exception:
            platform = Platform.getPlatformByName('CPU')
            prop = {}
            print("   ⚠️ CUDA not available. Falling back to CPU.")
            _mark(f"Selected OpenMM platform: {platform.getName()}")

        _mark("Creating Simulation object...")
        simulation = Simulation(modeller.topology, system, integrator, platform, prop)
        simulation.context.setPositions(modeller.positions)
        
        print("   📉 Minimizing Energy (강화된 최적화)...")
        _mark("Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=5000, tolerance=10*kilojoule_per_mole/nanometer)
        
        # Equilibration: 낮은 온도에서 시작하여 점진적으로 가열
        print("   🔥 Equilibration (50K -> 300K)...")
        _mark("Equilibration...")
        simulation.context.setVelocitiesToTemperature(50*kelvin)
        for temp in [50, 100, 150, 200, 250, 300]:
            integrator.setTemperature(temp*kelvin)
            simulation.step(1000)
        
        print("   🏃 Running MD Simulation (200,000 steps)...")
        _mark("Production MD...")
        
        simulation.reporters.append(PDBReporter(output_traj, 20000))
        simulation.reporters.append(StateDataReporter(sys.stdout, 20000, step=True, potentialEnergy=True, temperature=True))
        
        simulation.step(200000)
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
