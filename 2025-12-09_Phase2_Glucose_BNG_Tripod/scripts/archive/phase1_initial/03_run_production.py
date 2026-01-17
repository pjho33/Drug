import argparse
import os
import time

from rdkit import Chem
from rdkit.Chem import AllChem

from openmm import Platform
from openmm import MonteCarloBarostat
from openmm.app import DCDReporter, ForceField, Modeller, PDBFile, Simulation, StateDataReporter, CheckpointReporter
from openmm.app import PME, HBonds
from openmm.unit import bar, kelvin, molar, nanometer, picosecond, picoseconds
from openmm import LangevinMiddleIntegrator

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator


def _prepare_ligand(ligand_file: str):
    suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
    mol = next(suppl)
    if mol is None:
        raise ValueError("Ligand could not be loaded")

    if mol.GetNumAtoms() == sum([a.GetAtomicNum() != 1 for a in mol.GetAtoms()]):
        mol = Chem.AddHs(mol, addCoords=True)

    try:
        conf = mol.GetConformer()
        if not conf.Is3D():
            raise ValueError
    except Exception:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            AllChem.UFFOptimizeMolecule(mol)

    off_mol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    return mol, off_mol


def _merge_complex(protein_pdb: str, ligand_mol) -> str:
    temp_lig_pdb = "temp_lig.pdb"
    Chem.MolToPDBFile(ligand_mol, temp_lig_pdb)

    lines_prot = [l for l in open(protein_pdb).readlines() if not l.startswith("END") and not l.startswith("CONECT")]
    lines_lig = [l for l in open(temp_lig_pdb).readlines() if l.startswith("HETATM") or l.startswith("ATOM") or l.startswith("CONECT")]

    merged = "complex_merged.pdb"
    with open(merged, "w") as f:
        f.writelines(lines_prot)
        f.writelines(lines_lig)
        f.write("END\n")

    return merged


def run_production(protein_pdb: str, ligand_sdf: str, out_prefix: str, total_ns: float, temperature_k: float, pressure_bar: float,
                   dt_fs: float, save_ps: float, platform_name: str, device_index: str, precision: str,
                   padding_nm: float, ionic_strength_m: float, equil_ns: float):
    t0 = time.time()

    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    traj_dcd = out_prefix + ".dcd"
    final_pdb = out_prefix + "_final.pdb"
    log_csv = out_prefix + "_log.csv"
    chk_path = out_prefix + ".chk"

    ligand_rdkit, ligand_off = _prepare_ligand(ligand_sdf)

    cache_file = os.path.join(out_dir if out_dir else ".", "gaff_cache.json")
    gaff = GAFFTemplateGenerator(molecules=[ligand_off], cache=cache_file)

    merged_pdb = _merge_complex(protein_pdb, ligand_rdkit)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    forcefield.registerTemplateGenerator(gaff.generator)

    pdb = PDBFile(merged_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)

    modeller.addSolvent(
        forcefield,
        padding=padding_nm * nanometer,
        ionicStrength=ionic_strength_m * molar,
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )

    system.addForce(MonteCarloBarostat(pressure_bar * bar, temperature_k * kelvin))

    integrator = LangevinMiddleIntegrator(temperature_k * kelvin, 1 / picosecond, (dt_fs / 1000.0) * picoseconds)

    try:
        platform = Platform.getPlatformByName(platform_name)
        prop = {"DeviceIndex": device_index, "Precision": precision} if platform.getName() in ("CUDA", "OpenCL") else {}
    except Exception:
        platform = Platform.getPlatformByName("CPU")
        prop = {}

    simulation = Simulation(modeller.topology, system, integrator, platform, prop)
    simulation.context.setPositions(modeller.positions)

    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(temperature_k * kelvin)

    if equil_ns < 0:
        equil_ns = 0.0
    # 1 ns = 1,000,000 fs
    eq_steps = int((equil_ns * 1_000_000.0) / dt_fs)
    prod_steps = int((total_ns * 1_000_000.0) / dt_fs)
    # 1 ps = 1000 fs
    report_steps = max(1, int((save_ps * 1000.0) / dt_fs))

    simulation.reporters.append(DCDReporter(traj_dcd, report_steps))
    simulation.reporters.append(StateDataReporter(log_csv, report_steps, step=True, potentialEnergy=True, temperature=True, volume=True))
    simulation.reporters.append(CheckpointReporter(chk_path, report_steps))

    if os.path.exists(chk_path):
        simulation.loadCheckpoint(chk_path)
    else:
        if eq_steps > 0:
            simulation.step(eq_steps)

    simulation.step(prod_steps)

    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    for f in ("temp_lig.pdb", "complex_merged.pdb"):
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError:
                pass

    dt = time.time() - t0
    print(f"done_seconds {dt:.1f}")
    print(f"traj {traj_dcd}")
    print(f"final {final_pdb}")
    print(f"log {log_csv}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("protein_pdb")
    ap.add_argument("ligand_sdf")
    ap.add_argument("out_prefix")

    ap.add_argument("--total_ns", type=float, default=100.0)
    ap.add_argument("--equil_ns", type=float, default=1.0)
    ap.add_argument("--temperature_k", type=float, default=310.0)
    ap.add_argument("--pressure_bar", type=float, default=1.0)
    ap.add_argument("--dt_fs", type=float, default=2.0)
    ap.add_argument("--save_ps", type=float, default=10.0)

    ap.add_argument("--padding_nm", type=float, default=2.0)
    ap.add_argument("--ionic_strength_m", type=float, default=0.15)

    ap.add_argument("--platform", default=os.environ.get("OPENMM_PLATFORM", "CUDA"))
    ap.add_argument("--device_index", default=os.environ.get("OPENMM_DEVICE_INDEX", "0"))
    ap.add_argument("--precision", default="mixed")

    args = ap.parse_args()

    run_production(
        protein_pdb=args.protein_pdb,
        ligand_sdf=args.ligand_sdf,
        out_prefix=args.out_prefix,
        total_ns=args.total_ns,
        temperature_k=args.temperature_k,
        pressure_bar=args.pressure_bar,
        dt_fs=args.dt_fs,
        save_ps=args.save_ps,
        platform_name=args.platform,
        device_index=args.device_index,
        precision=args.precision,
        padding_nm=args.padding_nm,
        ionic_strength_m=args.ionic_strength_m,
        equil_ns=args.equil_ns,
    )


if __name__ == "__main__":
    main()
