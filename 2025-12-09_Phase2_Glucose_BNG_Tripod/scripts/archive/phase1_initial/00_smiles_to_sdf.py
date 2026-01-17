import argparse
from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_3d_sdf(smiles: str, out_sdf: str, name: str | None = None, seed: int = 42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxIterations = 5000
    params.useRandomCoords = True

    res = AllChem.EmbedMolecule(mol, params)
    if res != 0:
        params.maxIterations = 20000
        res = AllChem.EmbedMolecule(mol, params)
        if res != 0:
            raise RuntimeError("RDKit embedding failed")

    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=5000)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=5000)

    if name:
        mol.SetProp("_Name", name)

    w = Chem.SDWriter(out_sdf)
    w.write(mol)
    w.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", help="SMILES string")
    ap.add_argument("--smiles_file", help="File containing SMILES on first line")
    ap.add_argument("--out", required=True, help="Output SDF path")
    ap.add_argument("--name", default=None)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if (args.smiles is None) == (args.smiles_file is None):
        raise SystemExit("Provide exactly one of --smiles or --smiles_file")

    smiles = args.smiles
    if args.smiles_file is not None:
        with open(args.smiles_file, "r", encoding="utf-8") as f:
            smiles = f.readline().strip()

    smiles_to_3d_sdf(smiles, args.out, name=args.name, seed=args.seed)


if __name__ == "__main__":
    main()
