# scripts/02_run_docking.py
import os
import subprocess
import sys

def run_diffdock(protein_path, ligand_smiles, output_dir, n_samples=10, n_steps=20):
    if os.path.isfile(ligand_smiles):
        with open(ligand_smiles, 'r', encoding='utf-8') as f:
            ligand_smiles = f.readline().strip()

    protein_path = os.path.abspath(protein_path)
    output_dir = os.path.abspath(output_dir)

    print(f"🤖 [Step 2] Running DiffDock AI Prediction (High-Performance)...")
    print(f"   Target: {protein_path}")
    print(f"   Ligand: {ligand_smiles}")
    
    # DiffDock 경로 (워크스테이션 환경에 맞게 확인 필요!)
    diffdock_home = os.path.expanduser("~/projects/Drug/DiffDock")
    
    if not os.path.exists(diffdock_home):
        print(f"❌ Error: DiffDock not found at {diffdock_home}")
        return

    # 실행 명령어 구성
    cmd = [
        "python", "-m", "inference",
        "--protein_path", protein_path,
        "--ligand_description", ligand_smiles,
        "--out_dir", output_dir,
        "--inference_steps", str(n_steps),
        "--samples_per_complex", str(n_samples),
        # ✅ [복구 완료] Batch Size를 1 -> 5로 원상복구! (속도 향상)
        "--batch_size", "5"
    ]
    
    print("   🚀 Executing DiffDock...")
    try:
        subprocess.run(cmd, cwd=diffdock_home, check=True)
        print(f"   ✅ Docking Finished! Results in: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"   ❌ Docking Failed: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 02_run_docking.py <clean_pdb> <smiles> <output_dir>")
        sys.exit(1)
        
    run_diffdock(sys.argv[1], sys.argv[2], sys.argv[3])