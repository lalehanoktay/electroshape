import argparse
import gzip
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List

from rdkit import Chem, rdBase
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

from .electroshape import prepare_mol, process_record

def main():
    rdBase.DisableLog('rdApp.*')
    
    ap = argparse.ArgumentParser(description="ElectroShape 5D Generator (Armstrong 2011)")
    ap.add_argument("--sdf", required=True, help="Input SDF file (.sdf or .sdf.gz)")
    ap.add_argument("--out", required=True, help="Output CSV file (.csv or .csv.gz)")
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1), help="Number of CPU cores to use.")
    ap.add_argument("--chunk", type=int, default=5000, help="Number of molecules to process in a batch.")
    args = ap.parse_args()

    # Open file (handle gzip)
    if args.sdf.endswith(".gz"):
        suppl = Chem.ForwardSDMolSupplier(gzip.open(args.sdf, 'rb'))
    else:
        suppl = Chem.ForwardSDMolSupplier(args.sdf)

    is_gz = args.out.endswith(".gz")
    f_out = gzip.open(args.out, "wt") if is_gz else open(args.out, "wt")
    
    header = ["IDNUMBER", "Canonical_SMILES", "n_atoms", "charge_source"] + [f"ES_{i+1}" for i in range(18)]
    f_out.write(",".join(header) + "\n")

    tasks = []
    
    print(f"Reading molecules from {args.sdf}...")
    for mol in suppl:
        if mol is not None:
            res = prepare_mol(mol)
            if res:
                tasks.append((res[0], res[1]))
            
    print(f"Queued {len(tasks)} molecules for processing on {args.workers} cores...")
    
    buffer = []
    processed_count = 0
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_task = {executor.submit(process_record, t): t for t in tasks}
        
        pbar_kwargs = {"total": len(tasks), "unit": "mol", "desc": "Generating ElectroShape descriptors"}
        pbar = tqdm(**pbar_kwargs) if HAS_TQDM else None
        
        for future in as_completed(future_to_task):
            res = future.result()
            if pbar: pbar.update(1)
            
            if res:
                processed_count += 1
                buffer.append(res)
                
                if len(buffer) >= args.chunk:
                    for item in buffer:
                        vec_str = ",".join([f"{x:.6f}" for x in item["vector"]])
                        line = f"{item['mol_id']},{item['smiles']},{item['n_atoms']},{item['charge_source']},{vec_str}"
                        f_out.write(line + "\n")
                    buffer = []
        
        if pbar: pbar.close()

    if buffer:
        for item in buffer:
            vec_str = ",".join([f"{x:.6f}" for x in item["vector"]])
            line = f"{item['mol_id']},{item['smiles']},{item['n_atoms']},{item['charge_source']},{vec_str}"
            f_out.write(line + "\n")

    f_out.close()
    print(f"\n[SUCCESS] Processed {processed_count} of {len(tasks)} molecules. Results saved to {args.out}")

if __name__ == "__main__":
    main()
