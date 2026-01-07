import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen
from typing import List, Optional, Tuple
import logging


logging.basicConfig(
    filename='skipped_molecules.log',
    filemode='a', 
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.ERROR
)

def _moments(dists: np.ndarray) -> List[float]:
    """Calculate mean, std_dev, and cube root of 3rd central moment."""
    if len(dists) == 0:
        return [0.0, 0.0, 0.0]
    
    m1 = np.mean(dists)
    m2 = np.std(dists)
    m3 = np.cbrt(np.mean((dists - m1)**3))
    
    return [m1, m2, m3]

def electroshape_vector_5d(mol: Chem.Mol) -> np.ndarray:
    """
    Generates the 18D ElectroShape descriptor (5D space: x, y, z, Charge, LogP).
    Implements Armstrong et al. (2011) parameters.
    """
    # Scaling factors from Armstrong et al. (2011)
    MU_Q = 25.0  # Angstroms per electron unit
    MU_L = 5.0   # Angstroms per LogP unit

    conf = mol.GetConformer(0)
    N = mol.GetNumAtoms()
    coords_3d = conf.GetPositions()
    
    # Get atomic properties
    charges = np.array([
        mol.GetAtomWithIdx(i).GetDoubleProp("ESPCharge")
        for i in range(N)
    ])


    try:
        contribs = AllChem.CalcCrippenContribs(mol)
        
       
        if len(contribs) != N:
            raise ValueError(f"LogP atom count mismatch. Expected {N}, got {len(contribs)}.")
            
        logp_values = np.array([x[0] for x in contribs])
        
    except Exception as e:
        
        raise ValueError(f"LogP Calculation Failed: {str(e)}")

    
    # Build 5D Coordinates
    coords_5d = np.zeros((N, 5))
    coords_5d[:, :3] = coords_3d
    coords_5d[:, 3] = charges * MU_Q
    coords_5d[:, 4] = logp_values * MU_L

    # --- Define Centroids (C1 - C6) ---
    c1 = np.mean(coords_5d, axis=0)
    
    dists_c1 = np.linalg.norm(coords_5d - c1, axis=1)
    c2 = coords_5d[np.argmax(dists_c1)]
    
    dists_c2 = np.linalg.norm(coords_5d - c2, axis=1)
    c3 = coords_5d[np.argmax(dists_c2)]
    
    # Chirality calculations for C4/C5
    vec_a_spatial = (c2 - c1)[:3]
    vec_b_spatial = (c3 - c1)[:3]
    cross_spatial = np.cross(vec_a_spatial, vec_b_spatial)
    
    len_a_s = np.linalg.norm(vec_a_spatial)
    len_cross = np.linalg.norm(cross_spatial)
    
    vec_c_spatial = cross_spatial * (len_a_s / (2 * len_cross)) if len_cross > 1e-9 else np.zeros(3)
    base_spatial = c1[:3] + vec_c_spatial
    
    q_plus = np.max(coords_5d[:, 3])
    q_minus = np.min(coords_5d[:, 3])
    l_mean = np.mean(coords_5d[:, 4])
    
    c4 = np.zeros(5); c4[:3] = base_spatial; c4[3] = q_plus; c4[4] = l_mean
    c5 = np.zeros(5); c5[:3] = base_spatial; c5[3] = q_minus; c5[4] = l_mean
    
    # C6: Lipophilic Offset Centroid (Armstrong 2011 specific)
    dist_c1_c2 = np.linalg.norm(c1 - c2)
    c6 = c1.copy()
    c6[4] += dist_c1_c2
    
    centroids = [c1, c2, c3, c4, c5, c6]
    
    # Compute Moments
    descriptor = []
    for center in centroids:
        dists = np.linalg.norm(coords_5d - center, axis=1)
        descriptor.extend(_moments(dists))
        
    return np.array(descriptor)

def prepare_mol(mol: Chem.Mol) -> Optional[Tuple[str, str]]:
    """Adds H, generates 3D conformer, extracts ID."""
    try:
        if mol.HasProp('IDNUMBER'):
            unique_id = mol.GetProp('IDNUMBER')
        elif mol.HasProp('_Name') and mol.GetProp('_Name').strip():
            unique_id = mol.GetProp('_Name')
        else:
            unique_id = Chem.MolToSmiles(mol, canonical=True)

        mol = Chem.AddHs(mol, explicitOnly=True)
        
        if mol.GetNumConformers() == 0:
            cid = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if cid == -1:
                cid = AllChem.EmbedMolecule(mol, AllChem.SR(randomSeed=42))
            if cid < 0:
                return None
        
        return Chem.MolToMolBlock(mol), unique_id
    except:
        return None

def process_record(task: Tuple[str, str]) -> Optional[dict]:
    molblock, unique_id = task
    try:
        mol = Chem.MolFromMolBlock(molblock, removeHs=False)
        if mol is None: 
            logging.error(f"MOLECULE_INVALID | ID: {unique_id} | Could not parse MolBlock")
            return None
            
        try:
            mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            if mp:
                for a in mol.GetAtoms():
                    a.SetDoubleProp("ESPCharge", float(mp.GetMMFFPartialCharge(a.GetIdx())))
            else:
                raise ValueError("MMFF partial charge calculation failed for molecule.")
        except Exception as e:

            logging.error(f"CHARGE_FAIL | ID: {unique_id} | REASON: {e}")
            return None


        try:
            vec = electroshape_vector_5d(mol)
        except ValueError as ve:

            logging.error(f"DESCRIPTOR_FAIL | ID: {unique_id} | REASON: {ve}")
            return None

        smiles = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), canonical=True)
        
        return {
            "mol_id": unique_id,
            "smiles": smiles,
            "n_atoms": mol.GetNumAtoms(),
            "charge_source": "mmff",
            "vector": vec.tolist()
        }
    except Exception as e:

        logging.error(f"UNEXPECTED_FAIL | ID: {unique_id} | REASON: {e}")
        return None
