import os, sys
from io import StringIO
import json
from typing import Optional, Dict, List
from Bio.PDB import PDBIO, Residue
from rdkit import Chem
from rdkit.Chem import AllChem

from dimorphite_dl import dimorphite_dl as dl
from rcsb import get_smiles_from_rcsb

with open(os.path.join(os.path.dirname(__file__), 'manual_smiles.json')) as f:
    MANUAL_SMILES = json.load(f)


def get_reference_smi(pdb_id, ligand_id):
    if pdb_id in MANUAL_SMILES:
        smi = MANUAL_SMILES[pdb_id]
    else:
        smi = get_smiles_from_rcsb(ligand_id)
    return smi


def write_sdf(mol, sdf):
    writer = Chem.SDWriter(sdf)
    writer.write(mol)
    writer.close()


def get_formula(mol, as_dict=False):
    count = {}
    for at in mol.GetAtoms():
        if at.GetSymbol() == 'H':
            continue
        num = count.get(at.GetSymbol(), 0)
        num += 1
        count[at.GetSymbol()] = num
    if as_dict:
        return count 
    else:
        eles = list(count.keys())
        eles.sort()
        return ''.join([f'{e}{count[e]}' for e in eles])


def get_mol_info(mol, kekuleSmiles=True, isomericSmiles=False):
    formula = get_formula(mol)
    tmp_smi = Chem.MolToSmiles(mol, kekuleSmiles=kekuleSmiles, isomericSmiles=isomericSmiles)
    tmp_mol = Chem.MolFromSmiles(tmp_smi)
    if tmp_mol is None:
        smi = Chem.MolToSmiles(mol)
        key = Chem.MolToInchiKey(mol)
    else:
        smi = Chem.MolToSmiles(tmp_mol)
        key = Chem.MolToInchiKey(tmp_mol)
    mol.SetProp('LigandFixerSmiles', smi)
    mol.SetProp('LigandFixerInchiKey', key)
    mol.SetProp('LigandFixerFormula', formula)
    return smi, key, formula


def is_same_molecule(mol, ref_mol):
    '''
    Determine if the two provided molecules are the same

    Returns
    --------
    status: int
        - 0: Same molecule
        - 1: Different Protonation State / Tautomer
        - 2: Has same formula, may be can be fixed
        - 3: Different Molecule 
    '''
    if not mol.HasProp('LigandFixerSmiles'):
        get_mol_info(mol)
    if not ref_mol.HasProp('LigandFixerSmiles'):
        get_mol_info(ref_mol)

    key, ref_key = mol.GetProp('LigandFixerInchiKey'), ref_mol.GetProp('LigandFixerInchiKey')
    smi, ref_smi = mol.GetProp('LigandFixerSmiles'), ref_mol.GetProp('LigandFixerSmiles')
    formula, ref_formula = mol.GetProp('LigandFixerFormula'), ref_mol.GetProp('LigandFixerFormula')
    if key == ref_key:
        if smi == ref_smi:
            status = 0
        else:
            status = 1
    elif key[:-1] == ref_key[:-1]:
        status = 1
    else:
        if formula == ref_formula and mol.GetNumBonds() == ref_mol.GetNumBonds():
            status = 2
        else:
            status = 3
    return status


def reconstruct_mol(mol, all_bond_single=True):
    """
    Create a new molecule with all bonds SINGLE bond between heavy atoms
    """
    new_mol = Chem.RWMol()
    atom_map = {}
    num_atoms = 0
    pos = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            new_atom = Chem.Atom(atom.GetAtomicNum())
            new_atom.SetFormalCharge(atom.GetFormalCharge())
            new_mol.AddAtom(new_atom)
            atom_map[atom.GetIdx()] = num_atoms
            pos.append(mol.GetConformer().GetAtomPosition(atom.GetIdx()))
            num_atoms += 1

    for bond in mol.GetBonds():
        idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(idx1).GetSymbol() != 'H' and mol.GetAtomWithIdx(idx2).GetSymbol() != 'H':
            btype = Chem.BondType.SINGLE if all_bond_single else bond.GetBondType()
            new_mol.AddBond(atom_map[idx1], atom_map[idx2], btype)
    
    new_mol = new_mol.GetMol()
    for key, value in mol.GetPropsAsDict().items():
        if isinstance(value, int):
            new_mol.SetIntProp(key, value)
        elif isinstance(value, float):
            new_mol.SetDoubleProp(key, value)
        elif isinstance(value, bool):
            new_mol.SetBoolProp(key, value)
        elif isinstance(value, str):
            new_mol.SetProp(key, value)
        else:
            raise TypeError('Unrecognized type')
    
    conf = Chem.Conformer(num_atoms)
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, [float(x) for x in p])
    new_mol.AddConformer(conf)
    
    try:
        Chem.SanitizeMol(new_mol)
    except:
        new_mol = None
    return new_mol


def run_dl(mol):
    prev_out = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    mol_dl = dl.run_with_mol_list(
        [mol],
        min_ph=7.4,
        max_ph=7.4,
        pka_precision=0.0,
    )[0]
    sys.stdout = prev_out
    return mol_dl


def fix_valence(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'B' and len(atom.GetNeighbors()) == 4:
            atom.SetFormalCharge(-1)
        if atom.GetSymbol() == 'N' and len(atom.GetNeighbors()) == 4:
            atom.SetFormalCharge(1)
    return mol


def write_residue_to_pdb(residue: Residue, connect: Optional[Dict[int, List[int]]] = None, out_pdb: Optional[os.PathLike] = None):
    io = PDBIO()
    strio = StringIO()
    io.set_structure(residue)
    io.save(strio, preserve_atom_numbering=True)
    strio.seek(0)
    pdb_block = strio.read()
    if connect:
        record = []
        for atom in residue.get_atoms():
            atomId = atom.get_serial_number()
            if atomId in connect:
                record.append(f'CONECT {atomId:>5}' + ''.join(f'{a:>5}' for a in connect[atomId]))
        pdb_block = pdb_block.split('\n')[:-3] + record + ['END']
        pdb_block = '\n'.join(pdb_block)
    if out_pdb:
        with open(out_pdb, 'w') as f:
            f.write(pdb_block)
    return pdb_block


def read_pdbblock(block: str):
    try:
        mol = Chem.MolFromPDBBlock(block, sanitize=False)
        fix_valence(mol)
    except Exception as e:
        return None

    Chem.SanitizeMol(mol)
    mol.SetBoolProp('FromPDB', True)
    return mol



class LigandFixException(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def raise_exception(condition, msg):
    if not condition:
        raise LigandFixException(msg)
    

def fix_ligand(sdf, ref_smi, out_sdf, name=None):
    # Bad Ligand
    mol = Chem.SDMolSupplier(sdf, sanitize=False)[0]
    raise_exception(mol is not None, f'Bad ligand: {sdf}')

    # Exclude rare elements
    common = {
        'H', 'B', 'C', 'N', 'O', 'F', 
        'S', 'P', 'Cl', 'Br', 'I'
    }
    rares = [at.GetSymbol() for at in mol.GetAtoms() if at.GetSymbol() not in common]
    raise_exception(len(rares) == 0, f'Ligand has rare elements: {set(rares)}')

    # Sanitize Molecule
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        raise LigandFixException(f'Sanitize failed: {e}')
    
    # Get Reference    
    try:
        ref_mol = Chem.MolFromSmiles(ref_smi, sanitize=False)
        fix_valence(ref_mol)
        Chem.SanitizeMol(ref_mol)
        ref_smi = Chem.MolToSmiles(ref_mol, kekuleSmiles=True, isomericSmiles=False)
        ref_mol = Chem.MolFromSmiles(ref_smi)
        ref_mol = run_dl(ref_mol)
    except:
        ref_mol = None
    
    # Fix ligand
    fix_err = ""
    if ref_mol is not None:
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol)
        except:
            mol_flatten = reconstruct_mol(mol)
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol_flatten)
            except Exception as e:
                fix_err = f"Fix failed: {e}"
        
        if not fix_err:
            status = is_same_molecule(mol, ref_mol)
            fix_err = f'NOT same after fix. Error code: {status}' if status != 0 else ''

    mol_h = Chem.AddHs(mol, addCoords=True)
    for prop in mol_h.GetPropNames():
        mol_h.ClearProp(prop)
    mol_h.SetProp('_Name', os.path.basename(out_sdf) if name is None else name)
    write_sdf(mol_h, out_sdf)

    raise_exception(ref_mol is not None, "No reference found")
    raise_exception(not fix_err, fix_err)

    return mol_h