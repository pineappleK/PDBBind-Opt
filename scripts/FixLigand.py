import warnings
warnings.filterwarnings('ignore')
import os, glob, sys
import shutil
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, QED
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from dimorphite_dl import dimorphite_dl as dl
from openbabel import openbabel


# SUPERSEDES = {
#     '4pox': '5ovz', '5ab1': '6yhw',
#     '4jdf': '7oyz', '1qon': '6xyu'
# }

# with open('filter_ligands.txt') as f:
#     EXCLUDE = set([line.split()[0] for line in f.read().split('\n') if line])

# EXCLUDE.update([
#     "ZN", "MG", "NA", "CA", "CL", "I", "CO", 'NI', 'K', 'MN'
#     "DMS", "PO4", "SO4", "ACT", "MBO", "PAE", "EDO",
#     "TFA", "GOL", "Y1", 'OXL', 'DIO', 'PEG', 'MRD', 'CO3', 'CO2',
#     'PDO', 'AZI', 'NAG', 'PMP'
# ])


AMINO_ACIDS = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'HID': 'H', 'HIE': 'H'
}


def parse_mol2_sequence(file_path):
    seq = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            if line.startswith("@<TRIPOS>"):
                atom_section = False
            if atom_section:
                parts = line.split()
                if len(parts) > 7:  # Ensure there are enough parts in the line
                    atom_name = parts[1]
                    res_name = parts[7].strip()
                    if res_name not in AMINO_ACIDS:
                        seq = []
                        break
                    elif atom_name == 'N':
                        seq.append(AMINO_ACIDS[res_name])
    return ''.join(seq)


def parse_mol2_residues(file_path):
    residue_names = set()
    with open(file_path, 'r') as file:
        lines = file.readlines()
        atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            if line.startswith("@<TRIPOS>"):
                atom_section = False
            if atom_section:
                parts = line.split()
                if len(parts) > 7:  # Ensure there are enough parts in the line
                    res_name = parts[7]
                    residue_names.add(res_name.strip())
    return residue_names


def get_ligand_names(dataset_dir, index_file):
    '''
    Get name of the ligand. Ligand with name None is a 'polymer', e.g. polypeptide
    '''
    ligand_names = {}
    with open(index_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('//')
            code = parts[0].split()[0]
            name = parts[1].split()[1][1:-1]
            ligand_names[code] = name

    polymers = set()
    for code, name in ligand_names.items():
        if '-' in name or 'mer' in name:
            polymers.add(code)
            
    for code in tqdm(ligand_names, desc='Determine polymers from MOL2 file'):
        fpath = os.path.join(dataset_dir, f'{code}/{code}_ligand.mol2')
        residues = parse_mol2_residues(fpath)
        
        if len(residues) > 1:
            polymers.add(code)
        
    polymers = list(polymers)
    polymers.sort()

    for code in polymers:
        ligand_names[code] = None
    
    return ligand_names


def get_formula(mol):
    count = {}
    for at in mol.GetAtoms():
        if at.GetSymbol() == 'H':
            continue
        num = count.get(at.GetSymbol(), 0)
        num += 1
        count[at.GetSymbol()] = num
    return count 


SANITIZE_FLAG = Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY


def get_terminal_oxygens(atom):
    ter_os = []
    for nei in atom.GetNeighbors():
        if nei.GetSymbol() == 'O' and len(nei.GetNeighbors()) == 1:
            ter_os.append(nei)
    return ter_os


def fix_phosphate(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "P" and len(atom.GetNeighbors()) == 4:
            ter_os = get_terminal_oxygens(atom)
            atom.SetIsAromatic(False)
            atom.SetFormalCharge(0)
            atom.SetNumRadicalElectrons(0)
            ter_os[0].SetIsAromatic(False)
            bo1 = mol.GetBondBetweenAtoms(ter_os[0].GetIdx(), atom.GetIdx())
            bo1.SetIsAromatic(False)
            bo1.SetBondType(Chem.BondType.DOUBLE)
            for o in ter_os[1:]:
                o.SetIsAromatic(False)
                bo = mol.GetBondBetweenAtoms(o.GetIdx(), atom.GetIdx())
                bo.SetIsAromatic(False)
                bo.SetBondType(Chem.BondType.SINGLE)
                o.SetFormalCharge(-1)
        atom.UpdatePropertyCache()
    return mol


def fix_carboxylate(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and len(atom.GetNeighbors()) == 3:
            ter_os = get_terminal_oxygens(atom)
            if len(ter_os) == 2:
                atom.SetIsAromatic(False)
                ter_os[0].SetIsAromatic(False)
                ter_os[1].SetIsAromatic(False)
                bo1 = mol.GetBondBetweenAtoms(ter_os[0].GetIdx(), atom.GetIdx())
                bo1.SetIsAromatic(False)
                bo1.SetBondType(Chem.BondType.DOUBLE)
                
                ter_os[1].SetFormalCharge(-1)
                bo2 = mol.GetBondBetweenAtoms(ter_os[1].GetIdx(), atom.GetIdx())
                bo2.SetIsAromatic(False)
                bo2.SetBondType(Chem.BondType.SINGLE)
        atom.UpdatePropertyCache()
    return mol


def fix_guanidine(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 3:
            ns = [nei for nei in atom.GetNeighbors() if nei.GetSymbol() == 'N']
            
            num_ar_bonds = 0
            num_non_ring_ar_bonds = 0
            for i, n in enumerate(ns):
                bo = mol.GetBondBetweenAtoms(n.GetIdx(), atom.GetIdx())
                if bo.GetIsAromatic():
                    num_ar_bonds += 1
                    if not bo.IsInRing():
                        num_non_ring_ar_bonds += 1
            
            if num_ar_bonds >= 2 and num_non_ring_ar_bonds >= 1:
                atom.SetIsAromatic(False)
                flag = True
                for n in ns:
                    n.SetIsAromatic(False)
                    bo = mol.GetBondBetweenAtoms(n.GetIdx(), atom.GetIdx())
                    bo.SetIsAromatic(False)
                    if sum([a.GetSymbol() == 'H' for a in n.GetNeighbors()]) == 2 and flag:
                        n.SetFormalCharge(1)
                        bo.SetBondType(Chem.BondType.DOUBLE)
                        flag = False
                    else:
                        bo.SetBondType(Chem.BondType.SINGLE)
                        n.SetFormalCharge(0)
                assert (not flag), "Guandine is not fixed"
        atom.SetNumRadicalElectrons(0)
        atom.UpdatePropertyCache()
    return mol


def fix_sulfonate(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S' and len(atom.GetNeighbors()) == 4:
            ter_os = get_terminal_oxygens(atom)
            atom.SetIsAromatic(False)
            for o in ter_os[:2]:
                o.SetIsAromatic(False)
                bo = mol.GetBondBetweenAtoms(o.GetIdx(), atom.GetIdx())
                bo.SetBondType(Chem.BondType.DOUBLE)
                bo.SetIsAromatic(False)
                o.SetFormalCharge(0)
            for o in ter_os[2:]:
                o.SetIsAromatic(False)
                bo = mol.GetBondBetweenAtoms(o.GetIdx(), atom.GetIdx())
                bo.SetBondType(Chem.BondType.SINGLE)
                bo.SetIsAromatic(False)
                o.SetFormalCharge(-1)
        atom.UpdatePropertyCache()
    return mol


def fix_pyridine_oxide(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetIsAromatic() and atom.IsInRing():
            ter_os = get_terminal_oxygens(atom)
            if len(ter_os) == 1:
                o = ter_os[0]
                mol.GetBondBetweenAtoms(o.GetIdx(), atom.GetIdx()).SetBondType(Chem.BondType.SINGLE)
                o.SetFormalCharge(-1)
                atom.SetFormalCharge(1)
    return mol


def hard_fix_carbon(mol):
    for atom in mol.GetAtoms():
        numVal = int(sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()]))
        if numVal == 5:
            for b in atom.GetBonds():
                if b.GetBondType() == Chem.BondType.DOUBLE:
                    b.SetBondType(Chem.BondType.SINGLE)
                    break
    return mol


def fix_aromatic(mol):
    for atom in mol.GetAtoms():
        if not atom.IsInRing():
            atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        if not bond.IsInRing():
            bond.SetIsAromatic(False)
    return mol


def fix_boron(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'B' and len(atom.GetNeighbors()) == 4:
            atom.SetFormalCharge(-1)
    return mol


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
    conf.SetPositions(np.array(pos))
    new_mol.AddConformer(conf)
    
    try:
        Chem.SanitizeMol(new_mol)
    except:
        new_mol = None
    return new_mol


def AssignBondOrdersFromTemplate(refmol, mol):
    refmol2 = Chem.Mol(refmol)
    mol2 = Chem.Mol(mol)
    # do the molecules match already?
    matching = mol2.GetSubstructMatch(refmol2)
    if not matching:  # no, they don't match
        # check if bonds of mol are SINGLE
        for b in mol2.GetBonds():
            if b.GetBondType() != Chem.BondType.SINGLE:
                b.SetBondType(Chem.BondType.SINGLE)
                b.SetIsAromatic(False)
        # set the bonds of mol to SINGLE
        for b in refmol2.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        # set atom charges to zero;
        for a in refmol2.GetAtoms():
            a.SetFormalCharge(0)
        for a in mol2.GetAtoms():
            a.SetFormalCharge(0)

    matching = mol2.GetSubstructMatches(refmol2, uniquify=False)
    # do the molecules match now?
    if matching:
        if len(matching) > 1:
            warnings.warn("More than one matching pattern found - picking one")
        matching = matching[0]
        # apply matching: set bond properties
        for b in refmol.GetBonds():
            atom1 = matching[b.GetBeginAtomIdx()]
            atom2 = matching[b.GetEndAtomIdx()]
            b2 = mol2.GetBondBetweenAtoms(atom1, atom2)
            b2.SetBondType(b.GetBondType())
            b2.SetIsAromatic(b.GetIsAromatic())
        # apply matching: set atom properties
        for a in refmol.GetAtoms():
            a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
            a2.SetHybridization(a.GetHybridization())
            a2.SetIsAromatic(a.GetIsAromatic())
            a2.SetNumExplicitHs(a.GetNumExplicitHs())
            a2.SetFormalCharge(a.GetFormalCharge())
        # Don't sanitize. Will cause aromaticity errors
        # Chem.SanitizeMol(mol2)
        if hasattr(mol2, '__sssAtoms'):
            mol2.__sssAtoms = None  # we don't want all bonds highlighted
    else:
        raise ValueError("No matching found")
    
    # set num radicals to 0
    for atom in mol2.GetAtoms():
        atom.SetNumRadicalElectrons(0)
    mol2.UpdatePropertyCache()
    return mol2


def read_by_obabel(path, ob_fmt='xyz', secondary_conversion=True):
    # surpress warnings
    openbabel.obErrorLog.SetOutputLevel(0)

    suffix = path.split('.')[-1]
    in_fmt = 'mol' if suffix == 'sdf' else suffix

    ob_mol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(in_fmt, ob_fmt)
    obConversion.ReadFile(ob_mol, path)
    ob_mol.DeleteHydrogens()
    ob_str = obConversion.WriteString(ob_mol)

    # Secondary Conversion
    if secondary_conversion:
        ob_mol = openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(ob_fmt, 'mol')
        obConversion.ReadString(ob_mol, ob_str)
        ob_str = obConversion.WriteString(ob_mol)
        ob_fmt = 'mol'

    if ob_fmt == 'mol':
        rd_mol = Chem.MolFromMolBlock(ob_str)
    elif ob_fmt == 'mol2':
        rd_mol = Chem.MolFromMol2Block(ob_str)
    elif ob_fmt == 'xyz':
        rd_mol = Chem.MolFromXYZBlock(ob_str)
    elif ob_fmt == 'pdb':
        rd_mol = Chem.MolFromPDBBlock(ob_str)
    else:
        raise NotImplementedError(f'Unsupported format {ob_fmt}')

    return rd_mol


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


def mol_to_smiles_and_key(mol, kekuleSmiles=True):
    tmp_smi = Chem.MolToSmiles(mol, kekuleSmiles=kekuleSmiles, isomericSmiles=False)
    tmp_mol = Chem.MolFromSmiles(tmp_smi)
    if tmp_mol is None:
        smi = Chem.MolToSmiles(mol)
        key = Chem.MolToInchiKey(mol)
    else:
        smi = Chem.MolToSmiles(tmp_mol)
        key = Chem.MolToInchiKey(tmp_mol)
    return smi, key


COMMON_ELEMNTS = {
    'H', 'B', 'C', 'N', 'O', 'F', 'S', 'P', 'Cl', 'Br', 'I'
}


if __name__ == '__main__':
    # Covalent binders from https://yzhang.hpc.nyu.edu/CovBinderInPDB/
    covalents = set(pd.read_csv('CovBinderInPDB_2022Q4_AllRecords.csv')['pdb_id'])
    covalents = set([c.lower() for c in covalents])

    with open('PDBbind-v2020-lig-info.json') as f:
        infos = json.load(f)
    
    ligand_names = list(infos.keys())
    ligand_names.sort()
    # ligand_names = ['2ohp']

    do_log = True
    use_key = True
    
    results = []
    for code in tqdm(ligand_names, desc='Fixing ligands'):
        result = {
            'PDBID': code,
            'RareElement': False,
            'ReadStatus': -1,
            'CorrectStatus': -1,
            'QED': None,
            'NumHeavyAtoms': None,
            'Covalent': code in covalents,
            'Smi': '',
            'OriSmi': '',
            'RefSmi': '',
            'FixSmi': '',
            'FixErrMsg': '',
            'Peptide': False,
            'PeptideRefSmi': '',
            'PeptideFixSmi': '',
            'PeptideFixErrMsg': '',
        }
        # Read molecule
        f_mol2 = f'../raw/PDBbind-v2020/{code}/{code}_ligand.mol2'
        mol = Chem.MolFromMol2File(f_mol2, sanitize=False, removeHs=False, cleanupSubstructures=False)
        # Rare elements
        if not set(at.GetSymbol() for at in mol.GetAtoms()).issubset(COMMON_ELEMNTS):
            result['RareElement'] = True
            results.append(result)
            continue

        # Read Status:
        # 0 - Success with simple rdkit.Chem.Sanitize
        # 1 - Success with simple fix some funcional groups
        # 2 - Success with using OpenBabel convert to PDB then back to molblock
        # 3 - Success with reset all bond order to one. This completely destroy the molecule and need template matching

        # Try to sanitize
        try:
            Chem.SanitizeMol(mol)
            mol.SetIntProp('ReadStatus', 0)
        except Exception as e:
            # Santize failed, try to fix
            try:
                mol = fix_phosphate(mol)
                mol = fix_carboxylate(mol)
                mol = fix_guanidine(mol)
                mol = fix_sulfonate(mol)
                mol = fix_pyridine_oxide(mol)
                mol = fix_aromatic(mol)
                mol = fix_boron(mol)
                Chem.SanitizeMol(mol)
                mol.SetIntProp('ReadStatus', 1)
            except Exception as e:
                # Fix failed, try to parse with OpenBabel
                mol_by_obabel = read_by_obabel(f_mol2, ob_fmt='pdb', secondary_conversion=True)
                if (mol_by_obabel is not None) and (mol_by_obabel.GetNumBonds() > 0):
                    mol = mol_by_obabel
                    mol.SetIntProp('ReadStatus', 2)
                else:
                    # OpenBabel failed again, reset all bond order to 1
                    mol = reconstruct_mol(mol, True)
                    if mol is not None:
                        mol.SetIntProp('ReadStatus', 3)
                    else:
                        # All methods failed, continue
                        results.append(result)
                        continue
        
        mol = Chem.AddHs(mol, addCoords=True)
        read_status = mol.GetIntProp('ReadStatus')
        mol_block = Chem.MolToMolBlock(mol)
        mol = Chem.MolFromMolBlock(mol_block)
        mol.SetIntProp('ReadStatus', read_status)
        
        smi, key = mol_to_smiles_and_key(mol)
        result['OriSmi'] = smi
        formula = get_formula(mol)

        # Get reference
        ref_mols, ref_smis, ref_keys = [], [], []
        for info in infos[code]:
            ref_mol = None
            if info['smiles']:
                # some SMILES has issues
                ref_mol = Chem.MolFromSmiles(info['smiles'].replace('\n', "").replace(" ", ''))
            if ref_mol is None and info['inchi']:
                ref_mol = Chem.MolFromInchi(info['inchi'])
            
            # No reference or reference is an Ion
            if (ref_mol is None) or (ref_mol.GetNumHeavyAtoms() == 1):
                continue

            # DL predict Protonation State 
            ref_mol = run_dl(ref_mol)
            # Regularize Aromaticity
            ref_mol = Chem.MolFromSmiles(
                Chem.MolToSmiles(ref_mol, kekuleSmiles=True, isomericSmiles=False),
                sanitize=False
            )
            ref_mol = Chem.RemoveHs(ref_mol)
            ref_smi, ref_key = mol_to_smiles_and_key(ref_mol)
            ref_smis.append(ref_smi)
            ref_mols.append(ref_mol)
            ref_keys.append(ref_key)
        
        # Check the correctness of the molecule
        # 0 - Correct without template matching
        # 1 - Correct with template matching
        # 2 - Fail to found reference in PDB
        # 3 - Found references in PDB, but none of them matches the molecule
        # 4 - Found potential correct reference according to heavy atom formula matching, e.g. C3O3N4 == C3O3N4
        #     but template matching failes
        # 5 - Is a peptide, but fix according to sequence is failed, maybe is a cyclic peptide
        
        # No ref found
        if len(ref_mols) == 0:
            mol.SetIntProp('CorrectStatus', 2)
        else:
            # Ref found, need to match refs
            for ref_mol, ref_smi, ref_key in zip(ref_mols, ref_smis, ref_keys):
                condition = (key == ref_key) if use_key else smi == ref_smi
                if condition:
                    mol.SetIntProp('CorrectStatus', 0)
                    break
                # Found match
                if get_formula(ref_mol) == formula:
                    try:
                        fix_mol = AssignBondOrdersFromTemplate(ref_mol, mol)
                        fix_smi, fix_key = mol_to_smiles_and_key(fix_mol)
                        err_msg = 'AssignBondOrderFromTemplate Success'
                    except Exception as e:
                        err_msg = str(e)
                        fix_smi, fix_key = '', ''
                        fix_mol = None
                        # raise e
                    
                    condition = (fix_key == ref_key) if use_key else fix_smi == ref_smi
                    if condition:
                        fix_mol.SetIntProp('CorrectStatus', 1)
                        mol = fix_mol
                    else:
                        mol.SetIntProp('CorrectStatus', 4)
                    
                    result.update({'FixErrMsg': err_msg, 'FixSmi': fix_smi, 'RefSmi': ref_smi})
                    break
            else:
                mol.SetIntProp('CorrectStatus', 3)
    
        # Try if we can verify the correctness if the molecule is a non-cyclic peptide 
        if mol.GetIntProp('CorrectStatus') > 1:
            seq = parse_mol2_sequence(f_mol2)
            if seq:
                ref_mol = Chem.RemoveHs(run_dl(Chem.MolFromSequence(seq)))
                ref_smi, ref_key = mol_to_smiles_and_key(ref_mol)
                try:
                    fix_mol = AssignBondOrdersFromTemplate(ref_mol, mol)
                    fix_smi, fix_key = mol_to_smiles_and_key(fix_mol)
                    err_msg = 'AssignBondOrderFromTemplate Success'
                except Exception as e:
                    err_msg = str(e)
                    fix_smi, fix_key = '', ''
                    fix_mol = None
                
                condition = (fix_key == ref_key) if use_key else fix_smi == ref_smi
                if condition:
                    fix_mol.SetIntProp('CorrectStatus', 1)
                    mol = fix_mol
                else:
                    mol.SetIntProp('CorrectStatus', 5)
                
                result.update({
                    'Peptide': True,
                    'PeptideRefSmi': ref_smi,
                    'PeptideFixSmi': fix_smi,
                    'PeptideFixErrMsg': err_msg
                })
        
        result.update({
            'ReadStatus': mol.GetIntProp('ReadStatus'),
            'CorrectStatus': mol.GetIntProp('CorrectStatus'),
            'Smi': mol_to_smiles_and_key(mol)[0],
            'QED': QED.qed(mol),
            'NumHeavyAtoms': mol.GetNumHeavyAtoms(),
        })
        results.append(result)
        mol = Chem.AddHs(mol, addCoords=True)
        with Chem.SDWriter(f_mol2.replace('.mol2', '_processed.sdf')) as w:
            w.write(mol)

    df = pd.DataFrame(results)
    print('Success:', df.query('CorrectStatus == 0 | CorrectStatus == 1').shape[0])
    print('Rare Element:', df.query('RareElement').shape[0])
    print('Fail Read:', df.query('~RareElement & ReadStatus == -1').shape[0])
    print('Peptide:', df.query('Peptide').shape[0])
    print('Peptide Fail:', df.query('Peptide & CorrectStatus != 1').shape[0])
    if do_log:
        df.to_csv('FixLog.csv', index=None)