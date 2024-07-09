import warnings
warnings.filterwarnings('ignore')
import os, glob
import numpy as np
import pandas as pd
import shutil
import json
from scipy.spatial.distance import cdist
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import MDAnalysis as mda
from MDAnalysis.topology.tables import SYMB2Z


def steric_clash_filter():
    dataset_dir = '../raw/PDBbind-v2020'
    dirs = list(glob.glob(os.path.join(dataset_dir, '*')))

    result = []
    for dirname in tqdm(dirs):
        code = os.path.basename(dirname)

        # Get Ligand Positions
        mol2 = os.path.join(dirname, f'{code}_ligand.mol2')
        ligand = Chem.MolFromMol2File(mol2, sanitize=False, cleanupSubstructures=False)
        conf = ligand.GetConformer()
        ligand_pos = []
        for i in range(ligand.GetNumAtoms()):
            if ligand.GetAtomWithIdx(i).GetSymbol() != 'H':
                ligand_pos.append(conf.GetAtomPosition(i))
        ligand_pos = np.array(ligand_pos)
        
        # Get Protein Positions
        pdb = os.path.join(dirname, f'{code}_protein.pdb')
        protein = mda.Universe(pdb, pdb)
        try:
            protein_heavy = protein.select_atoms('not element H')
        except AttributeError as e:
            protein_heavy = protein.select_atoms('(not name H*) and (not name 1H*) and (not name 2H*) and (not name 3H*)')
        except Exception as e:
            print(dirname)
            raise e
        protein_pos = protein_heavy.positions
        
        dist_mat = cdist(ligand_pos, protein_pos)
        min_dist = np.min(dist_mat)
        min_loc = np.unravel_index(np.argmin(dist_mat), dist_mat.shape)

        try:
            prot_ele = protein_heavy.atoms[min_loc[1]].element
        except mda.exceptions.NoDataError as e:
            name: str = protein_heavy.atoms[min_loc[1]].name
            if name[:2].capitalize() in SYMB2Z:
                prot_ele = name[:2]
            else:
                prot_ele = name[0]
        except Exception as e:
            print(dirname)
            raise e

        result.append((
            code, min_dist, min_loc[0], min_loc[1], 
            ligand.GetAtomWithIdx(int(min_loc[0])).GetSymbol(),
            prot_ele
        ))
    
    df = pd.DataFrame(result, columns=['PDBID', 'min_dist', 'lig_idx', 'prot_idx', 'lig_atom', 'prot_atom'])
    df.sort_values(by=['PDBID'], inplace=True)
    df.to_csv('MinDist.csv', index=None)


def non_standard_pocket_filter():

    STD_RESIDUES = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
        'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO',
        'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HID', 'HIE', 'ACE', 'NME', 'HOH'
    }

    dataset_dir = '../raw/PDBbind-v2020'
    dirs = list(glob.glob(os.path.join(dataset_dir, '*')))
    
    log = []
    for dirname in tqdm(dirs):
        code = os.path.basename(dirname)
        f_pocket = os.path.join(dirname, f'{code}_pocket.pdb')
        # pocket = parmed.load_file()
        pocket = mda.Universe(f_pocket)
        diff = []
        for res in pocket.atoms.residues:
            if res.resname not in STD_RESIDUES:
                diff.append(res.resname)
        if len(diff) > 0:
            log.append((code, ','.join(diff)))
    
    log.sort(key=lambda x: x[0])
    with open('NonStandardPocket.txt', 'w') as f:
        f.write('\n'.join(f'{entry[0]} {entry[1]}' for entry in log))



if __name__ == '__main__':
    # steric_clash_filter()
    # non_standard_pocket_filter()
    df = pd.read_csv('MinDist.csv', index_col=0)
    codes = set(df.index)
    steric_clashes = set(df.query('min_dist > 4.5 | min_dist < 2.0').index)
    
    with open('NonStandardPocket.txt') as f:
        non_std = set([line.split()[0] for line in f.read().split('\n')])

    cov_info = pd.read_csv('CovBinderInPDB_2022Q4_AllRecords.csv')
    cov_ids = set([x.lower() for x in cov_info['pdb_id']])
    covalents = codes.intersection(cov_ids)
    
    lig_fix_info = pd.read_csv('FixLog.csv', index_col=0)
    bad_ligs = set(lig_fix_info.query('CorrectStatus != 0 & CorrectStatus != 1').index)
    # bad_ligs = set()

    small_qeds = set(lig_fix_info.dropna(subset=['QED']).query('QED < 0.2').index)
    
    print('Pocket with non-standard residues:', len(non_std))
    print('Steric Clashes:', len(steric_clashes))
    print('Bad Ligands:', len(bad_ligs))
    print('Small QEDs:', len(small_qeds))
    print('After filter:', len(codes - non_std - steric_clashes - cov_ids - bad_ligs - small_qeds))
    


    
    



        



