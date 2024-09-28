import os
import re
from copy import deepcopy
import json
import itertools
from typing import List, Optional, Union, Tuple
from collections import defaultdict
import numpy as np
import pandas as pd
from tqdm import tqdm

from rcsb import download_pdb_cif
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

import openmm.app as app
from pdbfixer import PDBFixer
import networkx as nx


building_blocks = pd.read_csv(os.path.join(os.path.dirname(__file__), 'building_blocks.csv'), index_col=0)


def mol_from_seq(seq):
    rxn_link = AllChem.ReactionFromSmarts('[N:1]-[C:2]-[C:3](=[O:4])[O:5].[N:6]-[C:7]-[C:8]=[O:9]>>[N:1]-[C:2]-[C:3](=[O:4])-[N:6]-[C:7]-[C:8]=[O:9].[O:5]')
    rxn_ncap = AllChem.ReactionFromSmarts('[C:1](=[O:2])[O:3].[N:4]-[C:5]-[C:6]=[O:7]>>[C:1](=[O:2])-[N:4]-[C:5]-[C:6]=[O:7].[O:3]')
    rxn_ccap = AllChem.ReactionFromSmarts('[N:1]-[C:2]-[C:3](=[O:4])[O:5].[#7H1,#7H2,#7H3:6]>>[N:1]-[C:2]-[C:3](=[O:4])-[N:6].[O:5]')
    # Generic amide rxn
    rxn_general = AllChem.ReactionFromSmarts('[C:1](=[O:2])[O:3].[#7H1,#7H2,#7H3:4]>>[C:1](=[O:2])[N:4].[O:3]')

    for i in range(len(seq)):
        assert seq[i] in building_blocks.index, f'Unsupported residue: {seq[i]}'
        assert building_blocks.loc[seq[i], 'Category'] != 'other', f'Not amino-acid: {seq[i]}'
        if i == 0:
            product = Chem.MolFromSmiles(building_blocks.loc[seq[i], 'Smiles'])
            continue
        
        if i == 1 and building_blocks.loc[seq[0], 'Category'] == 'ncap':
            if building_blocks.loc[seq[i], 'Category'] == 'ccap':
                # ACE-NME cases
                rxn = rxn_general
            else:
                rxn = rxn_ncap
        elif i == len(seq) - 1 and building_blocks.loc[seq[i], 'Category'] == 'ccap':
            rxn = rxn_ccap
        else:
            rxn = rxn_link

        rxn_res = rxn.RunReactants((product, Chem.MolFromSmiles(building_blocks.loc[seq[i], 'Smiles'])))
        product = rxn_res[0][0]
        Chem.SanitizeMol(product)
    return product


def form_disulfide_bond(m):
    mnoh = Chem.RemoveHs(m)
    thiols = [match[0] for match in mnoh.GetSubstructMatches(Chem.MolFromSmarts('[SH:1]-[C:2]'))]
    results = []
    for i, j in itertools.combinations(thiols, 2):
        new_mol = Chem.RWMol(mnoh)
        new_mol.AddBond(i, j, Chem.BondType.SINGLE)
        new_mol.GetAtomWithIdx(i).SetNumExplicitHs(0)
        new_mol.GetAtomWithIdx(j).SetNumExplicitHs(0)
        mol = new_mol.GetMol()
        Chem.SanitizeMol(mol)
        results.append(mol)
    return results



def process_everything(
    pdb_id: str, 
    ligand_id: Optional[str] = None, 
    ligand_info: Optional[List[Tuple[str, List[int]]]] = None, 
    dataset_dir: Optional[os.PathLike] = '../raw_data',
    binding_cutoff: float = 10.0, 
    hetatm_cutoff: float = 4.0,
    find_connected_ligand_residues: bool = True
):
    """
    Run process workflow
    """

    # Create folder
    if not os.path.isdir(dataset_dir):
        os.path.mkdir(dataset_dir)
    folder = os.path.join(dataset_dir, pdb_id)
    if not os.path.isdir(folder):
        os.mkdir(folder)

    download_pdb_cif(pdb_id, folder)


def create_residue_connection_graph(topology: app.Topology, positions=None):
    graph = nx.Graph()
    graph.add_nodes_from([res for res in topology.residues()])
    for bond in topology.bonds():
        res1, res2 = bond.atom1.residue, bond.atom2.residue
        # don't include metal bonds
        if len(res1) == 1 or len(res2) == 1:
            continue
        if not (res1 is res2):
            graph.add_edge(res1, res2)
    
    # In PDB topology, NH2 is not connected by default
    prev_residue = None
    for residue in topology.residues():
        if (residue.name == 'NH2') and (prev_residue is not None):
            n_index = [at.index for at in residue.atoms() if at.element.symbol != 'H'][0]
            vec = positions[n_index]
            nh2_pos = np.array([vec.x, vec.y, vec.z])
            prev_residue_pos = []
            for at in prev_residue.atoms():
                prev_residue_pos.append([positions[at.index].x, positions[at.index].y, positions[at.index].z])
            prev_residue_pos = np.array(prev_residue_pos)
            min_dist = np.min(np.linalg.norm(prev_residue_pos - nh2_pos, axis=1))
            if min_dist < 0.2:
                graph.add_edge(prev_residue, residue)
        prev_residue = residue
    return graph



def find_short_polymers(pdb_id, dataset_dir='../raw_data_pdbbind_poly'):
    datas = []
    try:
        pdb_file = os.path.join(dataset_dir, pdb_id, pdb_id + '.pdb')
        short_chains_from_struct = defaultdict(list)
        short_chains_from_seq = {}

        fixer = PDBFixer(pdb_file)
        for seq in fixer.sequences:
            if len(seq.residues) > 20 or len(seq.residues) < 2:
                continue
            short_chains_from_seq[seq.chainId] = seq.residues
        
        top = fixer.topology
        graph = create_residue_connection_graph(top, fixer.positions)
        
        for subgraph in nx.connected_components(graph):
            if len(subgraph) > 20 or len(subgraph) < 2:
                continue
            residues = list(subgraph)
            residues.sort(key=lambda x: x.index)
            chain_id = residues[0].chain.id
            short_chains_from_struct[chain_id] += residues
        
        for chain_id in set(chain.id for chain in top.chains()):
            pdb_seq = short_chains_from_seq.get(chain_id, [])
            
            residues = short_chains_from_struct[chain_id]
            for chain in top.chains():
                if chain.id == chain_id:
                    for residue in chain.residues():
                        if (residue.name in pdb_seq) and (residue not in residues):
                            residues.append(residue)
            residues.sort(key=lambda x: x.index)

            resnames = ' '.join(r.name for r in residues)
            resids = ' '.join(f'{r.id}{r.insertionCode}'.strip() for r in residues)
            if pdb_seq or resnames or resids:
                datas.append({
                    "PDBID": pdb_id, "chain": chain_id,
                    "pdb_seq": ' '.join(pdb_seq),
                    "resnames": resnames,
                    "resids": resids,
                    "note": None
                })
    except Exception as e:
        datas = [{'PDBID': pdb_id, 'chain': None, 'pdb_seq': None, 'resnames': None, 'resids': None, 'note': str(e)}]
    return datas


if __name__ == '__main__':
    import pandas as pd
    import multiprocessing as mp

    dataset_dir = '../raw_data_pdbbind_poly'
    if not os.path.isdir(dataset_dir):
        os.mkdir(dataset_dir)
    
    df = pd.read_csv('../biolip/PDBBind_poly.csv')
    # args = [[idx] for idx in df["PDBID"].unique()]
    args = df["PDBID"].unique().tolist()
    
    def wrap_process_wf(args):
        try:
            process_everything(args[0], dataset_dir=dataset_dir)
            return None
        except:
            return args[0]

    with mp.Pool(10) as p:
        results = list(tqdm(p.imap_unordered(find_short_polymers, args, chunksize=1), total=len(args)))
    
    infos = [data for datas in results for data in datas]
    infos = pd.DataFrame(infos)
    infos.sort_values(by=['PDBID', 'chain'])
    infos.to_csv('short_chain.csv', index=None)
