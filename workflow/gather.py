import os, sys
import argparse
from tqdm import tqdm
import pandas as pd
import numpy as np
import tarfile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED
import multiprocessing as mp


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', help='Input CSV containing metadata')
    parser.add_argument('-d', '--dir', dest='dir', help='Directory of the process.py output')
    parser.add_argument('-m', '--num_workers', dest='num_workers', type=int, default=64, help="Number of workers to run in parallel")
    parser.add_argument('-o', dest='output', help='Output metadata (csv) format')

    args = parser.parse_args()

    metadata = pd.read_csv(args.input)

    def process(pdbid):
        directory = os.path.join(args.dir, pdbid)
        if not os.path.isfile(os.path.join(directory, 'done.tag')):
            return []
        
        identifiers = []
        for dirname in os.listdir(directory):
            if (not dirname.startswith('.')) and os.path.isdir(os.path.join(directory, dirname)):
                identifiers.append(dirname)

        fcif = os.path.join(directory, f'{pdbid}.cif')
        info = MMCIF2Dict(fcif)
        date = info['_pdbx_database_status.recvd_initial_deposition_date'][0]
        year = date.split('-')[0]
        
        key = '_refine.ls_d_res_high'
        if key in info:
            resolution = info[key][0]
        else:
            resolution = "NMR"

        subdf = metadata.query(f'PDBID == "{pdbid}"')
        recs = []
        for idt in identifiers:
            _, name, chain, resnum = tuple(idt.split('_'))
            record = {
                "PDBID": pdbid,
                "Resolution": resolution,
                "Year": year,
                "Ligand Name": name,
                "Ligand Chain": chain,
                "Ligand Residue Number": resnum,
                "Binding Affinity Measurement": subdf.iloc[0]['measurement'],
                "Binding Affinity Sign": subdf.iloc[0]['sign'],
                "Binding Affinity Value": subdf.iloc[0]['value'],
                "Binding Affinity Unit": subdf.iloc[0]['unit'],
                "Log Binding Affinity": subdf.iloc[0]['logvalue'],
                "Binding Affinity Source": subdf.iloc[0]['source'],
                "Binding Affinity Annotation": subdf.iloc[0]['origin']
            }
            recs.append(record)
        return recs
    
    pdb_ids = list(metadata['PDBID'].unique())
    with mp.Pool(64) as p:
        results = list(tqdm(p.imap_unordered(process, pdb_ids, chunksize=1), total=len(pdb_ids)))
    
    records = []
    for recs in results:
        if len(recs) == 0:
            continue
        records += recs

    df = pd.DataFrame(records)
    df = df.sort_values(by=['PDBID', 'Ligand Chain'])
    df.to_csv(args.output, index=None)
    
    print(f"Number of PDB entries: {df['PDBID'].unique().shape[0]}")
    print(f"Number of protein-ligand structures: {df.shape[0]}")
    #     sdfs = []
    #     for idt in identifiers:
    #         sdfs.append(os.path.join(directory, idt, f'{idt}_ligand_refined.sdf'))
        
    #     mols = [Chem.SDMolSupplier(sdf, removeHs=False)[0] for sdf in sdfs]
    #     for m in mols:
    #         Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    #         continue

    #     mol = mols[0]

    #     ligand_props = dict(
    #         mol_weight = Descriptors.MolWt(mol),
    #         logp = Crippen.MolLogP(mol),
    #         tpsa = Descriptors.TPSA(mol),
    #         num_rotatable_bonds = Lipinski.NumRotatableBonds(mol),
    #         num_heavy_atoms = mol.GetNumHeavyAtoms(),
    #         smiles = Chem.MolToSmiles(Chem.RemoveHs(mol)),
    #         num_h_donors = Lipinski.NumHDonors(mol),
    #         num_h_acceptors = Lipinski.NumHAcceptors(mol),
    #         qed = QED.qed(mol)
    #     )
        
        

    #     for rec in subdf.to_dict('records'):
    #         rec.update(ligand_props)
    #         rec.update({"Year": year, "Resolution": resolution})
    #         records.append(rec)
        
    # df = pd.DataFrame(records)
    # columns = ['PDBID', "Year", 'Resolution'] + df.columns.tolist()[1:-2]
    # df = df[columns]
    # df.to_csv(args.metadata, index=None)

    # print(f"Number of protein-ligand structures: {num_structures}")