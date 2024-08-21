import gemmi, subprocess, glob, os, pickle
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from collections import defaultdict, Counter
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch, Select, PDBIO
from multiprocessing import Pool
from fix_protein import *

def mmcif_corrector(value, type = str):
    if value == '?': return None
    else:
        if type == str: return value.replace("'", "")
        elif type == int: return int(value)
        elif type == float: return float(value)
        elif type == bool: return True if value == 'yes' else False
        else: return value
  
# first step download both the pdb and cif files - pdb for final formate and cif for metadata
def download_pdb_cif(pdb_id, folder):
    """
    Download a CIF file using wget.

    Parameters:
    - pdb_id (str): The PDB ID of the file to download.
    - target_folder (str): The path to the folder where the file will be saved.
    """
    # URL for the PDB file (Replace with the base URL of your choice)
    url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
    url_pdb = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    # Construct the wget command
    command_cif = ["wget", "-q", "-P", folder, url_cif]
    command_pdb = ["wget", "-q", "-P", folder, url_pdb]

    # Execute the command
    try:
        subprocess.run(command_cif, check=True)
        subprocess.run(command_pdb, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to download PDB code: {pdb_id}. Error: {e}")
        return False

def download_ligand_sdf(pdbid, ligand_id, asym_id, auth_seq_id, folder):
    """
    Download a ligand SDF file using wget.

    Parameters:
    - pdbid (str): The PDB ID of the file to download.
    - ligand_id (str): The ligand ID associated with the PDB structure.
    - asym_id (str): The label_asym_id for the ligand in the structure.
    - auth_seq_id (int): The auth_seq_id of the ligand in the structure.
    - folder (str): The path to the folder where the file will be saved.
    """
    # URL for the SDF file
    url = f"https://models.rcsb.org/v1/{pdbid}/ligand?auth_seq_id={auth_seq_id}&label_asym_id={asym_id}&encoding=sdf&filename={pdbid}_{asym_id}_{ligand_id}.sdf"
    
    # Construct the wget command
    command = ["wget", "-q", "-P", folder, url]

    # Execute the command
    try:
        subprocess.run(command, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to download ligand {ligand_id} for PDB code: {pdbid}. Error: {e}")
        return False
  
def get_cif_header(cif_path, chains = None, lig_of_interest = None):
    """
    This function is for getting the cif header file for each individual pdb code for alignment purposes
    
    - the cif header contains rich information about the reference uniprot id, the sequence length, and the chain id
    """ 
    doc = gemmi.cif.read(cif_path)  # Using read_url to directly read from RCSB PDB
    block = doc.sole_block()
    if not block:
        raise ValueError("No Block Found")
    chain_meta_data = defaultdict(dict)
    uniprot_meta_data = {}
    entity_to_chain = defaultdict(list)

    # Prepare DataFrame for alignment data
    columns = ['pdb_id', 'chain_id', 'entity_id', 'pdbx_seq_one_letter_code', 'pdbx_seq_one_letter_code_can', 'uniprot_id', ' uniprot_name', 'ref_seq_one_letter',
               'pdb_seq_align_beg', 'pdb_seq_align_end', 'uniprot_seq_align_beg', 'uniprot_seq_align_end', 'nstd_linkage', 'nstd_monomer']
    alignment_data = []

    # Prepare DataFrame for sequence differences
    diff_columns = ['pdb_id', 'chain_id', 'pdb_res_num', 'pdb_res_name', 'uniprot_res_num', 'uniprot_res_name', 'details', 'in_align']
    sequence_differences = []
    pdb_id = block.find('_entry.', ['id'])[0][0]

    # Extracting _entity_poly information
    for item in block.find('_entity_poly.', ['entity_id', 'pdbx_strand_id', 'pdbx_seq_one_letter_code', 'pdbx_seq_one_letter_code_can', 'nstd_linkage', 'nstd_monomer']):
        chain_ids = item[1].split(",")
        for chain_id in chain_ids:
            entity_to_chain[item[0]].append(chain_id)
            if chains is None or chain_id in chains:
                entity_id = item[0]
                sequence = item[2].replace('\n', '').replace(';', '')  # Adjust sequence formatting
                canonical_sequence = item[3].replace('\n', '').replace(';', '')
                nstd_linkage = item[4]
                nstd_monomer = item[5]
                chain_meta_data[chain_id]['chain_id'] = mmcif_corrector(chain_id, str) 
                chain_meta_data[chain_id]['pdb_id'] = mmcif_corrector(pdb_id, str)
                chain_meta_data[chain_id]['entity_id'] = mmcif_corrector(entity_id, int)
                chain_meta_data[chain_id]['pdbx_seq_one_letter_code'] = mmcif_corrector(sequence, str)
                chain_meta_data[chain_id]['pdbx_seq_one_letter_code_can'] = mmcif_corrector(canonical_sequence, str)
                chain_meta_data[chain_id]['nstd_linkage'] = mmcif_corrector(nstd_linkage, bool)
                chain_meta_data[chain_id]['nstd_monomer'] = mmcif_corrector(nstd_monomer, bool)
        
    # Extracting reference sequence information
    struct_ref_table = block.find('_struct_ref.', ['db_name', 'db_code', 'entity_id', 'pdbx_db_accession', 'pdbx_seq_one_letter_code'])
    pdb_chains = struct_ref_table.find_column('entity_id')
    counter = Counter(pdb_chains)
    entity_to_skip = []
    for counter_key, counter_value in counter.items():
        if counter_value > 1:
            entity_to_skip.append(counter_key)
    entity_to_skip = set(entity_to_skip)
    chain_to_skip = []
    for key, value in entity_to_chain.items():
        if key in entity_to_skip:
            chain_to_skip += value
        
    for row in struct_ref_table:
        entity_id = row[2]
        if entity_id in entity_to_skip: continue
        if row[0] != 'UNP': continue
        db_code = row[1]
        uniprot_id = row[3]
        reference_sequence = row[4]
        uniprot_meta_data[uniprot_id] = {'db_code': db_code, 'reference_sequence': reference_sequence}

    # Extracting _struct_ref_seq category for sequence alignment
    for row in block.find('_struct_ref_seq.', ['pdbx_strand_id', 'pdbx_db_accession', 'db_align_beg', 'db_align_end', 'pdbx_auth_seq_align_beg', 'pdbx_auth_seq_align_end']):
        chain_ids = row[0].split(",")
        for chain_id in chain_ids:
            if chains is None or chain_id in chains:
                uniprot_id = row[1]
                if uniprot_id not in uniprot_meta_data:
                    reference_sequence  = None
                    uniprot_name = None
                    uniprot_id = None
                    uniprot_align_beg = None
                    uniprot_align_end = None
                    if chain_id in chain_to_skip:
                        pdb_align_beg = 1
                        pdb_align_end = len(chain_meta_data[chain_id]['pdbx_seq_one_letter_code'])
                    else:
                        pdb_align_beg = mmcif_corrector(row[4], int)
                        pdb_align_end = mmcif_corrector(row[5], int)
                else:
                    reference_sequence = mmcif_corrector(uniprot_meta_data.get(uniprot_id).get('reference_sequence').replace('\n', '').replace(';', ''), str)
                    uniprot_name = mmcif_corrector(uniprot_meta_data.get(uniprot_id).get('db_code'), str)
                    uniprot_id = mmcif_corrector(uniprot_id, str)
                    pdb_align_beg = mmcif_corrector(row[4], int)
                    pdb_align_end = mmcif_corrector(row[5], int)
                    uniprot_align_beg = mmcif_corrector(row[2], int)
                    uniprot_align_end = mmcif_corrector(row[3], int)
                chain_meta_data[chain_id]['uniprot_id'] = uniprot_id
                chain_meta_data[chain_id]['uniprot_name'] = uniprot_name
                chain_meta_data[chain_id]['ref_seq_one_letter'] = reference_sequence
                chain_meta_data[chain_id]['pdb_seq_align_beg'] = pdb_align_beg
                chain_meta_data[chain_id]['pdb_seq_align_end'] = pdb_align_end
                chain_meta_data[chain_id]['uniprot_seq_align_beg'] = uniprot_align_beg
                chain_meta_data[chain_id]['uniprot_seq_align_end'] = uniprot_align_end
    
    for key, value in chain_meta_data.items():
        alignment_data.append([value.get('pdb_id'), key, value.get('entity_id'), value.get('pdbx_seq_one_letter_code'), 
                               value.get('pdbx_seq_one_letter_code_can'), value.get('uniprot_id'), value.get('uniprot_name'), 
                               value.get('ref_seq_one_letter'), value.get('pdb_seq_align_beg'), value.get('pdb_seq_align_end'), 
                               value.get('uniprot_seq_align_beg'), value.get('uniprot_seq_align_end'), value.get('nstd_linkage'), value.get('nstd_monomer')])

    # Extracting _struct_ref_seq_dif category for sequence differences
    for row in block.find('_struct_ref_seq_dif.', ['pdbx_pdb_strand_id', 'mon_id', 'db_mon_id', 'pdbx_seq_db_seq_num', 'details', 'pdbx_auth_seq_num']):
        chain_id = row[0]
        if chains is None or chain_id in chains:
            pdb_res_name = mmcif_corrector(row[1], str)
            uniprot_res_name = mmcif_corrector(row[2], str)
            uniprot_res_num = mmcif_corrector(row[3], int)
            pdb_res_num = mmcif_corrector(row[5], int)
            details = mmcif_corrector(row[4], str)
            if pdb_res_num and chain_meta_data[chain_id]['pdb_seq_align_beg'] <= pdb_res_num and pdb_res_num <= chain_meta_data[chain_id]['pdb_seq_align_end']:
                sequence_differences.append([pdb_id, chain_id, pdb_res_num, pdb_res_name, uniprot_res_num, uniprot_res_name, details, True])
            else:
                sequence_differences.append([pdb_id, chain_id, pdb_res_num, pdb_res_name, uniprot_res_num, uniprot_res_name, details, False])
                
    # Extracting _pdbx_poly_seq_scheme category for sequence resid mapping
    residue_num_mapping = defaultdict(dict)
    for row in block.find('_pdbx_poly_seq_scheme.', ['seq_id', 'mon_id', 'pdb_seq_num', 'pdb_ins_code', 'pdb_strand_id']):
        chain_id = row[-1]
        if chains is None or chain_id in chains:
            res_name = mmcif_corrector(row[1], str)
            seq_id = mmcif_corrector(row[0], int)
            res_num = mmcif_corrector(row[2], str)
            ins_code = mmcif_corrector(row[3], str) if row[3]!='.' else ''
            resnum_total = f"{res_num}{ins_code}"
            residue_num_mapping[chain_id][seq_id] = resnum_total
    
    # Extracting _pdbx_nonpoly_scheme for asymid and strand id mapping
    lig_of_interest_mapping = {}
    for row in block.find('_pdbx_nonpoly_scheme.', ['asym_id', 'mon_id', 'pdb_strand_id', 'pdb_seq_num', 'auth_seq_num']):
        if lig_of_interest is not None and row[1] != lig_of_interest: continue
        asym_id = mmcif_corrector(row[0], str)
        pdb_strand_id = mmcif_corrector(row[2], str)
        pdb_seq_num = mmcif_corrector(row[3], int)
        auth_seq_num = mmcif_corrector(row[4], int)
        lig_of_interest_mapping[(pdb_strand_id, pdb_seq_num)] = (asym_id, auth_seq_num)
        
    # Create DataFrames
    differences_df = pd.DataFrame(sequence_differences, columns=diff_columns)
    alignment_df = pd.DataFrame(alignment_data, columns=columns)

    return alignment_df, differences_df, residue_num_mapping, lig_of_interest_mapping

def extract_pdb_information(pdb_file):
    # Define the headers we are interested in
    headers_of_interest = {"REMARK 465", "REMARK 470", "SEQRES", "MODRES", "SSBOND"}
    lines = defaultdict(list)
    interchain_ss = defaultdict(list)
    modres_list = []
    with open(pdb_file, 'r') as file:
        for line in file:
            # Check if the line starts with one of the headers of interest
            if any(line.startswith(header) for header in headers_of_interest):
                if line.startswith("REMARK"):
                    lines[line[:10]].append(line.strip())
                else:
                    lines[line[:6]].append(line.strip())
            if line.startswith("SSBOND") and len(line) > 29:
                chain1 = line[15]
                chain2 = line[29]
                if chain1 != chain2:
                    interchain_ss[chain1].append(chain2)
                    interchain_ss[chain2].append(chain1)
            if line.startswith("MODRES") and len(line) > 11:
                chain = line[16]
                res_name = line[12:15].strip()
                modres_list.append((chain, res_name))
    return lines, interchain_ss, modres_list

def extract_chain_specific_information(lines, chain_list):
    chain_info = []
    ssbond_lines = []
    for field, line in lines.items():
        if field == 'REMARK 465':
            chain_info.extend(line[:7])
            for l in line[7:]:
                if l[19] in chain_list:
                    chain_info.append(l)
        elif field == 'REMARK 470':
            chain_info.extend(line[:6])
            for l in line[6:]:
                if l[19] in chain_list:
                    chain_info.append(l)
        elif field == 'SEQRES':
            for l in line:
                if l[11] in chain_list:
                    chain_info.append(l)
        elif field == 'MODRES':
            for l in line:
                if l[16] in chain_list:
                    chain_info.append(l)
        elif field == 'SSBOND':
            for l in line:
                if l[15] in chain_list or l[29] in chain_list:
                    chain_info.append(l)
                    ssbond_lines.append(l)
    return chain_info, ssbond_lines


class CustomHetatmSelect(Select):
    def __init__(self, chain_resid_pairs):
        """
        Initialize with:
        - full_chain_id: the ID of the chain where all residues should be accepted.
        - chain_resid_pairs: a set of (chain_id, res_id) tuples for selective residue acceptance.
        """
        self.chain_resid_pairs = chain_resid_pairs

    def accept_residue(self, residue):
        """ Accept only residues in the chain_resid_pairs set. """
        chain_id = residue.get_parent().id
        res_id = residue.id[1]
        return (chain_id, res_id) in self.chain_resid_pairs

class CustomResidueSelect(Select):
    def __init__(self, full_chain_id, modified_residues):
        """
        Initialize with:
        - full_chain_id: the ID of the chain where all residues should be accepted.
        - chain_resid_pairs: a set of (chain_id, res_id) tuples for selective residue acceptance.
        """
        self.full_chain_id = full_chain_id
        modified_residues_to_include = []
        for mod in modified_residues:
            if mod.parent.id in full_chain_id:
                modified_residues_to_include.append(mod)
        self.modified_residues_to_include = modified_residues_to_include

    def accept_residue(self, residue):
        """ Accept all residues from one chain and specific residues from others. """
        if residue in self.modified_residues_to_include:
            return True
        chain_id = residue.get_parent().id
        if chain_id in self.full_chain_id and residue.id[0] == ' ':
            return True
        return False

class CustomAllSelect(Select):
    def __init__(self, full_chain_id, modified_residues, chain_resid_pairs):
        self.full_chain_id = full_chain_id
        modified_residues_to_include = []
        for mod in modified_residues:
            if mod.parent.id in full_chain_id:
                modified_residues_to_include.append(mod)
        self.modified_residues_to_include = modified_residues_to_include
        self.chain_resid_pairs = chain_resid_pairs
    
    def accept_residue(self, residue):
        """ Accept all residues from one chain and specific residues from others. """
        if residue in self.modified_residues_to_include:
            return True
        chain_id = residue.get_parent().id
        if chain_id in self.full_chain_id and residue.id[0] == ' ':
            return True
        res_id = residue.id[1]
        if (chain_id, res_id) in self.chain_resid_pairs:
            return True
        return False

def process_everything(pdb_id, ligand_id, folder, chain_dis_cutoff = 10, hetatm_distance_cutoff = 4, add_hydrogens=False, skip_nc_terminal=True):
    
    download_pdb_cif(pdb_id, folder)
    pdb_file = os.path.join(folder, f'{pdb_id}.pdb')
    # read in the key properties from the original pdb file
    key_properties, interchain_ss, modres_list = extract_pdb_information(pdb_file)
    alignment_info, diff_info, res_num_mapping, lig_of_interest_mapping = get_cif_header(pdb_file.replace('.pdb', '.cif'), lig_of_interest=ligand_id)
    pdbid = pdb_file.split('/')[-1].split('.')[0]
    structure = PDBParser().get_structure('pdb', pdb_file)[0]
    polymer_residues = []
    hetero_residues = []
    modified_residues = []
    ligand_instances = []
    
    for residue in structure.get_residues():
        if residue.id[0] != ' ':
            if ligand_id in residue.id[0]:
                ligand_instances.append(residue)
            else:
                if (residue.parent.id, residue.id[0].split("_")[-1]) in modres_list:
                    modified_residues.append(residue)
                else:
                    hetero_residues.append(residue)
        else:
            polymer_residues.append(residue)
    
    ns = NeighborSearch([atom for residue in polymer_residues for atom in residue.get_atoms()])
    all_chains_to_include = set()
    
    for ligand in ligand_instances:
        ligand_chain = ligand.parent.id
        ligand_res_id = ligand.id[1]
        #TODO: @Eric - process ligand here with the asym_id
        # asym_id, auth_seq_id = lig_of_interest_mapping[(ligand_chain, ligand_res_id)]
        
        close_chains = set([ligand_chain]) 
        
        for atom in ligand.get_atoms():
            close_chain_lig = ns.search(atom.coord, chain_dis_cutoff, level='C')
            if close_chain_lig:
                close_chains.update([cc.id for cc in close_chain_lig])
        
        for close_chain in list(close_chains)[:]:
            ss_bonded_chains = interchain_ss[close_chain]
            if len(ss_bonded_chains) > 0:
                close_chains.update(ss_bonded_chains)
        all_chains_to_include = all_chains_to_include.union(close_chains)
        
        # select the modified residues that are on the chain
        # the key thing is to match the residues with the chain id and be able to select them and write them before the TER
        # the protein files should be finished first and then the hetatms for all the other things should be written
        
        heteros_to_include = set()
        heteros_to_include.add(ligand)
        hetero_to_search = []
        for residue in hetero_residues:
            if residue.parent.id in close_chains:
                heteros_to_include.add(residue)
            else:
                hetero_to_search.append(residue)
        ligand_ns = NeighborSearch([atom for residue in polymer_residues for atom in residue.get_atoms()] + [atom for atom in ligand.get_atoms()])
        for residue in hetero_to_search:
            for atom in residue.get_atoms():
                close_chain_hetero = ligand_ns.search(atom.coord, hetatm_distance_cutoff, level='C')
                if close_chain_hetero and np.any([cc.id in close_chains for cc in close_chain_hetero]):
                    heteros_to_include.add(residue)
                    break
        # print(f'Found {len(heteros_to_include)} hetero residues close to ligand {ligand_id}')
        # extract just the chain specific information
        selection_res = CustomResidueSelect(close_chains, modified_residues)
        selection_all = CustomAllSelect(close_chains, modified_residues, [(residue.parent.id, residue.id[1]) for residue in heteros_to_include])
        selection_het_no_lig = CustomHetatmSelect( [(residue.parent.id, residue.id[1]) for residue in heteros_to_include if residue != ligand])
        extracted_pdb_file = os.path.join(folder, f'{pdbid}_{ligand_id}_{ligand_chain}_rcsb.pdb')
        protein_to_fix = os.path.join(folder, f'{pdbid}_{ligand_id}_{ligand_chain}_protein_temp.pdb')
        hetero_temp_path = os.path.join(folder, f'{pdbid}_{ligand_id}_{ligand_chain}_hetatm.pdb')
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(protein_to_fix, selection_res)
        io.save(extracted_pdb_file, selection_all)
        io.save(hetero_temp_path, selection_het_no_lig)
        chain_properties, ssbond_lines = extract_chain_specific_information(key_properties, close_chains)
        
        with open(extracted_pdb_file, "r") as f:
            pdb_lines = f.read()
        with open(extracted_pdb_file, "w") as f:
            f.write('\n'.join(chain_properties))
            f.write('\n')
            f.write(pdb_lines)
        
        # The complex file was fine from manual inspection of protein 1a4w; need to do the pdbfixer step and correct issues there.
        seqs = {}
        for _ , row in alignment_info.iterrows():
            if row['chain_id'] in close_chains:
                uniprot_id = row['uniprot_id']
                seq = row['pdbx_seq_one_letter_code']
                seqs[row['chain_id']] = (uniprot_id, seq)
        # print(f"Processing protein {extracted_pdb_file}")
        fixer, seqres_all = fix_protein(protein_to_fix, protein_to_fix.replace("_temp.pdb", ".pdb"), seqs, res_num_mapping, add_hydrogens, skip_nc_terminal)
        add_header(protein_to_fix.replace("_temp.pdb", ".pdb"), fixer, seqres_all, res_num_mapping, ssbond_lines)
        os.remove(protein_to_fix)
    
    # filter key info for the chains of interest
    alignment_info = alignment_info[alignment_info['chain_id'].isin(all_chains_to_include)]
    diff_info = diff_info[diff_info['chain_id'].isin(all_chains_to_include)]
        
    return alignment_info, diff_info, len(ligand_instances)                    

def process_wf(args):
    pdbid, ligid = args
    folder = f"/pscratch/sd/k/kysun/apo-holo-project/data_curation/CLP_PDBBIND/refined/{pdbid}"
    try:
        alignment_df, differences_df, num_lig_instances = process_everything(pdbid, ligid, folder)
    except Exception as e:
        print(f"Failed to extract cif information for {pdbid}")
        print(e)
        return None, None, pdbid
    if num_lig_instances == 0:
        return None, None, pdbid
    return alignment_df, differences_df, None

def batch_process_wf(refined_data_fp, num_processes = 64):
    """
    download and extract cif information for all the pdb files in a directory and concatenate the information into a single dataframe
    """
    df = pd.read_csv(refined_data_fp)
    total_alignment_df = []
    total_differences_df = []
    failed_pdbs = []
    args = [(row['PDB Code'], row['Ligand Name']) for _, row in df.iterrows()]
    
    with Pool(num_processes) as p:
        results = list(tqdm(p.imap(process_wf, args[:50]), total=len(args[:50])))

    for alignment_df, differences_df, failed_pdb in results:
        if alignment_df is not None:
            total_alignment_df.append(alignment_df)
        if differences_df is not None:
            total_differences_df.append(differences_df)
        if failed_pdb is not None:
            failed_pdbs.append(failed_pdb)
    
    total_alignment_df = pd.concat(total_alignment_df).reset_index(drop=True)
    total_differences_df = pd.concat(total_differences_df).reset_index(drop=True)
    return total_alignment_df, total_differences_df, failed_pdbs
  
if __name__ == "__main__":
    import warnings
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
    refined_data_fp = "/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/utils/refined_data.csv"
    total_alignment_df, total_differences_df, failed_pdbs = batch_process_wf(refined_data_fp)
    total_alignment_df.to_csv("/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/utils/alignment_data.csv", index=False)
    total_differences_df.to_csv("/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/utils/differences_data.csv", index=False)
    with open("/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/utils/failed_pdbs.pkl", "wb") as f:
        pickle.dump(failed_pdbs, f)
    