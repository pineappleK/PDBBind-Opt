from pdbfixer import PDBFixer
from openmm.app import PDBFile
from collections import defaultdict
from textwrap import wrap

aa_mapping = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

def convert_to_seqres(sequence, chain_id='A'):
    # Convert the one-letter code to three-letter code
    result = []
    i = 0
    while i < len(sequence):
        if sequence[i] == '(':
            # Find the closing parenthesis
            closing_index = sequence.find(')', i)
            if closing_index != -1:
                # Append the content inside the parentheses as is
                result.append(sequence[i+1:closing_index])
                i = closing_index + 1
            else:
                raise ValueError("Unmatched parenthesis in the sequence.")
        else:
            # Convert the single-letter code to three-letter code
            result.append(aa_mapping.get(sequence[i], 'XAA'))  # 'XAA' for unknown residues
            i += 1
    # Group the sequence into chunks of 13 residues
    lines = wrap(' '.join(result), 51)
    # Create the SEQRES lines
    seqres_lines = []
    for i, line in enumerate(lines):
        seqres_lines.append(f"SEQRES  {i+1: >2} {chain_id} {len(sequence): >4}  {line}")
    return '\n'.join(seqres_lines)

def fix_protein(input_pdb, output_pdb, seqs = None, res_num_mapping = None, add_hydrogens=False, skip_nc_terminal=True):
    """
    Fix the protein structure using PDBFixer.

    Args:
        input_pdb: input pdb file
        output_pdb: output pdb file
    
    Returns:
        a pdb file with the missing atoms and residues added
    """
    # Load the PDB file using PDBFixer
    seqres_all = []
    if seqs:
        for k, v in seqs.items():
            seqres = convert_to_seqres(v[1], chain_id=k)
            seqres_all.append(seqres)
        with open(input_pdb, 'r') as f:
            pdb = f.read().split('\n')
        # Add SEQRES records
        pdb = [line for line in pdb if not line.startswith('SEQRES')]
        seqres_lines = '\n'.join(seqres_all)
        pdbs = '\n'.join(pdb)
        final_pdb_string = seqres_lines + '\n\n' + pdbs
        with open(input_pdb, 'w') as f:
            f.write(final_pdb_string)
    if res_num_mapping:
        with open(input_pdb, 'r') as f:
            pdb = f.read().split('\n')
        with open(input_pdb, 'w') as f:
            for line in pdb:
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                    chain = line[21]
                    res_num = int(line[22:26].strip())
                    insert_code = line[26]
                    new_res_num = [k for k, v in res_num_mapping[chain].items() if v == f"{res_num}{insert_code}".strip()][0]
                    line = line[:22] + f"{new_res_num:>4} " + line[27:]
                if line.startswith("ATOM") and int(line[22:26].strip()) == 1:
                    if line[12:16].strip() == 'H':
                        # Change 'H' to 'H1' (right-aligned in a 4-character field)
                        line = line[:12] + ' H1 ' + line[16:]
                f.write(line + '\n')
    
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    
    if skip_nc_terminal and res_num_mapping:
        item_to_pop = []
        chains = list(fixer.topology.chains())
        missing_residues_mapping = defaultdict(list)
        for info, residues in fixer.missingResidues.items():
            chain_index, res_id_start = info
            if res_id_start == 0: item_to_pop.append(info)
            chain_id = chains[chain_index].id
            padding = len(missing_residues_mapping[chain_id])
            for i, res_name in enumerate(residues):
                missing_residues_mapping[chain_id].append(res_id_start+i+1+padding)
                if res_id_start+i+1+padding == len(res_num_mapping[chain_id]):
                    item_to_pop.append(info)
        for info in item_to_pop:
            fixer.missingResidues.pop(info)
        
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(7.4)

    # Write the modified structure to a new PDB file
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'), keepIds=True)
    return fixer, seqres_all

def add_header(pdb_file, fixer, seqres_all = None, res_num_mapping = None, ssbond_lines = []):
    """
    Add necessary metadata
    """
    with open(pdb_file, 'r') as file:
        pdb_lines = file.read().split('\n')
        
    if res_num_mapping:
        pdb_lines = []
        with open(pdb_file, 'r') as f:
            pdb = f.read().split('\n')
        for line in pdb:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                chain = line[21]
                res_num = int(line[22:26].strip())
                if res_num in res_num_mapping[chain]:
                    new_res_num = res_num_mapping[chain][res_num]
                    insert_code = ' '
                    if new_res_num[-1].isalpha():
                        insert_code = new_res_num[-1]
                        new_res_num = new_res_num[:-1]
                    line = line[:22] + f"{new_res_num:>4}" + f"{insert_code:>1}" + line[27:]
            pdb_lines.append(line)

    # Create lists to hold the new REMARK and MODRES lines
    remark_465_lines = [
        "REMARK 465",
        "REMARK 465 FIXED MISSING RESIDUES",
        "REMARK 465 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
        "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
        "REMARK 465",
        "REMARK 465   M RES C SSSEQI"
    ]
    
    remark_470_lines = [
        "REMARK 470",
        "REMARK 470 FIXED MISSING ATOM",
        "REMARK 470 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
        "REMARK 470 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
        "REMARK 470",
        "REMARK 470   M RES CSSEQI  ATOMS"
    ]

    modres_lines = []
    chains = list(fixer.topology.chains())
    
    # Add missing residues to REMARK 465
    missing_residues_mapping = defaultdict(dict)
    for info, residues in fixer.missingResidues.items():
        chain_index, res_id_start = info
        chain_id = chains[chain_index].id
        padding = len(missing_residues_mapping[chain_id])
        for i, res_name in enumerate(residues):
            missing_residues_mapping[chain_id][res_id_start+i+1+padding] = res_name
            
    for chain in missing_residues_mapping:
        chain_name = chain
        residues = missing_residues_mapping[chain]
        for res_id, res_name in residues.items():
            if not res_num_mapping:
                remark_465_lines.append(f"REMARK 465     {res_name:3} {chain_name:1}  {res_id:>4}")
            else:
                res_id_mapped = res_num_mapping[chain_name][res_id]
                insert_code = ' '
                if res_id_mapped[-1].isalpha():
                    insert_code = res_id_mapped[-1]
                    res_id_mapped = res_id_mapped[:-1]
                remark_465_lines.append(f"REMARK 465     {res_name:3} {chain_name:1}  {res_id_mapped:>4}{insert_code:>1}")

    # Add missing atoms to REMARK 470
    missing_atoms = fixer.missingAtoms
    for residue in missing_atoms:
        chain_name = residue.chain.id
        res_name = residue.name
        insert_code = ' '
        if not res_num_mapping:
            res_id = residue.id
        else:
            res_id = res_num_mapping[chain_name][int(residue.id)]
            if res_id[-1].isalpha():
                insert_code = res_id[-1]
                res_id = res_id[:-1]
        atoms = missing_atoms[residue]
        atom_list = " ".join([f"{atom.name:<4}" for atom in atoms])
        remark_470_lines.append(f"REMARK 470     {res_name:3} {chain_name:1} {res_id:>4}{insert_code:>1}  {atom_list}")
    
    # Add modified residues to MODRES
    for residue, std_res_name in fixer.nonstandardResidues:
        res_name = residue.name
        chain = residue.chain.id
        if not res_num_mapping:
            res_id = residue.id
        else:
            res_id = res_num_mapping[chain][int(residue.id)]
            if res_id[-1].isalpha():
                insert_code = res_id[-1]
                res_id = res_id[:-1]
        modres_lines.append(f"MODRES {res_name:3} {chain:1} {res_id:>4}   {std_res_name:3}   MODIFIED RESIDUE")

    if len(remark_465_lines) == 6: # No missing residues
        remark_465_lines = []
    if len(remark_470_lines) == 6: # No missing atoms
        remark_470_lines = []
    
    final_pdb_strings = [pdb_lines[0]] + remark_465_lines + remark_470_lines + seqres_all + modres_lines + ssbond_lines
    
    # Reopen the PDB file for writing and insert the REMARKs and MODRES
    with open(pdb_file, 'w') as file:
        file.write('\n'.join(final_pdb_strings))
        file.write('\n')
        file.write("\n".join(pdb_lines[1:]))