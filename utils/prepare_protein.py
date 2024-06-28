"""
Prepare the protein structure for the energy-free building process.

Author: Oliver Sun
Date: 01/15/2024
"""
from pdbfixer import PDBFixer
from Bio.PDB import PDBParser, Superimposer, PDBIO
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
sys.path.append('/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/')
from utils.alphafold import get_af2_protein

def process_protein(pdb_seq = None, pdb_file = None, af2_path = None, clean_path=None, keep_water=False, add_hydrogens=True, get_af2_structure = False):
    """
    Process the protein structure by fixing the missing atoms and residues.
    If pdb file is provided, the protein will be fixed and the sequence will be extracted.
    If pdb sequence is provided, the protein will be fixed and the sequence will be used to get the AF2 protein.
    When both pdb file and pdb sequence are provided, the pdb file will be used to extract sequence.
    
    return af2_path, clean_path, pdb_seq
    """
    if not pdb_seq and not pdb_file:
        raise ValueError("Either pdb_seq or pdb_files must be provided.")

    if pdb_file:
        fix_protein(pdb_file, clean_path, keep_water, add_hydrogens)
        pdb_seq = extract_sequence(clean_path)
        if not get_af2_structure:
            return None, clean_path, pdb_seq
    
    if pdb_seq:
        get_af2_protein(pdb_seq, af2_path, clean_path)
        fix_protein(af2_path, af2_path, keep_water, add_hydrogens)
        if not pdb_file:
            return af2_path, None, pdb_seq
        else:
            align_structure(af2_path, clean_path)
            return af2_path, clean_path, pdb_seq
    
def align_structure(af2_path, clean_path):
    """
    Align the AlphaFold2 structure to the fixed structure.
    
    Args:
        af2_path (str): Path to the AlphaFold2 structure file.
        clean_path (str): Path to the fixed structure file.
    
    Returns:
        str: Path to the aligned AlphaFold2 structure file.
    """
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure('AF2_structure', af2_path)
    structure2 = parser.get_structure('PDB_structure', clean_path)
    
    # Extract the C-alpha atoms
    ca_atoms1 = []
    ca_atoms2 = []
    
    for model in structure1:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    ca = residue['CA']
                    ca_atoms1.append(ca)
    
    for model in structure2:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    ca = residue['CA']
                    ca_atoms2.append(ca)
    
    # Check if we have the same number of C-alpha atoms in both structures
    if len(ca_atoms1) != len(ca_atoms2):
        raise ValueError("The number of C-alpha atoms does not match between the two structures.")
    
    # Align the two structures
    sup = Superimposer()
    sup.set_atoms(ca_atoms2, ca_atoms1)  # Reference (clean) structure first
    sup.apply(structure1)  # Apply transformation to the AlphaFold2 structure
    
    # Save the aligned structure
    io = PDBIO()
    io.set_structure(structure1)
    io.save(af2_path)
    
    return

# input should be fixed with hydrogens already
def extract_sequence(input_pdb):
    """
    Extract the sequence of the protein structure.

    Args:
        input_pdb: input pdb file
    
    Returns:
        str: the sequence of the protein structure
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', input_pdb)
    
    sequences = {}
    
    # Three-letter to one-letter amino acid code mapping
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = []
            
            for residue in chain:
                if residue.id[0] == ' ':  # Exclude heteroatoms
                    resname = residue.resname
                    if resname in three_to_one:
                        sequence.append(three_to_one[resname])
            
            sequences[chain_id] = "".join(sequence)
    
    return sequences

def fix_protein(input_pdb, output_pdb, keep_water=False, add_hydrogens=True):
    """
    Fix the protein structure using PDBFixer.

    Args:
        input_pdb: input pdb file
        output_pdb: output pdb file
    
    Returns:
        a pdb file with the missing atoms and residues added
    """
    # Load the PDB file using PDBFixer
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=keep_water)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(7.4)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
    clean_protein(output_pdb, output_pdb)
    return output_pdb

def clean_protein(input_pdb, output_pdb):
    """
    Clean the protein structure after fixing or md relaxation.

    Args:
        input_pdb: input pdb file
        output_pdb: output pdb file
    
    Returns:
        a pdb file with the HETATM lines removed and H1 fixed
    """
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
        lines = [line for line in lines if not line.startswith('HETATM')]
    
    with open(output_pdb, 'w') as f:
        for line in lines:
            if line.startswith("ATOM") and int(line[22:26].strip()) == 1:
                if line[12:16].strip() == 'H':
                    # Change 'H' to 'H1' (right-aligned in a 4-character field)
                    line = line[:12] + ' H1 ' + line[16:]
            f.write(line)
    return output_pdb
    
def get_potential_energy(protein, platform='CUDA'):
    """
    Get the potential energy of the protein structure.

    Args:
        protein (str): the protein structure file name

    Returns:
        float: the potential energy of the protein structure
    """
    pdb = PDBFile(protein)
    forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    platform = Platform.getPlatformByName(platform)
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    return potential_energy.value_in_unit(kilojoules_per_mole)

def relax_backbone(protein, clean_path_fixed = None, clean_path_relaxed = None, device='CUDA', fix = True):
    """
    Relax the backbone of the protein structure to remove the clashes between the backbone atoms.

    Args:
        protein (str): the protein structure file name

    Returns:
        str: the relaxed protein structure file name
        float: the reference potential of the relaxed protein structure
    """
    clean_path_fixed = clean_path_fixed if clean_path_fixed else protein[:-4] + '_fixed.pdb'
    clean_path_relaxed = clean_path_relaxed if clean_path_relaxed else protein[:-4] + '_relaxed.pdb'
    
    if fix:
        fix_protein(protein, clean_path_fixed)
        protein = clean_path_fixed
    
    pdb = PDBFile(protein)
    modeller = Modeller(pdb.topology, pdb.positions)
    padding = 1.0*nanometer

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller.addSolvent(forcefield, model='tip3p', padding=padding)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, constraints=HBonds, nonbondedCutoff=1*nanometer)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    platform = Platform.getPlatformByName(device)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(tolerance=1e-5, maxIterations=1000)
    positions = simulation.context.getState(getPositions=True).getPositions()
    
    with open(clean_path_relaxed, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)
    
    clean_protein(clean_path_relaxed, clean_path_relaxed)
    return clean_path_relaxed
