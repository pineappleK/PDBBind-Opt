"""
Prepare the ligand for energy minimization process

Author: Oliver Sun
Date: 05/10/2024
"""
import warnings
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import ToolkitRegistry, AmberToolsToolkitWrapper
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point3D

def read_molecule(molecule_file, sanitize=True, calc_charges=False):
    if molecule_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=True)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=True)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=True)
    else:
        raise ValueError('Expect the format of the molecule_file to be '
                         'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))
    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)
        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                warnings.warn('Unable to compute charges for the molecule.')
    except Exception as e:
        print(e)
        print("RDKit was unable to read the molecule.")
        return None
    return mol
    
def process_ligand(smiles=None, file_path=None, sdf_out_path=None, calc_charges=False, keep_coords = True, hydrogens=True):
    """
    Process ligand input and return canonical SMILES and save SDF file with conformation optimization.
    
    Args:
        smiles (str): Input SMILES string (optional).
        file_path (str): Input structural file path (optional).
        sdf_out_path (str): Output path for the SDF file.
        coord_file_path (str): Optional file path to get coordinates from.
    
    Returns:
        str: Canonical SMILES string.
    """
    mol_source = None
    if not smiles and not file_path:
        raise ValueError("At least one of SMILES or file path must be provided.")
    
    mol_from_smiles = mol_from_file = None

    if smiles:
        mol_from_smiles = Chem.MolFromSmiles(smiles)
        if mol_from_smiles is None:
            raise ValueError("Invalid SMILES string provided.")
        Chem.RemoveHs(mol_from_smiles)
    
    if file_path:
        mol_from_file = read_molecule(file_path, calc_charges=calc_charges)
        if mol_from_file is None:
            raise ValueError("Failed to read molecule from file.")

    if not mol_from_smiles and not mol_from_file:
        raise ValueError("Failed to read molecule from file or SMILES string.")
    
    elif mol_from_smiles and mol_from_file:
        smiles_from_file = Chem.MolToSmiles(mol_from_file, canonical = True)
        canonical_smiles = Chem.MolToSmiles(mol_from_smiles, canonical = True)
        if canonical_smiles != smiles_from_file:
            print("Conflict detected: Using SMILES input over file input.")
            mol_source = 'smiles'
            mol = mol_from_smiles
        else:
            mol = mol_from_file
            mol_source = 'file'
    
    elif mol_from_smiles:
        canonical_smiles = Chem.MolToSmiles(mol_from_smiles, canonical = True)
        mol = mol_from_smiles
        mol_source = 'smiles'

    else:
        canonical_smiles = Chem.MolToSmiles(mol_from_file, canonical = True)
        mol = mol_from_file
        mol_source = 'file'

    if keep_coords:
        if mol_from_file and mol_source == 'file':
            save_molecule_as_sdf(mol, sdf_out_path, skip_optimization=True, hydrogens=hydrogens)
    else:
        mol = Chem.MolFromSmiles(canonical_smiles)
        save_molecule_as_sdf(mol, sdf_out_path, hydrogens=hydrogens)

    return canonical_smiles, sdf_out_path

def save_molecule_as_sdf(mol, output_path, skip_optimization = False, hydrogens=True):
    """Save a molecule as an SDF file with conformation optimization."""
    if hydrogens:
        mol = Chem.AddHs(mol)
    w = Chem.SDWriter(output_path)
    if not skip_optimization:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Generate 3D coordinates
        AllChem.UFFOptimizeMolecule(mol)  # Optimize conformation
    w.write(mol)
    w.close()

def prepare_ligand_for_md(ligand_sdf_path, ligand_final_path = None, charge_model = 'am1bcc'):
    """
    Do the heavy calculation of partial charge first and reuse it for batch jobs.

    Args:
        ligand_sdf_path: path to the ligand sdf file
        ligand_final_path: path to write the final ligand pdb file
    """
    if not ligand_final_path:
        ligand_final_path = ligand_sdf_path
    # assuming here that the sdf file contains only the conformer that we are going to work with
    registry = ToolkitRegistry()
    registry.register_toolkit(AmberToolsToolkitWrapper())
    ligand = Molecule.from_file(ligand_sdf_path)
    ligand.assign_partial_charges(charge_model, toolkit_registry=registry)
    ligand.to_file(ligand_final_path, file_format='sdf')

def minimize_protein_ligand(protein_pdb_path, ligand_sdf_path, final_complex_path = None, platform = 'CUDA', fix_ligand = False):
    """
    Prepare and minimize the energy of a protein-ligand complex
    
    Args:
        protein_pdb_path: path to the protein pdb file
        ligand_sdf_path: path to the ligand sdf file
        final_complex_path: path to write the final complex pdb file
        platform: platform to run simulations on
    """
    pdb = PDBFile(protein_pdb_path)
    ligand = Molecule.from_file(ligand_sdf_path)
    ligand_topology = ligand.to_topology()
    
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.add(ligand_topology.to_openmm(), ligand_topology.get_positions().to_openmm())
    
    # Create a force field with GAFF for ligands and AMBER ff14SB for proteins
    forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    gaff = GAFFTemplateGenerator(molecules=ligand, forcefield = "gaff-2.11")
    forcefield.registerTemplateGenerator(gaff.generator)

    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    if fix_ligand:
        # Get ligand atom indices in the combined topology
        ligand_atom_indices = list(ligand_topology.atoms)
        ligand_indices = list(range(len(pdb.positions), len(pdb.positions) + len(ligand_atom_indices)))
        
        # Set the mass of each ligand atom to zero
        for idx in ligand_indices:
            system.setParticleMass(idx, 0.0)
    
    integrator = LangevinIntegrator(300*kelvin, 1/picoseconds, 0.002*picoseconds)

    platform = Platform.getPlatformByName(platform)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    
    simulation.minimizeEnergy(tolerance=1e-5, maxIterations=1000)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    print(f"Finish minimizing, the final potential energy is: {energy}")
    
    if final_complex_path:
        with open(final_complex_path, 'w') as outfile:
            PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), outfile)
    return energy.value_in_unit(kilojoules_per_mole)

# Example usage:
# energy = simulate_protein_ligand('path/to/protein.pdb', 'path/to/ligand.sdf')
if __name__ == "__main__":
    #prepare_ligand("data/cdd_1845.sdf")
    minimize_protein_ligand("/pscratch/sd/k/kysun/mcsce-precompute/data/mpro_1845_relaxed.pdb",
                            "/pscratch/sd/k/kysun/mcsce-precompute/data/cdd_1845.sdf")