import shutil
from typing import List
from collections import defaultdict
import numpy as np
from rdkit import Chem
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openff.toolkit import Molecule, Topology
from openmmforcefields.generators import GAFFTemplateGenerator



def parse_headers(pdb_file):
    headers = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                break
            headers.append(line)
    return headers


def parse_missing_residues_and_atoms(headers: List[str]):
    missing_residues = []
    missing_atoms = defaultdict(list)
    read_465, read_470 = False, False
    for line in headers:
        line = line.strip()
        if line.startswith('REMARK 465   M RES C SSSEQI'):
            read_465 = True
            continue
        if line.startswith('REMARK 470   M RES CSSEQI  ATOMS'):
            read_470 = True
            continue

        if read_470:
            if not line.startswith('REMARK 470'):
                read_470 = False
            else:
                content = line.split()
                resname, chain, resid, atoms = content[2], content[3], content[4], content[5:]
                missing_atoms[chain, resid, resname] = atoms

        if read_465:
            if not line.startswith('REMARK 465'):
                read_465 = False
            else:
                content = line.split()
                resname, chain, resid = content[2], content[3], content[4]
                missing_residues.append((chain, resid, resname))
    return missing_residues, missing_atoms


import shutil

def to_quantity(ndarray):
    value = [mm.Vec3(float(arr[0]), float(arr[1]), float(arr[2])) for arr in ndarray]
    quantity = unit.Quantity(value, unit=unit.nanometers)
    return quantity


def refinePositionsWithLigand(in_protein, in_ligand, out_protein=None, out_ligand=None, extra_ffs=None, num_opt_cycles=3):
    if out_protein is None:
        shutil.copyfile(in_protein, in_protein + '.backup')
        out_protein = in_protein
        in_protein = in_protein + '.backup'
    if out_ligand is None:
        shutil.copyfile(in_ligand, in_ligand + '.backup')
        out_ligand = in_ligand
        in_ligand = in_ligand + '.backup'
    
    headers = parse_headers(in_protein)
    missing_residues, missing_atoms = parse_missing_residues_and_atoms(headers)

    mol = Chem.SDMolSupplier(in_ligand, removeHs=False)[0]
    ligand_missing_atoms = [int(x) for x in mol.GetPropsAsDict().get('Missing Atoms', '')]
    
    # Build ligand force field
    off_mol = Molecule.from_rdkit(mol)
    off_mol.assign_partial_charges('gasteiger')
    ligand_pos = to_quantity(mol.GetConformer().GetPositions() / 10)
    ligand_top = Topology.from_molecules(off_mol).to_openmm()

    extra_ffs = extra_ffs if extra_ffs is not None else []
    ff = app.ForceField('amber14-all.xml', *extra_ffs)
    ff.registerTemplateGenerator(GAFFTemplateGenerator(molecules=[off_mol], forcefield='gaff-2.11').generator)

    pdb = app.PDBFile(in_protein)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.add(ligand_top, ligand_pos)

    top = modeller.getTopology()
    pos = modeller.getPositions()

    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=0.5 * unit.nanometers, constraints=None)
    constrIndices = []
    for residue in top.residues():
        if residue.index == top.getNumResidues() - 1:
            is_missing_residue = False
            missing_atoms_residue = [atom.name for i, atom in enumerate(residue.atoms()) if i in ligand_missing_atoms]
        else:
            resinfo = (residue.chain.id, f'{residue.id}{residue.insertionCode}'.strip(), residue.name)
            missing_atoms_residue = missing_atoms.get(resinfo, [])
            is_missing_residue = resinfo in missing_residues
        
        for atom in residue.atoms():
            if (atom.element.symbol == 'H') or is_missing_residue or (atom.name in missing_atoms_residue):
                continue
            constrIndices.append(atom.index)

    posConstr = np.array([[vec.x, vec.y, vec.z] for vec in pos])[constrIndices]
    masses = [system.getParticleMass(i) for i in range(system.getNumParticles())]

    integrator = mm.LangevinMiddleIntegrator(100*unit.kelvin, 1/unit.picosecond, 0.0005*unit.picoseconds)
    simulation = app.Simulation(top, system, integrator)

    # Do optimization
    # In one cycle: full opt -> reset posisitons of constr. atoms -> constr. opt
    # This will relax steric clashes
    for _ in range(num_opt_cycles):
        simulation.context.setPositions(pos)
        simulation.minimizeEnergy()
        simulation.step(50)
        posTmp = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value
        posTmp[constrIndices] = posConstr

        for i in constrIndices:
            system.setParticleMass(i, 0.0*unit.amu)

        simulation.context.setPositions(posTmp)
        simulation.minimizeEnergy()

        pos = simulation.context.getState(getPositions=True).getPositions()
        for i in range(system.getNumParticles()):
            system.setParticleMass(i, masses[i])
    
    # output 
    fp = open(out_protein, 'w')
    fp.write(''.join(headers))
    app.PDBFile.writeModel(pdb.topology, pos[:pdb.topology.getNumAtoms()], file=fp, keepIds=True)
    app.PDBFile.writeFooter(pdb.topology, file=fp)
    fp.close()

    writer = Chem.SDWriter(out_ligand)
    for i in range(mol.GetNumAtoms()):
        vec = pos[pdb.topology.getNumAtoms() + i]
        mol.GetConformer().SetAtomPosition(i, [vec.x * 10, vec.y * 10, vec.z * 10])
    writer.write(mol)
