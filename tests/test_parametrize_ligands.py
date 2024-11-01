
import openmm as mm
import openmm.app as app
import openmm.unit as unit

from openff.toolkit import Molecule, Topology
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from rdkit import Chem


def run(sdf):
    mol = Chem.SDMolSupplier(sdf, removeHs=False)[0]
    
    # Build ligand force field
    off_mol = Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    off_mol.assign_partial_charges('gasteiger')
    ligand_top = Topology.from_molecules(off_mol).to_openmm()
    try:
        ff = app.ForceField('amber14-all.xml')
        ff.registerTemplateGenerator(SMIRNOFFTemplateGenerator(molecules=[off_mol]).generator)
        system = ff.createSystem(ligand_top, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=0.5 * unit.nanometers, constraints=None)
    except Exception as e:
        ff = app.ForceField('amber14-all.xml')
        ff.registerTemplateGenerator(GAFFTemplateGenerator(molecules=[off_mol], forcefield='gaff-2.11').generator)
        system = ff.createSystem(ligand_top, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=0.5 * unit.nanometers, constraints=None)
    return None



if __name__ == '__main__':
    import glob
    import multiprocessing as mp
    from tqdm import tqdm

    sdfs = list(glob.glob('../raw_data_pdbbind_sm_ver2/*/*/*_ligand_fixed.sdf'))
    with mp.Pool(256) as p:
        results = list(tqdm(p.imap_unordered(run, sdfs, chunksize=1), total=len(sdfs)))