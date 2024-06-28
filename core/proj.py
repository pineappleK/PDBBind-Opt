"""
A general project class that contains the project information and the project directory.
"""

import os, shutil, sys, json
sys.path.append('/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/')
from utils.prepare_protein import process_protein
from utils.prepare_ligand import process_ligand

class Protein:
    """
    A class to store key information about a protein.
    """
    def __init__(self, identifier, pdb_seq, pdb_file, is_apo, af2_path, clean_path, keep_water, add_hydrogens, get_af2_structure):
        self.identifier = identifier
        self.is_apo = is_apo
        self.fail_to_process = True
        self.af2_path = None
        self.clean_path = None
        self.protein_seq = None
        result = self._process_protein(pdb_seq=pdb_seq, pdb_file=pdb_file, af2_path=af2_path, clean_path=clean_path, keep_water=keep_water, add_hydrogens=add_hydrogens, get_af2_structure=get_af2_structure)
        if result is not None:
            self.af2_path, self.clean_path, self.protein_seq = result
            self.fail_to_process = False
    
    def _process_protein(self, pdb_seq, pdb_file, af2_path, clean_path, keep_water, add_hydrogens, get_af2_structure):
        """
        clean up the protein by removing water molecules and other unwanted residues.
        """
        try:
            result = process_protein(pdb_seq = pdb_seq,
                        pdb_file = pdb_file,
                        af2_path = af2_path,
                        clean_path = clean_path,
                        keep_water = keep_water,
                        add_hydrogens = add_hydrogens,
                        get_af2_structure = get_af2_structure)
            af2_path, clean_path, seq = result
            return af2_path, clean_path, seq
        except Exception as e:
            print(e)
            print(f"Failed to process the protein {self.identifier}.")
            return None


class Ligand:
    """
    A class to store key information about a ligand.
    """
    def __init__(self, identifier, smiles=None, structural_file_path=None, sdf_out_path=None, calc_charges=False, keep_coords=True, hydrogens=True):
        self.identifier = identifier
        self.fail_to_process = True
        self.canonical_smiles = None
        self.sdf_out_path = None
        result = self._process_ligand(smiles, structural_file_path, sdf_out_path=sdf_out_path, calc_charges=calc_charges, keep_coords=keep_coords, hydrogens=hydrogens)
        if result is not None:
            self.canonical_smiles, self.sdf_out_path = result
            self.fail_to_process = False
    
    def _process_ligand(self, smiles, file_path, sdf_out_path=None, calc_charges=False, keep_coords=True, hydrogens=True):
        """
        clean up ligand input files
        """
        try:
            result = process_ligand(smiles=smiles, 
                                file_path=file_path, 
                                sdf_out_path=sdf_out_path, 
                                calc_charges=calc_charges, 
                                keep_coords=keep_coords, 
                                hydrogens=hydrogens)
            canonical_smiles, sdf_out_path = result
            return canonical_smiles, sdf_out_path
        except Exception as e:
            print(e)
            print(f"Failed to process the ligand {self.identifier}.")
            return None

class Project:
    """
    A project that safely stores the protein and ligand information for downstream tasks.
    """
    def __init__(self, proj_dir, proj_name):
        self.proj_dir = proj_dir
        self.create_proj_dir()
        self.proj_name = proj_name
        self.proj_info_file = os.path.join(self.proj_dir, self.proj_name + '.json')
        self.proj_info = self.load_proj_info()
        self.proteins = []
        self.ligands = []

    def load_proj_info(self):
        if os.path.exists(self.proj_info_file):
            with open(self.proj_info_file, 'r') as f:
                proj_info = json.load(f)
            return proj_info
        else:
            return {
                'proj_name': self.proj_name,
                'proj_dir': self.proj_dir,
                'proteins': {},
                'ligands': {}
            }
    
    def save_proj_info(self):
        with open(self.proj_info_file, 'w') as f:
            json.dump(self.proj_info, f, indent=4)
    
    def create_proj_dir(self):
        # store the proteins and ligands in separate directories
        if not os.path.exists(self.proj_dir):
            os.makedirs(self.proj_dir)
        if not os.path.exists(os.path.join(self.proj_dir, 'proteins')):
            os.makedirs(os.path.join(self.proj_dir, 'proteins'))
        if not os.path.exists(os.path.join(self.proj_dir, 'ligands')):
            os.makedirs(os.path.join(self.proj_dir, 'ligands'))
    
    def remove_proj_dir(self):
        if os.path.exists(self.proj_dir):
            shutil.rmtree(self.proj_dir)
    
    def add_info(self, category, key, value):
        self.proj_info[category][key] = value
        self.save_proj_info()
    
    def remove_info(self, category, key):
        if key in self.proj_info[category]:
            del self.proj_info[category][key]
            self.save_proj_info()
    
    def get_info(self, category, key):
        return self.proj_info[category].get(key, None)
    
    def get_all_info(self):
        return self.proj_info

    def add_protein(self, pdb_seq, pdb_path, is_apo, get_af2_structure, identifier=None, keep_water=False, add_hydrogens=True):
        """
        Adds a Protein instance to the project.
        """
        if identifier is None:
            identifier = f'protein_{len(self.proteins) + 1}'
            if pdb_path is not None:
                pdb_file_name = os.path.basename(pdb_path).split('.')[0]
                identifier = f'protein_{len(self.proteins) + 1}_{pdb_file_name}'
        
        af2_path = None
        clean_path = None
        
        if pdb_path is not None:
            clean_path = os.path.join(self.proj_dir, 'proteins', f'{identifier}.pdb')
            if get_af2_structure:
                af2_path = os.path.join(self.proj_dir, 'proteins', f'{identifier}_af2_aligned.pdb')
        
        elif pdb_seq is not None:
            af2_path = os.path.join(self.proj_dir, 'proteins', f'{identifier}_af2.pdb')
    
        protein = Protein(identifier, pdb_seq, pdb_path, is_apo, af2_path, clean_path, keep_water, add_hydrogens, get_af2_structure)
        self.proteins.append(protein)
        
        self.add_info('proteins', f'{protein.identifier}', {
            'failed': protein.fail_to_process,
            'is_apo': protein.is_apo,
            'af2_path': protein.af2_path,
            'clean_path': protein.clean_path,
            'protein_seq': protein.protein_seq
        })

    def get_protein_info(self, identifier):
        """
        Retrieves a Protein instance by its identifier.
        """
        return self.get_info('proteins', identifier)

    def add_ligand(self, smiles, file_path, identifier=None, calc_charges=False, keep_coords = True, hydrogens=True):
        """
        Adds a Ligand instance to the project.
        """
        if identifier is None:
            identifier = f'ligand_{len(self.ligands) + 1}'
            if file_path is not None:
                file_name = os.path.basename(file_path).split('.')[0]
                identifier = f'ligand_{len(self.ligands) + 1}_{file_name}'
        
        sdf_out_path = os.path.join(self.proj_dir, 'ligands', f'{identifier}.sdf')
        ligand = Ligand(identifier, smiles, file_path, sdf_out_path, calc_charges, keep_coords, hydrogens)
        self.ligands.append(ligand)
        
        self.add_info('ligands', f'{ligand.identifier}', {
            'failed': ligand.fail_to_process,
            'canonical_smiles': ligand.canonical_smiles,
            'sdf_out_path': ligand.sdf_out_path
        })
        
    def get_ligand_info(self, identifier):
        """
        Retrieves a Ligand instance by its identifier.
        """
        return self.get_info('ligands', identifier)

if __name__ == '__main__':
    proj = Project('test_proj', 'test_proj')
    proj.add_protein(pdb_seq = None, pdb_path = "/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/data/test/1q72_protein.pdb", is_apo = True, get_af2_structure = False)
    proj.add_ligand(smiles = None, file_path = "/pscratch/sd/k/kysun/apo-holo-project/EvoStruct/data/test/1q72_ligand.sdf")
    proj.add_ligand(smiles = "CCO", file_path = None)
    proj.add_ligand(smiles = "C1CCC", file_path = None)