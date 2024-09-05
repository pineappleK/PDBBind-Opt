# LP-PDBBind-Preprocess
Catching errors in ligand for LP-PDBBind

- Potential Errors in a Protein File
  - residue number not starting from 1 for a new chain [old pdbs]
  - seqres did not record missing residues - pdbfixer could not directly fix; need manual input [for instance, pdb: 6hqy]
  - pdbbind appears that the ligands interact with multiple chains whereas the ligand interacts with only one chain [6PHR]
  - pdbbind shows that the ligands interact with one chain but include multiple chains in their record [5CBM].
  - residue number not consistent across different pdbs - some starts from MET1; some with MET0
  - 


## Update 09/05

Steps to run workflow:
1. Run `biolip/create_dataset.ipynb` to download BioLiP entries and extract non-covalent small molecule binders with binding affinity annotation
2. Run `utils/process.py` to start the ligand/protein fix workflow

Problems
1. BioLiP may only annotate ONE residue with binding affinity even though the residue is part of a polymer, i.e. polypeptide(-like), polynucleotide, polysaccharide
2. The current workflow should be adapted: given a ligand id (CCD code) -> use the CONECT record to extract the polymer
3. Ligand fixer workflow should be adapted to handle polymers, mainly about get reference SMILES for polymer