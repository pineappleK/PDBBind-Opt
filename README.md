# LP-PDBBind-Preprocess
Catching errors in ligand for LP-PDBBind

- Potential Errors in a Protein File
  - residue number not starting from 1 for a new chain [old pdbs]
  - seqres did not record missing residues - pdbfixer could not directly fix; need manual input [for instance, pdb: 6hqy]
  - pdbbind appears that the ligands interact with multiple chains whereas the ligand interacts with only one chain [6PHR]
  - pdbbind shows that the ligands interact with one chain but include multiple chains in their record [5CBM].
  - residue number not consistent across different pdbs - some starts from MET1; some with MET0
  - 
