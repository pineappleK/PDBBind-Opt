# LP-PDBBind-Preprocess
Catching errors in ligand for LP-PDBBind

- Potential Errors in a Protein File
  - residue number not starting from 1 for a new chain [old pdbs]
  - seqres did not record missing residues - pdbfixer could not directly fix; need manual input [for instance, pdb: 6hqy]
  - pdbbind appears that the ligands interact with multiple chains whereas the ligand interact with only one chain [6PHR]
  - residue number not consistent across different pdbs - some starts from MET1; some with MET0
  - 
