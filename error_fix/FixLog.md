## Note 09/27
Eric: I manually check some errors

## Discard
+ 1d6s: L-Met covalently bind to the cofactor and Kd is measure on wt OASS (the structure is K41A)
+ 2xag, 2xaf, 2xah, 2xaq, 2xas: covalently bound to FAD
+ 3zju, 3zjt, 3zjv: Covalent to DNA + Boron contained
+ 1t7d: Missing atoms in M12 residue, should have aliphatic chains
+ 3krd: Missing atoms in HXD
+ 3m3r: Two 7-mer (7 GLC), not a typical protein-ligand binding
+ 4ofl: multiple ligands binding to metal
+ 5mby: multiple ligands binding to metal
+ 6a80: multiple cystidines binding
+ 3rzi: mutliple ligands, Kd=21+/-1uM (Phe only); Kd=4.7+/-0.1uM (Trp only); Kd=5.0+/-0.1uM (Phe to Trp-bound protein); Kd=1.08+/-0.03uM (Trp to Phe-bound protein)

## Manual fix
+ 3kyf, 3kyg: dna-like, two chains should be considered as one
+ 1ozv, 1ssq, 1xt8, 1yxd, 2nwl, 3f3a, 3qs4: Ligands labeled as the same chain
+ 2yln: disulfide bond is (partially) formed between to CYS
+ 2h2e: Wrong bonds in PDB
+ 3b3s: Wrong bonds in PDB
+ 4ykj: Wrong bonds in PDB, ligand is ALA
+ 4ykk: Wrong bonds in PDB, ligand is DSN
+ 5vdk: Wrong bonds in PDB, ligand is 8X7
+ 6mle: Wrong bonds in PDB, ligand is ARG
+ 6mlo: Wrong bonds in PDB, ligand is ARG
+ 6qpl: Wrong bonds in PDB, ligand is JC5
+ 4ay6: Wrong bonds in PDB, ligand is 12V
+ 4ay6: Wrong bonds in PDB, ligand is FYB
+ 6eiz: Wrong ligand name, ligand is FOF
+ 1icj: Wrong ligand name, should be 2PE

+ 2r1w: 2-mer (KDA-KDB)
+ 2ri9: 2-mer (NMA-LDY)
+ 2w87: 2-mer (GCU-UNL)
+ 2zg3: 3-mer (BGC-GAL-SIA)
+ 3x00: 3-mer (ZOF-EDN-ZOF)
+ 4nrl: 5-mer (BGC-SIA)
+ 5htb: n-mer (ARKKQT-NH2-6L5)
+ 5htc: n-mer (ARKKQT-66N-66M)
+ 5may: 2-mer (FUL-PK6)
+ 5mb1: 2-mer (FUL-7KT)
+ 5o4z: 2-mer (A-A)
+ 6hck: 2-mer (LEU-LEU)
+ 2k46: 2-mer (GLC-GLC)
+ 2r0h: 3-mer (NAG-NAG-NAG)
+ 3ehn: 2-mer (NDG-GAL)
+ 5gwy: 6-mer
+ 5yqw: 2-mer (NAG-NAG)

+ 1fls: Bad box
+ 1g42, 1gnm, 1gnn, 1gno, 1hps, 1htg, 1j1a, 4phv: ligand multiple occupancies
+ 1i6v, 3fi2, 3fv8, 5ayf: Bad PDB, UNK
+ 1yhm, 2hj4, 2hjb, 2q7q, 3au6, 5lub, 5lzh, 5ogl, 6h9v: Code problem, modres patch
+ 3ary, 3arz, 3as3, 5eie: Code problem, rdkit fail to generate smiles
+ 4hwb, 5ct1: Ligand not found - NHE?
+ 4mj5: No PDB
+ 5ghv, 5knj: rdkit fail
+ 6f4x: PDBFixer failed

