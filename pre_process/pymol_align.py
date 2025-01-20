import os
import json

from pymol import cmd
from pymol import util

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import AllChem


basedir = '/data02/venus/AlloBind/PDBBind-Opt/raw_data/allo_r2_1030'
savesdir = '/data02/venus/AlloBind/PDBBind-Opt/pre_process/allo_pro_lig_filter_formal/r6_clusters_ligand_same_align'

def pymol_align_receptor(prot1_name, prot2_name, savesdir):
    prot1_pdb = prot1_name.split('_')[0]
    prot1_rec = f"{basedir}/{prot1_pdb}/{prot1_name}/{prot1_name}_protein_refined.pdb"
    prot1_lig = f"{basedir}/{prot1_pdb}/{prot1_name}/{prot1_name}_ligand_refined.sdf"

    prot2_pdb = prot2_name.split('_')[0]
    prot2_rec = f"{basedir}/{prot2_pdb}/{prot2_name}/{prot2_name}_protein_refined.pdb"
    prot2_lig = f"{basedir}/{prot2_pdb}/{prot2_name}/{prot2_name}_ligand_refined.sdf"

    cmd.load(prot1_rec, f'{prot1_name}_rec')
    cmd.load(prot1_lig, f'{prot1_name}_lig')
    cmd.load(prot2_rec, f'{prot2_name}_rec')
    cmd.load(prot2_lig, f'{prot2_name}_lig')

    cmd.bg_color("white")
    cmd.set("cartoon_transparency", 0.3)
    cmd.color_deep("white", f'{prot1_name}_rec', 0)
    util.cba(17, f'{prot1_name}_lig', _self=cmd)
    cmd.color_deep("white", f'{prot2_name}_rec', 0)
    util.cba(9, f'{prot2_name}_lig', _self=cmd)

    cmd.align(f'{prot2_name}_rec', f'{prot1_name}_rec')
    cmd.matrix_copy(f'{prot2_name}_rec', f'{prot2_name}_lig', source_state=1, target_state=1)

    cmd.save(f'{savesdir}/{prot2_name}_aligned_ligand.sdf', f'{prot2_name}_lig')
    cmd.save(f'{savesdir}/{prot2_name}_aligned_{prot1_name}.pse')

    cmd.delete('all')
    return prot1_lig, f'{savesdir}/{prot2_name}_aligned_ligand.sdf'


with open(f'./allo_pro_lig_filter_formal/r6_clusters_ligand_same_134.json', 'r') as f:
    clusters_ligand_same = json.load(f)

for idx, cluster in clusters_ligand_same.items():
    for pair in cluster:
        prot1_lig, prot2_lig = pymol_align_receptor(pair[0], pair[1], savesdir)

        mol1 = Chem.SDMolSupplier(prot1_lig)[0]
        mol2 = Chem.SDMolSupplier(prot2_lig)[0]
        rmsd = rdMolAlign.CalcRMS(mol1, mol2)
        if rmsd > 5:
            print(pair[0], pair[1])