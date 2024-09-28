import os, glob
from collections import defaultdict


if __name__ == '__main__':
    dataset_dir = '../raw_data_pdbbind_poly'
    report_file = 'PDBBind-poly.md'
    
    err_types = defaultdict(list)
    for err in glob.glob(os.path.join(dataset_dir, '*/err')):
        pdbid = os.path.basename(os.path.dirname(err))
        with open(err) as f:
            last = f.read().strip().split('\n')[-1]
        if last.startswith('AssertionError: Number of ligand residues'):
            err_types['Covalent'].append(pdbid)
        elif last.startswith('RuntimeError: Rare element in ligand'):
            err_types['Ligand with rare elements'].append(pdbid)
        elif last.startswith('AssertionError: Steric clash'):
            err_types['Steric clash'].append(pdbid)
        elif last.startswith('fix_ligand.LigandFixException: Number of atoms not match'):
            err_types['Fail to match template SMILES'].append(pdbid)
        elif last.startswith('fix_ligand.LigandFixException: No reference found'):
            err_types['No reference'].append(pdbid)
        elif last.startswith('fix_ligand.LigandFixException:'):
            err_types['Fail to fix ligand'].append(pdbid)
        elif last.startswith('ValueError: No template found') or last.startswith('AssertionError: Not a modified amino acid'):
            err_types['Fail to fix protein'].append(pdbid)
        elif last.startswith('AssertionError: No ligands found'):
            err_types['Fail to find ligand'].append(pdbid)
        else:
            err_types['Others'].append(pdbid)
    
    for key, value in err_types.items():
        print(key, len(value))
    

    with open(report_file, 'w') as f:
        for key, value in err_types.items():
            value.sort()
            f.write(f'## {key} ({len(err_types[key])})\n')
            for val in value:
                f.write(f'+ {val}\n')
            f.write('\n')
