import os
import warnings
from collections import defaultdict
from typing import Dict, Tuple, List, Optional, Union
import xml.etree.ElementTree as ET 
import xml.dom.minidom

from textwrap import wrap
from Bio.PDB import Select

import numpy as np
import openmm as mm 
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer


aa_mapping = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR','V': 'VAL',
}


def convert_to_three_letter_seq(sequence: str):
    # Convert the one-letter code to three-letter code
    result = []
    i = 0
    while i < len(sequence):
        if sequence[i] == '(':
            # Find the closing parenthesis
            closing_index = sequence.find(')', i)
            if closing_index != -1:
                # Append the content inside the parentheses as is
                result.append(sequence[i+1:closing_index])
                i = closing_index + 1
            else:
                raise ValueError("Unmatched parenthesis in the sequence.")
        else:
            # Convert the single-letter code to three-letter code
            result.append(aa_mapping.get(sequence[i], 'XAA'))  # 'XAA' for unknown residues
            i += 1
    return result


def convert_to_seqres(sequence: List[str], chain_id: str):
    lines = wrap(' '.join([s.upper() for s in sequence]), 51)
    # Create the SEQRES lines
    seqres_lines = []
    for i, line in enumerate(lines):
        seqres_lines.append(f"SEQRES  {i+1: >2} {chain_id} {len(sequence): >4}  {line}")
    return '\n'.join(seqres_lines)


class Structure:
    def __init__(self, topology: app.Topology, positions: Union[unit.Quantity, np.ndarray], seqs: Dict[str, List[str]] = dict()):
        self.topology = topology
        self.positions = positions
        self.positions_numpy = np.array([[vec.x, vec.y, vec.z] for vec in self.positions], dtype=np.float32)
        self.seqs = seqs

    def get_positions_by_residues(self, residues):
        selected = []
        for residue in residues:
            for atom in residue.atoms():
                selected.append(atom.index)
        return self.positions_numpy[selected]
    
    def get_positions_by_atoms(self, atoms: List[app.Atom]):
        return self.positions_numpy[[atom.index for atom in atoms]]

    def select_residues(self, residues: List[app.Residue]):
        modeller = app.Modeller(self.topology, self.positions)
        to_delete = [res for res in self.topology.residues() if res not in residues]
        modeller.delete(to_delete)
        struct = Structure(modeller.topology, modeller.positions)
        return struct
    
    def select(self, select: Select):
        to_delete = []
        modeller = app.Modeller(self.topology, self.positions)
        for chain in modeller.topology.chains():
            if not select.accept_chain(chain):
                to_delete.append(chain)
                continue
            for residue in chain.residues():
                if not select.accept_residue(residue):
                    to_delete.append(residue)
                    continue
                for atom in residue.atoms():
                    if not select.accept_atom(atom):
                        to_delete.append(atom)
                        continue
        modeller.delete(to_delete)
        return Structure(modeller.topology, modeller.positions)

    def save(self, file, select: Optional[Select] = None, keepIds: bool = True, header: str = '', res_num_mapping: Optional[Dict[str, Dict[int, str]]] = None):
        if isinstance(file, str) or isinstance(file, os.PathLike):
            fp = open(file, 'w')
            self.save(fp, select, keepIds, header, res_num_mapping)
            fp.close()
        else:
            if select is not None:
                sub = self.select(select)
                top, pos = sub.topology, sub.positions
            else:
                top, pos = self.topology, self.positions

            if res_num_mapping:
                for residue in top.residues():
                    chain = residue.chain.id
                    new_res_num = [k for k, v in res_num_mapping[chain].items() if v == f"{residue.id}{residue.insertionCode}".strip()][0]
                    residue.id = str(new_res_num)
                    residue.insertionCode = " "
            
            app.PDBFile.writeHeader(top, file)
            if header:
                print(header[:-1] if header[-1] == '\n' else header, file=file)
            app.PDBFile.writeModel(top, pos, file, keepIds=keepIds)
            app.PDBFile.writeFooter(top, file)


def write_xml(element: ET.Element, fpath: os.PathLike):
    uglystr = ET.tostring(element, 'unicode', method='xml')
    pretxml = xml.dom.minidom.parseString(uglystr)
    pretstr = pretxml.toprettyxml(indent='  ').lstrip('<?xml version="1.0" ?>')
    with open(fpath, 'w') as f:
        f.write(pretstr)


def createPatch(residue, folder, amino_acid=True):
    atom_types = []
    atoms = []
    bonds = []
    nbforces = []

    for bond in residue.bonds():
        if (bond.atom1.residue is residue) and (bond.atom2.residue is residue):
            bonds.append({'atomName1': bond.atom1.name, 'atomName2': bond.atom2.name})
        elif amino_acid:
            tup = (bond.atom1.name, bond.atom2.name)
            assert tup == ('N', 'C') or tup == ('C', 'N'), f'Not a modified amino acid: {residue} with {bond}'
    
    for atom in residue.atoms():
        atype = f'{residue.name}-{atom.name}'
        ele = atom.element if atom.element is not None else app.Element.getBySymbol(atom.name[0])
        atom_types.append({'element': ele.symbol, 'name': atype, 'class': atype, 'mass': str(ele.mass._value)})
        atoms.append({'name': atom.name, 'type': atype, 'charge': '0.0'})
        nbforces.append({'type': atype, 'sigma': '0.2', 'epsilon': '0.2'})
    
    top_xml = os.path.join(folder, residue.name+'_top.xml')
    ff_xml = os.path.join(folder, residue.name+'.xml')
    
    ffElement = ET.Element('ForceField')
    atypesElement = ET.SubElement(ffElement, 'AtomTypes')
    for atype in atom_types:
        ET.SubElement(atypesElement, 'Type', attrib=atype)
    residuesElement = ET.SubElement(ffElement, 'Residues')

    for prefix in ['N', 'C', '']:
        residueElement = ET.SubElement(residuesElement, 'Residue', attrib={'name': prefix + residue.name})
        for atom in atoms:
            ET.SubElement(residueElement, 'Atom', attrib=atom)
        for bond in bonds:
            ET.SubElement(residueElement, 'Bond', attrib=bond)
        if prefix == 'N':
            ext = ['C']
        elif prefix == 'C':
            ext = ['N']
        else:
            ext = ['C', 'N']
        for name in ext:
            ET.SubElement(residueElement, 'ExternalBond', attrib={'atomName': name})

    nbElement = ET.SubElement(ffElement, 'NonbondedForce', attrib={'coulomb14scale': "0.8333333333333334", 'lj14scale': "0.5"})
    ET.SubElement(nbElement, 'UseAttributeFromResidue', attrib={'name': 'charge'})
    for nbf in nbforces:
        ET.SubElement(nbElement, 'Atom', attrib=nbf)
    write_xml(ffElement, ff_xml)

    topElement = ET.Element('Residues')
    resElement = ET.SubElement(topElement, 'Residue', attrib={'name': residue.name})
    if amino_acid:
        ET.SubElement(resElement, 'Bond', attrib={'from': '-C', 'to': 'N'})
    for bond in bonds:
        ET.SubElement(resElement, 'Bond', attrib={'from': bond['atomName1'], 'to': bond['atomName2']})
    write_xml(topElement, top_xml)
    return top_xml, ff_xml


class StandardizedPDBFixer(PDBFixer):
    """
    A class to fix standarized PDB file. Here the standardized means that:
    1. SEQRES record well documented in the PDB header
    2. PDB residues are numbered in `_pdbx_poly_seq_scheme.seq_id` with no insertion code
    """
    def __init__(self, filename: os.PathLike, verbose: bool = False, pdb_id: str = ""):
        super().__init__(filename=str(filename))
        self.mod_res_info = []
        self.missing_residues_added = []
        self.missing_atoms_added = []
        self.verbose = verbose
        self.pdb_id = pdb_id

    def log(self, *args):
        if self.verbose:
            print(*args)
    
    def findNonstandardResidues(self):
        # Find non-std residues
        super().findNonstandardResidues()
        for residue, std_res_name in self.nonstandardResidues:
            self.mod_res_info.append(
                (residue.chain.id, int(residue.id), residue.insertionCode, residue.name, std_res_name)
            )

    def findMissingResidues(self, skip_terminals: bool = True):
        sequences = {seq.chainId: seq.residues for seq in self.sequences}
        chains = list(self.topology.chains())
        
        # Find missing residues
        super().findMissingResidues()
        num_missing_residues = {chain.id: 0 for chain in chains}
        item_to_pop = []
        for info, residues in self.missingResidues.items():
            chain_index, res_id_start = info
            chain_id = chains[chain_index].id
            is_ter = (res_id_start == 0) or (num_missing_residues[chain_id] + len(residues) + res_id_start == len(sequences[chain_id]))
            if is_ter and skip_terminals:
                item_to_pop.append(info)
            else:
                offset = res_id_start + num_missing_residues[chain_id] + 1
                for i, residue in enumerate(residues):
                    self.missing_residues_added.append((chain_id, offset + i, residue))
            num_missing_residues[chain_id] += len(residues)

        for info in item_to_pop:
            self.missingResidues.pop(info)
        self.log('Missing residues:', self.missing_residues_added)
        
    def findMissingAtoms(self, skip_terminals: bool = False):
        # Find missing atoms
        super().findMissingAtoms()
        if skip_terminals:
            self.missingTerminals = {}
        
        missingAtomNames = {residue: [at.name for at in atoms] for residue, atoms in self.missingAtoms.items()}
        for residue in self.missingTerminals:
            if residue not in missingAtomNames:
                missingAtomNames[residue] = self.missingTerminals[residue]
            else:
                missingAtomNames[residue] += self.missingTerminals[residue]

        for residue, atom_names in missingAtomNames.items():            
            self.missing_atoms_added.append(
                (residue.chain.id, int(residue.id), residue.name, atom_names)
            )

        self.log('Missing atoms:', self.missing_atoms_added)

    def addMissingAtoms(self, seed=None):
        return super().addMissingAtoms(seed)
    
    @classmethod
    def _getMappedResnumWithIcode(cls, res_id: int, chain: str, res_num_mapping: Dict[str, Dict[int, str]]) -> Tuple[str, str]:
        res_id_mapped = res_num_mapping[chain][res_id]
        insert_code = ' '
        if res_id_mapped[-1].isalpha():
            insert_code = res_id_mapped[-1]
            res_id_mapped = res_id_mapped[:-1]
        return res_id_mapped, insert_code


    @classmethod
    def getFixedResidueRemarks(cls, missing_residues_info: List[Tuple[str, int, str]], res_num_mapping=None):
        if len(missing_residues_info) == 0:
            return []
        
        remark_465_lines = [
            "REMARK 465",
            "REMARK 465 FIXED MISSING RESIDUES",
            "REMARK 465 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
            "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
            "REMARK 465",
            "REMARK 465   M RES C SSSEQI"
        ]
        for chain, res_id, res_name in missing_residues_info:
            if not res_num_mapping:
                insert_code = ' '
            else:
                res_id, insert_code = cls._getMappedResnumWithIcode(res_id, chain, res_num_mapping)
            remark_465_lines.append(f"REMARK 465     {res_name:3} {chain:1}  {res_id:>4}{insert_code:>1}")
        return remark_465_lines
    
    @classmethod
    def getFixedAtomRemarks(cls, missing_atoms_info: List[Tuple[str, str, int, List[str]]], res_num_mapping=None):
        if len(missing_atoms_info) == 0:
            return []
        
        remark_470_lines = [
            "REMARK 470",
            "REMARK 470 FIXED MISSING ATOM",
            "REMARK 470 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
            "REMARK 470 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
            "REMARK 470",
            "REMARK 470   M RES CSSEQI  ATOMS"
        ]
        for chain, res_id, res_name, atoms in missing_atoms_info:
            if not res_num_mapping:
                insert_code = ' '
            else:
                res_id, insert_code = cls._getMappedResnumWithIcode(res_id, chain, res_num_mapping)
            atom_list = " ".join([f"{atom:<4}" for atom in atoms])
            remark_470_lines.append(f"REMARK 470     {res_name:3} {chain:1} {res_id:>4}{insert_code:>1}  {atom_list}")
        return remark_470_lines
    
    @classmethod
    def getModresRecords(cls, mod_res_info: List[Tuple[str, int, str, str, str]], res_num_mapping=None, pdb_id=""):
        pdb_id = 'XXXX' if not pdb_id else pdb_id.upper()
        modres_lines = []
        for chain, res_id, icode, res_name, std_res_name in mod_res_info:
            if res_num_mapping:
                res_id = res_num_mapping[chain][res_id]
                if res_id[-1].isalpha():
                    icode = res_id[-1]
                    res_id = res_id[:-1]
            # TODO: need to figure out how to deal with modfied residues with insertion code (Eric)
            modres_lines.append(f"MODRES {pdb_id:>4} {res_name:3} {chain:1} {res_id:>4}{icode:1} {std_res_name:3}   MODIFIED RESIDUE")
        return modres_lines
    
    def refineAddedAtomPositions(self, forcefield=None):
        if forcefield is None:
            forcefield = app.ForceField('amber14-all.xml', 'tip3p.xml')
        # Conver List missing atoms information to dictionary, for better indexing
        missing_atoms_info_as_dict = defaultdict(list)
        for chain, res_id, res_name, atoms in self.missing_atoms_added:
            missing_atoms_info_as_dict[(chain, res_id, res_name)] += atoms
        system = forcefield.createSystem(self.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=None, rigidWater=False)
        for residue in self.topology.residues():
            resdata = (residue.chain.id, int(residue.id), residue.name)
            if resdata in self.missing_residues_added:
                self.log(f'Found fixed residue: {residue}')
                continue
            
            for atom in residue.atoms():
                if (resdata in missing_atoms_info_as_dict) and (atom.name in missing_atoms_info_as_dict[resdata]):
                    self.log(f'Found fixed atom: {atom}')
                    continue
                if atom.element is app.element.hydrogen:
                    continue
                system.setParticleMass(atom.index, 0.0)

        integrator = mm.LangevinIntegrator(300*unit.kelvin, 10/unit.picosecond, 5*unit.femtosecond)
        context = mm.Context(system, integrator)
        context.setPositions(self.positions)
        mm.LocalEnergyMinimizer.minimize(context)
        self.positions = context.getState(getPositions=True).getPositions()
        return self.positions
    
    def runFixWorkflow(
        self, 
        output_file: os.PathLike, 
        skip_add_terminal_residues: bool = True, 
        skip_add_terminal_oxygens: bool = False, 
        add_hydrogens: bool = True, 
        refine_positions: bool = True, 
        res_num_mapping: Optional[Dict] = None
    ):
        """
        Parameters
        """
        self.findNonstandardResidues()
        self.findMissingResidues(skip_terminals=skip_add_terminal_residues)
        self.findMissingAtoms(skip_terminals=skip_add_terminal_oxygens)
        self.addMissingAtoms()
        if add_hydrogens or refine_positions:
            self.addMissingHydrogens(7.4)

        top_xmls, ff_xmls = [], []
        dirname = os.path.dirname(output_file)
        for residue, std_name in self.nonstandardResidues:
            top_xml, ff_xml = createPatch(residue, dirname)
            top_xmls.append(top_xml)
            ff_xmls.append(ff_xml)
        
        if refine_positions:
            for top_xml in top_xmls:
                app.Topology.loadBondDefinitions(top_xml)
            self.topology.createStandardBonds()
            ff = app.ForceField('amber14-all.xml', 'tip3p.xml', *list(set(ff_xmls)))
            self.refineAddedAtomPositions(ff)
        
        headers = []
        for seq in self.sequences:
            headers.append(convert_to_seqres(seq.residues, seq.chainId))
        headers += StandardizedPDBFixer.getFixedResidueRemarks(self.missing_residues_added, res_num_mapping)
        headers += StandardizedPDBFixer.getFixedAtomRemarks(self.missing_atoms_added, res_num_mapping)
        headers += StandardizedPDBFixer.getModresRecords(self.mod_res_info, res_num_mapping, self.pdb_id)

        fp = open(output_file, 'w')
        app.PDBFile.writeHeader(self.topology, fp)
        for line in headers:
            print(line, file=fp)
        # map back residue id and icode
        if res_num_mapping:
            for residue in self.topology.residues():
                res_id = res_num_mapping[residue.chain.id][int(residue.id)]
                if res_id[-1].isalpha():
                    insert_code = res_id[-1]
                    res_id = res_id[:-1]
                else:
                    insert_code = " "
                residue.id = res_id
                residue.insertionCode = insert_code

        app.PDBFile.writeModel(self.topology, self.positions, keepIds=True, file=fp)
        app.PDBFile.writeFooter(self.topology, file=fp)