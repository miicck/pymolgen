from pymolgen.molecule import Molecule, BondType
from typing import List, Tuple, TextIO
import networkx

def molecule_from_smiles(smiles: str) -> Molecule:
    """
    Loads a molecule from a smiles string.

    Parameters
    ----------
    smiles
        The smiles string to load.

    Returns
    -------
    A Molecule if successful, otherwise None
    """
    from rdkit import Chem

    rdmol: Chem.Mol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    if rdmol is None:
        return None

    graph = networkx.Graph()

    for atom in rdmol.GetAtoms():
        atom: Chem.Atom
        symbol = atom.GetSymbol()
        total_valence = atom.GetTotalValence() + atom.GetNumRadicalElectrons()
        graph.add_node(atom.GetIdx(), element=symbol, valence=total_valence)

    for bond in rdmol.GetBonds():
        bond: Chem.Bond
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond.GetBondTypeAsDouble())

    mol = Molecule()
    mol.graph = graph
    return mol


def molecule_to_smiles(mol: Molecule) -> str:
    """
    Converts a molecule to a smiles string
    Parameters
    ----------
    mol
        Molecule to convert

    Returns
    -------
    Smiles string of the molecule
    """
    from rdkit import Chem
    return Chem.MolToSmiles(molecule_to_rdkit(mol))


def graph_from_atoms_bonds(atoms: List[str], bonds: List[Tuple[int,int,int]]) -> 'Networkx graph':
    """
    Convert list of atoms and bonds to networkx graph

    Parameters
    ----------
    atoms
        list of atoms as strings
    bonds
        list of int, int, int with atom1, atom2 and bond order between the two

    Returns
    -------
    created graph if successful, otherwise None

    """

    graph = networkx.Graph()

    for n in range(len(atoms)):
        graph.add_node(n, element=atoms[n], valence=0)

    for bond in bonds:
        graph.add_edge(bond[0], bond[1], order=bond[2])

    return graph

def molecule_from_bru(bru: str) -> Molecule:
    """
    Loads this molecule from a bru format file.

    Parameters
    ----------
    asf
       Name of bru format file

    Returns
    -------
    created molecule if successful
    """

    brufilelines = bru.split("\n")

    atoms = []
    bonds = []

    for line in brufilelines:
        if 'atom' in line:
            atoms.append(line.split()[1])
        if 'bond' in line:
            atom1 = int(line.split()[1])
            atom2 = int(line.split()[2])
            order = int(line.split()[3])
            bonds.append((atom1, atom2, order))

    mol = Molecule()
    mol.graph = graph_from_atoms_bonds(atoms, bonds)
    mol.set_valence_from_bonds()

    return mol

def parse_sdf_lines(lines: List[str]) -> (List, Tuple[int, int, int]):
    """
    Parses atoms and bonds from lines of SDF file
    """

    atoms = []
    bonds = []

    for line in lines:
        if 'V2000' in line: 
            natoms = int(line[0:3].split()[0])
            nbonds = int(line[3:6].split()[0])
            break

    try:
        natoms
    except:
        raise Exception('natoms not set in parse_sdf_lines')

    try:
        nbonds
    except:
        raise Exception('nbonds not set in parse_sdf_lines')


    for line in lines[4:4+natoms]:
        atom = line.split()[3]
        atoms.append(atom)

    for line in lines[4+natoms:4+natoms+nbonds]:
        atom1 = int(line[0:3].split()[0]) - 1
        atom2 = int(line[3:6].split()[0]) - 1
        order = int(line[6:9].split()[0])
        bonds.append([atom1, atom2, order])

    return atoms, bonds    

def molecule_from_sdf(sdffilename: str) -> 'Molecule':
    """
    Loads this molecule from an SDF format file.

    Parameters
    ----------
    sdffilename
        Name of SDF format file

    Returns
    -------
    self if successful, otherwise None
    """

    atoms = []
    bonds = []

    sdffile = open(sdffilename)

    lines = sdffile.readlines()

    atoms, bonds = parse_sdf_lines(lines)

    mol = Molecule()
    mol.graph = graph_from_atoms_bonds(atoms, bonds)
    mol.set_valence_from_bonds()

    return mol

def molecule_from_sdf_lines(lines) -> 'Molecule':
    """
    Loads this molecule from an SDF format file.

    Parameters
    ----------
    lines
        Lines SDF format file

    Returns
    -------
    self if successful, otherwise None
    """

    atoms = []
    bonds = []

    atoms, bonds = parse_sdf_lines(lines)

    mol = Molecule()
    mol.graph = graph_from_atoms_bonds(atoms, bonds)
    mol.set_valence_from_bonds()

    return mol

def molecule_from_sdf_large(sdffilename: str, start_line: int) -> 'Molecule':
    """
    Loads molecule i from a large SDF file
    """
    lines = []

    with open(sdffilename) as infile:

        n = 0
        for line in infile:
            if n >= start_line:
                lines.append(line)

                if '$$$$' in line:
                    break
            n += 1

    atoms, bonds = parse_sdf_lines(lines)

    mol = Molecule()
    mol.graph = graph_from_atoms_bonds(atoms, bonds)
    mol.set_valence_from_bonds()

    return mol    

def molecule_to_atoms_bonds(molecule: Molecule) -> (List, Tuple[int, int, int]):
    """
    Writes this molecule to an SDF format file.

    Parameters
    ----------
    molecule
        Molecule object
    sdffilename
        Name of SDF format file

    Returns
    -------
    None, creates file
    """    

    atoms = []
    bonds = []

    ids = {i: n for n, i in enumerate(molecule.graph.nodes)}

    for i in molecule.graph.nodes:
        atom = molecule.graph.nodes[i]["element"]
        atoms.append(atom)

    n_edges = len(list(molecule.graph.edges))

    n = 0
    for i in molecule.graph.edges:
        atom1 = ids[list(molecule.graph.edges)[n][0]] + 1 
        atom2 = ids[list(molecule.graph.edges)[n][1]] + 1 
        order = molecule.graph.edges[i]["order"]
        bonds.append([atom1, atom2, order])
        n += 1

    return atoms, bonds

def molecule_to_sdf(molecule):
    atoms, bonds = molecule_to_atoms_bonds(molecule)
    lines = atoms_bonds_to_sdf(atoms, bonds)
    return lines

def atoms_bonds_to_sdf(atoms, bonds):
    lines = []
    lines.append('Molecule\n pymolgen\n\n')

    n_atoms = len(atoms)
    n_bonds = len(bonds)

    lines.append('{0: >3}{1: >3}  0  0  1  0  0  0  0  0999 V2000\n'.format(n_atoms, n_bonds))

    for atom in atoms:
        lines.append('    0.0000    0.0000    0.0000 {0: <3} 0  0  0  0  0  0  0  0  0  0  0  0\n'.format(atom))

    for bond in bonds:
        atom1 = bond[0]
        atom2 = bond[1]
        order = bond[2]
        lines.append('{0: >3}{1: >3}  {2}  0  0  0  0\n'.format(atom1, atom2, order))

    lines.append('M  END\n')

    return lines

def molecule_to_rdkit(molecule: Molecule) -> 'Chem.RWMol':
    """
    Converts to an RDkit molecule
    """
    from rdkit import Chem
    mol = Chem.RWMol()

    graph_to_export = molecule.graph.subgraph(i for i in molecule.graph if molecule.graph.nodes[i]["element"] != "H")
    if len(graph_to_export) == 0:
        graph_to_export = molecule.graph

    for i in graph_to_export:
        atom = Chem.Atom(graph_to_export.nodes[i]["element"])
        mol.AddAtom(atom)

    ids = {i: n for n, i in enumerate(graph_to_export.nodes)}
    for i in graph_to_export:
        for j in graph_to_export[i]:
            if i >= j:
                continue  # Avoid double counting

            bond_type = BondType.from_order(graph_to_export[i][j]["order"])

            mol.AddBond(ids[i], ids[j], {
                BondType.SINGLE: Chem.BondType.SINGLE,
                BondType.AROMATIC: Chem.BondType.AROMATIC,
                BondType.DOUBLE: Chem.BondType.DOUBLE,
                BondType.TRIPLE: Chem.BondType.TRIPLE
            }[bond_type])

    return mol

def sdf_to_asf(sdffilename: str, asffilename: str):
    """
    Converts multi-molecule SDF file to ASF file
    """
    if sdffilename == asffilename:
        raise Exception("sdffilename same as asffilename")

    infile = open(sdffilename)
    
    outfile = open(asffilename, 'w') 
    outfile.close()

    lines = []

    for line in infile:
        lines.append(line)
        if '$$$$' in line:
            atoms, bonds = parse_sdf_lines(lines)

            atoms_bonds_to_asf(atoms, bonds, asffilename)
            lines = []

    infile.close()

def atoms_bonds_to_asf(atoms, bonds, asffilename):

    with open(asffilename, 'a') as outfile:

        outfile.write('atoms')

        for atom in atoms:
            outfile.write(' %s' %atom)

        outfile.write('\n')

        outfile.write('bonds')

        for bond in bonds:
            atom1 = bond[0]
            atom2 = bond[1]
            order = bond[2]
            outfile.write(' %s %s %s' %(atom1, atom2, order))

        outfile.write('\n')

def read_asf_lines(lines: List[str]) -> Molecule:
    """
    Loads this molecule from am ASF format file.

    Parameters
    ----------
    asf
       Name of ASF format file

    Returns
    -------
    created molecule if successful
    """

    bonds = []

    for line in lines:
        if 'atoms' in line:
            atoms = line.split()[1:]
        if 'bonds' in line:

            n = 0
            for i in line.split()[1:]:
                if n%3 == 0:
                    atom1 = int(i)
                if n%3 == 1:
                    atom2 = int(i)
                if n%3 == 2:
                    order = int(i)
                    bonds.append((atom1, atom2, order))
                n += 1

    mol = Molecule()
    mol.graph = graph_from_atoms_bonds(atoms, bonds)
    mol.set_valence_from_bonds()

    return mol

def read_asf_file(asf):

    with open(asf) as infile:
        lines = infile.readlines()

    for i in range(len(lines)):
        if i%2 == 0:
            mol = read_asf_lines(lines[i:i+2])
            smi = molecule_to_smiles(mol)
            print(smi)