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
    brufilename
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

    for line in sdffile:
        if 'V2000' in line: 
            natoms = int(line[0:3].split()[0])
            break

    n = 1
    for line in sdffile:
        atom = line.split()[3]
        atoms.append(atom)
        if n == natoms: break
        n += 1
        
    for line in sdffile:
        if 'END' in line:break
        if 'CHG' in line: break
        atom1 = int(line[0:3].split()[0]) - 1
        atom2 = int(line[3:6].split()[0]) - 1
        order = int(line[6:9].split()[0])
        bonds.append([atom1, atom2, order])

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

def molecule_to_sdf(molecule, sdffilename):
    atoms, bonds = molecule_to_atoms_bonds(molecule)
    atoms_bonds_to_sdf(atoms, bonds, sdffilename)

def atoms_bonds_to_sdf(atoms, bonds, sdffilename):
    outfile = open(sdffilename, 'w')

    outfile.write('Molecule\n pymolgen\n\n')

    n_atoms = len(atoms)
    n_bonds = len(bonds)

    outfile.write(' %s %s  0  0  1  0  0  0  0  0999 V2000\n' %(n_atoms, n_bonds))

    for atom in atoms:
        outfile.write('    0.0000    0.0000    0.0000 {0: <3} 0  0  0  0  0  0  0  0  0  0  0  0\n'.format(atom))

    for bond in bonds:
        atom1 = bond[0]
        atom2 = bond[1]
        order = bond[2]
        outfile.write('{0: >3}{1: >3}  {2}  0  0  0  0\n'.format(atom1, atom2, order))

    outfile.write('M  END\n')

    outfile.close()

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
