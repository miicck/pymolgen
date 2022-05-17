from pymolgen.molecule import Molecule, BondType
from typing import List, Tuple
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

    atom_valences = {
    'H':1,
    'C':4,
    'N':3,
    'O':2,
    'S':2,
    'F':1,
    'Cl':1,
    'Br':1
    }

    for n in range(len(atoms)):
        atom = atoms[n]
        total_valence = atom_valences[atom]
        print(n, atom, total_valence)
        graph.add_node(n, element=atom, valence=total_valence)

    for bond in bonds:
        atom1 = bond[0]
        atom2 = bond[1]
        order = bond[2]
        graph.add_edge(atom1, atom2, order=order)

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
    
    return mol

def molecule_from_sdf(sdffilename: str) -> 'Molecule':
    """
    Loads this molecule from an SDF format file.

    Parameters
    ----------
    sdffilename
        Name of bru format file

    Returns
    -------
    self if successful, otherwise None
    """

    atoms = []
    bonds = []

    sdffile = open(sdffilename)

    for line in sdffile:
        if 'V2000' in line: 
            natoms = int(line.split()[0])
            break

    n = 1
    for line in sdffile:
        print(line)
        atom = line.split()[3]
        atoms.append(atom)
        if n == natoms: break
        n += 1
        
    for line in sdffile:
        if 'END' in line: break
        print(line)
        atom1 = int(line.split()[0]) - 1
        atom2 = int(line.split()[1]) - 1
        order = int(line.split()[2])
        bonds.append([atom1, atom2, order])

    mol = Molecule()
    mol.graph = graph_from_atoms_bonds(atoms, bonds)

    return mol


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
