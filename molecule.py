import random
import time
import networkx
from copy import deepcopy
from typing import Iterator, Dict, Optional, Tuple, Set, List
from rdkit import Chem
from rdkit.Chem import Draw
import multiprocessing
from enum import Enum
from pymolgen.bond_generator import BondGenerator


class BondType(Enum):
    SINGLE = 1
    AROMATIC = 2
    DOUBLE = 3
    TRIPLE = 4

    @staticmethod
    def from_order(order: float) -> 'BondType':
        half_orders = round(order * 2)
        return {
            2: BondType.SINGLE,
            3: BondType.AROMATIC,
            4: BondType.DOUBLE,
            6: BondType.TRIPLE
        }[half_orders]


class FractionalOrderException(Exception):
    def __init__(self, message="Fractional order detected"):
        super().__init__(message)


def is_integer_order(order: float) -> bool:
    return abs(order - round(order)) < 1e-3


def to_integer_order(order: float) -> int:
    if not is_integer_order(order):
        raise FractionalOrderException
    return round(order)


class Molecule:

    ##################
    # Dunder methods #
    ##################

    def __init__(self, allow_frac_order: bool = True):
        """
        Build an empty molecule
        """
        self._graph: Optional[networkx.Graph] = None
        self._allow_frac_order: bool = allow_frac_order

    def __str__(self):
        return Chem.MolToSmiles(self.to_rdkit())

    ##############
    # Properties #
    ##############

    @property
    def graph(self) -> networkx.Graph:
        if self._graph is None:
            self.graph = networkx.Graph()
        return self._graph

    @graph.setter
    def graph(self, val: networkx.Graph):

        if networkx.is_frozen(val):
            raise Exception("Tried to freeze graph!")

        # Check resulting graph has everything we need
        for i in val.nodes:
            assert val.nodes[i]["valence"] is not None
            assert val.nodes[i]["element"] is not None

        for i in val.nodes:
            for j in val[i]:
                assert val[i][j]["order"] > 0

        self._graph = val

    @property
    def atom_count(self) -> int:
        return len(self.graph)

    @property
    def total_free_valence(self) -> int:
        return sum(self.free_valence(i) for i in self.graph)

    @property
    def attach_points(self) -> List[int]:
        """
        Returns
        -------
        List of points where this molecule fragment
        could attach to another molecule fragment.
        """
        points = []
        for i in self.graph:
            v = self.free_valence(i)
            if is_integer_order(v) and to_integer_order(v) > 0:
                points.append(i)
        return points

    ###########
    # Methods #
    ###########

    def total_order_of_bonds(self, i: int) -> float:
        """
        Returns
        -------
        The total order of bonds connected to atom i
        """
        return sum(self.graph[i][j]["order"] for j in self.graph[i])

    def free_valence(self, i: int) -> float:
        """
        Returns
        -------
        The unbonded valence of the atom i
        """
        return max(self.graph.nodes[i]["valence"] - self.total_order_of_bonds(i), 0.0)

    def copy(self) -> 'Molecule':
        """
        Returns
        -------
        A deep copy of this molecule
        """
        return deepcopy(self)

    def load_smiles(self, smiles: str) -> 'Molecule':
        """
        Loads this molecule from a smiles string.

        Parameters
        ----------
        smiles
            The smiles string to load.

        Returns
        -------
        self if successful, otherwise None
        """

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

        self.graph = graph
        return self

    def get_fragment(self, nodes: Set[int], allow_breaking_aromatic=False) -> 'Molecule':
        """
        Returns the fragment of this molecule, containing
        only the specified atom indices. Will break bonds
        in order to create the fragment, remembering the
        valence that is freed up as a result.

        Parameters
        ----------
        nodes
            The nodes to keep in the fragment.

        Returns
        -------
            The resulting molecular fragment.
        """
        frag = Molecule()
        frag.graph = networkx.Graph(self.graph.subgraph(nodes))

        # Look for broken bonds
        for i in nodes:
            for n in self.graph[i]:
                if n in nodes:
                    continue

                if not allow_breaking_aromatic:
                    # Don't allow breaking of aromatic bonds
                    order = self.graph[i][n]["order"]
                    if abs(order - round(order)) > 0.1:
                        return None

        return frag

    def random_fragment(self, min_size: int = 1, max_size: int = None) -> 'Molecule':
        """
        Returns a random fragment (sub-molecule) of this molecule.
        Cleaves bonds and remembers the freed-up valence.

        Parameters
        ----------
        min_size
            The minimum size of the fragment to generate.
        max_size
            The maximum size of the fragment to gnerate
            (will be capped at the size of the molecule - 1).
        Returns
        -------
        fragment
            A fragment of this molecule, with apropritate
            free valence points along cleaved bonds
        """
        # Pick the size of the fragment
        if max_size is None or max_size > self.atom_count - 1:
            max_size = self.atom_count - 1
        size = random.randint(min_size, max_size)

        # Pick a random starting location
        i = random.choice(list(self.graph.nodes))
        nodes = {i}

        # Add a random neighbouring atom until
        # the fragment is of the desired size
        while len(nodes) < size:

            neighbours = set()
            for n in nodes:
                neighbours.update(self.graph[n])
            neighbours -= nodes

            if len(neighbours) == 0:
                break
            nodes.add(random.choice(list(neighbours)))

        frag = self.get_fragment(nodes)
        if frag is None:
            # Fragment was not sensible, try again
            return self.random_fragment(min_size=min_size, max_size=max_size)
        return frag

    def plot(self, timeout: float = None):
        """
        Show a plot of this molecule.
        """

        def plot_on_thread():
            Draw.ShowMol(self.to_rdkit(), size=(1024, 1024))

        if timeout is None:
            plot_on_thread()
        else:
            p = multiprocessing.Process(target=plot_on_thread)
            p.start()
            time.sleep(timeout)
            p.terminate()

    def make_disjoint(self, other: 'Molecule') -> None:
        """
        Relabel the atoms in this molecule so
        they don't conflict with the atoms in other.

        Parameters
        ----------
        other
            Molecule that we want to make disjoint to
        """
        i_new = max(list(other.graph.nodes) + list(self.graph.nodes)) + 1

        for i in self.graph.nodes:
            if i in other.graph.nodes:
                self.graph = networkx.relabel_nodes(self.graph, {i: i_new})
                i_new += 1

        assert len(set(self.graph.nodes).intersection(set(other.graph.nodes))) == 0

    def remove_random_atom(self, element: str = None) -> 'Molecule':
        """
        Remove a random leaf atom (only bonded to one other atom) from the molecule.

        Parameters
        ----------
        element
            Optionally constrain the element removed.

        Returns
        -------
            The molecular fragment, after removing the atom.
        """
        # Get the list of possible nodes to remove (without breaking the molecule apart)
        nodes = list(n for n in self.graph.nodes if len(self.graph[n]) == 1)
        if element is not None:
            nodes = [n for n in nodes if self.graph.nodes[n]["element"].upper() == element.upper()]
        if len(nodes) == 0:
            return False
        to_remove = random.choice(nodes)
        keep_nodes = set(self.graph.nodes) - {to_remove}
        return self.get_fragment(keep_nodes)

    def saturate_single_aromatic_fragment(self) -> bool:

        for i in self.graph.nodes:
            for j in self.graph[i]:

                order = self.graph[i][j]["order"]
                if not is_integer_order(order):

                    if self.free_valence(i) > 0.5 and self.free_valence(j) > 0.5:
                        self.graph[i][j]["order"] += 0.5
                        return True

        return False

    def hydrogenate(self, saturate_aromatic_fragments: bool = True) -> int:
        """
        Terminates all free valence in this molecule with hydrogen
        """

        if saturate_aromatic_fragments:
            while self.saturate_single_aromatic_fragment():
                pass

        added = 0
        offset = max(self.graph.nodes)

        for i in list(self.graph):

            v = self.free_valence(i)
            if v < 1e-3:
                continue  # Already saturated

            if not is_integer_order(v):
                print(f"Non-integer order {v} detected in hydrogenation!")
                v = float(int(v))

            for j in range(to_integer_order(v)):
                added += 1
                new_index = offset + added
                self.graph.add_node(new_index, element="H", valence=1.0)
                self.graph.add_edge(i, new_index, order=1.0)

        return added

    def to_rdkit(self) -> Chem.RWMol:
        """
        Converts to an RDkit molecule
        """
        mol = Chem.RWMol()

        graph_to_export = self.graph.subgraph(i for i in self.graph if self.graph.nodes[i]["element"] != "H")
        if len(graph_to_export) == 0:
            graph_to_export = self.graph

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

    @staticmethod
    def randomly_glue_together(
            molecule_a: 'Molecule',
            molecule_b: 'Molecule',
            bond_generator: BondGenerator) -> Optional['Molecule']:
        """
        Given two molecules, attempt to glue them together by creating a bond
        between free valence points. Returns None if no compatible valence points exist.

        Parameters
        ----------
        molecule_a
            First molecule
        molecule_b
            Second molecule
        Returns
        -------
        None if molecules could not be glued together, otherwise
        the newly-created, glued-together molecule.
        """

        if len(molecule_a.attach_points) == 0:
            return None
        if len(molecule_b.attach_points) == 0:
            return None

        # Create a copy of molecule b that does not share node ids with molecule a
        molecule_b = molecule_b.copy()
        molecule_b.make_disjoint(molecule_a)

        # Choose the pair of free valence points to glue together
        a_point = random.choice(molecule_a.attach_points)
        b_point = random.choice(molecule_b.attach_points)
        a_valence = molecule_a.free_valence(a_point)
        b_valence = molecule_b.free_valence(b_point)

        # The order of the new bond
        new_order = bond_generator.generate_and_check(
            molecule_a.graph.nodes[a_point]["element"],
            molecule_b.graph.nodes[b_point]["element"],
            a_valence, b_valence)

        if new_order < 0:
            return None

        # Build the glued-together molecule
        mol = Molecule()

        assert a_point in molecule_a.graph.nodes
        assert b_point in molecule_b.graph.nodes
        mol.graph = networkx.compose(molecule_b.graph, molecule_a.graph)
        assert a_point in mol.graph.nodes
        assert b_point in mol.graph.nodes

        # Create the bond
        mol.graph.add_edge(a_point, b_point, order=new_order)

        return mol
