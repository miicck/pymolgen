import random
import networkx
import pysmiles
from copy import deepcopy
from typing import Iterator, Dict, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Draw
from contextlib import redirect_stderr
from io import StringIO
import mendeleev


class FractionalOrderException(Exception):
    def __init__(self, message):
        super().__init__(message)


class Molecule:

    ##################
    # Dunder methods #
    ##################

    def __init__(self, allow_frac_order: bool = True):
        """
        Build an empty molecule
        """
        self._graph: networkx.Graph = networkx.Graph()
        self._allow_frac_order: bool = allow_frac_order
        self._valence_lookup: Dict[str, int] = dict()

    def __str__(self):
        return pysmiles.write_smiles(self.graph)

    ##############
    # Properties #
    ##############

    @property
    def graph(self) -> networkx.Graph:
        return self._graph

    @graph.setter
    def graph(self, val: networkx.Graph):
        self._graph = val
        for i in self._graph.nodes:
            self._graph.nodes[i]["stereo"] = None

    @property
    def atom_count(self) -> int:
        return len(self.graph)

    @property
    def total_free_valence(self) -> int:
        return sum(vp[1] for vp in self.free_valence_points)

    @property
    def free_valence_points(self) -> Iterator[Tuple[int, int]]:
        for i in self.graph.nodes:

            # Work out the total bond order
            # eminating from this atom
            total_order = 0
            for j in self.graph[i]:
                total_order += self.graph[i][j]["order"]

            # Work out the valence of the atom
            element = self.graph.nodes[i]["element"]
            if element not in self._valence_lookup:
                self._valence_lookup[element] = mendeleev.element(element).nvalence()
            valence = self._valence_lookup[element]

            # Check if the total order is an integer
            diff = abs(total_order - round(total_order))
            if diff > 10e-4:
                if self._allow_frac_order:
                    continue  # Skip fractional order

                # Error out if we're not allowing fractional orders
                bond_orders = [f"{self.graph.nodes[i]} - {self.graph[i][j]['order']} - {self.graph.nodes[j]}"
                               for j in self.graph[i]]
                raise FractionalOrderException(bond_orders)

            total_order = round(total_order)
            if valence - total_order > 0:
                yield i, valence - total_order

    ###########
    # Methods #
    ###########

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
        if Chem.MolFromSmiles(smiles) is None:
            return None

        io = StringIO()
        with redirect_stderr(io):
            self.graph = pysmiles.read_smiles(smiles, explicit_hydrogen=True)

        return self

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

            nodes.add(random.choice(list(neighbours)))

        frag = Molecule()
        frag.graph = self.graph.subgraph(nodes)

        # Look for broken bonds
        for i in nodes:
            for n in self.graph[i]:
                if n in nodes:
                    continue

                # Bond between i -> n has been broken
                if self.graph.nodes[i]["aromatic"] and self.graph.nodes[n]["aromatic"]:
                    return self.random_fragment(min_size=min_size, max_size=max_size)

        return frag

    @property
    def valid_smiles(self) -> bool:
        return Chem.MolFromSmiles(str(self)) is not None

    def plot(self):
        if not self.valid_smiles:
            raise Exception(f"Tried to plot structure with invalid smiles string: {self}")

        smiles = str(self)
        m = Chem.MolFromSmiles(smiles)
        Draw.ShowMol(m, size=(1024, 1024))

    def make_disjoint(self, other: 'Molecule') -> None:
        """
        Relabel the atoms in this molecule so
        they don't conflict with the atoms in other.

        Parameters
        ----------
        other
            Molecule that we want to make disjoint to
        """
        i_new = max(i for i in other.graph.nodes)
        i_new = max(i_new, max(i for i in self.graph.nodes))
        i_new += 1

        for i in self.graph.nodes:
            if i in other.graph.nodes:
                self.graph = networkx.relabel_nodes(self.graph, {i: i_new})
                i_new += 1

        for i in self.graph.nodes:
            assert i not in other.graph.nodes

    def remove_random_atom(self, element: str = None) -> bool:

        # Get the list of possible nodes to remove (without breaking the molecule apart)
        nodes = list(n for n in self.graph.nodes if len(self.graph[n]) == 1)
        if element is not None:
            nodes = [n for n in nodes if self.graph.nodes[n]["element"].upper() == element.upper()]

        if len(nodes) == 0:
            return False

        to_remove = random.choice(nodes)
        self.graph.remove_node(to_remove)
        return True

    def hydrogenate(self) -> int:
        """
        Terminates all free valence in this molecule with hydrogen
        """
        new_graph = self.graph.copy()

        added = 0
        for i, v in self.free_valence_points:

            for j in range(v):
                added += 1
                new_node = max(new_graph.nodes) + 1
                new_graph.add_node(new_node, element="H")
                new_graph.add_edge(i, new_node, order=1)

        self.graph = new_graph
        assert self.total_free_valence == 0
        return added

    @staticmethod
    def randomly_glue_together(molecule_a: 'Molecule', molecule_b: 'Molecule') -> Optional['Molecule']:
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

        # If molecules can't be glued together, return None
        if molecule_a.total_free_valence == 0:
            return None
        if molecule_b.total_free_valence == 0:
            return None

        # Create a copy of molecule b that does not share node ids with molecule a
        molecule_b = molecule_b.copy()
        molecule_b.make_disjoint(molecule_a)

        # Check that the above worked
        for n in molecule_b.graph.nodes:
            assert n not in molecule_a.graph.nodes

        # Choose the pair of free valence points to glue together
        a_point, a_valence = random.choice(list(molecule_a.free_valence_points))
        b_point, b_valence = random.choice(list(molecule_b.free_valence_points))

        # The order of the new bond
        new_order = random.randint(1, min(a_valence, b_valence))

        # Build the glued-together molecule
        mol = Molecule()
        mol.graph = networkx.compose(molecule_a.graph, molecule_b.graph)

        # Create the bond
        mol.graph.add_edge(a_point, b_point, order=new_order)

        return mol
