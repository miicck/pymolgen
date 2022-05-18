import random
import time
import networkx
from copy import deepcopy
from typing import Iterator, Dict, Optional, Tuple, Set, List
from enum import Enum
from pymolgen.bond_generator import BondGenerator
import numpy as np

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
        n = 0
        for i in val.nodes:
            n += 1
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

    def copy(self) -> 'Molecule':
        """
        Returns
        -------
        A deep copy of this molecule
        """
        return deepcopy(self)

    def set_valence_from_bonds(self):
        """
        Sets the valence of each atom according
        to the total bond order attached to that atom.
        """
        for n in self.graph.nodes:
            self.graph.nodes[n]["valence"] = self.total_order_of_bonds(n)

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

    def molecular_weight(self) -> float:
        d = mass_dict()
        mw = 0
        for i in self.graph.nodes:
            element = self.graph.nodes[i]["element"]
            mw += d[element]
        return mw

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
        return frag

    def is_cyclic(self, i: int) -> bool:
        """
        Returns true if the node at the given index is cyclic.
        Parameters
        ----------
        i
            The index of the node

        Returns
        -------
        True if the node is part of a cycle
        """
        for basis in networkx.cycle_basis(self.graph):
            if i in basis:
                return True
        return False

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

        while True:

            target_size = random.randint(min_size, max_size or len(self.graph))

            # Pick a random starting location
            i = random.choice(list(self.graph.nodes))
            nodes = {i}

            # Add a random neighbouring atom until
            # the fragment is of the desired size
            while len(nodes) < target_size:

                neighbours = set()
                for n in nodes:
                    neighbours.update(self.graph[n])
                neighbours -= nodes

                if len(neighbours) == 0:
                    break

                nodes.add(random.choice(list(neighbours)))

            frag = self.get_fragment(nodes)
            if frag is not None:
                return frag

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
            nodes = [n for n in nodes if self.graph.nodes[n]["element"].capitalize() == element.capitalize()]
        if len(nodes) == 0:
            return False
        to_remove = random.choice(nodes)
        keep_nodes = set(self.graph.nodes) - {to_remove}
        return self.get_fragment(keep_nodes)

    def remove_atom(self, index: int) -> 'Molecule':
        """
        Removes the atom with the given index, and returns the new molecule.
        """
        return self.get_fragment(set(self.graph.nodes) - {index})

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
                print(f"Non-integer free valence {v} detected in hydrogenation")
                v = float(int(v))

            for j in range(to_integer_order(v)):
                added += 1
                new_index = offset + added
                self.graph.add_node(new_index, element="H", valence=1.0)
                self.graph.add_edge(i, new_index, order=1.0)

        return added

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

    @staticmethod
    def glue_together_attachmentpoint(
            molecule_a: 'Molecule',
            molecule_b: 'Molecule',
            bond_generator: BondGenerator,
            attachment_point: 'int') -> Optional['Molecule']:
        """
        Given two molecules, attempt to glue them together by creating a bond
        between the attachment point of molecule_a and free valence points of molecule_b. 
        Returns None if no compatible valence points exist.

        Parameters
        ----------
        molecule_a
            First molecule
        molecule_b
            Second molecule
        attachment_point
            Attachment point of molecule_a
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
        a_point = attachment_point
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

def mass_dict():
    d = {
    'H':1.00797,
    'He':4.00260,
    'Li':6.941,
    'Be':9.01218,
    'B':10.81,
    'C':12.011,
    'N':14.0067,
    'O':15.9994,
    'F':18.998403,
    'Ne':20.179,
    'Na':22.98977,
    'Mg':24.305,
    'Al':26.98154,
    'Si':28.0855,
    'P':30.97376,
    'S':32.06,
    'Cl':35.453,
    'K':39.0983,
    'Ar':39.948,
    'Ca':40.08,
    'Sc':44.9559,
    'Ti':47.90,
    'V':50.9415,
    'Cr':51.996,
    'Mn':54.9380,
    'Fe':55.847,
    'Ni':58.70,
    'Co':58.9332,
    'Cu':63.546,
    'Zn':65.38,
    'Ga':69.72,
    'Ge':72.59,
    'As':74.9216,
    'Se':78.96,
    'Br':79.904,
    'Kr':83.80,
    'Rb':85.4678,
    'Sr':87.62,
    'Y':88.9059,
    'Zr':91.22,
    'Nb':92.9064,
    'Mo':95.94,
    'Ru':101.07,
    'Rh':102.9055,
    'Pd':106.4,
    'Ag':107.868,
    'Cd':112.41,
    'In':114.82,
    'Sn':118.69,
    'Sb':121.75,
    'I':126.9045,
    'Te':127.60,
    'Xe':131.30,
    'Cs':132.9054,
    'Ba':137.33,
    'La':138.9055,
    'Ce':140.12,
    'Pr':140.9077,
    'Nd':144.24,
    'Sm':150.4,
    'Eu':151.96,
    'Gd':157.25,
    'Tb':158.9254,
    'Dy':162.50,
    'Ho':164.9304,
    'Er':167.26,
    'Tm':168.9342,
    'Yb':173.04,
    'Lu':174.967,
    'Hf':178.49,
    'Ta':180.9479,
    'W':183.85,
    'Re':186.207,
    'Os':190.2,
    'Ir':192.22,
    'Pt':195.09,
    'Au':196.9665,
    'Hg':200.59,
    'Tl':204.37,
    'Pb':207.2,
    'Bi':208.9804,
    'Ra':226.0254,
    'Ac':227.0278,
    'Pa':231.0359,
    'Th':232.0381,
    'Np':237.0482,
    'U':238.029
    }
    return d
