from abc import ABC, abstractmethod
import random


class BondGenerator(ABC):

    @abstractmethod
    def generate(self, atom_from: str, atom_to: str, free_valence_from: float, free_valence_to: float) -> float:
        """
        Given a pair of atoms with the specified free valence on each,
        gives a bond order to generate between them.

        Parameters
        ----------
        atom_from
            Atom at start of new bond
        atom_to
            Atom at end of new bond
        free_valence_from
            Valence available at the start atom
        free_valence_to
            Valence available at the end atom

        Returns
        -------
        The order of the new bond to generate (or a negative number if generation was impossible).
        """
        raise NotImplementedError

    def generate_and_check(self, atom_from: str, atom_to: str,
                           free_valence_from: float, free_valence_to: float) -> float:
        max_order = BondGenerator.max_possible_bond_order(free_valence_from, free_valence_to)
        new_order = self.generate(atom_from, atom_to, free_valence_from, free_valence_to)
        assert new_order <= max_order
        return new_order

    @staticmethod
    def max_possible_bond_order(free_valence_from: float, free_valence_to: float) -> int:
        max_order = min(free_valence_from, free_valence_to)
        diff = abs(max_order - round(max_order))
        if diff > 10e-3:
            raise FractionalOrderException
        return round(max_order)


class RandomBondGenerator(BondGenerator):

    def generate(self, atom_from: str, atom_to: str, free_valence_from: float, free_valence_to: float) -> float:
        return random.randint(1, BondGenerator.max_possible_bond_order(free_valence_from, free_valence_to))


class MaxOrderBondGenerator(BondGenerator):

    def generate(self, atom_from: str, atom_to: str, free_valence_from: float, free_valence_to: float) -> float:
        return BondGenerator.max_possible_bond_order(free_valence_from, free_valence_to)
