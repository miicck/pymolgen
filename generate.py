from pymolgen.molecule import Molecule, FractionalOrderException
from pymolgen.bond_generator import BondGenerator
from typing import Iterable, Iterator, Callable, List, Dict
import random


class SmilesDataset:

    def __init__(self, smiles: Iterable[str]):
        self.smiles: List[str] = list(smiles)
        self.mols: Dict[int, Molecule] = dict()

    def random_molecule(self) -> Molecule:
        i = random.randint(0, len(self.smiles) - 1)

        # Retrieve existing molecule
        if i in self.mols:
            return self.mols[i]

        new_mol = Molecule().load_smiles(self.smiles[i])
        if new_mol is None:
            # Loading molecule failed, try again
            return self.random_molecule()

        # Store and return new molecule
        self.mols[i] = new_mol
        return new_mol


def generate_from_fragments(
        source_smiles: Iterable[str],
        accept: Callable[[Molecule], bool],
        bond_generator: BondGenerator,
        min_frag_size: int = 1,
        max_frag_size: int = None) -> Iterator[Molecule]:
    """
    Given a set of source smiles strings, generate molecules
    by gluing together fragments from the set.

    Parameters
    ----------
    source_smiles
        Smiles strings for source molecules
    accept
        Function to determine if an intermediate
        fragment structure is acceptable as a generated molecule
    bond_generator
        Object used to generate new bonds to glue fragments together
    min_frag_size
        The minimum size of fragments to generate
    max_frag_size
        The maximum size of fragments to generate (None => any size)

    Returns
    -------
    An iterator of generated molecules. Will iterate forever.
    """

    dataset = SmilesDataset(source_smiles)

    source_mols = []
    saved_mols = dict()

    while True:

        mol = None

        while True:

            # Get a random fragment from the dataset
            frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

            if mol is None:
                mol = frag
                continue

            # Yield acceptable molecule
            if accept(mol):
                yield mol
                break

            # Add a random fragment
            mol = Molecule.randomly_glue_together(mol, frag, bond_generator)
