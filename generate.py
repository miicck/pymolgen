from pymolgen.molecule import Molecule, FractionalOrderException
from pymolgen.bond_generator import BondGenerator
from typing import Iterable, Iterator, Callable, List, Dict, Tuple, Union
import random


class SmilesDataset:

    def __init__(self, smiles: Iterable[str]):
        self.smiles: List[str] = list(smiles)
        self.mols: Dict[int, Molecule] = dict()

    def __getitem__(self, i: int) -> Molecule:
        # Retrieve existing molecule
        if i in self.mols:
            return self.mols[i]

        # Store and return new molecule
        new_mol = Molecule().load_smiles(self.smiles[i])
        self.mols[i] = new_mol
        return new_mol

    def __len__(self):
        return len(self.smiles)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def get_smiles(self, i: int) -> str:
        return self.smiles[i]

    def random_molecule(self) -> Union[Molecule]:
        return self[random.randint(0, len(self.smiles) - 1)]


def generate_from_fragments(
        source_smiles: Iterable[str],
        accept: Callable[[Molecule], bool],
        bond_generator: BondGenerator,
        min_frag_size: int = 1,
        max_frag_size: int = None) -> Iterator[Tuple[Molecule, int]]:
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

    while True:

        mol = None
        fragments = 0

        while True:

            # Get a random fragment from the dataset
            frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

            if mol is None:
                mol = frag
                fragments = 1
                continue

            # Yield acceptable molecule
            if accept(mol):
                yield mol, fragments
                break

            # Add a random fragment
            mol = Molecule.randomly_glue_together(mol, frag, bond_generator)
            fragments += 1
