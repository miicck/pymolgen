from pymolgen.molecule import Molecule, FractionalOrderException
from typing import Iterable, Iterator, Callable
import random


def generate_from_fragments(
        source_smiles: Iterable[str],
        accept: Callable[[Molecule], bool],
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
    min_frag_size
        The minimum size of fragments to generate
    max_frag_size
        The maximum size of fragments to generate (None => any size)

    Returns
    -------
    An iterator of generated molecules. Will iterate forever.
    """
    source_mols = []
    saved_mols = dict()

    while True:

        mol = None

        while True:

            # Get a random fragment to add
            rand_smiles = random.choice(source_smiles)
            if rand_smiles in saved_mols:
                # Recover molecule from already-parsed smiles string
                mol_to_frag = saved_mols[rand_smiles]
            else:
                # Attempt to load a new molecule from a random smiles string
                try:
                    mol_to_frag = Molecule().load_smiles(rand_smiles)
                    if mol_to_frag is None:
                        continue
                except:
                    continue

                saved_mols[rand_smiles] = mol_to_frag

            frag = mol_to_frag.random_fragment(min_size=min_frag_size, max_size=max_frag_size)

            if mol is None:
                mol = frag
                continue

            # We've somehow generated an invalid molecule
            if not mol.valid_smiles:
                break

            # Yield acceptable molecule
            if accept(mol):
                yield mol

            try:
                # Add a random fragment
                mol = Molecule.randomly_glue_together(mol, frag)
            except FractionalOrderException:
                print(mol)
                print(frag)
                raise FractionalOrderException

            if mol == None:
                break  # Molecule has no free valence
