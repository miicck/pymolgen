from pymolgen.molecule import Molecule, FractionalOrderException
from pymolgen.molecule_formats import molecule_from_smiles, molecule_from_sdf, molecule_from_sdf_large, molecule_from_sdf_lines
from pymolgen.bond_generator import BondGenerator
from typing import Iterable, Iterator, Callable, List, Dict, Tuple, Union
import random
import glob
from abc import ABC, abstractmethod


class MoleculeDataset(ABC):

    def __init__(self):
        self.mols: Dict[int, Molecule] = dict()

    def __getitem__(self, i: int) -> Molecule:
        # Retrieve existing molecule
        if i in self.mols:
            return self.mols[i]

        # Store and return new molecule
        new_mol = self.load_molecule(i)
        self.mols[i] = new_mol
        return new_mol

    def random_molecule(self) -> Molecule:
        random_molecule_i = random.randint(0, len(self) - 1)
        #print('Random molecule i =', random_molecule_i)
        return self[random_molecule_i]

    @abstractmethod
    def load_molecule(self, i: int) -> Molecule:
        raise NotImplementedError()

    @abstractmethod
    def __len__(self):
        raise NotImplementedError()


class SmilesDataset(MoleculeDataset):

    def __init__(self, smiles: Iterable[str]):
        super().__init__()
        self.smiles: List[str] = list(smiles)

    def load_molecule(self, i: int) -> Molecule:
        return molecule_from_smiles(self.smiles[i])

    def __len__(self):
        return len(self.smiles)


class SDFDataset(MoleculeDataset):

    def __init__(self, sdf_file_abspaths: Iterable[str]):
        super().__init__()
        self.sdf_file_abspaths: List[str] = list(sdf_file_abspaths)

    def load_molecule(self, i: int) -> Molecule:
        return molecule_from_sdf(self.sdf_file_abspaths[i])

    def __len__(self):
        return len(self.sdf_file_abspaths)

class SDFDatasetLarge(MoleculeDataset):

    def __init__(self, sdf_file_abspath: str, max_n = None):
        super().__init__()
        self.sdf_file_abspath: str = sdf_file_abspath
        self.start_lines: List[int] = self.get_start_lines()

        if max_n is not None:
            self.start_lines = self.start_lines[:max_n]

    def load_molecule(self, i: int) -> Molecule:
        start_line = self.start_lines[i]
        return molecule_from_sdf_large(self.sdf_file_abspath, start_line)

    def __len__(self):
        return len(self.start_lines)

    def get_start_lines(self):
        
        start_lines = []
        
        with open(self.sdf_file_abspath) as infile:
            n = 1
            for line in infile:
                if 'V2000' in line:
                    start_lines.append(n-3)
                n += 1

        return start_lines

class SDFDatasetLargeRAM(MoleculeDataset):

    def __init__(self, sdf_file_abspath: str, max_n = None):
        super().__init__()
        self.sdf_file_abspath: str = sdf_file_abspath
        self.start_lines: List[int] = self.get_start_lines()

        if max_n is not None:
            self.start_lines = self.start_lines[:max_n]

        self.lines = []

        max_line = self.start_lines[-1]

        print('DATABASE MOLECULES =', len(self.start_lines))

        with open(self.sdf_file_abspath) as infile:
            n = 0
            for line in infile:
                self.lines.append(line)
                if n == max_line:
                    break

        if max_n is not None and max_n > len(self.lines):
            self.lines = self.lines[:start_lines[max_n+1]]

        print('DATABSE MOLECULES LINES =', len(self.lines))

    def load_molecule(self, i: int) -> Molecule:
        start_line = self.start_lines[i]

        if i < len(self.start_lines) - 1:
            end_line = self.start_lines[i+1]
            mol_lines = self.lines[start_line:end_line]
        else:
            mol_lines = self.lines[start_line:]

        return molecule_from_sdf_lines(mol_lines)

    def __len__(self):
        return len(self.start_lines)

    def get_start_lines(self):
        
        start_lines = []
        
        with open(self.sdf_file_abspath) as infile:
            n = 0
            for line in infile:
                if 'V2000' in line:
                    start_lines.append(n-3)
                n += 1

        return start_lines

def generate_from_molecules(
        dataset: Iterable['Molecule'],
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