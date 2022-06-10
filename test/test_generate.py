import random
from pymolgen.molecule import Molecule
from pymolgen.generate import *
from pymolgen.bond_generator import RandomBondGenerator

# Ensure same generation each time
random.seed(1024)


def test_gen():
    smiles = []
    with open("../datasets/test_set", "r") as smiles_chembl:
        for line in smiles_chembl:
            smiles.append(line)

    dataset = SmilesDataset(smiles)

    def accept_mol(mol: Molecule):
        return mol.atom_count > 40

    for i, (mol, fragments) in enumerate(generate_from_molecules(dataset, accept_mol, RandomBondGenerator())):
        assert mol is not None
        if i > 100:
            break


if __name__ == "__main__":
    test_gen()
