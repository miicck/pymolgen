import random

from pymolgen.molecule import Molecule
from pymolgen.generate import generate_from_fragments


def test_gen_chembl():
    smiles = []
    with open("../datasets/smiles_chembl", "r") as smiles_chembl:
        for line in smiles_chembl:
            smiles.append(line)

    def accept_mol(mol: Molecule):
        return mol.atom_count > 40

    for mol in generate_from_fragments(smiles, accept_mol):
        mol.plot()
