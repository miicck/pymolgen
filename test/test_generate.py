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
        smiles_before_h = str(mol)
        mol.hydrogenate()
        smiles_after_h = str(mol)
        if not mol.plot(timeout=0.5):
            raise Exception("Could not plot invalid smiles (before/after hydrogenation):\n"
                            f"{smiles_before_h}\n"
                            f"{smiles_after_h}")
