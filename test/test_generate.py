from pymolgen.molecule import Molecule
from pymolgen.generate import generate_from_fragments
from pymolgen.bond_generator import RandomBondGenerator


def test_gen():
    smiles = []
    with open("../datasets/test_set", "r") as smiles_chembl:
        for line in smiles_chembl:
            smiles.append(line)

    def accept_mol(mol: Molecule):
        return mol.atom_count > 40

    for i, mol in enumerate(generate_from_fragments(smiles, accept_mol, RandomBondGenerator())):
        assert mol is not None
        if i > 100:
            break
