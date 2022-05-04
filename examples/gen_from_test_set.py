from pymolgen.molecule import Molecule
from pymolgen.generate import generate_from_fragments

smiles = []
with open("../datasets/test_set", "r") as smiles_chembl:
    for line in smiles_chembl:
        smiles.append(line)


def accept_mol(mol: Molecule):
    return mol.atom_count > 40


for mol in generate_from_fragments(smiles, accept_mol):
    try:
        mol.plot(timeout=0.5)
    except:
        pass
