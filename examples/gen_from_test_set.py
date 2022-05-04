from pymolgen.molecule import Molecule
from pymolgen.generate import generate_from_fragments

with open("../datasets/test_set", "r") as test_set:
    smiles = list(line for line in test_set)

for mol in generate_from_fragments(smiles, lambda m: m.atom_count > 50, max_frag_size=20):
    try:
        mol.plot()
    except:
        pass
