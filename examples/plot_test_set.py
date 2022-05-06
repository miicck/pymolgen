from pymolgen.generate import SmilesDataset

with open("../datasets/test_set", "r") as test_set:
    smiles = list(line.strip() for line in test_set)

d = SmilesDataset(smiles)

for i, mol in enumerate(d):
    try:
        mol.plot(title=d.get_smiles(i))
    except:
        pass
