from pymolgen.generate import SmilesDataset

with open("../datasets/test_set", "r") as test_set:
    smiles = list(line.strip() for line in test_set)

d = SmilesDataset(smiles)

for i, mol in enumerate(d):
    try:
        print(d.get_smiles(i))
        mol.plot(title=d.get_smiles(i))
    except Exception as e:
        print("Failed: ", d.get_smiles(i), e)
        pass
