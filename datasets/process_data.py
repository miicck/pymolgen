import csv
import random


def convert_data():
    # Convert ChEMBL dataset
    with open("ChEMBL.csv", "r") as f:
        with open("smiles_chembl", "w") as f2:
            reader = csv.reader(f, delimiter=';')
            for i, row in enumerate(reader):
                if i == 0:
                    i_smiles = row.index("Smiles")
                else:
                    smiles = row[i_smiles].strip()
                    if len(smiles) > 0:
                        f2.write(row[i_smiles] + "\n")

    all_smiles = []
    with open("smiles_chembl", "r") as f:
        for line in f:
            all_smiles.append(line)

    subset = random.choices(all_smiles, k=1000)
    with open("test_set", "w") as f:
        for s in subset:
            f.write(s)


convert_data()
