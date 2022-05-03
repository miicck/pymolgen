import csv


def convert_data():
    # Convert ChEMBL dataset
    with open("datasets/ChEMBL.csv", "r") as f:
        with open("datasets/smiles_chembl", "w") as f2:
            reader = csv.reader(f, delimiter=';')
            for i, row in enumerate(reader):
                if i == 0:
                    i_smiles = row.index("Smiles")
                else:
                    smiles = row[i_smiles].strip()
                    if len(smiles) > 0:
                        f2.write(row[i_smiles] + "\n")


convert_data()
