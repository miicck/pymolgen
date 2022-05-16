import molecule
import generate
import bond_generator

min_frag_size = 1
max_frag_size = 50

mol = molecule.Molecule()

mol.load_bru('parent2.bru')

mol.remove_atom(30)

print(mol.attach_points)

smiles = []
with open("/home/pczbf/pymolgen/datasets/test_set", "r") as smiles_chembl:
    for line in smiles_chembl:
        smiles.append(line)

dataset = generate.SmilesDataset(smiles)

#generate a random fragment
frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

frag.plot()

print(frag.attach_points)

# Add generated random fragment to mol
mol = molecule.Molecule.randomly_glue_together(mol, frag, bond_generator.RandomBondGenerator())

mol.plot()
