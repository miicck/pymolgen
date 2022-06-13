import sys,os

from pymolgen.generate import SDFDatasetLarge
from pymolgen.molecule_formats import *

dataset = SDFDatasetLarge('../datasets/database1000/database1000.sdf')

#print(dataset.start_lines[:1000])

smiles = []

outfile = open('test_database.sdf', 'w')

for i in range(35):

	mol = dataset.load_molecule(i)

	smi = molecule_to_smiles(mol)

	smiles.append(smi)

	lines = molecule_to_sdf(mol)

	for line in lines:
		outfile.write(line)
	outfile.write('$$$$\n')

outfile.close()

os.system('obabel -xn -isdf test_database.sdf -ocan -O test_database.can -b')

with open('test_database.can') as infile1, open('../datasets/database1000/database1000.can') as infile2:

	n = 0
	for line1, line2 in zip(infile1, infile2):
		print(n)
		print(line1)
		print(line2)
		assert line1.strip() == line2.strip()
		n += 1

