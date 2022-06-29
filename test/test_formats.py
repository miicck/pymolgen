import sys,os

from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *
from pymolgen.generate import *

def test_bru():
	mol = molecule_from_bru("atom H\natom H\nbond 0 1 1")
	assert mol.atom_count == 2
	assert len(mol.graph.edges) == 1
	assert mol.graph.nodes[0]["valence"] == 1.0
	assert mol.graph.nodes[1]["valence"] == 1.0
	assert 0 in mol.graph[1]
	assert 1 in mol.graph[0]

def test_molecule_from_sdf():
	mol = molecule_from_sdf('mol-1.sdf')
	#plot_molecule(mol)
	#plot_molecule_graph(mol)

	smi = molecule_to_smiles(mol)
	print(smi)
	
	assert smi == 'CC1=CC(NC(=O)CS(=O)(=O)C2=CN(CC3=CC(C4=CNN=C4)=CC=C3)C3=C2C(F)=CC(F)=C3)=NO1'


def test_sdfdatasetlarge():
	dataset = SDFDatasetLarge('mol-1.sdf')

	mol = dataset.load_molecule(0)

	smi = molecule_to_smiles(mol)

	print(smi)

	assert smi == 'CC1=CC(NC(=O)CS(=O)(=O)C2=CN(CC3=CC(C4=CNN=C4)=CC=C3)C3=C2C(F)=CC(F)=C3)=NO1'

def test_sdfdatasetlarge2():

	dataset = SDFDatasetLarge('../datasets/database1000/database1000.sdf')

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

def test_sdf_to_asf():
	sdf_to_asf('../datasets/database1000/database1000.sdf', '../datasets/database1000/database1000.asf')
	read_asf_file('../datasets/database1000/database1000.asf')	

def test_molecule_to_atoms_bonds():

	mol = molecule_from_sdf('../datasets/sdf/mol-1.sdf')
	atoms, bonds, valences = molecule_to_atoms_bonds(mol)

	check = ['C', 'C', 'C', 'C', 'N', 'C', 'O', 'C', 'S', 'O', 'O', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'N', 'C', 'C', 'C', 'C', 'C', 'F', 'C', 'C', 'F', 'C', 'N', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

	assert atoms == check

	check = [4, 4, 4, 4, 3, 4, 2, 4, 6, 2, 2, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 1, 4, 4, 1, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

	assert valences == check

def test_molecule_to_sdf():

	mol = molecule_from_sdf('../datasets/sdf/mol-1.sdf')

	lines = molecule_to_sdf(mol)

	check = ['Molecule\n pymolgen\n\n', ' 55 59  0  0  1  0  0  0  0  0999 V2000\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n', '  1  2  1  0  0  0  0\n', '  1 37  1  0  0  0  0\n', '  1 38  1  0  0  0  0\n', '  1 39  1  0  0  0  0\n', '  2 36  1  0  0  0  0\n', '  2  3  2  0  0  0  0\n', '  3  4  1  0  0  0  0\n', '  3 40  1  0  0  0  0\n', '  4  5  1  0  0  0  0\n', '  4 35  2  0  0  0  0\n', '  5  6  1  0  0  0  0\n', '  5 41  1  0  0  0  0\n', '  6  7  2  0  0  0  0\n', '  6  8  1  0  0  0  0\n', '  8  9  1  0  0  0  0\n', '  8 42  1  0  0  0  0\n', '  8 43  1  0  0  0  0\n', '  9 10  2  0  0  0  0\n', '  9 11  2  0  0  0  0\n', '  9 12  1  0  0  0  0\n', ' 12 34  1  0  0  0  0\n', ' 12 13  2  0  0  0  0\n', ' 13 14  1  0  0  0  0\n', ' 13 44  1  0  0  0  0\n', ' 14 15  1  0  0  0  0\n', ' 14 27  1  0  0  0  0\n', ' 15 16  1  0  0  0  0\n', ' 15 45  1  0  0  0  0\n', ' 15 46  1  0  0  0  0\n', ' 16 26  2  0  0  0  0\n', ' 16 17  1  0  0  0  0\n', ' 17 18  2  0  0  0  0\n', ' 17 47  1  0  0  0  0\n', ' 18 19  1  0  0  0  0\n', ' 18 48  1  0  0  0  0\n', ' 19 20  2  0  0  0  0\n', ' 19 49  1  0  0  0  0\n', ' 20 21  1  0  0  0  0\n', ' 20 26  1  0  0  0  0\n', ' 21 25  2  0  0  0  0\n', ' 21 22  1  0  0  0  0\n', ' 22 23  2  0  0  0  0\n', ' 22 50  1  0  0  0  0\n', ' 23 24  1  0  0  0  0\n', ' 24 25  1  0  0  0  0\n', ' 24 51  1  0  0  0  0\n', ' 25 52  1  0  0  0  0\n', ' 26 53  1  0  0  0  0\n', ' 27 34  2  0  0  0  0\n', ' 27 28  1  0  0  0  0\n', ' 28 29  2  0  0  0  0\n', ' 28 54  1  0  0  0  0\n', ' 29 30  1  0  0  0  0\n', ' 29 31  1  0  0  0  0\n', ' 31 32  2  0  0  0  0\n', ' 31 55  1  0  0  0  0\n', ' 32 33  1  0  0  0  0\n', ' 32 34  1  0  0  0  0\n', ' 35 36  1  0  0  0  0\n', 'M  END\n', '> <valences>\n', '4 4 4 4 3 4 2 4 6 2 2 4 4 3 4 4 4 4 4 4 4 4 3 3 4 4 4 4 4 1 4 4 1 4 3 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n\n']

	assert lines == check

def test_read_from_sdf_valences():

	fragment_database = SDFDatasetLargeRAM('../datasets/fragments/fragments.sdf')

	elements = []
	valences = []

	for i in range(len(fragment_database)):
		mol = fragment_database.load_molecule(i)

		for i in mol.graph.nodes:
			elements.append(mol.graph.nodes[i]['element'])

		for i in mol.graph.nodes:
			valences.append(mol.graph.nodes[i]['valence'])

	check = ['C', 'H', 'H', 'H', 'C', 'C', 'C', 'N', 'O', 'H', 'H', 'N', 'C', 'O', 'H', 'H', 'C', 'S', 'O', 'O', 'C', 'C', 'C', 'N', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'C', 'N', 'N', 'C', 'F', 'H', 'O', 'C', 'C', 'N', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'H', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'N', 'C', 'N', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'S', 'H', 'S', 'C', 'H', 'C', 'H', 'H', 'H', 'N', 'C', 'O', 'C', 'N', 'C', 'O', 'C', 'H', 'N', 'C', 'O', 'C', 'N', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'N', 'C', 'C', 'C', 'C', 'C', 'H', 'N', 'C', 'H', 'H', 'N', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'C', 'N', 'O', 'C', 'C', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'H', 'H', 'C', 'N', 'C', 'C', 'H', 'N', 'N', 'C', 'H', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'N', 'C', 'O', 'C', 'C', 'C', 'H', 'Cl', 'H', 'H', 'H', 'H', 'H', 'N', 'C', 'C', 'C', 'N', 'N', 'C', 'C', 'C', 'H', 'H', 'O', 'C', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'H', 'H', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'H', 'C', 'N', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'N', 'C', 'C', 'H', 'H', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'N', 'C', 'H', 'H', 'C', 'H', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'N', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'C', 'N', 'C', 'H', 'H', 'H', 'C', 'C', 'C', 'C']

	assert elements == check

	check = [4, 1, 1, 1, 4, 4, 4, 3, 2, 1, 1, 3, 4, 2, 1, 1, 4, 6, 2, 2, 4, 4, 4, 3, 1, 1, 1, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 1, 1, 1, 4, 1, 1, 1, 4, 4, 3, 3, 4, 1, 1, 2, 4, 4, 3, 1, 4, 4, 4, 4, 4, 4, 4, 2, 1, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 3, 4, 3, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 1, 2, 4, 1, 4, 1, 1, 1, 3, 4, 2, 4, 3, 4, 2, 4, 1, 3, 4, 2, 4, 3, 4, 2, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 4, 4, 4, 4, 1, 3, 4, 1, 1, 3, 4, 4, 4, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 4, 3, 2, 4, 4, 1, 4, 4, 4, 4, 4, 4, 3, 1, 1, 4, 3, 4, 4, 1, 3, 3, 4, 1, 4, 1, 1, 1, 1, 1, 1, 3, 4, 2, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 3, 4, 4, 4, 3, 3, 4, 4, 4, 1, 1, 2, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 1, 1, 4, 4, 4, 4, 1, 4, 4, 1, 1, 4, 4, 4, 4, 4, 4, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 2, 4, 4, 4, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 4, 3, 4, 4, 4, 1, 1, 1, 1, 3, 4, 4, 1, 1, 4, 4, 2, 1, 1, 1, 1, 1, 3, 4, 1, 1, 4, 1, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 4, 4, 3, 4, 4, 1, 4, 4, 4, 1, 4, 3, 4, 1, 1, 1, 4, 4, 4, 4]

	assert valences == check

	print(elements)
	print(valences)

