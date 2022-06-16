import sys,os

from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *
from pymolgen.generate import SDFDatasetLarge

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
