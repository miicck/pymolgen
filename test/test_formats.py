from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *


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


