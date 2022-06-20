import sys,os

from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *
from pymolgen.molecule import *

def test_is_unsaturated():
	mol = molecule_from_sdf('mol-1.sdf')

	smi = molecule_to_smiles(mol)

	unsaturated_list = []

	for i in mol.graph.nodes:
		unsaturated_list.append(mol.is_unsaturated(i))

	check = [False, True, True, True, False, True, True, False, True, True, True, True, True, False, False, True, True, True, True, True, True, True, True, False, True, True, True, True, True, False, True, True, False, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False]

	for i in range(len(unsaturated_list)):
		assert check[i] == unsaturated_list[i]

def test_is_hydrogen():
	mol = molecule_from_sdf('mol-1.sdf')

	smi = molecule_to_smiles(mol)

	hydrogen_list = []

	for i in mol.graph.nodes:
		hydrogen_list.append(mol.is_hydrogen(i))

	check = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True]

	for i in range(len(hydrogen_list)):
		assert check[i] == hydrogen_list[i]

def test_get_hydrogen_neighbours():
	mol = molecule_from_sdf('mol-1.sdf')

	hydrogen_list = mol.get_hydrogen_neighbours(7)

	all_hydrogens_list = []

	for i in mol.graph.nodes:
		hydrogen_list = mol.get_hydrogen_neighbours(i)

		all_hydrogens_list.append(hydrogen_list)

	check = [[36, 37, 38], [], [39], [], [40], [], [], [41, 42], [], [], [], [], [43], [], [44, 45], [], [46], [47], [48], [], [], [49], [], [50], [51], [52], [], [53], [], [], [54], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]

	for i in range(len(all_hydrogens_list)):
		assert check[i] == all_hydrogens_list[i]

def test_get_noncyclic_unsaturated_fragments():

	mol = molecule_from_sdf('mol-1.sdf')

	fragments = mol.get_noncyclic_unsaturated_fragments()

	print(fragments)

test_get_noncyclic_unsaturated_fragments()