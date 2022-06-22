import sys,os

from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *
from pymolgen.molecule import *
from pymolgen.fragment_mol import *
from rdkit import Chem

dir = os.path.dirname(os.path.realpath(__file__))

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


def test_get_fragments():

	mol = molecule_from_sdf('mol-1.sdf')

	single_bonds = mol.get_single_bonds_not_h_not_c()

	new = split_mol(mol, single_bonds)

	mol2 = molecule_from_sdf('mol-1-can.sdf')

	bonds2 = mol2.get_single_bonds_not_h_not_c()

	new2 = split_mol(mol2, bonds2)

	n_equal = 0

	for i in new:
		for j in new2:
			
			if get_atom_list(i) == get_atom_list(j):

				n_equal += 1

				equal = networkx.is_isomorphic(i,j)

				assert equal is True

	assert n_equal == 16

def test_get_fragments_mol1():

	mol = molecule_from_sdf('mol-1.sdf')

	single_bonds = mol.get_single_bonds_not_h_not_c()

	new = split_mol(mol, single_bonds)

	new_nodes_view = []

	for fragment in new:
		nodes = []
		for i in fragment.nodes:
			nodes.append(i)
		new_nodes_view.append(nodes)

	print(new_nodes_view)


	print_fragments(new)

	saved_fragments = [[0, 36, 37, 38],[1, 2, 3, 34, 35, 39],[40, 4],[5, 6],[41, 42, 7],[8, 9, 10],[33, 11, 12, 13, 43, 53, 54, 26, 27, 28, 30, 31],[44, 45, 14],[46, 15, 16, 17, 18, 19, 52, 47, 48, 25],[49, 50, 51, 20, 21, 22, 23, 24],[29],[32]]

	for i in range(len(saved_fragments)):
		assert new_nodes_view[i] == saved_fragments[i]


def test_update_freq():

	frequencies = {}

	frequencies[(0,1,0,1)] = 1
	frequencies[(1,2,0,1)] = 1
	frequencies[(2,3,0,1)] = 1
	frequencies[(3,4,0,1)] = 1

	update_freq(frequencies, frag1_index = 0, frag2_index = 1, frag1_map = {10:0}, frag2_map = {20:1}, frag1_bond = 10, frag2_bond = 20)

	assert frequencies == {(0, 1, 0, 1): 2, (1, 2, 0, 1): 1, (2, 3, 0, 1): 1, (3, 4, 0, 1): 1}

	update_freq(frequencies, frag1_index = 4, frag2_index = 5, frag1_map = {10:0}, frag2_map = {20:1}, frag1_bond = 10, frag2_bond = 20)

	assert frequencies == {(0, 1, 0, 1): 2, (1, 2, 0, 1): 1, (2, 3, 0, 1): 1, (3, 4, 0, 1): 1, (4, 5, 0, 1): 1}

	update_freq(frequencies, frag1_index = 4, frag2_index = 5, frag1_map = {0:0}, frag2_map = {1:1}, frag1_bond = 0, frag2_bond = 1)

	assert frequencies == {(0, 1, 0, 1): 2, (1, 2, 0, 1): 1, (2, 3, 0, 1): 1, (3, 4, 0, 1): 1, (4, 5, 0, 1): 2}

	update_freq(frequencies, frag1_index = 4, frag2_index = 5, frag1_map = {2:2}, frag2_map = {3:3}, frag1_bond = 2, frag2_bond = 3)

	assert frequencies == {(0, 1, 0, 1): 2, (1, 2, 0, 1): 1, (2, 3, 0, 1): 1, (3, 4, 0, 1): 1, (4, 5, 0, 1): 2, (4, 5, 2, 3): 1}

def test_get_fragment_index():
	fragment_database = []
	frequencies = {}

	mol = molecule_from_sdf('mol-1.sdf')

	fragments, pairs, bonds = get_fragments_dataset(mol)

	for i in range(len(pairs)):
		update_database(pairs[i], bonds[i], fragment_database, fragments, frequencies)

	out = print_fragments(fragment_database)
	print(out)
	print(frequencies)

	mol = molecule_from_sdf('mol-1-can.sdf')

	for i in range(len(pairs)):
		update_database(pairs[i], bonds[i], fragment_database, fragments, frequencies)

	out = print_fragments(fragment_database)
	print(out)
	print(frequencies)

	check = {(0, 1, 0, 1): 2, (1, 2, 3, 4): 2, (2, 3, 4, 5): 2, (3, 4, 5, 7): 2, (4, 5, 7, 8): 2, (5, 6, 8, 11): 2, (4, 6, 7, 13): 2, (4, 7, 7, 15): 2, (7, 8, 19, 20): 2, (6, 9, 28, 29): 2, (6, 9, 31, 29): 2}

	assert frequencies ==  check

def test_make_fragment_database():

	make_fragment_database('../datasets/database1000/database10.sdf', 'outputs/fragments_sdf', 'outputs/fragments_txt', 'outputs/frequencies_txt')

test_make_fragment_database()