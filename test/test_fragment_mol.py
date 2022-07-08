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


def test_make_fragment_database():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database10.sdf', 'outputs/fragments.sdf', 'outputs/fragments.txt', 'outputs/frequencies.txt', 'outputs/frag_frequencies.txt')

def test_mol_bond_frequencies():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database10.sdf', max_n=1)

	print(print_fragments(fragment_database))
	print(frequencies)
	print(frag_frequencies)

	assert frag_frequencies == [1, 1, 1, 1, 2, 1, 1, 1, 1, 2]

def test_mol_bond_frequencies2():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database10.sdf', max_n=2)

	assert frag_frequencies == [1, 1, 1, 2, 2, 1, 1, 1, 1, 3, 1, 1, 1, 1]

	check = """[0, 36, 37, 38] ['C', 'H', 'H', 'H']
[1, 2, 3, 34, 35, 39] ['C', 'C', 'C', 'N', 'O', 'H']
[40, 4] ['H', 'N']
[5, 6] ['C', 'O']
[41, 42, 7] ['H', 'H', 'C']
[8, 9, 10] ['S', 'O', 'O']
[33, 11, 12, 13, 43, 53, 54, 26, 27, 28, 30, 31] ['C', 'C', 'C', 'N', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C']
[46, 15, 16, 17, 18, 19, 52, 47, 48, 25] ['H', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'C']
[49, 50, 51, 20, 21, 22, 23, 24] ['H', 'H', 'H', 'C', 'C', 'N', 'N', 'C']
[29] ['F']
[24, 2] ['H', 'O']
[3, 4, 5, 36, 9, 10, 11, 18, 20, 21, 22, 23, 25, 31] ['C', 'C', 'N', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'H', 'H']
[6, 7, 8, 26, 27, 28, 29, 30] ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H']
[32, 33, 34, 35, 12, 13, 14, 15, 16, 17] ['H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'N', 'C']
"""
	assert check == print_fragments(fragment_database)

def test_canonical_mapping():

	mol = molecule_from_sdf('mol-1.sdf')

	fragments, pairs, bonds = get_fragments_dataset(mol)

	all_canonical_mappings = []

	for fragment in fragments:
		all_canonical_mappings.append(get_canonical_mapping(fragment))

	print(print_fragments(fragments))
	print(all_canonical_mappings)

	assert all_canonical_mappings == [{0: 0, 36: 36, 37: 36, 38: 36}, {1: 1, 2: 2, 3: 3, 34: 34, 35: 35, 39: 39}, {40: 40, 4: 4}, {5: 5, 6: 6}, {41: 41, 7: 7, 42: 41}, {8: 8, 9: 9, 10: 9}, {33: 33, 11: 11, 12: 12, 13: 13, 43: 43, 26: 26, 27: 27, 53: 53, 28: 28, 30: 30, 54: 54, 31: 31}, {44: 44, 14: 14, 45: 44}, {46: 46, 16: 16, 15: 15, 17: 17, 18: 16, 19: 15, 47: 47, 48: 46, 25: 25, 52: 52}, {49: 49, 21: 21, 20: 20, 22: 22, 23: 23, 50: 50, 24: 24, 51: 51}, {29: 29}, {32: 32}]

def test_compound_dict():

	d1 = {0:1, 2:3, 4:5}

	d2 = {1:7, 3:9, 5:11}

	d1 = compound_dict(d1, d2)

	assert d1 == {0: 7, 2: 9, 4: 11}

def test_mol_bond_frequencies3():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database10.sdf', max_n=10)

	check = {(0, 1, 0, 1): 1, (1, 2, 3, 4): 1, (2, 3, 4, 5): 6, (3, 4, 5, 7): 4, (4, 5, 7, 8): 1, (5, 6, 8, 11): 1, (4, 6, 7, 13): 1, (4, 7, 7, 15): 1, (7, 8, 15, 20): 1, (6, 9, 28, 29): 1, (6, 9, 31, 29): 1, (3, 10, 5, 2): 2, (3, 11, 5, 3): 1, (11, 12, 5, 6): 1, (11, 13, 11, 12): 1, (9, 11, 29, 18): 1, (0, 14, 0, 1): 1, (3, 14, 5, 14): 1, (3, 15, 5, 17): 1, (3, 15, 5, 21): 1, (2, 16, 4, 25): 2, (4, 16, 7, 25): 5, (3, 16, 5, 25): 2, (4, 4, 7, 7): 5, (4, 17, 7, 29): 1, (17, 18, 30, 31): 2, (2, 4, 4, 7): 1, (4, 18, 7, 31): 1, (4, 14, 7, 42): 1, (3, 18, 5, 31): 3, (4, 14, 7, 50): 1, (4, 14, 7, 59): 1, (4, 19, 7, 61): 2, (4, 14, 7, 70): 1, (4, 20, 7, 72): 1, (10, 20, 2, 72): 1, (0, 21, 0, 1): 1, (0, 21, 0, 4): 1, (21, 22, 6, 7): 1, (22, 23, 12, 13): 1, (4, 23, 7, 16): 1, (22, 24, 23, 24): 1, (19, 24, 61, 29): 1, (3, 26, 5, 10): 1, (2, 7, 4, 15): 1, (7, 25, 15, 8): 1, (19, 26, 61, 16): 1, (0, 27, 0, 1): 3, (27, 28, 1, 2): 1, (25, 28, 8, 2): 1, (10, 28, 2, 8): 2, (27, 29, 1, 2): 1, (3, 29, 5, 5): 1, (2, 30, 4, 9): 1, (4, 30, 7, 16): 1, (4, 31, 7, 18): 1, (27, 29, 1, 24): 1, (0, 16, 0, 25): 2, (16, 27, 25, 1): 2, (16, 32, 25, 13): 1, (27, 32, 1, 13): 2, (4, 27, 7, 1): 1, (10, 16, 2, 25): 1, (4, 33, 7, 24): 1, (0, 34, 0, 1): 1, (3, 34, 5, 3): 1, (3, 35, 5, 6): 1, (27, 35, 1, 9): 1, (25, 32, 8, 13): 1, (3, 37, 5, 27): 1, (2, 20, 4, 72): 1, (5, 20, 8, 72): 1, (2, 5, 4, 8): 1, (4, 36, 7, 13): 1, (9, 37, 29, 30): 1}

	assert frequencies == check

def test_mol_bond_frequencies_11_20_1():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database11-20.sdf', max_n=1)

	for key, val in frequencies.items():
		print(key,val)

	print(frag_frequencies)

	print(print_fragments(fragment_database))

	print(frequencies)

	check = {(0, 1, 0, 1): 2, (1, 2, 1, 2): 1, (2, 3, 4, 5): 1, (3, 4, 5, 6): 2, (4, 4, 6, 6): 2, (3, 5, 5, 10): 1, (5, 6, 14, 15): 1, (1, 5, 1, 22): 1}

	assert frequencies == check

	check = '''[0, 34, 35, 36] ['C', 'H', 'H', 'H']
[1] ['O']
[32, 33, 2, 3, 4, 37, 55, 56, 57, 58, 27, 28, 29, 30, 31] ['C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'N', 'C', 'C', 'C']
[5, 38] ['N', 'H']
[40, 6, 39] ['H', 'C', 'H']
[10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 25, 26, 46, 47, 48, 49, 50, 54] ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
[15] ['Cl']
'''

	assert print_fragments(fragment_database) == check

	check = [2, 2, 1, 2, 3, 1, 1]

	assert frag_frequencies == check

def test_mol_bond_frequencies_11_20_2():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database11-20.sdf', max_n=2)

	for key, val in frequencies.items():
		print(key,val)

	print(frag_frequencies)

	print(print_fragments(fragment_database))

	print(frequencies)

	assert frag_frequencies ==[3, 2, 1, 5, 7, 1, 1, 1, 1, 2, 1, 5, 2, 3]

	assert frequencies == {(0, 1, 0, 1): 2, (1, 2, 1, 2): 1, (2, 3, 4, 5): 1, (3, 4, 5, 6): 3, (4, 4, 6, 6): 3, (3, 5, 5, 10): 1, (5, 6, 14, 15): 1, (1, 5, 1, 22): 1, (0, 7, 0, 1): 1, (4, 7, 6, 1): 1, (7, 10, 1, 15): 1, (4, 8, 6, 3): 1, (8, 9, 8, 9): 1, (8, 9, 11, 9): 1, (10, 11, 15, 19): 1, (3, 11, 5, 19): 3, (3, 12, 5, 22): 2, (4, 12, 6, 22): 2, (11, 12, 19, 22): 2, (4, 11, 6, 19): 1, (11, 13, 19, 33): 3}

	check = '''[0, 34, 35, 36] ['C', 'H', 'H', 'H']
[1] ['O']
[32, 33, 2, 3, 4, 37, 55, 56, 57, 58, 27, 28, 29, 30, 31] ['C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'N', 'C', 'C', 'C']
[5, 38] ['N', 'H']
[40, 6, 39] ['H', 'C', 'H']
[10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 25, 26, 46, 47, 48, 49, 50, 54] ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
[15] ['Cl']
[1] ['N']
[3, 4, 5, 6, 7, 8, 10, 11, 13, 14, 47] ['C', 'C', 'N', 'C', 'N', 'C', 'N', 'C', 'C', 'N', 'H']
[48, 9, 49] ['H', 'N', 'H']
[68, 69, 40, 41, 15, 16, 17, 18, 52, 53] ['H', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H']
[19, 20] ['C', 'O']
[22, 55] ['C', 'H']
[65, 33] ['H', 'O']
'''

	assert print_fragments(fragment_database) == check

def test_mol_bond_frequencies_11_20_3():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database11-20.sdf', max_n=3)

	for key, val in frequencies.items():
		print(key,val)

	print(frag_frequencies)

	print(print_fragments(fragment_database))

	print(frequencies)

	assert frag_frequencies == [4, 2, 1, 5, 9, 1, 2, 1, 1, 3, 1, 5, 2, 3, 1, 1, 1]

	assert frequencies == {(0, 1, 0, 1): 2, (1, 2, 1, 2): 1, (2, 3, 4, 5): 1, (3, 4, 5, 6): 3, (4, 4, 6, 6): 3, (3, 5, 5, 10): 1, (5, 6, 14, 15): 1, (1, 5, 1, 22): 1, (0, 7, 0, 1): 1, (4, 7, 6, 1): 1, (7, 10, 1, 15): 1, (4, 8, 6, 3): 1, (8, 9, 8, 9): 1, (8, 9, 11, 9): 1, (10, 11, 15, 19): 1, (3, 11, 5, 19): 3, (3, 12, 5, 22): 2, (4, 12, 6, 22): 2, (11, 12, 19, 22): 2, (4, 11, 6, 19): 1, (11, 13, 19, 33): 3, (0, 14, 0, 1): 1, (4, 14, 6, 1): 1, (4, 9, 6, 9): 1, (14, 15, 7, 8): 1, (6, 15, 15, 11): 1, (4, 15, 6, 14): 1, (4, 16, 6, 16): 1}

	check = '''[0, 34, 35, 36] ['C', 'H', 'H', 'H']
[1] ['O']
[32, 33, 2, 3, 4, 37, 55, 56, 57, 58, 27, 28, 29, 30, 31] ['C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'N', 'C', 'C', 'C']
[5, 38] ['N', 'H']
[40, 6, 39] ['H', 'C', 'H']
[10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 25, 26, 46, 47, 48, 49, 50, 54] ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
[15] ['Cl']
[1] ['N']
[3, 4, 5, 6, 7, 8, 10, 11, 13, 14, 47] ['C', 'C', 'N', 'C', 'N', 'C', 'N', 'C', 'C', 'N', 'H']
[48, 9, 49] ['H', 'N', 'H']
[68, 69, 40, 41, 15, 16, 17, 18, 52, 53] ['H', 'H', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H']
[19, 20] ['C', 'O']
[22, 55] ['C', 'H']
[65, 33] ['H', 'O']
[1, 2, 3, 4, 7] ['C', 'N', 'N', 'C', 'N']
[8, 9, 10, 11, 13, 14, 29, 30, 31] ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H']
[34, 35, 36, 37, 38, 16, 17, 18, 19, 20, 21] ['H', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C']
'''

	assert print_fragments(fragment_database) == check


def test_split_mol():

	mol = molecule_from_sdf('mol-1.sdf')

	fragments, pairs, single_bonds = get_fragments_dataset(mol)

	valences = []
	elements = []

	for frag in fragments:
		for i in frag.nodes:
			valences.append(frag.nodes[i]['valence'])
			elements.append(frag.nodes[i]['element'])

	check = [4, 1, 1, 1, 4, 4, 4, 3, 2, 1, 1, 3, 4, 2, 1, 1, 4, 6, 2, 2, 4, 4, 4, 3, 1, 1, 1, 4, 4, 4, 4, 4, 1, 1, 4, 1, 4, 4, 4, 4, 4, 1, 1, 1, 4, 1, 1, 1, 4, 4, 3, 3, 4, 1, 1]

	assert valences == check

	check = ['C', 'H', 'H', 'H', 'C', 'C', 'C', 'N', 'O', 'H', 'H', 'N', 'C', 'O', 'H', 'H', 'C', 'S', 'O', 'O', 'C', 'C', 'C', 'N', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'C', 'H', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'C', 'N', 'N', 'C', 'F', 'F']

	assert elements == check

def test_valence_make_fragment_database():

	fragment_database, frequencies, frag_frequencies = make_fragment_database('../datasets/database1000/database1.sdf')

	elements = []
	valences = []

	for i in range(len(fragment_database)):

		frag = fragment_database[i]

		for j in frag.nodes:
			elements.append(frag.nodes[j]['element'])
			valences.append(frag.nodes[j]['valence'])				


	check = [4, 1, 1, 1, 4, 4, 4, 3, 2, 1, 1, 3, 4, 2, 1, 1, 4, 6, 2, 2, 4, 4, 4, 3, 1, 1, 1, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 1, 1, 1, 4, 1, 1, 1, 4, 4, 3, 3, 4, 1]

	assert valences == check

	check = ['C', 'H', 'H', 'H', 'C', 'C', 'C', 'N', 'O', 'H', 'H', 'N', 'C', 'O', 'H', 'H', 'C', 'S', 'O', 'O', 'C', 'C', 'C', 'N', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'C', 'N', 'N', 'C', 'F']

	assert elements == check	


