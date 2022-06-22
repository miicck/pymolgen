from pymolgen.molecule_formats import *
from pymolgen.generate import *
from pymolgen.molecule import *
from networkx.algorithms import isomorphism
import networkx
import time

from os.path import expanduser
home = expanduser("~")

def get_fragments_dataset(mol):

	single_bonds = mol.get_single_bonds_not_h_not_c()

	fragments = split_mol(mol, single_bonds)

	pairs = get_pairs(single_bonds, fragments)

	return fragments, pairs, single_bonds

def print_fragments(fragments):

	for fragment in fragments:
		print(fragment.nodes)


def split_mol(mol, bonds):

	for bond in bonds:
		mol.graph.remove_edge(bond[0], bond[1])

	new = [mol.graph.subgraph(c) for c in networkx.connected_components(mol.graph)]

	return new

def get_pairs(bonds, fragments):

	pairs = []

	for bond in bonds:
		a = bond[0]
		b = bond[1]

		for i in range(len(fragments)):
			if a in fragments[i]:
				index_a = i
			if b in fragments[i]:
				index_b = i

		pairs.append([index_a, index_b])

	return pairs


def is_fragment_new(fragment, fragment_database):

	for i in fragment_database:
		if networkx.is_isomorphic(i,fragment):
			return False

	return True

def get_fragment_index(fragment, fragment_database, fragment_database_len, atom_list_all):

	is_new = True

	index = []

	map = {}

	for i in fragment.nodes:
		map[i] = i

	fragment_len = len(fragment)

	fragment_atom_list = get_atom_list(fragment)

	for i in range(len(fragment_database)):

		if fragment_len == fragment_database_len[i] and fragment_atom_list == atom_list_all[i]:

			gm = isomorphism.GraphMatcher(fragment, fragment_database[i])

			if gm.is_isomorphic():
				index.append(i)
				is_new = False
				map = gm.mapping

	if len(index) > 1:
		raise Exception('fragment in fragment_database more than once')

	if is_new is True:
		index.append(len(fragment_database))

	return is_new, index[0], map

def save_fragments_sdf(fragments, outfile_name):

	outfile = open(outfile_name, 'w')

	for fragment in fragments:
		mol = Molecule()
		mol.graph = fragment.copy()
		mol.set_valence_from_bonds()

		lines = molecule_to_sdf(mol)

		for line in lines:
			outfile.write(line)
		outfile.write('$$$$\n')

	outfile.close()

def make_fragment_database():

	

	dataset = SDFDatasetLarge(home + '/Downloads/chembl_30.sdf')

	fragment_database = []
	fragment_database_len = []

	atom_list_all = []
	bond_list_all = []

	frequencies = {}
	t0 = time.time()
	counter = 0
	for i in range(len(dataset)):

		counter += 1

		if counter % 100 == 0: 
			print(counter, time.time() - t0)
			t0 = time.time()

		mol = dataset.load_molecule(i)
		fragments, pairs, bonds = get_fragments_dataset(mol)
		
		for i in range(len(pairs)):

			frag1 = pairs[i][0]
			frag2 = pairs[i][1]

			frag1_bond = bonds[i][0]
			frag2_bond = bonds[i][1]

			frag1_is_new, frag1_index, frag1_map = get_fragment_index(fragments[frag1], fragment_database, fragment_database_len, atom_list_all)
			
			if frag1_is_new: 

				fragment_database.append(fragments[frag1])
				fragment_database_len.append(len(fragments[frag1]))
				atom_list_all.append(get_atom_list(fragments[frag1]))

			frag2_is_new, frag2_index, frag2_map = get_fragment_index(fragments[frag2], fragment_database, fragment_database_len, atom_list_all)

			if frag2_is_new: 
				fragment_database.append(fragments[frag2])
				fragment_database_len.append(len(fragments[frag2]))
				atom_list_all.append(get_atom_list(fragments[frag2]))

			update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond)


	print_fragments(fragment_database)

	print(frequencies)

	save_fragments_sdf(fragment_database, 'fragments3.sdf')

def get_atom_list(fragment):

	atom_list = []

	for i in fragment.nodes:
		atom_list.append(fragment.nodes[i]["element"])

	atom_list.sort()

	return atom_list

def to_np_matrix(fragments):
	for fragment in fragments:
		matrix = networkx.to_numpy_matrix(fragment)

def update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond):

	if frag1_index < frag2_index:
		i = frag1_index
		j = frag2_index
		k = frag1_map[frag1_bond]
		l = frag2_map[frag2_bond]

	else:
		i = frag2_index
		j = frag1_index
		k = frag2_map[frag2_bond]
		l = frag1_map[frag1_bond]

	for key in frequencies.keys():
		if (i,j,k,l) == key:
			frequencies[(i,j,k,l)] += 1	
			return

	frequencies[(i,j,k,l)] = 1

if __name__ == '__main__':
	make_fragment_database()