import sys,os
import networkx
import time

from pymolgen.molecule_formats import *
from pymolgen.generate import *
from pymolgen.molecule import *
from networkx.algorithms import isomorphism

home = os.path.expanduser("~")

def get_fragments_dataset(mol):

	single_bonds = mol.get_single_bonds_not_h_not_c()

	fragments = split_mol(mol, single_bonds)

	pairs = get_pairs(single_bonds, fragments)

	return fragments, pairs, single_bonds

def print_fragments(fragments):

	out = ''

	for fragment in fragments:

		atoms = []

		for i in fragment.nodes:
			atoms.append(fragment.nodes[i]['element'])

		out += str(fragment.nodes) + ' ' + str(atoms) + '\n'

	return out


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

def get_fragment_index(fragment, fragment_database, fragment_database_len=None, atom_list_all=None):

	is_new = True

	index = []

	map = {}

	for i in fragment.nodes:
		map[i] = i

	fragment_len = len(fragment)

	fragment_atom_list = get_atom_list(fragment)

	for i in range(len(fragment_database)):

		if fragment_database_len is not None:
			fragment_database_len_i = fragment_database_len[i]
		else:
			fragment_database_len_i = len(fragment_database[i])

		if atom_list_all is not None:
			atom_list_all_i = atom_list_all[i]
		else:
			atom_list_all_i = get_atom_list(fragment_database[i])

		if fragment_len == fragment_database_len_i and fragment_atom_list == atom_list_all_i:

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

def save_fragment_sdf(fragment, fragments_sdf):

	outfile = open(fragments_sdf, 'a')

	mol = Molecule()
	mol.graph = fragment.copy()
	mol.set_valence_from_bonds()

	lines = molecule_to_sdf(mol)

	for line in lines:
		outfile.write(line)
	outfile.write('$$$$\n')

	outfile.close()

def get_atom_list(fragment):

	atom_list = []

	for i in fragment.nodes:
		atom_list.append(fragment.nodes[i]["element"])

	atom_list.sort()

	return atom_list

def to_np_matrix(fragments):
	for fragment in fragments:
		matrix = networkx.to_numpy_matrix(fragment)

def update_database(pair, bond, fragment_database, fragments, frequencies, fragments_sdf=None, fragment_database_len=None, atom_list_all=None):

	frag1 = pair[0]
	frag2 = pair[1]

	frag1_bond = bond[0]
	frag2_bond = bond[1]

	frag1_is_new, frag1_index, frag1_map = get_fragment_index(fragments[frag1], fragment_database, fragment_database_len, atom_list_all)
	
	if frag1_is_new: 

		fragment_database.append(fragments[frag1])

		#fragment_database_len list to increase performance of get_fragment_index
		if fragment_database_len is not None:
			fragment_database_len.append(len(fragments[frag1]))
		#atom_list_all to increase performance of get_fragment_index
		if atom_list_all is not None:
			atom_list_all.append(get_atom_list(fragments[frag1]))

		if fragments_sdf is not None:
			save_fragment_sdf(fragments[frag1], fragments_sdf)

	frag2_is_new, frag2_index, frag2_map = get_fragment_index(fragments[frag2], fragment_database, fragment_database_len, atom_list_all)

	if frag2_is_new: 
		fragment_database.append(fragments[frag2])

		#fragment_database_len list to increase performance of get_fragment_index
		if fragment_database_len is not None:
			fragment_database_len.append(len(fragments[frag2]))
		#atom_list_all to increase performance of get_fragment_index
		if atom_list_all is not None:
			atom_list_all.append(get_atom_list(fragments[frag2]))

		if fragments_sdf is not None:
			save_fragment_sdf(fragments[frag2], fragments_sdf)

	update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond)


def update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond):

	#order fragments according to their indeces in the fragment database into i,j
	#then order atoms making bonds into k,l
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

	#increase frequencies of (i,j,k,l) if already in database
	for key in frequencies.keys():
		if (i,j,k,l) == key:
			frequencies[(i,j,k,l)] += 1	
			return

	#if (i,j,k,l) is new make new entry into dictionary
	frequencies[(i,j,k,l)] = 1

def make_fragment_database(database_file, fragments_sdf, fragments_txt, frequencies_txt):

	outfile = open(fragments_sdf, 'w')

	dataset = SDFDatasetLarge(database_file)

	fragment_database = []
	fragment_database_len = []

	atom_list_all = []

	frequencies = {}
	t0 = time.time()
	counter = 0

	#loop through every molecule in the dataset
	for i in range(len(dataset)):

		counter += 1

		if counter % 10 == 0: 
			print(counter, time.time() - t0)
			t0 = time.time()

		#load new molecule from database
		mol = dataset.load_molecule(i)

		#split molecule and get fragments, pairs means pairs of fragments bonded together, and bonds is bonds between atoms of each fragment
		fragments, pairs, bonds = get_fragments_dataset(mol)

		#update database in pairs of fragments by evaluating if 1. each fragment exist, 2. if a pair exists
		#then update fragments and/or bonds between fragments and frequencies accordingly
		for i in range(len(pairs)):
			update_database(pair=pairs[i], bond=bonds[i], fragment_database=fragment_database, frequencies=frequencies, 
				fragments=fragments, fragments_sdf=fragments_sdf, fragment_database_len=fragment_database_len, 
				atom_list_all=atom_list_all)

	with open(fragments_txt, 'w') as outfile:
		outfile.write(print_fragments(fragment_database))

	with open(frequencies_txt, 'w') as outfile:
		for key, val in frequencies.items():
			outfile.write(f"{str(key)}: {str(val)}\n")

if __name__ == '__main__':

	database_file = sys.argv[1]
	fragments_sdf = sys.argv[2]
	fragments_txt = sys.argv[3]
	frequencies_txt = sys.argv[4]

	make_fragment_database(database_file, fragments_sdf, fragments_txt, frequencies_txt)