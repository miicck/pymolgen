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

	if single_bonds is False:
		return False, False, False

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

def compound_dict(dict1, dict2):

	d = {}

	for key, val in dict1.items():
		d[key] = dict2[val]

	return d


def node_compare_element(node_1, node_2):
	return node_1["element"] == node_2["element"] and node_1["hybridization"] == node_2["hybridization"]

def get_fragment_index(fragment, fragment_database, fragment_database_len=None, atom_list_all=None):

	is_new = True

	index = []

	map = get_canonical_mapping(fragment)

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

			gm = isomorphism.GraphMatcher(fragment, fragment_database[i], node_match=node_compare_element)

			if gm.is_isomorphic():
				index.append(i)
				is_new = False
				newmap = gm.mapping

				map = get_canonical_mapping(fragment_database[i])

				map = compound_dict(newmap, map)

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

def update_bond_database(pair, bond, fragment_database, fragments, frequencies, frag_frequencies, fragments_sdf=None, fragment_database_len=None, atom_list_all=None):

	frag1 = pair[0]
	frag2 = pair[1]

	frag1_bond = bond[0]
	frag2_bond = bond[1]

	frag1_is_new, frag1_index, frag1_map = get_fragment_index(fragments[frag1], fragment_database, fragment_database_len, atom_list_all)
	
	if frag1_is_new: 

		frag_frequencies.append(1)

		fragment_database.append(fragments[frag1])

		#fragment_database_len list to increase performance of get_fragment_index
		if fragment_database_len is not None:
			fragment_database_len.append(len(fragments[frag1]))
		#atom_list_all to increase performance of get_fragment_index
		if atom_list_all is not None:
			atom_list_all.append(get_atom_list(fragments[frag1]))

		if fragments_sdf is not None:
			save_fragment_sdf(fragments[frag1], fragments_sdf)

	else:
		frag_frequencies[frag1_index] += 1

	frag2_is_new, frag2_index, frag2_map = get_fragment_index(fragments[frag2], fragment_database, fragment_database_len, atom_list_all)

	if frag2_is_new:

		frag_frequencies.append(1)

		fragment_database.append(fragments[frag2])

		#fragment_database_len list to increase performance of get_fragment_index
		if fragment_database_len is not None:
			fragment_database_len.append(len(fragments[frag2]))
		#atom_list_all to increase performance of get_fragment_index
		if atom_list_all is not None:
			atom_list_all.append(get_atom_list(fragments[frag2]))

		if fragments_sdf is not None:
			save_fragment_sdf(fragments[frag2], fragments_sdf)


	else:

		frag_frequencies[frag2_index] += 1

	update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond)


def update_fragment_database(fragment_database, fragment, frequencies, frag_frequencies, fragments_sdf=None, fragment_database_len=None, atom_list_all=None):

	frag1_is_new, frag1_index, frag1_map = get_fragment_index(fragment, fragment_database, fragment_database_len, atom_list_all)

	if frag1_is_new: 

		frag_frequencies.append(1)

		fragment_database.append(fragment)

		#fragment_database_len list to increase performance of get_fragment_index
		if fragment_database_len is not None:
			fragment_database_len.append(len(fragment))

		#atom_list_all to increase performance of get_fragment_index
		if atom_list_all is not None:
			atom_list_all.append(get_atom_list(fragment))

		if fragments_sdf is not None:
			save_fragment_sdf(fragment, fragments_sdf)

	else:
		frag_frequencies[frag1_index] += 1

	return frag1_index, frag1_map


def update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond, val=None):

	if val is None:
		val = 1


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
			frequencies[(i,j,k,l)] += val	
			return

	#if (i,j,k,l) is new make new entry into dictionary
	frequencies[(i,j,k,l)] = val

def make_fragment_database(database_file, fragments_sdf=None, fragments_txt=None, frequencies_txt=None, frag_frequencies_txt=None, max_n=None, verbose=False, fragment_database=None, frequencies=None, frag_frequencies=None):

	if fragments_sdf is not None:
		outfile = open(fragments_sdf, 'w')

	dataset = SDFDatasetLarge(database_file, max_n)

	if fragment_database is None:
		fragment_database = []

	fragment_database_len = []

	for i in fragment_database:
		fragment_database_len.append(len(i))

	atom_list_all = []

	for i in fragment_database:
		atom_list_all.append(get_atom_list(i))

	if frequencies is None:
		frequencies = {}
	
	t0 = time.time()
	counter = 0

	if frag_frequencies is None:
		frag_frequencies = []

	#loop through every molecule in the dataset
	for i in range(len(dataset)):

		counter += 1

		if verbose: print('%s ' %counter, end='', flush=True)

		if counter % 100 == 0: 
			if verbose: print()
			print('TIME', counter, time.time() - t0, flush=True)
			t0 = time.time()

		if counter % 5000 == 0:
			if verbose:
				print(print_fragments(fragment_database), flush=True)
				print(frequencies, flush=True)
				print(frag_frequencies, flush=True)

			save_fragments_txt(fragment_database, fragments_txt)
			save_frequencies_txt(frequencies, frequencies_txt)
			save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt)

		#load new molecule from database
		mol = dataset.load_molecule(i)

		#split molecule and get fragments, pairs means pairs of fragments bonded together, and bonds is bonds between atoms of each fragment
		fragments, pairs, bonds = get_fragments_dataset(mol)
		if fragments == False:
			continue
		frag1_index_list = []
		frag1_map_list = []

		#update database of fragments, if fragment exists increase frag_frequency, otherwise add fragment
		for fragment in fragments:
			frag1_index, frag1_map = update_fragment_database(fragment_database, fragment, frequencies, frag_frequencies, fragments_sdf, fragment_database_len, atom_list_all)

			frag1_index_list.append(frag1_index)
			frag1_map_list.append(frag1_map)

		for i in range(len(pairs)):

			frag1 = pairs[i][0]
			frag2 = pairs[i][1]

			frag1_index = frag1_index_list[frag1]
			frag1_map = frag1_map_list[frag1]

			frag1_bond = bonds[i][0]
			frag2_bond = bonds[i][1]

			frag2_is_new, frag2_index, frag2_map = get_fragment_index(fragments[frag2], fragment_database, fragment_database_len, atom_list_all)

			update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond)

	save_fragments_txt(fragment_database, fragments_txt)

	#frag_mapping = get_frag_mapping(fragments_txt)

	#frequencies = update_bond_frequencies(frequencies, frag_mapping)

	save_frequencies_txt(frequencies, frequencies_txt)
	save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt)

	return fragment_database, frequencies, frag_frequencies

def get_frag_mapping(fragments_txt):
    """
    Generates atom mapping as dictionary from original atom numbers to 0 - len(atoms)
    Returns list of dictionaries as mapping for each fragment
    """
    frag_mapping = []

    with open(fragments_txt) as infile:
        for line in infile:
            atoms = line.split(']')[0].strip('[').split(',')
            atoms = [int(x.strip()) for x in atoms]
            d = {}
            for i in range(len(atoms)):
                d[atoms[i]] = i
            frag_mapping.append(d)

    return frag_mapping

def update_bond_frequencies(bond_frequencies, frag_mapping):
    """
    Update bond frequencies for atom numbering in frag_mapping (typically numbering from 0)
    """
    d = {}

    for key, val in bond_frequencies.items():
        print(key, val)

        i = key[0]
        j = key[1]

        k = frag_mapping[i][key[2]]
        l = frag_mapping[j][key[3]]

        d[(i,j,k,l)] = val

    return d

def save_fragments_txt(fragment_database, fragments_txt):
	if fragments_txt is not None:
		with open(fragments_txt, 'w') as outfile:
			outfile.write(print_fragments(fragment_database))

def save_frequencies_txt(frequencies, frequencies_txt):
	if frequencies_txt is not None:
		with open(frequencies_txt, 'w') as outfile:
			for key, val in frequencies.items():
				outfile.write(f"{str(key)}: {str(val)}\n")

def save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt):
	if frag_frequencies_txt is not None:
		with open(frag_frequencies_txt, 'w') as outfile:
			for i in range(len(frag_frequencies)):
				outfile.write('%s ' %frag_frequencies[i])
				if i + 1 % 10 == 0: outfile.write('\n')
			outfile.write('\n')

def map_mols(mol1, mol2):

	gm = isomorphism.GraphMatcher(mol1, mol2, node_match=node_compare_element)

	all_mappings = []

	for i in gm.isomorphisms_iter():
		all_mappings.append(i)

	mapping = all_mappings[0]

	return mapping

def get_canonical_mapping(fragment):
	gm = isomorphism.GraphMatcher(fragment, fragment, node_match=node_compare_element)

	all_mappings = []

	for mapping in gm.isomorphisms_iter():
		all_mappings.append(mapping)

	canonical_mapping = all_mappings[0]

	for i in all_mappings:
		for key, val in i.items():
			if canonical_mapping[key] > val:
				canonical_mapping[key] = val

	return canonical_mapping

def renumber_fragment(fragment):

	fragment = networkx.convert_node_labels_to_integers(fragment,first_label=0)

	mapping = get_canonical_mapping(fragment)

	fragment = networkx.relabel_nodes(fragment, mapping)

	return fragment


if __name__ == '__main__':

	database_file = sys.argv[1]
	fragments_sdf = sys.argv[2]
	fragments_txt = sys.argv[3]
	frequencies_txt = sys.argv[4]
	frag_frequencies_txt = sys.argv[5]

	make_fragment_database(database_file, fragments_sdf, fragments_txt, frequencies_txt, frag_frequencies_txt)
	print('Normal termination')
