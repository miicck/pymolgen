import sys,os

from fragment_mol import *
from fragment_builder import *

def combine_fragment_databases(fragments_sdf_1, fragments_txt_1, frequencies_txt_1, frag_frequencies_txt_1, fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out):

	fragment_database_mol = get_fragment_database(fragments_sdf_1)

	fragment_database = []

	for i in fragment_database_mol:
		fragment_database.append(i.graph)

	frequencies = get_bond_frequencies(frequencies_txt_1)

	frag_frequencies = get_frag_frequencies(frag_frequencies_txt_1)

	frag_mapping = get_frag_mapping(fragments_txt_1)

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	fragment_database_mol_2 = get_fragment_database(fragments_sdf_2)

	fragment_database_2 = []

	for i in fragment_database_mol_2:
		fragment_database_2.append(i.graph)

	frequencies2 = get_bond_frequencies(frequencies_txt_2)

	frag_mapping2 = get_frag_mapping(fragments_txt_2)

	frequencies2 = update_bond_frequencies(frequencies2, frag_mapping2)

	frag_frequencies2 = get_frag_frequencies(frag_frequencies_txt_2)

	#mapping of fragment atom indeces from 2 to 1 (or 2 to 2 if new fragment)
	frag_mapping2to1 = []

	#map of fragment index from 2 to final database
	frag_index_mapping = []

	for i in range(len(fragment_database_2)):

		fragment = fragment_database_2[i]

		frag1_is_new, frag1_index, frag1_map = get_fragment_index(fragment, fragment_database)

		if frag1_is_new: 

			frag_frequencies.append(frag_frequencies2[i])

			frag_index_mapping.append(len(fragment_database))

			fragment_database.append(fragment)

		else:
			frag_frequencies[frag1_index] += frag_frequencies2[i]

			frag_index_mapping.append(frag1_index)

		frag_mapping2to1.append(frag1_map)

	for key, val in frequencies2.items():

		frag1_index = key[0]
		frag2_index = key[1]
		frag1_bond = key[2]
		frag2_bond = key[3]

		frag1_map = frag_mapping2to1[frag1_index]
		frag2_map = frag_mapping2to1[frag2_index]

		frag1_index = frag_index_mapping[key[0]]
		frag2_index = frag_index_mapping[key[1]]

		update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond, val)

	save_frequencies_txt(frequencies, frequencies_txt_out)

	save_fragments_sdf(fragment_database, fragments_sdf_out)

	save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt_out)


if __name__ == '__main__':

	fragments_sdf_1 = sys.argv[1]
	fragments_txt_1 = sys.argv[2]
	frequencies_txt_1 = sys.argv[3]
	frag_frequencies_txt_1 = sys.argv[4]
	fragments_sdf_2 = sys.argv[5]
	fragments_txt_2 = sys.argv[6]
	frequencies_txt_2 = sys.argv[7]
	frag_frequencies_txt_2 = sys.argv[8]
	fragments_sdf_out = sys.argv[9]
	fragments_txt_out = sys.argv[10]
	frequencies_txt_out = sys.argv[11]
	frag_frequencies_txt_out = sys.argv[12]


	combine_fragment_databases(fragments_sdf_1, fragments_txt_1, frequencies_txt_1, frag_frequencies_txt_1, fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out)
