import sys,os

from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *

def combine_fragment_databases(fragment_database, frequencies, frag_frequencies, frag_mapping, fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out):

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

		#get mapping for atom numbers
		frag1_map = frag_mapping2to1[frag1_index]
		frag2_map = frag_mapping2to1[frag2_index]

		#convert frag indeces 
		frag1_index = frag_index_mapping[key[0]]
		frag2_index = frag_index_mapping[key[1]]

		update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond, val)

	save_frequencies_txt(frequencies, frequencies_txt_out)

	save_fragments_sdf(fragment_database, fragments_sdf_out)

	save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt_out)

	save_fragments_txt(fragment_database, fragments_txt_out)

def loop(n, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out):

	fragment_database_mol = get_fragment_database('%s0.sdf' %(fragments_sdf_in))

	fragment_database = []

	for i in fragment_database_mol:
		fragment_database.append(i.graph)

	frequencies = get_bond_frequencies('%s0.txt' %frequencies_txt_in )

	frag_frequencies = get_frag_frequencies('%s0.txt' %frag_frequencies_txt_in)

	frag_mapping = get_frag_mapping('%s0.txt' %fragments_txt_in )

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	for i in range(1, n):

		fragments_sdf_2 = '%s%s.sdf' %(fragments_sdf_in, i)
		fragments_txt_2 = '%s%s.txt' %(fragments_txt_in, i) 
		frequencies_txt_2 = '%s%s.txt' %(frequencies_txt_in, i) 
		frag_frequencies_txt_2 = '%s%s.txt' %(frag_frequencies_txt_in, i)

		print(fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out)

		combine_fragment_databases(fragment_database, frequencies, frag_frequencies, frag_mapping, fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out)

def renumber_frequencies(fragments_txt_in, frequencies_txt_in, frequencies_txt_out):

	frequencies = get_bond_frequencies(frequencies_txt_in)

	frag_mapping = get_frag_mapping(fragments_txt_in)

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	save_frequencies_txt(frequencies, frequencies_txt_out)

if __name__ == '__main__':


	n = int(sys.argv[1])

	in_sub = sys.argv[2]
	out_sub = sys.argv[3]

	fragments_sdf_in = 'fragments%s_' %in_sub
	fragments_txt_in = 'fragments%s_' %in_sub
	frequencies_txt_in = 'frequencies%s_' %in_sub
	frag_frequencies_txt_in = 'frag_frequencies%s_' %in_sub
	fragments_sdf_out = 'fragments%s.sdf' %out_sub
	fragments_txt_out = 'fragments%s.txt' %out_sub
	frequencies_txt_out = 'frequencies%s.txt' %out_sub
	frag_frequencies_txt_out = 'frag_frequencies%s.txt' %out_sub

	loop(n, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out)


"""
	
	fragments_sdf_in = sys.argv[2]
	fragments_txt_in = sys.argv[3]
	frequencies_txt_in = sys.argv[4]
	frag_frequencies_txt_in = sys.argv[5]
	fragments_sdf_out = sys.argv[6]
	fragments_txt_out = sys.argv[7]
	frequencies_txt_out = sys.argv[8]
	frag_frequencies_txt_out = sys.argv[9]
"""