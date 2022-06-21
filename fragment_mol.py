from pymolgen.molecule_formats import *
from pymolgen.generate import *
from pymolgen.molecule import *
import networkx

from os.path import expanduser
home = expanduser("~")

def get_fragments():

	mol = molecule_from_sdf(home + '/pymolgen/test/mol-1.sdf')

	single_bonds = mol.get_single_bonds_not_h_not_c()

	new = split_mol(mol, single_bonds)


	mol2 = molecule_from_sdf(home + '/pymolgen/test/mol-1-can.sdf')

	bonds2 = mol2.get_single_bonds_not_h_not_c()

	new2 = split_mol(mol2, bonds2)

	for i in new:
		for j in new2:
			equal = networkx.is_isomorphic(i,j)

			if equal is True:
				print(i.nodes,j.nodes)

def get_fragments_dataset(mol):

	single_bonds = mol.get_single_bonds_not_h_not_c()

	fragments = split_mol(mol, single_bonds)

	return fragments

def print_fragments(fragments):

	for fragment in fragments:
		print(fragment.nodes)


def split_mol(mol, bonds):

	for bond in bonds:
		mol.graph.remove_edge(bond[0], bond[1])

	new = [mol.graph.subgraph(c) for c in networkx.connected_components(mol.graph)]

	return new

def is_fragment_new(fragment, fragment_database):

	for i in fragment_database:
		if networkx.is_isomorphic(i,fragment):
			return False

	return True

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

	dataset = SDFDatasetLarge(home + '/pymolgen/datasets/database1000/database1000.sdf', max_n = 1)

	fragment_database = []

	for i in range(len(dataset)):
		mol = dataset.load_molecule(i)
		fragments = get_fragments_dataset(mol)

		for fragment in fragments:

			print(type(fragment))
			if is_fragment_new(fragment, fragment_database):
				fragment_database.append(fragment)

	print('fragment_database')
	print_fragments(fragment_database)

	save_fragments_sdf(fragment_database, 'fragments.sdf')

make_fragment_database()