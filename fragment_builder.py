import sys,os

import networkx
from networkx.algorithms import isomorphism

from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *

def node_compare_element(node_1, node_2):
    return node_1["element"] == node_2["element"]

def get_frag_frequencies(frag_frequencies_txt):
    frag_frequencies = []

    with open(frag_frequencies_txt) as infile:
        for line in infile:
            frag_frequencies.extend(int(i) for i in line.split())

    return frag_frequencies

def get_bond_frequencies(bond_frequencies_txt):
    bond_frequencies = {}

    with open(bond_frequencies_txt) as infile:
        for line in infile:
            i = int(line.split()[0].strip('(').strip(','))
            j = int(line.split()[1].strip(','))
            k = int(line.split()[2].strip(','))
            l = int(line.split()[3].strip(':').strip(')').strip(','))
            f = int(line.split()[4])
            bond_frequencies[(i,j,k,l)] = f

    return bond_frequencies

def get_fragment_bond_frequencies(fragment_i, bond_frequencies):

    fragment_bond_frequencies = {}

    for key, val in bond_frequencies.items():
        if fragment_i == key[0] or fragment_i == key[1]:
            fragment_bond_frequencies[key] = val

    return fragment_bond_frequencies

def save_neighbours(fragment_i, fragment_bond_frequencies, fragment_database, outfile_name):

    neighbours = []

    with open(outfile_name, 'w') as outfile:
        pass

    for key, val in fragment_bond_frequencies.items():
        if fragment_i == key[0]: 
            neighbours.append(key[1])
        if fragment_i == key[1]:
            neighbours.append(key[0])

    for i in neighbours:

        mol = fragments_database.load_molecule(i)

        lines = molecule_to_sdf(mol)

        with open(outfile_name, 'a') as outfile:
            for line in lines:
                outfile.write(line)
            outfile.write('$$$$\n')

def make_fragment_database(fragments_sdf):

    fragment_database = SDFDatasetLargeRAM(fragments_sdf)

    return fragment_database

fragments_database = make_fragment_database('fragments50k.sdf')

def find_fragment(fragment, fragment_database):

    for i in range(len(fragment_database)):

        gm = isomorphism.GraphMatcher(fragment.graph, fragment_database[i].graph, node_match=node_compare_element)

        if gm.is_isomorphic():

            return i

    return False

def build_molecule(fragments_sdf, parent_file, parent_fragment_file, remove_hydrogens, remove_hydrogens_parent_fragment):

    fragment_database = make_fragment_database(fragments_sdf)

    mol = molecule_from_sdf(parent_file)

    parent_mw = Molecule.molecular_weight(mol)

    mw = parent_mw

    molecule_fragments = []
    frag_bonds = []

    parent_fragment = molecule_from_sdf(parent_fragment_file)

    for i in remove_hydrogens_parent_fragment:
        parent_fragment = parent_fragment.remove_atom(i)

    parent_fragment_i = find_fragment(parent_fragment, fragment_database)

    print(parent_fragment_i)

    bond_frequencies = get_bond_frequencies('frequencies50k.txt')

    fragment_bond_frequencies = get_fragment_bond_frequencies(parent_fragment_i, bond_frequencies)
    print(fragment_bond_frequencies)

    save_neighbours(parent_fragment_i, fragment_bond_frequencies, fragment_database, 'neighbours.sdf')

    if parent_fragment_i is False:
        sys.exit('Parent fragment not found')

    if remove_hydrogens is not None:
        for i in remove_hydrogens:
            mol = mol.remove_atom(i)

    #while mw < 500:
   #     fragment = fragments_database.random_molecule()


build_molecule('fragments50k.sdf', 'zgwhxzahbbyfix-26.sdf', 'zgwhxzahbbyfix-26-amide-6.sdf', [26], [6])














