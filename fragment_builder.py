import sys,os
import random

import networkx
from networkx.algorithms import isomorphism

from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_mol import print_fragments

def node_compare_element(node_1, node_2):
    return node_1["element"] == node_2["element"]

def get_frag_mapping(fragments_txt):

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

def get_fragment_bond_frequencies(fragment_i, atom_i, bond_frequencies):

    fragment_bond_frequencies = {}

    for key, val in bond_frequencies.items():
        if fragment_i == key[0] and atom_i == key[2]:
            fragment_bond_frequencies[key] = val
        if fragment_i == key[1] and atom_i == key[3]:
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

def get_fragment_database(fragments_sdf):

    fragment_database = SDFDatasetLargeRAM(fragments_sdf)

    return fragment_database

def find_fragment(fragment, fragment_database):

    for i in range(len(fragment_database)):

        gm = isomorphism.GraphMatcher(fragment.graph, fragment_database[i].graph, node_match=node_compare_element)

        if gm.is_isomorphic():

            return i

    return False

def get_random_neighbour(fragment_i, atom_i, fragment_bond_frequencies):

    keys = []
    vals = []

    for key, val in fragment_bond_frequencies.items():
        keys.append(key)
        vals.append(val)

    print('keys =', keys)
    print('vals =', vals)

    draw = random.choices(population=keys, weights=vals, k=1)[0]
    print('fragment_i', fragment_i)
    print('draw =', draw)
    if fragment_i == draw[0]:

        new_frag_i = draw[1]
        fragment_i_atom = draw[2]
        new_frag_i_atom = draw[3]

    if fragment_i == draw[1]:

        new_frag_i = draw[0]
        fragment_i_atom = draw[3]
        new_frag_i_atom = draw[2]

    return new_frag_i, fragment_i_atom, new_frag_i_atom

def update_bond_frequencies(bond_frequencies, frag_mapping):

    d = {}

    for key, val in bond_frequencies.items():

        i = key[0]
        j = key[1]
        k = frag_mapping[i][key[2]]
        l = frag_mapping[j][key[3]]

        d[(i,j,k,l)] = val

    return d

def get_length(list):

    length = 0

    for i in list:
        length += len(i)

    return length

def build_molecule(fragments_sdf, fragments_txt, frequencies_txt, parent_file, parent_fragment_file, remove_hydrogens, remove_hydrogens_parent_fragment):

    #make databases and update atom numberings
    fragment_database = get_fragment_database(fragments_sdf)
    frag_mapping = get_frag_mapping(fragments_txt)
    bond_frequencies = get_bond_frequencies(frequencies_txt)   
    bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)

    mol = molecule_from_sdf(parent_file)

    parent_mw = Molecule.molecular_weight(mol)

    mw = parent_mw

    #prepare parent
    frag_list = []
    frag_mol_list = []
    frag_bond_list = []
    frag_free_valence_list = []

    parent_fragment = molecule_from_sdf(parent_fragment_file)

    for i in remove_hydrogens_parent_fragment:
        parent_fragment = parent_fragment.remove_atom(i)

    parent_fragment_i = find_fragment(parent_fragment, fragment_database)

    parent_fragment = fragment_database[parent_fragment_i]

    frag_free_valence_list.append([])
    print('parent_fragment')
    print_molecule(parent_fragment)
    print(parent_fragment.free_valence_list)
    for i in parent_fragment.free_valence_list:
        frag_free_valence_list[0].append(i)

    if parent_fragment_i is False:
        sys.exit('Parent fragment not found')

    frag_list.append(parent_fragment_i)

    counter = 0
    while get_length(frag_free_valence_list) != 0:

        print('frag_free_valence_list =', frag_free_valence_list)
        counter += 1
        if counter == 100: break

        #choose random position of constituent fragment in molecule
        i = random.randrange(len(frag_list))

        #if chosen fragment has free valence points
        if len(frag_free_valence_list[i]) > 0:

            #get atom from fragment_i to build on
            atom_i = random.choice(frag_free_valence_list[i])

            #get fragment_i (index in fragment_database)
            fragment_i = frag_list[i]

            #get bond frequencies for fragment_i
            fragment_bond_frequencies = get_fragment_bond_frequencies(fragment_i, atom_i, bond_frequencies)
            
            #choose random neighbour
            new_frag_i, fragment_i_atom, new_frag_i_atom = get_random_neighbour(fragment_i, atom_i, fragment_bond_frequencies)
            print('neighbour =', new_frag_i, fragment_i_atom, new_frag_i_atom)

            #add neighbour index in fragment_database to frag_list
            frag_list.append(new_frag_i)

            #generate molecule object from new_frag_i
            new_frag = fragment_database[new_frag_i]

            #add new fragment molecule object to frag_mol_list
            frag_mol_list.append(new_frag)

            print('new mol =', end=''); print_molecule(new_frag)

            #get free valence points of new fragment
            new_free_valence_list = new_frag.free_valence_list

            print('new_frag_free_valence_list =', new_frag.free_valence_list)

            #add bond betweent current fragment and new fragment to list of bonds between fragments (frag_bond_list)
            frag_bond_list.append((fragment_i, new_frag_i, fragment_i_atom, new_frag_i_atom))

            #remove atom from current fragment making bond to new fragment from frag_free_valence_list[i]
            print('fragment_i_atom', fragment_i_atom)
            print('frag_free_valence_list[i]', frag_free_valence_list[i])
            frag_free_valence_list[i].remove(fragment_i_atom)

            #remove atom from new fragment making bond to current fragment from new fragment's list of free valence points
            print('new_frag_i_atom =',  new_frag_i_atom)
            new_free_valence_list.remove(new_frag_i_atom)

            #add new_free_valence_list to the list of available valence points in molecule being built
            frag_free_valence_list.append(new_free_valence_list)

if __name__ == '__main__':

    random.seed(100)

    build_molecule('fragments10.sdf', 'fragments10.txt', 'frequencies10.txt', 'zgwhxzahbbyfix-26.sdf', 'methane.sdf', [26], [3])














