import sys,os
import random
import numpy as np
import argparse
import subprocess

import networkx
from networkx.algorithms import isomorphism

from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_mol import print_fragments, get_canonical_mapping, map_mols, get_frag_mapping, update_bond_frequencies

def node_compare_element(node_1, node_2):
    return node_1["element"] == node_2["element"] and node_1["hybridization"] == node_2["hybridization"]

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

def get_random_neighbour(fragment_i, fragment_bond_frequencies):

    keys = []
    vals = []

    for key, val in fragment_bond_frequencies.items():
        keys.append(key)
        vals.append(val)

    if len(fragment_bond_frequencies) == 0:
        print('fragment bond frequencies =', fragment_bond_frequencies)
        print('fragment_i =', fragment_i)

    draw = random.choices(population=keys, weights=vals, k=1)[0]

    if fragment_i == draw[0]:

        new_frag_i = draw[1]
        fragment_i_atom = draw[2]
        new_frag_i_atom = draw[3]

    if fragment_i == draw[1]:

        new_frag_i = draw[0]
        fragment_i_atom = draw[3]
        new_frag_i_atom = draw[2]

    return new_frag_i, new_frag_i_atom

def get_length(list):

    length = 0

    for i in list:
        length += len(i)

    return length

def reverse_canonical_mapping(fragment):
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

def build_molecule(fragments_sdf, fragments_txt, frequencies_txt, parent_file, parent_fragment_file, remove_hydrogens, remove_hydrogens_parent_fragment, outfile_name, n_mol, filters=False, unique=False, figure=None, rules=False):

    # build pains_database if using filters
    if filters:
        from pymolgen.newmol import gen_pains_database
        try:
            pains_database = gen_pains_database()
        except:
            raise Exception("Could not generate pains database")

    #make databases and update atom numberings
    fragment_database = get_fragment_database(fragments_sdf)
    frag_mapping = get_frag_mapping(fragments_txt)
    bond_frequencies = get_bond_frequencies(frequencies_txt)   
    bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)

    if unique:
        candidate_list = []
        candidate_bond_list = []
    else:
        candidate_list = None
        candidate_bond_list = None

    parent_mol = molecule_from_sdf(parent_file)

    for i in remove_hydrogens:
        parent_mol = parent_mol.remove_atom(i)

    parent_mw = Molecule.molecular_weight(parent_mol)

    parent_fragment = molecule_from_sdf(parent_fragment_file)

    smi = molecule_to_smiles(parent_fragment)

    print('Parent fragment', smi)

    for i in remove_hydrogens_parent_fragment:
        parent_fragment = parent_fragment.remove_atom(i)

    parent_fragment_original = parent_fragment

    parent_fragment_i = find_fragment(parent_fragment, fragment_database)

    if parent_fragment_i is False:
        sys.exit('Parent fragment not found')

    parent_fragment = fragment_database[parent_fragment_i]

    print_molecule(parent_fragment_original)
    print_molecule(parent_fragment)

    parent_mapping = map_mols(parent_fragment_original.graph, parent_fragment.graph)

    print('parent_fragment')
    print_molecule(parent_fragment)
    print('parent_fragment.free_valence_list =', parent_fragment.free_valence_list)

    with open(outfile_name, 'w') as outfile:
        print('Writing to', outfile_name)

    if figure is not None:
        with open(figure, 'w') as outfile:
            print('Writing to figure', figure)

    n = 1
    while n <= n_mol:

        mol = build_mol_single(parent_mol, parent_fragment, parent_fragment_i, fragment_database, bond_frequencies, parent_mapping, filters, pains_database, candidate_list, candidate_bond_list, figure, rules)

        if mol is not None:

            smi = molecule_to_smiles(mol)
            mw = mol.molecular_weight()
            print('NEW_CANDIDATE %s %s %.1f' % (n, smi, mw))            

            lines = molecule_to_sdf(mol)

            with open(outfile_name, 'a') as outfile:
                for line in lines:
                    outfile.write(line)

                outfile.write('$$$$\n')

            if figure is not None:

                newatoms = []
                for i in mol.graph.nodes:
                    if i >= 44:
                        newatoms.append(i)

                fig = mol.get_fragment(newatoms)
                #fig.hydrogenate()
                smi = molecule_to_smiles(fig)
                #print('ATTACHED ', smi)
                print_molecule(fig)

                lines = molecule_to_sdf(fig)

                with open(figure, 'a') as outfile:
                    for line in lines:
                        outfile.write(line)

                    outfile.write('$$$$\n')               

            n += 1

def build_mol_single(parent_mol, parent_fragment, parent_fragment_i, fragment_database, bond_frequencies, parent_mapping, filters=False, pains_database=None, candidate_list=None, candidate_bond_list=None, figure=None, rules=False):

    #prepare parent fragment
    frag_list = []
    frag_mol_list = [parent_mol]
    frag_bond_list = []
    frag_free_valence_list = []

    frag_free_valence_list.append([])

    for i in parent_mol.free_valence_list:
        frag_free_valence_list[0].append(i)

    frag_list.append(-1)

    counter = 0
    while get_length(frag_free_valence_list) != 0:

        counter += 1
        if counter == 100:
            return None

        # choose random position of constituent fragments in molecule
        i = random.randrange(len(frag_list))

        # if chosen fragment has free valence points
        if len(frag_free_valence_list[i]) > 0:

            j = len(frag_list)

            # get atom from fragment_i to build on
            atom_i = random.choice(frag_free_valence_list[i])

            # get fragment_i (index in fragment_database)
            fragment_i = frag_list[i]

            if fragment_i == -1:
                fragment_i = parent_fragment_i
                
                # get mol for fragment_i
                fragment_i_mol = fragment_database[fragment_i]

                # get mapped atom_i since fragment_bond_frequencies are stored for canonical atoms
                atom_i_can = parent_mapping[atom_i]

                #print('fragment bond frequencies =', get_fragment_bond_frequencies(fragment_i, atom_i_can, bond_frequencies))


            else:
                # get mol for fragment_i
                fragment_i_mol = fragment_database[fragment_i]

                # get canonical mapping
                canonical_mapping = get_canonical_mapping(fragment_i_mol.graph)

                # get mapped atom_i since fragment_bond_frequencies are stored for canonical atoms
                atom_i_can = canonical_mapping[atom_i]

            # get bond frequencies for fragment_i
            fragment_bond_frequencies = get_fragment_bond_frequencies(fragment_i, atom_i_can, bond_frequencies)

            # return none molecule if fragment_bond_frequencies has length 0 (cannot build on fragment)
            # this shouldn't happen since all fragments come from molecules so they shuold all have bonds
            # but there could be errors in the database
            if len(fragment_bond_frequencies) == 0:
                return None

            if fragment_i == -1:
                print_molecule(fragment_database[fragment_i])

            # choose random neighbour
            new_frag_i, new_frag_i_atom = get_random_neighbour(fragment_i, fragment_bond_frequencies)

            # generate molecule object from new_frag_i
            new_frag = fragment_database[new_frag_i]

            # get free valence points of new fragment
            new_free_valence_list = new_frag.free_valence_list

            if not new_frag.is_fluorine() and new_frag_i_atom in new_free_valence_list:

                # add neighbour index in fragment_database to frag_list
                frag_list.append(new_frag_i)

                # add bond betweent current fragment and new fragment to list of bonds between fragments (frag_bond_list)
                frag_bond_list.append((i, j, atom_i, new_frag_i_atom))

                # remove atom from current fragment making bond to new fragment from frag_free_valence_list[i]
                frag_free_valence_list[i].remove(atom_i)

                # remove atom from new fragment making bond to current fragment from new fragment's list of free valence points
                try: new_free_valence_list.remove(new_frag_i_atom)
                except: 
                    #print_fragments([new_frag])
                    smi = molecule_to_smiles(new_frag)
                    print(smi)
                    print_molecule(new_frag)
                    lines = molecule_to_sdf(new_frag)
                    with open('new_frag.sdf', 'w') as outfile:
                        for line in lines:
                            outfile.write(line)
                    print(lines)
                    raise Exception('Could not remove atom', new_frag_i_atom, 'from', new_free_valence_list, 'fragment_i =', fragment_i, 'new_frag_i=', new_frag_i)
                    return None

                # add new_free_valence_list to the list of available valence points in molecule being built
                frag_free_valence_list.append(new_free_valence_list)

    for i in frag_list[1:]:
        frag_mol_list.append(fragment_database[i])

    if candidate_list is not None:
        if is_new_candidate(frag_list, frag_bond_list, candidate_list, candidate_bond_list) is True:
            candidate_list.append(frag_list)
            candidate_bond_list.append(frag_bond_list)
        else:
            print("Not unique")
            return None

    mol = combine_all_fragments(frag_mol_list, frag_list, frag_bond_list)

    if filters:
        from pymolgen.newmol import filters_final_mol
        
        try:
            filter_pass = filters_final_mol(mol, pains_database)
        except:
            filter_pass = False
            smi = molecule_to_smiles(mol)
            print('Could not run filters', smi)
        if filter_pass is False:
            return None

    if rules:
        smi = molecule_to_smiles(mol)
        with open('rules.smi', 'w') as outfile:
            outfile.write('%s\n' %smi)

        result = subprocess.run(['/home/pczbf/Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb rules.smi'], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

        if result == '':
            print('Failed rules', smi)
            return None

    return mol

def is_new_candidate(frag_list, frag_bond_list, candidate_list, candidate_bond_list):

    for i in range(len(candidate_list)):
        if frag_list == candidate_list[i]:
            if frag_bond_list == candidate_bond_list[i]:
                return False

    return True

def combine_all_fragments(frag_mol_list, frag_list, frag_bond_list):

    mol = Molecule()

    new_frag_bond_list = []

    frag_len_list = [len(i.graph.nodes) for i in frag_mol_list]

    added_frag_len_list = [0]

    for i in range(1,len(frag_len_list)):
        added_frag_len_list.append(sum(frag_len_list[:i]))

    for bond in frag_bond_list:
        i = bond[0]
        j = bond[1]
        k = bond[2]
        l = bond[3]

        k += added_frag_len_list[i]
        l += added_frag_len_list[j]

        new_frag_bond_list.append((i,j,k,l))

    graphs = [x.graph for x in frag_mol_list]

    mol.graph = networkx.disjoint_union_all(graphs)

    for bond in new_frag_bond_list:
        k = bond[2]
        l = bond[3]
        mol.graph.add_edge(k, l, order=1)        

    return mol

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file',required=True)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
    parser.add_argument('-x','--parent_fragment_file', help='Parent Fragment Structure File to search fragment database in SDF format',required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0',required=True)
    parser.add_argument('-R','--remove_hydrogens_parent_fragment', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points for the parent fragment in database, numbered from 0',required=True)
    parser.add_argument('-s','--seed', type=int, help='Seed for random number generator',required=False)
    parser.add_argument('-o','--outfile_name', help='Output File Name',required=True)
    parser.add_argument('-n','--n_mol', type=int, help='Number of molecules to generate',required=True)
    parser.add_argument('--unique', action='store_true', help='Generate unique set of molecules', required=False)
    parser.add_argument('--rules', action='store_true', help='Use rules to filter', required=False)

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    if args.unique:
        print('Unique not fully working since does not take symmetry into account')

    build_molecule(fragments_sdf=args.fragments_sdf, fragments_txt=args.fragments_txt, frequencies_txt=args.frequencies_txt, 
        parent_file=args.parent_file, parent_fragment_file=args.parent_fragment_file, remove_hydrogens=args.remove_hydrogens, 
        remove_hydrogens_parent_fragment=args.remove_hydrogens_parent_fragment, outfile_name=args.outfile_name, n_mol=args.n_mol, 
        unique=args.unique, rules=args.rules, filters=True)










