import sys,os

from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *
from pymolgen.molecule import *
from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *
from rdkit import Chem

def test_update_bond_frequencies():

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
    frag_mapping = get_frag_mapping('../datasets/database1000/fragments1.txt')
    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')   
    
    print(bond_frequencies)

    bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)

    print(bond_frequencies)

    check = {(0, 1, 0, 0): 1, (1, 2, 2, 1): 1, (2, 3, 1, 0): 1, (3, 4, 0, 2): 1, (4, 5, 2, 0): 1, (5, 6, 0, 1): 1, (4, 6, 2, 3): 1, (4, 7, 2, 1): 1, (7, 8, 1, 3): 1, (6, 9, 9, 0): 1, (6, 9, 11, 0): 1}

    assert bond_frequencies == check

def test_free_valence_list():

    fragment_database = get_fragment_database('../datasets/fragments/fragments10.sdf')

    free_valence_list_list = []

    for i in range(len(fragment_database)):

        free_valence_list_list.append(fragment_database[i].free_valence_list)

    check = [[0], [0, 2], [1, 1], [0, 0], [2, 2], [0, 0], [1, 3, 9, 11], [1, 5], [3], [0], [1], [0, 2, 6, 7], [0], [4], [0, 15, 22, 26, 31, 35], [7, 11], [0, 0, 0], [0, 1, 1], [2], [8], [0, 3], [0, 3, 4], [3, 10, 12], [0, 3], [6, 10], [0], [5, 11], [0, 0], [0, 1, 5, 8], [0, 3, 6], [0, 7], [3], [0, 5], [8], [0, 2], [2, 9], [4], [5, 8]]

    assert free_valence_list_list == check

