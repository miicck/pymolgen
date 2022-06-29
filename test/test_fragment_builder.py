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

