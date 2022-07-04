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

def test_filters_additive_mol():

    from openeye import oechem
    from pymolgen.newmol import filters_additive, gen_pains_database, filters_final
    pains_database = gen_pains_database()
    smi = 'COCCN1C(=O)C2(CN(C3C4=C(C=CC(C#N)=C4)OC(C)(C)C3O)C(=O)O2)C2=CC=C(C3=C(C)ON=C3C)C=C21'

    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, smi)

    oechem.OEAddExplicitHydrogens(oemol)

    print(filters_additive(oemol, smi))
    print(filters_final(oemol, smi, pains_database))

def test_filters_final_mol():

    mol = molecule_from_sdf('../datasets/sdf/mol-1.sdf')
    print(mol)
    from pymolgen.newmol import filters_additive, gen_pains_database, filters_final_mol
    pains_database = gen_pains_database()
    print(filters_final_mol(mol, pains_database))


test_filters_final_mol()
