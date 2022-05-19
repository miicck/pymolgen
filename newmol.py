import glob
import networkx
import sys,os
import random
from pymolgen.molecule_formats import molecule_from_sdf, molecule_to_smiles
from pymolgen.molecule_visualization import plot_molecule, plot_molecule_graph
from pymolgen.molecule import Molecule
from pymolgen.generate import SDFDataset
from pymolgen.bond_generator import RandomBondGenerator

def newmol(parent_file):
    min_frag_size = 1
    max_frag_size = 50

    mol = molecule_from_sdf(parent_file)

    mw = Molecule.molecular_weight(mol)
    print(mw)

    mol = mol.remove_atom(20)
    mol = mol.remove_atom(30)
    #mol = mol.remove_random_atom()

    mw = Molecule.molecular_weight(mol)
    print(mw)

    print(mol.attach_points)

    sdf_files = glob.glob('/home/pczbf/pymolgen/datasets/sdf/mol*.sdf') 

    sdf_files_abspath = [ os.path.abspath(sdf_file) for sdf_file in sdf_files]

    dataset = SDFDataset(sdf_files_abspath)

    #random.seed(100)

    #generate a random fragment
    frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

    smi = molecule_to_smiles(frag)
    mw = '%.1f' %Molecule.molecular_weight(frag)
    print('frag = ', smi, mw)

    # Add generated random fragment to mol
    mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol

    smi = molecule_to_smiles(mol)
    mw = '%.1f' %Molecule.molecular_weight(mol)
    print(smi, mw)

    return mol

def is_hydrogen(graph):
    if len(graph.nodes) != 1:
        return False

    for i in graph.nodes:
        return graph.nodes[i]["element"] == "H"

def newmol_mw_attachment_points(parent_file, remove_hydrogens, max_mw = 500):

    min_frag_size = 1
    max_frag_size = 50

    sdf_files = glob.glob('/home/pczbf/pymolgen/datasets/sdf/mol*.sdf') 

    sdf_files_abspath = [ os.path.abspath(sdf_file) for sdf_file in sdf_files]

    dataset = SDFDataset(sdf_files_abspath)

    mol = molecule_from_sdf(parent_file)

    parent_mw = Molecule.molecular_weight(mol)

    budget_mw = (max_mw - parent_mw) * random.random()

    for i in remove_hydrogens:
        mol = mol.remove_atom(i)

    mw = Molecule.molecular_weight(mol)

    attachment_points = mol.attach_points

#random.seed(100)

    random.shuffle(attachment_points)

    for attachment_point in attachment_points:

        #generate a random fragment
        n = 0
        while True:
            if n == 1000: 
                print('MAX LOOP')
                break
            n += 1
            frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)
            if len(frag.attach_points) > 0 and not is_hydrogen(frag.graph) and Molecule.molecular_weight(frag) + Molecule.molecular_weight(mol) + len(frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw - 14:
                

                smi = molecule_to_smiles(frag)
                mw = '%.1f' %Molecule.molecular_weight(frag)
                #print('frag = ', smi, mw)

                # Add generated random fragment to mol
                mol = Molecule.glue_together_attachmentpoint(mol, frag, RandomBondGenerator(), attachment_point) or mol

                break

    n = 0
    while True:
        n += 1
        if n == 100: 
            print('MAX LOOP, attachment_points =', attachment_points, 'budget_mw =', budget_mw)
            break
            frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)
            if len(frag.attach_points) > 0 and not is_hydrogen(frag.graph) and Molecule.molecular_weight(frag) + Molecule.molecular_weight(mol) + len(frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw - 14:
                smi = molecule_to_smiles(frag)
                mw = '%.1f' %Molecule.molecular_weight(frag)
                mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol
                break

    smi = molecule_to_smiles(mol)
    mw = '%.1f' %Molecule.molecular_weight(mol)
    print('mol = ', smi, mw)
    plot_molecule(mol)

    return mol

if __name__ == '__main__':
    parent_file = sys.argv[1]
    remove_hydrogens = [int(i) for i in sys.argv[2].split()] 
    newmol_mw_attachment_points(parent_file, remove_hydrogens, max_mw = 800)