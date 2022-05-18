import glob
import networkx
import sys,os
import random
from pymolgen.molecule_formats import molecule_from_sdf, molecule_to_smiles
from pymolgen.molecule_visualization import plot_molecule, plot_molecule_graph
from pymolgen.molecule import Molecule
from pymolgen.generate import SDFDataset
from pymolgen.bond_generator import RandomBondGenerator

def newmol2(parent_file):
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

    plot_molecule_graph(mol)

    sdf_files = glob.glob('/home/pczbf/pymolgen/datasets/sdf/mol*.sdf') 

    sdf_files_abspath = [ os.path.abspath(sdf_file) for sdf_file in sdf_files]

    dataset = SDFDataset(sdf_files_abspath)

    random.seed(100)

    #generate a random fragment
    frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

    smi = molecule_to_smiles(frag)
    mw = '%.1f' %Molecule.molecular_weight(frag)
    print(smi, mw)

    # Add generated random fragment to mol
    mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator())

    plot_molecule_graph(mol)

    plot_molecule(mol)

    smi = molecule_to_smiles(mol)
    mw = '%.1f' %Molecule.molecular_weight(mol)
    print(smi, mw)

def newmol1():
    min_frag_size = 1
    max_frag_size = 50

    mol = molecule.Molecule()

    mol.load_bru('parent.txt')

    #mol.remove_atom(30)

    print(mol.attach_points)

    smiles = []
    with open("/home/pczbf/pymolgen/datasets/test_set", "r") as smiles_chembl:
        for line in smiles_chembl:
            smiles.append(line)

    dataset = generate.SmilesDataset(smiles)

    #generate a random fragment
    frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

    networkx.draw(frag)
    frag.plot()

    print(frag.attach_points)

    # Add generated random fragment to mol
    mol = Molecule.randomly_glue_together(mol, frag, bond_generator.RandomBondGenerator())

    networkx.draw(mol)
    mol.plot()

if __name__ == '__main__':
    parent_file = sys.argv[1]
    newmol2(parent_file)
