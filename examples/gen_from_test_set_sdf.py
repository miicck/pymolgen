from pymolgen.molecule import Molecule
from pymolgen.generate import *
from pymolgen.molecule_visualization import plot_molecule, plot_molecule_graph
from pymolgen.bond_generator import RandomBondGenerator, MaxOrderBondGenerator
from molecule_formats import *
import matplotlib.pyplot as plt
import random
import glob
import os

sdf_files = glob.glob('../datasets/sdf/mol*.sdf') 

sdf_files_abspath = [ os.path.abspath(sdf_file) for sdf_file in sdf_files]

dataset = SDFDataset(sdf_files_abspath)

random.seed(100)

for mol, fragments in generate_from_molecules(
        dataset, lambda m: m.atom_count > 50, MaxOrderBondGenerator(), max_frag_size=20):
    smi = molecule_to_smiles(mol)
    print(smi)
    try:
        plot_molecule(mol, title=f"Molecule from {fragments} fragments")
        #plot_molecule_graph(mol)
    except Exception as e:
        print(f"Failed: {e}")
