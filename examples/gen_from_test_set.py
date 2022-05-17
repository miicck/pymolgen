from pymolgen.molecule import Molecule
from pymolgen.generate import generate_from_fragments
from pymolgen.molecule_visualization import plot_molecule, plot_molecule_graph
from pymolgen.bond_generator import RandomBondGenerator, MaxOrderBondGenerator
import matplotlib.pyplot as plt
import random

with open("../datasets/test_set", "r") as test_set:
    smiles = list(line for line in test_set)

random.seed(100)

for mol, fragments in generate_from_fragments(
        smiles, lambda m: m.atom_count > 50, MaxOrderBondGenerator(), max_frag_size=20):
    try:
        plot_molecule(mol, title=f"Molecule from {fragments} fragments")
        plot_molecule_graph(mol)
    except Exception as e:
        print(f"Failed: {e}")
