#!/usr/bin/python3

import glob
import networkx
import sys, os
import random
import argparse
from pymolgen.molecule_formats import molecule_from_sdf, molecule_to_smiles, molecule_to_atoms_bonds, molecule_to_sdf
from pymolgen.molecule_visualization import plot_molecule, plot_molecule_graph
from pymolgen.molecule import Molecule
from pymolgen.generate import SDFDataset, SDFDatasetLarge, SDFDatasetLargeRAM
from pymolgen.bond_generator import RandomBondGenerator

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import AllChem

from openeye import oechem

from pymolgen.properties_pymolgen import oeMolProp, num_atomatic_rings, num_chiral_centres, \
            num_lipinsky_donors, num_lipinsky_acceptors, molecular_weight, num_rot_bond, pains_filter, gen_pains_database

from auto_docker import PAINS_filter

from functools import partial
print = partial(print, flush=True)

# =============================================================================
# Define Thresholds
# =============================================================================

#   CHIRAL_THRESHOLD:   Maximum number of Chiral Centers (<= 2)
CHIRAL_THRESHOLD = 2
#   PSA_TRESHOLD:       Polar surface area (<= 140)
PSA_TRESHOLD = 140
#   PFI_TRESHOLD:       Property Forecast Index (LOGP + # of aromatic rings < 8)
PFI_TRESHOLD = 8
#   ROTBOND_THRESHOLD:  Number of rotatable bonds (<= 7)
ROTBOND_THRESHOLD = 7
#   WEIGHT_TRESHOLD:    Maximum MW in Daltons (<= 500)
WEIGHT_TRESHOLD = 500
#   H_DON_TRESHOLD:     Maximum number of hydrogen donors (<= 5)
H_DON_TRESHOLD = 5
#   H_ACC_TRESHOLD:     Maximum number of hydrogen acceptors (<= 10)
H_ACC_TRESHOLD = 10
#   LOGP_TRESHOLD:      Water/Octanol Partition Coefficient (0.5-5.0)
LOGP_TRESHOLD_UP = 5
LOGP_TRESHOLD_LOW = 0.5


def newmol(parent_file):
    min_frag_size = 1
    max_frag_size = 50

    mol = molecule_from_sdf(parent_file)

    mw = Molecule.molecular_weight(mol)
    print(mw)

    mol = mol.remove_atom(20)
    mol = mol.remove_atom(30)
    # mol = mol.remove_random_atom()

    mw = Molecule.molecular_weight(mol)
    print(mw)

    print(mol.attach_points)

    sdf_files = glob.glob('/home/pczbf/pymolgen/datasets/sdf/mol*.sdf')

    sdf_files_abspath = [os.path.abspath(sdf_file) for sdf_file in sdf_files]

    dataset = SDFDataset(sdf_files_abspath)

    # random.seed(100)

    # generate a random fragment
    frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)

    smi = molecule_to_smiles(frag)
    mw = '%.1f' % Molecule.molecular_weight(frag)
    print('frag = ', smi, mw)

    # Add generated random fragment to mol
    mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol

    smi = molecule_to_smiles(mol)
    mw = '%.1f' % Molecule.molecular_weight(mol)
    print(smi, mw)

    return mol


def is_hydrogen(molecule):
    if len(molecule.graph.nodes) != 1:
        return False

    for i in molecule.graph.nodes:
        return molecule.graph.nodes[i]["element"] == "H"


def newmol_mw_attachment_points(dataset_path, parent_file, remove_hydrogens, outfile_name, max_mw=500):
    min_frag_size = 1
    max_frag_size = 50

    sdf_files = glob.glob('%s/mol*.sdf' % dataset_path)

    sdf_files_abspath = [os.path.abspath(sdf_file) for sdf_file in sdf_files]

    dataset = SDFDataset(sdf_files_abspath)

    mol = molecule_from_sdf(parent_file)

    parent_mw = Molecule.molecular_weight(mol)

    budget_mw = (max_mw - parent_mw) * random.random()

    for i in remove_hydrogens:
        mol = mol.remove_atom(i)

    mw = Molecule.molecular_weight(mol)

    attachment_points = mol.attach_points

    random.seed(100)

    random.shuffle(attachment_points)

    for attachment_point in attachment_points:

        # generate a random fragment
        n = 0
        while True:
            if n == 100:
                print('MAX LOOP when attaching to attachment_point =', attachment_point)
                break
            n += 1
            frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)
            if len(frag.attach_points) > 0 and Molecule.molecular_weight(frag) + Molecule.molecular_weight(mol) + len(
                    frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw:
                smi = molecule_to_smiles(frag)
                plot_molecule(frag)
                plot_molecule_graph(frag)
                mw = '%.1f' % Molecule.molecular_weight(frag)
                mol = Molecule.glue_together_attachmentpoint(mol, frag, RandomBondGenerator(), attachment_point) or mol
                plot_molecule(mol)
                plot_molecule_graph(mol)

                break

    n = 0
    while True:
        n += 1
        if n == 100:
            print('MAX LOOP, attachment_points =', mol.attach_points, 'budget_mw =', budget_mw)
            break
        frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)
        if not is_hydrogen(frag) and len(frag.attach_points) > 0 and Molecule.molecular_weight(
                frag) + Molecule.molecular_weight(mol) + len(frag.attach_points) + len(
            mol.attach_points) - 2 <= parent_mw + budget_mw:
            smi = molecule_to_smiles(frag)
            mw = '%.1f' % Molecule.molecular_weight(frag)
            mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol
            plot_molecule(frag)
            plot_molecule_graph(frag)
            plot_molecule(mol)
            plot_molecule_graph(mol)

            print(molecule_to_smiles(mol))
        # break if already at the budget minus 10 (minus 10 because there are no non-hydrogen fragments with mass < 10)
        if Molecule.molecular_weight(mol) + len(mol.attach_points) >= parent_mw + budget_mw - 10:
            break
        if len(mol.attach_points) == 0:
            break

    smi = molecule_to_smiles(mol)
    mw = Molecule.molecular_weight(mol)
    left_budget = budget_mw + parent_mw - mw
    print('NEW_CANDIDATE %s %.1f %.1f' % (smi, mw, left_budget))
    molecule_to_sdf(mol, outfile_name)

    return mol


def newmol_mw_attachment_points_single(dataset, parent_mol, budget_mw, pains_database, remove_hydrogens=None, 
    remove_hydrogens_max=None, remove_hydrogens_max_n=None):
    
    min_frag_size = 1
    max_frag_size = 50

    mol = parent_mol

    parent_mw = mol.molecular_weight()

    #create attachment points from remove_hydrogens list
    if remove_hydrogens is not None:
        for i in remove_hydrogens:
            mol = mol.remove_atom(i)

    #create up to remove_hydrogens_max_n attachment points from remove_hydrogens_max list
    if remove_hydrogens_max is not None and remove_hydrogens_max_n is not None:
        random.shuffle(remove_hydrogens_max)
        remove_hydrogens_max = remove_hydrogens_max[:remove_hydrogens_max_n]
        for i in remove_hydrogens_max:
            mol = mol.remove_atom(i)

    protected_atoms = []

    for i in mol.graph:
        protected_atoms.append(i)

    mw = mol.molecular_weight()

    attachment_points = mol.attach_points
    if len(attachment_points) == 0:
        raise Exception("No attachment points in parent structure!")

    random.shuffle(attachment_points)
    frag_counter = 0
    n_rot_bonds = 0
    for attachment_point in attachment_points:

        # generate a random fragment
        n = 0
        if n_rot_bonds == ROTBOND_THRESHOLD:
            break
        while True:
            if n == 1000:
                print('MAX LOOP when attaching to attachment_point =', attachment_point)
                break
            n += 1
            mol_n_attach_points = len(mol.attach_points)
            current_budget = parent_mw + budget_mw - mol.molecular_weight() - mol_n_attach_points
            frag = dataset.random_molecule().random_fragment_keep_cycle_max_mass(min_mass=1.0, max_mass=current_budget,
                                                                                 counter=frag_counter)
            frag_counter += 1
            frag_n_attach_points = len(frag.attach_points)

            if frag_n_attach_points > 0 and frag.molecular_weight() + mol.molecular_weight() + frag_n_attach_points + mol_n_attach_points - 2 <= parent_mw + budget_mw:
                smi = molecule_to_smiles(frag)
                mw = '%.1f' % frag.molecular_weight()
                # print("budget = %.1f" %current_budget, "frag1 = ", smi)
                newmol = Molecule.glue_together_attachmentpoint(mol, frag, RandomBondGenerator(),
                                                                attachment_point) or mol
                filters_additive_pass, n_rot_bonds = filters_additive_mol(newmol)
                if filters_additive_pass:
                    mol = newmol
                    break

    n = 0
    while n_rot_bonds < ROTBOND_THRESHOLD:
        n += 1
        if n == 1000:
            smi = molecule_to_smiles(mol)
            print('MAX LOOP, attachment_points =', mol.attach_points, 'budget_mw =', budget_mw, smi)
            break
        frag = dataset.random_molecule().random_fragment_keep_cycle_max_mass(min_mass=1.0, max_mass=current_budget,
                                                                             counter=frag_counter)
        frag_counter += 1
        mol_n_attach_points = len(mol.attach_points)
        frag_n_attach_points = len(frag.attach_points)
        current_budget = parent_mw + budget_mw - mol.molecular_weight() - mol_n_attach_points
        if not is_hydrogen(frag) and frag_n_attach_points > 0 and frag.molecular_weight() + mol.molecular_weight() + frag_n_attach_points + mol_n_attach_points - 2 <= parent_mw + budget_mw:
            #smi = molecule_to_smiles(frag)
            #mw = '%.1f' % frag.molecular_weight()
            # print("budget = %.1f" %current_budget, "frag2 = ", smi)
            newmol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol
            filters_additive_pass, n_rot_bonds = filters_additive_mol(newmol)
            if filters_additive_pass:
                mol = newmol
                if n_rot_bonds == ROTBOND_THRESHOLD:
                    break
        # break if already at the budget minus 10 (minus 10 because there are no non-hydrogen fragments with mass < 10)
        mol_n_attach_points = len(mol.attach_points)
        if Molecule.molecular_weight(mol) + mol_n_attach_points >= parent_mw + budget_mw - 10:
            break
        if mol_n_attach_points == 0:
            remove_hydrogen_pass, mol = remove_hydrogen(mol, protected_atoms)
            if remove_hydrogen_pass == False:
                print("Cannot Generate new attachment point")
                break

    lines = molecule_to_sdf(mol)

    mw = mol.molecular_weight()
    left_budget = budget_mw + parent_mw - mw

    filter_pass = filters_final_mol(mol, pains_database)
    if filter_pass == False:
        return None

    #smi = molecule_to_smiles(mol)
    #print('NEW_CANDIDATE %s %.1f %.1f' % (smi, mw, left_budget))

    return mol


def remove_hydrogen(mol, protected_atoms):
    """
    Removes one hydrogen from molecule to create new attachment point at random
    without removing from protected atoms (original hydrogens in parent structure)
    """

    hydrogens = []
    for i in mol.graph:
        element = mol.graph.nodes[i]['element']
        if element == 'H' and i not in protected_atoms:
            hydrogens.append(i)

    if len(hydrogens) > 0:
        hydrogen = random.choice(hydrogens)
        mol = mol.remove_atom(hydrogen)
        return (True, mol)

    return (False, mol)


def filters_additive(oemol, smi):

    n_rot_bonds = num_rot_bond(oemol)

    if n_rot_bonds > ROTBOND_THRESHOLD:
        #print("Failed n_rot_bonds filter", smi)
        return (False, n_rot_bonds)

    n_chiral = num_chiral_centres(oemol)

    if n_chiral > 2:
        #print("Failed n_chiral filter", smi)
        return (False, n_rot_bonds)

    h_don = num_lipinsky_donors(oemol)

    if h_don > H_DON_TRESHOLD:
        #print("Failed h_don filter", smi)
        return (False, n_rot_bonds)

    h_acc = num_lipinsky_acceptors(oemol)

    if h_acc > H_ACC_TRESHOLD:
        #print("Failed h_acc filter", smi)
        return (False, n_rot_bonds)

    return (True, n_rot_bonds)


def filters_additive_mol(mol):
    """
    Runs filters for a molecule for quantities that
    are additive as the molecule increases in size.

    Parameters
    ----------
    mol:
        Molecule to evaluate filters for

    Returns
    -------
    passed: bool
        If the filters pass or not
    n_rot_bonds: int
        The number of rotatatable bonds in the molecule (if pass)
    """

    smi = molecule_to_smiles(mol)

    # generate openeye molecule and run filters on it
    try:
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smi)

        oechem.OEAddExplicitHydrogens(oemol)

        filters_additive_pass, n_rot_bonds = filters_additive(oemol, smi)
        #print("filters_additive_pass =", filters_additive_pass, "n_rot_bonds = ", n_rot_bonds)
        return (filters_additive_pass, n_rot_bonds)

    except:
        raise Exception("filters_additive_mol failed with ", smi)
        return (False, 0)


def filters_final(oemol, smi, pains_database):

    logp, PSA = oeMolProp(oemol)

    if logp > LOGP_TRESHOLD_UP or logp < LOGP_TRESHOLD_LOW:
        print("Failed logP filter", smi)
        return False

    if PSA > PSA_TRESHOLD:
        print("Failed PSA filter", smi)
        return False

    n_aromatic_rings = num_atomatic_rings(smi)

    PFI = n_aromatic_rings + logp

    if PFI > PFI_TRESHOLD:
        print("Failed PFI filter", smi)
        return False

    if pains_filter(oemol, pains_database) == False:
        print("Failed PAINS_filter", smi)
        return False

    return True


def filters_final_mol(mol, pains_database):

    smi = molecule_to_smiles(mol)

    # generate openeye molecule and run filters on it
    try:
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smi)

        oechem.OEAddExplicitHydrogens(oemol)

        #removed here because test carried out as you build molecule so it'll never fail
        #filters_additive_pass, n_not_bonds = filters_additive(oemol, smi, pains_database)
        #if filters_additive_pass == False:
        #    print("Failed filters_additive")
        #    return False

        filters_final_pass = filters_final(oemol, smi, pains_database)
        if filters_final_pass == False:
            return False

        return True

    except Exception as e:
        print("filters_final_mol failed with ", smi)
        print(e)
        return False


def newmol_mw_attachment_points_loop(dataset_path, parent_file, remove_hydrogens, outfile_name, n_mol, max_mw=500,
                                     seed=None):
    if seed is not None:
        random.seed(seed)

    sdf_files = glob.glob('%s/mol*.sdf' % dataset_path)
    if len(sdf_files) == 0:
        raise Exception("No molecules in dataset!")

    sdf_files_abspath = [os.path.abspath(sdf_file) for sdf_file in sdf_files]
    dataset = SDFDataset(sdf_files_abspath)

    with open(outfile_name, "w") as outfile:
        print("Generating ", outfile_name)

    parent_mol = molecule_from_sdf(parent_file)
    parent_mw = Molecule.molecular_weight(parent_mol)
    budget_mw = (max_mw - parent_mw) * random.random()

    try:
        pains_database = gen_pains_database()
    except:
        raise Exception("Could not generate pains database")

    counter = 0
    while counter < n_mol:
        mol = newmol_mw_attachment_points_single(dataset, parent_mol, remove_hydrogens, budget_mw, pains_database)
        if mol is None: continue

        counter += 1
        budget_mw = (max_mw - parent_mw) * random.random()
        lines = molecule_to_sdf(mol)

        with open(outfile_name, "a") as outfile:
            for line in lines:
                outfile.write(line)
            outfile.write('$$$$\n')


def count_generated_molecules(outfile_name):
    """
    Count the number of generated molecules in an SDF file
    """
    generated_molecules = 0
    with open(outfile_name) as f:
        for line in f:
            if 'V2000' in line:
                generated_molecules += 1

def newmol_mw_attachment_points_loop_large(dataset_file, parent_file, outfile_name, n_mol, remove_hydrogens=None, max_n = None, 
    max_mw=500, remove_hydrogens_max=None, remove_hydrogens_max_n=None, seed=None, logfilename=None, restart=False):
    
    if seed is not None:
        random.seed(seed)

    print('max_n =', max_n)
    dataset = SDFDatasetLargeRAM(dataset_file, max_n)

    print(sys.getsizeof(dataset.lines))

    if restart is False:
        with open(outfile_name, "w") as outfile:
            print("Generating ", outfile_name)

    parent_mol = molecule_from_sdf(parent_file)
    parent_mw = Molecule.molecular_weight(parent_mol)
    budget_mw = (max_mw - parent_mw) * random.random()

    try:
        pains_database = gen_pains_database()
    except:
        raise Exception("Could not generate pains database")

    if logfilename is not None and restart is False:
            open(logfilename, 'w').close()

    counter = 0

    if restart is True:
        counter = count_generated_molecules(outfile_name)

    while counter < n_mol:
        mol = newmol_mw_attachment_points_single(dataset, parent_mol, budget_mw, pains_database, 
            remove_hydrogens, remove_hydrogens_max, remove_hydrogens_max_n)
        if mol is None: continue

        counter += 1
        budget_mw = (max_mw - parent_mw) * random.random()
        lines = molecule_to_sdf(mol)

        with open(outfile_name, "a") as outfile:
            for line in lines:
                outfile.write(line)
            outfile.write('$$$$\n')

        smi = molecule_to_smiles(mol)
        printline = 'NEW_CANDIDATE %s %s' % (counter, smi)
        print(printline)
        if logfilename is not None:
            with open(logfilename, 'a') as logfile:
                logfile.write('%s\n' %printline)




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator')
    parser.add_argument('-d','--dataset_path', help='Dataset path from which fragments will be generated',required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
    parser.add_argument('-o','--outfile_name', help='Output File Name',required=True)
    parser.add_argument('-n','--n_mol', type=int, help='Number of molecules to generate',required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0',required=False)
    parser.add_argument('-R','--remove_hydrogens_max', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points up to remove_hydrogens_max_n simultaneously, numbered from 0',required=False)
    parser.add_argument('-m','--remove_hydrogens_max_n', type=int, help='Maximum number of hydrogen atoms that will be created as attachment points from remove_hydrogens_max simulaneouly',required=False)
    parser.add_argument('-s','--seed', type=int, help='Seed for random number generator',required=False)
    parser.add_argument('-l','--logfilename', help='Log file',required=False)
    parser.add_argument('--restart', action='store_true', help='Restart generation from previous run')
    args = parser.parse_args()

    if args.remove_hydrogens_max is not None and args.remove_hydrogens_max_n is None:
        sys.exit('Need to provide remove_hydrogens_max_n if remove_hydrogens_max set')

    if args.remove_hydrogens_max_n is not None and args.remove_hydrogens_max is None:
        sys.exit('Need to provide remove_hydrogens_max if remove_hydrogens_max_n set')

    if args.seed is not None and args.n_mol > 1000:
        sys.exit('Cannot run with seed and n_mol > 1000')

    newmol_mw_attachment_points_loop_large(dataset_file=args.dataset_path, parent_file=args.parent_file, 
        outfile_name=args.outfile_name, n_mol=args.n_mol, remove_hydrogens=args.remove_hydrogens, 
        remove_hydrogens_max=args.remove_hydrogens_max, remove_hydrogens_max_n=args.remove_hydrogens_max_n, 
        seed=args.seed, logfilename=args.logfilename, restart=args.restart)
