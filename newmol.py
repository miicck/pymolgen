import glob
import networkx
import sys,os
import random
from pymolgen.molecule_formats import molecule_from_sdf, molecule_to_smiles, molecule_to_atoms_bonds, molecule_to_sdf
from pymolgen.molecule_visualization import plot_molecule, plot_molecule_graph
from pymolgen.molecule import Molecule
from pymolgen.generate import SDFDataset
from pymolgen.bond_generator import RandomBondGenerator
from auto_docker import PAINS_filter, nonallowed_fragment_tautomers

# Import Openeye Modules
from openeye import oechem, oeomega, oedocking, oequacpac

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import AllChem

# Import properties modules
from properties_pymolgen import oeMolProp, num_atomatic_rings, num_chiral_centres, \
    num_lipinsky_donors, num_lipinsky_acceptors, molecular_weight, num_rot_bond

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


try:
    from line_profiler import LineProfiler

    def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner

except ImportError:
    def do_profile(follow=[]):
        "Helpful if you accidentally leave in production!"
        def inner(func):
            def nothing(*args, **kwargs):
                return func(*args, **kwargs)
            return nothing
        return inner

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

def is_hydrogen(molecule):
    if len(molecule.graph.nodes) != 1:
        return False

    for i in molecule.graph.nodes:
        return molecule.graph.nodes[i]["element"] == "H"

def newmol_mw_attachment_points(dataset_path, parent_file, remove_hydrogens, outfile_name, max_mw = 500):

    min_frag_size = 1
    max_frag_size = 50

    sdf_files = glob.glob('%s/mol*.sdf' %dataset_path) 

    sdf_files_abspath = [ os.path.abspath(sdf_file) for sdf_file in sdf_files]

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

        #generate a random fragment
        n = 0
        while True:
            if n == 100: 
                print('MAX LOOP when attaching to attachment_point =', attachment_point)
                break
            n += 1
            frag = dataset.random_molecule().random_fragment(min_size=min_frag_size, max_size=max_frag_size)
            if len(frag.attach_points) > 0 and Molecule.molecular_weight(frag) + Molecule.molecular_weight(mol) + len(frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw:
                smi = molecule_to_smiles(frag)
                plot_molecule(frag)
                plot_molecule_graph(frag)
                mw = '%.1f' %Molecule.molecular_weight(frag)
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
        if not is_hydrogen(frag) and len(frag.attach_points) > 0 and Molecule.molecular_weight(frag) + Molecule.molecular_weight(mol) + len(frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw:
            smi = molecule_to_smiles(frag)
            mw = '%.1f' %Molecule.molecular_weight(frag)
            mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol
            plot_molecule(frag)
            plot_molecule_graph(frag)
            plot_molecule(mol)
            plot_molecule_graph(mol)
            
            print(molecule_to_smiles(mol))
        #break if already at the budget minus 10 (minus 10 because there are no non-hydrogen fragments with mass < 10)
        if Molecule.molecular_weight(mol) + len(mol.attach_points) >= parent_mw + budget_mw - 10:    
            break
        if len(mol.attach_points) == 0:
            break

    smi = molecule_to_smiles(mol)
    mw = Molecule.molecular_weight(mol)
    left_budget = budget_mw + parent_mw - mw
    print('NEW_CANDIDATE %s %.1f %.1f' %(smi, mw, left_budget))
    molecule_to_sdf(mol, outfile_name)

    return mol

def newmol_mw_attachment_points_single(dataset, parent_mol, remove_hydrogens, budget_mw):

    min_frag_size = 1
    max_frag_size = 50

    mol = parent_mol

    parent_mw = Molecule.molecular_weight(mol)

    #budget_mw = (max_mw - parent_mw) * random.random()
    
    for i in remove_hydrogens:
        mol = mol.remove_atom(i)

    protected_atoms = []

    for i in mol.graph:
        protected_atoms.append(i)

    mw = Molecule.molecular_weight(mol)

    attachment_points = mol.attach_points

    random.shuffle(attachment_points)
    frag_counter = 0
    for attachment_point in attachment_points:

        #generate a random fragment
        n = 0
        while True:
            if n == 1000: 
                print('MAX LOOP when attaching to attachment_point =', attachment_point)
                break
            n += 1
            current_budget = parent_mw + budget_mw - Molecule.molecular_weight(mol) - len(mol.attach_points)
            frag = dataset.random_molecule().random_fragment_keep_cycle_max_mass(min_mass=1.0, max_mass=current_budget, counter=frag_counter)
            frag_counter += 1
            if len(frag.attach_points) > 0 and frag.molecular_weight() + mol.molecular_weight() + len(frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw:
                smi = molecule_to_smiles(frag)
                mw = '%.1f' %Molecule.molecular_weight(frag)
                mol = Molecule.glue_together_attachmentpoint(mol, frag, RandomBondGenerator(), attachment_point) or mol
                break

    n = 0
    while True:
        n += 1
        if n == 1000:
            smi = molecule_to_smiles(mol)         
            print('MAX LOOP, attachment_points =', mol.attach_points, 'budget_mw =', budget_mw, smi)
            break
        frag = dataset.random_molecule().random_fragment_keep_cycle_max_mass(min_mass=1.0, max_mass=current_budget, counter=frag_counter)
        frag_counter += 1
        if not is_hydrogen(frag) and len(frag.attach_points) > 0 and frag.molecular_weight() + mol.molecular_weight() + len(frag.attach_points) + len(mol.attach_points) - 2 <= parent_mw + budget_mw:
            smi = molecule_to_smiles(frag)
            mw = '%.1f' %frag.molecular_weight()
            mol = Molecule.randomly_glue_together(mol, frag, RandomBondGenerator()) or mol
        #break if already at the budget minus 10 (minus 10 because there are no non-hydrogen fragments with mass < 10)
        if Molecule.molecular_weight(mol) + len(mol.attach_points) >= parent_mw + budget_mw - 10:    
            break
        if len(mol.attach_points) == 0:
            remove_hydrogen_pass, mol = remove_hydrogen(mol, protected_atoms)
            if remove_hydrogen_pass == False:
                print("Cannot Generate new attachment point")
                break

    lines = molecule_to_sdf(mol)

    with open('pepe.sdf', 'w') as outfile:
        for line in lines:
            outfile.write(line)

    smi = molecule_to_smiles(mol)
    mw = mol.molecular_weight()
    left_budget = budget_mw + parent_mw - mw

    filter_pass = filters_final_mol(mol)
    if filter_pass == False:
        return None

    print('NEW_CANDIDATE %s %.1f %.1f' %(smi, mw, left_budget))

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

def filters_additive(oemol,smi):
    #filters:
    if PAINS_filter(oemol) == False: 
        print("Failed PAINS_filter", smi)
        return False

    n_chiral = num_chiral_centres(oemol)
    h_don = num_lipinsky_donors(oemol)
    h_acc = num_lipinsky_acceptors(oemol)
    n_rot_bonds = num_rot_bond(oemol)

    if n_chiral > 2:
        print("Failed n_chiral filter", smi)
        return False

    if h_don > H_DON_TRESHOLD:
        print("Failed h_don filter", smi)
        return False

    if h_acc > H_ACC_TRESHOLD:
        print("Failed h_acc filter", smi)
        return False

    if n_rot_bonds > ROTBOND_THRESHOLD:
        print("Failed n_rot_bonds filter", smi)
        return False   

    return True

def filters_additive_mol(mol):
    smi = molecule_to_smiles(mol)

    #generate openeye molecule and run filters on it
    try:
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smi)

        oechem.OEAddExplicitHydrogens(oemol)

        return filters_additive(oemol, smi) 

    except:
        print("filters_additive_mol failed with ", smi)
        return False

def filters_final(oemol, smi):

    logp, PSA = oeMolProp(oemol)
    n_aromatic_rings = num_atomatic_rings(smi)
    PFI = n_aromatic_rings + logp

    if logp > LOGP_TRESHOLD_UP or logp < LOGP_TRESHOLD_LOW:
        print("Failed logP filter", smi)
        return False

    if PSA > PSA_TRESHOLD:
        print("Failed PSA filter", smi)
        return False

    if PFI > PFI_TRESHOLD:
        print("Failed PFI filter", smi)
        return False

    return True    

def filters_final_mol(mol):
    smi = molecule_to_smiles(mol)

    #generate openeye molecule and run filters on it
    try:
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smi)

        oechem.OEAddExplicitHydrogens(oemol)

        filters_additive_pass = filters_additive(oemol, smi)
        if filters_additive_pass == False:
            return False

        filters_final_pass = filters_final(oemol, smi) 
        if filters_final_pass == False:
            return False

        return True

    except Exception as e: 
        print("filters_final_mol failed with ", smi)
        print(e)
        return False

#@do_profile(follow=[newmol_mw_attachment_points_single, Molecule.random_fragment_keep_cycle_max_mass])
def newmol_mw_attachment_points_loop(dataset_path, parent_file, remove_hydrogens, outfile_name, n_mol, max_mw = 500):

    sdf_files = glob.glob('%s/mol*.sdf' %dataset_path) 

    sdf_files_abspath = [ os.path.abspath(sdf_file) for sdf_file in sdf_files]

    dataset = SDFDataset(sdf_files_abspath)

    with open(outfile_name, "w") as outfile:
        print("Generating ", outfile_name)

    parent_mol = molecule_from_sdf(parent_file)

    parent_mw = Molecule.molecular_weight(parent_mol)

    budget_mw = (max_mw - parent_mw) * random.random()

    counter = 0
    while counter < n_mol:
        mol = newmol_mw_attachment_points_single(dataset, parent_mol, remove_hydrogens, budget_mw)
        if mol is None: continue

        counter += 1
        budget_mw = (max_mw - parent_mw) * random.random()

        lines = molecule_to_sdf(mol)

        with open(outfile_name, "a") as outfile:
            for line in lines:
                outfile.write(line)
            outfile.write('$$$$\n')


if __name__ == '__main__':
    dataset_path = sys.argv[1]
    parent_file = sys.argv[2]
    remove_hydrogens = [int(i) for i in sys.argv[3].split()]
    outfile_name = sys.argv[4]
    n_mol = int(sys.argv[5])

    random.seed(100)
    print("Random = ", random.random())
    


    newmol_mw_attachment_points_loop(dataset_path, parent_file, remove_hydrogens, outfile_name, n_mol)