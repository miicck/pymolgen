import random

from pymolgen.molecule import Molecule


def test_molecule():
    mol = Molecule().load_smiles("CSc1ncccc1C(=O)OCc1ccccc1F")


def test_cehmbl_molecules():
    with open("../datasets/smiles_chembl", "r") as smiles_chembl:
        for i, smiles_string in enumerate(smiles_chembl):
            m = Molecule().load_smiles(smiles_string)
            if i > 100:
                break


def test_order():
    m = Molecule().load_smiles("[H][H]")
    assert m.total_free_valence == 0


def test_fragments():
    m = Molecule().load_smiles("[H][H]")
    assert m.atom_count == 2

    f = m.random_fragment()
    assert len(list(f.free_valence_points)) == 1
    assert f.total_free_valence == 1


def test_glue_h2():
    # Create two free hydrogen atoms
    m1 = Molecule().load_smiles("[H][H]")
    f1 = m1.random_fragment()
    assert f1.total_free_valence == 1

    m2 = Molecule().load_smiles("[H][H]")
    f2 = m2.random_fragment()
    assert f2.total_free_valence == 1

    # Glue them together
    m3 = Molecule.randomly_glue_together(f1, f2)
    assert m3.atom_count == 2
    assert m3.total_free_valence == 0
    assert str(m3) == "[H][H]"


def test_glue_ch3s():
    # Create a CH3 by removing a hydrogen from methane
    m1 = Molecule().load_smiles("C")
    assert m1.atom_count == 5
    assert m1.remove_random_atom("H")
    assert m1.total_free_valence == 1

    m2 = m1.copy()
    assert m2.total_free_valence == 1

    m3 = Molecule.randomly_glue_together(m1, m2)
    assert m3.total_free_valence == 0
    assert str(m3) == "CC"


def test_hydrogenate():
    m1 = Molecule().load_smiles("C")
    m1.remove_random_atom("H")
    m1.hydrogenate()
    assert str(m1) == "C"


def test_hydrogenete_2():
    m = Molecule().load_smiles("[N]C([CH]Cc:1:c:c(:c(OC):[c]:[c]1)O[CH2])[C][C]")
    assert m.hydrogenate() == 13


def test_hydrogenete_3():
    m = Molecule().load_smiles("[C]:[C]OC")
    assert m.hydrogenate() == 5


def test_hydrogenete_4():
    m = Molecule().load_smiles("[O]c:1:c:c2:c([C]3[CH][CH][C](C([CH][C]3C(CC2)NC(=O)C)=O)[Br]):c(:c1[O])OC")
    assert m.hydrogenate() == 8


def test_gen_mol():
    source_mols = []
    with open("../datasets/smiles_chembl", "r") as smiles_chembl:
        for smiles_string in smiles_chembl:
            source_mols.append(Molecule().load_smiles(smiles_string))
            if len(source_mols) >= 100:
                break

    while True:

        target_atom_count = random.randint(20, 100)

        try:

            traj = []
            mol = Molecule().load_smiles("[H][H]")
            mol.remove_random_atom("H")
            traj.append(mol.copy())

            while mol.atom_count < target_atom_count:
                frag = random.choice(source_mols).random_fragment(max_size=16)
                mol = Molecule.randomly_glue_together(mol, frag)
                traj.append(mol.copy())

            traj[-1].plot()

        except:
            pass
