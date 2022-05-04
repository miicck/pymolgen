import random

from pymolgen.molecule import Molecule, FractionalOrderException


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
    assert m1.hydrogenate() == 1
    assert str(m1) == "C"


def test_hydrogenete_2():
    m = Molecule().load_smiles("[C]:[C]OC")
    assert m.hydrogenate() == 5


def test_hydrogenete_4():
    m = Molecule().load_smiles("[O]c:1:c:c2:c([C]3[CH][CH][C](C([CH][C]3C(CC2)NC(=O)C)=O)[Br]):c(:c1[O])OC")
    assert m.hydrogenate() == 8


def test_hydrogenete_5():
    # Before hydrogenation
    m1 = Molecule().load_smiles(
        "O=C1CN2[C](CSc:3:n:n:c(:[o]3)c:3:c:[c]:4:[c](:c(:n3)c:3:c:c:c:c(:c3)[N+]([O-])=O)"
        ":[nH]:[c]:3:c:c:c:c:[c]43)[C]N3CC(NC(C3C2C(N1)=O)=O)=O")

    assert m1.valid_smiles
    m2 = m1.copy()
    m2.hydrogenate()
    assert m2.valid_smiles


def test_ch4_from_bits():
    m1 = Molecule().load_smiles("[H]")
    m2 = Molecule().load_smiles("[CH]")
    m3 = Molecule.randomly_glue_together(m1, m2)
    m3.hydrogenate()
    assert str(m3) == "C"


def test_frac_valence():
    m1 = Molecule(allow_frac_order=True).load_smiles(
        "[N][C]([C]NC1CCCCC1)[C]")
    m2 = Molecule(allow_frac_order=True).load_smiles(
        "Cc:1:c:c(:c:c(:c1Oc:1:n:c(:n:[c]:2:c:c:[s]:[c]21)NC1CCN(C1)Cc:1:c:c:n:c:c1)[CH])[C]")
    m3 = Molecule.randomly_glue_together(m1, m2)
    m3.hydrogenate()
