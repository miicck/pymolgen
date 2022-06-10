from pymolgen.newmol import newmol_mw_attachment_points_loop


def test_gen():
    newmol_mw_attachment_points_loop("../datasets/sdf_test_set", "parent.sdf", [29, 30], "test_gen.sdf", 10)


if __name__ == "__main__":
    test_gen()
