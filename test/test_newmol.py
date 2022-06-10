import sys,os
import random

from pymolgen.newmol import newmol_mw_attachment_points_loop

dir_path = os.path.dirname(os.path.realpath(__file__))

def test_gen():
    newmol_mw_attachment_points_loop(dir_path + "/../datasets/sdf_test_set", "parent.sdf", [29, 30], "test_gen.sdf", 10, seed=100)


if __name__ == "__main__":
    test_gen()
