import sys,os
import random

from pymolgen.newmol import newmol_mw_attachment_points_loop, newmol_mw_attachment_points_loop_large

dir_path = os.path.dirname(os.path.realpath(__file__))

def test_gen():
    #newmol_mw_attachment_points_loop(dir_path + "/../datasets/sdf_test_set", "parent.sdf", [29, 30], "test_gen.sdf", 10, seed=100)
    newmol_mw_attachment_points_loop_large('/home/pczbf/Downloads/chembl_30.sdf', "parent.sdf", [29,30], 'test_gen.sdf', 100, max_n = 10000, seed=100)


if __name__ == "__main__":
    test_gen()
