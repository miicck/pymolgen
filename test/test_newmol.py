import sys,os
import random
import filecmp

from pymolgen.newmol import newmol_mw_attachment_points_loop, newmol_mw_attachment_points_loop_large

dir_path = os.path.dirname(os.path.realpath(__file__))

#def test_gen():
    #newmol_mw_attachment_points_loop(dir_path + "/../datasets/sdf_test_set", "parent.sdf", [29, 30], "test_gen.sdf", 10, seed=100)
    #newmol_mw_attachment_points_loop_large('/home/pczbf/Downloads/chembl_30.sdf', 'parent.sdf', [29,30], 'test_gen.sdf', 100, max_n = 100, seed=100)

def test_phenylisoxazole():
    newmol_mw_attachment_points_loop_large(dataset_file='../datasets/database1000/database1000.sdf', 
        parent_file='../datasets/phenylisoxazole/phenylisoxazole.sdf',
        outfile_name='test_gen.sdf', n_mol=100, remove_hydrogens=[17,20], 
        remove_hydrogens_max=[12,13,14,22,23], remove_hydrogens_max_n=2, 
        seed=100, logfilename='test_logfile.txt')

    compare = filecmp.cmp('test_logfile.txt', '../datasets/phenylisoxazole/phenylisoxazole.txt', shallow=False)
    assert compare is True

    compare = filecmp.cmp('test_gen.sdf', '../datasets/phenylisoxazole/pymolgen0.sdf', shallow=False)

if __name__ == '__main__':
    test_gen()
