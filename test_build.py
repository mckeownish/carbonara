
# import numpy as np
from glob import glob
import unittest


class TestCarbonara(unittest.TestCase):


    def test_files_written_simple(self):

        test_files = glob('newFitData/test_simple/fitdata/*')

        expected_docs =  [  'newFitData/test_simple/fitdata/fitLog1.dat',
                            'newFitData/test_simple/fitdata/fitLog2.dat',
                            'newFitData/test_simple/fitdata/fitLog3.dat',
                            'newFitData/test_simple/fitdata/mol1_EndScatter.dat',
                            'newFitData/test_simple/fitdata/mol1_InitialScatter.dat',
                            'newFitData/test_simple/fitdata/mol1_xyz_sub_0_step_1.dat',
                            'newFitData/test_simple/fitdata/mol2_EndScatter.dat',
                            'newFitData/test_simple/fitdata/mol2_InitialScatter.dat',
                            'newFitData/test_simple/fitdata/mol2_xyz_sub_0_step_1.dat',
                            'newFitData/test_simple/fitdata/mol3_EndScatter.dat',
                            'newFitData/test_simple/fitdata/mol3_InitialScatter.dat',
                            'newFitData/test_simple/fitdata/mol3_xyz_sub_0_step_1.dat']

        self.assertTrue(set(expected_docs).issubset(test_files), f"Files missing: {set(expected_docs) - set(test_files)}")


    # def test_initial_fit_simple(self):

    #     log_files = [ 'newFitData/test_simple/fitdata/fitLog1.dat',
    #                    'newFitData/test_simple/fitdata/fitLog2.dat',
    #                    'newFitData/test_simple/fitdata/fitLog3.dat' ]

    #     for file_name in log_files:

    #         log_name = file_name.split('/')[-1]
    #         f = open(file_name, "r")
    #         lst = f.readlines()

    #         self.assertTrue( len(lst) > 8, f"Number of lines in {log_name} expected greater than 8, got: {len(lst)}" )

    #         initial_fit = lst[7].split(' ')[2]

    #         diff = float(initial_fit) - float('0.00191271')

    #         with self.subTest(line=file_name):
    #             self.assertTrue( abs( diff ) < 1e-6, f"Initial fit in {log_name} not as expected, difference of {diff} in fit found!" ) 

    #         f.close()

    # def test_initial_scatter_simple(self):

    #     scatter_files = [ 'newFitData/test_simple/fitdata/mol1origScatter.dat',
    #                       'newFitData/test_simple/fitdata/mol2origScatter.dat',
    #                       'newFitData/test_simple/fitdata/mol3origScatter.dat' ]

    #     og_scatter = np.array([ [2.16667e-02, 2.68781e+04, 2.48196e+00, 2.45531e+00],
    #                             [4.50000e-02, 2.38818e+04, 2.36376e+00, 2.32630e+00],
    #                             [6.83333e-02, 1.95240e+04, 2.16229e+00, 2.19422e+00],
    #                             [9.16667e-02, 1.46964e+04, 1.87825e+00, 1.91043e+00],
    #                             [1.15000e-01, 1.02486e+04, 1.51778e+00, 1.55350e+00],
    #                             [1.38333e-01, 6.75520e+03, 1.10096e+00, 1.09895e+00] ])
        
    #     for file_name in scatter_files:

    #         scatter_name = file_name.split('/')[-1]

    #         diff = np.genfromtxt(file_name, skip_footer=1) - og_scatter

    #         with self.subTest(line=file_name):
    #             self.assertTrue( diff.sum() < 1e-6, f"Initial scattering in {scatter_name} not as expected")


if __name__ == '__main__':
    unittest.main()
