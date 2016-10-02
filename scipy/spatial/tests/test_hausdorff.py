import numpy as np
from numpy.testing import (TestCase,
                           assert_almost_equal)
import scipy
from scipy.spatial.distance import directed_hausdorff

class TestHausdorff(TestCase):
    '''Test various properties of the directed Hausdorff code.'''
   
    def setUp(self):
        self.random_angles = np.random.random((100,)) * np.pi * 2
        self.random_columns = np.column_stack((self.random_angles,
                                               self.random_angles,
                                               np.zeros((100,))))
        self.random_columns[...,0] = np.cos(self.random_columns[...,0])
        self.random_columns[...,1] = np.sin(self.random_columns[...,1])
        self.random_columns_2 = np.column_stack((self.random_angles,
                                                 self.random_angles,
                                                 np.zeros((100,))))
        self.random_columns_2[1:,0] = np.cos(self.random_columns_2[1:,0]) * 2.0
        self.random_columns_2[1:,1] = np.sin(self.random_columns_2[1:,1]) * 2.0
        # move one point farther out so we don't have two perfect circles
        self.random_columns_2[0,0] = np.cos(self.random_columns_2[0,0]) * 3.3
        self.random_columns_2[0,1] = np.sin(self.random_columns_2[0,1]) * 3.3
        self.path_1 = self.random_columns
        self.path_2 = self.random_columns_2
        self.path_1_4d = np.insert(self.path_1, 3, 5, axis = 1)
        self.path_2_4d = np.insert(self.path_2, 3, 27, axis = 1)

    def tearDown(self):
        del self.random_angles
        del self.random_columns
        del self.random_columns_2
        del self.path_1
        del self.path_2
        del self.path_1_4d
        del self.path_2_4d

    def test_symmetry(self):
        '''Ensure that the directed (asymmetric)
        Hausdorff distance is actually asymmetric
        '''
        forward = directed_hausdorff(self.path_1, self.path_2)
        reverse = directed_hausdorff(self.path_2, self.path_1)
        self.assertNotEqual(forward, reverse)

    def test_brute_force_comparison_forward(self):
        '''Ensure that the algorithm for 
        directed_hausdorff gives the same
        result as the simple / brute force
        approach in the forward direction.'''
        actual = directed_hausdorff(self.path_1, self.path_2)   
        # brute force over rows:
        expected = max(np.amin(scipy.spatial.distance.cdist(self.path_1,
                                                            self.path_2), 
                                                            axis = 1))
        assert_almost_equal(actual, expected, decimal=9)

    def test_brute_force_comparison_reverse(self):
        '''Ensure that the algorithm for 
        directed_hausdorff gives the same
        result as the simple / brute force
        approach in the reverse direction.'''
        actual = directed_hausdorff(self.path_2, self.path_1)
        # brute force over columns:
        expected = max(np.amin(scipy.spatial.distance.cdist(self.path_1,
                                                            self.path_2), 
                                                            axis = 0))
        assert_almost_equal(actual, expected, decimal=9)

    def test_degenerate_case(self):
        '''The directed Hausdorff distance must be zero
        if both input data arrays match.'''
        actual = directed_hausdorff(self.path_1, self.path_1)
        assert_almost_equal(actual, 0.0, decimal=9)

    def test_2d_data_forward(self):
        '''Ensure that 2D data is handled
        properly for a simple case relative
        to brute force approach.'''
        actual = directed_hausdorff(self.path_1[...,:2],
                                    self.path_2[...,:2])
        expected = max(np.amin(scipy.spatial.distance.cdist(
                                                      self.path_1[...,:2],
                                                      self.path_2[...,:2]), 
                                                      axis = 1))
        assert_almost_equal(actual, expected, decimal=9)

    def test_4d_data_reverse(self):
        '''Ensure that 4D data is handled
        properly for a simple case relative
        to brute force approach.'''
        actual = directed_hausdorff(self.path_2_4d, self.path_1_4d)
        # brute force over columns:
        expected = max(np.amin(scipy.spatial.distance.cdist(self.path_1_4d,
                                                            self.path_2_4d), 
                                                            axis = 0))
        assert_almost_equal(actual, expected, decimal=9)
