import pytest
import numpy as np
from numpy.testing import  assert_equal
from scipy.spatial.distance import continuous_frechet

class TestFrechet(object):

    def test_simple_case(self):
        # Test simple case            
        u = [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0)]   
        v = [(2.0, 2.0), (0.0, 1.0), (2.0, 4.0)]
        actual    = continuous_frechet(u,v)
        expected  = 2.0
        assert_equal(actual, expected)

    def test_invalid_dimensions(self):
        # Ensure that a ValueError is raised when the number of columns
        # is not the same.
        u = np.zeros((4, 2))
        v = np.zeros((4, 5))

        with pytest.raises(ValueError):
            continuous_frechet(u, v)
