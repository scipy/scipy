import math
from math import sin, cos, pi
import unittest

from numpy import array, arange, allclose

from mkufunc import Cfunc, genufunc, mkufunc


class Math_Tests(unittest.TestCase):
    
    def test_sin(self):
        @mkufunc([(float, float)])
        def u_sin(x):
            return sin(x)

        x = 1.23
        self.assert_(u_sin(x), sin(x))


if __name__ == '__main__':
    unittest.main()
