'''Try out simple problem that is failing.'''

import unittest

import numpy as np
from pyHiGHS import highs_wrapper

class TestHiGHS(unittest.TestCase):

    def test_simple_presolve(self):

        c = np.array([-1, 1], dtype='double')
        A = np.array([
            [1, 0],
            [0, 1],
        ], dtype='double')
        lhs = np.array([-np.inf, 2.], dtype='double')
        rhs = np.array([1., 2.], dtype='double')
        lb = np.array([1., 2.], dtype='double')
        ub = np.array([1., 2.], dtype='double')

        options = {'presolve': True, 'method': 'simplex'}
        with self.assertRaises(RuntimeError):
            res = highs_wrapper(c, A, rhs, lhs, lb, ub, options=options)

        options['presolve'] = False
        res = highs_wrapper(c, A, rhs, lhs, lb, ub, options=options)
        self.assertEqual(res['col_value'].tolist(), [1, 2])

if __name__ == '__main__':
    unittest.main()
