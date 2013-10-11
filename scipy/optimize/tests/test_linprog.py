"""
Unit test for Linear Programming via Simplex Algorithm.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_array_almost_equal, TestCase, \
                          assert_allclose, run_module_suite
import numpy as np

from scipy.optimize import linprog


class TestLinprog(TestCase):
    """
    Test SLSQP algorithm using Example 14.4 from Numerical Methods for
    Engineers by Steven Chapra and Raymond Canale.
    This example maximizes the function f(x) = 2*x*y + 2*x - x**2 - 2*y**2,
    which has a maximum at x=2, y=1.
    """
    def setUp(self):
        self.opts = {'disp': False}

    #def test_linprog_standard_form(self):
    #    """ Maximize linear function subject to linear upper-bound constraints, non-negative variables. """
    #    pass
    #
    #def test_linprog_n_equals_m(self):
    #    """ Maximize linear function where number of constraints equals the number of variables. """
    #    pass

    def test_linprog_minimization(self):
        """ Minimize linear function subject to linear, non-negative variables. """
        # Two-Phase example
        c = [6,3]
        A_lb = [[1, 1],
                [2,-1]]
        b_lb = [1,1]
        A_ub = [[0,3]]
        b_ub = [2]

        #http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
        linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,objtype='min',disp=False)

    #def test_linprog_mixed_constraints(self):
    #    """ Minimize linear function subject to linear upper, lower, and equality constraints, non-negative variables. """
    #    pass
    #
    #def test_linprog_cyclic_recovery(self):
    #    """ Test linprogs recovery from cycling using the Klee-Minty problem """
    #    pass
    #
    #def test_linprog_unbounded(self):
    #    """ Test linprog response to an unbounded problem """
    #    pass
    #
    #def test_linprog_infeasible(self):
    #    """ Test linrpog response to an infeasible problem """
    #    pass

if __name__ == "__main__":
    run_module_suite()
