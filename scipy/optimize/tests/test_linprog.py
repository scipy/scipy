"""
Unit test for Linear Programming via Simplex Algorithm.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_array_almost_equal, TestCase, \
                          assert_allclose, run_module_suite
import numpy as np

from scipy.optimize import linprog, Result


class TestLinprog(TestCase):
    """
    Test SLSQP algorithm using Example 14.4 from Numerical Methods for
    Engineers by Steven Chapra and Raymond Canale.
    This example maximizes the function f(x) = 2*x*y + 2*x - x**2 - 2*y**2,
    which has a maximum at x=2, y=1.
    """
    def setUp(self):
        self.opts = {'disp': False}

    def test_linprog_upper_bound_constraints(self):
        """ Maximize a linear function subject to only linear upper bound constraints. """
        #  http://www.dam.brown.edu/people/huiwang/classes/am121/Archive/simplex_121_c.pdf
        c = [3,2]
        b_ub = [10,8,4]
        A_ub = [[2,1],
                [1,1],
                [1,0]]

        res = (linprog(c,A_ub=A_ub,b_ub=b_ub,objtype='max',disp=False))

        assert_(res.status == 0,
                "Test of linprog upper bound constraints failed.  Expected status = 0, got {:d}.".format(res.status))

        assert_array_almost_equal(res.x,np.array([2.0,6.0]),
                                  err_msg="Test of linprog upper bound constraints failed with incorrect result.")

    def test_linprog_mixed_constraints(self):
        """ Minimize linear function subject to linear, non-negative variables. """
        #  http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
        c = [6,3]
        A_lb = [[1, 1],
                [2,-1]]
        b_lb = [1,1]
        A_ub = [[0,3]]
        b_ub = [2]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,objtype='min',disp=False)

        assert_(res.status == 0,
                "Test of linprog minimization failed.  Expected status = 0, got {:d}.".format(res.status))

        assert_array_almost_equal(res.x,[2/3,1/3],
                                  err_msg="Test of linprog minimization failed with incorrect result.")

    def test_linprog_cyclic_recovery(self):
        """ Test linprogs recovery from cycling using the Klee-Minty problem """
        #  Klee-Minty  http://www.math.ubc.ca/~israel/m340/kleemin3.pdf
        c = [100,10,1]
        A_ub = [[1, 0, 0],
                [20, 1, 0],
                [200,20, 1]]

        b_ub = [1,100,10000]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,objtype='max',disp=False)

        assert_(res.status == 0,
                "Test of linprog recovery from cycling failed.  Expected status = 0, got {:d}.".format(res.status))

        assert_array_almost_equal(res.x,[0,0,10000],
                                  err_msg="Test of linprog recovery from cycling failed with incorrect result.")

    def test_linprog_unbounded(self):
        """ Test linprog response to an unbounded problem """
        c = [1,1]
        A_lb = [[1,-1],
                [1,1]]
        b_lb = [1,2]

        res = linprog(c,A_lb=A_lb,b_lb=b_lb,objtype='max',disp=False)

        assert_(res.status == 3,"Test of linprog response to an unbounded problem failed.")

    def test_linprog_infeasible(self):
        """ Test linrpog response to an infeasible problem """
        c = [1,1]

        A_lb = [[1,1]]
        b_lb = [5]

        A_ub = [[1,0],
                [0,1]]
        b_ub = [2,2]

        res = linprog(c,A_lb=A_lb,b_lb=b_lb,A_ub=A_ub,b_ub=b_ub,objtype='max',disp=False)

        assert_(not res.success,"Test of linprog with an infeasible problem errantly ended with success")

        assert_(res.status == 2,"Test of linprog with an infeasible problem did not acknowledge its infeasibility")

    def test_nontrivial_problem(self):
        """ Test linprog for a problem involving all constraint types, negative resource limits, and rounding issues. """
        c = [-1,8,4,-6]

        A_ub = [[-7,-7,6,9]]
        b_ub = [-3]

        A_eq = [[-10,1,1,-8]]
        b_eq = [-4]

        A_lb = [[-1,1,3,0],
                [-10,10,7,-7],
                [-6,1,-3,-4]]
        b_lb = [-6,6,-6]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,A_eq=A_eq,b_eq=b_eq,objtype='min',disp=False)

        assert_(res.status == 0,
                "Test of linprog with nontrivial problem failed.  Expected status = 0, got {:d}.".format(res.status))

        assert_(res.fun == 7083/1391,
                "Test of linprog with nontrivial problem converged but yielded unexpected result (fun)")

        assert_array_almost_equal(res.x,[101/1391,1462/1391,0,752/1391],
                                  err_msg="Test of linprog with nontrivial problem converged but yielded unexpected result (x)")

if __name__ == "__main__":
    run_module_suite()
