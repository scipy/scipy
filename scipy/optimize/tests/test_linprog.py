"""
Unit test for Linear Programming via Simplex Algorithm.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_array_almost_equal, TestCase, \
                          assert_allclose, run_module_suite, assert_almost_equal
import numpy as np

from scipy.optimize import linprog, linprog_verbose_callback


def lpgen_2d(m,n):
    """ -> A b c LP test: m*n vars, m+n constraints
        row sums == n/m, col sums == 1
        https://gist.github.com/denis-bz/8647461
    """
    np.random.seed(0)
    c = - np.random.exponential(size=(m,n))
    Arow = np.zeros((m,m*n))
    brow = np.zeros(m)
    for j in range(m):
        j1 = j + 1
        Arow[j,j*n:j1*n] = 1
        brow[j] = n/m

    Acol = np.zeros((n,m*n))
    bcol = np.zeros(n)
    for j in range(n):
        j1 = j + 1
        Acol[j,j::n] = 1
        bcol[j] = 1

    A = np.vstack((Arow,Acol))
    b = np.hstack((brow,bcol))

    return A, b, c.ravel()


class TestLinprog(TestCase):
    """  Test the linprog routines using a variety of example problems.
    """
    def setUp(self):
        self.opts = {'disp': False}

    def test_linprog_upper_bound_constraints(self):
        """  Maximize a linear function subject to only linear upper bound
        constraints. """
        #  http://www.dam.brown.edu/people/huiwang/classes/am121/Archive/simplex_121_c.pdf
        c = np.array([3,2])*-1  # maximize
        A_ub = [[2,1],
                [1,1],
                [1,0]]
        b_ub = [10,8,4]

        res = (linprog(c,A_ub=A_ub,b_ub=b_ub))

        assert_(res.status == 0,
                "Test of linprog upper bound constraints failed.  "
                "Expected status = 0, got %d." % res.status)

        assert_array_almost_equal(res.x,np.array([2.0,6.0]),
                                  err_msg="Test of linprog upper bound "
                                          "constraints failed with incorrect "
                                          "result.")

        assert_almost_equal(-res.fun,18,err_msg="Test of linprog upper bound"
                            " constraints failed. "
                            "Expected f=18, got %f." % res.fun)

    def test_linprog_mixed_constraints(self):
        """ Minimize linear function subject to non-negative variables. """
        #  http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
        c = [6,3]
        A_ub = [[0, 3],
               [-1,-1],
               [-2, 1]]
        b_ub = [2,-1,-1]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_(res.status == 0,
                "Test of linprog with artificial variables failed.  "
                "Expected status = 0, got %d." % res.status)

        assert_array_almost_equal(res.x,[2/3,1/3],
                                  err_msg="Test of linprog with artificial "
                                          "variables failed with incorrect "
                                          "result.")

        assert_almost_equal(res.fun,5,err_msg="Test of linprog with artificial "
                                              "variables failed with incorrect "
                                              "objective value.")

    def test_linprog_cyclic_recovery(self):
        """  Test linprogs recovery from cycling using the Klee-Minty problem """
        #  Klee-Minty  http://www.math.ubc.ca/~israel/m340/kleemin3.pdf
        c = np.array([100,10,1])*-1  # maximize
        A_ub = [[1, 0, 0],
                [20, 1, 0],
                [200,20, 1]]

        b_ub = [1,100,10000]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_(res.status == 0,
                "Test of linprog recovery from cycling failed.  Expected status"
                " = 0, got %d." % res.status)

        assert_array_almost_equal(res.x,[0,0,10000],
                                  err_msg="Test of linprog recovery from "
                                          "cycling failed with incorrect "
                                          "result.")

    def test_linprog_unbounded(self):
        """  Test linprog response to an unbounded problem """
        c = np.array([1,1])*-1  # maximize
        A_ub = [[-1,1],
                [-1,-1]]
        b_ub = [-1,-2]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_(res.status == 3,"Test of linprog response to an "
                                "unbounded problem failed.")

    def test_linprog_infeasible(self):
        """  Test linrpog response to an infeasible problem """
        c = [-1,-1]

        A_ub = [[1,0],
                [0,1],
                [-1,-1]]
        b_ub = [2,2,-5]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_(not res.success,"Test of linprog with an infeasible problem "
                                "errantly ended with success")

        assert_(res.status == 2,"Test of linprog with an infeasible "
                                "problem did not acknowledge its infeasibility")

    def test_nontrivial_problem(self):
        """  Test linprog for a problem involving all constraint types,
        negative resource limits, and rounding issues. """
        c = [-1,8,4,-6]

        A_ub = [[-7,-7,6,9],
                [1,-1,-3,0],
                [10,-10,-7,7],
                [6,-1,3,4]]
        b_ub = [-3,6,-6,6]

        A_eq = [[-10,1,1,-8]]
        b_eq = [-4]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq)

        assert_(res.status == 0,
                "Test of linprog with nontrivial problem failed.  "
                "Expected status = 0, got %d." % res.status)

        assert_almost_equal(res.fun,7083/1391,9,
                "Test of linprog with nontrivial problem converged but yielded "
                "unexpected result (%f)" % res.fun)

        assert_array_almost_equal(res.x,[101/1391,1462/1391,0,752/1391],
                                  err_msg="Test of linprog with nontrivial "
                                          "problem converged but yielded "
                                          "unexpected result.")

    def test_negative_variable(self):
        """  Test linprog with a problem with one unbounded variable and
        another with a negative lower bound. """
        c = np.array([-1,4])*-1  # maximize

        A_ub = [[-3,1],
                [1,2]]

        b_ub = [6,4]

        x0_bounds = (-np.inf,np.inf)
        x1_bounds = (-3,np.inf)

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,bounds=(x0_bounds,x1_bounds))

        assert_(res.status == 0,
                "Test of linprog with negative variable failed.  "
                "Expected status = 0, got %d." % res.status)

        assert_allclose(-res.fun,80/7,err_msg="Test of linprog with negative "
                                              "variable converged but yielded "
                                              "unexpected result.")

        assert_array_almost_equal(res.x,[-8/7,18/7],
                                  err_msg="Test of linprog with negative "
                                          "variable converged but yielded "
                                          "unexpected result")

    def test_large_problem(self):
        """ Test linprog simplex with a rather large problem (400 variables,
        40 constraints) generated by https://gist.github.com/denis-bz/8647461
        """
        A,b,c = lpgen_2d(20,20)
        res = linprog(c,A_ub=A,b_ub=b)

        assert_(res.status == 0,
                "Test of linprog with large problem failed.  "
                "Expected status = 0, got %d." % res.status)

        assert_almost_equal(res.fun,-64.049494229,
                            err_msg="Test of linprog with 400 x 40 problem"
                                    "gave incorrect solution")


if __name__ == "__main__":
    run_module_suite()
