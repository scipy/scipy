"""
Unit test for Linear Programming via Simplex Algorithm.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import (assert_, assert_array_almost_equal, TestCase,
        assert_allclose, run_module_suite, assert_almost_equal, assert_raises,
        assert_equal)
import numpy as np

from scipy.optimize import linprog, OptimizeWarning
from scipy.lib._numpy_compat import _assert_warns


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
        # Maximize a linear function subject to only linear upper bound constraints.
        #  http://www.dam.brown.edu/people/huiwang/classes/am121/Archive/simplex_121_c.pdf
        c = np.array([3,2])*-1  # maximize
        A_ub = [[2,1],
                [1,1],
                [1,0]]
        b_ub = [10,8,4]

        res = (linprog(c,A_ub=A_ub,b_ub=b_ub))

        assert_equal(res.status, 0,
                err_msg="Test of linprog upper bound constraints failed.")

        assert_array_almost_equal(res.x,np.array([2.0,6.0]),
                                  err_msg="Test of linprog upper bound "
                                          "constraints failed with incorrect "
                                          "result.")

        assert_almost_equal(-res.fun, 18, err_msg="Test of linprog upper bound "
                "constraints converged with incorrect objective value.")

    def test_linprog_mixed_constraints(self):
        # Minimize linear function subject to non-negative variables.
        #  http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
        c = [6,3]
        A_ub = [[0, 3],
               [-1,-1],
               [-2, 1]]
        b_ub = [2,-1,-1]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_equal(res.status, 0,
                err_msg="Test of linprog with artificial variables failed.")

        assert_array_almost_equal(res.x,[2/3,1/3],
                                  err_msg="Test of linprog with artificial "
                                          "variables failed with incorrect "
                                          "result.")

        assert_almost_equal(res.fun,5,err_msg="Test of linprog with artificial "
                                              "variables failed with incorrect "
                                              "objective value.")

    def test_linprog_cyclic_recovery(self):
        # Test linprogs recovery from cycling using the Klee-Minty problem
        #  Klee-Minty  http://www.math.ubc.ca/~israel/m340/kleemin3.pdf
        c = np.array([100,10,1])*-1  # maximize
        A_ub = [[1, 0, 0],
                [20, 1, 0],
                [200,20, 1]]

        b_ub = [1,100,10000]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_equal(res.status, 0,
                err_msg="Test of linprog recovery from cycling failed.")

        assert_array_almost_equal(res.x,[0,0,10000],
                                  err_msg="Test of linprog recovery from "
                                          "cycling failed with incorrect "
                                          "result.")

    def test_linprog_cyclic_bland(self):
        # Test the effect of Bland's rule on a cycling problem
        c = np.array([-10, 57, 9, 24.])
        A_ub = np.array([[0.5, -5.5, -2.5, 9],
                         [0.5, -1.5, -0.5, 1],
                         [1, 0, 0, 0]])
        b_ub = [0, 0, 1]

        res = linprog(c, A_ub=A_ub, b_ub=b_ub,
                      options=dict(maxiter=100))
        assert_(not res.success)

        res = linprog(c, A_ub=A_ub, b_ub=b_ub,
                      options=dict(maxiter=100, bland=True,))
        assert_(res.success)
        assert_allclose(res.x, [1, 0, 1, 0])

    def test_linprog_unbounded(self):
        # Test linprog response to an unbounded problem
        c = np.array([1,1])*-1  # maximize
        A_ub = [[-1,1],
                [-1,-1]]
        b_ub = [-1,-2]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_equal(res.status, 3, err_msg="Test of linprog response to an "
                "unbounded problem failed.")

    def test_linprog_infeasible(self):
        # Test linrpog response to an infeasible problem
        c = [-1,-1]

        A_ub = [[1,0],
                [0,1],
                [-1,-1]]
        b_ub = [2,2,-5]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub)

        assert_(not res.success,"Test of linprog with an infeasible problem "
                                "errantly ended with success")

        assert_equal(res.status, 2, err_msg="Test of linprog with an "
                "infeasible problem did not acknowledge its infeasibility")

    def test_nontrivial_problem(self):
        # Test linprog for a problem involving all constraint types,
        # negative resource limits, and rounding issues.
        c = [-1,8,4,-6]

        A_ub = [[-7,-7,6,9],
                [1,-1,-3,0],
                [10,-10,-7,7],
                [6,-1,3,4]]
        b_ub = [-3,6,-6,6]

        A_eq = [[-10,1,1,-8]]
        b_eq = [-4]

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq)

        assert_equal(res.status, 0,
                err_msg="Test of linprog with nontrivial problem failed.")

        assert_almost_equal(res.fun, 7083/1391, 9,
                err_msg="Test of linprog with nontrivial problem converged "
                "but yielded unexpected result")

        assert_array_almost_equal(res.x,[101/1391,1462/1391,0,752/1391],
                                  err_msg="Test of linprog with nontrivial "
                                          "problem converged but yielded "
                                          "unexpected result.")

    def test_negative_variable(self):
        # Test linprog with a problem with one unbounded variable and
        # another with a negative lower bound.
        c = np.array([-1,4])*-1  # maximize

        A_ub = [[-3,1],
                [1,2]]

        b_ub = [6,4]

        x0_bounds = (-np.inf,np.inf)
        x1_bounds = (-3,np.inf)

        res = linprog(c,A_ub=A_ub,b_ub=b_ub,bounds=(x0_bounds,x1_bounds))

        assert_equal(res.status, 0,
                err_msg="Test of linprog with negative variable failed.")

        assert_allclose(-res.fun,80/7,err_msg="Test of linprog with negative "
                                              "variable converged but yielded "
                                              "unexpected result.")

        assert_array_almost_equal(res.x,[-8/7,18/7],
                                  err_msg="Test of linprog with negative "
                                          "variable converged but yielded "
                                          "unexpected result")

    def test_large_problem(self):
        # Test linprog simplex with a rather large problem (400 variables,
        # 40 constraints) generated by https://gist.github.com/denis-bz/8647461
        A,b,c = lpgen_2d(20,20)
        res = linprog(c,A_ub=A,b_ub=b)

        assert_equal(res.status, 0,
                err_msg="Test of linprog with large problem failed.")

        assert_almost_equal(res.fun,-64.049494229,
                            err_msg="Test of linprog with 400 x 40 problem"
                                    "gave incorrect solution")

    def test_network_flow(self):
        # A network flow problem with supply and demand at nodes
        # and with costs along directed edges.
        # https://www.princeton.edu/~rvdb/542/lectures/lec10.pdf
        c = [2, 4, 9, 11, 4, 3, 8, 7, 0, 15, 16, 18]
        n, p = -1, 1
        A_eq = [
                [n, n, p, 0, p, 0, 0, 0, 0, p, 0, 0],
                [p, 0, 0, p, 0, p, 0, 0, 0, 0, 0, 0],
                [0, 0, n, n, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, p, p, 0, 0, p, 0],
                [0, 0, 0, 0, n, n, n, 0, p, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, n, n, 0, 0, p],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, n, n, n]]
        b_eq = [0, 19, -16, 33, 0, 0, -36]
        res = linprog(c=c, A_eq=A_eq, b_eq=b_eq)

        assert_equal(res.status, 0,
                err_msg="Test of linprog solution of network flow failed.")

        assert_allclose(res.fun, 755,
                err_msg="Test of linprog solution of network flow "
                "converged but yielded unexpected total cost.")

    def test_network_flow_limited_capacity(self):
        # A network flow problem with supply and demand at nodes
        # and with costs and capacities along directed edges.
        # http://blog.sommer-forst.de/2013/04/10/
        cost = [2, 2, 1, 3, 1]
        bounds = [
                [0, 4],
                [0, 2],
                [0, 2],
                [0, 3],
                [0, 5]]
        n, p = -1, 1
        A_eq = [
                [n, n, 0, 0, 0],
                [p, 0, n, n, 0],
                [0, p, p, 0, n],
                [0, 0, 0, p, p]]
        b_eq = [-4, 0, 0, 4]
        res = linprog(c=cost, A_eq=A_eq, b_eq=b_eq, bounds=bounds)

        assert_equal(res.status, 0,
                err_msg="Test of linprog solution of network flow "
                "with limited capacity failed.")

        assert_allclose(res.fun, 14,
                err_msg="Test of linprog solution of network flow "
                "with limited capacity converged but yielded unexpected "
                "total cost.")

    def test_simplex_algorithm_wikipedia_example(self):
        # http://en.wikipedia.org/wiki/Simplex_algorithm#Example
        Z = [-2, -3, -4]
        A_ub = [
                [3, 2, 1],
                [2, 5, 3]]
        b_ub = [10, 15]
        res = linprog(c=Z, A_ub=A_ub, b_ub=b_ub)

        assert_equal(res.status, 0,
                err_msg="Test of linprog solution of Wikipedia example failed.")

        assert_allclose(res.fun, -20,
                err_msg="Test of linprog solution of Wikipedia example "
                "converged but yielded unexpected objective value.")

    def test_callback(self):
        # Check that callback is as advertised
        callback_complete = [False]
        last_xk = []

        def cb(xk, **kwargs):
            kwargs.pop('tableau')
            assert_(isinstance(kwargs.pop('phase'), int))
            assert_(isinstance(kwargs.pop('nit'), int))

            i, j = kwargs.pop('pivot')
            assert_(np.isscalar(i))
            assert_(np.isscalar(j))

            basis = kwargs.pop('basis')
            assert_(isinstance(basis, np.ndarray))
            assert_(basis.dtype == np.int_)

            complete = kwargs.pop('complete')
            assert_(isinstance(complete, bool))
            if complete:
                last_xk.append(xk)
                callback_complete[0] = True
            else:
                assert_(not callback_complete[0])

            # no more kwargs
            assert_(not kwargs)
        
        c = np.array([-3,-2])
        A_ub = [[2,1], [1,1], [1,0]]
        b_ub = [10,8,4]
        res = linprog(c,A_ub=A_ub,b_ub=b_ub, callback=cb)

        assert_(callback_complete[0])
        assert_allclose(last_xk[0], res.x)

    def test_unknown_options_or_solver(self):
        c = np.array([-3,-2])
        A_ub = [[2,1], [1,1], [1,0]]
        b_ub = [10,8,4]

        _assert_warns(OptimizeWarning, linprog,
                c, A_ub=A_ub, b_ub=b_ub, options=dict(spam='42'))

        assert_raises(ValueError, linprog,
                c, A_ub=A_ub, b_ub=b_ub, method='ekki-ekki-ekki')

    def test_no_constraints(self):
        res = linprog([-1, -2])
        assert_allclose(res.x, [0, 0])

    def test_simple_bounds(self):
        res = linprog([1, 2], bounds=(1, 2))
        assert_allclose(res.x, [1, 1])
        res = linprog([1, 2], bounds=[(1, 2), (1, 2)])
        assert_allclose(res.x, [1, 1])

    def test_invalid_inputs(self):
        for bad_bound in [[(5, 0), (1, 2), (3, 4)],
                          [(1, 2), (3, 4)],
                          [(1, 2), (3, 4), (3, 4, 5)],
                          [(1, 2), (np.inf, np.inf), (3, 4)],
                          [(1, 2), (-np.inf, -np.inf), (3, 4)],
                          ]:
            assert_raises(ValueError, linprog,
                          [1, 2, 3], bounds=bad_bound)

        assert_raises(ValueError, linprog, [1,2], A_ub=[[1,2]], b_ub=[1,2])
        assert_raises(ValueError, linprog, [1,2], A_ub=[[1]], b_ub=[1])
        assert_raises(ValueError, linprog, [1,2], A_eq=[[1,2]], b_eq=[1,2])
        assert_raises(ValueError, linprog, [1,2], A_eq=[[1]], b_eq=[1])
        assert_raises(ValueError, linprog, [1,2], A_eq=[1], b_eq=1)
        assert_raises(ValueError, linprog, [1,2], A_ub=np.zeros((1,1,3)), b_eq=1)

if __name__ == "__main__":
    run_module_suite()
