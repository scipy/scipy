"""
Unit test for Linear Programming via Simplex Algorithm.
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_, assert_array_almost_equal, assert_allclose,
        assert_almost_equal, assert_raises, assert_equal, run_module_suite)

from scipy.optimize import linprog, OptimizeWarning
from scipy._lib._numpy_compat import _assert_warns


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


def _assert_infeasible(res):
    # res: linprog result object
    assert_(not res.success, "incorrectly reported success")
    assert_equal(res.status, 2, "failed to report infeasible status")


def _assert_unbounded(res):
    # res: linprog result object
    assert_(not res.success, "incorrectly reported success")
    assert_equal(res.status, 3, "failed to report unbounded status")


def _assert_success(res, desired_fun=None, desired_x=None):
    # res: linprog result object
    # desired_fun: desired objective function value or None
    # desired_x: desired solution or None
    assert_(res.success)
    assert_equal(res.status, 0)
    if desired_fun is not None:
        assert_allclose(res.fun, desired_fun,
                        err_msg="converged to an unexpected objective value")
    if desired_x is not None:
        assert_allclose(res.x, desired_x,
                        err_msg="converged to an unexpected solution")


def test_linprog_upper_bound_constraints():
    # Maximize a linear function subject to only linear upper bound constraints.
    #  http://www.dam.brown.edu/people/huiwang/classes/am121/Archive/simplex_121_c.pdf
    c = np.array([3,2])*-1  # maximize
    A_ub = [[2,1],
            [1,1],
            [1,0]]
    b_ub = [10,8,4]
    res = (linprog(c,A_ub=A_ub,b_ub=b_ub))
    _assert_success(res, desired_fun=-18, desired_x=[2, 6])


def test_linprog_mixed_constraints():
    # Minimize linear function subject to non-negative variables.
    #  http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
    c = [6,3]
    A_ub = [[0, 3],
           [-1,-1],
           [-2, 1]]
    b_ub = [2,-1,-1]
    res = linprog(c,A_ub=A_ub,b_ub=b_ub)
    _assert_success(res, desired_fun=5, desired_x=[2/3, 1/3])


def test_linprog_cyclic_recovery():
    # Test linprogs recovery from cycling using the Klee-Minty problem
    #  Klee-Minty  http://www.math.ubc.ca/~israel/m340/kleemin3.pdf
    c = np.array([100,10,1])*-1  # maximize
    A_ub = [[1, 0, 0],
            [20, 1, 0],
            [200,20, 1]]
    b_ub = [1,100,10000]
    res = linprog(c,A_ub=A_ub,b_ub=b_ub)
    _assert_success(res, desired_x=[0, 0, 10000])


def test_linprog_cyclic_bland():
    # Test the effect of Bland's rule on a cycling problem
    c = np.array([-10, 57, 9, 24.])
    A_ub = np.array([[0.5, -5.5, -2.5, 9],
                     [0.5, -1.5, -0.5, 1],
                     [1, 0, 0, 0]])
    b_ub = [0, 0, 1]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, options=dict(maxiter=100))
    assert_(not res.success)
    res = linprog(c, A_ub=A_ub, b_ub=b_ub,
                  options=dict(maxiter=100, bland=True,))
    _assert_success(res, desired_x=[1, 0, 1, 0])


def test_linprog_unbounded():
    # Test linprog response to an unbounded problem
    c = np.array([1,1])*-1  # maximize
    A_ub = [[-1,1],
            [-1,-1]]
    b_ub = [-1,-2]
    res = linprog(c,A_ub=A_ub,b_ub=b_ub)
    _assert_unbounded(res)


def test_linprog_infeasible():
    # Test linrpog response to an infeasible problem
    c = [-1,-1]
    A_ub = [[1,0],
            [0,1],
            [-1,-1]]
    b_ub = [2,2,-5]
    res = linprog(c,A_ub=A_ub,b_ub=b_ub)
    _assert_infeasible(res)


def test_nontrivial_problem():
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
    _assert_success(res, desired_fun=7083/1391,
                    desired_x=[101/1391,1462/1391,0,752/1391])


def test_negative_variable():
    # Test linprog with a problem with one unbounded variable and
    # another with a negative lower bound.
    c = np.array([-1,4])*-1  # maximize
    A_ub = np.array([[-3,1],
                     [1, 2]], dtype=np.float64)
    A_ub_orig = A_ub.copy()
    b_ub = [6,4]
    x0_bounds = (-np.inf,np.inf)
    x1_bounds = (-3,np.inf)

    res = linprog(c,A_ub=A_ub,b_ub=b_ub,bounds=(x0_bounds,x1_bounds))

    assert_equal(A_ub, A_ub_orig)   # user input not overwritten
    _assert_success(res, desired_fun=-80/7, desired_x=[-8/7, 18/7])


def test_large_problem():
    # Test linprog simplex with a rather large problem (400 variables,
    # 40 constraints) generated by https://gist.github.com/denis-bz/8647461
    A,b,c = lpgen_2d(20,20)
    res = linprog(c,A_ub=A,b_ub=b)
    _assert_success(res, desired_fun=-64.049494229)


def test_network_flow():
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
    _assert_success(res, desired_fun=755)


def test_network_flow_limited_capacity():
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
    # Including the callback here ensures the solution can be
    # calculated correctly, even when phase 1 terminated
    # with some of the artificial variables as pivots
    # (i.e. basis[:m] contains elements corresponding to
    # the artificial variables)
    res = linprog(c=cost, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                  callback=lambda x, **kwargs: None)
    _assert_success(res, desired_fun=14)


def test_simplex_algorithm_wikipedia_example():
    # http://en.wikipedia.org/wiki/Simplex_algorithm#Example
    Z = [-2, -3, -4]
    A_ub = [
            [3, 2, 1],
            [2, 5, 3]]
    b_ub = [10, 15]
    res = linprog(c=Z, A_ub=A_ub, b_ub=b_ub)
    _assert_success(res, desired_fun=-20)


def test_enzo_example():
    # http://projects.scipy.org/scipy/attachment/ticket/1252/lp2.py
    #
    # Translated from Octave code at:
    # http://www.ecs.shimane-u.ac.jp/~kyoshida/lpeng.htm
    # and placed under MIT licence by Enzo Michelangeli
    # with permission explicitly granted by the original author,
    # Prof. Kazunobu Yoshida  
    c = [4, 8, 3, 0, 0, 0]
    A_eq = [
            [2, 5, 3, -1, 0, 0],
            [3, 2.5, 8, 0, -1, 0],
            [8, 10, 4, 0, 0, -1]]
    b_eq = [185, 155, 600]
    res = linprog(c=c, A_eq=A_eq, b_eq=b_eq)
    _assert_success(res, desired_fun=317.5,
                    desired_x=[66.25, 0, 17.5, 0, 183.75, 0])


def test_enzo_example_b():
    # rescued from https://github.com/scipy/scipy/pull/218
    c = [2.8, 6.3, 10.8, -2.8, -6.3, -10.8]
    A_eq = [[-1, -1, -1, 0, 0, 0],
            [0, 0, 0, 1, 1, 1],
            [1, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 1, 0],
            [0, 0, 1, 0, 0, 1]]
    b_eq = [-0.5, 0.4, 0.3, 0.3, 0.3]
    # Including the callback here ensures the solution can be
    # calculated correctly.
    res = linprog(c=c, A_eq=A_eq, b_eq=b_eq,
                  callback=lambda x, **kwargs: None)
    _assert_success(res, desired_fun=-1.77,
                    desired_x=[0.3, 0.2, 0.0, 0.0, 0.1, 0.3])


def test_enzo_example_c_with_degeneracy():
    # rescued from https://github.com/scipy/scipy/pull/218
    m = 20
    c = -np.ones(m)
    tmp = 2*np.pi*np.arange(1, m+1)/(m+1)
    A_eq = np.vstack((np.cos(tmp)-1, np.sin(tmp)))
    b_eq = [0, 0]
    res = linprog(c=c, A_eq=A_eq, b_eq=b_eq)
    _assert_success(res, desired_fun=0, desired_x=np.zeros(m))


def test_enzo_example_c_with_unboundedness():
    # rescued from https://github.com/scipy/scipy/pull/218
    m = 50
    c = -np.ones(m)
    tmp = 2*np.pi*np.arange(m)/(m+1)
    A_eq = np.vstack((np.cos(tmp)-1, np.sin(tmp)))
    b_eq = [0, 0]
    res = linprog(c=c, A_eq=A_eq, b_eq=b_eq)
    _assert_unbounded(res)


def test_enzo_example_c_with_infeasibility():
    # rescued from https://github.com/scipy/scipy/pull/218
    m = 50
    c = -np.ones(m)
    tmp = 2*np.pi*np.arange(m)/(m+1)
    A_eq = np.vstack((np.cos(tmp)-1, np.sin(tmp)))
    b_eq = [1, 1]
    res = linprog(c=c, A_eq=A_eq, b_eq=b_eq)
    _assert_infeasible(res)


def test_callback():
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


def test_unknown_options_or_solver():
    c = np.array([-3,-2])
    A_ub = [[2,1], [1,1], [1,0]]
    b_ub = [10,8,4]

    _assert_warns(OptimizeWarning, linprog,
                  c, A_ub=A_ub, b_ub=b_ub, options=dict(spam='42'))

    assert_raises(ValueError, linprog,
                  c, A_ub=A_ub, b_ub=b_ub, method='ekki-ekki-ekki')


def test_no_constraints():
    res = linprog([-1, -2])
    assert_equal(res.x, [0, 0])
    _assert_unbounded(res)


def test_simple_bounds():
    res = linprog([1, 2], bounds=(1, 2))
    _assert_success(res, desired_x=[1, 1])
    res = linprog([1, 2], bounds=[(1, 2), (1, 2)])
    _assert_success(res, desired_x=[1, 1])


def test_invalid_inputs():
    for bad_bound in [[(5, 0), (1, 2), (3, 4)],
                      [(1, 2), (3, 4)],
                      [(1, 2), (3, 4), (3, 4, 5)],
                      [(1, 2), (np.inf, np.inf), (3, 4)],
                      [(1, 2), (-np.inf, -np.inf), (3, 4)],
                      ]:
        assert_raises(ValueError, linprog, [1, 2, 3], bounds=bad_bound)

    assert_raises(ValueError, linprog, [1,2], A_ub=[[1,2]], b_ub=[1,2])
    assert_raises(ValueError, linprog, [1,2], A_ub=[[1]], b_ub=[1])
    assert_raises(ValueError, linprog, [1,2], A_eq=[[1,2]], b_eq=[1,2])
    assert_raises(ValueError, linprog, [1,2], A_eq=[[1]], b_eq=[1])
    assert_raises(ValueError, linprog, [1,2], A_eq=[1], b_eq=1)
    assert_raises(ValueError, linprog, [1,2], A_ub=np.zeros((1,1,3)), b_eq=1)


def test_basic_artificial_vars():
    # Test if linprog succeeds when at the end of Phase 1 some artificial
    # variables remain basic, and the row in T corresponding to the
    # artificial variables is not all zero.
    c = np.array([-0.1, -0.07, 0.004, 0.004, 0.004, 0.004])
    A_ub = np.array([[1.0, 0, 0, 0, 0, 0], [-1.0, 0, 0, 0, 0, 0],
                     [0, -1.0, 0, 0, 0, 0], [0, 1.0, 0, 0, 0, 0],
                     [1.0, 1.0, 0, 0, 0, 0]])
    b_ub = np.array([3.0, 3.0, 3.0, 3.0, 20.0])
    A_eq = np.array([[1.0, 0, -1, 1, -1, 1], [0, -1.0, -1, 1, -1, 1]])
    b_eq = np.array([0, 0])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  callback=lambda x, **kwargs: None)
    _assert_success(res, desired_fun=0, desired_x=np.zeros_like(c))


if __name__ == '__main__':
    run_module_suite()
