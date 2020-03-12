"""
Unit test for Linear Programming presolve routine.
"""

import numpy as np
from numpy.testing import assert_
from scipy.optimize._linprog_util import (
    _LPProblem,
    _presolve_infeasible_bounds,
    _presolve_infeasible_equality_constraints,
    _presolve_infeasible_inequality_constraints,
    _presolve_remove_fixed_variables,
    _presolve_remove_row_singletons)


def test_infeasible_bounds():
    """
    Test if upper bounds are strictly >= lower bounds.
    """
    lp = _LPProblem(
        c=None,
        A_ub=None,
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=None,
        x0=None
    )

    # (lp, rev, status) = _presolve_infeasible_bounds(lp)

    bounds = np.array([[0, 1], [0, 1]])
    lpp, rev, status = _presolve_infeasible_bounds(lp._replace(bounds=bounds))
    assert_(status == 6)

    bounds = np.array([[1, 1], [1, 1]])
    lpp, rev, status = _presolve_infeasible_bounds(lp._replace(bounds=bounds))
    assert_(status == 6, "Equal bounds expected to pass")

    bounds = np.array([[1, 0], [1, 1]])
    lpp, rev, status = _presolve_infeasible_bounds(lp._replace(bounds=bounds))
    assert_(status == 2)


def test_infeasible_equality_constraints():
    """
    Test if equality constraints are feasible, given upper and lower bounds.
    """
    lp = _LPProblem(
        c=None,
        A_ub=None,
        b_ub=None,
        A_eq=np.array([[1, 1], [2, -1]]),
        b_eq=None,
        bounds=np.array([[0, 1], [0, 1]]),
        x0=None
    )

    # (lp, rev, status) = _presolve_infeasible_equality_constraints(lp)

    b_eq = np.array([1, 1])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 6)

    b_eq = np.array([0, -1])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 6)

    b_eq = np.array([2, 2])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 6)

    b_eq = np.array([-1, 0])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 2)

    b_eq = np.array([3, 0])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 2)

    b_eq = np.array([1, -2])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 2)

    b_eq = np.array([1, 4])
    lpp, rev, status = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq))
    assert_(status == 2)


def test_infeasible_inequality_constraints():
    """
    Test if inequality constraints are feasible, given upper and lower bounds.
    """
    lp = _LPProblem(
        c=None,
        A_ub=np.array([[1, 1], [2, -1]]),
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=np.array([[0, 1], [0, 1]]),
        x0=None
    )

    # (lp, rev, status) = _presolve_infeasible_inequality_constraints(lp)

    b_ub = np.array([10, 10])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 6)

    b_ub = np.array([10, 10])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 6)

    b_ub = np.array([0, 10])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 6)

    b_ub = np.array([10, -1])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 6)

    b_ub = np.array([0, -1])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 6)

    b_ub = np.array([-2, 10])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 2)

    b_ub = np.array([10, -2])
    lpp, rev, status = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub))
    assert_(status == 2)


def test_remove_fixed_variables():
    """
    Test if fixed variables (lb == ub) are removed correctly.
    """
    lp = _LPProblem(
        c=np.array([1, 2, 3]),
        A_ub=np.array([[1, 3, 5], [-2, -2, 1]]),
        b_ub=np.array([1, -1]),
        A_eq=None,
        b_eq=None,
        bounds=None,
        x0=None,
    )
    tol = 1e-9

    # (lp, rev, status) = _presolve_remove_fixed_variables(lp)

    bounds = np.array([[0, 1], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), tol)
    assert_(status == 6)

    # Variable 0 fixed at 1
    bounds = np.array([[1, 1], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), tol)
    assert_(status == 5)
    assert_(np.all(lpp[0] == np.array([2, 3])))
    assert_(np.all(lpp[1] == np.array([[3, 5], [-2, 1]])))
    assert_(np.all(lpp[2] == np.array([0, 1])))

    cr, xr, x0r = rev(lpp[0], np.array([0, 0]), lpp[6])
    assert_(np.all(cr == np.array([1, 2, 3])))
    assert_(np.all(xr == np.array([1, 0, 0])))
    assert_(x0r is None)

    # Variable 1 fixed at -1
    bounds = np.array([[0, 1], [-1, -1], [0, 1]])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), tol)
    assert_(status == 5)
    assert_(np.all(lpp[0] == np.array([1, 3])))
    assert_(np.all(lpp[1] == np.array([[1, 5], [-2, 1]])))
    assert_(np.all(lpp[2] == np.array([4, -3])))

    cr, xr, x0r = rev(lpp[0], np.array([0, 0]), lpp[6])
    assert_(np.all(cr == np.array([1, 2, 3])))
    assert_(np.all(xr == np.array([0, -1, 0])))
    assert_(x0r is None)

    # Variable 2 fixed at 3
    bounds = np.array([[0, 1], [0, 1], [3, 3]])
    x0 = np.array([3, -2, 1])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds, x0=x0), tol)
    assert_(status == 5)
    assert_(np.all(lpp[0] == np.array([1, 2])))
    assert_(np.all(lpp[1] == np.array([[1, 3], [-2, -2]])))
    assert_(np.all(lpp[2] == np.array([-14, -4])))

    cr, xr, x0r = rev(lpp[0], np.array([0, 0]), lpp[6])
    assert_(np.all(cr == np.array([1, 2, 3])))
    assert_(np.all(xr == np.array([0, 0, 3])))
    assert_(np.all(x0r == x0))

    # Variables 0 and 2 fixed at 1 and 3
    bounds = np.array([[1, 1], [0, 1], [3, 3]])
    x0 = np.array([3, -2, 1])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds, x0=x0), tol)
    assert_(status == 5)
    assert_(np.all(lpp[0] == np.array([2])))
    assert_(np.all(lpp[1] == np.array([[3], [-2]])))
    assert_(np.all(lpp[2] == np.array([-15, -2])))

    cr, xr, x0r = rev(lpp[0], np.array([0]), lpp[6])
    assert_(np.all(cr == np.array([1, 2, 3])))
    assert_(np.all(xr == np.array([1, 0, 3])))
    assert_(np.all(x0r == x0))


def test_remove_fixed_variables_close_bounds():
    """
    Test tolerances with close bounds.
    """
    lp = _LPProblem(
        c=np.array([1, 2, 3]),
        A_ub=np.array([[1, 3, 5], [-2, -2, 1]]),
        b_ub=np.array([1, -1]),
        A_eq=None,
        b_eq=None,
        bounds=None,
        x0=None,
    )
    tol = 1e-2

    # (lp, rev, status) = _presolve_remove_fixed_variables(lp)
    # Variable 1 fixed at -1
    bounds = np.array([[0, 1], [-1, -1], [0, 1]])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), tol)
    assert_(status == 5)
    assert_(np.all(lpp[0] == np.array([1, 3])))
    assert_(np.all(lpp[1] == np.array([[1, 5], [-2, 1]])))
    assert_(np.all(lpp[2] == np.array([4, -3])))

    cr, xr, x0r = rev(lpp[0], np.array([0, 0]), lpp[6])
    assert_(np.all(cr == np.array([1, 2, 3])))
    assert_(np.all(xr == np.array([0, -1, 0])))
    assert_(x0r is None)

    # Variable 1 fixed at -1 within tolerance
    bounds = np.array([[0, 1], [-1 - 0.49 * tol, -1 + 0.49 * tol], [0, 1]])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), tol)
    assert_(status == 5)
    assert_(np.all(lpp[0] == np.array([1, 3])))
    assert_(np.all(lpp[1] == np.array([[1, 5], [-2, 1]])))
    assert_(np.all(lpp[2] == np.array([4, -3])))

    cr, xr, x0r = rev(lpp[0], np.array([0, 0]), lpp[6])
    assert_(np.all(cr == np.array([1, 2, 3])))
    assert_(np.all(xr == np.array([0, -1, 0])))
    assert_(x0r is None)

    # Variable 1 fixed at -1 outside tolerance
    bounds = np.array([[0, 1], [-1 - 0.51 * tol, -1 + 0.51 * tol], [0, 1]])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), tol)
    assert_(status == 6)


def test_remove_fixed_variables_solution():
    """
    Test result if all variables are fixed.
    Routine must conclude that problem is solved or that problem is infeasible.
    """
    lp = _LPProblem(
        c=np.array([1, 2, 3]),
        A_ub=np.array([[1, 3, 5], [-2, -2, 1]]),
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=np.array([[1, 1], [-2, -2], [3, 3]]),
        x0=None,
    )
    tol = 1e-2

    # A_ub * x = [10, 5]

    # b_ub exactly on border: solution
    b_ub = np.array([10, 5])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), tol)
    print("status:", status)
    assert_(status == 0)

    # b_ub too large: infeasible
    b_ub = np.array([11, 6])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), tol)
    assert_(status == 2)

    # b_ub small enough: solution
    b_ub = np.array([9, 4])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), tol)
    assert_(status == 0)

    # b_ub a little larger: infeasible since tolerance plays no role
    b_ub = np.array([10 + 0.5 * tol, 5])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), tol)
    assert_(status == 2)

    lp = _LPProblem(
        c=np.array([1, 2, 3]),
        A_ub=None,
        b_ub=None,
        A_eq=np.array([[1, 3, 5], [-2, -2, 1]]),
        b_eq=None,
        bounds=np.array([[1, 1], [-2, -2], [3, 3]]),
        x0=None,
    )
    tol = 1e-2

    # A_eq * x = [10, 5]
    # b_eq tolerance [5.91, 3]*1e-2

    # b_eq exact: solution
    b_eq = np.array([10, 5])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 0)

    # b_eq too large: infeasible
    b_eq = np.array([11, 6])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 2)

    # b_eq too small: infeasible
    b_eq = np.array([9, 4])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 2)

    # b_eq a little off, within tolerance: solution
    b_eq = np.array([10 + 1.5 * tol, 5 - 1.5 * tol])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 0)

    # b_eq a little off, within tolerance: solution
    b_eq = np.array([10 - 5.5 * tol, 5 + 2.5 * tol])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 0)

    # b_eq a little more off, outside tolerance: solution
    b_eq = np.array([10 - 5.5 * tol, 5 + 3.5 * tol])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 2)

    # b_eq a little more off, outside tolerance: solution
    b_eq = np.array([10 - 6.5 * tol, 5 + 2.5 * tol])
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), tol)
    assert_(status == 2)


def test_remove_singleton_rows():
    """
    Test singleton removal procedure.
    """
    tol = 1e-2
    lp = _LPProblem(
        c=np.array([1, 2, 3, 4, 5]),
        A_ub=None,
        b_ub=None,
        A_eq=np.array([
            [ 1,  3,  5,  7,  9],
            [ 0,  0,  0, -3,  0],
            [-2, -2,  1,  8,  1],
            [ 0,  2,  0,  2,  0],
            [ 0,  0,  0, -1,  0],
            [ 5,  0,  0,  0,  0]
            ]),
        b_eq=np.array([8, 6, -2, -3, 2 + 0.5 * tol, 4]),
        bounds=np.array([[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]]),
        x0=None,
    )

    # (lp, rev, status) = _presolve_remove_row_singletons(lp)
    lpp, rev, status = _presolve_remove_row_singletons(lp, tol)
    assert_(np.all(lpp[0] == np.array([2, 3, 5])))
    assert_(np.allclose(lpp[4], np.array([21.2, 15.6, 1.])))
    assert_(np.allclose(lpp[3], np.array([[3, 5, 9], [-2, 1, 1], [2, 0, 0]])))

    b_eq = np.array([8, 6, -2, -3, 2 + 1.5 * tol, 4])
    lpp, rev, status = _presolve_remove_row_singletons(lp._replace(b_eq=b_eq), tol)
    assert_(status == 2)
