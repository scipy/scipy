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
    _presolve_remove_equation_row_singletons,
    _presolve_remove_inequality_row_singletons,
    _presolve_remove_empty_rows,
    _presolve_remove_empty_columns)


def test_infeasible_bounds():
    """
    Test if upper bounds are strictly >= lower bounds.
    """

    # Problem parameters in preferred format
    lp = _LPProblem(
        c=np.zeros((2), dtype=float),  # dummy, irrelavant
        A_ub=np.empty((0, 2), dtype=float),
        b_ub=np.empty((0), dtype=float),
        A_eq=np.empty((0, 2), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.empty((0, 2), dtype=float),  # supplied below
        x0=None
    )
    tol = 1e-9

    # Valid finite bounds
    bounds = np.array([[0, 1], [0, 1]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(problem_feasible)

    # Valid finite bounds, first variable fixed
    bounds = np.array([[1, 1], [1, 1]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(problem_feasible, "Equal bounds expected to pass")

    # Invalid finite bounds, first lower bound lower than upper bound
    bounds = np.array([[1, 0], [1, 1]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(not problem_feasible)

    # Valid infinite bounds
    bounds = np.array([[1, np.inf], [1, 1]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(problem_feasible)

    # Valid infinite bounds
    bounds = np.array([[1, 5], [-np.inf, 1]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(problem_feasible)

    # Invalid infinite bounds
    bounds = np.array([[np.inf, np.inf], [1, 1]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(not problem_feasible)

    # Invalid infinite bounds
    bounds = np.array([[0, 2], [-np.inf, -np.inf]])
    problem_feasible = _presolve_infeasible_bounds(lp._replace(bounds=bounds), tol)
    assert_(not problem_feasible)


def test_infeasible_equality_constraints():
    """
    Test if equality constraints are feasible, given upper and lower bounds.
    """

    # Problem parameters in preferred format
    lp = _LPProblem(
        c=np.zeros((2), dtype=float),  # values irrelevant
        A_ub=np.empty((0, 2), dtype=float),
        b_ub=np.empty((0), dtype=float),
        A_eq=np.array([[1, 1], [2, -1]], dtype=float),
        b_eq=np.empty((0), dtype=float),  # dummy, values supplied below
        bounds=np.array([[0, 1], [0, 1]], dtype=float),
        x0=None
    )
    tol = 1e-9

    # A_eq * x = b_eq
    # Using bounds, A_eq * x is between [0, -1] and [2, 2]
    b_eq = np.array([1, 1])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(problem_feasible)

    b_eq = np.array([0, -1])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(problem_feasible)

    b_eq = np.array([2, 2])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(problem_feasible)

    b_eq = np.array([-1, 0])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(not problem_feasible)

    b_eq = np.array([3, 0])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(not problem_feasible)

    b_eq = np.array([1, -2])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(not problem_feasible)

    b_eq = np.array([1, 4])
    problem_feasible = _presolve_infeasible_equality_constraints(lp._replace(b_eq=b_eq), tol)
    assert_(not problem_feasible)


def test_infeasible_inequality_constraints():
    """
    Test if inequality constraints are feasible, given upper and lower bounds.
    """

    # Problem parameters in preferred format
    lp = _LPProblem(
        c=np.zeros((2), dtype=float),  # values irrelevant
        A_ub=np.array([[1, 1], [2, -1]], dtype=float),
        b_ub=np.empty((0), dtype=float),  # dummy, values supplied below,
        A_eq=np.empty((0, 2), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.array([[0, 1], [0, 1]], dtype=float),
        x0=None
    )
    tol = 1e-9

    # A_ub * x <= b_ub
    # Using bounds, A_ub * x has lower bound [0, -1]
    b_ub = np.array([10, 10])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(problem_feasible)

    b_ub = np.array([0, 10])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(problem_feasible)

    b_ub = np.array([10, -1])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(problem_feasible)

    b_ub = np.array([0, -1])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(problem_feasible)

    b_ub = np.array([-0.5*tol, -1])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(problem_feasible)

    b_ub = np.array([0, -1 - 1.5 * tol])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(not problem_feasible)

    b_ub = np.array([-2, 10])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(not problem_feasible)

    b_ub = np.array([10, -2])
    problem_feasible = _presolve_infeasible_inequality_constraints(lp._replace(b_ub=b_ub), tol)
    assert_(not problem_feasible)


def test_remove_fixed_variables():
    """
    Test if fixed variables (lb == ub) are removed correctly.
    """
    lp = _LPProblem(
        c=np.array([1, 2, 3], dtype=float),
        A_ub=np.array([[1, 3, 5], [-2, -2, 1]], dtype=float),
        b_ub=np.array([1, -1], dtype=float),
        A_eq=np.empty((0, 3), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.empty((0, 2), dtype=float),  # supplied below
        x0=None,
    )
    tol = 1e-9

    # No variables fixed
    bounds = np.array([[0, 1], [0, 1], [0, 1]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    assert_(rev is None)

    # Variable 0 fixed at 1
    bounds = np.array([[1, 1], [0, 1], [0, 1]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([2, 3], dtype=float)))
    assert_(np.all(lpp.A_ub == np.array([[3, 5], [-2, 1]], dtype=float)))
    assert_(np.all(lpp.b_ub == np.array([0, 1], dtype=float)))

    xr = rev(np.array([0., 0.], dtype=float))
    assert_(np.all(xr == np.array([1, 0, 0], dtype=float)))

    # Variable 1 fixed at -1
    bounds = np.array([[0, 1], [-1, -1], [0, 1]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([1, 3], dtype=float)))
    assert_(np.all(lpp.A_ub == np.array([[1, 5], [-2, 1]], dtype=float)))
    assert_(np.all(lpp.b_ub == np.array([4, -3], dtype=float)))

    xr = rev(np.array([0., 0.], dtype=float))
    assert_(np.all(xr == np.array([0, -1, 0], dtype=float)))

    # Variable 2 fixed at 3
    bounds = np.array([[0, 1], [0, 1], [3, 3]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([1, 2], dtype=float)))
    assert_(np.all(lpp.A_ub == np.array([[1, 3], [-2, -2]], dtype=float)))
    assert_(np.all(lpp.b_ub == np.array([-14, -4], dtype=float)))

    xr = rev(np.array([0., 0.], dtype=float))
    assert_(np.all(xr == np.array([0, 0, 3], dtype=float)))

    # Variables 0 and 2 fixed at 1 and 3
    bounds = np.array([[1, 1], [0, 1], [3, 3]], dtype=float)
    x0 = np.array([3, -2, 1])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([2], dtype=float)))
    assert_(np.all(lpp.A_ub == np.array([[3], [-2]], dtype=float)))
    assert_(np.all(lpp.b_ub == np.array([-15, -2], dtype=float)))

    xr = rev(np.array([0.]))
    assert_(np.all(xr == np.array([1, 0, 3])))

    # All variables fixed (at 1, -1 and 3)
    bounds = np.array([[1, 1], [-1, -1], [3, 3]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])  # Infeasibility not detected here
    assert_(status['bounded'])
    assert_(not status['loop'])
    assert_(lpp.c.size == 0)
    assert_(lpp.A_ub.size == 0)
    assert_(np.all(lpp.b_ub == np.array([-12, -4], dtype=float)))

    xr = rev(np.array([], dtype=float))
    assert_(np.all(xr == np.array([1, -1, 3], dtype=float)))


def test_remove_fixed_variables_close_bounds():
    """
    Test tolerances with close bounds.
    """
    lp = _LPProblem(
        c=np.array([1, 2, 3], dtype=float),
        A_ub=np.array([[1, 3, 5], [-2, -2, 1]], dtype=float),
        b_ub=np.array([1, -1], dtype=float),
        A_eq=np.empty((0, 3), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.empty((0, 2), dtype=float),  # supplied below
        x0=None,
    )
    tol = 1e-2

    # Variable 1 fixed at -1
    bounds = np.array([[0, 1], [-1, -1], [0, 1]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([1, 3], dtype=float)))
    assert_(np.all(lpp.A_ub == np.array([[1, 5], [-2, 1]], dtype=float)))
    assert_(np.all(lpp.b_ub == np.array([4, -3], dtype=float)))

    xr = rev(np.array([0., 0.], dtype=float))
    assert_(np.all(xr == np.array([0, -1, 0], dtype=float)))

    # Variable 1 fixed at -1 within tolerance
    bounds = np.array([[0, 1], [-1 - 0.49 * tol, -1 + 0.49 * tol], [0, 1]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([1, 3], dtype=float)))
    assert_(np.all(lpp.A_ub == np.array([[1, 5], [-2, 1]], dtype=float)))
    assert_(np.all(lpp.b_ub == np.array([4, -3], dtype=float)))

    xr = rev(np.array([0., 0.], dtype=float))
    assert_(np.all(xr == np.array([0, -1, 0], dtype=float)))

    # Variable 1 'fixed' at -1 outside tolerance
    bounds = np.array([[0, 1], [-1 - 0.51 * tol, -1 + 0.51 * tol], [0, 1]], dtype=float)
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(bounds=bounds), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])


def test_remove_fixed_variables_solution():
    """
    Test result if all variables are fixed.
    Routine must conclude that problem is solved or that problem is infeasible.
    """
    lp = _LPProblem(
        c=np.array([1, 2, 3], dtype=float),
        A_ub=np.array([[1, 3, 5], [-2, -2, 1]], dtype=float),
        b_ub=np.empty((0), dtype=float),  # supplied below
        A_eq=np.empty((0, 3), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.array([[1, 1], [-2, -2], [3, 3]], dtype=float),
        x0=None,
    )
    tol = 1e-2

    # A_ub * x = [10, 5] <= b_ub

    # b_ub exactly on border: solution
    b_ub = np.array([10, 5])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_inequality_constraints(lpp, tol)
    assert_(problem_feasible)

    # b_ub too small: infeasible
    b_ub = np.array([9, 4])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_inequality_constraints(lpp, tol)
    assert_(not problem_feasible)

    # b_ub sufficiently large: solution
    b_ub = np.array([11, 6])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_inequality_constraints(lpp, tol)
    assert_(problem_feasible)

    # b_ub a little smaller: infeasible since tolerance plays no role
    b_ub = np.array([10 - 0.5 * tol, 5])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_ub=b_ub), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_inequality_constraints(lpp, tol)
    assert_(not problem_feasible)

    lp = _LPProblem(
        c=np.array([1, 2, 3], dtype=float),
        A_ub=np.empty((0, 3), dtype=float),
        b_ub=np.empty((0), dtype=float),
        A_eq=np.array([[1, 3, 5], [-2, -2, 1]], dtype=float),
        b_eq=np.empty((0), dtype=float),  # supplied below
        bounds=np.array([[1, 1], [-2, -2], [3, 3]], dtype=float),
        x0=None,
    )
    tol = 1e-2

    # A_eq * x = [10, 5]

    # b_eq exact: solution
    b_eq = np.array([10, 5])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_equality_constraints(lpp, tol)
    assert_(problem_feasible)

    # b_eq too large: infeasible
    b_eq = np.array([11, 6])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_equality_constraints(lpp, tol)
    assert_(not problem_feasible)

    # b_eq too small: infeasible
    b_eq = np.array([9, 4])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_equality_constraints(lpp, tol)
    assert_(not problem_feasible)

    # b_eq a little off, within tolerance: solution
    b_eq = np.array([10 + 0.5 * tol, 5 - 0.5 * tol])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_equality_constraints(lpp, tol)
    assert_(problem_feasible)

    # b_eq a little more off, outside tolerance: solution
    b_eq = np.array([10 - 1.5 * tol, 5 + 0.5 * tol])
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_fixed_variables(lp._replace(b_eq=b_eq), status, tol)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    problem_feasible = _presolve_infeasible_equality_constraints(lpp, tol)
    assert_(not problem_feasible)


def test_remove_equation_singleton_rows():
    """
    Test removal singletons from equations (A_eq, b_eq).
    """
    tol = 1e-2
    lp = _LPProblem(
        c=np.array([1, 2, 3, 4, 5], dtype=float),
        A_ub=np.empty((0, 5), dtype=float),
        b_ub=np.empty((0), dtype=float),
        A_eq=np.array([
            [ 1,  3,  5,  7,  9],
            [ 0,  0,  0, -3,  0],
            [-2, -2,  1,  8,  1],
            [ 0,  2,  0,  2,  0],
            [ 0,  0,  0, -1,  0],
            [ 5,  0,  0,  0,  0]
            ], dtype=float),
        b_eq=np.array([8, 6, -2, -3, 2 + 0.5 * tol, 4], dtype=float),
        bounds=np.array([[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], dtype=float),
        x0=None,
    )

    # Row singletons at line row 1, 4 and 5
    # Two singletons refer to the same variable
    # -3 * x[3] = 6
    # -1 * x[3] = 2 + 0.5 * tol
    # 5 * x[0] = 4

    # Solution within tolerance
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_equation_row_singletons(lp, status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.c == np.array([2, 3, 5], dtype=float)))
    assert_(np.allclose(lpp.b_eq, np.array([21.2, 15.6, 1.], dtype=float)))
    assert_(np.allclose(lpp.A_eq, np.array([[3, 5, 9], [-2, 1, 1], [2, 0, 0]], dtype=float)))

    # Do not use an integer array to insert into:
    # the resulting array will also be integer.
    xr = rev(np.array([0., 0., 0.], dtype=float))
    assert_(np.allclose(xr, np.array([0.8, 0, 0, -2, 0], dtype=float)))

    # Solution outside tolerance
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_eq = np.array([8, 6, -2, -3, 2 + 1.5 * tol, 4], dtype=float)
    lpp, rev, status = _presolve_remove_equation_row_singletons(lp._replace(b_eq=b_eq), status, tol)
    assert_(not status['solved'])
    assert_(not status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])


def test_remove_inequality_singleton_rows():
    """
    Test removal singletons from inequalities (A_ub, b_ub).
    """
    tol = 1e-2
    lp = _LPProblem(
        c=np.array([1, 2, 3, 4, 5], dtype=float),
        A_ub=np.array([
            [ 1,  3,  5,  7,  9],
            [ 0,  0,  0, -3,  0],
            [-2, -2,  1,  8,  1],
            [ 0,  2,  0,  2,  0],
            [ 0,  0,  0, -1,  0],
            [ 5,  0,  0,  0,  0]
            ], dtype=float),
        b_ub=np.empty((0), dtype=float),  # values supplied below
        A_eq=np.empty((0, 5), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.array([[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], dtype=float),
        x0=None,
    )

    # Row singletons at line row 1, 4 and 5
    # Two singletons refer to the same variable

    # Bounds unchanged
    # -3 * x[3] <= 6, x[3] >= -2
    # -1 * x[3] <= 2, x[3] >= -2
    # 5 * x[0] <= 6, x[0] <= 1.2
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, 6, -2, -3, 2, 6], dtype=float)
    lpp, rev, status = _presolve_remove_inequality_row_singletons(lp._replace(b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.allclose(lpp.bounds, np.array([[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], dtype=float)))

    # Bounds reduced
    # -3 * x[3] <= 6, x[3] >= -2
    # -1 * x[3] <= 2, x[3] >= -2
    # 5 * x[0] <= 4, x[0] <= 0.8
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, 6, -2, -3, 2, 4], dtype=float)
    lpp, rev, status = _presolve_remove_inequality_row_singletons(lp._replace(b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.allclose(lpp.bounds, np.array([[0, 0.8], [0, 1], [0, 1], [0, 1], [0, 1]], dtype=float)))

    # Bounds fixed
    # -3 * x[3] <= -3, x[3] >= 1
    # -1 * x[3] <= 2, x[3] >= -2
    # 5 * x[0] <= 0, x[0] <= 0
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, -3, -2, -3, 2, 0], dtype=float)
    lpp, rev, status = _presolve_remove_inequality_row_singletons(lp._replace(b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.allclose(lpp.bounds, np.array([[0, 0], [0, 1], [0, 1], [1, 1], [0, 1]], dtype=float)))

    # Problem infeasible
    # -3 * x[3] <= -6, x[3] >= 2
    # -1 * x[3] <= 2, x[3] >= -2
    # 5 * x[0] <= 0, x[0] <= 0
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, -6, -2, -3, 2, 0], dtype=float)
    lpp, rev, status = _presolve_remove_inequality_row_singletons(lp._replace(b_ub=b_ub), status, tol)
    print(lpp.bounds)
    assert_(not status['solved'])
    assert_(not status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])

def test_remove_empty_rows():
    """
    Test empty row removal procedure.
    """
    tol = 1e-2
    A = np.array([
        [ 1,  3,  5,  7,  9],
        [ 0,  0,  0,  0,  0],
        [-2, -2,  1,  8,  1],
        [ 0,  2,  0,  2,  0],
        [ 0,  0,  0,  0,  0],
        [ 5,  0,  0,  0,  0]
        ], dtype=float)
    lp = _LPProblem(
        c=np.zeros((5), dtype=float),  # not relevant
        A_ub=np.empty((0, 5), dtype=float),
        b_ub=np.empty((0), dtype=float),
        A_eq=np.empty((0, 5), dtype=float),
        b_eq=np.empty((0), dtype=float),
        bounds=np.empty((0, 2), dtype=float),
        x0=None,
    )

    # Empty equations: no action
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_empty_rows(lp, status, tol)
    assert_(not status['solved'])  # There is no check on number of variables remaining
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    assert_(rev is None)

    # A_eq has zero rows at index 1 and 4
    # Corresponding elements in b_eq are not 0, so problem is infeasible
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_eq = np.array([8, 6, -2, -3, 2, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_eq=A, b_eq=b_eq), status, tol)
    assert_(not status['solved'])
    assert_(not status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])

    # A_eq has zero rows at index 1 and 4
    # Corresponding elements in b_eq are 0, so rows can be removed
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_eq = np.array([8, 0, -2, -3, 0, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_eq=A, b_eq=b_eq), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.A_eq.shape == (4, 5)))
    assert_(rev is None)

    # A_eq has zero rows at index 1 and 4
    # Corresponding elements in b_eq are 0 within tolerance, so rows can be removed
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_eq = np.array([8, 0.5 * tol, -2, -3, -0.9 * tol, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_eq=A, b_eq=b_eq), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.A_eq.shape == (4, 5)))
    assert_(rev is None)

    # A_eq has zero rows at index 1 and 4
    # Corresponding elements in b_eq are not 0 within tolerance, so problem is infeasible
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_eq = np.array([8, 1.5 * tol, -2, -3, -0.9 * tol, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_eq=A, b_eq=b_eq), status, tol)
    assert_(not status['solved'])
    assert_(not status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])

    # A_ub, b_ub
    # A_ub has zero rows at index 1 and 4
    # Corresponding elements in b_ub are positive, so rows can be removed
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, 6, -2, -3, 2, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_ub=A, b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.A_ub.shape == (4, 5)))
    assert_(rev is None)

    # A_ub has zero rows at index 1 and 4
    # Corresponding elements in b_ub are 0, so rows can be removed
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, 0, -2, -3, 0, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_ub=A, b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.A_ub.shape == (4, 5)))
    assert_(rev is None)

    # A_ub has zero rows at index 1 and 4
    # Corresponding elements in b_ub are 0 within tolerance, so rows can be removed
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, 0.5 * tol, -2, -3, -0.9 * tol, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_ub=A, b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp.A_ub.shape == (4, 5)))
    assert_(rev is None)

    # A_ub has zero rows at index 1 and 4
    # Corresponding elements in b_ub are not >= 0 within tolerance, so problem
    # is infeasible
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    b_ub = np.array([8, -1.5 * tol, -2, -3, -0.9 * tol, 4])
    lpp, rev, status = _presolve_remove_empty_rows(lp._replace(A_ub=A, b_ub=b_ub), status, tol)
    assert_(not status['solved'])
    assert_(not status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])


def test_remove_empty_columns():
    """
    Test empty column removal procedure.
    """
    lp = _LPProblem(
        c=np.array([1, -2, 3, -4, 5], dtype=float),
        A_ub=np.empty((0, 5), dtype=float),  # supplied below
        b_ub=np.empty((0), dtype=float),
        A_eq=np.empty((0, 5), dtype=float),  # supplied below
        b_eq=np.empty((0), dtype=float),
        bounds=np.array([[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], dtype=float),
        x0=None,
    )

    # One zero column
    A1 = np.array([
        [ 1,  3,  0,  7,  0],
        [ 1, -1,  0,  0,  0],
        [-2, -2,  0,  8,  0],
        [ 0,  2,  0,  2,  0],
        [ 0,  1,  0, -3,  2],
        [ 5,  0,  0,  0,  0]
        ], dtype=float)

    # Two zero columns
    A2 = np.array([
        [ 3,  3,  0,  0,  7],
        [ 1, -3,  0,  0,  0],
        [-2, -1,  0,  0,  2],
        ], dtype=float)

    # Empty equations: all variables will be removed and set to lower bound
    # for positive c and upper bound for negative c
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_empty_columns(lp, status)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    xr = rev(np.array([]))
    assert_(np.allclose(xr, np.array([0., 1., 0., 1., 0.])))

    # Empty equations with 0 in c
    # Value of variable x[4] does not affect minimum, arbitrary set midway interval
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    c = np.array([1, -2, 3, -4, 0], dtype=float)
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(c=c), status)
    assert_(status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(not status['loop'])
    xr = rev(np.array([]))
    assert_(np.allclose(xr, np.array([0., 1., 0., 1., 0.5])))

    # Zero column in A_eq
    # x[2] fixed on lower bound
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A1), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Zero column in A_eq
    # x[2] fixed on lower bound
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [0, np.inf], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A1, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Zero column in A_eq
    # x[2] lower bound = -inf, so problem unbounded
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [-np.inf, 1], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A1, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(not status['bounded'])
    assert_(not status['loop'])

    # Two zero columns in A_eq
    # c[2] > 0, x[2] lower bound = 0, OK
    # c[3] < 0, x[3] upper bound = inf, so unbounded
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [0, 1], [0, np.inf], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A2, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(not status['bounded'])
    assert_(not status['loop'])

    # Two zero columns in A_eq
    # c[2] > 0, x[2] lower bound = 0, OK
    # c[3] < 0, x[3] upper bound = 1, OK
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [0, 1], [-np.inf, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A2, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, 5])))

    xr = rev(np.array([0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 1., 0.])))

    # Columns in A_ub
    # x[2] fixed on lower bound
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_ub=A1), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Columns in A_ub
    # x[2] fixed on lower bound
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [0, np.inf], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_ub=A1, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Zero column in A_ub
    # x[2] lower bound = -inf, so problem unbounded
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [-np.inf, 1], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_ub=A1, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(not status['bounded'])
    assert_(not status['loop'])

    # Zero columns in A_eq and A_ub
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A1, A_ub=A2), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Zero columns in A_eq and A_ub
    # only column 2 occurs in both
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A2, A_ub=A1), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Zero columns in A_eq and A_ub
    # x[2] fixed on lower bound
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [0, np.inf], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A1, A_ub=A2, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(status['bounded'])
    assert_(status['loop'])
    assert_(np.all(lpp[0] == np.array([1, -2, -4, 5])))

    xr = rev(np.array([0., 0., 0., 0.]))  # arbitrary values
    assert_(np.allclose(xr, np.array([0., 0., 0., 0., 0.])))

    # Zero columns in A_eq and A_ub
    # x[2] lower bound = -inf, so problem unbounded
    status = {'solved': False, 'feasible': True, 'bounded': True, 'loop': False}
    bounds = np.array([[0, 1], [0, 1], [-np.inf, 1], [0, 1], [0, 1]])
    lpp, rev, status = _presolve_remove_empty_columns(lp._replace(A_eq=A1, A_ub=A2, bounds=bounds), status)
    assert_(not status['solved'])
    assert_(status['feasible'])
    assert_(not status['bounded'])
    assert_(not status['loop'])
