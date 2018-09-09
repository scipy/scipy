"""
Method agnostic utility functions for linear progamming
"""

import numpy as np


def _postsolve(x, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
               complete=False, undo=[], tol=1e-8):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1D array
        Solution vector to the standard-form problem.
    c : 1D array
        Original coefficients of the linear objective function to be minimized.
    A_ub : 2D array
        Original upper bound constraint matrix.
    b_ub : 1D array
        Original upper bound constraint vector.
    A_eq : 2D array
        Original equality constraint matrix.
    b_eq : 1D array
        Original equality constraint vector.
    bounds : sequence of tuples
        Bounds, as modified in presolve
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)
    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

    Returns
    -------
    x : 1D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack: 1D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub * x``
    con : 1D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq * x``
    """
    # note that all the inputs are the ORIGINAL, unmodified versions
    # no rows, columns have been removed
    # the only exception is bounds; it has been modified
    # we need these modified values to undo the variable substitutions
    # in retrospect, perhaps this could have been simplified if the "undo"
    # variable also contained information for undoing variable substitutions

    n_x = len(c)

    # we don't have to undo variable substitutions for fixed variables that
    # were removed from the problem
    no_adjust = set()

    # if there were variables removed from the problem, add them back into the
    # solution vector
    if len(undo) > 0:
        no_adjust = set(undo[0])
        x = x.tolist()
        for i, val in zip(undo[0], undo[1]):
            x.insert(i, val)
        x = np.array(x)

    # now undo variable substitutions
    # if "complete", problem was solved in presolve; don't do anything here
    if not complete and bounds is not None:  # bounds are never none, probably
        n_unbounded = 0
        for i, b in enumerate(bounds):
            if i in no_adjust:
                continue
            lb, ub = b
            if lb is None and ub is None:
                n_unbounded += 1
                x[i] = x[i] - x[n_x + n_unbounded - 1]
            else:
                if lb is None:
                    x[i] = ub - x[i]
                else:
                    x[i] += lb

    n_x = len(c)
    x = x[:n_x]  # all the rest of the variables were artificial
    fun = x.dot(c)
    slack = b_ub - A_ub.dot(x)  # report slack for ORIGINAL UB constraints
    # report residuals of ORIGINAL EQ constraints
    con = b_eq - A_eq.dot(x)

    # Patch for bug #8664. Detecting this sort of issue earlier
    # (via abnormalities in the indicators) would be better.
    bounds = np.array(bounds)  # again, this should have been the standard form
    lb = bounds[:, 0]
    ub = bounds[:, 1]
    lb[np.equal(lb, None)] = -np.inf
    ub[np.equal(ub, None)] = np.inf

    return x, fun, slack, con, lb, ub


def _check_result(x, fun, status, slack, con, lb, ub, tol, message):
    """
    Check the validity of the provided solution.

    A valid (optimal) solution satisfies all bounds, all slack variables are
    negative and all equality constraint residuals are strictly non-zero.
    Further, the lower-bounds, upper-bounds, slack and residuals contain
    no nan values.

    Parameters
    ----------
    x : 1D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties which could not resolved using
                 a more robust, albeit less efficient, solver encountered

    slack: 1D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub * x``
    con : 1D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq * x``
    message : str
        A string descriptor of the exit status of the optimization.
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

    Returns
    -------
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties which could not resolved using
                 a more robust, albeit less efficient, solver encountered

    message : str
        A string descriptor of the exit status of the optimization.
    """
    # Somewhat arbitrary, but status 5 is very unusual
    tol = np.sqrt(tol) * 10

    contains_nans = (
        np.isnan(x).any()
        or np.isnan(fun)
        or np.isnan(slack).any()
        or np.isnan(con).any()
    )

    if contains_nans:
        is_feasible = False
    else:
        invalid_bounds = (x < lb - tol).any() or (x > ub + tol).any()
        invalid_slack = status != 3 and (slack < -tol).any()
        invalid_con = status != 3 and (np.abs(con) > tol).any()
        is_feasible = not (invalid_bounds or invalid_slack or invalid_con)

    if status == 0 and not is_feasible:
        status = 4
        message = ("The solution does not satisfy the constraints, yet "
                   "no errors were raised and there is no certificate of "
                   "infeasibility or unboundedness. This is known to occur "
                   "if the `presolve` option is False and the problem is "
                   "infeasible. If you encounter this under different "
                   "circumstances, please submit a bug report. Otherwise, "
                   "please enable presolve.")
    elif status == 0 and contains_nans:
        status = 4
        message = ("Numerical difficulties were encountered but no errors "
                   "were raised. This is known to occur if the 'presolve' "
                   "option is False, 'sparse' is True, and A_eq includes "
                   "redundant rows. If you encounter this under different "
                   "circumstances, please submit a bug report. Otherwise, "
                   "remove linearly dependent equations from your equality "
                   "constraints or enable presolve.")
    elif status == 2 and is_feasible:
        # Occurs if the simplex method exits after phase one with a very
        # nearly basic feasible solution. Postsolving can then make the
        # solution basic. The solution is NOT optimal
        raise ValueError(message)

    return status, message


def _postprocess(x, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
                 complete=False, undo=[], status=0, message="", tol=1e-8):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1D array
        Solution vector to the standard-form problem.
    c : 1D array
        Original coefficients of the linear objective function to be minimized.
    A_ub : 2D array
        Original upper bound constraint matrix.
    b_ub : 1D array
        Original upper bound constraint vector.
    A_eq : 2D array
        Original equality constraint matrix.
    b_eq : 1D array
        Original equality constraint vector.
    bounds : sequence of tuples
        Bounds, as modified in presolve
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)
    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties which could not resolved using
                 a more robust, albeit less efficient, solver encountered

    message : str
        A string descriptor of the exit status of the optimization.
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

    Returns
    -------
    x : 1D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack: 1D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub * x``
    con : 1D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq * x``
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties which could not resolved using
                 a more robust, albeit less efficient, solver encountered

    message : str
        A string descriptor of the exit status of the optimization.

    """

    x, fun, slack, con, lb, ub = _postsolve(
        x, c, A_ub, b_ub, A_eq, b_eq,
        bounds, complete, undo, tol
    )

    status, message = _check_result(
        x, fun, status, slack, con,
        lb, ub, tol, message
    )

    return x, fun, slack, con, status, message
