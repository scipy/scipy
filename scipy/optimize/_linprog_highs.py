"""HiGHS Linear Optimization Methods

Interface to HiGHS linear optimization software.
https://www.maths.ed.ac.uk/hall/HiGHS/

.. versionadded:: 1.5.0

References
----------
.. [1] Q. Huangfu and J.A.J. Hall. "Parallelizing the dual revised simplex
           method." Mathematical Programming Computation, 10 (1), 119-142,
           2018. DOI: 10.1007/s12532-017-0130-5

"""

import numpy as np
from .optimize import _check_unknown_options
from ._highs.highs_wrapper import highs_wrapper, CONST_INF
from scipy.sparse import csc_matrix, vstack, issparse


def _replace_inf(x):
    # Replace `np.inf` with CONST_INF
    infs = (np.abs(x) == np.inf)
    x[infs] = np.sign(x[infs])*CONST_INF
    return x


def _linprog_highs(lp, solver, time_limit=1, presolve=True, parallel=False,
                   disp=False, maxiter=None, autoscale=False,
                   **unknown_options):
    r"""
    Solve the following linear programming problem using one of the HiGHS
    solvers:

    .. math::

        \min_x \ & c^T x \\
        \mbox{such that} \ & A_{ub} x \leq b_{ub},\\
        & A_{eq} x = b_{eq},\\
        & l \leq x \leq u ,

    where :math:`x` is a vector of decision variables; :math:`c`,
    :math:`b_{ub}`, :math:`b_{eq}`, :math:`l`, and :math:`u` are vectors; and
    :math:`A_{ub}` and :math:`A_{eq}` are matrices.

    Informally, that's:

    minimize::

        c @ x

    such that::

        A_ub @ x <= b_ub
        A_eq @ x == b_eq
        lb <= x <= ub

    Note that by default ``lb = 0`` and ``ub = None`` unless specified with
    ``bounds``.

    Parameters
    ----------
    lp :  _LPProblem
        A ``scipy.optimize._linprog_util._LPProblem`` ``namedtuple``.
    solver : "ipm" or "simplex"
        Which HiGHS solver to use.

    Options
    -------
    maxiter : int
        The maximum number of iterations to perform in either phase.
    disp : bool
        Set to ``True`` if indicators of optimization status are to be printed
        to the console each iteration; default ``False``.
    time_limit : float
        The maximum time allotted to solve the problem; default 1s.
    presolve : bool
        Set to ``False`` if presolve is to be disabled; default ``True``.
    parallel : bool
        Set to ``True`` to enable parallelization; default ``False``.
    autoscale : bool
        Set to ``True`` to automatically perform equilibration.
        Consider using this option if the numerical values in the
        constraints are separated by several orders of magnitude.
        Default: ``False``.
    unknown_options : dict
        Optional arguments not used by this particular solver. If
        `unknown_options` is non-empty, a warning is issued listing all
        unused options.

    Returns
    -------
    sol : dict
        A dictionary consisting of the fields:

            x : 1D array
                The values of the decision variables that minimizes the
                objective function while satisfying the constraints.
            fun : float
                The optimal value of the objective function ``c @ x``.
            slack : 1D array
                The (nominally positive) values of the slack variables,
                ``b_ub - A_ub @ x``.
            con : 1D array
                The (nominally zero) residuals of the equality constraints,
                ``b_eq - A_eq @ x``.
            success : bool
                ``True`` when the algorithm succeeds in finding an optimal
                solution.
            status : int
                An integer representing the exit status of the algorithm.

                ``0`` : Optimization terminated successfully.

                ``1`` : Iteration limit reached.

                ``2`` : Problem appears to be infeasible.

                ``3`` : Problem appears to be unbounded.

                ``4`` : Numerical difficulties encountered.

            nit : int
                The total number of iterations performed in all phases.
            message : str
                A string descriptor of the exit status of the algorithm.

    """

    _check_unknown_options(unknown_options)

    statuses = {
        12: (1, "Iteration limit reached."),
        9: (0, "Optimization terminated successfully."),
        8: (3, "The problem is unbounded."),
        7: (2, "The problem is infeasible."),
        12: (1, "Iteration limit reached."),
    }
    # "Iteration limit reached.",
    # "The problem appears infeasible, as the phase one auxiliary "
    # "problem terminated successfully with a residual of {0:.1e}, "
    # "greater than the tolerance {1} required for the solution to "
    # "be considered feasible. Consider increasing the tolerance to "
    # "be greater than {0:.1e}. If this tolerance is unnaceptably "
    # "large, the problem is likely infeasible.",
    # "The problem is unbounded, as the simplex algorithm found "
    # "a basic feasible solution from which there is a direction "
    # "with negative reduced cost in which all decision variables "
    # "increase.",
    # "Numerical difficulties encountered; consider trying "
    # "method='interior-point'.",

    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = lp

    lb, ub = bounds.T.copy()                # separate bounds, copy->C-cntgs
    # highs_wrapper solves LHS <= A*x <= RHS, not equality constraints
    lhs_ub = -np.ones_like(b_ub)*np.inf     # LHS of UB constraints is -inf
    rhs_ub = b_ub                           # RHS of UB constraints is b_ub
    lhs_eq = b_eq                           # Equality constaint is inequality
    rhs_eq = b_eq                           # constraint with LHS=RHS
    lhs = np.concatenate((lhs_ub, lhs_eq))
    rhs = np.concatenate((rhs_ub, rhs_eq))

    if issparse(A_ub) or issparse(A_eq):
        A = vstack((A_ub, A_eq))
    else:
        A = np.vstack((A_ub, A_eq))
    A = csc_matrix(A)

    options = {
        'presolve': presolve,
        'sense': 1,  # minimization
        'solver': solver,
        'parallel': parallel,
        'time_limit': time_limit,
        'message_level': disp * 1,  # 0 is none, 1 is some. Simplex only?
        'write_solution_to_file': False,
        'solution_file': 'test.sol',
        'write_solution_pretty': True,
    }

    # this doesn't seem to work, and I don't see an equivalent for IPM
    if maxiter is not None and solver == 'simplex':
        options['simplex_iteration_limit'] = maxiter

    rhs = _replace_inf(rhs)
    lhs = _replace_inf(lhs)
    lb = _replace_inf(lb)
    ub = _replace_inf(ub)

    res = highs_wrapper(c, A.indptr, A.indices, A.data, lhs, rhs,
                        lb, ub, options)

    sol = {'x': res.get('x', None),
           'slack': res.get('slack', None),
           'fun': res.get('fun', None),
           'con': res.get('con', None),
           'status': statuses[res['status']][0],
           'success': res['status'] == 9,
           'message': res['message'],
           'nit': (res.get('simplex_nit', 0) if solver == 'simplex'
                   else res.get('ipm_nit', 0))
           }
    if sol['x'] is not None:
        sol['x'] = np.array(sol['x'])
    return sol
