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

from collections import defaultdict
import numpy as np
from .optimize import _check_unknown_options, OptimizeWarning
from warnings import warn
from ._highs.highs_wrapper import highs_wrapper
from ._highs.constants import (
    CONST_I_INF,
    CONST_INF,
    MESSAGE_LEVEL_MINIMAL,

    MODEL_STATUS_NOTSET,
    MODEL_STATUS_LOAD_ERROR,
    MODEL_STATUS_MODEL_ERROR,
    MODEL_STATUS_MODEL_EMPTY,
    MODEL_STATUS_PRESOLVE_ERROR,
    MODEL_STATUS_SOLVE_ERROR,
    MODEL_STATUS_POSTSOLVE_ERROR,
    MODEL_STATUS_PRIMAL_INFEASIBLE,
    MODEL_STATUS_PRIMAL_UNBOUNDED,
    MODEL_STATUS_OPTIMAL,
    MODEL_STATUS_REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND
    as MODEL_STATUS_RDOVUB,
    MODEL_STATUS_REACHED_TIME_LIMIT,
    MODEL_STATUS_REACHED_ITERATION_LIMIT,

    HIGHS_SIMPLEX_STRATEGY_CHOOSE,
    HIGHS_SIMPLEX_STRATEGY_DUAL,
    HIGHS_SIMPLEX_STRATEGY_PRIMAL,

    HIGHS_SIMPLEX_CRASH_STRATEGY_OFF,
    HIGHS_SIMPLEX_CRASH_STRATEGY_BIXBY,
    HIGHS_SIMPLEX_CRASH_STRATEGY_LTSF,

    HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
    HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX,
    HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE,

    HIGHS_SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
    HIGHS_SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX,
)
from scipy.sparse import csc_matrix, vstack, issparse


def _replace_inf(x):
    # Replace `np.inf` with CONST_INF
    infs = np.isinf(x)
    x[infs] = np.sign(x[infs])*CONST_INF
    return x


def _convert_to_highs_enum(option, option_str, choices, default):
    # None indicates that we want the default value:
    if option is None:
        return default
    # If option is in the choices we can look it up, if not use
    # the default value and warn:
    try:
        return choices[option.lower()]
    except (KeyError, AttributeError):
        warn("Option {0} is {1}, but only values in {2} are allowed."
             " Using default: {3}".format(
                 option_str, option, set(choices.keys()), default),
             OptimizeWarning, stacklevel=3)
        return default


def _linprog_highs(lp, solver, time_limit=None, presolve=True,
                   disp=False, maxiter=None,
                   dual_feasibility_tolerance=None,
                   dual_objective_value_upper_bound=None,
                   ipm_optimality_tolerance=None,
                   primal_feasibility_tolerance=None,
                   simplex_crash_strategy=None,
                   simplex_dual_edge_weight_strategy=None,
                   simplex_primal_edge_weight_strategy=None,
                   simplex_strategy=None,
                   simplex_update_limit=None,
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
    solver : "ipm" or "simplex" or None
        Which HiGHS solver to use.  If ``None``, HiGHS will determine which
        solver to use based on the problem.

    Options
    -------
    maxiter : int
        The maximum number of iterations to perform in either phase. For
        ``solver='ipm'``, this does not include the number of crossover
        iterations.  Default is the largest possible value for an ``int``
        on the platform.
    disp : bool
        Set to ``True`` if indicators of optimization status are to be printed
        to the console each iteration; default ``False``.
    time_limit : float
        The maximum time in seconds allotted to solve the problem; default is
        the largest possible value for a ``double`` on the platform.
    presolve : bool
        Presolve attempts to identify trivial infeasibilities,
        identify trivial unboundedness, and simplify the problem before
        sending it to the main solver. It is generally recommended
        to keep the default setting ``True``; set to ``False`` if presolve is
        to be disabled.
    dual_feasibility_tolerance : double
        Dual feasibility tolerance.  Default is 1e-07.
        The minimum of this and ``primal_feasibility_tolerance``
        is used for the feasibility tolerance when ``solver='ipm'``.
    dual_objective_value_upper_bound : double
        Upper bound on objective value for dual simplex:
        algorithm terminates if reached.  Default is the largest possible
        value for a ``double`` on the platform.
        When ``solver='ipm'`` this value is ignored.
    ipm_optimality_tolerance : double
        Optimality tolerance for ``solver='ipm'``.  Default is 1e-08.
        Minimum possible value is 1e-12 and must be smaller than the largest
        possible value for a ``double`` on the platform.
    primal_feasibility_tolerance : double
        Primal feasibility tolerance.  Default is 1e-07.
        The minimum of this and ``dual_feasibility_tolerance``
        is used for the feasibility tolerance when ``solver='ipm'``.
    simplex_crash_strategy : str {'off', 'ltsf', 'bixby'}
        Strategy for simplex crash. ``'off'`` uses no simplex crash routine.
        ``'bixby'`` is based on the CPLEX simplex crash strategy from [2]_.
        ``'ltsf'`` is based on the LTSF simplex crash strategy found in [3]_
        Default is ``off``.
    simplex_dual_edge_weight_strategy : str {'dantzig', 'devex', 'steepest'}
        Strategy for simplex dual edge weights. ``'dantzig'`` uses Dantzig’s
        original strategy of choosing the most negative reduced cost.
        ``'devex'`` uses the strategy described in [4]_.  ``steepest`` uses
        the exact steepest edge strategy as described in [5]_.
        Default is ``'devex'``.
    simplex_primal_edge_weight_strategy : str {'dantzig', 'devex'}
        Strategy for simplex primal edge weights. These options correspond
        to those of the same for ``simplex_dual_edge_weight_strategy``.
        Default is ``devex``.
    simplex_strategy : str {'choose', 'dual', 'primal'}
        Strategy for simplex solver. ``'choose'`` chooses the strategy based
        on characteristics of the problem. ``'dual'`` uses the dual revised
        simplex method while ``'primal'`` uses the primal revised simplex
        method. Default: ``choose``.
    simplex_update_limit : int
        Limit on the number of updates made to the representation of the
        basis matrix inverse (i.e. basis matrix factorization) before a
        new representation is formed from scratch.
        If needed for efficiency or numerical stability, a new
        representation of the inverse may be formed before this limit is
        reached. Default is ``5000``.
    unknown_options : dict
        Optional arguments not used by this particular solver. If
        ``unknown_options`` is non-empty, a warning is issued listing all
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
                The (nominally positive) values of the slack,
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

                ``1`` : Iteration or time limit reached.

                ``2`` : Problem appears to be infeasible.

                ``3`` : Problem appears to be unbounded.

                ``4`` : The HiGHS solver ran into a problem.

            nit : int
                The total number of iterations performed.
                For ``solver='simplex'``, this includes iterations in all
                phases. For ``solver='ipm'``, this does not include
                crossover iterations.
            crossover_nit : int
                The number of primal/dual pushes performed during the
                crossover routine for ``solver='ipm'``.  This is ``0``
                for ``solver='simplex'``.
            message : str
                A string descriptor of the exit status of the algorithm.

    References
    ----------
    .. [2] Bixby, Robert E. "Implementing the simplex method: The initial
           basis." ORSA Journal on Computing 4.3 (1992): 267-284.
    .. [3] Maros, István, and Gautam Mitra. "Strategies for creating advanced
           bases for large-scale linear programming problems." INFORMS Journal
           on Computing 10.2 (1998): 248-260.
    .. [4] Harris, Paula MJ. "Pivot selection methods of the Devex LP code."
           Mathematical programming 5.1 (1973): 1-28.
    .. [5] Goldfarb, Donald, and John Ker Reid. "A practicable steepest-edge
           simplex algorithm." Mathematical Programming 12.1 (1977): 361-371.
    """

    _check_unknown_options(unknown_options)

    # Map options to HiGHS enum values
    simplex_strategy = _convert_to_highs_enum(
        simplex_strategy,
        'simplex_strategy',
        choices={'choose': HIGHS_SIMPLEX_STRATEGY_CHOOSE,
                 'dual': HIGHS_SIMPLEX_STRATEGY_DUAL,
                 'primal': HIGHS_SIMPLEX_STRATEGY_PRIMAL},
        default=HIGHS_SIMPLEX_STRATEGY_CHOOSE)

    simplex_crash_strategy = _convert_to_highs_enum(
        simplex_crash_strategy,
        'simplex_crash_strategy',
        choices={'off': HIGHS_SIMPLEX_CRASH_STRATEGY_OFF,
                 'bixby': HIGHS_SIMPLEX_CRASH_STRATEGY_BIXBY,
                 'ltsf': HIGHS_SIMPLEX_CRASH_STRATEGY_LTSF},
        default=HIGHS_SIMPLEX_CRASH_STRATEGY_OFF)

    simplex_dual_edge_weight_strategy = _convert_to_highs_enum(
        simplex_dual_edge_weight_strategy,
        'simplex_dual_edge_weight_strategy',
        choices={'dantzig': HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
                 'devex': HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX,
                 'steepest':
                 HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE},
        default=HIGHS_SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX)

    simplex_primal_edge_weight_strategy = _convert_to_highs_enum(
        simplex_primal_edge_weight_strategy,
        'simplex_primal_edge_weight_strategy',
        choices={'dantzig': HIGHS_SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
                 'devex': HIGHS_SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX},
        default=HIGHS_SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX)

    statuses = {
        MODEL_STATUS_NOTSET: (
            4,
            'HiGHS Status Code 0: HighsModelStatusNOTSET',
        ),
        MODEL_STATUS_LOAD_ERROR: (
            4,
            'HiGHS Status Code 1: HighsModelStatusLOAD_ERROR',
        ),
        MODEL_STATUS_MODEL_ERROR: (
            2,
            'HiGHS Status Code 2: HighsModelStatusMODEL_ERROR',
        ),
        MODEL_STATUS_MODEL_EMPTY: (
            4,
            'HiGHS Status Code 3: HighsModelStatusMODEL_EMPTY',
        ),
        MODEL_STATUS_PRESOLVE_ERROR: (
            4,
            'HiGHS Status Code 4: HighsModelStatusPRESOLVE_ERROR',
        ),
        MODEL_STATUS_SOLVE_ERROR: (
            4,
            'HiGHS Status Code 5: HighsModelStatusSOLVE_ERROR',
        ),
        MODEL_STATUS_POSTSOLVE_ERROR: (
            4,
            'HiGHS Status Code 6: HighsModelStatusPOSTSOLVE_ERROR',
        ),
        MODEL_STATUS_RDOVUB: (
            4,
            'HiGHS Status Code 10: '
            'HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND',
        ),
        MODEL_STATUS_PRIMAL_INFEASIBLE: (
            2,
            "The problem is infeasible.",
        ),
        MODEL_STATUS_PRIMAL_UNBOUNDED: (
            3,
            "The problem is unbounded.",
        ),
        MODEL_STATUS_OPTIMAL: (
            0,
            "Optimization terminated successfully.",
        ),
        MODEL_STATUS_REACHED_TIME_LIMIT: (
            1,
            "Time limit reached.",
        ),
        MODEL_STATUS_REACHED_ITERATION_LIMIT: (
            1,
            "Iteration limit reached.",
        ),
    }

    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = lp

    lb, ub = bounds.T.copy()  # separate bounds, copy->C-cntgs
    # highs_wrapper solves LHS <= A*x <= RHS, not equality constraints
    lhs_ub = -np.ones_like(b_ub)*np.inf  # LHS of UB constraints is -inf
    rhs_ub = b_ub  # RHS of UB constraints is b_ub
    lhs_eq = b_eq  # Equality constaint is inequality
    rhs_eq = b_eq  # constraint with LHS=RHS
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
        'time_limit': time_limit,
        'message_level': MESSAGE_LEVEL_MINIMAL * disp,
        'dual_feasibility_tolerance': dual_feasibility_tolerance,
        'dual_objective_value_upper_bound': dual_objective_value_upper_bound,
        'ipm_optimality_tolerance': ipm_optimality_tolerance,
        'primal_feasibility_tolerance': primal_feasibility_tolerance,
        'simplex_crash_strategy': simplex_crash_strategy,
        'simplex_dual_edge_weight_strategy': simplex_dual_edge_weight_strategy,
        'simplex_primal_edge_weight_strategy':
            simplex_primal_edge_weight_strategy,
        'simplex_strategy': simplex_strategy,
        'simplex_update_limit': simplex_update_limit,
    }

    options['ipm_iteration_limit'] = maxiter
    options['simplex_iteration_limit'] = maxiter

    # np.inf doesn't work; use very large constant
    rhs = _replace_inf(rhs)
    lhs = _replace_inf(lhs)
    lb = _replace_inf(lb)
    ub = _replace_inf(ub)

    res = highs_wrapper(c, A.indptr, A.indices, A.data, lhs, rhs,
                        lb, ub, options)

    # If we get a MODEL_ERROR, match behavior of other linprog implementations
    # and return infeasible status
    # if res['status'] == MODEL_STATUS_MODEL_ERROR:
    #     res['status'] == MODEL_STATUS_PRIMAL_INFEASIBLE

    # HiGHS represents constraints as lhs/rhs, so
    # Ax + s = b => Ax = b - s
    # and we need to split up s by A_ub and A_eq
    if 'slack' in res:
        slack = res['slack']
        con = np.array(slack[len(b_ub):])
        slack = np.array(slack[:len(b_ub)])
    else:
        slack, con = None, None

    sol = {'x': np.array(res['x']) if 'x' in res else None,
           'slack': slack,
           # TODO: Add/test dual info like:
           # 'lambda': res.get('lambda'),
           # 's': res.get('s'),
           'fun': res.get('fun'),
           'con': con,
           'status': statuses[res['status']][0],
           'success': res['status'] == MODEL_STATUS_OPTIMAL,
           'message': statuses[res['status']][1],
           'nit': res.get('simplex_nit', 0) or res.get('ipm_nit', 0),
           'crossover_nit': res.get('crossover_nit'),
           }
    return sol
