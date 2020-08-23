# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 19:49:17 2020

@author: matth
"""


def _linprog_highs_doc(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None,
                       bounds=None, method='highs', callback=None,
                       maxiter=None, disp=False, presolve=True,
                       time_limit=None,
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
    Linear programming: minimize a linear objective function subject to linear
    equality and inequality constraints using one of the HiGHS solvers.

    Linear programming solves problems of the following form:

    .. math::

        \min_x \ & c^T x \\
        \mbox{such that} \ & A_{ub} x \leq b_{ub},\\
        & A_{eq} x = b_{eq},\\
        & l \leq x \leq u ,

    where :math:`x` is a vector of decision variables; :math:`c`,
    :math:`b_{ub}`, :math:`b_{eq}`, :math:`l`, and :math:`u` are vectors; and
    :math:`A_{ub}` and :math:`A_{eq}` are matrices.

    Alternatively, that's:

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
    c : 1-D array
        The coefficients of the linear objective function to be minimized.
    A_ub : 2-D array, optional
        The inequality constraint matrix. Each row of ``A_ub`` specifies the
        coefficients of a linear inequality constraint on ``x``.
    b_ub : 1-D array, optional
        The inequality constraint vector. Each element represents an
        upper bound on the corresponding value of ``A_ub @ x``.
    A_eq : 2-D array, optional
        The equality constraint matrix. Each row of ``A_eq`` specifies the
        coefficients of a linear equality constraint on ``x``.
    b_eq : 1-D array, optional
        The equality constraint vector. Each element of ``A_eq @ x`` must equal
        the corresponding element of ``b_eq``.
    bounds : sequence, optional
        A sequence of ``(min, max)`` pairs for each element in ``x``, defining
        the minimum and maximum values of that decision variable. Use ``None``
        to indicate that there is no bound. By default, bounds are
        ``(0, None)`` (all decision variables are non-negative).
        If a single tuple ``(min, max)`` is provided, then ``min`` and
        ``max`` will serve as bounds for all decision variables.
    method : {'highs-simplex', 'highs-ipm', 'highs'}, 'interior-point',
    'revised simplex', 'simplex'}, optional

        The algorithm used to solve the standard form problem.
        :ref:`'highs-simplex' <optimize.linprog-highs-simplex>`,
        :ref:`'highs-ipm' <optimize.linprog-highs-ipm>`,
        :ref:`'highs' <optimize.linprog-highs>`,
        :ref:`'interior-point' <optimize.linprog-interior-point>` (default),
        :ref:`'revised simplex' <optimize.linprog-revised_simplex>`, and
        :ref:`'simplex' <optimize.linprog-simplex>` (legacy)
        are supported.

    options : dict, optional
        A dictionary of solver options. All methods accept the following
        options:

        maxiter : int
            The maximum number of iterations to perform in either phase.
            For ``solver='ipm'``, this does not include the number of
            crossover iterations.  Default is the largest possible value
            for an ``int`` on the platform.
        disp : bool (default: ``False``)
            Set to ``True`` if indicators of optimization status are to be
            printed to the console during optimization.
        presolve : bool (default: ``True``)
            Presolve attempts to identify trivial infeasibilities,
            identify trivial unboundedness, and simplify the problem before
            sending it to the main solver. It is generally recommended
            to keep the default setting ``True``; set to ``False`` if
            presolve is to be disabled.
        time_limit : float
            The maximum time in seconds allotted to solve the problem;
            default is the largest possible value for a ``double`` on the
            platform.
        dual_feasibility_tolerance : double (default: 1e-07)
            Dual feasibility tolerance.
            The minimum of this and ``primal_feasibility_tolerance``
            is used for the feasibility tolerance when ``solver='ipm'``.
        dual_objective_value_upper_bound : double
            Upper bound on objective value for dual simplex:
            algorithm terminates if reached.  Default is the largest
            possible value for a ``double`` on the platform.
            When ``solver='ipm'`` this value is ignored.
        ipm_optimality_tolerance : double
            Optimality tolerance for ``solver='ipm'``.  Default is 1e-08.
            Minimum possible value is 1e-12 and must be smaller than the
            largest possible value for a ``double`` on the platform.
        primal_feasibility_tolerance : double (default: 1e-07)
            Primal feasibility tolerance.
            The minimum of this and ``dual_feasibility_tolerance``
            is used for the feasibility tolerance when ``solver='ipm'``.
        simplex_crash_strategy : int {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
            Strategy for simplex crash: off / LTSSF / Bixby (0/1/2).
            Default is ``0``.  Corresponds to the following:

            ``0``: `SIMPLEX_CRASH_STRATEGY_OFF`

            ``1``: `SIMPLEX_CRASH_STRATEGY_LTSSF_K`

            ``2``: `SIMPLEX_CRASH_STRATEGY_BIXBY`

            ``3``: `SIMPLEX_CRASH_STRATEGY_LTSSF_PRI`

            ``4``: `SIMPLEX_CRASH_STRATEGY_LTSF_K`

            ``5``: `SIMPLEX_CRASH_STRATEGY_LTSF_PRI`

            ``6``: `SIMPLEX_CRASH_STRATEGY_LTSF`

            ``7``: `SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS`

            ``8``: `SIMPLEX_CRASH_STRATEGY_BASIC`

            ``9``: `SIMPLE_CRASH_STRATEGY_TEST_SING`

        ``SIMPLEX_CRASH_STRATEGY_*`` are defined as in HiGHS.

        simplex_dual_edge_weight_strategy : int {0, 1, 2, 3, 4}
            Strategy for simplex dual edge weights:
            Dantzig / Devex / Steepest Edge.  Default is ``2``.
            Corresponds to the following:

            ``0``: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG`

            ``1``: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX`

            ``2``: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH`

            ``3``: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE`

            ``4``: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL`

            ``SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_*`` are defined as in
            HiGHS.

        simplex_primal_edge_weight_strategy : int {0, 1}
            Strategy for simplex primal edge weights:
            Dantzig / Devex.  Default is ``0``.
            Corresponds to the following:

            ``0``: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG`

            ``1``: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX`

            ``SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_*`` are defined as in
            HiGHS.

        simplex_strategy : int {0, 1, 2, 3}
            Strategy for simplex solver. Default: ``1``.
            Corresponds to the following:

            ``0``: `SIMPLEX_STRATEGY_MIN`

            ``1``: `SIMPLEX_STRATEGY_DUAL`

            ``2``: `SIMPLEX_STRATEGY_DUAL_TASKS`

            ``3``: `SIMPLEX_STRATEGY_DUAL_MULTI`

            ``SIMPLEX_STRATEGY_*`` are defined as in HiGHS.

        simplex_update_limit : int (default: ``5000``)
            Limit on the number of updates made to the representation of
            the basis matrix inverse (e.g. LU factorization)
            before a new representation is formed from scratch.
            If needed for efficiency or numerical stability, a new
            representation of the inverse may be formed before this limit
            is reached. See [2]_  Secture 2.4 for more information about
            updating the representation of the basis matrix inverse.

        unknown_options : dict
            Optional arguments not used by this particular solver. If
            ``unknown_options`` is non-empty, a warning is issued listing
            all unused options.

    Returns
    -------
    res : OptimizeResult
        A :class:`scipy.optimize.OptimizeResult` consisting of the fields:

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

    Notes
    -----

    Method :ref:`'highs-simplex' <optimize.linprog-highs-simplex> is a wrapper
    of the C++ high performance dual revised simplex implementation (HSOL)
    [1]_, [2]_. Method :ref:`'highs-ipm' <optimize.linprog-highs-ipm>`
    is a wrapper of a C++ implementation of an **i**\ nterior-\ **p**\ oint
    **m**\ ethod [1]_; it features a crossover routine, so it is as accurate
    as a simplex solver. Method :ref:`'highs' <optimize.linprog-highs> chooses
    between the two automatically. For new code involving `linprog`, we
    recommend explicitly choosing one of these three method values instead of
    :ref:`'interior-point' <optimize.linprog-interior-point>` (default),
    :ref:`'revised simplex' <optimize.linprog-revised_simplex>`, and
    :ref:`'simplex' <optimize.linprog-simplex>` (legacy).

    References
    ----------
    .. [1] Huangfu, Q., Galabova, I., Feldmeier, M., and Hall, J. A. J.
           "HiGHS - high performance software for linear optimization."
           Accessed 4/16/2020 at https://www.maths.ed.ac.uk/hall/HiGHS/#guide
    .. [2] Huangfu, Q. and Hall, J. A. J. "Parallelizing the dual revised
           simplex method." Mathematical Programming Computation, 10 (1),
           119-142, 2018. DOI: 10.1007/s12532-017-0130-5

    """
    pass
