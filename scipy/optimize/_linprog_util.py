"""
Method agnostic utility functions for linear progamming
"""

import numpy as np
import scipy.sparse as sps
from warnings import warn
from .optimize import OptimizeWarning
from scipy.optimize._remove_redundancy import (
    _remove_redundancy, _remove_redundancy_sparse, _remove_redundancy_dense)
from collections import namedtuple

_LPProblem = namedtuple('_LPProblem', 'c A_ub b_ub A_eq b_eq bounds x0')
_LPProblem.__new__.__defaults__ = (None,) * 6  # make c the only required arg
_LPProblem.__doc__ = \
    """ Represents a linear-programming problem.

    Attributes
    ----------
    c : 1D array
        The coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        The inequality constraint matrix. Each row of ``A_ub`` specifies the
        coefficients of a linear inequality constraint on ``x``.
    b_ub : 1D array, optional
        The inequality constraint vector. Each element represents an
        upper bound on the corresponding value of ``A_ub @ x``.
    A_eq : 2D array, optional
        The equality constraint matrix. Each row of ``A_eq`` specifies the
        coefficients of a linear equality constraint on ``x``.
    b_eq : 1D array, optional
        The equality constraint vector. Each element of ``A_eq @ x`` must equal
        the corresponding element of ``b_eq``.
    bounds : various valid formats, optional
        The bounds of ``x``, as ``min`` and ``max`` pairs.
        If bounds are specified for all N variables separately, valid formats
        are:
        * a 2D array (N x 2);
        * a sequence of N sequences, each with 2 values.
        If all variables have the same bounds, the bounds can be specified as
        a 1-D or 2-D array or sequence with 2 scalar values.
        If all variables have a lower bound of 0 and no upper bound, the bounds
        parameter can be omitted (or given as None).
        Absent lower and/or upper bounds can be specified as -numpy.inf (no
        lower bound), numpy.inf (no upper bound) or None (both).
    x0 : 1D array, optional
        Guess values of the decision variables, which will be refined by
        the optimization algorithm. This argument is currently used only by the
        'revised simplex' method, and can only be used if `x0` represents a
        basic feasible solution.

    Notes
    -----
    This namedtuple supports 2 ways of initialization:
    >>> lp1 = _LPProblem(c=[-1, 4], A_ub=[[-3, 1], [1, 2]], b_ub=[6, 4])
    >>> lp2 = _LPProblem([-1, 4], [[-3, 1], [1, 2]], [6, 4])

    Note that only ``c`` is a required argument here, whereas all other arguments
    ``A_ub``, ``b_ub``, ``A_eq``, ``b_eq``, ``bounds``, ``x0`` are optional with
    default values of None.
    For example, ``A_eq`` and ``b_eq`` can be set without ``A_ub`` or ``b_ub``:
    >>> lp3 = _LPProblem(c=[-1, 4], A_eq=[[2, 1]], b_eq=[10])
    """


def _check_sparse_inputs(options, A_ub, A_eq):
    """
    Check the provided ``A_ub`` and ``A_eq`` matrices conform to the specified
    optional sparsity variables.

    Parameters
    ----------
    A_ub : 2-D array, optional
        2-D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    A_eq : 2-D array, optional
        2-D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    options : dict
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.

    Returns
    -------
    A_ub : 2-D array, optional
        2-D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    A_eq : 2-D array, optional
        2-D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    options : dict
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.
    """
    # This is an undocumented option for unit testing sparse presolve
    _sparse_presolve = options.pop('_sparse_presolve', False)
    if _sparse_presolve and A_eq is not None:
        A_eq = sps.coo_matrix(A_eq)
    if _sparse_presolve and A_ub is not None:
        A_ub = sps.coo_matrix(A_ub)

    sparse = options.get('sparse', False)
    if not sparse and (sps.issparse(A_eq) or sps.issparse(A_ub)):
        options['sparse'] = True
        warn("Sparse constraint matrix detected; setting 'sparse':True.",
             OptimizeWarning, stacklevel=4)
    return options, A_ub, A_eq


def _format_A_constraints(A, n_x, sparse_lhs=False):
    """Format the left hand side of the constraints to a 2-D array

    Parameters
    ----------
    A : 2-D array
        2-D array such that ``A @ x`` gives the values of the upper-bound
        (in)equality constraints at ``x``.
    n_x : int
        The number of variables in the linear programming problem.
    sparse_lhs : bool
        Whether either of `A_ub` or `A_eq` are sparse. If true return a
        coo_matrix instead of a numpy array.

    Returns
    -------
    np.ndarray or sparse.coo_matrix
        2-D array such that ``A @ x`` gives the values of the upper-bound
        (in)equality constraints at ``x``.

    """
    if sparse_lhs:
        return sps.coo_matrix(
            (0, n_x) if A is None else A, dtype=float, copy=True
        )
    elif A is None:
        return np.zeros((0, n_x), dtype=float)
    else:
        return np.array(A, dtype=float, copy=True)


def _format_b_constraints(b):
    """Format the upper bounds of the constraints to a 1-D array

    Parameters
    ----------
    b : 1-D array
        1-D array of values representing the upper-bound of each (in)equality
        constraint (row) in ``A``.

    Returns
    -------
    1-D np.array
        1-D array of values representing the upper-bound of each (in)equality
        constraint (row) in ``A``.

    """
    if b is None:
        return np.array([], dtype=float)
    b = np.array(b, dtype=float, copy=True).squeeze()
    return b if b.size != 1 else b.reshape((-1))


def _clean_inputs(lp):
    """
    Given user inputs for a linear programming problem, return the
    objective vector, upper bound constraints, equality constraints,
    and simple bounds in a preferred format.

    Parameters
    ----------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : various valid formats, optional
            The bounds of ``x``, as ``min`` and ``max`` pairs.
            If bounds are specified for all N variables separately, valid formats are:
            * a 2D array (2 x N or N x 2);
            * a sequence of N sequences, each with 2 values.
            If all variables have the same bounds, a single pair of values can
            be specified. Valid formats are:
            * a sequence with 2 scalar values;
            * a sequence with a single element containing 2 scalar values.
            If all variables have a lower bound of 0 and no upper bound, the bounds
            parameter can be omitted (or given as None).
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    Returns
    -------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, as ``min`` and ``max`` pairs, one for each of the N
            elements of ``x``. The N x 2 array contains lower bounds in the first
            column and upper bounds in the 2nd. Unbounded variables have lower
            bound -np.inf and/or upper bound np.inf.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    """
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = lp

    if c is None:
        raise TypeError

    try:
        c = np.array(c, dtype=np.float64, copy=True).squeeze()
    except ValueError as e:
        raise TypeError(
            "Invalid input for linprog: c must be a 1-D array of numerical "
            "coefficients") from e
    else:
        # If c is a single value, convert it to a 1-D array.
        if c.size == 1:
            c = c.reshape((-1))

        n_x = len(c)
        if n_x == 0 or len(c.shape) != 1:
            raise ValueError(
                "Invalid input for linprog: c must be a 1-D array and must "
                "not have more than one non-singleton dimension")
        if not(np.isfinite(c).all()):
            raise ValueError(
                "Invalid input for linprog: c must not contain values "
                "inf, nan, or None")

    sparse_lhs = sps.issparse(A_eq) or sps.issparse(A_ub)
    try:
        A_ub = _format_A_constraints(A_ub, n_x, sparse_lhs=sparse_lhs)
    except ValueError as e:
        raise TypeError(
            "Invalid input for linprog: A_ub must be a 2-D array "
            "of numerical values") from e
    else:
        n_ub = A_ub.shape[0]
        if len(A_ub.shape) != 2 or A_ub.shape[1] != n_x:
            raise ValueError(
                "Invalid input for linprog: A_ub must have exactly two "
                "dimensions, and the number of columns in A_ub must be "
                "equal to the size of c")
        if (sps.issparse(A_ub) and not np.isfinite(A_ub.data).all()
                or not sps.issparse(A_ub) and not np.isfinite(A_ub).all()):
            raise ValueError(
                "Invalid input for linprog: A_ub must not contain values "
                "inf, nan, or None")

    try:
        b_ub = _format_b_constraints(b_ub)
    except ValueError as e:
        raise TypeError(
            "Invalid input for linprog: b_ub must be a 1-D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_ub") from e
    else:
        if b_ub.shape != (n_ub,):
            raise ValueError(
                "Invalid input for linprog: b_ub must be a 1-D array; b_ub "
                "must not have more than one non-singleton dimension and "
                "the number of rows in A_ub must equal the number of values "
                "in b_ub")
        if not(np.isfinite(b_ub).all()):
            raise ValueError(
                "Invalid input for linprog: b_ub must not contain values "
                "inf, nan, or None")

    try:
        A_eq = _format_A_constraints(A_eq, n_x, sparse_lhs=sparse_lhs)
    except ValueError as e:
        raise TypeError(
            "Invalid input for linprog: A_eq must be a 2-D array "
            "of numerical values") from e
    else:
        n_eq = A_eq.shape[0]
        if len(A_eq.shape) != 2 or A_eq.shape[1] != n_x:
            raise ValueError(
                "Invalid input for linprog: A_eq must have exactly two "
                "dimensions, and the number of columns in A_eq must be "
                "equal to the size of c")

        if (sps.issparse(A_eq) and not np.isfinite(A_eq.data).all()
                or not sps.issparse(A_eq) and not np.isfinite(A_eq).all()):
            raise ValueError(
                "Invalid input for linprog: A_eq must not contain values "
                "inf, nan, or None")

    try:
        b_eq = _format_b_constraints(b_eq)
    except ValueError as e:
        raise TypeError(
            "Invalid input for linprog: b_eq must be a 1-D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_eq") from e
    else:
        if b_eq.shape != (n_eq,):
            raise ValueError(
                "Invalid input for linprog: b_eq must be a 1-D array; b_eq "
                "must not have more than one non-singleton dimension and "
                "the number of rows in A_eq must equal the number of values "
                "in b_eq")
        if not(np.isfinite(b_eq).all()):
            raise ValueError(
                "Invalid input for linprog: b_eq must not contain values "
                "inf, nan, or None")

    # x0 gives a (optional) starting solution to the solver. If x0 is None,
    # skip the checks. Initial solution will be generated automatically.
    if x0 is not None:
        try:
            x0 = np.array(x0, dtype=float, copy=True).squeeze()
        except ValueError as e:
            raise TypeError(
                "Invalid input for linprog: x0 must be a 1-D array of "
                "numerical coefficients") from e
        if x0.ndim == 0:
            x0 = x0.reshape((-1))
        if len(x0) == 0 or x0.ndim != 1:
            raise ValueError(
                "Invalid input for linprog: x0 should be a 1-D array; it "
                "must not have more than one non-singleton dimension")
        if not x0.size == c.size:
            raise ValueError(
                "Invalid input for linprog: x0 and c should contain the "
                "same number of elements")
        if not np.isfinite(x0).all():
            raise ValueError(
                "Invalid input for linprog: x0 must not contain values "
                "inf, nan, or None")

    # Bounds can be one of these formats:
    # (1) a 2-D array or sequence, with shape N x 2
    # (2) a 1-D or 2-D sequence or array with 2 scalars
    # (3) None (or an empty sequence or array)
    # Unspecified bounds can be represented by None or (-)np.inf.
    # All formats are converted into a N x 2 np.array with (-)np.inf where
    # bounds are unspecified.

    # Prepare clean bounds array
    bounds_clean = np.zeros((n_x, 2), dtype=float)

    # Convert to a numpy array.
    # np.array(..,dtype=float) raises an error if dimensions are inconsistent
    # or if there are invalid data types in bounds. Just add a linprog prefix
    # to the error and re-raise.
    # Creating at least a 2-D array simplifies the cases to distinguish below.
    if bounds is None or np.array_equal(bounds, []) or np.array_equal(bounds, [[]]):
        bounds = (0, np.inf)
    try:
        bounds_conv = np.atleast_2d(np.array(bounds, dtype=float))
    except ValueError as e:
        raise ValueError(
            "Invalid input for linprog: unable to interpret bounds, "
            "check values and dimensions: " + e.args[0]) from e
    except TypeError as e:
        raise TypeError(
            "Invalid input for linprog: unable to interpret bounds, "
            "check values and dimensions: " + e.args[0]) from e

    # Check bounds options
    bsh = bounds_conv.shape
    if len(bsh) > 2:
        # Do not try to handle multidimensional bounds input
        raise ValueError(
            "Invalid input for linprog: provide a 2-D array for bounds, "
            "not a {:d}-D array.".format(len(bsh)))
    elif np.all(bsh == (n_x, 2)):
        # Regular N x 2 array
        bounds_clean = bounds_conv
    elif (np.all(bsh == (2, 1)) or np.all(bsh == (1, 2))):
        # 2 values: interpret as overall lower and upper bound
        bounds_flat = bounds_conv.flatten()
        bounds_clean[:, 0] = bounds_flat[0]
        bounds_clean[:, 1] = bounds_flat[1]
    elif np.all(bsh == (2, n_x)):
        # Reject a 2 x N array
        raise ValueError(
            "Invalid input for linprog: provide a {:d} x 2 array for bounds, "
            "not a 2 x {:d} array.".format(n_x, n_x))
    else:
        raise ValueError(
            "Invalid input for linprog: unable to interpret bounds with this "
            "dimension tuple: {0}.".format(bsh))

    # The process above creates nan-s where the input specified None
    # Convert the nan-s in the 1st column to -np.inf and in the 2nd column
    # to np.inf
    i_none = np.isnan(bounds_clean[:, 0])
    bounds_clean[i_none, 0] = -np.inf
    i_none = np.isnan(bounds_clean[:, 1])
    bounds_clean[i_none, 1] = np.inf

    return _LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds_clean, x0)


def _presolve(lp, rr, tol=1e-9):
    """
    Given inputs for a linear programming problem in preferred format,
    presolve the problem: identify trivial infeasibilities, redundancies,
    and unboundedness, tighten bounds where possible, and eliminate fixed
    variables.

    Parameters
    ----------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, as ``min`` and ``max`` pairs, one for each of the N
            elements of ``x``. The N x 2 array contains lower bounds in the first
            column and upper bounds in the 2nd. Unbounded variables have lower
            bound -np.inf and/or upper bound np.inf.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    rr : bool
        If ``True`` attempts to eliminate any redundant rows in ``A_eq``.
        Set False if ``A_eq`` is known to be of full row rank, or if you are
        looking for a potential speedup (at the expense of reliability).
    tol : float
        The tolerance which determines when a solution is "close enough" to
        zero in Phase 1 to be considered a basic feasible solution or close
        enough to positive to serve as an optimal solution.

    Returns
    -------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, as ``min`` and ``max`` pairs, possibly tightened.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    c0 : 1D array
        Constant term in objective function due to fixed (and eliminated)
        variables.
    x : 1D array
        Solution vector (when the solution is trivial and can be determined
        in presolve)
    revstack: list of functions
        the functions in the list reverse the operations of _presolve()
        the function signature is x_org = f(x_mod), where x_mod is the result
        of a presolve step and x_org the value at the start of the step
        (currently, the revstack contains only one function)
    complete: bool
        Whether the solution is complete (solved or determined to be infeasible
        or unbounded in presolve)
    status : int
        An integer representing the exit status of the optimization::

         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded
         4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [5] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.
    .. [8] Andersen, Erling D., and Knud D. Andersen. "Presolving in linear
           programming." Mathematical Programming 71.2 (1995): 221-245.

    """
    # ideas from Reference [5] by Andersen and Andersen
    # however, unlike the reference, this is performed before converting
    # problem to standard form
    # There are a few advantages:
    #  * artificial variables have not been added, so matrices are smaller
    #  * bounds have not been converted to constraints yet. (It is better to
    #    do that after presolve because presolve may adjust the simple bounds.)
    # There are still improvements that can be made, namely:
    #  * remaining checks from [5] may be implemented, although they may prove
    #    to be too computationally involved
    #  * implement additional efficiency improvements in redundancy removal [2]
    #  * introcuce a presolve level (instead of just on/off)

    def _continue(p_stat):
        """
        Return True if status indicates that problem presolving should continue.
        The problem must not be solved, must be feasible, bounded and there
        must have been changes in a previous step.
        """
        return (not p_stat['solved'] and p_stat['feasible'] and p_stat['bounded'] and p_stat['loop'])

    # Presolve processing summary: initial situation
    pps = [len(lp.c), lp.A_eq.shape[0], lp.A_ub.shape[0]]

    if sps.issparse(lp.A_eq):
        A_eq = lp.A_eq.tocsr()
        A_ub = lp.A_ub.tocsr()
        lp = lp._replace(A_eq=A_eq, A_ub=A_ub)

    # Presolve status is a dictionary. Item 'loop' specifies whether the
    # presolve loop should be (re)entered.
    p_stat = {'solved': False, 'feasible': True, 'bounded': True, 'loop': True}

    # Reversal function list
    revstack = []

    # Initial feasibility check and solved check
    p_stat = _presolve_check_problem(lp, p_stat, tol)

    # Repeatedly preprocess the problem until no further changes
    while _continue(p_stat):

        # Preprocess the problem
        # This process may detect unboundedness (set 'bounded' to False).
        # If preprocessing changes the problem, 'loop' is set to True.
        p_stat['loop'] = False
        (lp, revstack, p_stat) = _presolve_preprocess_problem(lp, revstack, p_stat, tol)

        # If the problem has changed, recheck feasibility & solved
        if p_stat['loop']:
            p_stat = _presolve_check_problem(lp, p_stat, tol)

    # Presolve processing summary: final situation
    pps.extend([len(lp.c), lp.A_eq.shape[0], lp.A_ub.shape[0]])

    # Output: list containing
    # 0. number of variables
    # 1. number of variables eliminated
    # 2. number of equations
    # 3. number of equations eliminated
    # 4. number of redundant equations
    # 5. number of inequalities
    # 6. number of inequalities eliminated
    presolve_effect = [pps[0], pps[0] - pps[3], pps[1], pps[1] - pps[4],
                       0, pps[2], pps[2] - pps[5]]

    # Return the preprocessed problem with status numeral
    status = 0
    if not p_stat['bounded']:
        status = 3
    elif not p_stat['feasible']:
        status = 2
    complete = _presolve_complete(p_stat)
    message = ""
    if complete:
        return (lp, revstack, complete, status, message, presolve_effect)

    # If not complete: remove redundant (linearly dependent) rows from
    # equality constraints.
    A_eq = lp.A_eq
    b_eq = lp.b_eq
    n_rows_A = A_eq.shape[0]
    redundancy_warning = ("A_eq does not appear to be of full row rank. To "
                          "improve performance, check the problem formulation "
                          "for redundant equality constraints.")
    if (sps.issparse(A_eq)):
        if rr and A_eq.size > 0:  # TODO: Fast sparse rank check?
            A_eq, b_eq, status, message = _remove_redundancy_sparse(A_eq, b_eq)
            if A_eq.shape[0] < n_rows_A:
                warn(redundancy_warning, OptimizeWarning, stacklevel=1)
            if status != 0:
                complete = True
        presolve_effect[4] = pps[4] - A_eq.shape[0]
        return (lp._replace(A_eq=A_eq, b_eq=b_eq),
                revstack, complete, status, message, presolve_effect)

    # This is a wild guess for which redundancy removal algorithm will be
    # faster. More testing would be good.
    small_nullspace = 5
    if rr and A_eq.size > 0:
        try:  # TODO: instead use results of first SVD in _remove_redundancy
            rank = np.linalg.matrix_rank(A_eq)
        except Exception:  # oh well, we'll have to go with _remove_redundancy_dense
            rank = 0
    if rr and A_eq.size > 0 and rank < A_eq.shape[0]:
        warn(redundancy_warning, OptimizeWarning, stacklevel=3)
        dim_row_nullspace = A_eq.shape[0] - rank
        if dim_row_nullspace <= small_nullspace:
            A_eq, b_eq, status, message = _remove_redundancy(A_eq, b_eq)
        if dim_row_nullspace > small_nullspace or status == 4:
            A_eq, b_eq, status, message = _remove_redundancy_dense(A_eq, b_eq)
        if A_eq.shape[0] < rank:
            message = ("Due to numerical issues, redundant equality "
                       "constraints could not be removed automatically. "
                       "Try providing your constraint matrices as sparse "
                       "matrices to activate sparse presolve, try turning "
                       "off redundancy removal, or try turning off presolve "
                       "altogether.")
            status = 4
        if status != 0:
            complete = True
    presolve_effect[4] = pps[4] - A_eq.shape[0]
    return (lp._replace(A_eq=A_eq, b_eq=b_eq),
            revstack, complete, status, message, presolve_effect)


def _presolve_complete(p_stat):
    """
    Check if problem is complete. A problem is complete if presolve solved the
    problem, or if it detected it is not feasible or unbounded.

    Parameters
    ----------
    p_stat : scipy.optimize._linprog_util._LPProblem
        The problem parameters

    Returns
    -------
    Boolean
        True if status indicates that the problem needs no further processing.

    """
    return (p_stat['solved'] or not p_stat['feasible'] or not p_stat['bounded'])


def _presolve_check_problem(lp, p_stat, tol):
    """
    Carry out feasibility tests and check whether the problem is solved.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    p_stat : TYPE
        DESCRIPTION.

    """
    # Check bounds
    p_stat['feasible'] = _presolve_infeasible_bounds(lp, tol)

    # Check constraints implied by equalities
    if not _presolve_complete(p_stat):
        p_stat['feasible'] = _presolve_infeasible_equality_constraints(lp, tol)

    # Check constraints implied by inequalities
    if not _presolve_complete(p_stat):
        p_stat['feasible'] = _presolve_infeasible_inequality_constraints(lp, tol)

    # The problem is solved if there are no remaining variables
    p_stat['solved'] = np.size(lp.c) == 0

    # Return status
    return p_stat


def _presolve_preprocess_problem(lp, revstack, p_stat, tol):
    """
    Preprocess the problem in order to simplify it.
    The function will set p_stat['loop'] to True if the problem has been
    changed (variables have been removed).
    The function detects if the problem has been solved and may also detect
    infeasibility or unboundedness. If so, the corresponding status is set
    and the function exits.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    revstack : list of functions
        Functions list containing the functions to reverse the presolve
        preprocessing.
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    lp : scipy.optimize._linprog_util._LPProblem
        The modified problem parameters
    revstack : list of functions
        Functions list containing the functions to reverse the presolve
        preprocessing.
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status

    """
    # Remove variables which are fixed by bounds
    lp, rev, p_stat = _presolve_remove_fixed_variables(lp, p_stat, tol)
    if rev is not None:
        revstack.append(rev)

    # Remove variables which are fixed by singletons in A_eq
    if not _presolve_complete(p_stat):
        lp, rev, p_stat = _presolve_remove_equation_row_singletons(lp, p_stat, tol)
        if rev is not None:
            revstack.append(rev)

    # Remove singleton inequalities in A_ub and merge them into bounds
    if not _presolve_complete(p_stat):
        lp, rev, p_stat = _presolve_remove_inequality_row_singletons(lp, p_stat, tol)
        if rev is not None:
            revstack.append(rev)

    # Remove empty rows
    if not _presolve_complete(p_stat):
        lp, rev, p_stat = _presolve_remove_empty_rows(lp, p_stat, tol)
        if rev is not None:
            revstack.append(rev)

    # Remove empty columns
    if not _presolve_complete(p_stat):
        lp, rev, p_stat = _presolve_remove_empty_columns(lp, p_stat)
        if rev is not None:
            revstack.append(rev)

    # Return status w.r.t. solved, feasibility
    return (lp, revstack, p_stat)


def _presolve_infeasible_bounds(lp, tol):
    """
    Check bounds.
    Return True if bounds are feasible.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    Boolean
        True if bounds are feasible

    """
    lb = lp.bounds[:, 0]
    ub = lp.bounds[:, 1]
    return np.all(ub >= lb - tol) and not np.any(lb == np.inf) and not np.any(ub == -np.inf)


def _presolve_infeasible_equality_constraints(lp, tol):
    """
    Check for equalities which cannot be fulfilled given the bounds.
    Given bounds for x lead to bounds for product A_eq*x. If b_eq
    does not satisfy these bounds, a solution is infeasible.
    Return True if bounds are feasible.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    Boolean
        True if bounds are feasible

    """
    # Handle empty A_eq separately
    # In this case, feasibility depends on if b_eq == 0 within tolerance
    if np.size(lp.A_eq) == 0:
        if np.size(lp.b_eq) > 0:
            return np.all(np.abs(lp.b_eq) < tol)
        else:
            return True

    # Create two matrices G and H the shape of A_eq. G with upper bound
    # values where elements of A_eq are positive, lower bound values
    # where negative, and 0 where 0. H the opposite.
    # Note: treat zero-values in A_eq separately, to avoid final
    # mutiplication of 0 and inf.
    if sps.issparse(lp.A_eq):
        # Get row indices, column indices and (non-zero) values from A_eq
        iA, jA, valA = sps.find(lp.A_eq)
        # Masks for positive and negative values
        pos = valA > 0
        neg = np.logical_not(pos)
        # Get row and column indices of positive and negative values in A_eq
        iApos = iA[pos]
        jApos = jA[pos]
        iAneg = iA[neg]
        jAneg = jA[neg]
        # Create new sparse matrices G and H
        # Using dense matrices may be faster but that would ignore the user's
        # considerations to use a sparse A_eq.
        G = sps.lil_matrix(lp.A_eq.shape)
        H = sps.lil_matrix(lp.A_eq.shape)
        # For these locations, set elements in G and H
        # G[...] = ... and H[...] = ... issue a SparseEfficiencyWarning:
        # Changing the sparsity structure of a csr_matrix is expensive.
        # For that reson the more straightforward calls below
        #  G[iApos, jApos] = lp.bounds[jApos, 1]
        #  G[iAneg, jAneg] = lp.bounds[jAneg, 0]
        #  H[iApos, jApos] = lp.bounds[jApos, 0]
        #  H[iAneg, jAneg] = lp.bounds[jAneg, 1]
        # have been reduced to two:
        iApn = np.concatenate((iApos, iAneg))
        jApn = np.concatenate((jApos, jAneg))
        G[iApn, jApn] = np.concatenate((lp.bounds[jApos, 1], lp.bounds[jAneg, 0]))
        H[iApn, jApn] = np.concatenate((lp.bounds[jApos, 0], lp.bounds[jAneg, 1]))
        # Row sums of element-wise product gives range between which equations
        # can vary.
        u_eq = np.sum(lp.A_eq.multiply(G), 1).flatten()
        l_eq = np.sum(lp.A_eq.multiply(H), 1).flatten()
    else:
        G = np.zeros(lp.A_eq.shape, dtype=float)  # lp.A_eq.copy()
        H = np.zeros(lp.A_eq.shape, dtype=float)  # lp.A_eq.copy()
        pos = lp.A_eq > 0
        neg = lp.A_eq < 0
        if pos.any():
            G = np.where(pos, lp.bounds[:, 1], G)
            H = np.where(pos, lp.bounds[:, 0], H)
        if neg.any():
            G = np.where(neg, lp.bounds[:, 0], G)
            H = np.where(neg, lp.bounds[:, 1], H)
        # Row sums of element-wise product gives range between which equations
        # can vary.
        u_eq = np.sum(lp.A_eq * G, 1)
        l_eq = np.sum(lp.A_eq * H, 1)

    # Check whether b_eq is within this range
    return np.all(l_eq - tol <= lp.b_eq) and np.all(lp.b_eq <= u_eq + tol)


def _presolve_infeasible_inequality_constraints(lp, tol):
    """
    Check for inequalities which cannot be fulfilled given the bounds.
    Given bounds for x lead to bounds L and U for product A_ub*x. If
    these bounds do not overlap <= b_ub, i.e. if L > b_ub, then a
    solution is infeasible.
    Return True if bounds are feasible.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    Boolean
        Return True if bounds are feasible.

    """
    # Handle empty A_ub separately
    # In this case, feasibility depends on if 0 <= b_ub
    if np.size(lp.A_ub) == 0:
        if np.size(lp.b_ub) > 0:
            return np.all(lp.b_ub >= 0)
        else:
            return True

    # Create a matrix H the shape of A_ub, with upper bound
    # values where elements of A_ub are negative, lower bound values
    # where positive, and 0 where 0.
    # Note: treat zero-values in A_ub separately, to avoid final
    # mutiplication of 0 and inf.
    if sps.issparse(lp.A_ub):
        # Get row indices, column indices and (non-zero) values from A_ub
        iA, jA, valA = sps.find(lp.A_ub)
        # Masks for positive and negative values
        pos = valA > 0
        neg = np.logical_not(pos)
        # Get row and column indices of positive and negative values in A_ub
        iApos = iA[pos]
        jApos = jA[pos]
        iAneg = iA[neg]
        jAneg = jA[neg]
        # Create new sparse matrix H
        # Using dense matrices may be faster but that would ignore the user's
        # considerations to use a sparse A_eq.
        H = sps.lil_matrix(lp.A_ub.shape)
        # For these locations, set elements in H
        # H[...] = ... issue a SparseEfficiencyWarning:
        # Changing the sparsity structure of a csr_matrix is expensive.
        # For that reason, the more straightforward two calls below
        #  H[iApos, jApos] = lp.bounds[jApos, 0]
        #  H[iAneg, jAneg] = lp.bounds[jAneg, 1]
        # have been reduced to one:
        iApn = np.concatenate((iApos, iAneg))
        jApn = np.concatenate((jApos, jAneg))
        H[iApn, jApn] = np.concatenate((lp.bounds[jApos, 0], lp.bounds[jAneg, 1]))
        # Row sums of element-wise product gives range between which equations
        # can vary.
        l_ub = np.sum(lp.A_ub.multiply(H), 1).flatten()
    else:
        H = np.zeros(lp.A_ub.shape, dtype=float)  # lp.A_ub.copy()
        H = np.where(lp.A_ub > 0, lp.bounds[:, 0], H)
        H = np.where(lp.A_ub < 0, lp.bounds[:, 1], H)
        # Row sums of element-wise product gives minimum of range between
        # which equations can vary.
        l_ub = np.sum(lp.A_ub * H, 1)

    # Check whether this minimum is always less than b_ub
    # l_ub <= b_ub < b_ub + tol
    return np.all(l_ub <= lp.b_ub + tol)


def _presolve_remove_fixed_variables(pp, p_stat, tol=1e-9):
    """
    Check for lower bounds which are equal to upper bounds.
    Remove these variable from the problem matrices.
    Return an update problem and a reversal function. If there is nothing to
    revert, this function is None.
    There is no check if the resulting problem is feasible.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    lp : scipy.optimize._linprog_util._LPProblem
        The modified problem parameters
    revstack : list of functions
        Functions list containing the functions to reverse the presolve
        preprocessing.
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status

    """
    # TODO Perhaps use a tolerance vector to allow for tolerance per variable
    # Check which indices are fixed, given the tolerance
    fixed_indices = np.abs(pp.bounds[:, 0] - pp.bounds[:, 1]) < tol

    # Number of fixed variables
    n_fixed = np.sum(fixed_indices)

    # No fixed variables
    if n_fixed == 0:
        return (pp, None, p_stat)

    # Remove fixed variables from equations
    keep_indices = np.logical_not(fixed_indices)
    # Fix variables indicated by fixed_indices (midway bounds)
    # Make sure A_eq and A_ub have equal number of columns even if they are
    # empty
    x_fixed = np.sum(pp.bounds[fixed_indices, :], 1) / 2
    b_eq = pp.b_eq - (pp.A_eq[:, fixed_indices]).dot(x_fixed)
    A_eq = pp.A_eq[:, keep_indices]
    b_ub = pp.b_ub - (pp.A_ub[:, fixed_indices]).dot(x_fixed)
    A_ub = pp.A_ub[:, keep_indices]
    pp = pp._replace(b_eq=b_eq, A_eq=A_eq, b_ub=b_ub, A_ub=A_ub,
                     c=pp.c[keep_indices], bounds=pp.bounds[keep_indices, :])
    # x0 can be None, so check before removing elements
    if pp.x0 is not None:
        # pp.x0 = pp.x0[keep_indices]
        pp = pp._replace(x0=pp.x0[keep_indices])

    def rev(x_mod):
        # Reverse operation only relevant for x
        # When removing elements at positions k1, k2, k3, ...
        # these must be replaced at (after) positions k1-1, k2-2, k3-3, ...
        # in the modified array. If x_mod is empty, just return x_fixed
        if np.size(x_mod) == 0:
            return x_fixed
        # A non-empty x_mod is a little more work
        i = np.flatnonzero(fixed_indices)
        # Number of variables to restore
        N = len(i)
        index_offset = np.arange(N)
        # Create insert indices
        insert_indices = i - index_offset
        x_rev = np.insert(x_mod.astype(float), insert_indices, x_fixed)
        return x_rev

    p_stat['loop'] = True
    if np.size(pp.c) == 0:
        p_stat['solved'] = True
        p_stat['loop'] = False

    return (pp, rev, p_stat)


def _presolve_remove_equation_row_singletons(pp, p_stat, tol=1e-9):
    """
    Remove row singletons from A_eq * x = b_eq.

    Return an updated problem and a reversal function. If there is nothing to
    revert, this function is None. Update the presolve status if conflicting
    values for the same variable have been detected.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    lp : scipy.optimize._linprog_util._LPProblem
        The modified problem parameters
    revstack : list of functions
        Functions list containing the functions to reverse the presolve
        preprocessing.
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status

    """
    if np.size(pp.A_eq) == 0:
        # No singleton rows
        return (pp, None, p_stat)

    # Sparse
    is_A_eq_sparse = sps.issparse(pp.A_eq)
    # Detect which rows of A_eq have a single non-zero element.
    # Check exact zeros, not approximate ones.
    # Create a mask vector (1-D array) to index these rows
    A_eq_nonzero = pp.A_eq != 0
    sing_row_indices = np.sum(A_eq_nonzero, axis=1) == 1
    if is_A_eq_sparse:
        # sing_row_indices is a numpy matrix, need a 1-D-array
        sing_row_indices = np.array(sing_row_indices).flatten()

    # Mask for non-zero elements in singleton rows of A_eq
    A_sr_nonzero = A_eq_nonzero[sing_row_indices, :]
    # Singleton variable counts (a variable may appear in multiple singleton rows)
    fixed_counts = np.sum(A_sr_nonzero, axis=0)
    # Mask for singleton variables
    fixed_indices = fixed_counts >= 1
    if is_A_eq_sparse:
        # sing_row_indices is a numpy matrix, need a 1-D-array
        fixed_indices = np.array(fixed_indices).flatten()

    # Number of fixed variables
    n_fixed = np.sum(fixed_indices)

    # No fixed variables
    if n_fixed == 0:
        #TODO mod loop?
        return (pp, None, p_stat)

    # Determine singleton variable values by elimination from each singleton.
    # row. If several rows eliminate the same variable, check whether the
    # resulting values are identical.
    # The b-values of singleton rows
    b_sr = pp.b_eq[sing_row_indices]
    # The non-zero A-element at these rows
    A_sr = pp.A_eq[sing_row_indices, :]
    A_sr = A_sr[A_sr_nonzero]
    # The value of the eliminated variable
    x_sr = b_sr / A_sr
    # Corresponding variable indices: column indices of the nonzero elements in
    # A_sr_nonzero. Use i to index rows, j to index columns & variables.
    (_, j_sr) = np.nonzero(A_sr_nonzero)

    # Create empty array to store fixed variables
    x_fixed = np.zeros((np.sum(fixed_indices), ))
    # Loop through indices of variables which will be eliminated (fixed)
    k = 0
    for j in np.flatnonzero(fixed_indices):
        # get value(s) from x_sr for variable j
        if is_A_eq_sparse:
            x_j = x_sr[0, j_sr == j]
            # Need a 1-D array
            x_j = np.array(x_j).flatten()
        else:
            x_j = x_sr[j_sr == j]
        # If more than one, check if all values identical
        if len(x_j) > 1 and any(np.abs(x_j - x_j[0]) > tol):
            # Solution infeasible
            p_stat['feasible'] = False
            return (pp, None, p_stat)
        x_fixed[k] = x_j[0]
        k = k + 1

    # Remove fixed variables and singleton rows from equations
    keep_indices = np.logical_not(fixed_indices)

    # Column removal and b_eq update
    # b_eq = b_eq - A_eq[:,j_fixed] * x_fixed
    b_eq = pp.b_eq - pp.A_eq[:, fixed_indices].dot(x_fixed)
    A_eq = pp.A_eq[:, keep_indices]
    # Row removal
    keep_row_indices = np.logical_not(sing_row_indices)
    A_eq = A_eq[keep_row_indices, :]
    b_eq = b_eq[keep_row_indices]
    pp = pp._replace(b_eq=b_eq, A_eq=A_eq)
    # Make sure A_eq and A_ub have equal number of columns even if they are
    # empty
    b_ub = pp.b_ub - pp.A_ub[:, fixed_indices].dot(x_fixed)
    A_ub = pp.A_ub[:, keep_indices]
    pp = pp._replace(b_ub=b_ub, A_ub=A_ub, c=pp.c[keep_indices],
                     bounds=pp.bounds[keep_indices, :])
    if pp.x0 is not None:
        pp = pp._replace(x0=pp.x0[keep_indices])

    def rev(x_mod):
        # Reverse operation only relevant for x
        # When removing elements at positions k1, k2, k3, ...
        # these must be replaced at (after) positions k1-1, k2-2, k3-3, ...
        # in the modified array
        if np.size(x_mod) == 0:
            return x_fixed
        # A non-empty x_mod is a little more work
        i = np.flatnonzero(fixed_indices)
        # Number of variables to restore
        N = len(i)
        index_offset = np.arange(N)
        # Create insert indices
        insert_indices = i - index_offset
        x_rev = np.insert(x_mod.astype(float), insert_indices, x_fixed)
        return x_rev

    p_stat['loop'] = True
    if np.size(pp.c) == 0:
        p_stat['solved'] = True
        p_stat['loop'] = False
    return (pp, rev, p_stat)


def _presolve_remove_inequality_row_singletons(pp, p_stat, tol=1e-9):
    """
    Remove row singletons from A_ub * x <= b_ub.

    A row singleton in the inequalities represents a simple bound.
    For example if variable x[j] in row i is a singleton, then
    A_ub[i,j] * x[j] <= b_ub[i]. This gives x[j] <= b_ub[i]/A_ub[i,j] if
    A_ub[i,j] > 0 and  x[j] >= b_ub[i]/A_ub[i,j] if A_ub[i,j] < 0.

    Merge this bound with already given bounds and remove the row in A_ub
    and b_ub.
    The result may be:
    * bounds narrowing: signal that the problem has changed;
    * a fixed variable: signal that the problem has changed; the next
      preprocessing loop will remove the fixed variable;
    * infeasible bounds: signal infeasibility.

    Return an updated problem; no variables will be removed so the reversal
    function is None. Return p_stat['feasible'] = False if conflicting bounds
    have been detected.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    lp : scipy.optimize._linprog_util._LPProblem
        The modified problem parameters
    revstack : list of functions
        Functions list containing the functions to reverse the presolve
        preprocessing.
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status

    """
    if np.size(pp.A_ub) == 0:
        # No singleton rows
        return (pp, None, p_stat)

    # Sparse
    is_A_ub_sparse = sps.issparse(pp.A_ub)
    # Detect which rows of A_ub have a single non-zero element.
    # Check exact zeros, not approximate ones.
    A_ub_nonzero = pp.A_ub != 0
    # Mask 1 x N_u, N_sr True values
    sing_row_indices = np.sum(A_ub_nonzero, axis=1) == 1
    if is_A_ub_sparse:
        # sing_row_indices is a numpy matrix, need an 1-D-array
        sing_row_indices = np.array(sing_row_indices).flatten()

    # Loop to go through all singleton variables to adjust bounds
    # TODO try to remove loop
    for i in np.nonzero(sing_row_indices)[0]:
        # Which variable ivolved?
        if is_A_ub_sparse:
            # A_ub_nonzero[i, :] is a 2D-array of size 1 x N
            # We need the col-indices from nonzero()
            j = np.nonzero(A_ub_nonzero[i, :])[1][0]
        else:
            # A_ub_nonzero[i, :] is a 1D-array with N elements
            # We need the first (and only) element in the tuple returned
            j = np.nonzero(A_ub_nonzero[i, :])[0]
        # Lower or upper bound
        val = pp.b_ub[i] / pp.A_ub[i, j]
        # Check with data in bounds
        bnds = pp.bounds
        lbj = bnds[j, 0]
        ubj = bnds[j, 1]
        if pp.A_ub[i, j] > 0:
            # Merge upper bound with [lbj, ubj]
            if val < lbj - tol:
                # New upper bound infeasible with current bounds
                p_stat['feasible'] = False
                return (pp, None, p_stat)
            elif np.abs(val - lbj) < tol:
                # New fixed variable, new upper bound equals lower bound
                bnds[j, 1] = val
            elif val < ubj:
                # Use new upper bound instead of current
                bnds[j, 1] = val
        else:
            # Merge lower bound with [lbj, ubj]
            if val > ubj + tol:
                # New lower bound infeasible with current bounds
                p_stat['feasible'] = False
                return (pp, None, p_stat)
            elif np.abs(val - ubj) < tol:
                # New fixed variable, new lower bound equals upper bound
                bnds[j, 0] = val
            elif val > lbj:
                # Use new lower bound instead of current
                bnds[j, 0] = val

    # If there are changes, then save bounds, A_ub and b_ub
    if np.sum(sing_row_indices) > 0:
        # Remove the singleton rows from A_ub and b_ub
        keep_row_indices = np.logical_not(sing_row_indices)
        A_ub = pp.A_ub[keep_row_indices, :]
        b_ub = pp.b_ub[keep_row_indices]
        # Also include changed bounds, even if nothing changed
        pp = pp._replace(bounds=bnds, A_ub=A_ub, b_ub=b_ub)
        # Mark problem as changed
        p_stat['loop'] = True
        # No need to update 'solved', no variables removed

    return (pp, None, p_stat)


def _presolve_remove_empty_rows(pp, p_stat, tol=1e-9):
    """
    Check for empty rows in A_eq and A_ub.
    Remove these rows and the corresponding elements in b_eq and b_ub.

    Parameters
    ----------
    lp : scipy.optimize._linprog_util._LPProblem
        The problem parameters
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status
    tol : float
        Tolerance to be used for checking equality

    Returns
    -------
    lp : scipy.optimize._linprog_util._LPProblem
        The modified problem parameters
    revstack : list of functions
        Functions list containing the functions to reverse the presolve
        preprocessing.
    p_stat : dictionary with keys 'solved', 'feasible', 'bounded', 'loop'
        Presolve problem status

    """
    # Equations
    if np.size(pp.A_eq) == 0:
        zero_row_indices = []
    else:
        zero_row_indices = np.array(np.sum(pp.A_eq != 0, axis=1) == 0).flatten()
    if np.any(zero_row_indices):
        if np.any(np.abs(pp.b_eq[zero_row_indices]) > tol):
            # A_eq * x = b_eq: empty row makes feasible problem only if
            # corresponding element in b_eq is zero
            p_stat['feasible'] = False
            return (pp, None, p_stat)
        else:
            # if RHS is zero, we can eliminate this equation entirely.
            nonzero_row_indices = np.logical_not(zero_row_indices)
            A_eq = pp.A_eq[nonzero_row_indices, :]
            b_eq = pp.b_eq[nonzero_row_indices]
            pp = pp._replace(b_eq=b_eq, A_eq=A_eq)
            p_stat['loop'] = True

    # Inequalities
    if np.size(pp.A_ub) == 0:
        zero_row_indices = []
    else:
        zero_row_indices = np.array(np.sum(pp.A_ub != 0, axis=1) == 0).flatten()
    if np.any(zero_row_indices):
        if np.any(pp.b_ub[zero_row_indices] < -tol):
            # A_ub * x <= b_ub: empty row makes feasible problem only if
            # corresponding element in b_ub is bigger than zero.
            p_stat['feasible'] = False
            return (pp, None, p_stat)
        else:
            # if RHS is zero or below, we can eliminate this inequality.
            nonzero_row_indices = np.logical_not(zero_row_indices)
            A_ub = pp.A_ub[nonzero_row_indices, :]
            b_ub = pp.b_ub[nonzero_row_indices]
            pp = pp._replace(b_ub=b_ub, A_ub=A_ub)
            p_stat['loop'] = True

    return (pp, None, p_stat)


def _presolve_remove_empty_columns(pp, p_stat):
    """
    Check for empty columns in A_eq and A_ub.
    If a variable does not appear in any of the equations and
    inequalities, an optimal choice can be made for it using its bounds
    and the corrsponding value in c.
    This also holds if both A_eq and A_ub are empty.
    Remove these empty columns and the corresponding variables.
    """
    # Determine indices of zero columns (= indices of fixed variables)
    fixed_indices = np.full((len(pp.c), ), True)
    if np.size(pp.A_ub) > 0:
        zero_coll_indices_ub = np.array(np.sum(pp.A_ub != 0, axis=0) == 0)
        fixed_indices = np.logical_and(fixed_indices, zero_coll_indices_ub)
    if np.size(pp.A_eq) > 0:
        zero_coll_indices_eq = np.array(np.sum(pp.A_eq != 0, axis=0) == 0)
        fixed_indices = np.logical_and(fixed_indices, zero_coll_indices_eq)
    fixed_indices = fixed_indices.flatten()

    # If no empty columns, return
    if not np.any(fixed_indices):
        return (pp, None, p_stat)

    # Determine optimal values of fixed variables
    x_fixed = np.zeros((len(pp.c), ))
    fixed_indices_neg_c = np.logical_and(fixed_indices, pp.c < 0).flatten()
    x_fixed[fixed_indices_neg_c] = pp.bounds[fixed_indices_neg_c, 1]
    fixed_indices_pos_c = np.logical_and(fixed_indices, pp.c > 0).flatten()
    x_fixed[fixed_indices_pos_c] = pp.bounds[fixed_indices_pos_c, 0]

    # Also check c == 0...
    fixed_indices_zero_c = np.logical_and(fixed_indices, pp.c == 0)
    # For these indices select a value within the currrent bounds
    # (actually arbitrary: any value within the bounds is equally valid)
    # For twice unbounded choose 0, for once unbounded choose the other bound,
    # for bounded choose mid-way.
    if np.sum(fixed_indices_zero_c) > 0:
        lb_zero_c = pp.bounds[fixed_indices_zero_c, 0]
        ub_zero_c = pp.bounds[fixed_indices_zero_c, 1]
        lb_zero_c_inf = lb_zero_c == -np.inf
        ub_zero_c_inf = ub_zero_c == np.inf
        bothb_zero_c_inf = np.logical_and(lb_zero_c_inf, ub_zero_c_inf)
        lb_zero_c[bothb_zero_c_inf] = 0
        ub_zero_c[bothb_zero_c_inf] = 0
        lb_zero_c[lb_zero_c_inf] = ub_zero_c[lb_zero_c_inf]
        ub_zero_c[ub_zero_c_inf] = lb_zero_c[ub_zero_c_inf]
        x_fixed[fixed_indices_zero_c] = (ub_zero_c + lb_zero_c) / 2

    x_fixed = x_fixed[fixed_indices]
    if np.any(np.isinf(x_fixed)):
        p_stat['bounded'] = False
        return (pp, None, p_stat)

    # Modify c, A_ub, A_eq, bounds, x0
    keep_indices = np.logical_not(fixed_indices)
    # Column removal
    A_eq = pp.A_eq[:, keep_indices]
    A_ub = pp.A_ub[:, keep_indices]
    pp = pp._replace(A_eq=A_eq, A_ub=A_ub)
    # pp.c = pp.c[keep_indices]
    # pp.bounds = pp.bounds[keep_indices, :]
    pp = pp._replace(c=pp.c[keep_indices], bounds=pp.bounds[keep_indices, :])
    if pp.x0 is not None:
        # pp.x0 = pp.x0[keep_indices]
        pp = pp._replace(x0=pp.x0[keep_indices])

    def rev(x_mod):
        # Reverse operation only relevant for x
        # When removing elements at positions k1, k2, k3, ...
        # these must be replaced at (after) positions k1-1, k2-2, k3-3, ...
        # in the modified array
        if np.size(x_mod) == 0:
            return x_fixed
        # A non-empty x_mod is a little more work
        i = np.flatnonzero(fixed_indices)
        # Number of variables to restore
        N = len(i)
        index_offset = np.arange(N)
        # Create insert indices
        insert_indices = i - index_offset
        x_rev = np.insert(x_mod.astype(float), insert_indices, x_fixed)
        return x_rev

    p_stat['loop'] = True
    if np.size(pp.c) == 0:
        p_stat['solved'] = True
        p_stat['loop'] = False
    return (pp, rev, p_stat)


def _parse_linprog(lp, options):
    """
    Parse the provided linear programming problem

    ``_parse_linprog`` employs two main steps ``_check_sparse_inputs`` and
    ``_clean_inputs``. ``_check_sparse_inputs`` checks for sparsity in the
    provided constraints (``A_ub`` and ``A_eq) and if these match the provided
    sparsity optional values.

    ``_clean inputs`` checks of the provided inputs. If no violations are
    identified the objective vector, upper bound constraints, equality
    constraints, and simple bounds are returned in the expected format.

    Parameters
    ----------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : various valid formats, optional
            The bounds of ``x``, as ``min`` and ``max`` pairs.
            If bounds are specified for all N variables separately, valid formats are:
            * a 2D array (2 x N or N x 2);
            * a sequence of N sequences, each with 2 values.
            If all variables have the same bounds, a single pair of values can
            be specified. Valid formats are:
            * a sequence with 2 scalar values;
            * a sequence with a single element containing 2 scalar values.
            If all variables have a lower bound of 0 and no upper bound, the bounds
            parameter can be omitted (or given as None).
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    options : dict
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.

    Returns
    -------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, as ``min`` and ``max`` pairs, one for each of the N
            elements of ``x``. The N x 2 array contains lower bounds in the first
            column and upper bounds in the 2nd. Unbounded variables have lower
            bound -np.inf and/or upper bound np.inf.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    options : dict, optional
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.

    """
    if options is None:
        options = {}

    solver_options = {k: v for k, v in options.items()}
    solver_options, A_ub, A_eq = _check_sparse_inputs(solver_options, lp.A_ub, lp.A_eq)
    # Convert lists to numpy arrays, etc...
    lp = _clean_inputs(lp._replace(A_ub=A_ub, A_eq=A_eq))
    return lp, solver_options


def _get_Abc(lp, c0):
    """
    Given a linear programming problem of the form:

    Minimize::

        c @ x

    Subject to::

        A_ub @ x <= b_ub
        A_eq @ x == b_eq
         lb <= x <= ub

    where ``lb = 0`` and ``ub = None`` unless set in ``bounds``.

    Return the problem in standard form:

    Minimize::

        c @ x

    Subject to::

        A @ x == b
            x >= 0

    by adding slack variables and making variable substitutions as necessary.

    Parameters
    ----------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, lower bounds in the 1st column, upper
            bounds in the 2nd column. The bounds are possibly tightened
            by the presolve procedure.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables.

    Returns
    -------
    A : 2-D array
        2-D array such that ``A`` @ ``x``, gives the values of the equality
        constraints at ``x``.
    b : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in A (for standard form problem).
    c : 1-D array
        Coefficients of the linear objective function to be minimized (for
        standard form problem).
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables.
    x0 : 1-D array
        Starting values of the independent variables, which will be refined by
        the optimization algorithm

    References
    ----------
    .. [9] Bertsimas, Dimitris, and J. Tsitsiklis. "Introduction to linear
           programming." Athena Scientific 1 (1997): 997.

    """
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = lp

    if sps.issparse(A_eq):
        sparse = True
        A_eq = sps.csr_matrix(A_eq)
        A_ub = sps.csr_matrix(A_ub)

        def hstack(blocks):
            return sps.hstack(blocks, format="csr")

        def vstack(blocks):
            return sps.vstack(blocks, format="csr")

        zeros = sps.csr_matrix
        eye = sps.eye
    else:
        sparse = False
        hstack = np.hstack
        vstack = np.vstack
        zeros = np.zeros
        eye = np.eye

    # Variables lbs and ubs (see below) may be changed, which feeds back into
    # bounds, so copy.
    bounds = np.array(bounds, copy=True)

    # modify problem such that all variables have only non-negativity bounds
    lbs = bounds[:, 0]
    ubs = bounds[:, 1]
    m_ub, n_ub = A_ub.shape

    lb_none = np.equal(lbs, -np.inf)
    ub_none = np.equal(ubs, np.inf)
    lb_some = np.logical_not(lb_none)
    ub_some = np.logical_not(ub_none)

    # unbounded below: substitute xi = -xi' (unbounded above)
    # if -inf <= xi <= ub, then -ub <= -xi <= inf, so swap and invert bounds
    l_nolb_someub = np.logical_and(lb_none, ub_some)
    i_nolb = np.nonzero(l_nolb_someub)[0]
    lbs[l_nolb_someub], ubs[l_nolb_someub] = (
        -ubs[l_nolb_someub], -lbs[l_nolb_someub])
    lb_none = np.equal(lbs, -np.inf)
    ub_none = np.equal(ubs, np.inf)
    lb_some = np.logical_not(lb_none)
    ub_some = np.logical_not(ub_none)
    c[i_nolb] *= -1
    if x0 is not None:
        x0[i_nolb] *= -1
    if len(i_nolb) > 0:
        if A_ub.shape[0] > 0:  # sometimes needed for sparse arrays... weird
            A_ub[:, i_nolb] *= -1
        if A_eq.shape[0] > 0:
            A_eq[:, i_nolb] *= -1

    # upper bound: add inequality constraint
    i_newub, = ub_some.nonzero()
    ub_newub = ubs[ub_some]
    n_bounds = len(i_newub)
    if n_bounds > 0:
        shape = (n_bounds, A_ub.shape[1])
        if sparse:
            idxs = (np.arange(n_bounds), i_newub)
            A_ub = vstack((A_ub, sps.csr_matrix((np.ones(n_bounds), idxs),
                                                shape=shape)))
        else:
            A_ub = vstack((A_ub, np.zeros(shape)))
            A_ub[np.arange(m_ub, A_ub.shape[0]), i_newub] = 1
        b_ub = np.concatenate((b_ub, np.zeros(n_bounds)))
        b_ub[m_ub:] = ub_newub

    A1 = vstack((A_ub, A_eq))
    b = np.concatenate((b_ub, b_eq))
    c = np.concatenate((c, np.zeros((A_ub.shape[0],))))
    if x0 is not None:
        x0 = np.concatenate((x0, np.zeros((A_ub.shape[0],))))
    # unbounded: substitute xi = xi+ + xi-
    l_free = np.logical_and(lb_none, ub_none)
    i_free = np.nonzero(l_free)[0]
    n_free = len(i_free)
    c = np.concatenate((c, np.zeros(n_free)))
    if x0 is not None:
        x0 = np.concatenate((x0, np.zeros(n_free)))
    A1 = hstack((A1[:, :n_ub], -A1[:, i_free]))
    c[n_ub: n_ub + n_free] = -c[i_free]
    if x0 is not None:
        i_free_neg = x0[i_free] < 0
        x0[np.arange(n_ub, A1.shape[1])[i_free_neg]] = -x0[i_free[i_free_neg]]
        x0[i_free[i_free_neg]] = 0

    # add slack variables
    A2 = vstack([eye(A_ub.shape[0]), zeros((A_eq.shape[0], A_ub.shape[0]))])

    A = hstack([A1, A2])

    # lower bound: substitute xi = xi' + lb
    # now there is a constant term in objective
    i_shift = np.nonzero(lb_some)[0]
    lb_shift = lbs[lb_some].astype(float)
    c0 += np.sum(lb_shift * c[i_shift])
    if sparse:
        b = b.reshape(-1, 1)
        A = A.tocsc()
        b -= (A[:, i_shift] * sps.diags(lb_shift)).sum(axis=1)
        b = b.ravel()
    else:
        b -= (A[:, i_shift] * lb_shift).sum(axis=1)
    if x0 is not None:
        x0[i_shift] -= lb_shift

    return A, b, c, c0, x0


def _round_to_power_of_two(x):
    """
    Round elements of the array to the nearest power of two.
    """
    return 2**np.around(np.log2(x))


def _autoscale(A, b, c, x0):
    """
    Scales the problem according to equilibration from [12].
    Also normalizes the right hand side vector by its maximum element.
    """
    m, n = A.shape

    C = 1
    R = 1

    if A.size > 0:

        R = np.max(np.abs(A), axis=1)
        if sps.issparse(A):
            R = R.toarray().flatten()
        R[R == 0] = 1
        R = 1 / _round_to_power_of_two(R)
        A = sps.diags(R) * A if sps.issparse(A) else A * R.reshape(m, 1)
        b = b * R

        C = np.max(np.abs(A), axis=0)
        if sps.issparse(A):
            C = C.toarray().flatten()
        C[C == 0] = 1
        C = 1 / _round_to_power_of_two(C)
        A = A * sps.diags(C) if sps.issparse(A) else A * C
        c = c * C

    b_scale = np.max(np.abs(b)) if b.size > 0 else 1
    if b_scale == 0:
        b_scale = 1.
    b = b / b_scale

    if x0 is not None:
        x0 = x0 / b_scale * (1 / C)
    return A, b, c, x0, C, b_scale


def _unscale(x, C, b_scale):
    """
    Convert solution to _autoscale problem -> solution to original problem.
    """
    try:
        n = len(C)
        # fails if sparse or scalar; that's OK.
        # this is only needed for original simplex (never sparse)
    except TypeError:
        n = len(x)

    return x[:n] * b_scale * C


def _display_summary(message, status, fun, iteration):
    """
    Print the termination summary of the linear program

    Parameters
    ----------
    message : str
            A string descriptor of the exit status of the optimization.
    status : int
        An integer representing the exit status of the optimization::

                0 : Optimization terminated successfully
                1 : Iteration limit reached
                2 : Problem appears to be infeasible
                3 : Problem appears to be unbounded
                4 : Serious numerical difficulties encountered

    fun : float
        Value of the objective function.
    iteration : iteration
        The number of iterations performed.
    """
    print(message)
    if status in (0, 1):
        print("         Current function value: {0: <12.6f}".format(fun))
    print("         Iterations: {0:d}".format(iteration))


def _postsolve(x, postsolve_args, complete=False):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1-D array
        Solution vector to the standard-form problem.
    postsolve_args : tuple
        Data needed by _postsolve to convert the solution to the standard-form
        problem into the solution to the original problem, including:

    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, lower bounds in the 1st column, upper
            bounds in the 2nd column. The bounds are possibly tightened
            by the presolve procedure.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    revstack: list of functions
        the functions in the list reverse the operations of _presolve()
        the function signature is x_org = f(x_mod), where x_mod is the result
        of a presolve step and x_org the value at the start of the step
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)

    Returns
    -------
    x : 1-D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack : 1-D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub @ x``
    con : 1-D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq @ x``
    """
    # note that all the inputs are the ORIGINAL, unmodified versions
    # no rows, columns have been removed

    (c, A_ub, b_ub, A_eq, b_eq, bounds, x0), revstack, C, b_scale = postsolve_args

    x = _unscale(x, C, b_scale)

    # Undo variable substitutions of _get_Abc()
    # if "complete", problem was solved in presolve; don't do anything here
    n_x = bounds.shape[0]
    if not complete and bounds is not None:  # bounds are never none, probably
        n_unbounded = 0
        for i, bi in enumerate(bounds):
            lbi = bi[0]
            ubi = bi[1]
            if lbi == -np.inf and ubi == np.inf:
                n_unbounded += 1
                x[i] = x[i] - x[n_x + n_unbounded - 1]
            else:
                if lbi == -np.inf:
                    x[i] = ubi - x[i]
                else:
                    x[i] += lbi
    # all the rest of the variables were artificial
    x = x[:n_x]

    # If there were variables removed from the problem, add them back into the
    # solution vector
    # Apply the functions in revstack (reverse direction)
    for rev in reversed(revstack):
        x = rev(x)

    fun = x.dot(c)
    slack = b_ub - A_ub.dot(x)  # report slack for ORIGINAL UB constraints
    # report residuals of ORIGINAL EQ constraints
    con = b_eq - A_eq.dot(x)

    return x, fun, slack, con


def _check_result(x, fun, status, slack, con, bounds, tol, message):
    """
    Check the validity of the provided solution.

    A valid (optimal) solution satisfies all bounds, all slack variables are
    negative and all equality constraint residuals are strictly non-zero.
    Further, the lower-bounds, upper-bounds, slack and residuals contain
    no nan values.

    Parameters
    ----------
    x : 1-D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties encountered

    slack : 1-D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub @ x``
    con : 1-D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq @ x``
    bounds : 2D array
        The bounds on the original variables ``x``
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
             4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.
    """
    # Somewhat arbitrary, but status 5 is very unusual
    tol = np.sqrt(tol) * 10

    contains_nans = (
        np.isnan(x).any() or np.isnan(fun) or np.isnan(slack).any() or np.isnan(con).any()
    )

    if contains_nans:
        is_feasible = False
    else:
        invalid_bounds = (x < bounds[:, 0] - tol).any() or (x > bounds[:, 1] + tol).any()
        invalid_slack = status != 3 and (slack < -tol).any()
        invalid_con = status != 3 and (np.abs(con) > tol).any()
        is_feasible = not (invalid_bounds or invalid_slack or invalid_con)

    if status == 0 and not is_feasible:
        status = 4
        message = ("The solution does not satisfy the constraints within the "
                   "required tolerance of " + "{:.2E}".format(tol) + ", yet "
                   "no errors were raised and there is no certificate of "
                   "infeasibility or unboundedness. This is known to occur "
                   "if the `presolve` option is False and the problem is "
                   "infeasible. This can also occur due to the limited "
                   "accuracy of the `interior-point` method. Check whether "
                   "the slack and constraint residuals are acceptable; "
                   "if not, consider enabling presolve, reducing option "
                   "`tol`, and/or using method `revised simplex`. "
                   "If you encounter this message under different "
                   "circumstances, please submit a bug report.")
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
        # nearly basic feasible solution. Postsolving can make the solution
        # basic, however, this solution is NOT optimal
        raise ValueError(message)

    return status, message
