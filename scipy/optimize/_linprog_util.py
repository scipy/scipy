"""
Method agnostic utility functions for linear progamming
"""

import numpy as np
import scipy.sparse as sps
from warnings import warn
from .optimize import OptimizeWarning
from scipy.optimize._remove_redundancy import (
    _remove_redundancy, _remove_redundancy_sparse, _remove_redundancy_dense
    )
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
        If a specific variable has no lower bound and/or no upper bound,
        specify the respective bound as -numpy.inf, numpy.inf or None.
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
        c = np.array(c, dtype=np.float, copy=True).squeeze()
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: c must be a 1-D array of numerical "
            "coefficients")
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
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: A_ub must be a 2-D array "
            "of numerical values")
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
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: b_ub must be a 1-D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_ub")
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
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: A_eq must be a 2-D array "
            "of numerical values")
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
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: b_eq must be a 1-D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_eq")
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
        except ValueError:
            raise TypeError(
                "Invalid input for linprog: x0 must be a 1-D array of "
                "numerical coefficients")
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
    # Creating at least a 2-D array simplifies the cases to distinguish.
    try:
        bounds_conv = np.atleast_2d(np.array(bounds, dtype=float))
    except ValueError as e:
        raise ValueError(
            "Invalid input for linprog: unable to interpret bounds, "
            "check values and dimensions: " + e.args[0])
    except TypeError as e:
        raise TypeError(
            "Invalid input for linprog: unable to interpret bounds, "
            "check values and dimensions: " + e.args[0])

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
    elif bsh[0] == 1 and bsh[1] == 0:
        # [] converts to a (1,0) array (with no elements)
        bounds_clean[:, 1] = np.inf
    elif bsh[0] == 1 and bsh[1] == 1 and np.isnan(bounds_conv[0, 0]):
        # None converts to a (1,1) array with a nan-value
        bounds_clean[:, 1] = np.inf
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

lp, revstack, complete, status, message
    lp : see above
        The problem is modified by a number of presolve steps.
    revstack : list of functions
        List of functions, the application of which recreates the original
        ``c``, ``x`` and ``x0``.
        Each function corresponds to a presolve step taken. It takes ``c``,
        ``x`` and ``x0`` from after a step and returns the values from before
        the presolve step.
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
    # There are many improvements that can be made, namely:
    #  * implement remaining checks from [5]
    #  * loop presolve until no additional changes are made
    #  * implement additional efficiency improvements in redundancy removal [2]

    # Default return values
    complete = False
    message = None
    status = 0

    # List of functions for undoing presolve modifications
    revstack = []

    # Check infeasible bounds
    (lp, _, status) = _presolve_infeasible_bounds(lp)
    if status == 2:
        message = "The problem is (trivially) infeasible. One or more lower "
        "bounds are higher than corresponding upper bounds."
        return (lp, revstack, True, 2, message)

    # Check infeasible equality constraints
    (lp, _, status) = _presolve_infeasible_equality_constraints(lp)
    if status == 2:
        message = "The problem is (trivially) infeasible. One or more "
        "equalities cannot be fulfilled given the bounds."
        return (lp, revstack, True, 2, message)

    # Check infeasible inequality constraints
    (lp, _, status) = _presolve_infeasible_inequality_constraints(lp)
    if status == 2:
        message = "The problem is (trivially) infeasible. One or more "
        "inequalities cannot be fulfilled given the bounds."
        return (lp, revstack, True, 2, message)

    # Remove fixed variables
    (lp, rev, status) = _presolve_remove_fixed_variables(lp)
    if status == 0:
        message = "The solution was determined in presolve as there are "
        "no non-trivial constraints."
        return (lp, revstack, True, 0, message)
    elif status == 2:
        message = "The problem is (trivially) infeasible. One or more "
        "equalities cannot be fulfilled given the bounds."
        return (lp, revstack, True, 2, message)
    elif status == 5:
        revstack.append(rev)

    loop_modifications = True
    while loop_modifications:
        loop_modifications = False
        # Remove row singletons
        (lp, rev, status) = _presolve_remove_row_singletons(lp)
        if status == 0:
            message = "The solution was determined in presolve as there are "
            "no non-trivial constraints."
            return (lp, revstack, True, 0, message)
        elif status == 2:
            message = "The problem is (trivially) infeasible. One or more "
            "equalities cannot be fulfilled given the bounds."
            return (lp, revstack, True, 2, message)
        elif status == 5:
            revstack.append(rev)
            loop_modifications = True

    # Remove empty rows
    (lp, rev, status) = _presolve_remove_empty_rows(lp)
    if status == 2:
        message = "The problem is (trivially) infeasible due to a row of "
        "zeros in the (in)equality constraint matrix with a nonzero "
        "corresponding constraint value."
        return (lp, revstack, True, 2, message)

    # Remove empty columns
    (lp, rev, status) = _presolve_remove_empty_columns(lp)
    if status == 3:
        message = "The solution was determined in presolve as there are "
        "no non-trivial constraints."
        return (lp, revstack, True, 3, message)
    elif status == 5:
        revstack.append(rev)

    # Check if problem solved
    c = lp[0]
    if len(c) == 0:
        message = "The solution was determined in presolve as there are "
        "no non-trivial constraints."
        return (lp, revstack, True, 0, message)

    # At this point, complete = False, status cannot be anything else than 0
    # and message = None.
    return (lp, revstack, False, 0, None)


def _presolve_infeasible_bounds(pp):
    """
    Check for lower bounds which are larger than upper bounds
    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    """
    bounds = pp[5]
    infeasible_bounds = bounds[:, 0] > bounds[:, 1]
    if np.any(infeasible_bounds):
        status = 2  # infeasible bounds
    else:
        status = 6  # no change

    return (pp, None, status)


def _presolve_infeasible_equality_constraints(pp):
    """
    Check for equalities which cannot be fulfilled given the bounds.
    Given bounds for x lead to bounds for product A_eq*x. If b_eq
    does not satisfy these bounds, a solution is infeasible.
    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    """
    A_eq = pp[3]
    b_eq = pp[4]
    bounds = pp[5]

    if A_eq is None:
        status = 6  # no change
        return (pp, None, status)

    # Create two matrices G and H the shape of A_eq. G with upper bound
    # values where elements of A_eq are positive, lower bound values
    # where negative, and 0 where 0. H the opposite.
    # Note: treat zero-values in A_eq separately, to avoid final
    # mutiplication of 0 and inf.
    pos = A_eq > 0
    neg = A_eq < 0
    n_eq = len(b_eq)
    # Repeat lower bounds to get a matrix the shape of A_eq
    Ml = np.tile(bounds[:, 0], [n_eq, 1])
    # Repeat upper bounds to get a matrix the shape of A_eq
    Mu = np.tile(bounds[:, 1], [n_eq, 1])
    # Create G and H
    G = np.zeros(A_eq.shape)
    G[neg] = Ml[neg]
    G[pos] = Mu[pos]
    H = np.zeros(A_eq.shape)
    H[pos] = Ml[pos]
    H[neg] = Mu[neg]
    # Row sums of element-wise product gives range between which equations
    # can vary.
    u_eq = np.sum(A_eq * G, 1)
    l_eq = np.sum(A_eq * H, 1)
    # Check whether b_eq is within this range
    if np.any(b_eq < l_eq) or np.any(b_eq > u_eq):
        # Infeasible constraints
        status = 2  # infeasible bounds
    else:
        # No infeasible constraints detected
        status = 6  # no change
    return (pp, None, status)


def _presolve_infeasible_inequality_constraints(pp):
    """
    Check for inequalities which cannot be fulfilled given the bounds.
    Given bounds for x lead to bounds L and U for product A_ub*x. If
    these bounds do not overlap <= b_ub, i.e. if L > b_ub, then a
    solution is infeasible.
    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    """
    A_ub = pp[1]
    b_ub = pp[2]
    bounds = pp[5]

    if A_ub is None:
        status = 6  # no change
        return (pp, None, status)

    # Create a matrix H the shape of A_eq, with upper bound
    # values where elements of A_eq are negative, lower bound values
    # where positive, and 0 where 0.
    # Note: treat zero-values in A_eq separately, to avoid final
    # mutiplication of 0 and inf.
    pos = A_ub > 0
    neg = A_ub < 0
    n_eq = len(b_ub)
    # Repeat lower bounds to get a matrix the shape of A_eq
    Ml = np.tile(bounds[:, 0], [n_eq, 1])
    # Repeat upper bounds to get a matrix the shape of A_eq
    Mu = np.tile(bounds[:, 1], [n_eq, 1])
    # Create H
    H = np.zeros(A_ub.shape)
    H[pos] = Ml[pos]
    H[neg] = Mu[neg]
    # Row sums of element-wise product gives range between which equations
    # can vary.
    l_ub = np.sum(A_ub * H, 1)
    # Check whether b_eq is within this range
    if np.any(b_ub < l_ub):
        # Infeasible constraints
        status = 2  # infeasible bounds
    else:
        # No infeasible constraints detected
        status = 6  # no change
    return (pp, None, status)


def _presolve_remove_fixed_variables(pp, tol=1e-9):
    """
    Check for lower bounds which are equal to upper bounds.
    Remove these variable from the problem matrices.
    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    """
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = pp
    # TODO Use a tolerance vector to allow for tolerance per variable
    # Check which indices are fixed, given the tolerance
    fixed_indices = np.abs(bounds[:, 0] - bounds[:, 1]) <= tol

    # Number of fixed variables
    n_fixed = np.sum(fixed_indices)

    # No fixed variables
    if n_fixed == 0:
        status = 6  # no change
        return (pp, None, status)

    # All variables fixed
    if n_fixed == len(c):
        # If all variables have been fixed, the problem is either solved or
        # infeasible.
        # The fixed variables results in 0 >= b_ub and 0 == b_eq.
        # So check if the transformed b_eq equals 0 (within tolerance range).
        # For b_ub there is no need to use a tolerance, just check if b_ub <= 0.
        # The tolerance to use for element i in b_eq is:
        # tol_b[i]^2 = A_eq[i,1]^2 * tol_x[1]^2 + ... + A_eq[i,N] * tol_x[N]^2
        # if there are N variables in x.
        # All x-tolerances are the same: tol_x[j] = tol.
        x_fixed = np.sum(bounds[fixed_indices, :], 1) / 2
        if A_eq is not None:
            tol_b_eq = np.sqrt(np.sum(np.square(A_eq), 1)) * tol
            # Calculate residu A_eq * x - b_eq and check if zero
            r = np.dot(A_eq, x_fixed) - b_eq
            sol_eq = np.all(np.abs(r) <= tol_b_eq)
        else:
            sol_eq = True
        if A_ub is not None:
            # Calculate residu A_ub * x - b_ub and check if >= zero
            r = np.dot(A_ub, x_fixed) - b_ub
            sol_ub = np.all(r >= 0)
        else:
            sol_ub = True
        if sol_eq and sol_ub:
            # Solution determined in presolve
            status = 0
        else:
            # Solution infeasible
            status = 2
        return (pp, None, status)

    # Remove fixed variables from equations
    keep_indices = np.logical_not(fixed_indices)
    # Fix variables indicated by fixed_indices (midway bounds)
    x_fixed = np.sum(bounds[fixed_indices, :], 1) / 2
    if b_eq is not None:
        b_eq = b_eq - np.dot(A_eq[:, fixed_indices], x_fixed)
        A_eq = A_eq[:, keep_indices]
    if b_ub is not None:
        b_ub = b_ub - np.dot(A_ub[:, fixed_indices], x_fixed)
        A_ub = A_ub[:, keep_indices]
    c_fixed = c[fixed_indices]
    c = c[keep_indices]
    bounds = bounds[keep_indices, :]
    if x0 is not None:
        x0_fixed = x0[fixed_indices]
        x0 = x0[keep_indices]

    def rev(c_mod, x_mod, x0_mod):
        # Reverse operation only relevant for c, x and x0
        # When removing elements at positions k1, k2, k3, ...
        # these must be replaced at (after) positions k1-1, k2-2, k3-3, ...
        # in the modified array
        i = np.flatnonzero(fixed_indices)
        # Number of variables to restore
        N = len(i)
        index_offset = list(range(N))
        # Create insert indices
        insert_indices = np.subtract(i, index_offset).flatten()
        c_rev = np.insert(c_mod.astype(float), insert_indices, c_fixed)
        x_rev = np.insert(x_mod.astype(float), insert_indices, x_fixed)
        if x0 is not None:
            x0_rev = np.insert(x0_mod.astype(float), insert_indices, x0_fixed)
        else:
            x0_rev = None
        return (c_rev, x_rev, x0_rev)

    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)  # TODO check if necessary
    status = 5  # modifications
    return (pp, rev, status)


def _presolve_remove_row_singletons(pp, tol=1e-9):
    """
    Remove row singletons.
    """
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = pp

    if A_eq is None:
        # No singleton rows
        status = 6  # no change
        return (pp, None, status)

    # Detect which rows of A_eq have a single non-zero element.
    # Check exact zeros, not approximate ones.
    A_eq_nonzero = A_eq != 0
    sing_row_indices = np.sum(A_eq_nonzero, axis=1) == 1

    # Mask for non-zero elements in singleton rows of A_eq
    A_sr_nonzero = A_eq_nonzero[sing_row_indices, :]
    # Singleton variable counts (a variable may appear in multiple singleton rows)
    fixed_counts = np.sum(A_sr_nonzero, axis=0)
    # Mask for singleton variables
    fixed_indices = fixed_counts >= 1

    # Determine singleton variable values by elimination from each singleton row
    # Note: several rows may eliminate the same variable.
    b_sr = b_eq[sing_row_indices]
    A_sr = A_eq[sing_row_indices, :][A_sr_nonzero]
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
        x_j = x_sr[j_sr == j]
        # If more than one, check if all values identical
        if len(x_j) > 1 and any(np.abs(x_j - x_j[0]) > tol):
            # Solution infeasible
            status = 2
            return (pp, None, status)
        x_fixed[k] = x_j[0]
        k = k + 1

    # Number of fixed variables
    n_fixed = np.sum(fixed_indices)

    # No fixed variables
    if n_fixed == 0:
        status = 6  # no change
        return (pp, None, status)

    # All variables fixed
    if n_fixed == len(c):
        # If all variables have been fixed, the problem is either solved or
        # infeasible.
        # The fixed variables results in 0 >= b_ub and 0 == b_eq.
        # So check if the transformed b_eq equals 0 (within tolerance range).
        # For b_ub there is no need to use a tolerance, just check if b_ub <= 0.
        # The tolerance to use for element i in b_eq is:
        # tol_b[i]^2 = A_eq[i,1]^2 * tol_x[1]^2 + ... + A_eq[i,N] * tol_x[N]^2
        # if there are N variables in x.
        # All x-tolerances are the same: tol_x[j] = tol.
        tol_b_eq = np.sqrt(np.sum(np.square(A_eq), 1)) * tol
        # Calculate residu A_eq * x - b_eq and check if zero
        r = np.dot(A_eq, x_fixed) - b_eq
        sol_eq = np.all(np.abs(r) <= tol_b_eq)
        if A_ub is not None:
            # Calculate residu A_ub * x - b_ub and check if >= zero
            r = np.dot(A_ub, x_fixed) - b_ub
            sol_ub = np.all(r >= 0)
        else:
            sol_ub = True
        if sol_eq and sol_ub:
            # Solution determined in presolve
            status = 0
        else:
            # Solution infeasible
            status = 2
        return (pp, None, status)

    # Remove fixed variables and singleton rows from equations
    keep_indices = np.logical_not(fixed_indices)
    if b_eq is not None:
        # Column removal and b_eq update
        b_eq = b_eq - np.dot(A_eq[:, fixed_indices], x_fixed)
        A_eq = A_eq[:, keep_indices]
        # Row removal
        keep_row_indices = np.logical_not(sing_row_indices)
        A_eq = A_eq[keep_row_indices, :]
        b_eq = b_eq[keep_row_indices]
    if b_ub is not None:
        b_ub = b_ub - np.dot(A_ub[:, fixed_indices], x_fixed)
        A_ub = A_ub[:, keep_indices]
    c_fixed = c[fixed_indices]
    c = c[keep_indices]
    bounds = bounds[keep_indices, :]
    if x0 is not None:
        x0_fixed = x0[fixed_indices]
        x0 = x0[keep_indices]

    def rev(c_mod, x_mod, x0_mod):
        # Reverse operation only relevant for c, x and x0
        # When removing elements at positions k1, k2, k3, ...
        # these must be replaced at (after) positions k1-1, k2-2, k3-3, ...
        # in the modified array
        i = np.flatnonzero(fixed_indices)
        # Number of variables to restore
        N = len(i)
        index_offset = list(range(N))
        # Create insert indices
        insert_indices = np.subtract(i, index_offset).flatten()
        c_rev = np.insert(c_mod.astype(float), insert_indices, c_fixed)
        x_rev = np.insert(x_mod.astype(float), insert_indices, x_fixed)
        if x0 is not None:
            x0_rev = np.insert(x0_mod.astype(float), insert_indices, x0_fixed)
        else:
            x0_rev = None
        return (c_rev, x_rev, x0_rev)

    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)  # TODO check if necessary
    status = 5  # modifications
    return (pp, rev, status)


def _presolve_remove_empty_rows(pp, tol=1e-9):
    """
    Check for empty rows in A_eq and A_ub.
    Remove these rows and the corresponding elements in b_eq and b_ub.
    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    """
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = pp

    # Equations
    if A_eq is None:
        zero_row_indices = []
    else:
        zero_row_indices = np.array(np.sum(A_eq != 0, axis=1) == 0)
    if np.any(zero_row_indices):
        if np.any(np.abs(b_eq[zero_row_indices]) > tol):
            # infeasible if RHS is not zero
            status = 2
            return (pp, None, status)
        else:
            # if RHS is zero, we can eliminate this equation entirely.
            nonzero_row_indices = np.logical_not(zero_row_indices)
            A_eq = A_eq[nonzero_row_indices, :]
            b_eq = b_eq[nonzero_row_indices]

    # Inequalities
    if A_ub is None:
        zero_row_indices = []
    else:
        zero_row_indices = np.array(np.sum(A_ub != 0, axis=1) == 0)
    if np.any(zero_row_indices):
        if np.any(b_ub[zero_row_indices] > tol):
            # A_ub * x >= b_ub, with empty row causing LHS to be 0.
            # Infeasible if RHS bigger than zero.
            status = 2
            return (pp, None, status)
        else:
            # if RHS is zero or below, we can eliminate this inequality.
            nonzero_row_indices = np.logical_not(zero_row_indices)
            A_ub = A_ub[nonzero_row_indices, :]
            b_ub = b_ub[nonzero_row_indices]

    # Wrap up
    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)  # TODO check if necessary
    status = 6  # modifications made, but nothing to reverse, so return 6, not 5
    return (pp, None, status)


def _presolve_remove_empty_columns(pp):
    """
    Check for empty columns in A_eq and A_ub.
    If a variable does not appear in any of the equations and
    inequalities, an optimal choice can be made for it using its bounds
    and the corrsponding value in c.
    Remove these empty columns and the corresponding variables.
    """
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = pp

    # Determine indices of zero columns (= indices of fixed variables)
    fixed_indices = np.full((len(c), ), True)
    if A_ub is not None:
        zero_coll_indices_ub = np.array(np.sum(A_ub != 0, axis=0) == 0)
        fixed_indices = np.logical_and(fixed_indices, zero_coll_indices_ub)
    if A_eq is not None:
        zero_coll_indices_eq = np.array(np.sum(A_eq != 0, axis=0) == 0)
        fixed_indices = np.logical_and(fixed_indices, zero_coll_indices_eq)

    # If no empty columns, return
    if not np.any(fixed_indices):
        status = 6
        return (pp, None, status)

    # Determine optimal values of fixed variables
    x_fixed = np.zeros((len(c), ))
    fixed_indices_neg_c = np.logical_and(fixed_indices, c < 0)
    x_fixed[fixed_indices_neg_c] = bounds[fixed_indices_neg_c, 1]
    fixed_indices_pos_c = np.logical_and(fixed_indices, c > 0)
    x_fixed[fixed_indices_pos_c] = bounds[fixed_indices_pos_c, 0]
    x_fixed = x_fixed[fixed_indices]
    if np.any(np.isinf(x_fixed)):
        status = 3
        return (pp, None, status)

    # Modify c, A_ub, A_eq, bounds, x0
    keep_indices = np.logical_not(fixed_indices)
    if A_eq is not None:
        # Column removal
        A_eq = A_eq[:, keep_indices]
    if A_ub is not None:
        A_ub = A_ub[:, keep_indices]
    c_fixed = c[fixed_indices]
    c = c[keep_indices]
    bounds = bounds[keep_indices, :]
    if x0 is not None:
        x0_fixed = x0[fixed_indices]
        x0 = x0[keep_indices]

    def rev(c_mod, x_mod, x0_mod):
        # Reverse operation only relevant for c, x and x0
        # When removing elements at positions k1, k2, k3, ...
        # these must be replaced at (after) positions k1-1, k2-2, k3-3, ...
        # in the modified array
        i = np.flatnonzero(fixed_indices)
        # Number of variables to restore
        N = len(i)
        index_offset = list(range(N))
        # Create insert indices
        insert_indices = np.subtract(i, index_offset).flatten()
        c_rev = np.insert(c_mod.astype(float), insert_indices, c_fixed)
        x_rev = np.insert(x_mod.astype(float), insert_indices, x_fixed)
        if x0 is not None:
            x0_rev = np.insert(x0_mod.astype(float), insert_indices, x0_fixed)
        else:
            x0_rev = None
        return (c_rev, x_rev, x0_rev)

    pp = (c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    status = 5
    return (pp, rev, status)


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


def _get_Abc(lp, c0, undo=[]):
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

    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem

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

    # bounds will be modified, create a copy
    bounds = np.array(bounds, copy=True)
    # undo[0] contains indices of variables removed from the problem
    # however, their bounds are still part of the bounds list
    # they are needed elsewhere, but not here
    if undo is not None and undo != []:
        bounds = np.delete(bounds, undo[0], 0)

    # modify problem such that all variables have only non-negativity bounds
    lbs = bounds[:, 0]
    ubs = bounds[:, 1]
    m_ub, n_ub = A_ub.shape

    lb_none = np.equal(lbs, -np.inf)
    ub_none = np.equal(ubs, np.inf)
    lb_some = np.logical_not(lb_none)
    ub_some = np.logical_not(ub_none)

    # if preprocessing is on, lb == ub can't happen
    # if preprocessing is off, then it would be best to convert that
    # to an equality constraint, but it's tricky to make the other
    # required modifications from inside here.

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
    c[n_ub:n_ub+n_free] = -c[i_free]
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
        R = 1/_round_to_power_of_two(R)
        A = sps.diags(R)*A if sps.issparse(A) else A*R.reshape(m, 1)
        b = b*R

        C = np.max(np.abs(A), axis=0)
        if sps.issparse(A):
            C = C.toarray().flatten()
        C[C == 0] = 1
        C = 1/_round_to_power_of_two(C)
        A = A*sps.diags(C) if sps.issparse(A) else A*C
        c = c*C

    b_scale = np.max(np.abs(b)) if b.size > 0 else 1
    if b_scale == 0:
        b_scale = 1.
    b = b/b_scale

    if x0 is not None:
        x0 = x0/b_scale*(1/C)
    return A, b, c, x0, C, b_scale


def _unscale(x, C, b_scale):
    """
    Converts solution to _autoscale problem -> solution to original problem.
    """

    try:
        n = len(C)
        # fails if sparse or scalar; that's OK.
        # this is only needed for original simplex (never sparse)
    except TypeError:
        n = len(x)

    return x[:n]*b_scale*C


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


def _postsolve(x, postsolve_args, complete=False, tol=1e-8, copy=False):
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

    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

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
    bounds : 2D array
        The bounds on the original variables ``x``
    """
    # note that all the inputs are the ORIGINAL, unmodified versions
    # no rows, columns have been removed
    # the only exception is bounds; it has been modified
    # we need these modified values to undo the variable substitutions
    # in retrospect, perhaps this could have been simplified if the "undo"
    # variable also contained information for undoing variable substitutions

    (c, A_ub, b_ub, A_eq, b_eq, bounds, x0), undo, C, b_scale = postsolve_args
    x = _unscale(x, C, b_scale)

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
        copy = True
    if copy:
        x = np.array(x, copy=True)

    # now undo variable substitutions
    # if "complete", problem was solved in presolve; don't do anything here
    if not complete and bounds is not None:  # bounds are never none, probably
        n_unbounded = 0
        for i, bi in enumerate(bounds):
            if i in no_adjust:
                continue
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

    n_x = len(c)
    x = x[:n_x]  # all the rest of the variables were artificial
    fun = x.dot(c)
    slack = b_ub - A_ub.dot(x)  # report slack for ORIGINAL UB constraints
    # report residuals of ORIGINAL EQ constraints
    con = b_eq - A_eq.dot(x)

    return x, fun, slack, con, bounds


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
        np.isnan(x).any()
        or np.isnan(fun)
        or np.isnan(slack).any()
        or np.isnan(con).any()
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


def _postprocess(x, postsolve_args, complete=False, status=0, message="",
                 tol=1e-8, iteration=None, disp=False):
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
    c : 1-D array
        Original coefficients of the linear objective function to be minimized.
    A_ub : 2-D array, optional
        2-D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1-D array, optional
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2-D array, optional
        2-D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1-D array, optional
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : 2D array
        The bounds of ``x``, lower bounds in the 1st column, upper
        bounds in the 2nd column. The bounds are possibly tightened
        by the presolve procedure.
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
             4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

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

    x, fun, slack, con, bounds = _postsolve(
        x, postsolve_args, complete, tol
    )

    status, message = _check_result(
        x, fun, status, slack, con,
        bounds, tol, message
    )

    if disp:
        _display_summary(message, status, fun, iteration)

    return x, fun, slack, con, status, message
