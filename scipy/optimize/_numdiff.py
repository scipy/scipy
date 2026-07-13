"""Routines for numerical differentiation."""
import functools
import numpy as np

from scipy.sparse.linalg import LinearOperator
from ..sparse import issparse, spmatrix, find, csc_array, csr_array, csr_matrix
from ._group_columns import group_dense, group_sparse
from scipy._lib._array_api import (
    array_namespace, xp_result_type, xp_default_dtype, xp_copy, xp_capabilities,
    xp_size, xp_vector_norm
)
from scipy._lib._util import MapWrapper, xp_array_equal
from scipy._external.array_api_compat import is_numpy_namespace
from scipy._external import array_api_extra as xpx


def _adjust_scheme_to_bounds(x0, h, num_steps, scheme, lb, ub, xp=None):
    """Adjust final difference scheme to the presence of bounds.

    Parameters
    ----------
    x0 : ndarray, shape (n,)
        Point at which we wish to estimate derivative.
    h : ndarray, shape (n,)
        Desired absolute finite difference steps.
    num_steps : int
        Number of `h` steps in one direction required to implement finite
        difference scheme. For example, 2 means that we need to evaluate
        f(x0 + 2 * h) or f(x0 - 2 * h)
    scheme : {'1-sided', '2-sided'}
        Whether steps in one or both directions are required. In other
        words '1-sided' applies to forward and backward schemes, '2-sided'
        applies to center schemes.
    lb : ndarray, shape (n,)
        Lower bounds on independent variables.
    ub : ndarray, shape (n,)
        Upper bounds on independent variables.

    Returns
    -------
    h_adjusted : ndarray, shape (n,)
        Adjusted absolute step sizes. Step size decreases only if a sign flip
        or switching to one-sided scheme doesn't allow to take a full step.
    use_one_sided : ndarray of bool, shape (n,)
        Whether to switch to one-sided scheme. Informative only for
        ``scheme='2-sided'``.
    """
    xp = xp or array_namespace(x0)

    if scheme == '1-sided':
        use_one_sided = xp.ones_like(h, dtype=xp.bool)
    elif scheme == '2-sided':
        h = xp.abs(h)
        use_one_sided = xp.zeros_like(h, dtype=xp.bool)
    else:
        raise ValueError("`scheme` must be '1-sided' or '2-sided'.")

    if xp.all((lb == -xp.inf) & (ub == xp.inf)):
        return h, use_one_sided

    h_total = h * num_steps
    h_adjusted = xp_copy(h, xp=xp)

    lower_dist = x0 - lb
    upper_dist = ub - x0

    if scheme == '1-sided':
        x = x0 + h_total
        violated = (x < lb) | (x > ub)
        fitting = xp.abs(h_total) <= xp.maximum(lower_dist, upper_dist)
        h_adjusted[violated & fitting] *= -1

        forward = (upper_dist >= lower_dist) & ~fitting
        h_adjusted[forward] = upper_dist[forward] / num_steps
        backward = (upper_dist < lower_dist) & ~fitting
        h_adjusted[backward] = -lower_dist[backward] / num_steps
    elif scheme == '2-sided':
        central = (lower_dist >= h_total) & (upper_dist >= h_total)

        forward = (upper_dist >= lower_dist) & ~central
        h_adjusted[forward] = xp.minimum(
            h[forward], 0.5 * upper_dist[forward] / num_steps)
        use_one_sided[forward] = True

        backward = (upper_dist < lower_dist) & ~central
        h_adjusted[backward] = -xp.minimum(
            h[backward], 0.5 * lower_dist[backward] / num_steps)
        use_one_sided[backward] = True

        min_dist = xp.minimum(upper_dist, lower_dist) / num_steps
        adjusted_central = (~central & (xp.abs(h_adjusted) <= min_dist))
        h_adjusted[adjusted_central] = min_dist[adjusted_central]
        use_one_sided[adjusted_central] = False

    return h_adjusted, use_one_sided


@xp_capabilities(
    skip_backends=[
        (
            "jax.numpy",
            "needs jax control flow logic (if -> xp.where) and conversion to be "
            "functionally pure"
        ),
    ]
)
@functools.lru_cache
def _eps_for_method(x0_dtype, f0_dtype, method, xp):
    """
    Calculates relative EPS step to use for a given data type
    and numdiff step method.

    Progressively smaller steps are used for larger floating point types.

    Parameters
    ----------
    f0_dtype: dtype
        function evaluation

    x0_dtype: dtype
        parameter vector

    method: {'2-point', '3-point', 'cs'}

    xp: array-namespace

    Returns
    -------
    EPS: float
        relative step size. May be xp.float16, xp.float32, xp.float64

    Notes
    -----
    The default relative step will be xp.float64. However, if x0 or f0 are
    smaller floating point types (xp.float16, xp.float32), then the smallest
    floating point type is chosen.
    """
    # the default EPS value. Normally we want this to be 64 bit, but torch
    # has a 32 bit default.
    out_dtype = xp_default_dtype(xp)
    EPS = xp.finfo(out_dtype).eps

    bits = {}
    for _bits in (16, 32, 64):
        if hasattr(xp, f"float{_bits}"):
            bits[_bits] = getattr(xp, f"float{_bits}")

    x0_is_fp = False
    if xp.isdtype(x0_dtype, 'real floating'):
        # if you're a floating point type then over-ride the default EPS
        EPS = xp.finfo(x0_dtype).eps
        x0_size = xp.finfo(x0_dtype).bits
        out_dtype = bits[x0_size]
        x0_is_fp = True

    if xp.isdtype(f0_dtype, 'real floating'):
        f0_size = xp.finfo(f0_dtype).bits
        # choose the smallest itemsize between x0 and f0
        if x0_is_fp and f0_size < x0_size:
            EPS = xp.finfo(f0_dtype).eps
            out_dtype = bits[f0_size]

    if method in ["2-point", "cs"]:
        return xp.asarray(EPS**0.5, dtype=out_dtype)
    elif method in ["3-point"]:
        return xp.asarray(EPS**(1/3), dtype=out_dtype)
    else:
        raise RuntimeError("Unknown step method, should be one of "
                           "{'2-point', '3-point', 'cs'}")


@xp_capabilities(
    skip_backends=[
        (
            "jax.numpy",
            "needs jax control flow logic (if -> xp.where) and conversion to be "
            "functionally pure"
        ),
    ]
)
def _compute_absolute_step(rel_step, x0, f0, method, xp=None):
    """
    Computes an absolute step from a relative step for finite difference
    calculation.

    Parameters
    ----------
    rel_step: None or array-like
        Relative step for the finite difference calculation
    x0 : np.ndarray
        Parameter vector
    f0 : np.ndarray or scalar
    method : {'2-point', '3-point', 'cs'}
    xp : array-namespace

    Returns
    -------
    h : np.array
        The absolute step size, ``h.dtype==x0.dtype``.

    Notes
    -----
    `h` has the same dtype as `x0` because dx is later calculated
    as ``(x0 + h) - x0``, and problems would occur if ``x0.dtype==np.float16``
    with ``h.dtype==np.float64``.
    If `rel_step is None`, then a default relative step is calculated using the
    smallest floating point type of `x0` and `f0`, see _eps_for_method.

    """
    # this is used instead of np.sign(x0) because we need
    # sign_x0 to be 1 when x0 == 0.
    xp = xp or array_namespace(x0)
    gte0 = (x0 >= 0)
    sign_x0 = xp.astype(gte0, x0.dtype) * 2 - 1

    rstep = _eps_for_method(x0.dtype, f0.dtype, method, xp=xp)
    default_abs_step = xp.astype(
        rstep * sign_x0 * xp.maximum(xp.asarray(1.0, dtype=x0.dtype), xp.abs(x0)),
        x0.dtype
    )

    if rel_step is None:
        abs_step = default_abs_step
    else:
        # User has requested specific relative steps.
        # Don't multiply by max(1, abs(x0) because if x0 < 1 then their
        # requested step is not used.
        _abs_step = rel_step * sign_x0 * xp.abs(x0)
        abs_step = xp.asarray(_abs_step, dtype=x0.dtype)

        # however we don't want an abs_step of 0, which can happen if
        # rel_step is 0, or x0 is 0. Instead, substitute a realistic step
        dx = ((x0 + abs_step) - x0)
        abs_step = xp.where(
            dx == 0,
            default_abs_step,
            abs_step
        )

    return abs_step


def _prepare_bounds(bounds, x0, xp=None):
    """
    Prepares new-style bounds from a two-tuple specifying the lower and upper
    limits for values in x0. If a value is not bound then the lower/upper bound
    will be expected to be -np.inf/np.inf.

    Examples
    --------
    >>> _prepare_bounds([(0, 1, 2), (1, 2, np.inf)], [0.5, 1.5, 2.5])
    (array([0., 1., 2.]), array([ 1.,  2., inf]))
    """
    xp = xp or array_namespace(x0)

    lb, ub = (xp.asarray(b, dtype=x0.dtype) for b in bounds)
    lb = xp.broadcast_to(lb, x0.shape)
    ub = xp.broadcast_to(ub, x0.shape)

    return lb, ub


def group_columns(A, order=0):
    """Group columns of a 2-D matrix for sparse finite differencing [1]_.

    Two columns are in the same group if in each row at least one of them
    has zero. A greedy sequential algorithm is used to construct groups.

    Parameters
    ----------
    A : array_like or sparse array, shape (m, n)
        Matrix of which to group columns.
    order : int, iterable of int with shape (n,) or None
        Permutation array which defines the order of columns enumeration.
        If int or None, a random permutation is used with `order` used as
        a random seed. Default is 0, that is use a random permutation but
        guarantee repeatability.

    Returns
    -------
    groups : ndarray of int, shape (n,)
        Contains values from 0 to n_groups-1, where n_groups is the number
        of found groups. Each value ``groups[i]`` is an index of a group to
        which ith column assigned. The procedure was helpful only if
        n_groups is significantly less than n.

    References
    ----------
    .. [1] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
           sparse Jacobian matrices", Journal of the Institute of Mathematics
           and its Applications, 13 (1974), pp. 117-120.
    """
    if issparse(A):
        A = csc_array(A)
    else:
        A = np.atleast_2d(A)
        A = (A != 0).astype(np.int32)

    if A.ndim != 2:
        raise ValueError("`A` must be 2-dimensional.")

    m, n = A.shape

    if order is None or np.isscalar(order):
        rng = np.random.RandomState(order)
        order = rng.permutation(n)
    else:
        order = np.asarray(order)
        if order.shape != (n,):
            raise ValueError("`order` has incorrect shape.")

    A = A[:, order]

    if issparse(A):
        groups = group_sparse(m, n, A.indices, A.indptr)
    else:
        groups = group_dense(m, n, A)

    groups[order] = groups.copy()

    return groups


@xp_capabilities(
    skip_backends=[
        ("dask.array", "would need lazy_xp_function to avoid dask.compute"),
        (
            "jax.numpy",
            "needs jax control flow logic (if -> xp.where) and conversion to be "
            "functionally pure"
        ),
    ]
)
def approx_derivative(fun, x0, method='3-point', rel_step=None, abs_step=None,
                      f0=None, bounds=(-np.inf, np.inf), sparsity=None,
                      as_linear_operator=False, args=(), kwargs=None,
                      full_output=False, workers=None):
    """Compute finite difference approximation of the derivatives of a
    vector-valued function.

    If a function maps from R^n to R^m, its derivatives form m-by-n matrix
    called the Jacobian, where an element (i, j) is a partial derivative of
    f[i] with respect to x[j].

    Parameters
    ----------
    fun : callable
        Function of which to estimate the derivatives. The argument x
        passed to this function is ndarray of shape (n,) (never a scalar
        even if n=1). It must return 1-D array_like of shape (m,) or a scalar.
    x0 : array_like of shape (n,) or float
        Point at which to estimate the derivatives. Float will be converted
        to a 1-D array.
    method : {'3-point', '2-point', 'cs'}, optional
        Finite difference method to use:
            - '2-point' - use the first order accuracy forward or backward
                          difference.
            - '3-point' - use central difference in interior points and the
                          second order accuracy forward or backward difference
                          near the boundary.
            - 'cs' - use a complex-step finite difference scheme. This assumes
                     that the user function is real-valued and can be
                     analytically continued to the complex plane. Otherwise,
                     produces bogus results.
    rel_step : None or array_like, optional
        Relative step size to use. If None (default) the absolute step size is
        computed as ``h = (rel_step * sign(x0) * max(1, abs(x0))).astype(x0.dtype)``,
        with `rel_step` being selected automatically, see Notes. Otherwise
        ``h = (rel_step * sign(x0) * abs(x0)).astype(x0.dtype)``. For
        ``method='3-point'`` the sign of `h` is ignored. The calculated step size
        is possibly adjusted to fit into the bounds.
    abs_step : array_like, optional
        Absolute step size to use, possibly adjusted to fit into the bounds.
        For ``method='3-point'`` the sign of `abs_step` is ignored. By default
        relative steps are used, only if ``abs_step is not None`` are absolute
        steps used. `abs_step` is coerced to the dtype of `x0` during calculation.
    f0 : None or array_like, optional
        If not None it is assumed to be equal to ``fun(x0)``, in this case
        the ``fun(x0)`` is not called. Default is None.
    bounds : tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each bound must match the size of `x0` or be a scalar, in the latter
        case the bound will be the same for all variables. Use it to limit the
        range of function evaluation. Bounds checking is not implemented
        when `as_linear_operator` is True.
    sparsity : {None, array_like, sparse array, 2-tuple}, optional
        Defines a sparsity structure of the Jacobian matrix. If the Jacobian
        matrix is known to have only few non-zero elements in each row, then
        it's possible to estimate its several columns by a single function
        evaluation [3]_. To perform such economic computations two ingredients
        are required:

        * structure : array_like or sparse array of shape (m, n). A zero
          element means that a corresponding element of the Jacobian
          identically equals to zero.
        * groups : array_like of shape (n,). A column grouping for a given
          sparsity structure, use `group_columns` to obtain it.

        A single array or a sparse array is interpreted as a sparsity
        structure, and groups are computed inside the function. A tuple is
        interpreted as (structure, groups). If None (default), a standard
        dense differencing will be used.

        Note, that sparse differencing makes sense only for large Jacobian
        matrices where each row contains few non-zero elements.
    as_linear_operator : bool, optional
        When True the function returns an `scipy.sparse.linalg.LinearOperator`.
        Otherwise it returns a dense array or a sparse array depending on
        `sparsity`. The linear operator provides an efficient way of computing
        ``J.dot(p)`` for any vector ``p`` of shape (n,), but does not allow
        direct access to individual elements of the matrix. By default
        `as_linear_operator` is False.
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)``.
    full_output : bool, optional
        If True then the function also returns a dictionary with extra information
        about the calculation.
    workers : int or map-like callable, optional
        Supply a map-like callable, such as
        `multiprocessing.Pool.map` for evaluating the population in parallel.
        This evaluation is carried out as ``workers(fun, iterable)``.
        Alternatively, if `workers` is an int the task is subdivided into `workers`
        sections and the fun evaluated in parallel
        (uses `multiprocessing.Pool <multiprocessing>`).
        Supply -1 to use all available CPU cores.
        It is recommended that a map-like be used instead of int, as repeated
        calls to `approx_derivative` will incur large overhead from setting up
        new processes.

    Returns
    -------
    J : {ndarray, sparse array, LinearOperator}
        Finite difference approximation of the Jacobian matrix.
        If `as_linear_operator` is True returns a LinearOperator
        with shape (m, n). Otherwise it returns a dense array or sparse
        array depending on how `sparsity` is defined. If `sparsity`
        is None then an ndarray with shape (m, n) is returned. If
        `sparsity` is not None returns a csr_array or csr_matrix with
        shape (m, n) following the array/matrix type of the incoming structure.
        For sparse arrays and linear operators it is always returned as
        a 2-D structure. For ndarrays, if m=1 it is returned
        as a 1-D gradient array with shape (n,).
        `J.dtype` is the promoted type of ``fun(x0)`` and `x0`.

    info_dict : dict
        Dictionary containing extra information about the calculation. The
        keys include:

        - `nfev`, number of function evaluations. If `as_linear_operator` is True
           then `fun` is expected to track the number of evaluations itself.
           This is because multiple calls may be made to the linear operator which
           are not trackable here.

    See Also
    --------
    check_derivative : Check correctness of a function computing derivatives.

    Notes
    -----
    If `rel_step` is not provided, it is assigned as ``EPS**(1/s)``, where EPS is
    determined from the smallest floating point dtype of `x0` or ``fun(x0)``,
    ``np.finfo(x0.dtype).eps``, s=2 for '2-point' method and
    s=3 for '3-point' method. This relative step approximately minimizes a sum
    of truncation and round-off errors, see [1]_. Relative steps are used by
    default. However, absolute steps are used when ``abs_step is not None``.
    If any of the absolute or relative steps produces an indistinguishable
    difference from the original `x0`, ``(x0 + dx) - x0 == 0``, then a
    automatic step size is substituted for that particular entry.
    The calculated absolute step size, `h`, is coerced to have the same dtype as `x0`.

    The floating point precision of `J` is the promoted type of ``fun(x0)``
    and `x0`. If `rel_step` and `abs_step` are None, then the overall accuracy of
    `J` depends on the smallest floating point dtype of `x0` and ``fun(x0)``,
    as that dtype is used to calculate the default step size.

    A finite difference scheme for '3-point' method is selected automatically.
    The well-known central difference scheme is used for points sufficiently
    far from the boundary, and 3-point forward or backward scheme is used for
    points near the boundary. Both schemes have the second-order accuracy in
    terms of Taylor expansion. Refer to [2]_ for the formulas of 3-point
    forward and backward difference schemes.

    For dense differencing when m=1 Jacobian is returned with a shape (n,),
    on the other hand when n=1 Jacobian is returned with a shape (m, 1).
    Our motivation is the following: a) It handles a case of gradient
    computation (m=1) in a conventional way. b) It clearly separates these two
    different cases. b) In all cases np.atleast_2d can be called to get 2-D
    Jacobian with correct dimensions.

    References
    ----------
    .. [1] W. H. Press et. al. "Numerical Recipes. The Art of Scientific
           Computing. 3rd edition", sec. 5.7.

    .. [2] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
           sparse Jacobian matrices", Journal of the Institute of Mathematics
           and its Applications, 13 (1974), pp. 117-120.

    .. [3] B. Fornberg, "Generation of Finite Difference Formulas on
           Arbitrarily Spaced Grids", Mathematics of Computation 51, 1988.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize._numdiff import approx_derivative
    >>>
    >>> def f(x, c1, c2):
    ...     return np.array([x[0] * np.sin(c1 * x[1]),
    ...                      x[0] * np.cos(c2 * x[1])])
    ...
    >>> x0 = np.array([1.0, 0.5 * np.pi])
    >>> approx_derivative(f, x0, args=(1, 2))
    array([[ 1.,  0.],
           [-1.,  0.]])

    Bounds can be used to limit the region of function evaluation.
    In the example below we compute left and right derivative at point 1.0.

    >>> def g(x):
    ...     return x**2 if x >= 1 else x
    ...
    >>> x0 = 1.0
    >>> approx_derivative(g, x0, bounds=(-np.inf, 1.0))
    array([ 1.])
    >>> approx_derivative(g, x0, bounds=(1.0, np.inf))
    array([ 2.])

    We can also parallelize the derivative calculation using the workers
    keyword.

    >>> from multiprocessing import Pool
    >>> import time
    >>> def fun2(x):       # import from an external file for use with multiprocessing
    ...     time.sleep(0.002)
    ...     return rosen(x)

    >>> rng = np.random.default_rng()
    >>> x0 = rng.uniform(high=10, size=(2000,))
    >>> f0 = rosen(x0)

    >>> %timeit approx_derivative(fun2, x0, f0=f0)     # may vary
    10.5 s ± 5.91 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

    >>> elapsed = []
    >>> with Pool() as workers:
    ...     for i in range(10):
    ...         t = time.perf_counter()
    ...         approx_derivative(fun2, x0, workers=workers.map, f0=f0)
    ...         et = time.perf_counter()
    ...         elapsed.append(et - t)
    >>> np.mean(elapsed)    # may vary
    np.float64(1.442545195999901)

    Create a map-like vectorized version. `x` is a generator, so first of all
    a 2-D array, `xx`, is reconstituted. Here `xx` has shape `(Y, N)` where `Y`
    is the number of function evaluations to perform and `N` is the dimensionality
    of the objective function. The underlying objective function is `rosen`, which
    requires `xx` to have shape `(N, Y)`, so a transpose is required.

    >>> def fun(f, x, *args, **kwds):
    ...     xx = np.r_[[xs for xs in x]]
    ...     return f(xx.T)
    >>> %timeit approx_derivative(fun2, x0, workers=fun, f0=f0)    # may vary
    91.8 ms ± 755 μs per loop (mean ± std. dev. of 7 runs, 10 loops each)

    """
    if method not in ['2-point', '3-point', 'cs']:
        raise ValueError(f"Unknown method '{method}'. ")

    info_dict = {'nfev': None}

    xp = array_namespace(x0)
    _x = xpx.atleast_nd(xp.asarray(x0), ndim=1, xp=xp)
    _dtype = xp.float64
    if xp.isdtype(_x.dtype, "real floating"):
        _dtype = _x.dtype

    # promotes to floating
    x0 = xp.astype(_x, _dtype)

    if x0.ndim > 1:
        raise ValueError("`x0` must have at most 1 dimension.")

    # cast the bounds to the relevant array_namespace
    lb, ub = _prepare_bounds(bounds, x0, xp=xp)

    if lb.shape != x0.shape or ub.shape != x0.shape:
        raise ValueError("Inconsistent shapes between bounds and `x0`.")

    if not is_numpy_namespace(xp) and (sparsity is not None):
        raise RuntimeError("as_linear_operator and sparsity are not yet supported by"
                           "non-numpy array namespaces")

    if as_linear_operator and not (xp.all(xp.isinf(lb))
                                   and xp.all(xp.isinf(ub))):
        raise ValueError("Bounds not supported when "
                         "`as_linear_operator` is True.")

    if kwargs is None:
        kwargs = {}

    fun_wrapped = _Fun_Wrapper(fun, x0, args, kwargs)

    # Record how function evaluations are consumed by `approx_derivative`.
    # Historically this was done by upstream functions wrapping `fun`.
    # However, with parallelization via workers it was going to be impossible to
    # keep that counter updated across Processes. Counter synchronisation can
    # be achieved via multiprocessing.Value and a Pool. However, workers can be
    # any map-like, not necessarily a Pool, so initialization of the Value would
    # be difficult.
    nfev = _nfev = 0

    if f0 is None:
        f0 = fun_wrapped(x0)
        nfev = 1
    else:
        f0 = xpx.atleast_nd(f0, ndim=1, xp=xp)
        if f0.ndim > 1:
            raise ValueError("`f0` passed has more than 1 dimension.")

    if xp.any((x0 < lb) | (x0 > ub)):
        raise ValueError("`x0` violates bound constraints.")

    if as_linear_operator:
        if rel_step is None:
            rel_step = _eps_for_method(x0.dtype, f0.dtype, method, xp=xp)

        J, _ = _linear_operator_difference(fun_wrapped, x0,
                                           f0, rel_step, method)
    else:
        # by default we use rel_step
        if abs_step is None:
            h = _compute_absolute_step(rel_step, x0, f0, method, xp=xp)
        else:
            # user specifies an absolute step
            sign_x0 = xp.astype(x0 >= 0., x0.dtype) * 2 - 1
            h = abs_step

            # cannot have a zero step. This might happen if x0 is very large
            # or small. In which case fall back to relative step.
            dx = ((x0 + h) - x0)
            h = xp.where(dx == 0,
                         _eps_for_method(x0.dtype, f0.dtype, method, xp=xp) *
                         sign_x0 *
                         xp.maximum(xp.asarray(1.0, dtype=x0.dtype), xp.abs(x0)),
                         h)
            h = xp.astype(h, x0.dtype)

        if method == '2-point':
            h, use_one_sided = _adjust_scheme_to_bounds(
                x0, h, 1, '1-sided', lb, ub)
        elif method == '3-point':
            h, use_one_sided = _adjust_scheme_to_bounds(
                x0, h, 1, '2-sided', lb, ub)
        elif method == 'cs':
            use_one_sided = False

        # normalize workers
        workers = workers or map
        with MapWrapper(workers) as mf:
            if sparsity is None:
                J, _nfev = _dense_difference(fun_wrapped, x0, f0, h,
                                         use_one_sided, method,
                                         mf)
            else:
                if not issparse(sparsity) and len(sparsity) == 2:
                    structure, groups = sparsity
                else:
                    structure = sparsity
                    groups = group_columns(sparsity)

                if issparse(structure):
                    structure = structure.tocsc()
                else:
                    structure = np.atleast_2d(structure)
                groups = np.atleast_1d(groups)
                J, _nfev = _sparse_difference(fun_wrapped, x0, f0, h,
                                              use_one_sided, structure,
                                              groups, method, mf, xp=xp)

    if full_output:
        nfev += _nfev
        info_dict["nfev"] = nfev
        return J, info_dict
    else:
        return J


def _linear_operator_difference(fun, x0, f0, h, method, xp=None):
    m = xp_size(f0)
    n = xp_size(x0)

    xp = xp or array_namespace(x0)
    result_dtype = xp_result_type(x0, f0, force_floating=True, xp=xp)

    if method == '2-point':
        # nfev = 1
        def matvec(p):
            if xp_array_equal(p, xp.zeros_like(p), xp):
                return xp.zeros(m, dtype=result_dtype)
            dx = h / xp_vector_norm(p, xp=xp)
            x = x0 + dx*p
            df = fun(x) - f0
            return df / dx

    elif method == '3-point':
        # nfev = 2
        def matvec(p):
            if xp_array_equal(p, xp.zeros_like(p), xp):
                return xp.zeros(m, dtype=result_dtype)
            dx = 2*h / xp_vector_norm(p, xp=xp)
            x1 = x0 - (dx/2)*p
            x2 = x0 + (dx/2)*p
            f1 = fun(x1)
            f2 = fun(x2)
            df = f2 - f1
            return df / dx

    elif method == 'cs':
        # nfev = 1
        def matvec(p):
            if xp_array_equal(p, xp.zeros_like(p), xp):
                return xp.zeros(m, dtype=result_dtype)
            dx = h / xp_vector_norm(p, xp=xp)
            x = x0 + dx*p*1.j
            f1 = fun(x)
            df = xp.imag(f1)
            return df / dx
    else:
        raise RuntimeError("Never be here.")

    return LinearOperator(shape=(m, n), matvec=matvec, dtype=result_dtype, xp=xp), 0


def _dense_difference(fun, x0, f0, h, use_one_sided, method, workers, xp=None):
    m = xp_size(f0)
    n = xp_size(x0)
    nfev = 0

    xp = xp or array_namespace(x0)

    # h should have same dtype as x0
    result_type = xp_result_type(x0, f0, force_floating=True, xp=xp)
    # output dtype should be the same as df_dx
    J_transposed = xp.empty((n, m), dtype=result_type)

    if method == '2-point':
        def x_generator2(x0, h):
            for i in range(n):
                # If copying isn't done then it's possible for different workers
                # to see the same values of x1. (At least that's what happened
                # when I used `multiprocessing.dummy.Pool`).
                # I also considered creating all the vectors at once, but that
                # means assembling a very large N x N array. It's therefore a
                # trade-off between N array copies or creating an NxN array.
                x1 = xp_copy(x0, xp=xp)
                x1[i] = x0[i] + h[i]
                yield x1

        # only f_evals (numerator) needs parallelization, the denominator
        # (the step size) is fast to calculate.
        f_evals = workers(fun, x_generator2(x0, h))
        dx = [(x0[i] + h[i]) - x0[i] for i in range(n)]
        df = [f_eval - f0 for f_eval in f_evals]
        df_dx = [delf / delx for delf, delx in zip(df, dx)]
        nfev += len(df_dx)

    elif method == '3-point':
        def x_generator3(x0, h, use_one_sided):
            for i, one_sided in enumerate(use_one_sided):
                x1 = xp_copy(x0, xp=xp)
                x2 = xp_copy(x0, xp=xp)
                if one_sided:
                    x1[i] = x0[i] + h[i]
                    x2[i] = x0[i] + 2*h[i]
                else:
                    x1[i] = x0[i] - h[i]
                    x2[i] = x0[i] + h[i]
                yield x1
                yield x2

        # workers may return something like a list that needs to be turned
        # into an iterable (can't call `next` on a list)
        f_evals = iter(workers(fun, x_generator3(x0, h, use_one_sided)))
        gen = x_generator3(x0, h, use_one_sided)
        dx = list()
        df = list()
        for i, one_sided in enumerate(use_one_sided):
            l = next(gen)
            u = next(gen)

            f1 = next(f_evals)
            f2 = next(f_evals)
            if one_sided:
                dx.append(u[i] - x0[i])
                df.append(-3.0 * f0 + 4 * f1 - f2)
            else:
                dx.append(u[i] - l[i])
                df.append(f2 - f1)
        df_dx = [delf / delx for delf, delx in zip(df, dx)]
        nfev += 2 * len(df_dx)
    elif method == 'cs':
        bits = xp.finfo(x0).bits
        xc_dtype = getattr(xp, f"complex{bits*2}")

        def x_generator_cs(x0, h):
            for i in range(n):
                xc = xp.astype(x0, xc_dtype)
                xc[i] += h[i] * 1.0j
                yield xc

        f_evals = iter(workers(fun, x_generator_cs(x0, h)))
        df_dx = [xp.imag(f1) / hi for f1, hi in zip(f_evals, h)]
        nfev += len(df_dx)
    else:
        raise RuntimeError("Never be here.")

    for i, v in enumerate(df_dx):
        J_transposed[i, ...] = v

    if m == 1:
        # flatten the array
        return xp.reshape(J_transposed, (-1,)), nfev
    else:
        # fun(x0) has more than 1-dimension
        return J_transposed.T, nfev


def _sparse_difference(fun, x0, f0, h, use_one_sided,
                       structure, groups, method, workers, xp=None):
    xp = xp or array_namespace(x0)
    m = xp.size(f0)
    n = xp_size(x0)
    row_indices = []
    col_indices = []
    fractions = []
    result_type = xp_result_type(x0, f0, force_floating=True, xp=xp)

    n_groups = xp.max(groups) + 1
    nfev = 0

    def e_generator():
        # Perturb variables which are in the same group simultaneously.
        for group in range(n_groups):
            yield xp.equal(group, groups)

    def x_generator2():
        e_gen = e_generator()
        for e in e_gen:
            h_vec = h * e
            x = x0 + h_vec
            yield x

    def x_generator3():
        e_gen = e_generator()
        for e in e_gen:
            h_vec = h * e
            x1 = xp_copy(x0, xp=xp)
            x2 = xp_copy(x0, xp=xp)

            mask_1 = use_one_sided & e
            x1[mask_1] += h_vec[mask_1]
            x2[mask_1] += 2 * h_vec[mask_1]

            mask_2 = ~use_one_sided & e
            x1[mask_2] -= h_vec[mask_2]
            x2[mask_2] += h_vec[mask_2]
            yield x1
            yield x2

    def x_generator_cs():
        e_gen = e_generator()
        for e in e_gen:
            h_vec = h * e
            yield x0 + h_vec * 1.j

    # evaluate the function for each of the groups
    if method == '2-point':
        f_evals = iter(workers(fun, x_generator2()))
        xs = x_generator2()
    elif method == '3-point':
        f_evals = iter(workers(fun, x_generator3()))
        xs = x_generator3()
    elif method == 'cs':
        f_evals = iter(workers(fun, x_generator_cs()))

    for e in e_generator():
        # The result is written to columns which correspond to perturbed
        # variables.
        cols, = xp.nonzero(e)
        # Find all non-zero elements in selected columns of Jacobian.
        i, j, _ = find(structure[:, cols])
        # Restore column indices in the full array.
        j = cols[j]

        if method == '2-point':
            dx = next(xs) - x0
            df = next(f_evals) - f0
            nfev += 1
        elif method == '3-point':
            # Here we do conceptually the same but separate one-sided
            # and two-sided schemes.
            x1 = next(xs)
            x2 = next(xs)

            mask_1 = use_one_sided & e
            mask_2 = ~use_one_sided & e

            dx = xp.zeros(n, dtype=x0.dtype)
            dx[mask_1] = x2[mask_1] - x0[mask_1]
            dx[mask_2] = x2[mask_2] - x1[mask_2]

            f1 = next(f_evals)
            f2 = next(f_evals)
            nfev += 2

            mask = use_one_sided[j]
            df = xp.empty(m, dtype=f0.dtype)

            rows = i[mask]
            df[rows] = -3 * f0[rows] + 4 * f1[rows] - f2[rows]

            rows = i[~mask]
            df[rows] = f2[rows] - f1[rows]
        elif method == 'cs':
            f1 = next(f_evals)
            nfev += 1
            df = f1.imag
            dx = h * e
        else:
            raise ValueError("Never be here.")

        # All that's left is to compute the fraction. We store i, j and
        # fractions as separate arrays and later construct csr_array.
        row_indices.append(i)
        col_indices.append(j)
        fractions.append(df[i] / dx[j])

    row_indices = xp.concat(row_indices, axis=0)
    col_indices = xp.concat(col_indices, axis=0)
    fractions = xp.concat(fractions, axis=0)

    if isinstance(structure, spmatrix):
        return csr_matrix(
            (fractions, (row_indices, col_indices)),
            shape=(m, n),
            dtype=result_type
        ), nfev
    return csr_array(
        (fractions, (row_indices, col_indices)),
        shape=(m, n),
        dtype=result_type
    ), nfev


class _Fun_Wrapper:
    # Permits pickling of a wrapped function
    def __init__(self, fun, x0, args, kwargs):
        self.fun = fun
        self.x0 = x0
        self.args = args
        self.kwargs = kwargs
        self._xp = array_namespace(x0)

    def __call__(self, x):
        # send user function same fp type as x0. (but only if cs is not being
        # used
        if self._xp.isdtype(x.dtype, "real floating"):
            x = self._xp.astype(x, self.x0.dtype)
        f = xpx.atleast_nd(
            self.fun(x, *self.args, **self.kwargs),
            ndim=1,
            xp=self._xp
        )
        if f.ndim > 1:
            raise RuntimeError("`fun` return value has "
                               "more than 1 dimension.")
        return f

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_xp"] = state["_xp"].empty(0)
        return state

    def __setstate__(self, state):
        self._xp = array_namespace(state.pop("_xp"))
        self.__dict__.update(state)


def check_derivative(fun, jac, x0, bounds=(-np.inf, np.inf), args=(),
                     kwargs=None):
    """Check correctness of a function computing derivatives (Jacobian or
    gradient) by comparison with a finite difference approximation.

    Parameters
    ----------
    fun : callable
        Function of which to estimate the derivatives. The argument x
        passed to this function is ndarray of shape (n,) (never a scalar
        even if n=1). It must return 1-D array_like of shape (m,) or a scalar.
    jac : callable
        Function which computes Jacobian matrix of `fun`. It must work with
        argument x the same way as `fun`. The return value must be array_like
        or sparse array with an appropriate shape.
    x0 : array_like of shape (n,) or float
        Point at which to estimate the derivatives. Float will be converted
        to 1-D array.
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each bound must match the size of `x0` or be a scalar, in the latter
        case the bound will be the same for all variables. Use it to limit the
        range of function evaluation.
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun` and `jac`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)`` and the same
        for `jac`.

    Returns
    -------
    accuracy : float
        The maximum among all relative errors for elements with absolute values
        higher than 1 and absolute errors for elements with absolute values
        less or equal than 1. If `accuracy` is on the order of 1e-6 or lower,
        then it is likely that your `jac` implementation is correct.

    See Also
    --------
    approx_derivative : Compute finite difference approximation of derivative.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize._numdiff import check_derivative
    >>>
    >>>
    >>> def f(x, c1, c2):
    ...     return np.array([x[0] * np.sin(c1 * x[1]),
    ...                      x[0] * np.cos(c2 * x[1])])
    ...
    >>> def jac(x, c1, c2):
    ...     return np.array([
    ...         [np.sin(c1 * x[1]),  c1 * x[0] * np.cos(c1 * x[1])],
    ...         [np.cos(c2 * x[1]), -c2 * x[0] * np.sin(c2 * x[1])]
    ...     ])
    ...
    >>>
    >>> x0 = np.array([1.0, 0.5 * np.pi])
    >>> check_derivative(f, jac, x0, args=(1, 2))
    2.4492935982947064e-16
    """
    if kwargs is None:
        kwargs = {}
    J_to_test = jac(x0, *args, **kwargs)
    if issparse(J_to_test):
        J_diff = approx_derivative(fun, x0, bounds=bounds, sparsity=J_to_test,
                                   args=args, kwargs=kwargs)
        J_to_test = csr_array(J_to_test)
        abs_err = J_to_test - J_diff
        i, j, abs_err_data = find(abs_err)
        J_diff_data = np.asarray(J_diff[i, j]).ravel()
        return np.max(np.abs(abs_err_data) /
                      np.maximum(1, np.abs(J_diff_data)))
    else:
        J_diff = approx_derivative(fun, x0, bounds=bounds,
                                   args=args, kwargs=kwargs)
        abs_err = np.abs(J_to_test - J_diff)
        return np.max(abs_err / np.maximum(1, np.abs(J_diff)))
