"""Routines for numerical differentiation."""

from __future__ import division

import numpy as np

from ..sparse import issparse, csc_matrix, csr_matrix, coo_matrix, find
from ._group_columns import group_dense, group_sparse

EPS = np.finfo(np.float64).eps


def _adjust_scheme_to_bounds(x0, h, num_steps, scheme, lb, ub):
    """Adjust final difference scheme to the presence of bounds.

    Parameters
    ----------
    x0 : ndarray, shape (n,)
        Point at which we wish to estimate derivative.
    h : ndarray, shape (n,)
        Desired finite difference steps.
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
        Adjusted step sizes. Step size decreases only if a sign flip or
        switching to one-sided scheme doesn't allow to take a full step.
    use_one_sided : ndarray of bool, shape (n,)
        Whether to switch to one-sided scheme. Informative only for
        ``scheme='2-sided'``.
    """
    if scheme == '1-sided':
        use_one_sided = np.ones_like(h, dtype=bool)
    elif scheme == '2-sided':
        h = np.abs(h)
        use_one_sided = np.zeros_like(h, dtype=bool)
    else:
        raise ValueError("`scheme` must be '1-sided' or '2-sided'.")

    if np.all((lb == -np.inf) & (ub == np.inf)):
        return h, use_one_sided

    h_total = h * num_steps
    h_adjusted = h.copy()

    lower_dist = x0 - lb
    upper_dist = ub - x0

    if scheme == '1-sided':
        x = x0 + h_total
        violated = (x < lb) | (x > ub)
        fitting = np.abs(h_total) <= np.maximum(lower_dist, upper_dist)
        h_adjusted[violated & fitting] *= -1

        forward = (upper_dist >= lower_dist) & ~fitting
        h_adjusted[forward] = upper_dist[forward] / num_steps
        backward = (upper_dist < lower_dist) & ~fitting
        h_adjusted[backward] = -lower_dist[backward] / num_steps
    elif scheme == '2-sided':
        central = (lower_dist >= h_total) & (upper_dist >= h_total)

        forward = (upper_dist >= lower_dist) & ~central
        h_adjusted[forward] = np.minimum(
            h[forward], 0.5 * upper_dist[forward] / num_steps)
        use_one_sided[forward] = True

        backward = (upper_dist < lower_dist) & ~central
        h_adjusted[backward] = -np.minimum(
            h[backward], 0.5 * lower_dist[backward] / num_steps)
        use_one_sided[backward] = True

        min_dist = np.minimum(upper_dist, lower_dist) / num_steps
        adjusted_central = (~central & (np.abs(h_adjusted) <= min_dist))
        h_adjusted[adjusted_central] = min_dist[adjusted_central]
        use_one_sided[adjusted_central] = False

    return h_adjusted, use_one_sided


def _compute_absolute_step(rel_step, x0, method):
    if rel_step is None:
        if method == '2-point':
            rel_step = EPS**0.5
        elif method == '3-point':
            rel_step = EPS**(1 / 3)
        elif method == 'cs':
            rel_step = EPS**(0.5)
        else:
            raise ValueError("`method` must be '2-point' or '3-point'.")

    sign_x0 = (x0 >= 0).astype(float) * 2 - 1
    return rel_step * sign_x0 * np.maximum(1.0, np.abs(x0))


def _prepare_bounds(bounds, x0):
    lb, ub = [np.asarray(b, dtype=float) for b in bounds]
    if lb.ndim == 0:
        lb = np.resize(lb, x0.shape)

    if ub.ndim == 0:
        ub = np.resize(ub, x0.shape)

    return lb, ub


def group_columns(A, order=0):
    """Group columns of a 2-d matrix for sparse finite differencing [1]_.

    Two columns are in the same group if in each row at least one of them
    has zero. A greedy sequential algorithm is used to construct groups.

    Parameters
    ----------
    A : array_like or sparse matrix, shape (m, n)
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
        which i-th column assigned. The procedure was helpful only if
        n_groups is significantly less than n.

    References
    ----------
    .. [1] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
           sparse Jacobian matrices", Journal of the Institute of Mathematics
           and its Applications, 13 (1974), pp. 117-120.
    """
    if issparse(A):
        A = csc_matrix(A)
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


def approx_derivative(fun, x0, method='3-point', rel_step=None, f0=None,
                      bounds=(-np.inf, np.inf), vectorized=False,
                      sparsity=None, args=(), kwargs={}):
    """Compute finite difference approximation of the derivatives of a
    vector-valued function.

    If a function maps from R^n to R^m, its derivatives form m-by-n matrix
    called the Jacobian, where an element (i, j) is a partial derivative of
    f[i] with respect to x[j].

    Parameters
    ----------
    fun : callable
        Function of which to estimate the derivatives. It can be implemented
        in two different ways:

            * Non-vectorized version must accept ndarray of shape (n,) and
              return 1-d array_like of shape (m,).
            * Vectorized version must accept ndarray of shape (n, k) and
              return 2-d array_like of shape (m, k). In this way the Jacobian
              can be approximated using a single call of `fun`.

        Use `vectorized` parameter to select an appropriate variant.
    x0 : array_like with shape (n,) or (n, 1) or float
        Point at which to estimate the derivatives.
    method : {'3-point', '2-point'}, optional
        Finite difference method to use:
            - '2-point' - use the fist order accuracy forward or backward
                          difference.
            - '3-point' - use central difference in interior points and the
                          second order accuracy forward or backward difference
                          near the boundary.
            - 'cs' - use a complex-step finite difference scheme. This assumes
                     that the user function is real-valued and can be
                     analytically continued to the complex plane. Otherwise,
                     produces bogus results.
    rel_step : None or array_like, optional
        Relative step size to use. The absolute step size is computed as
        ``h = rel_step * sign(x0) * max(1, abs(x0))``, possibly adjusted to
        fit into the bounds. For ``method='3-point'`` the sign of `h` is
        ignored. If None (default) then step is selected automatically,
        see Notes.
    f0 : None or array_like with shape (m,) or (m, 1), optional
        If not None it is assumed to be equal to ``fun(x0)``, in  this case
        the ``fun(x0)`` is not called. Default is None.
    bounds : tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each bound must match the size of `x0` or be a scalar, in the latter
        case the bound will be the same for all variables. Use it to limit the
        range of function evaluation.
    vectorized : bool, optional
        Whether `fun` is implemented in a vectorized fashion. Default is False.
    sparsity : {None, array_like, sparse matrix, 2-tuple}, optional
        Defines a sparsity structure of the Jacobian matrix. If the Jacobian
        matrix is known to have only few non-zero elements in each row, then
        it's possible to estimate its several columns by a single function
        evaluation [3]_. To perform such economic computations two ingredients
        are required:

        * structure : array_like or sparse matrix of shape (m, n). A zero
          element means that a corresponding element of the Jacobian
          identically equals to zero.
        * groups : array_like of shape (n,). A column grouping for a given
          sparsity structure, use `group_columns` to obtain it.

        A single array or a sparse matrix is interpreted as a sparsity
        structure, and groups are computed inside the function. A tuple is
        interpreted as (structure, groups). If None (default), a standard
        dense differencing will be used.

        Note, that sparse differencing makes sense only for large Jacobian
        matrices where each row contains few non-zero elements.
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)``.

    Returns
    -------
    J : ndarray or csr_matrix
        Finite difference approximation of the Jacobian matrix. If `sparsity`
        is None then ndarray with shape (m, n) is returned. Although if m=1 it
        is returned as a gradient with shape (n,). If `sparsity` is not None,
        csr_matrix with shape (m, n) is returned.

    See Also
    --------
    check_derivative : Check correctness of a function computing derivatives.

    Notes
    -----
    If `rel_step` is not provided, it assigned to ``EPS**(1/s)``, where EPS is
    machine epsilon for float64 numbers, s=2 for '2-point' method and s=3 for
    '3-point' method. Such relative step approximately minimizes a sum of
    truncation and round-off errors, see [1]_.

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
    different cases. b) In all cases np.atleast_2d can be called to get 2-d
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
    >>> from scipy.optimize import approx_derivative
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
    """
    if method not in ['2-point', '3-point', 'cs']:
        raise ValueError("Unknown method '%s'. " % method)

    x0 = np.asarray(x0, dtype=float)
    if x0.ndim > 2 or x0.ndim == 2 and x0.shape[1] > 1:
        raise ValueError("`x0` must be float, 1-dimensional or 2-dimensional "
                         "with 1 column.")
    x0 = x0.ravel()

    lb, ub = _prepare_bounds(bounds, x0)

    if lb.shape != x0.shape or ub.shape != x0.shape:
        raise ValueError("Inconsistent shapes between bounds and `x0`.")

    if vectorized:
        def fun_wrapped(x):
            f = np.atleast_2d(fun(x, *args, **kwargs))
            if f.ndim > 2:
                raise ValueError("`fun` return has more than 2 dimensions.")
            return f
    else:
        def fun_wrapped(x):
            return np.atleast_1d(fun(x, *args, **kwargs))

    if f0 is not None:
        f0 = np.asarray(f0, dtype=float)
        if f0.ndim > 1 or f0.ndim == 2 and f0.shape[1] > 1:
            raise ValueError("`x0` must be float, 1-dimensional or"
                             "2-dimensional with 1 column.")
        f0 = f0.ravel()
    elif not vectorized:
        f0 = fun_wrapped(x0)
        if f0.ndim > 1:
            raise ValueError("`fun` return has more than 1 dimension, but "
                             "vectorized=False.")

    if np.any((x0 < lb) | (x0 > ub)):
        raise ValueError("`x0` violates bound constraints.")

    h = _compute_absolute_step(rel_step, x0, method)

    if method == '2-point':
        h, use_one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '1-sided', lb, ub)
    elif method == '3-point':
        h, use_one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', lb, ub)
    elif method == 'cs':
        use_one_sided = False

    if sparsity is None:
        if vectorized:
            return _dense_difference_vectorized(fun_wrapped, x0, f0, h,
                                                use_one_sided, method)
        else:
            return _dense_difference(fun_wrapped, x0, f0, h, use_one_sided,
                                     method)
    else:
        if not issparse(sparsity) and len(sparsity) == 2:
            structure, groups = sparsity
        else:
            structure = sparsity
            groups = group_columns(sparsity)

        if issparse(structure):
            structure = csc_matrix(structure)
        else:
            structure = np.atleast_2d(structure)

        groups = np.atleast_1d(groups)

        if vectorized:
            return _sparse_difference_vectorized(fun_wrapped, x0, f0, h,
                                                 use_one_sided, structure,
                                                 groups, method)
        else:
            return _sparse_difference(fun_wrapped, x0, f0, h, use_one_sided,
                                      structure, groups, method)


def _dense_difference(fun, x0, f0, h, use_one_sided, method):
    m = f0.size
    n = x0.size
    J_transposed = np.empty((n, m))
    h_vecs = np.diag(h)

    for i in range(h.size):
        if method == '2-point':
            x = x0 + h_vecs[i]
            dx = x[i] - x0[i]  # Recompute dx as exactly representable number.
            df = fun(x) - f0
        elif method == '3-point' and use_one_sided[i]:
            x1 = x0 + h_vecs[i]
            x2 = x0 + 2 * h_vecs[i]
            dx = x2[i] - x0[i]
            f1 = fun(x1)
            f2 = fun(x2)
            df = -3.0 * f0 + 4 * f1 - f2
        elif method == '3-point' and not use_one_sided[i]:
            x1 = x0 - h_vecs[i]
            x2 = x0 + h_vecs[i]
            dx = x2[i] - x1[i]
            f1 = fun(x1)
            f2 = fun(x2)
            df = f2 - f1
        elif method == 'cs':
            f1 = fun(x0 + h_vecs[i]*1.j)
            df = f1.imag
            dx = h_vecs[i, i]

        J_transposed[i] = df / dx

    if m == 1:
        J_transposed = np.ravel(J_transposed)

    return J_transposed.T


def _dense_difference_vectorized(fun, x0, f0, h, use_one_sided, method):
    h_vecs = np.diag(h)
    x0 = x0[:, None]
    if method == '2-point':
        X = x0 + h_vecs
        dx = np.diag(X) - x0.ravel()
        if f0 is None:
            X = np.hstack((X, x0))
            F = fun(X)
            f0 = F[:, -1]
            F = F[:, :-1]
        else:
            F = fun(X)
        F -= f0[:, None]
        F /= dx
    elif method == '3-point':
        n = x0.size
        one_sided = np.nonzero(use_one_sided)[0]
        two_sided = np.nonzero(~use_one_sided)[0]
        X = np.empty((n, 2 * n))
        X1 = X[:, :n]
        X2 = X[:, n:]
        X1[:] = x0 + h_vecs
        X2[:, one_sided] = x0 + 2 * h_vecs[:, one_sided]
        X2[:, two_sided] = x0 - h_vecs[:, two_sided]
        dx = np.empty_like(h)
        dx[one_sided] = X2[one_sided, one_sided] - x0[one_sided].ravel()
        dx[two_sided] = X1[two_sided, two_sided] - X2[two_sided, two_sided]

        if f0 is None:
            X = np.hstack((X, x0))
            F = fun(X)
            f0 = F[:, -1]
            F = F[:, :-1]
        else:
            F = fun(X)

        F1 = F[:, :n]
        F2 = F[:, n:]
        F1[:, one_sided] *= 4
        F1[:, one_sided] -= F2[:, one_sided]
        F1[:, one_sided] -= 3 * f0[:, None]

        F1[:, two_sided] -= F2[:, two_sided]
        F1 /= dx
        F = F1
    elif method == 'cs':
        X = x0 + h_vecs * 1j
        F = fun(X)
        F = F.imag / h

    if F.shape[0] == 1:
        F = np.ravel(F)

    return F


def _sparse_difference(fun, x0, f0, h, use_one_sided,
                       structure, groups, method):
    m = f0.size
    n = x0.size
    row_indices = []
    col_indices = []
    fractions = []

    n_groups = np.max(groups) + 1
    for group in range(n_groups):
        # Perturb variables which are in the same group simultaneously.
        e = np.equal(group, groups)
        h_vec = h * e
        if method == '2-point':
            x = x0 + h_vec
            dx = x - x0
            df = fun(x) - f0
            # The result is  written to columns which correspond to perturbed
            # variables.
            cols, = np.nonzero(e)
            # Find all non-zero elements in selected columns of Jacobian.
            i, j, _ = find(structure[:, cols])
            # Restore column indices in the full array.
            j = cols[j]
        elif method == '3-point':
            # Here we do conceptually the same but separate one-sided
            # and two-sided schemes.
            x1 = x0.copy()
            x2 = x0.copy()

            mask_1 = use_one_sided & e
            x1[mask_1] += h_vec[mask_1]
            x2[mask_1] += 2 * h_vec[mask_1]

            mask_2 = ~use_one_sided & e
            x1[mask_2] -= h_vec[mask_2]
            x2[mask_2] += h_vec[mask_2]

            dx = np.zeros(n)
            dx[mask_1] = x2[mask_1] - x0[mask_1]
            dx[mask_2] = x2[mask_2] - x1[mask_2]

            f1 = fun(x1)
            f2 = fun(x2)

            cols, = np.nonzero(e)
            i, j, _ = find(structure[:, cols])
            j = cols[j]

            mask = use_one_sided[j]
            df = np.empty(m)

            rows = i[mask]
            df[rows] = -3 * f0[rows] + 4 * f1[rows] - f2[rows]

            rows = i[~mask]
            df[rows] = f2[rows] - f1[rows]
        elif method == 'cs':
            f1 = fun(x0 + h_vec*1.j)
            df = f1.imag
            dx = h_vec
            cols, = np.nonzero(e)
            i, j, _ = find(structure[:, cols])
            j = cols[j]
        else:
            raise ValueError("Never be here.")

        # All that's left is to compute the fraction. We store i, j and
        # fractions as separate arrays and later construct coo_matrix.
        row_indices.append(i)
        col_indices.append(j)
        fractions.append(df[i] / dx[j])

    row_indices = np.hstack(row_indices)
    col_indices = np.hstack(col_indices)
    fractions = np.hstack(fractions)
    J = coo_matrix((fractions, (row_indices, col_indices)), shape=(m, n))
    return csr_matrix(J)


def _sparse_difference_vectorized(fun, x0, f0, h, use_one_sided, structure,
                                  groups, method):
    n = x0.shape[0]
    n_groups = np.max(groups) + 1
    h_vecs = np.empty((n_groups, n))
    cols_in_group = []
    for group in range(n_groups):
        e = np.equal(group, groups)
        cols_in_group.append(np.nonzero(e)[0])
        h_vecs[group] = h * e

    h_vecs = h_vecs.T
    x0 = x0[:, None]

    row_indices = []
    col_indices = []
    fractions = []

    if method == '2-point':
        X = x0 + h_vecs
        if f0 is None:
            X = np.hstack((X, x0))
            F = fun(X)
            X = X[:, :-1]
            f0 = F[:, -1]
            F = F[:, :-1]
        else:
            F = fun(X)

        m = f0.shape[0]

        F -= f0[:, None]
        X -= x0

        for cols, dx, df in zip(cols_in_group, X.T, F.T):
            i, j, _ = find(structure[:, cols])
            j = cols[j]
            row_indices.append(i)
            col_indices.append(j)
            fractions.append(df[i] / dx[j])

    # This case often occurs and can be implemented more efficiently.
    elif method == '3-point' and not np.any(use_one_sided):
        X = np.empty((n, 2 * n_groups))
        X1 = X[:, :n_groups]
        X2 = X[:, n_groups:]
        X1[:] = x0 - h_vecs
        X2[:] = x0 + h_vecs

        F = fun(X)
        m = F.shape[0]

        F1 = F[:, :n_groups]
        F2 = F[:, n_groups:]
        F2 -= F1

        X2 -= X1

        for cols, dx, df in zip(cols_in_group, X2.T, F2.T):
            i, j, _ = find(structure[:, cols])
            j = cols[j]
            row_indices.append(i)
            col_indices.append(j)
            fractions.append(df[i] / dx[j])

    elif method == '3-point':
        one_sided = np.nonzero(use_one_sided)[0]
        two_sided = np.nonzero(~use_one_sided)[0]
        X = np.empty((n, 2 * n_groups))
        X1 = X[:, :n_groups]
        X2 = X[:, n_groups:]
        X1[:] = x0 + h_vecs
        X2[one_sided] = x0[one_sided] + 2 * h_vecs[one_sided]
        X2[two_sided] = x0[two_sided] - h_vecs[two_sided]

        if f0 is None:
            X = np.hstack((X, x0))
            F = fun(X)
            f0 = F[:, -1]
            F = F[:, :-1]
        else:
            F = fun(X)

        m = f0.shape[0]

        F1 = F[:, :n_groups]
        F2 = F[:, n_groups:]

        X2[one_sided] -= x0[one_sided]
        X2[two_sided] *= -1
        X2[two_sided] += X1[two_sided]

        df = np.empty(m)
        for group, (cols, dx) in enumerate(zip(cols_in_group, X2.T)):
            i, j, _ = find(structure[:, cols])
            j = cols[j]

            mask = use_one_sided[j]

            rows = i[mask]
            df[rows] = 4 * F1[rows, group] - 3 * f0[rows] - F2[rows, group]

            rows = i[~mask]
            df[rows] = F1[rows, group] - F2[rows, group]

            row_indices.append(i)
            col_indices.append(j)
            fractions.append(df[i] / dx[j])

    elif method == 'cs':
        X = x0 + h_vecs * 1j
        F = fun(X)
        m = F.shape[0]
        df_all = F.imag.T
        dx_all = h_vecs.T

        for cols, dx, df in zip(cols_in_group, dx_all, df_all):
            i, j, _ = find(structure[:, cols])
            j = cols[j]
            row_indices.append(i)
            col_indices.append(j)
            fractions.append(df[i] / dx[j])

    row_indices = np.hstack(row_indices)
    col_indices = np.hstack(col_indices)
    fractions = np.hstack(fractions)
    J = coo_matrix((fractions, (row_indices, col_indices)), shape=(m, n))
    return csr_matrix(J)


def check_derivative(fun, jac, x0, bounds=(-np.inf, np.inf), vectorized=False,
                     args=(), kwargs={}):
    """Check correctness of a function computing derivatives (Jacobian or
    gradient) by comparison with a finite difference approximation.

    Parameters
    ----------
    fun : callable
        Function of which to estimate the derivatives. It can be implemented
        in two different ways:

            * Non-vectorized version must accept ndarray of shape (n,) and
              return 1-d array_like of shape (m,).
            * Vectorized version must accept ndarray of shape (n, k) and
              return 2-d array_like of shape (m, k). In this way the Jacobian
              can be approximated using a single call of `fun`.

        Use `vectorized` parameter to select an appropriate variant.
    jac : callable
        Function which computes Jacobian matrix of `fun`. Note that a raveled
        version of `x` will be passed to `jac`The return value must be
        array_like or sparse matrix with an appropriate shape.
    x0 : array_like with shape (n,) or (n, 1) or float
        Point at which to estimate the derivatives.
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each bound must match the size of `x0` or be a scalar, in the latter
        case the bound will be the same for all variables. Use it to limit the
        range of function evaluation.
    vectorized : bool, optional
        Whether `fun` is implemented in a vectorized fashion. Default is False.
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
    >>> from scipy.optimize import check_derivative
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
    x0 = np.asarray(x0, dtype=float)
    if x0.ndim > 2 or x0.ndim == 2 and x0.shape[1] > 1:
        raise ValueError("`x0` must be float, 1-dimensional or 2-dimensional "
                         "with 1 column.")
    if x0.ndim == 2:  # Keep float as float.
        x0 = x0.ravel()

    J_to_test = jac(x0, *args, **kwargs)
    if issparse(J_to_test):
        J_diff = approx_derivative(fun, x0, bounds=bounds, sparsity=J_to_test,
                                   vectorized=vectorized, args=args,
                                   kwargs=kwargs)
        J_to_test = csr_matrix(J_to_test)
        abs_err = J_to_test - J_diff
        i, j, abs_err_data = find(abs_err)
        J_diff_data = np.asarray(J_diff[i, j]).ravel()
        return np.max(np.abs(abs_err_data) /
                      np.maximum(1, np.abs(J_diff_data)))
    else:
        J_diff = approx_derivative(fun, x0, bounds=bounds,
                                   vectorized=vectorized, args=args,
                                   kwargs=kwargs)
        abs_err = np.abs(J_to_test - J_diff)
        return np.max(abs_err / np.maximum(1, np.abs(J_diff)))
