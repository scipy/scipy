"""Boundary value problem solver."""
from __future__ import division, print_function, absolute_import

from warnings import warn

import numpy as np
from numpy.linalg import norm, pinv

from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import splu
from scipy.optimize import OptimizeResult


EPS = np.finfo(float).eps


def compute_fun_jac(fun, x, y, p, f0=None):
    """Compute derivatives of and ODE system rhs with forward differences.

    Returns
    -------
    df_dy : ndarray, shape (n, n, m)
        Derivatives with respect to y. An element (i, j, k) corresponds to
        d f_i(x_k, y_k) / d y_j.
    df_dp : ndarray, shape (n, k, m) or (0,)
        Derivatives with respect to p. An element (i, j, k) corresponds to
        d f_i(x_k, y_k, p) / d p_j. If `p` is empty, None is returned.
    """
    n, m = y.shape
    if f0 is None:
        f0 = fun(x, y, p)

    df_dy = np.empty((n, n, m))
    h = EPS ** 0.5 * (1 + np.abs(y))
    for i in range(n):
        y_new = y.copy()
        y_new[i] += h[i]
        hi = y_new[i] - y[i]
        f_new = fun(x, y_new, p)
        df_dy[:, i, :] = (f_new - f0) / hi

    k = p.shape[0]
    if k == 0:
        df_dp = None
    else:
        df_dp = np.empty((n, k, m))
        h = EPS ** 0.5 * (1 + np.abs(p))
        for i in range(k):
            p_new = p.copy()
            p_new[i] += h[i]
            hi = p_new[i] - p[i]
            f_new = fun(x, y, p_new)
            df_dp[:, i, :] = (f_new - f0) / hi

    return df_dy, df_dp


def compute_bc_jac(bc, ya, yb, p, bc0=None):
    """Compute derivatives of boundary conditions with forward differences.

    Returns
    -------
    dbc_dya : ndarray, shape (n + k, n)
        Derivatives with respect to ya. An element (i, j) corresponds to
        d bc_i / d ya_j.
    dbc_dyb : ndarray, shape (n + k, n)
        Derivatives with respect to yb. An element (i, j) corresponds to
        d bc_i / d ya_j.
    dbc_dp : ndarray, shape (n + k, k)
        Derivatives with respect to p. An element (i, j) corresponds to
        d bc_i / d p_j. If `p` is empty, None is returned.
    """
    n = ya.shape[0]
    k = p.shape[0]

    if bc0 is None:
        bc0 = bc(ya, yb, p)

    dbc_dya = np.empty((n, n + k))
    h = EPS**0.5 * (1 + np.abs(ya))
    for i in range(n):
        ya_new = ya.copy()
        ya_new[i] += h[i]
        hi = ya_new[i] - ya[i]
        bc_new = bc(ya_new, yb, p)
        dbc_dya[i] = (bc_new - bc0) / hi
    dbc_dya = dbc_dya.T

    h = EPS**0.5 * (1 + np.abs(yb))
    dbc_dyb = np.empty((n, n + k))
    for i in range(n):
        yb_new = yb.copy()
        yb_new[i] += h[i]
        hi = yb_new[i] - yb[i]
        bc_new = bc(ya, yb_new, p)
        dbc_dyb[i] = (bc_new - bc0) / hi
    dbc_dyb = dbc_dyb.T

    if k == 0:
        dbc_dp = None
    else:
        h = EPS ** 0.5 * (1 + np.abs(p))
        dbc_dp = np.empty((k, n + k))
        for i in range(k):
            p_new = p.copy()
            p_new[i] += h[i]
            hi = p_new[i] - p[i]
            bc_new = bc(ya, yb, p_new)
            dbc_dp[i] = (bc_new - bc0) / hi
        dbc_dp = dbc_dp.T

    return dbc_dya, dbc_dyb, dbc_dp


def compute_jac_indices(n, m, k):
    """Compute indices for the global Jacobian construction."""
    i_col = np.repeat(np.arange((m - 1) * n), n)
    j_col = (np.tile(np.arange(n), n * (m - 1)) +
             np.repeat(np.arange(m - 1) * n, n**2))

    i_bc = np.repeat(np.arange((m - 1) * n, m * n + k), n)
    j_bc = np.tile(np.arange(n), n + k)

    i_p_col = np.repeat(np.arange((m - 1) * n), k)
    j_p_col = np.tile(np.arange(m * n, m * n + k), (m - 1) * n)

    i_p_bc = np.repeat(np.arange((m - 1) * n, m * n + k), k)
    j_p_bc = np.tile(np.arange(m * n, m * n + k), n + k)

    i = np.hstack((i_col, i_col, i_bc, i_bc, i_p_col, i_p_bc))
    j = np.hstack((j_col, j_col + n,
                   j_bc, j_bc + (m - 1) * n,
                   j_p_col, j_p_bc))

    return i, j


def construct_global_jac(n, m, k, i_jac, j_jac, h, df_dy, df_dy_middle, df_dp,
                         df_dp_middle, dbc_dya, dbc_dyb, dbc_dp):
    """Construct the Jacobian of the whole collocation system.

    Returns
    -------
    J : csc_matrix, shape (n * m + k, n * m + k)
        Global Jacobian in a sparse form.
    """
    df_dy = np.transpose(df_dy, (2, 0, 1))
    df_dy_middle = np.transpose(df_dy_middle, (2, 0, 1))

    h = h[:, np.newaxis, np.newaxis]

    dPhi_dy_0 = np.empty((m - 1, n, n))
    dPhi_dy_0[:] = -np.identity(n)
    dPhi_dy_0 -= h / 6 * (df_dy[:-1] + 2 * df_dy_middle)
    T = np.einsum('...ij,...jk->...ik', df_dy_middle, df_dy[:-1])
    dPhi_dy_0 -= h**2 / 12 * T

    dPhi_dy_1 = np.empty((m - 1, n, n))
    dPhi_dy_1[:] = np.identity(n)
    dPhi_dy_1 -= h / 6 * (df_dy[1:] + 2 * df_dy_middle)
    T = np.einsum('...ij,...jk->...ik', df_dy_middle, df_dy[1:])
    dPhi_dy_1 += h**2 / 12 * T

    values = np.hstack((dPhi_dy_0.ravel(), dPhi_dy_1.ravel(), dbc_dya.ravel(),
                        dbc_dyb.ravel()))

    if k > 0:
        df_dp = np.transpose(df_dp, (2, 0, 1))
        df_dp_middle = np.transpose(df_dp_middle, (2, 0, 1))
        T = np.einsum('...ij,...jk->...ik', df_dy_middle,
                      df_dp[:-1] - df_dp[1:])
        df_dp_middle += 0.125 * h * T
        dPhi_dp = -h/6 * (df_dp[:-1] + df_dp[1:] + 4 * df_dp_middle)
        values = np.hstack((values, dPhi_dp.ravel(), dbc_dp.ravel()))

    J = coo_matrix((values, (i_jac, j_jac)))
    return csc_matrix(J)


def collocation_fun(fun, y, p, x, h):
    """Evaluate collocation residuals.

    This function lies in the core of the method. The solution is sought
    as a cubic C^1 continuous spline with derivatives matching the ODE rhs
    at given nodes `x`. Collocation conditions are formed from the equality
    of our solution derivatives and rhs of the ODE system in the middle points
    between nodes.

    Such method is classified to Lobbato IIIA family in ODE literature.
    Refer to [1]_ for the formula and some discussion.

    Returns
    -------
    col_res : ndarray, shape (n, m - 1)
        Collocation residuals at the middle points of mesh intervals.
    y_middle : ndarray, shape (n, m - 1)
        Values of the cubic interpolant evaluated at the middle points of
        mesh intervals.
    f : ndarray, shape (n, m)
        RHS of the ODE system evaluated at mesh nodes.
    f_middle : ndarray, shape (n, m - 1)
        RHS of the ODE system evaluated at `y_middle.`

    References
    ----------
    .. [1] J. Kierzenka, L. F. Shampine, "A BVP Solver Based on Residual
           Control and the Maltab PSE", ACM Trans. Math. Softw., Vol. 27,
           Number 3, pp. 299-316, 2001.
    """
    f = fun(x, y, p)
    y_middle = (0.5 * (y[:, 1:] + y[:, :-1]) -
                0.125 * h * (f[:, 1:] - f[:, :-1]))
    f_middle = fun(x[:-1] + 0.5 * h, y_middle, p)
    col_res = y[:, 1:] - y[:, :-1] - h / 6 * (f[:, :-1] + f[:, 1:] +
                                              4 * f_middle)

    return col_res, y_middle, f, f_middle


def prepare_sys(n, m, k, fun, bc, x, h):
    """Create the function and the Jacobian for the collocation system."""
    x_middle = x[:-1] + 0.5 * h
    i_jac, j_jac = compute_jac_indices(n, m, k)

    def sys_fun(y, p):
        return collocation_fun(fun, y, p, x, h)

    def sys_jac(y, p, y_middle, f, f_middle, bc0):
        df_dy, df_dp = compute_fun_jac(fun, x, y, p, f)
        df_dy_middle, df_dp_middle = compute_fun_jac(fun, x_middle,
                                                     y_middle, p, f_middle)
        dbc_dya, dbc_dyb, dbc_dp = compute_bc_jac(bc, y[:, 0], y[:, -1],
                                                  p, bc0)
        return construct_global_jac(n, m, k, i_jac, j_jac, h, df_dy,
                                    df_dy_middle, df_dp, df_dp_middle, dbc_dya,
                                    dbc_dyb, dbc_dp)

    return sys_fun, sys_jac


def solve_newton(n, m, k, col_fun, bc, jac, y, p, B, tol, n_iter):
    """Solve the nonlinear collocation system by a Newton method.

    This is a simple Newton method with a backtracking line search. As
    advised in [1]_, an affine-invariant criterion function F = ||J^-1 r||^2
    is used, where J is the Jacobian matrix at the current iteration and r is
    the vector or residuals.

    There are other tricks proposed in [1]_, but they are not used as they
    don't seem to improve anything significantly, and even break the
    convergence on some test problems I tried.

    Parameters
    ----------
    n : int
        Number of equations in the ODE.
    m : int
        Number of nodes in the mesh.
    k : int
        Number of unknown parameters.
    col_fun : callable
        Function computing collocation residuals.
    bc : callable
        Function computing boundary condition residuals.
    jac : callable
        Function computing the Jacobian of the whole system (including
        collocation and boundary condition residuals). It is supposed to
        return csc_matrix.
    y : ndarray, shape (n, m)
        Initial guess for the function values at the mesh nodes.
    p : ndarray, shape (k,)
        Initial guess for the unknown parameters.
    B : ndarray with shape (n, n) or None
        Matrix to force the S y(a) = 0 condition for a problems with the
        singular term. If None, the singular term is assumed to be absent.
    tol : float
        Iterations are terminated if the uniform norm of the system residuals
        less than `tol`. The reasonable value of `tol` should be 1-2 orders
        less than the tolerance for our BVP solver.
    n_iter : int
        Max number of iterations. A small number should be used when solving a
        boundary value problem, because it's better to refine the mesh and
        start over than try hard on the current mesh.

    Returns
    -------
    y : ndarray, shape (n, m)
        Final iterate for the function values at the mesh nodes.
    p : ndarray, shape (k,)
        Final iterate for the unknown parameters.
    singular : bool
        True, if the LU decomposition failed because Jacobian turned out
        to be singular.

    References
    ----------
    .. [1]  U. Ascher, R. Mattheij and R. Russell "Numerical Solution of
       Boundary Value Problems for Ordinary Differential Equations"
    """
    res_col_new, y_middle, f, f_middle = col_fun(y, p)
    res_bc = bc(y[:, 0], y[:, -1], p)
    res = np.hstack((res_col_new.ravel(order='F'), res_bc))

    sigma = 0.05  # Minimum relative improvement (Armijo constant).
    tau = 0.5  # Step size decrease factor for backtracking.

    singular = False
    for iteration in range(n_iter):
        J = jac(y, p, y_middle, f, f_middle, res_bc)
        try:
            LU = splu(J)
        except RuntimeError:
            singular = True
            break

        step = -LU.solve(res)
        y_step = step[:m * n].reshape((n, m), order='F')
        p_step = step[m * n:]

        cost = np.dot(step, step)
        alpha = 1
        n_trial = 5
        for i in range(n_trial):
            y_new = y + alpha * y_step
            if B is not None:
                y_new[:, 0] = np.dot(B, y_new[:, 0])
            p_new = p + alpha * p_step

            res_col_new, y_middle, f, f_middle = col_fun(y_new, p_new)
            res_bc = bc(y_new[:, 0], y_new[:, -1], p_new)
            res_new = np.hstack((res_col_new.ravel(order='F'), res_bc))

            s = LU.solve(res_new)
            cost_new = np.dot(s, s)
            if cost_new <= (1 - 2 * alpha * sigma) * cost:
                break

            alpha *= tau

        y = y_new
        p = p_new
        res = res_new

        if norm(res, ord=np.inf) < tol:
            break

    return y, p, singular


def print_iteration_header():
    print("{:^15}{:^15}{:^15}{:^15}".format(
        "Iteration", "Max residual", "Total nodes", "Nodes added"))


def print_iteration_progress(iteration, residual, total_nodes, inserted_nodes):
    print("{:^15}{:^15.2e}{:^15}{:^15}".format(
        iteration, residual, total_nodes, inserted_nodes))


class BVPResult(OptimizeResult):
    pass


TERMINATION_MESSAGES = {
    0: "The algorithm converged to the desired accuracy.",
    1: "The maximum number of mesh nodes is exceeded.",
    2: "A singular Jacobian encountered when solving the collocation system."
}


def estimate_rms_residuals(fun, sol, x, h, p, res_middle, f_middle):
    """Estimate rms values of collocation residuals using Lobatto quadrature.

    The residuals are defined as the difference between the derivatives of
    our solution and the rhs of the ODE system. We use relative residuals, i.e.
    normalized by 1 + np.abs(f). RMS values are computed as sqrt
    from the integrals of the squared relative residuals over each interval.
    Integrals are estimated using 5-point Lobatto quadrature [1]_, we use the
    fact that residuals at the mesh nodes are identically zero.

    In [2] they don't normalize integrals by interval lengths, which gives
    a higher rate of convergence of the residuals by the factor of h**0.5.
    I chose to do such normalization for an ease of interpretation of return
    values as RMS estimates.

    Returns
    -------
    rms_res : ndarray, shape (m - 1,)
        Estimated rms values of the relative residuals over each interval.

    References
    ----------
    .. [1] http://mathworld.wolfram.com/LobattoQuadrature.html
    .. [2] J. Kierzenka, L. F. Shampine, "A BVP Solver Based on Residual
       Control and the Maltab PSE", ACM Trans. Math. Softw., Vol. 27,
       Number 3, pp. 299-316, 2001.
    """
    x_middle = x[:-1] + 0.5 * h
    s = 0.5 * h * (3/7)**0.5
    x1 = x_middle + s
    x2 = x_middle - s
    y1 = sol(x1)
    y2 = sol(x2)
    y1_prime = sol(x1, 1)
    y2_prime = sol(x2, 1)
    f1 = fun(x1, y1, p)
    f2 = fun(x2, y2, p)
    res_1 = y1_prime - f1
    res_2 = y2_prime - f2

    res_middle /= 1 + np.abs(f_middle)
    res_1 /= 1 + np.abs(f1)
    res_2 /= 1 + np.abs(f2)

    res_1 = np.sum(res_1**2, axis=0)
    res_2 = np.sum(res_2**2, axis=0)
    res_middle = np.sum(res_middle**2, axis=0)

    return (0.5 * (32/45 * res_middle + 49/90 * (res_1 + res_2))) ** 0.5


def create_spline(y, f, x, h):
    """Create a cubic spline given values and derivatives.

    Formulas for the coefficients are taken from interpolate.CubicSpline.

    Returns
    -------
    sol : PPoly
        Constructed spline as a PPoly instance.
    """
    from scipy.interpolate import PPoly

    n, m = y.shape
    c = np.empty((4, n, m - 1))
    slope = (y[:, 1:] - y[:, :-1]) / h
    t = (f[:, :-1] + f[:, 1:] - 2 * slope) / h
    c[0] = t / h
    c[1] = (slope - f[:, :-1]) / h - t
    c[2] = f[:, :-1]
    c[3] = y[:, :-1]
    c = np.rollaxis(c, 1)

    return PPoly(c, x, extrapolate=True, axis=1)


def modify_mesh(x, insert_1, insert_2):
    """Insert nodes to a mesh or remove nodes from a mesh.

    Currently only inserting is implemented. The need of knots deletion is
    questionable.

    Parameters
    ----------
    x : ndarray, shape (m,)
        Mesh nodes.
    insert_1 : ndarray
        Intervals to each insert 1 new node in the middle.
    insert_2 : ndarray
        Intervals to each insert 2 new nodes, such that divide an interval
        into 3 equal parts.

    Returns
    -------
    x_new : ndarray
        New mesh nodes.

    Notes
    -----
    `insert_1` and `insert_2` should not have common values.
    """
    nodes_1 = 0.5 * (x[insert_1] + x[insert_1 + 1])
    nodes_2 = np.hstack((2 * x[insert_2] + x[insert_2 + 1],
                         x[insert_2] + 2 * x[insert_2 + 1])) / 3

    ind = np.hstack((insert_1 + 1, np.tile(insert_2, 2) + 1))
    nodes = np.hstack((nodes_1, nodes_2))

    return np.insert(x, ind, nodes)


def solve_bvp(fun, bc, x, y, p=None, S=None, tol=1e-3, max_nodes=1000,
              verbose=0):
    """Solve a boundary-value problem for a system of ODEs.

    This function numerically solves a first order system of ODEs subject to
    two-point boundary conditions::

        dy / dx = f(x, y, p) + S * y / (x - a), a <= x <= b
        bc(y(a), y(b), p) = 0

    Here x is a 1-dimensional independent variable, y(x) is a n-dimensional
    vector-valued function and p is a k-dimensional vector of unknown
    parameters which is to be found along with y(x). For the problem to be
    determined there must be n + k boundary conditions, i.e. bc must be
    (n + k)-dimensional function.

    The last singular term in the right-hand side of the system is optional.
    It is defined by a n-by-n matrix S, such that the solution must satisfy
    S y(a) = 0. This condition will be forced during iterations, so it must not
    contradict boundary conditions. See [2] for the explanation how this term
    is handled when solving BVPs numerically.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is
        ``fun(x, y, p)``, where ``x`` is ndarray of shape (m,) and ``y`` is
        ndarray of shape (n, m), meaning that ``y[:, i]`` corresponds to
        ``x[i]``, ``p`` is optional with shape (k,). The return value must
        have shape (n, m) with the same layout as for ``y``.
    bc : callable
        Function evaluating residuals of the boundary conditions. The calling
        signature is ``bc(ya, yb, p)``, where ``ya`` and ``yb`` have shape
        (n,), ``p`` is optional with shape (k,). The return value must have
        shape (n + k,).
    x : array_like, shape (m,)
        Initial mesh. Must be a strictly increasing sequence with ``x[0]=a``
        and ``x[-1]=b``. It is recommended that m * n + k be less than
        `max_dense_size` such that the initial system can be solved with a
        robust dense solver.
    y : array_like, shape (n, m)
        Initial guess for the function values at the mesh nodes, i-th column
        corresponds to ``x[i]``.
    p : array_like with shape (k,) or None, optional
        Initial guess for the unknown parameters. If None (default), it is
        assumed that the problem doesn't depend on any parameters.
    S : array_like with shape (n, n) or None
        Matrix defining the singular term. If None (default), the problem is
        solved without the singular term.
    tol : float, optional
        Desired tolerance of the solution. If we define ``res = y' - f`` where
        ``y'`` is the derivative of the found solution, then the solver tries
        to achieve on each interval ``norm(res / (1 + abs(f)) < tol``, where
        ``norm`` is estimated in a root mean squared sense (using a numerical
        quadrature formula). Default is 1e-3.
    max_nodes : int, optional
        Maximum allowed number of the mesh nodes. If exceeded, the algorithm
        terminates. Default is 1000.
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

            * 0 (default) : work silently.
            * 1 : display a termination report.
            * 2 : display progress during iterations.

    Returns
    -------
    Bunch object with the following fields defined:
    sol : PPoly
        Found solution for y as `scipy.interpolate.PPoly` instance, a C^1
        continuous cubic spline.
    p : ndarray or None, shape (k,)
        Found parameters. None, if the parameters were not present in the
        problem.
    x : ndarray, shape (m,)
        Nodes of the final mesh.
    y : ndarray, shape (n, m)
        Function values evaluated at the mesh nodes.
    f : ndarray, shape (n, m)
        Function derivatives evaluated at the mesh nodes.
    res : ndarray, shape (m - 1,)
        Relative rms residuals for each mesh interval.
    niter : int
        Number of completed iterations.
    status : int
        The reason for algorithm termination:

            * 0: The algorithm converged to the desired accuracy.
            * 1: The maximum number of mesh nodes is exceeded.
            * 2: A singular Jacobian encountered when solving the collocation
                 system.

    message : string
        Verbal description of the termination reason.
    success : bool
        True if the algorithm converged to the desired accuracy (``status=0``).

    Notes
    -----
    This function implements a 4-th order collocation algorithm with the
    control of residuals similar to [1]_. A collocation system is solved
    by a damped Newton method with an affine-invariant criterion function as
    described in [3]_.

    Note that in [1]_  integral residuals are defined without normalization
    by interval lengths. So their definition is different by a multiplier of
    h**0.5 (h is an interval length) from the definition used here.

    .. versionadded:: 0.18.0

    References
    ----------
    .. [1] J. Kierzenka, L. F. Shampine, "A BVP Solver Based on Residual
       Control and the Maltab PSE", ACM Trans. Math. Softw., Vol. 27,
       Number 3, pp. 299-316, 2001.
    .. [2] L.F. Shampine, P. H. Muir and H. Xu, "A User-Friendly Fortran BVP
           Solver".
    .. [3]  U. Ascher, R. Mattheij and R. Russell "Numerical Solution of
            Boundary Value Problems for Ordinary Differential Equations"

    Examples
    --------
    We solve a simple Sturm-Liouville problem::

        y'' + k**2 * y = 0
        y(0) = y(1) = 0

    It is known that a non-trivial solution y = A * sin(k * x) is possible for
    k = pi * n, where n is an integer. To establish the normalization constant
    A = 1 we add one more boundary condition::

        y'(0) = k

    Now we rewrite our equation as a first order system and implement its
    right-hand side computation::

        y1' = y2
        y2' = -k**2 * y1

    >>> def fun(x, y, p):
    ...     k = p[0]
    ...     return np.vstack((y[1], -k**2 * y[0]))

    Note that parameters p are passed as a vector (with one element in our
    case).

    Implement the boundary conditions:

    >>> def bc(ya, yb, p):
    ...     k = p[0]
    ...     return np.array([ya[0], yb[0], ya[1] - k])

    We initially setup a coarse mesh with only 5 nodes. To avoid the trivial
    solution we set all values of y to 1 initially.

    >>> x = np.linspace(0, 1, 5)
    >>> y = np.ones((2, x.shape[0]))

    Now we are ready to run the solver. We choose 3 as an initial guess for
    k, aiming to find a solution for k = pi.

    >>> from scipy.integrate import solve_bvp
    >>> sol = solve_bvp(fun, bc, x, y, p=[3])

    We see that we found an approximately correct value for k:

    >>> sol.p[0]
    3.14169442647

    And finally we plot the solution. We take an advantage of having the
    solution in a spline form to produce a smooth plot.

    >>> x_plot = np.linspace(0, 1, 100)
    >>> y_plot = sol.sol(x_plot)[0]
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x_plot, y_plot)
    >>> plt.xlabel("x")
    >>> plt.ylabel("y")
    >>> plt.show()
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("`x` must be 1 dimensional.")
    h = np.diff(x)
    if np.any(h <= 0):
        raise ValueError("`x` must be strictly increasing.")
    a = x[0]

    y = np.asarray(y, dtype=float)
    if y.ndim != 2:
        raise ValueError("`y` must be 2 dimensional.")
    if y.shape[1] != x.shape[0]:
        raise ValueError("`y` is expected to have {} columns, but actually "
                         "has {}.".format(x.shape[0], y.shape[1]))

    if p is None:
        p = np.array([])
    else:
        p = np.asarray(p, dtype=float)
    if p.ndim != 1:
        raise ValueError("`p` must be 1 dimensional.")

    if tol < 100 * EPS:
        warn("`tol` is too low, setting to {:.2e}".format(100 * EPS))
        tol = 100 * EPS

    n = y.shape[0]
    k = p.shape[0]

    if S is not None:
        S = np.asarray(S)
        if S.shape != (n, n):
            raise ValueError("`S` is expected to have shape {}, "
                             "but actually has {}".format((n, n), S.shape))

        # Compute I - S^+ S to impose necessary boundary conditions.
        B = np.identity(n) - np.dot(pinv(S), S)

        y[:, 0] = np.dot(B, y[:, 0])

        # Compute (I - S)^+ to correct derivatives at x=a.
        D = pinv(np.identity(n) - S)
    else:
        B = None
        D = None

    if k == 0:
        if S is None:
            def fun_wrapped(x, y, _):
                return np.asarray(fun(x, y))
        else:
            def fun_wrapped(x, y, _):
                f = np.asarray(fun(x, y))
                if x[0] == a:
                    f[:, 0] = np.dot(D, f[:, 0])
                    f[:, 1:] += np.dot(S, y[:, 1:]) / (x[1:] - a)
                else:
                    f += np.dot(S, y) / (x - a)
                return f

        def bc_wrapped(ya, yb, _):
            return np.asarray(bc(ya, yb))
    else:
        if S is None:
            def fun_wrapped(x, y, p):
                return np.asarray(fun(x, y, p))
        else:
            def fun_wrapped(x, y, p):
                f = np.asarray(fun(x, y, p))
                if x[0] == a:
                    f[:, 0] = np.dot(D, f[:, 0])
                    f[:, 1:] += np.dot(S, y[:, 1:]) / (x[1:] - a)
                else:
                    f += np.dot(S, y) / (x - a)
                return f

        def bc_wrapped(ya, yb, p):
            return np.asarray(bc(ya, yb, p))

    f = fun_wrapped(x, y, p)
    if f.shape != y.shape:
        raise ValueError("`fun` return is expected to have shape {}, "
                         "but actually has {}.".format(y.shape, f.shape))

    bc_res = bc_wrapped(y[:, 0], y[:, -1], p)
    if bc_res.shape != (n + k,):
        raise ValueError("`bc` return is expected to have shape {}, "
                         "but actually has {}.".format((n + k,), bc_res.shape))

    status = 0
    iteration = 0
    if verbose == 2:
        print_iteration_header()

    while True:
        m = x.shape[0]

        col_fun, jac_sys = prepare_sys(n, m, k, fun_wrapped, bc_wrapped, x, h)
        y, p, singular = solve_newton(
            n, m, k, col_fun, bc_wrapped, jac_sys, y, p, B, tol * 1e-1, 4)
        iteration += 1

        col_res, y_middle, f, f_middle = collocation_fun(fun_wrapped, y,
                                                         p, x, h)
        # This relation is not trivial, but can be verified.
        res_middle = 1.5 * col_res / h
        sol = create_spline(y, f, x, h)
        rms_res = estimate_rms_residuals(
            fun_wrapped, sol, x, h, p, res_middle, f_middle)
        max_res = np.max(rms_res)

        if singular:
            status = 2
            break

        insert_1, = np.nonzero((rms_res > tol) &
                               (rms_res < 100 * tol))
        insert_2, = np.nonzero(rms_res >= 100 * tol)
        inserted = insert_1.shape[0] + 2 * insert_2.shape[0]

        if m + inserted > max_nodes:
            status = 1
            if verbose == 2:
                inserted = "({})".format(inserted)
                print_iteration_progress(iteration, max_res, m, inserted)
            break

        if verbose == 2:
            print_iteration_progress(iteration, max_res, m, inserted)

        if inserted > 0:
            x = modify_mesh(x, insert_1, insert_2)
            h = np.diff(x)
            y = sol(x)
        else:
            status = 0
            break

    if verbose > 0:
        if status == 0:
            print("Solved in {} iterations, number of nodes {}, "
                  "maximum relative residual {:.2e}."
                  .format(iteration, x.shape[0], max_res))
        elif status == 1:
            print("Number of nodes is exceeded after iteration {}, "
                  "maximum relative residual {:.2e}."
                  .format(iteration, max_res))
        elif status == 2:
            print("Singular Jacobian encountered when solving the collocation "
                  "system on iteration {}, maximum relative residual {:.2e}."
                  .format(iteration, max_res))

    if p.size == 0:
        p = None

    return BVPResult(sol=sol, p=p, x=x, y=y, f=f, res=rms_res, niter=iteration,
                     status=status, message=TERMINATION_MESSAGES[status],
                     success=status == 0)
