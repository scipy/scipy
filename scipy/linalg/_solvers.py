"""Matrix equation solver routines"""

# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# February 24, 2012

# Modified: Chad Fulton <ChadFulton@gmail.com>
# June 19, 2014

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import inv, LinAlgError

from .basic import solve
from .lapack import get_lapack_funcs
from .decomp_schur import schur
from .special_matrices import kron

__all__ = ['solve_sylvester', 'solve_lyapunov', 'solve_discrete_lyapunov',
           'solve_continuous_are', 'solve_discrete_are']


def solve_sylvester(a,b,q):
    """
    Computes a solution (X) to the Sylvester equation (AX + XB = Q).

    Parameters
    ----------
    a : (M, M) array_like
        Leading matrix of the Sylvester equation
    b : (N, N) array_like
        Trailing matrix of the Sylvester equation
    q : (M, N) array_like
        Right-hand side

    Returns
    -------
    x : (M, N) ndarray
        The solution to the Sylvester equation.

    Raises
    ------
    LinAlgError
        If solution was not found

    Notes
    -----
    Computes a solution to the Sylvester matrix equation via the Bartels-
    Stewart algorithm.  The A and B matrices first undergo Schur
    decompositions.  The resulting matrices are used to construct an
    alternative Sylvester equation (``RY + YS^T = F``) where the R and S
    matrices are in quasi-triangular form (or, when R, S or F are complex,
    triangular form).  The simplified equation is then solved using
    ``*TRSYL`` from LAPACK directly.

    .. versionadded:: 0.11.0

    """

    # Compute the Schur decomp form of a
    r,u = schur(a, output='real')

    # Compute the Schur decomp of b
    s,v = schur(b.conj().transpose(), output='real')

    # Construct f = u'*q*v
    f = np.dot(np.dot(u.conj().transpose(), q), v)

    # Call the Sylvester equation solver
    trsyl, = get_lapack_funcs(('trsyl',), (r,s,f))
    if trsyl is None:
        raise RuntimeError('LAPACK implementation does not contain a proper Sylvester equation solver (TRSYL)')
    y, scale, info = trsyl(r, s, f, tranb='C')

    y = scale*y

    if info < 0:
        raise LinAlgError("Illegal value encountered in the %d term" % (-info,))

    return np.dot(np.dot(u, y), v.conj().transpose())


def solve_lyapunov(a, q):
    """
    Solves the continuous Lyapunov equation (AX + XA^H = Q) given the values
    of A and Q using the Bartels-Stewart algorithm.

    Parameters
    ----------
    a : array_like
        A square matrix

    q : array_like
        Right-hand side square matrix

    Returns
    -------
    x : array_like
        Solution to the continuous Lyapunov equation

    See Also
    --------
    solve_sylvester : computes the solution to the Sylvester equation

    Notes
    -----
    Because the continuous Lyapunov equation is just a special form of the
    Sylvester equation, this solver relies entirely on solve_sylvester for a
    solution.

    .. versionadded:: 0.11.0

    """

    return solve_sylvester(a, a.conj().transpose(), q)


def _solve_discrete_lyapunov_direct(a, q):
    """
    Solves the discrete Lyapunov equation directly.

    This function is called by the `solve_discrete_lyapunov` function with
    `method=direct`. It is not supposed to be called directly.
    """

    lhs = kron(a, a.conj())
    lhs = np.eye(lhs.shape[0]) - lhs
    x = solve(lhs, q.flatten())

    return np.reshape(x, q.shape)


def _solve_discrete_lyapunov_bilinear(a, q):
    """
    Solves the discrete Lyapunov equation using a bilinear transformation.

    This function is called by the `solve_discrete_lyapunov` function with
    `method=bilinear`. It is not supposed to be called directly.
    """
    eye = np.eye(a.shape[0])
    aH = a.conj().transpose()
    aHI_inv = inv(aH + eye)
    b = np.dot(aH - eye, aHI_inv)
    c = 2*np.dot(np.dot(inv(a + eye), q), aHI_inv)
    return solve_lyapunov(b.conj().transpose(), -c)


def solve_discrete_lyapunov(a, q, method=None):
    """
    Solves the discrete Lyapunov equation :math:`(A'XA-X=-Q)`.

    Parameters
    ----------
    a : (M, M) array_like
        A square matrix

    q : (M, M) array_like
        Right-hand side square matrix

    method : {'direct', 'bilinear'}, optional
        Type of solver.

        If not given, chosen to be ``direct`` if ``M`` is less than 10 and
        ``bilinear`` otherwise.

    Returns
    -------
    x : ndarray
        Solution to the discrete Lyapunov equation

    See Also
    --------
    solve_lyapunov : computes the solution to the continuous Lyapunov equation

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *direct* if ``M`` is less than 10
    and ``bilinear`` otherwise.

    Method *direct* uses a direct analytical solution to the discrete Lyapunov
    equation. The algorithm is given in, for example, [1]_. However it requires
    the linear solution of a system with dimension :math:`M^2` so that
    performance degrades rapidly for even moderately sized matrices.

    Method *bilinear* uses a bilinear transformation to convert the discrete
    Lyapunov equation to a continuous Lyapunov equation :math:`(B'X+XB=-C)`
    where :math:`B=(A-I)(A+I)^{-1}` and
    :math:`C=2(A' + I)^{-1} Q (A + I)^{-1}`. The continuous equation can be
    efficiently solved since it is a special case of a Sylvester equation.
    The transformation algorithm is from Popov (1964) as described in [2]_.

    .. versionadded:: 0.11.0

    References
    ----------
    .. [1] Hamilton, James D. Time Series Analysis, Princeton: Princeton
       University Press, 1994.  265.  Print.
       http://www.scribd.com/doc/20577138/Hamilton-1994-Time-Series-Analysis
    .. [2] Gajic, Z., and M.T.J. Qureshi. 2008.
       Lyapunov Matrix Equation in System Stability and Control.
       Dover Books on Engineering Series. Dover Publications.

    """
    if method is None:
        # Select automatically based on size of matrices
        if a.shape[0] >= 10:
            method = 'bilinear'
        else:
            method = 'direct'

    meth = method.lower()

    if meth == 'direct':
        x = _solve_discrete_lyapunov_direct(a, q)
    elif meth == 'bilinear':
        x = _solve_discrete_lyapunov_bilinear(a, q)
    else:
        raise ValueError('Unknown solver %s' % method)

    return x


def solve_continuous_are(a, b, q, r):
    """
    Solves the continuous algebraic Riccati equation, or CARE, defined
    as (A'X + XA - XBR^-1B'X+Q=0) directly using a Schur decomposition
    method.

    Parameters
    ----------
    a : (M, M) array_like
        Input
    b : (M, N) array_like
        Input
    q : (M, M) array_like
        Input
    r : (N, N) array_like
        Non-singular, square matrix

    Returns
    -------
    x : (M, M) ndarray
        Solution to the continuous algebraic Riccati equation

    See Also
    --------
    solve_discrete_are : Solves the discrete algebraic Riccati equation

    Notes
    -----
    Method taken from:
    Laub, "A Schur Method for Solving Algebraic Riccati Equations."
    U.S. Energy Research and Development Agency under contract
    ERDA-E(49-18)-2087.
    http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf

    .. versionadded:: 0.11.0

    """

    try:
        g = inv(r)
    except LinAlgError:
        raise ValueError('Matrix R in the algebraic Riccati equation solver is ill-conditioned')

    g = np.dot(np.dot(b, g), b.conj().transpose())

    z11 = a
    z12 = -1.0*g
    z21 = -1.0*q
    z22 = -1.0*a.conj().transpose()

    z = np.vstack((np.hstack((z11, z12)), np.hstack((z21, z22))))

    # Note: we need to sort the upper left of s to have negative real parts,
    #       while the lower right is positive real components (Laub, p. 7)
    [s, u, sorted] = schur(z, sort='lhp')

    (m, n) = u.shape

    u11 = u[0:m//2, 0:n//2]
    u21 = u[m//2:m, 0:n//2]
    u11i = inv(u11)

    return np.dot(u21, u11i)


def solve_discrete_are(a, b, q, r):
    """
    Solves the disctrete algebraic Riccati equation, or DARE, defined as
    (X = A'XA-(A'XB)(R+B'XB)^-1(B'XA)+Q), directly using a Schur decomposition
    method.

    Parameters
    ----------
    a : (M, M) array_like
        Non-singular, square matrix
    b : (M, N) array_like
        Input
    q : (M, M) array_like
        Input
    r : (N, N) array_like
        Non-singular, square matrix

    Returns
    -------
    x : ndarray
        Solution to the continuous Lyapunov equation

    See Also
    --------
    solve_continuous_are : Solves the continuous algebraic Riccati equation

    Notes
    -----
    Method taken from:
    Laub, "A Schur Method for Solving Algebraic Riccati Equations."
    U.S. Energy Research and Development Agency under contract
    ERDA-E(49-18)-2087.
    http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf

    .. versionadded:: 0.11.0

    """

    try:
        g = inv(r)
    except LinAlgError:
        raise ValueError('Matrix R in the algebraic Riccati equation solver is ill-conditioned')

    g = np.dot(np.dot(b, g), b.conj().transpose())

    try:
        ait = inv(a).conj().transpose()  # ait is "A inverse transpose"
    except LinAlgError:
        raise ValueError('Matrix A in the algebraic Riccati equation solver is ill-conditioned')

    z11 = a+np.dot(np.dot(g, ait), q)
    z12 = -1.0*np.dot(g, ait)
    z21 = -1.0*np.dot(ait, q)
    z22 = ait

    z = np.vstack((np.hstack((z11, z12)), np.hstack((z21, z22))))

    # Note: we need to sort the upper left of s to lie within the unit circle,
    #       while the lower right is outside (Laub, p. 7)
    [s, u, sorted] = schur(z, sort='iuc')

    (m,n) = u.shape

    u11 = u[0:m//2, 0:n//2]
    u21 = u[m//2:m, 0:n//2]
    u11i = inv(u11)

    return np.dot(u21, u11i)
