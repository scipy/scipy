"""
Implementation of the Levinson-Durbin algorithm for solving Toeplitz systems.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg


__all__ = ['solve_toeplitz']


def solve_toeplitz(c, r=None, y=None):
    """
    Solve the matrix equation (T x = y) where T is a Toeplitz matrix.

    The square Toeplitz input matrix T is represented through its
    first column and optionally its first row.
    This representation, including the ordering of the first two args,
    is taken from the linalg.toeplitz function.

    Parameters
    ----------
    c : 1d array_like
        The first column of the Toeplitz matrix.
    r : 1d array_like, optional
        The first row of the Toeplitz matrix.
    y : array_like
        The rhs of the matrix equation.

    Returns
    -------
    x : ndarray
        The solution of the matrix equation (T x = y).

    Notes
    -----
    This is an implementation of the Levinson-Durbin algorithm which uses
    only O(N) memory and O(N^2) time, but which has less than stellar
    numerical stability.

    References
    ----------
    .. [1] Wikipedia, "Levinson recursion",
           http://en.wikipedia.org/wiki/Levinson_recursion

    """
    # This block has been copied from the linalg.toeplitz construction function.
    c = np.asarray(c).ravel()
    if r is None:
        r = c.conjugate()
    else:
        r = np.asarray(r).ravel()

    # Check that the rhs makes sense.
    if y is None:
        raise ValueError('missing rhs')
    y = np.asarray(y)
    N = y.shape[0]

    # Check that the Toeplitz representation makes sense
    # and is compatible with the rhs shape.
    if c.shape != r.shape:
        raise ValueError('expected the Toeplitz matrix to be square')
    if c.shape != (N,):
        raise ValueError('the rhs shape is incompatible with the matrix shape')

    # If the diagonal is zero, then the Levinson-Durbin implementation fails.
    if not c[0]:
        raise np.linalg.LinAlgError(
                'the scipy implementation of the Levinson-Durbin algorithm '
                'fails when the main diagonal is zero')

    # Key relating entries of the toeplitz matrix to entries of c, r,
    # assuming n is a positive integer less than N:
    # M[0, 0] == c[0]
    # M[n, :n] == c[n:0:-1]
    # M[0, 1:n+1] == r[1:n+1]

    # If any of the input arrays are complex then use complex dtype.
    # Otherwise use real dtype.
    if any(np.iscomplexobj(obj) for obj in (c, r, y)):
        mytype = np.complex128
    else:
        mytype = np.float64

    # Initialize the forward, backward, and solution vectors.
    f_prev = np.zeros(N, dtype=mytype)
    b_prev = np.zeros(N, dtype=mytype)
    x_prev = np.zeros(y.shape, dtype=mytype)
    f = np.zeros(N, dtype=mytype)
    b = np.zeros(N, dtype=mytype)
    x = np.zeros(y.shape, dtype=mytype)
    f[0] = 1 / c[0]
    b[0] = 1 / c[0]
    x[0] = y[0] / c[0]

    # Compute forward, backward, and solution vectors recursively.
    # The tensordot shenanigans is for multiple rhs.
    for n in range(1, N):
        f, f_prev = f_prev, f
        b, b_prev = b_prev, b
        x, x_prev = x_prev, x
        eps_f = np.dot(c[n:0:-1], f_prev[:n])
        eps_x = np.tensordot(c[n:0:-1], x_prev[:n], 1)
        eps_b = np.dot(r[1:n+1], b_prev[:n])
        f.fill(0)
        b.fill(0)
        x.fill(0)
        denom = 1 - eps_f * eps_b

        # Complain if the denominator is exactly zero.
        # For better numerical stability, maybe use a different algorithm.
        if not denom:
            raise np.linalg.LinAlgError(
                    'the Levinson-Durbin algorithm '
                    'failed to solve the matrix equation')

        coeff = 1 / denom
        f[:n] += coeff * f_prev[:n]
        f[1:n+1] -= coeff * eps_f * b_prev[:n]
        b[1:n+1] += coeff * b_prev[:n]
        b[:n] -= coeff * eps_b * f_prev[:n]
        x[:n+1] = x_prev[:n+1] + np.tensordot(b[:n+1], (y[n] - eps_x), 0)

    return x

