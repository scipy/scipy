# Author: Robert T. McGibbon, December 2014

from numpy import zeros, sum

# FIXME: from numpy.linalg import LinAlgError when it's supported by Pythran
LinAlgError = RuntimeError

# FIXME: use '@' once pythran supports '@' without linking to blas.
def _dot(x, y):
    return sum(x * y)


# pythran export levinson(float64 [], float64[])
# pythran export levinson(complex128 [], complex128[])
def levinson(a, b):
    """Solve a linear Toeplitz system using Levinson recursion.

    Parameters
    ----------
    a : array, dtype=double or complex128, shape=(2n-1,)
        The first column of the matrix in reverse order (without the diagonal)
        followed by the first (see below)
    b : array, dtype=double  or complex128, shape=(n,)
        The right hand side vector. Both a and b must have the same type
        (double or complex128).

    Notes
    -----
    For example, the 5x5 toeplitz matrix below should be represented as
    the linear array ``a`` on the right ::

        [ a0    a1   a2  a3  a4 ]
        [ a-1   a0   a1  a2  a3 ]
        [ a-2  a-1   a0  a1  a2 ] -> [a-4  a-3  a-2  a-1  a0  a1  a2  a3  a4]
        [ a-3  a-2  a-1  a0  a1 ]
        [ a-4  a-3  a-2  a-1 a0 ]

    Returns
    -------
    x : arrray, shape=(n,)
        The solution vector
    reflection_coeff : array, shape=(n+1,)
        Toeplitz reflection coefficients. When a is symmetric Toeplitz and
        ``b`` is ``a[n:]``, as in the solution of autoregressive systems,
        then ``reflection_coeff`` also correspond to the partial
        autocorrelation function.
    """
    # Adapted from toeplitz.f90 by Alan Miller, accessed at
    # http://jblevins.org/mirror/amiller/toeplitz.f90
    # Released under a Public domain declaration.

    dtype = a.dtype

    (n,) = b.shape
    x = zeros(n, dtype=dtype)  # result
    g = zeros(n, dtype=dtype)  # workspace
    h = zeros(n, dtype=dtype)  # workspace
    reflection_coeff = zeros(n + 1, dtype=dtype)  # history

    assert len(a) == (2 * n) - 1

    if a[n - 1] == 0:
        raise LinAlgError("Singular principal minor")

    x[0] = b[0] / a[n - 1]
    reflection_coeff[0] = 1
    reflection_coeff[1] = x[0]

    if n == 1:
        return x, reflection_coeff

    g[0] = a[n - 2] / a[n - 1]
    h[0] = a[n] / a[n - 1]

    for m in range(1, n):
        # Compute numerator and denominator of x[m]
        x_num = _dot(a[n + m - 1 : n - 1 : -1], x[:m]) - b[m]
        x_den = _dot(a[n + m - 1 : n - 1 : -1], g[m - 1 :: -1]) - a[n - 1]
        if x_den == 0:
            raise LinAlgError("Singular principal minor")
        x[m] = x_num / x_den
        reflection_coeff[m + 1] = x[m]

        # Compute x
        x[:m] -= x[m] * g[m - 1 :: -1]
        if m == n - 1:
            return x, reflection_coeff

        # Compute the numerator and denominator of g[m] and h[m]
        g_num = _dot(a[n - m - 1 : n - 1], g[:m]) - a[n - m - 2]
        h_num = _dot(a[n + m - 1 : n - 1 : -1], h[:m]) - a[n + m]
        g_den = _dot(a[n - m - 1 : n - 1], h[m - 1 :: -1]) - a[n - 1]

        if g_den == 0.0:
            raise LinAlgError("Singular principal minor")

        # Compute g and h
        g[m] = g_num / g_den
        h[m] = h_num / x_den
        k = m - 1
        m2 = (m + 1) // 2
        c1 = g[m]
        c2 = h[m]
        for j in range(m2):
            gj = g[j]
            gk = g[k]
            hj = h[j]
            hk = h[k]
            g[j] = gj - (c1 * hk)
            g[k] = gk - (c1 * hj)
            h[j] = hj - (c2 * gk)
            h[k] = hk - (c2 * gj)
            k -= 1
