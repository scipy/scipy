from numpy import zeros, asarray, polynomial, r_
from scipy import linalg

__all__ = ["pade"]

_Polynomial = polynomial.Polynomial


def pade(an, m, n=None):
    """
    Return Pade approximation to a polynomial as the ratio of two polynomials.

    Parameters
    ----------
    an : (N,) array_like
        Taylor series coefficients.
    m : int
        The order of the returned approximating polynomial `q`.
    n : int, optional
        The order of the returned approximating polynomial `p`. By default,
        the order is ``len(an)-1-m``.

    Returns
    -------
    p, q : Polynomial class
        The Pade approximation of the polynomial defined by `an` is
        ``p(x)/q(x)``.

    Examples
    --------
    >>> from scipy.interpolate import pade
    >>> e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]
    >>> p, q = pade(e_exp, 2)

    >>> e_poly = np.polynomial.Polynomial(e_exp)

    Compare ``e_poly(x)`` and the Pade approximation ``p(x)/q(x)``

    >>> e_poly(1)
    2.7166666666666668

    >>> p(1)/q(1)
    2.7179487179487181

    """
    an = asarray(an)
    if n is None:
        n = len(an) - 1 - m
        if n < 0:
            raise ValueError("Order of q <m> must be smaller than len(an)-1.")
    if n < 0:
        raise ValueError("Order of p <n> must be greater than 0.")
    N = m + n
    if N > len(an)-1:
        raise ValueError("Order of q+p <m+n> must be smaller than len(an).")
    an = an[:N+1]
    if m == 0:  # this is just the Taylor series
        return _Polynomial(an), _Polynomial([1])
    # first solve the Toeplitz system for q
    # first row might contain tailing zeros
    top = r_[an[n+1::-1][:m+1], [0]*(m-n-1)]
    an_mat = linalg.toeplitz(an[n+1:], top)
    # we set q[0] = 1 -> first column is -rhs
    q = r_[1.0, linalg.solve(an_mat[:, 1:], -an_mat[:, 0])]
    # substitute `q` to get `p`
    if m > 100:  # arbitrary threshold when to use dedicated Toeplitz function
        p = linalg.matmul_toeplitz((an[:n+1], zeros(m+1)), q)
    else:
        p = linalg.toeplitz(an[:n+1], zeros(m+1)) @ q
    return _Polynomial(p), _Polynomial(q)
