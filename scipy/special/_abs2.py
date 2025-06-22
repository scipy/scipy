import numpy as np

from ._ufuncs import _abs_sq


def abs_sq(z):
    """
    Absolute square of z.
    Equivalent to (but faster than) ``z.real**2 + z.imag**2``.

    Parameters
    ----------
    z : array_like
        Real or complex argument

    Returns
    -------
    scalar or ndarray
        Squared absolute of the argument
        (``x**2`` for reals, ``z.real**2 + z.imag**2`` for complex).

    Examples
    --------
    >>> from scipy import special
    >>> special.abs_sq(2.)
    4.0
    >>> special.abs_sq([1., 2., 3.])
    array([1., 4., 9.])
    >>> special.abs_sq(1+1j)
    2.0
    >>> special.abs_sq([1+1j, 2+2j, 3+3j])
    array([ 2.,  8., 18.])
    """
    if not isinstance(z, np.ndarray):
        z = np.asanyarray(z, order='K')
    return _abs_sq(z, z.flat[0].real)
