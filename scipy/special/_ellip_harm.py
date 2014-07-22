from __future__ import division, print_function, absolute_import
from ._ufuncs import _ellip_harm
from ._ellip_harm_2 import _ellipsoid, _ellipsoid_norm
import threading
import numpy as np
# the functions _ellipsoid, _ellipsoid_norm use global variables, the lock  
# protects them if the function is called from multiple threads simultaneously
_ellip_lock = threading.Lock()

def ellip_harm(h2, k2, n, p, s, signm=1, signn=1):
    r"""
    Ellipsoidal Harmonic functions E^p_n(l), also known as Lame Functions:The first kind

    Lame's Equation is the following Differential Equation:

    .. math:: (s^2 - h^2)(s^2 - k^2)E''(s) + s(2s^2 - h^2 - k^2)E'(s) + (p - qs)E(s) = 0

    Parameters
    ----------
    h2 : double
        :math:`h^2`
    k2 : double
        :math:`k^2`
    n : int
       degree
    p : int
       order, can range between [1,2n+1]
    signm : double, optional
           determines the sign of prefactor of functions. See Notes 
    signn : double, optional
           determines the sign of prefactor of functions. See Notes
    
    Returns
    -------
    ellip_harm : double
        the harmonic :math:`E^p_n(s)`

    See Also
    --------
    ellip_harm2, ellip_normal

    Notes
    -----
    Uses LAPACK subroutine DSTEVR
    The geometric intepretation is in accordance with [2]_,[3]_,[4]_
    signm and signn control the sign of prefactor for functions according to their type.
    K : +1
    L : signm
    M : signn
    N : signm*signn

    References
    ----------
    .. [1] Digital Libary of Mathematical Functions 29.12
       http://dlmf.nist.gov/29.12
    .. [2] Bardhan and Knepley, "Computational science and 
       re-discovery: open-source implementations of 
       ellipsoidal harmonics for problems in potential theory",
       Comput. Sci. Disc. 5, 014006 (2012)
       doi:10.1088/1749-4699/5/1/014006
    .. [3] David J.and Dechambre P, "Computation of Ellipsoidal
       Gravity Field Harmonics for small solar system bodies"
       pp. 30-36, 2000
    .. [4] George Dassios, "Ellipsoidal Harmonics: Theory and Applications"
       pp. 418, 2012
    
    Examples
    --------
    >>> from scipy.special import ellip_harm
    >>> w = ellip_harm(5,8,1,1,2.5)
    >>> w
    2.5

    """
    return _ellip_harm(h2, k2, n, p, s, signm, signn)

# np.vectorize does not work on Cython functions on Numpy < 1.8, so a wrapper is needed
def _ellip_harm_2_vec(h2, k2, n, p, s):
    return _ellipsoid(h2, k2, n, p, s)

_ellip_harm_2_vec = np.vectorize(_ellip_harm_2_vec, otypes='d')

def ellip_harm_2(h2, k2, n, p, s):
    r"""
    Ellipsoidal Harmonic functions F^p_n(l), also known as Lame Functions:The second kind

    .. math::

     F^p_n(s)=(2n + 1)E^p_n(s)\int_{0}^{1/s}\frac{du}{(E^p_n(1/u))^2\sqrt{(1-u^2k^2)(1-u^2h^2)}}

    Parameters
    ----------
    h2 : double
        :math:`h^2`
    k2 : double
        :math:`k^2`
    n : int
       degree
    p : int
       order, can range between [1,2n+1]
    
    Returns
    -------
    ellip_harm_2 : double
        the harmonic :math:`F^p_n(s)`

    See Also
    --------
    ellip_harm

    Examples
    --------
    >>> from scipy.special import ellip_harm_2
    >>> w = ellip_harm_2(5,8,2,1,10)
    >>> w
    0.00108056853382

    """
    with _ellip_lock:
        return _ellip_harm_2_vec(h2, k2, n, p, s)

def _ellip_normal_vec(h2, k2, n, p):
    return _ellipsoid_norm(h2, k2, n, p)

_ellip_normal_vec = np.vectorize(_ellip_normal_vec, otypes='d')

def ellip_normal(h2, k2, n, p):
    r"""
    Normalization constant for Ellipsoidal Harmonic Functions: the first kind

    .. math:: 

    \gamma^p_n=8\int_{0}^{h}\int_{h}^{k}\frac{(y^2-x^2)(E^p_n(y)E^p_n(x))^2}{\sqrt((k^2-y^2)(y^2-h^2)(h^2-x^2)(k^2-x^2)}dydx

    Parameters
    ----------
    h2: double
        :math:`h^2`
    k2: double
        :math:`k^2`
    n: int
       degree
    p: int
       order, can range between [1,2n+1]

    Returns
    -------
    ellip_normal : double
        the normalization constant :math:`\gamma^p_n`

    See Also
    --------
    ellip_harm

    Examples
    --------
    >>> from scipy.special import ellip_harm_2
    >>> w = ellip_normal(5,8,3,7)
    >>> w
    1723.38796997

    """
    with _ellip_lock:
        return _ellip_normal_vec(h2, k2, n, p)
