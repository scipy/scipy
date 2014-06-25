from __future__ import division, print_function, absolute_import
from ._ufuncs import _ellip_harm
from ._ellip_harm_2 import _ellipsoid, _ellipsoid_norm
import threading
import numpy as np

_ellip_lock = threading.Lock()

def ellip_harm(h2, k2, n, p, s, signm=1, signn=1):
    r"""
    Ellipsoidal Harmonic functions E^p_n(l), also known as Lame Functions:The first kind

    Lame's Equation is the following Differential Equation:

    .. math:: (s^2 - h^2)(s^2 - k^2)E''(s) + s(2s^2 - h^2 - k^2)E'(s) + (p - qs)E(s) = 0

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
    signm: double, optional
           determines the sign of prefactor of functions. See Notes 
    signn: double, optional
           determines the sign of prefactor of functions. See Notes
    
    Returns
    -------
    E^p_n(s) : double

    Notes
    -----
    Uses LAPACK subroutine DSTEVR
    The geometric intepretation is in accordance with [2],[3],[4]
    signm and signn control the sign of prefactor for functions according to their type.
    K: +1
    L: signm
    M: signn
    N: signm*signn

    References
    ----------
    .. [1] Digital Libary of Mathematical Functions 29.12
       http://dlmf.nist.gov/29.12
    .. [2] Bardhan and Knepley.Computational science and 
       re-discovery: open-source implementations of 
       ellipsoidal harmonics for problems in potential theory
       http://arxiv.org/abs/1204.0267
    .. [3] G. Romain and B. Jean-Pierre. Ellipsoidal harmonic expansions
       of the gravitational potential: theory and applications.
       http://link.springer.com/article/10.1023%2FA%3A1017555515763#close
    .. [4] David J.and Dechambre P. Computation of Ellipsoidal
       Gravity Field Harmonics for small solar system bodies
       http://ccar.colorado.edu/scheeres/scheeres/assets/Theses%20and%20Abstracts/dechambre_thesis.pdf
    
    Examples
    --------
    >>> from scipy.special import ellip_harm
    >>> w = ellip_harm(5,8,1,1,2.5)
    >>> w
    >>> 2.5

    """
    return _ellip_harm(h2, k2, n, p, s, signm, signn)

def ellip_harm_2(h2, k2, n, p, s):
    r"""
    Ellipsoidal Harmonic functions F^p_n(l), also known as Lame Functions:The second kind

    .. math::

     F^p_n(s)=(2n + 1)E^p_n(s)\int_{0}^{1/s}\frac{du}{(E^p_n(1/u))^2\sqrt{(1-u^2k^2)(1-u^2h^2)}}

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
    F^p_n(s) : double

    Notes
    -----
    The geometric intepretation is in accordance with [2],[3],[4]

    References
    ----------
    .. [1] Digital Libary of Mathematical Functions 29.12
       http://dlmf.nist.gov/29.12
    .. [2] Bardhan and Knepley.Computational science and 
       re-discovery: open-source implementations of 
       ellipsoidal harmonics for problems in potential theory
       http://arxiv.org/abs/1204.0267
    .. [3] G. Romain and B. Jean-Pierre. Ellipsoidal harmonic expansions
       of the gravitational potential: theory and applications.
       http://link.springer.com/article/10.1023%2FA%3A1017555515763#close
    .. [4] David J.and Dechambre P. Computation of Ellipsoidal
       Gravity Field Harmonics for small solar system bodies
       http://ccar.colorado.edu/scheeres/scheeres/assets/Theses%20and%20Abstracts/dechambre_thesis.pdf
    .. [5]George Dassios. Ellipsoidal Harmonics: Theory and Applications
    
    Examples
    --------
    >>> from scipy.special import ellip_harm_2
    >>> w = ellip_harm_2(5,8,2,1,10)
    >>> w
    >>> 0.00108056853382

    """
    with _ellip_lock:
        return _ellipsoid(h2, k2, n, p, s)

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
    \gamma^p_n(s) : double

    Notes
    -----
    The geometric intepretation is in accordance with [2],[3],[4]

    References
    ----------
    .. [1] Digital Libary of Mathematical Functions 29.12
       http://dlmf.nist.gov/29.12
    .. [2] Bardhan and Knepley.Computational science and 
       re-discovery: open-source implementations of 
       ellipsoidal harmonics for problems in potential theory
       http://arxiv.org/abs/1204.0267
    .. [3] G. Romain and B. Jean-Pierre. Ellipsoidal harmonic expansions
       of the gravitational potential: theory and applications.
       http://link.springer.com/article/10.1023%2FA%3A1017555515763#close
    .. [4] David J.and Dechambre P. Computation of Ellipsoidal
       Gravity Field Harmonics for small solar system bodies
       http://ccar.colorado.edu/scheeres/scheeres/assets/Theses%20and%20Abstracts/dechambre_thesis.pdf
    .. [5]George Dassios. Ellipsoidal Harmonics: Theory and Applications

    Examples
    --------
    >>> from scipy.special import ellip_harm_2
    >>> w = ellip_normal(5,8,3,7)
    >>> w
    >>> 1723.38796997

    """
    return _ellipsoid_norm(h2, k2, n, p)
