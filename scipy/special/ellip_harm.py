from __future__ import division, print_function, absolute_import

from ._ufuncs import _ellip_harm


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
