"""
Log of Bessel functions

References
----------
.. [dlmf] NIST Digital Library of Mathematical Functions
       https://dlmf.nist.gov/10.41#ii
.. [a144617] Sequence A144617, The On-Line Encyclopedia of Integer Sequences,
       https://oeis.org/A144617
"""
#
# Author:  Seth Troisi 2020

from numpy import log
from scipy.special import iv, ive
from ._logiv import log_iv_asym

__all__ = ['logiv', 'logive']

def logiv(v, z):
    """
    Log of modified Bessel function of the first kind of real order.

    See `iv` for Parameters and return
    """

    if abs(v) <= 50:
        return log(iv(v, z))

    return log_iv_asym(v, z)

def logive(v, z):
    """
    Log of Exponentially scaled modified Bessel function of the first kind of real order.

    See `ive` for Parameters and return
    """

    if abs(v) <= 50:
        return log(ive(v, z))

    # scaled by exp(-abs(z.real))
    return log_iv_asym(v, z) - abs(z)
