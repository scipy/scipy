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

import numpy as np
from scipy._lib._util import _lazywhere
from scipy.special import iv, ive
from ._ufuncs import _log_iv_asym

__all__ = ['logiv', 'logive']


def _log_value(val, a):
    # Handle negative, 0, inf, nan
    # inf => inf, 0 => -inf, negative, nan => nan
    return _lazywhere(val > 0, (val,),
        f=lambda v: np.log(v),
        f2=lambda v: np.where(v == 0, -np.inf, np.nan))


def logiv(v, z):
    """
    Log of modified Bessel function of the first kind of real order.

    See `iv` for Parameters and return
    """

    # TODO: if logkv is added can support -v via DLMF 10.27.2

    v = np.asarray(v)
    z = np.asarray(z)

    value = _lazywhere(
        (v > 50) & (z > 0),
        (v, z),
        f=_log_iv_asym,
        f2=lambda v, z: _log_value(val=iv(v, z), a=(v, z)))

    # HACK: converts np.array(<scalar>) to <scalar>
    return value + 0


def logive(v, z):
    """
    Log of Exponentially scaled modified Bessel function of the first kind of real order.

    See `ive` for Parameters and return
    """

    v = np.asarray(v)
    z = np.asarray(z)

    # scaled by exp(-abs(z.real))
    value = _lazywhere(
        (v > 50) & (z > 0),
        (v, z),
        f=lambda v, z: _log_iv_asym(v, z) - abs(z),
        f2=lambda v, z: _log_value(val=ive(v, z), a=(v, z)))

    # HACK: converts np.array(<scalar>) to <scalar>
    return value + 0
