
import numpy as np

from ._arraytools import axis_slice, axis_reverse
from ._signaltools import lfiltic, lfilter
from ._spline import symiirorder1_ic

__all__ = ['symiirorder1']


def symiirorder1(signal, c0, z1, precision=-1.0):
    """
    Implement a smoothing IIR filter with mirror-symmetric boundary conditions
    using a cascade of first-order sections.  The second section uses a
    reversed sequence.  This implements a system with the following
    transfer function and mirror-symmetric boundary conditions::

                           c0
           H(z) = ---------------------
                   (1-z1/z) (1 - z1 z)

    The resulting signal will have mirror symmetric boundary conditions
    as well.

    Parameters
    ----------
    input : ndarray
        The input signal.
    c0, z1 : scalar
        Parameters in the transfer function.
    precision :
        Specifies the precision for calculating initial conditions
        of the recursive filter based on mirror-symmetric input.

    Returns
    -------
    output : ndarray
        The filtered signal.
    """
    if np.abs(z1) >= 1:
        raise ValueError('|z1| must be less than 1.0')

    y0 = symiirorder1_ic(signal, z1, precision)

    # Apply first the system 1 / (1 - z1 * z^-1)
    b = np.ones(1, dtype=signal.dtype)
    a = np.r_[1, -z1]
    zi = lfiltic(b, a, y0)

    y1, _ = lfilter(b, a, axis_slice(signal, 1), zi=zi)
    y1 = np.r_[y0, y1]

    # Compute backward symmetric condition and apply the system
    # c0 / (1 - z1 * z)
    b = np.asarray([c0])
    out_last = -c0 / (z1 - 1.0) * y1[-1]

    zi = lfiltic(b, a, np.atleast_1d(out_last))
    out, _ = lfilter(b, a, axis_slice(y1, -2, step=-1), zi=zi)
    return np.r_[axis_reverse(out), out_last]
