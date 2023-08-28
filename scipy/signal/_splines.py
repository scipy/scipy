
import numpy as np

from ._arraytools import axis_slice, axis_reverse
from ._signaltools import lfiltic, lfilter, sosfilt
from ._spline import symiirorder1_ic, symiirorder2_ic_fwd, symiirorder2_ic_bwd

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


def symiirorder2(input, r, omega, precision=-1.0):
    """
    Implement a smoothing IIR filter with mirror-symmetric boundary conditions
    using a cascade of second-order sections.  The second section uses a
    reversed sequence.  This implements the following transfer function::

                                  cs^2
         H(z) = ---------------------------------------
                (1 - a2/z - a3/z^2) (1 - a2 z - a3 z^2 )

    where::
          a2 = 2 * r * cos(omega)
          a3 = - r ** 2
          cs = 1 - 2 * r * cos(omega) + r ** 2

    Parameters
    ----------
    input : ndarray
        The input signal.
    r, omega : float
        Parameters in the transfer function.
    precision : float
        Specifies the precision for calculating initial conditions
        of the recursive filter based on mirror-symmetric input.

    Returns
    -------
    output : ndarray
        The filtered signal.
    """
    if r >= 1.0:
        raise ValueError('r must be less than 1.0')

    if not input.flags.c_contiguous:
        input = input.copy()

    rsq = r * r
    a2 = 2 * r * np.cos(omega)
    a3 = -rsq
    cs = np.atleast_1d(1 - 2 * r * np.cos(omega) + rsq)

    # Find the starting (forward) conditions.
    ic_fwd = symiirorder2_ic_fwd(input, r, omega, precision)

    # Apply first the system cs / (1 - a2 * z^-1 - a3 * z^-2)
    b = cs
    a = np.r_[1, -a2, -a3]
    zi = lfiltic(b, a, ic_fwd)
    sos = np.r_[cs, 0, 0, 1, -a2, -a3]
    y_fwd, _ = sosfilt(sos, input[2:], zi=zi)
    y_fwd = np.r_[ic_fwd, y_fwd]

    # Then compute the symmetric backward starting conditions
    ic_bwd = symiirorder2_ic_bwd(input, r, omega, precision)

    # Apply the system cs / (1 - a2 * z^1 - a3 * z^2)
    zi = lfiltic(b, a, ic_bwd)
    y, _ = sosfilt(sos, axis_slice(y_fwd, -3, step=-1), zi=zi)
    return np.r_[axis_reverse(y), ic_bwd]
