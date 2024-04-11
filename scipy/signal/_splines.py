
import numpy as np

from ._arraytools import axis_slice, axis_reverse
from ._signaltools import lfilter, sosfilt
from ._spline import symiirorder1_ic, symiirorder2_ic_fwd, symiirorder2_ic_bwd

__all__ = ['symiirorder1', 'symiirorder2']


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
        The input signal. If 2D, then the filter will be applied in a batched
        fashion across the last axis.
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

    if signal.ndim > 2:
        raise ValueError('Input must be 1D or 2D')

    squeeze_dim = False
    if signal.ndim == 1:
        signal = signal[None, :]
        squeeze_dim = True

    if np.issubdtype(signal.dtype, np.integer):
        signal = signal.astype(np.float64)

    y0 = symiirorder1_ic(signal, z1, precision)

    # Apply first the system 1 / (1 - z1 * z^-1)
    b = np.ones(1, dtype=signal.dtype)
    a = np.r_[1, -z1]
    a = a.astype(signal.dtype)

    # Compute the initial state for lfilter.
    zii = y0 * z1

    y1, _ = lfilter(b, a, axis_slice(signal, 1), zi=zii)
    y1 = np.c_[y0, y1]

    # Compute backward symmetric condition and apply the system
    # c0 / (1 - z1 * z)
    b = np.asarray([c0], dtype=signal.dtype)
    out_last = -c0 / (z1 - 1.0) * axis_slice(y1, -1)

    # Compute the initial state for lfilter.
    zii = out_last * z1

    out, _ = lfilter(b, a, axis_slice(y1, -2, step=-1), zi=zii)
    out = np.c_[axis_reverse(out), out_last]

    if squeeze_dim:
        out = out[0]

    return out


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

    if input.ndim > 2:
        raise ValueError('Input must be 1D or 2D')

    if not input.flags.c_contiguous:
        input = input.copy()

    squeeze_dim = False
    if input.ndim == 1:
        input = input[None, :]
        squeeze_dim = True

    if np.issubdtype(input.dtype, np.integer):
        input = input.astype(np.float64)

    rsq = r * r
    a2 = 2 * r * np.cos(omega)
    a3 = -rsq
    cs = np.atleast_1d(1 - 2 * r * np.cos(omega) + rsq)
    sos = np.atleast_2d(np.r_[cs, 0, 0, 1, -a2, -a3]).astype(input.dtype)

    # Find the starting (forward) conditions.
    ic_fwd = symiirorder2_ic_fwd(input, r, omega, precision)

    # Apply first the system cs / (1 - a2 * z^-1 - a3 * z^-2)
    # Compute the initial conditions in the form expected by sosfilt
    # coef = np.asarray([[a3, a2], [0, a3]], dtype=input.dtype)
    coef = np.r_[a3, a2, 0, a3].reshape(2, 2).astype(input.dtype)
    zi = np.matmul(coef, ic_fwd[:, :, None])[:, :, 0]

    y_fwd, _ = sosfilt(sos, axis_slice(input, 2), zi=zi[None])
    y_fwd = np.c_[ic_fwd, y_fwd]

    # Then compute the symmetric backward starting conditions
    ic_bwd = symiirorder2_ic_bwd(input, r, omega, precision)

    # Apply the system cs / (1 - a2 * z^1 - a3 * z^2)
    # Compute the initial conditions in the form expected by sosfilt
    zi = np.matmul(coef, ic_bwd[:, :, None])[:, :, 0]
    y, _ = sosfilt(sos, axis_slice(y_fwd, -3, step=-1), zi=zi[None])
    out = np.c_[axis_reverse(y), axis_reverse(ic_bwd)]

    if squeeze_dim:
        out = out[0]

    return out
