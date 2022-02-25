# Author: Pim Schellart
# 2010 - 2011

"""Tools for spectral analysis of unequally sampled signals."""

import numpy as np
cimport numpy as np
cimport cython

np.import_array()

__all__ = ['_lombscargle']


cdef extern from "math.h":
    double cos(double)
    double sin(double)
    double atan2(double, double)

@cython.boundscheck(False)
def _lombscargle(np.ndarray[np.float64_t, ndim=1] x,
                np.ndarray[np.float64_t, ndim=1] y,
                np.ndarray[np.float64_t, ndim=1] freqs):
    """
    _lombscargle(x, y, freqs)

    Computes the Lomb-Scargle periodogram.

    Parameters
    ----------
    x : array_like
        Sample times.
    y : array_like
        Measurement values (must be registered so the mean is zero).
    freqs : array_like
        Angular frequencies for output periodogram.

    Returns
    -------
    pgram : array_like
        Lomb-Scargle periodogram.

    Raises
    ------
    ValueError
        If the input arrays `x` and `y` do not have the same shape.

    See also
    --------
    lombscargle

    """

    # Check input sizes
    if x.shape[0] != y.shape[0]:
        raise ValueError("Input arrays do not have the same size.")

    # Create empty array for output periodogram
    pgram = np.empty(freqs.shape[0], dtype=np.float64)

    # Local variables
    cdef Py_ssize_t i, j
    cdef double c, s, xc, xs, cc, ss, cs
    cdef double tau, c_tau, s_tau, c_tau2, s_tau2, cs_tau

    for i in range(freqs.shape[0]):

        xc = 0.
        xs = 0.
        cc = 0.
        ss = 0.
        cs = 0.

        for j in range(x.shape[0]):

            c = cos(freqs[i] * x[j])
            s = sin(freqs[i] * x[j])

            xc += y[j] * c
            xs += y[j] * s
            cc += c * c
            ss += s * s
            cs += c * s

        tau = atan2(2 * cs, cc - ss) / (2 * freqs[i])
        c_tau = cos(freqs[i] * tau)
        s_tau = sin(freqs[i] * tau)
        c_tau2 = c_tau * c_tau
        s_tau2 = s_tau * s_tau
        cs_tau = 2 * c_tau * s_tau

        pgram[i] = 0.5 * (((c_tau * xc + s_tau * xs)**2 / \
            (c_tau2 * cc + cs_tau * cs + s_tau2 * ss)) + \
            ((c_tau * xs - s_tau * xc)**2 / \
            (c_tau2 * ss - cs_tau * cs + s_tau2 * cc)))

    return pgram
