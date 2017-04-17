"""
FFTLog - fftlog.py
"""
# Created by Dieter Werthmüller, January 2017
from __future__ import division, print_function, absolute_import

import atexit
import numpy as np
from . import _fftl
from .basic import _asfarray

__all__ = ['fftlog', 'fftlogargs']

atexit.register(_fftl.destroy_fftl_x_cache)
atexit.register(_fftl.destroy_fftl_w_cache)
del atexit


def fftlog(x, dlogr, mu='sine', q=0.0, kr=1.0, rk=1.0):
    """Fourier transform of a logarithmically spaced periodic sequence.

    Fast Fourier transform of a real, discrete periodic sequence of
    logarithmically spaced points.

    `fftlog` computes a discrete version of the Fourier sine (if mu = 1/2) or
    cosine (if mu = -1/2) transform

    .. math::

        \tilde{A} = \sqrt{2/\pi} \int_0^\infty A(r) \sin(kr) dr,

        \tilde{A} = \sqrt{2/\pi} \int_0^\infty A(r) \cos(kr) dr


    by making the substitutions

    .. math::

        A(r) = a(r) r^{\mu-1/2},

        \tilde{A}(k) = \tilde{a}(k) k^{-\mu -1/2}

    and applying a biased Hankel transform to a(r).

    Parameters
    ----------
    x : array_like, real-valued
        Array A(r) to transform: a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr],
        where jc = (n+1)/2 = central index of array.

    dlogr : float, optional
        Separation between input-points (log10); may be positive or negative.
        Default is 0.01.

    mu : string, optional; {'sine', 'cosine'}
        Index of J_mu in Hankel transform: 0.5 for a sine transform and -0.5
        for a cosine transform. Default is 'sine' (0.5).

    q : float, optional
        Exponent of power law bias; q may be any real number, positive or
        negative.  If in doubt, use q = 0, for which case the Hankel transform
        is orthogonal, i.e. self-inverse, provided also that, for n even, kr is
        low-ringing.  Non-zero q may yield better approximations to the
        continuous Hankel transform for some functions. Default is 0
        (unbiased).

    kr : float, optional
        k_c r_c where c is central point of array
        = k_j r_(n+1-j) = k_(n+1-j) r_j .
        Normally one would choose kr to be about 1 (default) (or 2, or pi, to
        taste). Default is 1.

    rk : float, optional
        r_c/k_c = r_j/k_j (a constant, the same constant for any j); rk is not
        (necessarily) the same quantity as kr.  rk is used only to multiply the
        output array by sqrt(rk)^dir, so if you want to do the normalization
        later, or you don't care about the normalization, you can set rk = 1.
        Default is 1.

    n : int, optional
        Defines the length of the Fourier transform.  If `n` is not specified
        (the default) then ``n = x.shape[axis]``.  If ``n < x.shape[axis]``,
        `x` is truncated, if ``n > x.shape[axis]``, `x` is zero-padded.

    axis : int, optional
        The axis along which the transform is applied.  The default is the
        last axis.

    overwrite_x : bool, optional
        If set to true, the contents of `x` can be overwritten. Default is
        False.

    Returns
    -------
    y : real ndarray
        Transformed array Ã(k): a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].

    .. versionadded:: 1.0.0

    References
    ----------
    .. [1] 'Uncorrelated modes of the non-linear power spectrum', by A. J. S.
           Hamilton, `Monthly Notices of the Royal Astronomical Society` vol.
           312, pp. 257-284, http://dx.doy.org/10.1046/j.1365-8711.2000.03071.x
           (2000). Website of FFTLog: http://casa.colorado.edu/~ajsh/FFTLog.

    Examples
    --------
    >>> from scipy.fftpack import fftlog, fftlogargs
    >>> # Get fftlog-arguments
    >>> n, dlogr, logrc = 4, .1, 0
    >>> q = 0
    >>> mu = 'sine'
    >>> w, t, kr, rk = fftlogargs(n, dlogr, logrc, mu, q, 1, 1)
    >>> rk /= 2/np.pi    # Scale
    >>> # Analytical solution
    >>> fw = np.sqrt(np.pi/2/w)  # Frequency domain
    >>> ft = 1/np.sqrt(t)        # Time domain
    >>> # FFTLog
    >>> fftl = fftlog(fw, dlogr=dlogr, mu=mu, q=q, kr=kr, rk=rk)
    >>> fftl *= 2/np.pi  # Scale back
    >>> # Print result
    >>> print('Input      :', fw)
    >>> print('Analytical :', ft)
    >>> print('fftlog     :', fftl)
    Input      : [ 1.48956664  1.32757767  1.18320484  1.05453243]
    Analytical : [ 1.15380264  1.02832769  0.91649802  0.81682972]
    fftlog     : [ 1.15380264  1.02832769  0.91649802  0.81682972]

    """

    # Check that mu is {'sine', or 'cosine'}
    if mu not in ['sine', 'cosine']:
        raise ValueError("mu must be either 'sine' or 'cosine'.")
    if mu == 'sine':
        nmu = 0.5
    else:
        nmu = -0.5

    tmp = _asfarray(x)

    if len(tmp) < 1:
        raise ValueError("Invalid number of FFT data points "
                         "(%d) specified." % len(tmp))

    dlnr = dlogr*np.log(10.0)
    n = len(tmp)
    if np.iscomplexobj(tmp):  # Returns complex128
        y = (_fftl.drfftl(tmp.real, n, nmu, q, dlnr, kr, rk, 1) +
             1j*_fftl.drfftl(tmp.imag, n, nmu, q, dlnr, kr, rk, 1))
    else:  # Returns float64
        y = _fftl.drfftl(tmp, n, nmu, q, dlnr, kr, rk, 1)

    return y


def fftlogargs(n, dlogr=0.01, logrc=0.0, mu='sine', q=0, kr=1, kropt=0):
    """FFTLog input parameters (for usage with fftlog).

    Return the required input points and the corresponding output points, the
    (adjusted) kr and the corresponding rk for `fftlog`.

    Parameters
    ----------
    n : int
        Number of samples.

    dlogr : float, optional
        Separation between input-points (log10); may be positive or negative.
        Default is 0.01.

    logrc : float, optional
        Central point of periodic interval (log10). Default is 0.

    mu : string, optional; {'sine', 'cosine'}
        Index of J_mu in Hankel transform: 0.5 for a sine transform and -0.5
        for a cosine transform. Default is 'sine' (0.5).

    q : float, optional
        Exponent of power law bias; q may be any real number, positive or
        negative.  If in doubt, use q = 0, for which case the Hankel transform
        is orthogonal, i.e. self-inverse, provided also that, for n even, kr is
        low-ringing.  Non-zero q may yield better approximations to the
        continuous Hankel transform for some functions. Only used if kropt is
        1. Default is 0 (unbiased).

    kr : float, optional
        k_c r_c where c is central point of array
        = k_j r_(n+1-j) = k_(n+1-j) r_j .
        Normally one would choose kr to be about 1 (default) (or 2, or pi, to
        taste). Default is 1.

    kropt : int, optional; {0, 1}
        - 0 to use input kr as is (default);
        - 1 to change kr to nearest low-ringing kr, quietly.


    Returns
    -------
    inppts : ndarray
        Array of length `n`, containing the sample input-points.
    outpts : ndarray
        Array of length `n`, containing the corresponding output-points.
    kr : float
        Low-ringing kr if kropt=1; else same as input.
    rk : float
        Inverse of kr, shifted if logrc != 0 and kr != 1.

    .. versionadded:: 1.0.0

    Examples
    --------
    >>> from scipy import fftpack
    >>> n = 4
    >>> kr = 1
    >>> kropt = 1
    >>> outpts = fftpack.fftlogargs(n, kr=kr, kropt=kropt)
    >>> print('intpts :', outpts[0])
    >>> print('outpts :', outpts[1])
    >>> print('kr     :', outpts[2])
    >>> print('rk     :', outpts[3])
    intpts : [ 0.96605088  0.98855309  1.01157945  1.03514217]
    outpts : [ 0.97306236  0.9957279   1.01892138  1.04265511]
    kr     : 1.0072578812188107
    rk     : 0.992794416054

    """
    # Check that mu is {'sine', or 'cosine'}; get numeric mu (nmu)
    if mu not in ['sine', 'cosine']:
        raise ValueError("mu must be either 'sine' or 'cosine'.")
    if mu == 'sine':
        nmu = 0.5
    else:
        nmu = -0.5

    # Central index (1/2 integral if n is even)
    nc = (n + 1)/2.0

    # Input points (frequencies)
    inppts = 10**(logrc + (np.arange(n)+1 - nc)*dlogr)

    # Get low-ringing kr
    if kropt == 1:
        dlnr = dlogr*np.log(10.0)
        kr = _fftl.getkr(mu=nmu, q=q, dlnr=dlnr, kr=kr, kropt=kropt)

    # Central point log10(k_c) of periodic interval
    logkc = np.log10(kr) - logrc

    # rk = r_c/k_c
    rk = 10**(logrc - logkc)

    # Output points (times)
    outpts = 10**(logkc + (np.arange(n)+1 - nc)*dlogr)

    return inppts, outpts, kr, rk
