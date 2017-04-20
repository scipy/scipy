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


def fftlog(x, spacing, transform='sine', bias=0.0, kr=1.0, rk=1.0):
    """Fourier transform of a logarithmically spaced periodic sequence.

    Fast Fourier transform of a real, discrete periodic sequence of
    logarithmically spaced points.

    `fftlog` computes a discrete version of the Fourier sine or cosine
    transform

    .. math::

        G = \sqrt{2/\pi} \int_0^\infty F(r) \sin(kr) dr,

        G = \sqrt{2/\pi} \int_0^\infty F(r) \cos(kr) dr


    by making the substitutions

    .. math::

        F(r) = f(r) r^{ \mu - 1/2},

        G(k) = g(k) k^{-\mu - 1/2}

    and applying a biased Hankel transform to f(r);
    mu = 1/2 for the sine and -1/2 for the cosine transform.

    Parameters
    ----------
    x : array_like, real-valued
        Array F(r) to transform: f(j) is F(r_j) at r_j = r_c exp[(j-jc) dlnr],
        where jc = (n+1)/2 = central index of array.

    spacing : float, optional
        Separation between input-points (log10); may be positive or negative.
        Default is 0.01.

    transform : string, optional; {'sine', 'cosine'}
        Transform type to use, which defines index of J_mu in Hankel transform:
        mu is 0.5 for a sine transform and -0.5 for a cosine transform. Default
        is 'sine' (mu=0.5).

    bias : float, optional
        Exponent of power law bias; bias may be any real number, positive or
        negative. If in doubt, use bias = 0, for which case the Hankel transform
        is orthogonal, i.e. self-inverse, provided also that, for n even, kr is
        low-ringing. Non-zero bias may yield better approximations to the
        continuous Hankel transform for some functions. Default is 0
        (unbiased).

    kr : float, optional
        k_c r_c where c is central point of array
        = k_j r_(n+1-j) = k_(n+1-j) r_j .
        Normally one would choose kr to be about 1 (default) (or 2, or pi, to
        taste). Default is 1.

    rk : float, optional
        r_c/k_c = r_j/k_j (a constant, the same constant for any j); rk is not
        (necessarily) the same quantity as kr. rk is used only to multiply the
        output array by sqrt(rk)^dir, so if you want to do the normalization
        later, or you don't care about the normalization, you can set rk = 1.
        Default is 1.

    Returns
    -------
    y : real ndarray
        Transformed array G(k): g(j) is G(k_j) at k_j = k_c exp[(j-jc) dlnr].

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

    Get fftlog-arguments

    >>> n, spacing, center = 4, .1, 0
    >>> bias = 0
    >>> transform = 'sine'
    >>> w, t, kr, rk = fftlogargs(n, spacing, center, transform, bias, 1, 1)
    >>> rk /= 2/np.pi    # Scale

    Analytical solution

    >>> fw = np.sqrt(np.pi/2/w)  # Frequency domain
    >>> ft = 1/np.sqrt(t)        # Time domain

    FFTLog

    >>> fftl = fftlog(fw, spacing, transform, bias, kr, rk)
    >>> fftl *= 2/np.pi  # Scale back

    Print result

    >>> print('Input      :', fw)
    Input      : [ 1.48956664  1.32757767  1.18320484  1.05453243]
    >>> print('Analytical :', ft)
    Analytical : [ 1.15380264  1.02832769  0.91649802  0.81682972]
    >>> print('fftlog     :', fftl)
    fftlog     : [ 1.15380264  1.02832769  0.91649802  0.81682972]

    """

    # Check that transform is {'sine', or 'cosine'}
    if transform not in ['sine', 'cosine']:
        raise ValueError("transform must be either 'sine' or 'cosine'.")
    if transform == 'sine':
        mu = 0.5
    else:
        mu = -0.5

    tmp = _asfarray(x)

    if len(tmp) < 1:
        raise ValueError("Invalid number of FFT data points "
                         "(%d) specified." % len(tmp))

    dlnr = spacing*np.log(10.0)
    n = len(tmp)
    if np.iscomplexobj(tmp):  # Returns complex128
        y = (_fftl.drfftl(tmp.real, n, mu, bias, dlnr, kr, rk, 1) +
             1j*_fftl.drfftl(tmp.imag, n, mu, bias, dlnr, kr, rk, 1))
    else:  # Returns float64
        y = _fftl.drfftl(tmp, n, mu, bias, dlnr, kr, rk, 1)

    return y


def fftlogargs(n, spacing=0.01, center=0.0, transform='sine', bias=0, kr=1,
               kropt=False):
    """FFTLog input parameters (for usage with fftlog).

    Return the required input points and the corresponding output points, the
    (adjusted) kr and the corresponding rk for `fftlog`.

    Parameters
    ----------
    n : int
        Number of samples.

    spacing : float, optional
        Separation between input-points (log10); may be positive or negative.
        Default is 0.01.

    center : float, optional
        Central point of periodic interval (log10). Default is 0.

    transform : string, optional; {'sine', 'cosine'}
        Transform type to use, which defines index of J_mu in Hankel transform:
        mu is 0.5 for a sine transform and -0.5 for a cosine transform. Default
        is 'sine' (mu=0.5).

    bias : float, optional
        Exponent of power law bias; bias may be any real number, positive or
        negative. If in doubt, use bias = 0, for which case the Hankel transform
        is orthogonal, i.e. self-inverse, provided also that, for n even, kr is
        low-ringing. Non-zero bias may yield better approximations to the
        continuous Hankel transform for some functions. Only used if kropt is
        True. Default is 0 (unbiased).

    kr : float, optional
        k_c r_c where c is central point of array
        = k_j r_(n+1-j) = k_(n+1-j) r_j .
        Normally one would choose kr to be about 1 (default) (or 2, or pi, to
        taste). Default is 1.

    kropt : bool, optional
        - False to use input kr as is (default);
        - True to change kr to nearest low-ringing value.


    Returns
    -------
    inppts : ndarray
        Array of length `n`, containing the sample input-points.
    outpts : ndarray
        Array of length `n`, containing the corresponding output-points.
    kr : float
        Low-ringing kr if kropt=True; else same as input.
    rk : float
        Inverse of kr, shifted if center != 0 and kr != 1.

    .. versionadded:: 1.0.0

    Examples
    --------
    >>> from scipy.fftpack import fftlogargs
    >>> intpts, outpts, kr, rk = fftlogargs(n=4, kr=1, kropt=True)
    >>> print('intpts :', intpts)
    intpts : [ 0.96605088  0.98855309  1.01157945  1.03514217]
    >>> print('outpts :', outpts)
    outpts : [ 0.97306236  0.9957279   1.01892138  1.04265511]
    >>> print('kr     :', kr)
    kr     : 1.0072578812188107
    >>> print('rk     :', rk)
    rk     : 0.992794416054

    """
    # Check that transform is {'sine', or 'cosine'}; get mu
    if transform not in ['sine', 'cosine']:
        raise ValueError("transform must be either 'sine' or 'cosine'.")
    if transform == 'sine':
        mu = 0.5
    else:
        mu = -0.5

    # Central index (1/2 integral if n is even)
    nc = (n + 1)/2.0

    # Input points (frequencies)
    inppts = 10**(center + (np.arange(n)+1 - nc)*spacing)

    # Get low-ringing kr
    if kropt:
        dlnr = spacing*np.log(10.0)
        kr = _fftl.getkr(mu=mu, q=bias, dlnr=dlnr, kr=kr)

    # Central point log10(k_c) of periodic interval
    logkc = np.log10(kr) - center

    # rk = r_c/k_c
    rk = 10**(center - logkc)

    # Output points (times)
    outpts = 10**(logkc + (np.arange(n)+1 - nc)*spacing)

    return inppts, outpts, kr, rk
