import numpy as np
from warnings import warn
from ._basic import rfft, irfft
from ..special import loggamma, poch

__all__ = [
    'fht', 'ifht',
    'fhtoffset',
]


# constants
LN_2 = np.log(2)


def fht(a, dln, mu, offset=0.0, bias=0.0):

    # size of transform
    n = np.shape(a)[-1]

    # bias input array
    if bias != 0:
        # a_q(r) = a(r) (r/r_c)^{-q}
        j_c = (n-1)/2
        j = np.arange(n)
        a = a * np.exp(-bias*(j - j_c)*dln)

    # compute FHT coefficients
    u = fhtcoeff(n, dln, mu, offset=offset, bias=bias)

    # transform
    A = _fhtq(a, u)

    # bias output array
    if bias != 0:
        # A(k) = A_q(k) (k/k_c)^{-q} (k_c r_c)^{-q}
        A *= np.exp(-bias*((j - j_c)*dln + offset))

    return A


def ifht(A, dln, mu, offset=0.0, bias=0.0):

    # size of transform
    n = np.shape(A)[-1]

    # bias input array
    if bias != 0:
        # A_q(k) = A(k) (k/k_c)^{q} (k_c r_c)^{q}
        j_c = (n-1)/2
        j = np.arange(n)
        A = A * np.exp(bias*((j - j_c)*dln + offset))

    # compute FHT coefficients
    u = fhtcoeff(n, dln, mu, offset=offset, bias=bias)

    # transform
    a = _fhtq(A, u, inverse=True)

    # bias output array
    if bias != 0:
        # a(r) = a_q(r) (r/r_c)^{q}
        a /= np.exp(-bias*(j - j_c)*dln)

    return a


def fhtcoeff(n, dln, mu, offset=0.0, bias=0.0):
    '''Compute the coefficient array for a fast Hankel transform.
    '''

    lnkr, q = offset, bias

    # Hankel transform coefficients
    # u_m = (kr)^{-i 2m pi/(n dlnr)} U_mu(q + i 2m pi/(n dlnr))
    # with U_mu(x) = 2^x Gamma((mu+1+x)/2)/Gamma((mu+1-x)/2)
    xp = (mu+1+q)/2
    xm = (mu+1-q)/2
    y = np.linspace(0, np.pi*(n//2)/(n*dln), n//2+1)
    u = np.empty(n//2+1, dtype=complex)
    v = np.empty(n//2+1, dtype=complex)
    u.imag[:] = y
    u.real[:] = xm
    loggamma(u, out=v)
    u.real[:] = xp
    loggamma(u, out=u)
    y *= 2*(LN_2 - lnkr)
    u.real -= v.real
    u.real += LN_2*q
    u.imag += v.imag
    u.imag += y
    np.exp(u, out=u)

    # fix last coefficient to be real
    u.imag[-1] = 0

    # deal with special cases
    if not np.isfinite(u[0]):
        # write u_0 = 2^q Gamma(xp)/Gamma(xm) = 2^q poch(xm, xp-xm)
        # poch() handles special cases for negative integers correctly
        u[0] = 2**q * poch(xm, xp-xm)
        # the coefficient may be inf or 0, meaning the transform or the
        # inverse transform, respectively, is singular

    return u


def fhtoffset(dln, mu, initial=0.0, bias=0.0):
    '''Return optimal offset for a fast Hankel transform.

    Returns an offset close to `initial` that fulfils the low-ringing
    condition of [1]_ for the fast Hankel transform `fht` with logarithmic
    spacing `dln`, order `mu` and bias `bias`.

    Parameters
    ----------
    dln : float
        Uniform logarithmic spacing of the transform.
    mu : float
        Order of the Hankel transform, any positive or negative real number.
    initial : float, optional
        Initial value for the offset. Returns the closest value that fulfils
        the low-ringing condition.
    bias : float, optional
        Exponent of power law bias, any positive or negative real number.

    Returns
    -------
    offset : float
        Optimal offset of the uniform logarithmic spacing of the transform that
        fulfils a low-ringing condition.

    See Also
    --------
    fht : Definition of the fast Hankel transform.

    References
    ----------
    .. [1] Hamilton A. J. S., 2000, MNRAS, 312, 257 (astro-ph/9905191)

    '''

    lnkr, q = initial, bias

    xp = (mu+1+q)/2
    xm = (mu+1-q)/2
    y = np.pi/(2*dln)
    zp = loggamma(xp + 1j*y)
    zm = loggamma(xm + 1j*y)
    arg = (LN_2 - lnkr)/dln + (zp.imag + zm.imag)/np.pi
    return lnkr + (arg - np.round(arg))*dln


def _fhtq(a, u, inverse=False):
    '''Compute the biased fast Hankel transform.

    This is the basic FFTLog routine.
    '''

    # size of transform
    n = np.shape(a)[-1]

    # check for singular transform or singular inverse transform
    if np.isinf(u[0]) and not inverse:
        warn('singular transform; consider changing the bias')
        # fix coefficient to obtain (potentially correct) transform anyway
        u = u.copy()
        u[0] = 0
    elif u[0] == 0 and inverse:
        warn('singular inverse transform; consider changing the bias')
        # fix coefficient to obtain (potentially correct) inverse anyway
        u = u.copy()
        u[0] = np.inf

    # biased fast Hankel transform via real FFT
    A = rfft(a, axis=-1)
    if not inverse:
        # forward transform
        A *= u
    else:
        # backward transform
        A /= u.conj()
    A = irfft(A, n, axis=-1)
    A = A[..., ::-1]

    return A
