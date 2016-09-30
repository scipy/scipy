"""
Rational interpolation in barycentric form
"""
from __future__ import division, print_function, absolute_import

import sys
import warnings
import itertools

import numpy as np
from numpy import linspace, pi
import scipy.linalg

from scipy.interpolate import BarycentricInterpolator
from scipy._lib.six import xrange


__all__ = ['floater_hormann', 'lsq_rational',
           'RationalInterpolationWarning', 'RationalBarycentricInterpolation']


class RationalInterpolationWarning(Warning):
    pass


class RationalBarycentricInterpolation(object):
    r"""
    Rational function interpolator in barycentric form

    .. math::

       Q(z) = \frac{\sum_j \frac{u_j y_j}{x_j - z}}{\sum_j \frac{u_j}{x_j - z}}

    Parameters
    ----------
    x : array_like, shape (N,)
        Interpolation points
    y : array_like, shape (..., N)
        Data at interpolation points
    u : array_like, shape (..., N)
        Weights of interpolation points

    See Also
    --------
    floater_hormann, lsq_rational

    Methods
    -------
    __call__
    get_numerator
    get_denominator

    Attributes
    ----------
    x
    y
    u
    numerator_roots
    denominator_roots
    """

    def __init__(self, x, y, u):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.u = np.asarray(u)
        self._dtype = np.find_common_type([self.x.dtype, self.y.dtype, self.u.dtype, np.inexact], [])
        if self.x.ndim != 1 or self.x.size != self.y.shape[-1]:
            raise ValueError("x and y have incompatible shapes")
        if self.u.shape[-1] != self.y.shape[-1]:
            raise ValueError("y and u have incompatible shapes")
        self._numerator_roots = None
        self._denominator_roots = None

    def __call__(self, xp):
        xp = np.asarray(xp)
        xp_shape = xp.shape
        xp = xp.ravel()

        x = self.x
        y = self.y

        r = np.empty(y.shape[:-1] + (xp.size,), dtype=self._dtype)

        qs = xp - x[:,None]

        # Deal with exact zeros
        k = np.argmin(abs(qs), axis=0)
        j = np.arange(len(xp))
        m = (qs[k,j] == 0)
        r[...,m] = y[...,k[m]]

        # Other points: barycentric interpolation
        qs = qs[...,~m]
        np.reciprocal(qs, out=qs)
        numerator = np.einsum('...i,...i,...ij->...j', y, self.u, qs)
        denominator = np.einsum('...i,...ij->...j', self.u, qs)
        r[...,~m] = numerator/denominator

        return r.reshape(y.shape[:-1] + xp_shape)

    @staticmethod
    def _polyroots(x, u, dtype):
        r"""
        Find the roots of the polynomial

        .. math::

           p(z) = [\prod_n (z - x_n)] \sum_n \frac{u_n}{z - x_n}

        References
        ----------
        [1] R. M. Corless, Generalized companion matrices in the Lagrange basis,
            in Proceedings EACA, L. Gonzalez-Vega and T. Recio, eds., June 2004,
            pp. 317-322.
        """
        n = len(x)
        uv = u.reshape(-1, u.shape[-1])
        j = np.arange(1, n+1)
        z = np.empty((uv.shape[0], n-1), dtype=complex)

        for k in xrange(uv.shape[0]):
            A = np.zeros((n+1, n+1), dtype=dtype, order="F")
            A[0,1:] = -1
            A[1:,0] = u
            A[j,j] = x

            B = np.zeros((n+1, n+1), dtype=dtype, order="F")
            B[j,j] = 1

            with np.errstate(divide='ignore', invalid='ignore'):
                w = scipy.linalg.eig(A, B, right=False, overwrite_a=True, overwrite_b=True)

            # lexicographic sort puts infs to the end
            w[~np.isfinite(w)] = np.inf
            w.sort()
            z[k] = w[:-2]

        return z.reshape(u.shape[:-1] + (n-1,))

    def get_numerator(self, format='poly1d'):
        """
        Compute the numerator polynomial

        Parameters
        ----------
        format : {'poly1d', 'barycentric'}
            Format of returned polynomial, see below.

        Returns
        -------
        p
            Numerator polynomial

        Notes
        -----
        Format 'poly1d' polynomials are numerically unstable, and
        will return mostly noise at orders => 20.

        """
        if format == 'poly1d':
            return _lagpoly_poly1d(self.x, self.u * self.y / abs(self.u).max())
        elif format == 'barycentric':
            return _lagpoly_barycentric(self.x, self.u * self.y / abs(self.u).max())
        else:
            raise ValueError("Invalid polynomial format")

    def get_denominator(self, format='poly1d'):
        """
        Compute the denominator polynomial

        Parameters
        ----------
        format : {'poly1d', 'barycentric'}
            Format of returned polynomial, see below.

        Returns
        -------
        p
            Denominator polynomial

        Notes
        -----
        Format 'poly1d' polynomials are numerically unstable, and
        will return mostly noise at orders => 20.

        """
        if format == 'poly1d':
            return _lagpoly_poly1d(self.x, self.u / abs(self.u).max())
        elif format == 'barycentric':
            return _lagpoly_barycentric(self.x, self.u / abs(self.u).max())
        else:
            raise ValueError("Invalid polynomial format")

    @property
    def numerator_roots(self):
        """
        Roots of the numerator polynomial.

        Notes
        -----
        Note that if these roots coincide exactly with the roots of
        the denominator polynomial, the rational function may not
        have a root at that point.
        """
        if self._numerator_roots is None:
            self._numerator_roots = RationalBarycentricInterpolation._polyroots(
                self.x, self.y * self.u, self._dtype)
        return self._numerator_roots

    @property
    def denominator_roots(self):
        """
        Roots of the denominator polynomial.

        Notes
        -----
        Note that if these roots coincide exactly with the roots of
        the numerator polynomial, the rational function does not have
        a pole at that point.
        """
        if self._denominator_roots is None:
            self._denominator_roots = RationalBarycentricInterpolation._polyroots(
                self.x, self.u, self._dtype)
        return self._denominator_roots


def _lagpoly_poly1d(x, u):
    """
    poly1d power-basis representation of the polynomial

    .. math:: p(z) = [\prod_n (z - x_n)] \sum_n \frac{u_n}{z - x_n}
    """
    M = len(x)
    p = np.poly1d(0.0)
    for j in xrange(M):
        pt = np.poly1d(u[j])
        for k in xrange(M):
            if k == j:
                continue
            pt *= np.poly1d([1.0, -x[k]])
        p += pt
    return p


def _lagpoly_barycentric(x, u):
    """
    Barycentric representation of the polynomial

    .. math:: p(z) = [\prod_n (z - x_n)] \sum_n \frac{u_n}{z - x_n}
    """

    # XXX: this may underflow (accompanied by overflow inside BarycentricInterpolator)

    # w[i] = \prod_{j != i} (x[i] - x[j])
    w = np.zeros(len(x), dtype=np.find_common_type([x.dtype, np.inexact], []))
    w[0] = 1
    for j in xrange(1, len(x)):
        w[:j] *= (x[:j] - x[j])
        w[j] = np.multiply.reduce(x[j] - x[:j])

    return BarycentricInterpolator(x, u * w)


def floater_hormann(x, y, d=3):
    """
    Compute Floater-Hormann barycentric rational interpolant of data

    Parameters
    ----------
    x : array_like, shape (n,)
        Interpolation points
    y : array_like, shape (n,)
        Function values at interpolation points
    d : int, optional
        Floater-Hormann interpolant order

    Returns
    -------
    ip : RationalBarycentricInterpolation
        Rational barycentric interpolator of data

    Notes
    -----
    The Floater-Hormann interpolant is a rational function that
    interpolates the data with approximation order
    :math:`O(h^{d+1})`. The interpolant contains no poles on the real
    axis, unlike in typical rational interpolants.

    References
    ----------
    [1] M.S. Floater and K. Hormann, ''Barycentric rational interpolation
        with no poles and high rates of approximation'',
        Numer. Math. **107**, 315 (2007).
        http://dx.doi.org/10.1007/s00211-007-0093-y

    """

    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)

    if x.ndim != 1 or x.shape[0] != y.shape[-1]:
        raise ValueError("x and y do not have matching shape")

    d = min(x.size-1, d)
    if not d >= 0:
        raise ValueError("d not >= 0")

    x, y = _sort_and_average_duplicates(x, y)

    # Floater-Hormann weights
    u = np.zeros((n,), dtype=np.find_common_type([x.dtype, np.inexact], []))

    # Evaluate the start and end ranges with explicit for loops
    for i in itertools.chain(xrange(d), xrange(max(d, n - d), n)):
        for k in xrange(max(i-d, 0), min(i+1, n-d)):
            c = 1.0
            for p in xrange(k, k+d+1):
                if p != i:
                    c /= abs(x[i] - x[p])
            u[i] += c

    # Vectorize the evaluation for the center part
    if d < n-d:
        for k in xrange(-d, 1):
            c = 1.0
            for p in xrange(k, k+d+1):
                if p != 0:
                    c /= abs(x[d:n-d] - x[d+p:n+p-d])
            u[d:n-d] += c

    u[((d+1)%2)::2] *= -1

    return RationalBarycentricInterpolation(x, y, u)


def lsq_rational(x, y, m=None, tol=1e-14):
    """
    Compute least-squares classical barycentric rational interpolant of data

    Parameters
    ----------
    x : array_like, shape (N,)
        Interpolation points
    y : array_like, shape (..., N)
        Data at interpolation points
    m : int, optional
        Degree of the numerator polynomial, must be ``N-1 >= m >= (N-1)/2``.
        Default: max(N//2, (N-1)//2)

    Returns
    -------
    ip : RationalBarycentricInterpolation
        Rational barycentric interpolator of data

    Warns
    -----
    RationalInterpolationWarning
        Emitted if the computed interpolant does not pass through data
        points, to relative accuracy specified by `tol`.

    Notes
    -----
    The fitting does not constrain locations of the poles of the rational
    interpolant. The following caveats apply, as usual in classical rational
    interpolation:

    - The interpolant may contain divergences

    - The constraints on ``m`` may mean it is impossible for the
      interpolant to pass through all data points. If this occurs,
      this routine emits a warning.

    References
    ----------
    [1] J.-P. Berrut, R. Baltensperger, H.D. Mittelmann, ''Recent
        developments in barycentric rational interpolation'',
        ISNM International Series of Numerical Mathematics **151**, 27
        (2005). http://dx.doi.org/10.1007/3-7643-7356-3_3

    """

    x = np.asarray(x)
    y = np.asarray(y)
    N = len(x)

    if x.ndim != 1 or N != y.shape[-1]:
        raise ValueError("x and y do not have matching shape")

    if N <= 1:
        return RationalBarycentricInterpolation(x, y, np.array([1.0]))

    if m is None:
        m = max(N//2, (N-1)//2)

    if 2*m < N - 1:
        raise ValueError("m must be >= (N-1)/2")
    if m > N - 1:
        raise ValueError("m must be < N - 1")

    dtype = np.find_common_type([x.dtype, y.dtype, np.inexact], [])

    # Shift and scale the fitting problem
    x0 = x
    x = x - x.mean()
    x_scale = x.ptp()
    if x_scale != 0:
        x /= x_scale

    y0 = y
    vy = y.reshape(-1, N)
    vy = vy - vy.mean()
    y_scale = vy.ptp(axis=-1)
    y_scale[y_scale==0] = 1
    vy /= y_scale

    # Solve interpolation conditions
    u = np.empty((vy.shape[0], N), dtype=dtype)
    for k in xrange(vy.shape[0]):
        V = np.vander(x, m).T[::-1]
        F = vy[k][None,:] * np.vander(x, N - m - 1).T[::-1]
        A = np.concatenate((V, F), axis=0)

        # Appropriate weights are in the null space of A
        _, s, vh = scipy.linalg.svd(A, full_matrices=True, overwrite_a=True)
        u[k,:] = vh.T[:,-1].conj()

    u = u.reshape(y.shape)

    # Check interpolation
    if abs(u).min() <= tol*abs(u).max():
        warnings.warn("Rational interpolant may not fit all points",
                      RationalInterpolationWarning)

    return RationalBarycentricInterpolation(x0, y0, u)


def _sort_and_average_duplicates(x, y):
    x = np.atleast_1d(x)
    y = np.asarray(y)

    # Sort
    j = np.argsort(x)
    x = x[j]
    y = y[...,j]

    # Average duplicates
    m = np.flatnonzero((x[1:] != x[:-1]))
    if m.size < x.size:
        j = np.r_[0, 1 + m]
        y = np.add.reduceat(y, j, axis=-1) / np.diff(np.r_[j, x.size])
        x = x[j]

    return x, y
