# -*- coding: utf-8 -*-
#
# Author:  Travis Oliphant  2002-2011 with contributions from
#          SciPy Developers 2004-2011
#
import warnings
from collections.abc import Iterable
from functools import wraps, cached_property
import ctypes

import numpy as np
from numpy.polynomial import Polynomial
from scipy._lib.doccer import (extend_notes_in_docstring,
                               replace_notes_in_docstring,
                               inherit_docstring_from)
from scipy._lib._ccallback import LowLevelCallable
from scipy import optimize
from scipy import integrate
import scipy.special as sc

import scipy.special._ufuncs as scu
from scipy._lib._util import _lazyselect, _lazywhere
from . import _stats
from ._tukeylambda_stats import (tukeylambda_variance as _tlvar,
                                 tukeylambda_kurtosis as _tlkurt)
from ._distn_infrastructure import (
    get_distribution_names, _kurtosis,
    rv_continuous, _skew, _get_fixed_fit_value, _check_shape, _ShapeInfo)
from ._ksstats import kolmogn, kolmognp, kolmogni
from ._constants import (_XMIN, _EULER, _ZETA3, _SQRT_PI,
                         _SQRT_2_OVER_PI, _LOG_SQRT_2_OVER_PI)
import scipy.stats._boost as _boost
from scipy.optimize import root_scalar
from scipy.stats._warnings_errors import FitError
import scipy.stats as stats


def _remove_optimizer_parameters(kwds):
    """
    Remove the optimizer-related keyword arguments 'loc', 'scale' and
    'optimizer' from `kwds`.  Then check that `kwds` is empty, and
    raise `TypeError("Unknown arguments: %s." % kwds)` if it is not.

    This function is used in the fit method of distributions that override
    the default method and do not use the default optimization code.

    `kwds` is modified in-place.
    """
    kwds.pop('loc', None)
    kwds.pop('scale', None)
    kwds.pop('optimizer', None)
    kwds.pop('method', None)
    if kwds:
        raise TypeError("Unknown arguments: %s." % kwds)


def _call_super_mom(fun):
    # if fit method is overridden only for MLE and doesn't specify what to do
    # if method == 'mm', this decorator calls generic implementation
    @wraps(fun)
    def wrapper(self, *args, **kwds):
        method = kwds.get('method', 'mle').lower()
        if method == 'mm':
            return super(type(self), self).fit(*args, **kwds)
        else:
            return fun(self, *args, **kwds)
    return wrapper


class ksone_gen(rv_continuous):
    r"""Kolmogorov-Smirnov one-sided test statistic distribution.

    This is the distribution of the one-sided Kolmogorov-Smirnov (KS)
    statistics :math:`D_n^+` and :math:`D_n^-`
    for a finite sample size ``n >= 1`` (the shape parameter).

    %(before_notes)s

    See Also
    --------
    kstwobign, kstwo, kstest

    Notes
    -----
    :math:`D_n^+` and :math:`D_n^-` are given by

    .. math::

        D_n^+ &= \text{sup}_x (F_n(x) - F(x)),\\
        D_n^- &= \text{sup}_x (F(x) - F_n(x)),\\

    where :math:`F` is a continuous CDF and :math:`F_n` is an empirical CDF.
    `ksone` describes the distribution under the null hypothesis of the KS test
    that the empirical CDF corresponds to :math:`n` i.i.d. random variates
    with CDF :math:`F`.

    %(after_notes)s

    References
    ----------
    .. [1] Birnbaum, Z. W. and Tingey, F.H. "One-sided confidence contours
       for probability distribution functions", The Annals of Mathematical
       Statistics, 22(4), pp 592-596 (1951).

    %(example)s

    """
    def _argcheck(self, n):
        return (n >= 1) & (n == np.round(n))

    def _shape_info(self):
        return [_ShapeInfo("n", True, (1, np.inf), (True, False))]

    def _pdf(self, x, n):
        return -scu._smirnovp(n, x)

    def _cdf(self, x, n):
        return scu._smirnovc(n, x)

    def _sf(self, x, n):
        return sc.smirnov(n, x)

    def _ppf(self, q, n):
        return scu._smirnovci(n, q)

    def _isf(self, q, n):
        return sc.smirnovi(n, q)


ksone = ksone_gen(a=0.0, b=1.0, name='ksone')


class kstwo_gen(rv_continuous):
    r"""Kolmogorov-Smirnov two-sided test statistic distribution.

    This is the distribution of the two-sided Kolmogorov-Smirnov (KS)
    statistic :math:`D_n` for a finite sample size ``n >= 1``
    (the shape parameter).

    %(before_notes)s

    See Also
    --------
    kstwobign, ksone, kstest

    Notes
    -----
    :math:`D_n` is given by

    .. math::

        D_n = \text{sup}_x |F_n(x) - F(x)|

    where :math:`F` is a (continuous) CDF and :math:`F_n` is an empirical CDF.
    `kstwo` describes the distribution under the null hypothesis of the KS test
    that the empirical CDF corresponds to :math:`n` i.i.d. random variates
    with CDF :math:`F`.

    %(after_notes)s

    References
    ----------
    .. [1] Simard, R., L'Ecuyer, P. "Computing the Two-Sided
       Kolmogorov-Smirnov Distribution",  Journal of Statistical Software,
       Vol 39, 11, 1-18 (2011).

    %(example)s

    """
    def _argcheck(self, n):
        return (n >= 1) & (n == np.round(n))

    def _shape_info(self):
        return [_ShapeInfo("n", True, (1, np.inf), (True, False))]

    def _get_support(self, n):
        return (0.5/(n if not isinstance(n, Iterable) else np.asanyarray(n)),
                1.0)

    def _pdf(self, x, n):
        return kolmognp(n, x)

    def _cdf(self, x, n):
        return kolmogn(n, x)

    def _sf(self, x, n):
        return kolmogn(n, x, cdf=False)

    def _ppf(self, q, n):
        return kolmogni(n, q, cdf=True)

    def _isf(self, q, n):
        return kolmogni(n, q, cdf=False)


# Use the pdf, (not the ppf) to compute moments
kstwo = kstwo_gen(momtype=0, a=0.0, b=1.0, name='kstwo')


class kstwobign_gen(rv_continuous):
    r"""Limiting distribution of scaled Kolmogorov-Smirnov two-sided test statistic.

    This is the asymptotic distribution of the two-sided Kolmogorov-Smirnov
    statistic :math:`\sqrt{n} D_n` that measures the maximum absolute
    distance of the theoretical (continuous) CDF from the empirical CDF.
    (see `kstest`).

    %(before_notes)s

    See Also
    --------
    ksone, kstwo, kstest

    Notes
    -----
    :math:`\sqrt{n} D_n` is given by

    .. math::

        D_n = \text{sup}_x |F_n(x) - F(x)|

    where :math:`F` is a continuous CDF and :math:`F_n` is an empirical CDF.
    `kstwobign`  describes the asymptotic distribution (i.e. the limit of
    :math:`\sqrt{n} D_n`) under the null hypothesis of the KS test that the
    empirical CDF corresponds to i.i.d. random variates with CDF :math:`F`.

    %(after_notes)s

    References
    ----------
    .. [1] Feller, W. "On the Kolmogorov-Smirnov Limit Theorems for Empirical
       Distributions",  Ann. Math. Statist. Vol 19, 177-189 (1948).

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        return -scu._kolmogp(x)

    def _cdf(self, x):
        return scu._kolmogc(x)

    def _sf(self, x):
        return sc.kolmogorov(x)

    def _ppf(self, q):
        return scu._kolmogci(q)

    def _isf(self, q):
        return sc.kolmogi(q)


kstwobign = kstwobign_gen(a=0.0, name='kstwobign')


## Normal distribution

# loc = mu, scale = std
# Keep these implementations out of the class definition so they can be reused
# by other distributions.
_norm_pdf_C = np.sqrt(2*np.pi)
_norm_pdf_logC = np.log(_norm_pdf_C)


def _norm_pdf(x):
    return np.exp(-x**2/2.0) / _norm_pdf_C


def _norm_logpdf(x):
    return -x**2 / 2.0 - _norm_pdf_logC


def _norm_cdf(x):
    return sc.ndtr(x)


def _norm_logcdf(x):
    return sc.log_ndtr(x)


def _norm_ppf(q):
    return sc.ndtri(q)


def _norm_sf(x):
    return _norm_cdf(-x)


def _norm_logsf(x):
    return _norm_logcdf(-x)


def _norm_isf(q):
    return -_norm_ppf(q)


class norm_gen(rv_continuous):
    r"""A normal continuous random variable.

    The location (``loc``) keyword specifies the mean.
    The scale (``scale``) keyword specifies the standard deviation.

    %(before_notes)s

    Notes
    -----
    The probability density function for `norm` is:

    .. math::

        f(x) = \frac{\exp(-x^2/2)}{\sqrt{2\pi}}

    for a real number :math:`x`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return random_state.standard_normal(size)

    def _pdf(self, x):
        # norm.pdf(x) = exp(-x**2/2)/sqrt(2*pi)
        return _norm_pdf(x)

    def _logpdf(self, x):
        return _norm_logpdf(x)

    def _cdf(self, x):
        return _norm_cdf(x)

    def _logcdf(self, x):
        return _norm_logcdf(x)

    def _sf(self, x):
        return _norm_sf(x)

    def _logsf(self, x):
        return _norm_logsf(x)

    def _ppf(self, q):
        return _norm_ppf(q)

    def _isf(self, q):
        return _norm_isf(q)

    def _stats(self):
        return 0.0, 1.0, 0.0, 0.0

    def _entropy(self):
        return 0.5*(np.log(2*np.pi)+1)

    @_call_super_mom
    @replace_notes_in_docstring(rv_continuous, notes="""\
        For the normal distribution, method of moments and maximum likelihood
        estimation give identical fits, and explicit formulas for the estimates
        are available.
        This function uses these explicit formulas for the maximum likelihood
        estimation of the normal distribution parameters, so the
        `optimizer` and `method` arguments are ignored.\n\n""")
    def fit(self, data, **kwds):

        floc = kwds.pop('floc', None)
        fscale = kwds.pop('fscale', None)

        _remove_optimizer_parameters(kwds)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            # Without this check, this function would just return the
            # parameters that were given.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        if not np.isfinite(data).all():
            raise ValueError("The data contains non-finite values.")

        if floc is None:
            loc = data.mean()
        else:
            loc = floc

        if fscale is None:
            scale = np.sqrt(((data - loc)**2).mean())
        else:
            scale = fscale

        return loc, scale

    def _munp(self, n):
        """
        @returns Moments of standard normal distribution for integer n >= 0

        See eq. 16 of https://arxiv.org/abs/1209.4340v2
        """
        if n % 2 == 0:
            return sc.factorial2(n - 1)
        else:
            return 0.


norm = norm_gen(name='norm')


class alpha_gen(rv_continuous):
    r"""An alpha continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `alpha` ([1]_, [2]_) is:

    .. math::

        f(x, a) = \frac{1}{x^2 \Phi(a) \sqrt{2\pi}} *
                  \exp(-\frac{1}{2} (a-1/x)^2)

    where :math:`\Phi` is the normal CDF, :math:`x > 0`, and :math:`a > 0`.

    `alpha` takes ``a`` as a shape parameter.

    %(after_notes)s

    References
    ----------
    .. [1] Johnson, Kotz, and Balakrishnan, "Continuous Univariate
           Distributions, Volume 1", Second Edition, John Wiley and Sons,
           p. 173 (1994).
    .. [2] Anthony A. Salvia, "Reliability applications of the Alpha
           Distribution", IEEE Transactions on Reliability, Vol. R-34,
           No. 3, pp. 251-252 (1985).

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return [_ShapeInfo("a", False, (0, np.inf), (False, False))]

    def _pdf(self, x, a):
        # alpha.pdf(x, a) = 1/(x**2*Phi(a)*sqrt(2*pi)) * exp(-1/2 * (a-1/x)**2)
        return 1.0/(x**2)/_norm_cdf(a)*_norm_pdf(a-1.0/x)

    def _logpdf(self, x, a):
        return -2*np.log(x) + _norm_logpdf(a-1.0/x) - np.log(_norm_cdf(a))

    def _cdf(self, x, a):
        return _norm_cdf(a-1.0/x) / _norm_cdf(a)

    def _ppf(self, q, a):
        return 1.0/np.asarray(a-sc.ndtri(q*_norm_cdf(a)))

    def _stats(self, a):
        return [np.inf]*2 + [np.nan]*2


alpha = alpha_gen(a=0.0, name='alpha')


class anglit_gen(rv_continuous):
    r"""An anglit continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `anglit` is:

    .. math::

        f(x) = \sin(2x + \pi/2) = \cos(2x)

    for :math:`-\pi/4 \le x \le \pi/4`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # anglit.pdf(x) = sin(2*x + \pi/2) = cos(2*x)
        return np.cos(2*x)

    def _cdf(self, x):
        return np.sin(x+np.pi/4)**2.0

    def _ppf(self, q):
        return np.arcsin(np.sqrt(q))-np.pi/4

    def _stats(self):
        return 0.0, np.pi*np.pi/16-0.5, 0.0, -2*(np.pi**4 - 96)/(np.pi*np.pi-8)**2

    def _entropy(self):
        return 1-np.log(2)


anglit = anglit_gen(a=-np.pi/4, b=np.pi/4, name='anglit')


class arcsine_gen(rv_continuous):
    r"""An arcsine continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `arcsine` is:

    .. math::

        f(x) = \frac{1}{\pi \sqrt{x (1-x)}}

    for :math:`0 < x < 1`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # arcsine.pdf(x) = 1/(pi*sqrt(x*(1-x)))
        with np.errstate(divide='ignore'):
            return 1.0/np.pi/np.sqrt(x*(1-x))

    def _cdf(self, x):
        return 2.0/np.pi*np.arcsin(np.sqrt(x))

    def _ppf(self, q):
        return np.sin(np.pi/2.0*q)**2.0

    def _stats(self):
        mu = 0.5
        mu2 = 1.0/8
        g1 = 0
        g2 = -3.0/2.0
        return mu, mu2, g1, g2

    def _entropy(self):
        return -0.24156447527049044468


arcsine = arcsine_gen(a=0.0, b=1.0, name='arcsine')


class FitDataError(ValueError):
    """Raised when input data is inconsistent with fixed parameters."""
    # This exception is raised by, for example, beta_gen.fit when both floc
    # and fscale are fixed and there are values in the data not in the open
    # interval (floc, floc+fscale).
    def __init__(self, distr, lower, upper):
        self.args = (
            "Invalid values in `data`.  Maximum likelihood "
            "estimation with {distr!r} requires that {lower!r} < "
            "(x - loc)/scale  < {upper!r} for each x in `data`.".format(
                distr=distr, lower=lower, upper=upper),
        )


class FitSolverError(FitError):
    """
    Raised when a solver fails to converge while fitting a distribution.
    """
    # This exception is raised by, for example, beta_gen.fit when
    # optimize.fsolve returns with ier != 1.
    def __init__(self, mesg):
        emsg = "Solver for the MLE equations failed to converge: "
        emsg += mesg.replace('\n', '')
        self.args = (emsg,)


def _beta_mle_a(a, b, n, s1):
    # The zeros of this function give the MLE for `a`, with
    # `b`, `n` and `s1` given.  `s1` is the sum of the logs of
    # the data. `n` is the number of data points.
    psiab = sc.psi(a + b)
    func = s1 - n * (-psiab + sc.psi(a))
    return func


def _beta_mle_ab(theta, n, s1, s2):
    # Zeros of this function are critical points of
    # the maximum likelihood function.  Solving this system
    # for theta (which contains a and b) gives the MLE for a and b
    # given `n`, `s1` and `s2`.  `s1` is the sum of the logs of the data,
    # and `s2` is the sum of the logs of 1 - data.  `n` is the number
    # of data points.
    a, b = theta
    psiab = sc.psi(a + b)
    func = [s1 - n * (-psiab + sc.psi(a)),
            s2 - n * (-psiab + sc.psi(b))]
    return func


class beta_gen(rv_continuous):
    r"""A beta continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `beta` is:

    .. math::

        f(x, a, b) = \frac{\Gamma(a+b) x^{a-1} (1-x)^{b-1}}
                          {\Gamma(a) \Gamma(b)}

    for :math:`0 <= x <= 1`, :math:`a > 0`, :math:`b > 0`, where
    :math:`\Gamma` is the gamma function (`scipy.special.gamma`).

    `beta` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ia, ib]

    def _rvs(self, a, b, size=None, random_state=None):
        return random_state.beta(a, b, size)

    def _pdf(self, x, a, b):
        #                     gamma(a+b) * x**(a-1) * (1-x)**(b-1)
        # beta.pdf(x, a, b) = ------------------------------------
        #                              gamma(a)*gamma(b)
        return _boost._beta_pdf(x, a, b)

    def _logpdf(self, x, a, b):
        lPx = sc.xlog1py(b - 1.0, -x) + sc.xlogy(a - 1.0, x)
        lPx -= sc.betaln(a, b)
        return lPx

    def _cdf(self, x, a, b):
        return _boost._beta_cdf(x, a, b)

    def _sf(self, x, a, b):
        return _boost._beta_sf(x, a, b)

    def _isf(self, x, a, b):
        with warnings.catch_warnings():
            # See gh-14901
            message = "overflow encountered in _beta_isf"
            warnings.filterwarnings('ignore', message=message)
            return _boost._beta_isf(x, a, b)

    def _ppf(self, q, a, b):
        with warnings.catch_warnings():
            message = "overflow encountered in _beta_ppf"
            warnings.filterwarnings('ignore', message=message)
            return _boost._beta_ppf(q, a, b)

    def _stats(self, a, b):
        return(
            _boost._beta_mean(a, b),
            _boost._beta_variance(a, b),
            _boost._beta_skewness(a, b),
            _boost._beta_kurtosis_excess(a, b))

    def _fitstart(self, data):
        g1 = _skew(data)
        g2 = _kurtosis(data)

        def func(x):
            a, b = x
            sk = 2*(b-a)*np.sqrt(a + b + 1) / (a + b + 2) / np.sqrt(a*b)
            ku = a**3 - a**2*(2*b-1) + b**2*(b+1) - 2*a*b*(b+2)
            ku /= a*b*(a+b+2)*(a+b+3)
            ku *= 6
            return [sk-g1, ku-g2]
        a, b = optimize.fsolve(func, (1.0, 1.0))
        return super()._fitstart(data, args=(a, b))

    @_call_super_mom
    @extend_notes_in_docstring(rv_continuous, notes="""\
        In the special case where `method="MLE"` and
        both `floc` and `fscale` are given, a
        `ValueError` is raised if any value `x` in `data` does not satisfy
        `floc < x < floc + fscale`.\n\n""")
    def fit(self, data, *args, **kwds):
        # Override rv_continuous.fit, so we can more efficiently handle the
        # case where floc and fscale are given.

        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is None or fscale is None:
            # do general fit
            return super().fit(data, *args, **kwds)

        # We already got these from kwds, so just pop them.
        kwds.pop('floc', None)
        kwds.pop('fscale', None)

        f0 = _get_fixed_fit_value(kwds, ['f0', 'fa', 'fix_a'])
        f1 = _get_fixed_fit_value(kwds, ['f1', 'fb', 'fix_b'])

        _remove_optimizer_parameters(kwds)

        if f0 is not None and f1 is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        # Special case: loc and scale are constrained, so we are fitting
        # just the shape parameters.  This can be done much more efficiently
        # than the method used in `rv_continuous.fit`.  (See the subsection
        # "Two unknown parameters" in the section "Maximum likelihood" of
        # the Wikipedia article on the Beta distribution for the formulas.)

        if not np.isfinite(data).all():
            raise ValueError("The data contains non-finite values.")

        # Normalize the data to the interval [0, 1].
        data = (np.ravel(data) - floc) / fscale
        if np.any(data <= 0) or np.any(data >= 1):
            raise FitDataError("beta", lower=floc, upper=floc + fscale)

        xbar = data.mean()

        if f0 is not None or f1 is not None:
            # One of the shape parameters is fixed.

            if f0 is not None:
                # The shape parameter a is fixed, so swap the parameters
                # and flip the data.  We always solve for `a`.  The result
                # will be swapped back before returning.
                b = f0
                data = 1 - data
                xbar = 1 - xbar
            else:
                b = f1

            # Initial guess for a.  Use the formula for the mean of the beta
            # distribution, E[x] = a / (a + b), to generate a reasonable
            # starting point based on the mean of the data and the given
            # value of b.
            a = b * xbar / (1 - xbar)

            # Compute the MLE for `a` by solving _beta_mle_a.
            theta, info, ier, mesg = optimize.fsolve(
                _beta_mle_a, a,
                args=(b, len(data), np.log(data).sum()),
                full_output=True
            )
            if ier != 1:
                raise FitSolverError(mesg=mesg)
            a = theta[0]

            if f0 is not None:
                # The shape parameter a was fixed, so swap back the
                # parameters.
                a, b = b, a

        else:
            # Neither of the shape parameters is fixed.

            # s1 and s2 are used in the extra arguments passed to _beta_mle_ab
            # by optimize.fsolve.
            s1 = np.log(data).sum()
            s2 = sc.log1p(-data).sum()

            # Use the "method of moments" to estimate the initial
            # guess for a and b.
            fac = xbar * (1 - xbar) / data.var(ddof=0) - 1
            a = xbar * fac
            b = (1 - xbar) * fac

            # Compute the MLE for a and b by solving _beta_mle_ab.
            theta, info, ier, mesg = optimize.fsolve(
                _beta_mle_ab, [a, b],
                args=(len(data), s1, s2),
                full_output=True
            )
            if ier != 1:
                raise FitSolverError(mesg=mesg)
            a, b = theta

        return a, b, floc, fscale


beta = beta_gen(a=0.0, b=1.0, name='beta')


class betaprime_gen(rv_continuous):
    r"""A beta prime continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `betaprime` is:

    .. math::

        f(x, a, b) = \frac{x^{a-1} (1+x)^{-a-b}}{\beta(a, b)}

    for :math:`x >= 0`, :math:`a > 0`, :math:`b > 0`, where
    :math:`\beta(a, b)` is the beta function (see `scipy.special.beta`).

    `betaprime` takes ``a`` and ``b`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ia, ib]

    def _rvs(self, a, b, size=None, random_state=None):
        u1 = gamma.rvs(a, size=size, random_state=random_state)
        u2 = gamma.rvs(b, size=size, random_state=random_state)
        return u1 / u2

    def _pdf(self, x, a, b):
        # betaprime.pdf(x, a, b) = x**(a-1) * (1+x)**(-a-b) / beta(a, b)
        return np.exp(self._logpdf(x, a, b))

    def _logpdf(self, x, a, b):
        return sc.xlogy(a - 1.0, x) - sc.xlog1py(a + b, x) - sc.betaln(a, b)

    def _cdf(self, x, a, b):
        return sc.betainc(a, b, x/(1.+x))

    def _munp(self, n, a, b):
        if n == 1.0:
            return np.where(b > 1,
                            a/(b-1.0),
                            np.inf)
        elif n == 2.0:
            return np.where(b > 2,
                            a*(a+1.0)/((b-2.0)*(b-1.0)),
                            np.inf)
        elif n == 3.0:
            return np.where(b > 3,
                            a*(a+1.0)*(a+2.0)/((b-3.0)*(b-2.0)*(b-1.0)),
                            np.inf)
        elif n == 4.0:
            return np.where(b > 4,
                            (a*(a + 1.0)*(a + 2.0)*(a + 3.0) /
                             ((b - 4.0)*(b - 3.0)*(b - 2.0)*(b - 1.0))),
                            np.inf)
        else:
            raise NotImplementedError


betaprime = betaprime_gen(a=0.0, name='betaprime')


class bradford_gen(rv_continuous):
    r"""A Bradford continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `bradford` is:

    .. math::

        f(x, c) = \frac{c}{\log(1+c) (1+cx)}

    for :math:`0 <= x <= 1` and :math:`c > 0`.

    `bradford` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # bradford.pdf(x, c) = c / (k * (1+c*x))
        return c / (c*x + 1.0) / sc.log1p(c)

    def _cdf(self, x, c):
        return sc.log1p(c*x) / sc.log1p(c)

    def _ppf(self, q, c):
        return sc.expm1(q * sc.log1p(c)) / c

    def _stats(self, c, moments='mv'):
        k = np.log(1.0+c)
        mu = (c-k)/(c*k)
        mu2 = ((c+2.0)*k-2.0*c)/(2*c*k*k)
        g1 = None
        g2 = None
        if 's' in moments:
            g1 = np.sqrt(2)*(12*c*c-9*c*k*(c+2)+2*k*k*(c*(c+3)+3))
            g1 /= np.sqrt(c*(c*(k-2)+2*k))*(3*c*(k-2)+6*k)
        if 'k' in moments:
            g2 = (c**3*(k-3)*(k*(3*k-16)+24)+12*k*c*c*(k-4)*(k-3) +
                  6*c*k*k*(3*k-14) + 12*k**3)
            g2 /= 3*c*(c*(k-2)+2*k)**2
        return mu, mu2, g1, g2

    def _entropy(self, c):
        k = np.log(1+c)
        return k/2.0 - np.log(c/k)


bradford = bradford_gen(a=0.0, b=1.0, name='bradford')


class burr_gen(rv_continuous):
    r"""A Burr (Type III) continuous random variable.

    %(before_notes)s

    See Also
    --------
    fisk : a special case of either `burr` or `burr12` with ``d=1``
    burr12 : Burr Type XII distribution
    mielke : Mielke Beta-Kappa / Dagum distribution

    Notes
    -----
    The probability density function for `burr` is:

    .. math::

        f(x; c, d) = c d \frac{x^{-c - 1}}
                              {{(1 + x^{-c})}^{d + 1}}

    for :math:`x >= 0` and :math:`c, d > 0`.

    `burr` takes ``c`` and ``d`` as shape parameters for :math:`c` and
    :math:`d`.

    This is the PDF corresponding to the third CDF given in Burr's list;
    specifically, it is equation (11) in Burr's paper [1]_. The distribution
    is also commonly referred to as the Dagum distribution [2]_. If the
    parameter :math:`c < 1` then the mean of the distribution does not
    exist and if :math:`c < 2` the variance does not exist [2]_.
    The PDF is finite at the left endpoint :math:`x = 0` if :math:`c * d >= 1`.

    %(after_notes)s

    References
    ----------
    .. [1] Burr, I. W. "Cumulative frequency functions", Annals of
       Mathematical Statistics, 13(2), pp 215-232 (1942).
    .. [2] https://en.wikipedia.org/wiki/Dagum_distribution
    .. [3] Kleiber, Christian. "A guide to the Dagum distributions."
       Modeling Income Distributions and Lorenz Curves  pp 97-117 (2008).

    %(example)s

    """
    # Do not set _support_mask to rv_continuous._open_support_mask
    # Whether the left-hand endpoint is suitable for pdf evaluation is dependent
    # on the values of c and d: if c*d >= 1, the pdf is finite, otherwise infinite.

    def _shape_info(self):
        ic = _ShapeInfo("c", False, (0, np.inf), (False, False))
        id = _ShapeInfo("d", False, (0, np.inf), (False, False))
        return [ic, id]

    def _pdf(self, x, c, d):
        # burr.pdf(x, c, d) = c * d * x**(-c-1) * (1+x**(-c))**(-d-1)
        output = _lazywhere(x == 0, [x, c, d],
                   lambda x_, c_, d_: c_ * d_ * (x_**(c_*d_-1)) / (1 + x_**c_),
                   f2 = lambda x_, c_, d_: (c_ * d_ * (x_ ** (-c_ - 1.0)) /
                                            ((1 + x_ ** (-c_)) ** (d_ + 1.0))))
        if output.ndim == 0:
            return output[()]
        return output

    def _logpdf(self, x, c, d):
        output = _lazywhere(
            x == 0, [x, c, d],
            lambda x_, c_, d_: (np.log(c_) + np.log(d_) + sc.xlogy(c_*d_ - 1, x_)
                                - (d_+1) * sc.log1p(x_**(c_))),
            f2 = lambda x_, c_, d_: (np.log(c_) + np.log(d_)
                                     + sc.xlogy(-c_ - 1, x_)
                                     - sc.xlog1py(d_+1, x_**(-c_))))
        if output.ndim == 0:
            return output[()]
        return output

    def _cdf(self, x, c, d):
        return (1 + x**(-c))**(-d)

    def _logcdf(self, x, c, d):
        return sc.log1p(x**(-c)) * (-d)

    def _sf(self, x, c, d):
        return np.exp(self._logsf(x, c, d))

    def _logsf(self, x, c, d):
        return np.log1p(- (1 + x**(-c))**(-d))

    def _ppf(self, q, c, d):
        return (q**(-1.0/d) - 1)**(-1.0/c)

    def _stats(self, c, d):
        nc = np.arange(1, 5).reshape(4,1) / c
        #ek is the kth raw moment, e1 is the mean e2-e1**2 variance etc.
        e1, e2, e3, e4 = sc.beta(d + nc, 1. - nc) * d
        mu = np.where(c > 1.0, e1, np.nan)
        mu2_if_c = e2 - mu**2
        mu2 = np.where(c > 2.0, mu2_if_c, np.nan)
        g1 = _lazywhere(
            c > 3.0,
            (c, e1, e2, e3, mu2_if_c),
            lambda c, e1, e2, e3, mu2_if_c: (e3 - 3*e2*e1 + 2*e1**3) / np.sqrt((mu2_if_c)**3),
            fillvalue=np.nan)
        g2 = _lazywhere(
            c > 4.0,
            (c, e1, e2, e3, e4, mu2_if_c),
            lambda c, e1, e2, e3, e4, mu2_if_c: (
                ((e4 - 4*e3*e1 + 6*e2*e1**2 - 3*e1**4) / mu2_if_c**2) - 3),
            fillvalue=np.nan)
        if np.ndim(c) == 0:
            return mu.item(), mu2.item(), g1.item(), g2.item()
        return mu, mu2, g1, g2

    def _munp(self, n, c, d):
        def __munp(n, c, d):
            nc = 1. * n / c
            return d * sc.beta(1.0 - nc, d + nc)
        n, c, d = np.asarray(n), np.asarray(c), np.asarray(d)
        return _lazywhere((c > n) & (n == n) & (d == d), (c, d, n),
                          lambda c, d, n: __munp(n, c, d),
                          np.nan)


burr = burr_gen(a=0.0, name='burr')


class burr12_gen(rv_continuous):
    r"""A Burr (Type XII) continuous random variable.

    %(before_notes)s

    See Also
    --------
    fisk : a special case of either `burr` or `burr12` with ``d=1``
    burr : Burr Type III distribution

    Notes
    -----
    The probability density function for `burr12` is:

    .. math::

        f(x; c, d) = c d \frac{x^{c-1}}
                              {(1 + x^c)^{d + 1}}

    for :math:`x >= 0` and :math:`c, d > 0`.

    `burr12` takes ``c`` and ``d`` as shape parameters for :math:`c`
    and :math:`d`.

    This is the PDF corresponding to the twelfth CDF given in Burr's list;
    specifically, it is equation (20) in Burr's paper [1]_.

    %(after_notes)s

    The Burr type 12 distribution is also sometimes referred to as
    the Singh-Maddala distribution from NIST [2]_.

    References
    ----------
    .. [1] Burr, I. W. "Cumulative frequency functions", Annals of
       Mathematical Statistics, 13(2), pp 215-232 (1942).

    .. [2] https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/b12pdf.htm

    .. [3] "Burr distribution",
       https://en.wikipedia.org/wiki/Burr_distribution

    %(example)s

    """
    def _shape_info(self):
        ic = _ShapeInfo("c", False, (0, np.inf), (False, False))
        id = _ShapeInfo("d", False, (0, np.inf), (False, False))
        return [ic, id]

    def _pdf(self, x, c, d):
        # burr12.pdf(x, c, d) = c * d * x**(c-1) * (1+x**(c))**(-d-1)
        return np.exp(self._logpdf(x, c, d))

    def _logpdf(self, x, c, d):
        return np.log(c) + np.log(d) + sc.xlogy(c - 1, x) + sc.xlog1py(-d-1, x**c)

    def _cdf(self, x, c, d):
        return -sc.expm1(self._logsf(x, c, d))

    def _logcdf(self, x, c, d):
        return sc.log1p(-(1 + x**c)**(-d))

    def _sf(self, x, c, d):
        return np.exp(self._logsf(x, c, d))

    def _logsf(self, x, c, d):
        return sc.xlog1py(-d, x**c)

    def _ppf(self, q, c, d):
        # The following is an implementation of
        #   ((1 - q)**(-1.0/d) - 1)**(1.0/c)
        # that does a better job handling small values of q.
        return sc.expm1(-1/d * sc.log1p(-q))**(1/c)

    def _munp(self, n, c, d):
        nc = 1. * n / c
        return d * sc.beta(1.0 + nc, d - nc)


burr12 = burr12_gen(a=0.0, name='burr12')


class fisk_gen(burr_gen):
    r"""A Fisk continuous random variable.

    The Fisk distribution is also known as the log-logistic distribution.

    %(before_notes)s

    See Also
    --------
    burr

    Notes
    -----
    The probability density function for `fisk` is:

    .. math::

        f(x, c) = \frac{c x^{c-1}}
                       {(1 + x^c)^2}

    for :math:`x >= 0` and :math:`c > 0`.

    Please note that the above expression can be transformed into the following
    one, which is also commonly used:

    .. math::

        f(x, c) = \frac{c x^{-c-1}}
                       {(1 + x^{-c})^2}

    `fisk` takes ``c`` as a shape parameter for :math:`c`.

    `fisk` is a special case of `burr` or `burr12` with ``d=1``.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # fisk.pdf(x, c) = c * x**(-c-1) * (1 + x**(-c))**(-2)
        return burr._pdf(x, c, 1.0)

    def _cdf(self, x, c):
        return burr._cdf(x, c, 1.0)

    def _sf(self, x, c):
        return burr._sf(x, c, 1.0)

    def _logpdf(self, x, c):
        # fisk.pdf(x, c) = c * x**(-c-1) * (1 + x**(-c))**(-2)
        return burr._logpdf(x, c, 1.0)

    def _logcdf(self, x, c):
        return burr._logcdf(x, c, 1.0)

    def _logsf(self, x, c):
        return burr._logsf(x, c, 1.0)

    def _ppf(self, x, c):
        return burr._ppf(x, c, 1.0)

    def _munp(self, n, c):
        return burr._munp(n, c, 1.0)

    def _stats(self, c):
        return burr._stats(c, 1.0)

    def _entropy(self, c):
        return 2 - np.log(c)


fisk = fisk_gen(a=0.0, name='fisk')


class cauchy_gen(rv_continuous):
    r"""A Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `cauchy` is

    .. math::

        f(x) = \frac{1}{\pi (1 + x^2)}

    for a real number :math:`x`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # cauchy.pdf(x) = 1 / (pi * (1 + x**2))
        return 1.0/np.pi/(1.0+x*x)

    def _cdf(self, x):
        return 0.5 + 1.0/np.pi*np.arctan(x)

    def _ppf(self, q):
        return np.tan(np.pi*q-np.pi/2.0)

    def _sf(self, x):
        return 0.5 - 1.0/np.pi*np.arctan(x)

    def _isf(self, q):
        return np.tan(np.pi/2.0-np.pi*q)

    def _stats(self):
        return np.nan, np.nan, np.nan, np.nan

    def _entropy(self):
        return np.log(4*np.pi)

    def _fitstart(self, data, args=None):
        # Initialize ML guesses using quartiles instead of moments.
        p25, p50, p75 = np.percentile(data, [25, 50, 75])
        return p50, (p75 - p25)/2


cauchy = cauchy_gen(name='cauchy')


class chi_gen(rv_continuous):
    r"""A chi continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `chi` is:

    .. math::

        f(x, k) = \frac{1}{2^{k/2-1} \Gamma \left( k/2 \right)}
                   x^{k-1} \exp \left( -x^2/2 \right)

    for :math:`x >= 0` and :math:`k > 0` (degrees of freedom, denoted ``df``
    in the implementation). :math:`\Gamma` is the gamma function
    (`scipy.special.gamma`).

    Special cases of `chi` are:

        - ``chi(1, loc, scale)`` is equivalent to `halfnorm`
        - ``chi(2, 0, scale)`` is equivalent to `rayleigh`
        - ``chi(3, 0, scale)`` is equivalent to `maxwell`

    `chi` takes ``df`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("df", False, (0, np.inf), (False, False))]

    def _rvs(self, df, size=None, random_state=None):
        return np.sqrt(chi2.rvs(df, size=size, random_state=random_state))

    def _pdf(self, x, df):
        #                   x**(df-1) * exp(-x**2/2)
        # chi.pdf(x, df) =  -------------------------
        #                   2**(df/2-1) * gamma(df/2)
        return np.exp(self._logpdf(x, df))

    def _logpdf(self, x, df):
        l = np.log(2) - .5*np.log(2)*df - sc.gammaln(.5*df)
        return l + sc.xlogy(df - 1., x) - .5*x**2

    def _cdf(self, x, df):
        return sc.gammainc(.5*df, .5*x**2)

    def _sf(self, x, df):
        return sc.gammaincc(.5*df, .5*x**2)

    def _ppf(self, q, df):
        return np.sqrt(2*sc.gammaincinv(.5*df, q))

    def _isf(self, q, df):
        return np.sqrt(2*sc.gammainccinv(.5*df, q))

    def _stats(self, df):
        mu = np.sqrt(2)*np.exp(sc.gammaln(df/2.0+0.5)-sc.gammaln(df/2.0))
        mu2 = df - mu*mu
        g1 = (2*mu**3.0 + mu*(1-2*df))/np.asarray(np.power(mu2, 1.5))
        g2 = 2*df*(1.0-df)-6*mu**4 + 4*mu**2 * (2*df-1)
        g2 /= np.asarray(mu2**2.0)
        return mu, mu2, g1, g2


chi = chi_gen(a=0.0, name='chi')


class chi2_gen(rv_continuous):
    r"""A chi-squared continuous random variable.

    For the noncentral chi-square distribution, see `ncx2`.

    %(before_notes)s

    See Also
    --------
    ncx2

    Notes
    -----
    The probability density function for `chi2` is:

    .. math::

        f(x, k) = \frac{1}{2^{k/2} \Gamma \left( k/2 \right)}
                   x^{k/2-1} \exp \left( -x/2 \right)

    for :math:`x > 0`  and :math:`k > 0` (degrees of freedom, denoted ``df``
    in the implementation).

    `chi2` takes ``df`` as a shape parameter.

    The chi-squared distribution is a special case of the gamma
    distribution, with gamma parameters ``a = df/2``, ``loc = 0`` and
    ``scale = 2``.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("df", False, (0, np.inf), (False, False))]

    def _rvs(self, df, size=None, random_state=None):
        return random_state.chisquare(df, size)

    def _pdf(self, x, df):
        # chi2.pdf(x, df) = 1 / (2*gamma(df/2)) * (x/2)**(df/2-1) * exp(-x/2)
        return np.exp(self._logpdf(x, df))

    def _logpdf(self, x, df):
        return sc.xlogy(df/2.-1, x) - x/2. - sc.gammaln(df/2.) - (np.log(2)*df)/2.

    def _cdf(self, x, df):
        return sc.chdtr(df, x)

    def _sf(self, x, df):
        return sc.chdtrc(df, x)

    def _isf(self, p, df):
        return sc.chdtri(df, p)

    def _ppf(self, p, df):
        return 2*sc.gammaincinv(df/2, p)

    def _stats(self, df):
        mu = df
        mu2 = 2*df
        g1 = 2*np.sqrt(2.0/df)
        g2 = 12.0/df
        return mu, mu2, g1, g2


chi2 = chi2_gen(a=0.0, name='chi2')


class cosine_gen(rv_continuous):
    r"""A cosine continuous random variable.

    %(before_notes)s

    Notes
    -----
    The cosine distribution is an approximation to the normal distribution.
    The probability density function for `cosine` is:

    .. math::

        f(x) = \frac{1}{2\pi} (1+\cos(x))

    for :math:`-\pi \le x \le \pi`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # cosine.pdf(x) = 1/(2*pi) * (1+cos(x))
        return 1.0/2/np.pi*(1+np.cos(x))

    def _logpdf(self, x):
        c = np.cos(x)
        return _lazywhere(c != -1, (c,),
                          lambda c: np.log1p(c) - np.log(2*np.pi),
                          fillvalue=-np.inf)

    def _cdf(self, x):
        return scu._cosine_cdf(x)

    def _sf(self, x):
        return scu._cosine_cdf(-x)

    def _ppf(self, p):
        return scu._cosine_invcdf(p)

    def _isf(self, p):
        return -scu._cosine_invcdf(p)

    def _stats(self):
        return 0.0, np.pi*np.pi/3.0-2.0, 0.0, -6.0*(np.pi**4-90)/(5.0*(np.pi*np.pi-6)**2)

    def _entropy(self):
        return np.log(4*np.pi)-1.0


cosine = cosine_gen(a=-np.pi, b=np.pi, name='cosine')


class dgamma_gen(rv_continuous):
    r"""A double gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `dgamma` is:

    .. math::

        f(x, a) = \frac{1}{2\Gamma(a)} |x|^{a-1} \exp(-|x|)

    for a real number :math:`x` and :math:`a > 0`. :math:`\Gamma` is the
    gamma function (`scipy.special.gamma`).

    `dgamma` takes ``a`` as a shape parameter for :math:`a`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("a", False, (0, np.inf), (False, False))]

    def _rvs(self, a, size=None, random_state=None):
        u = random_state.uniform(size=size)
        gm = gamma.rvs(a, size=size, random_state=random_state)
        return gm * np.where(u >= 0.5, 1, -1)

    def _pdf(self, x, a):
        # dgamma.pdf(x, a) = 1 / (2*gamma(a)) * abs(x)**(a-1) * exp(-abs(x))
        ax = abs(x)
        return 1.0/(2*sc.gamma(a))*ax**(a-1.0) * np.exp(-ax)

    def _logpdf(self, x, a):
        ax = abs(x)
        return sc.xlogy(a - 1.0, ax) - ax - np.log(2) - sc.gammaln(a)

    def _cdf(self, x, a):
        fac = 0.5*sc.gammainc(a, abs(x))
        return np.where(x > 0, 0.5 + fac, 0.5 - fac)

    def _sf(self, x, a):
        fac = 0.5*sc.gammainc(a, abs(x))
        return np.where(x > 0, 0.5-fac, 0.5+fac)

    def _ppf(self, q, a):
        fac = sc.gammainccinv(a, 1-abs(2*q-1))
        return np.where(q > 0.5, fac, -fac)

    def _stats(self, a):
        mu2 = a*(a+1.0)
        return 0.0, mu2, 0.0, (a+2.0)*(a+3.0)/mu2-3.0


dgamma = dgamma_gen(name='dgamma')


class dweibull_gen(rv_continuous):
    r"""A double Weibull continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `dweibull` is given by

    .. math::

        f(x, c) = c / 2 |x|^{c-1} \exp(-|x|^c)

    for a real number :math:`x` and :math:`c > 0`.

    `dweibull` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _rvs(self, c, size=None, random_state=None):
        u = random_state.uniform(size=size)
        w = weibull_min.rvs(c, size=size, random_state=random_state)
        return w * (np.where(u >= 0.5, 1, -1))

    def _pdf(self, x, c):
        # dweibull.pdf(x, c) = c / 2 * abs(x)**(c-1) * exp(-abs(x)**c)
        ax = abs(x)
        Px = c / 2.0 * ax**(c-1.0) * np.exp(-ax**c)
        return Px

    def _logpdf(self, x, c):
        ax = abs(x)
        return np.log(c) - np.log(2.0) + sc.xlogy(c - 1.0, ax) - ax**c

    def _cdf(self, x, c):
        Cx1 = 0.5 * np.exp(-abs(x)**c)
        return np.where(x > 0, 1 - Cx1, Cx1)

    def _ppf(self, q, c):
        fac = 2. * np.where(q <= 0.5, q, 1. - q)
        fac = np.power(-np.log(fac), 1.0 / c)
        return np.where(q > 0.5, fac, -fac)

    def _munp(self, n, c):
        return (1 - (n % 2)) * sc.gamma(1.0 + 1.0 * n / c)

    # since we know that all odd moments are zeros, return them at once.
    # returning Nones from _stats makes the public stats call _munp
    # so overall we're saving one or two gamma function evaluations here.
    def _stats(self, c):
        return 0, None, 0, None


dweibull = dweibull_gen(name='dweibull')


class expon_gen(rv_continuous):
    r"""An exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `expon` is:

    .. math::

        f(x) = \exp(-x)

    for :math:`x \ge 0`.

    %(after_notes)s

    A common parameterization for `expon` is in terms of the rate parameter
    ``lambda``, such that ``pdf = lambda * exp(-lambda * x)``. This
    parameterization corresponds to using ``scale = 1 / lambda``.

    The exponential distribution is a special case of the gamma
    distributions, with gamma shape parameter ``a = 1``.

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return random_state.standard_exponential(size)

    def _pdf(self, x):
        # expon.pdf(x) = exp(-x)
        return np.exp(-x)

    def _logpdf(self, x):
        return -x

    def _cdf(self, x):
        return -sc.expm1(-x)

    def _ppf(self, q):
        return -sc.log1p(-q)

    def _sf(self, x):
        return np.exp(-x)

    def _logsf(self, x):
        return -x

    def _isf(self, q):
        return -np.log(q)

    def _stats(self):
        return 1.0, 1.0, 2.0, 6.0

    def _entropy(self):
        return 1.0

    @_call_super_mom
    @replace_notes_in_docstring(rv_continuous, notes="""\
        When `method='MLE'`,
        this function uses explicit formulas for the maximum likelihood
        estimation of the exponential distribution parameters, so the
        `optimizer`, `loc` and `scale` keyword arguments are
        ignored.\n\n""")
    def fit(self, data, *args, **kwds):
        if len(args) > 0:
            raise TypeError("Too many arguments.")

        floc = kwds.pop('floc', None)
        fscale = kwds.pop('fscale', None)

        _remove_optimizer_parameters(kwds)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        if not np.isfinite(data).all():
            raise ValueError("The data contains non-finite values.")

        data_min = data.min()

        if floc is None:
            # ML estimate of the location is the minimum of the data.
            loc = data_min
        else:
            loc = floc
            if data_min < loc:
                # There are values that are less than the specified loc.
                raise FitDataError("expon", lower=floc, upper=np.inf)

        if fscale is None:
            # ML estimate of the scale is the shifted mean.
            scale = data.mean() - loc
        else:
            scale = fscale

        # We expect the return values to be floating point, so ensure it
        # by explicitly converting to float.
        return float(loc), float(scale)


expon = expon_gen(a=0.0, name='expon')


class exponnorm_gen(rv_continuous):
    r"""An exponentially modified Normal continuous random variable.

    Also known as the exponentially modified Gaussian distribution [1]_.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponnorm` is:

    .. math::

        f(x, K) = \frac{1}{2K} \exp\left(\frac{1}{2 K^2} - x / K \right)
                  \text{erfc}\left(-\frac{x - 1/K}{\sqrt{2}}\right)

    where :math:`x` is a real number and :math:`K > 0`.

    It can be thought of as the sum of a standard normal random variable
    and an independent exponentially distributed random variable with rate
    ``1/K``.

    %(after_notes)s

    An alternative parameterization of this distribution (for example, in
    the Wikpedia article [1]_) involves three parameters, :math:`\mu`,
    :math:`\lambda` and :math:`\sigma`.

    In the present parameterization this corresponds to having ``loc`` and
    ``scale`` equal to :math:`\mu` and :math:`\sigma`, respectively, and
    shape parameter :math:`K = 1/(\sigma\lambda)`.

    .. versionadded:: 0.16.0

    References
    ----------
    .. [1] Exponentially modified Gaussian distribution, Wikipedia,
           https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("K", False, (0, np.inf), (False, False))]

    def _rvs(self, K, size=None, random_state=None):
        expval = random_state.standard_exponential(size) * K
        gval = random_state.standard_normal(size)
        return expval + gval

    def _pdf(self, x, K):
        return np.exp(self._logpdf(x, K))

    def _logpdf(self, x, K):
        invK = 1.0 / K
        exparg = invK * (0.5 * invK - x)
        return exparg + _norm_logcdf(x - invK) - np.log(K)

    def _cdf(self, x, K):
        invK = 1.0 / K
        expval = invK * (0.5 * invK - x)
        logprod = expval + _norm_logcdf(x - invK)
        return _norm_cdf(x) - np.exp(logprod)

    def _sf(self, x, K):
        invK = 1.0 / K
        expval = invK * (0.5 * invK - x)
        logprod = expval + _norm_logcdf(x - invK)
        return _norm_cdf(-x) + np.exp(logprod)

    def _stats(self, K):
        K2 = K * K
        opK2 = 1.0 + K2
        skw = 2 * K**3 * opK2**(-1.5)
        krt = 6.0 * K2 * K2 * opK2**(-2)
        return K, opK2, skw, krt


exponnorm = exponnorm_gen(name='exponnorm')


class exponweib_gen(rv_continuous):
    r"""An exponentiated Weibull continuous random variable.

    %(before_notes)s

    See Also
    --------
    weibull_min, numpy.random.Generator.weibull

    Notes
    -----
    The probability density function for `exponweib` is:

    .. math::

        f(x, a, c) = a c [1-\exp(-x^c)]^{a-1} \exp(-x^c) x^{c-1}

    and its cumulative distribution function is:

    .. math::

        F(x, a, c) = [1-\exp(-x^c)]^a

    for :math:`x > 0`, :math:`a > 0`, :math:`c > 0`.

    `exponweib` takes :math:`a` and :math:`c` as shape parameters:

    * :math:`a` is the exponentiation parameter,
      with the special case :math:`a=1` corresponding to the
      (non-exponentiated) Weibull distribution `weibull_min`.
    * :math:`c` is the shape parameter of the non-exponentiated Weibull law.

    %(after_notes)s

    References
    ----------
    https://en.wikipedia.org/wiki/Exponentiated_Weibull_distribution

    %(example)s

    """
    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ic = _ShapeInfo("c", False, (0, np.inf), (False, False))
        return [ia, ic]

    def _pdf(self, x, a, c):
        # exponweib.pdf(x, a, c) =
        #     a * c * (1-exp(-x**c))**(a-1) * exp(-x**c)*x**(c-1)
        return np.exp(self._logpdf(x, a, c))

    def _logpdf(self, x, a, c):
        negxc = -x**c
        exm1c = -sc.expm1(negxc)
        logp = (np.log(a) + np.log(c) + sc.xlogy(a - 1.0, exm1c) +
                negxc + sc.xlogy(c - 1.0, x))
        return logp

    def _cdf(self, x, a, c):
        exm1c = -sc.expm1(-x**c)
        return exm1c**a

    def _ppf(self, q, a, c):
        return (-sc.log1p(-q**(1.0/a)))**np.asarray(1.0/c)


exponweib = exponweib_gen(a=0.0, name='exponweib')


class exponpow_gen(rv_continuous):
    r"""An exponential power continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponpow` is:

    .. math::

        f(x, b) = b x^{b-1} \exp(1 + x^b - \exp(x^b))

    for :math:`x \ge 0`, :math:`b > 0`.  Note that this is a different
    distribution from the exponential power distribution that is also known
    under the names "generalized normal" or "generalized Gaussian".

    `exponpow` takes ``b`` as a shape parameter for :math:`b`.

    %(after_notes)s

    References
    ----------
    http://www.math.wm.edu/~leemis/chart/UDR/PDFs/Exponentialpower.pdf

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("b", False, (0, np.inf), (False, False))]

    def _pdf(self, x, b):
        # exponpow.pdf(x, b) = b * x**(b-1) * exp(1 + x**b - exp(x**b))
        return np.exp(self._logpdf(x, b))

    def _logpdf(self, x, b):
        xb = x**b
        f = 1 + np.log(b) + sc.xlogy(b - 1.0, x) + xb - np.exp(xb)
        return f

    def _cdf(self, x, b):
        return -sc.expm1(-sc.expm1(x**b))

    def _sf(self, x, b):
        return np.exp(-sc.expm1(x**b))

    def _isf(self, x, b):
        return (sc.log1p(-np.log(x)))**(1./b)

    def _ppf(self, q, b):
        return pow(sc.log1p(-sc.log1p(-q)), 1.0/b)


exponpow = exponpow_gen(a=0.0, name='exponpow')


class fatiguelife_gen(rv_continuous):
    r"""A fatigue-life (Birnbaum-Saunders) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `fatiguelife` is:

    .. math::

        f(x, c) = \frac{x+1}{2c\sqrt{2\pi x^3}} \exp(-\frac{(x-1)^2}{2x c^2})

    for :math:`x >= 0` and :math:`c > 0`.

    `fatiguelife` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    References
    ----------
    .. [1] "Birnbaum-Saunders distribution",
           https://en.wikipedia.org/wiki/Birnbaum-Saunders_distribution

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _rvs(self, c, size=None, random_state=None):
        z = random_state.standard_normal(size)
        x = 0.5*c*z
        x2 = x*x
        t = 1.0 + 2*x2 + 2*x*np.sqrt(1 + x2)
        return t

    def _pdf(self, x, c):
        # fatiguelife.pdf(x, c) =
        #     (x+1) / (2*c*sqrt(2*pi*x**3)) * exp(-(x-1)**2/(2*x*c**2))
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return (np.log(x+1) - (x-1)**2 / (2.0*x*c**2) - np.log(2*c) -
                0.5*(np.log(2*np.pi) + 3*np.log(x)))

    def _cdf(self, x, c):
        return _norm_cdf(1.0 / c * (np.sqrt(x) - 1.0/np.sqrt(x)))

    def _ppf(self, q, c):
        tmp = c*sc.ndtri(q)
        return 0.25 * (tmp + np.sqrt(tmp**2 + 4))**2

    def _sf(self, x, c):
        return _norm_sf(1.0 / c * (np.sqrt(x) - 1.0/np.sqrt(x)))

    def _isf(self, q, c):
        tmp = -c*sc.ndtri(q)
        return 0.25 * (tmp + np.sqrt(tmp**2 + 4))**2

    def _stats(self, c):
        # NB: the formula for kurtosis in wikipedia seems to have an error:
        # it's 40, not 41. At least it disagrees with the one from Wolfram
        # Alpha.  And the latter one, below, passes the tests, while the wiki
        # one doesn't So far I didn't have the guts to actually check the
        # coefficients from the expressions for the raw moments.
        c2 = c*c
        mu = c2 / 2.0 + 1.0
        den = 5.0 * c2 + 4.0
        mu2 = c2*den / 4.0
        g1 = 4 * c * (11*c2 + 6.0) / np.power(den, 1.5)
        g2 = 6 * c2 * (93*c2 + 40.0) / den**2.0
        return mu, mu2, g1, g2


fatiguelife = fatiguelife_gen(a=0.0, name='fatiguelife')


class foldcauchy_gen(rv_continuous):
    r"""A folded Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `foldcauchy` is:

    .. math::

        f(x, c) = \frac{1}{\pi (1+(x-c)^2)} + \frac{1}{\pi (1+(x+c)^2)}

    for :math:`x \ge 0` and :math:`c \ge 0`.

    `foldcauchy` takes ``c`` as a shape parameter for :math:`c`.

    %(example)s

    """
    def _argcheck(self, c):
        return c >= 0

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (True, False))]

    def _rvs(self, c, size=None, random_state=None):
        return abs(cauchy.rvs(loc=c, size=size,
                              random_state=random_state))

    def _pdf(self, x, c):
        # foldcauchy.pdf(x, c) = 1/(pi*(1+(x-c)**2)) + 1/(pi*(1+(x+c)**2))
        return 1.0/np.pi*(1.0/(1+(x-c)**2) + 1.0/(1+(x+c)**2))

    def _cdf(self, x, c):
        return 1.0/np.pi*(np.arctan(x-c) + np.arctan(x+c))

    def _stats(self, c):
        return np.inf, np.inf, np.nan, np.nan


foldcauchy = foldcauchy_gen(a=0.0, name='foldcauchy')


class f_gen(rv_continuous):
    r"""An F continuous random variable.

    For the noncentral F distribution, see `ncf`.

    %(before_notes)s

    See Also
    --------
    ncf

    Notes
    -----
    The probability density function for `f` is:

    .. math::

        f(x, df_1, df_2) = \frac{df_2^{df_2/2} df_1^{df_1/2} x^{df_1 / 2-1}}
                                {(df_2+df_1 x)^{(df_1+df_2)/2}
                                 B(df_1/2, df_2/2)}

    for :math:`x > 0` and parameters :math:`df_1, df_2 > 0` .

    `f` takes ``dfn`` and ``dfd`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        idfn = _ShapeInfo("dfn", False, (0, np.inf), (False, False))
        idfd = _ShapeInfo("dfd", False, (0, np.inf), (False, False))
        return [idfn, idfd]

    def _rvs(self, dfn, dfd, size=None, random_state=None):
        return random_state.f(dfn, dfd, size)

    def _pdf(self, x, dfn, dfd):
        #                      df2**(df2/2) * df1**(df1/2) * x**(df1/2-1)
        # F.pdf(x, df1, df2) = --------------------------------------------
        #                      (df2+df1*x)**((df1+df2)/2) * B(df1/2, df2/2)
        return np.exp(self._logpdf(x, dfn, dfd))

    def _logpdf(self, x, dfn, dfd):
        n = 1.0 * dfn
        m = 1.0 * dfd
        lPx = (m/2 * np.log(m) + n/2 * np.log(n) + sc.xlogy(n/2 - 1, x)
               - (((n+m)/2) * np.log(m + n*x) + sc.betaln(n/2, m/2)))
        return lPx

    def _cdf(self, x, dfn, dfd):
        return sc.fdtr(dfn, dfd, x)

    def _sf(self, x, dfn, dfd):
        return sc.fdtrc(dfn, dfd, x)

    def _ppf(self, q, dfn, dfd):
        return sc.fdtri(dfn, dfd, q)

    def _stats(self, dfn, dfd):
        v1, v2 = 1. * dfn, 1. * dfd
        v2_2, v2_4, v2_6, v2_8 = v2 - 2., v2 - 4., v2 - 6., v2 - 8.

        mu = _lazywhere(
            v2 > 2, (v2, v2_2),
            lambda v2, v2_2: v2 / v2_2,
            np.inf)

        mu2 = _lazywhere(
            v2 > 4, (v1, v2, v2_2, v2_4),
            lambda v1, v2, v2_2, v2_4:
            2 * v2 * v2 * (v1 + v2_2) / (v1 * v2_2**2 * v2_4),
            np.inf)

        g1 = _lazywhere(
            v2 > 6, (v1, v2_2, v2_4, v2_6),
            lambda v1, v2_2, v2_4, v2_6:
            (2 * v1 + v2_2) / v2_6 * np.sqrt(v2_4 / (v1 * (v1 + v2_2))),
            np.nan)
        g1 *= np.sqrt(8.)

        g2 = _lazywhere(
            v2 > 8, (g1, v2_6, v2_8),
            lambda g1, v2_6, v2_8: (8 + g1 * g1 * v2_6) / v2_8,
            np.nan)
        g2 *= 3. / 2.

        return mu, mu2, g1, g2


f = f_gen(a=0.0, name='f')


## Folded Normal
##   abs(Z) where (Z is normal with mu=L and std=S so that c=abs(L)/S)
##
##  note: regress docs have scale parameter correct, but first parameter
##    he gives is a shape parameter A = c * scale

##  Half-normal is folded normal with shape-parameter c=0.

class foldnorm_gen(rv_continuous):
    r"""A folded normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `foldnorm` is:

    .. math::

        f(x, c) = \sqrt{2/\pi} cosh(c x) \exp(-\frac{x^2+c^2}{2})

    for :math:`x \ge 0` and :math:`c \ge 0`.

    `foldnorm` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        return c >= 0

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (True, False))]

    def _rvs(self, c, size=None, random_state=None):
        return abs(random_state.standard_normal(size) + c)

    def _pdf(self, x, c):
        # foldnormal.pdf(x, c) = sqrt(2/pi) * cosh(c*x) * exp(-(x**2+c**2)/2)
        return _norm_pdf(x + c) + _norm_pdf(x-c)

    def _cdf(self, x, c):
        return _norm_cdf(x-c) + _norm_cdf(x+c) - 1.0

    def _stats(self, c):
        # Regina C. Elandt, Technometrics 3, 551 (1961)
        # https://www.jstor.org/stable/1266561
        #
        c2 = c*c
        expfac = np.exp(-0.5*c2) / np.sqrt(2.*np.pi)

        mu = 2.*expfac + c * sc.erf(c/np.sqrt(2))
        mu2 = c2 + 1 - mu*mu

        g1 = 2. * (mu*mu*mu - c2*mu - expfac)
        g1 /= np.power(mu2, 1.5)

        g2 = c2 * (c2 + 6.) + 3 + 8.*expfac*mu
        g2 += (2. * (c2 - 3.) - 3. * mu**2) * mu**2
        g2 = g2 / mu2**2.0 - 3.

        return mu, mu2, g1, g2


foldnorm = foldnorm_gen(a=0.0, name='foldnorm')


class weibull_min_gen(rv_continuous):
    r"""Weibull minimum continuous random variable.

    The Weibull Minimum Extreme Value distribution, from extreme value theory
    (Fisher-Gnedenko theorem), is also often simply called the Weibull
    distribution. It arises as the limiting distribution of the rescaled
    minimum of iid random variables.

    %(before_notes)s

    See Also
    --------
    weibull_max, numpy.random.Generator.weibull, exponweib

    Notes
    -----
    The probability density function for `weibull_min` is:

    .. math::

        f(x, c) = c x^{c-1} \exp(-x^c)

    for :math:`x > 0`, :math:`c > 0`.

    `weibull_min` takes ``c`` as a shape parameter for :math:`c`.
    (named :math:`k` in Wikipedia article and :math:`a` in
    ``numpy.random.weibull``).  Special shape values are :math:`c=1` and
    :math:`c=2` where Weibull distribution reduces to the `expon` and
    `rayleigh` distributions respectively.

    %(after_notes)s

    References
    ----------
    https://en.wikipedia.org/wiki/Weibull_distribution

    https://en.wikipedia.org/wiki/Fisher-Tippett-Gnedenko_theorem

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # weibull_min.pdf(x, c) = c * x**(c-1) * exp(-x**c)
        return c*pow(x, c-1)*np.exp(-pow(x, c))

    def _logpdf(self, x, c):
        return np.log(c) + sc.xlogy(c - 1, x) - pow(x, c)

    def _cdf(self, x, c):
        return -sc.expm1(-pow(x, c))

    def _sf(self, x, c):
        return np.exp(-pow(x, c))

    def _logsf(self, x, c):
        return -pow(x, c)

    def _ppf(self, q, c):
        return pow(-sc.log1p(-q), 1.0/c)

    def _munp(self, n, c):
        return sc.gamma(1.0+n*1.0/c)

    def _entropy(self, c):
        return -_EULER / c - np.log(c) + _EULER + 1

    @extend_notes_in_docstring(rv_continuous, notes="""\
        If ``method='mm'``, parameters fixed by the user are respected, and the
        remaining parameters are used to match distribution and sample moments
        where possible. For example, if the user fixes the location with
        ``floc``, the parameters will only match the distribution skewness and
        variance to the sample skewness and variance; no attempt will be made
        to match the means or minimize a norm of the errors.
        \n\n""")
    def fit(self, data, *args, **kwds):
        if kwds.pop('superfit', False):
            return super().fit(data, *args, **kwds)

        # this extracts fixed shape, location, and scale however they
        # are specified, and also leaves them in `kwds`
        data, fc, floc, fscale = _check_fit_input_parameters(self, data,
                                                             args, kwds)
        method = kwds.get("method", "mle").lower()

        # See https://en.wikipedia.org/wiki/Weibull_distribution#Moments for
        # moment formulas.
        def skew(c):
            gamma1 = sc.gamma(1+1/c)
            gamma2 = sc.gamma(1+2/c)
            gamma3 = sc.gamma(1+3/c)
            num = 2 * gamma1**3 - 3*gamma1*gamma2 + gamma3
            den = (gamma2 - gamma1**2)**(3/2)
            return num/den

        # For c in [1e2, 3e4], population skewness appears to approach
        # asymptote near -1.139, but past c > 3e4, skewness begins to vary
        # wildly, and MoM won't provide a good guess. Get out early.
        s = stats.skew(data)
        max_c = 1e4
        s_min = skew(max_c)
        if s < s_min and method == "mle" and fc is None and not args:
            return super().fit(data, *args, **kwds)

        # If method is method of moments, we don't need the user's guesses.
        # If method is "mle", extract the guesses from args and kwds.
        if method == "mm":
            c, loc, scale = None, None, None
        else:
            c = args[0] if len(args) else None
            loc = kwds.pop('loc', None)
            scale = kwds.pop('scale', None)

        if fc is None and c is None:  # not fixed and no guess: use MoM
            # Solve for c that matches sample distribution skewness to sample
            # skewness.
            # we start having numerical issues with `weibull_min` with
            # parameters outside this range - and not just in this method.
            # We could probably improve the situation by doing everything
            # in the log space, but that is for another time.
            c = root_scalar(lambda c: skew(c) - s, bracket=[0.02, max_c],
                            method='bisect').root
        elif fc is not None:  # fixed: use it
            c = fc

        if fscale is None and scale is None:
            v = np.var(data)
            scale = np.sqrt(v / (sc.gamma(1+2/c) - sc.gamma(1+1/c)**2))
        elif fscale is not None:
            scale = fscale

        if floc is None and loc is None:
            m = np.mean(data)
            loc = m - scale*sc.gamma(1 + 1/c)
        elif floc is not None:
            loc = floc

        if method == 'mm':
            return c, loc, scale
        else:
            # At this point, parameter "guesses" may equal the fixed parameters
            # in kwds. No harm in passing them as guesses, too.
            return super().fit(data, c, loc=loc, scale=scale, **kwds)


weibull_min = weibull_min_gen(a=0.0, name='weibull_min')


class truncweibull_min_gen(rv_continuous):
    r"""A doubly truncated Weibull minimum continuous random variable.

    %(before_notes)s

    See Also
    --------
    weibull_min, truncexpon

    Notes
    -----
    The probability density function for `truncweibull_min` is:

    .. math::

        f(x, a, b, c) = \frac{c x^{c-1} \exp(-x^c)}{\exp(-a^c) - \exp(-b^c)}

    for :math:`a < x <= b`, :math:`0 \le a < b` and :math:`c > 0`.

    `truncweibull_min` takes :math:`a`, :math:`b`, and :math:`c` as shape
    parameters.

    Notice that the truncation values, :math:`a` and :math:`b`, are defined in
    standardized form:

    .. math::

        a = (u_l - loc)/scale
        b = (u_r - loc)/scale

    where :math:`u_l` and :math:`u_r` are the specific left and right
    truncation values, respectively. In other words, the support of the
    distribution becomes :math:`(a*scale + loc) < x <= (b*scale + loc)` when
    :math:`loc` and/or :math:`scale` are provided.

    %(after_notes)s

    References
    ----------

    .. [1] Rinne, H. "The Weibull Distribution: A Handbook". CRC Press (2009).

    %(example)s

    """
    def _argcheck(self, c, a, b):
        return (a >= 0.) & (b > a) & (c > 0.)

    def _shape_info(self):
        ic = _ShapeInfo("c", False, (0, np.inf), (False, False))
        ia = _ShapeInfo("a", False, (0, np.inf), (True, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ic, ia, ib]

    def _fitstart(self, data):
        # Arbitrary, but default a=b=c=1 is not valid
        return super()._fitstart(data, args=(1, 0, 1))

    def _get_support(self, c, a, b):
        return a, b

    def _pdf(self, x, c, a, b):
        denum = (np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return (c * pow(x, c-1) * np.exp(-pow(x, c))) / denum

    def _logpdf(self, x, c, a, b):
        logdenum = np.log(np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return np.log(c) + sc.xlogy(c - 1, x) - pow(x, c) - logdenum

    def _cdf(self, x, c, a, b):
        num = (np.exp(-pow(a, c)) - np.exp(-pow(x, c)))
        denum = (np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return num / denum

    def _logcdf(self, x, c, a, b):
        lognum = np.log(np.exp(-pow(a, c)) - np.exp(-pow(x, c)))
        logdenum = np.log(np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return lognum - logdenum

    def _sf(self, x, c, a, b):
        num = (np.exp(-pow(x, c)) - np.exp(-pow(b, c)))
        denum = (np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return num / denum

    def _logsf(self, x, c, a, b):
        lognum = np.log(np.exp(-pow(x, c)) - np.exp(-pow(b, c)))
        logdenum = np.log(np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return lognum - logdenum

    def _isf(self, q, c, a, b):
        return pow(
            -np.log((1 - q) * np.exp(-pow(b, c)) + q * np.exp(-pow(a, c))), 1/c
            )

    def _ppf(self, q, c, a, b):
        return pow(
            -np.log((1 - q) * np.exp(-pow(a, c)) + q * np.exp(-pow(b, c))), 1/c
            )

    def _munp(self, n, c, a, b):
        gamma_fun = sc.gamma(n/c + 1.) * (
            sc.gammainc(n/c + 1., pow(b, c)) - sc.gammainc(n/c + 1., pow(a, c))
            )
        denum = (np.exp(-pow(a, c)) - np.exp(-pow(b, c)))
        return gamma_fun / denum


truncweibull_min = truncweibull_min_gen(name='truncweibull_min')


class weibull_max_gen(rv_continuous):
    r"""Weibull maximum continuous random variable.

    The Weibull Maximum Extreme Value distribution, from extreme value theory
    (Fisher-Gnedenko theorem), is the limiting distribution of rescaled
    maximum of iid random variables. This is the distribution of -X
    if X is from the `weibull_min` function.

    %(before_notes)s

    See Also
    --------
    weibull_min

    Notes
    -----
    The probability density function for `weibull_max` is:

    .. math::

        f(x, c) = c (-x)^{c-1} \exp(-(-x)^c)

    for :math:`x < 0`, :math:`c > 0`.

    `weibull_max` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    References
    ----------
    https://en.wikipedia.org/wiki/Weibull_distribution

    https://en.wikipedia.org/wiki/Fisher-Tippett-Gnedenko_theorem

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # weibull_max.pdf(x, c) = c * (-x)**(c-1) * exp(-(-x)**c)
        return c*pow(-x, c-1)*np.exp(-pow(-x, c))

    def _logpdf(self, x, c):
        return np.log(c) + sc.xlogy(c-1, -x) - pow(-x, c)

    def _cdf(self, x, c):
        return np.exp(-pow(-x, c))

    def _logcdf(self, x, c):
        return -pow(-x, c)

    def _sf(self, x, c):
        return -sc.expm1(-pow(-x, c))

    def _ppf(self, q, c):
        return -pow(-np.log(q), 1.0/c)

    def _munp(self, n, c):
        val = sc.gamma(1.0+n*1.0/c)
        if int(n) % 2:
            sgn = -1
        else:
            sgn = 1
        return sgn * val

    def _entropy(self, c):
        return -_EULER / c - np.log(c) + _EULER + 1


weibull_max = weibull_max_gen(b=0.0, name='weibull_max')


class genlogistic_gen(rv_continuous):
    r"""A generalized logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genlogistic` is:

    .. math::

        f(x, c) = c \frac{\exp(-x)}
                         {(1 + \exp(-x))^{c+1}}

    for :math:`x >= 0`, :math:`c > 0`.

    `genlogistic` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # genlogistic.pdf(x, c) = c * exp(-x) / (1 + exp(-x))**(c+1)
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        # Two mathematically equivalent expressions for log(pdf(x, c)):
        #     log(pdf(x, c)) = log(c) - x - (c + 1)*log(1 + exp(-x))
        #                    = log(c) + c*x - (c + 1)*log(1 + exp(x))
        mult = -(c - 1) * (x < 0) - 1
        absx = np.abs(x)
        return np.log(c) + mult*absx - (c+1) * sc.log1p(np.exp(-absx))

    def _cdf(self, x, c):
        Cx = (1+np.exp(-x))**(-c)
        return Cx

    def _ppf(self, q, c):
        vals = -np.log(pow(q, -1.0/c)-1)
        return vals

    def _stats(self, c):
        mu = _EULER + sc.psi(c)
        mu2 = np.pi*np.pi/6.0 + sc.zeta(2, c)
        g1 = -2*sc.zeta(3, c) + 2*_ZETA3
        g1 /= np.power(mu2, 1.5)
        g2 = np.pi**4/15.0 + 6*sc.zeta(4, c)
        g2 /= mu2**2.0
        return mu, mu2, g1, g2


genlogistic = genlogistic_gen(name='genlogistic')


class genpareto_gen(rv_continuous):
    r"""A generalized Pareto continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genpareto` is:

    .. math::

        f(x, c) = (1 + c x)^{-1 - 1/c}

    defined for :math:`x \ge 0` if :math:`c \ge 0`, and for
    :math:`0 \le x \le -1/c` if :math:`c < 0`.

    `genpareto` takes ``c`` as a shape parameter for :math:`c`.

    For :math:`c=0`, `genpareto` reduces to the exponential
    distribution, `expon`:

    .. math::

        f(x, 0) = \exp(-x)

    For :math:`c=-1`, `genpareto` is uniform on ``[0, 1]``:

    .. math::

        f(x, -1) = 1

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        return np.isfinite(c)

    def _shape_info(self):
        return [_ShapeInfo("c", False, (-np.inf, np.inf), (False, False))]

    def _get_support(self, c):
        c = np.asarray(c)
        b = _lazywhere(c < 0, (c,),
                       lambda c: -1. / c,
                       np.inf)
        a = np.where(c >= 0, self.a, self.a)
        return a, b

    def _pdf(self, x, c):
        # genpareto.pdf(x, c) = (1 + c * x)**(-1 - 1/c)
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.xlog1py(c + 1., c*x) / c,
                          -x)

    def _cdf(self, x, c):
        return -sc.inv_boxcox1p(-x, -c)

    def _sf(self, x, c):
        return sc.inv_boxcox(-x, -c)

    def _logsf(self, x, c):
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.log1p(c*x) / c,
                          -x)

    def _ppf(self, q, c):
        return -sc.boxcox1p(-q, -c)

    def _isf(self, q, c):
        return -sc.boxcox(q, -c)

    def _stats(self, c, moments='mv'):
        if 'm' not in moments:
            m = None
        else:
            m = _lazywhere(c < 1, (c,),
                           lambda xi: 1/(1 - xi),
                           np.inf)
        if 'v' not in moments:
            v = None
        else:
            v = _lazywhere(c < 1/2, (c,),
                           lambda xi: 1 / (1 - xi)**2 / (1 - 2*xi),
                           np.nan)
        if 's' not in moments:
            s = None
        else:
            s = _lazywhere(c < 1/3, (c,),
                           lambda xi: 2 * (1 + xi) * np.sqrt(1 - 2*xi) /
                                      (1 - 3*xi),
                           np.nan)
        if 'k' not in moments:
            k = None
        else:
            k = _lazywhere(c < 1/4, (c,),
                           lambda xi: 3 * (1 - 2*xi) * (2*xi**2 + xi + 3) /
                                      (1 - 3*xi) / (1 - 4*xi) - 3,
                           np.nan)
        return m, v, s, k

    def _munp(self, n, c):
        def __munp(n, c):
            val = 0.0
            k = np.arange(0, n + 1)
            for ki, cnk in zip(k, sc.comb(n, k)):
                val = val + cnk * (-1) ** ki / (1.0 - c * ki)
            return np.where(c * n < 1, val * (-1.0 / c) ** n, np.inf)
        return _lazywhere(c != 0, (c,),
                          lambda c: __munp(n, c),
                          sc.gamma(n + 1))

    def _entropy(self, c):
        return 1. + c


genpareto = genpareto_gen(a=0.0, name='genpareto')


class genexpon_gen(rv_continuous):
    r"""A generalized exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genexpon` is:

    .. math::

        f(x, a, b, c) = (a + b (1 - \exp(-c x)))
                        \exp(-a x - b x + \frac{b}{c}  (1-\exp(-c x)))

    for :math:`x \ge 0`, :math:`a, b, c > 0`.

    `genexpon` takes :math:`a`, :math:`b` and :math:`c` as shape parameters.

    %(after_notes)s

    References
    ----------
    H.K. Ryu, "An Extension of Marshall and Olkin's Bivariate Exponential
    Distribution", Journal of the American Statistical Association, 1993.

    N. Balakrishnan, "The Exponential Distribution: Theory, Methods and
    Applications", Asit P. Basu.

    %(example)s

    """
    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        ic = _ShapeInfo("c", False, (0, np.inf), (False, False))
        return [ia, ib, ic]

    def _pdf(self, x, a, b, c):
        # genexpon.pdf(x, a, b, c) = (a + b * (1 - exp(-c*x))) * \
        #                            exp(-a*x - b*x + b/c * (1-exp(-c*x)))
        return (a + b*(-sc.expm1(-c*x)))*np.exp((-a-b)*x +
                                                b*(-sc.expm1(-c*x))/c)

    def _logpdf(self, x, a, b, c):
        return np.log(a+b*(-sc.expm1(-c*x))) + (-a-b)*x+b*(-sc.expm1(-c*x))/c

    def _cdf(self, x, a, b, c):
        return -sc.expm1((-a-b)*x + b*(-sc.expm1(-c*x))/c)

    def _sf(self, x, a, b, c):
        return np.exp((-a-b)*x + b*(-sc.expm1(-c*x))/c)


genexpon = genexpon_gen(a=0.0, name='genexpon')


class genextreme_gen(rv_continuous):
    r"""A generalized extreme value continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_r

    Notes
    -----
    For :math:`c=0`, `genextreme` is equal to `gumbel_r` with
    probability density function

    .. math::

        f(x) = \exp(-\exp(-x)) \exp(-x),

    where :math:`-\infty < x < \infty`.

    For :math:`c \ne 0`, the probability density function for `genextreme` is:

    .. math::

        f(x, c) = \exp(-(1-c x)^{1/c}) (1-c x)^{1/c-1},

    where :math:`-\infty < x \le 1/c` if :math:`c > 0` and
    :math:`1/c \le x < \infty` if :math:`c < 0`.

    Note that several sources and software packages use the opposite
    convention for the sign of the shape parameter :math:`c`.

    `genextreme` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        return np.isfinite(c)

    def _shape_info(self):
        return [_ShapeInfo("c", False, (-np.inf, np.inf), (False, False))]

    def _get_support(self, c):
        _b = np.where(c > 0, 1.0 / np.maximum(c, _XMIN), np.inf)
        _a = np.where(c < 0, 1.0 / np.minimum(c, -_XMIN), -np.inf)
        return _a, _b

    def _loglogcdf(self, x, c):
        # Returns log(-log(cdf(x, c)))
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: sc.log1p(-c*x)/c, -x)

    def _pdf(self, x, c):
        # genextreme.pdf(x, c) =
        #     exp(-exp(-x))*exp(-x),                    for c==0
        #     exp(-(1-c*x)**(1/c))*(1-c*x)**(1/c-1),    for x \le 1/c, c > 0
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        cx = _lazywhere((x == x) & (c != 0), (x, c), lambda x, c: c*x, 0.0)
        logex2 = sc.log1p(-cx)
        logpex2 = self._loglogcdf(x, c)
        pex2 = np.exp(logpex2)
        # Handle special cases
        np.putmask(logpex2, (c == 0) & (x == -np.inf), 0.0)
        logpdf = _lazywhere(~((cx == 1) | (cx == -np.inf)),
                            (pex2, logpex2, logex2),
                            lambda pex2, lpex2, lex2: -pex2 + lpex2 - lex2,
                            fillvalue=-np.inf)
        np.putmask(logpdf, (c == 1) & (x == 1), 0.0)
        return logpdf

    def _logcdf(self, x, c):
        return -np.exp(self._loglogcdf(x, c))

    def _cdf(self, x, c):
        return np.exp(self._logcdf(x, c))

    def _sf(self, x, c):
        return -sc.expm1(self._logcdf(x, c))

    def _ppf(self, q, c):
        x = -np.log(-np.log(q))
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.expm1(-c * x) / c, x)

    def _isf(self, q, c):
        x = -np.log(-sc.log1p(-q))
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.expm1(-c * x) / c, x)

    def _stats(self, c):
        g = lambda n: sc.gamma(n*c + 1)
        g1 = g(1)
        g2 = g(2)
        g3 = g(3)
        g4 = g(4)
        g2mg12 = np.where(abs(c) < 1e-7, (c*np.pi)**2.0/6.0, g2-g1**2.0)
        gam2k = np.where(abs(c) < 1e-7, np.pi**2.0/6.0,
                         sc.expm1(sc.gammaln(2.0*c+1.0)-2*sc.gammaln(c + 1.0))/c**2.0)
        eps = 1e-14
        gamk = np.where(abs(c) < eps, -_EULER, sc.expm1(sc.gammaln(c + 1))/c)

        m = np.where(c < -1.0, np.nan, -gamk)
        v = np.where(c < -0.5, np.nan, g1**2.0*gam2k)

        # skewness
        sk1 = _lazywhere(c >= -1./3,
                         (c, g1, g2, g3, g2mg12),
                         lambda c, g1, g2, g3, g2gm12:
                             np.sign(c)*(-g3 + (g2 + 2*g2mg12)*g1)/g2mg12**1.5,
                         fillvalue=np.nan)
        sk = np.where(abs(c) <= eps**0.29, 12*np.sqrt(6)*_ZETA3/np.pi**3, sk1)

        # kurtosis
        ku1 = _lazywhere(c >= -1./4,
                         (g1, g2, g3, g4, g2mg12),
                         lambda g1, g2, g3, g4, g2mg12:
                             (g4 + (-4*g3 + 3*(g2 + g2mg12)*g1)*g1)/g2mg12**2,
                         fillvalue=np.nan)
        ku = np.where(abs(c) <= (eps)**0.23, 12.0/5.0, ku1-3.0)
        return m, v, sk, ku

    def _fitstart(self, data):
        # This is better than the default shape of (1,).
        g = _skew(data)
        if g < 0:
            a = 0.5
        else:
            a = -0.5
        return super()._fitstart(data, args=(a,))

    def _munp(self, n, c):
        k = np.arange(0, n+1)
        vals = 1.0/c**n * np.sum(
            sc.comb(n, k) * (-1)**k * sc.gamma(c*k + 1),
            axis=0)
        return np.where(c*n > -1, vals, np.inf)

    def _entropy(self, c):
        return _EULER*(1 - c) + 1


genextreme = genextreme_gen(name='genextreme')


def _digammainv(y):
    """Inverse of the digamma function (real positive arguments only).

    This function is used in the `fit` method of `gamma_gen`.
    The function uses either optimize.fsolve or optimize.newton
    to solve `sc.digamma(x) - y = 0`.  There is probably room for
    improvement, but currently it works over a wide range of y:

    >>> import numpy as np
    >>> rng = np.random.default_rng()
    >>> y = 64*rng.standard_normal(1000000)
    >>> y.min(), y.max()
    (-311.43592651416662, 351.77388222276869)
    >>> x = [_digammainv(t) for t in y]
    >>> np.abs(sc.digamma(x) - y).max()
    1.1368683772161603e-13

    """
    _em = 0.5772156649015328606065120
    func = lambda x: sc.digamma(x) - y
    if y > -0.125:
        x0 = np.exp(y) + 0.5
        if y < 10:
            # Some experimentation shows that newton reliably converges
            # must faster than fsolve in this y range.  For larger y,
            # newton sometimes fails to converge.
            value = optimize.newton(func, x0, tol=1e-10)
            return value
    elif y > -3:
        x0 = np.exp(y/2.332) + 0.08661
    else:
        x0 = 1.0 / (-y - _em)

    value, info, ier, mesg = optimize.fsolve(func, x0, xtol=1e-11,
                                             full_output=True)
    if ier != 1:
        raise RuntimeError("_digammainv: fsolve failed, y = %r" % y)

    return value[0]


## Gamma (Use MATLAB and MATHEMATICA (b=theta=scale, a=alpha=shape) definition)

## gamma(a, loc, scale)  with a an integer is the Erlang distribution
## gamma(1, loc, scale)  is the Exponential distribution
## gamma(df/2, 0, 2) is the chi2 distribution with df degrees of freedom.

class gamma_gen(rv_continuous):
    r"""A gamma continuous random variable.

    %(before_notes)s

    See Also
    --------
    erlang, expon

    Notes
    -----
    The probability density function for `gamma` is:

    .. math::

        f(x, a) = \frac{x^{a-1} e^{-x}}{\Gamma(a)}

    for :math:`x \ge 0`, :math:`a > 0`. Here :math:`\Gamma(a)` refers to the
    gamma function.

    `gamma` takes ``a`` as a shape parameter for :math:`a`.

    When :math:`a` is an integer, `gamma` reduces to the Erlang
    distribution, and when :math:`a=1` to the exponential distribution.

    Gamma distributions are sometimes parameterized with two variables,
    with a probability density function of:

    .. math::

        f(x, \alpha, \beta) = \frac{\beta^\alpha x^{\alpha - 1} e^{-\beta x }}{\Gamma(\alpha)}

    Note that this parameterization is equivalent to the above, with
    ``scale = 1 / beta``.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("a", False, (0, np.inf), (False, False))]

    def _rvs(self, a, size=None, random_state=None):
        return random_state.standard_gamma(a, size)

    def _pdf(self, x, a):
        # gamma.pdf(x, a) = x**(a-1) * exp(-x) / gamma(a)
        return np.exp(self._logpdf(x, a))

    def _logpdf(self, x, a):
        return sc.xlogy(a-1.0, x) - x - sc.gammaln(a)

    def _cdf(self, x, a):
        return sc.gammainc(a, x)

    def _sf(self, x, a):
        return sc.gammaincc(a, x)

    def _ppf(self, q, a):
        return sc.gammaincinv(a, q)

    def _isf(self, q, a):
        return sc.gammainccinv(a, q)

    def _stats(self, a):
        return a, a, 2.0/np.sqrt(a), 6.0/a

    def _entropy(self, a):
        return sc.psi(a)*(1-a) + a + sc.gammaln(a)

    def _fitstart(self, data):
        # The skewness of the gamma distribution is `2 / np.sqrt(a)`.
        # We invert that to estimate the shape `a` using the skewness
        # of the data.  The formula is regularized with 1e-8 in the
        # denominator to allow for degenerate data where the skewness
        # is close to 0.
        a = 4 / (1e-8 + _skew(data)**2)
        return super()._fitstart(data, args=(a,))

    @extend_notes_in_docstring(rv_continuous, notes="""\
        When the location is fixed by using the argument `floc`
        and `method='MLE'`, this
        function uses explicit formulas or solves a simpler numerical
        problem than the full ML optimization problem.  So in that case,
        the `optimizer`, `loc` and `scale` arguments are ignored.
        \n\n""")
    def fit(self, data, *args, **kwds):
        floc = kwds.get('floc', None)
        method = kwds.get('method', 'mle')

        if floc is None or method.lower() == 'mm':
            # loc is not fixed.  Use the default fit method.
            return super().fit(data, *args, **kwds)

        # We already have this value, so just pop it from kwds.
        kwds.pop('floc', None)

        f0 = _get_fixed_fit_value(kwds, ['f0', 'fa', 'fix_a'])
        fscale = kwds.pop('fscale', None)

        _remove_optimizer_parameters(kwds)

        # Special case: loc is fixed.

        if f0 is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            # Without this check, this function would just return the
            # parameters that were given.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        # Fixed location is handled by shifting the data.
        data = np.asarray(data)

        if not np.isfinite(data).all():
            raise ValueError("The data contains non-finite values.")

        if np.any(data <= floc):
            raise FitDataError("gamma", lower=floc, upper=np.inf)

        if floc != 0:
            # Don't do the subtraction in-place, because `data` might be a
            # view of the input array.
            data = data - floc
        xbar = data.mean()

        # Three cases to handle:
        # * shape and scale both free
        # * shape fixed, scale free
        # * shape free, scale fixed

        if fscale is None:
            # scale is free
            if f0 is not None:
                # shape is fixed
                a = f0
            else:
                # shape and scale are both free.
                # The MLE for the shape parameter `a` is the solution to:
                # np.log(a) - sc.digamma(a) - np.log(xbar) +
                #                             np.log(data).mean() = 0
                s = np.log(xbar) - np.log(data).mean()
                func = lambda a: np.log(a) - sc.digamma(a) - s
                aest = (3-s + np.sqrt((s-3)**2 + 24*s)) / (12*s)
                xa = aest*(1-0.4)
                xb = aest*(1+0.4)
                a = optimize.brentq(func, xa, xb, disp=0)

            # The MLE for the scale parameter is just the data mean
            # divided by the shape parameter.
            scale = xbar / a
        else:
            # scale is fixed, shape is free
            # The MLE for the shape parameter `a` is the solution to:
            # sc.digamma(a) - np.log(data).mean() + np.log(fscale) = 0
            c = np.log(data).mean() - np.log(fscale)
            a = _digammainv(c)
            scale = fscale

        return a, floc, scale


gamma = gamma_gen(a=0.0, name='gamma')


class erlang_gen(gamma_gen):
    """An Erlang continuous random variable.

    %(before_notes)s

    See Also
    --------
    gamma

    Notes
    -----
    The Erlang distribution is a special case of the Gamma distribution, with
    the shape parameter `a` an integer.  Note that this restriction is not
    enforced by `erlang`. It will, however, generate a warning the first time
    a non-integer value is used for the shape parameter.

    Refer to `gamma` for examples.

    """

    def _argcheck(self, a):
        allint = np.all(np.floor(a) == a)
        if not allint:
            # An Erlang distribution shouldn't really have a non-integer
            # shape parameter, so warn the user.
            warnings.warn(
                'The shape parameter of the erlang distribution '
                'has been given a non-integer value %r.' % (a,),
                RuntimeWarning)
        return a > 0

    def _shape_info(self):
        return [_ShapeInfo("a", True, (1, np.inf), (True, False))]

    def _fitstart(self, data):
        # Override gamma_gen_fitstart so that an integer initial value is
        # used.  (Also regularize the division, to avoid issues when
        # _skew(data) is 0 or close to 0.)
        a = int(4.0 / (1e-8 + _skew(data)**2))
        return super(gamma_gen, self)._fitstart(data, args=(a,))

    # Trivial override of the fit method, so we can monkey-patch its
    # docstring.
    @extend_notes_in_docstring(rv_continuous, notes="""\
        The Erlang distribution is generally defined to have integer values
        for the shape parameter.  This is not enforced by the `erlang` class.
        When fitting the distribution, it will generally return a non-integer
        value for the shape parameter.  By using the keyword argument
        `f0=<integer>`, the fit method can be constrained to fit the data to
        a specific integer shape parameter.""")
    def fit(self, data, *args, **kwds):
        return super().fit(data, *args, **kwds)


erlang = erlang_gen(a=0.0, name='erlang')


class gengamma_gen(rv_continuous):
    r"""A generalized gamma continuous random variable.

    %(before_notes)s

    See Also
    --------
    gamma, invgamma, weibull_min

    Notes
    -----
    The probability density function for `gengamma` is ([1]_):

    .. math::

        f(x, a, c) = \frac{|c| x^{c a-1} \exp(-x^c)}{\Gamma(a)}

    for :math:`x \ge 0`, :math:`a > 0`, and :math:`c \ne 0`.
    :math:`\Gamma` is the gamma function (`scipy.special.gamma`).

    `gengamma` takes :math:`a` and :math:`c` as shape parameters.

    %(after_notes)s

    References
    ----------
    .. [1] E.W. Stacy, "A Generalization of the Gamma Distribution",
       Annals of Mathematical Statistics, Vol 33(3), pp. 1187--1192.

    %(example)s

    """
    def _argcheck(self, a, c):
        return (a > 0) & (c != 0)

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ic = _ShapeInfo("c", False, (-np.inf, np.inf), (False, False))
        return [ia, ic]

    def _pdf(self, x, a, c):
        return np.exp(self._logpdf(x, a, c))

    def _logpdf(self, x, a, c):
        return _lazywhere((x != 0) | (c > 0), (x, c),
                          lambda x, c: (np.log(abs(c)) + sc.xlogy(c*a - 1, x)
                                        - x**c - sc.gammaln(a)),
                          fillvalue=-np.inf)

    def _cdf(self, x, a, c):
        xc = x**c
        val1 = sc.gammainc(a, xc)
        val2 = sc.gammaincc(a, xc)
        return np.where(c > 0, val1, val2)

    def _rvs(self, a, c, size=None, random_state=None):
        r = random_state.standard_gamma(a, size=size)
        return r**(1./c)

    def _sf(self, x, a, c):
        xc = x**c
        val1 = sc.gammainc(a, xc)
        val2 = sc.gammaincc(a, xc)
        return np.where(c > 0, val2, val1)

    def _ppf(self, q, a, c):
        val1 = sc.gammaincinv(a, q)
        val2 = sc.gammainccinv(a, q)
        return np.where(c > 0, val1, val2)**(1.0/c)

    def _isf(self, q, a, c):
        val1 = sc.gammaincinv(a, q)
        val2 = sc.gammainccinv(a, q)
        return np.where(c > 0, val2, val1)**(1.0/c)

    def _munp(self, n, a, c):
        # Pochhammer symbol: sc.pocha,n) = gamma(a+n)/gamma(a)
        return sc.poch(a, n*1.0/c)

    def _entropy(self, a, c):
        val = sc.psi(a)
        return a*(1-val) + 1.0/c*val + sc.gammaln(a) - np.log(abs(c))


gengamma = gengamma_gen(a=0.0, name='gengamma')


class genhalflogistic_gen(rv_continuous):
    r"""A generalized half-logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genhalflogistic` is:

    .. math::

        f(x, c) = \frac{2 (1 - c x)^{1/(c-1)}}{[1 + (1 - c x)^{1/c}]^2}

    for :math:`0 \le x \le 1/c`, and :math:`c > 0`.

    `genhalflogistic` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _get_support(self, c):
        return self.a, 1.0/c

    def _pdf(self, x, c):
        # genhalflogistic.pdf(x, c) =
        #    2 * (1-c*x)**(1/c-1) / (1+(1-c*x)**(1/c))**2
        limit = 1.0/c
        tmp = np.asarray(1-c*x)
        tmp0 = tmp**(limit-1)
        tmp2 = tmp0*tmp
        return 2*tmp0 / (1+tmp2)**2

    def _cdf(self, x, c):
        limit = 1.0/c
        tmp = np.asarray(1-c*x)
        tmp2 = tmp**(limit)
        return (1.0-tmp2) / (1+tmp2)

    def _ppf(self, q, c):
        return 1.0/c*(1-((1.0-q)/(1.0+q))**c)

    def _entropy(self, c):
        return 2 - (2*c+1)*np.log(2)


genhalflogistic = genhalflogistic_gen(a=0.0, name='genhalflogistic')


class genhyperbolic_gen(rv_continuous):
    r"""A generalized hyperbolic continuous random variable.

    %(before_notes)s

    See Also
    --------
    t, norminvgauss, geninvgauss, laplace, cauchy

    Notes
    -----
    The probability density function for `genhyperbolic` is:

    .. math::

        f(x, p, a, b) =
            \frac{(a^2 - b^2)^{p/2}}
            {\sqrt{2\pi}a^{p-0.5}
            K_p\Big(\sqrt{a^2 - b^2}\Big)}
            e^{bx} \times \frac{K_{p - 1/2}
            (a \sqrt{1 + x^2})}
            {(\sqrt{1 + x^2})^{1/2 - p}}

    for :math:`x, p \in ( - \infty; \infty)`,
    :math:`|b| < a` if :math:`p \ge 0`,
    :math:`|b| \le a` if :math:`p < 0`.
    :math:`K_{p}(.)` denotes the modified Bessel function of the second
    kind and order :math:`p` (`scipy.special.kn`)

    `genhyperbolic` takes ``p`` as a tail parameter,
    ``a`` as a shape parameter,
    ``b`` as a skewness parameter.

    %(after_notes)s

    The original parameterization of the Generalized Hyperbolic Distribution
    is found in [1]_ as follows

    .. math::

        f(x, \lambda, \alpha, \beta, \delta, \mu) =
           \frac{(\gamma/\delta)^\lambda}{\sqrt{2\pi}K_\lambda(\delta \gamma)}
           e^{\beta (x - \mu)} \times \frac{K_{\lambda - 1/2}
           (\alpha \sqrt{\delta^2 + (x - \mu)^2})}
           {(\sqrt{\delta^2 + (x - \mu)^2} / \alpha)^{1/2 - \lambda}}

    for :math:`x \in ( - \infty; \infty)`,
    :math:`\gamma := \sqrt{\alpha^2 - \beta^2}`,
    :math:`\lambda, \mu \in ( - \infty; \infty)`,
    :math:`\delta \ge 0, |\beta| < \alpha` if :math:`\lambda \ge 0`,
    :math:`\delta > 0, |\beta| \le \alpha` if :math:`\lambda < 0`.

    The location-scale-based parameterization implemented in
    SciPy is based on [2]_, where :math:`a = \alpha\delta`,
    :math:`b = \beta\delta`, :math:`p = \lambda`,
    :math:`scale=\delta` and :math:`loc=\mu`

    Moments are implemented based on [3]_ and [4]_.

    For the distributions that are a special case such as Student's t,
    it is not recommended to rely on the implementation of genhyperbolic.
    To avoid potential numerical problems and for performance reasons,
    the methods of the specific distributions should be used.

    References
    ----------
    .. [1] O. Barndorff-Nielsen, "Hyperbolic Distributions and Distributions
       on Hyperbolae", Scandinavian Journal of Statistics, Vol. 5(3),
       pp. 151-157, 1978. https://www.jstor.org/stable/4615705

    .. [2] Eberlein E., Prause K. (2002) The Generalized Hyperbolic Model:
        Financial Derivatives and Risk Measures. In: Geman H., Madan D.,
        Pliska S.R., Vorst T. (eds) Mathematical Finance - Bachelier
        Congress 2000. Springer Finance. Springer, Berlin, Heidelberg.
        :doi:`10.1007/978-3-662-12429-1_12`

    .. [3] Scott, David J, Würtz, Diethelm, Dong, Christine and Tran,
       Thanh Tam, (2009), Moments of the generalized hyperbolic
       distribution, MPRA Paper, University Library of Munich, Germany,
       https://EconPapers.repec.org/RePEc:pra:mprapa:19081.

    .. [4] E. Eberlein and E. A. von Hammerstein. Generalized hyperbolic
       and inverse Gaussian distributions: Limiting cases and approximation
       of processes. FDM Preprint 80, April 2003. University of Freiburg.
       https://freidok.uni-freiburg.de/fedora/objects/freidok:7974/datastreams/FILE1/content

    %(example)s

    """

    def _argcheck(self, p, a, b):
        return (np.logical_and(np.abs(b) < a, p >= 0)
                | np.logical_and(np.abs(b) <= a, p < 0))

    def _shape_info(self):
        ip = _ShapeInfo("p", False, (-np.inf, np.inf), (False, False))
        ia = _ShapeInfo("a", False, (0, np.inf), (True, False))
        ib = _ShapeInfo("b", False, (-np.inf, np.inf), (False, False))
        return [ip, ia, ib]

    def _fitstart(self, data):
        # Arbitrary, but the default a=b=1 is not valid
        return super()._fitstart(data, args=(1, 1, 0.5))

    def _logpdf(self, x, p, a, b):
        # kve instead of kv works better for large values of p
        # and smaller values of sqrt(a^2  - b^2)
        @np.vectorize
        def _logpdf_single(x, p, a, b):
            return _stats.genhyperbolic_logpdf(x, p, a, b)

        return _logpdf_single(x, p, a, b)

    def _pdf(self, x, p, a, b):
        # kve instead of kv works better for large values of p
        # and smaller values of sqrt(a^2  - b^2)
        @np.vectorize
        def _pdf_single(x, p, a, b):
            return _stats.genhyperbolic_pdf(x, p, a, b)

        return _pdf_single(x, p, a, b)

    def _cdf(self, x, p, a, b):

        @np.vectorize
        def _cdf_single(x, p, a, b):
            user_data = np.array(
                [p, a, b], float
                ).ctypes.data_as(ctypes.c_void_p)
            llc = LowLevelCallable.from_cython(
                _stats, '_genhyperbolic_pdf', user_data
                )

            t1 = integrate.quad(llc, -np.inf, x)[0]

            if np.isnan(t1):
                msg = ("Infinite values encountered in scipy.special.kve. "
                       "Values replaced by NaN to avoid incorrect results.")
                warnings.warn(msg, RuntimeWarning)

            return t1

        return _cdf_single(x, p, a, b)

    def _rvs(self, p, a, b, size=None, random_state=None):
        # note: X = b * V + sqrt(V) * X  has a
        # generalized hyperbolic distribution
        # if X is standard normal and V is
        # geninvgauss(p = p, b = t2, loc = loc, scale = t3)
        t1 = np.float_power(a, 2) - np.float_power(b, 2)
        # b in the GIG
        t2 = np.float_power(t1, 0.5)
        # scale in the GIG
        t3 = np.float_power(t1, - 0.5)
        gig = geninvgauss.rvs(
            p=p,
            b=t2,
            scale=t3,
            size=size,
            random_state=random_state
            )
        normst = norm.rvs(size=size, random_state=random_state)

        return b * gig + np.sqrt(gig) * normst

    def _stats(self, p, a, b):
        # https://mpra.ub.uni-muenchen.de/19081/1/MPRA_paper_19081.pdf
        # https://freidok.uni-freiburg.de/fedora/objects/freidok:7974/datastreams/FILE1/content
        # standardized moments
        p, a, b = np.broadcast_arrays(p, a, b)
        t1 = np.float_power(a, 2) - np.float_power(b, 2)
        t1 = np.float_power(t1, 0.5)
        t2 = np.float_power(1, 2) * np.float_power(t1, - 1)
        integers = np.linspace(0, 4, 5)
        # make integers perpendicular to existing dimensions
        integers = integers.reshape(integers.shape + (1,) * p.ndim)
        b0, b1, b2, b3, b4 = sc.kv(p + integers, t1)
        r1, r2, r3, r4 = [b / b0 for b in (b1, b2, b3, b4)]

        m = b * t2 * r1
        v = (
            t2 * r1 + np.float_power(b, 2) * np.float_power(t2, 2) *
            (r2 - np.float_power(r1, 2))
        )
        m3e = (
            np.float_power(b, 3) * np.float_power(t2, 3) *
            (r3 - 3 * b2 * b1 * np.float_power(b0, -2) +
             2 * np.float_power(r1, 3)) +
            3 * b * np.float_power(t2, 2) *
            (r2 - np.float_power(r1, 2))
        )
        s = m3e * np.float_power(v, - 3 / 2)
        m4e = (
            np.float_power(b, 4) * np.float_power(t2, 4) *
            (r4 - 4 * b3 * b1 * np.float_power(b0, - 2) +
             6 * b2 * np.float_power(b1, 2) * np.float_power(b0, - 3) -
             3 * np.float_power(r1, 4)) +
            np.float_power(b, 2) * np.float_power(t2, 3) *
            (6 * r3 - 12 * b2 * b1 * np.float_power(b0, - 2) +
             6 * np.float_power(r1, 3)) +
            3 * np.float_power(t2, 2) * r2
        )
        k = m4e * np.float_power(v, -2) - 3

        return m, v, s, k


genhyperbolic = genhyperbolic_gen(name='genhyperbolic')


class gompertz_gen(rv_continuous):
    r"""A Gompertz (or truncated Gumbel) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gompertz` is:

    .. math::

        f(x, c) = c \exp(x) \exp(-c (e^x-1))

    for :math:`x \ge 0`, :math:`c > 0`.

    `gompertz` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # gompertz.pdf(x, c) = c * exp(x) * exp(-c*(exp(x)-1))
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return np.log(c) + x - c * sc.expm1(x)

    def _cdf(self, x, c):
        return -sc.expm1(-c * sc.expm1(x))

    def _ppf(self, q, c):
        return sc.log1p(-1.0 / c * sc.log1p(-q))

    def _entropy(self, c):
        return 1.0 - np.log(c) - np.exp(c)*sc.expn(1, c)


gompertz = gompertz_gen(a=0.0, name='gompertz')


def _average_with_log_weights(x, logweights):
    x = np.asarray(x)
    logweights = np.asarray(logweights)
    maxlogw = logweights.max()
    weights = np.exp(logweights - maxlogw)
    return np.average(x, weights=weights)


class gumbel_r_gen(rv_continuous):
    r"""A right-skewed Gumbel continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_l, gompertz, genextreme

    Notes
    -----
    The probability density function for `gumbel_r` is:

    .. math::

        f(x) = \exp(-(x + e^{-x}))

    The Gumbel distribution is sometimes referred to as a type I Fisher-Tippett
    distribution.  It is also related to the extreme value distribution,
    log-Weibull and Gompertz distributions.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # gumbel_r.pdf(x) = exp(-(x + exp(-x)))
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return -x - np.exp(-x)

    def _cdf(self, x):
        return np.exp(-np.exp(-x))

    def _logcdf(self, x):
        return -np.exp(-x)

    def _ppf(self, q):
        return -np.log(-np.log(q))

    def _sf(self, x):
        return -sc.expm1(-np.exp(-x))

    def _isf(self, p):
        return -np.log(-np.log1p(-p))

    def _stats(self):
        return _EULER, np.pi*np.pi/6.0, 12*np.sqrt(6)/np.pi**3 * _ZETA3, 12.0/5

    def _entropy(self):
        # https://en.wikipedia.org/wiki/Gumbel_distribution
        return _EULER + 1.

    @_call_super_mom
    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        data, floc, fscale = _check_fit_input_parameters(self, data,
                                                         args, kwds)

        # By the method of maximum likelihood, the estimators of the
        # location and scale are the roots of the equations defined in
        # `func` and the value of the expression for `loc` that follows.
        # The first `func` is a first order derivative of the log-likelihood
        # equation and the second is from Source: Statistical Distributions,
        # 3rd Edition. Evans, Hastings, and Peacock (2000), Page 101.

        def get_loc_from_scale(scale):
            return -scale * (sc.logsumexp(-data / scale) - np.log(len(data)))

        if fscale is not None:
            # if the scale is fixed, the location can be analytically
            # determined.
            scale = fscale
            loc = get_loc_from_scale(scale)
        else:
            # A different function is solved depending on whether the location
            # is fixed.
            if floc is not None:
                loc = floc

                # equation to use if the location is fixed.
                # note that one cannot use the equation in Evans, Hastings,
                # and Peacock (2000) (since it assumes that the derivative
                # w.r.t. the log-likelihood is zero). however, it is easy to
                # derive the MLE condition directly if loc is fixed
                def func(scale):
                    term1 = (loc - data) * np.exp((loc - data) / scale) + data
                    term2 = len(data) * (loc + scale)
                    return term1.sum() - term2
            else:

                # equation to use if both location and scale are free
                def func(scale):
                    sdata = -data / scale
                    wavg = _average_with_log_weights(data, logweights=sdata)
                    return data.mean() - wavg - scale

            # set brackets for `root_scalar` to use when optimizing over the
            # scale such that a root is likely between them. Use user supplied
            # guess or default 1.
            brack_start = kwds.get('scale', 1)
            lbrack, rbrack = brack_start / 2, brack_start * 2

            # if a root is not between the brackets, iteratively expand them
            # until they include a sign change, checking after each bracket is
            # modified.
            def interval_contains_root(lbrack, rbrack):
                # return true if the signs disagree.
                return (np.sign(func(lbrack)) !=
                        np.sign(func(rbrack)))
            while (not interval_contains_root(lbrack, rbrack)
                   and (lbrack > 0 or rbrack < np.inf)):
                lbrack /= 2
                rbrack *= 2

            res = optimize.root_scalar(func, bracket=(lbrack, rbrack),
                                       rtol=1e-14, xtol=1e-14)
            scale = res.root
            loc = floc if floc is not None else get_loc_from_scale(scale)
        return loc, scale


gumbel_r = gumbel_r_gen(name='gumbel_r')


class gumbel_l_gen(rv_continuous):
    r"""A left-skewed Gumbel continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_r, gompertz, genextreme

    Notes
    -----
    The probability density function for `gumbel_l` is:

    .. math::

        f(x) = \exp(x - e^x)

    The Gumbel distribution is sometimes referred to as a type I Fisher-Tippett
    distribution.  It is also related to the extreme value distribution,
    log-Weibull and Gompertz distributions.

    %(after_notes)s

    %(example)s

    """

    def _shape_info(self):
        return []

    def _pdf(self, x):
        # gumbel_l.pdf(x) = exp(x - exp(x))
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return x - np.exp(x)

    def _cdf(self, x):
        return -sc.expm1(-np.exp(x))

    def _ppf(self, q):
        return np.log(-sc.log1p(-q))

    def _logsf(self, x):
        return -np.exp(x)

    def _sf(self, x):
        return np.exp(-np.exp(x))

    def _isf(self, x):
        return np.log(-np.log(x))

    def _stats(self):
        return -_EULER, np.pi*np.pi/6.0, \
               -12*np.sqrt(6)/np.pi**3 * _ZETA3, 12.0/5

    def _entropy(self):
        return _EULER + 1.

    @_call_super_mom
    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        # The fit method of `gumbel_r` can be used for this distribution with
        # small modifications. The process to do this is
        # 1. pass the sign negated data into `gumbel_r.fit`
        #    - if the location is fixed, it should also be negated.
        # 2. negate the sign of the resulting location, leaving the scale
        #    unmodified.
        # `gumbel_r.fit` holds necessary input checks.

        if kwds.get('floc') is not None:
            kwds['floc'] = -kwds['floc']
        loc_r, scale_r, = gumbel_r.fit(-np.asarray(data), *args, **kwds)
        return -loc_r, scale_r


gumbel_l = gumbel_l_gen(name='gumbel_l')


class halfcauchy_gen(rv_continuous):
    r"""A Half-Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfcauchy` is:

    .. math::

        f(x) = \frac{2}{\pi (1 + x^2)}

    for :math:`x \ge 0`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # halfcauchy.pdf(x) = 2 / (pi * (1 + x**2))
        return 2.0/np.pi/(1.0+x*x)

    def _logpdf(self, x):
        return np.log(2.0/np.pi) - sc.log1p(x*x)

    def _cdf(self, x):
        return 2.0/np.pi*np.arctan(x)

    def _ppf(self, q):
        return np.tan(np.pi/2*q)

    def _stats(self):
        return np.inf, np.inf, np.nan, np.nan

    def _entropy(self):
        return np.log(2*np.pi)


halfcauchy = halfcauchy_gen(a=0.0, name='halfcauchy')


class halflogistic_gen(rv_continuous):
    r"""A half-logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halflogistic` is:

    .. math::

        f(x) = \frac{ 2 e^{-x} }{ (1+e^{-x})^2 }
             = \frac{1}{2} \text{sech}(x/2)^2

    for :math:`x \ge 0`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # halflogistic.pdf(x) = 2 * exp(-x) / (1+exp(-x))**2
        #                     = 1/2 * sech(x/2)**2
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return np.log(2) - x - 2. * sc.log1p(np.exp(-x))

    def _cdf(self, x):
        return np.tanh(x/2.0)

    def _ppf(self, q):
        return 2*np.arctanh(q)

    def _munp(self, n):
        if n == 1:
            return 2*np.log(2)
        if n == 2:
            return np.pi*np.pi/3.0
        if n == 3:
            return 9*_ZETA3
        if n == 4:
            return 7*np.pi**4 / 15.0
        return 2*(1-pow(2.0, 1-n))*sc.gamma(n+1)*sc.zeta(n, 1)

    def _entropy(self):
        return 2-np.log(2)


halflogistic = halflogistic_gen(a=0.0, name='halflogistic')


class halfnorm_gen(rv_continuous):
    r"""A half-normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfnorm` is:

    .. math::

        f(x) = \sqrt{2/\pi} \exp(-x^2 / 2)

    for :math:`x >= 0`.

    `halfnorm` is a special case of `chi` with ``df=1``.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return abs(random_state.standard_normal(size=size))

    def _pdf(self, x):
        # halfnorm.pdf(x) = sqrt(2/pi) * exp(-x**2/2)
        return np.sqrt(2.0/np.pi)*np.exp(-x*x/2.0)

    def _logpdf(self, x):
        return 0.5 * np.log(2.0/np.pi) - x*x/2.0

    def _cdf(self, x):
        return _norm_cdf(x)*2-1.0

    def _ppf(self, q):
        return sc.ndtri((1+q)/2.0)

    def _stats(self):
        return (np.sqrt(2.0/np.pi),
                1-2.0/np.pi,
                np.sqrt(2)*(4-np.pi)/(np.pi-2)**1.5,
                8*(np.pi-3)/(np.pi-2)**2)

    def _entropy(self):
        return 0.5*np.log(np.pi/2.0)+0.5


halfnorm = halfnorm_gen(a=0.0, name='halfnorm')


class hypsecant_gen(rv_continuous):
    r"""A hyperbolic secant continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `hypsecant` is:

    .. math::

        f(x) = \frac{1}{\pi} \text{sech}(x)

    for a real number :math:`x`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        # hypsecant.pdf(x) = 1/pi * sech(x)
        return 1.0/(np.pi*np.cosh(x))

    def _cdf(self, x):
        return 2.0/np.pi*np.arctan(np.exp(x))

    def _ppf(self, q):
        return np.log(np.tan(np.pi*q/2.0))

    def _stats(self):
        return 0, np.pi*np.pi/4, 0, 2

    def _entropy(self):
        return np.log(2*np.pi)


hypsecant = hypsecant_gen(name='hypsecant')


class gausshyper_gen(rv_continuous):
    r"""A Gauss hypergeometric continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gausshyper` is:

    .. math::

        f(x, a, b, c, z) = C x^{a-1} (1-x)^{b-1} (1+zx)^{-c}

    for :math:`0 \le x \le 1`, :math:`a,b > 0`, :math:`c` a real number,
    :math:`z > -1`, and :math:`C = \frac{1}{B(a, b) F[2, 1](c, a; a+b; -z)}`.
    :math:`F[2, 1]` is the Gauss hypergeometric function
    `scipy.special.hyp2f1`.

    `gausshyper` takes :math:`a`, :math:`b`, :math:`c` and :math:`z` as shape
    parameters.

    %(after_notes)s

    References
    ----------
    .. [1] Armero, C., and M. J. Bayarri. "Prior Assessments for Prediction in
           Queues." *Journal of the Royal Statistical Society*. Series D (The
           Statistician) 43, no. 1 (1994): 139-53. doi:10.2307/2348939

    %(example)s

    """

    def _argcheck(self, a, b, c, z):
        # z > -1 per gh-10134
        return (a > 0) & (b > 0) & (c == c) & (z > -1)

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        ic = _ShapeInfo("c", False, (-np.inf, np.inf), (False, False))
        iz = _ShapeInfo("z", False, (-1, np.inf), (False, False))
        return [ia, ib, ic, iz]

    def _pdf(self, x, a, b, c, z):
        # gausshyper.pdf(x, a, b, c, z) =
        #   C * x**(a-1) * (1-x)**(b-1) * (1+z*x)**(-c)
        Cinv = sc.gamma(a)*sc.gamma(b)/sc.gamma(a+b)*sc.hyp2f1(c, a, a+b, -z)
        return 1.0/Cinv * x**(a-1.0) * (1.0-x)**(b-1.0) / (1.0+z*x)**c

    def _munp(self, n, a, b, c, z):
        fac = sc.beta(n+a, b) / sc.beta(a, b)
        num = sc.hyp2f1(c, a+n, a+b+n, -z)
        den = sc.hyp2f1(c, a, a+b, -z)
        return fac*num / den


gausshyper = gausshyper_gen(a=0.0, b=1.0, name='gausshyper')


class invgamma_gen(rv_continuous):
    r"""An inverted gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invgamma` is:

    .. math::

        f(x, a) = \frac{x^{-a-1}}{\Gamma(a)} \exp(-\frac{1}{x})

    for :math:`x >= 0`, :math:`a > 0`. :math:`\Gamma` is the gamma function
    (`scipy.special.gamma`).

    `invgamma` takes ``a`` as a shape parameter for :math:`a`.

    `invgamma` is a special case of `gengamma` with ``c=-1``, and it is a
    different parameterization of the scaled inverse chi-squared distribution.
    Specifically, if the scaled inverse chi-squared distribution is
    parameterized with degrees of freedom :math:`\nu` and scaling parameter
    :math:`\tau^2`, then it can be modeled using `invgamma` with
    ``a=`` :math:`\nu/2` and ``scale=`` :math:`\nu \tau^2/2`.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, a):
        # invgamma.pdf(x, a) = x**(-a-1) / gamma(a) * exp(-1/x)
        return np.exp(self._logpdf(x, a))

    def _logpdf(self, x, a):
        return -(a+1) * np.log(x) - sc.gammaln(a) - 1.0/x

    def _cdf(self, x, a):
        return sc.gammaincc(a, 1.0 / x)

    def _ppf(self, q, a):
        return 1.0 / sc.gammainccinv(a, q)

    def _sf(self, x, a):
        return sc.gammainc(a, 1.0 / x)

    def _isf(self, q, a):
        return 1.0 / sc.gammaincinv(a, q)

    def _stats(self, a, moments='mvsk'):
        m1 = _lazywhere(a > 1, (a,), lambda x: 1. / (x - 1.), np.inf)
        m2 = _lazywhere(a > 2, (a,), lambda x: 1. / (x - 1.)**2 / (x - 2.),
                        np.inf)

        g1, g2 = None, None
        if 's' in moments:
            g1 = _lazywhere(
                a > 3, (a,),
                lambda x: 4. * np.sqrt(x - 2.) / (x - 3.), np.nan)
        if 'k' in moments:
            g2 = _lazywhere(
                a > 4, (a,),
                lambda x: 6. * (5. * x - 11.) / (x - 3.) / (x - 4.), np.nan)
        return m1, m2, g1, g2

    def _entropy(self, a):
        return a - (a+1.0) * sc.psi(a) + sc.gammaln(a)


invgamma = invgamma_gen(a=0.0, name='invgamma')


class invgauss_gen(rv_continuous):
    r"""An inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invgauss` is:

    .. math::

        f(x, \mu) = \frac{1}{\sqrt{2 \pi x^3}}
                    \exp(-\frac{(x-\mu)^2}{2 x \mu^2})

    for :math:`x >= 0` and :math:`\mu > 0`.

    `invgauss` takes ``mu`` as a shape parameter for :math:`\mu`.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return [_ShapeInfo("mu", False, (0, np.inf), (False, False))]

    def _rvs(self, mu, size=None, random_state=None):
        return random_state.wald(mu, 1.0, size=size)

    def _pdf(self, x, mu):
        # invgauss.pdf(x, mu) =
        #                  1 / sqrt(2*pi*x**3) * exp(-(x-mu)**2/(2*x*mu**2))
        return 1.0/np.sqrt(2*np.pi*x**3.0)*np.exp(-1.0/(2*x)*((x-mu)/mu)**2)

    def _logpdf(self, x, mu):
        return -0.5*np.log(2*np.pi) - 1.5*np.log(x) - ((x-mu)/mu)**2/(2*x)

    # approach adapted from equations in
    # https://journal.r-project.org/archive/2016-1/giner-smyth.pdf,
    # not R code. see gh-13616

    def _logcdf(self, x, mu):
        fac = 1 / np.sqrt(x)
        a = _norm_logcdf(fac * ((x / mu) - 1))
        b = 2 / mu + _norm_logcdf(-fac * ((x / mu) + 1))
        return a + np.log1p(np.exp(b - a))

    def _logsf(self, x, mu):
        fac = 1 / np.sqrt(x)
        a = _norm_logsf(fac * ((x / mu) - 1))
        b = 2 / mu + _norm_logcdf(-fac * (x + mu) / mu)
        return a + np.log1p(-np.exp(b - a))

    def _sf(self, x, mu):
        return np.exp(self._logsf(x, mu))

    def _cdf(self, x, mu):
        return np.exp(self._logcdf(x, mu))

    def _ppf(self, x, mu):
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            x, mu = np.broadcast_arrays(x, mu)
            ppf = _boost._invgauss_ppf(x, mu, 1)
            i_wt = x > 0.5  # "wrong tail" - sometimes too inaccurate
            ppf[i_wt] = _boost._invgauss_isf(1-x[i_wt], mu[i_wt], 1)
            i_nan = np.isnan(ppf)
            ppf[i_nan] = super()._ppf(x[i_nan], mu[i_nan])
        return ppf

    def _isf(self, x, mu):
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            x, mu = np.broadcast_arrays(x, mu)
            isf = _boost._invgauss_isf(x, mu, 1)
            i_wt = x > 0.5  # "wrong tail" - sometimes too inaccurate
            isf[i_wt] = _boost._invgauss_ppf(1-x[i_wt], mu[i_wt], 1)
            i_nan = np.isnan(isf)
            isf[i_nan] = super()._isf(x[i_nan], mu[i_nan])
        return isf

    def _stats(self, mu):
        return mu, mu**3.0, 3*np.sqrt(mu), 15*mu

    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        method = kwds.get('method', 'mle')

        if type(self) == wald_gen or method.lower() == 'mm':
            return super().fit(data, *args, **kwds)

        data, fshape_s, floc, fscale = _check_fit_input_parameters(self, data,
                                                                   args, kwds)
        '''
        Source: Statistical Distributions, 3rd Edition. Evans, Hastings,
        and Peacock (2000), Page 121. Their shape parameter is equivilent to
        SciPy's with the conversion `fshape_s = fshape / scale`.

        MLE formulas are not used in 3 condtions:
        - `loc` is not fixed
        - `mu` is fixed
        These cases fall back on the superclass fit method.
        - `loc` is fixed but translation results in negative data raises
          a `FitDataError`.
        '''
        if floc is None or fshape_s is not None:
            return super().fit(data, *args, **kwds)
        elif np.any(data - floc < 0):
            raise FitDataError("invgauss", lower=0, upper=np.inf)
        else:
            data = data - floc
            fshape_n = np.mean(data)
            if fscale is None:
                fscale = len(data) / (np.sum(data ** -1 - fshape_n ** -1))
            fshape_s = fshape_n / fscale
        return fshape_s, floc, fscale


invgauss = invgauss_gen(a=0.0, name='invgauss')


class geninvgauss_gen(rv_continuous):
    r"""A Generalized Inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `geninvgauss` is:

    .. math::

        f(x, p, b) = x^{p-1} \exp(-b (x + 1/x) / 2) / (2 K_p(b))

    where `x > 0`, `p` is a real number and `b > 0`([1]_).
    :math:`K_p` is the modified Bessel function of second kind of order `p`
    (`scipy.special.kv`).

    %(after_notes)s

    The inverse Gaussian distribution `stats.invgauss(mu)` is a special case of
    `geninvgauss` with `p = -1/2`, `b = 1 / mu` and `scale = mu`.

    Generating random variates is challenging for this distribution. The
    implementation is based on [2]_.

    References
    ----------
    .. [1] O. Barndorff-Nielsen, P. Blaesild, C. Halgreen, "First hitting time
       models for the generalized inverse gaussian distribution",
       Stochastic Processes and their Applications 7, pp. 49--54, 1978.

    .. [2] W. Hoermann and J. Leydold, "Generating generalized inverse Gaussian
       random variates", Statistics and Computing, 24(4), p. 547--557, 2014.

    %(example)s

    """
    def _argcheck(self, p, b):
        return (p == p) & (b > 0)

    def _shape_info(self):
        ip = _ShapeInfo("p", False, (-np.inf, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ip, ib]

    def _logpdf(self, x, p, b):
        # kve instead of kv works better for large values of b
        # warn if kve produces infinite values and replace by nan
        # otherwise c = -inf and the results are often incorrect
        @np.vectorize
        def logpdf_single(x, p, b):
            return _stats.geninvgauss_logpdf(x, p, b)

        z = logpdf_single(x, p, b)
        if np.isnan(z).any():
            msg = ("Infinite values encountered in scipy.special.kve(p, b). "
                   "Values replaced by NaN to avoid incorrect results.")
            warnings.warn(msg, RuntimeWarning)
        return z

    def _pdf(self, x, p, b):
        # relying on logpdf avoids overflow of x**(p-1) for large x and p
        return np.exp(self._logpdf(x, p, b))

    def _cdf(self, x, *args):
        _a, _b = self._get_support(*args)

        @np.vectorize
        def _cdf_single(x, *args):
            p, b = args
            user_data = np.array([p, b], float).ctypes.data_as(ctypes.c_void_p)
            llc = LowLevelCallable.from_cython(_stats, '_geninvgauss_pdf',
                                               user_data)

            return integrate.quad(llc, _a, x)[0]

        return _cdf_single(x, *args)

    def _logquasipdf(self, x, p, b):
        # log of the quasi-density (w/o normalizing constant) used in _rvs
        return _lazywhere(x > 0, (x, p, b),
                          lambda x, p, b: (p - 1)*np.log(x) - b*(x + 1/x)/2,
                          -np.inf)

    def _rvs(self, p, b, size=None, random_state=None):
        # if p and b are scalar, use _rvs_scalar, otherwise need to create
        # output by iterating over parameters
        if np.isscalar(p) and np.isscalar(b):
            out = self._rvs_scalar(p, b, size, random_state)
        elif p.size == 1 and b.size == 1:
            out = self._rvs_scalar(p.item(), b.item(), size, random_state)
        else:
            # When this method is called, size will be a (possibly empty)
            # tuple of integers.  It will not be None; if `size=None` is passed
            # to `rvs()`, size will be the empty tuple ().

            p, b = np.broadcast_arrays(p, b)
            # p and b now have the same shape.

            # `shp` is the shape of the blocks of random variates that are
            # generated for each combination of parameters associated with
            # broadcasting p and b.
            # bc is a tuple the same lenth as size.  The values
            # in bc are bools.  If bc[j] is True, it means that
            # entire axis is filled in for a given combination of the
            # broadcast arguments.
            shp, bc = _check_shape(p.shape, size)

            # `numsamples` is the total number of variates to be generated
            # for each combination of the input arguments.
            numsamples = int(np.prod(shp))

            # `out` is the array to be returned.  It is filled in in the
            # loop below.
            out = np.empty(size)

            it = np.nditer([p, b],
                           flags=['multi_index'],
                           op_flags=[['readonly'], ['readonly']])
            while not it.finished:
                # Convert the iterator's multi_index into an index into the
                # `out` array where the call to _rvs_scalar() will be stored.
                # Where bc is True, we use a full slice; otherwise we use the
                # index value from it.multi_index.  len(it.multi_index) might
                # be less than len(bc), and in that case we want to align these
                # two sequences to the right, so the loop variable j runs from
                # -len(size) to 0.  This doesn't cause an IndexError, as
                # bc[j] will be True in those cases where it.multi_index[j]
                # would cause an IndexError.
                idx = tuple((it.multi_index[j] if not bc[j] else slice(None))
                            for j in range(-len(size), 0))
                out[idx] = self._rvs_scalar(it[0], it[1], numsamples,
                                            random_state).reshape(shp)
                it.iternext()

        if size == ():
            out = out.item()
        return out

    def _rvs_scalar(self, p, b, numsamples, random_state):
        # following [2], the quasi-pdf is used instead of the pdf for the
        # generation of rvs
        invert_res = False
        if not(numsamples):
            numsamples = 1
        if p < 0:
            # note: if X is geninvgauss(p, b), then 1/X is geninvgauss(-p, b)
            p = -p
            invert_res = True
        m = self._mode(p, b)

        # determine method to be used following [2]
        ratio_unif = True
        if p >= 1 or b > 1:
            # ratio of uniforms with mode shift below
            mode_shift = True
        elif b >= min(0.5, 2 * np.sqrt(1 - p) / 3):
            # ratio of uniforms without mode shift below
            mode_shift = False
        else:
            # new algorithm in [2]
            ratio_unif = False

        # prepare sampling of rvs
        size1d = tuple(np.atleast_1d(numsamples))
        N = np.prod(size1d)  # number of rvs needed, reshape upon return
        x = np.zeros(N)
        simulated = 0

        if ratio_unif:
            # use ratio of uniforms method
            if mode_shift:
                a2 = -2 * (p + 1) / b - m
                a1 = 2 * m * (p - 1) / b - 1
                # find roots of x**3 + a2*x**2 + a1*x + m (Cardano's formula)
                p1 = a1 - a2**2 / 3
                q1 = 2 * a2**3 / 27 - a2 * a1 / 3 + m
                phi = np.arccos(-q1 * np.sqrt(-27 / p1**3) / 2)
                s1 = -np.sqrt(-4 * p1 / 3)
                root1 = s1 * np.cos(phi / 3 + np.pi / 3) - a2 / 3
                root2 = -s1 * np.cos(phi / 3) - a2 / 3
                # root3 = s1 * np.cos(phi / 3 - np.pi / 3) - a2 / 3

                # if g is the quasipdf, rescale: g(x) / g(m) which we can write
                # as exp(log(g(x)) - log(g(m))). This is important
                # since for large values of p and b, g cannot be evaluated.
                # denote the rescaled quasipdf by h
                lm = self._logquasipdf(m, p, b)
                d1 = self._logquasipdf(root1, p, b) - lm
                d2 = self._logquasipdf(root2, p, b) - lm
                # compute the bounding rectangle w.r.t. h. Note that
                # np.exp(0.5*d1) = np.sqrt(g(root1)/g(m)) = np.sqrt(h(root1))
                vmin = (root1 - m) * np.exp(0.5 * d1)
                vmax = (root2 - m) * np.exp(0.5 * d2)
                umax = 1  # umax = sqrt(h(m)) = 1

                logqpdf = lambda x: self._logquasipdf(x, p, b) - lm
                c = m
            else:
                # ratio of uniforms without mode shift
                # compute np.sqrt(quasipdf(m))
                umax = np.exp(0.5*self._logquasipdf(m, p, b))
                xplus = ((1 + p) + np.sqrt((1 + p)**2 + b**2))/b
                vmin = 0
                # compute xplus * np.sqrt(quasipdf(xplus))
                vmax = xplus * np.exp(0.5 * self._logquasipdf(xplus, p, b))
                c = 0
                logqpdf = lambda x: self._logquasipdf(x, p, b)

            if vmin >= vmax:
                raise ValueError("vmin must be smaller than vmax.")
            if umax <= 0:
                raise ValueError("umax must be positive.")

            i = 1
            while simulated < N:
                k = N - simulated
                # simulate uniform rvs on [0, umax] and [vmin, vmax]
                u = umax * random_state.uniform(size=k)
                v = random_state.uniform(size=k)
                v = vmin + (vmax - vmin) * v
                rvs = v / u + c
                # rewrite acceptance condition u**2 <= pdf(rvs) by taking logs
                accept = (2*np.log(u) <= logqpdf(rvs))
                num_accept = np.sum(accept)
                if num_accept > 0:
                    x[simulated:(simulated + num_accept)] = rvs[accept]
                    simulated += num_accept

                if (simulated == 0) and (i*N >= 50000):
                    msg = ("Not a single random variate could be generated "
                           "in {} attempts. Sampling does not appear to "
                           "work for the provided parameters.".format(i*N))
                    raise RuntimeError(msg)
                i += 1
        else:
            # use new algorithm in [2]
            x0 = b / (1 - p)
            xs = np.max((x0, 2 / b))
            k1 = np.exp(self._logquasipdf(m, p, b))
            A1 = k1 * x0
            if x0 < 2 / b:
                k2 = np.exp(-b)
                if p > 0:
                    A2 = k2 * ((2 / b)**p - x0**p) / p
                else:
                    A2 = k2 * np.log(2 / b**2)
            else:
                k2, A2 = 0, 0
            k3 = xs**(p - 1)
            A3 = 2 * k3 * np.exp(-xs * b / 2) / b
            A = A1 + A2 + A3

            # [2]: rejection constant is < 2.73; so expected runtime is finite
            while simulated < N:
                k = N - simulated
                h, rvs = np.zeros(k), np.zeros(k)
                # simulate uniform rvs on [x1, x2] and [0, y2]
                u = random_state.uniform(size=k)
                v = A * random_state.uniform(size=k)
                cond1 = v <= A1
                cond2 = np.logical_not(cond1) & (v <= A1 + A2)
                cond3 = np.logical_not(cond1 | cond2)
                # subdomain (0, x0)
                rvs[cond1] = x0 * v[cond1] / A1
                h[cond1] = k1
                # subdomain (x0, 2 / b)
                if p > 0:
                    rvs[cond2] = (x0**p + (v[cond2] - A1) * p / k2)**(1 / p)
                else:
                    rvs[cond2] = b * np.exp((v[cond2] - A1) * np.exp(b))
                h[cond2] = k2 * rvs[cond2]**(p - 1)
                # subdomain (xs, infinity)
                z = np.exp(-xs * b / 2) - b * (v[cond3] - A1 - A2) / (2 * k3)
                rvs[cond3] = -2 / b * np.log(z)
                h[cond3] = k3 * np.exp(-rvs[cond3] * b / 2)
                # apply rejection method
                accept = (np.log(u * h) <= self._logquasipdf(rvs, p, b))
                num_accept = sum(accept)
                if num_accept > 0:
                    x[simulated:(simulated + num_accept)] = rvs[accept]
                    simulated += num_accept

        rvs = np.reshape(x, size1d)
        if invert_res:
            rvs = 1 / rvs
        return rvs

    def _mode(self, p, b):
        # distinguish cases to avoid catastrophic cancellation (see [2])
        if p < 1:
            return b / (np.sqrt((p - 1)**2 + b**2) + 1 - p)
        else:
            return (np.sqrt((1 - p)**2 + b**2) - (1 - p)) / b

    def _munp(self, n, p, b):
        num = sc.kve(p + n, b)
        denom = sc.kve(p, b)
        inf_vals = np.isinf(num) | np.isinf(denom)
        if inf_vals.any():
            msg = ("Infinite values encountered in the moment calculation "
                   "involving scipy.special.kve. Values replaced by NaN to "
                   "avoid incorrect results.")
            warnings.warn(msg, RuntimeWarning)
            m = np.full_like(num, np.nan, dtype=np.double)
            m[~inf_vals] = num[~inf_vals] / denom[~inf_vals]
        else:
            m = num / denom
        return m


geninvgauss = geninvgauss_gen(a=0.0, name="geninvgauss")


class norminvgauss_gen(rv_continuous):
    r"""A Normal Inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `norminvgauss` is:

    .. math::

        f(x, a, b) = \frac{a \, K_1(a \sqrt{1 + x^2})}{\pi \sqrt{1 + x^2}} \,
                     \exp(\sqrt{a^2 - b^2} + b x)

    where :math:`x` is a real number, the parameter :math:`a` is the tail
    heaviness and :math:`b` is the asymmetry parameter satisfying
    :math:`a > 0` and :math:`|b| <= a`.
    :math:`K_1` is the modified Bessel function of second kind
    (`scipy.special.k1`).

    %(after_notes)s

    A normal inverse Gaussian random variable `Y` with parameters `a` and `b`
    can be expressed as a normal mean-variance mixture:
    `Y = b * V + sqrt(V) * X` where `X` is `norm(0,1)` and `V` is
    `invgauss(mu=1/sqrt(a**2 - b**2))`. This representation is used
    to generate random variates.

    Another common parametrization of the distribution (see Equation 2.1 in
    [2]_) is given by the following expression of the pdf:

    .. math::

        g(x, \alpha, \beta, \delta, \mu) =
        \frac{\alpha\delta K_1\left(\alpha\sqrt{\delta^2 + (x - \mu)^2}\right)}
        {\pi \sqrt{\delta^2 + (x - \mu)^2}} \,
        e^{\delta \sqrt{\alpha^2 - \beta^2} + \beta (x - \mu)}

    In SciPy, this corresponds to
    `a = alpha * delta, b = beta * delta, loc = mu, scale=delta`.

    References
    ----------
    .. [1] O. Barndorff-Nielsen, "Hyperbolic Distributions and Distributions on
           Hyperbolae", Scandinavian Journal of Statistics, Vol. 5(3),
           pp. 151-157, 1978.

    .. [2] O. Barndorff-Nielsen, "Normal Inverse Gaussian Distributions and
           Stochastic Volatility Modelling", Scandinavian Journal of
           Statistics, Vol. 24, pp. 1-13, 1997.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _argcheck(self, a, b):
        return (a > 0) & (np.absolute(b) < a)

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (-np.inf, np.inf), (False, False))
        return [ia, ib]

    def _fitstart(self, data):
        # Arbitrary, but the default a=b=1 is not valid
        return super()._fitstart(data, args=(1, 0.5))

    def _pdf(self, x, a, b):
        gamma = np.sqrt(a**2 - b**2)
        fac1 = a / np.pi * np.exp(gamma)
        sq = np.hypot(1, x)  # reduce overflows
        return fac1 * sc.k1e(a * sq) * np.exp(b*x - a*sq) / sq

    def _sf(self, x, a, b):
        if np.isscalar(x):
            # If x is a scalar, then so are a and b.
            return integrate.quad(self._pdf, x, np.inf, args=(a, b))[0]
        else:
            result = []
            for (x0, a0, b0) in zip(x, a, b):
                result.append(integrate.quad(self._pdf, x0, np.inf,
                                             args=(a0, b0))[0])
            return np.array(result)

    def _isf(self, q, a, b):
        def _isf_scalar(q, a, b):

            def eq(x, a, b, q):
                # Solve eq(x, a, b, q) = 0 to obtain isf(x, a, b) = q.
                return self._sf(x, a, b) - q

            # Find a bracketing interval for the root.
            # Start at the mean, and grow the length of the interval
            # by 2 each iteration until there is a sign change in eq.
            xm = self.mean(a, b)
            em = eq(xm, a, b, q)
            if em == 0:
                # Unlikely, but might as well check.
                return xm
            if em > 0:
                delta = 1
                left = xm
                right = xm + delta
                while eq(right, a, b, q) > 0:
                    delta = 2*delta
                    right = xm + delta
            else:
                # em < 0
                delta = 1
                right = xm
                left = xm - delta
                while eq(left, a, b, q) < 0:
                    delta = 2*delta
                    left = xm - delta
            result = optimize.brentq(eq, left, right, args=(a, b, q),
                                     xtol=self.xtol)
            return result

        if np.isscalar(q):
            return _isf_scalar(q, a, b)
        else:
            result = []
            for (q0, a0, b0) in zip(q, a, b):
                result.append(_isf_scalar(q0, a0, b0))
            return np.array(result)

    def _rvs(self, a, b, size=None, random_state=None):
        # note: X = b * V + sqrt(V) * X is norminvgaus(a,b) if X is standard
        # normal and V is invgauss(mu=1/sqrt(a**2 - b**2))
        gamma = np.sqrt(a**2 - b**2)
        ig = invgauss.rvs(mu=1/gamma, size=size, random_state=random_state)
        return b * ig + np.sqrt(ig) * norm.rvs(size=size,
                                               random_state=random_state)

    def _stats(self, a, b):
        gamma = np.sqrt(a**2 - b**2)
        mean = b / gamma
        variance = a**2 / gamma**3
        skewness = 3.0 * b / (a * np.sqrt(gamma))
        kurtosis = 3.0 * (1 + 4 * b**2 / a**2) / gamma
        return mean, variance, skewness, kurtosis


norminvgauss = norminvgauss_gen(name="norminvgauss")


class invweibull_gen(rv_continuous):
    u"""An inverted Weibull continuous random variable.

    This distribution is also known as the Fréchet distribution or the
    type II extreme value distribution.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invweibull` is:

    .. math::

        f(x, c) = c x^{-c-1} \\exp(-x^{-c})

    for :math:`x > 0`, :math:`c > 0`.

    `invweibull` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    References
    ----------
    F.R.S. de Gusmao, E.M.M Ortega and G.M. Cordeiro, "The generalized inverse
    Weibull distribution", Stat. Papers, vol. 52, pp. 591-619, 2011.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # invweibull.pdf(x, c) = c * x**(-c-1) * exp(-x**(-c))
        xc1 = np.power(x, -c - 1.0)
        xc2 = np.power(x, -c)
        xc2 = np.exp(-xc2)
        return c * xc1 * xc2

    def _cdf(self, x, c):
        xc1 = np.power(x, -c)
        return np.exp(-xc1)

    def _ppf(self, q, c):
        return np.power(-np.log(q), -1.0/c)

    def _munp(self, n, c):
        return sc.gamma(1 - n / c)

    def _entropy(self, c):
        return 1+_EULER + _EULER / c - np.log(c)

    def _fitstart(self, data, args=None):
        # invweibull requires c > 1 for the first moment to exist, so use 2.0
        args = (2.0,) if args is None else args
        return super(invweibull_gen, self)._fitstart(data, args=args)


invweibull = invweibull_gen(a=0, name='invweibull')


class johnsonsb_gen(rv_continuous):
    r"""A Johnson SB continuous random variable.

    %(before_notes)s

    See Also
    --------
    johnsonsu

    Notes
    -----
    The probability density function for `johnsonsb` is:

    .. math::

        f(x, a, b) = \frac{b}{x(1-x)}  \phi(a + b \log \frac{x}{1-x} )

    where :math:`x`, :math:`a`, and :math:`b` are real scalars; :math:`b > 0`
    and :math:`x \in [0,1]`.  :math:`\phi` is the pdf of the normal
    distribution.

    `johnsonsb` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _argcheck(self, a, b):
        return (b > 0) & (a == a)

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (-np.inf, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ia, ib]

    def _pdf(self, x, a, b):
        # johnsonsb.pdf(x, a, b) = b / (x*(1-x)) * phi(a + b * log(x/(1-x)))
        trm = _norm_pdf(a + b*np.log(x/(1.0-x)))
        return b*1.0/(x*(1-x))*trm

    def _cdf(self, x, a, b):
        return _norm_cdf(a + b*np.log(x/(1.0-x)))

    def _ppf(self, q, a, b):
        return 1.0 / (1 + np.exp(-1.0 / b * (_norm_ppf(q) - a)))


johnsonsb = johnsonsb_gen(a=0.0, b=1.0, name='johnsonsb')


class johnsonsu_gen(rv_continuous):
    r"""A Johnson SU continuous random variable.

    %(before_notes)s

    See Also
    --------
    johnsonsb

    Notes
    -----
    The probability density function for `johnsonsu` is:

    .. math::

        f(x, a, b) = \frac{b}{\sqrt{x^2 + 1}}
                     \phi(a + b \log(x + \sqrt{x^2 + 1}))

    where :math:`x`, :math:`a`, and :math:`b` are real scalars; :math:`b > 0`.
    :math:`\phi` is the pdf of the normal distribution.

    `johnsonsu` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a, b):
        return (b > 0) & (a == a)

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (-np.inf, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ia, ib]

    def _pdf(self, x, a, b):
        # johnsonsu.pdf(x, a, b) = b / sqrt(x**2 + 1) *
        #                          phi(a + b * log(x + sqrt(x**2 + 1)))
        x2 = x*x
        trm = _norm_pdf(a + b * np.log(x + np.sqrt(x2+1)))
        return b*1.0/np.sqrt(x2+1.0)*trm

    def _cdf(self, x, a, b):
        return _norm_cdf(a + b * np.log(x + np.sqrt(x*x + 1)))

    def _ppf(self, q, a, b):
        return np.sinh((_norm_ppf(q) - a) / b)


johnsonsu = johnsonsu_gen(name='johnsonsu')


class laplace_gen(rv_continuous):
    r"""A Laplace continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `laplace` is

    .. math::

        f(x) = \frac{1}{2} \exp(-|x|)

    for a real number :math:`x`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return random_state.laplace(0, 1, size=size)

    def _pdf(self, x):
        # laplace.pdf(x) = 1/2 * exp(-abs(x))
        return 0.5*np.exp(-abs(x))

    def _cdf(self, x):
        with np.errstate(over='ignore'):
            return np.where(x > 0, 1.0 - 0.5*np.exp(-x), 0.5*np.exp(x))

    def _sf(self, x):
        # By symmetry...
        return self._cdf(-x)

    def _ppf(self, q):
        return np.where(q > 0.5, -np.log(2*(1-q)), np.log(2*q))

    def _isf(self, q):
        # By symmetry...
        return -self._ppf(q)

    def _stats(self):
        return 0, 2, 0, 3

    def _entropy(self):
        return np.log(2)+1

    @_call_super_mom
    @replace_notes_in_docstring(rv_continuous, notes="""\
        This function uses explicit formulas for the maximum likelihood
        estimation of the Laplace distribution parameters, so the keyword
        arguments `loc`, `scale`, and `optimizer` are ignored.\n\n""")
    def fit(self, data, *args, **kwds):
        data, floc, fscale = _check_fit_input_parameters(self, data,
                                                         args, kwds)

        # Source: Statistical Distributions, 3rd Edition. Evans, Hastings,
        # and Peacock (2000), Page 124

        if floc is None:
            floc = np.median(data)

        if fscale is None:
            fscale = (np.sum(np.abs(data - floc))) / len(data)

        return floc, fscale


laplace = laplace_gen(name='laplace')


class laplace_asymmetric_gen(rv_continuous):
    r"""An asymmetric Laplace continuous random variable.

    %(before_notes)s

    See Also
    --------
    laplace : Laplace distribution

    Notes
    -----
    The probability density function for `laplace_asymmetric` is

    .. math::

       f(x, \kappa) &= \frac{1}{\kappa+\kappa^{-1}}\exp(-x\kappa),\quad x\ge0\\
                    &= \frac{1}{\kappa+\kappa^{-1}}\exp(x/\kappa),\quad x<0\\

    for :math:`-\infty < x < \infty`, :math:`\kappa > 0`.

    `laplace_asymmetric` takes ``kappa`` as a shape parameter for
    :math:`\kappa`. For :math:`\kappa = 1`, it is identical to a
    Laplace distribution.

    %(after_notes)s

    References
    ----------
    .. [1] "Asymmetric Laplace distribution", Wikipedia
            https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution

    .. [2] Kozubowski TJ and Podgórski K. A Multivariate and
           Asymmetric Generalization of Laplace Distribution,
           Computational Statistics 15, 531--540 (2000).
           :doi:`10.1007/PL00022717`

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("kappa", False, (0, np.inf), (False, False))]

    def _pdf(self, x, kappa):
        return np.exp(self._logpdf(x, kappa))

    def _logpdf(self, x, kappa):
        kapinv = 1/kappa
        lPx = x * np.where(x >= 0, -kappa, kapinv)
        lPx -= np.log(kappa+kapinv)
        return lPx

    def _cdf(self, x, kappa):
        kapinv = 1/kappa
        kappkapinv = kappa+kapinv
        return np.where(x >= 0,
                        1 - np.exp(-x*kappa)*(kapinv/kappkapinv),
                        np.exp(x*kapinv)*(kappa/kappkapinv))

    def _sf(self, x, kappa):
        kapinv = 1/kappa
        kappkapinv = kappa+kapinv
        return np.where(x >= 0,
                        np.exp(-x*kappa)*(kapinv/kappkapinv),
                        1 - np.exp(x*kapinv)*(kappa/kappkapinv))

    def _ppf(self, q, kappa):
        kapinv = 1/kappa
        kappkapinv = kappa+kapinv
        return np.where(q >= kappa/kappkapinv,
                        -np.log((1 - q)*kappkapinv*kappa)*kapinv,
                        np.log(q*kappkapinv/kappa)*kappa)

    def _isf(self, q, kappa):
        kapinv = 1/kappa
        kappkapinv = kappa+kapinv
        return np.where(q <= kapinv/kappkapinv,
                        -np.log(q*kappkapinv*kappa)*kapinv,
                        np.log((1 - q)*kappkapinv/kappa)*kappa)

    def _stats(self, kappa):
        kapinv = 1/kappa
        mn = kapinv - kappa
        var = kapinv*kapinv + kappa*kappa
        g1 = 2.0*(1-np.power(kappa, 6))/np.power(1+np.power(kappa, 4), 1.5)
        g2 = 6.0*(1+np.power(kappa, 8))/np.power(1+np.power(kappa, 4), 2)
        return mn, var, g1, g2

    def _entropy(self, kappa):
        return 1 + np.log(kappa+1/kappa)


laplace_asymmetric = laplace_asymmetric_gen(name='laplace_asymmetric')


def _check_fit_input_parameters(dist, data, args, kwds):
    data = np.asarray(data)
    floc = kwds.get('floc', None)
    fscale = kwds.get('fscale', None)

    num_shapes = len(dist.shapes.split(",")) if dist.shapes else 0
    fshape_keys = []
    fshapes = []

    # user has many options for fixing the shape, so here we standardize it
    # into 'f' + the number of the shape.
    # Adapted from `_reduce_func` in `_distn_infrastructure.py`:
    if dist.shapes:
        shapes = dist.shapes.replace(',', ' ').split()
        for j, s in enumerate(shapes):
            key = 'f' + str(j)
            names = [key, 'f' + s, 'fix_' + s]
            val = _get_fixed_fit_value(kwds, names)
            fshape_keys.append(key)
            fshapes.append(val)
            if val is not None:
                kwds[key] = val

    # determine if there are any unknown arguments in kwds
    known_keys = {'loc', 'scale', 'optimizer', 'method',
                  'floc', 'fscale', *fshape_keys}
    unknown_keys = set(kwds).difference(known_keys)
    if unknown_keys:
        raise TypeError(f"Unknown keyword arguments: {unknown_keys}.")

    if len(args) > num_shapes:
        raise TypeError("Too many positional arguments.")

    if None not in {floc, fscale, *fshapes}:
        # This check is for consistency with `rv_continuous.fit`.
        # Without this check, this function would just return the
        # parameters that were given.
        raise RuntimeError("All parameters fixed. There is nothing to "
                           "optimize.")

    if not np.isfinite(data).all():
        raise ValueError("The data contains non-finite values.")

    return (data, *fshapes, floc, fscale)


class levy_gen(rv_continuous):
    r"""A Levy continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy_stable, levy_l

    Notes
    -----
    The probability density function for `levy` is:

    .. math::

        f(x) = \frac{1}{\sqrt{2\pi x^3}} \exp\left(-\frac{1}{2x}\right)

    for :math:`x >= 0`.

    This is the same as the Levy-stable distribution with :math:`a=1/2` and
    :math:`b=1`.

    %(after_notes)s

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import levy
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 1)

    Calculate the first four moments:

    >>> mean, var, skew, kurt = levy.stats(moments='mvsk')

    Display the probability density function (``pdf``):

    >>> # `levy` is very heavy-tailed.
    >>> # To show a nice plot, let's cut off the upper 40 percent.
    >>> a, b = levy.ppf(0), levy.ppf(0.6)
    >>> x = np.linspace(a, b, 100)
    >>> ax.plot(x, levy.pdf(x),
    ...        'r-', lw=5, alpha=0.6, label='levy pdf')

    Alternatively, the distribution object can be called (as a function)
    to fix the shape, location and scale parameters. This returns a "frozen"
    RV object holding the given parameters fixed.

    Freeze the distribution and display the frozen ``pdf``:

    >>> rv = levy()
    >>> ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

    Check accuracy of ``cdf`` and ``ppf``:

    >>> vals = levy.ppf([0.001, 0.5, 0.999])
    >>> np.allclose([0.001, 0.5, 0.999], levy.cdf(vals))
    True

    Generate random numbers:

    >>> r = levy.rvs(size=1000)

    And compare the histogram:

    >>> # manual binning to ignore the tail
    >>> bins = np.concatenate((np.linspace(a, b, 20), [np.max(r)]))
    >>> ax.hist(r, bins=bins, density=True, histtype='stepfilled', alpha=0.2)
    >>> ax.set_xlim([x[0], x[-1]])
    >>> ax.legend(loc='best', frameon=False)
    >>> plt.show()

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return []

    def _pdf(self, x):
        # levy.pdf(x) = 1 / (x * sqrt(2*pi*x)) * exp(-1/(2*x))
        return 1 / np.sqrt(2*np.pi*x) / x * np.exp(-1/(2*x))

    def _cdf(self, x):
        # Equivalent to 2*norm.sf(np.sqrt(1/x))
        return sc.erfc(np.sqrt(0.5 / x))

    def _sf(self, x):
        return sc.erf(np.sqrt(0.5 / x))

    def _ppf(self, q):
        # Equivalent to 1.0/(norm.isf(q/2)**2) or 0.5/(erfcinv(q)**2)
        val = -sc.ndtri(q/2)
        return 1.0 / (val * val)

    def _isf(self, p):
        return 1/(2*sc.erfinv(p)**2)

    def _stats(self):
        return np.inf, np.inf, np.nan, np.nan


levy = levy_gen(a=0.0, name="levy")


class levy_l_gen(rv_continuous):
    r"""A left-skewed Levy continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy, levy_stable

    Notes
    -----
    The probability density function for `levy_l` is:

    .. math::
        f(x) = \frac{1}{|x| \sqrt{2\pi |x|}} \exp{ \left(-\frac{1}{2|x|} \right)}

    for :math:`x <= 0`.

    This is the same as the Levy-stable distribution with :math:`a=1/2` and
    :math:`b=-1`.

    %(after_notes)s

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import levy_l
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 1)

    Calculate the first four moments:

    >>> mean, var, skew, kurt = levy_l.stats(moments='mvsk')

    Display the probability density function (``pdf``):

    >>> # `levy_l` is very heavy-tailed.
    >>> # To show a nice plot, let's cut off the lower 40 percent.
    >>> a, b = levy_l.ppf(0.4), levy_l.ppf(1)
    >>> x = np.linspace(a, b, 100)
    >>> ax.plot(x, levy_l.pdf(x),
    ...        'r-', lw=5, alpha=0.6, label='levy_l pdf')

    Alternatively, the distribution object can be called (as a function)
    to fix the shape, location and scale parameters. This returns a "frozen"
    RV object holding the given parameters fixed.

    Freeze the distribution and display the frozen ``pdf``:

    >>> rv = levy_l()
    >>> ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

    Check accuracy of ``cdf`` and ``ppf``:

    >>> vals = levy_l.ppf([0.001, 0.5, 0.999])
    >>> np.allclose([0.001, 0.5, 0.999], levy_l.cdf(vals))
    True

    Generate random numbers:

    >>> r = levy_l.rvs(size=1000)

    And compare the histogram:

    >>> # manual binning to ignore the tail
    >>> bins = np.concatenate(([np.min(r)], np.linspace(a, b, 20)))
    >>> ax.hist(r, bins=bins, density=True, histtype='stepfilled', alpha=0.2)
    >>> ax.set_xlim([x[0], x[-1]])
    >>> ax.legend(loc='best', frameon=False)
    >>> plt.show()

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return []

    def _pdf(self, x):
        # levy_l.pdf(x) = 1 / (abs(x) * sqrt(2*pi*abs(x))) * exp(-1/(2*abs(x)))
        ax = abs(x)
        return 1/np.sqrt(2*np.pi*ax)/ax*np.exp(-1/(2*ax))

    def _cdf(self, x):
        ax = abs(x)
        return 2 * _norm_cdf(1 / np.sqrt(ax)) - 1

    def _sf(self, x):
        ax = abs(x)
        return 2 * _norm_sf(1 / np.sqrt(ax))

    def _ppf(self, q):
        val = _norm_ppf((q + 1.0) / 2)
        return -1.0 / (val * val)

    def _isf(self, p):
        return -1/_norm_isf(p/2)**2

    def _stats(self):
        return np.inf, np.inf, np.nan, np.nan


levy_l = levy_l_gen(b=0.0, name="levy_l")


class logistic_gen(rv_continuous):
    r"""A logistic (or Sech-squared) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `logistic` is:

    .. math::

        f(x) = \frac{\exp(-x)}
                    {(1+\exp(-x))^2}

    `logistic` is a special case of `genlogistic` with ``c=1``.

    Remark that the survival function (``logistic.sf``) is equal to the
    Fermi-Dirac distribution describing fermionic statistics.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return random_state.logistic(size=size)

    def _pdf(self, x):
        # logistic.pdf(x) = exp(-x) / (1+exp(-x))**2
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        y = -np.abs(x)
        return y - 2. * sc.log1p(np.exp(y))

    def _cdf(self, x):
        return sc.expit(x)

    def _logcdf(self, x):
        return sc.log_expit(x)

    def _ppf(self, q):
        return sc.logit(q)

    def _sf(self, x):
        return sc.expit(-x)

    def _logsf(self, x):
        return sc.log_expit(-x)

    def _isf(self, q):
        return -sc.logit(q)

    def _stats(self):
        return 0, np.pi*np.pi/3.0, 0, 6.0/5.0

    def _entropy(self):
        # https://en.wikipedia.org/wiki/Logistic_distribution
        return 2.0

    @_call_super_mom
    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        data, floc, fscale = _check_fit_input_parameters(self, data,
                                                         args, kwds)

        # if user has provided `floc` or `fscale`, fall back on super fit
        # method. This scenario is not suitable for solving a system of
        # equations
        if floc is not None or fscale is not None:
            return super().fit(data, *args, **kwds)

        # rv_continuous provided guesses
        loc, scale = self._fitstart(data)
        # account for user provided guesses
        loc = kwds.pop('loc', loc)
        scale = kwds.pop('scale', scale)

        # the maximum likelihood estimators `a` and `b` of the location and
        # scale parameters are roots of the two equations described in `func`.
        # Source: Statistical Distributions, 3rd Edition. Evans, Hastings, and
        # Peacock (2000), Page 130
        def func(params, data):
            a, b = params
            n = len(data)
            c = (data - a) / b
            x1 = np.sum(sc.expit(c)) - n/2
            x2 = np.sum(c*np.tanh(c/2)) - n
            return x1, x2

        return tuple(optimize.root(func, (loc, scale), args=(data,)).x)


logistic = logistic_gen(name='logistic')


class loggamma_gen(rv_continuous):
    r"""A log gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `loggamma` is:

    .. math::

        f(x, c) = \frac{\exp(c x - \exp(x))}
                       {\Gamma(c)}

    for all :math:`x, c > 0`. Here, :math:`\Gamma` is the
    gamma function (`scipy.special.gamma`).

    `loggamma` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _rvs(self, c, size=None, random_state=None):
        # Use the property of the gamma distribution Gamma(c)
        #    Gamma(c) ~ Gamma(c + 1)*U**(1/c),
        # where U is uniform on [0, 1]. (See, e.g.,
        # G. Marsaglia and W.W. Tsang, "A simple method for generating gamma
        # variables", https://doi.org/10.1145/358407.358414)
        # So
        #    log(Gamma(c)) ~ log(Gamma(c + 1)) + log(U)/c
        # Generating a sample with this formulation is a bit slower
        # than the more obvious log(Gamma(c)), but it avoids loss
        # of precision when c << 1.
        return (np.log(random_state.gamma(c + 1, size=size))
                + np.log(random_state.uniform(size=size))/c)

    def _pdf(self, x, c):
        # loggamma.pdf(x, c) = exp(c*x-exp(x)) / gamma(c)
        return np.exp(c*x-np.exp(x)-sc.gammaln(c))

    def _logpdf(self, x, c):
        return c*x - np.exp(x) - sc.gammaln(c)

    def _cdf(self, x, c):
        return sc.gammainc(c, np.exp(x))

    def _ppf(self, q, c):
        return np.log(sc.gammaincinv(c, q))

    def _sf(self, x, c):
        return sc.gammaincc(c, np.exp(x))

    def _isf(self, q, c):
        return np.log(sc.gammainccinv(c, q))

    def _stats(self, c):
        # See, for example, "A Statistical Study of Log-Gamma Distribution", by
        # Ping Shing Chan (thesis, McMaster University, 1993).
        mean = sc.digamma(c)
        var = sc.polygamma(1, c)
        skewness = sc.polygamma(2, c) / np.power(var, 1.5)
        excess_kurtosis = sc.polygamma(3, c) / (var*var)
        return mean, var, skewness, excess_kurtosis


loggamma = loggamma_gen(name='loggamma')


class loglaplace_gen(rv_continuous):
    r"""A log-Laplace continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `loglaplace` is:

    .. math::

        f(x, c) = \begin{cases}\frac{c}{2} x^{ c-1}  &\text{for } 0 < x < 1\\
                               \frac{c}{2} x^{-c-1}  &\text{for } x \ge 1
                  \end{cases}

    for :math:`c > 0`.

    `loglaplace` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    References
    ----------
    T.J. Kozubowski and K. Podgorski, "A log-Laplace growth rate model",
    The Mathematical Scientist, vol. 28, pp. 49-60, 2003.

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # loglaplace.pdf(x, c) = c / 2 * x**(c-1),   for 0 < x < 1
        #                      = c / 2 * x**(-c-1),  for x >= 1
        cd2 = c/2.0
        c = np.where(x < 1, c, -c)
        return cd2*x**(c-1)

    def _cdf(self, x, c):
        return np.where(x < 1, 0.5*x**c, 1-0.5*x**(-c))

    def _ppf(self, q, c):
        return np.where(q < 0.5, (2.0*q)**(1.0/c), (2*(1.0-q))**(-1.0/c))

    def _munp(self, n, c):
        return c**2 / (c**2 - n**2)

    def _entropy(self, c):
        return np.log(2.0/c) + 1.0


loglaplace = loglaplace_gen(a=0.0, name='loglaplace')


def _lognorm_logpdf(x, s):
    return _lazywhere(x != 0, (x, s),
                      lambda x, s: -np.log(x)**2 / (2*s**2) - np.log(s*x*np.sqrt(2*np.pi)),
                      -np.inf)


class lognorm_gen(rv_continuous):
    r"""A lognormal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `lognorm` is:

    .. math::

        f(x, s) = \frac{1}{s x \sqrt{2\pi}}
                  \exp\left(-\frac{\log^2(x)}{2s^2}\right)

    for :math:`x > 0`, :math:`s > 0`.

    `lognorm` takes ``s`` as a shape parameter for :math:`s`.

    %(after_notes)s

    Suppose a normally distributed random variable ``X`` has  mean ``mu`` and
    standard deviation ``sigma``. Then ``Y = exp(X)`` is lognormally
    distributed with ``s = sigma`` and ``scale = exp(mu)``.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return [_ShapeInfo("s", False, (0, np.inf), (False, False))]

    def _rvs(self, s, size=None, random_state=None):
        return np.exp(s * random_state.standard_normal(size))

    def _pdf(self, x, s):
        # lognorm.pdf(x, s) = 1 / (s*x*sqrt(2*pi)) * exp(-1/2*(log(x)/s)**2)
        return np.exp(self._logpdf(x, s))

    def _logpdf(self, x, s):
        return _lognorm_logpdf(x, s)

    def _cdf(self, x, s):
        return _norm_cdf(np.log(x) / s)

    def _logcdf(self, x, s):
        return _norm_logcdf(np.log(x) / s)

    def _ppf(self, q, s):
        return np.exp(s * _norm_ppf(q))

    def _sf(self, x, s):
        return _norm_sf(np.log(x) / s)

    def _logsf(self, x, s):
        return _norm_logsf(np.log(x) / s)

    def _stats(self, s):
        p = np.exp(s*s)
        mu = np.sqrt(p)
        mu2 = p*(p-1)
        g1 = np.sqrt((p-1))*(2+p)
        g2 = np.polyval([1, 2, 3, 0, -6.0], p)
        return mu, mu2, g1, g2

    def _entropy(self, s):
        return 0.5 * (1 + np.log(2*np.pi) + 2 * np.log(s))

    @_call_super_mom
    @extend_notes_in_docstring(rv_continuous, notes="""\
        When `method='MLE'` and
        the location parameter is fixed by using the `floc` argument,
        this function uses explicit formulas for the maximum likelihood
        estimation of the log-normal shape and scale parameters, so the
        `optimizer`, `loc` and `scale` keyword arguments are ignored.
        \n\n""")
    def fit(self, data, *args, **kwds):
        floc = kwds.get('floc', None)
        if floc is None:
            # fall back on the default fit method.
            return super().fit(data, *args, **kwds)

        f0 = (kwds.get('f0', None) or kwds.get('fs', None) or
              kwds.get('fix_s', None))
        fscale = kwds.get('fscale', None)

        if len(args) > 1:
            raise TypeError("Too many input arguments.")
        for name in ['f0', 'fs', 'fix_s', 'floc', 'fscale', 'loc', 'scale',
                     'optimizer', 'method']:
            kwds.pop(name, None)
        if kwds:
            raise TypeError("Unknown arguments: %s." % kwds)

        # Special case: loc is fixed.  Use the maximum likelihood formulas
        # instead of the numerical solver.

        if f0 is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        if not np.isfinite(data).all():
            raise ValueError("The data contains non-finite values.")

        floc = float(floc)
        if floc != 0:
            # Shifting the data by floc. Don't do the subtraction in-place,
            # because `data` might be a view of the input array.
            data = data - floc
        if np.any(data <= 0):
            raise FitDataError("lognorm", lower=floc, upper=np.inf)
        lndata = np.log(data)

        # Three cases to handle:
        # * shape and scale both free
        # * shape fixed, scale free
        # * shape free, scale fixed

        if fscale is None:
            # scale is free.
            scale = np.exp(lndata.mean())
            if f0 is None:
                # shape is free.
                shape = lndata.std()
            else:
                # shape is fixed.
                shape = float(f0)
        else:
            # scale is fixed, shape is free
            scale = float(fscale)
            shape = np.sqrt(((lndata - np.log(scale))**2).mean())

        return shape, floc, scale


lognorm = lognorm_gen(a=0.0, name='lognorm')


class gibrat_gen(rv_continuous):
    r"""A Gibrat continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gibrat` is:

    .. math::

        f(x) = \frac{1}{x \sqrt{2\pi}} \exp(-\frac{1}{2} (\log(x))^2)

    `gibrat` is a special case of `lognorm` with ``s=1``.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return np.exp(random_state.standard_normal(size))

    def _pdf(self, x):
        # gibrat.pdf(x) = 1/(x*sqrt(2*pi)) * exp(-1/2*(log(x))**2)
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return _lognorm_logpdf(x, 1.0)

    def _cdf(self, x):
        return _norm_cdf(np.log(x))

    def _ppf(self, q):
        return np.exp(_norm_ppf(q))

    def _stats(self):
        p = np.e
        mu = np.sqrt(p)
        mu2 = p * (p - 1)
        g1 = np.sqrt((p - 1)) * (2 + p)
        g2 = np.polyval([1, 2, 3, 0, -6.0], p)
        return mu, mu2, g1, g2

    def _entropy(self):
        return 0.5 * np.log(2 * np.pi) + 0.5


# deprecation of gilbrat, see #15911
deprmsg = ("`gilbrat` is a misspelling of the correct name for the `gibrat` "
           "distribution, and will be removed in SciPy 1.11.")


class gilbrat_gen(gibrat_gen):
    # override __call__ protocol from rv_generic to also
    # deprecate instantiation of frozen distributions
    r"""

    .. deprecated:: 1.9.0
        `gilbrat` is deprecated, use `gibrat` instead!
        `gilbrat` is a misspelling of the correct name for the `gibrat`
        distribution, and will be removed in SciPy 1.11.

    """
    def __call__(self, *args, **kwds):
        # align with warning text from np.deprecated that's used for methods
        msg = "`gilbrat` is deprecated, use `gibrat` instead!\n" + deprmsg
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        return self.freeze(*args, **kwds)


gibrat = gibrat_gen(a=0.0, name='gibrat')
gilbrat = gilbrat_gen(a=0.0, name='gilbrat')


# since the deprecated class gets intantiated upon import (and we only want to
# warn upon use), add the deprecation to each (documented) class method, c.f.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gilbrat.html
_gibrat_method_names = [
    "cdf", "entropy", "expect", "fit", "interval", "isf", "logcdf", "logpdf",
    "logsf", "mean", "median", "moment", "pdf", "ppf", "rvs", "sf", "stats",
    "std", "var"
]
for m in _gibrat_method_names:
    wrapper = np.deprecate(getattr(gilbrat, m), f"gilbrat.{m}", f"gibrat.{m}",
                           deprmsg)
    setattr(gilbrat, m, wrapper)


class maxwell_gen(rv_continuous):
    r"""A Maxwell continuous random variable.

    %(before_notes)s

    Notes
    -----
    A special case of a `chi` distribution,  with ``df=3``, ``loc=0.0``,
    and given ``scale = a``, where ``a`` is the parameter used in the
    Mathworld description [1]_.

    The probability density function for `maxwell` is:

    .. math::

        f(x) = \sqrt{2/\pi}x^2 \exp(-x^2/2)

    for :math:`x >= 0`.

    %(after_notes)s

    References
    ----------
    .. [1] http://mathworld.wolfram.com/MaxwellDistribution.html

    %(example)s
    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return chi.rvs(3.0, size=size, random_state=random_state)

    def _pdf(self, x):
        # maxwell.pdf(x) = sqrt(2/pi)x**2 * exp(-x**2/2)
        return _SQRT_2_OVER_PI*x*x*np.exp(-x*x/2.0)

    def _logpdf(self, x):
        # Allow x=0 without 'divide by zero' warnings
        with np.errstate(divide='ignore'):
            return _LOG_SQRT_2_OVER_PI + 2*np.log(x) - 0.5*x*x

    def _cdf(self, x):
        return sc.gammainc(1.5, x*x/2.0)

    def _ppf(self, q):
        return np.sqrt(2*sc.gammaincinv(1.5, q))

    def _stats(self):
        val = 3*np.pi-8
        return (2*np.sqrt(2.0/np.pi),
                3-8/np.pi,
                np.sqrt(2)*(32-10*np.pi)/val**1.5,
                (-12*np.pi*np.pi + 160*np.pi - 384) / val**2.0)

    def _entropy(self):
        return _EULER + 0.5*np.log(2*np.pi)-0.5


maxwell = maxwell_gen(a=0.0, name='maxwell')


class mielke_gen(rv_continuous):
    r"""A Mielke Beta-Kappa / Dagum continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `mielke` is:

    .. math::

        f(x, k, s) = \frac{k x^{k-1}}{(1+x^s)^{1+k/s}}

    for :math:`x > 0` and :math:`k, s > 0`. The distribution is sometimes
    called Dagum distribution ([2]_). It was already defined in [3]_, called
    a Burr Type III distribution (`burr` with parameters ``c=s`` and
    ``d=k/s``).

    `mielke` takes ``k`` and ``s`` as shape parameters.

    %(after_notes)s

    References
    ----------
    .. [1] Mielke, P.W., 1973 "Another Family of Distributions for Describing
           and Analyzing Precipitation Data." J. Appl. Meteor., 12, 275-280
    .. [2] Dagum, C., 1977 "A new model for personal income distribution."
           Economie Appliquee, 33, 327-367.
    .. [3] Burr, I. W. "Cumulative frequency functions", Annals of
           Mathematical Statistics, 13(2), pp 215-232 (1942).

    %(example)s

    """
    def _shape_info(self):
        ik = _ShapeInfo("k", False, (0, np.inf), (False, False))
        i_s = _ShapeInfo("s", False, (0, np.inf), (False, False))
        return [ik, i_s]

    def _pdf(self, x, k, s):
        return k*x**(k-1.0) / (1.0+x**s)**(1.0+k*1.0/s)

    def _logpdf(self, x, k, s):
        # Allow x=0 without 'divide by zero' warnings.
        with np.errstate(divide='ignore'):
            return np.log(k) + np.log(x)*(k - 1) - np.log1p(x**s)*(1 + k/s)

    def _cdf(self, x, k, s):
        return x**k / (1.0+x**s)**(k*1.0/s)

    def _ppf(self, q, k, s):
        qsk = pow(q, s*1.0/k)
        return pow(qsk/(1.0-qsk), 1.0/s)

    def _munp(self, n, k, s):
        def nth_moment(n, k, s):
            # n-th moment is defined for -k < n < s
            return sc.gamma((k+n)/s)*sc.gamma(1-n/s)/sc.gamma(k/s)

        return _lazywhere(n < s, (n, k, s), nth_moment, np.inf)


mielke = mielke_gen(a=0.0, name='mielke')


class kappa4_gen(rv_continuous):
    r"""Kappa 4 parameter distribution.

    %(before_notes)s

    Notes
    -----
    The probability density function for kappa4 is:

    .. math::

        f(x, h, k) = (1 - k x)^{1/k - 1} (1 - h (1 - k x)^{1/k})^{1/h-1}

    if :math:`h` and :math:`k` are not equal to 0.

    If :math:`h` or :math:`k` are zero then the pdf can be simplified:

    h = 0 and k != 0::

        kappa4.pdf(x, h, k) = (1.0 - k*x)**(1.0/k - 1.0)*
                              exp(-(1.0 - k*x)**(1.0/k))

    h != 0 and k = 0::

        kappa4.pdf(x, h, k) = exp(-x)*(1.0 - h*exp(-x))**(1.0/h - 1.0)

    h = 0 and k = 0::

        kappa4.pdf(x, h, k) = exp(-x)*exp(-exp(-x))

    kappa4 takes :math:`h` and :math:`k` as shape parameters.

    The kappa4 distribution returns other distributions when certain
    :math:`h` and :math:`k` values are used.

    +------+-------------+----------------+------------------+
    | h    | k=0.0       | k=1.0          | -inf<=k<=inf     |
    +======+=============+================+==================+
    | -1.0 | Logistic    |                | Generalized      |
    |      |             |                | Logistic(1)      |
    |      |             |                |                  |
    |      | logistic(x) |                |                  |
    +------+-------------+----------------+------------------+
    |  0.0 | Gumbel      | Reverse        | Generalized      |
    |      |             | Exponential(2) | Extreme Value    |
    |      |             |                |                  |
    |      | gumbel_r(x) |                | genextreme(x, k) |
    +------+-------------+----------------+------------------+
    |  1.0 | Exponential | Uniform        | Generalized      |
    |      |             |                | Pareto           |
    |      |             |                |                  |
    |      | expon(x)    | uniform(x)     | genpareto(x, -k) |
    +------+-------------+----------------+------------------+

    (1) There are at least five generalized logistic distributions.
        Four are described here:
        https://en.wikipedia.org/wiki/Generalized_logistic_distribution
        The "fifth" one is the one kappa4 should match which currently
        isn't implemented in scipy:
        https://en.wikipedia.org/wiki/Talk:Generalized_logistic_distribution
        https://www.mathwave.com/help/easyfit/html/analyses/distributions/gen_logistic.html
    (2) This distribution is currently not in scipy.

    References
    ----------
    J.C. Finney, "Optimization of a Skewed Logistic Distribution With Respect
    to the Kolmogorov-Smirnov Test", A Dissertation Submitted to the Graduate
    Faculty of the Louisiana State University and Agricultural and Mechanical
    College, (August, 2004),
    https://digitalcommons.lsu.edu/gradschool_dissertations/3672

    J.R.M. Hosking, "The four-parameter kappa distribution". IBM J. Res.
    Develop. 38 (3), 25 1-258 (1994).

    B. Kumphon, A. Kaew-Man, P. Seenoi, "A Rainfall Distribution for the Lampao
    Site in the Chi River Basin, Thailand", Journal of Water Resource and
    Protection, vol. 4, 866-869, (2012).
    :doi:`10.4236/jwarp.2012.410101`

    C. Winchester, "On Estimation of the Four-Parameter Kappa Distribution", A
    Thesis Submitted to Dalhousie University, Halifax, Nova Scotia, (March
    2000).
    http://www.nlc-bnc.ca/obj/s4/f2/dsk2/ftp01/MQ57336.pdf

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, h, k):
        shape = np.broadcast_arrays(h, k)[0].shape
        return np.full(shape, fill_value=True)

    def _shape_info(self):
        ih = _ShapeInfo("h", False, (-np.inf, np.inf), (False, False))
        ik = _ShapeInfo("k", False, (-np.inf, np.inf), (False, False))
        return [ih, ik]

    def _get_support(self, h, k):
        condlist = [np.logical_and(h > 0, k > 0),
                    np.logical_and(h > 0, k == 0),
                    np.logical_and(h > 0, k < 0),
                    np.logical_and(h <= 0, k > 0),
                    np.logical_and(h <= 0, k == 0),
                    np.logical_and(h <= 0, k < 0)]

        def f0(h, k):
            return (1.0 - np.float_power(h, -k))/k

        def f1(h, k):
            return np.log(h)

        def f3(h, k):
            a = np.empty(np.shape(h))
            a[:] = -np.inf
            return a

        def f5(h, k):
            return 1.0/k

        _a = _lazyselect(condlist,
                             [f0, f1, f0, f3, f3, f5],
                             [h, k],
                             default=np.nan)

        def f0(h, k):
            return 1.0/k

        def f1(h, k):
            a = np.empty(np.shape(h))
            a[:] = np.inf
            return a

        _b = _lazyselect(condlist,
                             [f0, f1, f1, f0, f1, f1],
                             [h, k],
                             default=np.nan)
        return _a, _b

    def _pdf(self, x, h, k):
        # kappa4.pdf(x, h, k) = (1.0 - k*x)**(1.0/k - 1.0)*
        #                       (1.0 - h*(1.0 - k*x)**(1.0/k))**(1.0/h-1)
        return np.exp(self._logpdf(x, h, k))

    def _logpdf(self, x, h, k):
        condlist = [np.logical_and(h != 0, k != 0),
                    np.logical_and(h == 0, k != 0),
                    np.logical_and(h != 0, k == 0),
                    np.logical_and(h == 0, k == 0)]

        def f0(x, h, k):
            '''pdf = (1.0 - k*x)**(1.0/k - 1.0)*(
                      1.0 - h*(1.0 - k*x)**(1.0/k))**(1.0/h-1.0)
               logpdf = ...
            '''
            return (sc.xlog1py(1.0/k - 1.0, -k*x) +
                    sc.xlog1py(1.0/h - 1.0, -h*(1.0 - k*x)**(1.0/k)))

        def f1(x, h, k):
            '''pdf = (1.0 - k*x)**(1.0/k - 1.0)*np.exp(-(
                      1.0 - k*x)**(1.0/k))
               logpdf = ...
            '''
            return sc.xlog1py(1.0/k - 1.0, -k*x) - (1.0 - k*x)**(1.0/k)

        def f2(x, h, k):
            '''pdf = np.exp(-x)*(1.0 - h*np.exp(-x))**(1.0/h - 1.0)
               logpdf = ...
            '''
            return -x + sc.xlog1py(1.0/h - 1.0, -h*np.exp(-x))

        def f3(x, h, k):
            '''pdf = np.exp(-x-np.exp(-x))
               logpdf = ...
            '''
            return -x - np.exp(-x)

        return _lazyselect(condlist,
                           [f0, f1, f2, f3],
                           [x, h, k],
                           default=np.nan)

    def _cdf(self, x, h, k):
        return np.exp(self._logcdf(x, h, k))

    def _logcdf(self, x, h, k):
        condlist = [np.logical_and(h != 0, k != 0),
                    np.logical_and(h == 0, k != 0),
                    np.logical_and(h != 0, k == 0),
                    np.logical_and(h == 0, k == 0)]

        def f0(x, h, k):
            '''cdf = (1.0 - h*(1.0 - k*x)**(1.0/k))**(1.0/h)
               logcdf = ...
            '''
            return (1.0/h)*sc.log1p(-h*(1.0 - k*x)**(1.0/k))

        def f1(x, h, k):
            '''cdf = np.exp(-(1.0 - k*x)**(1.0/k))
               logcdf = ...
            '''
            return -(1.0 - k*x)**(1.0/k)

        def f2(x, h, k):
            '''cdf = (1.0 - h*np.exp(-x))**(1.0/h)
               logcdf = ...
            '''
            return (1.0/h)*sc.log1p(-h*np.exp(-x))

        def f3(x, h, k):
            '''cdf = np.exp(-np.exp(-x))
               logcdf = ...
            '''
            return -np.exp(-x)

        return _lazyselect(condlist,
                           [f0, f1, f2, f3],
                           [x, h, k],
                           default=np.nan)

    def _ppf(self, q, h, k):
        condlist = [np.logical_and(h != 0, k != 0),
                    np.logical_and(h == 0, k != 0),
                    np.logical_and(h != 0, k == 0),
                    np.logical_and(h == 0, k == 0)]

        def f0(q, h, k):
            return 1.0/k*(1.0 - ((1.0 - (q**h))/h)**k)

        def f1(q, h, k):
            return 1.0/k*(1.0 - (-np.log(q))**k)

        def f2(q, h, k):
            '''ppf = -np.log((1.0 - (q**h))/h)
            '''
            return -sc.log1p(-(q**h)) + np.log(h)

        def f3(q, h, k):
            return -np.log(-np.log(q))

        return _lazyselect(condlist,
                           [f0, f1, f2, f3],
                           [q, h, k],
                           default=np.nan)

    def _get_stats_info(self, h, k):
        condlist = [
            np.logical_and(h < 0, k >= 0),
            k < 0,
        ]

        def f0(h, k):
            return (-1.0/h*k).astype(int)

        def f1(h, k):
            return (-1.0/k).astype(int)

        return _lazyselect(condlist, [f0, f1], [h, k], default=5)

    def _stats(self, h, k):
        maxr = self._get_stats_info(h, k)
        outputs = [None if np.any(r < maxr) else np.nan for r in range(1, 5)]
        return outputs[:]

    def _mom1_sc(self, m, *args):
        maxr = self._get_stats_info(args[0], args[1])
        if m >= maxr:
            return np.nan
        return integrate.quad(self._mom_integ1, 0, 1, args=(m,)+args)[0]


kappa4 = kappa4_gen(name='kappa4')


class kappa3_gen(rv_continuous):
    r"""Kappa 3 parameter distribution.

    %(before_notes)s

    Notes
    -----
    The probability density function for `kappa3` is:

    .. math::

        f(x, a) = a (a + x^a)^{-(a + 1)/a}

    for :math:`x > 0` and :math:`a > 0`.

    `kappa3` takes ``a`` as a shape parameter for :math:`a`.

    References
    ----------
    P.W. Mielke and E.S. Johnson, "Three-Parameter Kappa Distribution Maximum
    Likelihood and Likelihood Ratio Tests", Methods in Weather Research,
    701-707, (September, 1973),
    :doi:`10.1175/1520-0493(1973)101<0701:TKDMLE>2.3.CO;2`

    B. Kumphon, "Maximum Entropy and Maximum Likelihood Estimation for the
    Three-Parameter Kappa Distribution", Open Journal of Statistics, vol 2,
    415-419 (2012), :doi:`10.4236/ojs.2012.24050`

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("a", False, (0, np.inf), (False, False))]

    def _pdf(self, x, a):
        # kappa3.pdf(x, a) = a*(a + x**a)**(-(a + 1)/a),     for x > 0
        return a*(a + x**a)**(-1.0/a-1)

    def _cdf(self, x, a):
        return x*(a + x**a)**(-1.0/a)

    def _ppf(self, q, a):
        return (a/(q**-a - 1.0))**(1.0/a)

    def _stats(self, a):
        outputs = [None if np.any(i < a) else np.nan for i in range(1, 5)]
        return outputs[:]

    def _mom1_sc(self, m, *args):
        if np.any(m >= args[0]):
            return np.nan
        return integrate.quad(self._mom_integ1, 0, 1, args=(m,)+args)[0]


kappa3 = kappa3_gen(a=0.0, name='kappa3')


class moyal_gen(rv_continuous):
    r"""A Moyal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `moyal` is:

    .. math::

        f(x) = \exp(-(x + \exp(-x))/2) / \sqrt{2\pi}

    for a real number :math:`x`.

    %(after_notes)s

    This distribution has utility in high-energy physics and radiation
    detection. It describes the energy loss of a charged relativistic
    particle due to ionization of the medium [1]_. It also provides an
    approximation for the Landau distribution. For an in depth description
    see [2]_. For additional description, see [3]_.

    References
    ----------
    .. [1] J.E. Moyal, "XXX. Theory of ionization fluctuations",
           The London, Edinburgh, and Dublin Philosophical Magazine
           and Journal of Science, vol 46, 263-280, (1955).
           :doi:`10.1080/14786440308521076` (gated)
    .. [2] G. Cordeiro et al., "The beta Moyal: a useful skew distribution",
           International Journal of Research and Reviews in Applied Sciences,
           vol 10, 171-192, (2012).
           http://www.arpapress.com/Volumes/Vol10Issue2/IJRRAS_10_2_02.pdf
    .. [3] C. Walck, "Handbook on Statistical Distributions for
           Experimentalists; International Report SUF-PFY/96-01", Chapter 26,
           University of Stockholm: Stockholm, Sweden, (2007).
           http://www.stat.rice.edu/~dobelman/textfiles/DistributionsHandbook.pdf

    .. versionadded:: 1.1.0

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        u1 = gamma.rvs(a=0.5, scale=2, size=size,
                       random_state=random_state)
        return -np.log(u1)

    def _pdf(self, x):
        return np.exp(-0.5 * (x + np.exp(-x))) / np.sqrt(2*np.pi)

    def _cdf(self, x):
        return sc.erfc(np.exp(-0.5 * x) / np.sqrt(2))

    def _sf(self, x):
        return sc.erf(np.exp(-0.5 * x) / np.sqrt(2))

    def _ppf(self, x):
        return -np.log(2 * sc.erfcinv(x)**2)

    def _stats(self):
        mu = np.log(2) + np.euler_gamma
        mu2 = np.pi**2 / 2
        g1 = 28 * np.sqrt(2) * sc.zeta(3) / np.pi**3
        g2 = 4.
        return mu, mu2, g1, g2

    def _munp(self, n):
        if n == 1.0:
            return np.log(2) + np.euler_gamma
        elif n == 2.0:
            return np.pi**2 / 2 + (np.log(2) + np.euler_gamma)**2
        elif n == 3.0:
            tmp1 = 1.5 * np.pi**2 * (np.log(2)+np.euler_gamma)
            tmp2 = (np.log(2)+np.euler_gamma)**3
            tmp3 = 14 * sc.zeta(3)
            return tmp1 + tmp2 + tmp3
        elif n == 4.0:
            tmp1 = 4 * 14 * sc.zeta(3) * (np.log(2) + np.euler_gamma)
            tmp2 = 3 * np.pi**2 * (np.log(2) + np.euler_gamma)**2
            tmp3 = (np.log(2) + np.euler_gamma)**4
            tmp4 = 7 * np.pi**4 / 4
            return tmp1 + tmp2 + tmp3 + tmp4
        else:
            # return generic for higher moments
            # return rv_continuous._mom1_sc(self, n, b)
            return self._mom1_sc(n)


moyal = moyal_gen(name="moyal")


class nakagami_gen(rv_continuous):
    r"""A Nakagami continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `nakagami` is:

    .. math::

        f(x, \nu) = \frac{2 \nu^\nu}{\Gamma(\nu)} x^{2\nu-1} \exp(-\nu x^2)

    for :math:`x >= 0`, :math:`\nu > 0`. The distribution was introduced in
    [2]_, see also [1]_ for further information.

    `nakagami` takes ``nu`` as a shape parameter for :math:`\nu`.

    %(after_notes)s

    References
    ----------
    .. [1] "Nakagami distribution", Wikipedia
           https://en.wikipedia.org/wiki/Nakagami_distribution
    .. [2] M. Nakagami, "The m-distribution - A general formula of intensity
           distribution of rapid fading", Statistical methods in radio wave
           propagation, Pergamon Press, 1960, 3-36.
           :doi:`10.1016/B978-0-08-009306-2.50005-4`

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("nu", False, (0, np.inf), (False, False))]

    def _pdf(self, x, nu):
        return np.exp(self._logpdf(x, nu))

    def _logpdf(self, x, nu):
        # nakagami.pdf(x, nu) = 2 * nu**nu / gamma(nu) *
        #                       x**(2*nu-1) * exp(-nu*x**2)
        return (np.log(2) + sc.xlogy(nu, nu) - sc.gammaln(nu) +
                sc.xlogy(2*nu - 1, x) - nu*x**2)

    def _cdf(self, x, nu):
        return sc.gammainc(nu, nu*x*x)

    def _ppf(self, q, nu):
        return np.sqrt(1.0/nu*sc.gammaincinv(nu, q))

    def _sf(self, x, nu):
        return sc.gammaincc(nu, nu*x*x)

    def _isf(self, p, nu):
        return np.sqrt(1/nu * sc.gammainccinv(nu, p))

    def _stats(self, nu):
        mu = sc.gamma(nu+0.5)/sc.gamma(nu)/np.sqrt(nu)
        mu2 = 1.0-mu*mu
        g1 = mu * (1 - 4*nu*mu2) / 2.0 / nu / np.power(mu2, 1.5)
        g2 = -6*mu**4*nu + (8*nu-2)*mu**2-2*nu + 1
        g2 /= nu*mu2**2.0
        return mu, mu2, g1, g2

    def _rvs(self, nu, size=None, random_state=None):
        # this relationship can be found in [1] or by a direct calculation
        return np.sqrt(random_state.standard_gamma(nu, size=size) / nu)

    def _fitstart(self, data, args=None):
        if args is None:
            args = (1.0,) * self.numargs
        # Analytical justified estimates
        # see: https://docs.scipy.org/doc/scipy/reference/tutorial/stats/continuous_nakagami.html
        loc = np.min(data)
        scale = np.sqrt(np.sum((data - loc)**2) / len(data))
        return args + (loc, scale)


nakagami = nakagami_gen(a=0.0, name="nakagami")


# The function name ncx2 is an abbreviation for noncentral chi squared.
def _ncx2_log_pdf(x, df, nc):
    # We use (xs**2 + ns**2)/2 = (xs - ns)**2/2  + xs*ns, and include the
    # factor of exp(-xs*ns) into the ive function to improve numerical
    # stability at large values of xs. See also `rice.pdf`.
    df2 = df/2.0 - 1.0
    xs, ns = np.sqrt(x), np.sqrt(nc)
    res = sc.xlogy(df2/2.0, x/nc) - 0.5*(xs - ns)**2
    corr = sc.ive(df2, xs*ns) / 2.0
    # Return res + np.log(corr) avoiding np.log(0)
    return _lazywhere(
        corr > 0,
        (res, corr),
        f=lambda r, c: r + np.log(c),
        fillvalue=-np.inf)


class ncx2_gen(rv_continuous):
    r"""A non-central chi-squared continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `ncx2` is:

    .. math::

        f(x, k, \lambda) = \frac{1}{2} \exp(-(\lambda+x)/2)
            (x/\lambda)^{(k-2)/4}  I_{(k-2)/2}(\sqrt{\lambda x})

    for :math:`x >= 0`, :math:`k > 0` and :math:`\lambda \ge 0`.
    :math:`k` specifies the degrees of freedom (denoted ``df`` in the
    implementation) and :math:`\lambda` is the non-centrality parameter
    (denoted ``nc`` in the implementation). :math:`I_\nu` denotes the
    modified Bessel function of first order of degree :math:`\nu`
    (`scipy.special.iv`).

    `ncx2` takes ``df`` and ``nc`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, df, nc):
        return (df > 0) & np.isfinite(df) & (nc >= 0)

    def _shape_info(self):
        idf = _ShapeInfo("df", False, (0, np.inf), (False, False))
        inc = _ShapeInfo("nc", False, (0, np.inf), (True, False))
        return [idf, inc]

    def _rvs(self, df, nc, size=None, random_state=None):
        return random_state.noncentral_chisquare(df, nc, size)

    def _logpdf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_ncx2_log_pdf,
                          f2=lambda x, df, _: chi2._logpdf(x, df))

    def _pdf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        with warnings.catch_warnings():
            message = "overflow encountered in _ncx2_pdf"
            warnings.filterwarnings("ignore", message=message)
            return _lazywhere(cond, (x, df, nc), f=_boost._ncx2_pdf,
                              f2=lambda x, df, _: chi2._pdf(x, df))

    def _cdf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_boost._ncx2_cdf,
                          f2=lambda x, df, _: chi2._cdf(x, df))

    def _ppf(self, q, df, nc):
        cond = np.ones_like(q, dtype=bool) & (nc != 0)
        with warnings.catch_warnings():
            message = "overflow encountered in _ncx2_ppf"
            warnings.filterwarnings("ignore", message=message)
            return _lazywhere(cond, (q, df, nc), f=_boost._ncx2_ppf,
                              f2=lambda x, df, _: chi2._ppf(x, df))

    def _sf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_boost._ncx2_sf,
                          f2=lambda x, df, _: chi2._sf(x, df))

    def _isf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        with warnings.catch_warnings():
            message = "overflow encountered in _ncx2_isf"
            warnings.filterwarnings("ignore", message=message)
            return _lazywhere(cond, (x, df, nc), f=_boost._ncx2_isf,
                              f2=lambda x, df, _: chi2._isf(x, df))

    def _stats(self, df, nc):
        return (
            _boost._ncx2_mean(df, nc),
            _boost._ncx2_variance(df, nc),
            _boost._ncx2_skewness(df, nc),
            _boost._ncx2_kurtosis_excess(df, nc),
        )


ncx2 = ncx2_gen(a=0.0, name='ncx2')


class ncf_gen(rv_continuous):
    r"""A non-central F distribution continuous random variable.

    %(before_notes)s

    See Also
    --------
    scipy.stats.f : Fisher distribution

    Notes
    -----
    The probability density function for `ncf` is:

    .. math::

        f(x, n_1, n_2, \lambda) =
            \exp\left(\frac{\lambda}{2} +
                      \lambda n_1 \frac{x}{2(n_1 x + n_2)}
                \right)
            n_1^{n_1/2} n_2^{n_2/2} x^{n_1/2 - 1} \\
            (n_2 + n_1 x)^{-(n_1 + n_2)/2}
            \gamma(n_1/2) \gamma(1 + n_2/2) \\
            \frac{L^{\frac{n_1}{2}-1}_{n_2/2}
                \left(-\lambda n_1 \frac{x}{2(n_1 x + n_2)}\right)}
            {B(n_1/2, n_2/2)
                \gamma\left(\frac{n_1 + n_2}{2}\right)}

    for :math:`n_1, n_2 > 0`, :math:`\lambda \ge 0`.  Here :math:`n_1` is the
    degrees of freedom in the numerator, :math:`n_2` the degrees of freedom in
    the denominator, :math:`\lambda` the non-centrality parameter,
    :math:`\gamma` is the logarithm of the Gamma function, :math:`L_n^k` is a
    generalized Laguerre polynomial and :math:`B` is the beta function.

    `ncf` takes ``df1``, ``df2`` and ``nc`` as shape parameters. If ``nc=0``,
    the distribution becomes equivalent to the Fisher distribution.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, df1, df2, nc):
        return (df1 > 0) & (df2 > 0) & (nc >= 0)

    def _shape_info(self):
        idf1 = _ShapeInfo("df1", False, (0, np.inf), (False, False))
        idf2 = _ShapeInfo("df2", False, (0, np.inf), (False, False))
        inc = _ShapeInfo("nc", False, (0, np.inf), (True, False))
        return [idf1, idf2, inc]

    def _rvs(self, dfn, dfd, nc, size=None, random_state=None):
        return random_state.noncentral_f(dfn, dfd, nc, size)

    def _pdf(self, x, dfn, dfd, nc):
        # ncf.pdf(x, df1, df2, nc) = exp(nc/2 + nc*df1*x/(2*(df1*x+df2))) *
        #             df1**(df1/2) * df2**(df2/2) * x**(df1/2-1) *
        #             (df2+df1*x)**(-(df1+df2)/2) *
        #             gamma(df1/2)*gamma(1+df2/2) *
        #             L^{v1/2-1}^{v2/2}(-nc*v1*x/(2*(v1*x+v2))) /
        #             (B(v1/2, v2/2) * gamma((v1+v2)/2))
        return _boost._ncf_pdf(x, dfn, dfd, nc)

    def _cdf(self, x, dfn, dfd, nc):
        return _boost._ncf_cdf(x, dfn, dfd, nc)

    def _ppf(self, q, dfn, dfd, nc):
        return _boost._ncf_ppf(q, dfn, dfd, nc)

    def _sf(self, x, dfn, dfd, nc):
        return _boost._ncf_sf(x, dfn, dfd, nc)

    def _isf(self, x, dfn, dfd, nc):
        return _boost._ncf_isf(x, dfn, dfd, nc)

    def _munp(self, n, dfn, dfd, nc):
        val = (dfn * 1.0/dfd)**n
        term = sc.gammaln(n+0.5*dfn) + sc.gammaln(0.5*dfd-n) - sc.gammaln(dfd*0.5)
        val *= np.exp(-nc / 2.0+term)
        val *= sc.hyp1f1(n+0.5*dfn, 0.5*dfn, 0.5*nc)
        return val

    def _stats(self, dfn, dfd, nc, moments='mv'):
        mu = _boost._ncf_mean(dfn, dfd, nc)
        mu2 = _boost._ncf_variance(dfn, dfd, nc)
        g1 = _boost._ncf_skewness(dfn, dfd, nc) if 's' in moments else None
        g2 = _boost._ncf_kurtosis_excess(
            dfn, dfd, nc) if 'k' in moments else None
        return mu, mu2, g1, g2


ncf = ncf_gen(a=0.0, name='ncf')


class t_gen(rv_continuous):
    r"""A Student's t continuous random variable.

    For the noncentral t distribution, see `nct`.

    %(before_notes)s

    See Also
    --------
    nct

    Notes
    -----
    The probability density function for `t` is:

    .. math::

        f(x, \nu) = \frac{\Gamma((\nu+1)/2)}
                        {\sqrt{\pi \nu} \Gamma(\nu/2)}
                    (1+x^2/\nu)^{-(\nu+1)/2}

    where :math:`x` is a real number and the degrees of freedom parameter
    :math:`\nu` (denoted ``df`` in the implementation) satisfies
    :math:`\nu > 0`. :math:`\Gamma` is the gamma function
    (`scipy.special.gamma`).

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("df", False, (0, np.inf), (False, False))]

    def _rvs(self, df, size=None, random_state=None):
        return random_state.standard_t(df, size=size)

    def _pdf(self, x, df):
        return _lazywhere(
            df == np.inf, (x, df),
            f=lambda x, df: norm._pdf(x),
            f2=lambda x, df: (
                np.exp(sc.gammaln((df+1)/2)-sc.gammaln(df/2))
                / (np.sqrt(df*np.pi)*(1+(x**2)/df)**((df+1)/2))
            )
        )

    def _logpdf(self, x, df):
        return _lazywhere(
            df == np.inf, (x, df),
            f=lambda x, df: norm._logpdf(x),
            f2=lambda x, df: (
                sc.gammaln((df+1)/2) - sc.gammaln(df/2)
                - (0.5*np.log(df*np.pi)
                   + (df+1)/2*np.log(1+(x**2)/df))
            )
        )

    def _cdf(self, x, df):
        return sc.stdtr(df, x)

    def _sf(self, x, df):
        return sc.stdtr(df, -x)

    def _ppf(self, q, df):
        return sc.stdtrit(df, q)

    def _isf(self, q, df):
        return -sc.stdtrit(df, q)

    def _stats(self, df):
        # infinite df -> normal distribution (0.0, 1.0, 0.0, 0.0)
        infinite_df = np.isposinf(df)

        mu = np.where(df > 1, 0.0, np.inf)

        condlist = ((df > 1) & (df <= 2),
                    (df > 2) & np.isfinite(df),
                    infinite_df)
        choicelist = (lambda df: np.broadcast_to(np.inf, df.shape),
                      lambda df: df / (df-2.0),
                      lambda df: np.broadcast_to(1, df.shape))
        mu2 = _lazyselect(condlist, choicelist, (df,), np.nan)

        g1 = np.where(df > 3, 0.0, np.nan)

        condlist = ((df > 2) & (df <= 4),
                    (df > 4) & np.isfinite(df),
                    infinite_df)
        choicelist = (lambda df: np.broadcast_to(np.inf, df.shape),
                      lambda df: 6.0 / (df-4.0),
                      lambda df: np.broadcast_to(0, df.shape))
        g2 = _lazyselect(condlist, choicelist, (df,), np.nan)

        return mu, mu2, g1, g2

    def _entropy(self, df):
        if df == np.inf:
            return norm._entropy()
        half = df/2
        half1 = (df + 1)/2
        return (half1*(sc.digamma(half1) - sc.digamma(half))
                + np.log(np.sqrt(df)*sc.beta(half, 0.5)))


t = t_gen(name='t')


class nct_gen(rv_continuous):
    r"""A non-central Student's t continuous random variable.

    %(before_notes)s

    Notes
    -----
    If :math:`Y` is a standard normal random variable and :math:`V` is
    an independent chi-square random variable (`chi2`) with :math:`k` degrees
    of freedom, then

    .. math::

        X = \frac{Y + c}{\sqrt{V/k}}

    has a non-central Student's t distribution on the real line.
    The degrees of freedom parameter :math:`k` (denoted ``df`` in the
    implementation) satisfies :math:`k > 0` and the noncentrality parameter
    :math:`c` (denoted ``nc`` in the implementation) is a real number.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, df, nc):
        return (df > 0) & (nc == nc)

    def _shape_info(self):
        idf = _ShapeInfo("df", False, (0, np.inf), (False, False))
        inc = _ShapeInfo("nc", False, (-np.inf, np.inf), (False, False))
        return [idf, inc]

    def _rvs(self, df, nc, size=None, random_state=None):
        n = norm.rvs(loc=nc, size=size, random_state=random_state)
        c2 = chi2.rvs(df, size=size, random_state=random_state)
        return n * np.sqrt(df) / np.sqrt(c2)

    def _pdf(self, x, df, nc):
        # Boost version has accuracy issues in left tail; see gh-16591
        n = df*1.0
        nc = nc*1.0
        x2 = x*x
        ncx2 = nc*nc*x2
        fac1 = n + x2
        trm1 = (n/2.*np.log(n) + sc.gammaln(n+1)
                - (n*np.log(2) + nc*nc/2 + (n/2)*np.log(fac1)
                   + sc.gammaln(n/2)))
        Px = np.exp(trm1)
        valF = ncx2 / (2*fac1)
        trm1 = (np.sqrt(2)*nc*x*sc.hyp1f1(n/2+1, 1.5, valF)
                / np.asarray(fac1*sc.gamma((n+1)/2)))
        trm2 = (sc.hyp1f1((n+1)/2, 0.5, valF)
                / np.asarray(np.sqrt(fac1)*sc.gamma(n/2+1)))
        Px *= trm1+trm2
        return np.clip(Px, 0, None)

    def _cdf(self, x, df, nc):
        return np.clip(_boost._nct_cdf(x, df, nc), 0, 1)

    def _ppf(self, q, df, nc):
        return _boost._nct_ppf(q, df, nc)

    def _sf(self, x, df, nc):
        return np.clip(_boost._nct_sf(x, df, nc), 0, 1)

    def _isf(self, x, df, nc):
        return _boost._nct_isf(x, df, nc)

    def _stats(self, df, nc, moments='mv'):
        mu = _boost._nct_mean(df, nc)
        mu2 = _boost._nct_variance(df, nc)
        g1 = _boost._nct_skewness(df, nc) if 's' in moments else None
        g2 = _boost._nct_kurtosis_excess(df, nc)-3 if 'k' in moments else None
        return mu, mu2, g1, g2


nct = nct_gen(name="nct")


class pareto_gen(rv_continuous):
    r"""A Pareto continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `pareto` is:

    .. math::

        f(x, b) = \frac{b}{x^{b+1}}

    for :math:`x \ge 1`, :math:`b > 0`.

    `pareto` takes ``b`` as a shape parameter for :math:`b`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("b", False, (0, np.inf), (False, False))]

    def _pdf(self, x, b):
        # pareto.pdf(x, b) = b / x**(b+1)
        return b * x**(-b-1)

    def _cdf(self, x, b):
        return 1 - x**(-b)

    def _ppf(self, q, b):
        return pow(1-q, -1.0/b)

    def _sf(self, x, b):
        return x**(-b)

    def _stats(self, b, moments='mv'):
        mu, mu2, g1, g2 = None, None, None, None
        if 'm' in moments:
            mask = b > 1
            bt = np.extract(mask, b)
            mu = np.full(np.shape(b), fill_value=np.inf)
            np.place(mu, mask, bt / (bt-1.0))
        if 'v' in moments:
            mask = b > 2
            bt = np.extract(mask, b)
            mu2 = np.full(np.shape(b), fill_value=np.inf)
            np.place(mu2, mask, bt / (bt-2.0) / (bt-1.0)**2)
        if 's' in moments:
            mask = b > 3
            bt = np.extract(mask, b)
            g1 = np.full(np.shape(b), fill_value=np.nan)
            vals = 2 * (bt + 1.0) * np.sqrt(bt - 2.0) / ((bt - 3.0) * np.sqrt(bt))
            np.place(g1, mask, vals)
        if 'k' in moments:
            mask = b > 4
            bt = np.extract(mask, b)
            g2 = np.full(np.shape(b), fill_value=np.nan)
            vals = (6.0*np.polyval([1.0, 1.0, -6, -2], bt) /
                    np.polyval([1.0, -7.0, 12.0, 0.0], bt))
            np.place(g2, mask, vals)
        return mu, mu2, g1, g2

    def _entropy(self, c):
        return 1 + 1.0/c - np.log(c)

    @_call_super_mom
    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        parameters = _check_fit_input_parameters(self, data, args, kwds)
        data, fshape, floc, fscale = parameters

        # ensure that any fixed parameters don't violate constraints of the
        # distribution before continuing.
        if floc is not None and np.min(data) - floc < (fscale or 0):
            raise FitDataError("pareto", lower=1, upper=np.inf)

        ndata = data.shape[0]

        def get_shape(scale, location):
            # The first-order necessary condition on `shape` can be solved in
            # closed form
            return ndata / np.sum(np.log((data - location) / scale))

        if floc is fscale is None:
            # The support of the distribution is `(x - loc)/scale > 0`.
            # The method of Lagrange multipliers turns this constraint
            # into an equation that can be solved numerically.
            # See gh-12545 for details.

            def dL_dScale(shape, scale):
                # The partial derivative of the log-likelihood function w.r.t.
                # the scale.
                return ndata * shape / scale

            def dL_dLocation(shape, location):
                # The partial derivative of the log-likelihood function w.r.t.
                # the location.
                return (shape + 1) * np.sum(1 / (data - location))

            def fun_to_solve(scale):
                # optimize the scale by setting the partial derivatives
                # w.r.t. to location and scale equal and solving.
                location = np.min(data) - scale
                shape = fshape or get_shape(scale, location)
                return dL_dLocation(shape, location) - dL_dScale(shape, scale)

            def interval_contains_root(lbrack, rbrack):
                # return true if the signs disagree.
                return (np.sign(fun_to_solve(lbrack)) !=
                        np.sign(fun_to_solve(rbrack)))

            # set brackets for `root_scalar` to use when optimizing over the
            # scale such that a root is likely between them. Use user supplied
            # guess or default 1.
            brack_start = kwds.get('scale', 1)
            lbrack, rbrack = brack_start / 2, brack_start * 2
            # if a root is not between the brackets, iteratively expand them
            # until they include a sign change, checking after each bracket is
            # modified.
            while (not interval_contains_root(lbrack, rbrack)
                   and (lbrack > 0 or rbrack < np.inf)):
                lbrack /= 2
                rbrack *= 2
            res = root_scalar(fun_to_solve, bracket=[lbrack, rbrack])
            if res.converged:
                scale = res.root
                loc = np.min(data) - scale
                shape = fshape or get_shape(scale, loc)

                # The Pareto distribution requires that its parameters satisfy
                # the condition `fscale + floc <= min(data)`. However, to
                # avoid numerical issues, we require that `fscale + floc`
                # is strictly less than `min(data)`. If this condition
                # is not satisfied, reduce the scale with `np.nextafter` to
                # ensure that data does not fall outside of the support.
                if not (scale + loc) < np.min(data):
                    scale = np.min(data) - loc
                    scale = np.nextafter(scale, 0)
                return shape, loc, scale
            else:
                return super().fit(data, **kwds)
        elif floc is None:
            loc = np.min(data) - fscale
        else:
            loc = floc
        # Source: Evans, Hastings, and Peacock (2000), Statistical
        # Distributions, 3rd. Ed., John Wiley and Sons. Page 149.
        scale = fscale or np.min(data) - loc
        shape = fshape or get_shape(scale, loc)
        return shape, loc, scale


pareto = pareto_gen(a=1.0, name="pareto")


class lomax_gen(rv_continuous):
    r"""A Lomax (Pareto of the second kind) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `lomax` is:

    .. math::

        f(x, c) = \frac{c}{(1+x)^{c+1}}

    for :math:`x \ge 0`, :math:`c > 0`.

    `lomax` takes ``c`` as a shape parameter for :math:`c`.

    `lomax` is a special case of `pareto` with ``loc=-1.0``.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # lomax.pdf(x, c) = c / (1+x)**(c+1)
        return c*1.0/(1.0+x)**(c+1.0)

    def _logpdf(self, x, c):
        return np.log(c) - (c+1)*sc.log1p(x)

    def _cdf(self, x, c):
        return -sc.expm1(-c*sc.log1p(x))

    def _sf(self, x, c):
        return np.exp(-c*sc.log1p(x))

    def _logsf(self, x, c):
        return -c*sc.log1p(x)

    def _ppf(self, q, c):
        return sc.expm1(-sc.log1p(-q)/c)

    def _stats(self, c):
        mu, mu2, g1, g2 = pareto.stats(c, loc=-1.0, moments='mvsk')
        return mu, mu2, g1, g2

    def _entropy(self, c):
        return 1+1.0/c-np.log(c)


lomax = lomax_gen(a=0.0, name="lomax")


class pearson3_gen(rv_continuous):
    r"""A pearson type III continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `pearson3` is:

    .. math::

        f(x, \kappa) = \frac{|\beta|}{\Gamma(\alpha)}
                       (\beta (x - \zeta))^{\alpha - 1}
                       \exp(-\beta (x - \zeta))

    where:

    .. math::

            \beta = \frac{2}{\kappa}

            \alpha = \beta^2 = \frac{4}{\kappa^2}

            \zeta = -\frac{\alpha}{\beta} = -\beta

    :math:`\Gamma` is the gamma function (`scipy.special.gamma`).
    Pass the skew :math:`\kappa` into `pearson3` as the shape parameter
    ``skew``.

    %(after_notes)s

    %(example)s

    References
    ----------
    R.W. Vogel and D.E. McMartin, "Probability Plot Goodness-of-Fit and
    Skewness Estimation Procedures for the Pearson Type 3 Distribution", Water
    Resources Research, Vol.27, 3149-3158 (1991).

    L.R. Salvosa, "Tables of Pearson's Type III Function", Ann. Math. Statist.,
    Vol.1, 191-198 (1930).

    "Using Modern Computing Tools to Fit the Pearson Type III Distribution to
    Aviation Loads Data", Office of Aviation Research (2003).

    """
    def _preprocess(self, x, skew):
        # The real 'loc' and 'scale' are handled in the calling pdf(...). The
        # local variables 'loc' and 'scale' within pearson3._pdf are set to
        # the defaults just to keep them as part of the equations for
        # documentation.
        loc = 0.0
        scale = 1.0

        # If skew is small, return _norm_pdf. The divide between pearson3
        # and norm was found by brute force and is approximately a skew of
        # 0.000016.  No one, I hope, would actually use a skew value even
        # close to this small.
        norm2pearson_transition = 0.000016

        ans, x, skew = np.broadcast_arrays(1.0, x, skew)
        ans = ans.copy()

        # mask is True where skew is small enough to use the normal approx.
        mask = np.absolute(skew) < norm2pearson_transition
        invmask = ~mask

        beta = 2.0 / (skew[invmask] * scale)
        alpha = (scale * beta)**2
        zeta = loc - alpha / beta

        transx = beta * (x[invmask] - zeta)
        return ans, x, transx, mask, invmask, beta, alpha, zeta

    def _argcheck(self, skew):
        # The _argcheck function in rv_continuous only allows positive
        # arguments.  The skew argument for pearson3 can be zero (which I want
        # to handle inside pearson3._pdf) or negative.  So just return True
        # for all skew args.
        return np.isfinite(skew)

    def _shape_info(self):
        return [_ShapeInfo("skew", False, (-np.inf, np.inf), (False, False))]

    def _stats(self, skew):
        m = 0.0
        v = 1.0
        s = skew
        k = 1.5*skew**2
        return m, v, s, k

    def _pdf(self, x, skew):
        # pearson3.pdf(x, skew) = abs(beta) / gamma(alpha) *
        #     (beta * (x - zeta))**(alpha - 1) * exp(-beta*(x - zeta))
        # Do the calculation in _logpdf since helps to limit
        # overflow/underflow problems
        ans = np.exp(self._logpdf(x, skew))
        if ans.ndim == 0:
            if np.isnan(ans):
                return 0.0
            return ans
        ans[np.isnan(ans)] = 0.0
        return ans

    def _logpdf(self, x, skew):
        #   PEARSON3 logpdf                           GAMMA logpdf
        #   np.log(abs(beta))
        # + (alpha - 1)*np.log(beta*(x - zeta))          + (a - 1)*np.log(x)
        # - beta*(x - zeta)                           - x
        # - sc.gammalnalpha)                              - sc.gammalna)
        ans, x, transx, mask, invmask, beta, alpha, _ = (
            self._preprocess(x, skew))

        ans[mask] = np.log(_norm_pdf(x[mask]))
        # use logpdf instead of _logpdf to fix issue mentioned in gh-12640
        # (_logpdf does not return correct result for alpha = 1)
        ans[invmask] = np.log(abs(beta)) + gamma.logpdf(transx, alpha)
        return ans

    def _cdf(self, x, skew):
        ans, x, transx, mask, invmask, _, alpha, _ = (
            self._preprocess(x, skew))

        ans[mask] = _norm_cdf(x[mask])

        skew = np.broadcast_to(skew, invmask.shape)
        invmask1a = np.logical_and(invmask, skew > 0)
        invmask1b = skew[invmask] > 0
        # use cdf instead of _cdf to fix issue mentioned in gh-12640
        # (_cdf produces NaNs for inputs outside support)
        ans[invmask1a] = gamma.cdf(transx[invmask1b], alpha[invmask1b])

        # The gamma._cdf approach wasn't working with negative skew.
        # Note that multiplying the skew by -1 reflects about x=0.
        # So instead of evaluating the CDF with negative skew at x,
        # evaluate the SF with positive skew at -x.
        invmask2a = np.logical_and(invmask, skew < 0)
        invmask2b = skew[invmask] < 0
        # gamma._sf produces NaNs when transx < 0, so use gamma.sf
        ans[invmask2a] = gamma.sf(transx[invmask2b], alpha[invmask2b])

        return ans

    def _rvs(self, skew, size=None, random_state=None):
        skew = np.broadcast_to(skew, size)
        ans, _, _, mask, invmask, beta, alpha, zeta = (
            self._preprocess([0], skew))

        nsmall = mask.sum()
        nbig = mask.size - nsmall
        ans[mask] = random_state.standard_normal(nsmall)
        ans[invmask] = random_state.standard_gamma(alpha, nbig)/beta + zeta

        if size == ():
            ans = ans[0]
        return ans

    def _ppf(self, q, skew):
        ans, q, _, mask, invmask, beta, alpha, zeta = (
            self._preprocess(q, skew))
        ans[mask] = _norm_ppf(q[mask])
        ans[invmask] = sc.gammaincinv(alpha, q[invmask])/beta + zeta
        return ans

    @_call_super_mom
    @extend_notes_in_docstring(rv_continuous, notes="""\
        Note that method of moments (`method='MM'`) is not
        available for this distribution.\n\n""")
    def fit(self, data, *args, **kwds):
        if kwds.get("method", None) == 'MM':
            raise NotImplementedError("Fit `method='MM'` is not available for "
                                      "the Pearson3 distribution. Please try "
                                      "the default `method='MLE'`.")
        else:
            return super(type(self), self).fit(data, *args, **kwds)


pearson3 = pearson3_gen(name="pearson3")


class powerlaw_gen(rv_continuous):
    r"""A power-function continuous random variable.

    %(before_notes)s

    See Also
    --------
    pareto

    Notes
    -----
    The probability density function for `powerlaw` is:

    .. math::

        f(x, a) = a x^{a-1}

    for :math:`0 \le x \le 1`, :math:`a > 0`.

    `powerlaw` takes ``a`` as a shape parameter for :math:`a`.

    %(after_notes)s

    For example, the support of `powerlaw` can be adjusted from the default
    interval ``[0, 1]`` to the interval ``[c, c+d]`` by setting ``loc=c`` and
    ``scale=d``. For a power-law distribution with infinite support, see
    `pareto`.

    `powerlaw` is a special case of `beta` with ``b=1``.

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("a", False, (0, np.inf), (False, False))]

    def _pdf(self, x, a):
        # powerlaw.pdf(x, a) = a * x**(a-1)
        return a*x**(a-1.0)

    def _logpdf(self, x, a):
        return np.log(a) + sc.xlogy(a - 1, x)

    def _cdf(self, x, a):
        return x**(a*1.0)

    def _logcdf(self, x, a):
        return a*np.log(x)

    def _ppf(self, q, a):
        return pow(q, 1.0/a)

    def _stats(self, a):
        return (a / (a + 1.0),
                a / (a + 2.0) / (a + 1.0) ** 2,
                -2.0 * ((a - 1.0) / (a + 3.0)) * np.sqrt((a + 2.0) / a),
                6 * np.polyval([1, -1, -6, 2], a) / (a * (a + 3.0) * (a + 4)))

    def _entropy(self, a):
        return 1 - 1.0/a - np.log(a)


powerlaw = powerlaw_gen(a=0.0, b=1.0, name="powerlaw")


class powerlognorm_gen(rv_continuous):
    r"""A power log-normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powerlognorm` is:

    .. math::

        f(x, c, s) = \frac{c}{x s} \phi(\log(x)/s)
                     (\Phi(-\log(x)/s))^{c-1}

    where :math:`\phi` is the normal pdf, and :math:`\Phi` is the normal cdf,
    and :math:`x > 0`, :math:`s, c > 0`.

    `powerlognorm` takes :math:`c` and :math:`s` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        ic = _ShapeInfo("c", False, (0, np.inf), (False, False))
        i_s = _ShapeInfo("s", False, (0, np.inf), (False, False))
        return [ic, i_s]

    def _pdf(self, x, c, s):
        # powerlognorm.pdf(x, c, s) = c / (x*s) * phi(log(x)/s) *
        #                                         (Phi(-log(x)/s))**(c-1),
        return (c/(x*s) * _norm_pdf(np.log(x)/s) *
                pow(_norm_cdf(-np.log(x)/s), c*1.0-1.0))

    def _cdf(self, x, c, s):
        return 1.0 - pow(_norm_cdf(-np.log(x)/s), c*1.0)

    def _ppf(self, q, c, s):
        return np.exp(-s * _norm_ppf(pow(1.0 - q, 1.0 / c)))


powerlognorm = powerlognorm_gen(a=0.0, name="powerlognorm")


class powernorm_gen(rv_continuous):
    r"""A power normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powernorm` is:

    .. math::

        f(x, c) = c \phi(x) (\Phi(-x))^{c-1}

    where :math:`\phi` is the normal pdf, and :math:`\Phi` is the normal cdf,
    and :math:`x >= 0`, :math:`c > 0`.

    `powernorm` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    def _pdf(self, x, c):
        # powernorm.pdf(x, c) = c * phi(x) * (Phi(-x))**(c-1)
        return c*_norm_pdf(x) * (_norm_cdf(-x)**(c-1.0))

    def _logpdf(self, x, c):
        return np.log(c) + _norm_logpdf(x) + (c-1)*_norm_logcdf(-x)

    def _cdf(self, x, c):
        return 1.0-_norm_cdf(-x)**(c*1.0)

    def _ppf(self, q, c):
        return -_norm_ppf(pow(1.0 - q, 1.0 / c))


powernorm = powernorm_gen(name='powernorm')


class rdist_gen(rv_continuous):
    r"""An R-distributed (symmetric beta) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rdist` is:

    .. math::

        f(x, c) = \frac{(1-x^2)^{c/2-1}}{B(1/2, c/2)}

    for :math:`-1 \le x \le 1`, :math:`c > 0`. `rdist` is also called the
    symmetric beta distribution: if B has a `beta` distribution with
    parameters (c/2, c/2), then X = 2*B - 1 follows a R-distribution with
    parameter c.

    `rdist` takes ``c`` as a shape parameter for :math:`c`.

    This distribution includes the following distribution kernels as
    special cases::

        c = 2:  uniform
        c = 3:  `semicircular`
        c = 4:  Epanechnikov (parabolic)
        c = 6:  quartic (biweight)
        c = 8:  triweight

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, np.inf), (False, False))]

    # use relation to the beta distribution for pdf, cdf, etc
    def _pdf(self, x, c):
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return -np.log(2) + beta._logpdf((x + 1)/2, c/2, c/2)

    def _cdf(self, x, c):
        return beta._cdf((x + 1)/2, c/2, c/2)

    def _ppf(self, q, c):
        return 2*beta._ppf(q, c/2, c/2) - 1

    def _rvs(self, c, size=None, random_state=None):
        return 2 * random_state.beta(c/2, c/2, size) - 1

    def _munp(self, n, c):
        numerator = (1 - (n % 2)) * sc.beta((n + 1.0) / 2, c / 2.0)
        return numerator / sc.beta(1. / 2, c / 2.)


rdist = rdist_gen(a=-1.0, b=1.0, name="rdist")


def _rayleigh_fit_check_error(ier, msg):
    if ier != 1:
        raise RuntimeError('rayleigh.fit: fsolve failed to find the root of '
                           'the first-order conditions of the log-likelihood '
                           f'function: {msg} (ier={ier})')


class rayleigh_gen(rv_continuous):
    r"""A Rayleigh continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rayleigh` is:

    .. math::

        f(x) = x \exp(-x^2/2)

    for :math:`x \ge 0`.

    `rayleigh` is a special case of `chi` with ``df=2``.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return chi.rvs(2, size=size, random_state=random_state)

    def _pdf(self, r):
        # rayleigh.pdf(r) = r * exp(-r**2/2)
        return np.exp(self._logpdf(r))

    def _logpdf(self, r):
        return np.log(r) - 0.5 * r * r

    def _cdf(self, r):
        return -sc.expm1(-0.5 * r**2)

    def _ppf(self, q):
        return np.sqrt(-2 * sc.log1p(-q))

    def _sf(self, r):
        return np.exp(self._logsf(r))

    def _logsf(self, r):
        return -0.5 * r * r

    def _isf(self, q):
        return np.sqrt(-2 * np.log(q))

    def _stats(self):
        val = 4 - np.pi
        return (np.sqrt(np.pi/2),
                val/2,
                2*(np.pi-3)*np.sqrt(np.pi)/val**1.5,
                6*np.pi/val-16/val**2)

    def _entropy(self):
        return _EULER/2.0 + 1 - 0.5*np.log(2)

    @_call_super_mom
    @extend_notes_in_docstring(rv_continuous, notes="""\
        Notes specifically for ``rayleigh.fit``: If the location is fixed with
        the `floc` parameter, this method uses an analytical formula to find
        the scale.  Otherwise, this function uses a numerical root finder on
        the first order conditions of the log-likelihood function to find the
        MLE.  Only the (optional) `loc` parameter is used as the initial guess
        for the root finder; the `scale` parameter and any other parameters
        for the optimizer are ignored.\n\n""")
    def fit(self, data, *args, **kwds):
        data, floc, fscale = _check_fit_input_parameters(self, data,
                                                         args, kwds)

        def scale_mle(loc, data):
            # Source: Statistical Distributions, 3rd Edition. Evans, Hastings,
            # and Peacock (2000), Page 175
            return (np.sum((data - loc) ** 2) / (2 * len(data))) ** .5

        def loc_mle(loc, data):
            # This implicit equation for `loc` is used when
            # both `loc` and `scale` are free.
            xm = data - loc
            s1 = xm.sum()
            s2 = (xm**2).sum()
            s3 = (1/xm).sum()
            return s1 - s2/(2*len(data))*s3

        def loc_mle_scale_fixed(loc, scale, data):
            # This implicit equation for `loc` is used when
            # `scale` is fixed but `loc` is not.
            xm = data - loc
            return xm.sum() - scale**2 * (1/xm).sum()

        if floc is not None:
            # `loc` is fixed, analytically determine `scale`.
            if np.any(data - floc <= 0):
                raise FitDataError("rayleigh", lower=1, upper=np.inf)
            else:
                return floc, scale_mle(floc, data)

        # Account for user provided guess of `loc`.
        loc0 = kwds.get('loc')
        if loc0 is None:
            # Use _fitstart to estimate loc; ignore the returned scale.
            loc0 = self._fitstart(data)[0]

        if fscale is not None:
            # `scale` is fixed
            x, info, ier, msg = optimize.fsolve(loc_mle_scale_fixed, x0=loc0,
                                                args=(fscale, data,),
                                                xtol=1e-10, full_output=True)
            _rayleigh_fit_check_error(ier, msg)
            return x[0], fscale
        else:
            # Neither `loc` nor `scale` are fixed.
            x, info, ier, msg = optimize.fsolve(loc_mle, x0=loc0, args=(data,),
                                                xtol=1e-10, full_output=True)
            _rayleigh_fit_check_error(ier, msg)
            return x[0], scale_mle(x[0], data)


rayleigh = rayleigh_gen(a=0.0, name="rayleigh")


class reciprocal_gen(rv_continuous):
    r"""A loguniform or reciprocal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for this class is:

    .. math::

        f(x, a, b) = \frac{1}{x \log(b/a)}

    for :math:`a \le x \le b`, :math:`b > a > 0`. This class takes
    :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    This doesn't show the equal probability of ``0.01``, ``0.1`` and
    ``1``. This is best when the x-axis is log-scaled:

    >>> import numpy as np
    >>> fig, ax = plt.subplots(1, 1)
    >>> ax.hist(np.log10(r))
    >>> ax.set_ylabel("Frequency")
    >>> ax.set_xlabel("Value of random variable")
    >>> ax.xaxis.set_major_locator(plt.FixedLocator([-2, -1, 0]))
    >>> ticks = ["$10^{{ {} }}$".format(i) for i in [-2, -1, 0]]
    >>> ax.set_xticklabels(ticks)  # doctest: +SKIP
    >>> plt.show()

    This random variable will be log-uniform regardless of the base chosen for
    ``a`` and ``b``. Let's specify with base ``2`` instead:

    >>> rvs = %(name)s(2**-2, 2**0).rvs(size=1000)

    Values of ``1/4``, ``1/2`` and ``1`` are equally likely with this random
    variable.  Here's the histogram:

    >>> fig, ax = plt.subplots(1, 1)
    >>> ax.hist(np.log2(rvs))
    >>> ax.set_ylabel("Frequency")
    >>> ax.set_xlabel("Value of random variable")
    >>> ax.xaxis.set_major_locator(plt.FixedLocator([-2, -1, 0]))
    >>> ticks = ["$2^{{ {} }}$".format(i) for i in [-2, -1, 0]]
    >>> ax.set_xticklabels(ticks)  # doctest: +SKIP
    >>> plt.show()

    """
    def _argcheck(self, a, b):
        return (a > 0) & (b > a)

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (0, np.inf), (False, False))
        ib = _ShapeInfo("b", False, (0, np.inf), (False, False))
        return [ia, ib]

    def _fitstart(self, data):
        # Reasonable, since support is [a, b]
        return super()._fitstart(data, args=(np.min(data), np.max(data)))

    def _get_support(self, a, b):
        return a, b

    def _pdf(self, x, a, b):
        # reciprocal.pdf(x, a, b) = 1 / (x*log(b/a))
        return 1.0 / (x * np.log(b * 1.0 / a))

    def _logpdf(self, x, a, b):
        return -np.log(x) - np.log(np.log(b * 1.0 / a))

    def _cdf(self, x, a, b):
        return (np.log(x)-np.log(a)) / np.log(b * 1.0 / a)

    def _ppf(self, q, a, b):
        return a*pow(b*1.0/a, q)

    def _munp(self, n, a, b):
        return 1.0/np.log(b*1.0/a) / n * (pow(b*1.0, n) - pow(a*1.0, n))

    def _entropy(self, a, b):
        return 0.5*np.log(a*b)+np.log(np.log(b*1.0/a))

    fit_note = """\
        `loguniform`/`reciprocal` is over-parameterized. `fit` automatically
         fixes `scale` to 1 unless `fscale` is provided by the user.\n\n"""

    @extend_notes_in_docstring(rv_continuous, notes=fit_note)
    def fit(self, data, *args, **kwds):
        fscale = kwds.pop('fscale', 1)
        return super().fit(data, *args, fscale=fscale, **kwds)


loguniform = reciprocal_gen(name="loguniform")
reciprocal = reciprocal_gen(name="reciprocal")


class rice_gen(rv_continuous):
    r"""A Rice continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rice` is:

    .. math::

        f(x, b) = x \exp(- \frac{x^2 + b^2}{2}) I_0(x b)

    for :math:`x >= 0`, :math:`b > 0`. :math:`I_0` is the modified Bessel
    function of order zero (`scipy.special.i0`).

    `rice` takes ``b`` as a shape parameter for :math:`b`.

    %(after_notes)s

    The Rice distribution describes the length, :math:`r`, of a 2-D vector with
    components :math:`(U+u, V+v)`, where :math:`U, V` are constant, :math:`u,
    v` are independent Gaussian random variables with standard deviation
    :math:`s`.  Let :math:`R = \sqrt{U^2 + V^2}`. Then the pdf of :math:`r` is
    ``rice.pdf(x, R/s, scale=s)``.

    %(example)s

    """
    def _argcheck(self, b):
        return b >= 0

    def _shape_info(self):
        return [_ShapeInfo("b", False, (0, np.inf), (True, False))]

    def _rvs(self, b, size=None, random_state=None):
        # https://en.wikipedia.org/wiki/Rice_distribution
        t = b/np.sqrt(2) + random_state.standard_normal(size=(2,) + size)
        return np.sqrt((t*t).sum(axis=0))

    def _cdf(self, x, b):
        return sc.chndtr(np.square(x), 2, np.square(b))

    def _ppf(self, q, b):
        return np.sqrt(sc.chndtrix(q, 2, np.square(b)))

    def _pdf(self, x, b):
        # rice.pdf(x, b) = x * exp(-(x**2+b**2)/2) * I[0](x*b)
        #
        # We use (x**2 + b**2)/2 = ((x-b)**2)/2 + xb.
        # The factor of np.exp(-xb) is then included in the i0e function
        # in place of the modified Bessel function, i0, improving
        # numerical stability for large values of xb.
        return x * np.exp(-(x-b)*(x-b)/2.0) * sc.i0e(x*b)

    def _munp(self, n, b):
        nd2 = n/2.0
        n1 = 1 + nd2
        b2 = b*b/2.0
        return (2.0**(nd2) * np.exp(-b2) * sc.gamma(n1) *
                sc.hyp1f1(n1, 1, b2))


rice = rice_gen(a=0.0, name="rice")


class recipinvgauss_gen(rv_continuous):
    r"""A reciprocal inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `recipinvgauss` is:

    .. math::

        f(x, \mu) = \frac{1}{\sqrt{2\pi x}}
                    \exp\left(\frac{-(1-\mu x)^2}{2\mu^2x}\right)

    for :math:`x \ge 0`.

    `recipinvgauss` takes ``mu`` as a shape parameter for :math:`\mu`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("mu", False, (0, np.inf), (False, False))]

    def _pdf(self, x, mu):
        # recipinvgauss.pdf(x, mu) =
        #                     1/sqrt(2*pi*x) * exp(-(1-mu*x)**2/(2*x*mu**2))
        return np.exp(self._logpdf(x, mu))

    def _logpdf(self, x, mu):
        return _lazywhere(x > 0, (x, mu),
                          lambda x, mu: (-(1 - mu*x)**2.0 / (2*x*mu**2.0)
                                         - 0.5*np.log(2*np.pi*x)),
                          fillvalue=-np.inf)

    def _cdf(self, x, mu):
        trm1 = 1.0/mu - x
        trm2 = 1.0/mu + x
        isqx = 1.0/np.sqrt(x)
        return _norm_cdf(-isqx*trm1) - np.exp(2.0/mu)*_norm_cdf(-isqx*trm2)

    def _sf(self, x, mu):
        trm1 = 1.0/mu - x
        trm2 = 1.0/mu + x
        isqx = 1.0/np.sqrt(x)
        return _norm_cdf(isqx*trm1) + np.exp(2.0/mu)*_norm_cdf(-isqx*trm2)

    def _rvs(self, mu, size=None, random_state=None):
        return 1.0/random_state.wald(mu, 1.0, size=size)


recipinvgauss = recipinvgauss_gen(a=0.0, name='recipinvgauss')


class semicircular_gen(rv_continuous):
    r"""A semicircular continuous random variable.

    %(before_notes)s

    See Also
    --------
    rdist

    Notes
    -----
    The probability density function for `semicircular` is:

    .. math::

        f(x) = \frac{2}{\pi} \sqrt{1-x^2}

    for :math:`-1 \le x \le 1`.

    The distribution is a special case of `rdist` with `c = 3`.

    %(after_notes)s

    References
    ----------
    .. [1] "Wigner semicircle distribution",
           https://en.wikipedia.org/wiki/Wigner_semicircle_distribution

    %(example)s

    """
    def _shape_info(self):
        return []

    def _pdf(self, x):
        return 2.0/np.pi*np.sqrt(1-x*x)

    def _logpdf(self, x):
        return np.log(2/np.pi) + 0.5*sc.log1p(-x*x)

    def _cdf(self, x):
        return 0.5+1.0/np.pi*(x*np.sqrt(1-x*x) + np.arcsin(x))

    def _ppf(self, q):
        return rdist._ppf(q, 3)

    def _rvs(self, size=None, random_state=None):
        # generate values uniformly distributed on the area under the pdf
        # (semi-circle) by randomly generating the radius and angle
        r = np.sqrt(random_state.uniform(size=size))
        a = np.cos(np.pi * random_state.uniform(size=size))
        return r * a

    def _stats(self):
        return 0, 0.25, 0, -1.0

    def _entropy(self):
        return 0.64472988584940017414


semicircular = semicircular_gen(a=-1.0, b=1.0, name="semicircular")


class skewcauchy_gen(rv_continuous):
    r"""A skewed Cauchy random variable.

    %(before_notes)s

    See Also
    --------
    cauchy : Cauchy distribution

    Notes
    -----

    The probability density function for `skewcauchy` is:

    .. math::

        f(x) = \frac{1}{\pi \left(\frac{x^2}{\left(a\, \text{sign}(x) + 1
                                                   \right)^2} + 1 \right)}

    for a real number :math:`x` and skewness parameter :math:`-1 < a < 1`.

    When :math:`a=0`, the distribution reduces to the usual Cauchy
    distribution.

    %(after_notes)s

    References
    ----------
    .. [1] "Skewed generalized *t* distribution", Wikipedia
       https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution#Skewed_Cauchy_distribution

    %(example)s

    """
    def _argcheck(self, a):
        return np.abs(a) < 1

    def _shape_info(self):
        return [_ShapeInfo("a", False, (-1.0, 1.0), (False, False))]

    def _pdf(self, x, a):
        return 1 / (np.pi * (x**2 / (a * np.sign(x) + 1)**2 + 1))

    def _cdf(self, x, a):
        return np.where(x <= 0,
                        (1 - a) / 2 + (1 - a) / np.pi * np.arctan(x / (1 - a)),
                        (1 - a) / 2 + (1 + a) / np.pi * np.arctan(x / (1 + a)))

    def _ppf(self, x, a):
        i = x < self._cdf(0, a)
        return np.where(i,
                        np.tan(np.pi / (1 - a) * (x - (1 - a) / 2)) * (1 - a),
                        np.tan(np.pi / (1 + a) * (x - (1 - a) / 2)) * (1 + a))

    def _stats(self, a, moments='mvsk'):
        return np.nan, np.nan, np.nan, np.nan

    def _fitstart(self, data):
        # Use 0 as the initial guess of the skewness shape parameter.
        # For the location and scale, estimate using the median and
        # quartiles.
        p25, p50, p75 = np.percentile(data, [25, 50, 75])
        return 0.0, p50, (p75 - p25)/2


skewcauchy = skewcauchy_gen(name='skewcauchy')


class skew_norm_gen(rv_continuous):
    r"""A skew-normal random variable.

    %(before_notes)s

    Notes
    -----
    The pdf is::

        skewnorm.pdf(x, a) = 2 * norm.pdf(x) * norm.cdf(a*x)

    `skewnorm` takes a real number :math:`a` as a skewness parameter
    When ``a = 0`` the distribution is identical to a normal distribution
    (`norm`). `rvs` implements the method of [1]_.

    %(after_notes)s

    %(example)s

    References
    ----------
    .. [1] A. Azzalini and A. Capitanio (1999). Statistical applications of
        the multivariate skew-normal distribution. J. Roy. Statist. Soc.,
        B 61, 579-602. :arxiv:`0911.2093`

    """
    def _argcheck(self, a):
        return np.isfinite(a)

    def _shape_info(self):
        return [_ShapeInfo("a", False, (-np.inf, np.inf), (False, False))]

    def _pdf(self, x, a):
        return _lazywhere(
            a == 0, (x, a), lambda x, a: _norm_pdf(x),
            f2=lambda x, a: 2.*_norm_pdf(x)*_norm_cdf(a*x)
        )

    def _cdf(self, x, a):
        cdf = _boost._skewnorm_cdf(x, 0, 1, a)
        # Boost is not accurate in left tail when a > 0
        i_small_cdf = (cdf < 1e-6) & (a > 0)
        cdf[i_small_cdf] = super()._cdf(x[i_small_cdf], a[i_small_cdf])
        return np.clip(cdf, 0, 1)

    def _ppf(self, x, a):
        return _boost._skewnorm_ppf(x, 0, 1, a)

    def _sf(self, x, a):
        # Boost's SF is implemented this way. Use whatever customizations
        # we made in the _cdf.
        return self._cdf(-x, -a)

    def _isf(self, x, a):
        return _boost._skewnorm_isf(x, 0, 1, a)

    def _rvs(self, a, size=None, random_state=None):
        u0 = random_state.normal(size=size)
        v = random_state.normal(size=size)
        d = a/np.sqrt(1 + a**2)
        u1 = d*u0 + v*np.sqrt(1 - d**2)
        return np.where(u0 >= 0, u1, -u1)

    def _stats(self, a, moments='mvsk'):
        output = [None, None, None, None]
        const = np.sqrt(2/np.pi) * a/np.sqrt(1 + a**2)

        if 'm' in moments:
            output[0] = const
        if 'v' in moments:
            output[1] = 1 - const**2
        if 's' in moments:
            output[2] = ((4 - np.pi)/2) * (const/np.sqrt(1 - const**2))**3
        if 'k' in moments:
            output[3] = (2*(np.pi - 3)) * (const**4/(1 - const**2)**2)

        return output

    # For odd order, the each noncentral moment of the skew-normal distribution
    # with location 0 and scale 1 can be expressed as a polynomial in delta,
    # where delta = a/sqrt(1 + a**2) and `a` is the skew-normal shape
    # parameter.  The dictionary _skewnorm_odd_moments defines those
    # polynomials for orders up to 19.  The dict is implemented as a cached
    # property to reduce the impact of the creation of the dict on import time.
    @cached_property
    def _skewnorm_odd_moments(self):
        skewnorm_odd_moments = {
            1: Polynomial([1]),
            3: Polynomial([3, -1]),
            5: Polynomial([15, -10, 3]),
            7: Polynomial([105, -105, 63, -15]),
            9: Polynomial([945, -1260, 1134, -540, 105]),
            11: Polynomial([10395, -17325, 20790, -14850, 5775, -945]),
            13: Polynomial([135135, -270270, 405405, -386100, 225225, -73710,
                            10395]),
            15: Polynomial([2027025, -4729725, 8513505, -10135125, 7882875,
                            -3869775, 1091475, -135135]),
            17: Polynomial([34459425, -91891800, 192972780, -275675400,
                            268017750, -175429800, 74220300, -18378360,
                            2027025]),
            19: Polynomial([654729075, -1964187225, 4714049340, -7856748900,
                            9166207050, -7499623950, 4230557100, -1571349780,
                            346621275, -34459425]),
        }
        return skewnorm_odd_moments

    def _munp(self, order, a):
        if order & 1:
            if order > 19:
                raise NotImplementedError("skewnorm noncentral moments not "
                                          "implemented for odd orders greater "
                                          "than 19.")
            # Use the precomputed polynomials that were derived from the
            # moment generating function.
            delta = a/np.sqrt(1 + a**2)
            return (delta * self._skewnorm_odd_moments[order](delta**2)
                    * _SQRT_2_OVER_PI)
        else:
            # For even order, the moment is just (order-1)!!, where !! is the
            # notation for the double factorial; for an odd integer m, m!! is
            # m*(m-2)*...*3*1.
            # We could use special.factorial2, but we know the argument is odd,
            # so avoid the overhead of that function and compute the result
            # directly here.
            return sc.gamma((order + 1)/2) * 2**(order/2) / _SQRT_PI

    @extend_notes_in_docstring(rv_continuous, notes="""\
        If ``method='mm'``, parameters fixed by the user are respected, and the
        remaining parameters are used to match distribution and sample moments
        where possible. For example, if the user fixes the location with
        ``floc``, the parameters will only match the distribution skewness and
        variance to the sample skewness and variance; no attempt will be made
        to match the means or minimize a norm of the errors.
        Note that the maximum possible skewness magnitude of a
        `scipy.stats.skewnorm` distribution is approximately 0.9952717; if the
        magnitude of the data's sample skewness exceeds this, the returned
        shape parameter ``a`` will be infinite.
        \n\n""")
    def fit(self, data, *args, **kwds):
        # this extracts fixed shape, location, and scale however they
        # are specified, and also leaves them in `kwds`
        data, fa, floc, fscale = _check_fit_input_parameters(self, data,
                                                             args, kwds)
        method = kwds.get("method", "mle").lower()

        # See https://en.wikipedia.org/wiki/Skew_normal_distribution for
        # moment formulas.
        def skew_d(d):  # skewness in terms of delta
            return (4-np.pi)/2 * ((d * np.sqrt(2 / np.pi))**3
                                  / (1 - 2*d**2 / np.pi)**(3/2))

        # If skewness of data is greater than max possible population skewness,
        # MoM won't provide a good guess. Get out early.
        s = stats.skew(data)
        s_max = skew_d(1)
        if abs(s) >= s_max and method == "mle" and fa is None and not args:
            return super().fit(data, *args, **kwds)

        # If method is method of moments, we don't need the user's guesses.
        # If method is "mle", extract the guesses from args and kwds.
        if method == "mm":
            a, loc, scale = None, None, None
        else:
            a = args[0] if len(args) else None
            loc = kwds.pop('loc', None)
            scale = kwds.pop('scale', None)

        if fa is None and a is None:  # not fixed and no guess: use MoM
            # Solve for a that matches sample distribution skewness to sample
            # skewness.
            s = np.clip(s, -s_max, s_max)
            d = root_scalar(lambda d: skew_d(d) - s, bracket=[-1, 1]).root
            with np.errstate(divide='ignore'):
                a = np.sqrt(np.divide(d**2, (1-d**2)))*np.sign(s)
        else:
            a = fa if fa is not None else a
            d = a / np.sqrt(1 + a**2)

        if fscale is None and scale is None:
            v = np.var(data)
            scale = np.sqrt(v / (1 - 2*d**2/np.pi))
        elif fscale is not None:
            scale = fscale

        if floc is None and loc is None:
            m = np.mean(data)
            loc = m - scale*d*np.sqrt(2/np.pi)
        elif floc is not None:
            loc = floc

        if method == 'mm':
            return a, loc, scale
        else:
            # At this point, parameter "guesses" may equal the fixed parameters
            # in kwds. No harm in passing them as guesses, too.
            return super().fit(data, a, loc=loc, scale=scale, **kwds)


skewnorm = skew_norm_gen(name='skewnorm')


class trapezoid_gen(rv_continuous):
    r"""A trapezoidal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The trapezoidal distribution can be represented with an up-sloping line
    from ``loc`` to ``(loc + c*scale)``, then constant to ``(loc + d*scale)``
    and then downsloping from ``(loc + d*scale)`` to ``(loc+scale)``.  This
    defines the trapezoid base from ``loc`` to ``(loc+scale)`` and the flat
    top from ``c`` to ``d`` proportional to the position along the base
    with ``0 <= c <= d <= 1``.  When ``c=d``, this is equivalent to `triang`
    with the same values for `loc`, `scale` and `c`.
    The method of [1]_ is used for computing moments.

    `trapezoid` takes :math:`c` and :math:`d` as shape parameters.

    %(after_notes)s

    The standard form is in the range [0, 1] with c the mode.
    The location parameter shifts the start to `loc`.
    The scale parameter changes the width from 1 to `scale`.

    %(example)s

    References
    ----------
    .. [1] Kacker, R.N. and Lawrence, J.F. (2007). Trapezoidal and triangular
       distributions for Type B evaluation of standard uncertainty.
       Metrologia 44, 117-127. :doi:`10.1088/0026-1394/44/2/003`


    """
    def _argcheck(self, c, d):
        return (c >= 0) & (c <= 1) & (d >= 0) & (d <= 1) & (d >= c)

    def _shape_info(self):
        ic = _ShapeInfo("c", False, (0, 1.0), (True, True))
        id = _ShapeInfo("d", False, (0, 1.0), (True, True))
        return [ic, id]

    def _pdf(self, x, c, d):
        u = 2 / (d-c+1)

        return _lazyselect([x < c,
                            (c <= x) & (x <= d),
                            x > d],
                           [lambda x, c, d, u: u * x / c,
                            lambda x, c, d, u: u,
                            lambda x, c, d, u: u * (1-x) / (1-d)],
                            (x, c, d, u))

    def _cdf(self, x, c, d):
        return _lazyselect([x < c,
                            (c <= x) & (x <= d),
                            x > d],
                           [lambda x, c, d: x**2 / c / (d-c+1),
                            lambda x, c, d: (c + 2 * (x-c)) / (d-c+1),
                            lambda x, c, d: 1-((1-x) ** 2
                                               / (d-c+1) / (1-d))],
                            (x, c, d))

    def _ppf(self, q, c, d):
        qc, qd = self._cdf(c, c, d), self._cdf(d, c, d)
        condlist = [q < qc, q <= qd, q > qd]
        choicelist = [np.sqrt(q * c * (1 + d - c)),
                      0.5 * q * (1 + d - c) + 0.5 * c,
                      1 - np.sqrt((1 - q) * (d - c + 1) * (1 - d))]
        return np.select(condlist, choicelist)

    def _munp(self, n, c, d):
        # Using the parameterization from Kacker, 2007, with
        # a=bottom left, c=top left, d=top right, b=bottom right, then
        #     E[X^n] = h/(n+1)/(n+2) [(b^{n+2}-d^{n+2})/(b-d)
        #                             - ((c^{n+2} - a^{n+2})/(c-a)]
        # with h = 2/((b-a) - (d-c)). The corresponding parameterization
        # in scipy, has a'=loc, c'=loc+c*scale, d'=loc+d*scale, b'=loc+scale,
        # which for standard form reduces to a'=0, b'=1, c'=c, d'=d.
        # Substituting into E[X^n] gives the bd' term as (1 - d^{n+2})/(1 - d)
        # and the ac' term as c^{n-1} for the standard form. The bd' term has
        # numerical difficulties near d=1, so replace (1 - d^{n+2})/(1-d)
        # with expm1((n+2)*log(d))/(d-1).
        # Testing with n=18 for c=(1e-30,1-eps) shows that this is stable.
        # We still require an explicit test for d=1 to prevent divide by zero,
        # and now a test for d=0 to prevent log(0).
        ab_term = c**(n+1)
        dc_term = _lazyselect(
            [d == 0.0, (0.0 < d) & (d < 1.0), d == 1.0],
            [lambda d: 1.0,
             lambda d: np.expm1((n+2) * np.log(d)) / (d-1.0),
             lambda d: n+2],
            [d])
        val = 2.0 / (1.0+d-c) * (dc_term - ab_term) / ((n+1) * (n+2))
        return val

    def _entropy(self, c, d):
        # Using the parameterization from Wikipedia (van Dorp, 2003)
        # with a=bottom left, c=top left, d=top right, b=bottom right
        # gives a'=loc, b'=loc+c*scale, c'=loc+d*scale, d'=loc+scale,
        # which for loc=0, scale=1 is a'=0, b'=c, c'=d, d'=1.
        # Substituting into the entropy formula from Wikipedia gives
        # the following result.
        return 0.5 * (1.0-d+c) / (1.0+d-c) + np.log(0.5 * (1.0+d-c))


trapezoid = trapezoid_gen(a=0.0, b=1.0, name="trapezoid")
# Note: alias kept for backwards compatibility. Rename was done
# because trapz is a slur in colloquial English (see gh-12924).
trapz = trapezoid_gen(a=0.0, b=1.0, name="trapz")
if trapz.__doc__:
    trapz.__doc__ = "trapz is an alias for `trapezoid`"


class triang_gen(rv_continuous):
    r"""A triangular continuous random variable.

    %(before_notes)s

    Notes
    -----
    The triangular distribution can be represented with an up-sloping line from
    ``loc`` to ``(loc + c*scale)`` and then downsloping for ``(loc + c*scale)``
    to ``(loc + scale)``.

    `triang` takes ``c`` as a shape parameter for :math:`0 \le c \le 1`.

    %(after_notes)s

    The standard form is in the range [0, 1] with c the mode.
    The location parameter shifts the start to `loc`.
    The scale parameter changes the width from 1 to `scale`.

    %(example)s

    """
    def _rvs(self, c, size=None, random_state=None):
        return random_state.triangular(0, c, 1, size)

    def _argcheck(self, c):
        return (c >= 0) & (c <= 1)

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, 1.0), (True, True))]

    def _pdf(self, x, c):
        # 0: edge case where c=0
        # 1: generalised case for x < c, don't use x <= c, as it doesn't cope
        #    with c = 0.
        # 2: generalised case for x >= c, but doesn't cope with c = 1
        # 3: edge case where c=1
        r = _lazyselect([c == 0,
                         x < c,
                         (x >= c) & (c != 1),
                         c == 1],
                        [lambda x, c: 2 - 2 * x,
                         lambda x, c: 2 * x / c,
                         lambda x, c: 2 * (1 - x) / (1 - c),
                         lambda x, c: 2 * x],
                        (x, c))
        return r

    def _cdf(self, x, c):
        r = _lazyselect([c == 0,
                         x < c,
                         (x >= c) & (c != 1),
                         c == 1],
                        [lambda x, c: 2*x - x*x,
                         lambda x, c: x * x / c,
                         lambda x, c: (x*x - 2*x + c) / (c-1),
                         lambda x, c: x * x],
                        (x, c))
        return r

    def _ppf(self, q, c):
        return np.where(q < c, np.sqrt(c * q), 1-np.sqrt((1-c) * (1-q)))

    def _stats(self, c):
        return ((c+1.0)/3.0,
                (1.0-c+c*c)/18,
                np.sqrt(2)*(2*c-1)*(c+1)*(c-2) / (5*np.power((1.0-c+c*c), 1.5)),
                -3.0/5.0)

    def _entropy(self, c):
        return 0.5-np.log(2)


triang = triang_gen(a=0.0, b=1.0, name="triang")


class truncexpon_gen(rv_continuous):
    r"""A truncated exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `truncexpon` is:

    .. math::

        f(x, b) = \frac{\exp(-x)}{1 - \exp(-b)}

    for :math:`0 <= x <= b`.

    `truncexpon` takes ``b`` as a shape parameter for :math:`b`.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("b", False, (0, np.inf), (False, False))]

    def _get_support(self, b):
        return self.a, b

    def _pdf(self, x, b):
        # truncexpon.pdf(x, b) = exp(-x) / (1-exp(-b))
        return np.exp(-x)/(-sc.expm1(-b))

    def _logpdf(self, x, b):
        return -x - np.log(-sc.expm1(-b))

    def _cdf(self, x, b):
        return sc.expm1(-x)/sc.expm1(-b)

    def _ppf(self, q, b):
        return -sc.log1p(q*sc.expm1(-b))

    def _munp(self, n, b):
        # wrong answer with formula, same as in continuous.pdf
        # return sc.gamman+1)-sc.gammainc1+n, b)
        if n == 1:
            return (1-(b+1)*np.exp(-b))/(-sc.expm1(-b))
        elif n == 2:
            return 2*(1-0.5*(b*b+2*b+2)*np.exp(-b))/(-sc.expm1(-b))
        else:
            # return generic for higher moments
            # return rv_continuous._mom1_sc(self, n, b)
            return self._mom1_sc(n, b)

    def _entropy(self, b):
        eB = np.exp(b)
        return np.log(eB-1)+(1+eB*(b-1.0))/(1.0-eB)


truncexpon = truncexpon_gen(a=0.0, name='truncexpon')


# logsumexp trick for log(p + q) with only log(p) and log(q)
def _log_sum(log_p, log_q):
    return sc.logsumexp([log_p, log_q], axis=0)


# same as above, but using -exp(x) = exp(x + πi)
def _log_diff(log_p, log_q):
    # need to broadcast in case a or b is 0.5; logsumexp doesn't
    log_p, log_q = np.broadcast_arrays(log_p, log_q)
    return sc.logsumexp([log_p, log_q+np.pi*1j], axis=0)


def _log_gauss_mass(a, b):
    """Log of Gaussian probability mass within an interval"""
    a, b = np.atleast_1d(a), np.atleast_1d(b)
    a, b = np.broadcast_arrays(a, b)

    # Calculations in right tail are inaccurate, so we'll exploit the
    # symmetry and work only in the left tail
    case_left = b <= 0.5
    case_right = a > 0.5
    case_central = ~(case_left | case_right)

    def mass_case_left(a, b):
        return _log_diff(sc.log_ndtr(b), sc.log_ndtr(a))

    def mass_case_right(a, b):
        return mass_case_left(-b, -a)

    def mass_case_central(a, b):
        left_mass = mass_case_left(a, 0.5)
        right_mass = mass_case_right(0.5, b)
        return _log_sum(left_mass, right_mass)

    # _lazyselect not working; don't care to debug it
    out = np.full_like(a, fill_value=np.nan, dtype=np.complex128)
    if a[case_left].size:
        out[case_left] = mass_case_left(a[case_left], b[case_left])
    if a[case_right].size:
        out[case_right] = mass_case_right(a[case_right], b[case_right])
    if a[case_central].size:
        out[case_central] = mass_case_central(a[case_central], b[case_central])
    return np.real(out)  # discard ~0j


class truncnorm_gen(rv_continuous):
    r"""A truncated normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    This distribution is the normal distribution centered on ``loc`` (default
    0), with standard deviation ``scale`` (default 1), and clipped at ``a``,
    ``b`` standard deviations to the left, right (respectively) from ``loc``.
    If ``myclip_a`` and ``myclip_b`` are clip values in the sample space (as
    opposed to the number of standard deviations) then they can be converted
    to the required form according to::

        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale

    %(example)s

    """

    def _argcheck(self, a, b):
        return a < b

    def _shape_info(self):
        ia = _ShapeInfo("a", False, (-np.inf, np.inf), (True, False))
        ib = _ShapeInfo("b", False, (-np.inf, np.inf), (False, True))
        return [ia, ib]

    def _fitstart(self, data):
        # Reasonable, since support is [a, b]
        return super()._fitstart(data, args=(np.min(data), np.max(data)))

    def _get_support(self, a, b):
        return a, b

    def _pdf(self, x, a, b):
        return np.exp(self._logpdf(x, a, b))

    def _logpdf(self, x, a, b):
        return _norm_logpdf(x) - _log_gauss_mass(a, b)

    def _cdf(self, x, a, b):
        return np.exp(self._logcdf(x, a, b))

    def _logcdf(self, x, a, b):
        return _log_gauss_mass(a, x) - _log_gauss_mass(a, b)

    def _sf(self, x, a, b):
        return np.exp(self._logsf(x, a, b))

    def _logsf(self, x, a, b):
        return _log_gauss_mass(x, b) - _log_gauss_mass(a, b)

    def _ppf(self, q, a, b):
        q, a, b = np.broadcast_arrays(q, a, b)

        case_left = a < 0
        case_right = ~case_left

        def ppf_left(q, a, b):
            log_Phi_x = _log_sum(sc.log_ndtr(a),
                                 np.log(q) + _log_gauss_mass(a, b))
            return sc.ndtri_exp(log_Phi_x)

        def ppf_right(q, a, b):
            log_Phi_x = _log_sum(sc.log_ndtr(-b),
                                 np.log1p(-q) + _log_gauss_mass(a, b))
            return -sc.ndtri_exp(log_Phi_x)

        out = np.empty_like(q)

        q_left = q[case_left]
        q_right = q[case_right]

        if q_left.size:
            out[case_left] = ppf_left(q_left, a[case_left], b[case_left])
        if q_right.size:
            out[case_right] = ppf_right(q_right, a[case_right], b[case_right])

        return out

    def _isf(self, q, a, b):
        # Mostly copy-paste of _ppf, but I think this is simpler than combining
        q, a, b = np.broadcast_arrays(q, a, b)

        case_left = b < 0
        case_right = ~case_left

        def isf_left(q, a, b):
            log_Phi_x = _log_diff(sc.log_ndtr(b),
                                  np.log(q) + _log_gauss_mass(a, b))
            return sc.ndtri_exp(np.real(log_Phi_x))

        def isf_right(q, a, b):
            log_Phi_x = _log_diff(sc.log_ndtr(-a),
                                  np.log1p(-q) + _log_gauss_mass(a, b))
            return -sc.ndtri_exp(np.real(log_Phi_x))

        out = np.empty_like(q)

        q_left = q[case_left]
        q_right = q[case_right]

        if q_left.size:
            out[case_left] = isf_left(q_left, a[case_left], b[case_left])
        if q_right.size:
            out[case_right] = isf_right(q_right, a[case_right], b[case_right])

        return out

    def _munp(self, n, a, b):
        def n_th_moment(n, a, b):
            """
            Returns n-th moment. Defined only if n >= 0.
            Function cannot broadcast due to the loop over n
            """
            pA, pB = self._pdf([a, b], a, b)
            probs = [pA, -pB]
            moments = [0, 1]
            for k in range(1, n+1):
                # a or b might be infinite, and the corresponding pdf value
                # is 0 in that case, but nan is returned for the
                # multiplication.  However, as b->infinity,  pdf(b)*b**k -> 0.
                # So it is safe to use _lazywhere to avoid the nan.
                vals = _lazywhere(probs, [probs, [a, b]],
                                  lambda x, y: x * y**(k-1), fillvalue=0)
                mk = np.sum(vals) + (k-1) * moments[-2]
                moments.append(mk)
            return moments[-1]

        return _lazywhere((n >= 0) & (a == a) & (b == b), (n, a, b),
                          np.vectorize(n_th_moment, otypes=[np.float64]),
                          np.nan)

    def _stats(self, a, b, moments='mv'):
        pA, pB = self.pdf(np.array([a, b]), a, b)

        def _truncnorm_stats_scalar(a, b, pA, pB, moments):
            m1 = pA - pB
            mu = m1
            # use _lazywhere to avoid nan (See detailed comment in _munp)
            probs = [pA, -pB]
            vals = _lazywhere(probs, [probs, [a, b]], lambda x, y: x*y,
                              fillvalue=0)
            m2 = 1 + np.sum(vals)
            vals = _lazywhere(probs, [probs, [a-mu, b-mu]], lambda x, y: x*y,
                              fillvalue=0)
            # mu2 = m2 - mu**2, but not as numerically stable as:
            # mu2 = (a-mu)*pA - (b-mu)*pB + 1
            mu2 = 1 + np.sum(vals)
            vals = _lazywhere(probs, [probs, [a, b]], lambda x, y: x*y**2,
                              fillvalue=0)
            m3 = 2*m1 + np.sum(vals)
            vals = _lazywhere(probs, [probs, [a, b]], lambda x, y: x*y**3,
                              fillvalue=0)
            m4 = 3*m2 + np.sum(vals)

            mu3 = m3 + m1 * (-3*m2 + 2*m1**2)
            g1 = mu3 / np.power(mu2, 1.5)
            mu4 = m4 + m1*(-4*m3 + 3*m1*(2*m2 - m1**2))
            g2 = mu4 / mu2**2 - 3
            return mu, mu2, g1, g2

        _truncnorm_stats = np.vectorize(_truncnorm_stats_scalar,
                                        excluded=('moments',))
        return _truncnorm_stats(a, b, pA, pB, moments)


truncnorm = truncnorm_gen(name='truncnorm', momtype=1)


class truncpareto_gen(rv_continuous):
    r"""An upper truncated Pareto continuous random variable.

    %(before_notes)s

    See Also
    --------
    pareto : Pareto distribution

    Notes
    -----
    The probability density function for `truncpareto` is:

    .. math::

        f(x, b, c) = \frac{b}{1 - c^{-b}} \frac{1}{x^{b+1}}

    for :math:`b > 0`, :math:`c > 1` and :math:`1 \le x \le c`.

    `truncpareto` takes `b` and `c` as shape parameters for :math:`b` and
    :math:`c`.

    Notice that the upper truncation value :math:`c` is defined in
    standardized form so that random values of an unscaled, unshifted variable
    are within the range ``[1, c]``.
    If ``u_r`` is the upper bound to a scaled and/or shifted variable,
    then ``c = (u_r - loc) / scale``. In other words, the support of the
    distribution becomes ``(scale + loc) <= x <= (c*scale + loc)`` when
    `scale` and/or `loc` are provided.

    %(after_notes)s

    References
    ----------
    .. [1] Burroughs, S. M., and Tebbens S. F.
        "Upper-truncated power laws in natural systems."
        Pure and Applied Geophysics 158.4 (2001): 741-757.

    %(example)s

    """

    def _shape_info(self):
        ib = _ShapeInfo("b", False, (0.0, np.inf), (False, False))
        ic = _ShapeInfo("c", False, (1.0, np.inf), (False, False))
        return [ib, ic]

    def _argcheck(self, b, c):
        return (b > 0) & (c > 1)

    def _get_support(self, b, c):
        return self.a, c

    def _pdf(self, x, b, c):
        return b * x**-(b+1) / (1 - c**-b)

    def _logpdf(self, x, b, c):
        return np.log(b) - np.log1p(-c**-b) - (b+1)*np.log(x)

    def _cdf(self, x, b, c):
        return (1 - x**-b) / (1 - c**-b)

    def _logcdf(self, x, b, c):
        return np.log1p(-x**-b) - np.log1p(-c**-b)

    def _ppf(self, q, b, c):
        return pow(1 - (1 - c**-b)*q, -1/b)

    def _sf(self, x, b, c):
        return (x**-b - c**-b) / (1 - c**-b)

    def _logsf(self, x, b, c):
        return np.log(x**-b - c**-b) - np.log1p(-c**-b)

    def _isf(self, q, b, c):
        return pow(c**-b + (1 - c**-b)*q, -1/b)

    def _entropy(self, b, c):
        return -(np.log(b/(1 - c**-b))
                 + (b+1)*(np.log(c)/(c**b - 1) - 1/b))

    def _munp(self, n, b, c):
        if n == b:
            return b*np.log(c) / (1 - c**-b)
        else:
            return b / (b-n) * (c**b - c**n) / (c**b - 1)

    def _fitstart(self, data):
        b, loc, scale = pareto.fit(data)
        c = (max(data) - loc)/scale
        return b, c, loc, scale


truncpareto = truncpareto_gen(a=1.0, name='truncpareto')


class tukeylambda_gen(rv_continuous):
    r"""A Tukey-Lamdba continuous random variable.

    %(before_notes)s

    Notes
    -----
    A flexible distribution, able to represent and interpolate between the
    following distributions:

    - Cauchy                (:math:`lambda = -1`)
    - logistic              (:math:`lambda = 0`)
    - approx Normal         (:math:`lambda = 0.14`)
    - uniform from -1 to 1  (:math:`lambda = 1`)

    `tukeylambda` takes a real number :math:`lambda` (denoted ``lam``
    in the implementation) as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, lam):
        return np.isfinite(lam)

    def _shape_info(self):
        return [_ShapeInfo("lam", False, (-np.inf, np.inf), (False, False))]

    def _pdf(self, x, lam):
        Fx = np.asarray(sc.tklmbda(x, lam))
        Px = Fx**(lam-1.0) + (np.asarray(1-Fx))**(lam-1.0)
        Px = 1.0/np.asarray(Px)
        return np.where((lam <= 0) | (abs(x) < 1.0/np.asarray(lam)), Px, 0.0)

    def _cdf(self, x, lam):
        return sc.tklmbda(x, lam)

    def _ppf(self, q, lam):
        return sc.boxcox(q, lam) - sc.boxcox1p(-q, lam)

    def _stats(self, lam):
        return 0, _tlvar(lam), 0, _tlkurt(lam)

    def _entropy(self, lam):
        def integ(p):
            return np.log(pow(p, lam-1)+pow(1-p, lam-1))
        return integrate.quad(integ, 0, 1)[0]


tukeylambda = tukeylambda_gen(name='tukeylambda')


class FitUniformFixedScaleDataError(FitDataError):
    def __init__(self, ptp, fscale):
        self.args = (
            "Invalid values in `data`.  Maximum likelihood estimation with "
            "the uniform distribution and fixed scale requires that "
            "data.ptp() <= fscale, but data.ptp() = %r and fscale = %r." %
            (ptp, fscale),
        )


class uniform_gen(rv_continuous):
    r"""A uniform continuous random variable.

    In the standard form, the distribution is uniform on ``[0, 1]``. Using
    the parameters ``loc`` and ``scale``, one obtains the uniform distribution
    on ``[loc, loc + scale]``.

    %(before_notes)s

    %(example)s

    """
    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return random_state.uniform(0.0, 1.0, size)

    def _pdf(self, x):
        return 1.0*(x == x)

    def _cdf(self, x):
        return x

    def _ppf(self, q):
        return q

    def _stats(self):
        return 0.5, 1.0/12, 0, -1.2

    def _entropy(self):
        return 0.0

    @_call_super_mom
    def fit(self, data, *args, **kwds):
        """
        Maximum likelihood estimate for the location and scale parameters.

        `uniform.fit` uses only the following parameters.  Because exact
        formulas are used, the parameters related to optimization that are
        available in the `fit` method of other distributions are ignored
        here.  The only positional argument accepted is `data`.

        Parameters
        ----------
        data : array_like
            Data to use in calculating the maximum likelihood estimate.
        floc : float, optional
            Hold the location parameter fixed to the specified value.
        fscale : float, optional
            Hold the scale parameter fixed to the specified value.

        Returns
        -------
        loc, scale : float
            Maximum likelihood estimates for the location and scale.

        Notes
        -----
        An error is raised if `floc` is given and any values in `data` are
        less than `floc`, or if `fscale` is given and `fscale` is less
        than ``data.max() - data.min()``.  An error is also raised if both
        `floc` and `fscale` are given.

        Examples
        --------
        >>> import numpy as np
        >>> from scipy.stats import uniform

        We'll fit the uniform distribution to `x`:

        >>> x = np.array([2, 2.5, 3.1, 9.5, 13.0])

        For a uniform distribution MLE, the location is the minimum of the
        data, and the scale is the maximum minus the minimum.

        >>> loc, scale = uniform.fit(x)
        >>> loc
        2.0
        >>> scale
        11.0

        If we know the data comes from a uniform distribution where the support
        starts at 0, we can use `floc=0`:

        >>> loc, scale = uniform.fit(x, floc=0)
        >>> loc
        0.0
        >>> scale
        13.0

        Alternatively, if we know the length of the support is 12, we can use
        `fscale=12`:

        >>> loc, scale = uniform.fit(x, fscale=12)
        >>> loc
        1.5
        >>> scale
        12.0

        In that last example, the support interval is [1.5, 13.5].  This
        solution is not unique.  For example, the distribution with ``loc=2``
        and ``scale=12`` has the same likelihood as the one above.  When
        `fscale` is given and it is larger than ``data.max() - data.min()``,
        the parameters returned by the `fit` method center the support over
        the interval ``[data.min(), data.max()]``.

        """
        if len(args) > 0:
            raise TypeError("Too many arguments.")

        floc = kwds.pop('floc', None)
        fscale = kwds.pop('fscale', None)

        _remove_optimizer_parameters(kwds)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        if not np.isfinite(data).all():
            raise ValueError("The data contains non-finite values.")

        # MLE for the uniform distribution
        # --------------------------------
        # The PDF is
        #
        #     f(x, loc, scale) = {1/scale  for loc <= x <= loc + scale
        #                        {0        otherwise}
        #
        # The likelihood function is
        #     L(x, loc, scale) = (1/scale)**n
        # where n is len(x), assuming loc <= x <= loc + scale for all x.
        # The log-likelihood is
        #     l(x, loc, scale) = -n*log(scale)
        # The log-likelihood is maximized by making scale as small as possible,
        # while keeping loc <= x <= loc + scale.   So if neither loc nor scale
        # are fixed, the log-likelihood is maximized by choosing
        #     loc = x.min()
        #     scale = x.ptp()
        # If loc is fixed, it must be less than or equal to x.min(), and then
        # the scale is
        #     scale = x.max() - loc
        # If scale is fixed, it must not be less than x.ptp().  If scale is
        # greater than x.ptp(), the solution is not unique.  Note that the
        # likelihood does not depend on loc, except for the requirement that
        # loc <= x <= loc + scale.  All choices of loc for which
        #     x.max() - scale <= loc <= x.min()
        # have the same log-likelihood.  In this case, we choose loc such that
        # the support is centered over the interval [data.min(), data.max()]:
        #     loc = x.min() = 0.5*(scale - x.ptp())

        if fscale is None:
            # scale is not fixed.
            if floc is None:
                # loc is not fixed, scale is not fixed.
                loc = data.min()
                scale = data.ptp()
            else:
                # loc is fixed, scale is not fixed.
                loc = floc
                scale = data.max() - loc
                if data.min() < loc:
                    raise FitDataError("uniform", lower=loc, upper=loc + scale)
        else:
            # loc is not fixed, scale is fixed.
            ptp = data.ptp()
            if ptp > fscale:
                raise FitUniformFixedScaleDataError(ptp=ptp, fscale=fscale)
            # If ptp < fscale, the ML estimate is not unique; see the comments
            # above.  We choose the distribution for which the support is
            # centered over the interval [data.min(), data.max()].
            loc = data.min() - 0.5*(fscale - ptp)
            scale = fscale

        # We expect the return values to be floating point, so ensure it
        # by explicitly converting to float.
        return float(loc), float(scale)


uniform = uniform_gen(a=0.0, b=1.0, name='uniform')


class vonmises_gen(rv_continuous):
    r"""A Von Mises continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `vonmises` and `vonmises_line` is:

    .. math::

        f(x, \kappa) = \frac{ \exp(\kappa \cos(x)) }{ 2 \pi I_0(\kappa) }

    for :math:`-\pi \le x \le \pi`, :math:`\kappa > 0`. :math:`I_0` is the
    modified Bessel function of order zero (`scipy.special.i0`).

    `vonmises` is a circular distribution which does not restrict the
    distribution to a fixed interval. Currently, there is no circular
    distribution framework in scipy. The ``cdf`` is implemented such that
    ``cdf(x + 2*np.pi) == cdf(x) + 1``.

    `vonmises_line` is the same distribution, defined on :math:`[-\pi, \pi]`
    on the real line. This is a regular (i.e. non-circular) distribution.

    `vonmises` and `vonmises_line` take ``kappa`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("kappa", False, (0, np.inf), (False, False))]

    def _rvs(self, kappa, size=None, random_state=None):
        return random_state.vonmises(0.0, kappa, size=size)

    @inherit_docstring_from(rv_continuous)
    def rvs(self, *args, **kwds):
        rvs = super().rvs(*args, **kwds)
        return np.mod(rvs + np.pi, 2*np.pi) - np.pi

    def _pdf(self, x, kappa):
        # vonmises.pdf(x, kappa) = exp(kappa * cos(x)) / (2*pi*I[0](kappa))
        #                        = exp(kappa * (cos(x) - 1)) /
        #                          (2*pi*exp(-kappa)*I[0](kappa))
        #                        = exp(kappa * cosm1(x)) / (2*pi*i0e(kappa))
        return np.exp(kappa*sc.cosm1(x)) / (2*np.pi*sc.i0e(kappa))

    def _logpdf(self, x, kappa):
        # vonmises.pdf(x, kappa) = exp(kappa * cosm1(x)) / (2*pi*i0e(kappa))
        return kappa * sc.cosm1(x) - np.log(2*np.pi) - np.log(sc.i0e(kappa))

    def _cdf(self, x, kappa):
        return _stats.von_mises_cdf(kappa, x)

    def _stats_skip(self, kappa):
        return 0, None, 0, None

    def _entropy(self, kappa):
        return (-kappa * sc.i1(kappa) / sc.i0(kappa) +
                np.log(2 * np.pi * sc.i0(kappa)))


vonmises = vonmises_gen(name='vonmises')
vonmises_line = vonmises_gen(a=-np.pi, b=np.pi, name='vonmises_line')


class wald_gen(invgauss_gen):
    r"""A Wald continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `wald` is:

    .. math::

        f(x) = \frac{1}{\sqrt{2\pi x^3}} \exp(- \frac{ (x-1)^2 }{ 2x })

    for :math:`x >= 0`.

    `wald` is a special case of `invgauss` with ``mu=1``.

    %(after_notes)s

    %(example)s
    """
    _support_mask = rv_continuous._open_support_mask

    def _shape_info(self):
        return []

    def _rvs(self, size=None, random_state=None):
        return random_state.wald(1.0, 1.0, size=size)

    def _pdf(self, x):
        # wald.pdf(x) = 1/sqrt(2*pi*x**3) * exp(-(x-1)**2/(2*x))
        return invgauss._pdf(x, 1.0)

    def _cdf(self, x):
        return invgauss._cdf(x, 1.0)

    def _sf(self, x):
        return invgauss._sf(x, 1.0)

    def _ppf(self, x):
        return invgauss._ppf(x, 1.0)

    def _isf(self, x):
        return invgauss._isf(x, 1.0)

    def _logpdf(self, x):
        return invgauss._logpdf(x, 1.0)

    def _logcdf(self, x):
        return invgauss._logcdf(x, 1.0)

    def _logsf(self, x):
        return invgauss._logsf(x, 1.0)

    def _stats(self):
        return 1.0, 1.0, 3.0, 15.0


wald = wald_gen(a=0.0, name="wald")


class wrapcauchy_gen(rv_continuous):
    r"""A wrapped Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `wrapcauchy` is:

    .. math::

        f(x, c) = \frac{1-c^2}{2\pi (1+c^2 - 2c \cos(x))}

    for :math:`0 \le x \le 2\pi`, :math:`0 < c < 1`.

    `wrapcauchy` takes ``c`` as a shape parameter for :math:`c`.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        return (c > 0) & (c < 1)

    def _shape_info(self):
        return [_ShapeInfo("c", False, (0, 1), (False, False))]

    def _pdf(self, x, c):
        # wrapcauchy.pdf(x, c) = (1-c**2) / (2*pi*(1+c**2-2*c*cos(x)))
        return (1.0-c*c)/(2*np.pi*(1+c*c-2*c*np.cos(x)))

    def _cdf(self, x, c):

        def f1(x, cr):
            # CDF for 0 <= x < pi
            return 1/np.pi * np.arctan(cr*np.tan(x/2))

        def f2(x, cr):
            # CDF for pi <= x <= 2*pi
            return 1 - 1/np.pi * np.arctan(cr*np.tan((2*np.pi - x)/2))

        cr = (1 + c)/(1 - c)
        return _lazywhere(x < np.pi, (x, cr), f=f1, f2=f2)

    def _ppf(self, q, c):
        val = (1.0-c)/(1.0+c)
        rcq = 2*np.arctan(val*np.tan(np.pi*q))
        rcmq = 2*np.pi-2*np.arctan(val*np.tan(np.pi*(1-q)))
        return np.where(q < 1.0/2, rcq, rcmq)

    def _entropy(self, c):
        return np.log(2*np.pi*(1-c*c))

    def _fitstart(self, data):
        # Use 0.5 as the initial guess of the shape parameter.
        # For the location and scale, use the minimum and
        # peak-to-peak/(2*pi), respectively.
        return 0.5, np.min(data), np.ptp(data)/(2*np.pi)


wrapcauchy = wrapcauchy_gen(a=0.0, b=2*np.pi, name='wrapcauchy')


class gennorm_gen(rv_continuous):
    r"""A generalized normal continuous random variable.

    %(before_notes)s

    See Also
    --------
    laplace : Laplace distribution
    norm : normal distribution

    Notes
    -----
    The probability density function for `gennorm` is [1]_:

    .. math::

        f(x, \beta) = \frac{\beta}{2 \Gamma(1/\beta)} \exp(-|x|^\beta),

    where :math:`x` is a real number, :math:`\beta > 0` and
    :math:`\Gamma` is the gamma function (`scipy.special.gamma`).

    `gennorm` takes ``beta`` as a shape parameter for :math:`\beta`.
    For :math:`\beta = 1`, it is identical to a Laplace distribution.
    For :math:`\beta = 2`, it is identical to a normal distribution
    (with ``scale=1/sqrt(2)``).

    References
    ----------

    .. [1] "Generalized normal distribution, Version 1",
           https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1

    .. [2] Nardon, Martina, and Paolo Pianca. "Simulation techniques for
           generalized Gaussian densities." Journal of Statistical
           Computation and Simulation 79.11 (2009): 1317-1329

    .. [3] Wicklin, Rick. "Simulate data from a generalized Gaussian
           distribution" in The DO Loop blog, September 21, 2016,
           https://blogs.sas.com/content/iml/2016/09/21/simulate-generalized-gaussian-sas.html

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("beta", False, (0, np.inf), (False, False))]

    def _pdf(self, x, beta):
        return np.exp(self._logpdf(x, beta))

    def _logpdf(self, x, beta):
        return np.log(0.5*beta) - sc.gammaln(1.0/beta) - abs(x)**beta

    def _cdf(self, x, beta):
        c = 0.5 * np.sign(x)
        # evaluating (.5 + c) first prevents numerical cancellation
        return (0.5 + c) - c * sc.gammaincc(1.0/beta, abs(x)**beta)

    def _ppf(self, x, beta):
        c = np.sign(x - 0.5)
        # evaluating (1. + c) first prevents numerical cancellation
        return c * sc.gammainccinv(1.0/beta, (1.0 + c) - 2.0*c*x)**(1.0/beta)

    def _sf(self, x, beta):
        return self._cdf(-x, beta)

    def _isf(self, x, beta):
        return -self._ppf(x, beta)

    def _stats(self, beta):
        c1, c3, c5 = sc.gammaln([1.0/beta, 3.0/beta, 5.0/beta])
        return 0., np.exp(c3 - c1), 0., np.exp(c5 + c1 - 2.0*c3) - 3.

    def _entropy(self, beta):
        return 1. / beta - np.log(.5 * beta) + sc.gammaln(1. / beta)

    def _rvs(self, beta, size=None, random_state=None):
        # see [2]_ for the algorithm
        # see [3]_ for reference implementation in SAS
        z = random_state.gamma(1/beta, size=size)
        y = z ** (1/beta)
        # convert y to array to ensure masking support
        y = np.asarray(y)
        mask = random_state.random(size=y.shape) < 0.5
        y[mask] = -y[mask]
        return y


gennorm = gennorm_gen(name='gennorm')


class halfgennorm_gen(rv_continuous):
    r"""The upper half of a generalized normal continuous random variable.

    %(before_notes)s

    See Also
    --------
    gennorm : generalized normal distribution
    expon : exponential distribution
    halfnorm : half normal distribution

    Notes
    -----
    The probability density function for `halfgennorm` is:

    .. math::

        f(x, \beta) = \frac{\beta}{\Gamma(1/\beta)} \exp(-|x|^\beta)

    for :math:`x, \beta > 0`. :math:`\Gamma` is the gamma function
    (`scipy.special.gamma`).

    `halfgennorm` takes ``beta`` as a shape parameter for :math:`\beta`.
    For :math:`\beta = 1`, it is identical to an exponential distribution.
    For :math:`\beta = 2`, it is identical to a half normal distribution
    (with ``scale=1/sqrt(2)``).

    References
    ----------

    .. [1] "Generalized normal distribution, Version 1",
           https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1

    %(example)s

    """
    def _shape_info(self):
        return [_ShapeInfo("beta", False, (0, np.inf), (False, False))]

    def _pdf(self, x, beta):
        #                                 beta
        # halfgennorm.pdf(x, beta) =  -------------  exp(-|x|**beta)
        #                             gamma(1/beta)
        return np.exp(self._logpdf(x, beta))

    def _logpdf(self, x, beta):
        return np.log(beta) - sc.gammaln(1.0/beta) - x**beta

    def _cdf(self, x, beta):
        return sc.gammainc(1.0/beta, x**beta)

    def _ppf(self, x, beta):
        return sc.gammaincinv(1.0/beta, x)**(1.0/beta)

    def _sf(self, x, beta):
        return sc.gammaincc(1.0/beta, x**beta)

    def _isf(self, x, beta):
        return sc.gammainccinv(1.0/beta, x)**(1.0/beta)

    def _entropy(self, beta):
        return 1.0/beta - np.log(beta) + sc.gammaln(1.0/beta)


halfgennorm = halfgennorm_gen(a=0, name='halfgennorm')


class crystalball_gen(rv_continuous):
    r"""
    Crystalball distribution

    %(before_notes)s

    Notes
    -----
    The probability density function for `crystalball` is:

    .. math::

        f(x, \beta, m) =  \begin{cases}
                            N \exp(-x^2 / 2),  &\text{for } x > -\beta\\
                            N A (B - x)^{-m}  &\text{for } x \le -\beta
                          \end{cases}

    where :math:`A = (m / |\beta|)^m  \exp(-\beta^2 / 2)`,
    :math:`B = m/|\beta| - |\beta|` and :math:`N` is a normalisation constant.

    `crystalball` takes :math:`\beta > 0` and :math:`m > 1` as shape
    parameters.  :math:`\beta` defines the point where the pdf changes
    from a power-law to a Gaussian distribution.  :math:`m` is the power
    of the power-law tail.

    References
    ----------
    .. [1] "Crystal Ball Function",
           https://en.wikipedia.org/wiki/Crystal_Ball_function

    %(after_notes)s

    .. versionadded:: 0.19.0

    %(example)s
    """
    def _argcheck(self, beta, m):
        """
        Shape parameter bounds are m > 1 and beta > 0.
        """
        return (m > 1) & (beta > 0)

    def _shape_info(self):
        ibeta = _ShapeInfo("beta", False, (0, np.inf), (False, False))
        im = _ShapeInfo("m", False, (1, np.inf), (False, False))
        return [ibeta, im]

    def _fitstart(self, data):
        # Arbitrary, but the default m=1 is not valid
        return super()._fitstart(data, args=(1, 1.5))

    def _pdf(self, x, beta, m):
        """
        Return PDF of the crystalball function.

                                            --
                                           | exp(-x**2 / 2),  for x > -beta
        crystalball.pdf(x, beta, m) =  N * |
                                           | A * (B - x)**(-m), for x <= -beta
                                            --
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) +
                   _norm_pdf_C * _norm_cdf(beta))

        def rhs(x, beta, m):
            return np.exp(-x**2 / 2)

        def lhs(x, beta, m):
            return ((m/beta)**m * np.exp(-beta**2 / 2.0) *
                    (m/beta - beta - x)**(-m))

        return N * _lazywhere(x > -beta, (x, beta, m), f=rhs, f2=lhs)

    def _logpdf(self, x, beta, m):
        """
        Return the log of the PDF of the crystalball function.
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) +
                   _norm_pdf_C * _norm_cdf(beta))

        def rhs(x, beta, m):
            return -x**2/2

        def lhs(x, beta, m):
            return m*np.log(m/beta) - beta**2/2 - m*np.log(m/beta - beta - x)

        return np.log(N) + _lazywhere(x > -beta, (x, beta, m), f=rhs, f2=lhs)

    def _cdf(self, x, beta, m):
        """
        Return CDF of the crystalball function
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) +
                   _norm_pdf_C * _norm_cdf(beta))

        def rhs(x, beta, m):
            return ((m/beta) * np.exp(-beta**2 / 2.0) / (m-1) +
                    _norm_pdf_C * (_norm_cdf(x) - _norm_cdf(-beta)))

        def lhs(x, beta, m):
            return ((m/beta)**m * np.exp(-beta**2 / 2.0) *
                    (m/beta - beta - x)**(-m+1) / (m-1))

        return N * _lazywhere(x > -beta, (x, beta, m), f=rhs, f2=lhs)

    def _ppf(self, p, beta, m):
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) +
                   _norm_pdf_C * _norm_cdf(beta))
        pbeta = N * (m/beta) * np.exp(-beta**2/2) / (m - 1)

        def ppf_less(p, beta, m):
            eb2 = np.exp(-beta**2/2)
            C = (m/beta) * eb2 / (m-1)
            N = 1/(C + _norm_pdf_C * _norm_cdf(beta))
            return (m/beta - beta -
                    ((m - 1)*(m/beta)**(-m)/eb2*p/N)**(1/(1-m)))

        def ppf_greater(p, beta, m):
            eb2 = np.exp(-beta**2/2)
            C = (m/beta) * eb2 / (m-1)
            N = 1/(C + _norm_pdf_C * _norm_cdf(beta))
            return _norm_ppf(_norm_cdf(-beta) + (1/_norm_pdf_C)*(p/N - C))

        return _lazywhere(p < pbeta, (p, beta, m), f=ppf_less, f2=ppf_greater)

    def _munp(self, n, beta, m):
        """
        Returns the n-th non-central moment of the crystalball function.
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) +
                   _norm_pdf_C * _norm_cdf(beta))

        def n_th_moment(n, beta, m):
            """
            Returns n-th moment. Defined only if n+1 < m
            Function cannot broadcast due to the loop over n
            """
            A = (m/beta)**m * np.exp(-beta**2 / 2.0)
            B = m/beta - beta
            rhs = (2**((n-1)/2.0) * sc.gamma((n+1)/2) *
                   (1.0 + (-1)**n * sc.gammainc((n+1)/2, beta**2 / 2)))
            lhs = np.zeros(rhs.shape)
            for k in range(n + 1):
                lhs += (sc.binom(n, k) * B**(n-k) * (-1)**k / (m - k - 1) *
                        (m/beta)**(-m + k + 1))
            return A * lhs + rhs

        return N * _lazywhere(n + 1 < m, (n, beta, m),
                              np.vectorize(n_th_moment, otypes=[np.float64]),
                              np.inf)


crystalball = crystalball_gen(name='crystalball', longname="A Crystalball Function")


def _argus_phi(chi):
    """
    Utility function for the argus distribution used in the pdf, sf and
    moment calculation.
    Note that for all x > 0:
    gammainc(1.5, x**2/2) = 2 * (_norm_cdf(x) - x * _norm_pdf(x) - 0.5).
    This can be verified directly by noting that the cdf of Gamma(1.5) can
    be written as erf(sqrt(x)) - 2*sqrt(x)*exp(-x)/sqrt(Pi).
    We use gammainc instead of the usual definition because it is more precise
    for small chi.
    """
    return sc.gammainc(1.5, chi**2/2) / 2


class argus_gen(rv_continuous):
    r"""
    Argus distribution

    %(before_notes)s

    Notes
    -----
    The probability density function for `argus` is:

    .. math::

        f(x, \chi) = \frac{\chi^3}{\sqrt{2\pi} \Psi(\chi)} x \sqrt{1-x^2}
                     \exp(-\chi^2 (1 - x^2)/2)

    for :math:`0 < x < 1` and :math:`\chi > 0`, where

    .. math::

        \Psi(\chi) = \Phi(\chi) - \chi \phi(\chi) - 1/2

    with :math:`\Phi` and :math:`\phi` being the CDF and PDF of a standard
    normal distribution, respectively.

    `argus` takes :math:`\chi` as shape a parameter.

    %(after_notes)s

    References
    ----------
    .. [1] "ARGUS distribution",
           https://en.wikipedia.org/wiki/ARGUS_distribution

    .. versionadded:: 0.19.0

    %(example)s
    """
    def _shape_info(self):
        return [_ShapeInfo("chi", False, (0, np.inf), (False, False))]

    def _logpdf(self, x, chi):
        # for x = 0 or 1, logpdf returns -np.inf
        with np.errstate(divide='ignore'):
            y = 1.0 - x*x
            A = 3*np.log(chi) - _norm_pdf_logC - np.log(_argus_phi(chi))
            return A + np.log(x) + 0.5*np.log1p(-x*x) - chi**2 * y / 2

    def _pdf(self, x, chi):
        return np.exp(self._logpdf(x, chi))

    def _cdf(self, x, chi):
        return 1.0 - self._sf(x, chi)

    def _sf(self, x, chi):
        return _argus_phi(chi * np.sqrt(1 - x**2)) / _argus_phi(chi)

    def _rvs(self, chi, size=None, random_state=None):
        chi = np.asarray(chi)
        if chi.size == 1:
            out = self._rvs_scalar(chi, numsamples=size,
                                   random_state=random_state)
        else:
            shp, bc = _check_shape(chi.shape, size)
            numsamples = int(np.prod(shp))
            out = np.empty(size)
            it = np.nditer([chi],
                           flags=['multi_index'],
                           op_flags=[['readonly']])
            while not it.finished:
                idx = tuple((it.multi_index[j] if not bc[j] else slice(None))
                            for j in range(-len(size), 0))
                r = self._rvs_scalar(it[0], numsamples=numsamples,
                                     random_state=random_state)
                out[idx] = r.reshape(shp)
                it.iternext()

        if size == ():
            out = out[()]
        return out

    def _rvs_scalar(self, chi, numsamples=None, random_state=None):
        # if chi <= 1.8:
        # use rejection method, see Devroye:
        # Non-Uniform Random Variate Generation, 1986, section II.3.2.
        # write: PDF f(x) = c * g(x) * h(x), where
        # h is [0,1]-valued and g is a density
        # we use two ways to write f
        #
        # Case 1:
        # write g(x) = 3*x*sqrt(1-x**2), h(x) = exp(-chi**2 (1-x**2) / 2)
        # If X has a distribution with density g its ppf G_inv is given by:
        # G_inv(u) = np.sqrt(1 - u**(2/3))
        #
        # Case 2:
        # g(x) = chi**2 * x * exp(-chi**2 * (1-x**2)/2) / (1 - exp(-chi**2 /2))
        # h(x) = sqrt(1 - x**2), 0 <= x <= 1
        # one can show that
        # G_inv(u) = np.sqrt(2*np.log(u*(np.exp(chi**2/2)-1)+1))/chi
        #          = np.sqrt(1 + 2*np.log(np.exp(-chi**2/2)*(1-u)+u)/chi**2)
        # the latter expression is used for precision with small chi
        #
        # In both cases, the inverse cdf of g can be written analytically, and
        # we can apply the rejection method:
        #
        # REPEAT
        #    Generate U uniformly distributed on [0, 1]
        #    Generate X with density g (e.g. via inverse transform sampling:
        #    X = G_inv(V) with V uniformly distributed on [0, 1])
        # UNTIL X <= h(X)
        # RETURN X
        #
        # We use case 1 for chi <= 0.5 as it maintains precision for small chi
        # and case 2 for 0.5 < chi <= 1.8 due to its speed for moderate chi.
        #
        # if chi > 1.8:
        # use relation to the Gamma distribution: if X is ARGUS with parameter
        # chi), then Y = chi**2 * (1 - X**2) / 2 has density proportional to
        # sqrt(u) * exp(-u) on [0, chi**2 / 2], i.e. a Gamma(3/2) distribution
        # conditioned on [0, chi**2 / 2]). Therefore, to sample X from the
        # ARGUS distribution, we sample Y from the gamma distribution, keeping
        # only samples on [0, chi**2 / 2], and apply the inverse
        # transformation X = (1 - 2*Y/chi**2)**(1/2). Since we only
        # look at chi > 1.8, gamma(1.5).cdf(chi**2/2) is large enough such
        # Y falls in the inteval [0, chi**2 / 2] with a high probability:
        # stats.gamma(1.5).cdf(1.8**2/2) = 0.644...
        #
        # The points to switch between the different methods are determined
        # by a comparison of the runtime of the different methods. However,
        # the runtime is platform-dependent. The implemented values should
        # ensure a good overall performance and are supported by an analysis
        # of the rejection constants of different methods.

        size1d = tuple(np.atleast_1d(numsamples))
        N = int(np.prod(size1d))
        x = np.zeros(N)
        simulated = 0
        chi2 = chi * chi
        if chi <= 0.5:
            d = -chi2 / 2
            while simulated < N:
                k = N - simulated
                u = random_state.uniform(size=k)
                v = random_state.uniform(size=k)
                z = v**(2/3)
                # acceptance condition: u <= h(G_inv(v)). This simplifies to
                accept = (np.log(u) <= d * z)
                num_accept = np.sum(accept)
                if num_accept > 0:
                    # we still need to transform z=v**(2/3) to X = G_inv(v)
                    rvs = np.sqrt(1 - z[accept])
                    x[simulated:(simulated + num_accept)] = rvs
                    simulated += num_accept
        elif chi <= 1.8:
            echi = np.exp(-chi2 / 2)
            while simulated < N:
                k = N - simulated
                u = random_state.uniform(size=k)
                v = random_state.uniform(size=k)
                z = 2 * np.log(echi * (1 - v) + v) / chi2
                # as in case one, simplify u <= h(G_inv(v)) and then transform
                # z to the target distribution X = G_inv(v)
                accept = (u**2 + z <= 0)
                num_accept = np.sum(accept)
                if num_accept > 0:
                    rvs = np.sqrt(1 + z[accept])
                    x[simulated:(simulated + num_accept)] = rvs
                    simulated += num_accept
        else:
            # conditional Gamma for chi > 1.8
            while simulated < N:
                k = N - simulated
                g = random_state.standard_gamma(1.5, size=k)
                accept = (g <= chi2 / 2)
                num_accept = np.sum(accept)
                if num_accept > 0:
                    x[simulated:(simulated + num_accept)] = g[accept]
                    simulated += num_accept
            x = np.sqrt(1 - 2 * x / chi2)

        return np.reshape(x, size1d)

    def _stats(self, chi):
        # need to ensure that dtype is float
        # otherwise the mask below does not work for integers
        chi = np.asarray(chi, dtype=float)
        phi = _argus_phi(chi)
        m = np.sqrt(np.pi/8) * chi * sc.ive(1, chi**2/4) / phi
        # compute second moment, use Taylor expansion for small chi (<= 0.1)
        mu2 = np.empty_like(chi)
        mask = chi > 0.1
        c = chi[mask]
        mu2[mask] = 1 - 3 / c**2 + c * _norm_pdf(c) / phi[mask]
        c = chi[~mask]
        coef = [-358/65690625, 0, -94/1010625, 0, 2/2625, 0, 6/175, 0, 0.4]
        mu2[~mask] = np.polyval(coef, c)
        return m, mu2 - m**2, None, None


argus = argus_gen(name='argus', longname="An Argus Function", a=0.0, b=1.0)


class rv_histogram(rv_continuous):
    """
    Generates a distribution given by a histogram.
    This is useful to generate a template distribution from a binned
    datasample.

    As a subclass of the `rv_continuous` class, `rv_histogram` inherits from it
    a collection of generic methods (see `rv_continuous` for the full list),
    and implements them based on the properties of the provided binned
    datasample.

    Parameters
    ----------
    histogram : tuple of array_like
        Tuple containing two array_like objects.
        The first containing the content of n bins,
        the second containing the (n+1) bin boundaries.
        In particular, the return value of `numpy.histogram` is accepted.

    density : bool, optional
        If False, assumes the histogram is proportional to counts per bin;
        otherwise, assumes it is proportional to a density.
        For constant bin widths, these are equivalent, but the distinction
        is important when bin widths vary (see Notes).
        If None (default), sets ``density=True`` for backwards compatibility,
        but warns if the bin widths are variable. Set `density` explicitly
        to silence the warning.

        .. versionadded:: 1.10.0

    Notes
    -----
    When a histogram has unequal bin widths, there is a distinction between
    histograms that are proportional to counts per bin and histograms that are
    proportional to probability density over a bin. If `numpy.histogram` is
    called with its default ``density=False``, the resulting histogram is the
    number of counts per bin, so ``density=False`` should be passed to
    `rv_histogram`. If `numpy.histogram` is called with ``density=True``, the
    resulting histogram is in terms of probability density, so ``density=True``
    should be passed to `rv_histogram`. To avoid warnings, always pass
    ``density`` explicitly when the input histogram has unequal bin widths.

    There are no additional shape parameters except for the loc and scale.
    The pdf is defined as a stepwise function from the provided histogram.
    The cdf is a linear interpolation of the pdf.

    .. versionadded:: 0.19.0

    Examples
    --------

    Create a scipy.stats distribution from a numpy histogram

    >>> import scipy.stats
    >>> import numpy as np
    >>> data = scipy.stats.norm.rvs(size=100000, loc=0, scale=1.5, random_state=123)
    >>> hist = np.histogram(data, bins=100)
    >>> hist_dist = scipy.stats.rv_histogram(hist, density=False)

    Behaves like an ordinary scipy rv_continuous distribution

    >>> hist_dist.pdf(1.0)
    0.20538577847618705
    >>> hist_dist.cdf(2.0)
    0.90818568543056499

    PDF is zero above (below) the highest (lowest) bin of the histogram,
    defined by the max (min) of the original dataset

    >>> hist_dist.pdf(np.max(data))
    0.0
    >>> hist_dist.cdf(np.max(data))
    1.0
    >>> hist_dist.pdf(np.min(data))
    7.7591907244498314e-05
    >>> hist_dist.cdf(np.min(data))
    0.0

    PDF and CDF follow the histogram

    >>> import matplotlib.pyplot as plt
    >>> X = np.linspace(-5.0, 5.0, 100)
    >>> fig, ax = plt.subplots()
    >>> ax.set_title("PDF from Template")
    >>> ax.hist(data, density=True, bins=100)
    >>> ax.plot(X, hist_dist.pdf(X), label='PDF')
    >>> ax.plot(X, hist_dist.cdf(X), label='CDF')
    >>> ax.legend()
    >>> fig.show()

    """
    _support_mask = rv_continuous._support_mask

    def __init__(self, histogram, *args, density=None, **kwargs):
        """
        Create a new distribution using the given histogram

        Parameters
        ----------
        histogram : tuple of array_like
            Tuple containing two array_like objects.
            The first containing the content of n bins,
            the second containing the (n+1) bin boundaries.
            In particular, the return value of np.histogram is accepted.
        density : bool, optional
            If False, assumes the histogram is proportional to counts per bin;
            otherwise, assumes it is proportional to a density.
            For constant bin widths, these are equivalent.
            If None (default), sets ``density=True`` for backward
            compatibility, but warns if the bin widths are variable. Set
            `density` explicitly to silence the warning.
        """
        self._histogram = histogram
        self._density = density
        if len(histogram) != 2:
            raise ValueError("Expected length 2 for parameter histogram")
        self._hpdf = np.asarray(histogram[0])
        self._hbins = np.asarray(histogram[1])
        if len(self._hpdf) + 1 != len(self._hbins):
            raise ValueError("Number of elements in histogram content "
                             "and histogram boundaries do not match, "
                             "expected n and n+1.")
        self._hbin_widths = self._hbins[1:] - self._hbins[:-1]
        bins_vary = not np.allclose(self._hbin_widths, self._hbin_widths[0])
        if density is None and bins_vary:
            message = ("Bin widths are not constant. Assuming `density=True`."
                       "Specify `density` explicitly to silence this warning.")
            warnings.warn(message, RuntimeWarning, stacklevel=2)
            density = True
        elif not density:
            self._hpdf = self._hpdf / self._hbin_widths

        self._hpdf = self._hpdf / float(np.sum(self._hpdf * self._hbin_widths))
        self._hcdf = np.cumsum(self._hpdf * self._hbin_widths)
        self._hpdf = np.hstack([0.0, self._hpdf, 0.0])
        self._hcdf = np.hstack([0.0, self._hcdf])
        # Set support
        kwargs['a'] = self.a = self._hbins[0]
        kwargs['b'] = self.b = self._hbins[-1]
        super().__init__(*args, **kwargs)

    def _pdf(self, x):
        """
        PDF of the histogram
        """
        return self._hpdf[np.searchsorted(self._hbins, x, side='right')]

    def _cdf(self, x):
        """
        CDF calculated from the histogram
        """
        return np.interp(x, self._hbins, self._hcdf)

    def _ppf(self, x):
        """
        Percentile function calculated from the histogram
        """
        return np.interp(x, self._hcdf, self._hbins)

    def _munp(self, n):
        """Compute the n-th non-central moment."""
        integrals = (self._hbins[1:]**(n+1) - self._hbins[:-1]**(n+1)) / (n+1)
        return np.sum(self._hpdf[1:-1] * integrals)

    def _entropy(self):
        """Compute entropy of distribution"""
        res = _lazywhere(self._hpdf[1:-1] > 0.0,
                         (self._hpdf[1:-1],),
                         np.log,
                         0.0)
        return -np.sum(self._hpdf[1:-1] * res * self._hbin_widths)

    def _updated_ctor_param(self):
        """
        Set the histogram as additional constructor argument
        """
        dct = super()._updated_ctor_param()
        dct['histogram'] = self._histogram
        dct['density'] = self._density
        return dct


class studentized_range_gen(rv_continuous):
    r"""A studentized range continuous random variable.

    %(before_notes)s

    See Also
    --------
    t: Student's t distribution

    Notes
    -----
    The probability density function for `studentized_range` is:

    .. math::

         f(x; k, \nu) = \frac{k(k-1)\nu^{\nu/2}}{\Gamma(\nu/2)
                        2^{\nu/2-1}} \int_{0}^{\infty} \int_{-\infty}^{\infty}
                        s^{\nu} e^{-\nu s^2/2} \phi(z) \phi(sx + z)
                        [\Phi(sx + z) - \Phi(z)]^{k-2} \,dz \,ds

    for :math:`x ≥ 0`, :math:`k > 1`, and :math:`\nu > 0`.

    `studentized_range` takes ``k`` for :math:`k` and ``df`` for :math:`\nu`
    as shape parameters.

    When :math:`\nu` exceeds 100,000, an asymptotic approximation (infinite
    degrees of freedom) is used to compute the cumulative distribution
    function [4]_ and probability distribution function.

    %(after_notes)s

    References
    ----------

    .. [1] "Studentized range distribution",
           https://en.wikipedia.org/wiki/Studentized_range_distribution
    .. [2] Batista, Ben Dêivide, et al. "Externally Studentized Normal Midrange
           Distribution." Ciência e Agrotecnologia, vol. 41, no. 4, 2017, pp.
           378-389., doi:10.1590/1413-70542017414047716.
    .. [3] Harter, H. Leon. "Tables of Range and Studentized Range." The Annals
           of Mathematical Statistics, vol. 31, no. 4, 1960, pp. 1122-1147.
           JSTOR, www.jstor.org/stable/2237810. Accessed 18 Feb. 2021.
    .. [4] Lund, R. E., and J. R. Lund. "Algorithm AS 190: Probabilities and
           Upper Quantiles for the Studentized Range." Journal of the Royal
           Statistical Society. Series C (Applied Statistics), vol. 32, no. 2,
           1983, pp. 204-210. JSTOR, www.jstor.org/stable/2347300. Accessed 18
           Feb. 2021.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import studentized_range
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 1)

    Calculate the first four moments:

    >>> k, df = 3, 10
    >>> mean, var, skew, kurt = studentized_range.stats(k, df, moments='mvsk')

    Display the probability density function (``pdf``):

    >>> x = np.linspace(studentized_range.ppf(0.01, k, df),
    ...                 studentized_range.ppf(0.99, k, df), 100)
    >>> ax.plot(x, studentized_range.pdf(x, k, df),
    ...         'r-', lw=5, alpha=0.6, label='studentized_range pdf')

    Alternatively, the distribution object can be called (as a function)
    to fix the shape, location and scale parameters. This returns a "frozen"
    RV object holding the given parameters fixed.

    Freeze the distribution and display the frozen ``pdf``:

    >>> rv = studentized_range(k, df)
    >>> ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

    Check accuracy of ``cdf`` and ``ppf``:

    >>> vals = studentized_range.ppf([0.001, 0.5, 0.999], k, df)
    >>> np.allclose([0.001, 0.5, 0.999], studentized_range.cdf(vals, k, df))
    True

    Rather than using (``studentized_range.rvs``) to generate random variates,
    which is very slow for this distribution, we can approximate the inverse
    CDF using an interpolator, and then perform inverse transform sampling
    with this approximate inverse CDF.

    This distribution has an infinite but thin right tail, so we focus our
    attention on the leftmost 99.9 percent.

    >>> a, b = studentized_range.ppf([0, .999], k, df)
    >>> a, b
    0, 7.41058083802274

    >>> from scipy.interpolate import interp1d
    >>> rng = np.random.default_rng()
    >>> xs = np.linspace(a, b, 50)
    >>> cdf = studentized_range.cdf(xs, k, df)
    # Create an interpolant of the inverse CDF
    >>> ppf = interp1d(cdf, xs, fill_value='extrapolate')
    # Perform inverse transform sampling using the interpolant
    >>> r = ppf(rng.uniform(size=1000))

    And compare the histogram:

    >>> ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
    >>> ax.legend(loc='best', frameon=False)
    >>> plt.show()

    """

    def _argcheck(self, k, df):
        return (k > 1) & (df > 0)

    def _shape_info(self):
        ik = _ShapeInfo("k", False, (1, np.inf), (False, False))
        idf = _ShapeInfo("df", False, (0, np.inf), (False, False))
        return [ik, idf]

    def _fitstart(self, data):
        # Default is k=1, but that is not a valid value of the parameter.
        return super(studentized_range_gen, self)._fitstart(data, args=(2, 1))

    def _munp(self, K, k, df):
        cython_symbol = '_studentized_range_moment'
        _a, _b = self._get_support()
        # all three of these are used to create a numpy array so they must
        # be the same shape.

        def _single_moment(K, k, df):
            log_const = _stats._studentized_range_pdf_logconst(k, df)
            arg = [K, k, df, log_const]
            usr_data = np.array(arg, float).ctypes.data_as(ctypes.c_void_p)

            llc = LowLevelCallable.from_cython(_stats, cython_symbol, usr_data)

            ranges = [(-np.inf, np.inf), (0, np.inf), (_a, _b)]
            opts = dict(epsabs=1e-11, epsrel=1e-12)

            return integrate.nquad(llc, ranges=ranges, opts=opts)[0]

        ufunc = np.frompyfunc(_single_moment, 3, 1)
        return np.float64(ufunc(K, k, df))

    def _pdf(self, x, k, df):

        def _single_pdf(q, k, df):
            # The infinite form of the PDF is derived from the infinite
            # CDF.
            if df < 100000:
                cython_symbol = '_studentized_range_pdf'
                log_const = _stats._studentized_range_pdf_logconst(k, df)
                arg = [q, k, df, log_const]
                usr_data = np.array(arg, float).ctypes.data_as(ctypes.c_void_p)
                ranges = [(-np.inf, np.inf), (0, np.inf)]

            else:
                cython_symbol = '_studentized_range_pdf_asymptotic'
                arg = [q, k]
                usr_data = np.array(arg, float).ctypes.data_as(ctypes.c_void_p)
                ranges = [(-np.inf, np.inf)]

            llc = LowLevelCallable.from_cython(_stats, cython_symbol, usr_data)
            opts = dict(epsabs=1e-11, epsrel=1e-12)
            return integrate.nquad(llc, ranges=ranges, opts=opts)[0]

        ufunc = np.frompyfunc(_single_pdf, 3, 1)
        return np.float64(ufunc(x, k, df))

    def _cdf(self, x, k, df):

        def _single_cdf(q, k, df):
            # "When the degrees of freedom V are infinite the probability
            # integral takes [on a] simpler form," and a single asymptotic
            # integral is evaluated rather than the standard double integral.
            # (Lund, Lund, page 205)
            if df < 100000:
                cython_symbol = '_studentized_range_cdf'
                log_const = _stats._studentized_range_cdf_logconst(k, df)
                arg = [q, k, df, log_const]
                usr_data = np.array(arg, float).ctypes.data_as(ctypes.c_void_p)
                ranges = [(-np.inf, np.inf), (0, np.inf)]

            else:
                cython_symbol = '_studentized_range_cdf_asymptotic'
                arg = [q, k]
                usr_data = np.array(arg, float).ctypes.data_as(ctypes.c_void_p)
                ranges = [(-np.inf, np.inf)]

            llc = LowLevelCallable.from_cython(_stats, cython_symbol, usr_data)
            opts = dict(epsabs=1e-11, epsrel=1e-12)
            return integrate.nquad(llc, ranges=ranges, opts=opts)[0]

        ufunc = np.frompyfunc(_single_cdf, 3, 1)

        # clip p-values to ensure they are in [0, 1].
        return np.clip(np.float64(ufunc(x, k, df)), 0, 1)


studentized_range = studentized_range_gen(name='studentized_range', a=0,
                                          b=np.inf)


# Collect names of classes and objects in this module.
pairs = list(globals().copy().items())
_distn_names, _distn_gen_names = get_distribution_names(pairs, rv_continuous)

__all__ = _distn_names + _distn_gen_names + ['rv_histogram']
