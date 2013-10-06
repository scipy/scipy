#
# Author:  Travis Oliphant  2002-2011 with contributions from
#          SciPy Developers 2004-2011
#
from __future__ import division, print_function, absolute_import

import math
import sys
import warnings

from scipy.lib.six import callable, string_types, get_method_function
from scipy.lib.six import exec_

from scipy.misc import comb, derivative
from scipy.misc.doccer import inherit_docstring_from
from scipy import special
from scipy import optimize
from scipy import integrate
from scipy.special import gammaln as gamln

import keyword
import re
import inspect
from numpy import all, where, arange, putmask, \
     ravel, take, ones, sum, shape, product, reshape, \
     zeros, floor, logical_and, log, sqrt, exp, arctanh, tan, sin, arcsin, \
     arctan, tanh, ndarray, cos, cosh, sinh, newaxis, log1p, expm1
from numpy import atleast_1d, polyval, ceil, place, extract, \
     any, argsort, argmax, vectorize, r_, asarray, nan, inf, pi, isinf, \
     NINF, empty
import numpy
import numpy as np
import numpy.random as mtrand
from numpy import flatnonzero as nonzero
from . import vonmises_cython
from ._tukeylambda_stats import tukeylambda_variance as _tlvar, \
                                tukeylambda_kurtosis as _tlkurt

__all__ = [
           'rv_continuous',
           'ksone', 'kstwobign', 'norm', 'alpha', 'anglit', 'arcsine',
           'beta', 'betaprime', 'bradford', 'burr', 'fisk', 'cauchy',
           'chi', 'chi2', 'cosine', 'dgamma', 'dweibull', 'erlang',
           'expon', 'exponweib', 'exponpow', 'fatiguelife', 'foldcauchy',
           'f', 'foldnorm', 'frechet_r', 'weibull_min', 'frechet_l',
           'weibull_max', 'genlogistic', 'genpareto', 'genexpon', 'genextreme',
           'gamma', 'gengamma', 'genhalflogistic', 'gompertz', 'gumbel_r',
           'gumbel_l', 'halfcauchy', 'halflogistic', 'halfnorm', 'hypsecant',
           'gausshyper', 'invgamma', 'invgauss', 'invweibull',
           'johnsonsb', 'johnsonsu', 'laplace', 'levy', 'levy_l',
           'levy_stable', 'logistic', 'loggamma', 'loglaplace', 'lognorm',
           'gilbrat', 'maxwell', 'mielke', 'nakagami', 'ncx2', 'ncf', 't',
           'nct', 'pareto', 'lomax', 'pearson3', 'powerlaw', 'powerlognorm',
           'powernorm', 'rdist', 'rayleigh', 'reciprocal', 'rice',
           'recipinvgauss', 'semicircular', 'triang', 'truncexpon',
           'truncnorm', 'tukeylambda', 'uniform', 'vonmises', 'wald',
           'wrapcauchy', 'entropy', 'rv_discrete', 'binom', 'bernoulli',
           'nbinom', 'geom', 'hypergeom', 'logser', 'poisson', 'planck',
           'boltzmann', 'randint', 'zipf', 'dlaplace', 'skellam'
          ]

floatinfo = numpy.finfo(float)
eps = numpy.finfo(float).eps

gam = special.gamma
random = mtrand.random_sample

import types
from scipy.misc import doccer

try:
    from new import instancemethod
except ImportError:
    # Python 3
    def instancemethod(func, obj, cls):
        return types.MethodType(func, obj)


# These are the docstring parts used for substitution in specific
# distribution docstrings

docheaders = {'methods':"""\nMethods\n-------\n""",
              'parameters':"""\nParameters\n---------\n""",
              'notes':"""\nNotes\n-----\n""",
              'examples':"""\nExamples\n--------\n"""}

_doc_rvs = \
"""rvs(%(shapes)s, loc=0, scale=1, size=1)
    Random variates.
"""
_doc_pdf = \
"""pdf(x, %(shapes)s, loc=0, scale=1)
    Probability density function.
"""
_doc_logpdf = \
"""logpdf(x, %(shapes)s, loc=0, scale=1)
    Log of the probability density function.
"""
_doc_pmf = \
"""pmf(x, %(shapes)s, loc=0, scale=1)
    Probability mass function.
"""
_doc_logpmf = \
"""logpmf(x, %(shapes)s, loc=0, scale=1)
    Log of the probability mass function.
"""
_doc_cdf = \
"""cdf(x, %(shapes)s, loc=0, scale=1)
    Cumulative density function.
"""
_doc_logcdf = \
"""logcdf(x, %(shapes)s, loc=0, scale=1)
    Log of the cumulative density function.
"""
_doc_sf = \
"""sf(x, %(shapes)s, loc=0, scale=1)
    Survival function (1-cdf --- sometimes more accurate).
"""
_doc_logsf = \
"""logsf(x, %(shapes)s, loc=0, scale=1)
    Log of the survival function.
"""
_doc_ppf = \
"""ppf(q, %(shapes)s, loc=0, scale=1)
    Percent point function (inverse of cdf --- percentiles).
"""
_doc_isf = \
"""isf(q, %(shapes)s, loc=0, scale=1)
    Inverse survival function (inverse of sf).
"""
_doc_moment = \
"""moment(n, %(shapes)s, loc=0, scale=1)
    Non-central moment of order n
"""
_doc_stats = \
"""stats(%(shapes)s, loc=0, scale=1, moments='mv')
    Mean('m'), variance('v'), skew('s'), and/or kurtosis('k').
"""
_doc_entropy = \
"""entropy(%(shapes)s, loc=0, scale=1)
    (Differential) entropy of the RV.
"""
_doc_fit = \
"""fit(data, %(shapes)s, loc=0, scale=1)
    Parameter estimates for generic data.
"""
_doc_expect = \
"""expect(func, %(shapes)s, loc=0, scale=1, lb=None, ub=None, conditional=False, **kwds)
    Expected value of a function (of one argument) with respect to the distribution.
"""
_doc_expect_discrete = \
"""expect(func, %(shapes)s, loc=0, lb=None, ub=None, conditional=False)
    Expected value of a function (of one argument) with respect to the distribution.
"""
_doc_median = \
"""median(%(shapes)s, loc=0, scale=1)
    Median of the distribution.
"""
_doc_mean = \
"""mean(%(shapes)s, loc=0, scale=1)
    Mean of the distribution.
"""
_doc_var = \
"""var(%(shapes)s, loc=0, scale=1)
    Variance of the distribution.
"""
_doc_std = \
"""std(%(shapes)s, loc=0, scale=1)
    Standard deviation of the distribution.
"""
_doc_interval = \
"""interval(alpha, %(shapes)s, loc=0, scale=1)
    Endpoints of the range that contains alpha percent of the distribution
"""
_doc_allmethods = ''.join([docheaders['methods'], _doc_rvs, _doc_pdf,
                           _doc_logpdf, _doc_cdf, _doc_logcdf, _doc_sf,
                           _doc_logsf, _doc_ppf, _doc_isf, _doc_moment,
                           _doc_stats, _doc_entropy, _doc_fit,
                           _doc_expect, _doc_median,
                           _doc_mean, _doc_var, _doc_std, _doc_interval])

# Note that the two lines for %(shapes) are searched for and replaced in
# rv_continuous and rv_discrete - update there if the exact string changes
_doc_default_callparams = \
"""
Parameters
----------
x : array_like
    quantiles
q : array_like
    lower or upper tail probability
%(shapes)s : array_like
    shape parameters
loc : array_like, optional
    location parameter (default=0)
scale : array_like, optional
    scale parameter (default=1)
size : int or tuple of ints, optional
    shape of random variates (default computed from input arguments )
moments : str, optional
    composed of letters ['mvsk'] specifying which moments to compute where
    'm' = mean, 'v' = variance, 's' = (Fisher's) skew and
    'k' = (Fisher's) kurtosis. (default='mv')
"""
_doc_default_longsummary = \
"""Continuous random variables are defined from a standard form and may
require some shape parameters to complete its specification.  Any
optional keyword parameters can be passed to the methods of the RV
object as given below:
"""
_doc_default_frozen_note = \
"""
Alternatively, the object may be called (as a function) to fix the shape,
location, and scale parameters returning a "frozen" continuous RV object:

rv = %(name)s(%(shapes)s, loc=0, scale=1)
    - Frozen RV object with the same methods but holding the given shape,
      location, and scale fixed.
"""
_doc_default_example = \
"""Examples
--------
>>> from scipy.stats import %(name)s
>>> numargs = %(name)s.numargs
>>> [ %(shapes)s ] = [0.9,] * numargs
>>> rv = %(name)s(%(shapes)s)

Display frozen pdf

>>> x = np.linspace(0, np.minimum(rv.dist.b, 3))
>>> h = plt.plot(x, rv.pdf(x))

Here, ``rv.dist.b`` is the right endpoint of the support of ``rv.dist``.

Check accuracy of cdf and ppf

>>> prb = %(name)s.cdf(x, %(shapes)s)
>>> h = plt.semilogy(np.abs(x - %(name)s.ppf(prb, %(shapes)s)) + 1e-20)

Random number generation

>>> R = %(name)s.rvs(%(shapes)s, size=100)
"""

_doc_default = ''.join([_doc_default_longsummary,
                        _doc_allmethods,
                        _doc_default_callparams,
                        _doc_default_frozen_note,
                        _doc_default_example])

_doc_default_before_notes = ''.join([_doc_default_longsummary,
                                     _doc_allmethods,
                                     _doc_default_callparams,
                                     _doc_default_frozen_note])

docdict = {'rvs':_doc_rvs,
           'pdf':_doc_pdf,
           'logpdf':_doc_logpdf,
           'cdf':_doc_cdf,
           'logcdf':_doc_logcdf,
           'sf':_doc_sf,
           'logsf':_doc_logsf,
           'ppf':_doc_ppf,
           'isf':_doc_isf,
           'stats':_doc_stats,
           'entropy':_doc_entropy,
           'fit':_doc_fit,
           'moment':_doc_moment,
           'expect':_doc_expect,
           'interval':_doc_interval,
           'mean':_doc_mean,
           'std':_doc_std,
           'var':_doc_var,
           'median':_doc_median,
           'allmethods':_doc_allmethods,
           'callparams':_doc_default_callparams,
           'longsummary':_doc_default_longsummary,
           'frozennote':_doc_default_frozen_note,
           'example':_doc_default_example,
           'default':_doc_default,
           'before_notes':_doc_default_before_notes}

# Reuse common content between continous and discrete docs, change some
# minor bits.
docdict_discrete = docdict.copy()

docdict_discrete['pmf'] = _doc_pmf
docdict_discrete['logpmf'] = _doc_logpmf
docdict_discrete['expect'] = _doc_expect_discrete
_doc_disc_methods = ['rvs', 'pmf', 'logpmf', 'cdf', 'logcdf', 'sf', 'logsf',
                     'ppf', 'isf', 'stats', 'entropy', 'expect', 'median',
                     'mean', 'var', 'std', 'interval']
for obj in _doc_disc_methods:
    docdict_discrete[obj] = docdict_discrete[obj].replace(', scale=1', '')
docdict_discrete.pop('pdf')
docdict_discrete.pop('logpdf')

_doc_allmethods = ''.join([docdict_discrete[obj] for obj in
                                              _doc_disc_methods])
docdict_discrete['allmethods'] = docheaders['methods'] + _doc_allmethods

docdict_discrete['longsummary'] = _doc_default_longsummary.replace(
                                      'Continuous', 'Discrete')
_doc_default_frozen_note = \
"""
Alternatively, the object may be called (as a function) to fix the shape and
location parameters returning a "frozen" discrete RV object:

rv = %(name)s(%(shapes)s, loc=0)
    - Frozen RV object with the same methods but holding the given shape and
      location fixed.
"""
docdict_discrete['frozennote'] = _doc_default_frozen_note

_doc_default_discrete_example = \
"""Examples
--------
>>> from scipy.stats import %(name)s
>>> [ %(shapes)s ] = [<Replace with reasonable values>]
>>> rv = %(name)s(%(shapes)s)

Display frozen pmf

>>> x = np.arange(0, np.minimum(rv.dist.b, 3))
>>> h = plt.vlines(x, 0, rv.pmf(x), lw=2)

Here, ``rv.dist.b`` is the right endpoint of the support of ``rv.dist``.

Check accuracy of cdf and ppf

>>> prb = %(name)s.cdf(x, %(shapes)s)
>>> h = plt.semilogy(np.abs(x - %(name)s.ppf(prb, %(shapes)s)) + 1e-20)

Random number generation

>>> R = %(name)s.rvs(%(shapes)s, size=100)

"""
docdict_discrete['example'] = _doc_default_discrete_example

_doc_default_before_notes = ''.join([docdict_discrete['longsummary'],
                                     docdict_discrete['allmethods'],
                                     docdict_discrete['callparams'],
                                     docdict_discrete['frozennote']])
docdict_discrete['before_notes'] = _doc_default_before_notes

_doc_default_disc = ''.join([docdict_discrete['longsummary'],
                             docdict_discrete['allmethods'],
                             docdict_discrete['frozennote'],
                             docdict_discrete['example']])
docdict_discrete['default'] = _doc_default_disc


# clean up all the separate docstring elements, we do not need them anymore
for obj in [s for s in dir() if s.startswith('_doc_')]:
    exec('del ' + obj)
del obj
try:
    del s
except NameError:
    # in Python 3, loop variables are not visible after the loop
    pass


def _moment(data, n, mu=None):
    if mu is None:
        mu = data.mean()
    return ((data - mu)**n).mean()


def _moment_from_stats(n, mu, mu2, g1, g2, moment_func, args):
    if (n == 0):
        return 1.0
    elif (n == 1):
        if mu is None:
            val = moment_func(1,*args)
        else:
            val = mu
    elif (n == 2):
        if mu2 is None or mu is None:
            val = moment_func(2,*args)
        else:
            val = mu2 + mu*mu
    elif (n == 3):
        if g1 is None or mu2 is None or mu is None:
            val = moment_func(3,*args)
        else:
            mu3 = g1 * np.power(mu2, 1.5)  # 3rd central moment
            val = mu3+3*mu*mu2+mu*mu*mu  # 3rd non-central moment
    elif (n == 4):
        if g1 is None or g2 is None or mu2 is None or mu is None:
            val = moment_func(4,*args)
        else:
            mu4 = (g2+3.0)*(mu2**2.0)  # 4th central moment
            mu3 = g1*np.power(mu2, 1.5)  # 3rd central moment
            val = mu4+4*mu*mu3+6*mu*mu*mu2+mu*mu*mu*mu
    else:
        val = moment_func(n, *args)

    return val


def _skew(data):
    """
    skew is third central moment / variance**(1.5)
    """
    data = np.ravel(data)
    mu = data.mean()
    m2 = ((data - mu)**2).mean()
    m3 = ((data - mu)**3).mean()
    return m3 / np.power(m2, 1.5)


def _kurtosis(data):
    """
    kurtosis is fourth central moment / variance**2 - 3
    """
    data = np.ravel(data)
    mu = data.mean()
    m2 = ((data - mu)**2).mean()
    m4 = ((data - mu)**4).mean()
    return m4 / m2**2 - 3


# Frozen RV class
class rv_frozen(object):

    def __init__(self, dist, *args, **kwds):
        self.args = args
        self.kwds = kwds
        self.dist = dist

    def pdf(self, x):    # raises AttributeError in frozen discrete distribution
        return self.dist.pdf(x, *self.args, **self.kwds)

    def logpdf(self, x):
        return self.dist.logpdf(x, *self.args, **self.kwds)

    def cdf(self, x):
        return self.dist.cdf(x, *self.args, **self.kwds)

    def logcdf(self, x):
        return self.dist.logcdf(x, *self.args, **self.kwds)

    def ppf(self, q):
        return self.dist.ppf(q, *self.args, **self.kwds)

    def isf(self, q):
        return self.dist.isf(q, *self.args, **self.kwds)

    def rvs(self, size=None):
        kwds = self.kwds.copy()
        kwds.update({'size':size})
        return self.dist.rvs(*self.args, **kwds)

    def sf(self, x):
        return self.dist.sf(x, *self.args, **self.kwds)

    def logsf(self, x):
        return self.dist.logsf(x, *self.args, **self.kwds)

    def stats(self, moments='mv'):
        kwds = self.kwds.copy()
        kwds.update({'moments':moments})
        return self.dist.stats(*self.args, **kwds)

    def median(self):
        return self.dist.median(*self.args, **self.kwds)

    def mean(self):
        return self.dist.mean(*self.args, **self.kwds)

    def var(self):
        return self.dist.var(*self.args, **self.kwds)

    def std(self):
        return self.dist.std(*self.args, **self.kwds)

    def moment(self, n):
        return self.dist.moment(n, *self.args, **self.kwds)

    def entropy(self):
        return self.dist.entropy(*self.args, **self.kwds)

    def pmf(self,k):
        return self.dist.pmf(k, *self.args, **self.kwds)

    def logpmf(self,k):
        return self.dist.logpmf(k, *self.args, **self.kwds)

    def interval(self, alpha):
        return self.dist.interval(alpha, *self.args, **self.kwds)


def valarray(shape,value=nan,typecode=None):
    """Return an array of all value.
    """

    out = ones(shape, dtype=bool) * value
    if typecode is not None:
        out = out.astype(typecode)
    if not isinstance(out, ndarray):
        out = asarray(out)
    return out


# This should be rewritten
def argsreduce(cond, *args):
    """Return the sequence of ravel(args[i]) where ravel(condition) is
    True in 1D.

    Examples
    --------
    >>> import numpy as np
    >>> rand = np.random.random_sample
    >>> A = rand((4,5))
    >>> B = 2
    >>> C = rand((1,5))
    >>> cond = np.ones(A.shape)
    >>> [A1,B1,C1] = argsreduce(cond,A,B,C)
    >>> B1.shape
    (20,)
    >>> cond[2,:] = 0
    >>> [A2,B2,C2] = argsreduce(cond,A,B,C)
    >>> B2.shape
    (15,)

    """
    newargs = atleast_1d(*args)
    if not isinstance(newargs, list):
        newargs = [newargs,]
    expand_arr = (cond == cond)
    return [extract(cond, arr1 * expand_arr) for arr1 in newargs]



parse_arg_template = """
def _parse_args(self, %(shape_arg_str)s %(locscale_in)s):
    return (%(shape_arg_str)s), %(locscale_out)s

def _parse_args_rvs(self, %(shape_arg_str)s %(locscale_in)s, size=None):
    return (%(shape_arg_str)s), %(locscale_out)s, size

def _parse_args_stats(self, %(shape_arg_str)s %(locscale_in)s, moments='mv'):
    return (%(shape_arg_str)s), %(locscale_out)s, moments
"""


class rv_generic(object):
    """Class which encapsulates common functionality between rv_discrete
    and rv_continuous.

    """
    def _construct_argparser(self, names_to_inspect, locscale_in, locscale_out):
        """Construct the parser for the shape arguments.

        Generates the argument-parsing functions dynamically.
        Modifies the calling class.
        Is supposed to be called in __init__ of a class for each distribution.

        If self.shapes is a non-empty string, interprets it as a comma-separated
        list of shape parameters.

        Otherwise inspects the call signatures of `names_to_inspect`
        and constructs the argument-parsing functions from these.
        In this case also sets `shapes` and `numargs`.
        """

        if self.shapes:
            # sanitize the user-supplied shapes
            if not isinstance(self.shapes, string_types):
                raise TypeError('shapes must be a string.')

            shapes = self.shapes.replace(',', ' ').split()

            for field in shapes:
                if keyword.iskeyword(field):
                    raise SyntaxError('keywords cannot be used as shapes.')
                if not re.match('^[_a-zA-Z][_a-zA-Z0-9]*$', field):
                    raise SyntaxError('shapes must be valid python identifiers')
        else:
            # find out the call signatures (_pdf, _cdf etc), deduce shape arguments
            shapes_list = []
            for name in names_to_inspect:
                # look for names in instance methods, then global namespace
                # the latter is needed for rv_discrete with explicit `values`
                try:
                    meth = get_method_function(getattr(self, name))
                except:
                    meth = globals()[name]
                shapes_args = inspect.getargspec(meth)
                shapes_list.append(shapes_args.args)

                # *args or **kwargs are not allowed w/automatic shapes
                # (generic methods have 'self, x' only)
                if len(shapes_args.args) > 2:
                    if shapes_args.varargs is not None:
                        raise TypeError('*args are not allowed w/out explicit shapes')
                    if shapes_args.keywords is not None:
                        raise TypeError('**kwds are not allowed w/out explicit shapes')
                    if shapes_args.defaults is not None:
                        raise TypeError('defaults are not allowed for shapes')

            shapes = max(shapes_list, key=lambda x: len(x))
            shapes = shapes[2:]  # remove self, x,

            # make sure the signatures are consistent
            # (generic methods have 'self, x' only)
            for item in shapes_list:
                if len(item) > 2 and item[2:] != shapes:
                    raise TypeError('Shape arguments are inconsistent.')

        # have the arguments, construct the method from template
        shapes_str = ', '.join(shapes) + ', ' if shapes else ''  # NB: not None
        dct = dict(shape_arg_str=shapes_str,
                   locscale_in=locscale_in,
                   locscale_out=locscale_out,
                  )
        ns = {}
        exec_(parse_arg_template % dct, ns)
        # NB: attach to the instance, not class
        for name in ['_parse_args', '_parse_args_stats', '_parse_args_rvs']:
            setattr(self, name,
                instancemethod(ns[name], self, self.__class__)
            )

        self.shapes = ', '.join(shapes) if shapes else None
        if not hasattr(self, 'numargs'):
            # allows more general subclassing with *args
            self.numargs = len(shapes)


    # These are actually called, and should not be overwritten if you
    # want to keep error checking.
    def rvs(self, *args, **kwds):
        """
        Random variates of given type.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).
        scale : array_like, optional
            Scale parameter (default=1).
        size : int or tuple of ints, optional
            Defining number of random variates (default=1).

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of given `size`.

        """
        discrete = kwds.pop('discrete', None)
        args, loc, scale, size = self._parse_args_rvs(*args, **kwds)
        cond = logical_and(self._argcheck(*args), (scale >= 0))
        if not all(cond):
            raise ValueError("Domain error in arguments.")

        # self._size is total size of all output values
        self._size = product(size, axis=0)
        if self._size is not None and self._size > 1:
            size = numpy.array(size, ndmin=1)

        if np.all(scale == 0):
            return loc*ones(size, 'd')

        vals = self._rvs(*args)
        if self._size is not None:
            vals = reshape(vals, size)

        vals = vals * scale + loc

        # Cast to int if discrete
        if discrete:
            if numpy.isscalar(vals):
                vals = int(vals)
            else:
                vals = vals.astype(int)

        return vals

    def median(self, *args, **kwds):
        """
        Median of the distribution.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            Location parameter, Default is 0.
        scale : array_like, optional
            Scale parameter, Default is 1.

        Returns
        -------
        median : float
            The median of the distribution.

        See Also
        --------
        stats.distributions.rv_discrete.ppf
            Inverse of the CDF

        """
        return self.ppf(0.5, *args, **kwds)

    def mean(self, *args, **kwds):
        """
        Mean of the distribution

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        mean : float
            the mean of the distribution
        """
        kwds['moments'] = 'm'
        res = self.stats(*args, **kwds)
        if isinstance(res, ndarray) and res.ndim == 0:
            return res[()]
        return res

    def var(self, *args, **kwds):
        """
        Variance of the distribution

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        var : float
            the variance of the distribution

        """
        kwds['moments'] = 'v'
        res = self.stats(*args, **kwds)
        if isinstance(res, ndarray) and res.ndim == 0:
            return res[()]
        return res

    def std(self, *args, **kwds):
        """
        Standard deviation of the distribution.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        std : float
            standard deviation of the distribution

        """
        kwds['moments'] = 'v'
        res = sqrt(self.stats(*args, **kwds))
        return res

    def interval(self, alpha, *args, **kwds):
        """
        Confidence interval with equal areas around the median.

        Parameters
        ----------
        alpha : array_like of float
            Probability that an rv will be drawn from the returned range.
            Each value should be in the range [0, 1].
        arg1, arg2, ... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            location parameter, Default is 0.
        scale : array_like, optional
            scale parameter, Default is 1.

        Returns
        -------
        a, b : ndarray of float
            end-points of range that contain ``100 * alpha %`` of the rv's possible
            values.

        """
        alpha = asarray(alpha)
        if any((alpha > 1) | (alpha < 0)):
            raise ValueError("alpha must be between 0 and 1 inclusive")
        q1 = (1.0-alpha)/2
        q2 = (1.0+alpha)/2
        a = self.ppf(q1, *args, **kwds)
        b = self.ppf(q2, *args, **kwds)
        return a, b


##  continuous random variables: implement maybe later
##
##  hf  --- Hazard Function (PDF / SF)
##  chf  --- Cumulative hazard function (-log(SF))
##  psf --- Probability sparsity function (reciprocal of the pdf) in
##                units of percent-point-function (as a function of q).
##                Also, the derivative of the percent-point function.

class rv_continuous(rv_generic):
    """
    A generic continuous random variable class meant for subclassing.

    `rv_continuous` is a base class to construct specific distribution classes
    and instances from for continuous random variables. It cannot be used
    directly as a distribution.

    Parameters
    ----------
    momtype : int, optional
        The type of generic moment calculation to use: 0 for pdf, 1 (default)
        for ppf.
    a : float, optional
        Lower bound of the support of the distribution, default is minus
        infinity.
    b : float, optional
        Upper bound of the support of the distribution, default is plus
        infinity.
    xtol : float, optional
        The tolerance for fixed point calculation for generic ppf.
    badvalue : object, optional
        The value in a result arrays that indicates a value that for which
        some argument restriction is violated, default is np.nan.
    name : str, optional
        The name of the instance. This string is used to construct the default
        example for distributions.
    longname : str, optional
        This string is used as part of the first line of the docstring returned
        when a subclass has no docstring of its own. Note: `longname` exists
        for backwards compatibility, do not use for new subclasses.
    shapes : str, optional
        The shape of the distribution. For example ``"m, n"`` for a
        distribution that takes two integers as the two shape arguments for all
        its methods.
    extradoc :  str, optional, deprecated
        This string is used as the last part of the docstring returned when a
        subclass has no docstring of its own. Note: `extradoc` exists for
        backwards compatibility, do not use for new subclasses.

    Methods
    -------
    rvs(<shape(s)>, loc=0, scale=1, size=1)
        random variates

    pdf(x, <shape(s)>, loc=0, scale=1)
        probability density function

    logpdf(x, <shape(s)>, loc=0, scale=1)
        log of the probability density function

    cdf(x, <shape(s)>, loc=0, scale=1)
        cumulative density function

    logcdf(x, <shape(s)>, loc=0, scale=1)
        log of the cumulative density function

    sf(x, <shape(s)>, loc=0, scale=1)
        survival function (1-cdf --- sometimes more accurate)

    logsf(x, <shape(s)>, loc=0, scale=1)
        log of the survival function

    ppf(q, <shape(s)>, loc=0, scale=1)
      percent point function (inverse of cdf --- quantiles)

    isf(q, <shape(s)>, loc=0, scale=1)
        inverse survival function (inverse of sf)

    moment(n, <shape(s)>, loc=0, scale=1)
        non-central n-th moment of the distribution.  May not work for array arguments.

    stats(<shape(s)>, loc=0, scale=1, moments='mv')
        mean('m'), variance('v'), skew('s'), and/or kurtosis('k')

    entropy(<shape(s)>, loc=0, scale=1)
        (differential) entropy of the RV.

    fit(data, <shape(s)>, loc=0, scale=1)
        Parameter estimates for generic data

    expect(func=None, args=(), loc=0, scale=1, lb=None, ub=None,
             conditional=False, **kwds)
        Expected value of a function with respect to the distribution.
        Additional kwd arguments passed to integrate.quad

    median(<shape(s)>, loc=0, scale=1)
        Median of the distribution.

    mean(<shape(s)>, loc=0, scale=1)
        Mean of the distribution.

    std(<shape(s)>, loc=0, scale=1)
        Standard deviation of the distribution.

    var(<shape(s)>, loc=0, scale=1)
        Variance of the distribution.

    interval(alpha, <shape(s)>, loc=0, scale=1)
        Interval that with `alpha` percent probability contains a random
        realization of this distribution.

    __call__(<shape(s)>, loc=0, scale=1)
        Calling a distribution instance creates a frozen RV object with the
        same methods but holding the given shape, location, and scale fixed.
        See Notes section.

    **Parameters for Methods**

    x : array_like
        quantiles
    q : array_like
        lower or upper tail probability
    <shape(s)> : array_like
        shape parameters
    loc : array_like, optional
        location parameter (default=0)
    scale : array_like, optional
        scale parameter (default=1)
    size : int or tuple of ints, optional
        shape of random variates (default computed from input arguments )
    moments : string, optional
        composed of letters ['mvsk'] specifying which moments to compute where
        'm' = mean, 'v' = variance, 's' = (Fisher's) skew and
        'k' = (Fisher's) kurtosis. (default='mv')
    n : int
        order of moment to calculate in method moments

    Notes
    -----

    **Methods that can be overwritten by subclasses**
    ::

      _rvs
      _pdf
      _cdf
      _sf
      _ppf
      _isf
      _stats
      _munp
      _entropy
      _argcheck

    There are additional (internal and private) generic methods that can
    be useful for cross-checking and for debugging, but might work in all
    cases when directly called.

    **Frozen Distribution**

    Alternatively, the object may be called (as a function) to fix the shape,
    location, and scale parameters returning a "frozen" continuous RV object:

    rv = generic(<shape(s)>, loc=0, scale=1)
        frozen RV object with the same methods but holding the given shape,
        location, and scale fixed

    **Subclassing**

    New random variables can be defined by subclassing rv_continuous class
    and re-defining at least the ``_pdf`` or the ``_cdf`` method (normalized
    to location 0 and scale 1) which will be given clean arguments (in between
    a and b) and passing the argument check method.

    If positive argument checking is not correct for your RV
    then you will also need to re-define the ``_argcheck`` method.

    Correct, but potentially slow defaults exist for the remaining
    methods but for speed and/or accuracy you can over-ride::

      _logpdf, _cdf, _logcdf, _ppf, _rvs, _isf, _sf, _logsf

    Rarely would you override ``_isf``, ``_sf`` or ``_logsf``, but you could.

    Statistics are computed using numerical integration by default.
    For speed you can redefine this using ``_stats``:

     - take shape parameters and return mu, mu2, g1, g2
     - If you can't compute one of these, return it as None
     - Can also be defined with a keyword argument ``moments=<str>``,
       where <str> is a string composed of 'm', 'v', 's',
       and/or 'k'.  Only the components appearing in string
       should be computed and returned in the order 'm', 'v',
       's', or 'k'  with missing values returned as None.

    Alternatively, you can override ``_munp``, which takes n and shape
    parameters and returns the nth non-central moment of the distribution.

    A note on ``shapes``: subclasses need not specify them explicitly. In this
    case, the `shapes` will be automatically deduced from the signatures of the
    overridden methods.
    If, for some reason, you prefer to avoid relying on introspection, you can
    specify ``shapes`` explicitly as an argument to the instance constructor.

    Examples
    --------
    To create a new Gaussian distribution, we would do the following::

        class gaussian_gen(rv_continuous):
            "Gaussian distribution"
            def _pdf(self, x):
                ...
            ...

    """

    def __init__(self, momtype=1, a=None, b=None, xtol=1e-14,
                 badvalue=None, name=None, longname=None,
                 shapes=None, extradoc=None):

        rv_generic.__init__(self)

        if badvalue is None:
            badvalue = nan
        if name is None:
            name = 'Distribution'
        self.badvalue = badvalue
        self.name = name
        self.a = a
        self.b = b
        if a is None:
            self.a = -inf
        if b is None:
            self.b = inf
        self.xtol = xtol
        self._size = 1
        self.m = 0.0
        self.moment_type = momtype

        self.expandarr = 1

        self.shapes = shapes
        self._construct_argparser(names_to_inspect=['_pdf', '_cdf'],
                                  locscale_in='loc=0, scale=1',
                                  locscale_out='loc, scale')

        # nin correction
        self.vecfunc = vectorize(self._ppf_single_call, otypes='d')
        self.vecfunc.nin = self.numargs + 1
        self.vecentropy = vectorize(self._entropy, otypes='d')
        self.vecentropy.nin = self.numargs + 1
        self.veccdf = vectorize(self._cdf_single_call, otypes='d')
        self.veccdf.nin = self.numargs + 1

        self.extradoc = extradoc
        if momtype == 0:
            self.generic_moment = vectorize(self._mom0_sc, otypes='d')
        else:
            self.generic_moment = vectorize(self._mom1_sc, otypes='d')
        self.generic_moment.nin = self.numargs+1  # Because of the *args argument
        # of _mom0_sc, vectorize cannot count the number of arguments correctly.

        if longname is None:
            if name[0] in ['aeiouAEIOU']:
                hstr = "An "
            else:
                hstr = "A "
            longname = hstr + name

        if sys.flags.optimize < 2:
            # Skip adding docstrings if interpreter is run with -OO
            if self.__doc__ is None:
                self._construct_default_doc(longname=longname, extradoc=extradoc)
            else:
                self._construct_doc()

        ## This only works for old-style classes...
        # self.__class__.__doc__ = self.__doc__


    def _construct_default_doc(self, longname=None, extradoc=None):
        """Construct instance docstring from the default template."""
        if longname is None:
            longname = 'A'
        if extradoc is None:
            extradoc = ''
        if extradoc.startswith('\n\n'):
            extradoc = extradoc[2:]
        self.__doc__ = ''.join(['%s continuous random variable.' % longname,
                                '\n\n%(before_notes)s\n', docheaders['notes'],
                                extradoc, '\n%(example)s'])
        self._construct_doc()

    def _construct_doc(self):
        """Construct the instance docstring with string substitutions."""
        tempdict = docdict.copy()
        tempdict['name'] = self.name or 'distname'
        tempdict['shapes'] = self.shapes or ''

        if self.shapes is None:
            # remove shapes from call parameters if there are none
            for item in ['callparams', 'default', 'before_notes']:
                tempdict[item] = tempdict[item].replace(
                        "\n%(shapes)s : array_like\n    shape parameters", "")
        for i in range(2):
            if self.shapes is None:
                # necessary because we use %(shapes)s in two forms (w w/o ", ")
                self.__doc__ = self.__doc__.replace("%(shapes)s, ", "")
            self.__doc__ = doccer.docformat(self.__doc__, tempdict)

    def _ppf_to_solve(self, x, q,*args):
        return self.cdf(*(x, )+args)-q

    def _ppf_single_call(self, q, *args):
        left = right = None
        if self.a > -np.inf:
            left = self.a
        if self.b < np.inf:
            right = self.b

        factor = 10.
        if not left:  # i.e. self.a = -inf
            left = -1.*factor
            while self._ppf_to_solve(left, q,*args) > 0.:
                right = left
                left *= factor
            # left is now such that cdf(left) < q
        if not right:  # i.e. self.b = inf
            right = factor
            while self._ppf_to_solve(right, q,*args) < 0.:
                left = right
                right *= factor
            # right is now such that cdf(right) > q

        return optimize.brentq(self._ppf_to_solve,
                               left, right, args=(q,)+args, xtol=self.xtol)

    # moment from definition
    def _mom_integ0(self, x,m,*args):
        return x**m * self.pdf(x,*args)

    def _mom0_sc(self, m,*args):
        return integrate.quad(self._mom_integ0, self.a,
                                    self.b, args=(m,)+args)[0]

    # moment calculated using ppf
    def _mom_integ1(self, q,m,*args):
        return (self.ppf(q,*args))**m

    def _mom1_sc(self, m,*args):
        return integrate.quad(self._mom_integ1, 0, 1,args=(m,)+args)[0]

    ## These are the methods you must define (standard form functions)
    def _argcheck(self, *args):
        # Default check for correct values on args and keywords.
        # Returns condition array of 1's where arguments are correct and
        #  0's where they are not.
        cond = 1
        for arg in args:
            cond = logical_and(cond,(asarray(arg) > 0))
        return cond

    def _pdf(self,x,*args):
        return derivative(self._cdf,x,dx=1e-5,args=args,order=5)

    ## Could also define any of these
    def _logpdf(self, x, *args):
        return log(self._pdf(x, *args))

    ##(return 1-d using self._size to get number)
    def _rvs(self, *args):
        ## Use basic inverse cdf algorithm for RV generation as default.
        U = mtrand.sample(self._size)
        Y = self._ppf(U,*args)
        return Y

    def _cdf_single_call(self, x, *args):
        return integrate.quad(self._pdf, self.a, x, args=args)[0]

    def _cdf(self, x, *args):
        return self.veccdf(x,*args)

    def _logcdf(self, x, *args):
        return log(self._cdf(x, *args))

    def _sf(self, x, *args):
        return 1.0-self._cdf(x,*args)

    def _logsf(self, x, *args):
        return log(self._sf(x, *args))

    def _ppf(self, q, *args):
        return self.vecfunc(q,*args)

    def _isf(self, q, *args):
        return self._ppf(1.0-q,*args)  # use correct _ppf for subclasses

    # The actual calculation functions (no basic checking need be done)
    # If these are defined, the others won't be looked at.
    # Otherwise, the other set can be defined.
    def _stats(self,*args, **kwds):
        return None, None, None, None

    #  Central moments
    def _munp(self,n,*args):
        return self.generic_moment(n,*args)

    def pdf(self,x,*args,**kwds):
        """
        Probability density function at x of the given RV.

        Parameters
        ----------
        x : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at x

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        x,loc,scale = map(asarray,(x,loc,scale))
        args = tuple(map(asarray,args))
        x = asarray((x-loc)*1.0/scale)
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x >= self.a) & (x <= self.b)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        putmask(output,(1-cond0)+np.isnan(x),self.badvalue)
        if any(cond):
            goodargs = argsreduce(cond, *((x,)+args+(scale,)))
            scale, goodargs = goodargs[-1], goodargs[:-1]
            place(output,cond,self._pdf(*goodargs) / scale)
        if output.ndim == 0:
            return output[()]
        return output

    def logpdf(self, x, *args, **kwds):
        """
        Log of the probability density function at x of the given RV.

        This uses a more numerically accurate calculation if available.

        Parameters
        ----------
        x : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        logpdf : array_like
            Log of the probability density function evaluated at x

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        x,loc,scale = map(asarray,(x,loc,scale))
        args = tuple(map(asarray,args))
        x = asarray((x-loc)*1.0/scale)
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x >= self.a) & (x <= self.b)
        cond = cond0 & cond1
        output = empty(shape(cond),'d')
        output.fill(NINF)
        putmask(output,(1-cond0)+np.isnan(x),self.badvalue)
        if any(cond):
            goodargs = argsreduce(cond, *((x,)+args+(scale,)))
            scale, goodargs = goodargs[-1], goodargs[:-1]
            place(output,cond,self._logpdf(*goodargs) - log(scale))
        if output.ndim == 0:
            return output[()]
        return output

    def cdf(self,x,*args,**kwds):
        """
        Cumulative distribution function of the given RV.

        Parameters
        ----------
        x : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        cdf : ndarray
            Cumulative distribution function evaluated at `x`

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        x,loc,scale = map(asarray,(x,loc,scale))
        args = tuple(map(asarray,args))
        x = (x-loc)*1.0/scale
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond2 = (x >= self.b) & cond0
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        place(output,(1-cond0)+np.isnan(x),self.badvalue)
        place(output,cond2,1.0)
        if any(cond):  # call only if at least 1 entry
            goodargs = argsreduce(cond, *((x,)+args))
            place(output,cond,self._cdf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logcdf(self,x,*args,**kwds):
        """
        Log of the cumulative distribution function at x of the given RV.

        Parameters
        ----------
        x : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        logcdf : array_like
            Log of the cumulative distribution function evaluated at x

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        x,loc,scale = map(asarray,(x,loc,scale))
        args = tuple(map(asarray,args))
        x = (x-loc)*1.0/scale
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond2 = (x >= self.b) & cond0
        cond = cond0 & cond1
        output = empty(shape(cond),'d')
        output.fill(NINF)
        place(output,(1-cond0)*(cond1 == cond1)+np.isnan(x),self.badvalue)
        place(output,cond2,0.0)
        if any(cond):  # call only if at least 1 entry
            goodargs = argsreduce(cond, *((x,)+args))
            place(output,cond,self._logcdf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def sf(self,x,*args,**kwds):
        """
        Survival function (1-cdf) at x of the given RV.

        Parameters
        ----------
        x : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        sf : array_like
            Survival function evaluated at x

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        x,loc,scale = map(asarray,(x,loc,scale))
        args = tuple(map(asarray,args))
        x = (x-loc)*1.0/scale
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond2 = cond0 & (x <= self.a)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        place(output,(1-cond0)+np.isnan(x),self.badvalue)
        place(output,cond2,1.0)
        if any(cond):
            goodargs = argsreduce(cond, *((x,)+args))
            place(output,cond,self._sf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logsf(self,x,*args,**kwds):
        """
        Log of the survival function of the given RV.

        Returns the log of the "survival function," defined as (1 - `cdf`),
        evaluated at `x`.

        Parameters
        ----------
        x : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        logsf : ndarray
            Log of the survival function evaluated at `x`.

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        x,loc,scale = map(asarray,(x,loc,scale))
        args = tuple(map(asarray,args))
        x = (x-loc)*1.0/scale
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond2 = cond0 & (x <= self.a)
        cond = cond0 & cond1
        output = empty(shape(cond),'d')
        output.fill(NINF)
        place(output,(1-cond0)+np.isnan(x),self.badvalue)
        place(output,cond2,0.0)
        if any(cond):
            goodargs = argsreduce(cond, *((x,)+args))
            place(output,cond,self._logsf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def ppf(self,q,*args,**kwds):
        """
        Percent point function (inverse of cdf) at q of the given RV.

        Parameters
        ----------
        q : array_like
            lower tail probability
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        x : array_like
            quantile corresponding to the lower tail probability q.

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        q, loc, scale = map(asarray,(q, loc, scale))
        args = tuple(map(asarray, args))
        cond0 = self._argcheck(*args) & (scale > 0) & (loc == loc)
        cond1 = (0 < q) & (q < 1)
        cond2 = cond0 & (q == 0)
        cond3 = cond0 & (q == 1)
        cond = cond0 & cond1
        output = valarray(shape(cond), value=self.badvalue)

        lower_bound = self.a * scale + loc
        upper_bound = self.b * scale + loc
        place(output, cond2, argsreduce(cond2, lower_bound)[0])
        place(output, cond3, argsreduce(cond3, upper_bound)[0])

        if any(cond):  # call only if at least 1 entry
            goodargs = argsreduce(cond, *((q,)+args+(scale,loc)))
            scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]
            place(output, cond, self._ppf(*goodargs) * scale + loc)
        if output.ndim == 0:
            return output[()]
        return output

    def isf(self, q, *args, **kwds):
        """
        Inverse survival function at q of the given RV.

        Parameters
        ----------
        q : array_like
            upper tail probability
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)

        Returns
        -------
        x : ndarray or scalar
            Quantile corresponding to the upper tail probability q.

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        q, loc, scale = map(asarray, (q, loc, scale))
        args = tuple(map(asarray, args))
        cond0 = self._argcheck(*args) & (scale > 0) & (loc == loc)
        cond1 = (0 < q) & (q < 1)
        cond2 = cond0 & (q == 1)
        cond3 = cond0 & (q == 0)
        cond = cond0 & cond1
        output = valarray(shape(cond), value=self.badvalue)

        lower_bound = self.a * scale + loc
        upper_bound = self.b * scale + loc
        place(output, cond2, argsreduce(cond2, lower_bound)[0])
        place(output, cond3, argsreduce(cond3, upper_bound)[0])

        if any(cond):
            goodargs = argsreduce(cond, *((q,)+args+(scale,loc)))
            scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]
            place(output, cond, self._isf(*goodargs) * scale + loc)
        if output.ndim == 0:
            return output[()]
        return output

    def stats(self,*args,**kwds):
        """
        Some statistics of the given RV

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)
        moments : str, optional
            composed of letters ['mvsk'] defining which moments to compute:
            'm' = mean,
            'v' = variance,
            's' = (Fisher's) skew,
            'k' = (Fisher's) kurtosis.
            (default='mv')

        Returns
        -------
        stats : sequence
            of requested moments.

        """
        args, loc, scale, moments = self._parse_args_stats(*args, **kwds)
        loc, scale = map(asarray, (loc, scale))
        args = tuple(map(asarray, args))
        cond = self._argcheck(*args) & (scale > 0) & (loc == loc)

        signature = inspect.getargspec(get_method_function(self._stats))
        if (signature[2] is not None) or ('moments' in signature[0]):
            mu, mu2, g1, g2 = self._stats(*args,**{'moments':moments})
        else:
            mu, mu2, g1, g2 = self._stats(*args)
        if g1 is None:
            mu3 = None
        else:
            mu3 = g1*np.power(mu2,1.5)  # (mu2**1.5) breaks down for nan and inf
        default = valarray(shape(cond), self.badvalue)
        output = []

        # Use only entries that are valid in calculation
        if any(cond):
            goodargs = argsreduce(cond, *(args+(scale,loc)))
            scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]
            if 'm' in moments:
                if mu is None:
                    mu = self._munp(1.0,*goodargs)
                out0 = default.copy()
                place(out0,cond,mu*scale+loc)
                output.append(out0)

            if 'v' in moments:
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    if mu is None:
                        mu = self._munp(1.0,*goodargs)
                    mu2 = mu2p - mu*mu
                    if np.isinf(mu):
                        #if mean is inf then var is also inf
                        mu2 = np.inf
                out0 = default.copy()
                place(out0,cond,mu2*scale*scale)
                output.append(out0)

            if 's' in moments:
                if g1 is None:
                    mu3p = self._munp(3.0,*goodargs)
                    if mu is None:
                        mu = self._munp(1.0,*goodargs)
                    if mu2 is None:
                        mu2p = self._munp(2.0,*goodargs)
                        mu2 = mu2p - mu*mu
                    mu3 = mu3p - 3*mu*mu2 - mu**3
                    g1 = mu3 / np.power(mu2, 1.5)
                out0 = default.copy()
                place(out0,cond,g1)
                output.append(out0)

            if 'k' in moments:
                if g2 is None:
                    mu4p = self._munp(4.0,*goodargs)
                    if mu is None:
                        mu = self._munp(1.0,*goodargs)
                    if mu2 is None:
                        mu2p = self._munp(2.0,*goodargs)
                        mu2 = mu2p - mu*mu
                    if mu3 is None:
                        mu3p = self._munp(3.0,*goodargs)
                        mu3 = mu3p - 3*mu*mu2 - mu**3
                    mu4 = mu4p - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4
                    g2 = mu4 / mu2**2.0 - 3.0
                out0 = default.copy()
                place(out0,cond,g2)
                output.append(out0)
        else:  # no valid args
            output = []
            for _ in moments:
                out0 = default.copy()
                output.append(out0)

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)

    def moment(self, n, *args, **kwds):
        """
        n'th order non-central moment of distribution.

        Parameters
        ----------
        n : int, n>=1
            Order of moment.
        arg1, arg2, arg3,... : float
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        kwds : keyword arguments, optional
            These can include "loc" and "scale", as well as other keyword
            arguments relevant for a given distribution.

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        if not (self._argcheck(*args) and (scale > 0)):
            return nan
        if (floor(n) != n):
            raise ValueError("Moment must be an integer.")
        if (n < 0):
            raise ValueError("Moment must be positive.")
        mu, mu2, g1, g2 = None, None, None, None
        if (n > 0) and (n < 5):
            signature = inspect.getargspec(get_method_function(self._stats))
            if (signature[2] is not None) or ('moments' in signature[0]):
                mdict = {'moments':{1:'m',2:'v',3:'vs',4:'vk'}[n]}
            else:
                mdict = {}
            mu, mu2, g1, g2 = self._stats(*args,**mdict)
        val = _moment_from_stats(n, mu, mu2, g1, g2, self._munp, args)

        # Convert to transformed  X = L + S*Y
        # so E[X^n] = E[(L+S*Y)^n] = L^n sum(comb(n,k)*(S/L)^k E[Y^k],k=0...n)
        if loc == 0:
            return scale**n * val
        else:
            result = 0
            fac = float(scale) / float(loc)
            for k in range(n):
                valk = _moment_from_stats(k, mu, mu2, g1, g2, self._munp, args)
                result += comb(n,k,exact=True)*(fac**k) * valk
            result += fac**n * val
            return result * loc**n

    def _nnlf(self, x, *args):
        return -sum(self._logpdf(x, *args),axis=0)

    def nnlf(self, theta, x):
        '''Return negative loglikelihood function

        Notes
        -----
        This is ``-sum(log pdf(x, theta), axis=0)`` where theta are the
        parameters (including loc and scale).
        '''
        try:
            loc = theta[-2]
            scale = theta[-1]
            args = tuple(theta[:-2])
        except IndexError:
            raise ValueError("Not enough input arguments.")
        if not self._argcheck(*args) or scale <= 0:
            return inf
        x = asarray((x-loc) / scale)
        cond0 = (x <= self.a) | (self.b <= x)
        if (any(cond0)):
            return inf
        else:
            N = len(x)
            return self._nnlf(x, *args) + N * log(scale)

    def _penalized_nnlf(self, theta, x):
        ''' Return negative loglikelihood function,
        i.e., - sum (log pdf(x, theta),axis=0)
           where theta are the parameters (including loc and scale)
        '''
        try:
            loc = theta[-2]
            scale = theta[-1]
            args = tuple(theta[:-2])
        except IndexError:
            raise ValueError("Not enough input arguments.")
        if not self._argcheck(*args) or scale <= 0:
            return inf
        x = asarray((x-loc) / scale)

        loginf = log(floatinfo.machar.xmax)

        if np.isneginf(self.a).all() and np.isinf(self.b).all():
            Nbad = 0
        else:
            cond0 = (x <= self.a) | (self.b <= x)
            Nbad = sum(cond0)
            if Nbad > 0:
                x = argsreduce(~cond0, x)[0]

        N = len(x)
        return self._nnlf(x, *args) + N*log(scale) + Nbad * 100.0 * loginf

    # return starting point for fit (shape arguments + loc + scale)
    def _fitstart(self, data, args=None):
        if args is None:
            args = (1.0,)*self.numargs
        return args + self.fit_loc_scale(data, *args)

    # Return the (possibly reduced) function to optimize in order to find MLE
    #  estimates for the .fit method
    def _reduce_func(self, args, kwds):
        args = list(args)
        Nargs = len(args)
        fixedn = []
        index = list(range(Nargs))
        names = ['f%d' % n for n in range(Nargs - 2)] + ['floc', 'fscale']
        x0 = []
        for n, key in zip(index, names):
            if key in kwds:
                fixedn.append(n)
                args[n] = kwds[key]
            else:
                x0.append(args[n])

        if len(fixedn) == 0:
            func = self._penalized_nnlf
            restore = None
        else:
            if len(fixedn) == len(index):
                raise ValueError("All parameters fixed. There is nothing to optimize.")

            def restore(args, theta):
                # Replace with theta for all numbers not in fixedn
                # This allows the non-fixed values to vary, but
                #  we still call self.nnlf with all parameters.
                i = 0
                for n in range(Nargs):
                    if n not in fixedn:
                        args[n] = theta[i]
                        i += 1
                return args

            def func(theta, x):
                newtheta = restore(args[:], theta)
                return self._penalized_nnlf(newtheta, x)

        return x0, func, restore, args

    def fit(self, data, *args, **kwds):
        """
        Return MLEs for shape, location, and scale parameters from data.

        MLE stands for Maximum Likelihood Estimate.  Starting estimates for
        the fit are given by input arguments; for any arguments not provided
        with starting estimates, ``self._fitstart(data)`` is called to generate
        such.

        One can hold some parameters fixed to specific values by passing in
        keyword arguments ``f0``, ``f1``, ..., ``fn`` (for shape parameters)
        and ``floc`` and ``fscale`` (for location and scale parameters,
        respectively).

        Parameters
        ----------
        data : array_like
            Data to use in calculating the MLEs.
        args : floats, optional
            Starting value(s) for any shape-characterizing arguments (those not
            provided will be determined by a call to ``_fitstart(data)``).
            No default value.
        kwds : floats, optional
            Starting values for the location and scale parameters; no default.
            Special keyword arguments are recognized as holding certain
            parameters fixed:

            f0...fn : hold respective shape parameters fixed.

            floc : hold location parameter fixed to specified value.

            fscale : hold scale parameter fixed to specified value.

            optimizer : The optimizer to use.  The optimizer must take func,
                        and starting position as the first two arguments,
                        plus args (for extra arguments to pass to the
                        function to be optimized) and disp=0 to suppress
                        output as keyword arguments.

        Returns
        -------
        shape, loc, scale : tuple of floats
            MLEs for any shape statistics, followed by those for location and
            scale.

        Notes
        -----
        This fit is computed by maximizing a log-likelihood function, with
        penalty applied for samples outside of range of the distribution. The
        returned answer is not guaranteed to be the globally optimal MLE, it
        may only be locally optimal, or the optimization may fail altogether.
        """
        Narg = len(args)
        if Narg > self.numargs:
            raise TypeError("Too many input arguments.")

        start = [None]*2
        if (Narg < self.numargs) or not ('loc' in kwds and
                                         'scale' in kwds):
            start = self._fitstart(data)  # get distribution specific starting locations
            args += start[Narg:-2]
        loc = kwds.get('loc', start[-2])
        scale = kwds.get('scale', start[-1])
        args += (loc, scale)
        x0, func, restore, args = self._reduce_func(args, kwds)

        optimizer = kwds.get('optimizer', optimize.fmin)
        # convert string to function in scipy.optimize
        if not callable(optimizer) and isinstance(optimizer, string_types):
            if not optimizer.startswith('fmin_'):
                optimizer = "fmin_"+optimizer
            if optimizer == 'fmin_':
                optimizer = 'fmin'
            try:
                optimizer = getattr(optimize, optimizer)
            except AttributeError:
                raise ValueError("%s is not a valid optimizer" % optimizer)
        vals = optimizer(func,x0,args=(ravel(data),),disp=0)
        if restore is not None:
            vals = restore(args, vals)
        vals = tuple(vals)
        return vals

    def fit_loc_scale(self, data, *args):
        """
        Estimate loc and scale parameters from data using 1st and 2nd moments.

        Parameters
        ----------
        data : array_like
            Data to fit.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).

        Returns
        -------
        Lhat : float
            Estimated location parameter for the data.
        Shat : float
            Estimated scale parameter for the data.

        """
        mu, mu2 = self.stats(*args,**{'moments':'mv'})
        tmp = asarray(data)
        muhat = tmp.mean()
        mu2hat = tmp.var()
        Shat = sqrt(mu2hat / mu2)
        Lhat = muhat - Shat*mu
        if not np.isfinite(Lhat):
            Lhat = 0
        if not (np.isfinite(Shat) and (0 < Shat)):
            Shat = 1
        return Lhat, Shat

    @np.deprecate
    def est_loc_scale(self, data, *args):
        """This function is deprecated, use self.fit_loc_scale(data) instead."""
        return self.fit_loc_scale(data, *args)

    def freeze(self,*args,**kwds):
        """Freeze the distribution for the given arguments.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution.  Should include all
            the non-optional arguments, may include ``loc`` and ``scale``.

        Returns
        -------
        rv_frozen : rv_frozen instance
            The frozen distribution.

        """
        return rv_frozen(self,*args,**kwds)

    def __call__(self, *args, **kwds):
        return self.freeze(*args, **kwds)

    def _entropy(self, *args):
        def integ(x):
            val = self._pdf(x, *args)
            return special.xlogy(val, val)

        entr = -integrate.quad(integ,self.a,self.b)[0]
        if not np.isnan(entr):
            return entr
        else:  # try with different limits if integration problems
            low,upp = self.ppf([0.001,0.999],*args)
            if np.isinf(self.b):
                upper = upp
            else:
                upper = self.b
            if np.isinf(self.a):
                lower = low
            else:
                lower = self.a
            return -integrate.quad(integ,lower,upper)[0]

    def entropy(self, *args, **kwds):
        """
        Differential entropy of the RV.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).
        scale : array_like, optional
            Scale parameter (default=1).

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        args = tuple(map(asarray,args))
        cond0 = self._argcheck(*args) & (scale > 0) & (loc == loc)
        output = zeros(shape(cond0),'d')
        place(output,(1-cond0),self.badvalue)
        goodargs = argsreduce(cond0, *args)
        # np.vectorize doesn't work when numargs == 0 in numpy 1.5.1
        if self.numargs == 0:
            place(output,cond0,self._entropy()+log(scale))
        else:
            place(output,cond0,self.vecentropy(*goodargs)+log(scale))

        return output

    def expect(self, func=None, args=(), loc=0, scale=1, lb=None, ub=None,
               conditional=False, **kwds):
        """Calculate expected value of a function with respect to the distribution

        The expected value of a function ``f(x)`` with respect to a
        distribution ``dist`` is defined as::

                    ubound
            E[x] = Integral(f(x) * dist.pdf(x))
                    lbound

        Parameters
        ----------
        func : callable, optional
            Function for which integral is calculated. Takes only one argument.
            The default is the identity mapping f(x) = x.
        args : tuple, optional
            Argument (parameters) of the distribution.
        lb, ub : scalar, optional
            Lower and upper bound for integration. default is set to the support
            of the distribution.
        conditional : bool, optional
            If True, the integral is corrected by the conditional probability
            of the integration interval.  The return value is the expectation
            of the function, conditional on being in the given interval.
            Default is False.

        Additional keyword arguments are passed to the integration routine.

        Returns
        -------
        expect : float
            The calculated expected value.

        Notes
        -----
        The integration behavior of this function is inherited from
        `integrate.quad`.

        """
        lockwds = {'loc': loc,
                   'scale':scale}
        self._argcheck(*args)
        if func is None:
            def fun(x, *args):
                return x * self.pdf(x, *args, **lockwds)
        else:
            def fun(x, *args):
                return func(x) * self.pdf(x, *args, **lockwds)
        if lb is None:
            lb = loc + self.a * scale
        if ub is None:
            ub = loc + self.b * scale
        if conditional:
            invfac = (self.sf(lb, *args, **lockwds)
                      - self.sf(ub, *args, **lockwds))
        else:
            invfac = 1.0
        kwds['args'] = args
        return integrate.quad(fun, lb, ub, **kwds)[0] / invfac


_EULER = 0.577215664901532860606512090082402431042  # -special.psi(1)
_ZETA3 = 1.202056903159594285399738161511449990765  # special.zeta(3,1)  Apery's constant


## Kolmogorov-Smirnov one-sided and two-sided test statistics
class ksone_gen(rv_continuous):
    """General Kolmogorov-Smirnov one-sided test.

    %(default)s

    """
    def _cdf(self, x, n):
        return 1.0 - special.smirnov(n, x)

    def _ppf(self, q, n):
        return special.smirnovi(n, 1.0 - q)
ksone = ksone_gen(a=0.0, name='ksone')


class kstwobign_gen(rv_continuous):
    """Kolmogorov-Smirnov two-sided test for large N.

    %(default)s

    """
    def _cdf(self,x):
        return 1.0-special.kolmogorov(x)

    def _sf(self,x):
        return special.kolmogorov(x)

    def _ppf(self,q):
        return special.kolmogi(1.0-q)
kstwobign = kstwobign_gen(a=0.0, name='kstwobign')


## Normal distribution

# loc = mu, scale = std
# Keep these implementations out of the class definition so they can be reused
# by other distributions.
_norm_pdf_C = math.sqrt(2*pi)
_norm_pdf_logC = math.log(_norm_pdf_C)


def _norm_pdf(x):
    return exp(-x**2/2.0) / _norm_pdf_C


def _norm_logpdf(x):
    return -x**2 / 2.0 - _norm_pdf_logC


def _norm_cdf(x):
    return special.ndtr(x)


def _norm_logcdf(x):
    return special.log_ndtr(x)


def _norm_ppf(q):
    return special.ndtri(q)


def _norm_sf(x):
    return special.ndtr(-x)


def _norm_logsf(x):
    return special.log_ndtr(-x)


def _norm_isf(q):
    return -special.ndtri(q)


class norm_gen(rv_continuous):
    """A normal continuous random variable.

    The location (loc) keyword specifies the mean.
    The scale (scale) keyword specifies the standard deviation.

    %(before_notes)s

    Notes
    -----
    The probability density function for `norm` is::

        norm.pdf(x) = exp(-x**2/2)/sqrt(2*pi)

    %(example)s

    """
    def _rvs(self):
        return mtrand.standard_normal(self._size)

    def _pdf(self,x):
        return _norm_pdf(x)

    def _logpdf(self, x):
        return _norm_logpdf(x)

    def _cdf(self,x):
        return _norm_cdf(x)

    def _logcdf(self, x):
        return _norm_logcdf(x)

    def _sf(self, x):
        return _norm_sf(x)

    def _logsf(self, x):
        return _norm_logsf(x)

    def _ppf(self,q):
        return _norm_ppf(q)

    def _isf(self,q):
        return _norm_isf(q)

    def _stats(self):
        return 0.0, 1.0, 0.0, 0.0

    def _entropy(self):
        return 0.5*(log(2*pi)+1)

    @inherit_docstring_from(rv_continuous)
    def fit(self, data, **kwds):
        """%(super)s
        This function (norm_gen.fit) uses explicit formulas for the maximum
        likelihood estimation of the parameters, so the `optimizer` argument
        is ignored.
        """
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            # Without this check, this function would just return the
            # parameters that were given.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        if floc is None:
            loc = data.mean()
        else:
            loc = floc

        if fscale is None:
            scale = np.sqrt(((data - loc)**2).mean())
        else:
            scale = fscale

        return loc, scale

norm = norm_gen(name='norm')


class alpha_gen(rv_continuous):
    """An alpha continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `alpha` is::

        alpha.pdf(x,a) = 1/(x**2*Phi(a)*sqrt(2*pi)) * exp(-1/2 * (a-1/x)**2),

    where ``Phi(alpha)`` is the normal CDF, ``x > 0``, and ``a > 0``.

    %(example)s

    """
    def _pdf(self, x, a):
        return 1.0/(x**2)/special.ndtr(a)*_norm_pdf(a-1.0/x)

    def _logpdf(self, x, a):
        return -2*log(x) + _norm_logpdf(a-1.0/x) - log(special.ndtr(a))

    def _cdf(self, x, a):
        return special.ndtr(a-1.0/x) / special.ndtr(a)

    def _ppf(self, q, a):
        return 1.0/asarray(a-special.ndtri(q*special.ndtr(a)))

    def _stats(self, a):
        return [inf]*2 + [nan]*2
alpha = alpha_gen(a=0.0, name='alpha')


class anglit_gen(rv_continuous):
    """An anglit continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `anglit` is::

        anglit.pdf(x) = sin(2*x + pi/2) = cos(2*x),

    for ``-pi/4 <= x <= pi/4``.

    %(example)s

    """
    def _pdf(self, x):
        return cos(2*x)

    def _cdf(self, x):
        return sin(x+pi/4)**2.0

    def _ppf(self, q):
        return (arcsin(sqrt(q))-pi/4)

    def _stats(self):
        return 0.0, pi*pi/16-0.5, 0.0, -2*(pi**4 - 96)/(pi*pi-8)**2

    def _entropy(self):
        return 1-log(2)
anglit = anglit_gen(a=-pi/4, b=pi/4, name='anglit')


class arcsine_gen(rv_continuous):
    """An arcsine continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `arcsine` is::

        arcsine.pdf(x) = 1/(pi*sqrt(x*(1-x)))
        for 0 < x < 1.

    %(example)s

    """
    def _pdf(self, x):
        return 1.0/pi/sqrt(x*(1-x))

    def _cdf(self, x):
        return 2.0/pi*arcsin(sqrt(x))

    def _ppf(self, q):
        return sin(pi/2.0*q)**2.0

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
    # This exception is raised by, for example, beta_gen.fit when both floc
    # and fscale  are fixed and there are values in the data not in the open
    # interval (floc, floc+fscale).
    def __init__(self, distr, lower, upper):
        self.args = ("Invalid values in `data`.  Maximum likelihood "
                      "estimation with {distr!r} requires that {lower!r} < x "
                      "< {upper!r} for each x in `data`.".format(distr=distr,
                      lower=lower, upper=upper),)


class FitSolverError(RuntimeError):
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
    psiab = special.psi(a + b)
    func = s1 - n * (-psiab + special.psi(a))
    return func


def _beta_mle_ab(theta, n, s1, s2):
    # Zeros of this function are critical points of
    # the maximum likelihood function.  Solving this system
    # for theta (which contains a and b) gives the MLE for a and b
    # given `n`, `s1` and `s2`.  `s1` is the sum of the logs of the data,
    # and `s2` is the sum of the logs of 1 - data.  `n` is the number
    # of data points.
    a, b = theta
    psiab = special.psi(a + b)
    func = [s1 - n * (-psiab + special.psi(a)),
            s2 - n * (-psiab + special.psi(b))]
    return func


class beta_gen(rv_continuous):
    """A beta continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `beta` is::

        beta.pdf(x, a, b) = gamma(a+b)/(gamma(a)*gamma(b)) * x**(a-1) *
        (1-x)**(b-1),

    for ``0 < x < 1``, ``a > 0``, ``b > 0``.

    %(example)s

    """
    def _rvs(self, a, b):
        return mtrand.beta(a,b,self._size)

    def _pdf(self, x, a, b):
        return np.exp(self._logpdf(x, a, b))

    def _logpdf(self, x, a, b):
        lPx = special.xlog1py(b-1.0, -x) + special.xlogy(a-1.0, x)
        lPx -= special.betaln(a,b)
        return lPx

    def _cdf(self, x, a, b):
        return special.btdtr(a,b,x)

    def _ppf(self, q, a, b):
        return special.btdtri(a,b,q)

    def _stats(self, a, b):
        mn = a*1.0 / (a + b)
        var = (a*b*1.0)/(a+b+1.0)/(a+b)**2.0
        g1 = 2.0*(b-a)*sqrt((1.0+a+b)/(a*b)) / (2+a+b)
        g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
        g2 /= a*b*(a+b+2)*(a+b+3)
        return mn, var, g1, g2

    def _fitstart(self, data):
        g1 = _skew(data)
        g2 = _kurtosis(data)

        def func(x):
            a, b = x
            sk = 2*(b-a)*sqrt(a + b + 1) / (a + b + 2) / sqrt(a*b)
            ku = a**3 - a**2*(2*b-1) + b**2*(b+1) - 2*a*b*(b+2)
            ku /= a*b*(a+b+2)*(a+b+3)
            ku *= 6
            return [sk-g1, ku-g2]
        a, b = optimize.fsolve(func, (1.0, 1.0))
        return super(beta_gen, self)._fitstart(data, args=(a,b))

    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        """%(super)s
        In the special case where both `floc` and `fscale` are given, a
        `ValueError` is raised if any value `x` in `data` does not satisfy
        `floc < x < floc + fscale`.
        """
        # Override rv_continuous.fit, so we can more efficiently handle the
        # case where floc and fscale are given.

        f0 = kwds.get('f0', None)
        f1 = kwds.get('f1', None)
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is None or fscale is None:
            # do general fit
            return super(beta_gen, self).fit(data, *args, **kwds)

        if f0 is not None and f1 is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        # Special case: loc and scale are constrained, so we are fitting
        # just the shape parameters.  This can be done much more efficiently
        # than the method used in `rv_continuous.fit`.  (See the subsection
        # "Two unknown parameters" in the section "Maximum likelihood" of
        # the Wikipedia article on the Beta distribution for the formulas.)

        # Normalize the data to the interval [0,1].
        data = (ravel(data) - floc) / fscale
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
            theta, info, ier, mesg = optimize.fsolve(_beta_mle_a, a,
                args=(b, len(data), np.log(data).sum()), full_output=True)
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
            s2 = np.log(1 - data).sum()

            # Use the "method of moments" to estimate the initial
            # guess for a and b.
            fac = xbar * (1 - xbar) / data.var(ddof=0) - 1
            a = xbar * fac
            b = (1 - xbar) * fac

            # Compute the MLE for a and b by solving _beta_mle_ab.
            theta, info, ier, mesg = optimize.fsolve(_beta_mle_ab, [a, b],
                args=(len(data), s1, s2), full_output=True)
            if ier != 1:
                raise FitSolverError(mesg=mesg)
            a, b = theta

        return a, b, floc, fscale

beta = beta_gen(a=0.0, b=1.0, name='beta')


class betaprime_gen(rv_continuous):
    """A beta prime continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `betaprime` is::

        betaprime.pdf(x, a, b) = x**(a-1) * (1+x)**(-a-b) / beta(a, b)

    for ``x > 0``, ``a > 0``, ``b > 0``, where ``beta(a, b)`` is the beta
    function (see `scipy.special.beta`).

    %(example)s

    """
    def _rvs(self, a, b):
        u1 = gamma.rvs(a,size=self._size)
        u2 = gamma.rvs(b,size=self._size)
        return (u1 / u2)

    def _pdf(self, x, a, b):
        return np.exp(self._logpdf(x, a, b))

    def _logpdf(self, x, a, b):
        return special.xlogy(a-1.0, x) - special.xlog1py(a+b, x) - special.betaln(a,b)

    def _cdf_skip(self, x, a, b):
        # remove for now: special.hyp2f1 is incorrect for large a
        x = where(x == 1.0, 1.0-1e-6,x)
        return pow(x,a)*special.hyp2f1(a+b,a,1+a,-x)/a/special.beta(a,b)

    def _munp(self, n, a, b):
        if (n == 1.0):
            return where(b > 1, a/(b-1.0), inf)
        elif (n == 2.0):
            return where(b > 2, a*(a+1.0)/((b-2.0)*(b-1.0)), inf)
        elif (n == 3.0):
            return where(b > 3, a*(a+1.0)*(a+2.0)/((b-3.0)*(b-2.0)*(b-1.0)),
                         inf)
        elif (n == 4.0):
            return where(b > 4,
                         a*(a+1.0)*(a+2.0)*(a+3.0)/((b-4.0)*(b-3.0)
                                                    * (b-2.0)*(b-1.0)), inf)
        else:
            raise NotImplementedError
betaprime = betaprime_gen(a=0.0, b=500.0, name='betaprime')


class bradford_gen(rv_continuous):
    """A Bradford continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `bradford` is::

        bradford.pdf(x, c) = c / (k * (1+c*x)),

    for ``0 < x < 1``, ``c > 0`` and ``k = log(1+c)``.

    %(example)s

    """
    def _pdf(self, x, c):
        return c / (c*x + 1.0) / log(1.0+c)

    def _cdf(self, x, c):
        return log(1.0+c*x) / log(c+1.0)

    def _ppf(self, q, c):
        return ((1.0+c)**q-1)/c

    def _stats(self, c, moments='mv'):
        k = log(1.0+c)
        mu = (c-k)/(c*k)
        mu2 = ((c+2.0)*k-2.0*c)/(2*c*k*k)
        g1 = None
        g2 = None
        if 's' in moments:
            g1 = sqrt(2)*(12*c*c-9*c*k*(c+2)+2*k*k*(c*(c+3)+3))
            g1 /= sqrt(c*(c*(k-2)+2*k))*(3*c*(k-2)+6*k)
        if 'k' in moments:
            g2 = c**3*(k-3)*(k*(3*k-16)+24)+12*k*c*c*(k-4)*(k-3) \
                 + 6*c*k*k*(3*k-14) + 12*k**3
            g2 /= 3*c*(c*(k-2)+2*k)**2
        return mu, mu2, g1, g2

    def _entropy(self, c):
        k = log(1+c)
        return k/2.0 - log(c/k)
bradford = bradford_gen(a=0.0, b=1.0, name='bradford')


class burr_gen(rv_continuous):
    """A Burr continuous random variable.

    %(before_notes)s

    See Also
    --------
    fisk : a special case of `burr` with ``d = 1``

    Notes
    -----
    The probability density function for `burr` is::

        burr.pdf(x, c, d) = c * d * x**(-c-1) * (1+x**(-c))**(-d-1)

    for ``x > 0``.

    %(example)s

    """
    def _pdf(self, x, c, d):
        return c*d*(x**(-c-1.0))*((1+x**(-c*1.0))**(-d-1.0))

    def _cdf(self, x, c, d):
        return (1+x**(-c*1.0))**(-d**1.0)

    def _ppf(self, q, c, d):
        return (q**(-1.0/d)-1)**(-1.0/c)

    def _stats(self, c, d, moments='mv'):
        g2c, g2cd = gam(1-2.0/c), gam(2.0/c+d)
        g1c, g1cd = gam(1-1.0/c), gam(1.0/c+d)
        gd = gam(d)
        k = gd*g2c*g2cd - g1c**2 * g1cd**2
        mu = g1c*g1cd / gd
        mu2 = k / gd**2.0
        g1, g2 = None, None
        g3c, g3cd = None, None
        if 's' in moments:
            g3c, g3cd = gam(1-3.0/c), gam(3.0/c+d)
            g1 = 2*g1c**3 * g1cd**3 + gd*gd*g3c*g3cd - 3*gd*g2c*g1c*g1cd*g2cd
            g1 /= sqrt(k**3)
        if 'k' in moments:
            if g3c is None:
                g3c = gam(1-3.0/c)
            if g3cd is None:
                g3cd = gam(3.0/c+d)
            g4c, g4cd = gam(1-4.0/c), gam(4.0/c+d)
            g2 = 6*gd*g2c*g2cd * g1c**2 * g1cd**2 + gd**3 * g4c*g4cd
            g2 -= 3*g1c**4 * g1cd**4 - 4*gd**2*g3c*g1c*g1cd*g3cd
        return mu, mu2, g1, g2
burr = burr_gen(a=0.0, name='burr')

#XXX: cf PR #2552
class fisk_gen(burr_gen):
    """A Fisk continuous random variable.

    The Fisk distribution is also known as the log-logistic distribution, and
    equals the Burr distribution with ``d == 1``.

    %(before_notes)s

    See Also
    --------
    burr

    %(example)s

    """
    def _pdf(self, x, c):
        return burr_gen._pdf(self, x, c, 1.0)

    def _cdf(self, x, c):
        return burr_gen._cdf(self, x, c, 1.0)

    def _ppf(self, x, c):
        return burr_gen._ppf(self, x, c, 1.0)

    def _stats(self, c):
        return burr_gen._stats(self, c, 1.0)

    def _entropy(self, c):
        return 2 - log(c)
fisk = fisk_gen(a=0.0, name='fisk')


# median = loc
class cauchy_gen(rv_continuous):
    """A Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `cauchy` is::

        cauchy.pdf(x) = 1 / (pi * (1 + x**2))

    %(example)s

    """
    def _pdf(self, x):
        return 1.0/pi/(1.0+x*x)

    def _cdf(self, x):
        return 0.5 + 1.0/pi*arctan(x)

    def _ppf(self, q):
        return tan(pi*q-pi/2.0)

    def _sf(self, x):
        return 0.5 - 1.0/pi*arctan(x)

    def _isf(self, q):
        return tan(pi/2.0-pi*q)

    def _stats(self):
        return inf, inf, nan, nan

    def _entropy(self):
        return log(4*pi)

    def _fitstart(self, data, args=None):
        return (0, 1)
cauchy = cauchy_gen(name='cauchy')


class chi_gen(rv_continuous):
    """A chi continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `chi` is::

        chi.pdf(x,df) = x**(df-1) * exp(-x**2/2) / (2**(df/2-1) * gamma(df/2))

    for ``x > 0``.

    Special cases of `chi` are:

        - ``chi(1, loc, scale) = `halfnormal`
        - ``chi(2, 0, scale) = `rayleigh`
        - ``chi(3, 0, scale) : `maxwell`

    %(example)s

    """
    def _rvs(self, df):
        return sqrt(chi2.rvs(df,size=self._size))

    def _pdf(self, x, df):
        return x**(df-1.)*exp(-x*x*0.5)/(2.0)**(df*0.5-1)/gam(df*0.5)

    def _cdf(self, x, df):
        return special.gammainc(df*0.5,0.5*x*x)

    def _ppf(self, q, df):
        return sqrt(2*special.gammaincinv(df*0.5,q))

    def _stats(self, df):
        mu = sqrt(2)*special.gamma(df/2.0+0.5)/special.gamma(df/2.0)
        mu2 = df - mu*mu
        g1 = (2*mu**3.0 + mu*(1-2*df))/asarray(np.power(mu2, 1.5))
        g2 = 2*df*(1.0-df)-6*mu**4 + 4*mu**2 * (2*df-1)
        g2 /= asarray(mu2**2.0)
        return mu, mu2, g1, g2
chi = chi_gen(a=0.0, name='chi')


## Chi-squared (gamma-distributed with loc=0 and scale=2 and shape=df/2)
class chi2_gen(rv_continuous):
    """A chi-squared continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `chi2` is::

        chi2.pdf(x,df) = 1 / (2*gamma(df/2)) * (x/2)**(df/2-1) * exp(-x/2)

    %(example)s

    """
    def _rvs(self, df):
        return mtrand.chisquare(df,self._size)

    def _pdf(self, x, df):
        return exp(self._logpdf(x, df))

    def _logpdf(self, x, df):
        # term1 = (df/2.-1)*log(x)
        # term1[(df==2)*(x==0)] = 0
        # avoid 0*log(0)==nan
        return (df/2.-1)*log(x+1e-300) - x/2. - gamln(df/2.) - (log(2)*df)/2.

    def _cdf(self, x, df):
        return special.chdtr(df, x)

    def _sf(self, x, df):
        return special.chdtrc(df, x)

    def _isf(self, p, df):
        return special.chdtri(df, p)

    def _ppf(self, p, df):
        return self._isf(1.0-p, df)

    def _stats(self, df):
        mu = df
        mu2 = 2*df
        g1 = 2*sqrt(2.0/df)
        g2 = 12.0/df
        return mu, mu2, g1, g2
chi2 = chi2_gen(a=0.0, name='chi2')


class cosine_gen(rv_continuous):
    """A cosine continuous random variable.

    %(before_notes)s

    Notes
    -----
    The cosine distribution is an approximation to the normal distribution.
    The probability density function for `cosine` is::

        cosine.pdf(x) = 1/(2*pi) * (1+cos(x))

    for ``-pi <= x <= pi``.

    %(example)s

    """
    def _pdf(self, x):
        return 1.0/2/pi*(1+cos(x))

    def _cdf(self, x):
        return 1.0/2/pi*(pi + x + sin(x))

    def _stats(self):
        return 0.0, pi*pi/3.0-2.0, 0.0, -6.0*(pi**4-90)/(5.0*(pi*pi-6)**2)

    def _entropy(self):
        return log(4*pi)-1.0
cosine = cosine_gen(a=-pi, b=pi, name='cosine')


class dgamma_gen(rv_continuous):
    """A double gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `dgamma` is::

        dgamma.pdf(x, a) = 1 / (2*gamma(a)) * abs(x)**(a-1) * exp(-abs(x))

    for ``a > 0``.

    %(example)s

    """
    def _rvs(self, a):
        u = random(size=self._size)
        return (gamma.rvs(a,size=self._size)*where(u >= 0.5,1,-1))

    def _pdf(self, x, a):
        ax = abs(x)
        return 1.0/(2*special.gamma(a))*ax**(a-1.0) * exp(-ax)

    def _logpdf(self, x, a):
        ax = abs(x)
        return (a-1.0)*log(ax) - ax - log(2) - gamln(a)

    def _cdf(self, x, a):
        fac = 0.5*special.gammainc(a,abs(x))
        return where(x > 0,0.5+fac,0.5-fac)

    def _sf(self, x, a):
        fac = 0.5*special.gammainc(a,abs(x))
        return where(x > 0,0.5-fac,0.5+fac)

    def _ppf(self, q, a):
        fac = special.gammainccinv(a,1-abs(2*q-1))
        return where(q > 0.5, fac, -fac)

    def _stats(self, a):
        mu2 = a*(a+1.0)
        return 0.0, mu2, 0.0, (a+2.0)*(a+3.0)/mu2-3.0
dgamma = dgamma_gen(name='dgamma')


class dweibull_gen(rv_continuous):
    """A double Weibull continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `dweibull` is::

        dweibull.pdf(x, c) = c / 2 * abs(x)**(c-1) * exp(-abs(x)**c)

    %(example)s

    """
    def _rvs(self, c):
        u = random(size=self._size)
        return weibull_min.rvs(c, size=self._size)*(where(u >= 0.5,1,-1))

    def _pdf(self, x, c):
        ax = abs(x)
        Px = c/2.0*ax**(c-1.0)*exp(-ax**c)
        return Px

    def _logpdf(self, x, c):
        ax = abs(x)
        return log(c) - log(2.0) + (c-1.0)*log(ax) - ax**c

    def _cdf(self, x, c):
        Cx1 = 0.5*exp(-abs(x)**c)
        return where(x > 0, 1-Cx1, Cx1)

    def _ppf_skip(self, q, c):
        fac = where(q <= 0.5,2*q,2*q-1)
        fac = pow(asarray(log(1.0/fac)),1.0/c)
        return where(q > 0.5,fac,-fac)

    def _stats(self, c):
        var = gam(1+2.0/c)
        return 0.0, var, 0.0, gam(1+4.0/c)/var
dweibull = dweibull_gen(name='dweibull')


## Exponential (gamma distributed with a=1.0, loc=loc and scale=scale)
class expon_gen(rv_continuous):
    """An exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `expon` is::

        expon.pdf(x) = lambda * exp(- lambda*x)

    for ``x >= 0``.

    The scale parameter is equal to ``scale = 1.0 / lambda``.

    `expon` does not have shape parameters.

    %(example)s

    """
    def _rvs(self):
        return mtrand.standard_exponential(self._size)

    def _pdf(self, x):
        return exp(-x)

    def _logpdf(self, x):
        return -x

    def _cdf(self, x):
        return -expm1(-x)

    def _ppf(self, q):
        return -log1p(-q)

    def _sf(self,x):
        return exp(-x)

    def _logsf(self, x):
        return -x

    def _isf(self,q):
        return -log(q)

    def _stats(self):
        return 1.0, 1.0, 2.0, 6.0

    def _entropy(self):
        return 1.0
expon = expon_gen(a=0.0, name='expon')


class exponweib_gen(rv_continuous):
    """An exponentiated Weibull continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponweib` is::

        exponweib.pdf(x, a, c) =
            a * c * (1-exp(-x**c))**(a-1) * exp(-x**c)*x**(c-1)

    for ``x > 0``, ``a > 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, a, c):
        exc = exp(-x**c)
        return a*c*(1-exc)**asarray(a-1) * exc * x**(c-1)

    def _logpdf(self, x, a, c):
        exc = exp(-x**c)
        return log(a) + log(c) + (a-1.)*log(1-exc) - x**c + (c-1.0)*log(x)

    def _cdf(self, x, a, c):
        exm1c = -expm1(-x**c)
        return (exm1c)**a

    def _ppf(self, q, a, c):
        return (-log1p(-q**(1.0/a)))**asarray(1.0/c)
exponweib = exponweib_gen(a=0.0, name='exponweib')


class exponpow_gen(rv_continuous):
    """An exponential power continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponpow` is::

        exponpow.pdf(x, b) = b * x**(b-1) * exp(1+x**b - exp(x**b))

    for ``x >= 0``, ``b > 0``.

    %(example)s

    """
    def _pdf(self, x, b):
        xbm1 = x**(b-1.0)
        xb = xbm1 * x
        return exp(1)*b*xbm1 * exp(xb - exp(xb))

    def _logpdf(self, x, b):
        xb = x**(b-1.0)*x
        return 1 + log(b) + (b-1.0)*log(x) + xb - exp(xb)

    def _cdf(self, x, b):
        return -expm1(-expm1(x**b))

    def _sf(self, x, b):
        return exp(-expm1(x**b))

    def _isf(self, x, b):
        return (log1p(-log(x)))**(1./b)

    def _ppf(self, q, b):
        return pow(log1p(-log1p(-q)), 1.0/b)
exponpow = exponpow_gen(a=0.0, name='exponpow')


class fatiguelife_gen(rv_continuous):
    """A fatigue-life (Birnbaum-Sanders) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `fatiguelife` is::

        fatiguelife.pdf(x,c) =
            (x+1) / (2*c*sqrt(2*pi*x**3)) * exp(-(x-1)**2/(2*x*c**2))

    for ``x > 0``.

    %(example)s

    """
    def _rvs(self, c):
        z = norm.rvs(size=self._size)
        x = 0.5*c*z
        x2 = x*x
        t = 1.0 + 2*x2 + 2*x*sqrt(1 + x2)
        return t

    def _pdf(self, x, c):
        return (x+1)/asarray(2*c*sqrt(2*pi*x**3))*exp(-(x-1)**2/asarray((2.0*x*c**2)))

    def _logpdf(self, x, c):
        return log(x+1) - (x-1)**2 / (2.0*x*c**2) - log(2*c) - 0.5*(log(2*pi) + 3*log(x))

    def _cdf(self, x, c):
        return special.ndtr(1.0/c*(sqrt(x)-1.0/asarray(sqrt(x))))

    def _ppf(self, q, c):
        tmp = c*special.ndtri(q)
        return 0.25*(tmp + sqrt(tmp**2 + 4))**2

    def _stats(self, c):
        c2 = c*c
        mu = c2 / 2.0 + 1
        den = 5*c2 + 4
        mu2 = c2*den / 4.0
        g1 = 4*c*sqrt(11*c2+6.0)/np.power(den, 1.5)
        g2 = 6*c2*(93*c2+41.0) / den**2.0
        return mu, mu2, g1, g2
fatiguelife = fatiguelife_gen(a=0.0, name='fatiguelife')


class foldcauchy_gen(rv_continuous):
    """A folded Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `foldcauchy` is::

        foldcauchy.pdf(x, c) = 1/(pi*(1+(x-c)**2)) + 1/(pi*(1+(x+c)**2))

    for ``x >= 0``.

    %(example)s

    """
    def _rvs(self, c):
        return abs(cauchy.rvs(loc=c,size=self._size))

    def _pdf(self, x, c):
        return 1.0/pi*(1.0/(1+(x-c)**2) + 1.0/(1+(x+c)**2))

    def _cdf(self, x, c):
        return 1.0/pi*(arctan(x-c) + arctan(x+c))

    def _stats(self, c):
        return inf, inf, nan, nan
foldcauchy = foldcauchy_gen(a=0.0, name='foldcauchy')


class f_gen(rv_continuous):
    """An F continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `f` is::

                             df2**(df2/2) * df1**(df1/2) * x**(df1/2-1)
        F.pdf(x, df1, df2) = --------------------------------------------
                             (df2+df1*x)**((df1+df2)/2) * B(df1/2, df2/2)

    for ``x > 0``.

    %(example)s

    """
    def _rvs(self, dfn, dfd):
        return mtrand.f(dfn, dfd, self._size)

    def _pdf(self, x, dfn, dfd):
        return exp(self._logpdf(x, dfn, dfd))

    def _logpdf(self, x, dfn, dfd):
        n = 1.0*dfn
        m = 1.0*dfd
        lPx = m/2*log(m) + n/2*log(n) + (n/2-1)*log(x)
        lPx -= ((n+m)/2)*log(m+n*x) + special.betaln(n/2,m/2)
        return lPx

    def _cdf(self, x, dfn, dfd):
        return special.fdtr(dfn, dfd, x)

    def _sf(self, x, dfn, dfd):
        return special.fdtrc(dfn, dfd, x)

    def _ppf(self, q, dfn, dfd):
        return special.fdtri(dfn, dfd, q)

    def _stats(self, dfn, dfd):
        v2 = asarray(dfd*1.0)
        v1 = asarray(dfn*1.0)
        mu = where(v2 > 2, v2 / asarray(v2 - 2), inf)
        mu2 = 2*v2*v2*(v2+v1-2)/(v1*(v2-2)**2 * (v2-4))
        mu2 = where(v2 > 4, mu2, inf)
        g1 = 2*(v2+2*v1-2)/(v2-6)*sqrt((2*v2-4)/(v1*(v2+v1-2)))
        g1 = where(v2 > 6, g1, nan)
        g2 = 3/(2*v2-16)*(8+g1*g1*(v2-6))
        g2 = where(v2 > 8, g2, nan)
        return mu, mu2, g1, g2
f = f_gen(a=0.0, name='f')


## Folded Normal
##   abs(Z) where (Z is normal with mu=L and std=S so that c=abs(L)/S)
##
##  note: regress docs have scale parameter correct, but first parameter
##    he gives is a shape parameter A = c * scale

##  Half-normal is folded normal with shape-parameter c=0.

class foldnorm_gen(rv_continuous):
    """A folded normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `foldnorm` is::

        foldnormal.pdf(x, c) = sqrt(2/pi) * cosh(c*x) * exp(-(x**2+c**2)/2)

    for ``c >= 0``.

    %(example)s

    """
    def _argcheck(self, c):
        return (c >= 0)

    def _rvs(self, c):
        return abs(norm.rvs(loc=c,size=self._size))

    def _pdf(self, x, c):
        return sqrt(2.0/pi)*cosh(c*x)*exp(-(x*x+c*c)/2.0)

    def _cdf(self, x, c):
        return special.ndtr(x-c) + special.ndtr(x+c) - 1.0

    def _stats(self, c):
        fac = special.erf(c/sqrt(2))
        mu = sqrt(2.0/pi)*exp(-0.5*c*c)+c*fac
        mu2 = c*c + 1 - mu*mu
        c2 = c*c
        g1 = sqrt(2/pi)*exp(-1.5*c2)*(4-pi*exp(c2)*(2*c2+1.0))
        g1 += 2*c*fac*(6*exp(-c2) + 3*sqrt(2*pi)*c*exp(-c2/2.0)*fac +
                       pi*c*(fac*fac-1))
        g1 /= pi*np.power(mu2, 1.5)

        g2 = c2*c2+6*c2+3+6*(c2+1)*mu*mu - 3*mu**4
        g2 -= 4*exp(-c2/2.0)*mu*(sqrt(2.0/pi)*(c2+2)+c*(c2+3)*exp(c2/2.0)*fac)
        g2 /= mu2**2.0
        return mu, mu2, g1, g2
foldnorm = foldnorm_gen(a=0.0, name='foldnorm')


## Extreme Value Type II or Frechet
## (defined in Regress+ documentation as Extreme LB) as
##   a limiting value distribution.
##
class frechet_r_gen(rv_continuous):
    """A Frechet right (or Weibull minimum) continuous random variable.

    %(before_notes)s

    See Also
    --------
    weibull_min : The same distribution as `frechet_r`.
    frechet_l, weibull_max

    Notes
    -----
    The probability density function for `frechet_r` is::

        frechet_r.pdf(x, c) = c * x**(c-1) * exp(-x**c)

    for ``x > 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        return c*pow(x,c-1)*exp(-pow(x,c))

    def _logpdf(self, x, c):
        return log(c) + (c-1)*log(x) - pow(x,c)

    def _cdf(self, x, c):
        return -expm1(-pow(x,c))

    def _ppf(self, q, c):
        return pow(-log1p(-q),1.0/c)

    def _munp(self, n, c):
        return special.gamma(1.0+n*1.0/c)

    def _entropy(self, c):
        return -_EULER / c - log(c) + _EULER + 1
frechet_r = frechet_r_gen(a=0.0, name='frechet_r')
weibull_min = frechet_r_gen(a=0.0, name='weibull_min')


class frechet_l_gen(rv_continuous):
    """A Frechet left (or Weibull maximum) continuous random variable.

    %(before_notes)s

    See Also
    --------
    weibull_max : The same distribution as `frechet_l`.
    frechet_r, weibull_min

    Notes
    -----
    The probability density function for `frechet_l` is::

        frechet_l.pdf(x, c) = c * (-x)**(c-1) * exp(-(-x)**c)

    for ``x < 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        return c*pow(-x,c-1)*exp(-pow(-x,c))

    def _cdf(self, x, c):
        return exp(-pow(-x,c))

    def _ppf(self, q, c):
        return -pow(-log(q),1.0/c)

    def _munp(self, n, c):
        val = special.gamma(1.0+n*1.0/c)
        if (int(n) % 2):
            sgn = -1
        else:
            sgn = 1
        return sgn * val

    def _entropy(self, c):
        return -_EULER / c - log(c) + _EULER + 1
frechet_l = frechet_l_gen(b=0.0, name='frechet_l')
weibull_max = frechet_l_gen(b=0.0, name='weibull_max')


class genlogistic_gen(rv_continuous):
    """A generalized logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genlogistic` is::

        genlogistic.pdf(x, c) = c * exp(-x) / (1 + exp(-x))**(c+1)

    for ``x > 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        Px = c*exp(-x)/(1+exp(-x))**(c+1.0)
        return Px

    def _logpdf(self, x, c):
        return log(c) - x - (c+1.0)*log1p(exp(-x))

    def _cdf(self, x, c):
        Cx = (1+exp(-x))**(-c)
        return Cx

    def _ppf(self, q, c):
        vals = -log(pow(q,-1.0/c)-1)
        return vals

    def _stats(self, c):
        zeta = special.zeta
        mu = _EULER + special.psi(c)
        mu2 = pi*pi/6.0 + zeta(2,c)
        g1 = -2*zeta(3,c) + 2*_ZETA3
        g1 /= np.power(mu2, 1.5)
        g2 = pi**4/15.0 + 6*zeta(4,c)
        g2 /= mu2**2.0
        return mu, mu2, g1, g2
genlogistic = genlogistic_gen(name='genlogistic')


class genpareto_gen(rv_continuous):
    """A generalized Pareto continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genpareto` is::

        genpareto.pdf(x, c) = (1 + c * x)**(-1 - 1/c)

    for ``c != 0``, and for ``x >= 0`` for all c,
    and ``x < 1/abs(c)`` for ``c < 0``.

    %(example)s

    """
    def _argcheck(self, c):
        c = asarray(c)
        self.b = where(c < 0, 1.0/abs(c), inf)
        return where(c == 0, 0, 1)

    def _pdf(self, x, c):
        Px = pow(1+c*x,asarray(-1.0-1.0/c))
        return Px

    def _logpdf(self, x, c):
        return (-1.0-1.0/c) * np.log1p(c*x)

    def _cdf(self, x, c):
        return 1.0 - pow(1+c*x,asarray(-1.0/c))

    def _ppf(self, q, c):
        vals = 1.0/c * (pow(1-q, -c)-1)
        return vals

    def _munp(self, n, c):
        k = arange(0,n+1)
        val = (-1.0/c)**n * sum(comb(n,k)*(-1)**k / (1.0-c*k),axis=0)
        return where(c*n < 1, val, inf)

    def _entropy(self, c):
        if (c > 0):
            return 1+c
        else:
            self.b = -1.0 / c
            return rv_continuous._entropy(self, c)
genpareto = genpareto_gen(a=0.0, name='genpareto')


class genexpon_gen(rv_continuous):
    """A generalized exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genexpon` is::

        genexpon.pdf(x, a, b, c) = (a + b * (1 - exp(-c*x))) * \
                                   exp(-a*x - b*x + b/c * (1-exp(-c*x)))

    for ``x >= 0``, ``a,b,c > 0``.

    References
    ----------
    H.K. Ryu, "An Extension of Marshall and Olkin's Bivariate Exponential
    Distribution", Journal of the American Statistical Association, 1993.

     N. Balakrishnan, "The Exponential Distribution: Theory, Methods and
     Applications", Asit P. Basu.

    %(example)s

    """
    def _pdf(self, x, a, b, c):
        return (a+b*(-expm1(-c*x)))*exp((-a-b)*x+b*(-expm1(-c*x))/c)

    def _cdf(self, x, a, b, c):
        return -expm1((-a-b)*x + b*(-expm1(-c*x))/c)

    def _logpdf(self, x, a, b, c):
        return np.log(a+b*(-expm1(-c*x))) + (-a-b)*x+b*(-expm1(-c*x))/c
genexpon = genexpon_gen(a=0.0, name='genexpon')


class genextreme_gen(rv_continuous):
    """A generalized extreme value continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_r

    Notes
    -----
    For ``c=0``, `genextreme` is equal to `gumbel_r`.
    The probability density function for `genextreme` is::

        genextreme.pdf(x, c) =
            exp(-exp(-x))*exp(-x),                    for c==0
            exp(-(1-c*x)**(1/c))*(1-c*x)**(1/c-1),    for x <= 1/c, c > 0

    %(example)s

    """
    def _argcheck(self, c):
        min = np.minimum
        max = np.maximum
        sml = floatinfo.machar.xmin
        self.b = where(c > 0, 1.0 / max(c, sml),inf)
        self.a = where(c < 0, 1.0 / min(c,-sml), -inf)
        return where(abs(c) == inf, 0, 1)

    def _pdf(self, x, c):
        cx = c*x
        logex2 = where((c == 0)*(x == x),0.0,log1p(-cx))
        logpex2 = where((c == 0)*(x == x),-x,logex2/c)
        pex2 = exp(logpex2)
        # Handle special cases
        logpdf = where((cx == 1) | (cx == -inf),-inf,-pex2+logpex2-logex2)
        putmask(logpdf, (c == 1) & (x == 1), 0.0)
        return exp(logpdf)

    def _cdf(self, x, c):
        loglogcdf = where((c == 0)*(x == x),-x,log1p(-c*x)/c)
        return exp(-exp(loglogcdf))

    def _ppf(self, q, c):
        x = -log(-log(q))
        return where((c == 0)*(x == x),x,-expm1(-c*x)/c)

    def _stats(self, c):
        g = lambda n: gam(n*c+1)
        g1 = g(1)
        g2 = g(2)
        g3 = g(3)
        g4 = g(4)
        g2mg12 = where(abs(c) < 1e-7,(c*pi)**2.0/6.0,g2-g1**2.0)
        gam2k = where(abs(c) < 1e-7,pi**2.0/6.0, expm1(gamln(2.0*c+1.0)-2*gamln(c+1.0))/c**2.0)
        eps = 1e-14
        gamk = where(abs(c) < eps,-_EULER,expm1(gamln(c+1))/c)

        m = where(c < -1.0,nan,-gamk)
        v = where(c < -0.5,nan,g1**2.0*gam2k)

        # skewness
        sk1 = where(c < -1./3,nan,np.sign(c)*(-g3+(g2+2*g2mg12)*g1)/((g2mg12)**(3./2.)))
        sk = where(abs(c) <= eps**0.29,12*sqrt(6)*_ZETA3/pi**3,sk1)

        # kurtosis
        ku1 = where(c < -1./4,nan,(g4+(-4*g3+3*(g2+g2mg12)*g1)*g1)/((g2mg12)**2))
        ku = where(abs(c) <= (eps)**0.23,12.0/5.0,ku1-3.0)
        return m,v,sk,ku

    def _munp(self, n, c):
        k = arange(0,n+1)
        vals = 1.0/c**n * sum(comb(n,k) * (-1)**k * special.gamma(c*k + 1),axis=0)
        return where(c*n > -1, vals, inf)
genextreme = genextreme_gen(name='genextreme')


def _digammainv(y):
    # Inverse of the digamma function (real positive arguments only).
    # This function is used in the `fit` method of `gamma_gen`.
    # The function uses either optimize.fsolve or optimize.newton
    # to solve `digamma(x) - y = 0`.  There is probably room for
    # improvement, but currently it works over a wide range of y:
    #    >>> y = 64*np.random.randn(1000000)
    #    >>> y.min(), y.max()
    #    (-311.43592651416662, 351.77388222276869)
    #    x = [_digammainv(t) for t in y]
    #    np.abs(digamma(x) - y).max()
    #    1.1368683772161603e-13
    #
    _em = 0.5772156649015328606065120
    func = lambda x: special.digamma(x) - y
    if y > -0.125:
        x0 = exp(y) + 0.5
        if y < 10:
            # Some experimentation shows that newton reliably converges
            # must faster than fsolve in this y range.  For larger y,
            # newton sometimes fails to converge.
            value = optimize.newton(func, x0, tol=1e-10)
            return value
    elif y > -3:
        x0 = exp(y/2.332) + 0.08661
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
    """A gamma continuous random variable.

    %(before_notes)s

    See Also
    --------
    erlang, expon

    Notes
    -----
    The probability density function for `gamma` is::

        gamma.pdf(x, a) = lambda**a * x**(a-1) * exp(-lambda*x) / gamma(a)

    for ``x >= 0``, ``a > 0``. Here ``gamma(a)`` refers to the gamma function.

    The scale parameter is equal to ``scale = 1.0 / lambda``.

    `gamma` has a shape parameter `a` which needs to be set explicitly. For instance:

        >>> from scipy.stats import gamma
        >>> rv = gamma(3., loc = 0., scale = 2.)

    produces a frozen form of `gamma` with shape ``a = 3.``, ``loc =0.``
    and ``lambda = 1./scale = 1./2.``.

    When ``a`` is an integer, `gamma` reduces to the Erlang
    distribution, and when ``a=1`` to the exponential distribution.

    %(example)s

    """
    def _rvs(self, a):
        return mtrand.standard_gamma(a, self._size)

    def _pdf(self, x, a):
        return exp(self._logpdf(x, a))

    def _logpdf(self, x, a):
        return special.xlogy(a-1.0, x) - x - gamln(a)

    def _cdf(self, x, a):
        return special.gammainc(a, x)

    def _sf(self, x, a):
        return special.gammaincc(a, x)

    def _ppf(self, q, a):
        return special.gammaincinv(a,q)

    def _stats(self, a):
        return a, a, 2.0/sqrt(a), 6.0/a

    def _entropy(self, a):
        return special.psi(a)*(1-a) + 1 + gamln(a)

    def _fitstart(self, data):
        # The skewness of the gamma distribution is `4 / sqrt(a)`.
        # We invert that to estimate the shape `a` using the skewness
        # of the data.  The formula is regularized with 1e-8 in the
        # denominator to allow for degenerate data where the skewness
        # is close to 0.
        a = 4 / (1e-8 + _skew(data)**2)
        return super(gamma_gen, self)._fitstart(data, args=(a,))

    @inherit_docstring_from(rv_continuous)
    def fit(self, data, *args, **kwds):
        f0 = kwds.get('f0', None)
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is None:
            # loc is not fixed.  Use the default fit method.
            return super(gamma_gen, self).fit(data, *args, **kwds)

        # Special case: loc is fixed.

        if f0 is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            # Without this check, this function would just return the
            # parameters that were given.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        # Fixed location is handled by shifting the data.
        data = np.asarray(data)
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
                # log(a) - special.digamma(a) - log(xbar) + log(data.mean) = 0
                s = log(xbar) - log(data).mean()
                func = lambda a: log(a) - special.digamma(a) - s
                aest = (3-s + math.sqrt((s-3)**2 + 24*s)) / (12*s)
                xa = aest*(1-0.4)
                xb = aest*(1+0.4)
                a = optimize.brentq(func, xa, xb, disp=0)

            # The MLE for the scale parameter is just the data mean
            # divided by the shape parameter.
            scale = xbar / a
        else:
            # scale is fixed, shape is free
            # The MLE for the shape parameter `a` is the solution to:
            # special.digamma(a) - log(data).mean() + log(fscale) = 0
            c = log(data).mean() - log(fscale)
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
        allpos = np.all(a > 0)
        if not allint:
            # An Erlang distribution shouldn't really have a non-integer
            # shape parameter, so warn the user.
            warnings.warn('The shape parameter of the erlang distribution '
                'has been given a non-integer value %r.' % (a,),
                RuntimeWarning)
        return allpos

    def _fitstart(self, data):
        # Override gamma_gen_fitstart so that an integer initial value is
        # used.  (Also regularize the division, to avoid issues when
        # _skew(data) is 0 or close to 0.)
        a = int(4.0 / (1e-8 + _skew(data)**2))
        return super(gamma_gen, self)._fitstart(data, args=(a,))

    # Trivial override of the fit method, so we can monkey-patch its
    # docstring.
    def fit(self, data, *args, **kwds):
        return super(erlang_gen, self).fit(data, *args, **kwds)

    if fit.__doc__ is not None:
        fit.__doc__ = (rv_continuous.fit.__doc__ +
            """
            Notes
            -----
            The Erlang distribution is generally defined to have integer values
            for the shape parameter.  This is not enforced by the `erlang` class.
            When fitting the distribution, it will generally return a non-integer
            value for the shape parameter.  By using the keyword argument
            `f0=<integer>`, the fit method can be constrained to fit the data to
            a specific integer shape parameter.
            """)
erlang = erlang_gen(a=0.0, name='erlang')


class gengamma_gen(rv_continuous):
    """A generalized gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gengamma` is::

        gengamma.pdf(x, a, c) = abs(c) * x**(c*a-1) * exp(-x**c) / gamma(a)

    for ``x > 0``, ``a > 0``, and ``c != 0``.

    %(example)s

    """
    def _argcheck(self, a, c):
        return (a > 0) & (c != 0)

    def _pdf(self, x, a, c):
        return abs(c) * exp((c*a-1)*log(x)-x**c - gamln(a))

    def _cdf(self, x, a, c):
        val = special.gammainc(a,x**c)
        cond = c + 0*val
        return where(cond > 0,val,1-val)

    def _ppf(self, q, a, c):
        val1 = special.gammaincinv(a,q)
        val2 = special.gammaincinv(a,1.0-q)
        ic = 1.0/c
        cond = c+0*val1
        return where(cond > 0,val1**ic,val2**ic)

    def _munp(self, n, a, c):
        return special.gamma(a+n*1.0/c) / special.gamma(a)

    def _entropy(self, a,c):
        val = special.psi(a)
        return a*(1-val) + 1.0/c*val + gamln(a)-log(abs(c))
gengamma = gengamma_gen(a=0.0, name='gengamma')


class genhalflogistic_gen(rv_continuous):
    """A generalized half-logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genhalflogistic` is::

        genhalflogistic.pdf(x, c) = 2 * (1-c*x)**(1/c-1) / (1+(1-c*x)**(1/c))**2

    for ``0 <= x <= 1/c``, and ``c > 0``.

    %(example)s

    """
    def _argcheck(self, c):
        self.b = 1.0 / c
        return (c > 0)

    def _pdf(self, x, c):
        limit = 1.0/c
        tmp = asarray(1-c*x)
        tmp0 = tmp**(limit-1)
        tmp2 = tmp0*tmp
        return 2*tmp0 / (1+tmp2)**2

    def _cdf(self, x, c):
        limit = 1.0/c
        tmp = asarray(1-c*x)
        tmp2 = tmp**(limit)
        return (1.0-tmp2) / (1+tmp2)

    def _ppf(self, q, c):
        return 1.0/c*(1-((1.0-q)/(1.0+q))**c)

    def _entropy(self,c):
        return 2 - (2*c+1)*log(2)
genhalflogistic = genhalflogistic_gen(a=0.0, name='genhalflogistic')


class gompertz_gen(rv_continuous):
    """A Gompertz (or truncated Gumbel) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gompertz` is::

        gompertz.pdf(x, c) = c * exp(x) * exp(-c*(exp(x)-1))

    for ``x >= 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        ex = exp(x)
        return c*ex*exp(-c*(ex-1))

    def _cdf(self, x, c):
        return 1.0-exp(-c*(exp(x)-1))

    def _ppf(self, q, c):
        return log(1-1.0/c*log(1-q))

    def _entropy(self, c):
        return 1.0 - log(c) - exp(c)*special.expn(1,c)
gompertz = gompertz_gen(a=0.0, name='gompertz')


class gumbel_r_gen(rv_continuous):
    """A right-skewed Gumbel continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_l, gompertz, genextreme

    Notes
    -----
    The probability density function for `gumbel_r` is::

        gumbel_r.pdf(x) = exp(-(x + exp(-x)))

    The Gumbel distribution is sometimes referred to as a type I Fisher-Tippett
    distribution.  It is also related to the extreme value distribution,
    log-Weibull and Gompertz distributions.

    %(example)s

    """
    def _pdf(self, x):
        ex = exp(-x)
        return ex*exp(-ex)

    def _logpdf(self, x):
        return -x - exp(-x)

    def _cdf(self, x):
        return exp(-exp(-x))

    def _logcdf(self, x):
        return -exp(-x)

    def _ppf(self, q):
        return -log(-log(q))

    def _stats(self):
        return _EULER, pi*pi/6.0, \
               12*sqrt(6)/pi**3 * _ZETA3, 12.0/5

    def _entropy(self):
        return 1.0608407169541684911
gumbel_r = gumbel_r_gen(name='gumbel_r')


class gumbel_l_gen(rv_continuous):
    """A left-skewed Gumbel continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_r, gompertz, genextreme

    Notes
    -----
    The probability density function for `gumbel_l` is::

        gumbel_l.pdf(x) = exp(x - exp(x))

    The Gumbel distribution is sometimes referred to as a type I Fisher-Tippett
    distribution.  It is also related to the extreme value distribution,
    log-Weibull and Gompertz distributions.

    %(example)s

    """
    def _pdf(self, x):
        ex = exp(x)
        return ex*exp(-ex)

    def _logpdf(self, x):
        return x - exp(x)

    def _cdf(self, x):
        return 1.0-exp(-exp(x))

    def _ppf(self, q):
        return log(-log(1-q))

    def _stats(self):
        return -_EULER, pi*pi/6.0, \
               -12*sqrt(6)/pi**3 * _ZETA3, 12.0/5

    def _entropy(self):
        return 1.0608407169541684911
gumbel_l = gumbel_l_gen(name='gumbel_l')


class halfcauchy_gen(rv_continuous):
    """A Half-Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfcauchy` is::

        halfcauchy.pdf(x) = 2 / (pi * (1 + x**2))

    for ``x >= 0``.

    %(example)s

    """
    def _pdf(self, x):
        return 2.0/pi/(1.0+x*x)

    def _logpdf(self, x):
        return np.log(2.0/pi) - np.log1p(x*x)

    def _cdf(self, x):
        return 2.0/pi*arctan(x)

    def _ppf(self, q):
        return tan(pi/2*q)

    def _stats(self):
        return inf, inf, nan, nan

    def _entropy(self):
        return log(2*pi)
halfcauchy = halfcauchy_gen(a=0.0, name='halfcauchy')


class halflogistic_gen(rv_continuous):
    """A half-logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halflogistic` is::

        halflogistic.pdf(x) = 2 * exp(-x) / (1+exp(-x))**2 = 1/2 * sech(x/2)**2

    for ``x >= 0``.

    %(example)s

    """
    def _pdf(self, x):
        return 0.5/(cosh(x/2.0))**2.0

    def _cdf(self, x):
        return tanh(x/2.0)

    def _ppf(self, q):
        return 2*arctanh(q)

    def _munp(self, n):
        if n == 1:
            return 2*log(2)
        if n == 2:
            return pi*pi/3.0
        if n == 3:
            return 9*_ZETA3
        if n == 4:
            return 7*pi**4 / 15.0
        return 2*(1-pow(2.0,1-n))*special.gamma(n+1)*special.zeta(n,1)

    def _entropy(self):
        return 2-log(2)
halflogistic = halflogistic_gen(a=0.0, name='halflogistic')


class halfnorm_gen(rv_continuous):
    """A half-normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfnorm` is::

        halfnorm.pdf(x) = sqrt(2/pi) * exp(-x**2/2)

    for ``x > 0``.

    `halfnorm` is a special case of `chi` with ``df == 1``.

    %(example)s

    """
    def _rvs(self):
        return abs(norm.rvs(size=self._size))

    def _pdf(self, x):
        return sqrt(2.0/pi)*exp(-x*x/2.0)

    def _logpdf(self, x):
        return 0.5 * np.log(2.0/pi) - x*x/2.0

    def _cdf(self, x):
        return special.ndtr(x)*2-1.0

    def _ppf(self, q):
        return special.ndtri((1+q)/2.0)

    def _stats(self):
        return sqrt(2.0/pi), 1-2.0/pi, sqrt(2)*(4-pi)/(pi-2)**1.5, \
               8*(pi-3)/(pi-2)**2

    def _entropy(self):
        return 0.5*log(pi/2.0)+0.5
halfnorm = halfnorm_gen(a=0.0, name='halfnorm')


class hypsecant_gen(rv_continuous):
    """A hyperbolic secant continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `hypsecant` is::

        hypsecant.pdf(x) = 1/pi * sech(x)

    %(example)s

    """
    def _pdf(self, x):
        return 1.0/(pi*cosh(x))

    def _cdf(self, x):
        return 2.0/pi*arctan(exp(x))

    def _ppf(self, q):
        return log(tan(pi*q/2.0))

    def _stats(self):
        return 0, pi*pi/4, 0, 2

    def _entropy(self):
        return log(2*pi)
hypsecant = hypsecant_gen(name='hypsecant')


class gausshyper_gen(rv_continuous):
    """A Gauss hypergeometric continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gausshyper` is::

        gausshyper.pdf(x, a, b, c, z) =
            C * x**(a-1) * (1-x)**(b-1) * (1+z*x)**(-c)

    for ``0 <= x <= 1``, ``a > 0``, ``b > 0``, and
    ``C = 1 / (B(a,b) F[2,1](c, a; a+b; -z))``

    %(example)s

    """
    def _argcheck(self, a, b, c, z):
        return (a > 0) & (b > 0) & (c == c) & (z == z)

    def _pdf(self, x, a, b, c, z):
        Cinv = gam(a)*gam(b)/gam(a+b)*special.hyp2f1(c,a,a+b,-z)
        return 1.0/Cinv * x**(a-1.0) * (1.0-x)**(b-1.0) / (1.0+z*x)**c

    def _munp(self, n, a, b, c, z):
        fac = special.beta(n+a,b) / special.beta(a,b)
        num = special.hyp2f1(c,a+n,a+b+n,-z)
        den = special.hyp2f1(c,a,a+b,-z)
        return fac*num / den
gausshyper = gausshyper_gen(a=0.0, b=1.0, name='gausshyper')


class invgamma_gen(rv_continuous):
    """An inverted gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invgamma` is::

        invgamma.pdf(x, a) = x**(-a-1) / gamma(a) * exp(-1/x)

    for x > 0, a > 0.

    `invgamma` is a special case of `gengamma` with ``c == -1``.

    %(example)s

    """
    def _pdf(self, x, a):
        return exp(self._logpdf(x,a))

    def _logpdf(self, x, a):
        return (-(a+1)*log(x)-gamln(a) - 1.0/x)

    def _cdf(self, x, a):
        return 1.0-special.gammainc(a, 1.0/x)

    def _ppf(self, q, a):
        return 1.0/special.gammaincinv(a,1-q)

    def _munp(self, n, a):
        return exp(gamln(a-n) - gamln(a))

    def _entropy(self, a):
        return a - (a+1.0)*special.psi(a) + gamln(a)
invgamma = invgamma_gen(a=0.0, name='invgamma')


# scale is gamma from DATAPLOT and B from Regress
class invgauss_gen(rv_continuous):
    """An inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invgauss` is::

        invgauss.pdf(x, mu) = 1 / sqrt(2*pi*x**3) * exp(-(x-mu)**2/(2*x*mu**2))

    for ``x > 0``.

    When `mu` is too small, evaluating the cumulative density function will be
    inaccurate due to ``cdf(mu -> 0) = inf * 0``.
    NaNs are returned for ``mu <= 0.0028``.

    %(example)s

    """
    def _rvs(self, mu):
        return mtrand.wald(mu, 1.0, size=self._size)

    def _pdf(self, x, mu):
        return 1.0/sqrt(2*pi*x**3.0)*exp(-1.0/(2*x)*((x-mu)/mu)**2)

    def _logpdf(self, x, mu):
        return -0.5*log(2*pi) - 1.5*log(x) - ((x-mu)/mu)**2/(2*x)

    def _cdf(self, x, mu):
        fac = sqrt(1.0/x)
        # Numerical accuracy for small `mu` is bad.  See #869.
        C1 = norm.cdf(fac*(x-mu)/mu)
        C1 += exp(1.0/mu) * norm.cdf(-fac*(x+mu)/mu) * exp(1.0/mu)
        return C1

    def _stats(self, mu):
        return mu, mu**3.0, 3*sqrt(mu), 15*mu
invgauss = invgauss_gen(a=0.0, name='invgauss')


class invweibull_gen(rv_continuous):
    """An inverted Weibull continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invweibull` is::

        invweibull.pdf(x, c) = c * x**(-c-1) * exp(-x**(-c))

    for ``x > 0``, ``c > 0``.

    References
    ----------
    F.R.S. de Gusmao, E.M.M Ortega and G.M. Cordeiro, "The generalized inverse
    Weibull distribution", Stat. Papers, vol. 52, pp. 591-619, 2011.

    %(example)s

    """
    def _pdf(self, x, c):
        xc1 = x**(-c-1.0)
        xc2 = x**(-c)
        xc2 = exp(-xc2)
        return c*xc1*xc2

    def _cdf(self, x, c):
        xc1 = x**(-c)
        return exp(-xc1)

    def _ppf(self, q, c):
        return pow(-log(q),asarray(-1.0/c))

    def _munp(self, n, c):
        return special.gamma(1 - n / c)

    def _entropy(self, c):
        return 1+_EULER + _EULER / c - log(c)
invweibull = invweibull_gen(a=0, name='invweibull')


class johnsonsb_gen(rv_continuous):
    """A Johnson SB continuous random variable.

    %(before_notes)s

    See Also
    --------
    johnsonsu

    Notes
    -----
    The probability density function for `johnsonsb` is::

        johnsonsb.pdf(x, a, b) = b / (x*(1-x)) * phi(a + b * log(x/(1-x)))

    for ``0 < x < 1`` and ``a,b > 0``, and ``phi`` is the normal pdf.

    %(example)s

    """
    def _argcheck(self, a, b):
        return (b > 0) & (a == a)

    def _pdf(self, x, a, b):
        trm = norm.pdf(a+b*log(x/(1.0-x)))
        return b*1.0/(x*(1-x))*trm

    def _cdf(self, x, a, b):
        return norm.cdf(a+b*log(x/(1.0-x)))

    def _ppf(self, q, a, b):
        return 1.0/(1+exp(-1.0/b*(norm.ppf(q)-a)))
johnsonsb = johnsonsb_gen(a=0.0, b=1.0, name='johnsonb')


class johnsonsu_gen(rv_continuous):
    """A Johnson SU continuous random variable.

    %(before_notes)s

    See Also
    --------
    johnsonsb

    Notes
    -----
    The probability density function for `johnsonsu` is::

        johnsonsu.pdf(x, a, b) = b / sqrt(x**2 + 1) *
                                 phi(a + b * log(x + sqrt(x**2 + 1)))

    for all ``x, a, b > 0``, and `phi` is the normal pdf.

    %(example)s

    """
    def _argcheck(self, a, b):
        return (b > 0) & (a == a)

    def _pdf(self, x, a, b):
        x2 = x*x
        trm = norm.pdf(a+b*log(x+sqrt(x2+1)))
        return b*1.0/sqrt(x2+1.0)*trm

    def _cdf(self, x, a, b):
        return norm.cdf(a+b*log(x+sqrt(x*x+1)))

    def _ppf(self, q, a, b):
        return sinh((norm.ppf(q)-a)/b)
johnsonsu = johnsonsu_gen(name='johnsonsu')


class laplace_gen(rv_continuous):
    """A Laplace continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `laplace` is::

        laplace.pdf(x) = 1/2 * exp(-abs(x))

    %(example)s

    """
    def _rvs(self):
        return mtrand.laplace(0, 1, size=self._size)

    def _pdf(self, x):
        return 0.5*exp(-abs(x))

    def _cdf(self, x):
        return where(x > 0, 1.0-0.5*exp(-x), 0.5*exp(x))

    def _ppf(self, q):
        return where(q > 0.5, -log(2*(1-q)), log(2*q))

    def _stats(self):
        return 0, 2, 0, 3

    def _entropy(self):
        return log(2)+1
laplace = laplace_gen(name='laplace')


class levy_gen(rv_continuous):
    """A Levy continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy_stable, levy_l

    Notes
    -----
    The probability density function for `levy` is::

        levy.pdf(x) = 1 / (x * sqrt(2*pi*x)) * exp(-1/(2*x))

    for ``x > 0``.

    This is the same as the Levy-stable distribution with a=1/2 and b=1.

    %(example)s

    """
    def _pdf(self, x):
        return 1/sqrt(2*pi*x)/x*exp(-1/(2*x))

    def _cdf(self, x):
        return 2*(1-norm._cdf(1/sqrt(x)))

    def _ppf(self, q):
        val = norm._ppf(1-q/2.0)
        return 1.0/(val*val)

    def _stats(self):
        return inf, inf, nan, nan
levy = levy_gen(a=0.0,name="levy")


class levy_l_gen(rv_continuous):
    """A left-skewed Levy continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy, levy_stable

    Notes
    -----
    The probability density function for `levy_l` is::

        levy_l.pdf(x) = 1 / (abs(x) * sqrt(2*pi*abs(x))) * exp(-1/(2*abs(x)))

    for ``x < 0``.

    This is the same as the Levy-stable distribution with a=1/2 and b=-1.

    %(example)s

    """
    def _pdf(self, x):
        ax = abs(x)
        return 1/sqrt(2*pi*ax)/ax*exp(-1/(2*ax))

    def _cdf(self, x):
        ax = abs(x)
        return 2*norm._cdf(1/sqrt(ax))-1

    def _ppf(self, q):
        val = norm._ppf((q+1.0)/2)
        return -1.0/(val*val)

    def _stats(self):
        return inf, inf, nan, nan
levy_l = levy_l_gen(b=0.0, name="levy_l")


class levy_stable_gen(rv_continuous):
    """A Levy-stable continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy, levy_l

    Notes
    -----
    Levy-stable distribution (only random variates available -- ignore other
    docs)

    %(example)s

    """
    def _rvs(self, alpha, beta):
        sz = self._size
        TH = uniform.rvs(loc=-pi/2.0,scale=pi,size=sz)
        W = expon.rvs(size=sz)
        if alpha == 1:
            return 2/pi*(pi/2+beta*TH)*tan(TH)-beta*log((pi/2*W*cos(TH))/(pi/2+beta*TH))

        ialpha = 1.0/alpha
        aTH = alpha*TH
        if beta == 0:
            return W/(cos(TH)/tan(aTH)+sin(TH))*((cos(aTH)+sin(aTH)*tan(TH))/W)**ialpha

        val0 = beta*tan(pi*alpha/2)
        th0 = arctan(val0)/alpha
        val3 = W/(cos(TH)/tan(alpha*(th0+TH))+sin(TH))
        res3 = val3*((cos(aTH)+sin(aTH)*tan(TH)-val0*(sin(aTH)-cos(aTH)*tan(TH)))/W)**ialpha
        return res3

    def _argcheck(self, alpha, beta):
        if beta == -1:
            self.b = 0.0
        elif beta == 1:
            self.a = 0.0
        return (alpha > 0) & (alpha <= 2) & (beta <= 1) & (beta >= -1)

    def _pdf(self, x, alpha, beta):
        raise NotImplementedError
levy_stable = levy_stable_gen(name='levy_stable')


class logistic_gen(rv_continuous):
    """A logistic (or Sech-squared) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `logistic` is::

        logistic.pdf(x) = exp(-x) / (1+exp(-x))**2

    `logistic` is a special case of `genlogistic` with ``c == 1``.

    %(example)s

    """
    def _rvs(self):
        return mtrand.logistic(size=self._size)

    def _pdf(self, x):
        ex = exp(-x)
        return ex / (1+ex)**2.0

    def _cdf(self, x):
        return 1.0/(1+exp(-x))

    def _ppf(self, q):
        return -log(1.0/q-1)

    def _stats(self):
        return 0, pi*pi/3.0, 0, 6.0/5.0

    def _entropy(self):
        return 1.0
logistic = logistic_gen(name='logistic')


class loggamma_gen(rv_continuous):
    """A log gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `loggamma` is::

        loggamma.pdf(x, c) = exp(c*x-exp(x)) / gamma(c)

    for all ``x, c > 0``.

    %(example)s

    """
    def _rvs(self, c):
        return log(mtrand.gamma(c, size=self._size))

    def _pdf(self, x, c):
        return exp(c*x-exp(x)-gamln(c))

    def _cdf(self, x, c):
        return special.gammainc(c, exp(x))

    def _ppf(self, q, c):
        return log(special.gammaincinv(c,q))

    def _stats(self, c):
        # See, for example, "A Statistical Study of Log-Gamma Distribution", by
        # Ping Shing Chan (thesis, McMaster University, 1993).
        mean = special.digamma(c)
        var = special.polygamma(1, c)
        skewness = special.polygamma(2, c) / np.power(var, 1.5)
        excess_kurtosis = special.polygamma(3, c) / (var*var)
        return mean, var, skewness, excess_kurtosis

loggamma = loggamma_gen(name='loggamma')


class loglaplace_gen(rv_continuous):
    """A log-Laplace continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `loglaplace` is::

    loglaplace.pdf(x, c) = c / 2 * x**(c-1),   for 0 < x < 1
                         = c / 2 * x**(-c-1),  for x >= 1

    for ``c > 0``.

    References
    ----------
    T.J. Kozubowski and K. Podgorski, "A log-Laplace growth rate model",
    The Mathematical Scientist, vol. 28, pp. 49-60, 2003.

    %(example)s

    """
    def _pdf(self, x, c):
        cd2 = c/2.0
        c = where(x < 1, c, -c)
        return cd2*x**(c-1)

    def _cdf(self, x, c):
        return where(x < 1, 0.5*x**c, 1-0.5*x**(-c))

    def _ppf(self, q, c):
        return where(q < 0.5, (2.0*q)**(1.0/c), (2*(1.0-q))**(-1.0/c))

    def _munp(self, n, c):
        return c**2 / (c**2 - n**2)

    def _entropy(self, c):
        return log(2.0/c) + 1.0
loglaplace = loglaplace_gen(a=0.0, name='loglaplace')


def _lognorm_logpdf(x, s):
    return -log(x)**2 / (2*s**2) + np.where(x == 0, 0, -log(s*x*sqrt(2*pi)))


class lognorm_gen(rv_continuous):
    """A lognormal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `lognorm` is::

        lognorm.pdf(x, s) = 1 / (s*x*sqrt(2*pi)) * exp(-1/2*(log(x)/s)**2)

    for ``x > 0``, ``s > 0``.

    If ``log(x)`` is normally distributed with mean ``mu`` and variance ``sigma**2``,
    then ``x`` is log-normally distributed with shape parameter sigma and scale
    parameter ``exp(mu)``.

    %(example)s

    """
    def _rvs(self, s):
        return exp(s * mtrand.standard_normal(self._size))

    def _pdf(self, x, s):
        return exp(self._logpdf(x, s))

    def _logpdf(self, x, s):
        return _lognorm_logpdf(x, s)

    def _cdf(self, x, s):
        return _norm_cdf(log(x) / s)

    def _ppf(self, q, s):
        return exp(s * _norm_ppf(q))

    def _stats(self, s):
        p = exp(s*s)
        mu = sqrt(p)
        mu2 = p*(p-1)
        g1 = sqrt((p-1))*(2+p)
        g2 = numpy.polyval([1,2,3,0,-6.0],p)
        return mu, mu2, g1, g2

    def _entropy(self, s):
        return 0.5 * (1 + log(2*pi) + 2 * log(s))
lognorm = lognorm_gen(a=0.0, name='lognorm')


class gilbrat_gen(rv_continuous):
    """A Gilbrat continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gilbrat` is::

        gilbrat.pdf(x) = 1/(x*sqrt(2*pi)) * exp(-1/2*(log(x))**2)

    `gilbrat` is a special case of `lognorm` with ``s = 1``.

    %(example)s

    """
    def _rvs(self):
        return exp(mtrand.standard_normal(self._size))

    def _pdf(self, x):
        return exp(self._logpdf(x))

    def _logpdf(self, x):
        return _lognorm_logpdf(x, 1.0)

    def _cdf(self, x):
        return _norm_cdf(log(x))

    def _ppf(self, q):
        return exp(_norm_ppf(q))

    def _stats(self):
        p = np.e
        mu = sqrt(p)
        mu2 = p * (p - 1)
        g1 = sqrt((p - 1)) * (2 + p)
        g2 = numpy.polyval([1, 2, 3, 0, -6.0], p)
        return mu, mu2, g1, g2

    def _entropy(self):
        return 0.5 * log(2 * pi) + 0.5
gilbrat = gilbrat_gen(a=0.0, name='gilbrat')


class maxwell_gen(rv_continuous):
    """A Maxwell continuous random variable.

    %(before_notes)s

    Notes
    -----
    A special case of a `chi` distribution,  with ``df = 3``, ``loc = 0.0``,
    and given ``scale = a``, where ``a`` is the parameter used in the
    Mathworld description [1]_.

    The probability density function for `maxwell` is::

        maxwell.pdf(x) = sqrt(2/pi)x**2 * exp(-x**2/2)

    for ``x > 0``.

    References
    ----------
    .. [1] http://mathworld.wolfram.com/MaxwellDistribution.html

    %(example)s
    """
    def _rvs(self):
        return chi.rvs(3.0,size=self._size)

    def _pdf(self, x):
        return sqrt(2.0/pi)*x*x*exp(-x*x/2.0)

    def _cdf(self, x):
        return special.gammainc(1.5,x*x/2.0)

    def _ppf(self, q):
        return sqrt(2*special.gammaincinv(1.5,q))

    def _stats(self):
        val = 3*pi-8
        return 2*sqrt(2.0/pi), 3-8/pi, sqrt(2)*(32-10*pi)/val**1.5, \
               (-12*pi*pi + 160*pi - 384) / val**2.0

    def _entropy(self):
        return _EULER + 0.5*log(2*pi)-0.5
maxwell = maxwell_gen(a=0.0, name='maxwell')


class mielke_gen(rv_continuous):
    """A Mielke's Beta-Kappa continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `mielke` is::

        mielke.pdf(x, k, s) = k * x**(k-1) / (1+x**s)**(1+k/s)

    for ``x > 0``.

    %(example)s

    """
    def _pdf(self, x, k, s):
        return k*x**(k-1.0) / (1.0+x**s)**(1.0+k*1.0/s)

    def _cdf(self, x, k, s):
        return x**k / (1.0+x**s)**(k*1.0/s)

    def _ppf(self, q, k, s):
        qsk = pow(q,s*1.0/k)
        return pow(qsk/(1.0-qsk),1.0/s)
mielke = mielke_gen(a=0.0, name='mielke')


class nakagami_gen(rv_continuous):
    """A Nakagami continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `nakagami` is::

        nakagami.pdf(x, nu) = 2 * nu**nu / gamma(nu) *
                              x**(2*nu-1) * exp(-nu*x**2)

    for ``x > 0``, ``nu > 0``.

    %(example)s

    """
    def _pdf(self, x, nu):
        return 2*nu**nu/gam(nu)*(x**(2*nu-1.0))*exp(-nu*x*x)

    def _cdf(self, x, nu):
        return special.gammainc(nu,nu*x*x)

    def _ppf(self, q, nu):
        return sqrt(1.0/nu*special.gammaincinv(nu,q))

    def _stats(self, nu):
        mu = gam(nu+0.5)/gam(nu)/sqrt(nu)
        mu2 = 1.0-mu*mu
        g1 = mu * (1 - 4*nu*mu2) / 2.0 / nu / np.power(mu2, 1.5)
        g2 = -6*mu**4*nu + (8*nu-2)*mu**2-2*nu + 1
        g2 /= nu*mu2**2.0
        return mu, mu2, g1, g2
nakagami = nakagami_gen(a=0.0, name="nakagami")


class ncx2_gen(rv_continuous):
    """A non-central chi-squared continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `ncx2` is::

        ncx2.pdf(x, df, nc) = exp(-(nc+df)/2) * 1/2 * (x/nc)**((df-2)/4)
                              * I[(df-2)/2](sqrt(nc*x))

    for ``x > 0``.

    %(example)s

    """
    def _rvs(self, df, nc):
        return mtrand.noncentral_chisquare(df,nc,self._size)

    def _logpdf(self, x, df, nc):
        a = asarray(df/2.0)
        fac = -nc/2.0 - x/2.0 + (a-1)*np.log(x) - a*np.log(2) - special.gammaln(a)
        return fac + np.nan_to_num(np.log(special.hyp0f1(a, nc * x/4.0)))

    def _pdf(self, x, df, nc):
        return np.exp(self._logpdf(x, df, nc))

    def _cdf(self, x, df, nc):
        return special.chndtr(x,df,nc)

    def _ppf(self, q, df, nc):
        return special.chndtrix(q,df,nc)

    def _stats(self, df, nc):
        val = df + 2.0*nc
        return df + nc, 2*val, sqrt(8)*(val+nc)/val**1.5, \
               12.0*(val+2*nc)/val**2.0
ncx2 = ncx2_gen(a=0.0, name='ncx2')


class ncf_gen(rv_continuous):
    """A non-central F distribution continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `ncf` is::

    ncf.pdf(x, df1, df2, nc) = exp(nc/2 + nc*df1*x/(2*(df1*x+df2)))
                    * df1**(df1/2) * df2**(df2/2) * x**(df1/2-1)
                    * (df2+df1*x)**(-(df1+df2)/2)
                    * gamma(df1/2)*gamma(1+df2/2)
                    * L^{v1/2-1}^{v2/2}(-nc*v1*x/(2*(v1*x+v2)))
                    / (B(v1/2, v2/2) * gamma((v1+v2)/2))

    for ``df1, df2, nc > 0``.

    %(example)s

    """
    def _rvs(self, dfn, dfd, nc):
        return mtrand.noncentral_f(dfn,dfd,nc,self._size)

    def _pdf_skip(self, x, dfn, dfd, nc):
        n1,n2 = dfn, dfd
        term = -nc/2+nc*n1*x/(2*(n2+n1*x)) + gamln(n1/2.)+gamln(1+n2/2.)
        term -= gamln((n1+n2)/2.0)
        Px = exp(term)
        Px *= n1**(n1/2) * n2**(n2/2) * x**(n1/2-1)
        Px *= (n2+n1*x)**(-(n1+n2)/2)
        Px *= special.assoc_laguerre(-nc*n1*x/(2.0*(n2+n1*x)),n2/2,n1/2-1)
        Px /= special.beta(n1/2,n2/2)
         # this function does not have a return
         #   drop it for now, the generic function seems to work ok

    def _cdf(self, x, dfn, dfd, nc):
        return special.ncfdtr(dfn,dfd,nc,x)

    def _ppf(self, q, dfn, dfd, nc):
        return special.ncfdtri(dfn, dfd, nc, q)

    def _munp(self, n, dfn, dfd, nc):
        val = (dfn * 1.0/dfd)**n
        term = gamln(n+0.5*dfn) + gamln(0.5*dfd-n) - gamln(dfd*0.5)
        val *= exp(-nc / 2.0+term)
        val *= special.hyp1f1(n+0.5*dfn, 0.5*dfn, 0.5*nc)
        return val

    def _stats(self, dfn, dfd, nc):
        mu = where(dfd <= 2, inf, dfd / (dfd-2.0)*(1+nc*1.0/dfn))
        mu2 = where(dfd <= 4, inf, 2*(dfd*1.0/dfn)**2.0 *
                    ((dfn+nc/2.0)**2.0 + (dfn+nc)*(dfd-2.0)) /
                    ((dfd-2.0)**2.0 * (dfd-4.0)))
        return mu, mu2, None, None
ncf = ncf_gen(a=0.0, name='ncf')


class t_gen(rv_continuous):
    """A Student's T continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `t` is::

                                       gamma((df+1)/2)
        t.pdf(x, df) = ---------------------------------------------------
                       sqrt(pi*df) * gamma(df/2) * (1+x**2/df)**((df+1)/2)

    for ``df > 0``.

    %(example)s

    """
    def _rvs(self, df):
        return mtrand.standard_t(df, size=self._size)

    def _pdf(self, x, df):
        r = asarray(df*1.0)
        Px = exp(gamln((r+1)/2)-gamln(r/2))
        Px /= sqrt(r*pi)*(1+(x**2)/r)**((r+1)/2)
        return Px

    def _logpdf(self, x, df):
        r = df*1.0
        lPx = gamln((r+1)/2)-gamln(r/2)
        lPx -= 0.5*log(r*pi) + (r+1)/2*log(1+(x**2)/r)
        return lPx

    def _cdf(self, x, df):
        return special.stdtr(df, x)

    def _sf(self, x, df):
        return special.stdtr(df, -x)

    def _ppf(self, q, df):
        return special.stdtrit(df, q)

    def _isf(self, q, df):
        return -special.stdtrit(df, q)

    def _stats(self, df):
        mu2 = where(df > 2, df / (df-2.0), inf)
        g1 = where(df > 3, 0.0, nan)
        g2 = where(df > 4, 6.0/(df-4.0), nan)
        return 0, mu2, g1, g2
t = t_gen(name='t')


class nct_gen(rv_continuous):
    """A non-central Student's T continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `nct` is::

                                            df**(df/2) * gamma(df+1)
        nct.pdf(x, df, nc) = ----------------------------------------------------
                             2**df*exp(nc**2/2) * (df+x**2)**(df/2) * gamma(df/2)

    for ``df > 0``.

    %(example)s

    """
    def _argcheck(self, df, nc):
        return (df > 0) & (nc == nc)

    def _rvs(self, df, nc):
        return norm.rvs(loc=nc,size=self._size)*sqrt(df) / sqrt(chi2.rvs(df,size=self._size))

    def _pdf(self, x, df, nc):
        n = df*1.0
        nc = nc*1.0
        x2 = x*x
        ncx2 = nc*nc*x2
        fac1 = n + x2
        trm1 = n/2.*log(n) + gamln(n+1)
        trm1 -= n*log(2)+nc*nc/2.+(n/2.)*log(fac1)+gamln(n/2.)
        Px = exp(trm1)
        valF = ncx2 / (2*fac1)
        trm1 = sqrt(2)*nc*x*special.hyp1f1(n/2+1,1.5,valF)
        trm1 /= asarray(fac1*special.gamma((n+1)/2))
        trm2 = special.hyp1f1((n+1)/2,0.5,valF)
        trm2 /= asarray(sqrt(fac1)*special.gamma(n/2+1))
        Px *= trm1+trm2
        return Px

    def _cdf(self, x, df, nc):
        return special.nctdtr(df, nc, x)

    def _ppf(self, q, df, nc):
        return special.nctdtrit(df, nc, q)

    def _stats(self, df, nc, moments='mv'):
        mu, mu2, g1, g2 = None, None, None, None
        val1 = gam((df-1.0)/2.0)
        val2 = gam(df/2.0)
        if 'm' in moments:
            mu = nc*sqrt(df/2.0)*val1/val2
        if 'v' in moments:
            var = (nc*nc+1.0)*df/(df-2.0)
            var -= nc*nc*df * val1**2 / 2.0 / val2**2
            mu2 = var
        if 's' in moments:
            g1n = 2*nc*sqrt(df)*val1*((nc*nc*(2*df-7)-3)*val2**2
                                      - nc*nc*(df-2)*(df-3)*val1**2)
            g1d = (df-3)*sqrt(2*df*(nc*nc+1)/(df-2) -
                              nc*nc*df*(val1/val2)**2) * val2 * \
                              (nc*nc*(df-2)*val1**2 -
                               2*(nc*nc+1)*val2**2)
            g1 = g1n/g1d
        if 'k' in moments:
            g2n = 2*(-3*nc**4*(df-2)**2 * (df-3) * (df-4)*val1**4 +
                     2**(6-2*df) * nc*nc*(df-2)*(df-4) *
                     (nc*nc*(2*df-7)-3)*pi*gam(df+1)**2 -
                     4*(nc**4*(df-5)-6*nc*nc-3)*(df-3)*val2**4)
            g2d = (df-3)*(df-4)*(nc*nc*(df-2)*val1**2 -
                                 2*(nc*nc+1)*val2)**2
            g2 = g2n / g2d
        return mu, mu2, g1, g2
nct = nct_gen(name="nct")


class pareto_gen(rv_continuous):
    """A Pareto continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `pareto` is::

        pareto.pdf(x, b) = b / x**(b+1)

    for ``x >= 1``, ``b > 0``.

    %(example)s

    """
    def _pdf(self, x, b):
        return b * x**(-b-1)

    def _cdf(self, x, b):
        return 1 - x**(-b)

    def _ppf(self, q, b):
        return pow(1-q, -1.0/b)

    def _stats(self, b, moments='mv'):
        mu, mu2, g1, g2 = None, None, None, None
        if 'm' in moments:
            mask = b > 1
            bt = extract(mask,b)
            mu = valarray(shape(b),value=inf)
            place(mu, mask, bt / (bt-1.0))
        if 'v' in moments:
            mask = b > 2
            bt = extract(mask,b)
            mu2 = valarray(shape(b), value=inf)
            place(mu2, mask, bt / (bt-2.0) / (bt-1.0)**2)
        if 's' in moments:
            mask = b > 3
            bt = extract(mask,b)
            g1 = valarray(shape(b), value=nan)
            vals = 2 * (bt + 1.0) * sqrt(bt - 2.0) / ((bt - 3.0) * sqrt(bt))
            place(g1, mask, vals)
        if 'k' in moments:
            mask = b > 4
            bt = extract(mask,b)
            g2 = valarray(shape(b), value=nan)
            vals = 6.0*polyval([1.0,1.0,-6,-2],bt) / \
                   polyval([1.0,-7.0,12.0,0.0],bt)
            place(g2, mask, vals)
        return mu, mu2, g1, g2

    def _entropy(self, c):
        return 1 + 1.0/c - log(c)
pareto = pareto_gen(a=1.0, name="pareto")


class lomax_gen(rv_continuous):
    """A Lomax (Pareto of the second kind) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The Lomax distribution is a special case of the Pareto distribution, with
    (loc=-1.0).

    The probability density function for `lomax` is::

        lomax.pdf(x, c) = c / (1+x)**(c+1)

    for ``x >= 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        return c*1.0/(1.0+x)**(c+1.0)

    def _logpdf(self, x, c):
        return log(c) - (c+1)*log(1+x)

    def _cdf(self, x, c):
        return 1.0-1.0/(1.0+x)**c

    def _sf(self, x, c):
        return 1.0/(1.0+x)**c

    def _logsf(self, x, c):
        return -c*log(1+x)

    def _ppf(self, q, c):
        return pow(1.0-q,-1.0/c)-1

    def _stats(self, c):
        mu, mu2, g1, g2 = pareto.stats(c, loc=-1.0, moments='mvsk')
        return mu, mu2, g1, g2

    def _entropy(self, c):
        return 1+1.0/c-log(c)
lomax = lomax_gen(a=0.0, name="lomax")


class pearson3_gen(rv_continuous):
    """A pearson type III continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `pearson3` is::

        pearson3.pdf(x, skew) = abs(beta) / gamma(alpha) *
            (beta * (x - zeta))**(alpha - 1) * exp(-beta*(x - zeta))

    where::

            beta = 2 / (skew * stddev)
            alpha = (stddev * beta)**2
            zeta = loc - alpha / beta

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

        ans, x, skew = np.broadcast_arrays([1.0], x, skew)
        ans = ans.copy()

        mask = np.absolute(skew) < norm2pearson_transition
        invmask = ~mask

        beta = 2.0 / (skew[invmask] * scale)
        alpha = (scale * beta)**2
        zeta = loc - alpha / beta

        transx = beta * (x[invmask] - zeta)
        return ans, x, transx, skew, mask, invmask, beta, alpha, zeta

    def _argcheck(self, skew):
        # The _argcheck function in rv_continuous only allows positive
        # arguments.  The skew argument for pearson3 can be zero (which I want
        # to handle inside pearson3._pdf) or negative.  So just return True
        # for all skew args.
        return np.ones(np.shape(skew), dtype=bool)

    def _stats(self, skew):
        ans, x, transx, skew, mask, invmask, beta, alpha, zeta = self._preprocess([1], skew)
        m = zeta + alpha / beta
        v = alpha / (beta**2)
        s = 2.0 / (alpha**0.5) * np.sign(beta)
        k = 6.0 / alpha
        return m, v, s, k

    def _pdf(self, x, skew):
        # Do the calculation in _logpdf since helps to limit
        # overflow/underflow problems
        ans = exp(self._logpdf(x, skew))
        if ans.ndim == 0:
            if np.isnan(ans):
                return 0.0
            return ans
        ans[np.isnan(ans)] = 0.0
        return ans

    def _logpdf(self, x, skew):
        #   PEARSON3 logpdf                           GAMMA logpdf
        #   np.log(abs(beta))
        # + (alpha - 1)*log(beta*(x - zeta))          + (a - 1)*log(x)
        # - beta*(x - zeta)                           - x
        # - gamln(alpha)                              - gamln(a)
        ans, x, transx, skew, mask, invmask, beta, alpha, zeta = self._preprocess(x, skew)

        ans[mask] = np.log(_norm_pdf(x[mask]))
        ans[invmask] = log(abs(beta)) + gamma._logpdf(transx, alpha)
        return ans

    def _cdf(self, x, skew):
        ans, x, transx, skew, mask, invmask, beta, alpha, zeta = self._preprocess(x, skew)

        ans[mask] = _norm_cdf(x[mask])
        ans[invmask] = gamma._cdf(transx, alpha)
        return ans

    def _rvs(self, skew):
        ans, x, transx, skew, mask, invmask, beta, alpha, zeta = self._preprocess([0], skew)
        if mask[0]:
            return mtrand.standard_normal(self._size)
        ans = mtrand.standard_gamma(alpha, self._size)/beta + zeta
        if ans.size == 1:
            return ans[0]
        return ans

    def _ppf(self, q, skew):
        ans, q, transq, skew, mask, invmask, beta, alpha, zeta = self._preprocess(q, skew)
        ans[mask] = _norm_ppf(q[mask])
        ans[invmask] = special.gammaincinv(alpha,q[invmask])/beta + zeta
        return ans
pearson3 = pearson3_gen(name="pearson3")


class powerlaw_gen(rv_continuous):
    """A power-function continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powerlaw` is::

        powerlaw.pdf(x, a) = a * x**(a-1)

    for ``0 <= x <= 1``, ``a > 0``.

    `powerlaw` is a special case of `beta` with ``d == 1``.

    %(example)s

    """
    def _pdf(self, x, a):
        return a*x**(a-1.0)

    def _logpdf(self, x, a):
        return log(a) + (a-1)*log(x)

    def _cdf(self, x, a):
        return x**(a*1.0)

    def _logcdf(self, x, a):
        return a*log(x)

    def _ppf(self, q, a):
        return pow(q, 1.0/a)

    def _stats(self, a):
        return (a / (a + 1.0),
                a / (a + 2.0) / (a + 1.0) ** 2,
                -2.0 * ((a - 1.0) / (a + 3.0)) * sqrt((a + 2.0) / a),
                6 * polyval([1, -1, -6, 2], a) / (a * (a + 3.0) * (a + 4)))

    def _entropy(self, a):
        return 1 - 1.0/a - log(a)
powerlaw = powerlaw_gen(a=0.0, b=1.0, name="powerlaw")


class powerlognorm_gen(rv_continuous):
    """A power log-normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powerlognorm` is::

        powerlognorm.pdf(x, c, s) = c / (x*s) * phi(log(x)/s) *
                                                (Phi(-log(x)/s))**(c-1),

    where ``phi`` is the normal pdf, and ``Phi`` is the normal cdf,
    and ``x > 0``, ``s, c > 0``.

    %(example)s

    """
    def _pdf(self, x, c, s):
        return c/(x*s)*norm.pdf(log(x)/s)*pow(norm.cdf(-log(x)/s),c*1.0-1.0)

    def _cdf(self, x, c, s):
        return 1.0 - pow(norm.cdf(-log(x)/s),c*1.0)

    def _ppf(self, q, c, s):
        return exp(-s*norm.ppf(pow(1.0-q,1.0/c)))
powerlognorm = powerlognorm_gen(a=0.0, name="powerlognorm")


class powernorm_gen(rv_continuous):
    """A power normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powernorm` is::

        powernorm.pdf(x, c) = c * phi(x) * (Phi(-x))**(c-1)

    where ``phi`` is the normal pdf, and ``Phi`` is the normal cdf,
    and ``x > 0``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        return c*_norm_pdf(x) * \
               (_norm_cdf(-x)**(c-1.0))

    def _logpdf(self, x, c):
        return log(c) + _norm_logpdf(x) + (c-1)*_norm_logcdf(-x)

    def _cdf(self, x, c):
        return 1.0-_norm_cdf(-x)**(c*1.0)

    def _ppf(self, q, c):
        return -norm.ppf(pow(1.0-q,1.0/c))
powernorm = powernorm_gen(name='powernorm')


class rdist_gen(rv_continuous):
    """An R-distributed continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rdist` is::

        rdist.pdf(x, c) = (1-x**2)**(c/2-1) / B(1/2, c/2)

    for ``-1 <= x <= 1``, ``c > 0``.

    %(example)s

    """
    def _pdf(self, x, c):
        return np.power((1.0 - x**2), c / 2.0 - 1) / special.beta(0.5, c / 2.0)

    def _cdf(self, x, c):
        term1 = x / special.beta(0.5, c / 2.0)
        res = 0.5 + term1 * special.hyp2f1(0.5, 1 - c / 2.0, 1.5, x**2)
        # There's an issue with hyp2f1, it returns nans near x = +-1, c > 100.
        # Use the generic implementation in that case.  See gh-1285 for
        # background.
        if any(np.isnan(res)):
            return rv_continuous._cdf(self, x, c)

        return res

    def _munp(self, n, c):
        return (1 - (n % 2)) * special.beta((n + 1.0) / 2, c / 2.0)
rdist = rdist_gen(a=-1.0, b=1.0, name="rdist")


class rayleigh_gen(rv_continuous):
    """A Rayleigh continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rayleigh` is::

        rayleigh.pdf(r) = r * exp(-r**2/2)

    for ``x >= 0``.

    `rayleigh` is a special case of `chi` with ``df == 2``.

    %(example)s

    """
    def _rvs(self):
        return chi.rvs(2, size=self._size)

    def _pdf(self, r):
        return r * exp(-0.5 * r**2)

    def _cdf(self, r):
        return 1 - exp(-0.5 * r**2)

    def _ppf(self, q):
        return sqrt(-2 * log(1 - q))

    def _stats(self):
        val = 4 - pi
        return np.sqrt(pi/2), val/2, 2*(pi-3)*sqrt(pi)/val**1.5, \
               6*pi/val-16/val**2

    def _entropy(self):
        return _EULER/2.0 + 1 - 0.5*log(2)
rayleigh = rayleigh_gen(a=0.0, name="rayleigh")


class reciprocal_gen(rv_continuous):
    """A reciprocal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `reciprocal` is::

        reciprocal.pdf(x, a, b) = 1 / (x*log(b/a))

    for ``a <= x <= b``, ``a, b > 0``.

    %(example)s

    """
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self.d = log(b*1.0 / a)
        return (a > 0) & (b > 0) & (b > a)

    def _pdf(self, x, a, b):
        return 1.0 / (x * self.d)

    def _logpdf(self, x, a, b):
        return -log(x) - log(self.d)

    def _cdf(self, x, a, b):
        return (log(x)-log(a)) / self.d

    def _ppf(self, q, a, b):
        return a*pow(b*1.0/a,q)

    def _munp(self, n, a, b):
        return 1.0/self.d / n * (pow(b*1.0,n) - pow(a*1.0,n))

    def _entropy(self,a,b):
        return 0.5*log(a*b)+log(log(b/a))
reciprocal = reciprocal_gen(name="reciprocal")


# FIXME: PPF does not work.
class rice_gen(rv_continuous):
    """A Rice continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rice` is::

        rice.pdf(x, b) = x * exp(-(x**2+b**2)/2) * I[0](x*b)

    for ``x > 0``, ``b > 0``.

    %(example)s

    """
    def _pdf(self, x, b):
        return x*exp(-(x*x+b*b)/2.0)*special.i0(x*b)

    def _logpdf(self, x, b):
        return log(x) - (x*x + b*b)/2.0 + log(special.i0(x*b))

    def _munp(self, n, b):
        nd2 = n/2.0
        n1 = 1+nd2
        b2 = b*b/2.0
        return 2.0**(nd2)*exp(-b2)*special.gamma(n1) * \
               special.hyp1f1(n1,1,b2)
rice = rice_gen(a=0.0, name="rice")


# FIXME: PPF does not work.
class recipinvgauss_gen(rv_continuous):
    """A reciprocal inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `recipinvgauss` is::

        recipinvgauss.pdf(x, mu) = 1/sqrt(2*pi*x) * exp(-(1-mu*x)**2/(2*x*mu**2))

    for ``x >= 0``.

    %(example)s

    """
    def _rvs(self, mu):
        return 1.0/mtrand.wald(mu, 1.0, size=self._size)

    def _pdf(self, x, mu):
        return 1.0/sqrt(2*pi*x)*exp(-(1-mu*x)**2.0 / (2*x*mu**2.0))

    def _logpdf(self, x, mu):
        return -(1-mu*x)**2.0 / (2*x*mu**2.0) - 0.5*log(2*pi*x)

    def _cdf(self, x, mu):
        trm1 = 1.0/mu - x
        trm2 = 1.0/mu + x
        isqx = 1.0/sqrt(x)
        return 1.0-_norm_cdf(isqx*trm1)-exp(2.0/mu)*_norm_cdf(-isqx*trm2)
recipinvgauss = recipinvgauss_gen(a=0.0, name='recipinvgauss')


class semicircular_gen(rv_continuous):
    """A semicircular continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `semicircular` is::

        semicircular.pdf(x) = 2/pi * sqrt(1-x**2)

    for ``-1 <= x <= 1``.

    %(example)s

    """
    def _pdf(self, x):
        return 2.0/pi*sqrt(1-x*x)

    def _cdf(self, x):
        return 0.5+1.0/pi*(x*sqrt(1-x*x) + arcsin(x))

    def _stats(self):
        return 0, 0.25, 0, -1.0

    def _entropy(self):
        return 0.64472988584940017414
semicircular = semicircular_gen(a=-1.0, b=1.0, name="semicircular")


class triang_gen(rv_continuous):
    """A triangular continuous random variable.

    %(before_notes)s

    Notes
    -----
    The triangular distribution can be represented with an up-sloping line from
    ``loc`` to ``(loc + c*scale)`` and then downsloping for ``(loc + c*scale)``
    to ``(loc+scale)``.

    The standard form is in the range [0, 1] with c the mode.
    The location parameter shifts the start to `loc`.
    The scale parameter changes the width from 1 to `scale`.

    %(example)s

    """
    def _rvs(self, c):
        return mtrand.triangular(0, c, 1, self._size)

    def _argcheck(self, c):
        return (c >= 0) & (c <= 1)

    def _pdf(self, x, c):
        return where(x < c, 2*x/c, 2*(1-x)/(1-c))

    def _cdf(self, x, c):
        return where(x < c, x*x/c, (x*x-2*x+c)/(c-1))

    def _ppf(self, q, c):
        return where(q < c, sqrt(c*q), 1-sqrt((1-c)*(1-q)))

    def _stats(self, c):
        return (c+1.0)/3.0, (1.0-c+c*c)/18, sqrt(2)*(2*c-1)*(c+1)*(c-2) / \
               (5 * np.power((1.0-c+c*c), 1.5)), -3.0/5.0

    def _entropy(self,c):
        return 0.5-log(2)
triang = triang_gen(a=0.0, b=1.0, name="triang")


class truncexpon_gen(rv_continuous):
    """A truncated exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `truncexpon` is::

        truncexpon.pdf(x, b) = exp(-x) / (1-exp(-b))

    for ``0 < x < b``.

    %(example)s

    """
    def _argcheck(self, b):
        self.b = b
        return (b > 0)

    def _pdf(self, x, b):
        return exp(-x)/(1-exp(-b))

    def _logpdf(self, x, b):
        return -x - log(1-exp(-b))

    def _cdf(self, x, b):
        return (1.0-exp(-x))/(1-exp(-b))

    def _ppf(self, q, b):
        return -log(1-q+q*exp(-b))

    def _munp(self, n, b):
        # wrong answer with formula, same as in continuous.pdf
        # return gam(n+1)-special.gammainc(1+n,b)
        if n == 1:
            return (1-(b+1)*exp(-b))/(-expm1(-b))
        elif n == 2:
            return 2*(1-0.5*(b*b+2*b+2)*exp(-b))/(-expm1(-b))
        else:
            # return generic for higher moments
            # return rv_continuous._mom1_sc(self,n, b)
            return self._mom1_sc(n, b)

    def _entropy(self, b):
        eB = exp(b)
        return log(eB-1)+(1+eB*(b-1.0))/(1.0-eB)
truncexpon = truncexpon_gen(a=0.0, name='truncexpon')


class truncnorm_gen(rv_continuous):
    """A truncated normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The standard form of this distribution is a standard normal truncated to
    the range [a,b] --- notice that a and b are defined over the domain of the
    standard normal.  To convert clip values for a specific mean and standard
    deviation, use::

        a, b = (myclip_a - my_mean) / my_std, (myclip_b - my_mean) / my_std

    %(example)s

    """
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self._nb = _norm_cdf(b)
        self._na = _norm_cdf(a)
        self._sb = _norm_sf(b)
        self._sa = _norm_sf(a)
        if self.a > 0:
            self._delta = -(self._sb - self._sa)
        else:
            self._delta = self._nb - self._na
        self._logdelta = log(self._delta)
        return (a != b)

    def _pdf(self, x, a, b):
        return _norm_pdf(x) / self._delta

    def _logpdf(self, x, a, b):
        return _norm_logpdf(x) - self._logdelta

    def _cdf(self, x, a, b):
        return (_norm_cdf(x) - self._na) / self._delta

    def _ppf(self, q, a, b):
        if self.a > 0:
            return _norm_isf(q*self._sb + self._sa*(1.0-q))
        else:
            return _norm_ppf(q*self._nb + self._na*(1.0-q))

    def _stats(self, a, b):
        nA, nB = self._na, self._nb
        d = nB - nA
        pA, pB = _norm_pdf(a), _norm_pdf(b)
        mu = (pA - pB) / d   # correction sign
        mu2 = 1 + (a*pA - b*pB) / d - mu*mu
        return mu, mu2, None, None
truncnorm = truncnorm_gen(name='truncnorm')


# FIXME: RVS does not work.
class tukeylambda_gen(rv_continuous):
    """A Tukey-Lamdba continuous random variable.

    %(before_notes)s

    Notes
    -----
    A flexible distribution, able to represent and interpolate between the
    following distributions:

        - Cauchy                (lam=-1)
        - logistic              (lam=0.0)
        - approx Normal         (lam=0.14)
        - u-shape               (lam = 0.5)
        - uniform from -1 to 1  (lam = 1)

    %(example)s

    """
    def _argcheck(self, lam):
        return np.ones(np.shape(lam), dtype=bool)

    def _pdf(self, x, lam):
        Fx = asarray(special.tklmbda(x,lam))
        Px = Fx**(lam-1.0) + (asarray(1-Fx))**(lam-1.0)
        Px = 1.0/asarray(Px)
        return where((lam <= 0) | (abs(x) < 1.0/asarray(lam)), Px, 0.0)

    def _cdf(self, x, lam):
        return special.tklmbda(x, lam)

    def _ppf(self, q, lam):
        q = q*1.0
        vals1 = (q**lam - (1-q)**lam)/lam
        vals2 = log(q/(1-q))
        return where((lam == 0) & (q == q), vals2, vals1)

    def _stats(self, lam):
        return 0, _tlvar(lam), 0, _tlkurt(lam)

    def _entropy(self, lam):
        def integ(p):
            return log(pow(p,lam-1)+pow(1-p,lam-1))
        return integrate.quad(integ,0,1)[0]
tukeylambda = tukeylambda_gen(name='tukeylambda')


class uniform_gen(rv_continuous):
    """A uniform continuous random variable.

    This distribution is constant between `loc` and ``loc + scale``.

    %(before_notes)s

    %(example)s

    """
    def _rvs(self):
        return mtrand.uniform(0.0,1.0,self._size)

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
uniform = uniform_gen(a=0.0, b=1.0, name='uniform')


class vonmises_gen(rv_continuous):
    """A Von Mises continuous random variable.

    %(before_notes)s

    Notes
    -----
    If `x` is not in range or `loc` is not in range it assumes they are angles
    and converts them to [-pi, pi] equivalents.

    The probability density function for `vonmises` is::

        vonmises.pdf(x, kappa) = exp(kappa * cos(x)) / (2*pi*I[0](kappa))

    for ``-pi <= x <= pi``, ``kappa > 0``.

    %(example)s

    """
    def _rvs(self, kappa):
        return mtrand.vonmises(0.0, kappa, size=self._size)

    def _pdf(self, x, kappa):
        return exp(kappa * cos(x)) / (2*pi*special.i0(kappa))

    def _cdf(self, x, kappa):
        return vonmises_cython.von_mises_cdf(kappa, x)

    def _stats_skip(self, kappa):
        return 0, None, 0, None
vonmises = vonmises_gen(name='vonmises')


class wald_gen(invgauss_gen):
    """A Wald continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `wald` is::

        wald.pdf(x, a) = 1/sqrt(2*pi*x**3) * exp(-(x-1)**2/(2*x))

    for ``x > 0``.

    `wald` is a special case of `invgauss` with ``mu == 1``.

    %(example)s
    """
    def _rvs(self):
        return mtrand.wald(1.0, 1.0, size=self._size)

    def _pdf(self, x):
        return invgauss._pdf(x, 1.0)

    def _logpdf(self, x):
        return invgauss._logpdf(x, 1.0)

    def _cdf(self, x):
        return invgauss._cdf(x, 1.0)

    def _stats(self):
        return 1.0, 1.0, 3.0, 15.0
wald = wald_gen(a=0.0, name="wald")


class wrapcauchy_gen(rv_continuous):
    """A wrapped Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `wrapcauchy` is::

        wrapcauchy.pdf(x, c) = (1-c**2) / (2*pi*(1+c**2-2*c*cos(x)))

    for ``0 <= x <= 2*pi``, ``0 < c < 1``.

    %(example)s

    """
    def _argcheck(self, c):
        return (c > 0) & (c < 1)

    def _pdf(self, x, c):
        return (1.0-c*c)/(2*pi*(1+c*c-2*c*cos(x)))

    def _cdf(self, x, c):
        output = 0.0*x
        val = (1.0+c)/(1.0-c)
        c1 = x < pi
        c2 = 1-c1
        xp = extract(c1,x)
        xn = extract(c2,x)
        if (any(xn)):
            valn = extract(c2, np.ones_like(x)*val)
            xn = 2*pi - xn
            yn = tan(xn/2.0)
            on = 1.0-1.0/pi*arctan(valn*yn)
            place(output, c2, on)
        if (any(xp)):
            valp = extract(c1, np.ones_like(x)*val)
            yp = tan(xp/2.0)
            op = 1.0/pi*arctan(valp*yp)
            place(output, c1, op)
        return output

    def _ppf(self, q, c):
        val = (1.0-c)/(1.0+c)
        rcq = 2*arctan(val*tan(pi*q))
        rcmq = 2*pi-2*arctan(val*tan(pi*(1-q)))
        return where(q < 1.0/2, rcq, rcmq)

    def _entropy(self, c):
        return log(2*pi*(1-c*c))
wrapcauchy = wrapcauchy_gen(a=0.0, b=2*pi, name='wrapcauchy')


# DISCRETE DISTRIBUTIONS

def entropy(pk, qk=None, base=None):
    """Calculate the entropy of a distribution for given probability values.

    If only probabilities `pk` are given, the entropy is calculated as
    ``S = -sum(pk * log(pk), axis=0)``.

    If `qk` is not None, then compute a relative entropy (also known as
    Kullback-Leibler divergence or Kullback-Leibler distance)
    ``S = sum(pk * log(pk / qk), axis=0)``.

    This routine will normalize `pk` and `qk` if they don't sum to 1.

    Parameters
    ----------
    pk : sequence
        Defines the (discrete) distribution. ``pk[i]`` is the (possibly
        unnormalized) probability of event ``i``.
    qk : sequence, optional
        Sequence against which the relative entropy is computed. Should be in
        the same format as `pk`.
    base : float, optional
        The logarithmic base to use, defaults to ``e`` (natural logarithm).

    Returns
    -------
    S : float
        The calculated entropy.

    """
    pk = asarray(pk)
    pk = 1.0*pk / sum(pk, axis=0)
    if qk is None:
        vec = special.xlogy(pk, pk)
    else:
        qk = asarray(qk)
        if len(qk) != len(pk):
            raise ValueError("qk and pk must have same length.")
        qk = 1.0*qk / sum(qk, axis=0)
        # If qk is zero anywhere, then unless pk is zero at those places
        #   too, the relative entropy is infinite.
        if any(take(pk, nonzero(qk == 0.0), axis=0) != 0.0, 0):
            return inf
        vec = -special.xlogy(pk, pk / qk)
    S = -sum(vec, axis=0)
    if base is not None:
        S /= log(base)
    return S


## Handlers for generic case where xk and pk are given

def _drv_pmf(self, xk, *args):
    try:
        return self.P[xk]
    except KeyError:
        return 0.0


def _drv_cdf(self, xk, *args):
    indx = argmax((self.xk > xk),axis=-1)-1
    return self.F[self.xk[indx]]


def _drv_ppf(self, q, *args):
    indx = argmax((self.qvals >= q),axis=-1)
    return self.Finv[self.qvals[indx]]


def _drv_nonzero(self, k, *args):
    return 1


def _drv_moment(self, n, *args):
    n = asarray(n)
    return sum(self.xk**n[newaxis,...] * self.pk, axis=0)


def _drv_moment_gen(self, t, *args):
    t = asarray(t)
    return sum(exp(self.xk * t[newaxis,...]) * self.pk, axis=0)


def _drv2_moment(self, n, *args):
    """Non-central moment of discrete distribution."""
    # many changes, originally not even a return
    tot = 0.0
    diff = 1e100
    # pos = self.a
    pos = max(0.0, 1.0*self.a)
    count = 0
    # handle cases with infinite support
    ulimit = max(1000, (min(self.b,1000) + max(self.a,-1000))/2.0)
    llimit = min(-1000, (min(self.b,1000) + max(self.a,-1000))/2.0)

    while (pos <= self.b) and ((pos <= ulimit) or
                               (diff > self.moment_tol)):
        diff = np.power(pos, n) * self.pmf(pos,*args)
        # use pmf because _pmf does not check support in randint
        #     and there might be problems ? with correct self.a, self.b at this stage
        tot += diff
        pos += self.inc
        count += 1

    if self.a < 0:  # handle case when self.a = -inf
        diff = 1e100
        pos = -self.inc
        while (pos >= self.a) and ((pos >= llimit) or
                                   (diff > self.moment_tol)):
            diff = np.power(pos, n) * self.pmf(pos,*args)
            # using pmf instead of _pmf, see above
            tot += diff
            pos -= self.inc
            count += 1
    return tot


def _drv2_ppfsingle(self, q, *args):  # Use basic bisection algorithm
    b = self.b
    a = self.a
    if isinf(b):            # Be sure ending point is > q
        b = int(max(100*q,10))
        while 1:
            if b >= self.b:
                qb = 1.0
                break
            qb = self._cdf(b,*args)
            if (qb < q):
                b += 10
            else:
                break
    else:
        qb = 1.0
    if isinf(a):    # be sure starting point < q
        a = int(min(-100*q,-10))
        while 1:
            if a <= self.a:
                qb = 0.0
                break
            qa = self._cdf(a,*args)
            if (qa > q):
                a -= 10
            else:
                break
    else:
        qa = self._cdf(a, *args)

    while 1:
        if (qa == q):
            return a
        if (qb == q):
            return b
        if b <= a+1:
    # testcase: return wrong number at lower index
    # python -c "from scipy.stats import zipf;print zipf.ppf(0.01,2)" wrong
    # python -c "from scipy.stats import zipf;print zipf.ppf([0.01,0.61,0.77,0.83],2)"
    # python -c "from scipy.stats import logser;print logser.ppf([0.1,0.66, 0.86,0.93],0.6)"
            if qa > q:
                return a
            else:
                return b
        c = int((a+b)/2.0)
        qc = self._cdf(c, *args)
        if (qc < q):
            if a != c:
                a = c
            else:
                raise RuntimeError('updating stopped, endless loop')
            qa = qc
        elif (qc > q):
            if b != c:
                b = c
            else:
                raise RuntimeError('updating stopped, endless loop')
            qb = qc
        else:
            return c


def reverse_dict(dict):
    newdict = {}
    sorted_keys = list(dict.keys())
    sorted_keys.sort()
    for key in sorted_keys[::-1]:
        newdict[dict[key]] = key
    return newdict


def make_dict(keys, values):
    d = {}
    for key, value in zip(keys, values):
        d[key] = value
    return d


# Must over-ride one of _pmf or _cdf or pass in
#  x_k, p(x_k) lists in initialization

class rv_discrete(rv_generic):
    """
    A generic discrete random variable class meant for subclassing.

    `rv_discrete` is a base class to construct specific distribution classes
    and instances from for discrete random variables. rv_discrete can be used
    to construct an arbitrary distribution with defined by a list of support
    points and the corresponding probabilities.

    Parameters
    ----------
    a : float, optional
        Lower bound of the support of the distribution, default: 0
    b : float, optional
        Upper bound of the support of the distribution, default: plus infinity
    moment_tol : float, optional
        The tolerance for the generic calculation of moments
    values : tuple of two array_like
        (xk, pk) where xk are points (integers) with positive probability pk
        with sum(pk) = 1
    inc : integer
        increment for the support of the distribution, default: 1
        other values have not been tested
    badvalue : object, optional
        The value in (masked) arrays that indicates a value that should be
        ignored.
    name : str, optional
        The name of the instance. This string is used to construct the default
        example for distributions.
    longname : str, optional
        This string is used as part of the first line of the docstring returned
        when a subclass has no docstring of its own. Note: `longname` exists
        for backwards compatibility, do not use for new subclasses.
    shapes : str, optional
        The shape of the distribution. For example ``"m, n"`` for a
        distribution that takes two integers as the first two arguments for all
        its methods.
    extradoc :  str, optional
        This string is used as the last part of the docstring returned when a
        subclass has no docstring of its own. Note: `extradoc` exists for
        backwards compatibility, do not use for new subclasses.

    Methods
    -------
    generic.rvs(<shape(s)>, loc=0, size=1)
        random variates

    generic.pmf(x, <shape(s)>, loc=0)
        probability mass function

    logpmf(x, <shape(s)>, loc=0)
        log of the probability density function

    generic.cdf(x, <shape(s)>, loc=0)
        cumulative density function

    generic.logcdf(x, <shape(s)>, loc=0)
        log of the cumulative density function

    generic.sf(x, <shape(s)>, loc=0)
        survival function (1-cdf --- sometimes more accurate)

    generic.logsf(x, <shape(s)>, loc=0, scale=1)
        log of the survival function

    generic.ppf(q, <shape(s)>, loc=0)
        percent point function (inverse of cdf --- percentiles)

    generic.isf(q, <shape(s)>, loc=0)
        inverse survival function (inverse of sf)

    generic.moment(n, <shape(s)>, loc=0)
        non-central n-th moment of the distribution.  May not work for array arguments.

    generic.stats(<shape(s)>, loc=0, moments='mv')
        mean('m', axis=0), variance('v'), skew('s'), and/or kurtosis('k')

    generic.entropy(<shape(s)>, loc=0)
        entropy of the RV

    generic.expect(func=None, args=(), loc=0, lb=None, ub=None, conditional=False)
        Expected value of a function with respect to the distribution.
        Additional kwd arguments passed to integrate.quad

    generic.median(<shape(s)>, loc=0)
        Median of the distribution.

    generic.mean(<shape(s)>, loc=0)
        Mean of the distribution.

    generic.std(<shape(s)>, loc=0)
        Standard deviation of the distribution.

    generic.var(<shape(s)>, loc=0)
        Variance of the distribution.

    generic.interval(alpha, <shape(s)>, loc=0)
        Interval that with `alpha` percent probability contains a random
        realization of this distribution.

    generic(<shape(s)>, loc=0)
        calling a distribution instance returns a frozen distribution

    Notes
    -----

    You can construct an arbitrary discrete rv where ``P{X=xk} = pk``
    by passing to the rv_discrete initialization method (through the
    values=keyword) a tuple of sequences (xk, pk) which describes only those
    values of X (xk) that occur with nonzero probability (pk).

    To create a new discrete distribution, we would do the following::

        class poisson_gen(rv_discrete):
            #"Poisson distribution"
            def _pmf(self, k, mu):
                ...

    and create an instance::

        poisson = poisson_gen(name="poisson",
                              longname='A Poisson')

    The docstring can be created from a template.

    Alternatively, the object may be called (as a function) to fix the shape
    and location parameters returning a "frozen" discrete RV object::

        myrv = generic(<shape(s)>, loc=0)
            - frozen RV object with the same methods but holding the given
              shape and location fixed.

    A note on ``shapes``: subclasses need not specify them explicitly. In this
    case, the `shapes` will be automatically deduced from the signatures of the
    overridden methods.
    If, for some reason, you prefer to avoid relying on introspection, you can
    specify ``shapes`` explicitly as an argument to the instance constructor.


    Examples
    --------

    Custom made discrete distribution:

    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> xk = np.arange(7)
    >>> pk = (0.1, 0.2, 0.3, 0.1, 0.1, 0.1, 0.1)
    >>> custm = stats.rv_discrete(name='custm', values=(xk, pk))
    >>> h = plt.plot(xk, custm.pmf(xk))

    Random number generation:

    >>> R = custm.rvs(size=100)

    Display frozen pmf:

    >>> numargs = generic.numargs
    >>> [ <shape(s)> ] = ['Replace with resonable value', ]*numargs
    >>> rv = generic(<shape(s)>)
    >>> x = np.arange(0, np.min(rv.dist.b, 3)+1)
    >>> h = plt.plot(x, rv.pmf(x))

    Here, ``rv.dist.b`` is the right endpoint of the support of ``rv.dist``.

    Check accuracy of cdf and ppf:

    >>> prb = generic.cdf(x, <shape(s)>)
    >>> h = plt.semilogy(np.abs(x-generic.ppf(prb, <shape(s)>))+1e-20)

    """

    def __init__(self, a=0, b=inf, name=None, badvalue=None,
                 moment_tol=1e-8,values=None,inc=1,longname=None,
                 shapes=None, extradoc=None):

        super(rv_generic,self).__init__()

        if badvalue is None:
            badvalue = nan
        if name is None:
            name = 'Distribution'
        self.badvalue = badvalue
        self.a = a
        self.b = b
        self.name = name
        self.moment_tol = moment_tol
        self.inc = inc
        self._cdfvec = vectorize(self._cdfsingle, otypes='d')
        self.return_integers = 1
        self.vecentropy = vectorize(self._entropy)
        self.shapes = shapes
        self.extradoc = extradoc

        if values is not None:
            self.xk, self.pk = values
            self.return_integers = 0
            indx = argsort(ravel(self.xk))
            self.xk = take(ravel(self.xk),indx, 0)
            self.pk = take(ravel(self.pk),indx, 0)
            self.a = self.xk[0]
            self.b = self.xk[-1]
            self.P = make_dict(self.xk, self.pk)
            self.qvals = numpy.cumsum(self.pk,axis=0)
            self.F = make_dict(self.xk, self.qvals)
            self.Finv = reverse_dict(self.F)
            self._ppf = instancemethod(vectorize(_drv_ppf, otypes='d'),
                                       self, rv_discrete)
            self._pmf = instancemethod(vectorize(_drv_pmf, otypes='d'),
                                       self, rv_discrete)
            self._cdf = instancemethod(vectorize(_drv_cdf, otypes='d'),
                                       self, rv_discrete)
            self._nonzero = instancemethod(_drv_nonzero, self, rv_discrete)
            self.generic_moment = instancemethod(_drv_moment,
                                                 self, rv_discrete)
            self.moment_gen = instancemethod(_drv_moment_gen,
                                             self, rv_discrete)
            self._construct_argparser(names_to_inspect=['_drv_pmf'],
                                  locscale_in='loc=0',
                                  locscale_out='loc, 1')  # scale=1 for discrete RVs
        else:
            self._construct_argparser(names_to_inspect=['_pmf', '_cdf'],
                                  locscale_in='loc=0',
                                  locscale_out='loc, 1')  # scale=1 for discrete RVs

            # nin correction needs to be after we know numargs
            # correct nin for generic moment vectorization
            self.vec_generic_moment = vectorize(_drv2_moment, otypes='d')
            self.vec_generic_moment.nin = self.numargs + 2
            self.generic_moment = instancemethod(self.vec_generic_moment,
                                                 self, rv_discrete)

            # correct nin for ppf vectorization
            _vppf = vectorize(_drv2_ppfsingle, otypes='d')
            _vppf.nin = self.numargs + 2  # +1 is for self
            self._vecppf = instancemethod(_vppf,
                                          self, rv_discrete)

        # now that self.numargs is defined, we can adjust nin
        self._cdfvec.nin = self.numargs + 1

        # generate docstring for subclass instances
        if longname is None:
            if name[0] in ['aeiouAEIOU']:
                hstr = "An "
            else:
                hstr = "A "
            longname = hstr + name

        if sys.flags.optimize < 2:
            # Skip adding docstrings if interpreter is run with -OO
            if self.__doc__ is None:
                self._construct_default_doc(longname=longname, extradoc=extradoc)
            else:
                self._construct_doc()

            #discrete RV do not have the scale parameter, remove it
            self.__doc__ = self.__doc__.replace('\n    scale : array_like, '
                            'optional\n        scale parameter (default=1)', '')


    def _construct_default_doc(self, longname=None, extradoc=None):
        """Construct instance docstring from the rv_discrete template."""
        if extradoc is None:
            extradoc = ''
        if extradoc.startswith('\n\n'):
            extradoc = extradoc[2:]
        self.__doc__ = ''.join(['%s discrete random variable.' % longname,
                                '\n\n%(before_notes)s\n', docheaders['notes'],
                                extradoc, '\n%(example)s'])
        self._construct_doc()

    def _construct_doc(self):
        """Construct the instance docstring with string substitutions."""
        tempdict = docdict_discrete.copy()
        tempdict['name'] = self.name or 'distname'
        tempdict['shapes'] = self.shapes or ''

        if self.shapes is None:
            # remove shapes from call parameters if there are none
            for item in ['callparams', 'default', 'before_notes']:
                tempdict[item] = tempdict[item].replace(
                        "\n%(shapes)s : array_like\n    shape parameters", "")
        for i in range(2):
            if self.shapes is None:
                # necessary because we use %(shapes)s in two forms (w w/o ", ")
                self.__doc__ = self.__doc__.replace("%(shapes)s, ", "")
            self.__doc__ = doccer.docformat(self.__doc__, tempdict)

    def _rvs(self, *args):
        return self._ppf(mtrand.random_sample(self._size),*args)

    def _nonzero(self, k, *args):
        return floor(k) == k

    def _argcheck(self, *args):
        cond = 1
        for arg in args:
            cond &= (arg > 0)
        return cond

    def _pmf(self, k, *args):
        return self._cdf(k,*args) - self._cdf(k-1,*args)

    def _logpmf(self, k, *args):
        return log(self._pmf(k, *args))

    def _cdfsingle(self, k, *args):
        m = arange(int(self.a),k+1)
        return sum(self._pmf(m,*args),axis=0)

    def _cdf(self, x, *args):
        k = floor(x)
        return self._cdfvec(k,*args)

    def _logcdf(self, x, *args):
        return log(self._cdf(x, *args))

    def _sf(self, x, *args):
        return 1.0-self._cdf(x,*args)

    def _logsf(self, x, *args):
        return log(self._sf(x, *args))

    def _ppf(self, q, *args):
        return self._vecppf(q, *args)

    def _isf(self, q, *args):
        return self._ppf(1-q,*args)

    def _stats(self, *args):
        return None, None, None, None

    def _munp(self, n, *args):
        return self.generic_moment(n, *args)

    def rvs(self, *args, **kwargs):
        """
        Random variates of given type.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).
        size : int or tuple of ints, optional
            Defining number of random variates (default=1).  Note that `size`
            has to be given as keyword, not as positional argument.

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of given `size`.

        """
        kwargs['discrete'] = True
        return super(rv_discrete, self).rvs(*args, **kwargs)

    def pmf(self, k,*args, **kwds):
        """
        Probability mass function at k of the given RV.

        Parameters
        ----------
        k : array_like
            quantiles
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            Location parameter (default=0).

        Returns
        -------
        pmf : array_like
            Probability mass function evaluated at k

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k,loc = map(asarray,(k,loc))
        args = tuple(map(asarray,args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b) & self._nonzero(k,*args)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        place(output,(1-cond0) + np.isnan(k),self.badvalue)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output,cond,self._pmf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logpmf(self, k,*args, **kwds):
        """
        Log of the probability mass function at k of the given RV.

        Parameters
        ----------
        k : array_like
            Quantiles.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter. Default is 0.

        Returns
        -------
        logpmf : array_like
            Log of the probability mass function evaluated at k.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k,loc = map(asarray,(k,loc))
        args = tuple(map(asarray,args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b) & self._nonzero(k,*args)
        cond = cond0 & cond1
        output = empty(shape(cond),'d')
        output.fill(NINF)
        place(output,(1-cond0) + np.isnan(k),self.badvalue)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output,cond,self._logpmf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def cdf(self, k, *args, **kwds):
        """
        Cumulative distribution function of the given RV.

        Parameters
        ----------
        k : array_like, int
            Quantiles.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).

        Returns
        -------
        cdf : ndarray
            Cumulative distribution function evaluated at `k`.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k,loc = map(asarray,(k,loc))
        args = tuple(map(asarray,args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k < self.b)
        cond2 = (k >= self.b)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        place(output,(1-cond0) + np.isnan(k),self.badvalue)
        place(output,cond2*(cond0 == cond0), 1.0)

        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output,cond,self._cdf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logcdf(self, k, *args, **kwds):
        """
        Log of the cumulative distribution function at k of the given RV

        Parameters
        ----------
        k : array_like, int
            Quantiles.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).

        Returns
        -------
        logcdf : array_like
            Log of the cumulative distribution function evaluated at k.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k,loc = map(asarray,(k,loc))
        args = tuple(map(asarray,args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k < self.b)
        cond2 = (k >= self.b)
        cond = cond0 & cond1
        output = empty(shape(cond),'d')
        output.fill(NINF)
        place(output,(1-cond0) + np.isnan(k),self.badvalue)
        place(output,cond2*(cond0 == cond0), 0.0)

        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output,cond,self._logcdf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def sf(self,k,*args,**kwds):
        """
        Survival function (1-cdf) at k of the given RV.

        Parameters
        ----------
        k : array_like
            Quantiles.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).

        Returns
        -------
        sf : array_like
            Survival function evaluated at k.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k,loc = map(asarray,(k,loc))
        args = tuple(map(asarray,args))
        k = asarray(k-loc)
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b)
        cond2 = (k < self.a) & cond0
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        place(output,(1-cond0) + np.isnan(k),self.badvalue)
        place(output,cond2,1.0)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output,cond,self._sf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logsf(self,k,*args,**kwds):
        """
        Log of the survival function of the given RV.

        Returns the log of the "survival function," defined as ``1 - cdf``,
        evaluated at `k`.

        Parameters
        ----------
        k : array_like
            Quantiles.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).

        Returns
        -------
        sf : ndarray
            Survival function evaluated at `k`.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k,loc = map(asarray,(k,loc))
        args = tuple(map(asarray,args))
        k = asarray(k-loc)
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b)
        cond2 = (k < self.a) & cond0
        cond = cond0 & cond1
        output = empty(shape(cond),'d')
        output.fill(NINF)
        place(output,(1-cond0) + np.isnan(k),self.badvalue)
        place(output,cond2,0.0)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output,cond,self._logsf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def ppf(self,q,*args,**kwds):
        """
        Percent point function (inverse of cdf) at q of the given RV

        Parameters
        ----------
        q : array_like
            Lower tail probability.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).
        scale : array_like, optional
            Scale parameter (default=1).

        Returns
        -------
        k : array_like
            Quantile corresponding to the lower tail probability, q.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        q,loc = map(asarray,(q,loc))
        args = tuple(map(asarray,args))
        cond0 = self._argcheck(*args) & (loc == loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q == 1) & cond0
        cond = cond0 & cond1
        output = valarray(shape(cond),value=self.badvalue,typecode='d')
        # output type 'd' to handle nin and inf
        place(output,(q == 0)*(cond == cond), self.a-1)
        place(output,cond2,self.b)
        if any(cond):
            goodargs = argsreduce(cond, *((q,)+args+(loc,)))
            loc, goodargs = goodargs[-1], goodargs[:-1]
            place(output,cond,self._ppf(*goodargs) + loc)

        if output.ndim == 0:
            return output[()]
        return output

    def isf(self,q,*args,**kwds):
        """
        Inverse survival function (1-sf) at q of the given RV.

        Parameters
        ----------
        q : array_like
            Upper tail probability.
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).

        Returns
        -------
        k : ndarray or scalar
            Quantile corresponding to the upper tail probability, q.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        q,loc = map(asarray,(q,loc))
        args = tuple(map(asarray,args))
        cond0 = self._argcheck(*args) & (loc == loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q == 1) & cond0
        cond = cond0 & cond1

        # same problem as with ppf; copied from ppf and changed
        output = valarray(shape(cond),value=self.badvalue,typecode='d')
        # output type 'd' to handle nin and inf
        place(output,(q == 0)*(cond == cond), self.b)
        place(output,cond2,self.a-1)

        # call place only if at least 1 valid argument
        if any(cond):
            goodargs = argsreduce(cond, *((q,)+args+(loc,)))
            loc, goodargs = goodargs[-1], goodargs[:-1]
            place(output,cond,self._isf(*goodargs) + loc)  # PB same as ticket 766

        if output.ndim == 0:
            return output[()]
        return output

    def stats(self, *args, **kwds):
        """
        Some statistics of the given discrete RV.

        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).
        moments : string, optional
            Composed of letters ['mvsk'] defining which moments to compute:

              - 'm' = mean,
              - 'v' = variance,
              - 's' = (Fisher's) skew,
              - 'k' = (Fisher's) kurtosis.

            The default is'mv'.

        Returns
        -------
        stats : sequence
            of requested moments.

        """
        try:
            kwds["moments"] = kwds.pop("moment") # test suite is full of these; a feature?
        except KeyError:
            pass
        args, loc, _, moments = self._parse_args_stats(*args, **kwds)
        loc = asarray(loc)
        args = tuple(map(asarray,args))
        cond = self._argcheck(*args) & (loc == loc)

        signature = inspect.getargspec(get_method_function(self._stats))
        if (signature[2] is not None) or ('moments' in signature[0]):
            mu, mu2, g1, g2 = self._stats(*args,**{'moments':moments})
        else:
            mu, mu2, g1, g2 = self._stats(*args)
        if g1 is None:
            mu3 = None
        else:
            mu3 = g1 * np.power(mu2, 1.5)
        default = valarray(shape(cond), self.badvalue)
        output = []

        # Use only entries that are valid in calculation
        goodargs = argsreduce(cond, *(args+(loc,)))
        loc, goodargs = goodargs[-1], goodargs[:-1]

        if 'm' in moments:
            if mu is None:
                mu = self._munp(1.0,*goodargs)
            out0 = default.copy()
            place(out0,cond,mu+loc)
            output.append(out0)

        if 'v' in moments:
            if mu2 is None:
                mu2p = self._munp(2.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)
                mu2 = mu2p - mu*mu
            out0 = default.copy()
            place(out0,cond,mu2)
            output.append(out0)

        if 's' in moments:
            if g1 is None:
                mu3p = self._munp(3.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    mu2 = mu2p - mu*mu
                mu3 = mu3p - 3*mu*mu2 - mu**3
                g1 = mu3 / np.power(mu2, 1.5)
            out0 = default.copy()
            place(out0,cond,g1)
            output.append(out0)

        if 'k' in moments:
            if g2 is None:
                mu4p = self._munp(4.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    mu2 = mu2p - mu*mu
                if mu3 is None:
                    mu3p = self._munp(3.0,*goodargs)
                    mu3 = mu3p - 3*mu*mu2 - mu**3
                mu4 = mu4p - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4
                g2 = mu4 / mu2**2.0 - 3.0
            out0 = default.copy()
            place(out0,cond,g2)
            output.append(out0)

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)

    def moment(self, n, *args, **kwds):
        """
        n'th non-central moment of the distribution

        Parameters
        ----------
        n : int, n>=1
            order of moment
        arg1, arg2, arg3,... : float
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : float, optional
            location parameter (default=0)
        scale : float, optional
            scale parameter (default=1)

        """
        args, loc, scale = self._parse_args(*args, **kwds)
        if not (self._argcheck(*args) and (scale > 0)):
            return nan
        if (floor(n) != n):
            raise ValueError("Moment must be an integer.")
        if (n < 0):
            raise ValueError("Moment must be positive.")
        mu, mu2, g1, g2 = None, None, None, None
        if (n > 0) and (n < 5):
            signature = inspect.getargspec(get_method_function(self._stats))
            if (signature[2] is not None) or ('moments' in signature[0]):
                dict = {'moments':{1:'m',2:'v',3:'vs',4:'vk'}[n]}
            else:
                dict = {}
            mu, mu2, g1, g2 = self._stats(*args,**dict)
        val = _moment_from_stats(n, mu, mu2, g1, g2, self._munp, args)

        # Convert to transformed  X = L + S*Y
        # so E[X^n] = E[(L+S*Y)^n] = L^n sum(comb(n,k)*(S/L)^k E[Y^k],k=0...n)
        if loc == 0:
            return scale**n * val
        else:
            result = 0
            fac = float(scale) / float(loc)
            for k in range(n):
                valk = _moment_from_stats(k, mu, mu2, g1, g2, self._munp, args)
                result += comb(n,k,exact=True)*(fac**k) * valk
            result += fac**n * val
            return result * loc**n

    def freeze(self, *args, **kwds):
        return rv_frozen(self, *args, **kwds)

    def _entropy(self, *args):
        if hasattr(self,'pk'):
            return entropy(self.pk)
        else:
            mu = int(self.stats(*args, **{'moments':'m'}))
            val = self.pmf(mu,*args)
            ent = -special.xlogy(val, val)
            k = 1
            term = 1.0
            while (abs(term) > eps):
                val = self.pmf(mu+k,*args)
                term = -special.xlogy(val, val)
                val = self.pmf(mu-k,*args)
                term -= special.xlogy(val, val)
                k += 1
                ent += term
            return ent

    def entropy(self, *args, **kwds):
        args, loc, _ = self._parse_args(*args, **kwds)
        loc = asarray(loc)
        args = list(map(asarray,args))
        cond0 = self._argcheck(*args) & (loc == loc)
        output = zeros(shape(cond0),'d')
        place(output,(1-cond0),self.badvalue)
        goodargs = argsreduce(cond0, *args)
        # np.vectorize doesn't work when numargs == 0 in numpy 1.5.1
        if self.numargs == 0:
            place(output, cond0, self._entropy())
        else:
            place(output, cond0, self.vecentropy(*goodargs))

        return output

    def __call__(self, *args, **kwds):
        return self.freeze(*args,**kwds)

    def expect(self, func=None, args=(), loc=0, lb=None, ub=None, conditional=False):
        """
        Calculate expected value of a function with respect to the distribution
        for discrete distribution

        Parameters
        ----------
        fn : function (default: identity mapping)
            Function for which sum is calculated. Takes only one argument.
        args : tuple
            argument (parameters) of the distribution
        lb, ub : numbers, optional
            lower and upper bound for integration, default is set to the support
            of the distribution, lb and ub are inclusive (ul<=k<=ub)
        conditional : bool, optional
            Default is False.
            If true then the expectation is corrected by the conditional
            probability of the integration interval. The return value is the
            expectation of the function, conditional on being in the given
            interval (k such that ul<=k<=ub).

        Returns
        -------
        expect : float
            Expected value.

        Notes
        -----
        * function is not vectorized
        * accuracy: uses self.moment_tol as stopping criterium
          for heavy tailed distribution e.g. zipf(4), accuracy for
          mean, variance in example is only 1e-5,
          increasing precision (moment_tol) makes zipf very slow
        * suppnmin=100 internal parameter for minimum number of points to evaluate
          could be added as keyword parameter, to evaluate functions with
          non-monotonic shapes, points include integers in (-suppnmin, suppnmin)
        * uses maxcount=1000 limits the number of points that are evaluated
          to break loop for infinite sums
          (a maximum of suppnmin+1000 positive plus suppnmin+1000 negative
          integers are evaluated)

        """

        # moment_tol = 1e-12 # increase compared to self.moment_tol,
        # too slow for only small gain in precision for zipf

        # avoid endless loop with unbound integral, eg. var of zipf(2)
        maxcount = 1000
        suppnmin = 100  # minimum number of points to evaluate (+ and -)

        if func is None:
            def fun(x):
                # loc and args from outer scope
                return (x+loc)*self._pmf(x, *args)
        else:
            def fun(x):
                # loc and args from outer scope
                return func(x+loc)*self._pmf(x, *args)
        # used pmf because _pmf does not check support in randint
        # and there might be problems(?) with correct self.a, self.b at this stage
        # maybe not anymore, seems to work now with _pmf

        self._argcheck(*args)  # (re)generate scalar self.a and self.b
        if lb is None:
            lb = (self.a)
        else:
            lb = lb - loc   # convert bound for standardized distribution
        if ub is None:
            ub = (self.b)
        else:
            ub = ub - loc   # convert bound for standardized distribution
        if conditional:
            if np.isposinf(ub)[()]:
                # work around bug: stats.poisson.sf(stats.poisson.b, 2) is nan
                invfac = 1 - self.cdf(lb-1,*args)
            else:
                invfac = 1 - self.cdf(lb-1,*args) - self.sf(ub,*args)
        else:
            invfac = 1.0

        tot = 0.0
        low, upp = self._ppf(0.001, *args), self._ppf(0.999, *args)
        low = max(min(-suppnmin, low), lb)
        upp = min(max(suppnmin, upp), ub)
        supp = np.arange(low, upp+1, self.inc)  # check limits
        # print 'low, upp', low, upp
        tot = np.sum(fun(supp))
        diff = 1e100
        pos = upp + self.inc
        count = 0

        # handle cases with infinite support

        while (pos <= ub) and (diff > self.moment_tol) and count <= maxcount:
            diff = fun(pos)
            tot += diff
            pos += self.inc
            count += 1

        if self.a < 0:  # handle case when self.a = -inf
            diff = 1e100
            pos = low - self.inc
            while (pos >= lb) and (diff > self.moment_tol) and count <= maxcount:
                diff = fun(pos)
                tot += diff
                pos -= self.inc
                count += 1
        if count > maxcount:
            warnings.warn('expect(): sum did not converge', RuntimeWarning)
        return tot/invfac


class binom_gen(rv_discrete):
    """A binomial discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `binom` is::

       binom.pmf(k) = choose(n,k) * p**k * (1-p)**(n-k)

    for ``k`` in ``{0,1,...,n}``.

    `binom` takes ``n`` and ``p`` as shape parameters.

    %(example)s

    """
    def _rvs(self, n, p):
        return mtrand.binomial(n,p,self._size)

    def _argcheck(self, n, p):
        self.b = n
        return (n >= 0) & (p >= 0) & (p <= 1)

    def _logpmf(self, x, n, p):
        k = floor(x)
        combiln = (gamln(n+1) - (gamln(k+1) +
                                           gamln(n-k+1)))
        return combiln + special.xlogy(k,p) + special.xlog1py(n-k, -p)

    def _pmf(self, x, n, p):
        return exp(self._logpmf(x, n, p))

    def _cdf(self, x, n, p):
        k = floor(x)
        vals = special.bdtr(k,n,p)
        return vals

    def _sf(self, x, n, p):
        k = floor(x)
        return special.bdtrc(k,n,p)

    def _ppf(self, q, n, p):
        vals = ceil(special.bdtrik(q,n,p))
        vals1 = vals-1
        temp = special.bdtr(vals1,n,p)
        return where(temp >= q, vals1, vals)

    def _stats(self, n, p):
        q = 1.0-p
        mu = n * p
        var = n * p * q
        g1 = (q-p) / sqrt(n*p*q)
        g2 = (1.0-6*p*q)/(n*p*q)
        return mu, var, g1, g2

    def _entropy(self, n, p):
        k = r_[0:n + 1]
        vals = self._pmf(k, n, p)
        h = -sum(special.xlogy(vals, vals), axis=0)
        return h
binom = binom_gen(name='binom')


class bernoulli_gen(binom_gen):
    """A Bernoulli discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `bernoulli` is::

       bernoulli.pmf(k) = 1-p  if k = 0
                        = p    if k = 1

    for ``k`` in ``{0,1}``.

    `bernoulli` takes ``p`` as shape parameter.

    %(example)s

    """
    def _rvs(self, p):
        return binom_gen._rvs(self, 1, p)

    def _argcheck(self, p):
        return (p >= 0) & (p <= 1)

    def _logpmf(self, x, p):
        return binom._logpmf(x, 1, p)

    def _pmf(self, x, p):
        return binom._pmf(x, 1, p)

    def _cdf(self, x, p):
        return binom._cdf(x, 1, p)

    def _sf(self, x, p):
        return binom._sf(x, 1, p)

    def _ppf(self, q, p):
        return binom._ppf(q, 1, p)

    def _stats(self, p):
        return binom._stats(1, p)

    def _entropy(self, p):
        h = -special.xlogy(p, p) - special.xlogy(1 - p, 1 - p)
        return h
bernoulli = bernoulli_gen(b=1,name='bernoulli')


class nbinom_gen(rv_discrete):
    """A negative binomial discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `nbinom` is::

         nbinom.pmf(k) = choose(k+n-1, n-1) * p**n * (1-p)**k

    for ``k >= 0``.

    `nbinom` takes ``n`` and ``p`` as shape parameters.

    %(example)s

    """
    def _rvs(self, n, p):
        return mtrand.negative_binomial(n, p, self._size)

    def _argcheck(self, n, p):
        return (n >= 0) & (p >= 0) & (p <= 1)

    def _pmf(self, x, n, p):
        return exp(self._logpmf(x, n, p))

    def _logpmf(self, x, n, p):
        coeff = gamln(n+x) - gamln(x+1) - gamln(n)
        return coeff + n*log(p) + x*log(1-p)

    def _cdf(self, x, n, p):
        k = floor(x)
        return special.betainc(n, k+1, p)

    def _sf_skip(self, x, n, p):
        # skip because special.nbdtrc doesn't work for 0<n<1
        k = floor(x)
        return special.nbdtrc(k,n,p)

    def _ppf(self, q, n, p):
        vals = ceil(special.nbdtrik(q,n,p))
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1,n,p)
        return where(temp >= q, vals1, vals)

    def _stats(self, n, p):
        Q = 1.0 / p
        P = Q - 1.0
        mu = n*P
        var = n*P*Q
        g1 = (Q+P)/sqrt(n*P*Q)
        g2 = (1.0 + 6*P*Q) / (n*P*Q)
        return mu, var, g1, g2
nbinom = nbinom_gen(name='nbinom')


class geom_gen(rv_discrete):
    """A geometric discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `geom` is::

        geom.pmf(k) = (1-p)**(k-1)*p

    for ``k >= 1``.

    `geom` takes ``p`` as shape parameter.

    %(example)s

    """
    def _rvs(self, p):
        return mtrand.geometric(p, size=self._size)

    def _argcheck(self, p):
        return (p <= 1) & (p >= 0)

    def _pmf(self, k, p):
        return (1-p)**(k-1) * p

    def _logpmf(self, k, p):
        return (k-1)*log(1-p) + log(p)

    def _cdf(self, x, p):
        k = floor(x)
        return (1.0-(1.0-p)**k)

    def _sf(self, x, p):
        k = floor(x)
        return (1.0-p)**k

    def _ppf(self, q, p):
        vals = ceil(log(1.0-q)/log(1-p))
        temp = 1.0-(1.0-p)**(vals-1)
        return where((temp >= q) & (vals > 0), vals-1, vals)

    def _stats(self, p):
        mu = 1.0/p
        qr = 1.0-p
        var = qr / p / p
        g1 = (2.0-p) / sqrt(qr)
        g2 = numpy.polyval([1,-6,6],p)/(1.0-p)
        return mu, var, g1, g2
geom = geom_gen(a=1,name='geom', longname="A geometric")


class hypergeom_gen(rv_discrete):
    """A hypergeometric discrete random variable.

    The hypergeometric distribution models drawing objects from a bin.
    M is the total number of objects, n is total number of Type I objects.
    The random variate represents the number of Type I objects in N drawn
    without replacement from the total population.

    %(before_notes)s

    Notes
    -----
    The probability mass function is defined as::

        pmf(k, M, n, N) = choose(n, k) * choose(M - n, N - k) / choose(M, N),
                                               for N - (M-n) <= k <= min(m,N)

    Examples
    --------
    >>> from scipy.stats import hypergeom

    Suppose we have a collection of 20 animals, of which 7 are dogs.  Then if
    we want to know the probability of finding a given number of dogs if we
    choose at random 12 of the 20 animals, we can initialize a frozen
    distribution and plot the probability mass function:

    >>> [M, n, N] = [20, 7, 12]
    >>> rv = hypergeom(M, n, N)
    >>> x = np.arange(0, n+1)
    >>> pmf_dogs = rv.pmf(x)

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, pmf_dogs, 'bo')
    >>> ax.vlines(x, 0, pmf_dogs, lw=2)
    >>> ax.set_xlabel('# of dogs in our group of chosen animals')
    >>> ax.set_ylabel('hypergeom PMF')
    >>> plt.show()

    Instead of using a frozen distribution we can also use `hypergeom`
    methods directly.  To for example obtain the cumulative distribution
    function, use:

    >>> prb = hypergeom.cdf(x, M, n, N)

    And to generate random numbers:

    >>> R = hypergeom.rvs(M, n, N, size=10)

    """
    def _rvs(self, M, n, N):
        return mtrand.hypergeometric(n,M-n,N,size=self._size)

    def _argcheck(self, M, n, N):
        cond = rv_discrete._argcheck(self,M,n,N)
        cond &= (n <= M) & (N <= M)
        self.a = max(N-(M-n), 0)
        self.b = min(n,N)
        return cond

    def _logpmf(self, k, M, n, N):
        tot, good = M, n
        bad = tot - good
        return gamln(good+1) - gamln(good-k+1) - gamln(k+1) + gamln(bad+1) \
            - gamln(bad-N+k+1) - gamln(N-k+1) - gamln(tot+1) + gamln(tot-N+1) \
            + gamln(N+1)

    def _pmf(self, k, M, n, N):
        # same as the following but numerically more precise
        # return comb(good,k) * comb(bad,N-k) / comb(tot,N)
        return exp(self._logpmf(k, M, n, N))

    def _stats(self, M, n, N):
        tot, good = M, n
        n = good*1.0
        m = (tot-good)*1.0
        N = N*1.0
        tot = m+n
        p = n/tot
        mu = N*p
        var = m*n*N*(tot-N)*1.0/(tot*tot*(tot-1))
        g1 = (m - n)*(tot-2*N) / (tot-2.0)*sqrt((tot-1.0)/(m*n*N*(tot-N)))
        m2, m3, m4, m5 = m**2, m**3, m**4, m**5
        n2, n3, n4, n5 = n**2, n**2, n**4, n**5
        g2 = m3 - m5 + n*(3*m2-6*m3+m4) + 3*m*n2 - 12*m2*n2 + 8*m3*n2 + n3 \
             - 6*m*n3 + 8*m2*n3 + m*n4 - n5 - 6*m3*N + 6*m4*N + 18*m2*n*N \
             - 6*m3*n*N + 18*m*n2*N - 24*m2*n2*N - 6*n3*N - 6*m*n3*N \
             + 6*n4*N + N*N*(6*m2 - 6*m3 - 24*m*n + 12*m2*n + 6*n2 +
                             12*m*n2 - 6*n3)
        return mu, var, g1, g2

    def _entropy(self, M, n, N):
        k = r_[N - (M - n):min(n, N) + 1]
        vals = self.pmf(k, M, n, N)
        h = -sum(special.xlogy(vals, vals), axis=0)
        return h

    def _sf(self, k, M, n, N):
        """More precise calculation, 1 - cdf doesn't cut it."""
        # This for loop is needed because `k` can be an array. If that's the
        # case, the sf() method makes M, n and N arrays of the same shape. We
        # therefore unpack all inputs args, so we can do the manual integration.
        res = []
        for quant, tot, good, draw in zip(k, M, n, N):
            # Manual integration over probability mass function. More accurate
            # than integrate.quad.
            k2 = np.arange(quant + 1, draw + 1)
            res.append(np.sum(self._pmf(k2, tot, good, draw)))
        return np.asarray(res)
hypergeom = hypergeom_gen(name='hypergeom')


# FIXME: Fails _cdfvec
class logser_gen(rv_discrete):
    """A Logarithmic (Log-Series, Series) discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `logser` is::

        logser.pmf(k) = - p**k / (k*log(1-p))

    for ``k >= 1``.

    `logser` takes ``p`` as shape parameter.

    %(example)s

    """
    def _rvs(self, p):
        # looks wrong for p>0.5, too few k=1
        # trying to use generic is worse, no k=1 at all
        return mtrand.logseries(p, size=self._size)

    def _argcheck(self, p):
        return (p > 0) & (p < 1)

    def _pmf(self, k, p):
        return -p**k * 1.0 / k / log(1 - p)

    def _stats(self, p):
        r = log(1 - p)
        mu = p / (p - 1.0) / r
        mu2p = -p / r / (p - 1.0)**2
        var = mu2p - mu*mu
        mu3p = -p / r * (1.0+p) / (1.0 - p)**3
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / np.power(var, 1.5)

        mu4p = -p / r * (1.0 / (p-1)**2 - 6*p / (p - 1)**3 +
                          6*p*p / (p-1)**4)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / var**2 - 3.0
        return mu, var, g1, g2
logser = logser_gen(a=1,name='logser', longname='A logarithmic')


class poisson_gen(rv_discrete):
    """A Poisson discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `poisson` is::

        poisson.pmf(k) = exp(-mu) * mu**k / k!

    for ``k >= 0``.

    `poisson` takes ``mu`` as shape parameter.

    %(example)s

    """
    def _rvs(self, mu):
        return mtrand.poisson(mu, self._size)

    def _logpmf(self, k, mu):
        Pk = k*log(mu)-gamln(k+1) - mu
        return Pk

    def _pmf(self, k, mu):
        return exp(self._logpmf(k, mu))

    def _cdf(self, x, mu):
        k = floor(x)
        return special.pdtr(k,mu)

    def _sf(self, x, mu):
        k = floor(x)
        return special.pdtrc(k,mu)

    def _ppf(self, q, mu):
        vals = ceil(special.pdtrik(q,mu))
        vals1 = vals-1
        temp = special.pdtr(vals1,mu)
        return where((temp >= q), vals1, vals)

    def _stats(self, mu):
        var = mu
        tmp = asarray(mu)
        g1 = sqrt(1.0 / tmp)
        g2 = 1.0 / tmp
        return mu, var, g1, g2
poisson = poisson_gen(name="poisson", longname='A Poisson')


class planck_gen(rv_discrete):
    """A Planck discrete exponential random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `planck` is::

        planck.pmf(k) = (1-exp(-lambda_))*exp(-lambda_*k)

    for ``k*lambda_ >= 0``.

    `planck` takes ``lambda_`` as shape parameter.

    %(example)s

    """
    def _argcheck(self, lambda_):
        if (lambda_ > 0):
            self.a = 0
            self.b = inf
            return 1
        elif (lambda_ < 0):
            self.a = -inf
            self.b = 0
            return 1
        else:
            return 0

    def _pmf(self, k, lambda_):
        fact = (1-exp(-lambda_))
        return fact*exp(-lambda_*k)

    def _cdf(self, x, lambda_):
        k = floor(x)
        return 1-exp(-lambda_*(k+1))

    def _ppf(self, q, lambda_):
        vals = ceil(-1.0/lambda_ * log1p(-q)-1)
        vals1 = (vals-1).clip(self.a, np.inf)
        temp = self._cdf(vals1, lambda_)
        return where(temp >= q, vals1, vals)

    def _stats(self, lambda_):
        mu = 1/(exp(lambda_)-1)
        var = exp(-lambda_)/(expm1(-lambda_))**2
        g1 = 2*cosh(lambda_/2.0)
        g2 = 4+2*cosh(lambda_)
        return mu, var, g1, g2

    def _entropy(self, lambda_):
        l = lambda_
        C = (1-exp(-l))
        return l*exp(-l)/C - log(C)
planck = planck_gen(name='planck',longname='A discrete exponential ')


class boltzmann_gen(rv_discrete):
    """A Boltzmann (Truncated Discrete Exponential) random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `boltzmann` is::

        boltzmann.pmf(k) = (1-exp(-lambda_)*exp(-lambda_*k)/(1-exp(-lambda_*N))

    for ``k = 0,...,N-1``.

    `boltzmann` takes ``lambda_`` and ``N`` as shape parameters.

    %(example)s

    """
    def _pmf(self, k, lambda_, N):
        fact = (1-exp(-lambda_))/(1-exp(-lambda_*N))
        return fact*exp(-lambda_*k)

    def _cdf(self, x, lambda_, N):
        k = floor(x)
        return (1-exp(-lambda_*(k+1)))/(1-exp(-lambda_*N))

    def _ppf(self, q, lambda_, N):
        qnew = q*(1-exp(-lambda_*N))
        vals = ceil(-1.0/lambda_ * log(1-qnew)-1)
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1, lambda_, N)
        return where(temp >= q, vals1, vals)

    def _stats(self, lambda_, N):
        z = exp(-lambda_)
        zN = exp(-lambda_*N)
        mu = z/(1.0-z)-N*zN/(1-zN)
        var = z/(1.0-z)**2 - N*N*zN/(1-zN)**2
        trm = (1-zN)/(1-z)
        trm2 = (z*trm**2 - N*N*zN)
        g1 = z*(1+z)*trm**3 - N**3*zN*(1+zN)
        g1 = g1 / trm2**(1.5)
        g2 = z*(1+4*z+z*z)*trm**4 - N**4 * zN*(1+4*zN+zN*zN)
        g2 = g2 / trm2 / trm2
        return mu, var, g1, g2
boltzmann = boltzmann_gen(name='boltzmann',
        longname='A truncated discrete exponential ')


class randint_gen(rv_discrete):
    """A uniform discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `randint` is::

        randint.pmf(k) = 1./(high - low)

    for ``k = low, ..., high - 1``.

    `randint` takes ``low`` and ``high`` as shape parameters.

    Note the difference to the numpy ``random_integers`` which
    returns integers on a *closed* interval ``[low, high]``.

    %(example)s

    """
    def _argcheck(self, low, high):
        self.a = low
        self.b = high - 1
        return (high > low)

    def _pmf(self, k, low, high):
        fact = 1.0 / (high - low)
        return fact

    def _cdf(self, x, low, high):
        k = floor(x)
        return (k-low + 1) * 1.0 / (high - low)

    def _ppf(self, q, low, high):
        vals = ceil(q*(high - low) + low) - 1
        vals1 = (vals-1).clip(low, high)
        temp = self._cdf(vals1, low, high)
        return where(temp >= q, vals1, vals)

    def _stats(self, low, high):
        m2, m1 = asarray(high), asarray(low)
        mu = (m2 + m1 - 1.0) / 2
        d = m2 - m1
        var = (d-1)*(d+1.0)/12.0
        g1 = 0.0
        g2 = -6.0/5.0*(d*d+1.0)/(d-1.0)*(d+1.0)
        return mu, var, g1, g2

    def _rvs(self, low, high=None):
        """An array of *size* random integers >= ``low`` and < ``high``.

        If ``high`` is ``None``, then range is >=0  and < low
        """
        return mtrand.randint(low, high, self._size)

    def _entropy(self, low, high):
        return log(high - low)
randint = randint_gen(name='randint',longname='A discrete uniform '
                      '(random integer)')


# FIXME: problems sampling.
class zipf_gen(rv_discrete):
    """A Zipf discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `zipf` is::

        zipf.pmf(k) = 1/(zeta(a)*k**a)

    for ``k >= 1``.

    `zipf` takes ``a`` as shape parameter.

    %(example)s

    """
    def _rvs(self, a):
        return mtrand.zipf(a, size=self._size)

    def _argcheck(self, a):
        return a > 1

    def _pmf(self, k, a):
        Pk = 1.0 / asarray(special.zeta(a,1) * k**a)
        return Pk

    def _munp(self, n, a):
        return special.zeta(a-n,1) / special.zeta(a,1)

    def _stats(self, a):
        sv = special.errprint(0)
        fac = asarray(special.zeta(a,1))
        mu = special.zeta(a-1.0,1)/fac
        mu2p = special.zeta(a-2.0,1)/fac
        var = mu2p - mu*mu
        mu3p = special.zeta(a-3.0,1)/fac
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / asarray(np.power(var, 1.5))

        mu4p = special.zeta(a-4.0,1)/fac
        sv = special.errprint(sv)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / asarray(var**2) - 3.0
        return mu, var, g1, g2
zipf = zipf_gen(a=1,name='zipf', longname='A Zipf')


class dlaplace_gen(rv_discrete):
    """A  Laplacian discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `dlaplace` is::

        dlaplace.pmf(k) = tanh(a/2) * exp(-a*abs(k))

    for ``a >0``.

    `dlaplace` takes ``a`` as shape parameter.

    %(example)s

    """
    def _pmf(self, k, a):
        return tanh(a/2.0)*exp(-a*abs(k))

    def _cdf(self, x, a):
        k = floor(x)
        ind = (k >= 0)
        const = exp(a)+1
        return where(ind, 1.0-exp(-a*k)/const, exp(a*(k+1))/const)

    def _ppf(self, q, a):
        const = 1.0/(1+exp(-a))
        cons2 = 1+exp(a)
        ind = q < const
        vals = ceil(where(ind, log(q*cons2)/a-1, -log((1-q)*cons2)/a))
        vals1 = (vals-1)
        temp = self._cdf(vals1, a)
        return where(temp >= q, vals1, vals)

    def _stats(self, a):
        ea = exp(a)
        mu2 = 2.*ea/(ea-1.)**2
        mu4 = 2.*ea*(ea**2+10.*ea+1.) / (ea-1.)**4
        return 0., mu2, 0., mu4/mu2**2 - 3.

    def _entropy(self, a):
        return a / sinh(a) - log(tanh(a/2.0))
dlaplace = dlaplace_gen(a=-inf,
                        name='dlaplace', longname='A discrete Laplacian')


class skellam_gen(rv_discrete):
    """A  Skellam discrete random variable.

    %(before_notes)s

    Notes
    -----
    Probability distribution of the difference of two correlated or
    uncorrelated Poisson random variables.

    Let k1 and k2 be two Poisson-distributed r.v. with expected values
    lam1 and lam2. Then, ``k1 - k2`` follows a Skellam distribution with
    parameters ``mu1 = lam1 - rho*sqrt(lam1*lam2)`` and
    ``mu2 = lam2 - rho*sqrt(lam1*lam2)``, where rho is the correlation
    coefficient between k1 and k2. If the two Poisson-distributed r.v.
    are independent then ``rho = 0``.

    Parameters mu1 and mu2 must be strictly positive.

    For details see: http://en.wikipedia.org/wiki/Skellam_distribution

    `skellam` takes ``mu1`` and ``mu2`` as shape parameters.

    %(example)s

    """
    def _rvs(self, mu1, mu2):
        n = self._size
        return np.random.poisson(mu1, n)-np.random.poisson(mu2, n)

    def _pmf(self, x, mu1, mu2):
        px = np.where(x < 0, ncx2.pdf(2*mu2, 2*(1-x), 2*mu1)*2,
                         ncx2.pdf(2*mu1, 2*(x+1), 2*mu2)*2)
        # ncx2.pdf() returns nan's for extremely low probabilities
        return px

    def _cdf(self, x, mu1, mu2):
        x = np.floor(x)
        px = np.where(x < 0, ncx2.cdf(2*mu2, -2*x, 2*mu1),
                         1-ncx2.cdf(2*mu1, 2*(x+1), 2*mu2))
        return px

    def _stats(self, mu1, mu2):
        mean = mu1 - mu2
        var = mu1 + mu2
        g1 = mean / np.sqrt((var)**3)
        g2 = 1 / var
        return mean, var, g1, g2
skellam = skellam_gen(a=-np.inf, name="skellam", longname='A Skellam')
