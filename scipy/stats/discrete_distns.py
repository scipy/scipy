#
# Author:  Travis Oliphant  2002-2011 with contributions from
#          SciPy Developers 2004-2011
#
from __future__ import division, print_function, absolute_import

import sys
import warnings

from scipy.misc import derivative
from scipy import special
from scipy.special import gammaln as gamln

from numpy import (all, where, arange, putmask, ravel, take, ones, sum, shape,
                   product, reshape, zeros, floor, logical_and, log, sqrt, exp,
                   tanh, ndarray, cosh,
                   sinh, newaxis, log1p, expm1)

from numpy import (atleast_1d, polyval, ceil, place, extract, any, argsort,
                   argmax, vectorize, r_, asarray, nan, inf, pi, isinf, NINF,
                   empty)

import numpy as np
import numpy.random as mtrand

from ._distn_infrastructure import (
        instancemethod,
        rv_generic, docdict_discrete, argsreduce, valarray,
        _lazywhere,
        _ncx2_pdf, _ncx2_cdf,
        )

__all__ = [
    'entropy', 'rv_discrete', 'binom', 'bernoulli',
    'nbinom', 'geom', 'hypergeom', 'logser', 'poisson', 'planck',
    'boltzmann', 'randint', 'zipf', 'dlaplace', 'skellam'
]

eps = np.finfo(float).eps

import types
from scipy.misc import doccer


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
        mask = qk == 0.0
        qk[mask] = 1.0  # Avoid the divide-by-zero warning
        quotient = pk / qk
        vec = -special.xlogy(pk, quotient)
        vec[mask & (pk != 0.0)] = -inf
        vec[mask & (pk == 0.0)] = 0.0
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
    indx = argmax((self.xk > xk), axis=-1)-1
    return self.F[self.xk[indx]]


def _drv_ppf(self, q, *args):
    indx = argmax((self.qvals >= q), axis=-1)
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
    ulimit = max(1000, (min(self.b, 1000) + max(self.a, -1000))/2.0)
    llimit = min(-1000, (min(self.b, 1000) + max(self.a, -1000))/2.0)

    while (pos <= self.b) and ((pos <= ulimit) or
                               (diff > self.moment_tol)):
        diff = np.power(pos, n) * self.pmf(pos, *args)
        # use pmf because _pmf does not check support in randint and there
        # might be problems ? with correct self.a, self.b at this stage
        tot += diff
        pos += self.inc
        count += 1

    if self.a < 0:  # handle case when self.a = -inf
        diff = 1e100
        pos = -self.inc
        while (pos >= self.a) and ((pos >= llimit) or
                                   (diff > self.moment_tol)):
            diff = np.power(pos, n) * self.pmf(pos, *args)
            # using pmf instead of _pmf, see above
            tot += diff
            pos -= self.inc
            count += 1
    return tot


def _drv2_ppfsingle(self, q, *args):  # Use basic bisection algorithm
    b = self.b
    a = self.a
    if isinf(b):            # Be sure ending point is > q
        b = int(max(100*q, 10))
        while 1:
            if b >= self.b:
                qb = 1.0
                break
            qb = self._cdf(b, *args)
            if (qb < q):
                b += 10
            else:
                break
    else:
        qb = 1.0
    if isinf(a):    # be sure starting point < q
        a = int(min(-100*q, -10))
        while 1:
            if a <= self.a:
                qb = 0.0
                break
            qa = self._cdf(a, *args)
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
    # python -c "from scipy.stats import zipf;print zipf.ppf(0.01, 2)" wrong
    # python -c "from scipy.stats import zipf;print zipf.ppf([0.01, 0.61, 0.77, 0.83], 2)"
    # python -c "from scipy.stats import logser;print logser.ppf([0.1, 0.66, 0.86, 0.93], 0.6)"
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
        non-central n-th moment of the distribution.  May not work for array
        arguments.

    generic.stats(<shape(s)>, loc=0, moments='mv')
        mean('m', axis=0), variance('v'), skew('s'), and/or kurtosis('k')

    generic.entropy(<shape(s)>, loc=0)
        entropy of the RV

    generic.expect(func=None, args=(), loc=0, lb=None, ub=None,
            conditional=False)
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
                 moment_tol=1e-8, values=None, inc=1, longname=None,
                 shapes=None, extradoc=None):

        super(rv_discrete, self).__init__()

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
        self._cdfvec = vectorize(self._cdf_single, otypes='d')
        self.return_integers = 1
        self.vecentropy = vectorize(self._entropy)
        self.shapes = shapes
        self.extradoc = extradoc

        if values is not None:
            self.xk, self.pk = values
            self.return_integers = 0
            indx = argsort(ravel(self.xk))
            self.xk = take(ravel(self.xk), indx, 0)
            self.pk = take(ravel(self.pk), indx, 0)
            self.a = self.xk[0]
            self.b = self.xk[-1]
            self.P = make_dict(self.xk, self.pk)
            self.qvals = np.cumsum(self.pk, axis=0)
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
                                      # scale=1 for discrete RVs
                                      locscale_out='loc, 1',
                                      morevars=globals())
        else:
            self._construct_argparser(names_to_inspect=['_pmf', '_cdf'],
                                      locscale_in='loc=0',
                                      # scale=1 for discrete RVs
                                      locscale_out='loc, 1',
                                      morevars=globals())

            # nin correction needs to be after we know numargs
            # correct nin for generic moment vectorization
            _vec_generic_moment = vectorize(_drv2_moment, otypes='d')
            _vec_generic_moment.nin = self.numargs + 2
            self.generic_moment = instancemethod(_vec_generic_moment,
                                                 self, rv_discrete)

            # backwards compatibility
            self.vec_generic_moment = _vec_generic_moment

            # correct nin for ppf vectorization
            _vppf = vectorize(_drv2_ppfsingle, otypes='d')
            _vppf.nin = self.numargs + 2  # +1 is for self
            self._ppfvec = instancemethod(_vppf,
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
                self._construct_default_doc(longname=longname,
                                            extradoc=extradoc)
            else:
                self._construct_doc()

            #discrete RV do not have the scale parameter, remove it
            self.__doc__ = self.__doc__.replace(
                '\n    scale : array_like, '
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

    def _nonzero(self, k, *args):
        return floor(k) == k

    def _pmf(self, k, *args):
        return self._cdf(k, *args) - self._cdf(k-1, *args)

    def _logpmf(self, k, *args):
        return log(self._pmf(k, *args))

    def _cdf_single(self, k, *args):
        m = arange(int(self.a), k+1)
        return sum(self._pmf(m, *args), axis=0)

    def _cdf(self, x, *args):
        k = floor(x)
        return self._cdfvec(k, *args)

    # generic _logcdf, _sf, _logsf, _ppf, _isf, _rvs defined in rv_generic

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

    def pmf(self, k, *args, **kwds):
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
        k, loc = map(asarray, (k, loc))
        args = tuple(map(asarray, args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b) & self._nonzero(k, *args)
        cond = cond0 & cond1
        output = zeros(shape(cond), 'd')
        place(output, (1-cond0) + np.isnan(k), self.badvalue)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output, cond, self._pmf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logpmf(self, k, *args, **kwds):
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
        k, loc = map(asarray, (k, loc))
        args = tuple(map(asarray, args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b) & self._nonzero(k, *args)
        cond = cond0 & cond1
        output = empty(shape(cond), 'd')
        output.fill(NINF)
        place(output, (1-cond0) + np.isnan(k), self.badvalue)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output, cond, self._logpmf(*goodargs))
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
        k, loc = map(asarray, (k, loc))
        args = tuple(map(asarray, args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k < self.b)
        cond2 = (k >= self.b)
        cond = cond0 & cond1
        output = zeros(shape(cond), 'd')
        place(output, (1-cond0) + np.isnan(k), self.badvalue)
        place(output, cond2*(cond0 == cond0), 1.0)

        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output, cond, self._cdf(*goodargs))
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
        k, loc = map(asarray, (k, loc))
        args = tuple(map(asarray, args))
        k = asarray((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k < self.b)
        cond2 = (k >= self.b)
        cond = cond0 & cond1
        output = empty(shape(cond), 'd')
        output.fill(NINF)
        place(output, (1-cond0) + np.isnan(k), self.badvalue)
        place(output, cond2*(cond0 == cond0), 0.0)

        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output, cond, self._logcdf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def sf(self, k, *args, **kwds):
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
        k, loc = map(asarray, (k, loc))
        args = tuple(map(asarray, args))
        k = asarray(k-loc)
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b)
        cond2 = (k < self.a) & cond0
        cond = cond0 & cond1
        output = zeros(shape(cond), 'd')
        place(output, (1-cond0) + np.isnan(k), self.badvalue)
        place(output, cond2, 1.0)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output, cond, self._sf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def logsf(self, k, *args, **kwds):
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
        logsf : ndarray
            Log of the survival function evaluated at `k`.

        """
        args, loc, _ = self._parse_args(*args, **kwds)
        k, loc = map(asarray, (k, loc))
        args = tuple(map(asarray, args))
        k = asarray(k-loc)
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b)
        cond2 = (k < self.a) & cond0
        cond = cond0 & cond1
        output = empty(shape(cond), 'd')
        output.fill(NINF)
        place(output, (1-cond0) + np.isnan(k), self.badvalue)
        place(output, cond2, 0.0)
        if any(cond):
            goodargs = argsreduce(cond, *((k,)+args))
            place(output, cond, self._logsf(*goodargs))
        if output.ndim == 0:
            return output[()]
        return output

    def ppf(self, q, *args, **kwds):
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
        q, loc = map(asarray, (q, loc))
        args = tuple(map(asarray, args))
        cond0 = self._argcheck(*args) & (loc == loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q == 1) & cond0
        cond = cond0 & cond1
        output = valarray(shape(cond), value=self.badvalue, typecode='d')
        # output type 'd' to handle nin and inf
        place(output, (q == 0)*(cond == cond), self.a-1)
        place(output, cond2, self.b)
        if any(cond):
            goodargs = argsreduce(cond, *((q,)+args+(loc,)))
            loc, goodargs = goodargs[-1], goodargs[:-1]
            place(output, cond, self._ppf(*goodargs) + loc)

        if output.ndim == 0:
            return output[()]
        return output

    def isf(self, q, *args, **kwds):
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
        q, loc = map(asarray, (q, loc))
        args = tuple(map(asarray, args))
        cond0 = self._argcheck(*args) & (loc == loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q == 1) & cond0
        cond = cond0 & cond1

        # same problem as with ppf; copied from ppf and changed
        output = valarray(shape(cond), value=self.badvalue, typecode='d')
        # output type 'd' to handle nin and inf
        place(output, (q == 0)*(cond == cond), self.b)
        place(output, cond2, self.a-1)

        # call place only if at least 1 valid argument
        if any(cond):
            goodargs = argsreduce(cond, *((q,)+args+(loc,)))
            loc, goodargs = goodargs[-1], goodargs[:-1]
            # PB same as ticket 766
            place(output, cond, self._isf(*goodargs) + loc)

        if output.ndim == 0:
            return output[()]
        return output

    def _entropy(self, *args):
        if hasattr(self, 'pk'):
            return entropy(self.pk)
        else:
            mu = int(self.stats(*args, **{'moments': 'm'}))
            val = self.pmf(mu, *args)
            ent = -special.xlogy(val, val)
            k = 1
            term = 1.0
            while (abs(term) > eps):
                val = self.pmf(mu+k, *args)
                term = -special.xlogy(val, val)
                val = self.pmf(mu-k, *args)
                term -= special.xlogy(val, val)
                k += 1
                ent += term
            return ent

    def expect(self, func=None, args=(), loc=0, lb=None, ub=None,
               conditional=False):
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
            lower and upper bound for integration, default is set to the
            support of the distribution, lb and ub are inclusive (ul<=k<=ub)
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
        * suppnmin=100 internal parameter for minimum number of points to
          evaluate could be added as keyword parameter, to evaluate functions
          with non-monotonic shapes, points include integers in (-suppnmin,
          suppnmin)
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
        # used pmf because _pmf does not check support in randint and there
        # might be problems(?) with correct self.a, self.b at this stage maybe
        # not anymore, seems to work now with _pmf

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
                invfac = 1 - self.cdf(lb-1, *args)
            else:
                invfac = 1 - self.cdf(lb-1, *args) - self.sf(ub, *args)
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
            while ((pos >= lb) and (diff > self.moment_tol) and
                   count <= maxcount):
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

       binom.pmf(k) = choose(n, k) * p**k * (1-p)**(n-k)

    for ``k`` in ``{0, 1,..., n}``.

    `binom` takes ``n`` and ``p`` as shape parameters.

    %(example)s

    """
    def _rvs(self, n, p):
        return mtrand.binomial(n, p, self._size)

    def _argcheck(self, n, p):
        self.b = n
        return (n >= 0) & (p >= 0) & (p <= 1)

    def _logpmf(self, x, n, p):
        k = floor(x)
        combiln = (gamln(n+1) - (gamln(k+1) + gamln(n-k+1)))
        return combiln + special.xlogy(k, p) + special.xlog1py(n-k, -p)

    def _pmf(self, x, n, p):
        return exp(self._logpmf(x, n, p))

    def _cdf(self, x, n, p):
        k = floor(x)
        vals = special.bdtr(k, n, p)
        return vals

    def _sf(self, x, n, p):
        k = floor(x)
        return special.bdtrc(k, n, p)

    def _ppf(self, q, n, p):
        vals = ceil(special.bdtrik(q, n, p))
        vals1 = vals-1
        temp = special.bdtr(vals1, n, p)
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

    for ``k`` in ``{0, 1}``.

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
bernoulli = bernoulli_gen(b=1, name='bernoulli')


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
        return special.nbdtrc(k, n, p)

    def _ppf(self, q, n, p):
        vals = ceil(special.nbdtrik(q, n, p))
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1, n, p)
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
        return np.power(1-p, k-1) * p

    def _logpmf(self, k, p):
        return (k-1) * log(1-p) + log(p)

    def _cdf(self, x, p):
        k = floor(x)
        return -np.expm1(np.log1p(-p)*k)

    def _sf(self, x, p):
        return np.exp(self._logsf(x, p))

    def _logsf(self, x, p):
        k = floor(x)
        return k*np.log1p(-p)

    def _ppf(self, q, p):
        vals = ceil(log(1.0-q)/log(1-p))
        temp = self._cdf(vals-1, p)
        return where((temp >= q) & (vals > 0), vals-1, vals)

    def _stats(self, p):
        mu = 1.0/p
        qr = 1.0-p
        var = qr / p / p
        g1 = (2.0-p) / sqrt(qr)
        g2 = np.polyval([1, -6, 6], p)/(1.0-p)
        return mu, var, g1, g2
geom = geom_gen(a=1, name='geom', longname="A geometric")


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
                                       for max(0, N - (M-n)) <= k <= min(n, N)

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
        return mtrand.hypergeometric(n, M-n, N, size=self._size)

    def _argcheck(self, M, n, N):
        cond = rv_discrete._argcheck(self, M, n, N)
        cond &= (n <= M) & (N <= M)
        self.a = max(N-(M-n), 0)
        self.b = min(n, N)
        return cond

    def _logpmf(self, k, M, n, N):
        tot, good = M, n
        bad = tot - good
        return gamln(good+1) - gamln(good-k+1) - gamln(k+1) + gamln(bad+1) \
            - gamln(bad-N+k+1) - gamln(N-k+1) - gamln(tot+1) + gamln(tot-N+1) \
            + gamln(N+1)

    def _pmf(self, k, M, n, N):
        # same as the following but numerically more precise
        # return comb(good, k) * comb(bad, N-k) / comb(tot, N)
        return exp(self._logpmf(k, M, n, N))

    def _stats(self, M, n, N):
        # tot, good, sample_size = M, n, N
        # "wikipedia".replace('N', 'M').replace('n', 'N').replace('K', 'n')
        M, n, N = 1.*M, 1.*n, 1.*N
        m = M - n
        p = n/M
        mu = N*p

        var = m*n*N*(M - N)*1.0/(M*M*(M-1))
        g1 = (m - n)*(M-2*N) / (M-2.0) * sqrt((M-1.0) / (m*n*N*(M-N)))

        g2 = M*(M+1) - 6.*N*(M-N) - 6.*n*m
        g2 *= (M-1)*M*M
        g2 += 6.*n*N*(M-N)*m*(5.*M-6)
        g2 /= n * N * (M-N) * m * (M-2.) * (M-3.)
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
        # therefore unpack all inputs args, so we can do the manual
        # integration.
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
        return -np.power(p, k) * 1.0 / k / log(1 - p)

    def _stats(self, p):
        r = log(1 - p)
        mu = p / (p - 1.0) / r
        mu2p = -p / r / (p - 1.0)**2
        var = mu2p - mu*mu
        mu3p = -p / r * (1.0+p) / (1.0 - p)**3
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / np.power(var, 1.5)

        mu4p = -p / r * (
            1.0 / (p-1)**2 - 6*p / (p - 1)**3 + 6*p*p / (p-1)**4)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / var**2 - 3.0
        return mu, var, g1, g2
logser = logser_gen(a=1, name='logser', longname='A logarithmic')


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
        return special.pdtr(k, mu)

    def _sf(self, x, mu):
        k = floor(x)
        return special.pdtrc(k, mu)

    def _ppf(self, q, mu):
        vals = ceil(special.pdtrik(q, mu))
        vals1 = vals-1
        temp = special.pdtr(vals1, mu)
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
planck = planck_gen(name='planck', longname='A discrete exponential ')


class boltzmann_gen(rv_discrete):
    """A Boltzmann (Truncated Discrete Exponential) random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `boltzmann` is::

        boltzmann.pmf(k) = (1-exp(-lambda_)*exp(-lambda_*k)/(1-exp(-lambda_*N))

    for ``k = 0,..., N-1``.

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
        p = np.ones_like(k) / (high - low)
        return np.where((k >= low) & (k < high), p, 0.)

    def _cdf(self, x, low, high):
        k = floor(x)
        return (k - low + 1.) / (high - low)

    def _ppf(self, q, low, high):
        vals = ceil(q * (high - low) + low) - 1
        vals1 = (vals - 1).clip(low, high)
        temp = self._cdf(vals1, low, high)
        return where(temp >= q, vals1, vals)

    def _stats(self, low, high):
        m2, m1 = asarray(high), asarray(low)
        mu = (m2 + m1 - 1.0) / 2
        d = m2 - m1
        var = (d*d - 1) / 12.0
        g1 = 0.0
        g2 = -6.0/5.0 * (d*d + 1.0) / (d*d - 1.0)
        return mu, var, g1, g2

    def _rvs(self, low, high=None):
        """An array of *size* random integers >= ``low`` and < ``high``.

        If ``high`` is ``None``, then range is >=0  and < low
        """
        return mtrand.randint(low, high, self._size)

    def _entropy(self, low, high):
        return log(high - low)
randint = randint_gen(name='randint', longname='A discrete uniform '
                      '(random integer)')


# FIXME: problems sampling.
class zipf_gen(rv_discrete):
    """A Zipf discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `zipf` is::

        zipf.pmf(k, a) = 1/(zeta(a) * k**a)

    for ``k >= 1``.

    `zipf` takes ``a`` as shape parameter.

    %(example)s

    """
    def _rvs(self, a):
        return mtrand.zipf(a, size=self._size)

    def _argcheck(self, a):
        return a > 1

    def _pmf(self, k, a):
        Pk = 1.0 / special.zeta(a, 1) / k**a
        return Pk

    def _munp(self, n, a):
        return _lazywhere(
            a > n + 1, (a, n),
            lambda a, n: special.zeta(a - n, 1) / special.zeta(a, 1),
            np.inf)
zipf = zipf_gen(a=1, name='zipf', longname='A Zipf')


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
        return mtrand.poisson(mu1, n) - mtrand.poisson(mu2, n)

    def _pmf(self, x, mu1, mu2):
        px = np.where(x < 0,
                _ncx2_pdf(2*mu2, 2*(1-x), 2*mu1)*2,
                _ncx2_pdf(2*mu1, 2*(1+x), 2*mu2)*2)
        # ncx2.pdf() returns nan's for extremely low probabilities
        return px

    def _cdf(self, x, mu1, mu2):
        x = np.floor(x)
        px = np.where(x < 0,
                _ncx2_cdf(2*mu2, -2*x, 2*mu1),
                1-_ncx2_cdf(2*mu1, 2*(x+1), 2*mu2))
        return px

    def _stats(self, mu1, mu2):
        mean = mu1 - mu2
        var = mu1 + mu2
        g1 = mean / np.sqrt((var)**3)
        g2 = 1 / var
        return mean, var, g1, g2
skellam = skellam_gen(a=-np.inf, name="skellam", longname='A Skellam')
