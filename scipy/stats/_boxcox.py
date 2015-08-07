"""
"""
from __future__ import division, print_function, absolute_import

import numpy as np

__all__ = ['boxcox_llf', 'boxcox', 'boxcox_normmax']


def _asarray_1d_positive(x):
    x = np.asarray(x)
    if len(x.shape) != 1:
        raise ValueError('Array must be 1d.')
    if not x.size:
        return x
    if np.any(x <= 0):
        raise ValueError("Array elements must be positive.")
    return x


def boxcox_llf(lmb, data):
    r"""The boxcox log-likelihood function.

    Parameters
    ----------
    lmb : scalar
        Parameter for Box-Cox transformation.  See `boxcox` for details.
    data : array_like
        Data to calculate Box-Cox log-likelihood for.  If `data` is
        multi-dimensional, the log-likelihood is calculated along the first
        axis.

    Returns
    -------
    llf : float or ndarray
        Box-Cox log-likelihood of `data` given `lmb`.  A float for 1-D `data`,
        an array otherwise.

    See Also
    --------
    boxcox, probplot, boxcox_normplot, boxcox_normmax

    Notes
    -----
    The Box-Cox log-likelihood function is defined here as

    .. math::

        llf = (\lambda - 1) \sum_i(\log(x_i)) -
              N/2 \log(\sum_i (y_i - \bar{y})^2 / N),

    where ``y`` is the Box-Cox transformed input data ``x``.

    Examples
    --------
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt
    >>> from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    >>> np.random.seed(1245)

    Generate some random variates and calculate Box-Cox log-likelihood values
    for them for a range of ``lmbda`` values:

    >>> x = stats.loggamma.rvs(5, loc=10, size=1000)
    >>> lmbdas = np.linspace(-2, 10)
    >>> llf = np.zeros(lmbdas.shape, dtype=float)
    >>> for ii, lmbda in enumerate(lmbdas):
    ...     llf[ii] = stats.boxcox_llf(lmbda, x)

    Also find the optimal lmbda value with `boxcox`:

    >>> x_most_normal, lmbda_optimal = stats.boxcox(x)

    Plot the log-likelihood as function of lmbda.  Add the optimal lmbda as a
    horizontal line to check that that's really the optimum:

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(lmbdas, llf, 'b.-')
    >>> ax.axhline(stats.boxcox_llf(lmbda_optimal, x), color='r')
    >>> ax.set_xlabel('lmbda parameter')
    >>> ax.set_ylabel('Box-Cox log-likelihood')

    Now add some probability plots to show that where the log-likelihood is
    maximized the data transformed with `boxcox` looks closest to normal:

    >>> locs = [3, 10, 4]  # 'lower left', 'center', 'lower right'
    >>> for lmbda, loc in zip([-1, lmbda_optimal, 9], locs):
    ...     xt = stats.boxcox(x, lmbda=lmbda)
    ...     (osm, osr), (slope, intercept, r_sq) = stats.probplot(xt)
    ...     ax_inset = inset_axes(ax, width="20%", height="20%", loc=loc)
    ...     ax_inset.plot(osm, osr, 'c.', osm, slope*osm + intercept, 'k-')
    ...     ax_inset.set_xticklabels([])
    ...     ax_inset.set_yticklabels([])
    ...     ax_inset.set_title('$\lambda=%1.2f$' % lmbda)

    >>> plt.show()

    """
    data = np.asarray(data)
    N = data.shape[0]
    if N == 0:
        return np.nan

    y = boxcox(data, lmb)
    y_mean = np.mean(y, axis=0)
    llf = (lmb - 1) * np.sum(np.log(data), axis=0)
    llf -= N / 2.0 * np.log(np.sum((y - y_mean)**2. / N, axis=0))
    return llf


def _boxcox_conf_interval(x, lmax, alpha, method='lrt'):
    methods = {'lrt': _boxcox_conf_interval_lrt,
               'asy': _boxcox_conf_interval_asy}
    if method not in methods.keys():
        raise ValueError("Method %s not recognized." % method)
    return methods[method](x, lmax, alpha)


def _boxcox_conf_interval_asy(x, lmax, alpha):
    zstar = -special.ndtri(alpha / 2)
    x = _asarray_1d_positive(x)
    d1, d2 = _boxcox_llf_derivs(lmax, np.log(x))
    if d2 >= 0:
        raise RuntimeError('the second derivative of the log likelihood '
                           'at the mle should be negative')
    delta = zstar / np.sqrt(-d2)
    return (lmax - delta, lmax + delta)


def _boxcox_conf_interval_lrt(x, lmax, alpha):
    # Need to find the lambda for which
    #  f(x,lmbda) >= f(x,lmax) - 0.5*chi^2_alpha;1
    fac = 0.5 * distributions.chi2.ppf(1 - alpha, 1)
    target = boxcox_llf(lmax, x) - fac

    def rootfunc(lmbda, data, target):
        return boxcox_llf(lmbda, data) - target

    # Find positive endpoint of interval in which answer is to be found
    newlm = lmax + 0.5
    N = 0
    while (rootfunc(newlm, x, target) > 0.0) and (N < 500):
        newlm += 0.1
        N += 1

    if N == 500:
        raise RuntimeError("Could not find endpoint.")

    lmplus = optimize.brentq(rootfunc, lmax, newlm, args=(x, target))

    # Now find negative interval in the same way
    newlm = lmax - 0.5
    N = 0
    while (rootfunc(newlm, x, target) > 0.0) and (N < 500):
        newlm -= 0.1
        N += 1

    if N == 500:
        raise RuntimeError("Could not find endpoint.")

    lmminus = optimize.brentq(rootfunc, newlm, lmax, args=(x, target))
    return lmminus, lmplus


def boxcox(x, lmbda=None, alpha=None, lmbda_estimation=None, conf_method='lrt'):
    r"""
    Return a positive dataset transformed by a Box-Cox power transformation.

    Parameters
    ----------
    x : ndarray
        Input array.  Should be 1-dimensional.
    lmbda : {None, scalar}, optional
        If `lmbda` is not None, do the transformation for that value.

        If `lmbda` is None, find the lambda that maximizes the log-likelihood
        function and return it as the second output argument.
    alpha : {None, float}, optional
        If ``alpha`` is not None, return the ``100 * (1-alpha)%`` confidence
        interval for `lmbda` as the third output argument.
        Must be between 0.0 and 1.0.
    lmbda_estimation : str, optional
        One of {'pearsonr', 'mle', 'mle-newton'}.
    conf_method : str, optional
        One of {'lrt', 'asy'}.  This is only used if ``alpha`` is not None.
        Default: 'lrt'

    Returns
    -------
    boxcox : ndarray
        Box-Cox power transformed array.
    maxlog : float, optional
        If the `lmbda` parameter is None, the second returned argument is
        the lambda that maximizes the log-likelihood function.
    (min_ci, max_ci) : tuple of float, optional
        If `lmbda` parameter is None and ``alpha`` is not None, this returned
        tuple of floats represents the minimum and maximum confidence limits
        given ``alpha``.

    See Also
    --------
    probplot, boxcox_normplot, boxcox_normmax, boxcox_llf

    Notes
    -----
    The Box-Cox transform is given by::

        y = (x**lmbda - 1) / lmbda,  for lmbda > 0
            log(x),                  for lmbda = 0

    `boxcox` requires the input data to be positive.  Sometimes a Box-Cox
    transformation provides a shift parameter to achieve this; `boxcox` does
    not.  Such a shift parameter is equivalent to adding a positive constant to
    `x` before calling `boxcox`.

    The confidence limits returned when ``alpha`` is provided give the interval
    where:

    .. math::

        llf(\hat{\lambda}) - llf(\lambda) < \frac{1}{2}\chi^2(1 - \alpha, 1),

    with ``llf`` the log-likelihood function and :math:`\chi^2` the chi-squared
    function.

    References
    ----------
    G.E.P. Box and D.R. Cox, "An Analysis of Transformations", Journal of the
    Royal Statistical Society B, 26, 211-252 (1964).

    Examples
    --------
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    We generate some random variates from a non-normal distribution and make a
    probability plot for it, to show it is non-normal in the tails:

    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(211)
    >>> x = stats.loggamma.rvs(5, size=500) + 5
    >>> prob = stats.probplot(x, dist=stats.norm, plot=ax1)
    >>> ax1.set_xlabel('')
    >>> ax1.set_title('Probplot against normal distribution')

    We now use `boxcox` to transform the data so it's closest to normal:

    >>> ax2 = fig.add_subplot(212)
    >>> xt, _ = stats.boxcox(x)
    >>> prob = stats.probplot(xt, dist=stats.norm, plot=ax2)
    >>> ax2.set_title('Probplot after Box-Cox transformation')

    >>> plt.show()

    """
    x = np.asarray(x)
    if x.size == 0:
        return x

    if any(x <= 0):
        raise ValueError("Data must be positive.")

    if lmbda is not None and lmbda_estimation is None:
        return special.boxcox(x, lmbda)

    # If requested, find the lmbda that maximizes the log-likelihood function.
    brack = (-2, 2)
    if lmbda_estimation is not None:
        method = lmbda_estimation
        if lmbda is not None:
            brack = (lmbda-2, lmbda+2)
    else:
        method = 'mle'
    lmax = boxcox_normmax(x, brack=brack, method=method)
    y = boxcox(x, lmax)

    if alpha is None:
        return y, lmax
    else:
        # Find confidence interval
        interval = _boxcox_conf_interval(x, lmax, alpha, method=conf_method)
        return y, lmax, interval


def _boxcox_llf_derivs(lam, logy):
    # equations (6, 15) in the masters thesis
    # Note that z, u, v are centered here but not in the paper.
    # Here v is negated relative to v in the paper.
    #
    # See the differentiation rule DLMF 13.3.15:
    # d/dx 1f1(a, b, x) = a/b * 1f1(a+1, b+1, x)
    #
    # 1f1(1, 2, x) = (e^x - 1) / x (= exprel(x))
    # 1f1(2, 3, x) = 2*(x*e^x - e^x + 1) / x^2
    # 1f1(3, 4, x) = 3*(x^2*e^x -2*x*e^x + 2*e^x - 2) / x^3
    #
    # Note that exprel and the hypergeometric functions
    # remove the singularity at 0.

    def _boxcox_confluent_hypergeometric_helper(k):
        x = lam * logy
        r = special.hyp1f1(k, k+1, x) * (logy ** k) / k
        if not np.isfinite(r).all():
            msg = 'boxcox derivatives instability k=%d x=%s' % (k, x)
            warnings.warn(msg, RuntimeWarning)
        return r - r.mean()

    n = logy.size
    z = _boxcox_confluent_hypergeometric_helper(1)
    u = _boxcox_confluent_hypergeometric_helper(2)
    v = _boxcox_confluent_hypergeometric_helper(3)
    var = np.var(z)
    fprime = logy.sum() - u.dot(z) / var
    fprime2 = (2/n) * (u.dot(z) / var)**2 - (v.dot(z) + u.dot(u)) / var
    return fprime, fprime2


def _minimize_scalar_newton(x0, fprimes, xatol=1.48e-8, maxiter=50):
    # Related to scipy.optimize.zeros.newton().
    if xatol <= 0:
        raise ValueError('xatol too small (%g <= 0)' % xatol)
    p0 = x0
    for i in range(maxiter):
        d1, d2 = fprimes(p0)
        if d2 <= 0:
            msg = 'second derivative is not positive'
            warnings.warn(msg, RuntimeWarning)
        # Newton step.
        p = p0 - d1 / d2
        if abs(p - p0) < xatol:
            return p
        p0 = p
    msg = 'Failed to converge after %d iterations, value is %s' % (maxiter, p)
    raise RuntimeError(msg)



def boxcox_normmax(x, brack=(-2.0, 2.0), method='pearsonr'):
    """Compute optimal Box-Cox transform parameter for input data.

    Parameters
    ----------
    x : array_like
        Input array.
    brack : 2-tuple, optional
        The starting interval for a downhill bracket search with
        `optimize.brent`.  Note that this is in most cases not critical; the
        final result is allowed to be outside this bracket.
    method : str, optional
        The method to determine the optimal transform parameter (`boxcox`
        ``lmbda`` parameter). Options are:

        'pearsonr'  (default)
            Maximizes the Pearson correlation coefficient between
            ``y = boxcox(x)`` and the expected values for ``y`` if `x` would be
            normally-distributed.

        'mle'
            Maximizes the log-likelihood `boxcox_llf` using direct function
            evaluation.  This is the method used in `boxcox`.

        'mle-newton'
            Maximizes the log-likelihood `boxcox_llf` with the help of its
            first and second derivatives.

        'all'
            Use all optimization methods available, and return all results.
            Useful to compare different methods.

    Returns
    -------
    maxlog : float or ndarray
        The optimal transform parameter found.  An array instead of a scalar
        for ``method='all'``.

    See Also
    --------
    boxcox, boxcox_llf, boxcox_normplot

    Examples
    --------
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt
    >>> np.random.seed(1234)  # make this example reproducible

    Generate some data and determine optimal ``lmbda`` in various ways:

    >>> x = stats.loggamma.rvs(5, size=30) + 5
    >>> y, lmax_mle = stats.boxcox(x)
    >>> lmax_pearsonr = stats.boxcox_normmax(x)

    >>> lmax_mle
    7.177...
    >>> lmax_pearsonr
    7.916...
    >>> stats.boxcox_normmax(x, method='all')
    array([ 7.91667384,  7.17718692])

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> prob = stats.boxcox_normplot(x, -10, 10, plot=ax)
    >>> ax.axvline(lmax_mle, color='r')
    >>> ax.axvline(lmax_pearsonr, color='g', ls='--')

    >>> plt.show()

    """
    x = _asarray_1d_positive(x)

    def _pearsonr(x, brack):
        osm_uniform = _calc_uniform_order_statistic_medians(x)
        xvals = distributions.norm.ppf(osm_uniform)

        def _eval_pearsonr(lmbda, xvals, samps):
            # This function computes the x-axis values of the probability plot
            # and computes a linear regression (including the correlation) and
            # returns ``1 - r`` so that a minimization function maximizes the
            # correlation.
            y = boxcox(samps, lmbda)
            yvals = np.sort(y)
            r, prob = stats.pearsonr(xvals, yvals)
            return 1 - r

        return optimize.brent(_eval_pearsonr, brack=brack, args=(xvals, x))

    def _mle(x, brack):
        def _eval_mle(lmb, data):
            # function to minimize
            return -boxcox_llf(lmb, data)

        return optimize.brent(_eval_mle, brack=brack, args=(x,))

    def _all(x, brack):
        maxlog = np.zeros(3, dtype=float)
        maxlog[0] = _pearsonr(x, brack)
        maxlog[1] = _mle(x, brack)
        maxlog[2] = _mle_newton(x, brack)
        return maxlog

    def _mle_newton(x, brack):
        lam0 = np.mean(brack)
        logx = np.log(x)

        def _fprimes(lam):
            nd1, nd2 = _boxcox_llf_derivs(lam, logx)
            return -nd1, -nd2

        return _minimize_scalar_newton(lam0, _fprimes)

    methods = {'pearsonr': _pearsonr,
               'mle': _mle,
               'mle-newton': _mle_newton,
               'all': _all}
    if method not in methods.keys():
        raise ValueError("Method %s not recognized." % method)

    optimfunc = methods[method]
    return optimfunc(x, brack)
