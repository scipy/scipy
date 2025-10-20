import warnings
import numpy as np
from . import distributions
from .._lib._array_api import xp_capabilities, xp_promote
from .._lib._bunch import _make_tuple_bunch
from ._axis_nan_policy import _axis_nan_policy_factory
from ._stats_pythran import siegelslopes as siegelslopes_pythran
import scipy.stats._stats_py as _stats_py

__all__ = ['_find_repeats', 'theilslopes', 'siegelslopes']

# This is not a namedtuple for backwards compatibility. See PR #12983
TheilslopesResult = _make_tuple_bunch('TheilslopesResult',
                                      ['slope', 'intercept',
                                       'low_slope', 'high_slope'])
SiegelslopesResult = _make_tuple_bunch('SiegelslopesResult',
                                       ['slope', 'intercept'])


def _n_samples_optional_x(kwargs):
    return 2 if kwargs.get('x', None) is not None else 1


@xp_capabilities(np_only=True)
@_axis_nan_policy_factory(TheilslopesResult, default_axis=None, n_outputs=4,
                          n_samples=_n_samples_optional_x,
                          result_to_tuple=lambda x, _: tuple(x), paired=True,
                          too_small=1)
def theilslopes(y, x=None, alpha=0.95, method='separate', *, axis=None):
    r"""
    Computes the Theil-Sen estimator for a set of points (x, y).

    `theilslopes` implements a method for robust linear regression.  It
    computes the slope as the median of all slopes between paired values.

    Parameters
    ----------
    y : array_like
        Dependent variable.
    x : array_like or None, optional
        Independent variable. If None, use ``arange(len(y))`` instead.
    alpha : float, optional
        Confidence degree between 0 and 1. Default is 95% confidence.
        Note that `alpha` is symmetric around 0.5, i.e. both 0.1 and 0.9 are
        interpreted as "find the 90% confidence interval".
    method : {'joint', 'separate'}, optional
        Method to be used for computing estimate for intercept.
        Following methods are supported,

            * 'joint': Uses np.median(y - slope * x) as intercept.
            * 'separate': Uses np.median(y) - slope * np.median(x)
                          as intercept.

        The default is 'separate'.

        .. versionadded:: 1.8.0

    axis : int or tuple of ints, default: None
        If an int or tuple of ints, the axis or axes of the input along which
        to compute the statistic. The statistic of each axis-slice (e.g. row)
        of the input will appear in a corresponding element of the output.
        If ``None``, the input will be raveled before computing the statistic.

    Returns
    -------
    result : ``TheilslopesResult`` instance
        The return value is an object with the following attributes:

        slope : float
            Theil slope.
        intercept : float
            Intercept of the Theil line.
        low_slope : float
            Lower bound of the confidence interval on `slope`.
        high_slope : float
            Upper bound of the confidence interval on `slope`.

    See Also
    --------
    siegelslopes : a similar technique using repeated medians

    Notes
    -----
    The implementation of `theilslopes` follows [1]_. The intercept is
    not defined in [1]_, and here it is defined as ``median(y) -
    slope*median(x)``, which is given in [3]_. Other definitions of
    the intercept exist in the literature such as  ``median(y - slope*x)``
    in [4]_. The approach to compute the intercept can be determined by the
    parameter ``method``. A confidence interval for the intercept is not
    given as this question is not addressed in [1]_.

    For compatibility with older versions of SciPy, the return value acts
    like a ``namedtuple`` of length 4, with fields ``slope``, ``intercept``,
    ``low_slope``, and ``high_slope``, so one can continue to write::

        slope, intercept, low_slope, high_slope = theilslopes(y, x)

    References
    ----------
    .. [1] P.K. Sen, "Estimates of the regression coefficient based on
           Kendall's tau", J. Am. Stat. Assoc., Vol. 63, pp. 1379-1389, 1968.
    .. [2] H. Theil, "A rank-invariant method of linear and polynomial
           regression analysis I, II and III",  Nederl. Akad. Wetensch., Proc.
           53:, pp. 386-392, pp. 521-525, pp. 1397-1412, 1950.
    .. [3] W.L. Conover, "Practical nonparametric statistics", 2nd ed.,
           John Wiley and Sons, New York, pp. 493.
    .. [4] https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-5, 5, num=150)
    >>> y = x + np.random.normal(size=x.size)
    >>> y[11:15] += 10  # add outliers
    >>> y[-5:] -= 7

    Compute the slope, intercept and 90% confidence interval.  For comparison,
    also compute the least-squares fit with `linregress`:

    >>> res = stats.theilslopes(y, x, 0.90, method='separate')
    >>> lsq_res = stats.linregress(x, y)

    Plot the results. The Theil-Sen regression line is shown in red, with the
    dashed red lines illustrating the confidence interval of the slope (note
    that the dashed red lines are not the confidence interval of the regression
    as the confidence interval of the intercept is not included). The green
    line shows the least-squares fit for comparison.

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, y, 'b.')
    >>> ax.plot(x, res[1] + res[0] * x, 'r-')
    >>> ax.plot(x, res[1] + res[2] * x, 'r--')
    >>> ax.plot(x, res[1] + res[3] * x, 'r--')
    >>> ax.plot(x, lsq_res[1] + lsq_res[0] * x, 'g-')
    >>> plt.show()

    """
    if method not in ['joint', 'separate']:
        raise ValueError("method must be either 'joint' or 'separate'."
                         f"'{method}' is invalid.")

    y, x = xp_promote(y, x, force_floating=True, xp=np)
    x = np.arange(y.shape[-1], dtype=y.dtype) if x is None else x
    y, x = np.broadcast_arrays(y, x)

    # Compute sorted slopes only when deltax > 0
    deltax = x[..., :, np.newaxis] - x[..., np.newaxis, :]
    deltay = y[..., :, np.newaxis] - y[..., np.newaxis, :]
    i = np.triu(np.ones(deltax.shape[-2:], dtype=bool), k=1)
    deltax = np.reshape(deltax, deltax.shape[:-2] + (-1,))
    deltay = np.reshape(deltay, deltay.shape[:-2] + (-1,))
    i = np.ravel(i)
    deltax, deltay = deltax[..., i], deltay[..., i]
    deltax[deltax == 0] = np.nan
    slopes = deltay / deltax
    slopes = np.sort(slopes, axis=-1)
    medslope = np.nanmedian(slopes, axis=-1)
    if method == 'joint':
        medinter = np.nanmedian(y - medslope * x, axis=-1)
    else:
        medinter = np.nanmedian(y, axis=-1) - medslope * np.nanmedian(x, axis=-1)
    # Now compute confidence intervals
    if alpha > 0.5:
        alpha = 1. - alpha

    z = distributions.norm.ppf(alpha / 2.)
    # This implements (2.6) from Sen (1968)
    _, nxreps = _stats_py._rankdata(x, method='average', return_ties=True)
    _, nyreps = _stats_py._rankdata(y, method='average', return_ties=True)
    nt = np.count_nonzero(np.isfinite(slopes), axis=-1, keepdims=True)  # N in Sen (1968)
    ny = y.shape[-1]                                                    # n in Sen (1968)
    # Equation 2.6 in Sen (1968):
    sigsq = 1/18. * (
        ny * (ny-1) * (2*ny+5)
        - np.sum(nxreps * (nxreps-1) * (2*nxreps + 5), axis=-1, keepdims=True)
        - np.sum(nyreps * (nyreps-1) * (2*nyreps + 5), axis=-1, keepdims=True))
    # Find the confidence interval indices in `slopes`
    sigma = np.sqrt(sigsq)
    Ru = np.minimum(np.astype(np.round((nt - z*sigma)/2.), int), nt-1)
    Rl = np.maximum(np.astype(np.round((nt + z*sigma)/2.), int) - 1, 0)
    R = np.concatenate((np.atleast_1d(Rl), np.atleast_1d(Ru)), axis=-1)
    delta = np.take_along_axis(slopes, R, axis=-1)
    i_nan = np.broadcast_to(sigsq < 0, delta.shape)
    delta[i_nan] = np.nan

    return TheilslopesResult(slope=medslope[()], intercept=medinter[()],
                             low_slope=delta[..., 0][()], high_slope=delta[..., 1][()])


def _find_repeats(arr):
    # This function assumes it may clobber its input.
    if len(arr) == 0:
        return np.array(0, np.float64), np.array(0, np.intp)

    # XXX This cast was previously needed for the Fortran implementation,
    # should we ditch it?
    arr = np.asarray(arr, np.float64).ravel()
    arr.sort()

    # Taken from NumPy 1.9's np.unique.
    change = np.concatenate(([True], arr[1:] != arr[:-1]))
    unique = arr[change]
    change_idx = np.concatenate(np.nonzero(change) + ([arr.size],))
    freq = np.diff(change_idx)
    atleast2 = freq > 1
    return unique[atleast2], freq[atleast2]


@xp_capabilities(np_only=True)
@_axis_nan_policy_factory(SiegelslopesResult, default_axis=None, n_outputs=2,
                          n_samples=_n_samples_optional_x,
                          result_to_tuple=lambda x, _: tuple(x), paired=True,
                          too_small=1)
def siegelslopes(y, x=None, method="hierarchical"):
    r"""
    Computes the Siegel estimator for a set of points (x, y).

    `siegelslopes` implements a method for robust linear regression
    using repeated medians (see [1]_) to fit a line to the points (x, y).
    The method is robust to outliers with an asymptotic breakdown point
    of 50%.

    Parameters
    ----------
    y : array_like
        Dependent variable.
    x : array_like or None, optional
        Independent variable. If None, use ``arange(len(y))`` instead.
    method : {'hierarchical', 'separate'}
        If 'hierarchical', estimate the intercept using the estimated
        slope ``slope`` (default option).
        If 'separate', estimate the intercept independent of the estimated
        slope. See Notes for details.

    Returns
    -------
    result : ``SiegelslopesResult`` instance
        The return value is an object with the following attributes:

        slope : float
            Estimate of the slope of the regression line.
        intercept : float
            Estimate of the intercept of the regression line.

    See Also
    --------
    theilslopes : a similar technique without repeated medians

    Notes
    -----
    With ``n = len(y)``, compute ``m_j`` as the median of
    the slopes from the point ``(x[j], y[j])`` to all other `n-1` points.
    ``slope`` is then the median of all slopes ``m_j``.
    Two ways are given to estimate the intercept in [1]_ which can be chosen
    via the parameter ``method``.
    The hierarchical approach uses the estimated slope ``slope``
    and computes ``intercept`` as the median of ``y - slope*x``.
    The other approach estimates the intercept separately as follows: for
    each point ``(x[j], y[j])``, compute the intercepts of all the `n-1`
    lines through the remaining points and take the median ``i_j``.
    ``intercept`` is the median of the ``i_j``.

    The implementation computes `n` times the median of a vector of size `n`
    which can be slow for large vectors. There are more efficient algorithms
    (see [2]_) which are not implemented here.

    For compatibility with older versions of SciPy, the return value acts
    like a ``namedtuple`` of length 2, with fields ``slope`` and
    ``intercept``, so one can continue to write::

        slope, intercept = siegelslopes(y, x)

    References
    ----------
    .. [1] A. Siegel, "Robust Regression Using Repeated Medians",
           Biometrika, Vol. 69, pp. 242-244, 1982.

    .. [2] A. Stein and M. Werman, "Finding the repeated median regression
           line", Proceedings of the Third Annual ACM-SIAM Symposium on
           Discrete Algorithms, pp. 409-413, 1992.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-5, 5, num=150)
    >>> y = x + np.random.normal(size=x.size)
    >>> y[11:15] += 10  # add outliers
    >>> y[-5:] -= 7

    Compute the slope and intercept.  For comparison, also compute the
    least-squares fit with `linregress`:

    >>> res = stats.siegelslopes(y, x)
    >>> lsq_res = stats.linregress(x, y)

    Plot the results. The Siegel regression line is shown in red. The green
    line shows the least-squares fit for comparison.

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, y, 'b.')
    >>> ax.plot(x, res[1] + res[0] * x, 'r-')
    >>> ax.plot(x, lsq_res[1] + lsq_res[0] * x, 'g-')
    >>> plt.show()

    """
    if method not in ['hierarchical', 'separate']:
        raise ValueError("method can only be 'hierarchical' or 'separate'")
    y = np.asarray(y).ravel()
    if x is None:
        x = np.arange(len(y), dtype=float)
    else:
        x = np.asarray(x, dtype=float).ravel()
        if len(x) != len(y):
            raise ValueError("Array shapes are incompatible for broadcasting.")
    if len(x) < 2:
        raise ValueError("`x` and `y` must have length at least 2.")

    dtype = np.result_type(x, y, np.float32)  # use at least float32
    y, x = y.astype(dtype), x.astype(dtype)
    medslope, medinter = siegelslopes_pythran(y, x, method)
    medslope, medinter = np.asarray(medslope)[()], np.asarray(medinter)[()]
    return SiegelslopesResult(slope=medslope, intercept=medinter)
