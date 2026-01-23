import numpy as np
from .._lib._array_api import xp_capabilities
from .._lib._bunch import _make_tuple_bunch
from ._axis_nan_policy import _axis_nan_policy_factory
from ._stats_pythran import siegelslopes as siegelslopes_pythran

__all__ = ['_find_repeats', 'siegelslopes']

# This is not a namedtuple for backwards compatibility. See PR #12983
TheilslopesResult = _make_tuple_bunch('TheilslopesResult',
                                      ['slope', 'intercept',
                                       'low_slope', 'high_slope'])
SiegelslopesResult = _make_tuple_bunch('SiegelslopesResult',
                                       ['slope', 'intercept'])


def _n_samples_optional_x(kwargs):
    return 2 if kwargs.get('x', None) is not None else 1


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
