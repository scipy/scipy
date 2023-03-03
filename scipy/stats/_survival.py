from dataclasses import make_dataclass
import numpy as np
from ._censored_data import CensoredData


__all__ = ['ecdf']


ECDFResult = make_dataclass('ECDFResult', ['x', 'cdf', 'sf'])


def ecdf(sample):
    """Empirical cumulative distribution function of a sample.

    The empirical cumulative distribution function (ECDF) is a step function
    estimate of the CDF of the distribution underlying a sample.

    Parameters
    ----------
    sample : 1D array_like or `stats.CensoredData`
        Besides array_like, instances of `stats.CensoredData` containing
        uncensored observations are supported. Currently, instances of
        `stats.CensoredData` with censored data will result in a
        ``NotImplementedError``, but future support for left-censored,
        right-centered, and interval-censored data is planned.

    Returns
    -------
    An object with the following attributes.

    x : ndarray
        The unique values at which the ECDF changes.
    cdf : ndarray
        The values of the ECDF corresponding with `x`.
    sf : ndarray
        The empirical survival function, the complement of the ECDF.

    Notes
    -----
    When each observation of the sample is a precise measurement, the ECDF
    steps up by ``1/len(sample)`` at each of the observations.

    When observations are lower bounds, upper bounds, or both upper and lower
    bounds, the data is said to be "censored", and `sample` may be provided as
    an instance of `stats.CensoredData`.

    For right-censored data, the ECDF is given by the Kaplan-Meier estimator
    [1]_; other forms of censoring are not supported at this time.

    References
    ----------
    .. [1] Conover, William Jay. Practical nonparametric statistics. Vol. 350.
           John Wiley & Sons, 1999.

    .. [2] Kaplan, Edward L., and Paul Meier. "Nonparametric estimation from
           incomplete observations." Journal of the American statistical
           association 53.282 (1958): 457-481.

    Examples
    --------
    As in the example from [1]_ page 79, five boys were selected at random from
    those in a single high school. Their one-mile run times were recorded as
    follows.

    >>> sample = [6.23, 5.58, 7.06, 6.42, 5.20]  # one-mile run times (minutes)

    The empirical distribution function, which approximates the distribution
    function of one-mile run times of the population from which the boys were
    sampled, is calculated as follows.

    >>> from scipy import stats
    >>> res = stats.ecdf(sample)
    >>> res.x
    array([5.2 , 5.58, 6.23, 6.42, 7.06])
    >>> res.cdf
    array([0.2, 0.4, 0.6, 0.8, 1. ])

    To plot the result as a step function:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> ax = plt.subplot()
    >>> ax.step(np.insert(res.x, 0, 4), np.insert(res.cdf, 0, 0), where='post')
    >>> ax.set_xlabel('One-Mile Run Time (minutes)')
    >>> ax.set_ylabel('Empirical CDF')
    >>> plt.show()

    """
    if not isinstance(sample, CensoredData):
        try:  # takes care of input standardization/validation
            sample = CensoredData(uncensored=sample)
        except ValueError as e:
            message = str(e).replace('uncensored', 'sample')
            raise type(e)(message) from e

    if sample.num_censored() == 0:
        res = _ecdf_uncensored(sample._uncensor())
    else:
        # Support censoring in follow-up PRs
        message = ("Currently, only uncensored data is supported.")
        raise NotImplementedError(message)
    return res


def _ecdf_uncensored(sample):
    sample = np.sort(sample)
    x, counts = np.unique(sample, return_counts=True)

    # [1].81 "the fraction of [observations] that are less than or equal to x
    cdf = np.cumsum(counts) / sample.size

    # [1].89 "the relative frequency of the sample that exceeds x in value"
    sf = 1 - cdf

    return ECDFResult(x, cdf, sf)
