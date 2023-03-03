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
        uncensored and right-censored observations are supported. Currently,
        other instances of `stats.CensoredData` will result in a
        ``NotImplementedError``.

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
    steps up by ``1/len(sample)`` at each of the observations [1]_.

    When observations are lower bounds, upper bounds, or both upper and lower
    bounds, the data is said to be "censored", and `sample` may be provided as
    an instance of `stats.CensoredData`.

    For right-censored data, the ECDF is given by the Kaplan-Meier estimator
    [2]_; other forms of censoring are not supported at this time.

    References
    ----------
    .. [1] Conover, William Jay. Practical nonparametric statistics. Vol. 350.
           John Wiley & Sons, 1999.

    .. [2] Kaplan, Edward L., and Paul Meier. "Nonparametric estimation from
           incomplete observations." Journal of the American statistical
           association 53.282 (1958): 457-481.

    .. [3] Goel, Manish Kumar, Pardeep Khanna, and Jugal Kishore.
           "Understanding survival analysis: Kaplan-Meier estimate."
           International journal of Ayurveda research 1.4 (2010): 274.

    Examples
    --------
    **Uncensored Data**

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

    **Right-censored Data**

    As in the example from [1]_ page 91, the lives of ten car fanbelts were
    tested. At the end of the test, five fanbelts had broken, but five were
    still unbroken. Their survival times were recorded as follows.

    >>> broken = [77, 47, 81, 56, 80]  # in thousands of miles
    >>> unbroken = [62, 60, 43, 71, 37]

    The empirical survival function is calculated as follows.

    >>> sample = stats.CensoredData(uncensored=broken, right=unbroken)
    >>> res = stats.ecdf(sample)
    >>> res.x
    array([37., 43., 47., 56., 60., 62., 71., 77., 80., 81.])
    >>> res.sf
    array([1.   , 1.   , 0.875, 0.75 , 0.75 , 0.75 , 0.75 , 0.5  , 0.25 , 0.   ])

    To plot the result as a step function:

    >>> ax = plt.subplot()
    >>> ax.step(np.insert(res.x, 0, 30), np.insert(res.sf, 0, 1), where='post')
    >>> ax.set_xlabel('Fanbelt Survival Time (thousands of miles)')
    >>> ax.set_ylabel('Empirical SF')
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
    elif sample.num_censored() == sample._right.size:
        res = _ecdf_right_censored(sample)
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


def _ecdf_right_censored(sample):
    # It is conventional to discuss right-censored data in terms of
    # "survival time", "death", and "loss" (e.g. [2]). We'll use that
    # terminology here.
    # This implementation was influenced by the references cited and also
    # https://www.youtube.com/watch?v=lxoWsVco_iM
    # https://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator
    # In retrospect it is probably most easily compared against [3].
    # Ultimately, the data needs to be sorted, so this implementation is
    # written to avoid a separate call to `unique` after sorting. In hope of
    # better performance on large datasets, it also computes survival
    # probabilities at unique times only rather than at each observation.
    tod = sample._uncensored  # time of "death"
    tol = sample._right  # time of "loss"
    times = np.concatenate((tod, tol))
    died = np.asarray([1]*tod.size + [0]*tol.size)

    # sort by times
    i = np.argsort(times)
    times = times[i]
    died = died[i]
    at_risk = np.arange(times.size, 0, -1)

    # logical indices of unique times
    j = np.diff(times, prepend=-np.inf, append=np.inf) > 0
    j_l = j[:-1]  # first instances of unique times
    j_r = j[1:]  # last instances of unique times

    # get number at risk and deaths at each unique time
    t = times[j_l]  # unique times
    n = at_risk[j_l]  # number at risk at each unique time
    cd = np.cumsum(died)[j_r]  # cumulative deaths up to/including unique times
    d = np.diff(cd, prepend=0)  # deaths at each unique time

    # compute survival function
    sf = np.cumprod((n - d) / n)
    cdf = 1 - sf
    return ECDFResult(t, cdf, sf)
