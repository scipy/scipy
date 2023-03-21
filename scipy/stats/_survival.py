import warnings
from dataclasses import dataclass, field
import numpy as np
from scipy import special
from scipy.stats._censored_data import CensoredData
from scipy.stats._common import ConfidenceInterval


__all__ = ['ecdf']


@dataclass
class EmpiricalDistributionFunction:
    """An empirical distribution function produced by `scipy.stats.ecdf`

    Attributes
    ----------
    points : ndarray
        The point estimate of the cumulative distribution function (CDF) or its
        complement, the survival function (SF), at unique values of the sample.
    """
    points: np.ndarray
    # Exclude these from __str__
    _x: np.ndarray = field(repr=False)  # points at which function is estimated
    _n: np.ndarray = field(repr=False)  # number "at risk"
    _d: np.ndarray = field(repr=False)  # number of "deaths"
    _sf: np.ndarray = field(repr=False)  # survival function for var estimate
    _kind: str = field(repr=False)  # type of function: "cdf" or "sf"

    def __init__(self, x, points, n, d, kind):
        self.points = points
        self._x = x
        self._n = n
        self._d = d
        self._sf = points if kind == 'sf' else 1 - points
        self._kind = kind

    def confidence_interval(self, confidence_level=0.95, *, method='linear'):
        """Compute a confidence interval around the CDF/SF point estimate

        Parameters
        ----------
        confidence_level : float, default: 0.95
            Confidence level for the computed confidence interval

        method : str, {"linear", "log-log"}
            Method used to compute the confidence interval. Options are
            "linear" for the conventional Greenwood confidence interval
            (default)  and "log-log" for the "exponential Greenwood",
            log-negative-log-transformed confidence interval.

        Returns
        -------
        ci : ``ConfidenceInterval``
            An object with attributes ``low`` and ``high``: arrays of the
            lower and upper bounds of the confidence interval at unique values
            of the sample.

        Notes
        -----
        Confidence intervals are computed according to the Greenwood formula
        (``method='linear'``) or the more recent "exponential Greenwood"
        formula (``method='log-log'``) as described in [1]_. The conventional
        Greenwood formula can result in lower confidence limits less than 0
        and upper confidence limits greater than 1; these are clipped to the
        unit interval. NaNs may be produced by either method; these are
        features of the formulas.

        References
        ----------
        .. [1] Sawyer, Stanley. "The Greenwood and Exponential Greenwood
               Confidence Intervals in Survival Analysis."
               https://www.math.wustl.edu/~sawyer/handouts/greenwood.pdf

        """
        methods = {'linear': self._linear_ci,
                   'log-log': self._loglog_ci}

        message = f"`method` must be one of {set(methods)}."
        if method.lower() not in methods:
            raise ValueError(message)

        message = "`confidence_level` must be a scalar between 0 and 1."
        confidence_level = np.asarray(confidence_level)[()]
        if confidence_level.shape or not (0 <= confidence_level <= 1):
            raise ValueError(message)

        method_fun = methods[method.lower()]
        low, high = method_fun(confidence_level)

        message = ("The confidence interval is undefined at some observations."
                   " This is a feature of the mathematical formula used, not"
                   " an error in its implementation.")
        if np.any(np.isnan(low) | np.isnan(high)):
            warnings.warn(message, RuntimeWarning, stacklevel=2)

        return ConfidenceInterval(np.clip(low, 0, 1), np.clip(high, 0, 1))

    def _linear_ci(self, confidence_level):
        sf, d, n = self._sf, self._d, self._n
        # When n == d, Greenwood's formula divides by zero.
        # When s != 0, this can be ignored: var == inf, and CI is [0, 1]
        # When s == 0, this results in NaNs. Produce an informative warning.
        with np.errstate(divide='ignore', invalid='ignore'):
            var = sf ** 2 * np.cumsum(d / (n * (n - d)))

        se = np.sqrt(var)
        z = special.ndtri(1 / 2 + confidence_level / 2)

        z_se = z * se
        low = self.points - z_se
        high = self.points + z_se

        return low, high

    def _loglog_ci(self, confidence_level):
        sf, d, n = self._sf, self._d, self._n

        with np.errstate(divide='ignore', invalid='ignore'):
            var = 1 / np.log(sf) ** 2 * np.cumsum(d / (n * (n - d)))

        se = np.sqrt(var)
        z = special.ndtri(1 / 2 + confidence_level / 2)

        with np.errstate(divide='ignore'):
            lnl_points = np.log(-np.log(sf))

        z_se = z * se
        low = np.exp(-np.exp(lnl_points + z_se))
        high = np.exp(-np.exp(lnl_points - z_se))
        if self._kind == "cdf":
            low, high = 1-high, 1-low

        return low, high


@dataclass
class ECDFResult:
    """ Result object returned by `scipy.stats.ecdf`

    Attributes
    ----------
    x : ndarray
        The unique values of the sample processed by `scipy.stats.ecdf`.
    cdf : `~scipy.stats._result_classes.EmpiricalDistributionFunction`
        An object representing the empirical cumulative distribution function.
    sf : `~scipy.stats._result_classes.EmpiricalDistributionFunction`
        An object representing the complement of the empirical cumulative
        distribution function.
    """
    x: np.ndarray
    cdf: EmpiricalDistributionFunction
    sf: EmpiricalDistributionFunction

    def __init__(self, x, cdf, sf, n, d):
        self.x = x
        # Both CDF and SF results need SF for variance est.
        self.cdf = EmpiricalDistributionFunction(x, cdf, n, d, "cdf")
        self.sf = EmpiricalDistributionFunction(x, sf, n, d, "sf")


def ecdf(sample):
    """Empirical cumulative distribution function of a sample.

    The empirical cumulative distribution function (ECDF) is a step function
    estimate of the CDF of the distribution underlying a sample. This function
    returns objects representing both the empirical distribution function and
    its complement, the empirical survival function.

    Parameters
    ----------
    sample : 1D array_like or `scipy.stats.CensoredData`
        Besides array_like, instances of `scipy.stats.CensoredData` containing
        uncensored and right-censored observations are supported. Currently,
        other instances of `scipy.stats.CensoredData` will result in a
        ``NotImplementedError``.

    Returns
    -------
    res : `~scipy.stats._result_classes.ECDFResult`
        An object with the following attributes.

        x : ndarray
            The unique values in the sample.
        cdf : `~scipy.stats._result_classes.EmpiricalDistributionFunction`
            An object representing the empirical cumulative distribution function.
        sf : `~scipy.stats._result_classes.EmpiricalDistributionFunction`
            An object representing the empirical survival function.

        The `cdf` and `sf` attributes themselves have the following attributes.

        points : ndarray
            The point estimate of the CDF/SF at the values in `x`.

        And the following method:

        confidence_interval(confidence_level=0.95) :
            Compute the confidence interval around the CDF/SF at the values in
            `x`.

    Notes
    -----
    When each observation of the sample is a precise measurement, the ECDF
    steps up by ``1/len(sample)`` at each of the observations [1]_.

    When observations are lower bounds, upper bounds, or both upper and lower
    bounds, the data is said to be "censored", and `sample` may be provided as
    an instance of `scipy.stats.CensoredData`.

    For right-censored data, the ECDF is given by the Kaplan-Meier estimator
    [2]_; other forms of censoring are not supported at this time.

    Confidence intervals are computed according to the Greenwood formula or the
    more recent "Exponential Greenwood" formula as described in [4]_.

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

    .. [4] Sawyer, Stanley. "The Greenwood and Exponential Greenwood Confidence
           Intervals in Survival Analysis."
           https://www.math.wustl.edu/~sawyer/handouts/greenwood.pdf

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
    >>> res.cdf.points
    array([0.2, 0.4, 0.6, 0.8, 1. ])

    To plot the result as a step function:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> ax = plt.subplot()
    >>> ax.step(np.insert(res.x, 0, 4), np.insert(res.cdf.points, 0, 0),
    ...         where='post')
    >>> ax.set_xlabel('One-Mile Run Time (minutes)')
    >>> ax.set_ylabel('Empirical CDF')
    >>> plt.show()

    **Right-censored Data**

    As in the example from [1]_ page 91, the lives of ten car fanbelts were
    tested. Five tests concluded because the fanbelt being tested broke, but
    the remaining tests concluded for other reasons (e.g. the study ran out of
    funding, but the fanbelt was still functional). The mileage driven
    with the fanbelts were recorded as follows.

    >>> broken = [77, 47, 81, 56, 80]  # in thousands of miles driven
    >>> unbroken = [62, 60, 43, 71, 37]

    Precise survival times of the fanbelts that were still functional at the
    end of the tests are unknown, but they are known to exceed the values
    recorded in ``unbroken``. Therefore, these observations are said to be
    "right-censored", and the data is represented using
    `scipy.stats.CensoredData`.

    >>> sample = stats.CensoredData(uncensored=broken, right=unbroken)

    The empirical survival function is calculated as follows.

    >>> res = stats.ecdf(sample)
    >>> res.x
    array([37., 43., 47., 56., 60., 62., 71., 77., 80., 81.])
    >>> res.sf.points
    array([1.   , 1.   , 0.875, 0.75 , 0.75 , 0.75 , 0.75 , 0.5  , 0.25 , 0.   ])

    To plot the result as a step function:

    >>> ax = plt.subplot()
    >>> ax.step(np.insert(res.x, 0, 30), np.insert(res.sf.points, 0, 1),
    ...         where='post')
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
        # Support additional censoring options in follow-up PRs
        message = ("Currently, only uncensored and right-censored data is "
                   "supported.")
        raise NotImplementedError(message)

    t, cdf, sf, n, d = res
    return ECDFResult(t, cdf, sf, n, d)


def _ecdf_uncensored(sample):
    sample = np.sort(sample)
    x, counts = np.unique(sample, return_counts=True)

    # [1].81 "the fraction of [observations] that are less than or equal to x
    cdf = np.cumsum(counts) / sample.size

    # [1].89 "the relative frequency of the sample that exceeds x in value"
    sf = 1 - cdf

    return x, cdf, sf, sf * sample.size, counts


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
    return t, cdf, sf, n, d
