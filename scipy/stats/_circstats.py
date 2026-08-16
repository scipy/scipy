import math
import scipy._external.array_api_extra as xpx
from scipy._lib._array_api import array_namespace, xp_promote, xp_capabilities
from scipy.stats._axis_nan_policy import _axis_nan_policy_factory, _broadcast_arrays
from scipy.stats._stats_py import _xp_mean


__all__ = ['circmedian']


def _circmeandev(sample, alpha, *, high=2*math.pi, low=0):
    r"""Compute the circular mean (absolute) deviation of an angular sample.

    Given :math:`n` angle observations :math:`x_1, \cdots, x_n` measured in
    radians, the *circular mean (absolute) deviation* about a given angle
    :math:`\alpha` is defined by ([1]_, Eq. 2.3.14)

    .. math::

       d_0(\alpha) = \frac{1}{n} \sum_{i=1}^n
                                \min \left(x_i - \alpha, 2 \pi - (x_i - \alpha) \right).

    For greater generality, this function computes:

    .. math::

       d_0(\alpha) = \frac{1}{n} \sum_{i=1}^n \min(y_i,  T - y_i),

    where :math:`y_i = (x_i - \alpha) \mod T` and :math:`T` is the period of a full
    revolution (e.g. :math:`360` when working in degrees).

    Parameters
    ----------
    sample : array_like
        Input array of angle observations.
    alpha : array_like
        The reference angle with respect to which the mean deviation is computed.
    high : float, optional
        Upper boundary of the principal value of an angle.  Default is ``2*pi``.
    low : float, optional
        Lower boundary of the principal value of an angle.  Default is ``0``.

    Returns
    -------
    circmeandev : float
        Circular mean (absolute) deviation, restricted to the range ``[low, high]``.

    See Also
    --------
    circmean : Circular mean.
    circstd : Circular standard deviation.
    circvar : Circular variance.

    References
    ----------
    .. [1] Mardia, K. V. and Jupp, P. E. *Directional Statistics*.
           John Wiley & Sons, 1999.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats._circstats import _circmeandev
    >>> import matplotlib.pyplot as plt
    >>> sample = np.array([20, 30, 330])
    >>> dev = _circmeandev(sample, sample, high=360, low=0)
    array([20.        , 23.33333333, 36.66666667])

    Among all the observations, the circular mean deviation is minimized about 20
    degrees, so 20 degrees is the circular median.

    """
    xp = array_namespace(sample)
    # assumes we're working along the last axis
    sample, alpha = _broadcast_arrays((sample, alpha), axis=-1, xp=xp)
    sample = sample[..., xp.newaxis, :]
    alpha = alpha[..., :, xp.newaxis]
    T = high - low
    y = (sample - alpha) % T
    return _xp_mean(xp.minimum(y, T - y), axis=-1)


@xp_capabilities()
@_axis_nan_policy_factory(lambda x: x, n_outputs=1, default_axis=None,
                          result_to_tuple=lambda x, _: (x,))
def circmedian(sample, *, convention='arc-distance', high=2*math.pi, low=0, axis=0):
    r"""Compute the circular median of an angular sample.

    According to [1]_ and [2]_, a *circular median* is an angle that bisects the data:
    half of the observations lie within 180 degrees clockwise and half lie within 180
    degrees counter-clockwise. The implementation of `circmedian` agrees with these
    references but follows the details of [3]_, which defines this "bisecting property"
    more rigorously with consideration for edge cases (e.g. ties, symmetry).

    Parameters
    ----------
    sample : array_like
        Input array of angle observations.
    convention : {'arc-distance', 'bisecting'}
        Definition of the circular median, following the terminology of [3]_.

        - An ``'arc-distance'`` median minimizes the circular mean deviation: the
          average arc distance between the observations and a common reference angle.
          Any arc-distance median also has the "bisecting property": the chord between
          a median and its antipode divide the circle in two such that at least
          half the observations lie in each semicircle.
        - A ``'bisecting'`` median has the bisecting property, but does not necessarily
          minimize the circular mean deviation. Instead, the number of observations
          within 90 degrees of a bisecting median is greater than the number of
          observations within 90 degrees of its antipode.

        See [3]_ for precise mathematical definitions.

    high : float, optional
        Upper boundary of the principal value of an angle.  Default is ``2*pi``.
    low : float, optional
        Lower boundary of the principal value of an angle.  Default is ``0``.
    axis : int, default: 0
        Axis along which the circular median is computed.

    Returns
    -------
    circmedian : float
        Circular median, restricted to the range ``[low, high]``.

    See Also
    --------
    circmean : Circular mean.
    circstd : Circular standard deviation.
    circvar : Circular variance.

    Notes
    -----
    There are several definitions of "circular median" in the literature. Seminal
    references for circular statistics [1]_ and [2]_ mention the bisecting property
    of the circular median and its connection with minimization of the circular mean
    deviation. However, [1]_ prioritizes minimization of the circular mean deviation
    as the defining property, whereas [2]_ requires that "the majority of the data
    points are nearer to [the median] than [its antipode]".

    Although the two definitions produce the same circular median in many cases
    (unimodal distribution, no ties), [3]_ demonstrates that the two definitions can
    lead to different - and in fact, entirely opposite - results. [3]_ also compares the
    arc-distance and bisecting medians with two other definitions adopted from
    vector-space literature, and it concludes that the arc-distance median is the only
    one that provides all four properties under consideration. This, and the fact that
    the arc-distance definition is used by default in other circular statistics software
    (e.g. [4]_, [5]_, [6]_), contributed to the choice of ``'arc-distance'`` as the
    default convention.

    [7]_ addresses the fact that many points - even continuous arcs of the unit circle -
    may satisfy the definition of a circular median. It proposes that "the estimate of
    the population circular median be the average (circular mean) of all angles
    satisfying the definition of median", and that "for odd samples, the candidate
    values are the observations themselves". However, it suggests that "for even
    samples, the candidate values are the midpoints of all neighboring observation".
    `circmedian` always finds all points *among the observations* (and their antipodes,
    in the ``'bisecting'`` case) that satisfy the chosen definition of a circular median
    and returns their circular mean. If the circular mean is poorly defined (i.e. the
    circular variance is with a small tolerance of ``1.0``), `circmedian` returns NaN.

    References
    ----------
    .. [1] Fisher, Nicholas I. *Statistical Analysis of Circular Data*. Cambridge
           University Press, 1995.
    .. [2] Mardia, K. V. and Jupp, P. E. *Directional Statistics*.
           John Wiley & Sons, 1999.
    .. [3] Storath, Martin, and Andreas Weinmann. "Fast median filtering for phase or
           orientation data." *IEEE Transactions on Pattern Analysis and Machine
           Intelligence* 40.3 (2017): 639-652.
    .. [4] Berens, Philipp. "CircStat: a MATLAB toolbox for circular statistics."
           *Journal of Statistical Software* 31 (2009): 1-21.
    .. [5] Lund, Ulric, Claudio Agostinelli, and Maintainer Claudio Agostinelli.
          "Package 'circular'." Repository CRAN 775.5 (2017): 20-135.
    .. [6] Huang, Ziwei. "PyCircStat2: Circular statistics with Python".
           https://github.com/circstat/pycircstat2.
    .. [7] Otieno, B., and Christine M. Anderson-Cook. "A more efficient way of
           obtaining a unique median estimate for circular data." *Journal of Modern
           Applied Statistical Methods* 2.1 (2003): 15.

    Examples
    --------
    Consider Example 1 from [3]_, expressed in degrees for readability.

    >>> import numpy as np
    >>> from scipy import stats
    >>> sample = np.array([101.25, 101.25, 0, -101.25, -101.25])

    The unique arc-distance median is 0 because it minimizes the mean arc distance to
    the other observations.

    >>> displacements = (sample[:, np.newaxis] - sample[np.newaxis, :]) % 360
    >>> distances = np.minimum(displacements, 360 - displacements)
    >>> mean_deviations = np.mean(distances, axis=-1)
    >>> sample[mean_deviations == np.min(mean_deviations)]
    array([0.])

    As expected,

    >>> stats.circmedian(sample, low=-180, high=180)
    np.float64(0.0)

    Note that the chord from *any* observation to its antipode bisects the observations
    in this example. Furthermore, each antipode is within 90 degrees of three or four
    observations, whereas each observation has either two or zero other observations
    within 90 degrees. Therefore, the antipodes ``[-78.75, -78.75, 180, 78.75, 78.75]``
    all satisfy the definition of a bisecting median. The circular mean of these is 180
    degrees,

    >>> stats.circmean([-78.75, -78.75, 180, 78.75, 78.75], low=-180, high=180)
    np.float64(179.99999999999994)

    and so:

    >>> stats.circmedian(sample, low=-180, high=180, convention='bisecting')
    np.float64(179.99999999999994)

    The discrepancy between the two conventions is extreme in this constructed example,
    but for data from a unimodal, continuous distribution, the two definitions tend to
    produce similar or identical results.

    >>> rng = np.random.default_rng(8180124976)
    >>> sample = rng.vonmises(mu=0, kappa=1, size=10)
    >>> stats.circmedian(sample), stats.circmedian(sample, convention='bisecting')
    (np.float64(0.44401177698735417), np.float64(0.44401177698735417))

    """
    xp = array_namespace(sample)
    sample = xp_promote(sample, force_floating=True, xp=xp)
    tol = 10*xp.finfo(sample).eps
    T = high - low
    two_pi = 2 * math.pi
    n = sample.shape[-1]

    # axis=-1 is guaranteed by the _axis_nan_policy decorator
    if convention == 'arc-distance':
        mad = _circmeandev(sample, sample, high=high, low=low)
        min_mad = xp.min(mad, axis=-1, keepdims=True)
        i = (mad - min_mad) < tol*min_mad
    elif convention == 'bisecting':
        displacements = (sample[..., xp.newaxis, :] - sample[..., :, xp.newaxis]) % T
        displacements = xp.where(displacements > T/2, displacements-T, displacements)
        distances = xp.abs(displacements)

        count_close = xp.count_nonzero(distances < T/4 , axis=-1) - 1
        count_far = xp.count_nonzero(distances > T/4, axis=-1)
        sample = xp.where(count_far > count_close, (sample + T/2) % T, sample)

        i_antipodal = (distances == T/2)
        i_left = (displacements <= 0) | i_antipodal
        i_right = (displacements >= 0) | i_antipodal

        count_left = xp.astype(xp.count_nonzero(i_left, axis=-1), sample.dtype)
        count_right = xp.astype(xp.count_nonzero(i_right, axis=-1), sample.dtype)
        i = (count_left >= n/2) & (count_right >= n/2)
    else:
        raise ValueError("`convention` must be either 'arc-distance' or 'bisecting'.")

    phi = xp.exp(1j * (sample - low) * two_pi / T)
    median = _xp_mean(phi, weights=xp.astype(i, phi.dtype), axis=-1)
    median = xpx.at(median)[xp.abs(median) < tol].set(xp.nan)
    median = xp.atan2(xp.imag(median), xp.real(median)) % (2*math.pi)
    median = median * T / two_pi + low
    return median[()]
