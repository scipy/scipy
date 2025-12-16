import math
from scipy._lib import array_api_extra as xpx
from scipy._lib._array_api import (
    array_namespace, xp_promote, xp_capabilities, xp_size, xp_vector_norm, is_marray,
)
from scipy.stats._axis_nan_policy import _axis_nan_policy_factory, _broadcast_arrays
from scipy.stats._stats_py import _xp_mean


__all__ = ['circmean', 'circvar', 'circstd', 'circmedian', 'directional_stats']


def _circfuncs_common(samples, period, xp=None):
    xp = array_namespace(samples) if xp is None else xp

    samples = xp_promote(samples, force_floating=True, xp=xp)

    # Recast samples as radians that range between 0 and 2 pi and calculate
    # the sine and cosine
    scaled_samples = samples * ((2.0 * math.pi) / period)
    sin_samp = xp.sin(scaled_samples)
    cos_samp = xp.cos(scaled_samples)

    return samples, sin_samp, cos_samp


@xp_capabilities()
@_axis_nan_policy_factory(
    lambda x: x, n_outputs=1, default_axis=None,
    result_to_tuple=lambda x, _: (x,)
)
def circmean(samples, high=2*math.pi, low=0, axis=None, nan_policy='propagate'):
    r"""Compute the circular mean of a sample of angle observations.

    Given :math:`n` angle observations :math:`x_1, \cdots, x_n` measured in
    radians, their *circular mean* is defined by ([1]_, Eq. 2.2.4)

    .. math::

       \mathrm{Arg} \left( \frac{1}{n} \sum_{k=1}^n e^{i x_k} \right)

    where :math:`i` is the imaginary unit and :math:`\mathop{\mathrm{Arg}} z`
    gives the principal value of the argument of complex number :math:`z`,
    restricted to the range :math:`[0,2\pi]` by default.  :math:`z` in the
    above expression is known as the `mean resultant vector`.

    Parameters
    ----------
    samples : array_like
        Input array of angle observations.  The value of a full angle is
        equal to ``(high - low)``.
    high : float, optional
        Upper boundary of the principal value of an angle.  Default is ``2*pi``.
    low : float, optional
        Lower boundary of the principal value of an angle.  Default is ``0``.

    Returns
    -------
    circmean : float
        Circular mean, restricted to the range ``[low, high]``.

        If the mean resultant vector is zero, an input-dependent,
        implementation-defined number between ``[low, high]`` is returned.
        If the input array is empty, ``np.nan`` is returned.

    See Also
    --------
    circstd : Circular standard deviation.
    circvar : Circular variance.

    References
    ----------
    .. [1] Mardia, K. V. and Jupp, P. E. *Directional Statistics*.
           John Wiley & Sons, 1999.

    Examples
    --------
    For readability, all angles are printed out in degrees.

    >>> import numpy as np
    >>> from scipy.stats import circmean
    >>> import matplotlib.pyplot as plt
    >>> angles = np.deg2rad(np.array([20, 30, 330]))
    >>> circmean = circmean(angles)
    >>> np.rad2deg(circmean)
    7.294976657784009

    >>> mean = angles.mean()
    >>> np.rad2deg(mean)
    126.66666666666666

    Plot and compare the circular mean against the arithmetic mean.

    >>> plt.plot(np.cos(np.linspace(0, 2*np.pi, 500)),
    ...          np.sin(np.linspace(0, 2*np.pi, 500)),
    ...          c='k')
    >>> plt.scatter(np.cos(angles), np.sin(angles), c='k')
    >>> plt.scatter(np.cos(circmean), np.sin(circmean), c='b',
    ...             label='circmean')
    >>> plt.scatter(np.cos(mean), np.sin(mean), c='r', label='mean')
    >>> plt.legend()
    >>> plt.axis('equal')
    >>> plt.show()

    """
    xp = array_namespace(samples)
    # Needed for non-NumPy arrays to get appropriate NaN result
    # Apparently atan2(0, 0) is 0, even though it is mathematically undefined
    if xp_size(samples) == 0:
        return xp.mean(samples, axis=axis)
    period = high - low
    samples, sin_samp, cos_samp = _circfuncs_common(samples, period, xp=xp)
    sin_sum = xp.sum(sin_samp, axis=axis)
    cos_sum = xp.sum(cos_samp, axis=axis)
    res = xp.atan2(sin_sum, cos_sum)

    res = res[()] if res.ndim == 0 else res
    return (res * (period / (2.0 * math.pi)) - low) % period + low


@xp_capabilities()
@_axis_nan_policy_factory(
    lambda x: x, n_outputs=1, default_axis=None,
    result_to_tuple=lambda x, _: (x,)
)
def circvar(samples, high=2*math.pi, low=0, axis=None, nan_policy='propagate'):
    r"""Compute the circular variance of a sample of angle observations.

    Given :math:`n` angle observations :math:`x_1, \cdots, x_n` measured in
    radians, their *circular variance* is defined by ([2]_, Eq. 2.3.3)

    .. math::

       1 - \left| \frac{1}{n} \sum_{k=1}^n e^{i x_k} \right|

    where :math:`i` is the imaginary unit and :math:`|z|` gives the length
    of the complex number :math:`z`.  :math:`|z|` in the above expression
    is known as the `mean resultant length`.

    Parameters
    ----------
    samples : array_like
        Input array of angle observations.  The value of a full angle is
        equal to ``(high - low)``.
    high : float, optional
        Upper boundary of the principal value of an angle.  Default is ``2*pi``.
    low : float, optional
        Lower boundary of the principal value of an angle.  Default is ``0``.

    Returns
    -------
    circvar : float
        Circular variance.  The returned value is in the range ``[0, 1]``,
        where ``0`` indicates no variance and ``1`` indicates large variance.

        If the input array is empty, ``np.nan`` is returned.

    See Also
    --------
    circmean : Circular mean.
    circstd : Circular standard deviation.

    Notes
    -----
    In the limit of small angles, the circular variance is close to
    half the 'linear' variance if measured in radians.

    References
    ----------
    .. [1] Fisher, N.I. *Statistical analysis of circular data*. Cambridge
           University Press, 1993.
    .. [2] Mardia, K. V. and Jupp, P. E. *Directional Statistics*.
           John Wiley & Sons, 1999.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import circvar
    >>> import matplotlib.pyplot as plt
    >>> samples_1 = np.array([0.072, -0.158, 0.077, 0.108, 0.286,
    ...                       0.133, -0.473, -0.001, -0.348, 0.131])
    >>> samples_2 = np.array([0.111, -0.879, 0.078, 0.733, 0.421,
    ...                       0.104, -0.136, -0.867,  0.012,  0.105])
    >>> circvar_1 = circvar(samples_1)
    >>> circvar_2 = circvar(samples_2)

    Plot the samples.

    >>> fig, (left, right) = plt.subplots(ncols=2)
    >>> for image in (left, right):
    ...     image.plot(np.cos(np.linspace(0, 2*np.pi, 500)),
    ...                np.sin(np.linspace(0, 2*np.pi, 500)),
    ...                c='k')
    ...     image.axis('equal')
    ...     image.axis('off')
    >>> left.scatter(np.cos(samples_1), np.sin(samples_1), c='k', s=15)
    >>> left.set_title(f"circular variance: {np.round(circvar_1, 2)!r}")
    >>> right.scatter(np.cos(samples_2), np.sin(samples_2), c='k', s=15)
    >>> right.set_title(f"circular variance: {np.round(circvar_2, 2)!r}")
    >>> plt.show()

    """
    xp = array_namespace(samples)
    period = high - low
    samples, sin_samp, cos_samp = _circfuncs_common(samples, period, xp=xp)
    sin_mean = xp.mean(sin_samp, axis=axis)
    cos_mean = xp.mean(cos_samp, axis=axis)
    hypotenuse = (sin_mean**2. + cos_mean**2.)**0.5
    # hypotenuse can go slightly above 1 due to rounding errors
    R = xp.clip(hypotenuse, max=1.)

    res = 1. - R
    return res


@xp_capabilities()
@_axis_nan_policy_factory(
    lambda x: x, n_outputs=1, default_axis=None,
    result_to_tuple=lambda x, _: (x,)
)
def circstd(samples, high=2*math.pi, low=0, axis=None, nan_policy='propagate', *,
            normalize=False):
    r"""
    Compute the circular standard deviation of a sample of angle observations.

    Given :math:`n` angle observations :math:`x_1, \cdots, x_n` measured in
    radians, their `circular standard deviation` is defined by
    ([2]_, Eq. 2.3.11)

    .. math::

       \sqrt{ -2 \log \left| \frac{1}{n} \sum_{k=1}^n e^{i x_k} \right| }

    where :math:`i` is the imaginary unit and :math:`|z|` gives the length
    of the complex number :math:`z`.  :math:`|z|` in the above expression
    is known as the `mean resultant length`.

    Parameters
    ----------
    samples : array_like
        Input array of angle observations.  The value of a full angle is
        equal to ``(high - low)``.
    high : float, optional
        Upper boundary of the principal value of an angle.  Default is ``2*pi``.
    low : float, optional
        Lower boundary of the principal value of an angle.  Default is ``0``.
    normalize : boolean, optional
        If ``False`` (the default), the return value is computed from the
        above formula with the input scaled by ``(2*pi)/(high-low)`` and
        the output scaled (back) by ``(high-low)/(2*pi)``.  If ``True``,
        the output is not scaled and is returned directly.

    Returns
    -------
    circstd : float
        Circular standard deviation, optionally normalized.

        If the input array is empty, ``np.nan`` is returned.

    See Also
    --------
    circmean : Circular mean.
    circvar : Circular variance.

    Notes
    -----
    In the limit of small angles, the circular standard deviation is close
    to the 'linear' standard deviation if ``normalize`` is ``False``.

    References
    ----------
    .. [1] Mardia, K. V. (1972). 2. In *Statistics of Directional Data*
       (pp. 18-24). Academic Press. :doi:`10.1016/C2013-0-07425-7`.
    .. [2] Mardia, K. V. and Jupp, P. E. *Directional Statistics*.
           John Wiley & Sons, 1999.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import circstd
    >>> import matplotlib.pyplot as plt
    >>> samples_1 = np.array([0.072, -0.158, 0.077, 0.108, 0.286,
    ...                       0.133, -0.473, -0.001, -0.348, 0.131])
    >>> samples_2 = np.array([0.111, -0.879, 0.078, 0.733, 0.421,
    ...                       0.104, -0.136, -0.867,  0.012,  0.105])
    >>> circstd_1 = circstd(samples_1)
    >>> circstd_2 = circstd(samples_2)

    Plot the samples.

    >>> fig, (left, right) = plt.subplots(ncols=2)
    >>> for image in (left, right):
    ...     image.plot(np.cos(np.linspace(0, 2*np.pi, 500)),
    ...                np.sin(np.linspace(0, 2*np.pi, 500)),
    ...                c='k')
    ...     image.axis('equal')
    ...     image.axis('off')
    >>> left.scatter(np.cos(samples_1), np.sin(samples_1), c='k', s=15)
    >>> left.set_title(f"circular std: {np.round(circstd_1, 2)!r}")
    >>> right.plot(np.cos(np.linspace(0, 2*np.pi, 500)),
    ...            np.sin(np.linspace(0, 2*np.pi, 500)),
    ...            c='k')
    >>> right.scatter(np.cos(samples_2), np.sin(samples_2), c='k', s=15)
    >>> right.set_title(f"circular std: {np.round(circstd_2, 2)!r}")
    >>> plt.show()

    """
    xp = array_namespace(samples)
    period = high - low
    samples, sin_samp, cos_samp = _circfuncs_common(samples, period, xp=xp)
    sin_mean = xp.mean(sin_samp, axis=axis)  # [1] (2.2.3)
    cos_mean = xp.mean(cos_samp, axis=axis)  # [1] (2.2.3)
    hypotenuse = (sin_mean**2. + cos_mean**2.)**0.5
    # hypotenuse can go slightly above 1 due to rounding errors
    R = xp.clip(hypotenuse, max=1.)  # [1] (2.2.4)

    res = (-2*xp.log(R))**0.5+0.0  # torch.pow returns -0.0 if R==1
    if not normalize:
        res *= (high-low)/(2.*math.pi)  # [1] (2.3.14) w/ (2.3.7)
    return res


@xp_capabilities()
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


class DirectionalStats:
    def __init__(self, mean_direction, mean_resultant_length):
        self.mean_direction = mean_direction
        self.mean_resultant_length = mean_resultant_length

    def __repr__(self):
        return (f"DirectionalStats(mean_direction={self.mean_direction},"
                f" mean_resultant_length={self.mean_resultant_length})")


@xp_capabilities()
def directional_stats(samples, *, axis=0, normalize=True):
    """
    Computes sample statistics for directional data.

    Computes the directional mean (also called the mean direction vector) and
    mean resultant length of a sample of vectors.

    The directional mean is a measure of "preferred direction" of vector data.
    It is analogous to the sample mean, but it is for use when the length of
    the data is irrelevant (e.g. unit vectors).

    The mean resultant length is a value between 0 and 1 used to quantify the
    dispersion of directional data: the smaller the mean resultant length, the
    greater the dispersion. Several definitions of directional variance
    involving the mean resultant length are given in [1]_ and [2]_.

    Parameters
    ----------
    samples : array_like
        Input array. Must be at least two-dimensional, and the last axis of the
        input must correspond with the dimensionality of the vector space.
        When the input is exactly two dimensional, this means that each row
        of the data is a vector observation.
    axis : int, default: 0
        Axis along which the directional mean is computed.
    normalize: boolean, default: True
        If True, normalize the input to ensure that each observation is a
        unit vector. It the observations are already unit vectors, consider
        setting this to False to avoid unnecessary computation.

    Returns
    -------
    res : DirectionalStats
        An object containing attributes:

        mean_direction : ndarray
            Directional mean.
        mean_resultant_length : ndarray
            The mean resultant length [1]_.

    See Also
    --------
    circmean: circular mean; i.e. directional mean for 2D *angles*
    circvar: circular variance; i.e. directional variance for 2D *angles*

    Notes
    -----
    This uses a definition of directional mean from [1]_.
    Assuming the observations are unit vectors, the calculation is as follows.

    .. code-block:: python

        mean = samples.mean(axis=0)
        mean_resultant_length = np.linalg.norm(mean)
        mean_direction = mean / mean_resultant_length

    This definition is appropriate for *directional* data (i.e. vector data
    for which the magnitude of each observation is irrelevant) but not
    for *axial* data (i.e. vector data for which the magnitude and *sign* of
    each observation is irrelevant).

    Several definitions of directional variance involving the mean resultant
    length ``R`` have been proposed, including ``1 - R`` [1]_, ``1 - R**2``
    [2]_, and ``2 * (1 - R)`` [2]_. Rather than choosing one, this function
    returns ``R`` as attribute `mean_resultant_length` so the user can compute
    their preferred measure of dispersion.

    References
    ----------
    .. [1] Mardia, Jupp. (2000). *Directional Statistics*
       (p. 163). Wiley.

    .. [2] https://en.wikipedia.org/wiki/Directional_statistics

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.stats import directional_stats
    >>> data = np.array([[3, 4],    # first observation, 2D vector space
    ...                  [6, -8]])  # second observation
    >>> dirstats = directional_stats(data)
    >>> dirstats.mean_direction
    array([1., 0.])

    In contrast, the regular sample mean of the vectors would be influenced
    by the magnitude of each observation. Furthermore, the result would not be
    a unit vector.

    >>> data.mean(axis=0)
    array([4.5, -2.])

    An exemplary use case for `directional_stats` is to find a *meaningful*
    center for a set of observations on a sphere, e.g. geographical locations.

    >>> data = np.array([[0.8660254, 0.5, 0.],
    ...                  [0.8660254, -0.5, 0.]])
    >>> dirstats = directional_stats(data)
    >>> dirstats.mean_direction
    array([1., 0., 0.])

    The regular sample mean on the other hand yields a result which does not
    lie on the surface of the sphere.

    >>> data.mean(axis=0)
    array([0.8660254, 0., 0.])

    The function also returns the mean resultant length, which
    can be used to calculate a directional variance. For example, using the
    definition ``Var(z) = 1 - R`` from [2]_ where ``R`` is the
    mean resultant length, we can calculate the directional variance of the
    vectors in the above example as:

    >>> 1 - dirstats.mean_resultant_length
    0.13397459716167093
    """
    xp = array_namespace(samples)
    samples = xp.asarray(samples)

    if samples.ndim < 2:
        raise ValueError("samples must at least be two-dimensional. "
                         f"Instead samples has shape: {tuple(samples.shape)}")
    samples = xp.moveaxis(samples, axis, 0)

    if is_marray(xp):
        _xp = array_namespace(samples.mask)
        mask = _xp.any(samples.mask, axis=-1, keepdims=True)
        samples = xp.asarray(samples.data, mask=mask)

    if normalize:
        vectornorms = xp_vector_norm(samples, axis=-1, keepdims=True, xp=xp)
        samples = samples/vectornorms
    mean = xp.mean(samples, axis=0)
    mean_resultant_length = xp_vector_norm(mean, axis=-1, keepdims=True, xp=xp)
    mean_direction = mean / mean_resultant_length
    mrl = xp.squeeze(mean_resultant_length, axis=-1)
    mean_resultant_length = mrl[()] if mrl.ndim == 0 else mrl
    return DirectionalStats(mean_direction, mean_resultant_length)
