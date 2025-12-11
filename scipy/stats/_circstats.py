import math
from scipy._lib._array_api import (
    array_namespace, xp_promote, xp_capabilities, xp_size, xp_vector_norm, is_marray,
)
from scipy.stats._axis_nan_policy import _axis_nan_policy_factory


__all__ = ['circmean', 'circvar', 'circstd', 'directional_stats']


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
