import numpy as np


def _validate_lower_upper(lower, upper):
    if lower.ndim != 1 or upper.ndim != 1:
        raise ValueError('The input arrays must be one-dimensional.')
    if np.shape(lower) != np.shape(upper):
        raise ValueError('lower and upper must have the same length.')
    if np.isnan(lower).any() or np.isnan(upper).any():
        raise ValueError('lower and upper must not contain nan.')
    if (np.isinf(lower) & np.isinf(upper)).any():
        raise ValueError('lower and upper must not both be infinite.')
    if (lower > upper).any():
        raise ValueError('Elements of lower must not be greater than the '
                         'corresponding elements of upper.')


def _validate_x_censored(x, censored):
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('x must be one-dimensional.')
    censored = np.asarray(censored)
    if censored.ndim != 1:
        raise ValueError('censored must be one-dimensional')
    if (~np.isfinite(x)).any():
        raise ValueError('x must not contain nan or inf.')
    if censored.size != x.size:
        raise ValueError('x and censored must have the same length.')
    return x, censored.astype(bool)


class CensoredData:
    """
    Instances of this class represent censored data.

    Instances may be passed to the ``fit`` method of continuous
    univariate SciPy distributions for maximum likelihood estimation.

    Left-, right-, and interval-censored data can be represented.

    For interval-censored data, pass the lower and upper limits of
    the intervals as one-dimensional arrays to the ``CensoredData``
    constuctor.  By using the values `np.inf` and `-np.inf`, left-
    and right-censored data may also be created with the constructor.

    For convenience, the class methods ``left_censored`` and
    ``right_censored`` are provided to create a ``CensoredData``
    instance from a single one-dimensional array of measurements
    and a corresponding boolean array to indicate which measurements
    are censored.

    Parameters
    ----------
    lower, upper : array_like
        The data, a possibly mixed collection of censored and uncensored
        measurements.  The arrays must be one-dimensional and have the same
        length.  The values are interpreted as follows:

        * ``lower[k] == upper[k]``: the common value is an uncensored
          measurement.
        * ``upper[k] == inf``: the measurement is right-censored; the
          value in ``lower[k]`` is the lower bound for the true unknown
          value.
        * ``lower[k] == -inf``: the measurement is left-censored; the value
          in ``upper[k]`` is the upper bound for the true unknown value.
        * ``lower[k] < upper[k]``, and both values are finite: the
          measurement is interval-censored; the true unknown value is
          between ``lower[k]`` and ``upper[k]``.

    Examples
    --------
    In the most general case, a censored data set may contain values that
    are left-censored, right-censored, interval-censored and not censored.
    For example, here we create a data set with length five:

    >>> from scipy.stats import CensoredData
    >>> data = CensoredData([-np.inf, 1, 1.5, 2, 10],
    ...                     [0,       1, 1.5, 3, np.inf])
    >>> print(data)
    CensoredData(5 values: 2 not censored, 1 left-censored,
    1 right-censored, 1 interval-censored)

    The first value represents a left-censored measurement; the true (but
    unknown) value is less than or equal to 0.  The second and third
    values, 1 and 1.5, are not censored.  The fourth value represents an
    interval-censored measurement; the true value is in the interval [2, 3].
    The last value represents a right-censored measurement; the true value
    is greater than or equal to 10.

    A common case is to have a mix of uncensored measurements and censored
    measurements that are all right-censored (or all left-censored). For
    example, consider an experiment in which six devices are started at
    various times and left running until they fail.  Assume that time is
    measured in hours, and the experiment is stopped after 30 hours, even
    if all the devices have not failed by that time.  We might end up with
    data such as this::

        Device  Start-time  Fail-time  Time-to-failure
           1         0         13           13
           2         2         24           22
           3         5         22           17
           4         8         23           15
           5        10        ***          >20
           6        12        ***          >18

    Two of the devices had not failed when the experiment was stopped;
    the measurements of the time-to-failure for these two devices are
    right-censored.  We'll use the method `CensoredData.right_censored` to
    create a representation of this data.  The time-to-failure measurements
    are put the list ``ttf``.  The ``censored`` list indicates which values
    in ``ttf`` are censored.

    >>> ttf = [13, 22, 17, 15, 20, 18]
    >>> censored = [False, False, False, False, True, True]

    Pass these lists to `CensoredData.right_censored` to create an
    instance of `CensoredData`.

    >>> data = CensoredData.right_censored(ttf, censored)
    >>> print(data)
    CensoredData(6 values: 4 not censored, 2 right-censored)
    """

    def __init__(self, lower, upper):
        lower = np.asarray(lower)
        upper = np.asarray(upper)
        _validate_lower_upper(lower, upper)
        self._lower = lower
        self._upper = upper

        # Precompute masks for the various types of censoring.
        self._not_censored = lower == upper
        self._left_censored = np.isneginf(lower)
        self._right_censored = np.isposinf(upper)
        self._interval_censored = ~(self._not_censored | self._left_censored
                                    | self._right_censored)

    def __repr__(self):
        lower_str = " ".join(np.array_repr(self._lower).split())
        upper_str = " ".join(np.array_repr(self._upper).split())
        return (f"CensoredData({lower_str}, {upper_str})")

    def __str__(self):
        n = len(self._lower)
        num_nc = np.count_nonzero(self._not_censored)
        num_lc = np.count_nonzero(self._left_censored)
        num_rc = np.count_nonzero(self._right_censored)
        num_ic = np.count_nonzero(self._interval_censored)
        parts = [f'{num_nc} not censored']
        if num_lc > 0:
            parts.append(f'{num_lc} left-censored')
        if num_rc > 0:
            parts.append(f'{num_rc} right-censored')
        if num_ic > 0:
            parts.append(f'{num_ic} interval-censored')
        return f'CensoredData({n} values: ' + ', '.join(parts) + ')'

    def __len__(self):
        """
        The number of values (censored and not censored).
        """
        return len(self._lower)

    def num_censored(self):
        return len(self._lower) - np.count_nonzero(self._not_censored)

    @classmethod
    def right_censored(cls, x, censored):
        """
        Create a CensoredData instance of right-censored data.

        Parameters
        ----------
        x : array_like
            `x` is the array of observed data or measurements.
            `x` must be a one-dimensional sequence of numbers.
        censored : array_like of bool
            `censored` must be a one-dimensional sequence of boolean
            values.  If ``censored[k]`` is True, the corresponding value
            in `x` is right-censored.  That is, the value ``x[k]``
            is the lower bound of the true (but unknown) value.
        """
        x, censored = _validate_x_censored(x, censored)
        lower = 1.0*x  # Copy x while ensuring lower is floating point.
        upper = lower.copy()
        upper[censored] = np.inf
        return cls(lower=lower, upper=upper)

    @classmethod
    def left_censored(cls, x, censored):
        """
        Create a CensoredData instance of left-censored data.

        Parameters
        ----------
        x : array_like
            `x` is the array of observed data or measurements.
            `x` must be a one-dimensional sequence of numbers.
        censored : array_like of bool
            `censored` must be a one-dimensional sequence of boolean
            values.  If ``censored[k]`` is True, the corresponding value
            in `x` is left-censored.  That is, the value ``x[k]``
            is the upper bound of the true (but unknown) value.
        """
        x, censored = _validate_x_censored(x, censored)
        upper = 1.0*x  # Copy x while ensuring upper is floating point.
        lower = upper.copy()
        lower[censored] = -np.inf
        return cls(lower=lower, upper=upper)

    def _uncensor(self):
        """
        This function is used when a non-censored version of the data
        is needed to create a rough estimate of the parameters of a
        distribution via the method of moments or some similar method.
        The data is "uncensored" by taking the given endpoints as the
        data for the left- or right-censored data, and the mean for the
        interval-censored data.
        """
        data = self._lower.copy()
        data[self._left_censored] = self._upper[self._left_censored]
        ic = self._interval_censored
        data[ic] = 0.5*(self._lower[ic] + self._upper[ic])
        return data
