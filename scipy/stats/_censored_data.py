import numpy as np


def _validate_1d(a, name):
    if np.ndim(a) != 1:
        raise ValueError(f'`{name}` must be a one-dimensional sequence.')
    if np.isnan(a).any():
        raise ValueError(f'`{name}` must not contain nan.')
    if np.isinf(a).any():
        raise ValueError(f'`{name}` must contain only finite values.')


def _validate_intervals(intervals):
    intervals = np.asarray(intervals)
    if intervals.shape == (0,):
        # The input was a sequence with length 0.
        intervals = intervals.reshape((0, 2))
    if intervals.ndim != 2 or intervals.shape[-1] != 2:
        raise ValueError('`intervals` must be a two-dimensional array'
                         ' with shape (m, 2).')

    if np.isnan(intervals).any():
        raise ValueError('`intervals` must not contain nan.')
    if np.isinf(intervals).all(axis=1).any():
        raise ValueError('In each row in `intervals`, both values must not'
                         ' be infinite.')
    if (intervals[:, 0] > intervals[:, 1]).any():
        raise ValueError('In each row of `intervals`, the left value must not'
                         ' exceed the right value.')

    uncensored = intervals[:, 0] == intervals[:, 1]
    left_censored = np.isinf(intervals[:, 0])
    right_censored = np.isinf(intervals[:, 1])
    interval_censored = np.isfinite(intervals).all(axis=1) & ~uncensored

    x2 = intervals[uncensored, 0]
    left2 = intervals[left_censored, 1]
    right2 = intervals[right_censored, 0]
    intervals2 = intervals[interval_censored]

    return x2, left2, right2, intervals2


def _validate_x_censored(x, censored):
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('`x` must be one-dimensional.')
    censored = np.asarray(censored)
    if censored.ndim != 1:
        raise ValueError('`censored` must be one-dimensional.')
    if (~np.isfinite(x)).any():
        raise ValueError('`x` must not contain nan or inf.')
    if censored.size != x.size:
        raise ValueError('`x` and `censored` must have the same length.')
    return x, censored.astype(bool)


class CensoredData:
    """
    Instances of this class represent censored data.

    Instances may be passed to the ``fit`` method of continuous
    univariate SciPy distributions for maximum likelihood estimation.
    The *only* method of the univariate continuous distributions that
    understands `CensoredData` is the ``fit`` method.  An instance of
    `CensoredData` can not be passed to methods such as ``pdf`` and
    ``cdf``.

    An observation is said to be *censored* when the precise value is unknown,
    but it has a known upper and/or lower bound.  The conventional terminology
    is:

    * left-censored: an observation is below a certain value but it is
      unknown by how much.
    * right-censored: an observation is above a certain value but it is
      unknown by how much.
    * interval-censored: an observation lies somewhere on an interval between
      two values.

    Left-, right-, and interval-censored data can be represented by
    `CensoredData`.

    For interval-censored data, pass the lower and upper limits of
    the intervals as one-dimensional arrays to the `CensoredData`
    constuctor.  By using the values `np.inf` and `-np.inf`, left-
    and right-censored data may also be created with the constructor.

    For convenience, the class methods ``left_censored`` and
    ``right_censored`` are provided to create a `CensoredData`
    instance from a single one-dimensional array of measurements
    and a corresponding boolean array to indicate which measurements
    are censored.

    Parameters
    ----------
    x : array_like, 1D
        Uncensored observations.
    left : array_like, 1D
        Left-censored observations.
    right : array_like, 1D
        Right-censored observations.
    intervals : array_like, 2D, with shape (m, 2)
        Interval-censored observations.  Each row ``intervals[k, :]``
        represents the interval for the kth interval-censored observation.

    Notes
    -----
    In the input array `intervals`, the lower bound of the interval may
    be ``-inf``, and the upper bound may be ``inf``, but at least one must be
    finite). When the lower bound is ``-inf``, the row represents a left
    censored observation, and when the upper bound is ``inf``, the row
    represents a right-censored observation.  If the length of an interval
    is 0 (i.e. ``intervals[k, 0] == intervals[k, 1]``, the observation is
    treated as uncensored.  Consequently, one can represent all the types
    of censored and uncensored data in ``intervals``, but it is generally
    more convenient to use `x`, `left` and `right` for uncensored,
    left-censored and right-censored observations, respectively.

    Examples
    --------
    In the most general case, a censored data set may contain values that
    are left-censored, right-censored, interval-censored, and uncensored.
    For example, here we create a data set with five observations.  Two
    are uncensored (values 1 and 1.5), one is a left-censored observation
    of 0, one is a right-censored observation of 10 and one is
    interval-censored in the interval [2, 3].

    >>> from scipy.stats import CensoredData
    >>> data = CensoredData(x=[1, 1.5], left=[0], right=[10],
    ...                     intervals=[[2, 3]])
    >>> print(data)
    CensoredData(5 values: 2 not censored, 1 left-censored,
    1 right-censored, 1 interval-censored)

    Equivalently,
    >>> data = CensoredData(intervals=[[1, 1],
    ...                                [1.5, 1.5],
    ...                                [-np.inf, 0],
    ...                                [10, np.inf],
    ...                                [2, 3]])
    >>> print(data)
    CensoredData(5 values: 2 not censored, 1 left-censored,
    1 right-censored, 1 interval-censored)

    A common case is to have a mix of uncensored observations and censored
    observations that are all right-censored (or all left-censored). For
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
    the observations of the time-to-failure for these two devices are
    right-censored.  We can represent this data with

    >>> data = CensoredData(x=[13, 22, 17, 15], right=[20, 18])
    >>> print(data)
    CensoredData(6 values: 4 not censored, 2 right-censored)

    Alternatively, we can use the method `CensoredData.right_censored` to
    create a representation of this data.  The time-to-failure observations
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

    def __init__(self, x=None, *, left=None, right=None, intervals=None):
        if x is None:
            x = []
        if left is None:
            left = []
        if right is None:
            right = []
        if intervals is None:
            intervals = np.empty((0, 2))

        _validate_1d(x, 'x')
        _validate_1d(left, 'left')
        _validate_1d(right, 'right')
        x2, left2, right2, intervals2 = _validate_intervals(intervals)

        self._x = np.concatenate((x, x2))
        self._left = np.concatenate((left, left2))
        self._right = np.concatenate((right, right2))
        # Note that by construction, the private attribute _intervals
        # will be a 2D array that contains only finite values representing
        # intervals with nonzero but finite length.
        self._intervals = intervals2

    def __repr__(self):
        x_str = " ".join(np.array_repr(self._x).split())
        left_str = " ".join(np.array_repr(self._left).split())
        right_str = " ".join(np.array_repr(self._right).split())
        intervals_str = " ".join(np.array_repr(self._intervals).split())
        return (f"CensoredData(x={x_str}, left={left_str}, right={right_str},"
                f" intervals={intervals_str})")

    def __str__(self):
        num_nc = len(self._x)
        num_lc = len(self._left)
        num_rc = len(self._right)
        num_ic = len(self._intervals)
        n = num_nc + num_lc + num_rc + num_ic
        parts = [f'{num_nc} not censored']
        if num_lc > 0:
            parts.append(f'{num_lc} left-censored')
        if num_rc > 0:
            parts.append(f'{num_rc} right-censored')
        if num_ic > 0:
            parts.append(f'{num_ic} interval-censored')
        return f'CensoredData({n} values: ' + ', '.join(parts) + ')'

    # This is not a complete implementation of the arithmetic operators.
    # All we need is subtracting a scalar and dividing by a scalar.

    def __sub__(self, other):
        return CensoredData(x=self._x - other,
                            left=self._left - other,
                            right=self._right - other,
                            intervals=self._intervals - other)

    def __truediv__(self, other):
        return CensoredData(x=self._x / other,
                            left=self._left / other,
                            right=self._right / other,
                            intervals=self._intervals / other)

    def __len__(self):
        """
        The number of values (censored and not censored).
        """
        return (len(self._x) + len(self._left) + len(self._right)
                + len(self._intervals))

    def num_censored(self):
        """
        Number of censored values.
        """
        return len(self._left) + len(self._right) + len(self._intervals)

    @classmethod
    def right_censored(cls, x, censored):
        """
        Create a `CensoredData` instance of right-censored data.

        Parameters
        ----------
        x : array_like
            `x` is the array of observed data or measurements.
            `x` must be a one-dimensional sequence of finite numbers.
        censored : array_like of bool
            `censored` must be a one-dimensional sequence of boolean
            values.  If ``censored[k]`` is True, the corresponding value
            in `x` is right-censored.  That is, the value ``x[k]``
            is the lower bound of the true (but unknown) value.

        Returns
        -------
        data : `CensoredData`
            An instance of `CensoredData` that represents the
            collection of uncensored and right-censored values.

        Examples
        --------
        >>> from scipy.stats import CensoredData

        Two uncensored values (4 and 10) and two right-censored values
        (24 and 25).

        >>> data = CensoredData.right_censored([4, 10, 24, 25],
        ...                                    [False, False, True, True])
        >>> data
        CensoredData(x=array([ 4., 10.]), left=array([], dtype=float64),
        right=array([24., 25.]), intervals=array([], shape=(0, 2),
        dtype=float64))
        >>> print(data)
        CensoredData(4 values: 2 not censored, 2 right-censored)
        """
        x, censored = _validate_x_censored(x, censored)
        return cls(x=x[~censored], right=x[censored])

    @classmethod
    def left_censored(cls, x, censored):
        """
        Create a `CensoredData` instance of left-censored data.

        Parameters
        ----------
        x : array_like
            `x` is the array of observed data or measurements.
            `x` must be a one-dimensional sequence of finite numbers.
        censored : array_like of bool
            `censored` must be a one-dimensional sequence of boolean
            values.  If ``censored[k]`` is True, the corresponding value
            in `x` is left-censored.  That is, the value ``x[k]``
            is the upper bound of the true (but unknown) value.

        Returns
        -------
        data : `CensoredData`
            An instance of `CensoredData` that represents the
            collection of uncensored and left-censored values.

        Examples
        --------
        >>> from scipy.stats import CensoredData

        Two uncensored values (0.12 and 0.033) and two left-censored values
        (both 1e-3).

        >>> data = CensoredData.left_censored([0.12, 0.033, 1e-3, 1e-3],
        ...                                    [False, False, True, True])
        >>> data
        CensoredData(x=array([0.12 , 0.033]), left=array([0.001, 0.001]),
        right=array([], dtype=float64), intervals=array([], shape=(0, 2),
        dtype=float64))
        >>> print(data)
        CensoredData(4 values: 2 not censored, 2 left-censored)
        """
        x, censored = _validate_x_censored(x, censored)
        return cls(x=x[~censored], left=x[censored])

    def _uncensor(self):
        """
        This function is used when a non-censored version of the data
        is needed to create a rough estimate of the parameters of a
        distribution via the method of moments or some similar method.
        The data is "uncensored" by taking the given endpoints as the
        data for the left- or right-censored data, and the mean for the
        interval-censored data.
        """
        data = np.concatenate((self._x, self._left, self._right,
                               self._intervals.mean(axis=1)))
        return data

    def _supported(self, a, b):
        """
        Return a subset of self containing the values that are in
        (or overlap with) the interval (a, b).
        """
        x = self._x
        x = x[(a < x) & (x < b)]
        left = self._left
        left = left[a < left]
        right = self._right
        right = right[right < b]
        intervals = self._intervals
        intervals = intervals[(a < intervals[:, 1]) & (intervals[:, 0] < b)]
        return CensoredData(x, left=left, right=right, intervals=intervals)
