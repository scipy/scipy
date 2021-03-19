import numpy as np


def _validate_lower_upper(lower, upper):
    if np.shape(lower) != np.shape(upper):
        raise ValueError('lower and upper must have the same shape.')
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
        raise ValueError('x must be one-dimensinal.')
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
    univariate SciPy distributions for MLE fitting.

    Left-, right-, and interval-censored data can be represented.

    For interval-censored data, pass the lower and upper limits of
    the intervals as one-dimensional arrays to the ``CensoredData``
    constuctor.  By using the values `np.inf` and `-np.inf`, left-
    and right-censored data may also be created with the constructor.

    For convenience, the class methods ``left_censored`` and
    ``right_censored`` are provided to create a ``CensoredData``
    instance from a single one-dimensional array of observations
    and a corresponding boolean array to indicate which observations
    are censored.

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

    The first value is left-censored; the true (but unknown) value is
    less than or equal to 0.  The second and third values, 1 and 1.5,
    are not censored.  The fourth value is interval-censored; the true
    value is in the interval [2, 3].  The last value is right-censored;
    the true value is greater than or equal to 10.

    A common case is to have data that is a mix of uncensored values and
    censored values that are all right-censored (or all left-censored)
    with the same bound on the censored values.  For example, a measuring
    device might have an upper limit of 100, so a reading of 100 means
    the true value is unknown, but the value is not less than 100.
    Suppose we have nine readings from such a device, with six being
    [65, 51, 88, 93, 96, 89] and three being 100.  To create an instance
    of `CensoredData` to represent this data, we can use the class method
    `CensoredData.right_censored` as follows:

    >>> values = [65, 51, 88, 93, 96, 89, 100, 100, 100]
    >>> censored = [0, 0, 0, 0, 0, 0, 1, 1, 1]

    A 1 (or any True value) in the sequence ``censored`` indicates
    that the corresponding value in ``values`` is right-censored.

    >>> data = CensoredData.right_censored(values, censored)
    >>> print(data)
    CensoredData(9 values: 6 not censored, 3 right-censored)
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
