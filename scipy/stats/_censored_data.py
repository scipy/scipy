import numpy as np


def _validate_lower_upper(lower, upper):
    if np.shape(lower) != np.shape(upper):
        raise ValueError('lower and upper must have the same shape.')
    if np.isnan(lower).any() or np.isnan(upper).any():
        raise ValueError('lower and upper must not contain nan.')
    if (np.isinf(lower) & np.isinf(upper)).any():
        raise ValueError('lower and upper must not both be infinity.')
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
            is the lower bound of the measurement.
        """
        x, censored = _validate_x_censored(x, censored)
        lower = 1.0*x  # Copy x while ensuring lower is floating point.
        upper = np.empty_like(lower)
        not_censored = ~censored
        upper[not_censored] = lower[not_censored]
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
            is the upper bound of the measurement.
        """
        x, censored = _validate_x_censored(x, censored)
        upper = 1.0*x  # Copy x while ensuring upper is floating point.
        lower = np.empty_like(upper)
        not_censored = ~censored
        lower[not_censored] = upper[not_censored]
        lower[censored] = -np.inf
        return cls(lower=lower, upper=upper)


def _uncensor(cdata):
    # This function is used when a non-censored version of the data
    # is needed to create a rough estimate of the parameters of a
    # distribution via the method of moments or some similar method.
    # The data is "uncensored" by taking the
    # given endpoints as the data for the left- or right-censored
    # data, and the mean for the interval-censored data.
    data = np.empty_like(cdata._lower)
    data[cdata._not_censored] = cdata._lower[cdata._not_censored]
    data[cdata._left_censored] = cdata._upper[cdata._left_censored]
    data[cdata._right_censored] = cdata._lower[cdata._right_censored]
    ic = cdata._interval_censored
    data[ic] = 0.5*(cdata._lower[ic] + cdata._upper[ic])
    return data
