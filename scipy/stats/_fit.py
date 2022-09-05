import warnings
from collections import namedtuple
import numpy as np
from scipy import optimize, stats
from scipy._lib._util import check_random_state
from scipy.stats import _stats_py

def _combine_bounds(name, user_bounds, shape_domain, integral):
    """Intersection of user-defined bounds and distribution PDF/PMF domain"""

    user_bounds = np.atleast_1d(user_bounds)

    if user_bounds[0] > user_bounds[1]:
        message = (f"There are no values for `{name}` on the interval "
                   f"{list(user_bounds)}.")
        raise ValueError(message)

    bounds = (max(user_bounds[0], shape_domain[0]),
              min(user_bounds[1], shape_domain[1]))

    if integral and (np.ceil(bounds[0]) > np.floor(bounds[1])):
        message = (f"There are no integer values for `{name}` on the interval "
                   f"defined by the user-provided bounds and the domain "
                   "of the distribution.")
        raise ValueError(message)
    elif not integral and (bounds[0] > bounds[1]):
        message = (f"There are no values for `{name}` on the interval "
                   f"defined by the user-provided bounds and the domain "
                   "of the distribution.")
        raise ValueError(message)

    if not np.all(np.isfinite(bounds)):
        message = (f"The intersection of user-provided bounds for `{name}` "
                   f"and the domain of the distribution is not finite. Please "
                   f"provide finite bounds for shape `{name}` in `bounds`.")
        raise ValueError(message)

    return bounds


class FitResult:
    r"""Result of fitting a discrete or continuous distribution to data

    Attributes
    ----------
    params : namedtuple
        A namedtuple containing the maximum likelihood estimates of the
        shape parameters, location, and (if applicable) scale of the
        distribution.
    success : bool or None
        Whether the optimizer considered the optimization to terminate
        successfully or not.
    message : str or None
        Any status message provided by the optimizer.

    """

    def __init__(self, dist, data, discrete, res):
        self._dist = dist
        self._data = data
        self.discrete = discrete
        self.pxf = getattr(dist, "pmf", None) or getattr(dist, "pdf", None)

        shape_names = [] if dist.shapes is None else dist.shapes.split(", ")
        if not discrete:
            FitParams = namedtuple('FitParams', shape_names + ['loc', 'scale'])
        else:
            FitParams = namedtuple('FitParams', shape_names + ['loc'])

        self.params = FitParams(*res.x)

        # Optimizer can report success even when nllf is infinite
        if res.success and not np.isfinite(self.nllf()):
            res.success = False
            res.message = ("Optimization converged to parameter values that "
                           "are inconsistent with the data.")
        self.success = getattr(res, "success", None)
        self.message = getattr(res, "message", None)

    def __repr__(self):
        keys = ["params", "success", "message"]
        m = max(map(len, keys)) + 1
        return '\n'.join([key.rjust(m) + ': ' + repr(getattr(self, key))
                          for key in keys if getattr(self, key) is not None])

    def nllf(self, params=None, data=None):
        """Negative log-likelihood function

        Evaluates the negative of the log-likelihood function of the provided
        data at the provided parameters.

        Parameters
        ----------
        params : tuple, optional
            The shape parameters, location, and (if applicable) scale of the
            distribution as a single tuple. Default is the maximum likelihood
            estimates (``self.params``).
        data : array_like, optional
            The data for which the log-likelihood function is to be evaluated.
            Default is the data to which the distribution was fit.

        Returns
        -------
        nllf : float
            The negative of the log-likelihood function.

        """
        params = params if params is not None else self.params
        data = data if data is not None else self._data
        return self._dist.nnlf(theta=params, x=data)

    def plot(self, ax=None):
        """Visualize the fit result.

        Superposes the PDF/PMF of the fitted distribution over a normalized
        histogram of the data.

        Available only if ``matplotlib`` is installed.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to draw the plot onto, otherwise uses the current Axes.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The matplotlib Axes object on which the plot was drawn.
        """
        try:
            from matplotlib.ticker import MaxNLocator
        except ModuleNotFoundError as exc:
            message = "matplotlib must be installed to use method `plot`."
            raise ValueError(message) from exc

        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()

        fit_params = np.atleast_1d(self.params)
        support = self._dist.support(*fit_params)
        lb = support[0] if np.isfinite(support[0]) else min(self._data)
        ub = support[1] if np.isfinite(support[1]) else max(self._data)

        if self.discrete:
            x = np.arange(lb, ub + 2)
            y = self.pxf(x, *fit_params)
            ax.vlines(x[:-1], 0, y[:-1], label='Fit Distribution PMF',
                      color='C0')
            options = dict(density=True, bins=x, align='left', color='C1')
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_xlabel('k')
            ax.set_ylabel('PMF')
        else:
            x = np.linspace(lb, ub, 200)
            y = self.pxf(x, *fit_params)
            ax.plot(x, y, '--', label='Fit Distribution PDF', color='C0')
            options = dict(density=True, bins=50, align='mid', color='C1')
            ax.set_xlabel('x')
            ax.set_ylabel('PDF')

        if len(self._data) > 50 or self.discrete:
            ax.hist(self._data, label="Histogram of Data", **options)
        else:
            ax.plot(self._data, np.zeros_like(self._data), "*",
                    label='Data', color='C1')

        ax.set_title(f"{self._dist.name} Fit")
        ax.legend(*ax.get_legend_handles_labels())
        return ax


def fit(dist, data, bounds=None, *, guess=None,
        optimizer=optimize.differential_evolution):
    r"""Fit a discrete or continuous distribution to data

    Given a distribution, data, and bounds on the parameters of the
    distribution, return maximum likelihood estimates of the parameters.

    Parameters
    ----------
    dist : `scipy.stats.rv_continuous` or `scipy.stats.rv_discrete`
        The object representing the distribution to be fit to the data.
    data : 1D array_like
        The data to which the distribution is to be fit. If the data contain
        any of ``np.nan``, ``np.inf``, or -``np.inf``, the fit method will
        raise a ``ValueError``.
    bounds : dict or sequence of tuples, optional
        If a dictionary, each key is the name of a parameter of the
        distribution, and the corresponding value is a tuple containing the
        lower and upper bound on that parameter.  If the distribution is
        defined only for a finite range of values of that parameter, no entry
        for that parameter is required; e.g., some distributions have
        parameters which must be on the interval [0, 1]. Bounds for parameters
        location (``loc``) and scale (``scale``) are optional; by default,
        they are fixed to 0 and 1, respectively.

        If a sequence, element *i* is a tuple containing the lower and upper
        bound on the *i*\ th parameter of the distribution. In this case,
        bounds for *all* distribution shape parameters must be provided.
        Optionally, bounds for location and scale may follow the
        distribution shape parameters.

        If a shape is to be held fixed (e.g. if it is known), the
        lower and upper bounds may be equal. If a user-provided lower or upper
        bound is beyond a bound of the domain for which the distribution is
        defined, the bound of the distribution's domain will replace the
        user-provided value. Similarly, parameters which must be integral
        will be constrained to integral values within the user-provided bounds.
    guess : dict or array_like, optional
        If a dictionary, each key is the name of a parameter of the
        distribution, and the corresponding value is a guess for the value
        of the parameter.

        If a sequence, element *i* is a guess for the *i*\ th parameter of the
        distribution. In this case, guesses for *all* distribution shape
        parameters must be provided.

        If `guess` is not provided, guesses for the decision variables will
        not be passed to the optimizer. If `guess` is provided, guesses for
        any missing parameters will be set at the mean of the lower and
        upper bounds. Guesses for parameters which must be integral will be
        rounded to integral values, and guesses that lie outside the
        intersection of the user-provided bounds and the domain of the
        distribution will be clipped.
    optimizer : callable, optional
        `optimizer` is a callable that accepts the following positional
        argument.

        fun : callable
            The objective function to be optimized. `fun` accepts one argument
            ``x``, candidate shape parameters of the distribution, and returns
            the negative log-likelihood function given ``x``, `dist`, and the
            provided `data`.
            The job of `optimizer` is to find values of the decision variables
            that minimizes `fun`.

        `optimizer` must also accepts the following keyword argument.

        bounds : sequence of tuples
            The bounds on values of the decision variables; each element will
            be a tuple containing the lower and upper bound on a decision
            variable.

        If `guess` is provided, `optimizer` must also accept the following
        keyword argument.

        x0 : array_like
            The guesses for each decision variable.

        If the distribution has any shape parameters that must be integral or
        if the distribution is discrete and the location parameter is not
        fixed, `optimizer` must also accept the following keyword argument.

        integrality : array_like of bools
            For each decision variable, True if the decision variable is
            must be constrained to integer values and False if the decision
            variable is continuous.

        `optimizer` must return an object, such as an instance of
        `scipy.optimize.OptimizeResult`, which holds the optimal values of
        the decision variables in an attribute ``x``. If attributes
        ``fun``, ``status``, or ``message`` are provided, they will be
        included in the result object returned by `fit`.

    Returns
    -------
    result : `~scipy.stats._result_classes.FitResult`
        An object with the following fields.

        params : namedtuple
            A namedtuple containing the maximum likelihood estimates of the
            shape parameters, location, and (if applicable) scale of the
            distribution.
        success : bool or None
            Whether the optimizer considered the optimization to terminate
            successfully or not.
        message : str or None
            Any status message provided by the optimizer.

        The object has the following method:

        nllf(params=None, data=None)
            By default, the negative log-likehood function at the fitted
            `params` for the given `data`. Accepts a tuple containing
            alternative shapes, location, and scale of the distribution and
            an array of alternative data.

        plot(ax=None)
            Superposes the PDF/PMF of the fitted distribution over a normalized
            histogram of the data.

    See Also
    --------
    rv_continuous,  rv_discrete

    Notes
    -----
    Optimization is more likely to converge to the maximum likelihood estimate
    when the user provides tight bounds containing the maximum likelihood
    estimate. For example, when fitting a binomial distribution to data, the
    number of experiments underlying each sample may be known, in which case
    the corresponding shape parameter ``n`` can be fixed.

    Examples
    --------
    Suppose we wish to fit a distribution to the following data.

    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng()
    >>> dist = stats.nbinom
    >>> shapes = (5, 0.5)
    >>> data = dist.rvs(*shapes, size=1000, random_state=rng)

    Suppose we do not know how the data were generated, but we suspect that
    it follows a negative binomial distribution with parameters *n* and *p*\.
    (See `scipy.stats.nbinom`.) We believe that the parameter *n* was fewer
    than 30, and we know that the parameter *p* must lie on the interval
    [0, 1]. We record this information in a variable `bounds` and pass
    this information to `fit`.

    >>> bounds = [(0, 30), (0, 1)]
    >>> res = stats.fit(dist, data, bounds)

    `fit` searches within the user-specified `bounds` for the
    values that best match the data (in the sense of maximum likelihood
    estimation). In this case, it found shape values similar to those
    from which the data were actually generated.

    >>> res.params
    FitParams(n=5.0, p=0.5028157644634368, loc=0.0)  # may vary

    We can visualize the results by superposing the probability mass function
    of the distribution (with the shapes fit to the data) over a normalized
    histogram of the data.

    >>> import matplotlib.pyplot as plt  # matplotlib must be installed to plot
    >>> res.plot()
    >>> plt.show()

    Note that the estimate for *n* was exactly integral; this is because
    the domain of the `nbinom` PMF includes only integral *n*, and the `nbinom`
    object "knows" that. `nbinom` also knows that the shape *p* must be a
    value between 0 and 1. In such a case - when the domain of the distribution
    with respect to a parameter is finite - we are not required to specify
    bounds for the parameter.

    >>> bounds = {'n': (0, 30)}  # omit parameter p using a `dict`
    >>> res2 = stats.fit(dist, data, bounds)
    >>> res2.params
    FitParams(n=5.0, p=0.5016492009232932, loc=0.0)  # may vary

    If we wish to force the distribution to be fit with *n* fixed at 6, we can
    set both the lower and upper bounds on *n* to 6. Note, however, that the
    value of the objective function being optimized is typically worse (higher)
    in this case.

    >>> bounds = {'n': (6, 6)}  # fix parameter `n`
    >>> res3 = stats.fit(dist, data, bounds)
    >>> res3.params
    FitParams(n=6.0, p=0.5486556076755706, loc=0.0)  # may vary
    >>> res3.nllf() > res.nllf()
    True  # may vary

    Note that the numerical results of the previous examples are typical, but
    they may vary because the default optimizer used by `fit`,
    `scipy.optimize.differential_evolution`, is stochastic. However, we can
    customize the settings used by the optimizer to ensure reproducibility -
    or even use a different optimizer entirely - using the `optimizer`
    parameter.

    >>> from scipy.optimize import differential_evolution
    >>> rng = np.random.default_rng(767585560716548)
    >>> def optimizer(fun, bounds, *, integrality):
    ...     return differential_evolution(fun, bounds, strategy='best2bin',
    ...                                   seed=rng, integrality=integrality)
    >>> bounds = [(0, 30), (0, 1)]
    >>> res4 = stats.fit(dist, data, bounds, optimizer=optimizer)
    >>> res4.params
    FitParams(n=5.0, p=0.5015183149259951, loc=0.0)

    """
    # --- Input Validation / Standardization --- #
    user_bounds = bounds
    user_guess = guess

    # distribution input validation and information collection
    if hasattr(dist, "pdf"):  # can't use isinstance for types
        default_bounds = {'loc': (0, 0), 'scale': (1, 1)}
        discrete = False
    elif hasattr(dist, "pmf"):
        default_bounds = {'loc': (0, 0)}
        discrete = True
    else:
        message = ("`dist` must be an instance of `rv_continuous` "
                   "or `rv_discrete.`")
        raise ValueError(message)

    try:
        param_info = dist._param_info()
    except AttributeError as e:
        message = (f"Distribution `{dist.name}` is not yet supported by "
                   "`scipy.stats.fit` because shape information has "
                   "not been defined.")
        raise ValueError(message) from e

    # data input validation
    data = np.asarray(data)
    if data.ndim != 1:
        message = "`data` must be exactly one-dimensional."
        raise ValueError(message)
    if not (np.issubdtype(data.dtype, np.number)
            and np.all(np.isfinite(data))):
        message = "All elements of `data` must be finite numbers."
        raise ValueError(message)

    # bounds input validation and information collection
    n_params = len(param_info)
    n_shapes = n_params - (1 if discrete else 2)
    param_list = [param.name for param in param_info]
    param_names = ", ".join(param_list)
    shape_names = ", ".join(param_list[:n_shapes])

    if user_bounds is None:
        user_bounds = {}

    if isinstance(user_bounds, dict):
        default_bounds.update(user_bounds)
        user_bounds = default_bounds
        user_bounds_array = np.empty((n_params, 2))
        for i in range(n_params):
            param_name = param_info[i].name
            user_bound = user_bounds.pop(param_name, None)
            if user_bound is None:
                user_bound = param_info[i].domain
            user_bounds_array[i] = user_bound
        if user_bounds:
            message = ("Bounds provided for the following unrecognized "
                       f"parameters will be ignored: {set(user_bounds)}")
            warnings.warn(message, RuntimeWarning, stacklevel=2)

    else:
        try:
            user_bounds = np.asarray(user_bounds, dtype=float)
            if user_bounds.size == 0:
                user_bounds = np.empty((0, 2))
        except ValueError as e:
            message = ("Each element of a `bounds` sequence must be a tuple "
                       "containing two elements: the lower and upper bound of "
                       "a distribution parameter.")
            raise ValueError(message) from e
        if (user_bounds.ndim != 2 or user_bounds.shape[1] != 2):
            message = ("Each element of `bounds` must be a tuple specifying "
                       "the lower and upper bounds of a shape parameter")
            raise ValueError(message)
        if user_bounds.shape[0] < n_shapes:
            message = (f"A `bounds` sequence must contain at least {n_shapes} "
                       "elements: tuples specifying the lower and upper "
                       f"bounds of all shape parameters {shape_names}.")
            raise ValueError(message)
        if user_bounds.shape[0] > n_params:
            message = ("A `bounds` sequence may not contain more than "
                       f"{n_params} elements: tuples specifying the lower and "
                       "upper bounds of distribution parameters "
                       f"{param_names}.")
            raise ValueError(message)

        user_bounds_array = np.empty((n_params, 2))
        user_bounds_array[n_shapes:] = list(default_bounds.values())
        user_bounds_array[:len(user_bounds)] = user_bounds

    user_bounds = user_bounds_array
    validated_bounds = []
    for i in range(n_params):
        name = param_info[i].name
        user_bound = user_bounds_array[i]
        param_domain = param_info[i].domain
        integral = param_info[i].integrality
        combined = _combine_bounds(name, user_bound, param_domain, integral)
        validated_bounds.append(combined)

    bounds = np.asarray(validated_bounds)
    integrality = [param.integrality for param in param_info]

    # guess input validation

    if user_guess is None:
        guess_array = None
    elif isinstance(user_guess, dict):
        default_guess = {param.name: np.mean(bound)
                         for param, bound in zip(param_info, bounds)}
        unrecognized = set(user_guess) - set(default_guess)
        if unrecognized:
            message = ("Guesses provided for the following unrecognized "
                       f"parameters will be ignored: {unrecognized}")
            warnings.warn(message, RuntimeWarning, stacklevel=2)
        default_guess.update(user_guess)

        message = ("Each element of `guess` must be a scalar "
                   "guess for a distribution parameter.")
        try:
            guess_array = np.asarray([default_guess[param.name]
                                      for param in param_info], dtype=float)
        except ValueError as e:
            raise ValueError(message) from e

    else:
        message = ("Each element of `guess` must be a scalar "
                   "guess for a distribution parameter.")
        try:
            user_guess = np.asarray(user_guess, dtype=float)
        except ValueError as e:
            raise ValueError(message) from e
        if user_guess.ndim != 1:
            raise ValueError(message)
        if user_guess.shape[0] < n_shapes:
            message = (f"A `guess` sequence must contain at least {n_shapes} "
                       "elements: scalar guesses for the distribution shape "
                       f"parameters {shape_names}.")
            raise ValueError(message)
        if user_guess.shape[0] > n_params:
            message = ("A `guess` sequence may not contain more than "
                       f"{n_params} elements: scalar guesses for the "
                       f"distribution parameters {param_names}.")
            raise ValueError(message)

        guess_array = np.mean(bounds, axis=1)
        guess_array[:len(user_guess)] = user_guess

    if guess_array is not None:
        guess_rounded = guess_array.copy()

        guess_rounded[integrality] = np.round(guess_rounded[integrality])
        rounded = np.where(guess_rounded != guess_array)[0]
        for i in rounded:
            message = (f"Guess for parameter `{param_info[i].name}` "
                       f"rounded from {guess_array[i]} to {guess_rounded[i]}.")
            warnings.warn(message, RuntimeWarning, stacklevel=2)

        guess_clipped = np.clip(guess_rounded, bounds[:, 0], bounds[:, 1])
        clipped = np.where(guess_clipped != guess_rounded)[0]
        for i in clipped:
            message = (f"Guess for parameter `{param_info[i].name}` "
                       f"clipped from {guess_rounded[i]} to "
                       f"{guess_clipped[i]}.")
            warnings.warn(message, RuntimeWarning, stacklevel=2)

        guess = guess_clipped
    else:
        guess = None

    # --- MLE Fitting --- #
    def nllf(free_params, data=data):  # bind data NOW
        with np.errstate(invalid='ignore', divide='ignore'):
            return dist._penalized_nnlf(free_params, data)

    with np.errstate(invalid='ignore', divide='ignore'):
        kwds = {}
        if bounds is not None:
            kwds['bounds'] = bounds
        if np.any(integrality):
            kwds['integrality'] = integrality
        if guess is not None:
            kwds['x0'] = guess
        res = optimizer(nllf, **kwds)

    return FitResult(dist, data, discrete, res)


GoodnessOfFitResult = namedtuple('GoodnessOfFitResult',
                                 ('fit_result', 'statistic', 'pvalue',
                                  'null_distribution'))


def goodness_of_fit(dist, data, *, known_params=None, fit_params=None,
                    guessed_params=None, statistic='ad', fit_method='mle',
                    n_resamples=9999, random_state=None):
    r"""
    Perform a goodness of fit test comparing data to a distribution family.

    Given a distribution family and data, perform a test of the null hypothesis
    that the data were drawn from a distribution in that family. Any known
    parameters of the distribution may be specified. Remaining parameters of
    the distribution will be fit to the data, and the p-value of the test
    is adjusted accordingly. Several statistics for comparing the distribution
    to data are available.

    Parameters
    ----------
    dist : `scipy.stats.rv_continuous`
        The object representing the distribution family to be fit to the data.
    data : 1D array_like
        Finite, uncensored data to which the distribution is to be fit.
    known_params : dict, optional
        A dictionary containing distribution name-value pairs of distribution
        parameters that are known. Data is resampled from the null-hypothesized
        distribution with these values of the parameters. Before the statistic
        is evaluated for each resample, only remaining unknown parameters of
        the  null-hypothesized distribution family are fit to the resampled
        data; these parameters are held fixed. (If all parameters of the
        distribution family are known, then the step of fitting the
        distribution family to each resample is omitted.)
    fit_params : dict, optional
        A dictionary containing name-value pairs of distribution parameters
        that have already been fit to the data, e.g. using `scipy.stats.fit`
        or the ``fit`` method of `dist`. Data is resampled from the
        null-hypothesized distribution with these specified values of the
        parameters, but these and all other unknown parameters of the
        null-hypothesized distribution family are fit to the resampled data
        before the statistic is evaluated.
    guessed_params : dict, optional
        A dictionary containing name-value pairs of distribution parameters
        which are guessed. Before resampling from the null-hypothesized
        distribution, the null-hypothesized distribution family is fit to the
        data with these parameters free, but the provided values will be used
        as guesses if needed by the fitting procedure. The null-hypothesized
        distribution family is fit to the resampled data with these parameters
        free before the statistic is evaluated.
    statistic : {"ad", "ks"}, optional
        The statistic used to compare data to a distribution after fitting
        unknown parameters of the distribution family to the data. The
        Anderson-Darling ("ad") [1]_ and Kolmogorov-Smirnov ("ks") [2]_
        statistics are available.
    fit_method : {"mle", "mm"}, optional
        The method used to fit null-hypothesized distribution family to
        the resampled data. Maximum likelihood estimation ("mle") is the
        default; method of moments ("mm") is also available
    n_resamples : int
        The number of resamples performed to form the null distribution
        of the statistic.
    random_state : {None, int, `numpy.random.Generator`,
                    `numpy.random.RandomState`}, optional

        Pseudorandom number generator state used to generate permutations.

        If `random_state` is ``None`` (default), the
        `numpy.random.RandomState` singleton is used.
        If `random_state` is an int, a new ``RandomState`` instance is used,
        seeded with `random_state`.
        If `random_state` is already a ``Generator`` or ``RandomState``
        instance then that instance is used.

    Returns
    -------
    res : GoodnessOfFitResult
        An object with the following attributes.

        fit_result : ~`scipy.stats._result_classes.FitResult`
            An object representing the fit of the provided `dist` to `data`.
            These are the parameters of the distribution family that define
            the specific distribution used to resample data under the null
            hypothesis.
        statistic : float
            The value of the statistic comparing provided `data` to the
            fitted distribution.
        null_distribution : ndarray
            The value of the statistic comparing data resampled under the null
            hypothesis to the fitted distribution.
        pvalue : float
            The proportion of elements in the null distribution with
            statistic values at least as extreme as the statistic value of the
            provided `data`.

    Notes
    -----
    This is a generalized Monte Carlo goodness-of-fit procedure, special cases
    of which correspond with various Anderson-Darling tests, Lilliefors' test,
    etc.

    Traditionally, critical values corresponding with a set a fixed set of
    signficance levels are pre-calculated using Monte Carlo methods. Users
    perform the test by calculating the value of the test statistic only for
    their observed `data` and comparing this value to the tabulated critical
    values. This practice is not very flexible, as tables are not available for
    all distributions and combinations of known and unknown parameter values,
    and results can be inaccurate when critical values are interpolated from
    limited tabulated data to correspond with the user's sample size and
    fitted parameter values. To overcome these shortcomings, this function
    allows the user to perform the Monte Carlo trials adapted to their
    particular data.

    First, any unknown parameters of the distribution family specified by
    `dist` are fit to the provided `data`. These values of the parameters
    specify a particular member of the distribution family referred to as the
    "null-hypothesized distribution", that is, the distribution from which the
    data were sampled under the null hypothesis. The `statistic`, which
    compares data to a distribution, is computed between `data` and the
    null-hypothesized distribution.

    Next, many (specifically `n_resamples`) new samples, each containing the
    same number of observations as `data`, are generated from the
    null-hypothesized distribution. Any unknown parameters of the distribution
    family `dist` are fit to *each resample*, and the `statistic` is computed
    between each resample and its corresponding fitted distribution. These
    values of the statistic form the Monte Carlo null distribution (not to be
    confused with the "null-hypothesized distribution" above).

    The p-value of the test is the proportion of statistic values in the Monte
    Carlo null distribution that are at least as extreme as the statistic value
    of the provided `data`. More precisely, the p-value is given by

    .. math::

        p = \frac{b + 1}
                 {m + 1}

    where :math:`b` is the number of statistic values in the Monte Carlo null
    distribution that are greater than or equal to the the statistic value for
    calculated for `data`, and :math:`m=` is the number of elements in the
    Monte Carlo null distribution (`n_resamples`). The addition of :math:`1`
    to the numerator and denominator can be thought of as including the
    value of the statistic corresponding with `data` in the null distribution,
    but a more formal explanation is given in [3]_.

    The test can be very slow for some distribution families because unknown
    parameters of the distribution family must be fit to each of the resamples,
    and for most distributions in SciPy, distribution fitting is performed via
    numerical optimization. For this reason, it may be tempting to treat
    parameters of the distribution pre-fit to `data` (by the user) as though
    they were `known_params`, as specification of all parameters of the
    distribution precludes the need to fit the distribution to each resample.
    (This is essentially how the original Kilmogorov-Smirnov test is
    performed.) Although such a test can provide evidence against the null
    hypothesis, the power of the test is low (that is, it is less likely to
    reject the null hypothesis even when the null hypoothesis is false) because
    the resampled data is less likely to agree with the null-hypothesized
    distribution as well as `data`. This tends to increase the values of the
    statistic recorded in the null distribution, so that a larger number of
    them exceed the value of statistic for `data`, inflating the p-value.

    References
    ----------
    .. [1] M. A. Stephens (1974). "EDF Statistics for Goodness of Fit and
           Some Comparisons." Journal of the American Statistical Association,
           Vol. 69, pp. 730-737.
    .. [2] F. J. Massey, Jr. (1951). "The Kolmogorov-Smirnov test for goodness
           of fit." Journal of the American Statistical Association 46.253:
           68-78.
    .. [3] B. Phipson and G. K. Smyth (2010). "Permutation P-values Should
           Never Be Zero: Calculating Exact P-values When Permutations Are
           Randomly Drawn." Statistical Applications in Genetics and Molecular
           Biology 9.1.
    .. [4] H. W. Lilliefors (1967). "On the Kolmogorov-Smirnov test for
           normality with mean and variance unknown." Journal of the American
           statistical Association 62.318: 399-402.

    Examples
    --------
    A well-known test of the null hypothesis that data were drawn from a
    given distribution is the Kolmogorov-Smirnov test, available in SciPy
    as `scipy.stats.ks_1samp`. Suppose we wish to test whether the following
    data:

    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng()
    >>> x = stats.uniform.rvs(size=75, random_state=rng)

    were sampled from a normal distribution. To perform a KS-test, the
    empirical distribution function of the observed data will be compared
    against the (theoretical) cumulative distribution function of a normal
    distribution. Of course, to do this, the normal distribution under the null
    hypothesis must be fully specified. This is commonly done by first fitting
    the ``loc`` and ``scale`` parameters of the distribution to the observed
    data, then performing the test.

    >>> loc, scale = np.mean(x), np.std(x, ddof=1)
    >>> cdf = stats.norm(loc, scale).cdf
    >>> stats.ks_1samp(x, cdf)
    KstestResult(statistic=0.1119257570456813, pvalue=0.2827756409939257)

    An advantage of the KS-test is that the p-value - the probability of
    obtaining a value of the test statistic under the null hypothesis as
    extreme as the value obtained from the observed data - can be calculated
    exactly (and quite quickly). `goodness_of_fit` can only approximate these
    results.

    >>> known_params = {'loc': loc, 'scale': scale}
    >>> res = stats.goodness_of_fit(stats.norm, x, known_params=known_params,
                                    statistic='ks', random_state=rng)
    >>> res.statistic, res.pvalue
    (0.1119257570456813, 0.2788)

    The statistic matches exactly, but the p-value is estimated by forming
    a "Monte Carlo null distribution", that is, by explicitly drawing random
    samples from `scipy.stats.norm` with the provided parameters and
    calculating the stastic for each. The fraction of these statistic values
    at least as extreme as ``res.statistic`` approximates the exact p-value
    calculated by `scipy.stats.ks_1samp`.

    However, in many cases, we would prefer to test only that the data were
    sampled from one of *any* member of the normal distribution family, not
    specifically from the normal distribution with the location and scale
    fitted to the observed sample. In this case, Lilliefors [4]_ argued that
    the KS test was far too conservative and thus lacked power - the ability
    to reject the null hypothesis when the null hypothesis is actually false.
    Indeed, our p-value above is approximately 0.28, which is far too large
    to reject the null hypothesis at any common significance level.

    Consider why this might be. Note that in the KS test above, the statistic
    always compares data against the CDF of a normal distribution fitted to the
    *observed data*. This tends to reduce the value of the statistic for the
    observed data, but it is "unfair" when computing the statistic for other
    samples, such as those we randomly draw to form the Monte Carlo null
    distribution. It is easy to correct for this: whenever we compute the KS
    statistic of a sample, we use the CDF of a normal distribution fitted
    to *that sample*. The null distribution in this case has not been
    calculated exactly and is tyically approximated using Monte Carlo methods
    as described above. This is where `goodness_of_fit` shines.

    >>> res = stats.goodness_of_fit(stats.norm, x, statistic='ks',
    ...                             random_state=rng)
    >>> res.statistic, res.pvalue
    (0.1119257570456813, 0.0196)

    Indeed, this p-value is much smaller, and small enough to reject the
    null hypothesis at common signficance levels, including 5% and 2.5%.

    However, the KS statistic is not very sensitive to all deviations from
    normality. The original advantage of the KS statistic was the ability
    to compute the null distribution theoeretically, but a more sensitive
    statistic can be used if we are willing to approximate the null
    distribution computationally. The Anderson-Darling statistic [1]_ tends
    to be more sensitive, and critical values of the this statistic have been
    tabulated for various significance levels and sample sizes using Monte
    Carlo methods.

    >>> res = stats.anderson(x, 'norm')
    >>> print(res.statistic)
    >>> print(res.critical_values)
    >>> print(res.significance_level)
    1.2139573337497467
    [0.549 0.625 0.75  0.875 1.041]
    [15.  10.   5.   2.5  1. ]

    Here, the observed value of the statistic exceeds the critical value
    corresponding with a 1% significance level. This tells us that the p-value
    of the observed data is less than 1%, but what is it? We could interpolate
    from these (altready interpolated) values, but `goodness_of_fit` can
    estimate it directly.

    >>> res = stats.goodness_of_fit(stats.norm, x, statistic='ad',
    ...                             random_state=rng)
    >>> res.statistic, res.pvalue
    (1.2139573337497467, 0.0034)

    A further advantage is that use of `goodness_of_fit` is not limited to
    a particular set of distributions or conditions on which parameters
    are known vs which must be estimated from data. Instead, `goodness_of_fit`
    can estimate p-values relatively quickly for any distribution with
    a sufficiently fast and reliable ``fit`` method. For instance,
    perform a goodness of fit test using the Cramer-von Mises statistic
    against the Rayleigh distribution with known location and unknown
    scale.

    >>> rng = np.random.default_rng()
    >>> x = stats.chi(df=2.2, loc=0, scale=2).rvs(size=1000, random_state=rng)
    >>> res = stats.goodness_of_fit(stats.rayleigh, x, statistic='cvm',
    ...                             known_params={'loc': 0}, random_state=rng)

    This executes fairly quickly, but to check the reliability of the ``fit``
    method, we should inspect the fit result.

    >>> import matplotlib.pyplot as plt  # matplotlib must be installed to plot
    >>> res.fit_result.plot()
    >>> plt.show()

    If the distribution is not fit to the observed data as well as possible,
    the test may not control the type I error rate, that is, the chance of
    rejecting the null hypothesis even when it is true.

    We should also look for extreme outliers in the null distribution that
    may be caused by unreliable fitting. These do not necessarily invalidate
    the result, but they tend to reduce the test's power.

    >>> _, ax = plt.subplots()
    >>> ax.hist(np.log10(res.null_distribution))
    >>> ax.set_xlabel("log10 of CVM statistic under the null hypothesis")
    >>> ax.set_ylabel("Frequency")
    >>> ax.set_title("Histogram of the Monte Carlo null distribution")
    >>> plt.show()

    If ``fit`` method is working reliably, and if the distribution of the test
    statistic is not particularly sensitive to the values of the fitted
    parameters, then the p-value provided by `goodness_of_fit` is expected to
    be a good approximation.
    >>> res.statistic, res.pvalue
    (0.06548669128345683, 0.5759)

    """
    args = gof_iv(dist, data, known_params, fit_params, guessed_params,
                  statistic, fit_method, n_resamples, random_state)
    (dist, data, fixed_nhd_params, fixed_rfd_params, guessed_nhd_params,
     guessed_rfd_params, statistic, fit_method, n_resamples_int,
     random_state) = args

    # Fit null hypothesis distribution to data
    nhd_fit_fun = _get_fit_fun(dist, data, fit_method,
                               guessed_nhd_params, fixed_nhd_params)
    nhd_vals = nhd_fit_fun(data)
    nhd_dist = dist(*nhd_vals)

    def rvs(size):
        return nhd_dist.rvs(size=size, random_state=random_state)

    # Define statistic
    fit_fun = _get_fit_fun(dist, data, fit_method,
                           guessed_rfd_params, fixed_rfd_params)
    compare_fun = _compare_dict[statistic]

    def statistic_fun(data, axis=-1):
        # Make things simple by always working along the last axis.
        data = np.moveaxis(data, axis, -1)
        rfd_vals = fit_fun(data)
        rfd_dist = dist(*rfd_vals)
        return compare_fun(rfd_dist, data)

    res = stats.monte_carlo_test(data, rvs, statistic_fun, vectorized=True,
                                 n_resamples=n_resamples, axis=-1,
                                 alternative='greater')
    opt_res = optimize.OptimizeResult()
    opt_res.success = True
    opt_res.message = "The fit was performed successfully."
    opt_res.x = nhd_vals
    # Only continuous distributions for now, hence discrete=False
    # There's no fundamental limitation; it's just that we're not using
    # stats.fit, discrete distributions don't have `fit` method, and
    # we haven't written any vectorized fit functions for a discrete
    # distribution yet.
    return GoodnessOfFitResult(FitResult(dist, data, False, opt_res),
                               res.statistic, res.pvalue,
                               res.null_distribution)


def _get_fit_fun(dist, data, fit_method, guessed_params, fixed_params):

    shape_names = [] if dist.shapes is None else dist.shapes.split(", ")
    param_names = shape_names + ['loc', 'scale']
    fparam_names = ['f'+name for name in param_names]
    all_fixed = not set(fparam_names).difference(fixed_params)
    guessed_shapes = [guessed_params.pop(x, None)
                      for x in shape_names if x in guessed_params]

    # Define statistic, including fitting distribution to data
    if dist in _fit_funs:
        def fit_fun(data):
            params = _fit_funs[dist](data, **fixed_params,
                                     fit_method=fit_method)
            params = np.asarray(np.broadcast_arrays(*params))
            if params.ndim > 1:
                params = params[..., np.newaxis]
            return params

    elif all_fixed:
        def fit_fun(data):
            return [fixed_params[name] for name in fparam_names]

    else:
        def fit_fun_1d(data):
            return dist.fit(data, *guessed_shapes, **guessed_params,
                            **fixed_params, method=fit_method)

        def fit_fun(data):
            params = np.apply_along_axis(fit_fun_1d, axis=-1, arr=data)
            if params.ndim > 1:
                params = params.T[..., np.newaxis]
            return params

    return fit_fun


# Vectorized fitting functions. These are to accept ND `data` in which each
# row (slice along last axis) is a sample to fit, scalar fixed parameters,
# and `fit_method`. They return a tuple of shape parameter arrays, each of
# shape data.shape[:-1].
def _fit_norm(data, floc=None, fscale=None, fit_method="mle"):
    loc = floc
    scale = fscale
    if loc is None and scale is None:
        loc = np.mean(data, axis=-1)
        scale = np.std(data, ddof=1, axis=-1)
    elif loc is None:
        loc = np.mean(data, axis=-1)
    elif scale is None:
        scale = np.sqrt(((data - loc)**2).mean(axis=-1))
    return loc, scale


_fit_funs = {stats.norm: _fit_norm}


# Vectorized goodness of fit statistic functions. These accept a frozen
# distribution object and `data` in which each row (slice along last axis) is
# a sample.


def _anderson_darling(dist, data):
    x = np.sort(data, axis=-1)
    n = data.shape[-1]
    i = np.arange(1, n+1)
    Si = (2*i - 1)/n * (dist.logcdf(x) + dist.logsf(x[..., ::-1]))
    S = np.sum(Si, axis=-1)
    return -n - S


def _kolmogorov_smirnov(dist, data):
    x = np.sort(data, axis=-1)
    cdfvals = dist.cdf(x)
    Dplus = _stats_py._compute_dplus(cdfvals, axis=-1)
    Dminus = _stats_py._compute_dminus(cdfvals, axis=-1)
    return np.maximum(Dplus, Dminus)


def _cramer_von_mises(dist, data):
    x = np.sort(data, axis=-1)
    n = data.shape[-1]
    cdfvals = dist.cdf(x)
    u = (2*np.arange(1, n+1) - 1)/(2*n)
    w = 1 / (12*n) + np.sum((u - cdfvals)**2, axis=-1)
    return w


_compare_dict = {"ad": _anderson_darling, "ks": _kolmogorov_smirnov,
                 "cvm": _cramer_von_mises}


def gof_iv(dist, data, known_params, fit_params, guessed_params,
           statistic, fit_method, n_resamples, random_state):

    if not isinstance(dist, stats.rv_continuous):
        message = ("`dist` must be a (non-frozen) instance of "
                   "`stats.rv_continuous`.")
        raise TypeError(message)

    data = np.asarray(data, dtype=float)
    if not data.ndim == 1:
        message = "`data` must be a one-dimensional array of numbers."
        raise ValueError(message)

    # Leave validation of these key/value pairs to the `fit` method,
    # but collect these into dictionaries that will be used
    known_params = known_params or dict()
    fit_params = fit_params or dict()
    guessed_params = guessed_params or dict()

    known_params_f = {("f"+key): val for key, val in known_params.items()}
    fit_params_f = {("f"+key): val for key, val in fit_params.items()}

    # These the the values of parameters of the null distribution family
    # with which resamples are drawn
    fixed_nhd_params = known_params_f.copy()
    fixed_nhd_params.update(fit_params_f)

    # These are fixed when fitting the distribution family to resamples
    fixed_rfd_params = known_params_f.copy()

    # These are used as guesses when fitting the distribution family to
    # the original data
    guessed_nhd_params = guessed_params.copy()

    # These are used as guesses when fitting the distribution family to
    # resamples
    guessed_rfd_params = fit_params.copy()
    guessed_rfd_params.update(guessed_params)

    statistics = {'ad', 'ks', 'cvm'}
    if statistic.lower() not in statistics:
        message = f"`statistic` must be one of {statistics}."
        raise ValueError(message)

    fit_methods = {'mle', 'mm'}
    if fit_method.lower() not in fit_methods:
        message = f"`fit_method` must be one of {fit_methods}."
        raise ValueError(message)

    n_resamples_int = int(n_resamples)
    if n_resamples_int != n_resamples:
        message = "`n_resamples` must be an integer."
        raise TypeError(message)

    random_state = check_random_state(random_state)

    return (dist, data, fixed_nhd_params, fixed_rfd_params, guessed_nhd_params,
            guessed_rfd_params, statistic, fit_method, n_resamples_int,
            random_state)
