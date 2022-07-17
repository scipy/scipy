import warnings
from collections import namedtuple
import numpy as np
from scipy import optimize


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
