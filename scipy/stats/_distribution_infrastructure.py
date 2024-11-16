import numpy as np

from scipy import special
from scipy.integrate._tanhsinh import _tanhsinh
from scipy.optimize._bracket import _bracket_minimum  # , _bracket_root
from scipy.optimize._chandrupatla import _chandrupatla_minimize  # , _chandrupatla
from scipy.stats._probability_distribution import (
    _ProbabilityDistribution, _isnull,
    _logexpxmexpy, _kwargs2args, _log_real_standardize,

)


__all__ = ['ContinuousDistribution']


_doc = r""" Class that represents a {continuous} statistical distribution.

Parameters
----------
tol : positive float, optional
    The desired relative tolerance of calculations. Left unspecified,
    calculations may be faster; when provided, calculations may be
    more likely to meet the desired accuracy.
validation_policy : {None, "skip_all"}
    Specifies the level of input validation to perform. Left unspecified,
    input validation is performed to ensure appropriate behavior in edge
    case (e.g. parameters out of domain, argument outside of distribution
    support, etc.) and improve consistency of output dtype, shape, etc.
    Pass ``'skip_all'`` to avoid the computational overhead of these
    checks when rough edges are acceptable.
cache_policy : {None, "no_cache"}
    Specifies the extent to which intermediate results are cached. Left
    unspecified, intermediate results of some calculations (e.g. distribution
    support, moments, etc.) are cached to improve performance of future
    calculations. Pass ``'no_cache'`` to reduce memory reserved by the class
    instance.

Attributes
----------
All parameters are available as attributes.

Methods
-------
support

plot

sample

moment

mean
median
mode

variance
standard_deviation

skewness
kurtosis

pdf
logpdf

cdf
icdf
ccdf
iccdf

logcdf
ilogcdf
logccdf
ilogccdf

entropy
logentropy

Notes
-----
The following abbreviations are used throughout the documentation.

- PDF: probability density function
- CDF: cumulative distribution function
- CCDF: complementary CDF
- entropy: differential entropy
- log-*F*: logarithm of *F* (e.g. log-CDF)
- inverse *F*: inverse function of *F* (e.g. inverse CDF)

The API documentation is written to describe the API, not to serve as
a statistical reference. Effort is made to be correct at the level
required to use the functionality, not to be mathematically rigorous.
For example, continuity and differentiability may be implicitly assumed.
For precise mathematical definitions, consult your preferred mathematical
text.

"""

class ContinuousDistribution(_ProbabilityDistribution):
    ## Algorithms

    def _quadrature(self, integrand, limits=None, args=None,
                    params=None, log=False):
        # Performs numerical integration of an integrand between limits.
        # Much of this should be added to `_tanhsinh`.
        a, b = self._support(**params) if limits is None else limits
        a, b = np.broadcast_arrays(a, b)
        if not a.size:
            # maybe need to figure out result type from a, b
            return np.empty(a.shape, dtype=self._dtype)
        args = [] if args is None else args
        params = {} if params is None else params
        f, args = _kwargs2args(integrand, args=args, kwargs=params)
        args = np.broadcast_arrays(*args)
        # If we know the median or mean, consider breaking up the interval
        rtol = None if _isnull(self.tol) else self.tol
        res = _tanhsinh(f, a, b, args=args, log=log, rtol=rtol)
        # For now, we ignore the status, but I want to return the error
        # estimate - see question 5 at the top.
        return res.integral

    ## Other

    def _overrides(self, method_name):
        # Determines whether a class overrides a specified method.
        # Returns True if the method implementation exists and is the same as
        # that of the `ContinuousDistribution` class; otherwise returns False.

        # Sometimes we use `_overrides` to check whether a certain method is overridden
        # and if so, call it. This begs the questions of why we don't do the more
        # obvious thing: restructure so that if the private method is overridden,
        # Python will call it instead of the inherited version automatically. The short
        # answer is that there are multiple ways a use might wish to evaluate a method,
        # and simply overriding the method with a formula is not always the best option.
        # For more complete discussion of the considerations, see:
        # https://github.com/scipy/scipy/pull/21050#discussion_r1707798901
        method = getattr(self.__class__, method_name, None)
        super_method = getattr(ContinuousDistribution, method_name, None)
        return method is not super_method

    ### Distribution Properties

    def _logentropy_quadrature(self, **params):
        def logintegrand(x, **params):
            logpdf = self._logpdf_dispatch(x, **params)
            return logpdf + np.log(0j+logpdf)
        res = self._quadrature(logintegrand, params=params, log=True)
        return _log_real_standardize(res + np.pi*1j)

    def _entropy_quadrature(self, **params):
        def integrand(x, **params):
            pdf = self._pdf_dispatch(x, **params)
            logpdf = self._logpdf_dispatch(x, **params)
            return logpdf * pdf
        return -self._quadrature(integrand, params=params)

    def _mode_optimization(self, **params):
        if not self._size:
            return np.empty(self._shape, dtype=self._dtype)

        a, b = self._support(**params)
        m = self._median_dispatch(**params)

        f, args = _kwargs2args(lambda x, **params: -self._pdf_dispatch(x, **params),
                               args=(), kwargs=params)
        res_b = _bracket_minimum(f, m, xmin=a, xmax=b, args=args)
        res = _chandrupatla_minimize(f, res_b.xl, res_b.xm, res_b.xr, args=args)
        mode = np.asarray(res.x)
        mode_at_boundary = res_b.status == -1
        mode_at_left = mode_at_boundary & (res_b.fl <= res_b.fm)
        mode_at_right = mode_at_boundary & (res_b.fr < res_b.fm)
        mode[mode_at_left] = a[mode_at_left]
        mode[mode_at_right] = b[mode_at_right]
        return mode[()]

    ## Cumulative Distribution Functions

    def _logcdf2_quadrature(self, x, y, **params):
        logres = self._quadrature(self._logpdf_dispatch, limits=(x, y),
                                  log=True, params=params)
        return logres

    def _logcdf_quadrature(self, x, **params):
        a, _ = self._support(**params)
        return self._quadrature(self._logpdf_dispatch, limits=(a, x),
                                params=params, log=True)

    def _cdf2_quadrature(self, x, y, **params):
        return self._quadrature(self._pdf_dispatch, limits=(x, y), params=params)

    def _cdf_quadrature(self, x, **params):
        a, _ = self._support(**params)
        return self._quadrature(self._pdf_dispatch, limits=(a, x),
                                params=params)

    def _logccdf2_addition(self, x, y, **params):
        logcdf_x = self._logcdf_dispatch(x, **params)
        logccdf_y = self._logccdf_dispatch(y, **params)
        return special.logsumexp([logcdf_x, logccdf_y], axis=0)

    def _logccdf_quadrature(self, x, **params):
        _, b = self._support(**params)
        return self._quadrature(self._logpdf_dispatch, limits=(x, b),
                                params=params, log=True)

    def _ccdf2_addition(self, x, y, **params):
        cdf_x = self._cdf_dispatch(x, **params)
        ccdf_y = self._ccdf_dispatch(y, **params)
        # even if x > y, cdf(x, y) + ccdf(x,y) sums to 1
        return cdf_x + ccdf_y

    def _ccdf_quadrature(self, x, **params):
        _, b = self._support(**params)
        return self._quadrature(self._pdf_dispatch, limits=(x, b),
                                params=params)

    ## Inverse cumulative distribution functions

    def _ilogcdf_inversion(self, x, **params):
        return self._solve_bounded(self._logcdf_dispatch, x, params=params)

    def _icdf_inversion(self, x, **params):
        return self._solve_bounded(self._cdf_dispatch, x, params=params)

    def _ilogccdf_inversion(self, x, **params):
        return self._solve_bounded(self._logccdf_dispatch, x, params=params)

    def _iccdf_inversion(self, x, **params):
        return self._solve_bounded(self._ccdf_dispatch, x, params=params)

    ### Sampling Functions

    def _sample_inverse_transform(self, sample_shape, full_shape, *, rng, **params):
        uniform = rng.random(size=full_shape, dtype=self._dtype)
        return self._icdf_dispatch(uniform, **params)

    ### Moments

    def _moment_integrate_pdf(self, order, center, **params):
        def integrand(x, order, center, **params):
            pdf = self._pdf_dispatch(x, **params)
            return pdf*(x-center)**order
        return self._quadrature(integrand, args=(order, center), params=params)

    def _moment_integrate_icdf(self, order, center, **params):
        def integrand(x, order, center, **params):
            x = self._icdf_dispatch(x, **params)
            return (x-center)**order
        return self._quadrature(integrand, limits=(0., 1.),
                                args=(order, center), params=params)

    def _logmoment(self, order=1, *, logcenter=None, standardized=False):
        # make this private until it is worked into moment
        if logcenter is None or standardized is True:
            logmean = self._logmoment_quad(self._one, -np.inf, **self._parameters)
        else:
            logmean = None

        logcenter = logmean if logcenter is None else logcenter
        res = self._logmoment_quad(order, logcenter, **self._parameters)
        if standardized:
            logvar = self._logmoment_quad(2, logmean, **self._parameters)
            res = res - logvar * (order/2)
        return res

    def _logmoment_quad(self, order, logcenter, **params):
        def logintegrand(x, order, logcenter, **params):
            logpdf = self._logpdf_dispatch(x, **params)
            return logpdf + order*_logexpxmexpy(np.log(x+0j), logcenter)
        return self._quadrature(logintegrand, args=(order, logcenter),
                                params=params, log=True)

    ### Convenience

    def plot(self, x='x', y='pdf', *, t=('cdf', 0.0005, 0.9995), ax=None):
        r"""Plot a function of the distribution.

        Convenience function for quick visualization of the distribution
        underlying the random variable.

        Parameters
        ----------
        x, y : str, optional
            String indicating the quantities to be used as the abscissa and
            ordinate (horizontal and vertical coordinates), respectively.
            Defaults are ``'x'`` (the domain of the random variable) and
            ``'pdf'`` (the probability density function). Valid values are:
            'x', 'pdf', 'cdf', 'ccdf', 'icdf', 'iccdf', 'logpdf', 'logcdf',
            'logccdf', 'ilogcdf', 'ilogccdf'.
        t : 3-tuple of (str, float, float), optional
            Tuple indicating the limits within which the quantities are plotted.
            Default is ``('cdf', 0.001, 0.999)`` indicating that the central
            99.9% of the distribution is to be shown. Valid values are:
            'x', 'cdf', 'ccdf', 'icdf', 'iccdf', 'logcdf', 'logccdf',
            'ilogcdf', 'ilogccdf'.
        ax : `matplotlib.axes`, optional
            Axes on which to generate the plot. If not provided, use the
            current axes.

        Returns
        -------
        ax : `matplotlib.axes`
            Axes on which the plot was generated.
            The plot can be customized by manipulating this object.

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Plot the PDF over the central 99.9% of the distribution.
        Compare against a histogram of a random sample.

        >>> ax = X.plot()
        >>> sample = X.sample(10000)
        >>> ax.hist(sample, density=True, bins=50, alpha=0.5)
        >>> plt.show()

        Plot ``logpdf(x)`` as a function of ``x`` in the left tail,
        where the log of the CDF is between -10 and ``np.log(0.5)``.

        >>> X.plot('x', 'logpdf', t=('logcdf', -10, np.log(0.5)))
        >>> plt.show()

        Plot the PDF of the normal distribution as a function of the
        CDF for various values of the scale parameter.

        >>> X = stats.Normal(mu=0., sigma=[0.5, 1., 2])
        >>> X.plot('cdf', 'pdf')
        >>> plt.show()

        """

        # Strategy: given t limits, get quantile limits. Form grid of
        # quantiles, compute requested x and y at quantiles, and plot.
        # Currently, the grid of quantiles is always linearly spaced.
        # Instead of always computing linearly-spaced quantiles, it
        # would be better to choose:
        # a) quantiles or probabilities
        # b) linearly or logarithmically spaced
        # based on the specified `t`.
        # TODO:
        # - smart spacing of points
        # - when the parameters of the distribution are an array,
        #   use the full range of abscissae for all curves

        t_is_quantile = {'x', 'icdf', 'iccdf', 'ilogcdf', 'ilogccdf'}
        t_is_probability = {'cdf', 'ccdf', 'logcdf', 'logccdf'}
        valid_t = t_is_quantile.union(t_is_probability)
        valid_xy =  valid_t.union({'pdf', 'logpdf'})

        ndim = self._ndim
        x_name, y_name = x, y
        t_name, tlim = t[0], np.asarray(t[1:])
        tlim = tlim[:, np.newaxis] if ndim else tlim

        # pdf/logpdf are not valid for `t` because we can't easily invert them
        message = (f'Argument `t` of `{self.__class__.__name__}.plot` "'
                   f'must be one of {valid_t}')
        if t_name not in valid_t:
            raise ValueError(message)

        message = (f'Argument `x` of `{self.__class__.__name__}.plot` "'
                   f'must be one of {valid_xy}')
        if x_name not in valid_xy:
            raise ValueError(message)

        message = (f'Argument `y` of `{self.__class__.__name__}.plot` "'
                   f'must be one of {valid_xy}')
        if t_name not in valid_xy:
            raise ValueError(message)

        # This could just be a warning
        message = (f'`{self.__class__.__name__}.plot` was called on a random '
                   'variable with at least one invalid shape parameters. When '
                   'a parameter is invalid, no plot can be shown.')
        if self._any_invalid:
            raise ValueError(message)

        # We could automatically ravel, but do we want to? For now, raise.
        message = ("To use `plot`, distribution parameters must be "
                   "scalars or arrays with one or fewer dimensions.")
        if ndim > 1:
            raise ValueError(message)

        try:
            import matplotlib.pyplot as plt  # noqa: F401, E402
        except ModuleNotFoundError as exc:
            message = ("`matplotlib` must be installed to use "
                       f"`{self.__class__.__name__}.plot`.")
            raise ModuleNotFoundError(message) from exc
        ax = plt.gca() if ax is None else ax

        # get quantile limits given t limits
        qlim = tlim if t_name in t_is_quantile else getattr(self, 'i'+t_name)(tlim)

        message = (f"`{self.__class__.__name__}.plot` received invalid input for `t`: "
                   f"calling {'i'+t_name}({tlim}) produced {qlim}.")
        if not np.all(np.isfinite(qlim)):
            raise ValueError(message)

        # form quantile grid
        grid = np.linspace(0, 1, 300)
        grid = grid[:, np.newaxis] if ndim else grid
        q = qlim[0] + (qlim[1] - qlim[0]) * grid

        # compute requested x and y at quantile grid
        x = q if x_name in t_is_quantile else getattr(self, x_name)(q)
        y = q if y_name in t_is_quantile else getattr(self, y_name)(q)

        # make plot
        ax.plot(x, y)
        ax.set_xlabel(f"${x_name}$")
        ax.set_ylabel(f"${y_name}$")
        ax.set_title(str(self))

        # only need a legend if distribution has parameters
        if len(self._parameters):
            label = []
            parameters = self._parameterization.parameters
            param_names = list(parameters)
            param_arrays = [np.atleast_1d(self._parameters[pname])
                            for pname in param_names]
            for param_vals in zip(*param_arrays):
                assignments = [f"${parameters[name].symbol}$ = {val:.4g}"
                               for name, val in zip(param_names, param_vals)]
                label.append(", ".join(assignments))
            ax.legend(label)

        return ax

    def llf(self, sample, /, *, axis=-1):
        r"""Log-likelihood function

        Given a sample :math:`x`, the log-likelihood function (LLF) is the logarithm
        of the joint probability density of the observed data. It is typically
        viewed as a function of the parameters :math:`\theta` of a statistical
        distribution:

        .. math::

            \mathcal{L}(\theta | x) = \log \left( \prod_i f_\theta(x_i) \right) = \sum_{i} \log(f_\theta(x_i))

        where :math:`f_\theta` is the probability density function with
        parameters :math:`\theta`.

        As a method of `ContinuousDistribution`, the parameter values are specified
        during instantiation; `llf` accepts only the sample :math:`x` as `sample`.

        Parameters
        ----------
        sample : array_like
            The given sample for which to calculate the LLF.
        axis : int or tuple of ints
            The axis over which the reducing operation (sum of logarithms) is performed.

        Notes
        -----
        The LLF is often viewed as a function of the parameters with the sample fixed;
        see the Notes for an example of a function with this signature.

        References
        ----------
        .. [1] Likelihood function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Likelihood_function

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from scipy import stats
        >>> X = stats.Normal(mu=0., sigma=1.)

        Evaluate the LLF with the given sample:

        >>> sample = [1., 2., 3.]
        >>> X.llf(sample)
        -9.756815599614018
        >>> np.allclose(X.llf(sample), np.sum(X.logpdf(sample)))
        True

        To generate a function that accepts only the parameters and
        holds the data fixed:

        >>> def llf(mu, sigma):
        ...     return stats.Normal(mu=mu, sigma=sigma).llf(sample)
        >>> llf(0., 1.)
        -9.756815599614018

        """ # noqa: E501
        return np.sum(self.logpdf(sample), axis=axis)


ContinuousDistribution.__doc__ = _doc.replace("{continuous}", 'continuous')
