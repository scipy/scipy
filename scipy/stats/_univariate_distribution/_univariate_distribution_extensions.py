"""Public functionality built on top of core univariate distribution classes.

May also include private helpers for the public functionality.
"""
import functools
import numpy as np
from math import inf

from scipy._lib._array_api import xp_promote, xp_copy, is_cupy, xp_device
from scipy._external import array_api_extra as xpx
from scipy import special, stats
from scipy.integrate import tanhsinh as _tanhsinh
from scipy.optimize import elementwise
from scipy.stats._univariate_distribution._probability_distribution import _ProbabilityDistribution
from scipy.stats._univariate_distribution._univariate_distribution_util import (
     _RealInterval, _IntegerInterval, _Parameterization, _guess_bracket,
     _RealParameter, _log_real_standardize,  rv_capabilities)
from scipy.stats._univariate_distribution._univariate_distribution import (
    UnivariateDistribution, ContinuousDistribution, DiscreteDistribution,
    CircularDistribution, TransformedDistribution, MonotonicTransformedDistribution,
    FoldedDistribution, _combine_docs)


__all__ = ['make_distribution', 'Mixture', 'order_statistic',
           'truncate', 'abs', 'exp', 'log']


# Special case the names of some new-style distributions in `make_distribution`
_distribution_names = {
    # Continuous
    'argus': 'ARGUS',
    'betaprime': 'BetaPrime',
    'chi2': 'ChiSquared',
    'crystalball': 'CrystalBall',
    'dgamma': 'DoubleGamma',
    'dweibull': 'DoubleWeibull',
    'expon': 'Exponential',
    'exponnorm': 'ExponentiallyModifiedNormal',
    'exponweib': 'ExponentialWeibull',
    'exponpow': 'ExponentialPower',
    'fatiguelife': 'FatigueLife',
    'foldcauchy': 'FoldedCauchy',
    'foldnorm': 'FoldedNormal',
    'genlogistic': 'GeneralizedLogistic',
    'gennorm': 'GeneralizedNormal',
    'genpareto': 'GeneralizedPareto',
    'genexpon': 'GeneralizedExponential',
    'genextreme': 'GeneralizedExtremeValue',
    'gausshyper': 'GaussHypergeometric',
    'gengamma': 'GeneralizedGamma',
    'genhalflogistic': 'GeneralizedHalfLogistic',
    'geninvgauss': 'GeneralizedInverseGaussian',
    'gumbel_r': 'Gumbel',
    'gumbel_l': 'ReflectedGumbel',
    'halfcauchy': 'HalfCauchy',
    'halflogistic': 'HalfLogistic',
    'halfnorm': 'HalfNormal',
    'halfgennorm': 'HalfGeneralizedNormal',
    'hypsecant': 'HyperbolicSecant',
    'invgamma': 'InverseGammma',
    'invgauss': 'InverseGaussian',
    'invweibull': 'InverseWeibull',
    'irwinhall': 'IrwinHall',
    'jf_skew_t': 'JonesFaddySkewT',
    'johnsonsb': 'JohnsonSB',
    'johnsonsu': 'JohnsonSU',
    'ksone': 'KSOneSided',
    'kstwo': 'KSTwoSided',
    'kstwobign': 'KSTwoSidedAsymptotic',
    'laplace_asymmetric': 'LaplaceAsymmetric',
    'levy_l': 'LevyLeft',
    'levy_stable': 'LevyStable',
    'loggamma': 'ExpGamma',  # really the Exponential Gamma Distribution
    'loglaplace': 'LogLaplace',
    'lognorm': 'LogNormal',
    'loguniform': 'LogUniform',
    'ncx2': 'NoncentralChiSquared',
    'nct': 'NoncentralT',
    'norm': 'Normal',
    'norminvgauss': 'NormalInverseGaussian',
    'powerlaw': 'PowerLaw',
    'powernorm': 'PowerNormal',
    'rdist': 'R',
    'rel_breitwigner': 'RelativisticBreitWigner',
    'recipinvgauss': 'ReciprocalInverseGaussian',
    'reciprocal': 'LogUniform',
    'semicircular': 'SemiCircular',
    'skewcauchy': 'SkewCauchy',
    'skewnorm': 'SkewNormal',
    'studentized_range': 'StudentizedRange',
    't': 'StudentT',
    'trapezoid': 'Trapezoidal',
    'triang': 'Triangular',
    'truncexpon': 'TruncatedExponential',
    'truncnorm': 'TruncatedNormal',
    'truncpareto': 'TruncatedPareto',
    'truncweibull_min': 'TruncatedWeibull',
    'tukeylambda': 'TukeyLambda',
    'vonmises_line': 'VonMisesLine',
    'weibull_min': 'Weibull',
    'weibull_max': 'ReflectedWeibull',
    'wrapcauchy': 'WrappedCauchyLine',
    # Discrete
    'betabinom': 'BetaBinomial',
    'betanbinom': 'BetaNegativeBinomial',
    'dlaplace': 'LaplaceDiscrete',
    'geom': 'Geometric',
    'hypergeom': 'Hypergeometric',
    'logser': 'LogarithmicSeries',
    'nbinom': 'NegativeBinomial',
    'nchypergeom_fisher': 'NoncentralHypergeometricFisher',
    'nchypergeom_wallenius': 'NoncentralHypergeometricWallenius',
    'nhypergeom': 'NegativeHypergeometric',
    'poisson_binom': 'PoissonBinomial',
    'randint': 'UniformDiscrete',
    'yulesimon': 'YuleSimon',
    'zipf': 'Zeta',
}


@rv_capabilities
def make_distribution(dist):
    """Generate a `UnivariateDistribution` class from a compatible object.

    The argument may be an instance of `rv_continuous`, `rv_discrete`, or an instance of
    another class that satisfies the interface described below.

    The returned value is a `ContinuousDistribution` subclass if the input defines a
    ``pdf`` or ``(i)(c)cdf`` method; it is a `DiscreteDistribution` subclass if
    the input defines a ``pmf`` method. Like any subclass of `UnivariateDistribution`,
    it must be instantiated (i.e. by passing all shape parameters as keyword arguments)
    before use. Once instantiated, the resulting object will have the same interface as
    any other instance of `UnivariateDistribution`; e.g., `scipy.stats.Normal`,
    `scipy.stats.Binomial`.

    .. note::

        `make_distribution` does not work perfectly with all instances of
        `rv_continuous` or `rv_discrete`. Known failures include `levy_stable`,
        `vonmises`, `hypergeom`, 'nchypergeom_fisher', 'nchypergeom_wallenius', and
        `poisson_binom`. Some methods of some distributions will not support
        array shape parameters.

    Parameters
    ----------
    dist : `rv_continuous`
        Instance of `rv_continuous`, `rv_discrete`, or an instance of any custom class
        with the following attributes:

        __make_distribution_version__ : str
            A string containing the version number of SciPy in which this interface
            is defined. The preferred interface may change in future SciPy versions,
            in which case support for an old interface version may be deprecated
            and eventually removed.
        parameters : dict or tuple
            If a dictionary, each key is the name of a parameter,
            and the corresponding value is either a dictionary or tuple.
            If the value is a dictionary, it may have the following items, with default
            values used for entries which aren't present.

            domain_type : {'continuous', 'discrete'}, default: 'continuous'
                A string identifying whether the domain of the parameter is continuous
                or discrete (i.e. accepts only integral values).

            endpoints : tuple, default: (-inf, inf)
                A tuple defining the lower and upper endpoints of the domain of the
                parameter; allowable values are floats, the name (string) of another
                parameter, or a callable taking parameters as keyword only
                arguments and returning the numerical value of an endpoint for
                given parameter values.

            inclusive : tuple of bool, default: (False, False)
                A tuple specifying whether the endpoints are included within the domain
                of the parameter.

            typical : tuple, default: ``endpoints``
                Defining endpoints of a typical range of values of a parameter. Can be
                used for sampling parameter values for testing. Behaves like the
                ``endpoints`` tuple above, and should define a subinterval of the
                domain given by ``endpoints``.

            A tuple value ``(a, b)`` associated to a key in the ``parameters``
            dictionary is equivalent to ``{endpoints: (a, b)}``.

            Custom distributions with multiple parameterizations can be defined by
            having the ``parameters`` attribute be a tuple of dictionaries with
            the structure described above. In this case, ``dist``\'s class must also
            define a method ``process_parameters`` to map between the different
            parameterizations. It must take all parameters from all parameterizations
            as optional keyword arguments and return a dictionary mapping parameters to
            values, filling in values from other parameterizations using values from
            the supplied parameterization. See example.

        support : dict or tuple
            A dictionary describing the support of the distribution or a tuple
            describing the endpoints of the support. This behaves identically to
            the values of the parameters dict described above. (``domain_type`` is
            inferred which methods are defined.)

        The class **must** also define at least one of the following methods:
        ``pdf``, ``pmf``, ``cdf``, ``ccdf``, ``icdf``,  and ``iccdf``.

        It **may not** not define both ``pdf`` and ``pmf``. Subject to that restriction,
        if the class defines ``pmf``, the support of the distribution is discrete (i.e.
        accepts only integral values); otherwise, the support is continuous.

        The class **may** also define methods
        ``logentropy``, ``entropy``, ``median``, ``mode``,
        ``logpdf``, ``logpmf``,
        ``logcdf``, ``logccdf``,
        ``ilogcdf``, ``ilogccdf``,
        ``moment``, ``lmoment``, and ``sample``.

        Defined methods must accept the parameters of the distribution as
        keyword arguments and also accept any positional-only arguments accepted by
        the corresponding method of `ContinuousDistribution`/`DiscreteDistribution`.
        When multiple parameterizations are defined, these methods must accept
        all parameters from all parameterizations. The ``moment`` method
        must accept the ``order`` and ``kind`` arguments by position or keyword, but
        may return ``None`` if a formula is not available for the arguments; in this
        case, the infrastructure will fall back to a default implementation.
        The ``lmoment`` method must accept the ``order`` argument by position or
        keyword and return the specified L-moment (not the L-moment ratio), but may
        return ``None`` if a formula is not available for the arguments. The
        ``sample`` method must accept ``shape`` by position or keyword, but contrary
        to the public method of the same name, the argument it receives will be the
        *full* shape of the output array - that is, the shape passed to the public
        method prepended to the broadcasted shape of random variable parameters.

    Returns
    -------
    CustomDistribution : `UnivariateDistribution`
        A subclass of `UnivariateDistribution` corresponding with `dist`. The
        initializer requires all shape parameters to be passed as keyword arguments
        (using the same names as the instance of `rv_continuous`/`rv_discrete`).

    Notes
    -----
    The documentation of `UnivariateDistribution` is not rendered. See below for
    an example of how to instantiate the class (i.e. pass all shape parameters of
    `dist` to the initializer as keyword arguments). Documentation of all methods
    is identical to that of `scipy.stats.Normal`. Use ``help`` on the returned
    class or its methods for more information.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats, special

    Create a `ContinuousDistribution` from `scipy.stats.loguniform`.

    >>> LogUniform = stats.make_distribution(stats.loguniform)
    >>> X = LogUniform(a=1.0, b=3.0)
    >>> np.isclose((X + 0.25).median(), stats.loguniform.ppf(0.5, 1, 3, loc=0.25))
    np.True_
    >>> X.plot()
    >>> sample = X.sample(10000, rng=np.random.default_rng())
    >>> plt.hist(sample, density=True, bins=30)
    >>> plt.legend(('pdf', 'histogram'))
    >>> plt.show()

    Create a custom continuous distribution.

    >>> class MyLogUniform:
    ...     @property
    ...     def __make_distribution_version__(self):
    ...         return "1.16.0"
    ...
    ...     @property
    ...     def parameters(self):
    ...         return {'a': {'endpoints': (0, np.inf),
    ...                       'inclusive': (False, False)},
    ...                 'b': {'endpoints': ('a', np.inf),
    ...                       'inclusive': (False, False)}}
    ...
    ...     @property
    ...     def support(self):
    ...         return {'endpoints': ('a', 'b'), 'inclusive': (True, True)}
    ...
    ...     def pdf(self, x, a, b):
    ...         return 1 / (x * (np.log(b)- np.log(a)))
    >>>
    >>> MyLogUniform = stats.make_distribution(MyLogUniform())
    >>> Y = MyLogUniform(a=1.0, b=3.0)
    >>> np.isclose(Y.cdf(2.), X.cdf(2.))
    np.True_

    Create a custom distribution with variable support.

    >>> class MyUniformCube:
    ...     @property
    ...     def __make_distribution_version__(self):
    ...         return "1.16.0"
    ...
    ...     @property
    ...     def parameters(self):
    ...         return {"a": (-np.inf, np.inf),
    ...                 "b": {'endpoints':('a', np.inf), 'inclusive':(True, False)}}
    ...
    ...     @property
    ...     def support(self):
    ...         def left(*, a, b):
    ...             return a**3
    ...
    ...         def right(*, a, b):
    ...             return b**3
    ...         return (left, right)
    ...
    ...     def pdf(self, x, *, a, b):
    ...         return 1 / (3*(b - a)*np.cbrt(x)**2)
    ...
    ...     def cdf(self, x, *, a, b):
    ...         return (np.cbrt(x) - a) / (b - a)
    >>>
    >>> MyUniformCube = stats.make_distribution(MyUniformCube())
    >>> X = MyUniformCube(a=-2, b=2)
    >>> Y = stats.Uniform(a=-2, b=2)**3
    >>> X.support()
    (-8.0, 8.0)
    >>> np.isclose(X.cdf(2.1), Y.cdf(2.1))
    np.True_

    Create a custom distribution with multiple parameterizations. Here we create a
    custom version of the beta distribution that has an alternative parameterization
    in terms of the mean ``mu`` and a dispersion parameter ``nu``.

    >>> class MyBeta:
    ...     @property
    ...     def __make_distribution_version__(self):
    ...         return "1.16.0"
    ...
    ...     @property
    ...     def parameters(self):
    ...         return ({"a": (0, np.inf), "b": (0, np.inf)},
    ...                 {"mu": (0, 1), "nu": (0, np.inf)})
    ...
    ...     def process_parameters(self, a=None, b=None, mu=None, nu=None):
    ...         if a is not None and b is not None:
    ...             nu = a + b
    ...             mu = a / nu
    ...         else:
    ...             a = mu * nu
    ...             b = nu - a
    ...         return dict(a=a, b=b, mu=mu, nu=nu)
    ...
    ...     @property
    ...     def support(self):
    ...         return {'endpoints': (0, 1)}
    ...
    ...     def pdf(self, x, a, b, mu, nu):
    ...         return special._ufuncs._beta_pdf(x, a, b)
    ...
    ...     def cdf(self, x, a, b, mu, nu):
    ...         return special.betainc(a, b, x)
    >>>
    >>> MyBeta = stats.make_distribution(MyBeta())
    >>> X = MyBeta(a=2.0, b=2.0)
    >>> Y = MyBeta(mu=0.5, nu=4.0)
    >>> np.isclose(X.pdf(0.3), Y.pdf(0.3))
    np.True_

    Create a custom discrete distribution with one continuous shape parameter and one
    discrete shape parameter.

    >>> class MyBinomial:
    ...     @property
    ...     def __make_distribution_version__(self):
    ...         return "1.16.0"
    ...
    ...     @property
    ...     def parameters(self):
    ...         return {'n': {'domain_type': 'discrete',
    ...                       'endpoints': (0., np.inf),
    ...                       'inclusive': (True, False)},
    ...                 'p': {'domain_type': 'continuous',
    ...                       'endpoints': (0., 1.),
    ...                       'inclusive': (True, True)}}
    ...
    ...     @property
    ...     def support(self):
    ...         return {'endpoints': (0, 'n'), 'inclusive': (True, True)}
    ...
    ...     def pmf(self, x, n, p):
    ...         return special.binom(n, x) * p**x * (1 - p)**(n-x)
    >>>
    >>> MyBinomial = stats.make_distribution(MyBinomial())
    >>> X = stats.Binomial(n=10, p=0.4)
    >>> Y = MyBinomial(n=10, p=0.4)
    >>> np.isclose(Y.cdf(8.), X.cdf(8.))
    np.True_

    Create a quantile-defined distribution; i.e., one defined by the inverse-CDF rather
    then PDF or PMF.

    >>> class MyTukeyLambda:
    ...     __make_distribution_version__ = "1.16.0"
    ...
    ...     @property
    ...     def parameters(self):
    ...         return {'lam': {'endpoints': (-np.inf, np.inf),
    ...                        'inclusive': (False, False)}}
    ...
    ...     @property
    ...     def support(self):
    ...         return {'endpoints': (lambda *, lam: np.where(lam > 0, -1/lam, -np.inf),
    ...                               lambda *, lam: np.where(lam > 0, 1/lam, np.inf)),
    ...                'inclusive': (True, True)}
    ...
    ...     def icdf(self, p, lam):
    ...         return 1/lam * (p**lam - (1 - p)**lam)
    >>>
    >>> MyTukeyLambda = stats.make_distribution(MyTukeyLambda())
    >>> X = stats.tukeylambda(lam=0.14)
    >>> Y = MyTukeyLambda(lam=0.14)
    >>> np.isclose(Y.pdf(-0.314), X.pdf(-0.314))
    True

    """
    if dist in {stats.levy_stable, stats.vonmises, stats.hypergeom,
                stats.nchypergeom_fisher, stats.nchypergeom_wallenius,
                stats.poisson_binom}:
        raise NotImplementedError(f"`{dist.name}` is not supported.")

    if isinstance(dist, stats.rv_continuous | stats.rv_discrete):
        return _make_distribution_rv_generic(dist)
    elif getattr(dist, "__make_distribution_version__", "0.0.0") >= "1.16.0":
        return _make_distribution_custom(dist)
    else:
        message = ("The argument must be an instance of `rv_continuous`, "
                   "`rv_discrete`, or an instance of a class with attribute "
                   "`__make_distribution_version__ >= 1.16`.")
        raise ValueError(message)


def _make_distribution_rv_generic(dist):
    parameters = []
    names = []
    support = getattr(dist, '_support', (dist.a, dist.b))
    for shape_info in dist._shape_info():
        domain = _RealInterval(endpoints=shape_info.endpoints,
                               inclusive=shape_info.inclusive)
        param = _RealParameter(shape_info.name, domain=domain)
        parameters.append(param)
        names.append(shape_info.name)

    repr_str = _distribution_names.get(dist.name, dist.name.capitalize())
    if isinstance(dist, stats.rv_continuous):
        old_class, new_class = stats.rv_continuous, ContinuousDistribution
    else:
        old_class, new_class = stats.rv_discrete, DiscreteDistribution

    def _overrides(method_name):
        return (getattr(dist.__class__, method_name, None)
                is not getattr(old_class, method_name, None))

    if _overrides("_get_support"):
        def left(**parameter_values):
            a, _ = dist._get_support(**parameter_values)
            return np.asarray(a)[()]

        def right(**parameter_values):
            _, b = dist._get_support(**parameter_values)
            return np.asarray(b)[()]

        endpoints = (left, right)
    else:
        endpoints = support

    _x_support = _RealInterval(endpoints=endpoints, inclusive=(True, True))
    _x_param = _RealParameter('x', domain=_x_support, typical=(-1, 1))

    class CustomDistribution(new_class):
        _parameterizations = ([_Parameterization(*parameters)] if parameters
                              else [])
        _variable = _x_param

        __class_getitem__ = None

        def __repr__(self):
            s = super().__repr__()
            return s.replace('CustomDistribution', repr_str)

        def __str__(self):
            s = super().__str__()
            return s.replace('CustomDistribution', repr_str)

    def _sample_formula(self, full_shape=(), *, rng=None, **kwargs):
        return dist._rvs(size=full_shape, random_state=rng, **kwargs)

    def _moment_raw_formula(self, order, **kwargs):
        return dist._munp(int(order), **kwargs)

    def _moment_raw_formula_1(self, order, **kwargs):
        if order != 1:
            return None
        return dist._stats(**kwargs)[0]

    def _moment_central_formula(self, order, **kwargs):
        if order != 2:
            return None
        return dist._stats(**kwargs)[1]

    def _moment_standard_formula(self, order, **kwargs):
        if order == 3:
            if dist._stats_has_moments:
                kwargs['moments'] = 's'
            return dist._stats(**kwargs)[int(order - 1)]
        elif order == 4:
            if dist._stats_has_moments:
                kwargs['moments'] = 'k'
            k = dist._stats(**kwargs)[int(order - 1)]
            return k if k is None else k + 3
        else:
            return None

    methods = {'_logpdf': '_logpdf_formula',
               '_pdf': '_pdf_formula',
               '_logpmf': '_logpmf_formula',
               '_pmf': '_pmf_formula',
               '_logcdf': '_logcdf_formula',
               '_cdf': '_cdf_formula',
               '_logsf': '_logccdf_formula',
               '_sf': '_ccdf_formula',
               '_ppf': '_icdf_formula',
               '_isf': '_iccdf_formula',
               '_entropy': '_entropy_formula',
               '_median': '_median_formula'}

    # These are not desirable overrides for the new infrastructure
    skip_override = {'norminvgauss': {'_sf', '_isf'}}

    for old_method, new_method in methods.items():
        if dist.name in skip_override and old_method in skip_override[dist.name]:
            continue
        # If method of old distribution overrides generic implementation...
        method = getattr(dist.__class__, old_method, None)
        super_method = getattr(old_class, old_method, None)
        if method is not super_method:
            # Make it an attribute of the new object with the new name
            setattr(CustomDistribution, new_method, getattr(dist, old_method))

    if _overrides('_munp'):
        CustomDistribution._moment_raw_formula = _moment_raw_formula

    if _overrides('_rvs'):
        CustomDistribution._sample_formula = _sample_formula

    if _overrides('_stats'):
        CustomDistribution._moment_standardized_formula = _moment_standard_formula
        if not _overrides('_munp'):
            CustomDistribution._moment_raw_formula = _moment_raw_formula_1
            CustomDistribution._moment_central_formula = _moment_central_formula

    support_etc = _combine_docs(CustomDistribution, include_examples=False).lstrip()
    docs = [
        f"This class represents `scipy.stats.{dist.name}` as a subclass of "
        f"`{new_class}`.",
        f"The `repr`/`str` of class instances is `{repr_str}`.",
        f"The PDF of the distribution is defined {support_etc}"
    ]
    CustomDistribution.__doc__ = ("\n".join(docs))

    return CustomDistribution


def _get_domain_info(info):
    domain_info = {"endpoints": info} if isinstance(info, tuple) else info
    domain_type_str = domain_info.pop("domain_type", "continuous")
    if domain_type_str not in {"discrete", "continuous"}:
        message = ("If specified, the `domain_type` value of a parameter "
                   "must be either 'continuous' or 'discrete'.")
        raise ValueError(message)
    domain_type = _IntegerInterval if domain_type_str == 'discrete' else _RealInterval
    typical = domain_info.pop("typical", None)
    return domain_type, domain_info, typical


def _make_distribution_custom(dist):
    dist_parameters = (
        dist.parameters if isinstance(dist.parameters, tuple) else (dist.parameters, )
    )
    parameterizations = []
    for parameterization in dist_parameters:
        # The attribute name ``parameters`` appears reasonable from a user facing
        # perspective, but there is a little tension here with the internal. It's
        # important to keep in mind that the ``parameters`` attribute in a
        # user-created custom distribution specifies ``_parameterizations`` within
        # the infrastructure.
        parameters = []

        for name, info in parameterization.items():
            domain_type, domain_info, typical = _get_domain_info(info)
            domain = domain_type(**domain_info)
            param = _RealParameter(name, domain=domain, typical=typical)
            parameters.append(param)

        if parameters:
            parameterizations.append(_Parameterization(*parameters))

    domain_type, domain_info, typical = _get_domain_info(dist.support)
    _x_support = domain_type(**domain_info)

    if not hasattr(dist, 'pmf') and (hasattr(dist, 'pdf') or
                                     hasattr(dist, 'cdf') or
                                     hasattr(dist, 'ccdf') or
                                     hasattr(dist, 'icdf') or
                                     hasattr(dist, 'iccdf')):
        pxf = 'PDF'
        distribution_subclass = ContinuousDistribution
    elif hasattr(dist, 'pmf') and not hasattr(dist, 'pdf'):
        pxf = 'PMF'
        distribution_subclass = DiscreteDistribution
    else:
        message = ("The argument of `make_distribution` must implement at least one "
                   "of `pdf`/`pmf`/`cdf`/`ccdf`/`icdf`/`iccdf`, and may not implement "
                   "both `pdf` and `pmf`.")
        raise ValueError(message)

    _x_param = _RealParameter('x', domain=_x_support, typical=typical)
    repr_str = dist.__class__.__name__

    class CustomDistribution(distribution_subclass):
        _parameterizations = parameterizations
        _variable = _x_param

        def __repr__(self):
            s = super().__repr__()
            return s.replace('CustomDistribution', repr_str)

        def __str__(self):
            s = super().__str__()
            return s.replace('CustomDistribution', repr_str)

    methods = {'sample', 'logentropy', 'entropy',
               'median', 'mode', 'logpdf', 'pdf', 'logpmf', 'pmf',
               'logcdf2', 'logcdf', 'cdf2', 'cdf',
               'logccdf2', 'logccdf', 'ccdf2', 'ccdf',
               'ilogcdf', 'icdf', 'ilogccdf', 'iccdf',
               'lmoment'}

    for method in methods:
        if hasattr(dist, method):
            # Make it an attribute of the new object with the new name
            new_method = f"_{method}_formula"
            setattr(CustomDistribution, new_method, getattr(dist, method))

    if hasattr(dist, 'moment'):
        def _moment_raw_formula(self, order, **kwargs):
            return dist.moment(order, kind='raw', **kwargs)

        def _moment_central_formula(self, order, **kwargs):
            return dist.moment(order, kind='central', **kwargs)

        def _moment_standardized_formula(self, order, **kwargs):
            return dist.moment(order, kind='standardized', **kwargs)

        CustomDistribution._moment_raw_formula = _moment_raw_formula
        CustomDistribution._moment_central_formula = _moment_central_formula
        CustomDistribution._moment_standardized_formula = _moment_standardized_formula

    if hasattr(dist, 'process_parameters'):
        setattr(
            CustomDistribution,
            "_process_parameters",
            getattr(dist, "process_parameters")
        )

    support_etc = _combine_docs(CustomDistribution, include_examples=False).lstrip()
    docs = [
        f"This class represents `{repr_str}` as a subclass of "
        f"`{str(distribution_subclass)}`.",
        f"The {pxf} of the distribution is defined {support_etc}"
    ]
    CustomDistribution.__doc__ = ("\n".join(docs))

    return CustomDistribution


class TruncatedDistribution(TransformedDistribution):
    """Truncated distribution."""
    # TODO:
    # - consider avoiding catastrophic cancellation by using appropriate tail
    # - if the mode of `_dist` is within the support, it's still the mode
    # - rejection sampling might be more efficient than inverse transform

    _lb_domain = _RealInterval(endpoints=(-inf, 'ub'), inclusive=(True, False))
    _lb_param = _RealParameter('lb', symbol=r'b_l',
                                domain=_lb_domain, typical=(0.1, 0.2))

    _ub_domain = _RealInterval(endpoints=('lb', inf), inclusive=(False, True))
    _ub_param = _RealParameter('ub', symbol=r'b_u',
                                  domain=_ub_domain, typical=(0.8, 0.9))

    _parameterizations = [_Parameterization(_lb_param, _ub_param),
                          _Parameterization(_lb_param),
                          _Parameterization(_ub_param)]

    def __init__(self, X, /, *args, lb=-inf, ub=inf, **kwargs):
        return super().__init__(X, *args, lb=lb, ub=ub, **kwargs)

    def _process_parameters(self, lb=None, ub=None, **params):
        xp = self._xp
        lb = lb if lb is not None else xp.full_like(lb, -inf)[()]
        ub = ub if ub is not None else xp.full_like(ub, inf)[()]
        parameters = self._dist._process_parameters(**params)
        a, b = self._support(lb=lb, ub=ub, **parameters)
        logmass = self._dist._logcdf2_dispatch(a, b, **parameters)
        parameters.update(dict(lb=lb, ub=ub, _a=a, _b=b, logmass=logmass))
        return parameters

    def _support(self, lb, ub, **params):
        xp = self._xp
        a, b = self._dist._support(**params)
        return xp.maximum(a, lb), xp.minimum(b, ub)

    def _overrides(self, method_name):
        return method_name == '_logpdf_dispatch'

    def _logpdf_dispatch(self, x, *args, lb, ub, _a, _b, logmass, **params):
        logpdf = self._dist._logpdf_dispatch(x, *args, **params)
        return logpdf - logmass

    def _logcdf_dispatch(self, x, *args, lb, ub, _a, _b, logmass, **params):
        logcdf = self._dist._logcdf2_dispatch(_a, x, *args, **params)
        # of course, if this result is small we could compute with the other tail
        return logcdf - logmass

    def _logccdf_dispatch(self, x, *args, lb, ub, _a, _b, logmass, **params):
        logccdf = self._dist._logcdf2_dispatch(x, _b, *args, **params)
        return logccdf - logmass

    def _logcdf2_dispatch(self, x, y, *args, lb, ub, _a, _b, logmass, **params):
        logcdf2 = self._dist._logcdf2_dispatch(x, y, *args, **params)
        return logcdf2 - logmass

    def _ilogcdf_dispatch(self, logp, *args, lb, ub, _a, _b, logmass, **params):
        xp = self._xp
        log_Fa = self._dist._logcdf_dispatch(_a, *args, **params)
        logp_adjusted = xp.logaddexp(log_Fa, logp + logmass)
        return self._dist._ilogcdf_dispatch(logp_adjusted, *args, **params)

    def _ilogccdf_dispatch(self, logp, *args, lb, ub, _a, _b, logmass, **params):
        xp = self._xp
        log_cFb = self._dist._logccdf_dispatch(_b, *args, **params)
        logp_adjusted = xp.logaddexp(log_cFb, logp + logmass)
        return self._dist._ilogccdf_dispatch(logp_adjusted, *args, **params)

    def _icdf_dispatch(self, p, *args, lb, ub, _a, _b, logmass, **params):
        xp = self._xp
        Fa = self._dist._cdf_dispatch(_a, *args, **params)
        p_adjusted = Fa + p*xp.exp(logmass)
        return self._dist._icdf_dispatch(p_adjusted, *args, **params)

    def _iccdf_dispatch(self, p, *args, lb, ub, _a, _b, logmass, **params):
        xp = self._xp
        cFb = self._dist._ccdf_dispatch(_b, *args, **params)
        p_adjusted = cFb + p*xp.exp(logmass)
        return self._dist._iccdf_dispatch(p_adjusted, *args, **params)

    def __repr__(self):
        with np.printoptions(threshold=10):
            return (f"truncate({repr(self._dist)}, "
                    f"lb={repr(self.lb)}, ub={repr(self.ub)})")

    def __str__(self):
        with np.printoptions(threshold=10):
            return (f"truncate({str(self._dist)}, "
                    f"lb={str(self.lb)}, ub={str(self.ub)})")


@rv_capabilities
def truncate(X, lb=-inf, ub=inf):
    """Truncate the support of a random variable.

    Given a random variable `X`, `truncate` returns a random variable with
    support truncated to the interval between `lb` and `ub`. The underlying
    probability density function is normalized accordingly.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable to be truncated.
    lb, ub : float array-like
        The lower and upper truncation points, respectively. Must be
        broadcastable with one another and the shape of `X`.

    Returns
    -------
    X : `ContinuousDistribution`
        The truncated random variable.

    References
    ----------
    .. [1] "Truncated Distribution". *Wikipedia*.
           https://en.wikipedia.org/wiki/Truncated_distribution

    Examples
    --------
    Compare against `scipy.stats.truncnorm`, which truncates a standard normal,
    *then* shifts and scales it.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> loc, scale, lb, ub = 1, 2, -2, 2
    >>> X = stats.truncnorm(lb, ub, loc, scale)
    >>> Y = scale * stats.truncate(stats.Normal(), lb, ub) + loc
    >>> x = np.linspace(-3, 5, 300)
    >>> plt.plot(x, X.pdf(x), '-', label='X')
    >>> plt.plot(x, Y.pdf(x), '--', label='Y')
    >>> plt.xlabel('x')
    >>> plt.ylabel('PDF')
    >>> plt.title('Truncated, then Shifted/Scaled Normal')
    >>> plt.legend()
    >>> plt.show()

    However, suppose we wish to shift and scale a normal random variable,
    then truncate its support to given values. This is straightforward with
    `truncate`.

    >>> Z = stats.truncate(scale * stats.Normal() + loc, lb, ub)
    >>> Z.plot()
    >>> plt.show()

    Furthermore, `truncate` can be applied to any random variable:

    >>> Rayleigh = stats.make_distribution(stats.rayleigh)
    >>> W = stats.truncate(Rayleigh(), lb=0.5, ub=3)
    >>> W.plot()
    >>> plt.show()

    """
    return TruncatedDistribution(X, lb=lb, ub=ub)


class OrderStatisticDistribution(TransformedDistribution):
    r"""Probability distribution of an order statistic.

    An instance of this class represents a random variable that follows the
    distribution of the :math:`r^{\text{th}}` order statistic of a
    sample of :math:`n` observations of a random variable :math:`X`.

    Parameters
    ----------
    dist : `ContinuousDistribution`
        The random variable :math:`X`.
    n : array_like
        The (positive integer) sample size :math:`n`.
    r : array_like
        The (positive integer) rank of the order statistic :math:`r`,
        satisfying :math:`1 \leq r \leq n`.


    Notes
    -----
    If we make :math:`n` observations of a continuous random variable
    :math:`X` and sort them in increasing order
    :math:`X_{(1)}, \dots, X_{(r)}, \dots, X_{(n)}`, then :math:`X_{(r)}`
    is known as the :math:`r^{\text{th}}` order statistic.

    If the PDF, CDF, and CCDF of :math:`X` are denoted by :math:`f`,
    :math:`F`, and :math:`G = 1 - F`, respectively, then the PDF of
    :math:`X_{(r)}` is given by:

    .. math::

        f_r(x) = \frac{n!}{(r-1)! (n-r)!} f(x) F(x)^{r-1} G(x)^{n - r}

    The CDF and other methods of the distribution of :math:`X_{(r)}`
    are calculated using the fact that :math:`X = F^{-1}(U)`, where :math:`U` is a
    standard uniform random variable, together with the fact that the order statistics
    of i.i.d. uniform random variables follow a beta distribution
    :math:`B(r, n - r + 1)`.

    References
    ----------
    .. [1] Order statistic. *Wikipedia*. https://en.wikipedia.org/wiki/Order_statistic

    Examples
    --------
    Suppose we are interested in order statistics of samples of size five drawn
    from the standard normal distribution. Plot the PDF of the fourth
    order statistic and compare with a normalized histogram from simulation.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> from scipy.stats._univariate_distribution import OrderStatisticDistribution
    >>>
    >>> X = stats.Normal()
    >>> data = X.sample(shape=(10000, 5))
    >>> sorted_data = np.sort(data, axis=1)
    >>> Y = OrderStatisticDistribution(X, r=4, n=5)
    >>>
    >>> ax = plt.gca()
    >>> Y.plot(ax=ax)
    >>> ax.hist(sorted_data[:, 3], density=True, bins=30)
    >>> plt.show()

    """

    _r_domain = _IntegerInterval(endpoints=(1, 'n'), inclusive=(True, True))
    _r_param = _RealParameter('r', domain=_r_domain, typical=(1, 2))

    _n_domain = _IntegerInterval(endpoints=(1, inf), inclusive=(True, True))
    _n_param = _RealParameter('n', domain=_n_domain, typical=(1, 4))

    _r_domain.define_parameters(_n_param)

    _parameterizations = [_Parameterization(_r_param, _n_param)]

    def __init__(self, dist, /, *args, r, n, **kwargs):
        super().__init__(dist, *args, r=r, n=n, **kwargs)

    def _support(self, *args, r, n, **kwargs):
        return self._dist._support(*args, **kwargs)

    def _process_parameters(self, r=None, n=None, **params):
        parameters = self._dist._process_parameters(**params)
        parameters.update(dict(r=r, n=n))
        return parameters

    def _overrides(self, method_name):
        return method_name in {'_logpdf_formula', '_pdf_formula',
                               '_cdf_formula', '_ccdf_formula',
                               '_icdf_formula', '_iccdf_formula'}

    def _logpdf_formula(self, x, r, n, **kwargs):
        xp = self._xp
        log_factor = special.betaln(r, n - r + 1)
        log_fX = self._dist._logpdf_dispatch(x, **kwargs)
        # log-methods sometimes use complex dtype with 0 imaginary component,
        # but `_tanhsinh` doesn't accept complex limits of integration; take `real`.
        log_FX = self._dist._logcdf_dispatch(xp.real(x), **kwargs)
        log_cFX = self._dist._logccdf_dispatch(xp.real(x), **kwargs)
        # This can be problematic when (r - 1)|(n-r) = 0 and `log_FX`|log_cFX = -inf
        # The PDF in these cases is 0^0, so these should be replaced with log(1)=0
        # return log_fX + (r-1)*log_FX + (n-r)*log_cFX - log_factor
        rm1_log_FX = xp.where((r - 1 == 0) & xp.isneginf(log_FX), 0, (r-1)*log_FX)
        nmr_log_cFX = xp.where((n - r == 0) & xp.isneginf(log_cFX), 0, (n-r)*log_cFX)
        return log_fX + rm1_log_FX + nmr_log_cFX - log_factor

    def _pdf_formula(self, x, r, n, **kwargs):
        # 1 / factor = factorial(n) / (factorial(r-1) * factorial(n-r))
        factor = special.beta(r, n - r + 1)
        fX = self._dist._pdf_dispatch(x, **kwargs)
        FX = self._dist._cdf_dispatch(x, **kwargs)
        cFX = self._dist._ccdf_dispatch(x, **kwargs)
        return fX * FX**(r-1) * cFX**(n-r) / factor

    def _cdf_formula(self, x, r, n, **kwargs):
        x_ = self._dist._cdf_dispatch(x, **kwargs)
        return special.betainc(r, n-r+1, x_)

    def _ccdf_formula(self, x, r, n, **kwargs):
        x_ = self._dist._cdf_dispatch(x, **kwargs)
        return special.betaincc(r, n-r+1, x_)

    def _icdf_formula(self, p, r, n, **kwargs):
        p_ = special.betaincinv(r, n-r+1, p)
        return self._dist._icdf_dispatch(p_, **kwargs)

    def _iccdf_formula(self, p, r, n, **kwargs):
        p_ = special.betainccinv(r, n-r+1, p)
        return self._dist._icdf_dispatch(p_, **kwargs)

    def __repr__(self):
        with np.printoptions(threshold=10):
            return (f"order_statistic({repr(self._dist)}, r={repr(self.r)}, "
                    f"n={repr(self.n)})")

    def __str__(self):
        with np.printoptions(threshold=10):
            return (f"order_statistic({str(self._dist)}, r={str(self.r)}, "
                    f"n={str(self.n)})")


@rv_capabilities
def order_statistic(X, /, *, r, n):
    r"""Probability distribution of an order statistic.

    Returns a random variable that follows the distribution of the
    :math:`r^{\text{th}}` order statistic of a sample of :math:`n`
    observations of a random variable :math:`X`.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable :math:`X`.
    r : array_like
        The (positive integer) rank of the order statistic :math:`r`,
        satisfying ``1 <= r <= n``.
    n : array_like
        The (positive integer) sample size :math:`n`.

    Returns
    -------
    Y : `ContinuousDistribution`
        A random variable that follows the distribution of the prescribed
        order statistic.

    Notes
    -----
    If we make :math:`n` observations of a continuous random variable
    :math:`X` and sort them in increasing order
    :math:`X_{(1)}, \dots, X_{(r)}, \dots, X_{(n)}`, then :math:`X_{(r)}`
    is known as the :math:`r^{\text{th}}` order statistic.

    If the PDF, CDF, and CCDF of :math:`X` are denoted by :math:`f`,
    :math:`F`, and :math:`G = 1 - F`, respectively, then the PDF of
    :math:`X_{(r)}` is given by:

    .. math::

        f_r(x) = \frac{n!}{(r-1)! (n-r)!} f(x) F(x)^{r-1} G(x)^{n - r}

    The CDF and other methods of the distribution of :math:`X_{(r)}`
    are calculated using the fact that :math:`X = F^{-1}(U)`, where :math:`U` is a
    standard uniform random variable, together with the fact that the order statistics
    of i.i.d. uniform random variables follow a beta distribution
    :math:`B(r, n - r + 1)`.

    References
    ----------
    .. [1] Order statistic. *Wikipedia*. https://en.wikipedia.org/wiki/Order_statistic

    Examples
    --------
    Suppose we are interested in order statistics of samples of size five drawn
    from the standard normal distribution. Plot the PDF of each
    order statistic and compare with a normalized histogram from simulation.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>>
    >>> X = stats.Normal()
    >>> data = X.sample(shape=(10000, 5))
    >>> sorted_data = np.sort(data, axis=1)
    >>> Y = stats.order_statistic(X, r=[1, 2, 3, 4, 5], n=5)
    >>>
    >>> ax = plt.gca()
    >>> colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    >>> for i in range(5):
    ...     y = sorted_data[:, i]
    ...     ax.hist(y, density=True, bins=30, alpha=0.1, color=colors[i])
    >>> Y.plot(ax=ax)
    >>> plt.show()

    """
    xp = X._xp
    r, n = xp.asarray(r), xp.asarray(n)
    if xp.any((r != xp.floor(r)) | (r < 0)) or xp.any((n != xp.floor(n)) | (n < 0)):
        message = "`r` and `n` must contain only positive integers."
        raise ValueError(message)
    return OrderStatisticDistribution(X, r=r, n=n, xp=xp)


def _raise_if_not_continuous(f):
    @functools.wraps(f)
    def wrapped(self, *args, **kwargs):
        if not self._continuous:
            raise NotImplementedError(f"`{f.__name__}` is implemented only for "
                                      "purely continuous `Mixture`s.")
        return f(self, *args, **kwargs)
    return wrapped


@rv_capabilities
class Mixture(_ProbabilityDistribution):
    r"""Representation of a mixture distribution.

    A mixture distribution is the distribution of a random variable
    defined in the following way: first, a random variable is selected
    from `components` according to the probabilities given by `weights`, then
    the selected random variable is realized.

    Parameters
    ----------
    components : sequence of `UnivariateDistribution`
        The underlying instances of `UnivariateDistribution`.
        All must have scalar shape parameters (if any); e.g., the `mean`
        must return a scalar.
    weights : sequence of floats, optional
        The corresponding probabilities of selecting each random variable.
        Must be non-negative and sum to one. The default behavior is to weight
        all components equally.

    Attributes
    ----------
    components : sequence of `UnivariateDistribution`
        The underlying instances of `UnivariateDistribution`.
    weights : ndarray
        The corresponding probabilities of selecting each random variable.

    Methods
    -------
    support

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

    Notes
    -----
    Methods `median`, `mode`, `entropy`, `logentropy`, `icdf`, `iccdf`, `ilogcdf`,
    and `ilogccdf` are currently implemented only when all components are continuous
    distributions.

    The following abbreviations are used throughout the documentation.

    - PDF: probability density function
    - CDF: cumulative distribution function
    - CCDF: complementary CDF
    - entropy: differential entropy
    - log-*F*: logarithm of *F* (e.g. log-CDF)
    - inverse *F*: inverse function of *F* (e.g. inverse CDF)

    References
    ----------
    .. [1] Mixture distribution, *Wikipedia*,
           https://en.wikipedia.org/wiki/Mixture_distribution
    .. [2] Zero-inflated model, *Wikipedia*,
           https://en.wikipedia.org/wiki/Zero-inflated_model

    Examples
    --------
    A mixture of normal distributions:

    >>> import numpy as np
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt
    >>> X1 = stats.Normal(mu=-2, sigma=1)
    >>> X2 = stats.Normal(mu=2, sigma=1)
    >>> mixture = stats.Mixture([X1, X2], weights=[0.4, 0.6])
    >>> print(f'mean: {mixture.mean():.2f}, '
    ...       f'median: {mixture.median():.2f}, '
    ...       f'mode: {mixture.mode():.2f}')
    mean: 0.40, median: 1.04, mode: 2.00
    >>> x = np.linspace(-10, 10, 300)
    >>> plt.plot(x, mixture.pdf(x))
    >>> plt.title('PDF of normal distribution mixture')
    >>> plt.show()

    A zero-inflated Poisson model.

    >>> Poisson = stats.make_distribution(stats.poisson)
    >>> lmb = 10  # expected Poisson count
    >>> pi = 0.4  # probability of extra zeros
    >>> X = stats.Binomial(n=0, p=0)  # degenerate distribution; all zero
    >>> Y = Poisson(mu=lmb)
    >>> ZIP = stats.Mixture((X, Y), weights=(pi, 1-pi))
    >>> true_mean = (1 - pi) * lmb  # see e.g. [2]
    >>> print(np.allclose(ZIP.mean(), true_mean))
    True
    >>> true_variance = lmb*(1 - pi)*(1 + pi*lmb)  # see e.g. [2]
    >>> print(np.allclose(ZIP.variance(), true_variance))
    True

    """
    # Todo:
    # Add support for array shapes, weights

    def _input_validation(self, components, weights):
        if len(components) == 0:
            message = ("`components` must contain at least one random variable.")
            raise ValueError(message)

        continuous = True
        for var in components:
            if not isinstance(var, UnivariateDistribution):
                message = ("Each element of `components` must be an instance of "
                           "`UnivariateDistribution`.")
                raise ValueError(message)
            if not var._shape == ():
                message = "All elements of `components` must have scalar shapes."
                raise ValueError(message)
            if isinstance(var, CircularDistribution):
                message = "`Mixture` does not currently support circular components."
                raise NotImplementedError(message)
            continuous = continuous and isinstance(var, ContinuousDistribution)

        xp = components[0]._xp
        self._xp = xp
        self._device = components[0]._device

        if weights is None:
            return components, weights, continuous

        weights = xp.asarray(weights)
        if weights.shape != (len(components),):
            message = "`components` and `weights` must have the same length."
            raise ValueError(message)

        if not xp.isdtype(weights.dtype, "real floating"):
            message = "`weights` must have floating point dtype."
            raise ValueError(message)

        if not xpx.isclose(xp.sum(weights), 1.0, xp=xp):
            message = "`weights` must sum to 1.0."
            raise ValueError(message)

        if not xp.all(weights >= 0):
            message = "All `weights` must be non-negative."
            raise ValueError(message)

        return components, weights, continuous

    def __init__(self, components, *, weights=None):
        components, weights, continuous = self._input_validation(components, weights)
        xp = self._xp
        n = len(components)
        dtype = xp.result_type(*(var._dtype for var in components))
        self._shape = xp.broadcast_shapes(*(var._shape for var in components))
        self._dtype, self._components = dtype, components
        self._weights = (xp.full(n, 1/n, dtype=dtype, device=self._device)
                         if weights is None else weights)
        self._continuous = continuous
        self.validation_policy = None

    @property
    def components(self):
        return list(self._components)

    @property
    def weights(self):
        return xp_copy(self._weights)

    def _full(self, val, *args):
        xp = self._xp
        args = [xp.asarray(arg, device=self._device) for arg in args]
        dtype = xp.result_type(self._dtype, *(arg.dtype for arg in args))
        shape = xp.broadcast_shapes(self._shape, *(arg.shape for arg in args))
        return xp.full(shape, val, dtype=dtype, device=self._device)

    def _sum(self, fun, *args):
        out = self._full(0, *args)
        for var, weight in zip(self._components, self._weights):
            out += getattr(var, fun)(*args) * weight
        return out[()]

    def _logsum(self, fun, *args):
        xp = self._xp
        out = self._full(-xp.inf, *args)
        for var, log_weight in zip(self._components, xp.log(self._weights)):
            xp.logaddexp(out, getattr(var, fun)(*args) + log_weight, out=out)
        return out[()]

    def support(self):
        xp = self._xp
        a = self._full(xp.inf)
        b = self._full(-xp.inf)
        for var in self._components:
            a = xp.minimum(a, var.support()[0])
            b = xp.maximum(b, var.support()[1])
        return a, b

    def _raise_if_method(self, method):
        if method is not None:
            raise NotImplementedError("`method` not implemented for this distribution.")

    @_raise_if_not_continuous
    def logentropy(self, *, method=None):
        xp = self._xp
        self._raise_if_method(method)
        def log_integrand(x):
            # `x` passed by `_tanhsinh` will be of complex dtype because
            # `log_integrand` returns complex values, but the imaginary
            # component is always zero. Extract the real part because
            # `logpdf` uses `logaddexp`, which fails for complex input.
            return self.logpdf(xp.real(x)) + xp.log(self.logpdf(xp.real(x)) + 0j)

        res = _tanhsinh(log_integrand, *self.support(), log=True).integral
        return _log_real_standardize(res + xp.pi*1j, xp=xp)

    @_raise_if_not_continuous
    def entropy(self, *, method=None):
        self._raise_if_method(method)
        return _tanhsinh(lambda x: -self.pdf(x) * self.logpdf(x),
                         *self.support()).integral

    @_raise_if_not_continuous
    def mode(self, *, method=None):
        self._raise_if_method(method)
        a, b = self.support()
        def f(x): return -self.pdf(x)
        res = elementwise.bracket_minimum(f, 1., xmin=a, xmax=b)
        res = elementwise.find_minimum(f, res.bracket)
        return res.x

    @_raise_if_not_continuous
    def median(self, *, method=None):
        self._raise_if_method(method)
        return self.icdf(0.5)

    def mean(self, *, method=None):
        self._raise_if_method(method)
        return self._sum('mean')

    def variance(self, *, method=None):
        self._raise_if_method(method)
        return self._moment_central(2)

    def standard_deviation(self, *, method=None):
        self._raise_if_method(method)
        return self.variance()**0.5

    def skewness(self, *, method=None):
        self._raise_if_method(method)
        return self._moment_standardized(3)

    def kurtosis(self, *, method=None):
        self._raise_if_method(method)
        return self._moment_standardized(4)

    def moment(self, order=1, kind='raw', *, method=None):
        self._raise_if_method(method)
        kinds = {'raw': self._moment_raw,
                 'central': self._moment_central,
                 'standardized': self._moment_standardized}
        order = ContinuousDistribution._validate_order(self, order)
        kind = ContinuousDistribution._validate_kind(self, kind, kinds)
        moment_kind = kinds[kind]
        return moment_kind(order)

    def _moment_raw(self, order):
        out = self._full(0)
        for var, weight in zip(self._components, self._weights):
            out += var.moment(order, kind='raw') * weight
        return out[()]

    def _moment_central(self, order):
        order = int(order)
        out = self._full(0)
        for var, weight in zip(self._components, self._weights):
            moment_as = [var.moment(order, kind='central')
                         for order in range(order + 1)]
            a, b = var.mean(), self.mean()
            moment = var._moment_transform_center(order, moment_as, a, b)
            out += moment * weight
        return out[()]

    def lmoment(self, order=1, *, standardize=True, method=None):
        message = "L-moments are not currently available for mixture distributions."
        raise NotImplementedError(message)

    def _moment_standardized(self, order):
        return self._moment_central(order) / self.standard_deviation()**order

    def pdf(self, x, /, *, method=None):
        self._raise_if_method(method)
        return self._sum('pdf', x)

    def logpdf(self, x, /, *, method=None):
        self._raise_if_method(method)
        return self._logsum('logpdf', x)

    def pmf(self, x, /, *, method=None):
        self._raise_if_method(method)
        return self._sum('pmf', x)

    def logpmf(self, x, /, *, method=None):
        self._raise_if_method(method)
        return self._logsum('logpmf', x)

    def cdf(self, x, y=None, /, *, method=None):
        self._raise_if_method(method)
        args = (x,) if y is None else (x, y)
        return self._sum('cdf', *args)

    def logcdf(self, x, y=None, /, *, method=None):
        self._raise_if_method(method)
        args = (x,) if y is None else (x, y)
        return self._logsum('logcdf', *args)

    def ccdf(self, x, y=None, /, *, method=None):
        self._raise_if_method(method)
        args = (x,) if y is None else (x, y)
        return self._sum('ccdf', *args)

    def logccdf(self, x, y=None, /, *, method=None):
        self._raise_if_method(method)
        args = (x,) if y is None else (x, y)
        return self._logsum('logccdf', *args)

    def _invert(self, fun, p):
        xp = self._xp
        xmin, xmax = self.support()
        fun = getattr(self, fun)
        f = lambda x, p: fun(x) - p  # noqa: E731 is silly
        xl0, xr0 = _guess_bracket(xmin, xmax, xp=xp)
        xl0, xr0, xmin, xmax, p = xp_promote(xl0, xr0, xmin, xmax, p, xp=xp)
        res = elementwise.bracket_root(f, xl0=xl0, xr0=xr0,
                                       xmin=xmin, xmax=xmax, args=(p,))
        return elementwise.find_root(f, res.bracket, args=(p,)).x

    @_raise_if_not_continuous
    def icdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('cdf', p)

    @_raise_if_not_continuous
    def iccdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('ccdf', p)

    @_raise_if_not_continuous
    def ilogcdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('logcdf', p)

    @_raise_if_not_continuous
    def ilogccdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('logccdf', p)

    def sample(self, shape=(), *, rng=None, method=None):
        xp = self._xp
        self._raise_if_method(method)
        rng = xp.random.default_rng(rng) if is_cupy(xp) else np.random.default_rng(rng)
        size = xp.prod(xpx.atleast_nd(xp.asarray(shape), ndim=1, xp=xp))
        ns = rng.multinomial(size, self._weights)
        x = xp.concat([var.sample(shape=n, rng=rng)
                       for n, var in zip(ns, self._components)])
        i = xp.asarray(rng.permutation(x.shape[0]), device=xp_device(x))
        shape = (shape,) if np.isscalar(shape) else shape
        x = xp.reshape(x[i], shape)
        return x[()]

    def __repr__(self):
        result = "Mixture(\n"
        result += "    [\n"
        with np.printoptions(threshold=10):
            for component in self.components:
                result += f"        {repr(component)},\n"
            result += "    ],\n"
            result += f"    weights={repr(self.weights)},\n"
        result += ")"
        return result

    def __str__(self):
        result = "Mixture(\n"
        result += "    [\n"
        with np.printoptions(threshold=10):
            for component in self.components:
                result += f"        {str(component)},\n"
            result += "    ],\n"
            result += f"    weights={str(self.weights)},\n"
        result += ")"
        return result

@rv_capabilities
def abs(X, /):
    r"""Absolute value of a random variable.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable :math:`X`.

    Returns
    -------
    Y : `ContinuousDistribution`
        A random variable :math:`Y = |X|`.

    Examples
    --------
    Suppose we have a normally distributed random variable :math:`X`:

    >>> import numpy as np
    >>> from scipy import stats
    >>> X = stats.Normal()

    We wish to have a random variable :math:`Y` distributed according to
    the folded normal distribution; that is, a random variable :math:`|X|`.

    >>> Y = stats.abs(X)

    The PDF of the distribution in the left half plane is "folded" over to
    the right half plane. Because the normal PDF is symmetric, the resulting
    PDF is zero for negative arguments and doubled for positive arguments.

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(0, 5, 300)
    >>> ax = plt.gca()
    >>> Y.plot(x='x', y='pdf', t=('x', -1, 5), ax=ax)
    >>> plt.plot(x, 2 * X.pdf(x), '--')
    >>> plt.legend(('PDF of `Y`', 'Doubled PDF of `X`'))
    >>> plt.show()

    """
    return FoldedDistribution(X)


@rv_capabilities
def exp(X, /):
    r"""Natural exponential of a random variable.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable :math:`X`.

    Returns
    -------
    Y : `ContinuousDistribution`
        A random variable :math:`Y = \exp(X)`.

    Examples
    --------
    Suppose we have a normally distributed random variable :math:`X`:

    >>> import numpy as np
    >>> from scipy import stats
    >>> X = stats.Normal()

    We wish to have a lognormally distributed random variable :math:`Y`,
    a random variable whose natural logarithm is :math:`X`.
    If :math:`X` is to be the natural logarithm of :math:`Y`, then we
    must take :math:`Y` to be the natural exponential of :math:`X`.

    >>> Y = stats.exp(X)

    To demonstrate that ``X`` represents the logarithm of ``Y``,
    we plot a normalized histogram of the logarithm of observations of
    ``Y`` against the PDF underlying ``X``.

    >>> import matplotlib.pyplot as plt
    >>> rng = np.random.default_rng(435383595582522)
    >>> y = Y.sample(shape=10000, rng=rng)
    >>> ax = plt.gca()
    >>> ax.hist(np.log(y), bins=50, density=True)
    >>> X.plot(ax=ax)
    >>> plt.legend(('PDF of `X`', 'histogram of `log(y)`'))
    >>> plt.show()

    """
    xp = X._xp
    return MonotonicTransformedDistribution(X, g=xp.exp, h=xp.log, dh=lambda u: 1 / u,
                                            logdh=lambda u: -xp.log(u), xp=xp)


@rv_capabilities
def log(X, /):
    r"""Natural logarithm of a non-negative random variable.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable :math:`X` with positive support.

    Returns
    -------
    Y : `ContinuousDistribution`
        A random variable :math:`Y = \log(X)`.

    Examples
    --------
    Suppose we have a gamma distributed random variable :math:`X`:

    >>> import numpy as np
    >>> from scipy import stats
    >>> Gamma = stats.make_distribution(stats.gamma)
    >>> X = Gamma(a=1.0)

    We wish to have an exp-gamma distributed random variable :math:`Y`,
    a random variable whose natural exponential is :math:`X`.
    If :math:`X` is to be the natural exponential of :math:`Y`, then we
    must take :math:`Y` to be the natural logarithm of :math:`X`.

    >>> Y = stats.log(X)

    To demonstrate that ``X`` represents the exponential of ``Y``,
    we plot a normalized histogram of the exponential of observations of
    ``Y`` against the PDF underlying ``X``.

    >>> import matplotlib.pyplot as plt
    >>> rng = np.random.default_rng(435383595582522)
    >>> y = Y.sample(shape=10000, rng=rng)
    >>> ax = plt.gca()
    >>> ax.hist(np.exp(y), bins=50, density=True)
    >>> X.plot(ax=ax)
    >>> plt.legend(('PDF of `X`', 'histogram of `exp(y)`'))
    >>> plt.show()

    """
    xp = X._xp
    if xp.any(X.support()[0] < 0):
        message = ("The logarithm of a random variable is only implemented when the "
                   "support is non-negative.")
        raise NotImplementedError(message)
    return MonotonicTransformedDistribution(X, g=xp.log, h=xp.exp, dh=xp.exp,
                                            logdh=lambda u: u, xp=xp)
