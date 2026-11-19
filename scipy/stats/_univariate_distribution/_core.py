"""Core univariate distribution classes.

Currently also contains some private helpers (e.g. _set_invalid_nan) and
somewhat specific distribution classes (e.g. MonotonicTransformedDistribution)
to avoid circular imports; ideally these move to `_util.py` and `_extensions.py`
over time.
"""
import functools
from functools import cached_property
import math
from math import inf

import numpy as np
from scipy._lib._array_api import (
    array_namespace, xp_promote, xp_size, xp_copy, xp_spacing,
    is_torch, xp_device, xp_result_device, xp_result_type,
    xp_compat_namespace
)
from scipy._lib._docscrape import ClassDoc, NumpyDocString
from scipy._lib._util import _rng_spawn
from scipy._external import array_api_extra as xpx
from scipy import special
from scipy.differentiate import derivative
from scipy.integrate import tanhsinh as _tanhsinh, nsum
from scipy.optimize import elementwise
from scipy.stats._univariate_distribution._docs import _ProbabilityDistribution
from scipy.stats._stats_py import _log1mexp
from scipy.stats import qmc
from scipy.stats._univariate_distribution._util import (
    _RealInterval, _RealParameter, _Parameterization,
    _set_invalid_nan_property, _dispatch, _cdf2_input_validation, _SKIP_ALL, _NO_CACHE,
    _guess_bracket, _log_real_standardize, _generate_domain_support, _logexpxmexpy,
    rv_capabilities, _shift_scale_distribution_function_2arg,
    _shift_scale_distribution_function, _shift_scale_inverse_function
)


# in case we need to distinguish between None and not specified
# Typically this is used to determine whether the tolerance has been set by the
# user and make a decision about which method to use to evaluate a distribution
# function. Sometimes, the logic does not consider the value of the tolerance,
# only whether this has been defined or not. This is not intended to be the
# best possible logic; the intent is to establish the structure, which can
# be refined in follow-up work.
# See https://github.com/scipy/scipy/pull/21050#discussion_r1714195433.
_null = object()
def _isnull(x):
    return type(x) is object or x is None


# TODO:
#  Test dtype behavior - specified dtype is only the base dtype for promotion with
#    parameters and arguments. Test sample dtypes.
#  Add MPArray support
#  Add distribution-specific checks for finite vs undefined vs infinite moments.
#  `_tanhsinh` special case when there are no abscissae between the limits
#     example: cdf of uniform between 1.0 and np.nextafter(1.0, np.inf)
#  implement symmetric distribution (take advantage of symmetry)
#  implement composite distribution (for using different implementations depending on
#    shape parameters - e.g. t distribution with df = np.inf delegates to normal)
#  Investigate use `median` information throughout, e.g. to improve integration of
#    CDF or to provide initial bracket for ICDF. (Only if `_median_formula` is
#    provided.)
#  Document user tips for faster execution:
#  - pass NumPy arrays
#  - pass inputs of floating point type (not integers)
#  - prefer NumPy scalars or 0d arrays over other size 1 arrays
#  - pass no invalid parameters and disable invalid parameter checks (validation_policy)
#  - provide a Generator if you're going to do sampling
#  Add options for drawing parameters: log-spacing
#  Accuracy benchmark suite. Ideally, when all code is array API compatible and mparray
#    provides an array-API compatible arbitrary precision array, we can just compare
#    results when using NumPy backend against mparray backend.
#  Should caches be attributes so we can more easily ensure that they are not
#    modified when caching is turned off?
#  Make ShiftedScaledDistribution more efficient - only process underlying
#    distribution parameters as necessary.
#  Reconsider `all_inclusive` - see comment in `contains` method of `_interval`
# `validation_policy` needs testing
# `tol` does not offer very fine-grained control; consider improving
# report accuracy estimates? How?
# reconsider distribution mutability

# Originally, I planned to filter out invalid distribution parameters;
# distribution implementation functions would always work with "compressed",
# 1D arrays containing only valid distribution parameters. There are two
# problems with this:
# - This essentially requires copying all arrays, even if there is only a
#   single invalid parameter combination. This is expensive. Then, to output
#   the original size data to the user, we need to "decompress" the arrays
#   and fill in the NaNs, so more copying. Unless we branch the code when
#   there are no invalid data, these copies happen even in the normal case,
#   where there are no invalid parameter combinations. We should not incur
#   all this overhead in the normal case.
# - For methods that accept arguments other than distribution parameters, the
#   user will pass in arrays that are broadcastable with the original arrays,
#   not the compressed arrays. This means that this same sort of invalid
#   value detection needs to be repeated every time one of these methods is
#   called.
# The much simpler solution is to keep the data uncompressed but to replace
# the invalid parameters and arguments with NaNs (and only if some are
# invalid). With this approach, the copying happens only if/when it is
# needed. Most functions involved in stats distribution calculations don't
# mind NaNs; they just return NaN. The behavior "If x_i is NaN, the result
# is NaN" is explicit in the array API. So this should be fine.
#
# Currently, I am still leaving the parameters and function arguments
# in their broadcasted shapes rather than, say, raveling. The intent
# is to avoid back and forth reshaping. If authors of distributions have
# trouble dealing with N-D arrays, we can reconsider this.
#
# Another important decision is that the *private* methods must accept
# the distribution parameters as inputs rather than relying on these
# cached properties directly (although the public methods typically pass
# the cached values to the private methods). This is because the elementwise
# algorithms for quadrature, differentiation, root-finding, and minimization
# prefer that the input functions are strictly elementwise in the sense
# that the value output for a given input element does not depend on the
# shape of the input or that element's location within the input array.
# When the computation has converged for an element, it is removed from
# the computation entirely. As a result, the shape of the arrays passed to
# the function will almost never be broadcastable with the shape of the
# cached parameter arrays.
#
# I've sprinkled in some optimizations for scalars and same-shape/type arrays
# throughout. The biggest time sinks before were:
# - broadcast_arrays
# - result_dtype
# - is_subdtype
# It is much faster to check whether these are necessary than to do them.


# ideally we would move this back to _util.py but currently it needs
# `CircularDistribution` and `DiscreteDistribution`.
def _set_invalid_nan(f):
    # Wrapper for input / output validation and standardization of distribution
    # functions that accept either the quantile or percentile as an argument:
    # logpdf, pdf
    # logpmf, pmf
    # logcdf, cdf
    # logccdf, ccdf
    # ilogcdf, icdf
    # ilogccdf, iccdf
    # Arguments that are outside the required range are replaced by NaN before
    # passing them into the underlying function. The corresponding outputs
    # are replaced by the appropriate value before being returned to the user.
    # For example, when the argument of `cdf` exceeds the right end of the
    # distribution's support, the wrapper replaces the argument with NaN,
    # ignores the output of the underlying function, and returns 1.0. It also
    # ensures that output is of the appropriate shape and dtype.

    endpoints = {'icdf': (0, 1), 'iccdf': (0, 1),
                 'ilogcdf': (-inf, 0), 'ilogccdf': (-inf, 0)}
    replacements = {'logpdf': (-inf, -inf), 'pdf': (0, 0),
                    'logpmf': (-inf, -inf), 'pmf': (0, 0),
                    '_logcdf1': (-inf, 0), '_logccdf1': (0, -inf),
                    '_cdf1': (0, 1), '_ccdf1': (1, 0)}
    replace_strict = {'pdf', 'logpdf', 'pmf', 'logpmf'}
    replace_exact = {'icdf', 'iccdf', 'ilogcdf', 'ilogccdf'}
    clip = {'_cdf1', '_ccdf1'}
    clip_log = {'_logcdf1', '_logccdf1'}
    # relevant to discrete distributions only
    replace_non_integral = {'pmf', 'logpmf', 'pdf', 'logpdf'}
    cdflike = {'_cdf1', '_logcdf1', '_ccdf1', '_logccdf1'}
    # relevant to circular distributions only
    wrap_unit = {'icdf', 'iccdf'}
    no_unwrap = {'pdf', 'logpdf', 'pmf', 'logpmf'}

    @functools.wraps(f)
    def filtered(self, x, *args, **kwargs):
        xp = self._xp

        if self.validation_policy == _SKIP_ALL:
            return f(self, x, *args, **kwargs)

        method_name = f.__name__
        device = self._device
        dtype = self._dtype
        shape = self._shape
        circular = isinstance(self, CircularDistribution)
        discrete = isinstance(self, DiscreteDistribution)
        discrete_cdflike = discrete and method_name in cdflike

        # Ensure that argument is at least as precise as distribution
        # parameters, which are already at least floats. This will avoid issues
        # with raising integers to negative integer powers and failure to replace
        # invalid integers with NaNs.
        if hasattr(x, 'dtype') and x.dtype != dtype:
            dtype = xp.result_type(x.dtype, dtype)

        x = xp.asarray(x, dtype=dtype, device=device)

        # Broadcasting is slow. Do it only if necessary.
        if not x.shape == shape:
            try:
                shape = xp.broadcast_shapes(x.shape, shape)
                x = xp.broadcast_to(x, shape)
                # Should we broadcast the distribution parameters to this shape, too?
            except (ValueError, RuntimeError) as e:
                message = (
                    f"The argument provided to `{self.__class__.__name__}"
                    f".{method_name}` cannot be broadcast to the same "
                    "shape as the distribution parameters.")
                raise ValueError(message) from e

        low, high = endpoints.get(method_name, self.support())

        if circular:
            x = xp.asarray(x, dtype=dtype, copy=True)  # ideally, avoid multiple copies
            x[xp.isinf(x)] = xp.nan
            a, b = self.support()  # not the same as low, high for inverse methods
            period = b - a
            if method_name in wrap_unit:
                turn = x // 1
                x = x % 1
            else:
                turn = (x - a) // period
                x = (x - a) % period + a
            # much of the code below that deals with high/low values could be skipped
            # leaving it untouched while adding circular distributions; optimize later

        # Check for arguments outside of domain. They'll be replaced with NaNs,
        # and the result will be set to the appropriate value.
        left_inc, right_inc = self._variable.domain.inclusive
        mask_low = (x < low if (method_name in replace_strict and left_inc)
                    or discrete_cdflike else x <= low)
        mask_high = (x > high if (method_name in replace_strict and right_inc)
                     else x >= high)
        mask_invalid = (mask_low | mask_high)
        any_invalid = (mask_invalid if mask_invalid.shape == ()
                       else xp.any(mask_invalid))

        # Check for arguments at domain endpoints, whether they
        # are part of the domain or not.
        any_endpoint = False
        if method_name in replace_exact:
            mask_low_endpoint = (x == low)
            mask_high_endpoint = (x == high)
            mask_endpoint = (mask_low_endpoint | mask_high_endpoint)
            any_endpoint = (mask_endpoint if mask_endpoint.shape == ()
                            else xp.any(mask_endpoint))

        # Check for non-integral arguments to CDF-like method of discrete distribution
        if discrete_cdflike:
            x = xp.floor(x)

        # Check for non-integral arguments to PMF method
        # or PDF of a discrete distribution.
        any_non_integral = False
        if discrete and method_name in replace_non_integral:
            mask_non_integral = (x != xp.floor(x))
            any_non_integral = (mask_non_integral if mask_non_integral.shape == ()
                                else xp.any(mask_non_integral))

        # Set out-of-domain arguments to NaN. The result will be set to the
        # appropriate value later.
        if any_invalid:
            x = xp.asarray(x, dtype=dtype, copy=True)
            x[mask_invalid] = xp.nan

        res = xp.asarray(f(self, x, *args, **kwargs))

        # Ensure that the result is the correct dtype and shape,
        # copying (only once) if necessary.
        res_needs_copy = False
        if res.dtype != dtype:
            dtype = xp.result_type(dtype, self._dtype)
            res_needs_copy = True

        if res.shape != shape:  # faster to check first
            res = xp.broadcast_to(res, self._shape)
            res_needs_copy = (res_needs_copy or any_invalid
                              or any_endpoint or any_non_integral)

        if res_needs_copy:
            res = xp.asarray(res, dtype=dtype, copy=True)

        # For non-integral arguments to PMF (and PDF of discrete distribution)
        # replace with zero.
        if any_non_integral:
            zero = -xp.inf if method_name in {'logpmf', 'logpdf'} else 0
            res[mask_non_integral & ~xp.isnan(res)] = zero

        # For arguments outside the function domain, replace results
        if any_invalid:
            replace_low, replace_high = (
                replacements.get(method_name, (math.nan, math.nan)))
            res[mask_low] = replace_low
            res[mask_high] = replace_high

        # For arguments at the endpoints of the domain, replace results
        if any_endpoint:
            a, b = self.support()
            if a.shape != shape:
                a = xp.asarray(xp.broadcast_to(a, shape), copy=True)
                b = xp.asarray(xp.broadcast_to(b, shape), copy=True)

            replace_low_endpoint = (
                b[mask_low_endpoint] if method_name.endswith('ccdf')
                else a[mask_low_endpoint])
            replace_high_endpoint = (
                a[mask_high_endpoint] if method_name.endswith('ccdf')
                else b[mask_high_endpoint])

            if not discrete_cdflike:
                res[mask_low_endpoint] = replace_low_endpoint
            res[mask_high_endpoint] = replace_high_endpoint

        # Clip probabilities to [0, 1]
        if not circular and method_name in clip:
            res = xp.clip(res, 0., 1.)
        elif not circular and method_name in clip_log:
            res = xp.real(res)  # exp(res) > 0
            res = xp.clip(res, None, 0.)  # exp(res) < 1

        if circular and method_name not in no_unwrap:
            turn = -turn if method_name in {'_ccdf1', 'iccdf'} else turn
            res += turn * period if method_name in wrap_unit else turn

        return res[()]

    return filtered


class UnivariateDistribution(_ProbabilityDistribution):
    r""" Class that represents a continuous statistical distribution.

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
    xp : module, optional
        The backend to use for calculations - either NumPy or CuPy. Default is
        inferred from the parameters if possible or NumPy otherwise.
    dtype : dtype, optional
        The base dtype of the distribution. The output dtype of methods that
        accept arguments will follow standard promotion rules. Default is inferred
        from the parameters if possible and the default floating point dtype of
        the backend otherwise.
    device : device, optional
        The base device of the distribution. Must be compatible with method arguments.
        Default is inferred from the parameters if possible and the default device of
        the backend otherwise.

    Attributes
    ----------
    Parameters `tol`, `validation_policy`, and `cache_policy` are available
    as attributes.

    Methods
    -------
    support

    plot

    sample

    moment
    lmoment

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

    See Also
    --------
    :ref:`rv_infrastructure` : Tutorial

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
    __array_priority__ = 1
    _parameterizations = []

    ### Initialization

    def __init__(self, *, tol=_null, validation_policy=None, cache_policy=None,
                 xp=None, dtype=None, device=None, **parameters):
        self.validation_policy = validation_policy
        self.cache_policy = cache_policy
        self._not_implemented = (
            f"`{self.__class__.__name__}` does not provide an accurate "
            "implementation of the required method. Consider leaving "
            "`method` and `tol` unspecified to use another implementation."
        )
        self._original_parameters = {}
        # We may want to override the `__init__` method with parameters so
        # IDEs can suggest parameter names. If there are multiple parameterizations,
        # we'll need the default values of parameters to be None; this will
        # filter out the parameters that were not actually specified by the user.
        parameters = {key: val for key, val in
                      sorted(parameters.items()) if val is not None}

        # There is a bit of a chicken and egg problem when placing the array_namespace
        # determination in the existing code - we need to know the namespace to validate
        # the parameters, but we want to validate the parameters before array_namespace
        # takes them and produces a less specific error message. I think the best
        # compromise without turning this into a project is to try-except: if success,
        # we have our `xp`; if not, we raise a more specific message later.
        try:
            if xp is None and getattr(self, "_xp", None) is None:
                self._xp = array_namespace(*parameters.values())
            else:
                self._xp = xp or self._xp
                self._xp = (xp_compat_namespace(np) if self._xp == np
                            else self._xp)
            device = device or xp_result_device(*parameters.values())
            dtype = dtype or xp_result_type(*parameters.values(),
                                            getattr(self, '_dtype', None),
                                            force_floating=True, xp=self._xp)
        except TypeError:
            self._xp = xp_compat_namespace(np)

        self.tol = tol
        self._dtype = dtype
        self._device = device
        self._update_parameters(**parameters)

    def _update_parameters(self, *, validation_policy=None, **params):
        r""" Update the numerical values of distribution parameters.

        Parameters
        ----------
        validation_policy : {None, "skip_all"}
            Specifies the level of input validation to perform. Left unspecified,
            input validation is performed to ensure appropriate behavior in edge
            case (e.g. parameters out of domain, argument outside of distribution
            support, etc.) and improve consistency of output dtype, shape, etc.
            Pass ``'skip_all'`` to avoid the computational overhead of these
            checks when rough edges are acceptable.

        **params : array_like
            Desired numerical values of the distribution parameters. Any or all
            of the parameters initially used to instantiate the distribution
            may be modified. Parameters used in alternative parameterizations
            are not accepted.
        """
        xp = self._xp
        parameters = original_parameters = self._original_parameters.copy()
        parameters.update(**params)
        parameterization = None
        self._invalid = False
        self._any_invalid = False
        self._shape = tuple()
        self._ndim = 0
        self._size = 1

        if not len(self._parameterizations):
            xp_info = self._xp.__array_namespace_info__()
            inferred_dtype = xp_info.default_dtypes()['real floating']
            # see data-apis/array-api-compat/issues/489
            inferred_device = (xp_device(xp.asarray(0)) if is_torch(xp)  # skip device check # noqa: E501
                               else xp_info.default_device())

            if parameters:
                message = (f"The `{self.__class__.__name__}` distribution "
                           "family does not accept parameters, but parameters "
                           f"`{set(parameters)}` were provided.")
                raise ValueError(message)

        elif (validation_policy or self.validation_policy) == _SKIP_ALL:
            parameters = self._process_parameters(**parameters)
            inferred_dtype = xp_result_type(*parameters.values(), force_floating=True,
                                            xp=xp)
            inferred_device = xp_result_device(*parameters.values())

        else:
            # This is default behavior, which re-runs all parameter validations
            # even when only a single parameter is modified. For many
            # distributions, the domain of a parameter doesn't depend on other
            # parameters, so parameters could safely be modified without
            # re-validating all other parameters. To handle these cases more
            # efficiently, we could allow the developer to override this
            # behavior.

            # Currently the user can only update the original parameterization.
            # Even though that parameterization is already known,
            # `_identify_parameterization` is called to produce a nice error
            # message if the user passes other values. To be a little more
            # efficient, we could detect whether the values passed are
            # consistent with the original parameterization rather than finding
            # it from scratch. However, we might want other parameterizations
            # to be accepted, which would require other changes, so I didn't
            # optimize this.

            parameterization = self._identify_parameterization(parameters)
            parameters, shape, size, ndim = self._broadcast(
                parameters, dtype=self._dtype, device=self._device)
            parameters, invalid, any_invalid, dtype, device = (
                self._validate(parameterization, parameters))
            parameters = self._process_parameters(**parameters)

            self._invalid = invalid
            self._any_invalid = any_invalid
            self._shape = shape
            self._size = size
            self._ndim = ndim
            inferred_dtype = dtype
            inferred_device = device

        self._dtype = self._dtype or inferred_dtype
        self._device = self._device or inferred_device

        if self._dtype is not None:
            if not xp.isdtype(self._dtype, 'real floating'):
                message = "`dtype` must be either `float32` or `float64`."
                raise ValueError(message)
            if self._dtype != inferred_dtype:
                values = (xp.asarray(val, dtype=self._dtype, device=self._device)
                          for val in parameters.values())
                parameters = dict(zip(parameters.keys(), values))

        self.reset_cache()
        self._parameters = parameters
        self._parameterization = parameterization
        self._original_parameters = original_parameters
        for name in self._parameters.keys():
            # Make parameters properties of the class; return values from the instance
            if hasattr(self.__class__, name):
                continue
            setattr(self.__class__, name, property(lambda self_, name_=name:
                                                   xp_copy(self_._parameters[name_])[()]))

    def __getstate__(self):
        state = self.__dict__.copy()
        if "_xp" in state:
            state["_xp"] = state["_xp"].asarray(1.)
        return state

    def __setstate__(self, state):
        if "_xp" in state:
            state["_xp"] = array_namespace(state["_xp"])
        self.__dict__.update(state)

    def reset_cache(self):
        r""" Clear all cached values.

        To improve the speed of some calculations, the distribution's support
        and moments are cached.

        This function is called automatically whenever the distribution
        parameters are updated.

        """
        # We could offer finer control over what is cleared.
        # For simplicity, these will still exist even if cache_policy is
        # NO_CACHE; they just won't be populated. This allows caching to be
        # turned on and off easily.
        self._moment_raw_cache = {}
        self._moment_central_cache = {}
        self._moment_standardized_cache = {}
        self._lmoment_cache = {}
        self._support_cache = None
        self._method_cache = {}
        self._constant_cache = None

    def _identify_parameterization(self, parameters):
        # Determine whether a `parameters` dictionary matches is consistent
        # with one of the parameterizations of the distribution. If so,
        # return that parameterization object; if not, raise an error.
        #
        # I've come back to this a few times wanting to avoid this explicit
        # loop. I've considered several possibilities, but they've all been a
        # little unusual. For example, we could override `_eq_` so we can
        # use _parameterizations.index() to retrieve the parameterization,
        # or the user could put the parameterizations in a dictionary so we
        # could look them up with a key (e.g. frozenset of parameter names).
        # I haven't been sure enough of these approaches to implement them.
        parameter_names_set = set(parameters)

        for parameterization in self._parameterizations:
            if parameterization.matches(parameter_names_set):
                break
        else:
            if not parameter_names_set:
                message = (f"The `{self.__class__.__name__}` distribution "
                           "family requires parameters, but none were "
                           "provided.")
            else:
                parameter_names = self._get_parameter_str(parameters)
                message = (f"The provided parameters `{parameter_names}` "
                           "do not match a supported parameterization of the "
                           f"`{self.__class__.__name__}` distribution family.")
            raise ValueError(message)

        return parameterization

    def _broadcast(self, parameters, *, dtype=None, device=None):
        # Broadcast the distribution parameters to the same shape. If the
        # arrays are not broadcastable, raise a meaningful error.
        #
        # We always make sure that the parameters *are* the same shape
        # and not just broadcastable. Users can access parameters as
        # attributes, and I think they should see the arrays as the same shape.
        # More importantly, arrays should be the same shape before logical
        # indexing operations, which are needed in infrastructure code when
        # there are invalid parameters, and may be needed in
        # distribution-specific code. We don't want developers to need to
        # broadcast in implementation functions.

        # It's much faster to check whether broadcasting is necessary than to
        # broadcast when it's not necessary.
        xp = self._xp
        parameter_vals = [xp.asarray(parameter, dtype=dtype, device=device)
                          for parameter in parameters.values()]
        parameter_shapes = set(parameter.shape for parameter in parameter_vals)
        if len(parameter_shapes) == 1:
            parameters = dict(zip(parameters.keys(), parameter_vals))
            return (parameters, parameter_vals[0].shape,
                    xp_size(parameter_vals[0]), parameter_vals[0].ndim)

        try:
            parameter_vals = xp.broadcast_arrays(*parameter_vals)
        except (ValueError, RuntimeError) as e:
            parameter_names = self._get_parameter_str(parameters)
            message = (f"The parameters `{parameter_names}` provided to the "
                       f"`{self.__class__.__name__}` distribution family "
                       "cannot be broadcast to the same shape.")
            raise ValueError(message) from e
        return (dict(zip(parameters.keys(), parameter_vals)),
                parameter_vals[0].shape,
                parameter_vals[0].size,
                parameter_vals[0].ndim)

    def _validate(self, parameterization, parameters):
        # Broadcasts distribution parameter arrays and converts them to a
        # consistent dtype. Replaces invalid parameters with `np.nan`.
        # Returns the validated parameters, a boolean mask indicated *which*
        # elements are invalid, a boolean scalar indicating whether *any*
        # are invalid (to skip special treatments if none are invalid), and
        # the common dtype.
        xp = self._xp
        values = xp_promote(*parameters.values(), force_floating=True, xp=xp)
        values = (values,) if not isinstance(values, tuple) else values
        device = xp_result_device(*parameters.values())
        parameters = dict(zip(parameters.keys(), values))
        valid = parameterization.validation(parameters, xp=xp)
        invalid = ~valid
        any_invalid = invalid if invalid.shape == () else xp.any(invalid)
        # If necessary, make the arrays contiguous and replace invalid with NaN
        if any_invalid:
            for parameter_name in parameters:
                parameters[parameter_name] = xp_copy(parameters[parameter_name])
                parameters[parameter_name][invalid] = xp.nan

        return parameters, invalid, any_invalid, values[0].dtype, device

    def _process_parameters(self, **params):
        r""" Process and cache distribution parameters for reuse.

        This is intended to be overridden by subclasses. It allows distribution
        authors to pre-process parameters for re-use. For instance, when a user
        parameterizes a LogUniform distribution with `a` and `b`, it makes
        sense to calculate `log(a)` and `log(b)` because these values will be
        used in almost all distribution methods. The dictionary returned by
        this method is passed to all private methods that calculate functions
        of the distribution.
        """
        return params

    def _get_parameter_str(self, parameters):
        # Get a string representation of the parameters like "{a, b, c}".
        return f"{{{', '.join(parameters.keys())}}}"

    def _copy_parameterization(self):
        self._parameterizations = self._parameterizations.copy()
        for i in range(len(self._parameterizations)):
            self._parameterizations[i] = self._parameterizations[i].copy()

    ### Attributes

    @property
    def tol(self):
        r"""positive float:
        The desired relative tolerance of calculations. Left unspecified,
        calculations may be faster; when provided, calculations may be
        more likely to meet the desired accuracy.
        """
        return self._tol

    @tol.setter
    def tol(self, tol):
        xp = self._xp

        if _isnull(tol):
            self._tol = tol
            return

        tol = xp.asarray(tol)
        if (tol.shape != () or not tol > 0 or  # catches NaNs
                not xp.isdtype(tol.dtype, "real floating")):
            message = (f"Attribute `tol` of `{self.__class__.__name__}` must "
                       "be a positive float, if specified.")
            raise ValueError(message)
        self._tol = float(tol[()])

    @property
    def cache_policy(self):
        r"""{None, "no_cache"}:
        Specifies the extent to which intermediate results are cached. Left
        unspecified, intermediate results of some calculations (e.g. distribution
        support, moments, etc.) are cached to improve performance of future
        calculations. Pass ``'no_cache'`` to reduce memory reserved by the class
        instance.
        """
        return self._cache_policy

    @cache_policy.setter
    def cache_policy(self, cache_policy):
        cache_policy = str(cache_policy).lower() if cache_policy is not None else None
        cache_policies = {None, 'no_cache'}
        if cache_policy not in cache_policies:
            message = (f"Attribute `cache_policy` of `{self.__class__.__name__}` "
                       f"must be one of {cache_policies}, if specified.")
            raise ValueError(message)
        self._cache_policy = cache_policy

    @property
    def validation_policy(self):
        r"""{None, "skip_all"}:
        Specifies the level of input validation to perform. Left unspecified,
        input validation is performed to ensure appropriate behavior in edge
        case (e.g. parameters out of domain, argument outside of distribution
        support, etc.) and improve consistency of output dtype, shape, etc.
        Use ``'skip_all'`` to avoid the computational overhead of these
        checks when rough edges are acceptable.
        """
        return self._validation_policy

    @validation_policy.setter
    def validation_policy(self, validation_policy):
        validation_policy = (str(validation_policy).lower()
                             if validation_policy is not None else None)
        iv_policies = {None, 'skip_all'}
        if validation_policy not in iv_policies:
            message = (f"Attribute `validation_policy` of `{self.__class__.__name__}` "
                       f"must be one of {iv_policies}, if specified.")
            raise ValueError(message)
        self._validation_policy = validation_policy

    ### Other magic methods

    def __repr__(self):
        r""" Returns a string representation of the distribution.

        Includes the name of the distribution family, the names of the
        parameters and the `repr` of each of their values.


        """
        class_name = self.__class__.__name__
        parameters = list(self._original_parameters.items())
        info = []
        with np.printoptions(threshold=10):
            str_parameters = [f"{symbol}={repr(value)}" for symbol, value in parameters]
        str_parameters = f"{', '.join(str_parameters)}"
        info.append(str_parameters)
        return f"{class_name}({', '.join(info)})"

    def __str__(self):
        class_name = self.__class__.__name__
        parameters = list(self._original_parameters.items())
        info = []
        with np.printoptions(threshold=10):
            str_parameters = [f"{symbol}={str(value)}" for symbol, value in parameters]
        str_parameters = f"{', '.join(str_parameters)}"
        info.append(str_parameters)
        return f"{class_name}({', '.join(info)})"

    def __add__(self, loc):
        return ShiftedScaledDistribution(self, loc=loc)

    def __sub__(self, loc):
        return ShiftedScaledDistribution(self, loc=-loc)

    def __mul__(self, scale):
        return ShiftedScaledDistribution(self, scale=scale)

    def __truediv__(self, scale):
        return ShiftedScaledDistribution(self, scale=1/scale)

    def __pow__(self, other):
        xp = self._xp
        if not np.isscalar(other) or other <= 0 or other != int(other):
            message = ("Raising a random variable to the power of an argument is only "
                       "implemented when the argument is a positive integer.")
            raise NotImplementedError(message)

        # Fill in repr_pattern with the repr of self before taking abs.
        # Avoids having unnecessary abs in the repr.
        with np.printoptions(threshold=10):
            repr_pattern = f"({repr(self)})**{repr(other)}"
            str_pattern = f"({str(self)})**{str(other)}"
        X = abs(self) if other % 2 == 0 else self

        funcs = dict(g=lambda u: u**other, repr_pattern=repr_pattern,
                     str_pattern=str_pattern,
                     h=lambda u: xp.sign(u) * xp.abs(u)**(1 / other),
                     dh=lambda u: 1/other * xp.abs(u)**(1/other - 1),
                     xp=xp)

        return MonotonicTransformedDistribution(X, **funcs, increasing=True)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rtruediv__(self, other):
        xp = self._xp
        a, b = self.support()
        with np.printoptions(threshold=10):
            funcs = dict(g=lambda u: 1 / u,
                         repr_pattern=f"{repr(other)}/({repr(self)})",
                         str_pattern=f"{str(other)}/({str(self)})",
                         h=lambda u: 1 / u, dh=lambda u: 1 / u ** 2,
                         xp=xp)
        if xp.all(a >= 0) or xp.all(b <= 0):
            out = MonotonicTransformedDistribution(self, **funcs, increasing=False)
        else:
            message = ("Division by a random variable is only implemented "
                       "when the support is either non-negative or non-positive.")
            raise NotImplementedError(message)
        if xp.all(xp.asarray(other) == 1):
            return out
        else:
            return out * other

    def __rpow__(self, other):
        xp = self._xp
        with np.printoptions(threshold=10):
            funcs = dict(g=lambda u: other**u,
                         h=lambda u: xp.log(u) / math.log(other),
                         dh=lambda u: 1 / xp.abs(u * math.log(other)),
                         repr_pattern=f"{repr(other)}**({repr(self)})",
                         str_pattern=f"{str(other)}**({str(self)})",
                         xp=xp)

        if not np.isscalar(other) or other <= 0 or other == 1:
            message = ("Raising an argument to the power of a random variable is only "
                       "implemented when the argument is a positive scalar other than "
                       "1.")
            raise NotImplementedError(message)

        if other > 1:
            return MonotonicTransformedDistribution(self, **funcs, increasing=True)
        else:
            return MonotonicTransformedDistribution(self, **funcs, increasing=False)

    def __neg__(self):
        return self * -1

    def __abs__(self):
        return FoldedDistribution(self)

    ### Utilities

    ## Input validation

    def _validate_order(self, order, fname='moment', min_order=0):
        # Yet another integer validating function. Unlike others in SciPy, it
        # Is quite flexible about what is allowed as an integer, and it
        # raises a distribution-specific error message to facilitate
        # identification of the source of the error.
        if self.validation_policy == _SKIP_ALL:
            return order

        message = (f"Argument `order` of `{self.__class__.__name__}.{fname}` "
                   f"must be a finite integer greater than or equal to {min_order}.")
        try:
            order_int = round(order)
            # If this fails for any reason (e.g. it's an array, it's infinite)
            # it's not a valid `order`.
        except Exception as e:
            raise ValueError(message) from e

        if order_int < min_order or order_int != order:
            raise ValueError(message)

        return float(order)

    def _validate_kind(self, kind, kinds):
        message = (f"Argument `kind` of `{self.__class__.__name__}.moment` "
                   f"must be one of {set(kinds)}.")
        kind = kind.lower()
        if kind not in kinds:
            raise ValueError(message)
        return kind

    def _preserve_type(self, x):
        xp = self._xp
        x = xp.asarray(x, device=self._device)
        if x.dtype != self._dtype:
            x = xp.astype(x, self._dtype)
        return x[()]

    ## Testing

    @classmethod
    def _draw(cls, sizes=None, rng=None, i_parameterization=None,
              proportions=None, *, xp):
        r""" Draw a specific (fully-defined) distribution from the family.

        See _Parameterization.draw for documentation details.
        """
        rng = np.random.default_rng(rng)
        if len(cls._parameterizations) == 0:
            # xpx.default_device may not be the device on which arrays are created;
            # see data-apis/array-api-compat/issues/489
            return cls(xp=xp, device=xp_device(xp.asarray(0)))  # skip device check
        if i_parameterization is None:
            n = cls._num_parameterizations()
            i_parameterization = rng.integers(0, max(0, n - 1), endpoint=True)

        parameterization = cls._parameterizations[i_parameterization]
        parameters = parameterization.draw(sizes, rng, proportions=proportions,
                                           region='typical', xp=xp)
        return cls(**parameters, xp=xp)

    @classmethod
    def _num_parameterizations(cls):
        # Returns the number of parameterizations accepted by the family.
        return len(cls._parameterizations)

    @classmethod
    def _num_parameters(cls, i_parameterization=0):
        # Returns the number of parameters used in the specified
        # parameterization.
        return (0 if not cls._num_parameterizations()
                else len(cls._parameterizations[i_parameterization]))

    ## Algorithms

    def _differentiation(self, f, x, bounds=None, args=None, params=None):
        xp = self._xp
        a, b = self._support(**params) if bounds is None else bounds
        x = xp.real(x)  # logentropy makes dtype complex
        half = xp.asarray(0.5, dtype=self._dtype, device=self._device)
        step = xp.minimum(half, (b - a)/2)
        direction = -(xp.sign(x - a - step) + xp.sign(x - b + step))
        args = [] if args is None else args
        params = {} if params is None else params
        args = xp.broadcast_arrays(*args)
        rtol = None if _isnull(self.tol) else self.tol
        res = derivative(f, x, initial_step=step, step_direction=xp.sign(direction),
                         args=args, kwargs=params, tolerances={'rtol': rtol})
        return res.df

    def _quadrature(self, integrand, limits=None, args=(), params=None, log=False):
        # Performs numerical integration/summation between limits or over support.
        xp = self._xp
        if limits is None:
            a, b = self._support(**params)
        else:
            # ideally, a and b could be left as floats, but `_tanhsinh` doesn't
            # know to get device from args/kwargs only
            a = xp.asarray(limits[0], dtype=self._dtype, device=self._device)
            b = xp.asarray(limits[1], dtype=self._dtype, device=self._device)
        rtol = None if _isnull(self.tol) else self.tol

        # cupy.broadcast_arrays accepts non-arrays but doesn't broadcast them.
        # `order` may be passed around in `args` as an `int`, so it wouldn't
        # be broadcast. That turns out to be a problem during integration.
        args = tuple(self._preserve_type(arg) for arg in args)

        # For now, we ignore the status, but I want to return the error estimate
        if isinstance(self, ContinuousDistribution):
            return _tanhsinh(integrand, a, b, args=args, kwargs=params,
                             log=log, rtol=rtol).integral
        else:
            res = nsum(integrand, a, b, args=args, kwargs=params,
                       log=log, tolerances=dict(rtol=rtol)).sum

            # b > a -> sum is zero... unless the parameters are nan.
            # (Resolving gh-22321 should fix this in nsum.)
            _true = xp.asarray(True, device=self._device)
            cond = xp.isnan(params.popitem()[1]) if params else _true
            return xp.where(cond, xp.nan, res)[()]

    def _solve_bounded(self, f, p, *, bounds=None, params=None, xatol=None):
        # Finds the argument of a function that produces the desired output.
        xp = self._xp
        p = xp.real(p)

        xmin, xmax = self._support(**params) if bounds is None else bounds
        xmin = xp.asarray(xmin, dtype=self._dtype)
        xmax = xp.asarray(xmax, dtype=self._dtype)

        def f2(x, _p, **kwargs):  # named `_p` to avoid conflict with shape `p`
            return f(x, **kwargs) - _p

        xl0, xr0 = _guess_bracket(xmin, xmax, xp=xp)

        res = elementwise.bracket_root(f2, xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax,
                                       args=(p,), kwargs=params)
        # For now, we ignore the status, but I want to use the bracket width
        # as an error estimate.

        xrtol = None if _isnull(self.tol) else self.tol
        xatol = None if xatol is None else xatol
        tolerances = dict(xrtol=xrtol, xatol=xatol, fatol=0, frtol=0)
        return elementwise.find_root(f2, res.bracket, args=(p,),
                                     kwargs=params, tolerances=tolerances)

    def _optimization(self, f, x0, xatol, params):
        xp = self._xp
        if not self._size:
            return xp.empty(self._shape, dtype=self._dtype, device=self._device)

        a, b = self._support(**params)

        res_b = elementwise.bracket_minimum(f, x0, xmin=a, xmax=b, kwargs=params)
        res = elementwise.find_minimum(f, res_b.bracket, kwargs=params,
                                       tolerances=dict(xatol=xatol))
        x = xp.asarray(res.x)

        # If the optimum is at an endpoint, `_bracket_minimum` cannot produce a valid
        # bracket; it may terminate with `fl < fm < fr` (and, e.g. `xl < xm < xr` but
        # all very close to `a`. In this case, we assume the function is unimodal, and
        # the optimum is at the endpoint.
        # TODO: see https://github.com/scipy/scipy/pull/25210#discussion_r3295327858
        x_at_boundary = res.status == -1
        fl, fm, fr = res_b.f_bracket
        x_at_left = x_at_boundary & (fl <= fm)
        x_at_right = x_at_boundary & (fr < fm)
        x[x_at_left] = a[x_at_left]
        x[x_at_right] = b[x_at_right]
        return x[()]

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
        super_method = getattr(UnivariateDistribution, method_name, None)
        return method is not super_method

    ### Distribution properties
    # The following "distribution properties" are exposed via a public method
    # that accepts only options (not distribution parameters or quantile/
    # percentile argument).
    # support
    # logentropy, entropy,
    # median, mode, mean,
    # variance, standard_deviation
    # skewness, kurtosis
    # Common options are:
    # method - a string that indicates which method should be used to compute
    #          the quantity (e.g. a formula or numerical integration).
    # Input/output validation is provided by the `_set_invalid_nan_property`
    # decorator. These are the methods meant to be called by users.
    #
    # Each public method calls a private "dispatch" method that
    # determines which "method" (strategy for calculating the desired quantity)
    # to use by default and, via the `@_dispatch` decorator, calls the
    # method and computes the result.
    # Dispatch methods always accept:
    # method - as passed from the public method
    # params - a dictionary of distribution shape parameters passed by
    #          the public method.
    # Dispatch methods accept `params` rather than relying on the state of the
    # object because iterative algorithms like `tanhsinh` and `find_root`
    # need their callable to follow a strict elementwise protocol: each element
    # of the output is determined solely by the values of the inputs at the
    # corresponding location. The public methods do not satisfy this protocol
    # because they do not accept the parameters as arguments, producing an
    # output that generally has a different shape than that of the input. Also,
    # by calling "dispatch" methods rather than the public methods, the
    # iterative algorithms avoid the overhead of input validation.
    #
    # Each dispatch method can designate the responsibility of computing
    # the required value to any of several "implementation" methods. These
    # methods accept only `**params`, the parameter dictionary passed from
    # the public method via the dispatch method. We separate the implementation
    # methods from the dispatch methods for the sake of simplicity (via
    # compartmentalization) and to allow subclasses to override certain
    # implementation methods (typically only the "formula" methods). The names
    # of implementation methods are combinations of the public method name and
    # the name of the "method" (strategy for calculating the desired quantity)
    # string. (In fact, the name of the implementation method is calculated
    # from these two strings in the `_dispatch` decorator.) Common method
    # strings are:
    # formula - distribution-specific analytical expressions to be implemented
    #           by subclasses.
    # log/exp - Compute the log of a number and then exponentiate it or vice
    #           versa.
    # quadrature - Compute the value via numerical integration.
    #
    # The default method (strategy) is determined based on what implementation
    # methods are available and the error tolerance of the user. Typically,
    # a formula is always used if available. We fall back to "log/exp" if a
    # formula for the logarithm or exponential of the quantity is available,
    # and we use quadrature otherwise.

    def support(self):
        # If this were a `cached_property`, we couldn't update the value
        # when the distribution parameters change.
        # Caching is important, though, because calls to _support take a few
        # microseconds even when `a` and `b` are already the same shape.
        xp = self._xp

        if self._support_cache is not None:
            return self._support_cache

        a, b = self._support(**self._parameters)
        if a.shape != self._shape:
            a = xp.broadcast_to(a, self._shape)
        if b.shape != self._shape:
            b = xp.broadcast_to(b, self._shape)

        if self._any_invalid:
            a, b = xp.asarray(a, copy=True), xp.asarray(b, copy=True)
            a[self._invalid], b[self._invalid] = xp.nan, xp.nan
            a, b = a[()], b[()]

        support = (a, b)

        if self.cache_policy != _NO_CACHE:
            self._support_cache = support

        return support

    def _support(self, **params):
        # Computes the support given distribution parameters
        xp = self._xp
        a, b = self._variable.domain.get_numerical_endpoints(
            params, xp=xp, device=self._device)
        if len(params):
            # the parameters should all be of the same dtype and shape at this point
            vals = list(params.values())
            shape = vals[0].shape
            a = xp.broadcast_to(a, shape) if a.shape != shape else a
            b = xp.broadcast_to(b, shape) if b.shape != shape else b
        return self._preserve_type(a), self._preserve_type(b)

    @_set_invalid_nan_property
    def logentropy(self, *, method=None):
        return self._logentropy_dispatch(method=method, **self._parameters) + 0j

    @_dispatch
    def _logentropy_dispatch(self, method=None, **params):
        if self._overrides('_logentropy_formula'):
            method = self._logentropy_formula
        elif self._overrides('_entropy_formula'):
            method = self._logentropy_logexp_safe
        else:
            method = self._logentropy_quadrature
        return method

    def _logentropy_formula(self, **params):
        raise NotImplementedError(self._not_implemented)

    def _logentropy_logexp(self, **params):
        xp = self._xp
        with np.errstate(invalid='ignore', divide='ignore'):
            # np.log warns with complex NaN argument
            entropy = xp.asarray(self._entropy_dispatch(**params),
                                 dtype=self._dtype, device=self._device)
            res = xp.log(entropy+0j)
            return _log_real_standardize(res, xp=xp)

    def _logentropy_logexp_safe(self, **params):
        xp = self._xp
        out = self._logentropy_logexp(**params)
        mask = xp.isinf(xp.real(out))
        if xp.any(mask):
            params_mask = {key:val[mask] for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._logentropy_quadrature(**params_mask)
        return out[()]

    def _logentropy_quadrature(self, **params):
        xp = self._xp
        def logintegrand(x, **params):
            logpxf = self._logpxf_dispatch(x, **params)
            return logpxf + xp.log(0j+logpxf)
        res = self._quadrature(logintegrand, params=params, log=True)
        return _log_real_standardize(res + xp.pi*1j, xp=xp)

    @_set_invalid_nan_property
    def entropy(self, *, method=None):
        return self._entropy_dispatch(method=method, **self._parameters)

    @_dispatch
    def _entropy_dispatch(self, method=None, **params):
        if self._overrides('_entropy_formula'):
            method = self._entropy_formula
        elif self._overrides('_logentropy_formula'):
            method = self._entropy_logexp
        else:
            method = self._entropy_quadrature
        return method

    def _entropy_formula(self, **params):
        raise NotImplementedError(self._not_implemented)

    def _entropy_logexp(self, **params):
        xp = self._xp
        with np.errstate(invalid='ignore'):
            # np.exp(np.nan) raises on some platforms?
            dtype = xp.result_type(self._dtype, 0j)
            logentropy = xp.asarray(self._logentropy_dispatch(**params),
                                    dtype=dtype, device=self._device)
            return xp.real(xp.exp(logentropy))

    def _entropy_quadrature(self, **params):
        xp = self._xp
        def integrand(x, **params):
            pxf = self._pxf_dispatch(x, **params)
            logpxf = self._logpxf_dispatch(x, **params)
            temp = xp.asarray(pxf)
            i = (pxf != 0)  # 0 * inf -> nan; should be 0
            temp[i] = -pxf[i]*logpxf[i]
            return temp
        return self._quadrature(integrand, params=params)

    @_set_invalid_nan_property
    def median(self, *, method=None):
        return self._median_dispatch(method=method, **self._parameters)

    @_dispatch
    def _median_dispatch(self, method=None, **params):
        if self._overrides('_median_formula'):
            method = self._median_formula
        elif isinstance(self, CircularDistribution):
            method = self._median_optimization
        else:
            method = self._median_icdf
        return method

    def _median_formula(self, **params):
        raise NotImplementedError(self._not_implemented)

    def _median_icdf(self, **params):
        xp = self._xp
        half = xp.asarray(0.5, dtype=self._dtype, device=self._device)
        return self._icdf_dispatch(half, **params)

    @_set_invalid_nan_property
    def mode(self, *, method=None):
        return self._mode_dispatch(method=method, **self._parameters)

    @_dispatch
    def _mode_dispatch(self, method=None, **params):
        # We could add a method that looks for a critical point with
        # differentiation and the root finder
        if self._overrides('_mode_formula'):
            method = self._mode_formula
        else:
            method = self._mode_optimization
        return method

    def _mode_formula(self, **params):
        raise NotImplementedError(self._not_implemented)

    def _mode_optimization(self, xatol=None, **params):
        m = self._median_dispatch(**params)
        def f(x, **params): return -self._pxf_dispatch(x, **params)
        return self._optimization(f, m, xatol, params)

    def mean(self, *, method=None):
        return self.moment(1, kind='raw', method=method)

    def variance(self, *, method=None):
        return self.moment(2, kind='central', method=method)

    def standard_deviation(self, *, method=None):
        return self.variance(method=method)**0.5

    def skewness(self, *, method=None):
        return self.moment(3, kind='standardized', method=method)

    def kurtosis(self, *, method=None, convention=None):
        convention = 'non-excess' if convention is None else convention
        conventions = {'non-excess', 'excess'}
        message = (f'Parameter `convention` of `{self.__class__.__name__}.kurtosis` '
                   f"must be one of {conventions}.")
        convention = convention.lower()
        if convention not in conventions:
            raise ValueError(message)
        k = self.moment(4, kind='standardized', method=method)
        return k - 3 if convention == 'excess' else k

    ### Distribution functions
    # The following functions related to the distribution PDF and CDF are
    # exposed via a public method that accepts one positional argument - the
    # quantile - and keyword options (but not distribution parameters).
    # logpdf, pdf
    # logcdf, cdf
    # logccdf, ccdf
    # The `logcdf` and `cdf` functions can also be called with two positional
    # arguments - lower and upper quantiles - and they return the probability
    # mass (integral of the PDF) between them. The 2-arg versions of `logccdf`
    # and `ccdf` return the complement of this quantity.
    # All the (1-arg) cumulative distribution functions have inverse
    # functions, which accept one positional argument - the percentile.
    # ilogcdf, icdf
    # ilogccdf, iccdf
    # Common keyword options include:
    # method - a string that indicates which method should be used to compute
    #          the quantity (e.g. a formula or numerical integration).
    # Tolerance options should be added.
    # Input/output validation is provided by the `_set_invalid_nan`
    # decorator. These are the methods meant to be called by users.
    #
    # Each public method calls a private "dispatch" method that
    # determines which "method" (strategy for calculating the desired quantity)
    # to use by default and, via the `@_dispatch` decorator, calls the
    # method and computes the result.
    # Each dispatch method can designate the responsibility of computing
    # the required value to any of several "implementation" methods. These
    # methods accept only `**params`, the parameter dictionary passed from
    # the public method via the dispatch method.
    # See the note corresponding with the "Distribution Parameters" for more
    # information.

    ## Probability Density/Mass Functions

    @_set_invalid_nan
    def logpdf(self, x, /, *, method=None):
        return self._logpdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _logpdf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_logpdf_formula'):
            method = self._logpdf_formula
        elif _isnull(self.tol):  # ensure that developers override _logpdf
            method = self._logpdf_logexp
        return method

    def _logpdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _logpdf_logexp(self, x, **params):
        xp = self._xp
        with np.errstate(divide='ignore'):
            return xp.log(self._pdf_dispatch(x, **params))

    @_set_invalid_nan
    def pdf(self, x, /, *, method=None):
        return self._pdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _pdf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_pdf_formula'):
            method = self._pdf_formula
        elif self._overrides('_logpdf_formula') or self._overrides('_logpdf_dispatch'):
            method = self._pdf_logexp
        elif isinstance(self, ContinuousDistribution):
            if self._overrides('_cdf_formula'):
                method = self._pdf_differentiation_cdf
            elif self._overrides('_ccdf_formula'):
                method = self._pdf_differentiation_ccdf
            elif self._overrides('_icdf_formula'):
                method = self._pdf_differentiation_icdf
            elif self._overrides('_iccdf_formula'):
                method = self._pdf_differentiation_iccdf

        return method

    def _pdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _pdf_differentiation_cdf(self, x, **params):
        return self._differentiation(self._cdf_dispatch, x, params=params)

    def _pdf_differentiation_ccdf(self, x, **params):
        return -self._differentiation(self._ccdf_dispatch, x, params=params)

    def _pdf_differentiation_icdf(self, x, **params):
        p = self._cdf_dispatch(x, **params)
        return 1 / self._differentiation(self._icdf_dispatch, p,
                                         bounds=(0, 1), params=params)

    def _pdf_differentiation_iccdf(self, x, **params):
        p = self._ccdf_dispatch(x, **params)
        return -1 / self._differentiation(self._iccdf_dispatch, p,
                                          bounds=(0, 1), params=params)

    def _pdf_logexp(self, x, **params):
        xp = self._xp
        return xp.exp(self._logpdf_dispatch(x, **params))

    @_set_invalid_nan
    def logpmf(self, x, /, *, method=None):
        return self._logpmf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _logpmf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_logpmf_formula'):
            method = self._logpmf_formula
        elif _isnull(self.tol):  # ensure that developers override _logpmf
            method = self._logpmf_logexp
        return method

    def _logpmf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _logpmf_logexp(self, x, **params):
        xp = self._xp
        with np.errstate(divide='ignore'):
            return xp.log(self._pmf_dispatch(x, **params))

    @_set_invalid_nan
    def pmf(self, x, /, *, method=None):
        return self._pmf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _pmf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_pmf_formula'):
            method = self._pmf_formula
        else:
            method = self._pmf_logexp
        return method

    def _pmf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _pmf_logexp(self, x, **params):
        xp = self._xp
        return xp.exp(self._logpmf_dispatch(x, **params))

    ## Cumulative Distribution Functions

    def logcdf(self, x, y=None, /, *, method=None):
        if y is None:
            return self._logcdf1(x, method=method)
        else:
            return self._logcdf2(x, y, method=method)

    @_cdf2_input_validation
    def _logcdf2(self, x, y, *, method):
        xp = self._xp
        out = self._logcdf2_dispatch(x, y, method=method, **self._parameters)
        return (out + 0j) if not xp.isdtype(out.dtype, "complex floating") else out

    @_dispatch
    def _logcdf2_dispatch(self, x, y, *, method=None, **params):
        if self._overrides('_logcdf2_formula'):
            method = self._logcdf2_formula
        elif (self._overrides('_logcdf_formula')
              or self._overrides('_logccdf_formula')):
            method = self._logcdf2_subtraction
        elif (self._overrides('_cdf_formula')
              or self._overrides('_ccdf_formula')):
            method = self._logcdf2_logexp_safe
        else:
            method = self._logcdf2_quadrature
        return method

    def _logcdf2_formula(self, x, y, **params):
        raise NotImplementedError(self._not_implemented)

    def _logcdf2_subtraction(self, x, y, **params):
        xp = self._xp
        flip_sign = x > y  # some results will be negative
        x, y = xp.minimum(x, y), xp.maximum(x, y)
        logcdf_x = self._logcdf_dispatch(x, **params)
        logcdf_y = self._logcdf_dispatch(y, **params)
        logccdf_x = self._logccdf_dispatch(x, **params)
        logccdf_y = self._logccdf_dispatch(y, **params)
        case_left = (logcdf_x < -1) & (logcdf_y < -1)
        case_right = (logccdf_x < -1) & (logccdf_y < -1)
        case_central = ~(case_left | case_right)
        log_mass = _logexpxmexpy(logcdf_y, logcdf_x, xp=xp)
        log_mass[case_right] = _logexpxmexpy(logccdf_x, logccdf_y, xp=xp)[case_right]
        with np.errstate(invalid='ignore'):
            # np.logaddexp warns with NaN argument
            log_tail = xp.logaddexp(logcdf_x, logccdf_y)[case_central]
        log_mass[case_central] = xp.astype(_log1mexp(log_tail, xp=xp), log_mass.dtype)
        log_mass[flip_sign] += xp.pi * 1j
        return log_mass[()] if xp.any(flip_sign) else xp.real(log_mass)[()]

    def _logcdf2_logexp(self, x, y, **params):
        xp = self._xp
        expres = self._cdf2_dispatch(x, y, **params)
        expres = expres + 0j if xp.any(x > y) else expres
        with np.errstate(divide='ignore'):
            return xp.log(expres)

    def _logcdf2_logexp_safe(self, x, y, **params):
        xp = self._xp
        out = self._logcdf2_logexp(x, y, **params)
        mask = xp.isinf(xp.real(out))
        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._logcdf2_quadrature(x[mask], y[mask], **params_mask)
        return out[()]

    def _logcdf2_quadrature(self, x, y, **params):
        logres = self._quadrature(self._logpxf_dispatch, limits=(x, y),
                                  log=True, params=params)
        return logres

    @_set_invalid_nan
    def _logcdf1(self, x, *, method=None):
        return self._logcdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _logcdf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_logcdf_formula'):
            method = self._logcdf_formula
        elif self._overrides('_logccdf_formula'):
            method = self._logcdf_complement
        elif self._overrides('_cdf_formula'):
            method = self._logcdf_logexp_safe
        elif (self._overrides('_ilogcdf_formula')
              and isinstance(self, ContinuousDistribution)):
            method = self._logcdf_inversion
        else:
            method = self._logcdf_quadrature
        return method

    def _logcdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _logcdf_complement(self, x, **params):
        xp = self._xp
        return _log1mexp(self._logccdf_dispatch(x, **params), xp=xp)

    def _logcdf_logexp(self, x, **params):
        xp = self._xp
        with np.errstate(divide='ignore'):
            return xp.log(self._cdf_dispatch(x, **params))

    def _logcdf_logexp_safe(self, x, **params):
        xp = self._xp
        out = self._logcdf_logexp(x, **params)
        mask = xp.isinf(out)
        if xp.any(mask):
            params_mask = {key:xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._logcdf_quadrature(x[mask], **params_mask)
        return out[()]

    def _logcdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._ilogcdf_dispatch, x,
                                              bounds=(-inf, 0), params=params)

    def _logcdf_quadrature(self, x, **params):
        a, _ = self._support(**params)
        return self._quadrature(self._logpxf_dispatch, limits=(a, x),
                                params=params, log=True)

    def cdf(self, x, y=None, /, *, method=None):
        if y is None:
            return self._cdf1(x, method=method)
        else:
            return self._cdf2(x, y, method=method)

    @_cdf2_input_validation
    def _cdf2(self, x, y, *, method):
        return self._cdf2_dispatch(x, y, method=method, **self._parameters)

    @_dispatch
    def _cdf2_dispatch(self, x, y, *, method=None, **params):
        if self._overrides('_cdf2_formula'):
            method = self._cdf2_formula
        elif (self._overrides('_logcdf_formula')
              or self._overrides('_logccdf_formula')):
            method = self._cdf2_logexp
        elif self._overrides('_cdf_formula') or self._overrides('_ccdf_formula'):
            method = self._cdf2_subtraction_safe
        else:
            method = self._cdf2_quadrature
        return method

    def _cdf2_formula(self, x, y, **params):
        raise NotImplementedError(self._not_implemented)

    def _cdf2_logexp(self, x, y, **params):
        xp = self._xp
        return xp.real(xp.exp(self._logcdf2_dispatch(x, y, **params)))

    def _cdf2_subtraction(self, x, y, **params):
        # Improvements:
        # Lazy evaluation of cdf/ccdf only where needed
        # Stack x and y to reduce function calls?
        xp = self._xp
        cdf_x = self._cdf_dispatch(x, **params)
        cdf_y = self._cdf_dispatch(y, **params)
        ccdf_x = self._ccdf_dispatch(x, **params)
        ccdf_y = self._ccdf_dispatch(y, **params)
        i = (ccdf_x < 0.5) & (ccdf_y < 0.5)
        return xp.where(i, ccdf_x-ccdf_y, cdf_y-cdf_x)

    def _cdf2_subtraction_safe(self, x, y, **params):
        # The "safe" version of this function may be used if `method` is unspecified,
        # but if `method='subtraction'`, the regular version  above is used.
        xp = self._xp
        cdf_x = self._cdf_dispatch(x, **params)
        cdf_y = self._cdf_dispatch(y, **params)
        ccdf_x = self._ccdf_dispatch(x, **params)
        ccdf_y = self._ccdf_dispatch(y, **params)
        i = (ccdf_x < 0.5) & (ccdf_y < 0.5)
        out = xp.where(i, ccdf_x-ccdf_y, cdf_y-cdf_x)
        # Can't just call `out = _cdf2_subtraction(self, x, y, **params)` here
        # because we need the partial results below. Could refactor, but we'll leave
        # that to future work, say if the improvements to _cdf2_subtraction are made.

        eps = xp.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else eps**0.5

        cdf_max = xp.maximum(cdf_x, cdf_y)
        ccdf_max = xp.maximum(ccdf_x, ccdf_y)
        with np.errstate(invalid='ignore'):
            spacing = xp_spacing(xp.where(i, ccdf_max, cdf_max))
        mask = xp.abs(tol * out) < spacing

        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._cdf2_quadrature(x[mask], y[mask], **params_mask)
        return out[()]

    def _cdf2_quadrature(self, x, y, **params):
        return self._quadrature(self._pxf_dispatch, limits=(x, y), params=params)

    @_set_invalid_nan
    def _cdf1(self, x, *, method):
        return self._cdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _cdf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_cdf_formula'):
            method = self._cdf_formula
        elif self._overrides('_logcdf_formula'):
            method = self._cdf_logexp
        elif self._overrides('_ccdf_formula'):
            method = self._cdf_complement_safe
        elif (self._overrides('_icdf_formula')
              and isinstance(self, ContinuousDistribution)):
            method = self._cdf_inversion
        else:
            method = self._cdf_quadrature
        return method

    def _cdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _cdf_logexp(self, x, **params):
        xp = self._xp
        return xp.exp(self._logcdf_dispatch(x, **params))

    def _cdf_complement(self, x, **params):
        return 1 - self._ccdf_dispatch(x, **params)

    def _cdf_complement_safe(self, x, **params):
        xp = self._xp
        ccdf = self._ccdf_dispatch(x, **params)
        out = 1 - ccdf
        eps = xp.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else eps**0.5
        with np.errstate(invalid='ignore'):
            mask = tol * out < xp_spacing(ccdf)
        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._cdf_quadrature(x[mask], **params_mask)
        return out[()]

    def _cdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._icdf_dispatch, x,
                                              bounds=(0, 1), params=params)

    def _cdf_quadrature(self, x, **params):
        a, _ = self._support(**params)
        return self._quadrature(self._pxf_dispatch, limits=(a, x),
                                params=params)

    def logccdf(self, x, y=None, /, *, method=None):
        if y is None:
            return self._logccdf1(x, method=method)
        else:
            return self._logccdf2(x, y, method=method)

    @_cdf2_input_validation
    def _logccdf2(self, x, y, *, method):
        return self._logccdf2_dispatch(x, y, method=method, **self._parameters)

    @_dispatch
    def _logccdf2_dispatch(self, x, y, *, method=None, **params):
        # if _logccdf2_formula exists, we could use the complement
        # if _ccdf2_formula exists, we could use log/exp
        if self._overrides('_logccdf2_formula'):
            method = self._logccdf2_formula
        else:
            method = self._logccdf2_addition
        return method

    def _logccdf2_formula(self, x, y, **params):
        raise NotImplementedError(self._not_implemented)

    def _logccdf2_addition(self, x, y, **params):
        xp = self._xp
        logcdf_x = self._logcdf_dispatch(x, **params)
        logccdf_y = self._logccdf_dispatch(y, **params)
        return special.logsumexp(xp.stack((logcdf_x, logccdf_y)), axis=0)

    @_set_invalid_nan
    def _logccdf1(self, x, *, method=None):
        return self._logccdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _logccdf_dispatch(self, x, method=None, **params):
        if self._overrides('_logccdf_formula'):
            method = self._logccdf_formula
        elif self._overrides('_logcdf_formula'):
            method = self._logccdf_complement
        elif self._overrides('_ccdf_formula'):
            method = self._logccdf_logexp_safe
        elif (self._overrides('_ilogccdf_formula')
              and isinstance(self, ContinuousDistribution)):
            method = self._logccdf_inversion
        else:
            method = self._logccdf_quadrature
        return method

    def _logccdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _logccdf_complement(self, x, **params):
        return _log1mexp(self._logcdf_dispatch(x, **params))

    def _logccdf_logexp(self, x, **params):
        xp = self._xp
        return xp.log(self._ccdf_dispatch(x, **params))

    def _logccdf_logexp_safe(self, x, **params):
        xp = self._xp
        out = self._logccdf_logexp(x, **params)
        mask = xp.isinf(out)
        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._logccdf_quadrature(x[mask], **params_mask)
        return out[()]

    def _logccdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._ilogccdf_dispatch, x,
                                              bounds=(-inf, 0), params=params)

    def _logccdf_quadrature(self, x, **params):
        _, b = self._support(**params)
        return self._quadrature(self._logpxf_dispatch, limits=(x, b),
                                params=params, log=True)

    def ccdf(self, x, y=None, /, *, method=None):
        if y is None:
            return self._ccdf1(x, method=method)
        else:
            return self._ccdf2(x, y, method=method)

    @_cdf2_input_validation
    def _ccdf2(self, x, y, *, method):
        return self._ccdf2_dispatch(x, y, method=method, **self._parameters)

    @_dispatch
    def _ccdf2_dispatch(self, x, y, *, method=None, **params):
        if self._overrides('_ccdf2_formula'):
            method = self._ccdf2_formula
        else:
            method = self._ccdf2_addition
        return method

    def _ccdf2_formula(self, x, y, **params):
        raise NotImplementedError(self._not_implemented)

    def _ccdf2_addition(self, x, y, **params):
        cdf_x = self._cdf_dispatch(x, **params)
        ccdf_y = self._ccdf_dispatch(y, **params)
        # even if x > y, cdf(x, y) + ccdf(x,y) sums to 1
        return cdf_x + ccdf_y

    @_set_invalid_nan
    def _ccdf1(self, x, *, method):
        return self._ccdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _ccdf_dispatch(self, x, method=None, **params):
        if self._overrides('_ccdf_formula'):
            method = self._ccdf_formula
        elif self._overrides('_logccdf_formula'):
            method = self._ccdf_logexp
        elif self._overrides('_cdf_formula'):
            method = self._ccdf_complement_safe
        elif (self._overrides('_iccdf_formula')
              and isinstance(self, ContinuousDistribution)):
            method = self._ccdf_inversion
        else:
            method = self._ccdf_quadrature
        return method

    def _ccdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _ccdf_logexp(self, x, **params):
        xp = self._xp
        return xp.exp(self._logccdf_dispatch(x, **params))

    def _ccdf_complement(self, x, **params):
        return 1 - self._cdf_dispatch(x, **params)

    def _ccdf_complement_safe(self, x, **params):
        xp = self._xp
        cdf = self._cdf_dispatch(x, **params)
        out = 1 - cdf
        eps = xp.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else eps**0.5
        with np.errstate(invalid='ignore'):
            mask = tol * out < xp_spacing(cdf)
        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._ccdf_quadrature(x[mask], **params_mask)
        return out[()]

    def _ccdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._iccdf_dispatch, x,
                                              bounds=(0, 1), params=params)

    def _ccdf_quadrature(self, x, **params):
        _, b = self._support(**params)
        return self._quadrature(self._pxf_dispatch, limits=(x, b),
                                params=params)

    ## Inverse cumulative distribution functions

    @_set_invalid_nan
    def ilogcdf(self, logp, /, *, method=None):
        return self._ilogcdf_dispatch(logp, method=method, **self._parameters)

    @_dispatch
    def _ilogcdf_dispatch(self, x, method=None, **params):
        if self._overrides('_ilogcdf_formula'):
            method = self._ilogcdf_formula
        elif self._overrides('_ilogccdf_formula'):
            method = self._ilogcdf_complement
        else:
            method = self._ilogcdf_inversion
        return method

    def _ilogcdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _ilogcdf_complement(self, x, **params):
        return self._ilogccdf_dispatch(_log1mexp(x), **params)

    def _ilogcdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._logcdf_dispatch, x, params=params)

    @_set_invalid_nan
    def icdf(self, p, /, *, method=None):
        return self._icdf_dispatch(p, method=method, **self._parameters)

    @_dispatch
    def _icdf_dispatch(self, x, method=None, **params):
        if self._overrides('_icdf_formula'):
            method = self._icdf_formula
        elif self._overrides('_iccdf_formula'):
            method = self._icdf_complement_safe
        else:
            method = self._icdf_inversion
        return method

    def _icdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _icdf_complement(self, x, **params):
        return self._iccdf_dispatch(1 - x, **params)

    def _icdf_complement_safe(self, x, **params):
        xp = self._xp
        out = self._icdf_complement(x, **params)
        eps = xp.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else eps**0.5
        with np.errstate(invalid='ignore'):
            mask = tol * x < xp_spacing(1 - x)
        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._icdf_inversion(x[mask], **params_mask)
        return out[()]

    def _icdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._cdf_dispatch, x, params=params)

    @_set_invalid_nan
    def ilogccdf(self, logp, /, *, method=None):
        return self._ilogccdf_dispatch(logp, method=method, **self._parameters)

    @_dispatch
    def _ilogccdf_dispatch(self, x, method=None, **params):
        if self._overrides('_ilogccdf_formula'):
            method = self._ilogccdf_formula
        elif self._overrides('_ilogcdf_formula'):
            method = self._ilogccdf_complement
        else:
            method = self._ilogccdf_inversion
        return method

    def _ilogccdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _ilogccdf_complement(self, x, **params):
        return self._ilogcdf_dispatch(_log1mexp(x), **params)

    def _ilogccdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._logccdf_dispatch, x, params=params)

    @_set_invalid_nan
    def iccdf(self, p, /, *, method=None):
        return self._iccdf_dispatch(p, method=method, **self._parameters)

    @_dispatch
    def _iccdf_dispatch(self, x, method=None, **params):
        if self._overrides('_iccdf_formula'):
            method = self._iccdf_formula
        elif self._overrides('_icdf_formula'):
            method = self._iccdf_complement_safe
        else:
            method = self._iccdf_inversion
        return method

    def _iccdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _iccdf_complement(self, x, **params):
        return self._icdf_dispatch(1 - x, **params)

    def _iccdf_complement_safe(self, x, **params):
        xp = self._xp
        out = self._iccdf_complement(x, **params)
        eps = xp.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else eps**0.5
        with np.errstate(invalid='ignore'):
            mask = tol * x < xp_spacing(1 - x)
        if xp.any(mask):
            params_mask = {key: xp.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = xp.asarray(out)
            out[mask] = self._iccdf_inversion(x[mask], **params_mask)
        return out[()]

    def _iccdf_inversion(self, x, **params):
        return self._solve_bounded_continuous(self._ccdf_dispatch, x, params=params)

    ### Sampling Functions
    # The following functions for drawing samples from the distribution are
    # exposed via a public method that accepts one positional argument - the
    # shape of the sample - and keyword options (but not distribution
    # parameters).
    # sample
    # ~~qmc_sample~~ built into sample now
    #
    # Common keyword options include:
    # method - a string that indicates which method should be used to compute
    #          the quantity (e.g. a formula or numerical integration).
    # rng - the NumPy Generator/SciPy QMCEngine object to used for drawing numbers.
    #
    # Input/output validation is included in each function, since there is
    # little code to be shared.
    # These are the methods meant to be called by users.
    #
    # Each public method calls a private "dispatch" method that
    # determines which "method" (strategy for calculating the desired quantity)
    # to use by default and, via the `@_dispatch` decorator, calls the
    # method and computes the result.
    # Each dispatch method can designate the responsibility of sampling to any
    # of several "implementation" methods. These methods accept only
    # `**params`, the parameter dictionary passed from the public method via
    # the "dispatch" method.
    # See the note corresponding with the "Distribution Parameters" for more
    # information.

    def sample(self, shape=(), *, method=None, rng=None):
        # needs output validation to ensure that developer returns correct
        # dtype and shape
        xp = self._xp
        sample_shape = (shape,) if not np.iterable(shape) else tuple(shape)
        full_shape = sample_shape + self._shape
        if isinstance(rng, qmc.QMCEngine):
            pass
        # elif is_cupy(xp):
        #     rng = xp.random.default_rng(rng)
        else:
            rng = np.random.default_rng(rng)
        res = self._sample_dispatch(full_shape, method=method, rng=rng,
                                    **self._parameters)

        return xp.astype(res, self._dtype, copy=False)

    @_dispatch
    def _sample_dispatch(self, full_shape, *, method, rng, **params):
        # make sure that tests catch if sample is 0d array
        if self._overrides('_sample_formula') and not isinstance(rng, qmc.QMCEngine):
            method = self._sample_formula
        else:
            method = self._sample_inverse_transform
        return method

    def _sample_formula(self, full_shape, *, rng, **params):
        raise NotImplementedError(self._not_implemented)

    def _sample_inverse_transform(self, full_shape, *, rng, **params):
        xp = self._xp
        if isinstance(rng, qmc.QMCEngine):
            uniform = self._qmc_uniform(full_shape, qrng=rng, **params)
        else:
            dtype = np.float32 if 'float32' in str(self._dtype) else np.float64
            uniform = rng.random(size=full_shape, dtype=dtype)
        uniform = xp.asarray(uniform, dtype=self._dtype, device=self._device)
        return self._icdf_dispatch(uniform, **params)

    def _qmc_uniform(self, full_shape, *, qrng, **params):
        # Generate QMC uniform sample(s) on unit interval with specified shape;
        # if `sample_shape != ()`, then each slice along axis 0 is independent.
        xp = self._xp

        sample_shape = full_shape[:len(full_shape)-len(self._shape)]
        # Determine the number of independent sequences and the length of each.
        n_low_discrepancy = sample_shape[0] if sample_shape else 1
        n_independent = math.prod(full_shape[1:] if sample_shape else full_shape)

        # For each independent sequence, we'll need a new QRNG of the appropriate class
        # with its own RNG. (If scramble=False, we don't really need all the separate
        # rngs, but I'm not going to add a special code path right now.)
        rngs = _rng_spawn(qrng.rng, n_independent)
        qrng_class = qrng.__class__
        kwargs = dict(d=1, scramble=qrng.scramble, optimization=qrng._optimization)
        if isinstance(qrng, qmc.Sobol):
            kwargs['bits'] = qrng.bits

        # Draw uniform low-discrepancy sequences scrambled with each RNG
        uniforms = []
        for rng in rngs:
            qrng = qrng_class(rng=rng, **kwargs)
            uniform = qrng.random(n_low_discrepancy)
            uniform = uniform.reshape(n_low_discrepancy if sample_shape else ())[()]
            uniforms.append(xp.asarray(uniform))

        # Reorder the axes and ensure that the shape is correct
        empty = xp.asarray([], device = self._device)
        uniform = xp.moveaxis(xp.stack(uniforms), -1, 0) if uniforms else empty
        return xp.reshape(uniform, full_shape)

    ### Moments
    # The `moment` method accepts two positional arguments - the order and kind
    # (raw, central, or standard) of the moment - and a keyword option:
    # method - a string that indicates which method should be used to compute
    #          the quantity (e.g. a formula or numerical integration).
    # Like the distribution properties, input/output validation is provided by
    # the `_set_invalid_nan_property` decorator.
    #
    # Unlike most public methods above, `moment` dispatches to one of three
    # private methods - one for each 'kind'. Like most *public* methods above,
    # each of these private methods calls a private "dispatch" method that
    # determines which "method" (strategy for calculating the desired quantity)
    # to use. Also, each dispatch method can designate the responsibility
    # computing the moment to one of several "implementation" methods.
    # Unlike the dispatch methods above, however, the `@_dispatch` decorator
    # is not used, and both logic and method calls are included in the function
    # itself.
    # Instead of determining which method will be used based solely on the
    # implementation methods available and calling only the corresponding
    # implementation method, *all* the implementation methods are called
    # in sequence until one returns the desired information. When an
    # implementation methods cannot provide the requested information, it
    # returns the object None (which is distinct from arrays with NaNs or infs,
    # which are valid values of moments).
    # The reason for this approach is that although formulae for the first
    # few moments of a distribution may be found, general formulae that work
    # for all orders are not always easy to find. This approach allows the
    # developer to write "formula" implementation functions that return the
    # desired moment when it is available and None otherwise.
    #
    # Note that the first implementation method called is a cache. This is
    # important because lower-order moments are often needed to compute
    # higher moments from formulae, so we eliminate redundant calculations
    # when moments of several orders are needed.

    @cached_property
    def _moment_methods(self):
        return {'cache', 'formula', 'transform',
                'normalize', 'general', 'quadrature'}

    @cached_property
    def _lmoment_methods(self):
        return {'cache', 'formula', 'general', 'order_statistics', 'quadrature_icdf'}

    @property
    def _zero(self):
        return self._constants()[0]

    @property
    def _one(self):
        return self._constants()[1]

    def _constants(self):
        if self._constant_cache is not None:
            return self._constant_cache

        constants = self._preserve_type([0, 1])

        if self.cache_policy != _NO_CACHE:
            self._constant_cache = constants

        return constants

    @_set_invalid_nan_property
    def moment(self, order=1, kind='raw', *, method=None):
        kinds = {'raw': self._moment_raw,
                 'central': self._moment_central}
        if not isinstance(self, CircularDistribution):
            kinds['standardized'] = self._moment_standardized
        order = self._validate_order(order)
        kind = self._validate_kind(kind, kinds)
        moment_kind = kinds[kind]
        return moment_kind(order, method=method)

    def _moment_raw(self, order=1, *, method=None):
        """Raw distribution moment about the origin."""
        # Consider exposing the point about which moments are taken as an
        # option. This is easy to support, since `_moment_transform_center`
        # does all the work.
        methods = self._moment_methods if method is None else {method}
        return self._moment_raw_dispatch(order, methods=methods, **self._parameters)

    def _moment_raw_dispatch(self, order, *, methods, **params):
        moment = None

        if 'cache' in methods:
            moment = self._moment_raw_cache.get(order, None)

        if moment is None and 'formula' in methods:
            moment = self._moment_raw_formula(order, **params)

        if moment is None and 'transform' in methods:
            moment = self._moment_raw_transform(order, **params)

        if moment is None and 'general' in methods:
            moment = self._moment_raw_general(order, **params)

        if moment is None and 'quadrature' in methods:
            moment = self._moment_from_pxf(order, center=self._zero, **params)

        if moment is None and 'quadrature_icdf' in methods:
            moment = self._moment_integrate_icdf(order, center=self._zero, **params)

        if moment is not None and self.cache_policy != _NO_CACHE:
            self._moment_raw_cache[order] = moment

        return moment

    def _moment_raw_formula(self, order, **params):
        return None

    def _moment_raw_transform(self, order, **params):
        xp = self._xp
        central_moments = []
        for i in range(int(order) + 1):
            methods = {'cache', 'formula', 'normalize', 'general'}
            moment_i = self._moment_central_dispatch(order=i, methods=methods, **params)
            if moment_i is None:
                return None
            moment_i = xp.asarray(moment_i, dtype=self._dtype, device=self._device)
            central_moments.append(moment_i)

        # Doesn't make sense to get the mean by "transform", since that's
        # how we got here. Questionable whether 'quadrature' should be here.
        mean_methods = {'cache', 'formula', 'quadrature'}
        mean = self._moment_raw_dispatch(1, methods=mean_methods, **params)
        if mean is None:
            return None

        moment = self._moment_transform_center(order, central_moments, mean, self._zero)
        return moment

    def _moment_raw_general(self, order, **params):
        # This is the only general formula for a raw moment of a probability
        # distribution
        return self._one if order == 0 else None

    def _moment_central(self, order=1, *, method=None):
        """Distribution moment about the mean."""
        methods = self._moment_methods if method is None else {method}
        return self._moment_central_dispatch(order, methods=methods, **self._parameters)

    def _moment_central_dispatch(self, order, *, methods, **params):
        xp = self._xp
        moment = None

        if 'cache' in methods:
            moment = self._moment_central_cache.get(order, None)

        if moment is None and 'formula' in methods:
            moment = self._moment_central_formula(order, **params)

        if moment is None and 'transform' in methods:
            moment = self._moment_central_transform(order, **params)

        if moment is None and 'normalize' in methods and order > 2:
            moment = self._moment_central_normalize(order, **params)

        if moment is None and 'general' in methods:
            moment = self._moment_central_general(order, **params)

        if moment is None and 'quadrature' in methods:
            mean = self._moment_raw_dispatch(1, **params, methods=self._moment_methods)
            if isinstance(self, CircularDistribution):
                mean = xp.angle(mean)
            moment = self._moment_from_pxf(order, center=mean, **params)

        if moment is None and 'quadrature_icdf' in methods:
            mean = self._moment_raw_dispatch(1, **params, methods=self._moment_methods)
            moment = self._moment_integrate_icdf(order, center=mean, **params)

        if moment is not None and self.cache_policy != _NO_CACHE:
            self._moment_central_cache[order] = moment

        return moment

    def _moment_central_formula(self, order, **params):
        return None

    def _moment_central_transform(self, order, **params):
        xp = self._xp
        raw_moments = []
        for i in range(int(order) + 1):
            methods = {'cache', 'formula', 'general'}
            moment_i = self._moment_raw_dispatch(order=i, methods=methods, **params)
            if moment_i is None:
                return None
            moment_i = xp.asarray(moment_i, dtype=self._dtype, device=self._device)
            raw_moments.append(moment_i)

        mean_methods = self._moment_methods
        mean = self._moment_raw_dispatch(1, methods=mean_methods, **params)

        moment = self._moment_transform_center(order, raw_moments, self._zero, mean)
        return moment

    def _moment_central_normalize(self, order, **params):
        methods = {'cache', 'formula', 'general'}
        standard_moment = self._moment_standardized_dispatch(order, **params,
                                                             methods=methods)
        if standard_moment is None:
            return None
        var = self._moment_central_dispatch(2, methods=self._moment_methods, **params)
        return standard_moment*var**(order/2)

    def _moment_central_general(self, order, **params):
        general_central_moments = {0: self._one, 1: self._zero}
        return general_central_moments.get(order, None)

    def _moment_standardized(self, order=1, *, method=None):
        """Standardized distribution moment."""
        methods = self._moment_methods if method is None else {method}
        return self._moment_standardized_dispatch(order, methods=methods,
                                                  **self._parameters)

    def _moment_standardized_dispatch(self, order, *, methods, **params):
        moment = None

        if 'cache' in methods:
            moment = self._moment_standardized_cache.get(order, None)

        if moment is None and 'formula' in methods:
            moment = self._moment_standardized_formula(order, **params)

        if moment is None and 'normalize' in methods:
            moment = self._moment_standardized_normalize(order, False, **params)

        if moment is None and 'general' in methods:
            moment = self._moment_standardized_general(order, **params)

        if moment is None and 'normalize' in methods:
            moment = self._moment_standardized_normalize(order, True, **params)

        if moment is not None and self.cache_policy != _NO_CACHE:
            self._moment_standardized_cache[order] = moment

        return moment

    def _moment_standardized_formula(self, order, **params):
        return None

    def _moment_standardized_normalize(self, order, use_quadrature, **params):
        methods = ({'quadrature'} if use_quadrature
                   else {'cache', 'formula', 'transform'})
        central_moment = self._moment_central_dispatch(order, **params,
                                                       methods=methods)
        if central_moment is None:
            return None
        var = self._moment_central_dispatch(2, methods=self._moment_methods,
                                            **params)
        return central_moment/var**(order/2)

    def _moment_standardized_general(self, order, **params):
        general_standard_moments = {0: self._one, 1: self._zero, 2: self._one}
        return general_standard_moments.get(order, None)

    def _moment_from_pxf(self, order, center, **params):
        def integrand(x, order, center, **params):
            pxf = self._pxf_dispatch(x, **params)
            return pxf*(x-center)**order
        return self._quadrature(integrand, args=(order, center), params=params)

    def _moment_integrate_icdf(self, order, center, **params):
        def integrand(x, order, center, **params):
            x = self._icdf_dispatch(x, **params)
            return (x-center)**order
        return self._quadrature(integrand, limits=(0., 1.),
                                args=(order, center), params=params)

    def _moment_transform_center(self, order, moment_as, a, b):
        xp = self._xp
        a, b = self._preserve_type(a), self._preserve_type(b)
        a, b, *moment_as = xp.broadcast_arrays(a, b, *moment_as)
        n = order
        # orthogonal to other axes
        i = xp.arange(n+1, device=self._device).reshape([-1]+[1]*a.ndim)
        i = self._preserve_type(i)
        n_choose_i = special.binom(n, i)
        with np.errstate(invalid='ignore'):  # can happen with infinite moment
            moment_as = self._preserve_type(xp.stack(tuple(moment_as)))
            moment_b = xp.sum(n_choose_i*moment_as*(a-b)**(n-i), axis=0)
        return moment_b

    ### L-Moments

    @_set_invalid_nan_property
    def lmoment(self, order=1, *, standardize=True, method=None):
        order = self._validate_order(order, fname='lmoment', min_order=1)
        return self._lmoment(order, standardize=standardize, method=method)

    def _lmoment(self, order, *, standardize, method):
        methods = self._lmoment_methods if method is None else {method}

        lmoment = self._lmoment_dispatch(order, methods=methods, **self._parameters)
        if lmoment is None:
            return None

        if standardize and order >= 3 and lmoment is not None:
            lscale = self._lmoment_dispatch(2, methods=methods, **self._parameters)
            if lscale is None:
                return None
            lmoment = lmoment / lscale

        return lmoment

    def _lmoment_dispatch(self, order, *, methods, **params):
        lmoment = None

        if 'cache' in methods:
            lmoment = self._lmoment_cache.get(order, None)

        if lmoment is None and 'formula' in methods:
            lmoment = self._lmoment_formula(order, **params)

        if lmoment is None and 'general' in methods:
            lmoment = self._lmoment_general(order, **params)

        if lmoment is None and 'quadrature_icdf' in methods and (
                self._overrides('_icdf_formula') or self._overrides('_iccdf_formula')):
            lmoment = self._lmoment_integrate_icdf(order, **params)

        if lmoment is None and 'order_statistics' in methods:
            lmoment = self._lmoment_from_order_statistics(order, **params)

        if lmoment is None and 'quadrature_icdf' in methods:
            lmoment = self._lmoment_integrate_icdf(order, **params)

        if lmoment is not None and self.cache_policy != _NO_CACHE:
            self._lmoment_cache[order] = lmoment

        return lmoment

    def _lmoment_formula(self, order, **params):
        return None

    def _lmoment_general(self, order, **params):
        return self.mean() if order == 1 else None

    def _lmoment_from_order_statistics(self, order, **params):
        from scipy.stats._univariate_distribution._extensions import order_statistic
        xp = self._xp
        k = xp.arange(order, dtype=self._dtype, device=self._device)
        k = xp.moveaxis(xpx.atleast_nd(k, ndim=self._ndim + 1), -1, 0)
        E = order_statistic(self, r=order-k, n=order).mean()
        bc = xp.astype(special.binom(order-1, k), self._dtype)
        return xp.sum((-1)**k * bc * E, axis=0) / order

    def _lmoment_integrate_icdf(self, order, **params):
        xp = self._xp
        def integrand(p, **params):
            x = self._icdf_dispatch(p, **params)
            P = special.eval_sh_legendre(order - 1, p)  # NOT array API compatible
            P = xp.astype(P, x.dtype)
            return x * P
        return self._quadrature(integrand, limits=(0., 1.), params=params)

    ### Convenience

    def plot(self, x='x', y=None, *, t=None, ax=None):
        r"""Plot a function of the distribution.

        Convenience function for quick visualization of the distribution
        underlying the random variable.

        Parameters
        ----------
        x, y : str, optional
            String indicating the quantities to be used as the abscissa and
            ordinate (horizontal and vertical coordinates), respectively.
            Defaults are ``'x'`` (the domain of the random variable) and either
            ``'pdf'`` (the probability density function) (continuous) or
            ``'pdf'`` (the probability density function) (discrete).
            Valid values are:
            'x', 'pdf', 'pmf', 'cdf', 'ccdf', 'icdf', 'iccdf', 'logpdf', 'logpmf',
            'logcdf', 'logccdf', 'ilogcdf', 'ilogccdf'.
        t : 3-tuple of (str, float, float), optional
            Tuple indicating the limits within which the quantities are plotted.
            The default is ``('cdf', 0.0005, 0.9995)`` if the domain is infinite,
            indicating that the central 99.9% of the distribution is to be shown;
            otherwise, endpoints of the support are used where they are finite.
            Valid values are:
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

        # NumPy-only, for now?

        discrete = isinstance(self, DiscreteDistribution)
        t_is_quantile = {'x', 'icdf', 'iccdf', 'ilogcdf', 'ilogccdf'}
        t_is_probability = {'cdf', 'ccdf', 'logcdf', 'logccdf'}
        valid_t = t_is_quantile.union(t_is_probability)
        valid_xy =  valid_t.union({'pdf', 'logpdf', 'pmf', 'logpmf'})
        y_default = 'pmf' if discrete else 'pdf'
        y = y_default if y is None else y

        ndim = self._ndim
        x_name, y_name = x, y
        t_name = 'cdf' if t is None else t[0]

        a, b = self.support()
        tliml_default = 0 if np.all(np.isfinite(a)) else 0.0005
        tliml = tliml_default if t is None else t[1]
        tlimr_default = 1 if np.all(np.isfinite(b)) else 0.9995
        tlimr = tlimr_default if t is None else t[2]
        tlim = np.asarray([tliml, tlimr])
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
        if y_name not in valid_xy:
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
            import matplotlib.pyplot as plt
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
        if discrete and x_name in t_is_quantile:
            # should probably aggregate for large ranges
            q = np.arange(np.min(qlim[0]), np.max(qlim[1]) + 1)
            q = q[:, np.newaxis] if ndim else q
        else:
            grid = np.linspace(0, 1, 300)
            grid = grid[:, np.newaxis] if ndim else grid
            q = qlim[0] + (qlim[1] - qlim[0]) * grid
            q = np.round(q) if discrete else q

        # compute requested x and y at quantile grid
        x = q if x_name in t_is_quantile else getattr(self, x_name)(q)
        y = q if y_name in t_is_quantile else getattr(self, y_name)(q)

        # make plot
        x, y = np.broadcast_arrays(x.T, np.atleast_2d(y.T))
        for xi, yi in zip(x, y):  # plot is vectorized, but bar/step don't seem to be
            if discrete and x_name in t_is_quantile and y_name == 'pmf':
                # should this just be a step plot, too?
                ax.bar(xi, yi, alpha=np.sqrt(1/y.shape[0]))  # alpha heuristic
            elif discrete and x_name in t_is_quantile:
                values = yi
                edges = np.concatenate((xi, [xi[-1]+1]))
                ax.stairs(values, edges, baseline=None)
            else:
                ax.plot(xi, yi)
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


    ### Fitting
    # All methods above treat the distribution parameters as fixed, and the
    # variable argument may be a quantile or probability. The fitting functions
    # are fundamentally different because the quantiles (often observations)
    # are considered to be fixed, and the distribution parameters are the
    # variables. In a sense, they are like an inverse of the sampling
    # functions.
    #
    # At first glance, it would seem ideal for `fit` to be a classmethod,
    # called like `LogUniform.fit(sample=sample)`.
    # I tried this. I insisted on it for a while. But if `fit` is a
    # classmethod, it cannot call instance methods. If we want to support MLE,
    # MPS, MoM, MoLM, then we end up with most of the distribution functions
    # above needing to be classmethods, too. All state information, such as
    # tolerances and the underlying distribution of `ShiftedScaledDistribution`
    # and `OrderStatisticDistribution`, would need to be passed into all
    # methods. And I'm not really sure how we would call `fit` as a
    # classmethod of a transformed distribution - maybe
    # ShiftedScaledDistribution.fit would accept the class of the
    # shifted/scaled distribution as an argument?
    #
    # In any case, it was a conscious decision for the infrastructure to
    # treat the parameters as "fixed" and the quantile/percentile arguments
    # as "variable". There are a lot of advantages to this structure, and I
    # don't think the fact that a few methods reverse the fixed and variable
    # quantities should make us question that choice. It can still accommodate
    # these methods reasonably efficiently.


class ContinuousDistribution(UnivariateDistribution):
    def _overrides(self, method_name):
        if method_name in {'_logpmf_formula', '_pmf_formula'}:
            return True
        return super()._overrides(method_name)

    def _pmf_formula(self, x, **params):
        # using idiom from DiscreteDistribution `_pdf_formula` for now;
        # this will be removed if gh-24582 is merged with item #1 included,
        # and will be replaced with something more concise/readable otherwise.
        xp = self._xp
        if params:
            p = next(iter(params.values()))
            nan_result = xp.isnan(x) | xp.isnan(p)
        else:
            nan_result = xp.isnan(x)
        return xp.where(nan_result, xp.nan, 0.)

    def _logpmf_formula(self, x, **params):
        xp = self._xp
        if params:
            p = next(iter(params.values()))
            nan_result = xp.isnan(x) | xp.isnan(p)
        else:
            nan_result = xp.isnan(x)
        return xp.where(nan_result, xp.nan, -xp.inf)

    def _pxf_dispatch(self, x, *, method=None, **params):
        return self._pdf_dispatch(x, method=method, **params)

    def _logpxf_dispatch(self, x, *, method=None, **params):
        return self._logpdf_dispatch(x, method=method, **params)

    def _solve_bounded_continuous(self, func, p, params, bounds=None, xatol=None):
        return self._solve_bounded(func, p, params=params, bounds=bounds, xatol=xatol).x


class DiscreteDistribution(UnivariateDistribution):
    def _overrides(self, method_name):
        if method_name in {'_logpdf_formula', '_pdf_formula'}:
            return True
        return super()._overrides(method_name)

    def _logpdf_formula(self, x, **params):
        xp = self._xp
        if params:
            p = next(iter(params.values()))
            nan_result = xp.isnan(x) | xp.isnan(p)
        else:
            nan_result = xp.isnan(x)
        return xp.where(nan_result, xp.nan, xp.inf)

    def _pdf_formula(self, x, **params):
        xp = self._xp
        if params:
            p = next(iter(params.values()))
            nan_result = xp.isnan(x) | xp.isnan(p)
        else:
            nan_result = xp.isnan(x)
        return xp.where(nan_result, xp.nan, xp.inf)

    def _pxf_dispatch(self, x, *, method=None, **params):
        return self._pmf_dispatch(x, method=method, **params)

    def _logpxf_dispatch(self, x, *, method=None, **params):
        return self._logpmf_dispatch(x, method=method, **params)

    def _cdf_quadrature(self, x, **params):
        xp = self._xp
        return super()._cdf_quadrature(xp.floor(x), **params)

    def _logcdf_quadrature(self, x, **params):
        xp = self._xp
        return super()._logcdf_quadrature(xp.floor(x), **params)

    def _ccdf_quadrature(self, x, **params):
        xp = self._xp
        return super()._ccdf_quadrature(xp.floor(x + 1), **params)

    def _logccdf_quadrature(self, x, **params):
        xp = self._xp
        return super()._logccdf_quadrature(xp.floor(x + 1), **params)

    def _cdf2(self, x, y, *, method):
        raise NotImplementedError(
            "Two argument cdf functions are currently only supported for "
            "continuous distributions.")

    def _ccdf2(self, x, y, *, method):
        raise NotImplementedError(
            "Two argument cdf functions are currently only supported for "
            "continuous distributions.")

    def _logcdf2(self, x, y, *, method):
        raise NotImplementedError(
            "Two argument cdf functions are currently only supported for "
            "continuous distributions.")

    def _logccdf2(self, x, y, *, method):
        raise NotImplementedError(
            "Two argument cdf functions are currently only supported for "
            "continuous distributions.")

    def _lmoment(self, order, *, standardize, method):
        raise NotImplementedError(
            "L-moments are currently available only "
            "for continuous distributions.")

    def _solve_bounded_discrete(self, func, p, params, comp):
        xp = self._xp
        # We're trying to solve one of these two problems:
        # a) find the smallest integer x* within the support s.t. F(x*) >= p
        # b) find the smallest integer x* within the support s.t. G(x*) = 1 - F(x*) <= p
        # Our approach is to solve a continuous version of the problem that narrows the
        # solution down to an integer x s.t. either x* = x or x* = x + 1. At the end,
        # we'll choose between them.

        # First, solve func(x) == p where func is a continuous, monotone interpolant
        # of either the monotone increasing F or monotone decreasing G.
        res = self._solve_bounded(func, p, params=params, xatol=0.9)
        # Here, `_solve_bounded` can terminate for one of three reasons:
        # 1. `func(res.x) == p` (`fatol = 0` is satisfied),
        # 2. `res.xl` and `res.xr` bracket the root and `|res.xr - res.xl| <= xatol`, or
        # 3. There is no solution within the support.
        # There are several possible strategies for using `res.xl`, `res.x`, and/or
        # `res.xr` to find a solution to the original, discrete problem. Here is ours.

        # Consider case 2a. Because F is an increasing function, we know
        # that F(xr) >= p (and F(xl) <= p), so F(floor(xr) + 1) >= p.
        # F(floor(xr)) *may* be >= p, but we can't know until we evaluate it.
        # F(floor(xr) - 1) < p (strictly) because floor(xr) - 1 < xl and F decreases
        # monotonically as the argument decreases. So we choose x = floor(xr), and
        # later we'll choose between x* = x and x* = x + 1.
        x = xp.asarray(xp.floor(res.bracket[-1]))
        # This is also suitable for case 2b. Because G is a *decreasing* function, we
        # know that G(xr) <= p (and G(xl) >= p), so G(floor(xr) + 1) <= p.
        # G(floor(xr)) *may* be <= p, but we can't know until we evaluate it.
        # G(floor(xr) - 1) > p (strictly) because floor(xr) - 1 < xl and G increases
        # as the argument decreases. So we would still want to choose x = floor(xr), and
        # later we'll choose between x* = x and x* = x + 1.

        # Now we consider case 1a/b. In this case, `res.x` solved the equation
        # *exactly*, so the algorithm may have terminated before the bracket is tight
        # enough to rely on `res.xr`. If `res.x` happens to be integral, `res.x` is
        # the solution to the discrete problem, and floor(res.x) == res.x, so
        # floor(res.x) is the solution to the discrete problem. If not:
        # a) F(floor(res.x)) < p (strictly) and F(floor(res.x) + 1) > p (strictly). So
        #    floor(res.x) + 1 is the solution to the discrete problem.
        # b) G(floor(res.x)) > p (strictly) and G(floor(res.x) + 1) < p (strictly). So
        #    floor(res.x) + 1 is again the solution to the discrete problem.
        # Either way, we can choose x = res.x, and at the end we'll choose between
        # x* = x and x* = x + 1.
        mask = res.f_x == 0
        x[mask] = xp.floor(res.x[mask])

        # For case 3, let xmin be the left endpoint of the support, and note that in
        # general, F(xmin) > 0 and G(xmin) < 1. Therefore it is possible that:
        # a) F(x) > p for all x in the support (e.g. because p ~ 0)
        # a) G(x) < p for all x in the support (e.g. because p ~ 1)
        # In these cases, `_solve_bounded` would fail to find a root of the continuous
        # equation above, but the solution to the original, discrete problem is the left
        # endpoint of the support.
        # This case is handled before we get to this function; otherwise,
        # `_solve_bounded` may spin its wheels for a long time in vain.

        # Now, we choose between x* = x and x* = x + 1: if func(x) satisfies the
        # comparison `comp` (>= for cdf, <= for ccdf), the solution is x* = x;
        # otherwise the solution must be x* = x + 1.
        f = func(x, **params)
        x = xp.where(comp(f, p), x, x + 1.0)
        x[xp.isnan(f)] = xp.nan  # needed? why would func(x) be NaN within support?

        return x

    def _base_discrete_inversion(self, p, func, comp, /, **params):
        # For discrete distributions, icdf(p) is defined as the minimum integer x*
        # within the support such that F(x*) >= p; iccdf(p) is the minimum integer x*
        # within the support such that G(x*) <= p.

        # Identify where the solution is xmin.
        # (See rationale in `_solve_bounded_discrete`.)
        xp = self._xp
        xmin, xmax = self._support(**params)
        p, xmin, _ = xp.broadcast_arrays(p, xmin, xmax)
        mask = comp(func(xmin, **params), p)

        # Use `apply_where` to perform the inversion only when necessary.
        def f1(p, *args):
            return self._solve_bounded_discrete(
                func, p, params=dict(zip(params.keys(), args)), comp=comp)

        x = xpx.apply_where(~mask, (p, *params.values()), f1, fill_value=xmin, xp=xp)

        # x above may be a finite value even when p is NaN, so the returned value
        # should be NaN. We need to handle this as a special case.
        x[xp.isnan(p)] = xp.nan
        return x[()]

    def _icdf_inversion(self, x, **params):
        xp = self._xp
        return self._base_discrete_inversion(x, self._cdf_dispatch,
                                             xp.greater_equal, **params)

    def _ilogcdf_inversion(self, x, **params):
        xp = self._xp
        return self._base_discrete_inversion(x, self._logcdf_dispatch,
                                             xp.greater_equal, **params)

    def _iccdf_inversion(self, x, **params):
        xp = self._xp
        return self._base_discrete_inversion(x, self._ccdf_dispatch,
                                             xp.less_equal, **params)

    def _ilogccdf_inversion(self, x, **params):
        xp = self._xp
        return self._base_discrete_inversion(x, self._logccdf_dispatch,
                                             xp.less_equal, **params)

    def _mode_optimization(self, **params):
        # If `x` is the true mode of a unimodal continuous function, we can find
        # the mode among integers by rounding in each direction and checking
        # which is better. If the difference between `x` and the nearest integer
        # is less than `xatol`, the computed value of `x` may end up on the wrong
        # side of the nearest integer. Setting `xatol=0.5` guarantees that at most
        # three integers need to be checked, the two nearest integers, ``floor(x)``
        # and ``round(x)`` and the nearest integer other than these.
        xp = self._xp
        x = super()._mode_optimization(xatol=0.5, **params)
        low, high = self.support()
        xl, xr = xp.floor(x), xp.ceil(x)
        nearest = xp.round(x)
        # Clip to stay within support. There will be redundant calculation
        # when clipping since `xo` will be one of `xl` or `xr`, but let's
        # keep the implementation simple for now.
        xo = xp.clip(nearest + xp.copysign(1, nearest - x), low, high)
        x = xp.stack([xl, xo, xr])
        idx = xp.argmax(self._pmf_dispatch(x, **params), axis=0)
        return np.choose(idx, [xl, xo, xr])  # TODO: fix this

    def _logentropy_quadrature(self, **params):
        def logintegrand(x, **params):
            xp = self._xp
            logpmf = self._logpmf_dispatch(x, **params)
            # Entropy summand is -pmf*log(pmf), so log-entropy summand is
            # logpmf + log(logpmf) + pi*j. But pmf is always between 0 and 1,
            # so logpmf is always negative, and so log(logpmf) = log(-logpmf) + pi*j.
            # The two imaginary components "cancel" each other out (which we would
            # expect because each term of the entropy summand is positive).
            return xp.where(xp.isfinite(logpmf), logpmf + xp.log(-logpmf), -xp.inf)

        with np.errstate(invalid='ignore', divide='ignore'):
            # Zeros and infinities in the integrand are noisy. It is slow and cumbersome
            # to use apply_where to avoid the warnings, so silence them.
            return self._quadrature(logintegrand, params=params, log=True)


class CircularDistribution(UnivariateDistribution):

    @cached_property
    def _moment_methods(self):
        return {'cache', 'formula', 'transform', 'general', 'quadrature'}

    def _cdf2(self, x, y, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not "
                                  "support two-argument `cdf`.")

    def _ccdf2(self, x, y, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not "
                                  "support two-argument `ccdf`.")

    def _lmoment(self, order, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `lmoment`.")

    def _logcdf1(self, x, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `logcdf`.")

    def _logcdf2(self, x, y, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `logcdf`.")

    def _logccdf1(self, x, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `logccdf`.")

    def _logccdf2(self, x, y, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `logccdf`.")

    def _ilogcdf_dispatch(self, p, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `ilogcdf`.")

    def _ilogccdf_dispatch(self, x, *, method, **kwargs):
        raise NotImplementedError("Circular distributions do not support `ilogccdf`.")

    def _mode_optimization(self, xatol=None, **params):
        # super()._mode_optimization guess is the median, which can be expensive
        a, b = self._support(**params)
        m = (a + b) / 2
        def f(x, **params): return -self._pxf_dispatch(x, **params)
        return self._optimization(f, m, xatol, params)

    def _median_optimization(self, **params):

        def f(m, **params):

            def integrand1(x, m, **params):
                a, b = self._support(**params)
                x_wrapped = (x - a) % (b - a) + a
                return (x - m) * self._pdf_dispatch(x_wrapped, **params)

            def integrand2(x, m, **params):
                a, b = self._support(**params)
                x_wrapped = (x - a) % (b - a) + a
                return (m - x) * self._pdf_dispatch(x_wrapped, **params)

            a, b = self._support(**params)
            half_period = (b - a) / 2
            res1 = self._quadrature(integrand1, limits=(m, m + half_period),
                                    args=(m,), params=params)
            res2 = self._quadrature(integrand2, limits=(m - half_period, m),
                                    args=(m,), params=params)
            return res1 + res2

        a, b = self._support(**params)
        m0 = (a + b) / 2
        return self._optimization(f, m0, None, params)

    def _moment_from_pxf(self, order, center, **params):
        xp = self._xp
        def integrand(x, order, center, **params):
            xp = self._xp
            a, b = self._support(**params)
            period = b - a
            scale = 2*xp.pi / period
            # All internal moment calculations are in radians with the origin at the
            # left endpoint of the support, so we are integrating over [0, 2*pi].
            # Therefore, x needs to be scaled and shifted to the support of the
            # distribution. But the returned PDF values need to be scaled back
            # because we are integrating over [0, 2*pi], not the support. This is the
            # same as evaluation of a PDF under generic shift/scale of a distribution;
            # It just looks slightly different.
            pdf = self._pdf_dispatch(x/scale + a, **params) / scale
            return xp.exp(1j * order * (x - center)) * pdf
        return self._quadrature(integrand, limits=(0, 2*xp.pi),
                                args=(order, center), params=params)

    def mean(self, *, method=None):
        xp = self._xp
        a, b = self.support()
        period = b - a
        scale = 2*xp.pi / period
        phi = self.moment(1, kind='raw', method=method)
        return xp.angle(phi) / scale + a

    def variance(self, *, method=None):
        xp = self._xp
        rho = xp.real(self.moment(1, kind='central', method=method))
        return 1 - rho

    def standard_deviation(self, *, method=None):
        xp = self._xp
        rho = xp.real(self.moment(1, kind='central', method=method))
        return xp.sqrt(-2*xp.log(rho))

    def skewness(self, *, method=None):
        xp = self._xp
        b2 = xp.imag(self.moment(2, kind='central', method=method))
        rho = xp.real(self.moment(1, kind='central', method=method))
        return b2 / (1 - rho)**1.5

    def kurtosis(self, *, method=None, convention=None):
        xp = self._xp
        message = (f'`{self.__class__.__name__}.kurtosis` supports only the default '
                   f"value of `convention`.")
        if convention is not None:
            raise ValueError(message)

        # This is the most common definition
        a2 = xp.real(self.moment(2, kind='central', method=method))
        rho = xp.real(self.moment(1, kind='central', method=method))
        return (a2 - rho**4) / (1 - rho)**2

    def _moment_raw_general(self, order, **params):
        general_raw_moments = {0: self._one + 0j}
        return general_raw_moments.get(order, None)

    def _moment_central_general(self, order, **params):
        general_central_moments = {0: self._one + 0j}
        return general_central_moments.get(order, None)

    def _moment_central_transform(self, order, **params):
        xp = self._xp
        methods = {'cache', 'formula', 'general'}
        moment = self._moment_raw_dispatch(order=order, methods=methods, **params)
        phi1 = self._moment_raw_dispatch(1, methods=methods, **params)
        if moment is None or phi1 is None:
            return None
        mu = xp.angle(phi1)
        moment = self._moment_transform_center(order, moment, self._zero, mu)
        return moment

    def _moment_raw_transform(self, order, **params):
        xp = self._xp
        methods = {'cache', 'formula', 'general'}
        moment = self._moment_central_dispatch(order=order, methods=methods, **params)
        phi1 = self._moment_raw_dispatch(1, methods=methods, **params)
        if moment is None or phi1 is None:
            return None
        mu = xp.angle(phi1)
        moment = self._moment_transform_center(order, moment, mu, self._zero)
        return moment

    def _moment_transform_center(self, order, moment, a, b):
        xp = self._xp
        a, b = self._preserve_type(a), self._preserve_type(b)
        a, b, moment = xp.broadcast_arrays(a, b, moment)
        n = order
        moment_b = moment * xp.exp(1j*n*(a - b))
        return moment_b


@rv_capabilities
class TransformedDistribution(ContinuousDistribution):
    def __init__(self, X, /, *args, **kwargs):
        if (isinstance(X, CircularDistribution)
                or not isinstance(X, ContinuousDistribution)):
            message = ("Transformations are currently supported only for continuous, "
                       "linear RVs.")
            raise NotImplementedError(message)
        self._copy_parameterization()
        self._variable = X._variable
        self._dist = X
        self._xp = X._xp
        self._dtype = X._dtype
        self._device = X._device
        if X._parameterization:
            # Add standard distribution parameters to our parameterization
            dist_parameters = X._parameterization.parameters
            set_params = set(dist_parameters)
            if not self._parameterizations:
                self._parameterizations.append(_Parameterization())
            for parameterization in self._parameterizations:
                if set_params.intersection(parameterization.parameters):
                    message = (f"One or more of the parameters of {X} has "
                               "the same name as a parameter of "
                               f"{self.__class__.__name__}. Name collisions "
                               "create ambiguities and are not supported.")
                    raise ValueError(message)
                parameterization.parameters.update(dist_parameters)
        super().__init__(*args, **kwargs)

    def _overrides(self, method_name):
        return (self._dist._overrides(method_name)
                or super()._overrides(method_name))

    def reset_cache(self):
        super().reset_cache()

    def _update_parameters(self, *, validation_policy=None, **params):
        # maybe broadcast everything before processing?
        parameters = {}
        # There may be some issues with _original_parameters
        # We only want to update with _dist._original_parameters during
        # initialization. Afterward that, we want to start with
        # self._original_parameters.
        parameters.update(self._dist._original_parameters)
        parameters.update(params)
        super()._update_parameters(validation_policy=validation_policy, **parameters)

    def _process_parameters(self, **params):
        return self._dist._process_parameters(**params)

    def __repr__(self):
        raise NotImplementedError()

    def __str__(self):
        raise NotImplementedError()


@rv_capabilities
class ShiftedScaledDistribution(TransformedDistribution):
    """Distribution with a standard shift/scale transformation."""
    # Unclear whether infinite loc/scale will work reasonably in all cases
    _loc_domain = _RealInterval(endpoints=(-inf, inf), inclusive=(True, True))
    _loc_param = _RealParameter('loc', symbol=r'\mu',
                                domain=_loc_domain, typical=(1, 2))

    _scale_domain = _RealInterval(endpoints=(-inf, inf), inclusive=(True, True))
    _scale_param = _RealParameter('scale', symbol=r'\sigma',
                                  domain=_scale_domain, typical=(0.1, 10))

    _parameterizations = [_Parameterization(_loc_param, _scale_param),
                          _Parameterization(_loc_param),
                          _Parameterization(_scale_param)]

    def _process_parameters(self, loc=None, scale=None, **params):
        xp = self._xp
        loc = loc if loc is not None else xp.zeros_like(scale)[()]
        scale = scale if scale is not None else xp.ones_like(loc)[()]
        sign = scale > 0
        parameters = self._dist._process_parameters(**params)
        parameters.update(dict(loc=loc, scale=scale, sign=sign))
        return parameters

    def _transform(self, x, loc, scale, **kwargs):
        return (x - loc)/scale

    def _itransform(self, x, loc, scale, **kwargs):
        return x * scale + loc

    def _support(self, loc, scale, sign, **params):
        # Add shortcut for infinite support?
        xp = self._xp
        a, b = self._dist._support(**params)
        a, b = self._itransform(a, loc, scale), self._itransform(b, loc, scale)
        return xp.where(sign, a, b)[()], xp.where(sign, b, a)[()]

    def __repr__(self):
        xp = self._xp
        with np.printoptions(threshold=10):
            result =  f"{repr(self.scale)}*{repr(self._dist)}"
            if not self.loc.ndim and self.loc < 0:
                result += f" - {repr(-self.loc)}"
            elif (xp.any(self.loc != 0)
                  or xp.result_type(self.loc.dtype, self.scale.dtype)
                      != self.scale.dtype):
                # We don't want to hide a zero array loc if it can cause
                # a type promotion.
                result += f" + {repr(self.loc)}"
        return result

    def __str__(self):
        xp = self._xp
        with np.printoptions(threshold=10):
            result =  f"{str(self.scale)}*{str(self._dist)}"
            if not self.loc.ndim and self.loc < 0:
                result += f" - {str(-self.loc)}"
            elif (xp.any(self.loc != 0)
                  or xp.result_type(self.loc.dtype, self.scale.dtype)
                      != self.scale.dtype):
                # We don't want to hide a zero array loc if it can cause
                # a type promotion.
                result += f" + {str(self.loc)}"
        return result

    # Here, we override all the `_dispatch` methods rather than the public
    # methods or _formula methods. Why not the public methods?
    # If we were to override the public methods, then other
    # TransformedDistribution classes (which could transform a
    # ShiftedScaledDistribution) would need to call the public methods of
    # ShiftedScaledDistribution, which would run the input validation again.
    # Why not the _formula methods? For distributions that rely on the
    # default implementation of methods (e.g. `quadrature`, `inversion`),
    # the implementation would "see" the location and scale like other
    # distribution parameters, so they could affect the accuracy of the
    # calculations. I think it is cleaner if `loc` and `scale` do not affect
    # the underlying calculations at all.

    def _entropy_dispatch(self, *args, loc, scale, sign, **params):
        xp = self._xp
        return (self._dist._entropy_dispatch(*args, **params)
                + xp.log(xp.abs(scale)))

    def _logentropy_dispatch(self, *args, loc, scale, sign, **params):
        xp = self._xp
        lH0 = self._dist._logentropy_dispatch(*args, **params)
        lls = xp.log(xp.log(xp.abs(scale))+0j)
        return special.logsumexp(xp.stack(xp.broadcast_arrays(lH0, lls)), axis=0)

    def _median_dispatch(self, *, method, loc, scale, sign, **params):
        raw = self._dist._median_dispatch(method=method, **params)
        return self._itransform(raw, loc, scale)

    def _mode_dispatch(self, *, method, loc, scale, sign, **params):
        raw = self._dist._mode_dispatch(method=method, **params)
        return self._itransform(raw, loc, scale)

    def _logpdf_dispatch(self, x, *args, loc, scale, sign, **params):
        xp = self._xp
        x = self._transform(x, loc, scale)
        logpdf = self._dist._logpdf_dispatch(x, *args, **params)
        return logpdf - xp.log(xp.abs(scale))

    def _pdf_dispatch(self, x, *args, loc, scale, sign, **params):
        xp = self._xp
        x = self._transform(x, loc, scale)
        pdf = self._dist._pdf_dispatch(x, *args, **params)
        return pdf / xp.abs(scale)

    def _logpmf_dispatch(self, x, *args, loc, scale, sign, **params):
        xp = self._xp
        x = self._transform(x, loc, scale)
        logpmf = self._dist._logpmf_dispatch(x, *args, **params)
        return logpmf - xp.log(xp.abs(scale))

    def _pmf_dispatch(self, x, *args, loc, scale, sign, **params):
        xp = self._xp
        x = self._transform(x, loc, scale)
        pmf = self._dist._pmf_dispatch(x, *args, **params)
        return pmf / xp.abs(scale)

    def _logpxf_dispatch(self, x, *args, loc, scale, sign, **params):
        xp = self._xp
        x = self._transform(x, loc, scale)
        logpxf = self._dist._logpxf_dispatch(x, *args, **params)
        return logpxf - xp.log(xp.abs(scale))

    def _pxf_dispatch(self, x, *args, loc, scale, sign, **params):
        xp = self._xp
        x = self._transform(x, loc, scale)
        pxf = self._dist._pxf_dispatch(x, *args, **params)
        return pxf / xp.abs(scale)

    @_shift_scale_distribution_function
    def _logcdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_distribution_function
    def _cdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_distribution_function
    def _logccdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_distribution_function
    def _ccdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_distribution_function_2arg
    def _logcdf2_dispatch(self, x, y, *, method=None, **params):
        pass

    @_shift_scale_distribution_function_2arg
    def _cdf2_dispatch(self, x, y, *, method=None, **params):
        pass

    @_shift_scale_distribution_function_2arg
    def _logccdf2_dispatch(self, x, y, *, method=None, **params):
        pass

    @_shift_scale_distribution_function_2arg
    def _ccdf2_dispatch(self, x, y, *, method=None, **params):
        pass

    @_shift_scale_inverse_function
    def _ilogcdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_inverse_function
    def _icdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_inverse_function
    def _ilogccdf_dispatch(self, x, *, method=None, **params):
        pass

    @_shift_scale_inverse_function
    def _iccdf_dispatch(self, x, *, method=None, **params):
        pass

    def _moment_standardized_dispatch(self, order, *, loc, scale, sign, methods,
                                      **params):
        xp = self._xp
        res = (self._dist._moment_standardized_dispatch(
            order, methods=methods, **params))
        return None if res is None else res * xp.sign(scale)**order

    def _moment_central_dispatch(self, order, *, loc, scale, sign, methods,
                                 **params):
        res = (self._dist._moment_central_dispatch(
            order, methods=methods, **params))
        return None if res is None else res * scale**order

    def _moment_raw_dispatch(self, order, *, loc, scale, sign, methods,
                             **params):
        raw_moments = []
        methods_highest_order = methods
        for i in range(int(order) + 1):
            methods = (self._moment_methods if i < order
                       else methods_highest_order)
            raw = self._dist._moment_raw_dispatch(i, methods=methods, **params)
            if raw is None:
                return None
            moment_i = raw * scale**i
            raw_moments.append(moment_i)

        return self._moment_transform_center(
            order, raw_moments, loc, self._zero)

    def _lmoment_dispatch(self, order, *, loc, scale, sign, methods, **params):
        xp = self._xp
        res = self._dist._lmoment_dispatch(order, methods=methods, **params)
        if res is None:  # if a specific method is requested but not available
            return None
        res = res * xp.abs(scale) * xp.sign(scale)**order
        return res + loc if order == 1 else res

    def _sample_dispatch(self, full_shape, *,
                         rng, loc, scale, sign, method, **params):
        rvs = self._dist._sample_dispatch(full_shape, method=method, rng=rng, **params)
        return self._itransform(rvs, loc=loc, scale=scale, sign=sign, **params)

    def __add__(self, loc):
        return ShiftedScaledDistribution(self._dist, loc=self.loc + loc,
                                         scale=self.scale)

    def __sub__(self, loc):
        return ShiftedScaledDistribution(self._dist, loc=self.loc - loc,
                                         scale=self.scale)

    def __mul__(self, scale):
        return ShiftedScaledDistribution(self._dist,
                                         loc=self.loc * scale,
                                         scale=self.scale * scale)

    def __truediv__(self, scale):
        return ShiftedScaledDistribution(self._dist,
                                         loc=self.loc / scale,
                                         scale=self.scale / scale)


@rv_capabilities
class MonotonicTransformedDistribution(TransformedDistribution):
    r"""Distribution underlying a strictly monotonic function of a random variable

    Given a random variable :math:`X`; a strictly monotonic function
    :math:`g(u)`, its inverse :math:`h(u) = g^{-1}(u)`, and the derivative magnitude
    :math: `|h'(u)| = \left| \frac{dh(u)}{du} \right|`, define the distribution
    underlying the random variable :math:`Y = g(X)`.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable :math:`X`.
    g, h, dh : callable
        Elementwise functions representing the mathematical functions
        :math:`g(u)`, :math:`h(u)`, and :math:`|h'(u)|`
    logdh : callable, optional
        Elementwise function representing :math:`\log(h'(u))`.
        The default is ``lambda u: np.log(dh(u))``, but providing
        a custom implementation may avoid over/underflow.
    increasing : bool, optional
        Whether the function is strictly increasing (True, default)
        or strictly decreasing (False).
    repr_pattern : str, optional
        A string pattern for determining the __repr__. The __repr__
        for X will be substituted into the position where `***` appears.
        For example:
            ``"exp(***)"`` for the repr of an exponentially transformed
            distribution
        The default is ``f"{g.__name__}(***)"``.
    str_pattern : str, optional
        A string pattern for determining `__str__`. The `__str__`
        for X will be substituted into the position where `***` appears.
        For example:
            ``"exp(***)"`` for the repr of an exponentially transformed
            distribution
        The default is the value `repr_pattern` takes.
    """

    def __init__(self, X, /, *args, g, h, dh, logdh=None,
                 increasing=True, repr_pattern=None,
                 str_pattern=None, **kwargs):
        super().__init__(X, *args, **kwargs)
        xp = self._xp
        self._g = g
        self._h = h
        self._dh = dh
        self._logdh = (logdh if logdh is not None
                       else lambda u: xp.log(dh(u)))
        if increasing:
            self._xdf = self._dist._cdf_dispatch
            self._cxdf = self._dist._ccdf_dispatch
            self._ixdf = self._dist._icdf_dispatch
            self._icxdf = self._dist._iccdf_dispatch
            self._logxdf = self._dist._logcdf_dispatch
            self._logcxdf = self._dist._logccdf_dispatch
            self._ilogxdf = self._dist._ilogcdf_dispatch
            self._ilogcxdf = self._dist._ilogccdf_dispatch
        else:
            self._xdf = self._dist._ccdf_dispatch
            self._cxdf = self._dist._cdf_dispatch
            self._ixdf = self._dist._iccdf_dispatch
            self._icxdf = self._dist._icdf_dispatch
            self._logxdf = self._dist._logccdf_dispatch
            self._logcxdf = self._dist._logcdf_dispatch
            self._ilogxdf = self._dist._ilogccdf_dispatch
            self._ilogcxdf = self._dist._ilogcdf_dispatch
        self._increasing = increasing
        self._repr_pattern = repr_pattern or f"{g.__name__}(***)"
        self._str_pattern = str_pattern or self._repr_pattern

    def __repr__(self):
        with np.printoptions(threshold=10):
            return self._repr_pattern.replace("***", repr(self._dist))

    def __str__(self):
        with np.printoptions(threshold=10):
            return self._str_pattern.replace("***", str(self._dist))

    def _overrides(self, method_name):
        # Do not use the generic overrides of TransformedDistribution
        return False

    def _support(self, **params):
        xp = self._xp
        a, b = self._dist._support(**params)
        # For reciprocal transformation, we want this zero to become -inf
        b = xp.where(b==0, -xp.zeros_like(b), b)
        with np.errstate(divide='ignore'):
            if self._increasing:
                return self._g(a), self._g(b)
            else:
                return self._g(b), self._g(a)

    def _logpdf_dispatch(self, x, *args, **params):
        return self._dist._logpdf_dispatch(self._h(x), *args, **params) + self._logdh(x)

    def _pdf_dispatch(self, x, *args, **params):
        return self._dist._pdf_dispatch(self._h(x), *args, **params) * self._dh(x)

    def _logcdf_dispatch(self, x, *args, **params):
        return self._logxdf(self._h(x), *args, **params)

    def _cdf_dispatch(self, x, *args, **params):
        return self._xdf(self._h(x), *args, **params)

    def _logccdf_dispatch(self, x, *args, **params):
        return self._logcxdf(self._h(x), *args, **params)

    def _ccdf_dispatch(self, x, *args, **params):
        return self._cxdf(self._h(x), *args, **params)

    def _ilogcdf_dispatch(self, p, *args, **params):
        return self._g(self._ilogxdf(p, *args, **params))

    def _icdf_dispatch(self, p, *args, **params):
        return self._g(self._ixdf(p, *args, **params))

    def _ilogccdf_dispatch(self, p, *args, **params):
        return self._g(self._ilogcxdf(p, *args, **params))

    def _iccdf_dispatch(self, p, *args, **params):
        return self._g(self._icxdf(p, *args, **params))

    def _sample_dispatch(self, full_shape, *, method, rng, **params):
        rvs = self._dist._sample_dispatch(full_shape, method=method, rng=rng, **params)
        return self._g(rvs)


class FoldedDistribution(TransformedDistribution):
    r"""Distribution underlying the absolute value of a random variable

    Given a random variable :math:`X`; define the distribution
    underlying the random variable :math:`Y = |X|`.

    Parameters
    ----------
    X : `ContinuousDistribution`
        The random variable :math:`X`.

    Returns
    -------
    Y : `ContinuousDistribution`
        The random variable :math:`Y = |X|`

    """
    # Many enhancements are possible if distribution is symmetric. Start
    # with the general case; enhance later.

    def __init__(self, X, /, *args, **kwargs):
        super().__init__(X, *args, **kwargs, dtype=X._dtype)
        # I think we need to allow `_support` to define whether the endpoints
        # are inclusive or not. In the meantime, it's best to ensure that the lower
        # endpoint (typically 0 for folded distribution) is inclusive so PDF evaluates
        # correctly at that point.
        self._variable.domain.inclusive = (True, self._variable.domain.inclusive[1])

    def _overrides(self, method_name):
        # Do not use the generic overrides of TransformedDistribution
        return False

    def _support(self, **params):
        xp = self._xp
        a, b = self._dist._support(**params)
        a_, b_ = xp.abs(a), xp.abs(b)
        a_, b_ = xp.minimum(a_, b_), xp.maximum(a_, b_)
        i = (a < 0) & (b > 0)
        a_ = xp.asarray(a_)
        a_[i] = 0
        return a_[()], b_[()]

    def _logpdf_dispatch(self, x, *args, method=None, **params):
        xp = self._xp
        x = xp.abs(x)
        right = self._dist._logpdf_dispatch(x, *args, method=method, **params)
        left = self._dist._logpdf_dispatch(-x, *args, method=method, **params)
        left = xp.asarray(left)
        right = xp.asarray(right)
        a, b = self._dist._support(**params)
        left[-x < a] = -xp.inf
        right[x > b] = -xp.inf
        logpdfs = xp.stack([left, right])
        return special.logsumexp(logpdfs, axis=0)

    def _pdf_dispatch(self, x, *args, method=None, **params):
        xp = self._xp
        x = xp.abs(x)
        right = self._dist._pdf_dispatch(x, *args, method=method, **params)
        left = self._dist._pdf_dispatch(-x, *args, method=method, **params)
        left = xp.asarray(left)
        right = xp.asarray(right)
        a, b = self._dist._support(**params)
        left[-x < a] = 0
        right[x > b] = 0
        return left + right

    def _logcdf_dispatch(self, x, *args, method=None, **params):
        xp = self._xp
        x = xp.abs(x)
        a, b = self._dist._support(**params)
        xl = xp.maximum(-x, a)
        xr = xp.minimum(x, b)
        res = self._dist._logcdf2_dispatch(xl, xr, *args, method=method, **params)
        return xp.real(res)

    def _cdf_dispatch(self, x, *args, method=None, **params):
        xp = self._xp
        x = xp.abs(x)
        a, b = self._dist._support(**params)
        xl = xp.maximum(-x, a)
        xr = xp.minimum(x, b)
        return self._dist._cdf2_dispatch(xl, xr, *args, **params)

    def _logccdf_dispatch(self, x, *args, method=None, **params):
        xp = self._xp
        x = xp.abs(x)
        a, b = self._dist._support(**params)
        xl = xp.maximum(-x, a)
        xr = xp.minimum(x, b)
        res = self._dist._logccdf2_dispatch(xl, xr, *args, method=method, **params)
        return xp.real(res)

    def _ccdf_dispatch(self, x, *args, method=None, **params):
        xp = self._xp
        x = xp.abs(x)
        a, b = self._dist._support(**params)
        xl = xp.maximum(-x, a)
        xr = xp.minimum(x, b)
        return self._dist._ccdf2_dispatch(xl, xr, *args, method=method, **params)

    def _sample_dispatch(self, full_shape, *, method, rng, **params):
        xp = self._xp
        rvs = self._dist._sample_dispatch(full_shape, method=method, rng=rng, **params)
        return xp.abs(rvs)

    def __repr__(self):
        with np.printoptions(threshold=10):
            return f"abs({repr(self._dist)})"

    def __str__(self):
        with np.printoptions(threshold=10):
            return f"abs({str(self._dist)})"


# These should probably go back to `_util.py`` but `_combine_docs` needs
# `UnivariateDistribution`
def _combine_docs(dist_family, *, include_examples=True):
    fields = set(NumpyDocString.sections)
    fields.remove('index')
    if not include_examples:
        fields.remove('Examples')

    doc = ClassDoc(dist_family)
    superdoc = ClassDoc(UnivariateDistribution)
    for field in fields:
        if field in {"Methods", "Attributes"}:
            doc[field] = superdoc[field]
        elif field in {"Summary"}:
            pass
        elif field == "Extended Summary":
            doc[field].append(_generate_domain_support(dist_family))
        elif field == 'Examples':
            doc[field] = [_generate_example(dist_family)]
        elif field == 'Notes':
            # Otherwise the (indented) array API support table doesn't
            # have ReST-required blank line after it
            doc[field] += ['\n'] + superdoc[field]
        else:
            doc[field] += superdoc[field]
    return str(doc)


def _generate_example(dist_family):
    n_parameters = dist_family._num_parameters(0)
    shapes = [()] * n_parameters
    rng = np.random.default_rng(615681484984984)
    i = 0
    xp = xp_compat_namespace(np)
    dist = dist_family._draw(shapes, rng=rng, i_parameterization=i, xp=xp)

    rng = np.random.default_rng(2354873452)
    name = dist_family.__name__
    if n_parameters:
        parameter_names = list(dist._parameterizations[i].parameters)
        parameter_values = [round(getattr(dist, name), 2) for name in
                            parameter_names]
        name_values = [f"{name}={value}" for name, value in
                       zip(parameter_names, parameter_values)]
        instantiation = f"{name}({', '.join(name_values)})"
        attributes = ", ".join([f"X.{param}" for param in dist._parameters])
        X = dist_family(**dict(zip(parameter_names, parameter_values)))
    else:
        instantiation = f"{name}()"
        X = dist

    p = 0.32
    x = round(X.icdf(p), 2)
    y = round(X.icdf(2 * p), 2)

    example = f"""
    To use the distribution class, it must be instantiated using keyword
    parameters corresponding with one of the accepted parameterizations.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> from scipy.stats import {name}
    >>> X = {instantiation}

    For convenience, the ``plot`` method can be used to visualize the density
    and other functions of the distribution.

    >>> X.plot()
    >>> plt.show()

    The support of the underlying distribution is available using the ``support``
    method.

    >>> X.support()
    {X.support()}
    """

    if n_parameters:
        example += f"""
        The numerical values of parameters associated with all parameterizations
        are available as attributes.

        >>> {attributes}
        {tuple(X._parameters.values())}
        """

    example += f"""
    To evaluate the probability density/mass function of the underlying distribution
    at argument ``x={x}``:

    >>> x = {x}
    >>> X.pdf(x), X.pmf(x)
    {X.pdf(x), X.pmf(x)}
    """

    if not issubclass(dist_family, CircularDistribution):
        example += """
    The cumulative distribution function, its complement, and the logarithm
    of these functions are evaluated similarly.

    >>> np.allclose(np.exp(X.logccdf(x)), 1 - X.cdf(x))
    True
        """

    # When two-arg CDF is implemented for circular and discrete distributions, consider
    # removing the special-casing here.
    if (issubclass(dist_family, ContinuousDistribution)
            and not issubclass(dist_family, CircularDistribution)):
        example_continuous = f"""
    The inverse of these functions with respect to the argument ``x`` is also
    available.

    >>> logp = np.log(1 - X.ccdf(x))
    >>> np.allclose(X.ilogcdf(logp), x)
    True

    Note that distribution functions and their logarithms also have two-argument
    versions for working with the probability mass between two arguments. The
    result tends to be more accurate than the naive implementation because it avoids
    subtractive cancellation.

    >>> y = {y}
    >>> np.allclose(X.ccdf(x, y), 1 - (X.cdf(y) - X.cdf(x)))
    True
        """
        example += example_continuous

    example += f"""
    There are methods for computing measures of central tendency,
    dispersion, higher moments, and entropy.

    >>> X.mean(), X.median(), X.mode()
    {X.mean(), X.median(), X.mode()}

    >>> X.variance(), X.standard_deviation()
    {X.variance(), X.standard_deviation()}

    >>> X.skewness(), X.kurtosis()
    {X.skewness(), X.kurtosis()}
    """

    if not issubclass(dist_family, CircularDistribution):
        example += """
    >>> np.allclose(X.moment(order=6, kind='standardized'),
    ...             X.moment(order=6, kind='central') / X.variance()**3)
    True
        """

    # When logentropy is implemented for DiscreteDistribution, remove special-casing
    if issubclass(dist_family, ContinuousDistribution):
        example += """
    >>> np.allclose(np.exp(X.logentropy()), X.entropy())
    True
        """
    else:
        example += f"""
        >>> X.entropy()
        {X.entropy()}
        """

    example += f"""
    Pseudo-random samples can be drawn from
    the underlying distribution using ``sample``.

    >>> X.sample(shape=(4,))
    {repr(X.sample(shape=(4,)))}  # may vary
    """
    # remove the indentation due to use of block quote within function;
    # eliminate blank first line
    example = "\n".join([line.lstrip() for line in example.split("\n")][1:])
    return example
