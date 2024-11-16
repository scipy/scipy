from abc import ABC, abstractmethod
import functools

import numpy as np
from numpy import inf

from scipy._lib._docscrape import ClassDoc, NumpyDocString
from scipy._lib._util import _lazywhere
from scipy import special
from scipy.integrate._tanhsinh import _tanhsinh
from scipy.optimize._bracket import _bracket_root, _bracket_minimum
from scipy.optimize._chandrupatla import _chandrupatla, _chandrupatla_minimize

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

# Could add other policies for broadcasting and edge/out-of-bounds case handling
# For instance, when edge case handling is known not to be needed, it's much
# faster to turn it off, but it might still be nice to have array conversion
# and shaping done so the user doesn't need to be so careful.
_SKIP_ALL = "skip_all"
# Other cache policies would be useful, too.
_NO_CACHE = "no_cache"

#### Domains

class _Domain(ABC):
    r""" Representation of the applicable domain of a parameter or variable.

    A `_Domain` object is responsible for storing information about the
    domain of a parameter or variable, determining whether a value is within
    the domain (`contains`), and providing a text/mathematical representation
    of itself (`__str__`). Because the domain of a parameter/variable can have
    a complicated relationship with other parameters and variables of a
    distribution, `_Domain` itself does not try to represent all possibilities;
    in fact, it has no implementation and is meant for subclassing.

    Attributes
    ----------
    symbols : dict
        A map from special numerical values to symbols for use in `__str__`

    Methods
    -------
    contains(x)
        Determine whether the argument is contained within the domain (True)
        or not (False). Used for input validation.
    get_numerical_endpoints()
        Gets the numerical values of the domain endpoints, which may have been
        defined symbolically.
    __str__()
        Returns a text representation of the domain (e.g. ``[0, b)``).
        Used for generating documentation.

    """
    symbols = {np.inf: r"\infty", -np.inf: r"-\infty", np.pi: r"\pi", -np.pi: r"-\pi"}

    @abstractmethod
    def contains(self, x):
        raise NotImplementedError()

    @abstractmethod
    def draw(self, n):
        raise NotImplementedError()

    @abstractmethod
    def get_numerical_endpoints(self, x):
        raise NotImplementedError()

    @abstractmethod
    def __str__(self):
        raise NotImplementedError()


class _SimpleDomain(_Domain):
    r""" Representation of a simply-connected domain defined by two endpoints.

    Each endpoint may be a finite scalar, positive or negative infinity, or
    be given by a single parameter. The domain may include the endpoints or
    not.

    This class still does not provide an implementation of the __str__ method,
    so it is meant for subclassing (e.g. a subclass for domains on the real
    line).

    Attributes
    ----------
    symbols : dict
        Inherited. A map from special values to symbols for use in `__str__`.
    endpoints : 2-tuple of float(s) and/or str(s)
        A tuple with two values. Each may be either a float (the numerical
        value of the endpoints of the domain) or a string (the name of the
        parameters that will define the endpoint).
    inclusive : 2-tuple of bools
        A tuple with two boolean values; each indicates whether the
        corresponding endpoint is included within the domain or not.

    Methods
    -------
    define_parameters(*parameters)
        Records any parameters used to define the endpoints of the domain
    get_numerical_endpoints(parameter_values)
        Gets the numerical values of the domain endpoints, which may have been
        defined symbolically.
    contains(item, parameter_values)
        Determines whether the argument is contained within the domain

    """
    def __init__(self, endpoints=(-inf, inf), inclusive=(False, False)):
        a, b = endpoints
        self.endpoints = np.asarray(a)[()], np.asarray(b)[()]
        self.inclusive = inclusive
        # self.all_inclusive = (endpoints == (-inf, inf)
        #                       and inclusive == (True, True))

    def define_parameters(self, *parameters):
        r""" Records any parameters used to define the endpoints of the domain.

        Adds the keyword name of each parameter and its text representation
        to the  `symbols` attribute as key:value pairs.
        For instance, a parameter may be passed into to a distribution's
        initializer using the keyword `log_a`, and the corresponding
        string representation may be '\log(a)'. To form the text
        representation of the domain for use in documentation, the
        _Domain object needs to map from the keyword name used in the code
        to the string representation.

        Returns None, but updates the `symbols` attribute.

        Parameters
        ----------
        *parameters : _Parameter objects
            Parameters that may define the endpoints of the domain.

        """
        new_symbols = {param.name: param.symbol for param in parameters}
        self.symbols.update(new_symbols)

    def get_numerical_endpoints(self, parameter_values):
        r""" Get the numerical values of the domain endpoints.

        Domain endpoints may be defined symbolically. This returns numerical
        values of the endpoints given numerical values for any variables.

        Parameters
        ----------
        parameter_values : dict
            A dictionary that maps between string variable names and numerical
            values of parameters, which may define the endpoints.

        Returns
        -------
        a, b : ndarray
            Numerical values of the endpoints

        """
        # TODO: ensure outputs are floats
        a, b = self.endpoints
        # If `a` (`b`) is a string - the name of the parameter that defines
        # the endpoint of the domain - then corresponding numerical values
        # will be found in the `parameter_values` dictionary. Otherwise, it is
        # itself the array of numerical values of the endpoint.
        try:
            a = np.asarray(parameter_values.get(a, a))
            b = np.asarray(parameter_values.get(b, b))
        except TypeError as e:
            message = ("The endpoints of the distribution are defined by "
                       "parameters, but their values were not provided. When "
                       f"using a private method of {self.__class__}, pass "
                       "all required distribution parameters as keyword "
                       "arguments.")
            raise TypeError(message) from e

        return a, b

    def contains(self, item, parameter_values=None):
        r"""Determine whether the argument is contained within the domain.

        Parameters
        ----------
        item : ndarray
            The argument
        parameter_values : dict
            A dictionary that maps between string variable names and numerical
            values of parameters, which may define the endpoints.

        Returns
        -------
        out : bool
            True if `item` is within the domain; False otherwise.

        """
        parameter_values = parameter_values or {}
        # if self.all_inclusive:
        #     # Returning a 0d value here makes things much faster.
        #     # I'm not sure if it's safe, though. If it causes a bug someday,
        #     # I guess it wasn't.
        #     # Even if there is no bug because of the shape, it is incorrect for
        #     # `contains` to return True when there are invalid (e.g. NaN)
        #     # parameters.
        #     return np.asarray(True)

        a, b = self.get_numerical_endpoints(parameter_values)
        left_inclusive, right_inclusive = self.inclusive

        in_left = item >= a if left_inclusive else item > a
        in_right = item <= b if right_inclusive else item < b
        return in_left & in_right


class _RealDomain(_SimpleDomain):
    r""" Represents a simply-connected subset of the real line; i.e., an interval

    Completes the implementation of the `_SimpleDomain` class for simple
    domains on the real line.

    Methods
    -------
    define_parameters(*parameters)
        (Inherited) Records any parameters used to define the endpoints of the
        domain.
    get_numerical_endpoints(parameter_values)
        (Inherited) Gets the numerical values of the domain endpoints, which
        may have been defined symbolically.
    contains(item, parameter_values)
        (Inherited) Determines whether the argument is contained within the
        domain
    __str__()
        Returns a string representation of the domain, e.g. "[a, b)".
    draw(size, rng, proportions, parameter_values)
        Draws random values based on the domain. Proportions of values within
        the domain, on the endpoints of the domain, outside the domain,
        and having value NaN are specified by `proportions`.

    """

    def __str__(self):
        a, b = self.endpoints
        left_inclusive, right_inclusive = self.inclusive

        left = "[" if left_inclusive else "("
        a = self.symbols.get(a, f"{a}")
        right = "]" if right_inclusive else ")"
        b = self.symbols.get(b, f"{b}")

        return f"{left}{a}, {b}{right}"

    def draw(self, n, type_, min, max, squeezed_base_shape, rng=None):
        r""" Draw random values from the domain.

        Parameters
        ----------
        n : int
            The number of values to be drawn from the domain.
        type_ : str
            A string indicating whether the values are

            - strictly within the domain ('in'),
            - at one of the two endpoints ('on'),
            - strictly outside the domain ('out'), or
            - NaN ('nan').
        min, max : ndarray
            The endpoints of the domain.
        squeezed_based_shape : tuple of ints
            See _RealParameter.draw.
        rng : np.Generator
            The Generator used for drawing random values.

        """
        rng = np.random.default_rng(rng)

        # get copies of min and max with no nans so that uniform doesn't fail
        min_nn, max_nn = min.copy(), max.copy()
        i = np.isnan(min_nn) | np.isnan(max_nn)
        min_nn[i] = 0
        max_nn[i] = 1

        shape = (n,) + squeezed_base_shape

        if type_ == 'in':
            z = rng.uniform(min_nn, max_nn, size=shape)

        elif type_ == 'on':
            z_on_shape = shape
            z = np.ones(z_on_shape)
            i = rng.random(size=n) < 0.5
            z[i] = min
            z[~i] = max

        elif type_ == 'out':
            # make this work for infinite bounds
            z = min_nn - rng.uniform(size=shape)
            zr = max_nn + rng.uniform(size=shape)
            i = rng.random(size=n) < 0.5
            z[i] = zr[i]

        elif type_ == 'nan':
            z = np.full(shape, np.nan)

        return z


class _IntegerDomain(_SimpleDomain):
    r""" Representation of a domain of consecutive integers.

    Completes the implementation of the `_SimpleDomain` class for domains
    composed of consecutive integer values.

    To be completed when needed.
    """
    def __init__(self):
        raise NotImplementedError


#### Parameters

class _Parameter(ABC):
    r""" Representation of a distribution parameter or variable.

    A `_Parameter` object is responsible for storing information about a
    parameter or variable, providing input validation/standardization of
    values passed for that parameter, providing a text/mathematical
    representation of the parameter for the documentation (`__str__`), and
    drawing random values of itself for testing and benchmarking. It does
    not provide a complete implementation of this functionality and is meant
    for subclassing.

    Attributes
    ----------
    name : str
        The keyword used to pass numerical values of the parameter into the
        initializer of the distribution
    symbol : str
        The text representation of the variable in the documentation. May
        include LaTeX.
    domain : _Domain
        The domain of the parameter for which the distribution is valid.
    typical : 2-tuple of floats or strings (consider making a _Domain)
        Defines the endpoints of a typical range of values of the parameter.
        Used for sampling.

    Methods
    -------
    __str__():
        Returns a string description of the variable for use in documentation,
        including the keyword used to represent it in code, the symbol used to
        represent it mathemtatically, and a description of the valid domain.
    draw(size, *, rng, domain, proportions)
        Draws random values of the parameter. Proportions of values within
        the valid domain, on the endpoints of the domain, outside the domain,
        and having value NaN are specified by `proportions`.
    validate(x):
        Validates and standardizes the argument for use as numerical values
        of the parameter.

   """
    def __init__(self, name, *, domain, symbol=None, typical=None):
        self.name = name
        self.symbol = symbol or name
        self.domain = domain
        if typical is not None and not isinstance(typical, _Domain):
            typical = _RealDomain(typical)
        self.typical = typical or domain

    def __str__(self):
        r""" String representation of the parameter for use in documentation."""
        return f"`{self.name}` for :math:`{self.symbol} \\in {str(self.domain)}`"

    def draw(self, size=None, *, rng=None, region='domain', proportions=None,
             parameter_values=None):
        r""" Draw random values of the parameter for use in testing.

        Parameters
        ----------
        size : tuple of ints
            The shape of the array of valid values to be drawn.
        rng : np.Generator
            The Generator used for drawing random values.
        region : str
            The region of the `_Parameter` from which to draw. Default is
            "domain" (the *full* domain); alternative is "typical". An
            enhancement would give a way to interpolate between the two.
        proportions : tuple of numbers
            A tuple of four non-negative numbers that indicate the expected
            relative proportion of elements that:

            - are strictly within the domain,
            - are at one of the two endpoints,
            - are strictly outside the domain, and
            - are NaN,

            respectively. Default is (1, 0, 0, 0). The number of elements in
            each category is drawn from the multinomial distribution with
            `np.prod(size)` as the number of trials and `proportions` as the
            event probabilities. The values in `proportions` are automatically
            normalized to sum to 1.
        parameter_values : dict
            Map between the names of parameters (that define the endpoints of
            `typical`) and numerical values (arrays).

        """
        parameter_values = parameter_values or {}
        domain = self.domain
        proportions = (1, 0, 0, 0) if proportions is None else proportions

        pvals = proportions / np.sum(proportions)

        a, b = domain.get_numerical_endpoints(parameter_values)
        a, b = np.broadcast_arrays(a, b)

        base_shape = a.shape
        extended_shape = np.broadcast_shapes(size, base_shape)
        n_extended = np.prod(extended_shape)
        n_base = np.prod(base_shape)
        n = int(n_extended / n_base) if n_extended else 0

        rng = np.random.default_rng(rng)
        n_in, n_on, n_out, n_nan = rng.multinomial(n, pvals)

        # `min` and `max` can have singleton dimensions that correspond with
        # non-singleton dimensions in `size`. We need to be careful to avoid
        # shuffling results (e.g. a value that was generated for the domain
        # [min[i], max[i]] ends up at index j). To avoid this:
        # - Squeeze the singleton dimensions out of `min`/`max`. Squeezing is
        #   often not the right thing to do, but here is equivalent to moving
        #   all the dimensions that are singleton in `min`/`max` (which may be
        #   non-singleton in the result) to the left. This is what we want.
        # - Now all the non-singleton dimensions of the result are on the left.
        #   Ravel them to a single dimension of length `n`, which is now along
        #   the 0th axis.
        # - Reshape the 0th axis back to the required dimensions, and move
        #   these axes back to their original places.
        base_shape_padded = ((1,)*(len(extended_shape) - len(base_shape))
                             + base_shape)
        base_singletons = np.where(np.asarray(base_shape_padded)==1)[0]
        new_base_singletons = tuple(range(len(base_singletons)))
        # Base singleton dimensions are going to get expanded to these lengths
        shape_expansion = np.asarray(extended_shape)[base_singletons]

        # assert(np.prod(shape_expansion) == n)  # check understanding
        # min = np.reshape(min, base_shape_padded)
        # max = np.reshape(max, base_shape_padded)
        # min = np.moveaxis(min, base_singletons, new_base_singletons)
        # max = np.moveaxis(max, base_singletons, new_base_singletons)
        # squeezed_base_shape = max.shape[len(base_singletons):]
        # assert np.all(min.reshape(squeezed_base_shape) == min.squeeze())
        # assert np.all(max.reshape(squeezed_base_shape) == max.squeeze())

        # min = np.maximum(a, _fiinfo(a).min/10) if np.any(np.isinf(a)) else a
        # max = np.minimum(b, _fiinfo(b).max/10) if np.any(np.isinf(b)) else b
        min = np.asarray(a.squeeze())
        max = np.asarray(b.squeeze())
        squeezed_base_shape = max.shape

        if region == 'typical':
            typical = self.typical
            a, b = typical.get_numerical_endpoints(parameter_values)
            a, b = np.broadcast_arrays(a, b)
            min_here = np.asarray(a.squeeze())
            max_here = np.asarray(b.squeeze())
            z_in = typical.draw(n_in, 'in', min_here, max_here, squeezed_base_shape,
                                rng=rng)
        else:
            z_in = domain.draw(n_in, 'in', min, max, squeezed_base_shape, rng=rng)
        z_on = domain.draw(n_on, 'on', min, max, squeezed_base_shape, rng=rng)
        z_out = domain.draw(n_out, 'out', min, max, squeezed_base_shape, rng=rng)
        z_nan= domain.draw(n_nan, 'nan', min, max, squeezed_base_shape, rng=rng)

        z = np.concatenate((z_in, z_on, z_out, z_nan), axis=0)
        z = rng.permuted(z, axis=0)

        z = np.reshape(z, tuple(shape_expansion) + squeezed_base_shape)
        z = np.moveaxis(z, new_base_singletons, base_singletons)
        return z

    @abstractmethod
    def validate(self, arr):
        raise NotImplementedError()


class _RealParameter(_Parameter):
    r""" Represents a real-valued parameter.

    Implements the remaining methods of _Parameter for real parameters.
    All attributes are inherited.

    """
    def validate(self, arr, parameter_values):
        r""" Input validation/standardization of numerical values of a parameter.

        Checks whether elements of the argument `arr` are reals, ensuring that
        the dtype reflects this. Also produces a logical array that indicates
        which elements meet the requirements.

        Parameters
        ----------
        arr : ndarray
            The argument array to be validated and standardized.
        parameter_values : dict
            Map of parameter names to parameter value arrays.

        Returns
        -------
        arr : ndarray
            The argument array that has been validated and standardized
            (converted to an appropriate dtype, if necessary).
        dtype : NumPy dtype
            The appropriate floating point dtype of the parameter.
        valid : boolean ndarray
            Logical array indicating which elements are valid (True) and
            which are not (False). The arrays of all distribution parameters
            will be broadcasted, and elements for which any parameter value
            does not meet the requirements will be replaced with NaN.

        """
        arr = np.asarray(arr)

        valid_dtype = None
        # minor optimization - fast track the most common types to avoid
        # overhead of np.issubdtype. Checking for `in {...}` doesn't work : /
        if arr.dtype == np.float64 or arr.dtype == np.float32:
            pass
        elif arr.dtype == np.int32 or arr.dtype == np.int64:
            arr = np.asarray(arr, dtype=np.float64)
        elif np.issubdtype(arr.dtype, np.floating):
            pass
        elif np.issubdtype(arr.dtype, np.integer):
            arr = np.asarray(arr, dtype=np.float64)
        else:
            message = f"Parameter `{self.name}` must be of real dtype."
            raise TypeError(message)

        valid = self.domain.contains(arr, parameter_values)
        valid = valid & valid_dtype if valid_dtype is not None else valid

        return arr[()], arr.dtype, valid


#### Parameterizations

class _Parameterization:
    r""" Represents a parameterization of a distribution.

    Distributions can have multiple parameterizations. A `_Parameterization`
    object is responsible for recording the parameters used by the
    parameterization, checking whether keyword arguments passed to the
    distribution match the parameterization, and performing input validation
    of the numerical values of these parameters.

    Attributes
    ----------
    parameters : dict
        String names (of keyword arguments) and the corresponding _Parameters.

    Methods
    -------
    __len__()
        Returns the number of parameters in the parameterization.
    __str__()
        Returns a string representation of the parameterization.
    copy
        Returns a copy of the parameterization. This is needed for transformed
        distributions that add parameters to the parameterization.
    matches(parameters)
        Checks whether the keyword arguments match the parameterization.
    validation(parameter_values)
        Input validation / standardization of parameterization. Validates the
        numerical values of all parameters.
    draw(sizes, rng, proportions)
        Draw random values of all parameters of the parameterization for use
        in testing.
    """
    def __init__(self, *parameters):
        self.parameters = {param.name: param for param in parameters}

    def __len__(self):
        return len(self.parameters)

    def copy(self):
        return _Parameterization(*self.parameters.values())

    def matches(self, parameters):
        r""" Checks whether the keyword arguments match the parameterization.

        Parameters
        ----------
        parameters : set
            Set of names of parameters passed into the distribution as keyword
            arguments.

        Returns
        -------
        out : bool
            True if the keyword arguments names match the names of the
            parameters of this parameterization.
        """
        return parameters == set(self.parameters.keys())

    def validation(self, parameter_values):
        r""" Input validation / standardization of parameterization.

        Parameters
        ----------
        parameter_values : dict
            The keyword arguments passed as parameter values to the
            distribution.

        Returns
        -------
        all_valid : ndarray
            Logical array indicating the elements of the broadcasted arrays
            for which all parameter values are valid.
        dtype : dtype
            The common dtype of the parameter arrays. This will determine
            the dtype of the output of distribution methods.
        """
        all_valid = True
        dtypes = set()  # avoid np.result_type if there's only one type
        for name, arr in parameter_values.items():
            parameter = self.parameters[name]
            arr, dtype, valid = parameter.validate(arr, parameter_values)
            dtypes.add(dtype)
            all_valid = all_valid & valid
            parameter_values[name] = arr
        dtype = arr.dtype if len(dtypes)==1 else np.result_type(*list(dtypes))

        return all_valid, dtype

    def __str__(self):
        r"""Returns a string representation of the parameterization."""
        messages = [str(param) for name, param in self.parameters.items()]
        return ", ".join(messages)

    def draw(self, sizes=None, rng=None, proportions=None, region='domain'):
        r"""Draw random values of all parameters for use in testing.

        Parameters
        ----------
        sizes : iterable of shape tuples
            The size of the array to be generated for each parameter in the
            parameterization. Note that the order of sizes is arbitary; the
            size of the array generated for a specific parameter is not
            controlled individually as written.
        rng : NumPy Generator
            The generator used to draw random values.
        proportions : tuple
            A tuple of four non-negative numbers that indicate the expected
            relative proportion of elements that are within the parameter's
            domain, are on the boundary of the parameter's domain, are outside
            the parameter's domain, and have value NaN. For more information,
            see the `draw` method of the _Parameter subclasses.
        domain : str
            The domain of the `_Parameter` from which to draw. Default is
            "domain" (the *full* domain); alternative is "typical".

        Returns
        -------
        parameter_values : dict (string: array)
            A dictionary of parameter name/value pairs.
        """
        # ENH: be smart about the order. The domains of some parameters
        # depend on others. If the relationshp is simple (e.g. a < b < c),
        # we can draw values in order a, b, c.
        parameter_values = {}

        if not len(sizes) or not np.iterable(sizes[0]):
            sizes = [sizes]*len(self.parameters)

        for size, param in zip(sizes, self.parameters.values()):
            parameter_values[param.name] = param.draw(
                size, rng=rng, proportions=proportions,
                parameter_values=parameter_values,
                region=region
            )

        return parameter_values


#### Decorators

def _set_invalid_nan(f):
    # Wrapper for input / output validation and standardization of distribution
    # functions that accept either the quantile or percentile as an argument:
    # logpdf, pdf
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
                 'ilogcdf': (-np.inf, 0), 'ilogccdf': (-np.inf, 0)}
    replacements = {'logpdf': (-inf, -inf), 'pdf': (0, 0),
                    '_logcdf1': (-inf, 0), '_logccdf1': (0, -inf),
                    '_cdf1': (0, 1), '_ccdf1': (1, 0)}
    replace_strict = {'pdf', 'logpdf'}
    replace_exact = {'icdf', 'iccdf', 'ilogcdf', 'ilogccdf'}
    clip = {'_cdf1', '_ccdf1'}
    clip_log = {'_logcdf1', '_logccdf1'}

    @functools.wraps(f)
    def filtered(self, x, *args, **kwargs):
        if self.validation_policy == _SKIP_ALL:
            return f(self, x, *args, **kwargs)

        method_name = f.__name__
        x = np.asarray(x)
        dtype = self._dtype
        shape = self._shape

        # Ensure that argument is at least as precise as distribution
        # parameters, which are already at least floats. This will avoid issues
        # with raising integers to negative integer powers and failure to replace
        # invalid integers with NaNs.
        if x.dtype != dtype:
            dtype = np.result_type(x.dtype, dtype)
            x = np.asarray(x, dtype=dtype)

        # Broadcasting is slow. Do it only if necessary.
        if not x.shape == shape:
            try:
                shape = np.broadcast_shapes(x.shape, shape)
                x = np.broadcast_to(x, shape)
                # Should we broadcast the distribution parameters to this shape, too?
            except ValueError as e:
                message = (
                    f"The argument provided to `{self.__class__.__name__}"
                    f".{method_name}` cannot be be broadcast to the same "
                    "shape as the distribution parameters.")
                raise ValueError(message) from e

        low, high = endpoints.get(method_name, self.support())

        # Check for arguments outside of domain. They'll be replaced with NaNs,
        # and the result will be set to the appropriate value.
        left_inc, right_inc = self._variable.domain.inclusive
        mask_low = (x < low if (method_name in replace_strict and left_inc)
                    else x <= low)
        mask_high = (x > high if (method_name in replace_strict and right_inc)
                     else x >= high)
        mask_invalid = (mask_low | mask_high)
        any_invalid = (mask_invalid if mask_invalid.shape == ()
                       else np.any(mask_invalid))

        # Check for arguments at domain endpoints, whether they
        # are part of the domain or not.
        any_endpoint = False
        if method_name in replace_exact:
            mask_low_endpoint = (x == low)
            mask_high_endpoint = (x == high)
            mask_endpoint = (mask_low_endpoint | mask_high_endpoint)
            any_endpoint = (mask_endpoint if mask_endpoint.shape == ()
                            else np.any(mask_endpoint))

        # Set out-of-domain arguments to NaN. The result will be set to the
        # appropriate value later.
        if any_invalid:
            x = np.array(x, dtype=dtype, copy=True)
            x[mask_invalid] = np.nan

        res = np.asarray(f(self, x, *args, **kwargs))

        # Ensure that the result is the correct dtype and shape,
        # copying (only once) if necessary.
        res_needs_copy = False
        if res.dtype != dtype:
            dtype = np.result_type(dtype, self._dtype)
            res_needs_copy = True

        if res.shape != shape:  # faster to check first
            res = np.broadcast_to(res, self._shape)
            res_needs_copy = res_needs_copy or any_invalid or any_endpoint

        if res_needs_copy:
            res = np.array(res, dtype=dtype, copy=True)

        #  For arguments outside the function domain, replace results
        if any_invalid:
            replace_low, replace_high = (
                replacements.get(method_name, (np.nan, np.nan)))
            res[mask_low] = replace_low
            res[mask_high] = replace_high

        # For arguments at the endpoints of the domain, replace results
        if any_endpoint:
            a, b = self.support()
            if a.shape != shape:
                a = np.array(np.broadcast_to(a, shape), copy=True)
                b = np.array(np.broadcast_to(b, shape), copy=True)

            replace_low_endpoint = (
                b[mask_low_endpoint] if method_name.endswith('ccdf')
                else a[mask_low_endpoint])
            replace_high_endpoint = (
                a[mask_high_endpoint] if method_name.endswith('ccdf')
                else b[mask_high_endpoint])

            res[mask_low_endpoint] = replace_low_endpoint
            res[mask_high_endpoint] = replace_high_endpoint

        # Clip probabilities to [0, 1]
        if method_name in clip:
            res = np.clip(res, 0., 1.)
        elif method_name in clip_log:
            res = res.real  # exp(res) > 0
            res = np.clip(res, None, 0.)  # exp(res) < 1

        return res[()]

    return filtered


def _set_invalid_nan_property(f):
    # Wrapper for input / output validation and standardization of distribution
    # functions that represent properties of the distribution itself:
    # logentropy, entropy
    # median, mode
    # moment
    # It ensures that the output is of the correct shape and dtype and that
    # there are NaNs wherever the distribution parameters were invalid.

    @functools.wraps(f)
    def filtered(self, *args, **kwargs):
        if self.validation_policy == _SKIP_ALL:
            return f(self, *args, **kwargs)

        res = f(self, *args, **kwargs)
        if res is None:
            # message could be more appropriate
            raise NotImplementedError(self._not_implemented)

        res = np.asarray(res)
        needs_copy = False
        dtype = res.dtype

        if dtype != self._dtype:  # this won't work for logmoments (complex)
            dtype = np.result_type(dtype, self._dtype)
            needs_copy = True

        if res.shape != self._shape:  # faster to check first
            res = np.broadcast_to(res, self._shape)
            needs_copy = needs_copy or self._any_invalid

        if needs_copy:
            res = res.astype(dtype=dtype, copy=True)

        if self._any_invalid:
            # may be redundant when quadrature is used, but not necessarily
            # when formulas are used.
            res[self._invalid] = np.nan

        return res[()]

    return filtered


def _cdf2_input_validation(f):
    # Wrapper that does the job of `_set_invalid_nan` when `cdf` or `logcdf`
    # is called with two quantile arguments.
    # Let's keep it simple; no special cases for speed right now.
    # The strategy is a bit different than for 1-arg `cdf` (and other methods
    # covered by `_set_invalid_nan`). For 1-arg `cdf`, elements of `x` that
    # are outside (or at the edge of) the support get replaced by `nan`,
    # and then the results get replaced by the appropriate value (0 or 1).
    # We *could* do something similar, dispatching to `_cdf1` in these
    # cases. That would be a bit more robust, but it would also be quite
    # a bit more complex, since we'd have to do different things when
    # `x` and `y` are both out of bounds, when just `x` is out of bounds,
    # when just `y` is out of bounds, and when both are out of bounds.
    # I'm not going to do that right now. Instead, simply replace values
    # outside the support by those at the edge of the support. Here, we also
    # omit some of the optimizations that make `_set_invalid_nan` faster for
    # simple arguments (e.g. float64 scalars).

    @functools.wraps(f)
    def wrapped(self, x, y, *args, **kwargs):
        func_name = f.__name__

        low, high = self.support()
        x, y, low, high = np.broadcast_arrays(x, y, low, high)
        dtype = np.result_type(x.dtype, y.dtype, self._dtype)
        # yes, copy to avoid modifying input arrays
        x, y = x.astype(dtype, copy=True), y.astype(dtype, copy=True)

        # Swap arguments to ensure that x < y, and replace
        # out-of domain arguments with domain endpoints. We'll
        # transform the result later.
        i_swap = y < x
        x[i_swap], y[i_swap] = y[i_swap], x[i_swap]
        i = x < low
        x[i] = low[i]
        i = y < low
        y[i] = low[i]
        i = x > high
        x[i] = high[i]
        i = y > high
        y[i] = high[i]

        res = f(self, x, y, *args, **kwargs)

        # Clipping probability to [0, 1]
        if func_name in {'_cdf2', '_ccdf2'}:
            res = np.clip(res, 0., 1.)
        else:
            res = np.clip(res, None, 0.)  # exp(res) < 1

        # Transform the result to account for swapped argument order
        res = np.asarray(res)
        if func_name == '_cdf2':
            res[i_swap] *= -1.
        elif func_name == '_ccdf2':
            res[i_swap] *= -1
            res[i_swap] += 2.
        elif func_name == '_logcdf2':
            res = np.asarray(res + 0j) if np.any(i_swap) else res
            res[i_swap] = res[i_swap] + np.pi*1j
        else:
            # res[i_swap] is always positive and less than 1, so it's
            # safe to ensure that the result is real
            res[i_swap] = _logexpxmexpy(np.log(2), res[i_swap]).real
        return res[()]

    return wrapped


def _dispatch(f):
    # For each public method (instance function) of a distribution (e.g. ccdf),
    # there may be several ways ("method"s) that it can be computed (e.g. a
    # formula, as the complement of the CDF, or via numerical integration).
    # Each "method" is implemented by a different private method (instance
    # function).
    # This wrapper calls the appropriate private method based on the public
    # method and any specified `method` keyword option.
    # - If `method` is specified as a string (by the user), the appropriate
    #   private method is called.
    # - If `method` is None:
    #   - The appropriate private method for the public method is looked up
    #     in a cache.
    #   - If the cache does not have an entry for the public method, the
    #     appropriate "dispatch " function is called to determine which method
    #     is most appropriate given the available private methods and
    #     settings (e.g. tolerance).

    @functools.wraps(f)
    def wrapped(self, *args, method=None, **kwargs):
        func_name = f.__name__
        method = method or self._method_cache.get(func_name, None)
        if callable(method):
            pass
        elif method is not None:
            method = 'logexp' if method == 'log/exp' else method
            method_name = func_name.replace('dispatch', method)
            method = getattr(self, method_name)
        else:
            method = f(self, *args, method=method, **kwargs)
            if self.cache_policy != _NO_CACHE:
                self._method_cache[func_name] = method

        try:
            return method(*args, **kwargs)
        except KeyError as e:
            raise NotImplementedError(self._not_implemented) from e

    return wrapped


#### Utilities

def _fiinfo(x):
    if np.issubdtype(x.dtype, np.inexact):
        return np.finfo(x.dtype)
    else:
        return np.iinfo(x)


def _kwargs2args(f, args=None, kwargs=None):
    # Wraps a function that accepts a primary argument `x`, secondary
    # arguments `args`, and secondary keyward arguments `kwargs` such that the
    # wrapper accepts only `x` and `args`. The keyword arguments are extracted
    # from `args` passed into the wrapper, and these are passed to the
    # underlying function as `kwargs`.
    # This is a temporary workaround until the scalar algorithms `_tanhsinh`,
    # `_chandrupatla`, etc., support `kwargs` or can operate with compressing
    # arguments to the callable.
    args = args or []
    kwargs = kwargs or {}
    names = list(kwargs.keys())
    n_args = len(args)

    def wrapped(x, *args):
        return f(x, *args[:n_args], **dict(zip(names, args[n_args:])))

    args = list(args) + list(kwargs.values())

    return wrapped, args


def _log1mexp(x):
    r"""Compute the log of the complement of the exponential.

    This function is equivalent to::

        log1mexp(x) = np.log(1-np.exp(x))

    but avoids loss of precision when ``np.exp(x)`` is nearly 0 or 1.

    Parameters
    ----------
    x : array_like
        Input array.

    Returns
    -------
    y : ndarray
        An array of the same shape as `x`.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import log1m
    >>> x = 1e-300  # log of a number very close to 1
    >>> _log1mexp(x)  # log of the complement of a number very close to 1
    -690.7755278982137
    >>> # p.log(1 - np.exp(x))  # -inf; emits warning

    """
    def f1(x):
        # good for exp(x) close to 0
        return np.log1p(-np.exp(x))

    def f2(x):
        # good for exp(x) close to 1
        return np.real(np.log(-special.expm1(x + 0j)))

    return _lazywhere(x < -1, (x,), f=f1, f2=f2)[()]


def _logexpxmexpy(x, y):
    """ Compute the log of the difference of the exponentials of two arguments.

    Avoids over/underflow, but does not prevent loss of precision otherwise.
    """
    # TODO: properly avoid NaN when y is negative infinity
    # TODO: silence warning with taking log of complex nan
    # TODO: deal with x == y better
    i = np.isneginf(np.real(y))
    if np.any(i):
        y = np.asarray(y.copy())
        y[i] = np.finfo(y.dtype).min
    x, y = np.broadcast_arrays(x, y)
    res = np.asarray(special.logsumexp([x, y+np.pi*1j], axis=0))
    i = (x == y)
    res[i] = -np.inf
    return res


def _log_real_standardize(x):
    """Standardizes the (complex) logarithm of a real number.

    The logarithm of a real number may be represented by a complex number with
    imaginary part that is a multiple of pi*1j. Even multiples correspond with
    a positive real and odd multiples correspond with a negative real.

    Given a logarithm of a real number `x`, this function returns an equivalent
    representation in a standard form: the log of a positive real has imaginary
    part `0` and the log of a negative real has imaginary part `pi`.

    """
    shape = x.shape
    x = np.atleast_1d(x)
    real = np.real(x).astype(x.dtype)
    complex = np.imag(x)
    y = real
    negative = np.exp(complex*1j) < 0.5
    y[negative] = y[negative] + np.pi * 1j
    return y.reshape(shape)[()]


#### Documentation

def _combine_docs(dist_family):
    fields = set(NumpyDocString.sections)
    fields.remove('index')

    doc = ClassDoc(dist_family)
    superdoc = ClassDoc(_ProbabilityDistribution)
    for field in fields:
        if field in {"Methods", "Attributes"}:
            doc[field] = superdoc[field]
        elif field in {"Summary"}:
            pass
        elif field == "Extended Summary":
            doc[field].append(_generate_domain_support(dist_family))
        elif field == 'Examples':
            doc[field] = [dist_family._generate_example()]
        else:
            doc[field] += superdoc[field]
    return str(doc)


def _generate_domain_support(dist_family):
    n_parameterizations = len(dist_family._parameterizations)

    domain = f"\nfor :math:`x` in {dist_family._variable.domain}.\n"

    if n_parameterizations == 0:
        support = """
        This class accepts no distribution parameters.
        """
    elif n_parameterizations == 1:
        support = f"""
        This class accepts one parameterization:
        {str(dist_family._parameterizations[0])}.
        """
    else:
        number = {2: 'two', 3: 'three', 4: 'four', 5: 'five'}[
            n_parameterizations]
        parameterizations = [f"- {str(p)}" for p in
                             dist_family._parameterizations]
        parameterizations = "\n".join(parameterizations)
        support = f"""
        This class accepts {number} parameterizations:

        {parameterizations}
        """
    support = "\n".join([line.lstrip() for line in support.split("\n")][1:])
    return domain + support


#### Univariate Distribution Base Class

class _ProbabilityDistribution(ABC):
    _parameterizations = []  # type: ignore[var-annotated]

    ### Initialization

    def __init__(self, *, tol=_null, validation_policy=None, cache_policy=None,
                 **parameters):
        self.tol = tol
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
        parameters = {key: val for key, val in parameters.items() if val is not None}
        self._update_parameters(**parameters)

    def _update_parameters(self, *, validation_policy=None, **params):
        r""" Update the numerical values of distribution parameters.

        Parameters
        ----------
        **params : array_like
            Desired numerical values of the distribution parameters. Any or all
            of the parameters initially used to instantiate the distribution
            may be modified. Parameters used in alternative parameterizations
            are not accepted.

        validation_policy : str
            To be documented. See Question 3 at the top.
        """

        parameters = original_parameters = self._original_parameters.copy()
        parameters.update(**params)
        parameterization = None
        self._invalid = np.asarray(False)
        self._any_invalid = False
        self._shape = tuple()
        self._ndim = 0
        self._size = 1
        self._dtype = np.float64

        if (validation_policy or self.validation_policy) == _SKIP_ALL:
            parameters = self._process_parameters(**parameters)
        elif not len(self._parameterizations):
            if parameters:
                message = (f"The `{self.__class__.__name__}` distribution "
                           "family does not accept parameters, but parameters "
                           f"`{set(parameters)}` were provided.")
                raise ValueError(message)
        else:
            # This is default behavior, which re-runs all parameter validations
            # even when only a single parameter is modified. For many
            # distributions, the domain of a parameter doesn't depend on other
            # parameters, so parameters could safely be modified without
            # re-validating all other parameters. To handle these cases more
            # efficiently, we could allow the developer  to override this
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
            parameters, shape, size, ndim = self._broadcast(parameters)
            parameters, invalid, any_invalid, dtype = (
                self._validate(parameterization, parameters))
            parameters = self._process_parameters(**parameters)

            self._invalid = invalid
            self._any_invalid = any_invalid
            self._shape = shape
            self._size = size
            self._ndim = ndim
            self._dtype = dtype

        self.reset_cache()
        self._parameters = parameters
        self._parameterization = parameterization
        self._original_parameters = original_parameters
        for name in self._parameters.keys():
            # Make parameters properties of the class; return values from the instance
            if hasattr(self.__class__, name):
                continue
            setattr(self.__class__, name, property(lambda self_, name_=name:
                                                   self_._parameters[name_].copy()[()]))

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

    def _broadcast(self, parameters):
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
        parameter_vals = [np.asarray(parameter)
                          for parameter in parameters.values()]
        parameter_shapes = set(parameter.shape for parameter in parameter_vals)
        if len(parameter_shapes) == 1:
            return (parameters, parameter_vals[0].shape,
                    parameter_vals[0].size, parameter_vals[0].ndim)

        try:
            parameter_vals = np.broadcast_arrays(*parameter_vals)
        except ValueError as e:
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
        valid, dtype = parameterization.validation(parameters)
        invalid = ~valid
        any_invalid = invalid if invalid.shape == () else np.any(invalid)
        # If necessary, make the arrays contiguous and replace invalid with NaN
        if any_invalid:
            for parameter_name in parameters:
                parameters[parameter_name] = np.copy(
                    parameters[parameter_name])
                parameters[parameter_name][invalid] = np.nan

        return parameters, invalid, any_invalid, dtype

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
        parameter_names_list = list(parameters.keys())
        parameter_names_list.sort()
        return f"{{{', '.join(parameter_names_list)}}}"

    def _copy_parameterization(self):
        self._parameterizations = self._parameterizations.copy()
        for i in range(len(self._parameterizations)):
            self._parameterizations[i] = self._parameterizations[i].copy()

    ### Documentation

    @classmethod
    def _generate_example(cls):
        n_parameters = cls._num_parameters(0)
        shapes = [()] * n_parameters
        rng = np.random.default_rng(615681484984984)
        i = 0
        dist = cls._draw(shapes, rng=rng, i_parameterization=i)

        rng = np.random.default_rng(2354873452)
        name = cls.__name__
        if n_parameters:
            parameter_names = list(dist._parameterizations[i].parameters)
            parameter_values = [round(getattr(dist, name), 2) for name in
                                parameter_names]
            name_values = [f"{name}={value}" for name, value in
                           zip(parameter_names, parameter_values)]
            instantiation = f"{name}({', '.join(name_values)})"
            attributes = ", ".join([f"X.{param}" for param in dist._parameters])
            X = cls(**dict(zip(parameter_names, parameter_values)))
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
        To evaluate the probability density function of the underlying distribution
        at argument ``x={x}``:

        >>> x = {x}
        >>> X.pdf(x)
        {X.pdf(x)}

        The cumulative distribution function, its complement, and the logarithm
        of these functions are evaluated similarly.

        >>> np.allclose(np.exp(X.logccdf(x)), 1 - X.cdf(x))
        True

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

        There are methods for computing measures of central tendency,
        dispersion, higher moments, and entropy.

        >>> X.mean(), X.median(), X.mode()
        {X.mean(), X.median(), X.mode()}
        >>> X.variance(), X.standard_deviation()
        {X.variance(), X.standard_deviation()}
        >>> X.skewness(), X.kurtosis()
        {X.skewness(), X.kurtosis()}
        >>> np.allclose(X.moment(order=6, kind='standardized'),
        ...             X.moment(order=6, kind='central') / X.variance()**3)
        True
        >>> np.allclose(np.exp(X.logentropy()), X.entropy())
        True

        Pseudo-random samples can be drawn from
        the underlying distribution using ``sample``.

        >>> X.sample(shape=(4,))
        {repr(X.sample(shape=(4,)))}  # may vary
        """
        # remove the indentation due to use of block quote within function;
        # eliminate blank first line
        example = "\n".join([line.lstrip() for line in example.split("\n")][1:])
        return example

    ### Attributes

    # `tol` attribute is just notional right now. See Question 4 above.
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
        if _isnull(tol):
            self._tol = tol
            return

        tol = np.asarray(tol)
        if (tol.shape != () or not tol > 0 or  # catches NaNs
                not np.issubdtype(tol.dtype, np.floating)):
            message = (f"Attribute `tol` of `{self.__class__.__name__}` must "
                       "be a positive float, if specified.")
            raise ValueError(message)
        self._tol = tol[()]

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

    ### Magic methods

    def __repr__(self):
        r""" Returns a string representation of the distribution.

        Includes the name of the distribution family, the names of the
        parameters, and the broadcasted shape and result dtype of the
        parameters.

        """
        class_name = self.__class__.__name__
        parameters = list(self._original_parameters.items())
        info = []
        if parameters:
            parameters.sort()
            if self._size <= 3:
                str_parameters = [f"{symbol}={value}" for symbol, value in parameters]
                str_parameters = f"{', '.join(str_parameters)}"
            else:
                str_parameters = f"{', '.join([symbol for symbol, _ in parameters])}"
            info.append(str_parameters)
        if self._shape:
            info.append(f"shape={self._shape}")
        if self._dtype != np.float64:
            info.append(f"dtype={self._dtype}")
        return f"{class_name}({', '.join(info)})"

    def __add__(self, loc):
        return ShiftedScaledDistribution(self, loc=loc)

    def __sub__(self, loc):
        return ShiftedScaledDistribution(self, loc=-loc)

    def __mul__(self, scale):
        return ShiftedScaledDistribution(self, scale=scale)

    def __truediv__(self, scale):
        return ShiftedScaledDistribution(self, scale=1/scale)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rtruediv__(self, other):
        message = "Division by a random variable is not yet implemented."
        raise NotImplementedError(message)

    def __neg__(self):
        return self * -1

    ## Input validation

    def _validate_order_kind(self, order, kind, kinds):
        # Yet another integer validating function. Unlike others in SciPy, it
        # Is quite flexible about what is allowed as an integer, and it
        # raises a distribution-specific error message to facilitate
        # identification of the source of the error.
        if self.validation_policy == _SKIP_ALL:
            return order

        order = np.asarray(order, dtype=self._dtype)[()]
        message = (f"Argument `order` of `{self.__class__.__name__}.moment` "
                   "must be a finite, positive integer.")
        try:
            order_int = round(order.item())
            # If this fails for any reason (e.g. it's an array, it's infinite)
            # it's not a valid `order`.
        except Exception as e:
            raise ValueError(message) from e

        if order_int <0 or order_int != order:
            raise ValueError(message)

        message = (f"Argument `kind` of `{self.__class__.__name__}.moment` "
                   f"must be one of {set(kinds)}.")
        if kind.lower() not in kinds:
            raise ValueError(message)

        return order

    def _preserve_type(self, x):
        x = np.asarray(x)
        if x.dtype != self._dtype:
            x = x.astype(self._dtype)
        return x[()]

    ## Algorithms

    def _solve_bounded(self, f, p, *, bounds=None, params=None):
        # Finds the argument of a function that produces the desired output.
        # Much of this should be added to _bracket_root / _chandrupatla.
        xmin, xmax = self._support(**params) if bounds is None else bounds
        params = {} if params is None else params

        p, xmin, xmax = np.broadcast_arrays(p, xmin, xmax)
        if not p.size:
            # might need to figure out result type based on p
            return np.empty(p.shape, dtype=self._dtype)

        def f2(x, p, **kwargs):
            return f(x, **kwargs) - p

        f3, args = _kwargs2args(f2, args=[p], kwargs=params)
        # If we know the median or mean, should use it

        # Any operations between 0d array and a scalar produces a scalar, so...
        shape = xmin.shape
        xmin, xmax = np.atleast_1d(xmin, xmax)

        a = -np.ones_like(xmin)
        b = np.ones_like(xmax)

        i = np.isfinite(xmin) & np.isfinite(xmax)
        a[i] = xmin[i]
        b[i] = xmax[i]

        i = np.isfinite(xmin) & ~np.isfinite(xmax)
        a[i] = xmin[i]
        b[i] = xmin[i] + 1

        i = np.isfinite(xmax) & ~np.isfinite(xmin)
        a[i] = xmax[i] - 1
        b[i] = xmax[i]

        xmin = xmin.reshape(shape)
        xmax = xmax.reshape(shape)
        a = a.reshape(shape)
        b = b.reshape(shape)

        res = _bracket_root(f3, xl0=a, xr0=b, xmin=xmin, xmax=xmax, args=args)
        # For now, we ignore the status, but I want to use the bracket width
        # as an error estimate - see question 5 at the top.
        xrtol = None if _isnull(self.tol) else self.tol
        return _chandrupatla(f3, a=res.xl, b=res.xr, args=args, xrtol=xrtol).x

    ## Testing

    @classmethod
    def _draw(cls, sizes=None, rng=None, i_parameterization=None,
              proportions=None):
        r""" Draw a specific (fully-defined) distribution from the family.

        See _Parameterization.draw for documentation details.
        """
        rng = np.random.default_rng(rng)
        if len(cls._parameterizations) == 0:
            return cls()
        if i_parameterization is None:
            n = cls._num_parameterizations()
            i_parameterization = rng.integers(0, max(0, n - 1), endpoint=True)

        parameterization = cls._parameterizations[i_parameterization]
        parameters = parameterization.draw(sizes, rng, proportions=proportions,
                                           region='typical')
        return cls(**parameters)

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

    ## Public API and helpers

    def support(self):
        r"""Support of the random variable

        The support of a random variable is set of all possible outcomes;
        i.e., the subset of the domain of argument :math:`x` for which
        the probability density function :math:`f(x)` is nonzero.

        This function returns lower and upper bounds of the support.

        Returns
        -------
        out : tuple of Array
            The lower and upper bounds of the support.

        See Also
        --------
        pdf

        References
        ----------
        .. [1] Support (mathematics), *Wikipedia*,
               https://en.wikipedia.org/wiki/Support_(mathematics)

        Notes
        -----
        Suppose a continuous probability distribution has support ``(l, r)``.
        The following table summarizes the value returned by methods
        of ``ContinuousDistribution`` for arguments outside the support.

        +----------------+---------------------+---------------------+
        | Method         | Value for ``x < l`` | Value for ``x > r`` |
        +================+=====================+=====================+
        | ``pdf(x)``     | 0                   | 0                   |
        +----------------+---------------------+---------------------+
        | ``logpdf(x)``  | -inf                | -inf                |
        +----------------+---------------------+---------------------+
        | ``cdf(x)``     | 0                   | 1                   |
        +----------------+---------------------+---------------------+
        | ``logcdf(x)``  | -inf                | 0                   |
        +----------------+---------------------+---------------------+
        | ``ccdf(x)``    | 1                   | 0                   |
        +----------------+---------------------+---------------------+
        | ``logccdf(x)`` | 0                   | -inf                |
        +----------------+---------------------+---------------------+

        For the ``cdf`` and related methods, the inequality need not be
        strict; i.e. the tabulated value is returned when the method is
        evaluated *at* the corresponding boundary.

        The following table summarizes the value returned by the inverse
        methods of ``ContinuousDistribution`` for arguments at the boundaries
        of the domain ``0`` to ``1``.

        +-------------+-----------+-----------+
        | Method      | ``x = 0`` | ``x = 1`` |
        +=============+===========+===========+
        | ``icdf(x)`` | ``l``     | ``r``     |
        +-------------+-----------+-----------+
        | ``icdf(x)`` | ``r``     | ``l``     |
        +-------------+-----------+-----------+

        For the inverse log-functions, the same values are returned for
        for ``x = log(0)`` and ``x = log(1)``. All inverse functions return
        ``nan`` when evaluated at an argument outside the domain ``0`` to ``1``.

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Retrieve the support of the distribution:

        >>> X.support()
        (-0.5, 0.5)

        For a distribution with infinite support,

        >>> X = stats.Normal()
        >>> X.support()
        (-inf, inf)

        Due to underflow, the numerical value returned by the PDF may be zero
        even for arguments within the support, even if the true value is
        nonzero. In such cases, the log-PDF may be useful.

        >>> X.pdf([-100., 100.])
        array([0., 0.])
        >>> X.logpdf([-100., 100.])
        array([-5000.91893853, -5000.91893853])

        Use cases for the log-CDF and related methods are analogous.

        """
        # If this were a `cached_property`, we couldn't update the value
        # when the distribution parameters change.
        # Caching is important, though, because calls to _support take a few
        # microseconds even when `a` and `b` are already the same shape.
        if self._support_cache is not None:
            return self._support_cache

        a, b = self._support(**self._parameters)
        if a.shape != self._shape:
            a = np.broadcast_to(a, self._shape)
        if b.shape != self._shape:
            b = np.broadcast_to(b, self._shape)

        if self._any_invalid:
            a, b = np.asarray(a).copy(), np.asarray(b).copy()
            a[self._invalid], b[self._invalid] = np.nan, np.nan
            a, b = a[()], b[()]

        support = (a, b)

        if self.cache_policy != _NO_CACHE:
            self._support_cache = support

        return support

    def _support(self, **params):
        # Computes the support given distribution parameters
        a, b = self._variable.domain.get_numerical_endpoints(params)
        if len(params):
            # the parameters should all be of the same dtype and shape at this point
            vals = list(params.values())
            shape = vals[0].shape
            a = np.broadcast_to(a, shape) if a.shape != shape else a
            b = np.broadcast_to(b, shape) if b.shape != shape else b
        return self._preserve_type(a), self._preserve_type(b)

    def sample(self, shape=(), *, method=None, rng=None):
        r"""Random sample from the distribution.

        Parameters
        ----------
        shape : tuple of ints, default: ()
            The shape of the sample to draw. If the parameters of the distribution
            underlying the random variable are arrays of shape ``param_shape``,
            the output array will be of shape ``shape + param_shape``.
        method : {None, 'formula', 'inverse_transform'}
            The strategy used to produce the sample. By default (``None``),
            the infrastructure chooses between the following options,
            listed in order of precedence.

            - ``'formula'``: an implementation specific to the distribution
            - ``'inverse_transform'``: generate a uniformly distributed sample and
              return the inverse CDF at these arguments.

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a `NotImplementedError``
            will be raised.
        rng : `numpy.random.Generator`, optional
            Pseudorandom number generator state. When `rng` is None, a new
            `numpy.random.Generator` is created using entropy from the
            operating system. Types other than `numpy.random.Generator` are
            passed to `numpy.random.default_rng` to instantiate a `Generator`.

        References
        ----------
        .. [1] Sampling (statistics), *Wikipedia*,
               https://en.wikipedia.org/wiki/Sampling_(statistics)

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=0., b=1.)

        Generate a pseudorandom sample:

        >>> x = X.sample((1000, 1))
        >>> octiles = (np.arange(8) + 1) / 8
        >>> np.count_nonzero(x <= octiles, axis=0)
        array([ 148,  263,  387,  516,  636,  751,  865, 1000])  # may vary

        >>> X = stats.Uniform(a=np.zeros((3, 1)), b=np.ones(2))
        >>> X.a.shape,
        (3, 2)
        >>> x = X.sample(shape=(5, 4))
        >>> x.shape
        (5, 4, 3, 2)

        """
        # needs output validation to ensure that developer returns correct
        # dtype and shape
        sample_shape = (shape,) if not np.iterable(shape) else tuple(shape)
        full_shape = sample_shape + self._shape
        rng = np.random.default_rng(rng)
        res = self._sample_dispatch(sample_shape, full_shape, method=method,
                                    rng=rng, **self._parameters)

        return res.astype(self._dtype, copy=False)

    @_dispatch
    def _sample_dispatch(self, sample_shape, full_shape, *, method, rng, **params):
        # make sure that tests catch if sample is 0d array
        if self._overrides('_sample_formula'):
            method = self._sample_formula
        else:
            method = self._sample_inverse_transform
        return method

    def _sample_formula(self, sample_shape, full_shape, *, rng, **params):
        raise NotImplementedError(self._not_implemented)

    @_set_invalid_nan_property
    def moment(self, order, kind='raw', *, method=None):
        r"""Raw, central, or standard moment of positive integer order.

        In terms of probability density function :math:`f(x)` and support
        :math:`\chi`, the "raw" moment (about the origin) of order :math:`n` of
        a random variable :math:`X` is:

        .. math::

            \mu'_n(X) = \int_{\chi} x^n f(x) dx

        The "central" moment is the raw moment taken about the mean,
        :math:`\mu = \mu'_1`:

        .. math::

            \mu_n(X) = \int_{\chi} (x - \mu) ^n f(x) dx

        The "standardized" moment is the central moment normalized by the
        :math:`n^\text{th}` power of the standard deviation
        :math:`\sigma = \sqrt{\mu_2}` to produce a scale invariant quantity:

        .. math::

            \tilde{\mu}_n(X) = \frac{\mu_n(X)}
                                    {\sigma^n}

        Parameters
        ----------
        order : int
            The integer order of the moment; i.e. :math:`n` in the formulae above.
        kind : {'raw', 'central', 'standardized'}
            Whether to return the raw (default), central, or standardized moment
            defined above.
        method : {None, 'formula', 'general', 'transform', 'normalize', 'quadrature', 'cache'}
            The strategy used to evaluate the moment. By default (``None``),
            the infrastructure chooses between the following options,
            listed in order of precedence.

            - ``'cache'``: use the value of the moment most recently calculated
              via another method
            - ``'formula'``: use a formula for the moment itself
            - ``'general'``: use a general result that is true for all distributions
              with finite moments; for instance, the zeroth raw moment is
              identically 1
            - ``'transform'``: transform a raw moment to a central moment or
              vice versa (see Notes)
            - ``'normalize'``: normalize a central moment to get a standardized
              or vice versa
            - ``'quadrature'``: numerically integrate according to the definition

            Not all `method` options are available for all orders, kinds, and
            distributions. If the selected `method` is not available, a
            ``NotImplementedError`` will be raised.

        Returns
        -------
        out : array
            The moment of the random variable of the specified order and kind.

        See Also
        --------
        pdf
        mean
        variance
        standard_deviation
        skewness
        kurtosis

        Notes
        -----
        Not all distributions have finite moments of all orders; moments of some
        orders may be undefined or infinite. If a formula for the moment is not
        specifically implemented for the chosen distribution, SciPy will attempt
        to compute the moment via a generic method, which may yield a finite
        result where none exists. This is not a critical bug, but an opportunity
        for an enhancement.

        The definition of a raw moment in the summary is specific to the raw moment
        about the origin. The raw moment about any point :math:`a` is:

        .. math::

            E[(X-a)^n] = \int_{\chi} (x-a)^n f(x) dx

        In this notation, a raw moment about the origin is :math:`\mu'_n = E[x^n]`,
        and a central moment is :math:`\mu_n = E[(x-\mu)^n]`, where :math:`\mu`
        is the first raw moment; i.e. the mean.

        The ``'transform'`` method takes advantage of the following relationships
        between moments taken about different points :math:`a` and :math:`b`.

        .. math::

            E[(X-b)^n] =  \sum_{i=0}^n E[(X-a)^i] {n \choose i} (a - b)^{n-i}

        For instance, to transform the raw moment to the central moment, we let
        :math:`b = \mu` and :math:`a = 0`.

        The distribution infrastructure provides flexibility for distribution
        authors to implement separate formulas for raw moments, central moments,
        and standardized moments of any order. By default, the moment of the
        desired order and kind is evaluated from the formula if such a formula
        is available; if not, the infrastructure uses any formulas that are
        available rather than resorting directly to numerical integration.
        For instance, if formulas for the first three raw moments are
        available and the third standardized moments is desired, the
        infrastructure will evaluate the raw moments and perform the transforms
        and standardization required. The decision tree is somewhat complex,
        but the strategy for obtaining a moment of a given order and kind
        (possibly as an intermediate step due to the recursive nature of the
        transform formula above) roughly follows this order of priority:

        #. Use cache (if order of same moment and kind has been calculated)
        #. Use formula (if available)
        #. Transform between raw and central moment and/or normalize to convert
           between central and standardized moments (if efficient)
        #. Use a generic result true for most distributions (if available)
        #. Use quadrature

        References
        ----------
        .. [1] Moment, *Wikipedia*,
               https://en.wikipedia.org/wiki/Moment_(mathematics)

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the first raw moment:

        >>> X.moment(order=1, kind='raw')
        1.0
        >>> X.moment(order=1, kind='raw') == X.mean() == X.mu
        True

        Evaluate the second central moment:

        >>> X.moment(order=2, kind='central')
        4.0
        >>> X.moment(order=2, kind='central') == X.variance() == X.sigma**2
        True

        Evaluate the fourth standardized moment:

        >>> X.moment(order=4, kind='standardized')
        3.0
        >>> X.moment(order=4, kind='standardized') == X.kurtosis(convention='non-excess')
        True

        """  # noqa:E501
        kinds = {'raw': self._moment_raw,
                 'central': self._moment_central,
                 'standardized': self._moment_standardized}
        order = self._validate_order_kind(order, kind, kinds)
        moment_kind = kinds[kind]
        return moment_kind(order, method=method)

    @functools.cached_property
    def _moment_methods(self):
        return {'cache', 'formula', 'transform',
                'normalize', 'general', 'quadrature'}

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

        if moment is None and 'transform' in methods and order > 1:
            moment = self._moment_raw_transform(order, **params)

        if moment is None and 'general' in methods:
            moment = self._moment_raw_general(order, **params)

        if moment is None and 'quadrature' in methods:
            moment = self._moment_integrate_pdf(order, center=self._zero, **params)

        if moment is None and 'quadrature_icdf' in methods:
            moment = self._moment_integrate_icdf(order, center=self._zero, **params)

        if moment is not None and self.cache_policy != _NO_CACHE:
            self._moment_raw_cache[order] = moment

        return moment

    def _moment_raw_formula(self, order, **params):
        return None

    def _moment_raw_transform(self, order, **params):
        central_moments = []
        for i in range(int(order) + 1):
            methods = {'cache', 'formula', 'normalize', 'general'}
            moment_i = self._moment_central_dispatch(order=i, methods=methods, **params)
            if moment_i is None:
                return None
            central_moments.append(moment_i)

        # Doesn't make sense to get the mean by "transform", since that's
        # how we got here. Questionable whether 'quadrature' should be here.
        mean_methods = {'cache', 'formula', 'quadrature'}
        mean = self._moment_raw_dispatch(self._one, methods=mean_methods, **params)
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
            mean = self._moment_raw_dispatch(self._one, **params,
                                             methods=self._moment_methods)
            moment = self._moment_integrate_pdf(order, center=mean, **params)

        if moment is None and 'quadrature_icdf' in methods:
            mean = self._moment_raw_dispatch(self._one, **params,
                                             methods=self._moment_methods)
            moment = self._moment_integrate_icdf(order, center=mean, **params)

        if moment is not None and self.cache_policy != _NO_CACHE:
            self._moment_central_cache[order] = moment

        return moment

    def _moment_central_formula(self, order, **params):
        return None

    def _moment_central_transform(self, order, **params):
        raw_moments = []
        for i in range(int(order) + 1):
            methods = {'cache', 'formula', 'general'}
            moment_i = self._moment_raw_dispatch(order=i, methods=methods, **params)
            if moment_i is None:
                return None
            raw_moments.append(moment_i)

        mean_methods = self._moment_methods
        mean = self._moment_raw_dispatch(self._one, methods=mean_methods, **params)

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

    def _moment_transform_center(self, order, moment_as, a, b):
        a, b, *moment_as = np.broadcast_arrays(a, b, *moment_as)
        n = order
        i = np.arange(n+1).reshape([-1]+[1]*a.ndim)  # orthogonal to other axes
        i = self._preserve_type(i)
        n_choose_i = special.binom(n, i)
        moment_b = np.sum(n_choose_i*moment_as*(a-b)**(n-i), axis=0)
        return moment_b

    def mean(self, *, method=None):
        r"""Mean (raw first moment about the origin)

        Parameters
        ----------
        method : {None, 'formula', 'transform', 'quadrature', 'cache'}
            Method used to calculate the raw first moment. Not
            all methods are available for all distributions. See
            `moment` for details.

        See Also
        --------
        moment
        median
        mode

        References
        ----------
        .. [1] Mean, *Wikipedia*,
               https://en.wikipedia.org/wiki/Mean#Mean_of_a_probability_distribution

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the variance:

        >>> X.mean()
        1.0
        >>> X.mean() == X.moment(order=1, kind='raw') == X.mu
        True

        """
        return self.moment(1, kind='raw', method=method)

    @_set_invalid_nan_property
    def median(self, *, method=None):
        r"""Median (50th percentil)

        If a continuous random variable :math:`X` has probability :math:`0.5` of
        taking on a value less than :math:`m`, then :math:`m` is the median.
        That is, the median is the value :math:`m` for which:

        .. math::

            P(X ≤ m) = 0.5 = P(X ≥ m)

        Parameters
        ----------
        method : {None, 'formula', 'icdf'}
            The strategy used to evaluate the median.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the median
            - ``'icdf'``: evaluate the inverse CDF of 0.5

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The median

        See Also
        --------
        mean
        mode
        icdf

        References
        ----------
        .. [1] Median, *Wikipedia*,
               https://en.wikipedia.org/wiki/Median#Probability_distributions

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Uniform(a=0., b=10.)

        Compute the median:

        >>> X.median()
        5
        >>> X.median() == X.icdf(0.5) == X.iccdf(0.5)
        True

        """
        return self._median_dispatch(method=method, **self._parameters)

    @_dispatch
    def _median_dispatch(self, method=None, **params):
        if self._overrides('_median_formula'):
            method = self._median_formula
        else:
            method = self._median_icdf
        return method

    def _median_formula(self, **params):
        raise NotImplementedError(self._not_implemented)

    def _median_icdf(self, **params):
        return self._icdf_dispatch(0.5, **params)

    @_set_invalid_nan_property
    def mode(self, *, method=None):
        r"""Mode (most likely value)

        Informally, the mode is a value that a random variable has the highest
        probability (density) of assuming. That is, the mode is the element of
        the support :math:`\chi` that maximizes the probability density
        function :math:`f(x)`:

        .. math::

            \text{mode} = \arg\max_{x \in \chi} f(x)

        Parameters
        ----------
        method : {None, 'formula', 'optimization'}
            The strategy used to evaluate the mode.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the median
            - ``'optimization'``: numerically maximize the PDF

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The mode

        See Also
        --------
        mean
        median
        pdf

        Notes
        -----
        For some distributions

        #. the mode is not unique (e.g. the uniform distribution);
        #. the PDF has one or more singularities, and it is debateable whether
           a singularity is considered to be in the domain and called the mode
           (e.g. the gamma distribution with shape parameter less than 1); and/or
        #. the probability density function may have one or more local maxima
           that are not a global maximum (e.g. mixture distributions).

        In such cases, `mode` will

        #. return a single value,
        #. consider the mode to occur at a singularity, and/or
        #. return a local maximum which may or may not be a global maximum.

        If a formula for the mode is not specifically implemented for the
        chosen distribution, SciPy will attempt to compute the mode
        numerically, which may not meet the user's preferred definition of a
        mode. In such cases, the user is encouraged to subclass the
        distribution and override ``mode``.

        References
        ----------
        .. [1] Mode (statistics), *Wikipedia*,
               https://en.wikipedia.org/wiki/Mode_(statistics)

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the mode:

        >>> X.mode()
        1.0

        If the mode is not uniquely defined, ``mode`` nonetheless returns a
        single value.

        >>> X = stats.Uniform(a=0., b=1.)
        >>> X.mode()
        0.5

        If this choice does not satisfy your requirements, subclass the
        distribution and override ``mode``:

        >>> class BetterUniform(stats.Uniform):
        ...     def mode(self):
        ...         return self.b
        >>> X = BetterUniform(a=0., b=1.)
        >>> X.mode()
        1.0

        """
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

    def variance(self, *, method=None):
        r"""Variance (central second moment)

        Parameters
        ----------
        method : {None, 'formula', 'transform', 'normalize', 'quadrature', 'cache'}
            Method used to calculate the central second moment. Not
            all methods are available for all distributions. See
            `moment` for details.

        See Also
        --------
        moment
        standard_deviation
        mean

        References
        ----------
        .. [1] Variance, *Wikipedia*,
               https://en.wikipedia.org/wiki/Variance#Absolutely_continuous_random_variable

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the variance:

        >>> X.variance()
        4.0
        >>> X.variance() == X.moment(order=2, kind='central') == X.sigma**2
        True

        """
        return self.moment(2, kind='central', method=method)

    def standard_deviation(self, *, method=None):
        r"""Standard deviation (square root of the second central moment)

        Parameters
        ----------
        method : {None, 'formula', 'transform', 'normalize', 'quadrature', 'cache'}
            Method used to calculate the central second moment. Not
            all methods are available for all distributions. See
            `moment` for details.

        See Also
        --------
        variance
        mean
        moment

        References
        ----------
        .. [1] Standard deviation, *Wikipedia*,
               https://en.wikipedia.org/wiki/Standard_deviation#Definition_of_population_values

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the standard deviation:

        >>> X.standard_deviation()
        2.0
        >>> X.standard_deviation() == X.moment(order=2, kind='central')**0.5 == X.sigma
        True

        """
        return np.sqrt(self.variance(method=method))

    def skewness(self, *, method=None):
        r"""Skewness (standardized third moment)

        Parameters
        ----------
        method : {None, 'formula', 'general', 'transform', 'normalize', 'cache'}
            Method used to calculate the standardized third moment. Not
            all methods are available for all distributions. See
            `moment` for details.

        See Also
        --------
        moment
        mean
        variance

        References
        ----------
        .. [1] Skewness, *Wikipedia*,
               https://en.wikipedia.org/wiki/Skewness

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the skewness:

        >>> X.skewness()
        0.0
        >>> X.skewness() == X.moment(order=3, kind='standardized')
        True

        """
        return self.moment(3, kind='standardized', method=method)

    def kurtosis(self, *, method=None, convention='non-excess'):
        r"""Kurtosis (standardized fourth moment)

        By default, this is the standardized fourth moment, also known as the
        "non-excess" or "Pearson" kurtosis (e.g. the kurtosis of the normal
        distribution is 3). The "excess" or "Fisher" kurtosis (the standardized
        fourth moment minus 3) is available via the `convention` parameter.

        Parameters
        ----------
        method : {None, 'formula', 'general', 'transform', 'normalize', 'cache'}
            Method used to calculate the standardized fourth moment. Not
            all methods are available for all distributions. See
            `moment` for details.
        convention : {'non-excess', 'excess'}
            Two distinct conventions are available:

            - ``'non-excess'``: the standardized fourth moment (Pearson's kurtosis)
            - ``'excess'``: the standardized fourth moment minus 3 (Fisher's kurtosis)

            The default is ``'non-excess'``.

        See Also
        --------
        moment
        mean
        variance

        References
        ----------
        .. [1] Kurtosis, *Wikipedia*,
               https://en.wikipedia.org/wiki/Kurtosis

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Normal(mu=1., sigma=2.)

        Evaluate the kurtosis:

        >>> X.kurtosis()
        3.0
        >>> (X.kurtosis()
        ...  == X.kurtosis(convention='excess') + 3.
        ...  == X.moment(order=4, kind='standardized'))
        True

        """
        conventions = {'non-excess', 'excess'}
        message = (f'Parameter `convention` of `{self.__class__.__name__}.kurtosis` '
                   f"must be one of {conventions}.")
        convention = convention.lower()
        if convention not in conventions:
            raise ValueError(message)
        k = self.moment(4, kind='standardized', method=method)
        return k - 3 if convention == 'excess' else k

    @_set_invalid_nan
    def pdf(self, x, /, *, method=None):
        r"""Probability density function

        The probability density function ("PDF"), denoted :math:`f(x)`, is the
        probability *per unit length* that the random variable will assume the
        value :math:`x`. Mathematically, it can be defined as the derivative
        of the cumulative distribution function :math:`F(x)`:

        .. math::

            f(x) = \frac{d}{dx} F(x)

        `pdf` accepts `x` for :math:`x`.

        Parameters
        ----------
        x : array_like
            The argument of the PDF.
        method : {None, 'formula', 'logexp'}
            The strategy used to evaluate the PDF. By default (``None``), the
            infrastructure chooses between the following options, listed in
            order of precedence.

            - ``'formula'``: use a formula for the PDF itself
            - ``'logexp'``: evaluate the log-PDF and exponentiate

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The PDF evaluated at the argument `x`.

        See Also
        --------
        cdf
        logpdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        By definition of the support, the PDF evaluates to its minimum value
        of :math:`0` outside the support; i.e. for :math:`x < l` or
        :math:`x > r`. The maximum of the PDF may be less than or greater than
        :math:`1`; since the valus is a probability *density*, only its integral
        over the support must equal :math:`1`.

        References
        ----------
        .. [1] Probability density function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Probability_density_function

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Uniform(a=-1., b=1.)

        Evaluate the PDF at the desired argument:

        >>> X.pdf(0.25)
        0.5

        """
        return self._pdf_dispatch(x, method=method, **self._parameters)

    @_dispatch
    def _pdf_dispatch(self, x, *, method=None, **params):
        if self._overrides('_pdf_formula'):
            method = self._pdf_formula
        else:
            method = self._pdf_logexp
        return method

    def _pdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _pdf_logexp(self, x, **params):
        return np.exp(self._logpdf_dispatch(x, **params))

    @_set_invalid_nan
    def logpdf(self, x, /, *, method=None):
        r"""Log of the probability density function

        The probability density function ("PDF"), denoted :math:`f(x)`, is the
        probability *per unit length* that the random variable will assume the
        value :math:`x`. Mathematically, it can be defined as the derivative
        of the cumulative distribution function :math:`F(x)`:

        .. math::

            f(x) = \frac{d}{dx} F(x)

        `logpdf` computes the logarithm of the probability density function
        ("log-PDF"), :math:`\log(f(x))`, but it may be numerically favorable
        compared to the naive implementation (computing :math:`f(x)` and
        taking the logarithm).

        `logpdf` accepts `x` for :math:`x`.

        Parameters
        ----------
        x : array_like
            The argument of the log-PDF.
        method : {None, 'formula', 'logexp'}
            The strategy used to evaluate the log-PDF. By default (``None``), the
            infrastructure chooses between the following options, listed in order
            of precedence.

            - ``'formula'``: use a formula for the log-PDF itself
            - ``'logexp'``: evaluate the PDF and takes its logarithm

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The log-PDF evaluated at the argument `x`.

        See Also
        --------
        pdf
        logcdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        By definition of the support, the log-PDF evaluates to its minimum value
        of :math:`-\infty` (i.e. :math:`\log(0)`) outside the support; i.e. for
        :math:`x < l` or :math:`x > r`. The maximum of the log-PDF may be less
        than or greater than :math:`\log(1) = 0` because the maximum of the PDF
        can be any positive real.

        For distributions with infinite support, it is common for `pdf` to return
        a value of ``0`` when the argument is theoretically within the support;
        this can occur because the true value of the PDF is too small to be
        represented by the chosen dtype. The log-PDF, however, will often be finite
        (not ``-inf``) over a much larger domain. Consequently, it may be preferred
        to work with the logarithms of probabilities and probability densities to
        avoid underflow.

        References
        ----------
        .. [1] Probability density function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Probability_density_function

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-1.0, b=1.0)

        Evaluate the log-PDF at the desired argument:

        >>> X.logpdf(0.5)
        -0.6931471805599453
        >>> np.allclose(X.logpdf(0.5), np.log(X.pdf(0.5)))
        True

        """
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
        return np.log(self._pdf_dispatch(x, **params))

    def cdf(self, x, y=None, /, *, method=None):
        r"""Cumulative distribution function

        The cumulative distribution function ("CDF"), denoted :math:`F(x)`, is
        the probability the random variable :math:`X` will assume a value
        less than or equal to :math:`x`:

        .. math::

            F(x) = P(X ≤ x)

        A two-argument variant of this function is also defined as the
        probability the random variable :math:`X` will assume a value between
        :math:`x` and :math:`y`.

        .. math::

            F(x, y) = P(x ≤ X ≤ y)

        `cdf` accepts `x` for :math:`x` and `y` for :math:`y`.

        Parameters
        ----------
        x, y : array_like
            The arguments of the CDF. `x` is required; `y` is optional.
        method : {None, 'formula', 'logexp', 'complement', 'quadrature', 'subtraction'}
            The strategy used to evaluate the CDF.
            By default (``None``), the one-argument form of the function
            chooses between the following options, listed in order of precedence.

            - ``'formula'``: use a formula for the CDF itself
            - ``'logexp'``: evaluate the log-CDF and exponentiate
            - ``'complement'``: evaluate the CCDF and take the complement
            - ``'quadrature'``: numerically integrate the PDF

            In place of ``'complement'``, the two-argument form accepts:

            - ``'subtraction'``: compute the CDF at each argument and take
              the difference.

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The CDF evaluated at the provided argument(s).

        See Also
        --------
        logcdf
        ccdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        The CDF :math:`F(x)` is related to the probability density function
        :math:`f(x)` by:

        .. math::

            F(x) = \int_l^x f(u) du

        The two argument version is:

        .. math::

            F(x, y) = \int_x^y f(u) du = F(y) - F(x)

        The CDF evaluates to its minimum value of :math:`0` for :math:`x ≤ l`
        and its maximum value of :math:`1` for :math:`x ≥ r`.

        The CDF is also known simply as the "distribution function".

        References
        ----------
        .. [1] Cumulative distribution function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Cumulative_distribution_function

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the CDF at the desired argument:

        >>> X.cdf(0.25)
        0.75

        Evaluate the cumulative probability between two arguments:

        >>> X.cdf(-0.25, 0.25) == X.cdf(0.25) - X.cdf(-0.25)
        True

        """  # noqa: E501
        if y is None:
            return self._cdf1(x, method=method)
        else:
            return self._cdf2(x, y, method=method)

    @_cdf2_input_validation
    def _cdf2(self, x, y, *, method):
        return self._cdf2_dispatch(x, y, method=method, **self._parameters)

    @_dispatch
    def _cdf2_dispatch(self, x, y, *, method=None, **params):
        # Should revisit this logic.
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
        return np.real(np.exp(self._logcdf2_dispatch(x, y, **params)))

    def _cdf2_subtraction(self, x, y, **params):
        # Improvements:
        # Lazy evaluation of cdf/ccdf only where needed
        # Stack x and y to reduce function calls?
        cdf_x = self._cdf_dispatch(x, **params)
        cdf_y = self._cdf_dispatch(y, **params)
        ccdf_x = self._ccdf_dispatch(x, **params)
        ccdf_y = self._ccdf_dispatch(y, **params)
        i = (ccdf_x < 0.5) & (ccdf_y < 0.5)
        return np.where(i, ccdf_x-ccdf_y, cdf_y-cdf_x)

    def _cdf2_subtraction_safe(self, x, y, **params):
        cdf_x = self._cdf_dispatch(x, **params)
        cdf_y = self._cdf_dispatch(y, **params)
        ccdf_x = self._ccdf_dispatch(x, **params)
        ccdf_y = self._ccdf_dispatch(y, **params)
        i = (ccdf_x < 0.5) & (ccdf_y < 0.5)
        out = np.where(i, ccdf_x-ccdf_y, cdf_y-cdf_x)

        eps = np.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else np.sqrt(eps)

        cdf_max = np.maximum(cdf_x, cdf_y)
        ccdf_max = np.maximum(ccdf_x, ccdf_y)
        spacing = np.spacing(np.where(i, ccdf_max, cdf_max))
        mask = abs(tol * out) < spacing

        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._cdf2_quadrature(x[mask], y[mask], *params_mask)
        return out[()]

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
        else:
            method = self._cdf_quadrature
        return method

    def _cdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _cdf_logexp(self, x, **params):
        return np.exp(self._logcdf_dispatch(x, **params))

    def _cdf_complement(self, x, **params):
        return 1 - self._ccdf_dispatch(x, **params)

    def _cdf_complement_safe(self, x, **params):
        ccdf = self._ccdf_dispatch(x, **params)
        out = 1 - ccdf
        eps = np.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else np.sqrt(eps)
        mask = tol * out < np.spacing(ccdf)
        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._cdf_quadrature(x[mask], *params_mask)
        return out[()]

    @_set_invalid_nan
    def icdf(self, p, /, *, method=None):
        r"""Inverse of the cumulative distribution function.

        The inverse of the cumulative distribution function ("inverse CDF"),
        denoted :math:`F^{-1}(p)`, is the argument :math:`x` for which the
        cumulative distribution function :math:`F(x)` evaluates to :math:`p`.

        .. math::

            F^{-1}(p) = x \quad \text{s.t.} \quad F(x) = p

        `icdf` accepts `p` for :math:`p \in [0, 1]`.

        Parameters
        ----------
        p : array_like
            The argument of the inverse CDF.
        method : {None, 'formula', 'complement', 'inversion'}
            The strategy used to evaluate the inverse CDF.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the inverse CDF itself
            - ``'complement'``: evaluate the inverse CCDF at the
              complement of `p`
            - ``'inversion'``: solve numerically for the argument at which the
              CDF is equal to `p`

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The inverse CDF evaluated at the provided argument.

        See Also
        --------
        cdf
        ilogcdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`. The
        inverse CDF returns its minimum value of :math:`l` at :math:`p = 0`
        and its maximum value of :math:`r` at :math:`p = 1`. Because the CDF
        has range :math:`[0, 1]`, the inverse CDF is only defined on the
        domain :math:`[0, 1]`; for :math:`p < 0` and :math:`p > 1`, `icdf`
        returns ``nan``.

        The inverse CDF is also known as the quantile function, percentile function,
        and percent-point function.

        References
        ----------
        .. [1] Quantile function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Quantile_function

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the inverse CDF at the desired argument:

        >>> X.icdf(0.25)
        -0.25
        >>> np.allclose(X.cdf(X.icdf(0.25)), 0.25)
        True

        This function returns NaN when the argument is outside the domain.

        >>> X.icdf([-0.1, 0, 1, 1.1])
        array([ nan, -0.5,  0.5,  nan])

        """
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
        out = self._icdf_complement(x, **params)
        eps = np.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else np.sqrt(eps)
        mask = tol * x < np.spacing(1 - x)
        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._icdf_inversion(x[mask], *params_mask)
        return out[()]

    def ccdf(self, x, y=None, /, *, method=None):
        r"""Complementary cumulative distribution function

        The complementary cumulative distribution function ("CCDF"), denoted
        :math:`G(x)`, is the complement of the cumulative distribution function
        :math:`F(x)`; i.e., probability the random variable :math:`X` will
        assume a value greater than :math:`x`:

        .. math::

            G(x) = 1 - F(x) = P(X > x)

        A two-argument variant of this function is:

        .. math::

            G(x, y) = 1 - F(x, y) = P(X < x \text{ or } X > y)

        `ccdf` accepts `x` for :math:`x` and `y` for :math:`y`.

        Parameters
        ----------
        x, y : array_like
            The arguments of the CCDF. `x` is required; `y` is optional.
        method : {None, 'formula', 'logexp', 'complement', 'quadrature', 'addition'}
            The strategy used to evaluate the CCDF.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the CCDF itself
            - ``'logexp'``: evaluate the log-CCDF and exponentiate
            - ``'complement'``: evaluate the CDF and take the complement
            - ``'quadrature'``: numerically integrate the PDF

            The two-argument form chooses between:

            - ``'formula'``: use a formula for the CCDF itself
            - ``'addition'``: compute the CDF at `x` and the CCDF at `y`, then add

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The CCDF evaluated at the provided argument(s).

        See Also
        --------
        cdf
        logccdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        The CCDF :math:`G(x)` is related to the probability density function
        :math:`f(x)` by:

        .. math::

            G(x) = \int_x^r f(u) du

        The two argument version is:

        .. math::

            G(x, y) = \int_l^x f(u) du + \int_y^r f(u) du

        The CCDF returns its minimum value of :math:`0` for :math:`x ≥ r`
        and its maximum value of :math:`1` for :math:`x ≤ l`.

        The CCDF is also known as the "survival function".

        References
        ----------
        .. [1] Cumulative distribution function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Cumulative_distribution_function#Derived_functions

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the CCDF at the desired argument:

        >>> X.ccdf(0.25)
        0.25
        >>> np.allclose(X.ccdf(0.25), 1-X.cdf(0.25))
        True

        Evaluate the complement of the cumulative probability between two arguments:

        >>> X.ccdf(-0.25, 0.25) == X.cdf(-0.25) + X.ccdf(0.25)
        True

        """  # noqa: E501
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
        else:
            method = self._ccdf_quadrature
        return method

    def _ccdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _ccdf_logexp(self, x, **params):
        return np.exp(self._logccdf_dispatch(x, **params))

    def _ccdf_complement(self, x, **params):
        return 1 - self._cdf_dispatch(x, **params)

    def _ccdf_complement_safe(self, x, **params):
        cdf = self._cdf_dispatch(x, **params)
        out = 1 - cdf
        eps = np.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else np.sqrt(eps)
        mask = tol * out < np.spacing(cdf)
        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._ccdf_quadrature(x[mask], *params_mask)
        return out[()]

    @_set_invalid_nan
    def iccdf(self, p, /, *, method=None):
        r"""Inverse complementary cumulative distribution function.

        The inverse complementary cumulative distribution function ("inverse CCDF"),
        denoted :math:`G^{-1}(p)`, is the argument :math:`x` for which the
        complementary cumulative distribution function :math:`G(x)` evaluates to
        :math:`p`.

        .. math::

            G^{-1}(p) = x \quad \text{s.t.} \quad G(x) = p

        `iccdf` accepts `p` for :math:`p \in [0, 1]`.

        Parameters
        ----------
        p : array_like
            The argument of the inverse CCDF.
        method : {None, 'formula', 'complement', 'inversion'}
            The strategy used to evaluate the inverse CCDF.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the inverse CCDF itself
            - ``'complement'``: evaluate the inverse CDF at the
              complement of `p`
            - ``'inversion'``: solve numerically for the argument at which the
              CCDF is equal to `p`

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The inverse CCDF evaluated at the provided argument.

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`. The
        inverse CCDF returns its minimum value of :math:`l` at :math:`p = 1`
        and its maximum value of :math:`r` at :math:`p = 0`. Because the CCDF
        has range :math:`[0, 1]`, the inverse CCDF is only defined on the
        domain :math:`[0, 1]`; for :math:`p < 0` and :math:`p > 1`, ``iccdf``
        returns ``nan``.

        See Also
        --------
        icdf
        ilogccdf

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the inverse CCDF at the desired argument:

        >>> X.iccdf(0.25)
        0.25
        >>> np.allclose(X.iccdf(0.25), X.icdf(1-0.25))
        True

        This function returns NaN when the argument is outside the domain.

        >>> X.iccdf([-0.1, 0, 1, 1.1])
        array([ nan,  0.5, -0.5,  nan])

        """
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
        out = self._iccdf_complement(x, **params)
        eps = np.finfo(self._dtype).eps
        tol = self.tol if not _isnull(self.tol) else np.sqrt(eps)
        mask = tol * x < np.spacing(1 - x)
        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._iccdf_inversion(x[mask], *params_mask)
        return out[()]

    def logcdf(self, x, y=None, /, *, method=None):
        r"""Log of the cumulative distribution function

        The cumulative distribution function ("CDF"), denoted :math:`F(x)`, is
        the probability the random variable :math:`X` will assume a value
        less than or equal to :math:`x`:

        .. math::

            F(x) = P(X ≤ x)

        A two-argument variant of this function is also defined as the
        probability the random variable :math:`X` will assume a value between
        :math:`x` and :math:`y`.

        .. math::

            F(x, y) = P(x ≤ X ≤ y)

        `logcdf` computes the logarithm of the cumulative distribution function
        ("log-CDF"), :math:`\log(F(x))`/:math:`\log(F(x, y))`, but it may be
        numerically favorable compared to the naive implementation (computing
        the CDF and taking the logarithm).

        `logcdf` accepts `x` for :math:`x` and `y` for :math:`y`.

        Parameters
        ----------
        x, y : array_like
            The arguments of the log-CDF. `x` is required; `y` is optional.
        method : {None, 'formula', 'logexp', 'complement', 'quadrature', 'subtraction'}
            The strategy used to evaluate the log-CDF.
            By default (``None``), the one-argument form of the function
            chooses between the following options, listed in order of precedence.

            - ``'formula'``: use a formula for the log-CDF itself
            - ``'logexp'``: evaluate the CDF and take the logarithm
            - ``'complement'``: evaluate the log-CCDF and take the
              logarithmic complement (see Notes)
            - ``'quadrature'``: numerically log-integrate the log-PDF

            In place of ``'complement'``, the two-argument form accepts:

            - ``'subtraction'``: compute the log-CDF at each argument and take
              the logarithmic difference (see Notes)

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The log-CDF evaluated at the provided argument(s).

        See Also
        --------
        cdf
        logccdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        The log-CDF evaluates to its minimum value of :math:`\log(0) = -\infty`
        for :math:`x ≤ l` and its maximum value of :math:`\log(1) = 0` for
        :math:`x ≥ r`.

        For distributions with infinite support, it is common for
        `cdf` to return a value of ``0`` when the argument
        is theoretically within the support; this can occur because the true value
        of the CDF is too small to be represented by the chosen dtype. `logcdf`,
        however, will often return a finite (not ``-inf``) result over a much larger
        domain. Similarly, `logcdf` may provided a strictly negative result with
        arguments for which `cdf` would return ``1.0``. Consequently, it may be
        preferred to work with the logarithms of probabilities to avoid underflow
        and related limitations of floating point numbers.

        The "logarithmic complement" of a number :math:`z` is mathematically
        equivalent to :math:`\log(1-\exp(z))`, but it is computed to avoid loss
        of precision when :math:`\exp(z)` is nearly :math:`0` or :math:`1`.
        Similarly, the term "logarithmic difference" of :math:`w` and :math:`z`
        is used here to mean :math:`\log(\exp(w)-\exp(z))`.

        If ``y < x``, the CDF is negative, and therefore the log-CCDF
        is complex with imaginary part :math:`\pi`. For
        consistency, the result of this function always has complex dtype
        when `y` is provided, regardless of the value of the imaginary part.

        References
        ----------
        .. [1] Cumulative distribution function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Cumulative_distribution_function

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the log-CDF at the desired argument:

        >>> X.logcdf(0.25)
        -0.287682072451781
        >>> np.allclose(X.logcdf(0.), np.log(X.cdf(0.)))
        True

        """  # noqa: E501
        if y is None:
            return self._logcdf1(x, method=method)
        else:
            return self._logcdf2(x, y, method=method)

    @_dispatch
    def _logcdf2_dispatch(self, x, y, *, method=None, **params):
        # dtype is complex if any x > y
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
        flip_sign = x > y
        x, y = np.minimum(x, y), np.maximum(x, y)
        logcdf_x = self._logcdf_dispatch(x, **params)
        logcdf_y = self._logcdf_dispatch(y, **params)
        logccdf_x = self._logccdf_dispatch(x, **params)
        logccdf_y = self._logccdf_dispatch(y, **params)
        case_left = (logcdf_x < -1) & (logcdf_y < -1)
        case_right = (logccdf_x < -1) & (logccdf_y < -1)
        case_central = ~(case_left | case_right)
        log_mass = _logexpxmexpy(logcdf_y, logcdf_x)
        log_mass[case_right] = _logexpxmexpy(logccdf_x, logccdf_y)[case_right]
        log_tail = np.logaddexp(logcdf_x, logccdf_y)[case_central]
        log_mass[case_central] = _log1mexp(log_tail)
        log_mass[flip_sign] += np.pi * 1j
        return log_mass[()]

    def _logcdf2_logexp(self, x, y, **params):
        expres = self._cdf2_dispatch(x, y, **params)
        expres = expres + 0j if np.any(x > y) else expres
        return np.log(expres)

    def _logcdf2_logexp_safe(self, x, y, **params):
        out = self._logcdf2_logexp(x, y, **params)
        mask = np.isinf(out.real)
        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._logcdf2_quadrature(x[mask], y[mask], **params_mask)
        return out[()]

    @_cdf2_input_validation
    def _logcdf2(self, x, y, *, method):
        out = self._logcdf2_dispatch(x, y, method=method, **self._parameters)
        return (out + 0j) if not np.issubdtype(out.dtype, np.complexfloating) else out

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
        else:
            method = self._logcdf_quadrature
        return method

    def _logcdf_formula(self, x, **params):
        raise NotImplementedError(self._not_implemented)

    def _logcdf_complement(self, x, **params):
        return _log1mexp(self._logccdf_dispatch(x, **params))

    def _logcdf_logexp(self, x, **params):
        return np.log(self._cdf_dispatch(x, **params))

    def _logcdf_logexp_safe(self, x, **params):
        out = self._logcdf_logexp(x, **params)
        mask = np.isinf(out)
        if np.any(mask):
            params_mask = {key:np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._logcdf_quadrature(x[mask], **params_mask)
        return out[()]

    @_set_invalid_nan
    def ilogcdf(self, logp, /, *, method=None):
        r"""Inverse of the logarithm of the cumulative distribution function.

        The inverse of the logarithm of the cumulative distribution function
        ("inverse log-CDF") is the argument :math:`x` for which the logarithm
        of the cumulative distribution function :math:`\log(F(x))` evaluates
        to :math:`\log(p)`.

        Mathematically, it is equivalent to :math:`F^{-1}(\exp(y))`, where
        :math:`y = \log(p)`, but it may be numerically favorable compared to
        the naive implementation (computing :math:`p = \exp(y)`, then
        :math:`F^{-1}(p)`).

        `ilogcdf` accepts `logp` for :math:`\log(p) ≤ 0`.

        Parameters
        ----------
        logp : array_like
            The argument of the inverse log-CDF.
        method : {None, 'formula', 'complement', 'inversion'}
            The strategy used to evaluate the inverse log-CDF.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the inverse log-CDF itself
            - ``'complement'``: evaluate the inverse log-CCDF at the
              logarithmic complement of `logp` (see Notes)
            - ``'inversion'``: solve numerically for the argument at which the
              log-CDF is equal to `logp`

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The inverse log-CDF evaluated at the provided argument.

        See Also
        --------
        icdf
        logcdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        The inverse log-CDF returns its minimum value of :math:`l` at
        :math:`\log(p) = \log(0) = -\infty` and its maximum value of :math:`r` at
        :math:`\log(p) = \log(1) = 0`. Because the log-CDF has range
        :math:`[-\infty, 0]`, the inverse log-CDF is only defined on the
        negative reals; for :math:`\log(p) > 0`, `ilogcdf` returns ``nan``.

        Occasionally, it is needed to find the argument of the CDF for which
        the resulting probability is very close to ``0`` or ``1`` - too close to
        represent accurately with floating point arithmetic. In many cases,
        however, the *logarithm* of this resulting probability may be
        represented in floating point arithmetic, in which case this function
        may be used to find the argument of the CDF for which the *logarithm*
        of the resulting probability is :math:`y = \log(p)`.

        The "logarithmic complement" of a number :math:`z` is mathematically
        equivalent to :math:`\log(1-\exp(z))`, but it is computed to avoid loss
        of precision when :math:`\exp(z)` is nearly :math:`0` or :math:`1`.

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the inverse log-CDF at the desired argument:

        >>> X.ilogcdf(-0.25)
        0.2788007830714034
        >>> np.allclose(X.ilogcdf(-0.25), X.icdf(np.exp(-0.25)))
        True

        """
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

    def logccdf(self, x, y=None, /, *, method=None):
        r"""Log of the complementary cumulative distribution function

        The complementary cumulative distribution function ("CCDF"), denoted
        :math:`G(x)` is the complement of the cumulative distribution function
        :math:`F(x)`; i.e., probability the random variable :math:`X` will
        assume a value greater than :math:`x`:

        .. math::

            G(x) = 1 - F(x) = P(X > x)

        A two-argument variant of this function is:

        .. math::

            G(x, y) = 1 - F(x, y) = P(X < x \quad \text{or} \quad X > y)

        `logccdf` computes the logarithm of the complementary cumulative
        distribution function ("log-CCDF"), :math:`\log(G(x))`/:math:`\log(G(x, y))`,
        but it may be numerically favorable compared to the naive implementation
        (computing the CDF and taking the logarithm).

        `logccdf` accepts `x` for :math:`x` and `y` for :math:`y`.

        Parameters
        ----------
        x, y : array_like
            The arguments of the log-CCDF. `x` is required; `y` is optional.
        method : {None, 'formula', 'logexp', 'complement', 'quadrature', 'addition'}
            The strategy used to evaluate the log-CCDF.
            By default (``None``), the one-argument form of the function
            chooses between the following options, listed in order of precedence.

            - ``'formula'``: use a formula for the log CCDF itself
            - ``'logexp'``: evaluate the CCDF and take the logarithm
            - ``'complement'``: evaluate the log-CDF and take the
              logarithmic complement (see Notes)
            - ``'quadrature'``: numerically log-integrate the log-PDF

            The two-argument form chooses between:

            - ``'formula'``: use a formula for the log CCDF itself
            - ``'addition'``: compute the log-CDF at `x` and the log-CCDF at `y`,
              then take the logarithmic sum (see Notes)

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The log-CCDF evaluated at the provided argument(s).

        See Also
        --------
        ccdf
        logcdf

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`.
        The log-CCDF returns its minimum value of :math:`\log(0)=-\infty` for
        :math:`x ≥ r` and its maximum value of :math:`\log(1) = 0` for
        :math:`x ≤ l`.

        For distributions with infinite support, it is common for
        `ccdf` to return a value of ``0`` when the argument
        is theoretically within the support; this can occur because the true value
        of the CCDF is too small to be represented by the chosen dtype. The log
        of the CCDF, however, will often be finite (not ``-inf``) over a much larger
        domain. Similarly, `logccdf` may provided a strictly negative result with
        arguments for which `ccdf` would return ``1.0``. Consequently, it may be
        preferred to work with the logarithms of probabilities to avoid underflow
        and related limitations of floating point numbers.

        The "logarithmic complement" of a number :math:`z` is mathematically
        equivalent to :math:`\log(1-\exp(z))`, but it is computed to avoid loss
        of precision when :math:`\exp(z)` is nearly :math:`0` or :math:`1`.
        Similarly, the term "logarithmic sum" of :math:`w` and :math:`z`
        is used here to mean the :math:`\log(\exp(w)+\exp(z))`, AKA
        :math:`\text{LogSumExp}(w, z)`.

        References
        ----------
        .. [1] Cumulative distribution function, *Wikipedia*,
               https://en.wikipedia.org/wiki/Cumulative_distribution_function#Derived_functions

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the log-CCDF at the desired argument:

        >>> X.logccdf(0.25)
        -1.3862943611198906
        >>> np.allclose(X.logccdf(0.), np.log(X.ccdf(0.)))
        True

        """  # noqa: E501
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
        else:
            method = self._logccdf_quadrature
        return method

    def _logccdf_formula(self):
        raise NotImplementedError(self._not_implemented)

    def _logccdf_complement(self, x, **params):
        return _log1mexp(self._logcdf_dispatch(x, **params))

    def _logccdf_logexp(self, x, **params):
        return np.log(self._ccdf_dispatch(x, **params))

    def _logccdf_logexp_safe(self, x, **params):
        out = self._logccdf_logexp(x, **params)
        mask = np.isinf(out)
        if np.any(mask):
            params_mask = {key: np.broadcast_to(val, mask.shape)[mask]
                           for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._logccdf_quadrature(x[mask], **params_mask)
        return out[()]

    @_set_invalid_nan
    def ilogccdf(self, logp, /, *, method=None):
        r"""Inverse of the log of the complementary cumulative distribution function.

        The inverse of the logarithm of the complementary cumulative distribution
        function ("inverse log-CCDF") is the argument :math:`x` for which the logarithm
        of the complementary cumulative distribution function :math:`\log(G(x))`
        evaluates to :math:`\log(p)`.

        Mathematically, it is equivalent to :math:`G^{-1}(\exp(y))`, where
        :math:`y = \log(p)`, but it may be numerically favorable compared to the naive
        implementation (computing :math:`p = \exp(y)`, then :math:`G^{-1}(p)`).

        `ilogccdf` accepts `logp` for :math:`\log(p) ≤ 0`.

        Parameters
        ----------
        x : array_like
            The argument of the inverse log-CCDF.
        method : {None, 'formula', 'complement', 'inversion'}
            The strategy used to evaluate the inverse log-CCDF.
            By default (``None``), the infrastructure chooses between the
            following options, listed in order of precedence.

            - ``'formula'``: use a formula for the inverse log-CCDF itself
            - ``'complement'``: evaluate the inverse log-CDF at the
              logarithmic complement of `x` (see Notes)
            - ``'inversion'``: solve numerically for the argument at which the
              log-CCDF is equal to `x`

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The inverse log-CCDF evaluated at the provided argument.

        Notes
        -----
        Suppose a continuous probability distribution has support :math:`[l, r]`. The
        inverse log-CCDF returns its minimum value of :math:`l` at
        :math:`\log(p) = \log(1) = 0` and its maximum value of :math:`r` at
        :math:`\log(p) = \log(0) = -\infty`. Because the log-CCDF has range
        :math:`[-\infty, 0]`, the inverse log-CDF is only defined on the
        negative reals; for :math:`\log(p) > 0`, `ilogccdf` returns ``nan``.

        Occasionally, it is needed to find the argument of the CCDF for which
        the resulting probability is very close to ``0`` or ``1`` - too close to
        represent accurately with floating point arithmetic. In many cases,
        however, the *logarithm* of this resulting probability may be
        represented in floating point arithmetic, in which case this function
        may be used to find the argument of the CCDF for which the *logarithm*
        of the resulting probability is `y = \log(p)`.

        The "logarithmic complement" of a number :math:`z` is mathematically
        equivalent to :math:`\log(1-\exp(z))`, but it is computed to avoid loss
        of precision when :math:`\exp(z)` is nearly :math:`0` or :math:`1`.

        See Also
        --------
        iccdf
        ilogccdf

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-0.5, b=0.5)

        Evaluate the inverse log-CCDF at the desired argument:

        >>> X.ilogccdf(-0.25)
        -0.2788007830714034
        >>> np.allclose(X.ilogccdf(-0.25), X.iccdf(np.exp(-0.25)))
        True

        """
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

    @_set_invalid_nan_property
    def logentropy(self, *, method=None):
        r"""Logarithm of the differential entropy

        In terms of probability density function :math:`f(x)` and support
        :math:`\chi`, the differential entropy (or simply "entropy") of a random
        variable :math:`X` is:

        .. math::

            h(X) = - \int_{\chi} f(x) \log f(x) dx

        `logentropy` computes the logarithm of the differential entropy
        ("log-entropy"), :math:`log(h(X))`, but it may be numerically favorable
        compared to the naive implementation (computing :math:`h(X)` then
        taking the logarithm).

        Parameters
        ----------
        method : {None, 'formula', 'logexp', 'quadrature}
            The strategy used to evaluate the log-entropy. By default
            (``None``), the infrastructure chooses between the following options,
            listed in order of precedence.

            - ``'formula'``: use a formula for the log-entropy itself
            - ``'logexp'``: evaluate the entropy and take the logarithm
            - ``'quadrature'``: numerically log-integrate the logarithm of the
              entropy integrand

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The log-entropy.

        See Also
        --------
        entropy
        logpdf

        Notes
        -----
        If the entropy of a distribution is negative, then the log-entropy
        is complex with imaginary part :math:`\pi`. For
        consistency, the result of this function always has complex dtype,
        regardless of the value of the imaginary part.

        References
        ----------
        .. [1] Differential entropy, *Wikipedia*,
               https://en.wikipedia.org/wiki/Differential_entropy

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> import numpy as np
        >>> from scipy import stats
        >>> X = stats.Uniform(a=-1., b=1.)

        Evaluate the log-entropy:

        >>> X.logentropy()
        (-0.3665129205816642+0j)
        >>> np.allclose(np.exp(X.logentropy()), X.entropy())
        True

        For a random variable with negative entropy, the log-entropy has an
        imaginary part equal to `np.pi`.

        >>> X = stats.Uniform(a=-.1, b=.1)
        >>> X.entropy(), X.logentropy()
        (-1.6094379124341007, (0.4758849953271105+3.141592653589793j))

        """
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
        res = np.log(self._entropy_dispatch(**params) + 0j)
        return _log_real_standardize(res)

    def _logentropy_logexp_safe(self, **params):
        out = self._logentropy_logexp(**params)
        mask = np.isinf(out.real)
        if np.any(mask):
            params_mask = {key:val[mask] for key, val in params.items()}
            out = np.asarray(out)
            out[mask] = self._logentropy_quadrature(**params_mask)
        return out[()]

    @_set_invalid_nan_property
    def entropy(self, *, method=None):
        r"""Differential entropy

        In terms of probability density function :math:`f(x)` and support
        :math:`\chi`, the differential entropy (or simply "entropy") of a
        continuous random variable :math:`X` is:

        .. math::

            h(X) = - \int_{\chi} f(x) \log f(x) dx

        Parameters
        ----------
        method : {None, 'formula', 'logexp', 'quadrature'}
            The strategy used to evaluate the entropy. By default (``None``),
            the infrastructure chooses between the following options, listed
            in order of precedence.

            - ``'formula'``: use a formula for the entropy itself
            - ``'logexp'``: evaluate the log-entropy and exponentiate
            - ``'quadrature'``: use numerical integration

            Not all `method` options are available for all distributions.
            If the selected `method` is not available, a ``NotImplementedError``
            will be raised.

        Returns
        -------
        out : array
            The entropy of the random variable.

        See Also
        --------
        logentropy
        pdf

        Notes
        -----
        This function calculates the entropy using the natural logarithm; i.e.
        the logarithm with base :math:`e`. Consequently, the value is expressed
        in (dimensionless) "units" of nats. To convert the entropy to different
        units (i.e. corresponding with a different base), divide the result by
        the natural logarithm of the desired base.

        References
        ----------
        .. [1] Differential entropy, *Wikipedia*,
               https://en.wikipedia.org/wiki/Differential_entropy

        Examples
        --------
        Instantiate a distribution with the desired parameters:

        >>> from scipy import stats
        >>> X = stats.Uniform(a=-1., b=1.)

        Evaluate the entropy:

        >>> X.entropy()
        0.6931471805599454

        """
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
        return np.real(np.exp(self._logentropy_dispatch(**params)))

#### Commentary

# TODO:
#  Test sample dtypes
#  Add dtype kwarg (especially for distributions with no parameters)
#  When drawing endpoint/out-of-bounds values of a parameter, draw them from
#   the endpoints/out-of-bounds region of the full `domain`, not `typical`.
#  Distributions without shape parameters probably need to accept a `dtype` parameter;
#    right now they default to float64. If we have them default to float16, they will
#    need to determine result_type when input is not float16 (overhead).
#  Test _solve_bounded bracket logic, and decide what to do about warnings
#  Get test coverage to 100%
#  Raise when distribution method returns wrong shape/dtype?
#  Consider ensuring everything is at least 1D for calculations? Would avoid needing
#    to sprinkle `np.asarray` throughout due to indescriminate conversion of 0D arrays
#    to scalars
#  Break up `test_basic`: test each method separately
#  Fix `sample` for QMCEngine (implementation does not match documentation)
#  When a parameter is invalid, set only the offending parameter to NaN (if possible)?
#  `_tanhsinh` special case when there are no abscissae between the limits
#    example: cdf of uniform betweeen 1.0 and np.nextafter(1.0, np.inf)
#  check behavior of moment methods when moments are undefined/infinite -
#    basically OK but needs tests
#  investigate use of median
#  implement symmetric distribution
#  implement composite distribution
#  implement wrapped distribution
#  implement folded distribution
#  implement double distribution
#  profile/optimize
#  general cleanup (choose keyword-only parameters)
#  compare old/new distribution timing
#  make video
#  PR
#  add array API support
#  why does dist.ilogcdf(-100) not converge to bound? Check solver response to inf
#  _chandrupatla_minimize should not report xm = fm = NaN when it fails
#  integrate `logmoment` into `moment`? (Not hard, but enough time and code
#   complexity to wait for reviewer feedback before adding.)
#  Eliminate bracket_root error "`min <= a < b <= max` must be True"
#  Test repr?
#  use `median` information to improve integration? In some cases this will
#   speed things up. If it's not needed, it may be about twice as slow. I think
#   it should depend on the accuracy setting.
#  in tests, check reference value against that produced using np.vectorize?
#  add `axis` to `ks_1samp`
#  User tips for faster execution:
#  - pass NumPy arrays
#  - pass inputs of floating point type (not integers)
#  - prefer NumPy scalars or 0d arrays over other size 1 arrays
#  - pass no invalid parameters and disable invalid parameter checks with iv_profile
#  - provide a Generator if you're going to do sampling
#  add options for drawing parameters: log-spacing
#  accuracy benchmark suite
#  Should caches be attributes so we can more easily ensure that they are not
#   modified when caching is turned off?
#  Make ShiftedScaledDistribution more efficient - only process underlying
#   distribution parameters as necessary.
#  Reconsider `all_inclusive`
#  Should process_parameters update kwargs rather than returning? Should we
#   update parameters rather than setting to what process_parameters returns?

# Questions:
# 1.  I override `__getattr__` so that distribution parameters can be read as
#     attributes. We don't want uses to try to change them.
#     - To prevent replacements (dist.a = b), I could override `__setattr__`.
#     - To prevent in-place modifications, `__getattr__` could return a copy,
#       or it could set the WRITEABLE flag of the array to false.
#     Which should I do?
# 2.  `cache_policy` is supported in several methods where I imagine it being
#     useful, but it needs to be tested. Before doing that:
#     - What should the default value be?
#     - What should the other values be?
#     Or should we just eliminate this policy?
# 3.  `validation_policy` is supported in a few places, but it should be checked for
#     consistency. I have the same questions as for `cache_policy`.
# 4.  `tol` is currently notional. I think there needs to be way to set
#     separate `atol` and `rtol`. Some ways I imagine it being used:
#     - Values can be passed to iterative functions (quadrature, root-finder).
#     - To control which "method" of a distribution function is used. For
#       example, if `atol` is set to `1e-12`, it may be acceptable to compute
#       the complementary CDF as 1 - CDF even when CDF is nearly 1; otherwise,
#       a (potentially more time-consuming) method would need to be used.
#     I'm looking for unified suggestions for the interface, not ad hoc ideas
#     for using tolerances. Suppose the user wants to have more control over
#     the tolerances used for each method - how do they specify it? It would
#     probably be easiest for the user if they could pass tolerances into each
#     method, but it's easiest for us if they can only set it as a property of
#     the class. Perhaps a dictionary of tolerance settings?
# 5.  I also envision that accuracy estimates should be reported to the user
#     somehow. I think my preference would be to return a subclass of an array
#     with an `error` attribute - yes, really. But this is unlikely to be
#     popular, so what are other ideas? Again, we need a unified vision here,
#     not just pointing out difficulties (not all errors are known or easy
#     to estimate, what to do when errors could compound, etc.).
# 6.  The term "method" is used to refer to public instance functions,
#     private instance functions, the "method" string argument, and the means
#     of calculating the desired quantity (represented by the string argument).
#     For the sake of disambiguation, shall I rename the "method" string to
#     "strategy" and refer to the means of calculating the quantity as the
#     "strategy"?

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
# object because iterative algorithms like `_tanhsinh` and `_chandrupatla`
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
# rng - the NumPy Generator object to used for drawing random numbers.
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
# quantities should make us question that choice. It can still accomodate
# these methods reasonably efficiently.


#### Transformed Distributions and Similar

class TransformedDistribution(_ProbabilityDistribution):
    # TODO: This may need some sort of default `_parameterizations` with a
    #       single `_Parameterization` that has no parameters. The reason is
    #       that `dist`'s parameters need to get added to it. If they're not
    #       added, then those parameter kwargs are not recognized in
    #       `_update_parameters`.
    def __init__(self, dist, *args, **kwargs):
        self._copy_parameterization()
        self._variable = dist._variable
        self._dist = dist
        if dist._parameterization:
            # Add standard distribution parameters to our parameterization
            dist_parameters = dist._parameterization.parameters
            set_params = set(dist_parameters)
            for parameterization in self._parameterizations:
                if set_params.intersection(parameterization.parameters):
                    message = (f"One or more of the parameters of {dist} has "
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
        self._dist.reset_cache()
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
        s = super().__repr__()
        return s.replace("Distribution",
                         self._dist.__class__.__name__)


# Shift/scale distributions. The purpose of
# making it a separate class is for
# a) simplicity of the ContinuousDistribution class and
# b) avoiding the requirement that every distribution accept loc/scale.
# The simplicity of ContinuousDistribution is important, because there are
# several other distribution transformations to be supported; e.g., truncation,
# wrapping, folding, and doubling. We wouldn't want to cram all of this
# into the `ContinuousDistribution` class. Also, the order of the composition
# matters (e.g. truncate then shift/scale or vice versa). It's easier to
# accommodate different orders if the transformation is built up from
# components rather than all built into `ContinuousDistribution`.

def _shift_scale_distribution_function_2arg(func):
    def wrapped(self, x, y, *args, loc, scale, sign, **kwargs):
        item = func.__name__

        f = getattr(self._dist, item)

        # Obviously it's possible to get away with half of the work here.
        # Let's focus on correct results first and optimize later.
        xt = self._transform(x, loc, scale)
        yt = self._transform(y, loc, scale)
        fxy = f(xt, yt, *args, **kwargs)
        fyx = f(yt, xt, *args, **kwargs)
        return np.where(sign, fxy, fyx)[()]

    return wrapped


def _shift_scale_distribution_function(func):
    # c is for complementary
    citem = {'_logcdf_dispatch': '_logccdf_dispatch',
             '_cdf_dispatch': '_ccdf_dispatch',
             '_logccdf_dispatch': '_logcdf_dispatch',
             '_ccdf_dispatch': '_cdf_dispatch'}
    def wrapped(self, x, *args, loc, scale, sign, **kwargs):
        item = func.__name__

        f = getattr(self._dist, item)
        cf = getattr(self._dist, citem[item])

        # Obviously it's possible to get away with half of the work here.
        # Let's focus on correct results first and optimize later.
        xt = self._transform(x, loc, scale)
        fx = f(xt, *args, **kwargs)
        cfx = cf(xt, *args, **kwargs)
        return np.where(sign, fx, cfx)[()]

    return wrapped


def _shift_scale_inverse_function(func):
    citem = {'_ilogcdf_dispatch': '_ilogccdf_dispatch',
             '_icdf_dispatch': '_iccdf_dispatch',
             '_ilogccdf_dispatch': '_ilogcdf_dispatch',
             '_iccdf_dispatch': '_icdf_dispatch'}
    def wrapped(self, p, *args, loc, scale, sign, **kwargs):
        item = func.__name__

        f = getattr(self._dist, item)
        cf = getattr(self._dist, citem[item])

        # Obviously it's possible to get away with half of the work here.
        # Let's focus on correct results first and optimize later.
        fx =  self._itransform(f(p, *args, **kwargs), loc, scale)
        cfx = self._itransform(cf(p, *args, **kwargs), loc, scale)
        return np.where(sign, fx, cfx)[()]

    return wrapped


class ShiftedScaledDistribution(TransformedDistribution):
    """Distribution with a standard shift/scale transformation."""
    # Unclear whether infinite loc/scale will work reasonably in all cases
    _loc_domain = _RealDomain(endpoints=(-inf, inf), inclusive=(True, True))
    _loc_param = _RealParameter('loc', symbol=r'\mu',
                                domain=_loc_domain, typical=(1, 2))

    _scale_domain = _RealDomain(endpoints=(-inf, inf), inclusive=(True, True))
    _scale_param = _RealParameter('scale', symbol=r'\sigma',
                                  domain=_scale_domain, typical=(0.1, 10))

    _parameterizations = [_Parameterization(_loc_param, _scale_param),
                          _Parameterization(_loc_param),
                          _Parameterization(_scale_param)]

    def _process_parameters(self, loc=None, scale=None, **params):
        loc = loc if loc is not None else np.zeros_like(scale)[()]
        scale = scale if scale is not None else np.ones_like(loc)[()]
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
        a, b = self._dist._support(**params)
        a, b = self._itransform(a, loc, scale), self._itransform(b, loc, scale)
        return np.where(sign, a, b)[()], np.where(sign, b, a)[()]

    # Here, we override all the `_dispatch` methods rather than the public
    # methods or _function methods. Why not the public methods?
    # If we were to override the public methods, then other
    # TransformedDistribution classes (which could transform a
    # ShiftedScaledDistribution) would need to call the public methods of
    # ShiftedScaledDistribution, which would run the input validation again.
    # Why not the _function methods? For distributions that rely on the
    # default implementation of methods (e.g. `quadrature`, `inversion`),
    # the implementation would "see" the location and scale like other
    # distribution parameters, so they could affect the accuracy of the
    # calculations. I think it is cleaner if `loc` and `scale` do not affect
    # the underlying calculations at all.

    def _entropy_dispatch(self, *args, loc, scale, sign, **params):
        return (self._dist._entropy_dispatch(*args, **params)
                + np.log(abs(scale)))

    def _logentropy_dispatch(self, *args, loc, scale, sign, **params):
        lH0 = self._dist._logentropy_dispatch(*args, **params)
        lls = np.log(np.log(abs(scale))+0j)
        return special.logsumexp(np.broadcast_arrays(lH0, lls), axis=0)

    def _median_dispatch(self, *, method, loc, scale, sign, **params):
        raw = self._dist._median_dispatch(method=method, **params)
        return self._itransform(raw, loc, scale)

    def _mode_dispatch(self, *, method, loc, scale, sign, **params):
        raw = self._dist._mode_dispatch(method=method, **params)
        return self._itransform(raw, loc, scale)

    def _logpdf_dispatch(self, x, *args, loc, scale, sign, **params):
        x = self._transform(x, loc, scale)
        logpdf = self._dist._logpdf_dispatch(x, *args, **params)
        return logpdf - np.log(abs(scale))

    def _pdf_dispatch(self, x, *args, loc, scale, sign, **params):
        x = self._transform(x, loc, scale)
        pdf = self._dist._pdf_dispatch(x, *args, **params)
        return pdf / abs(scale)

    # Sorry about the magic. This is just a draft to show the behavior.
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
        res = (self._dist._moment_standardized_dispatch(
            order, methods=methods, **params))
        return None if res is None else res * np.sign(scale)**order

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

    def _sample_dispatch(self, sample_shape, full_shape, *,
                         method, rng, **params):
        rvs = self._dist._sample_dispatch(
            sample_shape, full_shape, method=method, rng=rng, **params)
        return self._itransform(rvs, **params)

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


class Mixture(_ProbabilityDistribution):
    r"""Representation of a mixture distribution.

    A mixture distribution is the distribution of a random variable
    defined in the following way: first, a random variable is selected
    from `components` according to the probabilities given by `weights`, then
    the selected random variable is realized.

    Parameters
    ----------
    components : sequence of `ContinuousDistribution`
        The underlying instances of `ContinuousDistribution`.
        All must have scalar shape parameters (if any); e.g., the `pdf` evaluated
        at a scalar argument must return a scalar.
    weights : sequence of floats
        The corresponding probabilities of selecting each random variable.
        Must be non-negative and sum to one.

    Attributes
    ----------
    components : sequence of `ContinuousDistribution`
        The underlying instances of `ContinuousDistribution`.
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

    """
    # Todo:
    # Fix Normal(mu=0.5).logentropy() runtime warning
    # Add support for array shapes, weights

    def _input_validation(self, components, weights):
        if len(components) == 0:
            message = ("`components` must contain at least one random variable.")
            raise ValueError(message)

        for var in components:
            # will generalize to other kinds of distributions when there
            # *are* other kinds of distributions
            if not isinstance(var, _ProbabilityDistribution):
                message = ("Each element of `components` must be an instance of "
                           "`ContinuousDistribution`.")
                raise ValueError(message)
            if not var._shape == ():
                message = "All elements of `components` must have scalar shapes."
                raise ValueError(message)

        if weights is None:
            return

        weights = np.asarray(weights)
        if weights.shape != (len(components),):
            message = "`components` and `weights` must have the same length."
            raise ValueError(message)

        if not np.issubdtype(weights.dtype, np.inexact):
            message = "`weights` must have floating point dtype."
            raise ValueError(message)

        if not np.isclose(np.sum(weights), 1.0):
            message = "`weights` must sum to 1.0."
            raise ValueError(message)

        if not np.all(weights >= 0):
            message = "All `weights` must be non-negative."
            raise ValueError(message)

        return components, weights

    def __init__(self, components, *, weights=None):
        components, weights = self._input_validation(components, weights)
        n = len(components)
        dtype = np.result_type(*(var._dtype for var in components))
        self._shape = np.broadcast_shapes(*(var._shape for var in components))
        self._dtype, self._components = dtype, components
        self._weights = np.full(n, 1/n, dtype=dtype) if weights is None else weights
        self.validation_policy = None

    @property
    def components(self):
        return list(self._components)

    @property
    def weights(self):
        return self._weights.copy()

    def _full(self, val, *args):
        args = [np.asarray(arg) for arg in args]
        dtype = np.result_type(self._dtype, *(arg.dtype for arg in args))
        shape = np.broadcast_shapes(self._shape, *(arg.shape for arg in args))
        return np.full(shape, val, dtype=dtype)

    def _sum(self, fun, *args):
        out = self._full(0, *args)
        for var, weight in zip(self._components, self._weights):
            out += getattr(var, fun)(*args) * weight
        return out

    def _logsum(self, fun, *args):
        out = self._full(-np.inf, *args)
        for var, log_weight in zip(self._components, np.log(self._weights)):
            np.logaddexp(out, getattr(var, fun)(*args) + log_weight, out=out)
        return out

    def support(self):
        a = self._full(np.inf)
        b = self._full(-np.inf)
        for var in self._components:
            a = np.minimum(a, var.support()[0])
            b = np.maximum(b, var.support()[1])
        return a, b

    def _raise_if_method(self, method):
        if method is not None:
            raise NotImplementedError("`method` not implemented for this distribution.")

    def logentropy(self, *, method=None):
        self._raise_if_method(method)
        def log_integrand(x):
            # `x` passed by `_tanhsinh` will be of complex dtype because
            # `log_integrand` returns complex values, but the imaginary
            # component is always zero. Extract the real part because
            # `logpdf` uses `logaddexp`, which fails for complex input.
            return self.logpdf(x.real) + np.log(self.logpdf(x.real) + 0j)

        res = _tanhsinh(log_integrand, *self.support(), log=True).integral
        return _log_real_standardize(res + np.pi*1j)

    def entropy(self, *, method=None):
        self._raise_if_method(method)
        return _tanhsinh(lambda x: -self.pdf(x) * self.logpdf(x),
                         *self.support()).integral

    def mode(self, *, method=None):
        self._raise_if_method(method)
        a, b = self.support()
        def f(x): return -self.pdf(x)
        res = _bracket_minimum(f, 1., xmin=a, xmax=b)
        res = _chandrupatla_minimize(f, res.xl, res.xm, res.xr)
        return res.x

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
        order = _ProbabilityDistribution._validate_order_kind(self, order, kind, kinds)
        moment_kind = kinds[kind]
        return moment_kind(order)

    def _moment_raw(self, order):
        out = self._full(0)
        for var, weight in zip(self._components, self._weights):
            out += var.moment(order, kind='raw') * weight
        return out

    def _moment_central(self, order):
        order = int(order)
        out = self._full(0)
        for var, weight in zip(self._components, self._weights):
            moment_as = [var.moment(order, kind='central')
                         for order in range(order + 1)]
            a, b = var.mean(), self.mean()
            moment = var._moment_transform_center(order, moment_as, a, b)
            out += moment * weight
        return out

    def _moment_standardized(self, order):
        return self._moment_central(order) / self.standard_deviation()**order

    def pdf(self, x, /, *, method=None):
        self._raise_if_method(method)
        return self._sum('pdf', x)

    def logpdf(self, x, /, *, method=None):
        self._raise_if_method(method)
        return self._logsum('logpdf', x)

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
        xmin, xmax = self.support()
        fun = getattr(self, fun)
        f = lambda x, p: fun(x) - p  # noqa: E731 is silly
        res = _bracket_root(f, xl0=self.mean(), xmin=xmin, xmax=xmax, args=(p,))
        return _chandrupatla(f, a=res.xl, b=res.xr, args=(p,)).x

    def icdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('cdf', p)

    def iccdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('ccdf', p)

    def ilogcdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('logcdf', p)

    def ilogccdf(self, p, /, *, method=None):
        self._raise_if_method(method)
        return self._invert('logccdf', p)

    def sample(self, shape=(), *, rng=None, method=None):
        self._raise_if_method(method)
        rng = np.random.default_rng(rng)
        size = np.prod(np.atleast_1d(shape))
        ns = rng.multinomial(size, self._weights)
        x = [var.sample(shape=n, rng=rng) for n, var in zip(ns, self._components)]
        x = np.reshape(rng.permuted(np.concatenate(x)), shape)
        return x[()]
