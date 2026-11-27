"""Private helper functions and classes used by the univariate distribution classes.

Ideally, all private functions and classes used by the main distribution classes
go here.
"""
from abc import ABC, abstractmethod
import functools
import inspect
import math
from math import inf
from types import GenericAlias

import numpy as np
from scipy._lib._array_api import (xp_capabilities, xp_promote, xp_copy,  is_numpy,
                                   xp_device, xp_result_device, xp_compat_namespace)
from scipy._external import array_api_extra as xpx
from scipy import special


# Could add other policies for broadcasting and edge/out-of-bounds case handling
# For instance, when edge case handling is known not to be needed, it's much
# faster to turn it off, but it might still be nice to have array conversion
# and shaping done so the user doesn't need to be so careful.
_SKIP_ALL = "skip_all"
# Other cache policies would be useful, too.
_NO_CACHE = "no_cache"


rv_capabilities = xp_capabilities(
    cpu_only=True,
    exceptions=['cupy'],
    skip_backends=[('array_api_strict', "fancy indexing in EIM"),
                   ("jax.numpy", 'mutation'),
                   ('dask.array', 'no attempt')],
    method_capabilities={
            "plot": dict(np_only=True),
    }
)


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
        defined symbolically or through a callable.
    __str__()
        Returns a text representation of the domain (e.g. ``[0, b)``).
        Used for generating documentation.

    """
    symbols = {math.inf: r"\infty", -math.inf: r"-\infty",
               math.pi: r"\pi", -math.pi: r"-\pi", 2*math.pi: r"2\pi"}

    # generic type compatibility with scipy-stubs
    __class_getitem__: classmethod = classmethod(GenericAlias)

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


class _Interval(_Domain):
    r""" Representation of an interval defined by two endpoints.

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
    endpoints : 2-tuple of float(s) and/or str(s) and/or callable(s).
        A tuple with two values. Each may be either a float (the numerical
        value of the endpoints of the domain), a string (the name of the
        parameters that will define the endpoint), or a callable taking the
        parameters used to define the endpoints of the domain as keyword only
        arguments and returning a numerical value for the endpoint.
    inclusive : 2-tuple of bools
        A tuple with two boolean values; each indicates whether the
        corresponding endpoint is included within the domain or not.

    Methods
    -------
    define_parameters(*parameters)
        Records any parameters used to define the endpoints of the domain
    get_numerical_endpoints(parameter_values)
        Gets the numerical values of the domain endpoints, which may have been
        defined symbolically or through a callable.
    contains(item, parameter_values)
        Determines whether the argument is contained within the domain
    draw(size, rng, proportions, parameter_values)
        Draws random values based on the domain.

    """
    def __init__(self, endpoints=(-inf, inf), inclusive=(False, False)):
        self.symbols = super().symbols.copy()
        self.endpoints = endpoints
        self.inclusive = inclusive

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

    def get_numerical_endpoints(self, parameter_values, xp, device=None):
        r""" Get the numerical values of the domain endpoints.

        Domain endpoints may be defined symbolically or through a callable.
        This returns numerical values of the endpoints given numerical values for
        any variables.

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
        a, b = self.endpoints
        # If `a` (`b`) is a string - the name of the parameter that defines
        # the endpoint of the domain - then corresponding numerical values
        # will be found in the `parameter_values` dictionary.
        # If a callable, it will be executed with `parameter_values` passed as
        # keyword arguments, and it will return the numerical values.
        # Otherwise, it is itself the array of numerical values of the endpoint.
        try:
            if callable(a):
                a = a(**parameter_values)
            else:
                a = parameter_values.get(a, a)
            if callable(b):
                b = b(**parameter_values)
            else:
                b = parameter_values.get(b, b)
        except TypeError as e:
            message = ("The endpoints of the distribution are defined by "
                       "parameters, but their values were not provided. When "
                       f"using a private method of {self.__class__}, pass "
                       "all required distribution parameters as keyword "
                       "arguments.")
            raise TypeError(message) from e
        # Floating point types are used for even integer parameters.
        # Convert to float here to ensure consistency throughout framework.
        a, b = xp.asarray(a, device=device), xp.asarray(b, device=device)
        a, b = xp_promote(a, b, force_floating=True, xp=xp)
        return a, b

    def contains(self, item, parameter_values=None, *, xp):
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
        device = xp_device(item)
        a, b = self.get_numerical_endpoints(parameter_values, xp=xp, device=device)
        left_inclusive, right_inclusive = self.inclusive

        in_left = item >= a if left_inclusive else item > a
        in_right = item <= b if right_inclusive else item < b
        return in_left & in_right

    def draw(self, n, type_, min, max, squeezed_base_shape, rng=None, *, xp):
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
        squeezed_base_shape : tuple of ints
            See _RealParameter.draw.
        rng : np.Generator
            The Generator used for drawing random values.

        """
        rng = np.random.default_rng(rng)
        dtype = xpx.default_dtype(xp)
        device = xp_result_device(min, max)
        ddd = dict(dtype=dtype, device=device)

        # def ints(*args, **kwargs): return rng.integers(*args, **kwargs, endpoint=True)
        # uniform = rng.uniform if isinstance(self, _RealInterval) else ints
        def uniform(a, b, size):
            is_integer = isinstance(self, _IntegerInterval)
            offset = 1 if is_integer else 0
            raw = a + (b - a + offset) * xp.asarray(rng.random(size=size), **ddd)
            return xp.floor(raw) if is_integer else raw

        # get copies of min and max with no nans so that uniform doesn't fail
        min_nn, max_nn = xp_copy(min), xp_copy(max)
        i = xp.isnan(min_nn) | xp.isnan(max_nn)
        min_nn[i] = 0
        max_nn[i] = 1

        shape = (n,) + squeezed_base_shape
        if is_numpy(xp) and np.__version__ < "2.1":
            xp = xp_compat_namespace(xp)

        if type_ == 'in':
            z = uniform(min_nn, max_nn, size=shape)

        elif type_ == 'on':
            z_on_shape = shape
            z = xp.ones(z_on_shape, **ddd)
            i = xp.asarray(rng.random(size=n) < 0.5, device=device)
            z[i] = min
            z[~i] = max

        elif type_ == 'out':
            z = min_nn - uniform(1, 5, size=shape)   # 1, 5 is arbitrary; we just want
            zr = max_nn + uniform(1, 5, size=shape)  # some numbers outside domain
            i = xp.asarray(rng.random(size=n) < 0.5, device=device)
            z[i] = zr[i]

        elif type_ == 'nan':
            z = xp.full(shape, xp.nan, **ddd)

        return z


class _RealInterval(_Interval):
    r""" Represents a simply-connected subset of the real line; i.e., an interval

    Completes the implementation of the `_Interval` class for intervals
    on the real line.

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
    """

    def __str__(self):
        a, b = self.endpoints
        a, b = self._get_endpoint_str(a, "f1"), self._get_endpoint_str(b, "f2")
        left_inclusive, right_inclusive = self.inclusive
        left = "[" if left_inclusive else "("
        right = "]" if right_inclusive else ")"

        return f"{left}{a}, {b}{right}"

    def _get_endpoint_str(self, endpoint, funcname):
        if callable(endpoint):
            if endpoint.__doc__ is not None:
                return endpoint.__doc__
            params = inspect.signature(endpoint).parameters.values()
            params = [
                p.name for p in params if p.kind == inspect.Parameter.KEYWORD_ONLY
            ]
            return f"{funcname}({','.join(params)})"
        return self.symbols.get(endpoint, f"{endpoint}")


class _IntegerInterval(_Interval):
    r""" Represents an interval of integers

    Completes the implementation of the `_Interval` class for simple
    domains on the integers.

    Methods
    -------
    define_parameters(*parameters)
        (Inherited) Records any parameters used to define the endpoints of the
        domain.
    get_numerical_endpoints(parameter_values)
        (Inherited) Gets the numerical values of the domain endpoints, which
        may have been defined symbolically.
    contains(item, parameter_values)
        (Overridden) Determines whether the argument is contained within the
        domain
    draw(n, type_, min, max, squeezed_base_shape, rng=None)
        (Inherited) Draws random values based on the domain.
    __str__()
        Returns a string representation of the domain, e.g. "{a, a+1, ..., b-1, b}".

    """
    def contains(self, item, parameter_values=None, *, xp):
        super_contains = super().contains(item, parameter_values, xp=xp)
        integral = (item == xp.round(item))
        return super_contains & integral

    def __str__(self):
        a, b = self.endpoints
        a = self.symbols.get(a, a)
        b = self.symbols.get(b, b)

        a_str, b_str = isinstance(a, str), isinstance(b, str)
        a_inf = a == r"-\infty" if a_str else a == -math.inf
        b_inf = b == r"\infty" if b_str else b == math.inf

        # This doesn't work well for cases where ``a`` is floating point
        # number large enough that ``nextafter(a, inf) > a + 1``, and
        # similarly for ``b`` and nextafter(b, -inf). There may not be any
        # distributions fit for SciPy where we would actually need to handle these
        # cases though.
        ap1 = f"{a} + 1" if a_str else f"{a + 1}"
        bm1 = f"{b} - 1" if b_str else f"{b - 1}"

        if not a_str and not b_str:
            gap = b - a
            if gap == 3:
                return f"\\{{{a}, {ap1}, {bm1}, {b}\\}}"
            if gap == 2:
                return f"\\{{{a}, {ap1}, {b}\\}}"
            if gap == 1:
                return f"\\{{{a}, {b}\\}}"
            if gap == 0:
                return f"\\{{{a}\\}}"

        if not a_inf and b_inf:
            ap2 = f"{a} + 2" if a_str else f"{a + 2}"
            return f"\\{{{a}, {ap1}, {ap2}, ...\\}}"
        if a_inf and not b_inf:
            bm2 = f"{b} - 2" if b_str else f"{b - 2}"
            return f"\\{{{b}, {bm1}, {bm2}, ...\\}}"
        if a_inf and b_inf:
            return "\\{..., -2, -1, 0, 1, 2, ...\\}"

        return f"\\{{{a}, {ap1}, ..., {bm1}, {b}\\}}"


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
    docstring: str
        A string that describes the parameter.

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

    # generic type compatibility with scipy-stubs
    __class_getitem__: classmethod = classmethod(GenericAlias)

    def __init__(self, name, *, domain, symbol=None, typical=None, docstring=None):
        self.name = name
        self.symbol = symbol or name
        self.domain = domain
        if typical is not None and not isinstance(typical, _Domain):
            typical = domain.__class__(typical)
        self.typical = typical or domain
        self.__doc__ = self.docstring = docstring

    def __str__(self):
        r""" String representation of the parameter for use in documentation."""
        return f"`{self.name}` for :math:`{self.symbol} \\in {str(self.domain)}`"

    def draw(self, size=None, *, rng=None, region='domain', proportions=None,
             parameter_values=None, xp):
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

        a, b = domain.get_numerical_endpoints(parameter_values, xp=xp)
        a, b = xp.broadcast_arrays(a, b)

        # work with shapes in NumPy
        base_shape = a.shape
        extended_shape = np.broadcast_shapes(size, base_shape)
        n_extended = np.prod(extended_shape)
        n_base = np.prod(base_shape)
        n = int(n_extended / n_base) if n_extended else 0

        rng = np.random.default_rng(rng)
        pvals = proportions / np.sum(proportions)
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

        min = xp.asarray(a.squeeze())  # This is not array-API compliant, but
        max = xp.asarray(b.squeeze())  # it will work with supported backends
        squeezed_base_shape = max.shape

        if region == 'typical':
            typical = self.typical
            a, b = typical.get_numerical_endpoints(parameter_values, xp=xp)
            a, b = xp.broadcast_arrays(a, b)
            min_here = xp.asarray(a.squeeze())
            max_here = xp.asarray(b.squeeze())
            z_in = typical.draw(n_in, 'in', min_here, max_here, squeezed_base_shape,
                                rng=rng, xp=xp)
        else:
            z_in = domain.draw(n_in, 'in', min, max, squeezed_base_shape,
                               rng=rng, xp=xp)

        z_on = domain.draw(n_on, 'on', min, max, squeezed_base_shape, rng=rng, xp=xp)
        z_out = domain.draw(n_out, 'out', min, max, squeezed_base_shape, rng=rng, xp=xp)
        z_nan= domain.draw(n_nan, 'nan', min, max, squeezed_base_shape, rng=rng, xp=xp)

        z = xp.concat((z_in, z_on, z_out, z_nan), axis=0)
        # z = xp.asarray(rng.permuted(z, axis=0))  # TODO: fix permutation

        z = xp.reshape(z, tuple(shape_expansion) + squeezed_base_shape)
        z = xp.moveaxis(z, new_base_singletons, tuple(base_singletons))
        return z

    @abstractmethod
    def validate(self, arr):
        raise NotImplementedError()


class _RealParameter(_Parameter):
    r""" Represents a real-valued parameter.

    Implements the remaining methods of _Parameter for real parameters.
    All attributes are inherited.

    """
    def validate(self, arr, parameter_values, xp):
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
        valid : bool ndarray
            Logical array indicating which elements are valid (True) and
            which are not (False). The arrays of all distribution parameters
            will be broadcasted, and elements for which any parameter value
            does not meet the requirements will be replaced with NaN.

        """
        if not xp.isdtype(arr.dtype, 'real floating'):
            message = f"Parameter `{self.name}` must be of real dtype."
            raise TypeError(message)
        return self.domain.contains(arr, parameter_values, xp=xp)


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

    def get_param_doc(self, param):
        return self.parameters[param].__doc__ if param in self.parameters else None

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

    def validation(self, parameter_values, xp):
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
        """
        all_valid = True
        for name, arr in parameter_values.items():
            parameter = self.parameters[name]
            valid = parameter.validate(arr, parameter_values, xp=xp)
            all_valid = all_valid & valid
        return all_valid

    def __str__(self):
        r"""Returns a string representation of the parameterization."""
        messages = [str(param) for name, param in self.parameters.items()]
        return ", ".join(messages)

    def draw(self, sizes=None, rng=None, proportions=None, region='domain', *, xp):
        r"""Draw random values of all parameters for use in testing.

        Parameters
        ----------
        sizes : iterable of shape tuples
            The size of the array to be generated for each parameter in the
            parameterization. Note that the order of sizes is arbitrary; the
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
        # depend on others. If the relationships is simple (e.g. a < b < c),
        # we can draw values in order a, b, c.
        parameter_values = {}

        if sizes is None or not len(sizes) or not np.iterable(sizes[0]):
            sizes = [sizes]*len(self.parameters)

        for size, param in zip(sizes, self.parameters.values()):
            parameter_values[param.name] = param.draw(
                size, rng=rng, proportions=proportions,
                parameter_values=parameter_values,
                region=region, xp=xp
            )

        return parameter_values


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
        xp = self._xp
        if self.validation_policy == _SKIP_ALL:
            return f(self, *args, **kwargs)

        res = f(self, *args, **kwargs)
        if res is None:
            # message could be more appropriate
            raise NotImplementedError(self._not_implemented)

        # TODO: this should convert to desired dtype, and only once
        # In principle, this could lose precision if default dtype is < self._dtype
        res = xp.asarray(res, device=self._device)
        needs_copy = False
        dtype = res.dtype

        if dtype != self._dtype:  # this won't work for logmoments (complex)
            dtype = xp.result_type(dtype, self._dtype)
            needs_copy = True

        if res.shape != self._shape:  # faster to check first
            res = xp.broadcast_to(res, self._shape)
            needs_copy = needs_copy or self._any_invalid

        if needs_copy:
            res = xp.astype(res, dtype, copy=True)

        # Ensure invalid shape parameters produce NaN result
        if self._any_invalid:
            res[self._invalid] = xp.nan

        return res[()]

    return filtered


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
            if func_name != '_sample_dispatch' and self.cache_policy != _NO_CACHE:
                self._method_cache[func_name] = method

        try:
            return method(*args, **kwargs)
        except KeyError as e:
            raise NotImplementedError(self._not_implemented) from e

    return wrapped


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
        xp = self._xp
        func_name = f.__name__

        low, high = self.support()
        x, y, low, high = xp.broadcast_arrays(x, y, low, high)
        dtype = xp.result_type(x.dtype, y.dtype, self._dtype)
        # yes, copy to avoid modifying input arrays
        x, y = xp.astype(x, dtype, copy=True), xp.astype(y, dtype, copy=True)

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
            res = xp.clip(res, 0., 1.)
        else:
            res = xp.where(xp.real(res) > 0, 0., res)

        # Transform the result to account for swapped argument order
        res = xp.asarray(res)
        if func_name == '_cdf2':
            res[i_swap] *= -1.
        elif func_name == '_ccdf2':
            res[i_swap] *= -1
            res[i_swap] += 2.
        elif func_name == '_logcdf2':
            res = xp.asarray(res + 0j) if xp.any(i_swap) else res
            res[i_swap] = res[i_swap] + xp.pi*1j
        else:
            # res[i_swap] is always positive and less than 1, so it's
            # safe to ensure that the result is real
            two = xp.asarray(2, dtype=res.dtype, device=xp_device(res))
            res[i_swap] = xp.real(_logexpxmexpy(xp.log(two), res[i_swap], xp=xp))
        return res[()]

    return wrapped


@xp_capabilities(skip_backends=[('dask.array', "dtype inference fails")])
def _logexpxmexpy(x, y, *, xp):
    """ Compute the log of the difference of the exponentials of two arguments.

    Avoids over/underflow, but does not prevent loss of precision otherwise.
    """
    return xpx.apply_where(
        x != y,
        (x, y),
        lambda x, y: special.logsumexp(xp.stack((x, y+xp.pi*1j)), axis=0),
        fill_value=-xp.inf, xp=xp
    )


def _guess_bracket(xmin, xmax, *, xp):
    a = xp.full_like(xmin, -1.0)
    b = xp.ones_like(xmax)

    i = xp.isfinite(xmin) & xp.isfinite(xmax)
    a[i] = xmin[i]
    b[i] = xmax[i]

    i = xp.isfinite(xmin) & ~xp.isfinite(xmax)
    a[i] = xmin[i]
    b[i] = xmin[i] + 1

    i = xp.isfinite(xmax) & ~xp.isfinite(xmin)
    a[i] = xmax[i] - 1
    b[i] = xmax[i]

    return a, b


def _log_real_standardize(x, *, xp):
    """Standardizes the (complex) logarithm of a real number.

    The logarithm of a real number may be represented by a complex number with
    imaginary part that is a multiple of pi*1j. Even multiples correspond with
    a positive real and odd multiples correspond with a negative real.

    Given a logarithm of a real number `x`, this function returns an equivalent
    representation in a standard form: the log of a positive real has imaginary
    part `0` and the log of a negative real has imaginary part `pi`.

    """
    shape = x.shape
    x = xpx.atleast_nd(x, ndim=2, xp=xp)
    real = xp.astype(xp.real(x), x.dtype)
    complex = xp.imag(x)
    y = real
    negative = xp.real(xp.exp(complex*1j)) < 0
    y[negative] = y[negative] + xp.pi * 1j
    return xp.reshape(y, shape)[()]


def _generate_domain_support(dist_family):
    n_parameterizations = len(dist_family._parameterizations)

    domain = f"\nfor :math:`x \\in {dist_family._variable.domain}`.\n"

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


def _shift_scale_distribution_function_2arg(func):
    def wrapped(self, x, y, *args, loc, scale, sign, method=None, **kwargs):
        item = func.__name__

        f = getattr(self._dist, item)

        xt = self._transform(x, loc, scale)
        yt = self._transform(y, loc, scale)
        res = xpx.apply_where(
            sign,
            (xt, yt, *args),
            lambda xt, yt, *args, **kwargs: f(xt, yt, *args, method=method, **kwargs),
            lambda xt, yt, *args, **kwargs: f(yt, xt, *args, method=method, **kwargs),
            kwargs=kwargs
        )
        return res[()]

    return wrapped


def _shift_scale_distribution_function(func):
    # c is for complementary
    citem = {'_logcdf_dispatch': '_logccdf_dispatch',
             '_cdf_dispatch': '_ccdf_dispatch',
             '_logccdf_dispatch': '_logcdf_dispatch',
             '_ccdf_dispatch': '_cdf_dispatch'}
    def wrapped(self, x, *args, loc, scale, sign, method=None, **kwargs):
        item = func.__name__

        f = getattr(self._dist, item)
        cf = getattr(self._dist, citem[item])

        xt = self._transform(x, loc, scale)
        res = xpx.apply_where(
            sign,
            (xt, *args),
            lambda xt, *args, **kwargs: f(xt, *args, method=method, **kwargs),
            lambda xt, *args, **kwargs: cf(xt, *args, method=method, **kwargs),
            kwargs=kwargs
        )
        return res[()]

    return wrapped


def _shift_scale_inverse_function(func):
    citem = {'_ilogcdf_dispatch': '_ilogccdf_dispatch',
             '_icdf_dispatch': '_iccdf_dispatch',
             '_ilogccdf_dispatch': '_ilogcdf_dispatch',
             '_iccdf_dispatch': '_icdf_dispatch'}
    def wrapped(self, p, *args, loc, scale, sign, method=None, **kwargs):
        item = func.__name__

        f = getattr(self._dist, item)
        cf = getattr(self._dist, citem[item])

        res = xpx.apply_where(
            sign,
            (p, *args),
            lambda p, *args, **kwargs: f(p, *args, method=method, **kwargs),
            lambda p, *args, **kwargs: cf(p, *args, method=method, **kwargs),
            kwargs=kwargs
        )
        return self._itransform(res, loc, scale)[()]

    return wrapped

