from typing import (
    Callable,
    Iterable,
    Dict,
    Tuple,
    Any,
    Set,
    Optional,
    Type,
    Union,
    List,
)
import inspect
import functools
from . import _uarray  # type: ignore
import copyreg
import atexit
import pickle

ArgumentExtractorType = Callable[..., Tuple["Dispatchable", ...]]
ArgumentReplacerType = Callable[[Tuple, Dict, Tuple], Tuple[Tuple, Dict]]

from ._uarray import (
    BackendNotImplementedError,
    _Function,
    _SkipBackendContext,
    _SetBackendContext,
)


def unpickle_function(mod_name, qname):
    import importlib

    try:
        module = importlib.import_module(mod_name)
        func = getattr(module, qname)
        return func
    except (ImportError, AttributeError) as e:
        from pickle import UnpicklingError

        raise UnpicklingError from e


def pickle_function(func):
    mod_name = getattr(func, "__module__", None)
    qname = getattr(func, "__qualname__", None)

    try:
        test = unpickle_function(mod_name, qname)
    except pickle.UnpicklingError:
        test = None

    if test is not func:
        raise pickle.PicklingError(
            "Can't pickle {}: it's not the same object as {}".format(func, test)
        )

    return unpickle_function, (mod_name, qname)


copyreg.pickle(_Function, pickle_function)
atexit.register(_uarray.clear_all_globals)


def create_multimethod(*args, **kwargs):
    """
    Creates a decorator for generating multimethods.

    This function creates a decorator that can be used with an argument
    extractor in order to generate a multimethod. Other than for the
    argument extractor, all arguments are passed on to
    :obj:`generate_multimethod`.

    See Also
    --------
    generate_multimethod
        Generates a multimethod.
    """

    def wrapper(a):
        return generate_multimethod(a, *args, **kwargs)

    return wrapper


def generate_multimethod(
    argument_extractor: ArgumentExtractorType,
    argument_replacer: ArgumentReplacerType,
    domain: str,
    default: Optional[Callable] = None,
):
    """
    Generates a multimethod.

    Parameters
    ----------
    argument_extractor : ArgumentExtractorType
        A callable which extracts the dispatchable arguments. Extracted arguments
        should be marked by the :obj:`Dispatchable` class. It has the same signature
        as the desired multimethod.
    argument_replacer : ArgumentReplacerType
        A callable with the signature (args, kwargs, dispatchables), which should also
        return an (args, kwargs) pair with the dispatchables replaced inside the args/kwargs.
    domain : str
        A string value indicating the domain of this multimethod.
    default: Optional[Callable], optional
        The default implementation of this multimethod, where ``None`` (the default) specifies
        there is no default implementation.

    Examples
    --------
    In this example, ``a`` is to be dispatched over, so we return it, while marking it as an ``int``.
    The trailing comma is needed because the args have to be returned as an iterable.

    >>> def override_me(a, b):
    ...   return Dispatchable(a, int),

    Next, we define the argument replacer that replaces the dispatchables inside args/kwargs with the
    supplied ones.

    >>> def override_replacer(args, kwargs, dispatchables):
    ...     return (dispatchables[0], args[1]), {}

    Next, we define the multimethod.

    >>> overridden_me = generate_multimethod(
    ...     override_me, override_replacer, "ua_examples"
    ... )

    Notice that there's no default implementation, unless you supply one.

    >>> overridden_me(1, "a")
    Traceback (most recent call last):
        ...
    uarray.backend.BackendNotImplementedError: ...
    >>> overridden_me2 = generate_multimethod(
    ...     override_me, override_replacer, "ua_examples", default=lambda x, y: (x, y)
    ... )
    >>> overridden_me2(1, "a")
    (1, 'a')

    See Also
    --------
    uarray
        See the module documentation for how to override the method by creating backends.
    """
    kw_defaults, arg_defaults, opts = get_defaults(argument_extractor)
    ua_func = _Function(
        argument_extractor,
        argument_replacer,
        domain,
        arg_defaults,
        kw_defaults,
        default,
    )

    return functools.update_wrapper(ua_func, argument_extractor)


def set_backend(backend, *args, **kwargs):
    """
    A context manager that sets the preferred backend. Uses :obj:`BackendOptions` to create
    the options, and sets the backend with the preferred options.

    Parameters
    ----------
    backend
        The backend to set.

    See Also
    --------
    BackendOptions: The backend plus options.
    skip_backend: A context manager that allows skipping of backends.
    """
    return _SetBackendContext(backend, *args, **kwargs)


def skip_backend(backend):
    """
    A context manager that allows one to skip a given backend from processing
    entirely. This allows one to use another backend's code in a library that
    is also a consumer of the same backend.

    Parameters
    ----------
    backend
        The backend to skip.

    See Also
    --------
    set_backend: A context manager that allows setting of backends.
    """
    return _SkipBackendContext(backend)


def get_defaults(f):
    sig = inspect.signature(f)
    kw_defaults = {}
    arg_defaults = []
    opts = set()
    for k, v in sig.parameters.items():
        if v.default is not inspect.Parameter.empty:
            kw_defaults[k] = v.default
        if v.kind in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        ):
            arg_defaults.append(v.default)
        opts.add(k)

    return kw_defaults, tuple(arg_defaults), opts


def set_global_backend(backend):
    """
    This utility method replaces the default backend for permanent use. It
    will be tried in the list of backends automatically, unless the
    ``only`` flag is set on a backend. This will be the first tried
    backend outside the :obj:`set_backend` context manager.

    Note that this method is not thread-safe.

    .. warning::
        We caution library authors against using this function in
        their code. We do *not* support this use-case. This function
        is meant to be used only by users themselves, or by a reference
        implementation, if one exists.

    Parameters
    ----------
    backend
        The backend to register.
    """
    _uarray.set_global_backend(backend)


def register_backend(backend):
    """
    This utility method sets registers backend for permanent use. It
    will be tried in the list of backends automatically, unless the
    ``only`` flag is set on a backend.

    Note that this method is not thread-safe.

    Parameters
    ----------
    backend
        The backend to register.
    """
    _uarray.register_backend(backend)


class Dispatchable:
    """
    A utility class which marks an argument with a specific dispatch type.


    Attributes
    ----------
    value
        The value of the Dispatchable.

    type
        The type of the Dispatchable.

    Examples
    --------
    >>> x = Dispatchable(1, str)
    >>> x
    <Dispatchable: type=<class 'str'>, value=1>

    See Also
    --------
    all_of_type
        Marks all unmarked parameters of a function.

    mark_as
        Allows one to create a utility function to mark as a given type.
    """

    def __init__(self, value, dispatch_type, coercible=True):
        self.value = value
        self.type = dispatch_type
        self.coercible = coercible

    def __getitem__(self, index):
        return (self.type, self.value)[index]

    def __str__(self):
        return (
            "<{0}: type={1!r}, value={2!r}>".format(
                type(self).__name__, self.type, self.value
            )
        )

    __repr__ = __str__


def mark_as(dispatch_type):
    """
    Creates a utility function to mark something as a specific type.

    Examples
    --------
    >>> mark_int = mark_as(int)
    >>> mark_int(1)
    <Dispatchable: type=<class 'int'>, value=1>
    """
    return functools.partial(Dispatchable, dispatch_type=dispatch_type)


def all_of_type(arg_type):
    """
    Marks all unmarked arguments as a given type.

    Examples
    --------
    >>> @all_of_type(str)
    ... def f(a, b):
    ...     return a, Dispatchable(b, int)
    >>> f('a', 1)
    (<Dispatchable: type=<class 'str'>, value='a'>, <Dispatchable: type=<class 'int'>, value=1>)
    """

    def outer(func):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            extracted_args = func(*args, **kwargs)
            return tuple(
                Dispatchable(arg, arg_type)
                if not isinstance(arg, Dispatchable)
                else arg
                for arg in extracted_args
            )

        return inner

    return outer


def wrap_single_convertor(convert_single):
    """
    Wraps a ``__ua_convert__`` defined for a single element to all elements.
    If any of them return ``NotImplemented``, the operation is assumed to be
    undefined.

    Accepts a signature of (value, type, coerce).
    """

    @functools.wraps(convert_single)
    def __ua_convert__(dispatchables, coerce):
        converted = []
        for d in dispatchables:
            c = convert_single(d.value, d.type, coerce and d.coercible)

            if c is NotImplemented:
                return NotImplemented

            converted.append(c)

        return converted

    return __ua_convert__
