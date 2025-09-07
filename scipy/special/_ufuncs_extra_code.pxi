cimport scipy.special._ufuncs_cxx
cimport scipy.special._ellip_harm_2
import scipy.special._special_ufuncs
import scipy.special._gufuncs
import numpy as np


_sf_error_code_map = {
    # skip 'ok'
    'singular': 1,
    'underflow': 2,
    'overflow': 3,
    'slow': 4,
    'loss': 5,
    'no_result': 6,
    'domain': 7,
    'arg': 8,
    'other': 9,
    'memory': 10
}

_sf_error_action_map = {
    'ignore': 0,
    'warn': 1,
    'raise': 2,
    0: 'ignore',
    1: 'warn',
    2: 'raise'
}


def geterr():
    """Get the current way of handling special-function errors.

    Returns
    -------
    err : dict
        A dictionary with keys "singular", "underflow", "overflow",
        "slow", "loss", "no_result", "domain", "arg", and "other",
        whose values are from the strings "ignore", "warn", and
        "raise". The keys represent possible special-function errors,
        and the values define how these errors are handled.

    See Also
    --------
    seterr : set how special-function errors are handled
    errstate : context manager for special-function error handling
    numpy.geterr : similar numpy function for floating-point errors

    Notes
    -----
    For complete documentation of the types of special-function errors
    and treatment options, see `seterr`.

    Examples
    --------
    By default all errors are ignored.

    >>> import scipy.special as sc
    >>> for key, value in sorted(sc.geterr().items()):
    ...     print(f'{key}: {value}')
    ...
    arg: ignore
    domain: ignore
    loss: ignore
    memory: raise
    no_result: ignore
    other: ignore
    overflow: ignore
    singular: ignore
    slow: ignore
    underflow: ignore

    """
    err = {}
    for key, code in _sf_error_code_map.items():
        action = sf_error.get_action(code)
        err[key] = _sf_error_action_map[action]
    return err


def seterr(**kwargs):
    """Set how special-function errors are handled.

    Parameters
    ----------
    all : {'ignore', 'warn' 'raise'}, optional
        Set treatment for all type of special-function errors at
        once. The options are:

        - 'ignore' Take no action when the error occurs
        - 'warn' Print a `SpecialFunctionWarning` when the error
          occurs (via the Python `warnings` module)
        - 'raise' Raise a `SpecialFunctionError` when the error
          occurs.

        The default is to not change the current behavior. If
        behaviors for additional categories of special-function errors
        are specified, then ``all`` is applied first, followed by the
        additional categories.
    singular : {'ignore', 'warn', 'raise'}, optional
        Treatment for singularities.
    underflow : {'ignore', 'warn', 'raise'}, optional
        Treatment for underflow.
    overflow : {'ignore', 'warn', 'raise'}, optional
        Treatment for overflow.
    slow : {'ignore', 'warn', 'raise'}, optional
        Treatment for slow convergence.
    loss : {'ignore', 'warn', 'raise'}, optional
        Treatment for loss of accuracy.
    no_result : {'ignore', 'warn', 'raise'}, optional
        Treatment for failing to find a result.
    domain : {'ignore', 'warn', 'raise'}, optional
        Treatment for an invalid argument to a function.
    arg : {'ignore', 'warn', 'raise'}, optional
        Treatment for an invalid parameter to a function.
    other : {'ignore', 'warn', 'raise'}, optional
        Treatment for an unknown error.

    Returns
    -------
    olderr : dict
        Dictionary containing the old settings.

    See Also
    --------
    geterr : get the current way of handling special-function errors
    errstate : context manager for special-function error handling
    numpy.seterr : similar numpy function for floating-point errors

    Examples
    --------
    >>> import scipy.special as sc
    >>> from pytest import raises
    >>> sc.gammaln(0)
    inf
    >>> olderr = sc.seterr(singular='raise')
    >>> with raises(sc.SpecialFunctionError):
    ...     sc.gammaln(0)
    ...
    >>> _ = sc.seterr(**olderr)

    We can also raise for every category except one.

    >>> olderr = sc.seterr(all='raise', singular='ignore')
    >>> sc.gammaln(0)
    inf
    >>> with raises(sc.SpecialFunctionError):
    ...     sc.spence(-1)
    ...
    >>> _ = sc.seterr(**olderr)

    """
    olderr = geterr()

    if 'all' in kwargs.keys():
        action = kwargs.pop('all')
        newkwargs = {key: action for key in _sf_error_code_map.keys()}
        for key, value in kwargs.items():
            newkwargs[key] = value
        kwargs = newkwargs

    for error, action in kwargs.items():
        action = _sf_error_action_map[action]
        code = _sf_error_code_map[error]
        # Error handling state must be set for all relevant
        # extension modules in synchrony, since each carries
        # a separate copy of this state.
        _set_action(code, action)
        scipy.special._ufuncs_cxx._set_action(code, action)
        scipy.special._special_ufuncs._set_action(code, action)
        scipy.special._gufuncs._set_action(code, action)
        scipy.special._ellip_harm_2._set_action(code, action)

    return olderr


class errstate:
    """Context manager for special-function error handling.

    Using an instance of `errstate` as a context manager allows
    statements in that context to execute with a known error handling
    behavior. Upon entering the context the error handling is set with
    `seterr`, and upon exiting it is restored to what it was before.

    Parameters
    ----------
    kwargs : {all, singular, underflow, overflow, slow, loss, no_result, domain, arg, other}
        Keyword arguments. The valid keywords are possible
        special-function errors. Each keyword should have a string
        value that defines the treatment for the particular type of
        error. Values must be 'ignore', 'warn', or 'other'. See
        `seterr` for details.

    See Also
    --------
    geterr : get the current way of handling special-function errors
    seterr : set how special-function errors are handled
    numpy.errstate : similar numpy function for floating-point errors

    Examples
    --------
    >>> import scipy.special as sc
    >>> from pytest import raises
    >>> sc.gammaln(0)
    inf
    >>> with sc.errstate(singular='raise'):
    ...     with raises(sc.SpecialFunctionError):
    ...         sc.gammaln(0)
    ...
    >>> sc.gammaln(0)
    inf

    We can also raise on every category except one.

    >>> with sc.errstate(all='raise', singular='ignore'):
    ...     sc.gammaln(0)
    ...     with raises(sc.SpecialFunctionError):
    ...         sc.spence(-1)
    ...
    inf

    """
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def __enter__(self):
        self.oldstate = seterr(**self.kwargs)

    def __exit__(self, exc_type, exc_value, traceback):
        seterr(**self.oldstate)
