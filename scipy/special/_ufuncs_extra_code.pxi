cimport scipy.special._ufuncs_cxx

def errprint(inflag=None):
    """
    errprint(inflag=None)

    Set or return the error printing flag for special functions.

    Parameters
    ----------
    inflag : bool, optional
        Whether warnings concerning evaluation of special functions in
        ``scipy.special`` are shown. If omitted, no change is made to
        the current setting.

    Returns
    -------
    old_flag : bool
        Previous value of the error flag

    Examples
    --------
    Turn on error printing.

    >>> import warnings
    >>> import scipy.special as sc
    >>> sc.bdtr(-1, 10, 0.3)
    nan
    >>> sc.errprint(True)
    False
    >>> with warnings.catch_warnings(record=True) as w:
    ...     sc.bdtr(-1, 10, 0.3)
    ...
    nan
    >>> len(w)
    1
    >>> w[0].message
    SpecialFunctionWarning('scipy.special/bdtr: domain error',)

    """
    if inflag is not None:
        scipy.special._ufuncs_cxx._set_errprint(int(bool(inflag)))
        return bool(sf_error.set_print(int(bool(inflag))))
    else:
        return bool(sf_error.get_print())
