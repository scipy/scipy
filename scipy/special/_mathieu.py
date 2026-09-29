from scipy._lib.deprecation import _deprecated
from scipy.special._ufuncs import _mathieu_ce, _mathieu_se, _mathieu_cem, _mathieu_sem
from scipy.special._ufunc_tools import _make_ufunc_wrapper, _with_cache_optimization

__all__ = ["mathieu_cem", "mathieu_sem", "mathieu_ce", "mathieu_se"]

_mathieu_cem_doc = (
    r"""mathieu_cem(m, q, x, out=None)

    Even Mathieu function and its derivative.

    .. deprecated:: 2.0.0
        `mathieu_cem` is deprecated and will be removed in SciPy 2.2.0.
        Use `mathieu_ce` instead.

    Returns the even Mathieu function, ``ce_m(x, q)``, of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` expressed in radians.

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Derivative with respect to `x` expressed in radians.

    See Also
    --------
    mathieu_ce, mathieu_a, mathieu_b, mathieu_sem

    Notes
    -----
    The even Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (a_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`a_m` (calculated with `mathieu_a`)
    results in an even, periodic solution :math:`y(x)` with period 180 degrees
    (for even :math:`m`) or 360 degrees (for odd :math:`m`).

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function
    .. [2] Stuart Brorson, A New Implementation of the Mathieu Functions for SciPy.
           https://github.com/brorson/ScipyMathieuPaper

    Examples
    --------
    Plot even Mathieu functions of orders ``2`` and ``4``.

    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> m = np.asarray([2, 4])
    >>> q = 50
    >>> x = np.linspace(-180, 180, 300)[:, np.newaxis]
    >>> y, _ = special.mathieu_cem(m, q, x)
    >>> plt.plot(x, y)
    >>> plt.xlabel('x (degrees)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 2', 'm = 4'))

    Because the orders ``2`` and
    ``4`` are even, the period of each function is 180 degrees.

    """
)


mathieu_cem = _with_cache_optimization(
    name="mathieu_cem",
    arg_names=["m", "q", "x"],
    docstring=_mathieu_cem_doc,
    ufunc=_mathieu_cem,
    cache_arg_indices=[0, 1],
    module="scipy.special._mathieu",
)


_mathieu_sem_doc = (
    r"""mathieu_sem(m, q, x, out=None)

    Odd Mathieu function and its derivative.

    .. deprecated:: 2.0.0
        `mathieu_sem` is deprecated and will be removed in SciPy 2.2.0.
        Use `mathieu_se` instead.

    Returns the odd Mathieu function, se_m(x, q), of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` expressed in radians.

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*.
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Derivative with respect to `x` expressed in radians.

    See Also
    --------
    mathieu_se, mathieu_a, mathieu_b, mathieu_cem

    Notes
    -----
    Odd Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (b_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`b_m` (calculated with `mathieu_b`)
    results in an odd, periodic solution :math:`y(x)` with period 180 degrees
    (for even :math:`m`) or 360 degrees (for odd :math:`m`).

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function
    .. [2] Stuart Brorson, A New Implementation of the Mathieu Functions for SciPy.
           https://github.com/brorson/ScipyMathieuPaper

    Examples
    --------
    Plot odd Mathieu functions of orders ``2`` and ``4``.

    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> m = np.asarray([2, 4])
    >>> q = 50
    >>> x = np.linspace(-180, 180, 300)[:, np.newaxis]
    >>> y, _ = special.mathieu_sem(m, q, x)
    >>> plt.plot(x, y)
    >>> plt.xlabel('x (degrees)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 2', 'm = 4'))

    Because the orders ``2`` and
    ``4`` are even, the period of each function is 180 degrees.

    """
)

mathieu_sem = _with_cache_optimization(
    name="mathieu_sem",
    arg_names=["m", "q", "x"],
    docstring=_mathieu_sem_doc,
    ufunc=_mathieu_sem,
    cache_arg_indices=[0, 1],
    module="scipy.special._mathieu",
)


_mathieu_ce_doc = (
    r"""mathieu_ce(m, q, x, out=None)

    Even Mathieu function and its derivative.

    Returns the even Mathieu function, :math:`\mathrm{ce}_m(x, q)`, of order
    :math:`m` and parameter :math:`q` evaluated at :math:`x` (given in radians).
    Also returns its derivative with respect to :math:`x`.

    Parameters
    ----------
    m : array_like
        Order of the function. Must be a non-negative integer.
    q : array_like
        Parameter of the function.
    x : array_like
        Argument of the function, *given in radians*.
    out : tuple of ndarray, optional
        Optional output arrays for the function results.

    Returns
    -------
    y : scalar or ndarray
        Value of the function.
    yp : scalar or ndarray
        Derivative with respect to `x` (per radian).

    See Also
    --------
    mathieu_a, mathieu_b, mathieu_se

    Notes
    -----
    The even Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (a_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`a_m` (calculated with `mathieu_a`)
    results in an even, periodic solution :math:`y(x)` with period :math:`\pi`
    (for even :math:`m`) or :math:`2\pi` (for odd :math:`m`).

    .. versionadded:: 2.0.0

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function
    .. [2] Stuart Brorson, A New Implementation of the Mathieu Functions for SciPy.
           https://github.com/brorson/ScipyMathieuPaper

    Examples
    --------
    Plot even Mathieu functions of orders :math:`m = 2` and :math:`m = 4`.

    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> m = np.asarray([2, 4])
    >>> q = 50
    >>> x = np.linspace(-np.pi, np.pi, 300)[:, np.newaxis]
    >>> y, _ = special.mathieu_ce(m, q, x)
    >>> plt.plot(x, y)
    >>> plt.xticks(
    ...     [-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
    ...     [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    ... )
    >>> plt.xlabel('x (radians)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 2', 'm = 4'))
    >>> plt.show()

    Because the orders :math:`2` and :math:`4` are even, the period of each
    function is :math:`\pi`.

    Now plot the functions of odd orders :math:`m = 1` and :math:`m = 3`.

    >>> m = np.asarray([1, 3])
    >>> y, _ = special.mathieu_ce(m, q, x)
    >>> plt.figure()
    >>> plt.plot(x, y)
    >>> plt.xticks(
    ...     [-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
    ...     [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    ... )
    >>> plt.xlabel('x (radians)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 1', 'm = 3'))
    >>> plt.show()

    Because the orders :math:`1` and :math:`3` are odd, the period of each
    function is :math:`2\pi`.

    """
)

mathieu_ce = _with_cache_optimization(
    name="mathieu_ce",
    arg_names=["m", "q", "x"],
    docstring=_mathieu_ce_doc,
    ufunc=_mathieu_ce,
    cache_arg_indices=[0, 1],
    module="scipy.special._mathieu",
)


_mathieu_se_doc = (
    r"""mathieu_se(m, q, x, out=None)

    Odd Mathieu function and its derivative.

    Returns the odd Mathieu function, :math:`\mathrm{se}_m(x, q)`, of order
    :math:`m` and parameter :math:`q` evaluated at :math:`x` (given in radians).
    Also returns its derivative with respect to :math:`x`.

    Parameters
    ----------
    m : array_like
        Order of the function. Must be a non-negative integer.
    q : array_like
        Parameter of the function.
    x : array_like
        Argument of the function, *given in radians*.
    out : tuple of ndarray, optional
        Optional output arrays for the function results.

    Returns
    -------
    y : scalar or ndarray
        Value of the function.
    yp : scalar or ndarray
        Derivative with respect to `x` (per radian).

    See Also
    --------
    mathieu_a, mathieu_b, mathieu_ce

    Notes
    -----
    Odd Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (b_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`b_m` (calculated with `mathieu_b`)
    results in an odd, periodic solution :math:`y(x)` with period :math:`\pi`
    (for even :math:`m`) or :math:`2\pi` (for odd :math:`m`).

    For :math:`m = 0`, both outputs are zero.

    .. versionadded:: 2.0.0

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function
    .. [2] Stuart Brorson, A New Implementation of the Mathieu Functions for SciPy.
           https://github.com/brorson/ScipyMathieuPaper

    Examples
    --------
    Plot odd Mathieu functions of orders :math:`m = 2` and :math:`m = 4`.

    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> m = np.asarray([2, 4])
    >>> q = 50
    >>> x = np.linspace(-np.pi, np.pi, 300)[:, np.newaxis]
    >>> y, _ = special.mathieu_se(m, q, x)
    >>> plt.plot(x, y)
    >>> plt.xticks(
    ...     [-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
    ...     [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    ... )
    >>> plt.xlabel('x (radians)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 2', 'm = 4'))
    >>> plt.show()

    Because the orders :math:`2` and :math:`4` are even, the period of each
    function is :math:`\pi`.

    Now plot the functions of odd orders :math:`m = 1` and :math:`m = 3`.

    >>> m = np.asarray([1, 3])
    >>> y, _ = special.mathieu_se(m, q, x)
    >>> plt.figure()
    >>> plt.plot(x, y)
    >>> plt.xticks(
    ...     [-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
    ...     [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    ... )
    >>> plt.xlabel('x (radians)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 1', 'm = 3'))
    >>> plt.show()

    Because the orders :math:`1` and :math:`3` are odd, the period of each
    function is :math:`2\pi`.

    """
)

mathieu_se = _with_cache_optimization(
    name="mathieu_se",
    arg_names=["m", "q", "x"],
    docstring=_mathieu_se_doc,
    ufunc=_mathieu_se,
    cache_arg_indices=[0, 1],
    module="scipy.special._mathieu",
)


def _deprecated_mathieu(func, replacement):
    msg = (f"`scipy.special.{func.__name__}` is deprecated as of SciPy 2.0.0 "
           f"and will be removed in SciPy 2.2.0. Use "
           f"`scipy.special.{replacement}` instead, converting x from degrees "
           "to radians.")
    return _make_ufunc_wrapper(
        _deprecated(msg, stacklevel=3)(func), func, func.__name__,
        ["m", "q", "x"], func.__doc__, module="scipy.special",
    )
