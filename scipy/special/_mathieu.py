from scipy.special._ufuncs import _mathieu_cem, _mathieu_sem
from scipy.special._ufunc_tools import _with_cache_optimization

__all__ = ["mathieu_cem", "mathieu_sem"]

_mathieu_cem_doc = (
    r"""mathieu_cem(m, q, x, out=None)

    Even Mathieu function and its derivative.

    Returns the even Mathieu function, ``ce_m(x, q)``, of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` of ce_m(x, q)

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
        Value of the derivative vs x

    See Also
    --------
    mathieu_a, mathieu_b, mathieu_sem

    Notes
    -----
    The even Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (a_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`a_m` (calculated with `mathieu_a`)
    results in an odd, periodic solution :math:`y(x)` with period 180 degrees
    (for even :math:`m`) or 360 degrees (for odd :math:`m`).

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function
    .. [2] https://github.com/brorson/ScipyMathieuPaper

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
)


_mathieu_sem_doc = (
    r"""mathieu_sem(m, q, x, out=None)

    Odd Mathieu function and its derivative.

    Returns the odd Mathieu function, se_m(x, q), of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` of se_m(x, q).

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
        Value of the derivative vs x

    See Also
    --------
    mathieu_a, mathieu_b, mathieu_cem

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
    .. [2] https://github.com/brorson/ScipyMathieuPaper

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
)
