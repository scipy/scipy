from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.linalg import solve_banded
from . import _bspl
from . import fitpack

__all__ = ["BSpline"]


class BSpline(object):
    r"""Univariate spline in the B-spline basis.

    .. math::

        S(x) = \sum_{j=0}^{n-1} c_j  B_{j, k; t}(x)

    where :math:`B_{j, k; t}` are B-spline basis functions of degree `k`
    and knots `t`.

    Parameters
    ----------
    t : ndarray, shape (n+k+1,)
        knots
    c : ndarray, shape (>=n, ...)
        spline coefficients
    k : int
        B-spline order
    extrapolate : bool, optional
        whether to extrapolate beyond the base interval, ``t[k] .. t[n]``,
        or to return nans. 
        If True, extrapolates the first and last polynomial pieces of b-spline
        functions active on the base interval.
        Default is True.

    Methods
    -------
    __call__
    basis_element
    derivative
    antiderivative
    integrate

    Notes
    -----
    B-spline basis elements are defined via

    .. math::

        S(x) = \sum_{j=0}^{n-1} c_{j} B_{j, k; t}(x)

        B_{i, 0}(x) = 1, \textrm{if $t_i \le x < t_{i+1}$, otherwise $0$,}

        B_{i, k}(x) = \frac{x - t_i}{t_{i+k} - t_i} B_{i, k-1}(x)
                 + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_{i+1}} B_{i+1, k-1}(x)

    **Implementation details**

    - At least ``k+1`` coefficients are required for a spline of degree `k`,
      so that ``n >= k+1``. Additional coefficients, ``c[j]`` with
      ``j > n``, are ignored.

    - B-spline basis elements of degree `k` form a partition of unity on the
      *base interval*, ``t[k] <= x <= t[n]``.

    Examples
    --------

    Translating the recursive definition of B-splines into Python code, we have:

    >>> def B(x, k, i, t):
    ...    if k == 0:
    ...       return 1.0 if t[i] <= x < t[i+1] else 0.0
    ...    if t[i+k] == t[i]:
    ...       c1 = 0.0
    ...    else:
    ...       c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
    ...    if t[i+k+1] == t[i+1]:
    ...       c2 = 0.0
    ...    else:
    ...       c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
    ...    return c1 + c2

    >>> def bspline(x, t, c, k):
    ...    n = len(t) - k - 1
    ...    assert (n >= k+1) and (len(c) >= n)
    ...    return sum(c[i] * B(x, k, i, t) for i in range(n))

    Note that this is an inefficient (if straightforward) way to
    evaluate B-splines --- this spline class does it in an equivalent,
    but much more efficient way.

    Here we construct a quadratic spline function on the base interval
    ``2 <= x <= 4`` and compare with the naive way of evaluating the spline:

    >>> from scipy.interpolate import BSpline
    >>> k = 2
    >>> t = [0, 1, 2, 3, 4, 5, 6]
    >>> c = [-1, 2, 0, -1]
    >>> spl = BSpline(t, c, k)
    >>> spl(2.5)
    array(1.375)
    >>> bspline(2.5, t, c, k)
    1.375

    Note that outside of the base interval results differ. This is because
    `BSpline` extrapolates the first and last polynomial pieces of b-spline
    functions active on the base interval.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> xx = np.linspace(1.5, 4.5, 50)
    >>> ax.plot(xx, [bspline(x, t, c ,k) for x in xx], 'r-', lw=3)
    >>> ax.plot(xx, spl(xx), 'b-', lw=4, alpha=0.7)
    >>> ax.grid(True)
    >>> plt.show()


    References
    ----------
    .. [1] Tom Lyche and Knut Morken, Spline methods,
        http://www.uio.no/studier/emner/matnat/ifi/INF-MAT5340/v05/undervisningsmateriale/
    .. [2] Carl de Boor, A practical guide to splines, Springer, 2001.

    """
    def __init__(self, t, c, k, extrapolate=True):
        super(BSpline, self).__init__()

        self.k = int(k)
        self.c = np.asarray(c)
        self.t = np.ascontiguousarray(t, dtype=np.float64)
        self.extrapolate = bool(extrapolate)

        n = self.t.shape[0] - self.k - 1

        if k < 0:
            raise ValueError("Spline order cannot be negative.")
        if int(k) != k:
            raise ValueError("Spline order must be integer.")
        if self.t.ndim != 1:
            raise ValueError("Knot vector must be one-dimensional.")
        if n < self.k + 1:
            raise ValueError("Need at least %d knots for degree %d" %
                    (2*k + 2, k))
        if (np.diff(self.t) < 0).any():
            raise ValueError("Knots must be in a non-decreasing order.")
        if len(np.unique(self.t[k:n+1])) < 2:
            raise ValueError("Need at least two internal knots.")
        if not np.isfinite(self.t).all():
            raise ValueError("Knots should not have nans or infs.")
        if self.c.ndim < 1:
            raise ValueError("Coefficients must be at least 1-dimensional.")
        if self.c.shape[0] < n:
            raise ValueError("Knots, coefficients and degree are inconsistent.")

        dt = self._get_dtype(self.c.dtype)
        self.c = np.ascontiguousarray(self.c, dtype=dt)

    def _get_dtype(self, dtype):
        if np.issubdtype(dtype, np.complexfloating):
            return np.complex_
        else:
            return np.float_

    @classmethod
    def _construct_fast(cls, t, c, k, extrapolate=True):
        """Construct a spline without making checks.

        Accepts same parameters as the regular constructor. Input arrays
        `t` and `c` must of correct shape and dtype.
        """
        self = object.__new__(cls)
        self.t, self.c, self.k = t, c, k
        self.extrapolate = extrapolate
        return self

    @classmethod
    def basis_element(cls, t, extrapolate=True):
        """Return a B-spline basis element ``B(x | t[0], ..., t[k+1])``.

        Parameters
        ----------
        t : ndarray, shape (k+1,)
            internal knots
        extrapolate : bool, optional
            whether to extrapolate beyond the base interval, ``t[0] .. t[k+1]``,
            or to return nans. Default is True.

        Returns
        -------
        basis_element : callable
            A callable representing a B-spline basis element for the knot
            vector `t`.

        Notes
        -----
        The order of the b-spline, `k`, is inferred from the length of `t` as
        ``len(t)-2``. The knot vector is constructed by appending and prepending
        ``k+1`` elements to internal knots `t`.

        Examples
        --------

        Construct a cubic b-spline:

        >>> from scipy.interpolate import BSpline
        >>> b = BSpline.basis_element([0, 1, 2, 3, 4])
        >>> k = b.k
        >>> b.t[k:-k]
        array([ 0.,  1.,  2.,  3.,  4.])
        >>> k
        3

        Construct a second order b-spline on ``[0, 1, 1, 2]``, and compare
        to its explicit form:

        >>> t = [-1, 0, 1, 1, 2]
        >>> b = BSpline.basis_element(t[1:])
        >>> def f(x):
        ...     return np.where(x < 1, x*x, (2. - x)**2)
        >>>
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots()
        >>> x = np.linspace(0, 2, 51)
        >>> ax.plot(x, b(x), 'g', lw=3)
        >>> ax.plot(x, f(x), 'r', lw=8, alpha=0.4)
        >>> ax.grid(True)
        >>> plt.show()

        """
        k = len(t) - 2
        t = np.r_[(t[0]-1,) * k, t, (t[-1]+1,) * k]
        c = np.zeros_like(t)
        c[k] = 1.
        return cls(t, c, k, extrapolate)

    def __call__(self, x, nu=0, extrapolate=None):
        """
        Evaluate a spline function.

        Parameters
        ----------
        x : array_like
            points to evaluate the spline at.
        nu: int, optional
            derivative to evaluate (default = 0).
        extrapolate : bool, optional
            whether to extrapolate based on the first and last intervals
            or return nans. Default is `self.extrapolate`.

        Returns
        -------
        y : array_like
            Shape is determined by replacing the interpolation axis
            in the coefficient array with the shape of `x`.

        """
        if extrapolate is None:
            extrapolate = self.extrapolate
        x = np.asarray(x)
        x_shape = x.shape
        x = np.ascontiguousarray(x.ravel(), dtype=np.float_)
        out = np.empty((len(x), int(np.prod(self.c.shape[1:]))),
                dtype=self.c.dtype)
        self._ensure_c_contiguous()
        self._evaluate(x, nu, extrapolate, out)
        return out.reshape(x_shape + self.c.shape[1:])

    def _evaluate(self, xp, nu, extrapolate, out):
        _bspl.evaluate_spline(self.t, self.c.reshape(self.c.shape[0], -1),
                self.k, xp, nu, extrapolate, out)

    def _ensure_c_contiguous(self):
        """
        c and t may be modified by the user. The Cython code expects
        that they are C contiguous.

        """
        if not self.t.flags.c_contiguous:
            self.x = self.x.copy()
        if not self.c.flags.c_contiguous:
            self.c = self.c.copy()

    def derivative(self, nu=1):
        """Return a b-spline representing the derivative.

        Parameters
        ----------
        nu : int, optional
            Derivative order.
            Default is 1.

        Returns
        -------
        b : BSpline object
            A new instance representing the derivative.

        See Also
        --------
        splder, splatinder

        """
        c = self.c
        # pad the c array if needed
        ct = len(self.t) - len(c)
        if ct > 0:
            c = np.r_[c, np.zeros((ct,) + c.shape[1:])]
        tck = fitpack.splder((self.t, c, self.k), nu)
        return self._construct_fast(*tck, extrapolate=self.extrapolate)

    def antiderivative(self, nu=1):
        """Return a b-spline representing the antiderivative.

        Parameters
        ----------
        nu : int, optional
            Antiderivative order. Default is 1.

        Returns
        -------
        b : BSpline object
            A new instance representing the antiderivative.

        See Also
        --------
        splder, splantider

        """
        c = self.c
        # pad the c array if needed
        ct = len(self.t) - len(c)
        if ct > 0:
            c = np.r_[c, [0]*ct]
        tck = fitpack.splantider((self.t, c, self.k), nu)
        return self._construct_fast(*tck, extrapolate=self.extrapolate)

    def integrate(self, a, b, extrapolate=None):
        """Compute a definite integral of the spline.

        Parameters
        ----------
        a : float
            Lower limit of integration.
        b : float
            Upper limit of integration.
        extrapolate : bool, optional
            whether to extrapolate beyond the base interval, ``t[k] .. t[-k-1]``,
            or take the spline to be zero outside of the base interval.
            Default is True.

        Returns
        -------
        I : array_like
            Definite integral of the spline over the interval ``[a, b]``.

        Examples
        --------
        Construct the linear spline ``x if x < 1 else 2 - x`` on the base
        interval :math:`[0, 2]`, and integrate it

        >>> from scipy.interpolate import BSpline
        >>> b = BSpline.basis_element([0, 1, 2])
        >>> b.integrate(0, 1)
        array(0.5)


        If the integration limits are outside of the base interval, the result
        is controlled by the `extrapolate` parameter

        >>> b.integrate(-1, 1)
        array(0.0)
        >>> b.integrate(-1, 1, extrapolate=False)
        array(1.0)

        """
        if extrapolate is None:
            extrapolate = self.extrapolate

        if not extrapolate:
            # shrink the integration interval, if needed
            a = max(a, self.t[self.k])
            b = min(b, self.t[-self.k - 1])

        # compute the antiderivative
        c = self.c
        ct = len(self.t) - len(c)
        if ct > 0:
            c = np.r_[c, [0]*ct]
        t, c, k = fitpack.splantider((self.t, c, self.k), 1)

        # prepare x & c
        if not c.flags.c_contiguous:
            c = c.copy()

        # evaluate the diff of antiderivatives
        x = np.asarray([a, b], dtype=np.float_)
        out = np.empty((2, int(np.prod(c.shape[1:]))),
                dtype=c.dtype)
        _bspl.evaluate_spline(t, c.reshape(c.shape[0], -1),
                k, x, 0, extrapolate, out)
        out = out[1] - out[0]
        return out.reshape(c.shape[1:])
