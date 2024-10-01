import math
import heapq
import itertools

from dataclasses import dataclass, field
from types import ModuleType
from typing import Any, TYPE_CHECKING

from scipy._lib._array_api import array_namespace, xp_size, xp_copy
from scipy._lib._util import MapWrapper

from scipy.integrate._rules import (
    ProductNestedFixed,
    GaussKronrodQuadrature,
    GenzMalikCubature,
)
from scipy.integrate._rules._base import _split_patch

__all__ = ['cubature']

if TYPE_CHECKING:
    Array = Any  # To be changed to a Protocol later (see array-api#589)
else:
    Array = object


@dataclass
class CubaturePatch:
    estimate: Array
    error: Array
    a: Array
    b: Array
    _xp: ModuleType = field(repr=False)

    def __lt__(self, other):
        # Consider patches with higher error estimates as being "less than" patches with
        # lower order estimates, so that patches with high error estimates are placed at
        # the top of the heap.

        this_err = self._xp.max(self._xp.abs(self.error))
        other_err = self._xp.max(self._xp.abs(other.error))

        return this_err > other_err


@dataclass
class CubatureResult:
    estimate: Array
    error: Array
    status: str
    patches: list[CubaturePatch]
    subdivisions: int
    atol: float
    rtol: float


def cubature(f, a, b, rule="gk21", rtol=1e-8, atol=0, max_subdivisions=10000,
             args=(), workers=1, points=None, region=None):
    r"""
    Adaptive cubature of multidimensional array-valued function.

    Given an arbitrary integration rule, this function returns an estimate of the
    integral to the requested tolerance over the region defined by the arrays `a` and
    `b` specifying the corners of a hypercube, or over a more complex region specified
    with `region`.

    Convergence is not guaranteed for all integrals.

    Parameters
    ----------
    f : callable
        Function to integrate. `f` must have the signature::

            f(x : ndarray, *args) -> ndarray

        `f` should accept arrays ``x`` of shape::

            (npoints, ndim)

        and output arrays of shape::

            (npoints, output_dim_1, ..., output_dim_n)

        In this case, `cubature` will return arrays of shape::

            (output_dim_1, ..., output_dim_n)
    a, b : array_like
        Lower and upper limits of integration as 1D arrays specifying the left and right
        endpoints of the intervals being integrated over. Can specify infinite limits.
    rule : str, optional
        Rule used to estimate the integral. If passing a string, the options are
        "gauss-kronrod" (21 node), or "genz-malik" (degree 7). If a rule like
        "gauss-kronrod" is specified for an ``n``-dim integrand, the corresponding
        Cartesian product rule is used. "gk21", "gk15" are also supported for
        compatibility with `quad_vec`. See Notes.
    rtol, atol : float, optional
        Relative and absolute tolerances. Iterations are performed until the error is
        estimated to be less than ``atol + rtol * abs(est)``. Here `rtol` controls
        relative accuracy (number of correct digits), while `atol` controls absolute
        accuracy (number of correct decimal places). To achieve the desired `rtol`, set
        `atol` to be smaller than the smallest value that can be expected from
        ``rtol * abs(y)`` so that `rtol` dominates the allowable error. If `atol` is
        larger than ``rtol * abs(y)`` the number of correct digits is not guaranteed.
        Conversely, to achieve the desired `atol`, set `rtol` such that
        ``rtol * abs(y)`` is always smaller than `atol`. Default values are 1e-8 for
        `rtol` and 0 for `atol`.
    max_subdivisions : int, optional
        Upper bound on the number of subdivisions to perform. Default is 10,000.
    args : tuple, optional
        Additional positional args passed to `f`, if any.
    workers : int or map-like callable, optional
        If `workers` is an integer, part of the computation is done in parallel
        subdivided to this many tasks (using :class:`python:multiprocessing.pool.Pool`).
        Supply `-1` to use all cores available to the Process. Alternatively, supply a
        map-like callable, such as :meth:`python:multiprocessing.pool.Pool.map` for
        evaluating the population in parallel. This evaluation is carried out as
        ``workers(func, iterable)``.
    points : list of array_like, optional
        List of points by which to split up the initial region of integration. If the
        rule being used does not evaluate `f` on the boundary of a region (which is the
        case for all Genz-Malik and Gauss-Kronrod rules) then `f` will not be evaluated
        at these points. This should be a list of array-likes where each element has
        length ``ndim``. Default is empty. See Examples.
    region : list of callable, optional
        List of callables parameterising a potentially non-rectangular region of
        integration.

        To perform an integral like:

        .. math::

            \int^{b_1}_{a_1}
            \cdots
            \int^{b_n}_{a_n}
            \int^{B_1(x_1, \ldots, x_n)}_{A_1(x_1, \ldots, x_n)}
            \cdots
            \int^{B_m(x_1, \ldots, x_{n+m})}_{A_m(x_1, \ldots, x_{n+m})}
            f(x_1, \ldots, x_{n+m})
            dx_{n+m} \cdots dx_1

        where `f` has ndim ``n+m``, ``n`` variables of which range over the intervals
        ``[a_1, b_1], ..., [a_n, b_n]`` and ``m`` variables of which are described by
        the functions ``B_i(x_1, ..., x_{n+i})``, this can be achieved by setting `a` to
        ``[a_1, ..., a_n]``, `b` to ``[b_1, ..., b_n]``, and passing `region` as::

            region = [
                region_func_1,
                ...,
                region_func_m,
            ]

        where each ``region_func_i`` is a function which calculates::

            region_func_i = lambda x: (
                # Upper limits given x_1, ..., x_{i-1}
                [ A_i(x[:, 0], ..., x[:, i-1]) ],

                # Lower limits given x_1, ..., x_{i-1}
                [ B_i(x[:, 0], ..., x[:, i-1]) ],
            )

        Each ``region_func_i`` should accept arrays of shape
        ``(npoints, num_preceding variables)`` and
        return a tuple of two arrays with shape ``(npoints, num_new_variables)``. In the
        above case ``num_new_variables`` is 1, but it is also possible to return the
        limits for multiple variables at once. See Examples.

        `a` and `b` may contain infinite limits, but none of the functions in ``region``
        can return infinity.

        More details can be found in the Notes section.

    Returns
    -------
    res : object
        Object containing the results of the estimation. It has the following
        attributes:

        estimate : ndarray
            Estimate of the value of the integral over the overall region specified.
        error : ndarray
            Estimate of the error of the approximation over the overall region
            specified.
        status : str
            Whether the estimation was successful. Can be either: "converged",
            "not_converged".
        subdivisions : int
            Number of subdivisions performed.
        atol, rtol : float
            Requested tolerances for the approximation.
        patches: list of object
            List of objects containing the estimates of the integral over smaller
            patches of the domain.


        Each object in ``patches`` has the following attributes:

        a, b : ndarray
            Points describing the corners of the patch. If the original integral
            contained infinite limits or was over a region described by `region`,
            then `a` and `b` are in the transformed coordinates.
        estimate : ndarray
            Estimate of the value of the integral over this patch.
        error : ndarray
            Estimate of the error of the approximation over this patch.


    Notes
    -----
    The algorithm uses a similar algorithm to `quad_vec`, which itself is based on the
    implementation of QUADPACK's DQAG* algorithms, implementing global error control and
    adaptive subdivision.

    The source of the nodes and weights used for Gauss-Kronrod quadrature can be found
    in [1]_, and the algorithm for calculating the nodes and weights in Genz-Malik
    cubature can be found in [2]_.

    The rules currently supported via the `rule` argument are:

    - "gauss-kronrod", 21-node Gauss-Kronrod
    - "genz-malik", ``n``-node Genz-Malik

    If using Gauss-Kronrod for an ``n``-dim integrand where ``n > 2``, then the
    corresponding Cartesian product rule will be found by taking the Cartesian product
    of the nodes in the 1D case. This means that the number of nodes scales
    exponentially as ``21^n`` in the Gauss-Kronrod case, which may be problematic in a
    moderate number of dimensions. Although Genz-Malik is typically less accurate than
    Gauss-Kronrod it has much fewer nodes, so in this situation using "genz-malik" might
    be preferable.

    Infinite limits are handeled with an appropriate variable transformation. Assuming
    ``a = [a_1, ..., a_n]`` and ``b = [b_1, ..., b_n]``:

    If :math:`a_i` and :math:`b_i` range over :math:`x \in (-\infty, \infty)`, the i-th
    integration variable will use the transformation :math:`x = \frac{1-|t|}{t}` and
    :math:`t \in (-1, 1)`.

    If :math:`a_i` and :math:`b_i` range over :math:`x \in [a_i, \infty)`, the i-th
    integration variable will use the transformation :math:`x = a_i + \frac{1-t}{t}` and
    :math:`t \in (0, 1)`.

    If :math:`a_i` and :math:`b_i` range over :math:`x \in (-\infty, b_i]`, the i-th
    integration variable will use the transformation :math:`x = b_i - \frac{1+t}{t}` and
    :math:`t \in (0, 1)`.

    In all three of these cases, the Jacobian of the transformation is
    :math:`J(t) = t^{-2}`.

    When `region` is specified, the variable transformation described in [3]_ is used.
    If the original integral has the form

    .. math::

        \int^{b_1}_{a_1}
        \cdots
        \int^{b_n}_{a_n}
        \int^{B_{n+1}(x_1, \ldots, x_n)}_{A_{n+1}(x_1, \ldots, x_n)}
        \cdots
        \int^{B_{n+m}(x_1, \ldots, x_{n+m-1})}_{A_{n+m}(x_1, \ldots, x_{n+m-1})}
            f(x_1, \ldots, x_{n+m})
        dx_{n+m} \cdots dx_1

    :math:`x_1, \ldots, x_n` are referred to as the "outer variables" and their limits
    are specified by the constant intervals :math:`(a_i, b_i)` for
    :math:`1 \le i \le n`.

    :math:`x_{n+1}, \ldots, x_{n+m}` are referred to as the "inner variables" and their
    limits given the values of the preceding variables range over the interval
    :math:`(A_i(x_1, \ldots, x_{i-1}), B_i(x_1, \ldots, x_{i-1})`.

    Integrals of this form will be transformed into an integral over constant limits:

    .. math::

        \int^{b_1}_{a_1}
        \cdots
        \int^{b_n}_{a_n}
        \int^{1}_{-1}
        \cdots
        \int^{1}_{-1}
            g(t_1, \ldots, t_{n+m})
        dt_{n+m} \cdots dt_1

    The outer variables are unchanged, :math:`x_i = t_i` for :math:`1 \le x_i \le t_i`.

    The inner variables are mapped via the transformation:

    .. math::

        x_i = \frac{
            B_{i}(x_1, \ldots, x_{i-1}) + A_{i}(x_1, \ldots, x_{i-1})
        }{2} + t_i \frac{
            B_{i}(x_1, \ldots, x_{i-1}) - A_{i}(x_1, \ldots, x_{i-1})
        }{2}

    for :math:`i > n`. This transformation has Jacobian determinant:

    .. math::

        J(x_1, \ldots, x_{n+m})
        =
        \prod^{n+m}_{i = n+1}
        \frac{B_{i}(x_1, \ldots, x_{i-1}) - A_{i}(x_1, \ldots, x_{i-1})}{2}

    To specify this type of integral programatically, the limits should be specified as
    two arrays, ``a`` and ``b`` for the constant limits, and a list of callables
    ``region`` describing the function limits::

        a = [a_1, ..., a_n]
        b = [b_1, ..., b_n]

        region = [
            region_func_1,
            ...,
            region_func_k,
        ]

    Each ``region_func_i`` should accepts arrays of shape
    ``(npoints, num_preceding_variables)`` and return a tuple of two arrays, one for
    the next set of lower limits given the values of the preceding variables, and one
    for the next set of upper limits given the values of the preceding variables. Both
    of the returned arrays should have shape ``(npoints, num_new_variables)``.

    References
    ----------
    .. [1] R. Piessens, E. de Doncker, Quadpack: A Subroutine Package for Automatic
        Integration, files: dqk21.f, dqk15.f (1983).

    .. [2] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.

    .. [3] Philip J. Davis, Philip Rabinowitz, Methods of Numerical Integration,
        Second Edition, Academic Press, Section 5.4 (1984).

    Examples
    --------
    **1D integral with vector output**:

    .. math::

        \int^1_0 \mathbf f(x) \text dx

    Where ``f(x) = x^n`` and ``n = np.arange(10)`` is a vector. Since no rule is
    specified, the default "gk21" is used, which corresponds to Gauss-Kronrod
    integration with 21 nodes.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> def f(x, n):
    ...    # Make sure x and n are broadcastable
    ...    return x[:, np.newaxis]**n[np.newaxis, :]
    >>> res = cubature(
    ...     f,
    ...     a=[0],
    ...     b=[1],
    ...     args=(np.arange(10),),
    ... )
    >>> res.estimate
     array([1.        , 0.5       , 0.33333333, 0.25      , 0.2       ,
            0.16666667, 0.14285714, 0.125     , 0.11111111, 0.1       ])

    **7D integral with arbitrary-shaped array output**::

        f(x) = cos(2*pi*r + alphas @ x)

    for some ``r`` and ``alphas``, and the integral is performed over the unit
    hybercube, :math:`[0, 1]^7`. Since the integral is in a moderate number of
    dimensions, "genz-malik" is used rather than the default "gauss-kronrod" to
    avoid constructing a product rule with :math:`21^7 \approx 2 \times 10^9` nodes.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> def f(x, r, alphas):
    ...     # f(x) = cos(2*pi*r + alphas @ x)
    ...     # Need to allow r and alphas to be arbitrary shape
    ...     npoints, ndim = x.shape[0], x.shape[-1]
    ...     alphas = alphas[np.newaxis, ...]
    ...     x = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)
    ...     return np.cos(2*np.pi*r + np.sum(alphas * x, axis=-1))
    >>> rng = np.random.default_rng()
    >>> r, alphas = rng.random((2, 3)), rng.random((2, 3, 7))
    >>> res = cubature(
    ...     f=f,
    ...     a=np.array([0, 0, 0, 0, 0, 0, 0]),
    ...     b=np.array([1, 1, 1, 1, 1, 1, 1]),
    ...     rtol=1e-5,
    ...     rule="genz-malik",
    ...     args=(r, alphas),
    ... )
    >>> res.estimate
     array([[-0.79812452,  0.35246913, -0.52273628],
            [ 0.88392779,  0.59139899,  0.41895111]])

    **Parallel computation with** `workers`:

    >>> from concurrent.futures import ThreadPoolExecutor
    >>> with ThreadPoolExecutor() as executor:
    ...     res = cubature(
    ...         f=f,
    ...         a=np.array([0, 0, 0, 0, 0, 0, 0]),
    ...         b=np.array([1, 1, 1, 1, 1, 1, 1]),
    ...         rtol=1e-5,
    ...         rule="genz-malik",
    ...         args=(r, alphas),
    ...         workers=executor.map,
    ...      )
    >>> res.estimate
     array([[-0.79812452,  0.35246913, -0.52273628],
            [ 0.88392779,  0.59139899,  0.41895111]])

    When this is done with process-based parallelization (as would be the case passing
    `workers` as an integer) you should ensure the main module is import-safe.

    **2D integral with infinite limits**:

    .. math::

        \int^{ \infty }_{ -\infty }
        \int^{ \infty }_{ -\infty }
            e^{-x^2-y^2}
        \text dy
        \text dx

    >>> def gaussian(x):
    ...     return np.exp(-np.sum(x**2, axis=-1))
    >>> res = cubature(gaussian, [-np.inf, -np.inf], [np.inf, np.inf])
    >>> res.estimate
     3.1415926

    **1D integral with singularities avoided using** `points`:

    .. math::

        \int^{ 1 }_{ -1 }
          \frac{\sin(x)}{x}
        \text dx

    It is necessary to use the `points` parameter to avoid evaluating `f` at the origin.

    >>> def sinc(x):
    ...     return np.sin(x)/x
    >>> res = cubature(sinc, [-1], [1], points=[[0]])
    >>> res.estimate
     1.8921661

    **2D integral over a circle using** `region`:

    Calculate the mass of a series of discs with different densities:

    .. math::

        \int^{ 1 }_{ -1 }
        \int^{ \sqrt{1-x^2} }_{ -\sqrt{1-x^2} }
          \mathbf \rho(x_1, x_2)
        \text dx_2 \text dx_1

    >>> def density(x_arr):
    ...     x_0, x_1 = x_arr[:, 0], x_arr[:, 1]
    ...     # For example, 4 different densities:
    ...     return (x_0**2 + x_1**2)[:, np.newaxis] / np.arange(1, 5)[np.newaxis, :]
    >>> # density takes arrays of shape (npoints, 2) and returns arrays of shape
    >>> # (npoints, 4)
    >>> density(np.array([[-0.5, 0.5], [1, 1]]))
     array([[0.5       , 0.25      , 0.16666667, 0.125     ],
            [2.        , 1.        , 0.66666667, 0.5       ]])
    >>> res = cubature(
    ...     density,
    ...     # Outer limits:
    ...     a=[-1],
    ...     b=[1],
    ...     # Inner limits. In the problem there is only one pair of function limits, so
    ...     # region is a list of length 1.
    ...     region=[
    ...           # This function takes the value of the outer variables, which in this
    ...           # case is the value of x_0. It then returns the new limits for the
    ...           # inner variables, which in this case is just x_1.
    ...           #
    ...           # In general, there might be more outer variables and more inner
    ...           # variables, and new limits might need to be calculated for several
    ...           # different values of outer variables at once. This means x_arr is of
    ...           # shape (npoints, num_outer) and it needs to return a tuple of arrays
    ...           # of shape (npoints, num_inner).
    ...           lambda x_arr: (
    ...               -np.sqrt(1-x_arr[:, 0]**2)[:, np.newaxis],
    ...                np.sqrt(1-x_arr[:, 0]**2)[:, np.newaxis],
    ...           ),
    ...     ],
    ... )
    >>> res.estimate
     array([1.57079633, 0.78539816, 0.52359878, 0.39269908])

    **3D integral over a sphere with** `region`:

    .. math::

        \int^{ 1 }_{ -1 }
        \int^{ \sqrt{1-x_0^2} }_{ -\sqrt{1-x_0^2} }
        \int^{ \sqrt{1-x_0^2-x_1^2} }_{ -\sqrt{1-x_0^2-x_1^2} }
            1
        \text dx_2
        \text dx_1
        \text dx_0

    >>> res = cubature(
    ...     lambda x_arr: np.ones(x_arr.shape[0]),
    ...     a=[-1],
    ...     b=[1],
    ...     rtol=1e-5,  # Reduce tolerance, can be slow otherwise
    ...     region=[
    ...         lambda x_arr: (
    ...             -np.sqrt(1-x_arr[:, 0]**2)[:, np.newaxis],
    ...              np.sqrt(1-x_arr[:, 0]**2)[:, np.newaxis],
    ...         ),
    ...         lambda x_arr: (
    ...             -np.sqrt(1-x_arr[:, 0]**2-x_arr[:, 1]**2)[:, np.newaxis],
    ...              np.sqrt(1-x_arr[:, 0]**2-x_arr[:, 1]**2)[:, np.newaxis],
    ...         ),
    ...     ],
    ... )
    >>> res.estimate  # True value 4*pi/3
     4.188792

    **Returning multiple limits at once using** `region`:

    In the above uses of `region`, each function in the list has returned an array of
    shape ``(npoints, 1)``. In the previous example, it is **not** possible to combine
    these two functions into one:

    >>> # Won't work:
    >>> res = cubature(
    ...     lambda x_arr: np.ones(x_arr.shape[0]),
    ...     a=[-1],
    ...     b=[1],
    ...     rtol=1e-5,
    ...     region=[
    ...         lambda x_arr: (
    ...             [
    ...                 -np.sqrt(1-x_arr[:, 0]**2),
    ...                 -np.sqrt(1-x_arr[:, 0]**2-x_arr[:, 1]**2),
    ...                 #                         ^^^^^^^^^^^
    ...                 # This will cause an error since x_arr[:, 1] depends on the
    ...                 # value of the previous limit, which has not yet been returned.
    ...             ],
    ...             [
    ...                 np.sqrt(1-x_arr[:, 0]**2),
    ...                 np.sqrt(1-x_arr[:, 0]**2-x_arr[:, 1]**2),
    ...                 # Also a problem here:   ^^^^^^^^^^^
    ...             ],
    ...         ),
    ...     ],
    ... ) # doctest: +SKIP
     Traceback (most recent call last):
      ...
     IndexError: index 1 is out of bounds for axis 1 with size 1

    This is because the expression ``np.sqrt(1-x_arr[:, 0]**2-x_arr[:, 1]**2)`` depends
    on the value of ``x[:, 1]`` which will not be known when this function is called.
    This is because the limits for ``x[:, 1]`` cannot be determined without calling
    this function.

    If you have an integral where a group of limits do not depend on one another, you
    can return them together. For example:

    .. math::

        \int^{ 1 }_{ 0 }
        \int^{ x_0 }_{ 0 }
        \int^{ x_0 }_{ 0 } \text dx_2 \text dx_1 \text dx_0

    Here the limits for :math:`x_1` and :math:`x_2` depend only on :math:`x_0`. This
    means it is possible to specify `region` like so:

    >>> res = cubature(
    ...     lambda x_arr: np.ones(x_arr.shape[0]),
    ...     a=[-1],
    ...     b=[1],
    ...     region=[
    ...         lambda x_arr: (
    ...             np.zeros((x_arr.shape[0], 2)),            # [0, 0]
    ...             np.repeat(x_arr[:, 0], 2).reshape(-1, 2), # [x_0, x_0]
    ...         ),
    ...     ],
    ... )
    >>> res.estimate
     0.666666
    """

    # It is also possible to use a custom rule, but this is not yet part of the public
    # API. An example of this can be found in the class scipy.integrate._rules.Rule.

    xp = array_namespace(a, b)
    max_subdivisions = float("inf") if max_subdivisions is None else max_subdivisions
    points = [] if points is None else [xp.asarray(p, dtype=xp.float64) for p in points]

    # Convert a and b to arrays
    a = xp.asarray(a, dtype=xp.float64)
    b = xp.asarray(b, dtype=xp.float64)

    if xp_size(a) == 0 or xp_size(b) == 0:
        raise ValueError("`a` and `b` must be nonempty")

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("`a` and `b` must be 1D arrays")

    # If `region` is specified, apply a transformation to reduce to integration over
    # a rectangular region.
    if region is not None:
        f = _FuncLimitsTransform(f, a, b, region)
        a, b = f.transformed_limits

        # Map break points specified in original coordinates to the new coordinates
        points = [f.inv(point) for point in points]

        # Add any new problematic points introduced by the transformation
        points.extend(f.points)

    # If any of limits are the wrong way around (a > b), flip them and keep track of
    # the sign.
    sign = (-1.0) ** xp.sum(xp.astype(a > b, xp.float64))

    a_flipped = xp.min(xp.stack([a, b]), axis=0)
    b_flipped = xp.max(xp.stack([a, b]), axis=0)

    a, b = a_flipped, b_flipped

    # If any of the limits are infinite, apply a transformation
    if xp.any(xp.isinf(a)) or xp.any(xp.isinf(b)):
        f = _InfiniteLimitsTransform(f, a, b)
        a, b = f.transformed_limits

        points = [f.inv(point) for point in points]
        points.extend(f.points)

    # If any problematic points are specified, divide the initial patch so that these
    # points lie on the edge of a subpatch.
    #
    # This means ``f`` won't be evaluated there if the rule being used has no evaluation
    # points on the boundary.

    if points == []:
        initial_patch = [(a, b)]
    else:
        initial_patch = _split_patch_at_points(a, b, points)

    # If the rule is a string, convert to a corresponding product rule
    if isinstance(rule, str):
        ndim = xp_size(a)

        if rule == "genz-malik":
            rule = GenzMalikCubature(ndim, xp=xp)
        else:
            quadratues = {
                "gauss-kronrod": GaussKronrodQuadrature(21, xp=xp),

                # Also allow names quad_vec uses:
                "gk21": GaussKronrodQuadrature(21, xp=xp),
                "gk15": GaussKronrodQuadrature(15, xp=xp),
            }

            base_rule = quadratues.get(rule)

            if base_rule is None:
                raise ValueError(f"unknown rule {rule}")

            rule = ProductNestedFixed([base_rule] * ndim)

    patches = []
    est = 0.0
    err = 0.0

    for a_k, b_k in initial_patch:
        # If any of the initial patches have zero width in one dimension, we can
        # ignore this as the integral will be 0 there.
        est_k = rule.estimate(f, a_k, b_k, args)
        err_k = rule.estimate_error(f, a_k, b_k, args)
        patches.append(CubaturePatch(est_k, err_k, a_k, b_k, xp))

        est += est_k
        err += err_k

    subdivisions = 0
    success = True

    with MapWrapper(workers) as mapwrapper:
        while xp.any(err > atol + rtol * xp.abs(est)):
            # patch_k is the patch with highest estimated error.
            patch_k = heapq.heappop(patches)

            est_k = patch_k.estimate
            err_k = patch_k.error

            a_k, b_k = patch_k.a, patch_k.b

            # Subtract the estimate of the integral and its error over this patch from
            # the current global estimates, since these will be refined in the loop over
            # all subpatches.
            est -= est_k
            err -= err_k

            # Find all 2^ndim subpatches formed by splitting patch_k along each axis,
            # e.g. for 1D integrals this splits an estimate over an interval into an
            # estimate over two subintervals, for 3D integrals this splits an estimate
            # over a cube into 8 subcubes.
            #
            # For each of the new subpatches, calculate an estimate for the integral and
            # the error there, and push these patches onto the heap for potential
            # further subdividing.

            executor_args = zip(
                itertools.repeat(f),
                itertools.repeat(rule),
                itertools.repeat(args),
                _split_patch(a_k, b_k),
            )

            for subdivision_result in mapwrapper(_process_subpatch, executor_args):
                a_k_sub, b_k_sub, est_sub, err_sub = subdivision_result

                est += est_sub
                err += err_sub

                new_patch = CubaturePatch(est_sub, err_sub, a_k_sub, b_k_sub, xp)

                heapq.heappush(patches, new_patch)

            subdivisions += 1

            if subdivisions >= max_subdivisions:
                success = False
                break

        status = "converged" if success else "not_converged"

        # Apply sign change to handle any limits which were initially flipped.
        est = sign * est

        return CubatureResult(
            estimate=est,
            error=err,
            status=status,
            subdivisions=subdivisions,
            patches=patches,
            atol=atol,
            rtol=rtol,
        )


def _process_subpatch(data):
    f, rule, args, coord = data
    a_k_sub, b_k_sub = coord

    est_sub = rule.estimate(f, a_k_sub, b_k_sub, args)
    err_sub = rule.estimate_error(f, a_k_sub, b_k_sub, args)

    return a_k_sub, b_k_sub, est_sub, err_sub


def _is_strictly_in_patch(point, a, b):
    xp = array_namespace(point, a, b)

    if xp.all(point == a) or xp.all(point == b):
        return False

    return xp.all(a <= point) and xp.all(point <= b)


def _split_patch_at_points(a, b, points):
    """
    Given the integration limits `a` and `b` describing a rectangular patch and a list
    of `points`, find the list of ``[(a_1, b_1), ..., (a_l, b_l)]`` which breaks up the
    initial patch into smaller subpatches such that no `points` lie strictly inside
    any of the subpatches.
    """

    xp = array_namespace(a, b)
    patches = [(a, b)]

    for point in points:
        if xp.any(xp.isinf(point)):
            # If a point is specified at infinity, ignore.
            #
            # This case occurs when points are given by the user to avoid, but after
            # applying a transformation, they are removed.
            continue

        new_subpatches = []

        for a_k, b_k in patches:
            if _is_strictly_in_patch(point, a_k, b_k):
                subpatches = _split_patch(a_k, b_k, point)

                for left, right in subpatches:
                    # Skip any zero-width patches.
                    if xp.any(left == right):
                        continue
                    else:
                        new_subpatches.append((left, right))

                new_subpatches.extend(subpatches)

            else:
                new_subpatches.append((a_k, b_k))

        patches = new_subpatches

    return patches


class _VariableTransform:
    """
    A transformation that can be applied to an integral.
    """

    @property
    def transformed_limits(self):
        """
        New limits of integration after applying the transformation.
        """

        raise NotImplementedError

    @property
    def points(self):
        """
        Any problematic points introduced by the transformation.

        These should be specified as points where ``_VariableTransform(f)(self, point)``
        would be problematic.

        For example, if the transformation ``x = 1/((1-t)(1+t))`` is applied to a
        univariate integral, then points should return ``[ [1], [-1] ]``.
        """

        return []

    def inv(self, x):
        """
        Map points ``x`` to ``t`` such that if ``f`` is the original function and ``g``
        is the function after the transformation is applied, then::

            f(x) = g(self.inv(x))
        """

        raise NotImplementedError

    def __call__(self, t, *args, **kwargs):
        """
        Apply the transformation to ``f`` and multiply by the Jacobian determinant.
        This should be the new integrand after the transformation has been applied so
        that the following is satisfied::

            f_transformed = _VariableTransform(f)

            cubature(f, a, b) == cubature(
                f_transformed,
                *f_transformed.transformed_limits(a, b),
            )
        """

        raise NotImplementedError


class _InfiniteLimitsTransform(_VariableTransform):
    r"""
    Transformation for handling infinite limits.

    Assuming ``a = [a_1, ..., a_n]`` and ``b = [b_1, ..., b_n]``:

    If :math:`a_i` and :math:`b_i` range over :math:`x \in (-\infty, \infty)`, the i-th
    integration variable will use the transformation :math:`x = \frac{1-|t|}{t}` and
    :math:`t \in (-1, 1)`.

    If :math:`a_i` and :math:`b_i` range over :math:`x \in [a_i, \infty)`, the i-th
    integration variable will use the transformation :math:`x = a_i + \frac{1-t}{t}` and
    :math:`t \in (0, 1)`.

    If :math:`a_i` and :math:`b_i` range over :math:`x \in (-\infty, b_i]`, the i-th
    integration variable will use the transformation :math:`x = b_i - \frac{1+t}{t}` and
    :math:`t \in (0, 1)`.

    In all three of these cases, the Jacobian of the transformation is
    :math:`J(t) = t^{-2}`.
    """

    def __init__(self, f, a, b):
        self._xp = array_namespace(a, b)

        self._f = f
        self._orig_a = a
        self._orig_b = b

        self._negate_pos = []
        self._semi_inf_pos = []
        self._double_inf_pos = []

        ndim = xp_size(a)

        for i in range(ndim):
            if a[i] == -math.inf and b[i] == math.inf:
                # (-oo, oo) will be mapped to (-1, 1).
                self._double_inf_pos.append(i)

            elif a[i] != -math.inf and b[i] == math.inf:
                # (start, oo) will be mapped to (0, 1).
                self._semi_inf_pos.append(i)

            elif a[i] == -math.inf and b[i] != math.inf:
                # (-oo, end) will be mapped to (0, 1).
                #
                # This is handled by making the transformation t = -x and reducing it to
                # the other semi-infinite case.
                self._negate_pos.append(i)
                self._semi_inf_pos.append(i)

                # Since we flip the limits, we don't need to separately multiply the
                # integrand by -1.
                self._orig_a[i] = -b[i]
                self._orig_b[i] = -a[i]

    @property
    def transformed_limits(self):
        a = xp_copy(self._orig_a)
        b = xp_copy(self._orig_b)

        for index in self._double_inf_pos:
            a[index] = -1
            b[index] = 1

        for index in self._semi_inf_pos:
            a[index] = 0
            b[index] = 1

        return a, b

    @property
    def points(self):
        # If there are infinite limits, then the origin becomes a problematic point
        # due to a division by zero there.
        if len(self._double_inf_pos) != 0 or len(self._semi_inf_pos) != 0:
            return [self._xp.zeros(self._orig_a.shape)]
        else:
            return []

    def inv(self, x):
        t = xp_copy(x)

        for i in self._negate_pos:
            t[..., i] *= -1

        for i in self._double_inf_pos:
            t[..., i] = 1/(x[..., i] + self._xp.sign(x[..., i]))

        for i in self._semi_inf_pos:
            t[..., i] = 1/(x[..., i] - self._orig_a[i] + 1)

        return t

    def __call__(self, t, *args, **kwargs):
        x = xp_copy(t)
        jacobian = 1.0

        for i in self._negate_pos:
            x[..., i] *= -1

        for i in self._double_inf_pos:
            # For (-oo, oo) -> (-1, 1), use the transformation x = (1-|t|)/t.
            x[..., i] = (1 - self._xp.abs(t[..., i])) / t[..., i]
            jacobian *= 1/(t[..., i] ** 2)

        for i in self._semi_inf_pos:
            # For (start, oo) -> (0, 1), use the transformation x = start + (1-t)/t.
            x[..., i] = self._orig_a[i] + (1 - t[..., i]) / t[..., i]
            jacobian *= 1/(t[..., i] ** 2)

        f_x = self._f(x, *args, **kwargs)
        jacobian = self._xp.reshape(jacobian, (-1, *([1]*(len(f_x.shape)-1))))

        return f_x * jacobian


class _FuncLimitsTransform(_VariableTransform):
    r"""
    This transformation handles integrals of the form:

    .. math::

        \int^{b_1}_{a_1}
        \cdots
        \int^{b_n}_{a_n}
        \int^{B_{n+1}(x_1, \ldots, x_n)}_{A_{n+1}(x_1, \ldots, x_n)}
        \cdots
        \int^{B_{n+m}(x_1, \ldots, x_{n+m-1})}_{A_{n+m}(x_1, \ldots, x_{n+m-1})}
            f(x_1, \ldots, x_{n+m})
        dx_{n+m} \cdots dx_1

    :math:`x_1, \ldots, x_n` are referred to as the "outer variables" and their limits
    are specified by the constant intervals :math:`(a_i, b_i)` for
    :math:`1 \le i \le n`.

    :math:`x_{n+1}, \ldots, x_{n+m}` are referred to as the "inner variables" and their
    limits given the values of the preceding variables range over the interval
    :math:`(A_i(x_1, \ldots, x_{i-1}), B_i(x_1, \ldots, x_{i-1})`.

    Integrals of this form will be transformed into an integral over constant limits:

    .. math::

        \int^{b_1}_{a_1}
        \cdots
        \int^{b_n}_{a_n}
        \int^{1}_{-1}
        \cdots
        \int^{1}_{-1}
            g(t_1, \ldots, t_{n+m})
        dt_{n+m} \cdots dt_1

    The outer variables are unchanged, :math:`x_i = t_i` for :math:`1 \le x_i \le t_i`.

    The inner variables are mapped via the transformation:

    .. math::

        x_i = \frac{
            B_{i}(x_1, \ldots, x_{i-1}) + A_{i}(x_1, \ldots, x_{i-1})
        }{2} + t_i \frac{
            B_{i}(x_1, \ldots, x_{i-1}) - A_{i}(x_1, \ldots, x_{i-1})
        }{2}

    for :math:`i > n`. This transformation has Jacobian determinant:

    .. math::

        J(x_1, \ldots, x_{n+m})
        =
        \prod^{n+m}_{i = n+1}
        \frac{B_{i}(x_1, \ldots, x_{i-1}) - A_{i}(x_1, \ldots, x_{i-1})}{2}

    To specify this type of integral programatically, the limits should be specified as
    two arrays, ``a`` and ``b`` for the constant limits, and a list of callables
    ``region`` describing the function limits::

        a = [a_1, ..., a_n]
        b = [b_1, ..., b_n]

        region = [
            region_func_1,
            ...,
            region_func_k,
        ]

    Each ``region_func_i`` should accepts arrays of shape
    ``(npoints, num_preceding_variables)`` and return a tuple of two arrays, one for
    the next set of lower limits given the values of the preceding variables, and one
    for the next set of upper limits given the values of the preceding variables. Both
    of the returned arrays should have shape ``(npoints, num_new_variables)``.
    """

    def __init__(self, f, a, b, region):
        self._xp = array_namespace(a, b)

        self._f = f
        self._a_outer = a
        self._b_outer = b

        self._region = region

        # The "outer dimension" here is the number of non-function limits.
        self._outer_ndim = xp_size(self._a_outer)

        # The "inner dimension" here is the number of limits returned by the functions
        # in `region`.
        #
        # Without evaluating each of the functions in `region` at least once, it's
        # impossible to know the number of inner variables. Knowing this is required to
        # return the transformed limits.
        #
        # To evaluate the functions in `region` at least once, it's necessary to have an
        # initial assignment to the outer variables. Here `a` is picked.
        full_limits = self._xp.reshape(self._a_outer, (1, -1))

        for region_func in self._region:
            full_limits = self._xp.concat(
                [
                    full_limits,

                    # Evaluate the next `region_func` on all of the variables so far
                    # and take the lower limit of the result as the value of the next
                    # set of inner variables.
                    region_func(full_limits)[0],
                ],
                axis=-1,
            )

        self._inner_ndim = full_limits.shape[-1] - xp_size(self._a_outer)

    @property
    def transformed_limits(self):
        return (
            self._xp.concat([self._a_outer, -self._xp.ones(self._inner_ndim)]),
            self._xp.concat([self._b_outer,  self._xp.ones(self._inner_ndim)]),
        )

    def _inner_limits(self, t):
        # Given this value of `t`, calculate the values of the functions describing the
        # inner limits.
        x = xp_copy(t)

        a_inner = []
        b_inner = []

        # Each region_func depends on the values of variables preceding it, this keeps
        # track of which variables need to be sent to the next region_func to return
        # a new group of limits.
        num_preceding_vars = self._outer_ndim

        for region_func in self._region:
            # region_func will return a tuple (A_group, B_group), where A_group and
            # B_group are arrays of shape (npoints, num_new_vars) specifying the limits
            # for the next group of variables.

            A_group, B_group = region_func(x[..., :num_preceding_vars])
            num_new_vars = A_group.shape[-1]

            t_group = t[:, num_preceding_vars:num_preceding_vars+num_new_vars]

            # Apply transformation to this group of variables.
            x[..., num_preceding_vars:num_preceding_vars+num_new_vars] = (
                (B_group + A_group)/2 + t_group * (B_group - A_group)/2
            )

            a_inner.append(A_group)
            b_inner.append(B_group)

            # The next region_func will need to know the value of the variables that
            # have just been returned.
            num_preceding_vars += num_new_vars

        return self._xp.concat(a_inner, axis=-1), self._xp.concat(b_inner, axis=-1)

    def inv(self, x):
        t = xp_copy(x)

        # Each region_func depends on the values of variables preceding it, this keeps
        # track of which variables need to be sent to the next region_func to return
        # a new group of limits.
        num_preceding_vars = self._outer_ndim

        for region_func in self._region:
            # region_func will return a tuple (A_group, B_group), where A_group and
            # B_group are arrays of shape (npoints, num_new_vars) specifying the limits
            # for the next group of variables.

            A_group, B_group = region_func(x[..., :num_preceding_vars])
            num_new_vars = A_group.shape[-1]

            x_group = x[..., num_preceding_vars:num_preceding_vars+num_new_vars]

            # Apply inverse transformation to this group of variables.
            t[..., num_preceding_vars:num_preceding_vars+num_new_vars] = (
                (2*x_group - B_group - A_group)/(B_group - A_group)
            )

            # The next region_func will need to know the value of the variables that
            # have just been returned.
            num_preceding_vars += num_new_vars

        return t

    def __call__(self, t, *args, **kwargs):
        # Inner limits are the values of the function limits given this value of `t`.
        a_inner, b_inner = self._inner_limits(t)

        # Allow returning `a_inner` and `b_inner` as array_like rather than as ndarrays
        # since `a` and `b` can also be array_like.
        a_inner = self._xp.asarray(a_inner)
        b_inner = self._xp.asarray(b_inner)

        # Prepare the input to the original integrand f.
        # `x[i] = t[i]` where `i <= self._inner_ndim`, but we need to map but for the
        # rest, use the transformation described in the class-level docstring.
        x = xp_copy(t)

        for i in range(self._inner_ndim):
            a_i = a_inner[:, i]
            b_i = b_inner[:, i]

            x[:, self._outer_ndim + i] = (
                (b_i + a_i)/2 + (b_i - a_i)/2 * t[:, self._outer_ndim + i]
            )

        f_x = self._f(x, *args, **kwargs)
        jacobian = self._xp.prod(b_inner - a_inner, axis=-1) / 2**self._inner_ndim
        jacobian = self._xp.reshape(jacobian, (-1, *([1]*(len(f_x.shape) - 1))))

        return f_x * jacobian
