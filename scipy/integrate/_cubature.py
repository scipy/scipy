import heapq
import itertools

import numpy as np  # TODO: remove

from dataclasses import dataclass
from types import ModuleType
from typing import Any, TYPE_CHECKING

from scipy._lib._array_api import array_namespace, xp_size
from scipy._lib._util import MapWrapper

from scipy.integrate._rules import (
    ProductNestedFixed,
    GaussKronrodQuadrature,
    GenzMalikCubature,
)
from scipy.integrate._rules._base import _subregion_coordinates

__all__ = ['cubature']

if TYPE_CHECKING:
    Array = Any  # To be changed to a Protocol later (see array-api#589)
else:
    Array = object


@dataclass
class CubatureRegion:
    estimate: Array
    error: Array
    a: Array
    b: Array
    _xp: ModuleType

    def __lt__(self, other):
        # Consider regions with higher error estimates as being "less than" regions with
        # lower order estimates, so that regions with high error estimates are placed at
        # the top of the heap.

        this_err = self._xp.max(self._xp.abs(self.error))
        other_err = self._xp.max(self._xp.abs(other.error))

        return this_err > other_err


@dataclass
class CubatureResult:
    estimate: Array
    error: Array
    status: str
    regions: list[CubatureRegion]
    subdivisions: int
    atol: float
    rtol: float


def cubature(f, a, b, rule="gk21", rtol=1e-8, atol=0, max_subdivisions=10000,
             args=(), workers=1, points=None, region=None):
    r"""
    Adaptive cubature of multidimensional array-valued function.

    Given an arbitrary integration rule, this function returns an estimate of the
    integral to the requested tolerance over the region defined by the arrays `a` and
    `b` specifying the corners of a hypercube.

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
        endpoints of the intervals being integrated over. Infinite limits are currently
        not supported.
    rule : str, optional
        Rule used to estimate the integral. If passing a string, the options are
        "gauss-kronrod" (21 node), or "genz-malik" (degree 7). If a rule like
        "gauss-kronrod" is specified for an ``n``-dim integrand, the corresponding
        Cartesian product rule is used. See Notes.

        "gk21", "gk15" are also supported for compatibility with `quad_vec`.
    rtol, atol : float, optional
        Relative and absolute tolerances. Iterations are performed until the error is
        estimated to be less than ``atol + rtol * abs(est)``. Here `rtol` controls
        relative accuracy (number of correct digits), while `atol` controls absolute
        accuracy (number of correct decimal places). To achieve the desired `rtol`, set
        `atol` to be smaller than the smallest value that can be expected from
        ``rtol * abs(y)`` so that rtol dominates the allowable error. If `atol` is
        larger than ``rtol * abs(y)`` the number of correct digits is not guaranteed.
        Conversely, to achieve the desired `atol` set `rtol` such that ``rtol * abs(y)``
        is always smaller than `atol`. Default values are 1e-8 for `rtol` and 0 for
        `atol`.
    max_subdivisions : int, optional
        Upper bound on the number of subdivisions to perform to improve the estimate
        over a subregion. Default is 10,000.
    args : tuple, optional
        Additional positional args passed to `f`, if any.
    workers : int or map-like callable, optional
        If `workers` is an integer, part of the computation is done in parallel
        subdivided to this many tasks (using :class:`python:multiprocessing.pool.Pool`).
        Supply `-1` to use all cores available to the Process. Alternatively, supply a
        map-like callable, such as :meth:`python:multiprocessing.pool.Pool.map` for
        evaluating the population in parallel. This evaluation is carried out as
        ``workers(func, iterable)``.

    Returns
    -------
    res : CubatureResult
        Result of estimation. The estimate of the integral is `res.estimate`, and the
        estimated error is `res.error`. If the integral converges within
        `max_subdivions`, then `res.status` will be ``"converged"``, otherwise it will
        be ``"not_converged"``. See `CubatureResult`.

    Examples
    --------
    A simple 1D integral with vector output. Here ``f(x) = x^n`` is integrated over the
    interval ``[0, 1]``. Since no rule is specified, the default "gk21" is used, which
    corresponds to `GaussKronrod` rule with 21 nodes.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> def f(x, n):
    ...    return x.reshape(-1, 1)**n  # Make sure x and n are broadcastable
    >>> res = cubature(
    ...     f,
    ...     a=[0],
    ...     b=[1],
    ...     args=(
    ...         # Since f accepts arrays of shape (npoints, ndim) we need to
    ...         # make sure n is the right shape
    ...         np.arange(10).reshape(1, -1),
    ...     )
    ... )
    >>> res.estimate
     array([1.        , 0.5       , 0.33333333, 0.25      , 0.2       ,
            0.16666667, 0.14285714, 0.125     , 0.11111111, 0.1       ])

    A 7D integral with arbitrary-shaped array output. Here::

        f(x) = cos(2*pi*r + alphas @ x)

    for some ``r`` and ``alphas``, and the integral is performed over the unit
    hybercube, :math:`[0, 1]^7`. Since the integral is in a moderate number of
    dimensions, "genz-malik" is used rather than the default "gauss-kronrod" to avoid
    constructing a product rule with :math:`21^7 \approx 2 \times 10^9` nodes.

    >>> import numpy as np
    >>> from scipy.integrate import cubature
    >>> def f(x, r, alphas):
    ...     # f(x) = cos(2*pi*r + alpha @ x)
    ...     # Need to allow r and alphas to be arbitrary shape
    ...     npoints, ndim = x.shape[0], x.shape[-1]
    ...     alphas_reshaped = alphas[np.newaxis, :]
    ...     x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)
    ...     return np.cos(2*np.pi*r + np.sum(alphas_reshaped * x_reshaped, axis=-1))
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

    To compute in parallel, it is possible to use the argument `workers`, for example:

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

    Notes
    -----
    The algorithm uses a similar algorithm to `quad_vec`, which itself is based on the
    implementation of QUADPACK's DQAG* algorithms, implementing global error control and
    adaptive subdivision.

    The source of the nodes and weights used for Gauss-Kronrod quadrature can be found
    in [1]_, and the algorithm for calculating the nodes and weights in Genz-Malik
    cubature can be found in [2]_.

    The rules currently supported via the `rule` argument are:

    - ``"gauss-kronrod"``, 21-node Gauss-Kronrod
    - ``"genz-malik"``, n-node Genz-Malik

    If using Gauss-Kronrod for an ``n``-dim integrand where ``n > 2``, then the
    corresponding Cartesian product rule will be found by taking the Cartesian product
    of the nodes in the 1D case. This means that the number of nodes scales
    exponentially as ``21^n`` in the Gauss-Kronrod case, which may be problematic in a
    moderate number of dimensions.

    Genz-Malik is typically less accurate than Gauss-Kronrod but has much fewer nodes,
    so in this situation using "genz-malik" might be preferable.

    References
    ----------
    .. [1] R. Piessens, E. de Doncker, Quadpack: A Subroutine Package for Automatic
        Integration, files: dqk21.f, dqk15.f (1983).

    .. [2] A.C. Genz, A.A. Malik, Remarks on algorithm 006: An adaptive algorithm for
        numerical integration over an N-dimensional rectangular region, Journal of
        Computational and Applied Mathematics, Volume 6, Issue 4, 1980, Pages 295-302,
        ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90039-X.
    """

    # It is also possible to use a custom rule, but this is not yet part of the public
    # API. An example of this can be found in the class scipy.integrate._rules.Rule.

    xp = array_namespace(a, b)
    max_subdivisions = float("inf") if max_subdivisions is None else max_subdivisions
    points = [] if points is None else [xp.asarray(p, dtype=xp.float64) for p in points]

    # Convert a and b to arrays
    a = xp.asarray(a, dtype=xp.float64)
    b = xp.asarray(b, dtype=xp.float64)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("`a` and `b` must be 1D arrays")

    # If `region` is specified, apply a transformation to reduce to integration over
    # a rectangular region.
    if region is not None:
        f = _FuncLimitsTransform(f, a, b, region)
        a, b = f.limits

    # If any of limits are the wrong way around (a > b), flip them and keep track of
    # the sign.
    sign = (-1.0) ** xp.sum(xp.astype(a > b, xp.float64))
    a_flipped = xp.min(xp.asarray([a, b]), axis=0)
    b_flipped = xp.max(xp.asarray([a, b]), axis=0)

    a, b = a_flipped, b_flipped

    # If any of the limits are infinite, apply a transformation to handle these
    if np.any(np.isinf(a_flipped)) or np.any(np.isinf(b_flipped)):
        f = _InfiniteLimitsTransform(f, a_flipped, b_flipped)
        a, b = f.limits
        points.extend(f.points)

    # If any problematic points are specified, divide the initial region so that
    # these points lie on the edge of a subregion.
    # This means ``f`` won't be evaluated there, if the rule being used has no
    # evaluation points on the boundary.
    if points == []:
        initial_regions = [(a, b)]
    else:
        initial_regions = _split_at_points(a, b, points)

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

    regions = []
    est = 0.0
    err = 0.0

    # Compute an estimate over each of the initial regions
    for a_k, b_k in initial_regions:
        # If any of the initial regions have zero width in one dimension, we can
        # ignore this as the integral will be 0 there.
        if not np.any(a_k == b_k):
            est_k = rule.estimate(f, a_k, b_k, args)
            err_k = rule.estimate_error(f, a_k, b_k, args)
            regions.append(CubatureRegion(est_k, err_k, a_k, b_k, xp))

            est += est_k
            err += err_k

    subdivisions = 0
    success = True

    with MapWrapper(workers) as mapwrapper:
        while xp.any(err > atol + rtol * xp.abs(est)):
            # region_k is the region with highest estimated error
            region_k = heapq.heappop(regions)

            est_k = region_k.estimate
            err_k = region_k.error

            a_k, b_k = region_k.a, region_k.b

            # Subtract the estimate of the integral and its error over this region from
            # the current global estimates, since these will be refined in the loop over
            # all subregions.
            est -= est_k
            err -= err_k

            # Find all 2^ndim subregions formed by splitting region_k along each axis,
            # e.g. for 1D integrals this splits an estimate over an interval into an
            # estimate over two subintervals, for 3D integrals this splits an estimate
            # over a cube into 8 subcubes.
            #
            # For each of the new subregions, calculate an estimate for the integral and
            # the error there, and push these regions onto the heap for potential
            # further subdividing.

            executor_args = zip(
                itertools.repeat(f),
                itertools.repeat(rule),
                itertools.repeat(args),
                _subregion_coordinates(a_k, b_k),
            )

            for subdivision_result in mapwrapper(_process_subregion, executor_args):
                a_k_sub, b_k_sub, est_sub, err_sub = subdivision_result

                est += est_sub
                err += err_sub

                new_region = CubatureRegion(est_sub, err_sub, a_k_sub, b_k_sub, xp)

                heapq.heappush(regions, new_region)

            subdivisions += 1

            if subdivisions >= max_subdivisions:
                success = False
                break

        status = "converged" if success else "not_converged"

        # Apply sign change if any of the limits were initially flipped
        est = sign * est

        return CubatureResult(
            estimate=est,
            error=err,
            status=status,
            subdivisions=subdivisions,
            regions=regions,
            atol=atol,
            rtol=rtol,
        )


def _process_subregion(data):
    f, rule, args, coord = data
    a_k_sub, b_k_sub = coord

    est_sub = rule.estimate(f, a_k_sub, b_k_sub, args)
    err_sub = rule.estimate_error(f, a_k_sub, b_k_sub, args)

    return a_k_sub, b_k_sub, est_sub, err_sub


def _is_strictly_in_region(point, a, b):
    xp = array_namespace(point, a, b)

    if xp.all(point == a) or xp.all(point == b):
        return False

    return xp.all(a <= point) and xp.all(point <= b)


def _split_at_points(a, b, points):
    """
    Given the integration limits `a` and `b` describing a rectangular region and a list
    of `points`, find the list of ``[(a_1, b_1), ..., (a_l, b_l)]`` which breaks up the
    initial region into smaller subregions such that no `points` lie strictly inside
    any of the subregions.
    """

    regions = [(a, b)]

    for point in points:
        new_regions = []

        for a_k, b_k in regions:
            if _is_strictly_in_region(point, a_k, b_k):
                subregions = _subregion_coordinates(a_k, b_k, point)
                new_regions.extend(subregions)
            else:
                new_regions.append((a_k, b_k))

        regions = new_regions

    return regions


class _VariableTransform:
    @property
    def limits(self):
        """
        New limits after applying the transformation.
        """

        raise NotImplementedError

    @property
    def points(self):
        """
        Any problematic points introduced by the transformation to avoid evaluating.
        """

        return []

    def __call__(self, t, *args, **kwargs):
        """
        `f` after the transformation.
        """

        raise NotImplementedError


class _InfiniteLimitsTransform(_VariableTransform):
    r"""
    Given ``f`` and the limits ``a`` and ``b``, apply a transformation to ``f`` to map
    any infinite limits to a finite interval.

    `math:`[-\infty, \infty]` is mapped to `math:`(-1, 1)` using the transformation
   `math:`x = (1-|t|)/t`.

    `math:`[a, \infty]` is mapped to `math:`(0, 1)` using the transformation
    `math:`x = a + (1-t)/t`.

    `math:`[-\infty, b]` is mapped to `math:`(0, 1)` using the transformation
    `math:`x = a - (1+t)/t`.
    """

    def __init__(self, f, a, b):
        self._f = f
        self._orig_a = a
        self._orig_b = b

        self._negate_pos = []
        self._semi_inf_pos = []
        self._double_inf_pos = []

        for i in range(len(a)):
            if a[i] == -np.inf and b[i] == np.inf:
                # (-oo, oo) will be mapped to (-1, 1).
                self._double_inf_pos.append(i)

            elif a[i] != -np.inf and b[i] == np.inf:
                # (start, oo) will be mapped to (0, 1).
                self._semi_inf_pos.append(i)

            elif a[i] == -np.inf and b[i] != np.inf:
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

        self._semi_inf_pos = np.array(self._semi_inf_pos)
        self._double_inf_pos = np.array(self._double_inf_pos)
        self._negate_pos = np.array(self._negate_pos)

    @property
    def limits(self):
        a, b = np.copy(self._orig_a), np.copy(self._orig_b)

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
        if self._double_inf_pos.size != 0 or self._semi_inf_pos.size != 0:
            return [np.zeros(self._orig_a.shape)]
        else:
            return []

    def __call__(self, t, *args, **kwargs):
        x = np.copy(t)

        if self._negate_pos.size != 0:
            x[..., self._negate_pos] *= -1

        if self._double_inf_pos.size != 0:
            # For (-oo, oo) -> (-1, 1), use the transformation x = (1-|t|)/t.

            x[..., self._double_inf_pos] = (
                1 - np.abs(t[..., self._double_inf_pos])
            ) / t[..., self._double_inf_pos]

        if self._semi_inf_pos.size != 0:
            # For (start, oo) -> (0, 1), use the transformation x = start + (1 - t)/t.

            # Need to expand start so it is broadcastable, and transpose to flip the
            # axis order.
            start = self._orig_a[self._semi_inf_pos][:, np.newaxis].T

            x[..., self._semi_inf_pos] = start + (
                1 - t[..., self._semi_inf_pos]
            ) / t[..., self._semi_inf_pos]

        f_x = self._f(x, *args, **kwargs)
        jacobian = 1

        if self._double_inf_pos.size != 0:
            jacobian *= 1/np.prod(
                t[..., self._double_inf_pos] ** 2,
                axis=-1,
            )

        if self._semi_inf_pos.size != 0:
            jacobian *= 1/np.prod(
                t[..., self._semi_inf_pos] ** 2,
                axis=-1,
            )

        jacobian = jacobian.reshape(-1, *([1]*(len(f_x.shape)-1)))

        return f_x * jacobian


class _FuncLimitsTransform(_VariableTransform):
    r"""
    Transform an integral with functions as limits to an integral with constant limits.

    Given an integral of the form:

    ..math ::

        \int^{b_1}_{a_1}
        \cdots
        \int^{b_n}_{a_n}
        \int^{B_1(x_1, \ldots, x_n)}_{A_1(x_1, \ldots, x_n)}
        \cdots
        \int^{B_m(x_1, \ldots, x_{n+m})}_{A_m(x_1, \ldots, x_{n+m})}
        f(x_1, \ldots, x_{n+m})
        dx_{n+m} \cdots dx_1

    an integral with :math:`n` outer non-function limits, and :math:`m` inner function
    limits, this will transform it into an integral over

    ..math::

        \int^{b_1}_{a_1}
        \cdots
        \int^{b_n}_{a_n}
        \int^{1}_{-1}
        \cdots
        \int^{1}_{-1}
        g(x_1, \ldots, x_n, y_1, \cdots, y_m)
        dy_m \cdots dy_1 dx_n \cdots dx_1

    Which is an integral over the original outer non-function limits and where a
    transformation has been applied so that the original function limits become [-1, 1].
    """

    def __init__(self, f, a, b, region):
        self._f = f
        self._a_outer = a
        self._b_outer = b

        self._region = region

        # The "outer dimension" here is the number of the outer non-function limits.
        self._outer_ndim = len(self._a_outer)

        # Without evaluating `region` at least once, it's impossible to know the
        # number of inner variables, which is required to return new limits.
        # TODO: don't evaluate at boundary
        limits = self._a_outer.reshape(1, -1)

        for region_func in self._region:
            limits = np.concat([limits, region_func(limits)[0]], axis=-1)

        self._inner_ndim = np.array(limits).shape[-1] - len(self._a_outer)

    # TODO: rename
    def _get_inner_limits(self, x_and_t):
        x = x_and_t[:, :self._outer_ndim]
        t = x_and_t[:, self._outer_ndim:]

        a_inner, b_inner = np.empty((x.shape[0], 0)), np.empty((x.shape[0], 0))

        A_i, B_i = None, None
        x_and_y = x
        y_i = np.array([])

        for i, region_func in enumerate(self._region):
            A_i, B_i = region_func(x_and_y)

            outer_ndim = a_inner.shape[-1]
            inner_ndim = A_i.shape[-1]

            y_i = (B_i + A_i)/2 + t[:, outer_ndim:outer_ndim+inner_ndim] * (B_i - A_i)/2

            x_and_y = np.concat([x_and_y, y_i], axis=-1)
            a_inner = np.concat([a_inner, A_i], axis=-1)
            b_inner = np.concat([b_inner, B_i], axis=-1)

        return a_inner, b_inner

    @property
    def limits(self):
        return (
            np.concatenate([self._a_outer, -np.ones(self._inner_ndim)]),
            np.concatenate([self._b_outer,  np.ones(self._inner_ndim)]),
        )

    def __call__(self, x_and_t, *args, **kwargs):
        # x_and_t consists of the outer variables x_1, ... x_n, which don't need
        # changing, and then inner variables t_1, ..., t_n which are in the range
        # [-1, 1].

        a_inner, b_inner = self._get_inner_limits(x_and_t)

        # Allow returning `a_inner` and `b_inner` as array_like rather than as ndarrays
        # since `a` and `b` can also be array_like.
        a_inner = np.array(a_inner)
        b_inner = np.array(b_inner)

        npoints = x_and_t.shape[0]

        # x_and_y should be the input to the original integrand f.
        # No change needed to x, but we need to map the t section of x_and_t back to y.
        x_and_y = np.concatenate(
            [
                x_and_t[:, :self._outer_ndim],
                np.zeros((npoints, self._inner_ndim)),
            ],
            axis=-1,
        )

        for i in range(self._inner_ndim):
            a_i = a_inner[:, i]
            b_i = b_inner[:, i]

            x_and_y[:, self._outer_ndim + i] = (
                (b_i + a_i)/2 + (b_i - a_i)/2 * x_and_t[:, self._outer_ndim + i]
            )

        f_x = self._f(x_and_y, *args, **kwargs)
        jacobian = np.prod(b_inner - a_inner, axis=-1) / 2**self._inner_ndim
        jacobian = jacobian.reshape(-1, *([1]*(len(f_x.shape) - 1)))

        return f_x * jacobian
