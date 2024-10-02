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
from scipy.integrate._rules._base import _split_subregion

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
    _xp: ModuleType = field(repr=False)

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


def cubature(f, a, b, *, rule="gk21", rtol=1e-8, atol=0, max_subdivisions=10000,
             args=(), workers=1, points=None):
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
        regions: list of object
            List of objects containing the estimates of the integral over smaller
            regions of the domain.

        Each object in ``regions`` has the following attributes:

        a, b : ndarray
            Points describing the corners of the region. If the original integral
            contained infinite limits or was over a region described by `region`,
            then `a` and `b` are in the transformed coordinates.
        estimate : ndarray
            Estimate of the value of the integral over this region.
        error : ndarray
            Estimate of the error of the approximation over this region.

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
        ISSN 0377-0427
        :doi:`10.1016/0771-050X(80)90039-X`

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

    # If any problematic points are specified, divide the initial region so that these
    # points lie on the edge of a subregion.
    #
    # This means ``f`` won't be evaluated there if the rule being used has no evaluation
    # points on the boundary.
    if points == []:
        initial_regions = [(a, b)]
    else:
        initial_regions = _split_region_at_points(a, b, points)

    regions = []
    est = 0.0
    err = 0.0

    for a_k, b_k in initial_regions:
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
                _split_subregion(a_k, b_k),
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

        # Apply sign change to handle any limits which were initially flipped.
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


def _split_region_at_points(a, b, points):
    """
    Given the integration limits `a` and `b` describing a rectangular region and a list
    of `points`, find the list of ``[(a_1, b_1), ..., (a_l, b_l)]`` which breaks up the
    initial region into smaller subregion such that no `points` lie strictly inside
    any of the subregions.
    """

    xp = array_namespace(a, b)
    regions = [(a, b)]

    for point in points:
        if xp.any(xp.isinf(point)):
            # If a point is specified at infinity, ignore.
            #
            # This case occurs when points are given by the user to avoid, but after
            # applying a transformation, they are removed.
            continue

        new_subregions = []

        for a_k, b_k in regions:
            if _is_strictly_in_region(point, a_k, b_k):
                subregions = _split_subregion(a_k, b_k, point)

                for left, right in subregions:
                    # Skip any zero-width regions.
                    if xp.any(left == right):
                        continue
                    else:
                        new_subregions.append((left, right))

                new_subregions.extend(subregions)

            else:
                new_subregions.append((a_k, b_k))

        regions = new_subregions

    return regions


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
