import heapq
import itertools

from dataclasses import dataclass
from scipy._lib._util import MapWrapper

import numpy as np

from scipy.integrate._rules import (
    ProductNestedFixed,
    GaussKronrodQuadrature,
    GenzMalikCubature,
)
from scipy.integrate._rules._base import _subregion_coordinates


__all__ = ['cubature']


@dataclass
class CubatureRegion:
    estimate: np.ndarray
    error: np.ndarray
    a: np.ndarray
    b: np.ndarray

    def __lt__(self, other):
        # Consider regions with higher error estimates as being "less than" regions with
        # lower order estimates, so that regions with high error estimates are placed at
        # the top of the heap.

        return np.max(np.abs(self.error)) > np.max(np.abs(other.error))


@dataclass
class CubatureResult:
    estimate: np.ndarray
    error: np.ndarray
    status: str
    regions: list[CubatureRegion]
    subdivisions: int
    atol: float
    rtol: float


def cubature(f, a, b, rule="gk21", rtol=1e-8, atol=0, max_subdivisions=10000,
             args=(), workers=1):
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
            f(x : ndarray, \*args) -> ndarray

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
    """

    # It is also possible to use a custom rule, but this is not yet part of the public
    # API. An example of this can be found in the class scipy.integrate._rules.Rule.

    max_subdivisions = np.inf if max_subdivisions is None else max_subdivisions

    # Convert a and b to arrays of at least 1D
    a = np.array(a)
    b = np.array(b)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("`a` and `b` must be 1D arrays")

    # If the rule is a string, convert to a corresponding product rule
    if isinstance(rule, str):
        ndim = len(a)

        if rule == "genz-malik":
            rule = GenzMalikCubature(ndim)
        else:
            quadratues = {
                "gauss-kronrod": GaussKronrodQuadrature(21),

                # Also allow names quad_vec uses:
                "gk21": GaussKronrodQuadrature(21),
                "gk15": GaussKronrodQuadrature(15),
            }

            base_rule = quadratues.get(rule)

            if base_rule is None:
                raise ValueError(f"unknown rule {rule}")

            rule = ProductNestedFixed([base_rule] * ndim)

    est = rule.estimate(f, a, b, args)
    err = rule.estimate_error(f, a, b, args)

    regions = [CubatureRegion(est, err, a, b)]
    subdivisions = 0
    success = True

    with MapWrapper(workers) as mapwrapper:
        while np.any(err > atol + rtol * np.abs(est)):
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

                new_region = CubatureRegion(est_sub, err_sub, a_k_sub, b_k_sub)

                heapq.heappush(regions, new_region)

            subdivisions += 1

            if subdivisions >= max_subdivisions:
                success = False
                break

        status = "converged" if success else "not_converged"

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
