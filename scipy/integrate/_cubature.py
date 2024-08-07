import heapq
import itertools

from dataclasses import dataclass

import numpy as np

from scipy.integrate._rules import (
    ProductNestedFixed, NestedFixedRule,
    GaussKronrodQuad, NewtonCotesQuad,
)


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
        return _max_norm(self.error) > _max_norm(other.error)


@dataclass
class CubatureResult:
    estimate: np.ndarray
    error: np.ndarray
    success: bool
    status: str
    regions: list[CubatureRegion]
    subdivisions: int
    atol: float
    rtol: float


def cub(f, a, b, rule="gk21", rtol=1e-05, atol=1e-08, max_subdivisions=10000,
        args=(), kwargs=None):
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
            f(x : ndarray, \*args, \*\*kwargs) -> ndarray

        `f` should accept arrays ``x`` of shape::
            (npoints, ndim)

        and output arrays of shape::
            (npoints, output_dim_1, ..., output_dim_n)

        In this case, `cub` will return arrays of shape::
            (output_dim_1, ..., output_dim_n)
    a, b : array_like or float
        Lower and upper limits of integration as rank-1 arrays specifying the left and
        right endpoints of the intervals being integrated over. If a float is passed,
        these will be converted to singleton arrays. Infinite limits are currently not
        supported.
    rule : str or instance of ``Rule``, optional
        Cubature rule to use when estimating the integral. If passing a string, the
        options are "gk21", "gk15", or "trapezoid". Default is "gk21". If passing an
        instance of a ``Rule``, then the options can be found in
        ``scipy.integrate._rules``.
    rtol : float, optional
        Relative tolerance. Default is 1e-05.
    atol : float, optional
        Absolute tolerance. Default is 1e-08.
    max_subdivisions : int, optional
        Upper bound on the number of subdivisions to perform to improve the estimate
        over a subregion. Default is 10,000.
    args : tuple, optional
        Additional positional args passed to `f`, if any.
    kwargs : tuple, optional
        Additional keyword args passed to `f`, if any.

    Returns
    -------
    res : CubatureResult
        Result of estimation. See `CubatureResult`.

    Examples
    --------
    A simple 1D integral with vector output. Here ``f(x) = x^n`` is integrated over the
    interval ``[0, 1]``. Since no rule is specified, the default "gk21" is used, which
    corresponds to `GaussKronrod` rule with 21 nodes.

    >>> import numpy as np
    >>> from scipy.integrate import cub
    >>> def f(x, n):
    ...    return x.reshape(-1, 1)**n  # Make sure x and n are broadcastable
    >>> res = cub(
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
    hybercube, ``[0, 1]^7``.

    >>> import numpy as np
    >>> from scipy.integrate import cub
    >>> from scipy.integrate._rules import GenzMalikCub
    >>> def f(x, r, alphas):
    ...     # f(x) = cos(2*pi*r + alpha @ x)
    ...     # Need to allow r and alphas to be arbitrary shape
    ...     npoints, ndim = x.shape[0], x.shape[-1]
    ...     alphas_reshaped = alphas[np.newaxis, :]
    ...     x_reshaped = x.reshape(npoints, *([1]*(len(alphas.shape) - 1)), ndim)
    ...     return np.cos(2*np.pi*r + np.sum(alphas_reshaped * x_reshaped, axis=-1))
    >>> np.random.seed(1)
    >>> res = cub(
    ...     f=f,
    ...     a=np.array([0, 0, 0, 0, 0, 0, 0]),
    ...     b=np.array([1, 1, 1, 1, 1, 1, 1]),
    ...     rule=GenzMalikCub(7),
    ...     kwargs={
    ...         "r": np.random.rand(2, 3),
    ...         "alphas": np.random.rand(2, 3, 7),
    ...     }
    ... )
    >>> res.estimate
     array([[-0.61336595,  0.88388877, -0.57997549],
            [-0.86968418, -0.86877137, -0.64602074]])

    It is also possible to use a custom rule. In the following, a custom rule is created
    which uses 3D Genz-Malik cubature for the estimate of the integral, and the
    difference between this estimate and a less accurate estimate using 5-node
    Gauss-Legendre quadrature as an estimate for the error.

    >>> import numpy as np
    >>> from scipy.integrate import cub
    >>> from scipy.integrate._rules import (
    ...     Rule, ProductFixed, GenzMalikCub, GaussLegendreQuad
    ... )
    >>> genz = GenzMalikCub(3)
    >>> gauss = GaussLegendreQuad(5)
    >>> # Gauss-Legendre is 1D, so we find the 3D product rule:
    >>> gauss_3d = ProductFixed([gauss, gauss, gauss])
    >>> class CustomRule(Rule):
    ...     def estimate(self, f, a, b, args=(), kwargs=None):
    ...         kwargs = kwargs or {}
    ...         return genz.estimate(f, a, b, args, kwargs)
    ...     def estimate_error(self, f, a, b, args=(), kwargs=None):
    ...         kwargs = kwargs or {}
    ...         return np.abs(
    ...             genz.estimate(f, a, b, args, kwargs)
    ...             - gauss_3d.estimate(f, a, b, args, kwargs)
    ...         )
    >>> np.random.seed(1)
    >>> cub(
    ...     f=f,
    ...     a=np.array([0, 0, 0]),
    ...     b=np.array([1, 1, 1]),
    ...     rule=CustomRule(),
    ...     kwargs={
    ...         "r": np.random.rand(2),
    ...         "alphas": np.random.rand(3, 2, 3),
    ...     }
    ... ).estimate
     array([[-0.95179502,  0.12444608],
            [-0.96247411,  0.60866385],
            [-0.97360014,  0.25515587]])

    This particular example estimates the error using the difference between two
    approximations, one more accurate than the other. These are called nested rules
    and can be created using ``NestedRule``:

    >>> from scipy.integrate._rules import NestedRule
    >>> rule = NestedRule(
    ...     higher=genz,
    ...     lower=gauss_3d,
    ... ) # Equivalent to CustomRule()
    """

    kwargs = kwargs or {}
    max_subdivisions = np.inf if max_subdivisions is None else max_subdivisions

    # Convert a and b to arrays of at least 1D
    a = np.atleast_1d(a)
    b = np.atleast_1d(b)

    # If the rule is a string, convert to a corresponding product rule
    if isinstance(rule, str):
        quadratues = {
            "gk21": GaussKronrodQuad(21),
            "gk15": GaussKronrodQuad(15),
            "trapezoid": NestedFixedRule(NewtonCotesQuad(5), NewtonCotesQuad(3)),
        }

        base_quadrature = quadratues.get(rule)

        if base_quadrature is None:
            raise ValueError(f"unknown rule {rule}")

        rule = ProductNestedFixed([base_quadrature] * len(a))

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("a and b should be 1D arrays")

    est = rule.estimate(f, a, b, args, kwargs)

    try:
        err = rule.estimate_error(f, a, b, args, kwargs)
    except NotImplementedError:
        raise ValueError("attempting cubature with a rule that doesn't implement error"
                         "estimation.")

    regions = [CubatureRegion(est, err, a, b)]
    subdivisions = 0
    success = True

    while np.any(err > atol + rtol * np.abs(est)):
        # region_k is the region with highest estimated error
        region_k = heapq.heappop(regions)

        est_k = region_k.estimate
        err_k = region_k.error

        a_k, b_k = region_k.a, region_k.b

        # Subtract the estimate of the integral and its error over this region from the
        # current global estimates, since these will be refined in the loop over all
        # subregions.
        est -= est_k
        err -= err_k

        # Find all 2^ndim subregions formed by splitting region_k along each axis, e.g.
        # for 1D integrals this splits an estimate over an interval into an estimate
        # over two subintervals, for 3D integrals this splits an estimate over a cube
        # into 8 subcubes.
        #
        # For each of the new subregions, calculate an estimate for the integral and
        # the error there, and push these regions onto the heap for potential further
        # subdividing.
        for a_k_sub, b_k_sub in _subregion_coordinates(a_k, b_k):
            est_sub = rule.estimate(f, a_k_sub, b_k_sub, args, kwargs)
            err_sub = rule.estimate_error(f, a_k_sub, b_k_sub, args, kwargs)

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
        success=success,
        status=status,
        subdivisions=subdivisions,
        regions=regions,
        atol=atol,
        rtol=rtol,
    )


def _subregion_coordinates(a, b):
    """
    Given the coordinates of a region like a=[0, 0] and b=[1, 1], yield the coordinates
    of all subregions, which in this case would be::

        ([0, 0], [1/2, 1/2]),
        ([0, 1/2], [1/2, 1]),
        ([1/2, 0], [1, 1/2]),
        ([1/2, 1/2], [1, 1])
    """

    m = (a + b)/2

    for a_sub, b_sub in zip(
        itertools.product(*np.array([a, m]).T),
        itertools.product(*np.array([m, b]).T)
    ):
        yield np.array(a_sub), np.array(b_sub)


def _max_norm(x):
    return np.max(np.abs(x))
