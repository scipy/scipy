import heapq
import itertools

from dataclasses import dataclass

import numpy as np

from scipy.integrate._rules import (
    FixedProductErrorFromDifferenceCub, ErrorFromDifference,
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

    Given an arbitrary cubature rule, this function returns an estimate of the integral
    over the defined by the arrays ``a`` and ``b`` specifying the corners of a hypercube
    to the required tolerance.

    Convergence is not guaranteed for all integrals.

    Parameters
    ----------
    f : callable
        Function to integrate. ``f`` must have the signature::
            f(x : ndarray, \*args, \*\*kwargs) -> ndarray

        If ``f`` accepts arrays ``x`` of shape ``(input_dim_1, ..., input_dim_n,
        num_eval_points)`` and outputs arrays of shape ``(output_dim_1, ...,
        output_dim_n, num_eval_points)``, then ``cub`` will return arrays of shape
        ``(output_dim_1, ..., output_dim_n)``.
    a, b : array_like or float
        Lower and upper limits of integration as rank-1 arrays specifying the left and
        right endpoints of the intervals being integrated over. If a float is passed,
        these will be automatically converted to singleton arrays. Infinite limits are
        currently not supported.
    rule : Cubature or str, optional
        Cubature rule to use when estimating the integral. If passing a string, the
        options are "gk21", "gk15", or "trapezoid". Default is "gk21".
    rtol : float, optional
        Relative tolerance. Default is 1e-05.
    atol : float, optional
        Absolute tolerance. Default is 1e-08.
    max_subdivisions : int, optional
        Upper bound on the number of subdivisions to perform to improve the estimate
        over a subregion. Default is 10,000.
    args : tuple, optional
        Additional positional args passed to ``f``, if any.
    kwargs : tuple, optional
        Additional keyword args passed to ``f``, if any.

    Returns
    -------
    res : CubatureResult
        Result of estimation. See CubatureResult.

    Examples
    --------
    A simple 1D integral with vector output:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import GaussKronrodQuad
    >>> def f(x, n):
    ...    # Make sure x and n are broadcastable
    ...    return x.reshape(1, -1)**n
    >>> res = cub(
    ...     f,
    ...     [0],
    ...     [1],
    ...     GaussKronrodQuad(21), # 1D rule
    ...     args=(
    ...         # Since f accepts arrays of shape (ndim, num_eval_points) we need to
    ...         # make sure n is the right shape
    ...         np.arange(10).reshape(-1, 1),
    ...     )
    ... )
    >>> res.estimate
     array([1.        , 0.5       , 0.33333333, 0.25      , 0.2       ,
            0.16666667, 0.14285714, 0.125     , 0.11111111, 0.1       ])


    A 7D integral with arbitrary-shaped array output:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import GenzMalikCub
    >>> def f(x, r, alphas):
    ...     # f(x) = cos(2pi*r + alpha @ x)
    ...     ndim = x.shape[0]
    ...     num_eval_points = x.shape[-1]
    ...     r_reshaped = np.expand_dims(r, -1)
    ...     alphas_reshaped = np.expand_dims(alphas, -1)
    ...     x_reshaped = x.reshape(
    ...         ndim,
    ...         *([1]*(len(alphas.shape) - 1)),
    ...         num_eval_points
    ...     )
    ...     return np.cos(
    ...         2*np.pi*r_reshaped + np.sum(alphas_reshaped * x_reshaped, axis=0)
    ...     )
    >>> np.random.seed(1)
    >>> res = cub(
    ...     f=f,
    ...     a=np.array([0, 0, 0, 0, 0, 0, 0]),
    ...     b=np.array([1, 1, 1, 1, 1, 1, 1]),
    ...     rule=GenzMalikCub(7),
    ...     kwargs={
    ...         "r": np.random.rand(2, 3),
    ...         "alphas": np.random.rand(7, 2, 3),
    ...     }
    ... )
    >>> res.estimate
     array([[ 0.69578554, -0.87958878, -0.90278537],
            [ 0.53436481, -0.34633208, -0.16061907]])

    Calculating an integral using a custom cubature rule, in this case 3D Genz-Malik
    cubature for the estimate and then the difference between this and 3D 21-node
    Gauss-Kronrod for the error estimate.

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub, FixedProductErrorFromDifferenceCub
    >>> from scipy.integrate._rules import Cub, GenzMalikCub, GaussKronrodQuad
    >>> def f(x, r, alphas):
    ...     # f(x) = cos(2pi*r + alpha @ x)
    ...     ndim = x.shape[0]
    ...     num_eval_points = x.shape[-1]
    ...     r_reshaped = np.expand_dims(r, -1)
    ...     alphas_reshaped = np.expand_dims(alphas, -1)
    ...     x_reshaped = x.reshape(
    ...         ndim,
    ...         *([1]*(len(alphas.shape) - 1)),
    ...         num_eval_points
    ...     )
    ...     return np.cos(
    ...         2*np.pi*r_reshaped + np.sum(alphas_reshaped * x_reshaped, axis=0)
    ...     )
    >>> genz = GenzMalikCub(3)
    >>> kronrod = FixedProductErrorFromDifferenceCub([GaussKronrodQuad(21)] * 3)
    >>> class CustomRule(Cub):
    ...     def estimate(self, f, a, b, args=(), kwargs=None):
    ...         if kwargs is None: kwargs = dict()
    ...         return genz.estimate(f, a, b, args, kwargs)
    ...     def estimate_error(self, f, a, b, args=(), kwargs=None):
    ...         if kwargs is None: kwargs = dict()
    ...         return np.abs(
    ...             genz.estimate(f, a, b, args, kwargs)
    ...             - kronrod.estimate(f, a, b, args, kwargs)
    ...         )
    >>> np.random.seed(1)
    >>> r, alphas = np.random.rand(2, 3), np.random.rand(3, 2, 3)
    >>> cub(
    ...     f=f,
    ...     a=np.array([0, 0, 0]),
    ...     b=np.array([1, 1, 1]),
    ...     rule=CustomRule(),
    ...     kwargs={
    ...         "r": r,
    ...         "alphas": alphas,
    ...     }
    ... ).estimate
     array([[ 0.6369219 , -0.90656298, -0.82183961],
            [ 0.48873062, -0.87426057,  0.50205197]])
    """

    if max_subdivisions is None:
        max_subdivisions = np.inf

    if kwargs is None:
        kwargs = dict()

    # Convert a and b to arrays of at least 1D
    a = np.atleast_1d(a)
    b = np.atleast_1d(b)

    # If the rule is a string, convert to a corresponding product rule
    if isinstance(rule, str):
        quadratues = {
            "gk21": GaussKronrodQuad(21),
            "gk15": GaussKronrodQuad(15),
            "trapezoid": ErrorFromDifference(NewtonCotesQuad(5), NewtonCotesQuad(3)),
        }

        base_quadrature = quadratues.get(rule)

        if base_quadrature is None:
            raise ValueError(f"unknown rule {rule}")

        rule = FixedProductErrorFromDifferenceCub([base_quadrature] * len(a))

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("a and b should be 1D arrays")

    est = rule.estimate(f, a, b, args, kwargs)

    try:
        err = rule.estimate_error(f, a, b, args, kwargs)
    except NotImplementedError:
        raise ValueError("attempting cubature with a rule that doesn't implement error \
estimation.")

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
