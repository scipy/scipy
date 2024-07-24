import math
import heapq
import itertools

from dataclasses import dataclass
from functools import cached_property

import numpy as np


def cub(f, a, b, rule, rtol=1e-05, atol=1e-08, max_subdivisions=10000, args=()):
    """
    Adaptive cubature of vector- or tensor-valued function.

    Given an arbitrary cubature rule with error estimation, this function returns the
    estimate of the integral over the region described by corners ``a`` and ``b`` to
    the required tolerance.

    Parameters
    ----------
    f : callable
        Function to integrate. ``f`` must have the signature::
            f(x : ndarray, *args) -> ndarray
        If ``f(x, *args)`` accepts arrays ``x`` of shape ``(input_dim_1, ...,
        input_dim_n, num_eval_points)`` and outputs arrays of shape ``(output_dim_1,
        ..., output_dim_n, num_eval_points)``, then ``cub`` will return arrays of shape
        ``(output_dim_1, ..., output_dim_n)``.
    a, b : ndarray
        Lower and upper limits of integration. If the integral is being performed over
        intervals x_1 in [a_1, b_1], x_2 in [a_2, b_2], ..., x_n in [a_n, b_n], then
        a and b will be [a_1, ..., a_n] and [b_1, ..., b_n] correspondingly.
    rule : Cubature
        Cubature rule to use when estimating the integral. The cubature rule specified
        must implement ``error_estimate``. See Examples.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.
    max_subdivisions : int
        Maximum number of times to subdivide the region of integration in order to
        improve the estimate over a region with high error.
    args : tuple, optional
        Additional positional args passed to ``f``. See Examples.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x, alphas):
    ...     # Example function to integrate
    ...     # f(x) = exp(sum(alpha * x))
    ...
    ...     # Support arbitrary shaped x and alphas
    ...     ndim = x.shape[0]
    ...     num_eval_points = x.shape[-1]
    ...
    ...     alphas_reshaped = np.expand_dims(alphas, -1)
    ...     x_reshaped = x.reshape(
    ...         ndim, *([1]*(len(alphas.shape) - 1)), num_eval_points
    ...     )
    ...
    ...     return np.exp(np.sum(alphas_reshaped * x_reshaped, axis=0))
    >>> rule_1d = GaussKronrod(21) # Use Gauss-Kronrod, which has error estimate
    >>> alphas_1d = np.array([0.5])
    >>> res = cub(
    ...     f,
    ...     np.array([0]),
    ...     np.array([1]),
    ...     rule_1d,
    ...     args=(alphas_1d,)
    ... ); res.estimate
     np.float64(1.2974425414002562)
    >>> # Now calculate 3D integral, constructing 3D rule from product of 1D rules
    >>> product_rule = Product([rule_1d, rule_1d, rule_1d])
    >>> # New values of alpha, now f calcualtes exp(0.2*x_1 + 0.5*x_2 + 0.8*x_3)
    >>> alphas_3d = np.array([0.2, 0.5, 0.8])
    >>> res = cub(
    ...     f,
    ...     np.array([0, 0, 0]),
    ...     np.array([1, 1, 1]),
    ...     product_rule,
    ...     args=(alphas_3d,)
    ... ); res.estimate
     np.float64(2.2002853017758057)
    >>> # Evaluate for many values of alpha simultaneously
    >>> alphas_many = np.array([[0.1, 0.2, 2], [0.25, 0.5, 5], [0.4, 0.8, 8]])
    >>> res = cub(
    ...     f,
    ...     np.array([0, 0, 0]),
    ...     np.array([1, 1, 1]),
    ...     rule_3d,
    ...     args=(alphas_many,)
    ... ); res.estimate
     array([1.46914007e+00, 2.20028530e+00, 3.50827109e+04])
    """

    if max_subdivisions is None:
        max_subdivisions = np.inf

    est = rule.estimate(f, a, b, args)
    err = rule.error_estimate(f, a, b, args)

    if err is None:
        raise Exception("attempting cubature with a rule that doesn't implement error \
estimation.")

    regions = [CubatureRegion(est, err, a, b)]
    subdivisions = 0
    success = False

    while np.any(err > atol + rtol * np.abs(est)):
        # region_k is the region with highest estimated error
        region_k = heapq.heappop(regions)

        est_k = region_k.estimate
        err_k = region_k.error

        a_k, b_k = region_k.a, region_k.b
        m_k = (a_k + b_k) / 2

        # Find all 2^ndim subregions formed by splitting region_k along each axis
        # e.g. for 1D quadrature this splits an interval in two, for 3D cubature this
        # splits the cube under consideration into 8 subcubes.
        subregion_coordinates = zip(
            _cartesian_product(np.array([a_k, m_k]).T).T,
            _cartesian_product(np.array([m_k, b_k]).T).T
        )

        # Subtract the estimate of the integral and its error from the current global
        # estimates, since these will be refined in the loop over all subregions.
        est -= est_k
        err -= err_k

        # For each of the new subregions, calculate an estimate for the integral and
        # the error there, and push these regions onto the heap for potential further
        # subdividing.
        for a_k_sub, b_k_sub in subregion_coordinates:
            est_sub = rule.estimate(f, a_k_sub, b_k_sub, args)
            err_sub = rule.error_estimate(f, a_k_sub, b_k_sub, args)

            est += est_sub
            err += err_sub

            new_region = CubatureRegion(est_sub, err_sub, a_k_sub, b_k_sub)

            heapq.heappush(regions, new_region)

        subdivisions += 1

        if subdivisions >= max_subdivisions:
            break
    else:
        success = True

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


class Cubature:
    """
    A generic interface for numerical cubature algorithms.

    Finds an estimate for the integral of ``f`` over a (hyper)rectangular region
    described by two points ``a`` and ``b``, and optionally finds an estimate for the
    error of this approximation.

    If the cubature rule doesn't support error estimation, error_estimate will return
    None.
    """

    def estimate(self, f, a, b, args=()):
        """
        Calculate estimate of integral of ``f`` in region described by corners ``a``
        and ``b``.

        Parameters
        ----------
        f : callable
            Function ``f(x, *args)`` to integrate. ``f`` should accept arrays ``x`` of
            shape ``(ndim, eval_points)`` and return arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``.

        a, b : ndarray
            Lower and upper limits of integration.

        args : tuple, optional
            Extra arguments to pass to ``f``, if any.

        Returns
        -------
        est : ndarray
            Result of estimation. If ``f`` returns arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``, then ``err_est`` will be of shape
            ``(output_dim_1, ..., output_dim_n)``.
        """
        raise NotImplementedError

    def error_estimate(self, f, a, b, args=()):
        """
        Calculate the error estimate of this cubature rule for the integral of ``f`` in
        the region described by corners ``a`` and ``b``.

        If the cubature rule doesn't support error estimation, this will be ``None``.

        Parameters
        ----------
        f : callable
            Function ``f(x, *args)`` to integrate. ``f`` should accept arrays ``x`` of
            shape ``(ndim, eval_points)`` and return arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``.

        a, b : ndarray
            Lower and upper limits of integration.

        args : tuple, optional
            Extra arguments to pass to ``f``, if any.

        Returns
        -------
        err_est : ndarray
            Estimate of the error. If the cubature rule doesn't support error
            estimation, then this will be ``None``. If error estimation is supported and
            ``f`` returns arrays of shape ``(output_dim_1, ..., output_dim_n,
            eval_points)``, then ``err_est`` will be of shape ``(output_dim_1, ...,
            output_dim_n)``.
        """
        return None


class FixedCubature(Cubature):
    """
    A numerical cubature rule implemented as the weighted sum of function evaluations.

    See Also
    --------
    NewtonCotes, DualEstimateCubature
    """

    @property
    def nodes(self):
        """
        Nodes at which to evaluate the function being integrated. For an n-dim function,
        these will be of shape ``(ndim, num_nodes)``.

        Nodes should be for integrals over the region ``[-1, 1]^n``.
        """

        raise NotImplementedError

    @property
    def weights(self):
        """
        Weights to multiply function evaluations at nodes in the estimate. For an n-dim
        function, these will be of shape ``(num_nodes,)``.
        """

        raise NotImplementedError

    def estimate(self, f, a, b, args=()):
        """
        Calculate estimate of integral of ``f`` in region described by corners ``a``
        and ``b`` as ``sum(weights * f(nodes))``. Nodes and weights will automatically
        be adjusted from calculating integrals over [-1, 1]^n to [a, b].

        Parameters
        ----------
        f : callable
            Function ``f(x, *args)`` to integrate. ``f`` should accept arrays ``x`` of
            shape ``(ndim, eval_points)`` and return arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``.

        a, b : ndarray
            Lower and upper limits of integration, specified as arrays of shape
            ``(ndim,)``. Note that ``f`` accepts arrays of shape ``(ndim, eval_points)``
            so in order to calculate ``f(a)``, it is necessary to add an extra axis
            ``f(a[:, np.newaxis])``.

        args : tuple, optional
            Extra arguments to pass to ``f``, if any.

        Returns
        -------
        est : ndarray
            Result of estimation. If ``f`` returns arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``, then ``err_est`` will be of shape
            ``(output_dim_1, ..., output_dim_n)``.
        """

        # Since f accepts arrays of shape (ndim, eval_points), it is necessary to
        # add an extra axis to a and b so that ``f`` can be evaluated there.
        a = a[:, np.newaxis]
        b = b[:, np.newaxis]
        lengths = b - a

        # The underlying cubature rule is for the hypercube [-1, 1]^n.
        #
        # To handle arbitrary regions of integration, it's necessary to apply a linear
        # change of coordinates to map each interval [a[i], b[i]] to [-1, 1].
        nodes = (self.nodes + 1) * (lengths / 2) + a

        # Also need to multiply the weights by a scale factor equal to the determinant
        # of the Jacobian for this coordinate change.
        weight_scale_factor = math.prod(lengths / 2)
        weights = self.weights * weight_scale_factor

        # f(nodes) will have shape (output_dim_1, ..., output_dim_n, num_nodes)
        # Summing along the last axis means estimate will shape (output_dim_1, ...,
        # output_dim_n)
        est = np.sum(
            weights * f(nodes, *args),
            axis=-1
        )

        return est

    def error_estimate(self, f, a, b, args=()):
        # A fixed cubature rule has no built-in error estimation in general.
        return None


class DualEstimateCubature(Cubature):
    """
    A cubature rule with a higher-order and lower-order estimate, and where the
    difference between the two estimates is used as an estimate for the error.

    Attributes
    ----------
    higher : Cubature
        Higher-order cubature to use as the estimate.

    lower : Cubature
        Lower-order cubature. The difference between ``higher.estimate`` and
        ``lower.estimate`` is used as an estimate for the error.

    See Also
    --------
    GaussKronrod

    Examples
    --------
    Newton-Cotes quadrature doesn't come with any error estimation built in. It is
    possible to give it an error estimator by creating a DualEstimateCubature which
    takes the difference between 10-node Newton-Cotes and 8-node Newton-Cotes:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     return np.sin(np.sqrt(x))
    >>> rule = NewtonCotes(10)
    >>> a, b = np.array([0]), np.array([np.pi * 2])
    >>> rule.error_estimate(f, a, b) is None
     True
    >>> rule_with_err_est = DualEstimateCubature(NewtonCotes(10), NewtonCotes(8))
    >>> rule_with_err_est.error_estimate(f, a, b)
     array([6.21045267])
    """

    def __init__(self, higher, lower):
        self.higher = higher
        self.lower = lower

    def estimate(self, f, a, b, args=()):
        return self.higher.estimate(f, a, b, args)

    def error_estimate(self, f, a, b, args=()):
        return np.abs(
            self.higher.estimate(f, a, b, args) - self.lower.estimate(f, a, b, args)
        )


class GaussKronrod(DualEstimateCubature):
    """
    Gauss-Kronrod cubature. Gauss-Kronrod rules consist of two cubature rules, one
    higher-order and one lower-order. The higher-order rule is used as the estimate of
    the integral and the difference between them is used as an estimate for the error.

    Gauss-Kronrod is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the Product of multiple Gauss-Kronrod rules. See Examples.

    For n-node Gauss-Kronrod, the lower-order rule has ``n//2`` nodes, which are the
    ordinary Gauss-Legendre nodes with corresponding weights. The higher-order rule has
    ``n`` nodes, ``n//2`` of which are the same as the lower-order rule and the
    remaining nodes are the Kronrod extension of those nodes.

    Parameters
    ----------
    npoints : int
        Number of nodes for the higher-order rule.

    Attributes
    ----------
    higher : Cubature
        Higher-order rule.

    lower : Cubature
        Lower-order rule.

    Examples
    --------
    Evaluate a 1D integral. Note in this example that ``f`` returns an array, so the
    estimates will also be arrays, despite the fact that this is a 1D problem.

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     return np.cos(x)
    >>> rule = GaussKronrod(21) # Use 21-point GaussKronrod
    >>> a, b = np.array([0]), np.array([1])
    >>> rule.estimate(f, a, b) # True value sin(1), approximately 0.84147
     array([0.84147098])
    >>> rule.error_estimate(f, a, b)
     array([1.11022302e-16])

    Evaluate a 2D integral. Note that in this example ``f`` returns a float, so the
    estimates will also be floats.

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule = Product([GaussKronrod(15), GaussKronrod(15)]) # Use 15-point GaussKronrod
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.error_estimate(f, a, b)
     np.float64(2.220446049250313e-16)
    """

    def __init__(self, npoints):
        if npoints != 15 and npoints != 21:
            raise Exception("Gauss-Kronrod quadrature is currently only supported for \
15 or 21 nodes")

        self.higher = self._Higher(npoints)
        self.lower = self._Lower(npoints//2)

    class _Higher(FixedCubature):
        def __init__(self, npoints):
            self.npoints = npoints

        @cached_property
        def nodes(self):
            if self.npoints == 21:
                return np.array([[
                    0.995657163025808080735527280689003,
                    0.973906528517171720077964012084452,
                    0.930157491355708226001207180059508,
                    0.865063366688984510732096688423493,
                    0.780817726586416897063717578345042,
                    0.679409568299024406234327365114874,
                    0.562757134668604683339000099272694,
                    0.433395394129247190799265943165784,
                    0.294392862701460198131126603103866,
                    0.148874338981631210884826001129720,
                    0,
                    -0.148874338981631210884826001129720,
                    -0.294392862701460198131126603103866,
                    -0.433395394129247190799265943165784,
                    -0.562757134668604683339000099272694,
                    -0.679409568299024406234327365114874,
                    -0.780817726586416897063717578345042,
                    -0.865063366688984510732096688423493,
                    -0.930157491355708226001207180059508,
                    -0.973906528517171720077964012084452,
                    -0.995657163025808080735527280689003
                ]])
            elif self.npoints == 15:
                return np.array([[
                    0.991455371120812639206854697526329,
                    0.949107912342758524526189684047851,
                    0.864864423359769072789712788640926,
                    0.741531185599394439863864773280788,
                    0.586087235467691130294144838258730,
                    0.405845151377397166906606412076961,
                    0.207784955007898467600689403773245,
                    0.000000000000000000000000000000000,
                    -0.207784955007898467600689403773245,
                    -0.405845151377397166906606412076961,
                    -0.586087235467691130294144838258730,
                    -0.741531185599394439863864773280788,
                    -0.864864423359769072789712788640926,
                    -0.949107912342758524526189684047851,
                    -0.991455371120812639206854697526329
                ]])

        @cached_property
        def weights(self):
            if self.npoints == 21:
                return np.array([
                    0.011694638867371874278064396062192,
                    0.032558162307964727478818972459390,
                    0.054755896574351996031381300244580,
                    0.075039674810919952767043140916190,
                    0.093125454583697605535065465083366,
                    0.109387158802297641899210590325805,
                    0.123491976262065851077958109831074,
                    0.134709217311473325928054001771707,
                    0.142775938577060080797094273138717,
                    0.147739104901338491374841515972068,
                    0.149445554002916905664936468389821,
                    0.147739104901338491374841515972068,
                    0.142775938577060080797094273138717,
                    0.134709217311473325928054001771707,
                    0.123491976262065851077958109831074,
                    0.109387158802297641899210590325805,
                    0.093125454583697605535065465083366,
                    0.075039674810919952767043140916190,
                    0.054755896574351996031381300244580,
                    0.032558162307964727478818972459390,
                    0.011694638867371874278064396062192,
                ])
            elif self.npoints == 15:
                return np.array([
                    0.022935322010529224963732008058970,
                    0.063092092629978553290700663189204,
                    0.104790010322250183839876322541518,
                    0.140653259715525918745189590510238,
                    0.169004726639267902826583426598550,
                    0.190350578064785409913256402421014,
                    0.204432940075298892414161999234649,
                    0.209482141084727828012999174891714,
                    0.204432940075298892414161999234649,
                    0.190350578064785409913256402421014,
                    0.169004726639267902826583426598550,
                    0.140653259715525918745189590510238,
                    0.104790010322250183839876322541518,
                    0.063092092629978553290700663189204,
                    0.022935322010529224963732008058970
                ])

    class _Lower(FixedCubature):
        def __init__(self, npoints):
            self.npoints = npoints

        @cached_property
        def nodes(self):
            if self.npoints == 10:
                return np.array([[
                    0.973906528517171720077964012084452,
                    0.865063366688984510732096688423493,
                    0.679409568299024406234327365114874,
                    0.433395394129247190799265943165784,
                    0.148874338981631210884826001129720,
                    -0.148874338981631210884826001129720,
                    -0.433395394129247190799265943165784,
                    -0.679409568299024406234327365114874,
                    -0.865063366688984510732096688423493,
                    -0.973906528517171720077964012084452,
                ]])
            elif self.npoints == 7:
                return np.array([[
                    0.949107912342758524526189684047851,
                    0.741531185599394439863864773280788,
                    0.405845151377397166906606412076961,
                    0.000000000000000000000000000000000,
                    -0.405845151377397166906606412076961,
                    -0.741531185599394439863864773280788,
                    -0.949107912342758524526189684047851,
                ]])

        @cached_property
        def weights(self):
            if self.npoints == 10:
                return np.array([
                    0.066671344308688137593568809893332,
                    0.149451349150580593145776339657697,
                    0.219086362515982043995534934228163,
                    0.269266719309996355091226921569469,
                    0.295524224714752870173892994651338,
                    0.295524224714752870173892994651338,
                    0.269266719309996355091226921569469,
                    0.219086362515982043995534934228163,
                    0.149451349150580593145776339657697,
                    0.066671344308688137593568809893332,
                ])
            elif self.npoints == 7:
                return np.array([
                    0.129484966168869693270611432679082,
                    0.279705391489276667901467771423780,
                    0.381830050505118944950369775488975,
                    0.417959183673469387755102040816327,
                    0.381830050505118944950369775488975,
                    0.279705391489276667901467771423780,
                    0.129484966168869693270611432679082
                ])


class NewtonCotes(FixedCubature):
    """
    Newton-Cotes cubature. Newton-Cotes rules consist of function evaluations at equally
    spaced nodes.

    Newton-Cotes has no error estimator. It is possible to give it an error estimator
    by using DualEstimateCubature to estimate the error as the difference between two
    Newton-Cotes rules. See Examples.

    Newton-Cotes is a 1D rule. To use it for multidimensional integrals, it will be
    necessary to take the product of multiple Newton-Cotes rules. See Examples.

    Parameters
    ----------
    npoints : int
        Number of equally spaced evaluation points.

    open : bool, default=False
        Whether this should be an open rule. Open rules do not include the endpoints
        as nodes.

    Attributes
    ----------
    npoints : int
        Number of equally spaced evaluation points.

    open : bool, default=False
        Whether this is an open rule. Open rules do not include the endpoints as nodes.

    Examples
    --------
    Evaluating a simple integral, first without error estimation and then with error
    estimation:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     return np.sin(np.sqrt(x))
    >>> rule = NewtonCotes(10)
    >>> a, b = np.array([0]), np.array([np.pi * 2])
    >>> rule.estimate(f, a, b) # Very inaccuracte
     array([0.60003617])
    >>> rule.error_estimate(f, a, b) is None
     True
    >>> rule_with_err_est = DualEstimateCubature(NewtonCotes(10), NewtonCotes(8))
    >>> rule_with_err_est.error_estimate(f, a, b) # Error is high
     array([6.21045267])

    Evaluating a 2D integral, using the product of NewtonCotes with an error estimator:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule_with_err_est = DualEstimateCubature(NewtonCotes(10), NewtonCotes(8))
    >>> rule_2d = Product([rule_with_err_est, rule_with_err_est])
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule_2d.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.6829419696151342)
    >>> rule_2d.error_estimate(f, a, b)
     np.float64(6.823492881835591e-10)
    """

    def __init__(self, npoints, open=False):
        if npoints < 2:
            raise Exception(
                "At least 2 points required for Newton-Cotes cubature"
            )

        self.npoints = npoints
        self.open = open

    @cached_property
    def nodes(self):
        if self.open:
            h = 2/self.npoints
            return np.linspace(-1 + h, 1 - h, num=self.npoints).reshape(1, -1)
        else:
            return np.linspace(-1, 1, num=self.npoints).reshape(1, -1)

    @cached_property
    def weights(self):
        return _newton_cotes_weights(self.nodes.reshape(-1))


class Product(DualEstimateCubature):
    """
    Find the n-dimensional cubature rule constructed from cubature rules of lower
    dimension.

    Given a list of base dual-estimate cubature rules of dimension d_1, ..., d_n, this
    will find the (d_1 + ... d_n)-dimensional cubature rule formed from taking the
    Cartesian product of their weights.

    Parameters
    ----------
    base_rules : list
        List of base dual estimate cubatures used to find the product rule.

    Attributes
    ----------
    base_rules : list of DualEstimateCubature
        List of base dual estimate cubatures used to find the product rule.

    higher : FixedCubature
        Higher-order rule given by the product of the higher order rules given in the
        base rules.

    lower : FixedCubature
        Lower-order rule given by the product of the lower order rules given in the base
        rules.

    Example
    -------

    Evaluate a 2D integral by taking the product of two 1D rules:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule = Product([GaussKronrod(15), GaussKronrod(15)]) # Use 15-point GaussKronrod
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.error_estimate(f, a, b)
     np.float64(2.220446049250313e-16)
    """

    def __init__(self, base_rules):
        self.base_rules = base_rules
        self.higher = self._Higher(base_rules)
        self.lower = self._Lower(base_rules)

    class _Higher(FixedCubature):
        def __init__(self, base_rules):
            self.base_rules = base_rules

        @cached_property
        def nodes(self):
            return _cartesian_product([rule.higher.nodes for rule in self.base_rules])

        @cached_property
        def weights(self):
            return np.prod(
                _cartesian_product([rule.higher.weights for rule in self.base_rules]),
                axis=0
            )

    class _Lower(FixedCubature):
        def __init__(self, base_rules):
            self.base_rules = base_rules

        @cached_property
        def nodes(self):
            return _cartesian_product([rule.lower.nodes for rule in self.base_rules])

        @cached_property
        def weights(self):
            return np.prod(
                _cartesian_product([rule.lower.weights for rule in self.base_rules]),
                axis=0
            )


class GenzMalik(DualEstimateCubature):
    """
    Genz-Malik cubature. Genz-Malik cubature is a true cubature rule in that it is not
    constructed as the product of 1D rules.

    Genz-Malik is only defined for integrals of dimension >= 2.

    Parameters
    ----------
    ndim : int
        The spacial dimension of the integrand.

    Attributes
    ----------
    higher : Cubature
        Higher-order rule.

    lower : Cubature
        Lower-order rule.

    Examples
    --------
    Evaluate a 3D integral:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import *
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2) + cos(x_3)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule = GenzMalik(3) # Use 3D Genz-Malik
    >>> a, b = np.array([0, 0, 0]), np.array([1, 1, 1])
    >>> rule.estimate(f, a, b) # True value 3*sin(1), approximately 2.5244
     np.float64(2.5244129547230862)
    >>> rule.error_estimate(f, a, b)
     np.float64(1.378269656626685e-06)
    """

    def __init__(self, ndim):
        if ndim < 2:
            raise Exception("Genz-Malik cubature is only defined for ndim >= 2")

        self.higher = self._Higher(ndim)
        self.lower = self._Lower(ndim)

    class _Higher(FixedCubature):
        def __init__(self, ndim):
            self.ndim = ndim

        @cached_property
        def nodes(self):
            l_2 = math.sqrt(9/70)
            l_3 = math.sqrt(9/10)
            l_4 = math.sqrt(9/10)
            l_5 = math.sqrt(9/19)

            its = itertools.chain(
                [(0,) * self.ndim],
                _distinct_permutations((l_2,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((-l_2,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((l_3,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((-l_3,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((l_4, l_4) + (0,) * (self.ndim - 2)),
                _distinct_permutations((l_4, -l_4) + (0,) * (self.ndim - 2)),
                _distinct_permutations((-l_4, -l_4) + (0,) * (self.ndim - 2)),
                itertools.product((l_5, -l_5), repeat=self.ndim)
            )

            out_size = 1 + 2 * (self.ndim + 1) * self.ndim + 2**self.ndim

            out = np.fromiter(
                itertools.chain.from_iterable(zip(*its)),
                dtype=float,
                count=self.ndim * out_size
            )

            out.shape = (self.ndim, out_size)
            return out

        @cached_property
        def weights(self):
            w_1 = (2**self.ndim) * (12824 - 9120 * self.ndim + 400 * self.ndim**2) \
                / 19683
            w_2 = (2**self.ndim) * 980/6561
            w_3 = (2**self.ndim) * (1820 - 400 * self.ndim) / 19683
            w_4 = (2**self.ndim) * (200 / 19683)
            w_5 = 6859 / 19683

            return np.repeat(
                [w_1, w_2, w_3, w_4, w_5],
                [
                    1,
                    2 * self.ndim,
                    2*self.ndim,
                    2*(self.ndim - 1)*self.ndim,
                    2**self.ndim,
                ]
            )

    class _Lower(FixedCubature):
        def __init__(self, ndim):
            self.ndim = ndim

        @cached_property
        def nodes(self):
            l_2 = math.sqrt(9/70)
            l_3 = math.sqrt(9/10)
            l_4 = math.sqrt(9/10)

            its = itertools.chain(
                [(0,) * self.ndim],
                _distinct_permutations((l_2,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((-l_2,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((l_3,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((-l_3,) + (0,) * (self.ndim - 1)),
                _distinct_permutations((l_4, l_4) + (0,) * (self.ndim - 2)),
                _distinct_permutations((l_4, -l_4) + (0,) * (self.ndim - 2)),
                _distinct_permutations((-l_4, -l_4) + (0,) * (self.ndim - 2)),
            )

            out_size = 1 + 2 * (self.ndim + 1) * self.ndim

            out = np.fromiter(
                itertools.chain.from_iterable(zip(*its)),
                dtype=float,
                count=self.ndim * out_size
            )

            out.shape = (self.ndim, out_size)
            return out

        @cached_property
        def weights(self):
            w_1 = (2**self.ndim) * (729 - 950*self.ndim + 50*self.ndim**2) / 729
            w_2 = (2**self.ndim) * (245 / 486)
            w_3 = (2**self.ndim) * (265 - 100*self.ndim) / 1458
            w_4 = (2**self.ndim) * (25 / 729)

            return np.repeat(
                [w_1, w_2, w_3, w_4],
                [1, 2 * self.ndim, 2*self.ndim, 2*(self.ndim - 1)*self.ndim]
            )


def _cartesian_product(points):
    """
    Takes a list of arrays such as `[ [[x_1, x_2]], [[y_1, y_2]] ]` and
    returns all combinations of these points, such as
    `[[x_1, x_1, x_2, x_2], [y_1, y_2, y_1, y_2]]`

    Note that this is assuming that the spatial dimension is the first axis, as opposed
    to the last axis.
    """
    out = np.stack(np.meshgrid(*points, indexing='ij'), axis=0)
    out = out.reshape(len(points), -1)
    return out


def _max_norm(x):
    return np.max(np.abs(x))


def _distinct_permutations(iterable, r=None):
    """
    Find the number of distinct permutations of `r` elements of `iterable`.
    """

    # Algorithm: https://w.wiki/Qai
    def _full(A):
        while True:
            # Yield the permutation we have
            yield tuple(A)

            # Find the largest index i such that A[i] < A[i + 1]
            for i in range(size - 2, -1, -1):
                if A[i] < A[i + 1]:
                    break

            #  If no such index exists, this permutation is the last one
            else:
                return

            # Find the largest index j greater than j such that A[i] < A[j]
            for j in range(size - 1, i, -1):
                if A[i] < A[j]:
                    break

            # Swap the value of A[i] with that of A[j], then reverse the
            # sequence from A[i + 1] to form the new permutation
            A[i], A[j] = A[j], A[i]
            A[i+1:] = A[:i-size:-1]  # A[i + 1:][::-1]

    # Algorithm: modified from the above
    def _partial(A, r):
        # Split A into the first r items and the last r items
        head, tail = A[:r], A[r:]
        right_head_indexes = range(r - 1, -1, -1)
        left_tail_indexes = range(len(tail))

        while True:
            # Yield the permutation we have
            yield tuple(head)

            # Starting from the right, find the first index of the head with
            # value smaller than the maximum value of the tail - call it i.
            pivot = tail[-1]

            for i in right_head_indexes:
                if head[i] < pivot:
                    break

                pivot = head[i]

            else:
                return

            # Starting from the left, find the first value of the tail
            # with a value greater than head[i] and swap.
            for j in left_tail_indexes:
                if tail[j] > head[i]:
                    head[i], tail[j] = tail[j], head[i]
                    break

            # If we didn't find one, start from the right and find the first
            # index of the head with a value greater than head[i] and swap.
            else:
                for j in right_head_indexes:
                    if head[j] > head[i]:
                        head[i], head[j] = head[j], head[i]
                        break

            # Reverse head[i + 1:] and swap it with tail[:r - (i + 1)]
            tail += head[:i-r:-1]  # head[i + 1:][::-1]
            i += 1

            head[i:], tail[:] = tail[:r-i], tail[r-i:]

    items = sorted(iterable)
    size = len(items)

    if r is None:
        r = size

    if 0 < r <= size:
        return _full(items) if (r == size) else _partial(items, r)

    return iter(() if r else ((),))


def _newton_cotes_weights(points):
    order = len(points) - 1

    a = np.vander(points, increasing=True)
    a = np.transpose(a)

    i = np.arange(order + 1)
    b = (1 - np.power(-1, i + 1)) / (2 * (i + 1))

    # TODO: figure out why this 2x is necessary
    return 2*np.linalg.solve(a, b)
