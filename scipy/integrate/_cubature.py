import math
import heapq
import functools

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


def cub(f, a, b, rule, rtol=1e-05, atol=1e-08, limit=10000, args=()):
    if limit is None:
        limit = np.inf

    est = rule.estimate(f, a, b, args)
    err = rule.error_estimate(f, a, b, args)

    regions = [CubatureRegion(est, err, a, b)]
    subdivisions = 0
    success = False

    while np.any(err > atol + rtol * np.abs(est)):
        region = heapq.heappop(regions)

        est_k = region.estimate
        err_k = region.error

        a_k, b_k = region.a, region.b
        m_k = (a_k + b_k) / 2

        subregion_coordinates = zip(
            _cartesian_product(np.array([a_k, m_k]).T).T,
            _cartesian_product(np.array([m_k, b_k]).T).T
        )

        est -= est_k
        err -= err_k

        for a_k_sub, b_k_sub in subregion_coordinates:
            est_sub = rule.estimate(f, a_k_sub, b_k_sub, args)
            err_sub = rule.error_estimate(f, a_k_sub, b_k_sub, args)

            est += est_sub
            err += err_sub

            new_region = CubatureRegion(est_sub, err_sub, a_k_sub, b_k_sub)

            heapq.heappush(regions, new_region)

        subdivisions += 1

        if subdivisions >= limit:
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
        # Compare negative error so that regions with higher error are placed closer to
        # top of the heap.
        return -_max_norm(self.error) < -_max_norm(other.error)


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


class CubatureRule(ABC):
    def __init__(self, nodes, weights):
        self.nodes = nodes
        self.weights = weights

    @property
    def ndim(self):
        return self.nodes.shape[0]

    @abstractmethod
    def estimate(self, f, a, b, args):
        """
        Calculate estimate of integral of f in region described by corners a and b.

        f should accept arrays of shape (input_dim, num_eval_points) and return arrays
        of shape (output_dim_1, ..., output_dim_n, num_eval_points).
        """
        pass

    @abstractmethod
    def error_estimate(self, f, a, b, args):
        """
        Calculate error estimate for the integral of f in region described by corners a
        and b.

        f should accept arrays of shape (input_dim, num_eval_points) and return arrays
        of shape (output_dim_1, ..., output_dim_n, num_eval_points).
        """
        pass


class NestedCubatureRule(CubatureRule):
    """
    A cubature rule with a higher-order and lower-order estimate, and where the
    difference between the two is used to estimate the error.
    """
    def __init__(self, higher_nodes, higher_weights, lower_nodes, lower_weights):
        super().__init__(higher_nodes, higher_weights)

        # TODO: not sure this is 100% right, since in e.g. GaussKronrod21 we never call
        # init, and just define the higher_nodes etc as properties.
        self.higher_nodes = higher_nodes
        self.higher_weights = higher_weights
        self.lower_nodes = lower_nodes
        self.lower_weights = lower_weights

    def higher_estimate(self, f, a, b, args):
        # The underlying cubature rule is presumed to be for the hypercube [0, 1]^n.
        #
        # To handle arbitrary regions of integration, it's necessary to apply a linear
        # change of coordinates to map each interval [a[i], b[i]] to [-1, 1].
        a = a[:, np.newaxis]
        b = b[:, np.newaxis]

        lengths = b - a

        # nodes have shape (input_dim, eval_points)
        higher_nodes = (self.higher_nodes + 1) * (lengths / 2) + a

        # Also need to multiply the weights by a scale factor equal to the determinant
        # of the Jacobian for this coordinate change.
        weight_scale_factor = math.prod(lengths / 2)

        # weights have shape (eval_points,)
        higher_weights = self.higher_weights * weight_scale_factor

        # f(nodes) will have shape (eval_points, output_dim)
        # integral_estimate should have shape (output_dim,)
        higher_integral_estimate = np.sum(
            higher_weights * f(higher_nodes, *args),
            axis=-1
        )

        return higher_integral_estimate

    def lower_estimate(self, f, a, b, args):
        a = a[:, np.newaxis]
        b = b[:, np.newaxis]

        lengths = b - a

        lower_nodes = (self.lower_nodes + 1) * (lengths / 2) + a

        weight_scale_factor = math.prod(lengths / 2)
        lower_weights = self.lower_weights * weight_scale_factor

        lower_integral_estimate = np.sum(
            lower_weights * f(lower_nodes, *args),
            axis=-1
        )

        return lower_integral_estimate

    def estimate(self, f, a, b, args):
        return self.higher_estimate(f, a, b, args)

    def error_estimate(self, f, a, b, args):
        return np.abs(
            self.higher_estimate(f, a, b, args) - self.lower_estimate(f, a, b, args),
        )


class GaussKronrod21(NestedCubatureRule):
    def __init__(self):
        pass

    @functools.cached_property
    def higher_nodes(self):
        # 21-point nodes
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

    @functools.cached_property
    def higher_weights(self):
        # 21-point weights
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

    @functools.cached_property
    def lower_nodes(self):
        # 10-point nodes
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

    @functools.cached_property
    def lower_weights(self):
        # 10-point weights
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


class GaussKronrod15(NestedCubatureRule):
    def __init__(self):
        pass

    @functools.cached_property
    def higher_nodes(self):
        # 15-point nodes
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

    @functools.cached_property
    def higher_weights(self):
        # 15-point weights
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

    @functools.cached_property
    def lower_nodes(self):
        # 7-point nodes
        return np.array([[
            0.949107912342758524526189684047851,
            0.741531185599394439863864773280788,
            0.405845151377397166906606412076961,
            0.000000000000000000000000000000000,
            -0.405845151377397166906606412076961,
            -0.741531185599394439863864773280788,
            -0.949107912342758524526189684047851,
        ]])

    @functools.cached_property
    def lower_weights(self):
        # 7-point weights
        return np.array([
            0.129484966168869693270611432679082,
            0.279705391489276667901467771423780,
            0.381830050505118944950369775488975,
            0.417959183673469387755102040816327,
            0.381830050505118944950369775488975,
            0.279705391489276667901467771423780,
            0.129484966168869693270611432679082
        ])


class Trapezoid(NestedCubatureRule):
    def __init__(self):
        pass

    @functools.cached_property
    def higher_nodes(self):
        return np.array([[
            -1,
            0,
            1,
        ]])

    @functools.cached_property
    def higher_weights(self):
        return np.array([
            0.5,
            1,
            0.5
        ])

    @functools.cached_property
    def lower_nodes(self):
        return np.array([[
            -1,
            1
        ]])

    @functools.cached_property
    def error_weights(self):
        return np.array([
            1/2,
            1/2,
        ])


class ProductRule(NestedCubatureRule):
    def __init__(self, base_rules):
        self.base_rules = base_rules

    @functools.cached_property
    def higher_nodes(self):
        return _cartesian_product([rule.higher_nodes for rule in self.base_rules])

    @functools.cached_property
    def higher_weights(self):
        return np.prod(
            _cartesian_product([rule.higher_weights for rule in self.base_rules]),
            axis=0
        )

    @functools.cached_property
    def lower_nodes(self):
        return _cartesian_product([rule.lower_nodes for rule in self.base_rules])

    @functools.cached_property
    def lower_weights(self):
        return np.prod(
            _cartesian_product([rule.lower_weights for rule in self.base_rules]),
            axis=0
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
