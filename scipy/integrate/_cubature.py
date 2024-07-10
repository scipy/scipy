import math
import heapq

import numpy as np


class CubatureRule:
    def __init__(self, nodes, weights):
        self.nodes = nodes
        self.weights = weights

    # returns the input dimension that the cubature rule is expecting `f` to have
    def dimension(self):
        return self.nodes.shape[-1]

    # returns integral_estimate, error_estimate
    def estimate(self, f, a, b):
        raise NotImplementedError


class CubatureRuleLinearError(CubatureRule):
    def __init__(self, nodes, weights, error_nodes, error_weights):
        super().__init__(nodes, weights)

        self.error_nodes = error_nodes
        self.error_weights = error_weights

    def estimate(self, f, a, b):
        # TODO: need to output some sort of error when the rule and f have inconsistent
        # dimensions

        # The underlying cubature rule is presumed to be for the hypercube [0, 1]^n.
        #
        # To handle arbitrary regions of integration, it's necessary to apply a linear
        # change of coordinates to map each interval [a[i], b[i]] to [-1, 1].
        a = a[:, np.newaxis]
        b = b[:, np.newaxis]

        lengths = b - a

        nodes = (self.nodes + 1) * (lengths / 2) + a
        error_nodes = (self.error_nodes + 1) * (lengths / 2) + a

        # Also need to multiply the weights by a scale factor equal to the determinant
        # of the Jacobian for this coordinate change.
        weight_scale_factor = math.prod(lengths / 2)

        weights = self.weights * weight_scale_factor
        error_weights = self.error_weights * weight_scale_factor

        # Need to resize the weights so they have the right shape
        integral_estimate = np.sum(weights * f(nodes), axis=0)

        # print(integral_estimate)

        error_estimate = abs(np.sum(
            error_weights * f(error_nodes),
            axis=0
        ))

        # print(error_estimate)

        # TODO: not yet worrying about rounding error

        return integral_estimate, error_estimate


# For now, assuming that the underlying quadrature rules all have linear error
# estimators so that we can combine the error nodes and error weights
class ProductRule(CubatureRuleLinearError):
    def __init__(self, base_rules):
        self.nodes = _cartesian_product([rule.nodes for rule in base_rules])
        self.error_nodes = _cartesian_product([rule.error_nodes for rule in base_rules])

        self.weights = np.prod(
            _cartesian_product([rule.weights for rule in base_rules]),
            axis=0
        )
        self.error_weights = np.prod(
            _cartesian_product([rule.error_weights for rule in base_rules]),
            axis=0
        )


class GaussKronrod21(CubatureRuleLinearError):
    def __init__(self):
        # 21-point nodes
        nodes_21 = np.array([[
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

        # 21-point weights
        weights_21 = np.array([
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

        # 10-point nodes come from taking the 21-point nodes at the positions:
        #   1, 3, 5, 7, 9, 11, 13, 15, 17, 19
        # (indexed from 0). The 10-point weights are
        weights_10 = np.array([
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

        # In a Gauss-Kronrod rule, the full set of weights is used for the estimate of
        # the integral, and then the difference between the estimate using the higher
        # and lower order rule is used as the estimation of the error.

        self.nodes = nodes_21
        self.weights = weights_21

        # Error nodes
        self.error_nodes = nodes_21
        self.error_weights = np.array([
            weights_21[0],
            weights_21[1] - weights_10[0],
            weights_21[2],
            weights_21[3] - weights_10[1],
            weights_21[4],
            weights_21[5] - weights_10[2],
            weights_21[6],
            weights_21[7] - weights_10[3],
            weights_21[8],
            weights_21[9] - weights_10[4],
            weights_21[10],
            weights_21[11] - weights_10[5],
            weights_21[12],
            weights_21[13] - weights_10[6],
            weights_21[14],
            weights_21[15] - weights_10[7],
            weights_21[16],
            weights_21[17] - weights_10[8],
            weights_21[18],
            weights_21[19] - weights_10[9],
            weights_21[20],
        ])


class GaussKronrod15(CubatureRuleLinearError):
    def __init__(self):
        # 15-point nodes
        nodes_15 = np.array([[
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

        # 15-point weights
        weights_15 = np.array([
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

        # 7-point nodes come from taking the 15-point nodes at the positions:
        #   1, 3, 5, 7, 9, 11, 13
        # (indexed from 0). The 7-point weights are
        weights_7 = np.array([
            0.129484966168869693270611432679082,
            0.279705391489276667901467771423780,
            0.381830050505118944950369775488975,
            0.417959183673469387755102040816327,
            0.381830050505118944950369775488975,
            0.279705391489276667901467771423780,
            0.129484966168869693270611432679082
        ])

        self.nodes = nodes_15
        self.weights = weights_15

        self.error_nodes = nodes_15
        self.error_weights = np.array([
            weights_15[0],
            weights_15[1] - weights_7[0],
            weights_15[2],
            weights_15[3] - weights_7[1],
            weights_15[4],
            weights_15[5] - weights_7[2],
            weights_15[6],
            weights_15[7] - weights_7[3],
            weights_15[8],
            weights_15[9] - weights_7[4],
            weights_15[10],
            weights_15[11] - weights_7[5],
            weights_15[12],
            weights_15[13] - weights_7[6],
            weights_15[14],
        ])


class Trapezoid(CubatureRuleLinearError):
    def __init__(self):
        self.nodes = np.array([
            [-1],
            [0],
            [1]
        ])
        self.weights = np.array([
            0.5,
            1,
            0.5
        ])

        # These come from the difference using 2 trapezoidal subregions versus three
        # 1 trapezoidal subregions as an estimate for the error.
        self.error_nodes = np.array([
            [-1],
            [0],
            [1]
        ])
        self.error_weights = 1/3 * np.array([
            -0.5,
            1,
            -0.5,
        ])


# TODO: is this a sensible default for maxiter?
def cub(f, a, b, rule, rtol=1e-05, atol=1e-08, maxiter=10_000):
    est, err = rule.estimate(f, a, b)

    regions = [RegionInfo(est, err, a, b)]
    iterations = 1

    while _max_norm(err) > max(atol, rtol * np.abs(est)) and iterations < maxiter:
        region = heapq.heappop(regions)
        print(region)

        est_k = region.integral_estimate
        err_k = region.error_estimate

        a_k, b_k = region.a, region.b
        m_k = (a_k + b_k) / 2

        # TODO: find a cleaner way of finding coordinates of subregions
        subregion_coordinates = list(zip(
            _cartesian_product(np.array([a_k, m_k]).T).T,
            _cartesian_product(np.array([m_k, b_k]).T).T
        ))
        print("subregion coordinates", subregion_coordinates)

        est -= est_k
        err -= err_k

        for a_k_sub, b_k_sub in subregion_coordinates:
            est_sub, err_sub = rule.estimate(f, a_k_sub, b_k_sub)

            est += est_sub
            err += err_sub

            new_region = RegionInfo(est_sub, err_sub, a_k_sub, b_k_sub)

            heapq.heappush(regions, new_region)

        iterations += 1

    if err > max(atol, rtol * np.abs(est)):
        status = "not_converged"
    else:
        status = "converged"

    return CubatureResult(
        estimate=est,
        error=err,
        success=err <= max(atol, rtol * np.abs(est)),
        status=status,
        subdivisions=iterations,
        regions=regions,
    )


class RegionInfo:
    def __init__(self, integral_estimate, error_estimate, a, b):
        self.integral_estimate = integral_estimate
        self.error_estimate = error_estimate
        self.a = a
        self.b = b

    def neg_err(self):
        return -np.linalg.norm(abs(self.error_estimate))

    def __le__(self, other):
        return self.neg_err() < other.neg_err()

    def __lt__(self, other):
        return self.__le__(other)

    def __repr__(self):
        return f"""RegionInfo(est=({self.integral_estimate}) \
err={self.error_estimate}, \
max_err={-self.neg_err()}, \
a={self.a}
b={self.b}"""

    def __str__(self):
        return self.__repr__()


class CubatureResult:
    # TODO: how to handle reporting number of evals for vectorized functions?
    def __init__(self, estimate, error, success, status, subdivisions,
                 regions):
        self.estimate = estimate
        self.error = error
        self.success = success
        self.status = status
        self.regions = regions
        self.subdivisions = subdivisions


def _cartesian_product(points):
    out = np.stack(np.meshgrid(*points, indexing='ij'), axis=0)
    out = out.reshape(len(points), -1)
    return out


def _max_norm(x):
    return np.amax(abs(x))
