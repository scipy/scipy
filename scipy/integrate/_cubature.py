import heapq
import itertools

import numpy as np


class CubatureRule:
    def __init__(self, nodes, weights):
        self.nodes = nodes
        self.weights = weights

    # returns the dimension of the cubature rule
    # TODO: need to be more precise about what this means
    def dimension(self):
        return self.nodes.shape[-1]

    # returns integral_estimate, error_estimate
    def estimate(self):
        pass


class CubatureRuleLinearError(CubatureRule):
    def __init__(self, nodes, weights, error_nodes, error_weights):
        super().__init__(nodes, weights)

        self.error_nodes = error_nodes
        self.error_weights = error_weights

    # ranges looks like np.array([[0, 1], [2, 3]]) for
    # \int^0_1 \int_2^3 f(x, y) dydx
    # f is assumed to be a function from a function from an array of shape
    # (eval_points, dim_input) where dim_input is the dimension of the input
    # and should return an array of shape (eval_points, dim_output) where dim_output
    # is the dimension of the output, e.g. for f(x, y) = (x + y)^alpha where
    # alpha is one of [.2, .5, .8] would have dim_input=2 and dim_output=3
    # TODO: consider better name than ranges, e.g. look at nquad
    def estimate(self, f, ranges):
        # TODO: need to output some sort of error when the rule and f have inconsistent
        # dimensions

        range_start_points = ranges[:, 0]
        range_end_points = ranges[:, 1]
        range_lengths = range_end_points - range_start_points

        # The underlying cubature rule is presumed to work on the interval (-1, 1)
        # We need to apply a linear change of coordinates to map each interval (a, b)
        # to (-1, 1).
        # We also need to multiply the weights by a scale factor equal to the
        # determinant of the Jacobian for this change.
        weight_scale_factor = np.prod(range_lengths / 2)

        nodes_scaled = (self.nodes + 1) * (range_lengths / 2) + range_start_points
        error_nodes_scaled = (
            (self.error_nodes + 1) * (range_lengths / 2) + range_start_points
        )

        weights_scaled = self.weights * weight_scale_factor
        error_weights_scaled = self.error_weights * weight_scale_factor

        # Need to resize the weights so they have the right shape
        integral_estimate = np.sum(
            weights_scaled[:, np.newaxis] * f(nodes_scaled),
            axis=0
        )
        error_estimate = abs(np.sum(
            error_weights_scaled[:, np.newaxis] * f(error_nodes_scaled),
            axis=0
        ))

        # TODO: not yet worrying about rounding error

        return integral_estimate, error_estimate


# For now, assuming that the underlying quadrature rules all have linear error
# estimators so that we can combine the error nodes and error weights
class ProductRule(CubatureRuleLinearError):
    def __init__(self, base_rules):
        # TODO: using meshgrid like this probably isn't right
        self.nodes = _cartesian_product(
            [np.meshgrid(rule.nodes)[0] for rule in base_rules]
        )
        self.error_nodes = _cartesian_product(
            [np.meshgrid(rule.error_nodes)[0] for rule in base_rules]
        )

        self.weights = np.prod(
            _cartesian_product([rule.weights for rule in base_rules]),
            axis=1
        )
        self.error_weights = np.prod(
            _cartesian_product([rule.error_weights for rule in base_rules]),
            axis=1
        )


# Maybe in the future it would be better if we could specify the degree and then it
# would calculate the weights.
class GaussKronrod21(CubatureRuleLinearError):
    def __init__(self):
        # 21-point nodes
        nodes_21 = np.array([
            [0.995657163025808080735527280689003],
            [0.973906528517171720077964012084452],
            [0.930157491355708226001207180059508],
            [0.865063366688984510732096688423493],
            [0.780817726586416897063717578345042],
            [0.679409568299024406234327365114874],
            [0.562757134668604683339000099272694],
            [0.433395394129247190799265943165784],
            [0.294392862701460198131126603103866],
            [0.148874338981631210884826001129720],
            [0],
            [-0.148874338981631210884826001129720],
            [-0.294392862701460198131126603103866],
            [-0.433395394129247190799265943165784],
            [-0.562757134668604683339000099272694],
            [-0.679409568299024406234327365114874],
            [-0.780817726586416897063717578345042],
            [-0.865063366688984510732096688423493],
            [-0.930157491355708226001207180059508],
            [-0.973906528517171720077964012084452],
            [-0.995657163025808080735527280689003]
        ])

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
        # (indexed from 0). At the moment I'll just hard code this
        nodes_10 = np.array([
            nodes_21[1],
            nodes_21[3],
            nodes_21[5],
            nodes_21[7],
            nodes_21[9],
            nodes_21[11],
            nodes_21[13],
            nodes_21[15],
            nodes_21[17],
            nodes_21[19],
        ])

        # 10-point weights
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

        # TODO: avoid duplicates, this is not good since the sample points increase
        # exponentially with dimension
        self.error_nodes = np.concat([nodes_21, nodes_10])
        self.error_weights = np.concat([weights_21, -weights_10])


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


def split_range(range):
    a, b = range
    mid = (a + b) / 2

    return np.array([[a, mid], [mid, b]])


# TODO: also need to consider relative tolerance
# ranges is like [[a_x, b_x], [a_y, b_y]]
def cub(f, ranges, rule, tol):
    global_est, global_err = rule.estimate(f, ranges)

    regions = [RegionInfo(global_est, global_err, ranges)]

    # print(global_err)

    while _max_norm(global_err) > tol:
        region = heapq.heappop(regions)
        print(_max_norm(global_err), global_err, region.ranges)
        # print("Considering region", region)

        est_k = region.integral_estimate
        err_k = region.error_estimate
        ranges_k = region.ranges

        split_ranges = []

        for rng in ranges_k:
            split_ranges.append(split_range(rng))

        # print("split ranges:", split_ranges)

        # this hurts my eyes
        # sub_ranges is now like
        #     array([[[-1.,  0.],
        #     [-1.,  0.]],

        #    [[-1.,  0.],
        #     [ 0.,  1.]],

        #    [[ 0.,  1.],
        #     [-1.,  0.]],

        #    [[ 0.,  1.],
        #     [ 0.,  1.]]])
        all_sub_ranges = np.array(list(itertools.product(*split_ranges)))

        # print("all sub ranges:", all_sub_ranges)

        global_est -= est_k
        global_err -= err_k

        for sub_ranges in all_sub_ranges:
            # print(sub_ranges)
            est_sub, err_sub = rule.estimate(f, sub_ranges)

            global_est += est_sub
            global_err += err_sub

            new_region = RegionInfo(est_sub, err_sub, sub_ranges)
            # print("attempting to push", new_region, "onto heap")
            # print("regions is currently", regions)

            heapq.heappush(regions, new_region)

    return global_est, global_err


# Not sure this is really needed, was running into issues comparing very small floats
# in the heap so made this to make it more explicit
class RegionInfo:
    def __init__(self, integral_estimate, error_estimate, ranges):
        self.integral_estimate = integral_estimate
        self.error_estimate = error_estimate
        self.ranges = ranges

    def key(self):
        return -np.linalg.norm(abs(self.error_estimate))

    def __le__(self, other):
        # this is comparing -self.error_estimate against -other.error_estimate
        return self.key() < other.key()

    def __lt__(self, other):
        return self.__le__(other)

    def __repr__(self):
        return f"""Estimate: {self.integral_estimate}
Errors: {self.error_estimate} (max {-self.key()})
Ranges: {self.ranges}"""

    def __str__(self):
        return self.__repr__()


# Slight modification of https://stackoverflow.com/a/11146645
# TODO: Is there a more canonical way of finding the Cartesian product of two numpy
# arrays?
def _cartesian_product(arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[..., i] = a
    return arr.reshape(-1, la)


def _max_norm(x):
    return np.amax(abs(x))
