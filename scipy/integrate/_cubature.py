import math
import heapq
import itertools
import functools

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from numpy.polynomial.polynomial import Polynomial


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
        # Consider regions with higher error as smaller so that they are placed at the
        # top of the min-heap.
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


class CubatureRule(ABC):
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
    def _estimate_from_nodes(self, nodes, weights, f, a, b, args=()):
        # The underlying cubature rule is presumed to be for the hypercube [0, 1]^n.
        #
        # To handle arbitrary regions of integration, it's necessary to apply a linear
        # change of coordinates to map each interval [a[i], b[i]] to [-1, 1].
        a = a[:, np.newaxis]
        b = b[:, np.newaxis]

        lengths = b - a

        # nodes have shape (input_dim, eval_points)
        nodes = (nodes + 1) * (lengths / 2) + a

        # Also need to multiply the weights by a scale factor equal to the determinant
        # of the Jacobian for this coordinate change.
        weight_scale_factor = math.prod(lengths / 2)

        # weights have shape (eval_points,)
        weights = weights * weight_scale_factor

        # f(nodes) will have shape (eval_points, output_dim)
        # integral_estimate should have shape (output_dim,)
        estimate = np.sum(
            weights * f(nodes, *args),
            axis=-1
        )

        return estimate

    def higher_estimate(self, f, a, b, args=()):
        return self._estimate_from_nodes(
            self.higher_nodes,
            self.higher_weights,
            f,
            a,
            b,
            args,
        )

    def lower_estimate(self, f, a, b, args=()):
        return self._estimate_from_nodes(
            self.lower_nodes,
            self.lower_weights,
            f,
            a,
            b,
            args,
        )

    def estimate(self, f, a, b, args=()):
        return self.higher_estimate(f, a, b, args)

    def error_estimate(self, f, a, b, args):
        # Take the difference between the higher and lower estimate to obtain estimate
        # for the error.
        return np.abs(self._estimate_from_nodes(
            np.concat((self.higher_nodes, self.lower_nodes), axis=-1),
            np.concat((self.higher_weights, -self.lower_weights), axis=-1),
            f,
            a,
            b,
            args
        ))


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


class NewtonCotes(NestedCubatureRule):
    def __init__(self, npoints, open=False):
        if npoints <= 2:
            raise Exception(
                "At least 3 points required for trapezoid rule with error estimate"
            )

        self.npoints = npoints
        self.open = open

    @functools.cached_property
    def higher_nodes(self):
        if self.open:
            h = 2/self.npoints
            return np.linspace(-1 + h, 1 - h, num=self.npoints).reshape(1, -1)
        else:
            return np.linspace(-1, 1, num=self.npoints).reshape(1, -1)

    @functools.cached_property
    def higher_weights(self):
        return self._weights_from_nodes(self.higher_nodes.reshape(-1))

    @functools.cached_property
    def lower_nodes(self):
        if self.open:
            h = 2/(self.npoints-1)
            return np.linspace(-1 + h, 1 - h, num=self.npoints-1).reshape(1, -1)
        else:
            return np.linspace(-1, 1, num=self.npoints-1).reshape(1, -1)

    @functools.cached_property
    def lower_weights(self):
        return self._weights_from_nodes(self.lower_weights.reshape(-1))

    def _weights_from_nodes(self, nodes):
        """
        Calculates the weight array from a list of nodes by forming the Lagrange basis
        polynomial for each node x_i and then integrating over the interval [-1, 1].
        """

        weights = []

        for node in nodes:
            poly = _lagrange_basis_polynomial(node, nodes)
            indef = poly.integ()
            weights.append(indef(1) - indef(-1))

        return np.array(weights)


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


class GenzMalik(NestedCubatureRule):
    def __init__(self, ndim):
        if ndim < 2:
            raise Exception("Genz-Malik cubature is only defined for ndim >= 2")

        self.ndim = ndim

    @functools.cached_property
    def higher_nodes(self):
        l_2 = np.sqrt(9/70)
        l_3 = np.sqrt(9/10)
        l_4 = np.sqrt(9/10)
        l_5 = np.sqrt(9/19)

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

    @functools.cached_property
    def higher_weights(self):
        w_1 = (2**self.ndim) * (12824 - 9120 * self.ndim + 400 * self.ndim**2) / 19683
        w_2 = (2**self.ndim) * 980/6561
        w_3 = (2**self.ndim) * (1820 - 400 * self.ndim) / 19683
        w_4 = (2**self.ndim) * (200 / 19683)
        w_5 = 6859 / 19683

        return np.repeat(
            [w_1, w_2, w_3, w_4, w_5],
            [1, 2 * self.ndim, 2*self.ndim, 2*(self.ndim - 1)*self.ndim, 2**self.ndim]
        )

    @functools.cached_property
    def lower_nodes(self):
        out_size = 1 + 2 * (self.ndim + 1) * self.ndim
        out = self.higher_nodes[:, :out_size]

        return out

    @functools.cached_property
    def lower_weights(self):
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
    Find the number of distinct permutations of `r` elements of the iterable.
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


def _lagrange_basis_polynomial(x_i, xs):
    poly = Polynomial(1)

    for x_j in xs:
        if x_i == x_j:
            continue

        poly *= Polynomial([-x_j, 1])/(x_i - x_j)

    return poly
