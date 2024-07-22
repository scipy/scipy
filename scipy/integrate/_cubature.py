import math
import heapq
import itertools

from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property

import numpy as np


def cub(f, a, b, rule, rtol=1e-05, atol=1e-08, max_subdivisions=10000, args=()):
    if max_subdivisions is None:
        max_subdivisions = np.inf

    est = rule.estimate(f, a, b, args)
    err = rule.error_estimate(f, a, b, args)

    if err is None:
        # TODO: more descriptive error
        raise Exception("Attempting cubature with a rule that doesn't implement error \
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

        # Find all 2**ndim subregions formed by splitting region_k along each axis
        # e.g. for 1D quadrature this splits an interval in two, for 3D cubature this
        # splits the current cube under consideration into 8 subcubes.
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


class Cubature(ABC):
    @abstractmethod
    def estimate(self, f, a, b, args=()):
        """
        Calculate estimate of integral of f in region described by corners a and b.

        f should accept arrays of shape (input_dim, num_eval_points) and return arrays
        of shape (output_dim_1, ..., output_dim_n, num_eval_points).
        """
        pass

    @abstractmethod
    def error_estimate(self, f, a, b, args=()):
        """
        Calculate error estimate for the integral of f in region described by corners a
        and b.

        f should accept arrays of shape (input_dim, num_eval_points) and return arrays
        of shape (output_dim_1, ..., output_dim_n, num_eval_points).
        """
        pass


class FixedCubature(Cubature):
    @property
    def nodes(self):
        raise NotImplementedError

    @property
    def weights(self):
        raise NotImplementedError

    def estimate(self, f, a, b, args=()):
        # The underlying cubature rule is presumed to be for the hypercube [0, 1]^n.
        #
        # To handle arbitrary regions of integration, it's necessary to apply a linear
        # change of coordinates to map each interval [a[i], b[i]] to [-1, 1].
        a = a[:, np.newaxis]
        b = b[:, np.newaxis]

        lengths = b - a

        # nodes have shape (input_dim, eval_points)
        nodes = (self.nodes + 1) * (lengths / 2) + a

        # Also need to multiply the weights by a scale factor equal to the determinant
        # of the Jacobian for this coordinate change.
        weight_scale_factor = math.prod(lengths / 2)

        # weights have shape (eval_points,)
        weights = self.weights * weight_scale_factor

        # f(nodes) will have shape (eval_points, output_dim)
        # integral_estimate should have shape (output_dim,)
        estimate = np.sum(
            weights * f(nodes, *args),
            axis=-1
        )

        return estimate

    # TODO: evaluate whether this is correct
    def error_estimate(self, f, a, b, args=()):
        return None


class DualEstimateCubature(Cubature):
    """
    A cubature rule with a higher-order and lower-order estimate, and where the
    difference between the two is used to estimate the error.
    """

    def __init__(self, higher, lower):
        self.higher = higher
        self.lower = lower

    def estimate(self, f, a, b, args=()):
        return self.higher.estimate(f, a, b, args)

    def error_estimate(self, f, a, b, args=()):
        # Take the difference between the higher and lower estimate to obtain estimate
        # for the error.
        return np.abs(
            self.higher.estimate(f, a, b, args) - self.lower.estimate(f, a, b, args)
        )


class GaussKronrod(DualEstimateCubature):
    def __init__(self, npoints):
        if npoints != 15 and npoints != 21:
            raise Exception("Gauss-Kronrod quadrature is currently only supported for \
                            15 or 21 nodes")

        self.higher = self._Higher(npoints)
        self.lower = self._Lower(npoints)

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
            if self.npoints == 21:
                # 21//2 = 10-point nodes
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
            elif self.npoints == 15:
                # 15//2 = 7-point nodes
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
            if self.npoints == 21:
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
            else:
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


class NewtonCotes(FixedCubature):
    def __init__(self, npoints, open=False):
        if npoints < 2:
            # TODO: check if there really is a 1-point Newton-Cotes rule?
            raise Exception(
                "At least 2 points required for Newton Cotes"
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
