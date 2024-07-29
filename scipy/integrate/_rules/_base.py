import math
import itertools

import numpy as np

from functools import cached_property


class Cub:
    """
    A generic interface for numerical cubature algorithms.

    Finds an estimate for the integral of ``f`` over hypercube described by two arrays
    ``a`` and ``b`` via ``estimate``, and may also find an estimate for the error of
    this approximation via ``error_estimate``.

    If error estimation is not supported by default (as in the case of Gauss-Legendre
    or Newton-Cotes), then ``error_estimate`` will raise a NotImplementedError.
    """

    def estimate(self, f, a, b, args=(), kwargs=None):
        """
        Calculate estimate of integral of ``f`` in hypercube described by corners ``a``
        and ``b``.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, *args, **kwargs) -> ndarray
            If ``f`` accepts arrays ``x`` of shape ``(input_dim_1, ..., input_dim_n,
            num_eval_points)`` and outputs arrays of shape ``(output_dim_1, ...,
            output_dim_n, num_eval_points)``, then ``cub`` will return arrays of shape
            ``(output_dim_1, ..., output_dim_n)``.
        a, b : ndarray
            Lower and upper limits of integration as rank-1 arrays specifying the left
            and right endpoints of the intervals being integrated over. Infinite limits
            are currently not supported.
        args : tuple, optional
            Additional positional args passed to ``f``, if any.
        kwargs : tuple, optional
            Additional keyword args passed to ``f``, if any.

        Returns
        -------
        est : ndarray
            Result of estimation. If ``f`` returns arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``, then ``err_est`` will be of shape
            ``(output_dim_1, ..., output_dim_n)``.
        """
        raise NotImplementedError

    def error_estimate(self, f, a, b, args=(), kwargs=None):
        """
        Calculate the error estimate of this cubature rule for the integral of ``f`` in
        the hypercube described by corners ``a`` and ``b``.

        If the cubature rule doesn't support error estimation, this will raise a
        NotImplementedError.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, *args, **kwargs) -> ndarray
            If ``f`` accepts arrays ``x`` of shape ``(input_dim_1, ..., input_dim_n,
            num_eval_points)`` and outputs arrays of shape ``(output_dim_1, ...,
            output_dim_n, num_eval_points)``, then ``cub`` will return arrays of shape
            ``(output_dim_1, ..., output_dim_n)``.
        a, b : ndarray
            Lower and upper limits of integration as rank-1 arrays specifying the left
            and right endpoints of the intervals being integrated over. Infinite limits
            are currently not supported.
        args : tuple, optional
            Additional positional args passed to ``f``, if any.
        kwargs : tuple, optional
            Additional keyword args passed to ``f``, if any.

        Returns
        -------
        err_est : ndarray
            Estimate of the error. If the cubature rule doesn't support error
           estimation, then a NotImplementedError will be raised instead. If error
           estimation is supported and ``f`` returns arrays of shape ``(output_dim_1,
           ..., output_dim_n, eval_points)``, then ``err_est`` will be of shape
           ``(output_dim_1, ..., output_dim_n)``.
        """

        est = self.estimate(f, a, b, args, kwargs)
        refined_est = 0

        for a_k, b_k in _subregion_coordinates(a, b):
            refined_est += self.estimate(f, a_k, b_k, args, kwargs)

        return np.abs(est - refined_est)


class FixedCub(Cub):
    """
    A cubature rule implemented as the weighted sum of function evaluations.

    Attributes
    ----------
    TODO: rule attribute

    See Also
    --------
    NewtonCotes, GaussLegendre, ErrorFromDifference
    """

    @property
    def rule(self):
        raise NotImplementedError

    def estimate(self, f, a, b, args=(), kwargs=None):
        """
        Calculate estimate of integral of ``f`` in hypercube described by corners ``a``
        and ``b`` as ``sum(weights * f(nodes))``. Nodes and weights will automatically
        be adjusted from calculating integrals over :math:`[-1, 1]^n` to
        :math:`[a, b]^n`.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, *args, **kwargs) -> ndarray
            If ``f`` accepts arrays ``x`` of shape ``(input_dim_1, ..., input_dim_n,
            num_eval_points)`` and outputs arrays of shape ``(output_dim_1, ...,
            output_dim_n, num_eval_points)``, then ``cub`` will return arrays of shape
            ``(output_dim_1, ..., output_dim_n)``.
        a, b : ndarray
            Lower and upper limits of integration as rank-1 arrays specifying the left
            and right endpoints of the intervals being integrated over. Infinite limits
            are currently not supported.
        args : tuple, optional
            Additional positional args passed to ``f``, if any.
        kwargs : tuple, optional
            Additional keyword args passed to ``f``, if any.

        Returns
        -------
        est : ndarray
            Result of estimation. If ``f`` returns arrays of shape ``(output_dim_1,
            ..., output_dim_n, eval_points)``, then ``err_est`` will be of shape
            ``(output_dim_1, ..., output_dim_n)``.
        """
        return _apply_rule(f, a, b, self.rule, args, kwargs)


class ErrorFromDifference(FixedCub):
    """
    A fixed cubature rule with error estimate given by the difference between two
    underlying fixed cubature rules.

    This can be used to give error estimation to cubature rules which don't have default
    error estimation, such as NewtonCotes or GaussLegendre. This is in contrast to rules
    like GaussKronrod which have default error estimation given by a special property of
    that cubature rule. See Examples.

    Attributes
    ----------
    TODO: rule attribute and lower_rule attribute

    See Also
    --------
    GaussKronrod, NewtonCotes

    Examples
    --------
    TODO: examples for ErrorFromDifference
    """

    def __init__(self, higher, lower):
        self.higher = higher
        self.lower = lower

    @property
    def rule(self):
        return self.higher.rule

    @property
    def lower_rule(self):
        return self.lower.rule

    def error_estimate(self, f, a, b, args=(), kwargs=None):
        """
        TODO: docstring for error_estimate in ErrorFromDifference
        """
        nodes, weights = self.rule
        lower_nodes, lower_weights = self.lower_rule

        error_nodes = np.concat([nodes, lower_nodes], axis=-1)
        error_weights = np.concat([weights, -lower_weights], axis=-1)
        error_rule = (error_nodes, error_weights)

        return np.abs(
            _apply_rule(f, a, b, error_rule, args, kwargs)
        )


class FixedProductCub(ErrorFromDifference):
    """
    Find the n-dimensional cubature rule constructed from the Cartesian product of 1D
    cubature rules.

    Given a list of N 1-dimensional cubature rules which support error estimation using
    ErrorFromDifference, this will find the N-dimensional ErrorFromDifference cubature
    rule obtained by taking the Cartesian product of their nodes, and estimating the
    error by taking the difference with a lower-accuracy N-dimensional cubature rule
    obtained using the `.lower` rule in each of the base 1-dimensional rules.

    Parameters
    ----------
    base_rules : list of ErrorFromDifference
        List of base 1-dimensional ErrorFromDifference cubature rules.

    Attributes
    ----------
    base_rules : list of DualEstimateCubature
        List of base 1-dimensional ErrorFromDifference cubature rules.

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

    def __init__(self, base_cubatures):
        self.base_cubatures = base_cubatures

    @cached_property
    def rule(self):
        underlying_nodes = [cubature.rule[0] for cubature in self.base_cubatures]

        to_product = []

        for node_group in underlying_nodes:
            if len(node_group.shape) == 1:
                to_product.append(node_group)
            else:
                to_product.extend(node_group)

        nodes = _cartesian_product(
            [cubature.rule[0] for cubature in self.base_cubatures]
        )

        weights = np.prod(
            _cartesian_product(
                [cubature.rule[1] for cubature in self.base_cubatures]
            ),
            axis=0
        )

        return nodes, weights

    @cached_property
    def lower_rule(self):
        nodes = _cartesian_product(
            [cubature.lower_rule[0] for cubature in self.base_cubatures]
        )

        weights = np.prod(
            _cartesian_product(
                [cubature.lower_rule[1] for cubature in self.base_cubatures]
            ),
            axis=0
        )

        return nodes, weights


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


def _apply_rule(f, a, b, rule, args=(), kwargs=None):
    if kwargs is None:
        kwargs = dict()

    orig_nodes, orig_weights = rule

    rule_ndim = orig_nodes.shape[0]
    a_ndim = len(a)
    b_ndim = len(b)

    if rule_ndim != a_ndim or rule_ndim != b_ndim:
        raise ValueError(f"cubature rule and function are of incompatible dimension, \
nodes have ndim {rule_ndim}, while limit of integration have ndim \
a_ndim={a_ndim}, b_ndim={b_ndim}")

    # Since f accepts arrays of shape (ndim, eval_points), it is necessary to
    # add an extra axis to a and b so that ``f`` can be evaluated there.
    a = a[:, np.newaxis]
    b = b[:, np.newaxis]
    lengths = b - a

    # The underlying cubature rule is for the hypercube [-1, 1]^n.
    #
    # To handle arbitrary regions of integration, it's necessary to apply a linear
    # change of coordinates to map each interval [a[i], b[i]] to [-1, 1].
    nodes = (orig_nodes + 1) * (lengths / 2) + a

    # Also need to multiply the weights by a scale factor equal to the determinant
    # of the Jacobian for this coordinate change.
    weight_scale_factor = math.prod(lengths / 2)
    weights = orig_weights * weight_scale_factor

    # f(nodes) will have shape (output_dim_1, ..., output_dim_n, num_nodes)
    # Summing along the last axis means estimate will shape (output_dim_1, ...,
    # output_dim_n)
    est = np.sum(
        weights * f(nodes, *args, **kwargs),
        axis=-1
    )

    return est
