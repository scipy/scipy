import math
import itertools

import numpy as np

from functools import cached_property


class Cub:
    """
    Base class for numerical integration algorithms (cubatures).

    Finds an estimate for the integral of ``f`` over hypercube described by two arrays
    ``a`` and ``b`` via `estimate`, and find an estimate for the error of this
    approximation via `error_estimate`.

    If a subclass does not implement its own `error_estimate`, then it will use a
    default error estimate based on the difference between the estimate over the whole
    region and the sum of the estimates over each subregion.

    See Also
    --------
    FixedCub

    Examples
    --------
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
    ...     def error_estimate(self, f, a, b, args=(), kwargs=None):
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
     array([[-0.9635857 ,  0.48159229,  0.79091717],
            [-0.92119516,  0.07230874,  0.02134233]])
    """

    def estimate(self, f, a, b, args=(), kwargs=None):
        r"""
        Calculate estimate of integral of ``f`` in hypercube described by corners ``a``
        and ``b``.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, \*args, \*\*kwargs) -> ndarray

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
        r"""
        Calculate the error estimate of this cubature rule for the integral of ``f`` in
        the hypercube described by corners ``a`` and ``b``.

        If a subclass does not override this method, then a default error estimator is
        used. This estimates the error as ``|est - refined_est|`` where ``est`` is
        ``estimate(f, a, b)`` and ``refined_est`` is the sum of
        ``estimate(f, a_k, b_k)`` where ``a_k, b_k`` are the coordinates of each
        subregion of the hypercube described by ``a`` and ``b``. In the 1D case, this
        is equivalent to comparing the integral over an entire interval ``[a, b]`` to
        the sum of the integrals over the left and right subintervals, ``[a, (a+b)/2]``
        and ``[(a+b)/2, b]``.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, \*args, \*\*kwargs) -> ndarray

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
    A cubature rule implemented as the weighted sum of function evaluations at fixed
    nodes.

    Attributes
    ----------
    nodes_and_weights : (ndarray, ndarray)
        A tuple ``(nodes, weights)`` of nodes at which to evaluate ``f`` and the
        corresponding weights. ``nodes`` should be of shape ``(num_nodes,)`` for 1D
        cubature rules (quadratures) and more generally for N-D cubature rules, it
        should be of shape ``(ndim, num_nodes)``. ``weights`` should be of shape
        `(num_nodes,)`. The nodes and weights should be for integrals over `[-1, 1]^n`.

    See Also
    --------
    NewtonCotesQuad, GaussLegendreQuad, ErrorFromDifference
    """

    @property
    def nodes_and_weights(self):
        raise NotImplementedError

    def estimate(self, f, a, b, args=(), kwargs=None):
        r"""
        Calculate estimate of integral of ``f`` in hypercube described by corners ``a``
        and ``b`` as ``sum(weights * f(nodes))``. Nodes and weights will automatically
        be adjusted from calculating integrals over :math:`[-1, 1]^n` to
        :math:`[a, b]^n`.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, \*args, \*\*kwargs) -> ndarray

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
        nodes, weights = self.nodes_and_weights

        return _apply_fixed_rule(f, a, b, nodes, weights, args, kwargs)


class ErrorFromDifference(FixedCub):
    """
    A fixed cubature rule with error estimate given by the difference between two
    underlying fixed cubature rules.

    This can be used to give error estimation to cubature rules which don't have default
    error estimation, such as NewtonCotesQuad or GaussLegendreQuad. This is in contrast
    to rules like GaussKronrod which have default error estimation given by a special
    property of that cubature rule. See Examples.

    Attributes
    ----------
    higher : FixedCub
        Higher accuracy fixed cubature rule.

    lower : FixedCub
        Lower accuracy fixed cubature rule.

    nodes_and_weights : (ndarray, ndarray)
        Nodes and weights. These are mirrored from ``higher.nodes_and_weights``.

    lower_nodes_and_weights : (ndarray, ndarray)
        Nodes and weights for the lower rule. These are mirrored from
        ``lower.nodes_and_weights``.

    See Also
    --------
    GaussKronrodQuad, NewtonCotesQuad

    Examples
    --------
    Gauss-Legendre rules have no intrinsic error estimation. We can build a cubature
    rule with error estimation by taking the difference between an accurate cubature
    rule and a less accurate cubature rule:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import GaussLegendreQuad, ErrorFromDifference
    >>> higher_accuracy_rule = GaussLegendreQuad(20)
    >>> lower_accuracy_rule = GaussLegendreQuad(10)
    >>> custom_rule = ErrorFromDifference(
    ...     higher_accuracy_rule,
    ...     lower_accuracy_rule
    ... )
    >>> cub(lambda x: np.sin(x), 0, 1, custom_rule).estimate
     array([0.45969769])
    """

    def __init__(self, higher, lower):
        self.higher = higher
        self.lower = lower

    @property
    def nodes_and_weights(self):
        return self.higher.nodes_and_weights

    @property
    def lower_nodes_and_weights(self):
        return self.lower.nodes_and_weights

    def error_estimate(self, f, a, b, args=(), kwargs=None):
        r"""
        Calculate the error estimate of this cubature rule for the integral of ``f`` in
        the hypercube described by corners ``a`` and ``b``.

        The error estimate used is the difference between the approximation given by
        ``lower`` and ``higher``.

        Parameters
        ----------
        f : callable
            Function to integrate. ``f`` must have the signature::
                f(x : ndarray, \*args, \*\*kwargs) -> ndarray

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

        nodes, weights = self.nodes_and_weights
        lower_nodes, lower_weights = self.lower_nodes_and_weights

        error_nodes = np.concat([nodes, lower_nodes], axis=-1)
        error_weights = np.concat([weights, -lower_weights], axis=-1)

        return np.abs(
            _apply_fixed_rule(f, a, b, error_nodes, error_weights, args, kwargs)
        )


class FixedProductCub(FixedCub):
    """
    Find the n-dimensional cubature rule constructed from the Cartesian product of 1D
    fixed cubature rules.

    Parameters
    ----------
    base_rules : list of FixedCub
        List of base 1-dimensional FixedCub cubature rules.

    Attributes
    ----------
    base_rules : list of FixedCub
        List of base 1-dimensional FixedCub cubature rules.

    Examples
    --------

    Evaluate a 2D integral by taking the product of two 1D rules:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import (
    ...  FixedProductCub, NewtonCotesQuad
    ... )
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule = FixedProductCub(
    ...     [NewtonCotesQuad(10), NewtonCotesQuad(10)]
    ... ) # Use 10-point NewtonCotesQuad
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.error_estimate(f, a, b)
     np.float64(2.220446049250313e-16)
    """

    def __init__(self, base_rules):
        self.base_rules = base_rules

    @cached_property
    def nodes_and_weights(self):
        nodes = _cartesian_product(
            [rule.nodes_and_weights[0] for rule in self.base_rules]
        )

        weights = np.prod(_cartesian_product(
            [rule.nodes_and_weights[1] for rule in self.base_rules]
        ), axis=0)

        return nodes, weights


class FixedProductErrorFromDifferenceCub(ErrorFromDifference):
    """
    Find the n-dimensional cubature rule constructed from the Cartesian product of 1D
    `ErrorFromDifference` cubature rules, and estimate the error as the difference
    between the product of the higher rules and the product of the lower rules.

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
    base_rules : list of ErrorFromDifference
        List of base 1-dimensional ErrorFromDifference cubature rules.

    Examples
    --------

    Evaluate a 2D integral by taking the product of two 1D rules:

    >>> import numpy as np
    >>> from scipy.integrate._cubature import cub
    >>> from scipy.integrate._rules import (
    ...  FixedProductErrorFromDifferenceCub, GaussKronrodQuad
    ... )
    >>> def f(x):
    ...     # f(x) = cos(x_1) + cos(x_2)
    ...     return np.sum(np.cos(x), axis=0)
    >>> rule = FixedProductErrorFromDifferenceCub(
    ...     [GaussKronrodQuad(15), GaussKronrodQuad(15)]
    ... ) # Use 15-point GaussKronrodQuad, which implements ErrorFromDifference
    >>> a, b = np.array([0, 0]), np.array([1, 1])
    >>> rule.estimate(f, a, b) # True value 2*sin(1), approximately 1.6829
     np.float64(1.682941969615793)
    >>> rule.error_estimate(f, a, b)
     np.float64(2.220446049250313e-16)
    """

    def __init__(self, base_rules):
        self.base_rules = base_rules

    @cached_property
    def nodes_and_weights(self):
        nodes = _cartesian_product(
            [rule.nodes_and_weights[0] for rule in self.base_rules]
        )

        weights = np.prod(_cartesian_product(
            [rule.nodes_and_weights[1] for rule in self.base_rules]
        ), axis=0)

        return nodes, weights

    @cached_property
    def lower_nodes_and_weights(self):
        nodes = _cartesian_product(
            [cubature.lower_nodes_and_weights[0] for cubature in self.base_rules]
        )

        weights = np.prod(_cartesian_product(
            [cubature.lower_nodes_and_weights[1] for cubature in self.base_rules]
        ), axis=0)

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


def _apply_fixed_rule(f, a, b, orig_nodes, orig_weights, args=(), kwargs=None):
    if kwargs is None:
        kwargs = dict()

    # Ensure orig_nodes are at least 2D, since 1D cubature methods return arrays of
    # shape (num_eval_points,) rather than (1,num_eval_points)
    orig_nodes = np.atleast_2d(orig_nodes)

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
