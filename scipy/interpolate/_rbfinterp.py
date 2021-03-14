"""Module for RBF interpolation"""
import warnings
from functools import lru_cache
from itertools import combinations_with_replacement

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import binom

from ._rbfinterp_pythran import _build_system, _evaluate


__all__ = ['RBFInterpolator', 'KNearestRBFInterpolator']


@lru_cache()
def _monomial_powers(ndim, degree):
    """
    Returns the powers for each monomial in a polynomial with the specified
    number of dimensions and degree
    """
    out = []
    for deg in range(degree + 1):
        for itm in combinations_with_replacement(np.eye(ndim, dtype=int), deg):
            out.append(sum(itm, np.zeros(ndim, dtype=int)))

    if not out:
        out = np.zeros((0, ndim), dtype=int)
    else:
        out = np.array(out)

    return out


# The RBFs that are implemented
_AVAILABLE = {'linear', 'tps', 'cubic', 'quintic', 'mq', 'imq', 'iq', 'ga'}


# The shape parameter does not need to be specified when using these RBFs
_SCALE_INVARIANT = {'linear', 'tps', 'cubic', 'quintic'}


# For RBFs that are conditionally positive definite of order m, the interpolant
# should include polynomial terms with degree >= m - 1. Define the minimum
# degrees here. These values are from Chapter 8 of Fasshauer's "Meshfree
# Approximation Methods with MATLAB". The RBFs that are not in this dictionary
# are positive definite and do not need polynomial terms
_NAME_TO_MIN_DEGREE = {
    'mq': 0,
    'linear': 0,
    'tps': 1,
    'cubic': 1,
    'quintic': 2
    }


def _sanitize_init_args(y, d, smoothing, kernel, epsilon, degree, k):
    """
    Sanitize __init__ arguments for RBFInterpolator and KNearestRBFInterpolator
    """
    y = np.asarray(y, dtype=float, order='C')
    if y.ndim != 2:
        raise ValueError('`y` must be a 2-dimensional array')

    ny, ndim = y.shape

    dtype = complex if np.iscomplexobj(d) else float
    d = np.asarray(d, dtype=dtype, order='C')
    if d.shape[0] != ny:
        raise ValueError(
            'Expected the first axis of `d` to have length %d' % ny
            )

    if np.isscalar(smoothing):
        smoothing = np.full(ny, smoothing, dtype=float)
    else:
        smoothing = np.asarray(smoothing, dtype=float, order='C')
        if smoothing.shape != (ny,):
            raise ValueError(
                'Expected `smoothing` to be a scalar or have shape (%d,)' % ny
                )

    if kernel not in _AVAILABLE:
        raise ValueError('`kernel` must be one of %s' % _AVAILABLE)

    if epsilon is None:
        if kernel in _SCALE_INVARIANT:
            epsilon = 1.0
        else:
            raise ValueError(
                '`epsilon` must be specified if `kernel` is not one of %s' %
                _SCALE_INVARIANT
                )
    else:
        epsilon = float(epsilon)

    min_degree = _NAME_TO_MIN_DEGREE.get(kernel, -1)
    if degree is None:
        degree = max(min_degree, 0)
    else:
        degree = int(degree)
        if degree < -1:
            raise ValueError('`degree` must be at least -1')
        elif degree < min_degree:
            warnings.warn(
                '`degree` should not be below %d when `kernel` is "%s". The '
                'interpolant may not be uniquely solvable, and the smoothing '
                'parameter may have an unintuitive effect.' %
                (min_degree, kernel),
                UserWarning
                )

    if k is None:
        nobs = ny
    else:
        # make sure the number of nearest neighbors used for interpolation does
        # not exceed the number of observations
        k = int(min(k, ny))
        nobs = k

    # The polynomial matrix must have full column rank in order for the
    # interpolant to be well-posed, which is not possible if there are fewer
    # observations than monomials
    nmonos = int(binom(degree + ndim, ndim))
    if nmonos > nobs:
        raise ValueError(
            'At least %d data points are required when `degree` is %d and the '
            'number of dimensions is %d' % (nmonos, degree, ndim)
            )

    return y, d, smoothing, kernel, epsilon, degree, k


class RBFInterpolator:
    """
    Radial basis function (RBF) interpolation in N dimensions

    Parameters
    ----------
    y : (P, N) array_like
        Data point coordinates
    d : (P, ...) array_like
        Data values at `y`
    smoothing : float or (P,) array_like, optional
        Smoothing parameter. The interpolant perfectly fits the data when this
        is set to 0. For large values, the interpolant approaches a least
        squares fit of a polynomial with the specified degree.
    kernel : str, optional
        Type of RBF. This should be one of:

            - 'linear'                       : ``-r``
            - 'tps' (thin plate spline)      : ``r**2 * log(r)``
            - 'cubic'                        : ``r**3``
            - 'quintic'                      : ``-r**5``
            - 'mq' (multiquadratic)          : ``-sqrt(1 + r**2)``
            - 'imq' (inverse multiquadratic) : ``1/sqrt(1 + r**2)``
            - 'iq' (inverse quadratic)       : ``1/(1 + r**2)``
            - 'ga' (Gaussian)                : ``exp(-r**2)``

    epsilon : float, optional
        Shape parameter that scales the input to the RBF. This can be ignored
        if `kernel` is 'linear', 'tps', 'cubic', or 'quintic' because it has
        the same effect as scaling the smoothing parameter. This must be
        specified if `kernel` is 'mq', 'imq', 'iq', or 'ga'.
    degree : int, optional
        Degree of the added polynomial. Some RBFs have a minimum polynomial
        degree that is needed for the interpolant to be well-posed. Those RBFs
        and their corresponding minimum degrees are:

            - 'mq'      : 0
            - 'linear'  : 0
            - 'tps'     : 1
            - 'cubic'   : 1
            - 'quintic' : 2

        The default value is the minimum degree required for `kernel` or 0 if
        there is no minimum required degree. Set this to -1 for no added
        polynomial.

    Notes
    -----
    An RBF is a scalar valued function in N-dimensional space whose value at
    :math:`x` can be expressed in terms of :math:`r=||x - c||`, where :math:`c`
    is the center of the RBF.

    An RBF interpolant for the vector of observations :math:`d`, which are made
    at the locations :math:`y`, is a linear combination of RBFs centered at
    :math:`y` plus a polynomial with a specified degree. The RBF interpolant is
    written as

    .. math::
        f(x) = K(x, y) a + P(x) b

    where :math:`K(x, y)` is a matrix of RBFs with centers at :math:`y`
    evaluated at the interpolation points :math:`x`, and :math:`P(x)` is a
    matrix of monomials, which span polynomials with the specified degree,
    evaluated at :math:`x`. The coefficients :math:`a` and :math:`b` are the
    solution to the linear equations

    .. math::
        (K(y, y) + \\lambda I) a + P(y) b = d

    and

    .. math::
        P(y)^T a = 0,

    where :math:`\\lambda` is a positive smoothing parameter that controls how
    well we want to fit the observations. The observations are fit exactly when
    the smoothing parameter is zero.

    For the RBFs 'ga', 'imq', and 'iq', the solution for :math:`a` and
    :math:`b` is analytically unique if :math:`P(y)` has full column rank. As
    an example, :math:`P(y)` would not have full column rank if the
    observations are collinear in two-dimensional space and the degree of the
    added polynomial is 1. For the RBFs 'mq', 'linear', 'tps', 'cubic', and
    'quintic', the solution for  :math:`a` and :math:`b` is analytically unique
    if :math:`P(y)` has full column rank and the degree of the added polynomial
    is not lower than the minimum value listed above (see Chapter 7 of [1]_ or
    [2]_).

    When using an RBF that is not scale invariant ('mq', 'imq', 'iq', and
    'ga'), an appropriate shape parameter must be chosen (e.g., through cross
    validation). Smaller values for the shape parameter correspond to wider
    RBFs. The problem can become ill-conditioned or singular when the shape
    parameter is too small.

    See Also
    --------
    KNearestRBFInterpolator

    References
    ----------
    .. [1] Fasshauer, G., 2007. Meshfree Approximation Methods with Matlab.
        World Scientific Publishing Co.

    .. [2] http://amadeus.math.iit.edu/~fass/603_ch3.pdf

    .. [3] Wahba, G., 1990. Spline Models for Observational Data. SIAM.

    .. [4] http://pages.stat.wisc.edu/~wahba/stat860public/lect/lect8/lect8.pdf

    Examples
    --------
    Demonstrate interpolating scattered data to a grid in 2-D

    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import RBFInterpolator
    >>> np.random.seed(0)

    >>> xobs = np.random.uniform(-1, 1, (100, 2))
    >>> yobs = np.sum(xobs, axis=1)*np.exp(-6*np.sum(xobs**2, axis=1))

    >>> xgrid = np.mgrid[-1:1:50j, -1:1:50j]
    >>> xflat = xgrid.reshape(2, -1).T
    >>> yflat = RBFInterpolator(xobs, yobs)(xflat)
    >>> ygrid = yflat.reshape(50, 50)

    >>> plt.pcolormesh(*xgrid, ygrid, vmin=-0.25, vmax=0.25, shading='gouraud')
    >>> plt.scatter(*xobs.T, c=yobs, s=50, ec='k', vmin=-0.25, vmax=0.25)
    >>> plt.colorbar()
    >>> plt.show()

    """
    def __init__(self, y, d,
                 smoothing=0.0,
                 kernel='tps',
                 epsilon=None,
                 degree=None):
        y, d, smoothing, kernel, epsilon, degree, _ = _sanitize_init_args(
            y, d, smoothing, kernel, epsilon, degree, None
            )

        ny, ndim = y.shape
        data_shape = d.shape[1:]
        d = d.reshape((ny, -1))
        powers = _monomial_powers(ndim, degree)
        lhs, rhs, shift, scale = _build_system(
            y, d, smoothing, kernel, epsilon, powers
            )

        coeffs = np.linalg.solve(lhs, rhs)

        self.y = y
        self.kernel = kernel
        self.epsilon = epsilon
        self.powers = powers
        self.shift = shift
        self.scale = scale
        self.coeffs = coeffs
        self.data_shape = data_shape

    def __call__(self, x):
        """
        Evaluates the interpolant at `x`

        Parameters
        ----------
        x : (Q, N) array_like
            Interpolation point coordinates

        Returns
        -------
        (Q, ...) ndarray
            Values of the interpolant at `x`

        """
        x = np.asarray(x, dtype=float, order='C')
        if x.ndim != 2:
            raise ValueError('Expected `x` to be a 2-dimensional array')

        nx, ndim = x.shape
        if ndim != self.y.shape[1]:
            raise ValueError(
                'Expected the second axis of `x` to have length %d' %
                self.y.shape[1]
                )

        out = _evaluate(
            x, self.y, self.kernel, self.epsilon, self.powers, self.shift,
            self.scale, self.coeffs
            )
        out = out.reshape((nx,) + self.data_shape)
        return out


class KNearestRBFInterpolator:
    """
    RBF interpolation using the k nearest neighbors

    Parameters
    ----------
    y : (P, N) array_like
        Data point coordinates
    d : (P, ...) array_like
        Data values at `y`
    k : int, optional
        Number of nearest data points to use for each interpolation point
    smoothing : float or (P,) array_like, optional
        Smoothing parameter. The interpolant perfectly fits the data when this
        is set to 0.
    kernel : str, optional
        Type of RBF. This should be one of:

            - 'linear'                       : ``-r``
            - 'tps' (thin plate spline)      : ``r**2 * log(r)``
            - 'cubic'                        : ``r**3``
            - 'quintic'                      : ``-r**5``
            - 'mq' (multiquadratic)          : ``-sqrt(1 + r**2)``
            - 'imq' (inverse multiquadratic) : ``1/sqrt(1 + r**2)``
            - 'iq' (inverse quadratic)       : ``1/(1 + r**2)``
            - 'ga' (Gaussian)                : ``exp(-r**2)``

    epsilon : float, optional
        Shape parameter that scales the input to the RBF. This can be ignored
        if `kernel` is 'linear', 'tps', 'cubic', or 'quintic' because it has
        the same effect as scaling the smoothing parameter. This must be
        specified if `kernel` is 'mq', 'imq', 'iq', or 'ga'.
    degree : int, optional
        Degree of the added polynomial. Some RBFs have a minimum polynomial
        degree that is needed for the interpolant to be well-posed. Those RBFs
        and their corresponding minimum degrees are:

            - 'mq'      : 0
            - 'linear'  : 0
            - 'tps'     : 1
            - 'cubic'   : 1
            - 'quintic' : 2

        The default value is the minimum degree required for `kernel` or 0 if
        there is no minimum required degree. Set this to -1 for no added
        polynomial.

    See Also
    --------
    RBFInterpolator

    """
    def __init__(self, y, d,
                 k=50,
                 smoothing=0.0,
                 kernel='tps',
                 epsilon=None,
                 degree=None):
        y, d, smoothing, kernel, epsilon, degree, k = _sanitize_init_args(
            y, d, smoothing, kernel, epsilon, degree, k
            )

        ny, ndim = y.shape
        data_shape = d.shape[1:]
        d = d.reshape((ny, -1))
        powers = _monomial_powers(ndim, degree)
        tree = cKDTree(y)

        self.y = y
        self.d = d
        self.k = k
        self.smoothing = smoothing
        self.kernel = kernel
        self.epsilon = epsilon
        self.powers = powers
        self.tree = tree
        self.data_shape = data_shape

    def __call__(self, x):
        """
        Evaluates the interpolant at `x`

        Parameters
        ----------
        x : (Q, N) array_like
            Interpolation point coordinates

        Returns
        -------
        (Q, ...) ndarray
            Values of the interpolant at `x`

        """
        x = np.asarray(x, dtype=float, order='C')
        if x.ndim != 2:
            raise ValueError('Expected `x` to be a 2-dimensional array')

        nx, ndim = x.shape
        if ndim != self.y.shape[1]:
            raise ValueError(
                'Expected the second axis of `x` to have length %d' %
                self.y.shape[1]
                )

        dtype = complex if np.iscomplexobj(self.d) else float
        out = np.zeros((nx,) + self.data_shape, dtype=dtype)

        # get the indices of the k nearest observation points to each
        # interpolation point
        _, yindices = self.tree.query(x, self.k)
        if self.k == 1:
            # cKDTree squeezes the output when k=1
            yindices = yindices[:, None]

        # multiple interpolation points may have the same neighborhood of
        # observation points. Make the neighborhoods unique so that we only
        # compute the interpolation coefficients once for each neighborhood
        yindices = np.sort(yindices, axis=1)
        yindices, inv = np.unique(yindices, return_inverse=True, axis=0)
        # `inv` tells us which neighborhood will be used by each interpolation
        # point. Now we find which interpolation points will be using each
        # neighborhood
        xindices = [[] for _ in range(len(yindices))]
        for i, j in enumerate(inv):
            xindices[j].append(i)

        for xidx, yidx in zip(xindices, yindices):
            # `yidx` are the indices of the observations in this neighborhood.
            # `xidx` are the indices of the interpolation points that are using
            # this neighborood
            xnbr = x[xidx]
            ynbr = self.y[yidx]
            dnbr = self.d[yidx]
            snbr = self.smoothing[yidx]
            lhs, rhs, shift, scale = _build_system(
                ynbr, dnbr, snbr, self.kernel, self.epsilon, self.powers,
                )

            coeffs = np.linalg.solve(lhs, rhs)

            out[xidx] = _evaluate(
                xnbr, ynbr, self.kernel, self.epsilon, self.powers, shift,
                scale, coeffs
                ).reshape((-1,) + self.data_shape)

        return out
