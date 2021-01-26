"""Module for RBF interpolation"""
import warnings
from functools import lru_cache
from itertools import combinations_with_replacement

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from scipy.special import xlogy


__all__ = ['RBFInterpolator', 'KNearestRBFInterpolator']


def _distance(x, y):
    """
    Returns a distance matrix between `x` and `y`

    Parameters
    ----------
    x : (..., n, d) ndarray

    y : (..., m, d) ndarray

    Returns
    -------
    (..., n, m) ndarray

    """
    if (x.ndim == 2) & (y.ndim == 2):
        # if possible, use the faster function
        return cdist(x, y)
    else:
        return np.linalg.norm(x[..., None, :] - y[..., None, :, :], axis=-1)


def _linear(r):
    """linear / 1st order polyharmonic spline"""
    return -r


def _tps(r):
    """thin plate spline / 2nd order polyharmonic spline"""
    return xlogy(r**2, r)


def _cubic(r):
    """cubic / 3rd order polyharmonic spline"""
    return r*r*r # faster than r**3


def _quintic(r):
    """quintic / 5th order polyharmonic spline"""
    return -r*r*r*r*r # faster than r**5


def _mq(r):
    """multiquadratic"""
    return -np.sqrt(r**2 + 1)


def _imq(r):
    """inverse multiquadratic"""
    return 1/np.sqrt(r**2 + 1)


def _iq(r):
    """inverse quadratic"""
    return 1/(r**2 + 1)


def _ga(r):
    """gaussian"""
    return np.exp(-r**2)


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


def _vandermonde(x, degree):
    """
    Returns monomials evaluated at `x`. The monomials span the space of
    polynomials with the specified degree

    Parameters
    ----------
    x : (..., d) float array

    degree : int

    Returns
    -------
    (..., p) float array

    """
    pwr = _monomial_powers(x.shape[-1], degree)
    out = np.product(x[..., None, :]**pwr, axis=-1)
    return out


# Interpolation with conditionally positive definite RBFs has no assurances of
# being well-posed when the order of the added polynomial is not high enough.
# Define that minimum polynomial order here. These values are from Chapter 8 of
# Fasshauer's "Meshfree Approximation Methods with MATLAB". The RBFs that are
# not in this dictionary are positive definite and do not need a polynomial
_NAME_TO_MIN_DEGREE = {
    'mq': 0,
    'linear': 0,
    'tps': 1,
    'cubic': 1,
    'quintic': 2
    }


_NAME_TO_FUNC = {
    'mq': _mq,
    'linear': _linear,
    'tps': _tps,
    'cubic': _cubic,
    'quintic': _quintic,
    'imq': _imq,
    'iq': _iq,
    'ga': _ga
    }


class RBFInterpolator:
    """
    Radial basis function (RBF) interpolation in N dimensions

    Parameters
    ----------
    y : (P, N) array_like
        Data point coordinates

    d : (P,) or (P, M) array_like
        Data values at `y`. The data can be scalar or vector valued. If the
        data are vector valued, an RBF interpolant will be fit to each
        component.

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
        Shape parameter that scales the input to the RBF. Unless using a
        polyharmonic spline ('linear', 'tps', 'cubic', or 'quintic'), this
        should be tuned for each application.

    degree : int, optional
        Degree of the added polynomial. Set this to -1 for no added polynomial.
        For some choices of `kernel`, the interpolant may not be well-posed if
        `degree` is below some value. Those RBFs and their corresponding
        minimum degrees are:

            - 'mq'      : 0
            - 'linear'  : 0
            - 'tps'     : 1
            - 'cubic'   : 1
            - 'quintic' : 2

        The default value is the minimum degree required for `kernel` or 0 if
        there is no minimum required degree.

    Notes
    -----
    An RBF is a scalar valued function in N-dimensional space whose value at
    :math:`x` can be expressed in terms of :math:`r=||x - c||`, where :math:`c`
    is the center of the RBF.

    An RBF interpolant for the vector of observations :math:`d`, which are made
    at locations :math:`y`, is a linear combination of RBFs centered at
    :math:`y` plus a polynomial with a specified degree. The RBF interpolant is
    written as

    .. math::
        f(x) = K(x, y) a + P(x) b

    where :math:`K(x, y)` is a matrix of RBFs with centers at :math:`y`
    evaluated at the interpolation points :math:`x`, and :math:`P(x)` is a
    matrix of monomials evaluated at :math:`x`. The coefficients :math:`a` and
    :math:`b` are the solution to the linear equations

    .. math::
        (K(y, y) + \lambda I) a + P(y) b = d

    and

    .. math::
        P(y)^T a = 0,

    where :math:`\lambda` is a positive smoothing parameter that controls how
    well we want to fit the observations. The observations are fit exactly when
    the smoothing parameter is zero.

    For the RBFs 'ga', 'imq', and 'iq', the solution for :math:`a` and
    :math:`b` is unique as long as :math:`P(y)` has full rank. As an example,
    :math:`P(y)` would not have full rank if the observations are collinear in
    two-dimensional space and the degree of the added polynomial is 1. For the
    remaining RBFs, the solution for :math:`a` and :math:`b` is unique as long
    as :math:`P(y)` has full rank and the degree of the added polynomial is not
    lower than the minimum value listed above (see Chapter 7 of [1]_ or [2]_).

    References
    ----------
    .. [1] Fasshauer, G., 2007. Meshfree Approximation Methods with Matlab.
        World Scientific Publishing Co.

    .. [2] http://amadeus.math.iit.edu/~fass/603_ch3.pdf

    .. [3] Wahba, G., 1990. Spline Models for Observational Data. SIAM.

    .. [4] http://pages.stat.wisc.edu/~wahba/stat860public/lect/lect8/lect8.pdf

    """
    def __init__(self, y, d,
                 smoothing=0.0,
                 kernel='cubic',
                 epsilon=1.0,
                 degree=None):
        y = np.asarray(y, dtype=float)
        if y.ndim != 2:
            raise ValueError('Expected `y` to be a 2-dimensional array')

        ny = y.shape[0]
        d = np.asarray(d, dtype=float)
        if (d.ndim != 1) & (d.ndim != 2):
            raise ValueError('Expected `d` to be a 1 or 2-dimensional array')

        if d.shape[0] != ny:
            raise ValueError(
                'Expected the first axis of `d` to have length %d' % ny
                )

        data_shape = d.shape[1:]

        if np.isscalar(smoothing):
            smoothing = np.full(ny, smoothing, dtype=float)
        else:
            smoothing = np.asarray(smoothing, dtype=float)
            if smoothing.shape != (ny,):
                raise ValueError(
                    'Expected `smoothing` to be a scalar or have shape (%d,)'
                    % ny
                    )

        try:
            kernel_func = _NAME_TO_FUNC[kernel]
        except KeyError:
            raise ValueError(
                '`kernel` is %s but it should be one of {%s}' %
                (kernel, ', '.join(_NAME_TO_FUNC.keys()))
                )

        if not np.isscalar(epsilon):
            raise ValueError('Expected `epsilon` to be a scalar')

        min_degree = _NAME_TO_MIN_DEGREE.get(kernel, -1)
        if degree is None:
            degree = max(min_degree, 0)
        elif degree < min_degree:
            warnings.warn(
                'The polynomial degree should not be below %d for %s in order '
                'for the interpolant to be uniquely solvable' %
                (min_degree, kernel)
                )

        degree = int(degree)
        # For improved numerical stability, shift the observations so that
        # their centroid is zero
        center = y.mean(axis=0)
        y = y - center
        # Build the system of equations and solve for the RBF and monomial
        # coefficients
        Kyy = kernel_func(_distance(y, y)*epsilon)
        Kyy[range(ny), range(ny)] += smoothing
        Py = _vandermonde(y, degree)
        nmonos = Py.shape[1]
        if nmonos > ny:
            raise ValueError(
                'The polynomial degree is too high. The number of monomials, '
                '%d, exceeds the number of observations, %d' % (nmonos, ny)
                )

        Z = np.zeros((nmonos, nmonos), dtype=float)
        LHS = np.block([[Kyy, Py], [Py.T, Z]])
        z = np.zeros((nmonos,) + data_shape, dtype=float)
        rhs = np.concatenate((d, z), axis=0)
        coeff = np.linalg.solve(LHS, rhs)
        kernel_coeff, poly_coeff = coeff[:ny], coeff[ny:]

        self.y = y
        self.kernel = kernel
        self.kernel_func = kernel_func
        self.epsilon = epsilon
        self.degree = degree
        self.center = center
        self.kernel_coeff = kernel_coeff
        self.poly_coeff = poly_coeff
        self.data_shape = data_shape

    def __call__(self, x, chunk_size=1000):
        """
        Evaluates the interpolant at `x`

        Parameters
        ----------
        x : (Q, N) array_like
            Interpolation point coordinates

        chunk_size : int, optional
            Break `x` into chunks with this size and evaluate the interpolant
            for each chunk

        Returns
        -------
        (Q,) or (Q, M) ndarray
            Values of the interpolant at `x`

        """
        x = np.asarray(x, dtype=float)
        if x.ndim != 2:
            raise ValueError('Expected `x` to be a 2-dimensional array')

        nx = x.shape[0]
        if chunk_size is not None:
            out = np.zeros((nx,) + self.data_shape, dtype=float)
            for start in range(0, nx, chunk_size):
                stop = start + chunk_size
                out[start:stop] = self(x[start:stop], chunk_size=None)

            return out

        x = x - self.center
        Kxy = self.kernel_func(_distance(x, self.y)*self.epsilon)
        Px = _vandermonde(x, self.degree)
        out = Kxy.dot(self.kernel_coeff) + Px.dot(self.poly_coeff)
        return out


class KNearestRBFInterpolator:
    """
    RBF interpolation using the k nearest neighbors

    Parameters
    ----------
    y : (P, N) array_like
        Data point coordinates

    d : (P,) or (P, M) array_like
        Data values at `y`. The data can be scalar or vector valued. If the
        data are vector valued, an RBF interpolant will be fit to each
        component.

    k : int, optional
        Number of nearest neighbors to use for each interpolation point

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
        Shape parameter that scales the input to the RBF. Unless using a
        polyharmonic spline ('linear', 'tps', 'cubic', or 'quintic'), this
        should be tuned for each application.

    degree : int, optional
        Degree of the added polynomial. Set this to -1 for no added polynomial.
        For some choices of `kernel`, the interpolant may not be well-posed if
        `degree` is below some value. Those RBFs and their corresponding
        minimum degrees are:

            - 'mq'      : 0
            - 'linear'  : 0
            - 'tps'     : 1
            - 'cubic'   : 1
            - 'quintic' : 2

        The default value is the minimum degree required for `kernel` or 0 if
        there is no minimum required degree.

    """
    def __init__(self, y, d,
                 smoothing=0.0,
                 k=50,
                 kernel='cubic',
                 epsilon=1.0,
                 degree=None):
        y = np.asarray(y, dtype=float)
        if y.ndim != 2:
            raise ValueError('Expected `y` to be a 2-dimensional array')

        ny = y.shape[0]
        d = np.asarray(d, dtype=float)
        if (d.ndim != 1) & (d.ndim != 2):
            raise ValueError('Expected `d` to be a 1 or 2-dimensional array')

        if d.shape[0] != ny:
            raise ValueError(
                'Expected the first axis of `d` to have length %d' % ny
                )

        data_shape = d.shape[1:]

        if np.isscalar(smoothing):
            smoothing = np.full(ny, smoothing, dtype=float)
        else:
            smoothing = np.asarray(smoothing, dtype=float)
            if smoothing.shape != (ny,):
                raise ValueError(
                    'Expected `smoothing` to be a scalar or have shape (%d,)'
                    % ny
                    )

        # make sure the number of nearest neighbors used for interpolation does
        # not exceed the number of observations
        k = min(int(k), ny)
        try:
            kernel_func = _NAME_TO_FUNC[kernel]
        except KeyError:
            raise ValueError(
                '`kernel` is %s but it should be one of {%s}' %
                (kernel, ', '.join(_NAME_TO_FUNC.keys()))
                )

        if not np.isscalar(epsilon):
            raise ValueError('Expected `epsilon` to be a scalar')

        min_degree = _NAME_TO_MIN_DEGREE.get(kernel, -1)
        if degree is None:
            degree = max(min_degree, 0)
        elif degree < min_degree:
            warnings.warn(
                'The polynomial degree should not be below %d for %s in order '
                'for the interpolant to be uniquely solvable' %
                (min_degree, kernel)
                )

        degree = int(degree)
        tree = cKDTree(y)

        self.y = y
        self.d = d
        self.smoothing = smoothing
        self.k = k
        self.kernel = kernel
        self.kernel_func = kernel_func
        self.epsilon = epsilon
        self.degree = degree
        self.tree = tree
        self.data_shape = data_shape

    def __call__(self, x, chunk_size=1000):
        """
        Evaluates the interpolant at `x`

        Parameters
        ----------
        x : (Q, N) array_like
            Interpolation point coordinates

        chunk_size : int, optional
            Break `x` into chunks with this size and evaluate the interpolant
            for each chunk

        Returns
        -------
        (Q,) or (Q, M) ndarray
            Values of the interpolant at `x`

        """
        x = np.asarray(x, dtype=float)
        if x.ndim != 2:
            raise ValueError('Expected `x` to be a 2-dimensional array')

        nx = x.shape[0]
        if chunk_size is not None:
            out = np.zeros((nx,) + self.data_shape, dtype=float)
            for start in range(0, nx, chunk_size):
                stop = start + chunk_size
                out[start:stop] = self(x[start:stop], chunk_size=None)

            return out

        # get the indices of the k nearest observations for each interpolation
        # point
        _, nbr = self.tree.query(x, self.k)
        if self.k == 1:
            # cKDTree squeezes the output when k=1
            nbr = nbr[:, None]

        # multiple interpolation points may have the same neighborhood. Make
        # the neighborhoods unique so that we only compute the interpolation
        # coefficients once for each neighborhood
        nbr, inv = np.unique(np.sort(nbr, axis=1), return_inverse=True, axis=0)
        nnbr = nbr.shape[0]
        # Get the observation data for each neighborhood
        y, d, smoothing = self.y[nbr], self.d[nbr], self.smoothing[nbr]
        # Shift the centers of the neighborhoods to zero for improved numerical
        # stability
        centers = y.mean(axis=1)
        y = y - centers[:, None]
        # build the left-hand-side interpolation matrix consisting of the RBF
        # and monomials evaluated at each neighborhood
        Kyy = self.kernel_func(_distance(y, y)*self.epsilon)
        Kyy[:, range(self.k), range(self.k)] += smoothing
        Py = _vandermonde(y, self.degree)
        nmonos = Py.shape[2]
        if nmonos > self.k:
            raise ValueError(
                'The polynomial degree is too high. The number of monomials, '
                '%d, exceeds the number of neighbors used for interpolation, '
                '%d' % (nmonos, self.k)
                )

        PyT = np.transpose(Py, (0, 2, 1))
        Z = np.zeros((nnbr, nmonos, nmonos), dtype=float)
        LHS = np.block([[Kyy, Py], [PyT, Z]])
        # build the right-hand-side data vector consisting of the observations
        # for each neighborhood and extra zeros
        z = np.zeros((nnbr, nmonos) + self.data_shape, dtype=float)
        rhs = np.concatenate((d, z), axis=1)
        # solve for the RBF and polynomial coefficients for each neighborhood
        coeff = np.linalg.solve(LHS, rhs)
        # expand the arrays from having one entry per neighborhood to one entry
        # per interpolation point
        coeff = coeff[inv]
        y = y[inv]
        centers = centers[inv]
        # evaluate at the interpolation points
        x = x - centers
        kernel_coeff = coeff[:, :self.k]
        poly_coeff = coeff[:, self.k:]
        Kxy = self.kernel_func(_distance(x[:, None], y)*self.epsilon)[:, 0]
        Px = _vandermonde(x, self.degree)
        if self.data_shape:
            # if the data are vector valued, add an extra axis to Kxy and Px to
            # allow for broadcasting
            Kxy = Kxy[:, :, None]
            Px = Px[:, :, None]

        out = (Kxy*kernel_coeff).sum(axis=1) + (Px*poly_coeff).sum(axis=1)
        return out
